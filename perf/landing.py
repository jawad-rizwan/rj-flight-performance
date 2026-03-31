"""
Landing Analysis (Raymer 17.9).

Covers approach, flare, free roll, and braking ground roll.
"""

import numpy as np
from .atmosphere import G, isa_density, RHO_SL
from .level_flight import V_stall


# =========================================================================
# Approach
# =========================================================================

def approach_speed(W, S, rho, CL_max_L, factor=1.3):
    """Approach speed [ft/s].

    V_a = factor * V_stall (landing config).
    factor = 1.3 for civil (FAR), 1.2 for military.
    """
    return factor * V_stall(W, S, rho, CL_max_L)


def approach_angle(TW_idle, LD):
    """Approach angle [rad].

    Raymer Eq. (17.108) applied to approach:
        sin(gamma_a) = T/W - 1/(L/D)

    For idle thrust, T/W ~ 0.03-0.05 for jets.
    Transport approach angle should be <= 3 deg.
    """
    sin_g = TW_idle - 1.0 / LD
    sin_g = np.clip(sin_g, -1.0, 0.0)  # approach is descending
    return np.arcsin(sin_g)


def approach_distance(h_obstacle, gamma_a):
    """Horizontal approach distance from obstacle to flare [ft].

    S_a = h_obstacle / tan(|gamma_a|)

    Uses flare height h_f as the reference if needed.
    """
    return h_obstacle / np.tan(abs(gamma_a))


# =========================================================================
# Flare
# =========================================================================

def flare_parameters(W, S, rho, CL_max_L, approach_factor=1.3,
                     glideslope_deg=3.0):
    """Flare segment parameters.

    Touchdown speed V_TD = 1.15 * V_stall (civil) or 1.1 * V_stall (military).
    Average flare velocity V_f = 1.23 * V_stall.
    Flare load factor n = 1.2.
    Flare radius R = V_f^2 / (g*(n-1)) = V_f^2 / (0.2*g).

    Raymer Eq. (17.107) applied to flare.

    Parameters
    ----------
    glideslope_deg : float, approach glideslope angle [deg].
        Standard = 3.0.  Steep approaches use 4.5-5.5 deg (e.g. London City 5.5).

    Returns
    -------
    dict with V_TD, V_f, R, h_f, S_f, glideslope_deg
    """
    V_s = V_stall(W, S, rho, CL_max_L)
    V_a = approach_factor * V_s
    V_TD = 1.15 * V_s                # civil touchdown speed
    V_f = 1.23 * V_s                 # average flare velocity (Raymer 17.9.2)

    n = 1.2
    R = V_f**2 / (G * (n - 1.0))     # Eq. (17.107) for flare

    gamma_a = np.radians(glideslope_deg)

    # Flare height
    h_f = R * (1.0 - np.cos(gamma_a))    # Eq. (17.110) for flare

    # Flare distance
    S_f = R * np.sin(gamma_a)              # Eq. (17.109) for flare

    return {
        "V_stall": V_s,
        "V_approach": V_a,
        "V_TD": V_TD,
        "V_f": V_f,
        "R": R,
        "h_f": h_f,
        "S_f": S_f,
        "glideslope_deg": glideslope_deg,
    }


# =========================================================================
# Free Roll
# =========================================================================

def free_roll_distance(V_TD, t_free=3.0):
    """Free-roll distance after touchdown before braking [ft].

    S_free = V_TD * t_free
    Typical t_free = 1-3 s.
    """
    return V_TD * t_free


# =========================================================================
# Braking Ground Roll
# =========================================================================

def braking_distance(W, S, V_TD, mu_brake, rho, CD0, CL_ground, K,
                     T_idle=0.0, T_reverse=0.0):
    """Braking ground-roll distance [ft].

    Uses Eq. (17.102) run in reverse: initial V = V_TD, final V = 0.

    S_B = (1/(2g*K_A)) * ln( (K_T + K_A*0) / (K_T + K_A*V_TD^2) )

    K_T = (T_idle/W - mu_brake) for no thrust reverser
    K_A = rho/(2*W/S) * (mu_brake*CL_ground - CD0 - K*CL_ground^2)

    With thrust reversers: T_idle is replaced by -T_reverse (negative = retarding).

    Parameters
    ----------
    W          : float, landing weight [lb]
    S          : float, wing area [ft^2]
    V_TD       : float, touchdown speed [ft/s]
    mu_brake   : float, braking friction coefficient
    rho        : float, air density [slug/ft^3]
    CD0        : float, zero-lift drag (landing config)
    CL_ground  : float, ground CL
    K          : float, induced drag factor
    T_idle     : float, idle thrust [lb] (positive forward)
    T_reverse  : float, reverse thrust magnitude [lb] (positive = retarding)

    Returns
    -------
    S_B : float, braking distance [ft]
    """
    T_net = T_idle - T_reverse   # net forward thrust (negative if reversers)
    K_T = (T_net / W) - mu_brake                    # Eq. (17.103) modified
    K_A = (rho / (2.0 * W / S)) * (
        mu_brake * CL_ground - CD0 - K * CL_ground**2)  # Eq. (17.104)

    if abs(K_A) < 1e-10:
        return -V_TD**2 / (2.0 * G * K_T)

    S_B = (1.0 / (2.0 * G * K_A)) * np.log(K_T / (K_T + K_A * V_TD**2))

    return abs(S_B)


# =========================================================================
# Total Landing Distance
# =========================================================================

def total_landing_distance(W, S, rho, CL_max_L, CD0, CL_ground, K,
                           mu_brake, h_obstacle=50.0, T_idle=0.0,
                           T_reverse=0.0, approach_factor=1.3,
                           t_free=3.0, FAR_factor=True,
                           glideslope_deg=3.0):
    """Total landing distance [ft].

    Sum of: approach + flare + free roll + braking

    Parameters
    ----------
    FAR_factor     : bool, if True multiply by 1/0.6 for FAR field length
    glideslope_deg : float, approach glideslope angle [deg].
        Standard = 3.0.  Steep approaches use 4.5-5.5 deg.

    Returns
    -------
    dict with all segment distances and total
    """
    flare = flare_parameters(W, S, rho, CL_max_L, approach_factor,
                             glideslope_deg)

    # Approach
    gamma_a = np.radians(glideslope_deg)
    S_a = (h_obstacle - flare["h_f"]) / np.tan(gamma_a) if flare["h_f"] < h_obstacle else 0.0

    # Flare
    S_f = flare["S_f"]

    # Free roll
    S_free = free_roll_distance(flare["V_TD"], t_free)

    # Braking
    S_B = braking_distance(W, S, flare["V_TD"], mu_brake, rho,
                           CD0, CL_ground, K, T_idle, T_reverse)

    S_total = S_a + S_f + S_free + S_B

    # FAR field length = actual / 0.6  (FAR 25.125)
    S_FAR = S_total / 0.6 if FAR_factor else S_total

    # Sink rate on approach
    V_sink = flare["V_approach"] * np.sin(gamma_a)

    return {
        "V_stall": flare["V_stall"],
        "V_approach": flare["V_approach"],
        "V_TD": flare["V_TD"],
        "glideslope_deg": glideslope_deg,
        "V_sink_approach": V_sink,
        "S_approach": S_a,
        "S_flare": S_f,
        "S_free_roll": S_free,
        "S_braking": S_B,
        "S_total_actual": S_total,
        "S_FAR_field": S_FAR,
    }


# =========================================================================
# Steep Approach (EASA CS 25.1420 / FAR 25.125)
# =========================================================================

def approach_sink_rate(V_approach, glideslope_deg):
    """Sink rate during approach [ft/s].

    V_sink = V_approach * sin(gamma_a)

    Typical limits:
        Normal (3 deg):   ~12 ft/s  (~720 fpm)
        Steep (5.5 deg):  ~22 ft/s  (~1300 fpm)

    Parameters
    ----------
    V_approach     : float, approach speed [ft/s]
    glideslope_deg : float, glideslope angle [deg]

    Returns
    -------
    V_sink : float, sink rate [ft/s] (positive downward)
    """
    return V_approach * np.sin(np.radians(glideslope_deg))


def max_glideslope_angle(W, S, rho, CL_max_L, h_obstacle=50.0,
                         approach_factor=1.3):
    """Maximum feasible glideslope angle [deg].

    Limited by the condition that flare height h_f must not exceed
    the obstacle clearance height h_obstacle.  Beyond this angle the
    aircraft would need to begin its flare above the screen height,
    making the approach geometry infeasible.

    h_f = R * (1 - cos(gamma_a))  <=  h_obstacle
    => gamma_max = arccos(1 - h_obstacle / R)

    Parameters
    ----------
    h_obstacle : float, screen height [ft] (default 50)

    Returns
    -------
    gamma_max_deg : float, maximum glideslope angle [deg]
    """
    V_s = V_stall(W, S, rho, CL_max_L)
    V_f = 1.23 * V_s
    n = 1.2
    R = V_f**2 / (G * (n - 1.0))

    ratio = h_obstacle / R
    if ratio >= 2.0:
        return 90.0  # practically unlimited
    return np.degrees(np.arccos(1.0 - ratio))


def steep_approach_analysis(W, S, rho, CL_max_L, CD0, CL_ground, K,
                            mu_brake, h_obstacle=50.0, T_idle=0.0,
                            T_reverse=0.0, approach_factor=1.3,
                            t_free=3.0, angles_deg=None):
    """Landing distances for standard and steep approach angles.

    Computes full landing distance breakdown for each glideslope angle
    and flags angles that exceed the flare-height feasibility limit.

    Steep approach operations (glideslope > 3 deg) are used at airports
    with terrain or noise constraints (e.g. London City at 5.5 deg,
    Lugano at 6.65 deg).  Certification basis: EASA CS 25.1420 /
    FAR 25.125(b)(2)(iii).

    Parameters
    ----------
    angles_deg : list of float, glideslope angles to evaluate [deg].
        Default [3.0, 4.5, 5.5].

    Returns
    -------
    dict keyed by angle (float) with landing distance dicts,
    plus 'gamma_max_deg' for the feasibility limit.
    """
    if angles_deg is None:
        angles_deg = [3.0, 4.5, 5.5]

    gamma_max = max_glideslope_angle(W, S, rho, CL_max_L, h_obstacle,
                                     approach_factor)

    results = {"gamma_max_deg": gamma_max}

    for angle in angles_deg:
        ld = total_landing_distance(
            W, S, rho, CL_max_L, CD0, CL_ground, K,
            mu_brake, h_obstacle, T_idle, T_reverse,
            approach_factor, t_free, FAR_factor=True,
            glideslope_deg=angle)
        ld["feasible"] = angle <= gamma_max
        results[angle] = ld

    return results
