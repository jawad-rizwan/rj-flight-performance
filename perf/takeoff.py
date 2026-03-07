"""
Takeoff Analysis (Raymer 17.8).

Covers ground roll, rotation, transition, climb,
and balanced field length.
"""

import numpy as np
from .atmosphere import isa_density, G, RHO_SL
from .level_flight import V_stall


# =========================================================================
# Ground Roll
# =========================================================================

def ground_roll_distance(W, S, T, CD0, CL_ground, K, mu, rho, CL_max_TO):
    """Takeoff ground-roll distance [ft].

    Raymer Eq. (17.102):
        S_G = (1/(2g*K_A)) * ln( (K_T + K_A*V_TO^2) / (K_T + K_A*V_i^2) )

    where V_i = 0, V_TO = 1.1*V_stall, and:
        K_T = (T/W) - mu                          Eq. (17.103)
        K_A = (rho/(2*W/S)) * (mu*CL - CD0 - K*CL^2)   Eq. (17.104)

    Note K_A is typically negative (drag > friction relief from lift).

    Parameters
    ----------
    W          : float, takeoff weight [lb]
    S          : float, wing area [ft^2]
    T          : float, average thrust during ground roll [lb]
    CD0        : float, zero-lift drag coefficient (ground config)
    CL_ground  : float, ground-roll lift coefficient (~0.1)
    K          : float, induced drag factor
    mu         : float, rolling friction coefficient
    rho        : float, air density [slug/ft^3]
    CL_max_TO  : float, max CL in takeoff configuration

    Returns
    -------
    S_G : float, ground roll distance [ft]
    """
    V_s = V_stall(W, S, rho, CL_max_TO)
    V_TO = 1.1 * V_s                        # takeoff speed = 1.1 * V_stall

    K_T = (T / W) - mu                      # Eq. (17.103)
    K_A = (rho / (2.0 * W / S)) * (mu * CL_ground - CD0 - K * CL_ground**2)  # Eq. (17.104)

    if abs(K_A) < 1e-10:
        # No aero contribution, simple kinematics
        return V_TO**2 / (2.0 * G * K_T)

    S_G = (1.0 / (2.0 * G * K_A)) * np.log(
        (K_T + K_A * V_TO**2) / (K_T + K_A * 0.0**2))  # Eq. (17.102)

    return S_G


# =========================================================================
# Rotation
# =========================================================================

def rotation_distance(V_TO, t_rotate=3.0):
    """Rotation distance [ft].

    S_R ~ V_TO * t_rotate
    Typical t_rotate = 3 s for large aircraft, 1 s for small.

    Parameters
    ----------
    V_TO      : float, takeoff speed [ft/s]
    t_rotate  : float, rotation time [s]
    """
    return V_TO * t_rotate


# =========================================================================
# Transition
# =========================================================================

def transition_segment(W, S, rho, CL_max_TO, TW):
    """Transition arc parameters after rotation.

    Raymer Eqs. (17.105)-(17.111).

    During transition, V increases from V_TO to V_climb (1.2*V_stall).
    Average velocity ~ 1.15 * V_stall.
    Average CL ~ 0.9 * CL_max_TO.

    Returns
    -------
    dict with keys:
        V_TR    : transition velocity [ft/s]
        R       : transition arc radius [ft]
        gamma   : climb angle at end of transition [rad]
        S_TR    : horizontal distance during transition [ft]
        h_TR    : altitude gained during transition [ft]
    """
    V_s = V_stall(W, S, rho, CL_max_TO)
    V_TR = 1.15 * V_s                      # average transition velocity

    # Load factor during transition
    # Eq. (17.105)-(17.106): n ~ 1.2
    n = 1.2

    # Transition arc radius
    # Eq. (17.107): R = V_TR^2 / (g*(n-1)) = V_TR^2 / (0.2*g)
    R = V_TR**2 / (G * (n - 1.0))          # Eq. (17.107)

    # Climb angle at end of transition
    # Eq. (17.108): sin(gamma) = (T-D)/W ~ T/W - 1/(L/D)
    LD_climb = 0.9 * CL_max_TO / (0.9 * CL_max_TO * 0.05 + 0.022)  # rough estimate
    gamma = np.arcsin(np.clip(TW - 1.0 / LD_climb, -1, 1))

    # Use simplified: sin(gamma) ~ T/W - 1/(L/D)
    sin_gamma = TW - 1.0 / max(LD_climb, 1.0)
    sin_gamma = np.clip(sin_gamma, 0.0, 1.0)
    gamma = np.arcsin(sin_gamma)                # Eq. (17.108)

    # Horizontal distance
    S_TR = R * sin_gamma                         # Eq. (17.109)

    # Altitude gained
    h_TR = R * (1.0 - np.cos(gamma))            # Eq. (17.110)

    return {
        "V_TR": V_TR,
        "R": R,
        "gamma_rad": gamma,
        "gamma_deg": np.degrees(gamma),
        "S_TR": S_TR,
        "h_TR": h_TR,
    }


# =========================================================================
# Climb segment (to clear obstacle)
# =========================================================================

def climb_distance(h_obstacle, h_TR, gamma):
    """Horizontal distance during climb to clear obstacle [ft].

    Raymer Eq. (17.112):  S_c = (h_obstacle - h_TR) / tan(gamma)

    If h_TR >= h_obstacle, obstacle cleared during transition and S_c = 0.
    """
    if h_TR >= h_obstacle:
        return 0.0
    return (h_obstacle - h_TR) / np.tan(gamma)   # Eq. (17.112)


# =========================================================================
# Total takeoff distance
# =========================================================================

def total_takeoff_distance(W, S, T, CD0, CL_ground, K, mu, rho,
                           CL_max_TO, TW, h_obstacle=35.0, t_rotate=3.0):
    """Total takeoff distance over obstacle [ft].

    Sum of: ground roll + rotation + transition + climb

    Parameters
    ----------
    W, S, T, CD0, CL_ground, K, mu, rho : see ground_roll_distance
    CL_max_TO : float, max CL in takeoff config
    TW        : float, T/W at takeoff
    h_obstacle: float, obstacle height [ft] (35 ft FAR 25)
    t_rotate  : float, rotation time [s]

    Returns
    -------
    dict with all segment distances and total
    """
    V_s = V_stall(W, S, rho, CL_max_TO)
    V_TO = 1.1 * V_s

    S_G = ground_roll_distance(W, S, T, CD0, CL_ground, K, mu, rho, CL_max_TO)
    S_R = rotation_distance(V_TO, t_rotate)

    trans = transition_segment(W, S, rho, CL_max_TO, TW)
    S_TR = trans["S_TR"]
    h_TR = trans["h_TR"]
    gamma = trans["gamma_rad"]

    S_C = climb_distance(h_obstacle, h_TR, gamma) if gamma > 0.01 else 0.0

    S_total = S_G + S_R + S_TR + S_C

    return {
        "V_stall": V_s,
        "V_TO": V_TO,
        "S_ground_roll": S_G,
        "S_rotation": S_R,
        "S_transition": S_TR,
        "h_transition": h_TR,
        "S_climb": S_C,
        "gamma_deg": trans["gamma_deg"],
        "S_total": S_total,
    }


# =========================================================================
# Balanced Field Length
# =========================================================================

def balanced_field_length(W_over_S, CL_climb, h_obstacle, TW, U,
                          G_climb, BPR=None, N_e=2, D_p=0.0, bhp=0.0,
                          rho=RHO_SL, is_prop=False):
    """Balanced field length (empirical method).

    Raymer Eq. (17.113):
        BFL = 0.863/(1+2.3*G) * (W/S/(rho*g*CL_climb) + h_obstacle)
              * (1/(T_av/W - U) + 2.7) + (655/sqrt(rho/rho_SL))

    Parameters
    ----------
    W_over_S  : float, takeoff wing loading [lb/ft^2]
    CL_climb  : float, CL at climb speed (CL at 1.2*V_stall)
    h_obstacle: float, obstacle height [ft]
    TW        : float, all-engine T/W at takeoff
    U         : float, = 0.01*CL_max_TO + 0.02 (flaps in TO position)
    G_climb   : float, climb gradient = gamma_climb - gamma_min
                gamma_min = 0.024 (2-engine), 0.027 (3), 0.030 (4)
    BPR       : float, bypass ratio (jet only)
    N_e       : int, number of engines
    D_p       : float, propeller diameter [ft] (prop only)
    bhp       : float, total brake horsepower (prop only)
    rho       : float, air density [slug/ft^3]
    is_prop   : bool, True if propeller aircraft

    Returns
    -------
    BFL : float, balanced field length [ft]
    """
    sig = rho / RHO_SL

    # Average thrust during takeoff
    if is_prop:
        # Eq. (17.115)
        T_av = 5.75 * bhp * ((rho / RHO_SL) * N_e * D_p**2 / bhp)**(1.0/3.0)
        TW_av = T_av / (W_over_S * 1.0)  # placeholder; need actual W
    else:
        # Eq. (17.114): T_av = 0.75 * T_static * (5+BPR)/(4+BPR)
        T_ratio = 0.75 * (5.0 + BPR) / (4.0 + BPR) if BPR else 0.75
        TW_av = TW * T_ratio

    BFL = (0.863 / (1.0 + 2.3 * G_climb)) * (
        W_over_S / (rho * G * CL_climb) + h_obstacle
    ) * (1.0 / (TW_av - U) + 2.7) + 655.0 / np.sqrt(sig)   # Eq. (17.113)

    return BFL
