"""
Takeoff Analysis (Raymer 17.8).

Covers ground roll, rotation, transition, climb,
balanced field length, accelerate-stop distance, and V-speed estimation.
"""

import numpy as np
from .atmosphere import isa_density, G, RHO_SL
from .level_flight import V_stall


# =========================================================================
# Ground Roll (to arbitrary speed)
# =========================================================================

def _ground_roll_to_V(W, S, T, CD0, CL_ground, K, mu, rho, V_i, V_f):
    """Ground-roll distance from V_i to V_f [ft].

    Raymer Eq. (17.102):
        S_G = (1/(2g*K_A)) * ln( (K_T + K_A*V_f^2) / (K_T + K_A*V_i^2) )

    K_T = (T/W) - mu                          Eq. (17.103)
    K_A = (rho/(2*W/S)) * (mu*CL - CD0 - K*CL^2)   Eq. (17.104)
    """
    K_T = (T / W) - mu
    K_A = (rho / (2.0 * W / S)) * (mu * CL_ground - CD0 - K * CL_ground**2)

    if abs(K_A) < 1e-10:
        return (V_f**2 - V_i**2) / (2.0 * G * K_T)

    return (1.0 / (2.0 * G * K_A)) * np.log(
        (K_T + K_A * V_f**2) / (K_T + K_A * V_i**2))


def ground_roll_distance(W, S, T, CD0, CL_ground, K, mu, rho, CL_max_TO):
    """Takeoff ground-roll distance from 0 to V_TO (1.1*V_stall) [ft].

    Raymer Eq. (17.102)-(17.104).
    """
    V_s = V_stall(W, S, rho, CL_max_TO)
    V_TO = 1.1 * V_s
    return _ground_roll_to_V(W, S, T, CD0, CL_ground, K, mu, rho, 0.0, V_TO)


# =========================================================================
# Braking Distance (from V to 0)
# =========================================================================

def _braking_distance_from_V(W, S, V, mu_brake, rho, CD0, CL_ground, K,
                              T_idle=0.0, T_reverse=0.0):
    """Braking distance from speed V to standstill [ft].

    Uses Eq. (17.102)-(17.104) run in reverse with braking friction.
    """
    T_net = T_idle - T_reverse
    K_T = (T_net / W) - mu_brake
    K_A = (rho / (2.0 * W / S)) * (
        mu_brake * CL_ground - CD0 - K * CL_ground**2)

    if abs(K_A) < 1e-10:
        return -V**2 / (2.0 * G * K_T)

    return abs((1.0 / (2.0 * G * K_A)) * np.log(K_T / (K_T + K_A * V**2)))


# =========================================================================
# Rotation
# =========================================================================

def rotation_distance(V_TO, t_rotate=3.0):
    """Rotation distance [ft].

    S_R ~ V_TO * t_rotate
    Typical t_rotate = 3 s for large aircraft, 1 s for small.
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
# Total Takeoff Distance (TODR — all engines operating)
# =========================================================================

def total_takeoff_distance(W, S, T, CD0, CL_ground, K, mu, rho,
                           CL_max_TO, TW, h_obstacle=35.0, t_rotate=3.0):
    """All-engine takeoff distance over obstacle (TODR) [ft].

    Sum of: ground roll + rotation + transition + climb

    Returns
    -------
    dict with V-speeds, segment distances, TODR, and FAR-factored TODR
    """
    V_s = V_stall(W, S, rho, CL_max_TO)
    V_R = 1.1 * V_s      # rotation speed
    V_2 = 1.2 * V_s      # takeoff safety speed

    S_G = ground_roll_distance(W, S, T, CD0, CL_ground, K, mu, rho, CL_max_TO)
    S_R = rotation_distance(V_R, t_rotate)

    trans = transition_segment(W, S, rho, CL_max_TO, TW)
    h_TR = trans["h_TR"]
    gamma = trans["gamma_rad"]
    R = trans["R"]

    if h_TR >= h_obstacle:
        # Obstacle cleared during transition — truncate arc at h_obstacle
        theta = np.arccos(np.clip(1.0 - h_obstacle / R, -1, 1))
        S_TR = R * np.sin(theta)
        S_C = 0.0
    else:
        S_TR = trans["S_TR"]
        S_C = climb_distance(h_obstacle, h_TR, gamma) if gamma > 0.01 else 0.0

    TODR = S_G + S_R + S_TR + S_C
    TODR_factored = 1.15 * TODR   # FAR 25.113(a): 115% of all-engine distance

    return {
        "V_stall_TO": V_s,
        "V_R": V_R,
        "V_2": V_2,
        "S_ground_roll": S_G,
        "S_rotation": S_R,
        "S_transition": S_TR,
        "h_transition": h_TR,
        "S_climb": S_C,
        "gamma_deg": trans["gamma_deg"],
        "TODR": TODR,
        "TODR_factored": TODR_factored,
    }


# =========================================================================
# Accelerate-Stop Distance (ASDR)
# =========================================================================

def accelerate_stop_distance(W, S, T, CD0, CL_ground, K, mu_roll, mu_brake,
                             rho, CL_max_TO, V_EF,
                             T_idle=0.0, T_reverse=0.0, t_react=2.0):
    """Accelerate-stop distance for a given engine-failure speed V_EF [ft].

    ASD = distance to accelerate from 0 to V_EF
        + reaction distance (V_EF * t_react)
        + braking distance from V_EF to 0

    Parameters
    ----------
    V_EF      : float or array, engine-failure speed [ft/s]
    t_react   : float, pilot reaction time [s] (FAR 25 allows ~2 s)

    Returns
    -------
    ASDR : float or array, accelerate-stop distance [ft]
    """
    S_accel = _ground_roll_to_V(W, S, T, CD0, CL_ground, K, mu_roll, rho,
                                0.0, V_EF)
    S_react = V_EF * t_react
    S_brake = _braking_distance_from_V(W, S, V_EF, mu_brake, rho,
                                        CD0, CL_ground, K, T_idle, T_reverse)
    return S_accel + S_react + S_brake


# =========================================================================
# Accelerate-Go Distance (OEI — one engine inoperative)
# =========================================================================

def accelerate_go_distance(W, S, T, CD0, CL_ground, K, mu_roll, rho,
                           CL_max_TO, TW, h_obstacle, V_EF,
                           n_engines=2, t_rotate=3.0):
    """Accelerate-go distance for a given engine-failure speed V_EF [ft].

    Distance = acceleration on all engines to V_EF
             + continued acceleration on (N-1) engines to V_R
             + rotation + transition + climb to obstacle (OEI)

    Parameters
    ----------
    V_EF      : float, engine-failure speed [ft/s]
    n_engines : int, total number of engines

    Returns
    -------
    AGDR : float, accelerate-go distance [ft]
    """
    V_s = V_stall(W, S, rho, CL_max_TO)
    V_R = 1.1 * V_s

    # Phase 1: all-engine acceleration from 0 to V_EF
    S1 = _ground_roll_to_V(W, S, T, CD0, CL_ground, K, mu_roll, rho,
                           0.0, V_EF)

    # Phase 2: OEI acceleration from V_EF to V_R
    T_OEI = T * (n_engines - 1) / n_engines
    CD0_OEI = CD0 + 0.005   # windmilling drag increment
    if V_EF < V_R:
        S2 = _ground_roll_to_V(W, S, T_OEI, CD0_OEI, CL_ground, K, mu_roll,
                               rho, V_EF, V_R)
    else:
        S2 = 0.0

    # Rotation
    S_R = rotation_distance(V_R, t_rotate)

    # Transition and climb (OEI)
    TW_OEI = T_OEI / W
    trans = transition_segment(W, S, rho, CL_max_TO, TW_OEI)
    h_TR = trans["h_TR"]
    gamma = trans["gamma_rad"]
    R = trans["R"]

    if h_TR >= h_obstacle:
        theta = np.arccos(np.clip(1.0 - h_obstacle / R, -1, 1))
        S_TR = R * np.sin(theta)
        S_C = 0.0
    else:
        S_TR = trans["S_TR"]
        S_C = climb_distance(h_obstacle, h_TR, gamma) if gamma > 0.01 else 0.0

    return S1 + S2 + S_R + S_TR + S_C


# =========================================================================
# V1 Decision Speed (iterative)
# =========================================================================

def find_V1(W, S, T, CD0, CL_ground, K, mu_roll, mu_brake, rho,
            CL_max_TO, TW, h_obstacle=35.0, n_engines=2,
            T_idle=0.0, T_reverse=0.0, t_react=2.0, t_rotate=3.0):
    """Find V1 decision speed where ASDR = accelerate-go distance.

    Iterates over candidate V_EF speeds between V_MCG (~0.7*V_stall)
    and V_R (1.1*V_stall) to find the intersection.

    Returns
    -------
    dict with V1, ASDR, accelerate-go distance, and BFL (the max of the two)
    """
    V_s = V_stall(W, S, rho, CL_max_TO)
    V_R = 1.1 * V_s

    V_lo = 0.7 * V_s   # approximate lower bound
    V_hi = V_R

    # Bisection: find V where ASDR = AGDR
    for _ in range(60):
        V_mid = 0.5 * (V_lo + V_hi)

        asd = accelerate_stop_distance(W, S, T, CD0, CL_ground, K,
                                       mu_roll, mu_brake, rho, CL_max_TO,
                                       V_mid, T_idle, T_reverse, t_react)
        agd = accelerate_go_distance(W, S, T, CD0, CL_ground, K,
                                     mu_roll, rho, CL_max_TO, TW,
                                     h_obstacle, V_mid, n_engines, t_rotate)

        if asd < agd:
            V_lo = V_mid
        else:
            V_hi = V_mid

        if abs(asd - agd) < 1.0:  # converged within 1 ft
            break

    V1 = 0.5 * (V_lo + V_hi)
    asd_final = accelerate_stop_distance(W, S, T, CD0, CL_ground, K,
                                         mu_roll, mu_brake, rho, CL_max_TO,
                                         V1, T_idle, T_reverse, t_react)
    agd_final = accelerate_go_distance(W, S, T, CD0, CL_ground, K,
                                       mu_roll, rho, CL_max_TO, TW,
                                       h_obstacle, V1, n_engines, t_rotate)
    BFL = max(asd_final, agd_final)

    return {
        "V1": V1,
        "ASDR_at_V1": asd_final,
        "AGDR_at_V1": agd_final,
        "BFL": BFL,
    }


# =========================================================================
# ASDR and TODR curves (for plotting)
# =========================================================================

def asdr_todr_curves(W, S, T, CD0, CL_ground, K, mu_roll, mu_brake, rho,
                     CL_max_TO, TW, h_obstacle=35.0, n_engines=2,
                     T_idle=0.0, T_reverse=0.0, t_react=2.0, t_rotate=3.0,
                     n_points=80):
    """Compute ASDR and accelerate-go distance vs engine-failure speed.

    Returns
    -------
    dict with arrays:
        V_EF  : engine-failure speeds [ft/s]
        ASDR  : accelerate-stop distance at each V_EF [ft]
        AGDR  : accelerate-go distance at each V_EF [ft]
    """
    V_s = V_stall(W, S, rho, CL_max_TO)
    V_R = 1.1 * V_s
    V_arr = np.linspace(0.5 * V_s, V_R, n_points)

    asd_arr = np.array([
        accelerate_stop_distance(W, S, T, CD0, CL_ground, K,
                                 mu_roll, mu_brake, rho, CL_max_TO,
                                 v, T_idle, T_reverse, t_react)
        for v in V_arr])
    agd_arr = np.array([
        accelerate_go_distance(W, S, T, CD0, CL_ground, K,
                               mu_roll, rho, CL_max_TO, TW,
                               h_obstacle, v, n_engines, t_rotate)
        for v in V_arr])

    return {"V_EF": V_arr, "ASDR": asd_arr, "AGDR": agd_arr}


# =========================================================================
# Balanced Field Length (empirical — Raymer Eq. 17.113)
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
