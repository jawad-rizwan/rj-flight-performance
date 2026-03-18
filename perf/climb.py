"""
Steady Climbing and Descending Flight (Raymer 17.3).

Covers climb angle, rate of climb, best angle / rate for jet & prop,
time to climb, and fuel to climb.
"""

import numpy as np
from .atmosphere import isa_density, G


# =========================================================================
# Basic climb equations
# =========================================================================

def climb_angle(T, D, W):
    """Steady climb angle [rad].

    Raymer Eq. (17.38):  sin(gamma) = (T - D) / W
    For small angles: gamma ~ (T - D) / W

    Also Eq. (17.38): gamma = arcsin[(T/W) - 1/(L/D)]
    """
    sin_gamma = (T - D) / W
    sin_gamma = np.clip(sin_gamma, -1.0, 1.0)
    return np.arcsin(sin_gamma)   # Eq. (17.38)


def rate_of_climb(V, T, D, W):
    """Rate of climb (vertical velocity) [ft/s].

    Raymer Eq. (17.39):  Vv = V * sin(gamma) = V * (T - D) / W

    Parameters
    ----------
    V : float or array, true airspeed [ft/s]
    T : float or array, thrust [lb]
    D : float or array, drag [lb]
    W : float, weight [lb]

    Returns
    -------
    Vv : float or array, rate of climb [ft/s]
    """
    return np.asarray(V, dtype=float) * (np.asarray(T) - np.asarray(D)) / W  # Eq. (17.39)


def climb_velocity(W, S, rho, CL, gamma):
    """Velocity in steady climb.

    Raymer Eq. (17.40):  V = sqrt( (2/rho*CL) * (W/S) * cos(gamma) )
    """
    return np.sqrt(2.0 * W * np.cos(gamma) / (rho * CL * S))   # Eq. (17.40)


def TW_for_climb(gamma, LD):
    """T/W required for a steady climb at angle gamma.

    Raymer Eq. (17.41):  T/W = 1/(L/D) + sin(gamma)
    For small gamma:      T/W ~ 1/(L/D) + Vv/V
    """
    return 1.0 / LD + np.sin(gamma)   # Eq. (17.41)


# =========================================================================
# Best climb -- Jet
# =========================================================================

def V_best_ROC_jet(W, S, rho, CD0, K, TW):
    """Velocity for best rate of climb -- jet aircraft.

    Raymer Eq. (17.43):
        V = sqrt( (W/S)/(3*rho*CD0) * [T/W + sqrt((T/W)^2 + 12*CD0*K)] )

    Parameters
    ----------
    TW  : float, thrust-to-weight ratio at flight condition
    """
    term = TW + np.sqrt(TW**2 + 12.0 * CD0 * K)
    return np.sqrt(W / (S * 3.0 * rho * CD0) * term)   # Eq. (17.43)


def ROC_jet(V, T, CD0, K, W, S, rho):
    """Rate of climb for a jet at given V and T.

    Uses Eq. (17.39) with drag expanded:
        D = qS*CD0 + K*W^2/(qS)  (level-flight approx for small gamma)
    """
    q = 0.5 * rho * V**2
    D = q * S * CD0 + K * W**2 / (q * S)
    return V * (T - D) / W   # Eq. (17.39)


# =========================================================================
# Best climb -- Prop
# =========================================================================

def ROC_prop(P_avail, W, S, rho, CD0, K):
    """Best rate of climb for propeller aircraft.

    Raymer Eq. (17.45):
        Vv = P*eta_p/W  - D*V/W

    Best ROC for props occurs at velocity for minimum power required.
    Vv_max = (P_avail - P_min_required) / W

    Parameters
    ----------
    P_avail : float, available power = P_shaft * eta_p [ft*lb/s]
    """
    # Minimum power required at V_min_power
    V_mp = np.sqrt(2.0 * W / (rho * S)) * (K / (3.0 * CD0))**0.25  # Eq. (17.19)
    q = 0.5 * rho * V_mp**2
    CL = W / (q * S)
    CD = CD0 + K * CL**2
    D = q * S * CD
    P_req = D * V_mp
    return (P_avail - P_req) / W   # Eq. (17.45) rearranged


# =========================================================================
# Time to climb and fuel to climb
# =========================================================================

def time_to_climb(h1, h2, Vv1, Vv2):
    """Time to climb between two altitudes using linear interpolation.

    Raymer Eqs. (17.48)-(17.50):
        Vv = Vv_i - a*(h - h_i)
        a  = (Vv2 - Vv1) / (h2 - h1)
        dt = (1/a) * ln(Vv_i / Vv_{i+1})

    Parameters
    ----------
    h1, h2   : float, altitudes [ft]
    Vv1, Vv2 : float, rates of climb at h1, h2 [ft/s]

    Returns
    -------
    dt : float, time to climb [seconds]
    """
    if abs(Vv1 - Vv2) < 1e-6:
        # Constant ROC
        return (h2 - h1) / Vv1

    a = (Vv2 - Vv1) / (h2 - h1)                    # Eq. (17.49)
    # Eq. (17.50): dt = (1/a)*ln(Vv_i / Vv_{i+1})
    # When ROC decreases with altitude, a < 0 and ln > 0, so take abs.
    return abs((1.0 / a) * np.log(Vv1 / Vv2))       # Eq. (17.50)


def fuel_to_climb(C, T_avg, dt):
    """Fuel burned during climb segment.

    Raymer Eq. (17.51):  dW_fuel = C * T_avg * dt

    Parameters
    ----------
    C     : float, TSFC [1/s]
    T_avg : float, average thrust during segment [lb]
    dt    : float, time for segment [seconds]

    Returns
    -------
    dW_fuel : float [lb]
    """
    return C * T_avg * dt   # Eq. (17.51)


# =========================================================================
# FAR 25 Climb Gradient Check (Table F.4)
# =========================================================================

# Minimum required gradients for 2-engine aircraft (Table F.4)
FAR25_GRADIENTS_2ENG = {
    "1st_seg":     0.000,   # >= 0.0%    (LOF, TO flaps, gear down, OEI)
    "2nd_seg":     0.024,   # >= 2.4%    (V2, TO flaps, gear up, OEI)
    "4th_seg":     0.012,   # >= 1.2%    (>=1.25Vs clean, clean, gear up, OEI)
    "GA_approach": 0.021,   # >= 2.1%    (<=1.4Vs, approach flaps, gear up, OEI)
    "GA_landing":  0.032,   # >= 3.2%    (<=1.23Vs, landing flaps, gear down, AEO)
}

# Typical drag increments (consistent with Raymer methodology)
_DCD0_FLAP_TO   = 0.010    # TO flap drag increment
_DCD0_FLAP_LAND = 0.020    # landing flap drag increment
_DCD0_FLAP_APP  = 0.013    # approach flap drag (between TO and land)
_DCD0_GEAR      = 0.015    # landing-gear drag increment
_DCD0_WINDMILL  = 0.005    # windmilling-engine drag increment


def far25_climb_gradients(W, S, CD0, K, CL_max_TO, CL_max_L, CL_max_clean,
                          T_max_SL, n_engines, rho):
    """Compute FAR 25 climb gradients for all Table F.4 segments.

    Returns dict of dicts, each with keys: gradient, required, speed_kts, pass.

    Parameters
    ----------
    W           : float, takeoff weight [lb]
    S           : float, wing area [ft^2]
    CD0         : float, clean zero-lift drag
    K           : float, induced drag factor
    CL_max_TO   : float, max CL in TO config
    CL_max_L    : float, max CL in landing config
    CL_max_clean: float, max CL in clean config
    T_max_SL    : float, total installed SL thrust [lb]
    n_engines   : int, number of engines
    rho         : float, air density [slug/ft^3]
    """
    from .level_flight import V_stall

    results = {}
    T_OEI = T_max_SL * (n_engines - 1) / n_engines

    def _gradient(V, T, CD0_eff):
        """Climb gradient = (T - D) / W."""
        q = 0.5 * rho * V**2
        CL = W / (q * S)
        CD = CD0_eff + K * CL**2
        D = q * S * CD
        return (T - D) / W

    # 1st segment: LOF speed (~1.1Vs TO), TO flaps, gear DOWN, OEI
    Vs_to = V_stall(W, S, rho, CL_max_TO)
    V_lof = 1.1 * Vs_to
    cd0_1 = CD0 + _DCD0_FLAP_TO + _DCD0_GEAR + _DCD0_WINDMILL
    g1 = _gradient(V_lof, T_OEI, cd0_1)
    results["1st_seg"] = {
        "label": "1st segment (LOF, TO flaps, gear dn, OEI)",
        "gradient": g1, "required": FAR25_GRADIENTS_2ENG["1st_seg"],
        "speed_kts": V_lof / 1.6878, "pass": g1 >= FAR25_GRADIENTS_2ENG["1st_seg"],
    }

    # 2nd segment: V2 = 1.2Vs, TO flaps, gear UP, OEI
    V2 = 1.2 * Vs_to
    cd0_2 = CD0 + _DCD0_FLAP_TO + _DCD0_WINDMILL
    g2 = _gradient(V2, T_OEI, cd0_2)
    results["2nd_seg"] = {
        "label": "2nd segment (V2, TO flaps, gear up, OEI)",
        "gradient": g2, "required": FAR25_GRADIENTS_2ENG["2nd_seg"],
        "speed_kts": V2 / 1.6878, "pass": g2 >= FAR25_GRADIENTS_2ENG["2nd_seg"],
    }

    # 4th segment: >= 1.25Vs clean, clean config, gear UP, OEI
    Vs_clean = V_stall(W, S, rho, CL_max_clean)
    V_4th = 1.25 * Vs_clean
    cd0_4 = CD0 + _DCD0_WINDMILL
    g4 = _gradient(V_4th, T_OEI, cd0_4)
    results["4th_seg"] = {
        "label": "4th segment (1.25Vs, clean, gear up, OEI)",
        "gradient": g4, "required": FAR25_GRADIENTS_2ENG["4th_seg"],
        "speed_kts": V_4th / 1.6878, "pass": g4 >= FAR25_GRADIENTS_2ENG["4th_seg"],
    }

    # Go-around approach: <= 1.4Vs_ref, approach flaps, gear UP, OEI
    W_land = 0.85 * W
    Vs_land = V_stall(W_land, S, rho, CL_max_L)
    V_ga_app = 1.4 * Vs_land
    cd0_ga_app = CD0 + _DCD0_FLAP_APP + _DCD0_WINDMILL
    T_OEI_ga = T_OEI  # OEI
    g_ga_app = _gradient_w(V_ga_app, T_OEI_ga, cd0_ga_app, W_land, S, K, rho)
    results["GA_approach"] = {
        "label": "Go-around approach (1.4Vs, app flaps, gear up, OEI)",
        "gradient": g_ga_app, "required": FAR25_GRADIENTS_2ENG["GA_approach"],
        "speed_kts": V_ga_app / 1.6878, "pass": g_ga_app >= FAR25_GRADIENTS_2ENG["GA_approach"],
    }

    # Go-around landing: <= 1.23Vs_ref, landing flaps, gear DOWN, AEO
    V_ga_land = 1.23 * Vs_land
    cd0_ga_land = CD0 + _DCD0_FLAP_LAND + _DCD0_GEAR
    g_ga_land = _gradient_w(V_ga_land, T_max_SL, cd0_ga_land, W_land, S, K, rho)
    results["GA_landing"] = {
        "label": "Go-around landing (1.23Vs, ldg flaps, gear dn, AEO)",
        "gradient": g_ga_land, "required": FAR25_GRADIENTS_2ENG["GA_landing"],
        "speed_kts": V_ga_land / 1.6878, "pass": g_ga_land >= FAR25_GRADIENTS_2ENG["GA_landing"],
    }

    return results


def _gradient_w(V, T, CD0_eff, W, S, K, rho):
    """Climb gradient helper using explicit weight (for landing-weight segments)."""
    q = 0.5 * rho * V**2
    CL = W / (q * S)
    CD = CD0_eff + K * CL**2
    D = q * S * CD
    return (T - D) / W


def time_to_climb_profile(altitudes, ROCs):
    """Integrate time to climb over multiple altitude segments.

    Parameters
    ----------
    altitudes : array, altitude breakpoints [ft]  (ascending)
    ROCs      : array, rate of climb at each altitude [ft/s]

    Returns
    -------
    times : array, cumulative time [s] at each altitude
    """
    n = len(altitudes)
    times = np.zeros(n)
    for i in range(n - 1):
        times[i + 1] = times[i] + time_to_climb(
            altitudes[i], altitudes[i + 1], ROCs[i], ROCs[i + 1])
    return times
