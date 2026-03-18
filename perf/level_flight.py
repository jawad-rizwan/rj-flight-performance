"""
Steady Level Flight performance (Raymer 17.2).

Covers thrust required, power required, T/W, L/D,
minimum-thrust and minimum-power conditions.
"""

import numpy as np
from .atmosphere import isa_density, dynamic_pressure


# ---------------------------------------------------------------------------
# Thrust & drag in level flight
# ---------------------------------------------------------------------------

def thrust_required(W, CD0, K, S, rho, V):
    """Thrust required for steady level flight.

    Raymer Eq. (17.8): T = D = qS(CD0 + K*CL^2)
    with CL from Eq. (17.9): L = W = qSCL  =>  CL = W/(qS)

    Parameters
    ----------
    W   : float, aircraft weight [lb]
    CD0 : float, zero-lift drag coefficient
    K   : float, induced drag factor
    S   : float, wing reference area [ft^2]
    rho : float, air density [slug/ft^3]
    V   : float or array, true airspeed [ft/s]

    Returns
    -------
    T_req : float or array, thrust required [lb]
    """
    V = np.asarray(V, dtype=float)
    q = 0.5 * rho * V**2
    CL = W / (q * S)
    CD = CD0 + K * CL**2                  # Eq. (17.8) drag polar
    return q * S * CD


def power_required(W, CD0, K, S, rho, V):
    """Power required for steady level flight.

    Raymer Eq. (17.16): P = DV = (1/2)*rho*V^3*S*CD0 + K*W^2/(0.5*rho*V*S)

    Returns
    -------
    P_req : float or array [ft*lb/s]
    """
    T = thrust_required(W, CD0, K, S, rho, V)
    return T * np.asarray(V, dtype=float)


# ---------------------------------------------------------------------------
# T/W and L/D
# ---------------------------------------------------------------------------

def TW_level(CD0, K, W_over_S, q):
    """Thrust-to-weight ratio in level flight.

    Raymer Eq. (17.11): T/W = qCD0/(W/S) + (W/S)*K/q

    Parameters
    ----------
    CD0      : float
    K        : float
    W_over_S : float, wing loading [lb/ft^2]
    q        : float or array, dynamic pressure [lb/ft^2]
    """
    q = np.asarray(q, dtype=float)
    return q * CD0 / W_over_S + W_over_S * K / q


def LD_from_CL(CL, CD0, K):
    """Lift-to-drag ratio.  L/D = CL / (CD0 + K*CL^2)."""
    CL = np.asarray(CL, dtype=float)
    return CL / (CD0 + K * CL**2)


def LD_max(CD0, K):
    """Maximum L/D.

    Raymer Eq. (17.67): (L/D)_max = 1 / (2*sqrt(CD0*K))
    Occurs at CL = sqrt(CD0/K)  [Eq. (17.14)]
    """
    return 1.0 / (2.0 * np.sqrt(CD0 * K))


def CL_max_LD(CD0, K):
    """CL for maximum L/D (minimum drag).

    Raymer Eq. (17.14): CL = sqrt(CD0 / K)
    """
    return np.sqrt(CD0 / K)


# ---------------------------------------------------------------------------
# Minimum thrust (max L/D) velocity
# ---------------------------------------------------------------------------

def V_min_thrust(W, S, rho, CD0, K):
    """Velocity for minimum thrust required (maximum L/D).

    Raymer Eq. (17.13): V = sqrt(2W/(rho*S)) * (K/CD0)^0.25

    Parameters
    ----------
    W   : float, weight [lb]
    S   : float, wing area [ft^2]
    rho : float, density [slug/ft^3]
    CD0, K : float, drag polar parameters

    Returns
    -------
    V : float [ft/s]
    """
    return np.sqrt(2.0 * W / (rho * S)) * (K / CD0)**0.25   # Eq. (17.13)


def D_min_thrust(W, CD0, K):
    """Minimum drag (at best L/D speed).

    Raymer Eq. (17.15): D_min = qS*(CD0 + CD0) = 2*q*S*CD0
    At this condition induced drag equals zero-lift drag, so
    D_min = W / (L/D)_max.
    """
    return W / LD_max(CD0, K)


# ---------------------------------------------------------------------------
# Minimum power velocity
# ---------------------------------------------------------------------------

def V_min_power(W, S, rho, CD0, K):
    """Velocity for minimum power required.

    Raymer Eq. (17.19): V_min_power = sqrt(2W/(rho*S)) * (K/(3*CD0))^0.25

    This is about 0.76 times V_min_thrust.
    """
    return np.sqrt(2.0 * W / (rho * S)) * (K / (3.0 * CD0))**0.25  # Eq. (17.19)


def CL_min_power(CD0, K):
    """CL for minimum power.

    Raymer Eq. (17.20): CL = sqrt(3*CD0 / K)
    At this condition, induced drag = 3 * zero-lift drag.
    """
    return np.sqrt(3.0 * CD0 / K)   # Eq. (17.20)


def D_min_power(CD0, K, S, rho, W):
    """Drag at minimum-power condition.

    Raymer Eq. (17.21): D = qS*(CD0 + 3*CD0) = 4*q*S*CD0
    """
    V = V_min_power(W, S, rho, CD0, K)
    q = 0.5 * rho * V**2
    return q * S * (CD0 + 3.0 * CD0)   # Eq. (17.21)


# ---------------------------------------------------------------------------
# Level-flight velocity from wing loading
# ---------------------------------------------------------------------------

def V_level(W, S, rho, CL):
    """Velocity in level flight for a given CL.

    Raymer Eq. (17.10): V = sqrt(2/(rho*CL) * (W/S))
    """
    return np.sqrt(2.0 * W / (rho * CL * S))   # Eq. (17.10)


def V_stall(W, S, rho, CL_max):
    """Stall speed [ft/s].

    From Eq. (17.10) with CL = CL_max.
    """
    return V_level(W, S, rho, CL_max)


def V_max(W, S, rho, CD0, K, T_avail):
    """Maximum level-flight speed [ft/s].

    Raymer Sec. 17.2.3: the highest velocity where T_avail >= T_req.
    Solves  T = qS*CD0 + K*W^2/(qS)  for V via the quartic in V^2.

    Parameters
    ----------
    T_avail : float, thrust available at flight condition [lb]

    Returns
    -------
    V_max : float [ft/s], or NaN if no solution (thrust insufficient)
    """
    # T = (rho/2)*V^2*S*CD0 + K*W^2 / ((rho/2)*V^2*S)
    # Let x = q = 0.5*rho*V^2:  T = x*S*CD0 + K*W^2/(x*S)
    # x^2*(S*CD0) - x*T + K*W^2/S = 0
    a_coeff = S * CD0
    b_coeff = -T_avail
    c_coeff = K * W**2 / S
    disc = b_coeff**2 - 4.0 * a_coeff * c_coeff
    if disc < 0:
        return float('nan')
    q_max = (-b_coeff + np.sqrt(disc)) / (2.0 * a_coeff)  # larger root = higher V
    return np.sqrt(2.0 * q_max / rho)
