"""
Level Turning Flight (Raymer 17.4).

Covers turn rate, turn radius, load factor, instantaneous
and sustained turn performance, and corner speed.
"""

import numpy as np
from .atmosphere import G


# =========================================================================
# Basic turn equations
# =========================================================================

def turn_rate(V, n):
    """Level turn rate [rad/s].

    Raymer Eq. (17.52):  psi_dot = g*sqrt(n^2 - 1) / V

    Parameters
    ----------
    V : float or array, true airspeed [ft/s]
    n : float or array, load factor

    Returns
    -------
    psi_dot : float or array [rad/s]
    """
    V = np.asarray(V, dtype=float)
    n = np.asarray(n, dtype=float)
    return G * np.sqrt(n**2 - 1.0) / V   # Eq. (17.52)


def turn_rate_deg(V, n):
    """Level turn rate [deg/s]."""
    return np.degrees(turn_rate(V, n))


def turn_radius(V, n):
    """Level turn radius [ft].

    Raymer Eq. (17.79):  R = V^2 / (g*sqrt(n^2 - 1))
    """
    V = np.asarray(V, dtype=float)
    n = np.asarray(n, dtype=float)
    return V**2 / (G * np.sqrt(n**2 - 1.0))   # Eq. (17.79)


def bank_angle(n):
    """Bank angle [rad] for a level coordinated turn.

    cos(phi) = 1/n  =>  phi = arccos(1/n)
    """
    return np.arccos(1.0 / np.asarray(n, dtype=float))


def load_factor_from_bank(phi_rad):
    """Load factor from bank angle.  n = 1/cos(phi)."""
    return 1.0 / np.cos(phi_rad)


# =========================================================================
# Sustained turn
# =========================================================================

def n_sustained(TW, q, W_over_S, CD0, K):
    """Maximum sustained load factor.

    Raymer Eq. (17.54):
        n = sqrt( (q / (K*(W/S))) * (T/W - q*CD0/(W/S)) )

    In a sustained turn, thrust = drag and altitude is maintained.

    Parameters
    ----------
    TW       : float, thrust-to-weight ratio at flight condition
    q        : float or array, dynamic pressure [lb/ft^2]
    W_over_S : float, wing loading [lb/ft^2]
    CD0      : float, zero-lift drag coefficient
    K        : float, induced drag factor
    """
    q = np.asarray(q, dtype=float)
    inner = (q / (K * W_over_S)) * (TW - q * CD0 / W_over_S)
    inner = np.maximum(inner, 0.0)
    return np.sqrt(inner)   # Eq. (17.54)


def n_sustained_simple(TW, LD):
    """Sustained load factor (simplified).

    Raymer Eq. (17.53):  n = (T/W) * (L/D)
    """
    return TW * LD   # Eq. (17.53)


# =========================================================================
# Instantaneous turn
# =========================================================================

def n_instantaneous_stall(CL_max, q, W_over_S):
    """Instantaneous load factor limited by stall (CL_max).

    n = q * CL_max / (W/S)
    """
    return q * CL_max / W_over_S


def n_instantaneous_structural(n_max):
    """Structural limit on load factor."""
    return n_max


def corner_speed(W_over_S, rho, CL_max, n_max):
    """Corner speed -- intersection of stall and structural limits.

    At corner speed, both CL_max and n_max are reached simultaneously.
    V_corner = sqrt(2*n_max*(W/S) / (rho*CL_max))

    Parameters
    ----------
    W_over_S : float, wing loading [lb/ft^2]
    rho      : float, air density [slug/ft^3]
    CL_max   : float, maximum lift coefficient
    n_max    : float, structural load factor limit

    Returns
    -------
    V_corner : float [ft/s]
    """
    return np.sqrt(2.0 * n_max * W_over_S / (rho * CL_max))


# =========================================================================
# Sustained turn rate envelope
# =========================================================================

def CL_sustained_turn(CD0, K):
    """CL for maximum sustained load factor (= CL for max L/D).

    Raymer Eq. (17.55):  L = nW = qS*sqrt(CD0/K)
    """
    return np.sqrt(CD0 / K)   # Eq. (17.55)


def sustained_turn_envelope(V_array, T_avail, W, S, CD0, K):
    """Compute sustained turn rate vs velocity.

    Parameters
    ----------
    V_array : array, velocities [ft/s]
    T_avail : float, available thrust [lb]
    W       : float, weight [lb]
    S       : float, wing area [ft^2]
    CD0, K  : float, drag polar

    Returns
    -------
    psi_dot : array, sustained turn rate [deg/s]
    n_sust  : array, sustained load factor
    """
    V = np.asarray(V_array, dtype=float)
    q = 0.5 * 0.002377 * V**2  # placeholder; caller should provide rho
    TW = T_avail / W
    WS = W / S

    n_arr = n_sustained(TW, q, WS, CD0, K)
    n_arr = np.maximum(n_arr, 1.0)
    psi = turn_rate_deg(V, n_arr)
    return psi, n_arr


def sustained_turn_envelope_rho(V_array, T_avail, W, S, CD0, K, rho):
    """Compute sustained turn rate vs velocity at a given density.

    Parameters
    ----------
    V_array : array, velocities [ft/s]
    T_avail : float, available thrust [lb]
    W       : float, weight [lb]
    S       : float, wing area [ft^2]
    CD0, K  : float, drag polar
    rho     : float, air density [slug/ft^3]

    Returns
    -------
    psi_dot : array, sustained turn rate [deg/s]
    n_sust  : array, sustained load factor
    """
    V = np.asarray(V_array, dtype=float)
    q = 0.5 * rho * V**2
    TW = T_avail / W
    WS = W / S

    n_arr = n_sustained(TW, q, WS, CD0, K)
    n_arr = np.maximum(n_arr, 1.0)
    psi = turn_rate_deg(V, n_arr)
    return psi, n_arr
