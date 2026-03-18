"""
Energy-Maneuverability Methods (Raymer 17.6).

Covers specific energy, specific excess power (Ps),
energy height, minimum time-to-climb trajectory,
and mission-segment weight fractions.
"""

import numpy as np
from .atmosphere import G, isa_density


# =========================================================================
# Energy equations
# =========================================================================

def total_energy(W, h, V):
    """Total energy of an aircraft [ft*lb].

    Raymer Eq. (17.84):  E = W*h + (1/2)*(W/g)*V^2
    """
    return W * h + 0.5 * (W / G) * np.asarray(V, dtype=float)**2  # Eq. (17.84)


def specific_energy(h, V):
    """Specific energy (energy height) h_e [ft].

    Raymer Eq. (17.85):  h_e = h + V^2 / (2g)

    Dividing total energy by weight gives units of distance.
    """
    return h + np.asarray(V, dtype=float)**2 / (2.0 * G)   # Eq. (17.85)


def Ps(T, D, W, V):
    """Specific excess power [ft/s].

    Raymer Eq. (17.88):  Ps = V*(T - D)/W = dh/dt + (V/g)*(dV/dt)

    At constant velocity (dV/dt = 0), Ps equals the rate of climb.
    For load-factor-dependent Ps, use Ps_expanded() with Eq. (17.89).

    Parameters
    ----------
    T : float or array, thrust [lb]
    D : float or array, drag [lb]
    W : float, weight [lb]
    V : float or array, velocity [ft/s]

    Returns
    -------
    Ps : float or array [ft/s]
    """
    V = np.asarray(V, dtype=float)
    T = np.asarray(T, dtype=float)
    D = np.asarray(D, dtype=float)
    return V * (T - D) / W   # Eq. (17.88) with dV/dt=0


def Ps_expanded(V, TW, W_over_S, CD0, K, rho, n=1.0):
    """Specific excess power using aerodynamic coefficients.

    Raymer Eq. (17.89):
        Ps = V * [T/W - (rho*V^2/(2*(W/S)))*CD0 - n^2*K*(W/S)/(0.5*rho*V^2)]

    Parameters
    ----------
    V        : float or array, velocity [ft/s]
    TW       : float, T/W at flight condition
    W_over_S : float, wing loading [lb/ft^2]
    CD0      : float
    K        : float
    rho      : float, density [slug/ft^3]
    n        : float, load factor
    """
    V = np.asarray(V, dtype=float)
    q = 0.5 * rho * V**2
    return V * (TW - q * CD0 / W_over_S - n**2 * K * W_over_S / q)  # Eq. (17.89)


# =========================================================================
# Time to change energy height
# =========================================================================

def time_energy_height(dh_e, Ps_avg):
    """Approximate time to change energy height.

    Raymer Eq. (17.93):  dt ~ dh_e / Ps_avg

    Parameters
    ----------
    dh_e   : float, change in energy height [ft]
    Ps_avg : float, average Ps over the segment [ft/s]
    """
    return dh_e / Ps_avg   # Eq. (17.93)


def time_to_climb_energy(he_array, Ps_array):
    """Integrate time to climb using energy-height method.

    Raymer Eq. (17.92):  t = integral(1/Ps) dh_e

    Uses trapezoidal integration.

    Parameters
    ----------
    he_array : array, energy heights [ft]
    Ps_array : array, Ps at each energy height [ft/s]

    Returns
    -------
    t_total : float, total time [s]
    """
    inv_Ps = 1.0 / np.asarray(Ps_array, dtype=float)
    return np.trapz(inv_Ps, np.asarray(he_array, dtype=float))  # Eq. (17.92)


# =========================================================================
# Fuel-specific energy
# =========================================================================

def fuel_specific_energy(Ps, C, T):
    """Fuel-specific energy f_s [ft].

    Raymer Eq. (17.94):  f_s = Ps / (C*T)

    Minimise fuel to climb by maximising f_s.

    Parameters
    ----------
    Ps : float or array [ft/s]
    C  : float, TSFC [1/s]
    T  : float or array, thrust [lb]
    """
    return Ps / (C * T)   # Eq. (17.94)


def fuel_to_climb_energy(dh_e, fs_avg):
    """Approximate fuel weight to change energy height.

    Raymer Eq. (17.96):  dW_f ~ dh_e / f_s_avg

    Parameters
    ----------
    dh_e   : float, change in energy height [ft]
    fs_avg : float, average fuel-specific energy [ft]
    """
    return dh_e / fs_avg   # Eq. (17.96)


# =========================================================================
# Mission-segment weight fraction
# =========================================================================

def weight_fraction_energy(C, dh_e, V, LD, TW):
    """Mission-segment weight fraction for energy maneuver.

    Raymer Eq. (17.97):
        W_i / W_{i-1} = exp[ -C*dh_e / (V * {1 - 1/(TW * LD)}) ]

    Parameters
    ----------
    C    : float, TSFC [1/s]
    dh_e : float, change in energy height [ft]
    V    : float, velocity [ft/s]
    LD   : float, L/D
    TW   : float, T/W

    Returns
    -------
    Wi_over_Wim1 : float, weight fraction
    """
    denom = V * (1.0 - 1.0 / (TW * LD))
    if abs(denom) < 1e-10:
        return 1.0
    return np.exp(-C * dh_e / denom)   # Eq. (17.97)
