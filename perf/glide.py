"""
Gliding Flight (Raymer 17.5).

Covers straight glide, best glide ratio, minimum sink rate,
and turning glide.
"""

import numpy as np
from .atmosphere import G


# =========================================================================
# Straight glide
# =========================================================================

def glide_angle(LD):
    """Glide angle [rad].

    Raymer Eq. (17.64):  L/D = 1/tan(gamma) ~ 1/gamma  (small angles)
    => gamma = arctan(1/(L/D))
    """
    return np.arctan(1.0 / np.asarray(LD, dtype=float))   # Eq. (17.64)


def glide_range(h, LD):
    """Horizontal glide range from altitude h [ft].

    R = h * (L/D)    (from Eq. 17.64, small angle)
    """
    return h * LD


def V_best_glide(W, S, rho, CD0, K):
    """Velocity for best glide ratio (= max L/D velocity).

    Raymer Eq. (17.65):  V = sqrt(2W/(rho*S)) * (K/CD0)^0.25

    Same as V_min_thrust for powered flight.
    """
    return np.sqrt(2.0 * W / (rho * S)) * (K / CD0)**0.25   # Eq. (17.65)


def CL_best_glide(CD0, K):
    """CL for best glide ratio (= CL for max L/D).

    Raymer Eq. (17.66):  CL = sqrt(CD0 / K)
    """
    return np.sqrt(CD0 / K)   # Eq. (17.66)


def LD_max_glide(CD0, K):
    """Maximum glide L/D.

    Raymer Eq. (17.67):  (L/D)_max = (1/2) * sqrt(pi*A*e / CD0)
                                   = 1 / (2*sqrt(CD0*K))
    """
    return 1.0 / (2.0 * np.sqrt(CD0 * K))   # Eq. (17.67)


# =========================================================================
# Sink rate
# =========================================================================

def sink_rate(W, S, rho, CL, CD):
    """Sink rate (vertical velocity, positive downward) [ft/s].

    Raymer Eq. (17.68):
        Vv = V * sin(gamma) = sqrt( (W/S) * (2/(rho*CL)) ) * (CD/CL)

    Rewritten: Vv = sqrt(2W/(rho*S)) * CD / CL^(3/2)
    """
    return np.sqrt(2.0 * W / (rho * S)) * CD / CL**1.5   # Eq. (17.68)


def V_min_sink(W, S, rho, CD0, K):
    """Velocity for minimum sink rate.

    Raymer Eq. (17.73):  V = sqrt(2W/(rho*S)) * (K/(3*CD0))^0.25

    Same as V_min_power for powered flight.
    """
    return np.sqrt(2.0 * W / (rho * S)) * (K / (3.0 * CD0))**0.25  # Eq. (17.73)


def CL_min_sink(CD0, K):
    """CL for minimum sink rate.

    Raymer Eq. (17.72):  CL = sqrt(3*CD0 / K)
    """
    return np.sqrt(3.0 * CD0 / K)   # Eq. (17.72)


def LD_min_sink(CD0, K):
    """L/D at minimum sink rate condition.

    Raymer Eq. (17.74):
        (L/D)_min_sink = sqrt(3/(16*K*CD0)) = sqrt(3*pi*A*e / (16*CD0))
    """
    return np.sqrt(3.0 / (16.0 * K * CD0))   # Eq. (17.74)


def min_sink_rate(W, S, rho, CD0, K):
    """Minimum achievable sink rate [ft/s].

    Evaluated at CL_min_sink, V_min_sink.
    """
    CL = CL_min_sink(CD0, K)
    CD = CD0 + K * CL**2        # = CD0 + 3*CD0 = 4*CD0
    return sink_rate(W, S, rho, CL, CD)


# =========================================================================
# Turning glide (Raymer 17.5.2)
# =========================================================================

def sink_rate_turning(W, S, rho, CL, CD, phi_rad):
    """Sink rate in a turning glide [ft/s].

    Raymer Eq. (17.80):
        Vv_turn = (1/cos^(3/2)(phi)) * Vv_straight

    The bank angle phi increases the effective wing loading.
    """
    Vv_straight = sink_rate(W, S, rho, CL, CD)
    return Vv_straight / np.cos(phi_rad)**1.5   # Eq. (17.80)


def turn_radius_glide(W, S, rho, CL, phi_rad):
    """Turn radius in a gliding turn [ft].

    Raymer Eq. (17.81):  R = 2W / (rho * S * CL * g * sin(phi))
    """
    return 2.0 * W / (rho * S * CL * G * np.sin(phi_rad))   # Eq. (17.81)
