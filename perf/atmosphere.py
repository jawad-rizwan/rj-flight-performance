"""
International Standard Atmosphere (ISA) model.
Provides density, pressure, temperature, and speed of sound
as functions of altitude in feet.
"""

import numpy as np

# --- Sea-level constants (Imperial) ---
RHO_SL = 0.002377       # slug/ft^3
T_SL   = 518.67         # Rankine  (59 F)
P_SL   = 2116.22        # lb/ft^2
A_SL   = 1116.45        # ft/s, speed of sound
G      = 32.174         # ft/s^2
GAMMA  = 1.4            # ratio of specific heats for air
R_AIR  = 1716.49        # ft*lb/(slug*R), gas constant for air

# Troposphere lapse rate
LAPSE  = 0.003566       # R/ft  (about -6.5 K/km)
TROPO  = 36_089.0       # ft, tropopause altitude
T_TROPO = T_SL - LAPSE * TROPO  # ~389.97 R


def isa_temperature(h_ft):
    """ISA temperature [R] at altitude h_ft [ft]."""
    h = np.asarray(h_ft, dtype=float)
    T = np.where(h <= TROPO,
                 T_SL - LAPSE * h,
                 T_TROPO)
    return float(T) if T.ndim == 0 else T


def isa_pressure(h_ft):
    """ISA pressure [lb/ft^2] at altitude h_ft [ft]."""
    h = np.asarray(h_ft, dtype=float)
    T = isa_temperature(h)
    P_tropo = np.where(
        h <= TROPO,
        P_SL * (T / T_SL) ** (G / (LAPSE * R_AIR)),
        P_SL * (T_TROPO / T_SL) ** (G / (LAPSE * R_AIR))
    )
    P = np.where(
        h <= TROPO,
        P_tropo,
        P_tropo * np.exp(-G * (h - TROPO) / (R_AIR * T_TROPO))
    )
    return float(P) if P.ndim == 0 else P


def isa_density(h_ft, dT_C=0.0):
    """ISA density [slug/ft^3] at altitude h_ft [ft].

    Parameters
    ----------
    h_ft : float or array, altitude [ft]
    dT_C : float, ISA temperature deviation [deg C / K].
           e.g. dT_C=20 for ISA+20.  Pressure is unchanged;
           only temperature (and thus density) is affected.
    """
    P = isa_pressure(h_ft)
    T = isa_temperature(h_ft) + dT_C * 1.8   # convert C -> Rankine offset
    rho = P / (R_AIR * T)
    return float(rho) if np.ndim(rho) == 0 else rho


def speed_of_sound(h_ft):
    """Speed of sound [ft/s] at altitude h_ft [ft]."""
    T = isa_temperature(h_ft)
    a = np.sqrt(GAMMA * R_AIR * T)
    return float(a) if np.ndim(a) == 0 else a


def dynamic_pressure(V, h_ft):
    """Dynamic pressure q = 0.5 * rho * V^2  [lb/ft^2].

    Parameters
    ----------
    V    : float or array, true airspeed [ft/s]
    h_ft : float, altitude [ft]
    """
    rho = isa_density(h_ft)
    return 0.5 * rho * np.asarray(V, dtype=float)**2


def sigma(h_ft, dT_C=0.0):
    """Density ratio rho/rho_SL (ISA SL reference)."""
    return isa_density(h_ft, dT_C) / RHO_SL


def mach_number(V, h_ft):
    """Mach number for true airspeed V [ft/s] at altitude h_ft [ft]."""
    return V / speed_of_sound(h_ft)


def TAS_from_mach(M, h_ft):
    """True airspeed [ft/s] from Mach number at altitude h_ft [ft]."""
    return M * speed_of_sound(h_ft)


def kts_to_fps(V_kts):
    """Convert knots to ft/s."""
    return V_kts * 1.68781


def fps_to_kts(V_fps):
    """Convert ft/s to knots."""
    return V_fps / 1.68781


def thrust_at_altitude(T_SL, h_ft, BPR=5.0, dT_C=0.0):
    """Thrust lapse model for turbofans.

    T/T_SL ~ sigma^n with exponent scaling by BPR:
      - BPR <= 3  : n ≈ 0.7  (low-bypass, military-class)
      - BPR ~ 5   : n ≈ 0.75 (CF34-class)
      - BPR ~ 9   : n ≈ 0.82 (geared turbofan, PW1000G-class)
      - BPR >= 12  : n ≈ 0.88

    Above the tropopause an additional 0.85 factor accounts for
    the isothermal stratosphere where ram recovery drops off.

    Parameters
    ----------
    dT_C : float, ISA temperature deviation [deg C].
           Hot day reduces density and thus available thrust.
    """
    sig = sigma(h_ft, dT_C)
    h = np.asarray(h_ft, dtype=float)
    # Linear interpolation: n = 0.70 at BPR=3, 0.90 at BPR=13
    exp = np.clip(0.70 + (BPR - 3.0) * 0.02, 0.70, 0.90)

    # Below tropopause
    ratio = sig**exp

    # Above tropopause: additional penalty for constant-temperature region
    if np.ndim(h) == 0:
        if h > TROPO:
            ratio = ratio * 0.85  # stratospheric penalty
    else:
        mask = h > TROPO
        ratio = np.where(mask, ratio * 0.85, ratio)

    return T_SL * ratio
