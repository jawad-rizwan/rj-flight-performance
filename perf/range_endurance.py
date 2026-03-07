"""
Range and Endurance calculations (Raymer 17.2.4 -- 17.2.9).

Covers Breguet range (jet & prop), loiter endurance,
and range/endurance optimization velocities.
"""

import numpy as np
from .level_flight import LD_max, CL_max_LD, V_min_thrust, V_min_power


# =========================================================================
# Range -- Jet
# =========================================================================

def breguet_range_jet(V, C, LD, Wi, Wf):
    """Breguet range equation for jet aircraft.

    Raymer Eq. (17.23):  R = (V/C) * (L/D) * ln(Wi/Wf)

    Parameters
    ----------
    V  : float, cruise velocity [ft/s]
    C  : float, TSFC [1/s]  (divide hourly TSFC by 3600)
    LD : float, cruise L/D
    Wi : float, initial weight [lb]
    Wf : float, final weight [lb]

    Returns
    -------
    R : float, range [ft]
    """
    return (V / C) * LD * np.log(Wi / Wf)   # Eq. (17.23)


def breguet_range_jet_nmi(V_kts, C_hr, LD, Wi, Wf):
    """Breguet range for jet, practical units.

    Parameters
    ----------
    V_kts : float, cruise velocity [knots]
    C_hr  : float, TSFC [1/hr]
    LD    : float, cruise L/D
    Wi    : float, initial weight [lb]
    Wf    : float, final weight [lb]

    Returns
    -------
    R : float, range [nautical miles]
    """
    return (V_kts / C_hr) * LD * np.log(Wi / Wf)


def V_best_range_jet(W, S, rho, CD0, K):
    """Velocity for maximum range -- jet aircraft.

    Raymer Eq. (17.25):  V = sqrt(2W/(rho*S)) * (3K/CD0)^0.25

    At best-range speed for a jet, CL = sqrt(CD0/(3K))  [Eq. (17.26)].
    The drag coefficient is CD0 + CD0/3 = (4/3)*CD0  [Eq. (17.27)].
    L/D at best range is ~86.6% of (L/D)_max.
    """
    return np.sqrt(2.0 * W / (rho * S)) * (3.0 * K / CD0)**0.25   # Eq. (17.25)


def CL_best_range_jet(CD0, K):
    """CL for best range -- jet.

    Raymer Eq. (17.26): CL = sqrt(CD0 / (3K))
    """
    return np.sqrt(CD0 / (3.0 * K))   # Eq. (17.26)


def D_best_range_jet(CD0, S, q):
    """Drag at best-range condition -- jet.

    Raymer Eq. (17.27): D = qS*(CD0 + CD0/3) = (4/3)*qS*CD0
    """
    return q * S * (CD0 + CD0 / 3.0)   # Eq. (17.27)


# =========================================================================
# Range -- Prop
# =========================================================================

def breguet_range_prop(eta_p, C_power, LD, Wi, Wf):
    """Breguet range equation for propeller aircraft.

    Raymer Eq. (17.28): R = (eta_p / C_power) * (L/D) * ln(Wi/Wf)

    For a prop, maximum range occurs at maximum L/D.

    Parameters
    ----------
    eta_p   : float, propeller efficiency
    C_power : float, power-specific fuel consumption [1/s]
              (C_power = C_bhp / 550  if C_bhp in lb/hr/bhp and want 1/s)
    LD      : float, L/D
    Wi, Wf  : float, initial / final weight [lb]

    Returns
    -------
    R : float, range [ft]
    """
    return (eta_p / C_power) * LD * np.log(Wi / Wf)   # Eq. (17.28)


def breguet_range_prop_nmi(eta_p, C_bhp_hr, LD, Wi, Wf):
    """Breguet range for prop, practical units.

    Raymer Eq. (17.28) rearranged:
        R [nmi] = (550 * eta_p / C_bhp) * (L/D) * ln(Wi/Wf) / 6076

    Parameters
    ----------
    eta_p    : float, propeller efficiency
    C_bhp_hr : float, BSFC [lb/hr/bhp]
    LD       : float, L/D
    Wi, Wf   : float, weights [lb]

    Returns
    -------
    R : float, range [nmi]
    """
    # Convert: 1 bhp = 550 ft*lb/s, 1 nmi = 6076 ft
    R_ft = (550.0 * eta_p / (C_bhp_hr / 3600.0)) * LD * np.log(Wi / Wf)
    return R_ft / 6076.0


# =========================================================================
# Endurance / Loiter -- Jet
# =========================================================================

def endurance_jet(C, LD, Wi, Wf):
    """Loiter endurance for jet aircraft.

    Raymer Eq. (17.30):  E = (1/C) * (L/D) * ln(Wi/Wf)

    Best endurance for a jet => maximize L/D
    => fly at V for max L/D [Eq. (17.13)].

    Parameters
    ----------
    C  : float, TSFC [1/s]
    LD : float, L/D
    Wi, Wf : float, weights [lb]

    Returns
    -------
    E : float, endurance [seconds]
    """
    return (1.0 / C) * LD * np.log(Wi / Wf)   # Eq. (17.30)


def endurance_jet_hr(C_hr, LD, Wi, Wf):
    """Loiter endurance for jet [hours].

    Parameters
    ----------
    C_hr : float, TSFC [1/hr]
    """
    return (1.0 / C_hr) * LD * np.log(Wi / Wf)


# =========================================================================
# Endurance / Loiter -- Prop
# =========================================================================

def endurance_prop(eta_p, C_power, LD, V, Wi, Wf):
    """Loiter endurance for propeller aircraft.

    Raymer Eq. (17.31):  E = (L/D) * (eta_p / (C_power * V)) * ln(Wi/Wf)

    Best endurance for prop => maximize (L/D) / V
    => fly at V for minimum power [Eq. (17.19)].

    Parameters
    ----------
    eta_p   : float, propeller efficiency
    C_power : float, BSFC [1/s]
    LD      : float, L/D at loiter condition
    V       : float, loiter velocity [ft/s]
    Wi, Wf  : float, weights [lb]

    Returns
    -------
    E : float, endurance [seconds]
    """
    return (1.0 / C_power) * (eta_p / V) * LD * np.log(Wi / Wf)


# =========================================================================
# Loiter-Cruise relationship
# =========================================================================

def loiter_from_cruise(R_cruise_nmi, V_cruise_kts):
    """Approximate loiter time from known cruise range and speed.

    Raymer Eq. (17.34):  E_loiter ~ 1.14 * (R_cruise / V_cruise)

    Returns
    -------
    E : float, approximate loiter endurance [hours]
    """
    return 1.14 * (R_cruise_nmi / V_cruise_kts)   # Eq. (17.34)


# =========================================================================
# Range parameter
# =========================================================================

def range_parameter_jet(V, C, CL, CD):
    """Range parameter for jet: (V/C)*(CL/CD).

    Raymer Eq. (17.24): (V/C)*(L/D) = (V/C)*(CL/(CD0 + K*CL^2))

    Maximising this gives the best-range condition.
    """
    return (V / C) * (CL / CD)   # Eq. (17.24)
