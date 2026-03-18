"""
Transonic wave drag estimation.
Compressibility drag rise model based on Raymer Chapter 12.

Provides:
    mdd_wing    — drag-divergence Mach from Korn equation
    cd0_at_mach — CD0 adjusted for compressibility drag rise
"""

import numpy as np


def mdd_wing(t_c, sweep_qc_deg, cl_design=0.0, supercritical=True):
    """Wing drag-divergence Mach number (Korn equation).

    Raymer Eq. 12.46 / Fig. 12.29 approximation.

    Parameters
    ----------
    t_c            : float, thickness-to-chord ratio
    sweep_qc_deg   : float, quarter-chord sweep [deg]
    cl_design      : float, design lift coefficient (default 0)
    supercritical  : bool, True for supercritical airfoils (kA=0.95)

    Returns
    -------
    float : M_DD (Boeing definition, ~20 counts drag rise)
    """
    cos_s = np.cos(np.radians(sweep_qc_deg))
    ka = 0.95 if supercritical else 0.87
    return ka / cos_s - t_c / cos_s**2 - cl_design / (10.0 * cos_s**3)


def cd0_at_mach(cd0_sub, mach, mdd):
    """CD0 with compressibility drag rise (Raymer Fig. 12.32 method).

    Piecewise model:
        M <= M_cr (M_DD - 0.08)  : cd0_sub (no change)
        M_cr < M <= M_DD         : smooth quadratic rise (+0.002 at M_DD)
        M > M_DD                 : steep linear rise (slope ~0.15/Mach)

    Parameters
    ----------
    cd0_sub : float, subsonic CD0
    mach    : float or array, Mach number
    mdd     : float, drag-divergence Mach number

    Returns
    -------
    float or array : total CD0 at given Mach
    """
    m = np.asarray(mach, dtype=float)
    mdd = np.asarray(mdd, dtype=float)
    m_cr = mdd - 0.08
    cd = np.full_like(m, cd0_sub)

    # Crest-critical to drag-divergence: quadratic rise to +20 counts
    mask1 = (m > m_cr) & (m <= mdd)
    denom = np.maximum(mdd - m_cr, 1e-6)
    frac1 = np.where(mask1, (m - m_cr) / denom, 0.0)
    cd = np.where(mask1, cd0_sub + 0.002 * frac1**2, cd)

    # Above M_DD: steep drag rise
    mask2 = m > mdd
    excess = np.where(mask2, m - mdd, 0.0)
    cd = np.where(mask2, cd0_sub + 0.002 + 0.15 * excess, cd)

    return float(cd) if cd.ndim == 0 else cd
