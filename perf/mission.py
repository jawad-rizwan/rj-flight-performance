"""
Mission Profile Analysis (Raymer 17.6.5, Breguet).

9-segment FAR 121.645 international mission profile with
segment-by-segment weight fractions and fuel accounting.

Segments
--------
0. Warmup & Takeoff         (historical)
1. Climb to cruise          (Eq. 17.97 energy method)
2. Cruise                   (Breguet range, Eq. 17.23)
3. Contingency loiter       (10% cruise time, Eq. 17.30)
4. Attempted landing        (historical)
--- trip fuel boundary ---
5. Go-around climb          (Eq. 17.97)
6. Divert to alternate      (Breguet)
7. Regulatory hold 30 min   (Eq. 17.30)
8. Land                     (historical)
"""

import numpy as np
from .atmosphere import (G, isa_density, RHO_SL, speed_of_sound,
                         TAS_from_mach, kts_to_fps, fps_to_kts,
                         thrust_at_altitude)
from .energy import specific_energy, weight_fraction_energy
from .level_flight import LD_max, V_stall, thrust_required
from .climb import ROC_jet, V_best_ROC_jet


def mission_profile(ac, contingency_pct=0.10, hold_time_hr=0.5):
    """Compute weight fractions for a 9-segment FAR 121.645 mission.

    Parameters
    ----------
    ac : AircraftData
    contingency_pct : float, contingency as fraction of cruise time
    hold_time_hr : float, regulatory hold time [hr]

    Returns
    -------
    dict with 'segments' list, fuel summary, and compliance metrics
    """
    W      = ac.W_TO
    S      = ac.S
    CD0    = ac.CD0
    K      = ac.K
    C_hr   = ac.TSFC
    C_s    = C_hr / 3600.0
    C_loiter_s = 0.80 * C_hr / 3600.0   # loiter TSFC ~80% cruise

    h_cr   = ac.h_cruise_ft
    rho_cr = isa_density(h_cr)
    V_cr   = TAS_from_mach(ac.M_cruise, h_cr)
    q_cr   = 0.5 * rho_cr * V_cr**2

    ld_max = LD_max(CD0, K)
    design_range = ac.design_range_nm
    alt_range    = ac.alternate_range_nm

    segments = []
    W_cur = W

    # ── 0  Warmup & Takeoff ──────────────────────────────────────
    wf = 0.970
    fuel = W_cur * (1.0 - wf)
    segments.append(_seg("0. Warmup & Takeoff", "Historical", wf, W_cur, fuel))
    W_cur -= fuel

    # ── 1  Climb to cruise ───────────────────────────────────────
    V_climb_init = kts_to_fps(250.0)
    he_0  = specific_energy(0, V_climb_init)
    he_cr = specific_energy(h_cr, V_cr)
    dh_e  = he_cr - he_0

    h_mid   = h_cr / 2.0
    rho_mid = isa_density(h_mid)
    T_mid   = thrust_at_altitude(ac.T_max_SL, h_mid, ac.BPR)
    TW_mid  = T_mid / W_cur
    V_mid   = 0.5 * (V_climb_init + V_cr)
    q_mid   = 0.5 * rho_mid * V_mid**2
    CL_mid  = W_cur / (q_mid * S)
    CD_mid  = CD0 + K * CL_mid**2
    LD_mid  = CL_mid / CD_mid

    wf = weight_fraction_energy(C_s, dh_e, V_mid, LD_mid, TW_mid)
    fuel = W_cur * (1.0 - wf)
    segments.append(_seg("1. Climb to FL350", "Eq. 17.97", wf, W_cur, fuel))
    W_cur -= fuel

    # Approximate climb distance (for subtracting from cruise range)
    avg_roc = 2500.0 / 60.0   # ~2500 fpm average in ft/s
    climb_time_s = h_cr / avg_roc
    climb_dist_nm = V_mid * climb_time_s / 6076.0

    # ── 2  Cruise ────────────────────────────────────────────────
    R_cruise_nm = max(design_range - climb_dist_nm, 0.0)
    R_cruise_ft = R_cruise_nm * 6076.0

    CL_cr = W_cur / (q_cr * S)
    CD_cr = CD0 + K * CL_cr**2
    LD_cr = CL_cr / CD_cr

    wf = np.exp(-R_cruise_ft * C_s / (V_cr * LD_cr))
    fuel = W_cur * (1.0 - wf)
    cruise_time_hr = R_cruise_nm / fps_to_kts(V_cr)
    segments.append(_seg(f"2. Cruise ({R_cruise_nm:.0f} nm)", "Breguet", wf, W_cur, fuel))
    W_cur -= fuel

    # ── 3  Contingency loiter (10% cruise time) ─────────────────
    loiter_s = contingency_pct * cruise_time_hr * 3600.0
    wf = np.exp(-loiter_s * C_loiter_s / ld_max)
    fuel = W_cur * (1.0 - wf)
    segments.append(_seg(f"3. Contingency ({contingency_pct*100:.0f}%)",
                         "Loiter", wf, W_cur, fuel))
    W_cur -= fuel

    # ── 4  Attempted landing ─────────────────────────────────────
    wf = 0.995
    fuel = W_cur * (1.0 - wf)
    segments.append(_seg("4. Attempted landing", "Historical", wf, W_cur, fuel))
    W_cur -= fuel

    # ── Trip fuel boundary ───────────────────────────────────────
    trip_fuel = W - W_cur

    # ── 5  Go-around climb ───────────────────────────────────────
    h_ga = 1500.0
    V_ga = kts_to_fps(180.0)
    dh_e_ga = specific_energy(h_ga, V_ga) - specific_energy(0, V_ga)
    T_ga  = ac.T_max_SL
    TW_ga = T_ga / W_cur
    q_ga  = 0.5 * RHO_SL * V_ga**2
    CL_ga = W_cur / (q_ga * S)
    CD_ga = CD0 + K * CL_ga**2
    LD_ga = CL_ga / CD_ga

    wf = weight_fraction_energy(C_s, dh_e_ga, V_ga, LD_ga, TW_ga)
    fuel = W_cur * (1.0 - wf)
    segments.append(_seg("5. Go-around climb", "Eq. 17.97", wf, W_cur, fuel))
    W_cur -= fuel

    # ── 6  Divert to alternate ───────────────────────────────────
    h_div   = 25_000.0
    rho_div = isa_density(h_div)
    V_div   = TAS_from_mach(ac.M_cruise * 0.90, h_div)
    q_div   = 0.5 * rho_div * V_div**2
    CL_div  = W_cur / (q_div * S)
    CD_div  = CD0 + K * CL_div**2
    LD_div  = CL_div / CD_div
    R_div_ft = alt_range * 6076.0

    wf = np.exp(-R_div_ft * C_s / (V_div * LD_div))
    fuel = W_cur * (1.0 - wf)
    segments.append(_seg(f"6. Divert ({alt_range:.0f} nm)", "Breguet", wf, W_cur, fuel))
    W_cur -= fuel

    # ── 7  Regulatory hold (30 min) ──────────────────────────────
    hold_s = hold_time_hr * 3600.0
    wf = np.exp(-hold_s * C_loiter_s / ld_max)
    fuel = W_cur * (1.0 - wf)
    segments.append(_seg(f"7. Hold ({hold_time_hr*60:.0f} min)", "Loiter", wf, W_cur, fuel))
    W_cur -= fuel

    # ── 8  Land ──────────────────────────────────────────────────
    wf = 0.995
    fuel = W_cur * (1.0 - wf)
    segments.append(_seg("8. Land", "Historical", wf, W_cur, fuel))
    W_cur -= fuel

    # ── Summary ──────────────────────────────────────────────────
    total_fuel   = W - W_cur
    reserve_fuel = total_fuel - trip_fuel
    wf_total     = W_cur / W

    # FAR 25 OEI climb gradients
    oei_gradients = _compute_oei_gradients(ac)

    # Time to climb (step estimate)
    ttc_min = _time_to_climb_estimate(ac)

    return {
        "segments":       segments,
        "W_TO":           W,
        "W_final":        W_cur,
        "wf_total":       wf_total,
        "total_fuel":     total_fuel,
        "trip_fuel":      trip_fuel,
        "reserve_fuel":   reserve_fuel,
        "fuel_available": ac.W_fuel_max,
        "design_range_nm": design_range,
        "cruise_range_nm": R_cruise_nm,
        "cruise_LD":      LD_cr,
        "climb_dist_nm":  climb_dist_nm,
        "oei":            oei_gradients,
        "ttc_min":        ttc_min,
    }


# ── Helpers ───────────────────────────────────────────────────

def _seg(name, method, wf, W_start, fuel):
    return {
        "name":    name,
        "method":  method,
        "wf":      float(wf),
        "W_start": float(W_start),
        "W_end":   float(W_start - fuel),
        "fuel":    float(fuel),
    }


def _compute_oei_gradients(ac):
    """FAR 25 OEI climb gradients (2-engine minimums).

    2nd segment : gear up, TO flaps, OEI, V2 = 1.2 Vs   -> >= 2.4%
    En-route    : clean config, OEI                       -> >= 1.2%
    Approach    : landing config, OEI, 1.3 Vs             -> >= 2.1%
    """
    W   = ac.W_TO
    S   = ac.S
    CD0 = ac.CD0
    K   = ac.K
    rho = RHO_SL

    T_oei = ac.T_max_SL / ac.n_engines   # one engine out

    # 2nd segment: TO flaps, gear up, V2 = 1.2*Vs_TO
    CD0_to = CD0 + 0.015   # flap drag, gear retracted
    Vs_to  = V_stall(W, S, rho, ac.CL_max_TO)
    V2     = 1.2 * Vs_to
    q2     = 0.5 * rho * V2**2
    CL2    = W / (q2 * S)
    CD2    = CD0_to + K * CL2**2
    D2     = q2 * S * CD2
    grad_2nd = (T_oei - D2) / W

    # En-route: clean config, OEI
    Vs_cl  = V_stall(W, S, rho, ac.CL_max_clean)
    V_enr  = 1.3 * Vs_cl
    q_enr  = 0.5 * rho * V_enr**2
    CL_enr = W / (q_enr * S)
    CD_enr = CD0 + K * CL_enr**2
    D_enr  = q_enr * S * CD_enr
    grad_enr = (T_oei - D_enr) / W

    # Approach: landing config, OEI, 1.3*Vs_L at landing weight
    W_land   = ac.W_landing
    CD0_land = CD0 + 0.020 + 0.015   # landing flaps + gear
    Vs_land  = V_stall(W_land, S, rho, ac.CL_max_L)
    Va       = 1.3 * Vs_land
    qa       = 0.5 * rho * Va**2
    CLa      = W_land / (qa * S)
    CDa      = CD0_land + K * CLa**2
    Da       = qa * S * CDa
    grad_app = (T_oei - Da) / W_land

    return {
        "2nd_segment":   {"gradient_pct": grad_2nd * 100, "min_pct": 2.4},
        "en_route":      {"gradient_pct": grad_enr * 100, "min_pct": 1.2},
        "approach_climb": {"gradient_pct": grad_app * 100, "min_pct": 2.1},
    }


def _time_to_climb_estimate(ac):
    """Estimate time to climb from SL to cruise altitude [minutes]."""
    W, S, CD0, K = ac.W_TO, ac.S, ac.CD0, ac.K
    h_cr = ac.h_cruise_ft
    n_steps = 8
    alts = np.linspace(0, h_cr, n_steps + 1)
    dt_total = 0.0

    for i in range(n_steps):
        h_lo, h_hi = alts[i], alts[i + 1]
        h_mid = 0.5 * (h_lo + h_hi)
        rho_h = isa_density(h_mid)
        T_h   = thrust_at_altitude(ac.T_max_SL, h_mid, ac.BPR)
        TW_h  = T_h / W
        V_h   = V_best_ROC_jet(W, S, rho_h, CD0, K, TW_h)
        roc_h = ROC_jet(V_h, T_h, CD0, K, W, S, rho_h)
        roc_h = max(roc_h, 1.0)
        dt_total += (h_hi - h_lo) / roc_h

    return dt_total / 60.0   # minutes
