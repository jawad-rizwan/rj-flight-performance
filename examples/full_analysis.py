#!/usr/bin/env python3
"""
Full performance analysis for regional jet variants.

Uses methods from Raymer Chapter 17: Performance and Flight Mechanics.
Runs CRJ and ZRJ families separately, saving output to examples/output/.
"""

import sys, os, io, contextlib
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

import numpy as np
from data import crj700, crj1000, zrj50, zrj70, zrj100
from perf import (
    # atmosphere
    isa_density, speed_of_sound, TAS_from_mach, fps_to_kts, kts_to_fps,
    thrust_at_altitude, sigma, dynamic_pressure, RHO_SL, G,
    # level flight
    thrust_required, power_required, TW_level,
    LD_max, CL_max_LD, V_min_thrust, V_min_power,
    V_stall, V_level, V_max, LD_from_CL,
    # range & endurance
    breguet_range_jet_nmi, V_best_range_jet, CL_best_range_jet,
    endurance_jet_hr, loiter_from_cruise,
    # climb
    climb_angle, rate_of_climb, V_best_ROC_jet, ROC_jet,
    time_to_climb, fuel_to_climb, time_to_climb_profile,
    far25_climb_gradients,
    # turning
    turn_rate_deg, turn_radius, corner_speed,
    n_sustained, bank_angle,
    sustained_turn_envelope_rho,
    # glide
    V_best_glide, LD_max_glide, min_sink_rate, glide_range,
    # energy
    specific_energy, Ps_expanded,
    # takeoff
    total_takeoff_distance, balanced_field_length, find_V1,
    # landing
    total_landing_distance,
    # wave drag
    mdd_wing, cd0_at_mach,
)

# =====================================================================
# Aircraft families
# =====================================================================
FAMILIES = {
    "CRJ": [crj700, crj1000],
    "ZRJ": [zrj50, zrj70, zrj100],
}

# =====================================================================
# Helper
# =====================================================================
def separator(title, width=70):
    print(f"\n{'=' * width}")
    print(f"  {title}")
    print(f"{'=' * width}")


def subsection(title):
    print(f"\n--- {title} ---")

# =====================================================================
# Main analysis
# =====================================================================
def analyse(ac):
    """Run complete performance analysis on one aircraft."""

    separator(f"{ac.name}  PERFORMANCE ANALYSIS  (Raymer Ch.17 Methods)")

    W   = ac.W_TO
    S   = ac.S
    CD0 = ac.CD0
    K   = ac.K

    # ---- Cruise conditions ----
    h_cr = ac.h_cruise_ft
    rho_cr = isa_density(h_cr)
    a_cr = speed_of_sound(h_cr)
    V_cr = TAS_from_mach(ac.M_cruise, h_cr)
    q_cr = dynamic_pressure(V_cr, h_cr)
    T_avail_cr = thrust_at_altitude(ac.T_max_SL, h_cr, ac.BPR)

    CL_cr = W / (q_cr * S)
    CD_cr = CD0 + K * CL_cr**2
    LD_cr = CL_cr / CD_cr

    print(f"\n  Aircraft       : {ac.name}")
    print(f"  MTOW           : {W:,.0f} lb")
    print(f"  Wing area      : {S:.1f} ft^2")
    print(f"  Wing loading   : {ac.wing_loading:.1f} lb/ft^2")
    print(f"  T/W (SL)       : {ac.thrust_to_weight:.3f}")
    print(f"  CD0            : {CD0:.4f}")
    print(f"  K              : {K:.4f}")
    print(f"  Oswald e       : {ac.e:.2f}")

    # =================================================================
    # 17.2  STEADY LEVEL FLIGHT
    # =================================================================
    subsection("17.2  STEADY LEVEL FLIGHT")

    # Cruise
    print(f"  Cruise Mach       : {ac.M_cruise}")
    print(f"  Cruise altitude   : {h_cr:,.0f} ft")
    print(f"  Cruise TAS        : {fps_to_kts(V_cr):.1f} kts  ({V_cr:.1f} ft/s)")
    print(f"  Cruise q          : {q_cr:.1f} lb/ft^2")
    print(f"  Cruise CL         : {CL_cr:.4f}")
    print(f"  Cruise CD         : {CD_cr:.5f}")
    print(f"  Cruise L/D        : {LD_cr:.2f}")
    print(f"  Thrust available  : {T_avail_cr:,.0f} lb")
    print(f"  Thrust required   : {thrust_required(W, CD0, K, S, rho_cr, V_cr):,.0f} lb")

    # Max L/D  (Eq. 17.67)
    ld_max = LD_max(CD0, K)
    cl_ldmax = CL_max_LD(CD0, K)
    v_ldmax = V_min_thrust(W, S, rho_cr, CD0, K)
    print(f"\n  (L/D)_max         : {ld_max:.2f}  at CL = {cl_ldmax:.4f}")
    print(f"  V for max L/D     : {fps_to_kts(v_ldmax):.1f} kts  (Eq. 17.13)")

    # Min power
    v_mp = V_min_power(W, S, rho_cr, CD0, K)
    print(f"  V for min power   : {fps_to_kts(v_mp):.1f} kts  (Eq. 17.19)")
    print(f"  Ratio V_mp/V_md   : {v_mp/v_ldmax:.3f}  (should be ~0.76)")

    # Stall speeds
    rho_sl = RHO_SL
    vs_clean = V_stall(W, S, rho_sl, ac.CL_max_clean)
    vs_to    = V_stall(W, S, rho_sl, ac.CL_max_TO)
    vs_land  = V_stall(ac.W_landing, S, rho_sl, ac.CL_max_L)
    print(f"\n  V_stall clean (SL): {fps_to_kts(vs_clean):.1f} kts")
    print(f"  V_stall TO   (SL): {fps_to_kts(vs_to):.1f} kts")
    print(f"  V_stall land (SL): {fps_to_kts(vs_land):.1f} kts  (at 85% MTOW)")

    # Maximum speed (Sec. 17.2.3)
    v_max_sl = V_max(W, S, rho_sl, CD0, K, ac.T_max_SL)
    v_max_cr = V_max(W, S, rho_cr, CD0, K, T_avail_cr)
    print(f"\n  V_max (SL)        : {fps_to_kts(v_max_sl):.1f} kts  (Sec. 17.2.3)")
    print(f"  V_max (cruise alt): {fps_to_kts(v_max_cr):.1f} kts")

    # =================================================================
    # 17.2.4  RANGE
    # =================================================================
    subsection("17.2.5  RANGE (Jet)")

    C_hr = ac.TSFC  # [1/hr]
    W_fuel_cruise = 0.90 * ac.W_fuel_max  # ~90% of fuel for cruise (rest for reserves)
    Wi = W
    Wf = W - W_fuel_cruise

    # Best-range velocity (Eq. 17.25)
    v_br = V_best_range_jet(W, S, rho_cr, CD0, K)
    cl_br = CL_best_range_jet(CD0, K)
    cd_br = CD0 + K * cl_br**2
    ld_br = cl_br / cd_br

    print(f"  V best range      : {fps_to_kts(v_br):.1f} kts  (Eq. 17.25)")
    print(f"  CL best range     : {cl_br:.4f}  (Eq. 17.26)")
    print(f"  L/D best range    : {ld_br:.2f}  (~86.6% of L/D_max = {0.866*ld_max:.2f})")

    # Breguet range at cruise
    R_cruise = breguet_range_jet_nmi(fps_to_kts(V_cr), C_hr, LD_cr, Wi, Wf)
    R_best   = breguet_range_jet_nmi(fps_to_kts(v_br), C_hr, ld_br, Wi, Wf)

    print(f"\n  Breguet range (cruise M={ac.M_cruise}): {R_cruise:,.0f} nmi  (Eq. 17.23)")
    print(f"  Breguet range (best-range V)    : {R_best:,.0f} nmi")
    print(f"  Fuel used         : {W_fuel_cruise:,.0f} lb  ({W_fuel_cruise/W*100:.1f}% MTOW)")

    # =================================================================
    # 17.2.7  ENDURANCE / LOITER
    # =================================================================
    subsection("17.2.8  ENDURANCE / LOITER (Jet)")

    W_fuel_loiter = 0.5 * ac.W_fuel_max  # use half fuel for loiter example
    E_max_ld = endurance_jet_hr(C_hr, ld_max, W, W - W_fuel_loiter)
    E_cruise = endurance_jet_hr(C_hr, LD_cr, W, W - W_fuel_loiter)

    print(f"  Max endurance (at L/D_max): {E_max_ld:.2f} hr  (Eq. 17.30)")
    print(f"  Endurance at cruise L/D  : {E_cruise:.2f} hr")
    print(f"  Loiter from cruise (est) : {loiter_from_cruise(R_cruise, fps_to_kts(V_cr)):.2f} hr  (Eq. 17.34)")

    # =================================================================
    # 17.3  CLIMB
    # =================================================================
    subsection("17.3  STEADY CLIMB")

    # Sea-level, max thrust
    T_sl = ac.T_max_SL
    D_sl = thrust_required(W, CD0, K, S, rho_sl, vs_clean * 1.3)  # at ~1.3 Vs
    gamma_sl = climb_angle(T_sl, D_sl, W)
    V_climb_sl = vs_clean * 1.3
    roc_sl = rate_of_climb(V_climb_sl, T_sl, D_sl, W)

    print(f"  Sea-level max thrust    : {T_sl:,.0f} lb")
    print(f"  Climb speed (1.3 Vs)    : {fps_to_kts(V_climb_sl):.1f} kts")
    print(f"  Climb angle (SL)        : {np.degrees(gamma_sl):.1f} deg  (Eq. 17.38)")
    print(f"  Rate of climb (SL)      : {roc_sl:.0f} ft/s  ({roc_sl*60:.0f} fpm)  (Eq. 17.39)")

    # Best ROC velocity (Eq. 17.43)
    TW_sl = T_sl / W
    V_broc = V_best_ROC_jet(W, S, rho_sl, CD0, K, TW_sl)
    roc_best = ROC_jet(V_broc, T_sl, CD0, K, W, S, rho_sl)

    print(f"\n  V best ROC (SL)         : {fps_to_kts(V_broc):.1f} kts  (Eq. 17.43)")
    print(f"  Best ROC (SL)           : {roc_best:.0f} ft/s  ({roc_best*60:.0f} fpm)")

    # Climb profile (SL to cruise alt)
    alt_steps = np.linspace(0, h_cr, 8)
    rocs = []
    for h in alt_steps:
        rho_h = isa_density(h)
        T_h = thrust_at_altitude(T_sl, h, ac.BPR)
        TW_h = T_h / W
        V_h = V_best_ROC_jet(W, S, rho_h, CD0, K, TW_h)
        roc_h = ROC_jet(V_h, T_h, CD0, K, W, S, rho_h)
        rocs.append(max(roc_h, 0.1))  # keep positive
    rocs = np.array(rocs)

    times = time_to_climb_profile(alt_steps, rocs)
    fuel_climb = fuel_to_climb(ac.TSFC / 3600, T_sl * 0.9, times[-1])

    print(f"\n  Climb profile (SL -> {h_cr/1000:.0f}k ft):")
    print(f"  {'Alt (ft)':>10}  {'ROC (fpm)':>10}  {'Time (min)':>10}")
    for h, r, t in zip(alt_steps, rocs, times):
        print(f"  {h:10,.0f}  {r*60:10,.0f}  {t/60:10.1f}")
    print(f"  Total time to climb   : {times[-1]/60:.1f} min  (Eqs. 17.48-17.50)")
    print(f"  Est. fuel to climb    : {fuel_climb:,.0f} lb  (Eq. 17.51)")

    # =================================================================
    # 17.4  TURNING FLIGHT
    # =================================================================
    subsection("17.4  LEVEL TURNING FLIGHT")

    # Corner speed at SL
    V_corner = corner_speed(ac.wing_loading, rho_sl, ac.CL_max_clean, ac.n_max)
    psi_corner = turn_rate_deg(V_corner, ac.n_max)
    R_corner = turn_radius(V_corner, ac.n_max)

    print(f"  n_max (structural)      : {ac.n_max}")
    print(f"  Corner speed (SL)       : {fps_to_kts(V_corner):.1f} kts")
    print(f"  Turn rate at corner     : {psi_corner:.1f} deg/s  (Eq. 17.52)")
    print(f"  Turn radius at corner   : {R_corner:,.0f} ft  (Eq. 17.79)")
    print(f"  Bank angle at n_max     : {np.degrees(bank_angle(ac.n_max)):.1f} deg")

    # Sustained turn at cruise
    n_sust = n_sustained(T_avail_cr / W, q_cr, ac.wing_loading, CD0, K)
    if n_sust > 1.0:
        psi_sust = turn_rate_deg(V_cr, n_sust)
        R_sust = turn_radius(V_cr, n_sust)
        print(f"\n  Sustained turn (cruise):")
        print(f"    n_sustained           : {n_sust:.2f}  (Eq. 17.54)")
        print(f"    Turn rate             : {psi_sust:.2f} deg/s")
        print(f"    Turn radius           : {R_sust:,.0f} ft")
        print(f"    Bank angle            : {np.degrees(bank_angle(n_sust)):.1f} deg")
    else:
        print(f"    n_sustained < 1 at cruise -- cannot sustain turn")

    # =================================================================
    # 17.5  GLIDE
    # =================================================================
    subsection("17.5  GLIDING FLIGHT")

    V_bg = V_best_glide(W, S, rho_cr, CD0, K)
    ld_glide = LD_max_glide(CD0, K)
    sr = min_sink_rate(W, S, rho_cr, CD0, K)
    R_glide = glide_range(h_cr, ld_glide)

    print(f"  Best glide L/D          : {ld_glide:.2f}  (Eq. 17.67)")
    print(f"  V best glide (at {h_cr/1000:.0f}k) : {fps_to_kts(V_bg):.1f} kts  (Eq. 17.65)")
    print(f"  Min sink rate           : {sr:.1f} ft/s  ({sr*60:.0f} fpm)")
    print(f"  Glide range from {h_cr/1000:.0f}k ft: {R_glide/6076:.0f} nmi  ({R_glide/5280:.0f} SM)")

    # =================================================================
    # 17.6  ENERGY MANEUVERABILITY
    # =================================================================
    subsection("17.6  ENERGY MANEUVERABILITY")

    he_cr = specific_energy(h_cr, V_cr)
    Ps_cr = Ps_expanded(V_cr, T_avail_cr/W, ac.wing_loading, CD0, K, rho_cr, n=1.0)

    print(f"  Energy height (cruise)  : {he_cr:,.0f} ft  (Eq. 17.85)")
    print(f"  Ps at cruise (n=1)      : {Ps_cr:.1f} ft/s  ({Ps_cr*60:.0f} fpm)  (Eq. 17.89)")

    # Ps at SL, various n
    for n_val in [1, 2, 3]:
        ps_val = Ps_expanded(V_climb_sl, TW_sl, ac.wing_loading, CD0, K, rho_sl, n=n_val)
        print(f"  Ps at SL, n={n_val}           : {ps_val:.1f} ft/s  ({ps_val*60:.0f} fpm)")

    # =================================================================
    # 17.7  SERVICE CEILING
    # =================================================================
    subsection("17.7  SERVICE CEILING")

    # Find altitude where ROC = 500 fpm for jet service ceiling
    # Cap speed at Mmo and account for compressibility drag rise
    svc_ceil = None
    abs_ceil = None
    for h_test in np.arange(0, 60001, 500):
        rho_t = isa_density(h_test)
        a_t = speed_of_sound(h_test)
        T_t = thrust_at_altitude(T_sl, h_test, ac.BPR)
        V_t = V_best_ROC_jet(W, S, rho_t, CD0, K, T_t / W)

        # Cap at Mmo
        if ac.M_mo > 0:
            V_mmo = ac.M_mo * a_t
            V_t = min(V_t, V_mmo)

        # Wave drag: adjust CD0 for compressibility
        M_t = V_t / a_t
        q_t = 0.5 * rho_t * V_t**2
        CL_t = W / (q_t * S) if q_t > 0 else 0
        mdd_t = mdd_wing(ac.t_c, ac.sweep_qc_deg, CL_t)
        CD0_t = cd0_at_mach(CD0, M_t, mdd_t)

        roc_t = ROC_jet(V_t, T_t, CD0_t, K, W, S, rho_t)
        if svc_ceil is None and roc_t * 60 <= 500.0:
            svc_ceil = h_test
        if roc_t <= 0.5:
            abs_ceil = h_test
            break

    if svc_ceil:
        print(f"  Service ceiling (~500 fpm): {svc_ceil:,.0f} ft")
    else:
        print(f"  Service ceiling           : > 60,000 ft")
    if abs_ceil:
        print(f"  Absolute ceiling (~0 fpm) : {abs_ceil:,.0f} ft")
    else:
        print(f"  Absolute ceiling          : > 60,000 ft")

    # =================================================================
    # 17.8  TAKEOFF
    # =================================================================
    subsection("17.8  TAKEOFF ANALYSIS (Sea Level, ISA)")

    CD0_TO = CD0 + 0.015   # gear + flap drag increment
    # Raymer Eq. 17.114: average thrust = T_static * 0.75*(5+BPR)/(4+BPR)
    T_ratio = 0.75 * (5.0 + ac.BPR) / (4.0 + ac.BPR)
    T_TO = T_ratio * T_sl

    # All-engine TODR
    to = total_takeoff_distance(
        W=W, S=S, T=T_TO, CD0=CD0_TO,
        CL_ground=ac.CL_ground, K=K, mu=ac.mu_roll,
        rho=rho_sl, CL_max_TO=ac.CL_max_TO, TW=ac.thrust_to_weight,
        h_obstacle=ac.h_obstacle_TO, t_rotate=3.0)

    # V1 / BFL (iterative)
    T_idle_to = 0.05 * T_sl
    v1_result = find_V1(
        W=W, S=S, T=T_TO, CD0=CD0_TO,
        CL_ground=ac.CL_ground, K=K,
        mu_roll=ac.mu_roll, mu_brake=ac.mu_brake,
        rho=rho_sl, CL_max_TO=ac.CL_max_TO, TW=ac.thrust_to_weight,
        h_obstacle=ac.h_obstacle_TO, n_engines=ac.n_engines,
        T_idle=T_idle_to, T_reverse=0.0, t_react=2.0, t_rotate=3.0)

    # Empirical BFL (Eq. 17.113)
    gamma_min = 0.024  # 2-engine
    U = 0.01 * ac.CL_max_TO + 0.02
    CL_climb = ac.CL_max_TO / (1.2**2)  # CL at 1.2*Vs
    G_climb = np.arcsin(np.clip(ac.thrust_to_weight * 0.5 - 1.0 /
              (CL_climb / (CD0 + K * CL_climb**2)), -1, 1)) - gamma_min

    bfl_empirical = balanced_field_length(
        W_over_S=ac.wing_loading, CL_climb=CL_climb,
        h_obstacle=ac.h_obstacle_TO, TW=ac.thrust_to_weight,
        U=U, G_climb=max(G_climb, 0.01), BPR=ac.BPR,
        N_e=ac.n_engines, rho=rho_sl, is_prop=False)

    print(f"\n  V-speeds:")
    print(f"    V_stall (TO config)   : {fps_to_kts(to['V_stall_TO']):.1f} kts")
    print(f"    V1  (decision)        : {fps_to_kts(v1_result['V1']):.1f} kts")
    print(f"    VR  (rotation, 1.1 Vs): {fps_to_kts(to['V_R']):.1f} kts")
    print(f"    V2  (safety, 1.2 Vs)  : {fps_to_kts(to['V_2']):.1f} kts")

    print(f"\n  TODR (all-engine):")
    print(f"    Ground roll           : {to['S_ground_roll']:,.0f} ft  (Eq. 17.102)")
    print(f"    Rotation              : {to['S_rotation']:,.0f} ft")
    print(f"    Transition            : {to['S_transition']:,.0f} ft")
    print(f"    Climb to {ac.h_obstacle_TO:.0f} ft        : {to['S_climb']:,.0f} ft  (Eq. 17.112)")
    print(f"    Climb angle           : {to['gamma_deg']:.1f} deg")
    print(f"    TODR (unfactored)     : {to['TODR']:,.0f} ft")
    print(f"    TODR (FAR x1.15)      : {to['TODR_factored']:,.0f} ft")

    print(f"\n  Balanced field / V1:")
    print(f"    ASDR at V1            : {v1_result['ASDR_at_V1']:,.0f} ft")
    print(f"    AGDR at V1            : {v1_result['AGDR_at_V1']:,.0f} ft")
    print(f"    BFL (iterative)       : {v1_result['BFL']:,.0f} ft")
    print(f"    BFL (Raymer Eq.17.113): {bfl_empirical:,.0f} ft")

    # TOFL = max(FAR-factored all-engine TODR, BFL)
    TOFL = max(to['TODR_factored'], v1_result['BFL'])
    print(f"\n  >> TOFL (FAR 25)        : {TOFL:,.0f} ft  = max(TODR×1.15, BFL)")

    # =================================================================
    # FAR 25 CLIMB GRADIENTS (Table F.4)
    # =================================================================
    subsection("FAR 25 CLIMB GRADIENTS (Table F.4, 2-engine)")

    far_grads = far25_climb_gradients(
        W=W, S=S, CD0=CD0, K=K,
        CL_max_TO=ac.CL_max_TO, CL_max_L=ac.CL_max_L,
        CL_max_clean=ac.CL_max_clean,
        T_max_SL=ac.T_max_SL, n_engines=ac.n_engines, rho=rho_sl)

    all_pass = True
    for seg_key in ["1st_seg", "2nd_seg", "4th_seg", "GA_approach", "GA_landing"]:
        seg = far_grads[seg_key]
        status = "PASS" if seg["pass"] else "** FAIL **"
        if not seg["pass"]:
            all_pass = False
        print(f"    {seg['label']}")
        print(f"      Gradient: {seg['gradient']*100:.1f}%  "
              f"(required >= {seg['required']*100:.1f}%)  "
              f"@ {seg['speed_kts']:.0f} kts  [{status}]")

    print(f"\n  >> FAR 25 climb compliance: {'ALL PASS' if all_pass else 'FAIL — see above'}")

    # =================================================================
    # 17.9  LANDING
    # =================================================================
    subsection("17.9  LANDING ANALYSIS (Sea Level, ISA)")

    W_land = ac.W_landing
    T_idle_land = 0.05 * T_sl   # idle thrust ~ 5%

    la = total_landing_distance(
        W=W_land, S=S, rho=rho_sl, CL_max_L=ac.CL_max_L,
        CD0=CD0 + 0.02, CL_ground=ac.CL_ground, K=K,
        mu_brake=ac.mu_brake, h_obstacle=ac.h_obstacle_L,
        T_idle=T_idle_land, T_reverse=0.0,
        approach_factor=1.3, t_free=3.0, FAR_factor=True)

    print(f"\n  V-speeds:")
    print(f"    V_stall (landing)     : {fps_to_kts(la['V_stall']):.1f} kts")
    print(f"    V_approach (1.3 Vs)   : {fps_to_kts(la['V_approach']):.1f} kts")
    print(f"    V_touchdown (1.15 Vs) : {fps_to_kts(la['V_TD']):.1f} kts")

    print(f"\n  Landing weight          : {W_land:,.0f} lb  (85% MTOW)")
    print(f"    Approach distance     : {la['S_approach']:,.0f} ft")
    print(f"    Flare distance        : {la['S_flare']:,.0f} ft")
    print(f"    Free roll             : {la['S_free_roll']:,.0f} ft")
    print(f"    Braking distance      : {la['S_braking']:,.0f} ft")
    ldr_wet = la['S_FAR_field'] * 1.15

    print(f"    LDR (unfactored)      : {la['S_total_actual']:,.0f} ft")
    print(f"    LDR (FAR /0.6)        : {la['S_FAR_field']:,.0f} ft")
    print(f"    LDR wet (FAR x1.15)   : {ldr_wet:,.0f} ft  (FAR 121.195d)")

    # =================================================================
    # Summary table
    # =================================================================
    separator(f"{ac.name}  PERFORMANCE SUMMARY")
    summary = [
        ("MTOW",                    f"{W:,.0f} lb"),
        ("Wing loading (W/S)",      f"{ac.wing_loading:.1f} lb/ft^2"),
        ("T/W (SL static)",         f"{ac.thrust_to_weight:.3f}"),
        ("(L/D)_max",               f"{ld_max:.2f}"),
        ("Cruise L/D",              f"{LD_cr:.2f}"),
        ("Cruise Mach",             f"{ac.M_cruise}"),
        ("V_stall clean (SL)",      f"{fps_to_kts(vs_clean):.0f} kts"),
        ("V1 / VR / V2",           f"{fps_to_kts(v1_result['V1']):.0f} / {fps_to_kts(to['V_R']):.0f} / {fps_to_kts(to['V_2']):.0f} kts"),
        ("Best range (cruise)",     f"{R_cruise:,.0f} nmi"),
        ("Best ROC (SL)",           f"{roc_best*60:,.0f} fpm"),
        ("Service ceiling",         "see above"),
        ("Corner speed (SL)",       f"{fps_to_kts(V_corner):.0f} kts"),
        ("TODR (unfactored)",       f"{to['TODR']:,.0f} ft"),
        ("TODR (FAR x1.15)",        f"{to['TODR_factored']:,.0f} ft"),
        ("ASDR at V1",              f"{v1_result['ASDR_at_V1']:,.0f} ft"),
        ("BFL (iterative)",         f"{v1_result['BFL']:,.0f} ft"),
        ("BFL (Raymer Eq.17.113)",  f"{bfl_empirical:,.0f} ft"),
        ("TOFL (FAR 25)",           f"{TOFL:,.0f} ft"),
        ("LDR (unfactored)",        f"{la['S_total_actual']:,.0f} ft"),
        ("LDR (FAR /0.6)",          f"{la['S_FAR_field']:,.0f} ft"),
        ("LDR wet (FAR x1.15)",     f"{ldr_wet:,.0f} ft"),
        ("Glide range from cruise", f"{R_glide/6076:,.0f} nmi"),
    ]
    for label, value in summary:
        print(f"  {label:<28s}: {value}")

    print()
    return {
        "ac": ac, "LD_max": ld_max, "LD_cr": LD_cr,
        "R_cruise": R_cruise, "ROC_best_sl": roc_best,
        "V_corner": V_corner, "TODR": to["TODR"],
        "TODR_factored": to["TODR_factored"],
        "ASDR": v1_result["ASDR_at_V1"],
        "BFL": v1_result["BFL"],
        "BFL_empirical": bfl_empirical,
        "TOFL": TOFL,
        "LDR": la["S_total_actual"],
        "LDR_FAR": la["S_FAR_field"],
        "LDR_wet": ldr_wet,
    }


# =====================================================================
def run_comparison(results):
    """Print comparison table for a set of analysis results."""
    separator("VARIANT COMPARISON")
    header = f"  {'Metric':<28s}"
    for name in results:
        header += f"  {name:>12s}"
    print(header)
    print("  " + "-" * (28 + 14 * len(results)))

    metrics = [
        ("(L/D)_max",          "LD_max",         ".2f",  None),
        ("Cruise L/D",         "LD_cr",          ".2f",  None),
        ("Range (nmi)",        "R_cruise",       ",.0f", None),
        ("Best ROC (fpm)",     "ROC_best_sl",    ",.0f", lambda v: v * 60),
        ("Corner spd (kts)",   "V_corner",       ",.0f", fps_to_kts),
        ("TODR (ft)",          "TODR",           ",.0f", None),
        ("TODR FAR (ft)",      "TODR_factored",  ",.0f", None),
        ("ASDR at V1 (ft)",    "ASDR",           ",.0f", None),
        ("BFL iterative (ft)", "BFL",            ",.0f", None),
        ("BFL Eq.17.113 (ft)", "BFL_empirical",  ",.0f", None),
        ("TOFL FAR 25 (ft)",   "TOFL",           ",.0f", None),
        ("LDR (ft)",           "LDR",            ",.0f", None),
        ("LDR FAR (ft)",       "LDR_FAR",        ",.0f", None),
        ("LDR wet (ft)",       "LDR_wet",        ",.0f", None),
    ]
    for label, key, fmt, transform in metrics:
        row = f"  {label:<28s}"
        for name in results:
            val = results[name][key]
            if transform:
                val = transform(val)
            row += f"  {val:>12{fmt}}"
        print(row)


# =====================================================================
# Raymer Table 17.1 — surface coefficients (midpoint values)
# =====================================================================
SURFACES = [
    # (name,              mu_roll, mu_brake)
    ("Dry concrete",       0.03,    0.40),
    ("Wet concrete",       0.05,    0.225),
    ("Icy concrete",       0.02,    0.08),
    ("Hard turf",          0.05,    0.40),
    ("Firm dirt",          0.04,    0.30),
    ("Soft turf",          0.07,    0.20),
    ("Wet grass",          0.08,    0.20),
]


def run_surface_performance(variants):
    """Airfield performance on different runway surfaces (Raymer Table 17.1)."""
    separator("AIRFIELD PERFORMANCE BY RUNWAY SURFACE  (Table 17.1)")

    for ac in variants:
        subsection(f"{ac.name}")
        W, S, CD0, K = ac.W_TO, ac.S, ac.CD0, ac.K
        rho = RHO_SL
        CD0_TO = CD0 + 0.015
        T_ratio = 0.75 * (5.0 + ac.BPR) / (4.0 + ac.BPR)
        T_TO = T_ratio * ac.T_max_SL
        T_idle = 0.05 * ac.T_max_SL
        W_land = ac.W_landing

        header = (f"  {'Surface':<18s} {'mu_r':>5s} {'mu_b':>5s}"
                  f" {'TOFL':>8s} {'BFL':>8s}"
                  f" {'LDR':>8s} {'LDR FAR':>8s} {'Wet FAR':>8s}")
        print(header)
        print("  " + "-" * len(header.strip()))

        for name, mu_r, mu_b in SURFACES:
            # Takeoff: TOFL = max(TODR × 1.15, BFL)
            to = total_takeoff_distance(
                W=W, S=S, T=T_TO, CD0=CD0_TO,
                CL_ground=ac.CL_ground, K=K, mu=mu_r,
                rho=rho, CL_max_TO=ac.CL_max_TO, TW=ac.thrust_to_weight,
                h_obstacle=ac.h_obstacle_TO, t_rotate=3.0)

            v1 = find_V1(
                W=W, S=S, T=T_TO, CD0=CD0_TO,
                CL_ground=ac.CL_ground, K=K,
                mu_roll=mu_r, mu_brake=mu_b,
                rho=rho, CL_max_TO=ac.CL_max_TO, TW=ac.thrust_to_weight,
                h_obstacle=ac.h_obstacle_TO, n_engines=ac.n_engines,
                T_idle=T_idle, T_reverse=0.0, t_react=2.0, t_rotate=3.0)

            tofl = max(to['TODR_factored'], v1['BFL'])

            # Landing
            la = total_landing_distance(
                W=W_land, S=S, rho=rho, CL_max_L=ac.CL_max_L,
                CD0=CD0 + 0.02, CL_ground=ac.CL_ground, K=K,
                mu_brake=mu_b, h_obstacle=ac.h_obstacle_L,
                T_idle=T_idle, T_reverse=0.0,
                approach_factor=1.3, t_free=3.0, FAR_factor=True)

            ldr = la['S_total_actual']
            ldr_far = la['S_FAR_field']
            # FAR 121.195(d): wet = dry FAR × 1.15
            wet_far = ldr_far * 1.15 if "Dry" in name else ldr_far

            wet_str = f"{wet_far:8,.0f}" if "Dry" in name else "       -"

            print(f"  {name:<18s} {mu_r:5.2f} {mu_b:5.2f}"
                  f" {tofl:8,.0f} {v1['BFL']:8,.0f}"
                  f" {ldr:8,.0f} {ldr_far:8,.0f} {wet_str}")

        print()
        print("  TOFL = max(TODR × 1.15, BFL)  |  "
              "LDR FAR = actual / 0.6  |  "
              "Wet FAR = dry FAR × 1.15 (FAR 121.195d)")
        print()


# =====================================================================
if __name__ == "__main__":
    OUTPUT_DIR = os.path.join(os.path.dirname(__file__), "output")

    for family_name, variants in FAMILIES.items():
        family_dir = os.path.join(OUTPUT_DIR, family_name)
        os.makedirs(family_dir, exist_ok=True)
        output_path = os.path.join(family_dir, "analysis.txt")

        # Capture output to both stdout and file
        buf = io.StringIO()

        class Tee:
            def __init__(self, *streams):
                self.streams = streams
            def write(self, data):
                for s in self.streams:
                    s.write(data)
            def flush(self):
                for s in self.streams:
                    s.flush()

        tee = Tee(sys.stdout, buf)
        old_stdout = sys.stdout
        sys.stdout = tee

        try:
            results = {}
            for ac in variants:
                results[ac.name] = analyse(ac)
            run_comparison(results)
            run_surface_performance(variants)
        finally:
            sys.stdout = old_stdout

        with open(output_path, "w") as f:
            f.write(buf.getvalue())
        print(f"\n  Output saved to {output_path}")
