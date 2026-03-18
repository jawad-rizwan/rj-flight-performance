#!/usr/bin/env python3
"""
Generate performance charts for regional jet variants.
Uses methods from Raymer Chapter 17.

Charts produced:
  1. Thrust required & available vs velocity
  2. Power required & available vs velocity
  3. L/D vs CL
  4. Rate of climb vs altitude
  5. Specific excess power (Ps) vs Mach
  6. Turn rate & load factor vs velocity
  7. Payload-range diagram
  8. Takeoff & landing distance breakdown
  9. ASDR & AGDR vs engine failure speed (BFL)
 10. Glide polar (sink rate vs velocity)
"""

import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from data import crj700, crj1000, zrj50, zrj70, zrj100
from perf import (
    isa_density, speed_of_sound, TAS_from_mach, fps_to_kts, kts_to_fps,
    thrust_at_altitude, dynamic_pressure, RHO_SL, G,
    thrust_required, power_required,
    LD_max, CL_max_LD, V_min_thrust, V_min_power, V_stall, LD_from_CL,
    breguet_range_jet_nmi, V_best_range_jet, CL_best_range_jet,
    ROC_jet, V_best_ROC_jet,
    turn_rate_deg, turn_radius, corner_speed,
    n_sustained, sustained_turn_envelope_rho,
    Ps_expanded, specific_energy,
    V_best_glide, min_sink_rate, sink_rate,
    asdr_todr_curves, find_V1,
    mission_profile,
)

FAMILIES = {
    "CRJ": {
        "variants": [crj700, crj1000],
        "colors": {"CRJ-700": "tab:blue", "CRJ-1000": "tab:red"},
    },
    "ZRJ": {
        "variants": [zrj50, zrj70, zrj100],
        "colors": {"ZRJ50": "tab:green", "ZRJ70": "tab:orange", "ZRJ100": "tab:purple"},
    },
}

# Module-level globals — set per family run
VARIANTS = []
COLORS = {}
CHART_DIR = ""

def save(fig, name):
    path = os.path.join(CHART_DIR, name)
    fig.savefig(path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved: {path}")


# =====================================================================
# 1. Thrust Required & Available vs Velocity
# =====================================================================
def plot_thrust_vs_velocity():
    fig, axes = plt.subplots(1, len(VARIANTS), figsize=(14, 5), sharey=True)
    if len(VARIANTS) == 1:
        axes = [axes]

    for ax, ac in zip(axes, VARIANTS):
        W   = ac.W_TO
        S   = ac.S
        CD0 = ac.CD0
        K   = ac.K

        for h_ft, ls in [(0, "-"), (ac.h_cruise_ft, "--")]:
            rho = isa_density(h_ft)
            Vs = V_stall(W, S, rho, ac.CL_max_clean)
            V_arr = np.linspace(Vs * 0.9, kts_to_fps(450), 300)
            T_req = thrust_required(W, CD0, K, S, rho, V_arr)
            T_av = thrust_at_altitude(ac.T_max_SL, h_ft, ac.BPR)

            label_alt = f"{h_ft/1000:.0f}k ft"
            ax.plot(fps_to_kts(V_arr), T_req, ls, color=COLORS[ac.name],
                    label=f"T_req @ {label_alt}")
            ax.axhline(T_av, color=COLORS[ac.name], ls=ls, alpha=0.5,
                       label=f"T_avail @ {label_alt}")

        # Mark min-drag speed at SL
        rho_sl = RHO_SL
        v_md = V_min_thrust(W, S, rho_sl, CD0, K)
        T_md = thrust_required(W, CD0, K, S, rho_sl, v_md)
        ax.plot(fps_to_kts(v_md), T_md, 'o', color=COLORS[ac.name], ms=8)
        ax.annotate(f"  min drag\n  {fps_to_kts(v_md):.0f} kts",
                    (fps_to_kts(v_md), T_md), fontsize=8)

        ax.set_xlabel("True Airspeed [kts]")
        ax.set_title(f"{ac.name}")
        ax.legend(fontsize=7, loc="upper right")
        ax.set_ylim(0, ac.T_max_SL * 1.3)
        ax.grid(True, alpha=0.3)

    axes[0].set_ylabel("Thrust [lb]")
    fig.suptitle("Thrust Required & Available vs Velocity  (Raymer Eq. 17.8)", fontsize=12)
    save(fig, "01_thrust_vs_velocity.png")


# =====================================================================
# 2. Power Required & Available vs Velocity
# =====================================================================
def plot_power_vs_velocity():
    fig, axes = plt.subplots(1, len(VARIANTS), figsize=(14, 5), sharey=True)
    if len(VARIANTS) == 1:
        axes = [axes]

    for ax, ac in zip(axes, VARIANTS):
        W, S, CD0, K = ac.W_TO, ac.S, ac.CD0, ac.K
        rho = RHO_SL
        Vs = V_stall(W, S, rho, ac.CL_max_clean)
        V_arr = np.linspace(Vs * 0.9, kts_to_fps(400), 300)

        P_req = power_required(W, CD0, K, S, rho, V_arr) / 550.0  # to hp
        P_av = ac.T_max_SL * V_arr / 550.0  # max thrust * V

        ax.plot(fps_to_kts(V_arr), P_req, color=COLORS[ac.name], label="P required")

        v_mp = V_min_power(W, S, rho, CD0, K)
        P_mp = power_required(W, CD0, K, S, rho, v_mp) / 550.0
        ax.plot(fps_to_kts(v_mp), P_mp, 'o', color=COLORS[ac.name], ms=8)
        ax.annotate(f"  min power\n  {fps_to_kts(v_mp):.0f} kts",
                    (fps_to_kts(v_mp), P_mp), fontsize=8)

        ax.set_xlabel("True Airspeed [kts]")
        ax.set_title(f"{ac.name} (Sea Level)")
        ax.legend(fontsize=8)
        ax.grid(True, alpha=0.3)

    axes[0].set_ylabel("Power [hp]")
    fig.suptitle("Power Required vs Velocity  (Raymer Eq. 17.16-17.17)", fontsize=12)
    save(fig, "02_power_vs_velocity.png")


# =====================================================================
# 3. L/D vs CL
# =====================================================================
def plot_LD_vs_CL():
    fig, ax = plt.subplots(figsize=(8, 5))

    for ac in VARIANTS:
        CL_arr = np.linspace(0.05, 1.5, 200)
        LD_arr = LD_from_CL(CL_arr, ac.CD0, ac.K)
        ax.plot(CL_arr, LD_arr, color=COLORS[ac.name], label=ac.name, lw=2)

        cl_opt = CL_max_LD(ac.CD0, ac.K)
        ld_opt = LD_max(ac.CD0, ac.K)
        ax.plot(cl_opt, ld_opt, 'o', color=COLORS[ac.name], ms=8)
        ax.annotate(f"  L/D_max={ld_opt:.1f}\n  CL={cl_opt:.3f}",
                    (cl_opt, ld_opt), fontsize=8)

    ax.set_xlabel("$C_L$")
    ax.set_ylabel("L/D")
    ax.set_title("Lift-to-Drag Ratio vs $C_L$  (Raymer Eq. 17.14, 17.67)")
    ax.legend()
    ax.grid(True, alpha=0.3)
    save(fig, "03_LD_vs_CL.png")


# =====================================================================
# 4. Rate of Climb vs Altitude
# =====================================================================
def plot_ROC_vs_altitude():
    fig, ax = plt.subplots(figsize=(8, 6))

    for ac in VARIANTS:
        W, S, CD0, K = ac.W_TO, ac.S, ac.CD0, ac.K
        alts = np.linspace(0, 45000, 100)
        rocs = []
        for h in alts:
            rho_h = isa_density(h)
            T_h = thrust_at_altitude(ac.T_max_SL, h, ac.BPR)
            TW_h = T_h / W
            V_h = V_best_ROC_jet(W, S, rho_h, CD0, K, TW_h)
            roc = ROC_jet(V_h, T_h, CD0, K, W, S, rho_h)
            rocs.append(max(roc * 60, 0))  # fpm
        ax.plot(rocs, alts / 1000, color=COLORS[ac.name], label=ac.name, lw=2)

    ax.axvline(500, color='gray', ls='--', alpha=0.5, label="Service ceiling (500 fpm)")
    ax.set_xlabel("Best Rate of Climb [fpm]")
    ax.set_ylabel("Altitude [1000 ft]")
    ax.set_title("Rate of Climb vs Altitude  (Raymer Eq. 17.43, 17.39)")
    ax.legend()
    ax.grid(True, alpha=0.3)
    save(fig, "04_ROC_vs_altitude.png")


# =====================================================================
# 5. Specific Excess Power vs Mach
# =====================================================================
def plot_Ps_vs_Mach():
    fig, axes = plt.subplots(1, len(VARIANTS), figsize=(14, 5), sharey=True)
    if len(VARIANTS) == 1:
        axes = [axes]

    for ax, ac in zip(axes, VARIANTS):
        W, S, CD0, K = ac.W_TO, ac.S, ac.CD0, ac.K

        for h_ft in [0, 10000, 20000, 30000, ac.h_cruise_ft]:
            rho = isa_density(h_ft)
            a = speed_of_sound(h_ft)
            T_h = thrust_at_altitude(ac.T_max_SL, h_ft, ac.BPR)
            TW = T_h / W

            M_arr = np.linspace(0.2, 0.85, 200)
            V_arr = M_arr * a
            ps_arr = Ps_expanded(V_arr, TW, ac.wing_loading, CD0, K, rho, n=1.0)

            ax.plot(M_arr, ps_arr * 60, label=f"{h_ft/1000:.0f}k ft")

        ax.axhline(0, color='k', lw=0.5)
        ax.set_xlabel("Mach Number")
        ax.set_title(f"{ac.name}")
        ax.legend(fontsize=7, title="Altitude")
        ax.grid(True, alpha=0.3)
        ax.set_ylim(-2000, None)

    axes[0].set_ylabel("$P_s$  [ft/min]")
    fig.suptitle("Specific Excess Power vs Mach  (Raymer Eq. 17.89)", fontsize=12)
    save(fig, "05_Ps_vs_Mach.png")


# =====================================================================
# 6. Turn Rate & Load Factor vs Velocity
# =====================================================================
def plot_turn_performance():
    fig, axes = plt.subplots(2, len(VARIANTS), figsize=(14, 9), sharex="col")
    if len(VARIANTS) == 1:
        axes = axes.reshape(-1, 1)

    for col, ac in enumerate(VARIANTS):
        W, S, CD0, K = ac.W_TO, ac.S, ac.CD0, ac.K
        rho = RHO_SL
        T_h = ac.T_max_SL

        V_arr = np.linspace(kts_to_fps(100), kts_to_fps(500), 300)
        V_kts = fps_to_kts(V_arr)

        # Sustained turn
        psi_sust, n_sust = sustained_turn_envelope_rho(V_arr, T_h, W, S, CD0, K, rho)

        # Stall limit (instantaneous)
        q_arr = 0.5 * rho * V_arr**2
        n_stall = q_arr * ac.CL_max_clean / (W / S)
        n_stall = np.minimum(n_stall, ac.n_max)
        valid = n_stall > 1.0
        psi_stall = np.where(valid, turn_rate_deg(V_arr, np.where(valid, n_stall, 1.001)), 0)

        # Structural limit
        n_struct = ac.n_max * np.ones_like(V_arr)
        psi_struct = turn_rate_deg(V_arr, n_struct)

        # Top plot: turn rate
        ax_tr = axes[0, col] if axes.ndim == 2 else axes[0]
        ax_tr.plot(V_kts, psi_sust, color=COLORS[ac.name], lw=2, label="Sustained")
        ax_tr.plot(V_kts, psi_stall, '--', color=COLORS[ac.name], alpha=0.7, label="Stall limit")
        ax_tr.plot(V_kts, psi_struct, ':', color='gray', alpha=0.5, label=f"n={ac.n_max}")
        V_c = corner_speed(ac.wing_loading, rho, ac.CL_max_clean, ac.n_max)
        ax_tr.axvline(fps_to_kts(V_c), color='green', ls='--', alpha=0.5, label="Corner speed")
        ax_tr.set_ylabel("Turn Rate [deg/s]")
        ax_tr.set_title(f"{ac.name} (Sea Level)")
        ax_tr.legend(fontsize=7)
        ax_tr.set_ylim(0, 25)
        ax_tr.grid(True, alpha=0.3)

        # Bottom plot: load factor
        ax_n = axes[1, col] if axes.ndim == 2 else axes[1]
        ax_n.plot(V_kts, n_sust, color=COLORS[ac.name], lw=2, label="Sustained n")
        ax_n.plot(V_kts, n_stall, '--', color=COLORS[ac.name], alpha=0.7, label="Stall limit n")
        ax_n.axhline(ac.n_max, color='gray', ls=':', alpha=0.5, label=f"n_max={ac.n_max}")
        ax_n.set_xlabel("True Airspeed [kts]")
        ax_n.set_ylabel("Load Factor n")
        ax_n.legend(fontsize=7)
        ax_n.set_ylim(0, ac.n_max + 1)
        ax_n.grid(True, alpha=0.3)

    fig.suptitle("Turn Performance  (Raymer Eq. 17.52, 17.54)", fontsize=12)
    save(fig, "06_turn_performance.png")


# =====================================================================
# 7. Payload-Range Diagram
# =====================================================================
def plot_payload_range():
    fig, ax = plt.subplots(figsize=(9, 6))

    for ac in VARIANTS:
        W, S, CD0, K = ac.W_TO, ac.S, ac.CD0, ac.K
        h = ac.h_cruise_ft
        rho = isa_density(h)
        V = TAS_from_mach(ac.M_cruise, h)
        q = 0.5 * rho * V**2
        C_hr = ac.TSFC

        # Point A: max payload, max fuel (MTOW limit may reduce fuel)
        W_fuel_A = ac.W_TO - ac.W_empty - ac.W_payload
        Wi_A = ac.W_TO
        Wf_A = Wi_A - W_fuel_A
        CL_A = Wi_A / (q * S)
        CD_A = CD0 + K * CL_A**2
        LD_A = CL_A / CD_A
        R_A = breguet_range_jet_nmi(fps_to_kts(V), C_hr, LD_A, Wi_A, Wf_A)
        PL_A = ac.W_payload

        # Point B: max payload, max fuel (MTOW)
        # Same as A if fuel is limited by MTOW
        W_fuel_B = ac.W_fuel_max
        Wi_B = ac.W_empty + ac.W_fuel_max + ac.W_payload
        if Wi_B > ac.W_TO:
            Wi_B = ac.W_TO
            W_fuel_B = ac.W_TO - ac.W_empty - ac.W_payload

        PL_B = ac.W_payload
        Wf_B = Wi_B - W_fuel_B
        CL_B = Wi_B / (q * S)
        CD_B = CD0 + K * CL_B**2
        LD_B = CL_B / CD_B
        R_B = breguet_range_jet_nmi(fps_to_kts(V), C_hr, LD_B, Wi_B, Wf_B)

        # Point C: max fuel, reduced payload (trade payload for fuel)
        W_fuel_C = ac.W_fuel_max
        PL_C = ac.W_TO - ac.W_empty - W_fuel_C
        PL_C = max(PL_C, 0)
        Wi_C = ac.W_empty + W_fuel_C + PL_C
        Wf_C = Wi_C - W_fuel_C
        CL_C = Wi_C / (q * S)
        CD_C = CD0 + K * CL_C**2
        LD_C = CL_C / CD_C
        R_C = breguet_range_jet_nmi(fps_to_kts(V), C_hr, LD_C, Wi_C, Wf_C)

        # Point D: ferry range (zero payload, max fuel)
        W_fuel_D = ac.W_fuel_max
        Wi_D = ac.W_empty + W_fuel_D
        Wf_D = ac.W_empty
        CL_D = Wi_D / (q * S)
        CD_D = CD0 + K * CL_D**2
        LD_D = CL_D / CD_D
        R_D = breguet_range_jet_nmi(fps_to_kts(V), C_hr, LD_D, Wi_D, Wf_D)

        ranges  = [0, R_B, R_C, R_D]
        payloads = [PL_B, PL_B, PL_C, 0]

        ax.plot(ranges, np.array(payloads) / 1000, '-o', color=COLORS[ac.name],
                label=ac.name, lw=2, ms=6)

        # Annotate points
        ax.annotate(f"  A/B", (R_B, PL_B / 1000), fontsize=8)
        ax.annotate(f"  C", (R_C, PL_C / 1000), fontsize=8)
        ax.annotate(f"  D (ferry)", (R_D, 0), fontsize=8)

    ax.set_xlabel("Range [nmi]")
    ax.set_ylabel("Payload [1000 lb]")
    ax.set_title("Payload-Range Diagram  (Breguet, Raymer Eq. 17.23)")
    ax.legend()
    ax.grid(True, alpha=0.3)
    ax.set_ylim(bottom=-0.5)
    save(fig, "07_payload_range.png")


# =====================================================================
# 8. Takeoff & Landing Bar Chart
# =====================================================================
def plot_TO_landing_bars():
    from perf import total_takeoff_distance, total_landing_distance

    fig, axes = plt.subplots(1, 2, figsize=(12, 5))

    names = []
    to_segments = {"Ground Roll": [], "Rotation": [], "Transition": [], "Climb": []}
    la_segments = {"Approach": [], "Flare": [], "Free Roll": [], "Braking": []}

    for ac in VARIANTS:
        W, S, CD0, K = ac.W_TO, ac.S, ac.CD0, ac.K
        rho = RHO_SL
        T_TO = 0.75 * (5.0 + ac.BPR) / (4.0 + ac.BPR) * ac.T_max_SL

        to = total_takeoff_distance(W, S, T_TO, CD0+0.015, ac.CL_ground, K,
                                    ac.mu_roll, rho, ac.CL_max_TO,
                                    ac.thrust_to_weight, ac.h_obstacle_TO, 3.0)

        la = total_landing_distance(ac.W_landing, S, rho, ac.CL_max_L,
                                    CD0+0.02, ac.CL_ground, K, ac.mu_brake,
                                    ac.h_obstacle_L, 0.05*ac.T_max_SL, 0.0,
                                    1.3, 3.0, False)

        names.append(ac.name)
        to_segments["Ground Roll"].append(to["S_ground_roll"])
        to_segments["Rotation"].append(to["S_rotation"])
        to_segments["Transition"].append(to["S_transition"])
        to_segments["Climb"].append(to["S_climb"])

        la_segments["Approach"].append(la["S_approach"])
        la_segments["Flare"].append(la["S_flare"])
        la_segments["Free Roll"].append(la["S_free_roll"])
        la_segments["Braking"].append(la["S_braking"])

    x = np.arange(len(names))
    width = 0.5

    # Takeoff
    ax = axes[0]
    bottom = np.zeros(len(names))
    for seg, vals in to_segments.items():
        ax.bar(x, vals, width, bottom=bottom, label=seg)
        bottom += np.array(vals)
    ax.set_xticks(x)
    ax.set_xticklabels(names)
    ax.set_ylabel("Distance [ft]")
    ax.set_title("Takeoff Distance Breakdown\n(Raymer 17.8)")
    ax.legend(fontsize=8)
    ax.grid(True, alpha=0.3, axis='y')

    # Landing
    ax = axes[1]
    bottom = np.zeros(len(names))
    for seg, vals in la_segments.items():
        ax.bar(x, vals, width, bottom=bottom, label=seg)
        bottom += np.array(vals)
    ax.set_xticks(x)
    ax.set_xticklabels(names)
    ax.set_ylabel("Distance [ft]")
    ax.set_title("Landing Distance Breakdown\n(Raymer 17.9)")
    ax.legend(fontsize=8)
    ax.grid(True, alpha=0.3, axis='y')

    fig.tight_layout()
    save(fig, "08_TO_landing_breakdown.png")


# =====================================================================
# 9. ASDR & TODR vs Engine Failure Speed
# =====================================================================
def plot_asdr_todr():
    fig, axes = plt.subplots(1, len(VARIANTS), figsize=(14, 6), sharey=True)
    if len(VARIANTS) == 1:
        axes = [axes]

    for ax, ac in zip(axes, VARIANTS):
        W, S, CD0, K = ac.W_TO, ac.S, ac.CD0, ac.K
        rho = RHO_SL
        CD0_TO = CD0 + 0.015
        T_TO = 0.75 * (5.0 + ac.BPR) / (4.0 + ac.BPR) * ac.T_max_SL
        T_idle = 0.05 * ac.T_max_SL

        curves = asdr_todr_curves(
            W=W, S=S, T=T_TO, CD0=CD0_TO,
            CL_ground=ac.CL_ground, K=K,
            mu_roll=ac.mu_roll, mu_brake=ac.mu_brake,
            rho=rho, CL_max_TO=ac.CL_max_TO, TW=ac.thrust_to_weight,
            h_obstacle=ac.h_obstacle_TO, n_engines=ac.n_engines,
            T_idle=T_idle, T_reverse=0.0, t_react=2.0, t_rotate=3.0)

        v1_result = find_V1(
            W=W, S=S, T=T_TO, CD0=CD0_TO,
            CL_ground=ac.CL_ground, K=K,
            mu_roll=ac.mu_roll, mu_brake=ac.mu_brake,
            rho=rho, CL_max_TO=ac.CL_max_TO, TW=ac.thrust_to_weight,
            h_obstacle=ac.h_obstacle_TO, n_engines=ac.n_engines,
            T_idle=T_idle, T_reverse=0.0, t_react=2.0, t_rotate=3.0)

        V_kts = fps_to_kts(curves["V_EF"])
        V1_kts = fps_to_kts(v1_result["V1"])
        BFL = v1_result["BFL"]

        ax.plot(V_kts, curves["ASDR"], 'r-', lw=2, label="ASDR (accelerate-stop)")
        ax.plot(V_kts, curves["AGDR"], 'b-', lw=2, label="AGDR (accelerate-go)")

        # Mark V1 / BFL intersection
        ax.plot(V1_kts, BFL, 'ko', ms=10, zorder=5)
        ax.annotate(f"V1 = {V1_kts:.0f} kts\nBFL = {BFL:,.0f} ft",
                    (V1_kts, BFL), fontsize=9, fontweight='bold',
                    xytext=(-80, 30), textcoords='offset points',
                    arrowprops=dict(arrowstyle='->', color='black'))

        # Mark VR
        V_s = V_stall(W, S, rho, ac.CL_max_TO)
        VR_kts = fps_to_kts(1.1 * V_s)
        ax.axvline(VR_kts, color='green', ls='--', alpha=0.6, label=f"VR = {VR_kts:.0f} kts")

        ax.set_xlabel("Engine Failure Speed, $V_{EF}$ [kts]")
        ax.set_title(f"{ac.name} (Sea Level, ISA)")
        ax.legend(fontsize=8, loc="upper left")
        ax.grid(True, alpha=0.3)

    axes[0].set_ylabel("Distance [ft]")
    fig.suptitle("ASDR & AGDR vs Engine Failure Speed — Balanced Field Length",
                 fontsize=12)
    save(fig, "09_ASDR_TODR.png")


# =====================================================================
# 10. Glide Polar (Sink Rate vs Velocity)
# =====================================================================
def plot_glide_polar():
    fig, ax = plt.subplots(figsize=(8, 6))

    for ac in VARIANTS:
        W, S, CD0, K = ac.W_TO, ac.S, ac.CD0, ac.K
        h = ac.h_cruise_ft
        rho = isa_density(h)

        V_arr = np.linspace(kts_to_fps(100), kts_to_fps(350), 200)
        q_arr = 0.5 * rho * V_arr**2
        CL_arr = W / (q_arr * S)
        CD_arr = CD0 + K * CL_arr**2
        sr_arr = sink_rate(W, S, rho, CL_arr, CD_arr)

        ax.plot(fps_to_kts(V_arr), sr_arr * 60, color=COLORS[ac.name],
                label=ac.name, lw=2)

        # Min sink
        sr_min = min_sink_rate(W, S, rho, CD0, K)
        from perf import V_min_sink
        v_ms = V_min_sink(W, S, rho, CD0, K)
        ax.plot(fps_to_kts(v_ms), sr_min * 60, 'o', color=COLORS[ac.name], ms=8)
        ax.annotate(f"  min sink\n  {sr_min*60:.0f} fpm",
                    (fps_to_kts(v_ms), sr_min * 60), fontsize=8)

        # Best glide
        v_bg = V_best_glide(W, S, rho, CD0, K)
        cl_bg = W / (0.5 * rho * v_bg**2 * S)
        cd_bg = CD0 + K * cl_bg**2
        sr_bg = sink_rate(W, S, rho, cl_bg, cd_bg)
        ax.plot(fps_to_kts(v_bg), sr_bg * 60, 's', color=COLORS[ac.name], ms=8)
        ax.annotate(f"  best glide\n  {fps_to_kts(v_bg):.0f} kts",
                    (fps_to_kts(v_bg), sr_bg * 60), fontsize=8)

    ax.set_xlabel("True Airspeed [kts]")
    ax.set_ylabel("Sink Rate [fpm]")
    ax.set_title(f"Glide Polar at {ac.h_cruise_ft/1000:.0f}k ft  (Raymer Eq. 17.68)")
    ax.invert_yaxis()  # convention: sink is positive downward
    ax.legend()
    ax.grid(True, alpha=0.3)
    save(fig, "10_glide_polar.png")


# =====================================================================
# 11. Mission Profile (Weight vs Segment)
# =====================================================================
def plot_mission_profile():
    # Only plot for variants with design_range_nm set
    variants_with_range = [ac for ac in VARIANTS if ac.design_range_nm > 0]
    if not variants_with_range:
        return

    fig, axes = plt.subplots(1, 2, figsize=(15, 6))

    # Left: weight through mission
    ax = axes[0]
    for ac in variants_with_range:
        mp = mission_profile(ac)
        segs = mp["segments"]
        names = [s["name"].split(". ", 1)[-1] for s in segs]
        W_trace = [s["W_start"] for s in segs] + [segs[-1]["W_end"]]
        x = np.arange(len(W_trace))
        labels = ["Start"] + names

        ax.step(x, np.array(W_trace) / 1000, where="post",
                color=COLORS[ac.name], lw=2, label=ac.name)
        # Mark trip fuel boundary (after segment 4)
        ax.axvline(5, color='gray', ls=':', alpha=0.4)

    ax.set_xticks(np.arange(len(labels)))
    ax.set_xticklabels(labels, rotation=45, ha="right", fontsize=7)
    ax.set_ylabel("Weight [1000 lb]")
    ax.set_title("Mission Profile — Weight vs Segment")
    ax.legend(fontsize=8)
    ax.grid(True, alpha=0.3)
    ax.annotate("trip fuel | reserves", xy=(5, ax.get_ylim()[0]),
                fontsize=7, color='gray', ha='center', va='bottom')

    # Right: fuel per segment (stacked bar)
    ax = axes[1]
    seg_names = [s["name"].split(". ", 1)[-1]
                 for s in mission_profile(variants_with_range[0])["segments"]]
    x = np.arange(len(seg_names))
    width = 0.8 / len(variants_with_range)

    for i, ac in enumerate(variants_with_range):
        mp = mission_profile(ac)
        fuels = [s["fuel"] for s in mp["segments"]]
        ax.bar(x + i * width, fuels, width, color=COLORS[ac.name], label=ac.name)

    ax.set_xticks(x + width * (len(variants_with_range) - 1) / 2)
    ax.set_xticklabels(seg_names, rotation=45, ha="right", fontsize=7)
    ax.set_ylabel("Fuel Burned [lb]")
    ax.set_title("Fuel Burn per Mission Segment")
    ax.legend(fontsize=8)
    ax.grid(True, alpha=0.3, axis='y')

    fig.tight_layout()
    save(fig, "11_mission_profile.png")


# =====================================================================
# 12. FAR 25 OEI Climb Gradients
# =====================================================================
def plot_oei_gradients():
    variants_with_range = [ac for ac in VARIANTS if ac.design_range_nm > 0]
    if not variants_with_range:
        return

    fig, ax = plt.subplots(figsize=(10, 5))

    categories = ["2nd segment\n(TO flaps)", "En-route\n(clean)", "Approach\n(landing)"]
    keys = ["2nd_segment", "en_route", "approach_climb"]
    minimums = [2.4, 1.2, 2.1]
    x = np.arange(len(categories))
    width = 0.8 / len(variants_with_range)

    for i, ac in enumerate(variants_with_range):
        mp = mission_profile(ac)
        vals = [mp["oei"][k]["gradient_pct"] for k in keys]
        ax.bar(x + i * width, vals, width, color=COLORS[ac.name], label=ac.name)

    # Minimum lines
    for j, (cat_x, mn) in enumerate(zip(x, minimums)):
        ax.plot([cat_x - 0.1, cat_x + width * len(variants_with_range) + 0.1],
                [mn, mn], 'r--', lw=1.5, alpha=0.7)
        ax.annotate(f"min {mn}%", (cat_x + width * len(variants_with_range), mn),
                    fontsize=8, color='red', va='bottom')

    ax.set_xticks(x + width * (len(variants_with_range) - 1) / 2)
    ax.set_xticklabels(categories)
    ax.set_ylabel("Climb Gradient [%]")
    ax.set_title("FAR 25 OEI Climb Gradients  (2-engine)")
    ax.legend(fontsize=8)
    ax.grid(True, alpha=0.3, axis='y')
    fig.tight_layout()
    save(fig, "12_oei_gradients.png")


# =====================================================================
ALL_PLOTS = [
    plot_thrust_vs_velocity,
    plot_power_vs_velocity,
    plot_LD_vs_CL,
    plot_ROC_vs_altitude,
    plot_Ps_vs_Mach,
    plot_turn_performance,
    plot_payload_range,
    plot_TO_landing_bars,
    plot_asdr_todr,
    plot_glide_polar,
    plot_mission_profile,
    plot_oei_gradients,
]

if __name__ == "__main__":
    import __main__
    base = os.path.join(os.path.dirname(__file__), "charts")

    for family_name, family in FAMILIES.items():
        __main__.VARIANTS = family["variants"]
        __main__.COLORS = family["colors"]
        __main__.CHART_DIR = os.path.join(base, family_name)
        os.makedirs(CHART_DIR, exist_ok=True)

        print(f"\n{'=' * 50}")
        print(f"  Generating {family_name} charts...")
        print(f"{'=' * 50}")
        for fn in ALL_PLOTS:
            fn()

    print(f"\nAll charts saved to examples/charts/CRJ/ and examples/charts/ZRJ/")
