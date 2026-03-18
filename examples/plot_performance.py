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
  7. Takeoff & landing distance breakdown
  8. Glide polar (sink rate vs velocity)
  9. Operating envelope (altitude vs Mach)
 10. Climb hodograph (Vx vs Vy)
 11. Ps contour plot (altitude vs Mach)
 12. Airfield performance by runway surface
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
    ROC_jet, V_best_ROC_jet, climb_angle, rate_of_climb,
    turn_rate_deg, turn_radius, corner_speed,
    n_sustained, sustained_turn_envelope_rho,
    Ps_expanded, specific_energy,
    V_best_glide, min_sink_rate, sink_rate,
    asdr_todr_curves,
    mdd_wing, cd0_at_mach,
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
        alts = np.linspace(0, 50000, 100)
        rocs = []
        for h in alts:
            rho_h = isa_density(h)
            a_h = speed_of_sound(h)
            T_h = thrust_at_altitude(ac.T_max_SL, h, ac.BPR)
            TW_h = T_h / W
            V_h = V_best_ROC_jet(W, S, rho_h, CD0, K, TW_h)
            # Cap at Mmo
            if ac.M_mo > 0:
                V_h = min(V_h, ac.M_mo * a_h)
            # Wave drag
            M_h = V_h / a_h
            q_h = 0.5 * rho_h * V_h**2
            CL_h = W / (q_h * S) if q_h > 0 else 0
            mdd_h = mdd_wing(ac.t_c, ac.sweep_qc_deg, CL_h)
            CD0_h = cd0_at_mach(CD0, M_h, mdd_h)
            roc = ROC_jet(V_h, T_h, CD0_h, K, W, S, rho_h)
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
            # Wave drag: adjust CD0 for each Mach
            q_arr = 0.5 * rho * V_arr**2
            CL_arr = W / (q_arr * S)
            mdd_arr = mdd_wing(ac.t_c, ac.sweep_qc_deg, CL_arr)
            CD0_arr = cd0_at_mach(CD0, M_arr, mdd_arr)
            ps_arr = Ps_expanded(V_arr, TW, ac.wing_loading, CD0_arr, K, rho, n=1.0)

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
# 7. Takeoff & Landing Bar Chart
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
    save(fig, "07_TO_landing_breakdown.png")


# =====================================================================
# 8. Glide Polar (Sink Rate vs Velocity)
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
    save(fig, "08_glide_polar.png")


# =====================================================================
# 9. Operating Envelope (Altitude vs Mach)
# =====================================================================
def plot_operating_envelope():
    fig, ax = plt.subplots(figsize=(10, 7))

    for ac in VARIANTS:
        W, S, CD0, K = ac.W_TO, ac.S, ac.CD0, ac.K

        # Sweep altitudes
        alts = np.linspace(0, 55000, 200)
        M_stall = []      # low-speed stall boundary
        M_max_thrust = []  # high-speed thrust boundary (Ps=0)

        for h in alts:
            rho_h = isa_density(h)
            a_h = speed_of_sound(h)
            T_h = thrust_at_altitude(ac.T_max_SL, h, ac.BPR)

            # Stall boundary: 1g stall Mach at this altitude
            Vs = V_stall(W, S, rho_h, ac.CL_max_clean)
            M_stall.append(Vs / a_h)

            # Thrust boundary: find max Mach where T >= D (Ps >= 0)
            M_test = np.linspace(0.3, 0.95, 300)
            V_test = M_test * a_h
            # Wave drag: adjust CD0 for each Mach
            q_test = 0.5 * rho_h * V_test**2
            CL_test = W / (q_test * S)
            mdd_test = mdd_wing(ac.t_c, ac.sweep_qc_deg, CL_test)
            CD0_test = cd0_at_mach(CD0, M_test, mdd_test)
            T_req = np.array([thrust_required(W, cd0_m, K, S, rho_h, v)
                              for cd0_m, v in zip(CD0_test, V_test)])
            feasible = T_req <= T_h
            # Also cap at Mmo
            if ac.M_mo > 0:
                feasible = feasible & (M_test <= ac.M_mo)
            if np.any(feasible):
                M_max_thrust.append(M_test[feasible][-1])
            else:
                M_max_thrust.append(np.nan)

        M_stall = np.array(M_stall)
        M_max_thrust = np.array(M_max_thrust)
        alts_k = alts / 1000

        # Plot stall boundary (left side)
        ax.plot(M_stall, alts_k, '--', color=COLORS[ac.name], alpha=0.7, lw=1.5)
        # Plot thrust boundary (right side)
        ax.plot(M_max_thrust, alts_k, '-', color=COLORS[ac.name], lw=2,
                label=ac.name)
        # Fill the envelope
        ax.fill_betweenx(alts_k, M_stall, M_max_thrust,
                         color=COLORS[ac.name], alpha=0.08)

    # Mmo limit
    ax.axvline(0.85, color='red', ls=':', lw=1.5, alpha=0.7, label="$M_{mo}$ = 0.85")

    # Cruise point
    ax.axhline(35, color='gray', ls=':', alpha=0.4)
    ax.plot(0.78, 35, 'k*', ms=12, zorder=5, label="Design cruise (M 0.78, FL350)")

    ax.set_xlabel("Mach Number")
    ax.set_ylabel("Altitude [1000 ft]")
    ax.set_title("Operating Envelope — Altitude vs Mach")
    ax.set_xlim(0.15, 0.95)
    ax.set_ylim(0, 55)
    ax.legend(fontsize=8, loc="upper left")
    ax.grid(True, alpha=0.3)
    save(fig, "09_operating_envelope.png")


# =====================================================================
# 10. Climb Hodograph (Raymer 17.3.2)
# =====================================================================
def plot_climb_hodograph():
    fig, axes = plt.subplots(1, len(VARIANTS), figsize=(14, 5), sharey=True)
    if len(VARIANTS) == 1:
        axes = [axes]

    for ax, ac in zip(axes, VARIANTS):
        W, S, CD0, K = ac.W_TO, ac.S, ac.CD0, ac.K
        rho = RHO_SL
        T = ac.T_max_SL

        Vs = V_stall(W, S, rho, ac.CL_max_clean)
        V_arr = np.linspace(Vs, kts_to_fps(450), 300)
        T_req = thrust_required(W, CD0, K, S, rho, V_arr)
        gamma_arr = np.arcsin(np.clip((T - T_req) / W, -1, 1))
        Vx = fps_to_kts(V_arr * np.cos(gamma_arr))  # horizontal
        Vy = V_arr * np.sin(gamma_arr) * 60           # vertical [fpm]

        ax.plot(Vx, Vy, color=COLORS[ac.name], lw=2)

        # Best angle of climb (max gamma)
        idx_gamma = np.argmax(gamma_arr)
        ax.plot(Vx[idx_gamma], Vy[idx_gamma], 'o', color=COLORS[ac.name], ms=8)
        ax.annotate(f"  best angle\n  {np.degrees(gamma_arr[idx_gamma]):.1f} deg",
                    (Vx[idx_gamma], Vy[idx_gamma]), fontsize=8)

        # Best rate of climb (max Vy)
        idx_roc = np.argmax(Vy)
        ax.plot(Vx[idx_roc], Vy[idx_roc], 's', color=COLORS[ac.name], ms=8)
        ax.annotate(f"  best ROC\n  {Vy[idx_roc]:,.0f} fpm",
                    (Vx[idx_roc], Vy[idx_roc]), fontsize=8)

        ax.axhline(0, color='k', lw=0.5)
        ax.set_xlabel("Horizontal Speed [kts]")
        ax.set_title(f"{ac.name} (Sea Level)")
        ax.grid(True, alpha=0.3)

    axes[0].set_ylabel("Rate of Climb [fpm]")
    fig.suptitle("Climb Hodograph  (Raymer Sec. 17.3.2)", fontsize=12)
    save(fig, "10_climb_hodograph.png")


# =====================================================================
# 11. Ps Contour Plot — Altitude vs Mach (Raymer 17.6.2)
# =====================================================================
def plot_Ps_contours():
    fig, axes = plt.subplots(1, len(VARIANTS), figsize=(14, 5), sharey=True)
    if len(VARIANTS) == 1:
        axes = [axes]

    for ax, ac in zip(axes, VARIANTS):
        W, S, CD0, K = ac.W_TO, ac.S, ac.CD0, ac.K

        M_arr = np.linspace(0.15, 0.90, 150)
        h_arr = np.linspace(0, 50000, 150)
        M_grid, H_grid = np.meshgrid(M_arr, h_arr)
        Ps_grid = np.zeros_like(M_grid)

        for i in range(len(h_arr)):
            h = h_arr[i]
            rho_h = isa_density(h)
            a_h = speed_of_sound(h)
            T_h = thrust_at_altitude(ac.T_max_SL, h, ac.BPR)
            TW_h = T_h / W
            V_row = M_arr * a_h
            # Wave drag: adjust CD0 for each Mach
            q_row = 0.5 * rho_h * V_row**2
            CL_row = np.where(q_row > 0, W / (q_row * S), 0)
            mdd_row = mdd_wing(ac.t_c, ac.sweep_qc_deg, CL_row)
            CD0_row = cd0_at_mach(CD0, M_arr, mdd_row)
            Ps_grid[i, :] = Ps_expanded(V_row, TW_h, ac.wing_loading,
                                         CD0_row, K, rho_h, n=1.0) * 60  # fpm

        levels = [0, 500, 1000, 2000, 5000, 10000, 15000]
        cs = ax.contour(M_grid, H_grid / 1000, Ps_grid, levels=levels,
                        colors='k', linewidths=0.8)
        ax.clabel(cs, inline=True, fontsize=7, fmt='%g fpm')
        cf = ax.contourf(M_grid, H_grid / 1000, Ps_grid, levels=levels,
                         cmap='RdYlGn', alpha=0.4)

        # Ps = 0 contour (flight envelope boundary)
        ax.contour(M_grid, H_grid / 1000, Ps_grid, levels=[0],
                   colors='red', linewidths=2)

        # Cruise point
        ax.plot(ac.M_cruise, ac.h_cruise_ft / 1000, 'k*', ms=10, zorder=5)

        ax.set_xlabel("Mach Number")
        ax.set_title(f"{ac.name}")
        ax.set_xlim(0.15, 0.90)
        ax.set_ylim(0, 50)
        ax.grid(True, alpha=0.2)

    axes[0].set_ylabel("Altitude [1000 ft]")
    fig.suptitle("Specific Excess Power Contours  (Raymer Fig. 17.9, Eq. 17.89)", fontsize=12)
    save(fig, "11_Ps_contours.png")


# =====================================================================
# 12. Airfield Performance by Runway Surface (Raymer Table 17.1)
# =====================================================================
def plot_surface_performance():
    from perf import total_takeoff_distance, total_landing_distance, find_V1

    SURFACES = [
        ("Dry\nconcrete",  0.03, 0.40),
        ("Wet\nconcrete",  0.05, 0.225),
        ("Icy\nconcrete",  0.02, 0.08),
        ("Hard\nturf",     0.05, 0.40),
        ("Firm\ndirt",     0.04, 0.30),
        ("Soft\nturf",     0.07, 0.20),
        ("Wet\ngrass",     0.08, 0.20),
    ]

    fig, axes = plt.subplots(1, 2, figsize=(15, 6))
    x = np.arange(len(SURFACES))
    width = 0.8 / len(VARIANTS)

    for i, ac in enumerate(VARIANTS):
        W, S, CD0, K = ac.W_TO, ac.S, ac.CD0, ac.K
        rho = RHO_SL
        CD0_TO = CD0 + 0.015
        T_TO = 0.75 * (5.0 + ac.BPR) / (4.0 + ac.BPR) * ac.T_max_SL
        T_idle = 0.05 * ac.T_max_SL

        tofls = []
        ldrs = []
        for _, mu_r, mu_b in SURFACES:
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
            tofls.append(max(to['TODR_factored'], v1['BFL']))

            la = total_landing_distance(
                W=ac.W_landing, S=S, rho=rho, CL_max_L=ac.CL_max_L,
                CD0=CD0 + 0.02, CL_ground=ac.CL_ground, K=K,
                mu_brake=mu_b, h_obstacle=ac.h_obstacle_L,
                T_idle=T_idle, T_reverse=0.0,
                approach_factor=1.3, t_free=3.0, FAR_factor=True)
            ldrs.append(la['S_FAR_field'])

        axes[0].bar(x + i * width, tofls, width, color=COLORS[ac.name], label=ac.name)
        axes[1].bar(x + i * width, ldrs, width, color=COLORS[ac.name], label=ac.name)

    surf_labels = [s[0] for s in SURFACES]
    for ax, title in zip(axes, ["TOFL (FAR 25)", "Landing Field Length (FAR 25)"]):
        ax.set_xticks(x + width * (len(VARIANTS) - 1) / 2)
        ax.set_xticklabels(surf_labels, fontsize=8)
        ax.set_ylabel("Distance [ft]")
        ax.set_title(title)
        ax.legend(fontsize=8)
        ax.grid(True, alpha=0.3, axis='y')

    fig.suptitle("Airfield Performance by Runway Surface  (Raymer Table 17.1)", fontsize=12)
    fig.tight_layout()
    save(fig, "12_surface_performance.png")


# =====================================================================
# 13. Per-Surface Takeoff/Landing Breakdown + BFL Charts
# =====================================================================
def plot_per_surface_charts():
    from perf import total_takeoff_distance, total_landing_distance, find_V1

    SURFACES = [
        ("Dry concrete",  "dry_concrete",  0.03, 0.40),
        ("Wet concrete",  "wet_concrete",  0.05, 0.225),
        ("Icy concrete",  "icy_concrete",  0.02, 0.08),
        ("Hard turf",     "hard_turf",     0.05, 0.40),
        ("Firm dirt",     "firm_dirt",      0.04, 0.30),
        ("Soft turf",     "soft_turf",     0.07, 0.20),
        ("Wet grass",     "wet_grass",     0.08, 0.20),
    ]

    surf_dir = os.path.join(CHART_DIR, "surfaces")
    os.makedirs(surf_dir, exist_ok=True)

    for surf_name, surf_key, mu_r, mu_b in SURFACES:
        # ---- TO/Landing breakdown ----
        fig, axes = plt.subplots(1, 2, figsize=(12, 5))
        names = []
        to_segs = {"Ground Roll": [], "Rotation": [], "Transition": [], "Climb": []}
        la_segs = {"Approach": [], "Flare": [], "Free Roll": [], "Braking": []}

        for ac in VARIANTS:
            W, S, CD0, K = ac.W_TO, ac.S, ac.CD0, ac.K
            rho = RHO_SL
            T_TO = 0.75 * (5.0 + ac.BPR) / (4.0 + ac.BPR) * ac.T_max_SL
            T_idle = 0.05 * ac.T_max_SL

            to = total_takeoff_distance(W, S, T_TO, CD0+0.015, ac.CL_ground, K,
                                        mu_r, rho, ac.CL_max_TO,
                                        ac.thrust_to_weight, ac.h_obstacle_TO, 3.0)
            la = total_landing_distance(ac.W_landing, S, rho, ac.CL_max_L,
                                        CD0+0.02, ac.CL_ground, K, mu_b,
                                        ac.h_obstacle_L, T_idle, 0.0,
                                        1.3, 3.0, False)
            names.append(ac.name)
            to_segs["Ground Roll"].append(to["S_ground_roll"])
            to_segs["Rotation"].append(to["S_rotation"])
            to_segs["Transition"].append(to["S_transition"])
            to_segs["Climb"].append(to["S_climb"])
            la_segs["Approach"].append(la["S_approach"])
            la_segs["Flare"].append(la["S_flare"])
            la_segs["Free Roll"].append(la["S_free_roll"])
            la_segs["Braking"].append(la["S_braking"])

        x = np.arange(len(names))
        width = 0.5

        ax = axes[0]
        bottom = np.zeros(len(names))
        for seg, vals in to_segs.items():
            ax.bar(x, vals, width, bottom=bottom, label=seg)
            bottom += np.array(vals)
        ax.set_xticks(x); ax.set_xticklabels(names)
        ax.set_ylabel("Distance [ft]")
        ax.set_title(f"Takeoff Distance Breakdown")
        ax.legend(fontsize=8); ax.grid(True, alpha=0.3, axis='y')

        ax = axes[1]
        bottom = np.zeros(len(names))
        for seg, vals in la_segs.items():
            ax.bar(x, vals, width, bottom=bottom, label=seg)
            bottom += np.array(vals)
        ax.set_xticks(x); ax.set_xticklabels(names)
        ax.set_ylabel("Distance [ft]")
        ax.set_title(f"Landing Distance Breakdown")
        ax.legend(fontsize=8); ax.grid(True, alpha=0.3, axis='y')

        fig.suptitle(f"{surf_name}  (mu_roll={mu_r}, mu_brake={mu_b})", fontsize=12)
        fig.tight_layout()
        path = os.path.join(surf_dir, f"{surf_key}_TO_landing.png")
        fig.savefig(path, dpi=150, bbox_inches="tight"); plt.close(fig)
        print(f"  Saved: {path}")

        # ---- BFL chart ----
        fig, axes_bfl = plt.subplots(1, len(VARIANTS), figsize=(14, 6), sharey=True)
        if len(VARIANTS) == 1:
            axes_bfl = [axes_bfl]

        for ax, ac in zip(axes_bfl, VARIANTS):
            W, S, CD0, K = ac.W_TO, ac.S, ac.CD0, ac.K
            rho = RHO_SL
            CD0_TO = CD0 + 0.015
            T_TO = 0.75 * (5.0 + ac.BPR) / (4.0 + ac.BPR) * ac.T_max_SL
            T_idle = 0.05 * ac.T_max_SL

            curves = asdr_todr_curves(
                W=W, S=S, T=T_TO, CD0=CD0_TO,
                CL_ground=ac.CL_ground, K=K,
                mu_roll=mu_r, mu_brake=mu_b,
                rho=rho, CL_max_TO=ac.CL_max_TO, TW=ac.thrust_to_weight,
                h_obstacle=ac.h_obstacle_TO, n_engines=ac.n_engines,
                T_idle=T_idle, T_reverse=0.0, t_react=2.0, t_rotate=3.0)
            v1_result = find_V1(
                W=W, S=S, T=T_TO, CD0=CD0_TO,
                CL_ground=ac.CL_ground, K=K,
                mu_roll=mu_r, mu_brake=mu_b,
                rho=rho, CL_max_TO=ac.CL_max_TO, TW=ac.thrust_to_weight,
                h_obstacle=ac.h_obstacle_TO, n_engines=ac.n_engines,
                T_idle=T_idle, T_reverse=0.0, t_react=2.0, t_rotate=3.0)

            V_kts = fps_to_kts(curves["V_EF"])
            V1_kts = fps_to_kts(v1_result["V1"])
            BFL = v1_result["BFL"]

            ax.plot(V_kts, curves["ASDR"], 'r-', lw=2, label="ASDR")
            ax.plot(V_kts, curves["AGDR"], 'b-', lw=2, label="AGDR")
            ax.plot(V1_kts, BFL, 'ko', ms=10, zorder=5)
            ax.annotate(f"V1 = {V1_kts:.0f} kts\nBFL = {BFL:,.0f} ft",
                        (V1_kts, BFL), fontsize=9, fontweight='bold',
                        xytext=(-80, 30), textcoords='offset points',
                        arrowprops=dict(arrowstyle='->', color='black'))

            V_s = V_stall(W, S, rho, ac.CL_max_TO)
            VR_kts = fps_to_kts(1.1 * V_s)
            ax.axvline(VR_kts, color='green', ls='--', alpha=0.6,
                       label=f"VR = {VR_kts:.0f} kts")

            ax.set_xlabel("$V_{EF}$ [kts]")
            ax.set_title(f"{ac.name}")
            ax.legend(fontsize=8, loc="upper left")
            ax.grid(True, alpha=0.3)

        axes_bfl[0].set_ylabel("Distance [ft]")
        fig.suptitle(f"Balanced Field Length — {surf_name}  (mu_roll={mu_r}, mu_brake={mu_b})",
                     fontsize=12)
        path = os.path.join(surf_dir, f"{surf_key}_BFL.png")
        fig.savefig(path, dpi=150, bbox_inches="tight"); plt.close(fig)
        print(f"  Saved: {path}")


# =====================================================================
ALL_PLOTS = [
    plot_thrust_vs_velocity,
    plot_power_vs_velocity,
    plot_LD_vs_CL,
    plot_ROC_vs_altitude,
    plot_Ps_vs_Mach,
    plot_turn_performance,
    plot_TO_landing_bars,
    plot_glide_polar,
    plot_operating_envelope,
    plot_climb_hodograph,
    plot_Ps_contours,
    plot_surface_performance,
    plot_per_surface_charts,
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
