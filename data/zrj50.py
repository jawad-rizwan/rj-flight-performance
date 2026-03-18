"""
ZRJ50 — 50-seat high-wing regional jet (NA variant, scope-clause constrained).

Same airframe as ZRJ70, reduced payload, 65,000 lb MTOW cap.
Geometry and weights from rj-basic-aerodynamics and rj-mission-analysis.
Aero coefficients identical to ZRJ70 (same airframe).
"""

from .base import AircraftData

zrj50 = AircraftData(
    name="ZRJ50",

    # ---- Weights ----
    W_TO=65_000.0,          # lb, MTOW (scope clause cap)
    W_empty=45_578.0,       # lb, OEW (We 44987 + crew 591, 3 crew)
    W_fuel_max=8_072.0,    # lb, max fuel (MTOW-limited, tank cap 21258)
    W_payload=11_350.0,     # lb, 50 pax × 197 lbs + 50 bags × 30 lbs

    # ---- Wing geometry (same airframe as ZRJ70) ----
    S=1016.58,              # ft^2, reference area (trapezoidal)
    b=89.05,                # ft, span
    AR=7.8,                 # aspect ratio
    sweep_qc_deg=22.9,      # deg, quarter-chord sweep
    t_c=0.123,              # chord-weighted avg t/c
    taper=0.33,             # taper ratio

    # ---- Aerodynamics (same airframe as ZRJ70) ----
    CD0=0.01843,            # zero-lift drag coefficient
    K=0.05418,              # induced drag factor (LE suction method)
    e=0.753,                # Oswald span efficiency (LE suction)
    CL_max_clean=1.16,      # clean
    CL_max_TO=1.99,         # takeoff flaps
    CL_max_L=2.62,          # landing flaps + slats

    # ---- Propulsion (2x PW1200G) ----
    n_engines=2,
    T_max_SL=38_380.0,      # lb total (2 x 19,190 lb)
    TSFC=0.46,              # 1/hr (cruise, refined GTF estimate)
    BPR=9.0,

    # ---- Design cruise ----
    M_cruise=0.78,
    h_cruise_ft=35_000.0,   # ft

    # ---- Mission ----
    design_range_nm=1800.0,     # nm (same as ZRJ70)
    alternate_range_nm=100.0,   # nm

    # ---- Structural limits ----
    n_max=2.5,
    V_NE_kts=340.0,         # KEAS (estimate for Mmo 0.85)
)
