"""
ZRJ100 — 96-seat high-wing regional jet (EU variant).

Geometry and weights from rj-basic-aerodynamics and rj-mission-analysis.
Aero coefficients synced via sync_aero.py from rj-basic-aerodynamics.
"""

from .base import AircraftData

zrj100 = AircraftData(
    name="ZRJ100",

    # ---- Weights ----
    W_TO=89_414.0,          # lb, MTOW (Ch.17 mission analysis, shared wing fuel)
    W_empty=48_341.0,       # lb, OEW (We 47553 + crew 788)
    W_fuel_max=17_693.0,    # lb, total mission fuel (shared wing with ZRJ70)
    W_payload=23_380.0,     # lb, 96 pax + cargo

    # ---- Wing geometry ----
    S=1016.58,              # ft^2, reference area (trapezoidal)
    b=89.05,                # ft, span
    AR=7.8,                 # aspect ratio
    sweep_qc_deg=22.9,      # deg, quarter-chord sweep
    t_c=0.123,              # chord-weighted avg t/c
    taper=0.33,             # taper ratio

    # ---- Aerodynamics (synced from rj-basic-aerodynamics) ----
    CD0=0.01910,            # zero-lift drag coefficient
    K=0.05418,              # induced drag factor (LE suction method)
    e=0.753,                # Oswald span efficiency (LE suction)
    CL_max_clean=1.16,      # clean
    CL_max_TO=1.99,         # takeoff flaps
    CL_max_L=2.62,          # landing flaps + slats

    # ---- Propulsion (2x PW1200G) ----
    n_engines=2,
    T_max_SL=38_380.0,      # lb total (2 x 19,190 lb)
    TSFC=0.50,              # 1/hr (cruise, refined GTF estimate)
    BPR=9.0,

    # ---- Design cruise ----
    M_cruise=0.78,
    h_cruise_ft=35_000.0,   # ft

    # ---- Structural limits ----
    n_max=2.5,
    V_NE_kts=340.0,         # KEAS (estimate for Mmo 0.85)
)
