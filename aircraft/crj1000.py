"""
Bombardier CRJ-1000 (CRJ-900NextGen stretch) performance data.

Sources
-------
- Bombardier CRJ-1000 Airport Planning Manual
- Jane's All the World's Aircraft
- Type Certificate Data Sheet A21EA
"""

from .base import AircraftData

crj1000 = AircraftData(
    name="CRJ-1000",

    # ---- Weights ----
    W_TO=91_800.0,          # lb, MTOW
    W_empty=53_250.0,       # lb, OEW (typical)
    W_fuel_max=19_594.0,    # lb, max usable fuel
    W_payload=23_100.0,     # lb, typical 100 pax @ 200 lb + bags

    # ---- Wing geometry ----
    S=568.0,                # ft^2, reference area (extended wing)
    b=85.9,                 # ft, span
    AR=13.0,                # aspect ratio
    sweep_qc_deg=25.0,      # deg, quarter-chord sweep
    t_c=0.11,               # thickness-to-chord ratio
    taper=0.32,             # taper ratio

    # ---- Aerodynamics ----
    CD0=0.024,              # zero-lift drag (cruise, slightly higher due to stretch)
    K=0.036,                # induced drag factor
    e=0.82,                 # Oswald span efficiency (higher AR helps)
    CL_max_clean=1.40,
    CL_max_TO=1.80,
    CL_max_L=2.40,

    # ---- Propulsion (2x GE CF34-8C5A1) ----
    n_engines=2,
    T_max_SL=29_020.0,      # lb total (2 x 14,510 lb)
    TSFC=0.69,              # 1/hr (cruise estimate)
    BPR=5.0,

    # ---- Landing gear ----
    mu_roll=0.03,
    mu_brake=0.40,

    # ---- Takeoff / landing ----
    h_obstacle_TO=35.0,
    h_obstacle_L=50.0,
    CL_ground=0.10,

    # ---- Design cruise ----
    M_cruise=0.78,
    h_cruise_ft=35_000.0,

    # ---- Structural limits ----
    n_max=2.5,
    V_NE_kts=335.0,
)
