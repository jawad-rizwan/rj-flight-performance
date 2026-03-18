"""
Bombardier CRJ-700 (CRJ-705 baseline) performance data.

Sources
-------
- Bombardier CRJ-700 Airport Planning Manual
- Jane's All the World's Aircraft
- Type Certificate Data Sheet A21EA
"""

from .base import AircraftData

crj700 = AircraftData(
    name="CRJ-700",

    # ---- Weights ----
    W_TO=75_000.0,          # lb, MTOW
    W_empty=44_245.0,       # lb, OEW (typical)
    W_fuel_max=19_594.0,    # lb, max usable fuel
    W_payload=17_500.0,     # lb, typical 70 pax @ 200 lb + bags

    # ---- Wing geometry ----
    S=520.0,                # ft^2, reference area
    b=76.3,                 # ft, span
    AR=11.2,                # aspect ratio
    sweep_qc_deg=25.0,      # deg, quarter-chord sweep
    t_c=0.11,               # thickness-to-chord ratio (root avg)
    taper=0.32,             # taper ratio

    # ---- Aerodynamics ----
    CD0=0.022,              # zero-lift drag (cruise estimate)
    K=0.040,                # induced drag factor  (1/(pi*e*AR))
    e=0.80,                 # Oswald span efficiency
    CL_max_clean=1.40,      # clean
    CL_max_TO=2.20,         # takeoff flaps + slats
    CL_max_L=2.40,          # landing flaps + slats

    # ---- Propulsion (2x GE CF34-8C5B1) ----
    n_engines=2,
    T_max_SL=26_760.0,      # lb total (2 x 13,360 lb)
    TSFC=0.68,              # 1/hr (cruise estimate)
    BPR=5.0,

    # ---- Landing gear ----
    mu_roll=0.03,
    mu_brake=0.50,

    # ---- Takeoff / landing ----
    h_obstacle_TO=35.0,     # ft (FAR 25)
    h_obstacle_L=50.0,      # ft
    CL_ground=0.10,

    # ---- Design cruise ----
    M_cruise=0.78,
    h_cruise_ft=35_000.0,   # ft

    # ---- Structural limits ----
    n_max=2.5,
    V_NE_kts=335.0,         # KEAS
    M_mo=0.85,
)
