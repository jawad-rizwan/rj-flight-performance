"""
Aircraft data container for performance calculations.
All units are Imperial (ft, lb, slug, sec) consistent with Raymer.
"""

from dataclasses import dataclass, field
from typing import Optional


@dataclass
class AircraftData:
    """Container for aircraft parameters used in performance analysis.

    Modify the values in crj700.py / crj1000.py (or create a new file)
    to analyse a different aircraft.

    Units
    -----
    Weight      : lb (pounds-force)
    Length/span : ft
    Area        : ft^2
    Thrust      : lb (pounds-force)
    Power       : hp (only for props; jets use thrust)
    SFC         : 1/hr  (jet TSFC, lb-fuel/hr per lb-thrust)
    Speed       : ft/s  (or kts where noted)
    """

    # ---- Identity ----
    name: str = ""

    # ---- Weights ----
    W_TO: float = 0.0          # max takeoff weight [lb]
    W_empty: float = 0.0       # operating empty weight [lb]
    W_fuel_max: float = 0.0    # max fuel capacity [lb]
    W_payload: float = 0.0     # design payload [lb]

    # ---- Wing geometry ----
    S: float = 0.0             # wing reference area [ft^2]
    b: float = 0.0             # wing span [ft]
    AR: float = 0.0            # aspect ratio  (b^2 / S)
    sweep_qc_deg: float = 0.0  # quarter-chord sweep [deg]
    t_c: float = 0.0           # thickness-to-chord ratio
    taper: float = 0.0         # taper ratio (lambda)

    # ---- Aerodynamics (parabolic polar: CD = CD0 + K*CL^2) ----
    CD0: float = 0.0           # zero-lift drag coefficient
    K: float = 0.0             # induced-drag factor  (1 / pi*e*AR)
    e: float = 0.0             # Oswald span efficiency
    CL_max_clean: float = 0.0  # max CL, clean config
    CL_max_TO: float = 0.0     # max CL, takeoff flaps
    CL_max_L: float = 0.0      # max CL, landing flaps

    # ---- Propulsion (jet) ----
    n_engines: int = 2
    T_max_SL: float = 0.0      # total installed static thrust at SL [lb]
    TSFC: float = 0.0          # thrust-specific fuel consumption [1/hr]
    BPR: float = 0.0           # bypass ratio

    # ---- Propulsion (prop) -- set if propeller aircraft ----
    is_prop: bool = False
    P_bhp: float = 0.0         # total shaft power [bhp]
    eta_p: float = 0.0         # propeller efficiency
    C_bhp: float = 0.0         # brake-specific fuel consumption [lb/hr/bhp]

    # ---- Landing gear ----
    mu_roll: float = 0.03      # rolling friction coefficient (hard runway)
    mu_brake: float = 0.50     # braking friction coefficient (Raymer Table 17.1, civil)

    # ---- Takeoff / landing ----
    h_obstacle_TO: float = 35.0  # obstacle height, takeoff [ft]
    h_obstacle_L: float = 50.0   # obstacle height, landing [ft]
    CL_ground: float = 0.10      # ground-roll lift coefficient

    # ---- Design cruise ----
    M_cruise: float = 0.0      # design cruise Mach number
    h_cruise_ft: float = 0.0   # design cruise altitude [ft]

    # ---- Structural / speed limits ----
    n_max: float = 2.5         # max positive load factor (transport)
    V_NE_kts: float = 0.0      # never-exceed speed [KEAS]
    M_mo: float = 0.0          # max operating Mach number
    h_max_ft: float = 0.0      # max pressurization altitude [ft] (0 = no limit)

    # ---- Derived helpers ----
    @property
    def W_landing(self):
        """Typical landing weight ~ 85% MTOW (Raymer Ch.17)."""
        return 0.85 * self.W_TO

    @property
    def wing_loading(self):
        """W/S at MTOW [lb/ft^2]."""
        return self.W_TO / self.S

    @property
    def thrust_to_weight(self):
        """T/W at MTOW, sea-level static."""
        return self.T_max_SL / self.W_TO

    def __post_init__(self):
        """Compute AR if not provided."""
        if self.AR == 0.0 and self.b > 0 and self.S > 0:
            self.AR = self.b**2 / self.S
