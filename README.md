# rj-flight-performance

Aircraft performance analysis toolkit based on **Raymer, Chapter 17: Performance and Flight Mechanics** (*Aircraft Design: A Conceptual Approach*, 7th Edition).

Currently configured for the **Bombardier CRJ-700** and **CRJ-1000**, but designed so any aircraft can be swapped in by editing a single data file.

---

## Repository Structure

```
rj-flight-performance/
├── perf/                          # Core performance library
│   ├── atmosphere.py              # ISA standard atmosphere model
│   ├── level_flight.py            # 17.2   Steady level flight
│   ├── range_endurance.py         # 17.2.4 Range & endurance (Breguet)
│   ├── climb.py                   # 17.3   Climbing & descending flight
│   ├── turning.py                 # 17.4   Level turning flight
│   ├── glide.py                   # 17.5   Gliding flight
│   ├── energy.py                  # 17.6   Energy-maneuverability methods
│   ├── takeoff.py                 # 17.8   Takeoff analysis
│   └── landing.py                 # 17.9   Landing analysis
├── aircraft/                      # Aircraft data files
│   ├── base.py                    # AircraftData dataclass definition
│   ├── crj700.py                  # Bombardier CRJ-700 parameters
│   └── crj1000.py                 # Bombardier CRJ-1000 parameters
├── examples/
│   ├── full_analysis.py           # Complete numerical performance analysis
│   ├── plot_performance.py        # Generate all performance charts
│   └── charts/                    # Generated PNG figures (9 charts)
├── requirements.txt               # Python dependencies
└── README.md
```

---

## Quick Start

```bash
# Install dependencies
pip install -r requirements.txt

# Run full numerical analysis (prints tables to console)
python examples/full_analysis.py

# Generate all performance charts (saved to examples/charts/)
python examples/plot_performance.py
```

---

## Changing Aircraft Data

All aircraft parameters live in a single `AircraftData` dataclass. To analyse a different aircraft:

1. **Copy** `aircraft/crj700.py` to a new file (e.g. `aircraft/my_aircraft.py`)
2. **Edit** the parameters:

```python
from .base import AircraftData

my_aircraft = AircraftData(
    name="My Aircraft",

    # Weights
    W_TO=75000.0,          # lb, max takeoff weight
    W_empty=44000.0,       # lb, operating empty weight
    W_fuel_max=19500.0,    # lb, max usable fuel
    W_payload=17500.0,     # lb, design payload

    # Wing geometry
    S=520.0,               # ft^2, reference area
    b=76.3,                # ft, span
    AR=11.2,               # aspect ratio
    sweep_qc_deg=25.0,     # deg, quarter-chord sweep
    t_c=0.11,              # thickness-to-chord ratio
    taper=0.32,            # taper ratio

    # Aerodynamics (parabolic drag polar: CD = CD0 + K*CL^2)
    CD0=0.022,             # zero-lift drag coefficient
    K=0.040,               # induced drag factor
    e=0.80,                # Oswald span efficiency
    CL_max_clean=1.40,     # max CL, clean
    CL_max_TO=1.80,        # max CL, takeoff flaps
    CL_max_L=2.40,         # max CL, landing flaps

    # Propulsion (jet)
    n_engines=2,
    T_max_SL=26760.0,      # lb, total sea-level static thrust
    TSFC=0.68,             # 1/hr, cruise TSFC
    BPR=5.0,               # bypass ratio

    # Design cruise
    M_cruise=0.78,
    h_cruise_ft=35000.0,   # ft

    # Structural limits
    n_max=2.5,             # max load factor
)
```

3. **Import** it in `aircraft/__init__.py`:
```python
from .my_aircraft import my_aircraft
```

4. **Add** it to the `VARIANTS` list in `examples/full_analysis.py` and `examples/plot_performance.py`.

### Propeller Aircraft

Set `is_prop=True` and fill in `P_bhp`, `eta_p`, and `C_bhp` instead of the jet thrust parameters. The library includes both jet and prop versions of the Breguet range/endurance equations.

---

## Methods Implemented

All equations reference Raymer's *Aircraft Design: A Conceptual Approach*, 7th Edition, Chapter 17.

| Section | Topic | Key Equations | Module |
|---------|-------|---------------|--------|
| 17.1 | Equations of motion | 17.1-17.7 | `level_flight.py` |
| 17.2 | Steady level flight | 17.8-17.15 | `level_flight.py` |
| 17.2.1 | Minimum thrust required | 17.12-17.14 | `level_flight.py` |
| 17.2.2 | Minimum power required | 17.16-17.21 | `level_flight.py` |
| 17.2.4 | Range (Breguet, jet) | 17.22-17.27 | `range_endurance.py` |
| 17.2.5 | Range optimization (jet) | 17.24-17.27 | `range_endurance.py` |
| 17.2.6 | Range optimization (prop) | 17.28 | `range_endurance.py` |
| 17.2.7 | Loiter endurance | 17.29-17.30 | `range_endurance.py` |
| 17.2.8 | Loiter optimization (jet) | 17.30 | `range_endurance.py` |
| 17.2.9 | Loiter optimization (prop) | 17.31-17.33 | `range_endurance.py` |
| 17.2.10 | Loiter-cruise relationship | 17.34 | `range_endurance.py` |
| 17.3.1 | Climb equations of motion | 17.36-17.41 | `climb.py` |
| 17.3.3 | Best angle/rate of climb (jet) | 17.42-17.43 | `climb.py` |
| 17.3.4 | Best angle/rate of climb (prop) | 17.44-17.45 | `climb.py` |
| 17.3.5 | Time to climb & fuel to climb | 17.46-17.51 | `climb.py` |
| 17.4 | Level turning flight | 17.52 | `turning.py` |
| 17.4.1 | Instantaneous turn rate | 17.52 | `turning.py` |
| 17.4.2 | Sustained turn rate | 17.53-17.55 | `turning.py` |
| 17.5.1 | Straight gliding flight | 17.62-17.74 | `glide.py` |
| 17.5.2 | Turning glide | 17.75-17.81 | `glide.py` |
| 17.6.1 | Energy equations | 17.84-17.89 | `energy.py` |
| 17.6.3 | Minimum time-to-climb | 17.91-17.93 | `energy.py` |
| 17.6.4 | Minimum fuel-to-climb | 17.94-17.96 | `energy.py` |
| 17.6.5 | Mission-segment weight fraction | 17.97 | `energy.py` |
| 17.8.1 | Takeoff ground roll | 17.100-17.104 | `takeoff.py` |
| 17.8.2 | Takeoff transition | 17.105-17.111 | `takeoff.py` |
| 17.8.3 | Takeoff climb | 17.112 | `takeoff.py` |
| 17.8.4 | Balanced field length | 17.113-17.115 | `takeoff.py` |
| 17.9.1 | Landing approach | -- | `landing.py` |
| 17.9.2 | Landing flare | 17.107, 17.110 | `landing.py` |
| 17.9.3 | Landing ground roll (braking) | 17.102-17.104 | `landing.py` |

---

## Performance Charts

Running `python examples/plot_performance.py` generates 10 charts in `examples/charts/`:

| # | Chart | Raymer Reference |
|---|-------|-----------------|
| 1 | Thrust required & available vs velocity | Fig. 17.2, Eq. 17.8 |
| 2 | Power required vs velocity | Eq. 17.16-17.17 |
| 3 | L/D vs CL | Eq. 17.14, 17.67 |
| 4 | Rate of climb vs altitude | Eq. 17.43, 17.39 |
| 5 | Specific excess power (Ps) vs Mach | Fig. 17.9, Eq. 17.89 |
| 6 | Turn rate & load factor vs velocity | Fig. 17.6, Eq. 17.52, 17.54 |
| 7 | Payload-range diagram | Eq. 17.23 |
| 8 | Takeoff & landing distance breakdown | Sec. 17.8, 17.9 |
| 9 | ASDR & AGDR vs engine failure speed (BFL) | Sec. 17.8.4, Eq. 17.102 |
| 10 | Glide polar (sink rate vs velocity) | Fig. 17.7, Eq. 17.68 |

---

## Sample Output (CRJ-700 vs CRJ-1000)

```
  Metric                             CRJ-700      CRJ-1000
  --------------------------------------------------------
  (L/D)_max                            16.85         17.01
  Cruise L/D                           16.79         16.97
  Range (nmi)                          2,976         2,359
  Best ROC (fpm)                      11,202         9,393
  Corner spd (kts)                       276           292
  TODR (ft)                            7,455         9,249
  TODR FAR (ft)                        8,573        10,636
  ASDR at V1 (ft)                      8,520        10,587
  BFL iterative (ft)                   8,521        10,587
  BFL Eq.17.113 (ft)                   7,933         9,808
  TOFL FAR 25 (ft)                     8,573        10,636
  LDR (ft)                             4,258         4,592
  LDR FAR (ft)                         7,099         7,656
```

---

## Units

All calculations use **Imperial units** (ft, lb, slugs, seconds), consistent with Raymer. Key conversions:

| Quantity | Unit |
|----------|------|
| Weight | lb (pounds-force) |
| Length / span | ft |
| Area | ft^2 |
| Velocity | ft/s (converted to kts for display) |
| Thrust | lb |
| SFC | 1/hr (TSFC for jets) |
| Power | ft-lb/s (converted to hp for display) |
| Density | slug/ft^3 |
| Dynamic pressure | lb/ft^2 (psf) |

---

## Data Checklist for Your Aircraft

To run a full performance analysis on your own design, you need to gather the values below. All units are **Imperial** (lb, ft, ft², etc.). Fields marked *(optional)* have sensible defaults or are only needed for propeller aircraft.

### Weights

| Parameter | Field | Unit | Description |
|-----------|-------|------|-------------|
| Max takeoff weight | `W_TO` | lb | Maximum takeoff gross weight |
| Operating empty weight | `W_empty` | lb | Empty weight including crew & unusable fuel |
| Max fuel capacity | `W_fuel_max` | lb | Maximum usable fuel weight |
| Design payload | `W_payload` | lb | Design payload (passengers + cargo) |

### Wing Geometry

| Parameter | Field | Unit | Description |
|-----------|-------|------|-------------|
| Reference area | `S` | ft² | Trapezoidal wing planform area |
| Wingspan | `b` | ft | Tip-to-tip span |
| Aspect ratio | `AR` | — | b²/S (auto-computed if `b` and `S` are given) |
| Quarter-chord sweep | `sweep_qc_deg` | deg | Sweep angle at the quarter-chord line |
| Thickness-to-chord ratio | `t_c` | — | Average t/c of the wing |
| Taper ratio | `taper` | — | Tip chord / root chord |

### Aerodynamics

| Parameter | Field | Unit | Description |
|-----------|-------|------|-------------|
| Zero-lift drag coefficient | `CD0` | — | Parasite drag at zero lift |
| Induced-drag factor | `K` | — | K in CD = CD0 + K·CL² (= 1/π·e·AR) |
| Oswald efficiency | `e` | — | Span efficiency factor |
| CL max (clean) | `CL_max_clean` | — | Maximum lift coefficient, clean configuration |
| CL max (takeoff) | `CL_max_TO` | — | Maximum lift coefficient, takeoff flaps |
| CL max (landing) | `CL_max_L` | — | Maximum lift coefficient, landing flaps |

### Propulsion — Jet

| Parameter | Field | Unit | Description |
|-----------|-------|------|-------------|
| Number of engines | `n_engines` | — | Total engine count |
| Total static thrust | `T_max_SL` | lb | Combined sea-level static thrust (all engines) |
| Thrust SFC | `TSFC` | 1/hr | Cruise thrust-specific fuel consumption |
| Bypass ratio | `BPR` | — | Engine bypass ratio |

### Propulsion — Propeller *(set `is_prop=True`)*

| Parameter | Field | Unit | Description |
|-----------|-------|------|-------------|
| Shaft power | `P_bhp` | bhp | Total shaft brake horsepower (all engines) |
| Propeller efficiency | `eta_p` | — | Cruise propeller efficiency (typically 0.7–0.85) |
| Brake SFC | `C_bhp` | lb/hr/bhp | Brake-specific fuel consumption |

### Landing Gear *(optional — defaults shown)*

| Parameter | Field | Default | Description |
|-----------|-------|---------|-------------|
| Rolling friction coeff. | `mu_roll` | 0.03 | Hard-surface runway rolling friction |
| Braking friction coeff. | `mu_brake` | 0.40 | Hard-surface runway braking friction |

### Takeoff / Landing *(optional — defaults shown)*

| Parameter | Field | Default | Description |
|-----------|-------|---------|-------------|
| Takeoff obstacle height | `h_obstacle_TO` | 35 ft | FAR 25 screen height for takeoff |
| Landing obstacle height | `h_obstacle_L` | 50 ft | FAR 25 screen height for landing |
| Ground-roll CL | `CL_ground` | 0.10 | Lift coefficient during ground roll |

### Design Cruise

| Parameter | Field | Unit | Description |
|-----------|-------|------|-------------|
| Cruise Mach number | `M_cruise` | — | Design cruise Mach |
| Cruise altitude | `h_cruise_ft` | ft | Design cruise altitude |

### Structural Limits

| Parameter | Field | Default | Description |
|-----------|-------|---------|-------------|
| Max load factor | `n_max` | 2.5 | Positive limit load factor (FAR 25 transport) |
| Never-exceed speed | `V_NE_kts` | — *(optional)* | VNE in KEAS; only needed for V-n diagrams |

---

## Dependencies

- Python 3.8+
- NumPy >= 1.20
- Matplotlib >= 3.5

---

## References

- Raymer, D. P. *Aircraft Design: A Conceptual Approach*, 7th Edition, AIAA, 2024. Chapter 17: Performance and Flight Mechanics.
