# rj-flight-performance

Aircraft performance analysis toolkit based on **Raymer, Chapter 17: Performance and Flight Mechanics** (*Aircraft Design: A Conceptual Approach*, 7th Edition).

Includes the **ZRJ family** (ZRJ50, ZRJ70, ZRJ100 high-wing regional jets) and **Bombardier CRJ-700 / CRJ-1000** as reference aircraft. Any aircraft can be added by creating a single data file.

---

## Repository Structure

```
rj-flight-performance/
├── perf/                              # Core performance library
│   ├── atmosphere.py                  # ISA standard atmosphere model
│   ├── level_flight.py                # 17.2   Steady level flight
│   ├── range_endurance.py             # 17.2.4 Range & endurance (Breguet)
│   ├── climb.py                       # 17.3   Climbing & descending flight
│   ├── turning.py                     # 17.4   Level turning flight
│   ├── glide.py                       # 17.5   Gliding flight
│   ├── energy.py                      # 17.6   Energy-maneuverability methods
│   ├── takeoff.py                     # 17.8   Takeoff analysis
│   └── landing.py                     # 17.9   Landing analysis
├── data/                              # Aircraft data files
│   ├── base.py                        # AircraftData dataclass definition
│   ├── crj700.py                      # Bombardier CRJ-700
│   ├── crj1000.py                     # Bombardier CRJ-1000
│   ├── zrj50.py                       # ZRJ50  (50-seat, scope-clause)
│   ├── zrj70.py                       # ZRJ70  (76-seat, NA variant)
│   └── zrj100.py                      # ZRJ100 (100-seat, EU variant)
├── examples/
│   ├── full_analysis.py               # Complete numerical performance analysis
│   ├── plot_performance.py            # Generate all performance charts
│   ├── charts/
│   │   ├── CRJ/                       # CRJ family charts (12 PNGs)
│   │   │   └── surfaces/              # Per-surface TO/landing + BFL charts
│   │   └── ZRJ/                       # ZRJ family charts (12 PNGs)
│   │       └── surfaces/              # Per-surface TO/landing + BFL charts
│   └── output/
│       ├── CRJ/analysis.txt           # CRJ family analysis output
│       └── ZRJ/analysis.txt           # ZRJ family analysis output
├── sync_aero.py                       # Pull aero coefficients from rj-basic-aerodynamics
├── sync_mission.py                    # Pull weights/propulsion from rj-mission-analysis
├── requirements.txt                   # Python dependencies
└── README.md
```

---

## Quick Start

```bash
# Install dependencies
pip install -r requirements.txt

# Run full numerical analysis (CRJ + ZRJ families, output to console + files)
python examples/full_analysis.py

# Generate all performance charts (saved to examples/charts/CRJ/ and ZRJ/)
python examples/plot_performance.py
```

---

## ZRJ Family

Three variants of a high-wing regional jet with PW1200G geared turbofan engines:

| Variant | Seats | MTOW (lb) | OEW (lb) | Wing Area (ft²) | T/W | Cruise Mach |
|---------|-------|-----------|----------|-----------------|-----|-------------|
| ZRJ50   | 50    | 65,000    | 45,578   | 1,016.6          | 0.590 | 0.78        |
| ZRJ70   | 76    | 81,523    | 45,775   | 1,016.6          | 0.471 | 0.78        |
| ZRJ100  | 96    | 89,414    | 48,341   | 1,016.6          | 0.429 | 0.78        |

All three share the same wing (AR 7.8, 22.9 deg sweep) and engine (2x PW1200G, 38,380 lb total thrust, BPR 9.0). The ZRJ50 uses the same airframe as the ZRJ70 but is MTOW-capped at 65,000 lb per scope clause.

### Sample Output (ZRJ Family)

```
  Metric                               ZRJ50         ZRJ70        ZRJ100
  ----------------------------------------------------------------------
  (L/D)_max                            15.82         15.82         15.54
  Cruise L/D                           12.91         14.45         14.59
  Range (nmi)                          1,376         2,823         2,574
  Best ROC (fpm)                      17,672        13,986        12,462
  Corner spd (kts)                       202           226           237
  TODR (ft)                            2,327         3,176         3,645
  TODR FAR (ft)                        2,677         3,653         4,192
  ASDR at V1 (ft)                      2,439         3,395         3,956
  BFL iterative (ft)                   2,439         3,396         3,956
  BFL Eq.17.113 (ft)                   2,589         3,662         4,266
  TOFL FAR 25 (ft)                     2,677         3,653         4,192
  LDR (ft)                             2,293         2,557         2,682
  LDR FAR (ft)                         3,822         4,262         4,469
```

### MR&O Performance Attributes

| Category | MR&O Attribute | ZRJ50 | ZRJ70 | ZRJ100 |
|----------|---------------|-------|-------|--------|
| **Airworthiness** | Transport Canada, EASA, FAA (Part 25) | Yes | Yes | Yes |
| **Performance** | Range, BOW, LRC/FAR121, 100 nm alt. | 4,056 nm | 7,875 nm | 7,929 nm |
| | Max cruise range | 4,479 nm | 4,860 nm | 4,604 nm |
| | Range, Full PAX, LRC | 1,531 nm | 3,722 nm | 3,009 nm |
| | Operational range, Full PAX (mission analysis) | 489 nm | 1,800 nm | 1,549 nm |
| | Normal cruise (Optimized Mach) | M 0.78 | M 0.78 | M 0.78 |
| | Maximum cruise (Optimized Mach) | M 0.85 | M 0.85 | M 0.85 |
| | Takeoff BFL, ft | 2,439 | 3,396 | 3,956 |
| | Initial cruise altitude, ft | 35,000 | 35,000 | 35,000 |
| | Maximum cruise altitude, ft | >65,000 | 63,000 | 59,500 |
| | Single engine climb | 19.6% | 13.6% | 11.5% |
| **Weight & Payload** | Max. payload, lbs | 11,350 | 18,055 | 23,380 |
| | Passengers (<100 PAX market segment) | 50 | 76 | 96 |
| | Max. landing weight (~85% MTOW) | 55,250 | 69,295 | 76,002 |
| **Airfield Performance** | TOFL, ISA, SL (MTOW) | 2,677 ft | 3,653 ft | 4,192 ft |
| | LD, ISA, SL (MLW) | 3,822 ft | 4,262 ft | 4,469 ft |
| | Steep approach capability (MLW) | Yes | Yes | Yes |
| | Rate of climb [SL, AEO / OEI] | 17,672 / 6,006 fpm | 13,986 / 4,643 fpm | 12,462 / 4,070 fpm |
| | Glideslope | 3.0 deg | 3.0 deg | 3.0 deg |

---

## Sync Scripts

Two scripts pull data from sibling repos to keep aircraft definitions up to date:

| Script | Source Repo | Fields Synced |
|--------|------------|---------------|
| `sync_aero.py` | `rj-basic-aerodynamics` | CD0, K, e, CL_max_clean, CL_max_TO, CL_max_L |
| `sync_mission.py` | `rj-mission-analysis` | W_TO, W_empty, W_fuel_max, W_payload, T_max_SL, TSFC, BPR, n_engines, M_cruise, h_cruise_ft |

```bash
python sync_aero.py      # update aero coefficients from component buildup
python sync_mission.py   # update weights/propulsion from mission analysis
```

Both scripts use regex-based file updaters that modify values in-place in `data/zrj*.py` while preserving comments and formatting.

---

## Changing Aircraft Data

All aircraft parameters live in a single `AircraftData` dataclass. To analyse a different aircraft:

1. **Copy** `data/crj700.py` to a new file (e.g. `data/my_aircraft.py`)
2. **Edit** the parameters (see Data Checklist below)
3. **Import** it in `data/__init__.py`:
```python
from .my_aircraft import my_aircraft
```
4. **Add** it to a family in `examples/full_analysis.py` and `examples/plot_performance.py`.

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
| 17.2.3 | Maximum speed | 17.8-17.11 | `level_flight.py` |
| 17.2.4 | Range (Breguet, jet) | 17.22-17.27 | `range_endurance.py` |
| 17.2.5 | Range optimization (jet) | 17.24-17.27 | `range_endurance.py` |
| 17.2.6 | Range optimization (prop) | 17.28 | `range_endurance.py` |
| 17.2.7 | Loiter endurance | 17.29-17.30 | `range_endurance.py` |
| 17.2.8 | Loiter optimization (jet) | 17.30 | `range_endurance.py` |
| 17.2.9 | Loiter optimization (prop) | 17.31-17.33 | `range_endurance.py` |
| 17.2.10 | Loiter-cruise relationship | 17.34 | `range_endurance.py` |
| 17.3.1 | Climb equations of motion | 17.36-17.41 | `climb.py` |
| 17.3.2 | Climb hodograph | 17.38-17.39 | `climb.py` |
| 17.3.3 | Best angle/rate of climb (jet) | 17.42-17.43 | `climb.py` |
| 17.3.4 | Best angle/rate of climb (prop) | 17.44-17.45 | `climb.py` |
| 17.3.5 | Time to climb & fuel to climb | 17.46-17.51 | `climb.py` |
| — | FAR 25 climb gradients (Table F.4) | FAR 25.111-25.121 | `climb.py` |
| 17.4 | Level turning flight | 17.52 | `turning.py` |
| 17.4.1 | Instantaneous turn rate | 17.52 | `turning.py` |
| 17.4.2 | Sustained turn rate | 17.53-17.55 | `turning.py` |
| 17.5.1 | Straight gliding flight | 17.62-17.74 | `glide.py` |
| 17.5.2 | Turning glide | 17.75-17.81 | `glide.py` |
| 17.6.1 | Energy equations | 17.84-17.89 | `energy.py` |
| 17.6.2 | Ps contours (energy maneuverability) | 17.89, Fig. 17.9 | `energy.py` |
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

Running `python examples/plot_performance.py` generates 12 charts per family in `examples/charts/{CRJ,ZRJ}/`, plus per-surface subcharts in `surfaces/`:

| # | Chart | Raymer Reference |
|---|-------|-----------------|
| 1 | Thrust required & available vs velocity | Fig. 17.2, Eq. 17.8 |
| 2 | Power required vs velocity | Eq. 17.16-17.17 |
| 3 | L/D vs CL | Eq. 17.14, 17.67 |
| 4 | Rate of climb vs altitude | Eq. 17.43, 17.39 |
| 5 | Specific excess power (Ps) vs Mach | Fig. 17.9, Eq. 17.89 |
| 6 | Turn rate & load factor vs velocity | Fig. 17.6, Eq. 17.52, 17.54 |
| 7 | Takeoff & landing distance breakdown | Sec. 17.8, 17.9 |
| 8 | Glide polar (sink rate vs velocity) | Fig. 17.7, Eq. 17.68 |
| 9 | Operating envelope (altitude vs Mach) | Eq. 17.8, 17.12 |
| 10 | Climb hodograph (Vx vs Vy) | Sec. 17.3.2, Eq. 17.38-17.39 |
| 11 | Ps contour plot (altitude vs Mach) | Fig. 17.9, Eq. 17.89 |
| 12 | Airfield performance by runway surface | Table 17.1, Sec. 17.8-17.9 |
| — | Per-surface TO/landing + BFL (in `surfaces/`) | Table 17.1, Eq. 17.102 |

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
| Braking friction coeff. | `mu_brake` | 0.50 | Hard-surface runway braking friction (Raymer Table 17.1) |

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
