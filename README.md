# rj-flight-performance

Aircraft performance analysis toolkit based on **Raymer, Chapter 17: Performance and Flight Mechanics** (*Aircraft Design: A Conceptual Approach*, 7th Edition).

## Structure

```
rj-flight-performance/
├── perf/                      # Performance calculation library
│   ├── atmosphere.py          # ISA standard atmosphere
│   ├── level_flight.py        # 17.2  Steady level flight
│   ├── range_endurance.py     # 17.2.4-17.2.9  Range & endurance
│   ├── climb.py               # 17.3  Climbing & descending
│   ├── turning.py             # 17.4  Level turning flight
│   ├── glide.py               # 17.5  Gliding flight
│   ├── energy.py              # 17.6  Energy-maneuverability
│   ├── takeoff.py             # 17.8  Takeoff analysis
│   └── landing.py             # 17.9  Landing analysis
├── aircraft/                  # Aircraft data files (easy to swap)
│   ├── base.py                # AircraftData dataclass
│   ├── crj700.py              # Bombardier CRJ-700
│   └── crj1000.py             # Bombardier CRJ-1000
├── examples/
│   ├── full_analysis.py       # Complete numerical analysis
│   └── plot_performance.py    # Generate all performance charts
└── requirements.txt           # numpy, matplotlib
```

## Quick Start

```bash
pip install -r requirements.txt
python examples/full_analysis.py       # Print performance tables
python examples/plot_performance.py    # Generate charts in examples/charts/
```

## Changing Aircraft

Edit `aircraft/crj700.py` or `aircraft/crj1000.py`, or create a new file following the same pattern. All parameters are defined in a single `AircraftData` dataclass in `aircraft/base.py`.

## Methods Implemented

| Section | Topic | Key Equations |
|---------|-------|---------------|
| 17.2 | Steady level flight | 17.8-17.21 |
| 17.2.4 | Range (Breguet) | 17.22-17.28 |
| 17.2.7 | Endurance / loiter | 17.29-17.34 |
| 17.3 | Climb | 17.36-17.51 |
| 17.4 | Turning flight | 17.52-17.55 |
| 17.5 | Gliding flight | 17.62-17.81 |
| 17.6 | Energy maneuverability | 17.84-17.97 |
| 17.8 | Takeoff | 17.100-17.115 |
| 17.9 | Landing | 17.102 (reverse) |

## Charts Generated

1. Thrust required & available vs velocity
2. Power required vs velocity
3. L/D vs CL
4. Rate of climb vs altitude
5. Specific excess power (Ps) vs Mach
6. Turn rate & load factor vs velocity
7. Payload-range diagram
8. Takeoff & landing distance breakdown
9. Glide polar (sink rate vs velocity)

## Units

Imperial throughout (ft, lb, slugs, sec), consistent with Raymer.
