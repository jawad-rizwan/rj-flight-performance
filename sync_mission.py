#!/usr/bin/env python3
"""
sync_mission.py — Pull latest weight/propulsion data from rj-mission-analysis.

Imports weight breakdown and engine data from the sibling mission-analysis repo
and updates data/zrj*.py with computed values:
  W_TO, W_empty, W_fuel_max, W_payload, n_engines, T_max_SL, TSFC, BPR,
  M_cruise, h_cruise_ft
"""

import sys
import os
import re

# ── Locate sibling mission-analysis repo ─────────────────────
MISSION_REPO = os.path.abspath(
    os.path.join(os.path.dirname(__file__), "..", "rj-mission-analysis")
)
if not os.path.isdir(MISSION_REPO):
    sys.exit(f"Error: mission-analysis repo not found at {MISSION_REPO}")
sys.path.insert(0, MISSION_REPO)

from data.ZRJ50 import AIRCRAFT as ZRJ50_DATA
from data.ZRJ70 import AIRCRAFT as ZRJ70_DATA
from data.ZRJ100 import AIRCRAFT as ZRJ100_DATA


# ── Formatting per field ──────────────────────────────────────
FIELD_FMT = {
    "W_TO":        ",.1f",
    "W_empty":     ",.1f",
    "W_fuel_max":  ",.1f",
    "W_payload":   ",.1f",
    "n_engines":   "d",
    "T_max_SL":    ",.1f",
    "TSFC":        ".2f",
    "BPR":         ".1f",
    "M_cruise":    ".2f",
    "h_cruise_ft": ",.1f",
}


def format_value(key, val):
    """Format a numeric value for writing into a Python source file.

    Uses underscore separators for thousands (e.g. 85_888.0) to match
    the existing style in data/zrj*.py.
    """
    fmt = FIELD_FMT.get(key, ".1f")
    formatted = format(val, fmt)
    # Convert comma grouping to underscores: 85,888.0 → 85_888.0
    formatted = formatted.replace(",", "_")
    return formatted


# ── Derive flight-performance fields from mission-analysis dict ──
def derive_fields(ac):
    """Convert mission-analysis dict to flight-performance AircraftData fields."""
    n_crew = ac["n_pilots"] + ac["n_flight_attendants"]
    W_empty = ac["empty_weight"] + n_crew * ac["person_weight"]  # OEW

    W_payload = ac["payload_weight"]
    T_max_SL = ac["max_thrust_per_engine"] * ac["num_engines"]

    # MTOW: use mtow_limit if set, else OEW + payload + fuel
    W_TO = ac.get("mtow_limit", None)
    if W_TO is None:
        # No cap — estimate MTOW (would need iterative sizing; use existing value)
        W_TO = W_empty + W_payload + 21_258  # placeholder, should match repo
    W_TO = float(W_TO)

    W_fuel_max = W_TO - W_empty - W_payload
    if W_fuel_max < 0:
        print(f"  WARNING: negative W_fuel_max ({W_fuel_max:.0f} lb) — check MTOW/weights")
        W_fuel_max = 0.0

    return {
        "W_TO":        W_TO,
        "W_empty":     float(W_empty),
        "W_fuel_max":  W_fuel_max,
        "W_payload":   float(W_payload),
        "n_engines":   ac["num_engines"],
        "T_max_SL":    float(T_max_SL),
        "TSFC":        ac["tsfc_cruise"],
        "BPR":         ac["bypass_ratio"],
        "M_cruise":    ac["cruise_mach"],
        "h_cruise_ft": ac["cruise_altitude_ft"],
    }


# ── File updater (dataclass keyword-argument format) ──────────
def update_aircraft_file(path, values):
    """Update numeric values in a dataclass-based aircraft .py file."""
    with open(path) as f:
        text = f.read()

    for key, val in values.items():
        formatted = format_value(key, val)
        # Match:  field=<number>,  (number may contain underscores)
        # Handles: W_TO=85_888.0, or T_max_SL=38_380.0, or TSFC=0.46,
        pattern = rf'(\b{key}\s*=\s*)-?[\d_]+(?:\.[\d_]+)?(\s*,)'
        replacement = rf"\g<1>{formatted}\2"
        text, n = re.subn(pattern, replacement, text)
        if n == 0:
            print(f"  WARNING: '{key}' not found in {os.path.basename(path)}")
        else:
            print(f"  {key:<14s} = {formatted}")

    with open(path, "w") as f:
        f.write(text)


# ── Main ──────────────────────────────────────────────────────
AIRCRAFT_MAP = {
    "zrj50":  ZRJ50_DATA,
    "zrj70":  ZRJ70_DATA,
    "zrj100": ZRJ100_DATA,
}

AIRCRAFT_DIR = os.path.join(os.path.dirname(__file__), "data")


if __name__ == "__main__":
    print("Syncing mission data from rj-mission-analysis...\n")
    for name, ac_data in AIRCRAFT_MAP.items():
        print(f"{name}:")
        vals = derive_fields(ac_data)
        update_aircraft_file(os.path.join(AIRCRAFT_DIR, f"{name}.py"), vals)
        print()
    print("Done.")
