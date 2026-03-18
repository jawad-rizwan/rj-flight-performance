#!/usr/bin/env python3
"""
sync_aero.py — Pull latest aerodynamic coefficients from rj-basic-aerodynamics.

Imports geometry from the sibling aero repo, runs the Raymer Ch.12
component buildup, and updates data/zrj*.py with computed values:
  CD0, K, e, CL_max_clean, CL_max_TO, CL_max_L
"""

import sys
import os
import re
import numpy as np

# ── Locate sibling aero repo ───────────────────────────────────
AERO_REPO = os.path.abspath(
    os.path.join(os.path.dirname(__file__), "..", "rj-basic-aerodynamics")
)
if not os.path.isdir(AERO_REPO):
    sys.exit(f"Error: aero repo not found at {AERO_REPO}")
sys.path.insert(0, AERO_REPO)

from aero.form_factors import ff_wing, ff_fuselage, ff_nacelle, ff_tail_with_hinge
from aero.parasite_drag import cd0_component_buildup
from aero.lift import (
    cl_alpha_subsonic,
    fuselage_lift_factor,
    cl_max_clean,
    cl_max_flaps,
    effective_aspect_ratio_winglet,
)
from aero.induced_drag import oswald_e, k_factor_leading_edge_suction

from data.ZRJ70 import AIRCRAFT as ZRJ70_GEOM
from data.ZRJ100 import AIRCRAFT as ZRJ100_GEOM


# ── Formatting per field ────────────────────────────────────────
FIELD_FMT = {
    "CD0":          ".5f",
    "K":            ".5f",
    "e":            ".3f",
    "CL_max_clean": ".2f",
    "CL_max_TO":    ".2f",
    "CL_max_L":     ".2f",
}


# ── Core computation ───────────────────────────────────────────
def compute_aero(ac):
    """Compute mission-relevant aero coefficients from geometry data."""
    S_ref       = ac["S_ref"]
    AR          = ac["AR"]
    b           = ac["b"]
    MAC         = ac["MAC"]
    cruise_mach = ac["cruise_mach"]
    cruise_alt  = ac["cruise_alt"]

    # Effective aspect ratio
    AR_eff = effective_aspect_ratio_winglet(AR, ac["winglet_h"], b)

    # Lift-curve slope
    F  = fuselage_lift_factor(ac["fuse_d"], b)
    sf = min((ac["S_exposed"] / S_ref) * F, 0.98)
    CLa = cl_alpha_subsonic(
        A=AR_eff, mach=cruise_mach,
        sweep_max_t=np.radians(ac["sweep_mt_deg"]),
        eta=0.95, s_exposed_ratio=sf, F=1.0,
    )

    # Wetted areas
    S_wet_wing  = 2.0 * ac["S_exposed"]        * (1 + 0.25 * ac["t_c_wing"])
    S_wet_htail = 2.0 * ac["S_htail_exposed"]  * (1 + 0.25 * ac["t_c_htail"])
    S_wet_vtail = 2.0 * ac["S_vtail_exposed"]  * (1 + 0.25 * ac["t_c_vtail"])
    S_wet_nac   = np.pi * ac["nacelle_d"] * ac["nacelle_length"] * 0.8
    f = ac["fuse_length"] / ac["fuse_d"]
    S_wet_fuse  = (np.pi * ac["fuse_d"] * ac["fuse_length"]
                   * (1 - 2 / f) ** (2 / 3) * (1 + 1 / f ** 2))

    # Form factors
    sweep_mt    = np.radians(ac["sweep_mt_deg"])
    sweep_mt_ht = np.radians(ac["sweep_mt_htail_deg"])
    sweep_mt_vt = np.radians(ac["sweep_mt_vtail_deg"])
    FF_wing  = ff_wing(ac["x_c_max_wing"], ac["t_c_wing"], cruise_mach, sweep_mt)
    FF_fuse  = ff_fuselage(ac["fuse_length"], ac["fuse_d"])
    FF_htail = ff_tail_with_hinge(ac["x_c_max_htail"], ac["t_c_htail"], cruise_mach, sweep_mt_ht)
    FF_vtail = ff_tail_with_hinge(ac["x_c_max_vtail"], ac["t_c_vtail"], cruise_mach, sweep_mt_vt)
    FF_nac   = ff_nacelle(ac["nacelle_length"], ac["nacelle_d"])

    # CD0 component buildup (Raymer Eq. 12.24)
    components = [
        {"name": "Wing",     "s_wet": S_wet_wing,  "length": MAC,
         "ff": FF_wing,  "Q": ac["Q_wing"],  "pct_laminar": ac["lam_wing"],  "k": ac["k_composite"]},
        {"name": "Fuselage", "s_wet": S_wet_fuse,  "length": ac["fuse_length"],
         "ff": FF_fuse,  "Q": ac["Q_fuse"],  "pct_laminar": ac["lam_fuse"],  "k": ac["k_metal"]},
        {"name": "H-tail",   "s_wet": S_wet_htail, "length": ac["MAC_htail"],
         "ff": FF_htail, "Q": ac["Q_htail"], "pct_laminar": ac["lam_tail"],  "k": ac["k_composite"]},
        {"name": "V-tail",   "s_wet": S_wet_vtail, "length": ac["MAC_vtail"],
         "ff": FF_vtail, "Q": ac["Q_vtail"], "pct_laminar": ac["lam_tail"],  "k": ac["k_composite"]},
    ]
    for i in range(ac["n_nacelles"]):
        components.append({
            "name": f"Nacelle {'L' if i == 0 else 'R'}",
            "s_wet": S_wet_nac, "length": ac["nacelle_length"],
            "ff": FF_nac, "Q": ac["Q_nac"], "pct_laminar": ac["lam_nac"], "k": ac["k_composite"],
        })

    result = cd0_component_buildup(
        components, S_ref, cruise_mach, cruise_alt,
        cd_misc=0.0, leak_pct=ac["leak_pct"],
    )
    cd0 = result["cd0_total"]

    # Oswald e — use LE suction method (recommended by aero repo)
    K_les = k_factor_leading_edge_suction(AR_eff, CLa, ac["S_suction"])
    e_les = 1.0 / (np.pi * AR_eff * K_les)

    # CLmax
    sweep_qc_rad = np.radians(ac["sweep_qc_deg"])
    cl_clean = cl_max_clean(ac["cl_max_airfoil"], sweep_qc_rad)

    sweep_hl_rad = np.radians(ac["sweep_hl_deg"])
    cl_to = cl_max_flaps(
        cl_clean,
        ac["delta_cl_te_factor"] * ac["te_flap_takeoff_fraction"],
        ac["s_flapped_ratio"], sweep_hl_rad,
        delta_cl_max_le=0.0,
    )
    cl_land = cl_max_flaps(
        cl_clean,
        ac["delta_cl_te_factor"], ac["s_flapped_ratio"], sweep_hl_rad,
        delta_cl_max_le=(0.9 * ac["delta_cl_le_factor"]
                         * ac["s_slat_ratio"] * np.cos(sweep_hl_rad)),
    )

    return {
        "CD0":          float(cd0),
        "K":            float(K_les),
        "e":            float(e_les),
        "CL_max_clean": float(cl_clean),
        "CL_max_TO":    float(cl_to),
        "CL_max_L":     float(cl_land),
    }


# ── File updater (dataclass keyword-argument format) ───────────
def update_aircraft_file(path, values):
    """Update numeric values in a dataclass-based aircraft .py file."""
    with open(path) as f:
        text = f.read()

    for key, val in values.items():
        fmt = FIELD_FMT.get(key, ".5f")
        formatted = format(val, fmt)
        # Match:  field=<number>,  (preserving trailing comment)
        pattern = rf'(\b{key}=)-?[\d.]+(\s*,)'
        replacement = rf"\g<1>{formatted}\2"
        text, n = re.subn(pattern, replacement, text)
        if n == 0:
            print(f"  WARNING: '{key}' not found in {os.path.basename(path)}")
        else:
            print(f"  {key:<14s} = {formatted}")

    with open(path, "w") as f:
        f.write(text)


# ── Main ───────────────────────────────────────────────────────
AIRCRAFT_MAP = {
    "zrj70":  ZRJ70_GEOM,
    "zrj100": ZRJ100_GEOM,
}

AIRCRAFT_DIR = os.path.join(os.path.dirname(__file__), "data")


if __name__ == "__main__":
    print("Syncing aerodynamic data from rj-basic-aerodynamics...\n")
    for name, geom in AIRCRAFT_MAP.items():
        print(f"{name}:")
        vals = compute_aero(geom)
        update_aircraft_file(os.path.join(AIRCRAFT_DIR, f"{name}.py"), vals)
        print()
    print("Done.")
