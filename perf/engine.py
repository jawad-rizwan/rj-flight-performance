"""
Workbook-backed engine performance interpolation.

Uses a small XLSX workbook with three slices:
  - altitude sweep at low Mach
  - altitude sweep at cruise Mach
  - speed sweep at 35,000 ft

The workbook is parsed with the Python standard library so the repo
does not require openpyxl/pandas just to use a compact engine deck.
"""

from __future__ import annotations

import os
import zipfile
from dataclasses import dataclass
from functools import lru_cache
from pathlib import Path
from xml.etree import ElementTree as ET

import numpy as np

from .atmosphere import fps_to_kts, speed_of_sound


DEFAULT_ENGINE_WORKBOOK = Path(
    os.environ.get(
        "RJ_ENGINE_PERFORMANCE_XLSX",
        "/home/jawad/Documents/Engine performance.xlsx",
    )
)

_MAIN_NS = "http://schemas.openxmlformats.org/spreadsheetml/2006/main"
_REL_NS = "http://schemas.openxmlformats.org/officeDocument/2006/relationships"


@dataclass(frozen=True)
class EngineDeck:
    low_mach: float
    high_mach: float
    altitude_ft: np.ndarray
    thrust_low: np.ndarray
    thrust_high: np.ndarray
    tsfc_low: np.ndarray
    tsfc_high: np.ndarray
    speed_sweep_kts: np.ndarray
    sweep_mach: np.ndarray
    sweep_thrust: np.ndarray
    sweep_tsfc: np.ndarray


def has_engine_deck(path: str | os.PathLike[str] | None = None) -> bool:
    """Whether the default or provided engine workbook exists."""
    return Path(path or DEFAULT_ENGINE_WORKBOOK).exists()


def _interp_clamped(x, xp, fp):
    return np.interp(np.asarray(x, dtype=float), xp, fp, left=fp[0], right=fp[-1])


def _interp_linear_extrap(x, xp, fp):
    """Linear interpolation with endpoint extrapolation."""
    x_arr = np.asarray(x, dtype=float)
    x_flat = np.atleast_1d(x_arr)
    y = np.interp(x_flat, xp, fp)

    left = x_flat < xp[0]
    if np.any(left):
        slope = (fp[1] - fp[0]) / (xp[1] - xp[0])
        y[left] = fp[0] + slope * (x_flat[left] - xp[0])

    right = x_flat > xp[-1]
    if np.any(right):
        slope = (fp[-1] - fp[-2]) / (xp[-1] - xp[-2])
        y[right] = fp[-1] + slope * (x_flat[right] - xp[-1])

    return float(y[0]) if np.ndim(x) == 0 else y


def _xlsx_shared_strings(archive: zipfile.ZipFile) -> list[str]:
    if "xl/sharedStrings.xml" not in archive.namelist():
        return []

    root = ET.fromstring(archive.read("xl/sharedStrings.xml"))
    strings = []
    for item in root:
        texts = [t.text or "" for t in item.iter(f"{{{_MAIN_NS}}}t")]
        strings.append("".join(texts))
    return strings


def _sheet_rows(path: Path) -> dict[str, list[list[str]]]:
    with zipfile.ZipFile(path) as archive:
        workbook = ET.fromstring(archive.read("xl/workbook.xml"))
        rels = ET.fromstring(archive.read("xl/_rels/workbook.xml.rels"))
        rel_map = {rel.attrib["Id"]: rel.attrib["Target"] for rel in rels}
        shared_strings = _xlsx_shared_strings(archive)

        sheets = {}
        for sheet in workbook.find(f"{{{_MAIN_NS}}}sheets"):
            name = sheet.attrib["name"]
            rel_id = sheet.attrib[f"{{{_REL_NS}}}id"]
            target = rel_map[rel_id]
            if not target.startswith("xl/"):
                target = f"xl/{target}"

            root = ET.fromstring(archive.read(target))
            rows = []
            for row in root.findall(f".//{{{_MAIN_NS}}}row"):
                values = []
                for cell in row.findall(f"{{{_MAIN_NS}}}c"):
                    raw = cell.find(f"{{{_MAIN_NS}}}v")
                    if raw is None:
                        values.append("")
                        continue
                    if cell.attrib.get("t") == "s":
                        values.append(shared_strings[int(raw.text)])
                    else:
                        values.append(raw.text or "")
                rows.append(values)
            sheets[name] = rows
        return sheets


def _numeric_table(rows: list[list[str]], cols: int = 3) -> np.ndarray:
    data = []
    for row in rows[1:]:
        if len(row) < cols:
            continue
        try:
            data.append([float(row[i]) for i in range(cols)])
        except ValueError:
            continue
    return np.asarray(data, dtype=float)


@lru_cache(maxsize=4)
def load_engine_deck(path: str | os.PathLike[str] | None = None) -> EngineDeck:
    """Load and cache the workbook-backed engine deck."""
    workbook_path = Path(path or DEFAULT_ENGINE_WORKBOOK)
    sheets = _sheet_rows(workbook_path)

    cruise_rows = _numeric_table(sheets["Cruise 35000ft"])
    high_rows = _numeric_table(sheets["0.78 Mach"])
    low_rows = _numeric_table(sheets["0.144 Mach"])

    a_35k_kts = fps_to_kts(speed_of_sound(35_000.0))
    sweep_mach = cruise_rows[:, 0] / a_35k_kts

    # Add exact low/high Mach anchors so the weighting uses all three sheets
    # and matches the altitude slices exactly at the reference Mach values.
    low_anchor = np.array(
        [[0.144, np.interp(35_000.0, low_rows[:, 0], low_rows[:, 1]),
          np.interp(35_000.0, low_rows[:, 0], low_rows[:, 2])]],
        dtype=float,
    )
    high_anchor = np.array(
        [[0.78, np.interp(35_000.0, high_rows[:, 0], high_rows[:, 1]),
          np.interp(35_000.0, high_rows[:, 0], high_rows[:, 2])]],
        dtype=float,
    )
    sweep_aug = np.vstack((low_anchor, np.column_stack((sweep_mach, cruise_rows[:, 1:3])), high_anchor))
    sweep_aug = sweep_aug[np.argsort(sweep_aug[:, 0])]

    return EngineDeck(
        low_mach=0.144,
        high_mach=0.78,
        altitude_ft=low_rows[:, 0],
        thrust_low=low_rows[:, 1],
        thrust_high=high_rows[:, 1],
        tsfc_low=low_rows[:, 2],
        tsfc_high=high_rows[:, 2],
        speed_sweep_kts=cruise_rows[:, 0],
        sweep_mach=sweep_aug[:, 0],
        sweep_thrust=sweep_aug[:, 1],
        sweep_tsfc=sweep_aug[:, 2],
    )


def _mach_weight(deck: EngineDeck, mach):
    mach = np.asarray(mach, dtype=float)
    thrust_at_mach = _interp_clamped(mach, deck.sweep_mach, deck.sweep_thrust)
    thrust_low_35 = np.interp(35_000.0, deck.altitude_ft, deck.thrust_low)
    thrust_high_35 = np.interp(35_000.0, deck.altitude_ft, deck.thrust_high)
    denom = thrust_high_35 - thrust_low_35
    if abs(denom) < 1e-9:
        return np.zeros_like(mach)
    return np.clip((thrust_at_mach - thrust_low_35) / denom, 0.0, 1.0)


def net_thrust_per_engine(alt_ft, mach, path: str | os.PathLike[str] | None = None):
    """Interpolated net thrust per engine [lb]."""
    deck = load_engine_deck(path)
    alt = np.asarray(alt_ft, dtype=float)
    m = np.asarray(mach, dtype=float)
    thrust_low = _interp_linear_extrap(alt, deck.altitude_ft, deck.thrust_low)
    thrust_high = _interp_linear_extrap(alt, deck.altitude_ft, deck.thrust_high)
    weight = _mach_weight(deck, m)
    thrust = np.maximum(thrust_low + weight * (thrust_high - thrust_low), 0.0)
    return float(thrust) if np.ndim(thrust) == 0 else thrust


def tsfc_from_deck(alt_ft, mach, path: str | os.PathLike[str] | None = None):
    """Interpolated TSFC [lb/lb/hr]."""
    deck = load_engine_deck(path)
    alt = np.asarray(alt_ft, dtype=float)
    m = np.asarray(mach, dtype=float)
    tsfc_low = _interp_linear_extrap(alt, deck.altitude_ft, deck.tsfc_low)
    tsfc_high = _interp_linear_extrap(alt, deck.altitude_ft, deck.tsfc_high)

    tsfc_at_mach = _interp_clamped(m, deck.sweep_mach, deck.sweep_tsfc)
    tsfc_low_35 = np.interp(35_000.0, deck.altitude_ft, deck.tsfc_low)
    tsfc_high_35 = np.interp(35_000.0, deck.altitude_ft, deck.tsfc_high)
    denom = tsfc_high_35 - tsfc_low_35
    if abs(denom) < 1e-9:
        weight = np.zeros_like(m)
    else:
        weight = np.clip((tsfc_at_mach - tsfc_low_35) / denom, 0.0, 1.0)

    tsfc = tsfc_low + weight * (tsfc_high - tsfc_low)
    return float(tsfc) if np.ndim(tsfc) == 0 else tsfc


def total_net_thrust(alt_ft, mach, n_engines=2, path: str | os.PathLike[str] | None = None):
    """Interpolated total net thrust [lb]."""
    return n_engines * net_thrust_per_engine(alt_ft, mach, path=path)
