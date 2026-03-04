"""Build unified MIAC master datasets from MIAC and MIAC_m sources.

Outputs:
1) Adds reciprocal conversion columns in source folders:
   - data/MIAC/**/*.csv gains/refreshes miac_m
   - data/MIAC_m/**/*.csv gains/refreshes miac
2) Writes canonical unified master datasets to data/MIAC/** with columns:
   - molality, miac, miac_m (+ composition columns for mixed solvents)
3) Writes comparison tables and plots to data/MIAC_combined/**.
"""

from __future__ import annotations

import csv
import math
from collections import defaultdict
from pathlib import Path
from typing import Dict, Iterable, List, Tuple

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


REPO_ROOT = Path(__file__).resolve().parents[2]
MIAC_ROOT = REPO_ROOT / "data" / "MIAC"
MIAC_M_ROOT = REPO_ROOT / "data" / "MIAC_m"
COMBINED_ROOT = REPO_ROOT / "data" / "MIAC_combined"

MW_SOLVENT = {
    "water": 18.01528e-3,
    "methanol": 32.04e-3,
    "ethanol": 46.068e-3,
}

CATION_CHARGE = {"Li": 1, "Na": 1, "K": 1, "NH4": 1, "H": 1}
ANION_CHARGE = {"Cl": -1, "Br": -1, "I": -1}

SALT_CANONICAL = {
    "LI": "LiI",
    "LII": "LiI",
    "LICL": "LiCl",
    "LIBR": "LiBr",
    "NACL": "NaCl",
    "NABR": "NaBr",
    "NAI": "NaI",
    "KCL": "KCl",
    "KBR": "KBr",
    "KI": "KI",
    "NH4CL": "NH4Cl",
    "NH4BR": "NH4Br",
    "NH4I": "NH4I",
}

COMP_COL_LABELS = {
    "water": "x_H2O",
    "methanol": "x_Methanol",
    "ethanol": "x_Ethanol",
}


def _to_float(value) -> float | None:
    if value is None:
        return None
    text = str(value).strip()
    if not text:
        return None
    try:
        val = float(text)
    except ValueError:
        return None
    if not math.isfinite(val):
        return None
    return val


def _canonical_salt(token: str) -> str:
    key = token.strip().replace("_", "").replace("-", "")
    if not key:
        return token
    upper = key.upper()
    if upper in SALT_CANONICAL:
        return SALT_CANONICAL[upper]

    for c in ("NH4", "Li", "Na", "K", "H"):
        if upper.startswith(c.upper()):
            a = upper[len(c) :]
            if a in {"CL", "BR", "I"}:
                a_fmt = {"CL": "Cl", "BR": "Br", "I": "I"}[a]
                return f"{c}{a_fmt}" if c != "NH4" else f"NH4{a_fmt}"
    return token


def _parse_salt_from_stem(stem: str) -> str:
    token = stem.split("-")[-1] if "-" in stem else stem
    return _canonical_salt(token)


def _sum_nu_from_salt(salt: str) -> int:
    cat = None
    an = None
    for c in sorted(CATION_CHARGE, key=len, reverse=True):
        if salt.startswith(c):
            cat = c
            an = salt[len(c) :]
            break
    if cat is None or an not in ANION_CHARGE:
        return 2
    zc = abs(CATION_CHARGE[cat])
    za = abs(ANION_CHARGE[an])
    g = math.gcd(zc, za)
    return int((za // g) + (zc // g))


def _read_rows(path: Path) -> Tuple[List[str], List[Dict[str, str]]]:
    with path.open("r", newline="", encoding="utf-8-sig") as f:
        reader = csv.DictReader(f)
        raw_fields = [h for h in (reader.fieldnames or []) if h]
        fields: List[str] = []
        seen = set()
        for h in raw_fields:
            hs = h.strip()
            if hs and hs not in seen:
                fields.append(hs)
                seen.add(hs)
        rows = []
        for row in reader:
            clean = {}
            for k, v in row.items():
                if not k:
                    continue
                ks = k.strip()
                if not ks:
                    continue
                clean[ks] = v.strip() if isinstance(v, str) else v
            rows.append(clean)
    return fields, rows


def _write_rows(path: Path, fieldnames: List[str], rows: List[Dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def _field_lookup(fields: Iterable[str]) -> Dict[str, str]:
    return {f.strip().lower(): f for f in fields if f and f.strip()}


def _molality_column(fields: Iterable[str]) -> str | None:
    lookup = _field_lookup(fields)
    for c in ("molality", "molality (kg/mol)", "m (mol/kg)", "m", "x"):
        if c in lookup:
            return lookup[c]
    return None


def _source_value_column(fields: Iterable[str], source: str) -> str | None:
    lookup = _field_lookup(fields)
    if source == "MIAC":
        for c in ("miac", "y", "gamma"):
            if c in lookup:
                return lookup[c]
    else:
        for c in ("miac_m", "gamma"):
            if c in lookup:
                return lookup[c]
    return None


def _insert_after(fieldnames: List[str], anchor: str, new_col: str) -> List[str]:
    out = [f for f in fieldnames if f and f.strip()]
    if new_col in out:
        return out
    if anchor in out:
        i = out.index(anchor)
        return out[: i + 1] + [new_col] + out[i + 1 :]
    return out + [new_col]


def _composition_from_row(row: Dict[str, str], solvent_system: str) -> Tuple[Dict[str, float], str]:
    solvents = [s for s in solvent_system.split("-") if s]
    if len(solvents) <= 1:
        return ({solvents[0]: 1.0} if solvents else {}), "pure"

    lookup = {k.strip().lower(): k for k in row.keys() if k and k.strip()}

    def get_x(solvent: str) -> float | None:
        cands = [f"x_{solvent}"]
        if solvent == "water":
            cands += ["x_h2o", "x_water"]
        if solvent == "methanol":
            cands += ["x_meoh"]
        if solvent == "ethanol":
            cands += ["x_etoh"]
        for c in cands:
            key = lookup.get(c.lower())
            if key is None:
                continue
            v = _to_float(row.get(key))
            if v is not None:
                return v
        return None

    known = {s: get_x(s) for s in solvents}
    present = {k: v for k, v in known.items() if v is not None}

    if len(solvents) == 2 and len(present) == 1:
        s_known, x_known = next(iter(present.items()))
        s_other = [s for s in solvents if s != s_known][0]
        known[s_other] = 1.0 - x_known
        present[s_other] = known[s_other]

    if len(present) == 0:
        comp = {s: 1.0 / len(solvents) for s in solvents}
        mode = "equal-fallback"
    elif len(present) < len(solvents):
        rem = max(0.0, 1.0 - sum(present.values()))
        missing = [s for s in solvents if s not in present]
        each = rem / len(missing) if missing else 0.0
        comp = {s: (present[s] if s in present else each) for s in solvents}
        mode = "partial-composition"
    else:
        denom = sum(present.values())
        if abs(denom) < 1e-12:
            comp = {s: 1.0 / len(solvents) for s in solvents}
            mode = "equal-fallback"
        else:
            comp = {s: present[s] / denom for s in solvents}
            mode = "composition-columns"
    return comp, mode


def _mw_from_composition(comp: Dict[str, float], solvent_system: str) -> float:
    solvents = [s for s in solvent_system.split("-") if s]
    if not solvents:
        raise ValueError("Empty solvent system")
    if len(solvents) == 1:
        return MW_SOLVENT[solvents[0]]
    mw = 0.0
    for s in solvents:
        mw += comp.get(s, 0.0) * MW_SOLVENT[s]
    return mw


def _signature(comp: Dict[str, float], solvent_system: str) -> Tuple[Tuple[str, float], ...]:
    solvents = [s for s in solvent_system.split("-") if s]
    if len(solvents) <= 1:
        return tuple()
    return tuple((s, round(float(comp.get(s, 0.0)), 12)) for s in solvents)


def _dedupe_points(points: List[Tuple[float, float]]) -> Tuple[np.ndarray, np.ndarray]:
    by_m: Dict[float, List[float]] = defaultdict(list)
    for m, y in points:
        by_m[round(float(m), 12)].append(float(y))
    m_sorted = sorted(by_m.keys())
    y_avg = [float(np.mean(by_m[m])) for m in m_sorted]
    return np.asarray(m_sorted, dtype=float), np.asarray(y_avg, dtype=float)


def _interp_to_grid(x_src: np.ndarray, y_src: np.ndarray, grid: np.ndarray) -> np.ndarray:
    if x_src.size == 0:
        return np.full_like(grid, np.nan, dtype=float)
    if x_src.size == 1:
        out = np.full_like(grid, np.nan, dtype=float)
        mask = np.isclose(grid, x_src[0], atol=1e-12)
        out[mask] = y_src[0]
        return out
    return np.interp(grid, x_src, y_src, left=np.nan, right=np.nan)


def _combo_output_path(solvent_system: str, salt: str) -> Path:
    return MIAC_ROOT / solvent_system / f"{solvent_system}-{salt}.csv"


def _composition_columns(solvent_system: str) -> List[str]:
    solvents = [s for s in solvent_system.split("-") if s]
    if len(solvents) <= 1:
        return []
    cols = []
    for s in solvents:
        cols.append(COMP_COL_LABELS.get(s, f"x_{s}"))
    return cols


def _build_source_records(path: Path, source: str, grouped: dict, write_failures: List[Path]) -> None:
    base_root = MIAC_ROOT if source == "MIAC" else MIAC_M_ROOT
    rel = path.relative_to(base_root)
    solvent_system = rel.parts[0].replace("_", "-").lower()
    salt = _parse_salt_from_stem(path.stem)
    sum_nu = _sum_nu_from_salt(salt)

    fields, rows = _read_rows(path)
    if not rows:
        return

    m_col = _molality_column(fields)
    value_col = _source_value_column(fields, source)
    if m_col is None or value_col is None:
        return

    out_col = "miac_m" if source == "MIAC" else "miac"
    out_fields = _insert_after(fields, value_col, out_col)

    for row in rows:
        m_val = _to_float(row.get(m_col))
        val = _to_float(row.get(value_col))
        if m_val is None or val is None:
            row[out_col] = ""
            continue

        comp, mw_mode = _composition_from_row(row, solvent_system)
        mw_mix = _mw_from_composition(comp, solvent_system)
        factor = 1.0 + mw_mix * m_val * sum_nu

        if source == "MIAC":
            miac = val
            miac_m = miac / factor
            row[out_col] = f"{miac_m:.12g}"
        else:
            miac_m = val
            miac = miac_m * factor
            row[out_col] = f"{miac:.12g}"

        sig = _signature(comp, solvent_system)
        key = (solvent_system, salt, sig)
        bucket = grouped.setdefault(
            key,
            {
                "solvent_system": solvent_system,
                "salt": salt,
                "sum_nu": sum_nu,
                "comp": comp,
                "sources": defaultdict(list),
                "rows_for_table": [],
            },
        )
        bucket["sources"][source].append((m_val, miac_m))
        bucket["rows_for_table"].append(
            {
                "source": source,
                "origin_file": str(path.relative_to(REPO_ROOT)).replace("\\", "/"),
                "molality": m_val,
                "miac": miac,
                "miac_m": miac_m,
                "M_solvent_mix": mw_mix,
                "sum_nu": sum_nu,
                "conversion_factor": factor,
                "mw_mode": mw_mode,
            }
        )

    try:
        _write_rows(path, out_fields, rows)
    except PermissionError:
        write_failures.append(path)


def _build_unified_rows(grouped: dict) -> Dict[Tuple[str, str], List[dict]]:
    unified_by_combo: Dict[Tuple[str, str], List[dict]] = defaultdict(list)

    for key, payload in grouped.items():
        solvent_system, salt, _ = key
        comp = payload["comp"]
        sum_nu = int(payload["sum_nu"])

        p_miac = payload["sources"].get("MIAC", [])
        p_miac_m = payload["sources"].get("MIAC_m", [])

        x1, y1 = _dedupe_points(p_miac)
        x2, y2 = _dedupe_points(p_miac_m)

        grid = np.array(sorted(set(np.round(np.concatenate([x1, x2]), 12))), dtype=float) if (x1.size + x2.size) else np.array([])
        if grid.size == 0:
            continue

        y1i = _interp_to_grid(x1, y1, grid)
        y2i = _interp_to_grid(x2, y2, grid)

        y_u = np.where(np.isfinite(y1i) & np.isfinite(y2i), 0.5 * (y1i + y2i), np.where(np.isfinite(y1i), y1i, y2i))

        mw_mix = _mw_from_composition(comp, solvent_system)
        factor = 1.0 + mw_mix * grid * sum_nu
        miac_u = y_u * factor

        for m, miac, miac_m in zip(grid, miac_u, y_u):
            row = {
                "source": "UNIFIED",
                "origin_file": "",
                "molality": float(m),
                "miac": float(miac),
                "miac_m": float(miac_m),
                "M_solvent_mix": float(mw_mix),
                "sum_nu": sum_nu,
                "conversion_factor": float(1.0 + mw_mix * float(m) * sum_nu),
                "mw_mode": "unified",
                "comp": comp,
            }
            payload["rows_for_table"].append(row)
            unified_by_combo[(solvent_system, salt)].append(row)

    return unified_by_combo


def _write_canonical_miac(unified_by_combo: Dict[Tuple[str, str], List[dict]]) -> List[Path]:
    written: List[Path] = []
    for (solvent_system, salt), rows in sorted(unified_by_combo.items()):
        comp_cols = _composition_columns(solvent_system)
        out_fields = comp_cols + ["molality", "miac", "miac_m"]

        def sort_key(r: dict):
            comp_tuple = tuple(round(float(r["comp"].get(s, 0.0)), 12) for s in [x.split("_")[-1].lower() if x != "x_H2O" else "water" for x in comp_cols])
            return comp_tuple + (float(r["molality"]),)

        rows_sorted = sorted(rows, key=sort_key)
        out_rows = []
        for r in rows_sorted:
            row = {}
            comp = r["comp"]
            for col in comp_cols:
                s = "water" if col == "x_H2O" else col.split("_", 1)[1].lower()
                row[col] = f"{float(comp.get(s, 0.0)):.12g}"
            row["molality"] = f"{float(r['molality']):.12g}"
            row["miac"] = f"{float(r['miac']):.12g}"
            row["miac_m"] = f"{float(r['miac_m']):.12g}"
            out_rows.append(row)

        out_path = _combo_output_path(solvent_system, salt)
        _write_rows(out_path, out_fields, out_rows)
        written.append(out_path)
    return written


def _write_combined_tables(grouped: dict) -> List[Path]:
    COMBINED_ROOT.mkdir(parents=True, exist_ok=True)
    written: List[Path] = []

    for key, payload in sorted(grouped.items(), key=lambda kv: (kv[0][0], kv[0][1], kv[0][2])):
        solvent_system, salt, sig = key
        out_dir = COMBINED_ROOT / solvent_system
        out_dir.mkdir(parents=True, exist_ok=True)
        sig_token = "base" if not sig else "_".join([f"{s}{v:.3f}" for s, v in sig])
        out_path = out_dir / f"{solvent_system}-{salt}-{sig_token}-combined.csv"

        comp_cols = _composition_columns(solvent_system)
        fields = ["source", "origin_file"] + comp_cols + [
            "molality",
            "miac",
            "miac_m",
            "M_solvent_mix",
            "sum_nu",
            "conversion_factor",
            "mw_mode",
        ]

        rows_sorted = sorted(payload["rows_for_table"], key=lambda r: (r["source"], float(r["molality"])))
        out_rows = []
        for r in rows_sorted:
            row = {"source": r["source"], "origin_file": r.get("origin_file", "")}
            comp = payload["comp"]
            for col in comp_cols:
                s = "water" if col == "x_H2O" else col.split("_", 1)[1].lower()
                row[col] = f"{float(comp.get(s, 0.0)):.12g}"
            row["molality"] = f"{float(r['molality']):.12g}"
            row["miac"] = f"{float(r['miac']):.12g}"
            row["miac_m"] = f"{float(r['miac_m']):.12g}"
            row["M_solvent_mix"] = f"{float(r['M_solvent_mix']):.12g}"
            row["sum_nu"] = str(int(r["sum_nu"]))
            row["conversion_factor"] = f"{float(r['conversion_factor']):.12g}"
            row["mw_mode"] = r.get("mw_mode", "")
            out_rows.append(row)

        _write_rows(out_path, fields, out_rows)
        written.append(out_path)

    return written


def _plot_comparisons(grouped: dict) -> List[Path]:
    plot_dir = COMBINED_ROOT / "plots"
    plot_dir.mkdir(parents=True, exist_ok=True)

    by_combo: Dict[Tuple[str, str], List[Tuple[Tuple[str, str, Tuple[Tuple[str, float], ...]], dict]]] = defaultdict(list)
    for key, payload in grouped.items():
        by_combo[(key[0], key[1])].append((key, payload))

    written: List[Path] = []
    for (solvent_system, salt), items in sorted(by_combo.items()):
        fig, ax = plt.subplots(figsize=(7.6, 5.0))

        for idx, (key, payload) in enumerate(sorted(items, key=lambda kv: kv[0][2])):
            _, _, sig = key
            label_suffix = ""
            if sig:
                label_suffix = " | " + ", ".join([f"{s}={v:.3f}" for s, v in sig])

            rows_miac = [r for r in payload["rows_for_table"] if r["source"] == "MIAC"]
            rows_miacm = [r for r in payload["rows_for_table"] if r["source"] == "MIAC_m"]
            rows_u = [r for r in payload["rows_for_table"] if r["source"] == "UNIFIED"]

            if rows_miac:
                x = np.array([float(r["molality"]) for r in sorted(rows_miac, key=lambda r: float(r["molality"]))])
                y = np.array([float(r["miac_m"]) for r in sorted(rows_miac, key=lambda r: float(r["molality"]))])
                ax.plot(x, y, "o", color="tab:blue", markersize=3.5, label=f"MIAC->MIAC_m{label_suffix}")

            if rows_miacm:
                x = np.array([float(r["molality"]) for r in sorted(rows_miacm, key=lambda r: float(r["molality"]))])
                y = np.array([float(r["miac_m"]) for r in sorted(rows_miacm, key=lambda r: float(r["molality"]))])
                ax.plot(x, y, "s", color="tab:orange", markersize=3.5, label=f"MIAC_m source{label_suffix}")

            if rows_u:
                x = np.array([float(r["molality"]) for r in sorted(rows_u, key=lambda r: float(r["molality"]))])
                y = np.array([float(r["miac_m"]) for r in sorted(rows_u, key=lambda r: float(r["molality"]))])
                ax.plot(x, y, "-", color="black", linewidth=1.4, label=f"Unified avg{label_suffix}")

        ax.set_xlabel(r"molality, $m$ / mol kg$^{-1}$")
        ax.set_ylabel(r"$\gamma_{\pm}^{m}$")
        ax.set_title(f"{salt} in {solvent_system}: source vs unified")
        ax.grid(True, alpha=0.3)
        ax.legend(fontsize=7)
        fig.tight_layout()

        out = plot_dir / f"{solvent_system}-{salt}-miac_m-compare.png"
        fig.savefig(out, dpi=220)
        plt.close(fig)
        written.append(out)

    return written



def _cleanup_noncanonical_miac(canonical_paths: List[Path]) -> List[Path]:
    canonical_set = {p.resolve() for p in canonical_paths}
    removed: List[Path] = []

    for path in sorted(MIAC_ROOT.rglob("*.csv")):
        p_res = path.resolve()
        if p_res in canonical_set:
            continue

        try:
            rel = path.relative_to(MIAC_ROOT)
        except ValueError:
            continue
        solvent_system = rel.parts[0].replace("_", "-").lower()
        canonical_salt = _parse_salt_from_stem(path.stem)
        expected = _combo_output_path(solvent_system, canonical_salt).resolve()

        # Remove alias/legacy files only if canonical target exists.
        if expected in canonical_set and expected != p_res:
            path.unlink(missing_ok=True)
            removed.append(path)

    return removed

def _validate_canonical() -> Tuple[int, int]:
    files_checked = 0
    bad_rows = 0

    for path in sorted(MIAC_ROOT.rglob("*.csv")):
        if "plot_" in str(path).lower():
            continue
        fields, rows = _read_rows(path)
        lookup = _field_lookup(fields)
        if "molality" not in lookup or "miac" not in lookup or "miac_m" not in lookup:
            continue

        m_col = lookup["molality"]
        x_col = lookup["miac"]
        mact_col = lookup["miac_m"]

        solvent_system = path.relative_to(MIAC_ROOT).parts[0].replace("_", "-").lower()
        salt = _parse_salt_from_stem(path.stem)
        sum_nu = _sum_nu_from_salt(salt)
        files_checked += 1

        for row in rows:
            m = _to_float(row.get(m_col))
            miac = _to_float(row.get(x_col))
            miac_m = _to_float(row.get(mact_col))
            if m is None or miac is None or miac_m is None:
                continue
            comp, _ = _composition_from_row(row, solvent_system)
            mw = _mw_from_composition(comp, solvent_system)
            recon = miac / (1.0 + mw * m * sum_nu)
            if abs(recon - miac_m) > 1e-8:
                bad_rows += 1

    return files_checked, bad_rows


def main() -> None:
    grouped = {}
    write_failures: List[Path] = []

    for path in sorted(MIAC_ROOT.rglob("*.csv")):
        _build_source_records(path, "MIAC", grouped, write_failures)

    for path in sorted(MIAC_M_ROOT.rglob("*.csv")):
        if "plot_fits" in path.parts or "plot_history" in path.parts:
            continue
        _build_source_records(path, "MIAC_m", grouped, write_failures)

    unified = _build_unified_rows(grouped)
    canonical = _write_canonical_miac(unified)
    removed_aliases = _cleanup_noncanonical_miac(canonical)
    combined = _write_combined_tables(grouped)
    plots = _plot_comparisons(grouped)
    files_checked, bad_rows = _validate_canonical()

    print(f"Updated source CSV files in {MIAC_ROOT} and {MIAC_M_ROOT}")
    if write_failures:
        print("Skipped locked source files:")
        for p in write_failures:
            print(f"- {p}")
    print(f"Wrote {len(canonical)} canonical unified MIAC files into {MIAC_ROOT}")
    if removed_aliases:
        print(f"Removed {len(removed_aliases)} non-canonical alias files from {MIAC_ROOT}")
    print(f"Wrote {len(combined)} comparison tables under {COMBINED_ROOT}")
    print(f"Wrote {len(plots)} comparison plots under {COMBINED_ROOT / 'plots'}")
    print(f"Validation: checked {files_checked} canonical files; mismatched rows={bad_rows}")


if __name__ == "__main__":
    main()
