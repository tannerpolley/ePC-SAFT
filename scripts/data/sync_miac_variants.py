"""Merge/sync MIAC CSV variants (molality, miac_m, mole_fraction, miac).

This script can:
1) Sync missing counterpart values within a canonical MIAC CSV.
2) Merge one or more incoming CSV files (that may contain only one variant) into
   the target CSV without overwriting existing rows.

Example:
  python scripts/data/sync_miac_variants.py \
      --csv data/MIAC/ethanol/ethanol-LiBr.csv \
      --merge C:/Users/Tanner/Downloads/ethanol-LiBr-data2.csv
"""

from __future__ import annotations

import argparse
import csv
import math
from pathlib import Path
from typing import Dict, Iterable, List, Tuple

REPO_ROOT = Path(__file__).resolve().parents[2]

SOLVENT_MW = {
    "water": 18.01528e-3,
    "methanol": 32.04e-3,
    "ethanol": 46.068e-3,
}

COMP_COL_LABELS = {
    "water": "x_H2O",
    "methanol": "x_Methanol",
    "ethanol": "x_Ethanol",
}

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

CATION_CHARGE = {"Li": 1, "Na": 1, "K": 1, "NH4": 1, "H": 1}
ANION_CHARGE = {"Cl": -1, "Br": -1, "I": -1}


def _to_float(value) -> float | None:
    if value is None:
        return None
    text = str(value).strip()
    if not text:
        return None
    try:
        out = float(text)
    except ValueError:
        return None
    if not math.isfinite(out):
        return None
    return out


def _canonical_salt(token: str) -> str:
    key = token.strip().replace("_", "").replace("-", "")
    if not key:
        return token
    upper = key.upper()
    if upper in SALT_CANONICAL:
        return SALT_CANONICAL[upper]

    for cation in ("NH4", "Li", "Na", "K", "H"):
        if upper.startswith(cation.upper()):
            anion = upper[len(cation) :]
            if anion in {"CL", "BR", "I"}:
                anion_fmt = {"CL": "Cl", "BR": "Br", "I": "I"}[anion]
                return f"{cation}{anion_fmt}" if cation != "NH4" else f"NH4{anion_fmt}"
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
        fields = [h.strip() for h in (reader.fieldnames or []) if h and h.strip()]
        rows: List[Dict[str, str]] = []
        for row in reader:
            clean: Dict[str, str] = {}
            for k, v in row.items():
                if not k:
                    continue
                ks = k.strip()
                if not ks:
                    continue
                clean[ks] = v.strip() if isinstance(v, str) else ""
            rows.append(clean)
    return fields, rows


def _write_rows(path: Path, fieldnames: List[str], rows: List[Dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def _lookup(fields: Iterable[str]) -> Dict[str, str]:
    return {f.strip().lower(): f for f in fields if f and f.strip()}


def _solvent_key_from_col(col: str) -> str | None:
    c = col.strip().lower()
    if c in {"x_h2o", "x_water"}:
        return "water"
    if c in {"x_methanol", "x_meoh"}:
        return "methanol"
    if c in {"x_ethanol", "x_etoh"}:
        return "ethanol"
    if c.startswith("x_"):
        tail = c[2:]
        if tail in SOLVENT_MW:
            return tail
    return None


def _extract_comp(row: Dict[str, str], solvent_system: str) -> Dict[str, float]:
    solvents = [s for s in solvent_system.split("-") if s]
    if len(solvents) <= 1:
        return {solvents[0]: 1.0} if solvents else {}

    comp: Dict[str, float] = {}
    for col, val in row.items():
        key = _solvent_key_from_col(col)
        if key is None or key not in solvents:
            continue
        f = _to_float(val)
        if f is not None:
            comp[key] = f

    if len(solvents) == 2 and len(comp) == 1:
        known = next(iter(comp.keys()))
        other = [s for s in solvents if s != known][0]
        comp[other] = 1.0 - comp[known]

    if not comp:
        return {s: 1.0 / len(solvents) for s in solvents}

    for s in solvents:
        comp.setdefault(s, 0.0)
    denom = sum(comp.values())
    if abs(denom) < 1e-12:
        return {s: 1.0 / len(solvents) for s in solvents}
    return {s: comp[s] / denom for s in solvents}


def _mw_mix(solvent_system: str, comp: Dict[str, float]) -> float:
    solvents = [s for s in solvent_system.split("-") if s]
    if not solvents:
        raise ValueError("Empty solvent system")
    if len(solvents) == 1:
        return SOLVENT_MW[solvents[0]]

    total = 0.0
    for s in solvents:
        total += comp.get(s, 0.0) * SOLVENT_MW[s]
    return total


def _composition_columns(solvent_system: str) -> List[str]:
    solvents = [s for s in solvent_system.split("-") if s]
    if len(solvents) <= 1:
        return []
    return [COMP_COL_LABELS.get(s, f"x_{s}") for s in solvents]


def _format(v: float | None) -> str:
    return "" if v is None else f"{float(v):.12g}"


def _row_key(row: Dict[str, str], comp_cols: List[str], include_source: bool) -> Tuple:
    def _f(name: str) -> float:
        val = _to_float(row.get(name, ""))
        return round(val if val is not None else float("nan"), 12)

    comp = tuple(_f(c) for c in comp_cols)
    key = comp + (
        _f("molality"),
        _f("miac_m"),
        _f("mole_fraction"),
        _f("miac"),
    )
    if include_source:
        src = str(row.get("source", "") or "").strip().lower()
        key = key + (src,)
    return key


def _solve_and_normalize_rows(
    rows: List[Dict[str, str]],
    source_fields: List[str],
    solvent_system: str,
    salt: str,
) -> Tuple[List[Dict[str, str]], int]:
    look = _lookup(source_fields)
    m_key = None
    for c in ("molality", "molality (kg/mol)", "m", "m (mol/kg)", "x"):
        if c in look:
            m_key = look[c]
            break
    miac_m_key = next((look[c] for c in ("miac_m", "gamma") if c in look), None)
    miac_key = next((look[c] for c in ("miac", "y") if c in look), None)
    x_key = next((look[c] for c in ("mole_fraction", "x_salt") if c in look), None)
    source_key = next((look[c] for c in ("source",) if c in look), None)

    normalized: List[Dict[str, str]] = []
    skipped = 0
    sum_nu = _sum_nu_from_salt(salt)

    for src in rows:
        comp = _extract_comp(src, solvent_system)
        mw = _mw_mix(solvent_system, comp)

        m = _to_float(src.get(m_key, "")) if m_key else None
        miac_m = _to_float(src.get(miac_m_key, "")) if miac_m_key else None
        miac = _to_float(src.get(miac_key, "")) if miac_key else None
        x = _to_float(src.get(x_key, "")) if x_key else None

        if m is None and x is not None and x < 1.0:
            m = x / (mw * max(1e-15, 1.0 - x))

        if m is None and miac is not None and miac_m is not None and miac_m != 0.0:
            factor = miac / miac_m
            m = (factor - 1.0) / (mw * sum_nu)

        if m is None:
            skipped += 1
            continue

        factor = 1.0 + mw * m * sum_nu

        if x is None:
            denom = m + (1.0 / mw)
            x = m / denom if denom > 0.0 else 0.0

        if miac is None and miac_m is not None:
            miac = miac_m * factor
        if miac_m is None and miac is not None and factor != 0.0:
            miac_m = miac / factor

        if miac is None and miac_m is None:
            skipped += 1
            continue

        out: Dict[str, str] = {}
        for comp_col in _composition_columns(solvent_system):
            s = "water" if comp_col == "x_H2O" else comp_col.split("_", 1)[1].lower()
            out[comp_col] = _format(comp.get(s, 0.0))

        out["molality"] = _format(m)
        out["miac_m"] = _format(miac_m)
        out["mole_fraction"] = _format(x)
        out["miac"] = _format(miac)
        out["source"] = str(src.get(source_key, "") if source_key else "").strip()
        normalized.append(out)

    return normalized, skipped


def merge_and_sync(
    target_csv: Path, merge_csvs: List[Path], solvent_system: str | None = None, salt: str | None = None
) -> Dict[str, int]:
    target_csv = target_csv.resolve()
    if not target_csv.exists():
        raise FileNotFoundError(f"Target CSV not found: {target_csv}")

    if solvent_system is None:
        solvent_system = target_csv.parent.name.lower().replace("_", "-")
    if salt is None:
        salt = _parse_salt_from_stem(target_csv.stem)

    if not salt:
        raise ValueError("Could not infer salt from target CSV stem; pass --salt.")
    if solvent_system not in {
        "water",
        "methanol",
        "ethanol",
        "water-methanol",
        "water-ethanol",
        "water-methanol-ethanol",
    }:
        raise ValueError("Unsupported solvent_system; pass --solvent-system explicitly.")

    target_fields, target_rows = _read_rows(target_csv)
    target_norm, target_skipped = _solve_and_normalize_rows(target_rows, target_fields, solvent_system, salt)

    source_fields_seen = [target_fields]
    merged_norm: List[Dict[str, str]] = list(target_norm)
    incoming_total = 0
    incoming_skipped = 0
    for mpath in merge_csvs:
        fields, rows = _read_rows(mpath)
        source_fields_seen.append(fields)
        incoming_total += len(rows)
        norm, skipped = _solve_and_normalize_rows(rows, fields, solvent_system, salt)
        incoming_skipped += skipped
        merged_norm.extend(norm)

    include_source = any("source" in _lookup(fields) for fields in source_fields_seen)

    comp_cols = _composition_columns(solvent_system)
    unique_rows: List[Dict[str, str]] = []
    seen = set()
    for row in merged_norm:
        key = _row_key(row, comp_cols, include_source=include_source)
        if key in seen:
            continue
        seen.add(key)
        unique_rows.append(row)

    if include_source:
        unique_rows.sort(key=lambda r: (float(r["molality"]), str(r.get("source", "") or "")))
    else:
        unique_rows.sort(key=lambda r: float(r["molality"]))

    out_fields = comp_cols + ["molality", "miac_m", "mole_fraction", "miac"]
    if include_source:
        out_fields.append("source")
    _write_rows(target_csv, out_fields, unique_rows)

    return {
        "target_rows_in": len(target_rows),
        "target_skipped": target_skipped,
        "incoming_rows_in": incoming_total,
        "incoming_skipped": incoming_skipped,
        "rows_out": len(unique_rows),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description="Merge and sync MIAC CSV variants.")
    parser.add_argument("--csv", required=True, help="Target canonical CSV path to update in-place.")
    parser.add_argument("--merge", action="append", default=[], help="Incoming CSV path(s) to merge.")
    parser.add_argument("--solvent-system", default=None, help="Optional solvent system override.")
    parser.add_argument("--salt", default=None, help="Optional salt override.")
    args = parser.parse_args()

    target = Path(args.csv)
    merges = [Path(p) for p in args.merge]

    stats = merge_and_sync(target, merges, solvent_system=args.solvent_system, salt=args.salt)

    print(f"Updated: {target}")
    print(f"- target_rows_in: {stats['target_rows_in']}")
    print(f"- target_skipped: {stats['target_skipped']}")
    print(f"- incoming_rows_in: {stats['incoming_rows_in']}")
    print(f"- incoming_skipped: {stats['incoming_skipped']}")
    print(f"- rows_out: {stats['rows_out']}")


if __name__ == "__main__":
    main()
