"""Scan alternative Figure 6b accounting methods against digitized contribution lines."""

from __future__ import annotations

import argparse
import csv
import math
import sys
from pathlib import Path
from typing import Dict, Iterable, List, Tuple

CONTRIBUTIONS = ["born", "dh", "hc", "disp", "assoc"]
REPO_ROOT = Path(__file__).resolve().parents[6]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from scripts.plot_outputs import paper_validation_path

OUTPUT_ROOT = paper_validation_path(Path(__file__).resolve().parent, "output")
OUTPUT_DATA_DIR = OUTPUT_ROOT / "data"


def _read_csv_rows(path: Path) -> list[list[str]]:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        return list(csv.reader(handle))


def _interp(x_grid: list[float], y_grid: list[float], x: float) -> float:
    if x <= x_grid[0]:
        return y_grid[0]
    if x >= x_grid[-1]:
        return y_grid[-1]
    lo = 0
    hi = len(x_grid) - 1
    while hi - lo > 1:
        mid = (lo + hi) // 2
        if x_grid[mid] <= x:
            lo = mid
        else:
            hi = mid
    t = (x - x_grid[lo]) / (x_grid[hi] - x_grid[lo])
    return y_grid[lo] + t * (y_grid[hi] - y_grid[lo])


def _load_bookkeeping(path: Path) -> tuple[list[float], dict[str, list[float]]]:
    rows = _read_csv_rows(path)
    header = rows[0]
    idx = {name: i for i, name in enumerate(header)}
    needed = ["x_salt", "z_correction"] + CONTRIBUTIONS
    for key in needed:
        if key not in idx:
            raise ValueError(f"Missing column '{key}' in bookkeeping CSV {path}")
    x = [float(row[idx["x_salt"]]) for row in rows[1:]]
    series = {key: [float(row[idx[key]]) for row in rows[1:]] for key in needed if key != "x_salt"}
    return x, series


def _load_digitized(path: Path) -> dict[str, tuple[list[float], list[float]]]:
    rows = _read_csv_rows(path)
    header = rows[0]
    aliases = {
        "born": "born",
        "dh": "dh",
        "hc": "hc",
        "hard-chain": "hc",
        "disp": "disp",
        "dispersion": "disp",
        "assoc": "assoc",
        "association": "assoc",
        "total": "total",
        "data": "data",
    }
    out: dict[str, tuple[list[float], list[float]]] = {}
    for col in range(0, len(header), 2):
        name = aliases[header[col + 1].strip().lower()]
        xs: list[float] = []
        ys: list[float] = []
        for row in rows[1:]:
            if col < len(row) and col + 1 < len(row) and row[col].strip() and row[col + 1].strip():
                xs.append(float(row[col]))
                ys.append(float(row[col + 1]))
        out[name] = (xs, ys)
    return out


def _rmse(x_model: list[float], y_model: list[float], x_data: list[float], y_data: list[float]) -> float:
    err2 = 0.0
    for x, y in zip(x_data, y_data):
        d = _interp(x_model, y_model, x) - y
        err2 += d * d
    return math.sqrt(err2 / len(x_data))


def _apply_weights(base: dict[str, list[float]], weights: dict[str, float]) -> dict[str, list[float]]:
    z = base["z_correction"]
    out: dict[str, list[float]] = {}
    for key in CONTRIBUTIONS:
        w = float(weights.get(key, 0.0))
        out[key] = [b + w * zz for b, zz in zip(base[key], z)]
    return out


def _best_scalar_weight(x_model: list[float], base: list[float], zcorr: list[float], x_data: list[float], y_data: list[float]) -> tuple[float, float, float]:
    num = 0.0
    den = 0.0
    for x, y in zip(x_data, y_data):
        b = _interp(x_model, base, x)
        z = _interp(x_model, zcorr, x)
        num += z * (y - b)
        den += z * z
    w = 0.0 if den == 0.0 else num / den
    y_model = [b + w * z for b, z in zip(base, zcorr)]
    return w, _rmse(x_model, base, x_data, y_data), _rmse(x_model, y_model, x_data, y_data)


def run_scan(bookkeeping_csv: Path, digitized_csv: Path, summary_csv: Path, best_weights_csv: Path) -> None:
    x_model, base = _load_bookkeeping(bookkeeping_csv)
    digitized = _load_digitized(digitized_csv)

    methods: list[tuple[str, dict[str, float]]] = [
        ("mu_baseline", {}),
        ("add_z_to_dh", {"dh": 1.0}),
        ("add_z_to_hc", {"hc": 1.0}),
        ("add_z_to_disp", {"disp": 1.0}),
        ("add_z_to_assoc", {"assoc": 1.0}),
        ("equal_split_nonborn", {"dh": 0.25, "hc": 0.25, "disp": 0.25, "assoc": 0.25}),
    ]

    # Coarse constrained search: one shared z-allocation across DH/HC/disp/assoc, sum(weights)=1.
    best_grid_name = "best_grid_sum1"
    best_grid_weights: dict[str, float] | None = None
    best_grid_score = math.inf
    best_grid_rmses: dict[str, float] | None = None
    for i in range(21):
        for j in range(21 - i):
            for k in range(21 - i - j):
                l = 20 - i - j - k
                weights = {"dh": i / 20.0, "hc": j / 20.0, "disp": k / 20.0, "assoc": l / 20.0}
                adjusted = _apply_weights(base, weights)
                rmses = {key: _rmse(x_model, adjusted[key], *digitized[key]) for key in CONTRIBUTIONS}
                score = sum(rmses.values())
                if score < best_grid_score:
                    best_grid_score = score
                    best_grid_weights = weights
                    best_grid_rmses = rmses

    summary_rows: list[dict[str, float | str]] = []
    for name, weights in methods:
        adjusted = _apply_weights(base, weights)
        row: dict[str, float | str] = {"method": name}
        for key in CONTRIBUTIONS:
            row[f"{key}_weight"] = float(weights.get(key, 0.0))
            row[f"{key}_rmse"] = _rmse(x_model, adjusted[key], *digitized[key])
        row["score_sum_rmse"] = sum(float(row[f"{key}_rmse"]) for key in CONTRIBUTIONS)
        summary_rows.append(row)

    if best_grid_weights is not None and best_grid_rmses is not None:
        row = {"method": best_grid_name}
        for key in CONTRIBUTIONS:
            row[f"{key}_weight"] = float(best_grid_weights.get(key, 0.0))
            row[f"{key}_rmse"] = float(best_grid_rmses[key])
        row["score_sum_rmse"] = sum(float(row[f"{key}_rmse"]) for key in CONTRIBUTIONS)
        summary_rows.append(row)

    best_weight_rows: list[dict[str, float | str]] = []
    for key in CONTRIBUTIONS:
        w, rmse0, rmse_best = _best_scalar_weight(x_model, base[key], base["z_correction"], *digitized[key])
        best_weight_rows.append(
            {
                "contribution": key,
                "best_scalar_z_weight": w,
                "baseline_rmse": rmse0,
                "best_rmse": rmse_best,
            }
        )

    summary_csv.parent.mkdir(parents=True, exist_ok=True)
    with summary_csv.open("w", newline="", encoding="utf-8") as handle:
        fieldnames = ["method"] + [f"{key}_weight" for key in CONTRIBUTIONS] + [f"{key}_rmse" for key in CONTRIBUTIONS] + ["score_sum_rmse"]
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(summary_rows)

    with best_weights_csv.open("w", newline="", encoding="utf-8") as handle:
        fieldnames = ["contribution", "best_scalar_z_weight", "baseline_rmse", "best_rmse"]
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(best_weight_rows)

    print(f"Wrote accounting summary: {summary_csv}")
    print(f"Wrote per-series best z-weight fits: {best_weights_csv}")


def _parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Scan alternative Figure 6b accounting methods")
    parser.add_argument(
        "--bookkeeping-csv",
        type=Path,
        default=OUTPUT_DATA_DIR / "figure6b_bookkeeping.csv",
    )
    parser.add_argument(
        "--digitized-csv",
        type=Path,
        default=Path(r"C:\Users\Tanner\Documents\git\ePC-SAFT\scripts\paper_validation\2020_Bulow_analysis\figure_6\figure_6b\data\Figure6b_curves.csv"),
    )
    parser.add_argument(
        "--summary-csv",
        type=Path,
        default=OUTPUT_DATA_DIR / "figure6b_accounting_method_scan.csv",
    )
    parser.add_argument(
        "--best-weights-csv",
        type=Path,
        default=OUTPUT_DATA_DIR / "figure6b_accounting_best_scalar_zfits.csv",
    )
    return parser.parse_args()


def main() -> None:
    args = _parse_args()
    run_scan(
        bookkeeping_csv=Path(args.bookkeeping_csv),
        digitized_csv=Path(args.digitized_csv),
        summary_csv=Path(args.summary_csv),
        best_weights_csv=Path(args.best_weights_csv),
    )


if __name__ == "__main__":
    main()

