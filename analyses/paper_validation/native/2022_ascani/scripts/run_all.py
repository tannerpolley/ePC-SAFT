from __future__ import annotations

import csv
import json
import sys
from pathlib import Path
from typing import Any

REPO_ROOT = Path(__file__).resolve().parents[5]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from scripts.dev.native_runtime_env import apply_to_current_process

ANALYSIS_DIR = REPO_ROOT / "analyses" / "paper_validation" / "native" / "2022_ascani"
SOURCE_CSV = REPO_ROOT / "data" / "reference" / "multiphase" / "ascani_case2_model_comparison.csv"
PROCESSED_DIR = ANALYSIS_DIR / "data" / "processed"
RESULTS_DIR = ANALYSIS_DIR / "results" / "electrolyte_lle"
NORMALIZED_SOURCE_CSV = PROCESSED_DIR / "source_expected_phase_compositions.csv"
SUMMARY_JSON = RESULTS_DIR / "summary.json"

SPECIES = ["H2O", "Butanol", "Na+", "K+", "Cl-"]
FEED = [
    0.940373242284748,
    0.04879624542603625,
    0.0019339313461782701,
    0.003481324798429627,
    0.005415256144607897,
]


def _rel(path: Path) -> str:
    return path.relative_to(REPO_ROOT).as_posix()


def _load_source_rows() -> list[dict[str, str]]:
    with SOURCE_CSV.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def _write_normalized_source(rows: list[dict[str, str]]) -> None:
    phase_rows: list[dict[str, str]] = []
    mapping = {
        "$x_{water}^{(org)}$": ("org", "H2O"),
        "$x_{butanol}^{(org)}$": ("org", "Butanol"),
        "$x_{NaCl}^{(org)}$": ("org", "NaCl"),
        "$x_{KCl}^{(org)}$": ("org", "KCl"),
        "$x_{water}^{(aq)}$": ("aq", "H2O"),
        "$x_{butanol}^{(aq)}$": ("aq", "Butanol"),
        "$x_{NaCl}^{(aq)}$": ("aq", "NaCl"),
        "$x_{KCl}^{(aq)}$": ("aq", "KCl"),
    }
    for row in rows:
        mapped = mapping.get(row["quantity"])
        if mapped is None:
            continue
        phase, component = mapped
        phase_rows.append(
            {
                "phase": phase,
                "component": component,
                "paper_mole_fraction": row["paper"],
                "model_2020": row["model_2020"],
                "model_2025_num": row["model_2025_num"],
            }
        )
    PROCESSED_DIR.mkdir(parents=True, exist_ok=True)
    with NORMALIZED_SOURCE_CSV.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=["phase", "component", "paper_mole_fraction", "model_2020", "model_2025_num"])
        writer.writeheader()
        writer.writerows(phase_rows)


def _attempt_public_solve() -> tuple[bool, dict[str, Any]]:
    apply_to_current_process()
    import epcsaft
    from epcsaft import ePCSAFTMixture

    mix = ePCSAFTMixture.from_dataset("2022_Ascani", SPECIES, FEED, 298.15)
    options = epcsaft.EquilibriumOptions(max_iterations=180, tolerance=1.0e-8, min_composition=1.0e-12)
    runtime_ipopt = epcsaft.runtime_build_info()["native_dependencies"]["ipopt"]
    try:
        result = mix.equilibrium(kind="electrolyte_lle", T=298.15, P=1.0e5, z=FEED, options=options)
    except epcsaft.SolutionError as exc:
        diagnostics = dict(getattr(exc, "diagnostics", {}) or {})
        return False, {
            "accepted": False,
            "runtime_ipopt": runtime_ipopt,
            "exception_type": type(exc).__name__,
            "exception_message": str(exc),
            "diagnostics": diagnostics,
            "blocker": {
                "kind": "native_ipopt_solver_rejected",
                "route_status": diagnostics.get("route_status"),
                "solver_status": diagnostics.get("solver_status"),
            },
        }
    diagnostics = dict(getattr(result, "diagnostics", {}) or {})
    return True, {
        "accepted": bool(diagnostics.get("accepted", True)),
        "runtime_ipopt": runtime_ipopt,
        "solver_backend": diagnostics.get("solver_backend", diagnostics.get("backend", "ipopt")),
        "diagnostics": diagnostics,
    }


def main() -> int:
    rows = _load_source_rows()
    _write_normalized_source(rows)
    accepted, solve_payload = _attempt_public_solve()
    summary = {
        "status": "accepted" if accepted else "blocked",
        "lane": "ascani_2022_distributed_ion_lle",
        "source_records": [_rel(SOURCE_CSV), _rel(NORMALIZED_SOURCE_CSV)],
        "feed": {"species": SPECIES, "mole_fractions": FEED, "temperature_K": 298.15, "pressure_Pa": 1.0e5},
        "expected": {
            "accepted": True,
            "solver_backend": "ipopt",
            "material_balance_abs": 1.0e-8,
            "charge_balance_abs": 1.0e-8,
        },
        "solve": solve_payload,
    }
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    SUMMARY_JSON.write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(summary, indent=2))
    return 0 if accepted else 1


if __name__ == "__main__":
    raise SystemExit(main())
