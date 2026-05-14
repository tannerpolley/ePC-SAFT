from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[3] / "src"))

from epcsaft.equilibrium_core.thermo_diagnostics import evaluate_khudaida_solver_gate


def main() -> int:
    parser = argparse.ArgumentParser(description="Print seeded Khudaida electrolyte LLE solver diagnostics.")
    parser.add_argument("--figure", type=int, default=2)
    parser.add_argument("--tie-line", type=int, default=1)
    parser.add_argument("--source", choices=("package", "experimental"), default="package")
    args = parser.parse_args()

    diagnostics = evaluate_khudaida_solver_gate(
        figure=args.figure,
        tie_line=args.tie_line,
        source=args.source,
        seeded=True,
    )
    solver = diagnostics["solver_diagnostics"]
    summary = {
        "fixed_phase_residual_norm": diagnostics["fixed_phase_residual_norm"],
        "solver_residual_norm": solver.get("solver_residual_norm"),
        "gibbs_delta": diagnostics["gibbs_delta"],
        "solver_gibbs_delta": solver.get("gibbs_delta"),
        "phase_distance": solver.get("phase_distance"),
        "decision": diagnostics["decision"],
        "phase_labels_swapped": solver.get("phase_labels_swapped"),
    }
    print(json.dumps(summary, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
