from __future__ import annotations

import argparse
import json
import sys


from pathlib import Path
import sys as _bootstrap_sys
from pathlib import Path as _BootstrapPath

for _candidate in _BootstrapPath(__file__).resolve().parents:
    if (_candidate / "scripts" / "plot_outputs.py").is_file():
        if str(_candidate) not in _bootstrap_sys.path:
            _bootstrap_sys.path.insert(0, str(_candidate))
        break
else:
    raise ModuleNotFoundError("Could not locate repo root containing scripts/plot_outputs.py")
from scripts.plot_outputs import REPO_ROOT

sys.path.insert(0, str(REPO_ROOT / "src"))

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
