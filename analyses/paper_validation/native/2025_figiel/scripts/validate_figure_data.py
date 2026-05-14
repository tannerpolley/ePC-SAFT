from __future__ import annotations

import argparse
import subprocess
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

ROOT = Path(__file__).resolve().parent
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))
if str(ROOT.parent) not in sys.path:
    sys.path.insert(0, str(ROOT.parent))

from shared import figure_data

PLOT_SCRIPTS = [
    ROOT.parent / "figures" / "figure_4" / "scripts" / "plot_figure_4.py",
    ROOT.parent / "figures" / "figure_5" / "scripts" / "plot_figure_5.py",
    ROOT.parent / "figures" / "figure_6" / "scripts" / "plot_figure_6.py",
    ROOT.parent / "figures" / "figure_7" / "scripts" / "plot_figure_7.py",
    ROOT.parent / "figures" / "figure_8" / "scripts" / "plot_figure_8.py",
    ROOT.parent / "figures" / "figure_9" / "scripts" / "plot_figure_9.py",
]


def _run(args: list[str]) -> None:
    print("+ " + " ".join(args))
    subprocess.run(args, cwd=REPO_ROOT, check=True)


def main() -> None:
    parser = argparse.ArgumentParser(description="Validate 2025 Figiel CSV-backed plot payloads.")
    parser.add_argument(
        "--skip-backend", action="store_true", help="Skip native build, doctor, and runtime smoke tests."
    )
    parser.add_argument("--skip-plots", action="store_true", help="Skip PNG/SVG regeneration.")
    parser.add_argument("--rtol", type=float, default=1e-9)
    parser.add_argument("--atol", type=float, default=1e-10)
    args = parser.parse_args()

    if not args.skip_backend:
        _run([sys.executable, "scripts/dev/build_epcsaft.py"])
        _run([sys.executable, "scripts/dev/doctor.py"])
        _run([sys.executable, "run_pytest.py", "tests/api/runtime/test_runtime_exports_and_metadata.py", "-q"])

    failures = figure_data.compare_all(rtol=args.rtol, atol=args.atol)
    if failures:
        print("2025 Figiel figure-data validation failed:")
        for failure in failures[:40]:
            print(f"  - {failure}")
        if len(failures) > 40:
            print(f"  - ... {len(failures) - 40} more")
        raise SystemExit(1)

    if not args.skip_plots:
        for script in PLOT_SCRIPTS:
            _run([sys.executable, str(script.relative_to(REPO_ROOT))])

    print("2025 Figiel CSV-backed figure data validated.")


if __name__ == "__main__":
    main()
