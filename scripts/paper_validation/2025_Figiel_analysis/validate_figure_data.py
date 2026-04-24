from __future__ import annotations

import argparse
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent
REPO_ROOT = ROOT.parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

import figure_data

PLOT_SCRIPTS = [
    ROOT / "figure_4" / "plot_figure_4.py",
    ROOT / "figure_5" / "plot_figure_5.py",
    ROOT / "figure_6" / "plot_figure_6.py",
    ROOT / "figure_7" / "plot_figure_7.py",
    ROOT / "figure_8" / "plot_figure_8.py",
    ROOT / "figure_9" / "plot_figure_9.py",
]


def _run(args: list[str]) -> None:
    print("+ " + " ".join(args))
    subprocess.run(args, cwd=REPO_ROOT, check=True)


def main() -> None:
    parser = argparse.ArgumentParser(description="Validate 2025 Figiel CSV-backed plot payloads.")
    parser.add_argument("--skip-backend", action="store_true", help="Skip native build, doctor, and runtime smoke tests.")
    parser.add_argument("--skip-plots", action="store_true", help="Skip PNG/gallery regeneration.")
    parser.add_argument("--rtol", type=float, default=1e-9)
    parser.add_argument("--atol", type=float, default=1e-10)
    args = parser.parse_args()

    if not args.skip_backend:
        _run([sys.executable, "scripts/build_epcsaft.py"])
        _run([sys.executable, "scripts/codex_doctor.py"])
        _run([sys.executable, "run_pytest.py", "tests/test_runtime.py", "-q"])

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
        _run([sys.executable, "scripts/paper_validation/tools/build_analysis_galleries.py"])

    print("2025 Figiel CSV-backed figure data validated.")


if __name__ == "__main__":
    main()
