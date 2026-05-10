from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys

REPO_ROOT = Path(__file__).resolve().parents[1]
SRC_ROOT = REPO_ROOT / "src"
if str(SRC_ROOT) not in sys.path:
    sys.path.insert(0, str(SRC_ROOT))

from epcsaft.benchmarks.neutral_equilibrium import CASE_BUILDERS
from epcsaft.benchmarks.neutral_equilibrium import render_benchmark_table
from epcsaft.benchmarks.neutral_equilibrium import run_neutral_equilibrium_benchmarks


def _parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Benchmark neutral equilibrium workflows without requiring FeOs.")
    parser.add_argument("--warmup", type=int, default=20, help="Number of warmup iterations before measurement.")
    parser.add_argument("--repeat", type=int, default=100, help="Number of measured iterations per case.")
    parser.add_argument("--case", choices=tuple(CASE_BUILDERS), help="Run only one named benchmark case.")
    parser.add_argument("--json", type=Path, dest="json_path", help="Write benchmark output as JSON.")
    parser.add_argument("--baseline-json", type=Path, help="Optional baseline JSON to compare medians against.")
    return parser.parse_args()


def main() -> int:
    args = _parse_args()
    payload = run_neutral_equilibrium_benchmarks(
        warmup=args.warmup,
        repeat=args.repeat,
        case=args.case,
        baseline_json=args.baseline_json,
    )
    print(render_benchmark_table(payload))
    if args.json_path is not None:
        args.json_path.parent.mkdir(parents=True, exist_ok=True)
        args.json_path.write_text(json.dumps(payload, indent=2), encoding="utf-8")
        print(f"\nWrote JSON: {args.json_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
