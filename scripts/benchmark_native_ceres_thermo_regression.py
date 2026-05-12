from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
SRC_ROOT = REPO_ROOT / "src"
if str(SRC_ROOT) not in sys.path:
    sys.path.insert(0, str(SRC_ROOT))

from epcsaft.benchmarks.native_ceres_thermo_regression import (
    CASE_BUILDERS,
    render_benchmark_table,
    run_native_ceres_thermo_regression_benchmark,
)


def _parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Benchmark the native Ceres thermodynamic regression slice.")
    parser.add_argument("--warmup", type=int, default=1, help="Number of warmup iterations before measurement.")
    parser.add_argument("--repeat", type=int, default=3, help="Number of measured iterations.")
    parser.add_argument("--case", choices=tuple(CASE_BUILDERS), default="reactive_speciation_logk_implicit")
    parser.add_argument("--json", type=Path, dest="json_path", help="Write benchmark output as JSON.")
    return parser.parse_args()


def main() -> int:
    args = _parse_args()
    payload = run_native_ceres_thermo_regression_benchmark(warmup=args.warmup, repeat=args.repeat, case=args.case)
    print(render_benchmark_table(payload))
    if args.json_path is not None:
        args.json_path.parent.mkdir(parents=True, exist_ok=True)
        args.json_path.write_text(json.dumps(payload, indent=2), encoding="utf-8")
        print(f"\nWrote JSON: {args.json_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
