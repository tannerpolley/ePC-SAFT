from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys

REPO_ROOT = Path(__file__).resolve().parents[2]
SRC_ROOT = REPO_ROOT / "src"
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))
if str(SRC_ROOT) not in sys.path:
    sys.path.insert(0, str(SRC_ROOT))

from scripts.benchmarks.helpers.literature import LITERATURE_CASES
from scripts.benchmarks.helpers.literature import render_literature_benchmark_table
from scripts.benchmarks.helpers.literature import run_literature_benchmarks


def _parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Render the issue #119 literature benchmark and downstream gate registry.")
    parser.add_argument("--case", choices=tuple(LITERATURE_CASES), help="Render only one named literature case.")
    parser.add_argument("--json", type=Path, dest="json_path", help="Write literature-suite output as JSON.")
    parser.add_argument(
        "--registry-only",
        action="store_true",
        help="Only render the issue #119 registry without executing executable benchmark commands.",
    )
    return parser.parse_args()


def main() -> int:
    args = _parse_args()
    payload = run_literature_benchmarks(case=args.case, execute_commands=not args.registry_only)
    print(render_literature_benchmark_table(payload))
    if not args.registry_only:
        summary = payload["execution_summary"]
        print(
            "\nExecution summary: "
            f"attempted={summary['attempted']} "
            f"passed={summary['passed']} "
            f"failed={summary['failed']} "
            f"blocked={summary['blocked']}"
        )
    if args.json_path is not None:
        args.json_path.parent.mkdir(parents=True, exist_ok=True)
        args.json_path.write_text(json.dumps(payload, indent=2), encoding="utf-8")
        print(f"\nWrote JSON: {args.json_path}")
    return 0 if args.registry_only or payload["all_executable_commands_passed"] else 1


if __name__ == "__main__":
    raise SystemExit(main())

