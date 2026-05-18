from __future__ import annotations

import sys
from pathlib import Path

if __package__ in {None, ""}:
    sys.path.insert(0, str(Path(__file__).resolve().parent))

import _common as common


ROOT = Path(__file__).resolve().parent


def main() -> int:
    rows = common.load_source_rows()
    common.write_normalized_source(rows)
    accepted, solve_payload, mix, result = common.solve_payload()
    summary = common.summary_payload(accepted, solve_payload, result)
    common.write_json(common.SUMMARY_JSON, summary)
    common.write_retained_outputs(summary, mix, result)
    print(common.SUMMARY_JSON.read_text(encoding="utf-8"))
    return 0 if accepted else 1


if __name__ == "__main__":
    raise SystemExit(main())
