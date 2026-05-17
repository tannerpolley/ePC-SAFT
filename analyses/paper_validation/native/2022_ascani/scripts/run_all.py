from __future__ import annotations

import subprocess
import sys
from pathlib import Path

if __package__ in {None, ""}:
    sys.path.insert(0, str(Path(__file__).resolve().parent))

import _common as common


ROOT = Path(__file__).resolve().parent


def _run(script: Path) -> None:
    print(f"[run] {script}")
    subprocess.run([sys.executable, str(script)], check=True)


def _run_figure_workflows() -> None:
    workflows = (
        ("figure_4b", "plot_figure_4b.py"),
        ("table_5", "plot_table_5.py"),
        ("gibbs_summary", "plot_gibbs_summary.py"),
    )
    figures_root = common.ANALYSIS_DIR / "figures"
    for figure_id, plot_script in workflows:
        figure_scripts = figures_root / figure_id / "scripts"
        _run(figure_scripts / "generate_data.py")
        _run(figure_scripts / plot_script)


def main() -> int:
    rows = common.load_source_rows()
    common.write_normalized_source(rows)
    accepted, solve_payload, _mix, _result = common.solve_payload()
    summary = common.summary_payload(accepted, solve_payload)
    common.write_json(common.SUMMARY_JSON, summary)
    print(common.SUMMARY_JSON.read_text(encoding="utf-8"))
    if accepted:
        _run_figure_workflows()
    return 0 if accepted else 1


if __name__ == "__main__":
    raise SystemExit(main())
