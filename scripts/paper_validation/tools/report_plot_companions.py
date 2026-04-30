from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[3]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from scripts.paper_validation.tools import backfill_plotly_html, build_analysis_galleries


DEFAULT_REPORT = build_analysis_galleries.PLOTS_ROOT / "plotly_companion_report.csv"


def _backfill_status(png_path: Path) -> str:
    data_path = png_path.parent / "data" / f"{png_path.stem}_plot_data.csv"
    if not data_path.exists():
        return "missing_csv"
    reason = backfill_plotly_html._skip_reason(backfill_plotly_html._read_rows(data_path))
    return reason or "backfillable"


def companion_rows(plots_root: Path = build_analysis_galleries.PLOTS_ROOT) -> list[dict[str, str]]:
    pngs = build_analysis_galleries.collect_pngs(plots_root)
    rows: list[dict[str, str]] = []
    for item in build_analysis_galleries.image_manifest(pngs):
        output_path = item["output_path"]
        png_path = plots_root / output_path
        html_path = png_path.with_suffix(".html")
        svg_path = png_path.with_suffix(".svg")
        data_path = png_path.parent / "data" / f"{png_path.stem}_plot_data.csv"
        interactive_source = item["interactive_source"]
        has_html = html_path.exists()
        has_svg = svg_path.exists()
        has_csv = data_path.exists()
        rows.append(
            {
                "output_path": output_path,
                "source_path": item["source_path"],
                "interactive_source": interactive_source,
                "html_status": "missing_html"
                if not has_html
                else ("static_html" if interactive_source == "static_png_wrapper" else "interactive"),
                "interactive_status": "interactive" if has_html and interactive_source != "static_png_wrapper" else "static_only",
                "html_path": item["html_path"],
                "svg_path": item["svg_path"],
                "data_path": item["data_path"],
                "has_html": str(has_html).lower(),
                "has_plotly_html": str(has_html and interactive_source != "static_png_wrapper").lower(),
                "has_svg": str(has_svg).lower(),
                "has_csv": str(has_csv).lower(),
                "bundle_complete": str(has_html and has_svg and has_csv).lower(),
                "static_reason": _backfill_status(png_path) if interactive_source == "static_png_wrapper" or not has_html else "",
            }
        )
    return rows


def write_report(output: Path, plots_root: Path = build_analysis_galleries.PLOTS_ROOT) -> Path:
    rows = companion_rows(plots_root)
    output.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = [
        "output_path",
        "source_path",
        "interactive_source",
        "html_status",
        "interactive_status",
        "html_path",
        "svg_path",
        "data_path",
        "has_html",
        "has_plotly_html",
        "has_svg",
        "has_csv",
        "bundle_complete",
        "static_reason",
    ]
    with output.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)
    return output


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Report Plotly/SVG/CSV companions for docs/plots PNG files.")
    parser.add_argument("--root", type=Path, default=build_analysis_galleries.PLOTS_ROOT, help="Plot gallery root.")
    parser.add_argument("--output", type=Path, default=DEFAULT_REPORT, help="CSV report path.")
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    rows = companion_rows(args.root)
    write_report(args.output, args.root)
    interactive = sum(1 for row in rows if row["interactive_status"] == "interactive")
    static_html = sum(1 for row in rows if row["html_status"] == "static_html")
    svg = sum(1 for row in rows if row["has_svg"] == "true")
    complete = sum(1 for row in rows if row["bundle_complete"] == "true")
    print(
        f"Wrote {args.output} with {len(rows)} PNG plot(s): "
        f"{interactive} interactive HTML, {static_html} static HTML, {svg} SVG companions, "
        f"{complete} complete bundles."
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
