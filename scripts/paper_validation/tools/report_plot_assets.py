from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[3]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from scripts.paper_validation.tools import build_analysis_galleries

DEFAULT_REPORT = REPO_ROOT / "build" / "plot_gallery" / "plot_asset_report.csv"


def asset_rows(plots_root: Path = build_analysis_galleries.PLOTS_ROOT) -> list[dict[str, str]]:
    rows: list[dict[str, str]] = []
    for item in build_analysis_galleries.gallery_entries(plots_root):
        output_path = item["output_path"]
        png_path = build_analysis_galleries.gallery_repo_root() / output_path
        svg_path = png_path.with_suffix(".svg")
        data_path = png_path.parent / f"{png_path.stem}_plot_data.csv"
        has_svg = svg_path.exists()
        has_csv = data_path.exists()
        rows.append(
            {
                "output_path": output_path,
                "source_path": item["source_path"],
                "svg_path": item["svg_path"],
                "data_path": item["data_path"],
                "has_svg": str(has_svg).lower(),
                "has_csv": str(has_csv).lower(),
                "bundle_complete": str(has_svg and has_csv).lower(),
            }
        )
    return rows


def write_report(output: Path, plots_root: Path = build_analysis_galleries.PLOTS_ROOT) -> Path:
    rows = asset_rows(plots_root)
    output.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = [
        "output_path",
        "source_path",
        "svg_path",
        "data_path",
        "has_svg",
        "has_csv",
        "bundle_complete",
    ]
    with output.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)
    return output


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Report PNG/SVG/CSV static plot assets under source-local out folders."
    )
    parser.add_argument("--root", type=Path, default=build_analysis_galleries.PLOTS_ROOT, help="Plot gallery root.")
    parser.add_argument("--output", type=Path, default=DEFAULT_REPORT, help="CSV report path.")
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    rows = asset_rows(args.root)
    write_report(args.output, args.root)
    svg = sum(1 for row in rows if row["has_svg"] == "true")
    csv_count = sum(1 for row in rows if row["has_csv"] == "true")
    complete = sum(1 for row in rows if row["bundle_complete"] == "true")
    print(f"Wrote {args.output} with {len(rows)} PNG plot(s): {svg} SVG, {csv_count} CSV, {complete} complete bundles.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
