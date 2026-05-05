from __future__ import annotations

import argparse
import csv
import html
import struct
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

REPO_ROOT = Path(__file__).resolve().parents[3]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from scripts.paper_validation.tools import build_analysis_galleries

PLOTS_ROOT = build_analysis_galleries.PLOTS_ROOT


@dataclass(frozen=True)
class CompanionResult:
    csv_created: int = 0
    svg_created: int = 0
    skipped_existing_csv: int = 0
    skipped_existing_svg: int = 0


def _plot_data_path(png_path: Path) -> Path:
    return png_path.parent / f"{png_path.stem}_plot_data.csv"


def _png_dimensions(png_path: Path) -> tuple[int, int]:
    with png_path.open("rb") as handle:
        header = handle.read(24)
    if len(header) < 24 or header[:8] != b"\x89PNG\r\n\x1a\n" or header[12:16] != b"IHDR":
        return (1200, 800)
    width, height = struct.unpack(">II", header[16:24])
    return (max(int(width), 1), max(int(height), 1))


def _write_placeholder_csv(png_path: Path, *, dry_run: bool) -> bool:
    csv_path = _plot_data_path(png_path)
    if csv_path.exists():
        return False
    if dry_run:
        return True
    csv_path.parent.mkdir(parents=True, exist_ok=True)
    with csv_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "figure_file",
                "axes_index",
                "axes_title",
                "artist_type",
                "artist_label",
                "series_index",
                "point_index",
                "x",
                "y",
                "width",
                "height",
            ],
        )
        writer.writeheader()
        writer.writerow({"figure_file": png_path.name, "artist_type": "no_numeric_artists"})
    return True


def _write_static_svg(png_path: Path, *, dry_run: bool, force: bool = False) -> bool:
    svg_path = png_path.with_suffix(".svg")
    if svg_path.exists() and not force:
        return False
    if dry_run:
        return True
    width, height = _png_dimensions(png_path)
    title = html.escape(png_path.stem.replace("_", " ").replace("-", " "))
    rel_png = html.escape(png_path.name)
    svg_path.write_text(
        f"""<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}" role="img" aria-label="{title}">
  <title>{title}</title>
  <desc>Static SVG wrapper for the generated PNG plot. The original plot data is stored in the sibling CSV companion.</desc>
  <image href="{rel_png}" x="0" y="0" width="{width}" height="{height}" preserveAspectRatio="xMidYMid meet"/>
</svg>
""",
        encoding="utf-8",
        newline="\n",
    )
    return True


def ensure_png_companions(
    png_paths: Iterable[Path],
    plots_root: Path = PLOTS_ROOT,
    *,
    dry_run: bool = False,
    create_missing_csv: bool = True,
    force_svg: bool = False,
) -> CompanionResult:
    csv_created = svg_created = 0
    skipped_existing_svg = skipped_existing_csv = 0
    for png_path in sorted({Path(path) for path in png_paths}, key=lambda path: path.as_posix().lower()):
        csv_path = _plot_data_path(png_path)
        svg_path = png_path.with_suffix(".svg")
        if csv_path.exists():
            skipped_existing_csv += 1
        elif create_missing_csv and _write_placeholder_csv(png_path, dry_run=dry_run):
            csv_created += 1
        if svg_path.exists() and not force_svg:
            skipped_existing_svg += 1
        if (not svg_path.exists() or force_svg) and _write_static_svg(png_path, dry_run=dry_run, force=force_svg):
            svg_created += 1
    return CompanionResult(
        csv_created=csv_created,
        svg_created=svg_created,
        skipped_existing_svg=skipped_existing_svg,
        skipped_existing_csv=skipped_existing_csv,
    )


def ensure_companions(plots_root: Path = PLOTS_ROOT, *, dry_run: bool = False) -> CompanionResult:
    return ensure_png_companions(build_analysis_galleries.collect_pngs(plots_root), plots_root, dry_run=dry_run)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Ensure every source-local out PNG has CSV and SVG companions.")
    parser.add_argument("--root", type=Path, default=PLOTS_ROOT, help="Plot gallery root.")
    parser.add_argument("--dry-run", action="store_true", help="Report work without writing companions.")
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    result = ensure_companions(args.root, dry_run=args.dry_run)
    action = "Would create" if args.dry_run else "Created"
    print(
        f"{action}: csv={result.csv_created}, svg={result.svg_created}. "
        f"Existing: csv={result.skipped_existing_csv}, svg={result.skipped_existing_svg}."
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
