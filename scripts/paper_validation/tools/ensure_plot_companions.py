from __future__ import annotations

import argparse
import csv
import html
import struct
import sys
from dataclasses import dataclass
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[3]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from scripts.paper_validation.tools import build_analysis_galleries


PLOTS_ROOT = build_analysis_galleries.PLOTS_ROOT
STATIC_WRAPPER_MARKER = 'epcsaft-interactive-source="static_png_wrapper"'


@dataclass(frozen=True)
class CompanionResult:
    csv_created: int = 0
    html_created: int = 0
    svg_created: int = 0
    skipped_existing_html: int = 0
    skipped_existing_svg: int = 0
    skipped_existing_csv: int = 0


def _plot_data_path(png_path: Path) -> Path:
    return png_path.parent / "data" / f"{png_path.stem}_plot_data.csv"


def _png_dimensions(png_path: Path) -> tuple[int, int]:
    with png_path.open("rb") as handle:
        header = handle.read(24)
    if len(header) < 24 or header[:8] != b"\x89PNG\r\n\x1a\n" or header[12:16] != b"IHDR":
        return (1200, 800)
    width, height = struct.unpack(">II", header[16:24])
    return (max(int(width), 1), max(int(height), 1))


def _relative_posix(path: Path, start: Path) -> str:
    return path.relative_to(start).as_posix()


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
        writer.writerow(
            {
                "figure_file": png_path.name,
                "artist_type": "no_numeric_artists",
            }
        )
    return True


def _write_static_html(png_path: Path, *, plots_root: Path, dry_run: bool) -> bool:
    html_path = png_path.with_suffix(".html")
    if html_path.exists():
        return False
    if dry_run:
        return True
    rel_png = html.escape(png_path.name)
    csv_path = _plot_data_path(png_path)
    rel_csv = html.escape(_relative_posix(csv_path, png_path.parent)) if csv_path.exists() else ""
    rel_svg = html.escape(png_path.with_suffix(".svg").name)
    title = html.escape(png_path.stem.replace("_", " ").replace("-", " "))
    source_path = html.escape(_relative_posix(png_path, plots_root))
    csv_link = f'<a href="{rel_csv}">CSV</a>' if rel_csv else ""
    html_path.write_text(
        f"""<!doctype html>
<html lang="en">
<head>
  <meta charset="utf-8">
  <meta name="epcsaft-interactive-source" content="static_png_wrapper">
  <meta name="viewport" content="width=device-width, initial-scale=1">
  <title>{title}</title>
  <style>
    body {{ margin: 0; font-family: system-ui, -apple-system, Segoe UI, sans-serif; background: #f8fafc; color: #162033; }}
    main {{ box-sizing: border-box; min-height: 100vh; padding: 18px; display: grid; gap: 12px; }}
    header {{ display: flex; flex-wrap: wrap; gap: 10px; align-items: baseline; }}
    h1 {{ margin: 0; font-size: 18px; line-height: 1.2; }}
    p {{ margin: 0; color: #53627a; font-size: 13px; }}
    nav {{ display: flex; gap: 8px; flex-wrap: wrap; }}
    a {{ color: #0b64c0; text-decoration: none; border: 1px solid #cfe0f2; border-radius: 6px; padding: 4px 8px; background: white; }}
    .frame {{ min-height: 0; display: grid; place-items: center; }}
    img {{ max-width: 100%; max-height: calc(100vh - 112px); object-fit: contain; background: white; border: 1px solid #d8e2ee; border-radius: 6px; }}
  </style>
</head>
<body>
  <!-- {STATIC_WRAPPER_MARKER} -->
  <main>
    <header>
      <h1>{title}</h1>
      <p>{source_path}</p>
      <nav><a href="{rel_png}">PNG</a><a href="{rel_svg}">SVG</a>{csv_link}</nav>
    </header>
    <div class="frame"><img src="{rel_png}" alt="{title}"></div>
  </main>
</body>
</html>
""",
        encoding="utf-8",
        newline="\n",
    )
    return True


def _write_static_svg(png_path: Path, *, dry_run: bool) -> bool:
    svg_path = png_path.with_suffix(".svg")
    if svg_path.exists():
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


def ensure_companions(plots_root: Path = PLOTS_ROOT, *, dry_run: bool = False) -> CompanionResult:
    csv_created = html_created = svg_created = 0
    skipped_existing_html = skipped_existing_svg = skipped_existing_csv = 0
    for png_path in build_analysis_galleries.collect_pngs(plots_root):
        csv_path = _plot_data_path(png_path)
        html_path = png_path.with_suffix(".html")
        svg_path = png_path.with_suffix(".svg")
        if csv_path.exists():
            skipped_existing_csv += 1
        elif _write_placeholder_csv(png_path, dry_run=dry_run):
            csv_created += 1
        if html_path.exists():
            skipped_existing_html += 1
        elif _write_static_html(png_path, plots_root=plots_root, dry_run=dry_run):
            html_created += 1
        if svg_path.exists():
            skipped_existing_svg += 1
        elif _write_static_svg(png_path, dry_run=dry_run):
            svg_created += 1
    return CompanionResult(
        csv_created=csv_created,
        html_created=html_created,
        svg_created=svg_created,
        skipped_existing_html=skipped_existing_html,
        skipped_existing_svg=skipped_existing_svg,
        skipped_existing_csv=skipped_existing_csv,
    )


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Ensure every docs/plots PNG has CSV, HTML, and SVG companions.")
    parser.add_argument("--root", type=Path, default=PLOTS_ROOT, help="Plot gallery root.")
    parser.add_argument("--dry-run", action="store_true", help="Report work without writing companions.")
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    result = ensure_companions(args.root, dry_run=args.dry_run)
    action = "Would create" if args.dry_run else "Created"
    print(
        f"{action}: csv={result.csv_created}, html={result.html_created}, svg={result.svg_created}. "
        f"Existing: csv={result.skipped_existing_csv}, html={result.skipped_existing_html}, svg={result.skipped_existing_svg}."
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
