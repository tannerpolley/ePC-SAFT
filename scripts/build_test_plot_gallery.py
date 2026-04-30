from __future__ import annotations

import argparse
import subprocess
import sys
from contextlib import contextmanager
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from scripts.paper_validation.tools import backfill_plotly_html
from scripts.paper_validation.tools import build_analysis_galleries
from scripts.paper_validation.tools import ensure_plot_companions
from scripts.paper_validation.tools import report_plot_companions
from tests.plots import plot_registry


PLOTS_ROOT = REPO_ROOT / "docs" / "plots"


class PlotGalleryBuildError(RuntimeError):
    pass


@dataclass(frozen=True)
class GalleryBuildResult:
    recipe_count: int
    plot_targets: tuple[str, ...]
    output_dirs: tuple[str, ...]
    png_count: int
    plotly_html_created: int = 0
    static_html_created: int = 0
    svg_created: int = 0
    index_path: Path | None = None
    report_path: Path | None = None
    dry_run: bool = False


@contextmanager
def _gallery_root(root: Path):
    old_root = build_analysis_galleries.PLOTS_ROOT
    build_analysis_galleries.PLOTS_ROOT = root
    try:
        yield
    finally:
        build_analysis_galleries.PLOTS_ROOT = old_root


def _plot_data_path(png_path: Path) -> Path:
    return png_path.parent / "data" / f"{png_path.stem}_plot_data.csv"


def _candidate_pngs(plots_root: Path, output_dirs: Iterable[str]) -> tuple[Path, ...]:
    tests_root = plots_root / "tests"
    pngs: list[Path] = []
    for output_dir in output_dirs:
        root = tests_root / output_dir
        if root.exists():
            pngs.extend(path for path in root.rglob("*.png") if "__pycache__" not in path.parts)
    return tuple(sorted(set(pngs), key=lambda path: path.relative_to(plots_root).as_posix().lower()))


def _missing_csvs(pngs: Iterable[Path]) -> list[Path]:
    return [png_path for png_path in pngs if not _plot_data_path(png_path).exists()]


def _can_backfill_plotly(png_path: Path) -> bool:
    return backfill_plotly_html.figure_from_plot_csv(_plot_data_path(png_path), png_path) is not None


def _write_gallery_index(plots_root: Path, *, dry_run: bool) -> Path:
    index_path = plots_root / "index.html"
    if dry_run:
        return index_path
    plots_root.mkdir(parents=True, exist_ok=True)
    with _gallery_root(plots_root):
        pngs = build_analysis_galleries.collect_pngs(plots_root)
        for nested_index in plots_root.rglob("index.html"):
            if nested_index != index_path:
                nested_index.unlink()
        index_path.write_text(build_analysis_galleries.render_gallery_page(plots_root, pngs), encoding="utf-8")
    return index_path


def _write_report(plots_root: Path, *, dry_run: bool) -> Path:
    report_path = plots_root / "plotly_companion_report.csv"
    if dry_run:
        return report_path
    with _gallery_root(plots_root):
        report_plot_companions.write_report(report_path, plots_root)
    return report_path


def _run_plot_producers(plot_targets: tuple[str, ...], *, repo_root: Path, dry_run: bool, skip_pytest: bool) -> None:
    if skip_pytest or dry_run:
        return
    cmd = [sys.executable, "run_pytest.py", *plot_targets, "-q"]
    completed = subprocess.run(cmd, cwd=repo_root, check=False)
    if completed.returncode != 0:
        raise PlotGalleryBuildError(f"Plot producer pytest failed with exit code {completed.returncode}: {' '.join(cmd)}")


def build_gallery(
    recipes: Iterable[plot_registry.PlotRecipe],
    *,
    repo_root: Path = REPO_ROOT,
    plots_root: Path = PLOTS_ROOT,
    dry_run: bool = False,
    force_html: bool = False,
    skip_pytest: bool = False,
) -> GalleryBuildResult:
    selected_recipes = tuple(recipes)
    plot_targets = plot_registry.producer_targets(selected_recipes)
    output_dirs = plot_registry.output_dirs(selected_recipes)

    if not selected_recipes:
        raise PlotGalleryBuildError("No plot recipes were selected.")

    _run_plot_producers(plot_targets, repo_root=repo_root, dry_run=dry_run, skip_pytest=skip_pytest)

    pngs = _candidate_pngs(plots_root, output_dirs)
    missing_csv = _missing_csvs(pngs)
    if missing_csv:
        missing = "\n".join(f"  - {path.as_posix()}" for path in missing_csv)
        raise PlotGalleryBuildError(f"Missing CSV companions for generated PNG plot(s):\n{missing}")

    backfill = backfill_plotly_html.backfill_plotly_html_for_pngs(pngs, force=force_html, dry_run=dry_run)
    static_candidates = tuple(
        path
        for path in pngs
        if (force_html or not path.with_suffix(".html").exists()) and not _can_backfill_plotly(path)
    )
    static = ensure_plot_companions.ensure_png_companions(
        list(static_candidates),
        plots_root,
        dry_run=dry_run,
        create_missing_csv=False,
        force_html=force_html,
    )
    svg_candidates = tuple(path for path in pngs if path not in static_candidates and not path.with_suffix(".svg").exists())
    svg = ensure_plot_companions.ensure_png_companions(
        list(svg_candidates),
        plots_root,
        dry_run=dry_run,
        create_missing_csv=False,
    )

    index_path = _write_gallery_index(plots_root, dry_run=dry_run)
    report_path = _write_report(plots_root, dry_run=dry_run)

    return GalleryBuildResult(
        recipe_count=len(selected_recipes),
        plot_targets=plot_targets,
        output_dirs=output_dirs,
        png_count=len(pngs),
        plotly_html_created=backfill.created,
        static_html_created=static.html_created,
        svg_created=static.svg_created + svg.svg_created,
        index_path=index_path,
        report_path=report_path,
        dry_run=dry_run,
    )


def _format_result(result: GalleryBuildResult) -> str:
    action = "Would build" if result.dry_run else "Built"
    targets = ", ".join(result.plot_targets)
    folders = ", ".join(f"tests/{folder}" for folder in result.output_dirs)
    return (
        f"{action} gallery artifacts for {result.recipe_count} recipe(s).\n"
        f"Plot producers: {targets}\n"
        f"Output folders: {folders}\n"
        f"PNG files: {result.png_count}; Plotly HTML created: {result.plotly_html_created}; "
        f"static HTML created: {result.static_html_created}; SVG created: {result.svg_created}\n"
        f"Gallery index: {result.index_path}\n"
        f"Companion report: {result.report_path}"
    )


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Build gallery-ready plot artifacts for registered test plot recipes.",
    )
    parser.add_argument("targets", nargs="*", help="Source test file, test directory, or pytest node id.")
    parser.add_argument("--all", action="store_true", help="Build every registered test plot recipe.")
    parser.add_argument("--dry-run", action="store_true", help="Report intended work without running pytest or writing files.")
    parser.add_argument("--force-html", action="store_true", help="Regenerate HTML companions even when they already exist.")
    parser.add_argument("--skip-pytest", action="store_true", help="Use existing PNG/CSV outputs without running plot producer tests.")
    parser.add_argument("--plots-root", type=Path, default=PLOTS_ROOT, help=argparse.SUPPRESS)
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    if not args.all and not args.targets:
        print("Provide at least one test target, or use --all.", file=sys.stderr)
        return 2
    try:
        recipes = plot_registry.resolve_plot_recipes(args.targets, all_recipes=args.all)
        result = build_gallery(
            recipes,
            repo_root=REPO_ROOT,
            plots_root=args.plots_root,
            dry_run=args.dry_run,
            force_html=args.force_html,
            skip_pytest=args.skip_pytest,
        )
    except (plot_registry.PlotRecipeLookupError, PlotGalleryBuildError) as exc:
        print(exc, file=sys.stderr)
        return 1
    print(_format_result(result))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
