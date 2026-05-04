from __future__ import annotations

import argparse
import json
import subprocess
import sys
from contextlib import contextmanager
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from scripts import build_plot_manifest
from scripts.paper_validation.tools import build_analysis_galleries
from scripts.paper_validation.tools import ensure_plot_companions
from scripts.paper_validation.tools import render_plot_data_csv
from scripts.paper_validation.tools import report_plot_assets
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
    rendered_from_csv: int = 0
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
    return png_path.parent / f"{png_path.stem}_plot_data.csv"


def _candidate_pngs(plots_root: Path, output_dirs: Iterable[str], *, repo_root: Path = REPO_ROOT) -> tuple[Path, ...]:
    tests_root = repo_root / "tests" / "plots" / "out"
    pngs: list[Path] = []
    for output_dir in output_dirs:
        root = tests_root / output_dir
        if root.exists():
            pngs.extend(path for path in root.rglob("*.png") if "__pycache__" not in path.parts)
    return tuple(sorted(set(pngs), key=lambda path: path.relative_to(repo_root).as_posix().lower()))


def _missing_csvs(pngs: Iterable[Path]) -> list[Path]:
    return [png_path for png_path in pngs if not _plot_data_path(png_path).exists()]


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


def _write_report(plots_root: Path, *, repo_root: Path = REPO_ROOT, dry_run: bool) -> Path:
    report_path = repo_root / "build" / "plot_gallery" / "plot_asset_report.csv"
    if dry_run:
        return report_path
    with _gallery_root(plots_root):
        report_plot_assets.write_report(report_path, plots_root)
    return report_path


def _write_manifest(plots_root: Path, *, dry_run: bool) -> Path:
    manifest_path = plots_root / "manifest.json"
    if dry_run:
        return manifest_path
    with _gallery_root(plots_root):
        payload = build_plot_manifest.manifest_payload(
            build_analysis_galleries.image_manifest(build_analysis_galleries.collect_pngs(plots_root))
        )
    manifest_path.parent.mkdir(parents=True, exist_ok=True)
    manifest_path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return manifest_path


def _run_plot_producers(plot_targets: tuple[str, ...], *, repo_root: Path, dry_run: bool, skip_pytest: bool) -> None:
    if skip_pytest or dry_run:
        return
    cmd = [sys.executable, "run_pytest.py", *plot_targets, "-q"]
    completed = subprocess.run(cmd, cwd=repo_root, check=False)
    if completed.returncode != 0:
        raise PlotGalleryBuildError(
            f"Plot producer pytest failed with exit code {completed.returncode}: {' '.join(cmd)}"
        )


def _render_missing_static_assets(pngs: Iterable[Path], *, dry_run: bool, force_render: bool) -> int:
    rendered = 0
    for png_path in pngs:
        csv_path = _plot_data_path(png_path)
        if force_render or not png_path.exists() or not png_path.with_suffix(".svg").exists():
            rendered += 1
            if not dry_run:
                render_plot_data_csv.render_csv_to_static_assets(csv_path, png_path)
    return rendered


def build_gallery(
    recipes: Iterable[plot_registry.PlotRecipe],
    *,
    repo_root: Path = REPO_ROOT,
    plots_root: Path = PLOTS_ROOT,
    dry_run: bool = False,
    force_render: bool = False,
    skip_pytest: bool = False,
) -> GalleryBuildResult:
    selected_recipes = tuple(recipes)
    plot_targets = plot_registry.producer_targets(selected_recipes)
    output_dirs = plot_registry.output_dirs(selected_recipes)

    if not selected_recipes:
        raise PlotGalleryBuildError("No plot recipes were selected.")

    _run_plot_producers(plot_targets, repo_root=repo_root, dry_run=dry_run, skip_pytest=skip_pytest)

    pngs = _candidate_pngs(plots_root, output_dirs, repo_root=repo_root)
    missing_csv = _missing_csvs(pngs)
    if missing_csv:
        missing = "\n".join(f"  - {path.as_posix()}" for path in missing_csv)
        raise PlotGalleryBuildError(f"Missing CSV companions for generated PNG plot(s):\n{missing}")

    rendered = _render_missing_static_assets(pngs, dry_run=dry_run, force_render=force_render)
    svg_candidates = tuple(path for path in pngs if not path.with_suffix(".svg").exists())
    svg = ensure_plot_companions.ensure_png_companions(
        list(svg_candidates),
        plots_root,
        dry_run=dry_run,
        create_missing_csv=False,
    )

    _write_manifest(plots_root, dry_run=dry_run)
    index_path = _write_gallery_index(plots_root, dry_run=dry_run)
    report_path = _write_report(plots_root, repo_root=repo_root, dry_run=dry_run)

    return GalleryBuildResult(
        recipe_count=len(selected_recipes),
        plot_targets=plot_targets,
        output_dirs=output_dirs,
        png_count=len(pngs),
        rendered_from_csv=rendered,
        svg_created=svg.svg_created,
        index_path=index_path,
        report_path=report_path,
        dry_run=dry_run,
    )


def _format_result(result: GalleryBuildResult) -> str:
    action = "Would build" if result.dry_run else "Built"
    targets = ", ".join(result.plot_targets)
    folders = ", ".join(f"tests/{folder}" for folder in result.output_dirs)
    return (
        f"{action} static gallery artifacts for {result.recipe_count} recipe(s).\n"
        f"Plot producers: {targets}\n"
        f"Output folders: {folders}\n"
        f"PNG files: {result.png_count}; rendered from CSV: {result.rendered_from_csv}; SVG created: {result.svg_created}\n"
        f"Gallery index: {result.index_path}\n"
        f"Asset report: {result.report_path}"
    )


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Build static PNG/SVG/CSV gallery artifacts for registered test plot recipes.",
    )
    parser.add_argument("targets", nargs="*", help="Source test file, test directory, or pytest node id.")
    parser.add_argument("--all", action="store_true", help="Build every registered test plot recipe.")
    parser.add_argument(
        "--dry-run", action="store_true", help="Report intended work without running pytest or writing files."
    )
    parser.add_argument(
        "--force-render", action="store_true", help="Regenerate PNG and SVG assets from CSV companions."
    )
    parser.add_argument(
        "--skip-pytest", action="store_true", help="Use existing PNG/CSV outputs without running plot producer tests."
    )
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
            force_render=args.force_render,
            skip_pytest=args.skip_pytest,
        )
    except (
        plot_registry.PlotRecipeLookupError,
        PlotGalleryBuildError,
        render_plot_data_csv.PlotDataRenderError,
    ) as exc:
        print(exc, file=sys.stderr)
        return 1
    print(_format_result(result))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
