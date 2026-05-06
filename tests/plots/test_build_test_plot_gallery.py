from __future__ import annotations

import csv
from pathlib import Path

import pytest

from scripts import build_test_plot_gallery
from tests.plots import plot_registry

FIELDNAMES = [
    "figure_file",
    "axes_index",
    "axes_title",
    "x_label",
    "y_label",
    "artist_type",
    "artist_label",
    "color",
    "linestyle",
    "marker",
    "linewidth",
    "series_index",
    "point_index",
    "x",
    "y",
    "width",
    "height",
]


def _write_png(path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_bytes(
        b"\x89PNG\r\n\x1a\n"
        b"\x00\x00\x00\rIHDR"
        b"\x00\x00\x00\x10\x00\x00\x00\x08"
        b"\x08\x02\x00\x00\x00"
        b"\x00\x00\x00\x00"
    )


def _write_plot_csv(png_path: Path, rows: list[dict[str, object]]) -> Path:
    csv_path = png_path.parent / f"{png_path.stem}_plot_data.csv"
    csv_path.parent.mkdir(parents=True, exist_ok=True)
    with csv_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=FIELDNAMES, lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow({field: row.get(field, "") for field in FIELDNAMES})
    return csv_path


def _recipe(output_dir: str) -> plot_registry.PlotRecipe:
    return plot_registry.PlotRecipe(
        source_tests=("tests/example/test_feature.py",),
        plot_targets=("tests/plots/test_example_plot_outputs.py",),
        output_dirs=(output_dir,),
    )


def test_registry_resolves_exact_test_file_to_plot_producer() -> None:
    recipes = plot_registry.resolve_plot_recipes(["tests/equilibrium/test_electrolyte_lle.py"])

    targets = {target for recipe in recipes for target in recipe.plot_targets}
    output_dirs = {output_dir for recipe in recipes for output_dir in recipe.output_dirs}
    assert "tests/plots/test_equilibrium_plot_outputs.py" in targets
    assert "equilibrium/electrolyte_lle" in output_dirs


def test_registry_directory_target_resolves_registered_children() -> None:
    recipes = plot_registry.resolve_plot_recipes(["tests/equilibrium"])

    sources = {source for recipe in recipes for source in recipe.source_tests}
    assert "tests/equilibrium/test_electrolyte_lle.py" in sources
    assert "tests/equilibrium/test_vle.py" in sources


def test_registry_unregistered_test_reports_actionable_recipe_pattern() -> None:
    with pytest.raises(plot_registry.PlotRecipeLookupError) as excinfo:
        plot_registry.resolve_plot_recipes(["tests/does_not_exist/test_new_feature.py"])

    message = str(excinfo.value)
    assert "No plot recipe is registered" in message
    assert "PlotRecipe(" in message
    assert "tests/does_not_exist/test_new_feature.py" in message


def test_build_gallery_from_numeric_csv_creates_static_svg_and_report(tmp_path: Path) -> None:
    png_path = tmp_path / "tests" / "plots" / "out" / "example" / "numeric" / "numeric_case.png"
    _write_png(png_path)
    _write_plot_csv(
        png_path,
        [
            {
                "figure_file": png_path.name,
                "axes_index": 0,
                "axes_title": "Numeric case",
                "x_label": "x",
                "y_label": "value",
                "artist_type": "line",
                "series_index": 0,
                "point_index": 0,
                "x": 0.0,
                "y": 1.0,
            },
            {
                "figure_file": png_path.name,
                "axes_index": 0,
                "axes_title": "Numeric case",
                "x_label": "x",
                "y_label": "value",
                "artist_type": "line",
                "series_index": 0,
                "point_index": 1,
                "x": 1.0,
                "y": 2.0,
            },
        ],
    )

    result = build_test_plot_gallery.build_gallery(
        [_recipe("example/numeric")],
        repo_root=tmp_path,
        skip_pytest=True,
        force_render=True,
    )

    assert result.png_count == 1
    assert result.rendered_from_csv == 1
    assert not png_path.with_suffix(".html").exists()
    assert png_path.with_suffix(".svg").exists()
    assert not any(path.name == "manifest" + ".json" for path in tmp_path.rglob("*"))
    assert not any(path.suffix == ".html" for path in tmp_path.rglob("*"))
    assert (tmp_path / "build" / "plot_assets" / "plot_asset_report.csv").exists()


def test_build_gallery_rejects_non_numeric_csv_without_static_html(tmp_path: Path) -> None:
    png_path = tmp_path / "tests" / "plots" / "out" / "example" / "static" / "static_case.png"
    _write_png(png_path)
    _write_plot_csv(png_path, [{"figure_file": png_path.name, "artist_type": "no_numeric_artists"}])

    with pytest.raises(Exception, match="no numeric artists"):
        build_test_plot_gallery.build_gallery(
            [_recipe("example/static")],
            repo_root=tmp_path,
            skip_pytest=True,
            force_render=True,
        )

    assert not png_path.with_suffix(".html").exists()


def test_build_gallery_dry_run_reports_work_without_writing_outputs(tmp_path: Path) -> None:
    png_path = tmp_path / "tests" / "plots" / "out" / "example" / "dry_run" / "dry_run_case.png"
    _write_png(png_path)
    _write_plot_csv(
        png_path,
        [
            {
                "figure_file": png_path.name,
                "artist_type": "line",
                "series_index": 0,
                "point_index": 0,
                "x": 0.0,
                "y": 1.0,
            }
        ],
    )

    result = build_test_plot_gallery.build_gallery(
        [_recipe("example/dry_run")],
        repo_root=tmp_path,
        dry_run=True,
        skip_pytest=True,
        force_render=True,
    )

    assert result.png_count == 1
    assert result.rendered_from_csv == 1
    assert not png_path.with_suffix(".html").exists()
    assert not png_path.with_suffix(".svg").exists()
    assert not any(path.name == "manifest" + ".json" for path in tmp_path.rglob("*"))
    assert not (tmp_path / "build" / "plot_assets" / "plot_asset_report.csv").exists()


def test_build_gallery_requires_csv_companion_for_each_png(tmp_path: Path) -> None:
    png_path = tmp_path / "tests" / "plots" / "out" / "example" / "missing_csv" / "missing_csv_case.png"
    _write_png(png_path)

    with pytest.raises(build_test_plot_gallery.PlotGalleryBuildError) as excinfo:
        build_test_plot_gallery.build_gallery(
            [_recipe("example/missing_csv")],
            repo_root=tmp_path,
            skip_pytest=True,
        )

    assert "Missing CSV companions" in str(excinfo.value)
    assert png_path.as_posix() in str(excinfo.value)
