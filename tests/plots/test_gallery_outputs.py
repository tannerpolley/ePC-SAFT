from __future__ import annotations

import csv
import importlib
import json
import shutil
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pytest

from scripts import plot_outputs
from scripts.paper_validation.tools import ensure_plot_companions
from scripts.paper_validation.tools import plot_asset_index
from scripts.paper_validation.tools import render_plot_data_csv
from scripts.paper_validation.tools import report_plot_assets
from tests.plots.plot_helpers import assert_figure_text_is_inside_canvas
from tests.plots.plot_helpers import save_comparison_plot


def test_save_plot_figure_writes_csv_backing_data(tmp_path: Path) -> None:
    output_path = tmp_path / "figure.png"
    fig, ax = plt.subplots()
    ax.set_xlabel("composition")
    ax.set_ylabel("value")
    ax.plot([0.0, 1.0], [2.0, 3.0], color="#123456", linestyle="--", marker="o", label="line")
    ax.scatter([0.5], [2.5], color="#abcdef", label="point")

    try:
        plot_outputs.save_plot_figure(fig, output_path, dpi=72, svg_companion=True)
    finally:
        plt.close(fig)

    csv_path = tmp_path / "figure_plot_data.csv"
    svg_path = tmp_path / "figure.svg"
    assert output_path.exists()
    assert svg_path.exists()
    assert csv_path.exists()
    with csv_path.open("r", newline="", encoding="utf-8-sig") as handle:
        rows = list(csv.DictReader(handle))
    assert {row["artist_type"] for row in rows} >= {"line", "scatter"}
    assert any(row["artist_label"] == "line" and row["x"] == "1" and row["y"] == "3" for row in rows)
    line_row = next(row for row in rows if row["artist_type"] == "line")
    scatter_row = next(row for row in rows if row["artist_type"] == "scatter")
    assert line_row["x_label"] == "composition"
    assert line_row["y_label"] == "value"
    assert line_row["color"] == "#123456"
    assert line_row["linestyle"] == "--"
    assert line_row["marker"] == "o"
    assert scatter_row["color"] == "#abcdef"


def test_test_plot_path_maps_test_file_to_docked_plots_folder(tmp_path: Path, monkeypatch) -> None:
    root = tmp_path / "plots"
    monkeypatch.setattr(plot_outputs, "PLOTS_ROOT", root)
    monkeypatch.setattr(plot_outputs, "TEST_PLOTS_ROOT", root / "tests" / "out")

    output_path = plot_outputs.test_plot_path("tests/plots/test_reference_comparison_outputs.py", "x.png")
    semantic_path = plot_outputs.test_plot_path(
        "tests/plots/test_property_plot_outputs.py",
        "residual.png",
        category=("properties", "residual_energy"),
    )

    assert output_path == root / "tests" / "out" / "plots" / "reference_comparison_outputs" / "x.png"
    assert semantic_path == root / "tests" / "out" / "properties" / "residual_energy" / "residual.png"
    assert output_path.parent.is_dir()
    assert semantic_path.parent.is_dir()


def test_readable_plot_helper_wraps_dense_labels_without_clipping(tmp_path: Path, monkeypatch) -> None:
    root = tmp_path / "plots"
    monkeypatch.setattr(plot_outputs, "PLOTS_ROOT", root)
    monkeypatch.setattr(plot_outputs, "TEST_PLOTS_ROOT", root / "tests" / "out")

    labels = [
        "very long residual Helmholtz label",
        "temperature derivative residual Helmholtz",
        "mean ionic activity coefficient molality",
        "solvation free energy sodium ion",
        "osmotic coefficient solution diagnostic",
        "fugacity coefficient chloride ion",
    ]
    output_path = save_comparison_plot(
        "dense_readability.png",
        "Dense label readability regression",
        labels,
        np.asarray([1.0, -2.0e5, 3.0, -4.0e-4, 5.0, -6.0], dtype=float),
        np.asarray([1.0, -2.0e5, 3.0, -4.0e-4, 5.0, -6.0], dtype=float),
        category=("tests", "readability"),
    )

    assert output_path.exists()
    assert output_path.with_suffix(".svg").exists()
    assert not output_path.with_suffix(".html").exists()


def test_text_canvas_check_detects_clipped_labels() -> None:
    fig, ax = plt.subplots(figsize=(3.0, 2.0))
    ax.text(-0.35, 0.5, "clipped", transform=ax.transAxes)

    try:
        try:
            assert_figure_text_is_inside_canvas(fig)
        except AssertionError:
            pass
        else:
            raise AssertionError("Expected clipped text to fail the plot readability check.")
    finally:
        plt.close(fig)


def test_render_plot_data_csv_writes_png_and_svg(tmp_path: Path) -> None:
    csv_path = tmp_path / "example_plot_data.csv"
    csv_path.write_text(
        "\n".join(
            [
                "figure_file,axes_index,axes_title,x_label,y_label,artist_type,artist_label,color,linestyle,marker,linewidth,series_index,point_index,x,y,width,height",
                "example.png,0,Example,x,y,line,line,#123456,--,o,1,0,0,0,1,,",
                "example.png,0,Example,x,y,line,line,#123456,--,o,1,0,1,1,2,,",
                "example.png,0,Example,x,y,scatter,points,#abcdef,,,1,1,0,0.5,1.5,,",
                "example.png,0,Example,x,y,bar,bars,#ff0000,,,1,2,0,2,0,0.4,3",
            ]
        )
        + "\n",
        encoding="utf-8",
    )

    output = render_plot_data_csv.render_csv_to_static_assets(csv_path)

    assert output == tmp_path / "example.png"
    assert output.exists()
    assert output.with_suffix(".svg").exists()


def test_render_plot_data_csv_rejects_placeholder_rows(tmp_path: Path) -> None:
    csv_path = tmp_path / "placeholder_plot_data.csv"
    csv_path.write_text("figure_file,artist_type\nplaceholder.png,no_numeric_artists\n", encoding="utf-8")

    with pytest.raises(render_plot_data_csv.PlotDataRenderError, match="no numeric artists"):
        render_plot_data_csv.render_csv_to_static_assets(csv_path)


def test_plot_asset_index_describes_source_owned_static_assets(tmp_path: Path) -> None:
    child = tmp_path / "scripts" / "paper_validation" / "Example_analysis" / "figure_1" / "out"
    child.mkdir(parents=True)
    (child / "figure_1.png").write_bytes(b"png")
    (child / "figure_1.svg").write_text("<svg></svg>", encoding="utf-8")
    (child / "figure_1_plot_data.csv").write_text("x,y\n0,0\n", encoding="utf-8")
    for test_child in (
        tmp_path / "tests" / "plots" / "out" / "equilibrium" / "vle",
        tmp_path / "tests" / "plots" / "out" / "properties" / "residual_energy",
        tmp_path / "tests" / "plots" / "out" / "properties" / "activity_fugacity",
        tmp_path / "tests" / "plots" / "out" / "contributions" / "neutral",
        tmp_path / "tests" / "plots" / "out" / "regression" / "hydrocarbon",
        tmp_path / "tests" / "plots" / "out" / "native" / "derivatives",
    ):
        test_child.mkdir(parents=True)
        (test_child / f"{test_child.name}.png").write_bytes(b"png")

    assets = plot_asset_index.asset_rows(plot_asset_index.collect_pngs(tmp_path), repo_root=tmp_path)
    by_output = {item["output_path"]: item for item in assets}

    assert by_output["scripts/paper_validation/Example_analysis/figure_1/out/figure_1.png"]["svg_path"] == (
        "scripts/paper_validation/Example_analysis/figure_1/out/figure_1.svg"
    )
    assert by_output["scripts/paper_validation/Example_analysis/figure_1/out/figure_1.png"]["data_path"] == (
        "scripts/paper_validation/Example_analysis/figure_1/out/figure_1_plot_data.csv"
    )
    assert by_output["scripts/paper_validation/Example_analysis/figure_1/out/figure_1.png"]["source_path"] == (
        "scripts/paper_validation/Example_analysis/figure_1/figure_1.png"
    )
    assert by_output["tests/plots/out/equilibrium/vle/vle.png"]["source_path"] == (
        "tests/plots/equilibrium/vle/vle.png"
    )
    assert by_output["tests/plots/out/properties/residual_energy/residual_energy.png"]["source_folder"] == (
        "tests/plots/properties/residual_energy"
    )
    assert all(("html" + "_path") not in item for item in assets)
    assert all(("native_" + "plot" + "ly") not in json.dumps(item) for item in assets)


def test_khudaida_lle_plots_include_model_paper_experiment_and_feed_series(tmp_path: Path) -> None:
    common = importlib.import_module("scripts.paper_validation.2026_Khudaida_analysis._common")
    source_dir = Path("scripts/paper_validation/2026_Khudaida_analysis/figure_2").resolve()
    figure_dir = tmp_path / "figure_2"
    shutil.copytree(source_dir / "data", figure_dir / "data")
    (figure_dir / "out").mkdir()
    shutil.copy2(figure_dir / "data" / "model_tielines.csv", figure_dir / "out" / "model_tielines.csv")

    common.plot_lle_figure(figure_dir, 2, 293.15, 0.05)

    with (figure_dir / "out" / "figure_2_plot_data.csv").open("r", newline="", encoding="utf-8-sig") as handle:
        full_labels = {row["artist_label"] for row in csv.DictReader(handle)}
    with (figure_dir / "out" / "figure_2_scaled_plot_data.csv").open("r", newline="", encoding="utf-8-sig") as handle:
        scaled_labels = {row["artist_label"] for row in csv.DictReader(handle)}

    assert {"Exp.", "model ePC-SAFT", "paper ePC-SAFT", "Feed"} <= full_labels
    assert {
        "Exp. organic",
        "model ePC-SAFT organic",
        "paper ePC-SAFT organic",
    } <= scaled_labels


def test_ensure_plot_companions_adds_svg_and_csv(tmp_path: Path) -> None:
    png = tmp_path / "scripts" / "paper_validation" / "Example_analysis" / "figure_1" / "out" / "figure_1.png"
    png.parent.mkdir(parents=True)
    png.write_bytes(
        b"\x89PNG\r\n\x1a\n"
        b"\x00\x00\x00\rIHDR"
        b"\x00\x00\x00\x10\x00\x00\x00\x08"
        b"\x08\x02\x00\x00\x00"
        b"\x00\x00\x00\x00"
    )

    result = ensure_plot_companions.ensure_companions(tmp_path)

    assert result.csv_created == 1
    assert result.svg_created == 1
    assert (png.parent / "figure_1_plot_data.csv").exists()
    assert not png.with_suffix(".html").exists()
    assert 'href="figure_1.png"' in png.with_suffix(".svg").read_text(encoding="utf-8")


def test_plot_asset_report_lists_static_png_svg_csv_bundles(tmp_path: Path) -> None:
    complete = tmp_path / "scripts" / "paper_validation" / "Example_analysis" / "figure_1" / "out" / "figure_1.png"
    incomplete = tmp_path / "scripts" / "paper_validation" / "Example_analysis" / "figure_2" / "out" / "figure_2.png"
    complete.parent.mkdir(parents=True)
    incomplete.parent.mkdir(parents=True)
    complete.write_bytes(b"png")
    incomplete.write_bytes(b"png")
    complete.with_suffix(".svg").write_text("<svg></svg>", encoding="utf-8")
    (complete.parent / "figure_1_plot_data.csv").write_text("x,y\n0,0\n", encoding="utf-8")

    output = tmp_path / "report.csv"
    rows = report_plot_assets.asset_rows(tmp_path)
    report_plot_assets.write_report(output, tmp_path)
    report = output.read_text(encoding="utf-8")

    by_output = {row["output_path"]: row for row in rows}
    assert by_output["scripts/paper_validation/Example_analysis/figure_1/out/figure_1.png"]["bundle_complete"] == "true"
    assert (
        by_output["scripts/paper_validation/Example_analysis/figure_2/out/figure_2.png"]["bundle_complete"] == "false"
    )
    assert "bundle_complete" in report
    assert "has_svg" in report
    assert "has_csv" in report
    assert "has_html" not in report


def test_plot_asset_index_keeps_output_and_source_paths(tmp_path: Path) -> None:
    paths = [
        tmp_path / "scripts" / "paper_validation" / "Example_analysis" / "figure_1" / "out" / "figure_1.png",
        tmp_path / "scripts" / "fits" / "out" / "miac" / "water" / "miac.png",
        tmp_path / "tests" / "plots" / "out" / "equilibrium" / "vle" / "equilibrium_vle_compositions.png",
    ]
    for path in paths:
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_bytes(b"png")
    paths[2].with_suffix(".svg").write_text("<svg></svg>", encoding="utf-8")
    data_path = paths[2].parent / "equilibrium_vle_compositions_plot_data.csv"
    data_path.write_text("x,y\n0,0\n", encoding="utf-8")

    assets = plot_asset_index.asset_rows(plot_asset_index.collect_pngs(tmp_path), repo_root=tmp_path)

    by_output = {item["output_path"]: item for item in assets}
    test_output = "tests/plots/out/equilibrium/vle/equilibrium_vle_compositions.png"
    assert by_output["scripts/paper_validation/Example_analysis/figure_1/out/figure_1.png"]["source_path"] == (
        "scripts/paper_validation/Example_analysis/figure_1/figure_1.png"
    )
    assert by_output["scripts/fits/out/miac/water/miac.png"]["source_path"] == "scripts/fits/miac/water/miac.png"
    assert by_output[test_output]["source_path"] == "tests/plots/equilibrium/vle/equilibrium_vle_compositions.png"
    assert by_output[test_output]["svg_path"] == "tests/plots/out/equilibrium/vle/equilibrium_vle_compositions.svg"
    assert ("html" + "_path") not in by_output[test_output]
    assert by_output[test_output]["data_path"] == (
        "tests/plots/out/equilibrium/vle/equilibrium_vle_compositions_plot_data.csv"
    )


def test_docs_plots_has_no_repo_owned_html_or_catalog_app() -> None:
    docs_plots = Path("docs") / "plots"
    html_files = sorted(path.as_posix() for path in docs_plots.rglob("*.html"))
    assert html_files == []
    assert not (docs_plots / ("manifest" + ".json")).exists()
