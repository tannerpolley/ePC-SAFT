from __future__ import annotations

import csv
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pytest

from scripts import plot_outputs
from scripts.paper_validation.tools import build_analysis_galleries
from scripts.paper_validation.tools import ensure_plot_companions
from scripts.paper_validation.tools import render_plot_data_csv
from scripts.paper_validation.tools import report_plot_assets
from scripts.paper_validation.tools import serve_plot_gallery
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

    csv_path = tmp_path / "data" / "figure_plot_data.csv"
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
    monkeypatch.setattr(plot_outputs, "TEST_PLOTS_ROOT", root / "tests")

    output_path = plot_outputs.test_plot_path("tests/plots/test_reference_comparison_outputs.py", "x.png")
    semantic_path = plot_outputs.test_plot_path(
        "tests/plots/test_property_plot_outputs.py",
        "residual.png",
        category=("properties", "residual_energy"),
    )

    assert output_path == root / "tests" / "plots" / "reference_comparison_outputs" / "x.png"
    assert semantic_path == root / "tests" / "properties" / "residual_energy" / "residual.png"
    assert output_path.parent.is_dir()
    assert semantic_path.parent.is_dir()


def test_readable_plot_helper_wraps_dense_labels_without_clipping(tmp_path: Path, monkeypatch) -> None:
    root = tmp_path / "plots"
    monkeypatch.setattr(plot_outputs, "PLOTS_ROOT", root)
    monkeypatch.setattr(plot_outputs, "TEST_PLOTS_ROOT", root / "tests")

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
    csv_path = tmp_path / "data" / "example_plot_data.csv"
    csv_path.parent.mkdir()
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
    csv_path = tmp_path / "data" / "placeholder_plot_data.csv"
    csv_path.parent.mkdir()
    csv_path.write_text("figure_file,artist_type\nplaceholder.png,no_numeric_artists\n", encoding="utf-8")

    with pytest.raises(render_plot_data_csv.PlotDataRenderError, match="no numeric artists"):
        render_plot_data_csv.render_csv_to_static_assets(csv_path)


def test_root_gallery_embeds_static_single_page_manifest(tmp_path: Path, monkeypatch) -> None:
    root = tmp_path / "plots"
    child = root / "paper_validation" / "Example" / "figure_1"
    child.mkdir(parents=True)
    (child / "figure_1.png").write_bytes(b"png")
    (child / "figure_1.svg").write_text("<svg></svg>", encoding="utf-8")
    (child / "data").mkdir()
    (child / "data" / "figure_1_plot_data.csv").write_text("x,y\n0,0\n", encoding="utf-8")
    for test_child in (
        root / "tests" / "equilibrium" / "vle",
        root / "tests" / "properties" / "residual_energy",
        root / "tests" / "properties" / "activity_fugacity",
        root / "tests" / "contributions" / "neutral",
        root / "tests" / "regression" / "hydrocarbon",
        root / "tests" / "native" / "derivatives",
    ):
        test_child.mkdir(parents=True)
        (test_child / f"{test_child.name}.png").write_bytes(b"png")
    monkeypatch.setattr(build_analysis_galleries, "PLOTS_ROOT", root)

    html = build_analysis_galleries.render_gallery_page(root, build_analysis_galleries.collect_pngs(root))

    assert 'data-testid="folder-tree"' in html
    assert 'data-testid="gallery-grid"' in html
    assert "paper_validation/index.html" not in html
    assert '"output_path":"paper_validation/Example/figure_1/figure_1.png"' in html
    assert '"svg_path":"paper_validation/Example/figure_1/figure_1.svg"' in html
    assert '"data_path":"paper_validation/Example/figure_1/data/figure_1_plot_data.csv"' in html
    assert '"source_path":"scripts/paper_validation/Example_analysis/figure_1/figure_1.png"' in html
    assert '"source_path":"tests/equilibrium/vle/vle.png"' in html
    assert '"source_folder":"tests/properties/residual_energy"' in html
    assert '"source_folder":"tests/properties/activity_fugacity"' in html
    assert '"source_folder":"tests/contributions/neutral"' in html
    assert '"source_folder":"tests/regression/hydrocarbon"' in html
    assert '"source_folder":"tests/native/derivatives"' in html
    assert "tests/plots/reference_comparison_outputs" not in html
    assert "Source tree" in html
    assert "Output tree" in html
    assert "plotGallerySidebarWidth" in html
    assert "plotGalleryStaticStateV1" in html
    assert "function saveGalleryState" in html
    assert "function restoreGalleryState" in html
    assert "function expandAncestors" in html
    assert "localStorage.setItem(galleryStateStorageKey" in html
    assert 'id="sidebar-resizer"' in html
    assert 'role="separator"' in html
    assert 'aria-label="Resize plot folder sidebar"' in html
    assert "function setSidebarWidth" in html
    assert "startSidebarResize" in html
    assert "finishSidebarResize" in html
    assert "col-resize" in html
    assert 'id="data-modal"' in html
    assert 'id="asset-modal"' in html
    assert 'img.src = image.output_path;' in html
    assert "function showDataTable" in html
    assert "function showAssetPreview" in html
    assert "function makeAssetButton" in html
    assert "parseCsv" in html
    assert 'button.textContent = "Data";' in html
    assert "Select folders on the left" in html
    assert ("html" + "_path") not in html
    assert ("native_" + "plot" + "ly") not in html
    assert ("csv_" + "backfill") not in html
    assert ("static_png_" + "wrapper") not in html


def test_ensure_plot_companions_adds_svg_and_csv(tmp_path: Path, monkeypatch) -> None:
    root = tmp_path / "plots"
    png = root / "paper_validation" / "Example" / "figure_1" / "figure_1.png"
    png.parent.mkdir(parents=True)
    png.write_bytes(
        b"\x89PNG\r\n\x1a\n"
        b"\x00\x00\x00\rIHDR"
        b"\x00\x00\x00\x10\x00\x00\x00\x08"
        b"\x08\x02\x00\x00\x00"
        b"\x00\x00\x00\x00"
    )
    monkeypatch.setattr(build_analysis_galleries, "PLOTS_ROOT", root)

    result = ensure_plot_companions.ensure_companions(root)

    assert result.csv_created == 1
    assert result.svg_created == 1
    assert (png.parent / "data" / "figure_1_plot_data.csv").exists()
    assert not png.with_suffix(".html").exists()
    assert 'href="figure_1.png"' in png.with_suffix(".svg").read_text(encoding="utf-8")


def test_plot_asset_report_lists_static_png_svg_csv_bundles(tmp_path: Path, monkeypatch) -> None:
    root = tmp_path / "plots"
    complete = root / "paper_validation" / "Example" / "figure_1" / "figure_1.png"
    incomplete = root / "paper_validation" / "Example" / "figure_2" / "figure_2.png"
    complete.parent.mkdir(parents=True)
    incomplete.parent.mkdir(parents=True)
    complete.write_bytes(b"png")
    incomplete.write_bytes(b"png")
    complete.with_suffix(".svg").write_text("<svg></svg>", encoding="utf-8")
    (complete.parent / "data").mkdir()
    (complete.parent / "data" / "figure_1_plot_data.csv").write_text("x,y\n0,0\n", encoding="utf-8")
    monkeypatch.setattr(build_analysis_galleries, "PLOTS_ROOT", root)

    output = tmp_path / "report.csv"
    rows = report_plot_assets.asset_rows(root)
    report_plot_assets.write_report(output, root)
    report = output.read_text(encoding="utf-8")

    by_output = {row["output_path"]: row for row in rows}
    assert by_output["paper_validation/Example/figure_1/figure_1.png"]["bundle_complete"] == "true"
    assert by_output["paper_validation/Example/figure_2/figure_2.png"]["bundle_complete"] == "false"
    assert "bundle_complete" in report
    assert "has_svg" in report
    assert "has_csv" in report
    assert "has_html" not in report


def test_gallery_manifest_keeps_output_and_source_paths(tmp_path: Path, monkeypatch) -> None:
    root = tmp_path / "plots"
    paths = [
        root / "paper_validation" / "Example" / "figure_1" / "figure_1.png",
        root / "fits" / "miac" / "water" / "miac.png",
        root / "tests" / "equilibrium" / "vle" / "equilibrium_vle_compositions.png",
    ]
    for path in paths:
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_bytes(b"png")
    paths[2].with_suffix(".svg").write_text("<svg></svg>", encoding="utf-8")
    data_path = paths[2].parent / "data" / "equilibrium_vle_compositions_plot_data.csv"
    data_path.parent.mkdir(parents=True, exist_ok=True)
    data_path.write_text("x,y\n0,0\n", encoding="utf-8")
    monkeypatch.setattr(build_analysis_galleries, "PLOTS_ROOT", root)

    manifest = build_analysis_galleries.image_manifest(build_analysis_galleries.collect_pngs(root))

    by_output = {item["output_path"]: item for item in manifest}
    assert by_output["paper_validation/Example/figure_1/figure_1.png"]["source_path"] == (
        "scripts/paper_validation/Example_analysis/figure_1/figure_1.png"
    )
    assert by_output["fits/miac/water/miac.png"]["source_path"] == "scripts/fits/miac/water/miac.png"
    assert by_output["tests/equilibrium/vle/equilibrium_vle_compositions.png"]["source_path"] == (
        "tests/equilibrium/vle/equilibrium_vle_compositions.png"
    )
    assert by_output["tests/equilibrium/vle/equilibrium_vle_compositions.png"]["svg_path"] == (
        "tests/equilibrium/vle/equilibrium_vle_compositions.svg"
    )
    assert ("html" + "_path") not in by_output["tests/equilibrium/vle/equilibrium_vle_compositions.png"]
    assert by_output["tests/equilibrium/vle/equilibrium_vle_compositions.png"]["data_path"] == (
        "tests/equilibrium/vle/data/equilibrium_vle_compositions_plot_data.csv"
    )
    assert by_output["paper_validation/Example/figure_1/figure_1.png"]["folder"] == (
        "scripts/paper_validation/Example_analysis/figure_1"
    )
    assert by_output["paper_validation/Example/figure_1/figure_1.png"]["output_folder"] == (
        "paper_validation/Example/figure_1"
    )


def test_plot_gallery_main_writes_only_root_index_and_removes_nested_indexes(tmp_path: Path, monkeypatch) -> None:
    root = tmp_path / "plots"
    leaf = root / "paper_validation" / "Example" / "figure_1"
    nested = leaf / "diagnostics"
    nested.mkdir(parents=True)
    (leaf / "figure_1.png").write_bytes(b"png")
    (nested / "diagnostic.png").write_bytes(b"png")
    (leaf / "index.html").write_text("stale nested page", encoding="utf-8")
    monkeypatch.setattr(build_analysis_galleries, "PLOTS_ROOT", root)

    build_analysis_galleries.main()
    html = (root / "index.html").read_text(encoding="utf-8")

    assert (root / "index.html").exists()
    assert not (leaf / "index.html").exists()
    assert '"output_path":"paper_validation/Example/figure_1/figure_1.png"' in html
    assert '"output_path":"paper_validation/Example/figure_1/diagnostics/diagnostic.png"' in html
    assert "/index.html" not in html


def test_docs_plot_pngs_have_companion_csv_data() -> None:
    plots_root = Path("docs/plots")
    missing = []
    for png_path in plots_root.rglob("*.png"):
        if "__pycache__" in png_path.parts:
            continue
        csv_path = png_path.parent / "data" / f"{png_path.stem}_plot_data.csv"
        if not csv_path.exists():
            missing.append(png_path.as_posix())
    assert missing == []


def test_docs_plot_pngs_have_gallery_manifest_entries() -> None:
    plots_root = build_analysis_galleries.PLOTS_ROOT
    pngs = build_analysis_galleries.collect_pngs(plots_root)
    manifest = build_analysis_galleries.image_manifest(pngs)

    assert len(manifest) == len(pngs)
    for item in manifest:
        png_path = plots_root / item["output_path"]
        assert png_path.exists()
        assert item["path"] == item["output_path"]
        assert item["source_path"]
        assert ("html" + "_path") not in item
        assert ("interactive" + "_source") not in item
        if item["svg_path"]:
            assert (plots_root / item["svg_path"]).exists()
        if item["data_path"]:
            assert (plots_root / item["data_path"]).exists()


def test_docs_test_plot_pngs_have_companion_svg_data() -> None:
    plots_root = Path("docs/plots/tests")
    missing = []
    for png_path in plots_root.rglob("*.png"):
        if "__pycache__" in png_path.parts:
            continue
        svg_path = png_path.with_suffix(".svg")
        if not svg_path.exists():
            missing.append(png_path.as_posix())
    assert missing == []


def test_docs_plots_have_only_root_gallery_html() -> None:
    html_files = sorted(path.as_posix() for path in Path("docs/plots").rglob("*.html"))
    assert html_files == ["docs/plots/index.html"]


def test_plot_gallery_server_skips_busy_preferred_port() -> None:
    import socket

    with socket.socket(socket.AF_INET, socket.SOCK_STREAM) as occupied:
        occupied.bind(("127.0.0.1", 0))
        preferred_port = occupied.getsockname()[1]
        available_port = serve_plot_gallery.find_available_port("127.0.0.1", preferred_port, attempts=3)

    assert available_port != preferred_port
    assert preferred_port < available_port <= preferred_port + 2
