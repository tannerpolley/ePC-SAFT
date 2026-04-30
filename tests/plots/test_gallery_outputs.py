from __future__ import annotations

import csv
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from scripts import plot_outputs
from scripts.paper_validation.tools import build_analysis_galleries
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


def test_plot_html_path_uses_plot_stem() -> None:
    assert plot_outputs.plot_html_path("docs/plots/tests/api/parity/example.png").as_posix().endswith(
        "docs/plots/tests/api/parity/example.html"
    )


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


def test_plotly_companion_uses_browser_readable_math_labels(tmp_path: Path, monkeypatch) -> None:
    root = tmp_path / "plots"
    monkeypatch.setattr(plot_outputs, "PLOTS_ROOT", root)
    monkeypatch.setattr(plot_outputs, "TEST_PLOTS_ROOT", root / "tests")

    output_path = save_comparison_plot(
        "plotly_math_labels.png",
        "Plotly math label readability regression",
        ["rho", "ares", "lnphi water"],
        np.asarray([1.0, 2.0, 3.0], dtype=float),
        np.asarray([1.0, 2.0, 3.0], dtype=float),
        category=("tests", "readability"),
    )
    html = output_path.with_suffix(".html").read_text(encoding="utf-8")

    assert "$\\rho$" not in html
    assert "\\ln \\phi" not in html
    assert "rho (\\u03c1)" in html
    assert "ln phi (ln \\u03c6)_H2O" in html


def test_plotly_companion_translates_axis_math_labels(tmp_path: Path, monkeypatch) -> None:
    root = tmp_path / "plots"
    monkeypatch.setattr(plot_outputs, "PLOTS_ROOT", root)
    monkeypatch.setattr(plot_outputs, "TEST_PLOTS_ROOT", root / "tests")

    output_path = save_comparison_plot(
        "plotly_axis_math_labels.png",
        "Plotly axis label readability regression",
        ["rho"],
        np.asarray([1.0], dtype=float),
        np.asarray([1.0], dtype=float),
        category=("tests", "readability"),
        ylabel=r"Residual Helmholtz energy, $A^{res}$",
        relative_error=False,
    )
    html = output_path.with_suffix(".html").read_text(encoding="utf-8")

    assert "$A^{res}$" not in html
    assert "Residual Helmholtz energy, A^res" in html


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


def test_root_gallery_embeds_single_page_explorer_manifest(tmp_path: Path, monkeypatch) -> None:
    root = tmp_path / "plots"
    child = root / "paper_validation" / "Example" / "figure_1"
    child.mkdir(parents=True)
    (child / "figure_1.png").write_bytes(b"png")
    (child / "figure_1.svg").write_text("<svg></svg>", encoding="utf-8")
    (child / "figure_1.html").write_text(
        '<html><meta name="epcsaft-interactive-source" content="csv_backfill"></html>',
        encoding="utf-8",
    )
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
    assert '"html_path":"paper_validation/Example/figure_1/figure_1.html"' in html
    assert '"interactive_source":"csv_backfill"' in html
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
    assert "Interactive" in html
    assert "Static" in html
    assert "plotGalleryViewMode" in html
    assert "plotGallerySidebarWidth" in html
    assert 'id="sidebar-resizer"' in html
    assert 'role="separator"' in html
    assert 'aria-label="Resize plot folder sidebar"' in html
    assert "function setSidebarWidth" in html
    assert "startSidebarResize" in html
    assert "finishSidebarResize" in html
    assert "col-resize" in html
    assert "has-hidden-selection" in html
    assert "hidden-selection-pill" in html
    assert "function selectedDescendantCount" in html
    assert "function isDescendantFolder" in html
    assert "hidden inside this collapsed folder" in html
    assert "grid-template-columns: var(--sidebar) var(--sidebar-resizer) minmax(0, 1fr);" in html
    assert ".sidebar-resizer" in html
    assert "sidebar-resizing" in html
    assert 'id="data-modal"' in html
    assert 'id="asset-modal"' in html
    assert 'id="interactive-view"' in html
    assert 'id="static-view"' in html
    assert "interactive-frame" in html
    assert "--interactive-card-height: clamp(330px, 42vw, 560px);" in html
    assert "height: var(--interactive-card-height);" in html
    assert "min-height: var(--interactive-card-height);" in html
    assert 'frame.src = image.html_path;' in html
    assert 'img.src = image.output_path;' in html
    assert "function showDataTable" in html
    assert "function showAssetPreview" in html
    assert "function makeAssetButton" in html
    assert "function makeAssetLink" not in html
    assert 'button.addEventListener("click", () => showAssetPreview(label, path, image));' in html
    assert '.sort((a, b) =>' in html
    assert "parseCsv" in html
    assert 'button.textContent = "Data";' in html
    assert 'badge.textContent = "Static only";' in html
    assert 'badge.textContent = "CSV interactive";' in html
    assert 'badge.textContent = "Interactive";' in html
    assert "CSV-backed interactive reconstruction from plot data" in html
    assert "resource-link is-disabled" not in html
    assert 'target = "_blank"' not in html
    assert 'className = "image-link"' not in html
    assert "Select folders on the left" in html


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
    paths[2].with_suffix(".html").write_text("<html></html>", encoding="utf-8")
    paths[0].with_suffix(".html").write_text(
        '<html><meta name="epcsaft-interactive-source" content="csv_backfill"></html>',
        encoding="utf-8",
    )
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
    assert by_output["tests/equilibrium/vle/equilibrium_vle_compositions.png"]["html_path"] == (
        "tests/equilibrium/vle/equilibrium_vle_compositions.html"
    )
    assert by_output["tests/equilibrium/vle/equilibrium_vle_compositions.png"]["interactive_source"] == (
        "native_plotly"
    )
    assert by_output["paper_validation/Example/figure_1/figure_1.png"]["interactive_source"] == "csv_backfill"
    assert by_output["fits/miac/water/miac.png"]["interactive_source"] == "static_only"
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


def test_docs_test_plot_pngs_have_interactive_html_companions() -> None:
    plots_root = Path("docs/plots/tests")
    missing = []
    for png_path in plots_root.rglob("*.png"):
        if "__pycache__" in png_path.parts:
            continue
        html_path = png_path.with_suffix(".html")
        if not html_path.exists():
            missing.append(png_path.as_posix())
    assert missing == []


def test_plot_gallery_server_skips_busy_preferred_port() -> None:
    import socket

    with socket.socket(socket.AF_INET, socket.SOCK_STREAM) as occupied:
        occupied.bind(("127.0.0.1", 0))
        preferred_port = occupied.getsockname()[1]
        available_port = serve_plot_gallery.find_available_port("127.0.0.1", preferred_port, attempts=3)

    assert available_port != preferred_port
    assert preferred_port < available_port <= preferred_port + 2
