from __future__ import annotations

import csv
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

from scripts import plot_outputs
from scripts.paper_validation.tools import build_analysis_galleries
from scripts.paper_validation.tools import serve_plot_gallery


def test_save_plot_figure_writes_csv_backing_data(tmp_path: Path) -> None:
    output_path = tmp_path / "figure.png"
    fig, ax = plt.subplots()
    ax.plot([0.0, 1.0], [2.0, 3.0], label="line")
    ax.scatter([0.5], [2.5], label="point")

    try:
        plot_outputs.save_plot_figure(fig, output_path, dpi=72)
    finally:
        plt.close(fig)

    csv_path = tmp_path / "data" / "figure_plot_data.csv"
    assert output_path.exists()
    assert csv_path.exists()
    with csv_path.open("r", newline="", encoding="utf-8-sig") as handle:
        rows = list(csv.DictReader(handle))
    assert {row["artist_type"] for row in rows} >= {"line", "scatter"}
    assert any(row["artist_label"] == "line" and row["x"] == "1" and row["y"] == "3" for row in rows)


def test_root_gallery_embeds_single_page_explorer_manifest(tmp_path: Path, monkeypatch) -> None:
    root = tmp_path / "plots"
    child = root / "paper_validation" / "Example" / "figure_1"
    child.mkdir(parents=True)
    (child / "figure_1.png").write_bytes(b"png")
    monkeypatch.setattr(build_analysis_galleries, "PLOTS_ROOT", root)

    html = build_analysis_galleries.render_gallery_page(root, build_analysis_galleries.collect_pngs(root))

    assert 'data-testid="folder-tree"' in html
    assert 'data-testid="gallery-grid"' in html
    assert "paper_validation/index.html" not in html
    assert "paper_validation/Example/figure_1/figure_1.png" in html
    assert "Select folders on the left" in html


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
    assert "paper_validation/Example/figure_1/figure_1.png" in html
    assert "paper_validation/Example/figure_1/diagnostics/diagnostic.png" in html
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


def test_plot_gallery_server_skips_busy_preferred_port() -> None:
    import socket

    with socket.socket(socket.AF_INET, socket.SOCK_STREAM) as occupied:
        occupied.bind(("127.0.0.1", 0))
        preferred_port = occupied.getsockname()[1]
        available_port = serve_plot_gallery.find_available_port("127.0.0.1", preferred_port, attempts=3)

    assert available_port != preferred_port
    assert preferred_port < available_port <= preferred_port + 2
