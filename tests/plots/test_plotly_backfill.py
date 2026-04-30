from __future__ import annotations

import csv
from pathlib import Path

from scripts.paper_validation.tools import backfill_plotly_html


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


def _write_plot_csv(png_path: Path, rows: list[dict[str, object]]) -> Path:
    png_path.parent.mkdir(parents=True, exist_ok=True)
    png_path.write_bytes(b"png")
    csv_path = png_path.parent / "data" / f"{png_path.stem}_plot_data.csv"
    csv_path.parent.mkdir(parents=True, exist_ok=True)
    with csv_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=FIELDNAMES, lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow({field: row.get(field, "") for field in FIELDNAMES})
    return csv_path


def test_backfill_numeric_csv_generates_marked_plotly_html(tmp_path: Path) -> None:
    png_path = tmp_path / "docs" / "plots" / "fits" / "case" / "figure.png"
    _write_plot_csv(
        png_path,
        [
            {
                "figure_file": "figure.png",
                "axes_index": 0,
                "axes_title": "Line axis",
                "x_label": "mole fraction",
                "y_label": "property value",
                "artist_type": "line",
                "artist_label": "line series",
                "color": "#123456",
                "linestyle": "--",
                "marker": "o",
                "linewidth": 2.5,
                "series_index": 0,
                "point_index": 0,
                "x": 0.0,
                "y": 1.0,
            },
            {
                "figure_file": "figure.png",
                "axes_index": 0,
                "axes_title": "Line axis",
                "x_label": "mole fraction",
                "y_label": "property value",
                "artist_type": "line",
                "artist_label": "line series",
                "color": "#123456",
                "linestyle": "--",
                "marker": "o",
                "linewidth": 2.5,
                "series_index": 0,
                "point_index": 1,
                "x": 1.0,
                "y": 2.0,
            },
            {
                "figure_file": "figure.png",
                "axes_index": 1,
                "axes_title": "Bar axis",
                "x_label": "case",
                "y_label": "bar value",
                "artist_type": "bar",
                "artist_label": "bar series",
                "color": "#abcdef",
                "series_index": 1,
                "point_index": 0,
                "x": 2.0,
                "y": 0.0,
                "width": 0.8,
                "height": 3.0,
            },
            {
                "figure_file": "figure.png",
                "axes_index": 1,
                "axes_title": "Bar axis",
                "artist_type": "scatter",
                "artist_label": "points",
                "series_index": 2,
                "point_index": 0,
                "x": 3.0,
                "y": 4.0,
            },
        ],
    )

    result = backfill_plotly_html.backfill_plotly_html(tmp_path / "docs" / "plots", roots=("fits",))

    html_path = png_path.with_suffix(".html")
    html = html_path.read_text(encoding="utf-8")
    assert result.created == 1
    assert result.skipped == {}
    assert html_path.exists()
    assert "Plotly.newPlot" in html
    assert backfill_plotly_html.BACKFILL_MARKER in html
    assert "line series" in html
    assert "bar series" in html
    assert "points" in html
    assert "#123456" in html
    assert "#abcdef" in html
    assert "dash" in html
    assert "mole fraction" in html
    assert "property value" in html
    assert "Trace" not in html


def test_backfill_infers_readable_labels_when_csv_metadata_is_sparse(tmp_path: Path) -> None:
    png_path = tmp_path / "docs" / "plots" / "fits" / "miac" / "case" / "miac_m_ethanol_NaCl.png"
    csv_path = _write_plot_csv(
        png_path,
        [
            {
                "figure_file": "miac_m_ethanol_NaCl.png",
                "axes_index": 0,
                "axes_title": "Sodium chloride in ethanol",
                "artist_type": "line",
                "series_index": 0,
                "point_index": 0,
                "x": 0.0,
                "y": 1.0,
            },
            {
                "figure_file": "miac_m_ethanol_NaCl.png",
                "axes_index": 0,
                "axes_title": "Sodium chloride in ethanol",
                "artist_type": "line",
                "series_index": 0,
                "point_index": 1,
                "x": 1.0,
                "y": 0.8,
            },
            {
                "figure_file": "miac_m_ethanol_NaCl.png",
                "axes_index": 0,
                "axes_title": "Sodium chloride in ethanol",
                "artist_type": "scatter",
                "series_index": 0,
                "point_index": 0,
                "x": 0.5,
                "y": 0.9,
            },
        ],
    )

    result = backfill_plotly_html.backfill_plotly_html(tmp_path / "docs" / "plots", roots=("fits",))

    fig = backfill_plotly_html.figure_from_plot_csv(csv_path, png_path)
    html = png_path.with_suffix(".html").read_text(encoding="utf-8")
    assert result.created == 1
    assert fig is not None
    assert fig.layout.xaxis.title.text == "m / mol kg⁻¹"
    assert fig.layout.yaxis.title.text == "γ±"
    assert "\\u03b3\\u00b1" in html
    assert "gamma" not in html
    assert "NaCl fit" in html
    assert "NaCl data" in html
    assert "Curve 1" not in html
    assert "Data points 1" not in html
    assert "line 0" not in html
    assert "scatter 0" not in html
    assert fig.data[0].line.color == fig.data[1].marker.color


def test_backfill_infers_presentation_series_names_and_pairs_colors(tmp_path: Path) -> None:
    png_path = tmp_path / "docs" / "plots" / "fits" / "miac" / "case" / "miac_m_2025_present_ethanol_sodium_salts.png"
    csv_path = _write_plot_csv(
        png_path,
        [
            {
                "figure_file": png_path.name,
                "axes_index": 0,
                "axes_title": "Sodium salts in ethanol",
                "artist_type": "line",
                "series_index": 1,
                "point_index": 0,
                "x": 0.0,
                "y": 1.0,
            },
            {
                "figure_file": png_path.name,
                "axes_index": 0,
                "axes_title": "Sodium salts in ethanol",
                "artist_type": "line",
                "series_index": 1,
                "point_index": 1,
                "x": 1.0,
                "y": 0.8,
            },
            {
                "figure_file": png_path.name,
                "axes_index": 0,
                "axes_title": "Sodium salts in ethanol",
                "artist_type": "scatter",
                "series_index": 1,
                "point_index": 0,
                "x": 0.5,
                "y": 0.9,
            },
        ],
    )

    fig = backfill_plotly_html.figure_from_plot_csv(csv_path, png_path)

    assert fig is not None
    assert [trace.name for trace in fig.data] == ["NaBr fit", "NaBr data"]
    assert fig.data[0].line.color == fig.data[1].marker.color


def test_backfill_groups_sparse_grouped_bars_into_readable_series(tmp_path: Path) -> None:
    png_path = tmp_path / "docs" / "plots" / "paper_validation" / "2025_Figiel" / "figure_4" / "figure_4.png"
    rows: list[dict[str, object]] = []
    for category_index in range(3):
        for series_index, offset in enumerate((-0.22, 0.0, 0.22)):
            rows.append(
                {
                    "figure_file": png_path.name,
                    "axes_index": 0,
                    "axes_title": "Gibbs energy of solvation at infinite dilution in water",
                    "artist_type": "bar",
                    "series_index": category_index * 3 + series_index,
                    "point_index": 0,
                    "x": category_index + offset,
                    "width": 0.22,
                    "height": 100.0 + 10.0 * category_index + series_index,
                }
            )
    csv_path = _write_plot_csv(png_path, rows)

    fig = backfill_plotly_html.figure_from_plot_csv(csv_path, png_path)

    assert fig is not None
    assert [trace.name for trace in fig.data] == ["Literature", "ePC-SAFT 2025", "ePC-SAFT 2020"]
    assert all("Bar series" not in str(trace.name) for trace in fig.data)
    assert list(fig.data[0].x) == ["Li+", "Na+", "K+"]


def test_backfill_skips_placeholders_and_dry_run_does_not_write(tmp_path: Path) -> None:
    plots_root = tmp_path / "docs" / "plots"
    placeholder_png = plots_root / "paper_validation" / "placeholder.png"
    dry_run_png = plots_root / "fits" / "dry_run.png"
    _write_plot_csv(
        placeholder_png,
        [{"figure_file": "placeholder.png", "artist_type": "existing_png_backfill"}],
    )
    _write_plot_csv(
        dry_run_png,
        [{"figure_file": "dry_run.png", "axes_index": 0, "artist_type": "scatter", "x": 1.0, "y": 2.0}],
    )

    dry_run = backfill_plotly_html.backfill_plotly_html(plots_root, dry_run=True)
    assert dry_run.created == 1
    assert not dry_run_png.with_suffix(".html").exists()

    result = backfill_plotly_html.backfill_plotly_html(plots_root)

    assert result.created == 1
    assert result.skipped == {"no_numeric_artists": 1}
    assert dry_run_png.with_suffix(".html").exists()
    assert not placeholder_png.with_suffix(".html").exists()


def test_backfill_preserves_existing_html_unless_force_is_used(tmp_path: Path) -> None:
    plots_root = tmp_path / "docs" / "plots"
    png_path = plots_root / "fits" / "existing.png"
    html_path = png_path.with_suffix(".html")
    _write_plot_csv(
        png_path,
        [{"figure_file": "existing.png", "axes_index": 0, "artist_type": "scatter", "x": 1.0, "y": 2.0}],
    )
    html_path.write_text("existing html", encoding="utf-8")

    skipped = backfill_plotly_html.backfill_plotly_html(plots_root)
    assert skipped.created == 0
    assert skipped.skipped == {"existing_html": 1}
    assert html_path.read_text(encoding="utf-8") == "existing html"

    forced = backfill_plotly_html.backfill_plotly_html(plots_root, force=True)
    assert forced.created == 1
    assert backfill_plotly_html.BACKFILL_MARKER in html_path.read_text(encoding="utf-8")
