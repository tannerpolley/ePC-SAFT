from __future__ import annotations

from pathlib import Path

from epcsaft.equilibrium_core.confidence import run_confidence_suite
from analyses.package_plot_smokes.tests.plots.plot_helpers import assert_plot_with_data


def test_electrolyte_lle_confidence_plots_are_written_to_local_output(tmp_path: Path) -> None:
    report = run_confidence_suite("khudaida_2026", mode="smoke", output_root=tmp_path, write_gallery=True)

    assert report.output_dir == tmp_path / "khudaida_2026"
    for path in (
        report.residual_gate_plot,
        report.error_plot,
        report.all_tielines_plot,
        report.continuation_plot,
        report.sensitivity_plot,
    ):
        assert path.parent == report.output_dir
        assert_plot_with_data(path)
        assert not path.with_suffix(".html").exists()

    assert not (Path("docs") / "plots" / "tests").exists()
