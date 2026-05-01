from __future__ import annotations

from pathlib import Path

from epcsaft.equilibrium_core.confidence import run_confidence_suite
from tests.plots.plot_helpers import assert_plot_with_data


def test_electrolyte_lle_confidence_plots_are_written_to_gallery(tmp_path: Path) -> None:
    report = run_confidence_suite("khudaida_2026", mode="full", output_root=tmp_path, write_gallery=True)

    for path in (
        report.residual_gate_plot,
        report.error_plot,
        report.continuation_plot,
        report.sensitivity_plot,
    ):
        gallery_path = Path("docs") / "plots" / "tests" / "equilibrium" / "electrolyte_lle_confidence" / path.name
        assert_plot_with_data(gallery_path)
        assert gallery_path.with_suffix(".html").is_file()

