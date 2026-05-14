from __future__ import annotations

from tests.helpers.runtime_cases import _ionic_state
from tests.helpers.runtime_cases import _neutral_state
from analyses.package_validation.package_plot_smokes.tests.plots.plot_helpers import append_payload_rows
from analyses.package_validation.package_plot_smokes.tests.plots.plot_helpers import save_contribution_closure_plot
from analyses.package_validation.package_plot_smokes.tests.plots.plot_helpers import save_contribution_term_breakdown_plot


def _contribution_rows(state) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    append_payload_rows(rows, "ares", state.ares(return_contribution_terms=True))
    append_payload_rows(rows, "Z", state.z(return_contribution_terms=True))
    append_payload_rows(rows, "dadt", state.dadt(return_contribution_terms=True))
    append_payload_rows(rows, "mures", state.mures(return_contribution_terms=True))
    append_payload_rows(rows, "lnphi", state.fugcoef(return_contribution_terms=True))
    return rows


def test_runtime_neutral_contribution_closure_plot() -> None:
    state, _species = _neutral_state()
    save_contribution_closure_plot(
        "runtime_neutral_contribution_closure.png",
        "Neutral contribution closure: term sums vs reported totals",
        _contribution_rows(state),
        category=("contributions", "neutral"),
    )


def test_runtime_neutral_contribution_term_breakdown_plot() -> None:
    state, _species = _neutral_state()
    save_contribution_term_breakdown_plot(
        "runtime_neutral_contribution_term_breakdown.png",
        "Neutral contribution term breakdown by reported quantity",
        _contribution_rows(state),
        category=("contributions", "neutral"),
    )


def test_runtime_ionic_contribution_closure_plot() -> None:
    state, _species = _ionic_state()
    save_contribution_closure_plot(
        "runtime_ionic_contribution_closure.png",
        "Ionic contribution closure: term sums vs reported totals",
        _contribution_rows(state),
        category=("contributions", "ionic"),
    )


def test_runtime_ionic_contribution_term_breakdown_plot() -> None:
    state, _species = _ionic_state()
    save_contribution_term_breakdown_plot(
        "runtime_ionic_contribution_term_breakdown.png",
        "Ionic contribution term breakdown by reported quantity",
        _contribution_rows(state),
        category=("contributions", "ionic"),
    )
