from __future__ import annotations

from pathlib import Path

import numpy as np

import epcsaft
from epcsaft import _core
from epcsaft import ePCSAFTMixture


REPO_ROOT = Path(__file__).resolve().parents[2]


def _hydrocarbon_mixture() -> ePCSAFTMixture:
    params = {
        "m": np.asarray([1.0, 1.6069, 2.0020]),
        "s": np.asarray([3.7039, 3.5206, 3.6184]),
        "e": np.asarray([150.03, 191.42, 208.11]),
        "k_ij": np.asarray(
            [
                [0.0, 3.0e-4, 1.15e-2],
                [3.0e-4, 0.0, 5.10e-3],
                [1.15e-2, 5.10e-3, 0.0],
            ]
        ),
    }
    return ePCSAFTMixture.from_params(params, species=["Methane", "Ethane", "Propane"])


def test_native_equilibrium_entrypoint_is_exposed() -> None:
    assert hasattr(_core, "_solve_equilibrium_native")


def test_public_equilibrium_result_comes_from_native_backend() -> None:
    mix = _hydrocarbon_mixture()

    result = mix.equilibrium(kind="tp_flash", T=220.0, P=1.0e5, z=[0.1, 0.3, 0.6], backend="native")

    assert isinstance(result, epcsaft.EquilibriumResult)
    assert result.backend == "neutral_vle"
    assert result.diagnostics["solver_language"] == "c++"
    assert result.diagnostics["native_entrypoint"] == "_solve_equilibrium_native"


def test_equilibrium_runtime_does_not_import_scipy_optimizers() -> None:
    source = (REPO_ROOT / "src" / "epcsaft" / "equilibrium.py").read_text(encoding="utf-8")

    forbidden = ("scipy.optimize", "least_squares", "differential_evolution", "minimize_scalar")

    for token in forbidden:
        assert token not in source
