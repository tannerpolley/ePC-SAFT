"""State-view helpers for the public ePC-SAFT EOS Harness."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Mapping

import numpy as np

from ._types import InputError, SolutionError, vector_to_array

STATE_METHOD_ALIAS_MAP = {
    "pressure": "p",
    "density": "rho",
    "molar_density": "rho_molar",
    "mass_density": "rho_mass",
    "compressibility_factor": "z",
    "residual_helmholtz": "ares",
    "temperature_derivative_residual_helmholtz": "dadt",
    "composition_derivative_residual_helmholtz": "dadx",
    "residual_enthalpy": "hres",
    "residual_entropy": "sres",
    "residual_gibbs": "gres",
    "residual_chemical_potential": "mures",
    "activity_coefficient": "gamma",
    "fugacity_coefficient": "fugcoef",
    "relative_permittivity": "epsr",
    "osmotic_coefficient": "osmotic_coef",
    "state_diagnostics": "diag",
    "solvation_free_energy": "gsolv",
}
STATE_METHOD_ALIAS_LOOKUP = {alias: name for name, alias in STATE_METHOD_ALIAS_MAP.items()}
CONTRIBUTION_NAMES = ("hc", "disp", "assoc", "ion", "born")
CONTRIBUTION_PUBLIC_NAMES = {
    "hc": "hard_chain",
    "disp": "dispersion",
    "assoc": "association",
    "ion": "ionic",
    "born": "born",
}
GAS_CONSTANT = 8.31446261815324


@dataclass(frozen=True, slots=True)
class StateDiagnosticsView:
    """Read-only typed view over serialized EOS state diagnostics."""

    diagnostics: Mapping[str, Any]

    def __post_init__(self) -> None:
        object.__setattr__(self, "diagnostics", dict(self.diagnostics or {}))

    @property
    def temperature(self) -> float:
        """Return state temperature."""
        return float(self.diagnostics.get("T", float("nan")))

    @property
    def phase(self) -> int:
        """Return the native phase token."""
        return int(self.diagnostics.get("phase", -1))

    @property
    def composition(self) -> np.ndarray:
        """Return the state composition vector."""
        return np.asarray(self.diagnostics.get("x", []), dtype=float)

    @property
    def pressure(self) -> float:
        """Return state pressure."""
        return float(self.diagnostics.get("pressure", float("nan")))

    @property
    def density(self) -> float:
        """Return state molar density."""
        return float(self.diagnostics.get("density_molar", self.diagnostics.get("density", float("nan"))))

    @property
    def compressibility_factor(self) -> float:
        """Return the compressibility factor."""
        return float(self.diagnostics.get("compressibility_factor", float("nan")))

    @property
    def has_ionic_outputs(self) -> bool:
        """Return whether ionic-only state outputs are populated."""
        return self.diagnostics.get("relative_permittivity") is not None

    @property
    def fugacity_coefficient_terms(self) -> Mapping[str, Any]:
        """Return the serialized fugacity-coefficient contribution terms."""
        return dict(self.diagnostics.get("fugacity_coefficient_terms", {}))

    def as_dict(self) -> dict[str, Any]:
        """Return a shallow copy of the stable serialized diagnostics payload."""
        return dict(self.diagnostics)


def sum_vector_terms(terms: Mapping[str, Any]) -> np.ndarray:
    """Return the sum over contribution-family vector terms."""
    total = None
    for name in CONTRIBUTION_NAMES:
        value = np.asarray(terms[name], dtype=float)
        total = value.copy() if total is None else total + value
    return np.asarray(total, dtype=float)


def public_contribution_terms(terms: Mapping[str, Any]) -> dict[str, Any]:
    """Map internal contribution-family names to public names."""
    out = {}
    for internal, public in CONTRIBUTION_PUBLIC_NAMES.items():
        out[public] = terms[internal]
    if "ideal" in terms:
        out["ideal"] = terms["ideal"]
    return out


def scalar_terms_dict(terms: Any) -> dict[str, float]:
    """Convert native scalar contribution terms to a dict."""
    return {
        "hc": float(terms.hc),
        "disp": float(terms.disp),
        "assoc": float(terms.assoc),
        "ion": float(terms.ion),
        "born": float(terms.born),
        "total": float(terms.total),
    }


def vector_terms_dict(terms: Any, expected_size: int, label: str) -> dict[str, np.ndarray]:
    """Convert native vector contribution terms to arrays with shape validation."""
    out = {}
    blocks = {
        "hc": vector_to_array(terms.hc),
        "disp": vector_to_array(terms.disp),
        "assoc": vector_to_array(terms.assoc),
        "ion": vector_to_array(terms.ion),
        "born": vector_to_array(terms.born),
        "total": vector_to_array(terms.total),
    }
    for name, arr in blocks.items():
        arr = np.asarray(arr, dtype=float)
        if arr.size != expected_size:
            raise SolutionError(
                f"Unexpected {label} payload size for {name}: expected {int(expected_size)}, got {int(arr.size)}."
            )
        out[name] = arr
    return out


def derivative_result_payload(
    *,
    supported: bool,
    backend: str,
    message: str,
    value: Any,
    jacobian: Any,
    shape: list[int] | None = None,
    outputs: Any = None,
    variables: Any = None,
    **extra: Any,
) -> dict[str, Any]:
    """Build the public derivative-result payload shape."""
    value_arr = np.asarray(value, dtype=float)
    jac_arr = np.asarray(jacobian, dtype=float)
    if shape is None:
        if jac_arr.ndim == 2:
            shape = [int(jac_arr.shape[0]), int(jac_arr.shape[1])]
        elif value_arr.ndim == 0:
            shape = [1, int(jac_arr.size)]
        else:
            shape = [int(value_arr.size), int(jac_arr.size // max(int(value_arr.size), 1))]
    if outputs is None:
        outputs = extra.get("output_order", extra.get("component_order", ()))
    if variables is None:
        variables = extra.get("variable_order", extra.get("parameter_order", ()))
    payload = {
        "supported": bool(supported),
        "backend": str(backend),
        "derivative_backend": str(backend),
        "message": str(message),
        "value": value_arr,
        "jacobian": jac_arr.reshape(tuple(shape)) if bool(supported) and jac_arr.size else jac_arr,
        "outputs": [str(item) for item in (outputs or [])],
        "variables": [str(item) for item in (variables or [])],
        "shape": [int(shape[0]), int(shape[1])],
    }
    payload.update(extra)
    return payload


def backend_from_contribution_details(
    details: Mapping[str, Any], *, available: bool = True, reason: str = ""
) -> tuple[str, str]:
    """Collapse contribution-family derivative details into the public summary backend."""
    if not available:
        unsupported_derivative(reason or "composition derivative contribution backend is unsupported")
    labels = {str(v) for v in dict(details).values()}
    if "unsupported" in labels:
        unsupported_derivative(reason or "composition derivative contribution backend is unsupported")
    if "cppad" in labels:
        return "cppad", "mixed CppAD/analytic derivative contributions"
    if "cppad_implicit" in labels:
        return "cppad_implicit", "mixed implicit CppAD/analytic derivative contributions"
    if "analytic_implicit" in labels:
        return "analytic_implicit", "mixed implicit analytic/analytic derivative contributions"
    return "analytic", "analytic derivative contributions"


def unsupported_derivative(message: str) -> None:
    """Raise the public unsupported-derivative error shape."""
    text = str(message)
    if not text.startswith("unsupported:"):
        text = f"unsupported: {text}"
    raise InputError(text)


def composition_derivative_residual_helmholtz_result(state: Any) -> dict[str, Any]:
    """Return residual-Helmholtz composition derivative terms for an EOS state."""
    ncomp = int(state._x.size)
    result = state._native.composition_derivative_residual_helmholtz_result()
    dadx = vector_terms_dict(result.dadx, ncomp, "dadx")
    ares = scalar_terms_dict(result.ares)
    sum_x = scalar_terms_dict(result.sum_x_dadx)
    z_raw = scalar_terms_dict(result.z_raw)
    z_terms = scalar_terms_dict(result.z)
    terms = {name: np.asarray(dadx[name], dtype=float) for name in CONTRIBUTION_NAMES}
    return {
        "total": sum_vector_terms(terms),
        "terms": terms,
        "ares_terms": {name: float(ares[name]) for name in CONTRIBUTION_NAMES},
        "sum_x_terms": {name: float(sum_x[name]) for name in CONTRIBUTION_NAMES},
        "z_raw_terms": {name: float(z_raw[name]) for name in CONTRIBUTION_NAMES},
        "z_terms": {name: float(z_terms[name]) for name in CONTRIBUTION_NAMES},
        "z_total": float(z_terms["total"]),
        "derivative_backend": {str(k): str(v) for k, v in dict(result.derivative_backend).items()},
        "derivative_available": bool(result.derivative_available),
    }


def fugacity_coefficient_term_result(state: Any) -> dict[str, Any]:
    """Return fugacity-coefficient contribution terms for an EOS state."""
    ncomp = int(state._x.size)
    result = state._native.fugacity_coefficient_result()
    mu = vector_terms_dict(result.mu, ncomp, "mu")
    lnfug = vector_terms_dict(result.lnfugcoef, ncomp, "lnfugcoef")
    dadx = vector_terms_dict(result.composition.dadx, ncomp, "dadx")
    ares = scalar_terms_dict(result.composition.ares)
    sum_x = scalar_terms_dict(result.composition.sum_x_dadx)
    z_raw = scalar_terms_dict(result.composition.z_raw)
    z_terms = scalar_terms_dict(result.composition.z)
    return {
        "mu_hc": mu["hc"],
        "mu_disp": mu["disp"],
        "mu_assoc": mu["assoc"],
        "mu_ion": mu["ion"],
        "mu_born": mu["born"],
        "mu_total": mu["total"],
        "lnfugcoef_total": lnfug["total"],
        "lnfugcoef_hc": lnfug["hc"],
        "lnfugcoef_disp": lnfug["disp"],
        "lnfugcoef_assoc": lnfug["assoc"],
        "lnfugcoef_ion": lnfug["ion"],
        "lnfugcoef_born": lnfug["born"],
        "dadx_hc": dadx["hc"],
        "dadx_disp": dadx["disp"],
        "dadx_assoc": dadx["assoc"],
        "dadx_ion": dadx["ion"],
        "dadx_born": dadx["born"],
        "a_hc": ares["hc"],
        "a_disp": ares["disp"],
        "a_assoc": ares["assoc"],
        "a_ion": ares["ion"],
        "a_born": ares["born"],
        "sum_x_dadx_hc": sum_x["hc"],
        "sum_x_dadx_disp": sum_x["disp"],
        "sum_x_dadx_assoc": sum_x["assoc"],
        "sum_x_dadx_ion": sum_x["ion"],
        "sum_x_dadx_born": sum_x["born"],
        "z_raw_hc": z_raw["hc"],
        "z_raw_disp": z_raw["disp"],
        "z_raw_assoc": z_raw["assoc"],
        "z_raw_ion": z_raw["ion"],
        "z_raw_born": z_raw["born"],
        "z_hc": z_terms["hc"],
        "z_disp": z_terms["disp"],
        "z_assoc": z_terms["assoc"],
        "z_ion": z_terms["ion"],
        "z_born": z_terms["born"],
        "z_total": z_terms["total"],
    }


def state_diagnostics_payload(state: Any, species: Any = None) -> dict[str, Any]:
    """Return the stable serialized state diagnostics payload."""
    mix = state._mixture
    species = state._mixture.species if species is None else species
    z = np.asarray(mix._params.get("z", []), dtype=float).flatten()
    has_ions = bool(np.any(np.abs(z) > 1e-12))
    terms = fugacity_coefficient_term_result(state)
    fugacity_coefficient = state.fugacity_coefficient(natural_log=False)
    if has_ions:
        activity_coefficient = state.activity_coefficient(species=species, mean_ionic_form=False)
        relative_permittivity = state.relative_permittivity()
        osmotic_coefficient = state.osmotic_coefficient()
        mean_ionic_activity_coefficient_molality = state.activity_coefficient(
            species=species,
            mean_ionic_form=True,
            basis="molality",
        )
        mean_ionic_activity_coefficient_mole = state.activity_coefficient(
            species=species,
            mean_ionic_form=True,
            basis="mole",
        )
        solvation_free_energy = state.solvation_free_energy(species=species)
    else:
        activity_coefficient = {}
        relative_permittivity = None
        osmotic_coefficient = None
        mean_ionic_activity_coefficient_molality = {}
        mean_ionic_activity_coefficient_mole = {}
        solvation_free_energy = {}
    return {
        "T": state._T,
        "phase": state._phase,
        "x": np.asarray(state._x, dtype=float),
        "pressure": state.pressure(),
        "density": state.molar_density(),
        "density_molar": state.molar_density(),
        "mass_density": state.mass_density() if np.asarray(mix._params.get("MW", []), dtype=float).size else None,
        "compressibility_factor": state.compressibility_factor(),
        "residual_helmholtz": state.residual_helmholtz(),
        "residual_enthalpy": state.residual_enthalpy(),
        "residual_entropy": state.residual_entropy(),
        "residual_gibbs": state.residual_gibbs(),
        "residual_chemical_potential": state.residual_chemical_potential(),
        "fugacity_coefficient": fugacity_coefficient,
        "activity_coefficient": activity_coefficient,
        "fugacity_coefficient_terms": terms,
        "relative_permittivity": relative_permittivity,
        "osmotic_coefficient": osmotic_coefficient,
        "mean_ionic_activity_coefficient_molality": mean_ionic_activity_coefficient_molality,
        "mean_ionic_activity_coefficient_mole": mean_ionic_activity_coefficient_mole,
        "solvation_free_energy": solvation_free_energy,
    }
