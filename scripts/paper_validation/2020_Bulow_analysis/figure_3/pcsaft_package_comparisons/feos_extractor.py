from __future__ import annotations

import json
from dataclasses import dataclass
from pathlib import Path

import numpy as np

try:
    import feos
    import si_units as si
except ImportError as exc:  # pragma: no cover - exercised by runtime setup
    raise RuntimeError(
        "feos and si_units must be installed in the PC-SAFT conda environment. "
        "Expected local editable install from C:\\Users\\Tanner\\Documents\\git\\feos\\py-feos."
    ) from exc


R_GAS = 8.31446261815324
T_REF = 298.15
P_REF = 1.0e5
EPS = 1.0e-8
EPS_INF = 1.0e-12

TERMS = ("hc", "disp", "assoc", "dh", "born")
LABEL_TO_TERM = {
    "Hard Sphere": "hc",
    "Hard Chain": "hc",
    "Dispersion": "disp",
    "Association": "assoc",
    "Ionic": "dh",
    "Born": "born",
}
RT_KJMOL = R_GAS * T_REF / 1000.0

FEOS_ROOT = Path(r"C:\Users\Tanner\Documents\git\feos")
FEOS_PARAMETER_PATH = FEOS_ROOT / "parameters" / "epcsaft" / "held2014_w_permittivity_added.json"
FEOS_BINARY_PATH = FEOS_ROOT / "parameters" / "epcsaft" / "held2014_binary.json"


@dataclass(frozen=True)
class IonSetup:
    components: tuple[str, str, str]
    target_index: int


ION_SETUPS: dict[str, IonSetup] = {
    "Li+": IonSetup(("water", "lithium ion", "chloride ion"), 1),
    "Na+": IonSetup(("water", "sodium ion", "chloride ion"), 1),
    "K+": IonSetup(("water", "potassium ion", "chloride ion"), 1),
    "F-": IonSetup(("water", "sodium ion", "fluoride ion"), 2),
    "Cl-": IonSetup(("water", "sodium ion", "chloride ion"), 2),
    "Br-": IonSetup(("water", "sodium ion", "bromide ion"), 2),
    "I-": IonSetup(("water", "sodium ion", "iodide ion"), 2),
}


def _to_kj_per_mol(quantity) -> float:
    return float(quantity / (si.KILO * si.JOULE / si.MOL))


def _to_pa(quantity) -> float:
    return float(quantity / si.PASCAL)


def _build_reference_state(eos, target_index: int):
    x_bulk = np.asarray([1.0 - 2.0 * EPS, EPS, EPS], dtype=float)
    state_bulk = feos.State(
        eos,
        temperature=T_REF * si.KELVIN,
        pressure=P_REF * si.PASCAL,
        molefracs=x_bulk,
        total_moles=si.MOL,
    )
    density_bulk = state_bulk.density

    x_ref = np.asarray([1.0, 0.0, 0.0], dtype=float)
    state_ref = feos.State(
        eos,
        temperature=T_REF * si.KELVIN,
        density=density_bulk,
        molefracs=x_ref,
        total_moles=si.MOL,
    )
    p_ref = state_ref.pressure()

    x_inf = x_ref.copy()
    x_inf[target_index] = EPS_INF
    x_inf[0] = 1.0 - EPS_INF
    state_inf = feos.State(
        eos,
        temperature=T_REF * si.KELVIN,
        pressure=p_ref,
        molefracs=x_inf,
        total_moles=si.MOL,
    )
    return state_inf, float(p_ref / si.PASCAL)


def _map_contributions(pairs, converter, *, ignore_labels: set[str] | None = None) -> dict[str, float]:
    values = {term: 0.0 for term in TERMS}
    seen_terms: set[str] = set()
    ignore_labels = set() if ignore_labels is None else set(ignore_labels)
    for label, quantity in pairs:
        term = LABEL_TO_TERM.get(label)
        if term is None:
            if label in ignore_labels:
                continue
            raise ValueError(f"Unexpected feos contribution label: {label!r}")
        values[term] += converter(quantity)
        seen_terms.add(term)

    missing = [term for term in TERMS if term not in seen_terms]
    if missing:
        raise ValueError(f"Missing feos contribution buckets: {missing}")
    return values


def _compute_one(ion: str, setup: IonSetup) -> dict[str, float]:
    parameters = feos.Parameters.from_json(
        list(setup.components),
        str(FEOS_PARAMETER_PATH),
        str(FEOS_BINARY_PATH),
    )
    eos = feos.EquationOfState.epcsaft(parameters, epcsaft_variant="advanced")
    state_inf, p_ref_pa = _build_reference_state(eos, setup.target_index)

    mu_pairs = state_inf.chemical_potential_contributions(setup.target_index, feos.Contributions.Residual)
    terms = _map_contributions(mu_pairs, _to_kj_per_mol)
    mu_total = _to_kj_per_mol(state_inf.chemical_potential(feos.Contributions.Residual)[setup.target_index])
    total = float(R_GAS * T_REF * state_inf.ln_phi()[setup.target_index] / 1000.0)
    mu_sum = float(sum(terms[term] for term in TERMS))
    pressure_terms = _map_contributions(state_inf.pressure_contributions(), _to_pa, ignore_labels={"Ideal gas"})
    pressure_pa = _to_pa(state_inf.pressure())
    volume_m3 = float(state_inf.volume / (si.METER**3))
    z_total = pressure_pa * volume_m3 / (R_GAS * T_REF)
    z_residual = z_total - 1.0
    if abs(z_residual) <= 1.0e-14:
        raise ValueError(f"Degenerate feos residual compressibility for {ion}.")

    z_corrections = {
        term: float(RT_KJMOL * (-(pressure_terms[term] * volume_m3 / (R_GAS * T_REF)) / z_residual * np.log(z_total)))
        for term in TERMS
    }
    lnfug_terms = {term: float(terms[term] + z_corrections[term]) for term in TERMS}
    lnfug_sum = float(sum(lnfug_terms[term] for term in TERMS))

    if not np.isfinite(total) or not np.isfinite(mu_total) or not np.isfinite(mu_sum):
        raise ValueError(f"Non-finite feos result for {ion}.")
    if not np.isfinite(lnfug_sum):
        raise ValueError(f"Non-finite feos lnfug-sum result for {ion}.")

    return {
        **terms,
        **{f"{term}_z_correction_kj_mol": z_corrections[term] for term in TERMS},
        **{f"{term}_lnfug_kj_mol": lnfug_terms[term] for term in TERMS},
        "mu_total_kj_mol": mu_total,
        "mu_sum_kj_mol": mu_sum,
        "mu_gap_kj_mol": mu_total - mu_sum,
        "lnfug_sum_kj_mol": lnfug_sum,
        "lnfug_gap_kj_mol": total - lnfug_sum,
        "total_kj_mol": total,
        "reference_pressure_pa": p_ref_pa,
        "state_pressure_pa": pressure_pa,
        "state_volume_m3": volume_m3,
        "compressibility_factor": z_total,
    }


def compute_results() -> dict[str, dict[str, float]]:
    return {ion: _compute_one(ion, setup) for ion, setup in ION_SETUPS.items()}


def write_results(path: Path) -> Path:
    results = compute_results()
    payload = {
        "package": "feos",
        "variant": "advanced",
        "parameter_path": str(FEOS_PARAMETER_PATH),
        "binary_path": str(FEOS_BINARY_PATH),
        "module_path": str(Path(feos.__file__).resolve()),
        "results": results,
    }
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2, sort_keys=True), encoding="utf-8")
    return path


def main() -> None:
    output_path = Path(__file__).resolve().with_name("feos_raw.json")
    write_results(output_path)
    print(f"Wrote {output_path}", flush=True)


if __name__ == "__main__":
    main()
