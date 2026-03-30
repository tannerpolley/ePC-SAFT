from __future__ import annotations

import argparse
import csv
import json
import subprocess
import sys
import tempfile
from dataclasses import dataclass
from pathlib import Path

import numpy as np

REPO_ROOT = Path(__file__).resolve().parents[2]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from scripts.multiphase_model_analysis.ascani_case2_dataset_comparison import PAPER_TARGETS, _case2_feed


OUTPUT_DIR = Path(__file__).resolve().parent / "output"
CLAPEYRON_ROOT = Path(r"C:\Users\Tanner\Documents\git\Clapeyron.jl")
FEOS_ROOT = Path(r"C:\Users\Tanner\Documents\git\feos")
CLAPEYRON_RUNNER = Path(__file__).resolve().with_name("external_package_case2_clapeyron.jl")
CLAPEYRON_UNLIKE = CLAPEYRON_ROOT / "database" / "SAFT" / "PCSAFT" / "ePCSAFT" / "ePCSAFT_unlike.csv"
FEOS_EPC_PURE = FEOS_ROOT / "parameters" / "epcsaft" / "held2014_w_permittivity_added.json"
FEOS_EPC_BINARY = FEOS_ROOT / "parameters" / "epcsaft" / "held2014_binary.json"
FEOS_PCSAFT_PURE = FEOS_ROOT / "parameters" / "pcsaft" / "gross2002.json"
FEOS_PCSAFT_BINARY = FEOS_ROOT / "parameters" / "pcsaft" / "rehner2023_binary.json"


@dataclass(frozen=True)
class CaseSpec:
    key: str
    label: str
    neutral_components: tuple[str, str]
    ion_components: tuple[str, str, str]
    charges: tuple[int, int, int, int, int]
    salt_labels: tuple[str, str]
    feed_moles: tuple[float, ...] | None = None
    paper_targets: dict[str, float] | None = None

    @property
    def cation_indices(self) -> tuple[int, int]:
        return (2, 3)

    @property
    def chloride_index(self) -> int:
        return 4

    @property
    def component_labels(self) -> tuple[str, ...]:
        return (
            self.neutral_components[0],
            self.neutral_components[1],
            self.ion_components[0],
            self.ion_components[1],
            self.ion_components[2],
        )


_species, _feed_z, _mass_feed = _case2_feed()
WORKED_EXAMPLE = CaseSpec(
    key="worked_example_water_butanol_nacl_kcl",
    label="Water + 1-butanol + NaCl + KCl worked example",
    neutral_components=("water", "1-butanol"),
    ion_components=("sodium", "potassium", "chloride"),
    charges=(0, 0, 1, 1, -1),
    salt_labels=("NaCl", "KCl"),
    feed_moles=tuple(
        float(value)
        for value in np.array(
            [
                0.8094 / 18.01528,
                0.1728 / 74.123,
                0.0054 / 58.44,
                0.0124 / 74.5513,
                0.0054 / 58.44 + 0.0124 / 74.5513,
            ],
            dtype=float,
        )
    ),
    paper_targets=PAPER_TARGETS,
)

CASES: tuple[CaseSpec, ...] = (
    WORKED_EXAMPLE,
    CaseSpec(
        key="water_butanol_kcl_nh4cl",
        label="Water + 1-butanol + KCl + NH4Cl",
        neutral_components=("water", "1-butanol"),
        ion_components=("potassium", "ammonium", "chloride"),
        charges=(0, 0, 1, 1, -1),
        salt_labels=("KCl", "NH4Cl"),
    ),
    CaseSpec(
        key="water_propanol_nacl_kcl",
        label="Water + 1-propanol + NaCl + KCl",
        neutral_components=("water", "1-propanol"),
        ion_components=("sodium", "potassium", "chloride"),
        charges=(0, 0, 1, 1, -1),
        salt_labels=("NaCl", "KCl"),
    ),
    CaseSpec(
        key="water_propanol_kcl_nh4cl",
        label="Water + 1-propanol + KCl + NH4Cl",
        neutral_components=("water", "1-propanol"),
        ion_components=("potassium", "ammonium", "chloride"),
        charges=(0, 0, 1, 1, -1),
        salt_labels=("KCl", "NH4Cl"),
    ),
)


def _load_json(path: Path) -> object:
    return json.loads(path.read_text(encoding="utf-8"))


def _load_csv_rows(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle))


def _load_preamble_csv_rows(path: Path, header_index: int) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        rows = [row for row in csv.reader(handle) if row]
    header = rows[header_index]
    return [dict(zip(header, row, strict=False)) for row in rows[header_index + 1 :]]


def charge_residual(molefracs: np.ndarray, charges: tuple[int, ...]) -> float:
    return float(np.dot(np.asarray(molefracs, dtype=float), np.asarray(charges, dtype=float)))


def map_species_to_pseudo_salt_basis(
    molefracs: np.ndarray,
    case: CaseSpec,
    *,
    check_charge: bool = True,
    charge_tol: float = 1.0e-6,
) -> dict[str, float]:
    x = np.asarray(molefracs, dtype=float)
    if check_charge:
        residual = charge_residual(x, case.charges)
        if abs(residual) > charge_tol:
            raise ValueError(f"Charge residual {residual:+.3e} exceeds tolerance for {case.key}.")
    values = {
        "water": float(x[0]),
        "alcohol": float(x[1]),
        case.salt_labels[0]: float(x[case.cation_indices[0]]),
        case.salt_labels[1]: float(x[case.cation_indices[1]]),
    }
    total = sum(values.values())
    return {key: value / total for key, value in values.items()}


def mean_ionic_lnfugacities(lnfug_bar: np.ndarray, case: CaseSpec) -> dict[str, float]:
    lnf = np.asarray(lnfug_bar, dtype=float)
    chloride_value = float(lnf[case.chloride_index])
    return {
        salt: 0.5 * (float(lnf[idx]) + chloride_value)
        for salt, idx in zip(case.salt_labels, case.cation_indices, strict=True)
    }


def _pseudo_salt_comparison_rows(clap_result: dict[str, object], case: CaseSpec) -> list[dict[str, object]]:
    if case.paper_targets is None:
        return []

    organic = clap_result["organic"]
    aqueous = clap_result["aqueous"]
    lnf_mean = clap_result["mean_lnfugacity_bar"]
    targets = case.paper_targets
    mapping = [
        ("x_water_org", "organic", "water"),
        ("x_butanol_org", "organic", "alcohol"),
        ("x_nacl_org", "organic", "NaCl"),
        ("x_kcl_org", "organic", "KCl"),
        ("x_water_aq", "aqueous", "water"),
        ("x_butanol_aq", "aqueous", "alcohol"),
        ("x_nacl_aq", "aqueous", "NaCl"),
        ("x_kcl_aq", "aqueous", "KCl"),
    ]
    rows: list[dict[str, object]] = []
    for key, phase, component in mapping:
        package_value = (
            float(organic["pseudo_salt_basis"][component])
            if phase == "organic"
            else float(aqueous["pseudo_salt_basis"][component])
        )
        paper_value = float(targets[key])
        rows.append(
            {
                "quantity": key,
                "paper": paper_value,
                "package": package_value,
                "abs_delta": abs(package_value - paper_value),
                "rel_delta_pct": np.nan if abs(paper_value) <= 1.0e-16 else abs(package_value - paper_value) / abs(paper_value) * 100.0,
            }
        )

    for key, salt in (("lnf_water_bar", "water"), ("lnf_butanol_bar", "alcohol")):
        if salt == "water":
            package_value = float(clap_result["lnfugacity_bar_avg"]["water"])
        else:
            package_value = float(clap_result["lnfugacity_bar_avg"]["alcohol"])
        paper_value = float(targets[key])
        rows.append(
            {
                "quantity": key,
                "paper": paper_value,
                "package": package_value,
                "abs_delta": abs(package_value - paper_value),
                "rel_delta_pct": abs(package_value - paper_value) / abs(paper_value) * 100.0,
            }
        )

    return rows


def _clapeyron_method_configs() -> list[dict[str, object]]:
    return [
        {
            "label": "MichelsenTPFlash lle K0_A",
            "K0": [1.0e-3, 1.0e2, 1.0e3, 1.0e3, 1.0e3],
        },
        {
            "label": "MichelsenTPFlash lle K0_B",
            "K0": [1.0e-2, 1.0e1, 1.0e3, 1.0e3, 1.0e3],
        },
        {
            "label": "MichelsenTPFlash lle K0_C",
            "K0": [1.0e-1, 1.0e1, 1.0e2, 1.0e2, 1.0e2],
        },
    ]


def _run_clapeyron_case(case: CaseSpec, output_dir: Path) -> dict[str, object]:
    action = "flash" if case.feed_moles is not None else "build_only"
    input_payload = {
        "action": action,
        "case_key": case.key,
        "label": case.label,
        "neutral_components": list(case.neutral_components),
        "ion_components": list(case.ion_components),
        "charges": list(case.charges),
        "pressure_pa": 1.0e5,
        "temperature_k": 298.15,
        "method_configs": _clapeyron_method_configs(),
    }
    if case.feed_moles is not None:
        input_payload["feed_moles"] = list(case.feed_moles)

    output_path = output_dir / f"{case.key}_clapeyron_raw.json"
    with tempfile.TemporaryDirectory(prefix="pcsaft_clapeyron_case2_") as tmpdir:
        input_path = Path(tmpdir) / "input.json"
        input_path.write_text(json.dumps(input_payload, indent=2), encoding="utf-8")
        cmd = [
            "julia",
            f"--project={CLAPEYRON_ROOT}",
            "--startup-file=no",
            str(CLAPEYRON_RUNNER),
            str(input_path),
            str(output_path),
        ]
        completed = subprocess.run(
            cmd,
            cwd=REPO_ROOT,
            capture_output=True,
            text=True,
            check=True,
        )
    result = _load_json(output_path)
    if completed.stdout.strip():
        result["stdout"] = completed.stdout.strip()
    if completed.stderr.strip():
        result["stderr"] = completed.stderr.strip()
    return result


def _clapeyron_unlike_audit(case: CaseSpec) -> dict[str, object]:
    rows = _load_preamble_csv_rows(CLAPEYRON_UNLIKE, header_index=2)
    names = {case.neutral_components[1], *case.ion_components}
    relevant = []
    for row in rows:
        species1 = row["species1"].split("~|~")[0]
        species2 = row["species2"].split("~|~")[0]
        if species1 in names and species2 in names:
            relevant.append({"species1": species1, "species2": species2, "k": row["k"]})
    alcohol_ion_pairs = [
        row
        for row in relevant
        if case.neutral_components[1] in (row["species1"], row["species2"])
        and any(ion in (row["species1"], row["species2"]) for ion in case.ion_components)
    ]
    return {
        "relevant_unlike_rows": relevant,
        "alcohol_ion_rows": alcohol_ion_pairs,
        "alcohol_ion_row_count": len(alcohol_ion_pairs),
    }


def _postprocess_clapeyron_flash(raw: dict[str, object], case: CaseSpec) -> dict[str, object]:
    flash = raw.get("flash", {})
    if not flash or not flash.get("success"):
        return {"status": "flash_failed", "raw": raw}

    xmat = np.asarray(flash["phase_molefracs"], dtype=float)
    nmat = np.asarray(flash["phase_component_moles"], dtype=float)
    organic_row = int(flash["organic_phase_index"]) - 1
    aqueous_row = int(flash["aqueous_phase_index"]) - 1
    betas = np.asarray(flash["phase_betas"], dtype=float)
    org_lnf = np.asarray(flash["organic_lnfug_bar"], dtype=float)
    aq_lnf = np.asarray(flash["aqueous_lnfug_bar"], dtype=float)
    org_basis = map_species_to_pseudo_salt_basis(xmat[organic_row], case)
    aq_basis = map_species_to_pseudo_salt_basis(xmat[aqueous_row], case)
    mean_org = mean_ionic_lnfugacities(org_lnf, case)
    mean_aq = mean_ionic_lnfugacities(aq_lnf, case)
    neutral_gaps = {
        "water": float(org_lnf[0] - aq_lnf[0]),
        "alcohol": float(org_lnf[1] - aq_lnf[1]),
    }
    salt_gaps = {
        salt: float(mean_org[salt] - mean_aq[salt]) for salt in case.salt_labels
    }
    charge_residuals = {
        "organic": charge_residual(xmat[organic_row], case.charges),
        "aqueous": charge_residual(xmat[aqueous_row], case.charges),
    }
    comparison_basis = {
        "organic": {"pseudo_salt_basis": org_basis},
        "aqueous": {"pseudo_salt_basis": aq_basis},
        "mean_lnfugacity_bar": {
            salt: {"avg": float(0.5 * (mean_org[salt] + mean_aq[salt]))}
            for salt in case.salt_labels
        },
        "lnfugacity_bar_avg": {
            "water": float(0.5 * (org_lnf[0] + aq_lnf[0])),
            "alcohol": float(0.5 * (org_lnf[1] + aq_lnf[1])),
        },
    }

    result = {
        "status": "ok",
        "method": flash["method"],
        "components": flash["components"],
        "phase_betas": {
            "organic": float(betas[organic_row]),
            "aqueous": float(betas[aqueous_row]),
        },
        "organic": {
            "species_basis": dict(zip(case.component_labels, xmat[organic_row].tolist(), strict=True)),
            "component_moles": dict(zip(case.component_labels, nmat[organic_row].tolist(), strict=True)),
            "pseudo_salt_basis": org_basis,
            "lnfug_bar": {
                "water": float(org_lnf[0]),
                "alcohol": float(org_lnf[1]),
                **{salt: float(value) for salt, value in mean_org.items()},
            },
            "charge_residual": charge_residuals["organic"],
        },
        "aqueous": {
            "species_basis": dict(zip(case.component_labels, xmat[aqueous_row].tolist(), strict=True)),
            "component_moles": dict(zip(case.component_labels, nmat[aqueous_row].tolist(), strict=True)),
            "pseudo_salt_basis": aq_basis,
            "lnfug_bar": {
                "water": float(aq_lnf[0]),
                "alcohol": float(aq_lnf[1]),
                **{salt: float(value) for salt, value in mean_aq.items()},
            },
            "charge_residual": charge_residuals["aqueous"],
        },
        "lnfugacity_bar_avg": {
            "water": float(0.5 * (org_lnf[0] + aq_lnf[0])),
            "alcohol": float(0.5 * (org_lnf[1] + aq_lnf[1])),
        },
        "mean_lnfugacity_bar": {
            salt: {
                "organic": float(mean_org[salt]),
                "aqueous": float(mean_aq[salt]),
                "avg": float(0.5 * (mean_org[salt] + mean_aq[salt])),
                "gap": float(mean_org[salt] - mean_aq[salt]),
            }
            for salt in case.salt_labels
        },
        "max_neutral_gap_abs": max(abs(value) for value in neutral_gaps.values()),
        "max_mean_ionic_gap_abs": max(abs(value) for value in salt_gaps.values()),
        "paper_comparison": _pseudo_salt_comparison_rows(comparison_basis, case),
        "raw": raw,
    }
    return result


def _find_record(rows: list[dict[str, object]], name1: str, name2: str) -> dict[str, object] | None:
    for row in rows:
        first = row["id1"]["name"]
        second = row["id2"]["name"]
        if {first, second} == {name1, name2}:
            return row
    return None


def _feos_case_audit(case: CaseSpec) -> dict[str, object]:
    pure_rows_epc = _load_json(FEOS_EPC_PURE)
    binary_rows_epc = _load_json(FEOS_EPC_BINARY)
    pure_rows_pc = _load_json(FEOS_PCSAFT_PURE)
    binary_rows_pc = _load_json(FEOS_PCSAFT_BINARY)
    alcohol = case.neutral_components[1]

    def _has_pure(rows: list[dict[str, object]], name: str) -> bool:
        return any(row["identifier"]["name"] == name for row in rows)

    water_alcohol_binary = _find_record(binary_rows_pc, "water", alcohol)
    alcohol_ion_binary_hits = []
    for ion in case.ion_components:
        hit = _find_record(binary_rows_epc, alcohol, f"{ion} ion" if not ion.endswith("ion") else ion)
        if hit is not None:
            alcohol_ion_binary_hits.append(hit)

    return {
        "alcohol_pure_present_in_pcsaft": _has_pure(pure_rows_pc, alcohol),
        "water_pure_present_in_epcsaft": _has_pure(pure_rows_epc, "water"),
        "ion_pure_present_in_epcsaft": {
            ion: _has_pure(pure_rows_epc, f"{ion} ion")
            for ion in case.ion_components
        },
        "water_alcohol_binary_present_in_pcsaft": water_alcohol_binary is not None,
        "alcohol_ion_binary_present_in_epcsaft": len(alcohol_ion_binary_hits) > 0,
        "alcohol_ion_binary_hits": alcohol_ion_binary_hits,
        "alcohol_has_permittivity_record_in_epcsaft": _has_pure(pure_rows_epc, alcohol),
    }


def _run_feos_case(case: CaseSpec) -> dict[str, object]:
    import feos

    alcohol = case.neutral_components[1]
    ion_names = [f"{name} ion" for name in case.ion_components]
    pure_input = [
        (["water", *ion_names], str(FEOS_EPC_PURE)),
        ([alcohol], str(FEOS_PCSAFT_PURE)),
    ]
    audit = _feos_case_audit(case)
    results: dict[str, object] = {"audit": audit, "variants": {}}
    for variant in ("advanced", "revised"):
        try:
            params = feos.Parameters.from_multiple_json(pure_input, str(FEOS_EPC_BINARY))
            feos.EquationOfState.epcsaft(params, epcsaft_variant=variant)
            results["variants"][variant] = {"success": True}
        except Exception as exc:  # pragma: no cover - exercised in runtime integration
            results["variants"][variant] = {
                "success": False,
                "error": str(exc),
            }
    return results


def _build_markdown_report(summary: dict[str, object]) -> str:
    worked = summary["worked_example"]
    clap = worked["clapeyron"]
    feos = worked["feos"]
    lines = [
        "# External Package Ascani Case-2 Comparison",
        "",
        "## Worked example",
        "",
        f"- System: {WORKED_EXAMPLE.label}",
        f"- Feed mass fractions: water={_mass_feed['w_water']:.4f}, 1-butanol={_mass_feed['w_butanol']:.4f}, NaCl={_mass_feed['w_nacl']:.4f}, KCl={_mass_feed['w_kcl']:.4f}",
        f"- Clapeyron no-charge build success: {clap['raw']['build_without_charge']['success']}",
        f"- Clapeyron explicit-charge build success: {clap['raw']['build_with_charge']['success']}",
    ]
    if clap["status"] == "ok":
        lines.extend(
            [
                f"- Clapeyron flash method: `{clap['method']}`",
                f"- Clapeyron phase betas: organic={clap['phase_betas']['organic']:.6f}, aqueous={clap['phase_betas']['aqueous']:.6f}",
                f"- Clapeyron max neutral ln fugacity gap: {clap['max_neutral_gap_abs']:.3e}",
                f"- Clapeyron max mean-ionic ln fugacity gap: {clap['max_mean_ionic_gap_abs']:.3e}",
                f"- feos advanced build success: {feos['variants']['advanced']['success']}",
                f"- feos revised build success: {feos['variants']['revised']['success']}",
            ]
        )
    else:
        lines.append("- Clapeyron flash did not succeed.")
    lines.extend(
        [
            "",
            "## Worked example vs paper",
            "",
            "| Quantity | Paper | Package | |Δ| |",
            "|---|---:|---:|---:|",
        ]
    )
    for row in clap.get("paper_comparison", []):
        lines.append(
            f"| {row['quantity']} | {row['paper']:.6g} | {row['package']:.6g} | {row['abs_delta']:.6g} |"
        )
    lines.extend(
        [
            "",
            "## Capability audit",
            "",
            "| System | Clapeyron build no charge | Clapeyron build explicit charge | feos advanced | feos revised |",
            "|---|---:|---:|---:|---:|",
        ]
    )
    for row in summary["system_audit"]:
        lines.append(
            f"| {row['label']} | {row['clapeyron_no_charge']} | {row['clapeyron_with_charge']} | "
            f"{row['feos_advanced']} | {row['feos_revised']} |"
        )
    lines.extend(
        [
            "",
            "## Notes",
            "",
            f"- Clapeyron alcohol-ion unlike rows for worked example: {worked['clapeyron_unlike_audit']['alcohol_ion_row_count']}.",
            f"- feos alcohol permittivity present for worked example alcohol: {worked['feos']['audit']['alcohol_has_permittivity_record_in_epcsaft']}.",
            f"- feos water-alcohol binary present in shipped PC-SAFT data: {worked['feos']['audit']['water_alcohol_binary_present_in_pcsaft']}.",
            f"- feos alcohol-ion binary rows for worked example: {worked['feos']['audit']['alcohol_ion_binary_present_in_epcsaft']}.",
            (
                "- Clapeyron mean-ionic ln fugacities are captured in the raw JSON, but the package-accessible ionic "
                "reference basis does not appear to match the paper's tabulated mean-ionic fugacity basis closely "
                "enough for a direct numeric comparison."
            ),
        ]
    )
    return "\n".join(lines) + "\n"


def build_summary(output_dir: Path) -> dict[str, object]:
    output_dir.mkdir(parents=True, exist_ok=True)
    clapeyron_by_case = {case.key: _run_clapeyron_case(case, output_dir) for case in CASES}
    feos_by_case = {case.key: _run_feos_case(case) for case in CASES}

    worked_clap = _postprocess_clapeyron_flash(clapeyron_by_case[WORKED_EXAMPLE.key], WORKED_EXAMPLE)
    worked_feos = feos_by_case[WORKED_EXAMPLE.key]

    system_audit = []
    for case in CASES:
        clap_raw = clapeyron_by_case[case.key]
        feos_raw = feos_by_case[case.key]
        system_audit.append(
            {
                "key": case.key,
                "label": case.label,
                "clapeyron_no_charge": clap_raw["build_without_charge"]["success"],
                "clapeyron_with_charge": clap_raw["build_with_charge"]["success"],
                "feos_advanced": feos_raw["variants"]["advanced"]["success"],
                "feos_revised": feos_raw["variants"]["revised"]["success"],
            }
        )

    return {
        "worked_example": {
            "label": WORKED_EXAMPLE.label,
            "feed_species_basis": dict(zip(_species, _feed_z.tolist(), strict=True)),
            "feed_mass_fractions": _mass_feed,
            "paper_targets": PAPER_TARGETS,
            "clapeyron": worked_clap,
            "clapeyron_unlike_audit": _clapeyron_unlike_audit(WORKED_EXAMPLE),
            "feos": worked_feos,
        },
        "system_audit": system_audit,
        "raw": {
            "clapeyron": clapeyron_by_case,
            "feos": feos_by_case,
        },
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=OUTPUT_DIR,
        help="Directory for raw JSON and markdown summary output.",
    )
    args = parser.parse_args()

    summary = build_summary(args.output_dir)
    summary_json = args.output_dir / "external_package_case2_summary.json"
    summary_md = args.output_dir / "external_package_case2_summary.md"
    summary_json.write_text(json.dumps(summary, indent=2), encoding="utf-8")
    summary_md.write_text(_build_markdown_report(summary), encoding="utf-8")
    print(f"Wrote {summary_json}", flush=True)
    print(f"Wrote {summary_md}", flush=True)


if __name__ == "__main__":
    main()
