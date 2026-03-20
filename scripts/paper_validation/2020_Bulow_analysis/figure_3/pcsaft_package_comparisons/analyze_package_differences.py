from __future__ import annotations

import csv
import json
import sys
from pathlib import Path

import numpy as np

SCRIPT_DIR = Path(__file__).resolve().parent
FIGURE3_DIR = SCRIPT_DIR.parent
ANALYSIS_ROOT = FIGURE3_DIR.parent
REPO_ROOT = ANALYSIS_ROOT.parents[2]

if str(ANALYSIS_ROOT) not in sys.path:
    sys.path.insert(0, str(ANALYSIS_ROOT))
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from scripts._env import require_pcsaft_install

require_pcsaft_install()

import _model_overlay as overlay
import _plot_common as common
from pcsaft.parameters import get_prop_dict


DATA_PATH = FIGURE3_DIR / "data" / "water_contributions.csv"
FIGURE2_TOTALS_PATH = ANALYSIS_ROOT / "figure_2" / "data" / "water_comparisons.csv"
FEOS_RAW_PATH = SCRIPT_DIR / "feos_raw.json"
CLAPEYRON_RAW_PATH = SCRIPT_DIR / "clapeyron_raw.json"
OUTPUT_DIR = SCRIPT_DIR / "diagnostics"
BREAKDOWN_CSV = OUTPUT_DIR / "figure3_package_mu_breakdown_comparison.csv"
REPORT_MD = OUTPUT_DIR / "package_difference_report.md"

FEOS_PURE_PATH = Path(r"C:\Users\Tanner\Documents\git\feos\parameters\epcsaft\held2014_w_permittivity_added.json")
FEOS_BINARY_PATH = Path(r"C:\Users\Tanner\Documents\git\feos\parameters\epcsaft\held2014_binary.json")
CLAPEYRON_ROOT = Path(r"C:\Users\Tanner\Documents\git\Clapeyron.jl")

TARGET_IONS = ("Na+", "Cl-")
CONTRIBUTION_MAP = {
    "hc": {"paper_rows": ("hc avg", "hc"), "suffix": "hc"},
    "disp": {"paper_rows": ("disp avg", "disp"), "suffix": "disp"},
    "assoc": {"paper_rows": ("assoc avg", "assoc"), "suffix": "assoc"},
    "dh": {"paper_rows": (), "suffix": "ion"},
    "born": {"paper_rows": ("born avg", "born"), "suffix": "born"},
}
TERM_ORDER = (*CONTRIBUTION_MAP.keys(), "total")
RT_KJMOL = overlay.R_GAS * overlay.T_REF / 1000.0


def _paper_value(frame: common.Table, contribution: str, ion: str) -> float:
    if contribution == "dh":
        return 0.0
    for row_key in CONTRIBUTION_MAP[contribution]["paper_rows"]:
        if row_key in frame.index:
            return frame.scalar(row_key, ion)
    raise KeyError(f"Missing paper row for contribution {contribution!r}.")


def _load_json(path: Path) -> dict:
    with path.open("r", encoding="utf-8") as handle:
        return json.load(handle)


def _load_clapeyron_csv(path: Path) -> list[dict[str, str]]:
    with path.open("r", newline="", encoding="utf-8-sig") as handle:
        rows = [row for row in csv.reader(handle) if row]
    header = rows[2]
    return [dict(zip(header, row, strict=False)) for row in rows[3:]]


def _current_repo_term_row(frame: common.Table, totals: common.Table, ion: str, contribution: str) -> dict[str, float | str]:
    terms, idx = overlay._infinite_dilution_terms("advanced", ion, "water")
    if contribution == "total":
        mu_sum = sum(
            RT_KJMOL * float(np.asarray(terms[f"mu_{meta['suffix']}"], dtype=float)[idx])
            for meta in CONTRIBUTION_MAP.values()
        )
        lnfug_sum = sum(
            RT_KJMOL * float(np.asarray(terms[f"lnfugcoef_{meta['suffix']}"], dtype=float)[idx])
            for meta in CONTRIBUTION_MAP.values()
        )
        return {
            "paper_mu_contr": totals.scalar("advanced", ion),
            "pcsaft_mu_contr": mu_sum,
            "pcsaft_mu_manual_sum": mu_sum,
            "pcsaft_a_contr": "",
            "pcsaft_z_contr": "",
            "pcsaft_dadx_contr": "",
            "pcsaft_sum_xj_dadx_contr": "",
            "pcsaft_z_correction": lnfug_sum - mu_sum,
            "pcsaft_lnfug_contr": lnfug_sum,
            "pcsaft_total": overlay.gsolv_ion("advanced", ion, "water"),
        }

    suffix = str(CONTRIBUTION_MAP[contribution]["suffix"])
    mu_term = float(np.asarray(terms[f"mu_{suffix}"], dtype=float)[idx])
    lnfug_term = float(np.asarray(terms[f"lnfugcoef_{suffix}"], dtype=float)[idx])
    a_term = float(terms[f"a_{suffix}"])
    z_term = float(terms[f"z_{suffix}"])
    dadx_term = float(np.asarray(terms[f"dadx_{suffix}"], dtype=float)[idx])
    sum_x_dadx_term = float(terms[f"sum_x_dadx_{suffix}"])

    pcsaft_mu = RT_KJMOL * mu_term
    return {
        "paper_mu_contr": _paper_value(frame, contribution, ion),
        "pcsaft_mu_contr": pcsaft_mu,
        "pcsaft_mu_manual_sum": RT_KJMOL * (a_term + z_term + dadx_term - sum_x_dadx_term),
        "pcsaft_a_contr": RT_KJMOL * a_term,
        "pcsaft_z_contr": RT_KJMOL * z_term,
        "pcsaft_dadx_contr": RT_KJMOL * dadx_term,
        "pcsaft_sum_xj_dadx_contr": RT_KJMOL * sum_x_dadx_term,
        "pcsaft_z_correction": RT_KJMOL * (lnfug_term - mu_term),
        "pcsaft_lnfug_contr": RT_KJMOL * lnfug_term,
        "pcsaft_total": overlay.gsolv_ion("advanced", ion, "water"),
    }


def _package_term_values(payload: dict, ion: str, contribution: str, prefix: str) -> dict[str, float]:
    row = payload["results"][ion]
    if contribution == "total":
        return {
            f"{prefix}_mu_contr": float(row["mu_sum_kj_mol"]),
            f"{prefix}_z_correction": float(row["lnfug_sum_kj_mol"]) - float(row["mu_sum_kj_mol"]),
            f"{prefix}_lnfug_contr": float(row["lnfug_sum_kj_mol"]),
            f"{prefix}_mu_total": float(row["mu_total_kj_mol"]),
            f"{prefix}_total": float(row["total_kj_mol"]),
        }
    return {
        f"{prefix}_mu_contr": float(row[contribution]),
        f"{prefix}_z_correction": float(row[f"{contribution}_z_correction_kj_mol"]),
        f"{prefix}_lnfug_contr": float(row[f"{contribution}_lnfug_kj_mol"]),
        f"{prefix}_mu_total": float(row["mu_total_kj_mol"]),
        f"{prefix}_total": float(row["total_kj_mol"]),
    }


def _build_breakdown_rows() -> list[dict[str, object]]:
    frame = common.load_indexed_csv(DATA_PATH)
    totals = common.load_indexed_csv(FIGURE2_TOTALS_PATH)
    feos_payload = _load_json(FEOS_RAW_PATH)
    clapeyron_payload = _load_json(CLAPEYRON_RAW_PATH)

    rows: list[dict[str, object]] = []
    for ion in TARGET_IONS:
        for contribution in TERM_ORDER:
            row = {"ion": ion, "contr": contribution}
            row.update(_current_repo_term_row(frame, totals, ion, contribution))
            row.update(_package_term_values(feos_payload, ion, contribution, "feos"))
            row.update(_package_term_values(clapeyron_payload, ion, contribution, "clapeyron"))
            row["feos_minus_pcsaft_mu"] = float(row["feos_mu_contr"]) - float(row["pcsaft_mu_contr"])
            row["clapeyron_minus_pcsaft_mu"] = float(row["clapeyron_mu_contr"]) - float(row["pcsaft_mu_contr"])
            rows.append(row)
    return rows


def _write_breakdown_csv(rows: list[dict[str, object]]) -> Path:
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    fieldnames = [
        "ion",
        "contr",
        "paper_mu_contr",
        "pcsaft_mu_contr",
        "pcsaft_mu_manual_sum",
        "pcsaft_a_contr",
        "pcsaft_z_contr",
        "pcsaft_dadx_contr",
        "pcsaft_sum_xj_dadx_contr",
        "pcsaft_z_correction",
        "pcsaft_lnfug_contr",
        "pcsaft_total",
        "feos_mu_contr",
        "feos_z_correction",
        "feos_lnfug_contr",
        "feos_mu_total",
        "feos_total",
        "clapeyron_mu_contr",
        "clapeyron_z_correction",
        "clapeyron_lnfug_contr",
        "clapeyron_mu_total",
        "clapeyron_total",
        "feos_minus_pcsaft_mu",
        "clapeyron_minus_pcsaft_mu",
    ]
    with BREAKDOWN_CSV.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)
    return BREAKDOWN_CSV


def _find_feos_pure(name: str) -> dict:
    data = json.loads(FEOS_PURE_PATH.read_text(encoding="utf-8"))
    return next(record for record in data if record["identifier"]["name"] == name)


def _find_feos_binary(name1: str, name2: str) -> dict:
    data = json.loads(FEOS_BINARY_PATH.read_text(encoding="utf-8"))
    for record in data:
        pair = {record["id1"]["name"], record["id2"]["name"]}
        if pair == {name1, name2}:
            return record
    raise KeyError(f"feos binary record not found for {name1!r}, {name2!r}")


def _current_repo_param_snapshot() -> dict[str, float]:
    x = np.asarray([overlay.EPS, overlay.EPS, 1.0 - 2.0 * overlay.EPS], dtype=float)
    params = get_prop_dict("bulow_2020", ["Na+", "Cl-", "Water"], x, overlay.T_REF)
    return {
        "na_sigma": float(params["s"][0]),
        "na_epsilon": float(params["e"][0]),
        "na_d_born": float(params["d_born"][0]),
        "cl_sigma": float(params["s"][1]),
        "cl_epsilon": float(params["e"][1]),
        "cl_d_born": float(params["d_born"][1]),
        "water_na_kij": float(params["k_ij"][2][0]),
        "water_cl_kij": float(params["k_ij"][2][1]),
        "na_cl_kij": float(params["k_ij"][0][1]),
    }


def _clapeyron_param_snapshot() -> dict[str, float]:
    unlike_adv = _load_clapeyron_csv(CLAPEYRON_ROOT / "database" / "SAFT" / "PCSAFT" / "ePCSAFTAdv" / "ePCSAFTAdv_unlike.csv")
    born_like = _load_clapeyron_csv(CLAPEYRON_ROOT / "database" / "Electrolytes" / "Born" / "born_like.csv")

    def born_radius(species: str) -> float:
        row = next(item for item in born_like if item["species"].startswith(f"{species}~|~"))
        return float(row["sigma_born"])

    def unlike_k(species1: str, species2: str) -> tuple[float, float]:
        row = next(
            item
            for item in unlike_adv
            if {item["species1"].split("~|~")[0], item["species2"].split("~|~")[0]} == {species1, species2}
        )
        return float(row["k"]), float(row["kT"] or 0.0)

    water_na_k, water_na_kT = unlike_k("water", "sodium")
    water_cl_k, water_cl_kT = unlike_k("water", "chloride")
    na_cl_k, _ = unlike_k("sodium", "chloride")
    return {
        "na_d_born": born_radius("sodium"),
        "cl_d_born": born_radius("chloride"),
        "water_na_k": water_na_k,
        "water_na_kT": water_na_kT,
        "water_cl_k": water_cl_k,
        "water_cl_kT": water_cl_kT,
        "na_cl_k": na_cl_k,
    }


def _feos_param_snapshot() -> dict[str, object]:
    na = _find_feos_pure("sodium ion")
    cl = _find_feos_pure("chloride ion")
    water_na = _find_feos_binary("water", "sodium ion")
    water_cl = _find_feos_binary("water", "chloride ion")
    na_cl = _find_feos_binary("sodium ion", "chloride ion")
    return {
        "na_sigma": float(na["sigma"]),
        "na_epsilon": float(na["epsilon_k"]),
        "na_dielectric": float(na["permittivity_record"]["ExperimentalData"]["data"][0][1]),
        "cl_sigma": float(cl["sigma"]),
        "cl_epsilon": float(cl["epsilon_k"]),
        "cl_dielectric": float(cl["permittivity_record"]["ExperimentalData"]["data"][0][1]),
        "water_na_kij": list(water_na["k_ij"]),
        "water_cl_kij": list(water_cl["k_ij"]),
        "na_cl_kij": list(na_cl["k_ij"]),
    }


def _render_report(rows: list[dict[str, object]]) -> str:
    current_params = _current_repo_param_snapshot()
    feos_params = _feos_param_snapshot()
    clapeyron_params = _clapeyron_param_snapshot()

    def find(ion: str, contr: str) -> dict[str, object]:
        return next(row for row in rows if row["ion"] == ion and row["contr"] == contr)

    lines = [
        "# Package Difference Report",
        "",
        "This report compares the Figure 3 contribution bookkeeping for `Na+` and `Cl-` across the current repo `PC-SAFT`, `feos`, and `Clapeyron.jl` after correcting the package-side total check to use the fugacity-basis sum where it is reconstructible.",
        "",
        "## Basis Checks",
        "",
    ]
    for ion in TARGET_IONS:
        total = find(ion, "total")
        lines.append(
            f"- `{ion}`: `pcsaft total - adjusted_sum = {float(total['pcsaft_total']) - float(total['pcsaft_lnfug_contr']):+.6f}` kJ/mol, "
            f"`feos total - adjusted_sum = {float(total['feos_total']) - float(total['feos_lnfug_contr']):+.6f}` kJ/mol, "
            f"`clapeyron total - adjusted_sum = {float(total['clapeyron_total']) - float(total['clapeyron_lnfug_contr']):+.6f}` kJ/mol."
        )
    lines.extend(
        [
            "",
            "The `pcsaft` and `Clapeyron` corrected sums close their totals to roundoff. `feos` does not: even after applying the same $-\\frac{Z^\\alpha}{Z-1}\\ln Z$ correction using its pressure contributions, the reconstructed sum still misses the package total by about `-40.962` kJ/mol for both `Na+` and `Cl-`. That means the exposed `feos` contribution labels are not equivalent to the current repo's per-term $\\mu^\\alpha/Z^\\alpha$ bookkeeping.",
            "",
            "## Term-Level Observations",
            "",
        ]
    )
    for ion in TARGET_IONS:
        for contr in ("hc", "disp", "assoc", "born"):
            row = find(ion, contr)
            lines.append(
                f"- `{ion}` `{contr}` mu term: `pcsaft = {float(row['pcsaft_mu_contr']):+.6f}`, "
                f"`feos = {float(row['feos_mu_contr']):+.6f}`, "
                f"`clapeyron = {float(row['clapeyron_mu_contr']):+.6f}` kJ/mol."
            )
    lines.extend(
        [
            "",
            "For `feos`, the association and Born `mu` terms stay very close to the current repo, while the hard-chain term is much more positive and the dispersion term shifts enough that the package-side branch sum no longer matches either `mu_total` or the corrected fugacity sum. That pattern points to a decomposition mismatch more than a total-EOS mismatch, because the overall hydration totals remain very close to the current repo.",
            "",
            "For `Clapeyron`, the corrected fugacity sum closes exactly, so its remaining differences are genuine model/state differences rather than a missing basis conversion. Its Born term is much less negative than the current repo for both `Na+` and `Cl-`, and its association term is more negative than both the current repo and `feos`.",
            "",
            "## Likely Causes",
            "",
            "### 1. `feos` branch bookkeeping is not equivalent to the repo's `mu` decomposition",
            "",
            "- `feos` uses `chemical_potential_contributions(...)` for branch terms, but the exposed branch set does not close to either its `mu_total` or its fugacity-basis total after the standard $Z$ correction.",
            "- The reconstructed `feos` fugacity-basis sum still misses the package total by the same offset as the package `mu_total - mu_sum` gap. That makes it unlikely that the remaining discrepancy is just a missing $Z$ correction.",
            f"- `feos` Born uses the hard-sphere diameter directly in [born.rs](/Users/Tanner/Documents/git/feos/crates/feos/src/epcsaft/eos/born.rs), not an explicit `d_Born` table like the current repo or Clapeyron.",
            "",
            "### 2. `Clapeyron.jl` is a genuinely different advanced-like model setup",
            "",
            "- The extractor uses `ESElectrolyte + pharmaPCSAFT + DHBorn + LinMixRSP`, which is the closest composable analogue in-tree, but it is not the same implementation as the current repo's native `bulow_2020` path.",
            f"- Current repo explicit Born diameters: `Na+ = {current_params['na_d_born']:.3f}` A, `Cl- = {current_params['cl_d_born']:.3f}` A.",
            f"- Clapeyron explicit Born diameters from `born_like.csv`: `Na+ = {clapeyron_params['na_d_born']:.3f}` A, `Cl- = {clapeyron_params['cl_d_born']:.3f}` A.",
            "- Those Born-radius differences directly explain part of the Born-term gap. The remaining association/hc/disp differences are consistent with the different neutralmodel implementation and state point that follow from the composable Clapeyron setup.",
            "",
            "### 3. Binary interaction parameters are not the main remaining cause for `feos`",
            "",
            f"- Current repo `bulow_2020` uses `k_ij(H2O,Na+) = {current_params['water_na_kij']:.6f}`, `k_ij(H2O,Cl-) = {current_params['water_cl_kij']:.6f}`, `k_ij(Na+,Cl-) = {current_params['na_cl_kij']:.6f}`.",
            f"- `feos` Held-2014 binary JSON has `water/sodium ion k_ij = {feos_params['water_na_kij']}`, `water/chloride ion k_ij = {feos_params['water_cl_kij']}`, `sodium ion/chloride ion k_ij = {feos_params['na_cl_kij']}`.",
            f"- Clapeyron advanced-like unlike CSV has `water/sodium k = {clapeyron_params['water_na_k']:.6f}`, `water/chloride k = {clapeyron_params['water_cl_k']:.6f}`, `sodium/chloride k = {clapeyron_params['na_cl_k']:.6f}`.",
            "- So after fixing the missing binary files, the large residual `feos` branch discrepancy is still present even though the key aqueous `k` values now line up with the current repo. That reinforces that the remaining issue is the exposed contribution split, not simply missing unlike parameters.",
            "",
            "## Output Artifacts",
            "",
            f"- Detailed cross-package CSV: `{BREAKDOWN_CSV}`",
            f"- Main comparison totals/plots: `{SCRIPT_DIR}`",
        ]
    )
    return "\n".join(lines) + "\n"


def main() -> None:
    if not FEOS_RAW_PATH.exists() or not CLAPEYRON_RAW_PATH.exists():
        raise FileNotFoundError(
            "Expected feos_raw.json and clapeyron_raw.json in pcsaft_package_comparisons. "
            "Run run_package_comparisons.py first."
        )
    rows = _build_breakdown_rows()
    _write_breakdown_csv(rows)
    report = _render_report(rows)
    REPORT_MD.write_text(report, encoding="utf-8")
    print(f"Wrote {BREAKDOWN_CSV}", flush=True)
    print(f"Wrote {REPORT_MD}", flush=True)


if __name__ == "__main__":
    main()
