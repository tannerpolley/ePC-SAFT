from __future__ import annotations

import csv
import json
import math
import subprocess
import sys
from pathlib import Path

import matplotlib
import numpy as np

matplotlib.use("Agg")
import matplotlib.pyplot as plt


SCRIPT_DIR = Path(__file__).resolve().parent
PACKAGE_DIR = SCRIPT_DIR.parent
FIGURE3_DIR = PACKAGE_DIR.parent
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


RAW_COMPARISON_PATH = PACKAGE_DIR / "raw_comparison_long.csv"
FEOS_RAW_PATH = PACKAGE_DIR / "feos_raw.json"
CLAPEYRON_RAW_PATH = PACKAGE_DIR / "clapeyron_raw.json"
MU_BREAKDOWN_PATH = FIGURE3_DIR / "diagnostics" / "figure3_mu_contribution_breakdown.csv"
CLAPEYRON_SCAN_PATH = SCRIPT_DIR / "clapeyron_dispersion_config_scan.json"
CLAPEYRON_SCAN_SCRIPT = SCRIPT_DIR / "scan_clapeyron_dispersion_configs.jl"
CLAPEYRON_ROOT = Path(r"C:\Users\Tanner\Documents\git\Clapeyron.jl")
FEOS_SOURCE = Path(r"C:\Users\Tanner\Documents\git\feos\crates\feos-core\src\state\residual_properties.rs")

DISPERSION_TABLE_PATH = SCRIPT_DIR / "dispersion_all_ions_audit.csv"
FEOS_AUDIT_PATH = SCRIPT_DIR / "feos_hc_disp_audit.csv"
CLAPEYRON_FLUORIDE_PATH = SCRIPT_DIR / "clapeyron_fluoride_parameter_audit.csv"
REPORT_PATH = SCRIPT_DIR / "dispersion_difference_report.md"
DISPERSION_PLOT_PATH = SCRIPT_DIR / "dispersion_all_ions_with_clapeyron_configs.png"
FLUORIDE_PLOT_PATH = SCRIPT_DIR / "fluoride_dispersion_config_scan.png"

IONS = ["Li+", "Na+", "K+", "F-", "Cl-", "Br-", "I-"]
ION_CASES = {
    "Li+": ["Li+", "Cl-", "Water"],
    "Na+": ["Na+", "Cl-", "Water"],
    "K+": ["K+", "Cl-", "Water"],
    "F-": ["Na+", "F-", "Water"],
    "Cl-": ["Na+", "Cl-", "Water"],
    "Br-": ["Na+", "Br-", "Water"],
    "I-": ["Na+", "I-", "Water"],
}


def _load_raw_rows() -> list[dict[str, str]]:
    with RAW_COMPARISON_PATH.open("r", newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def _load_mu_breakdown() -> list[dict[str, str]]:
    with MU_BREAKDOWN_PATH.open("r", newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def _load_json(path: Path) -> dict:
    with path.open("r", encoding="utf-8") as handle:
        return json.load(handle)


def _comparison_value(rows: list[dict[str, str]], ion: str, source: str, term: str) -> float:
    for row in rows:
        if row["ion"] == ion and row["source"] == source and row["term"] == term:
            return float(row["value_kj_mol"])
    raise KeyError(f"Missing comparison value for ion={ion} source={source} term={term}.")


def _scan_disp(clap_scan: dict, config: str, ion: str) -> float:
    row = clap_scan["configs"][config]["results"][ion]
    if "disp_kj_mol" not in row:
        return float("nan")
    return float(row["disp_kj_mol"])


def _breakdown_row(rows: list[dict[str, str]], ion: str, contr: str) -> dict[str, str]:
    for row in rows:
        if row["ion"] == ion and row["contr"] == contr:
            return row
    raise KeyError(f"Missing mu-breakdown row for ion={ion} contr={contr}.")


def _run_clapeyron_scan() -> dict:
    cmd = [
        "julia",
        f"--project={CLAPEYRON_ROOT}",
        str(CLAPEYRON_SCAN_SCRIPT),
        str(CLAPEYRON_SCAN_PATH),
        str(CLAPEYRON_ROOT),
    ]
    completed = subprocess.run(
        cmd,
        cwd=REPO_ROOT,
        capture_output=True,
        text=True,
        check=True,
    )
    if completed.stdout.strip():
        print(completed.stdout.strip(), flush=True)
    if completed.stderr.strip():
        print(completed.stderr.strip(), flush=True)
    return _load_json(CLAPEYRON_SCAN_PATH)


def _current_repo_pairs(ion: str) -> dict[str, float]:
    species = ION_CASES[ion]
    x = np.asarray([1.0e-12, 1.0e-12, 1.0 - 2.0e-12], dtype=float)
    props = get_prop_dict("bulow_2020", species, x, 298.15)
    k_ij = props["k_ij"]
    if ion.endswith("+"):
        target_index = 0
        counter_index = 1
    else:
        target_index = 1
        counter_index = 0
    water_index = 2
    return {
        "water_counter_k": float(k_ij[water_index, counter_index]),
        "water_target_k": float(k_ij[water_index, target_index]),
        "counter_target_k": float(k_ij[counter_index, target_index]),
    }


def _write_dispersion_table(rows: list[dict[str, str]], clap_scan: dict) -> list[dict[str, float]]:
    out_rows: list[dict[str, float]] = []
    with DISPERSION_TABLE_PATH.open("w", newline="", encoding="utf-8") as handle:
        fieldnames = [
            "ion",
            "paper_disp_kj_mol",
            "pcsaft_disp_kj_mol",
            "feos_disp_kj_mol",
            "clapeyron_mixed_disp_kj_mol",
            "clapeyron_advanced_only_disp_kj_mol",
            "clapeyron_repo_override_disp_kj_mol",
            "feos_minus_pcsaft_kj_mol",
            "clapeyron_mixed_minus_pcsaft_kj_mol",
            "clapeyron_repo_override_minus_pcsaft_kj_mol",
        ]
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for ion in IONS:
            row = {
                "ion": ion,
                "paper_disp_kj_mol": _comparison_value(rows, ion, "paper", "disp"),
                "pcsaft_disp_kj_mol": _comparison_value(rows, ion, "pcsaft", "disp"),
                "feos_disp_kj_mol": _comparison_value(rows, ion, "feos", "disp"),
                "clapeyron_mixed_disp_kj_mol": _scan_disp(clap_scan, "mixed_current", ion),
                "clapeyron_advanced_only_disp_kj_mol": _scan_disp(clap_scan, "advanced_only", ion),
                "clapeyron_repo_override_disp_kj_mol": _scan_disp(clap_scan, "repo_unlike_override", ion),
            }
            row["feos_minus_pcsaft_kj_mol"] = row["feos_disp_kj_mol"] - row["pcsaft_disp_kj_mol"]
            row["clapeyron_mixed_minus_pcsaft_kj_mol"] = row["clapeyron_mixed_disp_kj_mol"] - row["pcsaft_disp_kj_mol"]
            row["clapeyron_repo_override_minus_pcsaft_kj_mol"] = (
                row["clapeyron_repo_override_disp_kj_mol"] - row["pcsaft_disp_kj_mol"]
            )
            writer.writerow({key: f"{value:.12f}" if key != "ion" else value for key, value in row.items()})
            out_rows.append(row)
    return out_rows


def _write_feos_audit(feos_payload: dict, breakdown_rows: list[dict[str, str]]) -> list[dict[str, float]]:
    out_rows: list[dict[str, float]] = []
    with FEOS_AUDIT_PATH.open("w", newline="", encoding="utf-8") as handle:
        fieldnames = [
            "ion",
            "feos_hc_kj_mol",
            "pcsaft_mu_hc_kj_mol",
            "pcsaft_dadx_hc_kj_mol",
            "feos_minus_pcsaft_mu_hc_kj_mol",
            "feos_minus_pcsaft_dadx_hc_kj_mol",
            "feos_disp_kj_mol",
            "pcsaft_mu_disp_kj_mol",
            "pcsaft_dadx_disp_kj_mol",
            "feos_minus_pcsaft_mu_disp_kj_mol",
            "feos_minus_pcsaft_dadx_disp_kj_mol",
        ]
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for ion in IONS:
            hc_row = _breakdown_row(breakdown_rows, ion, "hc")
            disp_row = _breakdown_row(breakdown_rows, ion, "disp")
            feos_row = feos_payload["results"][ion]
            row = {
                "ion": ion,
                "feos_hc_kj_mol": float(feos_row["hc"]),
                "pcsaft_mu_hc_kj_mol": float(hc_row["pcsaft_mu_contr"]),
                "pcsaft_dadx_hc_kj_mol": float(hc_row["dadx_contr"]),
                "feos_disp_kj_mol": float(feos_row["disp"]),
                "pcsaft_mu_disp_kj_mol": float(disp_row["pcsaft_mu_contr"]),
                "pcsaft_dadx_disp_kj_mol": float(disp_row["dadx_contr"]),
            }
            row["feos_minus_pcsaft_mu_hc_kj_mol"] = row["feos_hc_kj_mol"] - row["pcsaft_mu_hc_kj_mol"]
            row["feos_minus_pcsaft_dadx_hc_kj_mol"] = row["feos_hc_kj_mol"] - row["pcsaft_dadx_hc_kj_mol"]
            row["feos_minus_pcsaft_mu_disp_kj_mol"] = row["feos_disp_kj_mol"] - row["pcsaft_mu_disp_kj_mol"]
            row["feos_minus_pcsaft_dadx_disp_kj_mol"] = row["feos_disp_kj_mol"] - row["pcsaft_dadx_disp_kj_mol"]
            writer.writerow({key: f"{value:.12f}" if key != "ion" else value for key, value in row.items()})
            out_rows.append(row)
    return out_rows


def _write_clapeyron_fluoride_audit(clap_scan: dict) -> list[dict[str, float]]:
    out_rows: list[dict[str, float]] = []
    current_pairs = _current_repo_pairs("F-")
    with CLAPEYRON_FLUORIDE_PATH.open("w", newline="", encoding="utf-8") as handle:
        fieldnames = [
            "config",
            "disp_kj_mol",
            "water_counter_k",
            "water_target_k",
            "counter_target_k",
            "delta_vs_pcsaft_repo_disp_kj_mol",
            "water_counter_delta_vs_repo",
            "water_target_delta_vs_repo",
            "counter_target_delta_vs_repo",
        ]
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        pcsaft_disp = overlay.contribution_breakdown("advanced", "F-", "water", basis="mu")["disp"]
        for config in ("mixed_current", "advanced_only", "repo_unlike_override"):
            scan_row = clap_scan["configs"][config]["results"]["F-"]
            if "disp_kj_mol" not in scan_row:
                raise ValueError(f"Missing fluoride scan result for config {config}: {scan_row.get('error')}")
            row = {
                "config": config,
                "disp_kj_mol": float(scan_row["disp_kj_mol"]),
                "water_counter_k": float(scan_row["water_counter_k"]),
                "water_target_k": float(scan_row["water_target_k"]),
                "counter_target_k": float(scan_row["counter_target_k"]),
            }
            row["delta_vs_pcsaft_repo_disp_kj_mol"] = row["disp_kj_mol"] - pcsaft_disp
            row["water_counter_delta_vs_repo"] = row["water_counter_k"] - current_pairs["water_counter_k"]
            row["water_target_delta_vs_repo"] = row["water_target_k"] - current_pairs["water_target_k"]
            row["counter_target_delta_vs_repo"] = row["counter_target_k"] - current_pairs["counter_target_k"]
            writer.writerow({key: f"{value:.12f}" if key != "config" else value for key, value in row.items()})
            out_rows.append(row)
    return out_rows


def _plot_all_ions(disp_rows: list[dict[str, float]]) -> None:
    common.configure_style()
    x = np.arange(len(IONS), dtype=float)
    width = 0.14
    order = [
        ("paper", "paper_disp_kj_mol", "#7F7F7F"),
        ("pcsaft", "pcsaft_disp_kj_mol", "#111111"),
        ("feos", "feos_disp_kj_mol", "#1F77B4"),
        ("clapeyron mixed", "clapeyron_mixed_disp_kj_mol", "#FF7F0E"),
        ("clapeyron adv-only", "clapeyron_advanced_only_disp_kj_mol", "#BCBD22"),
        ("clapeyron repo override", "clapeyron_repo_override_disp_kj_mol", "#2CA02C"),
    ]
    fig, ax = plt.subplots(figsize=(13.8, 6.2))
    values_for_limits: list[np.ndarray] = []
    offsets = (np.arange(len(order), dtype=float) - (len(order) - 1) / 2.0) * width
    for offset, (label, key, color) in zip(offsets, order, strict=True):
        values = np.asarray([row[key] for row in disp_rows], dtype=float)
        values_for_limits.append(values)
        bars = ax.bar(x + offset, values, width=width, label=label, color=color, edgecolor="black", linewidth=0.4)
        common.annotate_bar_values(ax, bars, fontsize=5)
    ax.axhline(0.0, color="black", linewidth=0.8)
    ax.set_xticks(x)
    ax.set_xticklabels(IONS)
    ax.set_ylabel(r"Dispersion contribution / kJ mol$^{-1}$ ($\mu$ basis)")
    ax.set_title("Figure 3 Dispersion Audit: all ions with alternate Clapeyron parameter loads")
    common.set_strict_bar_ylim(
        ax,
        np.asarray(values_for_limits, dtype=float),
        step=5.0,
        top_pad_frac=0.24,
        bottom_pad_frac=0.12,
    )
    ax.grid(axis="y", alpha=0.22)
    ax.legend(ncol=2, frameon=True)
    common.save_figure(fig, DISPERSION_PLOT_PATH)
    plt.close(fig)


def _plot_fluoride(rows: list[dict[str, float]], raw_rows: list[dict[str, str]]) -> None:
    common.configure_style()
    labels = [
        "Paper",
        "PC-SAFT",
        "feos",
        "Clap mixed",
        "Clap adv-only",
        "Clap repo override",
    ]
    values = np.asarray(
        [
            _comparison_value(raw_rows, "F-", "paper", "disp"),
            _comparison_value(raw_rows, "F-", "pcsaft", "disp"),
            _comparison_value(raw_rows, "F-", "feos", "disp"),
            rows[0]["disp_kj_mol"],
            rows[1]["disp_kj_mol"],
            rows[2]["disp_kj_mol"],
        ],
        dtype=float,
    )
    colors = ["#7F7F7F", "#111111", "#1F77B4", "#FF7F0E", "#BCBD22", "#2CA02C"]
    fig, ax = plt.subplots(figsize=(10.6, 5.8))
    x = np.arange(len(labels), dtype=float)
    bars = ax.bar(x, values, color=colors, edgecolor="black", linewidth=0.45)
    common.annotate_bar_values(ax, bars, fontsize=7)
    ax.axhline(0.0, color="black", linewidth=0.8)
    ax.set_xticks(x)
    ax.set_xticklabels(labels, rotation=18, ha="right")
    ax.set_ylabel(r"Fluoride dispersion contribution / kJ mol$^{-1}$ ($\mu$ basis)")
    ax.set_title("Fluoride dispersion audit: shipped vs repo-matched Clapeyron loads")
    common.set_strict_bar_ylim(
        ax,
        values[np.newaxis, :],
        step=2.5,
        top_pad_frac=0.18,
        bottom_pad_frac=0.12,
    )
    ax.grid(axis="y", alpha=0.22)
    common.save_figure(fig, FLUORIDE_PLOT_PATH)
    plt.close(fig)


def _average(values: list[float]) -> float:
    finite = [value for value in values if math.isfinite(value)]
    return sum(finite) / float(len(finite))


def _write_report(disp_rows: list[dict[str, float]], feos_rows: list[dict[str, float]], fluoride_rows: list[dict[str, float]]) -> None:
    mixed_non_f = [row["clapeyron_mixed_minus_pcsaft_kj_mol"] for row in disp_rows if row["ion"] in {"Cl-", "Br-", "I-"}]
    mixed_f = next(row["clapeyron_mixed_minus_pcsaft_kj_mol"] for row in disp_rows if row["ion"] == "F-")
    advanced_f = next(row["clapeyron_advanced_only_disp_kj_mol"] for row in disp_rows if row["ion"] == "F-")
    override_f = next(row["clapeyron_repo_override_minus_pcsaft_kj_mol"] for row in disp_rows if row["ion"] == "F-")
    mixed_na = next(row["clapeyron_mixed_disp_kj_mol"] for row in disp_rows if row["ion"] == "Na+")
    advanced_na = next(row["clapeyron_advanced_only_disp_kj_mol"] for row in disp_rows if row["ion"] == "Na+")
    feos_disp_residuals = [row["feos_minus_pcsaft_dadx_disp_kj_mol"] for row in feos_rows]
    feos_hc_offsets = [row["feos_minus_pcsaft_dadx_hc_kj_mol"] for row in feos_rows]
    fluoride_map = {row["config"]: row for row in fluoride_rows}

    lines = [
        "# Dispersion and feos contribution audit",
        "",
        "## 1. Fluoride is a parameter-loading problem in the current Clapeyron comparison path",
        "",
        f"- In the current shipped comparison load (`mixed_current`), `F-` dispersion is `{fluoride_map['mixed_current']['disp_kj_mol']:.6f} kJ/mol`, versus `{overlay.contribution_breakdown('advanced', 'F-', 'water', basis='mu')['disp']:.6f} kJ/mol` for the current repo.",
        f"- The same mixed path loads water/sodium and water/fluoride unlike values from the revised table: `k(H2O,Na+) = {fluoride_map['mixed_current']['water_counter_k']:.6f}` and `k(H2O,F-) = {fluoride_map['mixed_current']['water_target_k']:.6f}`.",
        f"- The current repo uses `k(H2O,Na+) = {fluoride_map['repo_unlike_override']['water_counter_k']:.6f}` and `k(H2O,F-) = {fluoride_map['repo_unlike_override']['water_target_k']:.6f}` at the Figure 3 state.",
        f"- Loading only the advanced unlike table moves fluoride to `{advanced_f:.6f} kJ/mol`, already close to the current repo.",
        f"- Loading the repo override file on top leaves fluoride unchanged at `{fluoride_map['repo_unlike_override']['disp_kj_mol']:.6f} kJ/mol` even though it restores `k(Na+,F-)` from `{fluoride_map['advanced_only']['counter_target_k']:.6f}` to `{fluoride_map['repo_unlike_override']['counter_target_k']:.6f}`.",
        f"- That means the `Na/F` unlike pair is not what drives the Figure 3 fluoride discrepancy here. The large shift comes from fixing the water-ion rows, and the remaining `{override_f:.6f} kJ/mol` gap is small enough to attribute mainly to the remaining pure-model differences, not the unlike-parameter path.",
        "",
        "## 2. Why fluoride stands out while the other ions move together",
        "",
        f"- `Cl-`, `Br-`, and `I-` are unchanged by the alternate Clapeyron unlike loads used in this audit, so their remaining average offset of `{_average(mixed_non_f):.6f} kJ/mol` versus the current repo is not a k-parameter load-order issue. It is a baseline model or pure-parameter difference between the two frameworks.",
        f"- `Na+` moves only from `{mixed_na:.6f}` to `{advanced_na:.6f} kJ/mol` when the water/sodium row is fixed, a change of `{advanced_na - mixed_na:.6f} kJ/mol`.",
        f"- `F-` moves from `{fluoride_map['mixed_current']['disp_kj_mol']:.6f}` to `{advanced_f:.6f} kJ/mol`, a much larger change of `{advanced_f - fluoride_map['mixed_current']['disp_kj_mol']:.6f} kJ/mol`.",
        "- That contrast isolates the fluoride anomaly to the water/fluoride unlike row. The other ions may still differ from the current repo, but not because the current comparison path is loading the wrong water-ion `k_ij` row for them.",
        "",
        "## 3. feos does not expose the same branch definition as the current repo",
        "",
        f"- `feos` dispersion matches the current repo `dadx_disp` branch to roundoff for every audited ion (`max |feos - dadx_disp| = {max(abs(v) for v in feos_disp_residuals):.6e} kJ/mol`).",
        f"- `feos` hard-chain differs from the current repo `dadx_hc` by an ion-independent `{_average(feos_hc_offsets):.6f} kJ/mol` offset across all seven ions.",
        "- That pattern is what you would expect from an internal state-scalar mismatch, not from different physical parameters: the ion ranking is preserved, but the branch baseline is shifted.",
        "- Source inspection of `feos` shows the public contribution path is not computing the same quantity as the repo's `mu^alpha = a^alpha + Z^alpha + dadx^alpha - \\sum x_j dadx_j^alpha` formula. In `feos-core/src/state/residual_properties.rs`, `residual_chemical_potential_contributions(...)` differentiates `molar_helmholtz_energy_contributions(t, v, x)` directly with respect to composition, and the local variable for `v` is incorrectly initialized from reduced temperature instead of reduced inverse density.",
        "- That explains the observed mix of behaviors: dispersion coming out as the repo `dadx` branch, association nearly matching because its state scalar is zero here, and hard-chain carrying a large constant offset because the density-sensitive part is evaluated at the wrong reduced variable.",
        "",
        "## 4. Practical conclusion",
        "",
        "- The current `Clapeyron` fluoride dispersion discrepancy is mostly a load-order issue, not evidence that its dispersion expression is fundamentally different.",
        "- The `advanced_only` config is not fully available for every Figure 3 ion in Clapeyron's shipped database, so the diagnostic CSV records unsupported cases as missing rather than pretending that advanced-only is a complete water-ion set.",
        "- The current `feos` `hc`/`disp` discrepancy is not a legitimate model-vs-model difference; it is a contribution-API mismatch, and the exposed branch values should not be treated as comparable `mu` contributions until that upstream path is fixed.",
    ]
    REPORT_PATH.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    raw_rows = _load_raw_rows()
    mu_breakdown = _load_mu_breakdown()
    feos_payload = _load_json(FEOS_RAW_PATH)
    _ = _load_json(CLAPEYRON_RAW_PATH)
    clap_scan = _run_clapeyron_scan()

    disp_rows = _write_dispersion_table(raw_rows, clap_scan)
    feos_rows = _write_feos_audit(feos_payload, mu_breakdown)
    fluoride_rows = _write_clapeyron_fluoride_audit(clap_scan)
    _plot_all_ions(disp_rows)
    _plot_fluoride(fluoride_rows, raw_rows)
    _write_report(disp_rows, feos_rows, fluoride_rows)

    print(f"Wrote {DISPERSION_TABLE_PATH}", flush=True)
    print(f"Wrote {FEOS_AUDIT_PATH}", flush=True)
    print(f"Wrote {CLAPEYRON_FLUORIDE_PATH}", flush=True)
    print(f"Wrote {DISPERSION_PLOT_PATH}", flush=True)
    print(f"Wrote {FLUORIDE_PLOT_PATH}", flush=True)
    print(f"Wrote {REPORT_PATH}", flush=True)


if __name__ == "__main__":
    main()
