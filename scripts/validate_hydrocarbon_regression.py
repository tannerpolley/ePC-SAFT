from __future__ import annotations

import argparse

from epcsaft import fit_pure_neutral
from tests.helpers.regression_cases import (
    _load_workbook_reference_rows,
    _neutral_fixed_parameters,
    _real_saturation_records,
)


def _print_benchmark_table(rows: list[dict[str, float]]) -> None:
    print("\nHydrocarbon neutral-component regression benchmark")
    print("Source: Gross/Sadowski (2001) Table 2 targets with NIST saturation pressure and liquid-density data")
    print(
        f"{'Species':<10} {'m_fit':>12} {'m_ref':>12} {'dm':>12} "
        f"{'s_fit':>12} {'s_ref':>12} {'ds':>12} "
        f"{'e_fit':>12} {'e_ref':>12} {'de':>12} "
        f"{'rho_rms':>12} {'vle_rms':>12}"
    )
    for row in rows:
        print(
            f"{row['species']:<10} "
            f"{row['m_fit']:>12.6f} {row['m_ref']:>12.6f} {row['m_delta']:>12.3e} "
            f"{row['s_fit']:>12.6f} {row['s_ref']:>12.6f} {row['s_delta']:>12.3e} "
            f"{row['e_fit']:>12.6f} {row['e_ref']:>12.6f} {row['e_delta']:>12.3e} "
            f"{row['rho_rms']:>12.3e} {row['vle_rms']:>12.3e}"
        )


def run_full_validation() -> list[dict[str, float]]:
    csv_rows = _load_workbook_reference_rows()
    benchmark_rows: list[dict[str, float]] = []
    for component in ("Methane", "Ethane", "Propane"):
        records = _real_saturation_records(component)
        reference = csv_rows[component]
        result = fit_pure_neutral(
            records,
            component,
            assoc_scheme="",
            fixed_parameters=_neutral_fixed_parameters(component),
            initial_guess={
                "m": reference["m"] * 1.08,
                "s": reference["s"] * 0.96,
                "e": reference["e"] * 1.05,
            },
            bounds={
                "m": (0.5, 3.5),
                "s": (2.0, 5.0),
                "e": (50.0, 400.0),
            },
        )
        if not result.success:
            raise RuntimeError(result.message)
        if result.metrics_by_term["density"] >= 0.01:
            raise RuntimeError(f"{component} density RMS exceeded threshold: {result.metrics_by_term['density']}")
        if result.metrics_by_term["pure_vle_fugacity_balance"] >= 0.01:
            raise RuntimeError(
                f"{component} pure VLE RMS exceeded threshold: {result.metrics_by_term['pure_vle_fugacity_balance']}"
            )
        benchmark_rows.append(
            {
                "species": component,
                "m_fit": result.fitted_values["m"],
                "m_ref": reference["m"],
                "m_delta": result.fitted_values["m"] - reference["m"],
                "s_fit": result.fitted_values["s"],
                "s_ref": reference["s"],
                "s_delta": result.fitted_values["s"] - reference["s"],
                "e_fit": result.fitted_values["e"],
                "e_ref": reference["e"],
                "e_delta": result.fitted_values["e"] - reference["e"],
                "rho_rms": result.metrics_by_term["density"],
                "vle_rms": result.metrics_by_term["pure_vle_fugacity_balance"],
            }
        )
    return benchmark_rows


def main() -> int:
    parser = argparse.ArgumentParser(description="Run opt-in full hydrocarbon regression reproduction.")
    parser.add_argument("--full", action="store_true", help="Run the full methane/ethane/propane reproduction.")
    args = parser.parse_args()
    if not args.full:
        parser.error("Use --full to run the opt-in long-running reproduction.")
    _print_benchmark_table(run_full_validation())
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
