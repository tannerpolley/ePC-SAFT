# MEA Ion And Binary Regression Handoff For ePC-SAFT Agents

Date: 2026-05-03

This handoff is for a future Codex agent working in `C:\Users\Tanner\Documents\git\ePC-SAFT`. It supports the `MEA-Thermodynamics` paper project, which is moving to the `epcsaft` package as its only thermodynamic runtime.

The expected implementing agent may be a smaller model, such as 5.4 mini or Spark. Keep the implementation direct, test-driven, and aligned with the existing package API.

## Goal

Implement reusable ion and binary-interaction regression capability in the `epcsaft` package so `MEA-Thermodynamics` can regress parameters for the true-species `MEA-CO2-H2O` system without creating a one-off optimizer outside the package.

Current package state to verify before editing:

- `src/epcsaft/regression.py` exposes `fit_pure_neutral(...)`.
- `fit_pure_ion(...)` and `fit_binary_pair(...)` are present but currently deferred.
- `README.md` and `docs/pages/parameter_regression.rst` describe ion and binary regression as deferred.

## Required Package Behavior

Implement:

- `fit_pure_ion(...)`
- `fit_binary_pair(...)`
- `write_fit_result(...)` support for the new result shapes if needed
- tests for dataset reads/writes, bounds, multistart behavior, and user options

Use the existing public abstractions where possible:

- `FitProblem`
- `FitTerm`
- `FitBounds`
- `FitResult`
- `load_regression_records(...)`
- `get_prop_dict(...)`
- `create_parameter_template(...)`

Do not make MEA-specific assumptions in the package API. The package should support generic ion and binary-pair regression, while the MEA project supplies its own records, species, and objective choices.

## Modern Born Profile

The MEA project needs compatibility with the modern advanced Born SSM+DS route already used by `2025_Figiel` and `2026_Khudaida`.

Regression tests should include user options equivalent to:

```json
{
  "elec_model": {
    "rel_perm": {
      "rule": "empirical",
      "differential_mode": "numerical"
    },
    "born_model": {
      "d_Born_mode": 3,
      "solvation_shell_model": true,
      "dielectric_saturation": true,
      "mu_born_model": {
        "differential_mode": "numerical",
        "comp_dep_delta_d": true
      }
    }
  }
}
```

## Regression Scope

Ion regression should support fitting ion pure-parameter fields needed by the current dataset format. Start with the narrowest useful set for MEA and existing package datasets:

- `s`
- `e`
- optionally `d_born` when records and bounds request it

Binary-pair regression should support fitting matrix parameters in `mixed/binary_interaction`:

- `k_ij`
- `l_ij`
- `k_hb_ij` only if the existing runtime path can evaluate it consistently

The first implementation can be least-squares based and can reuse Python/SciPy orchestration if native support is not already present. If adding native kernels, preserve the existing Python-facing result types.

## Test Expectations

Add focused API tests before broad scientific fitting:

- `fit_pure_ion(...)` rejects missing required dataset fields with a clear `InputError`.
- `fit_pure_ion(...)` accepts a small synthetic or package-owned ion record set and returns a `FitResult`.
- `fit_binary_pair(...)` updates only the requested pair and respects symmetric matrix behavior expected by the loader.
- bounds and initial guesses are honored.
- `multistart` is deterministic enough for tests.
- advanced Born SSM+DS user options resolve and pass through objective construction.
- `write_fit_result(...)` updates only the intended CSV cells and protects existing values when `overwrite=False`.

Use small deterministic fixtures. Do not require the full MEA dataset to test generic package behavior.

## Validation Commands

Run from `C:\Users\Tanner\Documents\git\ePC-SAFT`:

```powershell
uv sync --no-install-project
uv run python scripts\build_epcsaft.py
uv run python run_pytest.py --confidence -q
```

If time is limited, run the focused regression/API tests first, then the full confidence command before finalizing.

## Coordination With MEA-Thermodynamics

After this package work exists, the MEA project should call the package APIs rather than maintaining local regression logic.

MEA target species:

- `CO2`
- `MEA`
- `H2O`
- `MEAH+`
- `MEACOO-`
- `HCO3-`
- `CO3^2-`
- `H3O+`
- `OH-`

MEA initial parameter stance:

- inherit/reference-seed: `H2O`, `CO2`, `H3O+`, `OH-`, `CO3^2-`
- regress or re-evaluate: `MEA`, `MEAH+`, `MEACOO-`, `HCO3-`
- fit binary interactions only where they materially improve VLE/speciation accuracy or stabilize the model

The MEA project handoff lives at:

`C:\Users\Tanner\Documents\git\MEA-Thermodynamics\docs\ePC-SAFT\mea-epcsaft-package-transition-handoff.md`

## Agent Notes

- Reload `AGENTS.md`, `.codex/config.toml`, and `.codex/agents/*.toml` before starting.
- Use the repo-local `epcsaft-equations` skill for equation/API/native tracing.
- Prefer owner agents for independent slices when available: `python_api_test_owner`, `native_equation_owner`, `build_packaging_owner`, and `command_runner`.
- Keep changes package-generic. The MEA paper project is the first consumer, not the only consumer.

