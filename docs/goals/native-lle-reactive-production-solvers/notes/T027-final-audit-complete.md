# T027 Final Audit

## Verdict

Complete. The board now satisfies the implementation and validation requirements from GitHub issues #116 and #117.

The GitHub issues were refreshed during this audit:

- #116: open, updated 2026-05-15 10:46:27 UTC.
- #117: open, updated 2026-05-15 10:46:29 UTC.

This receipt audits implementation completion only; it does not close the GitHub issues.

## Issue #116 Definition Of Done

- Accepted electrolyte LLE result is solved by Ceres trust-region residual solve: done. Native accepted diagnostics report `solver_backend = ceres` and `solver_method = ceres_trust_region_residual_solve`.
- Production residuals use explicit ePC-SAFT fugacity/chemical-potential evaluations: done. The native electrolyte residual surface and accepted Ceres route evaluate phase fugacity/chemical-potential residual blocks from ePC-SAFT phase states.
- Production Jacobian is analytic / CppAD / implicit, not manual numeric perturbation: done. Accepted diagnostics report `jacobian_backend = cppad_implicit` and `derivative_backend = cppad_implicit`; tests fail if `local_residual_slope` appears.
- Old hand-coded simplex route cannot produce an accepted production result: done. Old helpers remain only as seed/support paths, and accepted production diagnostics are guarded against the old labels.
- `newton_step` missing-sensitivity path is removed or no longer reachable for accepted production: done. Accepted production uses the Ceres residual solve and CppAD/implicit Jacobian diagnostics.
- Distributed-ion phase variables are explicit in public result and diagnostics: done.
- Each liquid phase is electroneutral: done and asserted in electrolyte benchmarks.
- Mixed-solvent and common-ion mixed-electrolyte cases are supported: done through the Ascani-style and mixed-salt tests.
- TPD/g-hat are stability and seed/acceptance tools, not the sole production solve: done.
- Ascani-style benchmark passes: done.
- Salting-out benchmark passes: done.
- Python API remains generic: done.
- All validation commands pass: done in T014 plus the additional #116 validation slices.

## Issue #117 Definition Of Done

- `ReactivePhaseEquilibriumProblem` production route is native coupled solve: done.
- Reaction and phase residuals are evaluated in one coupled solved state: done and asserted by native/API tests.
- `solver_backend` is Ceres for accepted production benchmarks: done.
- Jacobian is analytic / CppAD / implicit, not manual numeric perturbation: done, with `cppad_implicit` diagnostics.
- Staged reactive route is not the production result: done. The staged route remains explicitly named and production tests monkeypatch it to fail if called.
- #115 speciation and #116 LLE are used only as initialization/subcomponents: done.
- Neutral reactive LLE benchmark passes: done.
- Reactive electrolyte LLE benchmark passes: done.
- Phase charge balance is enforced for charged phases: done.
- Element/material balance is enforced: done.
- Reaction residual norm is below tolerance: done.
- Phase-equilibrium residual norm is below tolerance: done.
- Python API remains generic: done.
- No application-specific metric API is added: done.
- All validation commands pass: done.

## Validation Evidence

- `uv run python scripts/dev/build_epcsaft.py --clean --enable-ceres --enable-cppad`: pass.
- `uv run python run_pytest.py tests/native/equilibrium -q`: pass, 26 tests.
- `uv run python run_pytest.py tests/equilibrium/reactive -q`: pass, 4 tests.
- `uv run python run_pytest.py tests/api/reactive -q`: pass, 65 passed and 1 skipped.
- `uv run python scripts/dev/validate_project.py quick`: pass, 31 tests.
- `uv run python scripts/dev/validate_project.py docs`: pass.
- `uv run python run_pytest.py tests/equilibrium/electrolyte -q`: pass, 37 passed and 7 skipped.
- `uv run python run_pytest.py tests/equilibrium/core -q`: pass, 74 tests.
- `uv run python run_pytest.py tests/api/equilibrium -q`: pass, 1 test.
- `git diff --check`: pass.
- Issue #116 and #117 route guards ran and remaining hits were classified in `T014-validation-route-guards.md`.

## Completion State

All GoalBuddy tasks required for issues #116 and #117 are done or superseded by later repair tasks. No queued or active task remains required for the stated objective.
