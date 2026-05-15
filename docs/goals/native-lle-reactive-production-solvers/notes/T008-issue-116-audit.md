# T008 Issue 116 Audit

## Decision

Approved for issue #117. Issue #116 now has enough executable evidence to unlock the reactive coupled-solver work.

## Definition-of-Done Mapping

- Accepted electrolyte LLE result is solved by Ceres trust-region residual solve: done. Accepted diagnostics report `solver_backend = ceres`, `solver_method = ceres_trust_region_residual_solve`, and `acceptance_gate = ceres_residual_solve`.
- Production residuals use explicit ePC-SAFT fugacity/chemical-potential evaluations: done. Native residual rows are built from phase-state log fugacity and electroneutral salt-pair chemical-potential combinations.
- Production Jacobian is analytic / CppAD / implicit, not manual numeric perturbation: done. Accepted diagnostics report `jacobian_backend = cppad_implicit`, `derivative_backend = cppad_implicit`, `jacobian_available = true`, and the Ceres cost uses `electrolyte_residual_jacobian_row_major(...)`.
- Old hand-coded simplex route cannot produce an accepted electrolyte LLE production result: done. The electrolyte accepted route is Ceres; old derivative-free helpers remain only in neutral LLE or TPD seed-polishing roles.
- Missing-sensitivity Newton path is no longer the accepted production route: done. Public and native tests reject old accepted electrolyte labels.
- Distributed-ion variables are explicit in public result and diagnostics: done. Basis diagnostics report explicit species, charged/neutral species indices, salt-pair basis vectors, formula feed, and phase compositions.
- Each liquid phase is electroneutral: done. Benchmarks assert phase charge balance below tolerance.
- Mixed-solvent and common-ion mixed-electrolyte cases are supported: done. The mixed Na/K/Cl Ascani-style fixture passes through the public route.
- TPD/g-hat are stability and seed/acceptance tools, not sole production solve: done. Accepted solve diagnostics remain Ceres; TPD and Gibbs diagnostics are seed/stability/acceptance metadata.
- Ascani-style benchmark passes: done. `test_distributed_ion_lle_production_solver.py` passes.
- Salting-out benchmark passes: done. `test_salting_out_lle_benchmark.py` passes.
- Python API remains generic: done. `ElectrolyteLLEProblem` routes to the generic native Ceres production solver with no downstream-specific metric API.
- All validation commands pass: done. T007 recorded clean Ceres+CppAD build, focused pytest gates, quick validation, docs validation, diff check, and justified guard search.

## Residual Guard Matches

Remaining route-guard matches are not accepted electrolyte LLE production behavior:

- Neutral LLE still has its historical derivative-free route and tests.
- TPD polishing remains seed/stability machinery.
- IPOPT still names a historical native seed source.
- #117 staged reactive APIs and tests remain intentionally queued for the next issue.
- Goal notes and plans preserve historical intake evidence.

## Next Task

Start T009, the issue #117 route map, with #116 treated as available in this branch.
