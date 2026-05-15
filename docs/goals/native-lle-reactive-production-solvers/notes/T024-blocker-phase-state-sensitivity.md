# T024 Worker: Phase-State Sensitivity Blocker

## Outcome

blocked

The issue #116 accepted-solve derivative closure cannot be honestly completed in this slice.

## Evidence

The Ceres accepted electrolyte LLE solve converges for the distributed-ion diagnostic case when the Ceres cost uses `local_residual_slope`, but the accepted solve cannot yet be switched back to the T021 analytic transformed-variable Jacobian without reintroducing the distributed-ion collapse/failure.

Additional local probes showed:

- neutral phase-state sensitivities agree with projected perturbation behavior to about `1e-7`;
- the active association/ionic Ascani case shows large projected sensitivity mismatches, including order-one to order-ten errors in composition columns;
- the mismatch is upstream of the Ceres route and is consistent with the active association/ionic phase-state fugacity sensitivity path needing a source-backed derivative repair.

## Stop Condition Triggered

`The phase-state fugacity sensitivity formulas do not match source-evaluated perturbation behavior and cannot be repaired within this slice.`

## Next

The next task should be a read/derive/repair task centered on `src/epcsaft/native/epcsaft_ares.cpp` phase-state fugacity sensitivities for active association plus ionic/Born terms, with tests that compare projected solved-state sensitivities against source-evaluated perturbations for the Ascani distributed-ion fixture.
