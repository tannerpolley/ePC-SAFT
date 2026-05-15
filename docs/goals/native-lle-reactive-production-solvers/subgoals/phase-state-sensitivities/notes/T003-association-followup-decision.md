# T003 Judge: Association Follow-Up Decision

Date: 2026-05-15

## Decision

Proceed with an association-extension Worker. Parent T017 remains blocked until the phase-state sensitivity surface supports the associating electrolyte state used by the issue #116 fixtures or a different issue-valid fixture choice is documented.

## Evidence

- The T002 surface differentiates the actual pressure-based phase state for non-associating states.
- The repo-contained ionic helper uses active association for water: `assoc_num = [2, 0, 0]`.
- The same helper uses direct Born model 1 and linear dielectric mixing, so Born SSM/DS support is not required for this next fixture gate.
- Existing native code already has an association solved-state implicit helper for density response and an `association_implicit_terms_scalar_cpp(...)` recording surface that can be extended to composition variables and site fractions.

## Approved Next Worker

Implement association solved-state composition sensitivity inside the same phase-state surface by:

1. recording association outputs and residual equations with variables `[rho, x..., XA...]`;
2. solving the association residual sensitivity system for `dXA/drho` and `dXA/dx`;
3. adding association `dmu/dv` and `dzraw/dv` terms to the existing non-association phase-state chain rule;
4. keeping direct differentiation through density iteration out of scope;
5. proving the active-association fixture now returns a supported phase-state fugacity Jacobian.

Allowed files:

- `src/epcsaft/native/epcsaft_ares.cpp`
- `tests/native/cppad/test_phase_state_sensitivities.py`
- `docs/goals/native-lle-reactive-production-solvers/subgoals/phase-state-sensitivities/notes/**`

Verify:

- `uv run python scripts/dev/build_epcsaft.py --profile full`
- `uv run python run_pytest.py tests/native/cppad/test_phase_state_sensitivities.py tests/native/contracts/test_association_implicit_derivative_contract.py -q`
- `uv run python run_pytest.py tests/native/cppad tests/native/equilibrium -q`

Stop if:

- the association residual sensitivity matrix is singular for the fixture;
- the extension requires a broad association model rewrite;
- the proof would rely on a manually perturbed derivative oracle.
