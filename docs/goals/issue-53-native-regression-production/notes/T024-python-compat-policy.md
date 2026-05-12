# T024 Judge Decision: Python Compatibility Path Policy

`T024` decision after reviewing `reactive_regression.py` and the current runtime/test contracts:

- `fit_reactive_electrolyte_parameters(..., backend='python_compat')` remains intentionally supported.
- It is explicitly a compatibility/debug fallback and does not claim production readiness:
  - `diagnostics['backend'] == 'python_compat'`
  - `diagnostics['production_ready'] == False`
  - `diagnostics['Backend_unavailable_jacobian'] == True`
- Direct residual-evaluator access used by debug/pipeline diagnostics must stay available as a narrow escape hatch, but can never satisfy production contracts:
  - production capability checks continue to require native thermodynamic row solve and native production derivative semantics.
  - `issue53_native_production_ready` remains false until Ceres/CppAD derivative ownership is complete in supported rows.

No new implementation is introduced here; this is a policy-coverage decision to keep the branch honest without silently broadening production claims.
