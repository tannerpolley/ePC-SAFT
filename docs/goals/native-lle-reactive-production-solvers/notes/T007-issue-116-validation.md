# T007 Issue 116 Validation

## Result

Done. The issue #116 validation ladder passed after refreshing the quick-test slice name and equation registry line references.

## Validation

- `uv run python scripts/dev/build_epcsaft.py --clean --enable-ceres --enable-cppad`: pass
- `uv run python run_pytest.py tests/native/equilibrium -q`: pass, 22 tests
- `uv run python run_pytest.py tests/equilibrium/electrolyte -q`: pass, 37 passed, 7 skipped
- `uv run python run_pytest.py tests/equilibrium/core -q`: pass, 74 tests
- `uv run python run_pytest.py tests/api/equilibrium -q`: pass, 1 test
- `uv run python scripts/dev/validate_project.py quick`: pass, 31 tests
- `uv run python scripts/dev/validate_project.py docs`: pass
- `git diff --check`: pass

## Repairs During Validation

- `run_pytest.py` quick/equilibrium slices were updated from the old unavailable-derivative smoke name to the new production-derivative smoke name.
- `docs/equations_registry.yaml` and `docs/equations.md` were refreshed with `scripts/docs/sync_equation_registry.py` after native line numbers moved.

## Guard Search

Command:

```powershell
rg -n "native_derivative_free_nelder_mead|not_required.*phase split solve|newton_step.*missing sensitivity|reactive_staged_equilibrium" src tests docs
```

Remaining matches are not accepted electrolyte LLE production behavior:

- Issue #116/#117 planning text and GoalBuddy intake notes.
- Neutral LLE implementation/tests that are outside electrolyte LLE issue #116.
- #117 staged reactive APIs and tests, which remain the next gated issue.
- IPOPT seed naming and public docs for existing staged reactive APIs.
- Tests that assert old electrolyte route labels are absent from accepted diagnostics.

## Next Task

Run T008 audit for the issue #116 definition of done before unlocking issue #117 work.
