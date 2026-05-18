# T040 Worker Receipt

## Result

Done.

## Summary

`ParameterSet` now owns canonical-to-runtime parameter compilation through `to_runtime_dict()`. Dataset loading is adapted through `ParameterSet.from_dataset(...)`, and `ePCSAFTMixture.from_dataset(...)` now constructs a `ParameterSet` before building the native mixture.

The runtime payload preserves dataset runtime options such as electrolyte model controls while keeping neutral datasets on the existing empty-charge-vector behavior. Binary hydrogen-bond interaction records now compile to the native `k_hb` runtime key rather than a non-consumed `k_hb_ij` payload key.

## Changed Files

- `src/epcsaft/parameter_schema.py`
- `src/epcsaft/epcsaft.py`
- `tests/api/parameters/test_parameter_schema.py`

## Verification

- `uv run python run_pytest.py tests/api/parameters/test_parameter_schema.py -q` -> `8 passed`
- `uv run python run_pytest.py tests/api/parameters tests/api/runtime/test_runtime_exports_and_metadata.py tests/equilibrium/core/test_api.py -q` -> `58 passed`
- `uv run python scripts/dev/check_text_gates.py` -> passed
- `git diff --check` -> passed

## Boundaries

- No dataset numerical payloads were intentionally changed.
- The old neutral template charge-vector contract remains intact.
- No downstream-specific parameter interface was added.
