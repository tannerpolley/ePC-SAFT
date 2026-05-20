# T001 Roadmap Audit Receipt

Scope: audited the package against `docs/roadmaps/unified_equilibrium_core_algorithm.md` and the completion boundary in `docs/roadmaps/FULL_ROADMAP.md`.

## Coverage Matrix

| Roadmap requirement | Current package evidence | Alignment |
| --- | --- | --- |
| Public typed problems plus `mixture.equilibrium(kind=...)` facade normalize to route contracts | `src/epcsaft/equilibrium.py`, `src/epcsaft/epcsaft.py`, public problem-object exports in `src/epcsaft/__init__.py` | Mostly aligned |
| Native Ipopt production routes fail loudly when unavailable | Native bindings return `ipopt_dependency_required`; Python wrappers raise route-boundary `InputError` for gated routes | Aligned for implemented routes |
| `auto` Hessian must select exact or fail loudly | `src/epcsaft/native/equilibrium_nlp/ipopt_adapter.cpp` rejects `auto` when the NLP lacks a Hessian provider; tests cover exact and limited-memory modes | Aligned |
| Route metadata exposes `variable_model`, `density_backend`, `residual_families`, `constraint_families` | `src/epcsaft/native/equilibrium_nlp/route_metadata.h`, `src/epcsaft/equilibrium_core/native_results.py`, native-contract tests | Mostly aligned |
| Residual-family selectors activate only route-relevant rows | Neutral/electrolyte routes are mostly selected; neutral reactive LLE previously kept inactive zero phase-charge residual rows | Gap fixed in T003 |
| Phase eligibility mask is first-class generic core data | Reactive route payloads now expose `phase_eligibility_mask` and shape diagnostics; broader non-reactive route adoption remains pending | Partially aligned |
| Direct transfer and reaction residuals can coexist | Reactive LLE/electrolyte LLE residual surfaces include phase equilibrium plus reaction stationarity; phase-tagged cross-phase reaction residuals exist | Partially aligned, needs capability/benchmark proof |
| Shared stability prelayer and postsolve certification | Neutral/electrolyte stability routes exist; electrolyte LLE has TPDF postsolve certification; neutral/reactive acceptance does not use one shared stability/certification layer | Behind |
| EOS/state value contract and derivative layers | Native phase-state blocks expose pressure/fugacity/activity/density derivatives; derivative coverage, route-level Hessian evidence, and CppAD tests exist | Mostly aligned |
| Public production capability evidence is honest | `capabilities()` registers implemented Ipopt route evidence; reactive batch optimizer remains explicitly non-production | Honest, but reactive LLE/electrolyte LLE production status is not registered |
| Reactive stability | Public route is declared and gated before staged fallback | Behind: no native coupled reactive stability NLP |
| Reactive electrolyte LLE production benchmark proof | Rezaee lane exercises `mix.equilibrium(kind="reactive_electrolyte_lle")`, but current source-backed gate still fails holdout accuracy | Behind |

## Verified Gaps

1. **Fixed now: inactive phase-charge rows in neutral reactive LLE.** The neutral reactive residual surface included `phase_charge` rows even when all charges were zero. That violated the selector requirement that phase-charge rows be active only for relevant routes.
2. **Partially fixed: generic phase eligibility evidence.** Reactive route payloads now expose a flattened two-phase eligibility mask. Non-reactive routes and capability evidence still need the same shared contract if they start accepting phase-restricted species.
3. **Shared stability/certification is incomplete.** Electrolyte LLE has a hard TPDF certificate, but certification is not a common postsolve contract across neutral, reactive, and reactive electrolyte routes.
4. **Reactive stability is only a gate.** `kind="reactive_stability"` validates shape and raises at the native route boundary.
5. **Reactive LLE/electrolyte LLE are not capability-registered production routes.** Native route builders and public dispatch exist, but production capability evidence and benchmark gates are not strong enough to advertise them.
6. **Partially fixed: route-level derivative/Hessian evidence.** Production derivative coverage remains production-only, while `capabilities()["derivatives"]["equilibrium_route_evidence"]` now records exact-Hessian route-builder evidence without advertising reactive LLE as a production capability.
7. **Reactive electrolyte LLE benchmark closure is still behind.** The Rezaee package route solves source rows in the validation lane, but the published gate remains failed for accuracy/holdout behavior.

## Candidate Repair Slices

| Candidate | Allowed files | Verification |
| --- | --- | --- |
| T003 selector fix for neutral reactive phase-charge rows | `src/epcsaft/native/equilibrium/reactive_phase_equilibrium_problem.cpp`, `src/epcsaft/native/equilibrium_nlp/route_metadata.h`, focused tests | Completed |
| Add explicit phase eligibility mask to reactive route payloads and diagnostics | native reactive phase structs/bindings/tests | Native residual surface and phase-model public route tests |
| Promote route-level derivative/Hessian coverage into capability evidence | `src/epcsaft/capability_evidence.py`, `src/epcsaft/runtime.py`, runtime tests | `run_pytest.py --native-contracts -q`, runtime capability tests |
| Unify postsolve certification diagnostics before acceptance | `src/epcsaft/equilibrium.py`, native route result adapters/tests | focused neutral/electrolyte/reactive route tests |
| Build native reactive stability NLP | native equilibrium route builders, Python route gate/tests | focused native stability and reactive API tests |
| Upgrade reactive electrolyte LLE benchmark gate from failing evidence to accepted production proof | analysis fixtures, public route tests, capability evidence | Rezaee/Ascani workflow tests plus capability tests |

## T003 Repair Receipt

Changed native selector behavior so neutral reactive LLE no longer carries inactive phase-charge residual rows. Electrolyte reactive LLE still reports and enforces phase-charge rows.

Verification:

- `uv run python scripts/dev/build_epcsaft.py --build-only --parallel 10`
- `uv run python run_pytest.py tests/native/equilibrium/test_reactive_phase_equilibrium_residual_surface.py tests/native/equilibrium/test_reactive_phase_equilibrium_residual_jacobian.py tests/native/equilibrium/test_route_metadata_contracts.py tests/api/reactive/test_reactive_phase_equilibrium_problem_routes_native.py -q`
- `uv run python run_pytest.py tests/native/equilibrium/test_route_builders.py::test_reactive_lle_eos_route_builder_owns_canonical_initial_point tests/native/equilibrium/test_route_builders.py::test_reactive_lle_eos_route_uses_exact_hessian_when_requested tests/native/equilibrium/test_route_builders.py::test_reactive_electrolyte_lle_eos_route_builder_uses_liquid_root_residual_route tests/native/equilibrium/test_route_builders.py::test_reactive_electrolyte_lle_eos_route_uses_exact_hessian_when_requested -q`
- `uv run python run_pytest.py --native-contracts -q`

## T004 Repair Receipt

Added explicit phase-eligibility evidence to reactive phase residual and route-result payloads:

- `phase_eligibility_mask`: flattened 2 x species mask;
- `phase_eligibility_shape`: serialized shape for Python callers;
- diagnostics fields `phase_eligibility_mask_available`, `phase_eligibility_rows`, and `phase_eligibility_cols`.

Verification:

- `uv run python scripts/dev/build_epcsaft.py --build-only --parallel 10`
- `uv run python run_pytest.py tests/native/equilibrium/test_reactive_phase_equilibrium_residual_surface.py tests/equilibrium/reactive/test_reactive_electrolyte_lle_coupled_solver.py tests/api/reactive/test_reactive_phase_equilibrium_problem_routes_native.py -q`
- `uv run python run_pytest.py --native-contracts -q`

## T006 Repair Receipt

Added a separate equilibrium-route derivative evidence section under `epcsaft.capabilities()`:

- `coverage_matrix` remains production-only;
- `equilibrium_route_evidence` records exact Lagrangian Hessian route evidence for neutral, electrolyte, and reactive route builders;
- reactive LLE/electrolyte LLE remain explicitly capability-pending rather than advertised as implemented production equilibrium keys.

Verification:

- `uv run python run_pytest.py tests/api/runtime/test_runtime_capabilities_dependency_gates.py -q`
