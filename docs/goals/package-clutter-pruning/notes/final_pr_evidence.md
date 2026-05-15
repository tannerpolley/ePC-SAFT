# Final PR Evidence For Issue #120

This note is the Phase 10 evidence bundle for the package clutter pruning goal.
It summarizes the completed cleanup phases, the final validation gate, and the
scope boundary used for the PR.

## Completed Deliverables

- Phase 0 baseline: recorded branch, origin SHA, import/export surface,
  benchmark modules, optional modules, broad exception sites, generated-artifact
  scan, and baseline validation in `notes/intake.md`.
- Phase 0 prerequisite: added PyYAML to the test dependency group so
  `tests/workflows/repo -q` validates the cleanup branch.
- Phase 1 artifact scan: saved before/after scans in
  `generated_artifact_scan_before.txt` and `generated_artifact_scan_after.txt`.
  The remaining build-path matches are intentional workflow tests under
  `tests/workflows/build`.
- Phase 2 benchmark relocation: moved benchmark implementation helpers under
  `scripts/benchmarks/helpers`, kept package import shims under
  `src/epcsaft/benchmarks`, and updated workflow benchmark tests.
- Phase 3 diagnostics relocation: moved equilibrium confidence and
  thermodynamic diagnostic implementations under
  `scripts/validation/equilibrium_core`, kept package shims in
  `src/epcsaft/equilibrium_core`, and moved validation tests under
  `tests/workflows/validation/equilibrium_core`.
- Phase 4 native regression scaffold pruning: removed the unused
  `src/epcsaft/native/regression` scaffold, removed the CMake include/glob
  hooks, and kept the Ceres smoke binding tied directly to compiled Ceres
  availability.
- Phase 5 optional module cleanup: moved the IPOPT implementation under
  `src/epcsaft/_optional_backends`, kept `epcsaft.ipopt_backend` as a
  compatibility alias, and clarified sequential reactive workflow wording.
- Phase 6 ownership plan: added `docs/roadmaps/module_ownership_and_split_plan.md`
  for the large modules that remain too broad for this cleanup PR.
- Phase 7 broad exception cleanup: narrowed broad catches in `src/epcsaft` and
  documented the two retained conversion-boundary catches in
  `docs/roadmaps/exception_handling_cleanup_report.md`.
- Phase 8 wording cleanup: added `notes/ambiguous_wording_audit.md` and rewrote
  package-facing wording so optional, compatibility, and sequential routes are
  described with explicit boundaries.
- Phase 9 public API audit: added `docs/roadmaps/public_api_surface_report.md`
  and kept broad public imports intact where tests, docs, or downstream smokes
  exercise them.
- Final workflow repair: updated predefined pytest slices and benchmark metadata
  to the moved diagnostics validation paths, and added a workflow guard that
  verifies predefined pytest target files exist.

## Final Validation

- `uv run python scripts/dev/validate_project.py quick`: pass, 31 tests passed
  after doctor reported current install state and no native core errors.
- `uv run python scripts/dev/validate_project.py docs`: pass, Sphinx HTML build
  succeeded.
- `uv run python run_pytest.py tests/api/package -q`: pass, 6 tests passed.
- `uv run python run_pytest.py tests/api/runtime -q`: pass, 43 tests passed.
- `uv run python run_pytest.py tests/workflows/repo -q`: pass, 53 tests passed.
- `git diff --check`: pass, no whitespace errors.
- `git ls-files | rg "(__pycache__|\\.pyc$|\\.pyo$|\\.pyd$|\\.so$|\\.dll$|\\.dylib$|\\.egg-info/|(^|/)build/|(^|/)dist/)"`:
  pass with justified matches only in `tests/workflows/build/*`.
- Project-banned status-token scan over current diff: pass, no matches.

## Scope Boundary

This PR is a package-organization cleanup. It intentionally preserves scientific
behavior: no equation forms, solver algorithms, regression objectives,
published benchmark target values, or public runtime semantics were changed.
Compatibility shims remain where public imports already exist, and future API
shrinking is deferred to explicit deprecation or migration work.
