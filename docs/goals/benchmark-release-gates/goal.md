# Benchmark Release Gates Goal

Source authority: GitHub issue #119, "Convert benchmark inventory and downstream smokes into executable release gates"

Issue URL: https://github.com/tannerpolley/ePC-SAFT/issues/119

## Outcome

Execute issue #119 as a validation and release-gating goal. The package must move from benchmark inventory and downstream smokes to executable literature benchmarks, numeric tolerance checks, downstream install tests, real downstream workflow runs, and release-readiness reports.

This board setup does not complete issue #119. It creates the execution board so a later `/goal` run can follow the issue phases in order.

## Hard Constraints

- Follow issue #119 exactly. If this file and the issue disagree, the GitHub issue is authoritative.
- This issue must not implement core thermodynamic algorithms. If a benchmark or downstream workflow fails because package behavior is missing, keep issue #119 open and name the implementation issue that owns the gap.
- Keep core implementation paths out of scope except for import or path fixes that do not change package behavior.
- Do not close, mark ready, or call issue #119 complete from inventories, package-local synthetic smokes, documentation-only status, staged helpers, mocked outputs, or capability labels without executable passing checks.
- Keep `epcsaft` a general-purpose package. Do not add public APIs named after downstream applications or downstream metrics.
- Use only production derivative mechanisms allowed by the issue when a benchmark depends on derivatives: analytic, CppAD, analytic implicit sensitivity, or CppAD implicit sensitivity.
- If a downstream run exposes missing production behavior, report the exact missing capability and the owning implementation issue instead of hiding the gap inside the validation issue.
- Do not commit the old forbidden backend token, old forbidden numeric-derivative token, or the old numeric-derivative phrase as contiguous text. Construct guard searches from fragments when needed.

## Owned Paths

Issue #119 may edit only the validation, benchmark, reporting, and package-install surfaces:

- `src/epcsaft/benchmarks/`
- `scripts/benchmarks/`
- `scripts/validation/`
- `tests/workflows/benchmarks/`
- `tests/api/package/`
- `tests/regression/literature/`
- `tests/fixtures/literature/`
- `analyses/`
- `docs/pages/downstream_local_installs.rst`
- `docs/roadmaps/downstream_integration_report.md`
- `docs/roadmaps/release_benchmark_report.md`

Do not edit core implementation:

- `src/epcsaft/native/`
- `src/epcsaft/equilibrium.py`
- `src/epcsaft/regression.py`
- `src/epcsaft/reactive_speciation.py`

except for import or path fixes that do not change package behavior.

## Required Phase Order

1. Phase 0 - Intake and safety baseline: confirm the current benchmark registry, executable command surfaces, tolerance checks, downstream workflow candidates, report destinations, branch state, and scope blockers before any implementation in the issue-owned paths.
2. Phase 1 - Benchmark fixture registry: ensure every required benchmark family is represented as an executable case contract rather than an inventory-only note.
3. Phase 2 - Executable benchmark commands: make the literature benchmark suite and supporting workflow commands runnable from the repo-owned benchmark surfaces.
4. Phase 3 - Numeric tolerance checks: enforce absolute and relative tolerances plus the regression and equilibrium pass/fail rules required by issue #119.
5. Phase 4 - Downstream repository install tests: prove local `epcsaft` installation and package-owned API use through package-side tests and downstream-install guidance.
6. Phase 5 - Real downstream workflow runs: run one real workflow per downstream repository and record commands, outputs, capabilities consumed, and exact blockers.
7. Phase 6 - Release-readiness reports: write the benchmark and downstream reports with pass/fail status, tolerances, capability usage, and any remaining implementation gaps tied back to owning issues.
8. Required validation - run the exact suite named in issue #119.
9. Final audit - verify every issue #119 completion line is true before calling the issue complete; otherwise leave the issue open or the PR draft with a stopped-state report.

## Required Benchmark Families

Every named benchmark must be executable or the issue is not complete:

- Gross/Sadowski pure PC-SAFT nonassociating parameters
- Gross/Sadowski associating systems
- Baygi MEA association and MEA-water binary baseline
- Cameretti/Held aqueous electrolyte density and MIAC
- Held alcohol/salt mixed-solvent density, osmotic coefficient, and MIAC
- Bulow/Ascani concentration-dependent dielectric and Born behavior
- Figiel 2025 modified Born / SSM / DS
- Ascani 2022 distributed-ion electrolyte LLE
- Ascani 2023 reactive phase equilibrium
- Khudaida 2026 salting-out LLE
- Rezaee lithium extraction thermodynamic model inputs
- MEA true-species pressure/speciation workflow fixture

## Downstream Repositories

Run one real workflow per downstream repository:

- `tannerpolley/MEA-Thermodynamics`
- `tannerpolley/Lithium_Extraction`
- `tannerpolley/MEA-Absorption-Column`

Each downstream run must prove:

- local `epcsaft` was installed
- one real downstream command ran
- generic package APIs were used
- no copied EOS code was required for package-owned behavior
- the command, output, and result were recorded

## Required Validation

Run the exact validation named in issue #119:

```powershell
uv run python scripts/benchmarks/benchmark_literature_suite.py
uv run python run_pytest.py tests/workflows/benchmarks -q
uv run python run_pytest.py tests/regression/literature -q
uv run python run_pytest.py tests/api/package -q
uv run python scripts/dev/validate_project.py quick
uv run python scripts/dev/validate_project.py docs
git diff --check
```

## Required Reports

Write:

- `docs/roadmaps/release_benchmark_report.md`
- `docs/roadmaps/downstream_integration_report.md`

Each report must include:

- commands run
- pass/fail status
- data source
- tolerances
- package capabilities consumed
- remaining non-complete implementation areas
- the implementation issue that owns each remaining area

## Completion Proof

The goal is complete only when every issue #119 completion line is true:

- benchmarks execute rather than only listing inventory
- numeric tolerances exist
- downstream repositories run real workflows
- generic package APIs are used
- no application-specific public APIs are added
- missing production behavior is reported back to implementation issues
- no core package behavior is changed in this validation issue

The final audit must answer:

- Which executable benchmark commands now exist?
- What numeric tolerance checks guard them?
- Which real downstream workflows ran and what package capabilities did they consume?
- What package gaps remain and which implementation issues own them?
- Can the package make honest release-readiness claims from executable evidence rather than inventory?

## Starter Command

`/goal Follow docs/goals/benchmark-release-gates/goal.md.`
