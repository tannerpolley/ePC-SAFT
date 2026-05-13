# Issue Scope

Source: GitHub issue #95, "Literature benchmark suite"

Current issue scope:

- Inventory and then build generic literature benchmarks for package-level confidence.
- Keep `epcsaft` general-purpose.
- Do not introduce application-specific public APIs.
- Fixture inventory can start early.
- Implementation waits on relevant solver/regression issues.
- Benchmark groups named in scope:
  - MEA simple workflow benchmark
  - MDEA ePC-SAFT benchmark
  - Figiel 2025 SSM+DS Born benchmark
  - Held 2014 revised ePC-SAFT benchmark
  - non-electrolyte LLE benchmark
  - Ascani 2022 electrolyte LLE benchmark
  - Ascani 2023 reactive LLE benchmark
  - Khudaida salting-out LLE benchmark
  - Hubach/Yu lithium-related equilibrium benchmark

Constraints recorded from the issue body:

- Dependencies are listed as `none`.
- No finite difference.
- Use CppAD for explicit algebraic derivatives.
- Use analytic formulas where exact and validated.
- Use implicit sensitivities for solved states.
- No `backend_unavailable` for required workflows.
- `backend_unavailable` is allowed only for explicitly out-of-scope workflows.
- No application-specific public APIs.

Prepared during Goal Prep so the later `/goal` run can continue from the same issue scope without rereading the issue body.

