# Issue #91 Current Scope

Issue: https://github.com/tannerpolley/ePC-SAFT/issues/91

Current scope from the live issue:

- prove ordinary two-liquid-phase splitting before layering electrolyte constraints
- keep `epcsaft` general-purpose
- stay within the generic LLE benchmark and solver hardening tranche

In scope:

- phase split
- fugacity equality
- stability checks or anti-trivial-solution strategy
- clear phase diagnostics
- simple literature or repo benchmark

Out of scope:

- electrolyte accounting
- forcing a poor benchmark to pass
- application-specific public APIs
- finite difference

Dependency status:

- dependency C is closed
- PR #104 is merged
- commit `974e6a232025d9d305297590228c0d2131fbc4fe` is in `origin/main`
