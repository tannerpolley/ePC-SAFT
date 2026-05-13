# Roadmap Agent Scope Contract Policy

Roadmap agents must not silently narrow their tasks.

Every implementation task must either complete the in-scope work and open a focused PR, or stop with a concrete `BLOCKED_SCOPE_GAP` handoff.

## In-Scope Item Classification

For every in-scope parameter family, derivative family, solver route, API, benchmark, or workflow row, classify it as exactly one of:

```text
implemented
already_supported_with_tests
blocker_requires_followup
out_of_scope_by_roadmap
```

Avoid vague endings such as `recorded limitation`, `unsupported`, `maybe later`, or `too broad` unless the note also records:

- exact file/function
- exact missing derivative or behavior
- exact parameter family or route
- exact future owner
- coverage/capability row preserving the boundary
- test or handoff preventing silent omission or overclaiming

## Required End States

Allowed task end states:

```text
PR_OPENED
BLOCKED_SCOPE_GAP
PREPARED_WAITING
WATCHER_TIMEOUT
MERGED
```

After implementation starts, do not stop at `PREPARED_READY`.

## k_hb_ij Example

Correct Task B handling:

- `k_hb_ij` appears in coverage/capability.
- Direct explicit algebraic association-mixing effects are represented where possible.
- Full active-association property/regression derivatives are `blocker_requires_followup` until Task C provides implicit association-site-fraction sensitivity.
- A test proves `k_hb_ij` is not silently omitted or overclaimed.

Incorrect Task B handling:

```text
k_hb_ij recorded as limitation. No PR opened.
```
