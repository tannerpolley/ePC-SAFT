# T002: Worker Slice Decision

Task: `T002`
Kind: `judge`
Status: `current`

## Summary

The first safe implementation slice is to add a generic package-owned literature benchmark manifest and workflow surface that classifies every issue-scope anchor as `already_supported_with_tests`, `blocker_requires_followup`, or equivalent suite status without adding any application-specific public API. This moves the issue forward immediately, documents the real coverage state, and avoids pretending missing LLE/reactive benchmark routes are complete.

## Details

- Chosen worker objective:
  - add `src/epcsaft/benchmarks/literature.py` with a machine-readable manifest for the issue-scope benchmark anchors;
  - expose the manifest through `src/epcsaft/benchmarks/__init__.py`;
  - add a small CLI surface in `scripts/benchmark_literature_suite.py`;
  - add workflow tests that verify the manifest order, classifications, and CLI JSON behavior.
- Why this slice is safe:
  - it is package-owned and generic;
  - it does not require missing solver or regression derivatives;
  - it gives the issue a truthful literature-suite backbone without silently narrowing scope.
- Deferred to follow-up slices:
  - replacing application-specific benchmark helpers with generic regression APIs;
  - adding missing MDEA, Held 2014, Ascani 2022, and Ascani 2023 runnable benchmark cases;
  - promoting Hubach/Khudaida hard cases from mixed opt-in/test-only coverage to a unified suite with stable acceptance gates.

## Board Receipt Snippet

```yaml
receipt:
  result: done
  note: notes/T002-worker-slice-decision.md
  summary: "Selected a bounded first slice: build a package-owned literature-suite manifest and CLI/tests that classify the full issue scope without adding new application-specific APIs or pretending blocked benchmark routes are implemented."
```

