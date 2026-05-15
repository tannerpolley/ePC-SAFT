# Benchmark Module Inventory

Phase 2 source scope from issue #120:

- `src/epcsaft/benchmarks/literature.py`
- `src/epcsaft/benchmarks/neutral_equilibrium.py`
- `src/epcsaft/benchmarks/reactive_regression.py`

## Classification

| Original runtime path | Contents | New owner | Runtime package status |
| --- | --- | --- | --- |
| `src/epcsaft/benchmarks/literature.py` | Literature benchmark inventory, subprocess-backed validation path checks, table/JSON helpers. | `scripts/benchmarks/helpers/literature.py` | Thin compatibility shim only. |
| `src/epcsaft/benchmarks/neutral_equilibrium.py` | Neutral timing benchmark cases, runners, baseline comparison, table/JSON helpers. | `scripts/benchmarks/helpers/neutral_equilibrium.py` | Thin compatibility shim only. |
| `src/epcsaft/benchmarks/reactive_regression.py` | Reactive regression/speciation benchmark cases, runners, baseline comparison, table/JSON helpers. | `scripts/benchmarks/helpers/reactive_regression.py` | Thin compatibility shim only. |

## Import Updates

- `scripts/benchmarks/benchmark_literature_suite.py` imports the helper implementation directly.
- `scripts/benchmarks/benchmark_neutral_equilibrium.py` imports the helper implementation directly.
- `scripts/benchmarks/benchmark_reactive_regression.py` imports the helper implementation directly.
- `tests/workflows/benchmarks/*` import helper implementations directly.
- `docs/pages/package_architecture.rst` no longer lists `epcsaft.benchmarks` as a runtime API surface.

## Compatibility Decision

`src/epcsaft/benchmarks` remains only as import-forwarding shims for older callers during the cleanup window. The shims contain no benchmark execution logic and point to `scripts.benchmarks.helpers`.

## Acceptance Notes

Benchmark workflows are now owned by scripts/workflow-test infrastructure rather than runtime package source. Scientific equations, benchmark targets, tolerances, and benchmark calculations were not changed; only module ownership and import paths changed.
