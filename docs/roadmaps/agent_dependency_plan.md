# Agent Dependency Plan
## Task matrix
| Key | Title | Branch | Phase | Dependencies | Can prep now | Can code now |
|---|---|---|---|---|---:|---:|
| A | Derivative backend completion audit and coverage matrix hard gate | `codex/backend-coverage-hard-gate` | Backend | none | yes | yes |
| B | Explicit CppAD parameter derivatives for EOS/property APIs | `codex/cppad-explicit-parameter-derivatives` | Backend | A | yes | no |
| C | Generic implicit sensitivity framework for solved states | `codex/generic-implicit-sensitivity-framework` | Backend | A, B | yes | no |
| D | General reaction and equilibrium-constant convention layer | `codex/reaction-constant-conventions` | Equilibrium | A | yes | no |
| E | Generic target-row and dataset schema | `codex/generic-target-row-schema` | Regression | A | yes | no |
| F | Generic speciation solver using ePC-SAFT activities | `codex/generic-activity-speciation` | Equilibrium | C, D | yes | no |
| G | Generic VLE/fugacity-equilibrium solver for volatile neutral species | `codex/generic-vle-fugacity-equilibrium` | Equilibrium | C | yes | no |
| H | Generic non-electrolyte LLE benchmark and solver hardening | `codex/generic-non-electrolyte-lle` | LLE | C | yes | no |
| I | Generic electrolyte LLE with distributed ions | `codex/generic-electrolyte-lle` | LLE | H, C | yes | no |
| J | Generic reactive LLE and chemical phase equilibrium | `codex/generic-reactive-lle` | Reactive | I, D, F | yes | no |
| K | Generic regression row schema and native optimizer backend | `codex/generic-regression-backend` | Regression | B, C, E | yes | no |
| L | Literature benchmark suite | `codex/literature-benchmark-suite` | Benchmarks | none | yes | no |
| M | Downstream integration smoke tests | `codex/downstream-integration-smokes` | Downstream | F, G, I, J, K | yes | no |

## Parallelization rule

Only task A starts implementation immediately.

After A merges, B, D, and E can run in parallel. C waits for B. F waits for C and D. G and H wait for C. I waits for H and C. J waits for I, D, and F. K waits for B, C, and E. L can do inventory early but must wait for implementation dependencies before adding benchmark tests. M waits for the relevant generic APIs.

## Gate rule

Project status is not the true dependency gate.

Every agent must verify:

```text
issue dependencies closed
PR dependencies merged
origin/main contains dependency merge commits
current branch is correct
branch rebases cleanly on origin/main
```

If any check fails, the agent stops.


## Branch bootstrap update

Codex app worktrees may initially be created from `main`.

Every agent must switch/create its assigned branch before Goal Prep writes files.

The assigned branch is listed in the task matrix. Use:

```text
docs/roadmaps/branch_bootstrap_instructions.md
```

If branch bootstrap fails, the agent must stop and report `BRANCH_BOOTSTRAP_FAILED`.
