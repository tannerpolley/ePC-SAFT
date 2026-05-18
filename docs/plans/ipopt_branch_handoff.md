# Ipopt Branch Handoff

## Active Workspace

Use this checkout as authoritative:

`C:\Users\Tanner\Documents\git\ePC-SAFT`

The active local branch is:

`ipopt`

At the time this handoff was written, `ipopt` pointed at:

`197234dc Add Ipopt improvement plan`

Do not use `C:\Users\Tanner\Documents\git\ePC-SAFT-ipopt-profile` as the active workspace unless the user explicitly redirects there. That older IDE worktree was intentionally left in place and still checks out `codex/rezaee-reactive-electrolyte-lle` at the same commit.

## Branch Migration State

The local `ipopt` branch was created in `C:\Users\Tanner\Documents\git\ePC-SAFT` from `codex/rezaee-reactive-electrolyte-lle`.

The consolidated `codex/rezaee-reactive-electrolyte-lle` line already contained all local `codex/*` branches and the known remote-tracking `origin/codex/*` branches at the time of migration, so no content merge was needed.

Local stale branches deleted after creating `ipopt`:

- `codex/ipopt-profile-validation`
- `codex/native-ipopt-derivative-gates`
- `codex/prune`

Remaining local `codex/*` branch:

- `codex/rezaee-reactive-electrolyte-lle`, kept because it is checked out by `C:\Users\Tanner\Documents\git\ePC-SAFT-ipopt-profile`

Remote `origin/codex/*` branches were not deleted. Draft PR #129 was later closed because the Ipopt work is now a broader local `ipopt` overhaul lane and was not merge-ready.

## Current Dirty State To Preserve

The `ipopt` checkout had pre-existing untracked paper files at handoff time. Leave these alone unless the user explicitly asks:

- `docs/papers/md/Hubach et al. - 2024 - Li+ Extraction from Aqueous Medium Using Tetracyanoborate Ionic Liquids-Experiments and ePC-SAFT Mod.md`
- `docs/papers/md/Hubach et al. - 2024 - Supporting Information - Li+ Extraction from Aqueous Medium Using Tetracyanoborate Ionic Liquids-Experiments and ePC-SAFT Mod.md`
- `docs/papers/md/Khudaida et al. - 2026 - Supporting Information - Upgrading the Extraction of Ethanol and Isobutanol from an Aqueous Solution in the Presence of Sodiu.md`
- `docs/papers/md/Yu - 2024 - Highly efficient lithium extraction from magnesium-rich brines with ionic.md`
- `docs/papers/md/Yu - 2024 - Supplemental - Highly efficient lithium extraction from magnesium-rich brines with ionic.md`
- `docs/papers/pdf/Yu - 2024 - Supplemental - Highly efficient lithium extraction from magnesium-rich brines with ionic.pdf`

There were also separate local edits under `src/epcsaft/native/equilibrium_nlp/` when this handoff was committed. Treat those as separate work and do not revert them as part of handoff handling.

## Relevant Artifacts

Continue Ipopt improvement work from:

- `docs/plans/ipopt_improvement_plan.md`
- `analyses/paper_validation/application/2026_khudaida/scripts/diagnose_figures_2_7_live_ipopt.py`
- `analyses/paper_validation/application/2026_khudaida/diagnostics/figures_2_7_live_ipopt/`
- `analyses/paper_validation/native/2022_ascani/`

The first high-priority Ipopt improvement tranche in `docs/plans/ipopt_improvement_plan.md` is:

1. Real problem scaling.
2. Exact Hessian support.
3. Per-iteration diagnostics.
4. Warm starts with primal and multiplier state.

Recommended implementation order recorded in the plan:

1. Per-iteration diagnostics.
2. Scaling.
3. Warm starts.
4. Exact Hessian support.

## Known Ipopt Environment

Use this Windows Ipopt setup when an Ipopt-enabled proof is needed:

```powershell
$env:EPCSAFT_IPOPT_ROOT = "$env:USERPROFILE\Documents\deps\ipopt-msvc"
$env:PATH = "$env:EPCSAFT_IPOPT_ROOT\bin;$env:PATH"
$env:EPCSAFT_RUNTIME_DLL_DIRS = "$env:EPCSAFT_IPOPT_ROOT\bin"
uv run python scripts/dev/build_epcsaft.py --clean --profile ipopt --ipopt-root $env:EPCSAFT_IPOPT_ROOT --parallel 4
uv run python scripts/dev/doctor.py --require-ipopt
```

Do not treat any `Ipopt=OFF` build evidence as Ipopt proof.

## Khudaida Diagnostic Context

Figures 2-7 were the focus. Existing cached figure outputs showed good paper agreement, but live native Ipopt re-solves did not uniformly reproduce those cached splits.

Most recent 60 second live diagnostic first-row results from `diagnose_figures_2_7_live_ipopt.py`:

- Figure 2: rejected, `local_infeasibility`, residual about `0.0666`, max label `phase_equilibrium:Butanol`, organic fraction about `0.00098`.
- Figure 3: rejected, `max_iterations_exceeded`, residual about `8.83`, max label `phase_equilibrium:Butanol`.
- Figure 4: accepted, `success`, residual about `6.8e-09`.
- Figure 5: rejected, `local_infeasibility`, residual about `0.00483`, max label `phase_equilibrium:H2O`, organic fraction about `0.00095`.
- Figure 6: accepted, `acceptable_point`, residual about `4.16e-07`.
- Figure 7: accepted, `acceptable_point`, residual about `3.53e-07`.

Working interpretation: the live route is often landing near a disappearing organic phase for Figures 2 and 5. Cached paper-like splits are not enough to prove the current residual closes. This points to solver scaling, initialization, and residual parity issues rather than simply wrong digitized feed data.

## New Thread Startup

Start the new thread from:

`C:\Users\Tanner\Documents\git\ePC-SAFT`

Then verify:

```powershell
git status --short --branch
git log --oneline --decorate -5
git worktree list
```

Before implementation, read:

- project `AGENTS.md`
- `docs/.codex-journal/user_preferences.md`
- `docs/roadmaps/FULL_ROADMAP.md`
- `C:\Users\Tanner\.codex\PROJECT_ARCHITECTURE.md`
- relevant user-level Git, process, and command-context instructions under `C:\Users\Tanner\.codex\instructions\`

Suggested skills for the next session:

- `chemical-engineer`
- `diagnose` or `superpowers:systematic-debugging`
- `superpowers:test-driven-development`
- `superpowers:verification-before-completion`
- `jetbrains` when semantic C++ or Python tracing is useful

## Immediate Next Goal Candidate

If the user asks to continue the Ipopt improvement work, start with per-iteration diagnostics. That gives evidence for scaling, restoration, warm-start, and Hessian work without changing solver math first.

Keep the implementation generic for all equilibrium problems, not Khudaida-only. Public route proof and exact derivative/Jacobian requirements from `docs/roadmaps/FULL_ROADMAP.md` still apply for accepted solver claims.
