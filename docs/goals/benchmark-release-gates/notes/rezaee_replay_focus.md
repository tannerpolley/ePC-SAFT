# Rezaee Replay Focus

Date: 2026-05-15

Scope: user-directed narrowed Phase 5 focus on `Lithium_Extraction` only, centered on `analyses/rezaee_2026_pcsaft_epcsaft/scripts/rezaee_reactive_equilibrium_replay.py`.

## Local install used

Installed the active issue-119 worktree into the downstream repo-local virtual environment:

```powershell
$env:EPCSAFT_PEP517_BUILD_DIR = "C:\Users\Tanner\.codex\worktrees\e3ab\ePC-SAFT\build\phase5-package-install"
uv pip install --python "C:\Users\Tanner\Documents\git\Lithium_Extraction\.venv\Scripts\python.exe" --reinstall "epcsaft @ file:///C:/Users/Tanner/.codex/worktrees/e3ab/ePC-SAFT"
```

`verified`: the install completed and replaced the prior pinned-git package with `epcsaft==1.5.2 (from file:///C:/Users/Tanner/.codex/worktrees/e3ab/ePC-SAFT)`.

## Repo-local integration check

Command:

```powershell
& "C:\Users\Tanner\Documents\git\Lithium_Extraction\.venv\Scripts\python.exe" "C:\Users\Tanner\Documents\git\Lithium_Extraction\scripts\check_epcsaft_integration.py" --mode dev
```

Observed result:

- `epcsaft module path`: `C:\Users\Tanner\Documents\git\Lithium_Extraction\.venv\Lib\site-packages\epcsaft\__init__.py`
- `epcsaft source kind`: `local_file`
- `epcsaft source detail`: `C:\Users\Tanner\.codex\worktrees\e3ab\ePC-SAFT`
- adapter smoke module loaded: `scripts.epcsaft_compat`
- failure: `Resolved epcsaft source kind 'local_file' ... is not allowed for dev mode.`

`verified`: this is a downstream contract-path blocker, not a replay-runtime failure. The downstream contract still defines `dev_worktree_root` as `C:/Users/Tanner/.codex/worktrees/epcsaft-dev/ePC-SAFT`, and that path is missing on this machine.

## Real workflow run

Command:

```powershell
& "C:\Users\Tanner\Documents\git\Lithium_Extraction\.venv\Scripts\python.exe" "C:\Users\Tanner\Documents\git\Lithium_Extraction\analyses\rezaee_2026_pcsaft_epcsaft\scripts\rezaee_reactive_equilibrium_replay.py"
```

Observed result:

- exit code `0`
- printed summary `status: "published_mismatch_but_calibrated_actual_rows_consistent"`
- `row_count`: `26`
- published Table 8/9 median `lnQ - lnK`: Li `32.33363244970333`, Na `37.758477412972525`
- calibrated paper-`K` median `lnQ - lnK`: Li `-0.23055203253576018`, Na `0.20188825902796737`
- calibrated paper-`K` median absolute complex error: `RLi = 0.0024493390000924417`, `RNa = 0.007617085730106119`
- `max_abs_charge_residual`: `8.000012776017418e-07`
- replay summary now also carries figure-validation evidence for published Figs. 7, 8, 10, and 11
- figure-replication AARD before `k_ij`: Li extraction `15.076734918338135`, selectivity `17.871733495213796`
- figure-replication AARD after `k_ij`: Li extraction `9.204927345864292`, selectivity `10.475494176274777`

Generated artifacts:

- `C:\Users\Tanner\Documents\git\Lithium_Extraction\analyses\rezaee_2026_pcsaft_epcsaft\data\processed\rezaee_2026_reactive_equilibrium_replay.csv`
- `C:\Users\Tanner\Documents\git\Lithium_Extraction\analyses\rezaee_2026_pcsaft_epcsaft\data\processed\rezaee_2026_reactive_equilibrium_paper_k_calibration.json`
- `C:\Users\Tanner\Documents\git\Lithium_Extraction\analyses\rezaee_2026_pcsaft_epcsaft\data\processed\rezaee_2026_paper_figure_digitized_points.csv`
- `C:\Users\Tanner\Documents\git\Lithium_Extraction\analyses\rezaee_2026_pcsaft_epcsaft\results\reaction_equilibrium\rezaee_2026_reactive_equilibrium_replay_summary.json`
- `C:\Users\Tanner\Documents\git\Lithium_Extraction\analyses\rezaee_2026_pcsaft_epcsaft\results\reaction_equilibrium\rezaee_2026_reactive_equilibrium_replay.md`
- `C:\Users\Tanner\Documents\git\Lithium_Extraction\analyses\rezaee_2026_pcsaft_epcsaft\results\reaction_equilibrium\rezaee_2026_paper_figure_digitization_summary.json`
- `C:\Users\Tanner\Documents\git\Lithium_Extraction\analyses\rezaee_2026_pcsaft_epcsaft\results\reaction_equilibrium\rezaee_2026_section32_paper_figures_summary.json`

## Package surfaces consumed

`verified` from script inspection:

- `from epcsaft import ePCSAFTMixture`
- `ePCSAFTMixture.from_params(...)` for both aqueous and organic mixtures
- `mixture.state(T=..., x=..., P=...)`
- `state.activity_coefficient(species=...)` for the aqueous phase
- `state.fugacity_coefficient()` for pure-component and mixture fugacity terms in the organic phase

The replay does not call an upstream packaged reactive-LLE solver. It uses package-owned state/activity/fugacity evaluation inside the downstream paper-specific reaction-equilibrium closure.

`verified` from the latest replay-script ownership pass:

- `rezaee_reactive_equilibrium_replay.py` now triggers the published-figure digitization and rendering workflow directly
- the replay summary/report now include the before/after-`k_ij` figure AARD evidence rather than leaving those figures as standalone side artifacts

## Current result and remaining blocker

`verified`: the replay runtime itself is green and now produces a clean real-data result for the 26 SI rows, but only after loading the fixed-paper-`K` calibrated organic payload generated by `rezaee_reactive_equilibrium_fit.py`.

Interpretation constrained to current evidence:

- the package can evaluate the ePC-SAFT/PC-SAFT activity terms required by the Rezaee replay;
- the published Table 8/9 organic parameter payload still does not reproduce the SI `RLi` / `RNa` complex mole fractions under the current activity-reference convention;
- the fixed-paper-`K` calibrated actual-row replay closes the SI rows cleanly enough to serve as the current real-data downstream result;
- the replay workflow now also reproduces the published visual validation trend: the digitized before-`k_ij` figures sit near the paper `15.11%` / `16.73%` AARD level and the after-`k_ij` figures improve to near the paper `7.89%` / `8.63%` level;
- the remaining blocker is now narrower: the published Table 8/9 source/reference-state gap remains unresolved, but it is not preventing a real-data replay result in the downstream workflow.

## Next safe step from this narrowed focus

Keep the downstream focus on the Rezaee replay lane:

1. treat the missing `dev_worktree_root` contract path as a downstream integration-config blocker, not a package-core blocker;
2. use the calibrated replay summary/report artifacts as the Phase 5 downstream evidence for the narrowed user-directed scope;
3. keep the published Table 8/9 mismatch explicit in reports rather than presenting the calibrated payload as the published Rezaee model;
4. if the user wants deeper follow-up here, inspect the remaining published-parameter convention gap rather than broadening back out to the other downstream repos first.
