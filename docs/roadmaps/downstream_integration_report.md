# Downstream Integration Report

Status: in progress, not complete.

Issue scope: GitHub issue `#119`, "Convert benchmark inventory and downstream smokes into executable release gates".

Current report state:

- This report records the user-directed narrowed downstream focus on the `Lithium_Extraction` Rezaee replay lane.
- `MEA-Thermodynamics` and `MEA-Absorption-Column` remain pending for full issue `#119` completion.
- No package-core implementation gap has been proven from the Rezaee replay run itself.
- The currently observed blockers are a downstream integration-contract path mismatch and an unresolved published Table 8/9 source/reference-state gap.

## Summary

| Downstream repo | Local install proof | Repo-local contract check | Real workflow run | Current result | Package-core gap identified |
| --- | --- | --- | --- | --- | --- |
| `Lithium_Extraction` | complete | blocked by configured `dev_worktree_root` mismatch | complete | calibrated paper-`K` Rezaee replay is consistent on 26 SI rows; published Table 8/9 payload still mismatches | none proven from current evidence |
| `MEA-Thermodynamics` | pending in this narrowed slice | pending in this narrowed slice | pending in this narrowed slice | pending | pending |
| `MEA-Absorption-Column` | pending in this narrowed slice | pending in this narrowed slice | pending in this narrowed slice | pending | pending |

## Lithium_Extraction

### Local install used

Command run from the active `ePC-SAFT` worktree:

```powershell
$env:EPCSAFT_PEP517_BUILD_DIR = "C:\Users\Tanner\.codex\worktrees\e3ab\ePC-SAFT\build\phase5-package-install"
uv pip install --python "C:\Users\Tanner\Documents\git\Lithium_Extraction\.venv\Scripts\python.exe" --reinstall "epcsaft @ file:///C:/Users/Tanner/.codex/worktrees/e3ab/ePC-SAFT"
```

Observed result:

- install completed successfully
- downstream environment now resolves `epcsaft==1.5.2` from `file:///C:/Users/Tanner/.codex/worktrees/e3ab/ePC-SAFT`

### Repo-local integration check

Command:

```powershell
& "C:\Users\Tanner\Documents\git\Lithium_Extraction\.venv\Scripts\python.exe" "C:\Users\Tanner\Documents\git\Lithium_Extraction\scripts\check_epcsaft_integration.py" --mode dev
```

Observed result:

- `epcsaft source kind`: `local_file`
- `epcsaft source detail`: `C:\Users\Tanner\.codex\worktrees\e3ab\ePC-SAFT`
- adapter smoke module loaded: `scripts.epcsaft_compat`
- failure reason: repo contract only allows `live_worktree` for `--mode dev`

Interpretation:

- this is a downstream integration-contract configuration blocker
- the downstream contract still points at `C:/Users/Tanner/.codex/worktrees/epcsaft-dev/ePC-SAFT`
- that configured path is missing on this machine
- this is not evidence of missing package runtime behavior in the active `e3ab` worktree

### Real workflow run

Command:

```powershell
& "C:\Users\Tanner\Documents\git\Lithium_Extraction\.venv\Scripts\python.exe" "C:\Users\Tanner\Documents\git\Lithium_Extraction\analyses\rezaee_2026_pcsaft_epcsaft\scripts\rezaee_reactive_equilibrium_replay.py"
```

Observed result:

- exit code `0`
- printed `status: "published_mismatch_but_calibrated_actual_rows_consistent"`
- `row_count`: `26`
- published Table 8/9 median `lnQ - lnK`: Li `32.33363244970333`, Na `37.758477412972525`
- calibrated paper-`K` median `lnQ - lnK`: Li `-0.23055203253576018`, Na `0.20188825902796737`
- calibrated paper-`K` median absolute complex error: `RLi = 0.0024493390000924417`, `RNa = 0.007617085730106119`
- `max_abs_charge_residual`: `8.000012776017418e-07`
- replay-owned figure-validation AARD before `k_ij`: Li extraction `15.076734918338135`, selectivity `17.871733495213796`
- replay-owned figure-validation AARD after `k_ij`: Li extraction `9.204927345864292`, selectivity `10.475494176274777`

Generated downstream artifacts:

- `C:\Users\Tanner\Documents\git\Lithium_Extraction\analyses\rezaee_2026_pcsaft_epcsaft\data\processed\rezaee_2026_reactive_equilibrium_replay.csv`
- `C:\Users\Tanner\Documents\git\Lithium_Extraction\analyses\rezaee_2026_pcsaft_epcsaft\data\processed\rezaee_2026_reactive_equilibrium_paper_k_calibration.json`
- `C:\Users\Tanner\Documents\git\Lithium_Extraction\analyses\rezaee_2026_pcsaft_epcsaft\data\processed\rezaee_2026_paper_figure_digitized_points.csv`
- `C:\Users\Tanner\Documents\git\Lithium_Extraction\analyses\rezaee_2026_pcsaft_epcsaft\results\reaction_equilibrium\rezaee_2026_reactive_equilibrium_replay_summary.json`
- `C:\Users\Tanner\Documents\git\Lithium_Extraction\analyses\rezaee_2026_pcsaft_epcsaft\results\reaction_equilibrium\rezaee_2026_reactive_equilibrium_replay.md`
- `C:\Users\Tanner\Documents\git\Lithium_Extraction\analyses\rezaee_2026_pcsaft_epcsaft\results\reaction_equilibrium\rezaee_2026_paper_figure_digitization_summary.json`
- `C:\Users\Tanner\Documents\git\Lithium_Extraction\analyses\rezaee_2026_pcsaft_epcsaft\results\reaction_equilibrium\rezaee_2026_section32_paper_figures_summary.json`

### Generic package surfaces consumed

Verified from downstream script inspection:

- `epcsaft.ePCSAFTMixture.from_params(...)` for aqueous and organic mixtures
- `mixture.state(T=..., x=..., P=...)`
- `state.activity_coefficient(species=...)` for aqueous activity terms
- `state.fugacity_coefficient()` for pure-component and mixture fugacity terms in the organic phase

The replay does not use a package-owned reactive-LLE production solver. It uses generic package state/activity/fugacity calls inside the downstream paper-specific reaction-equilibrium closure.

The replay script now also owns the published-figure validation pass, so the Rezaee downstream workflow produces both real SI-row replay evidence and replicated Figs. 7, 8, 10, and 11 from the same top-level command.

### Current result and remaining blocker

The current blocker is not a failed import, install, or missing package call.

Current evidence supports this narrower, updated conclusion:

- the package can evaluate the ePC-SAFT/PC-SAFT activity terms needed by the Rezaee replay
- the published Table 2 constants together with the published Table 8/9 organic parameter payload do not reproduce the SI `RLi` / `RNa` complex mole fractions under the current source/reference-state convention
- a fixed-paper-`K` calibrated organic payload now yields a clean real-data replay result on the 26 SI rows without changing the aqueous ePC-SAFT calls or the paper equilibrium constants
- the replay-owned figure validation preserves the paper's before/after-`k_ij` improvement trend, with the no-`k_ij` digitized AARDs near the published `15.11%` / `16.73%` values and the with-`k_ij` digitized AARDs improved to near the published `7.89%` / `8.63%` values
- the remaining blocker is therefore the unresolved published Table 8/9 source/reference-state gap, not an inability to produce a real downstream replay result

Package-gap routing status:

- no package-owned missing capability has been identified from this run alone
- no upstream implementation issue can be named honestly yet from this evidence
- keep the published-parameter mismatch explicit until the source/reference-state gap is either traced to a package capability gap or ruled downstream-only

## Pending downstream lanes

The following required issue `#119` downstream proofs remain open:

- `MEA-Thermodynamics`: one real workflow run after local `epcsaft` install
- `MEA-Absorption-Column`: one real workflow run after local `epcsaft` install

Because the current user-directed focus was narrowed to the Rezaee replay lane, those repos were not advanced further in this report update.
