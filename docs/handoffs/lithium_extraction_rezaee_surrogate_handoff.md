# Lithium Extraction Rezaee Surrogate Handoff

This handoff tells the Lithium_Extraction agent how to consume the repaired
`epcsaft` package and the Rezaee 2025/2026 validation assets before building
any downstream surrogate. It is intentionally limited to package setup,
paper-validation reproduction, and pre-surrogate evidence gates.

## Required Package Version

Use the merged `tannerpolley/ePC-SAFT` `main` revision that contains this file
and the repair commit `276ad4e1` (`Repair reactive LLE validation`). Do not pin
to the old PR #126 merge commit `869e3354ddc0b52075ddc9efe687b34d6aa98316`.

From the package checkout:

```powershell
cd C:\Users\Tanner\Documents\git\ePC-SAFT
git fetch origin
git switch main
git pull --ff-only origin main
git log --oneline -- docs\handoffs\lithium_extraction_rezaee_surrogate_handoff.md
git log --oneline -- src\epcsaft\native\equilibrium\reactive_phase_equilibrium_problem.cpp
```

The first log must show this handoff on `main`; the second log must include the
reactive LLE repair history, including `276ad4e1` or its merge-equivalent.

## Build And Validate ePC-SAFT

Run from `C:\Users\Tanner\Documents\git\ePC-SAFT`:

```powershell
$env:OPENBLAS_NUM_THREADS = "1"
$env:OMP_NUM_THREADS = "1"
$env:MKL_NUM_THREADS = "1"
uv sync --no-install-project
uv run python scripts\dev\build_epcsaft.py --profile full
uv run python scripts\dev\doctor.py
uv run python run_pytest.py tests\workflows\paper_validation\test_rezaee_2026_paper_validation.py -q
```

For a wider confidence gate, also run:

```powershell
uv run python scripts\dev\validate_project.py confidence
```

Do not generate a surrogate if the focused Rezaee workflow test fails.

## Install Into Lithium_Extraction

For local development, install the same package checkout into the
Lithium_Extraction environment:

```powershell
cd C:\Users\Tanner\Documents\git\Lithium_Extraction
uv pip install -e C:\Users\Tanner\Documents\git\ePC-SAFT
uv run python -c "import epcsaft; print(epcsaft.__version__); print(epcsaft.capabilities())"
```

For a portable run, pin to the merged upstream `main` commit that contains this
handoff:

```powershell
uv pip install "epcsaft @ git+https://github.com/tannerpolley/ePC-SAFT.git@<MERGED_EPCSAFT_MAIN_COMMIT>"
```

Record the exact ePC-SAFT commit in any Lithium_Extraction surrogate metadata.

## Source Assets To Use

The package validation assets live here:

```text
C:\Users\Tanner\Documents\git\ePC-SAFT\analyses\paper_validation\application\2026_rezaee
```

To mirror the source-backed package-validation folder into Lithium_Extraction:

```powershell
robocopy C:\Users\Tanner\Documents\git\ePC-SAFT\analyses\paper_validation\application\2026_rezaee C:\Users\Tanner\Documents\git\Lithium_Extraction\analyses\rezaee_2026_pcsaft_epcsaft /E /XD __pycache__ /XF *.pyc
```

The key source-backed inputs are:

- `data\input\rezaee_2025_extraction_equilibrium_mole_fractions.csv`
- `data\input\rezaee_2025_doe_extraction_responses.csv`
- `data\input\rezaee_2025_headline_extraction_points.csv`
- `data\input\rezaee_2026_reaction_constants.csv`
- `data\input\rezaee_2026_epcsaft_species_dataset.json`
- `data\input\rezaee_2026_organic_pcsaft_parameters.csv`
- `data\input\rezaee_2026_organic_binary_interactions.csv`

The key pre-surrogate outputs are:

- `data\processed\rezaee_2025_extraction_target_summary.csv`
- `data\processed\rezaee_2025_extraction_equilibrium_summary.csv`
- `data\processed\rezaee_2026_reactive_equilibrium_replay.csv`
- `data\processed\rezaee_2026_reactive_equilibrium_fit.csv`
- `data\processed\rezaee_2026_section32_basis_inference_rows.csv`
- `data\processed\rezaee_2026_section32_equilibrium_replication_rows.csv`
- `results\reaction_equilibrium\rezaee_2026_reactive_equilibrium_replay_summary.json`
- `results\reaction_equilibrium\rezaee_2026_reactive_equilibrium_fit_summary.json`
- `results\reaction_equilibrium\rezaee_2026_section32_equilibrium_replication_summary.json`

## Reproduce The Package Validation Workflow

Run these from `C:\Users\Tanner\Documents\git\ePC-SAFT`:

```powershell
uv run python analyses\paper_validation\application\2026_rezaee\scripts\rezaee_des_epcsaft_parameter_smoke.py
uv run python analyses\paper_validation\application\2026_rezaee\scripts\rezaee_2025_target_summary.py
uv run python analyses\paper_validation\application\2026_rezaee\scripts\rezaee_reactive_equilibrium_replay.py
uv run python analyses\paper_validation\application\2026_rezaee\scripts\rezaee_reactive_equilibrium_fit.py
uv run python analyses\paper_validation\application\2026_rezaee\scripts\rezaee_reactive_convention_scan.py
uv run python analyses\paper_validation\application\2026_rezaee\scripts\rezaee_reactive_epcsaft_option_scan.py
uv run python analyses\paper_validation\application\2026_rezaee\scripts\rezaee_paper_basis_reaction_coordinate.py
uv run python analyses\paper_validation\application\2026_rezaee\scripts\rezaee_section32_basis_inference.py
uv run python analyses\paper_validation\application\2026_rezaee\scripts\rezaee_section32_equilibrium_replication.py
```

The package test
`tests\workflows\paper_validation\test_rezaee_2026_paper_validation.py` runs the
same strict source-backed workflow and checks the generated summaries.

## Supported Reactive LLE Route

The repaired package supports phase-tagged cross-phase reaction quotients
through `ReactionDefinition.phase_stoichiometry`. For Rezaee-style reactions,
build the reaction with an aggregate stoichiometry plus explicit per-phase
stoichiometry:

```python
import math
import epcsaft

reaction = epcsaft.ReactionDefinition.from_literature_constant(
    {"Li+": -1.0, "OH-": -1.0, "DES": -1.0, "RLi": 1.0, "H2O": 1.0},
    log_equilibrium_constant=math.log(3.2983e-9),
    name="Rezaee_Li_cross_phase",
    standard_state="mole_fraction_activity",
    phase_stoichiometry={
        "aq": {"Li+": -1.0, "OH-": -1.0, "H2O": 1.0},
        "org": {"DES": -1.0, "RLi": 1.0},
    },
    source="Rezaee2026_SI_Table2",
)
```

Use the analogous Na reaction from `rezaee_2026_reaction_constants.csv` when the
surrogate needs coupled Li/Na selectivity. In the native reactive LLE result,
`aq` maps to the first liquid phase and `org` maps to the second liquid phase.
Keep the species order identical across initial composition vectors, parameter
tables, and result parsing.

Expected diagnostics for the phase-tagged route include:

- `reaction_phase_scope == "phase_tagged_cross_phase"`
- `native_reaction_residual_size == 2` for the Li/Na coupled case
- `max_element_balance_norm <= 1e-10` for the source-backed replay
- `capabilities()["reactive_lle"]["supported"] is True`
- `capabilities()["reactive_lle"]["cross_phase_quotient"] == "phase_tagged_cross_phase"`

## Basis Limits To Preserve

The validation exposes a source/reference-state gap. The package can execute the
Rezaee-style cross-phase reactive LLE route, but the paper constants and local
source tables do not support claiming direct published-constant closure of the
reported extraction percentages.

Keep these labels in downstream surrogate work:

- Source-backed experimental target rows: valid for downstream target summaries.
- Replayed paper constants: diagnostic only unless the mismatch is resolved by a
  new source-backed reference-state convention.
- Bounded refit outputs: calibrated package-side evidence, not published
  constants.
- Section 3.2 replication rows: pre-surrogate diagnostic evidence with the
  source/reference-state gap stated plainly.

Do not replace this with synthetic tests, invented equilibrium rows, or default
fallback constants.

## Gate Before Surrogate Generation

Before creating surrogate inputs in Lithium_Extraction, verify and record:

- The exact merged ePC-SAFT commit.
- The focused Rezaee workflow test passes in ePC-SAFT.
- `rezaee_2026_reactive_equilibrium_replay_summary.json` reports 26 evaluated
  rows for the source-backed SI table.
- The phase-tagged package diagnostics are present.
- The 2025 extraction target summary is generated from source-backed rows.
- The source/reference-state mismatch is preserved in the surrogate notes.

Stop if any gate fails. Do not generate a surrogate from stale PR #126 code, the
old synthetic fixtures, or a locally patched Lithium_Extraction compatibility
layer that hides package failures.
