# Stage 0 Intake Gate

This note is the required issue #115 Stage 0 intake artifact. It must be completed by active task `T001` before any source edits happen.

## GitHub Source

- Issue: https://github.com/tannerpolley/ePC-SAFT/issues/115
- Title: Implement native coupled activity speciation and reactive VLE pressure benchmark
- Relevant downstream comment: https://github.com/tannerpolley/ePC-SAFT/issues/115#issuecomment-4457708894

## Current Reactive / Speciation Solver Files

To be verified by `T001` using read-only `rg` and JetBrains semantic navigation where useful.

## Current VLE / Bubble / Dew Files

To be verified by `T001` using read-only `rg` and JetBrains semantic navigation where useful.

## Current Python Wrapper Files

To be verified by `T001` using read-only `rg` and JetBrains semantic navigation where useful.

## Current Tests That Use Staged Or Fixed-Point Behavior

To be verified by `T001` using read-only test inspection and no source edits.

## Chosen Pressure / Speciation Benchmark

Default from issue #115 until `T001` verifies fixtures: generic CO2 + amine + water pressure/speciation case with repo-contained fixture data, preferably the MEA-style nine-species liquid basis with fixed literature equilibrium constants and ePC-SAFT activity evaluation.

## Intake Completion Checklist

- Current implementation files listed with exact paths.
- Current wrapper files listed with exact paths.
- Current staged or fixed-point tests listed with exact paths.
- Benchmark fixture path or fixture gap identified.
- Native-regression no-edit boundary confirmed.
- Banned-token hygiene checked for this goal directory.
