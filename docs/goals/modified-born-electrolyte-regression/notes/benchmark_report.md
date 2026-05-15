# Benchmark Report - Issue #118

Scope: liquid-electrolyte modified Born / SSM / DS regression and benchmark evidence. This report does not claim vapor electrolyte Born support.

## Figiel 2025 Liquid Regression Probe

Source:

- `tests/fixtures/literature/figiel_2025/miac_liquid_electrolyte.json`
- `docs/papers/md/Figiel, Yu, Held - 2025 - Predicting Thermodynamic Properties of Ions in Single Solvents and in Mixe.md`
- `data/reference/epcsaft_parameters/2025_Figiel`
- `data/reference/MIAC/water/water-NaBr.csv`

Target rows:

- density
- relative permittivity
- osmotic coefficient
- mean ionic activity coefficient

Regression evidence from the package-owned Figiel probe:

- optimizer backend: `ceres`
- derivative/Jacobian backend: `cppad_implicit`
- python objective used: `false`
- initial parameters: `d_born=3.4`, `f_solv=1.0`
- final parameters: `d_born=3.015191151543383`, `f_solv=1.8087053108712816`
- parameter movement: `d_born=-0.384808848456617`, `f_solv=0.8087053108712816`
- active bounds: none
- objective initial: `924592.7416940106`
- objective final: `924592.7192213314`
- target-family diagnostics: density `1359.8475749608403`, relative permittivity `0.10590098480272113`, osmotic coefficient `0.008893771672830161`, mean ionic activity `0.004654763490976747`

Tolerance and interpretation:

- This is a regression-path proof, not a final literature refit. The required tolerance is objective decrease with production Ceres/CppAD evidence, nonzero parameter movement, finite diagnostics for all liquid target families, and no Python objective loop.

## Held/Cameretti 2014 Aqueous NaCl

Source:

- `data/reference/epcsaft_parameters/2014_Held`
- `data/reference/osmotic/water/NaCl.csv`
- `data/reference/MIAC/water/water-NaCl.csv`
- `analyses/paper_validation/native/2014_held`

Benchmark state:

- species: `H2O`, `Na+`, `Cl-`
- molality: `0.498905908 mol/kg`
- temperature: `298.15 K`
- pressure: `101325 Pa`

Evidence:

- density calculation: `55689.15062695025 mol m^-3`, finite and positive
- osmotic coefficient calculation: `0.9034406217587285`
- osmotic source target: `0.920710059`
- osmotic tolerance: absolute `0.03`
- mean ionic activity calculation: `0.6403185776477144`
- mean ionic activity source target: `0.681`
- mean ionic activity tolerance: absolute `0.08`

## Held 2012 Water-Methanol NaCl

Source:

- `data/reference/epcsaft_parameters/2012_Held`
- `data/reference/MIAC/water-methanol/water-methanol-NaCl.csv`
- `analyses/paper_validation/native/2012_held`

Benchmark state:

- species: `Na+`, `Cl-`, `H2O`, `Methanol`
- salt-free solvent composition: `x_H2O=0.307777668482`, `x_Methanol=0.692222331518`
- molality: `0.01244813278 mol/kg`
- temperature: `298.15 K`
- pressure: `101325 Pa`

Evidence:

- density calculation: `30937.726992380474 mol m^-3`, finite and positive
- osmotic coefficient calculation: `0.7033718483555421`, finite and bounded between 0 and 2
- mean ionic activity calculation: `0.7599493097242703`
- mean ionic activity source target: `0.812312312312`
- mean ionic activity tolerance: absolute `0.08`

## Verification Commands

- `uv run python run_pytest.py tests/regression/literature/test_figiel_held_electrolyte_benchmarks.py -q`
- `uv run python run_pytest.py tests/regression/electrolyte/test_miac_liquid_electrolyte_regression.py -q`

Both commands passed after the benchmark implementation.
