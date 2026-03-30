---
name: multi-phase-equilibrium-algorithm
description: Multiphase electrolyte equilibrium algorithm notes and equation references from Ascani 2022 plus current PC-SAFT implementation mapping. Use when working on the multiphase equilibrium formulation, Lagrangian/constraint equations, E-matrix construction, TPDF seeding, or cross-checking the current two-phase solver against the paper.
---

# Multiphase Equilibrium Algorithm

Use this skill when the task involves the Ascani 2022 multiphase electrolyte algorithm or the current `pcsaft_multiphase_lle` implementation.

## Read Order

1. `docs/multiphase_algorithm_equations.md`
2. `docs/latex/multiphase_algorithm_equations.tex`
3. `docs/multiphase_lle.md`
4. `src/pcsaft/pcsaft.pyx`

## What This Skill Covers

- Ascani 2022 equations, especially eqs. 1-35
- Independent counterion-pair construction through the $\mathbf{E}$ matrix
- Transformed-variable $\boldsymbol{\xi}$ trial-phase generation
- TPDF-based stability seeding
- Current two-phase flash residuals and their departures from the paper
- Minimum EOS equations needed to connect $\mu_i$, $\ln \varphi_i$, and $\ln f_i$

## Current Implementation Notes

- The repo solver is two-liquid-phase only, not the paper's general $\pi$-phase repeated split workflow.
- The TPDF seed follows the Ascani-style $\mathbf{E}$ and $\boldsymbol{\xi}$ construction.
- The final solve is reformulated in $\mathbf{x}^{(1)}$ and $\beta$, not the paper's all-phase transformed-variable flash equations.
- Electroneutrality is enforced during the solve with an explicit residual and then checked on both phases after convergence.
- The current implementation uses `scipy.optimize.least_squares`.

## References

- `docs/multiphase_algorithm_equations.md`: primary paper-vs-code equation map
- `docs/latex/multiphase_algorithm_equations.tex`: TeX companion in the `equations_v2.tex` style
- `docs/papers/md/Ascani, Sadowski, Held - 2022 - Calculation of Multiphase Equilibria Containing Mixed Solvents and M.md`: local paper markdown
