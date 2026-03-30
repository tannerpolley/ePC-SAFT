# Multiphase Algorithm Equations

This note documents the equations that define the current `pcsaft_multiphase_lle` workflow and compares them to the source formulation in Ascani, Sadowski, and Held (2022), `Calculation of Multiphase Equilibria Containing Mixed Solvents and Mixed Electrolytes: General Formulation and Case Studies`.

## Sources

- Primary algorithm source: `docs/papers/md/Ascani, Sadowski, Held - 2022 - Calculation of Multiphase Equilibria Containing Mixed Solvents and M.md`
- Current solver implementation: `src/pcsaft/pcsaft.pyx`
- Current user-facing overview: `docs/multiphase_lle.md`
- Current case-study comparison script: `scripts/multiphase_model_analysis/ascani_case2_dataset_comparison.py`
- Current regression coverage: `tests/test_cython.py`
- EOS support equations: `docs/equations.md` and `docs/latex/equations_v2.tex`

## Status Legend

- `paper-identical`: same working equation as the Ascani 2022 paper
- `adapted`: same thermodynamic idea, but rearranged or reduced for the current solver
- `code-only`: introduced by the current implementation, not stated in the paper in this form
- `not implemented`: present in the paper but not carried into the current solver

## Main Paper-vs-Code Differences

- The paper presents a general $\pi$-phase liquid algorithm with repeated phase-splitting; the current repo implementation solves only a two-liquid-phase flash.
- The paper derives a transformed-variable flash in $\{n_i^{(j)}\}$ and $\{\xi_s^{(j)}\}$; the current solver uses the Ascani-style $\mathbf{E}$ matrix and $\xi$-based TPDF seeding, but reformulates the final solve in terms of one phase composition $\mathbf{x}^{(1)}$ and phase fraction $\beta$.
- The paper describes a stability-analysis plus equation-solving workflow implemented in FORTRAN with a modified Powell solver from IMSL; the current repo uses `scipy.optimize.least_squares`.
- The paper removes explicit phase-charge constraints through the transformed-variable formulation; the current solver keeps an explicit charge residual during optimization and then checks both phases after solving.
- The current repo adds practical controls that are not part of the paper derivation: `seed_x`, `force_seed_solve`, `split_tol`, `solver_accept_norm`, `charge_weight`, and `beta_bounds`.

## Coverage Map

| Paper equations | Covered here as | Current status |
| --- | --- | --- |
| eqs. (1)-(2) | `ASC-01`, `ASC-02` | paper-identical |
| eqs. (3)-(6) | `ASC-03` | paper-identical |
| eqs. (7)-(10) | `ASC-04` | paper-identical |
| eqs. (11)-(15) | `ASC-05`, `ASC-06` | paper-identical |
| eqs. (16)-(21) | `ASC-07` | paper-identical |
| eqs. (22)-(25) | `ASC-08` | adapted |
| eqs. (26)-(27) | `ASC-09` | paper-identical |
| eqs. (28)-(32) | `ASC-10`, `ASC-11` | adapted |
| eqs. (33)-(35) | `ASC-12`, `ASC-13` | adapted |

## Thermodynamic Derivation

### ASC-01

- Identifier/source: `ASC-01` / Ascani 2022 eq. (1)
- Equation:

$$
\mu_{i,\mathrm{el}}^{(j)} = \left(\frac{\partial G^{(j)}}{\partial n_i^{(j)}}\right)_{T,p,n_{k \ne i}^{(j)}} = \mu_i^{(j)} + z_i F \Phi^{(j)}
$$

- Use case in the algorithm: introduces electrostatic chemical potential and is the starting point for eliminating phase-potential differences from the working equilibrium conditions.
- Where used in code or docs: conceptual source for `docs/multiphase_lle.md`; not evaluated explicitly in `src/pcsaft/pcsaft.pyx`.
- Comparison status: `paper-identical`
- Note: the current solver never solves for $\Phi^{(j)}$ directly; it works with mean-ionic fugacity equalities instead.

### ASC-02

- Identifier/source: `ASC-02` / Ascani 2022 eq. (2)
- Equation:

$$
G = \sum_{j=1}^{\pi}\sum_{i=1}^{N_{\mathrm{neut}}} n_i^{(j)} \mu_i^{(j)} + \sum_{j=1}^{\pi}\sum_{i=1}^{N_{\mathrm{ch}}} n_i^{(j)} \left(\mu_i^{(j)} + z_i F \Phi^{(j)}\right)
$$

- Use case in the algorithm: defines the Gibbs-energy objective behind the flash problem.
- Where used in code or docs: paper derivation only; the current implementation does not minimize $G$ directly.
- Comparison status: `paper-identical`
- Note: the repo uses stability analysis plus nonlinear residual solving rather than direct Gibbs minimization.

### ASC-03

- Identifier/source: `ASC-03` / Ascani 2022 eqs. (3)-(6)
- Equation:

$$
\min_{\bar{\mathbf{n}}} G(\bar{\mathbf{n}})
$$

subject to

$$
N_i - \sum_{j=1}^{\pi} n_i^{(j)} = 0
$$

$$
\Delta_{\mathrm{ch}}^{(j)} - \sum_{i=1}^{N_{\mathrm{ch}}} z_i n_i^{(j)} = 0
$$

$$
n_i^{(j)} \ge 0
$$

- Use case in the algorithm: states the constrained equilibrium problem before the paper derives working equations.
- Where used in code or docs: conceptual source for the material-balance and charge checks in `src/pcsaft/pcsaft.pyx`.
- Comparison status: `paper-identical`
- Note: the current code keeps the mass-balance idea exactly, but handles the charge condition as a residual and post-solve validation rather than as the full paper optimization problem.

### ASC-04

- Identifier/source: `ASC-04` / Ascani 2022 eqs. (7)-(10)
- Equation:

$$
\mathcal{L}(\bar{\mathbf{n}},\bar{\mathbf{\lambda}},\bar{\mathbf{\delta}})
= \sum_{j=1}^{\pi}\sum_{i=1}^{N} n_i^{(j)}\bigl(\mu_i^{(j)} + z_i F\Phi^{(j)}\bigr)
+ \sum_{i=1}^{N} \lambda_i\left(N_i - \sum_{j=1}^{\pi} n_i^{(j)}\right)
+ \sum_{j=1}^{\pi-1} \delta_j\left(\Delta_{\mathrm{ch}}^{(j)} - \sum_{i=1}^{N_{\mathrm{ch}}} z_i n_i^{(j)}\right)
$$

with stationary conditions

$$
\mu_i^{(j)} - \lambda_i^* = 0 \quad (z_i = 0)
$$

$$
\mu_i^{(j)} - \lambda_i^* + \bigl(F\Phi^{(j)} - \delta_j^*\bigr) z_i = 0 \quad (z_i \ne 0)
$$

- Use case in the algorithm: derives the neutral and ionic equilibrium equalities used in the later flash equations.
- Where used in code or docs: derivation only; not assembled symbolically in `src/pcsaft/pcsaft.pyx`.
- Comparison status: `paper-identical`
- Note: the current solver starts from the already-reduced fugacity equalities and does not expose Lagrange multipliers.

### ASC-05

- Identifier/source: `ASC-05` / Ascani 2022 eq. (11)
- Equation:

$$
\mu_i^{(1)} = \mu_i^{(2)} = \cdots = \mu_i^{(\pi)} \qquad (z_i = 0)
$$

- Use case in the algorithm: neutral-component equilibrium condition.
- Where used in code or docs: implemented in residual form through neutral fugacity equality in `_phase_equilibrium_residual(...)`.
- Comparison status: `paper-identical`
- Note: at fixed $T$ and $p$, the code enforces equality through $\Delta \ln f_i = 0$ instead of writing $\Delta \mu_i = 0$ explicitly.

### ASC-06

- Identifier/source: `ASC-06` / Ascani 2022 eqs. (12)-(15)
- Equation:

$$
\frac{1}{z_i}\mu_i^{(j)} - \tilde{\lambda}_i^* + \tilde{\delta}_j^* = 0
$$

$$
\frac{1}{|z_{\mathrm{cat}}|}\mu_{\mathrm{cat}}^{(j)} + \frac{1}{|z_{\mathrm{an}}|}\mu_{\mathrm{an}}^{(j)}
= \text{phase-independent constant}
$$

$$
\frac{1}{|z_{\mathrm{cat}}|}\mu_{\mathrm{cat}}^{(1)} + \frac{1}{|z_{\mathrm{an}}|}\mu_{\mathrm{an}}^{(1)}
= \cdots =
\frac{1}{|z_{\mathrm{cat}}|}\mu_{\mathrm{cat}}^{(\pi)} + \frac{1}{|z_{\mathrm{an}}|}\mu_{\mathrm{an}}^{(\pi)}
$$

- Use case in the algorithm: eliminates electrostatic-potential differences and motivates mean-ionic equilibrium relations.
- Where used in code or docs: implemented in matrix form through `E.dot(delta_ch)` in `_phase_equilibrium_residual(...)`.
- Comparison status: `paper-identical`
- Note: the current code uses the log-fugacity equivalent rather than storing mean ionic chemical potentials explicitly.

### ASC-07

- Identifier/source: `ASC-07` / Ascani 2022 eqs. (16)-(21)
- Equation:

$$
\mu_{\pm}
= \frac{\dfrac{1}{|z_{\mathrm{cat}}|}\mu_{\mathrm{cat}} + \dfrac{1}{|z_{\mathrm{an}}|}\mu_{\mathrm{an}}}
{\dfrac{1}{|z_{\mathrm{cat}}|} + \dfrac{1}{|z_{\mathrm{an}}|}}
$$

$$
f_{\pm}
= \left(f_{\mathrm{cat}}^{1/|z_{\mathrm{cat}}|} f_{\mathrm{an}}^{1/|z_{\mathrm{an}}|}\right)^{
\left[\left(1/|z_{\mathrm{cat}}|\right)+\left(1/|z_{\mathrm{an}}|\right)\right]^{-1}}
$$

$$
\gamma_{\pm}
= \left(\gamma_{\mathrm{cat}}^{\nu_{\mathrm{cat}}}\gamma_{\mathrm{an}}^{\nu_{\mathrm{an}}}\right)^{1/(\nu_{\mathrm{cat}}+\nu_{\mathrm{an}})}
$$

- Use case in the algorithm: defines the mean-ionic quantities whose equality replaces direct single-ion equilibrium conditions.
- Where used in code or docs: conceptual basis for `ion_pair_rows`, `e_matrix`, and the mean-ionic checks reported by `scripts/multiphase_model_analysis/ascani_case2_dataset_comparison.py`.
- Comparison status: `paper-identical`
- Note: the current solver tracks the equivalent weighted log-fugacity difference rather than forming every mean-ionic quantity explicitly.

### ASC-08

- Identifier/source: `ASC-08` / Ascani 2022 eqs. (22)-(25)
- Equation:

$$
\mu_i^{(1)} = \mu_i^{(2)} = \cdots = \mu_i^{(\pi)} \qquad (z_i = 0)
$$

$$
\mu_{\pm,ik}^{(1)} = \mu_{\pm,ik}^{(2)} = \cdots = \mu_{\pm,ik}^{(\pi)}
$$

$$
N_i - \sum_{j=1}^{\pi} n_i^{(j)} = 0
$$

$$
\sum_{i=1}^{N_{\mathrm{ch}}} z_i n_i^{(j)} = 0
$$

- Use case in the algorithm: paper-level working flash equations.
- Where used in code or docs: adapted in `_phase_equilibrium_residual(...)`, `_residual_two_phase(...)`, and the convergence checks in `_solve_two_phase_lle(...)`.
- Comparison status: `adapted`
- Note: the current repo enforces the same neutral and mean-ionic equalities, but only for two phases and with an explicit charge residual rather than the full transformed-variable system across all phases.

## Algorithm Construction

### ASC-09

- Identifier/source: `ASC-09` / Ascani 2022 eqs. (26)-(27)
- Equation:

$$
\overline{M\mu_{\pm}} = \overline{\overline{E}} \cdot \overline{M\mu_s}
$$

with row entries

$$
e_{mn} =
\begin{cases}
1/|z_n|, & \text{if ion } n \text{ belongs to pair } m \\
0, & \text{otherwise}
\end{cases}
$$

- Use case in the algorithm: defines the independent counterion-pair matrix that maps single-ion chemical potentials into mean-ionic combinations.
- Where used in code or docs: implemented directly in `_build_e_matrix(...)`.
- Comparison status: `paper-identical`
- Note: the current code follows the paper’s high-concentration-first cation/anion pairing heuristic and checks `rank(E) = N^{\mathrm{ch}} - 1`.

### ASC-10

- Identifier/source: `ASC-10` / Ascani 2022 eqs. (28)-(30)
- Equation:

$$
\Delta n_i^{(j)} = \frac{1}{|z_i|}\xi_s^{(j)} \quad (z_i > 0),
\qquad
\Delta n_k^{(j)} = \frac{1}{|z_k|}\xi_s^{(j)} \quad (z_k < 0)
$$

$$
\bar{n}^{(j)} = \bar{n}_0^{(j)} + \bar{E}^{T}\bar{\xi}^{(j)}
$$

- Use case in the algorithm: transformed-variable idea that preserves electroneutrality while changing paired ionic mole numbers.
- Where used in code or docs: used in the TPDF seeding path through `_trial_x_from_n_xi(...)`.
- Comparison status: `adapted`
- Note: the current code uses this transformed-variable idea for trial-phase generation, but not as the final flash unknown set.

### ASC-11

- Identifier/source: `ASC-11` / Ascani 2022 eqs. (31)-(32)
- Equation:

$$
N_i - \sum_{j=1}^{\pi} n_i^{(j)} = 0 \qquad i = 1,\dots,N^{\mathrm{neut}}
$$

$$
N_k - \sum_{j=1}^{\pi} n_k^{(j)} = 0 \qquad k = 1,\dots,N^{\mathrm{ch}} - 1
$$

- Use case in the algorithm: reduced balance equations after eliminating one ionic degree of freedom per phase.
- Where used in code or docs: not solved in this exact form by the current two-phase implementation.
- Comparison status: `not implemented`
- Note: the repo bypasses this exact paper system by reconstructing the second phase from a two-phase material balance and then checking electroneutrality numerically.

### ASC-12

- Identifier/source: `ASC-12` / Ascani 2022 eq. (33)
- Equation:

$$
\operatorname{TPDF}(\mathbf{x}) = \sum_{i=1}^{N} x_i \left(\mu_i(\mathbf{x}) - \mu_i(\mathbf{z})\right)
$$

- Use case in the algorithm: phase-stability objective used to find an incipient second liquid phase.
- Where used in code or docs: implemented in `_tpdf_value(...)` through the log-fugacity equivalent.
- Comparison status: `adapted`
- Note: the code evaluates the equivalent fixed-$T,p$ form $\sum_i x_i(\ln f_i(\mathbf{x}) - \ln f_i(\mathbf{z}))$.

### ASC-13

- Identifier/source: `ASC-13` / Ascani 2022 eqs. (34)-(35)
- Equation:

$$
x_i = \frac{n_i}{\sum_{k=1}^{N^{\mathrm{neut}}} n_k + \sum_{n=1}^{N^{\mathrm{ch}}}\sum_{m=1}^{N^{\mathrm{ch}}-1} e_{nm}\xi_m}
$$

$$
x_s = \frac{\sum_{m=1}^{N^{\mathrm{ch}}-1} e_{sm}\xi_m}{\sum_{k=1}^{N^{\mathrm{neut}}} n_k + \sum_{n=1}^{N^{\mathrm{ch}}}\sum_{m=1}^{N^{\mathrm{ch}}-1} e_{nm}\xi_m}
$$

- Use case in the algorithm: converts trial neutral-component mole numbers and transformed ionic variables into a trial phase composition.
- Where used in code or docs: implemented directly in `_trial_x_from_n_xi(...)`.
- Comparison status: `paper-identical`
- Note: this is the part of the paper’s transformed-variable construction that survives verbatim in the current code.

## Current Solver Equations

### SOL-01

- Identifier/source: `SOL-01` / current code `_phase_state_liq(...)`
- Equation:

$$
\ln f_i = \ln \varphi_i + \ln x_i + \ln p
$$

- Use case in the algorithm: converts EOS fugacity-coefficient output into the log fugacity used everywhere in the flash residuals and TPDF objective.
- Where used in code or docs: `_phase_state_liq(...)`, `scripts/multiphase_model_analysis/ascani_case2_dataset_comparison.py`
- Comparison status: `code-only`
- Note: the paper stays on $\mu_i$ and $f_i$; the current implementation works directly on $\ln f_i$ because `pcsaft_lnfugcoef(...)` already returns $\ln \varphi_i$.

### SOL-02

- Identifier/source: `SOL-02` / current code `_trial_x_from_n_xi(...)`
- Equation:

$$
\mathbf{n}_{\mathrm{ch}} = \mathbf{E}^{T}\boldsymbol{\xi}, \qquad
\mathbf{x} = \frac{\bigl[\mathbf{n}_{\mathrm{neut}},\mathbf{n}_{\mathrm{ch}}\bigr]}
{\sum n_{\mathrm{neut}} + \sum n_{\mathrm{ch}}}
$$

- Use case in the algorithm: reconstructs a valid trial composition during the TPDF search.
- Where used in code or docs: `_trial_x_from_n_xi(...)`
- Comparison status: `adapted`
- Note: algebraically this is the paper construction, but the code clamps negative numerical noise with `np.maximum(..., 0.0)`.

### SOL-03

- Identifier/source: `SOL-03` / current code `_tpdf_value(...)`
- Equation:

$$
\operatorname{TPDF}_{\ln f}(\mathbf{x}) =
\sum_{i=1}^{N} x_i \left[\ln f_i(\mathbf{x}) - \ln f_i(\mathbf{z}_{\mathrm{feed}})\right]
$$

- Use case in the algorithm: stability-search objective minimized by the random global-plus-local search.
- Where used in code or docs: `_tpdf_value(...)`, `_find_tpdf_seed(...)`
- Comparison status: `adapted`
- Note: this is the fixed-$T,p$ fugacity form corresponding to the paper’s chemical-potential TPDF.

### SOL-04

- Identifier/source: `SOL-04` / current code `_residual_two_phase(...)`
- Equation:

$$
\mathbf{x}^{(2)} =
\frac{\mathbf{z}_{\mathrm{feed}} - \beta \mathbf{x}^{(1)}}{1 - \beta}
$$

- Use case in the algorithm: removes one phase composition from the unknown set by enforcing the two-phase material balance analytically.
- Where used in code or docs: `_residual_two_phase(...)`, `_solve_two_phase_lle(...)`
- Comparison status: `code-only`
- Note: this is the main place where the current solver departs from the paper’s all-phases transformed-variable flash.

### SOL-05

- Identifier/source: `SOL-05` / current code `_phase_equilibrium_residual(...)`
- Equation:

$$
\mathbf{r}_{\mathrm{neutral}} =
\Delta \ln \mathbf{f}_{\mathrm{neutral}}
= \ln \mathbf{f}_{\mathrm{neutral}}^{(1)} - \ln \mathbf{f}_{\mathrm{neutral}}^{(2)}
$$

- Use case in the algorithm: enforces neutral-component equilibrium.
- Where used in code or docs: `_phase_equilibrium_residual(...)`, `_residual_two_phase(...)`
- Comparison status: `adapted`
- Note: equivalent to the paper’s neutral $\mu$ equality at fixed $T$ and $p$.

### SOL-06

- Identifier/source: `SOL-06` / current code `_phase_equilibrium_residual(...)`
- Equation:

$$
\mathbf{r}_{\mathrm{ionic}} =
\mathbf{E}\,\Delta \ln \mathbf{f}_{\mathrm{ch}},
\qquad
\Delta \ln \mathbf{f}_{\mathrm{ch}}
= \ln \mathbf{f}_{\mathrm{ch}}^{(1)} - \ln \mathbf{f}_{\mathrm{ch}}^{(2)}
$$

- Use case in the algorithm: enforces mean-ionic equilibrium across the independent ionic pairs selected by $\mathbf{E}$.
- Where used in code or docs: `_phase_equilibrium_residual(...)`, `_residual_two_phase(...)`
- Comparison status: `adapted`
- Note: this is the direct code analogue of the paper’s mean-ionic equations, expressed on the log-fugacity basis.

### SOL-07

- Identifier/source: `SOL-07` / current code `_residual_two_phase(...)`
- Equation:

$$
r_{\mathrm{charge}} = w_{\mathrm{charge}} \, \mathbf{z}^{T}\mathbf{x}^{(1)}
$$

- Use case in the algorithm: penalizes non-electroneutral trial phase compositions during nonlinear least-squares solves.
- Where used in code or docs: `_residual_two_phase(...)`
- Comparison status: `code-only`
- Note: the paper removes explicit charge equations through transformed variables, while the current solver keeps this residual for robustness.

### SOL-08

- Identifier/source: `SOL-08` / current code `_solve_two_phase_lle(...)`
- Equation:

$$
\beta = \beta_{\mathrm{lo}} + \left(\beta_{\mathrm{hi}} - \beta_{\mathrm{lo}}\right)\sigma(y_{\beta}),
\qquad
\sigma(u) = \frac{1}{1 + e^{-u}}
$$

- Use case in the algorithm: keeps the phase fraction inside user-specified bounds during unconstrained least-squares optimization.
- Where used in code or docs: `_solve_two_phase_lle(...)`
- Comparison status: `code-only`
- Note: this bounded phase-fraction transform is a solver convenience and does not appear in the paper.

### SOL-09

- Identifier/source: `SOL-09` / current code `_solve_two_phase_lle(...)`
- Equation:

$$
x_i^{(1)} = \frac{e^{\ell_i}}{1 + \sum_{k=1}^{N-1} e^{\ell_k}} \quad (i=1,\dots,N-1),
\qquad
x_N^{(1)} = \frac{1}{1 + \sum_{k=1}^{N-1} e^{\ell_k}}
$$

- Use case in the algorithm: parameterizes the first phase composition by unconstrained logits while preserving positivity and normalization.
- Where used in code or docs: `_softmax(...)`, `_solve_two_phase_lle(...)`
- Comparison status: `code-only`
- Note: the paper never introduces a softmax parameterization because it stays in moles and transformed variables.

### SOL-10

- Identifier/source: `SOL-10` / current code `pcsaft_multiphase_lle(...)`
- Equation:

$$
\operatorname{TPDF}_{\min} \ge \mathrm{tpdf\_tol}
\;\Rightarrow\;
\text{single liquid phase}
$$

unless `force_seed_solve` and `seed_x` are both supplied.

- Use case in the algorithm: converts the stability-search result into the one-phase vs two-phase branch.
- Where used in code or docs: `pcsaft_multiphase_lle(...)`
- Comparison status: `code-only`
- Note: the paper continues with repeated split-and-restabilize logic for general multiphase systems; the current solver stops at one or two liquid phases.

## Direct EOS Dependency Equations

### EOS-01

- Identifier/source: `EOS-01` / `docs/equations.md`, `docs/latex/equations_v2.tex`
- Equation:

$$
\tilde{\mu}_k^{\mathrm{res}} = \frac{\mu_k^{\mathrm{res}}}{RT}
$$

- Use case in the algorithm: dimensionless residual chemical potential that underlies the package’s fugacity-coefficient calculations.
- Where used in code or docs: EOS documentation only; consumed indirectly through `pcsaft_lnfugcoef(...)`.
- Comparison status: `adapted`
- Note: the multiphase solver never asks for $\tilde{\mu}_k^{\mathrm{res}}$ directly, but it relies on the EOS path that converts this quantity into $\ln \varphi_k$.

### EOS-02

- Identifier/source: `EOS-02` / `docs/equations.md`, `docs/latex/equations_v2.tex`
- Equation:

$$
\ln \varphi_k = \tilde{\mu}_k^{\mathrm{res}} - \ln Z
$$

- Use case in the algorithm: bridge from residual chemical potential to fugacity coefficient.
- Where used in code or docs: EOS documentation only; consumed indirectly through `pcsaft_lnfugcoef(...)`.
- Comparison status: `adapted`
- Note: this is the minimal EOS relation needed to explain why the current flash can enforce equilibrium through $\ln f_i$.

### EOS-03

- Identifier/source: `EOS-03` / `docs/equations.md`, `docs/latex/equations_v2.tex`, current solver
- Equation:

$$
\mu_i^{(a)} = \mu_i^{(b)}
\quad \Longleftrightarrow \quad
f_i^{(a)} = f_i^{(b)}
\quad \Longleftrightarrow \quad
\ln f_i^{(a)} = \ln f_i^{(b)}
$$

- Use case in the algorithm: explains why the paper’s $\mu$-equality conditions can be enforced with log fugacities in the implementation.
- Where used in code or docs: conceptual bridge between Ascani 2022 and `src/pcsaft/pcsaft.pyx`
- Comparison status: `adapted`
- Note: this is the key thermodynamic equivalence that makes the current implementation faithful in meaning even where the algebra differs.

## Current Solver Scope Summary

- Implemented directly:
  - Ascani-style independent ionic-pair construction through $\mathbf{E}$
  - TPDF seeding on a transformed electroneutral trial space
  - neutral and mean-ionic phase-equilibrium conditions
  - two-phase material balance and electroneutrality validation
- Adapted materially:
  - final flash unknowns
  - charge handling during the solve
  - solver algorithm and bounded parameterization
- Not implemented from the paper in this form:
  - general $\pi$-phase repeated splitting loop
  - paper-style all-phases transformed-variable flash equations as the primary nonlinear solve

## Practical Reading Order

1. Read `ASC-08` through `ASC-13` to understand the paper algorithm.
2. Read `SOL-01` through `SOL-10` to see how the current repo reduces that algorithm to a practical two-phase solver.
3. Use `EOS-01` through `EOS-03` only as the minimum bridge between the paper’s chemical-potential equations and the package’s fugacity-coefficient API.
