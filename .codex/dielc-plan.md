**Goal**  
- Replace the scalar dielectric constant input with per-component dielectric constants and compute the mixture dielectric constant in C++ using the mole-fraction linear rule ($\varepsilon_r = \sum_j \varepsilon_{r,j} x_j$) from the Dielectric Constant section of `docs/equations.tex` (rule 1).

**Equation Basis**  
- Mixing rule: $\varepsilon_r(\mathbf x) = \sum_j \varepsilon_{r,j}\,x_j$ (mole-fraction rule in `docs/equations.tex`, Dielectric Constant subsection).  
- Derivative (mole-fraction rule): $\left(\partial \varepsilon_r / \partial x_i\right)_{x_{j\ne i}} = \varepsilon_{r,i}$.

**Interface Changes (Python/Cython/C++)**  
- `pcsaft.pxd` / `add_args` struct: change `dielc` from `double` to `vector<double>` to hold per-component $\varepsilon_r$.  
- `pcsaft.pyx::create_struct`: accept `params['dielc']` as an array; require length == ncomp; convert with `np_to_vector_double`; raise `InputError` otherwise.  
- Docstrings: update `dielc` parameter to “ndarray, shape (n,) dielectric constant per component; mixture $\varepsilon_r$ computed internally by mole-fraction rule.”  
- Keep `dielc_water(t)` helper; users can still pass identical values for all components if desired.

**Mixture Dielectric Evaluation (C++)**  
- Add helper `epsilon_mix(const vector<double>& x, const vector<double>& eps_r)` computing $\varepsilon_r = \sum_j \varepsilon_{r,j} x_j$.  
- Replace all uses of `cppargs.dielc` in electrolyte terms with a locally computed `eps_mix`:  
  - Debye–Hückel contributions: $\kappa$ and downstream uses in `pcsaft_*` functions (`pcsaft_den_cpp`, `pcsaft_fugcoef_cpp`, `pcsaft_ares_cpp`, `pcsaft_dadt_cpp`, etc.) where divisions by `cppargs.dielc` occur.  
  - Flash routines (`flashPQ_cpp`, `flashTQ_cpp`, `outerPQ`, `outerTQ`): remove overriding `cppargs.dielc = dielc_water(t)`; compute `eps_mix` from current phase composition before density/fugacity calls.  
  - `get_single_component`: carry the component’s dielectric constant when forming a single-component struct.

**Python-Level Validation/Behavior**  
- Require per-component dielectric arrays only; no scalar broadcast allowed.

**Tests to Update**  
- `tests/test_cython.py`: electrolyte cases currently set `dielc = dielc_water(t)` (scalar). Update to pass an array of length n (e.g., `np.full_like(m, dielc_water(t))` or explicit per-component values) so validation succeeds.  
- Sphinx docs likely unchanged (`docs/functions/pcsaft.dielc_water.rst`).

**Non-Goals / Deferred**  
- Born term changes remain out of scope.

**Decisions (per user)**  
- Require per-component dielectric input arrays only (no scalar broadcast).  
- Implement concentration-dependent $\varepsilon_r$ derivatives in Debye–Hückel chemical potentials using the “updated version for concentration dependent dielectric constant” (Bulow 2019) equations in `docs/equations.tex` (exclude Bjerrum treatment).

**Additional Scope for Derivatives (Bulow 2019 section)**  
- Use the provided expressions (from `docs/equations.tex`, “The updated version for concentration dependent dielectric constant”):  
  - $\kappa(\mathbf x) = \sqrt{ \dfrac{\rho_N e^2}{k_B T\,\varepsilon_0\,\varepsilon_r} \sum_j x_j z_j^2 }$  
  - $\left( \dfrac{\partial \kappa}{\partial x_i} \right)_{T,v,x_{j\ne i}}$ as given in the Bulow 2019 subsection (retain exact form).  
  - $\left( \dfrac{\partial \varepsilon_r}{\partial x_i} \right)_{x_{j\ne i}} = \varepsilon_{r,i}$  
  - $\left( \dfrac{\partial \chi_i}{\partial x_i} \right)_{T,v,x_{j\ne i}}$ using the expression with $\partial \kappa/\partial x_i$.  
  - $\left( \dfrac{\partial \tilde a^{DH}}{\partial x_i} \right)_{T,v,x_{j\ne i}}$ four-term bracket from Bulow 2019 (non-Bjerrum).  
- Implement these derivatives in C++ for $\mu^{DH}$ / $\ln \gamma$ to align with composition-dependent $\varepsilon_r$; keep Bjerrum variants disabled/not used.

**Note: 2019 vs. 2005 μ\_DH Equivalence (when ∂ε/∂x = 0)**
- Starting from Bulow 2019: $\mu^{DH}\_k = a^{DH} + Z^{DH} + \left(\partial a^{DH}/\partial x_k\right) - \sum_j x_j (\partial a^{DH}/\partial x_j)$ with the composition derivative terms driven by $\kappa(x)$ and $\chi(\kappa)$.
- If $\partial \varepsilon_r/\partial x_i = 0$ and $\kappa$ is composition-independent, $\partial a^{DH}/\partial x_k \to 0$ and $\mu^{DH}\_k \to a^{DH} + Z^{DH}$.
- In the legacy 2005 formulation, $a^{DH}$ and $Z^{DH}$ yield $\mu^{DH}\_k = -\dfrac{q_k^2 \kappa}{24\pi \varepsilon_0 \varepsilon_r k_B T}\left[2\chi_k + \dfrac{\sum_j x_j q_j^2 \sigma_j}{\sum_j x_j q_j^2}\right]$, which matches the Bulow form under the same assumptions (constant $\varepsilon_r$, no $\kappa$ composition dependence). This serves as a quick check when toggling between implementations.

**Derivation of $\partial a^{DH}/\partial x_k$ (Bulow 2019, scalar $\varepsilon_r$)**
- Starting point (scalar, composition-independent $\varepsilon_r$):
  - $a^{DH} = -\dfrac{1}{12\pi k_B T\,\varepsilon_0\varepsilon_r}\;\kappa \sum_j x_j z_j^{2}\,\chi_j$
  - $\kappa = \sqrt{\dfrac{\rho\,e^{2}}{k_B T\,\varepsilon_0\varepsilon_r}\;\sum_m x_m z_m^{2}}$
  - $\chi_j = \dfrac{3}{(\kappa a_j)^3}\left[\tfrac32 + \ln(1+\kappa a_j) - 2(1+\kappa a_j) + \tfrac12(1+\kappa a_j)^2\right]$.
- Composition derivative:
  \[
  \frac{\partial a^{DH}}{\partial x_k}
  = -\frac{1}{12\pi k_B T\,\varepsilon_0\varepsilon_r}
  \Bigg[
  \kappa\, z_k^2 \chi_k
  + \frac{\partial \kappa}{\partial x_k}\sum_j x_j z_j^2\chi_j
  + \kappa \sum_j x_j z_j^2 \frac{\partial \chi_j}{\partial \kappa}\frac{\partial \kappa}{\partial x_k}
  \Bigg].
  \]
- With $\partial \varepsilon_r/\partial x_k = 0$ and fixed $\rho$:
  - $\displaystyle \frac{\partial \kappa}{\partial x_k}
    = \frac{1}{2}\kappa \frac{z_k^2}{\sum_m x_m z_m^2}$.
  - $\displaystyle \frac{\partial \chi_j}{\partial \kappa}
    = 3\left[-\frac{\chi_j}{\kappa} + \frac{a_j}{(\kappa a_j)^3}\left(\frac{1}{1+\kappa a_j} - 1 + \kappa a_j\right)\right]$.
- Plugging these into the expression above yields the nonzero $\partial a^{DH}/\partial x_k$ used in the Bulow 2019 μ\_DH block. Only if $\kappa$ were (unphysically) treated composition‑independent **and** the explicit $x_k$ term were dropped would the derivative vanish.

**Planned Code Touchpoints (updated)**  
- Same as above for $\varepsilon_r$ mixing, plus:  
  - Debye chemical potential path: replace existing 2005-style direct $\mu_{\text{ion}}$ update with Bulow 2019 formulation using $\partial a^{DH}/\partial x_i$. Ensure $Z$ and $a_{\text{res}}$ derivatives remain consistent.  
  - Ensure no Bjerrum treatment is compiled/used in this pass.
