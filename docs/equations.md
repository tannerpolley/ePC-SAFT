# Equation Index

This file is generated from `docs/latex/equations.tex` by `scripts/sync_equation_registry.py`.
The LaTeX document remains the current source of truth; this Markdown view and `docs/equations_registry.yaml` stay aligned with it.

## Parameter Setup

### Base Mixture Quantities

#### `m_bar`
- Label: `eq:m_bar`
- Source: \cite{Gross2001}, Eq.~(A.5)
- Status: Close literature match
- Description: Provides a supporting relation used in hard-chain reference contribution.
- Change note: High textual similarity to a tagged equation in the cited local paper export.
- LaTeX: `docs/latex/equations.tex:27`
- C++: No `EqID` owner comment has been attached yet.

```tex
\bar{m}=\sum_{i} x_{i} m_{i}
```

#### `mw_bar`
- Label: `eq:mw_bar`
- Source: \cite{Figiel2025}, Eq.~(12) (derived helper)
- Status: Derived helper relation
- Description: Provides a supporting relation used in relative-permittivity and electrolyte reference calculations.
- Change note: This mean-molecular-weight helper is used to evaluate weighted relative-permittivity rules but is not a separately numbered equation in the paper.
- LaTeX: `docs/latex/equations.tex:38`
- C++: No `EqID` owner comment has been attached yet.

```tex
\overline{MW}=\sum_{j=1}^{N_{c}} x_{j}\,MW_{j}
```

#### `mw_solvent_bar`
- Label: `eq:mw_solvent_bar`
- Source: \cite{Figiel2025}, Eq.~(16)
- Status: Adapted implementation form
- Description: Provides a supporting relation used in relative-permittivity and electrolyte reference calculations.
- Change note: Lower similarity; likely algebraically adapted for implementation or combined terms.
- LaTeX: `docs/latex/equations.tex:49`
- C++: No `EqID` owner comment has been attached yet.

```tex
\overline{MW}_{sol}=\sum_{j\in\mathcal S}x_{j}\,MW_{j}
```

#### `solvent_ion_sets`
- Label: `eq:solvent_ion_sets`
- Source: \cite{Figiel2025}, Eq.~(11) (set-definition helper)
- Status: Derived helper relation
- Description: Provides a supporting relation used in relative-permittivity and electrolyte reference calculations.
- Change note: This solvent/ion index-set definition is implementation notation introduced in this documentation.
- LaTeX: `docs/latex/equations.tex:60`
- C++: No `EqID` owner comment has been attached yet.

```tex
\mathcal S=\{j:\,z_{j}=0\},\qquad \mathcal I=\{j:\,z_{j}\neq 0\}
```

#### `x_solvent_total`
- Label: `eq:x_solvent_total`
- Source: \cite{Figiel2025}, Eq.~(16)
- Status: Adapted implementation form
- Description: Provides a supporting relation used in relative-permittivity and electrolyte reference calculations.
- Change note: Lower similarity; likely algebraically adapted for implementation or combined terms.
- LaTeX: `docs/latex/equations.tex:71`
- C++: No `EqID` owner comment has been attached yet.

```tex
x_{sol}=\sum_{j\in\mathcal S}x_{j}
```

### Size and Diameter Rules

#### Segment Diameter

##### `d_segment`
- Label: `eq:d_segment`
- Source: \cite{Gross2001}, Eq.~(A.9)
- Status: Close literature match
- Description: Provides a supporting relation used in hard-chain reference contribution.
- Change note: High textual similarity to a tagged equation in the cited local paper export.
- LaTeX: `docs/latex/equations.tex:89`
- C++: No `EqID` owner comment has been attached yet.

```tex
d_{i}=\sigma_{i}\left[1-0.12 \exp \left(-3 \frac{\epsilon_{i}}{k T}\right)\right]
```

##### `d_segment_dT`
- Label: `eq:d_segment_dT`
- Source: \cite{Gross2001}, Eq.~(A.9) (derived temperature differential)
- Status: Manual literature match
- Description: Provides a differential relation needed for temperature differential calculations.
- Change note: Computed by differentiating Eq. (A.9) with respect to temperature.
- LaTeX: `docs/latex/equations.tex:100`
- C++: No `EqID` owner comment has been attached yet.

```tex
\left(\frac{\partial d_i}{\partial T}\right)_{\rho,x}
    \equiv d_{i,T}
    =\sigma_{i}\left(3\frac{\epsilon_{i}}{kT^2}\right)\left[-0.12\exp\left(-3\frac{\epsilon_{i}}{kT}\right)\right]
```

##### `half_d_identity`
- Label: `eq:half_d_identity`
- Source: \cite{Figiel2025}, Eq.~(20)
- Status: Adapted implementation form
- Description: Provides a supporting relation used in temperature differential equations.
- Change note: Moved into the upstream diameter setup block so this simple diameter identity is defined with the other diameter relations before it is reused downstream.
- LaTeX: `docs/latex/equations.tex:113`
- C++: No `EqID` owner comment has been attached yet.

```tex
\frac{1}{2}d_{i}=\left(\frac{d_{i}d_{i}}{d_{i}+d_{i}}\right)
```

#### Ion Diameter Rule

##### `d_ion_rule`
- Label: `eq:d_ion_rule`
- Source: \cite{Bulow2019}, Eq.~(3)
- Status: Adapted implementation form
- Description: Provides the grouped ion-diameter rule used in Debye-Huckel electrolyte calculations.
- Change note: Consolidates the documented sigma, constant-factor, and Barker-Henderson ion-diameter options into one visible case-set presentation.
- LaTeX: `docs/latex/equations.tex:127`
- C++: No `EqID` owner comment has been attached yet.

```tex
d_{\mathrm{ion},i}\equiv d_{i}=
    \begin{cases}
        \sigma_{i}, & \text{sigma rule}, \\[6pt]
        \left(1-0.12\right)\sigma_{i}=0.88\sigma_{i}, & \text{constant-factor rule}, \\[6pt]
        \left[1-0.12 \exp \left(-3 \frac{\epsilon_{i}}{k T}\right)\right]\sigma_{i}, & \text{Barker-Henderson rule},
    \end{cases}
    \qquad i\in\mathcal{I}
```

### Pair Mixing Rules

#### `d_ij`
- Label: `eq:d_ij`
- Source: \cite{Gross2001}, Eq.~(A.14)
- Status: Adapted notation
- Description: Provides a supporting relation used in association contribution.
- Change note: Mapped to the arithmetic combining rule form; this file writes it for effective diameters in association context.
- LaTeX: `docs/latex/equations.tex:147`
- C++: No `EqID` owner comment has been attached yet.

```tex
d_{i j}=\left(d_{i i}+d_{j j}\right) / 2 .
```

#### `sigma_ij`
- Label: `eq:sigma_ij`
- Source: \cite{Gross2001}, Eq.~(A.15)
- Status: Adapted notation
- Description: Provides a supporting relation used in dispersion contribution.
- Change note: Moderate-to-high similarity; notation/arrangement appears adapted from the cited equation.
- LaTeX: `docs/latex/equations.tex:158`
- C++: No `EqID` owner comment has been attached yet.

```tex
\sigma_{i j}=\frac{1}{2}\left(\sigma_{i}+\sigma_{j}\right) \cdot \left(1- l_{ij} \right)
```

#### `epsilon_ij_mixing`
- Label: `eq:epsilon_ij_mixing`
- Source: \cite{Gross2001}, Eq.~(A.15)
- Status: Adapted notation
- Description: Provides a supporting relation used in dispersion contribution.
- Change note: This is the base Berthelot-style pair-dispersion combining rule before any ePC-SAFT ionic override is applied.
- LaTeX: `docs/latex/equations.tex:169`
- C++: No `EqID` owner comment has been attached yet.

```tex
\epsilon_{ij}^{\mathrm{base}}
    =
    \sqrt{\epsilon_{i}\,\epsilon_{j}}\left(1-k_{ij}\right)
```

#### `epsilon_ij_ionic_zero`
- Label: `eq:epsilon_ij_ionic_zero`
- Source: ePC-SAFT implementation override in \texttt{pair\_epsilon\_cpp(...)} (not a standalone literature equation)
- Status: Project-specific modification
- Description: Provides a supporting relation used in dispersion contribution.
- Change note: The active implementation suppresses short-range dispersion for same-sign ionic pairs by overriding the base combining rule with zero; this is why the second equation is simply \(=0\).
- LaTeX: `docs/latex/equations.tex:182`
- C++: No `EqID` owner comment has been attached yet.

```tex
\epsilon_{ij}=0
    \qquad \text{for } z_{i}z_{j}>0,
```

#### `epsilon_assoc_mixing`
- Label: `eq:epsilon_assoc_mixing`
- Source: \cite{Gross2002}, Eq.~(2)
- Status: Manual literature match
- Description: Provides a supporting relation used in association contribution.
- Change note: Mapped manually to the Gross 2002 cross-association energy combining rule.
- LaTeX: `docs/latex/equations.tex:196`
- C++: No `EqID` owner comment has been attached yet.

```tex
\varepsilon^{A_{i}B_{j}}=\frac{1}{2}(\varepsilon^{A_{i}B_{i}}+\varepsilon^{A_{j}B_{j}})(1-k_{ij}^{\mathrm{hb}})
```

#### `kappa_assoc_mixing`
- Label: `eq:kappa_assoc_mixing`
- Source: \cite{Gross2002}, Eq.~(3)
- Status: Manual literature match
- Description: Provides a supporting relation used in association contribution.
- Change note: Mapped manually to the Gross 2002 cross-association volume combining rule.
- LaTeX: `docs/latex/equations.tex:207`
- C++: No `EqID` owner comment has been attached yet.

```tex
k^{A_{i}B_{j}}=\sqrt{k^{A_{i}B_{i}}k^{A_{j}B_{j}}}\quad\left(\frac{\sqrt{\sigma_{i}\sigma_{j}}}{1/2(\sigma_{i}+\sigma_{j})}\right)^3
```

### Relative Permittivity

#### `epsr_mix_rule`
- Label: `eq:epsr_mix_rule`
- Source: \cite{Ascani2021}, Eq.~(3), Eq.~(11), Eq.~(13); \cite{Figiel2025}, Eq.~(11)
- Status: Adapted implementation form
- Description: Provides the grouped relative-permittivity mixing rule used by the electrolyte contributions.
- Change note: Consolidates the documented upstream relative-permittivity options into one visible grouped display while keeping helper quantities upstream.
- LaTeX: `docs/latex/equations.tex:222`
- C++: No `EqID` owner comment has been attached yet.

```tex
\varepsilon_{r}=
    \begin{cases}
        \sum_{j=1}^{N_{c}}\varepsilon_{r,j}x_{j}, & \text{mole-fraction rule}, \\[10pt]
        \dfrac{\sum_{j=1}^{N_{c}}\varepsilon_{r,j}\,x_{j}\,MW_{j}}{\overline{MW}}, & \text{mass-fraction rule}, \\[12pt]
        x_{sol}\,\varepsilon_{r,sol}^{(w)} + \sum_{j\in\mathcal I}\varepsilon_{r,j}\,x_{j}, & \text{solvent-ion rule}, \\[10pt]
        \dfrac{\varepsilon_{r,\mathrm{solvent,mix}}^{\mathrm{salt-free}}}{1+7.01\,x_{\mathrm{ion}}}, & \text{ion-suppressed rule}.
    \end{cases}
```

#### Mole-Fraction Rule

##### `depsr_dxi_mole`
- Label: `eq:depsr_dxi_mole`
- Source: \cite{Ascani2021}, Eq.~(11) (derived differential)
- Status: Manual literature match
- Description: Specifies dielectric-property mixing or derivative form for debye and huckel electrolyte term contribution.
- Change note: Direct derivative of the mole-fraction dielectric mixing rule.
- LaTeX: `docs/latex/equations.tex:242`
- C++: No `EqID` owner comment has been attached yet.

```tex
\left(\frac{\partial \varepsilon_{r}}{\partial x_{i}}\right)_{x_{j\neq i}}
    = \varepsilon_{r, i}
```

#### Mass-Fraction Rule

##### `depsr_dxi_mass`
- Label: `eq:depsr_dxi_mass`
- Source: \cite{Ascani2021}, Eq.~(12) (derived differential)
- Status: Manual literature match
- Description: Specifies dielectric-property mixing or derivative form for debye and huckel electrolyte term contribution.
- Change note: Direct derivative of the mass-fraction dielectric mixing rule.
- LaTeX: `docs/latex/equations.tex:257`
- C++: No `EqID` owner comment has been attached yet.

```tex
\left(\frac{\partial \varepsilon_{r}}{\partial x_{i}}\right)_{x_{j\neq i}}
    =
    \frac{MW_{i}}{\overline{MW}}\left(\varepsilon_{r,i}-\varepsilon_{r}\right)
```

#### Solvent-Ion Rule

##### `depsr_dxi_combo`
- Label: `eq:depsr_dxi_combo`
- Source: \cite{Ascani2021}, Eq.~(13) (derived differential)
- Status: Manual literature match
- Description: Specifies dielectric-property mixing or derivative form for debye and huckel electrolyte term contribution.
- Change note: Direct derivative of the mixed solvent-plus-ion dielectric mixing rule.
- LaTeX: `docs/latex/equations.tex:273`
- C++: No `EqID` owner comment has been attached yet.

```tex
\left(\frac{\partial \varepsilon_{r}}{\partial x_{i}}\right)_{x_{j\neq i}}
    =
    \varepsilon_{r,sol}^{(w)}
    +
    x_{sol}\,\frac{MW_{i}}{\overline{MW}_{sol}}
    \left(\varepsilon_{r,i}-\varepsilon_{r,sol}^{(w)}\right),
    \qquad i\in\mathcal S
```

##### `epsr_solvent_mass`
- Label: `eq:epsr_solvent_mass`
- Source: \cite{Ascani2021}, Eq.~(13)
- Status: Adapted notation
- Description: Provides a supporting relation used in relative-permittivity calculations for electrolyte contributions.
- Change note: Algebraic expansion of Eq. (13) introducing explicit solvent-only weighted relative permittivity.
- LaTeX: `docs/latex/equations.tex:290`
- C++: No `EqID` owner comment has been attached yet.

```tex
\varepsilon_{r,sol}^{(w)} \equiv \sum_{j\in\mathcal S}\varepsilon_{r,j}\,w_{j}^{sol} = \frac{\sum_{j\in\mathcal S}\varepsilon_{r,j}\,x_{j}\,MW_{j}}{\sum_{j\in\mathcal S}x_{j}\,MW_{j}} = \frac{\sum_{j\in\mathcal S}\varepsilon_{r,j}\,x_{j}\,MW_{j}}{\overline{MW}_{sol}}
```

#### Ion-Suppressed Empirical Rule

##### `epsr_salt_free`
- Label: `eq:epsr_salt_free`
- Source: \cite{Figiel2025}, Eq.~(12)
- Status: New literature extension
- Description: Provides a supporting relation used in relative-permittivity calculations for electrolyte contributions.
- Change note: Salt-free solvent-mixture relative permittivity used by the ion-suppressed upstream rule.
- LaTeX: `docs/latex/equations.tex:304`
- C++: No `EqID` owner comment has been attached yet.

```tex
\varepsilon_{r,\mathrm{solvent,mix}}^{\mathrm{salt-free}}=\sum_{\mathrm{solvent}}w_{\mathrm{solvent}}^{\mathrm{salt-free}}\cdot\varepsilon_{r,\mathrm{solvent}}
```

##### `epsr_mix_suppressed`
- Label: `eq:epsr_mix_suppressed`
- Source: \cite{Bulow2021} (equation number unresolved in local corpus)
- Status: New literature extension
- Description: Provides a supporting relation used in debye and huckel electrolyte term contribution.
- Change note: Changed from earlier dielectric-mixing rules to the 2025 ion-suppression mixing form.
- LaTeX: `docs/latex/equations.tex:315`
- C++: No `EqID` owner comment has been attached yet.

```tex
\varepsilon_{r,\mathrm{mix}}(\mathbf{x})
    =
    \frac{\varepsilon_{sf}(\mathbf{x})}{1+7.01\,x_{\mathrm{ion}}(\mathbf{x})},
    \qquad
    \varepsilon_{sf}\equiv \varepsilon_{r,\mathrm{solvent,mix}}^{\mathrm{salt\text{-}free}} .
```

##### `epsr_sf`
- Label: `eq:epsr_sf`
- Source: \cite{Figiel2025}, Eq.~(12)
- Status: Manual literature match
- Description: Provides a supporting relation used in debye and huckel electrolyte term contribution.
- Change note: Mapped manually to the salt-free solvent dielectric mixture definition.
- LaTeX: `docs/latex/equations.tex:330`
- C++: No `EqID` owner comment has been attached yet.

```tex
\varepsilon_{sf}(\mathbf{x})
    =
    \sum_{s\in\mathcal{S}} w_{s}^{sf}\,\varepsilon_{r,s},
    \qquad
    w_{s}^{sf}
    =
    \frac{x_{s} M_{s}}{\displaystyle \sum_{m\in\mathcal{S}} x_{m} M_{m}},
    \qquad
    D\equiv \sum_{m\in\mathcal{S}} x_{m} M_{m} .
```

##### `x_ion_total`
- Label: `eq:x_ion_total`
- Source: \cite{Figiel2025}, Eq.~(11) (derived helper)
- Status: Manual literature match
- Description: Provides a differential relation needed for debye and huckel electrolyte term contribution calculations.
- Change note: Ion-fraction helper used in the Eq. (11) differential form.
- LaTeX: `docs/latex/equations.tex:349`
- C++: No `EqID` owner comment has been attached yet.

```tex
x_{\mathrm{ion}}(\mathbf{x})=\sum_{m\in\mathcal{I}} x_{m},
    \qquad
    \left(\frac{\partial x_{\mathrm{ion}}}{\partial x_{i}}\right)_{T,v,x_{j\neq i}}
    =
    \begin{cases}
        1, & i\in\mathcal{I},    \\
        0, & i\notin\mathcal{I}.
    \end{cases}
```

##### `depsr_sf_dxi`
- Label: `eq:depsr_sf_dxi`
- Source: \cite{Figiel2025}, Eq.~(12) (derived differential)
- Status: Manual literature match
- Description: Provides a differential relation needed for debye and huckel electrolyte term contribution calculations.
- Change note: Derivative of solvent-mixture dielectric expression in Eq. (12).
- LaTeX: `docs/latex/equations.tex:366`
- C++: No `EqID` owner comment has been attached yet.

```tex
\left(\frac{\partial \varepsilon_{sf}}{\partial x_{i}}\right)_{T,v,x_{j\neq i}}
    =
    \sum_{s\in\mathcal{S}} \varepsilon_{r,s}
    \left(\frac{\partial w_{s}^{sf}}{\partial x_{i}}\right)_{T,v,x_{j\neq i}}
    =
    \begin{cases}
        \dfrac{M_{i}}{D}\left(\varepsilon_{r,i}-\varepsilon_{sf}\right), & i\in\mathcal{S},    \\[10pt]
        0,                                                             & i\notin\mathcal{S}.
    \end{cases}
```

##### `depsr_mix_dxi`
- Label: `eq:depsr_mix_dxi`
- Source: \cite{Figiel2025}, Eq.~(11)-Eq.~(12) (derived differential)
- Status: Manual literature match
- Description: Provides a differential relation needed for debye and huckel electrolyte term contribution calculations.
- Change note: Derivative of the 2025 dielectric mixing rule via chain rule.
- LaTeX: `docs/latex/equations.tex:385`
- C++: No `EqID` owner comment has been attached yet.

```tex
\left(\frac{\partial \varepsilon_{r,\mathrm{mix}}}{\partial x_{i}}\right)_{T,v,x_{j\neq i}}
    =
    \frac{1}{1+7.01\,x_{\mathrm{ion}}}
    \left(\frac{\partial \varepsilon_{sf}}{\partial x_{i}}\right)_{T,v,x_{j\neq i}}
    -\frac{7.01\,\varepsilon_{sf}}{\left(1+7.01\,x_{\mathrm{ion}}\right)^2}
    \left(\frac{\partial x_{\mathrm{ion}}}{\partial x_{i}}\right)_{T,v,x_{j\neq i}} .
```

##### `depsr_mix_dxi_piecewise`
- Label: `eq:depsr_mix_dxi_piecewise`
- Source: \cite{Figiel2025}, Eq.~(11)-Eq.~(12) (derived closed form)
- Status: Manual literature match
- Description: Provides a differential relation needed for debye and huckel electrolyte term contribution calculations.
- Change note: Piecewise closed-form derivative obtained by combining Eq. (11) and Eq. (12).
- LaTeX: `docs/latex/equations.tex:401`
- C++: No `EqID` owner comment has been attached yet.

```tex
\left(\frac{\partial \varepsilon_{r,\mathrm{mix}}}{\partial x_{i}}\right)_{T,v,x_{j\neq i}}
    =
    \begin{cases}
        \dfrac{1}{1+7.01\,x_{\mathrm{ion}}}\;
        \dfrac{M_{i}}{\displaystyle \sum_{m\in\mathcal{S}} x_{m} M_{m}}
        \left(\varepsilon_{r,i}-\varepsilon_{sf}\right),
         & i\in\mathcal{S}, \\[14pt]
        -\dfrac{7.01\,\varepsilon_{sf}}{\left(1+7.01\,x_{\mathrm{ion}}\right)^{2}},
         & i\in\mathcal{I}.
    \end{cases}
```

## Density

### Packing-Fraction Moments

#### `zeta_n`
- Label: `eq:zeta_n`
- Source: \cite{Gross2001}, Eq.~(A.8)
- Status: Close literature match
- Description: Provides a supporting relation used in hard-chain reference contribution.
- Change note: High textual similarity to a tagged equation in the cited local paper export.
- LaTeX: `docs/latex/equations.tex:429`
- C++: No `EqID` owner comment has been attached yet.

```tex
\zeta_{n}=\frac{\pi}{6} \rho \sum_{i} x_{i} m_{i} d_{i}^n \quad n \in\{0,1,2,3\}
```

#### `zeta_n_xk`
- Label: `eq:zeta_n_xk`
- Source: \cite{Gross2001}, Eq.~(A.34)
- Status: Close literature match
- Description: Provides a differential relation needed for composition differential calculations.
- Change note: Ownership moved upstream so the composition derivatives are defined before the downstream property sections that use them.
- LaTeX: `docs/latex/equations.tex:440`
- C++: No `EqID` owner comment has been attached yet.

```tex
\zeta_{n,xk}=\left(\frac{\partial\zeta_{n}}{\partial x_{k}}\right)_{T,\boldsymbol{\rho},x_{j\neq k}}=\frac{\pi}{6}\rho m_{k}(d_{k})^{n}\quad n\in\{0,1,2,3\}
```

#### `zeta_n_dT`
- Label: `eq:zeta_n_dT`
- Source: \cite{Gross2001}, Eq.~(A.53)
- Status: Manual literature match
- Description: Provides a differential relation needed for temperature differential calculations.
- Change note: Mapped manually to the temperature derivative of zeta moments.
- LaTeX: `docs/latex/equations.tex:451`
- C++: No `EqID` owner comment has been attached yet.

```tex
\zeta_{n,T}=\frac{\partial\zeta_{n}}{\partial T}=\frac{\pi}{6}\rho\sum_{i}x_{i}m_{i}nd_{i,T}\left(d_{i}\right)^{n-1}\quad n\in\{1,2,3\}
```

### Packing Fraction

#### `zeta3_eta`
- Label: `eq:zeta3_eta`
- Source: \cite{Gross2001}, Eq.~(A.20)-Eq.~(A.22) (identity in notation)
- Status: Derived helper relation
- Description: Provides the canonical packing-fraction identity used before evaluating contribution expressions.
- Change note: Makes the \zeta_3 and \eta notation equivalence explicit in the new initial-density section.
- LaTeX: `docs/latex/equations.tex:465`
- C++: No `EqID` owner comment has been attached yet.

```tex
\eta=\zeta_{3}
```

### Density Conversion

#### `rho_from_eta`
- Label: `eq:rho_from_eta`
- Source: \cite{Gross2001}, Eq.~(A.20)
- Status: Close literature match
- Description: Provides a supporting relation used in dispersion contribution.
- Change note: High textual similarity to a tagged equation in the cited local paper export.
- LaTeX: `docs/latex/equations.tex:479`
- C++: No `EqID` owner comment has been attached yet.

```tex
\rho=\frac{6}{\pi} \eta\left(\sum_{i} x_{i} m_{i} d_{i}^3\right)^{-1}
```

#### `rho_reduced`
- Label: `eq:rho_reduced`
- Source: \cite{Gross2001}, Eq.~(A.21)
- Status: Close literature match
- Description: Provides a supporting relation used in dispersion contribution.
- Change note: High textual similarity to a tagged equation in the cited local paper export.
- LaTeX: `docs/latex/equations.tex:490`
- C++: No `EqID` owner comment has been attached yet.

```tex
\rho_{N}=N_{\mathrm{AV}}\rho,
    \qquad
    \hat{\rho}=\rho_{N}\cdot 10^{-30}
```

### Pressure-Density Closure

#### `pressure_from_z`
- Label: `eq:pressure_from_z`
- Source: \cite{Gross2001}, Eq.~(A.21)
- Status: Adapted implementation form
- Description: Provides the pressure closure used in the initial density solve.
- Change note: Moved upstream from the pressure section so the density-solver closure lives with the rest of the initial density-setup equations.
- LaTeX: `docs/latex/equations.tex:506`
- C++: `src/epcsaft/native/epcsaft_Z.cpp:76` (double p_cpp(double t, double rho, vector<double> x, const add_args &cppargs) {)

```tex
P = ZRT\rho = Zk_{B}T\rho_{N}
```

#### `z_from_rho`
- Label: `eq:z_from_rho`
- Source: \cite{Gross2001}, Eq.~(A.22)
- Status: Close literature match
- Description: Provides the combined compressibility-factor closure used in the initial density solve.
- Change note: The packing-fraction and density forms are now shown together in one equation so the closure reads as one relation instead of three separate displayed identities.
- LaTeX: `docs/latex/equations.tex:517`
- C++: `src/epcsaft/native/epcsaft_Z.cpp:63` (CompressibilityFactorResult compressibility_factor_result_cpp(double t, double rho, vector<double> x, const add_args &cppargs) {)

```tex
Z
    =
    1+\eta\left(\frac{\partial\tilde{a}^\mathrm{res}}{\partial\eta}\right)_{T,x_{i}}
    =
    1+\rho\left(\frac{\partial\tilde{a}^\mathrm{res}}{\partial\rho}\right)_{T,x_{i}}
```

#### `density_solve_residual`
- Label: `eq:density_solve_residual`
- Source: Derived from Eq.~\eqref{eq:pressure_from_z}
- Status: Derived solver relation
- Description: Provides the explicit scalar residual used to locate density roots at specified \(T\), \(P\), and \(\mathbf{x}\).
- Change note: Added as a solver-facing closeout equation so the initial-density section ends with the actual root-finding constraint rather than only helper relations.
- LaTeX: `docs/latex/equations.tex:534`
- C++: No `EqID` owner comment has been attached yet.

```tex
R_{\rho}\!\left(T,\rho,\mathbf{x};P^{\mathrm{spec}}\right)
    =
    P^{\mathrm{spec}}
    -
    Z(T,\rho,\mathbf{x})RT\rho
    =
    0
```

## Contribution Intermediates

### Hard-Chain

#### Hard-Sphere Helmholtz Term \texorpdfstring{\(\tilde{a

##### `ares_hs`
- Label: `eq:ares_hs`
- Source: \cite{Gross2001}, Eq.~(A.6)
- Status: Close literature match
- Description: Provides a residual Helmholtz-energy relation for hard-chain reference contribution.
- Change note: High textual similarity to a tagged equation in the cited local paper export.
- LaTeX: `docs/latex/equations.tex:562`
- C++: `src/epcsaft/native/epcsaft_ares.cpp:17` (static double ares_hs_cpp(const HardChainState &hc_state) {)

```tex
\begin{aligned}
        \tilde{a}^{\mathrm{hs}}
        &=
        \frac{1}{\zeta_{0}}
        \left[
            \frac{3 \zeta_{1} \zeta_{2}}{\left(1-\zeta_{3}\right)}
            +\frac{\zeta_{2}^3}{\zeta_{3}\left(1-\zeta_{3}\right)^2}
            +\left(\frac{\zeta_{2}^3}{\zeta_{3}^2}-\zeta_{0}\right) \ln \left(1-\zeta_{3}\right)
        \right]
    \end{aligned}
```

##### `hs_ares_dadrho`
- Label: `eq:hs_ares_dadrho`
- Source: \cite{Gross2001}, Eq.~(A.26)
- Status: Manual literature match
- Description: Provides the density differential of the hard-sphere Helmholtz term.
- Change note: Renamed and re-homed so the hard-sphere density differential sits directly below the owning hard-sphere expression.
- LaTeX: `docs/latex/equations.tex:582`
- C++: `src/epcsaft/native/contributions/epcsaft_contrib_hc.cpp:82` (static double dadrho_hs_cpp(const HardChainState &hc_state) {)

```tex
\begin{aligned}
        \rho\left(\frac{\partial\tilde{a}^{hs}}{\partial\rho}\right)_{T,x}
        &=
        \frac{\zeta_{3}}{(1-\zeta_{3})}
        +\frac{3\zeta_{1}\zeta_{2}}{\zeta_{0}(1-\zeta_{3})^{2}}
        +\frac{3\zeta_{2}^{3}-\zeta_{3}\zeta_{2}^{3}}{\zeta_{0}(1-\zeta_{3})^{3}}
    \end{aligned}
```

##### `hs_ares_dxk`
- Label: `eq:hs_ares_dxk`
- Source: \cite{Gross2001}, Eq.~(A.36)
- Status: Manual literature match
- Description: Provides the composition differential of the hard-sphere Helmholtz term.
- Change note: Renamed and re-homed so the hard-sphere composition differential sits directly below the owning hard-sphere expression.
- LaTeX: `docs/latex/equations.tex:599`
- C++: `src/epcsaft/native/contributions/epcsaft_contrib_hc.cpp:103` (static vector<double> dadx_hs_cpp(const MixtureState &thermo, const HardChainState &hc_state, const add_args &cppargs) {)

```tex
\begin{aligned}
        \left(\frac{\partial\tilde{a}^{\mathrm{hs}}}{\partial x_{k}}\right)_{T,\boldsymbol{\rho},x_{j\neq k}}
        &=- \frac{\zeta_{0,xk}}{\zeta_{0}}\tilde{a}^{\mathrm{hs}}+\frac{1}{\zeta_{0}}\biggl[
        \frac{3(\zeta_{1,xk}\zeta_{2}+\zeta_{1}\zeta_{2,xk})}{(1-\zeta_{3})}
        +\frac{3\zeta_{1}\zeta_{2}\zeta_{3,xk}}{\left(1-\zeta_{3}\right)^{2}}
        \\
        &\quad +\frac{3\zeta_{2}^{2}\zeta_{2,xk}}{\zeta_{3}(1-\zeta_{3})^{2}}
        +\frac{\zeta_{2}^{3}\zeta_{3,xk}(3\zeta_{3}-1)}{\zeta_{3}^{2}(1-\zeta_{3})^{3}}
        \\
        &\quad +\left(
        \frac{3{\zeta_{2}}^{2}{\zeta_{2,xk}}{\zeta_{3}}-2{\zeta_{2}}^{3}{\zeta_{3,xk}}}{{\zeta_{3}}^{3}}
        -{\zeta_{0,xk}}
        \right)\ln(1-{\zeta_{3}})
        \\
        &\quad +\left(\zeta_{0}-\frac{\zeta_{2}^{3}}{\zeta_{3}^{2}}\right)\frac{\zeta_{3,xk}}{(1-\zeta_{3})}
        \biggr]
    \end{aligned}
```

##### `hs_ares_dT`
- Label: `eq:hs_ares_dT`
- Source: \cite{Gross2001}, Eq.~(A.55)
- Status: Manual literature match
- Description: Provides the temperature differential of the hard-sphere Helmholtz term.
- Change note: Renamed and re-homed so the hard-sphere temperature differential sits directly below the owning hard-sphere expression.
- LaTeX: `docs/latex/equations.tex:626`
- C++: `src/epcsaft/native/contributions/epcsaft_contrib_hc.cpp:90` (static double dadt_hs_cpp(const HardChainState &hc_state, const vector<double> &dzeta_dt) {)

```tex
\begin{aligned}
        \left(\frac{\partial\tilde{a}^{hs}}{\partial T}\right)_{\rho,x_{i}}
        &=\frac{1}{\zeta_{0}}\Bigg[
        \frac{3(\zeta_{1,T}\zeta_{2}+\zeta_{1}\zeta_{2,T})}{(1-\zeta_{3})}
        +\frac{3\zeta_{1}\zeta_{2}\zeta_{3,T}}{(1-\zeta_{3})^{2}}
        \\
        &\quad
        +\frac{3{\zeta_{2}}^{2}{\zeta_{2,T}}}{\zeta_{3}{(1-\zeta_{3})}^{2}}
        +\frac{{\zeta_{2}}^{3}{\zeta_{3,T}}(3{\zeta_{3}}-1)}{{\zeta_{3}}^{2}{(1-\zeta_{3})}^{3}}
        \\
        &\quad
        +\left(\frac{3{\zeta_{2}}^{2}{\zeta_{2,T}}{\zeta_{3}}-2{\zeta_{2}}^{3}{\zeta_{3,T}}}{{\zeta_{3}}^{3}}\right)\ln(1-{\zeta_{3}})
        +\left({\zeta_{0}}-\frac{{\zeta_{2}}^{3}}{{\zeta_{3}}^{2}}\right)\frac{\zeta_{3,T}}{(1-\zeta_{3})}
        \Bigg]
    \end{aligned}
```

#### Hard-Sphere Contact Value \texorpdfstring{\(\mathrm{g

##### `g_hs_contact`
- Label: `eq:g_hs_contact`
- Source: \cite{Gross2001}, Eq.~(A.7)
- Status: Close literature match
- Description: Provides a supporting relation used in hard-chain reference contribution.
- Change note: High textual similarity to a tagged equation in the cited local paper export.
- LaTeX: `docs/latex/equations.tex:654`
- C++: No `EqID` owner comment has been attached yet.

```tex
\mathrm{g}_{i j}^{\mathrm{hs}}=\frac{1}{\left(1-\zeta_{3}\right)}+\left(\frac{d_{i} d_{j}}{d_{i}+d_{j}}\right) \frac{3 \zeta_{2}}{\left(1-\zeta_{3}\right)^2} + \left(\frac{d_{i} d_{j}}{d_{i}+d_{j}}\right)^2 \frac{2 \xi_{2}{ }^2}{\left(1-\xi_{3}\right)^3}
```

##### `ghs_contact_dadrho`
- Label: `eq:ghs_contact_dadrho`
- Source: \cite{Gross2001}, Eq.~(A.27)
- Status: Manual literature match
- Description: Provides the density differential of the hard-sphere contact value.
- Change note: Renamed and re-homed so the contact-value density differential sits directly below the owning contact-value relation.
- LaTeX: `docs/latex/equations.tex:665`
- C++: `src/epcsaft/native/contributions/epcsaft_contrib_hc.cpp:164` (double hs_contact_density_derivative_cpp(double pair_diameter, double zeta2, double zeta3) {)

```tex
\begin{aligned}
        \rho\left(\frac{\partial g_{ij}^{\mathrm{hs}}}{\partial\rho}\right)_{T,x}
        &=
        \frac{\zeta_{3}}{\left(1-\zeta_{3}\right)^{2}}
        +
        \left(\frac{d_{i}d_{j}}{d_{i}+d_{j}}\right)
        \left(
            \frac{3\zeta_{2}}{\left(1-\zeta_{3}\right)^{2}}
            +
            \frac{6\zeta_{2}\zeta_{3}}{\left(1-\zeta_{3}\right)^{3}}
        \right)
        \\
        &\quad +
        \left(\frac{d_{i}d_{j}}{d_{i}+d_{j}}\right)^{2}
        \left(
            \frac{4{\zeta_{2}}^{2}}{\left(1-{\zeta_{3}}\right)^{3}}
            +
            \frac{6{\zeta_{2}}^{2}\zeta_{3}}{\left(1-{\zeta_{3}}\right)^{4}}
        \right)
    \end{aligned}
```

##### `ghs_contact_dxk`
- Label: `eq:ghs_contact_dxk`
- Source: \cite{Gross2001}, Eq.~(A.37)
- Status: Close literature match
- Description: Provides the composition differential of the hard-sphere contact value.
- Change note: Renamed and re-homed so the contact-value composition differential sits directly below the owning contact-value relation.
- LaTeX: `docs/latex/equations.tex:695`
- C++: `src/epcsaft/native/contributions/epcsaft_contrib_hc.cpp:186` (double hs_contact_composition_derivative_cpp()

```tex
\begin{aligned}
        \left(\frac{\partial g_{ij}^{\mathrm{hs}}}{\partial x_{k}}\right)_{T,\boldsymbol{\rho},\boldsymbol{\rho}_{j\neq k}}
        &=\frac{\zeta_{3,x,k}}{\left(1-\zeta_{3}\right)^{2}}
        +\left(\frac{d_{i}d_{j}}{d_{i}+d_{j}}\right)
        \left(\frac{3\zeta_{2,x,k}}{(1-\zeta_{3})^{2}}+\frac{6\zeta_{2}\zeta_{3,x,k}}{(1-\zeta_{3})^{3}}\right)
        \\
        &\quad +\left(\frac{d_{i}d_{j}}{d_{i}+d_{j}}\right)^{2}
        \left(\frac{4\zeta_{2}\zeta_{2,x,k}}{(1-\zeta_{3})^{3}}+\frac{6\zeta_{2}^{2}\zeta_{3,x,k}}{(1-\zeta_{3})^{4}}\right)
    \end{aligned}
```

##### `ghs_contact_dT`
- Label: `eq:ghs_contact_dT`
- Source: \cite{Gross2001}, Eq.~(A.57)
- Status: Manual literature match
- Description: Provides the temperature differential of the hard-sphere contact value.
- Change note: Renamed and re-homed so the contact-value temperature differential sits directly below the owning contact-value relation.
- LaTeX: `docs/latex/equations.tex:714`
- C++: `src/epcsaft/native/contributions/epcsaft_contrib_hc.cpp:171` (double hs_contact_time_derivative_cpp()

```tex
\begin{aligned}
        \left(\frac{\partial g_{ii}^{\mathrm{hs}}}{\partial T}\right)_{\rho,x}
        &=\frac{\zeta_{3,T}}{\left(1-\zeta_{3}\right)^{2}}
        +\left(\frac{1}{2}d_{i,T}\right)\frac{3\zeta_{2}}{\left(1-\zeta_{3}\right)^{2}}
        \\
        &\quad +\left(\frac{1}{2}d_{i}\right)\left(\frac{3\zeta_{2,T}}{\left(1-\zeta_{3}\right)^{2}}+\frac{6\zeta_{2}\zeta_{3,T}}{\left(1-\zeta_{3}\right)^{3}}\right)
        +\left(\frac{1}{2}d_{i}d_{i,T}\right)\frac{2\zeta_{2}^{2}}{\left(1-\zeta_{3}\right)^{3}}
        \\
        &\quad +\left(\frac{1}{2}d_{i}\right)^{2}\left(\frac{4\xi_{2}\xi_{2,T}}{\left(1-\xi_{3}\right)^{3}}+\frac{6\xi_{2}{}^{2}\xi_{3,T}}{\left(1-\xi_{3}\right)^{4}}\right)
    \end{aligned}
```

### Dispersion

#### First-Order Compressibility Prefactor

##### `c1_disp`
- Label: `eq:c1_disp`
- Source: \cite{Gross2001}, Eq.~(A.11)
- Status: Adapted notation
- Description: Provides a residual Helmholtz-energy relation for dispersion contribution.
- Change note: Moderate-to-high similarity; notation/arrangement appears adapted from the cited equation.
- LaTeX: `docs/latex/equations.tex:742`
- C++: No `EqID` owner comment has been attached yet.

```tex
C_{1} = \left(1+\bar{m} \frac{8 \eta-2 \eta^2}{(1-\eta)^4}+\right.\left.\quad(1-\bar{m}) \frac{20 \eta-27 \eta^2+12 \eta^3-2 \eta^4}{[(1-\eta)(2-\eta)]^2}\right)
```

##### `c2_disp`
- Label: `eq:c2_disp`
- Source: \cite{Gross2001}, Eq.~(A.31)
- Status: Close literature match
- Description: Provides a differential relation needed for dispersion contribution calculations.
- Change note: High textual similarity to a tagged equation in the cited local paper export.
- LaTeX: `docs/latex/equations.tex:753`
- C++: No `EqID` owner comment has been attached yet.

```tex
C_{2}
    =\frac{\partial C_{1}}{\partial\eta}
    =- C_{1}^{2}\biggl(
    \bar{m}\frac{-4\eta^{2}+20\eta+8}{\left(1-\eta\right)^{5}}
    +(1-\bar{m})\frac{2\eta^{3}+12\eta^{2}-48\eta+40}{\left[(1-\eta)(2-\eta)\right]^{3}}
    \biggr)
```

##### `c1_xk`
- Label: `eq:c1_xk`
- Source: \cite{Gross2001}, Eq.~(A.41)
- Status: Manual literature match
- Description: Provides a residual Helmholtz-energy relation for dispersion contribution.
- Change note: Mapped manually to the composition derivative of \(C_1\).
- LaTeX: `docs/latex/equations.tex:769`
- C++: No `EqID` owner comment has been attached yet.

```tex
\left(\frac{\partial C_{1}}{\partial x_{k}}\right)_{T,\rho,x_{j\neq k}}
    \equiv C_{1,xk}
    =C_{2}\zeta_{3,xk}-C_{1}^{2}\left\{m_{k}\frac{8\eta-2\eta^{2}}{\left(1-\eta\right)^{4}}-m_{k}\frac{20\eta-27\eta^{2}+12\eta^{3}-2\eta^{4}}{\left[(1-\eta)(2-\eta)\right]^{2}}\right\}
```

#### First Mixed Dispersion Moment

##### `m2epssigma3_bar`
- Label: `eq:m2epssigma3_bar`
- Source: \cite{Gross2001}, Eq.~(A.13)
- Status: Adapted implementation form
- Description: Provides a residual Helmholtz-energy relation for dispersion contribution.
- Change note: Lower similarity; likely algebraically adapted for implementation or combined terms.
- LaTeX: `docs/latex/equations.tex:785`
- C++: No `EqID` owner comment has been attached yet.

```tex
\overline{m^2 \epsilon \sigma^3}=\sum_{i} \sum_{j} x_{i} x_{j} m_{i} m_{j}\left(\frac{\epsilon_{i j}}{k T}\right) \sigma_{i j}^3
```

##### `m2epssigma3_xk`
- Label: `eq:m2epssigma3_xk`
- Source: \cite{Gross2001}, Eq.~(A.39) (appendix sequence inference)
- Status: Manual literature match
- Description: Provides a residual Helmholtz-energy relation for dispersion contribution.
- Change note: Mapped by appendix sequence as the first dispersion composition derivative moment.
- LaTeX: `docs/latex/equations.tex:796`
- C++: No `EqID` owner comment has been attached yet.

```tex
\left(\frac{\partial\overline{m^{2}\epsilon\sigma^{3}}}{\partial x_{k}}\right)_{T,\rho,x_{j\neq k}}
    \equiv \overline{(m^{2}\epsilon\sigma^{3})}_{xk}
    =2m_{k}\sum_{j}x_{j}m_{j}\left(\frac{\epsilon_{kj}}{kT}\right)\sigma_{kj}^{3}
```

#### Second Mixed Dispersion Moment

##### `m2eps2sigma3_bar`
- Label: `eq:m2eps2sigma3_bar`
- Source: \cite{Gross2001}, Eq.~(A.13)
- Status: Adapted implementation form
- Description: Provides a residual Helmholtz-energy relation for dispersion contribution.
- Change note: Lower similarity; likely algebraically adapted for implementation or combined terms.
- LaTeX: `docs/latex/equations.tex:812`
- C++: No `EqID` owner comment has been attached yet.

```tex
\overline{m^2 \epsilon^2 \sigma^3}=\sum_{i} \sum_{j} x_{i} x_{j} m_{i} m_{j}\left(\frac{\epsilon_{i j}}{k T}\right)^2 \sigma_{i j}{ }^3
```

##### `m2eps2sigma3_xk`
- Label: `eq:m2eps2sigma3_xk`
- Source: \cite{Gross2001}, Eq.~(A.40) (appendix sequence inference)
- Status: Manual literature match
- Description: Provides a residual Helmholtz-energy relation for dispersion contribution.
- Change note: Mapped by appendix sequence as the second dispersion composition derivative moment.
- LaTeX: `docs/latex/equations.tex:823`
- C++: No `EqID` owner comment has been attached yet.

```tex
\left(\frac{\partial\overline{m^{2}\epsilon^{2}\sigma^{3}}}{\partial x_{k}}\right)_{T,\rho,x_{j\neq k}}
    \equiv (\overline{m^2\epsilon^2\sigma^3})_{xk}
    =2m_{k}\sum_{j}x_{j}m_{j}\left(\frac{\epsilon_{kj}}{kT}\right)^2\sigma_{kj}^3
```

#### First Packing Polynomial

##### `i1_disp`
- Label: `eq:i1_disp`
- Source: \cite{Gross2001}, Eq.~(A.17)
- Status: Adapted implementation form
- Description: Provides a residual Helmholtz-energy relation for dispersion contribution.
- Change note: Lower similarity; likely algebraically adapted for implementation or combined terms.
- LaTeX: `docs/latex/equations.tex:839`
- C++: No `EqID` owner comment has been attached yet.

```tex
I_{1}(\eta, \bar{m})=\sum_{i=0}^6 a_{i}(\bar{m}) \eta^i
```

##### `deta_i1_deta`
- Label: `eq:deta_i1_deta`
- Source: \cite{Gross2001}, Eq.~(A.30)
- Status: Adapted implementation form
- Description: Provides a differential relation needed for dispersion contribution calculations.
- Change note: Lower similarity; likely algebraically adapted for implementation or combined terms.
- LaTeX: `docs/latex/equations.tex:850`
- C++: No `EqID` owner comment has been attached yet.

```tex
\frac{\partial(\eta I_{1})}{\partial\eta}=\sum_{j=0}^6a_{j}(\bar{m})(j+1)\eta^j
```

##### `i1_xk`
- Label: `eq:i1_xk`
- Source: \cite{Gross2001}, Eq.~(A.42)
- Status: Close literature match
- Description: Provides a residual Helmholtz-energy relation for dispersion contribution.
- Change note: High textual similarity to a tagged equation in the cited local paper export.
- LaTeX: `docs/latex/equations.tex:861`
- C++: No `EqID` owner comment has been attached yet.

```tex
\left(\frac{\partial I_{1}}{\partial x_{k}}\right)_{T,\rho,x_{j\neq k}}
    \equiv I_{1,xk}
    =\sum_{i=0}^{6}[a_{\mathrm{i}}(\bar{m})i\zeta_{3,xk}\eta^{i-1}+a_{\mathrm{i,x}k}\eta^{i}]
```

#### Second Packing Polynomial

##### `i2_disp`
- Label: `eq:i2_disp`
- Source: \cite{Gross2001}, Eq.~(A.17)
- Status: Adapted implementation form
- Description: Provides a residual Helmholtz-energy relation for dispersion contribution.
- Change note: Lower similarity; likely algebraically adapted for implementation or combined terms.
- LaTeX: `docs/latex/equations.tex:877`
- C++: No `EqID` owner comment has been attached yet.

```tex
I_{2}(\eta, \bar{m})=\sum_{i=0}^6 b_{i}(\bar{m}) \eta^i
```

##### `deta_i2_deta`
- Label: `eq:deta_i2_deta`
- Source: \cite{Gross2001}, Eq.~(A.30)
- Status: Adapted implementation form
- Description: Provides a differential relation needed for dispersion contribution calculations.
- Change note: Lower similarity; likely algebraically adapted for implementation or combined terms.
- LaTeX: `docs/latex/equations.tex:888`
- C++: No `EqID` owner comment has been attached yet.

```tex
\frac{\partial(\eta I_{2})}{\partial\eta}=\sum_{i=0}^6b_{j}(\bar{m})(j+1)\eta^i
```

##### `i2_xk`
- Label: `eq:i2_xk`
- Source: \cite{Gross2001}, Eq.~(A.43)
- Status: Close literature match
- Description: Provides a residual Helmholtz-energy relation for dispersion contribution.
- Change note: High textual similarity to a tagged equation in the cited local paper export.
- LaTeX: `docs/latex/equations.tex:899`
- C++: No `EqID` owner comment has been attached yet.

```tex
\left(\frac{\partial I_{2}}{\partial x_{k}}\right)_{T,\rho,x_{j\neq k}}
    \equiv I_{2,xk}
    =\sum_{i=0}^{6}[b_{\mathrm{i}}(\bar{m})i\zeta_{3,xk}\eta^{i-1}+b_{\mathrm{i,xk}}\eta^{i}]
```

#### Dispersion Polynomial Coefficients \(a_i(\bar m)\)

##### `a_i_mbar`
- Label: `eq:a_i_mbar`
- Source: \cite{Gross2001}, Eq.~(A.19)
- Status: Adapted implementation form
- Description: Provides a supporting relation used in dispersion contribution.
- Change note: Lower similarity; likely algebraically adapted for implementation or combined terms.
- LaTeX: `docs/latex/equations.tex:915`
- C++: No `EqID` owner comment has been attached yet.

```tex
a_{i}(\bar{m})=a_{0 i}+\frac{\bar{m}-1}{\bar{m}} a_{1 i}+\frac{\bar{m}-1}{\bar{m}} \frac{\bar{m}-2}{\bar{m}} a_{2 i}
```

##### `a_i_xk`
- Label: `eq:a_i_xk`
- Source: \cite{Gross2001}, Eq.~(A.44)
- Status: Close literature match
- Description: Provides a supporting relation used in dispersion contribution.
- Change note: High textual similarity to a tagged equation in the cited local paper export.
- LaTeX: `docs/latex/equations.tex:926`
- C++: No `EqID` owner comment has been attached yet.

```tex
\left(\frac{\partial a_{i}}{\partial x_{k}}\right)_{T,\rho,x_{j\neq k}}
    \equiv a_{\mathrm{i},xk}
    =\frac{m_{k}}{\bar{m}^{2}}a_{1i}+\frac{m_{k}}{\bar{m}^{2}}\left(3-\frac{4}{\bar{m}}\right)a_{2i}
```

#### Dispersion Polynomial Coefficients \(b_i(\bar m)\)

##### `b_i_mbar`
- Label: `eq:b_i_mbar`
- Source: \cite{Gross2001}, Eq.~(A.19)
- Status: Adapted implementation form
- Description: Provides a supporting relation used in dispersion contribution.
- Change note: Lower similarity; likely algebraically adapted for implementation or combined terms.
- LaTeX: `docs/latex/equations.tex:942`
- C++: No `EqID` owner comment has been attached yet.

```tex
b_{i}(\bar{m})=b_{0 i}+\frac{\bar{m}-1}{\bar{m}} b_{1 i}+\frac{\bar{m}-1}{\bar{m}} \frac{\bar{m}-2}{\bar{m}} b_{2 i}
```

##### `b_i_xk`
- Label: `eq:b_i_xk`
- Source: \cite{Gross2001}, Eq.~(A.45)
- Status: Close literature match
- Description: Provides a supporting relation used in dispersion contribution.
- Change note: High textual similarity to a tagged equation in the cited local paper export.
- LaTeX: `docs/latex/equations.tex:953`
- C++: No `EqID` owner comment has been attached yet.

```tex
\left(\frac{\partial b_{i}}{\partial x_{k}}\right)_{T,\rho,x_{j\neq k}}
    \equiv b_{\mathrm{i},xk}
    =\frac{m_{k}}{\bar{m}^{2}}b_{1i}+\frac{m_{k}}{\bar{m}^{2}}\left(3-\frac{4}{\bar{m}}\right)b_{2i}
```

### Association

#### Association Site Fraction

##### `x_assoc_site`
- Label: `eq:x_assoc_site`
- Source: \cite{Chapman1990} (exact equation number not available in local PDFs)
- Status: No direct numbered source in local corpus
- Description: Provides a residual Helmholtz-energy relation for association contribution.
- Change note: Mass-action site-fraction relation from SAFT association theory; exact numbered equation unavailable in the local PDF set.
- LaTeX: `docs/latex/equations.tex:974`
- C++: No `EqID` owner comment has been attached yet.

```tex
\begin{aligned}
        X^{\mathrm{A}_{i}}
        &=\left[1+ \sum_{j} \sum_{\mathrm{B}_{j}} \rho x_{j}X^{\mathrm{B}_{j}} \Delta^{\mathrm{A}_{i} \mathrm{~B}_{j}}\right]^{-1},
        \\
        &\text{with }\sum_{\mathrm{B}_{j}}\text{ over all sites on molecule }j\ (\mathrm{A}_{j},\mathrm{B}_{j},\mathrm{C}_{j},\ldots),
        \\
        &\text{and }\sum_{j}\text{ over all components.}
    \end{aligned}
```

##### `dx_assoc_drho`
- Label: `eq:dx_assoc_drho`
- Source: \cite{Chapman1990} (exact equation number not available in local PDFs)
- Status: No direct numbered source in local corpus
- Description: Provides a differential relation needed for association contribution calculations.
- Change note: Site-fraction density derivative closure used in association compressibility calculations; exact numbered reference unavailable locally.
- LaTeX: `docs/latex/equations.tex:992`
- C++: No `EqID` owner comment has been attached yet.

```tex
\begin{aligned}
        \left(\frac{\partial X^{A_{i}}}{\partial\rho}\right)_{T,x}
        &=-(X^{A_{i}})^{2}\sum_{j}\sum_{B_{j}}x_{j}
        \Bigg[
        X^{B_{j}}\Delta^{A_{i}B_{j}}
        \\
        &\quad +\rho\Bigg(
        \Delta^{A_{i}B_{j}}\left(\frac{\partial X^{B_{j}}}{\partial\rho}\right)_{T,x}
        +X^{B_{j}}\left(\frac{\partial\Delta^{A_{i}B_{j}}}{\partial\rho}\right)_{T,x}
        \Bigg)
        \Bigg]
    \end{aligned}
```

##### `dx_assoc_dxk`
- Label: `eq:dx_assoc_dxk`
- Source: \cite{Chapman1990} (exact equation number not available in local PDFs)
- Status: No direct numbered source in local corpus
- Description: Provides a differential relation needed for association contribution calculations.
- Change note: Site-fraction composition derivative closure from SAFT association machinery; local progression PDFs do not provide this as a numbered equation.
- LaTeX: `docs/latex/equations.tex:1014`
- C++: No `EqID` owner comment has been attached yet.

```tex
\begin{aligned}
        \left(\frac{\partial X^{A_{i}}}{\partial x_{k}}\right)_{T,v,x_{j\neq k}}
        &=
        -(X^{A_{i}})^{2}\,\rho
        \Bigg[
        \sum_{B_{k}}X^{B_{k}}\Delta^{A_{i}B_{k}}
        \\
        &\quad +
        \sum_{j}x_{j}\sum_{B_{j}}
        \left(
        \Delta^{A_{i}B_{j}}\left(\frac{\partial X^{B_{j}}}{\partial x_{k}}\right)_{T,v,x_{j\neq k}}
        +
        X^{B_{j}}\left(\frac{\partial \Delta^{A_{i}B_{j}}}{\partial x_{k}}\right)_{T,v,x_{j\neq k}}
        \right)
        \Bigg]
    \end{aligned}
```

#### Component Site Density

##### `rho_j_assoc`
- Label: `eq:rho_j_assoc`
- Source: \cite{Chapman1990} (exact equation number not available in local PDFs)
- Status: No direct numbered source in local corpus
- Description: Provides a supporting relation used in association contribution.
- Change note: Component density relation used with association equations; no directly numbered counterpart found in local Chapman-accessible sources.
- LaTeX: `docs/latex/equations.tex:1043`
- C++: No `EqID` owner comment has been attached yet.

```tex
\rho_{j}=x_{j}\rho_{\mathrm{mixture}}
```

#### Association Strength

##### `delta_assoc`
- Label: `eq:delta_assoc`
- Source: \cite{Chapman1990} (exact equation number not available in local PDFs)
- Status: No direct numbered source in local corpus
- Description: Provides a residual Helmholtz-energy relation for association contribution.
- Change note: Association strength definition used in SAFT implementations; exact numbered Chapman equation not available in local PDFs.
- LaTeX: `docs/latex/equations.tex:1057`
- C++: No `EqID` owner comment has been attached yet.

```tex
\Delta^{\mathrm{A}, \mathrm{~B}_{j}}=d_{i j}{ }^3 \mathrm{g}_{i j}^{\mathrm{hs}}\left[\exp \left(\frac{\epsilon^{\mathrm{A}, \mathrm{~B}_{j}}}{k T}\right)-1\right]
```

##### `ddelta_assoc_drho`
- Label: `eq:ddelta_assoc_drho`
- Source: \cite{Chapman1990} (exact equation number not available in local PDFs)
- Status: No direct numbered source in local corpus
- Description: Defines the density differential of the association strength term.
- Change note: Re-homed under the association-strength definition so the \(\Delta^{A_iB_j}\) density derivative stays local to its owning baseline.
- LaTeX: `docs/latex/equations.tex:1068`
- C++: No `EqID` owner comment has been attached yet.

```tex
\left(\frac{\partial\Delta^{A_{i}B_{j}}}{\partial\rho}\right)_{T,x}=d_{ij}^3\kappa^{A_{i}B_{j}}\left[\exp(\epsilon^{A_{i}B_{j}}/kT)-1\right]\left(\frac{\partial g_{ij}^{hs}}{\partial\rho}\right)_{T,x}
```

##### `ddelta_assoc_dxk`
- Label: `eq:ddelta_assoc_dxk`
- Source: \cite{Chapman1990} (exact equation number not available in local PDFs)
- Status: No direct numbered source in local corpus
- Description: Defines the composition differential of the association strength term.
- Change note: Re-homed under the association-strength definition so the \(\Delta^{A_iB_j}\) composition derivative stays local to its owning baseline.
- LaTeX: `docs/latex/equations.tex:1079`
- C++: No `EqID` owner comment has been attached yet.

```tex
\left(\frac{\partial \Delta^{A_{i}B_{j}}}{\partial x_{k}}\right)_{T,v,x_{j\neq k}}
    =
    d_{ij}^{3}\kappa^{A_{i}B_{j}}
    \left[\exp\!\left(\frac{\epsilon^{A_{i}B_{j}}}{kT}\right)-1\right]
    \left(\frac{\partial g_{ij}^{hs}}{\partial x_{k}}\right)_{T,\rho,x_{j\neq k}}
```

### Debye Huckel

#### Screening Parameter \(\kappa\)

##### `kappa_dh`
- Label: `eq:kappa_dh`
- Source: \cite{Bulow2019}, Eq.~(2)
- Status: Adapted implementation form
- Description: Defines the Debye screening quantity used in debye and huckel electrolyte term contribution.
- Change note: Lower similarity; likely algebraically adapted for implementation or combined terms.
- LaTeX: `docs/latex/equations.tex:1102`
- C++: No `EqID` owner comment has been attached yet.

```tex
\kappa=\sqrt{\frac{\rho e^{2}}{k_{B}T\varepsilon_{0}\varepsilon_{r}}\sum_{j}x_{j}z_{j}^{2}}
```

##### `dkappa_dh_drho`
- Label: `eq:dkappa_dh_drho`
- Source: Derived from Eq.~\eqref{eq:kappa_dh}
- Status: Derived relation
- Description: Provides the explicit density derivative of the Debye screening parameter in debye and huckel electrolyte term contribution.
- Change note: Re-homed under the screening-parameter definition so the \(\kappa\) density derivative is owned locally.
- LaTeX: `docs/latex/equations.tex:1113`
- C++: No `EqID` owner comment has been attached yet.

```tex
\begin{aligned}
        \rho\left(\frac{\partial \kappa}{\partial \rho}\right)_{T,x}
        &=
        \frac{\rho}{2}
        \left[
            \frac{\rho e^{2}}{k_{B}T\varepsilon_{0}\varepsilon_{r}}
            \sum_{j}x_{j}z_{j}^{2}
        \right]^{-1/2}
        \Biggl[
            \frac{e^{2}}{k_{B}T\varepsilon_{0}\varepsilon_{r}}\sum_{j}x_{j}z_{j}^{2}
            -
            \frac{\rho e^{2}}{k_{B}T\varepsilon_{0}\varepsilon_{r}^{2}}
            \left(\frac{\partial \varepsilon_{r}}{\partial \rho}\right)_{T,x}
            \sum_{j}x_{j}z_{j}^{2}
        \Biggr]
        \\
        &=
        \frac{\kappa}{2}
        \left[
            1
            -
            \frac{\rho}{\varepsilon_{r}}
            \left(\frac{\partial \varepsilon_{r}}{\partial \rho}\right)_{T,x}
        \right]
        \\
        &=
        \frac{\kappa}{2}
        \qquad
        \text{when }\left(\frac{\partial \varepsilon_{r}}{\partial \rho}\right)_{T,x}=0
    \end{aligned}
```

##### `dkappa_dh_dxi`
- Label: `eq:dkappa_dh_dxi`
- Source: \cite{Bulow2019}, Eq.~(13)
- Status: Manual literature match
- Description: Defines the Debye screening quantity used in debye and huckel electrolyte term contribution.
- Change note: Re-homed under the screening-parameter definition so the \(\kappa\) composition differential is owned locally.
- LaTeX: `docs/latex/equations.tex:1153`
- C++: No `EqID` owner comment has been attached yet.

```tex
\begin{aligned}
        \left(\frac{\partial\kappa}{\partial x_{i}}\right)_{T,v,x_{j\neq i}}
        &=
        \frac12
        \left(\frac{\rho_{N}e^{2}}{k_{B}T\varepsilon_{0}\varepsilon_{r}}\sum_{j}x_{j}z_{j}^{2}\right)^{-\frac12}
        \Biggl\{
            -\frac{\rho_{N}e^2}{k_{B}T\left(\varepsilon_{0}\varepsilon_{r}\right)^2}
            \varepsilon_{0}\frac{\partial\left(\varepsilon_{r}\right)}{\partial x_{i}}
            \sum_{j}x_{j}z_{j}^2
            +
            \frac{\rho_{N}e^2}{k_{B}T\varepsilon_{0}\varepsilon_{r}}
            \alpha_{i}z_{i}^2
        \Biggr\}
        \\
        &=
        \left(\frac{\rho_{N}e^{2}}{k_{B}T}\right)^{\frac{1}{2}}
        \Biggl\{
            -\frac{1}{2}\left(\varepsilon_{0}\varepsilon_{r}\right)^{-\frac{3}{2}}
            \varepsilon_{0}\frac{\partial\left(\varepsilon_{r}\right)}{\partial x_{i}}
            \left[\sum_{j}x_{j}z_{j}^{2}\right]^{\frac{1}{2}}
            \\
            &\qquad\qquad
            +
            \frac{1}{2\sqrt{\varepsilon_{0}\varepsilon_{r}}}
            \Biggl[\sum_{j}x_{j}z_{j}^2\Biggr]^{-\frac{1}{2}}
            \alpha_{i}z_{i}^2
        \Biggr\}
    \end{aligned}
```

#### Finite-Size Correction \(\chi_i\)

##### `chi_dh`
- Label: `eq:chi_dh`
- Source: \cite{Bulow2019}, Eq.~(10)
- Status: Manual literature match
- Description: Defines the Debye screening quantity used in debye and huckel electrolyte term contribution.
- Change note: Mapped manually to the ion-size function definition in the concentration-dependent dielectric derivation.
- LaTeX: `docs/latex/equations.tex:1194`
- C++: No `EqID` owner comment has been attached yet.

```tex
\chi_{i}=\frac{3}{\left(\kappa d_{i}\right)^3}\left[\frac{3}{2}+\ln(1+\kappa d_{i})-2(1+\kappa d_{i})+\frac{1}{2}\left(1+\kappa d_{i}\right)^2\right]
```

##### `sigma_dh`
- Label: `eq:sigma_dh`
- Source: \cite{Cameretti2005}, Eq.~(20)
- Status: Close literature match
- Description: Defines the Debye screening quantity used in debye and huckel electrolyte term contribution.
- Change note: High textual similarity to a tagged equation in the cited local paper export.
- LaTeX: `docs/latex/equations.tex:1205`
- C++: No `EqID` owner comment has been attached yet.

```tex
\sigma_{i}=\left(\frac{\partial(\kappa\chi_{i})}{\partial\kappa}\right)_{T,\mathrm{x}}=-2\chi_{i}+\frac{3}{1+\kappa a_{i}}
```

##### `dchi_dh_drho`
- Label: `eq:dchi_dh_drho`
- Source: Derived from Eq.~\eqref{eq:chi_dh}, Eq.~\eqref{eq:sigma_dh}, and Eq.~\eqref{eq:dkappa_dh_drho}
- Status: Derived relation
- Description: Provides the explicit density derivative of the Debye-Huckel \(\chi_i\) function.
- Change note: Re-homed under the finite-size-correction definition so the \(\chi_i\) density derivative is owned locally.
- LaTeX: `docs/latex/equations.tex:1216`
- C++: No `EqID` owner comment has been attached yet.

```tex
\rho\left(\frac{\partial \chi_i}{\partial \rho}\right)_{T,x}
    =
    \frac{\sigma_i-\chi_i}{2}
```

##### `dchi_dh_dxi`
- Label: `eq:dchi_dh_dxi`
- Source: \cite{Bulow2019}, Eq.~(14)
- Status: Manual literature match
- Description: Defines the Debye screening quantity used in debye and huckel electrolyte term contribution.
- Change note: Re-homed under the finite-size-correction definition so the \(\chi_i\) composition differential is owned locally.
- LaTeX: `docs/latex/equations.tex:1229`
- C++: No `EqID` owner comment has been attached yet.

```tex
\begin{aligned}
        \left(\frac{\partial\chi_{i}}{\partial x_{i}}\right)_{T,v,x_{j\neq i}}
        &=
        -\frac{9}{\left(\kappa d_{i}\right)^{4}}
        \Biggl[
            \frac{3}{2}+\ln(1+\kappa d_{i})-2(1+\kappa d_{i})
            +\frac{1}{2}\left(1+\kappa d_{i}\right)^{2}
        \Biggr]
        d_{i}\frac{\partial\kappa}{\partial x_{i}}
        \\
        &\quad +
        \frac{3}{\left(\kappa d_{i}\right)^{3}}
        \Biggl[
            \frac{1}{1+\kappa d_{i}}-1+\kappa d_{i}
        \Biggr]
        d_{i}\frac{\partial\kappa}{\partial x_{i}}
        \\
        &=
        3\frac{\partial\kappa}{\partial x_{i}}
        \left\{
            -\frac{\chi_{i}}{\kappa}
            +
            \frac{d_{i}}{\left(\kappa d_{i}\right)^{3}}
            \Biggl[
                \frac{1}{1+\kappa d_{i}}-1+\kappa d_{i}
            \Biggr]
        \right\}
    \end{aligned}
```

### Born

#### Born Diameter Rule

##### `d_born_rule`
- Label: `eq:d_born_rule`
- Source: \cite{Bulow2021a}; \cite{Figiel2025}, Eq.~(6)-Eq.~(8) (parameterization choices)
- Status: Project-specific grouped presentation
- Description: Provides the selectable Born-diameter rule used in born electrolyte term contribution.
- Change note: Groups the repeated-left-hand-side Born-diameter setup equations into one visible case set so the selectable diameter definition is documented as one variable with alternative parameterization choices.
- LaTeX: `docs/latex/equations.tex:1275`
- C++: No `EqID` owner comment has been attached yet.

```tex
d^{\text{Born}}_{i}
    =
    \begin{cases}
        d_{i}, & \text{default Born-diameter rule}, \\
        d^{\text{Born, fitted}}_{i}, & \text{fitted Born-diameter rule},
    \end{cases}
    \qquad i\in\mathcal{I}
```

#### Born Mixing Function

##### `f_mix`
- Label: `eq:f_mix`
- Source: \cite{Figiel2025}, Eq.~(10)
- Status: Close literature match
- Description: Provides a supporting relation used in born electrolyte term contribution.
- Change note: High textual similarity to a tagged equation in the cited local paper export.
- LaTeX: `docs/latex/equations.tex:1295`
- C++: No `EqID` owner comment has been attached yet.

```tex
f_{mix}=\sum_{k}x_{k}f_{k}
```

#### Born Shell Shift

##### `delta_d_born`
- Label: `eq:delta_d_born`
- Source: \cite{Figiel2025}, Eq.~(8)
- Status: Close literature match
- Description: Provides a supporting relation used in born electrolyte term contribution.
- Change note: High textual similarity to a tagged equation in the cited local paper export.
- LaTeX: `docs/latex/equations.tex:1309`
- C++: No `EqID` owner comment has been attached yet.

```tex
\Delta d_{i}=\frac{(f_{mix}-1)}{|z_{i}|}\cdot d_{i}^{\mathrm{Born}}
```

##### `ddelta_d_dxi`
- Label: `eq:ddelta_d_dxi`
- Source: \cite{Figiel2025}, Eq.~(6)
- Status: Adapted implementation form
- Description: Provides a differential relation needed for born electrolyte term contribution calculations.
- Change note: Re-homed under the shell-shift definition so the \(\Delta d_i\) composition differential is owned locally.
- LaTeX: `docs/latex/equations.tex:1320`
- C++: No `EqID` owner comment has been attached yet.

```tex
\left(\frac{\partial \Delta d_{j}}{\partial x_{i}}\right)_{T,\rho,x_{k\neq i}}
    =\frac{d_{j}^{\text {Born }}}{\left|z_{j}\right|} \left(\frac{\partial f_{\text {mix }}}{\partial x_{i}}\right)_{T,\rho,x_{k\neq i}}
    =\frac{d_{j}^{\text {Born }}}{\left|z_{j}\right|} f_{i},
```

#### Base Inverse-Diameter Term

##### `dterm_born`
- Label: `eq:D_born`
- Source: \cite{Bulow2021a}, Eq.~(2) (base inverse-diameter term)
- Status: Adapted notation
- Description: Provides the base inverse-diameter term used in the modular Born Helmholtz expression.
- Change note: Names the base Born inverse-diameter contribution explicitly so the modular Born sum can reference a consistent symbol family.
- LaTeX: `docs/latex/equations.tex:1336`
- C++: No `EqID` owner comment has been attached yet.

```tex
D_{i,\mathrm{Born}}^{(\mathrm{bulk})}
=
\frac{1}{d_i^{\mathrm{Born}}},
```

#### Solvation-Shell Inverse-Diameter Term

##### `dterm_ssm`
- Label: `eq:D_ssm`
- Source: \cite{Figiel2025}, Eq.~(7)-Eq.~(8) (derived option term)
- Status: Project-specific modification
- Description: Provides the SSM inverse-diameter correction used in the modular Born Helmholtz expression.
- Change note: Restored as its own named option term because the SSM and DS contributions are distinct quantities rather than alternative definitions of one variable.
- LaTeX: `docs/latex/equations.tex:1352`
- C++: No `EqID` owner comment has been attached yet.

```tex
D_{i,\mathrm{SSM}}^{(\mathrm{bulk})}
=
\frac{1}{d_i^{\mathrm{Born}}+\Delta d_i}
-
\frac{1}{d_i^{\mathrm{Born}}},
```

#### Dielectric-Saturation Inverse-Diameter Term

##### `dterm_ds`
- Label: `eq:D_ds`
- Source: \cite{Figiel2025}, Eq.~(7)-Eq.~(8) (derived option term)
- Status: Project-specific modification
- Description: Provides the DS inverse-diameter correction used in the modular Born Helmholtz expression.
- Change note: Restored as its own named option term because the SSM and DS contributions are distinct quantities rather than alternative definitions of one variable.
- LaTeX: `docs/latex/equations.tex:1370`
- C++: No `EqID` owner comment has been attached yet.

```tex
D_{i,\mathrm{DS}}^{(\mathrm{ion})}
=
\frac{1}{d_i^{\mathrm{Born}}}
-
\frac{1}{d_i^{\mathrm{Born}}+\Delta d_i}.
```

## Residual Helmholtz Energy

### `ares_total`
- Label: `eq:ares_total`
- Source: No explicit citation in equations.tex context
- Status: Project summary relation
- Description: Provides a residual Helmholtz-energy relation for residual helmholz energy.
- Change note: No explicit citation on this equation block in the source file.
- LaTeX: `docs/latex/equations.tex:1390`
- C++: `src/epcsaft/native/epcsaft_ares.cpp:98` (AresContributions ares_contributions_cpp(double t, double rho, const vector<double> &x, const add_args &cppargs) {)

```tex
\tilde{a}^{\mathrm{res}}=\tilde{a}^{h c}+\tilde{a}^{\text {disp }}+\tilde{a}^{\text {assoc }}+\tilde{a}^{\text {DH }}+\tilde{a}^{\text {Born }}
```

### `dares_dT`
- Label: `eq:dares_dT`
- Source: \cite{Gross2001}, Eq.~(A.51)
- Status: Project-specific extension
- Description: Provides the total temperature differential of the residual Helmholtz energy.
- Change note: Moved here from the removed standalone temperature-differential section so one downstream total-\(d\tilde a^\mathrm{res}/dT\) summary relation remains available.
- LaTeX: `docs/latex/equations.tex:1402`
- C++: `src/epcsaft/native/epcsaft_ares.cpp:155` (ScalarContributionTerms temperature_derivative_residual_helmholtz_result_cpp(double t, double rho, vector<double> x, const add_args &cppargs) {)

```tex
\left(\frac{\partial\tilde{a}^\mathrm{res}}{\partial T}\right)_{\rho,x_{i}} =\left(\frac{\partial\tilde{a}^\mathrm{hc}}{\partial T}\right)_{\rho,x_{i}} +\left(\frac{\partial\tilde{a}^\mathrm{disp}}{\partial T}\right)_{\rho,x_{i}}
    +\left(\frac{\partial\tilde{a}^\mathrm{assoc}}{\partial T}\right)_{\rho,x_{i}}
    +\left(\frac{\partial\tilde{a}^\mathrm{DH}}{\partial T}\right)_{\rho,x_{i}}
    +\left(\frac{\partial\tilde{a}^\mathrm{Born}}{\partial T}\right)_{\rho,x_{i}}
```

### Hard-Chain Reference Contribution

#### `ares_hc`
- Label: `eq:ares_hc`
- Source: \cite{Gross2001}, Eq.~(A.4)
- Status: Close literature match
- Description: Provides a residual Helmholtz-energy relation for hard-chain reference contribution.
- Change note: High textual similarity to a tagged equation in the cited local paper export.
- LaTeX: `docs/latex/equations.tex:1418`
- C++: `src/epcsaft/native/epcsaft_ares.cpp:27` (static double ares_hc_cpp(const MixtureState &thermo, const HardChainState &hc_state, const vector<double> &x, const add_args &cppargs) {)

```tex
\begin{aligned}
        \mathmakebox[11.5em][l]{\tilde{a}^{\mathrm{hc}}}
        &=
        \bar{m} \tilde{a}^{\mathrm{hs}}
        -\sum_{i} x_{i}\left(m_{i}-1\right) \ln \mathrm{g}_{i i}^{\mathrm{hs}}\left(\sigma_{i i}\right)
    \end{aligned}
```

#### `hc_ares_dadrho`
- Label: `eq:hc_ares_dadrho`
- Source: \cite{Gross2001}, Eq.~(A.25)
- Status: Close literature match
- Description: Provides the density differential of the hard-chain Helmholtz contribution.
- Change note: Moved back under the hard-chain Helmholtz contribution so Section 4 owns all equations carrying the hard-chain contribution label.
- LaTeX: `docs/latex/equations.tex:1435`
- C++: `src/epcsaft/native/contributions/epcsaft_contrib_hc.cpp:242` (double dadrho_hc_cpp(const MixtureState &thermo, const HardChainState &hc_state, const vector<double> &x, const add_args &cppargs) {)

```tex
\begin{aligned}
        \mathmakebox[11.5em][l]{\rho\left(\frac{\partial\tilde{a}^{hc}}{\partial\rho}\right)_{T,x}}
        &=
        \bar{m}\rho\left(\frac{\partial\tilde{a}^{hs}}{\partial\rho}\right)_{T,x}
        -\sum_{i}x_{\mathrm{i}}(m_{i}-1)(g_{ii}^\mathrm{hs})^{-1}\rho\frac{\partial g_{ii}^\mathrm{hs}}{\partial\rho}
    \end{aligned}
```

#### `hc_ares_dxk`
- Label: `eq:hc_ares_dxk`
- Source: \cite{Gross2001}, Eq.~(A.35) (appendix sequence inference)
- Status: Manual literature match
- Description: Provides the composition differential of the hard-chain Helmholtz contribution.
- Change note: Moved back under the hard-chain Helmholtz contribution so Section 4 owns all equations carrying the hard-chain contribution label.
- LaTeX: `docs/latex/equations.tex:1451`
- C++: `src/epcsaft/native/contributions/epcsaft_contrib_hc.cpp:265` (ContributionDadxResult dadx_hc_cpp(const MixtureState &thermo, const HardChainState &hc_state, double t, double rho, const vector<double> &x, const add_args &cppargs) {)

```tex
\begin{aligned}
        \mathmakebox[11.5em][l]{\left(\frac{\partial\tilde{a}^{\mathrm{hc}}}{\partial x_{k}}\right)_{T,\boldsymbol{\rho},x_{j\neq k}}}
        &=m_{k}\tilde{a}^{\mathrm{hs}}
        +\bar{m}\left(\frac{\partial\tilde{a}^{\mathrm{hs}}}{\partial x_{k}}\right)_{T,\boldsymbol{\rho},x_{j\neq k}}
        \\
        &\quad -\sum_{i}x_{i}(m_{i}-1)(g_{ii}^{\mathrm{hs}})^{-1}
        \left(\frac{\partial g_{ii}^{\mathrm{hs}}}{\partial x_{k}}\right)_{T,\boldsymbol{\rho},x_{j\neq k}}
    \end{aligned}
```

#### `hc_ares_dT`
- Label: `eq:hc_ares_dT`
- Source: \cite{Gross2001}, Eq.~(A.54)
- Status: Manual literature match
- Description: Provides the temperature differential of the hard-chain Helmholtz contribution.
- Change note: Moved back under the hard-chain Helmholtz contribution so Section 4 owns all equations carrying the hard-chain contribution label.
- LaTeX: `docs/latex/equations.tex:1469`
- C++: `src/epcsaft/native/contributions/epcsaft_contrib_hc.cpp:254` (double dadt_hc_cpp(const MixtureState &thermo, const HardChainState &hc_state, const vector<double> &dzeta_dt, const vector<double> &x, const add_args &cppargs) {)

```tex
\begin{aligned}
        \mathmakebox[11.5em][l]{\left(\frac{\partial\tilde{a}^{\mathrm{hc}}}{\partial T}\right)_{\rho,x_{i}}}
        &=
        \bar{m}\left(\frac{\partial\tilde{a}^{\mathrm{hs}}}{\partial T}\right)_{\rho,x_{i}}
        -\sum_{i}x_{\mathrm{i}}(m_{\mathrm{i}}-1)(g_{ii}^{\mathrm{hs}})^{-1}\left(\frac{\partial g_{ii}^{\mathrm{hs}}}{\partial T}\right)_{\rho,x_{i}}
    \end{aligned}
```

### Dispersion Contribution

#### `ares_disp`
- Label: `eq:ares_disp`
- Source: \cite{Gross2001}, Eq.~(A.10)
- Status: Close literature match
- Description: Provides a residual Helmholtz-energy relation for dispersion contribution.
- Change note: High textual similarity to a tagged equation in the cited local paper export.
- LaTeX: `docs/latex/equations.tex:1488`
- C++: `src/epcsaft/native/epcsaft_ares.cpp:37` (static double ares_disp_cpp(const MixtureState &thermo, const DispersionPolynomialState &dispersion) {)

```tex
\tilde{a}^{\mathrm{disp}}=-2 \pi \rho I_{1}(\eta, \bar{m}) \overline{m^2 \epsilon \sigma^3}-\pi \rho \bar{m} C_{1} I_{2}(\eta, \bar{m}) \overline{m^2 \epsilon^2 \sigma^3}
```

#### `disp_ares_dadrho`
- Label: `eq:disp_ares_dadrho`
- Source: \cite{Gross2001}, Eq.~(A.28)
- Status: Close literature match
- Description: Provides the density differential of the dispersion Helmholtz contribution.
- Change note: Moved back under the dispersion Helmholtz contribution so Section 4 owns all equations carrying the dispersion contribution label.
- LaTeX: `docs/latex/equations.tex:1500`
- C++: `src/epcsaft/native/contributions/epcsaft_contrib_disp.cpp:33` (double dadrho_disp_cpp(const MixtureState &thermo, const HardChainState &hc_state, const DispersionPolynomialState &dispersion) {)

```tex
\begin{aligned}
        \rho\left(\frac{\partial\tilde{a}^{disp}}{\partial\rho}\right)_{T,x}
        &=-2\pi\rho\frac{\partial(\eta I_{1})}{\partial\eta}\overline{m^{2}\epsilon\sigma^{3}}
        \\
        &\quad -\pi\rho\bar{m}\left[C_{1}\frac{\partial(\eta I_{2})}{\partial\eta}+C_{2}\eta I_{2}\right]\overline{m^{2}\epsilon^{2}\sigma^{3}}
    \end{aligned}
```

#### `disp_ares_dxk`
- Label: `eq:disp_ares_dxk`
- Source: \cite{Gross2001}, Eq.~(A.38) (appendix sequence inference)
- Status: Manual literature match
- Description: Provides the composition differential of the dispersion Helmholtz contribution.
- Change note: Moved back under the dispersion Helmholtz contribution so Section 4 owns all equations carrying the dispersion contribution label.
- LaTeX: `docs/latex/equations.tex:1516`
- C++: `src/epcsaft/native/contributions/epcsaft_contrib_disp.cpp:48` (ContributionDadxResult dadx_disp_cpp(const MixtureState &thermo, const HardChainState &hc_state, const DispersionPolynomialState &dispersion, double t, double rho, const vector<double> &x, const add_args &cppargs) {)

```tex
\begin{aligned}
        \left(\frac{\partial \tilde{a}^{\mathrm{disp}}}{\partial x_k}\right)_{T,\boldsymbol{\rho},x_{j \neq k}}
        =& -2 \pi \rho\left[I_{1, x k} \overline{m^2 \epsilon \sigma^3}+I_1 \overline{\left(m^2 \epsilon \sigma^3\right)_{x k}}\right] \\
        &- \pi \rho\left\{\left[m_k C_1 I_2+\bar{m} C_{1, x k} I_2+\bar{m} C_1 I_{2, x k}\right] \overline{m^2 \epsilon^2 \sigma^3}+\right.
        \left.\bar{m} C_1 I_2\left(\bar{m}^2 \epsilon^2 \sigma^3\right)_{x k}\right\}
    \end{aligned}
```

#### `disp_ares_dT`
- Label: `eq:disp_ares_dT`
- Source: \cite{Gross2001}, Eq.~(A.58)
- Status: Manual literature match
- Description: Provides the temperature differential of the dispersion Helmholtz contribution.
- Change note: Moved back under the dispersion Helmholtz contribution so Section 4 owns all equations carrying the dispersion contribution label.
- LaTeX: `docs/latex/equations.tex:1532`
- C++: `src/epcsaft/native/contributions/epcsaft_contrib_disp.cpp:39` (double dadt_disp_cpp(const MixtureState &thermo, double deta_dt, double t, const DispersionPolynomialState &dispersion) {)

```tex
\left(\frac{\partial\tilde{a}^{disp}}{\partial T}\right)_{\rho,x_{i}}
    =
    -2\pi\rho\left(I_{1,T}\overline{m^{2}\epsilon\sigma^{3}}+I_{1}\overline{m^{2}\epsilon\sigma^{3}}_{,T}\right)
    -\pi\rho\bar{m}\left(C_{1,T}I_{2}+C_{1}I_{2,T}\right)\overline{m^{2}\epsilon^{2}\sigma^{3}}
    -\pi\rho\bar{m}C_{1}I_{2}\overline{m^{2}\epsilon^{2}\sigma^{3}}_{,T}
```

### Association Contribution

#### `ares_assoc`
- Label: `eq:ares_assoc`
- Source: \cite{Chapman1990} (exact equation number not available in local PDFs)
- Status: No direct numbered source in local corpus
- Description: Provides a residual Helmholtz-energy relation for association contribution.
- Change note: Association Helmholtz form is traced to Chapman/Wertheim SAFT association theory, but the exact numbered equation is not present in the local progression PDFs.
- LaTeX: `docs/latex/equations.tex:1550`
- C++: `src/epcsaft/native/epcsaft_ares.cpp:43` (static double ares_assoc_cpp(const AssociationIntermediateState &assoc_state, const vector<double> &x) {)

```tex
\tilde{a}^{\mathrm{assoc}}= \sum_{i} x_{i}\sum_{\mathrm{A}_{i}}\left(\ln X^{\mathrm{A}_{i}}-\frac{X^{\mathrm{A}_{i}}}{2}+\frac{1}{2}\right)
```

#### `assoc_ares_dadrho`
- Label: `eq:assoc_ares_dadrho`
- Source: \cite{Chapman1990} (exact equation number not available in local PDFs)
- Status: No direct numbered source in local corpus
- Description: Provides the density differential of the association Helmholtz contribution.
- Change note: Moved back under the association Helmholtz contribution so Section 4 owns all equations carrying the association contribution label.
- LaTeX: `docs/latex/equations.tex:1562`
- C++: `src/epcsaft/native/contributions/epcsaft_contrib_assoc.cpp:316` (double dadrho_assoc_cpp()

```tex
\rho\left(\frac{\partial\tilde{a}^{assoc}}{\partial\rho}\right)_{T,x}
    =\rho\sum_{i}x_{i}\sum_{A_{i}}\left(\frac{1}{X^{A_{i}}}-\frac{1}{2}\right)\left(\frac{\partial X^{A_{i}}}{\partial\rho}\right)_{T,x}.
```

#### `assoc_ares_dxk`
- Label: `eq:assoc_ares_dxk`
- Source: \cite{Chapman1990} (exact equation number not available in local PDFs)
- Status: No direct numbered source in local corpus
- Description: Provides the composition differential of the association Helmholtz contribution.
- Change note: Moved back under the association Helmholtz contribution so Section 4 owns all equations carrying the association contribution label.
- LaTeX: `docs/latex/equations.tex:1574`
- C++: `src/epcsaft/native/contributions/epcsaft_contrib_assoc.cpp:378` (ContributionDadxResult dadx_assoc_cpp(const MixtureState &thermo, const HardChainState &hc_state, const AssociationIntermediateState &assoc_state, double t, double rho, const vector<double> &x, const add_args &cppargs) {)

```tex
\left(\frac{\partial \tilde a^{assoc}}{\partial x_{k}}\right)_{T,v,x_{i\neq k}}
    =
    \sum_{A_{k}}\left(\ln X^{A_{k}}-\frac{X^{A_{k}}}{2}+\frac{1}{2}\right)
    +
    \sum_{i} x_{i}\sum_{A_{i}}
    \left(\frac{1}{X^{A_{i}}}-\frac{1}{2}\right)
    \left(\frac{\partial X^{A_{i}}}{\partial x_{k}}\right)_{T,v,x_{i\neq k}}
```

#### `assoc_ares_dT`
- Label: `eq:assoc_ares_dT`
- Source: Native implementation owner \texttt{dadt\_assoc\_cpp(...)}
- Status: Implementation-owned equation
- Description: Provides the analytical temperature differential of the association Helmholtz contribution.
- Change note: Moved back under the association Helmholtz contribution so Section 4 owns all equations carrying the association contribution label.
- LaTeX: `docs/latex/equations.tex:1591`
- C++: `src/epcsaft/native/contributions/epcsaft_contrib_assoc.cpp:366` (double dadt_assoc_cpp(const AssociationIntermediateState &assoc_state, const vector<double> &x) {)

```tex
\left(\frac{\partial\tilde{a}^{\mathrm{assoc}}}{\partial T}\right)_{\rho,x_i}
    =
    \sum_i x_i \sum_{A_i}
    \left(\frac{1}{X^{A_i}}-\frac{1}{2}\right)
    \left(\frac{\partial X^{A_i}}{\partial T}\right)_{\rho,x_i}
```

### Debye and Huckel Contribution

#### `ares_dh`
- Label: `eq:ares_dh`
- Source: \cite{Bulow2019}, Eq.~(2)
- Status: Adapted implementation form
- Description: Defines the Debye screening quantity used in debye and huckel electrolyte term contribution.
- Change note: Lower similarity; likely algebraically adapted for implementation or combined terms.
- LaTeX: `docs/latex/equations.tex:1609`
- C++: `src/epcsaft/native/epcsaft_ares.cpp:56` (static double ares_ion_cpp(double t, const IonIntermediateState &ion_state) {)

```tex
\tilde{a}^{DH}=-\frac{\kappa e^{2}}{12\pi\varepsilon_{0}\varepsilon_{r}k_{B}T}\sum_{i}x_{i}z_{i}^{2}\chi_{i}
```

#### `dh_ares_dadrho`
- Label: `eq:dh_ares_dadrho`
- Source: Derived from Eq.~\eqref{eq:ares_dh}, Eq.~\eqref{eq:dkappa_dh_drho}, and Eq.~\eqref{eq:dchi_dh_drho}
- Status: Derived relation
- Description: Provides the explicit density differential of the Debye-Huckel Helmholtz contribution in chain-rule form.
- Change note: Moved back under the Debye-Huckel Helmholtz contribution so Section 4 owns all equations carrying the Debye-Huckel contribution label.
- LaTeX: `docs/latex/equations.tex:1621`
- C++: No `EqID` owner comment has been attached yet.

```tex
\rho\left(\frac{\partial\tilde{a}^{DH}}{\partial\rho}\right)_{T,x}
    =
    -\frac{e^2}{12\pi\varepsilon_0\varepsilon_r k_B T}
    \left[
        \rho\left(\frac{\partial \kappa}{\partial \rho}\right)_{T,x}\sum_i x_i z_i^2 \chi_i
        +
        \kappa\sum_i x_i z_i^2 \rho\left(\frac{\partial \chi_i}{\partial \rho}\right)_{T,x}
    \right]
```

#### `dh_ares_dxi`
- Label: `eq:dh_ares_dxi`
- Source: \cite{Bulow2019} (equation number unresolved in local corpus)
- Status: Literature update
- Description: Defines the Debye screening quantity used in debye and huckel electrolyte term contribution.
- Change note: Moved back under the Debye-Huckel Helmholtz contribution so Section 4 owns all equations carrying the Debye-Huckel contribution label.
- LaTeX: `docs/latex/equations.tex:1639`
- C++: `src/epcsaft/native/contributions/epcsaft_contrib_ion.cpp:96` (ContributionDadxResult dadx_ion_cpp(const MixtureState &thermo, const IonIntermediateState &ion_state, double t, double rho, const vector<double> &x, const add_args &cppargs) {)

```tex
\begin{aligned}
        \left(\frac{\partial \tilde{a}^{DH}}{\partial x_{i}}\right)_{T,v,x_{j\neq i}}
        =
        -\frac{e^{2}}{12\pi\varepsilon_{0}k_{B}T}
        \Bigg[
          & \frac{1}{\varepsilon_{r}}\left(\frac{\partial\kappa}{\partial x_{i}}\right)_{T,v,x_{j\neq i}}
            \sum_{j}x_{j}z_{j}^{2}\chi_{j}
        \\
        - & \frac{\kappa}{\varepsilon_{r}^{2}}
            \left(\frac{\partial\varepsilon_{r}}{\partial x_{i}}\right)_{T,v,x_{j\neq i}}
            \sum_{j}x_{j}z_{j}^{2}\chi_{j}
        \\
        + & \frac{\kappa}{\varepsilon_{r}}
            \sum_{j}\chi_{j}z_{j}^{2}
            \left(\frac{\partial x_{j}}{\partial x_{i}}\right)_{T,v,x_{k\neq i}}
        \\
        + & \frac{\kappa}{\varepsilon_{r}}
            \sum_{j}x_{j}z_{j}^{2}
            \left(\frac{\partial\chi_{j}}{\partial x_{i}}\right)_{T,v,x_{k\neq i}}
            \Bigg]
    \end{aligned}
```

#### `dh_ares_dT`
- Label: `eq:dh_ares_dT`
- Source: \cite{Cameretti2005}, Eq.~(13)-Eq.~(14) (derived temperature differential)
- Status: Derived from literature equation
- Description: Provides the temperature differential of the Debye-Huckel Helmholtz contribution.
- Change note: Moved back under the Debye-Huckel Helmholtz contribution so Section 4 owns all equations carrying the Debye-Huckel contribution label.
- LaTeX: `docs/latex/equations.tex:1670`
- C++: `src/epcsaft/native/contributions/epcsaft_contrib_ion.cpp:79` (double dadt_ion_cpp(const IonIntermediateState &ion_state, double t, const vector<double> &x, const add_args &cppargs) {)

```tex
\left(\frac{\partial \tilde{a}^{DH}}{\partial T}\right)_{\rho,x_i}
    \approx
    -\frac{\tilde{a}^{DH}}{T}
```

### Born Contribution

#### `ares_born`
- Label: `eq:ares_born_modular`
- Source: \cite{Bulow2021a}, Eq.~(2)
- Status: Adapted notation
- Description: Provides a residual Helmholtz-energy relation for born electrolyte term contribution.
- Change note: Moderate-to-high similarity; notation/arrangement appears adapted from the cited equation.
- LaTeX: `docs/latex/equations.tex:1686`
- C++: `src/epcsaft/native/epcsaft_ares.cpp:65` (static double ares_born_cpp(double t, const BornIntermediateState &born_state) {)

```tex
\tilde{a}^{\mathrm{Born}}
=
-\frac{e^2}{4\pi\varepsilon_0 k_{\mathrm{B}} T}
\sum_i x_i z_i^2
\sum_{\delta \in \mathcal{D}}
\left(
1-\frac{1}{\varepsilon_{r,m_\delta}}
\right)
D_{i,\delta}^{(m_\delta)},
```

#### `born_ares_dadrho`
- Label: `eq:born_ares_dadrho`
- Source: \cite{Bulow2021a}, Eq.~(4)
- Status: Manual literature match
- Description: Provides the density differential of the Born Helmholtz contribution.
- Change note: Moved back under the Born Helmholtz contribution so Section 4 owns all equations carrying the Born contribution label.
- LaTeX: `docs/latex/equations.tex:1706`
- C++: `src/epcsaft/native/contributions/epcsaft_contrib_born.cpp:194` (double dadrho_born_cpp() {)

```tex
\rho\left(\frac{\partial\tilde{a}^{Born}}{\partial\rho}\right)_{T,x}= 0
```

#### `born_ares_dxi`
- Label: `eq:born_ares_dxi`
- Source: \cite{Bulow2021a}, Eq.~(3)
- Status: Project-specific modification
- Description: Provides a differential relation needed for born electrolyte term contribution calculations.
- Change note: Moved back under the Born Helmholtz contribution so Section 4 owns all equations carrying the Born contribution label.
- LaTeX: `docs/latex/equations.tex:1717`
- C++: No `EqID` owner comment has been attached yet.

```tex
\begin{aligned}
        \left(\frac{\partial \tilde{a}^{Born}}{\partial x_{i}}\right)_{T,v_{N},x_{j\neq i}}
        &=
        -\frac{e^{2}}{4\pi k_{B}T\varepsilon_{0}}
        \Biggl[
            \left(1-\frac{1}{\varepsilon_{r}}\right)\frac{z_{i}^{2}}{d^{\text{Born}}_{i}}
            \\
            &\qquad\qquad
            +
            \left(\frac{1}{\varepsilon_{r}^{2}}\right)\sum_{j} \frac{x_{j} z_{j}^2}{d^{\text{Born}}_{j}}
            \left(\frac{\partial\varepsilon_{r}}{\partial x_{i}}\right)_{T,v,x_{j\neq i}}
        \Biggr]
    \end{aligned}
```

#### `born_ares_ssmds_dxi`
- Label: `eq:born_ares_ssmds_dxi`
- Source: \cite{Figiel2025} (equation number unresolved in local corpus)
- Status: Project-specific modification
- Description: Provides a differential relation needed for born electrolyte term contribution calculations.
- Change note: Moved back under the Born Helmholtz contribution so Section 4 owns all equations carrying the Born contribution label.
- LaTeX: `docs/latex/equations.tex:1740`
- C++: No `EqID` owner comment has been attached yet.

```tex
\begin{aligned}
        \left(\frac{\partial \tilde{a}^{\mathrm{Born}}}{\partial x_{i}}\right)_{T,v,x_{j\neq i}}
         & =
        -\frac{e^2}{4\pi\varepsilon_{0} k_{B} T}
        \Bigg[
            z_{i}^2\!\left(
            \left(1-\frac{1}{\varepsilon_{r,\mathrm{bulk}}}\right)
            \left[D_{i,\mathrm{Born}}^{(\mathrm{bulk})}+D_{i,\mathrm{SSM}}^{(\mathrm{bulk})}\right]
            +
            \left(1-\frac{1}{\varepsilon_{r,\mathrm{ion}}}\right)D_{i,\mathrm{DS}}^{(\mathrm{ion})}
            \right)
        \\
         & \qquad\qquad
            +\sum_{j} x_{j} z_{j}^2
            \Bigg(
            \frac{1}{\varepsilon_{r,\mathrm{bulk}}^2}
            \frac{\partial \varepsilon_{r,\mathrm{bulk}}}{\partial x_{i}}
            \left[D_{j,\mathrm{Born}}^{(\mathrm{bulk})}+D_{j,\mathrm{SSM}}^{(\mathrm{bulk})}\right]
            \\
         & \qquad\qquad\qquad\qquad
            +
            \frac{1}{\varepsilon_{r,\mathrm{ion}}^2}
            \frac{\partial \varepsilon_{r,\mathrm{ion}}}{\partial x_{i}}
            D_{j,\mathrm{DS}}^{(\mathrm{ion})}
            \\
         & \qquad\qquad\qquad\qquad
            +
            \left(\frac{1}{\varepsilon_{r,\mathrm{ion}}}-\frac{1}{\varepsilon_{r,\mathrm{bulk}}}\right)
            \frac{1}{\left(d_{j}^{\mathrm{Born}}+\Delta d_{j}\right)^2}
            \frac{\partial \Delta d_{j}}{\partial x_{i}}
            \Bigg)
            \Bigg],
    \end{aligned}
```

#### `born_ares_dT`
- Label: `eq:born_ares_dT`
- Source: \cite{Bulow2021a}, Eq.~(2) (derived temperature differential)
- Status: Derived from literature equation
- Description: Provides the temperature differential of the Born Helmholtz contribution.
- Change note: Moved back under the Born Helmholtz contribution so Section 4 owns all equations carrying the Born contribution label.
- LaTeX: `docs/latex/equations.tex:1783`
- C++: `src/epcsaft/native/contributions/epcsaft_contrib_born.cpp:199` (double dadt_born_cpp(double t, const BornIntermediateState &born_state) {)

```tex
\begin{aligned}
        \left(\frac{\partial\tilde{a}^{\mathrm{Born}}}{\partial T}\right)_{\rho,x_i}
        &=
        -\frac{\tilde{a}^{\mathrm{Born}}}{T}
        \\
        &\quad -\frac{e^2}{4\pi\varepsilon_0 k_{\mathrm{B}} T}
        \sum_i x_i z_i^2
        \sum_{\delta \in \mathcal{D}}
        \left[
            \frac{D_{i,\delta}^{(m_\delta)}}{\varepsilon_{r,m_\delta}^{2}}
            \left(\frac{\partial \varepsilon_{r,m_\delta}}{\partial T}\right)_{\rho,x}
            +
            \left(1-\frac{1}{\varepsilon_{r,m_\delta}}\right)
            \left(\frac{\partial D_{i,\delta}^{(m_\delta)}}{\partial T}\right)_{\rho,x}
        \right]
    \end{aligned}
```

#### `born_mode_set`
- Label: `eq:born_mode_set`
- Source: \cite{Figiel2025}, Eq.~(7) (option-driven reformulation)
- Status: Project-specific modification
- Description: Defines the set of active Born subterms used in the modular Helmholtz expression.
- Change note: Documents the current option-driven Born-term set directly rather than the legacy version naming.
- LaTeX: `docs/latex/equations.tex:1812`
- C++: No `EqID` owner comment has been attached yet.

```tex
\mathcal{D}
=
\{\mathrm{Born}\}
\cup
\mathcal{D}_{\mathrm{add}},
\qquad
\mathcal{D}_{\mathrm{add}}
\subseteq
\{\mathrm{SSM},\mathrm{DS}\},
```

#### `born_mode_medium`
- Label: `eq:born_mode_medium`
- Source: \cite{Bulow2021a}, Eq.~(2); \cite{Figiel2025}, Eq.~(7) (option-driven reformulation)
- Status: Project-specific modification
- Description: Defines the medium assignment used by each active Born subterm in the modular Helmholtz expression.
- Change note: Makes the bulk-vs-ion medium choice explicit for the modular Born documentation.
- LaTeX: `docs/latex/equations.tex:1831`
- C++: No `EqID` owner comment has been attached yet.

```tex
m_\delta
    =
    \begin{cases}
        \mathrm{bulk}, & \delta\in\{\mathrm{Born},\mathrm{SSM}\}, \\
        \mathrm{ion},  & \delta=\mathrm{DS}.
    \end{cases}
```

## Compressibility Factor

### `z_total`
- Label: `eq:z_total`
- Source: \cite{Gross2001}, Eq.~(A.24)
- Status: Project-specific extension
- Description: Provides a supporting relation used in pressure (compressibility factor).
- Change note: Mapped to Eq. (A.24) and extended here by adding association, Debye-Huckel, and Born compressibility contributions.
- LaTeX: `docs/latex/equations.tex:1860`
- C++: `src/epcsaft/native/epcsaft_Z.cpp:63` (CompressibilityFactorResult compressibility_factor_result_cpp(double t, double rho, vector<double> x, const add_args &cppargs) {)

```tex
Z = 1
    + \rho\left(\frac{\partial\tilde{a}^{hc}}{\partial\rho}\right)_{T,x}
    + \rho\left(\frac{\partial\tilde{a}^{disp}}{\partial\rho}\right)_{T,x}
    + \rho\left(\frac{\partial\tilde{a}^{assoc}}{\partial\rho}\right)_{T,x}
    + \rho\left(\frac{\partial\tilde{a}^{DH}}{\partial\rho}\right)_{T,x}
    + \rho\left(\frac{\partial\tilde{a}^{Born}}{\partial\rho}\right)_{T,x}
```

### `z_alpha`
- Label: `eq:z_alpha`
- Source: Project bridge identity based on the contribution-resolved \(Z\) split
- Status: Project-specific organization
- Description: Gives the generic contribution compressibility bridge identity in compressibility factor.
- Change note: Replaces the five redundant contribution-specific bridge equations with one generic \(Z^\alpha\) form covering \(\alpha\in\{hc,disp,assoc,DH,Born\}\).
- LaTeX: `docs/latex/equations.tex:1876`
- C++: No `EqID` owner comment has been attached yet.

```tex
Z^{\alpha}
    =
    \rho\left(\frac{\partial\tilde{a}^{\alpha}}{\partial\rho}\right)_{T,x}
    ,
    \qquad
    \alpha\in\{hc,disp,assoc,DH,Born\}
```

## Chemical Potential

### `mu_res`
- Label: `eq:mu_res`
- Source: \cite{Gross2001}, Eq.~(A.1)
- Status: Adapted implementation form
- Description: Gives a chemical-potential contribution in chemical potential.
- Change note: Lower similarity; likely algebraically adapted for implementation or combined terms.
- LaTeX: `docs/latex/equations.tex:1897`
- C++: `src/epcsaft/native/epcsaft_mu.cpp:57` (ResidualChemicalPotentialResult residual_chemical_potential_result_cpp(double t, double rho, vector<double> x, const add_args &cppargs) {)

```tex
\tilde{\mu_{k}}^{\mathrm{res}}=\frac{\mu_{k}^\mathrm{res}}{kT} = \frac{\hat{\mu}_{k}^\mathrm{res}}{RT}
```

### `mu_res_from_ares`
- Label: `eq:mu_res_from_ares`
- Source: \cite{Gross2001}, Eq.~(A.33)
- Status: Close literature match
- Description: Gives a chemical-potential contribution in chemical potential.
- Change note: High textual similarity to a tagged equation in the cited local paper export.
- LaTeX: `docs/latex/equations.tex:1908`
- C++: `src/epcsaft/native/epcsaft_mu.cpp:6` (static vector<double> mu_contribution_cpp()

```tex
\begin{aligned}
        \tilde{\mu_{k}}^{\mathrm{res}}
        &=\tilde{a}^{\mathrm{res}}+(\mathrm{Z}-1)
        +\left(\frac{\partial\tilde{a}^{\mathrm{res}}}{\partial x_{k}}\right)_{T,\nu,x_{i\neq k}}
        \\
        &\quad -\sum_{j=1}^{N}\left[x_{j}\left(\frac{\partial\tilde{a}^{\mathrm{res}}}{\partial x_{j}}\right)_{T,\nu,x_{i\neq j}}\right]
    \end{aligned}
```

### `mu_res_sum`
- Label: `eq:mu_res_sum`
- Source: \cite{Gross2001}, Eq.~(A.33)
- Status: Manual literature match
- Description: Gives the explicit residual chemical-potential decomposition in chemical potential.
- Change note: Written in explicit non-summation form to match the style used for the other property families in this section.
- LaTeX: `docs/latex/equations.tex:1925`
- C++: `src/epcsaft/native/epcsaft_mu.cpp:57` (ResidualChemicalPotentialResult residual_chemical_potential_result_cpp(double t, double rho, vector<double> x, const add_args &cppargs) {)

```tex
\tilde{\mu_{k}}^{\mathrm{res}}=\tilde{\mu_{k}}^{\mathrm{hc}}+\tilde{\mu_{k}}^{\mathrm{disp}}+\tilde{\mu_{k}}^{\mathrm{assoc}}+\tilde{\mu_{k}}^{\mathrm{DH}}+\tilde{\mu_{k}}^{\mathrm{Born}}
```

### `mu_alpha`
- Label: `eq:mu_alpha`
- Source: \cite{Gross2001}, Eq.~(A.33)
- Status: Project-specific organization
- Description: Gives the generic contribution chemical-potential relation in chemical potential.
- Change note: Replaces the five redundant contribution-specific chemical-potential equations with one generic \(\alpha\)-form covering \(\alpha\in\{hc,disp,assoc,DH,Born\}\).
- LaTeX: `docs/latex/equations.tex:1936`
- C++: No `EqID` owner comment has been attached yet.

```tex
\begin{aligned}
        \tilde{\mu}^{\alpha}_{k}
        &=
        \tilde{a}^{\alpha}
        + Z^{\alpha}
        +\left(\frac{\partial\tilde{a}^{\alpha}}{\partial x_{k}}\right)_{T,\nu,x_{i\neq k}}
        \\
        &\quad
        -\sum_{j=1}^{N}\left[x_{j}\left(\frac{\partial\tilde{a}^{\alpha}}{\partial x_{j}}\right)_{T,\nu,x_{i\neq j}}\right],
        \qquad
        \alpha\in\{hc,disp,assoc,DH,Born\}
    \end{aligned}
```

## Fugacity Coefficient

### `lnphi_total`
- Label: `eq:lnphi_total`
- Source: \cite{Gross2001}, Eq.~(A.32)
- Status: Manual literature match
- Description: Gives the total fugacity-coefficient relation in fugacity coefficient.
- Change note: Mapped manually to the residual-chemical-potential form used in the PC-SAFT appendix.
- LaTeX: `docs/latex/equations.tex:1962`
- C++: `src/epcsaft/native/epcsaft_fugcoef.cpp:88` (FugacityContributionResult fugacity_coefficient_result_cpp(double t, double rho, vector<double> x, const add_args &cppargs) {)

```tex
\ln\varphi_{k}
    =
    \tilde{\mu}_{k}^{\mathrm{res}}
    -
    \ln Z
    =
    \frac{\mu_{k}^{\mathrm{res}}}{kT}
    -
    \ln Z
```

### `lnphi_total_sum`
- Label: `eq:lnphi_total_sum`
- Source: Project decomposition based on Eq.~\eqref{eq:lnphi_total}
- Status: Project-specific organization
- Description: Gives the explicit fugacity-coefficient decomposition in fugacity coefficient.
- Change note: Written in explicit non-summation form to match the contribution-by-contribution structure used throughout this section.
- LaTeX: `docs/latex/equations.tex:1983`
- C++: `src/epcsaft/native/epcsaft_fugcoef.cpp:88` (FugacityContributionResult fugacity_coefficient_result_cpp(double t, double rho, vector<double> x, const add_args &cppargs) {)

```tex
\ln\varphi_{k}
    =
    \ln\varphi_{k}^{hc}
    +
    \ln\varphi_{k}^{disp}
    +
    \ln\varphi_{k}^{assoc}
    +
    \ln\varphi_{k}^{DH}
    +
    \ln\varphi_{k}^{Born}
```

### `lnphi_alpha`
- Label: `eq:lnphi_alpha`
- Source: Project decomposition based on Eq.~\eqref{eq:lnphi_total}
- Status: Project-specific organization
- Description: Gives the generic contribution fugacity-coefficient relation in fugacity coefficient.
- Change note: Uses only the explicit $Z^\alpha$ allocation requested for the contribution-resolved fugacity-coefficient terms.
- LaTeX: `docs/latex/equations.tex:2004`
- C++: `src/epcsaft/native/epcsaft_fugcoef.cpp:38` (static vector<double> lnfug_contribution_cpp()

```tex
\ln\varphi_{k}^{\alpha}
    =
    \tilde{\mu}_{k}^{\alpha}
    -
    \frac{Z^{\alpha}}{Z-1}\ln Z,
    \qquad
    \alpha\in\{hc,disp,assoc,DH,Born\}
```

### `lnphi_alpha_near_ideal`
- Label: `eq:lnphi_alpha_near_ideal`
- Source: Derived from Eq.~\eqref{eq:lnphi_alpha}
- Status: Derived approximation
- Description: Gives the near-ideal approximation for a contribution fugacity coefficient in fugacity coefficient.
- Change note: Retained explicitly as an approximation only, using the $Z\rightarrow 1$ limit requested for documentation.
- LaTeX: `docs/latex/equations.tex:2021`
- C++: `src/epcsaft/native/epcsaft_fugcoef.cpp:14` (static double stable_logz_over_zminus1(double Z) {)

```tex
\ln\varphi_{k}^{\alpha}
    \approx
    \tilde{\mu}_{k}^{\alpha}
    -
    Z^{\alpha},
    \qquad
    \alpha\in\{hc,disp,assoc,DH,Born\}
```

## Activity Coefficient

### Symmetric Reference

#### `gamma_sym`
- Label: `eq:gamma_sym`
- Source: Standard thermodynamic definition
- Status: Project-specific organization
- Description: Gives the symmetric-reference activity-coefficient definition in activity coefficient.
- Change note: Added explicitly in non-log form so the activity section mirrors the fugacity-coefficient organization.
- LaTeX: `docs/latex/equations.tex:2047`
- C++: `src/epcsaft/native/epcsaft_activity.cpp:60` (vector<double> miac_gamma_vector_cpp()

```tex
\gamma_{i}
    =
    \frac{\varphi_{i}(T,p,\mathbf{x})}{\varphi_{0i}(T,p,x_{i}=1)}
```

#### `lngamma_sym`
- Label: `eq:lngamma_sym`
- Source: Derived from Eq.~\eqref{eq:gamma_sym}
- Status: Derived relation
- Description: Gives the symmetric-reference logarithmic activity-coefficient definition in activity coefficient.
- Change note: Expressed only in terms of fugacity coefficients, as requested for the activity section.
- LaTeX: `docs/latex/equations.tex:2060`
- C++: `src/epcsaft/native/epcsaft_activity.cpp:60` (vector<double> miac_gamma_vector_cpp()

```tex
\ln\gamma_{i}
    =
    \ln\varphi_{i}
    -
    \ln\varphi_{0i}
```

### Infinite-Dilution Reference

#### `gamma_asym_inf`
- Label: `eq:gamma_asym_inf`
- Source: Standard thermodynamic definition
- Status: Project-specific organization
- Description: Gives the infinite-dilution activity-coefficient definition in activity coefficient.
- Change note: Added explicitly in non-log form for ions and solutes referenced to infinite dilution.
- LaTeX: `docs/latex/equations.tex:2077`
- C++: `src/epcsaft/native/epcsaft_activity.cpp:456` (ActivityCoefficientNative activity_coefficient_values_cpp()

```tex
\gamma_{i}^{*}
    =
    \frac{\varphi_{i}(T,p,\mathbf{x})}{\varphi_{i}^{\infty}(T,p,x_{i}\to 0)}
```

#### `lngamma_asym_inf`
- Label: `eq:lngamma_asym_inf`
- Source: Derived from Eq.~\eqref{eq:gamma_asym_inf}
- Status: Derived relation
- Description: Gives the infinite-dilution logarithmic activity-coefficient definition in activity coefficient.
- Change note: Expressed only in terms of fugacity coefficients, as requested for the activity section.
- LaTeX: `docs/latex/equations.tex:2090`
- C++: `src/epcsaft/native/epcsaft_activity.cpp:456` (ActivityCoefficientNative activity_coefficient_values_cpp()

```tex
\ln\gamma_{i}^{*}
    =
    \ln\varphi_{i}
    -
    \ln\varphi_{i}^{\infty}
```

#### `lngamma_asym_sum`
- Label: `eq:lngamma_asym_sum`
- Source: Project decomposition based on Eq.~\eqref{eq:lngamma_asym_inf}
- Status: Project-specific organization
- Description: Gives the explicit infinite-dilution activity-coefficient decomposition in activity coefficient.
- Change note: Written in explicit non-summation form to mirror the fugacity-coefficient decomposition.
- LaTeX: `docs/latex/equations.tex:2107`
- C++: No `EqID` owner comment has been attached yet.

```tex
\ln\gamma_{k}^{*}
    =
    \ln\gamma_{k}^{hc,*}
    +
    \ln\gamma_{k}^{disp,*}
    +
    \ln\gamma_{k}^{assoc,*}
    +
    \ln\gamma_{k}^{DH,*}
    +
    \ln\gamma_{k}^{Born,*}
```

#### `gamma_alpha_asym`
- Label: `eq:gamma_alpha_asym`
- Source: Project decomposition based on Eq.~\eqref{eq:gamma_asym_inf}
- Status: Project-specific organization
- Description: Gives the generic contribution activity-coefficient definition in activity coefficient.
- Change note: Defined only from contribution fugacity coefficients, not from chemical-potential expressions.
- LaTeX: `docs/latex/equations.tex:2128`
- C++: No `EqID` owner comment has been attached yet.

```tex
\gamma_{k}^{\alpha,*}
    =
    \frac{\varphi_{k}^{\alpha}}{\varphi_{k}^{\alpha,\infty}}
```

#### `lngamma_alpha_asym`
- Label: `eq:lngamma_alpha_asym`
- Source: Derived from Eq.~\eqref{eq:gamma_alpha_asym}
- Status: Derived relation
- Description: Gives the generic logarithmic contribution activity-coefficient definition in activity coefficient.
- Change note: Written only in terms of contribution fugacity coefficients, consistent with the requested section logic.
- LaTeX: `docs/latex/equations.tex:2141`
- C++: No `EqID` owner comment has been attached yet.

```tex
\ln\gamma_{k}^{\alpha,*}
    =
    \ln\varphi_{k}^{\alpha}
    -
    \ln\varphi_{k}^{\alpha,\infty}
```

### Mean Ionic Activity Coefficient

#### `gamma_pm_asym`
- Label: `eq:gamma_pm_asym`
- Source: Standard thermodynamic definition
- Status: Project-specific organization
- Description: Gives the mean-ionic infinite-dilution activity coefficient in activity coefficient.
- Change note: Added in non-log form to parallel the logarithmic mean-ionic relation already used in electrolyte work.
- LaTeX: `docs/latex/equations.tex:2158`
- C++: No `EqID` owner comment has been attached yet.

```tex
\gamma_{\pm}^{*}
    =
    \left(
    \left(\gamma_{c}^{*}\right)^{\nu_{c}}
    \left(\gamma_{a}^{*}\right)^{\nu_{a}}
    \right)^{\frac{1}{\nu_{c}+\nu_{a}}}
```

#### `lngamma_pm_asym`
- Label: `eq:lngamma_pm_asym`
- Source: Derived from Eq.~\eqref{eq:gamma_pm_asym}
- Status: Derived relation
- Description: Gives the logarithmic mean-ionic infinite-dilution activity coefficient in activity coefficient.
- Change note: Written explicitly in the standard stoichiometric-weighted logarithmic form.
- LaTeX: `docs/latex/equations.tex:2174`
- C++: No `EqID` owner comment has been attached yet.

```tex
\ln\gamma_{\pm}^{*}
    =
    \frac{\nu_{c}\ln\gamma_{c}^{*}+\nu_{a}\ln\gamma_{a}^{*}}{\nu_{c}+\nu_{a}}
```

#### `lngamma_pm_alpha_asym`
- Label: `eq:lngamma_pm_alpha_asym`
- Source: Project decomposition based on Eq.~\eqref{eq:lngamma_pm_asym}
- Status: Project-specific organization
- Description: Gives the logarithmic contribution mean-ionic activity coefficient in activity coefficient.
- Change note: Contribution form written directly from the contribution activity-coefficient terms.
- LaTeX: `docs/latex/equations.tex:2187`
- C++: No `EqID` owner comment has been attached yet.

```tex
\ln\gamma_{\pm}^{\alpha,*}
    =
    \frac{\nu_{c}\ln\gamma_{c}^{\alpha,*}+\nu_{a}\ln\gamma_{a}^{\alpha,*}}{\nu_{c}+\nu_{a}}
```

## Enthalpy, Entropy, and Gibbs Free Energy

### Entropy

#### `h_res`
- Label: `eq:h_res`
- Source: \cite{Figiel2025}, Eq.~(13)
- Status: Adapted implementation form
- Description: Provides a differential relation needed for entropy calculations.
- Change note: Lower similarity; likely algebraically adapted for implementation or combined terms.
- LaTeX: `docs/latex/equations.tex:2208`
- C++: `src/epcsaft/native/epcsaft_thermo.cpp:4` (double hres_cpp(double t, double rho, vector<double> x, const add_args &cppargs) {)

```tex
\tilde{h}^{\mathrm{res}} = \frac{\hat{h}^{\mathrm{res}}}{RT}=-T\left(\frac{\partial\tilde{a}^{\mathrm{res}}}{\partial T}\right)_{\rho,x_{i}}+(Z-1)
```

#### `s_res_from_s_vol`
- Label: `eq:s_res_from_s_vol`
- Source: \cite{Gross2001}, Eq.~(A.47)
- Status: Manual literature match
- Description: Provides a supporting relation used in entropy.
- Change note: Mapped manually to the residual-entropy relation with logarithmic compressibility correction.
- LaTeX: `docs/latex/equations.tex:2221`
- C++: No `EqID` owner comment has been attached yet.

```tex
\tilde{s}^{\mathrm{res}} = \frac{\hat{s}^{\mathrm{res}}(P,T)}{R}=\frac{\hat{s}^{\mathrm{res}}(\nu,T)}{R}+\ln(Z)
```

#### `s_res`
- Label: `eq:s_res`
- Source: \cite{Gross2001}, Eq.~(A.48)
- Status: Manual literature match
- Description: Provides a differential relation needed for entropy calculations.
- Change note: Mapped manually to the residual-entropy temperature-derivative form.
- LaTeX: `docs/latex/equations.tex:2232`
- C++: `src/epcsaft/native/epcsaft_thermo.cpp:18` (double sres_cpp(double t, double rho, vector<double> x, const add_args &cppargs) {)

```tex
\tilde{s}^{\mathrm{res}} = \frac{\hat{s}^{\mathrm{res}}(P,T)}{R}=-T\left[\left(\frac{\partial\tilde{a}^{\mathrm{res}}}{\partial T}\right)_{\rho,x_{i}}+\frac{\tilde{a}^{\mathrm{res}}}{T}\right]+\ln(Z)
```

### Gibbs Free Energy

#### `g_res_from_hs`
- Label: `eq:g_res_from_hs`
- Source: \cite{Gross2001}, Eq.~(A.49)
- Status: Manual literature match
- Description: Provides a supporting relation used in gibbs free energy.
- Change note: Mapped manually to the residual Gibbs relation via enthalpy and entropy.
- LaTeX: `docs/latex/equations.tex:2245`
- C++: No `EqID` owner comment has been attached yet.

```tex
\tilde{g}^{\mathrm{res}}=\frac{\hat{g}^{\mathrm{res}}}{RT}=\frac{\hat{h}^{\mathrm{res}}}{RT}-\frac{\hat{s}^{\mathrm{res}}(P,T)}{R}
```

#### `g_res_from_ares`
- Label: `eq:g_res_from_ares`
- Source: \cite{Gross2001}, Eq.~(A.50)
- Status: Manual literature match
- Description: Provides a residual Helmholtz-energy relation for gibbs free energy.
- Change note: Mapped manually to the residual Gibbs relation in Helmholtz/compressibility form.
- LaTeX: `docs/latex/equations.tex:2256`
- C++: `src/epcsaft/native/epcsaft_thermo.cpp:11` (double gres_cpp(double t, double rho, vector<double> x, const add_args &cppargs) {)

```tex
\tilde{g}^{\mathrm{res}}=\tilde{a}^{\mathrm{res}}+(Z-1)-\ln(Z)
```

#### `delta_g_solv_inf_x`
- Label: `eq:delta_g_solv_inf_x`
- Source: \cite{Figiel2025}, Eq.~(14)
- Status: Adapted implementation form
- Description: Provides a residual Helmholtz-energy relation for gibbs free energy.
- Change note: Lower similarity; likely algebraically adapted for implementation or combined terms.
- LaTeX: `docs/latex/equations.tex:2269`
- C++: `src/epcsaft/native/epcsaft_activity.cpp:104` (vector<double> gsolv_values_cpp()

```tex
\Delta^{\mathrm{solv},x}G_{i}^{\infty}(T,p,x_{j\neq i},x_{i}\to0)=RT\ln(\varphi_{i}^{\infty}(T,p,x_{j\neq i},x_{i}\to0))
```

#### `delta_g_transfer_inf`
- Label: `eq:delta_g_transfer_inf`
- Source: \cite{Figiel2025}, Eq.~(14)
- Status: Adapted implementation form
- Description: Provides a residual Helmholtz-energy relation for gibbs free energy.
- Change note: Lower similarity; likely algebraically adapted for implementation or combined terms.
- LaTeX: `docs/latex/equations.tex:2281`
- C++: `src/epcsaft/native/epcsaft_activity.cpp:104` (vector<double> gsolv_values_cpp()

```tex
\Delta^{\mathrm{tr}}G_{i}^{\infty,\mathrm{S}1\to\mathrm{S}2}(T,p,x_{j\neq i},x_{i}\rightarrow0)=RT\ln\left(\frac{\varphi_{i}^{\infty,\mathrm{S}2}}{\varphi_{i}^{\infty,\mathrm{S}1}}\right)
```

## Conversions

### Mean Ionic Conversions

#### `mu_pm_charge`
- Label: `eq:mu_pm_charge`
- Source: \cite{Ascani2022}, Eq.~(16)
- Status: Close literature match
- Description: Gives the mean-ionic chemical-potential conversion from cation and anion chemical potentials.
- Change note: Added to gather the general mean-ionic property conversions in one dedicated section.
- LaTeX: `docs/latex/equations.tex:2298`
- C++: No `EqID` owner comment has been attached yet.

```tex
\mu_{\pm}
    =
    \frac{\dfrac{1}{|z_{\mathrm{cat}}|}\mu_{\mathrm{cat}}+\dfrac{1}{|z_{\mathrm{an}}|}\mu_{\mathrm{an}}}
    {\dfrac{1}{|z_{\mathrm{cat}}|}+\dfrac{1}{|z_{\mathrm{an}}|}}
```

#### `f_pm_charge`
- Label: `eq:f_pm_charge`
- Source: \cite{Ascani2022}, Eq.~(17)
- Status: Close literature match
- Description: Gives the mean-ionic fugacity conversion from cation and anion fugacities.
- Change note: Added to gather the general mean-ionic property conversions in one dedicated section.
- LaTeX: `docs/latex/equations.tex:2312`
- C++: No `EqID` owner comment has been attached yet.

```tex
f_{\pm}
    =
    \left[
        \left(f_{\mathrm{cat}}\right)^{1/|z_{\mathrm{cat}}|}
        \left(f_{\mathrm{an}}\right)^{1/|z_{\mathrm{an}}|}
    \right]^{\!\left[\left(1/|z_{\mathrm{cat}}|\right)+\left(1/|z_{\mathrm{an}}|\right)\right]^{-1}}
```

#### `a_pm_charge`
- Label: `eq:a_pm_charge`
- Source: \cite{Ascani2022}, Eq.~(18)
- Status: Close literature match
- Description: Gives the mean-ionic activity conversion from cation and anion activities.
- Change note: Added to gather the general mean-ionic property conversions in one dedicated section.
- LaTeX: `docs/latex/equations.tex:2328`
- C++: No `EqID` owner comment has been attached yet.

```tex
a_{\pm}
    =
    \left[
        \left(a_{\mathrm{cat}}\right)^{1/|z_{\mathrm{cat}}|}
        \left(a_{\mathrm{an}}\right)^{1/|z_{\mathrm{an}}|}
    \right]^{\!\left[\left(1/|z_{\mathrm{cat}}|\right)+\left(1/|z_{\mathrm{an}}|\right)\right]^{-1}}
```

#### `phi_pm_charge`
- Label: `eq:phi_pm_charge`
- Source: \cite{Ascani2022}, Eq.~(19)
- Status: Close literature match
- Description: Gives the mean-ionic fugacity-coefficient conversion from cation and anion fugacity coefficients.
- Change note: Added to gather the general mean-ionic property conversions in one dedicated section.
- LaTeX: `docs/latex/equations.tex:2344`
- C++: No `EqID` owner comment has been attached yet.

```tex
\varphi_{\pm}
    =
    \left[
        \left(\varphi_{\mathrm{cat}}\right)^{1/|z_{\mathrm{cat}}|}
        \left(\varphi_{\mathrm{an}}\right)^{1/|z_{\mathrm{an}}|}
    \right]^{\!\left[\left(1/|z_{\mathrm{cat}}|\right)+\left(1/|z_{\mathrm{an}}|\right)\right]^{-1}}
```

#### `gamma_pm_charge`
- Label: `eq:gamma_pm_charge`
- Source: \cite{Ascani2022}, Eq.~(20)
- Status: Close literature match
- Description: Gives the mean-ionic activity-coefficient conversion from cation and anion activity coefficients.
- Change note: Added to gather the general mean-ionic property conversions in one dedicated section.
- LaTeX: `docs/latex/equations.tex:2360`
- C++: No `EqID` owner comment has been attached yet.

```tex
\gamma_{\pm}
    =
    \left[
        \left(\gamma_{\mathrm{cat}}\right)^{1/|z_{\mathrm{cat}}|}
        \left(\gamma_{\mathrm{an}}\right)^{1/|z_{\mathrm{an}}|}
    \right]^{\!\left[\left(1/|z_{\mathrm{cat}}|\right)+\left(1/|z_{\mathrm{an}}|\right)\right]^{-1}}
```

### Mole-Fraction and Molality Basis

#### `gamma_pm_molality`
- Label: `eq:gamma_pm_molality`
- Source: \cite{Figiel2025}, Eq.~(22)
- Status: Close literature match
- Description: Gives the mean-ionic activity-coefficient conversion from mole-fraction to molality basis.
- Change note: Kept as the primary mole-fraction to molality conversion for mean-ionic activity coefficients.
- LaTeX: `docs/latex/equations.tex:2380`
- C++: No `EqID` owner comment has been attached yet.

```tex
\gamma_{\pm}^{*,m}
    =
    \frac{\gamma_{\pm}^{*,x}}{1+M_{\mathrm{solvent}}\cdot\tilde{m}_{\mathrm{solute}}\cdot\sum_{i}\nu_{i,\mathrm{ion}}}
```

#### `lngamma_pm_molality`
- Label: `eq:lngamma_pm_molality`
- Source: Derived from Eq.~\eqref{eq:gamma_pm_molality}
- Status: Derived relation
- Description: Gives the logarithmic mean-ionic activity-coefficient conversion from mole-fraction to molality basis.
- Change note: Added as the direct log-form companion to the existing salt-basis \(\gamma_{\pm}\) conversion.
- LaTeX: `docs/latex/equations.tex:2393`
- C++: No `EqID` owner comment has been attached yet.

```tex
\ln\gamma_{\pm}^{*,m}
    =
    \ln\gamma_{\pm}^{*,x}
    -
    \ln\left(1+M_{\mathrm{solvent}}\cdot\tilde{m}_{\mathrm{solute}}\cdot\sum_{i}\nu_{i,\mathrm{ion}}\right)
```

## Supplemental

### Legacy / Paper-Direct Expressions

#### Debye-Huckel Density Differential

##### `dadrho_dh_explicit`
- Label: `eq:dadrho_dh_explicit`
- Source: \cite{Cameretti2005}, Eq.~(10)
- Status: Original literature form
- Description: Gives the compact paper-direct Debye-Huckel density differential retained for traceability.
- Change note: Moved out of the main density-differential section because the active document now shows the explicit \(\kappa\)- and \(\chi_i\)-based chain-rule construction there.
- LaTeX: `docs/latex/equations.tex:2417`
- C++: `src/epcsaft/native/contributions/epcsaft_contrib_ion.cpp:71` (double dadrho_ion_cpp(double t, const IonIntermediateState &ion_state) {)

```tex
\rho\left(\frac{\partial\tilde{a}^{DH}}{\partial\rho}\right)_{T,x}
    =-\frac{\kappa e^2}{24\pi kT\epsilon}\sum_{i}x_{i}z_{i}{}^{2}\sigma_{i}
```

#### Cameretti 2005 Chemical-Potential Forms

##### `mu_dh_2005`
- Label: `eq:mu_dh_2005`
- Source: \cite{Cameretti2005}, Eq.~(11)
- Status: Original literature form
- Description: Defines the Debye screening quantity used in debye and huckel electrolyte term contribution.
- Change note: Baseline 2005 formulation retained for comparison to the reorganized active chemical-potential presentation.
- LaTeX: `docs/latex/equations.tex:2431`
- C++: No `EqID` owner comment has been attached yet.

```tex
\tilde{\mu}^{DH}_{i}= -\frac{q_{i}^2 \kappa}{24 \pi k T \epsilon}\left[2 \chi_{i}+\frac{\sum_{j} x_{j} q_{j}^2 \sigma_{k}}{\sum_{j} x_{j} q_{j}^2}\right]
```

##### `sigma_dh_2005`
- Label: `eq:sigma_dh_2005`
- Source: \cite{Cameretti2005}, Eq.~(20)
- Status: Close literature match
- Description: Defines the Debye screening quantity used in debye and huckel electrolyte term contribution.
- Change note: Retained beside Eq.~\eqref{eq:mu_dh_2005} so the original 2005 chemical-potential form stays self-contained in the legacy comparison block.
- LaTeX: `docs/latex/equations.tex:2442`
- C++: No `EqID` owner comment has been attached yet.

```tex
\sigma_{k}=\left(\frac{\partial\left(\kappa \chi_{k}\right)}{\partial \kappa}\right)_{T, \mathrm{~N}}=-2 \chi_{k}+\frac{3}{1+\kappa d_{k}}
```

### Bjerrum Treatment

#### State and Residual Helmholtz Relations

##### `r_ion_bjerrum`
- Label: `eq:r_ion_bjerrum`
- Source: \cite{Bulow2021}, Eq.~(16)
- Status: Manual literature match
- Description: Provides a supporting relation used in debye and huckel electrolyte term contribution.
- Change note: Mapped manually to the closest-approach radius selection between ion diameter and Bjerrum length.
- LaTeX: `docs/latex/equations.tex:2459`
- C++: No `EqID` owner comment has been attached yet.

```tex
R_{i}=
    \begin{cases}
        d_{\mathrm{ion},i}, & d_{\mathrm{ion},i}>l_{B}, \\
        l_{B},              & d_{\mathrm{ion},i}<l_{B}.
    \end{cases}
```

##### `bjerrum_length`
- Label: `eq:bjerrum_length`
- Source: \cite{Bulow2021}, Eq.~(10)
- Status: Close literature match
- Description: Specifies dielectric-property mixing or derivative form for debye and huckel electrolyte term contribution.
- Change note: High textual similarity to a tagged equation in the cited local paper export.
- LaTeX: `docs/latex/equations.tex:2474`
- C++: No `EqID` owner comment has been attached yet.

```tex
l_{B}=\frac{\left|z_{i}z_{j}\right|e^2}{8\pi\varepsilon_{0}\varepsilon_{r}k_{B}T}
```

##### `ares_dh_bjerrum`
- Label: `eq:ares_dh_bjerrum`
- Source: \cite{Bulow2021}, Eq.~(20)
- Status: Manual literature match
- Description: Defines the Debye screening quantity used in debye and huckel electrolyte term contribution.
- Change note: Mapped manually to the reduced Debye-Huckel Helmholtz form with dissociation degree factors.
- LaTeX: `docs/latex/equations.tex:2485`
- C++: No `EqID` owner comment has been attached yet.

```tex
\tilde{a}^{DH}=-\frac{\kappa e^{2}}{12\pi\varepsilon_{0}\varepsilon_{r}k_{B}T}\sum_{i}\alpha_{i}x_{i}z_{i}^{2}\chi_{i}
```

##### `kappa_dh_bjerrum`
- Label: `eq:kappa_dh_bjerrum`
- Source: \cite{Bulow2021}, Eq.~(18)
- Status: Manual literature match
- Description: Defines the Debye screening quantity used in debye and huckel electrolyte term contribution.
- Change note: Mapped manually to the Bjerrum-treatment Debye screening parameter definition.
- LaTeX: `docs/latex/equations.tex:2496`
- C++: No `EqID` owner comment has been attached yet.

```tex
\kappa=\sqrt{\frac{\rho e^{2}}{k_{B}T\varepsilon_{0}\varepsilon_{r}}\sum_{j}\alpha_{j}x_{j}z_{j}^{2}}
```

##### `chi_dh_bjerrum`
- Label: `eq:chi_dh_bjerrum`
- Source: \cite{Bulow2021}, Eq.~(19)
- Status: Manual literature match
- Description: Defines the Debye screening quantity used in debye and huckel electrolyte term contribution.
- Change note: Mapped manually to the Bjerrum-treatment chi-function definition.
- LaTeX: `docs/latex/equations.tex:2507`
- C++: No `EqID` owner comment has been attached yet.

```tex
\chi_{i}=\frac{3}{\left(\kappa R_{i}\right)^3}\left[\frac{3}{2}+\ln(1+\kappa R_{i})-2(1+\kappa R_{i})+\frac{1}{2}\left(1+\kappa R_{i}\right)^2\right]
```

#### Density Differential

##### `dadrho_dh_bjerrum`
- Label: `eq:dadrho_dh_bjerrum`
- Source: \cite{Bulow2021}, Eq.~(18)
- Status: Adapted implementation form
- Description: Defines the Debye screening quantity used in debye and huckel electrolyte term contribution.
- Change note: Grouped under the supplemental Bjerrum-treatment section so the main density-differential section stays focused on the base Debye-Huckel formulation.
- LaTeX: `docs/latex/equations.tex:2520`
- C++: No `EqID` owner comment has been attached yet.

```tex
\rho\left(\frac{\partial\tilde{a}^{DH}}{\partial\rho}\right)_{T,x}
    =-\frac{\kappa e^2}{24\pi kT\epsilon}\sum_{i}\alpha _{i}x_{i}z_{i}{}^{2}\sigma_{i}
```

##### `sigma_dh_bjerrum`
- Label: `eq:sigma_dh_bjerrum`
- Source: \cite{Bulow2021}, Eq.~(26)
- Status: Adapted implementation form
- Description: Defines the Debye screening quantity used in debye and huckel electrolyte term contribution.
- Change note: Grouped under the supplemental Bjerrum-treatment section so the Bjerrum-specific sigma relation stays beside the other extended Debye-Huckel helpers.
- LaTeX: `docs/latex/equations.tex:2532`
- C++: No `EqID` owner comment has been attached yet.

```tex
\sigma_{i}=\left(\frac{\partial(\kappa\chi_{i})}{\partial\kappa}\right)_{T,\mathrm{x}}=-2\chi_{i}+\frac{3}{1+\kappa R_{i}}
```

#### Composition Differential

##### `dares_dh_dxi_bjerrum`
- Label: `eq:dares_dh_dxi_bjerrum`
- Source: \cite{Bulow2021}, Eq.~(21)
- Status: Manual literature match
- Description: Defines the Debye screening quantity used in debye and huckel electrolyte term contribution.
- Change note: Grouped under the supplemental Bjerrum-treatment section so the composition-derivative variant stays with the rest of the extended Bjerrum relations.
- LaTeX: `docs/latex/equations.tex:2545`
- C++: No `EqID` owner comment has been attached yet.

```tex
\begin{aligned}
        \left(\frac{\partial \tilde{a}^{DH}}{\partial x_{i}}\right)_{T,v,x_{j\neq i}}
        =
        -\frac{e^{2}}{12\pi\varepsilon_{0}k_{B}T}
        \Bigg[
          & \frac{1}{\varepsilon_{r}}\left(\frac{\partial\kappa}{\partial x_{i}}\right)_{T,v,x_{j\neq i}}
            \sum_{j}\alpha_{j}x_{j}z_{j}^{2}\chi_{j}
        \\
        - & \frac{\kappa}{\varepsilon_{r}^{2}}
            \left(\frac{\partial\varepsilon_{r}}{\partial x_{i}}\right)_{T,v,x_{j\neq i}}
            \sum_{j}\alpha_{j}x_{j}z_{j}^{2}\chi_{j}
        \\
        + & \frac{\kappa}{\varepsilon_{r}}\alpha_{i}z_{i}^{2}\chi_{i}
        \\
        + & \frac{\kappa}{\varepsilon_{r}}
            \sum_{j}\alpha_{j}x_{j}z_{j}^{2}
            \left(\frac{\partial\chi_{j}}{\partial x_{i}}\right)_{T,v,x_{k\neq i}}
            \Bigg]
    \end{aligned}
```

##### `dkappa_dh_dxi_bjerrum`
- Label: `eq:dkappa_dh_dxi_bjerrum`
- Source: \cite{Bulow2021}, Eq.~(22)
- Status: Adapted notation
- Description: Defines the Debye screening quantity used in debye and huckel electrolyte term contribution.
- Change note: Grouped under the supplemental Bjerrum-treatment section so the auxiliary composition derivative for \kappa stays adjacent to the Bjerrum-specific Helmholtz relation.
- LaTeX: `docs/latex/equations.tex:2574`
- C++: No `EqID` owner comment has been attached yet.

```tex
\begin{aligned}
         & \left(\frac{\partial\kappa}{\partial x_{i}}\right)_{T,v,x_{j\neq i}}=\frac12\left(\frac{\rho_{N}e^{2}}{k_{B}T\varepsilon_{0}\varepsilon_{r}}\sum_{j}\alpha_{j}x_{j}z_{j}^{2}\right)^{-\frac12} \\&\left\{-\frac{\rho_{N}e^2}{k_{B}T\left(\varepsilon_{0}\varepsilon_{r}\right)^2}\varepsilon_{0}\frac{\partial\left(\varepsilon_{r}\right)}{\partial x_{i}}\sum_{j}\alpha_{j}x_{j}z_{j}^2+\frac{\rho_{N}e^2}{k_{B}T\varepsilon_{0}\varepsilon_{r}}\alpha_{i}z_{i}^2\right\}\\&=\left(\frac{\rho_{N}e^{2}}{k_{B}T}\right)^{\frac{1}{2}}\left\{-\frac{1}{2}\left(\varepsilon_{0}\varepsilon_{r}\right)^{-\frac{3}{2}}\varepsilon_{0}\frac{\partial\left(\varepsilon_{r}\right)}{\partial x_{i}}\left[\sum_{j}\alpha_{j}x_{j}z_{j}^{2}\right]^{\frac{1}{2}}\right.\\&+\frac{1}{2\sqrt{\varepsilon_{0}\varepsilon_{r}}}\Bigg[\sum_{j}\alpha_{j}x_{j}z_{j}^2\Bigg]^{-\frac{1}{2}}\alpha_{i}z_{i}^2\Bigg\}
    \end{aligned}
```

##### `dchi_dh_dxi_bjerrum`
- Label: `eq:dchi_dh_dxi_bjerrum`
- Source: \cite{Bulow2021}, Eq.~(23)
- Status: Manual literature match
- Description: Defines the Debye screening quantity used in debye and huckel electrolyte term contribution.
- Change note: Grouped under the supplemental Bjerrum-treatment section so the Bjerrum-specific chi derivative remains adjacent to the other extended Debye-Huckel derivatives.
- LaTeX: `docs/latex/equations.tex:2587`
- C++: No `EqID` owner comment has been attached yet.

```tex
\begin{aligned}
        \left(\frac{\partial\chi_{i}}{\partial x_{i}}\right)_{T,v,x_{j\neq i}} & =-\frac{9}{\left(\kappa R_{i}\right)^{4}}\biggl[\frac{3}{2}+\ln(1+\kappa R_{i})-2(1+\kappa R_{i}) \\&+\frac{1}{2}\left(1+\kappa R_{i}\right)^{2}\bigg]R_{i}\frac{\partial\kappa}{\partial x_{i}}+\frac{3}{\left(\kappa R_{i}\right)^{3}}\bigg[\frac{1}{1+\kappa R_{i}}-1+\kappa R_{i}\bigg]\\R_{i}\frac{\partial\kappa}{\partial x_{i}}&= 3\frac{\partial\kappa}{\partial x_{i}}\left\{-\frac{\chi_{i}}{\kappa}+\frac{R_{i}}{\left(\kappa R_{i}\right)^{3}}\biggl[\frac{1}{1+\kappa R_{i}}-1+\kappa R_{i}\biggr]\right\}
    \end{aligned}
```

#### Ion-Pairing Relations

##### `alpha_ion_pair`
- Label: `eq:alpha_ion_pair`
- Source: \cite{Bulow2021}, Eq.~(14)
- Status: Close literature match
- Description: Gives an activity or fugacity-coefficient relation in debye and huckel electrolyte term contribution.
- Change note: High textual similarity to a tagged equation in the cited local paper export.
- LaTeX: `docs/latex/equations.tex:2602`
- C++: No `EqID` owner comment has been attached yet.

```tex
\alpha=\frac{-1+\sqrt{1+4x_{\pm} K_{ip}\frac{\left(\gamma_{\pm}^*\left(x_{f}\right)\right)^2}{\gamma_{ip}^*}}}{2x_{\pm} K_{ip}\frac{\left(\gamma_{\pm}^*\left(x_{f}\right)\right)^2}{\gamma_{ip}^*}}
```

##### `k_ion_pair`
- Label: `eq:k_ion_pair`
- Source: \cite{Bulow2021}, Eq.~(9) and Eq.~(11)
- Status: Manual literature match
- Description: Specifies dielectric-property mixing or derivative form for debye and huckel electrolyte term contribution.
- Change note: This line combines the algebraic ion-pair equilibrium form and the configurational integral expression shown as separate equations in the paper.
- LaTeX: `docs/latex/equations.tex:2613`
- C++: No `EqID` owner comment has been attached yet.

```tex
K_{ip}(T)=\frac{(1-\alpha)\cdot\gamma_{ip}^{*}}{\alpha^{2}\cdot x_{\pm}\cdot\left(\gamma_{\pm}^{*}(x_{f}(\alpha))\right)^{2}} = 4\pi\rho_{N}\int_{a}^{l_{B}}\exp\left(\frac{\left|z_{i}z_{j}\right|e^2}{4\pi\varepsilon_{0}\varepsilon_{r}k_{B}T}\cdot\frac{1}{r}\right)r^2dr
```

##### `gamma_ion_pair_unity`
- Label: `eq:gamma_ion_pair_unity`
- Source: \cite{Bulow2021} (approximation not explicitly numbered)
- Status: Project-specific approximation
- Description: Provides a supporting relation used in debye and huckel electrolyte term contribution.
- Change note: Assumption $\gamma_{ip}^* \approx 1$ is used as a simplifying closure and is not a standalone numbered equation in the source paper.
- LaTeX: `docs/latex/equations.tex:2624`
- C++: No `EqID` owner comment has been attached yet.

```tex
\gamma_{ip}^* \approx 1
```

##### `x_pm_balance`
- Label: `eq:x_pm_balance`
- Source: \cite{Bulow2021}, Eq.~(15)
- Status: Adapted notation
- Description: Provides a supporting relation used in debye and huckel electrolyte term contribution.
- Change note: Moderate-to-high similarity; notation/arrangement appears adapted from the cited equation.
- LaTeX: `docs/latex/equations.tex:2635`
- C++: No `EqID` owner comment has been attached yet.

```tex
x_{\pm} = x_{f} + x_{ip}
```

##### `x_free_ion_from_alpha`
- Label: `eq:x_free_ion_from_alpha`
- Source: \cite{Bulow2021}, Eq.~(15)
- Status: Adapted implementation form
- Description: Provides a supporting relation used in debye and huckel electrolyte term contribution.
- Change note: Lower similarity; likely algebraically adapted for implementation or combined terms.
- LaTeX: `docs/latex/equations.tex:2646`
- C++: No `EqID` owner comment has been attached yet.

```tex
x_{\pm} = \alpha x_{f}
```

##### `x_ion_pair_from_alpha`
- Label: `eq:x_ion_pair_from_alpha`
- Source: \cite{Bulow2021}, Eq.~(9) (derived stoichiometric form)
- Status: Derived from literature equation
- Description: Provides a supporting relation used in debye and huckel electrolyte term contribution.
- Change note: Stoichiometric rearrangement used with Eq. (9) during alpha-based ion-pair splitting.
- LaTeX: `docs/latex/equations.tex:2657`
- C++: No `EqID` owner comment has been attached yet.

```tex
x_{\pm} = (1 - \alpha) x_{ip}
```

##### `x_pm_stoich`
- Label: `eq:x_pm_stoich`
- Source: \cite{Bulow2021}, Eq.~(7) (adapted variable form)
- Status: Derived from literature equation
- Description: Provides a supporting relation used in debye and huckel electrolyte term contribution.
- Change note: Geometric mean form mirrors the mean-ionic expression pattern and is written here for mole fractions.
- LaTeX: `docs/latex/equations.tex:2668`
- C++: No `EqID` owner comment has been attached yet.

```tex
x_{\pm}=(x_{c}^{\nu c}\cdot x_{a}^{\nu a})^{\frac{1}{\nu_{c}+\nu_{a}}}
```

##### `mu_dh_infinite_dilution`
- Label: `eq:mu_dh_infinite_dilution`
- Source: \cite{Bulow2021}, Eq.~(25)
- Status: Adapted implementation form
- Description: Defines the Debye screening quantity used in debye and huckel electrolyte term contribution.
- Change note: Lower similarity; likely algebraically adapted for implementation or combined terms.
- LaTeX: `docs/latex/equations.tex:2679`
- C++: No `EqID` owner comment has been attached yet.

```tex
\mu_{i}^{DH}\left(x_{f}\right)=\left(\frac{\partial A^{DH}}{\partial\rho_{i}\left(x_{f}\right)}\right)_{T,V,N_{j\neq i}}=-\frac{e^{2}z_{i}^{2}\kappa}{24\pi\varepsilon_{0}\varepsilon_{r}}\left[2\chi_{i}+\frac{\sum_{k}x_{k,f}z_{k}^{2}\sigma_{k}}{\sum_{k}x_{k,f}z_{k}^{2}}\right]
```

##### `lngamma_i_infinite_dilution`
- Label: `eq:lngamma_i_infinite_dilution`
- Source: \cite{Bulow2021}, Eq.~(13)
- Status: Manual literature match
- Description: Defines the Debye screening quantity used in debye and huckel electrolyte term contribution.
- Change note: Mapped manually to the infinite-dilution ionic activity-coefficient relation.
- LaTeX: `docs/latex/equations.tex:2690`
- C++: No `EqID` owner comment has been attached yet.

```tex
\ln\gamma_{i}^{*}\left(x_{f,i}\right)=-\frac{e^{2}z_{i}^{2}\kappa}{24\pi\varepsilon_{0}\varepsilon_{r}}\left[2\chi_{i}+\frac{\sum_{k}x_{k,f}z_{k}^{2}\sigma_{k}}{\sum_{k}x_{k,f}z_{k}^{2}}\right]
```

##### `gamma_pm_x`
- Label: `eq:gamma_pm_x`
- Source: \cite{Bulow2021}, Eq.~(26)
- Status: Adapted implementation form
- Description: Gives an activity or fugacity-coefficient relation in debye and huckel electrolyte term contribution.
- Change note: Lower similarity; likely algebraically adapted for implementation or combined terms.
- LaTeX: `docs/latex/equations.tex:2701`
- C++: No `EqID` owner comment has been attached yet.

```tex
\gamma_{\pm}^{*}=(\gamma_{c}^{*,\nu c}\cdot\gamma_{a}^{*,\nu a})^{\frac{1}{\nu_{c}+\nu_{a}}}
```
