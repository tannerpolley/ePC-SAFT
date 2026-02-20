## SUPPORTING INFORMATION

# Calculation of multiphase equilibria containing mixed solvents and mixed 

electrolytes: General formulation and Case studies.

Moreno Ascani, ${ }^{\text {a }}$ Gabriele Sadowski ${ }^{\text {a }}$, Christoph Held ${ }^{\text {a\# }}$

a Laboratory of Thermodynamics, Department of Biochemical and Chemical Engineering, TU Dortmund, Emil-Figge Str. 70, 44277 Dortmund, Germany
\# corresponding author: christoph.held@tu-dortmund.de

## Content: 1) ePC-SAFT framework

## 2) ePC-SAFT parameters

## 3) References

## 1) ePC-SAFT framework

The EoS ePC-SAFT advanced ${ }^{\mathrm{S} 1, \mathrm{~S} 2}$ was used in this work to model the phase behavior of the investigated systems. Like the previous electrolyte version ePC-SAFT revised ${ }^{\mathrm{S} 3}$ and the original PC-SAFT ${ }^{\mathrm{S} 4}$, the thermodynamic behavior of the system is expressed as the residual Helmholtz energy at given values of temperature T , molar volume v (or density rho) and composition $\bar{x}$ (Eq. (S1)).

$$
a^{r e s}=a^{r e s}(T, v, \bar{x})
$$

The residual Helmholtz energy is formulated as a sum of different contributions. In addition to the contribution due to original PC-SAFT, which are the hard-chain $a^{H C}$, the dispersion $a^{\text {disp }}$ and the association $a^{\text {assoc }}$ contribution, ePC-SAFT advanced includes the Debye-Hückel $a^{D H}$ and the Born term $a^{\text {Born }}$ in the final expression of the residual Helmholtz energy (Eq. (S2)).

$$
\begin{equation*}
a^{\text {res }}=a^{H C}+a^{\text {disp }}+a^{\text {assoc }}+a^{D H}+a^{\text {Born }} \tag{S2}
\end{equation*}
$$

In Eq. (S2), each term has a given physical meaning and aims at representing a specific molecular interaction. The hard-chain term $a^{H C}$ accounts for the repulsions forces due to the own volume of the molecules and for their non-spherical form. The dispersion term $a^{\text {disp }}$ accounts for short-range anisotrope attractive forces among the molecules, such as Van der Waals or weak dipolar forces. The association term $a^{\text {assoc }}$ accounts for short-range strong and highly directional forces such as hydrogen bonds. The Debye-Hückel $a^{D H}$ and the Born term $a^{\text {Born }}$ are needed to describe intermolecular interactions due to charged species in the system: the Debye-Hückel term $a^{D H}$ accounts for the longrange electrostatic forces among ions whereas the Born term $a^{\text {Born }}$ describes electrostatic ion-dipole interactions of the ions with the surrounding (charged and non-charged) components. The value of the fugacity coefficient $\varphi_{\mathrm{i}}$ of each component $i=1, \ldots, N$ can be obtained from the residual Helmholtz energy provided by ePC-SAFT advanced using Eq. (S3).

$$
\begin{equation*}
\ln \left(\varphi_{\mathrm{i}}\right)=\frac{\mu_{\mathrm{i}}^{\mathrm{res}}}{\mathrm{k}_{\mathrm{B}} \cdot \mathrm{~T}}-\ln \left(1+\left(\frac{\partial\left(\frac{\mathrm{a}^{\mathrm{res}}}{\mathrm{k}_{\mathrm{B}} \cdot \mathrm{~T}}\right)}{\partial \rho}\right)\right) \tag{S3}
\end{equation*}
$$

Like ePC-SAFT revised and the original PC-SAFT, ePC-SAFT advanced requires five pure-component parameters for each component: the segment diameter $\sigma_{i}$, the number of segments $m_{i}^{\text {seg }}$, the dispersionenergy parameter $u_{i}$, the association-energy parameter $\varepsilon^{A_{i} B_{i}}$ and association-volume parameter $\kappa^{A_{i} B_{i}}$. Furthermore, for each pair of components $i$ and $j$ a binary interaction parameter $k_{i j}$ is required. Formally, this corrects for deviations from the adopted combining rules for the dispersion energy between different components and is used to obtain quantitative agreement with experimental data. The used combining rules of Berthelot Lorentz are given by Eqs. (S4) - (S6).

| $\begin{aligned} \sigma_{i j} & =\frac{1}{2}\left(\sigma_{i}+\sigma_{j}\right)\left(1-l_{i j}\right) \\ u_{i j} & =\sqrt{u_{i} u_{j}}\left(1-k_{i j}(T)\right) \\ \varepsilon^{A_{i} B_{j}} & =\frac{\varepsilon^{A_{i} B_{i}}+\varepsilon^{A_{j} B_{j}}}{2}\left(1-k_{i j}^{h b}\right) \end{aligned}$ | (S4) <br> (S5) <br> (S6) |
| :--- | :--- |

The binary interaction parameter $k_{i j}(T)$ is described as a linear function of the temperature according to Eq (S7).

$$
\begin{equation*}
k_{i j}(T)=k_{i j, 298.15}+k_{i j \_T}(T-298.15) \tag{S7}
\end{equation*}
$$

The binary interaction parameters $l_{i j}$ and $k_{i j}^{h b}$ from Eqs. (S4) and (S6) are usually set to zero. In this work, they were used only between water and 1-butanol to improve the agreement between calculated and experimental binary LLE.

## 2) ePC-SAFT parameters

The pure-component parameters used in this work to model non-charged components are listed in Table S1.

Table S1: Pure-component parameters of the non-charged components used in this work.
| component | $m_{i}^{\text {seg }}$ | $\sigma_{i} / \AA$ | $u_{i} / k_{B} / \mathrm{K}$ | $\varepsilon^{A i B i} / \mathrm{K}$ | $\kappa^{\text {AiBi }}$ | Ref. |
| :--- | :--- | :--- | :--- | :--- | :--- | :--- |
| Water | 1.2047 | * | 353.95 | 2425.7 | 0.04509 | S5 |
| n-hexane | 3.0578 | 3.7983 | 236.77 | - | - | S4 |
| n-dodecane | 5.3060 | 3.8959 | 249.21 | - | - | S4 |
| 1-butanol | 2.7510 | 3.6139 | 259.59 | 2544.56 | 0.00669 | S6 |
| 1-propanol | 3.0000 | 3.2522 | 233.40 | 2276.78 | 0.01527 | S7 |
| $* \sigma=2.7927+\left(10.11 \cdot \mathrm{e}^{-0.01775 \mathrm{~T}}-1.417 \cdot \mathrm{e}^{-0.01146 \mathrm{~T}}\right)$ |  |  |  |  |  |  |


Pure-component parameters used to model ions are listed in Table S2.

Table S2: Pure-component parameters of the charged species used in this work.
| Ionic species | $m_{i}^{\text {seg }}$ | $\sigma_{i} / \AA$ | $u_{i} / k_{B} / \mathrm{K}$ | Ref. |
| :--- | :--- | :--- | :--- | :--- |


| $\left[\mathrm{P}_{66614}\right]^{+}$ | 1.4872 | 3.5926 | 206.49 | S8 |
| :--- | :--- | :--- | :--- | :--- |
| [DCA] ${ }^{-}$ | 3.7432 | 3.8771 | 509.31 | S9 |
| $\left[\mathrm{NH}_{4}\right]^{+}$ | 1 | 3.5740 | 230.00 | S10 |
| $[\mathrm{K}]^{+}$ | 1 | 3.3417 | 200.00 | S10 |
| $[\mathrm{Na}]^{+}$ | 1 | 2.8232 | 230.00 | S10 |
| $[\mathrm{Cl}]^{-}$ | 1 | 2.7560 | 170.00 | S10 |

The binary interaction parameters between water and the organic solvents used in this work are listed in Table S3.

Table S3: Binary interaction parameters between water and the organic solvents used in this work.
|  | 1-butanol | 1-pentanol |
| :--- | :--- | :--- |
| $k_{i j \_298.15 K}$ | -0.102 | -0.017 |
| $k_{i j_{-} T}$ | $2.94 \times 10^{-4}$ | 0 |
| $l_{i j}$ | -0.0044 | 0 |
| $k_{i j}^{h b}$ | 0.026 | 0 |
| Ref. | S6 | S11 |


Binary interaction parameters between ions were taken from Held et al. ${ }^{\mathrm{S} 3}$ and are shown in Table S 4 .

Table S4: Binary interaction parameters $k_{i j \_298.15 K}$ between the cations and the chloride anion used in this work.
| $k_{i j \_298.15 K}$ | $\mathrm{Na}^{+}$ | $\mathrm{K}^{+}$ | $\mathrm{NH}_{4}{ }^{+}$ |
| :---: | :---: | :---: | :---: |
| $\mathrm{Cl}^{-}$ | 0.317 | 0.064 | -0.566 |


## References

(S1) Bülow, M.; Ascani, M.; Held, C. ePC-SAFT advanced-Part I: Physical meaning of including a concentration-dependent dielectric constant in the born term and in the Debye-Hückel theory. Fluid Phase Equilibria 2021, 535, 112967.
(S2) Bülow, M.; Ascani, M.; Held, C. ePC-SAFT advanced-Part II: Application to Salt Solubility in Ionic and Organic Solvents and the Impact of Ion Pairing. Fluid Phase Equilibria 2021, 537, 112989.
(S3) Held, C.; Reschke, T.; Mohammad, S.; Luza, A.; Sadowski, G. ePC-SAFT revised. Chemical Engineering Research and Design 2014, 92, 2884-2897.
(S4) Gross, J.; Sadowski, G. Perturbed-chain SAFT: An equation of state based on a perturbation theory for chain molecules. Ind. Eng. Chem. Res. 2001, 40, 1244-1260.
(S5) Held, C.; Cameretti, L. F.; Sadowski, G. Modeling aqueous electrolyte solutions: Part 1. Fully dissociated electrolytes. Fluid Phase Equilibria 2008, 270, 87-96.
(S6) Nann, A.; Held, C.; Sadowski, G. Liquid-Liquid Equilibria of 1-Butanol/Water/IL Systems. Ind. Eng. Chem. Res. 2013, 52, 18472-18481.
(S7) Gross, J.; Sadowski, G. Application of the perturbed-chain SAFT equation of state to associating systems. Ind. Eng. Chem. Res. 2002, 41, 5510-5515.
(S8) Ji, X.; Held, C.; Sadowski, G. Modeling imidazolium-based ionic liquids with ePC-SAFT. Fluid Phase Equilib. 2012, 335, 64-73.
(S9) Ji, X.; Held, C. Modeling the density of ionic liquids with ePC-SAFT. Fluid Phase Equilib. 2016, 410, 9-22.
(S10) Held, C.; Reschke, T.; Mohammad, S.; Luza, A.; Sadowski, G. ePC-SAFT revised. Chemical Engineering Research and Design 2014, 92, 2884-2897.
(S11) Fuchs, D.; Fischer, J.; Tumakaka, F.; Sadowski, G. Solubility of amino acids: Influence of the pH value and the addition of alcoholic cosolvents on aqueous solubility. Ind. Eng. Chem. Res. 2006, 45, 6578-6584.

