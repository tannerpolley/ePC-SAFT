# Modeling aqueous electrolyte solutions Part 1. Fully dissociated electrolytes 

Christoph Held, Luca F. Cameretti, Gabriele Sadowski*<br>Laboratory of Thermodynamics, Faculty of Biochemical and Chemical Engineering, Dortmund University of Technology, Emil-Figge-Str. 70, 44227 Dortmund, Germany

## A R T I C L E I N F O

## Article history:

Received 12 March 2008
Received in revised form 17 June 2008
Accepted 19 June 2008
Available online 28 June 2008

## Keywords:

Thermodynamic properties
Density
Vapor-liquid equilibria
Mean ionic activity coefficient
Equation of state
PC-SAFT
Debye-Hückel
Aqueous solutions
Binary mixtures
Primitive model


#### Abstract

Liquid densities (pvT), vapor pressures (VLE), and mean ionic activity coefficients (MIAC) at $25^{\circ} \mathrm{C}$ of 115 single-salt electrolyte solutions containing univalent up to trivalent ions are modeled with the ePC-SAFT equation of state proposed by Cameretti et al. [L.F. Cameretti, G. Sadowski, J.M. Mollerup, Ind. Eng. Chem. Res. 44 (2005) 3355-3362; ibid., 8944]. For each ion, only two model parameters were adjusted to experimental density and MIAC data. Without using any additional binary parameters, ePC-SAFT is able to reproduce experimental data of the respective salt solutions up to high electrolyte molalities. Moreover, it is even able to describe the reversed MIAC series for alkali hydroxides and fluorides.


© 2008 Elsevier B.V. All rights reserved.

## 1. Introduction

Electrolyte solutions play an important role in biological and chemical engineering. Applications of these systems can be found in waste and drinking water treatment, fertilizer production, or enhanced oil recovery [2]. Other important technical processes based on thermodynamic properties of electrolytes are electrolysis [3], wet flue-gas scrubbing, osmosis and reverse osmosis of aqueous solutions, as well as reactive distillation with an electrolyte serving as entrainer. Furthermore, electrolytes are used to stabilize biomolecules or to salt them out.

However, the first step towards modeling of complex multicomponent biological solutions is the ability to accurately describe quasi-binary aqueous electrolyte systems. For that purpose, besides

[^0]the non-ionic short-range interactions resulting from repulsive as well as attractive forces, also the long-range Coulombic interactions of the charged species need to be accounted for. The first ones can be described by using either $\mathrm{g}^{\mathrm{E}}$ models or equations of state (EOS). The Coulombic interactions can be modeled, e.g. by the Debye-Hückel theory (DH) developed in 1923 [4] or by the mean spherical approximation (MSA) introduced by Waismann and Lebowitz [5] in 1970 who solved the Ornstein-Zernike equation for a fluid of charged spheres of equal size.

Examples for electrolyte $\mathrm{g}^{\mathrm{E}}$ models are the electrolyte NRTL model [6] or the Pitzer model [7]. Nasirzadeh et al. [8] used a MSANRTL model [9] as well as an extended Pitzer model of Archer [10] to describe osmotic coefficients of lithium hydroxide solutions. Both, MSA-NRTL and the Pitzer-Archer model turned out to be excellent models for the description of activity coefficients in electrolyte solutions. However, nine adjustable parameters are needed for the Pitzer-Archer model, five of them being temperature dependent. Six parameters have to be adjusted using the MSA-NRTL, one of them even being a function of concentration.

The second group of models is represented by equations of state. Myers et al. [11] developed an electrolyte model based on the Peng-Robinson EOS (PREOS). They used a Born-energy term for charging up the uncharged reference system in a continuous
medium as well as a restricted version of the MSA. The resulting EOS requires three adjustable salt parameters. Osmotic coefficients and salt activity coefficients of 138 aqueous single-salt solutions were modeled with an overall average relative deviation (ARD) of $1.95 \%$. However, the PREOS is known to fail when predicting liquid densities [1,11]. Fürst and Renon [12] proposed the combination of a modified Redlich-Kwong-Soave EOS with the MSA accounting for the ionic part. Using six correlation parameters, the osmotic coefficients of 28 halide systems and 35 non-halide systems could be calculated with a root mean square relative deviation of $2.9 \%$ and $2.2 \%$, respectively.

Since its development in 1989, the Statistical Associating Fluid Theory (SAFT) has been applied to many different systems, including electrolyte solutions. Paricaud et al. [13] gave an overview of the developments using the SAFT in order to examine the effect of added salt on the vapor-liquid equilibria (vapor pressure and density) of aqueous solutions. Liu and co-workers [14] combined SAFT with the MSA primitive model. They could describe mean ionic activity coefficients (MIAC) and solution densities (pvT) for 30 electrolytes with an overall ARD of $1.6 \%$. Galindo et al. [15] successfully extended the SAFT-VR EOS to electrolyte solutions by using an additive Coulombic (MSA as well as DH) contribution. Their model yields good results for vapor pressures and liquid densities of aqueous solutions of univalent ions. However, they only marginally considered the MIAC. Radosz et al. recently published their SAFT1 [16,17] and SAFT2 [18,19] EOS yielding excellent results concerning the properties of aqueous single-salt and multi-salt solutions. SAFT1 and SAFT2 are hybrid models that treat a salt as one molecule consisting of two segments, the cation and the anion. Consequently, three individual-ion parameters as well as one additional salt parameter have to be adjusted for each electrolyte solution. SAFT1 was used to describe six aqueous alkali halide solutions. It was modified in SAFT2 by a new dispersion term which was applied to single-salt aqueous systems. They could be modeled with absolute ARDs of $0.63 \%$ for the MIAC and $0.45 \%$ for the liquid densities, respectively.

In this work we use the ePC-SAFT model developed by Cameretti et al. [1]. It is a combination of the PC-SAFT EOS by Gross and Sadowski [20] and the Debye-Hückel contribution [4] which accounts for the Coulomb interactions. It considers the ionic species independent of the salt they are part of. Only two parameters are used to characterize each ion. The first one is the ionic diameter $\sigma_{j}$ which is actually the diameter of the hydrated ion. As the DH contribution only accounts for Coulombic forces among the ions, we describe the interactions between ion and water (hydration) by means of dispersion. This yields to a second ionic parameter: the dispersive-energy parameter $u_{j} / k_{\mathrm{B}}$, which reflects the strength of ionic hydration. The parameters in the previous work [1] were obtained by simultaneous fitting to density data and vapor pressures (VLE) of salt solutions. Calculations for 12 salts could be performed with an overall ARD of $0.7 \%$ (pvT) and $2.4 \%$ (VLE), respectively. However, using this parameter set, the MIAC can only be described in poor agreement with experimental data. Therefore, in this work we apply a new approach for parameter estimation and present a new parameter set that is able to reasonably describe the solution densities as well as the VLE and MIAC of about 115 electrolyte systems.

## 2. ePC-SAFT equation of state

The ePC-SAFT EOS is based on a perturbation theory where the hard-chain system is used as the reference system. All other contributions are considered as additive contributions that can be considered independently. Thus, the residual Helmholtz energy $a^{\text {res }}$
can be written as

$$
\begin{equation*}
\frac{A}{N}^{\text {res }}=a^{\text {res }}=a^{\text {hc }}+a^{\text {disp }}+a^{\text {assoc }}+a^{\text {ion }} \tag{1}
\end{equation*}
$$

where $N$ is the total number of molecules. $a^{\text {hc }}$ represents the hard-chain repulsion of the reference system. $a^{\text {disp }}, a^{\text {assoc }}$, and $a^{\text {ion }}$ account for the Helmholtz-energy contributions due to dispersive, associative, and Coulomb interactions, respectively. Whereas expressions for $a^{\text {disp }}$ and $a^{\text {assoc }}$ are used as in the original PCSAFT model [20], Cameretti et al. [1] used a Debye-Hückel term to describe the Helmholtz energy contribution $a^{\text {ion }}$ to a system that is caused by charging the species $j$ :

$$
\begin{equation*}
\frac{a^{\text {ion }}}{k_{\mathrm{B}} T}=-\frac{\kappa}{12 \pi k_{\mathrm{B}} T \varepsilon} \times \sum_{j} x_{j} q_{j}^{2} \chi_{j} \tag{2}
\end{equation*}
$$

Here, $x_{j}$ and $q_{j}$ are the mole fraction and the charge of ion $j$, respectively. $k_{\mathrm{B}}$ represents the Boltzmann constant and $T$ the temperature. The ions are treated as spherical species in a uniform dielectric continuum characterized by a dielectric constant $\varepsilon$. They can approach each other to the distance $a_{j}$ which is equivalent to the ion diameter $\sigma_{j}$. In contrast to some other groups [21,22], we use a concentrationindependent dielectric constant $\varepsilon$ for water. The quantity $\chi_{j}$ in Eq. (2) is defined as

$$
\begin{equation*}
\chi_{j}=\frac{3}{\left(\kappa a_{j}\right)^{3}} \times\left[\frac{3}{2}+\ln \left(1+\kappa a_{j}\right)-2\left(1+\kappa a_{j}\right)+\frac{1}{2}\left(1+\kappa a_{j}\right)^{2}\right] \tag{3}
\end{equation*}
$$

with $\kappa$ being the inverse Debye screening length given by

$$
\begin{equation*}
\kappa=\sqrt{\frac{N_{A}}{k_{\mathrm{B}} T \varepsilon} \times \sum_{j} q_{j}^{2} c_{j}}=\sqrt{\frac{\rho_{N} e^{2}}{k_{\mathrm{B}} T \varepsilon} \times \sum_{j} z_{j}^{2} x_{j}} \tag{4}
\end{equation*}
$$

$c_{j}$ is the molar concentration, $\rho_{N}$ the number density of the system and $N_{\mathrm{A}}$ is Avogadro's constant, respectively.

The repulsive interactions of the ions and the attractive interaction with water (hydration) are accounted for in $a^{\text {hc }}$ and $a^{\text {disp }}$, respectively.

Using ePC-SAFT, three pure-component parameters are used to describe the molecular properties of a molecule: the segment number $m_{\text {seg }}$, the segment diameter $\sigma$, and the dispersion energy $u / k_{\mathrm{B}}$. For associating components like water two additional parameters are required, namely the association energy $\varepsilon_{\mathrm{hb}}^{\mathrm{A} 1 \mathrm{~B} 1} / k_{\mathrm{B}}$ and the association volume $\kappa_{\mathrm{hb}}^{\mathrm{A} 1 \mathrm{~B} 1}$. The association scheme used here for water is the two-site 2B approach [23].

The segment number of ions is always set to one ( $m_{\text {seg }, j}=1$ ) yielding finally to two parameters for each ion: the diameter $\sigma_{j}$ and the dispersion energy $u_{j} / k_{\mathrm{B}}$ of the hydrated ion. Since $u_{j} / k_{\mathrm{B}}$ was determined from aqueous electrolyte solution data, it gives a direct hint to which extent the ion interacts with water.

To describe salt solutions, the conventional Berthelot-Lorenz combining rules are used:

$$
\begin{equation*}
\sigma_{i j}=\frac{1}{2}\left(\sigma_{i}+\sigma_{j}\right) \tag{5}
\end{equation*}
$$

$$
\begin{equation*}
u_{i j}=\sqrt{u_{i} u_{j}} \tag{6}
\end{equation*}
$$

Eq. (6) is applied for interactions between water and ions only. Dispersion between two ions is neglected in this work. Often, a binary interaction parameter $k_{i j}$ is introduced for the dispersive-energy parameter to account for deviations from the geometric-mean rule in the mixture. This parameter is usually fitted to the respective binary mixture's properties. However, since the parameter $u_{j} / k_{\mathrm{B}}$ had already been fitted to salt/water systems, a $k_{i j}$ is not applied here.

Table 1
PC-SAFT parameters for water used in this work
| Parameter | Unit | Abbreviation | Value |
| :--- | :--- | :--- | :--- |
| Segment number | - | $m_{\text {seg }, 1}$ | 1.204659 |
| Segment diameter | Å | $\sigma_{1}$ | 2.792700 |
|  | Å | $T_{\text {dep, } 1}$ | 10.1100 |
|  | 1/K | $T_{\text {dep }, 2}$ | -0.01775 |
|  | Å | $T_{\text {dep, } 3}$ | -1.41700 |
|  | 1/K | $T_{\text {dep }, 4}$ | -0.01146 |
| Dispersion energy | K | $u_{1} / k_{\mathrm{B}}$ | 353.9449 |
| Association sites | - | $\mathrm{N}_{1}$ | 2 |
| Association energy | K | $\varepsilon_{\mathrm{hb}}^{\mathrm{A} 1 \mathrm{~B} 1} / \mathrm{k}_{\mathrm{B}}$ | 2425.6714 |
| Association volume | - | $\kappa_{\mathrm{hb}}^{\mathrm{A} 1 \mathrm{~B} 1}$ | 0.0450989 |


## 3. Parameter estimation

The first step in calculating properties of aqueous electrolyte solutions is to model the solvent water as accurately as possible. The parameters of water described by the 2B association model [23] are resumed from Cameretti et al. [1] and are given in Table 1. To account for the density anomaly of water, this parameter set was readjusted between 0 and $100^{\circ} \mathrm{C}$ by introducing a temperaturedependent segment diameter $\sigma_{\mathrm{T}, \mathrm{W}}$ for water. The latter quantity is related to the temperature-independent one ( $\sigma_{\mathrm{w}}$ ) by

$$
\begin{equation*}
\sigma_{\mathrm{T}, \mathrm{~W}}=\sigma_{\mathrm{W}}+T_{\mathrm{dep}, 1} \times \exp \left\{T_{\mathrm{dep}, 2} T\right\}+T_{\mathrm{dep}, 3} \times \exp \left\{T_{\mathrm{dep}, 4} T\right\} \tag{7}
\end{equation*}
$$

The respective coefficients are given in Table 1. Using these parameters, the absolute average mean deviation for a temperature range between 0 and $100^{\circ} \mathrm{C}$ is $0.06 \%$ for the densities and $0.52 \%$ for the vapor pressures of water, respectively. Note that without introducing the temperature-dependence of $\sigma_{\mathrm{w}}$ in Eq. (7) it is not possible to accurately model the water density between 0 and $25^{\circ} \mathrm{C}$, respectively (compare Ref. [1]).

For the calculation of solution densities, VLE, and MIAC of electrolyte solutions, several assumptions have to be made. First, a reasonable approximation within the temperature range of this work is that the vapor phase above the solution consists of pure water only. Secondly, the considered electrolytes were regarded as

Table 2
ePC-SAFT parameters used in this work; only valid with parameter set of water in Table 1
| Univalent cations |  |  | Univalent anions |  |  |
| :--- | :--- | :--- | :--- | :--- | :--- |
| Ion | $\sigma_{j}(\AA)$ | $u_{j} / k_{\mathrm{B}}(\mathrm{K})$ | Ion | $\sigma_{j}(\AA)$ | $u_{j} / k_{\mathrm{B}}(\mathrm{K})$ |
| $\mathrm{H}^{+}$ | 2.2740 | 1616.4939 | $\mathrm{F}^{-}$ | 1.6132 | 648.3127 |
| $\mathrm{Li}^{+}$ | 1.8177 | 2697.2795 | $\mathrm{Cl}^{-}$ | 3.0575 | 47.2878 |
| $\mathrm{Na}^{+}$ | 2.4122 | 646.0504 | $\mathrm{Br}^{-}$ | 3.4573 | 60.2216 |
| $\mathrm{K}^{+}$ | 2.9698 | 271.0518 | $\mathrm{I}^{-}$ | 3.9319 | 80.4347 |
| $\mathrm{Rb}^{+}$ | 3.1443 | 215.3697 | $\mathrm{OH}^{-}$ | 1.6401 | 2444.7555 |
| $\mathrm{Cs}^{+}$ | 3.5606 | 175.9357 | $\mathrm{SCN}^{-}$ | 4.0715 | 69.6806 |
| $\mathrm{NH}_{4}{ }^{+}$ | 3.4755 | 212.3632 | $\mathrm{ClO}_{4}{ }^{-}$ | 4.0731 | 58.423 |
| Choline ${ }^{+}$ | 5.9216 | 220.4883 | $\mathrm{ClO}_{3}{ }^{-}$ | 3.8212 | 15.4978 |
|  |  |  | $\mathrm{H}_{2} \mathrm{PO}_{4}{ }^{-}$ | 3.7026 | 0.0000 |
|  |  |  | $\mathrm{BrO}_{3}{ }^{-}$ | 3.5765 | 0.0000 |
|  |  |  | $\mathrm{NO}_{3}{ }^{-}$ | 3.3805 | 0.0000 |
| Bi-/trivalent cations |  |  | Bivalent anions |  |  |
| Ion | $\sigma_{j}(\AA)$ | $u_{j} / k_{\mathrm{B}}(\mathrm{K})$ | Ion | $\sigma_{j}(\AA)$ | $u_{j} / k_{\mathrm{B}}(\mathrm{K})$ |
| $\mathrm{Mg}^{2+}$ | 2.3229 | 8145.3298 | $\mathrm{SO}_{4}{ }^{2-}$ | 2.4538 | 0.0000 |
| $\mathrm{Ca}^{2+}$ | 2.8889 | 2146.9794 | $\mathrm{HPO}_{4}{ }^{2-}$ | 4.4608 | 0.0000 |
| $\mathrm{Sr}^{2+}$ | 2.9882 | 1677.6061 |  |  |  |
| $\mathrm{Ba}^{2+}$ | 3.0982 | 1475.9880 |  |  |  |
| $\mathrm{Co}^{2+}$ | 2.4387 | 5837.6334 |  |  |  |
| $\mathrm{Cu}^{2+}$ | 2.6955 | 2545.1445 |  |  |  |
| $\mathrm{Fe}^{2+}$ | 2.4828 | 5495.6986 |  |  |  |
| $\mathrm{Cr}^{3+}$ | 3.2380 | 6905.3450 |  |  |  |


strong electrolytes that fully dissociate into their respective cations and anions. Dispersive interactions were only considered between water-water and water-ion pairs. Two ions were assumed to interact only by repulsive ( $a^{\text {hc }}$ ) and Coulomb forces ( $a^{\text {ion }}$ ).

The two ion parameters - diameter $\sigma_{j}$ and dispersion energy $u_{j} / k_{\mathrm{B}}$ - were determined from fitting them to experimental data of aqueous salt solutions. In the previous work [1] saturated vapor pressures and liquid densities were used for that purpose. However, applying the so-determined parameters leads to poor results when calculating MIAC values. Therefore, in this work in addition to density data (between 20 and $30^{\circ} \mathrm{C}$ ) also MIAC at $25^{\circ} \mathrm{C}$ were used for the parameter estimation, which are much more sensitive to the ionic parameters than VLE data.

The MIAC $\gamma_{ \pm}^{*, x}$ of an electrolyte is defined as the geometrical mean of the mole-fraction-based rational activity coefficients of the ions in solution [24]:

$$
\begin{equation*}
\gamma_{ \pm}^{*, x}=\left(\left(\gamma_{+}^{*, x}\right)^{v_{+}}\left(\gamma_{-}^{*, x}\right)^{v_{-}}\right)^{1 /\left(v_{+}+v_{-}\right)} \tag{8}
\end{equation*}
$$

Here, $v_{+}$and $v_{-}$are the stoichiometric coefficients of cation and anion in the salt [16] which add to $v$. The rational activity coefficients $\gamma_{j}^{*, x}$ of ions can be obtained by ePC-SAFT as the ion fugacity coefficient $\varphi_{j}$ at the actual concentration related to the one at infinite diluted state $\varphi_{j}^{\infty}$ :

$$
\begin{equation*}
\gamma_{j}^{*, x}=\frac{\varphi_{j}\left(T, p, x_{j}\right)}{\varphi_{j}^{\infty}\left(T, p, x_{j} \rightarrow 0\right)} \tag{9}
\end{equation*}
$$

$\gamma_{ \pm}^{*}$ can be determined directly by potentiometric methods or indirectly by isopiestic measurements and is available for many electrolytes in aqueous solutions at $25^{\circ} \mathrm{C}$. In the literature [25] MIAC are often given on a molal basis, whereas ePC-SAFT is mole-fraction based. Thus, the following conversion from mole fraction $(x)$-based to molality ( $m$ )-based values is applied:

$$
\begin{align*}
\gamma_{ \pm}^{*, m} & =\frac{\gamma_{ \pm}^{*, x}}{\left(1+0.001 \nu M_{\mathrm{W}} m_{ \pm}\right)} \\
& =\frac{1}{\left(1+0.001 \nu M_{\mathrm{W}} m_{ \pm}\right)} \times\left[\frac{\varphi_{ \pm}^{*, m}}{\varphi_{ \pm, x_{\text {solute }} \rightarrow 0}^{*, m}}\right] \tag{10}
\end{align*}
$$

$M_{\mathrm{W}}$ and $m_{ \pm}$are the molecular weight of water and the molality of the electrolyte in moles of salt per kg of water, respectively.

Table 3
Comparison of hydrated cation sizes: experimental values from X-ray and neutron diffraction measurements [26] versus ePC-SAFT parameters
| Ion | $\sigma_{\text {exp }}(\AA)$ | $\sigma_{\text {epC-SAFT }}(\AA)$ |
| :--- | :--- | :--- |
| $\mathrm{Li}^{+}$ | 1.86 | 1.81 |
| $\mathrm{Na}^{+}$ | 2.40 | 2.41 |
| $\mathrm{~K}^{+}$ | 3.10 | 2.97 |


Table 4
Number of data points (NP), maximum molality, and deviations of ePC-SAFT from experimental density, vapor pressure, and MIAC data
| Electrolyte | Density |  |  |  | Vapor pressure |  |  |  |  | Activity coefficient |  |  |  |
| :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- |
|  | N | $m_{\text {max }}$ (mol/kg) | ARD (\%) | AAD ( $\mathrm{kg} / \mathrm{m}^{3}$ ) | N | $m_{\text {max }}$ (mol/kg) | ARD (\%) | AAD (kPa) | Reference | N | $m_{\text {max }}$ (mol/kg) | ARD (\%) | AAD |
| Fluorides |  |  |  |  |  |  |  |  |  |  |  |  |  |
| LiF | 15 | 0.04 | 0.01 | 0.06 | - | - | - | - |  |  |  |  |  |
| NaF | 18 | 0.93 | 0.24 | 2.44 | - | - | - | - |  | 17 | 1 | 2.38 | 0.02 |
| KF | 18 | 8.59 | 0.43 | 4.83 | - | - | - | - |  | 17 | 5 | 3.26 | 0.02 |
| RbF | 9 | 6.38 | 0.30 | 3.89 | - | - | - | - |  | 24 | 3.5 | 6.70 | 0.05 |
| CsF | 10 | 5.39 | 0.21 | 2.60 | - | - | - | - |  | 10 | 3.5 | 6.39 | 0.05 |
| Chlorides |  |  |  |  |  |  |  |  |  |  |  |  |  |
| HCl | 15 | 6.86 | 0.16 | 1.70 | - | - | - | - |  | 17 | 5 | 9.80 | 0.14 |
| LiCl | 18 | 15.73 | 1.83 | 21.46 | 13 | 12.69 | 5.87 | 0.76 | [31,38] | 20 | 4.5 | 9.79 | 0.1 |
| NaCl | 18 | 5.82 | 0.74 | 8.31 | 6 | 6.22 | 0.77 | 0.06 | [39] | 16 | 3.2 | 3.43 | 0.02 |
| KCl | 17 | 3.93 | 0.49 | 5.45 | 8 | 4.29 | 1.13 | 0.04 | [30,40] | 20 | 4.5 | 2.38 | 0.01 |
| RbCl | 33 | 8.27 | 0.29 | 3.73 | 8 | 6.95 | 1.55 | 0.07 | [30,39] | 32 | 7.8 | 1.35 | 0.01 |
| CsCl | 18 | 3.96 | 0.22 | 2.79 | 40 | 8.59 | 3.52 | 0.43 | [30] | 18 | 5 | 1.91 | 0.01 |
| $\mathrm{NH}_{4} \mathrm{Cl}$ | 9 | 6.23 | 0.97 | 10.16 | 14 | 5.32 | 0.62 | 0.01 | [38] | 22 | 7 | 0.95 | 0.01 |
| ChCl | 7 | 0.24 | 0.48 | 4.90 | - | - | - | - |  | 23 | 6 | 15.48 | 0.08 |
| $\mathrm{MgCl}_{2}$ | 18 | 5.66 | 1.69 | 20.77 | 40 | 4.8 | 1.9 | 0.20 | [30] | 20 | 4.5 | 11.08 | 0.19 |
| $\mathrm{CaCl}_{2}$ | 19 | 7.37 | 2.79 | 35.38 | 40 | 7.89 | 12.52 | 1.02 | [30] | 19 | 5.5 | 26.11 | 0.8 |
| $\mathrm{SrCl}_{2}$ | 24 | 3.4 | 2.20 | 27.46 | 40 | 3.2 | 0.95 | 0.10 | [30,41] | 38 | 3.5 | 11.05 | 0.08 |
| $\mathrm{BaCl}_{2}$ | 18 | 1.6 | 0.99 | 11.58 | 25 | 1.39 | 0.74 | 0.12 | [30,42] | 36 | 1.79 | 5.68 | 0.03 |
| $\mathrm{FeCl}_{2}$ | - | - | - | - | - | - | - | - |  | 32 | 2 | 7.13 | 0.04 |
| $\mathrm{CuCl}_{2}$ | 10 | 1.86 | 1.98 | 22.49 | 25 | 3.8 | 2.52 | 0.04 | [43] | 20 | 2.8 | 19.56 | 0.09 |
| $\mathrm{CoCl}_{2}$ | - | - | - | - | - | - | - | - |  | 40 | 4 | 12.10 | 0.13 |
| $\mathrm{CrCl}_{3}$ | 13 | 2.48 | 1.42 | 16.24 | - | - | - | - |  | 11 | 1.2 | 13.60 | 0.05 |
| Bromides |  |  |  |  |  |  |  |  |  |  |  |  |  |
| HBr | 20 | 10.11 | 0.38 | 5.01 | - | - | - | - |  | 20 | 5 | 9.28 | 0.19 |
| LiBr | 18 | 7.68 | 0.86 | 10.83 | 45 | 15.97 | 30.27 | 0.92 | [31] | 20 | 4.5 | 3.52 | 0.04 |
| NaBr | 18 | 6.48 | 1.03 | 12.98 | 14 | 7.98 | 1.38 | 0.10 | [30,38] | 20 | 4.5 | 1.75 | 0.01 |
| KBr | 18 | 5.6 | 0.56 | 6.96 | 35 | 4.35 | 0.68 | 0.09 | [30,39,44] | 20 | 4.5 | 1.78 | 0.01 |
| RbBr | 36 | 7.39 | 0.31 | 4.38 | - | - | - | - |  | 27 | 5 | 1.89 | 0.01 |
| CsBr | 18 | 3.13 | 0.17 | 2.16 | 30 | 5.89 | 2.41 | 0.60 | [30] | 18 | 5 | 2.56 | 0.01 |
| $\mathrm{NH}_{4} \mathrm{Br}$ | 10 | 6.81 | 0.98 | 11.72 | - | - | - | - |  | 19 | 7.99 | 6.56 | 0.04 |
| ChBr | - | - | - | - | - | - | - | - |  | 24 | 7 | 17.39 | 0.07 |
| $\mathrm{MgBr}_{2}$ | 20 | 4.44 | 0.89 | 11.6 | - | - | - | - |  | 20 | 4.5 | 13.43 | 0.35 |
| $\mathrm{CaBr}_{2}$ | 22 | 5 | 2.06 | 28.62 | 40 | 4.6 | 6.25 | 0.84 | [30,45] | 22 | 3.5 | 18.31 | 0.24 |
| $\mathrm{SrBr}_{2}$ | 21 | 1.73 | 1.09 | 13.14 | 40 | 3.34 | 1.69 | 0.20 | [30] | 40 | 2.12 | 11.91 | 0.08 |
| $\mathrm{BaBr}_{2}$ | 27 | 2.24 | 0.70 | 9.03 | 40 | 3.4 | 1.33 | 0.19 | [30] | 42 | 2.32 | 6.47 | 0.04 |
| $\mathrm{CuBr}_{2}$ | - | - | - | - | - | - | - | - |  | 39 | 3.61 | 14.32 | 0.11 |
| $\mathrm{CoBr}_{2}$ | - | - | - | - | - | - | - | - |  | 44 | 5 | 8.55 | 0.18 |
| Iodides |  |  |  |  |  |  |  |  |  |  |  |  |  |
| HI | 20 | 6.4 | 0.12 | 1.43 | - | - | - | - |  | 19 | 3.5 | 2.56 | 0.03 |
| LiI | 18 | 4.98 | 0.85 | 10.77 | 21 | 10.13 | 6.06 | 0.53 | [31] | 17 | 3 | 4.49 | 0.05 |
| NaI | 18 | 4.45 | 0.52 | 6.55 | 36 | 8.4 | 2.9 | 0.37 | [30] | 19 | 4.5 | 1.37 | 0.01 |
| KI | 18 | 4.02 | 0.16 | 1.94 | 27 | 5.65 | 1.26 | 0.09 | [30,41] | 20 | 4.5 | 1.44 | 0.01 |
| RbI | 42 | 8.74 | 0.18 | 2.92 | - | - | - | - |  | 27 | 5 | 3.30 | 0.02 |
| CsI | 18 | 2.57 | 0.06 | 0.81 | 25 | 2.6 | 0.46 | 0.06 | [30] | 18 | 3 | 4.66 | 0.02 |
| $\mathrm{NH}_{4} \mathrm{I}$ | 9 | 3.72 | 0.38 | 4.74 | 20 | 13.89 | 2.12 | 0.12 | [41] | 17 | 1.1 | 2.36 | 0.02 |
| $\mathrm{MgI}_{2}$ | 18 | 2.4 | 0.76 | 9.96 | - | - | - | - |  | 20 | 4.5 | 47.18 | 3.29 |
| $\mathrm{CaI}_{2}$ | 18 | 2.27 | 1.36 | 17.73 | 25 | 2.92 | 3.66 | 0.57 | [30] | 15 | 2 | 10.15 | 0.09 |
| $\mathrm{SrI}_{2}$ | 30 | 2.4 | 1.30 | 17.64 | 40 | 4.16 | 3.21 | 0.33 | [30] | 37 | 1.9 | 11.26 | 0.08 |
| $\mathrm{BaI}_{2}$ | 39 | 3.84 | 0.88 | 14 | - | - | - | - |  | 38 | 2 | 9.46 | 0.07 |
| $\mathrm{CoI}_{2}$ | - | - | - | - | - | - | - | - |  | 44 | 5 | 10.62 | 0.76 |
| Hydroxides |  |  |  |  |  |  |  |  |  |  |  |  |  |
| LiOH | 15 | 3.79 | 0.20 | 2.17 | 115 | 4.77 | 0.95 | 0.34 | [8] | 19 | 5 | 3.01 | 0.02 |
| NaOH | 15 | 10.88 | 0.44 | 5.07 | 18 | 12.19 | 3.54 | 0.58 | [46] | 20 | 5.5 | 6.78 | 0.05 |
| KOH | 15 | 15.51 | 0.7 | 9.44 | - | - | - | - |  | 14 | 3 | 3.03 | 0.02 |
| RbOH | 26 | 3.25 | 0.36 | 4.53 | - | - | - | - |  | - | - | - | - |
| CsOH | 12 | 1.18 | 0.21 | 2.30 | - | - | - | - |  | 12 | 1.2 | 2.72 | 0.02 |
| $\mathrm{BaOH}_{2}$ | 6 | 0.31 | 0.38 | 3.98 | - | - | - | - |  | 6 | 0.2 | 2.26 | 0.01 |
| Nitrates |  |  |  |  |  |  |  |  |  |  |  |  |  |
| $\mathrm{HNO}_{3}$ | 16 | 8.55 | 0.50 | 5.75 | - | - | - | - |  | 17 | 5 | 6.87 | 0.06 |
| $\mathrm{LiNO}_{3}$ | 20 | 3.63 | 0.75 | 8.03 | 16 | 12.86 | 1.56 | 0.03 | [40] | 22 | 5.5 | 14.25 | 0.14 |
| $\mathrm{NaNO}_{3}$ | 20 | 2.94 | 0.33 | 3.54 | 20 | 13.35 | 2.71 | 0.10 | [38,39,45] | 24 | 7 | 1.99 | 0.01 |
| $\mathrm{KNO}_{3}$ | 20 | 2.47 | 0.19 | 1.99 | 18 | 8.31 | 2.3 | 0.14 | [41] | 18 | 3.5 | 2.28 | 0.01 |
| $\mathrm{NH}_{4} \mathrm{NO}_{3}$ | 22 | 12.49 | 1.85 | 21.74 | - | - | - | - |  | 22 | 9 | 29.15 | 0.1 |
| $\mathrm{RbNO}_{3}$ | 19 | 3.65 | 0.20 | 2.53 | - | - | - | - |  | 26 | 4.5 | 2.92 | 0.01 |
| $\mathrm{CsNO}_{3}$ | 15 | 1.28 | 0.03 | 0.27 | - | - | - | - |  | 14 | 1.5 | 3.68 | 0.02 |
| $\mathrm{Mg}\left(\mathrm{NO}_{3}\right)_{2}$ | 18 | 2.25 | 0.62 | 7.04 | 14 | 5.13 | 2.36 | 0.02 | [41,42] | 20 | 4.5 | 13.14 | 0.09 |
| $\mathrm{Ca}\left(\mathrm{NO}_{3}\right)_{2}$ | 22 | 6.09 | 1.38 | 17.72 | 12 | 10.42 | 2.60 | 0.05 | [42,45] | 20 | 4.5 | 13.56 | 0.05 |
| $\mathrm{Sr}\left(\mathrm{NO}_{3}\right)_{2}$ | 21 | 1.58 | 0.91 | 10.97 | - | - | - | - |  | 19 | 4 | 11.95 | 0.03 |


Table 4 (Continued)
| Electrolyte | Density |  |  |  | Vapor pressure |  |  |  |  | Activity coefficient |  |  |  |
| :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- |
|  | N | $m_{\text {max }}$ (mol/kg) | ARD (\%) | AAD ( $\mathrm{kg} / \mathrm{m}^{3}$ ) | $N$ | $m_{\text {max }}$ (mol/kg) | ARD (\%) | AAD (kPa) | Reference | N | $m_{\text {max }}$ (mol/kg) | ARD (\%) | AAD |
| $\mathrm{Ba}\left(\mathrm{NO}_{3}\right)_{2}$ | 5 | 0.2 | 0.23 | 2.43 | - | - | - | - |  | 10 | 0.4 | 3.55 | 0.01 |
| $\mathrm{Fe}\left(\mathrm{NO}_{3}\right)_{2}$ | 6 | 3 | 5.12 | 66.74 | - | - | - | - |  | - | - | - | - |
| $\mathrm{Cu}\left(\mathrm{NO}_{3}\right)_{2}$ | 15 | 4.19 | 1.32 | 16.39 | - | - | - | - |  | 44 | 5 | 32.96 | 0.32 |
| $\mathrm{Co}\left(\mathrm{NO}_{3}\right)_{2}$ | 12 | 2.34 | 0.93 | 10.87 | - | - | - | - |  | 44 | 5 | 15.92 | 0.19 |
| $\mathrm{Cr}\left(\mathrm{NO}_{3}\right)_{3}$ | 11 | 1.8 | 1.52 | 17.5 | - | - | - | - |  | 12 | 1.4 | 15.48 | 0.05 |
| Thiocyanates |  |  |  |  |  |  |  |  |  |  |  |  |  |
| NaSCN | 10 | 9.63 | 0.10 | 1.11 | - | - | - | - |  | 14 | 5 | 2.76 | 0.02 |
| KSCN | 12 | 3.3 | 0.21 | 2.34 | - | - | - | - |  | 14 | 5 | 1.07 | 0.01 |
| $\mathrm{NH}_{4} \mathrm{SCN}$ | - | - | - | - | - | - | - | - |  | 12 | 2 | 2.48 | 0.02 |
| Chlorates |  |  |  |  |  |  |  |  |  |  |  |  |  |
| $\mathrm{HClO}_{3}$ | 6 | 3.74 | 0.03 | 0.30 | - | - | - | - |  | - | - | - | - |
| $\mathrm{LiClO}_{3}$ | 18 | 33.19 | 0.83 | 11.85 | - | - | - | - |  | 14 | 10.27 | 11.60 | 0.17 |
| $\mathrm{NaClO}_{3}$ | 7 | 4.03 | 0.31 | 3.51 | - | - | - | - |  | 23 | 3 | 2.14 | 0.01 |
| $\mathrm{KClO}_{3}$ | 7 | 0.43 | 0.06 | 0.59 | - | - | - | - |  | 7 | 0.7 | 0.96 | 0.01 |
| $\mathrm{RbClO}_{3}$ | - | - | - | - | - | - | - | - |  | 6 | 0.3 | 0.98 | 0.01 |
| $\mathrm{CsClO}_{3}$ | - | - | - | - | - | - | - | - |  | 6 | 0.3 | 0.51 | 0 |
| $\mathrm{Mg}\left(\mathrm{ClO}_{3}\right)_{2}$ | 15 | 2.24 | 0.59 | 6.80 | - | - | - | - |  | - | - | - | - |
| $\mathrm{Ba}\left(\mathrm{ClO}_{3}\right)_{2}$ | 12 | 1.04 | 0.47 | 5.45 | - | - | - | - |  | - | - | - | - |
| Perchlorates |  |  |  |  |  |  |  |  |  |  |  |  |  |
| $\mathrm{HClO}_{4}$ | 30 | 23.23 | 1.05 | 13.99 | 16 | 11.91 | 11.64 | 0.12 | [47] | 27 | 5 | 9.02 | 0.12 |
| $\mathrm{LiClO}_{4}$ | 9 | 5.52 | 1.88 | 22.43 | - | - | - | - |  | 26 | 4.5 | 4.13 | 0.06 |
| $\mathrm{NaClO}_{4}$ | 26 | 9.98 | 0.50 | 6.63 | - | - | - | - |  | 29 | 6 | 14.14 | 0.09 |
| $\mathrm{KClO}_{4}$ | - | - | - | - | - | - | - | - |  | 6 | 0.3 | 3.71 | 0.03 |
| $\mathrm{RbClO}_{4}$ | - | - | - | - | - | - | - | - |  | 6 | 0.3 | 5.10 | 0.03 |
| $\mathrm{CsClO}_{4}$ | - | - | - | - | - | - | - | - |  | 6 | 0.3 | 5.62 | 0.04 |
| $\mathrm{NH}_{4} \mathrm{ClO}_{4}$ | 9 | 1.50 | 0.19 | 2.00 | - | - | - | - |  | 22 | 2.5 | 7.97 | 0.04 |
| $\mathrm{Mg}\left(\mathrm{ClO}_{4}\right)_{2}$ | 22 | 4.39 | 0.90 | 12.09 | - | - | - | - |  | 16 | 1 | 9.53 | 0.07 |
| $\mathrm{Ca}\left(\mathrm{ClO}_{4}\right)_{2}$ | 12 | 4.53 | 0.64 | 8.07 | - | - | - | - |  | 31 | 3 | 17.61 | 0.24 |
| $\mathrm{Sr}\left(\mathrm{ClO}_{4}\right)_{2}$ | 37 | 5.24 | 0.26 | 3.12 | - | - | - | - |  | 39 | 4 | 15.91 | 0.28 |
| $\mathrm{Ba}\left(\mathrm{ClO}_{4}\right)_{2}$ | 13 | 4.46 | 0.42 | 5.81 | - | - | - | - |  | 22 | 5.5 | 6.60 | 0.06 |
| $\mathrm{Cu}\left(\mathrm{ClO}_{4}\right)_{2}$ | - | - | - | - | - | - | - | - |  | 36 | 3 | 25.66 | 0.63 |
| $\mathrm{Fe}\left(\mathrm{ClO}_{4}\right)_{2}$ | 10 | 0.15 | 0.04 | 0.38 | - | - | - | - |  | - | - | - | - |
| $\mathrm{Co}\left(\mathrm{ClO}_{4}\right)_{2}$ | - | - | - | - | - | - | - | - |  | 40 | 4 | 12.65 | 1.68 |
| Bromates |  |  |  |  |  |  |  |  |  |  |  |  |  |
| $\mathrm{LiBrO}_{3}$ | 12 | 0.82 | 0.19 | 2.03 | - | - | - | - |  | - | - | - | - |
| $\mathrm{NaBrO}_{3}$ | 6 | 2.21 | 0.08 | 0.93 | - | - | - | - |  | 23 | 2.62 | 3.01 | 0.02 |
| $\mathrm{KBrO}_{3}$ | 9 | 0.32 | 0.05 | 0.45 | - | - | - | - |  | 11 | 0.5 | 1.32 | 0.01 |
| $\mathrm{RbBrO}_{3}$ | - | - | - | - | - | - | - | - |  | 6 | 0.3 | 1.95 | 0.01 |
| $\mathrm{CsBrO}_{3}$ | - | - | - | - | - | - | - | - |  | 6 | 0.3 | 1.72 | 0.01 |
| $\mathrm{Mg}\left(\mathrm{BrO}_{3}\right)_{2}$ | 15 | 1.53 | 0.15 | 1.74 | - | - | - | - |  | - | - | - | - |
| Sulfates |  |  |  |  |  |  |  |  |  |  |  |  |  |
| $\mathrm{Li}_{2} \mathrm{SO}_{4}$ | 6 | 3.03 | 2.23 | 26.26 | 20 | 2.99 | 2.70 | 0.08 | [41] | 20 | 3.17 | 20.29 | 0.06 |
| $\mathrm{Na}_{2} \mathrm{SO}_{4}$ | 6 | 2.35 | 0.41 | 4.86 | 6 | 3.18 | 2.55 | 0.23 | [45] | 23 | 4.25 | 32.63 | 0.05 |
| $\mathrm{K}_{2} \mathrm{SO}_{4}$ | 5 | 0.64 | 0.08 | 0.86 | 18 | 0.97 | 0.39 | 0.03 | [48] | 13 | 0.69 | 1.61 | 0.01 |
| $\mathrm{Rb}_{2} \mathrm{SO}_{4}$ | 21 | 1.61 | 0.68 | 8.35 | - | - | - | - |  | 14 | 1.8 | 11.94 | 0.03 |
| $\mathrm{Cs}_{2} \mathrm{SO}_{4}$ | 16 | 1.18 | 0.40 | 5.21 | - | - | - | - |  | 16 | 1.63 | 12.31 | 0.03 |
| $\left(\mathrm{NH}_{4}\right)_{2} \mathrm{SO}_{4}$ | 19 | 7.57 | 2.16 | 26.59 | 12 | 6.01 | 2.49 | 0.07 | [44] | 19 | 4 | 17.44 | 0.03 |
| $\mathrm{MgSO}_{4}$ | 18 | 2.77 | 0.54 | 6.55 | - | - | - | - |  | 19 | 3.6 | 19.36 | 0.01 |
| $\mathrm{FeSO}_{4}$ | 10 | 1.65 | 0.47 | 5.45 | - | - | - | - |  | - | - | - | - |
| $\mathrm{CuSO}_{4}$ | 9 | 1.38 | 0.66 | 7.48 | - | - | - | - |  | 13 | 1.75 | 12.76 | 0.01 |
| $\mathrm{CoSO}_{4}$ | 5 | 0.56 | 0.57 | 6.10 | - | - | - | - |  | 15 | 1.51 | 17.58 | 0.03 |
| $\mathrm{Cr}_{2}\left(\mathrm{SO}_{4}\right)_{3}$ | 13 | 1.70 | 1.26 | 16.91 | - | - | - | - |  | 11 | 1.2 | 35.18 | 0.01 |
| Phosphates |  |  |  |  |  |  |  |  |  |  |  |  |  |
| $\mathrm{Na}_{2} \mathrm{HPO}_{4}$ | 8 | 0.79 | 1.81 | 19.20 | - | - | - | - |  | 10 | 1 | 11.72 | 0.03 |
| $\mathrm{K}_{2} \mathrm{HPO}_{4}$ | - | - | - | - | - | - | - | - |  | 27 | 0.87 | 5.86 | 0.02 |
| $\left(\mathrm{NH}_{4}\right)_{2} \mathrm{HPO}_{4}$ | - | - | - | - | - | - | - | - |  | 37 | 3.11 | 13.91 | 0.03 |
| $\mathrm{NaH}_{2} \mathrm{PO}_{4}$ | 7 | 0.58 | 0.59 | 6.00 | - | - | - | - |  | 27 | 5 | 3.02 | 0.01 |
| $\mathrm{KH}_{2} \mathrm{PO}_{4}$ | 12 | 2.45 | 0.31 | 3.41 | - | - | - | - |  | 20 | 1.8 | 2.47 | 0.01 |
| $N_{\text {salts }}$ Average |  |  | 97 | 97 |  |  | 37 | 37 |  |  |  | 106 | 106 |
|  |  |  | 0.75 | 9.14 |  |  | 3.29 | 0.24 |  |  |  | 9.17 | 0.12 |


Experimental activity coefficients as well as solution densities are taken from Lobo et al. [24]. If not available there, the latter quantity is found in Ref. [25]. Vapor-pressure data is taken from Refs. [8,31,33-44].

As mentioned before, ion-specific instead of salt-specific parameters are used in ePC-SAFT. Thus, the ionic parameters determined, e.g. for $\mathrm{Li}^{+}$are applicable to all electrolytes containing this ion. Hence, $\sigma_{\mathrm{Li}+}$ and $u_{\mathrm{Li}+} / k_{\mathrm{B}}$ have the same values in $\mathrm{LiCl}, \mathrm{LiBr}, \mathrm{LiOH}$, etc.

Obtaining such a universal set of parameters requires a simultaneous regression of several electrolyte solutions. For that purpose, 14 electrolytes AnCat were selected with $\mathrm{Cat}^{+}=\left\{\mathrm{Na}^{+}, \mathrm{Li}^{+}, \mathrm{K}^{+}\right\}$and $\mathrm{An}^{-}=\left\{\mathrm{F}^{-}, \mathrm{Cl}^{-}, \mathrm{Br}^{-}, \mathrm{I}^{-}, \mathrm{OH}^{-}\right\}$. The ion parameters were adjusted to
the respective electrolyte solutions except to LiF where no values for activity coefficients were available. After having adjusted the segment diameter and dispersion energy of these eight ions, the parameters of other ions (e.g. $\mathrm{Mg}^{2+}$ ) were determined from fitting to the respective solution data (e.g. $\mathrm{MgCl}_{2}, \mathrm{MgBr}_{2}$, and $\mathrm{MgI}_{2}$ ). All parameters obtained this way are listed in Table 2. (Note that due to the parameter estimation using aqueous-solution data, these ion parameters are linked to the water parameters as listed in Table 1.)

For the cations $\mathrm{Li}^{+}, \mathrm{Na}^{+}$, and $\mathrm{K}^{+}$Collins et al. [26] reported X-ray and neutron diffraction measurements of aqueous salt solutions providing the distance between a central cation and the nearest water oxygen. This directly compares to the ePC-SAFT parameter $\sigma_{j}$, since this in fact is the size of the hydrated cation. Comparing the experimental values to the $\sigma_{j}$ parameters given in Table 2 shows an excellent agreement of experimental data and fitted parameters. This can be seen in Table 3.

It becomes apparent from Table 2 that anions containing three or more oxygen atoms (e.g. the nitrate or bromate anion) can be modeled with $u_{j} / k_{\mathrm{B}}=0$, i.e. without accounting for dispersion interactions. This can be explained by their structure: they are sufficiently 'padded' with oxygen atoms to prevent an extensive hydration (Ref. [24], p. 122). Robinson based this suggestion on conductivity measurements from which he concluded a high mobility of these ions compared to other anions and thus only weak hydration. Although bivalent anions are in general strongly hydrated, Robinson observed the described effect also for the anion $\mathrm{SO}_{4}{ }^{2-}$ which agrees also to the zero $\mathrm{SO}_{4}{ }^{2-}$ dispersion parameter in Table 2. On the other hand sulfate salts (e.g. ammonium sulfate) have a salting-out effect on amino acids which militates in favor of strong water-sulfate interactions. Additionally, Cannon et al. [27] investigated the structure of water and stated that a huge amount (12-13) of water molecules can be found in the first hydration shell of sulfate anions. These opposing conclusions drawn based on interpreting ePC-SAFT parameters as well as on analyzing different experimental data obviously need further investigations. Furthermore, the appearance of ion paring in sulfate systems makes an analysis even more difficult.

## 4. Modeling results

Using the parameters summarized in Table 2 liquid densities, vapor pressures (not included in the parameter fitting), and solute activity coefficients (MIAC) were modeled. Table 4 summarizes the absolute average deviations (AAD) and absolute relative deviations (ARD) of solution densities (pvT), vapor pressures (VLE) and MIAC for 115 aqueous electrolyte solutions. AAD and ARD are calculated by

$$
\begin{align*}
& \mathrm{AAD}=\frac{1}{\mathrm{NP}} \sum_{k=1}^{\mathrm{NP}}\left|\left(y_{k}^{\mathrm{calc}}-y_{k}^{\mathrm{exp}}\right)\right|  \tag{11}\\
& \mathrm{ARD}=100 \times \frac{1}{\mathrm{NP}} \sum_{k=1}^{\mathrm{NP}}\left|\left(1-\frac{y_{k}^{\mathrm{calc}}}{y_{k}^{\mathrm{exp}}}\right)\right|
\end{align*}
$$

The above-mentioned thermodynamic properties of 115 aqueous electrolyte systems can be modeled reasonably by ePC-SAFT with overall ARDs of $0.75 \%$ (pvT), $3.29 \%$ (VLE), and $9.17 \%$ (MIAC). Considerably high deviations of MIAC values were found for some sulfates, nitrates, and fluorides. One reason is that ion pairing, which is expected and experimentally proven for these systems, is not accounted for at the moment. This will be subject of the second part of this work [28].

![](https://cdn.mathpix.com/cropped/cce0ca6f-d0d1-49b6-828a-b0592f95584d-06.jpg?height=593&width=867&top_left_y=224&top_left_x=1053)
Fig. 1. Liquid densities of aqueous solutions of six cesium salts related to the density of pure water at $20^{\circ} \mathrm{C}$ as function of salt molality. The lines represent calculations by ePC-SAFT. The circles represent experimental data [29]. Parameters and deviations are shown in Tables 2 and 3.

### 4.1. Solution densities and vapor pressures

As a typical example, the liquid densities of six cesium-salt solutions are shown in Fig. 1. All solution densities are presented as density differences $\Delta \rho$ between the density of a salt solution and pure water. As to be seen, the experimental data are described with high accuracy even at high salt concentrations of up to $6 \mathrm{~mol} / \mathrm{kg}$.

Although the cesium parameters have been adjusted to cesiumsalt solution densities at $20^{\circ} \mathrm{C}$, ePC-SAFT correctly predicts the respective densities also at other temperatures. This is illustrated in Fig. 2 which shows predicted liquid densities of cesium-salts solutions at $40^{\circ} \mathrm{C}$. No experimental density data were available for CsF.

The ePC-SAFT parameters can be used to predict vapor pressures of salt solutions which is exemplarily shown in Fig. 3 for some cesium salts at concentrations of up to $10 \mathrm{~mol} / \mathrm{kg}$ at different temperatures. The predicted vapor pressures agree very well with the experimental data at different temperatures. In the case of CsCl solutions, applying the universal parameters leads to slight model inaccuracies for higher system temperatures at concentrations above $4 \mathrm{~mol} / \mathrm{kg}$. At these conditions, model deviations might be due to the fact that CsCl cannot be assumed as fully dissociated any more. Therefore, the strength of $\mathrm{Cs}^{+}$hydration might be over-

![](https://cdn.mathpix.com/cropped/cce0ca6f-d0d1-49b6-828a-b0592f95584d-06.jpg?height=611&width=854&top_left_y=1899&top_left_x=1055)
Fig. 2. Liquid densities of aqueous solutions of six cesium salts at $40^{\circ} \mathrm{C}$ as function of salt molality. No data available for CsF. Same notation as in Fig. 1.

![](https://cdn.mathpix.com/cropped/cce0ca6f-d0d1-49b6-828a-b0592f95584d-07.jpg?height=332&width=813&top_left_y=226&top_left_x=176)
Fig. 3. Vapor pressures of aqueous solutions of three cesium salts at $60^{\circ} \mathrm{C}$ (squares and lines), $50^{\circ} \mathrm{C}$ (crosses and dashed lines), and $30^{\circ} \mathrm{C}$ (circles and dashed dotted lines) as function of salt molality. Experimental data from Ref. [30].

estimated at these concentrations and temperatures by ePC-SAFT leading to overestimated vapor pressures.

Fig. 4a and b show the influence of halide anions (Fig. 4a) and alkali cations (Fig. 4b) on the water activity coefficients in aqueous solutions of alkali halides. As it can be seen, the predictions with ePC-SAFT are again in good agreement with the experimental data.

For all salt solutions, at very low salt concentrations first a slight increase of the experimental water activity coefficients is observed. After passing a maximum, the values decrease continuously. This behavior holds true for every strong electrolyte in aqueous solution. Moreover, the experimental data show decreasing water activity coefficients for increasing sizes of the anion but for decreasing sizes of the cations. This means that the ion-water interactions increase in the order $\mathrm{Cl}^{-}<\mathrm{Br}^{-}<\mathrm{I}^{-}$for the anions, but in the order $\mathrm{K}^{+}<\mathrm{Na}^{+}<\mathrm{Li}^{+}$for the cations. Since the smallest cation but the largest anion causes the lowest water activity coefficients, there are obviously two opposing effects in aqueous electrolyte solutions. The first effect is a sterical effect that allows more water molecules to be placed around the larger ions (Ref. [32], p. 209 and Ref. [33], pp. 32-33 and 55-57). This seems to be the important one when comparing the influence of different anions. The second one is an electrostatic effect: at the same charge, the smallest ion causes the strongest electrical field which results in the strongest interaction with the surrounding water molecules (Ref. [32], p. 209 and Ref. [33], pp. 32-33 and 55-57). The effect obviously dominates the alkali cation hydration. It can be further observed from Fig. 4 that the influence of the halide anions (Fig. 4a) on the water activity coefficient is much weaker than that of the alkali cations (Fig. 4b). This is due to the fact that the relative size difference in Pauling diameters [34] (Table 5) is much bigger for the alkali cations than it is for the halide anions.

![](https://cdn.mathpix.com/cropped/cce0ca6f-d0d1-49b6-828a-b0592f95584d-07.jpg?height=470&width=888&top_left_y=1921&top_left_x=140)
Fig. 4. Influence of salts on the activity coefficient of water at $30^{\circ} \mathrm{C}$. The symbols represent experimental data [25,30,31] from isopiestic or vapor pressure measurements. The lines are predictions with ePC-SAFT. (a) Influence of anions on the activity coefficient of water. Circles: LiI, squares: NaI, triangles: KI. (b) Influence of cations on the activity coefficient of water. Circles: NaI, squares: NaBr, triangles: NaCl. Largest anion and smallest cation cause the lowest water activity coefficient.

Table 5
Comparison of relative bare Pauling ionic diameters
| Cation | $\sigma(\AA)$ | Anion | $\sigma(\AA)$ |
| :--- | :--- | :--- | :--- |
| $\mathrm{Li}^{+}$ | $100 \%$ | $\mathrm{Cl}^{-}$ | $100 \%$ |
| $\mathrm{Na}^{+}$ | $160 \%$ | $\mathrm{Br}^{-}$ | $110 \%$ |
| $\mathrm{~K}^{+}$ | $220 \%$ | $\mathrm{I}^{-}$ | $120 \%$ |


The smallest cation as well as the smallest anion correspond to the $100 \%$ value.

As it can be seen from Table 2, also the ePC-SAFT diameters for the halide anions as well as for the alkali cations follow the same trend as the Pauling diameters. The size of $\mathrm{H}^{+}$is found between the lithium and sodium ion. This can be explained by the fact that $\mathrm{H}^{+}$always exists as hydronium ion in water. As already mentioned above, the size parameter $\sigma_{j}$ represent the diameter of the hydrated ions rather than the diameter of the bare ions. That's why all cation parameters are considerably larger than the Pauling radii. In contrast, the rather weakly hydrated anions possess more or less the same fitted diameter as the Pauling ones. Besides by their lower charge density, this can be explained by the structure of the water molecules: the hydrogen atoms of water can approach the anions about $0.85 \AA$ more closely than the oxygen atoms approaching the cations [26].

### 4.2. Mean ionic activity coefficients

A very important property characterizing electrolyte solutions is the mean ionic activity coefficient. As it deviates much more from unity than the activity coefficient of water, it is a much more sensitive quantity. The MIAC first decreases with increasing concentration (Fig. 5). After reaching a minimum it often increases with increasing salt concentration reaching in some cases very high values. The latter phenomenon might be an evidence for an extensive hydration [24].

Using the parameters from Table 2, ePC-SAFT performs well in modeling the MIAC for most electrolytes. Fig. 5 shows the modeled MIAC of some cesium electrolytes at $25^{\circ} \mathrm{C}$. Only CsF is calculated less accurately. At concentrations of above $0.2 \mathrm{~mol} / \mathrm{kg}$, the estimated value of $\gamma_{ \pm}^{*, m}$ for CsF is lower than that observed experimentally (Fig.5). This is observed also for rubidium fluoride, whereas sodium and potassium fluorides are modeled with considerably higher precision. Several authors report that fluoride salts form ion pairs even at room temperature [35,36]. In systems where extensive ion pairing is expected (e.g. fluorides) the model

![](https://cdn.mathpix.com/cropped/cce0ca6f-d0d1-49b6-828a-b0592f95584d-07.jpg?height=617&width=851&top_left_y=1860&top_left_x=1092)
Fig. 5. Mean ionic activity coefficients of aqueous solutions of six cesium salts at $25^{\circ} \mathrm{C}$ as function of salt molality. Experimental data from Lobo [25]. Same notation as in Fig. 1.

- based on the stoichiometric electrolyte concentration - underestimates the effects of non-ideality. The results in these systems can be improved by implementing a theory dealing with ion pair formation which will be pursued in the second part of this work [28].

Analyzing experimental mean ionic activity coefficients, e.g. of the alkali metal bromides (Fig. 6) again the hydration of ions can be compared. The following sequence is observed: $\gamma_{\text {LiBr }}^{*}>\gamma_{\text {NaBr }}^{*}> \gamma_{\mathrm{KBr}}^{*}>\gamma_{\mathrm{RbBr}}^{*}>\gamma_{\mathrm{CsBr}}^{*}$. As already concluded from Fig. 4, the smallest alkali cation (highest surface charge density) is most strongly hydrated and strong hydration(low water activity coefficient) leads to high MIAC values. The calculations for the MIAC series in Fig. 6 follow the same trend as the experimental data.

Additionally, it is generally known that $2: 1$ electrolytes (bivalent cation) like $\mathrm{CaCl}_{2}$ have much higher activity coefficients than $1: 2$ electrolytes (bivalent anion) such as, e.g. $\mathrm{Na}_{2} \mathrm{SO}_{4}$. This is again based on the extent of hydration: bivalent cations ( $\mathrm{Ca}^{2+}$ ) are much more hydrated than bivalent anions because of their higher charge density. This is also reflected by the dispersive-energy parameter $u_{j} / k_{\mathrm{B}}$ (e.g. $\mathrm{Ca}^{2+}$ has a much higher value than $\mathrm{SO}_{4}{ }^{2-}$ ).

Another interesting effect can be observed for hydroxides, fluorides, and acetates (the latter are not shown here). Electrolytes containing these anions show a reversed sequence of activity coefficients compared to alkali salts shown in Fig. 6. Here the sequence is $\gamma_{\mathrm{CSOH}}^{*}>\gamma_{\mathrm{RbOH}}^{*}>\gamma_{\mathrm{KOH}}^{*}>\gamma_{\mathrm{NaOH}}^{*}>\gamma_{\mathrm{LiOH}}^{*}$ which is illustrated for the alkali hydroxides in Fig. 7.

A possible explanation for the reversal of the activity coefficients is the so-call localized hydrolysis [37]. A hydrated cation $\mathrm{Cat}^{+}$can be regarded as a complex like $\mathrm{H}^{+}-\mathrm{OH}^{-} \ldots \mathrm{Cat}^{+} \ldots \mathrm{OH}^{-}-\mathrm{H}^{+}$with the cation placed in the center. A positive partial charge is built up at the outer surface of the hydration shell. An anion which is a strong proton acceptor (i.e. a derivative of a weak acid like $\mathrm{F}^{-}, \mathrm{OH}^{-}$, or the acetate ion), will "dock" at the hydration shell to build an acid with the respective proton. Consequently, the water molecule is split up to form the hydroxide of the cation and the acid of the anion:

$$
\mathrm{Cat}^{+} \ldots \mathrm{OH}^{-}-\mathrm{H}^{+}+\mathrm{F}^{-} \rightarrow \mathrm{Cat}^{+} \ldots \mathrm{OH}^{-} \| \mathrm{H}^{+} \ldots \mathrm{F}^{-}
$$

As a result, the effective number of the ions in the solution is decreased which also leads to a decrease of the MIAC. The localized hydrolysis depends at one side on how strong the cation polarizes the water molecules, and on the other-hand side on how strong the anion acts as proton acceptor. It is strongest in systems where

![](https://cdn.mathpix.com/cropped/cce0ca6f-d0d1-49b6-828a-b0592f95584d-08.jpg?height=582&width=647&top_left_y=1839&top_left_x=226)
Fig. 6. Mean ionic activity coefficients of five bromide salts in aqueous solution at $25^{\circ} \mathrm{C}$ as function of salt molality. Experimental data: squares ( LiBr ), stars $(\mathrm{NaBr})$, circles $(\mathrm{KBr})$, crosses $(\mathrm{RbBr})$, triangles $(\mathrm{CsBr})$. The dotted lines represent ePCSAFT calculations. Activity coefficients decrease with increasing size of the cation $\mathrm{Li}^{+}>\mathrm{Na}^{+}>\mathrm{K}^{+}>\mathrm{Rb}^{+}>\mathrm{Cs}^{+}$.

![](https://cdn.mathpix.com/cropped/cce0ca6f-d0d1-49b6-828a-b0592f95584d-08.jpg?height=593&width=650&top_left_y=228&top_left_x=1158)
Fig. 7. Mean ionic activity coefficients of aqueous solutions of four hydroxides at $25^{\circ} \mathrm{C}$ as function of salt molality. Experimental data: squares (LiOH), stars (NaOH), circles $(\mathrm{KOH})$, triangles $(\mathrm{CsOH})$. The dotted lines represent ePC-SAFT calculations. Activity coefficients increase with increasing atomic size of the cation in the sequence $\mathrm{Cs}^{+}>\mathrm{K}^{+}>\mathrm{Na}^{+}>\mathrm{Li}^{+}$.

![](https://cdn.mathpix.com/cropped/cce0ca6f-d0d1-49b6-828a-b0592f95584d-08.jpg?height=268&width=856&top_left_y=1038&top_left_x=1053)
Fig. 8. Thermodynamic properties of aqueous RbCl solutions predicted by ePC-SAFT without parameter fitting to RbCl . (a) Vapor pressures at $40^{\circ} \mathrm{C}$ (squares) and at $30^{\circ} \mathrm{C}$ (circles). Data from Ref. [30]. (b) Mean ionic activity coefficient at $25^{\circ} \mathrm{C}$. Data from Lobo et al. [25]. (c) Solution densities at $20^{\circ} \mathrm{C}$. Data from Ref. [29].

the cation has a high charge density and the anion is a strong proton acceptor. Thus, if the anion is a weak proton acceptor, as e.g. the bromide ion, the lithium salt will be found to have the highest activity coefficient compared to the other bromides (Fig. 6) because it is the most strongly hydrated cation. However, if the anion is a strong proton acceptor, e.g., the hydroxyl ion, the extent of localized hydrolysis increases in the order $\mathrm{Rb}^{+}<\mathrm{Cs}^{+}<\mathrm{K}^{+}<\mathrm{Na}^{+}<\mathrm{Li}^{+}$. This effect is large enough to reverse the order of the activity-coefficient values. It is worth mentioning, that ePC-SAFT can even describe this effect for the hydroxides and fluorides using the parameters given in Table 2. The description of the acetate systems requires accounting for ion pairing and will be considered in a subsequent publication [28].

To demonstrate the predictive power of ePC-SAFT, the parameters for the $\mathrm{Rb}^{+}$ion were adjusted to six rubidium salts except to RbCl . The parameters for the chloride ion had already been adjusted before to $\mathrm{LiCl}, \mathrm{NaCl}$, and KCl . In Fig. 8, thermodynamic properties of RbCl solutions are presented without explicitly having fitted the $\mathrm{Rb}^{+}$parameters to this salt. As it can be seen, vapor pressures, activity coefficients as well as solution densities of RbCl can be predicted precisely although these experimental data were not included in the parameter estimation.

## 5. Conclusions

In this study, the ePC-SAFT EOS proposed by Cameretti et al. [1] was applied to describe thermodynamic properties of numerous aqueous electrolyte solutions. Ion-specific parameters were used which are independent of the electrolyte the ions are part of. Two parameters, namely the diameter $\sigma_{j}$ and the dispersion energy $u_{j} / k_{\mathrm{B}}$
of the hydrated ions were adjusted to aqueous-solution densities and MIAC data. Both parameters possess a physical meaning and show reasonable trends within the ion series. Using this approach, liquid densities, vapor pressures, and solute activity coefficients of 115 aqueous electrolyte systems were modeled reasonably with overall ARDs of $0.75 \%$ (pvT), $3.29 \%$ (VLE), and $9.17 \%$ (MIAC). The predictability of the model has been proven by a quantitative description of $\Delta \rho, \Delta p$, and $\gamma_{ \pm}^{*}$ of RbCl without having adjusted the parameters of $\mathrm{Rb}^{+}$and $\mathrm{Cl}^{-}$to experimental data of that salt.

Electrolytes containing fluoride or hydroxyl anions show a reversed sequence of activity coefficients when compared to other anion series. Even this phenomenon is correctly described by ePCSAFT.

## List of symbols

| a | Helmholtz free energy per number of particles (J) |
| :--- | :--- |
| $a_{j}$ | distance at closest approach of two ions j (Å) |
| A | Helmholtz free energy (J) |
| $c_{i}$ | molarity (moles solute per liter solution) (mol/l) |
| e | elementary charge, $1.6022 \times 10^{-19}$ (C) |
| $k_{\mathrm{B}}$ | Boltzmann constant, $1.38065 \times 10^{-23}(\mathrm{~J} / \mathrm{K})$ |
| $k_{i j}$ | binary interaction parameter |
| $m$ | molality (moles solute $i$ per kg solvent) (mol/kg) |
| $m_{\text {seg }}$ | number of segments |
| M | molecular weight (g/mol) |
| N | total number of particles; number of association sites |
| NP | number of data points |
| $N_{\mathrm{A}}$ | Avogadro's constant, $6.023 \times 10^{-23}\left(\mathrm{~mol}^{-1}\right)$ |
| $p$ | pressure (kPa, bar) |
| $q_{j}$ | charge of ion $j$ (C) |
| T | temperature (K) |
| $T_{\text {dep, } 1.4}$ | parameters for the temperature-dependent segment diameter (Å, 1/K) |
| $u / k_{\mathrm{B}}$ | dispersion-energy parameter (K) |
| $x$ | mole fraction |
| $z$ | charge number |

## Greek letters

| $\gamma_{i}$ | symmetrical activity coefficient of component $i$ (related to pure component) |
| :--- | :--- |
| $\gamma_{i}^{*}$ | asymmetrical activity coefficient of component $i$ (related to infinite dilution) |
| $\varepsilon$ | dielectric constant of a medium, $\varepsilon_{\mathrm{r}} \varepsilon_{0}(\mathrm{C} / \mathrm{Vm})$ |
| $\varepsilon_{\mathrm{r}}$ | relative permittivity |
| $\varepsilon_{0}$ | permittivity in vacuum, $8.85416 \times 10^{-12}(\mathrm{C} / \mathrm{Vm})$ |
| $\varepsilon_{\mathrm{hb}}^{\mathrm{AiBi}} / k_{\mathrm{B}}$ | association-energy parameter (K) |
| $\kappa$ | Debye length, defined in Eq. (4) (1/Å) |
| $\kappa_{\mathrm{hb}}^{\mathrm{AiBi}} / k_{\mathrm{B}}$ | association-volume parameter |
| $v$ | stoichiometric factor |
| $\rho$ | density ( $\mathrm{kg} / \mathrm{m}^{3}$ ) |
| $\rho_{\mathrm{N}}$ | number density (number of particles per volume) ( $1 / \AA^{3}$ ) |
| $\sigma_{i}$ | temperature-independent segment diameter of molecule $i(\AA)$ |
| $\sigma_{\mathrm{T}, \mathrm{W}}$ | temperature-dependent segment diameter of water (Å) |
| $\varphi_{i}$ | fugacity coefficient of component $i$ |
| $\chi$ | abbreviation for expression in Eq. (3) |

## Subscripts

An, - anion
Cat, + cation
$i, j, k$ component indexes
seg segment

| W | water |
| :--- | :--- |
| $\pm$ | mean ionic |
| 0 | pure substance |
| 1 | water |

Superscripts

| assoc | association |
| :--- | :--- |
| calc | calculated |
| disp | dispersion |
| exp | experimental |
| hc | hard chain |
| hs | hard sphere |
| ion | ionic |
| $m$ | based on molality |
| res | residual |
| symm | symmetrical (related to pure component) |
| $x$ | based on mole fraction |
| ,+- | positive or negative charge |
| $\infty$ | infinitely diluted |

## Acknowledgements

The authors gratefully acknowledge the financial support by the German Society of Industrial Research (AiF) with Grant 14778N/3 and the German Science Foundation (DFG) with Grant SA 700/7. We also wish to thank our student Markus Arndt for his help in the parameter estimation.

## References

[1] L.F. Cameretti, G. Sadowski, J.M. Mollerup, Ind. Eng. Chem. Res. 44 (2005) 3355-3362, ibid., 8944.
[2] R.M. Enick, S.M. Klara, SPE Reserv. Eng. 7 (1992) 253-258.
[3] M. Luckas, J. Krissmann, Thermodynamik der Elektrolytlösungen: Eine einheitliche Darstellung der Berechnung komplexer Gleichgewichte, Springer, Berlin, 2001.
[4] P. Debye, E. Hückel, Phys. Z. 9 (1923) 185-206.
[5] E. Waisman, J.L. Lebowitz, J. Chem. Phys. 52 (1970) 4307-4311.
[6] C.C. Chen, H.I. Britt, J.F. Boston, L.B. Evans, AIChE J. 28 (1982) 588-596.
[7] K.S. Pitzer, J. Phys. Chem. 77 (1973) 268-277.
[8] K. Nasirzadeh, R. Neueder, W. Kunz, Ind. Eng. Chem. Res. 44 (2005) 3807-3814.
[9] N. Papaiconomou, J.P. Simonin, O. Bernard, W. Kunz, Phys. Chem. Chem. Phys. 4 (2002) 4435-4443.
[10] D.G. Archer, J. Phys. Chem. Ref. Data 20 (1991) 509-555.
[11] J.A. Myers, S.I. Sandler, R.H. Wood, Ind. Eng. Chem. Res. 41 (2002) 3282-3297.
[12] W. Fürst, H. Renon, AIChE J. 39 (1993) 335-343.
[13] P. Paricaud, A. Galindo, G. Jackson, Fluid Phase Equilib. 194 (2002) 87-96.
[14] W.B. Liu, Y.G. Li, J.F. Lu, Fluid Phase Equilib. 160 (1999) 595-606.
[15] A. Galindo, A. Gil-Villegas, G. Jackson, A.N. Burgess, J. Phys. Chem. B. 103 (1999) 10272-10281.
[16] S.P. Tan, H. Adidharma, M. Radosz, Ind. Eng. Chem. Res. 44 (2005) 4442-4452.
[17] X.Y. Ji, S.P. Tan, H. Adidharma, M. Radosz, Ind. Eng. Chem. Res. 44 (2005) 7584-7590.
[18] S.P. Tan, X.Y. Ji, H. Adidharma, M. Radosz, J. Phys. Chem. B. 110 (2006) 16694-16699.
[19] X.Y. Ji, H. Adidharma, Ind. Eng. Chem. Res. 45 (2006) 7719-7728.
[20] J. Gross, G. Sadowski, Ind. Eng. Chem. Res. 40 (2001) 1244-1260.
[21] R. Triolo, J.R. Grigera, L. Blum, J. Phys. Chem. 80 (1976) 1858-1861.
[22] F.X. Ball, H. Planche, W. Fürst, H. Renon, AIChE J. 31 (1985) 1233-1240.
[23] S.H. Huang, M. Radosz, Ind. Eng. Chem. Res. 29 (1990) 2284-2294.
[24] R.A. Robinson, R.H. Stokes, Electrolyte Solutions, 2nd ed., Butterworth, London, 1970.
[25] V.M.M. Lobo, J.L. Quaresma, Handbook of Electrolyte Solutions, Parts A and B, Elsevier, Amsterdam, 1989.
[26] K.D. Collins, G.W. Neilson, J.E. Enderby, Biophys. Chem. 128 (2007) 95-104.
[27] W.R. Cannon, B.M. Pettitt, J.A. Mccammon, J. Phys. Chem. 98 (1994) 6225-6230.
[28] C. Held, G. Sadowski, Fluid Phase Equilib., in preparation.
[29] J. D'Ans, H. Surawski, C. Synowietz, Landolt-Börnstein, Bd. IV/1b, Springer, Berlin, 1977.
[30] K.R. Patil, A.D. Tripathi, G. Pathak, S.S. Katti, J. Chem. Eng. Data 36 (1991) 225-230.
[31] K.R. Patil, A.D. Tripathi, G. Pathak, S.S. Katti, J. Chem. Eng. Data 35 (1990) 166-168.
[32] R. Steudel, Chemie der Nichtmetalle, 2. Auflage, W. de Gruyter, Berlin, 1998.
[33] C.H. Hamann, A. Hamnett, W. Vielstich, Electrochemistry, 2nd ed., Wiley-VCH, Weinheim, 2007.
[34] A.L. Horvath, Handbook of Aqueous Electrolyte Solutions: Physical Properties, Estimation and Correlation Methods, Ellis Horwood, Chichester, 1985.
[35] R. Buchner, G.T. Hefter, J. Barthel, J. Chem. Soc. Faraday Trans. 90 (1994) 2475-2479.
[36] A.D. Pethybridge, D.J. Spiers, J. Chem. Soc. Faraday Trans. 73 (1977) 768-775.
[37] R.A. Robinson, H.S. Harned, Chem. Rev. 28 (1941) 419-476.
[38] A. Apelblat, J. Chem. Therm. 25 (1993) 63-71.
[39] A. Apelblat, E. Korin, J. Chem. Therm. 30 (1998) 59-71.
[40] J.N. Pearce, A.F. Nelson, J. Am. Chem. Soc. 54 (1932) 3544-3555.
[41] A. Apelblat, E. Korin, J. Chem. Therm. 30 (1998) 459-471.
[42] A. Apelblat, J. Chem. Therm. 24 (1992) 619-626.
[43] T. Ewan, W.R. Ormandy, B. Berkeley, J. Chem. Soc. T. 61 (1892) 769-781.
[44] A. Apelblat, J. Chem. Therm. 25 (1993) 1513-1520.
[45] A. Apelblat, E. Korin, J. Chem. Therm. 34 (2002) 1621-1637.
[46] B. Behzadi, B.H. Patel, A. Galindo, C. Ghotbi, Fluid Phase Equilib. 236 (2005) 241-255.
[47] J.N. Pearce, A.F. Nelson, J. Am. Chem. Soc. 55 (1933) 3075-3081.
[48] H.G. Leopold, J. Johnston, J. Am. Chem. Soc. 49 (1927) 1974-1988.


[^0]:    Abbreviations: AnCat, electrolyte of the form anion-cation; AAD, absolute average deviation, defined in Eq. (11); ARD, absolute average relative deviation, defined in Eq. (11); DH, Debye-Hückel model; EOS, equation of state; $\mathrm{g}^{\mathrm{E}}$, excess Gibbs energy; MIAC, mean ionic activity coefficient; MSA, mean spherical approximation; OF , objective function; pvT, pressure-volume-temperature behavior (density data); VLE, vapor-liquid equilibrium.

    * Corresponding author. Tel.: +49 231 7552635; fax: +49 2317552572.

    E-mail address: g.sadowski@bci.tu-dortmund.de (G. Sadowski).

