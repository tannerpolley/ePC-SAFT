# Application of the perturbed chain-SAFT equation of state for modeling $\mathrm{CO}_{2}$ solubility in aqueous monoethanolamine solutions 

Sadjad Fakouri Baygi, Hassan Pahlavanzadeh*<br>Faculty of Chemical Engineering, Tarbiat Modares University, Tehran, Iran


#### Abstract

$\mathrm{CO}_{2}$ removal by treatment of acid gases by aqueous alkanolamines is a very significant operation from industrial and environmental point of view. To attain a comprehensive thermodynamic model of the $\mathrm{CO}_{2}-\mathrm{MEA}-\mathrm{H}_{2} \mathrm{O}$ in a wide range of temperature and CO2 partial pressures, Perturbed Chain-Statistical Associating Fluid Theory (PC-SAFT) EOS is applied to predict the absorption of carbon dioxide by MEA (MonoEthanolAmine). In order to find the best association scheme for MEA in PC-SAFT EOS, three pure parameter sets for MEA in the 2B, 3B and 4C association schemes are determined in temperature range $303.15-443.15 \mathrm{~K}$. Temperature independent binary interaction parameters have been adjusted in the VLE calculation for three schemes of MEA and two schemes of water. Binary VLE calculations show the 3B scheme for MEA and the 4C scheme for water indicate the best prediction in the MEA- $\mathrm{H}_{2} \mathrm{O}$ system. Excess enthalpy data for aqueous MEA are predicted by $k_{i j}$, which has been adjusted in VLE calculations. The 3B scheme for MEA and the 4 C scheme for water also are used to find $\mathrm{CO}_{2}$ solubility in the ternary system of $\mathrm{CO}_{2}-\mathrm{MEA}-\mathrm{H}_{2} \mathrm{O}$ system. Ideal Smith-Missen algorithm has been applied to find the concentration of all species in chemical equilibrium. Results show the 3B association scheme for MEA and the 4C association scheme for water in PC-SAFT EOS have better agreement with binary and ternary experimental data. PC-SAFT EOS is able to anticipate the $\mathrm{CO}_{2}$ solubility in the $\mathrm{CO}_{2}-\mathrm{MEA}-\mathrm{H}_{2} \mathrm{O}$ system without any regression in the ternary system. The $\mathrm{CO}_{2}$ solubility in ternary system is compared to e-NRTL as a common thermodynamics model. The average absolute partial pressure deviations for PC-SAFT and e-NRTL are calculated around $36 \%$ and $42 \%$, respectively.


© 2014 The Institution of Chemical Engineers. Published by Elsevier B.V. All rights reserved.

Keywords: PC-SAFT EOS; Association schemes; Monoethanolamine; Independent binary interaction; Smith-Missen Algorithm; $\mathrm{CO}_{2}$ capture

## 1. Introduction

Aqueous alkanolamines play a vital role in the removal of acid gases in many industries such as natural gas treatment and $\mathrm{CO}_{2}$ capture from flue gases. Consequently, a rigorous and predictive thermodynamics model of $\mathrm{CO}_{2}$ absorption is essential to model these types of systems. PC-SAFT EOS, developed by Gross and Sadowski $(2001,2002)$, is one of the most prevalent SAFT EOS versions, providing an applicable and reliable thermodynamics framework to model thermodynamics properties of many systems, especially asymmetric systems. The ability of PC-SAFT EOS to predict asymmetric
and associating systems was a primary consideration in its selection in this study to model $\mathrm{CO}_{2}$ solubility in aqueous MEA solutions. On the other hand, there are some practical issues about application of association term in PC-SAFT EOS, because the performance of the PC-SAFT EOS depends seriously on the selection of an association scheme due to its effect on parameter estimation. Different association schemes result in different parameter sets that yield similar pure component vapor pressure and liquid density. However, different mixture results, such as heat of mixing and binary VLE may be obtained. In addition, this issue is critical regarding alkanolamine molecules, because these molecules consist of

[^0]| Nomenclature |  |
| :--- | :--- |
| $a_{\mathrm{i}}$ | activity of species i |
| $\tilde{a}^{\text {res }}$ | residual Helmholtz energy |
| $\tilde{a}^{h c}$ | hard-chain Helmholtz energy |
| $\tilde{a}^{\text {disp }}$ | dispersion interactions Helmholtz energy |
| $\tilde{a}^{\text {assoc }}$ | association interactions Helmholtz energy |
| d | average segment diameter |
| $g_{i j}^{h s}$ | radial pair distribution function for segments of component $i$ in the hard sphere system |
| h | enthalpy |
| $k_{B}$ | Boltzmann constant ( $1.38066 \times 10^{-23} \mathrm{~J} . \mathrm{K}^{-1}$ ) |
| $\mathrm{K}_{\mathrm{X}}$ | The chemical equilibrium constants in mole fraction based |
| $K_{m}$ | The chemical equilibrium constants in molality based |
|  |  |
| m | pure number of segments in a molecule |
| $\mathrm{M}_{\mathrm{t}}$ | total moles of MEA |
| np | number of experimental data points |
| x | mole fraction in liquid phase |
| y | mole fraction in gas phase |
| $\mathrm{X}^{\mathrm{A}_{\mathrm{i}}}$ | fraction of molecules of species $i$ that are not bound at site A |
|  |  |
| $p$ | pressure (Pa) |
| T | temperature (K) |
| W | total moles of water in liquid phase |
| Greek symbols |  |
| $\alpha$ | $\mathrm{CO}_{2}$ liquid loading (mole $\mathrm{CO}_{2} /$ mole MEA) |
| $\varepsilon$ | segment dispersion interaction energy |
| $\varepsilon^{A_{i} B_{j}}$ | energy of hydrogen bonding between site $A$ at molecule $i$ and site $B$ at molecule $j$ |
| $\kappa^{\mathrm{A}_{i} \mathrm{~B}_{\mathrm{j}}}$ | volume of hydrogen bonding between site $A$ at molecule $i$ and site $B$ at molecule $j$ |
| $\rho$ | molar density |
| $\rho_{n}$ | molecular density (number of molecules/ $\AA^{3}$ ) |
| $\sigma$ | pure temperature independent segment diameter (Å) |
| $v_{i}$ | stoichiometric coefficient of component $i$ in reaction |
| $\Delta{ }^{\mathrm{A}_{\mathrm{i}} \mathrm{B}_{\mathrm{j}}}$ | strength of hydrogen bonding between the site $A$ at molecule $i$ and the site $B$ at molecule $j$ |
| Subscripts and superscripts |  |
| exp | experimental |
| cal | calculated |
| E | excess |
| L | liquid |
| min | minimum |
| max | maximum |
| i, j, k | indices |
| sat | saturation |

distinctive functional groups (hydroxyl and amine groups). Therefore, several association schemes should be checked in binary systems to find the suitable schemes (Avlund, 2011; Pahlavanzadeh and Fakouri Baygi, 2013).

In the past few decades, some $g^{E}$ thermodynamics models have been implemented for correlation and prediction of $\mathrm{CO}_{2}$ solubility in aqueous MEA solutions and other aqueous alkanolamines solutions, such as electrolyte NRTL by Zhang et al. (2011), Austgen et al. (1989) and Posey and Rochelle (1997)

![](https://cdn.mathpix.com/cropped/619bc8ff-8258-4e1f-b381-2302a83268ab-02.jpg?height=570&width=849&top_left_y=228&top_left_x=1047)
Fig. 1 - Schematic diagram of the model.

and extended UNIQUAC by Faramarzi et al. (2009). Common features shared by most of the $g^{E}$ thermodynamics models that are applied to model $\mathrm{CO}_{2}$-alkanolamine-water systems include a large number of complicated adjustable binary interaction parameters that must be correlated to multicomponent system experimental data. On the contrary, the models that use the equations of states, particularly association EOSs, do not require this amount of experimental data, except for modifying the interaction energies. However, these modifications do not require multicomponent (ternary, quaternary, etc.) system experimental data. Some association models have recently been applied in order to model these systems, such as CPA for $\mathrm{CO}_{2}$-MDEA- $\mathrm{H}_{2} \mathrm{O}$ system by Zoghi et al. (2012), SAFTVR for $\mathrm{CO}_{2}-$ MEA- $\mathrm{H}_{2} \mathrm{O}$ system and by Rodriguez et al. (2012) and Mac Dowell et al. (2009) and PC-SAFT by Nasrifar and Tafazzol (2010). Nasrifar and Tafazzol (2010) applied PC-SAFT EOS to describe solubility of $\mathrm{CO}_{2}$ and $\mathrm{H}_{2} \mathrm{~S}$ in aqueous solutions of MDEA, DEA, and MEA, but their model involves very large errors for systems involve $\mathrm{CO}_{2}$ especially $\mathrm{CO}_{2}$-MEA- $\mathrm{H}_{2} \mathrm{O}$ system.

Limitations in work of Nasrifar and Tafazzol (2010) are the motivating factors for this study. In this study, our previous work (Pahlavanzadeh and Fakouri Baygi, 2013) on thermodynamics modeling of $\mathrm{CO}_{2}$-MDEA- $\mathrm{H}_{2} \mathrm{O}$ system by PC-SAFT EOS is extended to the $\mathrm{CO}_{2}-$ MDEA- $\mathrm{H}_{2} \mathrm{O}$ system. PC-SAFT EOS was employed in this research, similarly to our earlier work, to find $\mathrm{CO}_{2}$ physical absorption, and an ideal Smith-Missen method was used to find chemical species concentration. A schematic diagram of the model is presented in Fig. 1.

## 2. PC-SAFT EOS

Gross and Sadowski $(2001,2002)$ introduced PC-SAFT EOS in terms of free Helmholtz energy.

$$
\begin{equation*}
\tilde{a}^{\text {res }}=\tilde{a}^{\text {hc }}+\tilde{a}^{\text {disp }}+\tilde{a}^{\text {assoc }} \tag{1}
\end{equation*}
$$

where $\tilde{a}^{\text {res }}, \tilde{a}^{\text {hc }}, \tilde{a}^{\text {disp }}$, and $\tilde{a}^{\text {assoc }}$ represent residual Helmholtz energy, the contribution of the hard-chain reference fluid consists of spherical segments Helmholtz energy, the contribution of dispersive attractions to the Helmholtz energy, and the contribution due to short-range association interactions (hydrogen bonding) Helmholtz energy, respectively. Gross and Sadowski (2001) presented the contribution of the hard-chain reference fluid Helmholtz energy and the contribution of dispersive attractions to the Helmholtz energy. The contribution
of the association Helmholtz energy is given by the following expression (Gross and Sadowski, 2002):

$$
\begin{equation*}
\tilde{a}^{\text {assoc }}=\sum_{i} x_{i} \sum_{A_{i}}\left[\ln X^{A_{i}}-\frac{1}{2} X^{A_{i}}+\frac{1}{2}\right] \tag{2}
\end{equation*}
$$

where $x_{i}$ is mole fraction of molecule $i$, and the $X^{A_{i}}$ is a fraction of molecules of species $i$ that are not bound at site $A$. $X^{A_{i}}$ is expressed by following equation (Gross and Sadowski, 2002):

$$
\begin{equation*}
X^{A_{i}}=\left(1+\rho_{n} \sum_{j} x_{j} \sum_{B_{j}} X^{B_{j}} \Delta^{A_{i} B_{j}}\right)^{-1} \tag{3}
\end{equation*}
$$

where $\rho_{n}$ and $\Delta^{A_{i} B_{j}}$ are molecular density and strength of hydrogen bonding between the site $A$ at molecule $i$ with the site $B$ at molecule $j$, respectively. Hydrogen bond strength is expressed by Eq. (4) (Gross and Sadowski, 2002)

$$
\begin{equation*}
\Delta^{A_{i} B_{j}}=d_{i j}^{3} g_{i j}^{h s}\left(d_{i j}\right) \kappa^{A_{i} B_{j}}\left[\exp \left(\frac{\varepsilon^{A_{i} B_{j}}}{k_{B} T}\right)-1\right] \tag{4}
\end{equation*}
$$

where $\varepsilon^{A_{i} B_{j}}$ and $\kappa^{A_{i} B_{j}}$ are energy and volume of hydrogen bonding between site $A$ at molecule $i$ and site $B$ at molecule $j$, respectively. $d_{i j}=\left(d_{i}+d_{j}\right) / 2$ is average temperature dependent segment diameter, $g_{i j}^{h s}$ is radial pair distribution function for segments of component $i$ in the hard sphere system which has been presented by Gross and Sadowski (2001), $k_{B}$ is Boltzmann constant and T represents temperature in Kelvin.

Temperature dependent segment diameter, $d_{i}$, is given by the following expression (Gross and Sadowski, 2001):

$$
\begin{equation*}
d_{i}=\sigma_{i}\left(1-0.12 \exp \left(\frac{-3 \varepsilon_{i}}{k_{\mathrm{B}} \mathrm{~T}}\right)\right) \tag{5}
\end{equation*}
$$

where $\sigma_{i}$ and $\varepsilon_{i}$ are pure temperature independent segment diameter and segment dispersion interaction energy of molecule i respectively.

Gross and Sadowski (2002) used two mixing rules for associating multicomponent systems that are expressed in Eqs. (6) and (7):

$$
\begin{equation*}
\kappa^{A_{i} B_{j}}=\sqrt{\kappa^{A_{i} B_{i}} \kappa^{A_{j} B_{j}}}\left(\frac{\sqrt{\sigma_{i} \sigma_{j}}}{\frac{1}{2}\left(\sigma_{i}+\sigma_{j}\right)}\right)^{3} \tag{6}
\end{equation*}
$$

$$
\begin{equation*}
\varepsilon^{A_{i} B_{j}}=\frac{\varepsilon^{A_{i} B_{i}}+\varepsilon^{A_{j} B_{j}}}{2} \tag{7}
\end{equation*}
$$

Specific terminology is used in the study that was introduced by Huang and Radosz (1990) to denote association structure. They have classified eight different association schemes. In their terminology, different molecules were characterized by different schemes and several possibilities for participation in association. The molecules with one proton donor site and one proton acceptor site correspond to the 2 B association scheme, the molecules with two proton donor sites and one proton acceptor site and vice versa correspond to the 3 B association scheme and the molecules with two proton donor sites and two proton acceptor sites correspond to the 4 C association scheme.

Table 1 - Gorrelations constants for the saturated vapor pressure ( $p^{\text {sat }}$ ) and saturated liquid density ( $\rho^{\text {L }}$,sat ) of MEA in Eqs. (9) and (10) (from Avlund et al. (2008)).
| Constants | $p^{\text {sat }}$ (Pa) | $\rho^{\mathrm{L}, \text { sat }}\left(\mathrm{mol} \mathrm{L}^{-1}\right)$ |
| :--- | :--- | :--- |
| A | 92.624 | 1.0011 |
| B | -10,367 | 0.22523 |
| C | -9.4699 | 678.2 |
| D | $1.9 \times 10^{-18}$ | 0.21515 |
| E | 6 | 0 |
| $\mathrm{T}_{\text {min }}(\mathrm{K})$ | 283.65 |  |
| $\mathrm{T}_{\text {max }}(\mathrm{K})$ | 678.2 |  |
| Maximum deviation | 8.1\% | 0.3\% |
| Average deviation | 2.1\% | 0.1\% |


## 3. Modeling results

### 3.1. Pure parameters calculation

Alkanolamines consist of hydroxyl groups and one or more amine groups, so sometimes selecting a precise association scheme is complicated. To address their distinctive functional groups, a variety association schemes can be proposed.

MEA is a primary alkanolamine consisting of one hydroxyl group and one primary amine group. There is some degree of ambiguity regarding the participation of primary amine group in association. Therefore, different association schemes are considered for MEA that do not specify the role of hydroxyl and primary amine groups in association; rather, these association schemes regard MEA as a molecule with some association sites that they can associate other molecules sites with same strength. For example, Nasrifar and Tafazzol (2010) considered the 4 C association scheme for MEA. On the other hand, some association schemes specify the role of hydroxyl and primary amine groups in association; for example, Mac Dowell et al. (2009) applied an association scheme that mediated three sites for hydroxyl groups and three sites for amine groups for SAFT-VR EOS.

In this study, parameter sets of MEA for PC-SAFT EOS are determined for 2B, 3B and 4C schemes. In addition, a parameter set for the 4 C association scheme are determined with less errors for saturated pressure and saturated density respect to results of Nasrifar and Tafazzol (2010).

To calculate parameter sets for MEA, the following objective function is applied

$$
\begin{equation*}
\mathrm{OF}=\sum_{i}^{n p} \frac{\left|p_{i}^{\text {sat, exp }}-p_{i}^{\text {sat, cal }}\right|}{p_{i}^{\text {sat, exp }}}+\sum_{i}^{n p} \frac{\left|\rho_{i}^{\text {sat, exp }}-\rho_{i}^{\text {sat, cal }}\right|}{\rho_{i}^{\text {sat, exp }}} \tag{8}
\end{equation*}
$$

where $p^{\text {sat }}$ and $\rho^{\text {sat }}$ indicate saturated vapor pressure and saturated liquid density of pure species, respectively. Experimental saturated vapor pressures and saturated liquid densities obtained from DIPPR correlations that Avlund et al. (2008) used for determining MEA parameters set for CPA EOS. In addition, Avlund et al. (2008) presented these correlations have reliable agreement with experimental data. The correlations are expressed by Eqs. (9) and (10) for the vapor pressure ( $p^{\text {sat }}$ ) and saturated liquid density ( $\rho^{\mathrm{L}, \text { sat }}$ ) of MEA, respectively.

$$
\begin{equation*}
p^{\text {sat }}(\mathrm{Pa})=\exp \left(\mathrm{A}+\frac{\mathrm{B}}{(\mathrm{~T} / \mathrm{K})}+\mathrm{C} \ln (\mathrm{~T} / \mathrm{K})+\mathrm{D}(\mathrm{~T} / \mathrm{K})^{E}\right) \tag{9}
\end{equation*}
$$

Table 2 - Pure component parameters used in this work.
| Species | Association scheme | $m$ | $\sigma$ | $\varepsilon / k_{\mathrm{B}}$ | $\kappa{ }^{\text {AB }}$ | $\varepsilon^{\mathrm{AB}} / \mathrm{k}_{\mathrm{B}}$ | AAD\% ${ }^{\text {a }}$ |  | Ref. |
| :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- |
|  |  |  |  |  |  |  | $p^{\text {sat }}$ | $\rho^{\text {L,sat }}$ |  |
| MEA | 4C | 4.5208 | 2.6574 | 237.6864 | 0.187533 | 989.8984 | 0.24 | 0.23 | This work |
|  | 3B | 4.5354 | 2.6019 | 204.0438 | 0.118488 | 2383.4744 | 1.75 | 0.43 | This work |
|  | 2B | 3.0353 | 3.0435 | 277.174 | 0.037470 | 2586.3 | 0.62 | 0.12 | This work |
| $\mathrm{H}_{2} \mathrm{O}$ | 2B | 1.9599 | 2.362 | 279.42 | 0.1750 | 2059.28 | 1.18 | 3.92 | Diamantonis and Economou (2011) |
|  | 4C | 2.1945 | 2.229 | 141.66 | 0.2039 | 1804.17 | 1.98 | 0.83 | Diamantonis and Economou (2011) |
| $\mathrm{CO}_{2}$ | Nonassociating | 2.0729 | 2.7852 | 169.21 | - | - | 2.78 | 2.73 | Gross and Sadowski (2001) |
| ${ }^{\text {a }}$ AAD $\%=100 / n p \times \Sigma\left\|\psi^{\text {exp }}-\psi^{\text {cal }}\right\| / \psi^{\text {exp }}, \psi$ is $p^{\text {sat }}$ and $\rho^{\text {L,sat }}$. |  |  |  |  |  |  |  |  |  |


$$
\begin{equation*}
\rho^{\mathrm{L}, \text { sat }}\left(\frac{\mathrm{mol}}{\mathrm{~L}}\right)=\frac{\mathrm{A}}{\left(\mathrm{~B}^{[1+(1-(\mathrm{T} / \mathrm{K}) / \mathrm{C})]^{\mathrm{D}}}\right)} \tag{10}
\end{equation*}
$$

where the correlations constants are presented in Table 1.
Pure parameter sets for MEA with $2 \mathrm{~B}, 3 \mathrm{~B}$ and 4 C association schemes in temperature range $303.15-443.15 \mathrm{~K}$ are determined and the results are presented in Table 2, and the diagrams of saturated vapor pressure versus saturated liquid density and saturated liquid density versus temperature of saturated pure MEA are presented in Figs. 2 and 3, respectively.

So far, a number of association schemes and parameter sets have been reported in the literature for water and carbon dioxide. Diamantonis and Economou (2011) used several saturation data for finding the accurate parameter sets for water in the 2B and 4C association schemes. So in this study, two association schemes for $\mathrm{H}_{2} \mathrm{O}$ from Diamantonis and Economou and one parameter set for $\mathrm{CO}_{2}$, are compared for selection of the best parameter sets for VLE calculations of multicomponent system. Meanwhile, $\mathrm{CO}_{2}$ molecule is not

![](https://cdn.mathpix.com/cropped/619bc8ff-8258-4e1f-b381-2302a83268ab-04.jpg?height=835&width=835&top_left_y=1758&top_left_x=164)
Fig. 2 - Saturated vapor pressure versus saturated liquid density for pure MEA at temperature from 303.15 K to $\mathbf{4 4 3 . 1 5 ~ K}$. Comparison of the 4C, 3B and 2B association schemes to DIPPR correlations.

considered as an associating molecule. These parameter sets are presented in Table 2. Binary VLE of $\mathrm{CO}_{2}-\mathrm{H}_{2} \mathrm{O}$ and MEA- $\mathrm{H}_{2} \mathrm{O}$ systems are calculated and compared by using two association schemes for $\mathrm{H}_{2} \mathrm{O}$.

![](https://cdn.mathpix.com/cropped/619bc8ff-8258-4e1f-b381-2302a83268ab-04.jpg?height=830&width=826&top_left_y=1174&top_left_x=1059)
Fig. 3 - Saturated liquid density versus temperature for pure MEA at temperature from 303.15 K to 443.15 K . Comparison of the 4C,3B and 2B association schemes to DIPPR correlations.

Table 3 - AAD\% for MEA-H2O system in different association schemes.
| Association schemes |  | AAD\% ${ }^{\text {a }}$ |  | $k_{i j}$ |
| :--- | :--- | :--- | :--- | :--- |
| MEA | $\mathrm{H}_{2} \mathrm{O}$ | $\mathrm{X}_{\text {MEA }}$ | УмеА |  |
| 2B | 2B | 5.66 | 2.38 | -0.0420 |
| 3B |  | 4.24 | 2.08 | -0.0146 |
| 4C |  | 5.41 | 2.28 | -0.0362 |
| 2B | 4C | 3.52 | 8.63 | -0.1305 |
| 3B |  | 3.54 | 1.54 | -0.0520 |
| 4C |  | 8.05 | 3.60 | -0.1245 |


[^1]![](https://cdn.mathpix.com/cropped/619bc8ff-8258-4e1f-b381-2302a83268ab-05.jpg?height=2187&width=1329&top_left_y=221&top_left_x=376)
Fig. 4 - Isobaric temperature-composition T -x slices of the vapor-liquid equilibrium of $\mathrm{MEA}-\mathrm{H}_{2} \mathrm{O}$ with three association scheme for MEA and 2B association scheme for $\mathrm{H}_{2} \mathrm{O}$. The symbols correspond to the experimental data at $\mathrm{P}=66.66 \mathrm{kPa}$ and P=101.33 kPa from Cai et al. (1996) and Park and Lee (1997), and the dotted and solid curves represent the PC-SAFT predictions $\left(k_{i j}=0\right)$ and correlations $\left(k_{i j}\right.$ from Table 3): (a) 2 B scheme for MEA at $P=66.66 \mathrm{kPa}$, (b) 2 B scheme for MEA at $P=101.33 \mathrm{kPa}$, (c) 3B scheme for MEA at $P=66.66 \mathrm{kPa}$, (d) 3B scheme for MEA at $P=101.33 \mathrm{kPa}$, (e) 4 C scheme for MEA at $P=66.66 \mathrm{kPa}$, (f) 4C scheme for MEA at $P=101.33 \mathrm{kPa}$.

![](https://cdn.mathpix.com/cropped/619bc8ff-8258-4e1f-b381-2302a83268ab-06.jpg?height=2187&width=1326&top_left_y=221&top_left_x=360)
Fig. 5 - Isobaric temperature-composition $\mathrm{T}-\mathrm{x}$ slices of the vapor-liquid equilibrium of $\mathrm{MEA}-\mathrm{H}_{2} \mathrm{O}$ with three association scheme for MEA and 4C association scheme for $\mathrm{H}_{2} \mathrm{O}$. The symbols correspond to the experimental data at $\mathrm{P}=66.66 \mathrm{kPa}$ and P = 101.33 kPa from Cai et al. (1996) and Park and Lee (1997), and the dotted and solid curves represent the PC-SAFT predictions $\left(k_{i j}=0\right)$ and correlations $\left(k_{i j}\right.$ from Table 3): (a) 2 B scheme for MEA at $P=66.66 \mathrm{kPa}$, (b) 2 B scheme for MEA at $P=101.33 \mathrm{kPa}$, (c) 3B scheme for MEA at $P=66.66 \mathrm{kPa}$, (d) 3B scheme for MEA at $P=101.33 \mathrm{kPa}$, (e) 4 C scheme for MEA at $P=66.66 \mathrm{kPa}$, (f) 4 C scheme for MEA at $P=101.33 \mathrm{kPa}$.

![](https://cdn.mathpix.com/cropped/619bc8ff-8258-4e1f-b381-2302a83268ab-07.jpg?height=773&width=844&top_left_y=223&top_left_x=173)
Fig. 6 - Excess enthalpy of MEA- $\mathrm{H}_{2} \mathrm{O}$ mixture at $\mathrm{T}=298.15 \mathrm{~K}$ and $\mathbf{T} \boldsymbol{=} \mathbf{3 4 2 . 1 5 ~ K}$. Comparison of PC-SAFT $\boldsymbol{(} \boldsymbol{k}_{\mathbf{i j}} \boldsymbol{=} \mathbf{- 0 . 0 5 2 0}$ ) calculations to experimental data from Posey (1996) and Touhara et al. (1982). $\boldsymbol{k}_{\boldsymbol{i j}}$ is adjusted to vapor-liquid equilibrium data (see Fig.5(c) and (d)). Dash lines and solid lines represent model predictions ( $k_{i j}=0$ ) and model correlations $\left(k_{i j}=-0.0520\right)$, respectively.

### 3.2. Binary systems

Finding proper parameter sets and schemes is the principal part of the modeling with association EOSs, because the pure parameter sets have to be adjusted with the saturation experimental data, and these parameter sets and association schemes exhibit same results in pure systems. Therefore, the best association schemes can be found in binary systems.

Binary interaction parameters are adjusted for binary systems to obtain accurate VLE calculations. The conventional Berthelot-Lorentz combining rule is expressed by Eq. (11) is applied in this study.

$$
\begin{equation*}
\varepsilon_{i j}=\sqrt{\varepsilon_{i} \varepsilon_{j}}\left(1-k_{i j}\right) \tag{11}
\end{equation*}
$$

where $k_{i j}$ is binary interaction parameter (Gross and Sadowski, 2001).

### 3.2.1. $\mathrm{MEA}-\mathrm{H}_{2} \mathrm{O}$ system

Binary VLE between MEA and water are shown in Figs. 4 and 5 with different association schemes for MEA and water, and binary interaction parameters are adjusted with the following objective function with experimental data from Cai et al. (1996)

$$
\begin{equation*}
\mathrm{OF}=\sum_{i=1}^{n p} \frac{\left|x_{i}^{\exp }-x_{i}^{c a l}\right|}{x_{i}^{\exp }}+\sum_{i=1}^{n p} \frac{\left|y_{i}^{\exp }-y_{i}^{c a l}\right|}{y_{i}^{\exp }} \tag{12}
\end{equation*}
$$

In order to use this objective function both Bubble-T and Dew-T calculations must be applied. This objective function is used to adjust general $k_{i j}$ even in small mole fractions in both phases. Binary interaction parameters between MEA and water are presented in Table 3. VLE results for experimental data from Cai et al. (1996) at pressures 66.66 kPa and 101.33 kPa are shown and compared in Figs. 4 and 5. In addition, experimental data from Park and Lee (1997) are predicted with PC-SAFT EOS. The AAD\% and $k_{i j}$ from Table 3 show that the 3B

![](https://cdn.mathpix.com/cropped/619bc8ff-8258-4e1f-b381-2302a83268ab-07.jpg?height=895&width=846&top_left_y=223&top_left_x=1064)
Fig. 7 - Isothermal pressure-composition $\mathrm{P}-\mathrm{x}$ slices of the vapor-liquid equilibrium of $\mathrm{CO}_{2}-\mathrm{H}_{2} \mathrm{O}$. The symbols correspond to the experimental data from Campos et al. (2009) at $\mathrm{T}=298.2 \mathrm{~K}, \mathrm{~T}=303.2 \mathrm{~K}, \mathrm{~T}=313.2 \mathrm{~K}$ and $\mathrm{T}=323.2 \mathrm{~K}$ and the solid curves represent to the PC-SAFT EOS predictions ( $k_{i j}=0$ ).

association scheme for MEA indicates the best agreement with experimental data for both schemes of water with respect to other schemes. Therefore, the 3B scheme for MEA is selected. The binary VLE results with both the 2 B and 4 C association schemes for water, and $2 \mathrm{~B}, 3 \mathrm{~B}$ and 4 C association schemes for MEA, are presented in Figs 4 and 5, respectively.

Experimental mixing enthalpies at 298.15 K and 342.15 K from Posey (1996) and Touhara et al. (1982) are predicted with

![](https://cdn.mathpix.com/cropped/619bc8ff-8258-4e1f-b381-2302a83268ab-07.jpg?height=778&width=846&top_left_y=1774&top_left_x=1064)
Fig. 8 - Comparison of the experimental data for species concentration in $\mathrm{CO}_{2}-\mathrm{MEA}-\mathrm{H}_{2} \mathrm{O}$ and the prediction of ideal Smith-Missen algorithm at $\mathbf{T} \boldsymbol{=} \mathbf{2 9 3 . 1 5 ~ K}$. MEA concentration is $30 \mathrm{wt} . \%$. Symbols represent experimental data from Jakobsen et al. (2005).

Table 4 - Equilibrium constants correlations used for modeling in Eq. (17).
| Reaction Equilibrium constant no. | A | B | C | D | Temperature range (K) | Ref. |
| :--- | :--- | :--- | :--- | :--- | :--- | :--- |
| (R1) | 132.899 | -13445.9 | 22.4773 | 0 | 273.15-498.15 | Edwards et al. (1978) |
| (R2) | 231.465 | -12092.10 | -36.7816 | 0 | 273.15-498.15 | Edwards et al. (1978) |
| (R3) | 216.049 | -12431.70 | -35.4819 | 0 | 273.15-498.15 | Edwards et al. (1978) |
| (R4) | -1.8652 | -1545.3 | 0 | 0 | 293.15-323.15 | Tong et al. (2012) |
| (R5) | 2.1211 | -8189.38 | 0 | -0.007484 | 273.15-323.15 | Bates and Pinching (1951) |
| ${ }^{\mathrm{a}}$ This equilibrium constant was converted to mole fraction based. |  |  |  |  |  |  |


the 3 B scheme for MEA and the 4 C scheme for water, and results are shown in Fig. 6.

### 3.2.2. $\quad \mathrm{CO}_{2}-\mathrm{H}_{2} \mathrm{O}$ system

Both 4 C and 2 B association schemes for water were applied in calculating the VLE $\mathrm{CO}_{2}-\mathrm{H}_{2} \mathrm{O}$ mixture in our previous work (Pahlavanzadeh and Fakouri Baygi, 2013), and the results show that the 4C association scheme for water includes better results respect with the 2 B association scheme in binary VLE of $\mathrm{CO}_{2}-\mathrm{H}_{2} \mathrm{O}$ system. In addition, the prediction of $\mathrm{CO}_{2}$ solubility in water at temperature from 298.2 K to 313.2 K is presented in Fig. 7, and is compared with experimental data from Campos et al. (2009) As one can see, PC-SAFT EOS shows great prediction for $\mathrm{CO}_{2}-\mathrm{H}_{2} \mathrm{O}$ system at very low $\mathrm{CO}_{2}$ mole fraction (less
than $10^{-3}$ ) in Fig. 7. These small errors in the range of $\mathrm{CO}_{2}$ solubility are proper to model $\mathrm{CO}_{2}-\mathrm{MEA}-\mathrm{H}_{2} \mathrm{O}$ system.

### 3.3. Ternary systems $\left(\mathrm{CO}_{2}-\mathrm{MEA}-\mathrm{H}_{2} \mathrm{O}\right)$

The procedure of finding $\mathrm{CO}_{2}$ solubility in this system is described in Fig. 1. This ternary system includes both chemical equilibrium and multicomponent system phase equilibrium. The liquid phase consists of both molecular species and ionic species, and the gas phase consists of only molecular species.

In this work, no regression is applied for VLE prediction of the $\mathrm{CO}_{2}$-MEA- $\mathrm{H}_{2} \mathrm{O}$ system by PC-SAFT EOS. This means the binary interaction parameters ( $k_{\mathrm{ij}}$ ) between ( $\mathrm{CO}_{2}-\mathrm{H}_{2} \mathrm{O}$ ) and ( $\mathrm{CO}_{2}-\mathrm{MEA}$ ) have been set to zero, and the binary

Table 5 - Comparison between experimental data and model predictions for $\mathrm{CO}_{2}$ partial pressure of the $\mathrm{CO}_{2}-\mathrm{MEA}-\mathrm{H}_{2} \mathrm{O}$ system.
| Source | Temperature, T (K) | $\mathrm{CO}_{2}$ partial pressure, $p$ (kPa) | MEA mole fraction | $\mathrm{CO}_{2}$ loading, $\alpha$ | e-NRTL by Zhang et al. (2011) |  | PC-SAFT EOS by this work |  |
| :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- |
|  |  |  |  |  | $n p$ | ${ }^{\mathrm{b}}$ AAD\% | np | ${ }^{\mathrm{b}}$ AAD\% |
| Jones et al. (1959) | 313-413 | 0.073-918 | 0.05 | 0.076-0.728 | 48 | 32.2 | 50 | 26.17 |
| Lee et al. (1974) | 313-373 | 1.2-6616 | 0.05 | 0.139-1.19 | 45 | 44.4 | 45 | 36.51 |
| Lawson and Garst (1976) | 313-413 | 1.3-2750 | 0.05-0.11 | 0.11-0.929 | 24 | 30.7 | 21 | 34.30 |
| Lee et al. (1976) | 298-393 | 0.1-10,000 | 0.02-0.11 | 0.065-2.152 | 256 | 50.5 | 253 | 34.76 |
| Isaacs et al. (1980) | 353-373 | 0.009-1.75 | 0.05 | 0.0368-0.315 | 19 | 112 | 16 | 22.66 |
| Austgen et al. (1991) | 313-353 | 0.093-228 | 0.05 | 0.266-0.698 | 8 | 42.6 | 8 | 11.97 |
| Shen and Li (1992) | 313 | 1.57-2550 | 0.05 | 0.561-1.049 | 13 | 35.5 | 13 | 30.92 |
| Dawodu and Meisen (1994) | 373 | 455-3863 | 0.09 | 0.541-0.723 | 5 | 43.2 | 5 | 17.05 |
| Song et al. (1996) | 312 | 3.1-2359 | 0.11 | 0.49-1.061 | 10 | 74.4 | 10 | 66.87 |
| Jane and Li (1997) | 353 | 3.57-121.8 | 0.05 | 0.363-0.58 | 7 | 37.7 | 7 | 30.91 |
| Park et al. (1997) | 313 | 3.5-2092 | 0.05 | 0.512-1.046 | 7 | 28.6 | 7 | 21.25 |
| Mathonat et al. (1998) | 313-393 | 2000-10000 | 0.11 | 0.55-1.07 | 9 | 36.0 | 9 | 35.62 |
| Park et al. (2002) | 313 | 2.6-2189 | 0.03-0.05 | 0.478-1.068 | 13 | 49.8 | 13 | 50.65 |
| Tong et al. (2012) | 313-393 | 3.95-408.17 | 0.11 | 0.211-0.748 | - | - | 21 | 53.53 |
| ${ }^{\mathrm{a}}$ Hilliard (2008) | 313-333 | 0.005-50.2 | 0.06-0.16 | 0.114-0.591 | 55 | 35.5 | 42 | 34.97 |
| ${ }^{\mathrm{a}}$ Jou et al. (1995) | 273-423 | 0.0012-19954 | 0.11 | 0.002-1.324 | 124 | 33.5 | 100 | 43.16 |
| ${ }^{\text {a }}$ Ma'mun et al. (2005) | 393 | 7.354-191.9 | 0.11 | 0.155-0.4182 | 19 | 13.5 | 19 | 21.03 |
| ${ }^{\mathrm{a}}$ Xu and Rochelle (2011) | 373-443 | 12-1626 | 0.11 | 0.303-0.52 | 63 | 28.0 | 52 | 46.88 |
| Overall |  |  |  |  | 725 | 42.29 | 691 | 36.42 |


[^2]![](https://cdn.mathpix.com/cropped/619bc8ff-8258-4e1f-b381-2302a83268ab-09.jpg?height=768&width=846&top_left_y=228&top_left_x=173)
Fig. 9 - Experimental data from Jones et al. (1959) in $w_{\text {MEA }}=15.3 \mathrm{wt} . \%$ at temperatures $313.15-413.15 \mathrm{~K}$ (symbols) and $\mathrm{CO}_{2}$ partial pressures prediction by PC-SAFT EOS (solid lines) of $\mathrm{CO}_{2}-\mathrm{MEA}-\mathrm{H}_{2} \mathrm{O}$ system.

interaction parameter between ( $\mathrm{MEA}-\mathrm{H}_{2} \mathrm{O}$ ) has been set to -0.0520 . Meanwhile, interaction energies between $\mathrm{CO}_{2}$ and MEA were considered in this study (see Table 2 and Eq. (11)), but no regression was applied to modify their interaction energies. In the other words, if $k_{i j}$ is set to zero between molecules $i$ and $j$ in PC-SAFT EOS, it does not mean that there is no interaction between molecules $i$ and $j$; rather, it means that PC-SAFT is able to predict interaction energies between molecules $i$ and $j$ without any correction.

The chemical reactions that take place in the liquid phase for $\mathrm{CO}_{2}-\mathrm{MEA}-\mathrm{H}_{2} \mathrm{O}$ system can be expressed as:

$$
\begin{equation*}
2 \mathrm{H}_{2} \mathrm{O} \leftrightarrow \mathrm{H}_{3} \mathrm{O}^{+}+\mathrm{OH}^{-} \tag{R1}
\end{equation*}
$$

![](https://cdn.mathpix.com/cropped/619bc8ff-8258-4e1f-b381-2302a83268ab-09.jpg?height=764&width=844&top_left_y=1827&top_left_x=173)
Fig. 10 - Experimental data from Lee et al. (1976) in $m_{\text {MEA }}=3.75$ at temperatures $313.15-413.15 \mathrm{~K}$ (symbols) and $\mathrm{CO}_{2}$ partial pressures prediction by PC-SAFT EOS (solid lines) of $\mathrm{CO}_{2}-\mathrm{MEA}-\mathrm{H}_{2} \mathrm{O}$ system.

![](https://cdn.mathpix.com/cropped/619bc8ff-8258-4e1f-b381-2302a83268ab-09.jpg?height=803&width=851&top_left_y=223&top_left_x=1064)
Fig. 11 - Comparison between experimental data from Jou et al. (1995), Tong et al. (2012) and Ma'mun et al. (2005) in $w_{\text {MEA }}=30 \mathrm{wt} . \%$ at temperatures 313 and 393 K (symbols) and $\mathrm{CO}_{2}$ partial pressures prediction by PC-SAFT EOS (solid lines) of $\mathrm{CO}_{2}$-MEA- $\mathrm{H}_{2} \mathrm{O}$ system.

$$
\begin{equation*}
\mathrm{CO}_{2}+2 \mathrm{H}_{2} \mathrm{O} \leftrightarrow \mathrm{H}_{3} \mathrm{O}^{+}+\mathrm{HCO}_{3}{ }^{-} \tag{R2}
\end{equation*}
$$

$$
\begin{equation*}
\mathrm{HCO}_{3}{ }^{-}+\mathrm{H}_{2} \mathrm{O} \leftrightarrow \mathrm{H}_{3} \mathrm{O}^{+}+\mathrm{CO}_{3}{ }^{2-} \tag{R3}
\end{equation*}
$$

$$
\begin{equation*}
\mathrm{MEACOO}^{-}+\mathrm{H}_{2} \mathrm{O} \leftrightarrow \mathrm{HCO}_{3}{ }^{-}+\mathrm{MEA} \tag{R4}
\end{equation*}
$$

$$
\begin{equation*}
\mathrm{MEAH}^{+}+\mathrm{H}_{2} \mathrm{O} \leftrightarrow \mathrm{H}_{3} \mathrm{O}^{+}+\mathrm{MEA} \tag{R5}
\end{equation*}
$$

Charge balance is expressed as:

$$
\begin{array}{r}
{\left[\mathrm{H}_{3} \mathrm{O}^{+}\right]+\left[\mathrm{MEAH}^{+}\right]-\left[\mathrm{MEACOO}^{-}\right]-\left[\mathrm{HCO}_{3}^{-}\right]-2\left[\mathrm{CO}_{3}^{2-}\right]-\left[\mathrm{OH}^{-}\right]} \\
=0 \tag{13}
\end{array}
$$

![](https://cdn.mathpix.com/cropped/619bc8ff-8258-4e1f-b381-2302a83268ab-09.jpg?height=766&width=847&top_left_y=1827&top_left_x=1066)
Fig. 12 - Comparison between experimental data from Park et al. (1997) and Shen and Li (1992) in $\mathbf{w}_{\text {MEA }}=\mathbf{1 5 . 3} \mathbf{~ w t} \boldsymbol{.} \boldsymbol{\%}$ at $\mathrm{T}=313.15 \mathrm{~K}$ (symbols) and $\mathrm{CO}_{2}$ partial pressures prediction by PC-SAFT EOS (solid lines) of $\mathrm{CO}_{2}$-MEA- $\mathrm{H}_{2} \mathrm{O}$ system.

![](https://cdn.mathpix.com/cropped/619bc8ff-8258-4e1f-b381-2302a83268ab-10.jpg?height=776&width=846&top_left_y=223&top_left_x=157)
Fig. 13 - Comparison between experimental data from Austgen et al. (1991) in $m_{\text {MEA }}=2.5 \mathrm{kmol} \mathrm{m}^{-3}$ at temperatures 313.15 and 353.15 K (symbols) and $\mathrm{CO}_{2}$ partial pressures prediction by PC-SAFT EOS (solid lines) of $\mathrm{CO}_{2}$-MEA- $\mathrm{H}_{2} \mathrm{O}$ system.

Three mass balances are also expressed as:

$$
\begin{equation*}
\mathrm{M}_{\mathrm{t}}=[\mathrm{MEA}]+\left[\mathrm{MEAH}^{+}\right]+\left[\mathrm{MEACOO}^{-}\right] \tag{14}
\end{equation*}
$$

$$
\begin{equation*}
\alpha \mathrm{M}_{\mathrm{t}}=\left[\mathrm{CO}_{2}\right]-\left[\mathrm{HCO}_{3}^{-}\right]+\left[\mathrm{CO}_{3}^{2-}\right]+\left[\mathrm{MEACOO}{ }^{-}\right] \tag{15}
\end{equation*}
$$

$$
\begin{equation*}
\mathrm{W}=\left[\mathrm{H}_{2} \mathrm{O}\right]+\left[\mathrm{HCO}_{3}^{-}\right]+\left[\mathrm{CO}_{3}^{2-}\right]+\left[\mathrm{H}_{3} \mathrm{O}^{+}\right]+\left[\mathrm{OH}^{-}\right] \tag{16}
\end{equation*}
$$

where $\alpha, \mathrm{M}_{\mathrm{t}}$, and W are loading of $\mathrm{CO}_{2}$ in aqueous MEA, total moles of MEA, and total amount of $\mathrm{H}_{2} \mathrm{O}$ in liquid phase, respectively.

Chemical reactions (R1)-(R5) and linear equations (13)-(16) are solved by the Lagrange method simultaneously, which was introduced by Smith-Missen (Austgen, 1989). In order to find the concentration of species in chemical equilibrium, it is assumed their activity equal to their mole fraction $\left(a_{i}=x_{i}\right)$ for all species, and the effects of the ion pairs are neglected in equilibrium. It is necessary to say that this assumption $\left(a_{i}=x_{i}\right)$ has been applied in several similar works (Nasrifar and Tafazzol, 2010; Gabrielsen et al., 2005). To check these assumed results with experimental data, the concentration of some molecular and ionic species are compared with NMR experimental speciation data from Jakobsen et al. (2005) in Fig. 8.

The chemical equilibrium constants for reactions (R1)-(R5) are calculated by Eq. (17):

$$
\begin{equation*}
\ln K_{X}=A+\frac{B}{(T / K)}+C \ln (T / K)+D(T / K) \tag{17}
\end{equation*}
$$

where the constants of Eq. (17) for reactions (R1)-(R5) are presented in Table 4.

In this work, the constants of Edwards et al. (1978) are used for reactions (R1)-(R3), and for achievement the best ideal chemical equilibrium (because of this assumption $\left(a_{i}=x_{i}\right)$ ) the constants of Tong et al. (2012) and Bates and Pinching (1951) are used for reactions (R4) and (R5), respectively. Chemical equilibrium for reaction (R4) was reported in molality based
by Edwards et al. (1978), and then these constants are converted to the mole fraction based on the following equation (Austgen, 1989).

$$
\begin{equation*}
\ln K_{x}=\ln K_{m}-\ln \left(\frac{1}{0.01802}\right) \sum_{i \neq \mathrm{H}_{2} \mathrm{O}} v_{i} \tag{18}
\end{equation*}
$$

where $K_{x}, K_{m}$ and $v$ are chemical equilibrium constant in mole fraction based, chemical equilibrium constant in molality based and stoichiometric coefficients in reaction, respectively.

Results of PC-SAFT EOS predictions compared to different experimental data are presented in Table 5. The results of the e-NRTL model as implemented by Zhang et al. (2011) for the similar systems have been illustrated in Table 5 for comparison reasons. It is worth it to point out that although the PC-SAFT EOS uses no regression in VLE calculation, the model predictions are acceptable within the average relative deviation in $\mathrm{CO}_{2}$ partial pressure for many experimental data are shown in Table 5.

Fig. 9 shows that the prediction of the $\mathrm{CO}_{2}$ partial pressure in $15.3 \mathrm{wt} . \%$ and at wide range of temperatures between 313.15 and 413.15 K , with experimental data from Jones et al. (1959) is reliable. Fig. 10 exhibits that PC-SAFT is able to anticipate the experimental $\mathrm{CO}_{2}$ partial pressure data from Lee et al. (1976) in molality 3.75 MEA and at temperatures from 313.15 to 393.15 K acceptably. Model predictions and experimental data from Jou et al. (1995), Tong et al. (2012) and Ma'mun et al. (2005) in $30 \mathrm{wt} . \%$ MEA at two temperatures 313 K and 393 K are compared in Fig. 11. One can see from Fig. 11 that PCSAFT EOS has acceptable prediction with experimental data at these temperatures. There are also some deviations between experimental data at 393.15 K . While experimental data from Tong et al. (2012) and Ma'mun et al. (2005) were reported in similar MEA concentration and temprature ranges. Systematic deviations can be observed at the high $\mathrm{CO}_{2}$ loading region (approximately more than 0.9 ) under the different MEA concentration and temperature ranges in Figs. 10-12. Comparison of experimental data from Austgen et al. (1991) and PC-SAFT EOS shows good performance at the average loading region (0.266-0.698) in Fig. 13.

## 4. Conclusion

In this work, PC-SAFT EOS is applied for the prediction of a $\mathrm{CO}_{2}$-MEA- $\mathrm{H}_{2} \mathrm{O}$ system devoid of using any regression in the ternary system in temperature range $298-443 \mathrm{~K}$ and CO2 partial pressure range $0.001-20,000 \mathrm{kPa}$. Because of the importance of association schemes in parameter estimation in associated EOSs, several association schemes are checked for MEA, and the results indicate that the 3 B association scheme for MEA shows more conformity with binary VLE experimental data more than other association scheme. Therefore, both the 3 B association scheme for MEA and the 4 C association scheme for water (due to its consistency in $\mathrm{CO}_{2}-\mathrm{H}_{2} \mathrm{O}$ system) are applied. Binary interaction parameters between MEA and $\mathrm{H}_{2} \mathrm{O}$ in different association schemes are adjusted to correct VLE prediction of this binary system. In addition, excess enthalpy of the MEA and the water system is predicted by the binary interaction parameter that was adjusted in VLE calculation. In summary, PC-SAFT EOS performs an acceptable prediction in the ternary system of $\mathrm{CO}_{2}-\mathrm{MEA}-\mathrm{H}_{2} \mathrm{O}$ in a wide range of temperatures and $\mathrm{CO}_{2}$ partial pressure without using any regression in multicomponent system, and its average
absolute partial pressure deviation is about $36.42 \%$ among 691 data points.

## Acknowledgments

Authors are thankful to Hossein Eghbali for providing VLE data set of $\mathrm{CO}_{2}-\mathrm{MEA}-\mathrm{H}_{2} \mathrm{O}$ system.

## References

Austgen Jr., D., 1989. A Model of Vapor-Liquid Equilibria for Acid Gas-Alkanolamine-Water Systems. Texas University, Austin, TX, USA.
Austgen, D.M., Rochelle, G.T., Peng, X., Chen, C.C., 1989. Ind. Eng. Chem. Res. 28, 1060-1073.
Austgen, D.M., Rochelle, G.T., Chen, C.C., 1991. Ind. Eng. Chem. Res. 30, 543-555.
Avlund, A.S., 2011. Extension of Association Models to Complex Chemicals. Department of Chemical and Biochemical Engineering, Technical University of Denmark.
Avlund, A.S., Kontogeorgis, G.M., Michelsen, M.L., 2008. Ind. Eng. Chem. Res. 47, 7441-7446.
Bates, R.G., Pinching, G.D., 1951. J. Res. Natl. Bur. Stand. 46, 349352.
Cai, Z., Xie, R., Wu, Z., 1996. J. Chem. Eng. Data 41, 1101-1103.
Campos, C.E.P.S., Villardi, H.G.D.A., Pessoa, F.L.P., Uller, A.M.C., 2009. J. Chem. Eng. Data 54, 2881-2886.

Dawodu, O.F., Meisen, A., 1994. J. Chem. Eng. Data 39, 548-552.
Diamantonis, N.I., Economou, I.G., 2011. Energy Fuels 25, 3334-3343.
Edwards, T.J., Maurer, G., Newman, J., Prausnitz, J.M., 1978. AIChE J. 24, 966-976.

Faramarzi, L., Kontogeorgis, G.M., Thomsen, K., Stenby, E.H., 2009. Fluid Phase Equilib. 282, 121-132.
Gabrielsen, J., Michelsen, M.L., Stenby, E.H., Kontogeorgis, G.M., 2005. Ind. Eng. Chem. Res. 44, 3348-3354.

Gross, J., Sadowski, G., 2001. Ind. Eng. Chem. Res. 40, 1244-1260.
Gross, J., Sadowski, G., 2002. Ind. Eng. Chem. Res. 41, 5510-5515.
Hilliard, M.D., 2008. A Predictive Thermodynamic Model for an Aqueous Blend of Potassium Carbonate, Piperazine, and Monoethanolamine for Carbon Dioxide Capture from Flue Gas. University of Texas at Austin.
Huang, S.H., Radosz, M., 1990. Ind. Eng. Chem. Res. 29, 2284-2294.
Isaacs, E.E., Otto, F.D., Mather, A.E., 1980. J. Chem. Eng. Data 25, 118-120.

Jakobsen, J.P., Krane, J., Svendsen, H.F., 2005. Ind. Eng. Chem. Res. 44, 9894-9903.
Jane, I.S., Li, M.-H., 1997. J. Chem. Eng. Data 42, 98-105.
Jones, J.H., Froning, H.R., Claytor, E.E., 1959. J. Chem. Eng. Data 4, 85-92.
Jou, F.-Y., Mather, A.E., Otto, F.D., 1995. Can. J. Chem. Eng. 73, 140-147.
Lawson, J.D., Garst, A.W., 1976. J. Chem. Eng. Data 21, 20-30.
Lee, J.I., Otto, F.D., Mather, A.E., 1974. Can. J. Chem. Eng. 52, 803-805.
Lee, J.I., Otto, F.D., Mather, A.E., 1976. J. Appl. Chem. Biotechnol. 26, 541-549.
Mac Dowell, N., Llovell, F., Adjiman, C.S., Jackson, G., Galindo, A., 2009. Ind. Eng. Chem. Res. 49, 1883-1899.

Ma'mun, S., Nilsen, R., Svendsen, H.F., Juliussen, O., 2005. J. Chem. Eng. Data 50, 630-634.
Mathonat, C., Majer, V., Mather, A.E., Grolier, J.P.E., 1998. Ind. Eng. Chem. Res. 37, 4136-4141.
Nasrifar, K., Tafazzol, A.H., 2010. Ind. Eng. Chem. Res. 49, 7620-7630.
Pahlavanzadeh, H., Fakouri Baygi, S., 2013. J. Chem. Thermodyn. 59, 214-221.
Park, S.-B., Lee, H., 1997. Korean J. Chem. Eng. 14, 146-148.
Park, S.-B., Shim, C.-S., Lee, H., Lee, K.-H., 1997. Fluid Phase Equilib. 134, 141-149.
Park, J.-Y., Yoon, S.J., Lee, H., Yoon, J.-H., Shim, J.-G., Lee, J.K., Min, B.-Y., Eum, H.-M., Kang, M.C., 2002. Fluid Phase Equilib. 202, 359-366.
Posey, M.L., 1996. Thermodynamic Model for Acid Gas Loaded Aqueous Alkanolamine Solutions. University of Texas at Austin.
Posey, M.L., Rochelle, G.T., 1997. Ind. Eng. Chem. Res. 36, 3944-3953.
Rodriguez, J., Mac Dowell, N., Llovell, F., Adjiman, C.S., Jackson, G., Galindo, A., 2012. Mol. Phys. 110, 1325-1348.
Shen, K.P., Li, M.H., 1992. J. Chem. Eng. Data 37, 96-100.
Song, J.-H., Yoon, J.-H., Lee, H., Lee, K.-H., 1996. J. Chem. Eng. Data 41, 497-499.
Tong, D., Trusler, J.P.M., Maitland, G.C., Gibbins, J., Fennell, P.S., 2012. Int. J. Greenh. Gas Control 6, 37-47.

Touhara, H., Okazaki, S., Okino, F., Tanaka, H., Ikari, K., Nakanishi, K., 1982. J. Chem. Thermodyn. 14, 145-156.

Xu, Q., Rochelle, G., 2011. Energy Proc. 4, 117-124.
Zhang, Y., Que, H., Chen, C.-C., 2011. Fluid Phase Equilib. 311, 67-75.
Zoghi, A.T., Feyzi, F., Dehghani, M.R., 2012. Ind. Eng. Chem. Res. 51, 9875-9885.


[^0]:    * Corresponding author. Tel.: +98 21 82883312; fax: +98 2182883381.

    E-mail address: pahlavzh@modares.ac.ir (H. Pahlavanzadeh).
    Available online 30 July 2014

[^1]:    ${ }^{\text {a }}$ AAD $\%=100 / n p \times \Sigma\left|\psi^{\text {exp }}-\psi^{\text {cal }}\right| / \psi^{\text {exp }}, \psi$ is $x$ and $y$.

[^2]:    ${ }^{\mathrm{a}}$ These data sets are used in e-NRTL regression by Zhang et al. (2011).
    ${ }^{\mathrm{b}} \mathrm{AAD} \%=100 / n p \times \sum\left|p_{\mathrm{CO}_{2}}^{\mathrm{exp}}-p_{\mathrm{CO}_{2}}^{\mathrm{cal}}\right| / p_{\mathrm{CO}_{2}}^{\mathrm{exp}}$.

