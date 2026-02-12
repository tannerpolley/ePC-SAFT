# Application of the Perturbed-Chain SAFT Equation of State to Associating Systems 

J oachim Gross ${ }^{\boldsymbol{\dagger} \boldsymbol{,} \boldsymbol{\ddagger}}$ and Gabriele Sadowski ${ }^{\boldsymbol{*}, \boldsymbol{\S}}$<br>Fachgebiet Thermodynamik und Thermische Verfahrenstechnik, Technische Universität Berlin, 10623 Berlin, Germany, and Lehrstuhl für Thermodynamik, Universität Dortmund, Emil-FiggeStrasse 70, 44227 Dortmund, Germany

Downloaded via BRIGHAM YOUNG UNIV on May 18, 2024 at 03:34:25 (UTC). See https://pubs.acs.org/sharingguidelines for options on how to legitimately share published articles.


#### Abstract

The perturbed-chain SAFT (PC-SAFT) equation of state is applied to pure associating components as well as to vapor-liquid and liquid-liquid equilibria of binary mixtures of associating substances. For these substances, the PC-SAFT equation of state requires five purecomponent parameters, two of which characterize the association. The pure-component parameters were identified for 18 associating substances by correlating vapor pressure and liquid density data. A comparison to an earlier version of SAFT confirms the good results for pure substances. When only one associating compound is present in a mixture, the PC-SAFT equation of state does not require mixing rules for the association term. Using one binary interaction parameter $\mathrm{k}_{\mathrm{ij}}$ for the dispersion term only, the model was applied to azeotropic and nonazeotropic vapor-liquid equilibria at low and at high pressures, as well as to liquid-liquid equilibria. Simple mixing and combining rules were adopted for mixtures with more than one associating compound, introducing no additional binary interaction parameter. The simultaneous description of liquidliquid and vapor-liquid equilibrium was also possible with a single $\mathrm{k}_{\mathrm{ij}}$ parameter.


## 1. Introduction

Modeling the phase equilibrium and thermodynamic properties of systems in which molecules exhibit associating interactions remains a challenging problem in chemical industry. The large number of recent studies on equations of state applied to associating systems is an indication of this circumstance. Detailed reviews are given, for example, by Müller and Gubbins, ${ }^{1}$ E conomou and Donohue, ${ }^{2}$ and Wei and Sadus. ${ }^{3}$

It is beyond dispute that associating interactions that contribute significantly to the overall intermolecular forces have to be explicitly taken into account in molecular theories and, further, that physically based equations of state have the potential for correlating and someday even predicting the behavior of such systems. Progress toward this goal has been achieved by applying principles of statistical mechanics. One step along this road was the development of the statistical association fluid theory (SAFT) by Chapman et al. ${ }^{4,5}$ The perturbation theory of Wertheim ${ }^{6-9}$ is suitable for highly directional interactions and provided the basis for the development of a fairly simple association term. ${ }^{10-12}$ Numerous subsequent investigations have aimed at

[^0]numerical aspects of solving the appropriate equations or at applications of the theory to real substances (see ref 1).

Huang and Radosz applied the SAFT framework developing a modified SAFT equation of state and determined pure-component parameters for an array of components of industrial relevance. ${ }^{13,14}$ This particular version was subsequently often used for academic and for industry applications. ${ }^{1-3}$

Many SAFT models account for the chainlike shape of the molecules only in the repulsive contribution of the equation of state (hard-chain contribution). In contrast, Gross and Sadowski devel oped the perturbedchain SAFT (PC-SAFT) equation of state, where the chain-length dependence of the attractive (dispersive) interactions is also taken into account. The PC-SAFT model was shown to be suitable for pure components and mixtures of solvents and gases, ${ }^{15}$ as well as for liquid-liquid and vapor-liquid equilibria of polymer systems. ${ }^{16}$ The PC-SAFT equation of state is implemented in a common process simulation tool (by AspenTech Inc.). It is further available as coded subroutines ${ }^{17}$ and is incorporated in the engineering software PE (Phase Equilibria, developed and maintained by Pfohl, Petkov, and Brunner ${ }^{18}$ ), which is freely available for scientific and industrial applications.

In this work, the PC-SAFT equation of state is applied to systems in which molecular association prevails. Pure components and mixtures are investigated. As in our earlier studies, we seek an evaluation of the PC-SAFT equation of state by comparing it to

Table 1. Pure-Component Parameters of the Perturbed-Chain SAFT Equation of State for Associating Components
| component $\mathrm{i}^{\mathrm{a}}$ | $\mathrm{M}_{\mathrm{i}}$ (g/mol) | $\mathrm{m}_{\mathrm{i}}$ | $\sigma_{\mathrm{i}}$ (Å) | $\epsilon_{\mathrm{i}} / \mathrm{k}$ (K) | $\kappa^{\mathrm{A}_{\mathrm{i}} \mathrm{B}_{\mathrm{i}}}$ | $\epsilon^{\mathrm{A}_{\mathrm{i}} \mathrm{B}_{\mathrm{i}}} / \mathrm{k}$ (K) | AAD\% |  | T range (K) | ref ${ }^{\text {b }}$ |
| :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- |
|  |  |  |  |  |  |  | psat | v |  |  |
| alkanols |  |  |  |  |  |  |  |  |  |  |
| methanol | 32.042 | 1.5255 | 3.2300 | 188.90 | 0.035176 | 2899.5 | 2.36 | 2.01 | 200-512 | 1, 2 |
| ethanol | 46.069 | 2.3827 | 3.1771 | 198.24 | 0.032384 | 2653.4 | 0.99 | 0.79 | 230-516 | 1, 2 |
| 1-propanol | 60.096 | 2.9997 | 3.2522 | 233.40 | 0.015268 | 2276.8 | 0.85 | 1.71 | 240-537 | 2 |
| 1-butanol | 74.123 | 2.7515 | 3.6139 | 259.59 | 0.006692 | 2544.6 | 1.78 | 1.63 | 184-563 | 1, 2 |
| 1-pentanol | 88.15 | 3.6260 | 3.4508 | 247.28 | 0.010319 | 2252.1 | 0.5 | 0.52 | 250-588 | 2 |
| 1-hexanol | 102.177 | 3.5146 | 3.6735 | 262.32 | 0.005747 | 2538.9 | 0.58 | 0.71 | 228-610 | 2 |
| 1-heptanol | 116.203 | 4.3985 | 3.5450 | 253.46 | 0.001155 | 2878.5 | 1.47 | 1.19 | 239-632 | 2 |
| 1-octanol | 130.23 | 4.3555 | 3.7145 | 262.74 | 0.002197 | 2754.8 | 1.28 | 0.92 | 258-652 | 2 |
| 1-nonanol | 144.257 | 4.6839 | 3.7292 | 263.64 | 0.001427 | 2941.9 | 0.86 | 1.23 | 268-670 | 2 |
| 2-propanol | 60.096 | 3.0929 | 3.2085 | 208.42 | 0.024675 | 2253.9 | 0.7 | 1.25 | 185-508 | 2 |
| 2-methyl-2-butanol | 88.15 | 2.5487 | 3.9053 | 266.01 | 0.001863 | 2618.8 | 0.19 | 0.44 | 264-545 | 2 |
| water |  |  |  |  |  |  |  |  |  |  |
| water | 18.015 | 1.0656 | 3.0007 | 366.51 | 0.034868 | 2500.7 | 1.88 | 6.83 | 273-647 | 1 |
| amines |  |  |  |  |  |  |  |  |  |  |
| methylamine | 31.06 | 2.3967 | 2.8906 | 214.94 | 0.095103 | 684.3 | 0.3 | 0.53 | 179-430 | 3 |
| ethylamine | 45.09 | 2.7046 | 3.1343 | 221.53 | 0.017275 | 854.7 | 0.41 | 1.15 | 210-456 | 3 |
| 1-propylamine | 59.11 | 2.4539 | 3.5347 | 250.52 | 0.022674 | 1028.1 | 0.32 | 0.7 | 220-497 | 3 |
| 2-propylamine | 59.11 | 2.5908 | 3.4777 | 231.80 | 0.021340 | 932.2 | 0.24 | 0.22 | 210-471 | 3 |
| aniline | 93.13 | 2.6607 | 3.7021 | 335.47 | 0.074883 | 1351.6 | 0.79 | 1.09 | 267-699 | 3 |
| acetic acid |  |  |  |  |  |  |  |  |  |  |
| acetic acid | 60.053 | 1.3403 | 3.8582 | 211.59 | 0.075550 | 3044.4 | 2.12 | 1.36 | 302-592 | 1, 2 |


${ }^{\mathrm{a}}$ Two association sites ( $N^{\text {site }}=2$ ) are assumed for all substances. ${ }^{\mathrm{b}}$ References: (1) VDI-Wärmeatlas; VDI-Gesellschaft Verfahrenstechnik und Chemieingenieurwesen (GVC): Düsseldorf, Germany, 1994. (2) Daubert, T. E.; Danner, R. P.; Sibul, H. M.; Stebbins, C. C. Physical and Thermodynamic Properties of Pure Chemicals: Data Compilation; Taylor \& Francis: Washington, DC, 1989. (3) Chao, J .; Gadalla, N. A. M.; Gammon, B. E. Marsh, K. N.; Rodgers, A. S.; Somayajulu, G. R.; Wilhoit, R. C. J. Phys. Chem. Ref. Data, 1990, 19 (6), 1547-1569.

Table 2. Binary Interaction Parameters $\mathbf{k}_{\mathbf{i j}}$ Correcting the Cross-Dispersive Interactions
| binary system | $\mathrm{k}_{\mathrm{ij}}{ }^{\mathrm{a}}$ |  |
| :--- | :--- | :--- |
|  | PC-SAFT | SAFT |
| methanol-isobutanol | 0.05 | 0.03 |
| methanol-cyclohexane | 0.051 | 0.044 |
| methanol-1-octanol | 0.020 | 0.034 |
| ethanol-n-butane | 0.028 | 0.021 |
| 1-propanol-benzene | 0.020 | 0.005 |
| 1-propanol-ethylbenzene | 0.023 | 0.011 |
| 2-propanol-benzene | 0.021 | 0.007 |
| 1-butanol-n-butane | 0.015 | 0.025 |
| 1-pentanol-benzene | 0.0135 | 0.0065 |
| water-1-pentanol | 0.016 | -0.069 |


${ }^{\mathrm{a}} \mathrm{k}_{\mathrm{ij}}=\mathrm{k}_{\mathrm{ji}}$, nonzero for $\mathrm{i} \neq \mathrm{j}$.
the SAFT version of Huang and Radosz ${ }^{13,14}$ (hereafter called SAFT for short).

## 2. Perturbed-Chain SAFT Equation of State

The PC-SAFT equation of state was derived and described in detail by Gross and Sadowski. ${ }^{15}$ In terms of the compressibility factor $Z$, the equation of state is given as the sum of the ideal gas contribution ( $Z^{\text {id }}=1$ ), the hard-chain term (hc), the dispersive part (disp), and the contribution due to association (assoc) according to

$$
\begin{equation*}
Z=1+Z^{\mathrm{hc}}+Z^{\mathrm{disp}}+Z^{\mathrm{assoc}} \tag{1}
\end{equation*}
$$

where $Z^{\text {assoc }}$ denotes the contribution of the associating interactions to the compressibility factor (details given in refs 5 and 14). The effect of multipole interactions (such as dipole-dipole forces) is not separately taken into account in eq 1 . This contribution will be the subject of further investigations.

Two pure-component parameters determine the associating interactions between the association site $\mathrm{A}_{\mathrm{i}}$ and $B_{i}$ of a pure component $i$ : the association energy
$\epsilon^{\mathrm{A}_{\mathrm{i}} \mathrm{B}_{\mathrm{i}}} / \mathrm{k}$ and the effective association volume $\kappa^{\mathrm{A}_{\mathrm{i}} \mathrm{B}_{\mathrm{i}}}$. It is generally far from trivial to obtain cross-association parameters between two different associating substances i and j. For many systems, however, approximate values for cross-association parameters can be determined from pure-component association parameters. Simple combining rules for cross-association were suggested by Wolbach and Sandler, ${ }^{19}$ from a consideration of gas-phase association in the low-pressure limit, as

$$
\begin{gather*}
\epsilon^{\mathrm{A}_{\mathrm{i}} \mathrm{~B}_{\mathrm{j}}}=\frac{1}{2}\left(\epsilon^{\mathrm{A}_{\mathrm{i}} \mathrm{~B}_{\mathrm{i}}}+\epsilon^{\mathrm{A}_{\mathrm{j}} \mathrm{~B}_{\mathrm{j}}}\right)  \tag{2}\\
\kappa^{\mathrm{A}_{\mathrm{i}} \mathrm{~B}_{\mathrm{j}} \mathrm{~A}}=\sqrt{\kappa^{\mathrm{A}_{\mathrm{i}} \mathrm{~B}_{\mathrm{i}}} \kappa^{\mathrm{A}_{\mathrm{j}} \mathrm{~B}_{\mathrm{j}}}}\left(\frac{\sqrt{\sigma_{\mathrm{ii}} \sigma_{\mathrm{jj}}}}{1 / 2\left(\sigma_{\mathrm{ii}}+\sigma_{\mathrm{jj}}\right)}\right)^{3} \tag{3}
\end{gather*}
$$

These rules will be used in this work. Within the scope of this investigation, no binary correction parameters are introduced in these combining rules. Hence, only one binary interaction parameter, $\mathrm{k}_{\mathrm{ij}}$, is permitted, which corrects the dispersive interactions ${ }^{15}$.

## 3. Pure-Component Parameters for Associating Compounds

All associating components are assigned two association sites (often referred to as the 2 B model ${ }^{13}$ ). Although this is a reasonable assumption for some species (such as alcohols), it is a considerable simplification for other compounds-in particular for water. A study of Economou and Tsonopoulos ${ }^{20}$ indicates that water is best represented with a four-site treatment, whereas Suresh and Elliott ${ }^{21}$ found the two-site model to perform at least as well. F or simplicity at this point, we follow the latter study. (Calculation results from the SAFT version of Huang and Radosz, which are shown later in this work

![](https://cdn.mathpix.com/cropped/8fadac5f-7817-40f6-806b-7bb45d749e64-3.jpg?height=770&width=819&top_left_y=163&top_left_x=193)
Figure 1. Saturated liquid and vapor densities for methanol, 1-pentanol, and 1-nonanol. Comparison of PC-SAFT (solid lines) and SAFT (dashed lines) to experimental data.

![](https://cdn.mathpix.com/cropped/8fadac5f-7817-40f6-806b-7bb45d749e64-3.jpg?height=652&width=828&top_left_y=1069&top_left_x=193)
Figure 2. Vapor-liquid equilibria of methanol-isobutane at $\mathrm{T}=100^{\circ} \mathrm{C}$. Comparison of experimental data ${ }^{23}$ to calculation results of PC-SAFT ( $\mathrm{k}_{\mathrm{ij}}=0.05$ ) and SAFT ( $\mathrm{k}_{\mathrm{ij}}=0.03$ ).

for comparison, were obtained using three association sites for water. ${ }^{13}$ )

The pure-component parameters of associating substances were identified by simultaneously fitting vapor pressure data and liquid density data. A total of five parameters must be adjusted for any associating component i, namely, the segment diameter $\sigma_{\mathrm{i}}$, the segment number $\mathrm{m}_{\mathrm{i}}$, the segment energy parameter $\epsilon_{\mathrm{i}} / \mathrm{k}$, the association energy $\epsilon^{\mathrm{A}_{\mathrm{i}} \mathrm{B}_{\mathrm{i}}} / \mathrm{k}$, and the effective association volume $\kappa^{\mathrm{A}} \mathrm{B}_{\mathrm{i}}$, Table 1 gives these parameters for some alkanols, amines, acetic acid, and water. Dimensionless absolute average deviations (AAD\%) of vapor pressures and liquid densities, along with the temperature range covered by the literature data, are also listed.

## 4. Results

This section presents correlation results for pure components as well as for vapor-liquid and liquidliquid equilibria of binary mixtures of associating substances obtained with the perturbed-chain SAFT

![](https://cdn.mathpix.com/cropped/8fadac5f-7817-40f6-806b-7bb45d749e64-3.jpg?height=654&width=826&top_left_y=161&top_left_x=1108)
Figure 3. Vapor-liquid equilibria of 1-propanol-ethylbenzene at two pressures. Comparison of experimental data ${ }^{24}$ to calculation results of PC-SAFT ( $\mathrm{k}_{\mathrm{ij}}=0.023$ ) and SAFT ( $\mathrm{k}_{\mathrm{ij}}=0.011$ ).

![](https://cdn.mathpix.com/cropped/8fadac5f-7817-40f6-806b-7bb45d749e64-3.jpg?height=676&width=834&top_left_y=947&top_left_x=1102)
Figure 4. Vapor-liquid equilibria of 1-pentanol-benzene at $\mathrm{T}=40^{\circ} \mathrm{C}$. Comparison of experimental data ${ }^{25}$ to calculation results of PC-SAFT ( $\mathrm{k}_{\mathrm{ij}}=0.0135$ ) and SAFT ( $\mathrm{k}_{\mathrm{ij}}=0.0065$ ).

equation of state. Mixtures containing one associating compound (self-association) are investigated in section 4.2; mixtures of two associating substances (crossassociation) are discussed in section 4.3. Table 1 gives the kij parameters of all mixtures investigated.
4.1. Pure Associating Substances. The absolute average deviations (AAD\%) of the PC-SAFT equation of state from experimental vapor pressure and density data in the majority of cases range below 2\% (Table 1). Water is an exception here. Average deviations of 6.8\% in liquid density and $1.88 \%$ in vapor pressure were obtained. Various assumptions inherent in the molecular model adopted here are specifically severe for water. Associating interactions of water lead to cyclic chainlike molecular topologies, which are currently not taken into account in the association term. ${ }^{2}$ Furthermore, the dipole-dipole interactions, are not separately taken into account in the present form of the PC-SAFT equation of state. These directional interactions are known to have a significant impact on the thermophysical behavior of water. ${ }^{22}$

Figure 1 shows the vapor-liquid equilibrium of three alkanols in a T- $\rho$ diagram. The densities of the coexisting phases are well-described by the PC-SAFT model.

![](https://cdn.mathpix.com/cropped/8fadac5f-7817-40f6-806b-7bb45d749e64-4.jpg?height=667&width=836&top_left_y=161&top_left_x=189)
Figure 5. Vapor-liquid equilibria of 2-propanol-benzene ( $\diamond$ ) and 1-propanol-benzene ( O ) at $\mathrm{T}=40^{\circ} \mathrm{C}$. Comparison of experimental data ${ }^{25}$ to calculation results of PC-SAFT (2-propanol-benzene, $\mathrm{k}_{\mathrm{ij}}=0.021$; 1-propanol-benzene, $\mathrm{k}_{\mathrm{ij}}=0.020$ ) and SAFT ( 2 -propanol-benzene, $\mathrm{k}_{\mathrm{ij}}=0.007$; 1-propanol-benzene, $\mathrm{k}_{\mathrm{ij}}=0.005$ ).

![](https://cdn.mathpix.com/cropped/8fadac5f-7817-40f6-806b-7bb45d749e64-4.jpg?height=680&width=836&top_left_y=1024&top_left_x=189)
Figure 6. Vapor-liquid equilibria of 1-butanol-n-butane at four temperatures. Comparison of experimental data ${ }^{26}$ to calculation results of PC-SAFT ( $\mathrm{k}_{\mathrm{ij}}=0.015$ ) and SAFT ( $\mathrm{k}_{\mathrm{ij}}=0.025$ ).

Results of the SAFT equation of state obtained with pure-component parameters determined by Huang and Radosz ${ }^{13}$ are given for comparison. It becomes apparent from Figure 1 that an improved description of dispersive interactions also leads to better results for associating compounds.
4.2. Mixtures Containing One Associating Substance. Mixtures with only one associating substance (self-association) do not require mixing rules for the association term. Calculations for such mixtures are suitable for investigating whether the association term contributes to the total Helmholtz energy (in eq 1) at the right order of magnitude.

Figures $2-5$ show vapor-liquid equilibria of binary, self-associating mixtures in which both compounds are below their critical point. Most of the systems displayed show azeotropic behavior. Mixtures in which one compound is above its critical point are given in Figures 6 and 7. Some of the improvement of the PC-SAFT equation of state over the SAFT model is actually due to an improved description of pure components. Figure 7, for example, shows deficiencies of the SAFT equation

![](https://cdn.mathpix.com/cropped/8fadac5f-7817-40f6-806b-7bb45d749e64-4.jpg?height=732&width=838&top_left_y=158&top_left_x=1100)
Figure 7. Vapor-liquid equilibria of ethanol-n-butane at four temperatures. Comparison of experimental data ${ }^{26}$ to calculation results of PC-SAFT ( $\mathrm{k}_{\mathrm{ij}}=0.028$ ) and SAFT ( $\mathrm{k}_{\mathrm{ij}}=0.021$ ).

![](https://cdn.mathpix.com/cropped/8fadac5f-7817-40f6-806b-7bb45d749e64-4.jpg?height=1017&width=838&top_left_y=1022&top_left_x=1100)
Figure 8. Isobaric vapor-liquid and liquid-liquid equilibria of methanol-cyclohexane at $\mathrm{P}=1.013$ bar. Comparison of experimental data ${ }^{27-30}$ to calculation results of PC-SAFT ( $\mathrm{k}_{\mathrm{ij}}=0.051$ ) and SAFT ( $\mathrm{k}_{\mathrm{ij}}=0.044$ ).

of state in reproducing the vapor pressure of butane at higher temperatures ( $\mathrm{T}=140^{\circ} \mathrm{C}$ ). Furthermore, a subcritical system is predicted at $160^{\circ} \mathrm{C}$, although the critical point of butane is about 8 K lower.

The phase behavior of a mixture of methanol and cyclohexane at normal pressure is displayed in Figure 8. This system exhibits an azeotropic vapor-liquid equilibrium at higher temperatures and shows a liquidliquid equilibrium at lower temperatures. The $\mathrm{k}_{\mathrm{ij}}$ parameters of the PC-SAFT equation of state and of the SAFT model were adjusted to the highest available liquid-liquid demixing temperature (close to the critical

![](https://cdn.mathpix.com/cropped/8fadac5f-7817-40f6-806b-7bb45d749e64-5.jpg?height=715&width=830&top_left_y=158&top_left_x=191)
Figure 9. Isobaric vapor-liquid equilibria of methanol-1-octanol at $\mathrm{P}=1.013$ bar. Comparison of experimental ${ }^{31}$ data to calculation results of PC-SAFT ( $\mathrm{k}_{\mathrm{ij}}=0.020$ ) and SAFT ( $\mathrm{k}_{\mathrm{ij}}=0.034$ ).

![](https://cdn.mathpix.com/cropped/8fadac5f-7817-40f6-806b-7bb45d749e64-5.jpg?height=708&width=841&top_left_y=1001&top_left_x=184)
Figure 10. Isobaric, heteroazeotropic vapor-liquid and liquidliquid equilibria of water-1-pentanol at $\mathrm{P}=1.013$ bar. Comparison of experimental data ${ }^{32,33}$ to calculation results of PC-SAFT ( $\mathrm{k}_{\mathrm{ij}}=0.016$ ) and SAFT ( $\mathrm{k}_{\mathrm{ij}}=-0.069$ ).

point of the liquid-liquid equilibrium at the given pressure). The vapor-liquid equilibrium is well-described using this $\mathrm{k}_{\mathrm{ij}}$ parameter. Figure 8 indicates that the association term of SAFT-type equations of state is a powerful expression that enables also the PC-SAFT model to correlate associating mixtures with only one binary interaction parameter covering wide temperature ranges.
4.3. Mixtures with Two Associating Substances. The modeling of mixtures in which associating interactions occur between different substances is particularly demanding for molecular theories. Also, the calculation results for cross-associating systems depend significantly on the combining rules (here, eqs 2 and 3 ).

A $\mathrm{T}-\mathrm{x}$ diagram of methanol and 1 -octanol at $\mathrm{P}=$ 1.013 bar is displayed in Figure 9. The two models perform equally well for this system.

The phase behavior of water and 1-pentanol at $\mathrm{P}=$ 1.013 bar is given in Figure 10 in a $\mathrm{T}-\mathrm{x}$ diagram. The system shows a liquid-liquid equilibrium at lower temperatures and a heteroazeotropic vapor-liquid equi-
librium. The $\mathrm{k}_{\mathrm{ij}}$ parameters of both the PC-SAFT equation of state and the SAFT model were fitted to the 1-pentanol-rich liquid concentration of the vapor-liquid-liquid equilibrium. The results of both models show some deviations from the vapor-liquid data. The PC-SAFT equation of state is found to agree well with the liquid-liquid data.

## 5. Conclusions

The perturbed-chain SAFT equation of state was used to model phase equilibria of associating components. The pure-component parameters of 18 associating substances were determined by simultaneously fitting vapor pressure and liquid density data. Good correlation results were obtained for alcohols, amines, and acetic acid, as a comparison to an earlier SAFT version proposed by Huang and Radosz confirms. Deviations of $6.83 \%$ in liquid density and $1.88 \%$ in vapor pressure were observed for water. Mixtures with one associating compound do not require mixing rules in the association term of the PC-SAFT equation of state. Both vaporliquid and liquid-liquid equilibria were calculated for mixtures of this type using one constant $\mathrm{k}_{\mathrm{ij}}$ parameter that corrects the dispersive interactions. Liquid-liquid and vapor liquid equilibria were also calculated for mixtures with two associating substances. Simple combination rules were adopted for the cross-associating interactions in this study; no additional binary interaction parameter was thereby introduced. A constant $\mathrm{k}_{\mathrm{ij}}$ parameter correcting the contribution due to dispersion was used in the correlation of vapor-liquid and liquidliquid demixing.

## Acknowledgment

The authors are grateful to the Deutsche Forschungsgemeinschaft for supporting this work with Grant SAD 700/3. We are indebted to Prof. Wolfgang Arlt for his continuous and energetic support of our work.

## Literature Cited

(1) Müller, E. A.; Gubbins, K. E. Molecular-Based Equations of State for Associating Fluids: A Review of SAFT and Related Approaches. Ind. Eng. Chem. Res. 2001, 40, 2193.
(2) Economou, I. G.; Donohue, M. D. Equations of state for hydrogen bonding systems. Fluid Phase Equilib. 1996, 116, 518.
(3) Wei, Y. S.; Sadus, R. J. Equations of State for the Calculation of Fluid-Phase Equilibria. AlChE J . 2000, 46, 169.
(4) Chapman, W. G.; J ackson, G.; Gubbins, K. E. Phase equilibria of associating fluids. Chain molecules with multiple bonding sites. Mol. Phys. 1988, 65, 1057.
(5) Chapman, W. G.; Gubbins, K. E.; J ackson, G.; Radosz, M. New Reference Equation of State for Associating Liquids. Ind. Eng. Chem. Res. 1990, 29, 1709.
(6) Wertheim, M. S. Fluids with highly directional attractive forces: I. Statistical thermodynamics. J. Stat. Phys. 1984, 35, 19.
(7) Wertheim, M. S. Fluids with highly directional attractive forces: II. Thermodynamic perturbation theory and integral equations. J. Stat. Phys. 1984, 35, 35.
(8) Wertheim, M. S. Fluids with highly directional attractive forces: III. Multiple attraction sites. J . Stat. Phys. 1986, 42, 459.
(9) Wertheim, M. S. Fluids with highly directional attractive forces: IV. Equilibrium polymerization. J . Stat. Phys. 1986, 42, 477.
(10) Chapman, W. G.; Gubbins, K. E.; J oslin, C. G.; Gray, C. G. Theory and simulation of associating liquid mixtures. Fluid Phase Equilib. 1986, 29, 337.
(11) J oslin, C. G.; Gray, C. G.; Chapman, W. G.; Gubbins, K. E. Theory and simulation of associating liquid mixtures. II. Mol. Phys. 1987, 62, 843.
(12) J ackson, G.; Chapman, W. G.; Gubbins, K. E. Phase equilibria of associating fluids. Spherical molecules with multiple bonding sites. Mol. Phys. 1988, 65, 1.
(13) Huang, S. H.; Radosz, M. Equation of State for Small, Large, Polydisperse, and Associating Molecules. Ind. Eng. Chem. Res. 1990, 29, 2284.
(14) Huang, S. H.; Radosz, M. Equation of State for Small, Large, Polydisperse, and Associating Molecules: Extensions to Fluid Mixtures. Ind. Eng. Chem. Res. 1991, 30, 1994.
(15) Gross, J .; Sadowski, G. Perturbed-Chain SAFT: An Equation of State Based on a Perturbation Theory for Chain Molecules. Ind. Eng. Chem. Res. 2001, 40, 1244.
(16) Gross, J .; Sadowski, G. Modeling Polymer Systems Using the Perturbed-Chain Statistical Associating Fluid Theory Equation of State. Ind. Eng. Chem. Res. 2002, 41, 1084-1093.
(17) These subroutines, written by J. Gross, G. Sadowski, and W. Arlt are available at http://thwww.chemietechnik. uni-dortmund.de.
(18) Pfohl, O.; Petkov, S.; Brunner, G. "PE" Quickly Makes Available the Newest Equations of State via the Internet. Ind. Eng. Chem. Res. 2000, 39, 4439.
(19) Wolbach, J. P.; Sandler, S. I. Using Molecular Orbital Calculations To Describe the Phase Behavior of Cross-Associating Mixtures. Ind. Eng. Chem. Res. 1998, 37, 2917.
(20) E conomou, I. G.; Tsonopoulos, C. Associating Models and Mixing Rules in Equations of State for Water/Hydrocarbon Mixtures. Chem. Eng. Sci. 1997, 52, 511.
(21) Suresh, S. J.; Elliott, J. R. Multiphase Equilibrium Analysis via a Generalized Equation of State for Associating Mixtures. Ind. Eng. Chem. Res. 1992, 31, 2783.
(22) Müller, E. A.; Gubbins, K. E. An Equation of State for Water from a Simplified Intermolecular Potential. Ind. Eng. Chem. Res. 1995, 34, 3662.
(23) Leu, A.-D.; Robinson, D. B. Equilibrium Phase Properties of the Methanol-Isobutane Binary System. J . Chem. Eng. Data 1992, 37, 10.
(24) Ellis, S. R. M.; Froome, B. A. Chem. Ind. 1954, 237. Data taken from: Gmehling, J .; Onken, U. Vapor-Liquid Equilibrium Data Collection; DECHEMA Chemistry Data Series; DECHEMA: Frankfurt, Germany, 1977; Vol. I, Part 2a.
(25) Rhodes, J. M.; Griffin, T. A.; Lazzaroni, M. J.; Bhethanabotla, V. R.; Campbell, S. W. Total pressure measurements for
benzene with 1-propanol, 2-propanol, 1-pentanol, 3-pentanol, and 2-methyl-2-butanol at 313.15 K. Fluid Phase Equilib. 2001, 179, 217.
(26) Deák, A.; Victorov, A. I.; de Loos, Th. W. High-pressure VLE in alkanol + alkane mixtures. Experimental results for n-butane + ethanol, + 1-propanol, + 1-butanol systems and calculations with three EOS methods. Fluid Phase Equilib. 1995, 107, 277.
(27) J ones, D. C.; Amstell, S. J . Chem. Soc. 1930, 1316. Data taken from: Sorensen, J . M.; Arlt, W. Liquid-Liquid Equilibrium Data Collection; DECHEMA Chemistry Data Series; DECHEMA: Frankfurt, Germany, 1979; Vol. V, Part 1.
(28) Madhavan, S.; Murti, P. S. Chem. Eng. Sci. 1966, 21, 465. Data taken from: Gmehling, J.; Onken, U. Vapor-Liquid Equilibrium Data Collection; DECHEMA Chemistry Data Series; DECHEMA: Frankfurt, Germany, 1977; Vol. I, Part 2a.
(29) Marinichev, A. N.; Susarev, M. P. Zh. Prikl. Khim. 1965, 38, 1619. Data taken from: Gmehling, J .; Onken, U. Vapor-Liquid Equilibrium Data Collection; DECHEMA Chemistry Data Series; DECHEMA: Frankfurt, Germany, 1977; Vol. I, Part 2a.
(30) Sorensen, J . M.; Arlt, W. Liquid-Liquid Equilibrium Data Collection; DECHEMA Chemistry Data Series; DECHEMA: Frankfurt, Germany, 1979; Vol. V, Part 1.
(31) Arce, A.; Blanco, A.; Soto, A.; Tojo, J . Isobaric VaporLiquid Equilibria of Methanol + 1-Octanol and Ethanol + 1-Octanol Mixtures. J . Chem. Eng. Data 1995, 40, 1011.
(32) Cho, T.-H.; Ochi, K.; Kojima, K. Kagaku Kogaku Ronbunshu 1984, 10, 181. Data taken from: Gmehling, J .; Onken, U.; Rarey-Nies, J. R. Vapor-Liquid Equilibrium Data Collection; DECHEMA Chemistry Data Series; DECHEMA: Frankfurt, Germany, 1988; Vol. I, Part 1b.
(33) Zhuravleva, I. K.; Zhuravlev, E. F. Izv. Vyssh. Ucheb. Zaved. Khim. Khim. Tekhnol. 1970, 13, 480. Data taken from: Sorensen, J . M.; Arlt, W. Liquid-Liquid Equilibrium Data Collection; DECHEMA Chemistry Data Series; DECHEMA: Frankfurt, Germany, 1979; Vol. V, Part 1.

Received for review November 28, 2001
Revised manuscript received April 22, 2002
Accepted April 25, 2002
IE010954D


[^0]:    * To whom correspondence should be addressed. Tel.: ++49-231-755 2635. Fax: ++49-231-755 2572. E-mail: G.Sadowski@ ct.uni-dortmund.de.
    ${ }^{\dagger}$ Technische Universität Berlin.
    ${ }^{\ddagger}$ Current address: BASF AG, Conceptual Process Design GIC/P-Q 920, 67056 Ludwigshafen, Germany. Fax: ++49-621-60 73488. E-mail: J oachim.Gross@alumni.tu-berlin.de.
    ${ }^{§}$ Universität Dortmund.

