# Measuring and modeling alcohol/salt systems 

Christoph Held, Axel Prinz, Viktoria Wallmeyer, Gabriele Sadowski*<br>Laboratory of Thermodynamics, Department of Biochemical and Chemical Engineering, Technische Universität Dortmund, Emil-Figge-Str. 70, 44227 Dortmund, Germany

## ARTICLE INFO

## Article history:

Received 12 May 2011
Received in revised form
12 September 2011
Accepted 20 September 2011
Available online 29 September 2011

## Keywords:

Phase equilibria
Equation of state
Ternary solutions
Electrolytes
Osmotic coefficient
Mean ionic activity coefficient


#### Abstract

Liquid densities, osmotic coefficients, and mean ionic activity coefficients (MIAC) at $25^{\circ} \mathrm{C}$ of single-salt alcohol (methanol and ethanol) solutions containing univalent ions were measured and modeled with the ePC-SAFT equation of state. In accordance with our previous work [Held, C., Cameretti, L.F., Sadowski, G., Fluid Phase Equlilib. 270 (2008) 87-96], only two solvent-specific ion parameters were adjusted to experimental solution densities and osmotic coefficients: the solvated ion diameter and the dispersion-energy parameter. ePC-SAFT was able to reproduce experimental data of the respective alcohol/salt systems with reasonable accuracy. Based on the solvent-specific ion-parameter sets, it is possible to predict densities and MIACs in ternary and quaternary water/alcohol(s)/salt solutions by introducing appropriate mixing rules that do not contain any additional fitting parameters.


© 2011 Elsevier Ltd. All rights reserved.

## 1. Introduction

Electrolyte solutions play an important role in biological and chemical engineering. These systems have been used in wastewater and drinking-water treatments, fertilizer production, and enhanced oil recovery (Enick and Klara, 1992). Water/salt systems have been widely investigated over the last century. The increasing interest in biological solutions requires a better understanding of how salt influences the thermodynamic behavior of aqueous solutions, which we have addressed in the past (Held et al., 2008; Held and Sadowski, 2009). Because biological and pharmaceutical systems often contain cosolvents, more complex solutions of electrolytes in mixed solvents (here water and alcohols) must also be considered. Despite the availability of well-established experimental methods, the data in the literature is rather sparse considering that there is an almost endless number of possible system conditions, for example, various types of solvent(s) and salt(s) and various temperatures and/or pressures, solvent compositions, and salt concentrations.

To complicate matters, the variety of experimental methods for determining thermodynamic properties (e.g., activity coefficients) in such systems is exceeded by the number of thermodynamic models developed for calculating these properties (Cameretti et al., 2005; Fürst and Renon, 1993; Huang et al., 2009; Ji and Adidharma, 2006; Ji et al., 2005; Liu et al., 1999; Papaiconomou et al., 2002; Salimi et al., 2005; Simonin et al.,

[^0]2006; Wang et al., 2011; Ye et al., 1994; Zuo et al., 2000). Most models were developed to describe very specific systems by fitting system-specific model parameters instead of being generally valid or predictive. This may be due to the absence or inconsistency of activity-coefficient data in binary alcohol/salt systems. For ethanol, for example, only the vapor-pressure depressions in $\mathrm{NaI}, \mathrm{LiCl}$, and LiBr solutions are available (Mato and Cocero, 1988). For methanol, the data in the literature for ammonium-salt solutions is highly unreliable. Therefore, it is impossible to determine universal ionic parameter sets for establishing predictive models based on the current binary solvent/salt data. This was the objective of this work. Here, we present an appropriate model that is able to calculate the properties of ternary (two solvents + salt) systems based on ion parameters that were fitted to binary (one solvent + salt) solution data using universal ionic parameters.

Due to the lack of or inconsistency in the literature data, new experimental activity-coefficient data were determined within this work for alcohol/salt systems. Solvent activity coefficients in binary solutions can be measured using a vapor-pressure osmometer, and such measurements have been performed for aqueous amino-acid systems (Held et al., 2011). These data can be converted into mean ionic activity coefficients (MIAC) by the Gibbs-Duhem relationship. Analogous to our previous work (Held et al., 2008), model parameters were adjusted to both these MIACs and solution-density data. As data from the literature were not always consistent, we also determined new liquid-density data for ethanol/salt systems. Thus, we present new experimental solution densities, osmotic coefficients, and MIAC data for eight salts in methanol and ethanol.

The first step toward modeling complex multicomponent solutions is accurately describing binary solvent/salt systems. For that purpose, the interactions between all species must be covered using a physically sound approach. In addition to the nonionic short-range interactions resulting from repulsive and attractive (van der Waals) forces, the long-range Coulombic interactions of the charged species must also be accounted for. The first forces can be described by using common $g^{\mathrm{E}}$ models or equations of state (EOS). The Coulombic interactions can be modeled by the Debye-Hückel theory (DH) developed in 1923 (Debye and Hückel, 1923) or by the mean spherical approximation (MSA), which solves the Ornstein-Zernike equation for a fluid of charged spheres of equal size (Waisman and Lebowitz, 1970). Coulomb and non-Coulomb expressions can be combined to yield state-of-the-art electrolyte models, which have been widely applied to binary salt/water systems (for references see (Held et al., 2008). To describe ternary systems (e.g., a salt dissolved in two solvents), such approaches can be extended to powerful multicomponent models by introducing appropriate mixing rules. Electrolyte $g^{\mathrm{E}}$ models that have been applied to solutions with more than one solvent or one solute include the modified Pitzer approach by Ye et al. (1994) and the model by Papaiconomou et al. (2002), which combined the NRTL theory with the MSA to calculate MIACs in solvent mixtures. However, the two groups applied ternary mixture parameters, which they fitted to the properties of their respective systems. Salimi et al. (2005) defined concentration-dependent "ion parameters" that vary based on the salt. For example, $\mathrm{Cl}^{-}$has different parameters in NaCl than in KCl . Applying this approach, MIACs in mixed-solvent and mixedsalt electrolyte solutions could be predicted without fitting ternary parameters.

System densities can also be calculated using EOS instead of $g^{\mathrm{E}}$ models. In the model of Fürst and Renon (1993), the Redlich-Kwong Soave EOS for short-range interactions was combined with the MSA term for long-range interactions. A second short-range contribution for ion-solvent forces was also included to predict osmotic coefficients in multi-solute aqueous solutions. Later, Zuo et al. (2000) extended this approach to mixed-solvent systems. However, as in the work of Salimi et al. (2005), the "ion parameters" had to be changed for different salts. Liu et al. (1999) could predict MIACs in two-salt solutions; however, they used salt-specific parameters. SAFT approaches were also applied to electrolyte solutions. For example, Ji et al. (2005) proposed the SAFT1-RPM to model osmotic coefficients in multi-solute aqueous systems with ion-specific parameters and subsequently extended the model (SAFT2) to also describe seawater solutions (Ji and Adidharma, 2006); however, they had to fit additional ion parameters.

In this work, we apply the ePC-SAFT model developed by Cameretti et al. (2005) to aqueous electrolyte solutions containing methanol and ethanol as cosolvents. The ePC-SAFT model is a combination of the PC-SAFT EOS by Gross and Sadowski (2001) and the Debye-Hückel contribution (Debye and Hückel, 1923), which accounts for the Coulomb interactions. Only two solventspecific parameters are used to characterize each ion: the solvated ion diameter $\sigma_{j}$ and the dispersive energy parameter $u_{j} / k_{B}$, which both reflect the strength of ionic solvation. Within ePCSAFT, the ionic species are considered independent of the salt they are part of. We present ion-parameter sets that are able to reasonably describe solution densities, MIACs, and vapor-liquid equilibria (VLE) of binary alcohol/electrolyte systems. This work differs from the above-mentioned studies in that these binary ion/solvent parameters are directly used to model ternary systems by applying standard mixing rules without introducing new parameters. This allows one to predict the phase behavior in electrolyte two-solvent systems without impairing the calculations of the binary solvent/salt solutions.

## 2. Experimental work

Although the literature contains the thermodynamic properties of numerous aqueous electrolyte solutions, experimental data for ethanol/salt systems is rather scarce. Therefore, in this work, we measured the solution densities and osmotic coefficients of several ethanol/electrolyte solutions. To determine solution densities, a vibrating-tube densimeter was used. This method has previously been applied to electrolyte solutions (Zhang and Han, 1996). However, determining osmotic/activity coefficients is much more complicated. The measurements were carried out using a vapor-pressure osmometer, which has proven its applicability in our earlier works (Held et al., 2011, 2010) in which we precisely determined osmotic coefficients in aqueous systems. Moreover, Mato and Cocero (1988) have shown that vaporpressure osmometry is also appropriate for measuring activity coefficients in ethanol/electrolyte solutions.

In the case of most methanol/electrolyte systems, experimental solution data is available in the literature. However, osmotic coefficients (1) appear to vary significantly, especially for the ammonium solutions and (2) are not available for lithium iodide in methanol. Thus, the vapor-pressure osmometer was also used to determine the osmotic coefficients in these systems.

### 2.1. Materials

LiChrosolv® ethanol from Merck with a purity of $99.9 \%$ and methanol (Merck, $>99.9 \%$ ) were used to prepare the electrolyte solutions. Lithium iodide (Merck, >98\%), lithium bromide (Sigma, $>99 \%$ ), lithium chloride (VWR Prolabo $>99 \%$ ), sodium iodide (Fluka $>99 \%$ ), sodium bromide (Merck, > 99\%), sodium chloride (Merck, $>=99.5 \%$ ), ammonium iodide (Acros, $>=99 \%$ ), ammonium bromide (Merck, $>=99.5 \%$ ), and ammonium chloride (Merck, $>99 \%$ ) were used without further purification. All solutions were prepared gravimetrically by weighing with an uncertainty of 0.01 mg .

### 2.2. Solution densities

Whereas the solution densities of methanol/salt systems are in the literature, those of ethanol/salt systems are either not available or inconsistent. Thus, densities of binary ethanol/salt solutions containing LiCl, LiBr, LiI, NaBr, NaI, NH4 Cl, NH4 Br, and NH4 I were determined using a vibrating-tube densimeter "DMA 602" from Anton Paar Germany GmbH (Ostfildern, Germany) at $25^{\circ} \mathrm{C}$ and ambient pressure. For these measurements, a U-tube is filled with the fluid of interest and set into oscillation by an electromagnetic field. Densities are obtained by measuring the Eigen frequency of the filled U-tube. The apparatus was calibrated with air and deionized water, and reference data were taken from (Wagner and Pruss, 2002). According to the manufacturer, the maximum uncertainty of this apparatus is within $\pm 1.5 \times 10^{-6} \mathrm{g} / \mathrm{cm}^{3}$. The results are summarized in Table 1. The density of pure ethanol at $25^{\circ} \mathrm{C}$ was determined prior to measuring each salt solution; its average density ( $\rho_{\text {ethanol }}=785.54 \mathrm{~kg} / \mathrm{m}^{3}$ ) is consistent with the published density of pure ethanol $\left(785.10 \mathrm{~kg} / \mathrm{m}^{3}\right.$; Lacmann and Synowietz, 1977). Based on our experience, knowing the initial slope of $\rho$ vs. $m$ enables one to extrapolate the mixture densities to higher salt concentrations. Thus, we particularly considered the regions of low concentration for measuring the solution densities of the ethanol/salt systems shown in Table 1.

### 2.3. Osmotic coefficients

Osmotic coefficients were measured with two different vaporpressure osmometers: the Gonotec Osmomat O70 for small

Table 1
Experimental densities of $\mathrm{LiCl}, \mathrm{LiBr}, \mathrm{LiI}, \mathrm{NaBr}, \mathrm{NaI}, \mathrm{NH}_{4} \mathrm{Cl}, \mathrm{NH}_{4} \mathrm{Br}$, and $\mathrm{NH}_{4} \mathrm{I}$ in ethanol at $25^{\circ} \mathrm{C}$ and ambient pressure.
| LiCl |  | LiBr |  | LiI |  | NaBr |  |
| :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- |
| $m_{\mathrm{LiCl}}$ [mol/kg] | $\rho_{\text {solution }}$ [ $\mathrm{kg} / \mathrm{m}^{3}$ ] | $m_{\mathrm{LiBr}}$ [mol/kg] | $\rho_{\text {solution }}$ [ $\mathrm{kg} / \mathrm{m}^{3}$ ] | $m_{\text {LiI }}$ [mol/kg] | $\rho_{\text {solution }}$ [ $\mathrm{kg} / \mathrm{m}^{3}$ ] | $m_{\mathrm{NaBr}}$ [mol/kg] | $\rho_{\text {solution }}$ [ $\mathrm{kg} / \mathrm{m}^{3}$ ] |
| 0.6044 | 805.52 | 0.2943 | 805.04 | 0.193 | 804.12 | 0.0884 | 792.23 |
| 1.2346 | 823.73 | 0.6049 | 823.87 | 0.3939 | 823.16 | 0.1422 | 795.96 |
|  |  |  |  |  |  | 0.2192 | 800.77 |
| NaI |  | $\mathrm{NH}_{4} \mathrm{Cl}$ |  | $\mathrm{NH}_{4} \mathrm{Br}$ |  | $\mathrm{NH}_{4} \mathrm{I}$ |  |
| $m_{\text {NaI }}$ [mol/kg] | $\rho_{\text {solution }}$ [ $\mathrm{kg} / \mathrm{m}^{3}$ ] | $m_{\mathrm{NH} 4 \mathrm{Cl}}$ [mol/kg] | $\rho_{\text {solution }}$ [ $\mathrm{kg} / \mathrm{m}^{3}$ ] | $m_{\mathrm{NH} 4 \mathrm{Br}}$ [mol/kg] | $\rho_{\text {solution }}$ [ $\mathrm{kg} / \mathrm{m}^{3}$ ] | $m_{\mathrm{NH} 4 \mathrm{I}}$ [mol/kg] | $\rho_{\text {solution }}$ [ $\mathrm{kg} / \mathrm{m}^{3}$ ] |
| 0.1715 | 803.38 | 0.0376 | 786.45 | 0.1046 | 791.64 | 0.1769 | 800.90 |
| 0.3511 | 820.99 | 0.0907 | 788.09 | 0.2064 | 797.58 | 0.3658 | 817.06 |
|  |  | 0.1377 | 789.06 | 0.2975 | 802.71 |  |  |


concentrations and the Knauer K-7000 for more concentrated solutions. The concentration limits of the osmometers were approximately 0.6 mol solute $/ \mathrm{kg}$ solvent $(0.6 \mathrm{~m})$ for the Gonotec and 2.5 m for the Knauer instruments, depending on the salt and solvent. Measurements that did not agree within $2 \%$ (related to osmotic coefficients) between both osmometers (e.g., at 0.4 m ) were discarded.

The measuring cell of a vapor-pressure osmometer contains two thermistors in a tempered, solvent-saturated atmosphere. Using a syringe, a drop of the solution of interest is placed on one thermistor while a drop of the pure solvent is placed on the second thermistor as a reference. The solute causes a vapor-pressure depression of the solvent. That is, the solvent condensates from the atmosphere into the solution drop. The resulting condensation enthalpy causes a temperature increase within the solution drop. The temperature difference between the solution and pure solvent is measured using a Wheatstone bridge; this difference is converted into osmotic coefficients (see (Held et al., 2010)).

Prior to the measurements, the instrument was calibrated with a system of known osmotic coefficients. Among the few systems with appropriate and available experimental data, alcoholic solutions (ethanol or methanol) containing sodium iodide are the most investigated. Thus, we used ethanol/NaI and methanol/NaI to calibrate our osmometers.

Barthel and Lauermann (1986) published the following relation for the concentration-dependent vapor-pressure depression of an ethanol/sodium iodide system:

$$
\begin{align*}
\Delta p^{\mathrm{LV}}\left(25^{\circ} \mathrm{C}\right)= & 0.008873+4.04961 m_{\mathrm{NaI}}-1.07993 m_{\mathrm{NaI}}^{2} \\
& +1.85744 m_{\mathrm{NaI}}^{3}-0.461216 m_{\mathrm{NaI}}^{4} \tag{1}
\end{align*}
$$

The $\Delta p^{\mathrm{LV}}$ values can be converted into osmotic coefficients ( $\phi$ ) or solvent activity coefficients ( $\gamma_{\text {solv }}$ ) by

$$
\begin{equation*}
a_{\mathrm{solv}}=\frac{p_{0 \mathrm{solv}}^{\mathrm{LV}}-\Delta p^{\mathrm{LV}}}{p_{0 \mathrm{solv}}^{\mathrm{LV}}} ; \quad \gamma_{\mathrm{solv}}=\frac{a_{\mathrm{solv}}}{x_{\mathrm{solv}}} ; \quad \phi=-\frac{\ln a_{\mathrm{solv}}}{v_{ \pm} m_{ \pm} M_{\mathrm{solv}}} \tag{2}
\end{equation*}
$$

where $a_{\text {solv }}$ is the activity of the solvent at the salt molality $m_{ \pm}$of the considered $1: 1$ electrolyte ( $v_{ \pm}=2$ ).

For methanol/sodium iodide solutions, we used data published by Tomasula et al. (1985) for the calibration; their experimental osmotic coefficients are best correlated by

$$
\begin{align*}
\phi\left(m_{\mathrm{Nal}}<0.2 \mathrm{~mol} / \mathrm{kg}\right)= & 417,726 m_{\mathrm{Nal}}^{6}-252,772 m_{\mathrm{Nal}}^{5}+59,746 m_{\mathrm{Nal}}^{4} \\
& -6,968 m_{\mathrm{Nal}}^{3}+418.02 m_{\mathrm{Nal}}^{2}-12.403 m_{\mathrm{Nal}}+1 \\
\phi\left(m_{\mathrm{Nal}} \geq 0.2 \mathrm{~mol} / \mathrm{kg}\right)= & 0.0002 m_{\mathrm{Nal}}^{6}-0.006 m_{\mathrm{Nal}}^{5}+0.0549 m_{\mathrm{Nal}}^{4} \\
- & 0.2182 m_{\mathrm{Nal}}^{3}+0.4128 m_{\mathrm{Nal}}^{2}-0.0826 m_{\mathrm{Nal}}+0.8427 \tag{3}
\end{align*}
$$

Calibrating the osmometer makes it possible to ascribe the signal of each measurement (i.e., the detected temperature difference between the two droplets) to a certain osmolality. Experimental osmotic coefficients and the osmometer signals ( $\Delta T$ ) are related via

$$
\begin{equation*}
\Delta T=K_{\mathrm{eb}} \cdot \phi \cdot v_{ \pm} m_{ \pm} \tag{4}
\end{equation*}
$$

where $K_{\mathrm{eb}}$ is the ebullioscopic constant of the respective solvent.
For the calibration and measurements, each experiment was repeated five times; these values were recorded and averaged. The average relative standard deviation of the experiments was approximately $2 \%$. The resulting osmotic coefficients are given in Tables 2 and 3 for the methanol and ethanol solutions, respectively.

In addition to the activity coefficient of the alcohols, the rational mean ionic activity coefficients (MIAC) $\gamma_{ \pm \pm}^{*}$ of the salts are also of interest in many applications. Applying the GibbsDuhem relation allows for the conversion of osmotic-coefficient data into MIACs. For that purpose, the experimental osmotic coefficients were first correlated by a power series

$$
\begin{equation*}
\phi-1=\sum_{i=1}^{n} A_{i} m^{i} \tag{5}
\end{equation*}
$$

and then converted into MIACs by applying the following equation:

$$
\begin{equation*}
\ln \gamma_{ \pm}^{*}=(\phi-1)+\int_{0}^{m} \frac{(\phi-1)}{m} d m \tag{6}
\end{equation*}
$$

In Eq. (5), the $A_{i}$ values are the coefficients of the power series, and $n$ refers to the number of parameters needed to represent the experimental osmotic coefficients. Note that we did not directly adjust the $A$ coefficients to the raw osmotic coefficients but rather to a best-fit correlation for our $\phi, m$-data. For example, for $\mathrm{NH}_{4} \mathrm{Br}$ at $25^{\circ} \mathrm{C}, A_{1}$ to $A_{5}$ were found to be $-15,216.8,-1289,3344$, and -3139, respectively. The MIACs for the measured systems can also be found in Tables 2 and 3.

As a typical example, the results of the osmometer measurements for ethanol/ $\mathrm{NH}_{4} \mathrm{Br}$ solutions are illustrated in Fig. 1. As for the aqueous salt systems, the MIACs are much smaller than the osmotic coefficients; both values reach unity at zero salt concentration.

## 3. Modeling with ePC-SAFT

## 3.1. ePC-SAFT equation of state

The ePC-SAFT EOS is based on a perturbation theory in which the hard-chain system is used as the reference system. All other

Table 2
Experimental osmotic coefficients of ethanol/salt solutions at $25^{\circ} \mathrm{C}$. MIACs were obtained by Eq. (6).
| $\mathbf{m}_{\text {LiCI }}$ [mol/kg] | $\phi$ [-] | $\gamma_{\text {LiCI }}^{*}$ [-] | $\boldsymbol{m}_{\text {LiBr }}$ [mol/kg] | $\phi$ [-] | $\gamma_{\text {LiBr }}^{*}$ [-] | $\boldsymbol{m}_{\text {LiI }}$ [mol/kg] | $\phi$ [-] | $\gamma_{\text {LiI }}^{*}$ [-] |
| :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- |
| 0.250 | 0.747 | 0.218 | 0.250 | 0.735 | 0.212 | 0.400 | 0.731 | 0.195 |
| 0.400 | 0.772 | 0.206 | 0.400 | 0.710 | 0.185 | 0.600 | 0.819 | 0.223 |
| 0.750 | 0.800 | 0.201 | 0.750 | 0.791 | 0.200 | 1.000 | 0.996 | 0.246 |
| 1.000 | 0.767 | 0.192 | 1.000 | 0.840 | 0.177 | 1.250 | 1.034 | 0.257 |
|  |  |  | 1.250 | 0.926 | 0.159 | 1.500 | 1.034 | 0.270 |
|  |  |  | 1.500 | 1.069 | 0.209 |  |  |  |
| $\mathbf{m}_{\mathrm{NH} 4 \mathrm{Cl}}$ [mol/kg] | $\phi$ [-] | $\gamma *_{\text {NH4Cl }}$ [-] | $\mathbf{m}_{\text {NH4Br }}$ [mol/kg] | $\phi$ [-] | $\gamma_{\mathrm{NH} 4 \mathrm{Br}}^{*}$ [-] | $\mathrm{m}_{\text {NH4I }}$ [mol/kg] | $\phi$ [-] | $\gamma_{\mathrm{NH} 4 \mathrm{I}}^{*}$ [-] |
|  |  |  |  |  |  |  |  |  |
| 0.038 | 0.860 | 0.471 | 0.105 | 0.728 | 0.337 | 0.200 | 0.593 | 0.251 |
| 0.091 | 0.697 | 0.419 | 0.206 | 0.711 | 0.277 | 0.400 | 0.627 | 0.201 |
|  |  |  | 0.297 | 0.717 | 0.239 | 0.600 | 0.670 | 0.176 |
|  |  |  |  |  |  | 1.000 | 0.638 | 0.097 |
| $\mathbf{m}_{\text {NaCl }}$ [mol/kg] |  | $\boldsymbol{\gamma}_{\mathbf{N a C l}}^{*}$ [-] |  | $\phi$ [-] | $\gamma_{\text {NaBr }}^{*}$ [-] |  |  |  |
|  | $\begin{aligned} & \phi \\ & {[-]} \end{aligned}$ |  | $\mathbf{m}_{\text {NaBr }}$ [mol/kg] |  |  |  |  |  |
| 0.009 | 0.924 | 0.698 | 0.088 | 0.790 | 0.371 |  |  |  |
|  |  |  | 0.142 | 0.779 | 0.332 |  |  |  |
|  |  |  | 0.228 | 0.734 | 0.289 |  |  |  |


Table 3
Experimental osmotic coefficients of methanol/salt solutions at $25^{\circ} \mathrm{C}$. MIACs were obtained by Eq. (6).
| $m_{\text {NH4Br }}$ [mol/kg] | $\phi$ [-] | $\gamma_{\mathrm{NH} 4 \mathrm{Br}}^{*}$ [-] | $\boldsymbol{m}_{\text {NH4I }}$ [mol/kg] | $\phi$ [-] | $\gamma_{\mathbf{N} \mathbf{H} \mathbf{4 I}}^{*}$ [-] | $\boldsymbol{m}_{\text {LiI }}$ [mol/kg] | $\phi$ [-] | $\gamma_{\text {LiI }}^{*}$ [-] |
| :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- |
| 0.102 | 0.796 | 0.362 | 0.103 | 0.864 | 0.361 | 0.102 | 0.945 | 0.409 |
| 0.201 | 0.791 | 0.355 | 0.203 | 0.779 | 0.356 | 0.202 | 0.891 | 0.410 |
| 0.297 | 0.773 | 0.377 | 0.300 | 0.808 | 0.377 | 0.298 | 0.908 | 0.420 |
| 0.399 | 0.774 | 0.328 | 0.400 | 0.800 | 0.331 | 0.403 | 0.889 | 0.410 |
| 0.606 | 0.806 | 0.224 | 0.596 | 0.851 | 0.238 | 0.598 | 1.037 | 0.414 |
| 0.791 | 0.840 | 0.270 | 0.799 | 0.889 | 0.326 |  |  |  |


![](https://cdn.mathpix.com/cropped/137993f5-9955-4aaa-862d-a06927e0cd5f-04.jpg?height=467&width=860&top_left_y=1545&top_left_x=155)
Fig. 1. Experimental osmotic coefficients (circles) in $\mathrm{NH}_{4} \mathrm{Br} /$ ethanol solutions at $25^{\circ} \mathrm{C}$ as function of salt molality and MIACs (squares) of $\mathrm{NH}_{4} \mathrm{Br}$ calculated with Gibbs-Duhem (Eq. (5)). Lines drawn to guide the readers' eyes.

contributions are assumed to be additive and can be considered independently. Thus, the residual Helmholtz energy, $a^{\text {res }}$, can be written as follows:

$$
\begin{equation*}
a^{\mathrm{res}}=a^{\mathrm{hc}}+a^{\mathrm{disp}}+a^{\mathrm{assOC}}+a^{\mathrm{ion}} \tag{7}
\end{equation*}
$$

where $N$ is the total number of molecules; $a^{\text {hc }}$ represents the hardchain repulsion of the reference system; and $a^{\text {disp }}, a^{\text {assoc }}$, and $a^{\text {ion }}$ account for the Helmholtz-energy contributions due to dispersive, associative, and Coulomb interactions, respectively. Whereas expressions for $a^{\text {disp }}$ and $a^{\text {assoc }}$ are used as in the original PC-SAFT model (Gross and Sadowski, 2001), Cameretti et al. (2005) used a Debye-Hückel term to describe the Helmholtz-energy contribution
$a^{\text {ion }}$ to a system that is caused by charging the species $j$

$$
\begin{equation*}
a^{\text {ion }}=-\frac{\kappa}{12 \pi \varepsilon} \sum_{j} x_{j} q_{j}^{2} \chi_{j} \tag{8}
\end{equation*}
$$

Here $x_{j}$ and $q_{j}$ are the mole fraction and charge of ion $j$, respectively. The ions are treated as spherical species in a uniform dielectric continuum characterized by a dielectric constant $\varepsilon$. They can approach each other to the distance $a_{j}$, which is equivalent to the ion diameter $\sigma_{j}$. The quantity $\chi_{j}$ in Eq. (8) is defined as

$$
\begin{equation*}
\chi_{j}=\frac{3}{\left(\kappa a_{j}\right)^{3}}\left[\frac{3}{2}+\ln \left(1+\kappa a_{j}\right)-2\left(1+\kappa a_{j}\right)+\frac{1}{2}\left(1+\kappa a_{j}\right)^{2}\right] \tag{9}
\end{equation*}
$$

with $\kappa$ being the inverse Debye screening length given by

$$
\begin{equation*}
\kappa=\sqrt{\frac{N_{A}}{k_{B} T \varepsilon} \sum_{j} q_{j}^{2} c_{j}}=\sqrt{\frac{\rho_{N} e^{2}}{k_{B} T \varepsilon} \sum_{j} z_{j}^{2} x_{j}} \tag{10}
\end{equation*}
$$

in which $c_{j}$ is the molar concentration and $N_{A}$ is Avogadro's constant.
In contrast to other groups (Ball et al., 1985, Triolo et al., 1976), we used a dielectric constant, $\varepsilon$, for solvents that neglects the influence of added salt. Table 4 contains empirical equations for the temperature dependence of $\varepsilon$ in pure methanol ("0M") and pure ethanol ("0E") as well as the concentration dependence of $\varepsilon$ in salt-free solvent mixtures at $25^{\circ} \mathrm{C}$.

In ePC-SAFT, two parameters describe the ion properties: the solvated ion diameter, $\sigma_{j}$, and the dispersion-energy parameter, $u_{j} / k_{B}$. For the associating nonspherical solvents in this work, three additional pure-component parameters are required, namely the segment number, $m_{\text {seg }}$, the association-energy parameter, $\varepsilon_{\mathrm{hb}}^{A i B i} / k_{B}$, and the association volume, $\kappa_{\mathrm{hb}}^{A i B i}$. The association scheme used here for all solvents is the two-site 2B approach (Huang and Radosz, 1990).

To describe the mixtures, the conventional Berthelot-Lorenz-combining rules were used

$$
\begin{equation*}
\sigma_{i j}=\frac{1}{2}\left(\sigma_{i}+\sigma_{j}\right)\left(1-l_{i j}\right) \tag{11}
\end{equation*}
$$

$$
\begin{equation*}
u_{i j}=\sqrt{u_{i} u_{j}}\left(1-k_{i j}\right) \tag{12}
\end{equation*}
$$

Eq. (12) is applied for solvent/solvent and solvent/ion interactions only. Dispersion between two ions is neglected in this work. Often, a binary interaction parameter ( $k_{i j}$ ) is introduced for the

Table 4
Relative permittivity of solvents and their mixtures, fitted to experimental data (Wohlfarth, 1991) with $T$ given in Kelvin.
| Methanol (M) | $\varepsilon_{r, 0 \mathrm{M}}(T)=-53.398 \ln (T)+336.170$ |
| :--- | :--- |
| Ethanol (E) | $\varepsilon_{r, 0 \mathrm{E}}(T)=-0.132 T+64.072$ |
| Methanol/water $\left(\mathrm{H}_{2} \mathrm{O}\right)$ | $\varepsilon_{r, \mathrm{M} / \mathrm{H}_{2} \mathrm{O}}\left(25^{\circ} \mathrm{C}\right)=9.3986 x_{\mathrm{H}_{2} \mathrm{O}}^{3}+13.622 x_{\mathrm{H}_{2} \mathrm{O}}^{2}+25.014 x_{\mathrm{H}_{2} \mathrm{O}}+\varepsilon_{r, \mathrm{OM}}\left(25^{\circ} \mathrm{C}\right)$ |
| Ethanol/water $\left(\mathrm{H}_{2} \mathrm{O}\right)$ | $\varepsilon_{r, \mathrm{E} / \mathrm{H}_{2} \mathrm{O}}\left(25^{\circ} \mathrm{C}\right)=50.738 x_{\mathrm{H}_{2} \mathrm{O}}^{3}-17.327 x_{\mathrm{H}_{2} \mathrm{O}}^{2}+21.874 x_{\mathrm{H}_{2} \mathrm{O}}+\varepsilon_{r, \mathrm{OE}}\left(25^{\circ} \mathrm{C}\right)$ |


Table 5
Pure-component PC-SAFT parameters for the solvents considered in this work and binary interaction parameters between the alcohols and water.
| Parameter | Abbreviation | Water ${ }^{\mathrm{a}}$ | Methanol ${ }^{\text {b }}$ | Ethanol ${ }^{\mathrm{b}}$ |
| :--- | :--- | :--- | :--- | :--- |
| Segment number | $m_{\text {seg }}$ | 1.2047 | 1.5255 | 2.3827 |
| Segment diameter | $\sigma$ | 2.7927 | 3.2300 | 3.1771 |
| Dispersion energy | $u / k_{B}$ | 353.94 | 188.90 | 198.24 |
| Association sites | $N$ | 2 | 2 | 2 |
| Association energy | $\varepsilon_{\mathrm{hb}}^{A i B i} / k_{B}$ | 2425.7 | 2899.5 | 2653.4 |
| Association volume | $\kappa_{\mathrm{hb}}^{\mathrm{AiBi}}$ | 0.045099 | 0.035176 | 0.032384 |
| Binary interaction parameter with water | $k_{i j}\left(\mathrm{H}_{2} \mathrm{O}\right)$ | - | $-0.085{ }^{\text {c }}$ | $-0.049{ }^{\text {c }}$ |
| Binary interaction parameter with water | $k_{i j}^{\mathrm{hb}}\left(\mathrm{H}_{2} \mathrm{O}\right)$ | - | - | $0.20^{\text {c }}$ |
| Binary interaction parameter with water | $l_{i j}\left(\mathrm{H}_{2} \mathrm{O}\right)$ | - | - | $-0.01{ }^{\text {c }}$ |


${ }^{\mathrm{a}}$ Cameretti and Sadowski (2008).
${ }^{\mathrm{b}}$ Gross and Sadowski (2002).
${ }^{\mathrm{c}}$ Binary interaction parameters determined in this work by data from (Hall et al., 1979).
dispersive-energy parameter to account for deviations from the geometric-mean rule in the mixture (Eq. (12)). This parameter is usually fitted to the respective binary mixture's properties. The cross-dispersion interaction between solvents is corrected with respective $k_{i j}$ values (given in Table 5), which were fitted to VLE data of the respective salt-free systems. All $k_{i j}$ parameters between solvents and ions were set to zero in this work.

To account for cross association between the two solvents (e.g., water and ethanol), Eqs. (13) and (14) proposed by Wolbach and Sandler (1998) were applied

$$
\begin{equation*}
\varepsilon_{\mathrm{hb}}^{A_{i} B_{j}}=\frac{\varepsilon_{\mathrm{hb}}^{A_{i} B_{j}}+\varepsilon_{\mathrm{hb}}^{A_{i} B_{j}}}{2} \tag{13}
\end{equation*}
$$

$$
\begin{equation*}
\kappa^{A i B j}=\sqrt{\kappa^{A i B i} \kappa^{A j B j}}\left(\frac{\sqrt{\sigma_{i i} \sigma_{j j}}}{1 / 2\left(\sigma_{i i}+\sigma_{j j}\right)}\right)^{3}\left(1-k_{i j}^{\mathrm{hb}}\right) \tag{14}
\end{equation*}
$$

The parameters $l_{i j}$ and $k_{i j}^{A i B j}$ in Eqs. (11) and (14) were only necessary for the mixture of ethanol and water and were fitted to the binary salt-free VLE at $25^{\circ} \mathrm{C}$ (data taken from Hall et al., 1979). Using this approach, the PC-SAFT equation of state allows for a quantitative description of the VLE in alcohol/water mixtures. All model parameters are given in Table 5. Note that without $l_{i j}$ (which is a rather unusual parameter), the VLE in water/ethanol at $25^{\circ} \mathrm{C}$ cannot be described quantitatively using the previously published PC-SAFT parameters (see e.g. Fuchs et al., 2006). This may be due to the high asymmetry of the water/ethanol system. For higher temperatures or other water/ alcohol systems, the parameter $l_{i j}$ is not necessary.

### 3.2. Modeling activity and osmotic coefficients

Calculating activity coefficients with an equation of state requires determining fugacity coefficients. The exact relationship between fugacity coefficients and the residual Helmholtz energy (Eq. (7)) is given elsewhere (Held et al., 2010). Solvent activity coefficients, $\gamma_{\text {solv }}$, are normalized to the pure solvent "Osolv" ( $x_{\text {solv }}=1$ )

$$
\begin{equation*}
\gamma_{\text {solv }}=\frac{\varphi_{\text {solv }}\left(T, p, x_{\text {solv }}\right)}{\varphi_{\text {Osolv }}\left(T, p, x_{\text {solv }}=1\right)} \tag{15}
\end{equation*}
$$

In contrast to solvent activity coefficients, $\gamma_{\text {solv }}$, solute activity coefficients, $\gamma_{j}{ }^{*}$, are related to infinite dilution ( $x_{j}=0$ ). The MIAC $\gamma_{ \pm}^{*, x}$ of an electrolyte is defined as the geometric mean of the mole-fraction-based rational activity coefficients of the ions in a solution (Robinson and Stokes, 1970)

$$
\begin{equation*}
\gamma_{ \pm}^{*, x}=\left(\left(\gamma_{+}^{*, x}\right)^{v_{+}} \cdot\left(\gamma_{+}^{*, x}\right)^{v_{-}}\right)^{1 /\left(v_{+}+v_{-}\right)} \tag{16}
\end{equation*}
$$

Here $v_{+}$and $v_{-}$are the stoichiometric coefficients of cations and anions in the salt (Tan et al., 2005), which add to $v$. In this work, the salts are assumed to be fully dissociated, i.e., the stoichiometric factor $v=2$ holds for all $1: 1$ electrolytes. The rational activity coefficients $\gamma_{j}^{*, x}$ of the ions are obtained by ePC-SAFT as the ion fugacity coefficient $\varphi_{j}$ at the actual concentration related to that at the infinite diluted state $\varphi_{j}^{\infty}$ at the same temperature and pressure

$$
\begin{equation*}
\gamma_{j}^{*, x}=\frac{\varphi_{j}\left(T, p, x_{j}\right)}{\varphi_{j}^{\infty}\left(T, p, x_{ \pm} \rightarrow 0\right)} \tag{17}
\end{equation*}
$$

The reference state "infinite dilution" in a binary salt/solvent solution is thus calculated for a system without salt ( $x_{ \pm}=0$ ), containing only pure solvent ( $x_{\text {solv }}=1$ ). In the case of two solvents, "infinite dilution" means the salt-free system ( $x_{ \pm}=0$ ) at the same solvent composition (e.g., $x_{\text {solv } 1} / x_{\text {solv } 2}=C$ ) as the considered solution

$$
\begin{equation*}
\gamma_{j}^{*, x, \text { ter }}=\frac{\varphi_{j}\left(T, p, x_{j}, x_{\text {solv }, 1} / x_{\text {solv }, 2}=C\right)}{\varphi_{j}^{\infty}\left(T, p, x_{ \pm} \rightarrow 0, x_{\text {solv }, 1} / x_{\text {solv }, 2}=C\right)} \tag{18}
\end{equation*}
$$

In the literature, MIACs are usually given on a molal basis whereas ePC-SAFT is mole-fraction based. Thus, the following conversion from mole fraction ( $x$ )-based to molality ( $m$ )-based values is applied to yield

$$
\begin{equation*}
\gamma_{ \pm}^{*, m}=\frac{\gamma_{ \pm}^{*, x}}{\left(1+\sum_{j} 0.001 v_{j} m_{j} M_{\text {solv }}\right)}=\frac{1}{\left(1+\sum_{j} 0.001 v_{j} m_{j} M_{\text {solv }}\right)}\left[\frac{\varphi_{ \pm}^{x}}{\varphi_{ \pm, x_{ \pm} \rightarrow 0}^{\infty, x}}\right] \tag{19}
\end{equation*}
$$

For systems composed of more than one solvent, the conventional mean molar mass of the solvents is used in Eq. (19)

$$
\begin{equation*}
M_{\text {solv }}=\sum_{i} \chi_{\text {solv }, i}^{\text {salt free }} \cdot M_{\text {solv }, i} \tag{20}
\end{equation*}
$$

where the solvent mole fractions in the salt-free composition $x_{\text {solv }, i}^{\text {salt free }}$ add up to unity.

## 4. Model-parameter estimation

The first step in calculating properties of multicomponent electrolyte solutions is to model each solvent and solvent mixture as accurately as possible. The parameters of water were obtained from Cameretti and Sadowski (2008) and are shown in Table 5 along with the alcohols' pure-component and solvent-mixture parameters.

Using these parameters, the average relative deviations for the liquid densities of pure water, pure methanol, and pure ethanol (Eq. (21)) are $0.06 \%(273-373 \mathrm{~K}), 2.01 \%(200-512 \mathrm{~K})$, and $0.79 \% (230-516 \mathrm{~K})$, respectively. For the vapor pressures of the pure solvents, the deviations are $0.52 \%, 2.36 \%$, and $0.99 \%$ for water, methanol, and ethanol, respectively, within the same temperature range as before.

To calculate the thermodynamic properties in binary alcohol/ electrolyte systems, the two-ion parameters, diameter ( $\sigma_{j}$ ) and dispersion energy $\left(u_{j} / k_{B}\right)$, were determined by simultaneous fitting to experimental solution densities and osmotic coefficients of various electrolyte solutions. For each of the alcohols (methanol and ethanol), solutions of electrolytes AnCat with $\mathrm{Cat}^{+}=\left[\mathrm{Na}^{+}, \mathrm{Li}^{+}, \mathrm{K}^{+}, \mathrm{NH}_{4}^{+}\right]$and $\mathrm{An}^{-}=\left[\mathrm{Cl}^{-}, \mathrm{Br}^{-}, \mathrm{I}^{-}\right]$were considered. All parameters obtained this way are listed in Table 6. Note that these parameter sets are only valid for the solvent parameters given in Table 5. It should be mentioned that using ion parameters which are non-specific to the solvent is not possible for accurate modeling results. This strategy would allow for applying $k_{i j}$ values between solvents and ions to correct for the dispersion energy without correcting for the ions' diameters. The failure of this procedure is most probably due to the DH-term which requires solvated-ion instead of bare-ion diameters.

In general, the ion parameters in the alcohol systems reflect the sequence of the ion parameters in aqueous solutions, which were already available (Held et al., 2008). The solvated-ion diameters correlate with the ions' molar masses, i.e., $\sigma_{\mathrm{K}+}>\sigma_{\mathrm{Na}+}>\sigma_{\mathrm{Li}+}$ for the cations and $\sigma_{\mathrm{I}-}>\sigma_{\mathrm{Br}-}>\sigma_{\mathrm{Cl}-}$ for the anions. However, the ammonium ion represents an exception, possibly because $\mathrm{NH}_{4}^{+}$is not a monoatomic ion. Generally, the ion diameters are larger in the alcohols than in water, which physically means that alcohols cannot approach the ions as close as water molecules can. Moreover, $\sigma_{j}$, which represents the solvated (ion plus solvation shell) diameter rather than the bare-ion diameter, is also larger due to the fact that alcohols are larger than water.

The ion dispersion-energy parameter reflects the strength of solvation (Held et al., 2008). For cations, small cations are the most

Table 6
Ion parameters used in this work; only valid for solvent parameters listed in Table 5.
| Ion | in water (taken from Held et al., 2008) |  | in methanol |  | in ethanol |  |
| :--- | :--- | :--- | :--- | :--- | :--- | :--- |
|  | $\sigma_{j}[\AA]$ | $u_{j} / k_{B}[\mathrm{~K}]$ | $\sigma_{j}[\AA]$ | $u_{j} / k_{B}[\mathrm{~K}]$ | $\sigma_{j}[\AA]$ | $u_{j} / k_{B}$ [K] |
| $\mathrm{Li}^{+}$ | 1.82 | 2697.3 | 3.60 | 1000.0 | 3.60 | 1000.0 |
| $\mathrm{Na}^{+}$ | 2.41 | 646.1 | 4.00 | 400.0 | 3.70 | 550.0 |
| $\mathrm{K}^{+}$ | 2.97 | 271.1 | 4.50 | 150.0 | $3.80{ }^{\text {a }}$ | $40.0{ }^{\text {a }}$ |
| $\mathrm{NH}_{4}^{+}$ | 3.48 | 176.4 | 3.00 | 400.0 | 3.00 | 400.0 |
| $\mathrm{Cl}^{-}$ | 3.06 | 47.3 | 3.70 | 50.0 | 3.30 | 50.0 |
| $\mathrm{Br}^{-}$ | 3.46 | 60.2 | 3.90 | 90.0 | 3.70 | 60.0 |
| $\mathrm{I}^{-}$ | 3.93 | 80.4 | 4.10 | 140.0 | 3.90 | 120.0 |


[^1]strongly solvated in all solvents because the charge density dominates the short-range, solvent-cation interactions. Thus, cation dispersionenergy parameters decrease in the following order: $u_{\mathrm{Li}+} / k_{B}>u_{\mathrm{Na}+} / k_{\mathrm{B}}>u_{\mathrm{K}+} / k_{\mathrm{B}}$. In contrast, steric effects dominate the solvent-anion interactions so that, for anions, large anions are the most strongly solvated in all solvents. Anion dispersion energy parameters follow the series $u_{\mathrm{I}-} / k_{B}>u_{\mathrm{Br}-} / k_{\mathrm{B}}>u_{\mathrm{Cl}-} / k_{\mathrm{B}}$. As for the aqueous electrolyte solutions, the cations in alcohol systems have a stronger influence on activity coefficients than anions do, which results in larger $u / k_{B}$ values for the cations than for the anions. Thus, the ion parameters differ more among the cation than the anion series. This can also be observed for the experimental osmotic coefficients, which show a more pronounced change upon varying the cation than the anion.

## 5. Modeling results

Using the parameters summarized in Tables 5 and 6, liquid densities, osmotic coefficients (corresponding to solvent activity coefficients), and mean ionic activity coefficients (MIAC) were modeled for binary solvent/salt systems. Furthermore, liquid densities and MIACs were predicted for some ternary (and one quaternary) systems containing water, alcohol(s), and salt. The average absolute deviations (AAD) and the average relative deviations (ARD) between the modeled and experimental thermodynamic properties are calculated by

$$
\begin{align*}
& \mathrm{AAD}=\frac{1}{N P} \sum_{k=1}^{N P}\left|\left(y_{k}^{\text {calc }}-y_{k}^{\exp }\right)\right| \\
& \mathrm{ARD}=100 \cdot \frac{1}{N P} \sum_{k=1}^{N P}\left|\left(1-\frac{y_{k}^{\text {calc }}}{y_{k}^{\text {exp }}}\right)\right| \tag{21}
\end{align*}
$$

### 5.1. Binary electrolyte solutions

The binary systems considered in this work are methanol/salt and ethanol/salt systems. Deviations between experiment and ePC-SAFT modeling are listed in Table 7. Moreover, the results are shown graphically for two systems in Fig. 2. The solution densities and osmotic coefficients are modeled with absolute relative deviations of $0.87 \%$ and $7.19 \%$, respectively, for nine binary ethanol/salt solutions. When considering ten methanol/salt systems, the deviations are slightly higher for the solution densities $(0.92 \%)$ and much better for the osmotic coefficients ( $4.26 \%$ ). Due to their very low solubility in the alcohols, no data were available for KCl or $\mathrm{NH}_{4} \mathrm{Cl}$. Note that the modeling results are slightly better for the alcohol/salt systems than for the binary water/salt solutions due to (1) the lower concentration range and (2) the lower number of considered electrolytes for the alcohol systems.

### 5.1.1. Modeling thermodynamic properties of electrolyte alcohol solutions

Fig. 2 graphically illustrates the characteristic salt influence on thermodynamic components/solution properties in ethanol and methanol electrolyte solutions for two systems. Solution densities (Fig. 2a), vapor pressures (b), osmotic coefficients (c), and MIACs (d) were modeled with a single ionic solvent-specific parameter set.

Fig. 2a exemplarily illustrates the influence of NaI on liquid alcohol solution densities. As seen in Fig. 2, the experimental data are described accurately up to salt concentrations of $2 \mathrm{~mol} / \mathrm{kg}$.

Also, the vapor pressures in these systems are modeled with high precision (Fig. 2b). However, such data is not very sensitive because the vapor pressure reflects the solvent activity, which is the product of solvent concentration and the solvent activity coefficient. Thus, osmotic coefficients and mean ionic activity

Table 7
Deviations of ePC-SAFT modeling from experimental solution density, osmotic coefficients, and MIAC data in the binary methanol and ethanol electrolyte solutions.
| Salt | Solution density |  |  |  |  | Mean ionic activity coefficient |  |  |  |  | Osmotic coefficient |  |  |  |  |
| :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- |
|  | NP | $m_{\text {max }}$ [mol/kg] | ARD [\%] | AAD [ $\mathrm{kg} / \mathrm{m}^{3}$ ] | Ref. | NP | $m_{\text {max }}$ [mol/kg] | ARD [\%] | AAD [-] | Ref. | NP | $m_{\text {max }}$ [mol/kg] | ARD [\%] | AAD [-] | Ref. |
| Methanol |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |
| LiI | - | - | - | - | - | 5 | 0.60 | 11.41 | 0.04 | a | 5 | 0.60 | 8.56 | 0.08 | a |
| LiBr | 18 | 2.03 | 1.10 | 9.40 | Pasztor and Criss (1978) | 15 | 2.05 | 8.06 | 0.04 | Zafarani-Moattar and Nasirzade (1998) | 17 | 0.20 | 1.99 | 0.02 | Zafarani-Moattar and Nasirzade (1998) |
| LiCl | 9 | 0.32 | 0.39 | 3.09 | Pasztor and Criss (1978) | 18 | 2.19 | 11.89 | 0.05 | Zafarani-Moattar and Nasirzade (1998) | 27 | 0.25 | 2.85 | 0.03 | Zafarani-Moattar and Nasirzade (1998) |
| NaI | 21 | 0.81 | 0.72 | 5.98 | Lankford et al. (1984) | 5 | 1.00 | 4.94 | 0.02 | Barthel et al. (1985) | 14 | 2.76 | 2.95 | 0.04 | Barthel et al. (1985) |
| NaBr | 13 | 1.10 | 1.32 | 11.21 | Pasztor and Criss (1978) | 38 | 0.65 | 5.11 | 0.02 | Barthel et al. (1985), Ye et al. (1994) | 18 | 0.65 | 4.50 | 0.04 | Barthel et al. (1985, Ye et al. (1994) |
| NaCl | 14 | 0.13 | 0.46 | 3.62 | Pasztor and Criss (1978) | 24 | 0.22 | 3.88 | 0.02 | Barthel et al. (1985) | 23 | 0.22 | 4.86 | 0.04 | Barthel et al. (1985) |
| $\mathrm{NH}_{4} \mathrm{I}$ | 3 | 0.55 | 1.01 | 8.48 | Gonzalez et al. (2005) | 6 | 0.80 | 12.93 | 0.04 | ${ }^{\text {a }}$ | 6 | 0.80 | 6.47 | 0.05 | a |
| $\mathrm{NH}_{4} \mathrm{Br}$ | 5 | 1.03 | 1.77 | 14.74 | Ortmaier (1989) | 6 | 0.79 | 11.64 | 0.04 | ${ }^{\text {a }}$ | 6 | 0.79 | 5.65 | 0.05 | ${ }^{\text {a }}$ |
| KI | 15 | 0.38 | 0.87 | 7.09 | Pasztor and Criss (1978) | 25 | 0.62 | 6.74 | 0.03 | Barthel et al. (1985) | 41 | 0.62 | 2.37 | 0.02 | Barthel et al. (1985) |
| KBr | 14 | 0.16 | 0.61 | 4.85 | Pasztor and Criss (1978) | 14 | 0.13 | 7.93 | 0.04 | Barthel et al. (1985) | 14 | 0.13 | 2.42 | 0.02 | Barthel et al. (1985) |
| $\emptyset$ |  |  | 0.92 | 7.61 |  |  |  | 8.45 | 0.03 |  |  |  | 4.26 | 0.04 |  |
| Ethanol |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |
| LiI | 8 | 1.22 | 0.95 | 8.11 | ${ }^{\mathrm{a}}$ | 5 | 1.25 | 9.63 | 0.02 | ${ }^{\text {a }}$ |  | 51.25 | 6.01 | 0.06 | ${ }^{\mathrm{a}}$ |
| LiBr | 7 | 1.57 | 1.89 | 15.99 | ${ }^{\mathrm{a}}$ | 6 | 1.50 | 12.35 | 0.02 | ${ }^{\text {a }}$ | 6 | 1.50 | 3.06 | 0.03 | ${ }^{\mathrm{a}}$ |
| LiCl | 5 | 2.05 | 2.37 | 19.59 | ${ }^{\mathrm{a}}$ | 4 | 1.00 | 10.37 | 0.02 | ${ }^{\mathrm{a}}$ | 4 | 1.00 | 8.13 | 0.06 | ${ }^{\mathrm{a}}$ |
| NaI | 8 | 1.09 | 0.98 | 8.31 | ${ }^{\mathrm{a}}$ | 7 | 2.00 | 13.92 | 0.06 | Barthel and Lauermann (1986) | 7 | 2.00 | 3.17 | 0.03 | Barthel and Lauermann (1986) |
| NaBr | 4 | 0.22 | 0.49 | 3.89 | ${ }^{\mathrm{a}}$, Lacmann and Synowietz (1977) | 3 | 0.23 | 17.11 | 0.06 | ${ }^{\text {a }}$ | 3 | 0.23 | 10.17 | 0.08 | ${ }^{\mathrm{a}}$ |
| NaCl | 5 | 0.01 | 0.03 | 0.22 | Lacmann and Synowietz (1977) | 1 | 0.01 | 13.22 | 0.09 | ${ }^{\text {a }}$ | 1 | 0.01 | 8.07 | 0.08 | ${ }^{\mathrm{a}}$ |
| $\mathrm{NH}_{4} \mathrm{I}$ | 19 | 1.12 | 0.70 | 6.04 | ${ }^{\mathrm{a}}$, Gonzalez et al. (2005) | 4 | 1.00 | 18.75 | 0.03 | a | 4 | 1.00 | 5.27 | 0.03 | ${ }^{\mathrm{a}}$ |
| $\mathrm{NH}_{4} \mathrm{Br}$ | 6 | 0.30 | 0.33 | 2.58 | ${ }^{\mathrm{a}}$ | 3 | 0.30 | 19.19 | 0.05 | ${ }^{\text {a }}$ | 3 | 0.30 | 11.84 | 0.09 | ${ }^{\mathrm{a}}$ |
| $\mathrm{NH}_{4} \mathrm{Cl}$ | 6 | 0.14 | 0.12 | 0.94 | ${ }^{\mathrm{a}}$ | 2 | 0.09 | 21.57 | 0.09 | ${ }^{\text {a }}$ | 2 | 0.09 | 8.97 | 0.07 | ${ }^{\mathrm{a}}$ |
| $\emptyset$ |  |  | 0.87 | 7.30 |  |  |  | 15.12 | 0.05 |  |  |  | 7.19 | 0.06 |  |


[^2]![](https://cdn.mathpix.com/cropped/137993f5-9955-4aaa-862d-a06927e0cd5f-08.jpg?height=817&width=1532&top_left_y=191&top_left_x=284)
Fig. 2. Thermodynamic properties of NaI alcohol solutions modeled with ePC-SAFT (lines) compared to the experimental data (Refs. see Table 7) in ethanol (open circles) and methanol solutions (full circles) at $25^{\circ} \mathrm{C}$. (a) Solution densities. (b) Vapor pressures. (c) Osmotic coefficients. (d) Mean ionic activity coefficients. Model parameters and deviations are given in Tables 6, 7.

![](https://cdn.mathpix.com/cropped/137993f5-9955-4aaa-862d-a06927e0cd5f-08.jpg?height=668&width=869&top_left_y=1214&top_left_x=146)
Fig. 3. Influence of NaI on solvent activity coefficients at $25^{\circ} \mathrm{C}$. The symbols (stars: water, squares: methanol, circles: ethanol) represent experimental data from Barthel and Lauermann (1986), Barthel et al. (1985), Lobo and Quaresma (1989), and Tomasula et al. (1985). The lines are calculations with ePC-SAFT.

coefficients (MIACs) were introduced because they distinctly deviate from unity (see also Fig. 1). Generally, the modeled osmotic coefficients and MIACs agree well with the experimental data, as seen in Fig. 2c and d for the ethanol/ NaI and methanol/NaI systems, respectively.

### 5.1.2. Comparing the influence of salt on different solvents

Based on the experimental data and ePC-SAFT modeling, it is possible to compare the influence of a salt on both alcohols considered in this work with the influence on water (modeled previously by Held et al., 2008).

Fig. 3 illustrates the influence of sodium iodide on the activity coefficients of water, methanol, and ethanol. In general, the influence of salt on alcohol activity coefficients is much stronger than on
water activity coefficients in both the diluted and the more concentrated regions. Moreover, the solvent activity coefficients are high for alcohol and low for water in the diluted region, which is not the case for the more concentrated solutions. Because of this phenomenon, the activity coefficients intersect. Up to a molality of approximately 0.9 m , the order is $\gamma_{\text {water }}<\gamma_{\text {methanol }}<\gamma_{\text {ethanol }}$, whereas the series is reversed at higher concentrations. Independent of the considered solvent, two effects appear to be at play in electrolyte solutions; the first effect leads to an increase in solvent activity coefficients, and the second effect causes a decrease in activity coefficients.

- At very small to moderate concentrations ( $<0.5 \mathrm{~m}$ ) the interionic Coulombic forces dominate the solution's behavior due to their long range. The strength of this interaction is determined by the capability of the solvent to shield ionic charges, expressed by the solvent's dielectric constant, $\varepsilon$. The lower the $\varepsilon$, the stronger the interionic interactions, i.e., the DebyeHückel effect becomes more pronounced, increasing the solvent activity coefficient (more solvent molecules are forced to pass into the vapor phase compared to pure solvent). The dielectric constants of the pure solvents are as follows: $\varepsilon_{\text {water }}>\varepsilon_{\text {methanol }}>\varepsilon_{\text {ethanol }}$, while the respective activity coefficients follow the reverse order.
- Further addition of salt overlays the Coulomb forces by shortrange solvation effects, which results in lower solvent activity coefficients (more solvent molecules are retained in the liquid phase than in pure solvent). The solvents considered here (water, alcohols) possess a strong hydrogen-bond network. To solvate the introduced ions, such associating solvents must break their hydrogen-bond networks. As the strength of these networks decreases from water to ethanol, the activity coefficients follow that sequence into the more concentrated regions ( $m>2$ ).

This complex behavior is well captured with ePC-SAFT; the model describes the experimental data in Fig. 3 with good accuracy.

### 5.2. Ternary mixed-solvent electrolyte solutions

The solvent-specific ion parameters obtained from estimating the parameters in binary solvent/electrolyte systems (Table 6) were used to directly predict properties in ternary mixed-solvent solutions. Although this approach may result in a less accurate description of the experimental data than directly fitting the ion parameters to the respective ternary properties of interest, this approach has three main advantages: (1) the number of model parameters is reduced because they are universal parameters; (2) ternary properties of any system can be predicted based on the binary modeling results; and (3) the binary solvent/salt systems can still be described in contrast to correlations, which are only valid for the considered ternary systems.

To describe ternary mixtures, the conventional BerthelotLorenz mixing rules (Eqs. (11) and (12)) were applied for the hard-chain and dispersion-energy contribution in Eq. (7). In addition to those standard combining rules, the hard-sphere and the Debye-Hückel contribution for ions are required for the solvated-ion diameters in the solvent mixture. Here, we average the ion diameters obtained for the different solvents (given in Table 6) using the salt-free solvent composition yielding the following for ion $j$ in a water/ethanol mixture:

$$
\begin{equation*}
\sigma_{j}^{\mathrm{H}_{2} \mathrm{O} / \mathrm{EtOH}}=x_{\mathrm{H}_{2} \mathrm{O}}^{\text {salt free }} \cdot \sigma_{j}^{\mathrm{H}_{2} \mathrm{O}}+x_{\mathrm{EtOH}}^{\text {salt free }} \cdot \sigma_{j}^{\mathrm{EtOH}} \tag{22}
\end{equation*}
$$

Note that such a mixing rule is not needed for the ion's dispersion energy parameter, as any dispersion interaction between ions is neglected within ePC-SAFT.

As for the binary solvent/salt systems (Table 7), the deviations between modeling and experiments for all considered ternary solutions are listed in Table 8.

### 5.2.1. Solution densities

In contrast to $g^{\mathrm{E}}$ models, equations of state make it possible to model solution densities. Takenaka et al. (1994) measured the influence of alkali chlorides on water/methanol solution densities at temperatures between 15 and $45^{\circ} \mathrm{C}$. As a typical example, the liquid densities of water/ethanol/Nal solutions between mole fractions of water of 0.41 and 0.79 in the salt-free system at $25^{\circ} \mathrm{C}$ are shown in Fig. 4. Fig. 4 shows that ePC-SAFT describes the experimental data with high accuracy at all solvent compositions and salt concentrations examined. Although the ion parameters were adjusted to binary solvent/salt densities at $25{ }^{\circ} \mathrm{C}$, ePC-SAFT correctly predicted the respective densities (1) in the ternary two-solvent systems at $25^{\circ} \mathrm{C}$ and even better (2) at other temperatures (results not shown). The deviation (ARD) between the model and experiment is $0.63 \%$ (or $5.6 \mathrm{~kg} / \mathrm{m}^{3}$ ) for water/ethanol/ NaI solutions (shown in Fig. 4) and $0.88 \%\left(8.1 \mathrm{~kg} / \mathrm{m}^{3}\right)$ for all considered systems (see Table 8).

### 5.2.2. Mean ionic activity coefficients

The most important property that characterizes an electrolyte in solution is the mean ionic activity coefficient. This property deviates much more from unity than the activity coefficient of the solvent. The MIAC in the binary solvent/salt solutions (see Fig. 2) first decreases with increasing concentration; after reaching a minimum, it often increases with increasing salt concentration, reaching very high values in some cases. This behavior is also observed for the ternary systems considered in this work. For NaCl in water/ethanol solutions at $25^{\circ} \mathrm{C}$, Fig. 5 illustrates a distinct decrease in the MIAC at small salt concentrations. However, at higher NaCl molalities ( $>0.7 \mathrm{~m}$ ), the MIAC remains almost constant. Moreover, the MIAC decreases as the ethanol concentration increases due to the lower permittivity of ethanol (i.e., weaker solvent/ion interactions) and is valid for NaCl within

Table 8
Deviations between ePC-SAFT modeling and experimental solution density and MIAC data in ternary aqueous methanol and ethanol electrolyte solutions at $25^{\circ} \mathrm{C}$. $x_{\mathrm{H} 2 \mathrm{O}}$ is the mole fraction of water in the salt-free system.
| Density |  |  | Mean ionic activity coefficient |  |  |
| :--- | :--- | :--- | :--- | :--- | :--- |
| $x_{\mathrm{H} 2 \mathrm{O}}$ [\%] | $m_{\text {max }}$ [mol/kg] | ARD [\%] | $x_{\mathrm{H} 2 \mathrm{O}}$ [\%] | $m_{\text {max }}$ [mol/kg] | ARD [\%] |
| Ethanol/Water/KCl (Galleguillos et al., 2003; Lopes et al., 1999) |  |  |  |  |  |
| 95 | 2.6 | 0.7 | 98 | 4.0 | 2.5 |
| 86 | 1.5 | 1.1 | 94 | 3.0 | 0.9 |
| 81 | 1.3 | 1.1 | 91 | 3.0 | 3.0 |
| 72 | 0.6 | 0.7 | - | - | - |
| 61 | 0.4 | 0.5 | - | - | - |
| Ethanol/Water/NaBr (Gonzalezdiaz et al., 1995) |  |  |  |  |  |
| - | - | - | 91 | 1.0 | 1.6 |
| - | - | - | 79 | 0.8 | 4.1 |
| - | - | - | 63 | 0.8 | 7.4 |
| - | - | - | 39 | 0.8 | 9.7 |
| Ethanol/Water/NaCl (Esteso et al., 1989; Galleguillos et al., 2003) |  |  |  |  |  |
| 95 | 4.3 | 1.2 | 91 | 2.0 | 1.2 |
| 86 | 3.0 | 1.4 | 79 | 1.5 | 1.6 |
| 81 | 2.4 | 1.2 | 63 | 1.0 | 2.2 |
| 72 | 1.5 | 0.5 | 39 | 0.1 | 5.0 |
| 60 | 0.5 | 0.4 | - | - | - |
| Ethanol/Water/NaI (Taniewska-Osinska and Chadzynski, 1976) |  |  |  |  |  |
| 90 | 1.0 | 0.5 | - | - | - |
| 79 | 0.8 | 0.6 | - | - | - |
| 63 | 0.7 | 0.3 | - | - | - |
| 41 | 0.6 | 0.4 | - | - | - |
| 11 | 0.5 | 1.1 | - | - | - |
| Methanol/Water/KCl (Basili et al., 1999; Takenaka et al., 1994) |  |  |  |  |  |
| 87 | 0.3 | 0.6 | 88 | 2.0 | 1.5 |
| 69 | 0.3 | 1.0 | 73 | 1.5 | 1.5 |
| 43 | 0.1 | 1.6 | 54 | 0.6 | 3.3 |
| 8 | 0.1 | 0.9 | 31 | 0.2 | 2.5 |
| Methanol/Water/NaCl (Takenaka et al., 1994; Yan et al., 1994) |  |  |  |  |  |
| 87 | 0.8 | 0.4 | 88 | 3.7 | 3.9 |
| 69 | 0.4 | 0.9 | 73 | 1.7 | 2.3 |
| 43 | 0.4 | 1.3 | 54 | 0.8 | 4.6 |
| 8 | 0.2 | 0.7 | 31 | 0.2 | 0.7 |
| - | - | - | 17 | 0.1 | 4.6 |
| Methanol/Water/NaBr (Ye et al., 1994) |  |  |  |  |  |
| 87 | 0.3 | 0.6 | 88 | 3.0 | 3.6 |
| 69 | 0.3 | 1.0 | 73 | 3.0 | 3.6 |
| 43 | 0.1 | 1.5 | 54 | 3.0 | 5.3 |
| 8 | 0.1 | 0.9 | 31 | 2.5 | 3.3 |
| Methanol/Water/LiCl (Mussini et al., 2000; Takenaka et al., 1994) |  |  |  |  |  |
| 87 | 0.8 | 0.7 | 73 | 9.0 | 6.3 |
| 69 | 1.2 | 1.0 | 54 | 5.0 | 6.4 |
| 43 | 0.5 | 1.4 | 31 | 2.0 | 3.7 |
| 8 | 0.4 | 0.6 | - | - | - |


the whole solubility range. This behavior is also observed for the MIACs of the other salts in this work (see Table 8). In general, the salt behaves in the solvent mixture as one would expect. The experimental MIACs in the mixture are between the MIACS of the binary $\mathrm{NaCl} /$ water and $\mathrm{NaCl} /$ ethanol solutions
ePC-SAFT is able to describe this behavior within the respective NaCl solubility limits with an average deviation (ARD) of $6.7 \%$ between 20 and $80 \mathrm{wt} \%$ ethanol in the salt-free solvent system. As for the modeled solution densities (Fig. 4), no additional adjustable parameters were introduced so that these results were also pure predictions. Moreover, due to the low solubility of NaCl in ethanol, the $\mathrm{Na}^{+}$and $\mathrm{Cl}^{-}$ePC-SAFT parameters were not adjusted to osmotic or activity coefficients of NaCl in ethanol, but rather to other salts containing $\mathrm{Na}^{+}$and $\mathrm{Cl}^{-}$(i.e., to $\mathrm{NaI}, \mathrm{LiCl}$, and $\mathrm{NH}_{4} \mathrm{Cl}$ ). This again reveals the predictive capability of ePC-SAFT and the advantage of applying ion-specific model parameters.

![](https://cdn.mathpix.com/cropped/137993f5-9955-4aaa-862d-a06927e0cd5f-10.jpg?height=502&width=873&top_left_y=215&top_left_x=146)
Fig. 4. Experimental and modeled solution densities of water/ethanol/NaI mixtures at $25^{\circ} \mathrm{C}$ and water fractions in the salt-free system between 41 and $79 \mathrm{~mol} \%$. The symbols represent experimental data from Taniewska-Osinska and Chadzynski (1976). The full lines are predictions with ePC-SAFT. The gray lines represent the density of the binary water/NaI and ethanol/NaI solutions, respectively.

![](https://cdn.mathpix.com/cropped/137993f5-9955-4aaa-862d-a06927e0cd5f-10.jpg?height=526&width=875&top_left_y=958&top_left_x=148)
Fig. 5. Experimental and modeled mean ionic activity coefficients of NaCl in water/ethanol mixtures at $25^{\circ} \mathrm{C}$ and salt-free water fractions between 39 and $91 \mathrm{~mol} \%$. The symbols represent experimental data taken from Esteso et al. (1989) The full lines are predictions with ePC-SAFT. The gray lines represent the MIACs of NaCl in the binary water/ NaCl and ethanol/ NaCl solutions, respectively.

Fig. 6 illustrates the concentration dependence of NaBr MIACs in different water/methanol mixtures at $25{ }^{\circ} \mathrm{C}$. In this case, $\mathrm{Na}^{+}$and $\mathrm{Br}^{-}$parameters were adjusted to the properties of binary salt/ solvent systems that also contained NaBr . The MIAC values in the solvent mixtures can be predicted with a deviation (ARD) of $3.6 \%$ between 0.31 and 0.54 mol fraction of water (salt-free solvent composition), which is slightly better than that shown in Fig. 5. This is expected as the $\mathrm{Na}^{+}$and $\mathrm{Br}^{-}$parameters in water-free salt/ methanol solutions were fitted to solutions also containing NaBr .

The properties of other aqueous ethanol and methanol electrolyte solutions were also investigated. ePC-SAFT was able to quantitatively predict densities and MIACs in these systems at $25^{\circ} \mathrm{C}$. Although not taken into account in the parameter estimation, densities can also be predicted accurately at temperatures other than $25^{\circ} \mathrm{C}$. Overall deviations between model predictions and experiments are given in Table 8. For systems with higher ARD values, the deviations may be caused by (1) model uncertainties stemming from the binary solvent/ salt system, (2) the fact that the absolute MIAC values are relatively small (slight errors cause high relative deviations), and (3) measurement uncertainties. Examples for which experimental data are doubtful are the measurements by Yang et al. (2007) and by Deyhimi et al. (2003) because their data are not within the MIACs of the respective binary salt/alcohol solutions. These systems (water/methanol/KI and water/methanol/ $\mathrm{NH}_{4} \mathrm{Cl}$ ) were not considered in this work.

![](https://cdn.mathpix.com/cropped/137993f5-9955-4aaa-862d-a06927e0cd5f-10.jpg?height=522&width=868&top_left_y=200&top_left_x=1085)
Fig. 6. Experimental and modeled mean ionic activity coefficients of NaBr in water/methanol mixtures at $25^{\circ} \mathrm{C}$ and salt-free water fractions between 31 and $54 \mathrm{~mol} \%$. The symbols represent experimental data taken from Ye et al. (1994). The full lines are predictions with ePC-SAFT. The gray lines represent the MIACs of NaBr in the binary water/ NaBr and methanol/ NaBr solutions, respectively.

![](https://cdn.mathpix.com/cropped/137993f5-9955-4aaa-862d-a06927e0cd5f-10.jpg?height=520&width=870&top_left_y=960&top_left_x=1083)
Fig. 7. Experimental and modeled mean ionic activity coefficients of NaCl in water/methanol/ethanol mixtures at $25^{\circ} \mathrm{C}$. The symbols represent experimental data (Hernandez-Hernandez, et al., 2007) at the following salt-free solvent compositions $(\mathrm{mol} \%): \quad x_{\text {water }}=0.918, \quad x_{\mathrm{MeOH}}=0.061 \quad($ circles $), \quad x_{\text {water }}=0.830$, $x_{\text {MeOH }}=0.100$ (squares), $x_{\text {water }}=0.776, x_{\text {MeOH }}=0.073$ (stars). The full lines are predictions with ePC-SAFT.

### 5.3. Quaternary mixed-solvent electrolyte solutions

Hernandez-Hernandez et al. (2007) recently measured the MIACs of NaCl in solutions containing three solvents (water, methanol, and ethanol) at $25^{\circ} \mathrm{C}$. They also presented values for the dielectric constant of the salt-free solvent mixtures, which is crucial for modeling the Coulombic interactions within ePC-SAFT (see Eq. (8)). Fig. 7 compares their experimental activity coefficients (as a function of the salt-free solvent composition as well as the salt concentration) with the ePC-SAFT predictions. The ePC-SAFT modeling also succeeds in systems with more than two solvents. Other than the parameters presented in Tables 5, 6 and the mixing rule (Eq. (22)), no additional parameters were applied for the model prediction.

The average deviation between model and experiments is $4.17 \%$ (ARD) for the three shown solvent compositions and salt concentrations. However, it is remarkable that the experimental data in Fig. 7 are lower than the ePC-SAFT predictions. In contrast, for the two-solvent systems (water/ethanol/ NaCl and water/methanol/ NaCl ), the experimental data were higher than the predicted MIACs (compare Fig. 5). This might be due to uncertainties in the dielectric constants, which strongly influence the results. Moreover, measurement uncertainties may also result in slight deviations.

In general, the model shows very good agreement with the results, which reveals that ePC-SAFT can be applied as a predictive
model to calculate densities and MIACs in multi-solvent electrolyte systems. Recently, Wang et al. (2011) presented the COSMOSAC model, which was used to predict MIACs in multi-solvent electrolyte solutions. Although both approaches (ePC-SAFT and COSMO-SAC) are predictive and provided promising results, ePCSAFT and the suggested mixing rules in this work appear to be more quantitative than COSMO-SAC. Moreover, Ingram et al. (2011) recently modeled MIACs with the COSMO-RS model. They showed for two alcohol/water/salt systems, that their approach yields also a good qualitative prediction with the experimental data. Based on this, they could predict the salt influence on liquidliquid equilibria of aqueous solvent systems.

## 6. Conclusions

In this study, the ePC-SAFT model was applied to describe thermodynamic properties of alcohol and mixed-solvent electrolyte solutions. Solvent-specific ion parameters were used, which are independent of the electrolyte the ions are part of. First, two parameters, namely the diameter $\sigma_{j}$ and the dispersion energy $u_{j} / k_{B}$ of the solvated ions, were adjusted to the alcoholic solution densities and osmotic-coefficient data. The two parameters possess a physical meaning and show reasonable trends within the ion series. Using this approach, thermodynamic properties of nine ethanol and ten methanol solutions were modeled with reasonable overall ARDs for solution densities (ethanol/salt: $0.87 \%$, methanol/salt: 0.92\%), osmotic coefficients (ethanol/salt: 7.19\%, methanol/salt: $4.24 \%$ ), and mean ionic activity coefficients (etha$\mathrm{nol} /$ salt: $15.19 \%$, methanol/salt: $8.49 \%$ ).

The experimental data used for the model parameterization in methanol/salt solutions were mostly taken from the literature. Due to the lack or inconsistency of available data, especially for the ethanol/ salt systems, a U-tube densimeter and vapor-pressure osmometry were used to measure the densities and MIACs in binary solutions.

Based on the binary solvent/salt systems, the densities and MIACs in ternary two-solvent electrolyte mixtures could be accurately modeled with ARDs of $0.12 \%$ for solution densities and $3.57 \%$ for MIACs. Because no additional parameters were applied for these ternary systems, quantitative predictions were possible with ePCSAFT, i.e., the calculations are based on the binary solvent/salt mixtures only. Once the solution's dielectric constant (see Table 4 or Wohlfarth, 1991) and the binary ion/solvent parameters (see Table 6) are available, these predictions can also be used for mixtures containing more than two solvents. In this work, we showed that the MIACs of NaCl in a solvent mixture of water, methanol, and ethanol can be quantitatively predicted using ePC-SAFT with an ARD of only $4.17 \%$.

## Nomenclature

| $a_{i}[-]$ | activity of component $i$ |
| :--- | :--- |
| $a$ [J] | Helmholtz energy per number of particles |
| $k_{B}[\mathrm{~J} / \mathrm{K}]$ | Boltzmann constant, $1.38065 \times 10^{-23} \mathrm{~J} / \mathrm{K}$ |
| $k_{i j}[-]$ | binary interaction parameter, dispersion |
| $k_{i j}^{\mathrm{hb}}$ [-] | binary interaction parameter, association |
| $l_{i j}$ [-] | binary interaction parameter, segment diameter |
| $m[\mathrm{~mol} / \mathrm{kg}]$ molality (moles solute $i$ per kg solvent) |  |
| $M[\mathrm{~g} / \mathrm{mol}]$ molecular weight |  |
| $m_{\text {seg }}$ [-] | number of segments |
| $N[-]$ | number of association sites |
| NP [-] | number of data points |
| $p$ [bar] | pressure |
| $q_{j}$ [C] | charge of ion $j$ |
| $T[\mathrm{~K}]$ | temperature |
| $u / k_{B}[\mathrm{~K}]$ | dispersion-energy parameter |
| $x$ [-] | mole fraction |

## Greek letters

| $\gamma_{i}[-]$ | symmetric activity coefficient of component $i$ (related to pure component) |
| :--- | :--- |
| $\gamma_{i}^{*}[-]$ | asymmetric activity coefficient of component $i$ (related to infinite dilution) |
| $\varphi_{i}[-]$ | fugacity coefficient of component $i$ |
| $\Phi$ [-] | osmotic coefficient |
| $\varepsilon[\mathrm{C} / \mathrm{Vm}]$ | dielectric constant of a medium, $\varepsilon_{\mathrm{r}} \cdot \varepsilon_{0}$ |
| $\varepsilon_{\mathrm{r}}$ [-] | relative permittivity |
| $\varepsilon_{0}$ [C/Vm] |  |
| $\varepsilon_{\mathrm{hb}}^{A i \mathrm{Bi}} / k_{B}[\mathrm{~K}]$ |  |
| $\kappa[1 / \AA]$ |  |
| $\kappa_{\mathrm{hb}}^{A i B i}$ [-] <br> $\rho\left[\mathrm{kg} / \mathrm{m}^{3}\right]$ | association-volume parameter |
|  |  |
| $v$ [-] | stoichiometric factor |
| $\sigma_{i}$ [A] | temperature-independent segment diameter of molecule $i$ |
| $\chi$ [-] | abbreviation for an expression used in Eq. (8) |

## Subscripts

| An, - | anion |
| :--- | :--- |
| Cat, + | cation |
| E, EtOH | ethanol |
| $i, j, k$ | component indexes |
| M | methanol |
| seg | segment |
| W, $\mathrm{H}_{2} \mathrm{O}$ | water |
| $\pm$ | mean ionic |
| 0 | pure substance |
| solv | solvent |

## Superscripts

| assoc | association |
| :--- | :--- |
| calc | calculated |
| disp | dispersion |
| exp | experimental |
| hc | hard chain |
| ion | ionic |
| m | based on molality |
| res | residual |
| $x$ | based on mole fraction |
| ,+- | positive or negative charge |
| $\infty$ | infinitely diluted |

## Abbreviations

| AnCat | electrolyte of the form anion-cation |
| :--- | :--- |
| ARD | average relative deviation, defined in Eq. (21) |
| DH | Debye-Hückel |
| EOS | equation of state |
| $g^{\text {E }}$ | excess Gibbs energy |
| MIAC | mean ionic activity coefficient |
| MSA | mean spherical approximation |
| VLE | vapor-liquid equilibrium |

## Acknowledgements

The authors gratefully acknowledge financial support from the German Society of Industrial Research (AiF) (Grant 16295N/1).

We also wish to thank our laboratory assistants, Susanne Richter and Anna Jurytko, for their help with the osmometer measurements.

## References

Ball, F.X., Planche, H., Fürst, W., Renon, H., 1985. Representation of deviation from ideality in concentrated aqueous-solutions of electrolytes using a mean spherical approximation molecular-model. AICHE Journal 31, 1233-1240.
Barthel, J., Lauermann, G., 1986. Vapor pressure measurements on non-aqueous electrolyte solutions. Part 3: solutions of sodium iodide in ethanol, 2-propanol, and acetonitrile. Journal of Solution Chemistry 15, 869-877.
Barthel, J., Neueder, R., Lauermann, G., 1985. Vapor-pressures of non-aqueous electrolyte-solutions.1. Alkali-metal salts in methanol. Journal of Solution Chemistry 14, 621-633.
Basili, A., Mussini, P.R., Mussini, T., Rondinini, S., Sala, B., Vertova, A., 1999. Transference numbers of alkali chlorides and characterization of salt bridges for use in methanol plus water mixed solvents. Journal of Chemical and Engineering Data 44, 1002-1008.
Cameretti, L.F., Sadowski, G., 2008. Modeling of aqueous amino acid and polypeptide solutions with PC-SAPT. Chemical Engineering and Processing 47, 1018-1025.
Cameretti, L.F., Sadowski, G., Mollerup, J.M., 2005. Modeling of aqueous electrolyte solutions with perturbed-chain statistical associated fluid theory. Industrial \& Engineering Chemistry Research 44, 3355-3362 ibid., 8944.
Debye, P., Hückel, E., 1923. Zur Theorie der Elektrolyte. I. Gefrierpunktserniedrigung und verwandte Erscheinungen. Physikalische Zeitschrift 9, 185-206.
Deyhimi, F., Salamat-Ahangari, R., Ghalami-tchoobar, B., 2003. Determination of activity coefficients of $\mathrm{NH}_{4} \mathrm{Cl}$ in methanol-water mixed solvents at 25 degrees C by electromotive force measurements. Physics and Chemistry of Liquids 41, 605-611.
Enick, R.M., Klara, S.M., 1992. Effects of $\mathrm{CO}_{2}$ solubility in brine on the compositional simulation of $\mathrm{CO}_{2}$ floods. SPE Reservoir Engineering 7, 253-258.
Esteso, M.A., Gonzalezdiaz, O.M., Hernandezluis, F.F., Fernandezmerida, L., 1989. Activity-coefficients for NaCl in ethanol-water mixtures at 25-Degrees-C. Journal of Solution Chemistry 18, 277-288.
Fuchs, D., Fischer, J., Tumakaka, F., Sadowski, G., 2006. Solubility of amino acids: influence of the pH value and the addition of alcoholic cosolvents on aqueous solubility. Industrial \& Engineering Chemistry Research 45, 6578-6584.
Fürst, W., Renon, H., 1993. Representation of excess properties of electrolytesolutions using a new equation of state. AICHE Journal 39, 335-343.
Galleguillos, H.R., Taboada, M.E., Graber, T.A., Bolado, S., 2003. Compositions, densities, and refractive indices of potassium chloride plus ethanol plus water and sodium chloride plus ethanol plus water solutions at ( 298.15 and 313.15) K. Journal of Chemical and Engineering Data 48, 405-410.

Gonzalez, B., Dominguez, A., Tojo, J., Esteves, M.J.C., Cardoso, A.J.E.D.M., Barcia, O.E., 2005. Dynamic viscosities of KI or $\mathrm{NH}_{4} \mathrm{I}$ in methanol and $\mathrm{NH}_{4} \mathrm{I}$ in ethanol at several temperatures and 0.1 MPa . Journal of Chemical and Engineering Data 50, 109-112.
Gonzalezdiaz, O.M., Fernandezmerida, L., Hernandezluis, F., Esteso, M.A., 1995. Activity-coefficients for NaBr in ethanol-water mixtures at 25-Degrees-C. Journal of Solution Chemistry 24, 551-563.
Gross, J., Sadowski, G., 2002. Application of the Perturbed-Chain SAFT equation of state to associating systems. Industrial \& Engineering Chemistry Research 41, 5510-5515.
Gross, J., Sadowski, G., 2001. Perturbed-Chain SAFT: an equation of state based on a perturbation theory for chain molecules. Industrial \& Engineering Chemistry Research 40, 1244-1260.
Hall, D.J., Mash, C.J., Pemberton, R.C., 1979. Vapour-liquid equilibrium for the systems water+methanol, water+ethanol, methanol+ethanol and water+ methanol+ethanol, National Physical Laboratory, NPL Report.
Held, C., Cameretti, L.F., Sadowski, G., 2011. Measuring and modeling activity coefficients in aqueous amino-acid solutions. Vol. 50, pp. 131-141.
Held, C., Cameretti, L.F., Sadowski, G., 2008. Modeling aqueous electrolyte solutions. Part1: fully dissociated electrolytes. Fluid Phase Equilibria 270, 87-96.
Held, C., Neuhaus, T., Sadowski, G., 2010. Thermodynamic properties of aqueous ectoine, proline, and urea solutions-measurement and modeling. Biophysical Chemistry 152, 28-39.
Held, C., Sadowski, G., 2009. Modeling aqueous electrolyte solutions. Part2: weak electrolytes. Fluid Phase Equilibria 279, 141-148.
Hernandez-Hernandez, F., Perez-Villasenor, F., Hernandez-Ruiz, V., Iglesias-Silva, G.A., 2007. Activity coefficients of NaCl in $\mathrm{H} 2 \mathrm{O}+\mathrm{MeOH}+\mathrm{EtOH}$ by electromotive force at 298.15 K . Journal of Chemical \& Engineering Data 52, 959-964.
Huang, J.Q., Li, J.D., Gmehling, J., 2009. Prediction of solubilities of salts, osmotic coefficients and vapor-liquid equilibria for single and mixed solvent electrolyte systems using the LIQUAC model. Fluid Phase Equilibria 275, 8-20.
Huang, S.H., Radosz, M., 1990. Equation of state for small, large, polydisperse, and associating molecules. Industrial \& Engineering Chemistry Research 29, 2284-2294.
Ingram, T., Gerlach, T., Mehling, T., Smirnova, I., Extension of COSMO-RS for monoatomic electrolytes: modeling of liquid-liquid equilibria in presence of salts. Fluid Phase Equilibria, doi:10.1016/j.fluid.2011.09.021. In press.
Ji, X.Y., Adidharma, H., 2006. Ion-based SAFT2 to represent aqueous single- and multiple-salt solutions at 298.15 K . Industrial \& Engineering Chemistry Research 45, 7719-7728.

Ji, X.Y., Tan, S.P., Adidharma, H., Radosz, M., 2005. Statistical associating fluid theory coupled with restricted primitive model to represent aqueous strong electrolytes: multiple-salt solutions. Industrial \& Engineering Chemistry Research 44, 7584-7590.
Lacmann, R., Synowietz, C., 1977. Landolt-Börnstein, Group IV Volume 1a: Densities of Liquid Systems: Nonaqueous Systems and Ternary Aqueous Systems. Springer, Berlin.
Lankford, J.I., Holladay, W.T., Criss, C.M., 1984. Isentropic compressibilities of univalent electrolytes in methanol at 25 -degrees-C. Journal of Solution Chemistry 13, 699-720.
Liu, W.B., Li, Y.G., Lu, J.F., 1999. A new equation of state for real aqueous ionic fluids based on electrolyte perturbation theory, mean spherical approximation and statistical associating fluid theory. Fluid Phase Equilibria 160, 595-606.
Lobo, V.M.M., Quaresma, J.L., 1989. Handbook of Electrolyte Solutions, Parts A and B. Elsevier, Amsterdam.

Lopes, A., Farelo, F., Ferra, M.I.A., 1999. Activity coefficients of potassium chloride in water-ethanol mixtures. Journal of Solution Chemistry 28, 117-131.
Mato, F., Cocero, M.J., 1988. Measurement of vapor-pressures of electrolytesolutions by vapor-pressure osmometry. Journal of Chemical and Engineering Data 33, 38-39.
Mussini, P.R., Mussini, T., Sala, B., 2000. Thermodynamics of the cell \{Li-Amalgam\} vertical bar LiX ( $m$ )vertical bar AgX vertical bar Ag ( $\mathrm{X}=\mathrm{Cl}, \mathrm{Br}$ ) and medium effects upon LiX in (acetonitrile plus water), (1,4-dioxane plus water), and (methanol plus water) solvent mixtures with related solvation parameters. Journal of Chemical Thermodynamics 32, 597-616.
Ortmaier, H., 1989. Studies of dielectric properties of methanolic electrolyte solutions in the range of 0.9 to 89 GHz . Diploma Thesis. University of Regensburg, Regensburg (Germany).
Papaiconomou, N., Simonin, J.P., Bernard, O., Kunz, W., 2002. MSA-NRTL model for the description of the thermodynamic properties of electrolyte solutions. Physical Chemistry Chemical Physics 4, 4435-4443.
Pasztor, A.J., Criss, C.M., 1978. Apparent molal volumes and heat-capacities of some 1-1 electrolytes in anhydrous methanol at 25-degrees-C. Journal of Solution Chemistry 7, 27-44.
Robinson, R.A., Stokes, R.H., 1970. Electrolyte Solutions, 2nd ed. Butterworth, London.
Salimi, H.R., Taghikhani, V., Ghotbi, C., 2005. Application of the GV-MSA model to the electrolyte solutions containing mixed salts and mixed solvents. Fluid Phase Equilibria 231, 67-76.
Simonin, J.P., Krebs, S., Kunz, W., 2006. Inclusion of ionic hydration and association in the MSA-NRTL model for a description of the thermodynamic properties of aqueous ionic solutions: Application to solutions of associating acids. Industrial \& Engineering Chemistry Research 45, 4345-4354.
Takenaka, N., Takemura, T., Sakurai, M., 1994. Partial molar volumes of uni-univalent electrolytes in methanol plus water.1. lithium-chloride, sodium-chloride, and potassium-chloride. Journal of Chemical and Engineering Data 39, 207-213.
Tan, S.P., Adidharma, H., Radosz, M., 2005. Statistical associating fluid theory coupled with restricted primitive model to represent aqueous strong electrolytes. Industrial \& Engineering Chemistry Research 44, 4442-4452.
Taniewska-Osinska, S., Chadzynski, P., 1976. Viscosity of NaI solutions in mixtures of water with methanol, ethanol and n-propanol. Acta Universitatis Lodziensis SeriesActa Univ. Lodz. Ser. II 6, 37-46.
Tomasula, P., Czerwienski, G.J., Tassios, D., 1985. Vapor pressures and osmotic coefficients: electrolyte solutions of methanol. Fluid Phase Equilibria 38, 129-135.
Triolo, R., Grigera, J.R., Blum, L., 1976. Simple electrolytes in mean spherical approximation. Journal of Physical Chemistry 80, 1858-1861.
Wagner, W., Pruss, A., 2002. The IAPWS formulation 1995 for the thermodynamic properties of ordinary water substance for general and scientific use. Journal of Physical and Chemical Reference Data 31, 387-535.
Waisman, E., Lebowitz, J.L., 1970. Exact solution of an integral equation for structure of a primitive model of electrolytes. Journal of Chemical Physics 52, 4307-4311.
Wang, S., Song, Y., Chen, C.-C., 2011. Extension of COSMO-SAC solvation model for electrolytes. Industrial \& Engineering Chemistry Research 50, 176-187.
Wohlfarth, C., 1991. Landolt-Börnstein, Group IV Volume 6: Static Dielectric Constant of Pure and Binary Liquid Mixtures. Springer, Berlin 221.
Wolbach, J.P., Sandler, S.I., 1998. Using molecular orbital calculations to describe the phase behavior of cross-associating mixtures. Industrial \& Engineering Chemistry Research 37, 2917-2928.
Yan, W.D., Xu, Y.J., Han, S.J., 1994. Activity coefficients of sodium chloride in methanol-water mixed solvents at 298.15 K. Acta Chimica Sinica 52, 937-946.
Yang, M.F., Leng, C.L., Li, S.C., Sun, R.Y., 2007. Activity coefficients of potassium iodide in methanol-water mixed solvents by electromotive force measurements. Chemical Research 18, 49-53.
Ye, S., Xans, P., Lagourette, B., 1994. Modification of the Pitzer model to calculate the mean activity-coefficients of electrolytes in a water-alcohol mixed-solvent solution. Journal of Solution Chemistry 23, 1301-1315.
Zafarani-Moattar, M.T., Nasirzade, K., 1998. Osmotic coefficient of methanol plus LiCl , plus LiBr , and plus LiCH 3 COO at 25 degrees C . Journal of Chemical and Engineering Data 43, 215-219.
Zhang, H.L., Han, S.J., 1996. Viscosity and density of water plus sodium chloride plus potassium chloride solutions at 298.15 K . Journal of Chemical and Engineering Data 41, 516-520.
Zuo, J.Y.X., Zhang, D., Fürst, W., 2000. Extension of the electrolyte EOS of Furst and Renon to mixed solvent electrolyte systems. Fluid Phase Equilibria 175, 285-310.


[^0]:    *Corresponding author. Tel.: +49 2317552635.
    E-mail address: g.sadowski@bci.tu-dortmund.de (G. Sadowski).

[^1]:    ${ }^{\mathrm{a}}$ Due to the low solubility of potassium salts in ethanol, the potassium parameters could not be adjusted to binary ethanol/salt solutions. Instead, they were estimated from ternary MIACs with the constraint to fit in the homologous parameter series of the cations (see text).

[^2]:    ${ }^{\mathrm{a}}$ This work.

