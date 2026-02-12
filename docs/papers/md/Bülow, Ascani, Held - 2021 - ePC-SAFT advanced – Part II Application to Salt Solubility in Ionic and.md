# ePC-SAFT advanced - Part II: Application to Salt Solubility in Ionic and Organic Solvents and the Impact of Ion Pairing 

Mark Bülow, Moreno Ascani, Christoph Held*<br>Laboratory of Thermodynamics, Department of Biochemical and Chemical Engineering, TU Dortmund, Emil-Figge Str. 70, 44277 Dortmund, Germany

## ARTICLE INFO

## Article history:

Received 5 January 2021
Revised 6 February 2021
Accepted 22 February 2021
Available online 2 March 2021

## KEYWORDS:

electrolyte thermodynamics
solubility product
ionic liquids
predictions
Bjerrum
Born


#### Abstract

The applications of electrolyte thermodynamic models to non-aqueous systems is of great value to reduce experimental effort and gain inside into molecular interactions. A large-scale application is for example the design of advanced battery electrolytes. For non-aqueous electrolyte systems, the Born term was found to be important, as it accounts for the transfer of ions from water into non-aqueous medium. In part one of this study [Bülow et al., Fluid Phase Equilibria 2021, 112967] the Born term was combined with a concentration-dependent dielectric constant within the ePC-SAFT framework (electrolyte Perturbed-Chain Statistical Associating Fluid Theory). In the present work, the Bjerrum treatment for ion pairing was included in the Debye-Hückel framework within ePC-SAFT. The approach was validated by experimental data for the dissociation of salts in organic solvents derived from conductivity measurements. Further, solubility was modeled of alkali halides in organic solvents and in ionic liquids. Modeling solubility required access to the solubility product $K_{S P}$, which does not depend on the solvent. The approach within this work was to first determine $K_{S P}$ using experimental solubility data in water and the respective ePC-SAFT predicted activity coefficients prior to predict activity coefficients in non-aqueous medium, finally yielding solubility. The so-determined solubility values were found to be in reasonable agreement with the experimental data without fitting model parameters to any data of the non-aqueous solutions. The solubility product requires the solid form of the precipitating salt to be equal for all solvents; as alkali salts precipitate from aqueous solutions as hydrates, the method cannot be applied. Therefore, a methodology is presented to extrapolate the high-temperature $K_{S P}$ of anhydrates to lower temperature. Using the so-extrapolated $K_{S P}$ allowed predicting solubility of non-solvates in other solvents.


© 2021 Elsevier B.V. All rights reserved.

## 1. Introduction

Besides use in oil and gas industries, electrolyte thermodynamic models such as present in ref. [1] potentially help reducing the experimental effort to design efficient batteries. Batteries support the changeover from fossil fuels that is the mobility turnaround and storage of energy from regenerative sources. One crucial factor for the battery efficiency is the ionic conductivity of the electrolyte medium. The conductivity depends on the electrolyte medium and the amount and kind of dissolved salts [2]. Classical electrolyte medium for batteries are of organic nature (e.g. carbonates); it provides a mean dielectric constant to limit ion pairing and allows high ionic mobility. Ionic liquids (ILs), consisting of ions in liquid state bringing intrinsic ionic conductivity, have the potential to substitute or dope the classical electrolyte medium. Research on ILs is nowadays a full-grown scientific tree, with applications in

[^0]chemical industry, downstream processing, and biotechnology. The interesting and unique properties of ILs e.g., the almost negligible vapor pressure, may be capitalized on for various applications, including batteries. The solubility of potential salts in an electrolyte medium dictates the success of new battery design [3]. Thus, understanding the molecular interactions governing the salt solubility facilitates the selection of promising salt-electrolyte pairs.

Facing the challenges for non-aqueous electrolyte solutions, electrolyte Perturbed-Chain Statistical Associating Fluid Theory (ePC-SAFT) was advanced with the Born term for ion solvation combined with a concentration-dependent dielectric constant, defined as altered Born contribution, c.f. part one of this research study [1]. The altered Born contribution describes the required Gibbs energy of solvation for the transfer of ions from vacuum into a respective solvent. With this, the advanced ePC-SAFT model considers the change from water to non-aqueous solvents by the drastic decrease of the dielectric constant of the solvent compared to water as well as the influence of ion concentration on the dielectric constant in the system. The advanced approach within the ePC-SAFT framework was intensively investigated by
predicting infinite-dilution properties as well as mean ionic activity coefficients (MIACs) up to high concentrations. The success of the MIAC prediction results was found to be caused by accounting for the salt-concentration dependence of the dielectric constant in the Debye-Hückel theory and the altered Born contribution. MIACs in water have been used for model validations similar to the ePC-SAFT advanced approach and the concentration dependence of the dielectric constant for ion-ion and ion-water interactions was addressed in works on aqueous electrolyte systems.[4-6] It is important to note that the Born term was implemented to other electrolyte EoS without considering the variation of the dielectric constant with composition, inter alia the EoS Peng-Robinson[7,8], Soave-Redlich-Kwong,[9] electrolyte Cubic-Plus-Association (eCPA)[9-11], ePPC-SAFT,[12,13] and SAFTVRE (SAFT variable range for electrolytes) [14]. Still, the use of the Born term upon model development is often discarded as its contribution to aqueous systems is negligible. Indeed, the impact of the Born term on excess properties is effectively zero as long as the dielectric constant is treated independent of composition. Consequently, ePC-SAFT advanced reduces to the original model ePCSAFT revised from $2014{ }^{15}$ if the concentration dependence of the dielectric constant is omitted.

This work additionally accounts for physico-chemical effects for ionic species when dissolved in medium of low dielectric constant. In aqueous solutions, the degree of dissociation for monovalent salts is usually high for conventional alkali halides, with degrees of dissociation commonly greater than 0.8 over a broad concentration range [16]. The assumption of complete dissociation for such alkali halides in aqueous systems is regarded eligible in the literature; assuming complete dissociation for such strong electrolytes has proven to give reasonable modeling results for properties of the electrolyte solutions.[17-21] The influence of minor species (ion pairs) on modeling results is usually compensated indirectly by the model parameters. Besides some exceptions form this rule (e.g., alkali fluorides), the formation of ion pairs and higher clusters is favored in solvents with a dielectric constant lower than that of water. The degree of dissociation of salts in such solvents e.g., alcohols might be lower than 0.2 even at low to moderate ion concentrations. Bjerrum provided an implicit theory based on the electrostatic potential for the formation of ion pairs, yielding a critical distance of counter ions [22]. The distance is the Bjerrum length $l_{B}$ that allows a differentiation to whether or not free ions will pair at a certain concentration [23]. The existence of strongly bounded, long-lived ion pairs is confirmed by spectroscopic and dielectric measurements [24,25]. It is nowadays widely accepted that ion pairs can exist not only as contact ion pairs, where cation and anions directly contact, but also as solvent-shared and solventseparated ion pairs [23]. The application of Bjerrum treatment has shown to be beneficial for the modeling of salt solubility in various solvents at low permittivity at 298 K with COSMO-RS-ES [26].

The present manuscript presents for the first time the successful application of an electrolyte EoS to predict salt solubility in non-aqueous systems. In a first step of this work, results for the ion-pairing behavior predicted with Bjerrum treatment are validated against experimental data on dissociation degrees, which have been derived from conductivity measurements. The approach was applied to electrolyte systems containing salt in organic solvents and salt in ionic liquids. Then, the impact of the Bjerrum treatment on solubility predictions with ePC-SAFT advanced was tested against experimental data of alkali halides in organic solvents. Predictions with the new methodology within the ePC-SAFT framework including the concentration-dependent dielectric constant[27] in the Debye-Hückel theory and the altered Born contribution[1] allow a qualitative representation of the solubility of monovalent alkali-halide salts in organic solvents and in ionic liquids. Modeling salt solubility was approached by a treatment of
pharmaceuticals with thermodynamic models [28]. This method allows reasonable predictions of salt solubility in organic solvents and ionic liquids.

## 2. Thermodynamic approaches: ePC-SAFT, Born theory, Bjerrum treatment and solubility product

## 2.1. ePC-SAFT advanced

ePC-SAFT advanced calculates the residual Helmholtz energy from five contributions now including the altered Born contribution as addition to the methodology of ePC-SAFT revised in 2014 [15] (Equation (1)).

$$
\begin{equation*}
a^{r e s}=a^{h c}+a^{d i s p}+a^{a s s o c}+a^{D H}(\alpha, \varepsilon(\bar{x}))+a^{B o r n}(\varepsilon(\bar{x})) \tag{1}
\end{equation*}
$$

Five pure-component parameters are assigned to associating compounds like water. The pure-component parameters are the segment number $m_{i}^{\text {seg }}$, the diameter $\sigma_{i}$, the dispersion energy $u_{i}$, and additionally the association energy $\varepsilon^{A i B i}$ and the association volume $\kappa^{A i B i}$. Only two pure-ion parameters characterize inorganic ions, $\sigma_{i}$ and $u_{i}$.

### 2.2. Modeling activity coefficients

The fugacity coefficient is derived from the residual Helmholtz energy by differentiation (Equation (2)).

$$
\begin{equation*}
\ln \left(\varphi_{\mathrm{i}}\right)=\frac{\mu_{\mathrm{i}}^{\mathrm{res}}}{\mathrm{k}_{\mathrm{B}} \cdot \mathrm{~T}}-\ln \left(1+\left(\frac{\partial\left(\frac{\mathrm{a}^{\mathrm{res}}}{\mathrm{k}_{\mathrm{B}} \cdot \mathrm{~T}}\right)}{\partial \rho}\right)\right) \tag{2}
\end{equation*}
$$

Modeling mixtures requires the combining rules of BerthelotLorentz Equations (3) and ((4)).

$$
\begin{equation*}
\sigma_{i j}=\frac{1}{2}\left(\sigma_{i}+\sigma_{j}\right) \tag{3}
\end{equation*}
$$

$$
\begin{equation*}
u_{i j}=\sqrt{u_{i} u_{j}}\left(1-k_{i j}\right) \tag{4}
\end{equation*}
$$

Equation (4) introduces the binary interaction parameter $\mathrm{k}_{\mathrm{ij}}$ that can be used to alter the dispersion energy $u_{i j}$ in the mixture. In this work, values were assigned for the pairs ion-ion and ion-solvent; these and the pure-component parameters of the ions were regressed previously [15]. For the prediction of salt solubility in the organic solvent or the IL, $k_{i j}$-values between ion and IL-ions or ion and organic solvent were set to zero. The generic activity coefficient of each ion $i$ in the mixture was calculated from the ratio of the fugacity coefficients in mixture to the pure-ion state.

$$
\begin{equation*}
\gamma_{i}=\frac{\varphi_{\mathrm{i}}(\mathrm{~T}, \mathrm{p}, \overrightarrow{\mathrm{x}})}{\varphi_{0 \mathrm{i}}\left(\mathrm{~T}, \mathrm{p}, \mathrm{x}_{i}=1\right)} \tag{5}
\end{equation*}
$$

The activity coefficient at infinite dilution $\gamma_{i}^{*}$ of each ion $i$ in the mixture was calculated from the ratio of the fugacity coefficients in the mixture and at infinite dilution according to Equation (6).

$$
\begin{equation*}
\gamma_{i}^{*}=\frac{\varphi_{\mathrm{i}}(\mathrm{~T}, \mathrm{p}, \overrightarrow{\mathrm{x}})}{\varphi_{0 \mathrm{i}}\left(\mathrm{~T}, \mathrm{p}, \mathrm{x}_{i} \rightarrow 0\right)} \tag{6}
\end{equation*}
$$

Salts were considered to be fully dissociated into cations $c$ and anions $a$ with the respective stoichiometric coefficients $v_{c}$ and $v_{a}$. Any ion-averaged activity coefficient was calculated from the ionbased activity coefficients of the cation and anion:

$$
\begin{equation*}
\gamma_{ \pm}^{*}=\left(\gamma_{c}^{*, v c} \cdot \gamma_{a}^{*, v a}\right)^{\frac{1}{v_{c}+v_{a}}} \tag{7}
\end{equation*}
$$

$$
\begin{equation*}
\gamma_{ \pm}=\left(\gamma_{c}^{\nu c} \cdot \gamma_{a}^{\nu a}\right)^{\frac{1}{v_{c}+v_{a}}} \tag{8}
\end{equation*}
$$

The MIAC $\gamma_{ \pm}^{*}$ and the generic gMIAC $\gamma_{ \pm}$are derived by inserting the respective activity coefficients into Equation (7) and (8).

### 2.3. Ion pairing and influence on $a^{D H}$

The Bjerrum treatment was developed to compensate for simplifications introduced by the linearization of the PoissonBoltzmann equation in deriving the Debye-Hückel theory. The assumption that the electrostatic energy $z_{j} e \psi_{i}(r)$ between any two ions $i$ and $j$ at distance $r$ between both ion centers is much lower than the mean thermal energy per degree of freedom $k_{B} T$ is not valid for small $r$.[29] In the Bjerrum treatment within the DebyeHückel framework, ionic species are defined as either ion pairs or free ions in the system based on the law of mass action. For a 1:1 electrolyte the law of mass action is given by the dissociation constant $K_{i p}$, the activity ratio between paired and unpaired ions (Equation (9)),

$$
\begin{equation*}
K_{i p}(T)=\frac{x_{i p} \cdot \gamma_{i p}^{*}}{x_{f} \cdot\left(\gamma_{ \pm}^{*}\left(x_{f}(\alpha)\right)^{2}\right.}=\frac{(1-\alpha) \cdot \gamma_{i p}^{*}}{\alpha^{2} \cdot x_{ \pm} \cdot\left(\gamma_{ \pm}^{*}\left(x_{f}(\alpha)\right)^{2}\right.} \tag{9}
\end{equation*}
$$

where $\alpha$ is the dissociation degree, $\gamma_{i p}^{*}$ the activity coefficient of the ion pair, $\gamma_{ \pm}^{*}\left(x_{f}(\alpha)\right)$ the MIAC at the mole fraction of the free ions $x_{f}$ (depending on the degree of dissociation) and $x_{ \pm}$is the sum of $x_{f}$ and $x_{i p}$. Within the Bjerrum treatment, the activity coefficient $\gamma_{i p}^{*}$ is taken as unity. Bjerrum derived an expression for $K_{i p}$ by integrating the probability function of finding a counter ion $j$ at distance $r$ from ion $i$ between $r=a$ and $r=l_{B}$. In first approximation, the closest approach distance $a$ between both ions is taken as the sum of both ionic radii. $l_{B}$ is the Bjerrum length (Equation (10)),

$$
\begin{equation*}
l_{B}=\frac{\left|z_{i} z_{j}\right| e^{2}}{8 \pi \varepsilon_{0} \varepsilon_{r} k_{B} T} \tag{10}
\end{equation*}
$$

defined by the ratio of the valence $z$ and the thermal energy $k_{B} T$ times the dielectric constant; it dictates the cut-off for ionpair building from electrostatic potentials. The overall mass law is then described with an expression that contains temperature, density and the apparent dielectric constant of the medium in Equation (11),

$$
\begin{equation*}
K_{i p}(T)=4 \pi \rho_{N} \int_{a}^{l_{B}} \exp \left(\frac{\left|z_{i} z_{j}\right| e^{2}}{4 \pi \varepsilon_{0} \varepsilon_{r} k_{B} T} \cdot \frac{1}{r}\right) r^{2} d r \tag{11}
\end{equation*}
$$

with $\rho_{N}$ as the number density. The equilibrium constant can be understood as the partition function of two generic ions $i$ and $j$ between the spatial domains for the ion pairs ( $a \leq r \leq l_{B}$ ) and for the individual free ion ( $l_{B} \leq r \leq \infty$ ). This holds at low ion concentration. Fig. 1 graphically explains the Bjerrum threshold and the distinction between free and paired ions from the radial distribution function.

The Bjerrum treatment was developed within the Debye-Hückel framework; the approximations behind both theories are equal. Only interionic electrostatic interactions with other free ions are considered for calculating the chemical potential of the free ions (Equation (12)).

$$
\begin{equation*}
\mu_{i}^{D H}\left(x_{f}\right)=\left(\frac{\partial A^{D H}}{\partial \rho_{i}\left(x_{f}\right)}\right)_{T, V, N_{j \neq i}}=-\frac{e^{2} z_{i}^{2} \kappa}{24 \pi \varepsilon_{0} \varepsilon_{r}}\left[2 \chi_{i}+\frac{\sum_{k} x_{k, f} z_{k}^{2} \sigma_{k}}{\sum_{k} x_{k, f} z_{k}^{2}}\right] \tag{12}
\end{equation*}
$$

The number of free ions $\rho_{i}\left(x_{f}\right)$ and $x_{f}$ is calculated from the related total amount of ion of type $i, \rho_{i}$ and $x_{i}$, by multiplying with the dissociation degree $\alpha$. The resulting activity coefficient of the free ion $i$ is given as in Equation (13).

$$
\begin{equation*}
\ln \gamma_{i}^{*}\left(x_{f}\right)=-\frac{e^{2} z_{i}^{2} \kappa}{24 \pi \varepsilon_{0} \varepsilon_{r}}\left[2 \chi_{i}+\frac{\sum_{k} x_{k, f} z_{k}^{2} \sigma_{k}}{\sum_{k} x_{k, f} z_{k}^{2}}\right] \tag{13}
\end{equation*}
$$

![](https://cdn.mathpix.com/cropped/225b84f4-84d1-4768-92e8-fb2e12773952-03.jpg?height=585&width=841&top_left_y=193&top_left_x=1079)
Fig. 1. Graphical representation of the distinction between free ions and ion pairs. The yellow area depicts the domain between the closest distance of two counter ions and the threshold that is the Bjerrum length. In this area, ion pairing occurs, while for the crosshatched domain free ions are present.

The activity coefficient of the ion pair is set to unity i.e., there is no Coulomb interactions induced by ion pairs with each other or with other ions; further, the dielectric constant of the system is not influenced by the formation of ion pairs at cost of free ions. In Equation (14), the dissociation degree is rewritten as function of $K_{i p}$ and $\gamma_{ \pm}^{*}$; Equation (14) was evaluated numerically using a bisection algorithm due to its implicit character.

$$
\begin{equation*}
\alpha=\frac{-1+\sqrt{1+4 x_{ \pm} K_{i p} \frac{\left(\gamma_{ \pm}^{*}\left(x_{f}\right)\right)^{2}}{\gamma_{i p}^{*}}}}{2 x_{ \pm} K_{i p} \frac{\left(\gamma_{ \pm}^{*}\left(x_{f}\right)\right)^{2}}{\gamma_{i p}^{*}}} \tag{14}
\end{equation*}
$$

The Bjerrum treatment was implemented to provide an implicit dissociation degree for monovalent salts in various media. From the dissociation degree, the thermodynamic properties of the components are re-evaluated to account for ion pairing in the ePC-SAFT framework. The implementation of Bjerrum treatment to ePC-SAFT advanced requires a recalculation of the Debye-Hückel theory to account for ion pairs and the concentration-dependent dielectric constant. In the current version of ePC-SAFT advanced, a mixing rule linear in mole fraction is used for calculating the dielectric constant of the system at the actual composition. Furthermore, the composition derivatives of the dielectric constant are accounted for in the calculation of the fugacity coefficients of all the ionic components. The procedure is explained in two previous works.[27,30] The incorporation of the dissociation degree into the routine requires reformulating the contribution of the Debye-Hückel theory to the free energy of the system at given conditions ( $T, V, \bar{x}$ ) with a degree of dissociation given by $\alpha_{i}$ and reapplying the procedure again. In the Debye-Hückel theory only free ions contribute to interionic electrostatic interactions. The total number density $\rho_{i}$ of ion $i$ to the number density of free ions $\rho_{i}^{f i}$ needs correction using the dissociation degree $\alpha_{i}$ in Equation (15).

$$
\begin{equation*}
\alpha_{i}=\frac{\rho_{i, f}}{\rho_{i}} \tag{15}
\end{equation*}
$$

Furthermore, the distance of closest approach $a_{i}$ was replaced with the Bjerrum length $l_{B}$. The distance of closest approach $R_{i}$ was approached according to Equation (16).

$$
R_{i}=\left\{\begin{array}{l}
a_{i} \text { if } a_{i}>l_{B}  \tag{16}\\
l_{B} \text { if } a_{i}<l_{B}
\end{array}\right.
$$

Inserting the dissociation degree and distance of closest approach into the Debye-Hückel contribution $A^{D H}$, the inverse

Debye screening length $\kappa$ and the auxiliary function $\chi_{i}$ in Equations (17) to (19) were derived.

$$
\begin{equation*}
A^{D H}=-\frac{1}{12 \pi \varepsilon_{0} \varepsilon_{r}} \sum_{j} \alpha_{j} \rho_{j} q_{j}^{2} \kappa \chi_{j}=-\frac{N_{A} e^{2}}{12 \pi \varepsilon_{0} \varepsilon_{r}} \sum_{j} \alpha_{j} c_{j} z_{j}^{2} \kappa \chi_{j} \tag{17}
\end{equation*}
$$

$$
\begin{equation*}
\kappa=\sqrt{\frac{\rho_{N} e^{2}}{k_{B} T \varepsilon_{0} \varepsilon_{r}} \sum_{j} \alpha_{j} x_{j} z_{j}^{2}} \tag{18}
\end{equation*}
$$

$$
\begin{equation*}
\chi_{i}=\frac{3}{\left(\kappa R_{i}\right)^{3}}\left[\frac{3}{2}+\ln \left(1+\kappa R_{i}\right)-2\left(1+\kappa R_{i}\right)+\frac{1}{2}\left(1+\kappa R_{i}\right)^{2}\right] \tag{19}
\end{equation*}
$$

For evaluating the fugacity coefficient of the ion $i$, the Debye-Hückel contribution to the reduced molar free energy $a^{D H}$ (Equation (20)) as well as its partial derivative $\left(\partial a^{D H} / \partial x_{i}\right)_{T, v, x_{j \neq i}}$ (Equation (21)) with respect to the molar fraction of ion $i$ need to be calculated.

$$
\begin{align*}
a^{D H}=\frac{A^{D H}}{k_{B} T \rho_{N}}= & -\frac{e^{2}}{12 \pi \varepsilon_{0} \varepsilon_{r} k_{B} T} \sum_{i} \alpha_{i} x_{i} z_{i}^{2} \kappa \chi_{i}  \tag{20}\\
\left(\frac{\partial a^{D H}}{\partial x_{i}}\right)_{T, v, x_{j \neq i}}= & -\frac{e^{2}}{12 \pi k_{B} T}\left[\frac{\partial \kappa}{\partial x_{i}} \frac{1}{\varepsilon_{0} \varepsilon_{r}} \sum_{j} \alpha_{j} x_{j} z_{j}^{2} \kappa \chi_{j}\right. \\
& -\frac{\kappa}{\left(\varepsilon_{0} \varepsilon_{r}\right)^{2}} \frac{\partial\left(\varepsilon_{0} \varepsilon_{r}\right)}{\partial x_{i}} \sum_{j} \alpha_{j} x_{j} z_{j}^{2} \kappa \chi_{j}+\frac{\kappa}{\varepsilon_{0} \varepsilon_{r}} \alpha_{i} z_{i}^{2} \chi_{i} \\
& \left.+\frac{\kappa}{\varepsilon_{0} \varepsilon_{r}} \sum_{j} \alpha_{j} x_{j} z_{j}^{2} \frac{\partial \chi_{j}}{\partial x_{i}}\right] \tag{21}
\end{align*}
$$

The partial derivatives of the inverse screening length $\kappa$ and the auxiliary function $\chi_{i}$ are given in Equation (22) and (23).

$$
\begin{align*}
& \left(\frac{\partial \kappa}{\partial x_{i}}\right)_{T, v, x_{j \neq i}}=\frac{1}{2}\left(\frac{\rho_{N} e^{2}}{k_{B} T \varepsilon_{0} \varepsilon_{r}} \sum_{j} \alpha_{j} x_{j} z_{j}^{2}\right)^{-\frac{1}{2}}  \tag{22}\\
& \left\{-\frac{\rho_{N} e^{2}}{k_{B} T\left(\varepsilon_{0} \varepsilon_{r}\right)^{2}} \frac{\partial\left(\varepsilon_{0} \varepsilon_{r}\right)}{\partial x_{i}} \sum_{j} \alpha_{j} x_{j} z_{j}^{2}+\frac{\rho_{N} e^{2}}{k_{B} T \varepsilon_{0} \varepsilon_{r}} \alpha_{i} z_{i}^{2}\right\} \\
& =\left(\frac{\rho_{N} e^{2}}{k_{B} T}\right)^{\frac{1}{2}}\left\{-\frac{1}{2}\left(\varepsilon_{0} \varepsilon_{r}\right)^{-\frac{3}{2}} \frac{\partial\left(\varepsilon_{0} \varepsilon_{r}\right)}{\partial x_{i}}\left[\sum_{j} \alpha_{j} x_{j} z_{j}^{2}\right]^{\frac{1}{2}}\right. \\
& \left.+\frac{1}{2 \sqrt{\varepsilon_{0} \varepsilon_{r}}}\left[\sum_{j} \alpha_{j} x_{j} z_{j}^{2}\right]^{-\frac{1}{2}} \alpha_{i} z_{i}^{2}\right\}  \tag{22}\\
& =-\frac{9}{\left(\kappa R_{i}\right)^{4}}\left[\frac{3}{2}+\ln \left(1+\kappa R_{i}\right)-2\left(1+\kappa R_{i}\right)\right. \\
& \left.\left(\frac{\partial \chi_{i}}{\partial x_{i}}\right)_{T, v, x_{j \neq i}}+\frac{1}{2}\left(1+\kappa R_{i}\right)^{2}\right] R_{i} \frac{\partial \kappa}{\partial x_{i}}+\frac{3}{\left(\kappa R_{i}\right)^{3}}\left[\frac{1}{1+\kappa R_{i}}-1+\kappa R_{i}\right] \\
& R_{i} \frac{\partial \kappa}{\partial x_{i}}=3 \frac{\partial \kappa}{\partial x_{i}}\left\{-\frac{\chi_{i}}{\kappa}+\frac{R_{i}}{\left(\kappa R_{i}\right)^{3}}\left[\frac{1}{1+\kappa R_{i}}-1+\kappa R_{i}\right]\right\} \tag{23}
\end{align*}
$$

It should be noted that the derivative of the dissociation degree with respect to the molar fraction was neglected within this work.

### 2.4. Ion solvation and $a^{\text {Born }}$

The altered Born contribution was derived in part one of this study;[1] it describes the transfer work of an initially charged
sphere from infinite distance to a first neutral species in the solvent that is afterwards recharged. The respective molar Helmholtz free energy contribution to ePC-SAFT is shown in Equation (24).

$$
\begin{equation*}
a^{\text {Born }}=\frac{A^{\text {Born }}}{k T N}=-\frac{e^{2}}{4 \pi \varepsilon_{0} k_{B} T}\left(1-\frac{1}{\varepsilon_{r}}\right) \sum_{i} \frac{x_{i} z_{i}^{2}}{a_{i}} \tag{24}
\end{equation*}
$$

The distance of closest approach $a_{i}$ is here set to the respective diameter $\sigma_{\text {ion }}$ i.e., to the original ePC-SAFT parameter. The partial derivative of the reduced molar Born contribution is derived in accordance with the previous work with a concentration-dependent dielectric constant, presented in Equation (25).

$$
\begin{align*}
\left(\frac{\partial a^{\text {Born }}}{\partial x_{i}}\right)_{T, v_{N}, x_{j \neq i}}= & -\frac{e^{2}}{4 \pi k_{B} T \varepsilon_{0}}\left(1-\frac{1}{\varepsilon_{r}}\right) \frac{z_{i}^{2}}{a_{i}} \\
& -\frac{e^{2}}{4 \pi k_{B} T \varepsilon_{0}}\left(\frac{1}{\varepsilon_{r}^{2}}\right)\left(\frac{\partial \varepsilon_{r}}{\partial x_{i}}\right)_{T, v, x_{j \neq i}} \tag{25}
\end{align*}
$$

The Gibbs energy of transfer is in direct connection with the successful implementation of the altered Born contribution and was used to validate the approach in part one of this research [1]. $\Delta^{T r} G_{i}$ is calculated from the ratio of the fugacity coefficients at infinite dilution in solvents $S 1$ and $S 2$.

$$
\begin{equation*}
\Delta^{T r} G_{i}=R T \cdot \ln \left(\frac{\varphi_{i}^{\infty, S 2}}{\varphi_{i}^{\infty, S 1}}\right) \tag{26}
\end{equation*}
$$

### 2.5. Framework for solubility prediction

The solubility of a component $i$ may be calculated from standard state properties as depicted in Equation (27) with $T_{i}^{T_{p}}$, $\Delta h_{i}^{S L}\left(T_{i}^{T_{p}}\right)$ and $\Delta c_{p, i}^{S L}$ being, respectively, the triple point temperature, the melting enthalpy of pure component $i$ at the triple point and the difference between the molar heat capacities of pure component $i$ in the solid and liquid state, $c_{p, i}^{0 L}-c_{p, i}^{0 S}$.

$$
\begin{equation*}
\ln x_{i}^{L} \gamma_{i}^{L}=\frac{\Delta h_{i}^{S L}\left(T_{i}^{T_{p}}\right)}{R T}\left(1-\frac{T_{i}^{T_{p}}}{T}\right)-\frac{\Delta c_{p, i}^{S L}}{R T}\left(T-T_{i}^{T_{p}}\right)+\frac{\Delta c_{p, i}^{S L}}{R} \ln \frac{T}{T_{i}^{T_{p}}} \tag{27}
\end{equation*}
$$

Equation (27) is valid under the assumption that the existing solid phase is an ideal pure crystal and the molar heat capacities in the solid and liquid state are independent of temperature. For alkali-halide salts, the standard state properties are scarce and far out of the region of interest for the prediction of the salt solubility in organic solvents and ILs. It would therefore require an extrapolation over a large temperature range.

The right-hand side of Equation (27) is a function of the system pressure p and temperature T, if the solid phase is assumed to consist of the pure salt without the formation of hydrates. For a single salt, Equation (27) can be rearranged using the solubility product $\ln K_{S P}$ (Equation (28)),

$$
\begin{equation*}
\ln K_{S P}=\sum_{i} v_{i} \ln x_{i} \gamma_{i} \tag{28}
\end{equation*}
$$

where $v_{i}$ is the stochiometric coefficient of ion $i$ forming the salt and the summation is over all ions of the specific salt. For monovalent alkali halides, the equation reduces to the summation of salt cation and anion. The solubility product is calculated once in a reference solvent; this value is independent of any solvent and only depends on temperature and pressure. The reference solvent in this work is water, allowing the calculation of exact $K_{S P}$ for all investigated salts over a broad temperature range. The applied ePC-SAFT pure-component parameters have been regressed to inter alia MIAC

![](https://cdn.mathpix.com/cropped/225b84f4-84d1-4768-92e8-fb2e12773952-05.jpg?height=579&width=665&top_left_y=193&top_left_x=234)
Fig. 2. Visualized calculation of the SLE cycle to obtain solubility of alkali-halide salts in organic or ionic solvents using the solubility product $K_{S P}$. Solubility is predicted through a hypothetical pure-liquid state of the salt (2).

data in water at 298 K . Thus, the starting point and $K_{S P}$ are considered to be as accurate as possible.

The solubility of salts in the organic solvents or ILs is afterwards predicted by minimizing the objective function (Equation (29)),

$$
\begin{equation*}
O F=\left(K_{S P}-\sum_{i} v_{i} \ln x_{i} \gamma_{ \pm}\right)^{2} \tag{29}
\end{equation*}
$$

where $x_{i}$ is the salt solubility and $\gamma_{ \pm}$is the generic gMIAC of the salt $i$. Fig. 2 visualizes the applied approach utilizing $K_{S P}$ to predict the solubility of alkali-halide salts in organic solvents and ILs. The cycle shows the starting point of dissolved salts in a solvent (water, 1), the transition via the hypothetical fluid pure salt (2) to the solubility in the second solvent (3).

### 2.6. Use of $K_{S P}$ and the hypothetical solubility of metastable anhydrates as access to solubility in non-aqueous systems

The application of the solubility product $K_{S P}$ for the prediction of solubility in different solvents is only allowed if the salt forms the same solid form in the respective solvents. A hydrate, in contrast to the pure salt, is a solid that includes $n$ molecules of water and occurs regularly in aqueous solutions in contrast to other solvents where the pure salt forms. Especially lithium salts form mono-, di, or higher hydrates at room temperature, while the precipitating solid in non-aqueous solvents can impossibly be hydrate forms. Lithium salts form stable anhydrates in water only at $\mathrm{T}>373 \mathrm{~K}, \mathrm{NaBr}$ at $\mathrm{T}>323 \mathrm{~K}$. To apply the $K_{S P}$ method to predict the solubility of these salts in non-aqueous systems at room temperature, a new methodology is presented here. In this method, the $K_{S P}$ values are accessed assuming that the non-solvate would theoretically form at 298 K ; this is possible by extrapolating the $K_{S P}$-values to lower temperatures. Alkali-halide salts that form hydrates in the investigated temperature range are per conventions not available for the application of $K_{S P}$. The van 't Hoff equation is the connection of equilibrium constants $K$, in this work the solubility product $K_{S P}$, and the temperature (Equation (30)).

$$
\begin{equation*}
\left(\frac{\partial \ln K_{S P}}{\partial T}\right)_{p}=\frac{\Delta h^{S L}}{R T^{2}} \Leftrightarrow K_{S P}(T)=K_{S P}\left(T_{0}\right) \cdot \exp \left(\frac{\Delta h^{S L}}{R T}\right) \tag{30}
\end{equation*}
$$

Within a certain temperature range, the van 't Hoff equation allows to extrapolate $K_{S P}$ if $c_{p}$ effects are neglected. Therefore, $K_{S P}$ was calculated with ePC-SAFT advanced for at least two temperatures where the anhydrate precipitates from water. Afterwards, $K_{S P}$

Table 1
Pure-component parameters of organic solvents used in this work. All compounds have a 2B association scheme.
| Solvent | $m_{i}^{\text {seg }}$ | $\sigma_{i} / \AA$ | $u_{i} / k_{B} / \mathrm{K}$ | $\varepsilon^{A i B i} / \mathrm{K}$ | $\kappa^{A i B i}$ | Ref. |
| :--- | :--- | :--- | :--- | :--- | :--- | :--- |
| Water | 1.2047 | 2.7927 | 353.95 | 2425.7 | 0.04509 | $[32]$ |
| Methanol | 1.5255 | 3.2300 | 188.90 | 2899.5 | 0.03518 | $[33]$ |
| Ethanol | 2.3827 | 3.1771 | 198.24 | 2653.4 | 0.03238 | $[33]$ |


Table 2
Pure-component parameters of IL-ions used in this work. All ILs are parameterized as non-associating compounds.
| IL-ion | $m_{i}^{\text {seg }}$ | $\sigma_{i} / \AA$ | $u_{i} / k_{B} / \mathrm{K}$ | Ref. |
| :--- | :--- | :--- | :--- | :--- |
| $\left[\mathrm{C}_{2} \mathrm{mim}\right]^{+}$ | 1.4872 | 3.5926 | 206.4900 | $[31]$ |
| $\left[\mathrm{C}_{4} \mathrm{mim}\right]^{+}$ | 2.4805 | 3.6371 | 218.1441 | $[31]$ |
| $[\mathrm{TfO}]^{-}$ | 3.7432 | 3.8771 | 509.3113 | $[34]$ |


Table 3
ePC-SAFT pure-component parameters and binary interaction parameters $k_{i j}$ with water for inorganic ions used in this work for the first calculation of $K_{S P}$ from water.[15] Segment number is unity for all ions.
| Ion | $\sigma_{i} / \AA$ | $u_{i} / k_{B} / \mathrm{K}$ | $k_{\text {ion,Water }}$ |
| :--- | :--- | :--- | :--- |
| $\mathrm{Li}^{+}$ | 2.8449 | 360.00 | -0.25 |
| $\mathrm{Na}^{+}$ | 2.8232 | 230.00 | $-7.981 \bullet 10^{-4} \mathrm{~T} / \mathrm{K}+2.38$ |
| $\mathrm{K}^{+}$ | 3.3417 | 200.00 | $-4.012 \bullet 10^{-4} T / K+1.396$ |
| $\mathrm{F}^{-}$ | 1.7712 | 275.00 | 0 |
| $\mathrm{Cl}^{-}$ | 2.7560 | 170.00 | -0.25 |
| $\mathrm{Br}^{-}$ | 3.0707 | 190.00 | -0.25 |
| $\mathrm{I}^{-}$ | 3.6672 | 200.00 | -0.25 |


was linearized at logarithmic scale over the inversed temperature. This yields the gradient $\Delta h^{S L}$, i.e. the melting enthalpy (c.f. righthand side of Equation (27)). Applying the van 't Hoff equation, $K_{S P}$ was extrapolated to lower temperatures (e.g., 298 K). Then, the solubility of the metastable anhydrate at the given lower temperature was predicted with ePC-SAFT advanced. The application of the so-derived $K_{S P}$ can indeed be applied to predict the salt solubility also in other solvents given that crystal solvates do not form. Certainly, for lithium salts, the temperature range and the higher order of formed hydrates circumvents the application of the van 't Hoff equation. However, it is comfortable to apply the methodology to salts with anhydrates available at temperatures close to the temperature of interest.

## 3. Results and Discussion

### 3.1. Pure-component and binary interaction parameters

ePC-SAFT advanced demands pure-component parameters for all components. Table 1 and Table 2 list pure-component parameters for the organic solvents, water and the IL-ions, respectively. All non-electrolyte solvents were modeled as associating fluids with a 2 B associating scheme with parameters from literature. The ILs were characterized by an ion-specific approach as non-associating ions with a non-unity segment number. Note, that the temperature dependency of the segment diameter for inorganic ions is neglected as suggested in ePC-SAFT revised from 2014 [15], while it was explicitly accounted for the IL-ions as suggested in the work of Ji et al [31]. The inorganic ions forming the alkali-halide salts are presented in Table 3 together with the binary interaction parameters $k_{i j}$ with water from ref [15]. In the same parameter estimation, interaction parameters between the inorganic ions have been regressed (Table 4). It is important to state that all ion-related parameter estimation was performed by using data of aqueous salt solutions. The transferability to non-aqueous systems was pre-

![](https://cdn.mathpix.com/cropped/225b84f4-84d1-4768-92e8-fb2e12773952-06.jpg?height=374&width=1102&top_left_y=193&top_left_x=482)
Fig. 3. a) Degree of dissociation $\alpha$ over salt mole fraction for the salt LiCl in methanol (triangles) and ethanol (stars) at 298 K and 1 bar. Symbols are data calculated as defined in Equation (31) and lines are predicted degrees of dissociation with the Bjerrum treatment. b) Degree of dissociation $\alpha$ over IL mole fraction for the $\left[\mathrm{C}_{2}\right.$ mims][TfO] in water at 298 K , symbols are experimental data from literature.[42] The minimum of IL dissociation is visible at low IL concentration before approaching a limiting value of $\alpha=0.75$ for the pure IL. Line depicts predictive results using the Bjerrum treatment.

Table 4
ePC-SAFT binary interaction parameters $k_{i j}$ between halide anions and alkali cations used in this work from ref.[15]
| $k_{\text {cat,an }}$ | $\mathrm{Li}^{+}$ | $\mathrm{Na}^{+}$ | $\mathrm{K}^{+}$ | $\mathrm{Cs}^{+}$ |
| :--- | :--- | :--- | :--- | :--- |
| $\mathrm{F}^{-}$ | 0 | 0.665 | 1.000 | 1.000 |
| $\mathrm{Cl}^{-}$ | 0.669 | 0.317 | 0.064 | -0.417 |
| $\mathrm{Br}^{-}$ | 0.591 | 0.290 | -0.102 | -0.670 |
| $\mathrm{I}^{-}$ | 0.002 | 0.018 | -0.312 | -1.000 |


a) All salts were modeled with a similar dielectric constant that is a mean of available experimental data.[40]

Table 5
Dielectric constants for solvents and salts applied in this work. All salts are modeled with a similar dielectric constant that is a mean of available experimental data.
| Component | Dielectric constant $/ \mathrm{C} \cdot \mathrm{Vm}^{-1}$ | Ref. |
| :--- | :--- | :--- |
| Water | $-105.2 \ln T+677.480$ | [37] |
| Methanol | $-0.192 T+90.09$ | This work, [38] |
| Ethanol | $-0.146 T+68.47$ | This work, [39] |
| $\left[\mathrm{C}_{4}\right.$ mima $][\mathrm{Tf} \mathrm{O}]$ | 13.2 | [36] |
| Salts | 8 | a) |


sented in part one this research series by predicting MIACs in organic solvents [1].

The applied dielectric constants of the solvents and salts are summarized in Table 5. The concentration dependence of the dielectric constant is a vital part of the new methodology within ePC-SAFT. A linearized connection in salt mole fraction of the dielectric constants was found sufficient for the application to liquidliquid equilibria[27] and MIACs [1]. All inorganic salts were assumed to have a similar value of $\varepsilon_{r}=8$ that represents most of the salts in good estimation and is almost unaltered by temperature [35]. The dielectric constant for solvents was taken from literature and was assumed to depend linearly on temperature. For ILs, the dielectric constant was measured only at 298 K in the literature [36].

### 3.2. Ion paring and the Bjerrum treatment

Information of ion pairing is accessible through experimental conductivity measurements. For the dissociation of ions in water, a great amount of literature is available. For non-aqueous systems, conductivity data and thus dissociation degrees are less available. Under the assumptions that ions are independent of each other in their motion and that the degree of dissociation $\alpha$ becomes unity for an infinitely diluted solution, $\alpha$ may be calculated from the ratio of the molar conductivity $\Lambda^{\mathrm{c}}$ and the limiting molar conductiv-

$$
\begin{align*}
& \text { ity } \Lambda^{0} \text { (Equation (31)). } \\
& \alpha=\frac{\Lambda^{c}}{\Lambda^{0}} \tag{31}
\end{align*}
$$

The molar conductivity $\Lambda^{c}$ of an electrolyte relates the specific conductivity x to its molar concentration ( $\Lambda^{c}=x / c$ ). $\Lambda^{0}$ therefore refers to an infinitely diluted solution, in which the motion of the ion is counterbalanced solely by the interactions with the solvent molecules [41]. Both conductivity data are available in literature for salts in various solvents. The Bjerrum treatment for ion pairing is validated exemplarily for the dissociation of LiCl in ethanol and methanol, respectively, as depicted in Fig. 3a. For both systems, the dissociation is correctly reproduced within the Bjerrum treatment. More details on calculations are presented in the ESI and Table S1. Also given are data for dissociation of salts in other solvents and the respective results for the Bjerrum treatment. For common organic solvents considered in this work (e.g., alcohols), the electrostatic forces are considerably larger than classical physical short-range forces, leading to the monotonic decrease in dissociation for LiCl . The cut-off for association within the Bjerrum treatment is solely related to electrostatic potential not taking classical physical short-range interactions into account, only calculating the threshold for ion pairing from electrostatic potentials. This is a fair consideration as long as other forces are not undercutting the electrostatic potentials. This suggests that the Bjerrum treatment can hardly be adopted to systems containing complex mixtures of organic compounds, deep eutectic solvents or ILs, where classical physical short-range forces play major roles.

Indeed, the failure of the Bjerrum approach to ILs can be observed in Fig. 3b). ILs are in liquid state around room temperature, caused by strong hydrogen bonding and dispersive energies. Such effects were not included into the Bjerrum approach. Bjerrum treatment predicts a characteristic decrease of the IL dissociation, while experimental data exhibit a minimum at very low concentration followed by approaching a value of almost 0.8 . Due to the assumption of an ideal-mixture behavior of the ion pair, the Bjerrum treatment predicts a monotone decrease of the number density of the free ions and an increase of the number of ion pairs with increasing electrolyte concentration [43,44]. There is still debate on the origin of this conductance minimum. Previous interpretations suggested the existence of an equilibrium between free ions and ion pairs to form charged clusters [44]. However, further theories developed from the Bjerrum treatment (or employing other expressions for the dissociation degree)[45] could predict the appearance of this conductance minimum without involving triplets or higher clusters [24,46]. Those theories mainly focused on accounting for the high dipole moment of ion pairs, which at high ion and dipole concentration causes non-negligible dipole-ion as well as dipole-dipole interactions,[46] and changes the relative dielectric constant of the system [47]. Two theories worth mention-

![](https://cdn.mathpix.com/cropped/225b84f4-84d1-4768-92e8-fb2e12773952-07.jpg?height=391&width=1100&top_left_y=193&top_left_x=484)
Fig. 4. Influence of salt concentration on the activity coefficient of methanol for the salts left: LiCl and right: LiBr at 298 K and 1 bar. Symbols are experimental data, lines are predictions with ePC-SAFT revised from $2014{ }^{15}$ (black line), ePC-SAFT advanced with the concentration dependent dielectric constant and the altered Born contribution, and ePC-SAFT advanced additionally including the Bjerrum treatment (orange line).

ing are the Fisher-Levin and the Weiss-Schröer theory [48]. To sum up, the Bjerrum treatment was not considered to model properties in ILs e.g., solubility of alkali-halides in ILs.

### 3.3. Solvent activity coefficients in presence of salts

Further, the Bjerrum treatment was validated against activity coefficients of the organic solvent with rising concentration of the salt. In part one of ePC-SAFT advanced[1] the MIACs of salts were quantitatively predicted. The secret behind was use of the altered Born contribution that contains a salt-concentration dependent dielectric constant same as the Debye-Hückel theory. Here, the ePCSAFT advanced approach with and without addition of the Bjerrum treatment into is compared to ePC-SAFT revised from 2014 [15]. Fig. 4 exemplifies the results for the three ePC-SAFT models for the activity coefficients of methanol upon addition of left: LiCl or right: LiBr. ePC-SAFT revised from $2014{ }^{15}$ yields a large overestimation of the activity coefficient. This behavior is explained by the missing contribution of the concentration-dependent dielectric constant included in the original Born term. Consequently, ePC-SAFT advanced with an altered Born contribution gives a more realistic expression of the concentration dependence of the activity coefficient. Including the Bjerrum treatment into ePC-SAFT advanced, the activity coefficient is met quantitatively over the whole concentration range. This is attributed to the reduced concentration of free ions that contribute to the Debye-Hückel energy. For that matter, the Bjerrum treatment is a valid extension to the model from a physical standpoint. Similar results with ePC-SAFT including the Bjerrum treatment outperforming the original ePC-SAFT were found for LiCl in ethanol (c.f. Figure S1 in the ESI).

## 3.4. $\mathbf{K}_{\mathbf{S P}}$ as access to the solubility of metastable anhydrates and salt solubility in non-aqueous systems

Accessing the solubility of the metastable anhydrate is important to predict salt solubility in non-aqueous solvents. $\mathrm{K}_{\mathrm{SP}}$ that is extrapolated to account for the solubility of the metastable anhydrate can therefore also be applied to other solvents. Exemplarily shown here is the prediction of the solubility of metastable anhydrates for $\mathrm{NaBr} . \mathrm{NaBr}$ forms a stable dihydrate below 323 K . Solubility data for both, the dihydrate and the anhydrate are available in literature [49]. The solubility product $K_{S P}$ was predicted with ePC-SAFT advanced for the anhydrate above 323 K and then linearized at logarithmic scale over the inverse temperature following the van 't Hoff approach. The respective plot for $K_{S P}$ of the hydrate (orange symbols) and the anhydrate (green symbols) is presented in Fig. 5. The solid line depicts the values for $K_{S P}$ calculated with ePC-SAFT advanced, which was obtained by using the experimental solubility data and the ePC-SAFT predicted gMIACs at the respective solubility concentration and temperature.

![](https://cdn.mathpix.com/cropped/225b84f4-84d1-4768-92e8-fb2e12773952-07.jpg?height=343&width=512&top_left_y=758&top_left_x=1244)
Fig. 5. van 't Hoff plot for the logarithmic solubility product $K_{S P}$ over the inversed temperature of the salt NaBr in water. The green symbols are $K_{S P}$-values calculated with ePC-SAFT for the anhydrate, orange symbols represent the $K_{S P}$-values for the dihydrate. The solid green line depicts the $K_{S P}$-values and the dashed green line the extrapolated values using the van 't Hoff equation.

![](https://cdn.mathpix.com/cropped/225b84f4-84d1-4768-92e8-fb2e12773952-07.jpg?height=365&width=507&top_left_y=1344&top_left_x=1249)
Fig. 6. T,x-diagram of the solubility of NaBr in water. Orange symbols are experimental solubilities of the dihydrate; green symbols are solubility data for the anhydrate.[49] The solid green line depicts the solubility of the anhydrate and the dashed green line depicts the hypothetical solubility of the metastable anhydrate predicted with ePC-SAFT using the extrapolated $K_{S P}$-values from the van 't Hoff plot.

After extrapolating the linearized $K_{S P}$ to lower temperatures (dashed green line), the solubility of the metastable salt is predicted by using $K_{S P}$ as input value into Equation (28). Fig. 6 gives the experimental solubility of the NaBr dihydrate at temperatures below 323 K and of the anhydrate above this temperature. ePCSAFT predictions of the anhydrate solubility (dashed green line) show a much higher solubility of the anhydrate compared to the dihydrate below 323 K . The dihydrate therefore represents the stabilized and precipitating solute consequently bringing a smaller solubility. Additionally, a more distinct temperature dependence is found for the solubility of the dihydrate.

The methodology and the derived $K_{S P}$-values are not only applicable to the solubility of the metastable salts in water. In fact, this procedure allows to predict the salt solubility in other solvents that also precipitate the non-solvated salts at lower temperatures. The approach is not validated for lithium salts caused by the very high temperature at which the anhydrous salt takes the stable form.

Table 6
Overview of the experimental data of solubilities of alkalihalide salts in water.
| Salt | Hydrate* | T range / K | NP | Reference |
| :--- | :--- | :--- | :--- | :--- |
| LiCl | yes | 273.15-373.15 | 33 | [50,51] |
| LiBr | yes | 273.15-373.15 | 12 | [52] |
| NaF | no | 273.15-373.15 | 12 | [52] |
| NaCl | no | 273.15-373.15 | 20 | [49,51] |
| NaBr | treated | 273.15-373.15 | 22 | [49,52] |
| NaI | treated | 273.15-373.15 | 12 | [52] |
| KF | no | 273.15-353.15 | 6 | [51] |
| KCl | no | 273.15-393.15 | 30 | [50,51] |
| KBr | no | 273.15-373.15 | 20 | [50,52,53] |
| KI | no | 273.15-373.15 | 12 | [52] |


Table 7
Overview of the experimental data of solubilities of alkali-halide salts in the investigated organic solvents.
| Salt | Organic solvent | $T$-range / K | NP | Ref. |
| :--- | :--- | :--- | :--- | :--- |
| LiCl | Methanol | 293.15-333.15 | 9 | [54] |
|  | Ethanol | 293.15-333.15 | 14 | [54,55] |
|  | Methanol | 293.15-338.15 | 18 | [54] |
|  | Ethanol | 293.15-333.15 | 14 | [55,56] |
| NaF | Methanol | 293.15-328.15 | 5 | [57] |
|  | Ethanol | 293.15-328.15 | 5 | [57] |
| NaCl | Methanol | 293.15-333.15 | 14 | [49,54] |
|  | Ethanol | 298.15-348.15 | 3 | [49] |
| NaBr | Methanol | 298.15-333.15 | 5 | [49] |
|  | Ethanol | 298.15-348.15 | 3 | [49] |
| NaI | Methanol | 283.15-333.15 | 9 | [58] |
|  | Ethanol | 283.15-373.15 | 5 | [59] |
| KF | Methanol | 293.15-328.15 | 5 | [59] |
|  | Ethanol | 293.15-328.15 | 5 | [59] |
| KCl | Methanol | 293.15-333.15 | 19 | [49,54,55] |
|  | Ethanol | 293.15-333.15 | 12 | [49,54] |
| KBr | Methanol | 298.15-333.15 | 5 | [53] |
|  | Ethanol | 293.15-353.15 | 3 | [53] |
| KI | Ethanol | 278.25-342.75 | 18 | [60,61] |


### 3.5. The influence of Bjerrum treatment on modeled solubility of alkali-halides salts in organic solvents

The approach suggested in Fig. 6 is henceforward applied to predict salt solubility in organic solvents including the Bjerrum treatment to account for ion pairing. Table 6 presents available data for salt solubility in water, which is used as input data to determine $K_{S P}$. The formation of hydrates is commented, and for lithium-based salts modeling is impossible due to the formation of hydrates over a broad temperature range. NaBr and NaI were treated with the methodology explained in Fig. 6. Table 7 gives an overview about solubility data for alkali-halide salts in the organic solvents methanol and ethanol. Experimental data are available in a larger variety of salts, organic solvents, and temperature range.

* Denotes if a hydrates is formed at any given temperature and if it was treated with the methodology for metastable anhydrates.

In general, the predictions of salt solubility with ePC-SAFT advanced including the altered Born contribution and $\varepsilon(\bar{x})$ are in the order of magnitude of the experimental data. The solubility increases with larger inorganic ions forming the salts are predicted qualitatively correct (c.f. Fig. 7). To the best of our knowledge, this work shows for the first time that the prediction of the solubility of salts in pure-organic solvents is achievable with an EoS. It should be noted that ePC-SAFT revised (2014) predicted highly overestimated solubility, and that even $k_{i j}$ parameter regressions were not possible. The reason behind is a missing contribution to the Gibbs energy of solvation that is required to decrease the activity coefficient of the solvent with increased ion concentration (c.f. Fig. 4).

![](https://cdn.mathpix.com/cropped/225b84f4-84d1-4768-92e8-fb2e12773952-08.jpg?height=387&width=593&top_left_y=193&top_left_x=1206)
Fig. 7. Predicted solubility of alkali-halide salts in the alcohols methanol and ethanol at 298 K and 1 bar. Gray bars are experimental data, green bars are predictions with ePC-SAFT advanced including the altered Born contribution and $\varepsilon(\overline{\mathrm{x}})$, orange bars are predictions additionally including the Bjerrum treatment.

Table 8
ePC-SAFT binary interaction parameters $k_{i j}$ between ions and organic solvents used to correlate solubility of alkalihalides in organic solvent.
|  | $\mathrm{Na}^{+}$ | $\mathrm{K}^{+}$ | $\mathrm{Cl}^{-}$ | $\mathrm{Br}^{-}$ | $\mathrm{I}^{-}$ |
| :--- | :--- | :--- | :--- | :--- | :--- |
| Methanol | -0.31 | 0.47 | -0.21 | -0.42 | - |
| Ethanol | 0.42 | -0.2 | -0.15 | -0.35 | -0.38 |


Figure 7 shows the solubility predictions at 298 K for monovalent salts and the respective alcohols methanol and ethanol. It can be seen that ePC-SAFT advanced cannot predict the extremely high solubility of NaBr . The solubility is underestimated for both approaches (with and without the Bjerrum treatment) without reaching a satisfying agreement with the experimental data. Still, the prediction results including Bjerrum are slightly better for the solubility at 298 K . In contrast to the sodium-based salts, the solubility of the potassium-based salts is predicted in better agreement. The Bjerrum treatment effectively reduces the number of free ions that are accounted for in the Debye-Hückel theory and therefore reducing Coulombic forces, intensifying favored interactions within the solvent and increasing the solubility. Generally, this assumes that the Bjerrum treatment will only give better agreement if the prediction with ePC-SAFT advanced is lower than the experimental data. However, the influence of the Bjerrum treatment and thus ion pairing on the solubility prediction results is not very pronounced for the investigated alcohol solvents. Other salt-solvent mixtures might benefit largely from the application of ion pairs. This was already implied by results for COSMO-RS [26].

In order to quantitatively model salt solubility in solvents interaction parameters are required. Table 8 lists binary interaction parameters $k_{i j}$ between ion and alcohol. No temperature dependence was assumed, and correlations were investigated at 298 K where the pure-component parameters for the ions are valid. Correlations were carried out for the model including the Bjerrum treatment and are further studied with the infinite dilution property Gibbs energy of transfer. Results for the solubility correlations are presented in Fig. 8a. Fig. 8b shows the respective impact on transfer properties for the applied ions.

Varying the dispersion energy between alcohol and ions by applying these $k_{i j}$ allowed quantitative agreement with the experimental data. The very large solubility of NaBr in methanol (and ethanol) was successfully modeled within this approach. Applying the $k_{i j}$-values also to Gibbs energy of transfer, the value for potassium in methanol represents the experimental data within the value of the assumed standard deviation from literature. The same was found for chloride and with some limitations also for sodium. Overall, the correlation to salt solubility also gives a better representation of $\Delta^{T r} G_{i}$, meaning that the correlation is of sound

![](https://cdn.mathpix.com/cropped/225b84f4-84d1-4768-92e8-fb2e12773952-09.jpg?height=398&width=1268&top_left_y=191&top_left_x=400)
Fig. 8. a) Solubility of alkali-halide salts in the alcohols methanol and ethanol at 298 K and 1 bar with ePC-SAFT advanced including the Bjerrum treatment using binary interaction parameters $k_{i j}$ from Table 8. b) Gibbs energy of transfer for the investigated ions of the alkali-halide salts. Gray: experimental data; orange: predictions using ePC-SAFT advanced including the Bjerrum treatment, red: correlations with $\mathrm{k}_{\mathrm{ij}}$ listed in Table 8.

Table 9
Overview of the experimental data of solubilities of alkali-halide salts in the investigated IL [ $\mathrm{C}_{4}$ mim][TfO][62] in the temperature range from 298 K to 378 K .
| Salt | Hydrate formation in water | NP |
| :--- | :--- | :--- |
| LiF | Yes | 5 |
| LiCl | Yes | 5 |
| LiBr | Yes | 5 |
| LiI | Yes | 5 |
| NaF | No | 5 |
| NaCl | No | 5 |
| NaBr | Treated | 5 |
| NaI | Treated | 5 |
| KF | No | 5 |
| KCl | No | 5 |
| KBr | No | 5 |
| KI | No | 5 |


physical foundation and not only a way to represent the experimental data.

### 3.6. Solubility of alkali-halides salts in ionic liquids

Experimental solubility data of alkali-halide salts in ILs is rather scarce compared to organic solvents. A holistic investigation was only available for systems with the IL [ $\mathrm{C}_{4}$ mim][TfO] [62]. The overview on experimental data used in this work is given in Table 9 and details on further work is given in Figure S2 in the ESI. A larger range of experimental data would benefit a more generalized investigation and would allow the influence of IL-anion or IL-cation to be accounted for. The formation of hydrates of higher order and anhydrates available only at very elevated temperatures rule out lithium-based salts for the prediction of salt solubility in the IL. Still, the salts are given in the overview for clarity. Predictions are performed applying pure-component and binary interaction parameters in Table 1 to Table 4.

A summary on predicted salt solubility in the ILs is given in Fig. 9 for 298 K . The results show a qualitative agreement with the experimental data for all salts investigated in terms of magnitude of concentration at 298 K using ePC-SAFT advanced without the Bjerrum treatment. Still, the results lag in quantitative agreement for iodide-based salts. Interestingly, compared to the experimental data in alcohols (Fig. 7), the unexpectedly high solubility is not found for NaBr . It is presumed that the larger radii of the iodideanion counteract the effective Coulomb forces that are predicted by ePC-SAFT advanced. Note, that ePC-SAFT revised did not allow obtaining results in the vicinity of the experimental data (not shown).

Certainly, based on the already good results for the predicted salt solubility at 298 K , correlations are possible to quantitatively meet the experimental data. Exemplarily, the solubility of sodium-

![](https://cdn.mathpix.com/cropped/225b84f4-84d1-4768-92e8-fb2e12773952-09.jpg?height=404&width=593&top_left_y=764&top_left_x=1206)
Fig. 9. Predicted solubility of alkali-halide salts in the investigated ILs [ $\mathrm{C}_{4}$ mim][TfO] at 298 K and 1 bar. Gray bars are experimental data,[62] green bars are predictions using ePC-SAFT advanced without the Bjerrum treatment.

![](https://cdn.mathpix.com/cropped/225b84f4-84d1-4768-92e8-fb2e12773952-09.jpg?height=386&width=593&top_left_y=1318&top_left_x=1206)
Fig. 10. Solubility of sodium-halide salts in the investigated IL [ $\mathrm{C}_{4}$ mim][TfO]] at 298 K and 1 bar using the binary interaction parameter $k_{\mathrm{Na}^{+},[\mathrm{Tf}]^{-}}=0.17$. Gray: experimental data,[62] red: correlations with ePC-SAFT advanced without the Bjerrum treatment.

based salts is correlated with a single binary interaction parameter between the sodium-cation and the IL-anion with a value of $k_{i j}=$ 0.17 The results are depicted in Fig. 10 giving satisfying agreement for all sodium-based salts.

## 4. Conclusion

In this work, ePC-SAFT advanced is combined with the Bjerrum treatment to predict and model the solubility of alkali-halide salts in pure-organic and pure-ionic systems. The Bjerrum treatment based on mass balance was introduced to account for ion pairing caused by electrostatic potentials. The predicted dissociation of salts in organic solvents was in quantitative agreement with experimental data. The Bjerrum treatment is suitable for organic solvents but fails to the admittedly unexpected dissociation behavior of ILs due to lack in physical short-range contributions to the model. For further validation, activity coefficients of solvents with increasing salt concentration were investigated and met with quantitative agreement. Inclusion of the Bjerrum treatment into
the Debye Hückel term did not significantly modify the prediction results.

A methodology to access salts forming hydrates is exemplarily shown utilizing the solubility product $K_{S P}$, giving the concentration of the hypothetical metastable anhydrate at lower temperatures. The derived $K_{S P}$-values are thereafter used as input for the prediction of salt solubility in alcohol solvents. After validation, ePC-SAFT advanced was tested against experimental solubility of monovalent alkali-halide salts in the organic solvents methanol and ethanol and in the IL [ $\mathrm{C}_{4}$ mim][TfO] again utilizing the solubility product $K_{S P}$. Again, inclusion of the Bjerrum treatment into the Debye Hückel term did not significantly modify the prediction results.

In sum, ePC-SAFT advanced allowed predicting reasonably salt solubility in ILs and organic solvents. The prediction of salt solubility in both solvent classes is satisfying given the complexity of the investigated systems, exhibiting ion pairing, ion solvation and the application of pure-component parameters estimated from data of aqueous electrolyte solutions. Correlations using binary interaction parameters were also successful to further improve the model accuracy. As a final statement, solubility results were only possible with ePC-SAFT advanced that includes the altered Born contribution and $\varepsilon(\bar{x})$, while ePC-SAFT revised (2014) without these failed.

## Acknowledgement

This work has been supported by the German Science Foundation (DFG) within the priority program SPP 1708 "Material Synthesis Near Room Temperature" (grant HE 7165/7-1). The authors would like to thank their student Eric Rosthoff for his much-valued work on literature and association of ions.

## Declaration of Competing Interest

The authors declare no conflict of interest.

## Supplementary materials

Supplementary material associated with this article can be found, in the online version, at doi:10.1016/j.fluid.2021.112989.

## References

[1] M. Bülow, M. Ascani, C. Held, ePC-SAFT advanced - Part I: Physical meaning of including a concentration-dependent dielectric constant in the born term and in the Debye-Hückel theory, Fluid Phase Equilibria 535 (2021) 112967.
[2] R. Black, B. Adams, L.F. Nazar, Non-Aqueous and Hybrid Li-O2 Batteries, Adv. Energy Mater. 2 (2012) 801-815.
[3] A. Ponrouch, D. Monti, A. Boschin, B. Steen, P. Johansson, M.R. Palacín, Non-aqueous electrolytes for sodium-ion batteries, J. Mater. Chem. A 3 (2015) 22-42.
[4] I.Y. Shilov, A.K. Lyashchenko, The Role of Concentration Dependent Static Permittivity of Electrolyte Solutions in the Debye-Hückel Theory, The journal of physical chemistry. B 119 (2015) 10087-10095.
[5] J. Vincze, M. Valiskó, D. Boda, The nonmonotonic concentration dependence of the mean activity coefficient of electrolytes is a result of a balance between solvation and ion-ion correlations, J. Chem. Phys. 133 (2010) 154507.
[6] Q. Lei, B. Peng, L. Sun, J. Luo, Y. Chen, G.M. Kontogeorgis, X. Liang, Predicting activity coefficients with the Debye-Hückel theory using concentration dependent static permittivity, AIChE J (2020) 66.
[7] J. Wu, J.M. Prausnitz, Phase Equilibria for Systems Containing Hydrocarbons, Water, and Salt: An Extended Peng-Robinson Equation of State, Ind. Eng. Chem. Res. 37 (1998) 1634-1643.
[8] J.A. Myers, S.I. Sandler, R.H. Wood, An Equation of State for Electrolyte Solutions Covering Wide Ranges of Temperature, Pressure, and Composition, Ind. Eng. Chem. Res. 41 (2002) 3282-3297.
[9] Y. Lin, K. Thomsen, J.-C.de Hemptinne, Multicomponent equations of state for electrolytes, AIChE J 53 (2007) 989-1005.
[10] R. Inchekel, J.-C.de Hemptinne, W. Fürst, The simultaneous representation of dielectric constant, volume and activity coefficients using an electrolyte equation of state, Fluid Phase Equilibria 271 (2008) 19-27.
[11] B. Maribo-Mogensen, K. Thomsen, G.M. Kontogeorgis, An electrolyte CPA equation of state for mixed solvent electrolytes, AIChE J 61 (2015) 2933-2950.
[12] S. Ahmed, N. Ferrando, J.-C. Hemptinne, J.-P. de; Simonin, O. Bernard, O. Baudouin, Modeling of mixed-solvent electrolyte systems, Fluid Phase Equilibria 459 (2018) 138-157.
[13] J. Rozmus, J.-C.de Hemptinne, A. Galindo, S. Dufal, P. Mougin, Modeling of Strong Electrolytes with ePPC-SAFT up to High Temperatures, Ind. Eng. Chem. Res. 52 (2013) 9979-9994.
[14] J.M.A. Schreckenberg, S. Dufal, A.J. Haslam, C.S. Adjiman, G. Jackson, A Galindo, Modelling of the thermodynamic and solvation properties of electrolyte solutions with the statistical associating fluid theory for potentials of variable range, Molecular Physics 112 (2014) 2339-2364.
[15] C. Held, T. Reschke, S. Mohammad, A. Luza, G. Sadowski, ePC-SAFT revised, Chem. Eng. Res. Des. 92 (2014) 2884-2897.
[16] R. Heyrovska, Degrees of Dissociation and Hydration Numbers of Alkali Halides in Aqueous Solutions at $25^{\circ} \mathrm{C}$ (Some up to Saturation), Croatica Chemica Acta (1997) 39-54.
[17] C. Held, L.F. Cameretti, G. Sadowski, Modeling aqueous electrolyte solutions, Fluid Phase Equilibria 270 (2008) 87-96.
[18] T. Reschke, C. Brandenbusch, G. Sadowski, Modeling aqueous two-phase systems: I. Polyethylene glycol and inorganic salts as ATPS former, Fluid Phase Equilibria 368 (2014) 91-103.
[19] C. Held, T. Reschke, R. Müller, W. Kunz, G. Sadowski, Measuring and modeling aqueous electrolyte/amino-acid solutions with ePC-SAFT, The Journal of Chemical Thermodynamics 68 (2014) 1-12.
[20][ A. Wangler, G. Sieder, T. Ingram, M. Heilig, C. Held, Prediction of CO2 and H2S solubility and enthalpy of absorption in reacting N-methyldiethanolamine /water systems with ePC-SAFT, Fluid Phase Equilibria 461 (2018) 15-27.
[21] D. Pabsch, C. Held, G. Sadowski, Modeling the CO 2 Solubility in Aqueous Electrolyte Solutions Using ePC-SAFT, J. Chem. Eng. Data 65 (2020) 5768-5777.
[22] N Bjerrum, Untersuchungen über Ionenassoziation. I. Der Einfluss der Ionenassoziation auf die Aktivität der Ionen bei Mittleren Assoziationsgraden, von Niels Bjerrum, B. Lunos, 1926.
[23] Y. Marcus, G. Hefter, Ion pairing, Chemical reviews 106 (2006) 4585-4621.
[24] Y. Jiang, H. Nadolny, S. Käshammer, S. Weibels, W. Schröer, H. Weingärtner, The ion speciation of ionic liquids in molecular solvents of low and medium polarity, Faraday discussions 154 (2012) 391-407 discussion 439-64, 465-71.
[25] W. Schröer, On the chemical and the physical approaches to ion association, Journal of Molecular Liquids 164 (2011) 3-10.
[26] S. Müller, A. González de Castilla, C. Taeschler, A. Klein, I. Smirnova, Calculation of thermodynamic equilibria with the predictive electrolyte model COS-MO-RS-ES: Improvements for low permittivity systems, Fluid Phase Equilibria 506 (2020) 112368.
[27] M. Bülow, X. Ji, C. Held, Incorporating a concentration-dependent dielectric constant into ePC-SAFT. An application to binary mixtures containing ionic liquids, Fluid Phase Equilib 492 (2019) 26-33.
[28] H. Veith, C. Luebbert, G. Sadowski, Correctly Measuring and Predicting Solubilities of Solvates, Hydrates, and Polymorphs, Crystal Growth \& Design 20 (2) (2020) 723-735 Crystal Growth \& Design20, 723-735.
[29] K.B. Lipkowitz, R. Larter, T.R. Cundari, Reviews in Computational Chemistry, John Wiley \& Sons, Inc, Hoboken, NJ, USA, 2003.
[30] Bülow, M.; Danzer, A.; Held, C. Liquid-Liquid Equilibria of Binary and Ternary Systems Containing Ionic Liquids. In Encyclopedia of Ionic Liquids; Zhang, S., Ed.; Springer Singapore: Singapore, 2020; 1-7.
[31] X. Ji, C. Held, G. Sadowski, Modeling imidazolium-based ionic liquids with ePC-SAFT, Fluid Phase Equilib 335 (2012) 64-73.
[32] D. Fuchs, J. Fischer, F. Tumakaka, G. Sadowski, Solubility of Amino Acids: Influence of the pH value and the Addition of Alcoholic Cosolvents on Aqueous Solubility, Ind. Eng. Chem. Res. 45 (2006) 6578-6584.
[33] J. Gross, G. Sadowski, Application of the Perturbed-Chain SAFT Equation of State to Associating Systems, Ind. Eng. Chem. Res. 41 (2002) 5510-5515.
[34] X. Ji, C. Held, Modeling the density of ionic liquids with ePC-SAFT, Fluid Phase Equilib 410 (2016) 9-22.
[35] R.P. Lowndes, D.H. Martin, Dielectric constants of ionic crystals and their variations with temperature and pressure, Proc. R. Soc. Lond. A 316 (1970) 351-375.
[36] H. Weingärtner, The Static Dielectric Constant of Ionic Liquids, Zeitschrift für Physikalische Chemie 220 (2006) 1395-1405.
[37] W.B. Floriano, M.A.C. Nascimento, Dielectric constant and density of water as a function of pressure at constant temperature, Braz. J. Phys. 34 (2004) 38-41.
[38] M.T. Khimenko, V.V. Aleksandrov, N.N. Gritsenko, Polarizability and radii of molecules of some pure liquids, Zh.Fiz.Khim (1973) 2914-2915.
[39] R.M. Shirke, A. Chaudhari, N.M. More, P.B. Patil, Temperature dependent dielectric relaxation study of ethyl acetate - Alcohol mixtures using time domain technique, Journal of Molecular Liquids 94 (2001) 27-36.
[40] C. Andeen, J. Fontanella, D. Schuele, Low-Frequency Dielectric Constant of LiF, $\mathrm{NaF}, \mathrm{NaCl}, \mathrm{NaBr}, \mathrm{KCl}$, and KBr by the Method of Substitution, Physical Review B (1970) 5068-5073.
[41] R.A. Robinson, R.H. Stokes, Electrolyte Solutions: The Measurement and Interpretation of Conductance, Chemical Potential and Diffusion in Solutions of Simple Electrolytes, 2nd ed., Butterworth \& Co: London, 1959.
[42] O. Nordness, L.D. Simoni, M.A. Stadtherr, J.F. Brennecke, Characterization of Aqueous 1-Ethyl-3-Methylimidazolium Ionic Liquids for Calculation of Ion Dissociation, The journal of physical chemistry. B 123 (2019) 1348-1358.
[43] Y. Levin, M.E. Fisher, Criticality in the hard-sphere ionic fluid, Physica A: Statistical Mechanics and its Applications 225 (1996) 164-220.
[44] R.M. Fuoss, C.A. Kraus, Properties of Electrolytic Solutions. IV. The Conductance Minimum and the Formation of Triple Ions Due to the Action of Coulomb Forces 1, J. Am. Chem. Soc. 55 (1933) 2387-2399.
[45] W. Ebeling, Zur Theorie der Bjerrumschen Ionenassoziation in Elektrolyten, Zeitschrift für Physikalische Chemie (1968) 2380.
[46] H. Weingärtner, V.C. Weiss, W. Schröer, Ion association and electrical conductance minimum in Debye-Hückel-based theories of the hard sphere ionic fluid, J. Chem. Phys. 113 (2000) 762-770.
[47] H. Weingärtner, H.G. Nadolny, S. Käshammer, Dielectric Properties of an Electrolyte Solution at Low Reduced Temperature, The journal of physical chemistry. B 103 (1999) 4738-4743.
[48] V.C. Weiss, W. Schröer, Macroscopic theory for equilibrium properties of ion-ic-dipolar mixtures and application to an ionic model fluid, J. Chem. Phys. 108 (1998) 7747-7757.
[49] S.P. Pinho, E.A. Macedo, Solubility of $\mathrm{NaCl}, \mathrm{NaBr}$, and KCl in Water, Methanol, Ethanol, and Their Mixed Solvents, J. Chem. Eng. Data 50 (2005) 29-32.
[50] J.W. Lorimer, M. Salomon, C.L. Young, Eds. Solubility Data Series: Volume 47: Alkali Metal and Ammonium Chlorides in Water and Heavy Water (Binary Systems), Pergamon Press, 1991.
[51] B. Scrosati, C.A. Vincent, Eds. Alkali metal, alkaline-earth metal and ammonium halides, amide solvents; Solubility data series 11, Pergamon Press, Oxford, 1980.
[52] O. Söhnel, P. Novotny, Densities of Aqueous Solutions of Inorganic Substances, Physical Science Data (2021).
[53] S.P. Pinho, E.A. Macedo, Experimental measurement and modelling of KBr solubility in water, methanol, ethanol, and its binary mixed solvents at different temperatures, The Journal of Chemical Thermodynamics 34 (2002) 337-360.
[54] M. Li, D. Constantinescu, L. Wang, A. Mohs, J. Gmehling, Solubilities of NaCl, $\mathrm{KCl}, \mathrm{LiCl}$, and LiBr in Methanol, Ethanol, Acetone, and Mixed Solvents and Correlation Using the LIQUAC Model, Ind. Eng. Chem. Res. 49 (2010) 4981-4988.
[55] N. Xin, Y. Sun, M. He, C.J. Radke, J.M. Prausnitz, Solubilities of six lithium salts in five non-aqueous solvents and in a few of their binary mixtures, Fluid Phase Equilibria 461 (2018) 1-7.
[56] M.-Y. Li, L.-S. Wang, K.-P. Wang, B. Jiang, J. Gmehling, Experimental measurement and modeling of solubility of LiBr and LiNO 3 in methanol, ethanol, 1-propanol, 2-propanol and 1-butanol, Fluid Phase Equilibria 307 (2011) 104-109.
[57] F.G. Germuth, The solubilities of alkali bromides and fluorides in anhydrous methanol, ethanol, and butanol, Journal of the Franklin Institute 212 (1931) 343-349.
[58] E. Lloyd, C.B. Brown, D.G.R. Bonnell, W.J Jones, XC.-Equilibrium between alcohols and salts. Part II, J. Chem. Soc. 0 (1928) 658-666.
[59] G.J. Janz, R.P.T. Tomkins, J. Ambrose, Nonaqueous electrolytes handbook; Nonaqueous electrolytes handbook v. 1, Academic Press, New York, 1972.
[60] R.R. Pawar, S.B. Nahire, M. Hasan, Solubility and Density of Potassium Iodide in Binary Ethanol-Water Solvent Mixture at (298.15, 303.15, 308.15, and 313.15) K, J. Chem. Eng. Data 54 (2009) 1935-1937.
[61] B. Long, D. Zhao, W. Liu, Thermodynamics Studies on the Solubility of Inorganic Salt in Organic Solvents: Application to KI in Organic Solvents and Wa-ter-Ethanol Mixtures, Ind. Eng. Chem. Res. 51 (2012) 9456-9467.
[62] O. Kuzmina, E. Bordes, J. Schmauck, P.A. Hunt, J.P. Hallett, T. Welton, Solubility of alkali metal halides in the ionic liquid C4C1imOTf, Physical chemistry chemical physics : PCCP 18 (2016) 16161-16168.


[^0]:    * Corresponding author.

    E-mail address: christoph.held@tu-dortmund.de (C. Held).

