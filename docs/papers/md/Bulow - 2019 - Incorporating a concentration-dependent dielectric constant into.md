# Incorporating a concentration-dependent dielectric constant into ePC-SAFT. An application to binary mixtures containing ionic liquids 

Mark Bülow ${ }^{\mathrm{a}}$, Xiaoyan Ji ${ }^{\mathrm{b}}$, Christoph Held ${ }^{\mathrm{a}, \text { * }}$<br>${ }^{\mathrm{a}}$ Laboratory of Thermodynamics, Department of Biochemical and Chemical Engineering, TU Dortmund, Emil-Figge Str. 70, 44277, Dortmund, Germany<br>${ }^{\mathrm{b}}$ Division of Energy Science/Energy Engineering, Lulea University of Technology, 97187, Lulea, Sweden

## ARTICLE INFO

## Article history:

Received 2 February 2019
Received in revised form
7 March 2019
Accepted 12 March 2019
Available online 15 March 2019

## Keywords:

Electrolytes
Liquid-liquid equilibria
Prediction
Thermodynamics


#### Abstract

Primitive thermodynamic models for electrolyte solutions require the dielectric constant $\varepsilon$. This property strongly depends on the concentration of the electrolytes in the mixture. Neglecting this dependency might be reasonable for modeling solutions at low electrolyte concentrations. However, in solutions containing ionic liquids (ILs) and especially for the calculation of liquid-liquid equilibria (LLE) of systems with ILs, liquid phases often contain high IL concentrations. At such conditions, neglecting the influence of concentration on $\varepsilon$ is an oversimplification. In this work, an approach to account for the concentrationdependent dielectric constant within the Debye-Hückel theory was implemented into electrolyte Perturbed-Chain Statistical Associating Fluid Theory (original ePC-SAFT). This new approach was then applied to model LLE of binary mixtures containing water and commonly used hydrophobic ILs. These common ILs are comprised of the IL-cations $\left[\mathrm{C}_{\mathrm{n}} \mathrm{mim}\right]^{+},\left[\mathrm{C}_{\mathrm{n}} \mathrm{py}\right]^{+},\left[\mathrm{C}_{\mathrm{n}} \mathrm{mpy}\right]^{+},\left[\mathrm{C}_{\mathrm{n}} \mathrm{mpyr}\right]^{+},\left[\mathrm{C}_{4} \mathrm{~m}_{4} \mathrm{py}\right]^{+}$and the IL-anions $\left[\mathrm{BF}_{4}\right]^{-},\left[\mathrm{NTf}_{2}\right]^{-},\left[\mathrm{PF}_{6}\right]^{-},[\mathrm{TFO}]^{-}$. The LLE of binary mixtures water + IL were modeled at ambient pressure and different temperatures with the new ePC-SAFT and with the original ePC-SAFT [Ji et al. DOI: 10.1016/j.fluid.2012.05.029] without the concentration-dependent $\varepsilon$. Overall, the new approach within ePC-SAFT shows superior modeling as well as correlation capability compared to original ePC-SAFT, which was concluded by comparing both models with LLE data from literature.


© 2019 Elsevier B.V. All rights reserved.

## 1. Introduction

In the past two decades, ionic liquids (ILs) have raised a great interest in science, and thus ILs have been intensively investigated [1]. The responsible property is their very small vapor pressure. They additionally usually exhibit low thermal degradability as well as low melting points usually below 373 K . This makes them favorable solvents for future process intensification [2]. There is an almost unlimited number of ILs as they are consisting of highly alterable IL-ions and are thus seen as "designer" solvents for different applications in the chemical industry (organic, inorganic, polymeric, catalytic) for which the conversion or selectivity can be enhanced [3-6]. IL-cations are mostly bulky organics and have different substitution patterns surrounding a specific center molecule (e.g. imidazolium). The IL-anion can be small ions like chloride or bulky as $\left[\mathrm{NTf}_{2}\right]^{-}$. The properties of an IL can thus be varied by small adjustments of either IL-cation and IL-anion, e.g. by

[^0]changing the alkyl chain length from $\left[\mathrm{C}_{2} \mathrm{mim}\right]^{+}$to $\left[\mathrm{C}_{10} \mathrm{mim}\right]^{+}$. Moreover, the IL-anion can easily be varied using ion exchange [7].

Besides the usage of pure ILs as solvents, many applications use aqueous IL solutions, which is mirrored in the variety of papers investigating the mutual solubility of water and ILs. For example, ILs are seen as potential co-solvents in liquid extraction processes [8,9]. Therefore, in a first step, the liquid-liquid equilibrium of the binary solvent systems and thereafter the influence of the IL on the solubility of the extract need to be investigated. A second possible benefit of ILs is the shift of biological or enzyme-catalyzed reactions to higher selectivity or yield through incorporating ILs into the aqueous solution, e.g. for introducing of a second (non-)aqueous phase [10,11]. The great number of possible ILs makes it impossible to determine all their capability experimentally, while the usage of thermodynamic models allows promising screening possibilities. Different contributions to this approach have been made so far. Great effort to investigate all kinds of phase equilibria (VLE, LLE and SLE) was made in the group of Domanska [12-15] for many different ILs and solvents, although, in the scope of this research, water as a solvent was under main investigation. Next to the experimental studies, the systems were mainly correlated with $\mathrm{g}^{\mathrm{E}}$
models, such as NRTL. Besides, Domanska et al. [16] were able to quantitatively calculate the LLE of the mixture water + an isoquinolinium IL with PC-SAFT in a semi-predictive manner by adjusting the cross association between water and IL to the activity coefficient of the molecular solvent at infinite dilution.

Likewise, Freire and co-workers investigated numerous different ILs and solvents and their mutual solubility [17-20]. The COSMO-RS method, a model based on unimolecular quantum calculations, was used by the group of Freire et al. [21] to predict the mutual solubility of water and a series of $\left[\mathrm{C}_{\mathrm{n}} \operatorname{mim}\right]\left[\mathrm{NTf}_{2}\right](\mathrm{n}=2-8)$. The results show acceptable agreement with the experimental data for the systems investigated. The compressible lattice model NonRandom Hydrogen Bonding (NRHB) was used by Tsioptsias et al. [22] to model binary water + IL and alcohol + IL systems. However, a temperature-independent interaction parameter was required to almost quantitatively correlate the composition of the two phases over the temperature range investigated. To the best of our knowledge, these are the only two papers issuing the mutual solubility of water and ILs modeled with an equation of state. Classical non-electrolyte PC-SAFT was used by Nann et al. [23] to investigate the LLE of four ternary water-butanol-IL systems. The model needed two binary interaction parameters (a temperaturedependent $k_{i j}(T)$ as well as an $l_{i j}$ for the segment diameter) for successfully modeling the systems for four different ILs. Classical PC-SAFT does not allow an IL-ion-specific modeling approach. Thus, the model used in this study is based on the electrolyte PC-SAFT equation of state (original ePC-SAFT) [24], which is the ion-based extension of PC-SAFT that is able to account for electrolytes. It has shown before to be an excellent choice for the prediction and modeling of different thermodynamic properties for ILs; amongst others, densities of pure ILs and mixtures of ILs and solvents [25,26], gas solubility [27], salt solubility [28], viscosity [29], speed of sound and other derivative properties [30]. Furthermore, as an example for an aqueous system, the influence of ILs on the yield of enzyme-catalyzed reactions was predicted with original ePC-SAFT [11,31]. In addition to these properties, original ePC-SAFT was recently used to qualitatively predict LLE of binary systems of selenium-containing ILs and organic solvents [32]. Although bringing satisfying results for the above-mentioned properties, quantitative representation of LLE of binary mixtures water + ILs was not possible with original ePC-SAFT. There are different opportunities to address this shortcoming, e.g. incorporating the Born term [33] or introducing ion association into the fitting routine of pure-component parameters. The Born term becomes quite important for mixtures with components that have very different dielectric constants. However, in the present work, the influence of the dielectric constant is addressed differently, motivated by the fact that the concentration of IL-ions can be very high in IL-water mixtures. Accordingly, neglecting the fact that the dielectric constant $\varepsilon$ depends strongly on IL concentration in water might be an invalid assumption. Therefore, original ePC-SAFT needs to be extended and the concentration-dependent dielectric constant $\varepsilon$ needs to be incorporated.

The goal of this work is to incorporate the concentrationdependent $\varepsilon$ into original ePC-SAFT (further on simply referred to as ePC-SAFT), to verify the ability of this new approach within ePCSAFT to predict the mutual solubility of water and ILs, and to compare with the results of the original ePC-SAFT model. The ILs under investigation contain the following most common IL-cations $\left[\mathrm{C}_{\mathrm{n}} \mathrm{mim}\right]^{+},\left[\mathrm{C}_{\mathrm{n}} \mathrm{py}\right]^{+},\left[\mathrm{C}_{\mathrm{n}} \mathrm{mpy}\right]^{+},\left[\mathrm{C}_{\mathrm{n}} \mathrm{mpyr}\right]^{+},\left[\mathrm{C}_{4} \mathrm{~m}_{4} \mathrm{py}\right]^{+}$and the ILanions $\left[\mathrm{BF}_{4}\right]^{-},\left[\mathrm{NTF}_{2}\right]^{-},\left[\mathrm{PF}_{6}\right]^{-},[\mathrm{TFO}]^{-}$. The abbreviations are summarized in Table 1. Small adjustments to the IL properties are possible by changing the alkyl chain length ( $\mathrm{C}_{\mathrm{n}}$ ) of the IL-cation or exchanging the IL-anion, and thus the mirroring of this effect with ePC-SAFT was also investigated.

Table 1
Full names of the IL-ions studied in this work.
| IL-ion | Full name |
| :--- | :--- |
| IL-cations |  |
| [ $\left.\mathrm{C}_{\mathrm{n}} \mathrm{mim}\right]^{+}$ | 1-n-alkyl-3-methylimidazolium |
| $\left[\mathrm{C}_{\mathrm{n}} \text { py }\right]^{+}$ | 1-n-alkyl-pyridinium |
| [ $\mathrm{C}_{\mathrm{n}}$ mpy] ${ }^{+}$ | 1-n-alkyl-3-methylpyridinium |
| $\left[\mathrm{C}_{\mathrm{n}} \mathrm{mpyr}\right]^{+}$ | 1-methyl-1-n-alkylpyrrolidinium |
| $\left[\mathrm{C}_{4} \mathrm{~m}_{4} \mathrm{py}\right]^{+}$ | 1-butyl-4-methylpyridinium |
| IL-anions |  |
| $\left[\mathrm{NTf}_{2}\right]^{-}$ | bis (trifluoromethyl-sulfonyl)amide |
| $\left[\mathrm{PF}_{6}\right]^{-}$ | hexafluorophosphate |
| $\left[\mathrm{BF}_{4}\right]^{-}$ | tetrafluoroborate |
| [TFO] ${ }^{-}$ | trifluoromethanesulfonate |


## 2. Modeling and theory

### 2.1. Original ePC-SAFT

In the family of equations of state (EOS), the application of the Statistical Associating Fluid Theory (SAFT) has been used and varied over the last three decades. It was first developed by Chapman et al. in 1989 [34] based on Wertheim's first order perturbation theory [35,36], and successfully applied to many applications. These EOS utilize the residual Helmholtz energy $a^{\text {res }}$, from which other residual properties can be easily derived by partial deviation and, since it is a function of temperature, volume, and mole numbers, it is easily accessible from experiments. One of the advances is Perturbed-Chain SAFT developed by Gross and Sadowski [37], which accounts for the perturbation of a hard chain instead of a hard sphere. The electrolyte extension in the framework of PCSAFT, referred to as original ePC-SAFT [24,38], was used in this work to model the thermodynamic properties for solutions containing ILs. A detailed description of the model can be found in Refs. [24,39], and here only a brief summary is given. The residual Helmholtz energy $a^{\text {res }}$ is calculated as shown in equation (1),

$$
\begin{equation*}
a^{\text {res }}=a^{\text {hc }}+a^{\text {disp }}+a^{\text {assoc }}+a^{\text {ion }} \tag{1}
\end{equation*}
$$

where $a^{h c}$ is the respective energy contribution of the hard chain, $a^{\text {disp }}$ denotes to the dispersion-energy contribution, and the association-energy contribution is calculated by $a^{\text {assoc }}$. The detailed calculation of these contributions can be found in the literature [37]. The contribution caused by ion-ion interactions (charge) to the residual Helmholtz energy was accounted for by the DebyeHückel theory [40] (equation (2)).

$$
\begin{equation*}
a^{\text {ion }}=-\frac{\kappa}{12 \pi \cdot \varepsilon} \sum_{j} \chi_{j} q_{j}^{2} x_{j} \tag{2}
\end{equation*}
$$

Here, $\kappa$ is the inverse Debye screening length, the quantity $\chi_{\mathrm{j}}$ helps simplifying the equation and is presented later on (see equation (10)), $\varepsilon$ is the dielectric constant of the solvent, which Debye and Hückel set to the value of pure water in their original work. The charge $q_{j}$ and the mole fraction $x_{j}$ refer to ionic properties, and thus each IL-anion and IL-cation are treated individually. In the previous work [26], different options for characterizing ILs have been investigated. As a result, ILs were determined as independent ions with three pure-component parameters. These parameters are the segment number $m_{i}^{\text {seg }}$, segment diameter $\sigma_{i}$ and the dispersion-energy parameter between two neighboring segments $u_{i} / k_{B}$, where $k_{B}$ is the Boltzmann constant. For water, additionally the association-energy parameter $\varepsilon^{A i B i}$ and associationvolume parameter $\kappa^{\text {AiBi }}$ were applied. For the calculation of mixtures, mixing and combining rules were used. In this framework,
the Lorenz-Berthelot combining rules were utilized (equation (3) and (4)),

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

where $k_{i j}$ is the binary interaction parameter between two components $i$ and $j$. In this work, the interaction parameters for the respective IL-ion and water were calculated with a linear molecu-lar-weight-dependent correlation for the IL-cation to model the LLE of water-IL systems over a wide temperature range.

### 2.2. Liquid-liquid equilibria

For the calculation of the LLE for water - IL systems the isofugacity criteria were applied according to equation (5),

$$
\begin{equation*}
\varphi_{i}^{I} x_{i}^{I}=\varphi_{i}^{I I} x_{i}^{I I} \tag{5}
\end{equation*}
$$

where $\varphi_{i}$ is the fugacity coefficient of component $i$ in either phase $I$ or II, and $x_{i}$ is the respective mole fraction. Within ePC-SAFT, the IL is characterized by independent IL-ions, and thus the isofugacity criterion was solved for a system containing three components. The fugacity coefficients can be derived from the aforementioned residual Helmholtz energy $a^{\text {res }}$ calculated with ePC-SAFT via the chemical potential $\mu_{i}^{\text {res }}$ as depicted in equation (6).

$$
\begin{equation*}
\ln \varphi_{i}=\frac{\mu_{i}^{r e s}(T, v)}{k T}-\ln Z \tag{6}
\end{equation*}
$$

Here, $Z$ is the compressibility factor.

### 2.3. Development of ePC-SAFT with a concentration-dependent dielectric constant

The incorporation of a concentration-dependent dielectric constant is physically meaningful for systems with high electrolyte concentrations, which is apparent in the investigated binary systems water-IL. The IL-rich phase of the equilibrated solutions consists of IL-ions and still a significant amount of water, while the water-rich phase does practically not contain IL. Additionally, the dielectric constant for common ILs was experimentally found to be $11 \pm 5$ for all ILs that have been studied so far [41], while for water the dielectric constant is about seven times larger at room temperature. Original ePC-SAFT treats the dielectric constant as a fixed value independent of the composition. This can be improved by different possibilities. First, the value of water or of the IL for $\varepsilon$ in the Debye-Hückel term might be used, and for the latter the experimentally found value or unity might be used. A second possibility is to use a compound-dedicated dielectric constant in the respective compound-rich phase, e.g. the dielectric constant of water in the water-rich phase and $\varepsilon_{\text {IL }}$ in the IL-rich phase. The third, most viable but also most complex option is to vary the dielectric constant with the IL-ion concentration in the respective phase. All these approaches are discussed in the results section. For a concentrationdependent dielectric constant, different aspects need to be discussed. First, a relationship between the varying IL-ion concentrations and $\varepsilon$ in the coexisting phases is needed. Koeberg et al. [42] showed that the dielectric constant of $\left[\mathrm{C}_{4} \mathrm{mim}\right]\left[\mathrm{BF}_{4}\right]$ exhibits a linear dependence on the volume fraction similar to other liquid mixtures. It was shown that this dependency is also valid in mole and weight fractions (see Fig. 1). Unfortunately, this data is the only available source for the concentration-dependency of the dielectric

![](https://cdn.mathpix.com/cropped/987fc63c-d692-452c-bca4-ac1ea8a573f8-3.jpg?height=494&width=609&top_left_y=223&top_left_x=1179)
Fig. 1. Dependency of the relative dielectric constant $\varepsilon_{\boldsymbol{r}}$ on the IL mole fraction. Experimental data from Koeberg et al. [42], converted to mole fractions in this work.

constant in water - IL systems. The possibility of a different behavior for other ILs than $\left[\mathrm{C}_{4} \mathrm{mim}\right]\left[\mathrm{BF}_{4}\right]$ is still given. However, in the present work, this dependence is seen as a suitable approach for the new development within ePC-SAFT. Thus, a simple mixing rule for the dielectric constant in dependence of the mole fraction as seen in equation (7) is introduced. This simple mixing rule can later on easily be changed for the application to other electrolyte systems

$$
\begin{equation*}
\varepsilon=\varepsilon_{\text {water }} \cdot X_{\text {water }}+\varepsilon_{I L} \cdot\left(1-X_{\text {water }}\right) \tag{7}
\end{equation*}
$$

where $\varepsilon_{\text {water }}$ (the value for pure water) was inherited from the original ePC-SAFT from Cameretti et al. while $\varepsilon_{I L}$ was set to 11 independent of IL and temperature, respectively. This linear correlation for the dielectric constant was then incorporated into the calculation of the ionic contribution to the residual Helmholtz energy $a^{i o n}$. The variation to $a^{i o n}$ is shown in equation (8).

$$
\begin{equation*}
a^{\text {ion }}=-\frac{\kappa(\varepsilon(x))}{12 \pi k_{B} T \cdot \varepsilon(x)} \sum_{j} q_{j}^{2} \cdot x_{j} \cdot \chi_{j}(\varepsilon(x)) \tag{8}
\end{equation*}
$$

As can be seen from equation (8), the concentration-dependent dielectric constant also influences the inverse Debye screening length that consequently affects the quantity $\chi_{i}$ (equation (9) and (10))

$$
\begin{equation*}
\kappa=\sqrt{\frac{\rho_{N} e^{2}}{k_{B} T \varepsilon(x)} \cdot \sum_{i} z_{i}^{2} \cdot x_{i}(\varepsilon(x))} \tag{9}
\end{equation*}
$$

$\chi_{i}=\frac{3}{\left(\kappa(\varepsilon(x)) a_{i}\right)^{3}} \cdot\left[\frac{3}{2}+\ln \left(1+\kappa(\varepsilon(x)) a_{i}\right)-2\left(1+\kappa(\varepsilon(x)) a_{i}\right)\right.$

$$
\begin{equation*}
\left.+\frac{1}{2}\left(1+\kappa(\varepsilon(x)) a_{i}\right)^{2}\right] \tag{10}
\end{equation*}
$$

Additionally, for the calculation of the chemical potential, the derivation of $a^{\text {ion }}$ with respect to the mole fraction of component $n$ needs to be calculated according to equation (11).

$$
\begin{align*}
\frac{d a^{i o n}}{d x_{n}}= & -\frac{1}{12 \pi \mathrm{k}_{\mathrm{B}} T} \cdot\left[\frac{d \kappa}{d x_{n}} \cdot \frac{1}{\varepsilon} \cdot \sum_{i} z_{i}^{2} \cdot x_{i} \cdot \chi_{i}+\kappa \cdot \frac{d \varepsilon^{-1}}{d x_{n}} \cdot \sum_{i} z_{i}^{2} \cdot x_{i} \cdot \chi_{i}\right. \\
& \left.+\frac{\kappa}{\varepsilon}\left[\sum_{i} z_{i}^{2} \cdot \frac{d x_{i}}{d x_{n}} \cdot \chi_{i}+\sum_{i} z_{i}^{2} \cdot x_{i} \cdot \frac{d \chi_{i}}{d x_{n}}\right]\right] \tag{11}
\end{align*}
$$

Furthermore, a distinction of cases for the derivation of the
dielectric constant itself has to be implemented (equation (12)).

$$
\frac{d \varepsilon}{d x_{n}}=\left\{\begin{array}{l}\varepsilon_{0} \cdot\left(\varepsilon_{1}\right), \text { if } n=1  \tag{12}\\ \varepsilon_{0} \cdot\left(\varepsilon_{2}\right), \text { if } n \neq 1\end{array}\right.
$$

Here, the $n=1$ denotes to water. Again, the inverse Debye screening length $\kappa$ and the quantity $\chi_{i}$ are influenced (equation (13) to (14)).

$$
\begin{align*}
\frac{d \kappa}{d x_{n}}= & \left(\frac{\rho_{N} e^{2}}{k_{B} T}\right)^{\frac{1}{2}} \cdot\left[-\frac{1}{2}[\varepsilon]^{-\frac{3}{2}} \cdot \frac{d \varepsilon}{d x_{n}}\left[\sum_{i} z_{i}^{2} \cdot x_{i}\right]^{\frac{1}{2}}\right. \\
& \left.+\frac{1}{2 \sqrt{\varepsilon}}\left[\sum_{i} z_{i}^{2} \cdot x_{i}\right]^{-\frac{1}{2}} \cdot\left[\sum_{i} z_{i}^{2} \cdot \frac{d x_{i}}{d x_{n}}\right]\right]  \tag{13}\\
\frac{d \chi_{i}}{d x_{n}}= & 3 \cdot \frac{d \kappa}{d x_{n}} \cdot\left[-\frac{\chi_{i}}{\kappa}+\frac{a_{i}}{\left(\kappa \cdot a_{i}\right)^{3}} \cdot\left[\frac{1}{1+\kappa a_{i}}-2+\left(1+\kappa a_{i}\right)\right]\right] \tag{14}
\end{align*}
$$

The resulting extension is henceforth addressed as ePC-SAFT. To illustrate the influence of the new approach, Fig. 2 shows the derivation of $a^{\text {ion }}$ with respect to the mole fraction of the IL (Fig. 2 a)). Fig. 2 b) depicts the resulting fugacity coefficients of an IL-ion for a constant dielectric constant ( $\varepsilon_{\text {water }}$ or $\varepsilon_{\mathrm{IL}}$ ), respectively, and with ePC-SAFT ( $\varepsilon\left(\mathrm{X}_{\mathrm{II}}\right)$ ). Concrete $a^{\text {ion }}$ values depend strongly on using either $\varepsilon_{\text {water }}$ or $\varepsilon_{\text {IL }}$. Especially for low IL concentrations, the dielectric constant of pure IL would cause a huge deviation in $a^{i o n}$. Subsequently, in the vice versa scenario ( $\varepsilon_{\text {water }}$ and high IL concentration) the derivative of $a^{\text {ion }}$ would be largely overestimated. A similar behavior is visible for the fugacity coefficient, but the difference becomes more apparent for high IL concentrations. Consequently, for the pure components ( $\mathrm{x}_{\mathrm{IL}}=1$ and $\mathrm{x}_{\mathrm{IL}}=0$ ), both properties, the derivative of $a^{i o n}$ and the fugacity coefficient, have to take the values for the case that $\varepsilon_{\text {water }}$ or $\varepsilon_{\mathrm{IL}}$ are used in the modeling.

## 3. Results and discussion

Pure-component parameters are used as an input to calculate the residual Helmholtz energy within the ePC-SAFT framework. They are needed to physically meaningful characterize the compounds. The pure-component parameters are depicted in Table 2 for water and Table 3 for the IL-ions. The pure-component parameters for IL-ions have been fitted to liquid density data at varying temperatures and pressures in our previous work with a
![](https://cdn.mathpix.com/cropped/987fc63c-d692-452c-bca4-ac1ea8a573f8-4.jpg?height=582&width=1416&top_left_y=1901&top_left_x=346)

Table 2
Pure-component parameters of water used for calculation of LLE [43].
| Parameter | Unit | value |
| :--- | :--- | :--- |
| Segment number $m_{\text {seg }}$ | [-] | 1.2047 |
| Segment diameter $\sigma_{i}$ | [Å] | 2.7927 |
| Dispersion-energy parameter $u_{i j} / k_{B}$ | [K] | 353.9449 |
| Association-energy parameter $\varepsilon^{A i B i} / k_{B}$ | [K] | 2425.6714 |
| Association-volume parameter $\kappa^{A i B i}$ | [-] | 0.04509 |
| Association scheme | [-] | 2B |


dielectric constant set to unity [26,27].
The experimental data used to validate the correlative and predictive abilities of original ePC-SAFT and of ePC-SAFT are summarized in Table 4. The systems under study are water - IL at 1 bar and a temperature ranging from 278.15 to 364.55 K . The ILs are comprised of the most common ions (see also Tables 3 and 5).

### 3.1. Verification of the previous pure-component parameters for ILions

In predecessor works on ILs with original ePC-SAFT, the dielectric constant was set to unity for fitting of the purecomponent parameters and also for the calculation of different thermodynamic properties as cited in the introduction. In contrast to this assumption, the dielectric constant of the investigated ILs varies around an average of 11 , and this value was used in this work for the calculation of the LLE. This new treatment would normally require a re-parametrization of the IL-ion pure-component parameters with a dielectric constant of 11 . As a test, the IL-ions $\left[\mathrm{C}_{\mathrm{n}} \mathrm{mim}\right]^{+}$and $\left[\mathrm{NTf}_{2}\right]^{-}$have been refitted to the liquid densities of the pure ILs and the impact of a re-parametrization on the LLE with water was investigated. The modeling results did not show significant differences for LLE water - IL systems. Thus, the IL-ion purecomponent parameters obtained in our previous work are still valid and there is no need to have a full re-parametrization of the IL-ion parameters to the new dielectric constant of 11 .

### 3.2. Modeling LLE using different treatments of the dielectric constant

In this work, four options to account for the dielectric constant in the ePC-SAFT framework have been investigated. Original ePCSAFT uses a dielectric constant with the value of water $\left(\varepsilon_{\mathrm{r}}=80\right)$

Table 3
Pure-component parameters of IL-ions used in this work [26,27].
| IL-ion | $M_{w}[\mathrm{~g} / \mathrm{mol}]$ | $m_{\text {seg }}$ | $\sigma_{i}[\AA]$ | $u_{i j} / k_{B}[K]$ |
| :--- | :--- | :--- | :--- | :--- |
| IL-anions |  |  |  |  |
| $\left[\mathrm{NTf}_{2}\right]^{-}$ | 280.145 | 6.0103 | 3.7469 | 375.6529 |
| $\left[\mathrm{BF}_{4}\right]^{-}$ | 86.809 | 3.8227 | 3.5088 | 496.1176 |
| $\left[\mathrm{PF}_{6}\right]^{-}$ | 144.973 | 4.2771 | 3.5889 | 492.2835 |
| [TFO] ${ }^{-}$ | 149.070 | 3.7432 | 3.8771 | 509.3113 |
| IL-cations |  |  |  |  |
| [ $\left.\mathrm{C}_{2} \mathrm{mim}\right]^{+}$ | 111.168 | 1.4872 | 3.5926 | 206.4924 |
| [ $\left.\mathrm{C}_{3} \mathrm{mim}\right]^{+}$ | 125.190 | 1.9857 | 3.6154 | 212.3200 |
| [ $\left.\mathrm{C}_{4} \mathrm{mim}\right]^{+}$ | 139.221 | 2.4805 | 3.6371 | 218.1441 |
| [ $\left.\mathrm{C}_{5} \mathrm{mim}\right]^{+}$ | 153.240 | 2.9226 | 3.6575 | 224.1573 |
| $\left[\mathrm{C}_{6} \mathrm{mim}\right]^{+}$ | 167.275 | 3.4131 | 3.6781 | 230.0000 |
| $\left[\mathrm{C}_{7} \mathrm{mim}\right]^{+}$ | 181.300 | 3.8598 | 3.6996 | 235.9986 |
| [ $\left.\mathrm{C}_{8} \mathrm{mim}\right]^{+}$ | 195.328 | 4.2977 | 3.7187 | 242.0000 |
| [ $\mathrm{C}_{10}$ mim] ${ }^{+}$ | 223.380 | 5.2653 | 3.7627 | 253.7600 |
| [ $\left.\mathrm{C}_{3} \mathrm{mpyr}\right]^{+}$ | 128.235 | 2.2941 | 3.6082 | 167.9026 |
| [ $\mathrm{C}_{4}$ mpyr] ${ }^{+}$ | 142.257 | 2.6364 | 3.7448 | 273.8825 |
| $\left[\mathrm{C}_{6} \mathrm{mpyr}\right]^{+}$ | 170.310 | 3.3212 | 4.0179 | 485.9104 |
| $\left[\mathrm{C}_{4} \mathrm{~m}_{4} \mathrm{py}\right]^{+}$ | 150.242 | 2.5504 | 3.8046 | 277.1743 |
| [ $\mathrm{C}_{4}$ py] ${ }^{+}$ | 136.210 | 3.1335 | 3.3131 | 235.6465 |
| [ $\mathrm{C}_{6}$ py] ${ }^{+}$ | 164.263 | 4.4916 | 3.3243 | 246.0532 |
| [ $\mathrm{C}_{8}$ py] ${ }^{+}$ | 192.316 | 5.8496 | 3.3354 | 256.4599 |
| [ $\mathrm{C}_{3}$ mpy] ${ }^{+}$ | 136.214 | 2.7283 | 3.4306 | 161.5381 |
| [C4mpy] ${ }^{+}$ | 150.242 | 2.5504 | 3.8046 | 277.1743 |


Table 4
Investigated LLE of binary water - IL systems at 1 bar used to validate the new model ePC-SAFT. NP is the number of data points available.
| IL | $\mathrm{T}_{\text {min }}$ | $\mathrm{T}_{\text {max }}$ | NP | Reference |
| :--- | :--- | :--- | :--- | :--- |
| [ $\mathrm{C}_{2}$ mim] [ $\mathrm{NTf}_{2}$ ] | 285.15 | 361.17 | 59 | [17,44-46] |
| [ $\mathrm{C}_{3}$ mim][ $\mathrm{NTf}_{2}$ ] | 288.15 | 318.15 | 14 | [17] |
| [ $\mathrm{C}_{4}$ mim] [ $\mathrm{NTf}_{2}$ ] | 288.15 | 360.05 | 51 | [17,44-48] |
| [ $\mathrm{C}_{5}$ mim] [ $\mathrm{NTf}_{2}$ ] | 288.15 | 318.15 | 14 | [17] |
| [ $\mathrm{C}_{6}$ mim] [ $\mathrm{NTf}_{2}$ ] | 288.15 | 333.15 | 44 | [17,46,49,50] |
| [ $\mathrm{C}_{7}$ mim] [ $\mathrm{NTf}_{2}$ ] | 288.15 | 318.15 | 14 | [17] |
| [ $\mathrm{C}_{8}$ mim][ $\mathrm{NTf}_{2}$ ] | 288.15 | 318.15 | 26 | [17,48,50] |
| [ $\mathrm{C}_{3}$ mima][ $\mathrm{PF}_{6}$ ] | 288.15 | 318.15 | 14 | [51] |
| [ $\mathrm{C}_{4}$ mima][ $\mathrm{PF}_{6}$ ] | 288.15 | 363.05 | 78 | [18,20,44,45,47,52-54] |
| [ $\mathrm{C}_{6}$ mima][ $\mathrm{PF}_{6}$ ] | 288.15 | 353.15 | 45 | [18,53] |
| [ $\mathrm{C}_{8}$ mima][ $\mathrm{PF}_{6}$ ] | 278.15 | 364.55 | 47 | [18,52,53] |
| [ $\mathrm{C}_{3}$ mpyr][ $\mathrm{NTf}_{2}$ ] | 278.15 | 343.15 | 46 | [18,20,55] |
| [ $\mathrm{C}_{4}$ mpyr][ $\mathrm{NTf}_{2}$ ] | 288.15 | 343.15 | 18 | [18,20,49] |
| [ $\mathrm{C}_{6}$ mpyr][ $\mathrm{NTf}_{2}$ ] | 288.15 | 338.15 | 12 | [49] |
| $\left[\mathrm{C}_{4} \mathrm{~m}_{4} \mathrm{py}\right]\left[\mathrm{NTf}_{2}\right]$ | 288.15 | 358.90 | 35 | [14,19,56] |
| [ $\mathrm{C}_{6} \mathrm{mim}$ ][ $\mathrm{BF}_{4}$ ] | 278.15 | 333.15 | 31 | [57,58] |
| [ $\mathrm{C}_{8}$ mima][ $\mathrm{BF}_{4}$ ] | 278.15 | 340.15 | 31 | [52,57] |
| [ $\mathrm{C}_{10} \mathrm{mim}$ ] $\left[\mathrm{BF}_{4}\right]$ | 278.15 | 325.15 | 24 | [59] |
| [ $\mathrm{C}_{8}$ mim][TFO] | 293.15 | 313.15 | 10 | [48] |
| [ $\mathrm{C}_{4}$ py][ $\mathrm{NTf}_{2}$ ] | 288.15 | 343.67 | 35 | [19,48,56] |
| [ $\mathrm{C}_{6}$ py][ $\mathrm{NTf}_{2}$ ] | 288.15 | 318.15 | 16 | [19,50] |
| [ $\mathrm{C}_{8}$ py][ $\mathrm{NTf}_{2}$ ] | 288.15 | 318.15 | 14 | [19] |
| [ $\mathrm{C}_{3}$ mpy] [ $\mathrm{PF}_{6}$ ] | 288.15 | 318.15 | 35 | [18,51] |
| [ $\mathrm{C}_{3}$ mpy][ $\mathrm{NTf}_{2}$ ] | 288.15 | 318.15 | 14 | [18] |
| [ $\mathrm{C}_{4}$ mpy][ $\mathrm{NTf}_{2}$ ] | 288.15 | 318.15 | 16 | [50,56] |


as intended by Debye and Hückel. For the binary mixtures of water and ILs, the use of a dielectric constant similar to the pure IL $\left(\varepsilon_{\mathrm{r}}=11\right)$ is the second option. A third option is to vary the value of the dielectric constant with the respective species-rich phase. Therefore, the water-rich phase was assigned with the value of pure water and the IL-rich phase with the value of the pure IL. The last option is to account for the concentration of the IL in both phases with ePC-SAFT. The difference of the four options is exemplarily discussed based on the system water $-\left[\mathrm{C}_{4} \mathrm{mim}\right]\left[\mathrm{NTf}_{2}\right]$. The results for the LLE are depicted in Fig. 3. The benchmark for the comparison is the prediction result of the original ePC-SAFT model (2014). The binary interaction parameter $k_{i j}$ was set to zero for all four options. The solubility of water in the IL is highly underestimated with an

Table 5
Comparison of original ePC-SAFT and ePC-SAFT for modeling the LLE of water - IL systems at 1 bar in terms of AAD and ARD for the IL-rich phase using the parameters from Tables 2, 3 and 6. AAD and ARD correspond to the same temperature ranges as shown in Table 4.
| IL | original ePC-SAFT |  | ePC-SAFT including concentrationdependent $\varepsilon$ |  |
| :--- | :--- | :--- | :--- | :--- |
|  | AAD [-] | ARD [\%] | AAD [-] | ARD [\%] |
| [ $\mathrm{C}_{2}$ mim] [ $\mathrm{NTf}_{2}$ ] | 0.1488 | 80.51 | 0.0295 | 5.09 |
| [ $\mathrm{C}_{3}$ mim] [ $\mathrm{NTf}_{2}$ ] | 0.1357 | 15.89 | 0.0178 | 2.54 |
| [ $\mathrm{C}_{4}$ mim][ $\mathrm{NTf}_{2}$ ] | 0.1120 | 13.16 | 0.0161 | 2.23 |
| [ $\mathrm{C}_{5}$ mim][ $\mathrm{NTf}_{2}$ ] | 0.1037 | 11.81 | 0.0325 | 4.43 |
| [ $\mathrm{C}_{6}$ mim] [ $\mathrm{NTf}_{2}$ ] | 0.0708 | 8.14 | 0.0234 | 3.09 |
| [ $\mathrm{C}_{7}$ mim][ $\mathrm{NTf}_{2}$ ] | 0.0775 | 8.81 | 0.0105 | 1.34 |
| [ $\mathrm{C}_{8}$ mim][ $\mathrm{NTf}_{2}$ ] | 0.0428 | 4.97 | 0.0013 | 0.16 |
| [ $\mathrm{C}_{3}$ mima][ $\mathrm{PF}_{6}$ ] | 0.0969 | 12.14 | 0.0287 | 4.26 |
| [ $\mathrm{C}_{4}$ mim] ] $\mathrm{PF}_{6}$ ] | 0.1128 | 13.52 | 0.0185 | 2.61 |
| [ $\mathrm{C}_{6} \mathrm{mim}$ ] $\left[\mathrm{PF}_{6}\right]$ | 0.1229 | 14.05 | 0.0180 | 2.40 |
| [ $\mathrm{C}_{8}$ mim] ] $\mathrm{PF}_{6}$ ] | 0.0654 | 11.31 | 0.0521 | 6.76 |
| [ $\mathrm{C}_{3} \mathrm{mpyr}$ ][ $\mathrm{NTf}_{2}$ ] | 0.0664 | 8.31 | 0.0279 | 3.84 |
| [ $\mathrm{C}_{4}$ mpyr][ $\mathrm{NTf}_{2}$ ] | 0.0753 | 13.70 | 0.0269 | 3.76 |
| [ $\mathrm{C}_{6}$ mpyr] $\left[\mathrm{NTf}_{2}\right]$ | 0.0695 | 8.42 | 0.0477 | 6.86 |
| [ $\mathrm{C}_{4} \mathrm{~m}_{4} \mathrm{py}$ ][ $\mathrm{NTf}_{2}$ ] | 0.0572 | 7.03 | 0.0066 | 0.91 |
| [ $\mathrm{C}_{6} \mathrm{mim}$ ][ $\mathrm{BF}_{4}$ ] | 0.3980 | 57.75 | 0.1732 | 42.85 |
| $\left[\mathrm{C}_{8} \mathrm{mim}\right]\left[\mathrm{BF}_{4}\right]$ | 0.3966 | 140.82 | 0.0179 | 12.36 |
| $\left[\mathrm{C}_{10} \mathrm{mim}\right]\left[\mathrm{BF}_{4}\right]$ | 0.4399 | 71.10 | 0.4687 | 72.80 |
| [ $\mathrm{C}_{8}$ mim][TFO] | 0.0896 | 12.21 | 0.0047 | 0.77 |
| [ $\mathrm{C}_{4}$ py][ $\mathrm{NTf}_{2}$ ] | 0.0561 | 6.85 | 0.0146 | 2.02 |
| [ $\mathrm{C}_{6}$ py][ $\mathrm{NTf}_{2}$ ] | 0.0476 | 5.72 | 0.0113 | 1.49 |
| [ $\mathrm{C}_{8}$ py][ $\mathrm{NTf}_{2}$ ] | 0.0325 | 4.12 | 0.0151 | 1.94 |
| [ $\mathrm{C}_{3}$ mpy][ $\mathrm{NTf}_{2}$ ] | 0.0890 | 10.58 | 0.0166 | 2.24 |
| [C4mpy][NTf ${ }_{2}$ ] | 0.0412 | 5.36 | 0.0263 | 3.52 |
| Sum | 2.9482 | 22.76 | 1.1059 | 7.93 |


Fig. 3. Comparison of the different options to treat the dielectric constant within the original ePC-SAFT framework for the system water - $\left[\mathrm{C}_{4} \mathrm{mim}\right]\left[\mathrm{NTf}_{2}\right]$ at 1 bar. Symbols are experimental data, lines are predictions: Blue: original ePC-SAFT ( $\varepsilon_{\mathrm{r}}=80$ ); orange: original ePC-SAFT ( $\varepsilon_{\mathrm{r}}=11$ ) and original ePC-SAFT ( $\varepsilon=$ phase-dependent); and green: ePC-SAFT. All $k_{i j}$ have been set to zero for the predictions. No major difference in prediction of IL solubility in water. Pure-component parameters are depicted in Tables 2 and 3

ARD of $22 \%$ using original ePC-SAFT. A dielectric constant of the IL for the binary mixture leads to an improved representation of the experimental data ( $8.75 \%$ ARD). However, the temperaturedependency of the binodal curve is not reproduced correctly. For the third option, the results for the solubility of water in the IL and the temperature-dependency only marginally change compared to the second option with a constant dielectric constant of the IL ( 8.76 $\%$ ARD). A concentration-dependent dielectric constant wihin ePCSAFT, in contrast to the prior methods, overestimates the solubility slightly, though giving the best reproduction of the
temperature-dependency and a still acceptable deviation to the experimental data ( 12.27 \%ARD). The results are listed in Table 5.

### 3.3. Modeling LLE with original ePC-SAFT and ePC-SAFT including the concentration-dependent dielectric constant

### 3.3.1. Prediction

The predictive capability of original ePC-SAFT and the new approach have already been compared for the system water [ $\left.\mathrm{C}_{4} \mathrm{mim}\right]\left[\mathrm{NTf}_{2}\right]$ in Fig. 3 over a wide temperature range. One major advantage of the new approach in comparison to original ePC-SAFT lies within the reproducibility of the influence of the alkyl chain length of the IL-cation. The comparison is exemplarily shown for the LLE of systems water and $\left[\mathrm{NTf}_{2}\right]$-ILs containing $\left[\mathrm{C}_{\mathrm{n}} \mathrm{mim}\right]^{+}$-cations ( $\mathrm{n}=2-8$ ) in Fig. 4. For original ePC-SAFT ( $\varepsilon_{\mathrm{r}}=80$ ), no considerable differentiation in the prediction of the water solubility for varying chain length is visible. In contrast, ePC-SAFT predicts the influence qualitatively, while giving a small overall overestimation.

### 3.3.2. Correlation

To increase the accuracy of ePC-SAFT (Fig. 4), a binary interaction parameter $k_{i j}$ (equation (4)) was applied. With this parameter, it is usually possible to model the desired properties in agreement with the experimental data. However, as can be seen in Fig. 5 for the system [ $\mathrm{C}_{4} \mathrm{mim}$ ][ $\mathrm{NTf}_{2}$ ], this approach did not yield good results for the IL-rich phase using original ePC-SAFT. In contrast, the new approach within ePC-SAFT allowed an almost quantitative correlation of the experimental data by using just one $k_{i j}$ between IL-ion and water.

To further compare the results obtained with ePC-SAFT including a concentration-dependent dielectric constant with those of original ePC-SAFT as well as with the experimental data, the absolute relative deviation (ARD) and the average absolute deviation (AAD) were calculated according to equations (15) and (16).

$$
\begin{equation*}
\mathrm{ARD}=\frac{1}{N P} \sum_{i=1}^{N P}\left|1-\frac{\chi_{I L}^{\text {calc }}}{\chi_{I L}^{\text {exp }}}\right| \cdot 100 \% \tag{15}
\end{equation*}
$$

$$
\begin{equation*}
\mathrm{AAD}=\frac{1}{N P} \sum_{i=1}^{N P}\left|x_{I L}^{\text {calc }}-x_{I L}^{\text {exp }}\right| \tag{16}
\end{equation*}
$$

Here, NP is the number of data points, $x_{I L}^{\text {exp }}$ is the IL mole fraction from experimental data, and $x_{I L}^{\text {calc }}$ is the IL mole fraction calculated with ePC-SAFT. Most of the ILs, and generally the ILs investigated in this work, have a very poor to almost negligible solubility in water.

![](https://cdn.mathpix.com/cropped/987fc63c-d692-452c-bca4-ac1ea8a573f8-6.jpg?height=451&width=609&top_left_y=1955&top_left_x=279)
Fig. 4. Predictions of the influence of the alkyl chain length of the IL-cation on the LLE of water - $\left[\mathrm{C}_{\mathrm{n}} \operatorname{mim}\right]\left[\mathrm{NTf}_{2}\right](\mathrm{n}=2-8)$ at ambient pressure and 298.15 K . Higher alkyl chain length reduces water solubility in the IL. Blue line is the original ePC-SAFT $\left(\varepsilon_{\mathrm{r}}=80\right)$ prediction, green bars are ePC-SAFT predictions and grey bars are experimental data. Note, that the y -axis develops from top $\left(x_{w}=0\right)$ to bottom $\left(x_{w}=1\right)$.

![](https://cdn.mathpix.com/cropped/987fc63c-d692-452c-bca4-ac1ea8a573f8-6.jpg?height=516&width=608&top_left_y=223&top_left_x=1212)
Fig. 5. Comparison of original ePC-SAFT and ePC-SAFT including a concentrationdependent dielectric constant to model the system water - [ $\mathrm{C}_{4} \mathrm{mim}$ ][NTf ${ }_{2}$ ] at 1 bar. Symbols are experimental data and lines are modeling results: blue, solid: original ePC-SAFT $\left(\varepsilon_{\mathrm{r}}=80\right) k_{i j}=0$; blue, dashed: original ePC-SAFT $\left(\varepsilon_{\mathrm{r}}=80\right)$ with $k_{i j}$; green, solid: ePC-SAFT (new appraoch) $k_{i j}=0$; and green, dashed: ePC-SAFT (new appraoch) with $k_{i j}$.

Thus, their solubility in the IL-rich phase is overall more interesting. Consequently, the deviations were given for the IL-rich phase, i.e. for the solubility of water in the IL phase. To model the LLE of water - IL systems, binary interaction parameters between the IL-ions and water were used as listed in Table 6. For both models, the interaction parameters between IL-cations and water were fitted assuming a linear correlation of $k_{i j}$ with the molecular weight of the IL, in accordance with the procedure for the pure-component parameters of the IL-ions. For the IL-anions, only one $k_{i j}$ value was adjusted for each pair water-ion. As already discussed earlier, the use of temperature-independent $k_{i j}$ for original ePC-SAFT leads to an underestimation of the upper critical solution temperature. In contrast, ePC-SAFT allows a much better representation of the temperature dependency of the experimental LLE without the need of a temperature-dependent $k_{i j}$.

For most of the systems investigated, the experimental data of different authors are equivalent within the error of the methods applied. However, some experimental data do vary significantly. Here, for the purpose of clarity, the ARD and AAD values were calculated from the experimental data with good consistency.

The developed model ePC-SAFT is suitable for calculating various LLE of the systems containing water and different ILs. The most common IL-cations under study are comprised of the imidazolium center molecule and varying alkyl chain length. This

Table 6
Correlation of binary interaction parameters $k_{i j}$ between water and the IL-ions used to model LLE systems between water and ILs with ePC-SAFT and original ePC-SAFT $\left(\varepsilon_{\mathrm{r}}=80\right)$. For the ILs with changing IL-cation chain length, a correlation with the molecular weight $\left(M w_{I L}\right)$ of the IL with $\left[\mathrm{NTf}_{2}\right]^{-}$is used for calculating $k_{i j}$-values.
| IL-Ions | ePC-SAFT including concentration-dependent $\varepsilon$ | original ePC-SAFT |
| :--- | :--- | :--- |
|  | $k_{i j}$ with water | $k_{i j}$ with water |
| IL-anions |  |  |
| [ $\left.\mathrm{NTf}_{2}\right]^{-}$ | 0.11 | -0.03 |
| $\left[\mathrm{BF}_{4}\right]^{-}$ | -0.17 | -0.01 |
| $\left[\mathrm{PF}_{6}\right]^{-}$ | 0.14 | -0.01 |
| [TFO] ${ }^{-}$ | -0.065 | -0.065 |
| IL-cations |  |  |
| [ $\left.\mathrm{C}_{\mathrm{n}} \mathrm{mim}\right]^{+}$ | $0.0032 \cdot M w_{I L}-1.3054$ | $-0.0006 \cdot M w_{I L}+0.2325$ |
| $\left[\mathrm{C}_{\mathrm{n}} \text { py }\right]^{+}$ | $0.0007 \cdot M w_{I L}-0.1968$ | $-0.0002 \cdot M w_{I L}+0.0292$ |
| [ $\mathrm{C}_{\mathrm{n}}$ mpy] ${ }^{+}$ | $-0.0043 \cdot M w_{I L}-1.6408$ | $-0.0014 \cdot M w_{I L}+0.6636$ |
| $\left[\mathrm{C}_{\mathrm{n}} \mathrm{mpyr}\right]^{+}$ | $-0.0005 \cdot M w_{I L}+0.2941$ | $-0.0016 \cdot M w_{I L}+0.7133$ |
| $\left[\mathrm{C}_{4} \mathrm{~m}_{4} \mathrm{py}\right]^{+}$ | 0.0718 | -0.05 |


![](https://cdn.mathpix.com/cropped/987fc63c-d692-452c-bca4-ac1ea8a573f8-7.jpg?height=559&width=1419&top_left_y=213&top_left_x=309)
Fig. 6. (a) Liquid-liquid equilibria of the systems water $-\left[\mathrm{C}_{\mathrm{n}} \operatorname{mim}\right]\left[\mathrm{NTf}_{2}\right](\mathrm{n}=2,4,6,8)$ at 1 bar. IL solubility in water is shown only exemplarily for $\left[\mathrm{C}_{8} \mathrm{mim}\right]\left[\mathrm{NTf}_{2}\right]$. Symbols depict experimental data (black: $C_{2}$, grey: $C_{4}$, orange: $C_{6}$, green: $C_{8}$ ), lines are ePC-SAFT calculations including a concentration-dependent dielectric constant. (b) Liquid-liquid equilibria of the systems water $-\left[\mathrm{C}_{8} \mathrm{mim}\right][\mathrm{x}]\left(\mathrm{x}=\left[\mathrm{BF}_{4}\right],[\mathrm{TFO}],\left[\mathrm{NTf}_{2}\right]\right.$, and $\left.\left[\mathrm{PF}_{6}\right]\right)$ at 1 bar. Symbols depict experimental data (black: $\left[\mathrm{NTf}_{2}\right]$, grey: $\left[\mathrm{BF}_{4}\right]$, orange: [TFO], green: [ $\left.\mathrm{PF}_{6}\right]$ ), lines are $\mathrm{ePC}^{-}$ SAFT calculations including concentration-dependent dielectric constant.

variation of the chain length and the resulting change in the LLE is generally covered with ePC-SAFT as depicted in Fig. 6a) for the system water $+\left[\mathrm{C}_{\mathrm{n}} \operatorname{mim}\right]\left[\mathrm{NTf}_{2}\right](\mathrm{n}=2,4,6,8)$ at one bar over the whole temperature range of existing data. The influence of the substitute IL-anion is likewise modeled within the new approach. Results for various IL-anions ( $\left[\mathrm{BF}_{4}\right]^{-},[\mathrm{TFO}]^{-},\left[\mathrm{NTf}_{2}\right]^{-}$, and $\left[\mathrm{PF}_{6}\right]^{-}$) with a fixed IL-cation ( $\left[\mathrm{C}_{8} \mathrm{mim}\right]^{+}$) are shown in Fig. 6b).

The modeling results follow the correct temperature dependency. The water-rich phase is presented only exemplarily for $\left[\mathrm{C}_{8} \mathrm{mim}\right]\left[\mathrm{NTf}_{2}\right]$ since all ILs are almost insoluble in water. This behavior is also predicted correctly with ePC-SAFT. To give a further insight in the accuracy of ePC-SAFT, the solubility of water in the ILrich phase of the LLE water + IL at 298 K and 1 bar is shown in Fig. 7. The overall modeling result using ePC-SAFT is extremely satisfactory considering the fact that only one $k_{i j}$ value for the pairs ILanion + water was used (Table 6) and a chain-length dependent correlation for the $k_{i j}$ of the pair IL-cation + water. Significant deviations between experimental data and ePC-SAFT were found only for the system water - $\left[\mathrm{C}_{\mathrm{n}} \operatorname{mim}\right]\left[\mathrm{BF}_{4}\right]$. The experimental data of these systems are unexpected since the water solubility in the ILrich phase increased in the order $\mathrm{C}_{10}<\mathrm{C}_{6}<\mathrm{C}_{8}$. This experimentally observed non-linear dependence on the chain length could not be modeled with both, original ePC-SAFT (not shown) and ePCSAFT.

![](https://cdn.mathpix.com/cropped/987fc63c-d692-452c-bca4-ac1ea8a573f8-7.jpg?height=557&width=707&top_left_y=1881&top_left_x=198)
Fig. 7. Modeling results of water solubility in IL-rich phase using ePC-SAFT for systems water + IL for 17 ILs at 1 bar and 298.15 K . For the system water - $\left[\mathrm{C}_{3} \mathrm{mpy}\right]\left[\mathrm{NTf}_{2}\right]$ the temperature was 303.15 K . Green bars are ePC-SAFT predictions including concentration-dependent dielectric constant and grey bars are experimental data.

## 4. Conclusion

In this work, ePC-SAFT is introduced as a new modeling approach for systems containing ILs, where the dielectric constant in the ionic contribution to the residual Helmholtz energy was treated as the concentration-dependent property. It is the extension to original ePC-SAFT to account for high ion concentrations ( $\mathrm{x}_{\mathrm{IL}}>0.5$ ), which is important when calculating LLE of systems with ionic liquids. The new approach was compared to the original ePCSAFT, showing significant improvement in both prediction of LLE and correlation of LLE using one binary temperature-independent interaction parameter $k_{i j}$ between water and IL-ion. The new model ePC-SAFT does not require new pure-component parameters for the IL-ions. The dielectric constant is indeed a very important property that needs to be taken into account for modeling physically meaningful thermodynamic behavior of solutions containing high ion concentrations.

## Acknowledgment

This work has been supported by the German Science Foundation (DFG), Germany within the priority program SPP 1708 "Material Synthesis Near Room Temperature" (grant HE 7165/7-1).

## References

[1] P. Wasserscheid, T. Welton, Ionic Liquids in Synthesis, John Wiley \& Sons, 2008.
[2] T. Welton, Room-temperature ionic liquids. Solvents for synthesis and catalysis, Chem. Rev. 99 (1999) 2071-2084, https://doi.org/10.1021/cr980032t.
[3] K.R. Seddon, Ionic liquids for clean technology, J. Chem. Technol. Biotechnol. 68 (1997) 351-356, https://doi.org/10.1002/(SICI)1097-4660(199704)68: 4<351:AID-JCTB613>3.0.CO;2-4.
[4] Ionic liquids as green solvents, in: R.D. Rogers, K.R. Seddon (Eds.), ACS Symposium Series, American Chemical Society, Washington, DC, 2003.
[5] J.D. Holbrey, M.B. Turner, R.D. Rogers, Selection of ionic liquids for green chemical applications. In ionic liquids as green solvents, in: R.D. Rogers, K.R. Seddon (Eds.), ACS Symposium Series, American Chemical Society, Washington, DC, 2003, pp. 2-12.
[6] M.J. Earle, K.R. Seddon, Ionic liquids. Green solvents for the future, Pure Appl. Chem. 72 (2000), https://doi.org/10.1351/pac200072071391.
[7] I. Dinarès, C. Garcia de Miguel, A. Ibáñez, N. Mesquida, E. Alcalde, Imidazolium ionic liquids: a simple anion exchange protocol, Green Chem. 11 (2009) 1507, https://doi.org/10.1039/b915743n.
[8] X. Han, D.W. Armstrong, Ionic liquids in separations, Acc. Chem. Res. 40 (2007) 1079-1086, https://doi.org/10.1021/ar700044y.
[9] H. Zhao, S. Xia, P. Ma, Use of ionic liquids as 'green' solvents for extractions, J. Chem. Technol. Biotechnol. 80 (2005) 1089-1096, https://doi.org/10.1002/ jctb. 1333.
[10] A. Schmid, J.S. Dordick, B. Hauer, A. Klener, M. Wubbolts, B. Witholt, Industrial biocatalysis today and tomorrow, Nature 2001 (2001) 258-268.
[11] M. Voges, I.V. Prikhodko, S. Prill, M. Hübner, G. Sadowski, C. Held, Influence of pH value and ionic liquids on the solubility of l -alanine and l -glutamic acid in aqueous solutions at $30{ }^{\circ} \mathrm{C}$, J. Chem. Eng. Data 62 (2017) 52-61, https:// doi.org/10.1021/acs.jced.6b00367.
[12] U. Domańska, M. Królikowski, A. Pobudkowska, P. Bocheńska, Solubility of ionic liquids in water and octan-1-ol and octan-1-ol/water, or 2phenylethanol/water partition coefficients, J. Chem. Thermodyn. 55 (2012) 225-233, https://doi.org/10.1016/j.jct.2012.06.003.
[13] U. Domańska, M. Królikowski, M. Zawadzki, A. Wróblewska, Phase equilibrium investigation with ionic liquids and selectivity in separation of 2phenylethanol from water, J. Chem. Thermodyn. 102 (2016) 357-366, https://doi.org/10.1016/j.jct.2016.07.025.
[14] U. Domańska, M. Królikowski, D. Ramjugernath, T.M. Letcher, K. Tumba, Phase equilibria and modeling of pyridinium-based ionic liquid solutions, J. Phys. Chem. B 114 (2010) 15011-15017, https://doi.org/10.1021/jp105825c.
[15] U. Domańska, E.V. Lukoshko, M. Królikowski, Phase behaviour of ionic liquid 1-butyl-1-methylpyrrolidinium tris(pentafluoroethyl)trifluorophosphate with alcohols, water and aromatic hydrocarbons, Fluid Phase Equilib. 345 (2013) 18-22, https://doi.org/10.1016/j.fluid.2013.01.027.
[16] N. Deenadayalu, K.C. Ngcongo, T.M. Letcher, D. Ramjugernath, Liquid-Liquid equilibria for ternary mixtures (an ionic liquid + Benzene + Heptane or Hexadecane) at $\mathrm{T}=298.2 \mathrm{~K}$ and atmospheric pressure, J. Chem. Eng. Data 51 (2006) 988-991, https://doi.org/10.1021/je0504941.
[17] M.G. Freire, P.J. Carvalho, R.L. Gardas, I.M. Marrucho, L.M.N.B.F. Santos, J.A.P. Coutinho, Mutual solubilities of water and the C(n)mimTf(2)N hydrophobic ionic liquids, J. Phys. Chem. B 112 (2008) 1604-1610, https://doi.org/ 10.1021/jp7097203.
[18] M.G. Freire, C.M.S.S. Neves, P.J. Carvalho, R.L. Gardas, A.M. Fernandes, I.M. Marrucho, L.M.N.B.F. Santos, J.A.P. Coutinho, Mutual solubilities of water and hydrophobic ionic liquids, J. Phys. Chem. B 111 (2007) 13082-13089, https://doi.org/10.1021/jp076271e.
[19] M.G. Freire, C.M.S.S. Neves, K. Shimizu, C.E.S. Bernardes, I.M. Marrucho, J.A.P. Coutinho, J.N. Canongia Lopes, L.P.N. Rebelo, Mutual solubility of water and structural/positional isomers of N-alkylpyridinium-based ionic liquids, J. Phys. Chem. B 114 (2010) 15925-15934, https://doi.org/10.1021/ jp1093788.
[20] M.G. Freire, C.M.S.S. Neves, S.P.M. Ventura, M.J. Pratas, I.M. Marrucho, J. Oliveira, J.A.P. Coutinho, A.M. Fernandes, Solubility of non-aromatic ionic liquids in water and correlation using a QSPR approach, Fluid Phase Equilib. 294 (2010) 234-240, https://doi.org/10.1016/j.fluid.2009.12.035.
[21] A.R. Ferreira, M.G. Freire, J.C. Ribeiro, F.M. Lopes, J.G. Crespo, J.A.P. Coutinho, Overview of the liquid-liquid equilibria of ternary systems composed of ionic liquid and aromatic and aliphatic hydrocarbons, and their modeling by COSMO-RS, Ind. Eng. Chem. Res. 51 (2012) 3483-3507, https://doi.org/ 10.1021/ie2025322.
[22] C. Tsioptsias, I. Tsivintzelis, C. Panayiotou, Equation-of-state modeling of mixtures with ionic liquids, Phys. Chem. Chem. Phys. : Phys. Chem. Chem. Phys. 12 (2010) 4843-4851, https://doi.org/10.1039/c000208a.
[23] A. Nann, C. Held, G. Sadowski, Liquid-liquid equilibria of 1-butanol/water/IL systems, Ind. Eng. Chem. Res. 52 (2013) 18472-18481, https://doi.org/ 10.1021/ie403246e.
[24] L.F. Cameretti, G. Sadowski, J.M. Mollerup, Modeling of aqueous electrolyte solutions with perturbed-chain statistical associated fluid theory, Ind. Eng. Chem. Res. 44 (2005) 3355-3362, https://doi.org/10.1021/ie0488142.
[25] X. Ji, C. Held, Modeling the density of ionic liquids with ePC-SAFT, Fluid Phase Equilib. 410 (2016) 9-22, https://doi.org/10.1016/j.fluid.2015.11.014.
[26] X. Ji, C. Held, G. Sadowski, Modeling imidazolium-based ionic liquids with ePC-SAFT, Fluid Phase Equilib. 335 (2012) 64-73, https://doi.org/10.1016/ j.fluid.2012.05.029.
[27] X. Ji, C. Held, G. Sadowski, Modeling imidazolium-based ionic liquids with ePC-SAFT. Part II. Application to H2S and synthesis-gas components, Fluid Phase Equilib. 363 (2014) 59-65, https://doi.org/10.1016/j.fluid.2013.11.019.
[28] C.M.S.S. Neves, C. Held, S. Mohammad, M. Schleinitz, J.A.P. Coutinhoa, M.G. Freire, Effect of salts on the solubility of ionic liquids in water: experimental and electrolyte perturbed-chain statistical associating fluid theory, Phys. Chem. Chem. Phys. : Phys. Chem. Chem. Phys. 17 (2015) 32044-32052, https://doi.org/10.1039/c5cp06166k.
[29] G. Shen, C. Held, J.-P. Mikkola, X. Lu, X. Ji, Modeling the viscosity of ionic liquids with the electrolyte perturbed-chain statistical association fluid theory, Ind. Eng. Chem. Res. 53 (2014) 20258-20268, https://doi.org/10.1021/ ie503485h.
[30] G. Shen, C. Held, X. Lu, X. Ji, Modeling thermodynamic derivative properties of ionic liquids with ePC-SAFT, Fluid Phase Equilib. 405 (2015) 73-82, https:// doi.org/10.1016/j.fluid.2015.07.018.
[31] M. Voges, F. Schmidt, D. Wolff, G. Sadowski, C. Held, Thermodynamics of the alanine aminotransferase reaction, Fluid Phase Equilib. 422 (2016) 87-98, https://doi.org/10.1016/j.fluid.2016.01.023.
[32] K. Klauke, D.H. Zaitsau, M. Bülow, L. He, M. Klopotowski, T.-O. Knedel, J. Barthel, C. Held, S.P. Verevkin, C. Janiak, Thermodynamic properties of selenoether-functionalized ionic liquids and their use for the synthesis of zinc selenide nanoparticles, Dalton Trans. 47 (2018) 5083-5097, https://doi.org/ 10.1039/C8DT00233A. Cambridge, England : 2003.
[33] M. Born, Volumen und Hydratationswärme der Ionen, Z. Phys. 1 (1920) 45-48, https://doi.org/10.1007/BF01881023.
[34] W.G. Chapman, K.E. Gubbins, G. Jackson, M. Radosz, SAFT: equation-of-state
solution model for associating fluids, Fluid Phase Equilib. 52 (1989) 31-38, https://doi.org/10.1016/0378-3812(89)80308-5.
[35] M.S. Wertheim, Fluids with highly directional attractive forces. I. Statistical thermodynamics, J. Stat. Phys. 35 (1984) 19-34, https://doi.org/10.1007/ BF01017362.
[36] M.S. Wertheim, Fluids with highly directional attractive forces. II. Thermodynamic perturbation theory and integral equations, J. Stat. Phys. 35 (1984) 35-47, https://doi.org/10.1007/BF01017363.
[37] J. Gross, G. Sadowski, Perturbed-chain SAFT: an equation of state based on a perturbation theory for chain molecules, Ind. Eng. Chem. Res. 40 (2001) 1244-1260, https://doi.org/10.1021/ie0003887.
[38] C. Held, T. Reschke, S. Mohammad, A. Luza, G. Sadowski, ePC-SAFT revised, Chem. Eng. Res. Des. 92 (2014) 2884-2897, https://doi.org/10.1016/ j.cherd.2014.05.017.
[39] J. Gross, G. Sadowski, Application of the perturbed-chain SAFT equation of state to associating systems, Ind. Eng. Chem. Res. 41 (2002) 5510-5515, https://doi.org/10.1021/ie010954d.
[40] P. Debye, E. Hückel, Zur Theorie der Elektrolyte: I. Gefrierpunktserniedrigung und verwandte Erscheinungen, Z. Phys. 1923 (1923) 185-206.
[41] H. Weingärtner, Understanding ionic liquids at the molecular level: facts, problems, and controversies, Angew. Chem. 47 (2008) 654-670, https:// doi.org/10.1002/anie.200604951.
[42] M. Koeberg, C.-C. Wu, D. Kim, M. Bonn, THz dielectric relaxation of ionic liquid: water mixtures, Chem. Phys. Lett. 439 (2007) 60-64, https://doi.org/ 10.1016/j.cplett.2007.03.075.
[43] C. Held, L.F. Cameretti, G. Sadowski, Modeling aqueous electrolyte solutions, Fluid Phase Equilib. 270 (2008) 87-96, https://doi.org/10.1016/ j.fluid.2008.06.010.
[44] S. Wang, J. Jacquemin, P. Husson, C. Hardacre, M.F. Costa Gomes, Liquid-liquid miscibility and volumetric properties of aqueous solutions of ionic liquids as a function of temperature, J. Chem. Thermodyn. 41 (2009) 1206-1214, https:// doi.org/10.1016/j.jct.2009.05.009.
[45] J.M. Crosthwaite, S.N.V.K. Aki, E.J. Maginn, J.F. Brennecke, Liquid phase behavior of imidazolium-based ionic liquids with alcohols, J. Phys. Chem. B 108 (2004) 5113-5119, https://doi.org/10.1021/jp037774x.
[46] C. Wertz, J.K. Lehmann, A. Heintz, Liquid-liquid phase equilibria of the binary ionic liquid systems [C n MIM][NTf 2 ] +n -butanol, +n -pentanol, +H 2 O using UV spectroscopic and densimetric analytical methods, J. Chem. Eng. Data 58 (2013) 2375-2380, https://doi.org/10.1021/je300672q.
[47] N.V. Shvedene, S.V. Borovskaya, V.V. Sviridov, E.R. Ismailova, I.V. Pletnev, Measuring the solubilities of ionic liquids in water using ion-selective electrodes, Anal. Bioanal. Chem. 381 (2005) 427-430, https://doi.org/10.1007/ s00216-004-3001-7.
[48] D. Santos, M. Góes, E. Franceschi, A. Santos, C. Dariva, M. Fortuny, S. Mattedi, PHASE EQUILIBRIA FOR BINARY SYSTEMS CONTAINING IONIC LIQUID WITH WATER OR HYDROCARBONS, Braz. J. Chem. Eng. 32 (2015) 967-974, https:// doi.org/10.1590/0104-6632.20150324s00003609.
[49] K. Rehák, P. Morávek, M. Strejc, Determination of mutual solubilities of ionic liquids and water, Fluid Phase Equilib. 316 (2012) 17-25, https://doi.org/ 10.1016/j.fluid.2011.12.008.
[50] A. Chapeaux, L.D. Simoni, M.A. Stadtherr, J.F. Brennecke, Liquid phase behavior of ionic liquids with water and 1-Octanol and modeling of 1-Octanol/water partition coefficients, J. Chem. Eng. Data 52 (2007) 2462-2467, https:// doi.org/10.1021/je7003935.
[51] C.M.S.S. Neves, M.L.S. Batista, A.F.M. Cláudio, L.M.N.B.F. Santos, I.M. Marrucho, M.G. Freire, J.A.P. Coutinho, Thermophysical properties and water saturation of [PF 6 ]-Based ionic liquids, J. Chem. Eng. Data 55 (2010) 5065-5073, https:// doi.org/10.1021/je100638g.
[52] J.L. Anthony, E.J. Maginn, J.F. Brennecke, Solution thermodynamics of imidazolium-based ionic liquids and water, J. Phys. Chem. B 105 (2001) 10942-10949, https://doi.org/10.1021/jp0112368.
[53] Y. Li, L.-S. Wang, S.-F. Cai, Mutual solubility of alkyl imidazolium Hexafluorophosphate ionic liquids and water, J. Chem. Eng. Data 55 (2010) 5289-5293, https://doi.org/10.1021/je1003059.
[54] D.S.H. Wong, J.P. Chen, J.M. Chang, C.H. Chou, Phase equilibria of water and ionic liquids [emim][PF6] and [bmim][PF6], Fluid Phase Equilib. 194-197 (2002) 1089-1095, https://doi.org/10.1016/S0378-3812(01)00790-7.
[55] F.M. Maia, O. Rodríguez, E.A. Macedo, Free energy of transfer of a Methylene group in biphasic systems of water and ionic liquids [C 3 mpip][NTf 2 ], [C 3 mpyrr][NTf 2 ], and [C 4 mpyrr][NTf 2 ], Ind. Eng. Chem. Res. 51 (2012) 8061-8068, https://doi.org/10.1021/ie300227f.
[56] E.J. González, E.A. Macedo, Influence of the number, position and length of the alkyl-substituents on the solubility of water in pyridinium-based ionic liquids, Fluid Phase Equilib. 383 (2014) 72-77, https://doi.org/10.1016/ j.fluid.2014.09.028.
[57] F.M. Maia, O. Rodríguez, E.A. Macedo, LLE for (water+ionic liquid) binary systems using [Cxmim][BF4] $(x=6,8)$ ionic liquids, Fluid Phase Equilib. 296 (2010) 184-191, https://doi.org/10.1016/j.fluid.2010.05.003.
[58] M. Wagner, O. Stanga, W. Schröer, Corresponding states analysis of the critical points in binary solutions of room temperature ionic liquids, Phys. Chem. Chem. Phys. 5 (2003) 3943-3950, https://doi.org/10.1039/ B305959F.
[59] F.M. Maia, O. Rodríguez, E.A. Macedo, Relative hydrophobicity of equilibrium phases in biphasic systems (ionic liquid+water), J. Chem. Thermodyn. 48 (2012) 221-228, https://doi.org/10.1016/j.jct.2011.12.025.


[^0]:    * Corresponding author.

    E-mail address: christoph.held@tu-dortmund.de (C. Held).

