# Predicting Thermodynamic Properties of lons in Single Solvents and in Mixed Solvents Using a Modified Born Term within the ePC-SAFT Framework 

Paul Figiel, Gangqiang Yu,* and Christoph Held*

Cite This: Ind. Eng. Chem. Res. 2025, 64, 9406-9418
Read Online
Downloaded via BRIGHAM YOUNG UNIV on May 15, 2025 at 17:06:32 (UTC). See https://pubs.acs.org/sharingguidelines for options on how to legitimately share published articles.


#### Abstract

In this work, the EoS ePC-SAFT was further developed by modifying the Born term and introducing a new mixing rule for the dielectric constant of electrolyte solutions. The modified Born term considers the solvent solvation-shell volume of ions and accounts for the radial dependence of the dielectric constant near a central ion. This required one additional parameter per ion and one per solvent, while other pure-component parameters were inherited from the literature. The performance ![](https://cdn.mathpix.com/cropped/d6b90f00-389e-416f-b97c-c8651d44a864-01.jpg?height=313&width=846&top_left_y=882&top_left_x=1090) of ePC-SAFT was validated against experimental literature data: Gibbs energies of solvation and of transfer and MIACs in pure water, methanol, ethanol, and in their mixtures. The new modeling approach improves the accuracy of Gibbs energies at infinite dilution while retaining the ability to model properties of concentrated electrolyte systems. Additionally, ePC-SAFT parameters for $\mathrm{V}^{3+}$ and $\mathrm{VO}^{2+}$ ions were fitted to density and osmotic-coefficient data of $\mathrm{VCl}_{3}$ and $\mathrm{VOSO}_{4}$ in water using the new approach.


## - INTRODUCTION

An increasing interest in environmentally benign technologies has been observed in recent years due to the threat of climate change. This was accompanied by a shift in the energy sector from fossil to renewable energy sources. ${ }^{1,2}$ However, renewable energy sources, such as solar or wind energy, are not able to provide electricity on demand. Even though great progress in energy-storage systems was achieved in the last years, their development and improvement are still ongoing. ${ }^{3,4}$ The most prominent representative of energy-storage systems is the battery, which usually consists of two electrodes and a connecting electrolyte. The composition of these electrolytes ranges from aqueous electrolyte solutions in the case of vanadium-flow batteries to organic electrolyte solutions in Li ion batteries and even solid-state electrolyte systems. ${ }^{5-7}$ A lot of research interest lies in developing better electrolyte solutions for these systems to reduce viscosity and increase power density, stability, safety, and voltage range of the battery. ${ }^{7,8}$ This includes both the development of better solvents or solvent mixtures and the use of different ions in the solution. One major goal is the replacement of lithium ions in Li-ion batteries with much cheaper and safer materials such as sodium ions. ${ }^{3,9}$ It is of great interest to understand how different ions behave in different solvents for these developments. Electrolyte solutions are also used in many other industrial fields. They are present in (reverse) osmosis processes, carbon-capture technologies, and electrolysis processes. ${ }^{10-12}$ Recent research focuses on replacing chemical processes (fossil-based reactants) with electrochemical
routes. ${ }^{13}$ One example is the production of ammonia using an electrolysis process with a LiOH solution. While originally only LiOH aqueous solution were tested, organic solvents and salt mixtures have also been used in more recent research in this process. ${ }^{13-15}$ Again, knowledge of the interactions between solvents and ions aids the development in these areas. A lot of thermodynamic data (e.g., density, infinitedilution solvation properties, activity coefficients) on saltsolvent systems can be found in the literature for common electrolytes (e.g., alkali halides) in common solvents (e.g., water, ethanol). However, most literature sources include data on only single-solvent single-electrolyte mixtures. However, data are scarce for an electrolyte in a solvent mixture. Electrolyte thermodynamic models can aid in these cases to describe (or even predict) the properties of more complex mixtures in order to reduce the high experimental effort to characterize multisolvent multielectrolyte mixtures. Many different models exist in the literature for electrolyte solutions. ${ }^{16-21}$ In general, such models were developed to describe aqueous electrolyte solutions, but such models have been subsequently improved to also describe nonaqueous electrolyte solutions. In this context, the implementation of the

[^0]Born term was intensely discussed in the literature, and many models included this term, which resulted in better modeling results in nonaqueous electrolyte solutions. ${ }^{16,22-25}$ The Born term calculates the transfer energy of an ion from a vacuum into any desired solvent, and this energy contribution is usually comparatively high. One of the models that include the Born term is electrolyte PC-SAFT (ePC-SAFT advanced ${ }^{16}$ ), which combines PC-SAFT with Coulomb forces as well as iondipolar interactions. It has already been used to model mean ionic activity coefficients (MIACs) and osmotic coefficients in aqueous and organic solvents. ${ }^{16,26}$ Further, ePC-SAFT was used to calculate salt solubilities in these solvents and to predict $\mathrm{CO}_{2}$ solubilities in aqueous and organic electrolyte solutions ${ }^{27,28}$ and even to predict the salt influence on esterification reaction equilibria and kinetics. ${ }^{29,30}$ All of these works focused on systems with moderate to high ion concentrations. However, the properties of electrolyte solutions at infinite dilution were not investigated in these works. Properties at infinite dilution (like the Gibbs energy of solvation or the Gibbs energy of transfer) are very important, e.g., for the calculation of battery cell voltages in different solvents. The aim of this work is to further improve ePC-SAFT advanced to better describe these infinite-dilution properties of electrolyte solutions while retaining the ability of ePC-SAFT to model properties of concentrated electrolyte solutions.

## - EPC-SAFT

ePC-SAFT is an extension of the PC-SAFT EoS in order to access the modeling of charged species. In general, PC-SAFT calculates the residual Helmholtz energy $a^{\text {res }}$ by evaluating different contributions. These contributions are the hard-chain contribution $a^{\text {hc }}$, the dispersion contribution $a^{\text {disp }}$, and the association contribution $a^{\text {assoc }}$. In ePC-SAFT, two additional contributions for charged species are included: the DebyeHückel contribution $a^{\mathrm{DH}}$ for ion-ion interactions and the Born contribution $a^{\text {Born }}$ for ion-dipolar interactions. The residual Helmholtz energy is obtained by summing these contributions (eq 1)

$$
\begin{equation*}
a^{\mathrm{res}}=a^{\mathrm{hc}}+a^{\mathrm{disp}}+a^{\mathrm{assoc}}+a^{\mathrm{DH}}+a^{\mathrm{Born}} \tag{1}
\end{equation*}
$$

Original ePC-SAFT requires two ion-specific parameters, while three pure-component parameters for nonassociating components and five pure-component parameters for associating components are required. These are the segment number $m_{i}^{\text {seg }}$, the segment diameter $\sigma_{\mathrm{i}}$, and the dispersion energy parameter $u_{i} / k_{\mathrm{B}}$ for nonassociating components, where $k_{\mathrm{B}}$ represents the Boltzmann constant. Associating components also require the association-energy parameter $\varepsilon^{\mathrm{AiBi}}$ and the association-volume parameter $\kappa^{\mathrm{AiBi}}$. Additionally, the binary interaction parameter $k_{i j}$ between every component pair might be adjusted to alter the dispersion energy in the mixture. The corresponding combining rule is given in eq 2

$$
\begin{equation*}
u_{i j}=\sqrt{u_{i} u_{j}} \cdot\left(1-k_{i j}\right) \tag{2}
\end{equation*}
$$

Note that the binary interaction parameter $k_{\mathrm{ij}}$ might depend on the temperature.

Both the expressions for $a^{\mathrm{DH}}$ and $a^{\text {Born }}$ require a diameter that characterizes the ionic species. In all former ePC-SAFT modeling approaches, the temperature-corrected segment diameter of the ion $d_{\text {ion }}$ was used in both expressions, $a^{\mathrm{DH}}$ and $a^{\text {Born }}$. eq 3 shows the general formula for temperaturecorrected segment diameter $d_{i}(T)$

$$
\begin{equation*}
d_{i}(T)=\sigma_{i} \cdot\left(1-0.12 \cdot \mathrm{e}^{-3 \cdot u_{i i} / k_{\mathrm{B}} T}\right) \tag{3}
\end{equation*}
$$

This equation can be simplified for ions since the dispersion energy between cation-cation pairs and anion-anion pairs is set to zero in ePC-SAFT ( $u_{i i}=0$ ), which leads to eq 4

$$
\begin{equation*}
d_{\mathrm{ion}}=\sigma_{\mathrm{ion}} \cdot 0.88 \tag{4}
\end{equation*}
$$

In this work, no changes were made in the expressions for $a^{\mathrm{DH}}$, and $d_{\text {ion }}$ is still used in $a^{\text {DH }}$. However, a new diameter parameter for the Born term is introduced in this work. This parameter was not tested at different temperatures and did not include a temperature dependency.

## - BORN TERM

The Born term is used in many thermodynamic models to calculate the transfer work of an ion from a vacuum into a pure solvent. ${ }^{16,23,31}$ In the Born term, the solvent is treated as a continuum that only interacts with the ion via its dielectric constant $\varepsilon_{r}{ }^{32}$ The energy connected to this process of a single ion can be summed for all ions in the system, yielding the reduced Helmholtz energy contribution from Born term $a^{\text {Born }}$, given in eq 5

$$
\begin{equation*}
a^{\mathrm{Born}}=-\frac{e^{2}}{4 \pi \varepsilon_{0} k_{\mathrm{B}} T}\left(1-\frac{1}{\varepsilon_{r}}\right) \sum_{i} \frac{x_{i} z_{i}^{2}}{d_{i}^{\mathrm{Born}}} \tag{5}
\end{equation*}
$$

Here, $d_{i}^{\text {Born }}$ is the Born diameter of ion $i$. This quantity is not well-defined, and many different possibilities exist for its determination. In all former ePC-SAFT modeling approaches, $d_{i}^{\text {Born }}$ was set to $\sigma_{\mathrm{i}}$ in order to reduce the number of the model parameters. However, this led to deviations between the modeled and experimental Gibbs energies of solvation at infinite dilution of the considered ion, $\Delta^{\text {solv }} G_{i}^{\infty} .^{16}$ Other authors used Pauling diameters or cavity diameters to estimate $d_{i}^{\text {Born }}$, which also did not yield optimal results for $\Delta^{\text {solv }} G_{i}^{\infty}$ when compared to literature data. ${ }^{33-35}$ Some authors suggested fitting $d_{i}^{\text {Born }}$ to literature $\Delta^{\text {solv }} G_{i}^{\infty}$ data. ${ }^{34,36}$ Another approach is to not only look at the diameter of the ion but also include the solvent solvation shell in the Born term. ${ }^{34,35,37-40}$ This is done by introducing an additional quantity into Born term $\Delta d_{i}$ that accounts for this effect. This quantity might depend on the kind of ion, the solvent, and the composition of the system. Approaches that include this quantity will be called "solvation shell models" (SSMs) in this work. The formula for the reduced Helmholtz energy contribution from the Born term with an SSM approach is given in eq 6

$$
\begin{equation*}
a_{\mathrm{SSM}}^{\mathrm{Born}}=-\frac{e^{2}}{4 \pi \varepsilon_{0} k_{\mathrm{B}} T}\left(1-\frac{1}{\varepsilon_{r}}\right) \sum_{\mathrm{i}} \frac{x_{i} z_{i}^{2}}{d_{i}^{\mathrm{Born}}+\Delta d_{i}} \tag{6}
\end{equation*}
$$

An SSM assumes that the volume between the surface of the ion and the first solvation shell does not contribute to the Helmholtz energy. This assumption does not hold when looking at the changes in the dielectric constant around an ion. The dielectric constant is much lower close to an ion in a solvent environment compared to the bulk phase since the solvent molecules do not affect the electric fields at those distances. ${ }^{41-43}$ This effect is called "dielectric saturation" (DS) and can be included in the Born term by introducing a function for the dielectric constant that depends on the distance from the ion. While it is possible to obtain a radiusdependent function for the dielectric constant, the effort is very high, especially if different solvents, solvent mixtures, and
concentrated electrolyte solutions are considered. An easier equation can be obtained by assuming a constant dielectric constant between the ion surface and the solvation shell that increases instantaneously to the bulk dielectric constant at distances greater than the radius of the solvation shell. ${ }^{39,40}$ In this work, the dielectric constant between the ion and the first solvation shell is assumed to be equal to the dielectric constant of the ion. The value of the relative dielectric constant of all ions is assumed to be $\varepsilon_{r}=8$ in the ePC-SAFT calculations. This value was determined by Bülow et al. as the average of the relative dielectric constant of different salts taken from the literature. ${ }^{16,44}$ The simplified dependence of the relative dielectric constant on the radius is schematically shown in Figure 1. Both the solvation shell and the dielectric saturation

![](https://cdn.mathpix.com/cropped/d6b90f00-389e-416f-b97c-c8651d44a864-03.jpg?height=464&width=644&top_left_y=775&top_left_x=273)
Figure 1. Relative dielectric constant $\varepsilon_{r}$ as a function of the radius $r$ as assumed in combined solvation shell and dielectric saturation models. $r_{i}$ is the radius of the ion, and $\Delta r_{i}$ is the distance between the ion surface and the first solvation shell.

effect are included in this combined approach (SSM + DS), which results in eq 7 for the reduced Helmholtz energy contribution of the Born term.

$$
\begin{align*}
a_{\mathrm{SSM}+\mathrm{DS}}^{\mathrm{Born}}= & -\frac{e^{2}}{4 \pi \varepsilon_{0} k_{\mathrm{B}} T}\left(1-\frac{1}{\varepsilon_{r, \text { bulk }}}\right)\left(\sum_{i} \frac{x_{i} z_{i}^{2}}{d_{i}^{\mathrm{Born}}+\Delta d_{i}}\right) \\
& +\left(1-\frac{1}{\varepsilon_{r, \text { bulk }}}\right)\left(\sum_{i} x_{i} z_{i}^{2}\left(\frac{1}{d_{i}^{\text {Born }}}-\frac{1}{d_{i}^{\text {Born }}+\Delta d_{i}}\right)\right) \tag{7}
\end{align*}
$$

For SSMs, the Born diameter and the solvation shell thickness have the same effect on the Helmholtz energy (eq 6), and it could be argued that the quantities cannot be considered independently. In contrast, for SSM+DS, a separate consideration for the Born diameter $d_{i}^{\text {Born }}$ and the combined diameter $d_{i}^{\text {Born }}+\Delta d_{i}$ is necessary since the consideration of the dielectric saturation introduces a contribution that only depends on the Born diameter (eq 7).

It is also possible in these kinds of models to explicitly consider more than one solvation shell layer. ${ }^{45}$ However, it was found that the performance does not significantly improve compared to the models that use only one solvation shell layer. ${ }^{46}$ This consideration might not even be necessary if parameters are fitted to literature data since these data include all solvation shell layers that might exist and therefore can only yield parameters that implicitly include these effects. A comparison of the performance of the three approaches (original, SSM, SSM+DS) is given in Chapter S1 and S2 in
the Supporting Information. The SSM + DS approach was chosen for all further calculations due to its superior performance compared with the other approaches.

In the literature, numerous approaches exist to calculate the solvation shell diameter. However, most of them only consider ions at infinite dilution, and such approaches yield poor results for concentrated solutions or are only applicable to water (for details, see Chapters S1 and S2 in the Supporting Information). Thus, a new parameter $f_{\text {mix }}$ was included in ePC-SAFT that depends on the solvent as well as on the ion concentration but not on the kind of ion. The solvation shell thickness was calculated using eq 8

$$
\begin{equation*}
\Delta d_{i}=\frac{\left(f_{\mathrm{mix}}-1\right)}{\left|z_{i}\right|} \cdot d_{i}^{\mathrm{Born}} \tag{8}
\end{equation*}
$$

The parameter $f_{\text {solv, mix }}$ was calculated via a mixing rule (eq 9)

$$
\begin{equation*}
f_{\operatorname{mix}}=\frac{\sum_{k} x_{k} f_{k}}{\sum_{k} x_{k}} \tag{9}
\end{equation*}
$$

Here, $k$ is the number of ions and solvents in the system, and $f_{k}$ is the "pure-component" adjustable parameter for each of these components. The $f_{k}$ parameter was set to $f_{k}=1$ for all ions since no solvation shell exists in a system consisting of pure ions. For now, it was chosen to include only ions and solvents in this mixing rule. This leads to a simplification of eq 9 since no other components are considered in this work (eq 10)

$$
\begin{equation*}
f_{\operatorname{mix}}=\sum_{k} x_{k} f_{k} \tag{10}
\end{equation*}
$$

The parameter for each solvent was adjusted to MIACs of NaBr in the respective solvent considered in this work (water, methanol, and ethanol). It should be noted that the $f_{k}$ value is constant for each solvent. In eq $8, \Delta d_{i}$ is normalized to the absolute charge of the respective ion $i$. This caused significantly better modeling results for multivalent ions.

## - DIELECTRIC CONSTANT

An important quantity for modeling electrolyte solutions with ePC-SAFT is the dielectric constant. The relative dielectric constants of all components modeled in this work are given in Table 1 at 298.15 K and 1 bar . All ions were characterized with

Table 1. Relative Dielectric Constants of All Solvents and Ions Used in This Work Are Shown at 298.15 K and 1 bar
| component | dielectric constant $(-)$ | refs |
| :--- | :---: | :--- |
| water | 78.09 | 47 |
| methanol | 33.05 | 48 |
| ethanol | 24.88 | 49 |
| salts (Ions) | 8 | $16^{a}$ |


${ }^{a}$ All salts were modeled with a similar dielectric constant that is a mean of available experimental data. ${ }^{44}$
the same value of 8 for the relative dielectric constant, which is the average of a number of different salts. This approach was also used in earlier versions of ePC-SAFT. ${ }^{16}$

Earlier modeling approaches using the ePC-SAFT framework assumed a linear dependency of the relative dielectric constant on the ion mole fraction. ${ }^{16,50,51}$ This assumption was validated against mixtures of ionic liquids with water. However, a linear dependency of the relative dielectric constant on the ion mole fraction cannot be observed in electrolyte solutions
that include salts, e.g., alkali halides. ${ }^{52-60}$ The relative dielectric constant decreases strongly at low ion mole fractions and gets less steep at higher ion mole fractions in solvents (in this work: water, methanol, and ethanol). A new mixing rule based on the ideas of Zuber et al. is introduced in this work to capture this behavior. ${ }^{52}$ For that purpose, relative dielectric constants in water and methanol were correlated over the ion mole fraction. This new mixing rule was also used to calculate the relative dielectric constants of salt solutions in ethanol to validate it in a different solvent that was not used in the correlation. The new mixing rule is given in eq 11

$$
\begin{equation*}
\varepsilon_{r, \text { mix }}=\frac{\varepsilon_{r, \text { solvent }, \text { mix }}^{\text {salt }- \text { free }}}{1+7.01 \cdot x_{\text {ion }}} \tag{11}
\end{equation*}
$$

Here, $\varepsilon_{r, \text { solvent, mix }}^{\text {salt-fre }}$ is the relative dielectric constant of the solvent mixture without ions, and $x_{\text {ion }}$ is the mole fraction of all ions in the system. In Figure 2, the results of this mixing rule are

![](https://cdn.mathpix.com/cropped/d6b90f00-389e-416f-b97c-c8651d44a864-04.jpg?height=568&width=644&top_left_y=916&top_left_x=271)
Figure 2. Relative dielectric constants of different aqueous salt solutions at 298.15 K and 1 bar. Solid line: linear mixing rule; dashed line: new mixing rule introduced in this work by eq 11 . Gray circles: Data of NaCl in water by Nörtemann et al. ${ }^{53}$ Green squares: Data of NaBr in water by Haggis et al. ${ }^{55}$ Blue triangles: Data of NaCl in water by Buchner et al. ${ }^{56}$ Orange upside-down triangles: Data of LiCl in water by Wei et al. ${ }^{54}$.

compared to experimental relative dielectric constants of salt solutions in pure water as well as to the linear mole fractionbased mixing rule.

The new mixing allows mimicking the experimental data of numerous salts in water sufficiently well, while the previously used linear mixing rule overestimates the experimental data in the considered mole fraction range. Figure 3 shows the performance of the new mixing for electrolyte solutions in water, methanol, and ethanol. The dielectric constants in all three solvents can be described without using solvent-specific parameters or ion-specific parameters. However, it should be mentioned that this new mixing rule was not tested in any other solvent, and the results might be poor for solvents with very low dielectric constants. ${ }^{61}$ Additionally, this new mixing rule was correlated using experimental data points up to about $x_{\text {ion }}=0.38$ and might deviate from experimental data points at higher ion concentrations. This is acceptable since the solubility of most salts is much lower than this limit. It is expected that the mixing rule is valid for temperatures other than 298.15 K only if the temperature influence on the pure solvents' dielectric constants is considered quantitatively

![](https://cdn.mathpix.com/cropped/d6b90f00-389e-416f-b97c-c8651d44a864-04.jpg?height=503&width=636&top_left_y=199&top_left_x=1194)
Figure 3. Relative dielectric constants of different salt solutions in different solvents over the relative dielectric constant of the pure solvent at 298.15 K and 1 bar. Solid line: Mixing rule introduced in this work by eq 11; gray circles: Data of different salts in water taken from refs 41,53-56; blue triangles: Data of different salts in methanol taken from refs 57,58; green squares: Data of different salts in ethanol taken from refs 59,60.

correct. However, this hypothesis was not tested in this work since all systems were investigated at 298.15 K .

Similar to earlier ePC-SAFT works, the dielectric constant of a solvent mixture $\varepsilon_{r, \text { solvent,mix }}^{\text {salt-fre }}$ was obtained by applying a weight fraction-based mixing rule (eq 12) to the pure-component dielectric constants ${ }^{50}$

$$
\begin{equation*}
\varepsilon_{r, \text { solvent }, \text { mix }}^{\text {salt }- \text { free }}=\sum_{\text {solvent }} w_{\text {solvent }}^{\text {salt }- \text { free }} \cdot \varepsilon_{r, \text { solvent }} \tag{12}
\end{equation*}
$$

Here, $w_{\text {solvent }}^{\text {salt-free }}$ is the weight fraction of each solvent in the system without salt. This equation (together with eq 11) allows the estimation of the relative dielectric constant in all multicomponent systems, given that the pure component dielectric constants are known (which in this work were taken from Table 5).

## - CALCULATION OF THERMODYNAMIC PROPERTIES

The fugacity coefficient of each component in a system can be calculated if the residual Helmholtz energy and its derivations are known by using eq 13

$$
\begin{align*}
\ln \left(\varphi_{i}\right)= & (Z-1)-\ln (Z)+a^{\mathrm{res}}+\left(\frac{\partial a^{\mathrm{res}}}{\partial x_{i}}\right) \\
& -\sum_{j} x_{j}\left(\frac{\partial a^{\mathrm{res}}}{\partial x_{j}}\right) \tag{13}
\end{align*}
$$

The fugacity coefficient can be used to calculate different quantities of the system.

Solvation describes the process of transferring a species from a pure, ideal gaseous phase (g) into a liquid phase (l). The corresponding Gibbs energy of solvation of a species $\Delta^{\text {solv }} G_{\mathrm{i}}\left(T, p, x_{j}\right)$ can be calculated as the difference of its chemical potential $\mu_{i}$ in both phases as given in eq 14 .

$$
\begin{equation*}
\Delta^{\text {solv }} G_{i}\left(T, p, x_{j}\right)=\mu_{i}^{\mathrm{l}}\left(T^{\mathrm{l}}, p^{\mathrm{l}}, x_{j}^{\mathrm{l}}\right)-\mu_{0 i}^{\mathrm{id}, \mathrm{~g}}\left(T^{\mathrm{g}}, p^{\mathrm{g}}\right) \tag{14}
\end{equation*}
$$

The temperature is usually the same for both phases, while the pressure is often different for both phases depending on the chosen standard process. The role of the different standard processes is explained later in this chapter. The chemical
potential of a species $i$ can be expressed in the following way (eq 15):

$$
\begin{equation*}
\mu_{i}\left(T, p, x_{j}\right)=\mu_{0 i}^{\mathrm{id}, \mathrm{~g}}(T, p)+R T \ln \left(\frac{f_{i}}{p}\right) \tag{15}
\end{equation*}
$$

Here, $f_{i}$ is the fugacity of component $i$ in the mixture, and $\mu_{0 i}^{\mathrm{id}, g}(T, p)$ is the chemical potential of a pure, ideal gas at pressure $p$ and temperature $T$. The fugacity can be expressed via fugacity coefficient $\varphi_{i}$ and mole fraction $x_{i}$ as shown in eq 16

$$
\begin{equation*}
f_{i}=\varphi_{i} x_{i} p \tag{16}
\end{equation*}
$$

In theory, the pressure of the gaseous phase and the concentration in the liquid phase can be arbitrary. However, different standard treatments have been established that yield different Gibbs energies. The three most common standard processes found in the literature are the $\rho$-treatment, the $x$ treatment, and the $m$-treatment. The $\rho$-treatment corresponds to a transfer from an ideal gas at a given concentration (usually $1 \mathrm{~mol} / \mathrm{L}$ ) to a hypothetical solution with the same solute concentration at 1 bar pressure and the properties of an infinitely diluted system. ${ }^{62}$ Note that the pressure in both phases is not equal in this treatment. The $x$-treatment describes the transfer "from an ideal gas at 1 atm ( 1 bar in this work) pressure to a hypothetical dilute-ideal solution (at the same pressure) in which the mole fraction (of the solute) $s$ is unity" (see ${ }^{62}$ page 16). The $m$-treatment corresponds to a transfer from an ideal gaseous phase at 1 bar pressure to a hypothetical solution at the same pressure with a solute concentration of $1 \mathrm{~mol} / \mathrm{kg}$ and the properties of an infinitely diluted system. ${ }^{62}$ Note that this concentration is given in the molality scale and refers to 1 mol solute per 1 kg pure solvent (or solvent mixture). All three treatments assume the same temperature in the gaseous and liquid phases. In this work, the Gibbs energies of solvation are given for the $x$-treatment, and literature data was converted to this treatment when necessary. The Gibbs energy of solvation for the $x$-process was obtained by combining eqs $14-16$ and is given in eq 17

$$
\begin{align*}
\Delta^{\mathrm{solv}, x} G_{i}^{\infty}\left(T, p, x_{j \neq i}, x_{i} \rightarrow\right. & 0)=\mu_{0 i}^{\mathrm{id}, \mathrm{~g}}(T, p) \\
& +R T \ln \left(\frac{\varphi_{i}^{\infty} x_{i} p}{p}\right)-\mu_{0 i}^{\mathrm{id}, \mathrm{~g}}(T, p) \\
& =R T \ln \\
& \left(\varphi_{i}^{\infty}\left(T, p, x_{j \neq i}, x_{i} \rightarrow 0\right)\right) \tag{17}
\end{align*}
$$

Here, $\varphi_{i}^{\infty}$ is the fugacity coefficient of the ion at infinite dilution, which corresponds to the properties at infinite dilution as mentioned in the definition of the $x$-treatment. Therefore, the superscript $\infty$ is added to the Gibbs energies in this work. To convert values of the Gibbs energy of solvation from the $\rho$-treatment or the $m$-treatment, eq 18 can be used.

$$
\begin{align*}
\Delta^{\text {solv }, x} G_{i}^{\infty}\left(T, p, x_{j \neq i}, x_{i} \rightarrow\right. & 0) \\
& =\Delta^{\text {solv }, \rho} G_{i}^{\infty}\left(T, p, x_{j \neq i}, x_{i} \rightarrow 0\right) \\
& +R T \ln \left(\frac{R T \rho_{\text {solvent }}}{p^{0} M_{\text {solvent }}}\right) \\
& =\Delta^{\text {solv }, m} G_{i}^{\infty}\left(T, p, x_{j \neq i}, x_{i} \rightarrow 0\right) \\
& +R T \ln \left(\frac{1}{\tilde{m}^{0} M_{\text {solvent }}}\right) \tag{18}
\end{align*}
$$

Here, $p^{0}$ is the chosen standard pressure for the x -process (usually 1 bar), $\tilde{m}^{0}$ is the standard molality ( $1 \mathrm{~mol} / \mathrm{kg}$ ), $\rho_{\text {solvent }}$ the density of the solute-free solvent (or solvent mixture), and $M_{\text {solvent }}$ is the molar mass of the solvent (or solvent mixture). A detailed deviation of the conversion is found in the Supporting Information in Chapter S3. The superscript $x$ for the $x$-process will be omitted in the following text since only values for the x process are presented.

The Gibbs energy of transfer at infinite dilution from a solvent S 1 to another solvent or solvent mixture S 2 $\Delta^{\mathrm{tr}} G_{i}^{\infty, S 1} \rightarrow S 2$ was calculated as the difference of the Gibbs energies of solvation in the two solvents (eq 19)

$$
\begin{align*}
\Delta^{\operatorname{tr}} G_{i}^{\infty, S 1 \rightarrow S 2}\left(T, p, x_{j \neq i}, x_{i}\right. & \rightarrow 0) \\
& =\Delta^{\text {solv }} G_{i}^{\infty, S 2}\left(T, p, x_{j \neq i}, x_{i} \rightarrow 0\right) \\
& -\Delta^{\text {solv }} G_{i}^{\infty, S 1}\left(T, p, x_{j \neq i}, x_{i} \rightarrow 0\right) \\
& =R T \ln \left(\frac{\varphi_{i}^{\infty, S 2}}{\varphi_{i}^{\infty, S 1}}\right) \tag{19}
\end{align*}
$$

The activity coefficients of ions were calculated by using eq 20

$$
\begin{equation*}
\gamma_{i}^{*}=\frac{\varphi_{i}^{\mathrm{S}}}{\varphi_{i}^{\infty, \mathrm{S}}} \tag{20}
\end{equation*}
$$

The "*" indicates that the reference state for the activity coefficient is the infinitely diluted solution, which is most commonly used for ions. The mean ionic activity coefficient (MIAC) $\gamma_{ \pm}^{*}$ is usually given in the literature. It can be calculated from the single ion activity coefficients using eq 21

$$
\begin{equation*}
\left.\gamma_{ \pm}^{*}=\left(\gamma_{\mathrm{cation}}^{*} \cdot \gamma_{\mathrm{anion}}^{*}\right)^{*, \nu_{\mathrm{an}}}\right)^{1 / \nu_{\mathrm{cat}}+\nu_{\mathrm{an}}} \tag{21}
\end{equation*}
$$

Here, $\nu_{\mathrm{cat}}$ and $\nu_{\mathrm{an}}$ are the stoichiometric coefficients of the ions. MIACs are usually given in the molality scale $\gamma_{ \pm}^{*, m}$ in the literature while ePC-SAFT calculates mole fraction-based MIACs $\gamma_{ \pm}^{*, x}$. They were converted using eq 22

$$
\begin{equation*}
\gamma_{ \pm}^{*, m}=\frac{\gamma_{ \pm}^{*, x}}{1+M_{\mathrm{solvent}} \cdot \tilde{m}_{\mathrm{solute}} \cdot \sum_{i} \nu_{i, \text { ion }}} \tag{22}
\end{equation*}
$$

Here, $M_{\text {solvent }}$ is the molar mass of the solvent (in SI units), and $\tilde{m}_{\text {solute }}$ is the molality of the nondissociated solute. In this work, all MIACs are given in the molality scale.

## - MODEL PARAMETERS

The pure-component parameters for the ions and solvents were taken from the literature, with the exception of the parameters for $\mathrm{V}^{3+}$ and $\mathrm{VO}^{2+}$. These were fitted in this work to our own density data of $\mathrm{VCl}_{3}$ and $\mathrm{VOSO}_{4}$ in water (see

Chapter S4 in the Supporting Information for additional details). The new modeling approach used in this work requires two additional parameters: The Born diameter $d_{i}^{\text {Born }}$ of ions and the solvation parameter $f_{k}$ for solvents (set to 1 for all ions). The Born diameters were adjusted to Gibbs energies of solvation at infinite dilution in water (details on the used Gibbs energies of solvation in water are found in Chapter S5 in the Supporting Information), and the solvent-specific (ionindependent) $f_{k}$ parameters were adjusted to MIAC data of NaBr in each respective solvent. It can be observed that $f_{k}$ is close to 1.5 for the considered solvents. Indeed, we found that $f_{k}=1.5$ for the considered solvents if only electrolyte solutions with low ion concentrations are considered (ca. 10 mM ). The influence of $f_{k}$ on $\Delta^{\text {solv, } x} G_{i}^{\infty}$ is stronger with the smaller $d_{i}^{\text {Born }}$ caused by the dependency of $a^{\text {Born }}$ on $1 / d_{i}^{\text {Born }}$, cf. eq 7. Thus, the influence of $f_{k}$ on MIACs is much more pronounced than that on $\Delta^{\text {solv }, x} G_{i}^{\infty}$.

New binary interaction parameters $k_{i j}$ between the ions and the solvents and between all ions were also adjusted in this work. The $k_{i j}$ between an ion and water and the $k_{i j}$ among ions were adjusted to experimental MIAC literature data of different salts in water. In the special case of $\mathrm{V}^{3+}$ and $\mathrm{VO}^{2+}$, the $k_{i j}$ parameters $\mathrm{V}^{3+}$-water, $\mathrm{VO}^{2+}$-water, $\mathrm{V}^{3+}-\mathrm{Cl}^{-}$, and $\mathrm{VO}^{2+}-\mathrm{SO}_{4}{ }^{2-}$ were fitted to new osmotic coefficients of $\mathrm{VCl}_{3}$ and $\mathrm{VOSO}_{4}$ in water.

This order for fitting the different parameters ( $d_{i}^{\text {Born }}, f_{k}$ and $k_{i j}$ between ions and water) was chosen for the following reasons: The Gibbs energy of solvation of an ion in water is strongly influenced by the $d_{i}^{\text {Born }}$ parameter and only slightly changes with changing $f_{k}$ and $k_{i j}$. The $d_{i}^{\text {Born }}$ parameter was therefore fit to this quantity, while MIACs in water were used to adjust the other parameters. An iteration was then performed to get an optimal fit for all quantities. However, the change in all parameters during the iteration was minimal.

Gibbs energies of solvation at infinite dilution in organic solvents were calculated using the literature Gibbs energies of solvation in water and the literature Gibbs energies of transfer from water to the solvent. These data in the pure organic solvents were used to determine the $k_{i j}$ between ion/methanol and ion/ethanol. $k_{i j}$ values for solvent-solvent pairs were taken from the literature since the new modeling approach does not affect systems that only consist of neutral components. All pure-component parameters are given in Tables 2 and 3, and all binary interaction parameters used in this work are presented in Tables 4 and 5 . The values of the $k_{i j}$ parameters between the ions and the solvents deviate strongly from 0 , especially between an ion and an organic solvent. This strong deviation can be explained by the low sensitivity of the dispersion interactions compared to ion-dipolar forces

Table 2. Pure-Component Parameters of All Solvents Used in This Work ${ }^{a}$
| organic solvent | water ${ }^{63}$ | methanol ${ }^{64}$ | ethanol ${ }^{64}$ |
| :--- | :--- | :--- | :--- |
| $m_{i}^{\text {seg }}$ | 1.2047 | 1.5255 | 2.3827 |
| $\sigma_{i} / \AA$ | $b$ | 3.2300 | 3.1771 |
| $u_{i} / k_{\mathrm{B}} / \mathrm{K}$ | 353.95 | 188.90 | 198.24 |
| $\varepsilon^{\mathrm{AiBi}} / \mathrm{K}$ | 2425.7 | 2899.5 | 2653.4 |
| $\kappa^{\mathrm{AiBi}}$ | 0.04509 | 0.03518 | 0.03238 |
| $f_{k}$ | 1.5 | 1.4 | 1.6 |


${ }^{a}$ All solvents have a 2B association scheme. ${ }^{b} \sigma=2.7927+(10.11 \cdot \mathrm{e}^{-0.01775 \mathrm{~T} / \mathrm{K}}-1.417 \cdot \mathrm{e}^{-0.01146 \mathrm{~T} / \mathrm{K}}$ )

Table 3. ePC-SAFT Pure-Ion Parameters, $\sigma_{i}$ and $u_{i} / k_{B}$, Taken from ${ }^{a 17}$
| ion | $\sigma_{i} / \AA$ | $u_{i} / k_{\mathrm{B}} / \mathrm{K}$ | $d_{i}^{\text {Born }}(\AA)$ |
| :--- | :--- | :--- | :--- |
| $\mathrm{H}^{+}$ | 3.4654 | 500.00 | 1.218 |
| $\mathrm{Li}^{+}$ | 2.8449 | 360.00 | 2.784 |
| $\mathrm{Na}^{+}$ | 2.8232 | 230.00 | 3.445 |
| $\mathrm{K}^{+}$ | 3.3417 | 200.00 | 4.150 |
| $\mathrm{VO}^{2+}$ | 3.4638 | 1488.28 | 2.930 |
| $\mathrm{V}^{3+}$ | 3.1401 | 1415.66 | 2.940 |
| $\mathrm{Cl}^{-}$ | 2.7560 | 170.00 | 4.100 |
| $\mathrm{Br}^{-}$ | 3.0707 | 190.00 | 4.480 |
| $\mathrm{I}^{-}$ | 3.6672 | 200.00 | 4.985 |
| $\mathrm{SO}_{4}{ }^{2-}$ | 2.6491 | 80.00 | 4.953 |


${ }^{a}$ Segment number is 1 for all ions.

Table 4. ePC-SAFT Binary Interaction Parameters $k_{i j}$ between Different Anions and Cations Fitted in This Work Are Valid at 298.15 K
| $k_{\text {cat,an }}$ | $\mathrm{H}^{+}$ | $\mathrm{Li}^{+}$ | $\mathrm{Na}^{+}$ | $\mathrm{K}^{+}$ | $\mathrm{VO}^{2+}$ | $\mathrm{V}^{3+}$ |
| :--- | :---: | :--- | :--- | :--- | :---: | :---: |
| $\mathrm{Cl}^{-}$ | -0.9 | 0.8 | 0.8 | 0 |  | 1 |
| $\mathrm{Br}^{-}$ | -0.7 | 0.5 | 0.65 | -0.35 |  |  |
| $\mathrm{I}^{-}$ |  |  | 0.45 | 0 |  |  |
| $\mathrm{SO}_{4}{ }^{2-}$ |  | 0.8 | 0.7 | -0.6 | 0 |  |


Table 5. ePC-SAFT Binary Interaction Parameters $\boldsymbol{k}_{\boldsymbol{i j}}$ between Different Components at 298.15 K Used in This Work ${ }^{\text {a }}$
| $k_{i j}$ | water | methanol | ethanol |
| :--- | :--- | :--- | :--- |
| $\mathrm{Cl}^{-}$ | -0.3 | 0.5 | 0.8 |
| $\mathrm{Br}^{-}$ | -0.3 | 0.15 | 0 |
| $\mathrm{I}^{-}$ | -0.05 | 0.37 | 0.18 |
| $\mathrm{SO}_{4}{ }^{2-}$ | -0.1 |  |  |
| $\mathrm{H}^{+}$ | 0 | -0.3 | -0.6 |
| $\mathrm{Li}^{+}$ | -0.4 | -0.9 | -0.8 |
| $\mathrm{Na}^{+}$ | -0.3 | -0.25 | 0.05 |
| $\mathrm{K}^{+}$ | -0.1 | 0.32 | 0.53 |
| $\mathrm{VO}^{2+}$ | -0.08 |  |  |
| $\mathrm{V}^{3+}$ | -0.3 |  |  |
| methanol | $-0.0878^{65}$ |  |  |
| ethanol | -0.0617 ${ }^{66}$ |  |  |


${ }^{a}$ All $\mathrm{k}_{i j} s$ between ions and solvents were fitted in this work, and the $\mathrm{k}_{i j} s$ between different solvents were inherited from the literature. Valid at 298.15 K
described by the Born term, which was also found in ref 16 Therefore, high $k_{i j}$ are needed to have an effect on the modeling results.

## - GIBBS ENERGY OF SOLVATION IN WATER

The Gibbs energy of solvation at infinite dilution in water $\Delta^{\text {solv }} G_{i}^{\infty, \text { water }}$ was calculated by using the parameters from Tables 2 and 1. It is clear that the different mixing rules for the dielectric constant have no influence on the Gibbs energies of solvation since that is a property at infinite dilution. The results of $\Delta^{\text {solv }} G_{i}^{\infty, \text { water }}$ using the earlier version ePC-SAFT advanced and the new ePC-SAFT modeling approach are shown in Figure 4 together with literature data. The Gibbs energy of solvation at infinite dilution in water becomes less negative from $\mathrm{Li}^{+}$to $\mathrm{K}^{+}$and also from $\mathrm{Cl}^{-}$to $\mathrm{I}^{-}$. Larger ions have a lower charge density at the surface, which leads to weaker

![](https://cdn.mathpix.com/cropped/d6b90f00-389e-416f-b97c-c8651d44a864-07.jpg?height=442&width=1270&top_left_y=197&top_left_x=419)
Figure 4. Gibbs energy of solvation at infinite dilution in water of different ions at 298.15 K and 1 bar. Gray: Average value of different literature sources ${ }^{36,38,40,67-69}$ for $\mathrm{H}^{+}$only the value of Tissandier et al. ${ }^{69}$ is shown as recommended by Klinksiek et al.; ${ }^{70}$ blue: New modeling approach introduced in this work using parameters from Tables 1-5; and brown: ePC-SAFT advanced results taken from ref 16 All values are given in the mole fraction scale.

![](https://cdn.mathpix.com/cropped/d6b90f00-389e-416f-b97c-c8651d44a864-07.jpg?height=529&width=1270&top_left_y=825&top_left_x=419)
Figure 5. Molality-based MIACs of different salts in water at 298.15 K and 1 bar. Left: $\mathrm{Cl}^{-}$salts; right: $\mathrm{Br}^{-}$salts. Symbols: Literature data taken from. ${ }^{71}$ Purple upside-down triangles: $\mathrm{H}^{+}$components; orange triangles: $\mathrm{Li}^{+}$salts; green squares: $\mathrm{Na}^{+}$salts; and gray circles: $\mathrm{K}^{+}$salts. Lines : Calculated using ePC-SAFT with the modeling approach introduced in this work using parameters from Tables 1-5.

interactions and consequently to Gibbs energies of solvation closer to zero.

The calculations with the new modeling approach are in good agreement with the literature values for $\Delta^{\text {solv }} G_{i}^{\infty, \text { water }}$. The reason behind this is that $\Delta^{\text {solv }, x} G_{i}^{\infty}$ is sensitive to the Born diameter. ePC-SAFT advanced predicts reasonable $\Delta^{\text {solv }} G_{i}^{\infty, \text { water }}$ values for $\mathrm{Li}^{+}$but is not able to accurately describe this quantity for other ions. This is related to the use of segment diameters in the Born term in ePC-SAFT advanced since those differ greatly from the fitted Born diameters (with the exception of $\mathrm{Li}^{+}$). Rashin et al. and Schmid et al. also used fitted Born diameters in the same range as the fitted ones in this work. ${ }^{34,36}$ The differences in the exact values of the Born diameters can be attributed to the used literature Gibbs energies of solvation, and the fact that in this work, other contributions than the Born term are included in calculating this quantity. Fawcett et al. calculated the Gibbs energies of solvation of different ions in water and used Shannon-Prewitt crystal radii in the Born term, which are smaller than the fitted values in this work and the works of Rashin et al. and Schmid et al. ${ }^{38}$ This led to Gibbs energies of solvation that were much more negative than the literature values. Other authors used covalent radii or ionic (Pauling) radii, which also resulted in too negative Gibbs energies of solvation. ${ }^{33,34}$ Therefore, fitting Born diameters to Gibbs energies of solvation is justified, especially since the values in this work are comparable to the values of other authors. However, this again shows the poor definition of the Born diameter and raises questions on its
physical interpretation. The modeled Gibbs energies of solvation at infinite dilution in water for all investigated ions are found in the Supporting Information in Table S4.

## - MIACS IN WATER

Most applications that use electrolyte solutions operate at higher electrolyte concentrations and cannot be described using only infinite-dilution properties. ePC-SAFT advanced showed great results in modeling concentrated aqueous salt solutions and represents a good benchmark for the new modeling approach. It should again be noted that the aim of this work is the description of properties at infinite dilution while retaining the ability of ePC-SAFT to describe concentrated electrolyte solutions quantitatively correct. Figure 5 shows the MIACs of different salts in water. The new modeling approach is able to describe the initial decrease and the later increase of the MIACs with higher salt concentrations very well for the investigated salts. However, a re-evaluation of the binary interaction parameters compared to earlier versions of ePC-SAFT was necessary to achieve the quantitative agreement. In contrast, the new treatment did not require refitting the pure-ion parameters. Note, the $f_{k}$ parameter for water was only fitted to MIAC data of NaBr in water, and it remains constant for all other calculations.

## - PROPERTIES IN ORGANIC AND MIXED SOLVENTS

The next step in this work was the evaluation of infinitedilution properties of electrolytes in organic solvents. Methanol

![](https://cdn.mathpix.com/cropped/d6b90f00-389e-416f-b97c-c8651d44a864-08.jpg?height=990&width=1274&top_left_y=195&top_left_x=417)
Figure 6. Gibbs energies of transfer at infinite dilution from water to water + organic solvent systems $\Delta^{\text {trans }} G_{i}^{\infty}$ of different ions at 298.15 K and 1 bar. Top left: $\mathrm{K}^{+}$in water-methanol mixtures; top right: $\mathrm{Br}^{-}$in water-methanol mixtures; bottom left: $\mathrm{Na}^{+}$in water-ethanol mixtures; and bottom right: $\mathrm{Cl}^{-}$in water-ethanol mixtures. Lines: Calculated using ePC-SAFT with the approach introduced in this work using parameters from Tables 1-5; circles: Literature data taken from refs 73,74. All values are given on the mole fraction scale.

and ethanol were chosen as they represent common solvents, and there is more literature data on their electrolyte solutions available compared to other organic solvents. The Gibbs energy of transfer at infinite dilution from water to organic solvents is a very important property for transferring ions from an aqueous to an organic system. It represents the difference in Gibbs energy of solvation of an ion in the two solvents.

The new treatment of the Born term was then applied in the ePC-SAFT framework to predict Gibbs energies of transfer from water to mixed solvents and to predict MIACs of all electrolytes (except NaBr , this was used to fit the $f_{k}$ parameter) in single solvents and solvent mixtures without fitting additional parameters. Figure 6 shows the Gibbs energies of transfer at infinite dilution from water to mixed solvent systems. ePC-SAFT correlates well with the experimental values in the pure solvents. Further, ePC-SAFT correctly predicts increasing Gibbs energy of transfer from water to mixed solvents with increasing alcohol content. The dependence of $\Delta^{\text {trans }} G_{i}^{\infty}$ on the solvent composition was discussed by Hefter et al., who attributed the observed behavior to changes in the three-dimensional structure of water molecules with an increase in the alcohol content, which results in complex enthalpy and entropy changes. ${ }^{72}$ The latter two partially compensate each other, which leads to an apparently simple behavior of $\Delta^{\text {trans }} G_{i}^{\infty}$ with varying solvent compositions. We also checked whether previous versions of ePC-SAFT could describe these data, but this could not be confirmed, which is again caused by the strong influence of the Born diameters on Gibbs energies of solvation, where differences between modeling and experimental data cannot be overcome completely by only fitting $k_{i j}$ parameters between ion and solvent.

The properties of concentrated solutions in organic and mixed solvents were also investigated to evaluate the performance of the new modeling approach in these systems. First, experimental MIACs of NaBr in methanol and in ethanol were used to adjust the $f_{k}$ parameter of methanol and ethanol. A comparison of the MIACs of NaBr in methanol with the modeling results is presented in Figure 7. It can be observed that ePC-SAFT modeling MIACs in methanol is quantitatively correct using the mixing rule for the dielectric constant according to eq 16. In contrast, modeling using the mole

![](https://cdn.mathpix.com/cropped/d6b90f00-389e-416f-b97c-c8651d44a864-08.jpg?height=527&width=636&top_left_y=1833&top_left_x=1194)
Figure 7. Molality-based MIACs of NaBr in methanol at 298.15 K and 1 bar. Circles: Literature data from ref 77; solid line: ePC-SAFT with the modeling approach introduced in this work using parameters from Tables 1-5; and dashed line: ePC-SAFT with the modeling approach in this work but with the mole fraction-based linear mixing rule for the dielectric constant.

![](https://cdn.mathpix.com/cropped/d6b90f00-389e-416f-b97c-c8651d44a864-09.jpg?height=1070&width=1291&top_left_y=191&top_left_x=406)
Figure 8. Molality-based MIACs of LiBr (top left), NaI (top right), and NaBr (bottom) in methanol and in ethanol at 298.15 K and 1 bar. Gray circles: Literature MIAC data in methanol; ${ }^{75-78}$ green squares: Literature MIAC data in ethanol; ${ }^{77,79-81}$ lines: ePC-SAFT with the modeling approach introduced in this work using parameters from Tables 1-5; gray line: Methanol; and black line: Ethanol. The MIACs were modeled up to the solubility limit of the salts in the respective solvent as measured in refs 82,83 .

fraction-based linear mixing rule for the dielectric constant yielded high deviations between ePC-SAFT and the experimental MIACs.

As a next step, the MIACs of all other salts in ethanol and methanol were modeled without fitting further parameters since the binary interaction parameters ion-methanol and ionethanol were fitted to Gibbs energies of solvation at infinite dilution. The MIACs of $\mathrm{LiBr}, \mathrm{NaI}$, and NaBr in pure methanol and pure ethanol are shown in Figure 8. Two different sources for LiBr MIACs in methanol were found in the literature. ${ }^{75,76}$ The ePC-SAFT-predicted MIACs agree well with these two data sources up to a molality of about $3.5 \mathrm{~mol} \mathrm{~kg}^{-1}$. At molalities greater than $3.5 \mathrm{~mol} \mathrm{~kg}^{-1}$, ePC-SAFT underpredicts the experimental MIACs, which is quite obvious for LiBr in ethanol. These deviations might be explained by the fact that the binary interaction parameters ion-methanol and ionethanol were fitted to infinite-dilution data. Note that electrolyte solubility in organic solvents is usually not very high, so a low-to-moderate concentration range is more important.

Further, the ability of ePC-SAFT to predict MIACs in mixed solvent systems was evaluated. Figure 9 shows selected MIACs in water-methanol and water-ethanol mixtures at different solvent compositions. Again, the new modeling approach yields good results for the investigated systems, especially at lower electrolyte concentrations. The modeling results for the system water-ethanol- NaBr show deviations at ion molalities greater than about $2 \mathrm{~mol} \mathrm{~kg}^{-1}$ at salt-free weight fractions of ethanol equal to 0.4 . ePC-SAFT underestimates the increase in the MIACs with electrolyte molality at molalities greater than
about $1 \mathrm{~mol} \mathrm{~kg}^{-1}$. This could be related to the fact that only data at low electrolyte concentrations were used in fitting the parameters for the NaBr -ethanol system (infinite-dilution properties were used for $k_{i j}$, and MIAC data were used up to $0.25 \mathrm{~mol} \mathrm{~kg}^{-1}$ for $f_{\text {ethanol }}$ ).

To sum up, the new modeling approach within the ePCSAFT framework allows accurate modeling of properties at infinite dilution and in concentrated electrolyte single-solvent systems and even prediction in mixed-solvent electrolyte systems. The pure-ion parameters could be inherited from earlier works, but the binary interaction parameters needed to be readjusted. The improved performance of ePC-SAFT using our new approach vs recent modeling approaches for infinitely diluted systems was found to be in the new treatment of Born diameters. The good performance of the model for concentrated systems can be attributed to the use of a new mixing rule for the dielectric constant in combination with this new treatment of Born diameter.

## - CONCLUSIONS

The goal of this work was to increase the accuracy of ePCSAFT toward modeling properties at infinite dilution in single solvents and in mixed solvents at 298.15 K while keeping the well-known accuracy of describing properties at more concentrated electrolyte solutions. To achieve this goal, ePCSAFT was further improved by accounting for the solvationshell thickness and the dielectric saturation effect within the Born term. These modifications required the adjustment of the Born diameter of each ion and one ion-independent solvation parameter per solvent as well as one binary interaction

![](https://cdn.mathpix.com/cropped/d6b90f00-389e-416f-b97c-c8651d44a864-10.jpg?height=1053&width=1276&top_left_y=195&top_left_x=417)
Figure 9. Molality-based MIACs of NaBr (top) and NaCl (bottom) in methanol-water (left) and ethanol-water (right) mixtures at 298.15 K and 1 bar. Orange triangles: Literature MIACs in $w_{\text {organic solvent }}^{\text {salt-fre }}=0.8$ taken from refs $77,84,85$; green squares: Literature MIACs in $w_{\text {organic solvent }}^{\text {salt-fre }}=0.4$ taken from refs $77,84,85$. Lines are modeling results with ePC-SAFT from this work using parameters from Tables $1-5$. The MIACs were modeled up to the solubility of the salts in the respective solvent mixture as measured in refs 82,83 .

parameter ion-solvent. It was shown that infinite-dilution properties can be used to adjust the Born diameter, while binary interaction parameters were used to correlate finiteconcentration properties (MIACs in this work). All other pureion parameters were taken from earlier versions of ePC-SAFT without change. New parameters were determined for vanadium ions $\mathrm{V}^{3+}$ and $\mathrm{VO}^{2+}$ by fitting to density data and to osmotic-coefficient data measured in this work. Additionally, a new mixing rule for the dielectric constant was introduced by correlating experimental dielectric constant data of single electrolytes in water and in methanol. The mixing rule was then validated and holds true also for single electrolytes in ethanol, which was not used for determining the correlation. The new treatment of the Born term allowed for a better description of properties of ions at infinite dilution in pure and mixed solvents, which was shown by comparing modeled Gibbs energies of solvation and transfer at infinite dilution with literature data. The new modeling approach, including the new mixing rule for the dielectric constant, retained the ability of earlier versions to describe properties of concentrated electrolyte solutions in pure solvents, which was shown on MIACs in water, methanol, and ethanol. Even more, MIACs in mixed solvents were predicted accurately up to molalities of about $2 \mathrm{~mol} \mathrm{~kg}^{-1}$. This is a significant advance over previous modeling works and should also be helpful for other electrolyte thermodynamic models.

- ASSOCIATED CONTENT


## (5) Supporting Information

The Supporting Information is available free of charge at https://pubs.acs.org/doi/10.1021/acs.iecr.5c00475.

Details on different approaches to calculate the solvation shell diameter; comparison of the performance of different Born contributions in ePC-SAFT together with different mixing rules for the dielectric constant; equations for the conversion of the Gibbs energy of solvation for different standard treatments; density and osmotic coefficient data for aqueous $\mathrm{VCl}_{3}$ and $\mathrm{VOSO}_{4}$ solutions; comparison of literature data of the Gibbs energy of solvation of different ions in water (PDF)

## - AUTHOR INFORMATION

## Corresponding Authors

Gangqiang Yu - Laboratory of Thermodynamics, Department of Biochemical and Chemical Engineering, TU Dortmund, Dortmund 44277, Germany; College of Environmental science and Engineering, Beijing University of Technology, Beijing 100124, China; © orcid.org/0000-0002-35956972; Email: yugq@bjut.edu.cn
Christoph Held - Laboratory of Thermodynamics, Department of Biochemical and Chemical Engineering, TU Dortmund, Dortmund 44277, Germany; © orcid.org/0000-0003-1074-177X; Email: christoph.held@tu-dortmund.de

## Author

Paul Figiel - Laboratory of Thermodynamics, Department of Biochemical and Chemical Engineering, TU Dortmund,

Dortmund 44277, Germany; orcid.org/0000-0002-59477247
Complete contact information is available at:
https://pubs.acs.org/10.1021/acs.iecr.5c00475

## Notes

The authors declare no competing financial interest.

## - ACKNOWLEDGMENTS

This work was sponsored by the Beijing Nova Program (No. 20240484727) and financially supported by the National Natural Science Foundation of China under Grants (No. 22378006). The authors acknowledge funding from the Alexander von Humboldt Foundation (G.Y.). The work was funded by the German Research Foundation (DFG), HE716515/1, Project nr. 525252957.

## REFERENCES

(1) Adelekan, O. A.; Ilugbusi, B. S.; Adisa, O.; Obi, O. C.; Awonuga, K. F.; Asuzu, O. F.; Ndubuisi, N. L. Energy transition policies: a global review of shifts towards renewable sources. Eng. Sci. Technol. J. 2024, 5 (2), 272-287.
(2) Hassan, Q.; Viktor, P.; Al-Musawi, T. J.; Ali, B. M.; Algburi, S.; Alzoubi, H. M.; Al-Jiboory, A. K.; Sameen, A. Z.; Salman, H. M.; Jaszczur, M. The renewable energy role in the global energy Transformations. Renewable Energy Focus 2024, 48, No. 100545.
(3) Singh, A. N.; Islam, M.; Meena, A.; Faizan, M.; Han, D.; Bathula, C.; Hajibabaei, A.; Anand, R.; Nam, K.-W. Unleashing the Potential of Sodium-Ion Batteries: Current State and Future Directions for Sustainable Energy Storage. Adv. Funct. Mater. 2023, 33 (46), No. 2304617.
(4) Wei, P.; Abid, M.; Adun, H.; Awoh, D. K.; Cai, D.; Zaini, J. H.; Bamisile, O. Progress in Energy Storage Technologies and Methods for Renewable Energy Systems Application. Appl. Sci. 2023, 13 (9), No. 5626.
(5) Kulova, T. L.; Fateev, V. N.; Seregina, E. A.; Grigoriev, A. S. A Brief Review of Post-Lithium-Ion Batteries. Int. J. Electrochem. Sci. 2020, 15 (8), 7242-7259.
(6) Chalamala, B. R.; Soundappan, T.; Fisher, G. R.; Anstey, M. R.; Viswanathan, V. V.; Perry, M. L. Redox Flow Batteries: An Engineering Perspective. Proc. IEEE 2014, 102 (6), 976-999.
(7) Khan, F. N. U.; Rasul, M. G.; Sayem, A.; Mandal, N. Maximizing energy density of lithium-ion batteries for electric vehicles: A critical review. Energy Rep. 2023, 9, 11-21.
(8) Khan, F. M. N. U.; Rasul, M. G.; Sayem, A.; Mandal, N. K. Design and optimization of lithium-ion battery as an efficient energy storage device for electric vehicles: A comprehensive review. J. Energy Storage 2023, 71, No. 108033.
(9) Hosaka, T.; Kubota, K.; Hameed, A. S.; Komaba, S. Research Development on K-Ion Batteries. Chem. Rev. 2020, 120 (14), 63586466.
(10) Haupt, A.; Lerch, A. Forward Osmosis Application in Manufacturing Industries: A Short Review. Membranes 2018, 8 (3), No. 47.
(11) Wang, Y.; Feric, T. G.; Tang, J.; Fang, C.; Hamilton, S. T.; Halat, D. M.; Wu, B.; Celik, H.; Rim, G.; DuBridge, T.; Oshiro, J.; Wang, R.; Park, A.-H. A.; Reimer, J. A. Carbon capture in polymerbased electrolytes. Sci. Adv. 2024, 10 (16), No. eadk2350.
(12) Liu, Y.; Wang, Y.; Zhao, S. Journey of electrochemical chlorine production: From brine to seawater. Curr. Opin. Electrochem 2023, 37, No. 101202.
(13) Xia, R.; Overa, S.; Jiao, F. Emerging Electrochemical Processes to Decarbonize the Chemical Industry. JACS Au 2022, 2 (5), 10541070.
(14) McEnaney, J. M.; Singh, A. R.; Schwalbe, J. A.; Kibsgaard, J.; Lin, J. C.; Cargnello, M.; Jaramillo, T. F.; Nørskov, J. K. Ammonia synthesis from N 2 and H 2 O using a lithium cycling electrification
strategy at atmospheric pressure. Energy Environ. Sci. 2017, 10 (7), 1621-1630.
(15) Andersen, S. Z.; Statt, M. J.; Bukas, V. J.; Shapel, S. G.; Pedersen, J. B.; Krempl, K.; Saccoccio, M.; Chakraborty, D.; Kibsgaard, J.; Vesborg, P. C. K.; Nørskov, J.; Chorkendorff, I. Increasing stability, efficiency, and fundamental understanding of lithium-mediated electrochemical nitrogen reduction. Energy Environ. Sci. 2020, 13 (11), 4291-4300.
(16) Bülow, M.; Ascani, M.; Held, C. ePC-SAFT advanced - Part I: Physical meaning of including a concentration-dependent dielectric constant in the born term and in the Debye-Hückel theory. Fluid Phase Equilib. 2021, 535, No. 112967.
(17) Held, C.; Reschke, T.; Mohammad, S.; Luza, A.; Sadowski, G. ePC-SAFT revised. Chem. Eng. Res. Des. 2014, 92 (12), 2884-2897.
(18) Silva, G. M.; Maribo-Mogensen, B.; Liang, X.; Kontogeorgis, G. M. Improving the Born equation: Origin of the Born radius and introducing dielectric saturation effects. Fluid Phase Equilib. 2024, 576, No. 113955.
(19) Friese, E.; Ebel, A. Temperature dependent thermodynamic model of the system $\mathrm{H}(+)-\mathrm{NH}_{4}(+)-\mathrm{Na}(+)-\mathrm{SO}_{4}{ }^{2-}-\mathrm{NO}_{3}{ }^{-}-\mathrm{Cl}^{-}-\mathrm{H}_{2} \mathrm{O} . J$. Phys. Chem. A 2010, 114 (43), 11595-11631.
(20) Pitzer, K. S. Thermodynamics of electrolytes. I. Theoretical basis and general equations. J. Phys. Chem. A 1973, 77 (2), 268-277.
(21) Chen, C.-C.; Britt, H. I.; Boston, J. F.; Evans, L. B. Local composition model for excess Gibbs energy of electrolyte systems. Part I: Single solvent, single completely dissociated electrolyte systems. AIChE J. 1982, 28 (4), 588-596.
(22) Kontogeorgis, G. M.; Maribo-Mogensen, B.; Thomsen, K. The Debye-Hückel theory and its importance in modeling electrolyte solutions. Fluid Phase Equilib. 2018, 462, 130-152.
(23) Kontogeorgis, G. M.; Schlaikjer, A.; Olsen, M. D.; MariboMogensen, B.; Thomsen, K.; von Solms, N.; Liang, X. A Review of Electrolyte Equations of State with Emphasis on Those Based on Cubic and Cubic-Plus-Association (CPA) Models. Int. J. Thermophys. 2022, 43 (4), No. 54.
(24) Chen, C.-C.; Song, Y. Generalized electrolyte-NRTL model for mixed-solvent electrolyte systems. AIChE J. 2004, 50 (8), 1928-1941.
(25) Simonin, J.-P. On the "Born" term used in thermodynamic models for electrolytes. J. Chem. Phys. 2019, 150 (24), No. 244503.
(26) Bülow, M.; Ascani, M.; Held, C. ePC-SAFT advanced - Part II: Application to Salt Solubility in Ionic and Organic Solvents and the Impact of Ion Pairing. Fluid Phase Equilib. 2021, 537, No. 112989.
(27) Pabsch, D.; Held, C.; Sadowski, G. Modeling the CO 2 Solubility in Aqueous Electrolyte Solutions Using ePC-SAFT. J. Chem. Eng. Data 2020, 65 (12), 5768-5777.
(28) Pabsch, D.; Figiel, P.; Sadowski, G.; Held, C. Solubility of Electrolytes in Organic Solvents: Solvent-Specific Effects and IonSpecific Effects. J. Chem. Eng. Data 2022, 67 (9), 2706-2718.
(29) Pabsch, D.; Lindfeld, J.; Schwalm, J.; Strangmann, A.; Figiel, P.; Sadowski, G.; Held, C. Influence of solvent and salt on kinetics and equilibrium of esterification reactions. Chem. Eng. Sci. 2022, 263, No. 118046.
(30) Klinksiek, M.; Baco, S.; Leveneur, S.; Legros, J.; Held, C. Activity-Based Models to Predict Kinetics of Levulinic Acid Esterification. ChemPhysChem 2023, 24 (4), No. e202200729.
(31) Maribo-Mogensen, B.; Thomsen, K.; Kontogeorgis, G. M. An electrolyte CPA equation of state for mixed solvent electrolytes. AIChE J. 2015, 61 (9), 2933-2950.
(32) Born, M. Volumen und Hydratationswärme der Ionen. Z. Phys. 1920, 1 (1), 45-48.
(33) Latimer, W. M.; Pitzer, K. S.; Slansky, C. M. The Free Energy of Hydration of Gaseous Ions, and the Absolute Potential of the Normal Calomel Electrode. J. Chem. Phys. 1939, 7 (2), 108-111.
(34) Rashin, A. A.; Honig, B. Reevaluation of the Born model of ion hydration. J. Phys. Chem. A 1985, 89 (26), 5588-5593.
(35) Latimer, W. M. Single Ion Free Energies and Entropies of Aqueous Ions. J. Chem. Phys. 1955, 23 (1), 90-92.
(36) Schmid, R.; Miah, A. M.; Sapunov, V. N. A new table of the thermodynamic quantities of ionic hydration: values and some
applications (enthalpy-entropy compensation and Born radii). Phys. Chem. Chem. Phys. 2000, 2 (1), 97-102.
(37) KJ Laidler, C. P. The influence of dielectric saturation on the thermodynamic properties of aqueous ions. Proc. R. Soc. London, Ser. A 1957, 241 (1224), 80-92.
(38) Fawcett, W. R. Thermodynamic Parameters for the Solvation of Monatomic Ions in Water. J. Phys. Chem. B 1999, 103 (50), 1118111185.
(39) Kumar, A. A Modified Born Equation for Solvation Energy of Ions. J. Phys. Soc. Jpn. 1992, 61 (11), 4247-4250.
(40) Marcus, Y. Thermodynamics of solvation of ions. Part 5.Gibbs free energy of hydration at 298.15 K. J. Chem. Soc., Faraday Trans. 1991, 87 (18), 2995-2999.
(41) Hasted, J. B.; Ritson, D. M.; Collie, C. H. Dielectric Properties of Aqueous Ionic Solutions. Parts I and II. J. Chem. Phys. 1948, 16 (1), 1-21.
(42) Duan, X.; Nakamura, I. A new lattice Monte Carlo simulation for dielectric saturation in ion-containing liquids. Soft Matter 2015, 11 (18), 3566-3571.
(43) Booth, F. Dielectric Constant of Polar Liquids at High Field Strengths. J. Chem. Phys. 1955, 23 (3), 453-457.
(44) Andeen, C.; Fontanella, J.; Schuele, D. Low-Frequency Dielectric Constant of $\mathrm{LiF}, \mathrm{NaF}, \mathrm{NaCl}, \mathrm{NaBr}, \mathrm{KCl}$, and KBr by the Method of Substitution. Phys. Rev. B 1970, 2 (12), No. 5068.
(45) Abraham, M. H.; Liszi, J.; Mészaros, L. Calculations on ionic solvation. III. The electrostatic free energy of solvation of ions, using a multilayered continuum model. J. Chem. Phys. 1979, 70 (5), 24912496.
(46) Abraham, M. H.; Liszi, J.; Kristof, E. Calculations on ionic solvation. VII. The free energy of solvation of ions calculated from various local solvent dielectric constant-distance functions. Aust. J. Chem. 1982, 35 (7), 1273-1279.
(47) Floriano, W. B.; Nascimento, M. A. C. Dielectric constant and density of water as a function of pressure at constant temperature. Braz. J. Phys. 2004, 34 (1), 38-41.
(48) Khimenko, M. T.; Aleksandrov, V. V.; Gritsenko, N. N. Polarizability and radii of molecules of some pure liquids. Zh. Fiz. Khim. 1973, No. 47, 2914-2915.
(49) Shirke, R. M.; Chaudhari, A.; More, N. M.; Patil, P. B. Temperature dependent dielectric relaxation study of ethyl acetate Alcohol mixtures using time domain technique. J. Mol. Liq. 2001, 94 (1), 27-36.
(50) Bülow, M.; Ji, X.; Held, C. Incorporating a concentrationdependent dielectric constant into ePC-SAFT. An application to binary mixtures containing ionic liquids. Fluid Phase Equilib. 2019, 492, 26-33.
(51) Ascani, M.; Held, C. Prediction of salting-out in liquid-liquid two-phase systems with ePC-SAFT: Effect of the Born term and of a concentration-dependent dielectric constant. Z. Anorg. Allg. Chem. 2021, 647 (12), 1305-1314.
(52) Zuber, A.; Cardozo-Filho, L.; Cabral, V. F.; Checoni, R. F.; Castier, M. An empirical equation for the dielectric constant in aqueous and nonaqueous electrolyte mixtures. Fluid Phase Equilib. 2014, 376, 116-123.
(53) Nörtemann, K.; Hilland, J.; Kaatze, U. Dielectric Properties of Aqueous NaCl Solutions at Microwave Frequencies. J. Phys. Chem. A 1997, 101 (37), 6864-6869.
(54) Wei, Y.-Z.; Sridhar, S. Dielectric spectroscopy up to 20 GHz of LiCl/H2O solutions. J. Chem. Phys. 1990, 92 (2), 923-928.
(55) Haggis, G. H.; Hasted, J. B.; Buchanan, T. J. The Dielectric Properties of Water in Solutions. J. Chem. Phys. 1952, 20 (9), 14521465.
(56) Buchner, R.; Hefter, G. T.; May, P. M. Dielectric Relaxation of Aqueous NaCl Solutions. J. Phys. Chem. A 1999, 103 (1), 1-9.
(57) Kaatze, U.; Adolph, D.; Gottlob, D.; Pottel, R. Static Permittivity and Dielectric Relaxation of Solutions of Ions in Methanol. Ber. Bunsenges. Phys. Chem. 1980, 84 (12), 1198-1203.
(58) Winsor, P.; Cole, R. H. Dielectric properties of electrolyte solutions. 2. Alkali halides in methanol. J. Phys. Chem. A 1982, 86 (13), 2491-2494.
(59) Barthel, J.; R, N.; Chemistry Data Series, 12,2a. Electrolyte Data Collection, 2a: Dielectric Properties of Nonaqueous Electrolyte Solutions; DECHEMA, 1996.
(60) Hasted, J. B.; Roderick, G. W. Dielectric Properties of Aqueous and Alcoholic Electrolytic Solutions. J. Chem. Phys. 1958, 29 (1), 1726.
(61) Weingärtner, H.; Nadolny, H. G.; Käshammer, S. Dielectric Properties of an Electrolyte Solution at Low Reduced Temperature. J. Phys. Chem. B 1999, 103 (22), 4738-4743.
(62) Ben-Naim, A. Y. Solvation Thermodynamics; Springer, 1987.
(63) Cameretti, L. F.; Sadowski, G. Modeling of aqueous amino acid and polypeptide solutions with PC-SAFT. Chem. Eng. Process.: Process Intensif. 2008, 47 (6), 1018-1025.
(64) Gross, J.; Sadowski, G. Application of the Perturbed-Chain SAFT Equation of State to Associating Systems. Ind. Eng. Chem. Res. 2002, 41 (22), 5510-5515.
(65) Held, C.; Sadowski, G. Manuel 'e-PC-SAFT 1.0'; TU Dortmund University, Laboratory of Thermodynamics, 2017.
(66) Dohrn, S.; Reimer, P.; Luebbert, C.; Lehmkemper, K.; Kyeremateng, S. O.; Degenhardt, M.; Sadowski, G. Thermodynamic Modeling of Solvent-Impact on Phase Separation in Amorphous Solid Dispersions during Drying. Mol. Pharmaceutics 2020, 17 (7), 27212733.
(67) Conway, B. E. The evaluation and use of properties of individual ions in slution. J. Solution Chem. 1978, 7 (10), 721-770.
(68) Pliego, J. R., Jr; Riveros, J. M. Gibbs energy of solvation of organic ions in aqueous and dimethyl sulfoxide solutions. Phys. Chem. Chem. Phys. 2002, 4 (9), 1622-1627.
(69) Tissandier, M. D.; Cowen, K. A.; Feng, W. Y.; Gundlach, E.; Cohen, M. H.; Earhart, A. D.; Coe, J. V.; Tuttle, T. R. The Proton's Absolute Aqueous Enthalpy and Gibbs Free Energy of Solvation from Cluster-Ion Solvation Data. J. Phys. Chem. A 1998, 102 (40), 77877794.
(70) Klinksiek, M.; Baco, S.; Leveneur, S.; Legros, J.; Held, C. Predicting Proton Activity and Acid Dissociation Equilibria in MixedSolvent Systems, and Their Impact on Esterification Kinetics of Levulinic Acid. Ind. Eng. Chem. Res. 2025, 64 (11), 6188-6202.
(71) Hamer, W. J.; Wu, Y.-C. Osmotic Coefficients and Mean Activity Coefficients of Uni-univalent Electrolytes in Water at $25^{\circ} \mathrm{C} . J$. Phys. Chem. Ref. Data 1972, 1 (4), 1047-1100.
(72) Hefter, G. T.; Mclay, P. J. Solvation of Fluoride Ions. II. Enthalpies and Entropies of Transfer From Water to Aqueous Methanol. Aust. J. Chem. 1988, 41 (12), No. 1971.
(73) Marcus, Y. Gibbs energies of transfer of anions from water to mixed aqueous organic solvents. Chem. Rev. 2007, 107 (9), 38803897.
(74) Kalidas, C.; Hefter, G.; Marcus, Y. Gibbs energies of transfer of cations from water to mixed aqueous organic solvents. Chem. Rev. 2000, 100 (3), 819-852.
(75) Barthel, J.; Lauermann, G.; Neueder, R. Vapor pressure measurements on non-aqueous electrolyte solutions. Part 2. Tetraalkylammonium salts in methanol. Activity coefficients of various $1-1$ electrolytes at high concentrations. J. Solution Chem. 1986, 15 (10), 851-867.
(76) Safarov, J. T. Study of thermodynamic properties of binary solutions of lithium bromide or lithium chloride with methanol. Fluid Phase Equilib. 2005, 236 (1-2), 87-95.
(77) Ye, S.; Xans, P.; Lagourette, B. Modification of the Pitzer model to calculate the mean activity coefficients of electrolytes in a wateralcohol mixed solvent solution. J. Solution Chem. 1994, 23 (12), 1301-1315.
(78) Nasirzadeh, K.; Papaiconomou, N.; Neueder, R.; Kunz, W. Vapor Pressures, Osmotic and Activity Coefficients of Electrolytes in Protic Solvents at Different Temperatures. 1. Lithium Bromide in Methanol. J. Solution Chem. 2004, 33 (3), 227-245.
(79) Barthel, J.; Lauermann, G. Vapor pressure measurements on non-aqueous electrolyte solutions. Part 3: Solutions of sodium lodide in ethanol, 2-propanol, and acetonitrile. J. Solution Chem. 1986, 15 (10), 869-877.
(80) Held, C.; Prinz, A.; Wallmeyer, V.; Sadowski, G. Measuring and modeling alcohol/salt systems. Chem. Eng. Sci. 2012, 68 (1), 328339.
(81) Safarov, J. T. Vapor pressures of lithium bromide or lithium chloride and ethanol solutions. Fluid Phase Equilib. 2006, 243 (1-2), 38-44.
(82) Pinho, S. P.; Macedo, E. A. Solubility of $\mathrm{NaCl}, \mathrm{NaBr}$, and KCl in Water, Methanol, Ethanol, and Their Mixed Solvents. J. Chem. Eng. Data 2005, 50 (1), 29-32.
(83) Li, M.; Constantinescu, D.; Wang, L.; Mohs, A.; Gmehling, J. Solubilities of $\mathrm{NaCl}, \mathrm{KCl}, \mathrm{LiCl}$, and LiBr in Methanol, Ethanol, Acetone, and Mixed Solvents and Correlation Using the LIQUAC Model. Ind. Eng. Chem. Res. 2010, 49 (10), 4981-4988.
(84) Yan, W.; Xu, Y.; Han, S. Activity coefficients of sodium chloride in methanol-water mixed solvents at 298.15 K. Acta Chim. Sin. 1994, 52 (10), 937-946.
(85) Esteso, M. A.; Gonzalez-Diaz, O. M.; Hernandez-Luis, F. F.; Fernandez-Merida, L. Activity coefficients for NaCl in ethanol-water mixtures at 25C. J. Solution Chem. 1989, 18 (3), 277-288.


[^0]:    Received: February 3, 2025
    Revised: April 10, 2025
    Accepted: April 14, 2025
    Published: April 28, 2025

