# Prediction of salting-out in liquid-liquid two-phase systems with ePC-SAFT: Effect of the Born term and of a concentration-dependent dielectric constant 

Moreno Ascani ${ }^{[\mathrm{a}]}$ and Christoph Held ${ }^{*[\mathrm{a}]}$<br>We dedicate this to honor the $60^{\text {th }}$ birthday of Prof. Christoph Janiak for outstanding networking activities and great research.

Knowledge on phase equilibria is of crucial importance in designing industrial processes. However, modeling phase equilibria in liquid-liquid two-phase systems (LLTPS) containing electrolytes is still a challenge for electrolyte thermodynamic models and modeling still requires a lot of experimental input data. Further, modeling electrolyte solutions requires accounting for different physical effects in the electrolyte theory, especially the change of the dielectric properties of the medium at different compositions and the related change of solvation free energy of the dissolved ions. In a previous work, the Born term was altered by combining it with a concentration-dependent dielectric constant within the framework of electrolyte Perturbed-Chain Statistical Associating Fluid Theory (ePC-SAFT),
and hence called 'ePC-SAFT advanced'. In the present work, ePC-SAFT advanced was validated against liquid-liquid equilibria (LLE) of LLTPS water + organic solvents + alkali halides as well as aqueous two-phase systems containing the phase formers poly (propylene glycol) and an ionic liquid. All the ePCSAFT parameters were used as published in the literature, and each binary interaction parameter between ion-solvent was set to zero. ePC-SAFT advanced allowed quantitatively predicting the salt effect on LLTPS without adjusting binary interaction parameters, while classical ePC-SAFT or meaningless mixing rules for the dielectric constant term failed in predicting the phase behavior of the LLTPS.

## Introduction

Mixed-solvent aqueous electrolyte solutions or even nonaqueous electrolyte solutions are gaining increasing importance for innovative industrial applications. Non-aqueous electrolyte solutions are used as working fluids in reverse electrodialysis, ${ }^{[1]}$ for the production of rechargeable lithium-based batteries ${ }^{[2]}$ and as substitute of aqueous salt solutions for the production of absorption refrigeration machines and heat pumps. ${ }^{[3]}$ Furthermore, feedstocks originating from fermentation processes (fermentation broth) for the production of biomassderived chemical products are aqueous-based; however, they usually contain a huge number of different organic components and salts, where the latter can be already present in the system or might be intentionally added to change the phase equilibrium. ${ }^{[4-6]}$ Knowledge of liquid-liquid phase equilibria (LLE) in presence of salts and mixed solvents is crucial for designing liquid-liquid two-phase systems (LLTPS) used for separation processes to extract the target component from the fermenta-

[^0]tion broth. Knowledge and modelling of the phase behavior of mixed-solvent electrolyte solutions is also important for the design of salt-containing extractive distillation ${ }^{[7,8]}$ and extractive crystallization ${ }^{[8,9]}$ processes.

Salting-in or salting-out effects on macromolecules in solution are historically assessed based on the so-called 'Hofmeister series ${ }^{[10]}$. According to Hofmeister ions with a very high charge density will tend to induce much stronger saltingout effects (i.e., they will displace macromolecules) than ions with a low charge density, whereas presence of some salts will result in salting-in (i.e. they show stabilization effects toward macromolecules in solution). The usual classification divides the salts in "kosmotropes" and "chaotropes". That is, there are salts that tend to strengthen structure in the water thereby inducing salting-out, whereas others tend to disrupt the water structure thereby inducing salting-in. ${ }^{[10,11]}$

Among the LLTPS, there are several types of possible mixtures such as aqueous/organic mixtures. Further, there are aqueous two-phase systems (ATPS), which contain a high amount of water (usually more than $40 \%$ ) in each of the two liquid phases. ${ }^{[12]}$ ATPS are usually formed by combining aqueous systems of two polymers, of one polymer and one salt, ${ }^{[13]}$ or even of two salts. ${ }^{[14,15]}$ Liquid-liquid phase separation can be induced, for instance, by introducing a salting-out agent into a water-polymer solution, or by mixing two aqueous solutions of incompatible polymers, or even adding an ionic liquid (IL) to such a solution. Because of the high water content within ATPS, they are in general biocompatible media for biologically active substances ${ }^{[14]}$ and therefore have been proposed as substitute of traditional aqueous/organic LLTPS for
the recovery and purification of biomolecules. ${ }^{[16-18]}$ Most of the used ATPS phase formers are nonvolatile and biodegradable, and thus environmental friendlier than traditional solvents. Poly (alkylene glycol)s such as polyethylene glycol (PEG) or polypropylene glycol (PPG) are widely used as phase formers due to their biodegradability, low volatility, low viscosity and relative low cost. ${ }^{[19,20]}$ Several advantages are reported upon substituting inorganic salts with IL as salting-out agents in polymer-salt based ATPS. Due to the possibility to design ILs with defined properties, polymer-IL based ATPS can be optimized to meet specific requirements. The polarity of the IL-rich phase can be adjusted to tune the partition coefficient of selected biomolecules. Addition of IL to ternary systems polymer-salt-water was also studied. ${ }^{[21]}$ Since many ILs are liquid at room temperature, such ATPS do not suffer from the risk of IL crystallization. ${ }^{[14]}$ Depending on the chosen IL and the type of polymer, ILs can act either as salting-in agents (they promote dissolution) or salting-out agents (they promote phase separation). Because of the large number of polymer-IL combinations, there is an urgent need for reliable predictive models which can help reducing the experimental effort needed to design such systems.

Unfortunately, while the development of processes mainly based on aqueous electrolyte systems can be considered as well-optimized, water-poor electrolyte systems still require a lot of data and of understanding thermo-physical behavior at different conditions (solvent composition, kind and concentration of dissolved electrolytes). Modeling the phase behavior of electrolyte systems has been performed by different groups; these already included the Born theory into thermodynamic models e.g., the electrolyte Cubic-Plus-Association (eCPA), ${ }^{[22-24]}$ Peng-Robinson ${ }^{[25,26]}$ and Soave-Redlich-Kwong ${ }^{[22]}$ models and electrolyte Perturbed-Chain Statistical Associating Fluid Theory (ePC-SAFT) ${ }^{[27-42]}$ as well as the most recent version, 'ePC-SAFT advanced'. ${ }^{[43]}$ In that work the Born term was altered by a concentration-dependent dielectric constant. This allowed for a more predictive treatment of binary data (activity and osmotic coefficients) and for calculation of transfer properties (Gibbs free energies of solvation and transfer) ${ }^{[43]}$.

In the present work, 'ePC-SAFT advanced' which includes an altered Born term ${ }^{[43]}$ was further investigated for predicting LLE of mixed solvents containing dissolved electrolytes. Furthermore, the role of an appropriate expression of the dielectric constant as function of the system composition was studied with respect to the impact on the quality of the predicted LLE. As pointed out by different authors ${ }^{[44-48]}$ the dielectric constant is among the most important factors for modeling multicomponent electrolyte systems at high salt concentration. Different attempts to model the dielectric constant of pure fluids and fluid mixtures have been made since the beginning of the last century. The most notable of them is the framework developed by Onsager ${ }^{[49]}$ and Kirkwood ${ }^{[50]}$ for predicting the relative dielectric constant of a compound from its molecular properties, and its extension to fluid mixtures ${ }^{[51,52]}$ and electrolytic systems. ${ }^{[53]}$ It is important to note that the aim of this work is to show the potential of modeling improvement by accounting for the change of the dielectric properties of the
medium under different conditions rather than to provide a phenomenological description of the dielectric constant of complex systems.

## Thermodynamic framework

Liquid-liquid equilibria. According to the necessary and sufficient criterion for thermodynamic equilibrium, a system at given pressure and temperature is stable only if the Gibbs energy of this system reaches its global minimum. This state might be reached by splitting into two phases, which results from the complex interplay of all intermolecular interactions as function of temperature, pressure, type of components, and type and concentration of additives such as salts.

In this work, phase equilibria calculation was accomplished by exploiting the isofugacity criterion (Eq. (1)), which states that at thermodynamic equilibrium the fugacity of each component $i=1, \ldots, N$ is equal in all the phases $j=1, . . \pi$ :

$$
\begin{equation*}
f_{i}^{1}=f_{i}^{2}=\ldots=f_{i}^{\pi} \tag{1}
\end{equation*}
$$

In order to avoid convergence to the trivial solution, minimization of the Gibbs energy was performed before imposing the isofugacity criterion. For the calculation of the fugacities, the so-called " $\varphi-\varphi$ " criterion was used. The fugacity of each component in both phases was expressed using the respective fugacity coefficient $\varphi$ according to Eq. (2).

$$
\begin{equation*}
\varphi_{i}^{L 1} \mathrm{x}_{i}^{L 1}=\varphi_{i}^{L 2} \mathrm{x}_{i}^{L 2} \tag{2}
\end{equation*}
$$

For dissociated salts, the fugacity coefficient was calculated from the fugacity coefficient of the anion and the cation according to Eq. (3).

$$
\begin{equation*}
\varphi_{ \pm}=\left(\varphi_{+}^{\nu_{+}} \varphi_{-}^{\nu_{-}}\right)^{\frac{1}{\nu_{+}+\nu_{-}}} \tag{3}
\end{equation*}
$$

One possibility to obtain fugacity coefficients is to derive them by Helmholtz-energy models (i.e., equations of state) using Eq. (4).

$$
\begin{equation*}
\ln \left(\varphi_{\mathrm{i}}\right)=\frac{\mu_{\mathrm{i}}^{\text {res }}}{\mathrm{k}_{\mathrm{B}} \cdot \mathrm{~T}}-\ln \left(1+\left(\frac{\partial\left(\frac{\mathrm{a}^{\text {res }}}{\mathrm{k}_{\mathrm{B}} \cdot \mathrm{~T}}\right)}{\partial \rho}\right)\right) \tag{4}
\end{equation*}
$$

It can be seen that residual properties are required for the modeling. Among such models that allow calculating the missing derivatives in Eq. (4) is ePC-SAFT. The most recent version is called 'ePC-SAFT advanced', and it contains the Born term altered by a concentration-dependent dielectric constant. ${ }^{[43]}$ The theory is briefly discussed in the following.
ePC-SAFT advanced. Among the different classes of thermodynamic models, SAFT-based models represent a family of equations of state (EoS) developed by Chapman in 1989. ${ }^{[54]}$ SAFT-based EoS are known to be able to represent the thermodynamic behavior of complex mixtures, such as those formed by chain-like, polar or associating molecules. PC-SAFT
was first published by Gross and Sadowski in $2001{ }^{[55]}$ and was developed by applying the second-order perturbation theory of Barker and Henderson ${ }^{[56,57]}$ directly to a hard-chain reference fluid. ${ }^{[58]}$ Thus, a more realistic dispersion potential has been developed, which implicitly takes into account the nonspherical shape of the molecules. Later, PC-SAFT was extended to include electrostatic interionic interactions ${ }^{[59]}$ and to account for the change of the dielectric constant at different concentrations. ${ }^{[34]}$ Recently, an altered Born term was introduced within ePC-SAFT, ${ }^{[43]}$ which accounts for electrostatic interactions between charged components and their surrounding medium. The resulting residual Helmholtz energy is given as the sum of all contributions according to Eq. (5).

$$
\begin{equation*}
a^{\text {res }}=a^{h c}+a^{\text {disp }}+a^{\text {assoc }}+a^{D H}+a^{\text {Born }} \tag{5}
\end{equation*}
$$

'ePC-SAFT advanced' falls back to 'ePC-SAFT revised' if $a^{\text {Born }}=0$ and if the dielectric constant is not treated as a function of the salt concentration. 'ePC-SAFT advanced' requires five pure-component parameters for each associating component $i$ : the segment diameter $\sigma_{i}$, the number of segments $m_{i}^{\text {seg }}$, the dispersion-energy parameter $u_{i}$, the association-energy parameter $\varepsilon^{\text {AiBi }}$ and association-volume parameter $\kappa^{\text {AiBi }}$. Furthermore, a binary interaction parameter $k_{i j}$ is required for each pair of components $i$ and $j$ to correct for the deviation of the dispersion energy from the chosen mixing rule (Eq. (7)), as well as the dielectric constant of the pure component at the given temperature $T$ and pressure $p$. Two more interaction parameters $I_{i j}$ and $k_{i j}^{h b}$ can also be introduced (Eq. (6) and (8)), for instance, to accurately describe phase equilibria of mixtures of highly asymmetric components. ${ }^{[42]}$ Mixture properties are calculated according to the combining rules of Berthelot-Lorentz (Eq. (68)).

$$
\begin{align*}
& \sigma_{i j}=\frac{1}{2}\left(\sigma_{i}+\sigma_{j}\right)\left(1-l_{i j}\right)  \tag{6}\\
& u_{i j}=\sqrt{u_{i} u_{j}}\left(1-k_{i j}(T)\right)  \tag{7}\\
& \varepsilon^{A_{i} B_{j}}=\frac{\varepsilon^{A_{i} B_{i}}+\varepsilon^{A_{j} B_{j}}}{2}\left(1-k_{i j}^{h b}\right) \tag{8}
\end{align*}
$$

The binary interaction parameter $k_{i j}(T)$ is described as a linear function of the temperature according to Eq. (9).

$$
\begin{equation*}
k_{i j}(T)=k_{i j, 298.15 k}+k_{i j, T}(T / K-298.15 K) \tag{9}
\end{equation*}
$$

Each of the contributions in Eq. (5) requires one or more pure-component parameters: the hard-chain and dispersion terms are function of segment diameter $\sigma_{i}$, segment number $m_{i}$ and dispersion-energy parameter $u_{i}$ as well as their binary interaction parameters (Eq. 6-7), the association term is function of segment diameter $\sigma_{i}$, association-energy $\varepsilon^{\text {AiBi }}$ (and the binary interaction parameter given by Eq. (8)) and associa-tion-volume parameter $\kappa^{\text {AiBi }}$ and the Debye-Hückel and Born terms are described by the segment diameter $\sigma_{i}$, ion valence $z_{i}$ and the relative dielectric constant $\varepsilon_{r}$. The parameters in Eq. (6-
9) were inherited from literature for the pairs water/organic solvent and anion/cation. No such binary parameters were applied for the pair solvent/ion since electrostatic ion-solvent interaction are now explicitly described by the Born term as explained in the next section.

The Born term and the concentration-dependent dielectric constant. The Born Term is a contribution to the residual Helmholtz energy that accounts for the electrostatic interactions of ionic compounds with their surrounding medium. Differently from neutral species, ionic compounds interact via ion-ion interactions with other ions but also through strong dipole-ion interactions with all the other compounds in the system. The latter interactions are more pronounced when ions are immersed in a highly dipolar solvent, such as water, rather than in a non-polar solvent such as hexane. The presence of ions causes water molecules to reorient around the ions and to shield their electrostatic field. Those dipole-ion interactions can be explicitly accounted for, using for instance the non-primitive MSA, ${ }^{[60-62]}$ or implicitly using the other physical terms of the equation of state. The latter strategy was adopted in the previous version of ePC-SAFT where the dispersion energy between ions and water were adjusted using binary interaction parameters (which are set to zero in the present work). The Born term accounts for those interaction by treating the surrounding medium as a continuum, where the relative dielectric constant of the system is the macroscopic property which accounts for the dipolar nature of the medium. The Born contribution to the residual Helmholtz energy can be derived as the work required to discharge the ions in the vacuum (with dielectric constant $\varepsilon=\varepsilon_{0}$ ) and recharging them in the system ( $\varepsilon=\varepsilon_{0} \varepsilon_{r}$ ). Each ion is treated as a charged sphere characterized by its diameter $\sigma_{i}$ and valence $z_{i} \cdot{ }^{[63]}$ The final expression $a^{\text {Born }}$ is given by Eq. (10):

$$
\begin{equation*}
a^{\text {Born }}=-\frac{e^{2}}{4 \pi \varepsilon_{0} k_{B} T}\left(1-\frac{1}{\varepsilon_{r}}\right) \sum_{i} \frac{x_{i} z_{i}^{2}}{a_{i}} \tag{10}
\end{equation*}
$$

In Eq. (10), $e$ represents the elementary charge. Within this work, the segment diameter $\sigma_{i}$ is used as the pure-component parameter in the hard-chain $a^{\text {hc }}$, Debye-Hückel $a^{D H}$ and Born $a^{\text {Born }}$ terms. Since the dipolar character of the medium is captured by its relative dielectric constant $\varepsilon_{r}$, an expression for the concentration-dependence of the relative dielectric constant on the concentration in Eq. (10) is crucial for the correct representation of the electrostatic contribution to the solvation energy by means of the Born term. Maribo-Mogensen et al. ${ }^{[64]}$ highlighted the importance of the static permittivity, and they proposed a physical model for its calculation of a solvent mixture, ${ }^{[47]}$ even in electrolyte solvent mixtures. ${ }^{[46]}$ Vincze et al. ${ }^{[45]}$ and Valiskó et al. ${ }^{[44]}$ showed the importance of the concentration dependence of the dielectric constant, coupled with the Born term and a term for ion-ion interactions, which explains the non-monotonic behavior of mean ionic activity coefficients (MIACs) of alkali halides in water. In the present work, the treatment of the dielectric constant was extended to multicomponent mixtures by introducing a mass-fraction based
mixing rule. More specifically, four strategies for phase equilibria calculations in LLTPS were compared:

- Strategy 1: without Born term including the concentrationdependent dielectric constant, which is treated by a molefraction based mixing rule (theory from ref. ${ }^{[34]}$ );
- Strategy 2: with altered Born term including the concen-tration-dependent dielectric constant, which is treated by a mole-fraction based mixing rule (theory from ref. ${ }^{[43]}$ );
- Strategy 3: with altered Born term including the concen-tration-dependent dielectric constant, which is treated by a mass fraction-based mixing rule (theory from this work);
- Strategy 4: with altered Born term including the concen-tration-dependent dielectric constant, which is treated by a combined mixing rule, where the dielectric constant of the solvent is calculated according to the mass-fraction mixing rule and the influence of dissolved salts according to the mole-fraction based mixing rule (theory from this work).
The three employed mixing rules are depicted respectively in Eq. (11-13).

$$
\begin{align*}
& \varepsilon_{r}=\sum_{j=1}^{N_{c}} \varepsilon_{r, j} x_{j}  \tag{11}\\
& \varepsilon_{r}=\sum_{j=1}^{N_{c}} \varepsilon_{r, j} w_{j}  \tag{12}\\
& \varepsilon_{r}=\left(\sum_{j=1}^{N_{c}^{s o l}} \varepsilon_{r, j} w_{j}^{s o l}\right) x_{s o l}+\sum_{j=1}^{N_{c}^{\prime}} \varepsilon_{r, j} x_{j} \tag{13}
\end{align*}
$$

In Eq. (11-13) $x_{j}$ and $w_{j}$ denote the mole fraction and mass fraction of the generic component j, respectively. $w_{j}^{\text {sol }}$ in Eq. (13) is the mass fraction of solvent component $j$ in the salt-free solvent, $x_{s o l}$ is the (overall) solvent mole fraction in Eq. (13). $N_{c}$, $N_{c}^{s o l}, N_{c}^{l}$ are the total number of components, the total number of components in the salt-free solvent and the total number of ionic components, respectively. The mass-fraction based mixing rule (Eq. (12)) is expected to give better results than the molefraction based rule (Eq.(11)) in describing the dielectric constant of multicomponent mixtures with components of very different molecular masses, which is also confirmed by experiments ${ }^{[65,66]}$ (see Figure 1).

The reason of using the third mixing rule (Eq. (13)) is the success of the mole-fraction based mixing rule to account for the salt-concentration dependence of the dielectric constant; this was observed in previous work. ${ }^{[34,43]}$ The implementation of a concentration-dependent dielectric constant requires calculating the partial derivative of the respective expression and inserting them in the Debye-Hückel and Born terms. The whole procedure is explained in Bülow et $a l^{[34]}$. Since the Born term is employed, electrostatic interactions of ions with their surrounding medium are considered explicitly. All the binary interaction parameters between ions and any solvent were set to zero, i.e. 'ePC-SAFT advanced' was used as a fully predictive tool. It was found that previously determined ion parameters using data

![](https://cdn.mathpix.com/cropped/34702bfc-4887-4319-921d-54948fd992c1-04.jpg?height=507&width=628&top_left_y=226&top_left_x=1166)
Figure 1. Relative dielectric constant of the mixture water-poly (ethylene glycol)400 as function of the water mole fraction (triangles) and of the water weight fraction (stars). Experimental data from Mali et al. ${ }^{[65]}$. The dashed line represents the calculation of the relative dielectric constant using the weight-fraction based mixing rule (Eq. (12)).

from aqueous solutions only ${ }^{[67,68]}$ gave excellent results when used with the Born term in non-aqueous media. Therefore, the literature pure-ion ePC-SAFT parameters from the version 'ePCSAFT revised' from $2014{ }^{[67]}$ were used in this work without reparametrization.

## Results and Discussion

Model parameters used in this work. Within this work, the LLE of different LLTPS containing salts was modeled using singleion parameters of the alkali-halide salts fitted exclusively to aqueous solutions as published in 'ePC-SAFT revised' from 2014, ${ }^{[67]}$ whereas single-ion parameters of the ionic liquid were fitted to pure-component density data in previous works. ${ }^{[69,70]}$ Water and the organic solvents were modeled as associating fluids with a "2B association scheme". According to this association scheme, the components are assigned two association sites, which model anisotropic short-range interactions (such as hydrogen bonding). The pure-solvent parameters are listed in Table 1, the single-ion parameters and all the binary

Table 1. Pure-component parameters of water and organic solvents used in this work. All components have a 2B association scheme.
| Organic Solvent | $m_{i}^{\text {seg }}$ | $\sigma_{i} / \AA$ | $u_{i} / k_{B} / \mathrm{K}$ | $\varepsilon^{A i B i} / \mathrm{K}$ | $\kappa^{A i B i}$ | Ref. |
| :--- | :--- | :--- | :--- | :--- | :--- | :--- |
| Water | 1.2047 | * | 353.95 | 2425.7 | 0.04509 | [40] |
| 1-Butanol | 2.751 | 3.6139 | 259.59 | 2544.56 | 0.0067 | [31] |
| 1-Pentanol | 3.6260 | 3.4508 | 247.28 | 2252.1 | 0.010319 | [71] |
| MEK | 3.0748 | 3.3932 | 252.27 | 0 | 0.01 | [72] |
| MIBK | 3.3628 | 3.6799 | 259.89 | 0 | 0.01 | [73] |
| PPG | $\mathrm{M}_{\text {PPG }}$. 0.0363 | 3.35 | 190.49 | 1749.0 | 0.0298 | [74] |
| ${ }^{*} \sigma=2.7927+\left(10.11 \cdot \mathrm{e}^{-0.01775 \mathrm{~T}}-1.417 \cdot \mathrm{e}^{-0.01146 \mathrm{~T}}\right)$ |  |  |  |  |  |  |


interaction parameters are listed, respectively, in Tables 2-4. The success of the transferability of the single-ion parameters to model MIACs and solubility in pure alcohol solvents was shown in previous work. The focus of the present work lies on solvent mixtures based on water + organic solvent, which additionally contain dissolved salt or ionic liquid. Aqueous solution properties (solution density, osmotic and activity coefficients) were used to adjust single-ion parameters in 'ePCSAFT revised' from 2014, ${ }^{[67]}$ which is a model without the Born term. In order to fit experimental data up to high ion concentrations, binary interaction parameters had to be included in the 'ePC-SAFT revised' model. These parameters were inherited in the present work. In contrast, interaction parameters between ions and water were not inherited from 'ePC-SAFT revised', but all ion-solvent interaction parameters were set to zero in the present work.

The relative dielectric constants of the pure components investigated in this work, are listed in Table 5. These were used to describe the mean relative dielectric constant of the mixtures under study according to the used mixing rules in Eq. (11-13).

Experimental data. As model systems to validate the new methodology, five binary systems water + organic solvent and their combination with different salts were chosen. Five solvents used as extraction agents in the chemical and pharmaceutical industry were chosen: the ketones methyl ethyl ketone (MEK)

| Table 2. ePC-SAFT pure-ion parameters from ref. ${ }^{[67]}$ Segment number is unity for all ions. $k_{i j}$ between ions and water is set to zero. |  |  |  |  |
| :--- | :--- | :--- | :--- | :--- |
| Ion | $m_{i}^{\text {seg }} /-$ | $\sigma_{i} / \AA$ | $u_{i} / k_{B} / \mathrm{K}$ | Ref. |
| $\mathrm{Na}^{+[\mathrm{a}]}$ | 1 | 2.8232 | 230.00 | [67] |
| $\mathrm{K}^{+[\mathrm{a}]}$ | 1 | 3.3417 | 200.00 | [67] |
| $\mathrm{Cl}^{-[\mathrm{a}]}$ | 1 | 2.7560 | 170.00 | [67] |
| $\mathrm{Br}^{-[\mathrm{a}]}$ | 1 | 3.0707 | 190.00 | [67] |
| $\mathrm{Ac}^{-[\mathrm{b}]}$ | 3.7266 | 3.5605 | 533.1138 | [70] |
| $\mathrm{C}_{4} \mathrm{mim}^{+[\mathrm{b}]}$ | 2.4805 | 3.6371 | 218.1441 | [69] |
| [a] ion-ion dispersion allowed only between cation-anion; diameter $\mathrm{d}=$ \sigma is independent of temperature; [b] ion-ion dispersion allowed; diameter $\mathrm{d}(\mathrm{T})$ is temperature-dependent. |  |  |  |  |


| Table 3. ePC-SAFT binary interaction parameters $k_{i j}$ between <br> halide anions and alkali cations. [67] |  |  |
| :--- | :--- | :--- |
| $k_{\text {cat, an }}$ | $\mathrm{Na}^{+}$ | $\mathrm{K}^{+}$ |
| $\mathrm{Cl}^{-}$ | 0.317 | 0.064 |
| $\mathrm{Br}^{-}$ | 0.290 | -0.102 |


| Table 5. Relative dielectric constants for solvents and salts used in this work. T in Kelvin. |  |  |
| :--- | :--- | :--- |
| Component | Relative dielectric constant | Ref. |
| Water | $-105.2 \ln T+677.480$ | [75] |
| 1-Butanol | $-0.1077 T+49.723$ | This work ${ }^{[76-78]}$ |
| 1-Pentanol | $73.397-0.28165 T+2.8427-10^{-4} \mathrm{~T}^{2}$ | [79] |
| MEK | $18.389-0.1025(T-298.15)$ | This work ${ }^{[80]}$ |
| MIBK | $36.341-0.09712 T+6.1896-10^{-5} T^{2}$ | [79] |
| PPG (298.15 K) | $92.483\left(M_{w}\right)^{-0.347}$ | This work ${ }^{[81]}$ |
| Alkali halides | 8 | a) |
| [ $\mathrm{C}_{4}$ mim] [Ac] | 11 | b) |
| a) All alkali halides were modeled with a relative dielectric constant that is a mean of available experimental data. ${ }^{[82]}$ <br> b) For the IL , a relative dielectric constant of 11 was chosen according to a previous work. ${ }^{[34]}$ |  |  |


| Table 6. Investigated LLE data (at 298.15 K and 1 bar) with the number of data points (NP) used in this work to test the different modeling strategies. |  |  |
| :--- | :--- | :--- |
| System | NP | Ref |
| Water + 1-butanol + KBr | 10 | [86] |
| Water + 1-butanol + NaCl | 10 | [86] |
| Water + 1-pentanol + KCl | 10 | [87] |
| Water $+\mathrm{MEK}+\mathrm{KBr}$ | 10 | [86] |
| Water + MEK + KCl | 10 | [86] |
| Water + MEK + NaCl | 10 | [86] |
| Water $+\mathrm{MIBK}+\mathrm{NaCl}$ | 10 | [32] |
| Water + PPG725 + NaCl | 8 | [88] |
| Water + PPG400 + [ $\mathrm{C}_{4}$ mim][Ac] | 10 | [89] |

and methyl isobutyl ketone (MIBK), the aliphatic alcohols 1butanol and 1-pentanol and the polymer PPG. Table 6 provides an overview of the studied systems, the temperature, and the reference of the used experimental data. In order to obtain good predictions in the phase equilibria calculation of multicomponent systems, it is important to accurately model each of the binary sub-system, as already demonstrated in previous works. ${ }^{[29,30,83-85]}$ Binary interaction parameters between organic solvent and water are available for most of the studied water + organic solvent mixtures, see Table 4. However, binary parameters for water + MEK and water + PPG were missing in the literature and were fitted in this work to LLE data of binary mixtures with water.

| Table 4. ePC-SAFT binary interaction parameters $k_{i j, 298.15 K}, k_{i j, T}, l_{i j}$ and $k_{i j}^{h b}$ between water and the investigated organic solvents. |  |  |  |  |  |
| :--- | :--- | :--- | :--- | :--- | :--- |
|  | 1-Butanol | 1-Pentanol | MEK | MIBK | PPG |
| $k_{i j, 298.15 k}$ | -0.102 | 0.016 | -0.2719 | -0.055 | -0.167 |
| $k_{i j, T}$ | $2.94 \times 10^{-4}$ | 0 | $4.27 \times 10^{-4}$ | 0 | $6.67 \times 10^{-4}$ |
| $l_{i j}$ | -0.0044 | 0 | 0.0579 | 0 | 0 |
| $k_{i j}^{h b}$ | 0.026 | 0 | 0 | 0.097 | 0 |
| Ref. | [31] | [71] | This work | [32] | This work |

Effect of the Born term and of different treatments of the dielectric constant in the mixture on the accuracy of the predicted LLE of LLTPS. Table 7 gives an overview of the performance of the different modeling strategies employed in this work. The modeling results are evaluated by two measures, AAD (Eq. (14)) and ARD (Eq. (15)) according to:

$$
\begin{align*}
& A R D=\frac{1}{N P \cdot N_{C}} \sum_{k=1}^{N P} \sum_{i=1}^{N_{C}}\left|1-\frac{w_{i}^{\text {calc }}}{w_{i}^{\text {exp }}}\right| \cdot 100 \%  \tag{14}\\
& A A D=\frac{1}{N P \cdot N_{C}} \sum_{k=1}^{N P} \sum_{i=1}^{N_{C}}\left|w_{i}^{\text {calc }}-w_{i}^{\text {exp }}\right| \tag{15}
\end{align*}
$$

In Eq. (14-15) both sums are made over the number of data points ( $N P$ ) and the number of components ( $N_{c}$ ) for each data point. It is important to state that the AAD and ARD refers to all data points in both phases for all components, respectively. As it can be seen, strategy 1 (without the Born term) predicts that four of the nine studied system are LLTPS, while the other five systems were predicted to be homogeneous mixture at the experimental conditions. Including the Born term with a molefraction based concentration dependence of the dielectric constant (strategy 2) partially improves the predictive capability. However, for three systems out of the nine systems under study strategy 2 did not predict phase separation, while still high AAD and ARD values remain. A remarkable improvement can be observed by combining the Born term with either the massfraction based mixing rule (strategy 3) or the combined mixing rule (strategy 4) to described the dielectric constant of the mixture. Strategies 3 and 4 allow reducing the ARD and AAD values, in average, by one order of magnitude compared to strategy 2. That is, not only the Born term itself is important for the successful modeling of electrolyte systems, also taking into account the dependence of the dielectric constant on salt concentration is a very sensitive property that has to be included consistently in any electrolyte model. This answers partially also the main questions from recent literature. ${ }^{[90,91]}$

Figure 2 compares experimental LLE data and modeling results for the system water $+\mathrm{MEK}+\mathrm{NaCl}$ and provides insights
into the physical meaning of each used strategies and mixing rules within the ePC-SAFT framework. The results shown are representative for all the investigated LLTPS water+organic solvent + salt. The binary system water+MEK shows partial miscibility at 298.15 K and 1 bar, with the organic phase containing a water mass-fraction of 0.12 and the aqueous phase containing a MEK mass-fraction of 0.25 . Addition of NaCl leads to a strong salting-out of MEK from the aqueous phase and of water from the organic phase. At the same time NaCl prevails in the aqueous phase and is almost absent in the organic phase. In the following, the results using the different modeling strategies 1-4 are discussed.

Strategy 1 without using the Born term (Figure 2a) predicted phase diagrams that are fundamentally wrong. Instead of an enlargement of the miscibility gap (salting-out), a salting-in effect is predicted manifested in an increasing mutual solubility between MEK and water with increasing salt concentration. Furthermore, the partition of salt is predicted qualitatively incorrect, i.e. strategy 1 predicts that salt partitions preferentially to the organic phase. This is in contradiction to the experimental data. This misprediction is due to incomplete theory of 'ePC-SAFT revised'. On the one hand, according to the Debye-Hückel theory, decreasing the dielectric constant leads to a weaker screening of the ions electrostatic fields and therefore to stronger ion-ion interactions. On the other hand, according to the Born term increasing the dielectric constant leads to a stronger electrostatic ion solvation, and the impact of this effect on the residual Helmholtz energy is in general much greater than non-electrostatic or ion-ion electrostatic contributions ${ }^{[92]}$ (at least for ions with high charge density). Thus, neglecting explicitly electrostatic ion-solvent interactions within a modeling framework is the cause of the wrong prediction observed in LLE calculations with dissolved salts using the previous 'ePC-SAFT revised' version. ${ }^{[32,33,93]}$ The only way to compensate for the misprediction is to adjust the dispersion energy between ion-solvent, which resulted in very high absolute numbers of the binary interaction parameters; alternatively, introducing solvent-dependent pure-ion parameters ${ }^{[42]}$ was applied, which is not the best (and also not a pragmatic) solution.

Table 7. Comparison of the predictive modeling results for the different strategies used in this work, evaluated as ARD and AAD (overall composition within both phases) to predict LLE of different LLTPS. No value means that the strategy within the ePC-SAFT framework did not predict a phase split into two liquid phases at the experimental conditions.
| System | Strategy 1 |  | Strategy 2 |  | Strategy 3 |  | Strategy 4 |  |
| :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- |
|  | ARD/\% | AAD/- | ARD/\% | AAD/- | ARD/\% | AAD/- | ARD/\% | AAD/- |
| Water + 1-butanol + KBr | 306 | 0.0648 | 196 | 0.0554 | 15 | 0.0073 | 12 | 0.0066 |
| Water + 1-butanol + NaCl | 1021 | 0.0974 | 494 | 0.1270 | 23 | 0.0089 | 19 | 0.0089 |
| Water + 1-pentanol + KCl | 3356 | 0.0665 | 455 | 0.0162 | 28 | 0.0088 | 29 | 0.0086 |
| Water $+\mathrm{MEK}+\mathrm{KBr}$ | - | - | - | - | 15 | 0.0074 | 9 | 0.0078 |
| Water + MEK + KCl | - | - | 1604 | 0.0897 | 14 | 0.0137 | 10 | 0.0097 |
| Water + MEK + NaCl | - | - | 269 | 0.0626 | 20 | 0.0073 | 15 | 0.0043 |
| Water + MIBK+ NaCl | 3764 | 0.0774 | 150 | 0.0054 | 28 | 0.0013 | 21 | 0.0010 |
| Water + PPG725 + NaCl | - | - | - | - | 74 | 0.1216 | 69 | 0.0710 |
| Water + PPG400 + [ $\mathrm{C}_{4}$ mim] [Ac] | - | - | - | - | 203 | 0.3761 | - | - |


![](https://cdn.mathpix.com/cropped/34702bfc-4887-4319-921d-54948fd992c1-07.jpg?height=1008&width=1375&top_left_y=209&top_left_x=346)
Figure 2. Comparison of the different modeling strategies used to predict LLE of the LLTPS water $+\mathrm{MEK}+\mathrm{NaCl}$ at $\mathrm{T}=298.15 \mathrm{~K}$ and 1 bar. Grey circles and solid lines are experimental tie-lines (refs. in Table 6), white circles with dashed lines are ePC-SAFT predicted tie-lines: a.) Strategy 1, b.) Strategy 2, c.) Strategy 3, d.) Strategy 4. All strategies were applied using the same parameters listed in Tables 1-4.

In other words, in all recent ePC-SAFT works, strong ionsolvent interactions were corrected by fitting binary interaction parameters to experimental LLE data. Certainly, this made ePCSAFT a strong correlative tool with high accuracy and the ability to extrapolation (to a certain extent). However, reliable predictive capability could not be achieved so far by ePC-SAFT modeling.

By introducing the Born term with a mole-fraction based mixing rule for the dielectric constant in mixtures (Figure 2b), LLE modeling could be improved, but the results are still far away from being in agreement with the experimental LLE data. Although the tie-line slope (i.e. the partitioning of NaCl between both phases) is predicted qualitatively correct, still the salting-out effect could not be predicted.

Qualitative and quantitative agreement of LLE prediction with experiments could only be achieved using strategy 3 and 4 (i.e. using the Born term with a mass-fraction based rule or a combined mixing rule for the dielectric constant), cf. Figure 2c and Figure 2d. These results emphasize the importance of using an accurate value for the relative dielectric constant for modeling complex electrolyte systems combined with the altered Born term. The weight-based mixing rule (in Eq. (12)) used to calculate the relative dielectric constant in strategy 3 represents a much more physically sound approach than using, for example, the pure-component relative dielectric constant of one of the solvent components or a mole-fraction based mixing rule (see Figure 1). Below a certain salt concentration, a mole-
fraction based mixing rule between salt and solvent ${ }^{[34,43]}$ allowed describing solvent-salt interactions meaningfully. Thus, a combined approach for taking into account the dependence of the dielectric constant of a mixture salt-solvent was also confirmed; that is, ARD and AAD values slightly decreased when strategy 4 is used instead of strategy 3 . The only exception remains the mixture water + PEG400 $+\left[\mathrm{C}_{4} \mathrm{mim}\right][\mathrm{Ac}]$. Only strategy 3 allowed predicting the phase equilibrium in this ATPS. It cannot be excluded that the mole-fraction based mixing rule fails to capture the real dependency of the relative dielectric constant on concentration, caused by a complex interaction interplay between the IL and the solvent components. Further, the used mixing rule in strategy 4 (Eq. (13)) is of empirical nature.

Figure 3 presents a comparison between experimental and calculated partition coefficients $K$ of salt, expressed as massfraction ratio (i.e., $w_{\text {salt }}$ in the aqueous phase related to $w_{\text {salt }}$ in the organic phase), for four of the investigated systems. As the figures show, the correct order of magnitude of the salt partition coefficient could be predicted only using strategies 3 and 4 , whereas the partition coefficients are heavily underestimated using the strategies 1 and 2. Partition coefficients of the salts KCl and KBr in water-MEK (Figure 3c.) and d.) could also be predicted, using strategies 3 and 4, in very close agreement to the experimental data. Results show that without the Born term (strategy 1) $K$ values are below one, i.e. strategy 1 predicts that the salt is mainly present in the organic phase.

![](https://cdn.mathpix.com/cropped/34702bfc-4887-4319-921d-54948fd992c1-08.jpg?height=1291&width=1511&top_left_y=202&top_left_x=277)
Figure 3. Partition coefficients of a salt in four different LLTPS at 298.15 K and 1 bar. a.) water +1 -butanol +NaCl , b.) water $+\mathrm{MEK}+\mathrm{NaCl}$, c.) water + MEK + KCI, d.) water + MEK + KBr. Symbols are partition coefficients calculated from experimental data (full circles), predictions using the ePC-SAFT framework and strategy 1 (triangles), strategy 2 (diamonds), strategy 3 (open circles) or strategy 4 (stars). Lines are intended to guide the readers eyes.

The conclusion from this is that the Born term is crucially required to obtain the correct partitioning of a salt between an aqueous phase and an organic phase without fitting unreasonable model parameters.

## Conclusion

Within this work 'ePC-SAFT advanced' was used to predict LLE of mixed-solvent systems containing dissolved salts, using different mixing rules for the dielectric constant into the Born and the Debye-Hückel contributions to the residual Helmholtz energy. The importance of an altered Born term was demonstrated in two previous works for calculating MIACs and solidliquid equilibria of alkali halides in alcoholic solvents. In the present work, first an appropriate expression of the concen-tration-dependence of the dielectric constant for mixtures containing water + organic solvent + salt was suggested. The vital importance of the dielectric constant was already discussed by Maribo-Mogensen et al., which proposed a new
modeling approach taking into account hydrogen bonding of associating fluids ${ }^{[47]}$ and phenomena like dielectric saturation. ${ }^{[46]}$ However, by using a simple mass-fraction based mixing rule for the dielectric constant quantitative prediction of ternary LLE was possible for almost all the studied systems, without using any binary interaction parameters. The combined mixing rule seems to perform slightly better than the pure mass-fraction based rule, except for the system containing an ionic liquid. Using the previous version 'ePC-SAFT revised from 2014' (without the Born term) did not allow any predictions, and LLE of electrolyte systems could only be correlated. The Born term and a concentration-dependent dielectric constant enhance the predictive capability of 'ePC-SAFT advanced' and improve its correlative capability. This work doesn't aim, however, at underemphasize the importance of rigorous modeling of the dielectric constant of electrolyte systems: because of the empirical nature of the used expression, correctness cannot be guaranteed for all the systems and at all existing conditions, which should be investigated in future works. Furthermore, the performance of the Born term with a concentration-dependent
dielectric constant should be tested for multivalent electrolytes $\left(\mathrm{Mg}^{2+}, \mathrm{Ca}^{2+}\right)$ and for ions which deviate from the spheroidal form $\left(\mathrm{CH} 3 \mathrm{COO}^{-}, \mathrm{SO}_{4}{ }^{2-}, \mathrm{PO}_{4}{ }^{3-}\right)$ assumed within the Born term.

## Acknowledgements

This work is dedicated to Christoph Janiak and his research group. We met within the priority program SPP 1708 "Material Synthesis Near Room Temperature" and from then on started successful and fruitful collaboration. Open access funding enabled and organized by Projekt DEAL.

Keywords: ATPS • ionic liquids • electrolyte thermodynamics • liquid-liquid equilibria • partition coefficient
[1] X. Wu, Y. Gong, S. Xu, Z. Yan, X. Zhang, S. Yang, J. Chem. Eng. Data 2019, 64, 4319-4329.
[2] B. Golembiewski, N. vom Stein, N. Sick, H.-D. Wiemhöfer, J. Cleaner Prod. 2015, 87, 800-810.
[3] J. T. Safarov, Fluid Phase Equilib. 2005, 236, 87-95.
[4] E. O. Eisen, J. Joffe, J. Chem. Eng. Data 1966, 11, 480-484.
[5] Y. Román-Leshkov, C. J. Barrett, Z. Y. Liu, J. A. Dumesic, Nature 2007, 447, 982-985.
[6] A. M. Lopes, V. C. Santos-Ebinuma, A. Pessoa Júnior, C. O. Rangel-Yagui, Braz. J. Chem. Eng. 2014, 31, 1057-1064.
[7] J. Fu, Ind. Eng. Chem. Res. 2004, 43, 1279-1283.
[8] Y. Chen, F. Mutelet, J.-N. Jaubert, J. Phys. Chem. B 2012, 116, 14375-14388.
[9] D. A. Weingaertner, S. Lynn, D. N. Hanson, Ind. Eng. Chem. Res. 1991, 30, 490-501.
[10] F. Hofmeister, Arch. Exp. Pathol. Pharmakol. 1888, 25, 1-30.
[11] K. D. Collins, M. W. Washabaugh, Q. Rev. Biophys. 1985, 18, 323-422.
[12] M.-K. Shahbazinasab, F. Rahimpour, J. Chem. Eng. Data 2012, 57, 1867-1874.
[13] R. Hatti-Kaul, M\&B Lab. Bull. 2001, 19, 269-278.
[14] M. G. Freire, A. F. M. Cláudio, J. M. M. Araújo, J. A. P. Coutinho, I. M. Marrucho, J. N. Canongia Lopes, L. P. N. Rebelo, Chem. Soc. Rev. 2012, 41, 4966-4995.
[15] N. J. Bridges, K. E. Gutowski, R. D. Rogers, Green Chem. 2007, 9, 177-183.
[16] P. G. Mazzola, A. M. Lopes, F. A. Hasmann, A. F. Jozala, T. C. V. Penna, P. O. Magalhaes, C. O. Rangel-Yagui, A. Pessoa Jr, J. Chem. Technol. Biotechnol. 2008, 83, 143-157.
[17] T. J. Peters, Cell Biochem. Funct. 1987, 5, 233-234.
[18] B. J. Zaslavskij, Aqueous two-phase partitioning: Physical chemistry and bioanalytical applications; Dekker: New York, 1995.
[19] H. D. Willauer, J. G. Huddleston, R. D. Rogers, Ind. Eng. Chem. Res. 2002, 41, 2591-2601.
[20] J. G. Huddleston, H. D. Willauer, S. T. Griffin, R. D. Rogers, Ind. Eng. Chem. Res. 1999, 38, 2523-2539.
[21] J. F. B. Pereira, Á. S. Lima, M. G. Freire, J. A. P. Coutinho, Green Chem. 2010, 12, 1661.
[22] Y. Lin, K. Thomsen, J.-C. de Hemptinne, AIChE J. 2007, 53, 9891005.
[23] R. Inchekel, J.-C. de Hemptinne, W. Fürst, Fluid Phase Equilib. 2008, 271, 19-27.
[24] B. Maribo-Mogensen, K. Thomsen, G. M. Kontogeorgis, AIChE J. 2015, 61, 2933-2950.
[25] J. Wu, J. M. Prausnitz, Ind. Eng. Chem. Res. 1998, 37, 1634-1643.
[26] J. A. Myers, S. I. Sandler, R. H. Wood, Ind. Eng. Chem. Res. 2002, 41, 3282-3297.
[27] K. Wysoczanska, H. T. Do, G. Sadowski, E. A. Macedo, C. Held, AIChE J. 2020, 66.
[28] K. Wysoczanska, B. Nierhauve, G. Sadowski, E. A. Macedo, C. Held, Fluid Phase Equilib. 2021, 527, 112830.
[29] T. Reschke, C. Brandenbusch, G. Sadowski, Fluid Phase Equilib. 2014, 368, 91-103.
[30] T. Reschke, C. Brandenbusch, G. Sadowski, Fluid Phase Equilib. 2014, 375, 306-315.
[31] A. Nann, C. Held, G. Sadowski, Ind. Eng. Chem. Res. 2013, 52, 18472-18481.
[32] S. Mohammad, C. Held, E. Altuntepe, T. Köse, T. Gerlach, I. Smirnova, G. Sadowski, Fluid Phase Equilib. 2016, 416, 83-93.
[33] S. Mohammad, C. Held, E. Altuntepe, T. Köse, G. Sadowski, J. Phys. Chem. B 2016, 120, 3797-3808.
[34] M. Bülow, X. Ji, C. Held, Fluid Phase Equilib. 2019, 492, 26-33.
[35] K. Wysoczanska, E. A. Macedo, G. Sadowski, C. Held, Ind. Eng. Chem. Res. 2019, 58, 21761-21771.
[36] M. Voges, I. V. Prikhodko, S. Prill, M. Hübner, G. Sadowski, C. Held, J. Chem. Eng. Data 2017, 62, 52-61.
[37] C. M. S. S. Neves, C. Held, S. Mohammad, M. Schleinitz, J. A. P. Coutinhoa, M. G. Freire, Phys. Chem. Chem. Phys. 2015, 17, 32044-32052.
[38] A. Nann, J. Mündges, C. Held, S. P. Verevkin, G. Sadowski, J. Phys. Chem. B 2013, 117, 3173-3185.
[39] X. Ji, C. Held, G. Sadowski, Fluid Phase Equilib. 2014, 363, 5965.
[40] C. Held, L. F. Cameretti, G. Sadowski, Fluid Phase Equilib. 2008, 270, 87-96.
[41] C. Held, G. Sadowski, Fluid Phase Equilib. 2009, 279, 141-148.
[42] C. Held, A. Prinz, V. Wallmeyer, G. Sadowski, Chem. Eng. Sci. 2012, 68, 328-339.
[43] M. Bülow, M. Ascani, C. Held, Fluid Phase Equilib. 2021, 112967.
[44] M. Valiskó, D. Boda, J. Chem. Phys. 2014, 140, 234508.
[45] J. Vincze, M. Valiskó, D. Boda, J. Chem. Phys. 2010, 133, 154507.
[46] B. Maribo-Mogensen, G. M. Kontogeorgis, K. Thomsen, J. Phys. Chem. B 2013, 117, 10523-10533.
[47] B. Maribo-Mogensen, G. M. Kontogeorgis, K. Thomsen, J. Phys. Chem. B 2013, 117, 3389-3397.
[48] I. Y. Shilov, A. K. Lyashchenko, J. Phys. Chem. B 2015, 119, 10087-10095.
[49] L. Onsager, J. Am. Chem. Soc. 1936, 58, 1486-1493.
[50] J. G. Kirkwood, J. Chem. Phys. 1939, 7, 911-919.
[51] G. Oster, J. Am. Chem. Soc. 1946, 68, 2036-2041.
[52] A. H. Harvey, J. M. Prausnitz, J. Solution Chem. 1987, 16, 857869.
[53] P. Wang, A. Anderko, Fluid Phase Equilib. 2001, 186, 103-122.
[54] W. G. Chapman, K. E. Gubbins, G. Jackson, M. Radosz, Fluid Phase Equilib. 1989, 52, 31-38.
[55] J. Gross, G. Sadowski, Ind. Eng. Chem. Res. 2001, 40, 1244-1260.
[56] J. A. Barker, D. Henderson, J. Chem. Phys. 1967, 47, 4714-4721.
[57] J. A. Barker, D. Henderson, J. Chem. Phys. 1967, 47, 2856-2861.
[58] W. G. Chapman, G. Jackson, K. E. Gubbins, Mol. Phys. 1988, 65, 1057-1079.
[59] L. F. Cameretti, G. Sadowski, J. M. Mollerup, Ind. Eng. Chem. Res. 2005, 44, 3355-3362.
[60] L. Blum, F. Vericat, W. R. Fawcett, J. Chem. Phys. 1992, 96, 3039-3044.
[61] L. Blum, D. Q. Wei, J. Chem. Phys. 1987, 87, 555-565.
[62] D. Wei, L. Blum, J. Chem. Phys. 1987, 87, 2999-3007.
[63] M. Born, Z. Phys. 1920, 1, 45-48.
[64] B. Maribo-Mogensen, G. M. Kontogeorgis, K. Thomsen, Ind. Eng. Chem. Res. 2012, 51, 5353-5363.
[65] C. S. Mali, S. D. Chavan, K. S. Kanse, A. C. Kumbharkhane, S. C. Mehrotra, Indian J. Pure Appl. Phys. 2007 45(5), 476-481.
[66] C. G. Malmberg, A. A. Maryott, J. Res. Natl. Bur. Stand. 1950, 45, 299-303.
[67] C. Held, T. Reschke, S. Mohammad, A. Luza, G. Sadowski, Chem. Eng. Res. Des. 2014, 92, 2884-2897.
[68] H. Niedermeyer, J. P. Hallett, I. J. Villar-Garcia, P. A. Hunt, T. Welton, Chem. Soc. Rev. 2012, 41, 7780-7802.
[69] X. Ji, C. Held, G. Sadowski, Fluid Phase Equilib. 2012, 335, 6473.
[70] X. Ji, C. Held, Fluid Phase Equilib. 2016, 410, 9-22.
[71] J. Gross, G. Sadowski, Ind. Eng. Chem. Res. 2002, 41, 5510-5515.
[72] A. Tihic, G. M. Kontogeorgis, N. von Solms, M. L. Michelsen, Fluid Phase Equilib. 2006, 248, 29-43.
[73] M. Knierbein, M. Voges, C. Held, Org. Process Res. Dev. 2020, 24, 1052-1062.
[74] I. Stoychev, J. Galy, B. Fournel, P. Lacroix-Desmazes, M. Kleiner, G. Sadowski, J. Chem. Eng. Data 2009, 54, 1551-1559.
[75] W. B. Floriano, M. A. C. Nascimento, Braz. J. Phys. 2004, 34, 3841.
[76] R. D. Heyding, C. A. Winkler, Can. J. Chem. 1951, 29, 790-803.
[77] W. Dannhauser, R. H. Cole, J. Chem. Phys. 1955, 23, 1762-1766.
[78] D.-Y. Chu, Q. Zhang, R.-L. Liu, J. Chem. Soc. Faraday Trans. 1 1987, 83, 635.
[79] H. Landolt, K.-H. Hellwege, O. Madelung, Zahlenwerte und Funktionen aus Naturwissenschaften und Technik, Neue Serie; Springer: Berlin, 1991.
[80] A. Ghanadzadeh Gilani, N. Paktinat, M. Moghadam, J. Chem. Thermodyn. 2011, 43, 569-575.
[81] A. V. Sarode, A. C. Kumbharkhane, J. Mol. Liq. 2011, 160, 109113.
[82] C. Andeen, J. Fontanella, D. Schuele, Phys. Rev. B 1970, 2, 5068-5073.
[83] D. Pabsch, C. Held, G. Sadowski, J. Chem. Eng. Data 2020, 65, 5768-5777.
[84] C. Held, G. Sadowski, A. Carneiro, O. Rodríguez, E. A. Macedo, AIChE J. 2013, 59, 4794-4805.
[85] T. Reschke, C. Brandenbusch, G. Sadowski, Fluid Phase Equilib. 2015, 387, 178-189.
[86] Z. Li, Y. Tang, Y. Liu, Y. Li, Fluid Phase Equilib. 1995, 103, 143153.
[87] N. Boluda, V. Gomis, F. Ruiz, M. D. Saquete, N. Barnes, Fluid Phase Equilib. 2001, 179, 269-276.
[88] E. L. Cheluget, S. Gelinas, J. H. Vera, M. E. Weber, J. Chem. Eng. Data 1994, 39, 127-130.
[89] C. Wu, J. Wang, Y. Pei, H. Wang, Z. Li, J. Chem. Eng. Data 2010, 55, 5004-5008.
[90] C. Held, J. Chem. Eng. Data 2020, 65, 5073-5082.
[91] G. M. Kontogeorgis, B. Maribo-Mogensen, K. Thomsen, Fluid Phase Equilib. 2018, 462, 130-152.
[92] S. Ahmed, N. Ferrando, J.-C. de Hemptinne, J.-P. Simonin, O. Bernard, O. Baudouin, Fluid Phase Equilib. 2018, 459, 138-157.
[93] S. Mohammad, G. Grundl, R. Müller, W. Kunz, G. Sadowski, C. Held, Fluid Phase Equilib. 2016, 428, 102-111.

Manuscript received: January 26, 2021
Revised manuscript received: March 24, 2021
Accepted manuscript online: March 26, 2021


[^0]:    [a] M. Ascani, C. Held
    Laboratory of Thermodynamics, Department of Biochemical and Chemical Engineering, TU Dortmund University, Emil-Figge Str. 70, 44277 Dortmund, Germany
    E-mail: christoph.held@tu-dortmund.de
    © 2021 The Authors. Zeitschrift für anorganische und allgemeine Chemie published by Wiley-VCH GmbH. This is an open access article under the terms of the Creative Commons Attribution License, which permits use, distribution and reproduction in any medium, provided the original work is properly cited.

