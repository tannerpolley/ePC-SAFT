# ePC-SAFT advanced - Part I: Physical meaning of including a concentration-dependent dielectric constant in the born term and in the Debye-Hückel theory 

Mark Bülow, Moreno Ascani, Christoph Held*<br>Laboratory of Thermodynamics, Department of Biochemical and Chemical Engineering, TU Dortmund, Emil-Figge Str. 70, 44277 Dortmund, Germany

## ARTICLE INFO

## Article history:

Received 4 November 2020
Revised 28 December 2020
Accepted 12 January 2021
Available online 24 January 2021

## Keywords:

Non-aqueous systems
Mean ionic activity coefficients
Electrolyte thermodynamics
Gibbs energy of transfer
Gibbs energy of solvation, Dielectric constant


#### Abstract

The transition from aqueous electrolyte systems to non-aqueous electrolyte systems is highly demanded in industrial applications and especially challenging for physics-based thermodynamic models. Electrolyte thermodynamics is a complex matter, and still not all physico-chemical effects are accounted for in state-of-the-art equations of state. The dielectric constant of non-aqueous electrolyte systems changes drastically compared to aqueous systems. One main consequence is that ions are very differently solvated in non-aqueous medium compared to aqueous medium. The Born term represents a methodology to account for the influence of solvation energies of ions, which is based on influences of solvent and salt on the dielectric constant. Utilizing the Born term in electrolyte models is extensively debated, and it is often reasonably neglected in predominantly aqueous systems. Yet, it has a significant influence on transferability from aqueous to non-aqueous media i.e., systems with a large difference in polarity or permittivity compared to aqueous systems. In this work, a modified Born term was combined with electrolyte Perturbed-Chain Statistical Associating Fluid Theory (ePC-SAFT) by introducing additionally a salt concentration-dependent dielectric constant, henceforth called altered Born contribution. The new methodology was validated against infinite dilution properties for ion-solvent interactions: Gibbs energy of hydration and Gibbs energy of transfer of alkali halides from water to alcoholic solvents. Further, mean ionic activity coefficients (MIACs) of alkali halides in alcoholic solvents were quantitatively correct predicted with the advanced ePC-SAFT approach. Original ePC-SAFT parameters were applied for all predictions, and no further binary parameters were adjusted. Based on the success of the model predictions, the transferability of pure-ion ePC-SAFT parameters to organic solvents was verified and the incorporation of concentration-dependent dielectric constant into the altered Born contribution and Debye-Hückel theory was proven to be meaningful methods for the transfer of electrolyte thermodynamic models from aqueous to non-aqueous systems.


© 2021 Elsevier B.V. All rights reserved.

## 1. Introduction

The significance of non-aqueous electrolyte systems and pureionic systems in industrial applications is intensifying. Usually, unit operations involving electrolytes are perfectly optimized for applications with water as the main solvent; however, intensifying upstream and downstream processes involving organic solvents or even ionic fluids are emerging. Non-aqueous solvents may enhance the efficiency of reactions and purifications to reduce energy demand and costs. Organic solvents in classical non-electrolyte thermodynamic models are handled frequently without large expenditure. Processes with electrolytes especially involve phase changes,

[^0]solubility issues and other important properties, which govern the design of unit operations. The design of such unit operations therefore is built on operating experience with a lack in general understanding of intrinsic and bulk interactions rather than on true knowledge.

Facing these challenges, non-aqueous systems are more and more in focus of methods from electrolyte thermodynamics. The transfer from aqueous to pure-organic systems requires redesigning state-of-the-art modeling theories. This is mainly caused by model parameters, which have usually been designed especially for aqueous systems. Two arising major aspects that need to be addressed are the change in polarity - evident in a drastic decrease in static dielectric constant of the solvent - and direct interactions of the continuum solvent with the electrolytes that is the change in solvation energies. Both aspects are related and
must be addressed in a combined extension to the classical theory derived for aqueous systems. Extending an already established model comes with the benefit of remaining reliability for aqueous systems while accessing the pure-organic and pure-ionic systems in one holistic methodology. The model must furthermore be reducible to the aqueous application without changing any parameters. Moreover, it allows the model to account for additional theories in the future, covering so far unemployed physics-based interactions. Kontogeorgis et al. [1] assembled 17 central questions for electrolyte systems that are not less important for modeling non-aqueous systems. We address and try to give answers to the most important (in our opinion) and still open questions concerning the changeovers required for organic systems, i.e.:

- Do we need the Born term in modeling?
- Should the dielectric constant be considered constant or should we account for the concentration dependency including all derivatives?
- How many and which terms should we use in electrolyte equations of state (EoS)? For engineering applications?
- Which ion diameters should we use?

In this work, the decrease in the dielectric constant of the continuum medium upon changing from water solvent to organic medium is accounted for. Ion-solvent interactions lower the dielectric constant as applied for the continuum in primitive models. This is by now accounted for by a concentration-dependent dielectric constant based on the ion concentration. [2, 3] The Born term is a continuous solvation model to calculate the solvation energy of ions, related to the transfer of an ionic species from ideal gas (vacuum) into a medium. Ions are considered as charged spheres in a primitive continuum solvent, expressed as dielectric constant $\varepsilon=\varepsilon_{0} \varepsilon_{r}$ i.e., product of vacuum dielectric constant multiplied by relative permittivity. In fact, the Born contribution and the Debye-Hückel theory were originally derived based on exactly the same physical framework, and the Born term is an auxiliary constraint in that theoretical framework. The thermodynamic model used in this work is electrolyte Perturbed-Chain Statistical Associating Fluid Theory (ePC-SAFT), the extension to PC-SAFT for electrolyte systems using the Debye-Hückel theory. ePC-SAFT was applied so far to optimize separations and reactions of mainly aqueous systems. [4-10] It turned out that the incorporation of the Born term to treat solvation energies was not necessary with purecomponent parameters for ions directly regressed to aqueous thermodynamic properties. [11] Additionally, if the dielectric constant is independent of the system composition, the Born term will not contribute to excess properties. Still, the utilization of the Born term is a widely discussed topic. Thermodynamic approaches such as ePC-SAFT and other models treated solvation classically by dispersion forces or association forces between ion and solvent. However, it was shown that the Born term has the main contribution to the total free energy of solvation for alkali halides. [12] Research from various groups applying EoS previously included the Born term into their methodology; inter alia the EoS Peng-Robinson [13, 14], Soave-Redlich-Kwong, [15] electrolyte Cubic-Plus-Association (eCPA) [15-17], polar ePPC-SAFT, [12, 18] and SAFT-VRE (variable range electrolyte). [19] Recent works addressed modeling electrolyte solutions containing low-permittivity solvents, [20] and the usefulness of the Born term was discussed also in a review by Held [21] and by the several other work. [1, 17]

Until now, the application of a continuous approach for a concentration-dependent dielectric constant included in the Debye-Hückel theory and in the Born term is not available in the open literature. Continuous means that the pure-ion parameters in ePC-SAFT revised are still applicable and that the model reduces to the original ePC-SAFT revised model if the concentration dependence of the dielectric constant is omitted. This manuscript will

![](https://cdn.mathpix.com/cropped/5ee97251-f053-43e8-9915-6210d4b452e4-2.jpg?height=512&width=686&top_left_y=193&top_left_x=1156)
Figure 1. The Born cycle for calculating the electrostatic contribution to Helmholtz energy of solvation. The total solvation Helmholtz energy $\Delta^{\text {Solv }} A_{i}(1-4)$ is given as sum of three contributions: discharging the ion in vacuum $\Delta^{\text {Disch }} A(1-2)$, transferring the neutral ion in the solvent $\Delta^{T r} A$ (2-3) and recharging the ion $\Delta A^{\text {Charge }}(3-4)$.

show that incorporating both of such effects is of vital importance towards successful modeling properties of non-aqueous electrolyte systems. The new methodology within the ePC-SAFT framework is validated against finite and infinite dilution properties: i) exclusive ion-solvent interactions are governed by the infinite-dilution Gibbs energy of hydration and Gibbs energy of transfer; ion-ion interactions do not perturb such systems. ii) Salt-concentration dependent activity coefficients (MIAC) give inside into the goodness of transferability of the developed theory to non-aqueous electrolyte systems; MIACs are of more industrial and applicative importance due to the dependence on salt and solvent concentration.

## 2. Thermodynamic theories: ePC-SAFT and the altered Born contribution

### 2.1. Ion solvation and thermodynamic properties

The Born term for ion solvation is derived as the transfer work of an initially charged sphere from infinite distance to a first neutral species in the solvent that is afterwards recharged. The work is then summed for each ion in the system to account for the salt concentration dependence. In fact, the Born term may be visualized as a cycle of energy needed to discharge, transfer and recharge an ion as depicted in Figure 1.

The main contribution to the Helmholtz energy in the cycle is the charge process according to Eq. (1).

$$
\begin{equation*}
\Delta A^{\text {Charge }}=\frac{N_{a} e^{2}}{4 \pi \varepsilon} \sum_{i} \frac{n_{i} z_{i}^{2}}{a_{i}} \tag{1}
\end{equation*}
$$

The energy needed for the formation of cavities in the solvent and the formation of solvation structures is generally small compared to the electrostatic energy [12, 22, 23] and is reasonably neglected. The reduced molar Helmholtz free energy contribution of the Born term is therefore obtained from the discharging in vacuum and recharging of the solvent as shown in Eq. (2).

$$
\begin{equation*}
a^{\text {Born }}=\frac{A^{\text {Born }}}{k T N}=-\frac{e^{2}}{4 \pi \varepsilon_{0} k_{B} T}\left(1-\frac{1}{\varepsilon_{r}}\right) \sum_{i} \frac{x_{i} z_{i}^{2}}{a_{i}} \tag{2}
\end{equation*}
$$

In this work, the Born radii $a_{i}$ were set to the respective diameters of the ions estimated within the original ePC-SAFT parametrization method. Thus, new parameters are not required. Using a concentration-dependent dielectric constant $\varepsilon_{r}(\bar{x})$, the derivative with respect to the composition of the system are changed, yielding the further called 'altered Born contribution'. The partial derivative of the reduced molar altered Born contribution is derived in this work and presented in Eq. (3).

$$
\begin{align*}
\left(\frac{\partial a^{\text {Born }}}{\partial x_{i}}\right)_{T, v_{N}, x_{j \neq i}}= & -\frac{e^{2}}{4 \pi k_{B} T \varepsilon_{0}}\left(1-\frac{1}{\varepsilon_{r}(\bar{x})}\right) \frac{z_{i}^{2}}{a_{i}} \\
& -\frac{e^{2}}{4 \pi k_{B} T \varepsilon_{0}}\left(\frac{1}{\varepsilon_{r}(\bar{x})^{2}}\right)\left(\frac{\partial \varepsilon_{r}(\bar{x})}{\partial x_{i}}\right)_{T, v, x_{j \neq i}} \tag{3}
\end{align*}
$$

An expression for the altered Born contribution with respect to the residual pressure and with respect to the particle volume is automatically set to zero as the pressure dependence of the dielectric constant and of the Born radii are neglected (Eq. (4)).

$$
\begin{equation*}
\frac{p^{\text {Born }}}{k_{B} T}=\left(\frac{\partial a^{\text {Born }}}{\partial v_{N}}\right)_{T, v_{N}, x_{j \neq i}}=0 \tag{4}
\end{equation*}
$$

In contrast, temperature dependence of the dielectric constant cannot be neglected. From fundamental thermodynamics (Eq. (5)) an expression for the enthalpy can be derived.

$$
\begin{equation*}
H=A+P V+T S \tag{5}
\end{equation*}
$$

Here, the Born contribution to pressure (PV) is equal to zero and the partial derivation of the Helmholtz energy over temperature is equal to the negative entropy. Eq. (6) expresses the reduced molar enthalpy caused by the ion solvation.

$$
\begin{equation*}
h^{\text {Born }}=\frac{H^{\text {Born }}}{k T N}=-\frac{e^{2}}{4 \pi \varepsilon_{0} k_{B} T}\left(1-\frac{1}{\varepsilon_{r}}-\frac{T}{\varepsilon_{r}^{2}} \frac{\partial \varepsilon_{r}}{\partial T}\right) \sum_{i} \frac{x_{i} z_{i}^{2}}{a_{i}} \tag{6}
\end{equation*}
$$

The Gibbs free energy of solvation is calculated from fugacity coefficients of the ion $i$ at infinite dilution in the respective solvent S1 (Eq. (7)).

$$
\begin{equation*}
\Delta^{H y d} G_{i}(T, p)=R T \cdot \ln \varphi_{i}^{\infty, S 1}(T, p) \tag{7}
\end{equation*}
$$

Consequently, the Gibbs energy of transfer is calculated from the ratio of the fugacity coefficients at infinite dilution in solvents $S 1$ and $S 2$.

$$
\begin{equation*}
\Delta^{T r} G_{i}=R T \cdot \ln \left(\frac{\varphi_{i}^{\infty, S 2}}{\varphi_{i}^{\infty, S 1}}\right) \tag{8}
\end{equation*}
$$

Further, the MIAC is available experimentally as finite dilution property. If a salt is assumed to be fully dissociated, the MIAC is calculated from the ion-based activity coefficients of the cation cat and anion an as depicted in Eq. (9).

$$
\begin{equation*}
\gamma_{ \pm}^{*}=\left(v_{c a t} \ln \gamma_{c a t}^{*}+v_{a n} \ln \gamma_{a n}^{*}\right)^{\frac{1}{v_{c a t}+v_{a n}}} \tag{9}
\end{equation*}
$$

The activity coefficient of each of the ions $i$ in the mixture is calculated from the ration of the fugacity coefficients in mixture and of the fugacity coefficient at infinite dilution of the ion.

$$
\begin{equation*}
\gamma_{i}^{*}=\frac{\varphi_{\mathrm{i}}(\mathrm{~T}, \mathrm{p}, \overrightarrow{\mathrm{x}})}{\varphi_{0 \mathrm{i}}\left(\mathrm{~T}, \mathrm{p}, \mathrm{x}_{\mathrm{i}} \rightarrow 0\right)} \tag{10}
\end{equation*}
$$

## 2.2. ePC-SAFT advanced

In this work, the infinite-dilution fugacity coefficient $\varphi_{i}^{\infty}$ and the MIAC $\gamma_{ \pm}^{*}$ of the salt in the respective solvent (water, organic solvent) were calculated with the advanced ePC-SAFT approach, see Eq. (11).

$$
\begin{equation*}
a^{\text {res }}=a^{h c}+a^{\text {disp }}+a^{\text {assoc }}+a^{\text {DH }}\left(\varepsilon_{r}(\bar{x})\right)+a^{\text {Born }}\left(\varepsilon_{r}(\bar{x})\right) \tag{11}
\end{equation*}
$$

The model including $\varepsilon_{r}(\bar{x})$ and the required deviations will from now on be denoted ePC-SAFT advanced. ePC-SAFT advanced requires five pure-component parameters for associating compounds like water: segment number $m_{i}^{\text {seg }}$, diameter $\sigma_{i}$, dispersion energy $u_{i}$, association energy $\varepsilon^{\text {AiBi }}$ and volume $\kappa^{\text {AiBi }}$, while only two purecomponent parameters for the ions are applied: diameter $\sigma_{i}$ and
dispersion energy $u_{i}$. To calculate mixture properties, the combining rules of Berthelot-Lorentz were applied Eqs. (12) and ((13)).

$$
\begin{equation*}
\sigma_{i j}=\frac{1}{2}\left(\sigma_{i}+\sigma_{j}\right) \tag{12}
\end{equation*}
$$

$$
\begin{equation*}
u_{i j}=\sqrt{u_{i} u_{j}}\left(1-k_{i j}\right) \tag{13}
\end{equation*}
$$

Eq. (13) introduces the binary interaction parameter $k_{i j}$ that can be used to alter the dispersion energy $u_{i j}$ in the mixture. For electrolyte solutions, $k_{\text {ion }}$, solvent as well as $k_{\text {cat }}$, an are required for systems salt + water with pure-component parameters of the ions regressed including the interaction parameters. [11] For the prediction of infinite dilution properties and MIACs in the organic solvent, $\mathrm{k}_{\mathrm{ij}}$-values between ion and solvent were set to zero. Thus, the results are denoted 'predictions'. By derivation of the residual Helmholtz energy with respect to density and mole fraction, the fugacity coefficient is obtained (Eq. (14)).

$$
\begin{equation*}
\ln \left(\varphi_{\mathrm{i}}\right)=\frac{\mu_{\mathrm{i}}^{\mathrm{res}}}{\mathrm{k}_{\mathrm{B}} \cdot \mathrm{~T}}-\ln \left(1+\left(\frac{\partial\left(\frac{\mathrm{a}^{\mathrm{res}}}{\mathrm{k}_{\mathrm{B}} \cdot \mathrm{~T}}\right)}{\partial \rho}\right)\right) \tag{14}
\end{equation*}
$$

## 3. Results and Discussion

### 3.1. Pure-component and binary interaction parameters

The advanced treatment including the altered Born contribution and the concentration-dependent dielectric constant within ePC-SAFT aims at using already published pure-component parameters for all components and ions. The pure-component parameters for the organic solvents and water are presented in Table 1. All solvents were modeled as associating fluids with a 2B associating scheme. The inorganic ions representing alkali-halide salts are presented in Table 2 together with the binary interaction parameters $k_{\text {ion, water }}[11]$ and interaction parameters between the inorganic ions Table 3). In ePC-SAFT revised [11], all these parameters were estimated exclusively to properties (e.g. MIACs) in water solvent.

A transferability study to non-aqueous systems will be presented in this work. The already existing interaction parameters were applied in the present work to predict infinite and finite properties in organic solvents. That is, further $k_{i j}$-values were not used i.e., $k_{i j}=0$ for ion-organic solvent pairs.

The dielectric constant is a fundamental part of the modelling within the new treatment for ePC-SAFT. The relative dielectric constants $\varepsilon_{r}$ applied in this work are summarized in Table 4. A linear concentration dependence of $\varepsilon_{r}$ with the mole-based salt concentration is assumed according to previous work. [2] In that work, the linear relation has proven to be feasible for modeling LLE of water and ionic liquids. Certainly, this is not a perfect representation of the actual dependence with concentration and is much more a practical simplification with the application to engineering models in mind. Beneficially, this approach does not rely on experimental data for dielectric constants in mixtures. Even for binary systems experimental data is rarely available. For multi-component systems the effort for the determination of the dielectric constant for mixtures is enormously in terms of time and effort. The expression linear in concentration shows another benefit of the new approach. The derivative of the mixing rule for the dielectric constant is always the value of the pure component for all species involved in the system, $\frac{\partial\left(\varepsilon_{r}(\bar{x})\right)}{\partial x_{i}}=\varepsilon_{r, i}$. For the majority of industrially important solvent classes and even for inorganic salts, experimental data for the dielectric constant of the pure compound are available in the literature. Bearing in mind the experimental lack for multi component systems, the application of a simple mixing rule is feasible.

Table 1
Pure-component parameters of water and organic solvents used in this work. All components have a 2B association scheme.
| Organic Solvent | $m_{i}^{\text {seg }}$ | $\sigma_{i} / \AA$ | $u_{i} / k_{B} / \mathrm{K}$ | $\varepsilon^{\text {AiBi }} / \mathrm{K}$ | $\kappa^{\text {AiBi }}$ | Ref. |
| :--- | :--- | :--- | :--- | :--- | :--- | :--- |
| Water | 1.2047 | $*$ | 353.95 | 2425.7 | 0.04509 | $[24]$ |
| Methanol | 1.5255 | 3.2300 | 188.90 | 2899.5 | 0.03518 | $[25]$ |
| Ethanol | 2.3827 | 3.1771 | 198.24 | 2653.4 | 0.03238 | $[25]$ |


* $\sigma=2.7927+\left(10.11 \cdot \mathrm{e}^{-0.01775 \mathrm{~T} / \mathrm{K}}-1.417 \cdot \mathrm{e}^{-0.01146 \mathrm{~T} / \mathrm{K}}\right)$

Table 2
ePC-SAFT pure-ion parameters and binary interaction parameters $k_{i j}$ between water and inorganic ions from ref. [11] Segment number is unity for all ions. $k_{i j}$ valid only at 298.15 K .
| Ion | $\sigma_{i} / \AA$ | $u_{i} / k_{B} / \mathrm{K}$ | $k_{\text {ion, water }}$ |
| :--- | :--- | :--- | :--- |
| $\mathrm{Li}^{+}$ | 2.8449 | 360.00 | -0.2500 |
| $\mathrm{Na}^{+}$ | 2.8232 | 230.00 | 0.0045 |
| $\mathrm{~K}^{+}$ | 3.3417 | 200.00 | 0.1997 |
| $\mathrm{~F}^{-}$ | 1.7712 | 275.00 | 0.0000 |
| $\mathrm{Cl}^{-}$ | 2.7560 | 170.00 | -0.2500 |
| $\mathrm{Br}^{-}$ | 3.0707 | 190.00 | -0.2500 |
| $\mathrm{I}^{-}$ | 3.6672 | 200.00 | -0.2500 |


Table 3
ePC-SAFT binary interaction parameters $k_{i j}$ between halide anions and alkali cations from ref. [11]
|  | $\mathrm{Li}^{+}$ | $\mathrm{Na}^{+}$ | $\mathrm{K}^{+}$ |
| :--- | :--- | :--- | :--- |
| $\mathrm{F}^{-}$ | 0 | 0.665 | 1.000 |
| $\mathrm{Cl}^{-}$ | 0.669 | 0.317 | 0.064 |
| $\mathrm{Br}^{-}$ | 0.591 | 0.290 | -0.102 |
| $\mathrm{I}^{-}$ | 0.002 | 0.018 | -0.312 |


Table 4
Dielectric constants for solvents and salts applied in this work at 298.15 K and 1 bar.
| Component | Dielectric constant $/ \mathrm{C} \cdot \mathrm{Vm}^{-1}$ | Ref. |
| :--- | :--- | :--- |
| Water | 78.09 | $[30]$ |
| Methanol | 33.05 | $[29]$ |
| Ethanol | 24.88 | $[28]$ |
| Salts | 8 | a) |


${ }^{\text {a) }}$ All salts were modeled with a similar dielectric constant that is a mean of available experimental data.[27]

It has to be mentioned that there are theoretical models able to predict the dielectric constant in a (binary) mixture with at least qualitative agreement. [26] With a future model performing on a quantitative level, the above-mentioned benefits might be dropped in favor of a more physical meaningful approach coupled with a numerical derivation with respect to the mole fraction of component $\frac{\partial\left(\varepsilon_{r}(\bar{x})\right)}{\partial x_{i}}$. All salts were modeled with the same value of $\varepsilon_{r}=8$ that represents the value for most of the salts ( $7<\varepsilon_{r}<9$ ) in good estimation. [27] The dielectric constant for solvents was taken from literature and are here given for the temperature 298 K and 1 bar. [28, 29]

### 3.2. Gibbs energy of solvation

The advanced ePC-SAFT model was tested against literature data on Gibbs energy of solvation and Gibbs energy of transfer. At infinite dilution, the only interactions involving ions are ionsolvent interactions, reducing the altered Born contribution to the original Born term. Both, $\Delta^{\text {Hyd }} G_{i}$ and $\Delta^{\text {Tr }} G_{i}$ are important properties; successfully modeling both such properties provides a meaningful validation of the new approach within ePC-SAFT including

![](https://cdn.mathpix.com/cropped/5ee97251-f053-43e8-9915-6210d4b452e4-4.jpg?height=393&width=675&top_left_y=551&top_left_x=1163)
Figure 2. Gibbs energy of hydration at 298 K and 1 bar at infinite dilution ePCSAFT revised ${ }^{11}$ : orange; ePC-SAFT advanced: green; SAFT-VR: violet [19]. Literature data: gray, from Fawcett et al. [31]

the assumption to replace the Born radii with the ePC-SAFT ion diameters $\sigma_{i} . \Delta^{\text {Hyd }} G_{i}$ values for different alkali-halide ions have been correlated by Fawcett et al. [31] from experimental data and used for validation of SAFT-VRE by Schreckenberg et al. [19] In Figure 2, the literature $\Delta^{\text {Hyd }} G_{i}$ values for the single alkali cations and halide anions are compared to SAFT-VRE and to ePC-SAFT advanced (green bars) from this work and to ePC-SAFT revised (orange bars) from 2014. [11]

The authors would like to stress that the literature $\Delta^{\text {Hyd }} G_{i}$ data are based on correlations of effects averaged from considering different salts. Certainly, $\Delta^{\text {Hyd }} G_{i}$ values for single ions cannot be experimentally measured. Thus, a conclusion on the quantitative accordance between such data and the SAFT modeling approaches is rather complicated. Nonetheless, the modeling results give a good indication of the correctness of the different modeling approaches. The data show that the Gibbs energy of hydration increases for the alkali cations from $\mathrm{Li}^{+}$to $\mathrm{K}^{+}$and for the halide anions from $\mathrm{F}^{-}$to $\mathrm{I}^{-}$. The increase correlates with the size of the ions - the smallest ions possess the highest surface charge density; thus, the smaller the ion the higher the water-ion interactions at infinite dilution. Obviously, ePC-SAFT revised does not allow predicting $\Delta^{\text {Hyd }} G_{i}$; ePC-SAFT revised yields results of the wrong sign for all ions except $\mathrm{Li}^{+}$. This behavior is the reason behind the use of large binary interaction parameters to correlate distribution coefficients of salts in water/MIBK mixed solvent [32] or in water/polymer two phase systems [9]. In contrast, ePC-SAFT advanced from this work and SAFT-VRE estimate much more meaningful $\Delta^{\text {Hyd }} G_{i}$ values; both predict more negative values than the literature data. $\Delta^{\text {Hyd }} G_{i}$ for the alkali cations is more negative for SAFT-VRE than for ePC-SAFT advanced, and vice versa for the halide anions. The resulting mean Gibbs energy of hydration (salt specific) predicted by both models, SAFT-VRE and ePC-SAFT advanced, are roughly equal for all alkali halides at the conditions under study (not graphically visualized).

Further, a detailed investigation into the different contributions to the Gibbs energy of hydration is presented in Figure 3. The predominant contribution to $\Delta^{\text {Hyd }} G_{i}$ is coming from the Born term. This has two implications. First, simplifications within deriving the Born term as discussed earlier (cf. Figure 1 and text) might be seen critical. The application of the diameter $\sigma_{i}$ instead of the Born radii

![](https://cdn.mathpix.com/cropped/5ee97251-f053-43e8-9915-6210d4b452e4-5.jpg?height=396&width=679&top_left_y=193&top_left_x=228)
Figure 3. Contributions to the Gibbs energy of hydration $\Delta^{\text {Hyd }} \boldsymbol{G}_{\boldsymbol{i}}$ of ions at infinite dilution within ePC-SAFT advanced from this work. For each alkali-halide ion, the Born term (orange) provides the largest contribution. Debye-Hückel does not contribute at infinite dilution. Gray: hard-chain contribution, dark gray: dispersion contribution, black: association contribution.

or any of the other ion radii available in the literature might be a reason for the deviations between ePC-SAFT advanced and the literature $\Delta^{\text {Hyd }} G_{i}$. Nevertheless, Schreckenberg et al. applied Pauling radii, which also did not allow to further minimize the deviations between SAFT-VRE and the literature $\Delta^{\text {Hyd }} G_{i}$. In fact, $\Delta^{\text {Hyd }} G_{i}$ data might be useful to regress more meaningful radii, which is not done in this work due to the desire to use a minimum number of model parameters.

Second, and more generally valid, neglecting the Born term is not recommended for future work with ePC-SAFT especially for modeling non-aqueous electrolyte systems. Upon neglecting the Born term, huge $k_{i j}$-values between ion and solvent or even solvent-dependent ion parameters were necessary in other works [33, 34]. It should be noted here again that $\Delta^{\text {Hyd }} G_{i}$ is an infinitedilution property. The concentration-dependency of the altered Born contribution and the dielectric constant does not play a role for modeling $\Delta^{\text {Hyd }} G_{i}$; the importance of the concentration dependence will be shown in the section "MIAC in organic solvents".

### 3.3. Gibbs energy of transfer

The Gibbs energy of transfer ( $\Delta^{T r} G_{i}$ ) represents an important property for the transfer from aqueous to non-aqueous electrolyte systems. $\Delta^{T r} G_{i}$ values are experimentally available e.g., for alkalihalide ions from water to methanol or to ethanol. The comparison of ePC-SAFT modeling results for $\Delta^{T r} G_{i}$ with vs. without the Born term is compared to experimental data in Figure 4. $\Delta^{T r} G_{i}$ values for the halide anions estimated with ePC-SAFT advanced are in good agreement with the experimental data while ePC-SAFT revised (without the Born term) largely underestimates and even qualitatively incorrect represents the experimental data. This indicates that the Born term allows correcting the ion-solvent interactions in a manner to bring the modeling results to the same order of magnitude compared to the experimental $\Delta^{T r} G_{i}$ values. Still, the variation of $\Delta^{T r} G_{i}$ with increasing ion diameter is not correctly represented with ePC-SAFT advanced for the cations. High deviations are still observed for $\mathrm{Li}^{+}$and $\mathrm{K}^{+}$, which might be attributed to the potentially too high dispersion energy between ion and water. However, this work does not aim at readjusting model parameters.

The single contributions to the Gibbs energy of transfer within ePC-SAFT advanced for $\mathrm{Na}^{+}, \mathrm{Cl}^{-}$and $\mathrm{I}^{-}$for the transfer from water to methanol and ethanol, respectively, is depicted in Figure 5. The most important aspect here is that the Born term contributes significantly for the systems under study.

The ePC-SAFT diameter $\sigma_{i}$ for $\mathrm{Na}^{+}$and $\mathrm{Cl}^{-}$are almost similar, which explains the similar values for $\Delta^{T r} G_{N a+}$ and $\Delta^{T r} G_{C l-}$, while $\Delta^{T r} G_{I-}$ is significantly smaller caused by the larger $\sigma_{i}$-value. The percental contribution of the Born term on $\Delta^{T r} G_{i}$ is much smaller
than on $\Delta^{\text {Hyd }} G_{i}$. This can be explained by the Born self-energy that balances out for the transfer from water to the alcohol. Still, neglecting the Born term would translate to an oversimplification, which might only be compensated by applying unphysically high binary interaction parameters $k_{i j}$ - a strategy that was used in previous works [32] but which is not further recommended any more.

### 3.4. MIAC in alcoholic solvents

Gibbs energies of hydration and Gibbs energy of transfer present infinite-dilution properties. Experimental MIACs for salts in different solvents [37, 38] are suitable to validate ePC-SAFT advanced and to proof the physical meaning of including a concentration-dependent dielectric constant in the altered Born contribution and in the Debye-Hückel theory (cf. Figure 6 and Figure 7). Originally in ePC-SAFT revised, [11] MIACs at 298.15 K were used as input data to the parameter estimation. The following results were achieved using these originally-determined parameters. The methods aim at transferring all the parameters to predict MIACs in non-aqueous solvents.

Figure 6a illustrates the MIACs of LiBr in ethanol for different modeling prediction approaches within the ePC-SAFT framework, all of them without the use of $k_{i j}$-parameters between ion and solvent: ePC-SAFT revised, ePC-SAFT revised with the Born term, ePCSAFT revised but included with a concentration-dependent dielectric constant in the Debye-Hückel theory and ePC-SAFT advanced. Prediction of MIACs fail by using ePC-SAFT without and with the Born term as long as the concentration-dependence $\varepsilon_{r}(\bar{x})$ is excluded. The reason is that excess properties such as the MIACs are not affected by the Born term as long as the dielectric constant is of constant value; the Born contribution cancels out e.g., in Eq. (10). Neglecting the Born term for modeling aqueous systems without a consideration of a concentration-dependent dielectric constant is therefore not only a simplification, but also thermodynamically and mathematically consistent. Predictions with quantitative agreement over the whole investigated concentration range of the salt were only achieved by incorporating the concentrationdependent dielectric constant into the altered Born contribution and Debye-Hückel theory (ePC-SAFT advanced). A detailed inspection of the different contributions to the MIAC within the ePCSAFT advanced model is exemplarily shown for LiBr in ethanol in Figure 6b. The short-range forces i.e., the classical physical forces described by original PC-SAFT, contribute with a linear rise to the MIACs with increasing concentration of salt; however, the contribution to the MIAC is surprisingly low. Certainly, the electrostatic forces are expected to play a more fundamental role in organic vs. aqueous systems. Debye and Hückel stated that their treatment is only applicable to a physical model that is fully developed and incorporates all necessary short-range forces. For MIACs of the salts $(\mathrm{LiBr}, \mathrm{LiCl}, \mathrm{NaBr})$ under study in ethanol and methanol, this statement might need to be revised as long as water is not present still a surprising conclusion. Further, this might even be an argument to not adjust ion-related model parameters to data of systems alcohol + salt. Nevertheless, more complex systems will certainly require reasonable models for the short-range interactions e.g., liquid-liquid equilibria or solubility of salts in solvents or solvent mixtures. The Debye-Hückel theory gives a monotone nonlinear decrease that is responsible for the negative value of the MIAC at small salt concentration ranges (usually $<0.1 \mathrm{~mol} / \mathrm{kg}$ ). Only utilizing the altered Born contribution allows describing the increase of MIACs with increasing salt concentration - the key towards quantitatively predicting the MIACs.

Similar accuracy was obtained also for MIACs of other salts in alcohol (c.f. Figure 7). This reassures that the pure-ion ePC-SAFT parameters are still applicable also to non-aqueous systems, as long as the advanced theory within ePC-SAFT is applied. This is

![](https://cdn.mathpix.com/cropped/5ee97251-f053-43e8-9915-6210d4b452e4-6.jpg?height=403&width=1438&top_left_y=211&top_left_x=316)
Figure 4. Gibbs energy of transfer $\Delta^{\boldsymbol{T r}} \boldsymbol{G}_{\boldsymbol{i}}$ of alkali-halide ions at infinite dilution for the transfer from water to alcohol at 298.15 K and 1 bar. Comparison between ePC-SAFT advanced: green; this work) and ePC-SAFT revised: [11] orange with experimental data from Marcus et al.: gray bars. [35, 36] Left: $\Delta^{\boldsymbol{T r}} \boldsymbol{G}_{\boldsymbol{i}}$ from water to methanol, right: from water to ethanol.

![](https://cdn.mathpix.com/cropped/5ee97251-f053-43e8-9915-6210d4b452e4-6.jpg?height=912&width=1448&top_left_y=766&top_left_x=312)
Figure 5. Different contributions within ePC-SAFT advanced (bars: hard chain, dispersion, association, Debye-Hückel, Born) to the Gibbs energy of transfer for ions from water to methanol and from water to ethanol at infinite dilution at 298.15 K and 1 bar. Blue: $\mathrm{Na}^{+}$, red: $\mathrm{Cl}^{-}$or $\mathrm{I}^{-}$. Top: ethanol; left: NaCl , right: NaI . Bottom: methanol; left: NaCl, right: NaI.The contribution is in the same order of absolute magnitude compared to contributions from hard-chain repulsion, dispersion, and association. While the contributions due to hard-chain repulsion and dispersion more or less balance out, the Born term forces $\Delta^{\boldsymbol{T r}} \boldsymbol{G}_{\boldsymbol{i}}$ to become positive - a realistic value. Further, the results show that using ePC-SAFT ion diameter as input to the Born term gives reasonable results.

![](https://cdn.mathpix.com/cropped/5ee97251-f053-43e8-9915-6210d4b452e4-6.jpg?height=533&width=1450&top_left_y=1899&top_left_x=310)
Figure 6. Mean ionic activity coefficients (MIACs) of LiBr in ethanol at 298 K and 1 bar. Circles are experimental data. [23]. Left: Black: ePC-SAFT revised, [11] orange: ePC-SAFT + Born without $\boldsymbol{\varepsilon}\left(\boldsymbol{x}_{\text {salt }}\right)$, gray: ePC-SAFT with $\boldsymbol{\varepsilon}\left(\boldsymbol{x}_{\text {salt }}\right)$ in DH without Born, green: ePC-SAFT advanced with Born and $\boldsymbol{\varepsilon}\left(\boldsymbol{x}_{\text {salt }}\right)$ in Born and in DH. Right: Single contributions within ePC-SAFT advanced to the predicted MIAC. Orange: Born contribution with $\boldsymbol{\varepsilon}\left(\boldsymbol{x}_{\text {salt }}\right)$; blue: Debye-Hückel contribution $\boldsymbol{\varepsilon}\left(\boldsymbol{x}_{\text {salt }}\right)$. Gray solid: hard-chain contribution; gray dashed: dispersion contribution and gray dashed-dotted: association contribution.

![](https://cdn.mathpix.com/cropped/5ee97251-f053-43e8-9915-6210d4b452e4-7.jpg?height=1022&width=1452&top_left_y=176&top_left_x=310)
Figure 7. Mean ionic activity coefficients (MIACs) of alkali-halide salts in alcohols methanol and ethanol at 298 K and 1 bar. Cycles are experimental data. [37, 38] Green: Predictions with ePC-SAFT advanced, black: ePC-SAFT revised. Top: ethanol; left: LiBr, right: NaBr. Bottom: methanol; left: LiCl, right: NaBrConclusion

an extremely important discovery, meaning the new treatment not only allows the transferability of pure-ion parameters, which is highly debated and that was at no time verified. Effectively, this work shows the successfully achieved transition of ePC-SAFT to non-aqueous electrolyte systems. To the best of our knowledge, this is the first possible transition of any (electrolyte) EoS to nonaqueous solutions without considering large interaction parameters or the need for re-estimating new pure-component parameters.

It is now possible to answer some of the questions of Kontogeorgis et al. [1] based on the modeling results within ePC-SAFT advanced gained in the present work:

- Yes, we need the Born term in modeling, but only combined with a concentration-dependent dielectric constant including all derivatives, c.f. altered Born contribution.
- We should use as short-range and long-range contributions in electrolyte equations of state (EoS), as engineering applications aim at complex systems involving a cocktail of diverse species.
- The use of unified ion diameters that reflect the distance of closest approach in Debye-Hückel theory as well as in the altered Born contribution is apparently feasible.

Further, including infinite-dilution properties of organic solvents into the parameter estimation might be an interesting option, already implemented by other models. [20] Still, such data scatter significantly (c.f. Figure 4). Nevertheless, including such data would increase the available data set. As a final statement, a re-estimation of the pure-ion parameters or an inclusion of nonunified ion parameters (different ion-size parameters for DebyeHückel and altered Born) within ePC-SAFT is not intended, given the success of ePC-SAFT advanced and the required effort for a reparametrization.

## 4. Conclusion

This work advances the ePC-SAFT framework by incorporating a concentration-dependent dielectric constant into the altered Born contribution and in the Debye-Hückel theory to the residual Helmholtz energy. By this ePC-SAFT is capable to describe physical phenomena in solvents with low dielectric constant. The altered Born contribution addresses increased ion solvation by low permittivity, which is induced by decreased dielectric constant upon increasing salt concentration. The impact of the Born contribution is evaluated against Gibbs energy of hydration and Gibbs energy of transfer for alkali cations and halide anions. Such data are meaningful indicators for the physical soundness of any electrolyte thermodynamic model. It was shown that the altered Born contribution is of great importance for the transfer of the ePC-SAFT framework and the ion parameters from aqueous to non-aqueous systems. This was proven by the large percentage of the Born contribution to both, Gibbs energy of hydration and Gibbs energy of transfer. Further, the advanced ePC-SAFT framework was validated against experimental MIACs in alcohol. The predictions for MIACs were found to be in very good agreement with experimental data. That is, i) pure-component parameters for ions, estimated originally using aqueous properties, are transferable non-aqueous systems, ii) the importance of the Born term is significant only if the concentration-dependent dielectric constant is used. Overall, the transition of ePC-SAFT to a model for non-aqueous systems is accomplished. A re-estimation of the applied pure-component parameters for ions is ruled out by the accurate prediction of MIACs and by the lack in precise experimental data for infinite-dilution properties. The consequences of the approach presented in this work is the need to rethink classical approaches upon model development in the future: more pronounced physical phenomena cannot be neglected in electrolyte solutions instead of compensating such effects and hiding them by fitting many model parameters.

Further, four central questions arisen in the work of Kontogeorgis et al. have been addressed and intensively debated. Further, including the altered Born contribution brings a first step towards universal modelling, solvent-independent, predictive, and extrapolative in salt concentration.

## Declaration of Competing Interest

The authors declare no conflict of interest.

## CRediT authorship contribution statement

Mark Bülow: Investigation, Validation, Formal analysis, Writing - original draft, Conceptualization, Methodology. Moreno Ascani: Investigation, Validation, Formal analysis. Christoph Held: Conceptualization, Methodology, Project administration, Writing - review \& editing, Supervision.

## Acknowledgement

This work has been supported by the German Science Foundation (DFG) within the priority program SPP 1708 "Material Synthesis Near Room Temperature" (grant HE 7165/7-1).

## References

[1] G.M. Kontogeorgis, B. Maribo-Mogensen, K. Thomsen, The Debye-Hückel theory and its importance in modeling electrolyte solutions, Fluid Phase Equilibria 462 (2018) 130-152.
[2] M. Bülow, X. Ji, C. Held, Incorporating a concentration-dependent dielectric constant into ePC-SAFT. An application to binary mixtures containing ionic liquids, Fluid Phase Equilib 492 (2019) 26-33.
[3] M. Bülow, A. Danzer, C. Held, Liquid-Liquid Equilibria of Binary and Ternary Systems Containing Ionic Liquids, in: S Zhang (Ed.), Encyclopedia of Ionic Liquids, Springer Singapore, Singapore, 2020, pp. 1-7.
[4] C. Held, G. Sadowski, Thermodynamics of Bioreactions, Annual review of chemical and biomolecular engineering 7 (2016) 395-414.
[5] C. Held, T. Reschke, R. Müller, W. Kunz, G. Sadowski, Measuring and modeling aqueous electrolyte/amino-acid solutions with ePC-SAFT, The Journal of Chemical Thermodynamics 68 (2014) 1-12.
[6] M. Bülow, A. Schmitz, T. Mahmoudi, D. Schmidt, F. Junglas, C. Janiak, C. Held, Odd-even effect for efficient bioreactions of chiral alcohols and boosted stability of the enzyme, RSC Adv 10 (2020) 28351-28354.
[7] K. Wysoczanska, E.A. Macedo, G. Sadowski, C. Held, Solubility Enhancement of Vitamins in Water in the Presence of Covitamins: Measurements and ePC-SAFT Predictions, Ind. Eng. Chem. Res. 58 (2019) 21761-21771.
[8] T. Reschke, C. Brandenbusch, G. Sadowski, Modeling aqueous two-phase systems: I. Polyethylene glycol and inorganic salts as ATPS former, Fluid Phase Equilibria 368 (2014) 91-103.
[9] M. Wessner, M. Nowaczyk, C. Brandenbusch, Rapid identification of tailormade aqueous two-phase systems for the extractive purification of high-value biomolecules, Journal of Molecular Liquids 314 (2020) 113655.
[10] M. Schleinitz, L. Nolte, C. Brandenbusch, Predicting protein-protein interactions using the ePC-SAFT equation-of-state, Journal of Molecular Liquids 298 (2020) 112011.
[11] C. Held, T. Reschke, S. Mohammad, A. Luza, G. Sadowski, ePC-SAFT revised, Chem. Eng. Res. Des. 92 (2014) 2884-2897.
[12] S. Ahmed, N. Ferrando, J.-C. Hemptinne, J.-P. de; Simonin, O. Bernard, O. Baudouin, Modeling of mixed-solvent electrolyte systems, Fluid Phase Equilibria 459 (2018) 138-157.
[13] J. Wu, J.M. Prausnitz, Phase Equilibria for Systems Containing Hydrocarbons, Water, and Salt: An Extended Peng-Robinson Equation of State, Ind. Eng. Chem. Res. 37 (1998) 1634-1643.
[14] J.A. Myers, S.I. Sandler, R.H. Wood, An Equation of State for Electrolyte Solutions Covering Wide Ranges of Temperature, Pressure, and Composition, Ind. Eng. Chem. Res. 41 (2002) 3282-3297.
[15] Y. Lin, K. Thomsen, J.-C.de Hemptinne, Multicomponent equations of state for electrolytes, AIChE J 53 (2007) 989-1005.
[16] de R. Inchekel, J.-C. Hemptinne, W. Fürst, The simultaneous representation of dielectric constant, volume and activity coefficients using an electrolyte equation of state, Fluid Phase Equilibria 271 (2008) 19-27.
[17] B. Maribo-Mogensen, K. Thomsen, G.M. Kontogeorgis, An electrolyte CPA equation of state for mixed solvent electrolytes, AIChE J 61 (2015) 2933-2950.
[18] de J. Rozmus, J.-C. Hemptinne, A. Galindo, S. Dufal, P. Mougin, Modeling of Strong Electrolytes with ePPC-SAFT up to High Temperatures, Ind. Eng. Chem. Res. 52 (2013) 9979-9994.
[19] J.M.A. Schreckenberg, S. Dufal, A.J. Haslam, C.S. Adjiman, G. Jackson, A Galindo, Modelling of the thermodynamic and solvation properties of electrolyte solutions with the statistical associating fluid theory for potentials of variable range, Molecular Physics 112 (2014) 2339-2364.
[20] S. Müller, A. González de Castilla, C. Taeschler, A. Klein, I. Smirnova, Calculation of thermodynamic equilibria with the predictive electrolyte model COS-MO-RS-ES: Improvements for low permittivity systems, Fluid Phase Equilibria 506 (2020) 112368.
[21] C. Held, Thermodynamic gE Models and Equations of State for Electrolytes in a Water-Poor Medium: A Review, J. Chem. Eng. Data 65 (2020) 5073-5082.
[22] J.N. Israelachvili, Intermolecular and Surface Forces, 3rd ed., Elsevier Science, Saint Louis, 2015.
[23] A.A. Rashin, B. Honig, Reevaluation of the Born model of ion hydration, J. Phys. Chem. 89 (1985) 5588-5593.
$[24]$ D. Fuchs, J. Fischer, F. Tumakaka, G. Sadowski, Solubility of Amino Acids: Influence of the pH value and the Addition of Alcoholic Cosolvents on Aqueous Solubility, Ind. Eng. Chem. Res. 45 (2006) 6578-6584.
[25] J. Gross, G. Sadowski, Application of the Perturbed-Chain SAFT Equation of State to Associating Systems, Ind. Eng. Chem. Res. 41 (2002) 5510-5515.
[26] B. Maribo-Mogensen, G.M. Kontogeorgis, K. Thomsen, Modeling of dielectric properties of aqueous salt solutions with an equation of state, The journal of physical chemistry. B 117 (2013) 10523-10533.
[27] C. Andeen, J. Fontanella, D. Schuele, Low-Frequency Dielectric Constant of LiF, $\mathrm{NaF}, \mathrm{NaCl}, \mathrm{NaBr}, \mathrm{KCl}$, and KBr by the Method of Substitution, Physical Review B (1970) 5068-5073.
[28] R.M. Shirke, A. Chaudhari, N.M. More, P.B. Patil, Temperature dependent dielectric relaxation study of ethyl acetate - Alcohol mixtures using time domain technique, Journal of Molecular Liquids 94 (2001) 27-36.
[29] M.T. Khimenko, V.V. Aleksandrov, N.N. Gritsenko, Polarizability and radii of molecules of some pure liquids, Zh.Fiz.Khim (1973) 2914-2915.
[30] W.B. Floriano, M.A.C. Nascimento, Dielectric constant and density of water as a function of pressure at constant temperature, Braz. J. Phys. 34 (2004) 38-41.
[31] W.R. Fawcett, Thermodynamic Parameters for the Solvation of Monatomic Ions in Water, The journal of physical chemistry. B 103 (1999) 11181-11185.
[32] S. Mohammad, C. Held, E. Altuntepe, T. Köse, G. Sadowski, Influence of Salts on the Partitioning of 5-Hydroxymethylfurfural in Water/MIBK, The journal of physical chemistry. B 120 (2016) 3797-3808.
[33] C. Held, A. Prinz, V. Wallmeyer, G. Sadowski, Measuring and modeling alcohol/salt systems, Chemical Engineering Science 68 (2012) 328-339.
[34] S. Mohammad, C. Held, E. Altuntepe, T. Köse, T. Gerlach, I. Smirnova, G. Sadowski, Salt influence on MIBK/water liquid-liquid equilibrium: Measuring and modeling with ePC-SAFT and COSMO-RS, Fluid Phase Equilibria 416 (2016) 83-93.
[35] C. Kalidas, G. Hefter, Y. Marcus, Gibbs energies of transfer of cations from water to mixed aqueous organic solvents, Chemical reviews 100 (2000) 819-852.
[36] Y. Marcus, Gibbs energies of transfer of anions from water to mixed aqueous organic solvents, Chemical reviews 107 (2007) 3880-3897.
[37] M.T. Zafarani-Moattar, K. Nasirzade, Osmotic Coefficient of Methanol $+\mathrm{LiCl},+\mathrm{LiBr}$, and +LiCH 3 COO at $25^{\circ} \mathrm{C}$, J. Chem. Eng. Data 43 (1998) 215-219.
[38] S. Han, H. Pan, Thermodynamics of the sodium bromide-methanol-water and sodium bromide-ethanol-water two ternary systems by the measurements of electromotive force at 298. 15K, Fluid Phase Equilibria 83 (1993) 261-270.


[^0]:    * Corresponding author:

    E-mail address: christoph.held@tu-dortmund.de (C. Held).

