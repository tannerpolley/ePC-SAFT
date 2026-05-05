# Predicting Vapor-Liquid Equilibria for Sour-Gas Absorption in Aqueous Mixtures of Chemical and Physical Solvents or Ionic Liquids with ePC-SAFT 

Mark Bülow, Nevin Gerek Ince, Seiya Hirohama, Gabriele Sadowski, and Christoph Held*

Cite This: Ind. Eng. Chem. Res. 2021, 60, 6327-6336
Read Online
Downloaded via BRIGHAM YOUNG UNIV on June 8, 2023 at 01:05:07 (UTC). See https://pubs.acs.org/sharingguidelines for options on how to legitimately share published articles.


#### Abstract

Sour-gas absorption is the main unit operation used in refineries and petrochemical and natural gas processing plants for the effective reduction of climate-wrecking gases, mainly $\mathrm{CO}_{2}$ and $\mathrm{H}_{2} \mathrm{~S}$. Absorption is typically accomplished in an aqueous solvent mixture. The solvent mixture is vastly dependent on the application range; it might contain chemical solvents (amines), activators, and physical solvents. In this work, the vapor-liquid equilibria for absorption of the sour gases $\mathrm{CO}_{2}$ and $\mathrm{H}_{2} \mathrm{~S}$ was investigated in systems containing the chemical solvent methyl diethanolamine (MDEA) and the physical solvents tetrahydrothiophene-1,1-dioxide (sulfolane) or the ionic liquid 1-butyl-3-methylimidazolium acetate. The solubilities of $\mathrm{CO}_{2}$ and $\mathrm{H}_{2} \mathrm{~S}$ were predicted and validated using experimental ![](https://cdn.mathpix.com/cropped/ab722c49-d22f-4bae-b719-9d92d8fdb4fd-01.jpg?height=454&width=809&top_left_y=912&top_left_x=1159) literature data in a broad range of temperature ( $313-373 \mathrm{~K}$ ), sour-gas loading (up to 2 moles gas per moles of MDEA), and pressure (up to 180 bar ) at constant MDEA weight fraction ( $20.9 \mathrm{wt} \%$ ) and sulfolane weight fraction ( $30.5 \mathrm{wt} \%$ ). The equation-of-state electrolyte perturbed-chain statistical associating fluid theory (ePCSAFT) was utilized in this work for the predictions combined with the Born term to physically correctly describe the Gibbs energy of solvation of ions in the aqueous mixture of chemical and physical solvents; this was introduced in a recent work [Bülow, M. et al. Fluid Phase Equilib. 2021, 535, 112967]. Using this approach allowed reducing the total number of binary interaction parameters in these systems of maximum 11 species to a minimum; these parameters were fitted exclusively to data of binary mixtures. The ePCSAFT predictions of the gas solubility were most accurate at low sour-gas loadings and high temperatures. This work provides a thermodynamic framework for the solvent selection for sour-gas absorption in a broad range of conditions. This enables a realistic decrease in experimental effort for solvent selection in sour-gas absorption.


## - INTRODUCTION

Nowadays, further reducing the emission of climate-wrecking gases in power plants and in the (petro-)chemical industry is an urgent demand. State-of-the-art sour-gas absorption units often use solvent blends, i.e., aqueous mixtures containing chemical and physical solvents. Chemical solvents comprise amines, which dissociate in aqueous sour-gas solutions and thus bind the sour gases chemically by amine protonation or carbamate formation. So far, primary, secondary, and tertiary amines have been investigated and industrially used for the absorption process. Absorption kinetics and the amount of absorbed sour gas varies depending on the amine. ${ }^{1}$ The performance of aqueous amine solutions regarding selective absorption can be increased further using mixed amine blends. ${ }^{2}$ Physical solvents often occur in such systems either by constraints of the embedded process or they are added intentionally, e.g., to overcome stoichiometric limitations of the chemical solvent, to increase absorption by physical interaction with the sour gases, and to improve absorption kinetics. For example, at high loadings, $\alpha$, the partial pressure
of $\mathrm{H}_{2} \mathrm{~S}$ in a system including tetrahydrothiophene-1,1-dioxide (hereafter called sulfolane) and the amine methyl diethanolamine (MDEA) is about 10 bar lower compared to systems without sulfolane. That is, physical solvents might strongly influence the overall absorption behavior of the system. A characterization of phase behavior of such systems typically requires expensive lab experiments. A thermodynamic framework capable of predicting the solvent's influence on the system behavior ("solvent screening") is highly desired. In the present work, the solvent effect is studied using sulfolane, which is widely used in industries, as well as ionic liquids (ILs). ILs are widely observed to have a high affinity to sour gases, resulting in superior solubility over other physical solvents. ${ }^{3}$

Received: January 18, 2021
Revised: April 13, 2021
Accepted: April 13, 2021
Published: April 23, 2021
![](https://cdn.mathpix.com/cropped/ab722c49-d22f-4bae-b719-9d92d8fdb4fd-01.jpg?height=216&width=158&top_left_y=2371&top_left_x=1807)

However, the high viscosity of ILs in their pure states is challenging for processing in industrial sour-gas absorption. An application of ILs as an additive to amine-driven sour-gas absorption could benefit the overall solubility. Furthermore, ILs may be specially designed, e.g., with amine-functionalized IL cations ${ }^{4}$ to take over both physical and chemical functions in the absorption process. In the present work, the IL 1-butyl-3-methylimidazolium acetate ( $\left[\mathrm{C}_{4} \mathrm{mim}\right][\mathrm{Ac}]$ ) is used as a model compound to study the influence of ILs on gas solubility. In a solvent screening methodology, a suitable solvent blend is selected by solving a multidimensional optimization problem. This includes, e.g., the choice of the amine and of the physical solvent, their concentration range, temperature, and required pressure to achieve desired gas loadings. This can vary for each application due to the composition of the gas stream that has to be purified. This also includes the presence of a mixture of various sour gases. The sour gases might interact with each other, additionally influencing the respective absorption behavior and the system pressure. The absorption of $\mathrm{CO}_{2}$ gets more complicated in the presence of additional $\mathrm{H}_{2} \mathrm{~S}$, as $\mathrm{H}_{2} \mathrm{~S}$ increases the required pressure to assure a high $\mathrm{CO}_{2}$ solubility. Further, the partial pressure of both sour gases differs largely under defined system conditions. This requires a suitable solvent blend that is advantageous to reach a demanded selectivity of sour gases. To include vapor-liquid equilibria (VLE) predictions for selective sour-gas absorption with a thermodynamic model to further reduce the experimental effort is strongly desired in solvent screening tools.

Different approaches to model VLE for sour-gas absorption have been applied over the last decades. Generally, these frameworks can be divided into activity coefficient models and equations of state (EOS). Activity coefficient models regress the reacting system and have been used early on. ${ }^{5-16}$ Beneficially, these models are easy to use; however, they usually require a huge number of adjustable parameters, i.e., they are highly correlative. In the last decade, EOS have been favored to model phase behavior for absorption processes. EOS have been applied to model multicomponent systems in broad ranges of temperatures and pressures. ${ }^{13,17-27}$ These EOS include nonelectrolyte and electrolyte models, as well as models with or without accounting for association forces due to hydrogen bonding. However, each of such different approaches have investigated very specific systems with a scope on a limited range of conditions, without providing very general and broadly applicable framework.

The present work suggests applying an electrolyte EOS with the focus on a solvent screening: The influence of physical solvents was predicted using a framework based on electrolyte perturbed-chain statistical associating fluid theory (ePCSAFT). The model was used to predict solubility of single sour gases in binary (MDEA) and ternary (MDEA + sulfolane; MDEA $+\left[\mathrm{C}_{4} \mathrm{mim}\right][\mathrm{Ac}]$ ) aqueous blends as well as the selective absorption of one sour gas over a second one in the ternary solvent blend water + MDEA + sulfolane containing both the chemical solvent and physical solvent. All this was studied in broad ranges of temperature, loading, pressure, and concentrations. In previous publications, ePC-SAFT has already been successfully applied to predict VLE for single sour-gas absorption in an aqueous single-amine system and in the presence of methane. ${ }^{25,27}$ In these previous publications, a large number of binary systems have been investigated. Binary interaction parameters were regressed to these systems prior to
the predictions of the actual reacting system. The new approach presented in this work allowed us to drastically reduce the total amount of interaction parameters while maintaining the order of magnitude of the accuracy of the predictions. This was possible by introducing the Born term into ePC-SAFT as suggested in previous works. ${ }^{28,29}$ The Born term accounts for the change in Gibbs energy of solvation for ions upon presence of non-aqueous solvents. In the present work, the original Born term was used that does not account for the ion-concentration dependence of the dielectric constant. This was necessary due to the variety of the ions present in the system under study (e.g., protons, organic ions, inorganic spherical ions, inorganic ions, nonspherical ions that could additionally form associates); it is still unknown how these different ions effect the dielectric constant of the system.

## THEORY

Modeling Reaction Equilibria and Phase Equilibria. Thermodynamic modeling of sour-gas absorption requires to simultaneously solve phase equilibria and reaction equilibria. The underlying reactions taking place in the liquid phase are summarized in eqs 1-5: autoprotolysis of water (1), protonation of the amine MDEA (2), hydration of $\mathrm{CO}_{2}$ (3) and subsequent deprotonation of bicarbonate (4), and deprotonation of $\mathrm{H}_{2} \mathrm{~S}$ (5). Consequently, all reactions are linked to pH in the liquid phase.

$$
\begin{equation*}
\mathrm{H}_{2} \mathrm{O} \rightleftharpoons \mathrm{OH}^{-}+\mathrm{H}^{+} \tag{1}
\end{equation*}
$$

$$
\begin{equation*}
\mathrm{MDEAH}^{+} \rightleftharpoons \mathrm{MDEA}+\mathrm{H}^{+} \tag{2}
\end{equation*}
$$

$$
\begin{equation*}
\mathrm{CO}_{2}+\mathrm{H}_{2} \mathrm{O} \rightleftharpoons \mathrm{HCO}_{3}^{-}+\mathrm{H}^{+} \tag{3}
\end{equation*}
$$

$$
\begin{equation*}
\mathrm{HCO}_{3}^{-} \rightleftharpoons \mathrm{CO}_{3}^{2-}+\mathrm{H}^{+} \tag{4}
\end{equation*}
$$

$$
\begin{equation*}
\mathrm{H}_{2} \mathrm{~S} \rightleftharpoons \mathrm{HS}^{-}+\mathrm{H}^{+} \tag{5}
\end{equation*}
$$

Reaction equilibria are characterized with the activity-based reaction constants $K_{\mathrm{a}} \cdot K_{\mathrm{a}}$ is independent of composition in the liquid phase and of occurring parallel reactions. The dependence of reaction equilibria on concentration or solvents is accounted for by the activity coefficients $\gamma_{i}$ of the reactants and products $i$. That said, $K_{\mathrm{a}}$ is only depending on temperature and pressure; it is the product of the activity coefficient and the mole fraction $x_{i}$ of each reactant and product $i$ (cf., eq 6), which can be expressed as ratios of the product's mole fraction or activity coefficients to the reactant's mole fraction or activity coefficient $\mathrm{K}_{x}$ and $\mathrm{K}_{\gamma}$.

$$
\begin{equation*}
K_{\mathrm{a}}(T, p)=K_{x}(T, p, x) \cdot K_{\gamma}(T, p, x)=\prod_{i}\left(x_{i} \cdot \gamma_{i}\right)^{\nu_{\mathrm{i}}} \tag{6}
\end{equation*}
$$

For the five reactions considered in this work, the $K_{\mathrm{a}}$ values are expressed according to the stoichiometric coefficients $\nu_{\mathrm{i}}$ as follows (eqs 7-11)

$$
\begin{align*}
& K_{\mathrm{a}, 1}=\frac{x_{\mathrm{H}^{+}} \cdot x_{\mathrm{OH}^{-}}}{x_{\mathrm{H}_{2} \mathrm{O}}} \cdot \frac{\gamma_{\mathrm{H}^{+}}^{*} \cdot \gamma_{\mathrm{OH}^{-}}^{*}}{\gamma_{\mathrm{H}_{2} \mathrm{O}}^{0}}  \tag{7}\\
& K_{\mathrm{a}, 2}=\frac{x_{\mathrm{H}^{+}} \cdot x_{\mathrm{MDEA}}}{x_{\mathrm{MDEAH}^{+}}} \cdot \frac{\gamma_{\mathrm{H}^{+}}^{*} \cdot \gamma_{\mathrm{MDEA}}^{*}}{\gamma_{\mathrm{MDEAH}^{+}}^{*}} \tag{8}
\end{align*}
$$

$$
\begin{align*}
& K_{\mathrm{a}, 3}=\frac{x_{\mathrm{H}^{+}} \cdot x_{\mathrm{HCO}_{3}^{-}}}{x_{\mathrm{H}_{2} \mathrm{O}} \cdot x_{\mathrm{CO}_{2}}} \cdot \frac{\gamma_{\mathrm{H}^{+}}^{*} \cdot \gamma_{\mathrm{HCO}_{3}^{-}}^{*}}{\gamma_{\mathrm{H}_{2} \mathrm{O}}^{0} \cdot \gamma_{\mathrm{CO}_{2}}^{*}}  \tag{9}\\
& K_{\mathrm{a}, 4}=\frac{x_{\mathrm{H}^{+}} \cdot x_{\mathrm{CO}_{3}^{2-}}}{x_{\mathrm{HCO}_{3}^{-}}} \cdot \frac{\gamma_{\mathrm{H}^{+}}^{*} \cdot \gamma_{\mathrm{CO}_{3}^{2-}}^{*}}{\gamma_{\mathrm{HCO}_{3}^{-}}^{*}}  \tag{10}\\
& K_{\mathrm{a}, \mathrm{~S}}=\frac{x_{\mathrm{H}^{+}} \cdot x_{\mathrm{HS}^{-}}}{x_{\mathrm{H}_{2} \mathrm{~S}}} \cdot \frac{\gamma_{\mathrm{H}^{+}}^{*} \cdot \gamma_{\mathrm{HS}^{-}}^{*}}{\gamma_{\mathrm{H}_{2} \mathrm{~S}}^{*}} \tag{11}
\end{align*}
$$

Two different reference states were used in the literature to establish correlations for these $K_{\mathrm{a}}$ values: the pure-component reference state and the reference state of an ideally diluted solution. The generic activity coefficient $\gamma_{i}^{0}$ (eq 12) refers to the pure component and is only used for water. In contrast, the rational activity coefficient $\gamma_{i}^{*}$ refers to the infinite dilution state in pure water. This is important to mention since it is exactly the way the $K_{\mathrm{a}}$ values were determined in the literature. If any other reference states are chosen (which is not forbidden), these must be converted accordingly. The ideal dilution in pure water was applied to all charged species, sour gases, and even the chemical solvent MDEA (eq 13). The latter is often misused in the literature. Activity coefficients can be obtained from EOS by a ratio of fugacity coefficients $\varphi_{\mathrm{i}}$ of the component $i$ at the considered composition of the liquid phase over the composition in the reference state

$$
\begin{align*}
& \gamma_{i}^{0}=\frac{\varphi_{i}(T, p, x)}{\varphi_{0 i}\left(T, p, x_{i}=1\right)}  \tag{12}\\
& \gamma_{i}^{*}=\frac{\varphi_{i}(T, p, x)}{\varphi_{i}^{\infty}\left(T, p, x_{i}=0\right)} \tag{13}
\end{align*}
$$

A thermodynamic framework for calculating the reaction equilibria requires reaction constants. These values are readily available in the literature (cf. Table S1 in the SI). To cover broad temperature ranges, the $K_{a}$ values are expressed as a function of temperature (eq 14).

$$
\begin{equation*}
\ln K_{\mathrm{a}}(T)=c_{1}+\frac{c_{2}}{T}+c_{3} \cdot \ln T \tag{14}
\end{equation*}
$$

Because of discussions on the validity of such constants, it should be mentioned here that the $K_{a}$ value for $\mathrm{MDEAH}^{+}$ dissociation was obtained in the literature from titration measurements and extrapolation of experimental data to infinite dilution. Thus, a thermodynamic model was not used to determine $K_{a, 2}$, and further, the reference state "ideally diluted solution" must be used to convert $K_{a}$ values to $K_{x}$ values (eq 6) under defined reaction conditions.

Throughout the ongoing calculations of the underlying reactions in the liquid phase, the mole balances for each component are accounted for by eq 15. The balance is depending on the reaction coordination number $\lambda_{r}$ for each reaction $r$. Certainly, electroneutrality is respected.

$$
\begin{equation*}
n_{i}=n_{i, 0}+\sum_{r=1}^{5} v_{i} \cdot \lambda_{r} \tag{15}
\end{equation*}
$$

The thermodynamic representation of the absorption of sour gases into aqueous solvent blends is realized by vapor-liquid equilibrium conditions. In the present work, VLE is modeled based on the isofugacity criterion, which must be fulfilled for
each component (eq 16) that is present in the two phases $L$ and $V$. In this work, only the components $i=$ water, MDEA, sulfolane, $\mathrm{CO}_{2}$, and $\mathrm{H}_{2} \mathrm{~S}$ are assumed to be volatile, i.e., only these are distributed between $L$ and $V$.

$$
\begin{equation*}
\varphi_{i}^{L} \cdot x_{i}=\varphi_{i}^{V} \cdot y_{i} \tag{16}
\end{equation*}
$$

Here, $y_{i}$ is the mole fraction of the component $i$ in the vapor phase. By definition within this work, charged species have a negligible vapor pressure and are only present in the liquid phase, and eq 16 does not need to be solved for charged species. However, the fugacity coefficients $\varphi_{i}^{L}$ in the liquid phase are influenced by the electrolyte species, which is explicitly accounted for in this work.
ePC-SAFT Framework. In this work, the activity coefficients were calculated with ePC-SAFT. ePC-SAFT accounts for four independent contributions to the residual Helmholtz energy $a^{\text {res }}$ (eq 17).

$$
\begin{equation*}
a^{\mathrm{res}}=a^{\mathrm{hc}}+a^{\mathrm{disp}}+a^{\mathrm{assoc}}+a^{\mathrm{ion}}+a^{\mathrm{born}} \tag{17}
\end{equation*}
$$

The independent contributions are $a^{\text {hc }}$ (hard-chain contribution), $a^{\text {disp }}$ (dispersion contribution), and $a^{\text {assoc }}$ (association contribution); these are available through abstracting the gases, solvents, reactants, and products with pure-component parameters. These are the segment number $m_{\text {seg }}$, segment diameter $\sigma_{i}$, dispersion-energy parameter $u_{i} / k_{\mathrm{B}}$, and two association parameters $\varepsilon^{\mathrm{A}_{i} \mathrm{~B}_{i}}$ (energy parameter) and $\kappa^{\mathrm{A}_{i} \mathrm{~B}_{i}}$ (volume parameter). For associating components, additionally, the association scheme is needed, depicting proton donators and proton acceptors, e.g., 2 B model for water. ${ }^{30}$ The charge of electrolyte species is accounted for by the Debye-Hückel contribution $a^{\text {ion }}$ (eq 18).

$$
\begin{equation*}
a^{\text {ion }}=-\frac{\kappa}{12 \pi \cdot \varepsilon} \sum_{j} \chi_{j} q_{j}^{2} x_{j} \tag{18}
\end{equation*}
$$

with the inverse Debye length $\kappa$, the auxiliary function $\chi_{j}$, and the charge $q_{j}$ (detailed in ref 31.). The dielectric constant $\varepsilon$ is the product of the vacuum dielectric constant $\varepsilon_{0}$ and the dielectric constant of the solvent $\varepsilon_{r}\left(\varepsilon=\varepsilon_{0} \cdot \varepsilon_{r}\right)$. Further, charged species are assumed to be nonassociative, and dispersive cation-cation and anion-anion interactions are neglected.

The dielectric constant $\varepsilon_{r}$ is an important variable in the Debye-Hückel theory and in the Born term, as it screens the interaction between a charged species and the surrounding medium. This indicates how strong the charged species is shielded by the solvent from interactions with other (charged) species. $\varepsilon_{r}$ depends on the solvent composition and on temperature, changing from about 78 for water to 44 for sulfolane and 22.4 for MDEA at $298 \mathrm{~K} .^{32}$ To account for possible variation in solvent blends of ternary and higher orders, $\varepsilon_{r}$ is derived by summing up the solvents' purecomponent constants $\varepsilon_{r, i}$ multiplied by the respective solvent mole fraction $x_{i}$ in the mixture (eq 19); such a mixing rule is required in consequence of a lack in experimental data. As it is still unknown how different ions effect the dielectric constant of the system, the ions are in this work neglected in the mixing rule; $\varepsilon_{r}$ depends on the solvents.

$$
\begin{equation*}
\varepsilon_{\mathrm{r}}(\bar{x}, T)=\sum_{i=\text { solvent }} \varepsilon_{\mathrm{r}, i}(T) \cdot x_{i} \tag{19}
\end{equation*}
$$

The Born term ${ }^{33}$ developed in 1920 is used to calculate the Gibbs energy of solvation needed to transfer the respective ion
from vacuum into the solvent (or the solvent blend in this work). The Born contribution to the residual Helmholtz energy is calculated via eq 20

$$
\begin{equation*}
a^{\text {born }}=-\frac{e^{2}}{4 \pi \varepsilon_{0} k T}\left(1-\frac{1}{\varepsilon_{r}}\right) \sum_{i} \frac{x_{i} z_{i}^{2}}{a_{i}} \tag{20}
\end{equation*}
$$

where $\alpha_{i}$ is the same diameter used in the Debye-Hückel theory that is similar to the ePC-SAFT segment diameter $\sigma_{i}$.

The fugacity coefficient is calculated by differentiating $a^{\text {res }}$ with respect to density and with respect to mole fraction (chemical potential $\mu_{i}^{\text {res }}$ ).

$$
\begin{equation*}
\ln \left(\varphi_{i}\right)=\frac{\mu_{\mathrm{i}}^{\text {res }}}{k_{\mathrm{B}} \cdot T}-\ln \left(1+\left(\frac{\partial\left(\frac{a^{\text {res }}}{k_{\mathrm{B}} \cdot T}\right)}{\partial \rho}\right)\right) \tag{21}
\end{equation*}
$$

VLE for sour-gas absorption needs prediction of properties in mixtures. Therefore, Berthelot-Lorenz ${ }^{34,35}$ and WolbachSandler ${ }^{36}$ combining and mixing rules are incorporated into ePC-SAFT (eqs 22-25). This is the conventional procedure in PC-SAFT modeling. ${ }^{38}$

$$
\begin{align*}
& \sigma_{i j}=\frac{1}{2}\left(\sigma_{i}+\sigma_{j}\right)  \tag{22}\\
& \varepsilon^{A_{i} B_{j}}=\frac{1}{2} \cdot\left(\varepsilon^{A_{i} B_{i}}+\varepsilon^{A_{j} B_{j}}\right)  \tag{23}\\
& u_{i j}=\sqrt{u_{i} u_{j}}\left(1-k_{i j}(T)\right)  \tag{24}\\
& \kappa^{A_{i} B_{j}}=\sqrt{\kappa^{A_{i} B_{i}} \cdot \kappa^{A_{j} B_{j}}}\left(\frac{\sqrt{\sigma_{i} \cdot \sigma_{j}}}{\sigma_{i j}}\right)^{3} \tag{25}
\end{align*}
$$

Eq 24 introduces the binary interaction parameter $k_{i j}$; see eq 26.

$$
\begin{align*}
k_{i j}(T)= & k_{i j, a}+k_{i j, T} \cdot\left(T-T^{\mathrm{ref}}\right)+k_{i j, T^{2}} \cdot\left(T-T^{\mathrm{ref}}\right)^{2} \\
& +k_{i j, T^{3}} \cdot\left(T-T^{\mathrm{ref}}\right)^{3} \tag{26}
\end{align*}
$$

For the temperature dependency, $T^{\text {ref }}=298.15 \mathrm{~K}$ denotes an arbitrarily chosen reference temperature. This procedure was introduced by Wangler et al. ${ }^{27}$ to describe the binary system water + MDEA as accurately as possible. A flowchart for the VLE predictions within the ePC-SAFT framework is shown in Figure S3 in the SI.

Deriving Pure-Component Parameters. ePC-SAFT has already been applied to water-MDEA-sour gas systems, and the majority of the pure-component parameters are available in the literature (cf. Table 1). In this work, the physical solvent sulfolane was additionally considered and required purecomponent parameters. The pure-component parameters for sulfolane were fitted in this work to pure-component vapor pressures and liquid densities from the literature (see Table S2 and Figures S1 and S2). Pure-component parameters were obtained by minimizing the following objective function (OF) (eq 27)

Table 1. ePC-SAFT Pure-Component Parameters for the Solvents and Sour Gases ${ }^{\boldsymbol{a}}$
| parameter | $\mathrm{H}_{2} \mathrm{O}$ | MDEA | sulfolane | $\mathrm{CO}_{2}$ | $\mathbf{H}_{\mathbf{2}} \mathbf{S}$ |
| :--- | :--- | :--- | :--- | :--- | :--- |
| $m_{\mathrm{i}}^{\text {seg }}$ | 1.2046 | 3.6750 | 2.7894 | 2.0729 | 1.6686 |
| $\sigma / \AA$ | $b$ | 3.5630 | 3.7020 | 2.7852 | 3.0349 |
| $u_{i} / k_{\mathrm{B}} / \mathrm{K}$ | 353.94 | 228.71 | 436.47 | 169.21 | 229.00 |
| $N^{\text {assoc }}$ | 1:1 | 2:2 | 1:1 | 1:1 |  |
| $\frac{\varepsilon^{\mathrm{A}_{\mathrm{i}} \mathrm{B}_{\mathrm{i}}}}{k_{B}} / \mathrm{K}$ | 2425.6 | 2046.6 | 0 | 0 |  |
| $\kappa^{\mathrm{A}_{\mathrm{i}} \mathrm{B}_{\mathrm{i}}}$ | 0.045 | 0.1238 | 0.045 | 0.045 |  |
| Ref | 30 | 24 | This work | 38 | 23 |


${ }^{a} N^{\text {assoc }}$ gives the number of associating sites. ${ }^{b} \sigma=2.7927+(10.11 \cdot \mathrm{e}^{-0.01775 \mathrm{~T}}-1.417 \cdot \mathrm{e}^{-0.01146 \mathrm{~T}}$ ).

$$
\begin{align*}
& \mathrm{OF}=\sum_{k}^{\mathrm{NP}\left(\rho_{0 i}\right)}\left|1-\left(\frac{\rho_{0 i}^{\text {pred }}}{\rho_{0 i}^{\exp }}\right)\right|+\sum_{m}^{\mathrm{NP}\left(p_{0 i}^{\mathrm{LV}}\right)}\left|1-\left(\frac{p_{0 i}^{\mathrm{LV}, \text { pred }}}{p_{0 i}^{\mathrm{LV}, \exp }}\right)\right| \\
& =\min ! \tag{27}
\end{align*}
$$

where $\rho_{0 i}^{\text {pred }}$ and $\rho_{0 i}^{\exp }$ are the predicted and experimental densities of pure sulfolane, respectively. Equally, $\rho_{0 i}^{\mathrm{LV} \text {,pred }}$ and $\rho_{0 i}^{\mathrm{LV}, \text { exp }}$ are predicted and experimental pure-component vapor pressures, respectively. NP is the number of available experimental data points used in the fitting process.

## - RESULTS AND DISCUSSION

ePC-SAFT predictions and modeling results, respectively, are compared to experimental data in the literature. For easy comparison, the absolute relative deviation (ARD\%) calculated as in eq 28 is applied.

$$
\begin{equation*}
\mathrm{ARD} \%=\frac{100 \%}{\mathrm{NP}} \cdot \sum_{i}^{\mathrm{NP}}\left|1-\frac{p_{i}^{\mathrm{pred}}}{p_{i}^{\exp }}\right| \tag{28}
\end{equation*}
$$

where $p_{i}^{\text {pred }}$ and $p_{i}^{\text {exp }}$ denote ePC-SAFT results and experimental pressures, respectively; "pred" denotes results that were achieved without fitting a parameter to the mixture that is under consideration. Modeling means that one or more parameters were fitted to the mixture under study.
ePC-SAFT Pure-Component Parameters and Binary Parameters. The pure-component parameters for all components involved in the VLE modeling for sour-gas absorption except sulfolane are already available from previous works. Pure-component parameters for the molecular compounds water, MDEA, and the sour gases $\mathrm{CO}_{2}$ and $\mathrm{H}_{2} \mathrm{~S}$ are listed in Table 1. In this work, $\mathrm{CO}_{2}$ was modeled with induced association toward water. Experimental vapor pressure for sulfolane is much scattered in the available temperature range. Therefore, vapor pressures data calculated from Antoine parameters (Table S2 in the SI) was used in the parameter estimation, and the results were then evaluated against the available experimental data. Liquid densities are available in the literature. ${ }^{37}$ Both data were used for fitting. Sulfolane was assumed to be a nonassociating component for the parameter regression. Similar to $\mathrm{CO}_{2}$, sulfolane is modeled with an induced association toward water. The derived purecomponent parameters are listed in Table 1. The comparisons of experimental and fitted liquid densities are depicted in Figure S1, with an ARD of $0.068 \%$. Vapor-pressure fitting is shown in Figure S2 in the SI (ARD $=0.598 \%$ ). Although only
regressed in the temperature range above 382 K (valid range for the employed Antoine parameters), ePC-SAFT reliably predicts the vapor pressure at lower temperatures; these temperatures are especially important for the prediction of the reacting systems for sour-gas absorption.

The products of the reactions taking place in the liquid phase are charged species. Pure-component parameters for these species have also been determined in previous publications and are listed in Table 2.

Table 2. ePC-SAFT Pure-Component Parameters for the Charged Species Including the Ions of the IL
| ion | $m_{i}^{\text {seg }}$ | $\sigma / \AA$ | $u_{i} / k_{\mathrm{B}} / \mathrm{K}$ | $z_{i}$ | ref |
| :--- | :--- | :--- | :--- | :--- | :--- |
| $\mathrm{H}^{+}$ | 1 | 3.4654 | 500.00 | +1 | 26 |
| $\mathrm{OH}^{-}$ | 1 | 2.0177 | 650.00 | -1 | 26 |
| $\mathrm{HCO}_{3}{ }^{-}$ | 1 | 2.9296 | 70.00 | -1 | 26 |
| $\mathrm{CO}_{3}{ }^{2-}$ | 1 | 2.4422 | 249.26 | -2 | 26 |
| $\mathrm{HS}^{-}$ | 1 | 3.0349 | 229.00 | -1 | 24 |
| $\mathrm{MDEAH}^{+}$ | 1 | 3.5630 | 228.71 | +1 | 22 |
| $\left[\mathrm{C}_{4} \mathrm{mim}\right]^{+}$ | 2.4805 | 3.6371 | 218.14 | +1 | 39 |
| $[\mathrm{AC}]^{-}$ | 3.7266 | 3.5605 | 533.11 | -1 | 40 |


The binary interaction parameters regressed in earlier publications on sour-gas absorption ${ }^{25,27}$ for molecular compounds and ion-solvent interactions are not applied in this work. In fact, due to the inclusion of the Born term, the number of binary interaction parameters could be significantly reduced, emphasizing the predictive character of ePC-SAFT. The sour gas $\mathrm{CO}_{2}$ was modeled with induced association, and new binary interaction parameters have been regressed to the respective VLE data with water. The updated binary interaction parameters are listed in Table 3.

For the system sulfolane $+\mathrm{H}_{2} \mathrm{~S}$, no fitting was needed. Data of the binary system sulfolane + MDEA were not available in the literature, and $k_{i j}$ was set to zero. Additionally, no interaction parameters were fitted between sulfolane and charged species. Further, VLE modeling of sulfolane with water showed that induced association covers the experimental data best applying the " 2 B " association model according to Wertheim.

Modeling the VLE of Sulfolane + Sour Gas and Sulfolane + Water. The VLEs of the binary system sulfolane + sour gas are depicted in Figure 1 for $\mathrm{CO}_{2}$ and in Figure 2 for $\mathrm{H}_{2} \mathrm{~S}$. For the VLE with $\mathrm{CO}_{2}$, a temperature-independent binary interaction parameter was adjusted to $k_{i j}=0.055$, while the system sulfolane $+\mathrm{H}_{2} \mathrm{~S}$ was accurately predicted ( $k_{i j}=0$ ) with ePC-SAFT. The modeled VLE was compared to experimental data from Jou et al. ${ }^{41}$ Isothermal data were

![](https://cdn.mathpix.com/cropped/ab722c49-d22f-4bae-b719-9d92d8fdb4fd-05.jpg?height=323&width=430&top_left_y=235&top_left_x=1327)
Figure 1. $P$, $x$-diagram of the binary system sulfolane $+\mathrm{CO}_{2}$ at four temperatures varying from 298 to 403 K . Lines are ePC-SAFT correlations; symbols are experimental data from Jou et al. ${ }^{41} k_{i j}=$ 0.055 .

![](https://cdn.mathpix.com/cropped/ab722c49-d22f-4bae-b719-9d92d8fdb4fd-05.jpg?height=327&width=430&top_left_y=759&top_left_x=1327)
Figure 2. $P$, $x$-diagram of the binary system sulfolane $+\mathrm{H}_{2} \mathrm{~S}$ at four temperatures varying from 298 to 403 K . Lines are ePC-SAFT predictions; symbols are experimental data from Jou et al. ${ }^{41}$

available for temperatures varying from 298 to 403 K for both systems sulfolane $+\mathrm{CO}_{2}$ and sulfolane $+\mathrm{H}_{2} \mathrm{~S}$.

The binary VLE of the system sulfolane + water was modeled with ePC-SAFT in the temperature range of 298-323 K (cf. Figure 3). For this purpose, the induced-association

![](https://cdn.mathpix.com/cropped/ab722c49-d22f-4bae-b719-9d92d8fdb4fd-05.jpg?height=330&width=430&top_left_y=1480&top_left_x=1327)
Figure 3. $P$, $x$-diagram of the binary system sulfolane + water at five temperatures varying from 298 to 323 K . Lines are ePC-SAFT modeling results with induced association between sulfolane and water; symbols are experimental data from Tommila et al. ${ }^{43} k_{i j}=$ 0.015 .

approach proposed by Kleiner et al. ${ }^{42}$ with water was found to best represent the experimental data. Induced association especially enhanced the prediction accuracy of the system

Table 3. Parameters for the Calculation of $k_{i j}$ via eq $26{ }^{a}$
| pairs | $k_{i j, a}$ | $k_{i j, T}$ | $k_{i j, T^{2}}$ | $k_{i j, T^{3}}$ | ref |
| :--- | :--- | :--- | :--- | :--- | :--- |
| water-MDEA | -0.118 | $-7.5910^{-4}$ | $-1.3210^{-5}$ | $1.2510^{-7}$ | this work |
| water- $\mathrm{CO}_{2}$ | $-2.110^{-3}$ | $4.5310^{-4}$ | 0 | 0 | this work |
| water- $\mathrm{H}^{+}$ | 0.25 | 0 | 0 | 0 | 31 |
| water $-\mathrm{OH}^{-}$ | -0.25 | 0 | 0 | 0 | 31 |
| water- $\mathrm{CO}_{3}^{2-}$ | -0.25 | 0 | 0 | 0 | 31 |
| sulfolane-water | 0.015 | 0 | 0 | 0 | this work |
| sulfolane- $\mathrm{CO}_{2}$ | 0.055 | 0 | 0 | 0 | this work |


[^0]pressure in regions of high sulfolane mole fractions. For an accurate modeling of the VLE, a temperature-independent $k_{i j}=$ 0.015 was adjusted.

## Predicting VLE for Sour-Gas Absorption in Aqueous

MDEA Blends. The VLE for absorption of the sour gases $\mathrm{CO}_{2}$ and $\mathrm{H}_{2} \mathrm{~S}$ into aqueous MDEA blends has been predicted with ePC-SAFT in two previous papers. ${ }^{25,27}$ In the following, these predictions were performed again including the Born term within ePC-SAFT prior to investigating the influence of sulfolane on the VLE. For clarity, the systems under study are of similar conditions with respect to the weight fraction of MDEA: $19 \mathrm{wt} \%$ in the aqueous MDEA system and $20.5 \mathrm{wt} \%$ in the system with MDEA and sulfolane. Figure 4 presents the

![](https://cdn.mathpix.com/cropped/ab722c49-d22f-4bae-b719-9d92d8fdb4fd-06.jpg?height=314&width=430&top_left_y=772&top_left_x=410)
Figure 4. Partial pressure of $\mathrm{CO}_{2}$ vs sour-gas loading in the system water + MDEA $(19 \mathrm{wt} \%)+\mathrm{CO}_{2}$. Lines are predictions with ePCSAFT using parameters from Tables 1-3; symbols are experimental data (points: 313 K , squares: 333 K , triangles: 373 K , stars: 393 K ). ${ }^{10}$ ARD\% are listed in Table S3 in the SI.

predicted system pressures for $\mathrm{CO}_{2}$ absorption in systems consisting of $19 \mathrm{wt} \%$ MDEA. Results are evaluated in terms of ARD\% in Table S3 in the SI. The results for higher temperatures ( 373 and 393 K ) are more accurate than for 313 and 333 K . The pressure at lower loadings is underestimated at higher temperatures, while it is overestimated at low temperatures. Bearing in mind the drastic reduction in the number of binary interaction parameters regressed to binary systems, the results are in good to very good agreement with the experimental data. At lower temperatures ( 313 and 333 K ), a distinct underestimation of the partial pressure of $\mathrm{CO}_{2}$ can be observed for high loadings ( $\alpha>1$ ). That is, ePC-SAFT overestimates the extent of the reactions to $\mathrm{HCO}_{3}{ }^{-}$and $\mathrm{CO}_{3}{ }^{2-}$ in the liquid phase, i.e., the solubility of $\mathrm{CO}_{2}$ is overestimated.

The results for the system water + MDEA $(19 \mathrm{wt} \%)+\mathrm{H}_{2} \mathrm{~S}$ are depicted in Figure 5 for temperatures varying between 313 and 393 K . At high loadings, the system pressure is underestimated.

![](https://cdn.mathpix.com/cropped/ab722c49-d22f-4bae-b719-9d92d8fdb4fd-06.jpg?height=318&width=430&top_left_y=2109&top_left_x=410)
Figure 5. Partial pressure of $\mathrm{H}_{2} \mathrm{~S}$ vs sour-gas loading in the system water + MDEA $(19 \mathrm{wt} \%)+\mathrm{H}_{2} \mathrm{~S}$. Lines are predictions with ePCSAFT using parameters from Tables 1-3; symbols are experimental data (points: 313 K , squares: 333 K , triangles: 373 K , stars: 393 K ). ${ }^{10}$ ARD\% are listed in Table S4 in the SI.

Predicting VLE for Sour-Gas Absorption in Aqueous MDEA Blends in the Presence of Sulfolane. VLE for sourgas absorption into aqueous MDEA solutions, also in the presence of an inert gas, has already been successfully predicted with ePC-SAFT in the literature. ${ }^{22,24}$ Results including the Born term are depicted in the previous chapter. The focus now is to account for the addition of sulfolane as an example for a physical solvent to predict VLE to absorb a single sour gas, $\mathrm{CO}_{2}$ or $\mathrm{H}_{2} \mathrm{~S}$. Interaction parameters have not been fitted to data of the reacting system, and all following results with ePC-SAFT are thus predictions. Experimental data for systems including sulfolane are rather scarce, and only a few references were found that contains VLE data of selective $\mathrm{CO}_{2}$ and $\mathrm{H}_{2} \mathrm{~S}$ absorption in MDEA solution containing sulfolane, e.g., in ref 34. These experimental data cover temperatures of 313 and 373 K with $20.9 \mathrm{wt} \%$ MDEA and $30.5 \mathrm{wt} \%$ sulfolane. These $\mathrm{wt} \%$ values are related to constant composition of the aqueous solvent blend.

The VLE of the reacting system water + MDEA + sulfolane $+\mathrm{H}_{2} \mathrm{~S}$ is illustrated in Figure 6. Predictions using ePC-SAFT

![](https://cdn.mathpix.com/cropped/ab722c49-d22f-4bae-b719-9d92d8fdb4fd-06.jpg?height=317&width=430&top_left_y=1036&top_left_x=1327)
Figure 6. Partial pressure of $\mathrm{H}_{2} \mathrm{~S}$ vs sour-gas loading in the system water + MDEA $(20.9 \mathrm{wt} \%)+$ sulfolane $(30.5 \mathrm{wt} \%)+\mathrm{H}_{2} \mathrm{~S}$. Lines are predictions with ePC-SAFT using parameters from Tables 1-3; symbols are experimental data (stars: 373 K , triangles: 313 K ). ${ }^{41} \operatorname{ARD}(373 \mathrm{~K})=27.26 \% ; \operatorname{ARD}(313 \mathrm{~K})=13.15 \%$.

are based on the pure-component parameters and binary interaction parameters listed in Table 1 to Table 3. For isothermal conditions at 313 K , the partial pressure of $\mathrm{H}_{2} \mathrm{~S}$ is underestimated compared to the experimental data with an ARD of $27.26 \%$. Nevertheless, the exponential dependency of the partial pressure with increased loadings is predicted qualitatively correct. At high loadings, the partial pressure of $\mathrm{H}_{2} \mathrm{~S}$ is slightly underestimated. This is a very satisfying result considering that no binary interaction parameters between sulfolane and $\mathrm{H}_{2} \mathrm{~S}$ were used. At 313 K , the prediction of the required partial pressure is in almost quantitative agreement with the experimental data, again with a small deviation at higher loadings, overall leading to ARD of only $13.15 \%$. At loadings $\alpha<1$, the predictions give ARD $<10 \%$. Similar to the results for the systems without sulfolane, ePC-SAFT gives better representation at higher temperatures.

The ePC-SAFT-predicted VLE of the reacting system water $+\mathrm{MDEA}+$ sulfolane $+\mathrm{CO}_{2}$ is also in good agreement with the experimental data at 313 and 373 K (Figure 7). The strong scattering of the experimental data at 373 K is included in the ARD. Due to the scattering, ePC-SAFT predictions are slightly deviating compared to the results for the $\mathrm{H}_{2} \mathrm{~S}$ systems with $\mathrm{ARD}=36.50 \%$. The partial pressure of $\mathrm{CO}_{2}$ is marginally underestimated over the whole loading range at 373 K . That is, the influence of sulfolane on the chemical reactions, the hydration of $\mathrm{CO}_{2}$, and deprotonation of the carbonate (eqs 3 and 4) are overestimated. For isothermal conditions at 313 K ,

![](https://cdn.mathpix.com/cropped/ab722c49-d22f-4bae-b719-9d92d8fdb4fd-07.jpg?height=312&width=430&top_left_y=235&top_left_x=412)
Figure 7. Partial pressure of $\mathrm{CO}_{2}$ vs sour-gas loading in the system water + MDEA ( $20.9 \mathrm{wt} \%$ ) + sulfolane ( $30.5 \mathrm{wt} \%$ ) + $\mathrm{CO}_{2}$. Lines are predictions with ePC-SAFT using parameters from Tables 1-3; symbols are experimental data (stars: 373 K , triangles: 313 K ). ${ }^{41}$ ARD $\%(373 \mathrm{~K})=36.50 \%$; ARD $\%(313 \mathrm{~K})=35.62 \%$.

the same overall effect of sulfolane results in an underestimation of the partial pressure in comparison to experimental pressures for $\mathrm{CO}_{2}$. Still the prediction shows satisfying agreement ( $\mathrm{ARD}=35.62 \%$ ). The small deviation can easily be overcome by introducing binary interaction parameters for sulfolane with the carbonates to reduce the influence of sulfolane on the reaction equilibrium. This is not the focus of this work.

Especially at lower temperatures, deviations to the experimental data become visible. More data on the solvent system water + MDEA + sulfolane are available for the sour gas $\mathrm{CO}_{2}$. Data with three compositions of the solvent blend are available varying from 30 to $40 \mathrm{wt} \%$ for MDEA and from 10 to $20 \mathrm{wt} \%$ for sulfolane. To understand the interactions at lower temperatures (the largest deviations for the systems without the physical solvent), Figure 8 shows predictive results for ePC-SAFT at 328.15 K , the lowest investigated temperature.

The overall results are satisfying with ARD\% values not higher than $28 \%$. Still, there is a pattern visible for loadings $\alpha$ near and above 1 mol of $\mathrm{CO}_{2}$ per mol of MDEA. This pattern was already visible in the binary solvent blend where a change in the predictions in comparison to the experimental data is visible. Near or above $\alpha=1$, the pressure increases immediately, indicating the limitation of the chemisorption to some extent. The model here shows a slower increase in pressure. Following the whole concentration range, the model performs almost quantitatively up to $\alpha=1$. Still, the unsatisfying deviation for higher loadings needs a careful investigation. It might be reasoned in the induced association of $\mathrm{CO}_{2}$ but given the fact that the model's performance for $\mathrm{H}_{2} \mathrm{~S}$ is comparable, other influences might be more important. At higher loadings, most of the systems' components are of ionic character. Representing these species is very important for a correct prediction. In the future, an extended version of the model ${ }^{29}$ that includes an altered Born contribution for a
composition-dependent approach also including ionic species might be applied.

Predicting VLE for Sour-Gas Absorption in Aqueous MDEA Blends in the Presence of the lonic Liquid [C4mim][Ac]. Promising experimental results for sour-gas solubility in pure ILs have suggested their usage as physical solvents in reactive systems as substitutes to classical physical solvents. The screening of the huge amount of possible ILs is only meaningful with predictive tools such as ePC-SAFT. Basically, ILs might be modeled using two approaches. One which accounts for the ionic character and utilizes ion-specific parameters, and one that, similar to conventional solvents, is based on a molecular parameterization. In a recent publication, the molecular and ion-specific approaches have been compared with emphasis on the prediction of sour-gas absorption into the pure IL. ${ }^{44}$ Although utilizing the extremely low vapor pressure of the IL in the parameter estimation in the molecular approach, the ion-specific parameter set was superior at higher pressures. Therefore, an ion-specific approach is also chosen here to model the IL. The effect of the IL $\left[\mathrm{C}_{4} \operatorname{mim}\right][\mathrm{Ac}]$ on VLE for sour-gas absorption is predicted with ePC-SAFT including the Born term as depicted in Figure 9. These results are compared in the same figure to the $\mathrm{CO}_{2}$ absorption in the aqueous MDEA blend ( $30 \mathrm{wt} \%$ ) without IL. The latter system can be predicted with ARD $=32.68 \%$, caused by the high deviations at $\alpha>1$. This is attributed to the overestimation of reaction eqs 3 and/or 4 that causes the concentration of $\mathrm{CO}_{2}$ to be underestimated and therefore reduces the partial pressure of $\mathrm{CO}_{2}$ stronger than expected. Although this problem was addressed by the introduction of the Born term (hindering the formation of ions due to the increased solvation energies), the effect is still pronounced at high loadings, especially at low temperatures. The underestimation of the pressure at higher loadings is transferred also to the system including the IL, giving similar deviations for the overall systems (ARD $= 32.44 \%$ ). Reducing the error in the IL-free MDEA system thus ultimately reduces the deviation also introduced to the system including the IL. That said, the overall effect of the IL on the sour gas partial pressure is covered reasonably at loadings up to 0.9.

Predicting VLE for Selective Sour-Gas Absorption in Aqueous MDEA Blends in the Presence of Sulfolane. Real processes involve increasingly more presence of solvent blends as well as a mixture of at least two sour gases that mutually influences the absorption behavior. The sour gases $\mathrm{CO}_{2}$ and $\mathrm{H}_{2} \mathrm{~S}$ are chemically absorbed into the liquid phase yielding protons available for MDEA protonation. If two such gases are simultaneously involved in the absorption process, the reactions share protons (or $\mathrm{MDEAH}^{+}$) as a common

![](https://cdn.mathpix.com/cropped/ab722c49-d22f-4bae-b719-9d92d8fdb4fd-07.jpg?height=339&width=1477&top_left_y=2162&top_left_x=346)
Figure 8. Partial pressure of $\mathrm{CO}_{2}$ vs sour-gas loading in the system water $+\mathrm{MDEA}+$ sulfolane $+\mathrm{CO}_{2}$ at $T=328 \mathrm{~K}$. Lines are predictions with ePCSAFT using parameters from Tables $1-3$; symbols are experimental data: (a) MDEA $40 \mathrm{wt} \%$, sulfolane $10 \mathrm{wt} \%$, (b) MDEA $30 \mathrm{wt} \%$, sulfolane 20 $\mathrm{wt} \%$, (c) MDEA $40 \mathrm{wt} \%$, sulfolane $20 \mathrm{wt} \%$. ARD\%: $(\mathrm{a})=26.38 \%$, $(\mathrm{b})=27.34 \%$, $(\mathrm{c})=19.26 \%$.

![](https://cdn.mathpix.com/cropped/ab722c49-d22f-4bae-b719-9d92d8fdb4fd-08.jpg?height=327&width=430&top_left_y=235&top_left_x=412)
Figure 9. Partial pressure of $\mathrm{CO}_{2}$ vs sour-gas loading in the system water $+\operatorname{MDEA}(30 \mathrm{wt} \%)+\left[\mathrm{C}_{4} \operatorname{mim}\right][\mathrm{Ac}](10 \mathrm{wt} \%)+\mathrm{CO}_{2}$ (stars) and in the system water + MDEA (30 wt \%) for comparison (triangles). Lines are predictions with ePC-SAFT using parameters from Tables 1-3; symbols are experimental data. ${ }^{45}$ ARD\% (MDEA + $\mathrm{IL})=32.44 \%$; ARD\% (MDEA) $=32.68 \%$.

product. The importance of accounting for the interactions with a versatile thermodynamic model to represent a realistic scenario is imminent. Recently, Cleeton et al. ${ }^{20}$ presented ePCSAFT VLE predictions for competitive $\mathrm{CO}_{2}$ and $\mathrm{H}_{2} \mathrm{~S}$ absorption in an aqueous MDEA blend. In this work, ePCSAFT predictions were performed for the first time toward VLE for the selective sour-gas absorption in a solvent blend of MDEA and additionally sulfolane; these results are presented in Figure 10.

![](https://cdn.mathpix.com/cropped/ab722c49-d22f-4bae-b719-9d92d8fdb4fd-08.jpg?height=327&width=432&top_left_y=1221&top_left_x=410)
Figure 10. Partial pressure of a sour gas vs. sour-gas loading in the system water + MDEA ( $20.9 \mathrm{wt} \%$ ) + sulfolane ( $30.5 \mathrm{wt} \%$ ) + $\mathrm{CO}_{2}+ \mathrm{H}_{2} \mathrm{~S}$. Lines are predictions with ePC-SAFT using parameters from Tables $1-3$; symbols are experimental data. ${ }^{41}$ Triangles: quasiconstant $\mathrm{CO}_{2}$ loading ( $\alpha_{\mathrm{CO}_{2}}=0.35 \pm 0.05$ ), varying concentration of $\mathrm{H}_{2} \mathrm{~S}$ at 313 K ; stars: quasi-constant $\mathrm{H}_{2} \mathrm{~S}$ loading ( $\alpha_{\mathrm{H}_{2} \mathrm{~S}}=0.5 \pm 0.05$ ), varying concentration of $\mathrm{CO}_{2}$ at $373 \mathrm{~K} . \mathrm{ARD} \%(373 \mathrm{~K})=39.95 \%$; $\mathrm{ARD} \%(313 \mathrm{~K})=41.72 \%$.

The absorption behavior was predicted for the aqueous solvent system water + MDEA ( $20.9 \mathrm{wt} \%$ ) + sulfolane ( 30.5 $\mathrm{wt} \%$ ) and variation of one sour-gas loading at a quasi-constant sour-gas loading of the second sour gas ( $\alpha_{\mathrm{CO}_{2}}=0.35 \pm 0.05$; $\left.\alpha_{\mathrm{H}_{2} \mathrm{~S}}=0.5 \pm 0.05\right)$ under isothermal conditions. The temperatures were chosen according to the availability of experimental data, i.e., 373 K at quasi-constant $\mathrm{H}_{2} \mathrm{~S}$ loading and 313 K at quasi-constant loading of $\mathrm{CO}_{2}$. Note that the literature data were measured originally for a scattered variety of loadings of both sour gases. To allow graphical comparison with the ePC-SAFT predictions, sour-gas absorption within a constant sour-gas loading $\pm 0.05$ for one sour gas was selected from the experimental data. Note that the ePC-SAFT predictions were performed for the exactly correct given loadings, i.e., only the graphical representation underlies clustering the data at quasi-constant loadings. Nevertheless, the predictions are in very good agreement with the experimental data in the range of the studied conditions.

The general underestimation of $\mathrm{CO}_{2}$ partial pressure is also visible for the selective absorption in the presence of a second sour gas.
In summary, it might be concluded that ePC-SAFT is an accurate modeling tool that captures the effect of physical solvent blends on the VLE for selective absorption of sour gases at least qualitatively correct. This is the crucial prerequisite to serve as a powerful screening tool.

## - CONCLUSIONS

Sour-gas absorption remains a hot topic in times of climate change. Understanding the interactions within aqueous solvent blends of physical and chemical solvents is a key step toward development of systems that allow an efficient reduction of emissions. Holistic and predictive thermodynamic models give information on interactions and help reduce time-consuming experiments. In this work, the thermodynamic model ePCSAFT extended with the original Born term to account for the change in the Gibbs energy of solvation of ions in multicomponent solvent blends was applied. The influence of ions was neglected in the mixing rule for the dielectric constant, which was expressed only as a function of the solvent concentration. The predictive strength of ePC-SAFT was tested against experimental VLE data for single sour-gas absorption in an aqueous MDEA solvent of various conditions and in a solvent blend composed of $20.9 \mathrm{wt} \%$ MDEA and 30.5 $\mathrm{wt} \%$ sulfolane between 313 and 373 K . The promising solvent class of ILs was also included in this work to additionally stress ePC-SAFT toward sour-gas absorption in more complex mixtures. Additionally, VLE for selective sour-gas absorption was investigated in solvent blends containing sulfolane.
ePC-SAFT-predicted VLE for single-gas absorption of $\mathrm{CO}_{2}$ and $\mathrm{H}_{2} \mathrm{~S}$ was found to be in reasonably good agreement with experimental VLE data. The approach to include the Born term significantly reduces the number of binary interaction parameters. This is an important step for the development of a screening tool that is applied to multicomponent solvent blends. These binary interaction parameters have not been fitted to VLE data of the reacting system. The approach was finally applied to account for the presence of a second sour gas, i.e., absorption of $\mathrm{CO}_{2}$ and a quasi-constant loading of $\mathrm{H}_{2} \mathrm{~S}$ ( $\alpha =0.5)$ and vice versa at $\alpha_{\mathrm{CO}_{2}}=0.35$. For all systems under study, gas solubility was overestimated at high loadings due to overpredicted ion-solvent interactions, which causes underestimation of partial pressure of the sour gas especially at high loadings and low temperatures.

Overall, the prediction results are satisfying and qualitatively correct, i.e., the influence of temperature, loading, physical solvent, and a second sour gas on the pressure was predicted with qualitative agreement with the experimental data. That shows the significance of the Born term for electrolyte systems containing non-aqueous solvents in the blend. The Born term accounts for the change in ion-solvent interactions in multicomponent blends. This work supports the effort to reduce $\mathrm{CO}_{2}$ contamination of the environment and takes a step toward holistic thermodynamic prediction of sour-gas absorption for solvent screening.

## - ASSOCIATED CONTENT

## (I) Supporting Information

The Supporting Information is available free of charge at https://pubs.acs.org/doi/10.1021/acs.iecr.1c00176.

Chemical equilibrium constants, Antoine constants for sulfolane, and ARD values between model and experiments (Tables S1-S4) and density and vapor pressure of sulfolane and flowchart of calculation procedure (Figure S1-S3) (PDF)

## AUTHOR INFORMATION

## Corresponding Author

Christoph Held - Laboratory of Thermodynamics, Department of Biochemical and Chemical Engineering, TU Dortmund University, 44277 Dortmund, Germany; © orcid.org/0000-0003-1074-177X; Email: christoph.held@tu-dortmund.de

## Authors

Mark Bülow - Laboratory of Thermodynamics, Department of Biochemical and Chemical Engineering, TU Dortmund University, 44277 Dortmund, Germany
Nevin Gerek Ince - AVEVA Group plc, Cambridge CB3 OHB, England
Seiya Hirohama - AVEVA Group plc, Cambridge CB3 OHB, England
Gabriele Sadowski - Laboratory of Thermodynamics, Department of Biochemical and Chemical Engineering, TU Dortmund University, 44277 Dortmund, Germany; © orcid.org/0000-0002-5038-9152
Complete contact information is available at:
https://pubs.acs.org/10.1021/acs.iecr.1c00176

## Notes

The authors declare no competing financial interest.

## - ACKNOWLEDGMENTS

This work was supported by the German Science Foundation (DFG) within the priority program SPP 1708 "Material Synthesis Near Room Temperature" (grant HE 7165/7-1). The authors thank AVEVA for the financial support from their University Research and Development Programme.

## - REFERENCES

(1) Kim, Y. E.; Lim, J. A.; Jeong, S. K.; Yoon, Y. I.; Bae, S. T.; Nam, S. C. Comparison of Carbon Dioxide Absorption in Aqueous MEA, DEA, TEA, and AMP Solutions. Bull. Korean Chem. Soc. 2013, 34, 783-787.
(2) Haghtalab, A.; Izadi, A. Simultaneous measurement solubility of carbon dioxide+hydrogen sulfide into aqueous blends of alkanolamines at high pressure. Fluid Phase Equilib. 2014, 375, 181-190.
(3) Bates, E. D.; Mayton, R. D.; Ntai, I.; Davis, J. H. CO(2) capture by a task-specific ionic liquid. J. Am. Chem. Soc. 2002, 124, 926-927.
(4) Sharma, P.; Park, S. D.; Park, K. T.; Nam, S. C.; Jeong, S. K.; Yoon, Y. I.; Baek, I. H. Solubility of carbon dioxide in aminefunctionalized ionic liquids: Role of the anions. Chem. Eng. J. 2012, 193-194, 267-275.
(5) Austgen, D. M.; Rochelle, G. T.; Chen, C. C. Model of vaporliquid equilibria for aqueous acid gas-alkanolamine systems. 2. Representation of hydrogen sulfide and carbon dioxide solubility in aqueous MDEA and carbon dioxide solubility in aqueous mixtures of MDEA with MEA or DEA. Ind. Eng. Chem. Res. 1991, 30, 543-555.
(6) Barreau, A.; Le Blanchon Bouhelec, E.; Habchi Tounsi, K. N.; Mougin, P.; Lecomte, F. Absorption of H 2 S and CO 2 in Alkanolamine Aqueous Solution: Experimental Data and Modelling with the Electrolyte-NRTL Model. Oil Gas Sci. Technol. - Rev. IFP 2006, 61, 345-361.
(7) Ermatchkov, V.; Pérez-Salado Kamps, Á.; Maurer, G. Solubility of Carbon Dioxide in Aqueous Solutions of N -Methyldiethanolamine in the Low Gas Loading Region. Ind. Eng. Chem. Res. 2006, 45, 60816091.
(8) Faramarzi, L.; Kontogeorgis, G. M.; Thomsen, K.; Stenby, E. H. Extended UNIQUAC model for thermodynamic modeling of CO2 absorption in aqueous alkanolamine solutions. Fluid Phase Equilib. 2009, 282, 121-132.
(9) Pérez-Salado Kamps, Á.; Balaban, A.; Jödecke, M.; Kuranov, G.; Smirnova, N. A.; Maurer, G. Solubility of Single Gases Carbon Dioxide and Hydrogen Sulfide in Aqueous Solutions of N -Methyldiethanolamine at Temperatures from 313 to 393 K and Pressures up to 7.6 MPa : New Experimental Data and Model Extension. Ind. Eng. Chem. Res. 2001, 40, 696-706.
(10) Kuranov, G.; Rumpf, B.; Smirnova, N. A.; Maurer, G. Solubility of Single Gases Carbon Dioxide and Hydrogen Sulfide in Aqueous Solutions of N -Methyldiethanolamine in the Temperature Range 313-413 K at Pressures up to 5 MPa . Ind. Eng. Chem. Res. 1996, 35, 1959-1966.
(11) Sadegh, N.; Stenby, E. H.; Thomsen, K. Thermodynamic modeling of hydrogen sulfide absorption by aqueous N -methyldiethanolamine using the Extended UNIQUAC model. Fluid Phase Equilib. 2015, 392, 24-32.
(12) Sadegh, N.; Stenby, E. H.; Thomsen, K. Thermodynamic modelling of acid gas removal from natural gas using the Extended UNIQUAC model. Fluid Phase Equilib. 2017, 442, 38-43.
(13) Zhang, Y.; Chen, C.-C. Modeling Gas Solubilities in the Aqueous Solution of Methyldiethanolamine. Ind. Eng. Chem. Res. 2011, 50, 6436-6446.
(14) Zhang, Y.; Chen, C.-C. Thermodynamic Modeling for CO 2 Absorption in Aqueous MDEA Solution with Electrolyte NRTL Model. Ind. Eng. Chem. Res. 2011, 50, 163-175.
(15) Zoghi, A. T.; Feyzi, F.; Dehghani, M. R. Modeling CO 2 Solubility in Aqueous N -methyldiethanolamine Solution by Electrolyte Modified Peng-Robinson Plus Association Equation of State. Ind. Eng. Chem. Res. 2012, 51, 9875-9885.
(16) Zong, L.; Chen, C.-C. Thermodynamic modeling of CO2 and H2S solubilities in aqueous DIPA solution, aqueous sulfolane-DIPA solution, and aqueous sulfolane-MDEA solution with electrolyte NRTL model. Fluid Phase Equilib. 2011, 306, 190-203.
(17) Alkhatib, I. I. I.; Pereira, L. M. C.; Vega, L. F. 110th Anniversary: Accurate Modeling of the Simultaneous Absorption of H 2 S and CO 2 in Aqueous Amine Solvents. Ind. Eng. Chem. Res. 2019, 58, 6870-6886.
(18) Avlund, A. S.; Kontogeorgis, G. M.; Michelsen, M. L. Modeling Systems Containing Alkanolamines with the CPA Equation of State. Ind. Eng. Chem. Res. 2008, 47, 7441-7446.
(19) Chunxi, L.; Fürst, W. Representation of CO 2 and H 2 S solubility in aqueous MDEA solutions using an electrolyte equation of state. Chem. Eng. Sci. 2000, 55, 2975-2988.
(20) Cleeton, C.; Kvam, O.; Rea, R.; Sarkisov, L.; Angelis, M. G. de. Competitive $\mathrm{H} 2 \mathrm{~S}-\mathrm{CO} 2$ absorption in reactive aqueous methyldiethanolamine solution: Prediction with ePC-SAFT. Fluid Phase Equilib. 2020, 511, No. 112453.
(21) Huttenhuis, P. J. G.; Agrawal, N. J.; Solbraa, E.; Versteeg, G. F. The solubility of carbon dioxide in aqueous N -methyldiethanolamine solutions. Fluid Phase Equilib. 2008, 264, 99-112.
(22) Mac Dowell, N.; Llovell, F.; Adjiman, C. S.; Jackson, G.; Galindo, A. Modeling the Fluid Phase Behavior of Carbon Dioxide in Aqueous Solutions of Monoethanolamine Using Transferable Parameters with the SAFT-VR Approach. Ind. Eng. Chem. Res. 2010, 49, 1883-1899.
(23) Nasrifar, K.; Tafazzol, A. H. Vapor-Liquid Equilibria of Acid Gas-Aqueous Ethanolamine Solutions Using the PC-SAFT Equation of State. Ind. Eng. Chem. Res. 2010, 49, 7620-7630.
(24) Pahlavanzadeh, H.; Fakouri Baygi, S. Modeling CO2 solubility in Aqueous Methyldiethanolamine Solutions by Perturbed ChainSAFT Equation of State. J. Chem. Thermodyn. 2013, 59, 214-221.
(25) Uyan, M.; Sieder, G.; Ingram, T.; Held, C. Predicting CO2 solubility in aqueous N-methyldiethanolamine solutions with ePCSAFT. Fluid Phase Equilib. 2015, 393, 91-100.
(26) Wang, T.; El Ahmar, E.; Coquelet, C.; Kontogeorgis, G. M. Improvement of the PR-CPA equation of state for modelling of acid gases solubilities in aqueous alkanolamine solutions. Fluid Phase Equilib. 2018, 471, 74-87.
(27) Wangler, A.; Sieder, G.; Ingram, T.; Heilig, M.; Held, C. Prediction of CO2 and H2S solubility and enthalpy of absorption in reacting N-methyldiethanolamine /water systems with ePC-SAFT. Fluid Phase Equilib. 2018, 461, 15-27.
(28) Bülow, M.; Ascani, M.; Held, C. ePC-SAFT advanced - Part II: Application to Salt Solubility in Ionic and Organic Solvents and the Impact of Ion Pairing. Fluid Phase Equilib. 2021, 537, No. 112989.
(29) Bülow, M.; Ascani, M.; Held, C. ePC-SAFTadvanced - Part I: Physical meaning of including a concentration-dependentdielectric constant in the born term and in the Debye-Hückel theory. Fluid Phase Equilib. 2021, 535, No. 112967.
(30) Cameretti, L. F.; Sadowski, G. Modeling of aqueous amino acid and polypeptide solutions with PC-SAFT. Chem. Eng. Process. Process Intensif. 2008, 47, 1018-1025.
(31) Held, C.; Reschke, T.; Mohammad, S.; Luza, A.; Sadowski, G. ePC-SAFT revised. Chem. Eng. Res. Des. 2014, 92, 2884-2897.
(32) Vahidi, M.; Moshtari, B. Dielectric data, densities, refractive indices, and their deviations of the binary mixtures of N methyldiethanolamine with sulfolane at temperatures 293.15328.15 K and atmospheric pressure. Thermochim. Acta 2013, 551, 1-6.
(33) Born, M. Volumen und Hydratationswärme der Ionen. $Z$. Physik 1920, 1, 45-48.
(34) Berthelot, D. Sur le mélange des gaz. Compt. Rendus 1898, 126, 1703-1706.
(35) Lorentz, H. A. Ueber die Anwendung des Satzes vom Virial in der kinetischen Theorie der Gase. Ann. Phys. 1881, 248, 127-136.
(36) Wolbach, J. P.; Sandler, S. I. Using Molecular Orbital Calculations To Describe the Phase Behavior of Hydrogen-Bonding Fluids. Ind. Eng. Chem. Res. 1997, 36, 4041-4051.
(37) Saleh, M. A.; Shamsuddin Ahmed, M.; Begum, S. K. Density, viscosity and thermodynamic activation for viscous flow of water +sulfolane. Phys. Chem. Liq. 2006, 44, 153-165.
(38) Gross, J.; Sadowski, G. Perturbed-Chain SAFT: An Equation of State Based on a Perturbation Theory for Chain Molecules. Ind. Eng. Chem. Res. 2001, 40, 1244-1260.
(39) Ji, X.; Held, C.; Sadowski, G. Modeling imidazolium-based ionic liquids with ePC-SAFT. Fluid Phase Equilib. 2012, 335, 64-73.
(40) Ji, X.; Held, C. Modeling the density of ionic liquids with ePCSAFT. Fluid Phase Equilib. 2016, 410, 9-22.
(41) Jou, F.-Y.; Deshmukh, R. D.; Otto, F. D.; Mather, A. E. Solubility of $\mathrm{H} 2 \mathrm{~S}, \mathrm{CO} 2, \mathrm{CH} 4$ and C 2 H 6 in sulfolane at elevated pressures. Fluid Phase Equilib. 1990, 56, 313-324.
(42) Kleiner, M.; Sadowski, G. Modeling of Polar Systems Using PCP-SAFT: An Approach to Account for Induced-Association Interactions. J. Phys. Chem. C 2007, 111, 15544-15553.
(43) Tommila, E.; Lindell, E.; Virtalai, M.; Laakso, R. Density, viscosities, surface tensions, dielectric constants, vapor pressures, activities and heats of mixing of sulphlane - water, sulpholane methanol, and sulpholane - ethanol mixtures. Suom. Kemistil. 1969, 42, 95.
(44) Bülow, M.; Greive, M.; Zaitsau, D. H.; Verevkin, S. P.; Held, C. Extremely Low Vapor-Pressure Data as Access to PC-SAFT Parameter Estimation for Ionic Liquids and Modeling of Precursor Solubility in Ionic Liquids. ChemistryOpen 2021, 10, 216-226.
(45) Shojaeian, A.; Haghtalab, A. Solubility and density of carbon dioxide in different aqueous alkanolamine solutions blended with 1 -butyl-3-methylimidazolium acetate ionic liquid at high pressure. J. Mol. Liq. 2013, 187, 218-225.


[^0]:    ${ }^{a} k_{i j}$ Calculated this way are valid in the temperature range $293-413 \mathrm{~K}$. Only valid with the pure-component parameters in Tables 1 and 2.

