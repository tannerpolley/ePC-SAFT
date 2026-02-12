\title{
ePC-SAFT revised
}

\author{
Christoph Held, Thomas Reschke, Sultan Mohammad, Armando Luza, Gabriele Sadowski* \\ Laboratory of Thermodynamics, Department of Biochemical and Chemical Engineering, Technische Universität Dortmund, Emil-Figge-Str. 70, 44227 Dortmund, Germany
}

\begin{abstract}
So far, the electrolyte PC-SAFT equation of state developed in Cameretti et al. (2005) has been applied to model solution densities, vapor-liquid equilibria (VLE), liquid-liquid equilibria (LLE), and solid-liquid equilibria (SLE) of solutions containing electrolytes. For that purpose, two ion-specific parameters were used to characterize any ion: the diameter of the solvated ion and the dispersion-energy parameter between ion and solvent. Dispersion was only considered between ions and solvents. Considering the small number of adjustable parameters, this approach yielded acceptable results especially for low and moderate electrolyte concentrations. However, for high salt concentrations, a distinct deviation between modeled and experimental data was observed.

In this work a new modeling approach is suggested that accounts explicitly also for dispersion interactions between anions and cations. This yields a much more precise description of electrolyte solutions at higher concentrations compared to original ePC-SAFT. With this new approach it is also possible to directly model weak electrolyte solutions without using an additional approach that accounts for ion-pair formation.

The new approach for applying ePC-SAFT is now able to model phase equilibria (VLE, LLE, SLE) of ternary electrolyte solutions containing water, organic solvents, salts, and amino acids even at high salt concentrations in good agreement with experimental data.
\end{abstract}
© 2014 The Institution of Chemical Engineers. Published by Elsevier B.V. All rights reserved.

Keywords: Electrolytes; Osmotic coefficients; Activity coefficients; Solubility; LLE; VLE

\section*{1. Introduction}

Electrolytes are present in many chemical and biochemical processes and systems, e.g. in electrolysis (Luckas and Krissmann, 2001) or in aqueous two-phase systems (He et al., 2005), and in biological media in which additionally biochemical reactions may take place. In technical processes electrolytes may be present invariably or added intentionally. Due to the remarkable impact of electrolytes on such processes and due to the high thermodynamic non-ideality of electrolyte solutions, modeling thermodynamic properties and phase equilibria in electrolyte solutions gains increasing consideration.

Modeling such systems requires considering short-range (SR) interactions resulting from repulsion and van der Waals attraction, polar, or association interactions. Further, the longrange (LR) interactions among the charged species need to be accounted for. Debye and Hückel proposed already in 1923 their approach for describing dilute electrolyte solutions (Debye and Hückel, 1923). They considered ions to be point charges that are shielded within a certain radius which is denoted as the reversed Debye screening length $\kappa$. Debye and Hückel assumed that ions only interact via LR forces due to their charges. The model of Debye and Hückel was shown to agree well (Maribo-Mogensen et al., 2012) with mean spherical approximations (Waisman and Lebowitz, 1970). For the LR

\footnotetext{
Abbreviations: $\mathrm{Ac}^{-}$, acetate anion; ARD, absolute average relative deviation; BuOH , 1-butanol; EOS, equation of state; $\mathrm{Fo}^{-}$, formate anion; $g^{\mathrm{E}}$, excess Gibbs energy; LLE, liquid-liquid equilibria; OF, objective function; SAFT, Statistical Associating Fluid Theory; SLE, solid-liquid equilibria; VLE, vapor-liquid equilibria.
* Corresponding author. Tel.: +49 2317552635.

E-mail address: g.sadowski@bci.tu-dortmund.de (G. Sadowski).
Received 4 February 2014; Received in revised form 20 May 2014; Accepted 23 May 2014
Available online 2 June 2014
http://dx.doi.org/10.1016/j.cherd.2014.05.017
0263-8762/© 2014 The Institution of Chemical Engineers. Published by Elsevier B.V. All rights reserved.
}

\section*{Nomenclature}

Variables
a Helmholtz free energy per number of molecules [-]
a activity [-]
$d_{i} \quad$ temperature-dependent segment diameter of molecule i [ ${ }^{\circ} \mathrm{A}$ ]
$f$ fugacity [Pa]
$\Delta h^{\mathrm{SL}} \quad$ melting enthalpy $[\mathrm{kJ} / \mathrm{kg}]$
$k_{\mathrm{B}} \quad$ Boltzmann constant, $1.38065 \times 10^{-23} \mathrm{~J} / \mathrm{K}[\mathrm{J} / \mathrm{K}]$
$k_{i j} \quad$ binary interaction parameter used in Eq. (16) [-]
$k_{i j, T} \quad$ parameter accounting for the temperature dependency of the $k_{\mathrm{ij}}$ parameter [ $1 / \mathrm{K}$ ]
$k_{\mathrm{ij}, 298.15 \mathrm{~K}}$ binary interaction parameter at 298.15 K [-]
$k_{i j}^{\mathrm{hb}} \quad$ binary interaction parameter used in Eq. (17) [-]
$K_{\text {acid }} \quad$ acid constant $[\mathrm{mol} / \mathrm{L}]$
$K_{\text {diss }} \quad$ dissociation constant [-]
$l_{i j} \quad$ binary interaction parameter used in Eq. (15) [-]
m molality (moles solute i per kg solvent) [mol/kg]
M molecular weight [ $\mathrm{g} / \mathrm{mol]}$
$m^{\text {seg }}$ number of segments [-]
$n_{i} \quad$ mole number of component $i[-]$
N number of association sites [-]
NP number of data points [-]
$p$ pressure [Pa]
$R \quad$ ideal gas constant $[\mathrm{J} / \mathrm{mol} / \mathrm{K}]$
T temperature [K]
$\mathrm{T}^{\text {SL }} \quad$ melting temperature [K]
$u / k_{\mathrm{B}}$ dispersion-energy parameter [K]
$x$ mole fraction [-]
Z compressibility factor [-]
Greek symbols

\begin{tabular}{|l|l|}
\hline $\gamma_{\mathrm{i}}$ & symmetrical activity coefficient of component i (related to pure component) [-] \\
\hline $\varphi_{\mathrm{i}}$ & fugacity coefficient of component i [-] \\
\hline $\varepsilon_{\mathrm{r}}$ & relative dielectric constant [-] \\
\hline $\varepsilon^{A_{i} B_{i} / k_{\mathrm{B}}}$ & association-energy parameter [K] \\
\hline $\kappa^{A_{i} B_{i}}$ & association-volume parameter [-] \\
\hline $\phi$ & osmotic coefficient [-] \\
\hline $\rho$ & density [ $\mathrm{kg} / \mathrm{m}^{3}$ ] \\
\hline $\nu$ & stoichiometric factor [-] \\
\hline $\sigma_{\mathrm{i}}$ & temperature-independent segment diameter of molecule $i[\AA]$ \\
\hline
\end{tabular}

\section*{Subscripts}
i,j component indices
T function of temperature
seg segment
w water
0 pure substance

\section*{Superscripts}
assoc association
disp dispersion
hc hard chain
L liquid phase
oc osmotic coefficient
res residual
+, - positive or negative charge
forces, it was shown that the dielectric constant of an electrolyte solution strongly influences the results with electrolyte models (Maribo-Mogensen et al., 2012, 2013).

Since the 1970s electrolyte models have been developed which consider also SR forces. These can be described by using either excess Gibbs energy models ( $g^{\mathrm{E}}$ models) or equations of state (EOS). Examples for electrolyte $g^{\mathrm{E}}$ models are the electrolyte NRTL model (Chen et al., 1982) or the Pitzer model (Pitzer, 1973). Examples for electrolyte EOS are an electrolyte Peng-Robinson EOS (Myers et al., 2002) and various electrolyte SAFT-based (Statistical Associating Fluid Theory) approaches (Held et al., 2008; Ji and Adidharma, 2006; Galindo et al., 1999; Tan et al., 2008). Besides, existing electrolyte models differ in their approach for describing the SR forces (e.g. Renon and Prausnitz, 1968; Abrams and Prausnitz, 1975; Chapman et al., 1989; Gross et al., 2001) as well as the LR forces (e.g. Debye and Hückel, 1923; Waisman and Lebowitz, 1970). Another difference among these models is the description of ion dissociation. Salts might be considered as fully dissociated into cations and anions, as being completely non-dissociated, or as partly dissociated. These modeling approaches require different types of component parameters, either ion-specific parameters (e.g. Cameretti et al., 2005), salt-specific parameters (Myers et al., 2002), or the combination of both (e.g. Chapman et al., 1990; Ji et al., 2005). Moreover, the description of ion solvation differs in different approaches. One possibility is applying a chemical theory (e.g. Myers et al., 2002). However, in most cases, solvation is described via SR forces without using a specific chemical theory. Thus, only within SAFT-based models ion-water interactions are described via dispersion SR forces (electrolyte Perturbed-Chain SAFT, ePCSAFT (Cameretti et al., 2005)), associative SR forces (SAFT for potentials of variable range and electrolytes, SAFT-VRE (Galindo et al., 1999; Gil-Villegas et al., 2001)), or even both (Galindo et al., 1999; Gil-Villegas et al., 2001; Herzog et al., 2010; Rozmus et al., 2013; Lee and Kim, 2009). Moreover, the existing models differ in the number and kind of model parameters. In some models the ion diameters are treated as adjustable parameters (Cameretti et al., 2005), whereas in others the experimental values for bare-ion radii (Pauling) (Rozmus et al., 2013) or solvated-ion (Stokes) radii (Papaiconomou et al., 2002) are used. In some models, the ion parameters are even defined as function of salt concentration (e.g. Papaiconomou et al., 2002). A review about the number of parameters applied within SAFT-based approaches was given by Tan et al. (2008), and further comparisons between the different electrolyte SAFT models can be found in the work of Rozmus et al. (2013).

Although many approaches to model electrolyte solutions have been developed so far, this research field still reveals potential for improvements. This is mainly due to the fact that most models were applied to binary salt/water solutions only, without validating their applicability to systems containing more components. Moreover, many published papers on electrolyte models contain parameters only for a few salts or ions, so that the range of applicability is limited. From scientific point of view and for the broad application in industry however, it is highly desirable to apply one single model to (1) many different systems in a broad temperature and concentration range, (2) different properties/equilibria calculations, (3) systems where no experimental data are available (model predictions), and (4) conditions at which parameters have not been adjusted (model extrapolations).

One model that has been shown to fulfill these criteria is ePC-SAFT (Cameretti et al., 2005). Within ePC-SAFT, the SR
forces are described by PC-SAFT (Gross et al., 2001). The LR forces among ions are described using the approach of Debye and Hückel (Held et al., 2008). Accounting for the formation of ion pairs, also weak electrolyte solutions could be modeled successfully (Held and Sadowski, 2009). Applying ePC-SAFT, 120 aqueous electrolyte solutions composed of 24 ions were described requiring 48 ion parameters only. Based on this, the salt influence on the phase behavior of biomolecule/water (Held et al., 2010, 2014; Sadeghi et al., 2012) systems could be predicted in reasonable agreement with experimental data.

Although ePC-SAFT has successfully been applied to many electrolyte systems, there are cases for which the model accuracy is still not satisfying. Especially solutions at high salt concentrations and solutions containing weakly-hydrated ions like phosphates or sulfates could not be described accurately with original ePC-SAFT (Held et al., 2008). The average relative deviation (ARD) between experimental and osmotic coefficients modeled using ePC-SAFT is larger than $6 \%$ (Held et al., 2008). Further, to describe weak electrolyte solutions, a chemical approach had to be combined with ePC-SAFT for an accurate data description (Held and Sadowski, 2009). Moreover, temperature effects described with original ePC-SAFT are only in poor agreement with experimental data (Held et al., 2008).

The aim of this work is to develop a modeling approach within the ePC-SAFT framework that corrects for these shortcomings. This will be done by considering SR dispersion forces also between anions and cations in addition to the interactions already accounted for by original ePC-SAFT (Held et al., 2008).

\section*{2. ePC-SAFT equation of state}

\subsection*{2.1. The model}

In this work electrolyte PC-SAFT as proposed by Cameretti et al. (2005) was used to model electrolyte solutions. In addition to original PC-SAFT (Gross et al., 2001), the LR forces are explicitly taken into account within ePC-SAFT. Thus, the residual Helmholtz energy $a^{\text {res }}$ of an electrolyte system is calculated by ePC-SAFT as:
$$
\begin{equation*}
a^{\mathrm{res}}=a^{\mathrm{hc}}+a^{\mathrm{disp}}+a^{\mathrm{assoc}}+a^{\mathrm{ion}} \tag{1}
\end{equation*}
$$
ePC-SAFT is based on a perturbation theory. The hardchain system represents the reference system and is described via the Helmholtz-energy contribution $a^{\text {hc }}$. Perturbations that ePC-SAFT explicitly accounts for are molecular attractions due to dispersive van der Waals ( $a^{\text {disp }}$ ) and associative hydrogenbonding forces ( $a^{\text {assoc }}$ ). For a system with charged species, a Debye-Hückel term is used to describe the Helmholtz-energy contribution $a^{\text {ion }}$. All contributions are treated as additive contributions that can be considered independently.

Within ePC-SAFT, ions are treated as spherical species in a uniform dielectric continuum which is characterized by a dielectric constant $\varepsilon$. They can approach each other to a certain distance which is assumed to be the solvated-ion diameter $\sigma_{\mathrm{j}}$. Repulsive interactions and the van der Waals attraction between ions and solvent (solvation) are accounted for in $a^{\text {hc }}$ and $a^{\text {disp }}$, respectively. The contribution $a^{\text {ion }}$ requires information about the dielectric constant of the solution but no additional adjustable parameters. According to our previous work (Held et al., 2008), an electrolyte-independent dielectric constant of the solvent was used again in this work.

\subsection*{2.2 Calculation of thermodynamic properties with ePC-SAFT}

The description of phase equilibria with EOS-type models requires fugacity coefficients $\varphi_{i}$ of the components $i$. These are accessible from the residual Helmholtz energy $a^{\text {res }}$ according to
$$
\begin{equation*}
\ln \varphi_{i}=(Z-1)-\ln Z+a^{\text {res }}+\frac{\partial a^{\text {res }}}{\partial x_{i}}-\sum_{j} X_{j}\left(\frac{\partial a^{\text {res }}}{\partial x_{j}}\right) \tag{2}
\end{equation*}
$$
where $x_{i}$ is the mole fraction of component $i$. The compressibility factor $Z$ is obtained by
$$
\begin{equation*}
Z(\rho)=1+\rho\left(\frac{\partial a^{\text {res }}}{\partial \rho}\right)_{T, x} \tag{3}
\end{equation*}
$$

Based on the fugacity coefficients, also activity coefficients $\gamma_{i}$ can be calculated by
$$
\begin{equation*}
\gamma_{i}=\frac{\varphi_{i}(\mathrm{~T}, p, \vec{x})}{\varphi_{0 i}\left(\mathrm{~T}, p, x_{i}=1\right)} \tag{4}
\end{equation*}
$$
where the subscript $0 i$ denotes the reference state pure component at same temperature $T$ and pressure $p$ as the actual solution of composition $\vec{x}$.

Another possibility for the characterization of solutions containing (various) solutes is using the osmotic coefficient $\phi$. This property applies to single-solvent solutions and is used in this work for the description of electrolyte solutions. Whereas water activity coefficients in such solutions often only marginally deviate from unity, osmotic coefficients describe the thermodynamic non-ideality of aqueous solutions with considerably higher sensitivity. The molal osmotic coefficient used in this work relates to the water activity coefficient $\gamma_{\mathrm{w}}$ according to
$$
\begin{equation*}
\phi=-\frac{1000 \cdot \ln \left(x_{w} \gamma_{w}\right)}{M_{w} \sum_{j \neq w} m_{j}} \tag{5}
\end{equation*}
$$
where $M_{W}$ is the molar mass of water $(g / \mathrm{mol})$. The molality $m_{j}$ of a solute $j$ is given in moles of species $j$ per kg water.

Calculating the solute solubility in a solvent requires an equilibrium condition between the liquid solution and the solid phase. Assuming a pure solid phase and neglecting the influence of different heat capacities of solid and liquid solute, the mole fraction of the solute in the liquid phase $x_{i}^{L}$ (its solubility) can be calculated by Prausnitz (1969).
$$
\begin{equation*}
x_{i}^{\mathrm{L}}=\frac{1}{\gamma_{i}} \cdot \exp \left\{-\frac{\Delta h_{0 i}^{\mathrm{SL}}}{\mathrm{RT}}\left(1-\frac{\mathrm{T}}{\mathrm{T}_{0 \mathrm{i}}^{\mathrm{SL}}}\right)\right\} \tag{6}
\end{equation*}
$$
$\Delta h_{0 i}^{\mathrm{SL}}$ and $\mathrm{T}_{0 i}^{\mathrm{SL}}$ are the melting enthalpy and the melting temperature of the solute, respectively.

For calculation of liquid-liquid equilibria, the isofugacity conditions were applied:
$$
\begin{equation*}
f_{i}^{\mathrm{L} 1}\left(\mathrm{~T}, p, \vec{x}^{\mathrm{L} 1}\right)=f_{i}^{\mathrm{L} 2}\left(\mathrm{~T}, p, \vec{x}^{\mathrm{L} 2}\right) \tag{7}
\end{equation*}
$$
where $f_{i}^{\mathrm{L} 1}$ and $f_{i}^{\mathrm{L} 2}$ denote the fugacity of component $i$ in the liquid phases L1 and L2, respectively. The fugacity in the liquid phases is calculated by
$$
\begin{equation*}
f_{i}^{\mathrm{L}}\left(\mathrm{T}, p, \vec{x}^{\mathrm{L}}\right)=p x_{i}^{\mathrm{L}} \varphi_{i}^{\mathrm{L}}\left(\mathrm{T}, p, \vec{x}^{\mathrm{L}}\right) \tag{8}
\end{equation*}
$$

Fugacity coefficients of volatile components required in Eq. (8) were calculated according to Eq. (2) whereas the fugacity of a salt or an acid ( $f_{ \pm}$) was calculated by
$$
\begin{equation*}
f_{ \pm}=\left(f_{+}^{v_{+}} f_{-}^{v_{-}}\right)^{\left(1 /\left(v_{+}+v_{-}\right)\right)} \tag{9}
\end{equation*}
$$
where $v_{+}$and $v_{-}$are the number of cations and anions of the salt or acid.

Considering phosphate solutions requires knowing the amount of water present at certain salt concentrations. The mole number of water $n_{\mathrm{w}}$ is reduced dramatically when a weak or medium-strong acid is added into the solution. For the case of an electrolyte containing $\mathrm{H}^{+}$ions (acids), it is assumed that all protons react to the hydronium ion $\mathrm{H}_{3} \mathrm{O}^{+}$. Consequently, the mole numbers of water $n_{\mathrm{w}}$ and acid $n_{\mathrm{H}_{\nu+} \mathrm{An}_{\nu-}}$ in the electrolyte solution are reduced by the mole number of formed hydronium ions $n_{\mathrm{H}_{3} \mathrm{O}^{+}}$according to:
$$
\begin{equation*}
n_{\mathrm{w}}=n_{0, \mathrm{w}}-n_{\mathrm{H}_{3} \mathrm{O}^{+}} \tag{10}
\end{equation*}
$$
$$
\begin{equation*}
n_{\mathrm{H}_{\nu+} \mathrm{An}_{\nu-}}=n_{0, \mathrm{H}_{\nu+}} \mathrm{An}_{\nu-}-n_{\mathrm{H}_{3} \mathrm{O}^{+}} \tag{11}
\end{equation*}
$$
where $n_{0, H_{v}+A n_{v-}}$ is the initial acid mole number, $n_{0, w}$ is $1000 / \mathrm{M}_{\mathrm{w}}$, and $n_{\mathrm{H}_{3} \mathrm{O}^{+}}$is obtained from the dissociation equilibrium
$$
\begin{equation*}
\mathrm{H}_{2} \mathrm{O}+\mathrm{H}_{v+} \mathrm{An}_{v-} \stackrel{K_{\text {diss }}}{\longleftrightarrow} v_{+} \mathrm{H}_{3} \mathrm{O}^{+}+v_{-} \mathrm{An} \tag{12}
\end{equation*}
$$
with $v_{+}$and $v_{-}$being the number of hydronium ions and anions in the acid $\mathrm{H}_{v+} \mathrm{An}_{v-}$. The dissociation constant $\mathrm{K}_{\text {diss }}$ is defined as
$$
\begin{equation*}
K_{\text {diss }}=\frac{n_{\mathrm{H}_{3} \mathrm{O}^{+}} n_{\mathrm{An}}^{v-}}{n_{w} n_{\mathrm{H}_{v+}} \mathrm{An}_{v-}} \tag{13}
\end{equation*}
$$
$K_{\text {diss }}$ is accessible via the molarity-based acid constant $K_{\text {acid }}^{c}$, given in the literature at infinite dilution. At infinite dilution, $K_{\text {diss }}$ in Eq. (13) becomes:
$$
\begin{equation*}
K_{\mathrm{diss}}=\frac{K_{\mathrm{acid}}^{\mathrm{c}}}{c_{0, \mathrm{w}}} \tag{14}
\end{equation*}
$$
using $c_{0, w}=55.366 \mathrm{~mol} / \mathrm{L}$. The acid constants $\mathrm{K}_{\text {acid }}^{\mathrm{c}}$ at 298.15 K used in this work for phosphoric acid were $10^{-2.16}$ and $10^{-7.21}$, respectively. Combining Eqs. (11) and (14) yields a non-linear expression for $n_{\mathrm{H}_{3} \mathrm{O}^{+}}$as function of $K_{\text {diss }}$, which was applied in this work in order to correct for the concentration effects based on the reduction of $n_{\mathrm{w}}$ upon phosphate-salt addition.

\subsection*{2.3 Modeling strategies}

Describing electrolyte solutions with ePC-SAFT requires a modeling approach for solvents and electrolytes. Two strategies in the ePC-SAFT framework were investigated in this work. Strategy 1 was published in (Held et al., 2008), whereas strategy 2 was developed in this work.

The reason for establishing strategy 2 was the higher deviation between ePC-SAFT (strategy 1) modeled and experimental osmotic coefficients at high electrolyte concentrations as well as for weak electrolyte solutions. These deviations can be explained by the fact that LR interactions among ions can be physically-meaningful described by the Debye-Hückel approach only at low electrolyte concentrations. At high electrolyte concentrations, ions cannot be
assumed to be completely solvated any more. SR interactions among ions become more important. A similar effect is present in weak-electrolyte solutions, in which ions do not dissociate completely but interact strongly with each other and even form ion pairs. For both applications (high electrolyte concentrations and weak electrolyte solutions) anion-cation interactions can obviously neither be neglected nor be described with Debye-Hückel alone. Considering additional SR anion-cation attractive forces should correct for this shortcoming and thus should allow for improving the modeling results for electrolyte solutions at high electrolyte concentrations and also for weak electrolyte solutions.

In both strategies, hard-chain repulsion, dispersive attraction, and association are considered for solvents. Further, in both strategies, the considered salts and acids are assumed to be completely dissociated into their ions. Both strategies consider repulsion among ions as well as between solvents and ions. The diameters of solvated ions are assumed to be temperature independent. In both strategies, dispersion interactions between ions and solvents are accounted for and the LR interactions among ions are described by a Debye-Hückel contribution $a^{\text {ion }}$ to the residual Helmholtz energy. The application of the Debye-Hückel contribution requires the use of the solvent's dielectric constant, which is assumed to depend only on the solvent and temperature, i.e. being independent of salt or acid for both strategies.

The two strategies differ in modeling the dispersion interactions among ions. Whereas strategy 1 does not account for the dispersion among ions, strategy 2 explicitly considers dispersion between anions and cations. Dispersion between two anions or two cations is not considered neither in strategy 1 nor strategy 2.

In both strategies, conventional Lorentz-Berthelot - combining rules are used for estimating the segment diameter $\sigma_{\mathrm{ij}}$ and the dispersive energy $u_{i j}$ between two components $i$ and $j$ (e.g. water and ion):
$$
\begin{align*}
& \sigma_{i j}=\frac{1}{2}\left(\sigma_{i}+\sigma_{j}\right) \cdot\left(1-l_{i j}\right)  \tag{15}\\
& u_{i j}=\sqrt{u_{i} u_{j}} \cdot\left(1-k_{i j}\right) \tag{16}
\end{align*}
$$

Whereas binary interaction parameters $k_{i j}$ between ions and solvents at 298.15 K were not used in strategy $1, k_{\mathrm{ij}}$ values were applied for correcting the cross-dispersive interactions between solvent and ions as well as between cation and anion in strategy 2. Modeling an electrolyte solution containing a salt or acid applying strategy 2 thus requires three binary $k_{i j}$ parameters (anion-solvent, cation-solvent, anion-cation) whereas $k_{i j}$ s for ion/solvent pairs at 298.15 K and for anion/cation pairs were not used in strategy 1.

Considering compounds that may form hydrogen bonds requires the use of combining rules for the association parameters. In this work, the Wolbach-Sandler rules were applied (Wolbach and Sandler, 1998):
$$
\begin{equation*}
\varepsilon^{A_{i} B_{j}}=\frac{1}{2}\left(\varepsilon^{A_{i} B_{i}}+\varepsilon^{A_{j} B_{j}}\right)\left(1-k_{i j}^{h b}\right) \tag{17}
\end{equation*}
$$
$$
\begin{equation*}
k^{A_{i} B_{j}}=\sqrt{k^{A_{i} B_{i} k^{A_{j} B_{j}}}} \quad\left(\frac{\sqrt{\sigma_{i} \sigma_{j}}}{1 / 2\left(\sigma_{i}+\sigma_{j}\right)}\right)^{3} \tag{18}
\end{equation*}
$$

In order to model water/alcohol mixtures very precisely with PC-SAFT, $k_{i j}^{\mathrm{hb}}$ and $l_{\mathrm{ij}}$ binary parameters were already introduced in previously (Held et al., 2012; Nann et al., 2013a).

With this, PC-SAFT allows for a quantitative modeling the binary 1-butanol/water LLE. The so-obtained accuracy is a prerequisite for modeling ternary 1-butanol/water/salt systems.

\subsection*{2.4. Pure-component and binary ePC-SAFT parameters}

Modeling electrolyte solutions with ePC-SAFT requires pure-component parameters for solvents and ions. The purecomponent parameters for the solvents and the $k_{i j}$ interaction parameters between an organic solvent and water are listed in Table 1. Table 1 further contains the pure-component parameters and $k_{i j} \mathrm{~s}$ with water for glycine and alanine, which were also considered in this work.

The pure-component parameters of the ions were fitted to solution densities and osmotic coefficients of single-solute aqueous electrolyte solutions in both strategies (strategy 1 (Held et al., 2008) and strategy 2 [this work]). Although solution-density data and osmotic-coefficient data are not the most important properties for biological and biopharmaceutical applications, such data has been shown to be very suitable for the pure-component parameter estimation of ions (Held et al., 2008) and other solutes (Held et al., 2014, 2013; Nann et al., 2013b; Grosse Daldrup et al., 2011; Passos et al., 2014; Hoffmann et al., 2013). The data used for the parameter estimation was generally limited to atmospheric pressure, 298.15 K , and only data with salt concentrations up to 5 m were considered. The objective function for the density OF ${ }^{\text {density }}$ and the osmotic coefficients $\mathrm{OF}^{\circ \mathrm{C}}$ in the parameter-estimation procedure used in this work were:

OF ${ }^{\text {density }}$
$$
\begin{equation*}
=\sum_{i=1}^{\mathrm{NP}}\left[\left(\rho_{i}^{\text {solution,calc }}-\rho_{0, \text { water }}^{\text {calc }}\right)-\left(\rho_{i}^{\text {solution, exp }}-\rho_{0, \text { water }}^{\exp }\right)\right]^{2} \tag{19}
\end{equation*}
$$
$$
\begin{equation*}
\mathrm{OF}^{\mathrm{OC}}=\sum_{i=1}^{\mathrm{NP}}\left(1-\frac{\phi^{\mathrm{calc}}}{\phi^{\mathrm{exp}}}\right)^{2} \tag{20}
\end{equation*}
$$
where exp and calc are experimental data and ePC-SAFT calculated properties and the sums are built over all number of data points NP. In Eq. (19), the experimental densities of pure water were obtained from the equation of Wagner and Pruss (2002). The parameter-estimation strategy was carried out as follows. First, a sensitivity analysis has been performed on the dispersion-energy parameters $u_{\text {ion }} / k_{B}$ of the halide anions and of the alkali cations. For all ions, this parameter was not allowed to exceed 400 K . Next to the fact that very high dispersion-energy values are physically not meaningful, such high values might cause inconsistencies (e.g. liquid-liquid demixing of aqueous salt solutions beyond the (solid) salt solubility). After this, $k_{i j}$ values ( $k_{i j}$ water-ion and $k_{i j}$ anion-cation) and the segment diameters $\sigma_{\text {ion }}$ of alkali cations and halide anions were adjusted using the OF in Eqs. (19) and (20). In the next steps, the parameters ( $u_{\text {ion }} / k_{\mathrm{B}}, k_{\mathrm{ij}}$ and $\sigma_{\text {ion }}$ ) of all other anions and cations listed in Table 2 were adjusted successively using the OF in Eqs. (19) and (20).

Table 2 contains the pure-component parameters for the ions determined in previous work (Held et al., 2008) (strategy 1: left part of the table) and determined in this work (strategy 2: right part of the table) as well as binary interaction parameters $k_{i j}$ between water and the ions (only required for

\begin{tabular}{|l|l|l|l|l|l|l|}
\hline \multicolumn{7}{|c|}{Table 1 - PC-SAFT parameters for non-charged components used in this work.} \\
\hline Parameter & Water (Fuchs et al., 2006) & 1-Butanol (Nann et al., 2013a) & Benzene (Gross et al., 2001) & Toluene (Gross et al., 2001) & Glycine (Held et al., 2011) & Alanine (Held et al., 2011) \\
\hline $m_{i}^{\text {seg }}$ & 1.2047 & 2.7515 & 2.4653 & 2.8149 & 4.8507 & 5.4647 \\
\hline $\sigma_{i} / \AA$ & a & 3.6139 & 3.6478 & 3.7169 & 2.3270 & 2.5222 \\
\hline $u_{i} \cdot k_{\mathrm{B}}^{-1} / \mathrm{K}$ & 353.95 & 259.59 & 287.35 & 285.69 & 216.96 & 287.59 \\
\hline $N_{i}$ & 2 & 2 & - & - & 2 & 2 \\
\hline $\varepsilon^{A_{i} B_{i}} \cdot k_{B}^{-1} / \mathrm{K}$ & 2425.7 & 2544.6 & - & - & 2598.1 & 3176.6 \\
\hline $k^{A_{i} B_{i}}$ & 0.0451 & 0.0067 & - & - & 0.0393 & 0.0819 \\
\hline $k_{i j}$ with water & - & T $2.94 \times 10^{-4}$ & $6.07 \times 10^{-4} \times(T-298.15)$ & $7.26 \times 10^{-4} \times(T-298.15)$ & -0.0612 & T $2.91 \times 10^{-4}$ \\
\hline & & $-0.102^{\mathrm{b}}$ & +0.0255 & +0.021 & & -0.14796 \\
\hline \multicolumn{7}{|r|}{\multirow{2}{*}{\begin{tabular}{l}
${ }^{\mathrm{a}}$ The expression $\sigma=2.7927+10.11 \exp (-0.01775 \mathrm{~T})-1.417 \exp (0.01146 \mathrm{~T})$ was used (Fuchs et al., 2006). \\
${ }^{\mathrm{b}}$ For 1-butanol/water, additional binary interaction parameters were used: $l_{\mathrm{ij}}=-0.0044$ and $k_{\mathrm{ij}}^{\mathrm{hb}}=0.026$ (Nann et al., 2013a).
\end{tabular}}} \\
\hline & & & & & & \\
\hline
\end{tabular}

\begin{table}
\captionsetup{labelformat=empty}
\caption{Table 2 - ePG-SAFT parameters for ions as determined in a previous work (Held et al., 2008) (strategy 1, left part of the table) and determined in this work (strategy 2 , right part of the table). The parameters $\sigma_{\text {ion }}$ and $u_{\text {ion }}$ / $k_{\mathrm{B}}$ of strategy 2 are independent of the solvent. The ion parameters of strategy 1 and the $k_{\text {ion-water }}$ values are only valid in combination with the PG-SAFT parameters of water listed in Table 1.}
\begin{tabular}{|l|l|l|l|l|l|}
\hline \multicolumn{3}{|l|}{Strategy 1 previous work (Held et al., 2008)} & \multicolumn{3}{|l|}{Strategy 2 this work} \\
\hline Ion & $\sigma_{\text {ion }} / \AA$ & $u_{\text {ion }} \cdot k_{\mathrm{B}}^{-1} / \mathrm{K}$ & $\sigma_{\text {ion }} / \AA$ & $u_{\text {ion }} \cdot k_{\mathrm{B}}^{-1} / \mathrm{K}$ & $k_{\text {ion-water }}$ \\
\hline \multicolumn{6}{|l|}{Univalent anions} \\
\hline $\mathrm{H}_{3} \mathrm{O}^{+}$ & 2.2740 & 1616.49 & 3.4654 & 500.00 & 0.25 \\
\hline $\mathrm{Li}^{+}$ & 1.8177 & 2697.28 & 2.8449 & 360.00 & -0.25 \\
\hline $\mathrm{Na}^{+}$ & 2.4122 & 646.05 & 2.8232 & 230.00 & a \\
\hline $\mathrm{K}^{+}$ & 2.9698 & 271.05 & 3.3417 & 200.00 & b \\
\hline $\mathrm{Cs}^{+}$ & 3.5606 & 175.94 & 3.9246 & 180.00 & 0.081 \\
\hline $\mathrm{NH}_{4}{ }^{+}$ & 3.4755 & 212.36 & 3.5740 & 230.00 & 0.064 \\
\hline \multicolumn{6}{|l|}{Univalent anions} \\
\hline $\mathrm{F}^{-}$ & 1.6132 & 648.31 & 1.7712 & 275.00 & \\
\hline $\mathrm{Cl}^{-}$ & 3.0575 & 47.29 & 2.7560 & 170.00 & -0.25 \\
\hline $\mathrm{Br}^{-}$ & 3.4573 & 60.22 & 3.0707 & 190.00 & -0.25 \\
\hline $\mathrm{I}^{-}$ & 3.9319 & 80.43 & 3.6672 & 200.00 & -0.25 \\
\hline $\mathrm{NO}_{3}{ }^{-}$ & 3.3805 & 0.00 & 3.2988 & 130.00 & 0.098 \\
\hline $\mathrm{OH}^{-}$ & 1.6401 & 2444.76 & 2.0177 & 650.00 & -0.25 \\
\hline $\mathrm{ClO}_{4}{ }^{-}$ & 4.0731 & 58.42 & 4.0186 & 104.16 & -0.25 \\
\hline $\mathrm{HCO}_{3}{ }^{-}$ & n.a. & n.a. & 2.9296 & 70.00 & 0.00 \\
\hline $\mathrm{H}_{2} \mathrm{PO}_{4}{ }^{-}$ & 3.7020 & 60.00 & 3.6505 & 95.00 & 0.25 \\
\hline $\mathrm{Fo}^{-}$ & n.a. & n.a. & 3.3077 & 190.00 & -0.23 \\
\hline $\mathrm{Ac}^{-}$ & 4.3000 & 300.00 & 3.9328 & 150.00 & -0.23 \\
\hline \multicolumn{6}{|l|}{Bivalent cations} \\
\hline $\mathrm{Mg}^{2+}$ & 2.3229 & 8145.33 & 3.1327 & 1500.00 & -0.25 \\
\hline $\mathrm{Ca}^{2+}$ & 2.8889 & 2146.98 & 3.2648 & 1060.00 & $4.1 \times 10^{-3}$ \\
\hline $\mathrm{Cu}^{2+}$ & 2.6955 & 2545.14 & 3.8379 & 1610.90 & 0.25 \\
\hline $\mathrm{Zn}^{2+}$ & n.a. & n.a. & 2.9798 & 1250.00 & -0.25 \\
\hline \multicolumn{6}{|l|}{Bi-/trivalent anions} \\
\hline $\mathrm{SO}_{4}{ }^{2-}$ & 2.4538 & 0.00 & 2.6491 & 80.00 & 0.25 \\
\hline $\mathrm{S}_{2} \mathrm{O}_{3}{ }^{2-}$ & n.a. & n.a. & 3.1877 & 80.00 & 0.25 \\
\hline $\mathrm{CO}_{3}{ }^{2-}$ & n.a. & n.a. & 2.4422 & 249.26 & -0.25 \\
\hline $\mathrm{HPO}_{4}{ }^{2-}$ & 4.4608 & 0.00 & 2.1621 & 146.02 & 0.25 \\
\hline $\mathrm{PO}_{4}{ }^{3-}$ & n.a. & n.a. & 2.5516 & 310.00 & -0.25 \\
\hline
\end{tabular}
\end{table}
${ }^{\mathrm{a}}$ For water $/ \mathrm{Na}^{+}$, a temperature-dependent $k_{\mathrm{ij}}$ was applied: $k_{\mathrm{ij}}=-0.007981 \mathrm{~T} / \mathrm{K}+2.37999$.
${ }^{\mathrm{b}}$ For water $/ \mathrm{K}^{+}$, a temperature-dependent $\mathrm{k}_{\mathrm{ij}}$ was applied: $\mathrm{k}_{\mathrm{ij}}=-0.004012 \mathrm{~T} / \mathrm{K}+1.3959$.
strategy 2). Table 3 contains $k_{i j}$ parameters between two ions (only required for strategy 2 ). The $k_{i j}$ parameters between ion and water and among two ions (Tables 2 and 3) were fitted together with $\sigma_{\text {ion }}$ and $u_{\text {ion }} / k_{\mathrm{B}}$ simultaneously to solutiondensity data and osmotic-coefficient data. It might appear that
fitting $u_{\mathrm{ion}} / k_{\mathrm{B}}$ and the $k_{\mathrm{ij}}$ parameters simultaneously points to a redundant parameter-optimization problem (two parameters in one equation, Eq. (16)). Nevertheless, this is not the case as $u_{\text {ion }} / k_{\mathrm{B}}$ is an input in two different cross-dispersion interactions, namely $u_{\text {ion-water }} / k_{\mathrm{B}}$ as well as $u_{\text {anion-cation }} / k_{\mathrm{B}}$.

\begin{table}
\captionsetup{labelformat=empty}
\caption{Table 3 - Binary $k_{i j}$ interaction parameters between ions used in this work; only valid with pure-component parameter set of ions determined in this work (Table 2 strategy 2).}
\begin{tabular}{|l|l|l|l|l|l|l|l|l|l|l|}
\hline $k_{\text {cation-anion }}$ & $\mathrm{H}_{3} \mathrm{O}^{+}$ & $\mathrm{Li}^{+}$ & $\mathrm{Na}^{+}$ & $\mathrm{K}^{+}$ & $\mathrm{Cs}^{+}$ & $\mathrm{NH}_{4}{ }^{+}$ & $\mathrm{Mg}^{2+}$ & $\mathrm{Ca}^{2+}$ & $\mathrm{Cu}^{2+}$ & $\mathrm{Zn}^{2+}$ \\
\hline $\mathrm{F}^{-}$ & n.a. & n.a. & 0.665 & 1.000 & 1.000 & n.a. & n.a. & n.a. & n.a. & n.a. \\
\hline $\mathrm{Cl}^{-}$ & 0.654 & 0.669 & 0.317 & 0.064 & -0.417 & -0.566 & 0.817 & 1.000 & -0.216 & -0.705 \\
\hline $\mathrm{Br}^{-}$ & 0.645 & 0.591 & 0.290 & -0.102 & -0.670 & -0.639 & 0.752 & 0.993 & n.a. & n.a. \\
\hline $\mathrm{I}^{-}$ & 0.497 & 0.002 & 0.018 & -0.312 & -1.000 & -0.787 & 0.317 & 0.229 & n.a. & n.a. \\
\hline $\mathrm{ClO}_{4}{ }^{-}$ & 0.861 & 0.406 & -0.118 & n.a. & n.a. & -1.000 & 0.122 & 0.674 & n.a. & n.a. \\
\hline $\mathrm{NO}_{3}{ }^{-}$ & n.a. & 0.364 & -0.300 & -0.585 & -0.855 & -0.419 & 0.285 & -0.101 & n.a. & n.a. \\
\hline $\mathrm{H}_{2} \mathrm{PO}_{4}{ }^{-}$ & n.a. & n.a. & -0.071 & 0.018 & n.a. & -1.000 & n.a. & n.a. & n.a. & n.a. \\
\hline $\mathrm{HPO}_{4}{ }^{2-}$ & n.a. & n.a. & -1.000 & 1.000 & n.a. & -0.556 & n.a. & n.a. & n.a. & n.a. \\
\hline $\mathrm{PO}_{4}{ }^{3-}$ & n.a. & n.a. & -1.000 & 1.000 & n.a. & n.a. & n.a. & n.a. & n.a. & n.a. \\
\hline $\mathrm{Ac}^{-}$ & n.a. & -0.998 & 0.246 & 1.000 & 0.785 & n.a. & n.a. & n.a. & n.a. & n.a. \\
\hline $\mathrm{Fo}^{-}$ & n.a. & n.a. & -0.370 & 0.340 & n.a. & n.a. & n.a. & n.a. & n.a. & n.a. \\
\hline $\mathrm{OH}^{-}$ & n.a. & -0.566 & 0.649 & 1.000 & 0.564 & n.a. & n.a. & n.a. & n.a. & n.a. \\
\hline $\mathrm{SO}_{4}{ }^{2-}$ & n.a. & -0.952 & -1.000 & 1.000 & -1.000 & -1.000 & -1.000 & -0.908 & -1.000 & -0.446 \\
\hline $\mathrm{CO}_{3}{ }^{2-}$ & n.a. & n.a. & -1.000 & 1.000 & n.a. & n.a. & n.a. & n.a. & n.a. & n.a. \\
\hline $\mathrm{HCO}_{3}{ }^{-}$ & n.a. & n.a. & -0.514 & -0.476 & n.a. & n.a. & n.a. & n.a. & n.a. & n.a. \\
\hline $\mathrm{S}_{2} \mathrm{O}_{3}{ }^{2-}$ & n.a. & n.a. & 0.184 & n.a. & n.a. & n.a. & n.a. & n.a. & n.a. & n.a. \\
\hline
\end{tabular}
\end{table}

\begin{table}
\captionsetup{labelformat=empty}
\caption{Table 4 - Binary interaction parameters between an ion and an organic solvent used in this work, only valid with parameter set of solvents in Table 1 and ions in Table 2 (strategy 2) and Table 3.}
\begin{tabular}{|l|l|l|l|}
\hline Ion & Solvent & $k_{i j}$ & $l_{i j}$ \\
\hline $\mathrm{NH}_{4}{ }^{+}$ & 1-Butanol & 0.29 & 0.140 \\
\hline $\mathrm{Cl}^{-}$ & 1-Butanol & 0.22 & 0.245 \\
\hline $\mathrm{Na}^{+}$ & Benzene & 0.35 & - \\
\hline $\mathrm{Cl}^{-}$ & Benzene & 0.35 & - \\
\hline $\mathrm{Br}^{-}$ & Benzene & 0.15 & - \\
\hline $\mathrm{I}^{-}$ & Benzene & 0.07 & - \\
\hline $\mathrm{SO}_{4}{ }^{2-}$ & Benzene & 1.00 & - \\
\hline $\mathrm{Na}^{+}$ & Toluene & 0.35 & - \\
\hline $\mathrm{Cl}^{-}$ & Toluene & 0.35 & - \\
\hline $\mathrm{Br}^{-}$ & Toluene & 0.15 & - \\
\hline
\end{tabular}
\end{table}

Table 4 lists binary interaction parameters ( $k_{\mathrm{ij}}$ and $l_{\mathrm{ij}}$ ) between organic solvents and ions (only required for strategy 2). The pure-component parameters for the ions $\sigma_{\text {ion }}$ and $u_{\text {ion }} / k_{\mathrm{B}}$ as well as the $k_{\mathrm{ij}}$ parameters between the ions determined in strategy 2 (listed in Tables 2 and 3) are universally valid independent of the kind of the electrolyte or solvent (ion parameters also used in organic solvent/salt solutions, see Section 3.3). The ion parameters determined in strategy 1 (Table 2) are limited to solutions containing water as the only solvent. These ion parameters (strategy 1) and the $k_{\text {ion-water }}$ parameters (strategy 2) are valid only in combination with the PC-SAFT parameters for water listed in Table 1.

Next to the pure-component and binary interaction parameters, expressions for the temperature-dependent diameter of a pure component are required for PC-SAFT calculations. In case of the cations and anions in this work, the expression
$$
\begin{equation*}
d_{\text {ion }}=\sigma_{\text {ion }}[1-0.12] \tag{21}
\end{equation*}
$$
was used as in all our previous papers on electrolytes (Held et al., 2008, 2010, 2012, 2014; Held and Sadowski, 2009; Sadeghi et al., 2012, 2014). For all other pure components considered in this work, the following expression was applied:
$$
\begin{equation*}
d_{i}=\sigma_{i}\left[1-0.12 \exp \left(-\frac{3 u_{i}}{k_{\mathrm{B}} \mathrm{~T}}\right)\right] \tag{22}
\end{equation*}
$$

The expression in Eq. (22) was used also in the original PC-SAFT publication (Gross et al., 2001), and was derived originally from Chen and Kreglewski (1977) in order to describe the soft repulsion of molecules. Eqs. (21) and (22) differ from each other as the pure-dispersion interactions among two equal ions were set to zero in ePC-SAFT resulting in a zero $u_{i} / k_{B}$ value in Eq. (22) finally yielding Eq. (21) for all ions.

In both strategies, the relative dielectric constant for the electrolyte solutions (containing water and 1-butanol in this work) was assumed to be independent of added salt or acid. A value for $\varepsilon_{\mathrm{r}, 1 \text {-butanol }}$ of 17.51 at 298.15 K (Chu et al., 1987) was used, and for $\varepsilon_{\text {r,water }}$ the relation given in (Held et al., 2008) was applied. For modeling 1-butanol/water mixtures containing salts, the dielectric constants were weight-averaged according to:
$$
\begin{equation*}
\varepsilon_{\mathrm{r}, 1 \text {-butanol } / \text { water }}=w_{1 \text {-butanol }} \cdot \varepsilon_{\mathrm{r}, 1 \text {-butanol }}+w_{\text {water }} \cdot \varepsilon_{\mathrm{r}, \text { water }} \tag{23}
\end{equation*}
$$

Comparing the pure-ion parameters of the two strategies (Table 2) let conclude the following facts:
- Dispersion-energy parameter of ions $u_{\text {ion }} / k_{B}$ : The cation parameter values are smaller for strategy 2 than for strategy 1. This is especially valid for small cations with high charge
density, for which the $u_{\text {ion }} / k_{\mathrm{B}}$ values could be reduced dramatically (e.g. by one order of magnitude for $\mathrm{Li}^{+}$). High $u_{\text {ion }} / k_{\mathrm{B}}$ values are not physically meaningful and might cause model inconsistencies (e.g. liquid-liquid demixing of aqueous salt solutions). This was e.g. the reason for refitting the $u_{\mathrm{Na}^{+}} / k_{\mathrm{B}}$ parameter in (Naeem and Sadowski, 2010) to be 250 K (compared to 646 K in original strategy 1 (Held et al., 2008)). In contrast, the dispersion-energy parameters of the anions increased compared to strategy 1.
- Despite the differences in their absolute values, the $u_{\text {ion }} / k_{\mathrm{B}}$ parameters of the ions are ranked in the same order in both, strategy 1 and 2. That is, they increase with cation charge density ( $\mathrm{K}^{+}<\mathrm{Na}^{+}<\mathrm{Li}^{+}$) and they increase with anion size $\mathrm{Cl}^{-}<\mathrm{Br}^{-}<\mathrm{I}^{-}$, which has already been discussed previously (Held et al., 2008).
- Some of the ions were modeled in strategy 1 by setting $u_{\text {ion }} / k_{\mathrm{B}}$ to zero (e.g. $\mathrm{SO}_{4}{ }^{2-}$ ). This strategy has not been followed for strategy 2 due to the fact that by applying mixing rules, zero values also cause dispersive interactions between these ions and any other component being zero, which is unreasonable in many cases and leads to bad descriptions, e.g. for ternary systems.
- Diameter of the solvated ions $\sigma_{\mathrm{ion}}$ : the diameter of the solvated ions increases in the same order as the Pauling diameters. This is physically meaningful and true for both strategies. The diameters of the solvated cations are larger for strategy 2 than for strategy 1, especially for $\mathrm{Li}^{+}$. The diameters of solvated $\mathrm{Na}^{+}$and solvated $\mathrm{Li}^{+}$do not differ much, which is in accordance to recent MD studies on electrolyte solutions (Deublein et al., 2012).

Note, that the ion parameters (fitted to solution-density data and osmotic-coefficient data) allow modeling the thermodynamic properties of the ions, such as the mean ionic activity coefficients of salts in water (e.g. shown also in (Held et al., 2008)) as well as salt solubility in water. These results are not shown in this work.

\section*{3. Results and discussion}

The main focus of this work was to accurately model phase equilibria of ternary electrolyte solutions with ePC-SAFT with strategy 2. In a first step, solution densities and osmotic coefficients of binary single-solute aqueous electrolyte systems were modeled. Based on this, the ternary mixtures biomolecule/salt/water and organic solvent/salt/water were modeled. For these mixtures, liquid-liquid equilibria (LLE), solid-liquid equilibria (SLE), and vapor-liquid equilibria (VLE) were considered.

All results were compared to original ePC-SAFT (strategy 1) in terms of the ARD in order to show the benefit of the new modeling approach:
$$
\begin{equation*}
\mathrm{ARD}=100 \cdot \frac{1}{\mathrm{NP}} \sum_{k=1}^{\mathrm{NP}}\left|\left(1-\frac{\mathrm{y}_{k}^{\text {calc }}}{\mathrm{y}_{k}^{\text {exp }}}\right)\right| \tag{24}
\end{equation*}
$$
where $y$ denotes the considered properties (e.g. density, osmotic coefficient, solubility).

\subsection*{3.1. Aqueous electrolyte solutions}

\subsection*{3.1.1. Strong-electrolyte solutions}

Aqueous solutions containing the ions ( $\mathrm{Li}^{+}, \mathrm{Na}^{+}, \mathrm{K}^{+}, \mathrm{Ca}^{2+}$, $\mathrm{Mg}^{2+}, \mathrm{NH}_{4}{ }^{+}, \mathrm{Cl}^{-}, \mathrm{Br}^{-}, \mathrm{I}^{-}$, acetate ${ }^{-}\left(\mathrm{Ac}^{-}\right)$, formate $\left.\mathrm{Fo}^{-}\right), \mathrm{F}^{-}$, $\mathrm{OH}^{-}, \mathrm{NO}_{3}{ }^{-}, \mathrm{SO}_{4}{ }^{2-}, \mathrm{HSO}_{4}{ }^{-}, \mathrm{PO}_{4}{ }^{3-}, \mathrm{H}_{2} \mathrm{PO}_{4}{ }^{-}, \mathrm{HPO}_{4}{ }^{2-}$ ) were considered. Solution densities and osmotic coefficients of almost 100 aqueous electrolyte solutions were investigated up to salt or acid concentrations of 5 m .

Table 5 compares the ARDs for aqueous electrolyte solutions with respect to solution densities and osmotic coefficients. Experimental data are compared to ePC-SAFT modeling results using strategy 1 and strategy 2 . Applying the original parameters from strategy 1 yields overall ARD values of $5.04 \%$ for osmotic coefficients of the univalent-cation electrolyte solutions up to 5 m at 298.15 K (Table 5). In contrast, the new parameters determined in this work (Table 2 (strategy 2) and Table 3) allow for much better modeling of these osmotic coefficients with overall ARD values of $1.57 \%$ only (Table 5). This pronounced improvement in the modeling is due to the better description of osmotic coefficients at moderate to high salt or acid concentrations ( $1-5 \mathrm{~m}$ ). Further, modeling osmotic coefficients of solutions containing bivalent anions ( $\mathrm{SO}_{4}{ }^{2-}$ and $\mathrm{HPO}_{4}{ }^{2-}$ ) is remarkably improved with strategy 2 (average ARD $=4 \%$, Table 5) compared to strategy 1 (average $\mathrm{ARD}=15 \%$, Table 5). Modeling solution densities of the univalent-cation electrolyte solutions up to 5 m at 298.15 K containing the ions considered in Table 5 can be modeled by both strategies with similar ARDs of $0.71 \%$ (strategy 1, Table 5) and $0.60 \%$ (strategy 2 , Table 5 ), respectively.

To sum up, the results obtained with strategy 2 are highly accurate, and of equal or better quality compared to other electrolyte SAFT models published in literature (see also Tan et al., 2008). Further, strategy 2 can also be applied to weak electrolyte solutions (ARD osmotic coefficient $<1 \%$ for acetate and formate solutions). Such a result could be obtained with strategy 1 only by accounting for the formation of ion pairs and additional introduction of ion-pair species requiring also ePC-SAFT parameters for the ion pairs.

Fig. 1 compares as example the modeling results of strategy 1 and strategy 2 for osmotic coefficients of NaCl /water and KBr /water solutions at 298.15 K . Both strategies describe the typical behavior of strong electrolyte solutions in a correct manner. That is, osmotic coefficients first decrease at low salt concentrations. After passing a minimum, the values increase continuously. By comparing both strategies it can be easily observed from Fig. 1 that modeling strategy 2 is more precise, especially at higher salt concentrations ( $>1 \mathrm{~m}$ ). The higher precision is at cost of using more model parameters. The accurate description applying strategy 2 is especially important for systems that contain additional components besides water and salt (shown in Section 3.2). In these systems, modeling results

\begin{figure}
\includegraphics[alt={},max width=\textwidth]{https://cdn.mathpix.com/cropped/4b629789-ae7b-4b68-92b5-97b9fc58b080-08.jpg?height=554&width=687&top_left_y=228&top_left_x=1146}
\captionsetup{labelformat=empty}
\caption{Fig. 1 - Osmotic coefficients of aqueous solutions of NaCl and KBr at 298.15 K . The symbols represent experimental data from Lobo and Quaresma (1989) (KBr: squares, NaCl: circles). The lines are modeling results with strategy 1 (thin lines and parameters from Held et al. (2008), Table 2 (strategy 1)) and with strategy 2 (bold lines and parameters from Table 2 (strategy 2) and Table 3).}
\end{figure}
strongly depend on the accuracy of the binary salt/solvent system.

Osmotic coefficients of electrolyte solutions are known to be only slightly temperature dependent. In this work, osmotic coefficients of $\mathrm{NaCl} /$ water and $\mathrm{KCl} /$ water solutions were modeled with strategy 2 at 273.15 K and at 298.15 K , illustrated in Fig. 2. It can be observed that experimental osmotic coefficients increase within the temperature range considered. In order to describe this behavior with ePC-SAFT, a temperature-dependent $k_{i j}$ between water and cation was applied (see Table 1), whereas the $k_{i j}$ value between water and the $\mathrm{Cl}^{-}$anion was considered to be independent of temperature. This strategy was chosen as alkali cations are generally stronger hydrated than halide anions; the temperature influence on ion hydration can thus be assumed to be stronger for the cation than for the anion. By completely neglecting

\begin{figure}
\includegraphics[alt={},max width=\textwidth]{https://cdn.mathpix.com/cropped/4b629789-ae7b-4b68-92b5-97b9fc58b080-08.jpg?height=618&width=769&top_left_y=1818&top_left_x=1105}
\captionsetup{labelformat=empty}
\caption{Fig. 2 - Osmotic coefficients of binary aqueous solutions of NaCl and KCl . The symbols represent experimental data (KCl: squares, NaCl: circles) from Lobo and Quaresma (1989) at 298.15 K (gray symbols) and at 273.15 K (black symbols), respectively. The lines are modeling results with strategy 2 (gray lines: 298.15 K , black lines: 273.15 K ) using the ion parameter set determined in this work (Tables 2 (strategy 2) and 3).}
\end{figure}

\begin{table}
\captionsetup{labelformat=empty}
\caption{Table 5 - ARDs (\%) between ePC-SAFT and experimental solution densities and osmotic coefficients of binary aqueous electrolyte solutions with univalent cations at 298.15 K for concentrations up to 5 m .}
\begin{tabular}{|l|l|l|l|l|l|l|l|l|l|l|l|l|}
\hline Salt or acid & $\mathrm{F}^{-}$ & $\mathrm{Cl}^{-}$ & $\mathrm{Br}^{-}$ & $\mathrm{I}^{-}$ & $\mathrm{ClO}_{4}{ }^{-}$ & $\mathrm{NO}_{3}{ }^{-}$ & $\mathrm{H}_{2} \mathrm{PO}_{4}{ }^{-}$ & $\mathrm{HPO}_{4}{ }^{2-}$ & $\mathrm{Ac}^{-}$ & $\mathrm{Fo}^{-}$ & $\mathrm{OH}^{-}$ & $\mathrm{SO}_{4}{ }^{2-}$ \\
\hline \multicolumn{13}{|l|}{ARD solution density (\%) for strategy 1} \\
\hline $\mathrm{H}^{+}$ & - & 0.20 & 0.39 & 0.16 & 0.54 & - & - & - & - & - & - & - \\
\hline $\mathrm{Li}^{+}$ & - & 0.58 & 0.37 & 1.07 & 1.88 & 0.53 & - & - & - & - & 0.20 & 2.29 \\
\hline $\mathrm{Na}^{+}$ & 0.24 & 0.52 & 1.06 & 0.97 & 0.49 & 0.38 & 0.27 & 1.82 & - & - & 0.34 & 0.22 \\
\hline $\mathrm{K}^{+}$ & 0.43 & 0.47 & 0.31 & 0.17 & - & 0.09 & 0.20 & 3.43 & - & - & 0.16 & 0.08 \\
\hline $\mathrm{Cs}^{+}$ & 0.21 & 0.22 & 0.17 & 0.06 & - & 0.03 & - & - & - & - & 0.21 & 0.40 \\
\hline $\mathrm{NH}_{4}{ }^{+}$ & - & 0.40 & 0.30 & 3.22 & 0.19 & 1.06 & - & - & - & - & - & 0.94 \\
\hline Average & 0.29 & 0.40 & 0.43 & 0.94 & 0.78 & 0.42 & 0.24 & 2.63 & - & - & 0.23 & 0.79 \\
\hline \multicolumn{13}{|l|}{ARD osmotic coefficient (\%) for strategy 1} \\
\hline $\mathrm{H}^{+}$ & - & 2.31 & 2.11 & 1.58 & 3.62 & - & - & - & - & - & - & - \\
\hline $\mathrm{Li}^{+}$ & - & 2.62 & 1.42 & 3.14 & 2.21 & 3.89 & - & - & - & - & 2.20 & 6.21 \\
\hline $\mathrm{Na}^{+}$ & 0.82 & 2.36 & 1.30 & 1.29 & 6.07 & 1.11 & 2.77 & 16.52 & - & - & 2.04 & 5.81 \\
\hline $\mathrm{K}^{+}$ & 1.81 & 2.31 & 1.64 & 0.60 & - & 1.17 & 2.25 & 30.66 & - & - & 1.22 & 0.78 \\
\hline $\mathrm{Cs}^{+}$ & 1.70 & 1.41 & 1.93 & 1.73 & - & 2.66 & - & - & - & - & 1.91 & 4.21 \\
\hline $\mathrm{NH}_{4}{ }^{+}$ & - & 0.80 & 0.99 & 1.74 & 3.27 & 19.15 & - & - & - & - & - & 15.14 \\
\hline Average & 1.44 & 1.97 & 1.57 & 1.68 & 3.79 & 5.60 & 2.51 & 23.59 & - & - & 1.84 & 6.43 \\
\hline \multicolumn{13}{|l|}{ARD solution density (\%) for strategy 2} \\
\hline $\mathrm{H}^{+}$ & - & 0.3 & 0.23 & 0.1 & 0.09 & - & - & - & - & - & - & - \\
\hline $\mathrm{Li}^{+}$ & - & 0.3 & 0.17 & 0.61 & 0.23 & 0.09 & - & - & - & - & 0.98 & 0.98 \\
\hline $\mathrm{Na}^{+}$ & 0.52 & 0.69 & 0.13 & 0.26 & 0.1 & 0.29 & 0.26 & 0.77 & 0.71 & 0.07 & 1.44 & 0.87 \\
\hline $\mathrm{K}^{+}$ & 2.67 & 0.44 & 0.04 & 0.37 & - & 0.14 & 0.18 & 0.86 & 0.23 & 0.37 & 1.25 & 0.58 \\
\hline $\mathrm{Cs}^{+}$ & 2.03 & 0.16 & 0.25 & 0.42 & - & 0.47 & - & - & - & - & 0.88 & 0.71 \\
\hline $\mathrm{NH}_{4}{ }^{+}$ & - & 0.45 & 0.18 & 1.38 & 0.45 & 1.13 & - & - & - & - & - & 0.63 \\
\hline Average & 1.74 & 0.39 & 0.18 & 0.52 & 0.22 & 0.42 & 0.22 & 0.82 & 0.47 & 0.22 & 1.14 & 0.75 \\
\hline \multicolumn{13}{|l|}{ARD osmotic coefficient (\%) for strategy 2} \\
\hline $\mathrm{H}^{+}$ & & 0.73 & 0.75 & 1.42 & 3.06 & - & - & - & - & - & - & - \\
\hline $\mathrm{Li}^{+}$ & - & 1.03 & 1.66 & 1.02 & 1.2 & 0.61 & - & - & 0.41 & - & 0.48 & 3.11 \\
\hline $\mathrm{Na}^{+}$ & 0.72 & 2.95 & 0.83 & 1.16 & 0.71 & 0.37 & 2.64 & 6.43 & 0.71 & 1.22 & 1.13 & 2.41 \\
\hline $\mathrm{K}^{+}$ & 2.03 & 0.51 & 0.47 & 0.36 & - & 0.89 & 0.75 & 4.95 & 0.54 & 1.1 & 0.77 & 0.42 \\
\hline $\mathrm{Cs}^{+}$ & 1.96 & 1.45 & 1.57 & 0.73 & - & 2.78 & - & - & 0.85 & - & 0.89 & 1.6 \\
\hline $\mathrm{NH}_{4}{ }^{+}$ & - & 0.64 & 0.6 & 0.63 & 1.07 & 0.54 & - & - & - & - & - & 3.22 \\
\hline Average & 1.57 & 1.22 & 0.22 & 0.89 & 1.51 & 1.04 & 1.7 & 5.69 & 0.63 & 1.16 & 0.82 & 2.15 \\
\hline
\end{tabular}
\end{table}
the temperature influence on anion hydration the number of adjustable parameters was further reduced. It can be seen that the $k_{i j}$ value for cation/water is only slightly temperature dependent. Nevertheless, this temperature dependence is required in order to describe osmotic coefficients at temperatures different from 298.15 K quantitatively. However, a temperature-dependent $k_{i j}$ is not necessary for strategy 2 for small temperature ranges ( $\Delta \mathrm{T}<25 \mathrm{~K}$ ) especially at small salt concentrations (<1 m).

Note, that precise osmotic-coefficient data are not available for most electrolytes at temperatures other than 298.15 K . Therefore, this work does not focus on determining temperature-dependent $k_{i j}$ values for all the ions considered in this work (Table 1). It should be emphasized that also strategy 1 can be applied to describe temperature-dependent osmotic coefficients of aqueous electrolyte solutions. For an accurate modeling result of $\mathrm{KCl} /$ water solutions at 273.15 K (results not shown), $k_{i j}$ values between $\mathrm{K}^{+}$and water as well as between $\mathrm{Cl}^{-}$and water of -0.1 were needed $\left(k_{i j}=0\right.$ at 298.15 K$)$ (Sadeghi et al., 2014). That is, the $k_{i j}$ between ions and water (which are zero at 298.15 K ) much stronger depends on temperature in strategy 1 compared to strategy 2.

\subsection*{3.1.2. Weak-electrolyte solutions}

The osmotic coefficients of salt/water solutions behave differently depending on the dissociation behavior of the salt.

Solutions in which the salt is completely dissociated (strong electrolytes) typically have high osmotic coefficients, caused by strong ion-solvent interactions. In contrast, ion-ion interactions are more pronounced in weak electrolyte solutions, which cause lower osmotic coefficients. This effect is e.g. known to reverse the series of osmotic coefficients: for aqueous solutions containing alkali chlorides, bromides, or iodides (X), experimental osmotic coefficients are sequenced in the order $\mathrm{LiX}>\mathrm{NaX}>\mathrm{KX}$. For the alkali acetates, the formation of ion pairs (strong ion-ion interactions) in sodium and lithium acetate solutions (Held and Sadowski, 2009) causes the reversal of these series, so that experimental osmotic coefficients are sequenced in the order $\mathrm{LiAc}<\mathrm{NaAc}<\mathrm{KAc}$ (Fig. 3).

Using the universal ion parameter set of strategy 1 (Table 2) did not even allow for a qualitative modeling of aqueous alkali acetate systems (Fig. 3a). In order to describe weak electrolyte solutions quantitatively with strategy 1 , a chemical model of ion pairing was combined with ePC-SAFT in a previous work (Held and Sadowski, 2009). Although this yielded meaningful results, the approach is not easy as a new species (the ion pair) has to be characterized, parameterized, and implemented into ePC-SAFT.

Thus, this work suggest a different approach: strategy 2 treats formation of ion pairs as strong anion-cation interaction which is accounted for by cross-dispersion interactions. The strength of this cross dispersion can be adjusted by the

\begin{figure}
\includegraphics[alt={},max width=\textwidth]{https://cdn.mathpix.com/cropped/4b629789-ae7b-4b68-92b5-97b9fc58b080-10.jpg?height=688&width=1115&top_left_y=219&top_left_x=482}
\captionsetup{labelformat=empty}
\caption{Fig. 3 - Osmotic coefficients of aqueous solutions of LiAc, NaAc, and KAc. The symbols represent experimental data from Lobo and Quaresma (1989) at 298.15 K (KAc: stars, NaAc: circles, LiAc: triangles). The lines are modeling results: KAc: dashed lines; NaAc: full lines; LiAc: thin lines. Experimental osmotic coefficients increase in the order LiAc<NaAc<KAc. (a) Calculations performed with the classical ePC-SAFT (parameters from Held and Sadowski (2009), Table 2 (strategy 1)). Modeled osmotic coefficients increase in the order $\mathrm{LiAc}>\mathrm{NaAc}>\mathrm{KAc}$. (b) Calculations performed with ePC-SAFT using strategy 2 and the ion parameter set determined in this work (Tables 2 (strategy 2) and 3). Modeled osmotic coefficients increase in the experimentally validated order $\mathrm{LiAc}<\mathrm{NaAc}<\mathrm{KAc}$.}
\end{figure}
$k_{\mathrm{ij}}$ value between anion and cation. This procedure is much easier to handle in phase-equilibrium calculations of multicomponent systems since it does not introduce additional species (ion pairs). Thus, the kind and number of model parameters for weak electrolytes are identical to those of strong electrolytes in strategy 2.

Applying strategy 2 to the aqueous acetate systems, osmotic coefficients of $\mathrm{Li}^{+}, \mathrm{Na}^{+}$, and $\mathrm{K}^{+}$acetate solutions can be modeled in almost quantitative agreement with the experimental data (Fig. 3b). Considering the $k_{i j}$ values between anions and cations, it can be observed that the order is $\mathrm{K}^{+}<\mathrm{Na}^{+}<\mathrm{Li}^{+}$for $k_{i j}$ values between these cations and chloride, bromide, or iodide. For the acetates, the $k_{i j}$ values are sequenced in the opposite order (see Table 3). This result is very reasonable as this means that cation-anion interactions increase in the order $\mathrm{LiX}<\mathrm{NaX}<\mathrm{KX}$ for $\mathrm{X}=$ halide anion but in the order $\mathrm{KAc}<\mathrm{NaAc}<\mathrm{LiAc}$. This finding is valid independent of the solvent, and agrees with investigations of ion interactions in aqueous (Lyklema, 1995) and non-aqueous (Lyklema, 2013) electrolyte solutions. It has to be noted here that KAc is not a weak electrolyte (only NaAc and LiAc are). Thus, $\mathrm{K}^{+}$and $\mathrm{Ac}^{-}$are dissociated completely in water, which might explain also the fact that the $k_{i j}$ between $\mathrm{K}^{+}$and $\mathrm{Ac}^{-}$is unity (Table 3) resulting in a zero dispersion energy between them.

\subsection*{3.2 Ternary biomolecule/salt/water solutions}

\subsection*{3.2.1. Amino acid/salt/water}

The (e)PC-SAFT parameters for amino acids and ions, which resulted from fitting to binary amino acid/water (parameters in Table 1) and binary salt/water solutions (parameters in Tables 1-3) were used for predicting properties of ternary amino acid/salt/water systems. For this purpose, osmotic coefficients at 298.15 K and solubilities of amino acids in electrolyte solutions were modeled.

Osmotic coefficients are appropriate quantities for understanding interactions in aqueous amino acid/salt solutions. This is due to the fact that many salts strongly influence
osmotic coefficients of biological solutions (Held et al., 2010, 2014; Sadeghi et al., 2012). Fig. 4a illustrates the influence of 1 m KCl on experimental osmotic coefficients of aminoacid solutions. Whereas experimental osmotic coefficients of alanine/ $\mathrm{KCl} /$ water solutions increase linearly with alanine molality at 1 m KCl , experimental osmotic coefficients of glycine/ $\mathrm{KCl} /$ water solutions at 1 m KCl pass a maximum. Applying strategy 2 (with parameters in Tables 1-3) to these systems allows for a reasonable prediction of the experimental data. That is, $k_{i j}$ values between ions and amino acids were not used. The average deviations between prediction and experiments are within the uncertainty of the experimental data $( \pm 1 \%)$. Slightly higher deviations for ePC-SAFT predicted osmotic coefficients were observed for the alanine/ KCl/water system compared to the water/glycine/ KCl system. Fig. 4b shows that the solubilities of glycine and alanine in water predicted with ePC-SAFT (strategy 2) decrease with increasing KCl molality. Also the experimental data show a salting-out tendency for both amino acids.

Comparing the results in this work (strategy 2 and bold lines in Fig. 4) with those of the original modeling (strategy 1 (Held et al., 2014) and thin lines in Fig. 4) it can be stated that both strategies can qualitatively predict osmotic coefficients in ternary amino acid/salt/water systems. In 36 systems containing one amino acid (glycine, alanine, valine, or proline), one salt ( $\mathrm{LiCl}, \mathrm{LiBr}, \mathrm{LiI}, \mathrm{NaCl}, \mathrm{NaBr}, \mathrm{NaI}, \mathrm{KCl}, \mathrm{KBr}, \mathrm{KI}, \mathrm{NaNO}_{3}, \mathrm{KNO}_{3}$, $\mathrm{Na}_{2} \mathrm{SO}_{4}$, $\left(\mathrm{NH}_{4}\right)_{2} \mathrm{SO}_{4}$ ) and water, the ARD values of osmotic coefficients are comparable (<4\% for strategy 1 (Held et al., 2014) and $<5 \%$ for strategy 2 ). For strategy 1 , modeling results at low salt concentrations are very accurate, whereas deviations are distinctly higher for salt molalities $>1 \mathrm{~m}$. This is different for strategy 2 . Strategy 2 yields deviations ( $\pm 2.5 \%$ ) that are independent of the salt molality. This constant deviation might be corrected with $k_{\mathrm{ij}}$ parameters between an ion and an amino acid, which, however, was not the focus of this work.

Fig. 4b compares the prediction of amino-acid solubilities in aqueous electrolyte systems using the parameters

\begin{figure}
\includegraphics[alt={},max width=\textwidth]{https://cdn.mathpix.com/cropped/4b629789-ae7b-4b68-92b5-97b9fc58b080-11.jpg?height=650&width=1116&top_left_y=219&top_left_x=468}
\captionsetup{labelformat=empty}
\caption{Fig. 4 - Thermodynamic properties of amino acid/salt/water solutions at 298.15 K . The symbols represent data from Held et al. (2014) (squares: glycine, circles: alanine), bold lines are predictions with ePC-SAFT using the parameters in Tables 1, 2 (strategy 2), and 3. Thin lines are predictions with ePC-SAFT using the parameters determined in previous work (Tables 1 and 2 (strategy 1)). (a) Osmotic coefficients at 1 m KCl . (b) Amino-acid solubility in $\mathrm{KCl} /$ water solutions.}
\end{figure}
from this work (strategy 2) with those of the original model (strategy 1). Fig. 4b shows that applying strategy 2 predicts a salting-out effect of KCl on glycine and alanine solubilities, which is in good accordance to literature data for alanine. Strategy 1 predicts a salting-in effect of KCl on glycine and alanine solubilities, which is in good accordance to literature data for glycine only at high salt concentrations ( $>3 \mathrm{~m}$ ).

To sum up, both strategies allow for qualitative predictions of osmotic coefficients and solubility in most ternary amino acid/salt/water systems. Aiming at a quantitative modeling of such systems would require adjusting binary interaction parameters $k_{i j}$ between ions and amino acids, which is not within the focus of this work. It should be moreover noted, that the introduction of such binary $k_{i j}$ values would have no effect for modeling electrolytes containing the ions $\mathrm{NO}_{3}{ }^{-}$and $\mathrm{SO}_{4}{ }^{2-}$ with strategy 1 as these ions were described in strategy 1 with $u / k_{\mathrm{B}}=0 \mathrm{~K}$, i.e. non-zero $k_{i j}$ values would not have any effect due to the validity of Eq. (16). At the other hand, this discovers the strongest advantage of strategy 2 , where special care was taken to obtain non-zero $u / k_{B}$ values for all ions.

\subsection*{3.3. Ternary solvent/salt/water solutions}
ePC-SAFT was applied to model the salt influence on the miscibility of organic solvents (benzene, toluene, 1-butanol) with water. Both, strategy 1 and strategy 2, were used for this purpose. It turned out that only strategy 2 allowed modeling the salt influence on the miscibility gap of an organic solvent with water with quantitative agreement to experimental data. The modeling with strategy 1 caused tie lines that are not in accordance with experimental data (tie-line length, tie-line slope). Only the results of strategy 2 are thus shown in this section.

\subsection*{3.3.1. Salt influence on benzene/water and toluene/water systems}

Literature data on the miscibility of benzene and toluene with water in the presence of salts is available only for the concentrations in the water-rich phase (Poulson et al., 1999; Wen-Hui et al., 1990; McDevit and Long, 1952; Šír et al., 1980). Independent of the kind and concentration of the salt, the experimental data reveals that the miscibility of toluene or
benzene with water strongly decreases upon salt addition. This is illustrated in Fig. 5 for the systems benzene/water/ NaCl and benzene/water/NaI. It can be observed from the experimental data that both salts very similarly decrease the solubility of benzene or toluene in water (salting-out effect), whereas NaCl has the strongest salting-out effect among the salts studied $\left(\mathrm{NaCl}, \mathrm{NaBr}, \mathrm{NaI}, \mathrm{Na}_{2} \mathrm{SO}_{4}\right)$.

Modeling ternary systems containing salt, water, and benzene or toluene requires pure-component and binary ePC-SAFT parameters. The pure-component parameters of benzene and toluene were taken from Gross et al. (2001), whereas the binary interaction parameters between benzene and water as well as between toluene and water were fitted in this work to benzene/water and toluene/water LLE data from Neely et al. (2008) and Jou and Mather (2003) between 297.8 and 490.8 K . This yielded mole-fraction based ARD values (Eq. (24)) of $48 \%$ (water-rich phase) and $1 \%$ (benzene-rich phase) for the water-benzene LLE and 58\% (water-rich phase) and 2\% (toluene-rich phase) for the water-toluene LLE. These high relative deviations are due to the very low absolute mole-fraction solubilities (absolute deviations are in the $10^{-4}$ order of

\begin{figure}
\includegraphics[alt={},max width=\textwidth]{https://cdn.mathpix.com/cropped/4b629789-ae7b-4b68-92b5-97b9fc58b080-11.jpg?height=511&width=704&top_left_y=1966&top_left_x=1119}
\captionsetup{labelformat=empty}
\caption{Fig. 5 - Solubility of benzene in aqueous salt solutions at 298.15 K and 1 bar as function of salt molality. The symbols represent experimental data (Poulson et al., 1999; Wen-Hui et al., 1990; McDevit and Long, 1952; Šír et al., 1980) (NaCl: circles, NaI: squares). Lines are modeling results with ePC-SAFT using the ion parameter set determined in this work (Tables 2 (strategy 2), 3, and 4).}
\end{figure}

\begin{table}
\captionsetup{labelformat=empty}
\caption{Table 6 - ARDs [\%] between ePC-SAFT and experimental benzene or toluene solubilities in aqueous electrolyte solutions at 298.15 K .}
\begin{tabular}{|l|l|l|l|}
\hline System & ARD (\%) & System & ARD (\%) \\
\hline \multirow[t]{2}{*}{Water/toluene/ NaCl} & \multirow[t]{2}{*}{0.47} & Water/benzene/NaCl & 7.54 \\
\hline & & Water/benzene/ NaBr & 1.25 \\
\hline \multirow[t]{2}{*}{Water/toluene/ NaBr} & \multirow[t]{2}{*}{1.52} & Water/benzene/NaI & 0.75 \\
\hline & & Water/benzene/ $\mathrm{Na}_{2} \mathrm{SO}_{4}$ & 5.43 \\
\hline
\end{tabular}
\end{table}
magnitude in mole fractions). The pure-component and binary ePC-SAFT parameters are given in Table 1. Applying these parameters (Table 1) and the ion-parameter sets in Table 2 (strategy 2) and 3 allowed for a qualitative prediction (binary interaction parameters between ions and benzene or ions and toluene set to zero) of the salting-out effect on benzene or toluene (results not shown).

In order to describe the salting-out effect on benzene or toluene with quantitative agreement to the experimental data, $k_{i j}$ parameters between ions and benzene or ions and toluene were adjusted in this work to benzene and toluene solubility data in aqueous salt solutions (Poulson et al., 1999; Wen-Hui et al., 1990; McDevit and Long, 1952; Šír et al., 1980). These $k_{i j}$ values are listed in Table 4. Due to the small concentrations of benzene and toluene in water at 298.15 K , the $k_{i j}$ parameters between ion and benzene and ion and toluene are not very sensitive to the results shown in Fig. 5, and large changes in the $k_{i j} \mathrm{~s}$ are needed to influence these modeling results. As it can be observed, all $k_{i j}$ values between ion and benzene or ion and toluene were set to positive values in order to reduce the dispersive attractions between them. This yields a more pronounced salting-out effect compared to the predictions ( $k_{i j}=0$ between ions and organic solvent). It can be observed from Table 4 that the $k_{i j}$ s between ions and benzene also apply to the respective toluene-ion pairs.

Table 6 lists the ARD values for the solubility of benzene or toluene in the considered aqueous salt solutions. The fact that the ARD values are higher than $5 \%$ in some cases can be explained by the very small absolute solubility values. Absolute deviations between experiments and ePC-SAFT were always in the order of $10^{-5} \mathrm{~mol} / \mathrm{kg}$ only.

The overall description of the miscibility of organic solvents with water in the presence of electrolytes leads to reasonable results and verifies the reasonability of the ion parameters (strategy 2) proposed in this work.

\subsection*{3.3.2. Salt influence on 1-butanol/water systems}

For the system water/1-butanol experimental data for the concentrations in both, water-rich phase and 1-butanol-rich phase, are available in the absence and presence of $\mathrm{NH}_{4} \mathrm{Cl}$ at 298.15 K (Pirahmadi et al., 2012).

The data in Fig. 6a shows that $\mathrm{NH}_{4} \mathrm{Cl}$ addition causes a salting out of 1-butanol in the aqueous phase. In turn, the 1 -butanol content in the organic phase is increased by $20 \mathrm{wt} \%$ upon addition of $20 \mathrm{wt} \% \mathrm{NH}_{4} \mathrm{Cl}$ to the salt-free water/1-butanol system (salting-out effect on water in the 1-butanol phase). Fig. 6b illustrates the mass-fraction based ternary liquid-liquid phase diagram for the water $/ \mathrm{NH}_{4} \mathrm{Cl} / 1$-butanol system including the experimental binodal curve and the tie lines in the two-phase region. It can be observed that the miscibility gap is increased upon salt addition. Modeling the ternary LLE water/salt/1-butanol first requires an accurate modeling of the binary LLE water/1-butanol. Pure-component parameters and binary interaction parameters for water and 1-butanol were taken from literature Nann et al. (2013a) and are given in Table 1. These parameters allow precisely modeling of the binary LLE water/1-butanol.

Moreover, modeling water/salt/organic solvent systems with strategy 2 requires binary interaction parameters between ions and the organic solvent. Thus, binary ion/1butanol interaction parameters were adjusted in this work to experimental LLE data of the water/ $\mathrm{NH}_{4} \mathrm{Cl} / 1$-butanol system at 298.15 K . These parameters are listed in Table 4 . The use of $l_{i j}$ parameters for modeling ion/alcohol solutions with ePC-SAFT is not new and has already been addressed in Held et al. (2012). The $l_{i j}$ parameter corrects for the diameter of the solvated ion according to Eq. (15). This binary parameter is required as $\sigma_{\text {ion }}$ was adjusted to properties in aqueous solutions. As the solvated ion diameter strongly depends on the kind of solvent, the simple mixing rule in Eq. (15) without accounting for $l_{i j}$ does not provide a correct picture. It was thus necessary to cor-

\begin{figure}
\includegraphics[alt={},max width=\textwidth]{https://cdn.mathpix.com/cropped/4b629789-ae7b-4b68-92b5-97b9fc58b080-12.jpg?height=481&width=1179&top_left_y=2083&top_left_x=443}
\captionsetup{labelformat=empty}
\caption{Fig. 6 - Liquid-liquid equilibrium of the system 1-butanol/ $\mathrm{NH}_{4} \mathrm{Cl} /$ water at 298.15 K at atmospheric pressure. The symbols represent experimental data (Pirahmadi et al., 2012) and lines are modeling results with ePC-SAFT and the ion parameter set determined in this work (Tables 1, 2 (strategy 2), 3, and 4). (a) wt.\% of 1-butanol in the organic and in the aqueous phase as function of $\mathrm{wt} . \% \mathrm{NH}_{4} \mathrm{Cl}$ in the aqueous phase (squares: 1-butanol in the aqueous phase, circles: 1-butanol in the 1-butanol phase). (b) Ternary phase diagram in mass fractions.}
\end{figure}
rect this by the $l_{i j}$ parameter. Due to the low solubility of salts in 1-butanol, the $l_{i j}$ parameters between ion and 1-butanol was adjusted to ternary systems with water.

Applying ePC-SAFT with strategy 2 (pure-component parameters in Tables 1 and 2 and the binary parameters in Tables 3 and 4) to the water $/ \mathrm{NH}_{4} \mathrm{Cl} / 1$-butanol LLE at 298.15 K yields a good agreement between modeled and experimental data. This is shown in Fig. 6. The weight fractions of 1-butanol in the aqueous phase and in the 1-butanol phase upon $\mathrm{NH}_{4} \mathrm{Cl}$ addition (Fig. 6a), the binodal curve (Fig. 6b) and the tie lines of the two-phase region (Fig. 6b) modeled with strategy 2 (parameters in Tables 1-4) are in good agreement with the experimental data. Although not shown, this is also valid for this system at a higher temperature $(308.15 \mathrm{~K})$, for which data is also available in the literature (Pirahmadi et al., 2012). For 308.15 K , also the parameters from Tables $1-4$ were applied. Temperature-dependent binary interaction parameters (ion/ion, ion/1-butanol, ion/water) were not applied.

To sum up, it can be stated that the new modeling strategy for ePC-SAFT allows modeling the salt influence on aqueous/organic LLE in good agreement with experimental data. At this point it should be mentioned again that the modeling with strategy 1 caused tie lines which were not in accordance with the experimental data (tie-line length, tie-line slope).

\section*{4. Conclusion}

Original ePC-SAFT from Held et al. (2008) has been applied to many electrolyte systems until now. However, the model accuracy was still not satisfying for modeling thermodynamic properties at high salt concentrations and for modeling weak electrolyte solutions (containing e.g. phosphates, sulfates, or acetates). The aim of this work was to quantitatively describe densities and phase equilibria of electrolyte solutions at high salt concentrations and weak-electrolyte solutions applying ePC-SAFT with a new modeling strategy for ions. Compared to original ePC-SAFT (strategy 1), this new strategy (strategy 2) explicitly accounts for dispersion interactions between anions and cations. The new strategy allowed modeling densities and osmotic coefficients of aqueous electrolyte solutions with ARD values of $0.6 \%$ and $1.57 \%$ at 298.15 K for more than 100 electrolyte solutions with molalities up to 5 m , including weak electrolyte solutions. These results are much better than those with strategy 1 (ARDs of $0.71 \%$ and $5.04 \%$ for solution densities and osmotic coefficients, respectively).

Based on the ion-parameter sets determined from binary salt/water solutions, ternary aqueous electrolyte solutions were modeled with both, strategy 1 and strategy 2 . These include systems containing water, salts, amino acids, and organic solvents. Modeling results were compared with experimental osmotic coefficients, solubility data and LLE data. Strategy 1 allowed for modeling the salt influence on osmotic coefficients and solubility in good agreement to the experimental data. However, the modeling results with strategy 1 did not show agreement to experimental salt/water/organic solvent LLE data. In contrast, modeling the salt influence on LLEs with strategy 2 did allow for quantitative results.

The results with strategy 2 were satisfactory at low and high electrolyte concentrations, for strong and weak electrolyte solutions, at different temperatures, for systems containing more than one solvent and more than one solute. Thus, the new strategy developed in this work can be applied to model different types of phase equilibria of electrolyte solutions in a broad range of components and conditions.

\section*{Acknowledgements}

The authors gratefully acknowledge the financial support of the IGF-project 17114N/1 of the DECHEMA e.V. that was funded by the Federal Ministry of Economics and Technology (grant no. 005-1009-0053) based on an enactment of the German Federal Parliament (BMWi). The project was supported within the program "Promoting the Industrial Collective Research (IGF) with the help of the German Federation of Industrial Research Associations (AiF). Further, Christoph Held acknowledges financial support from the European Union within the EFRE (European Fond for regional development) program "Ziel2NRW - Europe: Investments in our future".

\section*{References}

Abrams, D.S., Prausnitz, J.M., 1975. Statistical thermodynamics of liquid-mixtures - new expression for excess Gibbs energy of partly or completely miscible systems. AIChE J. 21, 116-128.
Cameretti, L.F., Sadowski, G., Mollerup, J.M., 2005. Modeling of aqueous electrolyte solutions with perturbed-chain statistical associated fluid theory. Ind. Eng. Chem. Res. 44, 3355-3362, ibid., 8944.
Chapman, W.G., Gubbins, K.E., Jackson, G., Radosz, M., 1989. SAFT - equation-of-state solution model for associating fluids. Fluid Phase Equilib. 52, 31-38.
Chapman, K.W., Johnson, W.C., Mclean, T.J., 1990. A high-speed statistical process-control application of machine vision to electronics manufacturing. Comput. Ind. Eng. 19, 234-238.
Chen, S.S., Kreglewski, A., 1977. Applications of augmented van der waals theory of fluids. 1. Pure fluids. Ber. Der Bunsen-Ges. - Phys. Chem. 81, 1048-1052.

Chen, C.C., Britt, H.I., Boston, J.F., Evans, L.B., 1982. Local composition model for excess Gibbs energy of electrolyte systems. Part I: single solvent, single completely dissociated electrolyte systems. AIChE J. 28, 588-596.
Chu, D.Y., Zhang, Q., Liu, R.L., 1987. Standard Gibbs free-energies of transfer of Nacl and Kcl from water to mixtures of the 4 isomers of butyl alcohol with water - the use of ion-selective electrodes to study the thermodynamics of solutions. J. Chem. Soc., Faraday Trans. I 83, 635-644.
Debye, P., Hückel, E., 1923. Zur Theorie der Elektrolyte. I. Gefrierpunktserniedrigung und verwandte Erscheinungen. Phys. Z. 24, 185-206.
Deublein, S., Vrabec, J., Hasse, H., 2012. A set of molecular models for alkali and halide ions in aqueous solution. J. Chem. Phys. 136, 084501.
Fuchs, D., Fischer, J., Tumakaka, F., Sadowski, G., 2006. Solubility of amino acids: influence of the pH value and the addition of alcoholic cosolvents on aqueous solubility. Ind. Eng. Chem. Res. 45, 6578-6584.
Galindo, A., Gil-Villegas, A., Jackson, G., Burgess, A.N., 1999. SAFT-VRE: phase behavior of electrolyte solutions with the statistical associating fluid theory for potentials of variable range. J. Phys. Chem. B 103, 10272-10281.
Gil-Villegas, A., Galindo, A., Jackson, G., 2001. A statistical associating fluid theory for electrolyte solutions (SAFT-VRE). Mol. Phys. 99, 531-546.
Gross, J., Sadowski, G., Perturbed-Chain, S.A.F.T., 2001. An equation of state based on a perturbation theory for chain molecules. Ind. Eng. Chem. Res. 40, 1244-1260.
Grosse Daldrup, J.B., Held, C., Sadowski, G., Schembecker, G., 2011. Modeling pH and solubilities in aqueous multisolute amino-acid solutions. Ind. Eng. Chem. Res. 50, 3503-3509.
He, C.Y., Li, S.H., Liu, H.W., Li, K., Liu, F., 2005. Extraction of testosterone and epitestosterone in human urine using aqueous two-phase systems of ionic liquid and salt. J. Chromatogr. A 1082, 143-149.

Held, C., Sadowski, G., 2009. Modeling aqueous electrolyte solutions. Part 2: weak electrolytes. Fluid Phase Equilib. 279, 141-148.
Held, C., Cameretti, L.F., Sadowski, G., 2008. Modeling aqueous electrolyte solutions: Part 1. Fully dissociated electrolytes. Fluid Phase Equilib. 270, 87-96.
Held, C., Neuhaus, T., Sadowski, G., 2010. Thermodynamic properties of aqueous ectoine, proline, and urea solutions measurement and modeling. Biophys. Chem. 152, 28-39.
Held, C., Cameretti, L.F., Sadowski, G., 2011. Measuring and modeling activity coefficients in aqueous amino-acid solutions. Ind. Eng. Chem. Res. 50, 131-141.
Held, C., Prinz, A., Wallmeyer, V., Sadowski, G., 2012. Measuring and modeling alcohol/salt systems. Chem. Eng. Sci. 68, 328-339.
Held, C., Carneiro, A., Macedo, E.A., Sadowski, G., 2013. Modeling thermodynamic properties of aqueous single-solute and multi-solute sugar solutions with PC-SAFT. AIChE J. 59, 4794-4805.
Held, C., Reschke, T., Müller, R., Kunz, W., Sadowski, G., 2014. Measuring and modeling aqueous electrolyte/amino-acid solutions with ePC-SAFT. J. Chem. Thermodyn. 68, 1-12.
Herzog, S., Gross, J., Arlt, W., 2010. Equation of state for aqueous electrolyte systems based on the semirestricted non-primitive mean spherical approximation. Fluid Phase Equilib. 297, 23-33.
Hoffmann, P., Voges, M., Held, C., Sadowski, G., 2013. The role of activity coefficients in bioreaction equilibria: thermodynamics of methyl ferulate hydrolysis. Biophys. Chem. 173, 21-30.
Ji, X.Y., Adidharma, H., 2006. Ion-based SAFT2 to represent aqueous single- and multiple-salt solutions at 298.15 K . Ind. Eng. Chem. Res. 45, 7719-7728.
Ji, X.Y., Tan, S.P., Adidharma, H., Radosz, M., 2005. Statistical associating fluid theory coupled with restricted primitive model to represent aqueous strong electrolytes: multiple-salt solutions. Ind. Eng. Chem. Res. 44, 7584-7590.
Jou, F.Y., Mather, A.E., 2003. Liquid-liquid equilibria for binary mixtures of water plus benzene, water plus toluene, and water plus p-xylene from 273 K to 458 K . J. Chem. Eng. Data 48, 750-752.
Lee, B.S., Kim, K.C., 2009. Modeling of aqueous electrolyte solutions based on perturbed-chain statistical associating fluid theory incorporated with primitive mean spherical approximation. Korean J. Chem. Eng. 26, 1733-1747.
Lobo, V.M.M., Quaresma, J.L., 1989. Handbook of Electrolyte Solutions Parts A and B. Elsevier, Amsterdam.
Luckas, M., Krissmann, J., 2001. Thermodynamik der Elektrolytlösungen: Eine einheitliche Darstellung der Berechnung komplexer Gleichgewichte. Springer, Berlin.
Lyklema, J., 1995. Fundamentals of Interface and Colloid Science: Solid-Liquid Interfaces. Academic Press, London.
Lyklema, J., 2013. Principles of interactions in non-aqueous electrolyte solutions. Curr. Opin. Colloid Interface Sci. 18, 116-128.
Maribo-Mogensen, B., Kontogeorgis, G.M., Thomsen, K., 2012. Comparison of the Debye-Hückel and the mean spherical approximation theories for electrolyte solutions. Ind. Eng. Chem. Res. 51, 5353-5363.
Maribo-Mogensen, B., Kontogeorgis, G.M., Thomsen, K., 2013. Modeling of dielectric properties of complex fluids with an equation of state. J. Phys. Chem. B 117, 3389-3397.
McDevit, W.F., Long, F.A., 1952. The activity coefficient of benzene in aqueous salt solutions. J. Am. Chem. Soc. 74, 1773-1777.
Myers, J.A., Sandler, S.I., Wood, R.H., 2002. An equation of state for electrolyte solutions covering wide ranges of temperature, pressure, and composition. Ind. Eng. Chem. Res. 41, 3282-3297.

Naeem, S., Sadowski, G., 2010. pePC-SAFT: modeling of polyelectrolyte systems. 1. Vapor-liquid equilibria. Fluid Phase Equilib. 299, 84-93.
Nann, A., Held, C., Sadowski, G., 2013a. Liquid-liquid equilibria of 1-butanol/water/IL systems. Ind. Eng. Chem. Res. 52, 18472-18481.
Nann, A., Mündges, J., Held, C., Verevkin, S., Sadowski, G., 2013b. Molecular interactions in 1-butanol + IL solutions by measuring and modeling activity coefficients. J. Phys. Chem. B 117, 3173-3185.
Neely, B.J., Wagner, J., Robinson, R.L., Gasem, K.A.M., 2008. Mutual solubility measurements of hydrocarbon-water systems containing benzene, toluene, and 3-methylpentane. J. Chem. Eng. Data 53, 165-174.
Papaiconomou, N., Simonin, J.P., Bernard, O., Kunz, W., 2002. MSA-NRTL model for the description of the thermodynamic properties of electrolyte solutions. Phys. Chem. Chem. Phys. 4, 4435-4443.
Passos, H., Khan, I., Mutelet, F., Oliveira, M.B., Carvalho, P.J., Santos, L.M.N.B.F., Held, C., Sadowski, G., Freire, M.G., Coutinho, J.A.P., 2014. Vapor-liquid equilibria of water plus alkylimidazolium-based ionic liquids: measurements and perturbed-chain statistical associating fluid theory modeling. Ind. Eng. Chem. Res. 53, 3737-3748.
Pirahmadi, F., Dehghani, M.R., Behzadi, B., 2012. Experimental and theoretical study on liquid-liquid equilibrium of 1-butanol+water+ $\mathrm{NH}_{4} \mathrm{Cl}$ at 298.15, 308.15 and 318.15 K . Fluid Phase Equilib. 325, 1-5.
Pitzer, K.S., 1973. Thermodynamics of electrolytes. I. Theoretical basis and general equations. J. Phys. Chem. 77, 268-277.
Poulson, S.R., Harrington, R.R., Drever, J.I., 1999. The solubility of toluene in aqueous salt solutions. Talanta 48, 633-641.
Prausnitz, J.M., 1969. Molecular Thermodynamics of Fluid-Phase Equilibria, 1st ed. Prentince Hall PTR, Upper Saddle River, NJ.
Renon, H., Prausnitz, J.M., 1968. Local compositions in thermodynamic excess functions for liquid mixtures. AIChE J. 14, 135-144.
Rozmus, J., de Hemptinne, J.C., Galindo, A., Dufal, S., Mougin, P., 2013. Modeling of strong electrolytes with ePPC-SAFT up to high temperatures. Ind. Eng. Chem. Res. 52, 9979-9994.
Sadeghi, M., Held, C., Samieenasab, A., Ghotbi, C., Abdekhodaie, M.J., Taghikhani, V., Sadowski, G., 2012. Thermodynamic properties of aqueous salt containing urea solutions. Fluid Phase Equilib. 325, 71-79.
Sadeghi, M., Held, C., Ghotbi, C., Abdekhodaie, M.J., Sadowski, G., 2014. Thermodynamic properties of aqueous glucose-urea-salt systems. J. Solution Chem., http://dx.doi.org/10.1007/s10953-014-0192-1, Accepted for publication.
Šír, Z., Ludmila, S., Rod, V., 1980. Contribution to the application of scaled particle theory to prediction of the salting coefficient. Collect. Czechoslov. Chem. Commun. 45, 679-689.
Tan, S.P., Adidharma, H., Radosz, M., 2008. Recent advances and applications of statistical associating fluid theory. Ind. Eng. Chem. Res. 47, 8063-8082.
Wagner, W., Pruss, A., 2002. The IAPWS formulation 1995 for the thermodynamic properties of ordinary water substance for general and scientific use. J. Phys. Chem. Ref. Data 31, 387-535.
Waisman, E., Lebowitz, J.L., 1970. Exact solution of an integral equation for structure of a primitive model of electrolytes. J. Chem. Phys. 52, 4307-4311.
Wen-Hui, X., Jing-Zhe, S., Xi-Ming, X., 1990. Studies on the activity coefficient of benzene and its derivatives in aqueous salt solutions. Thermochim. Acta 169, 271-286.
Wolbach, J.P., Sandler, S.I., 1998. Using molecular orbital calculations to describe the phase behavior of cross-associating mixtures. Ind. Eng. Chem. Res. 37, 2917-2928.