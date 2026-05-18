\title{
Highly efficient lithium extraction from magnesium-rich brines with ionic liquid-based collaborative extractants: Thermodynamics and molecular insights
}

\author{
Gangqiang Yu ${ }^{\mathrm{a}, \mathrm{b}}$, Xinhe Zhang ${ }^{\mathrm{a}}$, Tobias Hubach ${ }^{\mathrm{b}}$, Biaohua Chen ${ }^{\mathrm{a}, *}$, Christoph Held ${ }^{\mathrm{b}, *}$ \\ ${ }^{\mathrm{a}}$ Faculty of Environment and Life, Beijing University of Technology, 100 Ping Le Yuan, Chaoyang District, Beijing 100124, China \\ ${ }^{\mathrm{b}}$ Laboratory of Thermodynamics, Department of Biochemical and Chemical Engineering, TU Dortmund University, Emil-Figge-Str. 70, 44227 Dortmund, Germany
}

\section*{ARTICLE INFO}

\section*{Keywords:}

Extraction of lithium
Ionic liquid
ePC-SAFT
Quantum chemical calculation
Molecular mechanism

\begin{abstract}
Selective extraction of $\mathrm{Li}^{+}$from high $\mathrm{Mg}^{2+} / \mathrm{Li}^{+}$ratio brines with ionic liquid (IL) based collaborative extractants was investigated by experiments, thermodynamic analyses, and quantum chemical (QC) calculations. Effects of different IL cationic structures and organophosphorus ligands on extraction performances were studied. The results demonstrated that the system 1-(2-hydroxyethyl)-3-methylimidazolium bis(trifluoromethylsulfonyl) imide + trioctyl phosphate ([HOEMIM][Tf2N] + TOP) was considered as the best extractant, with the very high extraction efficiency of $\mathrm{Li}^{+}$( 83.16 \%) and separation selectivity of $\mathrm{Li}^{+} / \mathrm{Mg}^{2+}$ (742.11), which is higher than values reported in literature. The thermodynamic model ePC-SAFT was first extended to quantitatively predict the phase equilibria of the so-called "organic-inorganic complex strong electrolyte system" presented in this work. The molecular-level extraction mechanism was explored by QC calculation, indicating that the strong multi-site intermolecular interactions between $\mathrm{Li}^{+}$and [HOEMIM][Tf2N] + TOP break the $\mathrm{Li}^{+}$hydration. This work provides guidance to rationally design novel IL-based extractants for efficient extraction of $\mathrm{Li}^{+}$.
\end{abstract}

\section*{1. Introduction}

Lithium (Li) as a strategy resource, has been widely applied in the various fields of lithium-ion batteries, aerospace, nuclear energy, medicine and ceramics due to its unique physicochemical properties (Su et al., 2022; Tarascon, 2010; Zhou et al., 2021). In response to the fossil energy crisis, electric vehicles are being developed more and more rapidly around the world, which leads to an unprecedently increase in the demand for lithium resources (Weil et al., 2009; Zhang et al., 2018; Zhou et al., 2019). Li resources mainly exist in hard rocks and salt lake brines, and the latter accounts for over $70 \%$ of Li reserves worldwide (Choubey et al., 2016; Kesler et al., 2012). Currently, more than half of the $\mathrm{Li}_{2} \mathrm{CO}_{3}$ equivalents are obtained from salt lake brines because its cost of preparation is lower than that of extracting $\mathrm{Li}_{2} \mathrm{CO}_{3}$ from hard rock, as the brine is in a liquid state, allowing for a shorter process route (Kuang et al., 2018; Shi et al., 2018). The greatest challenge in extracting $\mathrm{Li}^{+}$ from brines lies in the selective separation of $\mathrm{Li}^{+} / \mathrm{Mg}^{2+}$ since $\mathrm{Mg}^{2+}$ is always accompanied by $\mathrm{Li}^{+}$, often in greater amounts than $\mathrm{Li}^{+}$, and to be separated in an efficient and economical manner. Furthermore, the $\mathrm{Mg}^{2+} / \mathrm{Li}^{+}$ratio is a crucial factor in determining the feasibility and
approach to $\mathrm{Li}^{+}$extraction from brines. Generally, South American brines have a relatively lower $\mathrm{Mg}^{2+} / \mathrm{Li}^{+}$ratio, allowing for the relatively simple and low-cost precipitation method of Li extraction, despite with some inherent drawbacks. $\mathrm{Mg}^{2+} / \mathrm{Li}^{+}$ratios are relatively high in most brines in China and in some brines in North America, and thus, it is a worldwide challenge to selectively extract $\mathrm{Li}^{+}$from the high $\mathrm{Mg}^{2+} / \mathrm{Li}^{+}$ ratio brines (Zhao et al., 2013).

Over the past decades, there are several methods for extracting lithium from brines, e.g., solvent extraction (Shi et al., 2020; Su et al., 2020; Zhou et al., 2012), precipitation (Li et al., 2015; Liu et al., 2018), adsorption (Gu et al., 2018; Jin et al., 2021), membrane separation (Li et al., 2019a; Xu et al., 2021) and electrochemical methods (Zhao et al., 2019). Among them, solvent extraction is the most promising for industrial applications due to its high extraction efficiency of $\mathrm{Li}^{+}$, high separation selectivity of $\mathrm{Li}^{+} / \mathrm{Mg}^{2+}$, simple operation, low cost (Keller et al., 2021; Shi et al., 2019). In terms of solvent extraction, the extractants mainly used are $\beta$-diketones, crown ethers and organophosphorus. $\beta$-Diketones are usually used together with organophosphorus as a synergist extractant ( Li and Binnemans, 2020), and suitable for separation $\mathrm{Li}^{+} / \mathrm{Na}^{+}$instead of $\mathrm{Li}^{+} / \mathrm{Mg}^{2+}$ (Wang et al., 2019; Zhang

\footnotetext{
* Corresponding authors.

E-mail addresses: chenbh@bjut.edu.cn (B. Chen), christoph.held@tu-dortmund.de (C. Held).
}
et al., 2021). Crown ethers are expensive due to the complex synthesis, which are usually used for separation of lithium isotopes ( $\mathrm{Li}^{6} / \mathrm{Li}^{7}$ ), but not suitable for $\mathrm{Li}^{+}$extraction from bines (Cui et al., 2021). It is worth noting that organophosphorus ligands like tributyl phosphate (TBP), trioctylphosphine oxide (TOPO), triisobutyl phosphate (TIBP), and trioctyl phosphate (TOP), are usually considered as a kind of more suitable extractants for $\mathrm{Li}^{+}$extraction from high $\mathrm{Mg}^{2+} / \mathrm{Li}^{+}$ratio brines (Pramanik et al., 2020; Yu et al., 2021b; Yu et al., 2019b). Particularly, TBP as one of the most representative organophosphorus extractants, is usually used in combination with $\mathrm{FeCl}_{3}$ as a co-extractant as well as volatile organic solvents such as dichloromethane, dichloroethane, methyl isobutyl ketone, and kerosene as diluents. However, this inevitably results in the loss of volatile solvents and harm to the environment (Shi et al., 2020; Sun et al., 2023) Furthermore, the addition of $\mathrm{FeCl}_{3}$ leads to the easy formation of a third phase in the extraction system (Zhou et al., 2012) which hinders mass transfer and is not conducive to continuously industrial processes. It is urgently required to develop and search a new class of co-extractants to replace $\mathrm{FeCl}_{3}$ to form diluent-free extraction systems for efficient $\mathrm{Li}^{+}$extraction from high $\mathrm{Mg}^{2+} / \mathrm{Li}^{+}$ratio brines.

Ionic liquids (ILs) as a new type of green solvents and materials have been widely studied in the fields of material synthesis (Chen et al., 2020), chemical reaction (Chen et al., 2011), extractive distillation (Lei et al., 2014), gas separation and absorption (Yu et al., 2021a; Yu et al., 2023; Yu et al., 2020; Yu et al., 2019a) liquid-liquid extraction like extraction of valuable chemicals (Jiang et al., 2021; Qin et al., 2016), aromatics extraction (Navarro et al., 2023; Song et al., 2016), extractive desulfurization (Song et al., 2015), and extractive denitrification (Chen et al., 2014), and extraction of metal ions (Cai et al., 2021; Li et al., 2019b; Zheng et al., 2021), owning to their unique characteristics like negligible vapor pressure, good hydrolysis and thermal stabilities, flammability, and tunable structures. Especially, ILs as extractants can be tuned through the design of cations and anions to improve the extraction efficiency and the selectivity of metal ions (Dupont and Binnemans, 2015; Li et al., 2021; Wang et al., 2017), For example, Rout and Binnemans have studied the solvent extraction of neodymium(III) from nitric acid medium by neutral extractant Cyanex 923 mixed with $\left[\mathrm{Tf}_{2} \mathrm{~N}\right]^{-}$-based ILs. They showed that the extraction efficiency can be related to the solubility of the ILs cation in the aqueous phase (Rout and Binnemans, 2015), So far, the reports regarding the extraction of $\mathrm{Li}^{+}$ using ILs as co-extractants are relatively rare from high $\mathrm{Mg}^{2+} / \mathrm{Li}^{+}$ratio brines especially. The extraction efficiency of $\mathrm{Li}^{+}$with TBP combing with ILs ( $\left[\mathrm{PF}_{6}\right]^{-}$and $\left[\mathrm{Tf}_{2} \mathrm{~N}\right]^{-}$as anions) as mixed extractants was investigated, and it is proved that the extraction mechanism of $\mathrm{Li}^{+}$follows the cation exchange, and the IL anion replaced $\mathrm{FeCl}_{4}^{-}$in the traditional TBP- $\mathrm{FeCl}_{3}$ system to form a complex with $\mathrm{Li}^{+}$and TBP (Shi et al., 2014; Shi et al., 2016; Shi et al., 2017). The researcher synthesized $[\mathrm{BMIM}]_{3}\left[\mathrm{PW}_{2} \mathrm{O}_{40}\right]$ to achieve a highly selective extraction of $\mathrm{Li}^{+}$from the high $\mathrm{Mg}^{2+} / \mathrm{Li}^{+}$ratio brine, and the single-stage extraction efficiency of $\mathrm{Li}^{+}$reached $69.18 \%$ under the optimal conditions (Wang et al., 2018). A comparative study on the extraction effect of $\mathrm{Li}^{+}$with different ILs based on C923 as a neutral ligand was achieved (Cui et al., 2019). However, the above-mentioned studies focus mainly on the experimental performances of the extraction of $\mathrm{Li}^{+}$, and there are relatively few studies on extraction mechanisms. Particularly, there is a lack of systematic understanding on the effect of the structure of ILs on the extraction performance, as well as the thermodynamic behaviors during the extraction process. Even there is almost no research on unraveling the essential characteristics (intermolecular interactions) for separation mechanism of $\mathrm{Li}^{+} / \mathrm{Mg}^{2+}$ at the molecular level related to IL-based extractants. These are crucial for the development and design of new IL-based extractants for highly efficient $\mathrm{Li}^{+} / \mathrm{Mg}^{2+}$ separation. Furthermore, trial-and-error experimental measurements to obtain $\mathrm{Li}^{+} / \mathrm{Mg}^{2+}$ separation performances are very time-consuming. Thus, it is needed to resort to molecular thermodynamic models, which can predict thermodynamic equilibria required for $\mathrm{Li}^{+}$extraction based on the molecular structures of components in mixed systems. This will further
provide a thermodynamic property method to build the equilibrium stage and non-equilibrium stage models in the continuous industrialscale $\mathrm{Li}^{+}$extraction process simulation and design. To the best our knowledge, there are no studies on the thermodynamic modelling of ILrelated $\mathrm{Li}^{+}$extraction systems, due to the complexity of this system, and we refer to this system as "organic-inorganic complex strong electrolyte system". Thus, the ePC-SAFT equation of state (Cameretti et al., 2005) will be employed to describe the thermodynamic phase behavior for $\mathrm{Li}^{+}$ extraction from bines with ILs.

In this work, we selected an aqueous solution of high $\mathrm{Mg}^{2+} / \mathrm{Li}^{+}$ molar ratio of $40: 1$ formulated with $\mathrm{MgCl}_{2}$ and LiCl as the model brine as referred by the previous research (Li and Binnemans, 2021; Shi et al., 2016; Yu et al., 2019b), four [ $\left.\mathrm{Tf}_{2} \mathrm{~N}\right]^{-}$-based ILs, and five organophosphorus ligands (corresponding structures and abbreviations are shown in Fig. 1). This work focuses on addressing the following issues: (i) identifying structure-property relationships between ILs (i.e., different cations) and organophosphorus ligand structures and extraction performance of $\mathrm{Li}^{+}$(extraction efficiency of $\mathrm{Li}^{+}$and separation selectivity of $\mathrm{Li}^{+} / \mathrm{Mg}^{2+}$ ) to determine the best IL-based extractants; (ii) exploring coordination mechanisms for $\mathrm{Li}^{+}$extraction by experimental coupled the advanced spectroscopy techniques; (iii) determining whether ePCSAFT can be extended to predict the so-called phase equilibria of organic-inorganic complex strong electrolyte systems, and analyzing the thermodynamic properties (i.e., Gibbs free energy, enthalpy, and entropic changes) of the extraction process; (iv) revealing the molecularlevel mechanism of $\mathrm{Li}^{+} / \mathrm{Mg}^{2+}$ separation with IL-based extractants by analyzing binding energies, electrostatic potential (ESP) on molecular van der Waals (vdW) surfaces, and independent gradient models (IGM) based on quantum chemical (QC) calculations; and (v) investigating the washing and stripping performances of $\mathrm{Li}^{+}$-rich organic phase, and recycling performance of collaborative extractants. This work provides a systematic understanding of IL-based extractants for Li extraction at the multiscale perspectives ranging from microscopic molecular mechanisms, thermodynamics and to extraction processes, aiming to provide theoretical guidance for developing and designing the novel IL-based extraction systems for high-efficiency $\mathrm{Li}^{+} / \mathrm{Mg}^{2+}$ separation.

\section*{2. Experimental sections}

\subsection*{2.1. Chemicals}
$\mathrm{LiCl}(>99 \mathrm{wt} \%), \mathrm{MgCl}_{2}(>99 \mathrm{wt} \%), \mathrm{NaCl}(99 \mathrm{wt} \%), \mathrm{NaOH}$ and HCl were purchased from Beijing Enokai Technology Co., Ltd., China. Organophosphorus ligands, Cyanex 923 (C923), dibutyl phosphate (DBP), di-(2-ethylhexyl)phosphoric acid (P204), tributyl phosphate (TBP) and trioctyl phosphate (TOP) were obtained from Shanghai Bide Medical Technology Co., Ltd., China. ILs, 1-ethyl-3-methylimidazolium bis(trifluoromethylsulfonyl)imide ([EMIM][ $\left.\mathrm{Tf}_{2} \mathrm{~N}\right]$ ), 1-butyl-3-methylimidazolium bis(trifluoromethylsulfonyl)imide ([BMIM][Tf2N]), 1-octyl-3-methylimidazolium bis(trifluoromethylsulfonyl)imide ([OMIM] $\left[\mathrm{Tf}_{2} \mathrm{~N}\right]$ ), 1-(2-hydroxyethyl)-3-methylimidazolium bis(trifluoromethyl sulfonyl)imide ([HOEMIM][Tf2N]) and 1-butylpyridinium bis (trifluoromethyl-sulfonyl)imide ([BPy][ $\left.\mathrm{Tf}_{2} \mathrm{~N}\right]$ ) were purchased from Shanghai Chengjie Chemical Co.,Ltd., China.

\subsection*{2.2. Extraction experiments}

All extraction experiments were carried out in a centrifuge tube of 50 mL . First, simulated brine $\left(\mathrm{Li}^{+}: 0.766 \mathrm{~g} / \mathrm{L}, \mathrm{Mg}^{2+}: 98.984 \mathrm{~g} / \mathrm{L}\right.$, the high $\mathrm{Mg}^{2+} / \mathrm{Li}^{+}$molar ratio of $40: 1$ ) was prepared by ultra-pure water as the aqueous phase. Then, a certain amount of IL was added to organophosphorus ligands to prepare the organic phase as an extractant system. The aqueous phase of 5 mL was taken and mixed with the organic phase at the given volume ratio. The mixture was shaken and subjected to centrifugal extraction at $300 \mathrm{r} / \mathrm{min}$ for 30 min to achieve extraction equilibrium. Subsequently, it was centrifuged at $8000 \mathrm{r} / \mathrm{min}$ for 5 min to

\begin{figure}
\includegraphics[alt={},max width=\textwidth]{https://cdn.mathpix.com/cropped/91e8b19b-5f41-4aba-88a0-e1b0958e5483-03.jpg?height=936&width=1801&top_left_y=176&top_left_x=139}
\captionsetup{labelformat=empty}
\caption{Fig. 1. Structures and abbreviations for organophosphorus ligands, cations and anions of ILs used in this work.}
\end{figure}
attain complete separation of the aqueous and organic phases. The concentrations of metal ions in the aqueous phase before and after extraction were measured by the inductively coupled plasma emission spectrometer (ICP-OES Optima 8000). The extraction efficiency ( $E_{\mathrm{i}}$ ) of metal ions, the distribution coefficient ( $D_{\mathrm{i}}$ ) of metal ions between the aqueous and organic phases, as well as the separation selectivity of $\mathrm{Li}^{+} / \mathrm{Mg}^{2+}\left(\beta_{\mathrm{Li}^{+} / \mathrm{Mg}^{2+}}\right)$, were determined by the following relation:
$$
\begin{equation*}
E_{i}=\frac{C_{0}-C_{1}}{C_{0}} \times 100 \% \tag{1}
\end{equation*}
$$
$$
\begin{equation*}
D_{i}=\frac{C_{0}-C_{1}}{C_{0}} \times\left(V_{\mathrm{A}} / V_{\mathrm{O}}\right) \tag{2}
\end{equation*}
$$
$$
\begin{equation*}
\beta_{\mathrm{Li}+/ \mathrm{Mg}^{2+}}=\frac{D_{\mathrm{Li}^{+}}}{D_{\mathrm{Mg}^{2+}}} \tag{3}
\end{equation*}
$$

Where $C_{0}$ and $C_{1}$ represents the concentration of given metal ions in the aqueous phase before and after extraction measured experimentally, respectively. $V_{\mathrm{A}}$ and $V_{\mathrm{O}}$ represent the volume of the aqueous and organic phases, respectively.

\subsection*{2.3. Measurements and characterization}

The concentration for cation of ILs in aqueous phase was evaluated on the UV-visible spectroscopic (SHIMADZU UV-2600, Japan) to
measure the number of cations involved in the reaction. The change of pH value in the aqueous phase was monitored via a pH meter (type SN B27310, Thermo Fisher Scientific Inc.) with an accuracy of 0.01 . FTIR spectrometer (Spectrum Two ${ }^{\text {TM }}$ FTIR Spectrometers, PerkinElmer) and ${ }^{1}$ H NMR spectra (Bruker AVANCE HD III) were used to analyze the organic phase after and before extraction.

\section*{3. Thermodynamic modeling and theoretical calculations}

\section*{3.1. ePC-SAFT model}

The ePC-SAFT framework (Cameretti et al., 2005) used in this work calculates the residual Helmholtz energy ( $a^{\text {res }}$ ) by incorporating various Helmholtz energy contributions. It extends the classical PC-SAFT model developed by Gross and Sadowski (Gross and Sadowski, 2001) which considers dispersive perturbations (represented by the contribution $a^{\text {disp }}$ ) and associating perturbations (represented by the contribution $a^{\text {assoc }}$ ) of a hard-chain reference fluid (described by $a^{\text {hc }}$ ). Additionally, ePC-SAFT considers the contribution of electrostatic interactions through the Debye-Hückel theory (Bülow et al., 2021a,b), namely the Helmholtz energy contribution $a^{\mathrm{DH}}$.
$$
\begin{equation*}
a^{\text {res }}=a^{\text {hc }}+a^{\text {disp }}+a^{\text {assoc }}+a^{\mathrm{DH}}\left(\varepsilon_{r}(x)\right) \tag{4}
\end{equation*}
$$

Eq. (4) will be further denoted ePC-SAFT in this work. The dielectric constant $\varepsilon_{r}(x)$ of the medium is utilized in the calculation of $a^{\mathrm{DH}}$. The

\begin{table}
\captionsetup{labelformat=empty}
\caption{Table 1
Dielectric constants for all components applied in this work.}
\begin{tabular}{|l|l|l|}
\hline Component & Dielectric constant/ $\mathrm{C} \mathrm{Vm}^{-1}$ & Ref. \\
\hline Water & $-105.2 \ln (T[\mathrm{~K}])+677.480$ & (Ascani and Held, 2021) \\
\hline TOP & 11 & Set to constant number in this work \\
\hline [HOEMIM][ $\mathrm{Tf}_{2} \mathrm{~N}$ ] & 11 & a) \\
\hline $\mathrm{Li}^{+}$ & 8 & b) \\
\hline $\mathrm{Mg}^{2+}$ & 8 & b) \\
\hline $\mathrm{Cl}^{-}$ & 8 & b) \\
\hline
\end{tabular}
\end{table}

For the IL [HOEMIM] $\left[\mathrm{Tf}_{2} \mathrm{~N}\right]$, a relative dielectric constant of 11 was chosen according to a previous work (Bülow et al., 2019). All salts were modeled with a similar dielectric constant that is a mean of available experimental data (Andeen et al., 1970).
required pure-component values $\varepsilon_{r, j}$ for the solvents and ions used in the calculation are summarized in Table 1 (Andeen et al., 1970; Ascani and Held, 2021; Bülow et al., 2019). All ions were assigned the same value of $\varepsilon_{r, i n n}=8$. The dielectric constant of the medium $\varepsilon_{r}(x)$ was then calculated from a linear combination in the solvent mass fraction and the ion mole fraction:
$$
\begin{equation*}
\varepsilon_{r}(x)=\left(\sum_{j=1}^{N^{\text {solv }}} \varepsilon_{r, j} w_{j}^{\text {solv }}\right) x_{\text {solv }}+\sum_{j=1}^{N^{\text {ion }}} \varepsilon_{r, j} w_{j}^{\text {ion }} \tag{5}
\end{equation*}
$$

In this equation, $N^{\text {solv }}$ expresses the total number of components in the salt-free solvent mixture and $N^{\text {ion }}$ the total number of ions. $x_{j}$ and $w_{j}$ denote the mole fraction and the mass fraction of component $j . w_{j}{ }^{\text {solv }}$ is the mass fraction of solvent $j$ in the salt-free solvent mixture, while $x_{\text {solv }}$ is the sum of the mole fraction of all solvents in the overall mixture. Associating components such as water are assigned five pure-component parameters, namely the segment number $m_{i}^{\text {seg }}$, segment diameter $\sigma_{\mathrm{i}}$, dispersion-energy parameter $u_{i} / k_{\mathrm{b}}$ association-energy parameter $\varepsilon^{A_{i} B_{i}} / k_{b}$, and association-volume parameter $\kappa^{A_{i} B_{i}}$. For the ions considered in this work, two pure-component parameters, $m_{i}^{\text {seg }}$ and $\sigma_{\text {i }}$, are utilized. Additionally, the 2 B association scheme is used for water, TOP, IL and $\mathrm{Li}^{+}$during calculating associating interactions (Cameretti and Sadowski, 2008). Table 2 presents the pure-component parameters employed in this study (Held et al., 2008; Held et al., 2014).

The modeling of mixtures necessitates the utilization of the Berthelot-Lorentz (Berthelot, 1898; Lorentz, 1881) and WolbachSandler (Wolbach and Sandler, 1997) combining rules:
$$
\begin{equation*}
\sigma_{i j}=\frac{1}{2}\left(\sigma_{i}+\sigma_{j}\right) \tag{6}
\end{equation*}
$$
$$
\begin{equation*}
u_{i j}=\sqrt{u_{i} u_{j}\left(1-k_{i j}\right)} \tag{7}
\end{equation*}
$$
$$
\begin{equation*}
\varepsilon^{A_{i} B_{j}}=\frac{1}{2}\left(\varepsilon^{A_{i} B_{i}}+\varepsilon^{A_{j} B_{j}}\right)\left(1-k_{i j}^{h b}\right) \tag{8}
\end{equation*}
$$
$$
\begin{equation*}
\kappa^{A_{i} B_{j}}=\sqrt{\kappa^{A_{i} B_{i}} \kappa^{A_{j} B_{j}}}\left(\frac{\sqrt{\sigma_{i} \sigma_{j}}}{1 / 2\left(\sigma_{i}+\sigma_{j}\right)}\right)^{3} \tag{9}
\end{equation*}
$$

In Eq. (8), the binary interaction parameter $k_{i j}$ is introduced, allowing for the modification of the dispersion energy of the pair $i j$. Binary interaction parameters were required for the following pairs: ion-ion, ionwater, ion-TOP and TOP-water. Particularly, the ePC-SAFT model has been shown to be predictive for phase equilibria in complex mixtures of up to 6 components (Ascani et al., 2022). That is, fit parameters usually have not to be adjusted to the multicomponent mixture of interest, but only to the subsystems. The parameter for the IL-water pair was fitted to the phase equilibrium between IL and water, while the remaining parameters were fitted to one experimental data point. In addition, the binary interaction parameter $k_{i j}^{h b}$ for the association parameter of mixtures for the TOP- $\mathrm{Li}^{+}$system was adjusted to account for the coordination of $\mathrm{Li}^{+}$to the $\mathrm{P}=\mathrm{O}$ group of TOP. The binary interaction parameters used in this study are listed in Table 3 (Chapeaux et al., 2007; Held et al., 2014). For $\mathrm{Mg}^{2+}$ binary interaction parameters were not applied. In this

\begin{table}
\captionsetup{labelformat=empty}
\caption{Table 3
Binary interaction parameters $k_{i j}$ used in this work.}
\begin{tabular}{|l|l|l|}
\hline Parameters & $k_{\mathrm{ij}}$ & $k_{i j}^{h b}$ \\
\hline Water - [HOEMIM][ $\mathrm{Tf}_{2} \mathrm{~N}$ ] & $0.007{ }^{\text {a) }}$ & - \\
\hline Li+-TOP & 0.3 & 0.3 \\
\hline Water - TOP & 1 & - \\
\hline $\mathrm{Li}^{+}$- [HOEMIM] [Tf $\left.{ }_{2} \mathrm{~N}\right]^{-}$ & 1 & 1 \\
\hline $\mathrm{Li}^{+}$- Water & - & 1 \\
\hline $\mathrm{Li}^{+}-\mathrm{Cl}^{-}$ & 0.669 from reference (Held et al., 2014) & - \\
\hline
\end{tabular}
\end{table}
a) Fitted to binary data water- $[\mathrm{HOEMIM}]\left[\mathrm{Tf}_{2} \mathrm{~N}\right]$ from the reference (Chapeaux et al., 2007).
study, we found that only one fit parameter was required to match the data quantitatively correct ( $\mathrm{Li}^{+}$-TOP, see Table 3), while other interaction parameters were not used or were used to switch off the dispersion interaction (i.e., value equal to one, see Table 3), or fitted to the subsystem (e.g., $\mathrm{Li}^{+}-\mathrm{Cl}^{-}$, or water-IL).

\subsection*{3.2. Quantum chemical (QC) calculations}

The quantum chemical (QC) calculation based on the density function theory (DFT) was performed to identify the molecular-level extraction separation mechanism. The initial geometrical configurations regarding metal ion hydrates, IL, TOP, and the complex of metal ions and IL/TOP were optimized by the Gaussian 09 program (Version D.01) (Frisch et al., 2009) at the B3LYP/6-31+G(d, p) level (M. et al., 2005), with the DFT-D3(BJ) dispersion correction (Grimme et al., 2011). The frequency calculation was carried out to ensure the stability of the obtained geometry with the lowest energy. Based on the above configuration, the binding energies of complexes ( $\Delta E_{\text {Complex }}$ ) between clusters were calculated by:
$$
\begin{equation*}
\Delta E_{\text {Complex }}=E_{A-B}-E_{A}-E_{B} \tag{10}
\end{equation*}
$$

Where $E_{\mathrm{A}-\mathrm{B}}, E_{\mathrm{A}}$ and $E_{\mathrm{B}}$ represent the electron energy (in kcal/mol) of complexes of A-B, cluster A and cluster B, respectively. Moreover, the ESP analysis of molecular vdW surfaces and IGM analysis were performed by the Multiwfn tool (version 3.8) as developed by Lu and Chen (2012) to visually identify the type, strength and location of intermolecular interactions formed and to deeply explore extraction mechanism.

\section*{4. Results and discussion}

\subsection*{4.1. Structural effects of ILs and organophosphorus ligands on extraction performance}

The extraction performance for the five different structural organophosphorus ligands (i.e., C923, DBP, P204, TBP and TOP) containing the $\mathrm{P}=\mathrm{O}$ group were experimentally measured and compared to explore the structure-property relationship under the specific IL $[\mathrm{BMIM}]\left[\mathrm{Tf}_{2} \mathrm{~N}\right]$ at the given operation conditions (e.g., molar ratio of $\mathrm{Mg}^{2+}$ to $\mathrm{Li}^{+}$in the aqueous phase $\mathrm{Mg}^{2+} / \mathrm{Li}^{+}=40: 1$, the volume ratio for organic to aqueous phase $\mathrm{O} / \mathrm{A}=3: 1, \mathrm{pH}=7$, and IL concentration in organic phase $C_{\mathrm{IL}}=0.1 \mathrm{~mol} / \mathrm{L}$ ), and the results are shown in Fig. 2 (see Table S1 in

\begin{table}
\captionsetup{labelformat=empty}
\caption{Table 2
The ePC-SAFT pure-component parameters used in this work.}
\begin{tabular}{|l|l|l|l|l|l|l|l|}
\hline Component & $m_{i}^{\text {seg }}$ & $\sigma_{\mathrm{i}}(\AA)$ & $u_{i} / k_{\mathrm{b}}(\mathrm{K})$ & $\varepsilon^{A_{i} B_{i}} / k_{b}$ (K) & $\kappa^{A_{i} B_{i}}$ & Association scheme & Refs \\
\hline Water & 1.2047 & a) & 353.95 & 2425.7 & 0.0451 & 2B & (Held et al., 2008) \\
\hline TOP & 4.2032 & 5.4506 & 280.4777 & 6393.5 & 0.0001 & 2B & This work \\
\hline [HOEMIM][ $\mathrm{Tf}_{2} \mathrm{~N}$ ] & 4.073 & 4.6432 & 434.6120 & 5000 & 0.1 & 2B & This work \\
\hline $\mathrm{Li}^{+}$ & 1 & 2.8449 & 360.00 & 0 & 100 & 2B & This work \\
\hline $\mathrm{Mg}^{2+}$ & 1 & 3.1327 & 1500 & - & - & - & (Held et al., 2014) \\
\hline $\mathrm{Cl}^{-}$ & 1 & 2.7560 & 170.00 & - & - & - & (Held et al., 2014) \\
\hline
\end{tabular}
\end{table}
a) $\sigma=2.7927+\left(10.11 \cdot e^{-0.01775 T}-1.417 \cdot e^{-0.01146 T}\right)$ with $T$ in K .

\begin{figure}
\includegraphics[alt={},max width=\textwidth]{https://cdn.mathpix.com/cropped/91e8b19b-5f41-4aba-88a0-e1b0958e5483-05.jpg?height=513&width=1405&top_left_y=183&top_left_x=331}
\captionsetup{labelformat=empty}
\caption{Fig. 2. Effects of different ligands on the extraction efficiency (a) and separation performances ( $D_{\mathrm{Li}^{+}}$and $\beta_{\mathrm{Li}^{+} / \mathrm{Mg}^{2+}}$ ) of $\mathrm{Li}^{+} / \mathrm{Mg}^{2+}$ (b) in the presence of [BMIM][Tf ${ }_{2} \mathrm{~N}$ ] (molar ratio of $\mathrm{Mg}^{2+}$ to $\mathrm{Li}^{+}$in the aqueous phase $\mathrm{Mg}^{2+} / \mathrm{Li}^{+}=40: 1$; the volume ratio for organic to aqueous phases $\mathrm{O} / \mathrm{A}=3: 1 ; \mathrm{pH}=7$; and IL concentration in organic phase $C_{\mathrm{IL}}=0.1 \mathrm{~mol} / \mathrm{L}$ ) at 298.15 K .}
\end{figure}

Supplementary material for more details). It can be seen form Fig. 2a that the $\mathrm{Li}^{+}$extraction efficiencies for ligands TBP (72.27 \%) and TOP (69.81\%) are very close to each other, which are far higher than those for any of others. Meanwhile, the similar phenomenon is also presented in the distribution coefficients of $\mathrm{Li}^{+}$(Fig. 2b). Noteworthily, the separation selectivity of $\mathrm{Li}^{+} / \mathrm{Mg}^{2+}$ of TOP is up to the highest 199.63 among all the ligands and much higher than that of TBP. This is because that TOP has the longer alkyl chains with the branch chains, which leads to greater steric hindrance to combine $\mathrm{Mg}^{2+}$ with larger radius. Therefore, TOP is a preferred candidate ligand to perform the experiment below in this work.

The structural effect of hydrophobic ILs containing different cations under the same anion $\left[\mathrm{Tf}_{2} \mathrm{~N}\right]^{-}$on the extraction performance of $\mathrm{Li}^{+}$with TOP as the ligand at the given operating conditions ( $\mathrm{Mg}^{2+} / \mathrm{Li}^{+}=40: 1$; $\mathrm{O} / \mathrm{A}=3: 1 ; \mathrm{pH}=7$ and $C_{\mathrm{IL}}=0.1 \mathrm{~mol} / \mathrm{L}$ ), and the results are shown in Fig. 3 (see Table S2 in Supplementary material for more details). Fig. 3a shows that the extraction efficiency of $\mathrm{Li}^{+}$in different $\left[\mathrm{Tf}_{2} \mathrm{~N}\right]^{-}$-based ILs follows the order of [OMIM][Tf ${ }_{2} \mathrm{~N}$ ] (39.51\%) < [BMIM][ $\mathrm{Tf}_{2} \mathrm{~N}$ ] (69.6 $\%)<[$ EMIM $]\left[\mathrm{Tf}_{2} \mathrm{~N}\right](73.69 \%)$, which evidently indicates that the IL with the longer alkyl chain is unfavorable for extracting $\mathrm{Li}^{+}$ions. This may be attributed to two reasons as follows: the longer alkyl chain to (i) increase the IL hydrophobicity resulting in weakening the exchange capacity with $\mathrm{Li}^{+}$, and (ii) increase the IL viscosity leading to weaken the $\mathrm{Li}^{+}$mass transfer from aqueous to organic phases. It should be noted that
the cationic structural effect on the extraction efficiency and distribution coefficient of $\mathrm{Li}^{+}$follows the order of [OMIM][ $\left.\mathrm{Tf}_{2} \mathrm{~N}\right]<[\mathrm{BMIM}]\left[\mathrm{Tf}_{2} \mathrm{~N}\right] <[$ BPy $]\left[\mathrm{Tf}_{2} \mathrm{~N}\right]<[$ EMIM $]\left[\mathrm{Tf}_{2} \mathrm{~N}\right]<[$ HOEMIM $]\left[\mathrm{Tf}_{2} \mathrm{~N}\right]$ (Fig. 3a). Generally, it seems to be not obvious for the cationic structural effect on the extraction performance of $\mathrm{Mg}^{2+}$. This leads to the result that the cationic structural effect on the separation selectivity of $\mathrm{Li}^{+} / \mathrm{Mg}^{2+}$ is similar to the effects on extraction efficiency and distribution coefficient of $\mathrm{Li}^{+}$ following the same trend (except for the inconsistent phenomenon of $[\mathrm{BMIM}]\left[\mathrm{Tf}_{2} \mathrm{~N}\right]>[\mathrm{BPy}]\left[\mathrm{Tf}_{2} \mathrm{~N}\right]$ ) (Fig. 3b). Thus, the cation with the shorter alkyl chain is preferred. Particularly, the extraction efficiency of $\mathrm{Li}^{+}$can be greatly enhanced from $73.69 \%$ to $80.27 \%$, and separation selectivity of $\mathrm{Li}^{+} / \mathrm{Mg}^{2+}$ much increased from 206.88 to 709.92 when a hydroxyl group is introduced into the alkyl chain of the cation of [EMIM][ $\mathrm{Tf}_{2} \mathrm{~N}$ ] to become [HOEMIM][ $\mathrm{Tf}_{2} \mathrm{~N}$ ]. This is resulted from the fact that the presence of hydroxyl group can intensify the affinity of cations to water, which made it easier to enter the aqueous phase for exchange with $\mathrm{Li}^{+}$. Therefore, [HOEMIM][ $\mathrm{Tf}_{2} \mathrm{~N}$ ] is selected as the IL candidate for the subsequent study. Meanwhile, The IL [HOEMIM] $\left[\mathrm{Tf}_{2} \mathrm{~N}\right]$ exhibits good thermal stability (the onset decomposition temperature is about $400^{\circ} \mathrm{C}$ ), and the relatively low viscosity ( $73 \mathrm{mPa} \cdot \mathrm{s}$ at $30^{\circ} \mathrm{C}$ ) (see Fig. S1 in Supplementary material), which can be further decrease due to the presence of water during $\mathrm{Li}^{+}$extraction processes.

\begin{figure}
\includegraphics[alt={},max width=\textwidth]{https://cdn.mathpix.com/cropped/91e8b19b-5f41-4aba-88a0-e1b0958e5483-05.jpg?height=662&width=1401&top_left_y=1832&top_left_x=335}
\captionsetup{labelformat=empty}
\caption{Fig. 3. Effects of different ILs on extraction efficiency (a) and separation performances ( $D_{\mathrm{Li}^{+}}$and $\beta_{\mathrm{Li}+} / \mathrm{Mg}^{2+}$ ) of $\mathrm{Li}^{+} / \mathrm{Mg}^{2+}$ (b) with TOP as a neutral ligand ( $\mathrm{Mg}^{2+} / \mathrm{Li}^{+}=$ 40:1; $\mathrm{O} / \mathrm{A}=3: 1 ; \mathrm{pH}=7$ and $C_{\mathrm{IL}}=0.1 \mathrm{~mol} / \mathrm{L}$ ) at 298.15 K .}
\end{figure}

\subsection*{4.2. Optimizing operating conditions for extraction processes}

The operating conditions (i.e., pH values, IL concentration $C_{\mathrm{IL}}$ and O/A ratio) for extraction of $\mathrm{Li}^{+}$were optimized, and the results are shown in Fig. 4 (see Tables S3-S5 in Supplementary material for more details). As shown in Fig. 4a and b, when the pH value increased from 1 to 4, the extraction efficiency of $\mathrm{Li}^{+}$shows a slight upward trend from $75.82 \%$ to $83.44 \%$, while there will be no significant change when the pH value is higher than 4 . The extraction efficiency of $\mathrm{Mg}^{2+}$ has a slight decrease trend over the range from 1 to 7 . When pH value is greater than 7, there are white $\mathrm{Mg}(\mathrm{OH})_{2}$ precipitates from the aqueous phase, so the pH value has to be controlled to be within 7 . The strong acidity environment is unbeneficial to extraction $\mathrm{Li}^{+}$process, which is attributed to the fact that the competitive effect between $\mathrm{H}^{+}$and $\mathrm{Li}^{+}$result in the lower separation selectivity of $\mathrm{Li}^{+} / \mathrm{Mg}^{2+}$. Particularly, the higher pH value is in favor of higher separation selectivity, which is due to gradually decreasing distribution coefficient of $\mathrm{Mg}^{2+}$ with increasing the pH value. Moreover, separation selectivity of $\mathrm{Li}^{+} / \mathrm{Mg}^{2+}$ increases from 67.85 to 709.01 rapidly as the pH value rises from 4 to 7 . Accordingly, the optimal pH value for the extraction system is 7 .

The effect of IL concentration on the extraction performance is shown in Fig. 4c and d. It is found that using only TOP has almost no
extraction ability for $\mathrm{Li}^{+}$and $\mathrm{Mg}^{2+}$. When the [HOEMIM][ $\mathrm{Tf}_{2} \mathrm{~N}$ ] concentration increases from 0 to $0.09 \mathrm{~mol} / \mathrm{L}$, the extraction efficiency of $\mathrm{Li}^{+}$sharply increases with up to $80 \%$ for the concentration of $0.09 \mathrm{~mol} /$ L. For IL concentrations greater than $0.09 \mathrm{~mol} / \mathrm{L}$, the extraction efficiency of $\mathrm{Li}^{+}$remains almost a constant. However, the separation selectivity of $\mathrm{Li}^{+} / \mathrm{Mg}^{2+}$ decreases for IL concentrations greater than $0.09 \mathrm{~mol} / \mathrm{L}$. This is because of the fact that the saturated extraction of $\mathrm{Li}^{+}$will be reached under the IL concentration of $0.09 \mathrm{~mol} / \mathrm{L}$, which leads to an excess of IL increasing the ability to extract $\mathrm{Mg}^{2+}$ and elevating the distribution coefficient of $\mathrm{Mg}^{2+}$. Therefore, the [HOEMIM] $\left[\mathrm{Tf}_{2} \mathrm{~N}\right]$ concentration of $0.09 \mathrm{~mol} / \mathrm{L}$ is selected an optimal value.

The effect of $\mathrm{O} / \mathrm{A}$ ratios on the extraction performance are shown in Fig. 4e and f. It is found from Fig. 4e that when the O/A ratio increases from $1: 1$ to $3: 1$, the extraction efficiency of $\mathrm{Li}^{+}$presents an obvious upward trend. Subsequently, when the O/A ratio increases from 3:1 to 6:1, the increase in extraction efficiency of $\mathrm{Li}^{+}$tends to level off. Fig. 4f shows that separation selectivity of $\mathrm{Li}^{+} / \mathrm{Mg}^{2+}$ increases significantly when the $\mathrm{O} / \mathrm{A}$ ratio increases from $1: 1$ to $3: 1$, with a maximum value of 742.11, while it will decrease when the O/A ratio is increased to values greater than 3:1. This is because the extraction of $\mathrm{Li}^{+}$by $\mathrm{IL}+$ TOP reached saturation at this ratio, and continuing to increase the amount of $\mathrm{IL}+$ TOP enhanced the extraction capacity of $\mathrm{Mg}^{2+}$. Thus, O/A ratio

\begin{figure}
\includegraphics[alt={},max width=\textwidth]{https://cdn.mathpix.com/cropped/91e8b19b-5f41-4aba-88a0-e1b0958e5483-06.jpg?height=1471&width=1239&top_left_y=1029&top_left_x=413}
\captionsetup{labelformat=empty}
\caption{Fig. 4. Effects of pH values $\left(\mathrm{a}, \mathrm{b} ; \mathrm{Mg}^{2+} / \mathrm{Li}^{+}=40: 1, \mathrm{O} / \mathrm{A}=3: 1\right.$ and $\left.C_{\mathrm{IL}}=0.1 \mathrm{~mol} / \mathrm{L}\right)$, IL concentrations $C_{\mathrm{IL}}\left(\mathrm{c}, \mathrm{d}^{-} \mathrm{Mg}^{2+} / \mathrm{Li}^{+}=40: 1, \mathrm{O} / \mathrm{A}=3: 1\right.$ and $\left.\mathrm{pH}=7\right)$ and $\mathrm{O} / \mathrm{A}$ ratios ( $\mathrm{e}, \mathrm{f} ; \mathrm{Mg}^{2+} / \mathrm{Li}^{+}=40: 1 ;$ and $\mathrm{pH}=7$ and $C_{\mathrm{IL}}=0.09 \mathrm{~mol} / \mathrm{L}$ ) on extraction performances with [HOEMIM][Tf ${ }_{2} \mathrm{~N}$ ] + TOP as extractant system at 298.15 K .}
\end{figure}

\begin{table}
\captionsetup{labelformat=empty}
\caption{Table 4
Comparison of extraction performances of $\mathrm{Li}^{+}$with different extractant systems at 298.15 K .}
\begin{tabular}{|l|l|l|l|l|l|l|l|}
\hline Ligands & ILs & $\mathrm{Mg}^{2+} / \mathrm{Li}^{+}$molar ratios & One-stage $\mathrm{Li}^{+}$extraction efficiency (\%) & $\beta_{\mathrm{Li}^{+}} / \mathrm{Mg}^{2+}$ & O/A & pH & Refs. \\
\hline TBP & [HOEMIM][ $\mathrm{Tf}_{2} \mathrm{~N}$ ] & 40 & 80 & 170 & 1:1 & 7 & (Li and Binnemans, 2021) \\
\hline TBP & [BMIM][ $\mathrm{PF}_{6}$ ] & 12.58 & 90.93 (92.5) & - (-) & 2:1(3:1) & 7 & (Shi et al., 2014) \\
\hline TBP & [BMIM][ $\mathrm{Tf}_{2} \mathrm{~N}$ ] & 13.05 & 92.37 (94.61) & 403.33 (-) & 2:1(3:1) & 7 & (Shi et al., 2016) \\
\hline TBP & [BMIM] ${ }_{3} \mathrm{PW}_{12} \mathrm{O}_{40}$ & 78.33 & 69.18 (72) & 283.06 (210) & 1:1(3:1) & 7 & (Wang et al., 2018) \\
\hline TIBP & [EMIM][ $\mathrm{Tf}_{2} \mathrm{~N}$ ] & 9.54 & 83.71 (91.66) & 71.24 (138.61) & 1:1(2.5:1) & 7 & (Gao et al., 2015) \\
\hline TBP & [MMIM] ${ }_{4} \mathrm{SiW}_{12} \mathrm{O}_{40}$ & 78.33 & 60 (65) & 35 (125) & 1:1(3:1) & 7 & (Wang et al., 2023) \\
\hline TOP & [HOEMIM][ $\mathrm{Tf}_{2} \mathrm{~N}$ ] & 40 & 83.16 & 742.11 & 3:1 & 7 & This work \\
\hline
\end{tabular}
\end{table}

Note: The O/A in the brackets corresponds to the $\mathrm{Li}^{+}$extraction efficiency and $\beta_{\mathrm{Li}^{+}} / \mathrm{Mg}^{2+}$.
of 3:1 is a suitable selection. In addition, a comparison of extraction performance with different extractant systems collected from literature are shown in Table 4 (Gao et al., 2015; Li and Binnemans, 2021; Shi et al., 2014; Shi et al., 2016; Wang et al., 2023; Wang et al., 2018). It is found that the [HOEMIM] [ $\mathrm{Tf}_{2} \mathrm{~N}$ ] + TOP extractant system in this work has the very highest max one-stage extraction efficiency of $\mathrm{Li}^{+}$(83.16), and separation selectivity of $\mathrm{Mg}^{2+} / \mathrm{Li}^{+}$(742.11) is highest among all extractant systems collected from literature. Besides, when compared with the extractant systems of [HOEMIM][Tf2N] + TBP, TBP$[\mathrm{BMIM}]_{3} \mathrm{PW}_{12} \mathrm{O}_{40}+\mathrm{TBP}$ and $[\mathrm{MMIM}]_{4} \mathrm{SiW}_{12} \mathrm{O}_{40}$, + TBP, the extraction performance using [HOEMIM][Tf ${ }_{2} \mathrm{~N}$ ] + TOP is significantly improved. In summary, the extractant system of [HOEMIM][ $\left.\mathrm{Tf}_{2} \mathrm{~N}\right]+$ TOP proposed in this study had great selective extraction ability for $\mathrm{Li}^{+}$.

The multistage extraction experiment was carried out to validate the extraction performance of $\mathrm{Li}^{+}$with the extractant of $\mathrm{IL}+$ TOP in the real industrial application at the above-mentioned optimized conditions, and the results are shown in Fig. S2 in Supplementary material. It can be seen that the extraction efficiency of $\mathrm{Li}^{+}$reached $99.07 \%$, and $\mathrm{Li}^{+}$concentration in the aqueous phase decreased from the original $760.6 \mathrm{mg} / \mathrm{L}$ to $7.1 \mathrm{mg} / \mathrm{L}$ after five stage extractions. The results show that the IL + TOP composition is a promising extractant applied for extraction of $\mathrm{Li}^{+}$from the brine with high $\mathrm{Mg}^{2+} / \mathrm{Li}^{+}$ratio.

\subsection*{4.3. Analyses of extraction mechanisms via spectroscopy techniques}

\subsection*{4.3.1. Analyses of ${ }^{1}$ H NMR, UV-visible, and FTIR}

Fig. $5 a$ and $b$ show that the ${ }^{1} \mathrm{H}$ NMR spectra of the organic phase does not change significantly before and after extraction. The chemical shift corresponding to the hydroxyl peak of [HOEMIM] ${ }^{+}$at 3.81 ppm , and the peak intensity is weakened after extraction, indicating that [HOEMIM] ${ }^{+}$ may be exchanged into the aqueous phase, resulting in a decrease in [HOEMIM] ${ }^{+}$concentration. In addition, the UV-visible spectrum was used to analyze the spectrum of raffinate (i.e., the aqueous phase) at 190 $\sim 250 \mathrm{~nm}$ (see Fig. S3), there is an obvious absorbance peak of [HOEMIM] ${ }^{+}$presented at 209 nm . It is found that the cation [HOEMIM] ${ }^{+}$species in the aqueous phase after extraction are determined although when the metal ions are not added into the aqueous phase. This is due to the certain solubility of IL in the aqueous solution with the [HOEMIM] ${ }^{+}$concentration of $0.7135 \mathrm{~g} / \mathrm{L}$ measured quantitatively. However, the absorbance strength of the aqueous phase with metal ions (the [HOEMIM] ${ }^{+}$concentration of $11.7477 \mathrm{~g} / \mathrm{L}$ measured quantitatively) is significantly higher than that without metal ions. This completely the presence of metal ions in the aqueous phase can enhance the [HOEMIM] ${ }^{+}$entered from the oil phase into the aqueous phase, which is the significant evidence for occurring the cationic exchange mechanism. Fig. 5c showed the FTIR spectra of [HOEMIM][ $\mathrm{Tf}_{2} \mathrm{~N}$ ], TOP and their mixing. It is found that all peaks of [HOEMIM][ $\mathrm{Tf}_{2} \mathrm{~N}$ ] and TOP still existed after mixing, indicating that [HOEMIM] and TOP were physically mixed. Fig. 5d and e show the characteristic information of FTIR spectra for the organic phase before and after extraction, and the wavelengths at $2959 \mathrm{~cm}^{-1}$ and $2872 \mathrm{~cm}^{-1}$ corresponded to the asymmetric and symmetric tensile vibration of $\mathrm{CH}_{3}$ group, respectively.

Moreover, the peaks at $1463 \mathrm{~cm}^{-1}$ and $1380 \mathrm{~cm}^{-1}$ are assigned to the asymmetric and symmetric bending vibrations of $\mathrm{CH}_{3}$ respectively, $1014 \mathrm{~cm}^{-1}$ is attributed to the tensile vibration of P-O-C, and all of these peaks are not shifted. However, the peak at $1269 \mathrm{~cm}^{-1}$ representing $\mathrm{P}=\mathrm{O}$ in TOP is shifted to $1259 \mathrm{~cm}^{-1}$ after extraction. This is because that the $\mathrm{P}=\mathrm{O}$ vibration frequency changed is resulted by combining TOP with $\mathrm{Li}^{+}$, representing the strong interaction is formed between TOP with $\mathrm{Li}^{+}$. The peaks at $1352 \mathrm{~cm}^{-1}, 1189 \mathrm{~cm}^{-1}$, and $1051 \mathrm{~cm}^{-1}$ are corresponding to characteristic information of $\mathrm{S}=\mathrm{O}, \mathrm{C}-\mathrm{F}, \mathrm{S}-\mathrm{N}-\mathrm{S}$ in $\left[\mathrm{Tf}_{2} \mathrm{~N}\right]^{-}$, respectively, and all of them have not obvious shift, indicating the anion of IL still stays in the organic phase after extraction instead of the aqueous phase. It is found from Fig. S4 in Supplementary material that the $\mathrm{C}=\mathrm{C}$ stretching vibration in imidazolium ring on [HOEMIM] ${ }^{+}$ at $1572 \mathrm{~cm}^{-1}$ disappeared after extraction, which is due to the ion exchange happened between [HOEMIM] ${ }^{+}$and $\mathrm{Li}^{+}$in the aqueous phase, resulting in the decrease of [ HOEMIM$]^{+}$concentration in organic phase. The results are in the agreement with those obtained by UV-visible and ${ }^{1} \mathrm{H}$ NMR spectra.

\subsection*{4.3.2. Extraction equilibrium parameters}

According to the analyses of FTIR, ${ }^{1} \mathrm{H}$ NMR and UV-visible spectra, the extraction process can be treated as the following relation:
$$
\begin{align*}
{[\mathrm{HOEMIM}]\left[\mathrm{Tf}_{2} \mathrm{~N}\right]_{\text {org }} } & +\mathrm{Li}_{\mathrm{aq}}^{+}+n \mathrm{TOP}_{\text {org }} \rightarrow[\mathrm{HOEMIM}]_{\mathrm{aq}}^{+}+[\mathrm{Li} \cdot n \mathrm{TOP}]_{\text {org }}^{+}  \tag{11}\\
& +\left[\mathrm{Tf}_{2} \mathrm{~N}\right]_{\text {org }}^{-}
\end{align*}
$$
where subscripts "aq" and "org" represent the aqueous and organic phases respectively, and $n$ represents the number of TOP molecules. The extraction equilibrium constant ( $K$ ) can be expressed as:
$$
\begin{equation*}
K=\frac{D_{\mathrm{Li}^{+}} \cdot[\mathrm{HOEMIM}]_{\mathrm{aq}}^{+}}{[\mathrm{HOEMIM}]_{\mathrm{aq}}^{+}[\mathrm{HOEMIM}]_{\mathrm{org}}^{n}}, D_{\mathrm{Li}^{+}}=\frac{[\mathrm{Li} \cdot n \mathrm{TOP}]_{\mathrm{org}}^{+}}{[\mathrm{Li}]_{\mathrm{aq}}^{+}} \tag{12}
\end{equation*}
$$

Since, it can be further represented by as follows:
$$
\begin{equation*}
\log D_{\mathrm{Li}^{+}}+\log [\mathrm{HOEMIM}]_{\mathrm{aq}}^{+}-\log [\mathrm{HOEMIM}]_{\mathrm{org}}^{+}=n \log [\mathrm{TOP}]_{\mathrm{org}}+\log K \tag{13}
\end{equation*}
$$

The $n=1$ was determined by calculating the slop of curve (see Fig. S5 in Supplementary material), indicating that one TOP molecule can bind with one $\mathrm{Li}^{+}$.

\subsection*{4.4. Extraction thermodynamic behaviors}

\subsection*{4.4.1. Phase equilibria behavior predicted by the ePC-SAFT model}

In this study, phase equilibria were calculated using the isofugacity criterion, which states that at thermodynamic equilibrium, the fugacity of each component is equal in all phases.
$$
\begin{equation*}
f_{i}^{1}=f_{i}^{2}=\cdots=f_{i}^{\pi} \tag{14}
\end{equation*}
$$

For calculation of the fugacities, the so-called " $\varphi-\varphi$ " criterion was used. The fugacity of each component in both phases was expressed using the respective fugacity coefficient as below:

\begin{figure}
\includegraphics[alt={},max width=\textwidth]{https://cdn.mathpix.com/cropped/91e8b19b-5f41-4aba-88a0-e1b0958e5483-08.jpg?height=1920&width=1439&top_left_y=180&top_left_x=314}
\captionsetup{labelformat=empty}
\caption{Fig. 5. ${ }^{1} \mathrm{H}$ NMR of ( $\mathrm{a}, \mathrm{b}$; b is the partial magnification of a ) and FTIR ( c -e; e is the partial magnification of d ) spectra for organic phase before and after extraction (10 $\mathrm{v} \% \mathrm{IL}+90 \mathrm{v} \% \mathrm{TOP}$, and $\mathrm{O} / \mathrm{A}=3: 1$ ) at 298.15 K .}
\end{figure}
$$
\begin{equation*}
\varphi_{i}^{L 1} x_{i}^{L 1}=\varphi_{i}^{L 2} x_{i}^{L 2} \tag{15}
\end{equation*}
$$

The fugacity coefficients $\varphi_{i}$ can be derived with ePC-SAFT from the residual Helmholtz energy $a^{\text {res }}$ via the chemical potential $\mu_{i}$ as expressed in Equation (16)
$$
\begin{equation*}
\ln \left(\varphi_{i}\right)=\frac{\mu_{i}^{\text {res }}}{k_{B} T}-\ln \left(1+\left(\frac{\partial\left(\frac{a^{\text {res }}}{k_{B} T}\right)}{\partial \rho}\right)\right) \tag{16}
\end{equation*}
$$

Therefore, with the phase equilibrium calculation, the extraction efficiency of $\mathrm{Li}^{+}$and $\mathrm{Mg}^{2+}$ were determined using:
$$
\begin{equation*}
E_{i}=\frac{\phi \cdot x_{i, \text { org,phase }}}{x_{i, \text { Feed }}} \tag{17}
\end{equation*}
$$
where $\phi$ accounts for the molar ratio between the organic and aqueous phase.

The three binary parameters between $\mathrm{Li}^{+}$and TOP and between TOP and the IL were optimized to enhance the accuracy of the $E_{\mathrm{Li}+}$ data point

\begin{figure}
\includegraphics[alt={},max width=\textwidth]{https://cdn.mathpix.com/cropped/91e8b19b-5f41-4aba-88a0-e1b0958e5483-09.jpg?height=674&width=813&top_left_y=183&top_left_x=159}
\captionsetup{labelformat=empty}
\caption{Fig. 6. Comparison of the calculated (lines) and experimental (points) on the effects of the volume ratio for the organic to aqueous phase (O/A) on $\mathrm{Li}^{+}$and $\mathrm{Mg}^{2+}$ extraction efficiency with [HOEMIM][Tf2 N ] + TOP as extraction system ( $\mathrm{Mg}^{2+} / \mathrm{Li}^{+}=40: 1$; and $\mathrm{pH}=7$ and $\mathrm{C}_{\mathrm{IL}}=0.09 \mathrm{~mol} / \mathrm{L}$ ). ePC-SAFT parameters are taken from Tables 1-3, and the parameters were fitted to $\mathrm{E}_{\mathrm{Li}+}$ at $\mathrm{O} / \mathrm{A}=2: 1$ and 298.15 K , i.e., all other modeling lines are predictions.}
\end{figure}
at $\mathrm{O} / \mathrm{A}=2: 1$. All the binary interaction parameters to $\mathrm{Mg}^{2+}$ were set to zero. Further, dispersion between water and TOP, as well as TOP and IL was not accounted for $\left(k_{i j}=1\right)$. Fig. 6 illustrates the modelled extraction efficiency for $\mathrm{Li}^{+}$and $\mathrm{Mg}^{2+}$ at various ratios between the organic and aqueous phases (see Table S6 in Supplementary material for more details). Not surprisingly, the ePC-SAFT modelling of $E_{\mathrm{Li}+}$ at $\mathrm{O} / \mathrm{A}=2: 1$ is in the quantitative agreement to the experimental data since this was the fitting point. Conversely, the calculations for the other O/A ratios are predictive. The results clearly demonstrate that a ratio $\mathrm{O} / \mathrm{A}=1: 1$ yields a notably lower extraction efficiency, which is captured by ePC-SAFT qualitatively. The ePC-SAFT overestimates the extraction efficiency at this ratio. Anyway, the prediction at ratios $\mathrm{O} / \mathrm{A}>1$ agrees quantitatively with the experimental data. As the $\mathrm{O} / \mathrm{A}$ ratios increase, $E_{\mathrm{Li}+}$ finally arrives at a limiting value of approximately $85 \%$. Additionally, the ePCSAFT correctly predicts that the extraction efficiency of $\mathrm{Mg}^{2+}$ is negligible across all ratios. Remarkably, this behavior was obtained without the need to adjust additionally binary interaction parameters. Therefore, the ePC-SAFT model demonstrates the powerful prediction capacity and is first successfully extended to describe the so-called "organ-ic-inorganic complex strong electrolyte system" in $\mathrm{Li}^{+} / \mathrm{Mg}^{2+}$ extraction separation with ILs.

\subsection*{4.4.2. Analyses of Gibbs free energy, enthalpy and entropy}

The thermodynamic analyses for Gibbs free energy ( $\Delta G$ ), enthalpy $(\Delta H)$ and entropy ( $\Delta S$ ) in the extraction process were investigated. As shown in Fig. 7a, both the extraction efficiency and distribution coefficient of $\mathrm{Li}^{+}$decrease with the increase in temperatures. Based on the Van't Hoff relation, the enthalpy change can be calculated by
$$
\begin{equation*}
\log K=\frac{\Delta H}{2.303 R T}+C \tag{18}
\end{equation*}
$$
where $R$ is the universal gas constant, and $C$ is the integral constant. Fig. 7b shows that the slope of the plot is $\Delta H / 2.303 R$, so the calculated $\Delta H$ is $-10.33 \mathrm{~kJ} / \mathrm{mol}$, indicating the extraction of $\mathrm{Li}^{+}$is an exothermic process, and low temperature is beneficial to this process. $\Delta G$ and the $\Delta S$ can be also obtained by
$$
\begin{equation*}
\Delta G=-2.303 R T \cdot \log K=\Delta H-T \cdot \Delta S \tag{19}
\end{equation*}
$$

Consequently, the values of $\Delta G$ and $\Delta S$ are $-0.90 \mathrm{~kJ} / \mathrm{mol}$ and $-31.65 \mathrm{~J} / \mathrm{K} / \mathrm{mol}$, respectively, indicating that the extraction reaction is orderly spontaneous.

\subsection*{4.5. Molecular insights into extraction processes}

\subsection*{4.5.1. Analysis of binding energy}

It is very significant to identify the extraction mechanism at the molecular level for extraction separation of $\mathrm{Li}^{+} / \mathrm{Mg}^{2+}$. The binding energies of between $\mathrm{Li}^{+} / \mathrm{Mg}^{2+}$ and $\mathrm{H}_{2} \mathrm{O}$, IL, and TOP were calculated by QC calculations, and the results are shown in Fig. 8. In the aqueous solution, $\mathrm{Li}^{+}$can be coordinated with four $\mathrm{H}_{2} \mathrm{O}$ molecules, while $\mathrm{Mg}^{2+}$ can be coordinated with six $\mathrm{H}_{2} \mathrm{O}$ molecules (Marcus, 1988; Olsher et al., 1991). By comparing Fig. 8b and e, it can be found that the binding energy of $\mathrm{Mg}^{2+}-6 \mathrm{H}_{2} \mathrm{O}(-337.41 \mathrm{kcal} / \mathrm{mol})$ is much greater than that of $\mathrm{Li}^{+}-4 \mathrm{H}_{2} \mathrm{O}(-114.39 \mathrm{kcal} / \mathrm{mol})$, indicating that $\mathrm{Mg}^{2+}$ has stronger hydration in the aqueous phase, which makes $\mathrm{Mg}^{2+}$ more difficult to extract from the aqueous phase than $\mathrm{Li}^{+}$. As shown in Fig. 8A-C, the magnitude of interaction energy is in the order of $\mathrm{Li}^{+}-\left[\mathrm{Tf}_{2} \mathrm{~N}\right]^{-}(-133.09 \mathrm{kcal} / \mathrm{mol})>\mathrm{Li}^{+}-4 \mathrm{H}_{2} \mathrm{O}(-114.39 \mathrm{kcal} / \mathrm{mol})>\mathrm{Li}^{+}-\mathrm{TOP}(-67.49 \mathrm{kcal} /$ mol ). It can be concluded that $\mathrm{Li}^{+}$is more inclined to bind to $\left[\mathrm{Tf}_{2} \mathrm{~N}\right]^{-}$in the extraction system, and the neutral ligand TOP coordinate with $\mathrm{Li}^{+}$to ensure the stability of the complex. This is consistent with the experimental results in Fig. 4, indicating that the IL plays a crucial role in the extraction process and only TOP cannot perform the extraction of $\mathrm{Li}^{+}$. In addition, the stronger binding energy of $\mathrm{Li}^{+}-\left[\mathrm{Tf}_{2} \mathrm{~N}\right]^{-}$than that of $\mathrm{Li}^{+}$$4 \mathrm{H}_{2} \mathrm{O}$ indicates that the IL can effectively extract $\mathrm{Li}^{+}$from the aqueous phase. Fig. 8d-f show that there is the order of $\mathrm{Mg}^{2+}-6 \mathrm{H}_{2} \mathrm{O}(-337.41 \mathrm{kcal} / \mathrm{mol})>\mathrm{Mg}^{2+}-\left[\mathrm{Tf}_{2} \mathrm{~N}\right]^{-}(-337.27 \mathrm{kcal} / \mathrm{mol})>\mathrm{Mg}^{2+}$-TOP ( -177.46 $\mathrm{kcal} / \mathrm{mol}$ ) in the magnitude of interaction energies. This indicates that the very strong interaction presented between $\mathrm{Mg}^{2+}$ and $\mathrm{H}_{2} \mathrm{O}$ lead to the

\begin{figure}
\includegraphics[alt={},max width=\textwidth]{https://cdn.mathpix.com/cropped/91e8b19b-5f41-4aba-88a0-e1b0958e5483-09.jpg?height=556&width=1387&top_left_y=1972&top_left_x=340}
\captionsetup{labelformat=empty}
\caption{Fig. 7. Effect of temperature (a) and plot of versus $1000 / T$ (b) for $\mathrm{Li}^{+}$extraction ( $\mathrm{O} / \mathrm{A}=3: 1$ and $\mathrm{pH}=7$ ) at the normal atmosphere.}
\end{figure}

\begin{figure}
\includegraphics[alt={},max width=\textwidth]{https://cdn.mathpix.com/cropped/91e8b19b-5f41-4aba-88a0-e1b0958e5483-10.jpg?height=1585&width=1777&top_left_y=183&top_left_x=144}
\captionsetup{labelformat=empty}
\caption{Fig. 8. Optimized structures and binding energy for clusters by QC calculations at B3LYP/6-31+g(d,p) with DFT-D3(BJ) dispersion correction.}
\end{figure}

IL failed to effectively extract $\mathrm{Mg}^{2+}$. This will cause to low extraction efficiency of $\mathrm{Mg}^{2+}$ and further increase the separation selectivity of $\mathrm{Li}^{+} / \mathrm{Mg}^{2+}$, which is also consistent with the experimental extraction results. In addition, it can be seen from by comparing Fig. $8 \mathrm{~d}, \mathrm{f}$ and h that the binding energy of [HOEMIM] ${ }^{+}-\left[\mathrm{Tf}_{2} \mathrm{~N}\right]^{-}(-81.95 \mathrm{kcal} / \mathrm{mol})$ is the weakest, indicating that metal ions ( $\mathrm{Li}^{+}$and $\mathrm{Mg}^{2+}$ ) bind the anion $\left[\mathrm{Tf}_{2} \mathrm{~N}\right]^{-}$more competitively than the cation [HOEMIM] ${ }^{+}$. Meanwhile, it can be found from Fig. 8g, i and j that the binding energy of $\mathrm{Li}^{+}$[HOEMIM] ${ }^{+}$and $\mathrm{Mg}^{2+}$ - [HOEMIM] ${ }^{+}$are positive, indicating that there is repulsive interaction presented between them. The negative binding energy of $\mathrm{H}_{2} \mathrm{O}-[\mathrm{HOEMIM}]^{+}(-12.07 \mathrm{kcal} / \mathrm{mol})$ can lead to the cation [HOEMIM] ${ }^{+}$entering the aqueous phase, while it can weaken the hydration of $\mathrm{Li}^{+}$, and further promote transfer of $\mathrm{Li}^{+}$from the aqueous to organic phases. It can be seen from Fig. 8k that the binding energy of $\mathrm{Li}^{+}$-TOP- $\left[\mathrm{Tf}_{2} \mathrm{~N}\right]^{-}(-183.60 \mathrm{kcal} / \mathrm{mol})$ is significantly stronger than those of $\mathrm{Li}^{+}$-TOP ( $-67.49 \mathrm{kcal} / \mathrm{mol}$ ) and $\mathrm{Li}^{+}-\left[\mathrm{Tf}_{2} \mathrm{~N}\right]^{-}(-133.09 \mathrm{kcal} /$ mol ), indicating that the TOP and anion $\left[\mathrm{Tf}_{2} \mathrm{~N}\right]^{-}$have the collaborative effect on extraction of $\mathrm{Li}^{+}$.

\subsection*{4.5.2. Electrostatic potential (ESP) analysis}

ESP is an effective method to determine the location and strength of possibly electrostatic interactions between molecules or clusters by analyzing the three-dimensional distribution map of ESP extreme points on the vdW surface of molecules. It can be seen from the Fig. 9a that the global minimum ESP point ( $-48.74 \mathrm{kcal} / \mathrm{mol}$ ) of TOP was distributed near the O of $\mathrm{P}=\mathrm{O}$ group, where will electrostatically interact with other molecules with positive ESP values. This is due to the oxygen atom have strong electronegativity due to the presence of lone-pair electrons, and the $\mathrm{P}=\mathrm{O}$ group greatly enhanced the concentrated distribution of $\pi$-electron leading to the concentration of negative ESP points near the $\mathrm{P}=\mathrm{O}$ group. Meanwhile, the dispersion of positive ESP points distributes around hydrogen atoms in the alkyl chain. Fig. 9b shows that the global minimum point of ESP ( $-114.38 \mathrm{kcal} / \mathrm{mol}$ ) of $\left[\mathrm{Tf}_{2} \mathrm{~N}\right]^{-}$distributes near the oxygen atom, and the maximum point ( $-50.04 \mathrm{kcal} / \mathrm{mol}$ ) appears above fluorine atoms. Notably, the global ESP minimum of $\left[\mathrm{Tf}_{2} \mathrm{~N}\right]^{-} (-114.38 \mathrm{kcal} / \mathrm{mol})$ is significantly more negative than that of TOP $(-48.74 \mathrm{kcal} / \mathrm{mol})$. As a result of the strong electrostatic attraction, $\mathrm{Li}^{+}$ with a positive charge, preferably coordinates with the oxygen atom at the $\mathrm{S}=\mathrm{O}$ group of $\left[\mathrm{Tf}_{2} \mathrm{~N}\right]^{-}$, resulting in the formation of complex $\mathrm{Li}^{+}$-

\begin{figure}
\includegraphics[alt={},max width=\textwidth]{https://cdn.mathpix.com/cropped/91e8b19b-5f41-4aba-88a0-e1b0958e5483-11.jpg?height=1362&width=1814&top_left_y=174&top_left_x=129}
\captionsetup{labelformat=empty}
\caption{Fig. 9. ESP of molecular vdW surfaces of TOP, complex, cations and anions of ILs (The ESP maximum and minimum values are represented by yellow and cyan balls, respectively, units in $\mathrm{kcal} / \mathrm{mol}$ ). (For interpretation of the references to color in this figure legend, the reader is referred to the web version of this article.)}
\end{figure}
$\left[\mathrm{Tf}_{2} \mathrm{~N}\right]^{-}$, as shown in Fig. 9d. The maximum ESP value of $\mathrm{Li}^{+}-\left[\mathrm{Tf}_{2} \mathrm{~N}\right]^{-}$is concentrated near $\mathrm{Li}^{+}$( $+207.18 \mathrm{kcal} / \mathrm{mol}$ ). Due to the strong electronegativity of the oxygen atom on TOP, the O of the $\mathrm{P}=\mathrm{O}$ group in TOP coordinates with $\mathrm{Li}^{+}$, thereby forming $\mathrm{Li}^{+}$-TOP complex as shown in Fig. 9E. In addition, the global maximum ESP point for $\mathrm{Li}^{+}$-TOP appears near $\mathrm{Li}^{+}$( $+274.11 \mathrm{kcal} / \mathrm{mol}$ ), resulting in a strong electrostatic attraction between $\mathrm{Li}^{+}$and other molecule with negative ESP values like $\left[\mathrm{Tf}_{2} \mathrm{~N}\right]^{-}$. Fig. 9c shows that the global ESP maximum point $(+122.39 \mathrm{kcal} / \mathrm{mol}$ ) for [HOEMIM] ${ }^{+}$is distributed around the acid hydrogen atom on the imidazolium ring of cation, which is significantly lower than that of $\mathrm{Li}^{+}$-TOP ( $+274.11 \mathrm{kcal} / \mathrm{mol}$ ). This interaction promotes cation exchange between $\mathrm{Li}^{+}$and [HOEMIM] ${ }^{+}$, causing the transfer of $\mathrm{Li}^{+}$from the aqueous to organic phase by the strong electrostatic interaction between $\mathrm{Li}^{+}$and TOP $/\left[\mathrm{Tf}_{2} \mathrm{~N}\right]^{-}$. These findings are consistent with the experimental, spectroscopic characterization and binding energy analyses.

\subsection*{4.5.3. Independent gradient model (IGM) analysis}

IGM analysis as proposed by Lefebvre et al., can intuitively identify the size, type and location of noncovalent interactions between molecules, the IGM analyses regarding the molecular-level mechanism for extraction of $\mathrm{Li}^{+}$are shown in Fig. 10. As shown in Fig. 10a, a blue-green isosurface and coordination bond ( $\mathrm{Li}^{+} \rightarrow \mathrm{O}=\mathrm{P}$ ) appears between the $\mathrm{Li}^{+}$ and oxygen atom of the $\mathrm{P}=\mathrm{O}$ group, which is due to the coordination interaction between the oxygen atom containing lone pair electrons and
$\mathrm{Li}^{+}$providing empty electron orbitals. Fig. 10b shows that there are four blue green isosurfaces presented between $\mathrm{Li}^{+}$and the oxygen atoms of the $\mathrm{H}_{2} \mathrm{O}$ molecules, indicating the strong $\mathrm{Li}^{+}$hydration formed in the H $\mathrm{O} \cdots \mathrm{Li}^{+}$bond similar to HBs. It can be found in Fig. 10c that there are two blue green isosurfaces located between $\mathrm{Li}^{+}$and the nitrogen with the oxygen of $\left[\mathrm{Tf}_{2} \mathrm{~N}\right]^{-}$, respectively, indicating the electrostatic attraction between $\mathrm{Li}^{+}$and $\left[\mathrm{Tf}_{2} \mathrm{~N}\right]$ in $\mathrm{S}-\mathrm{N} \cdots \mathrm{Li}^{+}$and $\mathrm{S}=\mathrm{O} \cdots \mathrm{Li}^{+}$interactions. Fig. 10d shows that there are two blue-green areas presented between $\mathrm{Mg}^{2+}$ and the oxygen atom on TOP with carbons alkyl chain, respectively. Compared with area of isosurface of $\mathrm{Li}^{+}$-TOP, that of $\mathrm{Mg}^{2+}$-TOP is bluer and larger, representing a stronger interaction formed in the latter. This is consistent with the binding energy and ESP analyses. As shown in Fig. 10e, a huge blue-green isosurface occupies almost all of the space between the $\mathrm{Mg}^{2+}$ and six $\mathrm{H}_{2} \mathrm{O}$ molecules, which represents the very strong $\mathrm{Mg}^{2+}$ hydration formed in the $\mathrm{H}-\mathrm{O} \cdots \mathrm{Mg}$ bonds and the its strength is obviously stronger than that of $\mathrm{Li}^{+}-4 \mathrm{H}_{2} \mathrm{O}$ due to the larger area of isosurfaces in the former when compared with that in the latter (see Fig. 10b). Fig. 10f shows that the larger bule green area is resented in $\mathrm{Mg}^{2+}-\left[\mathrm{Tf}_{2} \mathrm{~N}\right]^{-}$than that in $\mathrm{Li}^{+}-\left[\mathrm{Tf}_{2} \mathrm{~N}\right]^{-}$, indicating that the former has the stronger interaction than the latter. These findings are in the accordance with the binding energy and ESP analyses. Meanwhile, the $\mathrm{Mg}^{2+}-6 \mathrm{H}_{2} \mathrm{O}$ complex has the larger area of isosurfaces than $\mathrm{Mg}^{2+}$-TOP and $\mathrm{Mg}^{2+}-\left[\mathrm{Tf}_{2} \mathrm{~N}\right]^{-}$complexes, indicating that the hydration interaction between $\mathrm{Mg}^{2+}$ and $\mathrm{H}_{2} \mathrm{O}$ molecules is stronger than that between $\mathrm{Mg}^{2+}$ and TOP/IL. This validates that it is difficult to extract $\mathrm{Mg}^{2+}$ from the

\begin{figure}
\includegraphics[alt={},max width=\textwidth]{https://cdn.mathpix.com/cropped/91e8b19b-5f41-4aba-88a0-e1b0958e5483-12.jpg?height=1736&width=1685&top_left_y=183&top_left_x=189}
\captionsetup{labelformat=empty}
\caption{Fig. 10. IGM visualization for different complex. (The isovalue of IGM is 0.01 , and the color scale is shown in ranging from -0.05 to 0.05 au ).}
\end{figure}
aqueous to organic phases by using the extractant of IL + TOP, which further leads to the high separation selectivity for $\mathrm{Li}^{+} / \mathrm{Mg}^{2+}$. For [HOEMIM] $\left[\mathrm{Tf}_{2} \mathrm{~N}\right]$, it is found in Fig. 10g that an isosurface is formed between oxygen atoms of $\left[\mathrm{Tf}_{2} \mathrm{~N}\right]^{-}$and the hydrogen atom on hydroxyl of [HOEMIM] ${ }^{+}$, representing the formation of HB interaction of $\mathrm{O}-\mathrm{H} \cdots \mathrm{O}$. Similarly, hydrogen atoms on alky chain of [HOEMIM] ${ }^{+}$and oxygen atoms of $\left[\mathrm{Tf}_{2} \mathrm{~N}\right]^{-}$can form the electrostatic interaction of $\mathrm{C}-\mathrm{H} \cdots \mathrm{O}$ and vdW dispersion. Fig. 10h shows that the electrostatic interaction of $\mathrm{S}=\mathrm{O} \cdots \mathrm{Li}^{+}$can be formed between $\mathrm{Li}^{+}$and $\left[\mathrm{Tf}_{2} \mathrm{~N}\right]^{-}$except for the coordination interaction of $\mathrm{Li}^{+} \rightarrow \mathrm{O}=\mathrm{P}$ in the complex $\mathrm{Li}^{+}-\mathrm{TOP}-\left[\mathrm{Tf}_{2} \mathrm{~N}\right]^{-}$. Moreover, $\mathrm{C}-\mathrm{H} \cdots \mathrm{O}$ interactions can be also formed between the hydrogen atoms on alkyl chains of TOP and oxygen atoms of $\left[\mathrm{Tf}_{2} \mathrm{~N}\right]^{-}$, and a thin green isosurface exists between the fluorine atom of $\left[\mathrm{Tf}_{2} \mathrm{~N}\right]^{-}$ and hydrogen atoms of TOP, where the vdW dispersion is formed. Therefore, various interactions of $\mathrm{Li}^{+} \rightarrow \mathrm{O}=\mathrm{P}, \mathrm{S}-\mathrm{N} \cdots \mathrm{Li}^{+}, \mathrm{S}=\mathrm{O} \cdots \mathrm{Li}^{+}, \mathrm{C}-$
$\mathrm{H} \cdots \mathrm{O}$, and vdW dispersion are simultaneously existed in the complex of $\mathrm{Li}^{+}$-TOP- $\left[\mathrm{Tf}_{2} \mathrm{~N}\right]^{-}$to form the so-called multi-site interaction, demonstrating that $\mathrm{Li}^{+}$-TOP- $\left[\mathrm{Tf}_{2} \mathrm{~N}\right]^{-}$has stronger interaction than $\mathrm{Li}^{+}$-TOP and $\mathrm{Li}^{+}-\left[\mathrm{Tf}_{2} \mathrm{~N}\right]^{-}$. This also confirm that both TOP and IL contribute to the extraction of $\mathrm{Li}^{+}$, which is consistent with the experimental results of extraction equilibrium, that is, $\mathrm{Li}^{+}$is extracted from aqueous to organic phases by the way of the complex $\mathrm{Li}^{+}-\mathrm{TOP}-\left[\mathrm{Tf}_{2} \mathrm{~N}\right]$. Fig. 10I shows that there is a bule isosurface presented between the hydrogen atom on hydroxyl of [HOEMIM] ${ }^{+}$and the oxygen atom of $\mathrm{H}_{2} \mathrm{O}$, representing the formation of $\mathrm{HB}(\mathrm{H}-\mathrm{O} \cdot \cdot \mathrm{H})$ interaction, which can weaken the hydration of $\mathrm{Li}^{+}$in the aqueous phase and promote its exchange transfer into the organic phase. These above-mentioned findings are consistent with experimental extraction results, spectroscopic characterization, binding energy and ESP analyses. In short, $\mathrm{Li}^{+}$mainly interacts electrostatically with the anion $\left[\mathrm{Tf}_{2} \mathrm{~N}\right]^{-}$of IL in the form of $\mathrm{S}-\mathrm{N} \cdots \mathrm{Li}^{+}$and $\mathrm{S}=\mathrm{O} \cdots \mathrm{Li}^{+}$. The
cation mainly undergoes HB interaction with $\mathrm{H}_{2} \mathrm{O}$ and enters the aqueous phase undertaking ion exchange with $\mathrm{Li}^{+}$, so that $\mathrm{Li}^{+}$in the aqueous phase enters the organic phase of IL + TOP through the $\mathrm{Li}^{+}$-TOP- $\left[\mathrm{Tf}_{2} \mathrm{~N}\right]^{-}$form of the complex. In this case, the coordinate interaction $\left(\mathrm{Li}^{+} \rightarrow \mathrm{O}=\mathrm{P},\right)$ can be formed between $\mathrm{Li}^{+}$and TOP, the electrostatic interactions ( $\mathrm{S}-\mathrm{N} \cdots \mathrm{Li}^{+}$and $\mathrm{S}=\mathrm{O} \cdots \mathrm{Li}^{+}$) can be formed $\mathrm{Li}^{+}-\left[\mathrm{Tf}_{2} \mathrm{~N}\right]^{-}$, as well the $\mathrm{C}-\mathrm{H} \cdots \mathrm{O}$ and vdW interactions can be presented between TOP and $\left[\mathrm{Tf}_{2} \mathrm{~N}\right]^{-}$. Therefore, the IL can act as a co-extractant as well as a stabilizing complex in an electrically neutral form during the $\mathrm{Li}^{+}$ extraction process.

\subsection*{4.6. Washing and stripping experiment results}

During the extraction $\mathrm{Li}^{+}$process, despite a high separation selectivity of $\mathrm{Li}^{+} / \mathrm{Mg}^{2+}$, a small amount of $\mathrm{Mg}^{2+}$ are still extracted into the $\mathrm{Li}^{+}$-loaded organic phase. To remove these excess $\mathrm{Mg}^{2+}$ and obtain highpure $\mathrm{Li}^{+}$organic phase for preparing the product containing $\mathrm{Li}^{+}$, the washing experiment was conducted. According to the cation exchange theory, the binding order of TOP and metal ions is $\mathrm{H}^{+}>\mathrm{Li}^{+}>\mathrm{Na}^{+}> \mathrm{Mg}^{2+}$ (Wang et al., 2018). Therefore, NaCl and LiCl selected as detergents to remove $\mathrm{Mg}^{2+}$ in the organic phase was investigated experimentally. As shown in Fig. 11a (see Tables S7 in Supplementary material for more details), the washing efficiencies of both $\mathrm{Li}^{+}$and $\mathrm{Mg}^{2+}$ rise with the increase in the NaCl concentration. When the NaCl concentration is $0.5 \mathrm{~mol} / \mathrm{L}$, the washing efficiency of $\mathrm{Mg}^{2+}$ reaches $99.29 \%$. In this case, the washing efficiency of $\mathrm{Li}^{+}$can also reach $45.71 \%$, which leads to the severe loss of $\mathrm{Li}^{+}$. To overcome this problem, LiCl was also added to the detergent under the $0.5 \mathrm{~mol} / \mathrm{L} \mathrm{NaCl}$, and the results are shown in the

Fig. 11b (see Tables S8 in Supplementary material for more details). It is found that the washing efficiency of $\mathrm{Li}^{+}$significantly decreases with the increase in the LiCl concentration. Under the LiCl concentration of 0.1 $\mathrm{mol} / \mathrm{L}$, the washing efficiencies of $\mathrm{Li}^{+}$and $\mathrm{Mg}^{2+}$ in the organic phase after extraction are $1.77 \%$ and $99.7 \%$, respectively, demonstrating the effectiveness of LiCl against washing $\mathrm{Li}^{+}$. The $\mathrm{Mg}^{2+}$ concentration in the organic phase after washing is almost zero, and then the $\mathrm{Li}^{+}$and introduced $\mathrm{Na}^{+}$need to be stripped by HCl solution of $0.5 \mathrm{~mol} / \mathrm{L}$ from the organic phase into the aqueous phase to facilitate the production of $\mathrm{Li}_{2} \mathrm{CO}_{3}$ subsequently. The experimental results are shown in Fig. 11c (see Table S9 in Supplementary material more details). When the ratio of $\mathrm{O} / \mathrm{A}$ is $4: 1$, the $\mathrm{Li}^{+}$stripping rate is only $49.2 \%$ due to the small amount of HCl aqueous solution leading to insufficient mass transfer. The stripping efficiency of $\mathrm{Li}^{+}$significantly increase with the increase in O/A ratio ranging from $4: 1$ to $2: 1$. When the $\mathrm{O} / \mathrm{A}$ ratio over $2: 1$, the stripping efficiency of $\mathrm{Li}^{+}$is almost unchanged with increasing the O/A ratio, so the optimal O/A ratio of $1: 1$ is selected with the stripping efficiency of $\mathrm{Li}^{+}$of 99.90 \%.

\subsection*{4.7. Recycling performance of extractant}

The extractant-loaded organic phase is acidic after stripping experiment, so it is not suitable for repeated extraction of $\mathrm{Li}^{+}$. Therefore, ammonia solution of $1 \mathrm{~mol} / \mathrm{L}$ was used to neutralize the $\mathrm{H}^{+}$in the organic phase to regenerate the extractant. It can be seen from the Fig. 12a (see Tables S10 in Supplementary material for more details) that after 6 extractions with the regenerated extractant of IL + TOP, the extraction of $\mathrm{Li}^{+}$has not obvious decrease keep with around $80 \%$.

\begin{figure}
\includegraphics[alt={},max width=\textwidth]{https://cdn.mathpix.com/cropped/91e8b19b-5f41-4aba-88a0-e1b0958e5483-13.jpg?height=1244&width=1405&top_left_y=1250&top_left_x=331}
\captionsetup{labelformat=empty}
\caption{Fig. 11. Effects of $\mathrm{NaCl}(\mathrm{O} / \mathrm{A}=3: 1)(a)$ and $\mathrm{LiCl}\left(\mathrm{O} / \mathrm{A}=3: 1\right.$ and $\left.C_{\mathrm{NaCl}}=0.5 \mathrm{~mol} / \mathrm{L}\right)$ (b) concentrations on washing efficiency of metal ions, and effects of different $\mathrm{O} /$ A ratios on the stripping efficiency of $\mathrm{Li}^{+}$(c) at 298.15 K .}
\end{figure}

\begin{figure}
\includegraphics[alt={},max width=\textwidth]{https://cdn.mathpix.com/cropped/91e8b19b-5f41-4aba-88a0-e1b0958e5483-14.jpg?height=591&width=1407&top_left_y=177&top_left_x=330}
\captionsetup{labelformat=empty}
\caption{Fig. 12. The cycling extraction efficiency of $\mathrm{Li}^{+}$(a) and FTIR spectra of the regenerated extractant of $[\mathrm{HOEMIM}]\left[\mathrm{Tf}{ }_{2} \mathrm{~N}\right]+\mathrm{TOP}$ (b) at 298.15 K .}
\end{figure}

Fig. 12b shows that the FTIR spectra of extractant regenerated once and regenerated 6 times are the same, indicating that the collaborative extractant of IL + TOP has excellent regenerable and cycling extraction performance with the promising potential in the real industrial application.

\section*{5. Conclusions}

In this work, we systematically investigated selective extraction of $\mathrm{Li}^{+}$with IL-based collaborative extractant from brines with high $\mathrm{Mg}^{2+} / \mathrm{Li}^{+}$molar ratio. The [HOEMIM][ $\mathrm{Tf}_{2} \mathrm{~N}$ ] + TOP was selected as the best candidate with the extraction efficiency of $\mathrm{Li}^{+}$up to $83.16 \%$ and separation selectivity of $\mathrm{Li}^{+} / \mathrm{Mg}^{2+}$ up to 742.11 at the single-stage extraction, as well as the extraction efficiency of $\mathrm{Li}^{+}$up to $99.07 \%$ at the fivestage extraction. It was found that $\mathrm{Li}^{+}$can be extracted from the aqueous to organic phases in the form of $\mathrm{Li}^{+}$-TOP- $\left[\mathrm{Tf}_{2} \mathrm{~N}\right]^{-}$complex of $1: 1: 1$ in molar ratio. The thermodynamic model ePC-SAFT was first extended to quantitatively predict the phase equilibria of the so-called "organ-ic-inorganic complex strong electrolyte system" presented in this work, without adjusting additionally binary interaction parameters. The molecular-level extraction mechanisms were investigated by analyzing binding energies, ESP and IGM. The results showed that the $\mathrm{Li}^{+}$and $\mathrm{Mg}^{2+}$ were hydrated by strong interaction of "Metal-O-H", and the $\mathrm{Li}^{+}$ hydration is weaker than $\mathrm{Mg}^{2+}$ hydration. The $\mathrm{Li}^{+}$can be extracted from the aqueous to organic phases by producing the complex of $\mathrm{Li}^{+}$-TOP$\left[\mathrm{Tf}_{2} \mathrm{~N}\right]^{-}$dominated by the multi-site interaction consisting of $\mathrm{Li}^{+} \rightarrow \mathrm{O}=\mathrm{P}, \mathrm{S}-\mathrm{N} \cdots \mathrm{Li}^{+}, \mathrm{S}=\mathrm{O} \cdots \mathrm{Li}^{+} \mathrm{C}-\mathrm{H} \cdots \mathrm{O}$, and vdW dispersion interactions together, which are stronger than $\mathrm{Li}^{+}$hydration energy. In this case, TOP plays an extractant by the coordinate interaction ( $\mathrm{Li}^{+} \rightarrow \mathrm{O}=\mathrm{P}$ ) and the IL plays a co-extractant as well as a stabilizing complex in an electrically neutral form. The collaborative extractant of [HOEMIM][ $\mathrm{Tf}_{2} \mathrm{~N}$ ] + TOP has the good recycling capacity without significant decrease for the extraction performance after six cycles. This work aims to provide theoretical guidance for rational design of a novel IL-based extractant for the high-efficiency extraction of $\mathrm{Li}^{+}$from ${\mathrm{high} \mathrm{Mg}^{2+}} / \mathrm{Li}^{+}$ratio brines.

\section*{CRediT authorship contribution statement}

Gangqiang Yu: Conceptualization, Investigation, Methodology, Validation, Visualization, Writing - original draft. Xinhe Zhang: Data curation, Investigation, Methodology, Validation. Tobias Hubach: Conceptualization, Formal analysis, Investigation, Validation, Writing original draft. Biaohua Chen: Conceptualization, Project administration, Supervision. Christoph Held: Conceptualization, Methodology, Resources, Software, Supervision, Investigation, Validation, Writing review \& editing.

\section*{Declaration of competing interest}

The authors declare that they have no known competing financial interests or personal relationships that could have appeared to influence the work reported in this paper.

\section*{Data availability}

Data will be made available on request.

\section*{Acknowledgments}

This work was financially supported by the National Natural Science Foundation of China under Grants (Nos. 22378006 and 22008003). The authors acknowledge funding from the Alexander von Humboldt Foundation (G. Yu). We acknowledge financial support realized through the DEAL contract with Elsevier.

\section*{Appendix A. Supplementary data}

Supplementary data to this article can be found online at https://doi. org/10.1016/j.ces.2023.119682.

\section*{References}

Andeen, C., Fontanella, J., Schuele, D., 1970. Low-frequency dielectric constant of LiF, $\mathrm{NaF}, \mathrm{NaCl}, \mathrm{NaBr}, \mathrm{KCl}$, and KBr by the method of substitution. Phys. Rev. B 2, 5068. Andersson, M.P., Uvdal, A., 2005. New scale factors for harmonic vibrational frequencies using the B3LYP density functional method with the triple-ζ basis set 6-311+G(d, p). J. Phys. Chem. A 109, 2937-2941.

Ascani, M., Held, C., 2021. Prediction of salting-out in liquid-liquid two-phase systems with ePC-SAFT: Effect of the Born term and of a concentration-dependent dielectric constant. Z. Anorg. Allg. Chem. 647, 1305-1314.
Ascani, M., Sadowski, G., Held, C., 2022. Calculation of multiphase equilibria containing mixed solvents and mixed electrolytes: general formulation and case studies. J. Chem. Eng. Data 67, 1972-1984.

Berthelot, D., 1898. Sur le mélange des gaz. Compt. Rendus. 126, 15.
Bülow, M., Ji, X., Held, C., 2019. Incorporating a concentration-dependent dielectric constant into ePC-SAFT. An application to binary mixtures containing ionic liquids. Fluid Phase Equilib. 492, 26-33.
Bülow, M., Ascani, M., Held, C., 2021a. ePC-SAFT advanced - Part I: Physical meaning of including a concentration-dependent dielectric constant in the born term and in the Debye-Hückel theory. Fluid Phase Equilib. 535, 112967.
Bülow, M., Ascani, M., Held, C., 2021b. ePC-SAFT advanced - Part II: Application to salt solubility in ionic and organic solvents and the impact of ion pairing. Fluid Phase Equilib. 537, 112989.
Cai, C., Hanada, T., Fajar, A.T.N., Goto, M., 2021. An ionic liquid extractant dissolved in an ionic liquid diluent for selective extraction of $\mathrm{Li}(\mathrm{I})$ from salt lakes. Desalination 509.

Cameretti, L.F., Sadowski, G., 2008. Modeling of aqueous amino acid and polypeptide solutions with PC-SAFT. Chem. Eng. Process. Process Intensification 47, 1018-1025.

Cameretti, L.F., Sadowski, G., Mollerup, J.M., 2005. Modeling of aqueous electrolyte solutions with perturbed-chain statistical associated fluid theory. Ind. Eng. Chem. Res. 44, 3355-3362.
Chapeaux, A., Simoni, L.D., Stadtherr, M.A., Brennecke, J.F., 2007. Liquid phase behavior of ionic liquids with water and 1-octanol and modeling of 1-octanol/water partition coefficients. J. Chem. Eng. Data 52, 2462-2467.
Chen, X., Yuan, S., Abdeltawab, A.A., Al-Deyab, S.S., Zhang, J., Yu, L., Yu, G., 2014. Extractive desulfurization and denitrogenation of fuels using functional acidic ionic liquids. Sep. Purif. Technol. 133, 187-193.
Chen, L., Zhou, T., Chen, L., Ye, Y., Qi, Z., Freund, H., Sundmacher, K., 2011. Selective oxidation of cyclohexanol to cyclohexanone in the ionic liquid 1-octyl-3-methylimidazolium chloride. Chem. Comm. 47, 9354-9356.
Chen, L., Zhang, T., Cheng, H., Richards, R.M., Qi, Z., 2020. A microwave assisted ionic liquid route to prepare bivalent Mn 5 O 8 nanoplates for 5-hydroxymethylfurfural oxidation. Nanoscale 12, 17902-17914.
Choubey, P.K., Kim, M.-S., Srivastava, R.R., Lee, J.-C., Lee, J.-Y., 2016. Advance review on the exploitation of the prominent energy-storage element: Lithium. Part I: From mineral and brine resources. Miner. Eng. 89, 119-137.
Cui, L., Jiang, K., Wang, J., Dong, K., Zhang, X., Cheng, F., 2019. Role of ionic liquids in the efficient transfer of lithium by Cyanex 923 in solvent extraction system. AIChE J. 65, e16606.
Cui, L., Li, S., Kang, J., Yin, C., Guo, Y., He, H., Cheng, F., 2021. A novel ion-pair strategy for efficient separation of lithium isotopes using crown ethers. Sep. Purif. Technol. 274, 118989.
Dupont, D., Binnemans, K., 2015. Rare-earth recycling using a functionalized ionic liquid for the selective dissolution and revalorization of Y 2 O 3: Eu 3+ from lamp phosphor waste. Green Chem. 17, 856-868.
Frisch, M., Trucks, G., Schlegel, H., Scuseria, G., Robb, M., Cheeseman, J., Scalmani, G., Barone, V., Mennucci, B., Petersson, G., 2009. Gaussian 09 (Revision D.01).
Gao, D., Yu, X., Guo, Y., Wang, S., Liu, M., Deng, T., Chen, Y., Belzile, N., 2015. Extraction of lithium from salt lake brine with triisobutyl phosphate in ionic liquid and kerosene. Chem. Res. Chin. Univ. 31, 621-626.
Grimme, S., Ehrlich, S., Goerigk, L., 2011. Effect of the damping function in dispersion corrected density functional theory. J. Comput. Chem. 32, 1456-1465.
Gross, J., Sadowski, G., 2001. Perturbed-chain SAFT: An equation of state based on a perturbation theory for chain molecules. Ind. Eng. Chem. Res. 40, 1244-1260.
Gu, D., Sun, W., Han, G., Cui, Q., Wang, H., 2018. Lithium ion sieve synthesized via an improved solid state method and adsorption performance for West Taijinar Salt Lake brine. Chem. Eng. J. 350, 474-483.
Held, C., Cameretti, L.F., Sadowski, G., 2008. Modeling aqueous electrolyte solutions: Part 1. Fully dissociated electrolytes. Fluid Phase Equilib. 270, 87-96.
Held, C., Reschke, T., Mohammad, S., Luza, A., Sadowski, G., 2014. ePC-SAFT revised. Chemical Eng. Res. Des. 92, 2884-2897.
Jiang, C., Cheng, H., Qin, Z., Wang, R., Chen, L., Yang, C., Qi, Z., Liu, X., 2021. COSMORS prediction and experimental verification of 1, 5-pentanediamine extraction from aqueous solution by ionic liquids. Green Energy Environ. 6, 422-431.
Jin, W., Hu, M., Sun, Z., Huang, C.-H., Zhao, H., 2021. Simultaneous and precise recovery of lithium and boron from salt lake brine by capacitive deionization with oxygen vacancy-rich CoP/Co3O4-graphene aerogel. Chem. Eng. J. 420, 127661.
Keller, A., Hlawitschka, M., Bart, H.-J., 2021. Manganese recycling of spent lithium-ion batteries via solvent extraction. Sep. Purif. Technol. 275, 119166.
Kesler, S.E., Gruber, P.W., Medina, P.A., Keoleian, G.A., Everson, M.P., Wallington, T.J., 2012. Global lithium resources: Relative importance of pegmatite, brine and other deposits. Ore. Geol. Rev. 48, 55-69.
Kuang, G., Liu, Y., Li, H., Xing, S., Li, F., Guo, H., 2018. Extraction of lithium from $\beta$-spodumene using sodium sulfate solution. Hydrometallurgy 177, 49-56.
Lei, Z., Dai, C., Zhu, J., Chen, B., 2014. Extractive distillation with ionic liquids: a review. AIChE J. 60, 3312-3329.
Li, Z., Binnemans, K., 2020. Selective removal of magnesium from lithium-rich brine for lithium purification by synergic solvent extraction using $\beta$-diketones and Cyanex 923. AIChE J. 66, e16246.

Li, Z., Binnemans, K., 2021. Opposite selectivities of tri- n -butyl phosphate and Cyanex 923 in solvent extraction of lithium and magnesium. AIChE J. 67, e17219.
Li, Z., Mercken, J., Li, X., Riaño, S., Binnemans, K., 2019b. Efficient and sustainable removal of magnesium from brines for lithium/magnesium separation using binary extractants. ACS Sustainable Chem. Eng. 7, 19225-19234.
Li, Z., Zhang, Z., Onghena, B., Li, X., Binnemans, K., 2021. Ethylammonium nitrate enhances the extraction of transition metal nitrates by tri-n-butyl phosphate (TBP). AIChE J. 67, e17213.
Li, Y.-H., Zhao, Z.-W., Liu, X.-H., Chen, X.-Y., Zhong, M.-L., 2015. Extraction of lithium from salt lake brine by aluminum-based alloys. Trans. Nonferrous. Met. Soc. China 25, 3484-3489.
Li, Y., Zhao, Y., Wang, H., Wang, M., 2019a. The application of nanofiltration membrane for recovering lithium from salt lake brine. Desalination 468.
Liu, X., Zhong, M., Chen, X., Zhao, Z., 2018. Separating lithium and magnesium in brine by aluminum-based materials. Hydrometallurgy 176, 73-77.
Lorentz, H.A., 1881. Ueber die Anwendung des Satzes vom Virial in der kinetischen Theorie der Gase. Ann. Phys. 248, 127-136.
Lu, T., Chen, F., 2012. Multiwfn: a multifunctional wavefunction analyzer. J. Comput. Chem. 33, 580-592.
Marcus, Y., 1988. Ionic radii in aqueous solutions. Chem. Rev. 88, 1475-1498.

Navarro, P., Moreno, D., Larriba, M., García, J., Rodríguez, F., Canales, R.I., Palomar, J., 2023. An overview process analysis of the aromatic-aliphatic separation by liquidliquid extraction with ionic liquids. Sep. Purif. Technol. 316, 123848.
Olsher, U., Izatt, R.M., Bradshaw, J.S., Dalley, N.K., 1991. Coordination chemistry of lithium ion: a crystal and molecular structure review. Chemical Reviews 91, 137-164.
Pramanik, B.K., Nghiem, L.D., Hai, F.I., 2020. Extraction of strategically important elements from brines: Constraints and opportunities. Water Res. 168, 115149.
Qin, L., Zhang, J., Cheng, H., Chen, L., Qi, Z., Yuan, W., 2016. Selection of imidazoliumbased ionic liquids for vitamin E extraction from deodorizer distillate. ACS Sustainable Chem. Eng. 4, 583-590.
Rout, A., Binnemans, K., 2015. Influence of the ionic liquid cation on the solvent extraction of trivalent rare-earth ions by mixtures of Cyanex 923 and ionic liquids. Dalton. Trans. 44, 1379-1387.
Shi, D., Cui, B., Li, L., Peng, X., Zhang, L., Zhang, Y., 2019. Lithium extraction from lowgrade salt lake brine with ultrahigh $\mathrm{Mg} / \mathrm{Li}$ ratio using TBP - kerosene - $\mathrm{FeCl}_{3}$ system. Sep. Purif. Technol. 211, 303-309.
Shi, D., Cui, B., Li, L., Xu, M., Zhang, Y., Peng, X., Zhang, L., Song, F., Ji, L., 2020. Removal of calcium and magnesium from lithium concentrated solution by solvent extraction method using D2EHPA. Desalination 479, 114306.
Shi, C., Duan, D., Jia, Y., Jing, Y., 2014. A highly efficient solvent system containing ionic liquid in tributyl phosphate for lithium ion extraction. J. Mol. Liq. 200, 191-195.
Shi, C., Jing, Y., Jia, Y., 2016. Solvent extraction of lithium ions by tri-n-butyl phosphate using a room temperature ionic liquid. J. Mol. Liq. 215, 640-646.
Shi, C., Jing, Y., Xiao, J., Wang, X., Yao, Y., Jia, Y., 2017. Solvent extraction of lithium from aqueous solution using non-fluorinated functionalized ionic liquids as extraction agents. Sep. Purif. Technol. 172, 473-479.
Shi, D., Zhang, L., Peng, X., Li, L., Song, F., Nie, F., Ji, L., Zhang, Y., 2018. Extraction of lithium from salt lake brine containing boron using multistage centrifuge extractors. Desalination 441, 44-51.
Song, Z., Zhou, T., Zhang, J., Cheng, H., Chen, L., Qi, Z., 2015. Screening of ionic liquids for solvent-sensitive extraction-with deep desulfurization as an example. Chem. Eng. Sci. 129, 69-77.
Song, Z., Zeng, Q., Zhang, J., Cheng, H., Chen, L., Qi, Z., 2016. Solubility of imidazoliumbased ionic liquids in model fuel hydrocarbons: A COSMO-RS and experimental study. J. Mol. Liq. 224, 544-550.
Su, H., Li, Z., Zhang, J., Liu, W., Zhu, Z., Wang, L., Qi, T., 2020. Combining selective extraction and easy stripping of lithium using a ternary synergistic solvent extraction system through regulation of Fe3+ coordination. ACS Sustainable Chem. Eng. 8, 1971-1979.
Su, H., Tan, B., Zhang, J., Liu, W., Wang, L., Wang, Y., Zhu, Z., Qi, T., 2022. Modelling of lithium extraction with TBP/P507-FeCl3 system from salt-lake brine. Sep. Purif. Technol. 282, 120110.
Sun, X., Wang, X., Wan, Y., Guo, Y., Deng, T., Yu, X., 2023. Synthesis of functional ionic liquids with high extraction rate and electroconductivity for lithium-magnesium separation and metallic magnesium production from salt lake brine. Chem. Eng. J. 452, 139610.
Tarascon, J.-M., 2010. Is lithium the new gold? Nature Chem. 2, 510.
Wang, K., Adidharma, H., Radosz, M., Wan, P., Xu, X., Russell, C.K., Tian, H., Fan, M., Yu, J., 2017. Recovery of rare earth elements with ionic liquids. Green Chem. 19, 4469-4493.
Wang, Y., Liu, H., Fan, J., Liu, X., Hu, Y., Hu, Y., Zhou, Z., Ren, Z., 2018. Recovery of lithium ions from salt lake brine with a high magnesium/lithium ratio using heteropolyacid ionic liquid. ACS Sustainable Chem. Eng. 7, 3062-3072.
Wang, Y., Duan, W., Li, R., Zhang, F., Tian, S., Ren, Z., Zhou, Z., 2023. One-step synthesis of heteropolyacid ionic liquid as co-extraction agent for recovery of lithium from aqueous solution with high magnesium/lithium ratio. Desalination 548.
Wang, J., Yang, S., Bai, R., Chen, Y., Zhang, S., 2019. Lithium recovery from the mother liquor obtained in the process of Li2CO3 production. Ind. Eng. Chem. Res. 58, 1363-1372.
Weil, M., Ziemann, S., Schebek, L., 2009. How to assess the availability of resources for new technologies? Case study: lithium a strategic metal for emerging technologies. Metall. Res. Technol. 106, 554-558.
Wolbach, J.P., Sandler, S.I., 1997. Using molecular orbital calculations to describe the phase behavior of hydrogen-bonding fluids. Ind. Eng. Chem. Res. 36, 4041-4051.
Xu, S., Song, J., Bi, Q., Chen, Q., Zhang, W.-M., Qian, Z., Zhang, L., Xu, S., Tang, N., He, T., 2021. Extraction of lithium from Chinese salt-lake brines by membranes: Design and practice. J Membr Sci. 635, 119441.
Yu, X., Cui, J., Liu, C., Yuan, F., Guo, Y., Deng, T., 2021b. Separation of magnesium from high $\mathrm{Mg} / \mathrm{Li}$ ratio brine by extraction with an organic system containing ionic liquid. Chem. Eng. Sci. 229, 116019.
Yu, G., Sui, X., Lei, Z., Dai, C., Chen, B., 2019a. Air-drying with ionic liquids. AIChE J. 65, 479-482.
Yu, G., Mu, M., Li, J., Wu, B., Xu, R., Liu, N., Chen, B., Dai, C., 2020. Imidazolium-based ionic liquids introduced into $\pi$-electron donors: highly efficient toluene capture. ACS Sustainable Chem. Eng. 8, 9058-9069.
Yu, G., Dai, C., Wu, B., Liu, N., Chen, B., Xu, R., 2021a. Chlorine drying with hygroscopic ionic liquids. Green Energy Environ. 6, 350-362.
Yu, X., Fan, X., Guo, Y., Deng, T., 2019b. Recovery of lithium from underground brine by multistage centrifugal extraction using tri-isobutyl phosphate. Sep. Purif. Technol. 211, 790-798.

Yu, G., Gajardo-Parra, N.F., Chen, M., Chen, B., Sadowski, G., Held, C., 2023. Aromatic volatile organic compounds absorption with phenyl-based deep eutectic solvents: A molecular thermodynamics and dynamics study. AIChE J. 69, e18053.
Zhang, L., Li, L., Shi, D., Peng, X., Song, F., Nie, F., Han, W., 2018. Recovery of lithium from alkaline brine by solvent extraction with $\beta$-diketone. Hydrometallurgy 175, 35-42.
Zhang, J., Liu, Y., Liu, W., Wang, L., Chen, J., Zhu, Z., Qi, T., 2021. Mechanism study on the synergistic effect and emulsification formation of phosphine oxide with $\beta$-diketone for lithium extraction from alkaline systems. Sep. Purif. Technol. 279, 119648.

Zhao, A., Liu, J., Ai, X., Yang, H., Cao, Y., 2019. Highly selective and pollution-free electrochemical extraction of lithium by a polyaniline/Lix Mn2 O4 cell. ChemSusChem 12, 1361-1367.

Zhao, Z., Si, X., Liu, X., He, L., Liang, X., 2013. Li extraction from high Mg/Li ratio brine with $\mathrm{LiFePO}_{4} / \mathrm{FePO}_{4}$ as electrode materials. Hydrometallurgy 133, 75-83.
Zheng, H., Dong, T., Sha, Y., Jiang, D., Zhang, H., Zhang, S., 2021. Selective extraction of lithium from spent lithium batteries by functional ionic liquid. ACS Sustainable Chem. Eng. 9, 7022-7029.
Zhou, W., Li, Z., Xu, S., 2021. Extraction of lithium from magnesium-rich solution using tri-n-butyl phosphate and sodium hexafluorophosphate. J. Sustain. Metall. 7, 1368-1378.
Zhou, Z., Qin, W., Liang, S., Tan, Y., Fei, W., 2012. Recovery of lithium using tributyl phosphate in methyl isobutyl ketone and $\mathrm{FeCl3}$. Ind. Eng. Chem. Res. 51, 12926-12932.
Zhou, Z., Liu, H., Fan, J., Liu, X., Hu, Y., Hu, Y., Wang, Y., Ren, Z., 2019. Selective extraction of lithium ion from aqueous solution with sodium phosphomolybdate as a coextraction agent. ACS Sustainable Chem. Eng. 7, 8885-8892.