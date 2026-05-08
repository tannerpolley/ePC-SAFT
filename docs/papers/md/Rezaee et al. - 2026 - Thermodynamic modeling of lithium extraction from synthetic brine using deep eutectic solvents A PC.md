\title{
Thermodynamic modeling of lithium extraction from synthetic brine using deep eutectic solvents: A PC-SAFT approach
}

\author{
Mehran Rezaee ${ }^{\mathrm{a}}$ (D), Farzaneh Feyzi ${ }^{\mathrm{a}, *{ }^{(B)} \text {, Rohaldin Miri }{ }^{\mathrm{b}}}$ \\ ${ }^{\mathrm{a}}$ Thermodynamics Research Laboratory, School of Chemical, Petroleum and Gas Engineering, Iran University of Science and Technology, Tehran 16846-13114, Iran \\ ${ }^{\mathrm{b}}$ School of Chemical, Petroleum and Gas Engineering, Iran University of Science and Technology, Tehran 16765-163, Iran
}

\section*{A R T I C L E I N F O}

\section*{Keywords:}

Lithium extraction
Deep eutectic solvents
Thermodynamic modeling
Density measurement
PC-SAFT equation of state

\begin{abstract}
Renewable energy sources are essential for mitigating global warming, with lithium-ion batteries serving as a key technology for effective energy storage in electric vehicles and grid-scale applications. This research investigates the thermodynamic modeling of lithium extraction from synthetic brine using a hydrophobic deep eutectic solvent (DES). The extraction process involves ion exchange between aqueous ions ( $\mathrm{Li}^{+}, \mathrm{Na}^{+}$) and an organic phase containing DES (tetrabutylammonium chloride (TBAC) and Decanoic acid (DA) with $1: 2$ molar ratio) and Trioctylphosphine oxide (TOPO). Equilibrium reaction constants were calculated using Joback-Reid and Mostafa group contribution methods. Equilibrium calculations were performed using PC-SAFT and ePC-SAFT equations of state (EOS) for the organic phase and aqueous phase, respectively. EOS parameters for the DES and TOPO were derived from experimental density measurement data. Binary interaction parameters of the EOS, in the organic phase, were fitted to experimental equilibrium data. The average absolute relative deviation (AARD) for lithium extraction percentage and selectivity over sodium were $7.89 \%$ and $8.63 \%$, respectively.
\end{abstract}

\section*{1. Introduction}

Lithium extraction and energy storage in batteries are critical components of the global shift toward sustainable energy that significantly impact the economy and development of renewable technologies. The role of lithium in electric vehicles, renewable energy storage, and consumer electronics continues to grow in importance [1]. The efficiency and sustainability of lithium extraction methods are under scrutiny, as traditional methods often lead to environmental degradation and water scarcity, raising concerns about the balance between resource extraction and ecological preservation [2].

Lithium is primarily sourced from hard rock deposits and brine reserves. Each source presents unique environmental challenges and extraction efficiencies. Traditional techniques, such as saline extraction, involve high water consumption and land disruption, while hard-rock mining is energy-intensive and associated with greenhouse gas emissions [2,3].

In contrast, emerging technologies like Direct Lithium Extraction promise faster recovery rates and lower ecological impact by utilizing innovative processes that minimize land degradation and water consumption, addressing some of the critical challenges faced by the
industry today [4].
Organic solvents and Ionic Liquids (ILs) are widely used for lithium extraction. They have exhibited good performance in some cases, but organic solvents have environmental issues, and ILs face severe regeneration problems. Meanwhile, the high cost of ILs is also a challenge. The mentioned problems make these solvents not suitable for large scale and continuous processes and limit their application in an industrial scale [5, $6]$.

Lithium purification using liquid-liquid extraction methods by deep eutectic solvents (DESs) has emerged as a groundbreaking approach [7]. DESs are utilized in many fields of science and technology, such as desulfurization [8], metal extraction [9], $\mathrm{CO}_{2}$ absorption [10], and water treatment [11]. DESs, are composed of a hydrogen bond donor (HBD) and a hydrogen bond acceptor (HBA). Also, these solvents offer significant advantages over traditional solvents due to their low toxicity, biodegradability, and ability to operate under milder conditions, making them particularly suitable for environmentally friendly extraction processes [12].

The significance of utilizing DESs in lithium extraction lies not only in their eco-friendly characteristics, but also in their demonstrated efficiency in solubilizing lithium salts from complex mixtures [13].

\footnotetext{
* Corresponding author.

E-mail address: feyzi@iust.ac.ir (F. Feyzi).
}

Hanada et al. [14] used trioctylphosphine oxide (TOPO) and triphenyl phosphate (TPP) as HBA and thenoyltrifluoroacetone (HTTA) and benzoyltrifluoroacetone (HBTA) as HBD for DES synthesis. They reached up to $99 \%$ lithium extraction at the best operating condition. Zhang et al. [13] used Trialkylphosphine oxide, N,N-bis(2-ethylhexyl) acetamide and tributyl phosphate (TBP) as HBA and HBTA as HBD for lithium extraction and reached up to $87 \%$ extraction at the best. Xue et al. [15] used DL-Menthol as HBA and three different carboxylic acids, including octanoic acid, decanoic acid (DA) and lauric acid with 71 \% lithium extraction.

Yu et al. [16] used methyltrioctylammonium chloride as HBA and DA as HBD for DES synthesis and added Bis(2-EthylHexyl)Phosphate as co-extractant to increase the efficiency of lithium extraction up to $80 \%$. Chen et al. [7] used tetrabutylammonium chloride (TBAC) and oleic acid for DES synthesis and TBP as co-extractant and achieved $77 \%$ lithium extraction. In our previous work [17], TBAC as HBA and DA as HBD, were utilized for making a DES, and TOPO was used as a co-extractant.

Development of thermodynamic modeling tools is essential for optimizing the operational parameters in both hard rock and brine extraction processes. These tools facilitate the management of environmental impacts while ensuring effective resource utilization [18]. Also, thermodynamics plays a critical role in the extraction of lithium from brines, particularly in understanding the energy changes and the feasibility of various extraction processes. The extraction of lithium from brines is an exothermic process, which means that lower temperature is more beneficial for extraction [19,20]. Thermodynamic models can also be used to optimize extraction processes and effective parameters of the extraction performance [21].

The most critical step in conducting thermodynamic modeling of a system is a suitable model that aligns with the nature of the system and its constituents. Many liquid-liquid extraction systems are composed of an aqueous and an organic phase. Thus, an appropriate model with high accuracy for such systems which contain ions in aqueous phase and complex molecules like DES in organic phase is needed.

Activity coefficient models can serve as good option for the modeling of liquid-liquid extraction systems. However, the necessity of experimental binary data is one of the primary challenges of these models that limits their applicability [8,22,23].

Statistical thermodynamics equations of state (EOSs) are suitable for liquid-liquid extraction systems due to their high accuracy for both aqueous and organic phases. The perturbed chain statistical association fluid theory (PC-SAFT) EOS is a strong model for this purpose because of its capability for predicting the properties of long chain molecules like DESs. Also, adding the electrolyte term to PC-SAFT EOS can expand its application to ionic solutions. The statistical thermodynamics models like PC-SAFT possess predictive capabilities and their accuracy can be further enhanced, in specific cases, by fitting binary interaction parameters [24].

A notable limitation of statistical thermodynamics EOSs is their requirement for pure substance parameters. For complex materials derived from combination of other substances, such as DESs, it is necessary to adjust these parameters over density and vapor pressure data at various temperatures [22,24,25].

The thermodynamic behavior of lithium extraction is fundamentally influenced by the interactions between HBA, HBD and other molecules present in the solution. The application of the PC-SAFT EOS model facilitates the accurate prediction of phase equilibria and thermophysical properties in complex systems [26]. In several studies, PC-SAFT EOS was applied to complex systems containing a DES. However, PC-SAFT modeling of lithium extraction systems using DESs has not been reported in the literature. Recently, a research done by Hubach et al. [27] used electrolyte PC-SAFT (ePC-SAFT) EOS for modeling a lithium extraction system by ionic liquids. Since the aqueous phase contains ions, it is necessary to add the electrolyte term. As reported by Hubach et al. [27] ePC-SAFT EOS is suitable for complex systems including ILs with an acceptable average absolute relative deviation (AARD) of $0.1 \%$
for lithium extraction data. Also, Yu et al. [28] applied ePC-SAFT model for calculating activity coefficients in a lithium extraction system using ILs.

In this research, thermodynamic modeling of lithium extraction from brine using a hydrophobic DES was investigated. Besides lithium, the brine contains sodium as a competing cation. Due to the presence of ions in the aqueous phase and long chain organic compounds in the equilibrium system, PC-SAFT and ePC-SAFT EOSs were used for organic and aqueous phases, respectively. PC-SAFT parameters for the organic phase compounds were fitted to density experimental data obtained in this research. Phase equilibrium data were directly used from our previous publication [17].

\section*{2. Thermodynamic modeling}

\subsection*{2.1. Phase equilibrium system}

The liquid-liquid equilibrium system, as described in detail in our previous work [17], includes an aqueous phase and an organic phase. The aqueous phase mainly includes lithium chloride and sodium chloride solutions, while the organic phase contains the DES and TOPO. The DES made of TBAC as the HBA, and DA as the HBD, with a molar ratio of $1: 2$, and TOPO as the co-extractant, is considered for thermodynamic modeling in this work. The details of DES synthesis and structure, confirmed by analyses are presented in our previous work [17]. Due to the high hydrophobicity of the organic phase, its miscibility with the aqueous phase is negligible. The ions exist only in the aqueous phase while organic phase contains four organic components. Fig. 1 depicts the schematic of the phase equilibrium system. According to Fig. 1, due to the high hydrophobicity of the organic phase, mutual solubility between phases is negligible ( $<1 \mathrm{wt} \%$ based on experimental measurements in our previous work [17]). Lithium and sodium do not exist as identical chemical species in both phases. In the aqueous phase, they are present as solvated ions, whereas in the organic phase they are incorporated into neutral, DES-associated complexes formed through ion-exchange reactions. As a result, conventional fugacity equalities cannot be written in the standard form because no chemically identical species coexist in both phases. To address this situation, we adopt a reaction-equilibrium formulation that is thermodynamically equivalent to enforcing chemical potential balance, but expressed through ion-exchange reactions between phase-specific species. This approach does not decouple phase equilibrium from reaction equilibrium; rather, it provides a practical and internally consistent framework for enforcing equilibrium when direct phase-to-phase fugacity equalities are not applicable.

The dominant mechanism for lithium extraction in this system is the ion exchange reaction. Since higher pH values are more favorable for extraction performance, the presence of hydroxide ions in the aqueous phase is necessary. The ion exchange reactions are as follows:

$$
\begin{equation*}
\mathrm{Li}_{(\mathrm{aq})}^{+}+\mathrm{R}-\mathrm{H}_{(\mathrm{org})}+\mathrm{OH}_{(\mathrm{aq})}^{-} \xrightarrow{\mathrm{K} 1 \mathrm{TOPO}} \mathrm{R}-\mathrm{Li}_{(\mathrm{org})}+\mathrm{H}_{2} \mathrm{O}_{(\mathrm{aq})} \tag{A}
\end{equation*}
$$


$$
\begin{equation*}
\mathrm{Na}_{(\mathrm{aq})}^{+}+\mathrm{R}-\mathrm{H}_{(\mathrm{org})}+\mathrm{OH}_{(\mathrm{aq})}^{-} \xrightarrow{K 2 \mathrm{TOPO}} \mathrm{R}-\mathrm{Na}_{(\mathrm{org})}+\mathrm{H}_{2} \mathrm{O}_{(\mathrm{aq})} \tag{B}
\end{equation*}
$$


In these equilibrium reactions, R-H represents the DES. Lithium and sodium cations are exchanged with the hydrogen in the structure of the DES. Then, the released hydrogen ions pair with $\mathrm{OH}^{-}$in the aqueous phase and decrease pH of aqueous phase. The organic phase contains two new species after the reactions ( $\mathrm{R}-\mathrm{Li}_{\text {(org) }}$ and $\mathrm{R}-\mathrm{Na}_{\text {(org) }}$ ). In other words, the extracted lithium and sodium are now attached to DES molecules, forming new compounds that are called RLi and RNa.

Based on the FTIR spectrum reported in our previous work [17], lithium enters the organic phase through ion exchange with hydrogen and subsequently replaces hydrogen within the DES structure. The presence of hydrogen in the aqueous phase leads to a decrease in the aqueous phase pH ; after equilibrium is reached, the pH of the aqueous phase decreases from the initial range of 9-11 to approximately 7.5-8.

\begin{figure}
\includegraphics[alt={},max width=\textwidth]{https://cdn.mathpix.com/cropped/93f7cf5c-3cc6-42a8-8ff1-401396fd69cc-03.jpg?height=574&width=1435&top_left_y=187&top_left_x=316}
\captionsetup{labelformat=empty}
\caption{Fig. 1. Schematic of the equilibrium system.}
\end{figure}

According to the principle of charge balance, the number of moles of lithium and sodium transferred into the organic phase must equal the number of hydrogen ions released into the aqueous phase. TOPO acts as an assisting agent that facilitates the transfer of lithium and sodium into the organic phase by forming complexes with these cations. Once transferred, lithium and sodium bind to the DES structure, replacing the hydrogen that has been released, thereby completing the ion exchange process. Since the DES is hydrophobic and the organic phase contains negligible amount of water, ionic dissociation does not occur; therefore, lithium and sodium do not exist as free ions in the organic phase but are incorporated into the DES structure.

The reaction equilibrium constants ( $K_{1}$ and $K_{2}$ ) of reactions (A) and (B) are calculated thermodynamically using Eq. (1) [29]. In Eq. (1), $\Delta G^{0}$ stands for the standard Gibbs energy change of reaction, $R$ is the universal gas constant and $T$ is the temperature. $\Delta G^{0}$ is defined by Eq. (2) and is calculated by Eq. (3) [29].

$$
\begin{align*}
& K=\exp \left(-\frac{\Delta G^{0}}{R T}\right)  \tag{1}\\
& \Delta G^{0}=\sum_{i=1}^{N} \nu_{i} G_{i}^{0}  \tag{2}\\
& G^{0}=\Delta G_{f i}^{0} \tag{3}
\end{align*}
$$

where, $\Delta G_{f i}^{0}$ is the standard Gibbs energy of formation of component $i$. In this study, $\Delta G_{f}^{0}$ for the components in the aqueous phase are extracted from literature [30,31]. Since the organic phase contains components that their standard Gibbs energies of formation are not reported in literature, their values were calculated using group contribution methods [32,33]. The Gibbs energy of formation for all the species that participate in equilibrium reactions are presented in Table 1. The details

\begin{table}
\captionsetup{labelformat=empty}
\caption{Table 1
Standard Gibbs energy of formation for components participating in equilibrium reaction.}
\begin{tabular}{|l|l|l|}
\hline Component & $\Delta G_{f}^{0}\left(\frac{\mathrm{~kJ}}{\mathrm{~mol}}\right)$ & Reference \\
\hline $\mathrm{H}_{2} \mathrm{O}$ & -237.13 & [30] \\
\hline $\mathrm{Li}^{+}$ & -292.59 & [30,31] \\
\hline $\mathrm{Na}^{+}$ & -261.88 & [30,31] \\
\hline $\mathrm{OH}^{-}$ & -157.28 & [30,31] \\
\hline TBAC & 77.93 & This work \\
\hline DA & -310.59 & This work \\
\hline DES & -543.25 & This work \\
\hline RLi & -707.58 & This work \\
\hline RNa & -657.13 & This work \\
\hline
\end{tabular}
\end{table}
of group contribution calculations are presented in Supporting Information.

The standard Gibbs energy change of reactions along with their corresponding equilibrium constants, are presented in Table 2.

It should be mentioned that since the temperature range of lithium extraction data ( $293.2-313.2^{\circ} \mathrm{C}$ ) is narrow and 298 K lies in this range the reaction equilibrium constants have been assumed constant.

For reactions occurring in non-ideal liquid phase, by neglecting the Poynting factor, Eq. (4) is available [29].

$$
\begin{equation*}
K=\prod_{i}^{n}\left(x_{i} \gamma_{i}\right)^{v_{i}} \tag{4}
\end{equation*}
$$

where $x$ is the mol fraction, $\gamma$ is the activity coefficient, $v$ is the stoichiometric coefficient, and $i$ denotes the species contributing in the reaction. Applying Eq. (4) to the system under study results in Eqs. (5) and 6.

$$
\begin{equation*}
K_{1}=\frac{\left(x_{R L i} \gamma_{R L i}\right)\left(x_{\mathrm{H}_{2} \mathrm{O}} \gamma_{\mathrm{H}_{2} \mathrm{O}}\right)}{\left(x_{\mathrm{Li}^{+}} \gamma_{\mathrm{Li}^{+}}\right)\left(x_{\mathrm{OH}^{-}} \gamma_{\mathrm{OH}^{-}}\right)\left(x_{D E S} \gamma_{D E S}\right)} \tag{5}
\end{equation*}
$$


$$
\begin{equation*}
K_{2}=\frac{\left(x_{R N a} \gamma_{R N a}\right)\left(x_{\mathrm{H}_{2} \mathrm{O}} \gamma_{\mathrm{H}_{2} \mathrm{O}}\right)}{\left(x_{\mathrm{Na}^{+}} \gamma_{\mathrm{Na}^{+}}\right)\left(x_{\mathrm{OH}^{-}} \gamma_{\mathrm{OH}^{-}}\right)\left(x_{D E S} \gamma_{D E S}\right)} \tag{6}
\end{equation*}
$$


The activity coefficients in Eqs. (5) and 6, can be calculated by a suitable model according thermodynamic Eq. (7).

$$
\begin{equation*}
\gamma_{i}=\frac{\widehat{\varphi}_{i}}{\varphi_{0 i}} \tag{7}
\end{equation*}
$$


In Eq. (7), $\widehat{\varphi}_{i}$ stands for the fugacity coefficient of component $i$ in the mixture, $\varphi_{0 i}$ represents the fugacity coefficient of pure component $i$ for non-ionic components, or fugacity coefficient at infinite dilution for ionic components. The fugacity coefficients in Eq. (7) need to be calculated by an EOS.

Since ions are present in the aqueous phase, ePC-SAFT EOS, as developed by Held et al. [34], was applied to this phase. On the other hand, the organic phase involves four components: DES, TOPO, RLi and RNa which have long hydrocarbon chains and association capabilities. Thus, the PC-SAFT model, as developed by Gross and Sadowski [35], was applied to the organic phase. The challenge in the organic phase

\begin{table}
\captionsetup{labelformat=empty}
\caption{Table 2
Standard Gibbs energy change of reactions and equilibrium constants at 298 K .}
\begin{tabular}{lll}
\hline & $\Delta G^{0}\left(\frac{\mathrm{~kJ}}{\mathrm{~mol}}\right)$ & $K$ \\
\hline Reaction A & 48.41 & $3.2983 \times 10^{-9}$ \\
Reaction B & 68.15 & $1.1476 \times 10^{-12}$ \\
\hline
\end{tabular}
\end{table}
modeling is PC-SAFT pure component parameters which are not available for the substances under consideration in literature. Therefore, experimental data are needed for determining these parameters. The perfect data for determining PC-SAFT parameters are density and vapor pressure in a range of temperatures. Since the components in the organic phase have very low vapor pressure, only density data was used to adjust the parameters. The density of pure DES and solution of TOPO in DES in the range of 293.2 to 343.2 K was measured by Anton Paar DMA HPM instrument according to ASTM D4052 with $\pm 0.1 \mathrm{~kg} / \mathrm{m}^{3}$ uncertainty in measurement. Due to the very low amount of RLi and RNa in the organic phase and the difficulty in purifying them, their density could not be measured experimentally. Thus, the PC-SAFT parameters of RLi and RNa were fitted to experimental equilibrium data.

\subsection*{2.2. PC-SAFT EOS model}

The PC-SAFT EOS was first introduced by Gross and Sadowski [35]. This model predicts the deviation of the system from ideal gas behavior by calculating the residual Helmholtz free energy ( $A^{\text {res }}$ ). The model needs 5 parameters for a pure substance including segment number ( $m$ ), segment diameter ( $\sigma$ ), segment energy ( $u / k$ ), effective association volume ( $\kappa^{A_{i} B_{j}}$ ) and association energy ( $\varepsilon^{A_{i} B_{j}} / k$ ) as the input data and calculates $A^{\text {res }}$ by Eq. (8) ( $k$ is the Boltzmann constant).

$$
\begin{equation*}
A^{\text {res }}=A^{H C}+A^{\text {disp }}+A^{\text {Assoc }} \tag{8}
\end{equation*}
$$


In Eq. (8), $A^{H C}$ is the Helmholtz free energy of the hard chain term, $A^{\text {disp }}$ stands for the dispersion term and $A^{\text {Assoc }}$ for the association contribution that accounts for strong intermolecular interactions like hydrogen bonds. Because of the presence of water in the aqueous phase, and hydrogen donating and accepting sites on DES and TOPO molecules in the organic phase, calculation of association term is necessary. To calculate mixture properties, the combining rules of Berthelot-Lorentz, Eqs. (9) and 10 are applied in this study.

$$
\begin{equation*}
\sigma_{i j}=\frac{1}{2}\left(\sigma_{i}+\sigma_{j}\right) \tag{9}
\end{equation*}
$$


$$
\begin{equation*}
u_{i j}=\sqrt{u_{i} u_{j}}\left(1-k_{i j}\right) \tag{10}
\end{equation*}
$$


After determining $A^{\text {res }}$, other properties can be calculated by rigorous
thermodynamic relations.
Since the aqueous phase contains strong ions and there is long-range Columbic interactions, the electrolyte contribution is also needed for a more precise calculation of the activity coefficients. Therefore, in the aqueous phase the Debye-Hükel term was used to describe electrolyte contributions, following the approach of Held et al. [34] and residual Helmholtz energy was calculated using Eq. (11).

$$
\begin{equation*}
A^{\text {res }}=A^{H C}+A^{\text {disp }}+A^{\text {Assoc }}+A^{\text {elec }} \tag{11}
\end{equation*}
$$


The terms in Eq. (11) are defined in details by Held et al. [34] and were used in this work accordingly.

\section*{3. Results and discussion}

\subsection*{3.1. Density measurement}

The experimental density data of the DES (TBAC/DA 1:2) are obtained in this research. These data along with two more literature data for the same DES, are depicted in Fig. 2. The density data measured in this work are also tabulated in Tables S5 and S6 in supporting information along with their uncertainty values.

The slight differences between the three sets of experimental data arise from operating conditions and density measurement methods. To generate a consistent dataset from the three literature sources and remove experimental noise, a Tait correlation [37-40], as introduced by Eq. (12), was fit to all available data points $\left(R^{2}=0.9914\right)$. This smoothed dataset ( 50 points from $293-343 \mathrm{~K}$ ) was then used for PC-SAFT parameter regression.

Eq. (12) presents density ( $D$ ) at ambient pressure as a function of temperature.

$$
\begin{equation*}
D=A_{0}+A_{1} T+A_{2} T^{2}+A_{3} T^{3} \tag{12}
\end{equation*}
$$


The parameters of Eq. (12) are presented in Table 3.
After determining Tait constants, a range of temperatures from 293.2 to 343.2 K was selected. This range of temperature was divided into 1 K steps and the density of DES was calculated at every point by Eq. (12). These density values were used to fit the PC-SAFT parameters of pure DES. These Parameters are presented in Table 4.

The density values of the DES were then calculated by PC-SAFT EOS

\begin{figure}
\includegraphics[alt={},max width=\textwidth]{https://cdn.mathpix.com/cropped/93f7cf5c-3cc6-42a8-8ff1-401396fd69cc-04.jpg?height=840&width=1133&top_left_y=1688&top_left_x=465}
\captionsetup{labelformat=empty}
\caption{Fig. 2. Density of DES as a function of temperature, this work and from [25,36]. The data in this work have $\pm 0.1 \mathrm{~kg} / \mathrm{m}^{3}$ uncertainty in measurement.}
\end{figure}

\begin{table}
\captionsetup{labelformat=empty}
\caption{Table 3
Parameters of Tait correlation of the DES.}
\begin{tabular}{|l|l|l|}
\hline Parameter & Value & Unit \\
\hline $A_{0}$ & 1172.12 & $\frac{\mathrm{kg}}{\mathrm{m}^{3}}$ \\
\hline $A_{1}$ & -1.5749 & $\frac{\mathrm{kg}}{\mathrm{m}^{3} \mathrm{~K}}$ \\
\hline $A_{2}$ & 0.000404 & $\frac{\mathrm{kg}}{\mathrm{m}^{3} \mathrm{~K}^{2}}$ \\
\hline $A_{3}$ & $-5.5079 \times 10^{-6}$ & $\frac{\mathrm{kg}}{\mathrm{m}^{3} \mathrm{~K}^{3}}$ \\
\hline
\end{tabular}
\end{table}

\begin{table}
\captionsetup{labelformat=empty}
\caption{Table 4
PC-SAFT EOS parameters of the DES.}
\begin{tabular}{|l|l|l|}
\hline Parameter & value & unit \\
\hline $m$ & 10.9898 & - \\
\hline $\sigma$ & 3.1570 & Å \\
\hline $u / k$ & 305.09 & K \\
\hline Number of association sites & 2 & - \\
\hline $\kappa^{A_{i} B_{j}}$ & 0.15 & - \\
\hline $\varepsilon^{A_{i} B_{j}} / k$ & 6900 & K \\
\hline
\end{tabular}
\end{table}
and were compared to the experimental data and the Tait correlation in Fig. 3.

Another component in the organic phase is TOPO, which serves as a co-extractant. The PC-SAFT EOS parameters for TOPO are not available in literature because of its solid state at ambient conditions. Therefore, its parameters should be adjusted to the solution density data. The preferred organic phase composition reported in the experimental section of our previous work [17] was 90 \%wt DES and 10 \%wt TOPO. Thus, the density of this solution was measured the same as pure DES. Also, the Tait correlation according to Eq. (12) was used to correlate these experimental data. The obtained Tait correlation parameters ( $R^{2}=$ 0.9994 ) are presented in Table 5.

The PC-SAFT EOS parameters of TOPO were fitted to the 50 extended density values of the solution by means of Table 5 and Eq. (12). The obtained PC-SAFT parameters for TOPO are presented in Table 6.

The experimental density data of the solution, along with the corresponding values of Tait correlation and the prediction of PC-SAFT EOS

\begin{table}
\captionsetup{labelformat=empty}
\caption{Table 5
Parameters of Tait correlation of $10 \%$ TOPO solution in DES.}
\begin{tabular}{|l|l|l|}
\hline Parameter & Value & Unit \\
\hline $A_{0}$ & 890.42 & $\frac{\mathrm{kg}}{\mathrm{m}^{3}}$ \\
\hline $A_{1}$ & 1.3138 & $\frac{\mathrm{kg}}{\mathrm{m}^{3} \mathrm{~K}}$ \\
\hline $A_{2}$ & -0.006050 & $\frac{\mathrm{kg}}{\mathrm{m}^{3} \mathrm{~K}^{2}}$ \\
\hline $A_{3}$ & $6.0606 \times 10^{-6}$ & $\frac{\mathrm{kg}}{\mathrm{m}^{3} \mathrm{~K}^{3}}$ \\
\hline
\end{tabular}
\end{table}

\begin{table}
\captionsetup{labelformat=empty}
\caption{Table 6
PC-SAFT EOS parameters of TOPO.}
\begin{tabular}{|l|l|l|}
\hline Parameter & value & unit \\
\hline $m$ & 13.2090 & - \\
\hline $\sigma$ & 3.7550 & Å \\
\hline $u / k$ & 395.55 & K \\
\hline Number of association sites & 1 & - \\
\hline $\kappa^{A_{i} B_{j}}$ & 0.01 & - \\
\hline $\varepsilon^{A_{i} B_{j}} / k$ & 1520 & K \\
\hline
\end{tabular}
\end{table}
are depicted in Fig. 4.
Table 7 shows the AARD( \%) of Tait correlation and PC-SAFT EOS predicted density from experimental data, where $D$ is the density.

$$
\begin{equation*}
\operatorname{AARD}(\%)=\frac{100}{N} \times \sum_{i=1}^{N} \frac{\left|D_{i}^{c a l c}-D_{i}^{\exp }\right|}{D_{i}^{\exp }} \tag{13}
\end{equation*}
$$


The other two components in the organic phase are RNa and RLi. Their PC-SAFT parameters are fitted to phase equilibrium experimental data. The reason is that their concentration is very low to be purified for density measurements.

\subsection*{3.2. Reaction equilibrium calculation}

Experimental equilibrium compositions ( 36 equilibrium data points [17]) along with reaction equilibrium constants were used as input data for adjusting PC-SAFT parameters of RLi and RNa. Due to the absence of hydrogen in their structures the association term was not considered and

\begin{figure}
\includegraphics[alt={},max width=\textwidth]{https://cdn.mathpix.com/cropped/93f7cf5c-3cc6-42a8-8ff1-401396fd69cc-05.jpg?height=843&width=1133&top_left_y=1683&top_left_x=467}
\captionsetup{labelformat=empty}
\caption{Fig. 3. Comparison of the calculated density and experimental data for DES.}
\end{figure}

\begin{figure}
\includegraphics[alt={},max width=\textwidth]{https://cdn.mathpix.com/cropped/93f7cf5c-3cc6-42a8-8ff1-401396fd69cc-06.jpg?height=837&width=1133&top_left_y=193&top_left_x=467}
\captionsetup{labelformat=empty}
\caption{Fig. 4. Experimental and calculated density of $10 \%$ TOPO solution in DES. The data in this work have $\pm 0.1 \mathrm{~kg} / \mathrm{m}^{3}$ uncertainty in measurement.}
\end{figure}

\begin{table}
\captionsetup{labelformat=empty}
\caption{Table 7
AARD( \%) of Tait correlation and PC-SAFT EOS.}
\begin{tabular}{|l|l|l|}
\hline \multirow[t]{2}{*}{Component} & \multicolumn{2}{|l|}{AARD( \%)} \\
\hline & Tait & PC-SAFT \\
\hline DES & 0.087 & 0.106 \\
\hline TOPO-DES mixture & 0.103 & 0.261 \\
\hline
\end{tabular}
\end{table}
only 3 parameters ( $m, \sigma, u$ ) were fitted. A total of 6 parameters were estimated as the initial guess and the activity coefficients of all components were calculated. The new values of the composition of RLi and RNa were calculated by rearranging Eqs. (5) and 6 into Eqs. (14) and 15, respectively.

$$
\begin{equation*}
x_{R L i}=\frac{K_{1}\left(x_{L L^{+}} \gamma_{L i^{+}}\right)\left(x_{O H^{-}} \gamma_{O H^{-}}\right)\left(x_{D E S} \gamma_{D E S}\right)}{\left(\gamma_{R L i}\right)\left(x_{\mathrm{H}_{2} \mathrm{O}} \gamma_{\mathrm{H}_{2} \mathrm{O}}\right)} \tag{14}
\end{equation*}
$$


$$
\begin{equation*}
x_{R N a}=\frac{K_{2}\left(x_{N a^{+}} \gamma_{\mathrm{Na}^{+}}\right)\left(x_{\mathrm{OH}^{-}} \gamma_{\mathrm{OH}^{-}}\right)\left(x_{D E S} \gamma_{D E S}\right)}{\left(\gamma_{R N a}\right)\left(x_{\mathrm{H}_{2} \mathrm{O}} \gamma_{\mathrm{H}_{2} \mathrm{O}}\right)} \tag{15}
\end{equation*}
$$


The difference between the new calculated compositions and the experimental data was evaluated by an objective function (OF) defined by Eq. (16).

$$
\begin{equation*}
\mathrm{OF}=\sum_{i=1}^{N}\left|x_{i}^{\text {calc }}-x_{i}^{\exp }\right| \tag{16}
\end{equation*}
$$


The Levenberg-Marquardt algorithm was employed for minimizing the objective function. Because of the small values of RLi and RNa mole fractions, the threshold of objective function minimization was considered to be $1 \times 10^{-6}$. The algorithm of fitting PC-SAFT parameters for RLi and RNa is exhibited by means of a flowchart in Fig. 5. The corresponding fitted PC-SAFT parameters are presented in Table 8.

After determining PC-SAFT parameters for all the components, the reaction equilibrium calculation can be performed. The purpose of this calculation is to predict the equilibrium composition of the mixture and evaluate the lithium extraction percentage and selectivity over sodium.

Eq. (17) was used based on the reaction coordinate parameter $(\varepsilon)$ for determining the equilibrium composition [29].

$$
\begin{equation*}
x_{i}=\frac{n_{i 0}+\sum_{j} \nu_{i, j} \varepsilon_{j}}{n_{0}+\sum_{j} \nu_{j} \varepsilon_{j}} \tag{17}
\end{equation*}
$$

where $n_{i 0}$ stands for the initial number of mols of component $i, n_{0}$ is the total number of initial mols of all components in the system, $v$ is the stoichiometric coefficient and $j$ is the reaction counter.

The input parameters for equilibrium calculation are the number of initial mols of all the components, PC-SAFT pure parameters and equilibrium constants of reactions. The reaction coordinates were estimated at the beginning of the calculations. Substituting the required parameters in Eq. (17) will provide mol fractions of all the components. In the next step, the activity coefficients of the participating components in reactions are calculated by the ePC-SAFT model. New compositions of RLi and RNa were calculated by Eqs. (14) and 15. The obtained compositions were compared with the last step compositions and if their differences were negligible, the compositions were considered as equilibrium values. Otherwise, the reaction coordinates were modified using Eq. (17) using the previous values of mol fractions. The flowchart of equilibrium calculations is presented in Fig. 6.

The calculated equilibrium compositions were used to evaluate lithium extraction percentage and selectivity defined by Eqs. (18) and 19 [17].

$$
\begin{equation*}
E_{i}(\%)=\frac{M_{i i n i}-M_{e q}}{M_{i i n i}} \times 100 \tag{18}
\end{equation*}
$$


$$
\begin{equation*}
S=\frac{E_{L i}}{E_{N a}} \tag{19}
\end{equation*}
$$

where $E, S$, and $M$ stand for extraction percentage, selectivity and mass of lithium or sodium, respectively. Subscripts ini and eq represent initial and equilibrium values, respectively.

The AARD of the extraction percentage and selectivity from experimental data were calculated using Eq. (20).

$$
\begin{equation*}
\operatorname{AARD}(\%)=\frac{100}{N} \times \sum_{i=1}^{N} \frac{\left|Y_{i}^{c a l c}-Y_{i}^{\exp }\right|}{Y_{i}^{\exp }} \tag{20}
\end{equation*}
$$


Here $Y$ stands for the extraction percent or selectivity and $N$ is the number of experimental data. Superscripts calc and exp stand for

\begin{figure}
\includegraphics[alt={},max width=\textwidth]{https://cdn.mathpix.com/cropped/93f7cf5c-3cc6-42a8-8ff1-401396fd69cc-07.jpg?height=1715&width=889&top_left_y=191&top_left_x=589}
\captionsetup{labelformat=empty}
\caption{Fig. 5. Flowchart of fitting PC-SAFT parameters for RLi and RNa .}
\end{figure}

\begin{table}
\captionsetup{labelformat=empty}
\caption{Table 8
PC-SAFT parameters for RLi and RNa.}
\begin{tabular}{|l|l|l|l|}
\hline Parameter & RLi & RNa & Unit \\
\hline $m$ & 11.1223 & 11.1504 & - \\
\hline $\sigma$ & 3.8327 & 4.0254 & Å \\
\hline $u / k$ & 392.09 & 427.44 & K \\
\hline $\kappa^{A_{i} B_{j}}$ & 0 & 0 & - \\
\hline $\varepsilon^{A_{i} B_{j}} / k$ & 0 & 0 & K \\
\hline
\end{tabular}
\end{table}
calculated and experimental, respectively. The AARD was $15.11 \%$ for the lithium extraction percentage and $16.73 \%$ for the selectivity. Figs. 7 and 8 show the calculated equilibrium lithium extraction percentage and selectivity versus experimental data.

The obtained AARD is not satisfactory and need to be improved. To
improve the accuracy of the model, binary interaction parameters ( $k_{i j}$ ) were fitted to equilibrium data for the components in the organic phase. Fig. 9 shows the flowchart for fitting $k_{i j}$.

According to Fig. 9, the input parameters are initial mols of all the components, PC-SAFT pure component parameters and equilibrium constants of reactions. Also, an initial guess for binary interaction parameters is required. Activity coefficients of all the components were calculated using the initial $k_{i j}$ values. Then, the mol fractions of RLi and RNa were calculated using Eqs. (14) and 15. If the deviations from experimental mol fractions are negligible, the loop breaks and $k_{i j}$ values are determined. If not, with new $k_{i j}$ values, the calculation continues until the deviation of calculated mol fractions from experimental data is negligible. The determined binary interaction parameters are presented in Table 9.

The binary interaction parameters (Table 9) reveal:

\begin{figure}
\includegraphics[alt={},max width=\textwidth]{https://cdn.mathpix.com/cropped/93f7cf5c-3cc6-42a8-8ff1-401396fd69cc-08.jpg?height=1814&width=886&top_left_y=191&top_left_x=592}
\captionsetup{labelformat=empty}
\caption{Fig. 6. Flowchart of equilibrium calculations.}
\end{figure}
a) Positive $k_{i j}$ between DES-TOPO (0.0623) suggests weaker dispersion interactions than predicted by the Berthelot-Lorentz rule, possibly due to steric hindrance between bulky molecules.
b) Negative $k_{i j}$ values between TOPO- RNa ( -0.0139 ) and RLi-RNa (-0.0127) indicate stronger than expected interactions, likely due to ion-dipole forces between charged complexes and the phosphoryl group of TOPO.
c) The small magnitude of $k_{i j}(<0.07)$ validates the combining rules for this system.

Figs. 10 and 11 depict the deviation of calculated lithium extraction percentage and selectivity by using $k_{i j}$ from experimental data.

The AARD for lithium extraction percentage and selectivity after
determining binary interaction parameters decreased to 7.89 \% and 8.63 $\%$, respectively. The improvement of modeling accuracy is quite clear by comparison of AARD before and after using binary interaction parameters in the organic phase. The achieved AARD values of $7.89 \%$ (extraction) and $8.63 \%$ (selectivity) are within acceptable ranges for complex extraction systems. The main sources of uncertainty are:
1. Limited experimental data for parameter fitting: Only density data were available for DES and TOPO, lacking vapor pressure or heat capacity data that would further constrain parameters.
2. Group contribution approximations: Equilibrium constants calculated from group contribution methods have inherent uncertainties of $\pm 2-10 \%$ [32].

\begin{figure}
\includegraphics[alt={},max width=\textwidth]{https://cdn.mathpix.com/cropped/93f7cf5c-3cc6-42a8-8ff1-401396fd69cc-09.jpg?height=871&width=1147&top_left_y=185&top_left_x=460}
\captionsetup{labelformat=empty}
\caption{Fig. 7. Deviation of calculated extraction percentage from experimental data [17].}
\end{figure}

\begin{figure}
\includegraphics[alt={},max width=\textwidth]{https://cdn.mathpix.com/cropped/93f7cf5c-3cc6-42a8-8ff1-401396fd69cc-09.jpg?height=882&width=1140&top_left_y=1177&top_left_x=465}
\captionsetup{labelformat=empty}
\caption{Fig. 8. Deviation of calculated selectivity from experimental data [17].}
\end{figure}
3. Temperature-independent equilibrium constants: The assumption of constant K over $293-313 \mathrm{~K}$ also extends calculation deviation.
4. Activity coefficient calculation: The PC-SAFT model's accuracy depends on the association scheme chosen. More complex association schemes (e.g., allowing cross-association between DES and metal complexes) might improve predictions.

Despite these limitations, the model successfully captures the extraction behavior and can guide process optimization and scale-up.

\section*{4. Conclusion}

This study presents the first application of PC-SAFT and ePC-SAFT equations of state for modeling lithium extraction from brine using a hydrophobic deep eutectic solvent (TBAC/DA 1:2) with TOPO as coextractant. Key findings include: The predominant mechanism for lithium extraction in this context is ion exchange reactions. Conventional fugacity equalities cannot be written in the standard form because no chemically identical species coexist in both phases. Consequently, the reaction coordinate parameters of the ion exchange reactions are calculated to determine the equilibrium compositions. The equilibrium

\begin{figure}
\includegraphics[alt={},max width=\textwidth]{https://cdn.mathpix.com/cropped/93f7cf5c-3cc6-42a8-8ff1-401396fd69cc-10.jpg?height=1752&width=889&top_left_y=191&top_left_x=589}
\captionsetup{labelformat=empty}
\caption{Fig. 9. Flowchart for fitting $k_{i j}$ in organic phase.}
\end{figure}

\begin{table}
\captionsetup{labelformat=empty}
\caption{Table 9
Binary interaction parameters of the PC-SAFT EOS in organic phase.}
\begin{tabular}{|l|l|l|l|l|}
\hline & DES & TOPO & RLi & RNa \\
\hline DES & - & 0.0623 & 0.0104 & 0.0158 \\
\hline TOPO & 0.0623 & - & 0.0115 & -0.0139 \\
\hline RLi & 0.0104 & 0.0115 & - & -0.0127 \\
\hline RNa & 0.0158 & -0.0139 & -0.0127 & - \\
\hline
\end{tabular}
\end{table}
constants for these reactions are calculated based on the Gibbs energy changes by group contribution methods. Since the PC-SAFT parameters for pure components in the organic phase are not available in literature, they were fitted to experimentally measured density data for TOPO and DES, as well as to reaction equilibrium data of RLi and RNa. The fitted
parameters are subsequently applied to reaction equilibrium calculations. The AARD values obtained for lithium extraction percentage and its selectivity over sodium are found to be $15.11 \%$ and $16.73 \%$, respectively. Furthermore, the binary interaction parameters in the organic phase were adjusted to experimental equilibrium data, enhancing the accuracy of modeling and reducing AARD values for lithium extraction and selectivity to $7.89 \%$ and $8.63 \%$, respectively.

\section*{CRediT authorship contribution statement}

Mehran Rezaee: Writing - original draft, Resources, Methodology, Investigation, Formal analysis, Data curation, Conceptualization. Farzaneh Feyzi: Writing - review \& editing, Supervision, Resources, Project administration, Methodology, Funding acquisition,

\begin{figure}
\includegraphics[alt={},max width=\textwidth]{https://cdn.mathpix.com/cropped/93f7cf5c-3cc6-42a8-8ff1-401396fd69cc-11.jpg?height=864&width=1133&top_left_y=189&top_left_x=467}
\captionsetup{labelformat=empty}
\caption{Fig. 10. Deviation of calculated lithium extraction percentage from experimental data [17] using $k_{i j}$.}
\end{figure}

\begin{figure}
\includegraphics[alt={},max width=\textwidth]{https://cdn.mathpix.com/cropped/93f7cf5c-3cc6-42a8-8ff1-401396fd69cc-11.jpg?height=882&width=1140&top_left_y=1179&top_left_x=465}
\captionsetup{labelformat=empty}
\caption{Fig. 11. Deviation of calculated selectivity from experimental data [17] using $k_{i j}$.}
\end{figure}

Conceptualization. Rohaldin Miri: Writing - review \& editing, Visualization, Conceptualization.

\section*{Declaration of competing interest}

The authors declare that there is no competing financial interest or personal relationship that could influence the research reported in this article.

\section*{Acknowledgements}

This work is based upon research funded by Iran National Science Foundation (INSF) under project No. 4005204 and the authors are humbly grateful for their support.

\section*{Supplementary materials}

Supplementary material associated with this article can be found, in the online version, at doi:10.1016/j.fluid.2026.114737.

\section*{Data availability}

\section*{Data will be made available on request.}

\section*{References}
[1] L. Kölbel, et al., Lithium extraction from geothermal brines in the upper Rhine Graben: a case study of potential and current state of the art, Hydrometallurgy 221 (2023) 106131.
[2] M.L. Vera, W.R. Torres, C.I. Galli, A. Chagnes, V. Flexer, Environmental impact of direct lithium extraction from brines, Nat. Rev. Earth Environ. 4 (3) (2023) 149-165.
[3] E.S. Rentier, C. Hoorn, A.C. Seijmonsbergen, Lithium brine mining affects geodiversity and sustainable development goals, Renew. Sustain. Energy Rev. 202 (2024) 114642.
[4] C.B. Tabelin, et al., Towards a low-carbon society: a review of lithium resource availability, challenges and innovations in mining, extraction and recycling, and future perspectives, Miner. Eng. 163 (2021) 106743.
[5] D. Gao, Y. Guo, X. Yu, S. Wang, T. Deng, Extracting lithium from the high concentration ratio of magnesium and lithium brine using imidazolium-based ionic liquids with varying alkyl chain lengths, J. Chem. Eng. Jpn. 49 (2) (2018) 104-110.
[6] G. Liu, Z. Zhao, A. Ghahreman, Novel approaches for lithium extraction from saltlake brines: a review, Hydrometallurgy 187 (2019) 81-100.
[7] W. Chen, et al., Tailoring hydrophobic deep eutectic solvent for selective lithium recovery from the mother liquor of Li2CO3, Chem. Eng. J. 420 (2021) 127648.
[8] M. Rezaee, F. Feyzi, M.R. Dehghani, Extractive desulfurization of dibenzothiophene from normal octane using deep eutectic solvents as extracting agent, J. Mol. Liq. 333 (2021) 115991.
[9] Z. Yuan, H. Liu, W.F. Yong, Q. She, J. Esteban, Status and advances of deep eutectic solvents for metal separation and recovery, Green Chem. 24 (5) (2022) 1895-1929, https://doi.org/10.1039/D1GC03851F.
[10] A.A. Manafpour, F. Feyzi, M. Rezaee, An environmentally friendly deep eutectic solvent for CO2 capture, Sci. Rep. 14 (1) (2024) 19744.
[11] C.M. Chabib, J.K. Ali, M.A. Jaoude, E. Alhseinat, I.A. Adeyemi, I.M. Al Nashef, Application of deep eutectic solvents in water treatment processes: a review, J. Water Process Eng. 47 (2022) 102663.
[12] Q. Wen, J.-X. Chen, Y.-L. Tang, J. Wang, Z. Yang, Assessing the toxicity and biodegradability of deep eutectic solvents, Chemosphere 132 (2015) 63-69.
[13] L. Zhang, J. Li, L. Ji, L. Li, Separation of lithium from alkaline solutions with hydrophobic deep eutectic solvents based on $\beta$-diketone, J. Mol. Liq. 344 (2021) 117729.
[14] T. Hanada, M. Goto, Synergistic deep eutectic solvents for lithium extraction, ACS Sustain. Chem. Eng. 9 (5) (2021) 2152-2160.
[15] K. Xue, et al., Lithium extraction from aqueous medium using hydrophobic deep eutectic solvents, J. Environ. Chem. Eng. 11 (5) (2023) 110490.
[16] L.-Y. Yu, K.-J. Wu, C.-H. He, Tailoring hydrophobic deep eutectic solvent for selective lithium recovery from dilute aqueous solutions, Sep. Purif. Technol. 281 (2022) 119928.
[17] M. Rezaee, F. Feyzi, S. Javanshir, Application of response surface methodology for selective extraction of lithium using a hydrophobic deep eutectic solvent, Ind. Eng. Chem. Res. 64 (4) (2025) 2294-2308.
[18] S.S. Parker, M.J. Clifford, B.S. Cohen, Potential impacts of proposed lithium extraction on biodiversity and conservation in the contiguous United States, Sci. Total Environ. 911 (2024) 168639.
[19] W.T. Stringfellow, P.F. Dobson, Technology for the recovery of lithium from geothermal brines, Energies 14 (20) (2021) 6805.
[20] L. Zhang, et al., Lithium recovery from effluent of spent lithium battery recycling process using solvent extraction, J. Hazard. Mater. 398 (2020) 122840.
[21] G. Yu, X. Zhang, T. Hubach, B. Chen, C. Held, Highly efficient lithium extraction from magnesium-rich brines with ionic liquid-based collaborative extractants: thermodynamics and molecular insights, Chem. Eng. Sci. 286 (2024) 119682.
[22] E.A. Crespo, et al., Characterization and modeling of the liquid phase of deep eutectic solvents based on fatty acids/alcohols and choline chloride, Ind. Eng. Chem. Res. 56 (42) (2017) 12192-12202.
[23] M.A. Kareem, F.S. Mjalli, M.A. Hashim, M.K.O. Hadj-Kali, F.S. Ghareh Bagh, I. M. Alnashef, Phase equilibria of toluene/heptane with deep eutectic solvents based on ethyltriphenylphosphonium iodide for the potential use in the separation of aromatics from naphtha, J. Chem. Thermodyn. 65 (2013) 138-149.
[24] P.V.A. Pontes, et al., Measurement and PC-SAFT modeling of solid-liquid equilibrium of deep eutectic solvents of quaternary ammonium chlorides and carboxylic acids, Fluid Ph. Equilibria 448 (2017) 69-80.
[25] C.H.J.T. Dietz, et al., PC-SAFT modeling of CO2 solubilities in hydrophobic deep eutectic solvents, Fluid Ph. Equilibria 448 (2017) 94-98.
[26] N. Mgxadeni, et al., Application of PC-SAFT EoS and SCFT for the modeling of thermophysical properties comprising deep eutectic solvent and alcohols, Asia-Pac. J. Chem. Eng. 19 (3) (2024) e3031.
[27] T. Hubach, Z. Er, C. Held, Li+ extraction from aqueous medium using tetracyanoborate ionic liquids-experiments and ePC-SAFT modeling, J. Chem. Eng. Data 69 (6) (2024) 2255-2263.
[28] G. Yu, et al., Lithium extraction from salt-lake brines using dicationic ionic liquids: ePC-SAFT modeling and molecular mechanisms, AIChE J. 71 (9) (2025) e18907.
[29] J.M. Smith, H.C. Van Ness, M.M. Abbott, M.T. Swihart, Introduction to Chemical Engineering Thermodynamics, McGraw-Hill Singapore, 1949.
[30] M. Ejeian, A. Grant, H.K. Shon, A. Razmjou, Is lithium brine water? Desalination 518 (2021) 115169.
[31] E.H. Oelkers, H.C. Helgeson, E.L. Shock, D.A. Sverjensky, J.W. Johnson, V. A. Pokrovskii, Summary of the apparent standard partial molal gibbs free energies of formation of aqueous species, minerals, and gases at pressures 1 to 5000 bars and temperatures 25 to $1000^{\circ} \mathrm{C}$, J. Phys. Chem. Ref. Data 24 (4) (1995) 1401-1560.
[32] K.G. Joback, R.C. Reid, Estimation of pure-component properties from groupcontributions, Chem. Eng. Commun. 57 (1-6) (1987) 233-243.
[33] A.T.M.G. Mostafa, J.M. Eakman, S.L. Yarbro, Prediction of standard heats and Gibbs free energies of formation of solid inorganic salts from group contributions, Ind. Eng. Chem. Res. 34 (12) (1995) 4577-4582.
[34] C. Held, T. Reschke, S. Mohammad, A. Luza, G. Sadowski, ePC-SAFT revised, Chem. Eng. Res. Des. 92 (12) (2014) 2884-2897.
[35] J. Gross, G. Sadowski, Application of the perturbed-chain SAFT equation of state to associating systems, Ind. Eng. Chem. Res. 41 (22) (2002) 5510-5515.
[36] G. Xu, et al., Tuning the composition of deep eutectic solvents consisting of tetrabutylammonium chloride and n-decanoic acid for adjustable separation of ethylene and ethane, Sep. Purif. Technol. 298 (2022) 121680.
[37] J.H. Dymond, R. Malhotra, The Tait equation: 100 years on, Int. J. Thermophys. 9 (6) (1988) 941-951.
[38] G.A. Neece, D.R. Squire, Tait and related empirical equations of state, J. Phys. Chem. 72 (1) (1968) 128-136.
[39] R. Aitbelale, Y. Chhiti, F.E.M.H. Alaoui, A. Sahib Eddine, N. Muñoz Rujas, F. Aguilar, High-pressure soybean oil biodiesel density: experimental measurements, correlation by tait equation, and perturbed chain SAFT (PC-SAFT) modeling, J. Chem. Eng. Data 64 (9) (2019) 3994-4004.
[40] S. Bi, T. Jia, K. Zhao, X. Meng, J. Wu, Liquid density of 2-methoxyethyl acetate, 2ethylhexyl acetate, and diethyl succinate at temperatures from 283.15 K to 363.15 K and pressures up to 100 MPa , J. Chem. Eng. Data 60 (12) (2015) 3532-3538.