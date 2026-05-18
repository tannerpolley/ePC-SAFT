\title{
$\mathrm{Li}^{+}$Extraction from Aqueous Medium Using Tetracyanoborate lonic Liquids-Experiments and ePC-SAFT Modeling
}

\author{
Tobias Hubach, Zeynep Er, and Christoph Held*
}

Cite This: J. Chem. Eng. Data 2024, 69, 2255-2263
Read Online
Downloaded via BRIGHAM YOUNG UNIV on February 25, 2026 at 18:22:19 (UTC). See https://pubs.acs.org/sharingguidelines for options on how to legitimately share published articles.

\begin{abstract}
The global need for rechargeable batteries has increased significantly, leading to a corresponding rise in demand for lithium. This study explores the use of tetracyanoborate-based ionic liquids (IL) and tributyl phosphate (TBP) for the liquid-liquid extraction of lithium ions ( $\mathrm{Li}^{+}$) from aqueous sources. The investigation includes comprehensive experimental analyses and modeling with the electrolyte perturbed-chain statistical associating fluid theory (ePC-SAFT). The Experimental Section investigates the process conditions for the extraction efficiency of the extraction system. It includes investigations into the best TBP/ IL ratio of the TBP/IL extraction agent mixture and the influence of the phase ratio (mass of IL + TBP/mass of the aqueous salt solution) on extraction efficiency. The results proved the synergetic effects of the two extraction agents, IL and TBP, providing maximum efficiency at TBP/IL $=0.85 . \mathrm{Li}^{+}$was almost completely
![](https://cdn.mathpix.com/cropped/a59933d6-bca6-45db-b31f-2c9dd99cc78c-1.jpg?height=477&width=629&top_left_y=821&top_left_x=1303) extracted from an aqueous LiCl solution using optimized conditions. The ePCSAFT modeling approach accounted for the ion-pair-assisted extraction of $\mathrm{Li}^{+}$into the organic phase, enabling the description of the experimental behavior as quantitatively correct. This provided deeper understanding of the thermodynamic behavior within the liquid-liquid extraction system and paves the way for the screening of numerous IL systems for $\mathrm{Li}^{+}$extraction and the prediction of optimal process conditions in the future.
\end{abstract}

\section*{1. INTRODUCTION}

The extraction of lithium, a vital resource in battery production, has gained significant attention in light of the rapid expansion of portable electronic devices, electric vehicles, and renewable energy systems. The growing demand for these technologies requires a reliable and economically viable lithium supply. A striking example of the potential impact of limited lithium availability on the global mobility transition was highlighted by the German Raw Materials Agency. ${ }^{1}$ In 2010, global lithium production reached 28,100 tons, which nearly tripled to 82,500 tons within ten years by 2020 . Forecasts project a substantial increase to $316,000-558,800$ tons by 2030. ${ }^{1}$ On a global scale, lithium resources are estimated to be around 89 million tons, with approximately 22 million tons currently extractable at an economical level, categorized as reserves. ${ }^{2}$ However, the extraction of lithium from underground deposits poses significant challenges and limitations. Notably, the costs associated with lithium extraction from ores are $30-50 \%$ higher than those from brines, and these resources could potentially be depleted in the future. ${ }^{3-5}$

The largest share of lithium resources, totaling 50 million tons, is concentrated in the lithium-rich region known as the Lithium Triangle, spanning across South America, including Chile, Argentina, and Bolivia. These resources predominantly exist in the form of continental brines, primarily restricted to this specific area. In the prevalent evaporation technology employed in extensive ponds, water is evaporated over a period
of 10-24 months, followed by further steps to purify lithium compounds. ${ }^{6}$ However, in other regions of the world, there are lithium brine sources with lower lithium concentrations, such as geothermal brines, where the evaporative technology is not applicable. Consequently, various promising methods for diversified lithium production have been explored in recent times. ${ }^{6}$

Among the various techniques thoroughly investigated, collectively referred to as direct lithium extraction (DLE), are adsorption, ${ }^{7}$ precipitation, ${ }^{8}$ and membrane separation. ${ }^{9,10}$ Further, solvent extraction is promising due to the advantages of continuous operation, slow energy consumption, and noncomplex apparatus. ${ }^{11,12}$ An extraction agent is required, and many different chemical systems have been suggested to efficiently extract $\mathrm{Li}^{+}$from aqueous media, such as the tributyl phosphate (TBP) $-\mathrm{FeCl}_{3}$ extraction agent mixture. ${ }^{13,14}$ While it offers high extraction efficiency, a third phase might occur during the extraction, and a strongly acidic environment is required for re-extraction. ${ }^{11}$ This necessitates specialized

\footnotetext{
Received: February 9, 2024
Revised: April 27, 2024
Accepted: April 29, 2024
Published: May 13, 2024
}
![](https://cdn.mathpix.com/cropped/a59933d6-bca6-45db-b31f-2c9dd99cc78c-1.jpg?height=218&width=151&top_left_y=2342&top_left_x=1781)

\begin{table}
\captionsetup{labelformat=empty}
\caption{Table 1. Characteristics of Chemicals Used in This Study}
\begin{tabular}{|l|l|l|l|l|l|}
\hline component & abbreviation & $M / \mathrm{g} \mathrm{mol}^{-1}$ & CAS & supplier & purity ${ }^{a}$ /wt \% \\
\hline 1-ethyl-3-methylimidazolium tetracyanoborate & [emim][tcb] & 226.05 & 742099-80-5 & Merck KgaA & >99.5 \\
\hline 1-butyl-3-methylimidazolium tetracyanoborate & [bmim][tcb] & 282.16 & 742099-78-1 & Merck KgaA & >99.5 \\
\hline 1-hexyl-3-methylimidazolium tetracyanoborate & [hmim][tcb] & 254.1 & 1240857-50-4 & Merck KGaA & >98 \\
\hline tributyl phosphate & TBP & 266.32 & 126-73-8 & Merck KGaA & 97 \\
\hline lithium chloride monohydrate & $\mathrm{LiCl} \cdot \mathrm{H}_{2} \mathrm{O}$ & 60.39 & 16712-20-2 & Alfa Aesar & 99.95 \\
\hline water & $\mathrm{H}_{2} \mathrm{O}$ & 18.015 & 7732-18-5 & Milli-Q Integral System, Merck KGaA & \\
\hline
\end{tabular}
\end{table}
${ }^{a}$ Purity from the datasheet provided by the supplier.
equipment and large quantities of base to neutralize the organic phase. ${ }^{12}$ In contrast, a promising alternative to conventional organic solvents is the use of ionic liquids (ILs) in liquid-liquid extraction. ${ }^{15}$ Many of such ILs possess desirable properties such as low vapor pressure, high chemical and thermal stability, and a melting temperature below $100^{\circ} \mathrm{C}$. Implementing ILs in the process allows reducing solvent losses to the atmosphere.

The initial exploration of utilizing ILs for $\mathrm{Li}^{+}$extraction was undertaken by Shi et al. ${ }^{16}$ Their research involved employing 1-butyl-3-methyl-imidazolium hexafluorophosphate [bmim]$\left[\mathrm{bf}_{6}\right]$ in conjunction with TBP to extract $\mathrm{Li}^{+}$from an aqueous feed, achieving a $\mathrm{Li}^{+}$extraction efficiency of $90.93 \%$ with a volumetric $R^{\mathrm{V}}(\mathrm{O} / \mathrm{A})$ ratio of $2: 1$ and a volumetric ratio of TBP to IL of 9:1, whereby the extraction efficiency in this work is defined as follows
$$
\begin{equation*}
E_{\mathrm{Li}^{+}}=\frac{n_{\mathrm{Li}^{+}, \text {org }}}{n_{\mathrm{Li}^{+}, \text {init }}} \cdot 100 \% \tag{1}
\end{equation*}
$$

Subsequently, in a follow-up study, the IL 1-butyl-3methylimidazolium bis(trifluoromethylsulfonyl)imide ([bmim] $\left[\mathrm{NTf}_{2}\right]$ ) was also examined, resulting in a slightly higher extraction efficiency of $92.37 \%$ under comparable conditions. ${ }^{17}$

Bai et al. ${ }^{15}$ conducted experiments utilizing a system comprising trialkyl methylammonium di(2-ethylhexyl)orthophosphinate ([N1888] [P507]), TBP, and $\mathrm{FeCl}_{3}$ for $\mathrm{Li}^{+}$extraction, with a specific focus on reducing the required acidity for re-extraction. Notably, a concentration of $[\mathrm{HCl}]= 1.5 \mathrm{~mol} \mathrm{~L}^{-1}$ yielded a re-extraction efficiency of $90 \%$. In a separate study, Li and Binnemans ${ }^{18}$ proposed a novel approach by substituting the commonly used TBP with Cyanex 923, a commercial mixture of trialkyl phosphine oxides predominantly composed of $n$-hexyl and $n$-octyl chains, which demonstrated superior performance compared to TBP in the extraction of various metals. However, the selectivity with Cyanex $923 / \mathrm{FeCl}_{3}$ in their work was disadvantageous, with higher extraction efficiency observed for $\mathrm{Mg}^{2+}$ rather than $\mathrm{Li}^{+}$. Bai et al. ${ }^{19}$ conducted a study utilizing the hydroxylfunctionalized IL 1-hydroxyethyl-3-methyl imidazolium bis(trifluoromethylsulfonyl)imide ([hoemim][ $\left.\mathrm{NTf}_{2}\right]$ ) in combination with TBP for $\mathrm{Li}^{+}$extraction. Their investigation yielded an extraction efficiency of $95.01 \%$, and a $\mathrm{Li}^{+} / \mathrm{Mg}^{2+}$ separation factor of 405.76 was achieved. Zhou et al. ${ }^{11}$ also employed [hoemim] [ $\mathrm{NTf}_{2}$ ] as a coextraction agent in addition to TBP. Under optimal conditions, an extraction efficiency of $94.2 \%$ for $\mathrm{Li}^{+}$was reached.

The existing literature on modeling and predicting liquidliquid extraction is presently constrained, despite the paramount importance of having a rapid screening method for various ILs in process development. Moreover, it proves beneficial to possess the capability to model or even predict the
impact of process conditions, such as the $R^{\mathrm{w}}\left(\frac{\mathrm{O}}{\mathrm{A}}\right)$ ratio, which describes the phase ratio between TBP/IL and the aqueous phase and the composition of the TBP/IL extractant mixture on $\mathrm{Li}^{+}$extraction. In a study by Zhao et al., ${ }^{20}$ the authors employed the conductor-like screening model for real solvents (COSMO-RS) to provide an initial assessment and screening of various ILs for their potential to extract $\mathrm{Li}^{+}$ions. However, it is worth noting that their experimental validation solely involved the use of pure ILs as extractants. Surprisingly, none of the synthesized ILs tested in the experimental part of their work exhibited an extraction efficiency for $\mathrm{Li}^{+}$ions exceeding $34.6 \%$, revealing that ILs alone are not suitable for the extraction of $\mathrm{Li}^{+}$ions. ${ }^{20}$ Recently, Yu et al. ${ }^{21}$ successfully employed the EoS ePC-SAFT for the first time to model $\mathrm{Li}^{+}$ extraction using the promising extraction mixture ([hoemim]$\left.\left[\mathrm{NTf}_{2}\right]\right)+$ trioctyl phosphate. A high $\mathrm{Li}^{+}$-extraction efficiency $(83.16 \%)$ and separation factor between $\mathrm{Li}^{+}$and $\mathrm{Mg}^{2+}$ (742.11), outperforming all existing literature records, were modeled and validated through experimental investigations.

While the existing literature focuses mainly on extraction systems using the anion bis(trifluoromethylsulfonyl)imide $\mathrm{NTf}_{2},{ }^{11,18,22-24}$ this study aims to investigate the extraction of $\mathrm{Li}^{+}$using tetracyanoborate-based (tcb) IL anions, which, to our knowledge, have never been proposed for this application in literature before. Mota-Martinez et al. ${ }^{25}$ stated that [hmim][tcb] has a lower viscosity and lower eco-toxicity compared to fluorinated ILs and can therefore be beneficial in industrial applications. Furthermore, the tcb-based ILs employed in this study demonstrate water stability, in contrast to tetrafluoroborate and hexafluorophosphate-based ILs, which undergo hydrolysis, leading to the formation of hydrofluoric acid upon exposure to water. ${ }^{26,27}$

This study seeks to illustrate the potential of tcb-based ILs for $\mathrm{Li}^{+}$-extraction while bridging the existing gap in the field of modeling and predicting liquid-liquid extraction processes. The proposed study aims to provide insights into the performance and efficiency of tcb-based ILs as coextraction agents for $\mathrm{Li}^{+}$extraction. Moreover, the development of a thermodynamic model for this extraction system enables the prediction and optimization of the process, ultimately contributing to the advancement of novel $\mathrm{Li}^{+}$extraction methods.

\section*{2. EXPERIMENTAL SECTION}
2.1. Chemicals. All the reagents used were of analytical reagent grade. The extractants (shown in Table 1) 1-ethyl-3methylimidazolium tetracyanoborate [emim][tcb] and 1-hexyl-3-methylimidazolium tetracyanoborate [hmim][tcb], 1-butyl-3-methylimidazolium tetracyanoborate [bmim][tcb] and tributyl phoshate (TBP) were purchased from Merck, Germany. Lithium chloride monohydrate $\left(\mathrm{LiCl} \cdot 1 \mathrm{H}_{2} \mathrm{O}\right)$ was
purchased by Thermo Fisher Scientific, Germany. The aqueous solutions were prepared using deionized water (Milli-Q), produced by the Millipore system (Millipore System, France) with an electrical conductivity $<0.55 \mu \mathrm{~S} \mathrm{~cm}^{-1}$ and $\mathrm{pH}=7.34 \pm$ 0.2 .
2.2. Extraction Procedure and Analytics. All the experiments in this work were conducted at $T=(294.16 \pm 1 \mathrm{~K})$ at 1 atm , with at least two experiments per data point. The organic phase contained a certain concentration of TBP and IL. The aqueous electrolyte solution and the organic phase were mixed in Eppendorf tubes ( 2 mL ) using a pulsing vortex mixer (VWR, Germany) with a specified weight-fraction-based $R^{\mathrm{w}}\left(\frac{\mathrm{O}}{\mathrm{A}}\right)$ ratio at $21^{\circ} \mathrm{C}$ for 2 h . After stirring, the mixtures were centrifuged for 3 min at $13,400 \mathrm{rpm}$ for phase separation. Subsequently, organic or aqueous samples were taken from the bottom or top phases by using a syringe. The concentrations of $\mathrm{Li}^{+}$and other cations in the aqueous phase before and after extraction were determined by ion chromatography using a Dionex ICS-2100 after dilution with Millipore water. A precision balance (Mettler Toledo NewClassic MF) was used with an accuracy of $\pm 0.1 \mathrm{mg}$. The cations were analyzed with an analytical column of $2 \times 30 \mathrm{~mm}$ IonPac CS16-Fast-4 $\mu \mathrm{m}$ (Dionex), at $40^{\circ} \mathrm{C}$ with a $10 \mu \mathrm{~L}$ injection volume and a conductivity detector. Millipore water with a 30 mM concentration of methanesulfonic acid served as an eluent. For quantitative analysis, the cation concentrations were calibrated beforehand. The ion concentrations in the organic phase were calculated by mass balance. Besides the extraction efficiency ( $E_{i}$ ) of ions (cf. eq 1), the mass-based distribution coefficient of ions $D_{i}^{\mathrm{w}}$ between the aqueous and organic phases was determined by eq 2 .
$$
\begin{equation*}
D_{i}^{\mathrm{w}}=\frac{w_{i, \mathrm{org}}}{w_{i, \mathrm{aq}}} \tag{2}
\end{equation*}
$$

To determine the solubility of water in TBP, depending on the LiCl concentration, reversed-phase high-performance liquid chromatography (RP-HPLC) was used. A chromatograph of the Agilent Technologies Inc. type 2160 Series was used for this purpose. In this work, a C18 column (Agilent Poroshell 120 EC) was used. The eluent was pure acetonitrile, and the samples were diluted using butyl acetate. The wavelength of the UV/vis detector was set to 265 nm . A flow rate of 1 mL $\min ^{-1}$ through the dense, solid bed of the column was applied.

\section*{3. THERMODYNAMIC MODELING}
3.1. ePC-SAFT. The ePC-SAFT framework in this study was employed to compute the residual Helmholtz energy ( $a^{\text {res }}$ ) of the extraction systems. The framework incorporates various contributions to Helmholtz energy, building upon the classical PC-SAFT model initially proposed by Gross and Sadowski, ${ }^{28}$ which considers dispersive perturbations (represented by the contribution $a^{\text {disp }}$ ) and associating perturbations (represented by the contribution $a^{\text {assoc }}$ ) of a hard-chain reference fluid (described by $a^{\text {hc }}$ ). Moreover, the ePC-SAFT approach incorporates the contributions of electrostatic interactions through the Debye-Hückel theory and Born theory, represented by the Helmholtz energy contributions $a^{\mathrm{DH}}$ and a ${ }^{\text {Born }}$
$$
\begin{equation*}
a^{\text {res }}=a^{\text {hc }}+a^{\text {disp }}+a^{\text {assoc }}+a^{\text {DH }}(\epsilon(x))+a^{\text {Born }}(\epsilon(x)) \tag{3}
\end{equation*}
$$

Equation 3 will be further denoted ePC-SAFT advanced in this work. The dielectric constant $\epsilon_{\mathrm{r}}(x)$ of the medium is utilized in the computation of $a^{\mathrm{DH}}$ and $a^{\text {Born }}$. The required purecomponent values $\epsilon_{\mathrm{r}, j}$ for the solvents and ions used in the calculation are listed in Table 2. All ions were assigned the

\begin{table}
\captionsetup{labelformat=empty}
\caption{Table 2. Dielectric Constants for All Components Applied in This Work}
\begin{tabular}{|l|l|l|}
\hline component & dielectric constant/C, $\mathrm{V} \mathrm{m}^{-1}$ & ref \\
\hline water & -105.2 $\ln (T[\mathrm{~K}])+677.480$ & 29 \\
\hline TBP & 11 & set to a constant number in this work \\
\hline [emim] [tcb] & 11 & $a$ \\
\hline $\mathrm{Li}^{+}$ & 8 & $b$ \\
\hline $\mathrm{Mg}^{2+}$ & 8 & $b$ \\
\hline $\mathrm{Cl}^{-}$ & 8 & $b$ \\
\hline
\end{tabular}
\end{table}
${ }^{a}$ For the IL [emim][tcb], a relative dielectric constant of 11 was chosen according to a previous work. ${ }^{30}$ All salts were modeled with the same dielectric constant that is a mean of available experimental data. ${ }^{31}$
same value of $\epsilon_{\mathrm{r}, \text { ion }}=8$. The dielectric constant of the medium $\epsilon_{\mathrm{r}}(x)$ was then determined from a linear combination of the solvent mass fraction and the ion mole fraction
$$
\begin{equation*}
\epsilon_{\mathrm{r}}(x)=\left(\sum_{j=1}^{N^{\text {solv }}} \epsilon_{\mathrm{r}, j} w_{j}^{\text {solv }}\right) x_{\text {solv }}+\sum_{j=1}^{N^{\text {ion }}} \epsilon_{\mathrm{r}, j} x_{j}^{\text {ion }} \tag{4}
\end{equation*}
$$

In this equation, $N^{\text {solv }}$ expresses the total number of components in the salt-free solvent mixture, and $N^{\text {ion }}$ the total number of ions. $x_{j}$ and $w_{j}$ denote the mole fraction, the mass fraction of component $j$, and $w_{j}^{\text {solv }}$ the mass fraction of solvent $j$ in the salt-free solvent mixture, while $x_{\text {solv }}$ is the sum of the mole fraction of all solvents in the overall mixture.

Associating components such as water are assigned five pure-component parameters, namely the segment number $m_{i}^{\text {seg }}$, segment diameter $\sigma_{i}$, dispersion-energy parameter $\frac{u_{i}}{k_{\mathrm{b}}}$, associa-tion-energy parameter $\epsilon^{A_{i} B_{i}}$, and association-volume parameter $\kappa^{A_{i} B_{i}}$. For the ions considered in this work, two pure-component parameters, $m_{i}^{\text {seg }}$ and $\sigma_{i}$, are utilized. It is important to mention that the temperature-corrected ion diameter ( $d_{i}=0.88 \sigma_{i}$ ) was used in $a^{\mathrm{DH}}$ and $a^{\text {Born }}$ (see Acknowledgment). Table 3 presents the pure-component parameters used in this study. The purecomponent parameters for TBP and [emim][tcb] were obtained by fitting the parameters to the pure-component temperature-dependent TBP and [emim][tcb] experimental densities. ${ }^{32,33}$ The $\mathrm{Li}^{+}$parameters were inherited from Held et al., ${ }^{34}$ and association was only allowed between $\mathrm{Li}^{+}$and TBP according to the reaction scheme in eq 9.

The modeling of mixtures necessitates the utilization of the Berthelot-Lorentz and Wolbach-Sandler combining rules
$$
\begin{align*}
& \sigma_{i j}=\frac{1}{2}\left(\sigma_{i}+\sigma_{j}\right)  \tag{5}\\
& u_{i j}=\sqrt{u_{i} u_{j}\left(1-k_{i j}\right)}  \tag{6}\\
& \epsilon^{A_{i} B_{j}}=\frac{1}{2}\left(\epsilon^{A_{i} B_{i}}+\epsilon^{A_{j} B_{j}}\right)\left(1-k_{i j}^{\mathrm{hb}}\right) \tag{7}
\end{align*}
$$

\begin{table}
\captionsetup{labelformat=empty}
\caption{Table 3. ePC-SAFT Pure-Component Parameters Used in This Work}
\begin{tabular}{|l|l|l|l|l|l|l|l|}
\hline component & $m_{i}^{\text {seg }}$ & $\sigma_{i} / \AA$ & $\frac{u_{i}}{k_{\mathrm{b}}} / \mathrm{K}$ & $\frac{\epsilon^{A_{i} B_{i}}}{k_{\mathrm{b}}} / \mathrm{K}$ & $\kappa^{A B_{i}}$ & association scheme & ref \\
\hline water & 1.2047 & $a$ & 353.95 & 2425.7 & 0.0451 & 2B & 35 \\
\hline TBP & 10.796 & 3.2510 & 217.09 & 5000 & 0.01 & 2B & this work \\
\hline [emim][tcb] & 7.1807 & 3.5680 & 338.10 & 5000 & 0.1 & 2B & this work \\
\hline $\mathrm{Li}^{+}$ & 1 & $2.8449{ }^{b}$ & 360.00 & $0{ }^{c}$ & 1 & 2B & this work \\
\hline $\mathrm{Mg}^{2+}$ & 1 & $3.1327{ }^{b}$ & 1500 & & & & 34 \\
\hline $\mathrm{Cl}^{-}$ & 1 & $2.7560{ }^{b}$ & 170.00 & & & & 34 \\
\hline
\end{tabular}
\end{table}
${ }^{a} \sigma=2.7927+\left(10.11 \mathrm{e}^{-0.01775 T}-1.417 \cdot \mathrm{e}^{-0.01146 T}\right)$ with $T$ in $\mathrm{K} .{ }^{b}$ Consider the information in the Acknowledgment section. ${ }^{c}$ Only association between $\mathrm{Li}^{+}$and TBP allowed, all other cross association forbidden, see Table 4.
$$
\begin{equation*}
\kappa^{A_{i} B_{j}}=\sqrt{\kappa^{A_{i} B_{i}} \kappa^{A_{j} B_{j}}}\left(\frac{\sqrt{\sigma_{i} \sigma_{j}}}{\frac{1}{2}\left(\sigma_{i}+\sigma_{j}\right)}\right)^{3} \tag{8}
\end{equation*}
$$

In eq 6, the binary interaction parameter $k_{i j}$ is introduced, allowing for the modification of the dispersion energy of the pair $i j$. Binary interaction parameters were required for the following pairs: ion/ion, ion/water, ion/TBP, and TBP/water. The parameter for the $\mathrm{IL} /$ water pair was fitted to the phase equilibrium between IL and water, while the remaining parameters were fitted to experimental data points as described below. In addition, the binary interaction parameter $k_{i j}^{\mathrm{hb}}$ for the association parameter of mixtures for the TBP/Li ${ }^{+}$system (cf. eq 7) was adjusted to account for the coordination of $\mathrm{Li}^{+}$to the $\mathrm{P}=\mathrm{O}$ group of TBP according the following mechanism
$$
\begin{equation*}
[\mathrm{Cmim}]_{\text {org }}^{+}+n \cdot \mathrm{TBP}_{\text {org }}+\mathrm{Li}_{\mathrm{aq}}^{+} \rightarrow[\mathrm{Li} \cdot n \mathrm{TBP}]_{\text {org }}^{+}+\mathrm{Cmim}_{\mathrm{aq}}^{+} \tag{9}
\end{equation*}
$$

The binary interaction parameters used in this study are listed in Table 4.

\begin{table}
\captionsetup{labelformat=empty}
\caption{Table 4. Binary Interaction Parameters $\boldsymbol{k}_{\boldsymbol{i j}}$ Used in This Work}
\begin{tabular}{|l|l|l|l|}
\hline parameters & $k_{i j}$ & $k_{i j}^{\mathrm{hb}}$ & refs \\
\hline water-[emim][tcb] & -0.0175 & 0 & this work \\
\hline $\mathrm{Li}^{+}$-TBP & -0.01 & -0.68 & this work \\
\hline $\mathrm{Cl}^{-}$-TBP & -0.01 & n.a. & this work \\
\hline water-TBP & -0.0475 & 0 & this work \\
\hline $\mathrm{Li}^{+}$-water & 0 & $1{ }^{a}$ & this work \\
\hline $\mathrm{Cl}^{-}$-water & 0 & n.a. & this work \\
\hline $\mathrm{Li}^{+}-\mathrm{Cl}^{-}$ & 0.669 & n.a. & 34 \\
\hline $\mathrm{Li}^{+}$-[emim][tcb] & -0.035 & $1{ }^{a}$ & this work \\
\hline $\mathrm{Cl}^{-}-[\mathrm{emim}][\mathrm{tcb}]$ & -0.035 & n.a. & this work \\
\hline $\mathrm{Li}^{+}-[\mathrm{bmim}][\mathrm{tcb}]$ & 0.12 & $1{ }^{a}$ & this work \\
\hline $\mathrm{Cl}^{-}-[\mathrm{bmim}][\mathrm{tcb}]$ & 0.12 & n.a. & this work \\
\hline $\mathrm{Li}^{+}-[\mathrm{hmim}][\mathrm{tcb}]$ & 0.15 & $1{ }^{a}$ & this work \\
\hline $\mathrm{Cl}^{-}-[\mathrm{hmim}][\mathrm{tcb}]$ & 0.15 & n.a. & this work \\
\hline
\end{tabular}
\end{table}
${ }^{a}$ Set to 1.00 in order to forbid association between $\mathrm{Li}^{+} / \mathrm{IL}$ and $\mathrm{Li}^{+}$/ water.

In this study, the isofugacity criterion was utilized to calculate phase equilibria. This criterion asserts that, at thermodynamic equilibrium, the fugacity of each component is identical in all phases.
$$
\begin{equation*}
f_{i}^{1}=f_{i}^{2}=\ldots=f_{i}^{\pi} \tag{10}
\end{equation*}
$$

For the calculation of the fugacities, the so-called " $\varphi-\varphi$ " criterion was employed. The fugacity of each component in
both phases was expressed using the respective fugacity coefficient as shown below
$$
\begin{equation*}
\varphi_{i}^{\mathrm{L} 1} x_{i}^{\mathrm{L} 1}=\varphi_{i}^{\mathrm{L} 2} x_{i}^{\mathrm{L} 2} \tag{11}
\end{equation*}
$$

The fugacity coefficients $\varphi_{i}$ were with ePC-SAFT advanced from the residual Helmholtz energy $a^{\text {res }}$ via the chemical potential $\mu_{i}$, as expressed in eq 12
$$
\begin{equation*}
\ln \left(\varphi_{i}\right)=\frac{\mu_{i}^{\text {res }}}{k_{\mathrm{B}} T}-\ln \left(1+\left(\frac{\partial\left(\frac{a^{\text {res }}}{k_{\mathrm{B}} T}\right)}{\partial \rho}\right)\right) \tag{12}
\end{equation*}
$$

Consequently, through the phase equilibrium calculations, we determined the extraction efficiency of $\mathrm{Li}^{+}$using the following equation
$$
\begin{equation*}
E_{i}=\frac{\phi \cdot x_{i, \text { org phase }}}{x_{i, \text { feed }}} \tag{13}
\end{equation*}
$$
where $\phi$ accounts for the mole-based ratio between the organic and aqueous phase.

The evaluation of the discrepancy between modeled and experimental data was evaluated by the average relative deviation (ARD) and the average absolute deviation (AAD), as indicated by equation
$$
\begin{align*}
& \mathrm{ARD}=\frac{100}{\mathrm{NP}} \sum_{i}^{\mathrm{NP}}\left|1-\frac{z_{i}^{\mathrm{ePC}-\mathrm{SAFT}}}{z_{i}^{\exp , \text { mean }}}\right|  \tag{14}\\
& \mathrm{AAD}=\frac{1}{\mathrm{NP}} \sum_{i}^{\mathrm{NP}}\left|z_{i}^{\mathrm{ePC}-\mathrm{SAFT}}-z_{i}^{\exp , \text { mean }}\right| \tag{15}
\end{align*}
$$
where NP is the number of data points, $z_{i}^{\exp }$ denotes experimental data ( $E_{\mathrm{Li}^{+}}$or solubility) at a given temperature and pressure, and $z_{i}^{\text {ePC-SAFT }}$ denotes ePC-SAFT advanced modeling result ( $E_{\mathrm{Li}^{+}}$or solubility).

\section*{4. RESULTS AND DISCUSSION}
4.1. Modeling the Underlying Equilibria. The efficiency of $\mathrm{Li}^{+}$extraction depends on both the phase ratio of the organic to aqueous phase after extraction and the concentration of $\mathrm{Li}^{+}$in the organic phase. To predict the behavior of the five-component system comprising $\mathrm{Li}^{+}, \mathrm{Cl}^{-}, \mathrm{TBP}$, water, and [emim][tcb] using ePC-SAFT, the parameters for simplified four-component systems needed to be adjusted first, as detailed in this section.
4.1.1. Ionic Liquid + LiCl + Water. Quantitatively modeling the extraction efficiency of various ion IL necessitated the determination of binary interaction parameters $\left(k_{i j}\right)$ between
the IL and the salt ions. These $k_{i j}$ values were determined by fitting the experimental data of this study and are listed in Table 4. Upon analyzing Table 4, it becomes evident that small $k_{i j}$ values were required between IL and ion. In order to reduce the number of fitting parameters, the $k_{i j}$ between the IL and Li ${ }^{+}$ were set to the same values for the pair IL and $\mathrm{Cl}^{-}$, i.e., $k_{i j}$ (IL/ $\left.\mathrm{Li}^{+}\right)=k_{i j}\left(\mathrm{IL} / \mathrm{Cl}^{-}\right)$. This strategy was sufficient to accurately model the IL solubility in water. The specific equilibrium data used to fit these parameters can be found in Tables S1, S3, and S4 in the Supporting Information, along with the corresponding ARD and AAD for each binary system. In conclusion, applying these $k_{i j}$ values enables precise modeling of IL solubility in water.

Figure 1 compares the experimental solubility of three IL cations in water and in water +LiCl with the ePC-SAFT

\begin{figure}
\includegraphics[alt={},max width=\textwidth]{https://cdn.mathpix.com/cropped/a59933d6-bca6-45db-b31f-2c9dd99cc78c-5.jpg?height=638&width=844&top_left_y=816&top_left_x=173}
\captionsetup{labelformat=empty}
\caption{Figure 1. Solubility of IL cations [ emim$]^{+}$(red, squares, dashed), [bmim] ${ }^{+}$(blue, diamonds, solid), and $[\mathrm{hmim}]^{+}$(green, circles, dotted) in water at 1 bar and 294.15 K ; experimental data (symbols), ePC-SAFT advanced (lines) using the parameters from Tables 2-4. Exp. data are given in Tables S1, S3, and S4. Standard uncertainty is $u(T)=0.5 \mathrm{~K} ; u(m)=0.0001 \mathrm{~g}$. The combined expanded uncertainty is $U=0.0363$ with a coverage factor of $k=2$ ( $95 \%$ confidence interval).}
\end{figure}
advanced results. The order of increasing solubility of the IL under the conditions ( $T=294.15 \mathrm{~K}, p=1 \mathrm{~atm}$ ) shown in Figure 1 is as follows: [bmim][tcb] < [hmim][tcb] < [emim][tcb]. Obviously, there is no linear relationship between the chain lengths of the IL cation and the IL solubility. In comparison with the literature under similar conditions (temperature and pressure), the measured values are consistent (cf. Tables S2 and S5 in Supporting Information). Interestingly, LiCl causes a solubility increase for $[\mathrm{emim}]^{+}$, while the solubility of other IL cations is not much influenced by LiCl addition. Note that the data shown do not reflect the complete LLE but are limited to the influence of the LiCl concentration in the aqueous feed on the cation solubility in the aqueous phase. The ePC-SAFT modeling results are in quantitative agreement with the experimental data.
4.1.2. $\mathrm{TBP}+\mathrm{LiCl}+$ Water. Next, the liquid-liquid equilibrium between water and TBP and the influence of the electrolyte LiCl on the mole fraction solubility of TBP in the organic phase and the fraction of $\mathrm{Li}^{+}$in the aqueous phase after extraction were studied. In the absence of LiCl , where only
water and TBP are present, the equilibrium was accurately modeled by applying the parameter $k_{i_{\text {TBP }} \text { /water }}$ to -0.0475 , as demonstrated in Figure $2\left(w_{\text {LiCl }}=0\right)$. Various values have been

\begin{figure}
\includegraphics[alt={},max width=\textwidth]{https://cdn.mathpix.com/cropped/a59933d6-bca6-45db-b31f-2c9dd99cc78c-5.jpg?height=631&width=839&top_left_y=369&top_left_x=1095}
\captionsetup{labelformat=empty}
\caption{Figure 2. Solubility of TBP in the organic phase at 1 bar and 294.15 K. Line: ePC-SAFT advanced using the parameters from Tables 2-4. Symbols: Exp. data according to Table S6. Standard uncertainty is $u(T)=0.5 \mathrm{~K} ; u(m)=0.0001 \mathrm{~g}$. The combined expanded uncertainty is $U=0.0363$, with a coverage factor of $k=2$ ( $95 \%$ confidence interval).}
\end{figure}
reported for the solubility of TBP in the organic phase in the binary water/TBP mixture, but the value measured in this work is within the range of published values (see Table S7 in Supporting Information). The concentration of LiCl in the aqueous feed influences the solubility of TBP in the organic phase. Adjusting the interaction parameter $k_{i j \mathrm{TBP} / \mathrm{Li}}{ }^{+}$allowed for qualitative modeling of the experimental phase behavior, but the experimentally observed solubility minimum could not be modeled quantitatively. Extrema in salting-in/salting-out phenomena are known from the literature, e.g., for ion effects on organic molecules in water.
4.2. IL-Cation Effect on Extraction Performance. Using the same IL anion $[\mathrm{tcb}]^{-}$, we studied the effect of the IL cation on the extraction efficiency of $\mathrm{Li}^{+} E_{\mathrm{Li}^{+}}$using TPB as a neutral extractant. The findings are presented in Figure 3, indicating the $E_{\mathrm{Li}}{ }^{+}$in various imidazolium-based tetracyanoborate ILs as [hmim][tcb] $(48.06 \%)<[$ bmim $][$ tcb $](65.49 \%)<[$ emim[tcb] $(88.63 \%)$, all at $R^{\mathrm{w}}\left(\frac{\mathrm{O}}{\mathrm{A}}\right)=1$. This comparison clearly demonstrates the preference for ILs with shorter IL-cation chain lengths for extracting $\mathrm{Li}^{+}$. These results align with previous findings in the literature, such as those related to $\left[\mathrm{NTf}_{2}\right]$-based ILs. ${ }^{36}$ The longer the IL-cation alkyl chain, the higher the hydrophobicity, which reduces the intermolecular interactions between IL and $\mathrm{Li}^{+}$ions. Moreover, longer ILcation alkyl chains cause increased viscosity, which hampers the mass transport of $\mathrm{Li}^{+}$from the aqueous to the organic phase. Consequently, [emim][tcb] is selected as the IL candidate for subsequent studies.
4.3. Optimizing Operating Conditions for $\mathrm{Li}^{+}$Extraction Processes. The extraction efficiency was determined for different concentrations of $\mathrm{Li}^{+}$in the aqueous feed. The extracant mixture TBP and [emim][tcb] was used for the extraction experiments for $\mathrm{Li}^{+}$concentrations ranging from 0.2 ppm (seawater $\mathrm{Li}^{+}$concentration level) to 2923 ppm (salt lake

\begin{figure}
\includegraphics[alt={},max width=\textwidth]{https://cdn.mathpix.com/cropped/a59933d6-bca6-45db-b31f-2c9dd99cc78c-6.jpg?height=590&width=746&top_left_y=204&top_left_x=221}
\captionsetup{labelformat=empty}
\caption{Figure 3. Effects of different ILs on experimental extraction efficiency $E_{\mathrm{Li}}{ }^{+}$with TBP and IL as the extractant $\left(w_{\mathrm{Li}^{+}}=2000 \mathrm{ppm}\right), R^{\mathrm{w}}\left(\frac{\mathrm{O}}{\mathrm{A}}\right)=1$ , $T=294.15 \mathrm{~K}, p=1.013$ bar, mixing time 2 h at 2500 rpm . Exp. data are given in Table S8. Standard uncertainty is $u(T)=0.5 \mathrm{~K} ; u(m)=$ 0.0001 g . The combined expanded uncertainty is $U=0.001$ with a coverage factor of $k=2$ ( $95 \%$ confidence interval).}
\end{figure}
$\mathrm{Li}^{+}$concentration). Figure 4 shows that most of the $\mathrm{Li}^{+}$can be effectively extracted, regardless of the initial $\mathrm{Li}^{+}$concentration.

\begin{figure}
\includegraphics[alt={},max width=\textwidth]{https://cdn.mathpix.com/cropped/a59933d6-bca6-45db-b31f-2c9dd99cc78c-6.jpg?height=655&width=844&top_left_y=1203&top_left_x=173}
\captionsetup{labelformat=empty}
\caption{Figure 4. Effect of the $\mathrm{Li}^{+}$content in aqueous feed on extraction efficiency $E_{\mathrm{Li}^{+}}$with $[\mathrm{emim}][\mathrm{tcb}]$ and TBP $\left(w_{\mathrm{TBP}, \text { organic }}=0.80\right)$, $R^{\mathrm{w}}\left(\frac{\mathrm{O}}{\mathrm{A}}\right)=1, T=294.15 \mathrm{~K}, p=1.013$ bar, mixing time 2 h at 2500 rpm. Exp. data are given in Table S9. Standard uncertainty is $u(T)= 0.5 \mathrm{~K} ; u(m)=0.0001 \mathrm{~g}$. The combined expanded uncertainty is $U=$ 0.0067 with a coverage factor of $k=2$ ( $95 \%$ confidence interval).}
\end{figure}

However, as the $\mathrm{Li}^{+}$concentration increases, $E_{\mathrm{Li}^{+}}$tends to decrease. At a $\mathrm{Li}^{+}$concentration of 2923 ppm , only $70 \%$ of the initial $\mathrm{Li}^{+}$content can be extracted. Nevertheless, the combination of [emim][tcb] and TBP can effectively extract $\mathrm{Li}^{+}$over a wide concentration range.

The operating conditions, including the IL/TBP ratio and $R^{w}\left(\frac{\mathrm{O}}{\mathrm{A}}\right)$ ratio on the extraction of $\mathrm{Li}^{+}$, were further optimized, and the results are presented in Figure 5. Individually, TBP and IL do not allow for an efficient $\mathrm{Li}^{+}$extraction. However, when used in combination, a synergistic effect is observed, and the extraction efficiency becomes dependent on the TBP content

\begin{figure}
\includegraphics[alt={},max width=\textwidth]{https://cdn.mathpix.com/cropped/a59933d6-bca6-45db-b31f-2c9dd99cc78c-6.jpg?height=590&width=748&top_left_y=204&top_left_x=1136}
\captionsetup{labelformat=empty}
\caption{Figure 5. Effects of the TBP content in the organic phase on extraction efficiency $E_{\mathrm{Li}^{+}}$with [emim][tcb] and TBP ( $w_{\text {TBP, organic }}=$ 0.80 ) as the extractant ( $w_{\mathrm{Li}^{+}}=2000 \mathrm{ppm}$ ), $R^{\mathrm{w}}\left(\frac{\mathrm{O}}{\mathrm{A}}\right)=1, T=294.15 \mathrm{~K}$, $p=1.013$ bar, mixing time 2 h at 2500 rpm . Exp. data are given in Table S10. Standard uncertainty is $u(T)=0.5 \mathrm{~K} ; u(m)=0.0001 \mathrm{~g}$. The combined expanded uncertainty is $U=0.006$ with a coverage factor of $k=2$ ( $95 \%$ confidence interval).}
\end{figure}
in the organic phase. This is due to the difference in $\mathrm{Li}^{+}$ extraction mechanisms between using IL-based extraction system vs nonionic extraction systems. On the one hand, transfer of metal cations through organic extractants into an organic phase usually occurs via ion-pair extraction. On the other hand, cation extraction in IL-based systems underlies an ion-exchange mechanism, which shows superior extraction efficiency. ${ }^{37}$ At a TBP content of $60 \mathrm{wt} \%$ in the organic phase, the extraction efficiency reaches $74 \%$. This efficiency progressively increases with a higher TBP content, peaking at $91 \%$ when the TBP content is $85 \mathrm{wt} \%$ in the organic phase. This trend is attributed to an enhanced likelihood of the complexation reaction (see eq 9) between $\mathrm{Li}^{+}$ions and TBP. The increase in extraction efficiency could also be attributed to the lower concentration of ionic liquid. This decreases the viscosity of the extractant, promoting the diffusion of $\mathrm{Li}^{+}$ions within the organic phase. ${ }^{38}$ Nevertheless, exceeding a certain TBP content leads to a decrease in extraction efficiency. For instance, with a TBP content of $95 \mathrm{wt} \%$ in the organic feed, only a $70 \%$ extraction efficiency of $\mathrm{Li}^{+}$was achieved. The findings are consistent with other IL/TBP extraction mixtures for $\mathrm{Li}^{+}$extraction from the literature. ${ }^{19,39-41}$

The optimization of the organic-to-aqueous phase ratio $R^{\mathrm{w}}\left(\frac{\mathrm{O}}{\mathrm{A}}\right)$ is crucial to ensure resource-efficient extraction. Ideally, the use of organic feed should be minimized while extracting $\mathrm{Li}^{+}$ions from an aqueous source. Figure 6 illustrates the extraction efficiency of $\mathrm{Li}^{+}$at various mass ratios of organic to aqueous phase $R^{\mathrm{w}}\left(\frac{\mathrm{O}}{\mathrm{A}}\right)$ ranging from 0.52 to 3.93 with a TBP content in the organic feed phase of $w_{\text {TBP }}=0.85$.

The results reveal a significant increase in extraction efficiency from $75 \%$ to $89 \%$ as $R^{\mathrm{w}}\left(\frac{\mathrm{O}}{\mathrm{A}}\right)$ increases from 0.52 to 0.99 . Further increasing the percentage of organic feed leads to a gradual increase in extraction efficiency but with a less pronounced slope compared to lower values. The maximum extraction efficiency of $98.33 \%$ was achieved at an $R^{\mathrm{w}}\left(\frac{\mathrm{O}}{\mathrm{A}}\right)=$ 3.93. This is an expected result since TBP excess molecules

\begin{figure}
\includegraphics[alt={},max width=\textwidth]{https://cdn.mathpix.com/cropped/a59933d6-bca6-45db-b31f-2c9dd99cc78c-7.jpg?height=590&width=744&top_left_y=204&top_left_x=223}
\captionsetup{labelformat=empty}
\caption{Figure 6. Effects of $R^{\mathrm{w}}\left(\frac{\mathrm{O}}{\mathrm{A}}\right)$ on experimental extraction efficiency $E_{\mathrm{Li}^{+}}$ with [emim] [tcb] and TBP ( $w_{\text {TBP, organic }}=0.85$ ) as the extractant ( $w_{\mathrm{Li}}{ }^{+} =2000 \mathrm{ppm}$ ), $T=294.15 \mathrm{~K}, p=1.013$ bar, mixing time 2 h at 2500 rpm. Exp. data are given in Table S11. The combined expanded uncertainty is $U=0.01$ with a coverage factor of $k=2(95 \%$ confidence interval).}
\end{figure}
ensure that the reaction equilibrium between TBP and $\mathrm{Li}^{+}$ions (eq 9) is far on the right-hand side. The higher $R^{\mathrm{w}}\left(\frac{\mathrm{O}}{\mathrm{A}}\right)$, the lower the concentration of $\mathrm{Li}^{+}$ions in the organic phase, which makes re-extraction more complex. These findings are in agreement with previous studies using TBP and a different IL. ${ }^{36,37,42,43}$ For reasons of efficiency and economy, a balanced $R^{\mathrm{w}}\left(\frac{\mathrm{O}}{\mathrm{A}}\right)$ value of $1: 1$ is proposed to strike a good compromise between extraction efficiency and resource utilization.

Finally, temperature will also impact extraction efficiency. All results obtained here are valid at 294.15 K . We briefly checked the influence of a lower temperature on the extraction efficiency, and we found that extraction at 275.15 K yielded $E_{\mathrm{Li}^{+}}$of about $2 \%$ higher than at 294.15 K (data not shown). This is an expected result, which will be studied further in the future.
4.3.1. $T B P+[$ emim $][$ tcb $]+L i C I+$ Water. As explained in Section 3, complex formation (eq 9) is the underlying mechanism when TBP and an IL extract $\mathrm{Li}^{+}$ions, which are then transferred to the organic phase. To characterize this complex formation in the ePC-SAFT framework, an additional association site was assigned to $\mathrm{Li}^{+}$ions, allowing a strong interaction between TBP and $\mathrm{Li}^{+}$. We used one experimental data point $\left[E_{\mathrm{Li}^{+}}\right.$at $\left.R^{\mathrm{w}}\left(\frac{\mathrm{O}}{\mathrm{A}}\right)=2\right]$ to adjust the strength of this cross-association. This was done by fitting $k_{i j}^{\mathrm{hb}}$ parameter between TBP and $\mathrm{Li}^{+}$at $R^{\mathrm{w}}\left(\frac{\mathrm{O}}{\mathrm{A}}\right)=2$, cf. Table 4, while all other data points were then obtained by predictive modeling. Note that cross association involving $\mathrm{Li}^{+}$was only allowed for the pair $\mathrm{Li}^{+}$/TBP, while it was forbidden for $\mathrm{Li}^{+}$/water and $\mathrm{Li}^{+}$/IL. Figure 7 demonstrates that the predictions obtained using ePC-SAFT advanced are in qualitative agreement with the experimental data. Notably, ePC-SAFT advanced accurately predicts the increase in extraction efficiency as $R^{\mathrm{w}}\left(\frac{\mathrm{O}}{\mathrm{A}}\right)$ increases, converging to a maximum value of $98.3 \%$ at $R^{\mathrm{w}}\left(\frac{\mathrm{O}}{\mathrm{A}}\right)=3.93$. A comprehensive comparison between the predicted modeling results presented in this study and the

\begin{figure}
\includegraphics[alt={},max width=\textwidth]{https://cdn.mathpix.com/cropped/a59933d6-bca6-45db-b31f-2c9dd99cc78c-7.jpg?height=683&width=839&top_left_y=202&top_left_x=1095}
\captionsetup{labelformat=empty}
\caption{Figure 7. Effect of $R^{\mathrm{w}}(\mathrm{O} / \mathrm{A})$ on the extraction efficiency $E_{\mathrm{Li}^{+}}$with [emim][tcb] and TBP as the extractant ( $w_{\text {TBP, organic }}=0.85, w_{\mathrm{Li}}{ }^{+}=$ 2000 ppm ), $T=294.15 \mathrm{~K}, p=1.013$ bar, mixing time 2 h at 2500 rpm. Line: ePC-SAFT using the parameters from Tables 2-4. Symbols: Exp. data according to Table S11.}
\end{figure}
experimental data is provided in the Supporting Information Table S11. Remarkably, no adjustments were necessary for any binary interaction parameters between water and ions or between [emim][tcb] and TBP.

\section*{5. CONCLUSIONS}

This study focused on the capability of tetracyanoborate-based ILs for $\mathrm{Li}^{+}$extraction from aqueous solutions. The significance of this work is underscored by the increasing demand for lithium in various applications, including battery production for electronic devices and electric vehicles. With the global demand for lithium projected to continue rising, the development of efficient and environmentally friendly extraction methods is of paramount importance.

The study delves into the potential of tetracyanoboratebased ILs, particularly [emim][tcb], as promising coextraction agents for $\mathrm{Li}^{+}$extraction aided by TBP. Through a combination of experimental work and ePC-SAFT modeling, we systematically examined the various factors influencing the extraction process. These include the choice of IL, the TBP content in the organic phase, the organic-to-aqueous phase ratio $\left(R^{\mathrm{w}}\left(\frac{\mathrm{O}}{\mathrm{A}}\right)\right)$, and the $\mathrm{Li}^{+}$concentration in the aqueous feed. Notably, the results indicate that [emim][tcb] is a suitable candidate for lithium extraction, with an extraction efficiency exceeding 98\% achieved under optimized conditions, which is comparable to the most promising extractants proposed in the literature to date. ${ }^{44}$

The highest extraction efficiency for a one-step extraction process assisted by the combining of IL and TBP in the literature was measured with $\mathrm{TBP}+[\mathrm{Na}]\left[\mathrm{NTf}_{2}\right]$ with $E_{\mathrm{Li}^{+}}= 98 \% .^{12}$ This is a very good result considering that even $\mathrm{Mg}^{2+}$ was present in the aqueous feed, which is known to be disadvantageous for $\mathrm{Li}^{+}$extraction. Higher extraction efficiencies were only achieved with a four-step extraction in the work of Li et al. ${ }^{45}$ with a value of $E_{\mathrm{Li}^{+}}=99.47 \%$ using the extraction system TBP + [bmim][bph4]. Based on our results, the next steps must be to analyze the tetracyanoborate-based ILs for their selectivity in terms of $\mathrm{Mg}^{2+} / \mathrm{Li}^{+}$separation and temper-
ature-dependency and to investigate the regeneration of the extraction solvent, which is necessary to develop a sustainable and scalable process. The IL cation transferred to the aqueous phase during the ion exchange mechanism should also be considered.

This work also contributes to the field of thermodynamic modeling by providing insights into the phase equilibria of mixtures containing ILs, LiCl , water, and TBP and the underlying equilibria. The ePC-SAFT modeling results were in quantitative agreement with experimental data, further validating the reliability of the predictions.

In summary, this research offers an exploration of tetracyanoborate-based ILs as efficient coextraction agents for $\mathrm{Li}^{+}$. The findings provide insights into the optimization of extraction processes and the thermodynamic modeling of liquid-liquid extraction of $\mathrm{Li}^{+}$. These insights are crucial for the development of sustainable and economically viable $\mathrm{Li}^{+}$ extraction methods, which are essential for meeting the increasing demand for lithium in various technological applications.

\section*{ASSOCIATED CONTENT}

\section*{(I) Supporting Information}

The Supporting Information is available free of charge at https://pubs.acs.org/doi/10.1021/acs.jced.4c00074.

Experimental data sets (PDF)

\section*{- AUTHOR INFORMATION}

\section*{Corresponding Author}

Christoph Held - Laboratory of Thermodynamics, TU
Dortmund University, 44227 Dortmund, Germany;
© orcid.org/0000-0003-1074-177X; Phone: +49-
2317552086; Email: christoph.held@tu-dortmund.de

\section*{Authors}

Tobias Hubach - Laboratory of Thermodynamics, TU Dortmund University, 44227 Dortmund, Germany
Zeynep Er - Laboratory of Thermodynamics, TU Dortmund University, 44227 Dortmund, Germany
Complete contact information is available at:
https://pubs.acs.org/10.1021/acs.jced.4c00074

\section*{Author Contributions}
T.H.: methodology, software, validation, formal analysis, data curation, investigation, writing-original draft preparation, and visualization. Z.E.: methodology and investigation. C.H.: resources, supervision, writing-review and editing, project administration, and funding acquisition.

\section*{Funding}

This research received no external funding.

\section*{Notes}

The authors declare no competing financial interest.

\section*{- ACKNOWLEDGMENTS}
C.H. wants to thank Pierre Walker for pointing out the important notice that it is required to use the temperaturecorrected ion diameters in all ePC-SAFT publications in order to achieve the results obtained by our group. This information is crucial, and it was obvious to me since DH and Born theories describe the temperature-dependent interactions of ions, for which the temperature treatment within the PC-SAFT framework has also been applied to ions throughout all our
works. This information should be taken into account by other researchers who want to rebuild a code based on ePC-SAFT.

\section*{- REFERENCES}
(1) Schmidt, M. Rohstoffrisikobewertung-Lithium, 2022nd ed.; DERA: Berlin, 2023.
(2) Mineral Commodity Summaries 2022, 2022.
(3) Zhang, Y.; Sun, W.; Xu, R.; Wang, L.; Tang, H. Lithium extraction from water lithium resources through green electro-chemical-battery approaches: A comprehensive review. J. Clean. Prod. 2021, 285, 124905.
(4) Kumar, A.; Naidu, G.; Fukuda, H.; Du, F.; Vigneswaran, S.; Drioli, E.; Lienhard, J. H. Metals Recovery from Seawater Desalination Brines: Technologies, Opportunities, and Challenges. ACS Sustainable Chem. Eng. 2021, 9, 7704-7712.
(5) Li, L.; Deshmane, V. G.; Paranthaman, M. P.; Bhave, R.; Moyer, B. A.; Harrison, S. Lithium Recovery from Aqueous Resources and Batteries: A Brief Review. Johnson Matthey Technol. Rev. 2018, 62, 161-176.
(6) Vera, M. L.; Torres, W. R.; Galli, C. I.; Chagnes, A.; Flexer, V. Environmental impact of direct lithium extraction from brines. Nat. Rev. Earth Environ. 2023, 4, 149-165.
(7) Lai, X.; Yuan, Y.; Chen, Z.; Peng, J.; Sun, H.; Zhong, H. Adsorption-Desorption Properties of Granular EP/HMO Composite and Its Application in Lithium Recovery from Brine. Ind. Eng. Chem. Res. 2020, 59, 7913-7925.
(8) Battaglia, G.; Berkemeyer, L.; Cipollina, A.; Cortina, J. L.; Fernandez de Labastida, M.; Lopez Rodriguez, J.; Winter, D. Recovery of Lithium Carbonate from Dilute Li-Rich Brine via Homogenous and Heterogeneous Precipitation. Ind. Eng. Chem. Res. 2022, 61, 13589-13602.
(9) Li, X.; Mo, Y.; Qing, W.; Shao, S.; Tang, C. Y.; Li, J. Membranebased technologies for lithium recovery from water lithium resources: A review. J. Membr. Sci. 2019, 591, 117317.
(10) Hubach, T.; Pillath, M.; Knaup, C.; Schlüter, S.; Held, C. Li+ Separation from Multi-Ionic Mixtures by Nanofiltration Membranes: Experiments and Modeling. Modelling 2023, 4, 408-425.
(11) Zhou, W.; Xu, S.; Li, Z. Recovery of Lithium from Brine with a High Mg/Li Ratio Using Hydroxyl-Functionalized Ionic Liquid and Tri-n-butyl Phosphate. J. Sustain. Metall. 2021, 7, 256-265.
(12) Bai, R.; Wang, J.; Wang, D.; Cui, J.; Zhang, Y. Recovery of lithium from high $\mathrm{Mg} / \mathrm{Li}$ ratio salt-lake brines using ion-exchange with NaNTf2 and TBP. Hydrometallurgy 2022, 213, 105914.
(13) Zhou, Z.; Qin, W.; Fei, W. Extraction Equilibria of Lithium with Tributyl Phosphate in Three Diluents. J. Chem. Eng. Data 2011, 56, 3518-3522.
(14) Zhou, Z.; Qin, W.; Liu, Y.; Fei, W. Extraction Equilibria of Lithium with Tributyl Phosphate in Kerosene and FeCl 3. J. Chem. Eng. Data 2012, 57, 82-86.
(15) Bai, R.; Wang, J.; Wang, D.; Zhang, Y.; Cui, J. Selective separation of lithium from the high magnesium brine by the extraction system containing phosphate-based ionic liquids. Sep. Purif. Technol. 2021, 274, 119051.
(16) Shi, C.; Duan, D.; Jia, Y.; Jing, Y. A highly efficient solvent system containing ionic liquid in tributyl phosphate for lithium ion extraction. J. Mol. Liq. 2014, 200, 191-195.
(17) Shi, C.; Jing, Y.; Jia, Y. Solvent extraction of lithium ions by tri-n-butyl phosphate using a room temperature ionic liquid. J. Mol. Liq. 2016, 215, 640-646.
(18) Li, Z.; Binnemans, K. Opposite selectivities of tri- n -butyl phosphate and Cyanex 923 in solvent extraction of lithium and magnesium. AIChE J. 2021, 67, No. e17219.
(19) Bai, R.; Wang, J.; Cui, L.; Yang, S.; Qian, W.; Cui, P.; Zhang, Y. Efficient Extraction of Lithium Ions from High Mg/Li Ratio Brine through the Synergy of TBP and Hydroxyl Functional Ionic Liquids. Chin. J. Chem. 2020, 38, 1743-1751.
(20) Zhao, X.; Wu, H.; Duan, M.; Hao, X.; Yang, Q.; Zhang, Q.; Huang, X. Liquid-liquid extraction of lithium from aqueous solution
using novel ionic liquid extractants via COSMO-RS and experiments. Fluid Phase Equilib. 2018, 459, 129-137.
(21) Yu, G.; Zhang, X.; Hubach, T.; Chen, B.; Held, C. Highly efficient lithium extraction from magnesium-rich brines with ionic liquid-based collaborative extractants: Thermodynamics and molecular insights. Chem. Eng. Sci. 2024, 286, 119682.
(22) Bezdomnikov, A. A.; Demina, L. I.; Kuz'mina, L. G.; Kostikova, G. V.; Zhilov, V. I.; Tsivadze, A. Y. Study of Lithium-Extraction Systems Based on Benzo-15-Crown-5 Ether and AlkylimidazoliumBased Ionic Liquid. Molecules 2023, 28, 935.
(23) Masmoudi, A.; Zante, G.; Trébouet, D.; Barillon, R.; Boltoeva, M. Solvent extraction of lithium ions using benzoyltrifluoroacetone in new solvents. Sep. Purif. Technol. 2021, 255, 117653.
(24) Chen, S.; Gao, D.; Yu, X.; Guo, Y.; Deng, T. Thermokinetics of lithium extraction with the novel extraction systems (tri-isobutyl phosphate + ionic liquid + kerosene). J. Chem. Thermodyn. 2018, 123, 79-85.
(25) Mota-Martinez, M. T.; Althuluth, M.; Kroon, M. C.; Peters, C. J. Solubility of carbon dioxide in the low-viscosity ionic liquid 1-hexyl-3-methylimidazolium tetracyanoborate. Fluid Phase Equilib. 2012, 332, 35-39.
(26) Freire, M. G.; Neves, C. M. S. S.; Marrucho, I. M.; Coutinho, J. A. P.; Fernandes, A. M. Hydrolysis of tetrafluoroborate and hexafluorophosphate counter ions in imidazolium-based ionic liquids. J. Phys. Chem. A 2010, 114, 3744-3749.
(27) Binnemans, K.; Jones, P. T. Ionic Liquids and Deep-Eutectic Solvents in Extractive Metallurgy: Mismatch Between Academic Research and Industrial Applicability. J. Sustain. Metall. 2023, 9, 423438.
(28) Gross, J.; Sadowski, G. Perturbed-Chain SAFT: An Equation of State Based on a Perturbation Theory for Chain Molecules. Ind. Eng. Chem. Res. 2001, 40, 1244-1260.
(29) Ascani, M.; Held, C. Prediction of salting-out in liquid-liquid two-phase systems with ePC-SAFT: Effect of the Born term and of a concentration-dependent dielectric constant. Z. Anorg. Allg. Chem. 2021, 647, 1305-1314.
(30) Bülow, M.; Ji, X.; Held, C. Incorporating a concentrationdependent dielectric constant into ePC-SAFT, An application to binary mixtures containing ionic liquids. Fluid Phase Equilib. 2019, 492, 26-33.
(31) Andeen, C.; Fontanella, J.; Schuele, D. Low-Frequency Dielectric Constant of $\mathrm{LiF}, \mathrm{NaF}, \mathrm{NaCl}, \mathrm{NaBr}, \mathrm{KCl}$, and KBr by the Method of Substitution. Phys. Rev. B 1970, 2, 5068-5073.
(32) Tian, Q.; Liu, H. Densities and Viscosities of Binary Mixtures of Tributyl Phosphate with Hexane and Dodecane from (298.15 to 328.15) K. J. Chem. Eng. Data 2007, 52, 892-897.
(33) Hiraga, Y.; Hagiwara, S.; Sato, Y.; Smith, R. L. Measurement and Correlation of High-Pressure Densities and Atmospheric Viscosities of Ionic Liquids: 1-Butyl-1-methylpyrrolidinium Bis(trifluoromethylsulfonyl)imide, 1-Allyl-3-methylimidazolium Bis(trifluoromethylsulfonyl)imide, 1-Ethyl-3-methylimidazolium Tetracyanoborate, and 1-Hexyl-3-methylimidazolium Tetracyanoborate. J. Chem. Eng. Data 2018, 63, 972-980.
(34) Held, C.; Reschke, T.; Mohammad, S.; Luza, A.; Sadowski, G. ePC-SAFT revised. Chem. Eng. Res. Des. 2014, 92, 2884-2897.
(35) Held, C.; Cameretti, L. F.; Sadowski, G. Modeling aqueous electrolyte solutions. Fluid Phase Equilib. 2008, 270, 87-96.
(36) Wang, X.; Jing, Y.; Liu, H.; Yao, Y.; Shi, C.; Xiao, J.; Wang, S.; Jia, Y. Extraction of lithium from salt lake brines by bis-[(trifluoromethyl)sulfonyl]imide-based ionic liquids. Chem. Phys. Lett. 2018, 707, 8-12.
(37) Shi, C.; Jing, Y.; Jia, Y. Tri-n-butyl phosphate-ionic liquid mixtures for $\mathrm{Li}+$ extraction from $\mathrm{Mg} 2+$-containing brines at $303-343$ K. Russ. J. Phys. Chem. 2017, 91, 692-696.
(38) Zheng, H.; Dong, T.; Sha, Y.; Jiang, D.; Zhang, H.; Zhang, S. Selective Extraction of Lithium from Spent Lithium Batteries by Functional Ionic Liquid. ACS Sustainable Chem. Eng. 2021, 9, 70227029.
(39) Yang, S.; Liu, G.; Wang, J.; Cui, L.; Chen, Y. Recovery of lithium from alkaline brine by solvent extraction with functionalized ionic liquid. Fluid Phase Equilib. 2019, 493, 129-136.
(40) Cui, L.; Jiang, K.; Wang, J.; Dong, K.; Zhang, X.; Cheng, F. Role of ionic liquids in the efficient transfer of lithium by Cyanex 923 in solvent extraction system. AIChE J. 2019, 65, No. e16606.
(41) Masmoudi, A.; Zante, G.; Trébouet, D.; Barillon, R.; Boltoeva, M. Understanding the Mechanism of Lithium Ion Extraction Using Tributyl Phosphate in Room Temperature Ionic Liquid. Solvent Extr. Ion Exch. 2020, 38, 777-799.
(42) Gao, D.; Guo, Y.; Yu, X.; Wang, S.; Deng, T. Extracting Lithium from the High Concentration Ratio of Magnesium and Lithium Brine Using Imidazolium-Based Ionic Liquids with Varying Alkyl Chain Lengths. J. Chem. Eng. Jpn. 2016, 49, 104-110.
(43) Gao, D.; Yu, X.; Guo, Y.; Wang, S.; Liu, M.; Deng, T.; Chen, Y.; Belzile, N. Extraction of lithium from salt lake brine with triisobutyl phosphate in ionic liquid and kerosene. Chem. Res. Chin. Univ. 2015, 31, 621-626.
(44) Khalil, A.; Mohammed, S.; Hashaikeh, R.; Hilal, N. Lithium recovery from brine: Recent developments and challenges. Desalination 2022, 528, 115611.
(45) Li, R.; Wang, W.; Wang, Y.; Wei, X.; Cai, Z.; Zhou, Z. Novel ionic liquid as co-extractant for selective extraction of lithium ions from salt lake brines with high $\mathrm{Mg} / \mathrm{Li}$ ratio. Sep. Purif. Technol. 2021, 277, 119471.

\begin{figure}
\includegraphics[alt={},max width=\textwidth]{https://cdn.mathpix.com/cropped/a59933d6-bca6-45db-b31f-2c9dd99cc78c-9.jpg?height=1129&width=216&top_left_y=1503&top_left_x=1095}
\captionsetup{labelformat=empty}
\caption{CAS BIOFINDER DISCOVERY PLATFORM ${ }^{\text {TM }}$
STOP DIGGING THROUGH DATA -START MAKING DISCOVERIES
CAS BioFinder helps you find the right biological insights in seconds}
\end{figure}

\section*{Start your search}