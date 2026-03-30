# Upgrading the Extraction of Ethanol and Isobutanol from an Aqueous Solution in the Presence of Sodium Chloride with Molecular Interaction Studies 

Salal Hasan Khudaida, Serli Dwi Rahayu, Paul Figiel, Christoph Held,* and Ardila Hayu Tiwikrama*

Cite This: J. Chem. Eng. Data 2026, 71, 708-717
Read Online
Downloaded via BRIGHAM YOUNG UNIV on March 20, 2026 at 20:56:05 (UTC). See https://pubs.acs.org/sharingguidelines for options on how to legitimately share published articles.


#### Abstract

Salting-out extraction systems composed of alcohols and salts are widely used for separating valuable chemicals and fuels from aqueous solutions. Therefore, in this study, the phase equilibria of the water + ethanol + isobutyl alcohol system in the presence of NaCl were investigated. Liquidliquid equilibrium (LLE) for the quaternary system (water + ethanol + isobutanol +NaCl ) was measured at 293.15 and 313.15 K and atmospheric pressure using a jacketed equilibrium cell. The effect of NaCl concentration on LLE was studied at 5 and $10 \mathrm{wt} \%$. Increasing NaCl enhanced the heterogeneous region and promoted isobutanol separation from aqueous ethanol solutions, with higher concentrations intensifying the salting-out effect and phase splitting. The electrolyte Perturbed-Chain Statistical Associating ![](https://cdn.mathpix.com/cropped/61a9ddb4-14df-4553-8a5c-b608ccc464f6-01.jpg?height=409&width=675&top_left_y=886&top_left_x=1253) Fluid Theory (ePC-SAFT) and eNRTL models were used to satisfactorily correlate the experimental LLE data. Regression of binary interaction parameters (BIPs) showed that eNRTL required 13 parameters, whereas ePC-SAFT required only 7 parameters, demonstrating that ePC-SAFT efficiently predicts phase behavior with fewer parameters. Additionally, molecular interactions among water, ethanol, isobutyl alcohol, and NaCl were analyzed using sigma profiles. The results highlight the practical applicability of LLE data for designing separation processes and show that ePC-SAFT provides accurate predictions for electrolyte-aqueous alcohol systems when experimental data are available.


## 1. INTRODUCTION

Liquid-liquid equilibrium (LLE) and the salting-out effect have been broadly studied as required knowledge in chemical engineering processes, biofuels, wastewater treatment, and pharmaceuticals because of their industrial relevance and economic advantages. ${ }^{1-9}$ Salting-out extraction systems composed of alcohols and salts are also commonly employed in many industrial processes for the separation of useful chemicals and fuels from aqueous solutions. ${ }^{10,11}$

The addition of electrolytes to liquid mixtures can significantly influence the mutual solubility of components through the salting-out effect, thereby enhancing the phase separation and improving the extraction efficiency. For example, Pai and Rao ${ }^{12}$ investigated the salt effect on the liquid-liquid equilibria of the water + ethanol + ethyl acetate system in the presence of saturated solutions of potassium acetate (KAc), sodium acetate (NaAc), potassium chloride $(\mathrm{KCl})$, sodium chloride $(\mathrm{NaCl})$, potassium sulfate $\left(\mathrm{K}_{2} \mathrm{SO}_{4}\right)$, and sodium sulfate $\left(\mathrm{Na}_{2} \mathrm{SO}_{4}\right)$ at 303.15 K . Their results indicated that $\mathrm{K}_{2} \mathrm{SO}_{4}$ caused the greatest solubility depression of water in the organic phase. Lin et al. ${ }^{13}$ studied the influence of adding electrolytes such as KAc or calcium chloride ( $\mathrm{CaCl}_{2}$ ) on the phase behavior of the water + ethanol + ethyl acetate mixture at temperatures ranging from 283.15 to 313.15 K .

They found that liquid-liquid phase splitting was significantly enhanced by the addition of these salts. Xie et al. ${ }^{10}$ examined the effect of several electrolytes-potassium carbonate $\left(\mathrm{K}_{2} \mathrm{CO}_{3}\right)$, dipotassium phosphate $\left(\mathrm{K}_{2} \mathrm{HPO}_{4}\right)$, tripotassium phosphate $\left(\mathrm{K}_{3} \mathrm{PO}_{4}\right)$, and potassium pyrophosphate $\left(\mathrm{K}_{4} \mathrm{P}_{2} \mathrm{O}_{7}\right)$ on the water + ethanol system at 298 K . Their results demonstrated improved separation efficiency induced by the electrolytes: by salting-out, salts such as $\mathrm{K}_{4} \mathrm{P}_{2} \mathrm{O}_{7}$ promote an almost complete separation of ethanol from water, leaving the ethanol nearly salt-free and concentrating the salts in the water-rich phase.

Although phase equilibrium data for alcohol-water systems with and without salts have been reported in the literature ${ }^{10,14-19}$ as summarized in Table 1, experimental data for quaternary systems containing standard salts (here: NaCl ) remain limited. To the best of our knowledge, no liquid-liquid equilibrium (LLE) data are available for the

[^0]![](https://cdn.mathpix.com/cropped/61a9ddb4-14df-4553-8a5c-b608ccc464f6-01.jpg?height=216&width=151&top_left_y=2342&top_left_x=1781)

Table 1. Phase Equilibrium Data of Water/Alcohol Systems with and without Salts Reported in the Literature
| systems | phase behavior | T/K | references |
| :--- | :--- | :--- | :--- |
| water + ethanol | VLE | 303-363 | 14 |
| ethanol + isobutanol | VLE | 351-381 | 19 |
| water + ethanol $+\mathrm{K}_{2} \mathrm{CO}_{3}$ | LLE | 298.15 | 10 |
| water + ethanol $+\mathrm{K}_{2} \mathrm{HPO}_{4}$ | LLE | 298.15 | 10 |
| water + ethanol $+\mathrm{K}_{3} \mathrm{PO}_{4}$ | LLE | 298.15 | 10 |
| water + ethanol $+\mathrm{K}_{4} \mathrm{P2O}_{7}$ | LLE | 298.15 | 10 |
| water + ethanol $+\mathrm{CaCl}_{2}$ | VLE | 298.15 | 16 |
| water + ethanol +KAc | VLE | 354-380 | 17 |
| water + isobutanol +NaCl | LLE | 298.13 | 15 |
|  | LLE | 298.15 | 18 |


water + ethanol + isobutanol +NaCl system, which is investigated for the first time in the present work. In addition, for mixed-solvent electrolyte solutions, the experimental data from various sources might be inconsistent. Therefore, extensive data collection and new experimental measurements are necessary to establish a reliable database, improve thermodynamic models, and support the development of separation processes using salting-out effects. ${ }^{20}$ The objective of the present study is to investigate the LLE behavior of the water + ethanol + isobutanol system with 5 and $10 \mathrm{wt} \% \mathrm{NaCl}$, at temperatures ranging from 293.15 to 313.15 K under atmospheric pressure. The experimental LLE tie-line data for the quaternary system (water + ethanol + isobutanol +NaCl ) were correlated using the electrolyte Perturbed-Chain Statistical Associating Fluid Theory (ePC-SAFT) ${ }^{21}$ and eNRTL ${ }^{22}$ models. Additionally, the molecular interaction of water, ethanol, isobutanol, and NaCl was investigated by means of sigma profiles.

## 2. EXPERIMENTAL SECTION

2.1. Materials. Ethanol (purity $>0.998$ mass fraction) and sodium chloride ( NaCl , purity $>0.995$ mass fraction) were obtained from Fischer Chemical, while isobutanol (purity $>0.99$ mass fraction) was provided by Thermo Scientific. All chemicals were used as received from the suppliers. Deionized ultrapure water (resistivity: $18.3 \mathrm{M} \Omega \cdot \mathrm{cm}$ ) was prepared in our laboratory using a NANO Pure-Ultrapure water system. Table 2 summarizes the materials used in this study.
2.2. Experimental Procedure. The tie-lines for LLE were measured using a glass-made equilibrium cell with an internal volume of approximately $25 \mathrm{~cm}^{3}$ at temperatures ranging from
293.15 to 313.15 K and at atmospheric pressure ( $P=0.1 \mathrm{MPa}$ ) as described in our previous works. ${ }^{23,24}$ Thermostatic water from a circulating water bath was used to maintain the cell at a constant temperature with a standard uncertainty, $u(T)$, of $\pm 0.1 \mathrm{~K}$. The temperature inside the cell was monitored with a precision thermometer (Model 1560, Hart Scientific Co.) with an uncertainty of $\pm 0.02 \mathrm{~K}$. The liquid mixture was vigorously agitated with a magnetic stirrer. Different feed compositions of water + ethanol +NaCl were investigated. The aqueous feed was mixed with the extraction solvent (isobutanol) in a mass ratio of $1: 1$. The ethanol content in the aqueous feed was varied from 2 to $12 \mathrm{wt} \%$ at a fixed NaCl concentration of 5 wt $\%$ and from 1 to $5 \mathrm{wt} \%$ at a fixed NaCl concentration of 10 wt $\%$ (based on the total feed). The mixture of water + ethanol + isobutanol +NaCl was stirred in the equilibrium cell for 5 h and then allowed to settle for at least 12 h at 298.15 and 313.15 K until complete phase separation was achieved. Samples were then withdrawn from the two phases for analysis. The organic-rich phase was collected from the top opening of the cell using a syringe, while the water-rich phase was sampled from the bottom port to avoid cross-contamination. The equilibrium compositions were determined from three replicate samples of each phase.
2.3. Analytical Measurements. Samples from the organic (upper) and aqueous (lower) phases were analyzed using a gas chromatograph (GC, Model 1000, China Chromatography Co., Taiwan) equipped with a thermal conductivity detector (TCD). A stainless-steel column packed with Porapak Qs 80/ 100 ( 4 m length) was used, with high-purity helium ( $99.99 \%$ ) as the carrier gas, to determine the equilibrium compositions of water, ethanol, and isobutanol in both phases. To prevent salt from entering the GC column, a stainless-steel tube filled with glass wool was installed at the column inlet and replaced when noticeable salt accumulation occurred. The injector and detector temperatures were set to $210{ }^{\circ} \mathrm{C}$; the oven temperature was maintained at $200^{\circ} \mathrm{C}$, and the injection volume was $1 \mu \mathrm{~L}$. Calibration curves were established before analysis. The NaCl concentrations in the organic- and aqueous-rich phases were determined gravimetrically. Approximately 1 mL of the liquid from each phase was divided into three vials. Each sample was evaporated in an oven at $105^{\circ} \mathrm{C}$ until the NaCl mass remained constant. Three replicate measurements were performed for each phase.

Table 2. Description of Chemicals Used in This Study
| Materials | Structure | CAS Number | Source | Purity in mass fraction ${ }^{\mathrm{a}}$ | MW (g/mol) |
| :--- | :--- | :--- | :--- | :--- | :--- |
| Ethanol | <smiles>CCO</smiles> | 64-17-5 | Fisher Chemical | > 0.998 | $46.07{ }^{\text {a }}$ |
| Isobutanol | <smiles>CC(C)CO</smiles> | 78-83-1 | Thermo Scientific | > 0.99 | $74.12^{\mathrm{a}}$ |
| Sodium Chloride | $\mathrm{Na}^{+} \mathrm{Cl}^{-}$ | 7647-14-5 | Fisher Chemical | > 0.995 | $58.44^{\mathrm{a}}$ |
| Water ( $\mathrm{H}_{2} \mathrm{O}$ ) | <smiles>O</smiles> | 7732-18-5 | Prepared in our Lab. | Deionized ultrapure water ${ }^{\text {b }}$ | $18.015{ }^{\text {a }}$ |


[^1]Table 3. Experimental LLE Tie-Line Data for Water (1) + Ethanol (2) + Isobutanol (3) + $\mathbf{5 ~ w t} \boldsymbol{\% ~ N a C l}$ (4) System at Temperatures ( $T=293.15$ and 313.15 K ) and Atmospheric Pressure ( $P=0.1 \mathrm{MPa})^{a}$
| $T(\mathrm{~K})$ | organic phase ( $x_{i}$ ) |  |  | $x_{4}{ }^{\text {I }}$ | aqueous phase ( $x_{i}$ ) |  |  | $x_{4}{ }^{\mathrm{II}}$ |
| :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- |
|  | $x_{1}{ }^{\text {I }}$ | $x_{2}{ }^{\mathrm{I}}$ | $x_{3}{ }^{\mathrm{I}}$ |  | $x_{1}{ }^{\mathrm{II}}$ | $x_{2}{ }^{\mathrm{II}}$ | $x_{3}{ }^{\mathrm{II}}$ |  |
| 293.15 | 0.3705 | 0.1624 | 0.4650 | 0.0021 | 0.9331 | 0.0159 | 0.0028 | 0.0482 |
|  | 0.3703 | 0.1400 | 0.4880 | 0.0017 | 0.9404 | 0.0121 | 0.0027 | 0.0448 |
|  | 0.3691 | 0.1144 | 0.5153 | 0.0012 | 0.9446 | 0.0110 | 0.0022 | 0.0422 |
|  | 0.3695 | 0.0927 | 0.5366 | 0.0012 | 0.9486 | 0.0092 | 0.0037 | 0.0385 |
|  | 0.3692 | 0.0759 | 0.5537 | 0.0012 | 0.9509 | 0.0073 | 0.0034 | 0.0384 |
|  | 0.3685 | 0.0646 | 0.5658 | 0.0011 | 0.9539 | 0.0056 | 0.0034 | 0.0371 |
|  | 0.3676 | 0.0478 | 0.5836 | 0.0010 | 0.9567 | 0.0044 | 0.0036 | 0.0353 |
|  | 0.3641 | 0.0336 | 0.6015 | 0.0008 | 0.9580 | 0.0030 | 0.0037 | 0.0353 |
| 303.15 | 0.3906 | 0.1597 | 0.4480 | 0.0017 | 0.9345 | 0.0149 | 0.0025 | 0.0481 |
|  | 0.3874 | 0.1371 | 0.4740 | 0.0015 | 0.9397 | 0.0124 | 0.0029 | 0.0450 |
|  | 0.3887 | 0.1105 | 0.4995 | 0.0013 | 0.9437 | 0.0099 | 0.0036 | 0.0428 |
|  | 0.3849 | 0.0885 | 0.5258 | 0.0008 | 0.9474 | 0.0081 | 0.0044 | 0.0401 |
|  | 0.3821 | 0.0759 | 0.5410 | 0.0010 | 0.9489 | 0.0069 | 0.0049 | 0.0393 |
|  | 0.3768 | 0.0623 | 0.5601 | 0.0008 | 0.9514 | 0.0055 | 0.0050 | 0.0381 |
|  | 0.3753 | 0.0476 | 0.5767 | 0.0004 | 0.9550 | 0.0043 | 0.0041 | 0.0366 |
|  | 0.3740 | 0.0316 | 0.5940 | 0.0004 | 0.9572 | 0.0026 | 0.0046 | 0.0356 |
| 313.15 | 0.4154 | 0.1313 | 0.4513 | 0.0020 | 0.9406 | 0.0122 | 0.0037 | 0.0435 |
|  | 0.4130 | 0.1094 | 0.4759 | 0.0017 | 0.9426 | 0.0101 | 0.0041 | 0.0432 |
|  | 0.4077 | 0.0866 | 0.5045 | 0.0012 | 0.9457 | 0.0079 | 0.0043 | 0.0421 |
|  | 0.4091 | 0.0701 | 0.5192 | 0.0016 | 0.9497 | 0.0064 | 0.0042 | 0.0397 |
|  | 0.4079 | 0.0577 | 0.5333 | 0.0011 | 0.9520 | 0.0051 | 0.0044 | 0.0385 |
|  | 0.4013 | 0.0487 | 0.5489 | 0.0011 | 0.9536 | 0.0043 | 0.0045 | 0.0376 |
|  | 0.3954 | 0.0342 | 0.5694 | 0.0010 | 0.9560 | 0.0029 | 0.0046 | 0.0365 |


${ }^{a}$ Standard uncertainties are $u(T)=0.02 \mathrm{~K}, u(P)=0.001 \mathrm{MPa}, u\left(x_{i}\right)=0.0016$ (water), 0.0007 (ethanol), 0.0015 (isobutanol), and $0.0002(\mathrm{NaCl})$.

Table 4. Experimental LLE Tie-Line Data for Water (1) + Ethanol (2) + Isobutanol (3) + $\mathbf{1 0} \mathbf{~ w t} \boldsymbol{\%} \mathbf{~ N a C l}$ (4) System at Temperatures $(T=293.15-313.15 \mathrm{~K})$ and Atmospheric Pressure $(P=0.1 \mathrm{MPa})^{a}$
| T (K) | organic phase ( $x_{i}$ ) |  |  | $x_{4}{ }^{\text {I }}$ | aqueous phase ( $x_{i}$ ) |  |  | $x_{4}{ }^{\mathrm{II}}$ |
| :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- |
|  | $x_{1}{ }^{\mathrm{I}}$ | $x_{2}{ }^{\mathrm{I}}$ | $x_{3}{ }^{\mathrm{I}}$ |  | $x_{1}{ }^{\mathrm{II}}$ | $x_{2}{ }^{\text {II }}$ | $x_{3}{ }^{\mathrm{II}}$ |  |
| 293.15 | 0.2456 | 0.0919 | 0.6612 | 0.0013 | 0.9223 | 0.0044 | 0.0013 | 0.0720 |
|  | 0.2480 | 0.0767 | 0.6739 | 0.0014 | 0.9253 | 0.0033 | 0.0010 | 0.0704 |
|  | 0.2476 | 0.0615 | 0.6896 | 0.0013 | 0.9277 | 0.0025 | 0.0013 | 0.0685 |
|  | 0.2543 | 0.0383 | 0.7062 | 0.0012 | 0.9311 | 0.0014 | 0.0014 | 0.0661 |
|  | 0.2560 | 0.0182 | 0.7248 | 0.0010 | 0.9353 | 0.0003 | 0.0001 | 0.0643 |
| 303.15 | 0.2587 | 0.0885 | 0.6521 | 0.0007 | 0.9235 | 0.0049 | 0.0015 | 0.0701 |
|  | 0.2566 | 0.0735 | 0.6593 | 0.0006 | 0.9249 | 0.0038 | 0.0016 | 0.0697 |
|  | 0.2632 | 0.0554 | 0.6811 | 0.0003 | 0.9280 | 0.0025 | 0.0016 | 0.0679 |
|  | 0.2681 | 0.0359 | 0.6957 | 0.0003 | 0.9309 | 0.0014 | 0.0012 | 0.0665 |
|  | 0.2736 | 0.0179 | 0.7078 | 0.0007 | 0.9363 | 0.0005 | 0.0003 | 0.0629 |
| 313.15 | 0.2684 | 0.1031 | 0.6268 | 0.0017 | 0.9210 | 0.0057 | 0.0013 | 0.0720 |
|  | 0.2669 | 0.0954 | 0.6358 | 0.0019 | 0.9242 | 0.0042 | 0.0008 | 0.0708 |
|  | 0.2693 | 0.0759 | 0.6535 | 0.0013 | 0.9271 | 0.0028 | 0.0008 | 0.0694 |
|  | 0.2695 | 0.0615 | 0.6677 | 0.0013 | 0.9283 | 0.0022 | 0.0013 | 0.0682 |
|  | 0.2733 | 0.0408 | 0.6845 | 0.0014 | 0.9295 | 0.0015 | 0.0014 | 0.0676 |
|  | 0.2831 | 0.0176 | 0.6982 | 0.0011 | 0.9333 | 0.0005 | 0.0010 | 0.0652 |


${ }^{a}$ Standard uncertainties are $u(T)=0.02 \mathrm{~K}, u(P)=0.001 \mathrm{MPa}, u\left(x_{i}\right)=0.0016$ (water), 0.0007 (ethanol), 0.0015 (isobutanol), and $0.0002(\mathrm{NaCl})$.

## 3. RESULTS AND DISCUSSIONS

3.1. Tie-Line Data Measurement. The experimental LLE data of water (1) + ethanol (2) + isobutanol (3) with the presence of NaCl (4) were measured at temperature ranges from 293.15 to 313.15 K and at atmospheric pressure. To validate the experimental procedures, the LLE tie-line data of the ethanol + water + ethyl acetate + potassium acetate system were obtained at $T=313.15 \mathrm{~K}$ and at $P=0.1 \mathrm{MPa}$. The LLE tie-line data were compared with the results reported by Lin et
al. ${ }^{13}$ under the same conditions. The comparison shows satisfactory agreement, as illustrated in Figure S1 in the Supporting Information, confirming the reliability and precision of the measurements.

By using validated experimental procedures, the experimental LLE tie-line data of water (1) + ethanol (2) + isobutanol (3) with the presence of NaCl were investigated at temperature ranges from 293.15 to 313.15 K at $P=0.1$ MPa. As shown in Tables 3 and 4, the solubility of ethanol decreases with increasing temperature. Increasing the NaCl
concentration up to $10 \mathrm{wt} \%$ induces the salting-out effect on solid-liquid-liquid equilibria (SLLE) phase behavior. Therefore, the tie-line data for $10 \mathrm{wt} \%$ are more limited compared to those for $5 \mathrm{wt} \% \mathrm{NaCl}$. The pseudoternary water (1) + ethanol (2) + isobutanol (3) +NaCl (4) system exhibits type-1 LLE behavior, as only one of the binary pairs (water + isobutanol) is partially miscible. At 313.15 K and at a fixed concentration of $10 \mathrm{wt} \% \mathrm{NaCl}$ in the mixture, the solubility of water in isobutanol is 0.2831 in mole fraction, while the solubility of isobutanol in water is 0.0010 in mole fraction. Thus, the effect of adding NaCl on the mutual solubilities of water and isobutanol can provide useful insight into liquid-phase-splitting enhancement, ${ }^{13}$ as shown in Figure 1. The corresponding LLE

![](https://cdn.mathpix.com/cropped/61a9ddb4-14df-4553-8a5c-b608ccc464f6-04.jpg?height=594&width=635&top_left_y=738&top_left_x=278)
Figure 1. LLE for the system water + ethanol + isobutanol without and with NaCl at 293.15 K and atmospheric pressure expressed as salt-free composition: black (without NaCl ), red ( $5 \mathrm{wt} \% \mathrm{NaCl}$ ), and blue ( $10 \mathrm{wt} \% \mathrm{NaCl}$ ).

data for water (1) + ethanol (2) + isobutanol (3) + $5 \mathrm{wt} \%$ NaCl (4) and water (1) + ethanol (2) + isobutanol (3) +10 $\mathrm{wt} \% \mathrm{NaCl}$ (4) are presented in Figures 2-7.

### 3.2. Distribution Coefficient ( $D$ ) and Separation

Factor (S). The separation behavior between two alcohols (ethanol and isobutanol) in the aqueous system with present NaCl can be evaluated through their distribution coefficients

![](https://cdn.mathpix.com/cropped/61a9ddb4-14df-4553-8a5c-b608ccc464f6-04.jpg?height=577&width=637&top_left_y=1855&top_left_x=278)
Figure 2. LLE for the system water + ethanol + isobutanol $+5 \mathrm{wt} \%$ NaCl at 293.15 K and atmospheric pressure expressed as salt-free composition: black (exp), red (ePC-SAFT), blue (eNRTL), and green (feed compositions).

![](https://cdn.mathpix.com/cropped/61a9ddb4-14df-4553-8a5c-b608ccc464f6-04.jpg?height=592&width=642&top_left_y=202&top_left_x=1190)
Figure 3. LLE for the system water + ethanol + isobutanol + $5 \mathrm{wt} \%$ NaCl at 303.15 K and atmospheric pressure expressed as salt-free composition: black (exp), red (ePC-SAFT), blue (eNRTL), and green (feed compositions).

![](https://cdn.mathpix.com/cropped/61a9ddb4-14df-4553-8a5c-b608ccc464f6-04.jpg?height=584&width=642&top_left_y=992&top_left_x=1190)
Figure 4. LLE for the system water + ethanol + isobutanol + $5 \mathrm{wt} \%$ NaCl at 313.15 K and atmospheric pressure expressed as salt-free composition: black (exp), red (ePC-SAFT), blue (eNRTL), and green (feed compositions).

![](https://cdn.mathpix.com/cropped/61a9ddb4-14df-4553-8a5c-b608ccc464f6-04.jpg?height=577&width=640&top_left_y=1772&top_left_x=1192)
Figure 5. LLE for the system water + ethanol + isobutanol + $10 \mathrm{wt} \%$ NaCl at 293.15 K and atmospheric pressure expressed as salt-free composition: black (exp), red (ePC-SAFT), blue (eNRTL), and green (feed compositions).

![](https://cdn.mathpix.com/cropped/61a9ddb4-14df-4553-8a5c-b608ccc464f6-05.jpg?height=592&width=637&top_left_y=204&top_left_x=278)
Figure 6. LLE for the system water + ethanol + isobutanol + $10 \mathrm{wt} \%$ NaCl at 303.15 K and atmospheric pressure expressed as salt-free composition: black (exp), red (ePC-SAFT), blue (eNRTL), and green (feed composition).

![](https://cdn.mathpix.com/cropped/61a9ddb4-14df-4553-8a5c-b608ccc464f6-05.jpg?height=574&width=637&top_left_y=997&top_left_x=278)
Figure 7. LLE for the system water + ethanol + isobutanol + $10 \mathrm{wt} \%$ NaCl at 313.15 K and atmospheric pressure expressed as salt-free composition: black (exp), red (ePC-SAFT), blue (eNRTL), and green (feed compositions).

(D) and separation factor ( $S$ ). The separation factor is calculated from the ratio of the distribution coefficient of ethanol ( $D_{2}$ ) to water ( $D_{1}$ ), as expressed below

$$
\begin{equation*}
S=\frac{D_{2}}{D_{1}}=\frac{x_{2}^{\mathrm{I}} / x_{2}^{\mathrm{II}}}{x_{1}^{\mathrm{I}} / x_{1}^{\mathrm{II}}} \tag{1}
\end{equation*}
$$

Superscripts I and II denote the organic and aqueous phases, respectively, while $x_{1}$ and $x_{2}$ represent the mole fractions of water and ethanol, respectively. The values of the distribution coefficient and the separation factor for the water + ethanol + isobutanol +NaCl system at each temperature are presented in Tables S1 and S2. As can be seen, both the distribution coefficients and separation factors are significantly greater than one (separation factors varying between 24.47 and 242.70), confirming that alcohols can be effectively separated from an aqueous solution. The experimental results indicate that the mixture containing $10 \mathrm{wt} \% \mathrm{NaCl}$ has a high separation factor ( $S=61.87-242.71$ ), which shows the high ability of isobutanol to extract ethanol from water, as depicted in Figure 8. Increasing the NaCl concentration from 5 to $10 \mathrm{wt} \%$ strengthens the salting-out effect, thereby promoting liquid-

![](https://cdn.mathpix.com/cropped/61a9ddb4-14df-4553-8a5c-b608ccc464f6-05.jpg?height=637&width=844&top_left_y=204&top_left_x=1090)
Figure 8. Separation factor of water over ethanol for the mixture water + ethanol + isobutanol +NaCl ; black ( 293.15 K ), red ( 303.15 K ), blue $(313.15 \mathrm{~K})$, squares (without NaCl$)$, circles $(5 \mathrm{wt} \% \mathrm{NaCl})$, and triangles ( $10 \mathrm{wt} \% \mathrm{NaCl}$ ).

liquid phase splitting. The plots of distribution coefficients and separation factor versus the mole fraction of ethanol in the presence of NaCl at temperatures ranging from 293.15 to 313.15 K are shown in Figure 9.

![](https://cdn.mathpix.com/cropped/61a9ddb4-14df-4553-8a5c-b608ccc464f6-05.jpg?height=664&width=842&top_left_y=1238&top_left_x=1090)
Figure 9. Distribution coefficient of ethanol over water of the mixture water + ethanol + isobutanol +NaCl ; black $(293.15 \mathrm{~K})$, red ( 303.15 $\mathrm{K})$, blue $(313.15 \mathrm{~K})$, squares (without NaCl$)$, circles $(5 \mathrm{wt} \% \mathrm{NaCl})$, and triangles ( $10 \mathrm{wt} \% \mathrm{NaCl}$ ).

It can be observed from Figures 8 and 9 that the distribution coefficients and separation factor increase with temperature from 293.15 to 303.15 K . The results clearly show that the salting-out effect of NaCl is very strong and increases the distribution coefficients and separation factor. In general, the addition of NaCl decreases the solubility of alcoholic compounds in the aqueous phase. ${ }^{10,25}$ This is also concluded from the data in Figures 8 and 9. Upon NaCl addition to the aqueous phase, the strong ion-dipole interactions between $\mathrm{Na}^{+} / \mathrm{Cl}^{-}$ions and water molecules decrease the interaction between water and ethanol. ${ }^{26}$ Consequently, the ethanol molecule is preferentially distributed to the organic phase, leading to higher ethanol distribution coefficients and larger

Table 5. Pure-Component Parameters of All Solvents Used in This Work
| organic solvent | water ${ }^{33}$ | ethanol ${ }^{28}$ | isobutanol ${ }^{34}$ | $\mathrm{Na}^{+30}$ | $\mathrm{Cl}^{-30}$ |
| :--- | :--- | :--- | :--- | :--- | :--- |
| $m_{i}^{\text {seg }}$ | 1.2047 | 2.3827 | 2.8440 |  |  |
| $\sigma_{i} / \AA$ | a | 3.1771 | 3.5561 | 2.8232 | 2.7560 |
| $u_{i} / k_{\mathrm{B}} / \mathrm{K}$ | 353.95 | 198.24 | 257.13 | 230.00 | 170.00 |
| $\varepsilon^{\mathrm{AiBi}} / \mathrm{k}_{\mathrm{B}} / \mathrm{K}$ | 2425.7 | 2653.4 | 2514.3 |  |  |
| $\kappa^{\mathrm{AiBi}}$ | 0.04509 | 0.03238 | 0.00391 |  |  |
| $f_{k}$ | 1.5 | 1.6 | 1.5 |  |  |
| $d_{i}^{\text {Born }} / \AA$ |  |  |  | 3.445 | 4.100 |
| ${ }^{a} \sigma=2.7927+\left(10.11 \mathrm{e}^{-0.01775 T / K}-1.417 \mathrm{e}^{-0.01146 T / K}\right)$. |  |  |  |  |  |


separation factors. These results indicate that the recovery of ethanol from an aqueous solution by isobutanol in the presence of NaCl is highly favorable. The salting-out effect not only increases the extraction efficiency but also improves the selectivity. Higher temperatures are favorable for the thermodynamic efficiency, as well (but not necessarily for the energy requirements of the extraction-based process, which is not part of this work).
3.3. Correlation of LLE Tie-Line Data with ePC-SAFT and eNRTL. In this work, both ePC-SAFT and eNRTL were used to correlate the LLE data of the system water + ethanol +2-butanol in the presence of NaCl . To describe charged species, ePC-SAFT was applied as an extension of the PCSAFT equation of state ( EoS ). In PC-SAFT, the residual Helmholtz energy $a^{\text {res }}$ is expressed as the sum of the hard-chain contribution $a^{\text {hc }}$, the dispersion contribution $a^{\text {disp }}$, and the association contribution $a^{\text {assoc }}$. Compared to the original PCSAFT, ${ }^{27,28}$ ePC-SAFT ${ }^{29-32}$, as well as the most recently developed version of ePC-SAFT, ${ }^{21}$ introduces two additional contributions for charged species: the Debye-Hückel contribution $a^{\mathrm{DH}}$ for ion-ion interactions and the Born contribution $a^{\text {Born }}$ for ion-dipolar interactions. Therefore, the residual Helmholtz energy is given by the sum of these five contributions, as shown in eq 2

$$
\begin{equation*}
a^{\mathrm{res}}=a^{\mathrm{hc}}+a^{\mathrm{disp}}+a^{\mathrm{assoc}}+a^{\mathrm{DH}}+a^{\mathrm{Born}} \tag{2}
\end{equation*}
$$

The original ePC-SAFT parameters include the segment number $m_{i}^{\text {seg }}$, the segment diameter $\sigma_{i}$, and the dispersion energy parameter $u_{i} / k_{\mathrm{B}}$ for nonassociating components, where $k_{\mathrm{B}}$ is the Boltzmann constant. Associating components additionally requires the association-energy parameter $\varepsilon^{\mathrm{AiBi}}$ and the association-volume parameter $\kappa^{\mathrm{AiBi}}$. Furthermore, the binary interaction parameter (BIP) $k_{i j}$ between each component pair can be adjusted to modify the dispersion energy of a pair $i j$ in the mixture. The corresponding combining rule is given by eq 3

$$
\begin{equation*}
u_{i j}=\sqrt{u_{i} u_{j}} \cdot\left(1-k_{i j}\right) \tag{3}
\end{equation*}
$$

Additionally, the dielectric constant values of the solution are required to calculate $a^{\mathrm{DH}}$ and $a^{\text {Born }}$. The dielectric constant of electrolyte solutions depends on salt concentration, solvent composition, temperature, and pressure. In this work, the dielectric constant of solutions containing $\mathrm{Na}^{+}$and $\mathrm{Cl}^{-}$ions was represented as the average of the available experimental data, represented by the following mixing rule.

$$
\begin{equation*}
\varepsilon_{\mathrm{r}, \text { mix }}=\frac{\varepsilon_{\mathrm{r}, \text { solvent, mix }}^{\text {salt-free }}}{1+7.01 \cdot x_{\text {ion }}} \tag{4}
\end{equation*}
$$

To apply ePC-SAFT for correlating the LLE of the water + ethanol + isobutanol system in the presence of NaCl , pure-
component parameters, dielectric constants, and binary interaction parameters are required. The pure-component parameters are listed in Table 5, the relative dielectric constants of all components considered in this work are given in Table 6, and the ePC-SAFT binary interaction

Table 6. Relative Dielectric Constants for All the Components Used in This Work at 298.15 K and Atmospheric Pressure ( $\boldsymbol{P}=0.1 \mathrm{MPa}$ )
| component | dielectric constant | references |
| :--- | :---: | :---: |
| water | 78.09 | 35 |
| ethanol | 24.88 | 36 |
| isobutanol | 17.20 | 37 |
| NaCl (Ions) | 8 | 38 |


parameters between component pairs are presented in Table 7. The table also contains a parameter that is used in the most recent ePC-SAFT version, called $f_{k}$, cf. ePC-SAFT ref 21 .

Table 7. ePC-SAFT Binary Interaction Parameters $k_{i j}$ between Different Components Used in This Work ${ }^{a, b}$
| pair | $k_{i j}$ | references |
| :--- | :--- | :--- |
| ethanol/water | -0.06167 | 39 |
| ethanol/isobutanol | 0.00154 | 19 |
| ethanol/ $\mathrm{Na}^{+}$ | 0.05 | 31 |
| ethanol/ $\mathrm{Cl}^{-}$ | 0.8 | 31 |
| water/isobutanol | -0.0071 | 34 |
| water $/ \mathrm{Na}^{+}$ | -0.3 | 31 |
| water/ $\mathrm{Cl}^{-}$ | -0.3 | 31 |
| isobutanol $/ \mathrm{Na}^{+}$ | 0 |  |
| isobutanol $/ \mathrm{Cl}^{-}$ | 0 |  |
| $\mathrm{Na}^{+} / \mathrm{Cl}^{-}$ | 0.317 | 30 |


${ }^{a}$ All binary interaction parameters $k_{i j}$ for the different components were taken from the literature. ${ }^{b}$ Only valid with the ePC-SAFT version from ref 21 and with the pure-component parameters from Table 5.

The eNRTL model has also been used to calculate the activity coefficients of the various species in the solutions, as given by

$$
\begin{equation*}
\frac{G^{E}}{R T}=\frac{G^{E, \mathrm{PDH}}}{R T}+\frac{G^{E, \mathrm{Born}}}{R T}+\frac{G^{E, \mathrm{lc}}}{R T} \tag{5}
\end{equation*}
$$

The long-range interaction contribution is represented by the Pitzer-Debye-Hückel's formula, which is normalized to mole fractions of unity for solvents and zero for electrolytes

$$
\begin{align*}
& \frac{G^{E, \mathrm{PDH}}}{R T}=-\left(\sum_{k} x_{k}\right)\left(\frac{1000}{M_{\mathrm{B}}}\right)^{0.5}\left(\frac{4 A_{\varphi} I_{x}}{\rho}\right) \ln \left(1+\rho I_{x}^{0.5}\right)  \tag{6}\\
& \frac{G^{E, \text { Born }}}{R T}=\frac{Q_{\mathrm{e}}^{2}}{R T}\left(\frac{1}{\varepsilon}-\frac{1}{\varepsilon_{\mathrm{w}}}\right)\left(\frac{\sum_{i} x_{i} z_{i}^{2}}{r_{i}}\right) 10^{-2}  \tag{7}\\
& A_{\varphi}=\frac{1}{3}\left(\frac{2 \pi N_{\mathrm{A}} d}{1000}\right)^{0.5}\left(\frac{Q_{\mathrm{e}}^{2}}{\varepsilon_{\mathrm{w}}} k T\right)^{0.5} \tag{8}
\end{align*}
$$

where $x_{k}$ is the mole fraction of component $k, M_{\mathrm{B}}$ is the molecular weight of the solvent $\mathrm{B}, A_{\varphi}$ is Debye-Hückel's parameter, $I_{x}$ is the ionic strength in mole fraction, $Q_{\mathrm{e}}$ is the electron charge, $\varepsilon_{\mathrm{w}}$ is the dielectric constant of water, $d$ is the density of solvent, $k$ is the Boltzmann constant, $r_{i}$ is Born's radius, and $T$ is the temperature.

$$
\begin{align*}
\ln \gamma_{i}= & \ln \gamma_{i}^{\mathrm{PDH}}+\ln \gamma_{i}^{\mathrm{Born}}+\ln \gamma_{i}^{\mathrm{lc}, m}  \tag{9}\\
\ln \gamma_{i}^{\mathrm{PDH}}= & -\left(\frac{1000}{M_{\mathrm{B}}}\right)^{0.5} A_{\varphi}\left[\left(\frac{2 z_{i}^{2}}{\rho}\right) \ln \left(1+\rho I_{x}^{0.5}\right)\right. \\
& \left.+\frac{z_{i}^{2} I_{x}^{0.5}-2 I_{x}^{0.5}}{1+\rho I_{x}^{0.5}}\right]  \tag{10}\\
\ln \gamma_{i}^{\mathrm{Born}}= & \frac{Q_{\mathrm{e}}^{2}}{2 k T}\left(\frac{1}{\varepsilon}-\frac{1}{\varepsilon_{\mathrm{w}}}\right)\left(\frac{\sum_{i} x_{i} z_{i}^{2}}{r_{i}}\right) 10^{-2}  \tag{11}\\
\ln \gamma_{i}^{\mathrm{lc}, m}= & 2\left(\frac{y^{2} \pm G_{21} \tau_{21}}{\left(y \pm y_{2} G_{21}\right)^{2}}-\frac{y \pm y_{2} G_{12} \tau_{12}}{\left(2 y \pm y_{2} G_{12}\right)^{2}}\right. \\
& \left.+\frac{y \pm G_{12} \tau_{12}}{\left(2 y \pm y_{2} G_{12}\right)}\right) \tag{12}
\end{align*}
$$

Table 8 shows the optimized binary interaction parameters (BIPs) of eNRTL of water-ethanol, water-isobutanol, and

Table 8. Binary Interaction Parameters (BIPs) ${ }^{a}$ of eNRTL for the System Water (1) + Ethanol (2) + Isobutanol (3) + $\mathbf{N a C l}$ (4) at $293.15-313.15 \mathrm{~K}$ and Atmospheric Pressure ( $P =0.1 \mathrm{MPa}$ )
| component (ij) | $b_{\mathrm{ij}}(\mathrm{K})$ | $b_{\mathrm{ji}}(\mathrm{K})$ | $\alpha$ |
| :--- | :--- | :--- | :--- |
| 1-2 | -404.067 | 1823.31 | 0.2 |
| 1-3 | 1271.07 | -150.629 | 0.2 |
| 1-4 | -35,373.9 | -1471.08 | 0.2 |
| 2-3 | -658.804 | 92.8789 | 0.2 |
| 2-4 | 68,085.7 | 36,654.7 | 0.2 |
| 3-4 | 2345.56 | 44,854 | 0.2 |


${ }^{a}$ BIPs by experimental data regression.
water -NaCl obtained in this work. Figures $2-7$ correspond to the experimental LLE tie-line measurements with the correlated results from the ePC-SAFT and eNRTL models for the water + ethanol + isobutanol system at NaCl concentrations of 5 and $10 \mathrm{wt} \%$ over the temperature range of $293.15-313.15 \mathrm{~K}$ at atmospheric pressure. The average absolute deviations (AADs) for the ePC-SAFT and eNRTL models are listed in Tables 9 and 10. From Figures 2-7 and the AAD values, it can be observed that the agreement
between the experimental and the correlated results is satisfactory. However, the e-PC-SAFT model represents the experimental LLE tie-line measurements of the water + ethanol + isobutanol system more accurately at $5 \mathrm{wt} \% \mathrm{NaCl}$ than at 10 $\mathrm{wt} \% \mathrm{NaCl}$. For comparison purposes, the optimized BIPs were also compared with those from the Aspen Plus databank V. 14 (APV140 VLE-LIT) for water-ethanol, water-isobutanol, and ethanol-isobutanol binary systems. The comparison results are presented in Tables S3 and S4 and Figures S2 and S3. The regression of all parameters from the experimental data yields smaller AADs than available BIPs. A relatively high deviation was found for water + ethanol + isobutanol $+10 \mathrm{wt} \% \mathrm{NaCl}$.

Figure $5-7$ show the calculated values in good agreement with the experimental data compared to those in Figure S3.
3.4. Quantum Chemical Calculations. The surface charge distribution of the molecule is represented by the sigma profile, which can be employed to assess intermolecular interactions. Hydrogen bonding energies of NaCl , water, ethanol, and isobutyl alcohol, as well as electrostatic energy differences, are illustrated in Figure 10. The B3LYP functional and the $6-31 \mathrm{G}^{*}$ basis set were employed to optimize the molecular geometry using DFT theory. Subsequently, the molecular structure file was generated using the optimized structure. The $\sigma$-profile of ethanol exhibited four distinct peaks, with the majority of these peaks located in the nonpolar region ( $\sigma=-0.004,0.002 \mathrm{e} \cdot \AA^{-2}$ ), hydrogen bond acceptor region ( $\sigma=0.017 \mathrm{e} \cdot \AA^{-2}$ ), and hydrogen bond donor ( $\sigma= \left.-0.016 \mathrm{e} \cdot \AA^{-2}\right)$. The miscibility of water and ethanol has a highly strong interaction in these three regions. The ionic interaction from NaCl induces strong interactions in the hydrogen bond donor and the hydrogen bond acceptor regions. As a result of the polar nature of water, NaCl is more soluble in water, which facilitates ionization. Isobutanol is significantly less polar, which fails to solvate ions, compared to water.

## 4. CONCLUSIONS

LLE data for the pseudoternary system of water (1) + ethanol (2) + isobutanol (3) in the presence of NaCl (4) with two different concentrations, 5 and $10 \mathrm{wt} \%$, respectively, were measured at temperature ranges from 293.15 to 313.15 K at atmospheric pressure. The addition of NaCl to the mixtures increased the heterogeneous region. Increasing the NaCl concentration from 5 to $10 \mathrm{wt} \%$ intensified the salting-out effect and promoted liquid-liquid phase separation. For example, compared with the salt-free system, the presence of $10 \mathrm{wt} \% \mathrm{NaCl}$ increased the separation factor by up to 25 times and the distribution coefficient by up to 14 times. These findings demonstrate that the salting-out effect improves not only the extraction efficiency but also selectivity. The measured LLE data were further used to validate the applicability of ePCSAFT and eNRTL for modeling the investigated systems. The sigma profiles indicated the molecular interaction of water, ethanol, isobutyl alcohol, and NaCl . Based on the sigma profiles, a strong interaction of water and NaCl could be concluded, which induces the separation factor of the considered system water + ethanol + isobutanol, with improvements in phase-splitting of the liquid-liquid mixture. Regression of all BIPs using the eNRTL model requires 13 parameters, while the ePC-SAFT model predicts the electro-lyte-aqueous ethanol + isobutanol solution with 7 parameters. Using fewer BIPs, ePC-SAFT predicted electrolyte-aqueous water + isobutyl alcohol phase splitting. With available

Table 9. AAD of LLE Tie-Line Data Correlation with ePC-SAFT and eNRTL Models for the Quaternary System of Water (1) + Ethanol (2) + Isobutanol (3) + $\mathbf{5 ~ w t} \boldsymbol{\% ~ N a C l}$ (4) at Elevated Temperatures 293.15-313.15 K and Atmospheric Pressure ( $P=$ 0.1 MPa)
| T/K | model | phase | $\Delta x_{1}{ }^{a}$ | $\Delta x_{2}{ }^{a}$ | $\Delta x_{3}{ }^{a}$ | $\Delta x_{4}{ }^{a}$ | grand $\mathrm{AAD}{ }^{b}$ |
| :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- |
|  |  |  | AAD | AAD | AAD | AAD |  |
| 293.15 | ePC-SAFT | organic | 0.0287 | 0.0031 | 0.0313 | 0.0009 | 0.0161 |
|  |  | aqueous | 0.0482 | 0.0064 | 0.0021 | 0.0083 |  |
|  | eNRTL | organic | 0.0059 | 0.0048 | 0.0010 | 0.0011 | 0.0046 |
|  |  | aqueous | 0.0102 | 0.0025 | 0.0080 | 0.0029 |  |
| 303.15 | ePC-SAFT | organic | 0.0305 | 0.0027 | 0.0336 | 0.0005 | 0.0168 |
|  |  | aqueous | 0.0495 | 0.0061 | 0.0026 | 0.0088 |  |
|  | eNRTL | organic | 0.0053 | 0.0040 | 0.0015 | 0.0007 | 0.0047 |
|  |  | aqueous | 0.0120 | 0.0029 | 0.0097 | 0.0012 |  |
| 313.15 | ePC-SAFT | organic | 0.0245 | 0.0033 | 0.0276 | 0.0006 | 0.0203 |
|  |  | aqueous | 0.0487 | 0.0052 | 0.0025 | 0.0081 |  |
|  | eNRTL | organic | 0.0147 | 0.0117 | 0.0041 | 0.0011 | 0.0078 |
|  |  | aqueous | 0.0128 | 0.0030 | 0.0126 | 0.0027 |  |


${ }^{a} \Delta x_{i}=\sum_{k=1}^{n_{\mathrm{TL}}}\left|x_{i k}^{\text {calc }}-x_{i k}^{\text {expt }}\right| / n_{\mathrm{TL}}{ }^{b}$ grand AAD $=\left[\sum_{k=1}^{n_{\mathrm{TL}}} \sum_{j=1}^{2} \sum_{i=1}^{3}\left|\left(x_{i j k}^{\text {calc }}-x_{i j k}^{\text {expt }}\right)\right|\right] / 8 n_{\mathrm{TL}}$, where $n_{\mathrm{TL}}$ is the number of tie-lines.

Table 10. AAD of LLE Tie-Line Data Correlation with ePC-SAFT and eNRTL Models for the Quaternary System of Water (1) + Ethanol (2) + Isobutanol (3) + $\mathbf{1 0 ~ w t ~ \% ~ N a C l}$ (4) at Elevated Temperatures 293.15-313.15 K and Atmospheric Pressure ( $P =0.1 \mathrm{MPa}$ )
| T/K | model | phase | $\Delta x_{1}{ }^{a}$ | $\Delta x_{2}{ }^{a}$ | $\Delta x_{3}{ }^{a}$ | $\Delta x_{4}{ }^{a}$ | grand $\mathrm{AAD}{ }^{b}$ |
| :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- |
|  |  |  | AAD | AAD | AAD | AAD |  |
| 293.15 | ePC-SAFT | organic | 0.0729 | 0.0096 | 0.0739 | 0.0005 | 0.0342 |
|  |  | aqueous | 0.0975 | 0.0020 | 0.0007 | 0.0161 |  |
|  | eNRTL | organic | 0.0010 | 0.0026 | 0.0022 | 0.0007 | 0.0086 |
|  |  | aqueous | 0.0312 | 0.0008 | 0.0062 | 0.0241 |  |
| 303.15 | ePC-SAFT | organic | 0.0760 | 0.0021 | 0.0774 | 0.0005 | 0.0341 |
|  |  | aqueous | 0.0973 | 0.0023 | 0.0008 | 0.0165 |  |
|  | eNRTL | organic | 0.0011 | 0.0022 | 0.0033 | 0.0022 | 0.0083 |
|  |  | aqueous | 0.0289 | 0.0006 | 0.0079 | 0.0204 |  |
| 313.15 | ePC-SAFT | organic | 0.0801 | 0.0039 | 0.0783 | 0.0002 | 0.0357 |
|  |  | aqueous | 0.1021 | 0.0025 | 0.0003 | 0.0179 |  |
|  | eNRTL | organic | 0.0011 | 0.0035 | 0.0028 | 0.0004 | 0.0091 |
|  |  | aqueous | 0.0323 | 0.0020 | 0.0098 | 0.0205 |  |


${ }^{a} \Delta x_{i}=\sum_{k=1}^{n_{\mathrm{TL}}}\left|x_{i k}^{\text {calc }}-x_{i k}^{\text {expt }}\right| / n_{\mathrm{TL}}{ }^{b}$ grand AAD $=\left[\sum_{k=1}^{n_{\mathrm{TL}}} \sum_{j=1}^{2} \sum_{i=1}^{3}\left|\left(x_{i j k}^{\text {calc }}-x_{i j k}^{\text {expt }}\right)\right|\right] / 8 n_{\mathrm{TL}}$, where $n_{\mathrm{TL}}$ is the number of tie-lines.

![](https://cdn.mathpix.com/cropped/61a9ddb4-14df-4553-8a5c-b608ccc464f6-08.jpg?height=662&width=840&top_left_y=1772&top_left_x=175)
Figure 10. Sigma profile ( $\sigma$ ) with the segmented surface for water, ethanol, isobutanol, and NaCl .

experimental data of electrolyte-aqueous water + isobutanol, the ePC-SAFT successfully predicted the phase splitting, which improved the accuracy of the ePC-SAFT for electrolyte predictions. The LLE data obtained in this work can be applied in chemical engineering processes, such as ethanol extraction from fermentation broths or ethanol regeneration from industrial process solutions using salting-out effects. This information provides guidance for the process design and optimization.

## - ASSOCIATED CONTENT

## (5) Supporting Information

The Supporting Information is available free of charge at https://pubs.acs.org/doi/10.1021/acs.jced.5c00780.

> Distribution and separation factor of water + ethanol + isobutanol +NaCl , AAD of LLE tie-line data correlation for water + ethanol + isobutanol +NaCl , and plots of comparison experimental with literature for LLE of ethanol + water + ethyl acetate + potassium acetate and LLE of water + ethanol + isobutanol $+\mathrm{NaCl}(\mathrm{PDF})$

## - AUTHOR INFORMATION

## Corresponding Authors

Christoph Held - Laboratory of Thermodynamics, Department of Biochemical and Chemical Engineering, Technische Universität Dortmund, 44225 Dortmund, Germany; © orcid.org/0000-0003-1074-177X; Email: christoph.held@tu-dortmund.de
Ardila Hayu Tiwikrama - Chemical Engineering and Biotechnology Department, National Taipei University of Technology, Taipei 10608, Taiwan; © orcid.org/0000-0002-0393-6774; Email: ardilahayu@mail.ntut.edu.tw

## Authors

Salal Hasan Khudaida - Laboratory of Thermodynamics, Department of Biochemical and Chemical Engineering, Technische Universität Dortmund, 44225 Dortmund, Germany; Chemical Engineering and Biotechnology Department, National Taipei University of Technology, Taipei 10608, Taiwan
Serli Dwi Rahayu - Chemical Engineering and Biotechnology Department, National Taipei University of Technology, Taipei 10608, Taiwan
Paul Figiel - Laboratory of Thermodynamics, Department of Biochemical and Chemical Engineering, Technische Universität Dortmund, 44225 Dortmund, Germany; © orcid.org/0000-0002-5947-7247
Complete contact information is available at:
https://pubs.acs.org/10.1021/acs.jced.5c00780

## Author Contributions

S.H.K.: conceptualization, formal analysis, methodology, writing-original draft, writing-review and editing, validation, software, and supervision. S.D.R.: data curation, formal analysis, investigation, and software. P.F.: formal analysis and software. C.H. and A.H.T.: conceptualization, funding acquisition, methodology, validation, project administration, writing-review and editing, and supervision.

## Notes

The authors declare no competing financial interest.

## - ACKNOWLEDGMENTS

The authors are grateful for the financial support and the fellowship provided by the National of Science and Technology Council, Taiwan, through Grant No. 113-2221-E-027-007. The authors acknowledge funding from the Alexander von Humboldt Foundation (S.H.K.).

## - REFERENCES

(1) Wannachod, T.; Hronec, M.; Soták, T.; Fulajtárová, K.; Pancharoen, U.; Arpornwichanop, A. Effects of Salt on the LLE and Tie-Line Data for Furfuryl Alcohol - n-Butanol-Water at $\mathrm{T}=$ 298.15 K. J. Mol. Liq. 2016, 218, 50-58.
(2) Simić, Z. V.; Radović, I. R.; Kijevčanin, M. L. Measurement and Correlation of Liquid-Liquid Equilibrium Data for Ternary Systems Water + Acetone + Organic Solvents (Isoamyl Acetate, or Ethyl Butyrate, or 1-Octanol, or 1-Decanol) at 298.15 K and Atmospheric Pressure. J. Chem. Eng. Data 2025, 70, 3295.
(3) Kosolapova, S. M.; Efimov, I.; Grai, K. M.; Pyagay, I. N.; Rudko, V. A. Liquid-Liquid Equilibrium of FAEs - Glycerol Phase - Ethanol in the Ethanol Regeneration and FAE Separation Stage. Colloids Surf., A 2024, 703, No. 135412.
(4) Ma, B.; Tan, J.; Liu, S.; Sun, Z.; Xie, S. Enhancing Biofuel Dewatering Efficiency by Ionizing Food-Additives-Based Salting-out Agents. Sep. Purif. Technol. 2025, 354, No. 128906.
(5) Adeh, R. A.; Nemati-Kande, E.; Hashemi, V. H. Liquid-Liquid Equilibrium and Partitioning of Ciprofloxacin Drug in ATPSs Involve PVP10000 + $\mathrm{Li}_{2} \mathrm{SO} 4$ + Water: The Effect of Temperature and PH. J. Chem. Eng. Data 2024, 69 (8), 2793-2804.
(6) Nann, A.; Held, C.; Sadowski, G. Liquid-Liquid Equilibria of 1Butanol/Water/IL Systems. Ind. Eng. Chem. Res. 2013, 52 (51), 18472-18481.
(7) Khudaida, S. H.; Lee, M. J.; Tiwikrama, A. H. N-[Tris-(Hydroxymethyl)Methyl]-3-Aminopropanesulfonic Acid (TAPS) Effect on the Separation of Carboxylic Acids (Acetic Acid, Propionic Acid, or Lactic Acid) from Their Aqueous Solutions Using Methyl Isobutyl Ketone as Extraction Solvent. J. Mol. Liq. 2023, 383, No. 122151.
(8) Khudaida, S. H.; Liu, J.; Ahmed, M. T.; Chang, Y. C.; Tiwikrama, A. H. Toward Evaluation of Diethylene Glycol-Based Deep Eutectic Solvents for the Separation of Benzene and n -Hexane Using Pseudoternary Liquid-Liquid Equilibrium Data. J. Chem. Eng. Data 2025, 70, 3758.
(9) Mohammad, S.; Held, C.; Altuntepe, E.; Köse, T.; Gerlach, T.; Smirnova, I.; Sadowski, G. Salt Influence on MIBK/Water LiquidLiquid Equilibrium: Measuring and Modeling with EPC-SAFT and COSMO-RS. Fluid Phase Equilib. 2016, 416, 83-93.
(10) Xie, S.; Song, W.; Yi, C.; Qiu, X. Salting-out Extraction Systems of Ethanol and Water Induced by High-Solubility Inorganic Electrolytes. J. Ind. Eng. Chem. 2017, 56, 145-150.
(11) Li, F.; Li, Q.; Wu, S.; Tan, Z. Salting-out Extraction of Allicin from Garlic (Allium Sativum L.) Based on Ethanol/Ammonium Sulfate in Laboratory and Pilot Scale. Food Chem. 2017, 217, 91-97.
(12) Pai, M. U.; Rao, K. M. Salt Effect on Liquid-Liquid Equilibria in the Ethyl Acetate-Ethyl Alcohol-Water System. J. Chem. Eng. Data 1966, 11 (3), 353-356.
(13) Lin, H. M.; Yeh, C. E.; Hong, G. B.; Lee, M. J. Enhancement of Liquid Phase Splitting of Water + Ethanol + Ethyl Acetate Mixtures in the Presence of a Hydrophilic Agent or an Electrolyte Substance. Fluid Phase Equilib. 2005, 237 (1-2), 21-30.
(14) Gmehling, J.; Onken, U. Vapor-Liquid Equilibrium Data Collection. Aqueous-Organic Systems, 1977.
(15) Gomis, V.; Ruiz, F.; Asensi, J. C.; Saquete, M. D. Liquid-LiquidSolid Equilibria for the Ternary Systems Butanols + Water + Sodium Chloride or + Potassium Chloride. J. Chem. Eng. Data 1996, 41 (2), 188-191.
(16) Kumagae, Y.; Mishima, K.; Hongo, M.; Kusunoki, M.; Arai, Y. Effect of Calcium Chloride on Vapor-liquid Equilibria of Alcoholalcohol and Alcohol-water Binary Systems. Can. J. Chem. Eng. 1992, 70 (6), 1180-1185.
(17) Zemp, R. J.; Francesconi, A. Z. Salt Effect on Phase Equilibria by a Recirculating Still. J. Chem. Eng. Data 1992, 37 (3), 313-316.
(18) Fritzsche, R. H.; Stockton, D. L. Systems Containing Isobutanol and Tetrachloroethane. J. Ind. Eng. Chem. 1946, 38 (7), 737-740.
(19) Resa, J. M.; González, C.; Goenaga, J. M.; Iglesias, M. Density, Refractive Index, and Speed of Sound at 298.15 K and Vapor-Liquid Equilibria at 101.3 KPa for Binary Mixtures of Ethyl Acetate + 1Pentanol and Ethanol + 2-Methyl-1-Propanol. J. Chem. Eng. Data 2004, 49 (4), 804-808.
(20) Yang, F.; Ngo, T. D.; Kontogeorgis, G. M.; De Hemptinne, J. C. A Benchmark Database for Mixed-Solvent Electrolyte Solutions: Consistency Analysis Using E-NRTL. Ind. Eng. Chem. Res. 2022, 61 (42), 15576-15593.
(21) Figiel, P.; Yu, G.; Held, C. Predicting Thermodynamic Properties of Ions in Single Solvents and in Mixed Solvents Using a Modified Born Term within the EPC-SAFT Framework. Ind. Eng. Chem. Res. 2025, 64 (18), 9406-9418.
(22) Chen, C. C.; Britt, H. I.; Boston, J. F.; Evans, L. B. Local Composition Model for Excess Gibbs Energy of Electrolyte Systems. Part I: Single Solvent, Single Completely Dissociated Electrolyte Systems. AIChE J. 1982, 28 (4), 588-596.
(23) Tetrisyanda, R.; Khudaida, S. H.; Lee, H. Y.; Lee, M. J. LiquidLiquid Equilibria of Ternary Aqueous Mixtures Containing Alcohols and Alkyl Levulinates. J. Chem. Thermodyn. 2021, 156, No. 106384.
(24) Tiwikrama, A. H.; Altway, S. Effect of EPPS Buffer on the Liquid-Liquid Equilibria of Carboxylic Acids (Acetic Acid or Propionic Acid)-Water-Methyl Isobutyl Ketone at Elevated Temperatures. J. Chem. Eng. Data 2022, 67 (5), 1228-1236.
(25) Malinowski, J. J.; Daugulis, A. J. Salt Effects in Extraction of Ethanol, 1-butanol and Acetone from Aqueous Solutions. AIChE J. 1994, 40 (9), 1459-1465.
(26) Hasseine, A.; Kabouche, A.; Meniaic, A. H.; Korichi, M. Salting Effect of NaCl and KCl on the Liquid-Liquid Equilibria of Water + Ethyl Acetate + Ethanol System and Interaction Parameters Estimation Using the Genetic Algorithm. Desalin. Water Treat. 2011, 29 (1-3), 47-55.
(27) Gross, J.; Sadowski, G. Perturbed-Chain SAFT: An Equation of State Based on a Perturbation Theory for Chain Molecules. Ind. Eng. Chem. Res. 2001, 40 (4), 1244-1260.
(28) Gross, J.; Sadowski, G. Application of the Perturbed-Chain SAFT Equation of State to Associating Systems. Ind. Eng. Chem. Res. 2002, 41 (22), 5510-5515.
(29) Cameretti, L. F.; Sadowski, G.; Mollerup, J. M. Modeling of Aqueous Electrolyte Solutions with Perturbed-Chain Statistical Associated Fluid Theory. Ind. Eng. Chem. Res. 2005, 44 (9), 33553362.
(30) Held, C.; Reschke, T.; Mohammad, S.; Luza, A.; Sadowski, G. EPC-SAFT Revised. Chem. Eng. Res. Des. 2014, 92 (12), 2884-2897.
(31) Bülow, M.; Ascani, M.; Held, C. EPC-SAFT Advanced - Part I: Physical Meaning of Including a Concentration-Dependent Dielectric Constant in the Born Term and in the Debye-Hückel Theory. Fluid Phase Equilib. 2021, 535, No. 112967.
(32) Ascani, M.; Held, C. Prediction of Salting-out in Liquid-Liquid Two-Phase Systems with EPC-SAFT: Effect of the Born Term and of a Concentration-Dependent Dielectric Constant. Z. Anorg. Allg. Chem. 2021, 647 (12), 1305-1314.
(33) Cameretti, L. F.; Sadowski, G. Modeling of Aqueous Amino Acid and Polypeptide Solutions with PC-SAFT. Chem. Eng. Process.: Process Intensif. 2008, 47 (6), 1018-1025.
(34) Wu, Y. Y.; Pan, D. T.; Zhu, J. W.; Chen, K.; Wu, B.; Ji, L. J. Liquid-Liquid Equilibria of Water + 2,3-Butanediol + Isobutanol at Several Temperatures. Fluid Phase Equilib. 2012, 325, 100-104.
(35) Floriano, W. B.; Nascimento, M. A. C. Dielectric Constant and Density of Water as a Function of Pressure at Constant Temperature. Braz. J. Phys. 2004, 34 (1), 38-41.
(36) Shirke, R. M.; Chaudhari, A.; More, N. M.; Patil, P. B. Temperature Dependent Dielectric Relaxation Study of Ethyl Acetate - Alcohol Mixtures Using Time Domain Technique. J. Mol. Liq. 2001, 94 (1), 27-36.
(37) Sastrya, S.; Parvateesama, K.; Vishwamb, T.; Murthyc, V. R. K. Investigation of Intermolecular Interaction between Isobutanol and Methyl Benzoate Using Excess Dielectric and Thermodynamic Parameters. Int. J. Eng. Res. 2014, 3 (3), 1779-1787.
(38) Andeen, C.; Fontanella, J.; Schuele, D. Low-Frequency Dielectric Constant of $\mathrm{LiF}, \mathrm{NaF}, \mathrm{NaCl}, \mathrm{NaBr}, \mathrm{KCl}$, and KBr by the Method of Substitution. Phys. Rev. B 1970, 2 (12), No. 5068.
(39) Dohrn, S.; Reimer, P.; Luebbert, C.; Lehmkemper, K.; Kyeremateng, S. O.; Degenhardt, M.; Sadowski, G. Thermodynamic Modeling of Solvent-Impact on Phase Separation in Amorphous Solid Dispersions during Drying. Mol. Pharmaceutics 2020, 17 (7), 27212733.
![](https://cdn.mathpix.com/cropped/61a9ddb4-14df-4553-8a5c-b608ccc464f6-10.jpg?height=1136&width=841&top_left_y=1496&top_left_x=1095)


[^0]:    Received: November 24, 2025
    Revised: December 20, 2025
    Accepted: December 23, 2025
    Published: December 30, 2025

[^1]:    ${ }^{a}$ From supplier. ${ }^{b}$ Milli-Q purification system with a resistivity of $18.3 \mathrm{M} \Omega \cdot \mathrm{cm}$.

