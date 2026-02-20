# Calculation of Multiphase Equilibria Containing Mixed Solvents and Mixed Electrolytes: General Formulation and Case Studies 

Moreno Ascani, Gabriele Sadowski, and Christoph Held*

Cite This: J. Chem. Eng. Data 2022, 67, 1972-1984
Read Online
Downloaded via BRIGHAM YOUNG UNIV on February 13, 2026 at 18:36:22 (UTC). See https://pubs.acs.org/sharingguidelines for options on how to legitimately share published articles.


#### Abstract

Aqueous and non-aqueous electrolyte solutions are ubiquitous in chemical and biochemical applications, especially in innovative processes, and they play a major role in geochemistry, environmental science, and numerous other scientific fields. Despite the obvious importance of electrolyte systems, research success on electrolyte thermodynamics is still behind all of the advances on non-electrolyte thermodynamics. After decades of research, several issues of thermodynamic models for electrolytes remain the object of discussion in the thermodynamic community. Still today, only a few simulation packages offer a general approach to calculate phase equilibria of electrolyte systems in a broad application range regarding the kind of salts and solvent, the number of phases, and the number of non-ionic or ionic species involved. In this work, the general background and the ![](https://cdn.mathpix.com/cropped/790aed8a-6454-4c32-b8d0-0b109e9556ff-01.jpg?height=364&width=507&top_left_y=821&top_left_x=1425) assumptions behind the equilibrium conditions of multiphase electrolyte systems with distributed ions are reviewed. A general methodology is proposed, which can be used to determine the number of liquid phases and their composition at equilibrium of any electrolyte system independent of the number of components. The algorithm was implemented in a FORTRAN routine using the equation of state "ePC-SAFT advanced" (Bülow, M.; Ascani, M.; Held, C. ePC-SAFT advanced-Part I: Physical meaning of including a concentration-dependent dielectric constant in the born term and in the Debye-Hückel theory. Fluid Phase Equilib. 2021, 535, 112967) to estimate the fugacity of each species. The algorithm was successfully tested against experimental data using case studies including three-phase liquid-liquid-liquid equilibria with two ionic species and two-phase liquid-liquid equilibria with three ionic species.


## 1. INTRODUCTION

The importance of electrolytes in industrial and scientific research can hardly be overemphasized. Electrolyte systems are present in geochemistry, environmental science, medicine, and in many systems encountered in the chemical and biochemical industry. ${ }^{1}$ They play a major role in designing downstream processes for the extraction and purification of biomolecules, ${ }^{2,3}$ for sour-gas absorption in petrochemistry and for natural gas purification, ${ }^{4}$ for carbon capture and storage, ${ }^{5}$ and in some innovative separation processes (for instance, those based on ionic liquids (ILs) ${ }^{6}$ ). Modeling and calculation of phase equilibria containing electrolytes is a very active field of research, motivated by the enormous implication of electrolyte solutions in several industrial and scientific applications. ${ }^{1}$ ILs are of special importance, and Joan Brennecke's group ${ }^{7}$ has contributed significantly to the understanding of the dissociation behavior of pure ILs and mixtures containing ILs. Her group found that typical dissociation degrees of $75 \%$ can be found in many pure ILs, and the range goes from $45 \%$ to almost complete dissociation depending on the type of IL. Thus, molecular as well as ionic approaches exist for the modeling of IL solutions. The two most successful theories to describe the long-range (LR) ion-ion electrostatic interactions in electrolyte solutions are the Debye-Hückel theory, ${ }^{8}$ which is based on a solution of the linearized Poisson-Boltzmann differential equations, and the more theoretically founded
mean spherical approximation (MSA), ${ }^{9}$ which is based on a solution of the Ornstein-Zernike equation on the MSA. Both Debye-Hückel and MSA show a similar performance for calculating the Helmholtz energy and its partial derivatives. ${ }^{10}$ Caused by the increasing importance of electrolyte solutions in the chemical industry, classical thermodynamic models have been further developed to deal with complex electrolyte solutions. Almost all of the electrolyte models make use of an already existing expression to describe the short-range (SR) interactions and either the Debye-Hückel or MSA theory to describe the LR interactions. Popular electrolyte models are the Pitzer model, ${ }^{11,12}$ e-NRTL, ${ }^{13,14}$ e-UNIQUAC, ${ }^{15,16}$ and electrolyte models based on an equation of state (EoS) for the SR interactions, e.g., cubic EoS, ${ }^{17,18}$ SAFT-VR, ${ }^{19,20}$ and PCSAFT. ${ }^{21-24}$ Different modeling strategies can be applied: this includes considering association between ions and solvent or between cations and anions, treating the dielectric constant as salt-dependent or considering the Born contribution ${ }^{25}$ (see ref

[^0]26 for an extensive review) to the Helmholtz energy. Some models consider the ions as species in the system, ${ }^{27}$ whereas others treat the salt as an own component ${ }^{21,28}$ (ion-based vs salt-based approach). While both approaches are useful, ${ }^{26}$ only the ion-based approach guarantees flexibility and allows modeling of electrolyte system containing complex mixtures of multiple salts, especially if these contain a common ion. More recent works consider the behavior of electrolyte solutions containing organic solvents, which was reviewed, e.g., by Held. ${ }^{29}$ Besides the development of new thermodynamic models for electrolytes, great effort has been spent to further develop algorithmic approaches for the efficient computation of phase equilibria, with or without coupled chemical equilibria such as reaction equilibria or dissociation equilibria. Gautam et al. ${ }^{30}$ extended the RAND ${ }^{31}$ algorithm to solve the phase equilibria and chemical equilibria of electrolyte systems through global minimization of the Gibbs energy subject to mass balance and electroneutrality constraints. The group of Brennecke ${ }^{32,33}$ proposed an asymmetric framework to model the liquid-liquid equilibrium (LLE) of mixed solvents containing ILs by accounting for different degrees of dissociation in the different phases. Zuend et al. ${ }^{34}$ published a method to calculate the LLE of multicomponent organicwater electrolyte systems using physiochemical constraints to generate initial guesses. Focusing on geochemical systems, Tsanas et al. ${ }^{35}$ applied the first-order Lagrange multipliers and the RAND method to describe multiphase chemical equilibria in electrolyte systems. Besides, there are a lot of publications in the literature dealing with the simulation of chemical and phase equilibrium of multicomponent systems relevant to geochemical applications. ${ }^{36-44}$ However, these works mainly consider one liquid phase, one vapor phase, and one or more solid phases, which coexist at equilibrium while dissolved ions are restricted to be in the liquid phase. Very recently, Tsanas et al. ${ }^{45}$ extended their previously developed methodology: they included the electrochemical potential and the electroneutrality condition in the phase stability analysis and in the system of working equations. They tested their methodology to a 12 -component two-phase liquid-liquid reactive system containing multiple ions, which can distribute between the two liquid phases.

As it can be seen from the last paragraph, there is an increasing interest in phase-equilibrium calculations containing electrolytes; however, only a very limited number of publications exist on generic algorithmic approaches to treat complex electrolytic systems containing several ions that can distribute between different phases. The purpose of the present work is (1) to review and discuss the background and the assumptions behind the equilibrium conditions (thermodynamic framework) for multiphase systems with distributed ions, (2) to present a simple and general phase-equilibrium calculation methodology that includes multiple ions, and (3) to design and test an algorithm for performing multiphase flash calculations in systems containing mixed solvents and mixed electrolytes.

## 2. MATHEMATICAL DERIVATION OF EQUILIBRIUM CONDITIONS

Großmann et al. ${ }^{46}$ and Luckas et al. ${ }^{47}$ presented detailed treatments on the thermodynamics of phase equilibria containing electrolytes. Großmann et al. ${ }^{46}$ have shown that equilibrium conditions of electrolytes in a liquid-liquid twophase system (LLTPS) can be formulated, for a general
number of ionic components, without accounting for the electrostatic potential difference between the two phases. In the following, the mathematical derivation of Großmann et al. and Luckas et al. will be extended to include any number of phases and will be used in the following sections to derive a general approach for computing phase equilibria containing distributed ions. For a charged component $i$ in a phase $j$, the dependence of the Gibbs energy $G^{(j)}$ on mole number $n_{i}^{(j)}$ is described by the electrostatic chemical potential, defined by eq 1. ${ }^{46,47}$

$$
\begin{equation*}
\mu_{i}^{(j) \mathrm{el}}=\left(\frac{\partial G^{(j)}}{\partial n_{i}^{(j)}}\right)_{T, p, n_{j \neq i}}=\mu_{i}^{(j)}+z_{i} F \Phi^{(j)} \tag{1}
\end{equation*}
$$

In eq $1, \mu_{i}^{(j)}$ is the classical chemical potential, $F$ denotes the Faraday constant, and $\Phi^{(j)}$ is the electrostatic potential, which either arises from an external field or from net (positive or negative) charges within phase $j$. To derive a general equilibrium condition for the phase equilibrium of systems containing ionic and molecular components, we start with the formulation of the Gibbs energy of the initially electroneutral system. For a system containing $\pi$ phases and $N$ components comprising neutral species $N^{\text {neut }}$ and charged ions $N^{\text {ch }}$ ( $N= N^{\text {neut }}+N^{\text {ch }}$ ), the Gibbs energy of the system includes a "chemical" contribution ( $\mu_{i}^{(j)}$ ) and an "electrical" ( $z_{i} F \Phi^{(j)}$ ) contribution and is given by eq 2 .

$$
\begin{align*}
G & =G^{(1)}+G^{(2)}+\ldots+G^{(\pi)} \\
& =\sum_{i=1}^{N^{\text {neut }}} \sum_{j=1}^{\pi} n_{i}^{(j)} \mu_{i}^{(j)}+\sum_{i=1}^{N^{\mathrm{ch}}} \sum_{j=1}^{\pi} n_{i}^{(j)}\left(\mu_{i}^{(j)}+z_{i} F \Phi^{(j)}\right) \tag{2}
\end{align*}
$$

In a nonreactive system, all components must fulfill the mass-balance condition; i.e., the overall mole number of each component $N_{i}$ in the system remains constant. Furthermore, at equilibrium, the net charge of each phase ( $\sum_{i} n_{i}^{(j)} z_{i}$ ) attains a constant but unknown value, which depends on the final charge distribution $\Delta_{\mathrm{ch}}^{(j)}$. The necessary and sufficient condition for thermodynamic equilibrium of a generic isobaric isothermal system requires that the Gibbs energy (formulated, for an electrolyte system, by eq 2) is globally minimized. ${ }^{48}$ Thus, a constrained optimization problem can be formulated

$$
\begin{align*}
& \underset{\underline{n}}{\min } G=\sum_{i=1}^{N^{\text {neut }}} \sum_{j=1}^{\pi} n_{i}^{(j)} \mu_{i}^{(j)}+\sum_{i=1}^{N^{\mathrm{ch}}} \sum_{j=1}^{\pi} n_{i}^{(j)}\left(\mu_{i}^{(j)}+z_{i} F \Phi^{(j)}\right)  \tag{3}\\
& \text { s.t. } \quad N_{i}-\sum_{j=1}^{\pi} n_{i}^{(j)}=0 \quad i=1, \ldots, N  \tag{4}\\
& \Delta_{\mathrm{ch}}^{(j)}-\sum_{i=1}^{N^{\mathrm{ch}}} z_{i} n_{i}^{(j)}=0 \quad j=1, \ldots, \pi-1  \tag{5}\\
& n_{i}^{(j)} \geq 0 \quad i=1, \ldots, N \quad j=1, \ldots, \pi \tag{6}
\end{align*}
$$

with $\bar{n}=\left(n_{1}^{(1)}, n_{2}^{(1)}, \ldots, n_{i}^{(j)}, \ldots, n_{N}^{(\pi)}\right)$ as the array containing all mole numbers of all components in all phases, which build up the complete set of $N \cdot \pi$ variables of the problem. The condition given by eq 5 is defined only for $\pi-1$ phases due to the electroneutrality of the whole system. Based on the constrained optimization problem defined by eqs $3-5$, the Lagrangian function $\mathcal{L}$ can be defined by eq 7 .

$$
\begin{align*}
& \mathcal{L}(\bar{n}, \bar{\lambda}, \bar{\delta})=\sum_{i=1}^{N^{\text {neut }}} \sum_{j=1}^{\pi} n_{i}^{(j)} \mu_{i}^{(j)} \\
& \quad+\sum_{i=1}^{N^{\mathrm{ch}}} \sum_{j=1}^{\pi} n_{i}^{(j)}\left(\mu_{i}^{(j)}+z_{i} F \Phi^{(j)}\right)+\sum_{i=1}^{N} \lambda_{i}\left(N_{i}-\sum_{j=1}^{\pi} n_{i}^{(j)}\right) \\
& \quad+\sum_{j=1}^{\pi-1} \delta_{j}\left(\Delta_{\mathrm{ch}}^{(j)}-\sum_{i=1}^{N^{\mathrm{ch}}} z_{i} n_{i}^{(j)}\right) \tag{7}
\end{align*}
$$

The symbols $\bar{\lambda}=\left(\lambda_{1}, \lambda_{2}, \ldots, \lambda_{N}\right), \lambda_{i} \in \mathbb{R}$ and $\bar{\delta}=\left(\delta_{1}, \delta_{2}, \ldots, \delta_{\pi}\right), \delta_{j} \in \mathbb{R}$ represent the arrays of the

Lagrange multipliers, respectively, for the constraints, eq 4 and eq 5 . Necessary conditions of thermodynamic equilibrium are given by a set of positive mole numbers $\bar{n}^{*}$ and of the Lagrange multipliers $\bar{\lambda}^{*}, \bar{\delta}^{*}$ that satisfy eq 8 .

$$
\begin{equation*}
\nabla \mathcal{L}\left(\bar{n}^{*}, \bar{\lambda}^{*}, \bar{\delta}^{*}\right)=\overline{0} \tag{8}
\end{equation*}
$$

The partial derivatives in eq 5 with respect to one of the Lagrange multipliers ( $\lambda_{i}$ or $\delta_{j}$ ) result in one of the constraints given by eq 4 or 5 . The partial derivatives with respect to one of the mole numbers $n_{i}^{(j)}$ lead to either eq 9 or 10 , depending on the nature of component $i$ (if it is a molecule or a charged species).

$$
\begin{array}{lrl}
\frac{\partial \mathcal{L}}{\partial n_{i}^{(j)}}=\mu_{i}^{(j)}-\lambda_{i}^{*}=0 & \text { if } z_{i}=0 & \left(i=1, \ldots, N^{\text {neut }} ; j=1, \ldots, \pi\right) \\
\frac{\partial \mathcal{L}}{\partial n_{i}^{(j)}}=\mu_{i}^{(j)}-\lambda_{i}^{*}+\left(F \Phi^{(j)}-\delta_{j}^{*}\right) z_{i}=0 & \text { if } z_{i} \neq 0 & \left(i=1, \ldots, N^{\mathrm{ch}} ; j=1, \ldots, \pi\right) \tag{10}
\end{array}
$$

The Lagrange multiplier $\lambda_{i}^{*} \in \mathbb{R}$ in eqs 9 and 10 takes the same value for each component, and the Lagrange multiplier $\delta_{j}^{*} \in \mathbb{R}$ takes the same value for each phase. By eliminating the Lagrange multipliers from eq 10, the equality of the chemical potential in each phase for each neutral component is recovered (eq 11).

$$
\begin{equation*}
\mu_{i}^{(1)}=\mu_{i}^{(2)}=\ldots=\mu_{i}^{(\pi)} \tag{11}
\end{equation*}
$$

To derive an equivalent condition for charged species, eq 10 is divided by the ion valence $z_{i}$; then, the term $\left(F \phi^{(j)}-\delta_{j}^{*}\right)$ is grouped into a common real number $\tilde{\delta}_{j}^{*}$, finally obtaining eq 12,

$$
\begin{gather*}
\frac{1}{z_{i}} \mu_{i}^{(j)}-\tilde{\lambda}_{i}^{*}+\tilde{\delta}_{j}^{*}=0 \quad \text { if } z_{i} \neq 0 \\
\quad\left(i=1, \ldots, N^{\mathrm{ch}} ; j=1, \ldots, \pi\right) \tag{12}
\end{gather*}
$$

where $\tilde{\lambda}_{i}^{*}$ denotes the original Lagrange multiplier $\lambda_{i}^{*}$ divided by the ion valence $z_{i}$. Both $\tilde{\delta}_{j}^{*}$ and $\lambda_{i}^{*}$ are real numbers, like the original Lagrange multipliers $\delta_{j}^{*}$ and $\lambda_{i}^{*}$; these must become constant values at equilibrium, depending only on the given phase $j$ (for $\tilde{\delta}_{j}^{*}$ ) or the given component $i$ (for $\tilde{\lambda}_{i}^{*}$ ). In order to eliminate $\delta_{j}^{*}$, any of the equations given by eq 12 referring to any of the charged species $i=1, \ldots, N^{\text {ch }}$ in the phase $j$ can be linearly combined with another equation from eq 12 for another component $k=1, \ldots, N^{\mathrm{ch}}, k \neq i$, in the same phase $j$. The set of equations given by eq 12 can be combined within the same phase $j$ to eliminate $\delta_{j}^{*}$. For this, a set of independent couples of cations and anions can be chosen, and the respective eqs 12 are combined. The procedure is repeated for all other phases $j=1, \ldots, \pi$ for the same couple of independent cations and anions; thus, out of the original $N^{\text {ch }} \cdot \pi$ equilibrium conditions given by eq 12 for all of the phases, only ( $N^{\text {ch }}- 1) \cdot \pi$ equilibrium conditions remain, which are given by eq 13 for each independent pair of cations and anions,

$$
\begin{equation*}
\frac{1}{z_{\mathrm{cat}}} \mu_{\mathrm{cat}}^{(j)}-\frac{1}{z_{\mathrm{an}}} \mu_{\mathrm{an}}^{(j)}=\tilde{\lambda}_{\mathrm{cat}}^{*}-\tilde{\lambda}_{\mathrm{an}}^{*} \tag{13}
\end{equation*}
$$

or in equivalent form by eq 14 .

$$
\begin{equation*}
\frac{1}{\left|z_{\mathrm{cat}}\right|} \mu_{\mathrm{cat}}^{(j)}+\frac{1}{\left|z_{\mathrm{an}}\right|} \mu_{\mathrm{an}}^{(j)}=\tilde{\lambda}_{\mathrm{cat}}^{*}-\tilde{\lambda}_{\mathrm{an}}^{*} \tag{14}
\end{equation*}
$$

The term $\tilde{\lambda}_{\text {cat }}^{*}-\tilde{\lambda}_{\text {an }}^{*}$ is the same for all independent couples of cations and anions upon which the number $N^{\text {ch }}-1$ equations are defined; the term is independent of the considered phase. Thus, a similar condition as given by eq 11 involving the chemical potential of pairs of counterions must be fulfilled at equilibrium, which is given by eq 15 .

$$
\begin{align*}
& \frac{1}{\left|z_{\mathrm{cat}}\right|} \mu_{\mathrm{cat}}^{(1)}+\frac{1}{\left|z_{\mathrm{an}}\right|} \mu_{\mathrm{an}}^{(1)}=\frac{1}{\mid z_{\mathrm{cat}}} \mu_{\mathrm{cat}}^{(2)}+\frac{1}{\left|z_{\mathrm{an}}\right|} \mu_{\mathrm{an}}^{(2)}=\ldots \\
& \quad=\frac{1}{\mid z_{\mathrm{cat}}} \mu_{\mathrm{cat}}^{(\pi)}+\frac{1}{\left|z_{\mathrm{an}}\right|} \mu_{\mathrm{an}}^{(\pi)} \tag{15}
\end{align*}
$$

Equation 15 defines at the same time the so-called mean ionic chemical potential, often denoted by the symbol $\mu_{ \pm}$(eq 16 ). Likewise, also other properties, among which a mean ionic fugacity and a mean ionic activity as well as a mean ionic fugacity coefficient or a mean ionic activity coefficient (MIAC), can be defined (eqs 16-20).

$$
\begin{align*}
& \mu_{ \pm}=\frac{\frac{1}{\mid z_{\text {cat }}} \mu_{\text {cat }}+\frac{1}{\mid z_{\text {an }}} \mu_{\text {an }}}{\frac{1}{\left|z_{\text {cat }}\right|}+\frac{1}{\left|z_{\text {an }}\right|}}  \tag{16}\\
& f_{ \pm}=\left(\left(f_{\text {cat }}\right)^{1 /\left|z_{\text {cat }}\right|} \cdot\left(f_{\text {an }}\right)^{1 /\left|z_{\text {an }}\right|}\right)^{1 /\left[\left(1 /\left|z_{\text {cat }}\right|\right)+\left(1 /\left|z_{\text {an }}\right|\right)\right]}  \tag{17}\\
& a_{ \pm}=\left(\left(a_{\text {cat }}\right)^{1 /\left|z_{\text {cat }}\right|} \cdot\left(a_{\text {an }}\right)^{1 /\left|z_{\text {an }}\right|}\right)^{1 /\left[\left(1 /\left|z_{\text {cat }}\right|\right)+\left(1 /\left|z_{\text {an }}\right|\right)\right]}  \tag{18}\\
& \varphi_{ \pm}=\left(\left(\varphi_{\text {cat }}\right)^{1 /\left|z_{\text {cat }}\right|} \cdot\left(\varphi_{\text {an }}\right)^{1 /\left|z_{\text {an }}\right|}\right)^{1 /\left[\left(1 /\left|z_{\text {cat }}\right|\right)+\left(1 /\left|z_{\text {an }}\right|\right)\right]}  \tag{19}\\
& \gamma_{ \pm}=\left(\left(\gamma_{\text {cat }}\right)^{1 /\left|z_{\text {cat }}\right|} \cdot\left(\gamma_{\text {an }}\right)^{1 /\left|z_{\text {an }}\right|}\right)^{1 /\left[\left(1 /\left|z_{\text {cat }}\right|\right)+\left(1 /\left|z_{\text {an }}\right|\right)\right]} \tag{20}
\end{align*}
$$

Equation 20 is equivalent to the definition of MIAC encountered in the literature, ${ }^{49,50}$ which is based on the stoichiometric coefficients $\nu_{i}$ rather than the charge number $\left|z_{i}\right|$ of each ion in a salt (see eq 21). The salt CA dissolves into $\nu_{\text {cat }}$ moles of cation Cat and $\nu_{\text {an }}$ moles of anion An (according to $\mathrm{CA} \rightarrow \nu_{\mathrm{cat}} \mathrm{Cat}+\nu_{\mathrm{an}} \mathrm{An}$ ), and due to electroneutrality, it holds that $\nu_{\text {cat }}\left|z_{\text {cat }}\right|=\nu_{\text {an }}\left|z_{\text {an }}\right|$; thus, the definition in eq 21 is obtained.

$$
\begin{align*}
\gamma_{ \pm} & =\left(\left(\gamma_{\mathrm{cat}}\right)^{\nu_{\mathrm{cat}}} \cdot\left(\gamma_{\mathrm{an}}\right)^{\nu_{\mathrm{an}}}\right)^{1 /\left(\nu_{\mathrm{cat}}+\nu_{\mathrm{an}}\right)} \\
& =\left(\left(\gamma_{\mathrm{cat}}\right)^{1 / \mid z_{\mathrm{cat}}} \cdot\left(\gamma_{\mathrm{an}}\right)^{1 / \mid \mathrm{zan}_{\mathrm{an}}}\right)^{1 /\left[\left(1 /\left|\mathrm{z}_{\mathrm{cat}}\right|\right)+\left(1 /\left|\mathrm{z}_{\mathrm{an}}\right|\right)\right]} \tag{21}
\end{align*}
$$

Furthermore, for an initially neutral system, it can be assumed that the electroneutrality is always maintained within the bulk of each phase (in other words, the net charge $\Delta_{\mathrm{ch}}^{(j)}$ can be assumed to be zero, and thus, eq 5 becomes effectively $\left(\sum_{i} n_{i}^{(j)} z_{i}=0\right)$ ). This assumption arises from the fact that macroscopic charge separation requires a prohibitive amount of energy, ${ }^{51}$ while microscopic charge separation is possible, but is limited to a small region close to the interfacial region. ${ }^{52}$ Thus, for an initially electroneutral $N$-component system containing molecular as well as charged species, the ( $p, T, \underline{z}$ )flash calculation problem can be solved by computing a system of equations containing the phase-equilibrium conditions (eqs 22 and 23), the mass balance (eq 24), and the electroneutrality constraint of each phase (eq 25).

$$
\begin{align*}
& \mu_{i}^{(1)}=\mu_{i}^{(2)}=\ldots=\mu_{i}^{(\pi)} \quad z_{i}=0  \tag{22}\\
& \mu_{ \pm, i k}^{(1)}=\mu_{ \pm, i k}^{(2)}=\ldots=\mu_{ \pm, i k}^{(\pi)} \quad z_{i}, z_{k} \neq 0 ; i=\text { cat }, k=\text { an }  \tag{23}\\
& N_{i}-\sum_{j=1}^{\pi} n_{i}^{(j)}=0 \quad i=1, \ldots, N  \tag{24}\\
& \sum_{k=1}^{N^{\mathrm{ch}}} z_{i} n_{i}^{(j)}=0 \quad j=1, \ldots, \pi-1 \tag{25}
\end{align*}
$$

The equilibrium conditions build up a set of $(N-1) \cdot(\pi-$ 1) equations. Furthermore, the mass balance and the electroneutrality constraint yield, respectively, $N$ and $\pi-1$ equations. Therefore, a system of $N \cdot \pi$ equations with $N \cdot \pi$ variables (given by the mole numbers of all of the components in all of the phases) must be solved.

## 3. ALGORITHMIC IMPLEMENTATION

### 3.1. Definition of Pairs of Independent Counterions.

Computing the system of eqs $22-25$ in mixtures containing more than two different ions requires as a first step specifying the set of independent counterions needed to build up $N^{\text {ch }}-1$ mean ionic chemical potentials (or, equivalently, mean ionic fugacities) from the $N^{\text {ch }}$ chemical potentials of the single ions. The mean ionic chemical potentials can be represented from the single-ion chemical potentials in a vector-matrix form (eqs 26 and 27).

$$
\begin{align*}
& \overline{M \mu_{ \pm}}=\overline{\bar{E}} \cdot \overline{M \mu_{s}}  \tag{26}\\
& {\left[\begin{array}{c}
\hat{\mu}_{ \pm, 1}^{(j)} \\
\hat{\mu}_{ \pm, 2}^{(j)} \\
\vdots \\
\hat{\mu}_{ \pm, N^{\mathrm{ch}}-1}^{(j)}
\end{array}\right]=\left[\begin{array}{cccc}
e_{1,1} & e_{1,2} & \cdots & e_{1, N^{\mathrm{ch}}} \\
e_{2,1} & e_{2,2} & \cdots & e_{2, N^{\mathrm{ch}}} \\
\vdots & \vdots & \ddots & \vdots \\
e_{N^{\mathrm{ch}}-1,1} & e_{N^{\mathrm{ch}}-1,2} & \cdots & e_{N^{\mathrm{ch}}-1, N^{\mathrm{ch}}}
\end{array}\right] \cdot\left[\begin{array}{c}
\mu_{1}^{(j)} \\
\mu_{2}^{(j)} \\
\vdots \\
\mu_{N^{\mathrm{ch}}}^{(j)}
\end{array}\right]} \tag{27}
\end{align*}
$$

The symbol $\hat{\mu}$ denotes the mean ionic chemical potential in eq 17 with the sum in the denominator $\frac{1}{\left|z_{\text {cat }}\right|}+\frac{1}{\left|z_{\text {an }}\right|}$ set equal to 1 (eq 16). Since the denominator is a constant, its value does
not affect the validity of eq 23 if $\hat{\mu}_{ \pm}$instead of $\mu_{ \pm}$is used. Each row $m$ of the matrix $\overline{\bar{E}}$ has elements $e_{m n}$ which are equal to $1 /\left|z_{n}\right|$, if the component $n$ belongs to the pair of counterions upon which the mean ionic chemical potential is defined, and equal to zero otherwise. The mean ionic chemical potential $\mu_{ \pm, i}$ can be recovered from $\hat{\mu}_{ \pm, i}$ if the value of $\hat{\mu}_{ \pm, i}$ is divided by the sum of the elements $e_{i, j}$ of the row number $i$. The chosen $N^{\mathrm{ch}}$ -1 pair of counterions is an independent set (and, thus, can be used to define the $N^{\text {ch }}-1$ mean ionic chemical potentials) if the matrix $\overline{\bar{E}}$ has full rank, i.e., if $\operatorname{rank}(\overline{\bar{E}})=N^{\text {ch }}-1$. In this section, we propose a simple procedure to build up a matrix $\overline{\bar{E}}$ with $\operatorname{rank}(\overline{\bar{E}})=N^{\text {ch }}-1$. This procedure is a preprocessing step to the actual ( $p, T, \underline{z}$ )-flash calculation, which is operated only once at the beginning of the phase-equilibrium calculation. It consists of the following steps:

1. Initialize the matrix $\overline{\bar{E}}=\overline{\bar{E}}\left(N^{\text {ch }}-1, N^{\text {ch }}\right)$ with all elements $e_{m n}=0$ and determine the number of cations $N_{\text {cat }}$ and anions $N_{\text {an }}$ among all of the components in the feed $F$.
2. Prepare two arrays of integer $v_{\text {cat }}\left(1, \ldots, N_{\text {cat }}\right)$ and $v_{\text {an }}(1, \ldots$, $N_{\mathrm{an}}$ ) representing, respectively, the index of all of the cations and anions in the feed $F$ in the order from highest to lowest cation or anion concentration $z_{i}$. For instance, the value of $v_{\text {cat }}(1)$ represents the index of the cation with the highest mole fraction $z_{v_{\text {cat }}(1)}$, followed by $v_{\text {cat }}(2)$ (if a second cation is present) and so on.
3. If $N_{\mathrm{cat}} \leq N_{\mathrm{an}}$, set the two matrix elements $e_{k, \nu_{\mathrm{cat}}(1)}=1 / \mid z_{v_{\text {cat }}(1)} \mid$ and $e_{k, v_{\text {an }}(k)}=1 /\left|z_{v_{\text {an }}(k)}\right|$ for $k=1, \ldots, N_{\text {an }}$. If $N_{\text {cat }}> N_{\mathrm{an}}$, set the two matrix elements $e_{k, v_{\mathrm{an}}(1)}=1 /\left|z_{v_{\mathrm{an}}(1)}\right|$ and $e_{k, v_{\text {cat }}(k)}=1 /\left|z_{v_{\text {cat }}(k)}\right|$ for $k=1, \ldots, N_{\text {cat }}$.
4. If $N_{\mathrm{cat}} \leq N_{\mathrm{an}}$ and $N_{\mathrm{cat}}>1$, set the two matrix elements $e_{N_{\text {an }}+k, v_{\text {cat }}(k+1)}=1 /\left|z_{v_{\text {cat }}(k+1)}\right|$ and $e_{N_{\text {an }}+k, v_{\text {an }}(k)}=1 /\left|z_{v_{\text {an }}(k)}\right|$ for $k =1, \ldots, N_{\mathrm{cat}}-1$. If $N_{\mathrm{cat}}>N_{\mathrm{an}}$ and $N_{\mathrm{an}}<1$, set the two matrix elements $e_{N_{\text {cat }}+k, v_{\text {an }}(k+1)}=1 /\left|z_{v_{\text {an }}(k+1)}\right|$ and $e_{N_{\text {cat }}+k, v_{\text {cat }}(k)} =1 /\left|z_{v_{\text {cat }}(k)}\right|$ for $k=1, \ldots, N_{\text {an }}-1$.
The preprocessing-step procedure will be applied in the following to a concrete example. Consider a system containing water, 1 -butanol, KCl , and $\mathrm{Na}_{2} \mathrm{SO}_{4}$. The mole fractions of the respective components, given on a salt basis, are $\underline{z}_{\mathrm{F}}=(0.400$, $0.400,0.110,0.090)$. For phase-equilibrium calculation, the composition is translated on a single-ion basis. Considering the components in the order (water, 1 -butanol, $\mathrm{K}^{+}, \mathrm{Cl}^{-}, \mathrm{Na}^{+}$, $\mathrm{SO}_{4}{ }^{2-}$ ), the ion-based feed composition is $\underline{z}_{\mathrm{F}}=(0.310,0.310$, $0.085,0.085,0.140,0.070)$ in mole fraction. The preprocess-ing-step procedure will be explained step by step in the following:
5. A number of anions and cations equal to two are present in the mixture, $N_{\mathrm{cat}}=N_{\mathrm{an}}=2$, and the number of charged species is four, $N^{\text {ch }}=4$. The matrix $\overline{\bar{E}}$ is thus a 3 $\times 4$ matrix and is initialized with all elements equal to zero.

$$
\overline{\bar{E}}=\left[\begin{array}{llll}
0 & 0 & 0 & 0 \\
0 & 0 & 0 & 0 \\
0 & 0 & 0 & 0 \\
0 & 0 & 0 & 0
\end{array}\right]
$$

2. Recalling the order at which the ions (without considering non-charged components) are present in
the feed, $\left(\mathrm{K}^{+}, \mathrm{Cl}^{-}, \mathrm{Na}^{+}, \mathrm{SO}_{4}{ }^{2-}\right)=(1,2,3,4)$, two arrays $v_{\text {cat }}$ and $v_{\text {an }}$ are prepared, respectively, with a length equal to the number of cations $N_{\mathrm{cat}}$ and anions $N_{\mathrm{an}}$ present in the system. Thus, both arrays have a length equal to two. They contain the index of cations and anions in the order of highest to lowest ion mole fraction. The following arrays are recovered:

$$
v_{\mathrm{cat}}=(3,1) \quad v_{\mathrm{an}}=(2,4)
$$

If two cations or anions have the same mole fraction, the order at which they appear in the array is irrelevant and can be set arbitrarily.
3. Since $N_{\mathrm{cat}}=N_{\mathrm{an}}=2$, for $k=1, \ldots, N_{\mathrm{an}}=1,2$, we set the two matrix elements $e_{k, v_{\text {cat }}(1)}=1 /\left|z_{v_{\text {cat }}(1)}\right|$ and $e_{k, v_{\text {an }}(k)}=1 / \mid z_{v_{\mathrm{an}}(k)} \mid$. Thus, for the considered example, $e_{1,3}=e_{2,3}=1$, $e_{1,2}=1$, and $e_{2,4}=1 / 2$. The matrix $\overline{\bar{E}}$ is now given as follows:

$$
\overline{\bar{E}}=\left[\begin{array}{cccc}
0 & 1 & 1 & 0 \\
0 & 0 & 1 & 1 / 2 \\
0 & 0 & 0 & 0
\end{array}\right]
$$

4. In the last step, for $k=1, \ldots, N_{\mathrm{cat}}-1$, we set the two matrix elements $e_{N_{\mathrm{an}}+k, v_{\mathrm{cat}}(k+1)}=1 /\left|z_{v_{\mathrm{cat}}(k+1)}\right|$ and $e_{N_{\mathrm{an}}+k, v_{\mathrm{an}}(k)} =1 /\left|z_{v_{\mathrm{an}}(k)}\right|$. Thus, for the considered example, $e_{3,1}=1$ and $e_{3,2}=1$. The final matrix $\overline{\bar{E}}$, which has $\operatorname{rank}(\overline{\bar{E}})=3$, is given as follows:

$$
\overline{\bar{E}}=\left[\begin{array}{cccc}
0 & 1 & 1 & 0 \\
0 & 0 & 1 & 1 / 2 \\
1 & 1 & 0 & 0
\end{array}\right]
$$

Thus, three pairs of independent cations and anions ( $\left(\mathrm{Cl}^{-}, \mathrm{Na}^{+}\right),\left(\mathrm{Na}^{+}, \mathrm{SO}_{4}{ }^{2-}\right)$, and $\left(\mathrm{K}^{+}, \mathrm{Cl}^{-}\right)$) are defined for the given system using the proposed preprocessingstep methodology.
3.2. Variable Reformulation. The phase-equilibrium problem can be solved by computing the system of $N \cdot \pi$ equations with $N \cdot \pi$ variables (number of moles) defined by eqs $22-25$. However, the full system of equations can be reduced by an appropriate variable reformulation, which allows removing the electroneutrality condition given by eq 25 . The main idea is to use a linear transformation of the number of moles $n_{i}^{(j)}$ of all of the charged species $i=1, \ldots, N^{\text {ch }}$, which forces the phase composition to remain electroneutral by simultaneously changing the mole numbers of a cation and an anion. This requires introducing $N^{\text {ch }}-1$ transformed variables $\xi_{s}^{(j)}$, where $s=1, \ldots, N^{\mathrm{ch}}-1$. Each of the transformed variables can change the number of moles of a pair of an independent cation and anion, already defined in the previous section, by an amount proportional to the inverse of the respective ion valence (eq 28).

$$
\begin{equation*}
\Delta n_{i}^{(j)}=\frac{1}{\left|z_{i}\right|} \xi_{s}^{(j)}\left(z_{i}>0\right) ; \quad \Delta n_{k}^{(j)}=\frac{1}{\left|z_{k}\right|} \xi_{s}^{(j)}\left(z_{k}<0\right) \tag{28}
\end{equation*}
$$

Please note that the transformed variables proposed in eq 28 must not be confused with the reaction coordinates used to balance a system subject to chemical reactions, although both change the composition of the system by changing the composition of multiple components. According to eq 24, the initially electroneutral phase $j$ remains electroneutral after
changing the mole numbers by $\xi_{s}^{(j)}$, since $\Delta n_{i}^{(j)} \cdot z_{i}+\Delta n_{k}^{(j)} \cdot z_{k}$ will be zero. This concept can be generalized to multiple electrolytes and several transformed variables, each of them changing the mole number of an independent pair of counterions according to eq 28, by defining the number of moles of charged species in a vector-matrix form (eqs 29 and 30).

$$
\begin{align*}
& \bar{n}^{(j)}=\bar{n}_{0}^{(j)}+\bar{E}^{T} \cdot \bar{\xi}^{(j)}  \tag{29}\\
& {\left[\begin{array}{c}
n_{1}^{(j)} \\
n_{2}^{(j)} \\
\vdots \\
n_{N^{c h}}^{(j)}
\end{array}\right]=\left[\begin{array}{c}
n_{1,0}^{(j)} \\
n_{2,0}^{(j)} \\
\vdots \\
n_{N^{c h}, 0}^{(j)}
\end{array}\right]+\left[\begin{array}{cccc}
e_{1,1} & e_{2,1} & \cdots & e_{N^{c h}-1,1} \\
e_{1,2} & e_{2,2} & \cdots & e_{N^{c h}-1,2} \\
\vdots & \vdots & \ddots & \vdots \\
e_{1, N^{c h}} & e_{2, N^{c h}} & \cdots & e_{N^{c h}-1, N^{c h}}
\end{array}\right] \cdot\left[\begin{array}{c}
\xi_{1}^{(j)} \\
\xi_{2}^{(j)} \\
\vdots \\
\xi_{N^{c h}-1}^{(j)}
\end{array}\right]} \tag{30}
\end{align*}
$$

The matrix $\overline{\bar{E}}^{T}$ represents the transposed matrix of $\overline{\bar{E}}$ defined in eqs 26 and 27, $\bar{n}^{(j)}$ represents the vector of length $N^{\text {ch }}$ containing the number of moles of charged species in phase $j$, and $\bar{n}_{0}^{(j)}$ represents the vector of length $N^{\mathrm{ch}}$ containing the initial number of moles of the initially electroneutral phase $j$. Thus, by changing the phase amount and phase composition through the $N^{\text {ch }}-1$ variables $\xi_{s}^{(j)}$, the $N^{\text {ch }}$ number of moles of charged species $n_{i}^{(j)}$ are forced to remain in a subspace of $\mathbb{R}^{N^{\mathrm{ch}}}$ with dimension ( $N^{\text {ch }}-1$ ) where the electroneutrality of phase $j$ remains fulfilled. This variable reformulation is similar in spirit to the transformed coordinates proposed by Ung et al. ${ }^{53,54}$ for multireaction systems; the expressions were defined to permanently fulfill the stoichiometry of all of the chemical reactions in the system. The system of equations is thus reformulated, using the transformed variables $\xi_{i}^{(j)}$, by eqs 31 and 32 as the mass balance condition together with the equality of the "classical" and mean ionic chemical potentials (eqs 22 and 23). The electroneutrality condition (eq 25) is implicitly fulfilled.

$$
\begin{align*}
& N_{i}-\sum_{j=1}^{\pi} n_{i}^{(j)}=0 \quad i=1, \ldots, N^{\text {neut }}  \tag{31}\\
& N_{k}-\sum_{j=1}^{\pi} n_{k}^{(j)}=0 \quad k=1, \ldots, N^{\mathrm{ch}}-1 \tag{32}
\end{align*}
$$

Equations 22 and 23 and equations 31 and 32 build a system of $(N-1) \cdot \pi$ equations in $(N-1) \cdot \pi$ variables ( $N^{\text {neut }} \cdot \pi$ number of moles of neutral components $n_{i}^{(j)}$ and $\left(N^{\text {ch }}-1\right) \cdot \pi$ transformed variables $\xi_{s}^{(j)}$ in each phase $j$ ).

### 3.3. Reformulated Tangent-Plane Distance Function.

The tangent-plane distance function (TPDF) for phasestability analysis originally proposed by Michelsen ${ }^{55,56}$ aims at testing whether or not a homogeneous phase at composition $\bar{z}$ can reach a lower Gibbs energy $g$, which is true if an infinitesimal amount of a second phase of composition $\bar{x}$ separates from the original phase. The expression of the TPDF is given by eq 33 .

$$
\begin{equation*}
\operatorname{TPDF}(\bar{x})=\sum_{i=1}^{N} x_{i} \cdot\left(\mu_{i}(\bar{x})-\mu_{i}(\bar{z})\right) \tag{33}
\end{equation*}
$$

The stability criterion proposed by Michelsen states that, if the value of $\operatorname{TPDF}(\bar{x})$ becomes negative at some values of the
trial phase $\bar{x}$, then the homogeneous phase at composition $\bar{z}$ is unstable under the given conditions ( $T$ and $p$ ) and, according to the used model, separates in (at least) one additional phase. The estimate of the formed phase is given by the composition value $\bar{x}^{*}$, which minimizes eq 33 . Equation 33 can be reformulated, in the general case of a multicomponent system with several ionic components, using the transformed variables $\xi_{s}$. Based on an initially electroneutral feed $F$ of composition $\bar{z}$, the estimate of the composition $\bar{x}$ of the second phase can be calculated for a defined matrix $\overline{\bar{E}}$, using number of moles $n_{i}$ of neutral components and transformed variables $\xi_{s}$ (with $n_{i}, \xi_{s}>$ 0 ) as adjusting variables. The array of mole fractions $\bar{x}$ of the second phase results from the variables $\bar{n}, \bar{\xi}$ and the matrix $\overline{\bar{E}}$ using eqs 34 and 35.

$$
\begin{array}{ll}
x_{i}=\frac{n_{i}}{\sum_{k=1}^{N^{\text {neut }}} n_{k}+\sum_{n=1}^{N^{\mathrm{ch}}} \sum_{m=1}^{N^{\mathrm{ch}}-1} e_{n, m} \cdot \xi_{m}} & i=1, \ldots, N^{\text {neut }}  \tag{34}\\
x_{s}=\frac{\sum_{m=1}^{N^{\mathrm{ch}}-1} e_{s, m} \cdot \xi_{m}}{\sum_{k=1}^{N^{\text {neut }}} n_{k}+\sum_{n=1}^{N^{\mathrm{ch}} \sum_{m=1}^{N^{\mathrm{ch}}-1} e_{n, m} \cdot \xi_{m}}} & s=1, \ldots, N^{\mathrm{ch}}
\end{array}
$$

Please note that $\xi_{m}$ has the dimension of moles.

## 4. RESULTS

The thermodynamic framework based on the phase-equilibrium conditions (chapter 2) combined with the proposed variable reformulation (chapter 3) was used for the design of an algorithm to perform phase-equilibrium calculations. A FORTRAN code was programmed, which makes use of the EoS ePC-SAFT advanced (see the Supporting Information) to predict the fugacity coefficients of the charged species and neutral components in each phase. The code performs basically a ( $p, T, \underline{z}$ )-flash calculation involving only liquid phases (however, the extension to handle a vapor phase is straightforward); the code was used to perform LLE as well as LLLE calculations. The architecture of the algorithm is shown in Figure 1.

In a first step, a unimolar amount of feed of composition $\underline{z}_{\mathrm{F}}$ at given temperature $T$ and pressure $p$ is initialized. Phasestability analysis is performed by minimizing the TPDF of this feed composition $\underline{z}_{\mathrm{F}}$ (eqs 33-35) using a random-search algorithm to generate compositions $\bar{x}$ of the trial phase. Basically, the algorithm generates random values between 0 and 1 of each of the number of moles of $N^{\text {neut }}$ neutral components $\left[n_{1}, \ldots, n_{N^{\text {netu }}}\right]$ and $N^{\text {ch }}-1$ transformed variables $\left[\xi_{1}, \ldots, \xi_{N^{\text {ch }}-1}\right]$, from which a trial composition $\bar{x}^{\text {trial }}$ is calculated using eqs 34 and 35. The algorithm first performs a global search over the entire composition space while keeping track of the minimal value of $\operatorname{TPDF}\left(\bar{x}^{*}\right)$ and $\bar{x}^{*}$ during the search. After a fixed number of iterations, the optimal value $\bar{x}^{*}$ is improved by moving it using small random perturbations [ $\delta n_{1}$, $\ldots, \delta n_{N^{\text {net }}}$ ] and $\left[\delta \xi_{1}, \ldots, \delta \xi_{N^{\text {ch }}-1}\right]$ and accepting only new values of $\bar{x}$ if they lead to a decrease of the TPDF. The maximal allowed length $\left\|\delta n_{i}, \delta \xi_{j}\right\|$ of the perturbation is decreased, as the algorithm proceeds, according to a predefined reduction protocol. If instability is detected (confirmed by the negative value of $\operatorname{TPDF}\left(\bar{x}^{*}\right)$ when the search is terminated), then a liquid phase with composition $\bar{x}^{*}$ is added, and the amount of the new phase separated from the feed is calculated by minimizing the Gibbs energy of the system, using the phase

![](https://cdn.mathpix.com/cropped/790aed8a-6454-4c32-b8d0-0b109e9556ff-06.jpg?height=1139&width=644&top_left_y=202&top_left_x=1188)
Figure 1. Architecture of the proposed algorithm for the determination of the number and composition of coexisting liquid phases at equilibrium in a system of feed composition $\underline{Z}_{\mathrm{F}}$ at given temperature $T$ and pressure $p . \pi$ represents the number of phases (which is, at the same time, an output of the routine), and $j$ represents the phase which gives the largest split after phase stability analysis. This phase is split into two phases $\left(\underline{x}^{(j)} \rightarrow\left[\underline{x}^{(\pi)}, \underline{x}^{(\pi+1)}\right]\right)$ which are then added to the list of compositions of all of the phases.

amount as the optimization variable. Using an EoS like it was done in this work, the absolute Gibbs energy $g$ was not used; rather, the modeling relies on a value given by $\hat{g}= R T \sum_{k=1}^{\pi} \sum_{i=1}^{N} \alpha^{(k)} x_{i}^{(k)} \ln f_{i}^{(k)}=g-\sum_{i} z_{i} \mu_{0 i}^{i d}$, with $\alpha^{(k)}$ being the molar fraction of phase $k$ and $\mu_{0 i}^{i d}$ the chemical potential of pure component $i$ as hypothetical ideal gas at the temperature $T$ and pressure $p$ of the considered system; this shall just be mentioned for the sake of completeness. Using $\hat{g}$ and $g$ yields the same results. The composition of the remaining phase results from the optimal amount of the new phase. The concentration estimates provided by the phase-stability analysis are used as a starting value for the subsequent phase-equilibrium calculation. Within the phase-equilibria calculation, the system of equations given by eqs 22 and 23 and eqs 31 and 32 for a unimolar system is solved numerically by using nonlinear equations based on a modified Powell algorithm available in the IMSL FORTRAN libraries. Phasestability analysis is then applied to the light phase (denoted by the superscript (1)). If instability is detected, then the number of phases is increased by one, and the composition of one of the phases (in our approach, this is the phase which gives the highest split by minimizing the Gibbs energy $\hat{g}$ when the trial phase is separated) is substituted by the estimate of both phases provided by the phase-stability analysis. By adding this
new estimate, phase-equilibrium calculations are performed again, and this procedure is repeated until the system is stable. It should be noted that the employed procedure based on successive phase-stability analysis prior to solving the phaseequilibrium equations was originally proposed by Michelsen ${ }^{55,56}$ and followed by other works ${ }^{57,58}$ where it proved to be reliable for determining the exact number and composition of phases at equilibrium. However, the thermodynamic framework (section 2) and the variable reformulation (eqs 28-30) are not restricted to one calculation method and can be applied to all of the other methods already proposed in the literature for non-charged components, including hybrid ${ }^{59}$ and global methods. ${ }^{60}$ In the following, two case studies are presented to illustrate the application of the proposed methodology using the designed algorithm. Please note that this approach uses the EoS ePC-SAFT advanced to predict fugacities; i.e., the results of the case studies do depend on the specific model employed to estimate the fugacities, whereas the specific computational approach will determine whether the correct solution of the phase-equilibrium problem will be found.

### 4.1. Case Study 1: LLLE of Ternary Systems Containing

ILs. Three-phase liquid-liquid-liquid systems (LLLE) are extremely interesting from an industrial point of view due to their importance in several processes involving microemulsions. ${ }^{61}$ Such processes include enhanced oil recovery, ${ }^{62-67}$ microextraction of lipids and polymers, ${ }^{68}$ wastewater treatment, ${ }^{69}$ and phase-transfer catalytic processes. ${ }^{70,71}$ If the occurrence of LLLE is not known in the design of technical systems, this might lead to operational inconsistencies within the equipment. ${ }^{61}$ In this work, two different ternary systems containing ILs were chosen as liquid three-phase systems to test the proposed computational methodology; the considered systems and the references are shown in Table 1.

Table 1. Overview of the Systems Showing a Three-Phase LLLE Considered in This Work
| system | $T / \mathrm{K}$ | $p /$ bar | ref. |
| :---: | :---: | :---: | :---: |
| water $+\left[\mathrm{P}_{66614}\right][\mathrm{DCA}]+n$-dodecane | 298.15 | 1.00 | 66 |
| water $+\left[\mathrm{P}_{66614}\right][\mathrm{DCA}]+n$-hexane | 298.15 | 1.00 | 65 |


The investigated IL has been modeled as a fully dissociated salt using the EoS ePC-SAFT advanced, ${ }^{23,24}$ which includes a concentration-dependent dielectric constant and the Born term as a further contribution to the residual Helmholtz energy; this version has been used in a recent work to predict liquid-liquid equilibria containing dissolved ions. ${ }^{72}$ The pure-component parameters of the ionic species and the neutral components used in this work are listed in the Supporting Information. The relative dielectric constant $\varepsilon_{r}$ of the IL-cation and the IL-anion was assumed to be $\varepsilon_{r}=11$, according to a previous work. ${ }^{72}$ The $\varepsilon_{r}$ value of alkanes might be estimated to have a value of $\varepsilon_{r} =2$ (cf. ref 73), and thus, this constant value was assumed to be valid also for dodecane. According to a previous work, ${ }^{22}$ the value for $\varepsilon_{r}$ of pure water was calculated from the correlation of Floriano et al. ${ }^{74}$ The resulting mean dielectric constant of the mixture was calculated by a linear weight-fraction-based mixing rule from the dielectric constants of the pure components.

The complete calculation procedure, which leads to the calculation of the three-phase LLLE will be shown in detail in the following only for the ternary system water + $\left[\mathrm{P}_{66614}\right][\mathrm{DCA}]+n$-dodecane ("dodec"). This ternary system is classified as a Treybal type III (cf. ref 75), in which each
binary subsystem shows partial miscibility and a three-liquid phase miscibility gap is reported. To test the computational approach, a feed mass-fraction-based composition of $\underline{w}_{Z}= \left(w_{\text {water }}, w_{\text {dodec }} w_{\text {IL }}\right)=(0.0903,0.1220,0.7877)$ at 298.15 K and 1 bar was chosen in the present work. All binary interaction parameters were set to zero for the prediction of the fugacity coefficients with ePC-SAFT advanced. The homogeneous feed at the given conditions has a value of $\hat{g}=-92802.136 \mathrm{~J} \cdot \mathrm{~mol}^{-1}$, according to ePC-SAFT advanced. For this feed, the phasestability analysis predicts instability of the homogeneous mixture, and thus, two liquid phases were initialized in a first step. The phase-equilibrium calculation yielded two points $\underline{w}^{\mathrm{I}} =\left(w_{\text {water }}=0.00142, w_{\text {dodec }}=0.99858, w_{\text {IL }}=6.96 \times 10^{-32}\right)$ and $\underline{w}^{\mathrm{II}}=\left(w_{\text {water }}=0.18445, w_{\text {dodec }}=0.00835, w_{\text {IL }}=0.80720\right)$ with $\hat{g} =-93702.172 \mathrm{~J} \cdot \mathrm{~mol}^{-1}$ of the whole system and fugacities given by Table 2.

Table 2. Fugacities of the Neutral Components and Mean Ionic Fugacity of the IL in the Two Liquid Phases at Equilibrium, Obtained after the First Split of the Homogeneous Feed $\underline{w}_{Z}$ at 298.15 K and 1 bar, and Obtained with ePC-SAFT Advanced for the System Water + [ $\mathbf{P}_{66614}$ ][DCA] + Dodecane Using the Pure-Component Parameters according to Tables S1 and S2
|  | phase (I) | phase (II) |
| :--- | :--- | :--- |
| $\ln \left(f_{\text {water }} /\right.$ bar $)$ | -3.282 | -3.282 |
| $\ln \left(f_{\text {dodecane }} /\right.$ bar $)$ | -8.570 | -8.570 |
| $\ln \left(f_{P_{66614}+\mathrm{DCA}-} /\right.$ bar $)$ | -187.840 | -187.840 |


In a next step, phase-stability analysis was applied to phase (I) to ascertain the stability of the two-phase system. Instability was detected, while the second phase was split by removing the amount of the trial phase, which minimizes the Gibbs energy of the system. Thus, estimates of the two phases resulting from splitting the second phases are used, together with the composition of the first phase, as an initial estimate to a three-phase-equilibria calculation. The resulting three phases at equilibrium have the following compositions:

$$
\begin{aligned}
& \underline{w}^{\mathrm{I}}=\left(w_{\text {water }}=0.00109, w_{\text {dodec }}=0.99891, w_{\mathrm{IL}}\right. \\
& \left.\quad=4.43 \times 10^{-32}\right) \\
& \underline{w}^{\mathrm{II}}=\left(w_{\text {water }} \approx 1.00000, w_{\text {dodec }}=2.16 \times 10^{-8}, w_{\mathrm{IL}}\right. \\
& \left.\quad=6.68 \times 10^{-7}\right) \\
& \underline{w}^{\mathrm{III}}=\left(w_{\text {water }}=0.14000, w_{\text {dodec }}=0.00764, w_{\mathrm{IL}}\right. \\
& \quad=0.85236)
\end{aligned}
$$

The whole mixture has a value of the calculated $\hat{g}= -93746.040 \mathrm{~J} \cdot \mathrm{~mol}^{-1}$. Respective fugacities are listed in Table 3.

Due to the extremely low IL content in the aqueous phase ( $w_{\text {IL }}=6.68 \times 10^{-7} \mathrm{~g} / \mathrm{g}$ ) and in the organic phase ( $w_{\text {IL }}=4.43 \times 10^{-32} \mathrm{~g} / \mathrm{g}$ ), the mean ionic fugacity $f_{P_{66614} / \text { DCA- }}$ deviates slightly between the three phases after convergence of the algorithm. However, the same results could be obtained by simply ignoring the equality of the mean ionic chemical potential and forcing the IL to be entirely present in the IL-rich phase, without the composition of the three phases being even quantitatively different after convergence of the routine.

Table 3. Fugacities of the Neutral Components and Mean Ionic Fugacity of the IL in the Three Liquid Phases at Equilibrium at 298.15 K and 1 bar , Obtained with ePCSAFT Advanced for the System Water $+\left[\mathbf{P}_{66614}\right][$ DCA $]+$ Dodecane Using the Pure-Component Parameters according to Tables S1 and S2 ${ }^{\boldsymbol{a}}$
|  | phase (I) | phase (II) | phase (III) |
| :--- | :--- | :--- | :--- |
| $\ln \left(f_{\text {water }} /\right.$ bar $)$ | -3.458 | -3.458 | -3.458 |
| $\ln \left(f_{\text {dodecane }} /\right.$ bar $)$ | -8.567 | -8.567 | -8.567 |
| $\ln \left(f_{P_{66614}+\mathrm{DCA}-} /\right.$ bar $)$ | -187.323 | -187.330 | -187.323 |


${ }^{a}$ The equilibrium composition is given in the text.
Results of the prediction are compared with experimental data in the ternary diagram shown in Figure 2.

To test the possibility of correlating the experimental data using the algorithm, binary interaction parameters between ILcation and the neutral components were adjusted to LLE of the binary subsystems (results not shown). The respective prediction of the ternary systems is shown in Figure 3.

Using binary interaction parameters obtained from data of the binary subsystems, the tie-lines in the two-phase region close to the binary mixture water-IL could be modeled quantitatively correctly. However, according to the calculations, the length of the tie-lines in the second two-phase region (close to the binary mixture dodecane-IL) rapidly decreases upon adding water to the system, while the length of the experimental tie-lines remains almost constant. Due to this behavior, the calculated three-phase region is smaller than the experimentally expected one. Nevertheless, ePC-SAFT advanced shows satisfactory performance, being the only electrolyte version of PC-SAFT able to predict a three-phase region with the used pure-component parameters (without binary interaction parameters, see Figure 2).
4.2. Case Study 2: LLE of Water-Organic Solvents Containing Multiple Salts. Salting-out is of great interest in the industry for designing separation processes containing electrolytes. Mixtures of water and short-chain alcohols such as ethanol or propanol are known to form a miscibility gap upon addition of salts; ${ }^{76-80}$ thus, such systems are considered to be alternative extraction systems for conventional polymer-based ATPS due to their lower cost and much lower viscosity. ${ }^{81}$

Experimental LLE data of ternary water-alcohol-salt systems are available in the literature. However, experimental LLE data of quaternary systems water-alcohol-salt1-salt2 are relatively scarce, probably due to the increasing effort to characterize such complex systems experimentally. In this work, four quaternary systems were considered containing water, either 1butanol or 1 -propanol, and two additional salts (e.g., NaCl , KCl , and $\mathrm{NH}_{4} \mathrm{Cl}$ ). The systems and the experimental conditions are listed in Table 4.

The salts have been modeled as fully dissociated using ePCSAFT advanced. The pure-component parameters of the neutral components and of the ions are shown in the Supporting Information in Tables S1 and S2. All of the binary interaction parameters between ions and the neutral components were set to zero. Binary interaction parameters between water and alcohol as well as between cations and the chloride anion were inherited from previous work and are given in Table S3 and Table S4 in the Supporting Information. The dielectric constant $\varepsilon_{r}$ of the mixture was estimated from the dielectric constant of the pure components using a combined weight-fraction-based mixing rule for the salt-free system and a molar-fraction-based mixing rule between the salt-free systems and the ions, according to our previous work. ${ }^{72}$ The value of the relative dielectric constant $\varepsilon_{r}$ of the studied ions was assumed to be $\varepsilon_{r}=8$ according to a previous work. ${ }^{23}$ The relative dielectric constant $\varepsilon_{r}$ of 1 -butanol was inherited from existing experimental data at 1 bar $^{85-87}$ as used in a previous work, ${ }^{72}$ and $\varepsilon_{r}$ of 1 -propanol was taken as $\varepsilon_{r}=$ 20.47 according to experimental data at 1 bar and $298.15 \mathrm{~K} .^{88}$ In the following, the step-by-step procedure of the developed algorithm will be explained for the system water +1 -butanol + $\mathrm{NaCl}+\mathrm{KCl}$ with mass-fraction-based feed compositions of the respective components $\underline{z}=\left(w_{\text {water }}=0.8094, w_{\text {butanol }}=0.1728\right.$, $w_{\mathrm{NaCl}}=0.0054, w_{\mathrm{KCl}}=0.0124$ ). Since three ions are present in the system, the number of transformed variables $\xi_{i}$ as well as of mean ionic fugacities $f_{ \pm}$needed to define the phase-equilibrium problem is two. Those are defined between cation-anion pairs $\mathrm{Na}^{+} / \mathrm{Cl}^{-}$and $\mathrm{K}^{+} / \mathrm{Cl}^{-}$using the proposed methodology. According to ePC-SAFT advanced, the feed with the above composition has a value of $\hat{g}=-27361.317 \mathrm{~J} \cdot \mathrm{~mol}^{-1}$. The phase-stability analysis confirmed instability of the given feed and provides initial values to the phase-equilibrium calculation

![](https://cdn.mathpix.com/cropped/790aed8a-6454-4c32-b8d0-0b109e9556ff-08.jpg?height=553&width=1239&top_left_y=1853&top_left_x=434)
Figure 2. Gibbs energy surface $\hat{g}$ (left) of the mixture water $+\left[\mathrm{P}_{66614}\right][\mathrm{DCA}]+$ dodecane at 298.15 K and 1 bar and respective ternary phase diagram (right) showing the experimental data (gray symbols connected by solid lines) and the predictions (white symbols connected by dashed lines) of the three-phase region. The feed composition used as an initial step for the algorithm is represented by the triangle. The predictions are a priori predictions with ePC-SAFT advanced using the pure-component parameters according to Tables S1 and S2 and without the use of any binary interaction parameters.

![](https://cdn.mathpix.com/cropped/790aed8a-6454-4c32-b8d0-0b109e9556ff-09.jpg?height=540&width=1239&top_left_y=199&top_left_x=434)
Figure 3. Ternary phase diagram of the mixture water $+\left[\mathrm{P}_{66614}\right][\mathrm{DCA}]+$ dodecane (left) and water $+\left[\mathrm{P}_{66614}\right][\mathrm{DCA}]+$ hexane (right) at 298.15 K and 1 bar. Gray symbols connected by solid lines represent experimental tie-lines, while white symbols connected by dot-dashed lines are ePCSAFT advanced modeling results. For the ePC-SAFT advanced modeling, the pure-component parameters according to Tables S1 and S2 were used and the following binary interaction parameters were used: $k_{\text {water }-\mathrm{P}_{66614}^{+}}=0.2, k_{\text {hexane }-\mathrm{P}_{66614}^{+}}=-0.175$, and $k_{\text {dodec }-\mathrm{P}_{66614}^{+}}=-0.125$; these were fit to data of binary systems (no parameters were fitted to the ternary system).

Table 4. Overview about the Systems Considered in This Work Showing a Two-Phase LLE with Two Salts
| system | $T / \mathrm{K}$ | $p /$ bar | ref. |
| :--- | :---: | :---: | :---: |
| water + 1-butanol $+\mathrm{KCl}+\mathrm{NH}_{4} \mathrm{Cl}$ | 298.15 | 1.00 | 82 |
| water + 1-butanol $+\mathrm{NaCl}+\mathrm{KCl}$ | 298.15 | 1.00 | 83 |
| water + 1-propanol $+\mathrm{KCl}+\mathrm{NH}_{4} \mathrm{Cl}$ | 298.15 | 1.00 | 82 |
| water + 1-propanol $+\mathrm{NaCl}+\mathrm{KCl}$ | 298.15 | 1.00 | 84 |


with $\pi=2$ phases. The phase-equilibrium calculation converged to a solution with mole-fraction-based phase compositions $\underline{x}^{(1)}=\left(x_{\text {water }}=0.4426, x_{\text {butanol }}=0.5570, x_{\mathrm{NaCl}}=\right. \left.4.15 \times 10^{-5}, x_{\mathrm{KCl}}=4.20 \times 10^{-4}\right)$ and $\underline{x}^{(2)}=\left(x_{\text {water }}=0.9627\right.$, $\left.x_{\text {butanol }}=0.0122, x_{\mathrm{NaCl}}=0.0076, x_{\mathrm{KCl}}=0.0174\right)$ and a value of $\hat{g} =-27479.860 \mathrm{~J} \cdot \mathrm{~mol}^{-1}$, which is lower than that for the homogeneous feed $\left(\hat{g}=-27361.317 \mathrm{~J} \cdot \mathrm{~mol}^{-1}\right)$. Furthermore, the two-phase system is stable according to the phase-stability test; thus, the obtained compositions are the solution of the phase-equilibrium problem for this feed composition under the given conditions. The values of the fugacities in both phases are listed in Table 5.

Calculation results of the four systems considered in this case study are illustrated in Figure 4.

To conclude, ePC-SAFT advanced shows excellent accuracy in predicting the LLE of systems of water-organic solvents containing one or more distributed salts. Using the proposed

Table 5. Fugacities of the Neutral Components and Mean Ionic Fugacity of the Two Salts in the Two Liquid Phases at Equilibrium Obtained after the First Split of the Feed $\underline{w}_{Z}$ at 298.15 K and 1 bar, Obtained with ePC-SAFT Advanced for the System Water + 1-Butanol + NaCl + KCl Using the Pure-Component Parameters according to Tables S1 and S2 and the Binary Parameters between Water and Butanol (Table S3) and Those among Ions (Table S4) ${ }^{a}$
|  | organic phase (1) | aqueous phase (2) |
| :--- | :---: | :---: |
| $\ln \left(f_{\text {water }} /\right.$ bar $)$ | -3.521 | -3.521 |
| $\ln \left(f_{\text {butanol }} /\right.$ bar $)$ | -5.088 | -5.088 |
| $\ln \left(f_{ \pm, \mathrm{K}^{+}} / \mathrm{Cl}^{-} /\right.$bar $)$ | -206.733 | -206.733 |
| $\ln \left(f_{ \pm, \mathrm{Na}^{+}} / \mathrm{Cl}^{-} /\right.$bar $)$ | -224.891 | -224.891 |


[^1]methodology allows predicting correctly the appearance of the kind and number of multiple liquid phases. In the case of ionicliquid systems, quantitative agreement is only possible by using binary parameters, which were fitted in this work exclusively to LLE data of binary subsystems; nevertheless, this opens the possibility toward predictions of phase equilibria in higher systems using the binary parameters determined from experimental data of the binary subsystems. Ultimately, experimental data in systems with more than two components seem to not be required for accurately predicting the LL phase behavior of higher systems; these findings are valid for the systems shown in this work, and validation is required also for other systems in future work.

## - CONCLUSION

In this work, a general mathematical formulation of multiphase equilibria with multiple distributed ions has been derived with the goal to compute multiphase phase equilibria of mixed solvents and mixed electrolytes systems, and the methodology has been applied to design an algorithm for multiphase equilibrium calculation. The developed algorithm is an equation-solving approach coupled with the tangent-plane stability analysis, but any calculation methodology originally developed for the calculation of phase equilibria of noncharged components (such as global optimization or hybrid methods) can be extended to electrolyte systems. However, it must be pointed out that the calculation efficiency and the numerical robustness heavily rely on the chosen calculation methodology, which depends on the specific considered phaseequilibrium problem. The proposed algorithm was implemented in a FORTRAN routine using the equation of state ePC-SAFT advanced to estimate the fugacity of each species. We conclude with two major notes: First, the proposed algorithm allowed correctly predicting the appearance of the kind and number of multiple liquid phases in chosen test systems. Second, the applied thermodynamic model used within the algorithm (ePC-SAFT advanced) showed excellent accuracy in predicting the LLE of systems of water-organic solvents containing one or more inorganic salts or ionic liquids. In the case of mixtures with inorganic salts, quantitatively correct a-priori predictions without using any binary parameters were possible. In the case of mixtures with ionic

![](https://cdn.mathpix.com/cropped/790aed8a-6454-4c32-b8d0-0b109e9556ff-10.jpg?height=1261&width=1403&top_left_y=193&top_left_x=351)
Figure 4. Quaternary phase diagrams based on weight fraction showing the LLE of four different systems composed of water, an organic solvent, and two inorganic salts at 298.15 K and 1 bar. White symbols connected by dotted lines are experimental tie-lines, ${ }^{82-84}$ and black symbols connected by solid lines are predictions using ePC-SAFT advanced. Parameters were taken from the literature (see the Supporting Information), and no further binary interaction parameter was fitted to predict the phase compositions.

liquids, quantitative agreement is only possible by using binary parameters, which were fitted in this work exclusively to LLE data of binary subsystems. In the future, the developed algorithm will be used to model even more complex systems.

## - ASSOCIATED CONTENT

## (5) Supporting Information

The Supporting Information is available free of charge at https://pubs.acs.org/doi/10.1021/acs.jced.1c00866.

Details about the EoS used for phase equilibrium calculation as well as the pure component and binary interaction parameter used in this work (PDF)

## - AUTHOR INFORMATION

## Corresponding Author

Christoph Held - Laboratory of Thermodynamics, Department of Biochemical and Chemical Engineering, TU Dortmund, 44277 Dortmund, Germany; © orcid.org/0000-0003-1074-177X; Email: christoph.held@tu-dortmund.de

## Authors

Moreno Ascani - Laboratory of Thermodynamics, Department of Biochemical and Chemical Engineering, TU Dortmund, 44277 Dortmund, Germany

Gabriele Sadowski - Laboratory of Thermodynamics, Department of Biochemical and Chemical Engineering, TU Dortmund, 44277 Dortmund, Germany; orcid.org/0000-0002-5038-9152

Complete contact information is available at:
https://pubs.acs.org/10.1021/acs.jced.1c00866

## Notes

The authors declare no competing financial interest.

## - ACKNOWLEDGMENTS

Funded by the Deutsche Forschungsgemeinschaft (DFG, German Research Foundation) under Germany's Excellence Strategy - EXC 2033-390677874 - RESOLV. We dedicate this work to Joan Brennecke. Joan, you are a fantastic teacher and scientist. All the best for your future.

## - REFERENCES

(1) Kontogeorgis, G. M.; Folas, G. K. Thermodynamic models for industrial applications: from classical and advanced mixing rules to association theories; John Wiley \& Sons: 2009.
(2) Mohammad, S.; Grundl, G.; Müller, R.; Kunz, W.; Sadowski, G.; Held, C. Influence of electrolytes on liquid-liquid equilibria of water/ 1-butanol and on the partitioning of 5-hydroxymethylfurfural in water/1-butanol. Fluid Phase Equilib. 2016, 428, 102-111.
(3) Mohammad, S.; Held, C.; Altuntepe, E.; Köse, T.; Sadowski, G. Influence of Salts on the Partitioning of 5-Hydroxymethylfurfural in Water/MIBK. J. Phys. Chem. B 2016, 120, 3797-3808.
(4) Bülow, M.; Gerek Ince, N.; Hirohama, S.; Sadowski, G.; Held, C. Predicting Vapor-Liquid Equilibria for Sour-Gas Absorption in Aqueous Mixtures of Chemical and Physical Solvents or Ionic Liquids with ePC-SAFT. Ind. Eng. Chem. Res. 2021, 60, 6327-6336.
(5) Pabsch, D.; Held, C.; Sadowski, G. Modeling the $\mathrm{CO}_{2}$ Solubility in Aqueous Electrolyte Solutions Using ePC-SAFT. J. Chem. Eng. Data 2020, 65, 5768-5777.
(6) Buelow, M.; Ji, X.; Held, C. Incorporating a concentrationdependent dielectric constant into ePC-SAFT. An application to binary mixtures containing ionic liquids. Fluid Phase Equilib. 2019, 492, 26-33.
(7) Nordness, O.; Brennecke, J. F. Ion dissociation in ionic liquids and ionic liquid solutions. Chem. Rev. 2020, 120, 12873-12902.
(8) Huckel, E.; Debye, P. Zur Theorie der Elektrolyte. I. Gefrierpunktserniedrigung und verwandte Erscheinungen. Phys. Z. 1923, 24, 185-206.
(9) Blum, Lo Mean spherical model for asymmetric electrolytes: I. Method of solution. Mol. Phys. 1975, 30, 1529-1535.
(10) Maribo-Mogensen, B.; Kontogeorgis, G. M.; Thomsen, K. Comparison of the Debye-Hückel and the Mean Spherical Approximation Theories for Electrolyte Solutions. Ind. Eng. Chem. Res. 2012, 51, 5353-5363.
(11) Pitzer, K. S. Thermodynamics of electrolytes. I. Theoretical basis and general equations. J. Phys. Chem. 1973, 77, 268-277.
(12) Pitzer, K. S.; Mayorga, G. Thermodynamics of electrolytes. II. Activity and osmotic coefficients for strong electrolytes with one or both ions univalent. J. Phys. Chem. 1973, 77, 2300-2308.
(13) Chen, C.-C.; Britt, H. I.; Boston, J. F.; Evans, L. B. Local composition model for excess Gibbs energy of electrolyte systems. Part I: Single solvent, single completely dissociated electrolyte systems. AIChE J. 1982, 28, 588-596.
(14) Chen, C.-C.; Evans, L. B. A local composition model for the excess Gibbs energy of aqueous electrolyte systems. AIChE J. 1986, 32, 444-454.
(15) Thomsen, K.; Rasmussen, P.; Gani, R. Correlation and prediction of thermal properties and phase behaviour for a class of aqueous electrolyte systems. Chem. Eng. Sci. 1996, 51, 3675-3683.
(16) Thomsen, K.; Rasmussen, P. Modeling of vapor-liquid-solid equilibrium in gas-aqueous electrolyte systems. Chem. Eng. Sci. 1999, 54, 1787-1802.
(17) Wu, J.; Prausnitz, J. M. Phase equilibria for systems containing hydrocarbons, water, and salt: An extended Peng- Robinson equation of state. Ind. Eng. Chem. Res. 1998, 37, 1634-1643.
(18) Myers, J. A.; Sandler, S. I.; Wood, R. H. An equation of state for electrolyte solutions covering wide ranges of temperature, pressure, and composition. Ind. Eng. Chem. Res. 2002, 41, 3282-3297.
(19) Zhao, H.; Dos Ramos, M. C.; McCabe, C. Development of an equation of state for electrolyte solutions by combining the statistical associating fluid theory and the mean spherical approximation for the nonprimitive model. J. Chem. Phys. 2007, 126, 244503.
(20) Galindo, A.; Gil-Villegas, A.; Jackson, G.; Burgess, A. N. SAFTVRE: phase behavior of electrolyte solutions with the statistical associating fluid theory for potentials of variable range. J. Phys. Chem. B 1999, 103, 10272-10281.
(21) Cameretti, L. F.; Sadowski, G.; Mollerup, J. M. Modeling of aqueous electrolyte solutions with perturbed-chain statistical associated fluid theory. Ind. Eng. Chem. Res. 2005, 44, 3355-3362.
(22) Held, C.; Reschke, T.; Mohammad, S.; Luza, A.; Sadowski, G. ePC-SAFT revised. Chem. Eng. Res. Des. 2014, 92, 2884-2897.
(23) Bülow, M.; Ascani, M.; Held, C. ePC-SAFT advanced-Part I: Physical meaning of including a concentration-dependent dielectric constant in the born term and in the Debye-Hückel theory. Fluid Phase Equilib. 2021, 535, 112967.
(24) Bülow, M.; Ascani, M.; Held, C. ePC-SAFT advanced-Part II: Application to Salt Solubility in Ionic and Organic Solvents and the Impact of Ion Pairing. Fluid Phase Equilib. 2021, 537, 112989.
(25) Born, M. Volumen und hydratationswärme der ionen. Zeitschrift für physik 1920, 1, 45-48.
(26) Kontogeorgis, G. M.; Maribo-Mogensen, B.; Thomsen, K. The Debye-Hückel theory and its importance in modeling electrolyte solutions. Fluid Phase Equilib. 2018, 462, 130-152.
(27) Maribo-Mogensen, B.; Thomsen, K.; Kontogeorgis, G. M. An electrolyte CPA equation of state for mixed solvent electrolytes. AIChE J. 2015, 61, 2933-2950.
(28) Schlaikjer, A.; Thomsen, K.; Kontogeorgis, G. M. eCPA: An ion-specific approach to parametrization. Fluid Phase Equilib. 2018, 470, 176-187.
(29) Held, C. Thermodynamic g E Models and Equations of State for Electrolytes in a Water-Poor Medium: A Review. J. Chem. Eng. Data 2020, 65, 5073-5082.
(30) Gautam, R.; Seider, W. D. Computation of phase and chemical equilibrium: Part III. Electrolytic solutions. AIChE J. 1979, 25, 10061015.
(31) White, W. B.; Johnson, S. M.; Dantzig, G. B. Chemical Equilibrium in Complex Mixtures. J. Chem. Phys. 1958, 28, 751-755.
(32) Simoni, L. D.; Brennecke, J. F.; Stadtherr, M. A. Asymmetric Framework for Predicting Liquid- Liquid Equilibrium of Ionic Liquid- Mixed-Solvent Systems. 1. Theory, Phase Stability Analysis, and Parameter Estimation. Ind. Eng. Chem. Res. 2009, 48, 7246-7256.
(33) Simoni, L. D.; Chapeaux, A.; Brennecke, J. F.; Stadtherr, M. A. Asymmetric Framework for Predicting Liquid- Liquid Equilibrium of Ionic Liquid- Mixed-Solvent Systems. 2. Prediction of Ternary Systems. Ind. Eng. Chem. Res. 2009, 48, 7257-7265.
(34) Zuend, A.; Seinfeld, J. H. A practical method for the calculation of liquid-liquid equilibria in multicomponent organic-waterelectrolyte systems using physicochemical constraints. Fluid Phase Equilib. 2013, 337, 201-213.
(35) Tsanas, C.; Stenby, E. H.; Yan, W. Calculation of multiphase chemical equilibrium in electrolyte solutions with non-stoichiometric methods. Fluid Phase Equilib. 2019, 482, 81-98.
(36) Medeiros, F. d. A.; Stenby, E. H.; Yan, W. RAND-based geochemical equilibrium algorithms with applications to underground geological storage of CO2. Advances in Water Resources 2021, 152, 103918.
(37) Bethke, C. M. Geochemical and Biogeochemical Reaction Modeling; Cambridge University Press: Cambridge, U.K., 2007.
(38) Duan, Z.; Sun, R. An improved model calculating CO2 solubility in pure water and aqueous NaCl solutions from 273 to 533 K and from 0 to 2000 bar. Chem. Geol. 2003, 193, 257-271.
(39) Duan, Z.; Sun, R.; Zhu, C.; Chou, I.-M. An improved model for the calculation of CO 2 solubility in aqueous solutions containing Na +, $\mathrm{K}+, \mathrm{Ca} 2+, \mathrm{Mg} 2+, \mathrm{Cl}-$, and SO $42-$. Marine chemistry 2006, 98, 131-139.
(40) Kulik, D. A.; Wagner, T.; Dmytrieva, S. V.; Kosakowski, G.; Hingerl, F. F.; Chudnenko, K. V.; Berner, U. R. GEM-Selektor geochemical modeling package: revised algorithm and GEMS3K numerical kernel for coupled simulation codes. Computational Geosciences 2013, 17, 1-24.
(41) Leal, A. M.M.; Blunt, M. J.; LaForce, T. C. Efficient chemical equilibrium calculations for geochemical speciation and reactive transport modelling. Geochim. Cosmochim. Acta 2014, 131, 301-322.
(42) Li, D.; Duan, Z. The speciation equilibrium coupling with phase equilibrium in the $\mathrm{H} 2 \mathrm{O}-\mathrm{CO} 2-\mathrm{NaCl}$ system from 0 to 250 C , from 0 to 1000 bar, and from 0 to 5 molality of NaCl . Chem. Geol. 2007, 244, 730-751.
(43) Parkhurst, D. L.; Appelo, C. A. J. Description of input and examples for PHREEQC version 3: a computer program for speciation, batch-reaction, one-dimensional transport, and inverse geochemical calculations; Techniques and Methods 6-A43; Reston, VA, 2013.
(44) Xu, T.; Sonnenthal, E.; Spycher, N.; Pruess, K. TOUGHREACT User's Guide: A Simulation Program for Non-isothermal Multiphase Reactive Geochemical Transport in Variably Saturated Geologic Media, V1.2.1, 2008.
(45) Tsanas, C.; Hemptinne, J.-C. de; Mougin, P. Calculation of phase and chemical equilibrium for multiple ion-containing phases including stability analysis. Chem. Eng. Sci. 2022, 248, 117174.
(46) Großmann, C.; Maurer, G. On the calculation of phase equilibria in aqueous two-phase systems containing ionic solutes. Fluid Phase Equilib. 1995, 106, 17-25.
(47) Luckas, M.; Krissmann, J. Thermodynamik der Elektrolytlösungen: eine einheitliche Darstellung der Berechnung komplexer Gleichgewichte; Springer-Verlag: 2013.
(48) Gmehling, J.; Kolbe, B. Thermodynamik; VCH: Weinheim, Germany, 1992.
(49) Held, C.; Cameretti, L. F.; Sadowski, G. Modeling aqueous electrolyte solutions: Part 1. Fully dissociated electrolytes. Fluid Phase Equilib. 2008, 270, 87-96.
(50) Held, C.; Sadowski, G. Modeling aqueous electrolyte solutions. Part 2. Weak electrolytes. Fluid Phase Equilib. 2009, 279, 141-148.
(51) Haynes, C. A.; Carson, J.; Blanch, H. W.; Prausnitz, J. M. Electrostatic potentials and protein partitioning in aqueous two-phase systems. AIChE J. 1991, 37, 1401-1409.
(52) Pratt, L. R. Contact potentials of solution interfaces: phase equilibrium and interfacial electric fields. J. Phys. Chem. 1992, 96, 2533.
(53) Ung, S.; Doherty, M. F. Theory of phase equilibria in multireaction systems. Chem. Eng. Sci. 1995, 50, 3201-3216.
(54) Ung, S.; Doherty, M. F. Vapor-liquid phase equilibrium in systems with multiple chemical reactions. Chem. Eng. Sci. 1995, 50, 23-48.
(55) Michelsen, M. L. The isothermal flash problem. Part I. Stability. Fluid Phase Equilib. 1982, 9, 1-19.
(56) Michelsen, M. L. The isothermal flash problem. Part II. Phasesplit calculation. Fluid Phase Equilib. 1982, 9, 21-40.
(57) Stateva, R. P.; Cholakov, G. S.; Galushko, A. A.; Wakeham, W.
A. A powerful algorithm for liquid-liquid-liquid equilibria predictions and calculations. Chem. Eng. Sci. 2000, 55, 2121-2129.
(58) Lucia, A.; Padmanabhan, L.; Venkataraman, S. Multiphase equilibrium flash calculations. Comput. Chem. Eng. 2000, 24, 25572569.
(59) Bausa, J.; Marquardt, W. Quick and reliable phase stability test in VLLE flash calculations by homotopy continuation. Comput. Chem. Eng. 2000, 24, 2447-2456.
(60) McDonald, C. M.; Floudas, C. A. Global optimization for the phase and chemical equilibrium problem: application to the NRTL equation. Comput. Chem. Eng. 1995, 19, 1111-1139.
(61) Bairagya, P.; Kundu, D.; Banerjee, T. A priori prediction of complex liquid-liquid-liquid equilibria in organic systems using a continuum solvation model. Phys. Chem. Chem. Phys. 2020, 22, 22023-22034.
(62) Lago, S.; Francisco, M.; Arce, A.; Soto, A. Enhanced oil recovery with the ionic liquid trihexyl (tetradecyl) phosphonium chloride: a phase equilibria study at $75^{\circ}$ C. Energy Fuels 2013, 27, 5806-5810.
(63) Lago, S.; Rodriguez, H.; Khoshkbarchi, M. K.; Soto, A.; Arce, A. Enhanced oil recovery using the ionic liquid trihexyl (tetradecyl) phosphonium chloride: phase behaviour and properties. RSC Adv. 2012, 2, 9392-9397.
(64) Lago, S.; Rodríguez-Cabo, B.; Arce, A.; Soto, A. Water/oil/[P6, 6, 6, 14][NTf2] phase equilibria. J. Chem. Thermodyn. 2014, 75, 6368.
(65) Rodríguez-Escontrela, I.; Arce, A.; Soto, A.; Marcilla, A.; Olaya, M. d. M.; Reyes-Labarta, J. A. Correlation of three-liquid-phase equilibria involving ionic liquids. Phys. Chem. Chem. Phys. 2016, 18, 21610-21617.
(66) Rodriguez-Escontrela, I.; Rodriguez-Palmeiro, I.; Rodriguez, O.; Arce, A.; Soto, A. Liquid-liquid-liquid equilibria for water+[P6 6 6 14][DCA]+ dodecane ternary system. Fluid Phase Equilib. 2015, 405, 124-131.
(67) Rodríguez-Palmeiro, I.; Rodríguez, O.; Soto, A.; Held, C. Measurement and PC-SAFT modelling of three-phase behaviour. Phys. Chem. Chem. Phys. 2015, 17, 1800-1810.
(68) He, X.; Huang, K.; Yu, P.; Zhang, C.; XIie, K.; Li, P.; Wang, J.; AN, Z.; LIU, H. Liquid-liquid-liquid three phase extraction apparatus: operation strategy and influences on mass transfer efficiency. Chinese Journal of Chemical Engineering 2012, 20, 27-35.
(69) Zhang, C.; Huang, K.; Yu, P.; Liu, H. Ionic liquid based three-liquid-phase partitioning and one-step separation of Pt (IV), Pd (II) and Rh (III). Sep. Purif. Technol. 2013, 108, 166-173.
(70) Krueger, J. J.; Amiridis, M. D.; Ploehn, H. J. Modeling mass transfer and interfacial reactions in three liquid phase transfer catalysis. Ind. Eng. Chem. Res. 2001, 40, 3158-3163.
(71) Nardello-Rataj, V.; Caron, L.; Borde, C.; Aubry, J.-M. Oxidation in three-liquid-phase microemulsion systems using "Balanced Catalytic Surfactants. J. Am. Chem. Soc. 2008, 130, 1491414915.
(72) Ascani, M.; Held, C. Prediction of salting-out in liquid-liquid two-phase systems with ePC-SAFT: Effect of the Born term and of a concentration-dependent dielectric constant. Zeitschrift für anorganische und allgemeine Chemie 2021, 647, 1305.
(73) Sen, A. D.; Anicich, V. G.; Arakelian, T. Dielectric constant of liquid alkanes and hydrocarbon mixtures. J. Phys. D: Appl. Phys. 1992, 25, 516.
(74) Floriano, W. B.; Nascimento, M. A. C. Dielectric constant and density of water as a function of pressure at constant temperature. Brazilian Journal of Physics 2004, 34, 38-41.
(75) Treybal, R. E. Methods of Calculation II. Stagewise Contact, Multicomponent Systems. Liquid Extraction, 2nd ed.; McGraw-Hill Book Company: New York, 1963; pp 275-276.
(76) Arzideh, S. M.; Movagharnejad, K.; Pirdashti, M. Influence of the temperature, type of salt, and alcohol on phase diagrams of 2propanol+ inorganic salt aqueous two-phase systems: experimental determination and correlation. J. Chem. Eng. Data 2018, 63, 28132824.
(77) Guo, W.; Ma, J.; Wang, Y.; Han, J.; Li, Y.; Song, S. Liquidliquid equilibrium of aqueous two-phase systems composed of hydrophilic alcohols (ethanol/2-propanol/1-propanol) and MgSO4/ ZnSO 4 at (303.15 and 313.15) K and correlation. Thermochimica acta 2012, 546, 8-15.
(78) Katayama, H.; Kitagawa, K. Liquid- liquid phase equilibria of (1-propanol or 2-propanol+ water) containing dipotassium hydrogen phosphate. J. Chem. Eng. Data 2006, 51, 2103-2106.
(79) Nemati-Kande, E.; Shekaari, H.; Jafari, S. A. Liquid-liquid equilibrium of 1-propanol, 2-propanol, 2-methyl-2-propanol or 2butanol+ sodium sulfite+ water aqueous two phase systems. Fluid Phase Equilib. 2012, 329, 42-54.
(80) Wang, J.; Zhang, Y.; Wang, Y. Liquid- liquid equilibria for 1propanol (or 2-propanol)- water systems containing potassium fluoride. J. Chem. Eng. Data 2002, 47, 110-112.
(81) Garcia-Cano, J.; Gomis, A.; Font, A.; Gomis, V. Effect of temperature on the phase-separation ability of KCl in aqueous twophase systems composed of propanols: Determination of the critical temperature and extension of the results to other salts. J. Chem. Thermodyn. 2019, 136, 88-99.
(82) Chen, J.; Zhong, Y.; Han, J.; Su, M.; Shi, X. Liquid-liquid equilibria for water+ 1-propanol (or 1-butanol)+ potassium chloride+ ammonium chloride quaternary systems at 298. 15 K . Fluid Phase Equilibria 2015, 397, 50-57.
(83) Pedraza, R.; Ruiz, F.; Saquete, M. D.; Gomis, V. Liquid-liquid-solid equilibrium for the water+ sodium chloride+ potassium chloride +1 -butanol quaternary system at $25^{\circ}$ C. Fluid Phase Equilib. 2004, 216, 27-31.
(84) Pedraza, R.; Ruiz, F.; Saquete, M. D.; Gomis, V. Liquid—liquid and liquid-liquid-solid equilibrium for the water+ sodium chloride+ potassium chloride+ 1 -propanol quaternary system at $25^{\circ}$ C. Fluid Phase Equilib. 2004, 221, 97-101.
(85) Heyding, R. D.; Winkler, C. A. Solvent Effect on Iodide Exchange. Can. J. Chem. 1951, 29, 790-803.
(86) Dannhauser, W.; Cole, R. H. Dielectric Properties of Liquid Butyl Alcohols. J. Chem. Phys. 1955, 23, 1762-1766.
(87) Chu, D.-Y.; Zhang, Q.; Liu, R.-L. Standard Gibbs free energies of transfer of NaCl and KCl from water to mixtures of the four isomers of butyl alcohol with water. The use of ion-selective electrodes to study the thermodynamics of solutions. J. Chem. Soc., Faraday Trans. 1 1987, 83, 635.
(88) Johari, G. P. Dielectric constants, densities, and viscosities of acetone-1-propanol and acetone-n-hexane mixtures at 25 . deg. C. J. Chem. Eng. Data 1968, 13, 541-543.
![](https://cdn.mathpix.com/cropped/790aed8a-6454-4c32-b8d0-0b109e9556ff-13.jpg?height=1135&width=223&top_left_y=1499&top_left_x=1090)

## CAS BIOFINDER DISCOVERY PLATFORM ${ }^{\text {TM }}$

## STOP DIGGING THROUGH DATA -START MAKING DISCOVERIES

CAS BioFinder helps you find the right biological insights in seconds

## Start your search


[^0]:    Special Issue: In Honor of Joan F. Brennecke
    Received: November 14, 2021
    Accepted: April 15, 2022
    Published: May 2, 2022

[^1]:    ${ }^{a}$ The equilibrium composition is given in the text.

