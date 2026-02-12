# Perturbed-Chain SAFT: An Equation of State Based on a Perturbation Theory for Chain Molecules 

J oachim Gross and Gabriele Sadowski*<br>Lehrstuhl für Thermodynamik, Universität Dortmund, Emil-FiggeStrasse 70, 44227 Dortmund, Germany

Downloaded via BRIGHAM YOUNG UNIV on May 14, 2024 at 14:27:31 (UTC). See https://pubs.acs.org/sharingguidelines for options on how to legitimately share published articles.


#### Abstract

A modified SAFT equation of state is developed by applying the perturbation theory of Barker and Henderson to a hard-chain reference fluid. With conventional one-fluid mixing rules, the equation of state is applicabl e to mixtures of small spherical molecules such as gases, nonspherical solvents, and chainlike polymers. The three pure-component parameters required for nonassociating molecules were identified for 78 substances by correlating vapor pressures and liquid volumes. The equation of state gives good fits to these properties and agrees well with caloric properties. When applied to vapor-liquid equilibria of mixtures, the equation of state shows substantial predictive capabilities and good precision for correlating mixtures. Comparisons to the SAFT version of Huang and Radosz reveal a clear improvement of the proposed model. A brief comparison with the Peng-Robinson model is also given for vapor-liquid equilibria of binary systems, confirming the good performance of the suggested equation of state. The applicability of the proposed model to polymer systems was demonstrated for high-pressure liquid-liquid equilibria of a polyethylene mixture. The pure-component parameters of polyethylene were obtained by extrapolating pure-component parameters of the n-alkane series to high molecular weights.


## 1. Introduction

The prediction or correlation of thermodynamic properties and phase equilibria with equations of state remains an important goal in chemical and related industries. Although the use of equations of state has for a long time been restricted to systems of simple fluids, there is an increasing demand for models that are also suitable for complex and macromolecular compounds. Clearly, the most apparent progress toward equations of state with such capabilities was made by applying principles of statistical mechanics. Some early models derived from statistical mechanics assumed molecules to be arranged in a lattice, ${ }^{1}$ whereas many of the more recent theories picture molecules to be moving freely in continuous space. A detailed review of different lines of development is given by Wei and Sadus. ${ }^{2}$

During the past few years, many studies assumed nonspherical molecules to be chains of freely jointed spherical segments. Despite its simplicity, this molecular model accounts for size and shape effects of molecules and has successfully been applied to simple species as well as large polymeric fluids and their mixtures. The first widely applied equation of state based on this molecular view was the perturbed hard-chain theory (PHCT) equation of state developed by Beret and Prausnitz ${ }^{3}$ and Donohue and Prausnitz. ${ }^{4}$ Their work revealed the potential of molecular-based theories and has been the inspiration for further developments.

A more recent equation-of-state concept for chain molecules is based on Wertheim's thermodynamic perturbation theory of first order. ${ }^{5-8}$ By applying Wertheim's theory and extending it to mixtures, Chapman et al. ${ }^{9,10}$ derived the statistical associating fluid theory

[^0](SAFT) equation of state for chain mixtures. Many modifications of the SAFT model were suggested over the years, examples being LJ -SAFT versions, ${ }^{11-17}$ in which Lennard-J ones spheres served as a reference for the chain formation, and VR-SAFT, in which the attractive potentials are allowed to show variable widths. ${ }^{18}$ Despite many theoretical improvements, one of the most successful modification remains the SAFT model suggested by Huang and Radosz, ${ }^{19,20}$ who applied a dispersion term developed by Chen and Kreglewski ${ }^{21}$ in the framework of SAFT. This dispersion term was derived by fitting a perturbation expansion to the experimental data of argon. The nonspherical shape of molecules is not accounted for in their dispersion term.

In a previous study, ${ }^{22}$ the authors were concerned with developing a theory for chain molecules, applying the perturbation theory of Barker and Henderson ${ }^{23,24}$ to a hard-chain reference. This theory was compared to simulation data for square-well chains and was found to give good results.

In this work, we will derive a dispersion expression for chain molecules by applying a perturbation theory for chain molecules and adjusting the appropriate model constants to the pure-component properties of $n$-alkanes. The equation of state developed here uses the same chain term and association term as the earlier SAFT equations. Because a hard-chain fluid serves as a reference for the perturbation theory in the present work, rather than spherical molecules as in the former SAFT modifications, the proposed model is referred to as perturbed-chain SAFT (PC-SAFT).

Because a new dispersion term is derived here, this work is focused on nonassociating components, for which the total attraction is dominated by dispersive forces. The SAFT version of Huang and Radosz (hereafter SAFT for brevity) is certainly the most widely applied version in today's industries. To eval uate the PC-SAFT model, comparisons are made to this SAFT model.

## 2. Theory

In a previous work, an equation of state for squarewell chain fluids was derived. ${ }^{22}$ This theory will now be extended to real substances. The procedure leading to the equation of state for real substances is similar in spirit to the work of Chen and K reglewski. ${ }^{21}$ After Alder et al. ${ }^{25}$ developed an expression for square-well spheres, Chen and Kreglewski extended the theory to describe real fluids of approximately spherical shape. In this work, we proceed similarly: starting from a theory for square-well chain molecules, we will obtain a model for real chain molecules of any length, from spheres to pol ymers.
2.1. Molecular Model. I n the proposed equation of state, molecules are conceived to be chains composed of spherical segments. The pair potential for the segment of a chain is given by a modified square-well potential, which was suggested by Chen and Kreglewski ${ }^{21}$

$$
u(r)= \begin{cases}\infty & r<\left(\sigma-\mathrm{s}_{1}\right)  \tag{1}\\ 3 \epsilon & \left(\sigma-\mathrm{s}_{1}\right) \leq \mathrm{r}<\sigma \\ -\epsilon & \sigma \leq \mathrm{r}<\lambda \sigma \\ 0 & r \geq \lambda \sigma\end{cases}
$$

where $u(r)$ is the pair potential, $r$ is the radial distance between two segments, $\sigma$ is the temperature-independent segment diameter, $\epsilon$ denotes the depth of the potential well, and $\lambda$ is the reduced well width. As suggested by Chen and Kreglewski, a ratio of $\mathrm{S}_{\mathrm{I}} / \sigma=0.12$ is assumed. In contrast to the work of Chen and Kreglewski, no additional temperature correction for the potential depth is introduced. Accordingly, nonassociating molecules are characterized by three purecomponent parameters: the temperature-independent segment diameter $\sigma$, the depth of the potential $\epsilon$, and the number of segments per chain m .

Although this potential model is very simple, the step function in the pair potential at $\mathrm{r}<\sigma$ accounts for an essential feature of real molecule behavior, the soft repulsion. The soft repulsion is introduced, because molecules have a collision diameter of $\sigma$ only when they collide at infinitely slow speed (zero temperature limit). Increasing temperature will result in a lower collision diameter.
2.2. Equation of State. According to perturbation theories, the interactions of molecules can be divided into a repulsive part and a contribution due to the attractive part of the potential. To calculate the repulsive contribution, a reference fluid in which no attractions are present is defined. The attractive interactions are treated as a perturbation to the reference system.

In the framework of Barker and Henderson's perturbation theory, ${ }^{24}$ a reference fluid with hard repulsion and a temperature-dependent segment diameter d(T) can be used to describe the soft repulsion of molecules, where

$$
\begin{equation*}
\mathrm{d}(\mathrm{~T})=\int_{0}^{\sigma}\left[1-\exp \left(-\frac{\mathrm{u}(\mathrm{r})}{\mathrm{kT}}\right)\right] \mathrm{dr} \tag{2}
\end{equation*}
$$

The reference fluid is given here by the hard-chain fluid and $\mathrm{d}(\mathrm{T})$ is the effective collision diameter of the chain segments. For the potential function given in eq 1, integration leads to the temperature-dependent hard segment diameter $\mathrm{d}_{\mathrm{i}}(\mathrm{T})$ of component i , according to

$$
\begin{equation*}
\mathrm{d}_{\mathrm{i}}(\mathrm{~T})=\sigma_{\mathrm{i}}\left[1-0.12 \exp \left(-\frac{3 \epsilon_{\mathrm{i}}}{\mathrm{kT}}\right)\right] \tag{3}
\end{equation*}
$$

The complete equation of state is given as an ideal gas contribution (id), a hard-chain contribution (hc), and a perturbation contribution, which accounts for the attractive interactions (disp).

$$
\begin{equation*}
Z=Z^{i d}+Z^{h c}+Z^{d i s p} \tag{4}
\end{equation*}
$$

where $Z$ is the compressibility factor, with $Z=P v /(R T)$ and $\mathrm{Z}^{\mathrm{id}}=1 ; \mathrm{P}$ is the pressure; v is the molar volume; and R denotes the gas constant. At this point, only dispersive attractions are considered. Specific interactions, such as hydrogen bonding or multipole interactions can be treated separately ${ }^{9,15}$ and will be considered in a subsequent investigation.

Hard-Chain Reference Equation of State. Based on Wertheim's ${ }^{5-8}$ thermodynamic perturbation theory of first-order, Chapman et al. ${ }^{9,10}$ developed an equation of state, which, for homonuclear hard-sphere chains comprising $m$ segments, is given by

$$
\begin{gather*}
\mathrm{z}^{\mathrm{hc}}=\overline{\mathrm{m}} \mathrm{z}^{\mathrm{hs}}-\sum_{\mathrm{i}} \mathrm{x}_{\mathrm{i}}\left(\mathrm{~m}_{\mathrm{i}}-1\right) \rho \frac{\partial \mathrm{In} \mathrm{~g}_{\mathrm{ii}}^{\mathrm{hs}}}{\partial \rho}  \tag{5}\\
\overline{\mathrm{~m}}=\sum_{\mathrm{i}} \mathrm{x}_{\mathrm{i}} \mathrm{~m}_{\mathrm{i}} \tag{6}
\end{gather*}
$$

where $x_{i}$ is the mole fraction of chains of component $i$, $\mathrm{m}_{\mathrm{i}}$ is the number of segments in a chain of component $\mathrm{i}, \rho$ is the total number density of molecules, $\mathrm{g}_{\mathrm{ii}}^{\text {hs }}$ is the radial pair distribution function for segments of component $i$ in the hard sphere system, and the superscript hs indicates quantities of the hard-sphere system. Expressions of Boublik ${ }^{26}$ and Mansoori et al. ${ }^{27}$ are used for mixtures of the hard-sphere reference system in eq 5 , given by

$$
\begin{array}{r}
Z^{h c}=\frac{\zeta_{3}}{\left(1-\zeta_{3}\right)}+\frac{3 \zeta_{1} \zeta_{2}}{\zeta_{0}\left(1-\zeta_{3}\right)^{2}}+\frac{3 \zeta_{2}^{3}-\zeta_{3} \zeta_{2}^{3}}{\zeta_{0}\left(1-\zeta_{3}\right)^{3}} \\
g_{i j}^{h s}=\frac{1}{\left(1-\zeta_{3}\right)}+\left(\frac{d_{i} d_{j}}{d_{i}+d_{j}}\right) \frac{3 \zeta_{2}}{\left(1-\zeta_{3}\right)}+ \\
\left(\frac{d_{i} d_{j}}{d_{i}+d_{j}}\right)^{2} \frac{2 \zeta_{2}^{2}}{\left(1-\zeta_{3}\right)^{3}} \tag{8}
\end{array}
$$

where

$$
\begin{equation*}
\xi_{n}=\frac{\pi}{6} \rho \sum x_{i} m_{i} d_{i}^{n} \quad n \in\{0,1,2,3\} \tag{9}
\end{equation*}
$$

Note that the compressibility factors of both the hardchain fluid in eq 5 and the hard-sphere fluid in eq 7 are residual properties, whereas they are often written including the ideal gas contribution in other literature (e.g., our previous work ${ }^{22}$ ).

Perturbation Theory for Pure Chain Molecules. After the reference chain fluid has been defined (it is identical to the SAFT reference fluid), the perturbation theory of Barker and Henderson can be used to calculate the attractive part of the chain interactions. It is a theory of second order, where the Helmholtz free energy is given as a sum of first- and second-order contributions via

$$
\begin{equation*}
\frac{\mathrm{A}^{\text {disp }}}{\mathrm{kTN}}=\frac{\mathrm{A}_{1}}{\mathrm{kTN}}+\frac{\mathrm{A}_{2}}{\mathrm{kTN}} \tag{10}
\end{equation*}
$$

Barker and Henderson derived their theory for spherical molecules. This theory can be extended to chain molecules, as each segment of a considered chain is again of spherical shape. The total interaction between two chain molecules required in the perturbation theory is then given by the sum of all individual segmentsegment interactions. ${ }^{22}$ Chiew ${ }^{28}$ obtained expressions for the individual segment-segment radial distribution function $\mathrm{g}_{\alpha \beta}^{\mathrm{hc}}\left(\mathrm{m} ; \mathrm{r}_{\alpha \beta}, \rho\right)$, which represents the radial distribution function for a segment $\alpha$ of one chain and a segment $\beta$ of another chain separated by the radial distance $\mathrm{r}_{\alpha \beta}$. Chiew also introduced an average interchain segment-segment radial distribution function $\mathrm{g}^{\text {hc_ }}$ ( $m ; r, \rho$ ), where different segments in a chain are nondistinguishable. It is convenient to determine the total interaction between two chains by applying this average radial distribution function. Gross and Sadowski ${ }^{22}$ used the results of Chiew and tested the theory for squarewell chains. The appropriate equations can easily be written for any potential function as

$$
\begin{array}{r}
\frac{\mathrm{A}_{1}}{\mathrm{kTN}}=-2 \pi \rho \mathrm{~m}^{2}\left(\frac{\epsilon}{\mathrm{kT}}\right) \sigma^{3} \int_{1}^{\infty} \tilde{\mathrm{u}}(\mathrm{x}) \mathrm{g}^{\mathrm{hc}}\left(\mathrm{~m} ; \mathrm{x} \frac{\sigma}{\mathrm{~d}}\right) \mathrm{x}^{2} \mathrm{dx} \\
\frac{\mathrm{~A}_{2}}{\mathrm{kTN}}=-\pi \rho \mathrm{m}\left(1+\mathrm{Z}^{\mathrm{hc}}+\rho \frac{\partial \mathrm{Z}^{\mathrm{hc}}}{\partial \rho}\right)^{-1} \mathrm{~m}^{2}\left(\frac{\epsilon}{\mathrm{kT}}\right)^{2} \times \\
\sigma^{3} \frac{\partial}{\partial \rho}\left[\rho \int_{1}^{\infty} \tilde{\mathrm{u}}(\mathrm{x})^{2} \mathrm{~g}^{\mathrm{hc}}\left(\mathrm{~m} ; \mathrm{x} \frac{\sigma}{\mathrm{~d}}\right) \mathrm{x}^{2} \mathrm{dx}\right] \tag{12}
\end{array}
$$

where $x$ is the reduced radial distance around a segment $(\mathrm{x}=\mathrm{r} / \sigma), \mathrm{u}(\mathrm{x})=\mathrm{u}(\mathrm{x}) / \epsilon$ denotes the reduced potential function, and $\mathrm{g}^{\mathrm{hc}}(\mathrm{m} ; x \sigma / \mathrm{d})$ is the average segmentsegment radial distribution function of the hard-chain fluid with temperature-dependent segment diameter $\mathrm{d}(\mathrm{T})$. The compressibility term in eq 12 can be obtained from eq 5 in the form

$$
\begin{align*}
& \left(1+\mathrm{Z}^{\mathrm{hc}}+\rho \frac{\partial \mathrm{Z}^{\mathrm{hc}}}{\partial \rho}\right)= \\
& \left(1+\mathrm{m} \frac{8 \eta-2 \eta^{2}}{(1-\eta)^{4}}+(1-\mathrm{m}) \frac{20 \eta-27 \eta^{2}+12 \eta^{3}-2 \eta^{4}}{[(1-\eta)(2-\eta)]^{2}}\right) \tag{13}
\end{align*}
$$

where the packing fraction $\eta$ is equal to $\zeta_{3}$ defined in eq 9. The packing fraction $\eta$ represents a reduced segment density.

Expressions for the radial distribution function of the hard-chain system are available in analytic form; ${ }^{29,30}$ however, these expressions are lengthy and lead to tedious calculations here, as an integration over $g^{h c}(r)$ is required in eqs 11 and 12. It is desirable therefore, to simplify the equation of state, and to do so, we first introduce the following abbreviations for the integrals in eq 11 and 12:

$$
\begin{gather*}
\mathrm{I}_{1}=\int_{1}^{\infty} \tilde{\mathrm{u}}(\mathrm{x}) \mathrm{g}^{\mathrm{hc}}\left(\mathrm{~m} ; \mathrm{x} \frac{\sigma}{\mathrm{~d}}\right) \mathrm{x}^{2} \mathrm{dx}  \tag{14}\\
\mathrm{I}_{2}=\frac{\partial}{\partial \rho}\left[\rho \int_{1}^{\infty} \tilde{\mathrm{u}}(\mathrm{x})^{2} \mathrm{~g}^{\mathrm{hc}}\left(\mathrm{~m} ; \mathrm{x} \frac{\sigma}{\mathrm{~d}}\right) \mathrm{x}^{2} \mathrm{dx}\right] \tag{15}
\end{gather*}
$$

For square-well chains, those integrals are functions of density and segment number only. For molecules
exhibiting soft repulsion, $\mathrm{I}_{1}$ and $\mathrm{I}_{2}$ are also functions of temperature. However, the temperature dependence due to $\mathrm{g}^{\mathrm{hc}}(\mathrm{m} ; x \sigma / \mathrm{d})$ is moderate and will be neglected here. With this assumption, it is possible to substitute the integrals $\mathrm{I}_{1}$ and $\mathrm{I}_{2}$ by power series in density $\eta$, where the coefficients of the power series are functions of the chain length.

$$
\begin{align*}
& \mathrm{I}_{1}(\eta, \mathrm{~m})=\sum_{\mathrm{i}=0}^{6} \mathrm{a}_{\mathrm{i}}(\mathrm{~m}) \eta^{\mathrm{i}}  \tag{16}\\
& \mathrm{I}_{2}(\eta, \mathrm{~m})=\sum_{\mathrm{i}=0}^{6} \mathrm{~b}_{\mathrm{i}}(\mathrm{~m}) \eta^{\mathrm{i}} \tag{17}
\end{align*}
$$

It was shown earlier ${ }^{22}$ that an expression proposed by Liu and $\mathrm{Hu}^{31}$ captures the dependency of coefficients $\mathrm{a}_{\mathrm{i}}(\mathrm{m})$ and $\mathrm{b}_{\mathrm{i}}(\mathrm{m})$ upon segment number accurately. It is given by

$$
\begin{align*}
& a_{i}(m)=a_{0 i}+\frac{m-1}{m} a_{1 i}+\frac{m-1}{m} \frac{m-2}{m} a_{2 i}  \tag{18}\\
& b_{i}(m)=b_{0 i}+\frac{m-1}{m} b_{1 i}+\frac{m-1}{m} \frac{m-2}{m} b_{2 i} \tag{19}
\end{align*}
$$

These equations were derived from a perturbation theory (sticky-point model based on Cummings and Stell ${ }^{32,33}$ ) assuming a correlation of nearest-neighbor segments and next-nearest neighbors. Equations 18 and 19 thus account for the bonding of one segment to a nearest-neighbor segment and for the possible bonding of the neighbor segment to a next-nearest-neighbor segment.

Determining Mode Constants. Let us now be concerned with identifying the model constants $a_{0 i}, a_{1 i}$, and $a_{2 i}$ as well as $b_{0 i}, b_{1 i}$, and $b_{2 i}$ of eqs 18 and 19. In a previous work, ${ }^{22}$ these constants were obtained by fitting the eqs 16 and 17 to eqs 14 and 15 for a squarewell potential using the radial distribution function proposed by Chiew. ${ }^{28}$ The appropriate model constants were universal, as the enti re ranges for the parameters m and $\eta$ were covered, i.e., m varies between $\mathrm{m}=1$ for spherical molecules and $\mathrm{m} \rightarrow \infty$ for infinitely long chains, and the packing fraction ranges between 0 for an ideal gas and $\eta \leq 0.74$ for the closest packing of segments.

This procedure is, in principle, possible for the potential function given above (eq 1); however, it has proven to be of advantage for an equation of state to incorporate information of real substance behavior. Many of the most successful models derived from statistical mechanics adjusted model constants to purecomponent data of real substances. Pure-component data of argon was, for example, used to adjust model constants in Chen and Kreglewski's BACK equation of state. ${ }^{21}$ The same dispersion term was also used in the SAFT model of Huang and Radosz. In the case of PHCT, argon and methane served as the model substances ${ }^{3}$. The model constants were in all cases considered to be universal. The reason that this procedure leads to superior models is three-fold: First, there are uncertainties in the dispersion properties, namely, in the assumed perturbing potential $\mathrm{u}(\mathrm{x})$ as well as approximations in $g^{h c}(r)$. Second, errors introduced in the reference equation of state can be corrected to a certain extent. Last, the molecular model assuming molecules to be chains of spherical segments might be oversimpli-

Table 1. Universal Model Constants for Equations 18 and 19
| i | $\mathrm{a}_{0 \mathrm{i}}$ | $\mathrm{a}_{1 \mathrm{i}}$ | $\mathrm{a}_{2 \mathrm{i}}$ | $\mathrm{b}_{0 \mathrm{i}}$ | $\mathrm{b}_{1 \mathrm{i}}$ | $\mathrm{b}_{2 \mathrm{i}}$ |
| :--- | :--- | :--- | :--- | :--- | :--- | :--- |
| 0 | 0.9105631445 | -0.3084016918 | -0.0906148351 | 0.7240946941 | -0.5755498075 | 0.0976883116 |
| 1 | 0.6361281449 | 0.1860531159 | 0.4527842806 | 2.2382791861 | 0.6995095521 | -0.2557574982 |
| 2 | 2.6861347891 | -2.5030047259 | 0.5962700728 | -4.0025849485 | 3.8925673390 | -9.1558561530 |
| 3 | -26.547362491 | 21.419793629 | -1.7241829131 | -21.003576815 | -17.215471648 | 20.642075974 |
| 4 | 97.759208784 | -65.255885330 | -4.1302112531 | 26.855641363 | 192.67226447 | -38.804430052 |
| 5 | -159.59154087 | 83.318680481 | 13.776631870 | 206.55133841 | -161.82646165 | 93.626774077 |
| 6 | 91.297774084 | -33.746922930 | -8.6728470368 | -355.60235612 | -165.20769346 | -29.666905585 |


fied. To correct for these shortcomings, we adjust the model constants $a_{0 i}, a_{1 i}, a_{2 i}, b_{0 i}, b_{1 i}$, and $b_{2 i}$ to experimental pure-component data. Because the proposed model accounts for the chainlike shape of molecules in the dispersion term, it is essential to include elongated molecules in the fitting procedure, and the series of n -alkanes is best suited to serve as model substances here. Methane can be assumed to be of spherical shape and will be used to determine the boundary case of $m =1$, where only the constants $a_{0 i}$ and $b_{0 i}$ in eqs 18 and 19 are relevant.

Our objective is to fit the power-series coefficients of the first-order term $a_{0 i}, a_{1 i}$, and $a_{2 i}$ as well as those of the second-order term $\mathrm{b}_{0 \mathrm{i}}, \mathrm{b}_{1 \mathrm{i}}$, and $\mathrm{b}_{2 \mathrm{i}}$ for $\mathrm{i}=1, \ldots, 6$ to pure-component data of $n$-alkanes. To obtain purecomponent parameters for the n -al kane components, an intermediate step has to be taken. We have assumed a Lennard-J ones perturbing potential in eqs 14 and 15 and used an expression for the average radial distribution function $\mathrm{g}^{\mathrm{hc}}(\mathrm{r})$ for hard chains proposed by Chiew. ${ }^{28}$ The integral expressions $\mathrm{I}_{1}$ and $\mathrm{I}_{2}$ (eqs 14 and 15) were so determined for Lennard-J ones-like chains. The three pure-component parameters ( $\mathrm{m}, \sigma, \epsilon / \mathrm{k}$ ) of the n-alkanes were identified for this equation of state by fitting vapor pressures and PvT data.

In a subsequent step, the coefficients $a_{0 i}, a_{1 i}, a_{2 i}, b_{0 i}$, $b_{1 i}$, and $b_{2 i}$ were regressed using the pure-component parameters determined before. Vapor pressures and liquid, vapor, and supercritical volumes were used in the regression, applying a Levenberg-Marquardt algorithm for minimizing the objective function

$$
\begin{equation*}
\operatorname{Min}=\sum_{\mathrm{i}=1}^{\mathrm{N}^{\exp }}\left|\frac{\boldsymbol{\Omega}_{\mathrm{i}}^{\exp }-\boldsymbol{\Omega}_{\mathrm{i}}^{\text {calc }}}{\boldsymbol{\Omega}_{\mathrm{i}}^{\exp }}\right|^{2} \tag{20}
\end{equation*}
$$

where $\Omega \in\left(\mathrm{P}^{\text {sat }}, \mathrm{v}\right)$ is the vapor pressure or the molar volume and $\mathrm{N}^{\text {exp }}$ is the total number of experimental points. The results for the coefficients $a_{0 i}, a_{1 i}$, and $a_{2 i}$ and $\mathrm{b}_{0 \mathrm{i}}, \mathrm{b}_{1 \mathrm{i}}$, and $\mathrm{b}_{2 \mathrm{i}}$ are given in Table 1. As for the dispersion expression of the BACK equation of state, these values are subsequently treated as universal model constants.

Mixtures. The perturbation theory of Barker and Henderson, as proposed here, makes use of an average radial distribution function and thus treats the segments of a chain as indistinguishable. Within this concept, a rigorous application of the perturbation theory to mixtures is, in principle, possible. O'Lenick et al. ${ }^{34}$ have derived a set of equations for the average radial pair distribution function of mixtures. Unfortunately, these expressions are not available in analytic form. However, the equation of state can easily be extended to mixtures by applying one-fluid mixing rules. Comparisons with simulation data of short-chain mixtures showed that the chain structure does not introduce any significant additional error to the one-fluid mixing
rule. ${ }^{17,22}$ Applying the van der Waals one-fluid mixing rules to the perturbation terms gives

$$
\begin{gather*}
\frac{\mathrm{A}_{1}}{\mathrm{kTN}}=-2 \pi \rho \mathrm{I}_{1}(\eta, \overline{\mathrm{~m}}) \sum_{\mathrm{i}} \sum_{\mathrm{j}} \mathrm{x}_{\mathrm{i}} \mathrm{x}_{\mathrm{j}} \mathrm{~m}_{\mathrm{i}} \mathrm{~m}_{\mathrm{j}}\left(\frac{\epsilon_{\mathrm{ij}}}{\mathrm{kT}}\right) \sigma_{\mathrm{ij}}^{3}  \tag{21}\\
\frac{\mathrm{~A}_{2}}{\mathrm{kTN}}=-\pi \rho \overline{\mathrm{m}}\left(1+\mathrm{Z}^{\mathrm{hc}}+\rho \frac{\partial \mathrm{Z}^{\mathrm{hc}}}{\partial \rho}\right)^{-1} \times \\
\mathrm{I}_{2}(\eta, \overline{\mathrm{~m}}) \sum_{\mathrm{i}} \sum_{\mathrm{j}} \mathrm{x}_{\mathrm{i}} \mathrm{x}_{\mathrm{j}} \mathrm{~m}_{\mathrm{i}} \mathrm{~m}_{\mathrm{j}}\left(\frac{\epsilon_{\mathrm{ij}}}{\mathrm{kT}}\right)^{2} \sigma_{\mathrm{ij}}^{3} \tag{22}
\end{gather*}
$$

where the power series $\mathrm{I}_{1}$ and $\mathrm{I}_{2}$ (eqs 14 and 15) can now be evaluated for the mean segment number $\bar{m}$ of the mixture, which was given by eq 6 .

The parameters for a pair of unlike segments are obtained by conventional Berthelot-Lorentz combining rules

$$
\begin{gather*}
\sigma_{\mathrm{ij}}=\frac{1}{2}\left(\sigma_{\mathrm{i}}+\sigma_{\mathrm{j}}\right)  \tag{23}\\
\epsilon_{\mathrm{ij}}=\sqrt{\epsilon_{\mathrm{i}} \epsilon_{\mathrm{j}}}\left(1-\mathrm{k}_{\mathrm{ij}}\right) \tag{24}
\end{gather*}
$$

where one binary interaction parameter, $\mathrm{k}_{\mathrm{ij}}$, is introduced to correct the segment-segment interactions of unlike chains. We also apply the one-fluid mixing concept to the compressibility term in eq 22, i.e., similarly to eq 13, it is

$$
\begin{align*}
& \left(1+\mathrm{Z}^{\mathrm{hc}}+\rho \frac{\partial \mathrm{Z}^{\mathrm{hc}}}{\partial \rho}\right)= \\
& \left(1+\overline{\mathrm{m}} \frac{8 \eta-2 \eta^{2}}{(1-\eta)^{4}}+(1-\overline{\mathrm{m}}) \frac{20 \eta-27 \eta^{2}+12 \eta^{3}-2 \eta^{4}}{[(1-\eta)(2-\eta)]^{2}}\right) \tag{25}
\end{align*}
$$

Equations for pressure, fugacity coefficients, and caloric properties can be derived from the Helmholtz free energy by applying classical thermodynamic relations. A summary of equations of the PC-SAFT model is given in Appendix A .

## 3. Results

3.1. Volumetric Data and Phase Equilibrium of Pure Components. The three pure-component parameters of the PC-SAFT model (segment number m , segment diameter $\sigma$, and segment energy parameter $\epsilon / \mathrm{k}$ ) can be regressed by fitting pure-component data. We applied the PC-SAFT model to nonassociating fluids including normal alkanes, branched alkanes, cyclic alkanes, alkenes, aromatics, chlorinated hydrocarbons, permanent gases, ethers, and esters. Data for the vapor pressure and liquid molar volumes were considered in all cases; additional PvT data was included in some

Table 2. Pure-Component Parameters for Nonassociating Substances
| substance | M [g/mol] | m [-] | $\sigma$ [Å] | $\epsilon / \mathrm{k}$ [K] | AAD\% |  | T range [K] | ref ${ }^{\text {a }}$ |  |
| :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- |
|  |  |  |  |  | psat | v |  | psat | v |
| n-alkanes |  |  |  |  |  |  |  |  |  |
| methane | 16.043 | 1.0000 | 3.7039 | 150.03 | 0.36 | $0.67{ }^{\text {b }}$ | 97-300 | 1 | 1,1 |
| ethane | 30.07 | 1.6069 | 3.5206 | 191.42 | 0.3 | 0.57 | 90-305 | 10 | 4 |
| propane | 44.096 | 2.0020 | 3.6184 | 208.11 | 1.29 | $0.77{ }^{\text {b }}$ | 85-523 | 2, 5 | 4, 6, 3 |
| butane | 58.123 | 2.3316 | 3.7086 | 222.88 | 0.75 | $1.59{ }^{\text {b }}$ | 135-573 | 2, 12 | 4, 7, 3 |
| pentane | 72.146 | 2.6896 | 3.7729 | 231.20 | 1.45 | 0.78 | 143-469 | 2, 5, 9 | 4 |
| hexane | 86.177 | 3.0576 | 3.7983 | 236.77 | 0.31 | 0.76 | 177-503 | 2, 5, 9 | 4 |
| heptane | 100.203 | 3.4831 | 3.8049 | 238.40 | 0.34 | $2.1{ }^{\text {b }}$ | 182-623 | 5, 7, 9 | 4, 7, 3 |
| octane | 114.231 | 3.8176 | 3.8373 | 242.78 | 0.77 | $1.59{ }^{\text {b }}$ | 216-569 | 5, 9 | 4, 6, 2 |
| nonane | 128.25 | 4.2079 | 3.8448 | 244.51 | 0.89 | 0.32 | 219-595 | 2, 5, 9 | 4 |
| decane | 142.285 | 4.6627 | 3.8384 | 243.87 | 0.24 | $1.18{ }^{\text {b }}$ | 243-617 | 9, 11 | 4, 12, 2 |
| undecane | 156.312 | 4.9082 | 3.8893 | 248.82 | 2.02 | 0.69 | 247-639 | 2, 5, 9 | 4 |
| dodecane | 170.338 | 5.3060 | 3.8959 | 249.21 | 2.1 | 0.93 | 263-658 | 2, 5, 9 | 4 |
| tridecane | 184.365 | 5.6877 | 3.9143 | 249.78 | 3.15 | 1.77 | 267-675 | 2, 5, 9 | 4 |
| tetradecane | 198.392 | 5.9002 | 3.9396 | 254.21 | 4.8 | 1.28 | 279-693 | 2, 5, 9 | 4, 5 |
| pentadecane | 212.419 | 6.2855 | 3.9531 | 254.14 | 1.04 | 5.35 | 283-708 | 2, 5, 9 | 4, 5 |
| hexadecane | 226.446 | 6.6485 | 3.9552 | 254.70 | 4.88 | 0.75 | 291-723 | 2, 5, 9 | 4, 5 |
| heptadecane | 240.473 | 6.9809 | 3.9675 | 255.65 | 5.35 | 0.69 | 295-736 | 2, 5, 9 | 4, 5 |
| octadecane | 254.5 | 7.3271 | 3.9668 | 256.20 | 4.99 | 1.19 | 301-747 | 2, 5, 9 | 4, 5 |
| nonadecane | 268.527 | 7.7175 | 3.9721 | 256.00 | 5.08 | 0.64 | 305-758 | 2, 5, 9 | 4, 5 |
| eicosane | 282.553 | 7.9849 | 3.9869 | 257.75 | 7.65 | $1.13{ }^{\text {b }}$ | 309-775 | 2, 9, 11 | 4, 5, 13 |
| branched alkanes |  |  |  |  |  |  |  |  |  |
| isobutane | 58.123 | 2.2616 | 3.7574 | 216.53 | 0.55 | 1.47 | 113-407 | 5 | 5 |
| isopentane | 72.15 | 2.5620 | 3.8296 | 230.75 | 0.4 | 1.53 | 113-460 | 5 | 5 |
| neopentane | 72.15 | 2.3543 | 3.9550 | 225.69 | 0.04 | 0.38 | 256-433 | 5 | 5 |
| 2-methylpentane | 86.177 | 2.9317 | 3.8535 | 235.58 | 0.61 | 0.59 | 119-498 | 5 | 5 |
| 2,2-dimethylbutane | 86.177 | 2.6008 | 4.0042 | 243.51 | 0.32 | 0.21 | 174-488 | 5 | 5 |
| 2,3-dimethylbutane | 86.177 | 2.6853 | 3.9545 | 246.07 | 0.38 | 0.46 | 145-500 | 5 | 5 |
| 3-methylpentane | 86.177 | 2.8852 | 3.8605 | 240.48 | 0.33 | 0.65 | 110-504 | 5 | 5 |
| 2-methylhexane | 100.204 | 3.3478 | 3.8612 | 237.42 | 0.74 | 1.44 | 154-530 | 5 | 5 |
| cycloalkanes |  |  |  |  |  |  |  |  |  |
| cyclopentane | 70.13 | 2.3655 | 3.7114 | 265.83 | 0.69 | 0.2 | 193-503 | 2 | 2 |
| cyclohexane | 84.147 | 2.5303 | 3.8499 | 278.11 | 0.53 | 3.12 | 279-553 | 5 | 5, 6, 12 |
| methylcyclopentane | 84.156 | 2.6130 | 3.8253 | 265.12 | 0.88 | 0.37 | 183-532 | 2 | 2 |
| methylcyclohexane | 98.182 | 2.6637 | 3.9993 | 282.33 | 1.91 | 0.31 | 203-572 | 2 | 2 |
| ethylcyclopentane | 98.182 | 2.9062 | 3.8873 | 270.50 | 2.55 | 0.75 | 134-569 | 2, 5 | 4, 5 |
| ethylcyclohexane | 112.215 | 2.8256 | 4.1039 | 294.04 | 2.38 | 1.04 | 263-609 | 2, 5 | 4, 5 |
| alkenes |  |  |  |  |  |  |  |  |  |
| ethylene | 28.05 | 1.5930 | 3.4450 | 176.47 | 1.16 | $2.61{ }^{\text {b }}$ | 104-400 | 16 | 15, 15 |
| propylene | 42.081 | 1.9597 | 3.5356 | 207.19 | 0.66 | 1.41 | 87-364 | 2, 5, 12 | 5, 12 |
| 1-butene | 56.107 | 2.2864 | 3.6431 | 222.00 | 0.69 | 0.52 | 87-419 | 2, 5, 12 | 4, 5, 6 |
| 1-pentene | 70.134 | 2.6006 | 3.7399 | 231.99 | 0.31 | 1.04 | 108-464 | 2, 5, 14 | 5 |
| 1-hexene | 84.616 | 2.9853 | 3.7753 | 236.81 | 0.42 | 1.23 | 133-504 | 2, 5, 14 | 5 |
| 1-octene | 112.215 | 3.7424 | 3.8133 | 243.02 | 0.79 | 0.75 | 171-567 | 5, 14 | 5 |
| cyclopentene | 68.114 | 2.2934 | 3.6668 | 267.76 | 0.14 | 0.07 | 223-393 | 2 | 2 |
| benzene derivatives |  |  |  |  |  |  |  |  |  |
| benzene | 78.114 | 2.4653 | 3.6478 | 287.35 | 0.64 | 1.42 | 278-562 | 5, 12 | 5, 12 |
| toluene | 92.141 | 2.8149 | 3.7169 | 285.69 | 2.41 | 1.35 | 178-594 | 2, 5, 12 | 5, 12 |
| ethylbenzene | 106.167 | 3.0799 | 3.7974 | 287.35 | 0.41 | 1.05 | 178-617 | 5 | 5 |
| m-xylene | 106.167 | 3.1861 | 3.7563 | 283.98 | 1.08 | 1.05 | 225-619 | 2, 5, 12 | 5, 12 |
| o-xylene | 106.167 | 3.1362 | 3.7600 | 291.05 | 0.64 | 1.15 | 248-630 | 2, 5 | 2, 4, 5 |
| p-xylene | 106.167 | 3.1723 | 3.7781 | 283.77 | 2.09 | 0.58 | 286-616 | 2, 5 | 2, 4, 5 |
| n-propyl benzene | 120.194 | 3.3438 | 3.8438 | 288.13 | 1.29 | 1.19 | 173-638 | 5 | 5 |
| tetralin | 132.205 | 3.3131 | 3.8750 | 325.07 | 1.06 | 0.54 | 237-720 | 5 | 5 |
| n-butylbenzene | 134.221 | 3.7662 | 3.8727 | 283.07 | 0.6 | 1.32 | 185-660 | 5 | 5 |
| biphenyl | 154.211 | 3.8877 | 3.8151 | 327.42 | 0.94 | 1.14 | 342-773 | 5 | 5 |
| gases |  |  |  |  |  |  |  |  |  |
| carbon monoxide | 28.01 | 1.3097 | 3.2507 | 92.15 | 2.21 | 1.38 | 80-133 | 5 | 5 |
| nitrogen | 28.01 | 1.2053 | 3.3130 | 90.96 | 0.34 | 1.5 | 63-126 | 17 | 17 |
| argon | 39.948 | 0.9285 | 3.4784 | 122.23 | 0.32 | 0.68 | 84-151 | 5 | 5 |
| carbon dioxide | 44.01 | 2.0729 | 2.7852 | 169.21 | 2.78 | $2.73{ }^{\text {b }}$ | 216-304 | 8 | 8, 8 |
| sulfur dioxide | 64.065 | 2.8611 | 2.6826 | 205.35 | 1.82 | 1.33 | 198-431 | 5 | 5 |
| chlorine | 70.905 | 1.5514 | 3.3672 | 265.67 | 0.99 | 0.57 | 172-417 | 5 | 5 |
| carbon disulfide | 76.143 | 1.6919 | 3.6172 | 334.82 | 2.26 | 0.77 | 200-552 | 5 | 5 |
| halogenated hydrocarbons |  |  |  |  |  |  |  |  |  |
| methyl chloride | 50.488 | 1.9297 | 3.2293 | 240.56 | 0.48 | 1.55 | 175-416 | 5 | 5 |
| chloroethane | 64.514 | 2.2638 | 3.4160 | 245.43 | 1.63 | 1.77 | 134-460 | 5 | 5 |
| 2-chloropropane | 78.541 | 2.4151 | 3.6184 | 251.47 | 2.89 | 0.61 | 156-489 | 5 | 5 |
| 1-chlorobutane | 92.568 | 2.8585 | 3.6424 | 258.66 | 0.99 | 1.14 | 150-537 | 5 | 5 |
| chlorobenzene | 112.558 | 2.6485 | 3.7533 | 315.04 | 2.59 | 0.66 | 228-632 | 2, 5 | 5 |
| bromobenzene | 157.01 | 2.6456 | 3.8360 | 334.37 | 1.27 | 1.44 | 242-670 | 5 | 5 |


Table 2 (Continued)
| substance | M [g/mol] | m [-] | $\sigma$ [Å] | $\epsilon / \mathrm{k}$ [K] | AAD\% |  | T range [K] | ref ${ }^{a}$ |  |
| :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- |
|  |  |  |  |  | psat | v |  | psat | v |
| ethers |  |  |  |  |  |  |  |  |  |
| dimethyl ether | 46.069 | 2.3071 | 3.2528 | 211.06 | 0.25 | 0.8 | 200-400 | 5, 12 | 5, 12 |
| diethyl ether | 74.123 | 2.9686 | 3.5147 | 220.09 | 0.46 | 0.94 | 240-467 | 5, 12 | 5, 12 |
| methyl-n-propyl ether | 74.123 | 3.0087 | 3.4569 | 222.73 | 1.99 | 1.04 | 220-476 | 5 | 5 |
| esters |  |  |  |  |  |  |  |  |  |
| methyl methanoate | 60.053 | 2.6784 | 3.0875 | 242.63 | 1.17 | 1.48 | 174-487 | 5 | 5 |
| ethyl methanoate | 74.079 | 2.8876 | 3.3109 | 246.47 | 0.92 | 1.8 | 193-508 | 5 | 5 |
| n-propyl methanoate | 88.106 | 3.2088 | 3.4168 | 246.46 | 0.95 | 1.01 | 270-538 | 5 | 5 |
| ethyl ethanoate | 88.106 | 3.5375 | 3.3079 | 230.80 | 1.3 | 2.24 | 189-523 | 5 | 5 |
| n-propyl ethanoate | 102.133 | 3.7861 | 3.4227 | 235.76 | 1.22 | 1.07 | 290-549 | 5 | 5 |
| n-butyl ethanoate | 116.16 | 3.9808 | 3.5427 | 242.52 | 2.03 | 1.26 | 199-579 | 5 | 5 |
| isopropyl ethanoate | 102.133 | 3.6084 | 3.4818 | 231.91 | 4.82 | 1.09 | 199-532 | 5 | 5 |
| methyl propanoate | 88.106 | 3.4793 | 3.3142 | 234.96 | 1.26 | 1.87 | 185-530 | 5 | 5 |
| ethyl propanoate | 102.133 | 3.8371 | 3.4031 | 232.78 | 2.18 | 1.07 | 199-546 | 5 | 5 |
| n-propyl propanoate | 116.16 | 4.1155 | 3.4875 | 235.60 | 3.53 | 0.73 | 220-568 | 5 | 5 |
| methyl butanoate | 102.133 | 3.6758 | 3.4437 | 240.62 | 1.36 | 0.94 | 187-554 | 5 | 5 |


${ }^{\mathrm{a}}$ References: (1) International Union of Pure and Applied Chemistry. International ThermodynamicTables of Fluid State5: Methane; Pergamon Press: Oxford, U.K., 1978. (2) Vargaftik, N. B. Tables of Thermophysical Properties of Liquids and Gases; J ohn Wiley \& Sons: New York, 1975. (3) Landolt-Börnstein, sechste Auflage, II. Band, 1. Teil: Mechanisch-Thermische Zustandgrössen; Springer: Berlin, 1971. (4) American Petroleum Institute Research Project 44. Selected Values of Properties of Hydrocarbons and Related Compounds; Texas A\&M University: College Station, Texas, 1973; Table 23-2-(1.101)-d. (5) Daubert, T. E.; Danner, R. P.; Sibul, H. M.; Stebbins, C. C. Physical and Thermodynamic Properties of PureChemicals: Data Compilation; Taylor \& Francis: Washington, D.C., 1989. (6) LandoltBörnstein, II. Band, 2. Teil Gleichgewichte ausser Schmelzgleichgewichten; Springer: Berlin, 1960. (7) Recommended Data of Selected Compounds and Binary Mixtures; DECHEMA Chemistry Data Series; DECHEMA: Frankfurt/Main, Germany, 1987; Vol. VI, Parts 1 and 2. (8) International Union of Pure and Applied Chemistry. International Thermodynamic Tables of Fluid State: Carbon Dioxide; Pergamon Press: Oxford, U.K., 1973. (9) Ruzicka, K.; Majer, V. J . Phys. Chem. Ref. Data 1994, 23 (1), 1-27. (10) Smith, B. D.; Srivastava, R. Thermodynamic Data for PureCompounds: Part A. Hydrocarbons and Ketones; Physical Sciences Data Series; Elsevier: Amsterdam, The Netherlands, 1986; Vol. 25. (11) American Petroleum Institute Research Project 44. Selected Values of Properties of Hydrocarbons and Related Compounds; Texas A\&M University: College Station, Texas, 1974; Table 23-2-(1.101)-kb. (12) VDI-Wärmeatlas, VDIGesellschaft Verfahrenstechnik und Chemieingenieurwesen (GVC), Düsseldorf, Germany, 1994. (13) Doolittle, A. K. J . Chem. Eng. Data 1964, 9 (2), 275-279. (14) Steele, W. V.; Chirico, R. D. J . Phys. Chem. Ref. Data 1993, 22 (2), 377-421. (15) Nowak, P.; Kleinrahm, R.; Wagner, W. J. Chem. Thermodyn. 1996, 28, 1423-1439. (16) Nowak, P.; Kleinrahm, R.; Wagner, W. J. Chem. Thermodyn. 1996, 28, 1441-1460. (17) International Union of Pure and Applied Chemistry. International Thermodynamic Tables of Fluid State: Nitrogen; Pergamon Press: Oxford, U.K., 1979. ${ }^{\text {b }}$ Supercritical density data was also used.
cases. Table 2 presents correlation results for 78 nonassociating compounds. Listed are the molar mass $\mathrm{M}_{\mathrm{i}}$, the segment number $\mathrm{m}_{\mathrm{i}}$, the segment diameter $\sigma_{\mathrm{i}}$, and the segment energy parameter $\epsilon_{\mathrm{i}} / \mathrm{k}$; furthermore, the temperature range covered by the experimental data and the absolute average deviation (AAD\%) of the vapor pressures and of the volume data are also considered. Considering the wide temperature range chosen for the regression, the average deviations in vapor pressure and liquid volumes are remarkably small.

It is an essential attribute for an EOS to accurately correlate vapor-liquid equilibria of pure components. Any systematic error observed in such calculations will usually propagate into mixture calculations. Figures 1A and 2A show the deviations between the PC-SAFT correlation results and the experimental data for some arbitrarily chosen representatives of different homologous series. Results of SAFT are presented in Figures 1 B and 2 B using the same scale of the diagram axis. It is not the purpose of Figures 1 and 2 to illustrate the deviations for some specific components but rather to determine systematic errors of the theories. The purecomponent parameters of the SAFT equation of state were, for this comparison, refitted to the same purecomponent data as for the PC-SAFT model. It should be noted that, because of the wide temperature range considered here, more pronounced errors of the SAFT equation of state were observed when the parameters derived by Huang and Radosz were used. From Figures 1 to 2 it becomes apparent that systematic deviations prevail for densities and for vapor pressures as calculated by the SAFT model, whereas a much improved
correlation of pure-component behavior is observed for the PC-SAFT equation of state. The improvement can be attributed to the dispersion contribution, which now accounts for the nonspherical shape of molecules.

Figure 3 shows the vapor-liquid equilibria of methane, propylene, diethyl ether, and toluene in a $\mathrm{T}-\rho$ diagram. The densities of the coexisting phases are welldescribed by the PC-SAFT equation of state. For methane, which is close to spherical in shape, results of SAFT are very similar to those of PC-SAFT. For nonspherical molecules, though, the PC-SAFT equation of state gives considerably better results for the liquid phase at low temperature and in the vicinity of the critical point. An improvement in these areas prevails not only for the arbitrarily chosen components displayed in Figure 3, but also for other pure nonpolar components.

Figure 4 shows the complete fluid-phase behavior of carbon dioxide in a pressure versus density plot. Displayed is the vapor-liquid region of carbon dioxide, from triple point to critical point, and subcritical and supercritical isotherms up to pressures of $10^{3}$ bar and temperatures of 450 K . The PC-SAFT equation of state describes the densities of the fluid region of carbon dioxide with good precision up to elevated pressures.

Well-behaved pure-component parameters are of importance in many petroleum industry applications, because they allow pure-component parameters to be correlated and, on demand, extrapolated to unknown or weakly characterized components. Pure-component parameters for n -alkanes are plotted in Figure 5. The parameters show a smooth course, approaching constant values for increasing molar mass. The pure-component

![](https://cdn.mathpix.com/cropped/84ec81da-9023-4737-8283-e6c19c5cd224-07.jpg?height=1351&width=837&top_left_y=158&top_left_x=186)
Figure 1. Relative deviation of theory from experiment in vapor pressure. The deviation from each experimental point is displayed with a small point; lines are for optical guidance. (A) PC-SAFT, (B) SAFT (given on the same scale).

parameters can be correlated using a functional form similar to eq 18, as displayed in Figure 5 (refer to Appendix B for details).
3.2. Caloric Data of Pure Components. The calculation of caloric properties requires temperature derivatives of the Helmholtz free energy. Although the hard-chain reference equation of state [with segment diameter $\mathrm{d}(\mathrm{T})$ ] has a contribution to these temperature derivatives, it is the dispersion term that plays a dominant role. Hence, applying an equation of state to caloric data is a meaningful test for the dispersion term. This is particularly the case as the caloric data were not included in the regression of pure-component parameters. Figure 6 compares experimental heats of vaporization to calculations from the PC-SAFT model and the SAFT equation of state. The PC-SAFT equation of state is capable of predicting the enthalpy of vaporization of nonassociating components.

Table 3 compares the two equations of state in calculations of heat capacities of pure liquids. Listed are calculation results for the PC-SAFT and SAFT models, the average deviations in both cases, the literature sources for $c_{p}$, and the ideal gas contribution of $c_{p}$. The agreement between PC-SAFT and the experimental data is very good for a predicted property and confirms a proper temperature behavior of the model.

![](https://cdn.mathpix.com/cropped/84ec81da-9023-4737-8283-e6c19c5cd224-07.jpg?height=1340&width=836&top_left_y=158&top_left_x=1102)
Figure 2. Relative deviation of theory from experiment in saturated liquid density. The deviation from each experimental point is displayed with a small point; lines are for optical guidance. (A) PC-SAFT, (B) SAFT (given on the same scale).

3.3. Mixtures. The PC-SAFT model is applicable to mixtures of gases, solvents, and polymers. The work presented here is focused on gas and solvent mixtures of nonassociating systems. Any binary pair of components can be corrected with one binary interaction parameter, $\mathrm{k}_{\mathrm{ij}}$, which was introduced in eq 24 . Engineers in chemical and petroleum industries often depend on the predictive capabilities of an equation of state. As a first step, the PC-SAFT model will be applied to binary mixtures with $\mathrm{k}_{\mathrm{ij}}$ set to 0 in order to explore its ability to predict the phase behavior of asymmetric systems. Figures $7-10$ show $\mathrm{P}-\mathrm{x}$ diagrams of methane-butane mixtures, ethane-decane mixtures, a diethyl etherethane mixture, and propane-benzene mixtures, respectively, at various temperatures. The PC-SAFT equation of state is in close agreement with the experimental data in all cases. When compared to SAFT, the PC-SAFT equation of state seems to be more predictive for liquid- and vapor-phase compositions and seems superior in the vicinity of the critical point. It is a general trend that the SAFT model predicts mixture critical points at too-high concentrations of the light (supercritical) components, whereas the PC-SAFT equation of state is seen to be more reliable here.

Results for a ternary system of heptane-butaneethane are displayed in Figure 11 for two pressures. As before, the phase diagram is predicted with all $\mathrm{k}_{\mathrm{ij}}$

![](https://cdn.mathpix.com/cropped/84ec81da-9023-4737-8283-e6c19c5cd224-08.jpg?height=875&width=843&top_left_y=163&top_left_x=184)
Figure 3. Saturated liquid and vapor densities for methane, propylene, diethyl ether, and toluene. Comparison of PC-SAFT (solid lines) and SAFT (dashed lines) to experimental data.

![](https://cdn.mathpix.com/cropped/84ec81da-9023-4737-8283-e6c19c5cd224-08.jpg?height=650&width=819&top_left_y=1164&top_left_x=191)
Figure 4. $\mathrm{P} \rho \mathrm{T}$ behavior of carbon dioxide. Subcritical and supercritical isotherms and coexisting densities from triple point to critical point from experiment (symbols) and PC-SAFT equation of state (solid lines).

parameters set to 0 . PC-SAFT describes the phase boundary as well as the coexisting concentrations (tie lines) with good accuracy.

Equations of state derived from statistical mechanics have successfully been established for systems of associating or macromolecular systems. Despite those achievements, one has to concede that cubic equations of state are still very powerful for the case of vaporliquid equilibria of volatile nonpolar substances. Han et at. ${ }^{35}$ compared seven cubic equations of state by fitting $\mathrm{k}_{\mathrm{ij}}$ parameters to vapor-liquid equilibria of 87 binary mixtures. Their study revealed the Peng-Robinson (PR) equation of state to be among the most reliable models. We will give a brief comparison of the PC-SAFT, SAFT, and PR equations of state by adjusting $\mathrm{k}_{\mathrm{ij}}$ parameters to vapor-liquid equilibria of 24 mixtures. The objective function for fitting the $\mathrm{k}_{\mathrm{ij}}$ values is chosen to be the same as described in the work of Han et al., ${ }^{35}$

![](https://cdn.mathpix.com/cropped/84ec81da-9023-4737-8283-e6c19c5cd224-08.jpg?height=1808&width=856&top_left_y=154&top_left_x=1093)
Figure 5. Pure-component parameters of the PC-SAFT equation of state for the series of $n$-alkanes as a function of molar mass. (A) $\epsilon_{i} / k$, segment energy parameter. (B) ( $m_{i} / M_{i}$ ), segment number parameter. (C) $\sigma_{\mathrm{i}}$, segment diameter.

namely minimization of the sum of the squared relative deviations, $\sum_{i}\left[\left(K_{i}^{\text {exp }}-K_{i}^{\text {calc }}\right) / K_{i}^{\text {exp }}\right]^{2}$, for the $K$ values of the two components in all experimental data pairs. Table 4 presents the results for these mixtures. The PCSAFT model is clearly superior to the SAFT equation of state for these systems. The absol ute average deviations in $\mathrm{K}_{1}$ and $\mathrm{K}_{2}$ values obtained from PC-SAFT are approximately one-half of those obtained from the SAFT results. It is also noteworthy that the $\mathrm{k}_{\mathrm{ij}}$ parameters from PC-SAFT are in average approximately one-half of the absolute values of the parameters from SAFT, which might also indicate higher predictive capabilities of the PC-SAFT equation of state. Table 4 also presents results for the PR model. Correlation results of the PC-

![](https://cdn.mathpix.com/cropped/84ec81da-9023-4737-8283-e6c19c5cd224-09.jpg?height=626&width=819&top_left_y=163&top_left_x=193)
Figure 6. Heat of vaporization for benzene (open circles), isobutane (solid circles), and propylene (diamonds). Comparison of PCSAFT (solid lines) and SAFT (dashed lines) to experimental data.

Table 3. Percent Absolute Average Deviation (AAD\%) in Heat Capacities of the PC-SAFT and SAFT Equation of State
|  | AAD\% |  | T range [ ${ }^{\circ} \mathrm{C}$ ] | no. of points | ref ${ }^{\text {a }}$ |  |
| :--- | :--- | :--- | :--- | :--- | :--- | :--- |
|  | PCSAFT |  |  |  |  |  |
|  |  | SAFT |  |  | $c_{p}$ | $c_{p}^{\text {ideal }}$ |
| isobutane | 0.96 | 8.65 | -50 to -25 | 2 | 1 | 2 |
| 2,3-dimethylbutane | 2.99 | 4.59 | -50 to 20 | 4 | 1 | 2 |
| propylene | 3.78 | 7.62 | -50 to 50 | 4 | 1 | 2 |
| 1-octene | 0.66 | 8.55 | -50 to 100 | 5 | 2 | 2 |
| cyclopentane | 0.51 | 12.45 | -50 to 20 | 4 | 1 | 2 |
| benzene | 1.75 | 1.87 | 50 to 200 | 4 | 3 | 2 |
| n-butylbenzene | 1.28 | 6.40 | -50 to 150 | 5 | 1 | 2 |
| overall average ${ }^{\text {b }}$ | 1.7 | 7.1 |  |  |  |  |


${ }^{\mathrm{a}}$ (1) VDI-Wärmeatlas, VDI-Gesellschaft Verfahrenstechnik und Chemieingenieurwesen (GVC), Düsseldorf, Germany, 1994. (2) Daubert, T. E.; Danner, R. P.; Sibul, H. M.; Stebbins, C. C. Physical and Thermodynamic Properties of PureChemicals: Data Compilation; Taylor \& Francis: Washington, D.C., 1989. (3) Vargaftik, N. B. Tables of Thermophysical Properties of Liquids and Gases; J ohn Wiley \& Sons: New York, 1975. ${ }^{\mathrm{b}}$ Average weight by the number of data points.

![](https://cdn.mathpix.com/cropped/84ec81da-9023-4737-8283-e6c19c5cd224-09.jpg?height=654&width=828&top_left_y=1708&top_left_x=193)
Figure 7. Vapor-liquid equilibrium of $n$-butane-methane mixtures at two temperatures. Comparison of PC-SAFT and SAFT predictions ( $\mathrm{k}_{\mathrm{ij}}=0$ in both cases) to experimental data (Sage et al., 1940).

SAFT equation of state are slightly better on average than those of the PR model. We consider this an important result, as it is our belief that particular

![](https://cdn.mathpix.com/cropped/84ec81da-9023-4737-8283-e6c19c5cd224-09.jpg?height=652&width=824&top_left_y=158&top_left_x=1112)
Figure 8. Vapor-liquid equilibrium of n -decane-ethane mixtures at two temperatures. Comparison of PC-SAFT and SAFT predictions ( $\mathrm{k}_{\mathrm{ij}}=0$ in both cases) to experimental data (Reamer and Sage, 1962).

![](https://cdn.mathpix.com/cropped/84ec81da-9023-4737-8283-e6c19c5cd224-09.jpg?height=617&width=832&top_left_y=956&top_left_x=1102)
Figure 9. Vapor-liquid isotherm of a diethyl ether-ethane mixture at T = $25^{\circ} \mathrm{C}$. Comparison of PC-SAFT (solid lines) and SAFT (dashed lines) predictions ( $\mathrm{k}_{\mathrm{ij}}=0$ in both cases) to experimental data.

![](https://cdn.mathpix.com/cropped/84ec81da-9023-4737-8283-e6c19c5cd224-09.jpg?height=630&width=826&top_left_y=1723&top_left_x=1106)
Figure 10. Vapor-liquid equilibrium of benzene-propane mixtures at three temperatures. Comparison of PC-SAFT and SAFT predictions ( $\mathrm{k}_{\mathrm{ij}}=0$ in both cases) to experimental data (Glanville et al., 1950).

strengths of semiempirical cubic equations of state lie in the vapor-liquid equilibria of systems such as those investigated here. When densities of liquids or generally

![](https://cdn.mathpix.com/cropped/84ec81da-9023-4737-8283-e6c19c5cd224-10.jpg?height=693&width=832&top_left_y=158&top_left_x=191)
Figure 11. Vapor-liquid equilibrium of ternary n-heptane-n-butane-ethane mixtures at $\mathrm{T}=121^{\circ} \mathrm{C}$ and for pressures of $\mathrm{P}=$ 69 bar (circles) and $\mathrm{P}=48$ bar (triangles). Comparison of PCSAFT predictions (all $\mathrm{k}_{\mathrm{ij}}=0$ ) to experimental data (Mehra and Thodos, 1966).

more complex systems are considered, it is advantageous to use models with a strong physical background, and equations of state derived from statistical mechanics probably provide the highest degree of potential here.

The PC-SAFT equation of state will now be applied to more nonideal systems. Although carbon dioxide has

![](https://cdn.mathpix.com/cropped/84ec81da-9023-4737-8283-e6c19c5cd224-10.jpg?height=645&width=828&top_left_y=165&top_left_x=1108)
Figure 12. Vapor-liquid equilibrium of n-decane- $\mathrm{CO}_{2}$ mixtures at two temperatures. Comparison of PC-SAFT ( $\mathrm{k}_{\mathrm{ij}}=0.128$ ) and SAFT ( $\mathrm{k}_{\mathrm{ij}}=0.115$ ) correlation results to experimental data (Reamer and Sage, 1963).

a quadrupole moment, it is here assumed to exhibit dispersive attractions only. A $\mathrm{k}_{\mathrm{ij}}$ parameter is required for binary $\mathrm{CO}_{2}$ systems. Figures 12 and 13 display binary mixtures of $\mathrm{CO}_{2}$-decane and $\mathrm{CO}_{2}$-cyclohexane, respectively, at two temperatures. The PC-SAFT equation of state gives good correlation results for both systems. Figure 12 suggests that improvements of the PC-SAFT model observed for simple fluids in Figures 7

Table 4. Correlation Results for Binary Vapor-Liquid Equilibria of Simple Fluids with the PC-SAFT Equation, the SAFT Model, and the PR Equation of State
| system | PC-SAFT |  |  | SAFT |  |  | PR |  |  | T range $\left[{ }^{\circ} \mathrm{C}\right]$ | ref ${ }^{\text {a }}$ |
| :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- |
|  | AAD\% |  |  | AAD\% |  |  | AAD\% |  |  |  |  |
|  | $\mathrm{k}_{\mathrm{ij}}$ | $\mathrm{K}_{1}$ | $\mathrm{K}_{2}$ | $\mathrm{k}_{\mathrm{ij}}$ | $\mathrm{K}_{1}$ | $\mathrm{K}_{2}$ | $\mathrm{k}_{\mathrm{ij}}$ | $\mathrm{K}_{1}$ | $\mathrm{K}_{2}$ |  |  |
| methane-butane | 0.022 | 1.64 | 2.60 | 0.074 | 5.62 | 2.84 | 0.038 | 1.91 | 3.59 | 21-121 | 4 |
| methane-pentane | 0.024 | 1.56 | 1.65 | 0.077 | 6.89 | 4.39 | 0.033 | 1.67 | 1.67 | 71-171 | 5 |
| methane-hexane | 0.021 | 3.62 | 1.44 | 0.076 | 9.30 | 5.29 | 0.040 | 5.64 | 2.20 | 0.01-138 | 6, 7 |
| methane-heptane | 0.016 | 7.28 | 3.06 | 0.053 | 12.90 | 9.15 | 0.017 | 5.72 | 3.30 | 4-238 | 8 |
| methane-decane | 0.056 | 1.58 | 2.11 | 0.121 | 6.34 | 4.10 | 0.073 | 2.85 | 3.64 | 150-310 | 9 |
| methane-isobutane | 0.028 | 0.56 | 1.25 | 0.076 | 5.23 | 3.04 | 0.046 | 0.88 | 2.33 | 38-104 | 10 |
| methane-cyclohexane | 0.045 | 4.18 | 2.26 | 0.076 | 7.48 | 6.02 | 0.048 | 4.22 | 1.91 | 21-171 | 11 |
| methane-benzene | 0.037 | 5.19 | 2.85 | 0.091 | 9.93 | 4.47 | 0.064 | 5.75 | 3.85 | 148-228 | 9 |
| methane-toluene | 0.052 | 3.34 | 3.47 | 0.111 | 8.03 | 6.33 | 0.091 | 4.87 | 5.52 | 149-270 | 9 |
| methane-m-xylene | 0.045 | 1.95 | 1.21 | 0.122 | 6.49 | 4.52 | 0.092 | 3.59 | 3.96 | 188-309 | 12 |
| methane-tetralin | 0.075 | 3.95 | 2.58 | 0.125 | 7.09 | 5.04 | 0.140 | 6.52 | 4.22 | 189-391 | 13 |
| nitrogen-ethylene | 0.075 | 2.79 | 4.38 | 0.059 | 4.07 | 5.37 | 0.084 | 2.41 | 4.86 | -7 to -13 | 14 |
| nitrogen-hexane | 0.119 | 7.73 | 3.97 | 0.127 | 10.27 | 7.68 | 0.124 | 8.29 | 5.86 | 38-171 | 15 |
| ethylene-heptane | 0.018 | 4.33 | 1.45 | 0.051 | 11.79 | 2.10 | 0.019 | 4.43 | 2.37 | 60-243 | 16 |
| ethylene-propylene | -0.001 | 3.52 | 1.38 | 0.024 | 6.63 | 4.55 | -0.007 | 3.16 | 1.57 | 10-25 | 17 |
| propane-butane | 0.003 | 1.95 | 1.61 | 0.011 | 4.53 | 2.48 | -0.002 | 2.89 | 2.72 | 60-140 | 18 |
| propylene-isobutane | -0.006 | 3.41 | 3.25 | 0.007 | 4.92 | 2.68 | -0.019 | 5.31 | 4.56 | 42-124 | 19 |
| pentane-heptane | 0.011 | 2.95 | 2.35 | 0.017 | 4.17 | 2.65 | 0.014 | 3.78 | 1.72 | 131-249 | 20 |
| $\mathrm{CO}_{2}$-methane | 0.065 | 3.04 | 2.61 | 0.091 | 2.98 | 2.70 | 0.107 | 1.84 | 2.30 | -43 to -3 | 21 |
| $\mathrm{CO}_{2}$-propane | 0.109 | 3.92 | 5.45 | 0.108 | 4.80 | 6.02 | 0.128 | 3.90 | 4.79 | -33 to 57 | 22, 23 |
| $\mathrm{CO}_{2}$-butane | 0.120 | 5.91 | 5.56 | 0.116 | 6.72 | 6.19 | 0.124 | 8.52 | 7.62 | -45 to 145 | 1-3 |
| $\mathrm{CO}_{2}$-pentane | 0.143 | 6.99 | 5.04 | 0.149 | 11.18 | 6.24 | 0.141 | 11.59 | 3.28 | 4-104 | 24 |
| $\mathrm{CO}_{2}$-heptane | 0.129 | 8.68 | 2.65 | 0.129 | 16.85 | 2.00 | 0.107 | 9.90 | 1.87 | 37-204 | 25 |
| $\mathrm{CO}_{2}$-methylcyclohexane | 0.144 | 6.39 | 2.72 | 0.121 | 7.34 | 2.07 | 0.129 | 9.76 | 4.04 | 65-204 | 26 |
| average |  | 4.06 | 2.80 |  | 7.76 | 4.58 |  | 5.11 | 3.54 |  |  |


[^1]![](https://cdn.mathpix.com/cropped/84ec81da-9023-4737-8283-e6c19c5cd224-11.jpg?height=639&width=830&top_left_y=163&top_left_x=193)
Figure 13. Vapor-liquid equilibrium of cyclohexane- $\mathrm{CO}_{2}$ mixtures two temperatures. Comparison of PC-SAFT ( $\mathrm{k}_{\mathrm{ij}}=0.13$ ) and SAFT ( $\mathrm{k}_{\mathrm{ij}}=0.135$ ) correlation results to experimental data (Nagarajan and Robinson, 1987; Krichevskii and Sorina, 1960).

![](https://cdn.mathpix.com/cropped/84ec81da-9023-4737-8283-e6c19c5cd224-11.jpg?height=708&width=834&top_left_y=964&top_left_x=191)
Figure 14. Vapor-liquid equilibrium of a ternary $\mathrm{CO}_{2}$-ethanemethane mixture at $\mathrm{T}=-43.15^{\circ} \mathrm{C}$ and $\mathrm{P}=25.3$ bar. Comparison of PC-SAFT correlation results to experimental data (Wei et al., 1995).

and 8 also prevail for $\mathrm{CO}_{2}$-systems in a similar manner. This improvement seems to be of systematic nature. Figure 13 indicates that PC-SAFT is able to describe the phase equilibrium of the $\mathrm{CO}_{2}$-cyclohexane system over wide ranges of temperatures, whereas SAFT cannot be fitted well to both the temperatures with a constant $\mathrm{k}_{\mathrm{ij}}$ value.

Results for a ternary system of $\mathrm{CO}_{2}$-methaneethane are shown in Figure 14. The methane-ethane binary mixture behaves rather ideally, and one can safely assume a $\mathrm{k}_{\mathrm{ij}}$ equal to 0 . A $\mathrm{k}_{\mathrm{ij}}$ parameter for the $\mathrm{CO}_{2}$-methane binary system was given in Table 4. The $\mathrm{CO}_{2}$-ethane system, which shows an azeotropic point, was correl ated to experimental data at 250 and 270 K , resulting in a good fit with a binary interaction parameter of $\mathrm{k}_{\mathrm{ij}}=0.095$. Using only binary information, the phase boundaries in the ternary diagram are welldescribed, and the calculated tie lines match the experimental data.

Figure 15 shows the isobaric vapor-liquid equilibrium of a 1-chlorobutane-hexane mixture. The PC-SAFT equation of state can correlate the phase envel ope. The

![](https://cdn.mathpix.com/cropped/84ec81da-9023-4737-8283-e6c19c5cd224-11.jpg?height=639&width=822&top_left_y=163&top_left_x=1110)
Figure 15. Vapor-liquid isobar of a 1-chlorobutane-n-hexane mixture at $\mathrm{P}=0.944$ bar. Comparison of PC-SAFT ( $\mathrm{k}_{\mathrm{ij}}=0.017$ ) and SAFT ( $\mathrm{k}_{\mathrm{ij}}=0.016$ ) correlation results to experimental data.

![](https://cdn.mathpix.com/cropped/84ec81da-9023-4737-8283-e6c19c5cd224-11.jpg?height=620&width=824&top_left_y=930&top_left_x=1108)
Figure 16. Excess enthalpy of a 1-chlorobutane-n-hexane mixture at $\mathrm{T}=25^{\circ} \mathrm{C}$. Comparison of PC-SAFT ( $\mathrm{k}_{\mathrm{ij}}=0.017$ ) and SAFT ( $\mathrm{k}_{\mathrm{ij}}=0.016$ ) calculations to experimental data. Binary interaction parameters $\mathrm{k}_{\mathrm{ij}}$ were adjusted to vapor-liquid equilibrium data (see Figure 15).

SAFT model gives good results for the largest part of the binary system but calculates a slightly erroneous pure-component boiling point of hexane. The excess heat of mixing was subsequently calculated with the same $\mathrm{k}_{\mathrm{ij}}$ parameter; the results are displayed in Figure 16. The PC-SAFT model is in good agreement with the experimental data. Generally, an accurate excess heat of mixing indicates a reliable temperature behavior of the thermodynamic model.
3.4. Polymer Systems. Phase equilibria, and particularly liquid-liquid equilibria, of polymer systems are a challenging test for an equation of state. Only a cursory look at polymer systems will be taken here. A detailed investigation of polymer systems will be presented in a subsequent investigation.

Here, we will explore whether extrapolating the purecomponent parameters of the n -alkane series to molecular weights of polymers will result in a reasonable description for polyethylene. Correlation results for the pure-component parameters of small alkanes are given in Appendix B. Upon extrapolation to molecular weights of polymers, these correlations converge to values of $\sigma =4.072, \mathrm{~m} / \mathrm{M}=0.02434$, and $\epsilon / \mathrm{k}=269.67$ for polyethylene. The segment number is given as a ratio with

![](https://cdn.mathpix.com/cropped/84ec81da-9023-4737-8283-e6c19c5cd224-12.jpg?height=611&width=822&top_left_y=161&top_left_x=195)
Figure 17. Molar mass distribution of polyethylene as characterized by Kiran et al., 1994. ${ }^{36}$ Squares represent results from gelpermeation chromatography. Squares are used as pseudocomponents (representing the molar mass distribution of polyethylene) in a subsequent calculation.

![](https://cdn.mathpix.com/cropped/84ec81da-9023-4737-8283-e6c19c5cd224-12.jpg?height=682&width=830&top_left_y=964&top_left_x=191)
Figure 18. High-pressure phase equilibria of polyethylene-npentane at $186.85^{\circ} \mathrm{C}$. The molar mass distribution of the polyethylene is accounted for by assuming 10 pseudocomponents. Comparison of PC-SAFT (solid lines, $\mathrm{k}_{\mathrm{ij}}=-0.0195$ ) and SAFT (dashed lines, $\mathrm{k}_{\mathrm{ij}}=-0.003$ ) correlation results to experimental data.

the molar mass, as this is a convenient parameter for polymer systems, where the segment number is a linear function of the molar mass.

Kiran et al. ${ }^{36}$ have presented experimental liquidliquid equilibria measurements for a pentane-pol yethylene system of well-characterized polyethylene. The molecular weight distribution of the polyethylene was determined by gel-permeation chromatography, and as a result of this analysis, Kiran et al. have characterized the polyethylene by 11 components with specified molar masses and their appropriate weight fractions. These 11 polyethylene components represent the polydisperse polymer, and they were used as pseudocomponents in our calculation (Figure 17). The binary interaction parameter of the polyethylene-pentane system was adjusted to be $\mathrm{k}_{\mathrm{ij}}=-0.019$. The correlation result is given in Figure 18. The PC-SAFT results are in good agreement with the experimental data. It is particularly remarkable that the calculated demixing region is not as narrow as that of the SAFT equation of state, as this is a systematic deficiency of SAFT, ${ }^{37,38}$

An absolute average deviation of 6.6\% is obtained for the density data of pure polyethylene ${ }^{39}$ (the pressure range is from 1 to 1000 bar; the temperature range is from 153 to $200^{\circ} \mathrm{C}$ ). Because the pure-component parameters stem from extrapolations of the purecomponent parameters of $n$-alkanes, this is a very reasonable result.

Results of PC-SAFT for polyethylene, the purecomponent parameters of which were obtained by extrapolating pure-component parameters of Iow-mo-lecular-weight alkanes, are promising and reveal PCSAFT to be a highly robust model. One might argue, though, that, because of self-interactions and entanglement, the interactions in a polymer have to be different to a certain extent, from those of normal fluids. Obviously, pure-component parameters of polyethylene obtained from extrapolating the pure-component parameters of n -alkanes cannot account for these effects. In a subsequent investigation, we will therefore revise the pure-component parameters of polyethylene.

## 4. Conclusions

The perturbed-chain SAFT (PC-SAFT) equation of state presented in this work is based on the statistical associating fluid theory (SAFT) as a reference, applying Barker and Henderson's second-order perturbation theory. A hard-chain reference is considered in the perturbation theory, leading to a dispersion term that properly depends on the chain length of a molecule. The equation of state was simplified by substituting a power series in density for lengthy volume integrals of the Barker-Henderson theory. To incorporate information of real fluid behavior, the model constants of these power series were adjusted to vapor pressure and PvT data of the n-alkane series. This approach is in the spirit of earlier work by Chen and Kreglewski ${ }^{21}$ and Beret and Prausnitz, ${ }^{3}$ who chosen argon and methane as their model substances.

The PC-SAFT equation of state was applied to pure components and mixtures of nonassociating and weakly polar components. When compared to the SAFT version of Huang and Radosz, the PC-SAFT model was found to improve pure-component representation. Systematic errors, which occur when the dispersion term does not account for the nonspherical shape of molecules, were not found for the PC-SAFT equation of state. Caloric properties, which were not included in the parameter fitting, were calculated with good precision. A clear improvement over the earlier version of SAFT was also found for mixtures. The PC-SAFT model showed predictive strength and revealed precision in the correlation of asymmetric mixtures. The robustness of this equation of state was tested by extrapolating pure-component parameters of $n$-alkanes to high molecular weights of polymers, resulting in a set of polyethylene parameters. PC-SAFT proved to be in good agreement with experimental liquid-liquid data of a polyethylene-pentane mixture.

The goal of future work will be to demonstrate the applicability of the PC-SAFT equation of state to associating molecules by adopting the association expression of the SAFT framework (Chapman et al. ${ }^{9}$ ) and to apply this model to polymer systems.

## Acknowledgment

The authors are grateful to the Deutsche Forschungsgemeinschaft for supporting this work with Grant SAD 700/3. We especially thank Wolfgang Arlt, who encour-
aged us to focus on the improvement of equation-of-state modeling and who continuously accompanied the project with scientific interest and a high degree of commitment.

## List of Symbols

$\mathrm{A}=$ Helmholtz free energy, J
$\mathrm{A}_{1}=$ Helmholtz free energy of first-order perturbation term, J
$\mathrm{A}_{2}=$ Helmholtz free energy of second-order perturbation term, J
$\mathrm{a}_{01}, \mathrm{a}_{02}, \mathrm{a}_{03}=$ model constants; defined in eq 18
$a_{j}(m)=$ functions defined by eqs 18 and 19
$\mathrm{b}_{01}, \mathrm{~b}_{02}, \mathrm{~b}_{03}=$ model constants defined in eq 19
$\mathrm{d}=$ temperature-dependent segment diameter, Å
$\mathrm{g}^{\mathrm{hc}}=$ average radial distribution function of hard-chain fluid
$\mathrm{g}_{\alpha \beta}^{\text {hc }}=$ site-site radial distribution function of hard-chain fluid
$\mathrm{I}_{1}, \mathrm{I}_{2}=$ abbreviations defined by eqs 14-17
$\mathrm{k}=$ Boltzmann constant, J/K
$\mathrm{k}_{\mathrm{ij}}=$ binary interaction parameter
$\mathrm{K}=\mathrm{K}$ factor, $\mathrm{K}_{\mathrm{i}}=\mathrm{y}_{\mathrm{i}} / \mathrm{x}_{\mathrm{i}}$
$\mathrm{m}=$ number of segments per chain
$\overline{\mathrm{m}}=$ mean segment number in the system, defined in eq 6
$\mathrm{M}=$ molar mass, $\mathrm{g} / \mathrm{mol}$
$\mathrm{N}=$ total number of molecules
$\mathrm{P}=$ pressure, Pa
$\mathrm{R}=$ gas constant, $\mathrm{J} \mathrm{mol}^{-1} \mathrm{~K}^{-1}$
$\mathrm{r}=$ radial distance between two segments, Å
$\mathrm{s}_{1}=$ constant defining the pair potential, defined in eq 1 , Å
$\mathrm{T}=$ temperature, K
$\mathrm{u}(\mathrm{r})=$ pair potential function, J
$\mathrm{v}=$ molar volume, $\mathrm{m}^{3} / \mathrm{mol}$
$\mathrm{x}=$ reduced radial distance between two segments
$x_{i}=$ mole fraction of component $i$
$\mathrm{Z}=$ compressibility factor
Greek Letters
$\epsilon=$ depth of pair potential, J
$\eta=$ packing fraction, $\eta=\xi_{3}$ (see eq 9)
$\lambda=$ reduced well width of square-well potential
$\rho=$ total number density of molecules, $1 /$ Å $^{3}$
$\sigma=$ segment diameter, Å
$\xi_{\mathrm{n}}=$ abbreviation ( $\mathrm{n}=0, \ldots, 3$ ) defined by eq $9, \AA^{\mathrm{n}-3}$
Superscripts
calc $=$ calculated property
crit = critical property
disp = contribution due to dispersive attraction
$\exp =$ experimental property
hc = residual contribution of hard-chain system
hs = residual contribution of hard-sphere system
id = ideal gas contribution
sat $=$ property at saturation condition
Appendix A. Summary of Equations for Calculating Pressure, Density, Fugacity Coefficients, and Caloric Properties Using the Perturbed-Chain SAFT Model

This section provides a summary of equations for calculating thermophysical properties using the per-turbed-chain SAFT equation of state. The Helmholtz free energy $\mathrm{A}^{\text {res }}$ is the starting point in this paragraph, as all other properties can be obtained as derivatives of $\mathrm{A}^{\text {res }}$. In the following, a tilde ( $\sim$ ) will be used for reduced quantities, and caret symbols ( $\wedge$ ) will indicate molar
quantities. The reduced Helmholtz free energy, for example, is given by

$$
\begin{equation*}
\tilde{a}^{\text {res }}=\frac{A^{\text {res }}}{N k T} \tag{A.1}
\end{equation*}
$$

At the same time, one can write in terms of the molar quantity

$$
\begin{equation*}
\tilde{a}^{\text {res }}=\frac{\hat{a}^{\text {res }}}{R T} \tag{A.2}
\end{equation*}
$$

Helmholtz Free Energy. The residual Helmholtz free energy consists of the hard-chain reference contribution and the dispersion contribution.

$$
\begin{equation*}
\tilde{a}^{\text {res }}=\tilde{a}^{\text {hc }}+\tilde{a}^{\text {disp }} \tag{A.3}
\end{equation*}
$$

Hard-Chain Reference Contribution.

$$
\begin{equation*}
\tilde{a}^{h c}=\bar{m} \tilde{a n}^{h s}-\sum_{i} x_{i}\left(m_{i}-1\right) \ln g_{i i}^{h s}\left(\sigma_{i i}\right) \tag{A.4}
\end{equation*}
$$

where $\bar{m}$ is the mean segment number in the mixture.

$$
\begin{equation*}
\overline{\mathrm{m}}=\sum_{i} x_{i} \mathrm{~m}_{i} \tag{A.5}
\end{equation*}
$$

The Helmholtz free energy of the hard-sphere fluid is given on a per-segment basis

$$
\begin{align*}
\tilde{a}^{\text {hs }} & =\frac{A^{\text {hs }}}{N_{s} \mathrm{kT}} \\
& =\frac{1}{\zeta_{0}}\left[\frac{3 \zeta_{1} \zeta_{2}}{\left(1-\zeta_{3}\right)}+\frac{\zeta_{2}^{3}}{\zeta_{3}\left(1-\zeta_{3}\right)^{2}}+\left(\frac{\zeta_{2}^{3}}{\zeta_{3}^{2}}-\zeta_{0}\right) \ln \left(1-\zeta_{3}\right)\right] \tag{A.6}
\end{align*}
$$

and the radial distribution function of the hard-sphere fluid is

$$
\begin{align*}
g_{i j}^{h s}=\frac{1}{\left(1-\zeta_{3}\right)}+\left(\frac{d_{i} d_{j}}{d_{i}+d_{j}}\right) & \frac{3 \zeta_{2}}{\left(1-\zeta_{3}\right)^{2}}+ \\
& \left(\frac{d_{i} d_{j}}{d_{i}+d_{j}}\right)^{2} \frac{2 \zeta_{2}^{2}}{\left(1-\zeta_{3}\right)^{3}} \tag{A.7}
\end{align*}
$$

with $\zeta_{\mathrm{n}}$ defined as

$$
\begin{equation*}
\zeta_{n}=\frac{\pi}{6} \rho \sum x_{i} m_{i} d_{i}^{n} \quad n \in\{0,1,2,3\} \tag{A.8}
\end{equation*}
$$

The temperature-dependent segment diameter $d_{i}$ of component $i$ is given by

$$
\begin{equation*}
\mathrm{d}_{\mathrm{i}}=\sigma_{\mathrm{i}}\left[1-0.12 \exp \left(-3 \frac{\epsilon_{\mathrm{i}}}{\mathrm{kT}}\right)\right] \tag{A.9}
\end{equation*}
$$

Dispersion Contribution. The dispersion contribution to the Helmholtz free energy is given by

$$
\begin{equation*}
\tilde{\mathrm{a}}^{\mathrm{disp}}=-2 \pi \rho \mathrm{l}_{1}(\eta, \overline{\mathrm{~m}}) \mathrm{m}^{2} \epsilon \sigma^{3}-\pi \rho \overline{\mathrm{m}} \mathrm{C}_{1} \mathrm{l}_{2}(\eta, \overline{\mathrm{~m}}) \mathrm{m}^{2} \epsilon^{2} \sigma^{3} \tag{A.10}
\end{equation*}
$$

where we have introduced an abbreviation $\mathrm{C}_{1}$ for the compressibility expression, which is defined as

$$
\begin{align*}
\mathrm{C}_{1}= & \left(1+\mathrm{Z}^{\mathrm{hc}}+\rho \frac{\partial \mathrm{Z}^{\mathrm{hc}}}{\partial \rho}\right)^{-1} \\
= & \left(1+\overline{\mathrm{m}} \frac{8 \eta-2 \eta^{2}}{(1-\eta)^{4}}+\right. \\
& \left.\quad(1-\overline{\mathrm{m}}) \frac{20 \eta-27 \eta^{2}+12 \eta^{3}-2 \eta^{4}}{[(1-\eta)(2-\eta)]^{2}}\right) \tag{A.11}
\end{align*}
$$

We have also introduced the abbreviations

$$
\begin{align*}
& \overline{\mathrm{m}^{2} \epsilon \sigma^{3}}=\sum_{\mathrm{i}} \sum_{\mathrm{j}} \mathrm{x}_{\mathrm{i}} \mathrm{x}_{\mathrm{j}} \mathrm{~m}_{\mathrm{i}} \mathrm{~m}_{\mathrm{j}}\left(\frac{\epsilon_{\mathrm{ij}}}{\mathrm{kT}}\right) \sigma_{\mathrm{ij}}^{3}  \tag{A.12}\\
& \overline{\mathrm{~m}^{2} \epsilon^{2} \sigma^{3}}=\sum_{\mathrm{i}} \sum_{\mathrm{j}} \mathrm{x}_{\mathrm{i}} \mathrm{x}_{\mathrm{j}} \mathrm{~m}_{\mathrm{i}} \mathrm{~m}_{\mathrm{j}}\left(\frac{\epsilon_{\mathrm{ij}}}{\mathrm{kT}}\right)^{2} \sigma_{\mathrm{ij}}^{3} \tag{A.13}
\end{align*}
$$

Conventional combining rules are employed to determine the parameters for a pair of unlike segments.

$$
\begin{gather*}
\sigma_{\mathrm{ij}}=\frac{1}{2}\left(\sigma_{\mathrm{i}}+\sigma_{\mathrm{j}}\right)  \tag{A.14}\\
\epsilon_{\mathrm{ij}}=\sqrt{\epsilon_{\mathrm{i}} \epsilon_{\mathrm{j}}}\left(1-\mathrm{k}_{\mathrm{ij}}\right) \tag{A.15}
\end{gather*}
$$

The integrals of the perturbation theory are substituted by simple power series in density

$$
\begin{align*}
& \mathrm{I}_{1}(\eta, \overline{\mathrm{~m}})=\sum_{\mathrm{i}=0}^{6} \mathrm{a}_{\mathrm{i}}(\overline{\mathrm{~m}}) \eta^{\mathrm{i}}  \tag{A.16}\\
& \mathrm{I}_{2}(\eta, \overline{\mathrm{~m}})=\sum_{\mathrm{i}=0}^{6} \mathrm{~b}_{\mathrm{i}}(\overline{\mathrm{~m}}) \eta^{\mathrm{i}} \tag{A.17}
\end{align*}
$$

where the coefficients $a_{i}$ and $b_{i}$ depend on the chain length according to

$$
\begin{align*}
& a_{i}(\bar{m})=a_{0 i}+\frac{\bar{m}-1}{\bar{m}} a_{1 i}+\frac{\bar{m}-1}{\bar{m}} \frac{\bar{m}-2}{\bar{m}} a_{2 i}  \tag{A.18}\\
& b_{i}(\bar{m})=b_{0 i}+\frac{\bar{m}-1}{\bar{m}} b_{1 i}+\frac{\bar{m}-1}{\bar{m}} \frac{\bar{m}-2}{\bar{m}} b_{2 i} \tag{A.19}
\end{align*}
$$

The universal model constants for the $a_{0 i}, a_{1 i}, a_{2 i}, b_{0 i}$, $\mathrm{b}_{1 \mathrm{i}}$, and $\mathrm{b}_{2 \mathrm{i}}$ are given in Table 1.

Density. The density at a given system pressure $\mathrm{P}^{\text {sys }}$ must be determined iteratively by adjusting the reduced density $\eta$ until $\mathrm{P}^{\text {calc }}=\mathrm{P}^{\text {sys }}$. A suitable starting value for a liquid phase is $\eta=0.5$; for a vapor phase, $\eta=10^{-10}$. Values of $\eta>0.7405[=\pi /(3 \sqrt{ } 2)]$ are higher than the closest packing of segments and have no physical relevance. The number density of molecules $\rho$ is calculated from $\eta$ through

$$
\begin{equation*}
\rho=\frac{6}{\pi} \eta\left(\sum_{i} x_{i} m_{i} d_{i}^{3}\right)^{-1} \tag{A.20}
\end{equation*}
$$

The quantities $\zeta_{\mathrm{n}}$ given in eq A. 8 can now be cal culated. For a converged value of $\eta$, we obtain the molar density $\hat{\rho}$, in units of kmol/m³, from

$$
\begin{equation*}
\hat{\rho}=\frac{\rho}{\mathrm{N}_{\mathrm{AV}}}\left(10^{10} \frac{\AA}{\mathrm{~m}}\right)^{3}\left(10^{-3} \frac{\mathrm{kmol}}{\mathrm{~mol}}\right) \tag{A.21}
\end{equation*}
$$

where $\rho$ is, according to eq A.20, given in units of $\AA^{-3}$ and $\mathrm{N}_{\mathrm{AV}}=6.022 \times 10^{23} \mathrm{~mol}^{-1}$ denotes Avogadro's number.

Pressure. Equations for the compressibility factor will be derived using the thermodynamic relation

$$
\begin{equation*}
\mathrm{Z}=1+\eta\left(\frac{\partial \tilde{\mathrm{a}}^{\mathrm{res}}}{\partial \eta}\right)_{\mathrm{T}, \mathrm{x}_{\mathrm{i}}} \tag{A.22}
\end{equation*}
$$

The pressure can be calculated in units of $\mathrm{Pa}=\mathrm{N} / \mathrm{m}^{2}$ by applying the relation

$$
\begin{equation*}
\mathrm{P}=\mathrm{ZkT} \rho\left(10^{10} \frac{\AA \AA}{\mathrm{~m}}\right)^{3} \tag{A.23}
\end{equation*}
$$

From eqs A. 22 and A.3, it is

$$
\begin{equation*}
\mathrm{Z}=1+\mathrm{Z}^{\mathrm{hc}}+\mathrm{Z}^{\mathrm{disp}} \tag{A.24}
\end{equation*}
$$

Hard-Chain Reference Contribution. The residual hard-chain contribution to the compressibility factor is given by

$$
\begin{equation*}
z^{h c}=\bar{m} z^{h s}-\sum_{i} x_{i}\left(m_{i}-1\right)\left(g_{i i}^{h s}\right)^{-1} \rho \frac{\partial g_{i i}^{h s}}{\partial \rho} \tag{A.25}
\end{equation*}
$$

where $Z^{\text {hs }}$ is the residual contribution of the hard-sphere fluid, given by

$$
\begin{align*}
& \mathrm{Z}^{\mathrm{hs}}=\frac{\zeta_{3}}{\left(1-\xi_{3}\right)}+\frac{3 \zeta_{1} \zeta_{2}}{\zeta_{0}\left(1-\xi_{3}\right)^{2}}+\frac{3 \zeta_{2}^{3}-\xi_{3} \zeta_{2}^{3}}{\zeta_{0}\left(1-\xi_{3}\right)^{3}} \text { (A.26) }  \tag{A.26}\\
& \rho \frac{\partial \mathrm{g}_{\mathrm{ij}}^{\mathrm{hs}}}{\partial \rho}=\frac{\zeta_{3}}{\left(1-\xi_{3}\right)^{2}}+\left(\frac{\mathrm{d}_{\mathrm{i}} \mathrm{~d}_{\mathrm{j}}}{\mathrm{~d}_{\mathrm{i}}+\mathrm{d}_{\mathrm{j}}}\right)\left(\frac{3 \zeta_{2}}{\left(1-\xi_{3}\right)^{2}}+\frac{6 \xi_{2} \zeta_{3}}{\left(1-\xi_{3}\right)^{3}}\right)+ \\
&\left(\frac{\mathrm{d}_{\mathrm{i}} \mathrm{~d}_{\mathrm{j}}}{\mathrm{~d}_{\mathrm{i}}+\mathrm{d}_{\mathrm{j}}}\right)^{2}\left(\frac{4 \zeta_{2}^{2}}{\left(1-\xi_{3}\right)^{3}}+\frac{6 \zeta_{2}^{2} \xi_{3}}{\left(1-\xi_{3}\right)^{4}}\right) \text { (A.27) } \tag{A.27}
\end{align*}
$$

and $\mathrm{g}_{\mathrm{ij}}^{\text {hs }}$ was given in eq A.7.
Dispersion Contribution. The dispersion contribution to the compressibility factor can be written as

$$
\begin{align*}
& \mathrm{Z}^{\text {disp }}=-2 \pi \rho \frac{\partial\left(\eta \mathrm{I}_{1}\right)}{\partial \eta} \overline{\mathrm{m}^{2} \epsilon \sigma^{3}}- \\
& \pi \rho \overline{\mathrm{m}}\left[\mathrm{C}_{1} \frac{\partial\left(\eta \mathrm{I}_{2}\right)}{\partial \eta}+\mathrm{C}_{2} \eta \mathrm{I}_{2}\right] \overline{\mathrm{m}^{2} \epsilon^{2} \sigma^{3}} \tag{A.28}
\end{align*}
$$

where

$$
\begin{align*}
& \frac{\partial\left(\eta \mathrm{I}_{1}\right)}{\partial \eta}=\sum_{\mathrm{j}=0}^{6} \mathrm{a}_{\mathrm{j}}(\overline{\mathrm{~m}})(\mathrm{j}+1) \eta^{\mathrm{j}}  \tag{A.29}\\
& \frac{\partial\left(\eta \mathrm{I}_{2}\right)}{\partial \eta}=\sum_{\mathrm{j}=0}^{6} \mathrm{~b}_{\mathrm{j}}(\overline{\mathrm{~m}})(\mathrm{j}+1) \eta^{\mathrm{j}} \tag{A.30}
\end{align*}
$$

and where $\mathrm{C}_{2}$ is an abbreviation defined as

$$
\begin{align*}
\mathrm{C}_{2}=\frac{\partial \mathrm{C}_{1}}{\partial \eta}=- & \mathrm{C}_{1}^{2}\left(\frac{-4 \eta^{2}+20 \eta+8}{(1-\eta)^{5}}+\right. \\
& \left.(1-\overline{\mathrm{m}}) \frac{2 \eta^{3}+12 \eta^{2}-48 \eta+40}{[(1-\eta)(2-\eta)]^{3}}\right) \tag{A.31}
\end{align*}
$$

Fugacity Coefficient. The fugacity coefficient $\varphi_{\mathrm{k}}-$ ( $\mathrm{T}, \mathrm{P}$ ) is related to the residual chemical potential according to

$$
\begin{equation*}
\ln \varphi_{\mathrm{k}}=\frac{\mu_{\mathrm{k}}^{\mathrm{res}}(\mathrm{~T}, \mathrm{v})}{\mathrm{kT}}-\ln \mathrm{Z} \tag{A.32}
\end{equation*}
$$

The chemical potential can be obtained from

$$
\begin{align*}
\frac{\mu_{k}^{\text {res }}(T, v)}{k T}= & \tilde{a}^{\text {res }}+(Z-1)+ \\
& \left(\frac{\partial \tilde{a}^{\text {res }}}{\partial x_{k}}\right)_{T, v, x_{i \neq k}}-\sum_{j=1}^{N}\left[x_{j}\left(\frac{\partial \tilde{a}^{\text {res }}}{\partial x_{j}}\right)_{T, v, x_{i \neq j}}\right] \tag{A.33}
\end{align*}
$$

where derivatives with respect to mole fractions are calculated regardless of the summation relation $\sum_{\mathrm{j}} \mathrm{x}_{\mathrm{j}}=1$. For convenience, one can define abbreviations for derivatives of eq A. 8 with respect to mole fraction.

$$
\begin{equation*}
\zeta_{n, x k}=\left(\frac{\partial \zeta_{n}}{\partial x_{k}}\right)_{T, \rho, x_{j \neq k}}=\frac{\pi}{6} \rho m_{k}\left(d_{k}\right)^{n} \quad n \in\{0,1,2,3\} \tag{A.34}
\end{equation*}
$$

Hard-Chain Reference Contribution.

$$
\begin{align*}
\begin{aligned}
\left.\frac{\partial \tilde{a}^{h c}}{\partial x_{k}}\right)_{T, \rho, x_{j \neq k}}= & m_{k} \tilde{a}^{h s}+\bar{m}\left(\frac{\partial \tilde{a}^{h s}}{\partial x_{k}}\right)_{T, \rho, x_{j \neq k}}- \\
& \sum_{i} x_{i}\left(m_{i}-1\right)\left(g_{i i}^{h s}\right)^{-1}\left(\frac{\partial g_{i i}^{h s}}{\partial x_{k}}\right)_{T, \rho, x_{j \neq k}}
\end{aligned}
\end{align*}
$$

with

$$
\begin{array}{r}
\left(\frac{\partial \tilde{a}^{\text {hs }}}{\partial x_{\mathrm{k}}}\right)_{\mathrm{T}, \rho, \mathrm{x} \neq \mathrm{k}}=-\frac{\xi_{0, \mathrm{xk}}}{\zeta_{0}} \tilde{a}^{\text {hs }}+\frac{1}{\zeta_{0}}\left[\frac{3\left(\zeta_{1, \mathrm{xk}} \xi_{2}+\xi_{1} \zeta_{2, \mathrm{xk}}\right)}{\left(1-\zeta_{3}\right)}+\right. \\
\frac{3 \zeta_{1} \zeta_{2} \zeta_{3, \mathrm{xk}}}{\left(1-\xi_{3}\right)^{2}}+\frac{3 \zeta_{2}{ }^{2} \xi_{2, \mathrm{xk}}}{\xi_{3}\left(1-\xi_{3}\right)^{2}}+\frac{\zeta_{2}{ }^{3} \xi_{3, \mathrm{xk}}\left(3 \zeta_{3}-1\right)}{\zeta_{3}{ }^{2}\left(1-\xi_{3}\right)^{3}}+ \\
\left(\frac{3 \zeta_{2}{ }^{2} \xi_{2, \mathrm{xk}} \xi_{3}-2 \zeta_{2}{ }^{3} \xi_{3, \mathrm{xk}}}{\xi_{3}{ }^{3}}-\zeta_{0, \mathrm{xk}}\right) \ln \left(1-\xi_{3}\right)+ \\
\left.\left(\xi_{0}-\frac{\zeta_{2}{ }^{3}}{\zeta_{3}{ }^{2}}\right) \frac{\zeta_{3, \mathrm{xk}}}{\left(1-\xi_{3}\right)}\right] \tag{A.36}
\end{array}
$$

$$
\begin{align*}
\left(\frac{\partial g_{i j}^{h s}}{\partial x_{k}}\right)_{T, \rho, x_{j \neq k}} & =\frac{\zeta_{3, x k}}{\left(1-\zeta_{3}\right)^{2}}+ \\
& \left(\frac{d_{i} d_{j}}{d_{i}+d_{j}}\right)\left(\frac{3 \zeta_{2, x k}}{\left(1-\zeta_{3}\right)^{2}}+\frac{6 \zeta_{2} \zeta_{3, x k}}{\left(1-\zeta_{3}\right)^{3}}\right)+ \\
& \left(\frac{d_{i} d_{j}}{d_{i}+d_{j}}\right)^{2}\left(\frac{4 \zeta_{2} \zeta_{2, x k}}{\left(1-\zeta_{3}\right)^{3}}+\frac{6 \zeta_{2}^{2} \zeta_{3, x k}}{\left(1-\zeta_{3}\right)^{4}}\right) \tag{A.37}
\end{align*}
$$

Dispersion Contribution.

$$
\begin{array}{r}
\left(\frac{\partial \tilde{a}^{h \mathrm{~s}}}{\partial \mathrm{x}_{\mathrm{k}}}\right)_{\mathrm{T}, \rho, \mathrm{x}_{\mathrm{j} \neq \mathrm{k}}}=-2 \pi \rho\left[\mathrm{I}_{1, \mathrm{xk}} \overline{\mathrm{~m}^{2} \epsilon \sigma^{3}}+\mathrm{I}_{1}\left(\overline{\mathrm{~m}^{2} \epsilon \sigma^{3}}\right)_{\mathrm{xk}}\right]- \\
\pi \rho\left\{\left[\mathrm{m}_{\mathrm{k}} \mathrm{C}_{1} \mathrm{I}_{2}+\overline{\mathrm{m} \mathrm{C}_{1, \mathrm{xk}} \mathrm{I}_{2}}+\overline{\left.\mathrm{m} \mathrm{C}_{1} \mathrm{I}_{2, \mathrm{xk}}\right] \overline{\mathrm{m}^{2} \epsilon^{2} \sigma^{3}}+}\right.\right. \\
\overline{\left.\mathrm{m} \mathrm{C}_{1} \mathrm{I}_{2}\left(\mathrm{~m}^{2} \epsilon^{2} \sigma^{3}\right)_{\mathrm{xk}}\right\}}(\mathrm{A} \tag{A.38}
\end{array}
$$

with

$$
\begin{array}{r}
\left(\overline{\mathrm{m}^{2} \epsilon \sigma^{3}}\right)_{\mathrm{xk}}=2 \mathrm{~m}_{\mathrm{k}} \sum_{\mathrm{j}} \mathrm{x}_{\mathrm{j}} \mathrm{~m}_{\mathrm{j}}\left(\frac{\epsilon_{\mathrm{kj}}}{\mathrm{kT}}\right) \sigma_{\mathrm{kj}}^{3} \\
\left(\overline{\mathrm{~m}^{2} \epsilon^{2} \sigma^{3}}\right)_{\mathrm{xk}}=2 \mathrm{~m}_{\mathrm{k}} \sum_{\mathrm{j}} \mathrm{x}_{\mathrm{j}} \mathrm{~m}_{\mathrm{j}}\left(\frac{\epsilon_{\mathrm{kj}}}{\mathrm{kT}}\right)^{2} \sigma_{\mathrm{kj}}^{3} \\
\mathrm{C}_{1, \mathrm{xk}}=\mathrm{C}_{2} \zeta_{3, \mathrm{xk}}- \\
\mathrm{C}_{1}^{2}\left\{\mathrm{~m}_{\mathrm{k}} \frac{8 \eta-2 \eta^{2}}{(1-\eta)^{4}}-\mathrm{m}_{\mathrm{k}} \frac{20 \eta-27 \eta^{2}+12 \eta^{3}-2 \eta^{4}}{[(1-\eta)(2-\eta)]^{2}}\right\} \tag{A.41}
\end{array}
$$

$$
\begin{equation*}
\mathrm{I}_{1, \mathrm{xk}}=\sum_{\mathrm{i}=0}^{6}\left[\mathrm{a}_{\mathrm{i}}(\overline{\mathrm{~m}}) \mathrm{i} \xi_{3, \mathrm{xk}} \eta^{\mathrm{i}-1}+\mathrm{a}_{\mathrm{i}, \mathrm{xk}} \eta^{\mathrm{i}}\right] \tag{A.42}
\end{equation*}
$$

$$
\begin{equation*}
\mathrm{I}_{2, \mathrm{xk}}=\sum_{\mathrm{i}=0}^{6}\left[\mathrm{~b}_{\mathrm{i}}(\overline{\mathrm{~m}}) \mathrm{i} \zeta_{3, \mathrm{xk}} \eta^{\mathrm{i}-1}+\mathrm{b}_{\mathrm{i}, \mathrm{xk}} \eta^{\mathrm{i}}\right] \tag{A.43}
\end{equation*}
$$

$$
\begin{equation*}
a_{i, x k}=\frac{m_{k}}{\bar{m}^{2}} a_{1 i}+\frac{m_{k}}{\bar{m}^{2}}\left(3-\frac{4}{\bar{m}}\right) a_{2 i} \tag{A.44}
\end{equation*}
$$

$$
\begin{equation*}
\mathrm{b}_{\mathrm{i}, \mathrm{xk}}=\frac{\mathrm{m}_{\mathrm{k}}}{\overline{\mathrm{~m}}^{2}} \mathrm{~b}_{1 \mathrm{i}}+\frac{\mathrm{m}_{\mathrm{k}}}{\overline{\mathrm{~m}}^{2}}\left(3-\frac{4}{\overline{\mathrm{~m}}}\right) \mathrm{b}_{2 \mathrm{i}} \tag{A.45}
\end{equation*}
$$

Enthalpy and Entropy. The molar enthalpy $\hat{h}^{\text {res }}$ is obtained from a derivative of the Helmholtz free energy with respect to temperature, according to

$$
\begin{equation*}
\frac{\hat{h}^{\mathrm{res}}}{\mathrm{RT}}=-\mathrm{T}\left(\frac{\partial \tilde{a}^{\mathrm{res}}}{\partial \mathrm{~T}}\right)_{\rho, \mathrm{x}_{\mathrm{i}}}+(\mathrm{Z}-1) \tag{A.46}
\end{equation*}
$$

Unlike the enthalpy of an ideal gas, which is a function of temperature only, the entropy of an ideal gas is a function of both temperature and pressure (or density). Hence, the residual entropy in the variables P and T is different from the residual entropy for the specified conditions v and T . It is

$$
\begin{equation*}
\frac{\hat{s}^{\text {res }}(P, T)}{R}=\frac{\hat{s}^{\text {res }}(v, T)}{R}+\ln (Z) \tag{A.47}
\end{equation*}
$$

All of the equations for anes are given in the variables v and T , so that the residual entropy can be written as

$$
\begin{equation*}
\frac{\hat{\mathrm{s}}^{\text {res }}(\mathrm{P}, \mathrm{~T})}{\mathrm{R}}=-\mathrm{T}\left[\left(\frac{\partial \tilde{\mathrm{ar}}^{\text {res }}}{\partial \mathrm{T}}\right)_{\rho, \mathrm{x}_{\mathrm{i}}}+\frac{\tilde{\mathrm{a}}^{\text {res }}}{\mathrm{T}}\right]+\ln (\mathrm{Z}) \tag{A.48}
\end{equation*}
$$

The residual molar Gibbs free energy $\hat{g}^{\text {res }}(\mathrm{P}, \mathrm{T})$ is defined as

$$
\begin{equation*}
\frac{\hat{g}^{\text {res }}}{R T}=\frac{\hat{h}^{\text {res }}}{R T}-\frac{\hat{s}^{\text {res }}(P, T)}{R} \tag{A.49}
\end{equation*}
$$

or simply as

$$
\begin{equation*}
\frac{\hat{g}^{\text {res }}}{\mathrm{RT}}=\tilde{a}^{\text {res }}+(\mathrm{Z}-1)-\ln (\mathrm{Z}) \tag{A.50}
\end{equation*}
$$

The temperature derivative of ãres in eqs A. 46 and A. 48 is again the sum of two contributions.

$$
\begin{equation*}
\left(\frac{\partial \tilde{a}^{\text {res }}}{\partial \mathrm{T}}\right)_{\rho, \mathrm{x}_{\mathrm{i}}}=\left(\frac{\partial \tilde{a}^{\text {hc }}}{\partial \mathrm{T}}\right)_{\rho, \mathrm{x}_{\mathrm{i}}}+\left(\frac{\partial \tilde{a}^{\text {disp }}}{\partial \mathrm{T}}\right)_{\rho, \mathrm{x}_{\mathrm{i}}} \tag{A.51}
\end{equation*}
$$

With abbreviations for two temperature derivatives

$$
\begin{array}{r}
\mathrm{d}_{\mathrm{i}, \mathrm{~T}}=\frac{\partial \mathrm{d}_{\mathrm{i}}}{\partial \mathrm{~T}}=\sigma_{\mathrm{i}}\left(3 \frac{\epsilon_{\mathrm{i}}}{\mathrm{kT}^{2}}\right)\left[-0.12 \exp \left(-3 \frac{\epsilon_{\mathrm{i}}}{\mathrm{kT}}\right)\right] \text { (A.52 } \\
\zeta_{\mathrm{n}, \mathrm{~T}}=\frac{\partial \zeta_{\mathrm{n}}}{\partial \mathrm{~T}}=\frac{\pi}{6} \rho \sum_{\mathrm{i}} \mathrm{x}_{\mathrm{i}} \mathrm{~m}_{\mathrm{i}} \mathrm{nd}_{\mathrm{i}, \mathrm{~T}}\left(\mathrm{~d}_{\mathrm{i}}\right)^{\mathrm{n}-1} \quad \mathrm{n} \in\{1,2,3\} \tag{A.53}
\end{array}
$$

the hard-chain contribution and the dispersion contribution can conveniently be written.

Hard-Chain Reference Contribution.

$$
\begin{align*}
& \left.\left(\frac{\partial \tilde{a}^{\text {nc }}}{\partial \mathrm{T}}\right)_{\rho, \mathrm{x}_{\mathrm{i}}}=\overline{\mathrm{m}}\left(\frac{\partial \tilde{a}^{\text {ns }}}{\partial \mathrm{T}}\right)_{\rho, \mathrm{x}_{\mathrm{i}}}-\sum_{\mathrm{i}} \mathrm{x}_{\mathrm{i}}\left(\mathrm{~m}_{\mathrm{i}}-1\right)\left(\mathrm{g}_{\mathrm{ii}}^{\text {hs }}\right)^{-1}\left(\frac{\partial \mathrm{~g}_{\mathrm{ii}}^{\text {ns }}}{\partial \mathrm{T}}\right)_{\left(\mathrm{A} . \mathrm{x}_{\mathrm{i}}\right.} 4\right)  \tag{A.54}\\
& \left(\frac{\partial \tilde{a}^{\text {hs }}}{\partial \mathrm{T}}\right)_{\rho, \mathrm{x}_{\mathrm{i}}}=\frac{1}{\zeta_{0}}\left[\frac{3\left(\zeta_{1, \mathrm{~T}} \zeta_{2}+\zeta_{1} \zeta_{2, \mathrm{~T}}\right)}{\left(1-\xi_{3}\right)}+\frac{3 \zeta_{1} \zeta_{2} \zeta_{3, \mathrm{~T}}}{\left(1-\xi_{3}\right)^{2}}+\right. \\
& \frac{3 \zeta_{2}{ }^{2} \zeta_{2, \mathrm{~T}}}{\zeta_{3}\left(1-\zeta_{3}\right)^{2}}+\frac{\zeta_{2}{ }^{3} \zeta_{3, \mathrm{~T}}\left(3 \zeta_{3}-1\right)}{\zeta_{3}{ }^{2}\left(1-\zeta_{3}\right)^{3}}+ \\
& \left.\left(\frac{3 \zeta_{2}{ }^{2} \zeta_{2, \mathrm{~T}} \zeta_{3}-2 \zeta_{2}{ }^{3} \zeta_{3, \mathrm{~T}}}{\zeta_{3}{ }^{3}}\right) \ln \left(1-\zeta_{3}\right)+\left(\zeta_{0}-\frac{\zeta_{2}{ }^{3}}{\zeta_{3}{ }^{2}}\right) \frac{\zeta_{3, \mathrm{~T}}}{\left(1-\zeta_{3}\right)}\right] \tag{A.55}
\end{align*}
$$

Equation A. 54 requires only the i-i pairs in the temperature derivative of the radial pair distribution function $\mathrm{g}_{\mathrm{ij}}^{\mathrm{hs}}$. For simplicity, one can restrict oneself to the i-i pairs in eq A. 7 by equating

$$
\begin{equation*}
\frac{1}{2} d_{i}=\left(\frac{d_{i} d_{i}}{d_{i}+d_{i}}\right) \tag{A.56}
\end{equation*}
$$

The temperature derivative of the radial pair distribution function $\mathrm{g}_{\mathrm{ii}}^{\mathrm{hs}}$ is then

$$
\begin{align*}
& \frac{\partial g_{i i}^{\text {hs }}}{\partial T}=\frac{\zeta_{3, T}}{\left(1-\xi_{3}\right)^{2}}+\left(\frac{1}{2} d_{i, T}\right) \frac{3 \zeta_{2}}{\left(1-\xi_{3}\right)^{2}}+ \\
& \left(\frac{1}{2} d_{i}\right)\left(\frac{3 \zeta_{2, T}}{\left(1-\xi_{3}\right)^{2}}+\frac{6 \zeta_{2} \xi_{3, T}}{\left(1-\xi_{3}\right)^{3}}\right)+\left(\frac{1}{2} d_{i} d_{i, T}\right) \frac{2 \zeta_{2}^{2}}{\left(1-\xi_{3}\right)^{3}}+ \\
& \left(\frac{1}{2} d_{i}\right)^{2}\left(\frac{4 \zeta_{2} \xi_{2, T}}{\left(1-\xi_{3}\right)^{3}}+\frac{6 \xi_{2}^{2} \xi_{3, T}}{\left(1-\xi_{3}\right)^{4}}\right)^{2} \text { (A.5 } \tag{A.57}
\end{align*}
$$

Dispersion Contribution.

$$
\begin{align*}
\left(\frac{\partial \tilde{a}^{\text {disp }}}{\partial \mathrm{T}}\right)_{\rho, \mathrm{x}_{\mathrm{i}}} & =-2 \pi \rho\left(\frac{\partial \mathrm{I}_{1}}{\partial \mathrm{~T}}-\frac{\mathrm{I}_{1}}{\mathrm{~T}}\right) \overline{\mathrm{m}^{2} \epsilon \sigma^{3}}- \\
& \pi \rho \overline{\mathrm{m}}\left[\frac{\partial \mathrm{C}_{1}}{\partial \mathrm{~T}} \mathrm{I}_{2}+\mathrm{C}_{1} \frac{\partial \mathrm{I}_{2}}{\partial \mathrm{~T}}-2 \mathrm{C}_{1} \frac{\mathrm{I}_{2}}{\mathrm{~T}}\right] \overline{\mathrm{m}^{2} \epsilon^{2} \sigma^{3}} \tag{A.58}
\end{align*}
$$

with

$$
\begin{gather*}
\frac{\partial \mathrm{I}_{1}}{\partial \mathrm{~T}}=\sum_{\mathrm{I}=0}^{6} \mathrm{a}_{\mathrm{i}}(\overline{\mathrm{~m}}) \mathrm{i} \zeta_{3, \mathrm{~T}} \eta^{\mathrm{i}-1}  \tag{A.59}\\
\frac{\partial \mathrm{I}_{2}}{\partial \mathrm{~T}}=\sum_{\mathrm{I}=0}^{6} \mathrm{~b}_{\mathrm{i}}(\overline{\mathrm{~m}}) \mathrm{i} \zeta_{3, \mathrm{~T}} \eta^{\mathrm{i}-1}  \tag{A.60}\\
\frac{\partial \mathrm{C}_{1}}{\partial \mathrm{~T}}=\zeta_{3, \mathrm{~T}} \mathrm{C}_{2} \tag{A.61}
\end{gather*}
$$

## Appendix B. Correlations for Pure-Component Parameters of the Perturbed-Chain SAFT Equation of State.

Equation 18, which was given the purpose of model development, is a well suitable function for correlating pure-component parameters with varying segment number. It is suitable, because it allows for varying parameters of short chains but converges to constant values as segment number increases. Generally, this equation captures the effect of chain length on physical properties. It is convenient, though, to modify this equation to obtain a correlation of pure-component parameters with molar mass, rather than with varying segment number. For convenience, we choose the molar mass of one hydrocarbon unit to be equal to the molar mass of methane ( $\mathrm{M}_{\mathrm{CH}_{4}}=16.043 \mathrm{~g} / \mathrm{mol}$ ) and obtain the following relation for the segment diameter of the n -alkane series:

$$
\begin{equation*}
\sigma_{\mathrm{i}}=\mathrm{q}_{01}+\frac{\mathrm{M}_{\mathrm{i}}-\mathrm{M}_{\mathrm{CH}_{4}}}{\mathrm{M}_{\mathrm{i}}} \mathrm{q}_{11}+\frac{\mathrm{M}_{\mathrm{i}}-\mathrm{M}_{\mathrm{CH}_{4}}}{\mathrm{M}_{\mathrm{i}}} \frac{\mathrm{M}_{\mathrm{i}}-2 \mathrm{M}_{\mathrm{CH}_{4}}}{\mathrm{M}_{\mathrm{i}}} \mathrm{q}_{21} \tag{B.1}
\end{equation*}
$$

For the ratio of the segment number to the molar mass ( $\mathrm{m}_{\mathrm{i}} / \mathrm{M}_{\mathrm{i}}$ ) and the energy parameter $\epsilon_{\mathrm{i}} / \mathrm{k}$ of the n -alkanes, we obtain

$$
\begin{align*}
&\left(m_{i} / M_{i}\right)=q_{02}+\frac{M_{i}-M_{C H_{4}}}{M_{i}} q_{12}+ \\
& \frac{M_{i}-M_{C H_{4}}}{M_{i}} \frac{M_{i}-2 M_{C H_{4}}}{M_{i}} q_{22} \tag{B.2}
\end{align*}
$$

$$
\begin{align*}
&\left(\epsilon_{i} / k\right)=q_{03}+\frac{M_{i}-M_{C H_{4}}}{M_{i} q_{13}}+ \\
& \frac{M_{i}-M_{C H_{4}}}{M_{i}} \frac{M_{i}-2 M_{C H_{4}}}{M_{i}} q_{23} \tag{B.3}
\end{align*}
$$

where $\mathrm{q}_{\mathrm{jk}}$ are constants that can be fitted to the purecomponent parameter given in Table 2. For the n-alkane series, these constants are

| $j$ | units | 0 | 1 | 2 |
| :--- | :--- | :---: | :---: | :---: |
| $\mathrm{q}_{\mathrm{j} 1}$ | $\AA$ | 3.7039 | -0.3226 | 0.6907 |
| $\mathrm{q}_{\mathrm{j} 2}$ | $\mathrm{~mol} / \mathrm{g}$ | 0.06233 | -0.02236 | -0.01563 |
| $\mathrm{q}_{\mathrm{j} 3}$ | K | 150.03 | 80.68 | 38.96 |

## Literature Cited

[^2](4) Donohue, M. D.; Prausnitz, J. M. Perturbed Hard Chain Theory for Fluid Mixtures: Thermodynamic Properties for Mixtures in Natural Gas and Petroleum Technology. AIChE J . 1978, 24, 849.
(5) Wertheim, M. S. Fluids with highly directional attractive forces. I. Statistical thermodynamics. J. Stat. Phys. 1984, 35, 19.
(6) Wertheim, M. S. Fluids with highly directional attractive forces. II. Thermodynamic perturbation theory and integral equations. J . Stat. Phys. 1984, 35, 35.
(7) Wertheim, M. S. Fluids with highly directional attractive forces. III. Multiple attraction sites. J . Stat. Phys. 1986, 42, 459.
(8) Wertheim, M. S. Fluids with highly directional attractive forces. IV. Equilibrium polymerization. J . Stat. Phys. 1986, 42, 477.
(9) Chapman, W. G.; J ackson, G.; Gubbins, K. E. Phase equilibria of associating fluids. Chain molecules with multiple bonding sites. Mol. Phys. 1988, 65, 1057.
(10) Chapman, W. G.; Gubbins, K. E.; J ackson, G.; Radosz, M. New Reference Equation of State for Associating Liquids. Ind. Eng. Chem. Res. 1990, 29, 1709.
(11) Ghonasgi, D.; Chapman, W. G. Prediction of the Properties of Model Polymer Solutions and Blends. AIChE J . 1994, 40, 878.
(12) Banaszak, M.; Chiew, Y. C.; O'Lenick, R.; Radosz, M. Thermodynamic Perturbation Theory: Lennard-J ones Chains. J . Chem. Phys. 1994, 100, 3803.
(13) J ohnson, J . K.; Müller, E. A.; Gubbins, K. E. Equation of State for Lennard-J ones Chains. J . Phys. Chem. 1994, 98, 6413.
(14) Müller, E. A.; Vega, L. F.; Gubbins, K. E. Theory and simulation of associating fluids: Lennard-J ones chains with association sites. Mol. Phys. 1994, 83, 1209.
(15) Kraska, T.; Gubbins, K. E. Phase Equilibria Calculations with a Modified SAFT Equation of State. 1. Pure Alkanes, Alkanols, and Water. Ind. Eng. Chem. Res. 1996, 35, 4727.
(16) Kraska, T.; Gubbins, K. E. Phase Equilibria Calculations with a Modified SAFT Equation of State. 2. Binary Mixtures of n-Alkanes, 1-Alkanols, and Water. I nd. Eng. Chem. Res. 1996, 35, 4738.
(17) Blas, F. J.; Vega, L. F. Thermodynamic behaviour of homonuclear and heteronuclear Lennard-J ones chains with association sites from simulation and theory. Mol. Phys. 1997, 92, 135.
(18) Gil-Villegas, A.; Galindo, A.; Whitehead, P. J .; Mills, S. J .; J ackson, G.; Burgess, A. N. Statistical associating fluid theory for chain molecules with attractive potentials of variable range. J . Chem. Phys. 1997, 106, 4168.
(19) Huang, S. H.; Radosz, M. Equation of State for Small, Large, Polydisperse, and Associating Molecules. Ind. Eng. Chem. Res. 1990, 29, 2284.
(20) Huang, S. H.; Radosz, M. Equation of State for Small, Large, Polydisperse and Associating Molecules: Extensions to Fluid Mixtures. Ind. Eng. Chem. Res. 1991, 30, 1994.
(21) Chen, S. S.; Kreglewski, A. Applications of the Augmented van der Waals Theory of Fluids. I. Pure Fluids. Ber. Bunsen-Ges. 1977, 81 (10), 1048.
(22) Gross, J .; Sadowski, G. Application of perturbation theory to a hard-chain reference fluid: An equation of state for squarewell chains. Fluid Phase Equilib. 2000, 168, 183.
(23) Barker, J . A.; Henderson, D. Perturbation Theory and Equation of State for Fluids: The Square-Well Potential. J . Chem. Phys. 1967, 47, 2856.
(24) Barker, J . A.; Henderson, D. Perturbation Theory and Equation of State for Fluids. II. A Successful Theory of Liquids. J . Chem. Phys. 1967, 47, 4714.
(25) Alder, B. J .; Young, D. A.; Mark, M. A. Studies in Molecular Dynamics. X. Corrections to the Augmented van der Waals Theory for the Square-Well Fluid. J . Chem. Phys. 1972, 56, 3013.
(26) Boublik, T. Hard-Sphere Equation of State. J . Chem. Phys. 1970, 53, 471.
(27) Mansoori, G. A.; Carnahan, N. F.; Starling, K. E.; Leland, T. W. Equilibrium Thermodynamic Properties of the Mixture of Hard Spheres. J . Chem. Phys. 1971, 54, 1523.
(28) Chiew, Y. C. Percus-Yevick integral equation theory for athermal hard-sphere chains. II. Average intermolecular correlation functions. Mol. Phys. 1991, 73, 359.
(29) Tang, Y.; Lu, B. C.-Y. Direct calculations of radial distribution function for hard-sphere chains. J. Chem. Phys. 1996, 105, 8262.
(30) Chang, J .; Kim, H. Analytical expression for the correlation function of a hard sphere chain fluid. Mol. Phys. 1999, 96, 1789.
(31) Liu, H.; Hu, Y. Molecuar thermodynamic theory for polymer systems. II. Equation of state for chain fluids. Fluid Phase Equilib. 1996, 122, 75.
(32) Cummings, P. T.; Stell, G. Statistical mechanical models of chemical reactions, analytic solution of models in the PercusYevick approximation. Mol. Phys. 1984, 51, 253.
(33) Cummings, P. T.; Stell G. Statistical mechanical models of chemical reactions. II. Analytic solution of the Percus-Yevick approximation for a model of homogeneous association. Mol. Phys. 1985, 55, 33.
(34) O'Lenick, R.; Li, X. J .; Chiew, Y. C. Correlation functions of hard-sphere chain mixtures: Integral equation theory and simulation results. Mol. Phys. 1995, 86, 1123.
(35) Han, S. J.; Lin, H. M.; Chao, K. C. Vapor-Liquid Equilibrium of Molecular Fluid Mixtures by Equation of State. Chem. Eng. Sci. 1988, 43, 2327.
(36) Kiran, E.; Xiong, Y.; Zhuang, W. Effect of Polydispersity on the Demixing Pressure of Polyethylene in Near- or Supercritical Alkanes. J . Supercrit. Fluids 1994, 7, 283.
(37) Koak, N.; Visser, R. M.; de Loos, Th. W. High-pressure phase behavior of the systems polyethylene + ethylene and polybutene + 1-butene. Fluid Phase Equilib. 1999, 158-160, 835.
(38) Lee, S.-H.; Hasch, B. M.; McHugh, M. A. Calculating copolymer solution behavior with Statistical Associating Fluid Theory. Fluid Phase Equilib. 1996, 117, 61.
(39) Gross, J. Entwicklung einer Zustandsgleichung für einfache, assoziierende und makromolekulare Stoffe. Dissertation, Technische Universität Berlin, Germany, 2000.

Received for review April 24, 2000
Revised manuscript received October 24, 2000
Accepted November 1, 2000
IE0003887


[^0]:    * Author to whom correspondence should be addressed. E-mail: G.Sadowski@ct.uni-dortmund.de. Fax: ++49-2317552572.

[^1]:    ${ }^{\mathrm{a}}$ References: (1) Kalra et al. J . Chem. Eng. Data 1976, 21, 222-225. (2) Leu and Robinson. J . Chem. Eng. Data 1987, 32, 444. (3) Robinson and Kalra. Proceedings of the 53rd Annual Convention of the Gas Processors Association, 1977. (4) Sage et al. Ind. Eng. Chem. 1940, 32, 1085. (5) Sage et al. Ind. Eng. Chem. 1942, 34, 1108. (6) Gunn et al. AICHE J . 1974, 20, 347. (7) Lin et al. J. Chem. Eng. Data 1977, 22, 402. (8) Reamer et al. Ind. Eng. Chem. Data Ser. 1 1956, 29. (9) Lin et al. J. Chem. Eng. Data 1979, 24, 146. (10) Olds et al. Ind. Eng. Chem. 1942, 34, 1008. (11) Reamer et al. Ind. Eng. Chem. Data Ser. 3 1958, 240. (12) Sebastian et al. Fluid Phase Equilib. 1979, 3, 145. (13) Sebastian et al. J. Chem. Eng. Data 1979, 24, 149. (14) Grauso et al. Fluid Phase Equilib. 1977, 1, 13. (15) Poston et al. J. Chem. Eng. Data 1966, 11, 364. (16) Kay. Ind. Eng. Chem. 1948, 40, 1459. (17) Ohgaki et al. Fluid Phase Equilib. 1982, 8, 113. (18) Kay. J. Chem. Eng. Data 1970, 15, 46. (19) Gilliland and Scheeline. Ind. Eng. Chem. 1940, 32, 48. (20) Cummings et al. Ind. Eng. Chem. 1933, 25, 728. (21) Wei et al. J . Chem. Eng. Data 1995, 40, 726. (22) Nagahama et al. J . Chem. Eng. J pn. 1974, 7, 323. (23) Yucelen and Kidnay. J. Chem. Eng. Data 1999, 44, 926. (24) Besserer and Robinson. J. Chem. Eng. Data 1973, 18. (25) Kalra et al. J. Chem. Eng. Data 1978, 23, 317. (26) Ng and Robinson. Fluid Phase Equilib. 1979, 2, 283. (26) Rodrigues et al. J. Chem. Eng. Data 1968, 13, 165.

[^2]:    (1) Sanchez, I. C.; Lacombe, R. H. Theory of liquid-liquid and liquid-vapor equilibria. Nature (London) 1974, 252, 381.
    (2) Wei, Y. S.; Sadus, R. J. Equations of State for the Calculation of Fluid-Phase Equilibria. AIChE J . 2000, 46, 169.
    (3) Beret, S.; Prausnitz, J . M. Perturbed Hard-Chain Theory: An Equation of State for Fluids Containing Small or Large Molecules. AIChE J . 1975, 21, 1123.

