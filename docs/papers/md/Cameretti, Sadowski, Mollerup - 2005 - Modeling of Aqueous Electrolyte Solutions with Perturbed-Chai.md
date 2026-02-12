# Modeling of Aqueous Electrolyte Solutions with Perturbed-Chain Statistical Associated Fluid Theory 

Luca F. Cameretti and Gabriele Sadowski*<br>Department of Biochemical and Chemical Engineering, University of Dortmund, 44227 Dortmund, Germany<br>Jørgen M. Mollerup<br>IVC-SEP, Department of Chemical Engineering, Technical University of Denmark, 2800 Lyngby, Denmark

Downloaded via BRIGHAM YOUNG UNIV on May 18, 2024 at 02:56:36 (UTC). See https://pubs.acs.org/sharingguidelines for options on how to legitimately share published articles.


#### Abstract

The vapor pressures and liquid densities of single-salt electrolyte solutions containing NaCl, $\mathrm{LiCl}, \mathrm{KCl}, \mathrm{NaBr}, \mathrm{LiBr}, \mathrm{KBr}, \mathrm{NaI}, \mathrm{LiI}, \mathrm{KI}, \mathrm{Li}_{2} \mathrm{SO}_{4}, \mathrm{Na}_{2} \mathrm{SO}_{4}$, and $\mathrm{K}_{2} \mathrm{SO}_{4}$ were modeled with an equation of state based on perturbed-chain statistical associated fluid theory (PC-SAFT). The PC-SAFT model was extended to charged compounds using a Debye-Hückel term for the electrostatic interactions. Two model parameters for each ion were fitted to experimental $p V T$ and vapor-pressure data. The model is able to excellently reproduce the experimental data up to high salt molalities and even to predict vapor pressures in mixed-salt solutions.


## 1. Introduction

Electrolyte solutions are commonly encountered in many chemical processes including desalination, wastewater treatment, and extractive distillation. ${ }^{1}$ They also play an important role in biotechnology. For example, salts are used to precipitate proteins out of the filtered fermentation broth and enhance their crystallization. ${ }^{2}$ Conditions for this purification step, i.e., temperature, salt concentration, etc., are still found by trail-and-error methods. ${ }^{3}$ The knowledge of phase equilibria of protein solutions containing electrolytes is not only important for production processes. Precipitation of protein aggregates is suspected of being responsible for a couple of diseases such as Alzheimer's, Creuzfeldt-Jacob disease, and eye cataracts, just to name a few. ${ }^{4}$ Therefore, many research groups are interested in phase equilibria of protein solutions containing salts, on the one hand, to optimize process variables and, on the other hand, to find appropriate remedies against diseases. ${ }^{5,6}$

The first step toward modeling complex protein solutions is to be able to accurately model electrolyte systems. Since the beginning of the last century, many researchers have concentrated their efforts in correlating and predicting phase equilibria in electrolyte solutions. In 1923, Debye and Hückel ${ }^{7}$ published one of the pioneering papers dealing with aqueous electrolyte solutions. They considered a system of hard spheres in a dielectric continuum and calculated the contribution to the system energy for charging up the spheres. Their model is often treated as an excess Gibbs energy model; however, what they developed was an electrostatic contribution to the Helmholtz energy when charging hard spheres at constant temperature, volume, and chemical potential of the solvent. Their development includes the background potential as well as the selfpotential. The latter is also known as the Born potential, which is independent of the configuration as long as we are dealing with pure solvents.
$g^{\mathrm{E}}$ models such as the local composition nonrandom two-liquid activity coefficient model by Chen et al. ${ }^{8}$ and

[^0]its extensions ${ }^{9-11}$ or other group contribution models ${ }^{12,13}$ have been applied to correlate activity coefficients of single-salt and mixed-salt electrolyte solutions. However, the great disadvantage of excess Gibbs enthalpy models is that they are not able to predict densities. Equation of state (EOS) models circumvent this crucial disadvantage. Several models based on EOSs were extended to electrolyte systems. Fürst and Renon ${ }^{14}$ combined the nonelectrolyte EOS of Schwartzentruber et al. with a mean spherical approximation (MSA) longrange term to account for the electrostatic interactions. Wu and Prausnitz ${ }^{15}$ presented an electrolyte EOS based on the Peng-Robinson EOS (PREOS), which itself already accounts for hard-core and dispersion interactions. They added a Born energy term for charging up the uncharged reference system in a continuous medium of given permittivity, a Coulombic term to account for electrostatic interactions between the ions, and further an association term to consider hydrogen bonds between the water molecules. Myers et al. ${ }^{16}$ followed a similar approach using the ideal gas mixture as a reference system, considering also the Born energy to discharge the ions in a vacuum, and then recharging them in the dielectric medium. Both approaches yield good agreement with experimental vapor-pressure data. However, they are based on the semiempirical PREOS, which is known to fail when predicting liquid densities.

A series of electrolyte EOSs are based on perturbation theories. Significant progress in developing powerful EOSs based on perturbation theories was achieved especially over the last 2 decades. Starting from an appropriate reference system, these theories have been proven to successfully account for the nonspherical shape of the molecules, for nonspecific attractive interactions, and for specific interactions such as association as well as for polar or quadrupolar interactions.

A whole series of models of this kind is based on the statistical associated fluid theory (SAFT), ${ }^{17-20}$ which considers a molecule as a chain of tangent spherical segments. Starting from the Helmholtz energy of a hardsphere reference system $a^{\mathrm{hs}}$, different perturbation contributions are considered that are assumed to effect independently, namely, attractive interactions of the $m$
(nonbonded) segments ( $a^{\text {disp }}$ ), hard-sphere chain formation ( $a^{\text {chain }}$ ), and association ( $a^{\text {assoc }}$ ):

$$
\begin{equation*}
a^{\mathrm{res}}=A^{\mathrm{res}} / N=m a^{\mathrm{hs}}+m a^{\mathrm{disp}}+a^{\mathrm{chain}}+a^{\mathrm{assoc}} \tag{1}
\end{equation*}
$$

where $N$ is the total number of molecules. In SAFT, the Carnahan-Starling expression is used for the segment hard-sphere contribution $a^{\mathrm{hs}, 21}$ the segment-segment dispersion $a^{\text {disp }}$ is described by a fourth-order perturbation term; ${ }^{22,23}$ the contribution of chain formation as well as the association term is accounted for based on the work of Wertheim. ${ }^{24-27}$ Subsequently, various perturbation theories were suggested that differ in the use of specific perturbation expressions. Examples are the perturbed hard-sphere-chain theory ${ }^{28}$ as well as the models proposed by Chang and Sandler, ${ }^{29}$ Gil-Villegas et al., ${ }^{30}$ and Hino and Prausnitz. ${ }^{31}$ Each of these models considers the nonspherical shape of a molecule, on the one hand, and the attractive interaction, on the other hand, as independent perturbations of the reference system. Several attempts have been made to overcome this deficiency. Various models were suggested that use the square-well sphere (e.g., refs 30,32 , and 33 ) or the Lennard-Jones sphere (e.g., refs 34-36) rather than the hard sphere as the reference to modify the chain contribution $a^{\text {chain }}$. The recently proposed perturbedchain SAFT (PC-SAFT) model ${ }^{37,38}$ adopts the opposite idea: here, a second-order perturbation theory is applied to the reference system of hard chains instead of hard spheres to develop a dispersion term $a^{\text {disp }}$ for nonspherical molecules.

$$
\begin{equation*}
a^{\mathrm{res}}=a^{\mathrm{hc}}+a^{\mathrm{disp}}(m)+a^{\mathrm{assoc}} \tag{2}
\end{equation*}
$$

Whereas the contributions to describe the hard-chain system as well as association are identical with those of the original SAFT model, the dispersion term was modified to account for the influence of the nonspherical shape of the molecule on the number of intermolecular interactions and is therefore a function of segment number $m$. The PC-SAFT model was successfully applied to a wide variety of uncharged systems, demonstrating that the modeling results could not only be improved for chain molecules such as polymers. Even for smaller nonspherical substances, the modeling results could be improved considerably compared to the original SAFT model. ${ }^{36-40}$

To describe electrolyte solutions, Galindo et al. ${ }^{41,42}$ have successfully extended the SAFT-VR EOS ${ }^{30}$ to electrolyte solutions by using an additive electrostatic term for the Helmholtz free energy obtained from the solution of the Ornstein-Zernicke equation for the restricted primitive model (RPM) with the MSA closure. Their model yields good results for vapor pressures and liquid densities of aqueous solutions of monovalent ions.

In this paper, we apply the PC-SAFT model to electrolyte solutions. Whereas the hard-core, dispersion, and association interactions are already taken into account by the original PC-SAFT, the charging of the ions will be considered by the Debye-Hückel (DH) term. Therefore, the DH term will first be discussed in detail.

DH Theory of Electrolyte Solutions. Debye and Hückel ${ }^{7}$ considered diluted electrolyte solutions and regarded the solvent (water) as a dielectric continuum. This assumption is feasible because the amount of water molecules is much greater than the total amount of ions. The ions are treated as hard spheres that can approach each other to distance $a_{i}$. This value is equivalent to an
ion diameter. Debye and Hückel proposed a contribution to the Helmholtz free energy for charging up a hardsphere system. This kind of model, where the ions are considered as charged spheres and the solvent is implicitly regarded as a dielectric continuum, is referred to as primitive model (PM). It should be denoted that the original DH theory was developed as a PM and not as a RPM. Both types of models consider the size of the ions, but in the RPM, the diameters of all ion species are equal.

The starting point for the treatment of charged particles is the Poisson equation:

$$
\begin{equation*}
\nabla^{2} \Phi_{\text {out }}(r)=-\rho / \epsilon \tag{3}
\end{equation*}
$$

where $\Phi_{\text {out }}(r)$ is the electric potential dependent on the distance $r$ from the center of an arbitrarily chosen ion. The index "out" denotes that eq 3 is valid outside the ionic sphere of diameter $a_{i} . \rho$ is the volumetric charge density, and $\epsilon$ is the dielectric constant of the medium. Here $\epsilon=\epsilon_{0} \epsilon_{\mathrm{r}}$, where $\epsilon_{0}$ is the permittivity in vacuo and $\epsilon_{\mathrm{r}}$ is the relative permittivity of the medium.

Following the Boltzmann principle, $\rho$ itself can be expressed by

$$
\begin{equation*}
\rho=N_{\mathrm{A}} \sum_{i} q_{i} c_{i} \exp \left(-\frac{q_{i} \Phi_{\mathrm{out}}}{k T}\right) \tag{4}
\end{equation*}
$$

where $q_{i}$ and $c_{i}$ are the charge (in coulomb) and the number concentration of ion $i$, respectively. Combining eqs 3 and 4 yields the Poisson-Boltzmann differential equation, which, in order to be easily solved, is linearized as

$$
\begin{align*}
\nabla^{2} \Phi_{\text {out }} & =\frac{N_{\mathrm{A}}}{k T \epsilon} \sum_{i} q_{i}^{2} c_{i} \Phi_{\text {out }} \\
& =\kappa^{2} \Phi_{\text {out }}, \quad \text { with } \kappa^{2}=\frac{N_{\mathrm{A}}}{k T \epsilon} \sum_{i} q_{i}^{2} c_{i} \tag{5}
\end{align*}
$$

$\kappa$ is called the inverse Debye screening length and has units of reciprocal meter. The differential eq 5 leads to the solution

$$
\begin{equation*}
\Phi_{\text {out }}(r)=\frac{A}{r} \mathrm{e}^{-\kappa r}+\frac{A^{\prime}}{r} \mathrm{e}^{\kappa r} \quad \text { for } r \geq a_{i} \tag{6}
\end{equation*}
$$

where obviously $A^{\prime}=0$ because the potential $\Phi_{\text {out }}(r)$ has to vanish for infinite distance $r$. This solution is valid for $r \geq a_{i}$. The internal region of the ionic sphere (index "in") is regarded as a continuum of given permittivity with a point charge at the center. Here, the centers of the surrounding ions are excluded and the Poisson equation (eq 3) simplifies to Laplace's equation, namely

$$
\begin{equation*}
\nabla^{2} \Phi_{\mathrm{in}}=0 \tag{7}
\end{equation*}
$$

and thus the electric potential is given by

$$
\begin{equation*}
\Phi_{\mathrm{in}}(r)=\frac{q_{i}}{4 \pi \epsilon} \frac{1}{r}+B \quad \text { for } r \leq a_{i} \tag{8}
\end{equation*}
$$

The first term of eq 8 is the self-potential often denoted as the Born potential. The two constants $A$ and $B$ are obtained by taking into consideration the boundary conditions that at $r=a_{i}$ both $\Phi_{\text {in }}$ and $\Phi_{\text {out }}$ as well as
their gradients $\partial \Phi_{\text {in }} / \partial r$ and $\partial \Phi_{\text {out }} / \partial r$ must be identical, respectively. Hence,

$$
\begin{equation*}
A=\frac{q_{i}}{4 \pi \epsilon} \frac{\mathrm{e}^{\kappa a_{i}}}{1+\kappa a_{i}} \quad \text { and } \quad B=-\frac{q_{i}}{4 \pi \epsilon} \frac{\kappa}{1+\kappa a_{i}} \tag{9}
\end{equation*}
$$

Constant $B$ represents the potential of the point charge in the center of the ion sphere. Therefore, the potential energy of one ion relative to its environment is given by

$$
\begin{equation*}
u_{i}=q_{i} B=-\frac{q_{i}^{2}}{4 \pi \epsilon} \frac{\kappa}{1+\kappa a_{i}} \tag{10}
\end{equation*}
$$

and the potential energy of the whole system becomes

$$
\begin{align*}
\frac{U^{\text {elec }}}{k T} & =\sum_{i} \frac{N_{i}}{2} \frac{u_{i}}{k T} \\
& =-\sum_{i} \frac{N_{i}}{8 \pi \epsilon} \frac{q_{i}^{2}}{k T} \frac{\kappa}{1+\kappa a_{i}} \\
& =-\frac{\kappa}{8 \pi \epsilon k T} \sum_{i} \frac{N_{i} q_{i}^{2}}{1+\kappa a_{i}} \tag{11}
\end{align*}
$$

To obtain an expression for the Helmholtz free energy, the following standard thermodynamic relationship is used:

$$
\begin{equation*}
\mathrm{d}\left(\frac{A}{T}\right)=U \mathrm{~d}\left(\frac{1}{T}\right)-\frac{P}{T} \mathrm{~d} V+\frac{1}{T} \sum_{i} \mu_{i} \mathrm{~d} n_{i} \tag{12}
\end{equation*}
$$

At constant volume and composition, the molar electrostatic Helmholtz free energy becomes ${ }^{7,43}$

$$
\begin{align*}
\frac{a^{\text {elec }}}{k T}= & -\frac{1}{4 \pi \epsilon k T} \sum_{i} \frac{x_{i} q_{i}^{2}}{3} \times \\
& \kappa\left[\frac{3}{2}+\ln \left(1+\kappa a_{i}\right)-2\left(1+\kappa a_{i}\right)+\frac{1}{2}\left(1+\kappa a_{i}\right)^{2}\right] \\
\frac{a^{\text {elec }}}{k T}= & -\frac{1}{4 \pi \epsilon k T} \sum_{i} \frac{x_{i} q_{i}^{2}}{3} \kappa \chi_{i} \tag{13}
\end{align*}
$$

where

$$
\begin{equation*}
\chi_{i}=\frac{3}{\left(\kappa a_{i}\right)^{3}}\left[\frac{3}{2}+\ln \left(1+\kappa a_{i}\right)-2\left(1+\kappa a_{i}\right)+\frac{1}{2}\left(1+\kappa a_{i}\right)^{2}\right] \tag{14}
\end{equation*}
$$

and $x_{i}$ is the mole fraction of ion $i$.
Augmented DH Theories. The original DH theory can be enhanced by considering more accurate descriptions of both the nonionic and ionic interactions. The nonelectrostatic corrections may consist of adding an attractive interaction to the hard-sphere repulsion. The electrostatic interactions may be described more accurately, for example, by using a quadratic series expansion of eq 4. These kinds of extensions to both ionic and nonionic contributions assume implicitly that the respective interactions are decoupled. In reality, this is not the case. The use of integral equation theories
provides an alternative without assuming the decoupling of these effects. ${ }^{42}$ For example, the solution of the Ornstein-Zernicke equation for charged particles with the MSA closure for the RPM and PM can be found in refs 44 and 45.

## 2. The Model

In this work, we extend the PC-SAFT EOS to account also for the interaction of charged molecules, hereafter referred to as ePC-SAFT. We restrict the grade of complexity to the PM of Debye and Hückel, keeping in mind that our long-term objective is not to optimize the EOS for electrolyte solutions but rather to find a meaningful basis for the description of aqueous protein solutions containing electrolytes. The ePC-SAFT model considers the following contributions:

$$
\begin{equation*}
a^{\mathrm{res}}=a^{\mathrm{hc}}+a^{\mathrm{disp}}(m)+a^{\mathrm{assoc}}+a^{\mathrm{elec}} \tag{15}
\end{equation*}
$$

The respective equations for the residual Helmholtz free-energy contributions to hard-chain repulsion $a^{\mathrm{hc}}$, dispersion $a^{\text {disp }}$, and association $a^{\text {assoc }}$ are summarized elsewhere. ${ }^{37}$ The contribution due to charging up the system $a^{\text {elec }}$ is calculated by eq 13. The system pressure is obtained by the standard thermodynamic relationship

$$
\begin{equation*}
\frac{p}{k T}=\frac{p^{\mathrm{id}}}{k T}-\left(\frac{\partial a^{\mathrm{res}} / k T}{\partial v}\right)_{T, \mathrm{~N}} \tag{16}
\end{equation*}
$$

The chemical potential of each component $j$ is expressed by

$$
\begin{equation*}
\frac{\mu_{j}}{k T}=\frac{\mu^{\mathrm{id}}}{k T}+\left(\frac{\partial A^{\mathrm{res}} / k T}{\partial N_{j}}\right)_{T, V, N_{i} \neq N_{j}} \tag{17}
\end{equation*}
$$

For the electrostatic contributions (elec), eqs 16 and 17 become

$$
\begin{align*}
\frac{p^{\text {elec }}}{k T}=-\left(\frac{\partial a^{\text {elec }} / k T}{\partial v}\right)_{T, \mathrm{~N}} & =-\frac{\kappa \rho_{N}}{24 \pi k T \epsilon} \sum_{k} x_{k} q_{k}^{2} \sigma_{k}  \tag{18}\\
\frac{\mu_{j}^{\text {elec }}}{k T}=\left|\frac{\partial A^{\text {elec }} / k T}{\partial N_{j}}\right|_{T, V, N_{i} \neq N_{j}} & = \\
& \left.-\frac{q_{j}^{2} \kappa}{24 \pi k T \epsilon} \left\lvert\, 2 \chi_{j}+\frac{\sum_{k} x_{k} q_{k}^{2} \sigma_{k}}{\sum_{k} x_{k} q_{k}^{2}}\right.\right) \tag{19}
\end{align*}
$$

Here $\rho_{\mathrm{N}}$ is the number density of the system and

$$
\begin{equation*}
\sigma_{k}=\left(\frac{\partial\left(\kappa \chi_{k}\right)}{\partial \kappa}\right)_{T, \mathrm{~N}}=-2 \chi_{k}+\frac{3}{1+\kappa a_{k}} \tag{20}
\end{equation*}
$$

Here, we neglect the pressure and density dependencies of the dielectric constant of water because the values of $\epsilon_{\mathrm{r}}$ in the temperature range of $278-373 \mathrm{~K}$ and at $p=1$ kPa compared to the values at $p=1 \mathrm{MPa}$ differ only at the second decimal place.

The original PC-SAFT model for associating, uncharged molecules has five parameters, namely, the

![](https://cdn.mathpix.com/cropped/52ef64d1-fd2e-401b-a5cf-6d75761c4aff-4.jpg?height=632&width=1333&top_left_y=161&top_left_x=397)
Figure 1. Liquid-phase densities (left) and vapor pressures (right) of pure water calculated with PC-SAFT (line). The crosses represent experimental data. ${ }^{47}$

Table 1. PC-SAFT Parameters for Water
| segment number | $m_{1}=1.09528$ |
| :--- | :--- |
| segment diameter $[\AA]$ | $\sigma_{1}=2.88980$ |
| dispersion energy $[\mathrm{K}]$ | $u_{1}^{0}=365.956$ |
| association sites | $N_{1}=2$ |
| association energy $[\mathrm{K}]$ | $\epsilon^{\mathrm{A}_{i} \mathrm{~B}_{j}}=2515.6706$ |
| association volume | $\kappa^{\mathrm{A}_{i} \mathrm{~B}_{j}}=0.0348679836$ |


segment number, the segment diameter, the dispersion energy, the association energy, and the association volume. The DH term does not require any additional adjustable parameter because the charge of the ions is given by their respective valence.

## 3. Results

Water. To be able to calculate the vapor pressures and densities of electrolyte solutions, the main component, namely, water, has to be modeled as accurately as possible. Gross and Sadowski ${ }^{39}$ fitted PC-SAFT model parameters to vapor-pressure and density data in a temperature range from the triple point to the critical point.

Because our long-term objective is to describe phase equilibria in aqueous protein-electrolyte solutions, the temperature range of interest is restricted to $T=278-$ 393 K . To obtain better agreement between experimental and modeled data in this temperature range, the parameter set was refitted to experimental data in this temperature range using a nonlinear least-squares algorithm. ${ }^{46}$ The parameters are summarized in Table 1. Although a four-site model would reflect best the physics of water molecules, it was shown earlier that a two-site approach yields better agreement between model and reality. ${ }^{39}$ Another advantage of using only two association sites instead of four is the decrease in computational time.

Figure 1 shows the excellent agreement between experimental and calculated data. The average deviations for vapor pressures and liquid densities are less than $0.9 \%$. The deviation of the liquid density at temperatures around the triple point is due to the wellknown anomaly of water; i.e., it exhibits a density maximum at about $4^{\circ} \mathrm{C}$, which cannot be described by any state-of-the-art model.

Electrolyte Solutions. Because there are no vaporpressure or density data available for the pure ions, their parameters must be obtained from aqueous solu-
tions. For this purpose, vapor-pressure and density data of solutions containing only one of the respective salts were considered. The investigated salts were alkali halides AnCat with $\mathrm{An}=\left\{\mathrm{Na}^{+}, \mathrm{Li}^{+}, \mathrm{K}^{+}\right\}$and $\mathrm{Cat}= \left\{\mathrm{Cl}^{-}, \mathrm{Br}^{-}, \mathrm{I}^{-}\right\}$. Further, the sulfate ion that was part of $\mathrm{Li}_{2} \mathrm{SO}_{4}, \mathrm{Na}_{2} \mathrm{SO}_{4}$, and $\mathrm{K}_{2} \mathrm{SO}_{4}$ was chosen as an example of a bivalent ion.

For the calculation of the vapor pressures and densities of electrolyte solutions, several assumptions have to be made. A reasonable approximation within the temperature range of this and our future work ( $T<333$ K ) is that the vapor phase above the solution consists of pure water only because the dissolved inorganic salts considered here are nonvolatile. A proof for the validity of this assumption is given, e.g., by Parisod and Plattner. ${ }^{48}$ At $T=653 \mathrm{~K}$, they measured a NaCl concentration of less than $0.1 \%$ in the vapor phase of a $\mathrm{NaCl}- \mathrm{H}_{2} \mathrm{O}$ solution.

The considered salts are regarded as strong electrolytes; i.e., they fully dissociate into the respective cations and anions. The ions are treated as charged hard spheres with diameter $\sigma_{j}$ ( $j$ refers to ions) that interact among each other solely by electrostatic forces. The diameter is equivalent to the distance of closest approach $a_{j}$. Dispersive interactions reign exclusively between water-water ( $u_{11}^{0}$ ) and water-ion pairs ( $u_{1 j}^{0}$ ); hence, $u_{i j}^{0}=0$ for $i \neq 1$. Association is considered only among water molecules.

In the ePC-SAFT model, the water molecules are considered explicitly in the hard-chain, the association, and the dispersion terms. The interactions arising from the polarity and polarizability of the water molecules and their effect on the electrical potential in the solution are implicitly accounted for by the dielectric constant of water. The temperature-dependent dielectric constant of water was calculated by applying the correlation of Floriano and Nascimento ${ }^{49}$ and setting the pressure to 100 kPa .

The parameters to be fitted for each ion are its hydrated diameter $\sigma_{j}$ and its dispersion energy $u_{j}^{0}$ when interacting in water. Note that $u_{j}^{0}$ is a purecomponent value of the respective ion $j$. The dispersion energy between water ( $i=1$ ) and an ion $j$ is calculated using the standard van der Waals mixing rule $u_{i j}= \left(u_{i}^{0} u_{j}^{0}\right)^{0.5}\left(1-k_{i j}\right)$.

![](https://cdn.mathpix.com/cropped/52ef64d1-fd2e-401b-a5cf-6d75761c4aff-5.jpg?height=1124&width=1397&top_left_y=161&top_left_x=363)
Figure 2. Vapor pressures of aqueous solutions of nine alkali halides at different temperatures. The lines are calculated with ePCSAFT. The crosses represent the experimental data (see Table 3 for details and references).

Table 2. Optimized ePC-SAFT Parameters for Alkali Halide Ions ${ }^{a}$
|  | $\sigma_{j}^{\mathrm{P}}[\AA]$ | $\sigma_{j}[\AA]$ | $u_{j}^{0}[\mathrm{~K}]$ |
| :--- | :--- | :--- | :--- |
| $\mathrm{Li}^{+}$ | 1.20 | 1.8059 | 1110.9261 |
| $\mathrm{Na}^{+}$ | 1.90 | 1.6262 | 119.8060 |
| $\mathrm{K}^{+}$ | 2.66 | 2.7602 | 8.8773 |
| $\mathrm{Cl}^{-}$ | 3.62 | 3.5991 | 359.6604 |
| $\mathrm{Br}^{-}$ | 3.90 | 3.8225 | 524.0636 |
| $\mathrm{I}^{-}$ | 4.32 | 4.1766 | 413.0494 |
| $\mathrm{SO}_{2}{ }^{4-}$ | 4.601 | 4.3294 | 159.9912 |


${ }^{a}$ For comparison, the Pauling ionic diameters $\sigma_{j}^{\mathrm{P}}$ are also given.

The parameters for the considered monovalent ions were fitted simultaneously to nine salt solutions, each with one of the possible cation-anion combinations (192 data points for vapor pressure and 189 data points for density). The sulfate parameters were estimated in an analogous way. Thus, the diameter and dispersion energy are determined for each ion regardless of the salt it is part of. All binary parameters $k_{i j}$ were set to zero. Nevertheless, it should be emphasized that the obtained parameters are only valid for aqueous solutions because the dielectric constant of water was used in the $a^{\text {elec }}$ contribution. Adding or exchanging the solvent will lead to a different behavior of the ions. Table 2 summarizes the optimized ion parameters.

The optimization routine yields ionic diameters that are close to but partly smaller than the respective crystal ionic (Pauling) diameters (see, e.g., ref 50). Therefore, the physical meaning of the values is somewhat unclear. In a second parameter estimation, the ion diameters were set to the Pauling values and only the dispersion energies were fitted. Using the so-obtained

Table 3. Vapor-Pressure Data ${ }^{\boldsymbol{a}}$
| salt | $T[\mathrm{~K}]$ | ref |
| :--- | :--- | :--- |
| NaCl | 278.15, 283.15, 288.15, 293.15, 298.15, 303.15, 308.15, 313.15, 318.15, 323.15 | 52 |
| NaBr | 303.15, 313.15, 323.15 | 53 |
| NaI | 303.15, 313.15, 323.15, 333.15 | 53 |
| LiCl | 303.15, 313.15, 323.15, 343.15 | 54 |
| LiBr | 303.15, 313.15, 323.15 | 54 |
| LiBr, LiI | 303.15, 313.15, 323.15, 333.15, 343.15 | 54 |
| KCl | 298.15, 303.15, 313.15 | 53, 55 |
| KBr | 291.15, 303.15, 313.15 | 53 |
| KI | 298.15, 303.15, 313.15 | 53 |
| $\mathrm{Li}_{2} \mathrm{SO}_{4}$ | 278.15, 283.15, 288.15, 293.15, 298.15, 303.15, 308.15, 313.15, 318.15, 323.15 | 56 |
| $\mathrm{Na}_{2} \mathrm{SO}_{4}$ | 293.15, 293.98, 298.15, 298.17, 301.61, 305.53, 308.15, 308.15, 308.17, 313.15, 316.03, 318.15, 323.15, 323.70, 328.15, 333.15 | 57, 58 |
| $\mathrm{K}_{2} \mathrm{SO}_{4}$ | 292.05, 295.40, 297.88, 298.73, 304.70, 309.11, 316.57, 316.57, 321.56, 325.45 | 57 |


${ }^{a}$ Temperatures and references.
parameter set, the experimental data could not be correlated well (data not shown). Thus, the parameter set was discarded.

The dispersion energies of the cations follow the correct trend; i.e., they decrease from the smallest to the largest ion. This is consistent with the fact that $\mathrm{Li}^{+}$, having the highest charge density due to the smallest diameter, interacts with more water molecules than $\mathrm{Na}^{+}$ and $\mathrm{K}^{+}$. However, the dispersion energies differ by 2 orders of magnitude. This discrepancy is not found for the dispersion energy of the anions. Nevertheless, it would be logical if $\mathrm{Cl}^{-}$had the largest value because in analogy to the $\mathrm{Li}^{+}$cation it has the highest charge density. The same problems are encountered for the dispersion energies when using the SAFT-VR EOS

![](https://cdn.mathpix.com/cropped/52ef64d1-fd2e-401b-a5cf-6d75761c4aff-6.jpg?height=1116&width=1397&top_left_y=171&top_left_x=363)
Figure 3. Liquid-phase densities of aqueous solutions of nine alkali halides at different temperatures. The lines are calculated with ePC-SAFT. The crosses represent the experimental data (taken from ref 59) at $T=293.15$ and 313.15 K .

![](https://cdn.mathpix.com/cropped/52ef64d1-fd2e-401b-a5cf-6d75761c4aff-6.jpg?height=1118&width=1409&top_left_y=1398&top_left_x=358)
Figure 4. Vapor pressures (left) and liquid-phase densities (right) of aqueous solutions of $\mathrm{Li}_{2} \mathrm{SO}_{4}, \mathrm{Na}_{2} \mathrm{SO}_{4}$, and $\mathrm{K}_{2} \mathrm{SO}_{4}$. The lines are calculated with ePC-SAFT. The crosses represent the experimental data (densities taken from ref 59 at $T=293.15 \mathrm{~K}$ for all salts and additionally at $T=333.15 \mathrm{~K}$ for $\mathrm{Li}_{2} \mathrm{SO}_{4}$ ).

Table 4. Average and Maximum Relative Errors of Vapor Pressures and Densities of Aqueous Solutions of Alkali Halides Calculated by $\boldsymbol{\Delta}=\mathbf{1}$ - Model/Expt
| salt | $\Delta p[\%]$ |  |  | $\Delta \rho[\%]$ |  |  |
| :--- | :--- | :--- | :--- | :--- | :--- | :--- |
|  | $m_{\text {max }}{ }^{a}$ | average | max | $m_{\text {max }}{ }^{a}$ | average | max |
| NaCl | 6.2 | 0.5 | 1.7 | 5.7 | 0.4 | 0.9 |
| NaBr | 8.0 | 3.6 | 7.8 | 8.0 | 0.9 | 1.7 |
| NaI | 8.4 | 0.4 | 1.9 | 10.0 | 0.3 | 0.6 |
| LiCl | 12.7 | 4.1 | 8.0 | 19.3 | 1.1 | 2.6 |
| LiBr | 13.8 | 7.3 | 18.9 | 21.4 | 1.6 | 2.7 |
| LiI | 7.5 | 1.8 | 5.1 | 13.9 | 0.6 | 1.4 |
| KCl | 4.8 | 1.8 | 7.7 | 4.5 | 0.8 | 1.8 |
| KBr | 5.0 | 3.5 | 8.5 | 5.6 | 1.0 | 1.7 |
| KI | 8.5 | 0.9 | 6.0 | 9.0 | 0.4 | 0.9 |
| $\mathrm{Li}_{2} \mathrm{SO}_{4}$ | 3.2 | 3.8 | 10.5 | 3.0 | 0.3 | 1.5 |
| $\mathrm{Na}_{2} \mathrm{SO}_{4}$ | 3.5 | 0.7 | 4.6 | 2.3 | 0.9 | 0.8 |
| $\mathrm{K}_{2} \mathrm{SO}_{4}$ | 1.0 | 0.7 | 2.2 | 0.6 | 0.2 | 0.6 |


${ }^{a} m_{\text {max }}$ denotes the maximum molality of the electrolyte solution data points for vapor pressure and density.

Table 5. Comparison of Predicted Vapor-Pressure and Experimental Data of Aqueous $\mathrm{NaCl}-\mathrm{KBr}$ and $\mathrm{NaBr}-\mathrm{KCl}$ Solutions
| $T[\mathrm{~K}]$ | $\mathrm{NaCl}+\mathrm{KBr}$ |  |  |  |  | $\mathrm{NaBr}+\mathrm{KCl}$ |  |  |  |  |
| :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- |
|  | NaCl [ $m$ ] | KBr [ $m$ ] | $p_{\exp }$ [kPa] | $p_{\text {mod }}$ [kPa] | $\Delta p$ [\%] | NaBr [ $m$ ] | KBr [ $m$ ] | $p_{\text {exp }}$ [kPa] | $p_{\text {mod }}$ [kPa] | $\Delta p$ [\%] |
| 303.15 | 0.999 | 1.001 | 3.94 | 4.02 | 2.04 | 1.000 | 1.002 | 3.94 | 4.02 | 2.03 |
| 303.15 | 1.500 | 1.498 | 3.77 | 3.93 | 4.27 | 1.502 | 1.498 | 3.77 | 3.93 | 4.27 |
| 303.15 | 2.000 | 1.999 | 3.59 | 3.85 | 7.33 | 1.999 | 2.002 | 3.60 | 3.85 | 7.02 |
| 308.15 | 1.001 | 0.998 | 5.22 | 5.31 | 1.75 | 1.002 | 0.999 | 5.21 | 5.31 | 1.95 |
| 308.15 | 1.496 | 1.502 | 4.99 | 5.19 | 4.08 | 1.501 | 1.499 | 5.01 | 5.19 | 3.66 |
| 308.15 | 2.001 | 1.998 | 4.80 | 5.09 | 6.04 | 2.001 | 1.998 | 4.79 | 5.09 | 6.27 |
| 313.15 | 0.998 | 1.001 | 6.82 | 6.95 | 1.91 | 0.999 | 0.998 | 6.83 | 6.95 | 1.77 |
| 313.15 | 1.501 | 1.498 | 6.58 | 6.80 | 3.28 | 1.498 | 1.500 | 6.60 | 6.80 | 2.97 |
| 313.15 | 1.998 | 2.001 | 6.30 | 6.66 | 5.72 | 1.999 | 1.999 | 6.32 | 6.66 | 5.39 |
| 318.15 | 0.501 | 0.499 | 9.23 | 9.24 | 0.14 | 0.501 | 0.500 | 9.23 | 9.24 | 0.14 |
| 318.15 | 1.002 | 0.999 | 8.93 | 9.01 | 0.93 | 1.003 | 0.009 | 8.91 | 9.23 | 3.65 |
| 318.15 | 1.499 | 1.501 | 8.53 | 8.81 | 3.31 | 1.497 | 1.503 | 8.55 | 8.81 | 3.07 |
| 318.15 | 2.000 | 2.001 | 8.18 | 8.64 | 5.58 | 1.999 | 2.001 | 8.20 | 8.64 | 5.31 |
| 323.15 | 0.499 | 0.498 | 11.95 | 11.88 | 0.56 | 0.499 | 0.501 | 11.94 | 11.88 | 0.49 |
| 323.15 | 0.998 | 1.004 | 11.49 | 11.59 | 0.84 | 1.002 | 0.999 | 11.49 | 11.59 | 0.84 |
| 323.15 | 1.496 | 1.501 | 11.02 | 11.33 | 2.81 | 1.499 | 1.499 | 11.00 | 11.33 | 2.99 |
| 323.15 | 1.999 | 1.998 | 10.50 | 11.10 | 5.73 | 2.003 | 1.998 | 10.49 | 11.10 | 5.83 |
| 328.15 | 0.496 | 0.503 | 15.20 | 15.15 | 0.33 | 0.501 | 0.500 | 15.20 | 15.15 | 0.34 |
| 328.15 | 1.001 | 0.999 | 14.59 | 14.77 | 1.26 | 0.997 | 1.004 | 14.60 | 14.77 | 1.18 |
| 328.15 | 2.001 | 1.998 | 14.01 | 14.44 | 3.09 | 1.501 | 1.502 | 14.00 | 14.44 | 3.16 |
| 328.15 | 1.502 | 1.497 | 13.33 | 14.15 | 6.19 | 2.002 | 1.999 | 13.33 | 14.15 | 6.18 |
| 333.15 | 0.499 | 0.502 | 19.17 | 19.17 | 0.02 | 0.503 | 0.497 | 19.17 | 19.17 | 0.02 |
| 333.15 | 1.000 | 1.001 | 18.48 | 18.69 | 1.12 | 1.002 | 1.000 | 18.48 | 18.69 | 1.11 |
| 333.15 | 1.499 | 1.502 | 17.71 | 18.27 | 3.15 | 1.499 | 1.503 | 17.71 | 18.27 | 3.15 |
| 333.15 | 1.997 | 1.998 | 16.93 | 17.91 | 5.76 | 1.998 | 2.001 | 16.91 | 17.90 | 5.88 |


approach. ${ }^{41}$ One possible explanation for this physical inconsistency might be that both the diameter and dispersion energy parameters have to compensate for the deficiency of the ePC-SAFT and SAFT-VR EOSs, respectively. They do not account for effects accompanying ion solvation such as the disruption of the water structure in the hydration sheath around the ion. Further, the polarizability of ions and water molecules, which has been shown to play an important role in the modeling of electrolyte solution properties (see, e.g., ref 51), has been neglected in the two approaches.

Experimental vapor pressures (see Table 3 for temperatures and references) of aqueous electrolyte solutions with salt molar fractions up to $x_{\mathrm{s}} \approx 0.22$ are compared to those of ePC-SAFT calculations. Figures 2 and 4 (left) show that theory and experimental data are in very good agreement over a wide range of temperatures, albeit the optimized parameters are not temper-ature-dependent. The prediction of vapor pressures for $\mathrm{Li}_{2} \mathrm{SO}_{4}$ shows higher deviations at elevated temperatures. Solution densities are also well calculated by the ePC-SAFT EOS, as shown in Figures 3 and 4 (right).

The maximum relative errors vary from $0.4 \%$ for KI to $1.6 \%$ for LiBr , which is comparable to the results presented by Galindo et al. ${ }^{41}$ Average and maximum errors for all salt solutions are tabulated in Table 4.

Although the simple DH term is expected to give very good results only at low salt concentrations ( $m<0.1$ ), the combination of the DH term with the PC-SAFT EOS allows for the prediction of solution behavior also at high molalities (also given in Table 4).

Mixed Electrolyte Solutions. In a second step, the vapor pressures of solutions containing two salts were predicted and the calculated data were compared to experimental data. ${ }^{60}$ As can be seen in Table 5, the ePCSAFT EOS can easily handle such systems without any readjustment of the ion parameters and without introduction of any binary parameters $k_{i j}$.

## 4. Conclusions

The PC-SAFT EOS has been extended to electrolyte systems by taking into account the electrostatics of ion solutions by the fairly simple DH term. The obtained ePC-SAFT EOS is able to accurately correlate and predict vapor pressures and densities of single-salt solutions as well as their mixtures up to high salt concentrations. In contrast to other electrolyte EOSs (e.g., refs 15 and 16), for ePC-SAFT two parameters are fitted for the single ions and not for the complete salt, making the model more flexible. In this work, we have put our emphasis on monovalent and bivalent strong electrolytes. The extension to protein solutions is the goal of future work.

## Acknowledgment

This work was financed by Deutsche Forschungsgemeinschaft Grant SA700/7-1. The authors thank University of Dortmund for financing a research stay of Prof. J. M. Mollerup in Dortmund with the Gambrinus Fellowship.

## Literature Cited

(1) Loehe, J. R.; Donohue, M. D. AIChE J. 1997, 43 (1), 180195.
(2) Wu, J.; Prausnitz, J. M. Fluid Phase Equilib. 1999, 155, 139-154.
(3) George, A.; Wilson, W. W. Acta Crystallogr. 1994, D50, 361365.
(4) Young, L. R. D.; Fink, A. L.; Dill, K. A. Acc. Chem. Res. 1993, 26, 614-620.
(5) Broide, M. L.; Berland, C. R.; Pande, J.; Ogun, O.; Benedek, G. B. PNAS 1991, 88, 5660-5664.
(6) Solms, N. V.; Anderson, C. O.; Blanch, H. W.; Prausnitz, J. M. AIChE J. 2002, 48 (6), 1292-1300.
(7) Debye, P.; Hückel, E. Phys. Z. 1923, 9, 185-206.
(8) Chen, C.-C.; Britt, H. I.; Boston, J. F.; Evans, L. B. AIChE J. 1982, 28, 588-596.
(9) Chen, C.-C.; Evans, L. B. AIChE J. 1986, 32, 444-454.
(10) Mock, B.; Chen, C.-C.; Evans, L. B. AIChE J. 1986, 32, 1655-1664.
(11) Chen, C.-C.; Matthias, P. M.; Orbey, H. AIChE J. 1999, 45, 1576-1586.
(12) Pérez-Villaseñor, F.; Iglesias-Silva, G.; Hall, K. R. Ind. Eng. Chem. Res. 2003, 42, 1087-1092.
(13) Xu, X.; Macedo, E. A. Ind. Eng. Chem. Res. 2003, 42, 57025707.
(14) Fürst, W.; Renon, H. Fluid Phase Equilib. 1993, 39 (2), 335-343.
(15) Wu, J.; Prausnitz, J. M. Ind. Eng. Chem. Res. 1998, 37 (5), 1634-1643.
(16) Myers, J. A.; Sandler, S. I.; Wood, R. H. Ind. Eng. Chem. Res. 2002, 41, 3282-3297.
(17) Chapman, W. G.; Gubbins, K. E.; Jackson, G.; Radosz, M. Fluid Phase Equilib. 1989, 52, 31-38.
(18) Chapman, W. G.; Gubbins, K. E.; Jackson, G.; Radosz, M. Ind. Eng. Chem. Res. 1990, 29, 1709-1721.
(19) Huang, S. H.; Radosz, M. Ind. Eng. Chem. Res. 1990, 29, 2284-2294.
(20) Huang, S. H.; Radosz, M. Ind. Eng. Chem. Res. 1991, 30, 1994-2005.
(21) Carnahan, N. F.; Starling, K. E. J. Chem. Phys. 1969, 51 (2), 635-636.
(22) Alder, B. J.; Young, D. A.; Mark, M. A. J. Chem. Phys. 1972, 56, 3013-3029.
(23) Chen, S. S.; Kreglewski, A. Ber. Bunsen. Ges. 1977, 81, 1048-1052.
(24) Wertheim, M. S. J. Stat. Phys. 1984, 35, 19-34.
(25) Wertheim, M. S. J. Stat. Phys. 1984, 35, 35-47.
(26) Wertheim, M. S. J. Stat. Phys. 1986, 42, 459-476.
(27) Wertheim, M. S. J. Stat. Phys. 1986, 42, 477-492.
(28) Song, Y.; Lambert, S. M.; Prausnitz, J. M. Ind. Eng. Chem. Res. 1994, 33, 1047-1057.
(29) Chang, J.; Sandler, S. I. Mol. Phys. 1994, 81 (3), 735-744.
(30) Gil-Villegas, A.; Galindo, A.; Whitehead, P. J.; Mills, S. J.; Jackson, G.; Burgess, A. N. J. Chem. Phys. 1997, 106, 4168-4186.
(31) Hino, T.; Prausnitz, J. M. Fluid Phase Equilib. 1997, 138, 105-130.
(32) Banaszak, M.; Chiew, Y. C.; Radosz, M. Phys. Rev. E 1993, 48, 3760-3765.
(33) Tavares, F. W.; Chang, J.; Sandler, S. I. Mol. Phys. 1995, 86 (6), 1451-1471.
(34) Chapman, W. G. J. Chem. Phys. 1990, 93, 4299-4304.
(35) Müller, E. A.; Vega, L. F.; Gubbins, K. E. Mol. Phys. 1994, 83 (6), 1209-1222.
(36) Blas, F. J.; Vega, L. F. Mol. Phys. 1997, 92, 135-150.
(37) Gross, J.; Sadowski, G. Ind. Eng. Chem. Res. 2001, 40, 1244-1260.
(38) Gross, J.; Sadowski, G. Ind. Eng. Chem. Res. 2002, 41, 1084-1093.
(39) Gross, J.; Sadowski, G. Ind. Eng. Chem. Res. 2002, 41, 5510-5515.
(40) Tumakaka, F.; Sadowski, G. Fluid Phase Equilib. 2004, 217, 233-239.
(41) Galindo, A.; Gil-Villegas, A.; Jackson, G.; Burgess, A. N. J. Phys. Chem. B 1999, 103, 10272-10281.
(42) Gil-Villegas, A.; Galindo, A.; Jackson, G. Mol. Phys. 2001, 99 (6), 531-546.
(43) Breil, M. P. Thermodynamics, Experimental, and Modelling of Aqueous Electrolyte and Amino Acid Solutions. Thesis, DTU, Lyngby, Denmark, 2001.
(44) Lee, L. L. Molecular Thermodynamics of Nonideal Fluids; Butterworth Publishers: Stoneham, MA, 1988.
(45) Qin, Y.; Prausnitz, J. M. J. Chem. Phys. 2004, 121 (7), 3181-3183.
(46) Mathworks, T. Optimization Toolbox-For Use with MAT$L A B$; The MathWorks Inc.: Natick, MA, 2004.
(47) VDI-Wärmeatlas; VDI-Gesellschaft Verfahrenstechnik und Chemieingenieurwesen (GVC): Düsseldorf, Germany, 1994.
(48) Parisod, C. J.; Plattner, E. J. Chem. Eng. Data 1981, 26, 16-20.
(49) Floriano, W. B.; Nascimento, M. A. C. Braz. J. Phys. 2004, 34 (1), 38-41.
(50) Horvath, A. L. Handbook of Aqueous Electrolyte Solutions; Ellis Horwood: Chichester, England, 1985.
(51) Jungwirth, P.; Tobias, D. J. J. Phys. Chem. B 2001, 105, 10468-10472.
(52) Apelblat, A.; Korin, E. J. Chem. Thermodyn. 1998, 30, 5971.
(53) Patil, K. R.; Tripathi, A. D.; Pathak, G.; Katti, S. S. J. Chem. Eng. Data 1991, 36, 225-230.
(54) Patil, K. R.; Tripathi, A. D.; Pathak, G.; Katti, S. S. J. Chem. Eng. Data 1990, 35, 166-168.
(55) Timmermanns, J. The Physicochemical Constants of Binary Systems in Concetrated Solutions; Interscience: New York, 1990; Vol. 3.
(56) Apelblat, A.; Korin, E. J. Chem. Thermodyn. 1998, 30, 459-471.
(57) Leopold, H. G.; Johnston, J. J. Am. Chem. Soc. 1927, 49, 1974-1988.
(58) Apelblat, A.; Korin, E. J. Chem. Thermodyn. 2002, 34, 1621-1637.
(59) Landolt-Börnstein, Numerical Data and Functional Relationships in Science and Technology, New Series, Group IV; Springer-Verlag: Berlin, 1977; Vol. 1b.
(60) Hsu, H.; Wu, Y.; Lee, L. J. Chem. Eng. Data 2003, 48, 514518.

Received for review December 6, 2004
Revised manuscript received February 9, 2005
Accepted February 25, 2005
IE0488142


[^0]:    * To whom correspondence should be addressed. Tel.: +49-231-755-2635. Fax: +49-231-755-2572. E-mail: g.sadowski@ bci.uni-dortmund.de.

