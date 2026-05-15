# Simultaneous Predictions of Chemical and Phase Equilibria in Systems with an Esterification Reaction Using PC-SAFT 

Moreno Ascani, Gabriele Sadowski (D) and Christoph Held *(D)<br>Laboratory of Thermodynamics, Department of Biochemical and Chemical Engineering, TU Dortmund University, Emil-Figge Str. 70, 44277 Dortmund, Germany<br>* Correspondence: christoph.held@tu-dortmund.de; Tel.: +49-2317552086

Citation: Ascani, M.; Sadowski, G.; Held, C. Simultaneous Predictions of Chemical and Phase Equilibria in Systems with an Esterification Reaction Using PC-SAFT. Molecules 2023, 28, 1768. https://doi.org/ 10.3390/molecules28041768

Academic Editor: Erich A. Müller

Received: 30 December 2022
Revised: 10 February 2023
Accepted: 11 February 2023
Published: 13 February 2023

Copyright: © 2023 by the authors. Licensee MDPI, Basel, Switzerland. This article is an open access article distributed under the terms and conditions of the Creative Commons Attribution (CC BY) license (https:// creativecommons.org/licenses/by/ 4.0/).


#### Abstract

The study of chemical reactions in multiple liquid phase systems is becoming more and more relevant in industry and academia. The ability to predict combined chemical and phase equilibria is interesting from a scientific point of view but is also crucial to design innovative separation processes. In this work, an algorithm to perform the combined chemical and liquid-liquid phase equilibrium calculation was implemented in the PC-SAFT framework in order to predict the thermodynamic equilibrium behavior of two multicomponent esterification systems. Esterification reactions involve hydrophobic reacting agents and water, which might cause liquid-liquid phase separation along the reaction coordinate, especially if long-chain alcoholic reactants are used. As test systems, the two quaternary esterification systems starting from the reactants acetic acid +1 -pentanol and from the reactants acetic acid + 1-hexanol were chosen. It is known that both quaternary systems exhibit composition regions of overlapped chemical and liquid-liquid equilibrium. To the best of our knowledge, this is the first time that PC-SAFT was used to calculate simultaneous chemical and liquid-liquid equilibria. All the binary subsystems were studied prior to evaluating the predictive capability of PC-SAFT toward the simultaneous chemical equilibria and phase equilibria. Overall, PC-SAFT proved its excellent capabilities toward predicting chemical equilibrium composition in the homogeneous composition range of the investigated systems as well as liquid-liquid phase behavior. This study highlights the potential of a physical sound model to perform thermodynamic-based modeling of chemical reacting systems undergoing liquid-liquid phase separation.


Keywords: liquid-liquid equilibrium; reaction equilibrium; thermodynamics; equation of state; reactive separation

## 1. Introduction

The study of chemical reactions in multiphase systems is important from a purely scientific point of view and becomes crucial in technical applications involving chemical and phase equilibria (CPE) [1]. The presence of one or more chemical reactions adds internal degrees of freedom (i.e., the species cannot only move between phases but also be converted into different components) to the multicomponent system that must be considered for designing chemical processes, thus increasing the intrinsic complexity of the system. Even for kinetically controlled systems, CPE dictates the ultimate state of the system at the given conditions $(T, p, \bar{x})[2-4]$ since the molecules that take part in the reactive systems (i.e., the reacting agents) ultimately decide on the position of the reaction equilibrium and on the heterogeneity of the reaction system. Thermodynamic analysis, performed using accurate thermodynamic models, represents a powerful tool to enable the estimation of the reaction equilibrium (i.e., the maximum yield of a chemical reaction) and to assess the possible demixing into two or more phases over the reaction coordinate. The integration of chemical reaction and phase separation into one single unit, called reactive separation, has several applications in the industry and in academia, e.g., reactive distillation [5-8], reactive extraction [9-11], reactive crystallization [12,13] or reactive absorption [14-16]. Reactive
processes may offer several advantages over their non-reactive counterparts, such as an increase in reaction yield and selectivity, the overcoming of thermodynamic restrictions (e.g., azeotropes), energy saving or capital cost reduction [17-19]. Further, reactions in multiple phases systems are encountered in living systems such as coacervates [20,21], in downstream processes for biomolecules' purification [22,23] or in geological fluid systems [24]. In process engineering, several methods have been developed so far to generate feasible flowsheets for a given separation/reaction task or to dimension and optimize the required unit operations [17], with all the proposed methods relying on the availability of experimental data. While accurate thermodynamic and kinetic data are required for a rigorous apparatus design, more qualitative data reflecting the thermodynamic behavior of the system (such as the appearance of azeotropes or the presence of a heterogeneity at given conditions, see for instance [25]) must be known at the early stage of project development. In this sense, advanced thermodynamic models offer great potential for cost and time saving since they allow existing experimental data to be correlated and their values to be predicted at experimental conditions that were not investigated [26]. As well as this correlative/extrapolative purpose, other properties that are difficult to measure (such as the interface tension [27] or diffusion behavior [28]) can be estimated using a thermodynamic model.

Motivated by the aforementioned importance, reactive phase equilibria have been the subject of many scientific works, dealing with theory and/or with experimental studies. Several authors focused on a pure theoretical-mathematical description of the reactive phase equilibrium problem, ranging from exploring the conditions of uniqueness of the phase equilibrium solution and the chemical equilibrium (CE) solution [29,30] or on its properties [31-33] to the development of new mathematical formulations of the CPE problem [34-40]. Doherty and coworkers [32,33,41,42] proposed the use of transformed variables, a concept already developed in soviet literature [43,44], to represent reactive phase diagrams and formulate working equations for CPE calculations. Later, they adapted the transformed variable to concrete calculation problems, for instance, to phase stability analysis [45] and the prediction of reactive azeotropes in a reactive mixture [46,47]. Several algorithmic approaches to solve the combined CPE problem were developed [48-57]. One successful approach is, for instance, the RAND method by Gautam et al., based on the work of White et al. [58]. The RAND algorithm is nowadays employed in commercial software such as Aspen Plus [59]. Recently, Koulocheris et al. [60] published an algorithmic approach to perform CPE of reactive VLEs using traditional $\mathrm{G}^{\mathrm{E}}$-models and tested their approach to several systems containing azeotropes. For an extensive discussion on the different available algorithms, we refer to the recent literature [53,61-63]. The group of Toikka has worked extensively in the field of CPE with theoretical considerations [64-70] and experimental/modeling contributions [71-76]. They have also published some extensive reviews on existing data for CPE systems [1,67].

The purpose of the present work is to test the thermodynamic model Perturbed Chain Statistical Associating Fluid Theory (PC-SAFT) for the modeling of CE in systems exhibiting liquid-liquid equilibria (LLE) and to develop an algorithmic approach to perform this task. PC-SAFT was chosen as it has already been used to model the CE in an esterification system of 1-butanol with acetic acid by Grob et al. [77], and by Riechert et al. [78] to model the reaction equilibrium in the esterification of ethanol and 1-propanol with acetic acid. The electrolyte version ePC-SAFT was also successfully used to model the CE of enzyme-catalyzed reactions [79,80], and the Michaelis constant [81-83], which can also be interpreted as a reaction equilibrium between free and bounded enzyme. However, none of these works considered the combined CPE problem. The proposed algorithm is a stoichiometric, equation-solving method based on a double-nested procedure with a successive update of the fugacity coefficients of all the components. The number of phases and initial composition estimates are provided by the tangent plane stability analysis [84,85]. Although the general idea of the double-nested procedure is widely known (see, for instance, [54,57,60,86-88]), in this work, we proposed some new ideas to improve the
robustness of the CPE calculation. The algorithm was applied to predict the CPE of two quaternary systems with an esterification reaction and the modeling results were compared with experimental data from the literature. Section 2.1 provides a summary of the theory of chemical reactions in multiphase equilibria, the different approaches to treat the CPE and the consequences of the occurrence of chemical reactions on the topology of phase diagrams. Details about the derivation and the structure of the proposed algorithm can be found in Section 2.2. Sections 3.1-3.5 provide a description of the investigated systems and the modeling strategy, as well as the calculated phase diagrams including the reaction equilibria. Details about PC-SAFT are summarized in Appendix A.

## 2. Algorithmic Approach

### 2.1. Thermodynamics of Chemical Reactions and Multiple Liquid Phase Equilibria

For a system at given temperature $T$, pressure $p$ and respective total moles $\bar{n}^{F}=\left(n_{1}^{F}, n_{2}^{F}, \ldots, n_{N}^{F}\right)$ of each components $i=1, \ldots, N$ at an arbitrarily chosen feed composition $F$, the mathematical solution of the CPE problem is defined by the number of phases $\pi$ and number of moles $n_{i}^{(j)}$ of each component $i$ in each phase $j$ that minimizes the total Gibbs energy of the system (Equation (1)) [4,89]:

$$
\begin{equation*}
\min _{\overline{\bar{n}}} G=\sum_{j=1}^{\pi} \sum_{i=1}^{N} n_{i}^{(j)} \mu_{i}^{(j)}=\sum_{j=1}^{\pi} \bar{n}^{(j) T} \cdot \bar{\mu}^{(j)} \tag{1}
\end{equation*}
$$

This is subject to the element conservation (Equation (2)), mass balance (Equation (3)) of each component, and non-negativity of the number of moles of each component (Equation (4)) [90].

$$
\begin{gather*}
\overline{\bar{A}} \cdot \bar{n}^{F}-\bar{b}=\overline{0}  \tag{2}\\
\bar{n}^{F}-\sum_{j=1}^{\pi} \bar{n}^{(j)}=\overline{0}  \tag{3}\\
n_{i}^{(j)} \geq 0 i=1, \ldots, N j=1, \ldots, \pi \tag{4}
\end{gather*}
$$

In Equation (2), $\overline{\bar{A}}$ is the component-element matrix, which has dimension $N^{E} \times N$ ( $N^{E}$ being the number of elements that build up the $N$ components) and whose $i$-column $\left(a_{1 i}, a_{2 i}, \ldots, a_{N^{E}}\right)^{T}$ contains the number of each element present in component $i$. The $N^{E}$-dimensional array $\bar{b}$ contains the number of moles of each element present in the system (which remains unchanged). However, among the $N^{E}$ linear equations defined by Equation (2), only $N^{C} \leq N^{E}$ equations, which are given by the rank of the matrix $\overline{\bar{A}}$, $N^{C}=\operatorname{rank}(\overline{\bar{A}})$, are sufficient to uniquely represent the condition of element conservation. The difference $N^{R}=N^{E}-N^{C}$ is denoted number of key reactions in the literature.

Based on the employed strategy to solve the CPE problem, two classes of methods can be distinguished, called non-stoichiometric and stoichiometric [38,90]. Non-stoichiometric methods directly attempt at solving the optimization problem given by Equations (1)-(4), using for instance the Kuhn-Tucker necessary conditions for minimization. Stoichiometric methods avoid the element balance constraints (Equation (2)) by formulating the component balance as a linear combination of $N^{R}$ reaction coordinates $\lambda_{j}, j=1, \ldots, N^{E}$ according to Equation (5) [60].

$$
\begin{equation*}
n_{i}^{F}=n_{i, 0}^{F}+\sum_{k=1}^{N^{R}} v_{i k} \cdot \lambda_{k} i=1, \ldots, N \bar{n}^{F}=\bar{n}_{0}^{F}+\overline{\bar{v}} \cdot \bar{\lambda} \tag{5}
\end{equation*}
$$

Here, $\bar{n}_{0}^{F}$ is a set of molar fractions satisfying the elemental abundance condition (Equation (2)) and $\overline{\bar{v}}$ being the stoichiometric matrix, which is a matrix of real numbers of dimension $N \times N^{R}$ satisfying the condition given by Equation (6) [90].

$$
\begin{equation*}
\overline{\bar{A}} \cdot \overline{\bar{v}}=\overline{\overline{0}} \tag{6}
\end{equation*}
$$

$\overline{\overline{0}}$ is a zero matrix of dimension $N^{E} \times N^{R}$. Rules for the correct choice of the number of key reactions $N^{R}$ and of the stoichiometric matrix $\overline{\bar{v}}$ is extensively discussed in the literature [3,91] and will not be repeated here. From Equations (1), (3) and (5), the necessary condition for the CPE can be derived, representing the condition that holds at equilibrium and is expressed by Equations (7) and (8) [34].

$$
\begin{gather*}
\mu_{i}^{(1)}=\mu_{i}^{(2)}=\ldots .=\mu_{i}^{(\pi)} i=1, \ldots, N  \tag{7}\\
\overline{\bar{v}} \cdot \bar{\mu}^{(R)}=\overline{0} \tag{8}
\end{gather*}
$$

where $\bar{\mu}^{(R)}$ represents the array of the chemical potentials of all the components in one chosen reference phase. Equation (8) is equivalent to the analogue reformulation based on the activity-based equilibrium constant, which is written for a reaction $k$ in Equation (9) [2,90].

$$
\begin{equation*}
\sum_{i=1}^{N} v_{i k} \cdot \mu_{i}^{(R)}=\Delta^{R} g_{k}^{0}+R T \sum_{i=1}^{N} \ln \left(x_{i}^{(R)} \gamma_{i}^{(R)}\right)^{v_{i k}}=0 \tag{9}
\end{equation*}
$$

$\Delta^{R} g_{k}^{0}$ is the standard Gibbs energy of reaction. Equation (9) summarizes the wellknown expression of the activity-based equilibrium constant $K_{a, k}$ given by Equation (10) [2,90].

$$
\begin{equation*}
K_{a, k}=\exp \left(-\frac{\Delta^{R} g_{k}^{0}}{R T}\right)=\prod_{i=1}^{N}\left(x_{i}^{(R)} \gamma_{i}^{(R)}\right)^{v_{i k}} K_{a, k}=K_{x, k}^{(R)} K_{\gamma, k}^{(R)} \tag{10}
\end{equation*}
$$

The previous conditions (Equations (7) and (8)), together with the mass balance (Equations (4) and (5)) build up a system of $\pi \cdot N+N^{R}$ equations that must be solved on $\pi \cdot N+N^{R}$ variables (which are the $\pi \cdot N$ non-negative number of moles if each component in each phase and $N^{R}$ reaction coordinates). Equations (7) and (8) build the necessary equilibrium condition of an N -component multiphase system with $N_{R}$ key reactions. It must be pointed out that there could be more solutions of Equations (7) and (8) for a different number or even for the correct number of phases $\pi$ (for instance, at the trivial solution or when an LL solution is present above the boiling point of the mixture) since the described condition is necessary but not sufficient.

For the determination of the equilibrium constant $K_{a}$ of a chemical reaction, two approaches can be employed: the first is to estimate the standard Gibbs energy of reaction $\Delta^{R} g_{k}^{0}$ from the standard energy of formation $\Delta^{F} g_{i}^{0}$ of each component $i$ in a chosen reference state. The second approach relies on the availability of at least one experimental equilibrium composition $\bar{x}^{*, \text { exp }}$ and of a thermodynamic model to predict the set of activity coefficients $\bar{\gamma}\left(T, p, \bar{x}^{*, \text { exp }}\right)$, which can then be used to determine the $K_{a, k}$ based on Equation (10).

In a system with chemical reactions, the number of degrees of freedom $F$ according to the Gibbs phase rule involving $\pi$ phases and $N$ components is reduced by the number of key reactions $N_{R}$, according to Equation (11) [64,92].

$$
\begin{equation*}
F=2-\pi+N-N_{R} \tag{11}
\end{equation*}
$$

The consequences of chemical reactions in a multiphase system for the resulting phase diagrams have been discussed extensively, for instance, in previous publications of the Toikka group [64-66]. Each key reaction decreases by one the dimension of the allowable composition space that can be reached by the system at equilibrium. For example, in a
ternary homogeneous system with a key reaction of the form $A+B \rightleftharpoons C$ at constant $T$ and $p$, the number of degrees of freedom is $F=2-1+3-1-2=1$. That is, the dimension of the composition space is reduced from two (non-reactive) to one (reactive), and thus, all the equilibrium compositions will belong to a line, called chemical equilibrium curve (CE curve) [64]. The corresponding equilibrium composition of a quaternary system with one reaction $A+B \rightleftharpoons C+D$ will span a chemical equilibrium surface (CE surface) [64]. A manifold of CE curves for different values of the equilibrium constant $K_{a}$ in a ternary ideal system are given in the ternary diagram of Figure 1.

If the system under consideration shows high non-ideality up to miscibility gap, some of the curves of the manifold can show a strong deviation from the ideal hyperbolic form shown in Figure 1. Othmer et al. [30] showed that for strongly non-ideal systems, up to three solutions can be present for some curves, although only one or two (belonging to a tie-line) corresponds to the stable solution. Figure 2 shows some CE curves, calculated for the same reaction $A+B \rightleftharpoons C$ as of Figure 1, crossing the miscibility gap of a strongly non-ideal system. Calculations are performed using the algorithm developed in Section 2.2. The tie-line that connects two points at CE is called a reactive tie-line [32], with a further distinction of a unique reactive tie-line if it appears in a ternary system [1,32,64].

![](https://cdn.mathpix.com/cropped/5459cb46-7072-4e11-a55a-d4c2d12194ce-05.jpg?height=634&width=668&top_left_y=1123&top_left_x=584)
Figure 1. CE curves for different values of the CE constant $\mathrm{K}_{\mathrm{a}}$ in a ternary system $A+B \rightleftharpoons C$ that assumes ideal mixing behavior ( $\mathrm{K}_{\gamma}=1$ in Equation (10)).

![](https://cdn.mathpix.com/cropped/5459cb46-7072-4e11-a55a-d4c2d12194ce-05.jpg?height=643&width=671&top_left_y=1932&top_left_x=584)
Figure 2. Hypothetical CE curves for different values of the CE constant $\mathrm{K}_{\mathrm{a}}$ in a strongly non-ideal ternary system $A+B \rightleftharpoons C$ with a miscibility gap. Calculations were performed using PC-SAFT with the algorithm developed in Section 2.2 and parameters listed in Section 3.2.

Whether or not phase split occurs in the reaction system is dictated by the properties of the pure components and of the mixture. The kind of molecules involved in the reaction mixture determine the equilibrium constant $K_{a}$, which in turn dictates if the reaction path undergoes phase split(s). The ternary diagram (Figure 2) shows one CE curve passing through the homogeneous region ( $K_{a}=1.8$ ), three crossing the two-phase region one time (one, however, showing one tangent point to the binodal more than crossing it) and one crossing the two-phase region twice ( $K_{a}=2.25$ ). The composition points inside the binodal will not exist in a system in equilibrium and will split into two phases (given by two points connected by a tie-line). However, the crossing points of CE curve and binodal represents points belongs to the reactive tie-lines. An exemplary diagram of a reactive system with two reactive tie-lines is shown in Figure 3, representing the real solution of the system in Figure 2 with $K_{a}=2.25$.

![](https://cdn.mathpix.com/cropped/5459cb46-7072-4e11-a55a-d4c2d12194ce-06.jpg?height=823&width=853&top_left_y=909&top_left_x=584)
Figure 3. Resulting phase diagram of the hypothetical non-ideal ternary system $A+B \rightleftharpoons C$ from Figure 2 for $K_{a}=2.25$. Visible are the two disconnected CE curves passing through the homogeneous phase and the two unique reactive tie-lines I and II.

### 2.2. Algorithm Architecture

In this section, our implementation of the algorithm to perform numerical calculation of the CPE is presented. Although the implementation is general, i.e., not limited to a single key reaction and neutral components (see, for instance, our previous work [93,94]), its scope in this work is to calculate the coupled reaction equilibrium and LLE of the two investigated esterification systems. Thus, only systems containing four molecular (i.e., non-charged) components and described by one key reaction are considered.

The phase equilibrium calculation starts with a given number of phases and an initial guess composition, which in our approach is provided by the tangent plan stability analysis [84,85]. Initially, a unimolar feed amount is initialized ( $n_{i, 0}^{F}=x_{i, 0}^{F}$ ). The number of moles of each component in each phase is thus given by Equation (12).

$$
\begin{equation*}
n_{i, 0}^{(j)}=x_{i, 0}^{(j)} \alpha_{0}^{(j)} \tag{12}
\end{equation*}
$$

Instead of searching for a direct solution of Equation (7), the objective function is reformulated as given by Equation (13) for each component.

$$
\begin{equation*}
1-\frac{x_{i}^{(j)}}{x_{i}^{\left(R_{i}\right)}} \exp \left\{\ln \varphi_{i}^{(j)}-\ln \varphi_{i}^{\left(R_{i}\right)}\right\}=0 j=1, \ldots, N_{\text {neut }} \tag{13}
\end{equation*}
$$

In Equation (13), $\varphi_{i}^{(j)}$ is the fugacity coefficient of component $i$ in phase $(j)$ and $R_{i}$ represents a reference phase, which is chosen for component $i$. Thus, $\exp \left\{\ln \varphi_{i}^{(j)}-\ln \varphi_{i}^{\left(R_{i}\right)}\right\}$ represents the partition coefficient of component $i$ between phase $j$ and the reference phase $R_{i}$. For a reactive system with $N^{R}$ key reactions, $N^{R}$ further equations must be formulated, representing the CE condition given by Equation (8). In our work, the $N^{R}$ objective functions are reformulated by Equation (14).

$$
\begin{equation*}
1-\frac{K_{\gamma, m}^{(R)} K_{x, m}^{(R)}}{K_{a, m}}=0 \tag{14}
\end{equation*}
$$

The superscript $(R)$ in Equation (14) means that the activity coefficients in $K_{\gamma, m}^{(R)}$ and the concentrations in $K_{x, m}^{(R)}$ refer to the component in only one reference phase (which do not have to be necessarily the same for all the components). Equation (13) builds up a set of $(\pi-1) \cdot N$ equations, and if $N^{R}$ key reactions must be defined, the total number of equations becomes $(\pi-1) \cdot N+N^{R}$.

Thus, the roots of the system of equations given by Equations (13) and (14) can be found by changing the value of $(\pi-1) \cdot N+N^{R}$ variables. For each component $i$, a reference phase $R_{i}$ is defined, which at the same time is used to impose the isofugacity criterion by Equation (13). The reference phase of each component is chosen as the phase in which the highest initial number of moles of component $i$ (according to Equation (11)) is present at the beginning of the CPE calculation. The number of moles of each component $i$ in the respective reference phase $R_{i}$ is found by mass balance by imposing that the total number of moles must be equal to the number of moles of the feed (as given by Equation (3)). If $N^{R}$ key reaction must be defined, then the total feed composition is corrected using the stoichiometric coefficients and $N^{R}$ reaction coordinates, which are varied by the algorithm as well. The implemented equation for the number of moles of each component in its reference phase is given by Equation (15).

$$
\begin{equation*}
n_{i}^{(R)}=\left(n_{i}^{F}-\sum_{m=1}^{N_{R}} v_{i}^{(m)} \lambda_{m}\right)-\sum_{\substack{j=1 \\ j \neq R}}^{\pi} n_{i}^{(j)} \tag{15}
\end{equation*}
$$

Thus, $(\pi-1) \cdot N+N^{R}$ variables $\left((\pi-1) \cdot N\right.$ number of moles $n_{i}^{(j)}$ of each component and $N^{R}$ reaction coordinates $\lambda_{m}, m=1, \ldots, N^{R}$ ) are varied to find the root of the system of $(\pi-1) \cdot N+N^{R}$ equations given by Equations (13) and (14).

As noted by other authors, the computationally most expensive step in phase equilibrium calculation is the evaluation of the fugacity coefficients in Equations (13) and (14) [54,87,88]. Boston et al. [87] suggested to solve the working equations in an inner loop using constant values of the fugacity coefficients and to update them, after convergence, in an outer loop. The iteration was performed until the relative difference between two successive solutions fell below a certain value. Upon using this method in our algorithm, we found convergence for the investigated esterification systems in this study. However, we observed oscillation and ultimately divergence when treating systems of concentrated electrolytes, high-pressure VLE with supercritical components and concentrated polymer solutions. In order to guarantee general robustness, we modified the double-nested approach. Instead
of working with constant fugacity coefficients in the inner loop, we approximated them as linear functions of the compositions using partial derivatives, as given by Equation (16).

$$
\begin{equation*}
\ln \varphi_{i}^{(j)}\left(T, p, \bar{x}^{(j)}\right)=\ln \varphi_{i}^{(j)}\left(T, p, \bar{x}_{0}^{(j)}\right)+\sum_{k=1}^{N} \frac{\partial \ln \varphi_{i}^{(j)}}{\partial x_{k}}\left(x_{k, 0}^{(j)}-x_{k}^{(j)}\right)_{j}^{i}=1, \ldots, N \tag{16}
\end{equation*}
$$

The partial derivatives of the fugacity coefficients in Equation (16) were evaluated by numeric differencing at the beginning of the calculation, and then updated from the previous values using a Broyden estimation [95] after each calculation step. The resulting non-linear system of equations was solved using a Newton algorithm with variable step length $\alpha$ (Equation (17)).

$$
\begin{equation*}
\bar{X}_{k+1}=\bar{X}_{k}-\alpha \cdot \overline{\bar{J}}^{-1}\left(\bar{X}_{k}\right) \cdot \bar{F}\left(\bar{X}_{k}\right) \tag{17}
\end{equation*}
$$

In Equation (17), $\bar{X}$ represents the array of $(\pi-1) \cdot N+N^{R}$ variables, $\bar{F}\left(\bar{X}_{k}\right)$ represents the array of objective functions (Equations (13) and (14)) calculated at the point $\bar{X}_{k}$ and $\overline{\bar{J}}\left(\bar{X}_{k}\right)$ is the respective Jacobian matrix. The partial derivatives in the Jacobian matrix are calculated via automatic differentiation using dual numbers [96]. The step length $\alpha$ is reduced if the new estimate $\bar{X}_{k+1}$ leads to one or more negative concentrations, or if the norm of the new objective function array is greater than the previous (i.e., if $\bar{F}\left(\bar{X}_{k+1}\right)>\bar{F}\left(\bar{X}_{k}\right)$ ). After convergence of Equations (13) and (14), the value of the fugacity coefficients and their derivatives in Equation (16) is updated and the iteration of Equation (17) is started again. The double-nested procedure is repeated until the change in the calculated mole numbers between two updates falls below a tolerance $\delta_{\text {tol }}=10^{-8}$.

## Algorithmic Structure

In following, the calculation procedure is shown using the hypothetical mixture of Figures 2 and 3 as test systems. These systems serve also as a first validation of our approach, since they show characteristic topologies of reactive phase diagrams already reported in the literature [30,66]. In sum, the algorithmic procedure to perform a reactive flash calculation according to our strategy consists of the following steps:
1- First the feed composition $\bar{x}^{F}$, the temperature $p$ and the pressure $T$ is given. Initially, an homogeneous CE calculation is performed at these conditions, according to the stoichiometry of the defined key reactions $N^{R}$. This is equivalent to moving the composition point, along a trajectory imposed by the stoichiometry called stoichiometry line, to the (hyper-)surface (composition $\bar{x}^{F \prime}$ ) where the CE condition for each key reaction is fulfilled (Equation (8)). For a simple reaction $A+B \rightleftharpoons C$ and the corresponding ternary phase diagram, this chemical equilibration step can be visualized in Figure 4.
2- Secondly, phase stability analysis according to the (non-reactive) tangent plane distance function [84,85] is performed for the chemically equilibrated feed. If the equilibrated feed lies inside the miscibility gap, two estimates of both liquid phase concentrations are provided (Figure 5).
3- Third, CE is performed for each of the single phases provided by step 2. This is equivalent to moving each single phase, according to the reaction stoichiometry, to the chemical equilibrium (hyper-)surface. The overall feed composition will move as well; however, it will in general not lie to the chemical equilibrium (hyper-)surface as with the single phases. This third step will finally provide good initial point for the final reactive flash calculation.
4- Finally, rigorous reactive flash calculation according to the strategy proposed in the last section is applied. After final convergence, two equilibrium points that satisfy Equations (7) and (8) are returned (Figure 6).

![](https://cdn.mathpix.com/cropped/5459cb46-7072-4e11-a55a-d4c2d12194ce-09.jpg?height=595&width=618&top_left_y=333&top_left_x=584)
Figure 4. The initial feed point $\bar{x}^{F}$ is moved along a stoichiometry line to the CE curve (green line). In the depicted system the final composition $\bar{x}^{F \prime}$ lies inside the miscibility gap (grey area) for the chosen initial feed composition.

![](https://cdn.mathpix.com/cropped/5459cb46-7072-4e11-a55a-d4c2d12194ce-09.jpg?height=533&width=556&top_left_y=1151&top_left_x=580)
Figure 5. After reaching the CE curve, phase stability is performed for the chemical equilibrated feed $\bar{x}^{F \prime}$. The tangent plane criterion is implemented in this work, which provides estimate of the phase compositions (grey points) if instability of the liquid phase is detected.

![](https://cdn.mathpix.com/cropped/5459cb46-7072-4e11-a55a-d4c2d12194ce-09.jpg?height=531&width=552&top_left_y=1914&top_left_x=577)
Figure 6. After convergence of the CPE calculation, the concentration of the two phases $\bar{x}^{\prime}$ and $\bar{x}^{\prime \prime}$ belonging to the same tie-line and the chemical equilibrium curve (the so-called "unique reactive tie-line") is found. The resulting feed composition still remains on the stoichiometric line (NOT on the chemical equilibrium line).

Figure 7 summarizes, in a flowchart, the computational steps of the CPE procedure explained in this section and their calling order within the implemented algorithm.

![](https://cdn.mathpix.com/cropped/5459cb46-7072-4e11-a55a-d4c2d12194ce-10.jpg?height=1040&width=883&top_left_y=461&top_left_x=575)
Figure 7. Flowchart that illustrates the computational steps of the CPE procedure implemented in the proposed algorithm. $\underline{z}^{*}$ denotes a chemical equilibrated, homogeneous feed (relevant for a feed lying outside the miscibility gap), and $\underline{x}^{(1, \text { est })}$ and $\underline{x}^{(2, \text { est })}$ denote composition estimates of the heterogeneous feed after the phase stability test.

## 3. Results

### 3.1. The Reaction Systems Considered in This Work

In this work, two esterification systems were considered to test our approach and the performance of PC-SAFT to predict, simultaneously, the occurrence of LLE along the CE. Those are the quaternary system acetic acid + 1-pentanol + pentyl acetate + water (system 1) and the system acetic acid +1 -hexanol + hexyl acetate + water (system 2 ) according to the chemical Equations (18) and (19).

$$
\begin{align*}
& \mathrm{CH}_{3} \mathrm{COOH}+\mathrm{C}_{5} \mathrm{H}_{11} \mathrm{OH} \rightleftharpoons \mathrm{CH}_{3} \mathrm{COOC}_{5} \mathrm{H}_{11}+\mathrm{H}_{2} \mathrm{O}  \tag{18}\\
& \mathrm{CH}_{3} \mathrm{COOH}+\mathrm{C}_{6} \mathrm{H}_{13} \mathrm{OH} \rightleftharpoons \mathrm{CH}_{3} \mathrm{COOC}_{6} \mathrm{H}_{13}+\mathrm{H}_{2} \mathrm{O} \tag{19}
\end{align*}
$$

The first system was characterized by Senina et al. [72] at $T=318.15 \mathrm{~K}$ and $p=1.013$ bar, the second system was extensively studied by Schmitt et al. [97-99] in a larger temperature range ( $293.15-403.15 \mathrm{~K}$ with special focus on the range $353.15-393.15 \mathrm{~K}$ ) at $p=1.013$ bar. Within the investigated temperature and pressure range, all the binary subsystems given by the alcohols and respective acetate esters with water show partial miscibility [100,101]. Thus, both the quaternary and all the ternary subsystems show a miscibility gap [ $98,102,103$ ]; the only exception is the homogenous system acetic acid + alcohol (1-pentanol or 1-hexanol) + ester (pentyl acetate or hexyl acetate). Esterification is a catalytic process and is practically frozen without a catalyst [72]: Schmitt [97] investigated the autocatalytic esterification of 1-hexanol with acetic acid at 298.15 K , showing that the reaction did not approach the equilibrium even after weeks. The reaction can be carried out in the presence of an
inorganic acid (homogeneous catalysis) or using a solid catalyst (heterogeneous catalysis). Senina et al. [72] used $\mathrm{HCl}_{\mathrm{aq}}$ in concentrations less than $2 \mathrm{wt} \%$, whereas Schmitt [97] employed an ion-exchange resin (Amberlyst CSP2). Due to the relatively low catalyst concentration, the catalyst was not considered in our calculations since it only marginally affects the phase equilibrium. In both works [72,97], the measurement of the final equilibrium composition (homogeneous CE, LLE or simultaneous CE and LLE) was carried out via gas chromatography.

### 3.2. PC-SAFT Parameters for the Considered Reaction Systems

In order to model the subsystems of the esterification systems (18) and (19), the parameters of the applied model must be determined. All the pure-component PC-SAFT parameters used in this work were retrieved from the literature and are listed in Table 1.

Table 1. PC-SAFT pure-component parameters used in this work to model the CE and LLE in the investigated systems.
| Component | $m_{i}^{\text {seg }} /-$ | $\sigma_{i} /$ | $\boldsymbol{u}_{\boldsymbol{i}} \boldsymbol{k}_{\boldsymbol{B}}^{\boldsymbol{-} \mathbf{1}} \boldsymbol{/} \mathbf{K}$ | $N_{i}$ | $\varepsilon^{A_{i} B_{i}} k_{B}^{-1} / K$ | $\boldsymbol{\kappa}^{A_{i} B_{i} /-}$ | Ref. |
| :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- |
| Water | 1.2047 | * | 353.95 | 2 | 2425.7 | 0.04509 | [104] |
| Acetic acid | 1.3402 | 3.8582 | 311.59 | 2 | 3044.4 | 0.07555 | [105] |
| 1-Pentanol | 3.6260 | 3.4508 | 247.28 | 2 | 2252.1 | 0.01033 | [105] |
| 1-Hexanol | 3.5146 | 3.6735 | 262.32 | 2 | 2538.9 | 0.00575 | [105] |
| Pentyl Acetate | 4.7077 | 3.4729 | 234.57 | 2 | 0.0 | 0.04509 | [106] |
| Hexyl Acetate | 4.8847 | 3.5834 | 241.42 | 2 | 0.0 | 0.04509 | [107] |
| $* \sigma=2.7927+\left(10.11 e^{-0.01775 T}-1.417 e^{-0.01146 T}\right)$ |  |  |  |  |  |  |  |


The binary interaction parameters used in this work were in part retrieved from the literature and, if not available, were regressed from mixture properties (LLE data in binary or ternary systems and VLE data in binary systems, see Table 2).

Table 2. Binary interaction parameters used in this work to model multicomponent mixtures using PC-SAFT. Definition of the $\mathrm{k}_{\mathrm{ij}}$ values according to the Appendix A .
| Component 1 | Component 2 | $\boldsymbol{k}_{\boldsymbol{i j}, \mathbf{2 9 8 . 1 5}} \boldsymbol{/} \boldsymbol{-}$ | $k_{i j, T} / K$ | Property Used for Estimation | Ref. |
| :--- | :--- | :--- | :--- | :--- | :--- |
| Water | Acetic acid | -0.1247 | - | VLE-binary | [107] |
| Water | 1-Pentanol | 0.001604 | 0.00016 | LLE-binary | [108] |
| Water | Pentyl Acetate | -0.0228 | - | LLE-binary | This work (using data from [100]) |
| Water | 1-Hexanol | 0.010105 | 0.000404 | LLE-binary | [108] |
| Water | Hexyl Acetate | -0.01 | 0.0015 | LLE-binary | This work (using data from [100]) |
| Acetic acid | 1-Pentanol | -0.1 | - | LLE-ternary | This work (using data from [103]) |
| Acetic acid | 1-Hexanol | -0.033 | - | LLE-ternary | This work (using data from [103]) |
| Acetic acid | Pentyl Acetate | -0.1 | - | LLE-ternary | This work (using data from [102]) |
| Acetic acid | Hexyl Acetate | -0.08 | -0.0004 | LLE-ternary | This work (using data from [98]) |
| 1-Pentanol | Pentyl Acetate | -0.0095 | - | VLE-binary | This work (using data from [109]) |
| 1-Hexanol | Hexyl Acetate | -0.0042 | - | VLE-binary | This work (using data from [98]) |


For the calculation of the CE curves in Figures 2 and 3, pure-component parameters listed in Table 3 were used. The hypothetical components are called A, B, C, as used in the calculated ternary diagrams, and the used hypothetical binary interaction parameters were chosen to: $k_{A B}=-0.045, k_{B C}=-0.025, k_{A C}=0.045$.

Table 3. PC-SAFT pure-component parameters of hypothetical mixture $\mathrm{A}+\mathrm{B}+\mathrm{C}$ used to calculate the exemplary CE curves in Figures 2 and 3.
| Component | $m_{i}^{\text {seg }} /-$ | $\sigma_{i} / \AA$ | $\boldsymbol{u}_{\boldsymbol{i}} \boldsymbol{k}_{\boldsymbol{B}}^{\boldsymbol{-} \mathbf{1}} \boldsymbol{/} \boldsymbol{K}$ | $N_{i}$ | $\varepsilon^{A_{i} B_{i}} k_{B}^{-1} / K$ | $\kappa^{A_{i} B_{i}} \boldsymbol{/}$ |
| :--- | :--- | :--- | :--- | :--- | :--- | :--- |
| A | 2.4000 | 3.2000 | 200.00 | 2 | 2500.0 | 0.05 |
| B | 1.0800 | 3.0000 | 400.00 | 2 | 2500.0 | 0.05 |
| C | 2.8000 | 3.8000 | 280.00 | 0 | - | - |


### 3.3. The Reaction Equilibrium Constants $K_{a}$ of the Considered Chemical Reactions

For the determination of the equilibrium constant $K_{a}$ of both chemical reactions (Equations (18) and (19)), we used one experimental equilibrium composition $\bar{x}^{*, \text { exp }}$ and predicted a set of activity coefficients $\bar{\gamma}\left(T, p, \bar{x}^{*, \text { exp }}\right)$ using PC-SAFT and the parameters applied in Tables 1 and 2. This was then used to determine the $K_{a, k}$ based on Equation (10). This method circumvents the approximations made in the estimation of the standard energy of formations. Senina et al. [72] measured the CE composition in the homogeneous region of system 1 as well as nine quaternary tie-lines at the CE at the given $T$ and $p$ conditions. Other experimental data [110,111] were determined in the homogeneous liquid phase but at saturation condition, i.e., along the condition of liquid-vapor coexistence. Since $T$ and $p$ at saturation vary continuously with composition, only the data of Senina et al. (determined at fixed $T$ and $p$ ) were considered in this work. For the same reasons, only CE compositions at fixed $T$ and $p$, determined by Schmitt et al. [97], were considered in this work for the determination of $K_{a}$ in system 2. The resulting $K_{a}$ values determined based on the experimental data from the literature are listed in Table 4.

Table 4. Obtained $K_{a}$ values for both systems 1 and 2, as well as the respective conditions ( $T$ and $p$ ) and the according references for the experimental equilibrium compositions.
| System | $\boldsymbol{T} / \boldsymbol{K}$ | $\boldsymbol{p} / \boldsymbol{b} \boldsymbol{a r}$ | $\boldsymbol{K}_{\boldsymbol{a}} /-$ | Ref. (for the Data) |
| :--- | :--- | :--- | :--- | :--- |
| 1 | 318.15 | 1 | 43.99 | [72] |
| 2 | 353.6 | 1 | 22.92 | [97] |


### 3.4. Prediction Results of the CPE Problem for Both Reactions under Study

Figures 8 and 9 show the CE surface in the composition tetrahedron of both systems, including the heterogeneous region of CE ("unique chemical reactive surface", according to [72]). These results were obtained using our developed algorithm (Section 2.2) fed by PC-SAFT (see the Appendix A) and the used parameters (Tables 1 and 2) for the activity coefficients as well as the equilibrium constants (Table 4). Both, the reaction surface (CE surface) and the liquid-liquid miscibility gap (binodal curve) were predicted in good agreement with the experimental data.

![](https://cdn.mathpix.com/cropped/5459cb46-7072-4e11-a55a-d4c2d12194ce-12.jpg?height=426&width=440&top_left_y=2108&top_left_x=575)
Figure 8. Quaternary phase diagram of system 2 (Equation (19)) at 353.6 K and 1 bar showing the PC-SAFT-predicted CE surface (green surface) and the PC-SAFT-predicted binodal (black curve encompassing the grey area). Experimental CE compositions of Schmitt et al. [97] are represented as grey spheres. All PC-SAFT predictions using parameters in Tables 1 and 2.

![](https://cdn.mathpix.com/cropped/5459cb46-7072-4e11-a55a-d4c2d12194ce-13.jpg?height=588&width=1210&top_left_y=338&top_left_x=577)
Figure 9. Quaternary phase diagrams of system 1 (Equation (18)) at 318.15 K and 1 bar showing the PC-SAFT-predicted CE surface (green surface) and the calculated binodal (black curve encompassing the grey area). Left: Experimental CE compositions of Senina et al. [72] (grey spheres). Right: Experimental tie-lines [72] (black spheres connected by a dashed line) and PC-SAFT-predicted tie-lines (grey spheres connected by a solid line). All PC-SAFT predictions using parameters in Tables 1 and 2.

### 3.5. Discussion

Figures 8 and 9 show that the CE composition is predicted quantitatively correct by PCSAFT, for both systems at the investigated $T$ and $p$ conditions, and in the whole composition range. System 2 shows a much broader miscibility gap in the CE surface than system 1, even at a higher temperature ( 353.6 K , compared to 318.15 K of system 1). This is in accordance with the subsystems, i.e., the much greater miscibility gap of 1-hexanol and hexyl acetate with water compared to their homologues 1-pentanol and pentyl acetate. The absence of experimental data of the CE tie-lines for system 2 did not allow for a direct comparison with the two-phase CPE prediction results in this system. Experimental CE tie-lines are available for system 1, and thus they were compared with the PC-SAFT predictions. The predicted CE tie-lines show qualitative agreement with the experimental data. It can be observed from Figure 9 that deviations between PC-SAFT and the experiments occur when acetic acid is present in the system. This inaccurate behavior of PC-SAFT at high acetic acid concentrations is already knows from previous work [77] and is probably due to the lack of representation of the dimerization behavior of acetic acid and the cross-association with the other components present in the mixture. A more detailed investigation of phase equilibria with acetic acid should be carried out in the future, trying to better capture the real association behavior of acetic acid in complex mixtures. This may require the investigation of the more refined (and likely more phenomenological) parametrization strategies of acetic acid and binary mixtures containing acetic acid.

Nevertheless, in sum, it can be concluded that the mathematical algorithm that has been developed in this work allows a satisfying estimation of chemical equilibria as well as liquid-liquid phase separation in the chemical reaction space by using PC-SAFT as the input tool for the activity coefficients. The results shown for the CE and CPE of the quaternary esterification systems are pure predictions since the model was parametrized using only pure-component vapor pressure and density (to determine the pure-component parameters, see Table 1) as well as the VLE and LLE of the binary and ternary subsystems (to determine the binary interaction parameters in Equation (A5), see Table 2). This is an important contribution to the design of reactive systems that may undergo phase separation.

## 4. Conclusions

In this work, an algorithm was successfully designed and implemented to predict CPE in multiphase multicomponent systems. New ideas were proposed to improve the
robustness of the calculation procedure when calculating the CPE of strongly non-ideal systems. The algorithm uses PC-SAFT to describe the thermodynamic behavior of the system, i.e., the fugacity coefficients of the reacting agents. Prior to modeling, the related literature on the thermodynamics of multiphase reactive systems was reviewed, and the proposed algorithm was tested against a hypothetical ternary mixture with a chemical reaction, showing that the topologies of reactive phase diagrams that are reported in the literature are also predicted well by our approach. Using the implemented algorithm, the predictive capability of PC-SAFT on the CPE could be tested successfully for the first time, against the simultaneous CE and LLE in two quaternary esterification systems, formed respectively by esterification of acetic acid and 1-pentanol and of acetic acid and 1-hexanol. The CE composition in the homogeneous phase were predicted quantitatively correct by PC-SAFT in both systems and over the whole composition range and the investigated $T$ and $p$ condition. The prediction of simultaneous CE and LLE was qualitatively correct in the whole composition range, showing higher deviations from experimental data in the presence of acetic acid. This study suggests potential improvements, possibly in a new parametrization strategy for pure and binary mixtures of acetic acid, but more importantly suggests the use of PC-SAFT to design reactive systems that may undergo phase separation.

Author Contributions: Conceptualization, M.A.; methodology, M.A.; software, M.A. and C.H.; validation, M.A.; formal analysis, M.A.; investigation, M.A.; resources, G.S. and C.H.; data curation, M.A.; writing-original draft preparation, M.A.; writing-review and editing, C.H.; visualization, M.A.; supervision, G.S. and C.H.; project administration, G.S. and C.H.; funding acquisition, G.S. and C.H. All authors have read and agreed to the published version of the manuscript.

Funding: The authors acknowledge funding from the Deutsche Forschungsgemeinschaft (DFG, German Research Foundation) under Germany's Excellence Strategy-EXC 2033-project number 390677874. Translation into German required: "Gefördert durch die Deutsche Forschungsgemeinschaft (DFG) im Rahmen der Exzellenzstrategie des Bundes und der Länder-EXC 2033-Projektnummer 390677874-RESOLV".

Institutional Review Board Statement: Not applicable.
Informed Consent Statement: Not applicable.
Data Availability Statement: The data supporting the reported results are all given in this manuscript and in Appendix A.

Acknowledgments: The authors acknowledge the work of Rinesh Vennoli.
Conflicts of Interest: The funders had no role in the design of the study; in the collection, analyses, or interpretation of data; in the writing of the manuscript; or in the decision to publish the results.

## Appendix A. PC-SAFT

PC-SAFT is a thermodynamic model belonging to the SAFT-type equations of state, and it was developed by Groß and Sadowski in 2001 [105,112]. As with all the SAFT-like equations of state, PC-SAFT bases its description of intermolecular forces on the formalism of perturbation theory [113], which allows an exact separation of each contribution of the N -body intermolecular potential $U_{N}$ into two (or more) contributions in the final expression of the residual Helmholtz energy $A^{\text {res }}$ (Equation (A1)).

$$
\begin{equation*}
U_{N}=U_{N}^{(0)}+U_{N}^{(1)} \tag{A1}
\end{equation*}
$$

For a spherical symmetric fluid, the potential separation given in Equation (A1) translates (exactly) into the residual Helmholtz energy as given by Equation (A2).

$$
\begin{equation*}
A^{r e s}=-k_{B} T \ln Z_{N}=-k_{B} T \ln Z_{N}^{(0)}\left\langle\exp \left(-\frac{U_{N}^{(1)}}{k_{B} T}\right)\right\rangle 0=A^{(0)}+A^{(1)} \tag{A2}
\end{equation*}
$$

In Equation (A2), $Z_{N}^{(0)}$ denotes the configurational integral of the unperturbed system and ${ }_{0}$ indicates a canonical average in the unperturbed system [114]. From Equation (A2), the residual Helmholtz energy can be written as the sum of a different contribution. Within this work only three terms, the hard chain, the dispersion and the association contributions, are considered to model the investigated systems (Equation (A3)).

$$
\begin{equation*}
a^{\text {res }}=\frac{A^{\text {res }}}{N k_{B} T}=\left(a^{H S}+a^{\text {Ch,for }}\right)+a^{\text {Disp }}+a^{\text {Assoc }} \tag{A3}
\end{equation*}
$$

In Equation (A3), $a^{H S}$ denotes the hard-sphere contribution described using the Carnahan-Starling equation for a mixture of hard spheres, $a^{\mathrm{Ch}, \text { for }}$ is the contribution due to chain formation from unconnected hard-sphered described by the Wertheim association theory [115-118], $a^{\text {Disp }}$ is the contribution of short-ranged, anisotropy dispersion forces developed by Groß and Sadowski [112] and $a^{A s s o c}$ is the contribution due to strong, shortranged and highly isotropic forces (such as hydrogen bonds) described by the Wertheim association theory as well. The term in brackets denotes the contribution of the hard-chain reference fluid $a^{H C}=a^{H S}+a^{C h, f o r}$ over which the dispersion contribution is developed. The main difference between PC-SAFT and the other SAFT-based equations of state is the choice of the reference fluid to develop the dispersion term $a^{\text {Disp }}$ according to Equation (A2), which is the hard-chain fluid for PC-SAFT. Therefore, information about the deviation from the symmetrical, spherical shape in real molecules is accounted for in the dispersion term. The association term $a^{\text {Assoc }}$, on the other hand, is developed as a perturbation of the hard-sphere $a^{H S}$, as in all the SAFT-like equations of state. The model needs three pure-component parameters for non-associating components and five for associating components. Those are the segment number $m_{i}^{\text {seg }}$, the segment diameter $\sigma_{i}$ and the dispersion energy $u_{i} / k_{B}$, and for an associating component with $\mathrm{N}_{\mathrm{i}}$ association sites, the association energy $\varepsilon^{A_{i} B_{i}}$ and the association volume $\kappa^{A_{i} B_{i}}$ of the association sites $A_{i}$ and $\mathrm{B}_{\mathrm{i}}$ are further needed. PC-SAFT is extended to mixtures using the Berthelot-Lorentz combining rules, cf. Equations (A4)-(A7), to estimate the interaction between two different components $i$ and $j$.

$$
\begin{gather*}
\sigma_{i j}=\frac{1}{2}\left(\sigma_{i}+\sigma_{j}\right)  \tag{A4}\\
u_{i j}=\sqrt{u_{i} u_{j}}\left(1-k_{i j}(T)\right)  \tag{A5}\\
\varepsilon^{A_{i} B_{j}}=\frac{\varepsilon^{A_{i} B_{i}}+\varepsilon^{A_{j} B_{j}}}{2}  \tag{A6}\\
\kappa^{A_{i} B_{j}}=\sqrt{\kappa^{A_{i} B_{i}} \kappa^{A_{j} B_{j}}}\left(\frac{\sqrt{\sigma_{i} \sigma_{j}}}{\frac{1}{2}\left(\sigma_{i}+\sigma_{j}\right)}\right)^{2} \tag{A7}
\end{gather*}
$$

The binary interaction parameter $k_{i j}$ is described as a linear function of the temperature according to Equation (A8).

$$
\begin{equation*}
k_{i j}(T)=k_{i j, 298.15 K}+k_{i j, T}(T / K-298.15 K) \tag{A8}
\end{equation*}
$$

The fugacity coefficient of component $i$ (at $T, p, \bar{x}$ ) is calculated from $a^{\text {res }}$, (all) the derivatives $\left(\partial a^{\text {res }} / \partial x_{i}\right)_{T, v, x_{k \neq i}}$ and the compressibility factor $Z$ using Equation (A9).

$$
\begin{equation*}
\ln \varphi_{i}=a^{r e s}+(Z-1)+\left(\frac{\partial a^{r e s}}{\partial x_{i}}\right)_{T, \rho, x_{k \neq i}}-\sum_{j=1}^{N} x_{j}\left(\frac{\partial a^{r e s}}{\partial x_{j}}\right)_{T, \rho, x_{k \neq j}}-\ln Z \tag{A9}
\end{equation*}
$$

Since $a^{\text {res }}$ and $\left(\partial a^{\text {res }} / \partial x_{i}\right)_{T, v, x_{k \neq i}}$ are the explicit function of $T, \rho, \bar{x}$ (and not $T, p, \bar{x}$ ), one needs, in advance, to calculate the number density $\rho$ corresponding to the pressure $p$, $\rho=\rho(T, p, \bar{x})$, by iteratively solving Equation (A10).

$$
\begin{equation*}
p=k_{B} T \rho\left[1+\rho\left(\frac{\partial a^{\text {res }}}{\partial \rho}\right)_{T, \rho, x_{k \neq i}}\right] \tag{A10}
\end{equation*}
$$

The whole program was written in FORTRAN. All the derivative properties in Equations (A9) and (A10) were calculated by means of automatic differentiation using dual numbers [96].

## References

1. Toikka, A.M.; Samarov, A.A.; Toikka, M.A. Phase and chemical equilibria in multicomponent fluid systems with a chemical reaction. Russ. Chem. Rev. 2015, 84, 378-392. [CrossRef]
2. Gmehling, J.; Kolbe, B. Thermodynamik, 2nd überarbeitete Auflage ed; VCH: Weinheim, Germany, 1992; ISBN 3527285474.
3. Smith, J.M.; van Ness, H.C.; Abbott, M.M. Introduction to Chemical Engineering Thermodynamics, 7th ed.; McGraw-Hill: Boston, MA, USA, 2005; ISBN 0073104450.
4. Prausnitz, J.M.; de Azevedo, E.G.; Lichtenthaler, R.N. Molecular Thermodynamics of Fluid-Phase Equilibria, 3rd ed.; Prentice Hall PTR: Upper Saddle River, NJ, USA, 1999; ISBN 0139777458.
5. Sundmacher, K.; Kienle, A. Reactive Distillation: Status and Future Directions; Wiley-VCH: Weinheim, Germany, 2003; ISBN 3527305793.
6. Górak, A.; Sorensen, E. Distillation: Fundamentals and Principles; Academic Press: Cambridge, MA, USA, 2014; ISBN 0123865484.
7. Górak, A.; Olujic, Z. Distillation: Equipment and Processes; Academic Press: Cambridge, MA, USA, 2014; ISBN 0123868793.
8. Serafimov, L.A.; Pisarenko, Y.A.; Kulov, N.N. Coupling chemical reaction with distillation: Thermodynamic analysis and practical applications. Chem. Eng. Sci. 1999, 54, 1383-1388. [CrossRef]
9. Brouwer, T.; Blahusiak, M.; Babic, K.; Schuur, B. Reactive extraction and recovery of levulinic acid, formic acid and furfural from aqueous solutions containing sulphuric acid. Sep. Purif. Technol. 2017, 185, 186-195. [CrossRef]
10. Maurer, G. Modeling the liquid-liquid equilibrium for the recovery of carboxylic acids from aqueous solutions. Fluid Phase Equilibria 2006, 241, 86-95. [CrossRef]
11. Schulz, R.; Waluga, T. Reactive extraction. In Process Intensification by Reactive and Membrane-Assisted Separations, 2nd ed.; De Gruyter: Berlin, Germany; Boston, MA, USA, 2022; ISBN 978-3-11-072045-7.
12. Berry, D.A.; Ng, K.M. Synthesis of reactive crystallization processes. AIChE J. 1997, 43, 1737-1750. [CrossRef]
13. McDonald, M.A.; Salami, H.; Harris, P.R.; Lagerman, C.E.; Yang, X.; Bommarius, A.S.; Grover, M.A.; Rousseau, R.W. Reactive crystallization: A review. React. Chem. Eng. 2021, 6, 364-400. [CrossRef]
14. Kenig, E.Y.; Schneider, R.; Górak, A. Reactive absorption: Optimal process design via optimal modelling. Chem. Eng. Sci. 2001, 56, 343-350. [CrossRef]
15. Kenig, E.Y.; Górak, A. Reactive absorption. In Integrated Chemical Processes: Synthesis, Operation, Analysis, and Control; Wiley: Hoboken, NJ, USA, 2005; pp. 265-311.
16. Kunze, A.-K. Reactive absorption. In Process Intensification by Reactive and Membrane-assisted Separations, 2nd ed.; De Gruyter: Berlin, Germany; Boston, MA, USA, 2022; ISBN 978-3-11-072045-7.
17. Skiborowski, M.; Górak, A. Hybrid separation processes. In Process Intensification by Reactive and Membrane-assisted Separations, 2nd ed.; De Gruyter: Berlin, Germany; Boston, MA, USA, 2022; ISBN 978-3-11-072045-7.
18. Schembecker, G.; Tlatlik, S. Process synthesis for reactive separations. Chem. Eng. Process. 2003, 42, 179-189. [CrossRef]
19. Malone, M.F.; Huss, R.S.; Doherty, M.F. Green chemical engineering aspects of reactive distillation. Environ. Sci. Technol. 2003, 37, 5325-5329. [CrossRef]
20. Nakashima, K.K. Chemistry of Active Coacervate Droplets: Liquid Droplets as a Minimal Model of Life. Ph.D. Thesis, Radboud University Nijmegen, Nijmegen, The Netherlands, 2021.
21. Nakashima, K.K.; Baaij, J.F.; Spruijt, E. Reversible generation of coacervate droplets in an enzymatic network. Soft Matter 2018, 14,361-367. [CrossRef]
22. Kim, Y.; Park, K.; Lee, H.; Jang, S.; Song, H.-C.; Shin, H.-C.; Park, J.J.; Park, J.; Maken, S. Purification of native and modified enzymes using a reactive aqueous two-phase system. J. Ind. Eng. Chem. 2004, 10, 384-388.
23. Campos-García, V.R.; Benavides, J.; González-Valdez, J. Reactive aqueous two-phase systems for the production and purification of PEGylated proteins. Electron. J. Biotechnol. 2021, 54, 60-68. [CrossRef]
24. Schick, D.; Bierhaus, L.; Strangmann, A.; Figiel, P.; Sadowski, G.; Held, C. Predicting $\mathrm{CO}_{2}$ solubility in aqueous and organic electrolyte solutions with ePC-SAFT advanced. Fluid Phase Equilibria 2023, 567, 113714. [CrossRef]
25. NguyenHuynh, D.; Mai, C.T.Q.; Tran, S.T.K.; Nguyen, X.T.T.; Baudouin, O. Modelling of phase behavior of ammonia and its mixtures using the mg-SAFT. Fluid Phase Equilibria 2020, 523, 112689. [CrossRef]
26. Kontogeorgis, G.M.; Folas, G.K. Thermodynamic Models for Industrial Applications: From Classical and Advanced Mixing Rules to Association Theories; John Wiley \& Sons: Hoboken, NJ, USA, 2009; ISBN 0470747544.
27. Danzer, A.; Enders, S. Comparison of two modelling approaches for the interfacial tension of binary aqueous mixtures. J. Mol. Liq. 2018, 266, 309-320. [CrossRef]
28. Borrmann, D.; Danzer, A.; Sadowski, G. Generalized Diffusion-Relaxation Model for Solvent Sorption in Polymers. Ind. Eng. Chem. Res. 2021, 60, 15766-15781. [CrossRef]
29. Caram, H.S.; Scriven, L.E. Non-unique reaction equilibria in non-ideal systems. Chem. Eng. Sci. 1976, 31, 163-168. [CrossRef]
30. Othmer, H.G. Nonuniqueness of equilibria in closed reacting systems. Chem. Eng. Sci. 1976, 31, 993-1003. [CrossRef]
31. Heidemann, R.A. Non-uniqueness in phase and reaction equilibrium computations. Chem. Eng. Sci. 1978, 33, 1517-1528. [CrossRef]
32. Ung, S.; Doherty, M.F. Theory of phase equilibria in multireaction systems. Chem. Eng. Sci. 1995, 50, 3201-3216. [CrossRef]
33. Barbosa, D.; Doherty, M.F. A new set of composition variables for the representation of reactive-phase diagrams. Proc. R. Soc. London. Ser. A Math. Phys. Sci. 1987, 413, 459-464.
34. Jiang, Y.; Smith, W.R.; Chapman, G.R. Global optimality conditions and their geometric interpretation for the chemical and phase equilibrium problem. SIAM J. Optim. 1995, 5, 813-834. [CrossRef]
35. Jiang, Y.; Chapman, G.R.; Smith, W.R. On the geometry of chemical reaction and phase equilibria. Fluid Phase Equilibria 1996, 118, 77-102. [CrossRef]
36. Smith, W.R.; Missen, R.W. Strategies for solving the chemical equilibrium problem and an efficient microcomputer-based algorithm. Can. J. Chem. Eng. 1988, 66, 591-598. [CrossRef]
37. Smith, J.V.; Missen, R.W.; Smith, W.R. General optimality criteria for multiphase multireaction chemical equilibrium. AIChE J. 1993, 39, 707-710. [CrossRef]
38. Smith, W.R. The computation of chemical equilibria in complex systems. Ind. Eng. Chem. Fundam. 1980, 19, 1-10. [CrossRef]
39. Zeleznik, F.J.; Gordon, S. Calculation of complex chemical equilibria. Ind. Eng. Chem. 1968, 60, 27-57. [CrossRef]
40. Gautam, R.; Wareck, J.S. Computation of physical and chemical equilibria-Alternate specifications. Comput. Chem. Eng. 1986, 10, 143-151. [CrossRef]
41. Barbosa, D.; Doherty, M.F. Theory of phase diagrams and azeotropic conditions for two-phase reactive systems. Proc. R. Soc. London. Ser. A. Math. Phys. Sci. 1987, 413, 443-458.
42. Barbosa, D.; Doherty, M.F. The influence of equilibrium chemical reactions on vapor-Liquid phase diagrams. Chem. Eng. Sci. 1988, 43, 529-540. [CrossRef]
43. Zharov, V.T. Open evaporation of solutions of reacting substances. Zh. Fiz. Khim 1970, 44, 1967.
44. Zharov, V.T.; Pervukhin, O.K. Structure of the Vapor-liquid Equilibrium Diagrams of Reactive Systems: II. Methanol-Formic Acid-Methyl Formate-Water System. Zh. Fiz. Khim 1972, 46, 1970.
45. Wasylkiewicz, S.K.; Ung, S. Global phase stability analysis for heterogeneous reactive mixtures and calculation of reactive liquid-liquid and vapor-liquid-liquid equilibria. Fluid Phase Equilibria 2000, 175, 253-272. [CrossRef]
46. Okasinski, M.J.; Doherty, M.F. Thermodynamic behavior of reactive azeotropes. AIChE J. 1997, 43, 2227-2238. [CrossRef]
47. Ung, S.; Doherty, M.F. Necessary and sufficient conditions for reactive azeotropes in multireaction mixtures. AIChE J. 1995, 41, 2383-2392. [CrossRef]
48. Ung, S.; Doherty, M.F. Vapor-liquid phase equilibrium in systems with multiple chemical reactions. Chem. Eng. Sci. 1995, 50, 23-48. [CrossRef]
49. McDonald, C.M.; Floudas, C.A. Global optimization for the phase and chemical equilibrium problem: Application to the NRTL equation. Comput. Chem. Eng. 1995, 19, 1111-1139. [CrossRef]
50. McDonald, C.M.; Floudas, C.A. Global optimization for the phase stability problem. AIChE J. 1995, 41, 1798-1814. [CrossRef]
51. McDonald, C.M.; Floudas, C.A. GLOPEQ: A new computational tool for the phase and chemical equilibrium problem. Comput. Chem. Eng. 1997, 21, 1-23. [CrossRef]
52. Jalali-Farahani, F.; Seader, J.D. Use of homotopy-continuation method in stability analysis of multiphase, reacting systems. Comput. Chem. Eng. 2000, 24, 1997-2008. [CrossRef]
53. Tsanas, C.; Stenby, E.H.; Yan, W. Calculation of multiphase chemical equilibrium by the modified RAND method. Ind. Eng. Chem. Res. 2017, 56, 11983-11995. [CrossRef]
54. Stateva, R.P.; Wakeham, W.A. Phase equilibrium calculations for chemically reacting systems. Ind. Eng. Chem. Res. 1997, 36, 5474-5482. [CrossRef]
55. Sanderson, R.V.; Chien, H.H. Simultaneous chemical and phase equilibrium calculation. Ind. Eng. Chem. Process Des. Dev. 1973, 12, 81-85. [CrossRef]
56. Coatléven, J.; Michel, A. A successive substitution approach with embedded phase stability for simultaneous chemical and phase equilibrium calculations. Comput. Chem. Eng. 2022, 168, 108041. [CrossRef]
57. Gupta, A.K.; Bishnoi, P.R.; Kalogerakis, N. A method for the simultaneous phase equilibria and stability calculations for multiphase reacting and non-reacting systems. Fluid Ph. Equilibria 1991, 63, 65-89. [CrossRef]
58. White, W.B.; Johnson, S.M.; Dantzig, G.B. Chemical Equilibrium in Complex Mixtures. J. Chem. Phys. 1958, 28, 751-755. [CrossRef]
59. Liu, Q.; Proust, C.; Gomez, F.; Luart, D.; Len, C. The prediction multi-phase, multi reactant equilibria by minimizing the Gibbs energy of the system: Review of available techniques and proposal of a new method based on a Monte Carlo technique. Chem. Eng. Sci. 2020, 216, 115433. [CrossRef]
60. Koulocheris, V.; Panteli, M.; Petropoulou, E.; Louli, V.; Voutsas, E. Modeling of Simultaneous Chemical and Phase Equilibria in Systems Involving Non-reactive and Reactive Azeotropes. Ind. Eng. Chem. Res. 2020, 59, 8836-8847. [CrossRef]
61. Leal, A.M.M.; Kulik, D.A.; Smith, W.R.; Saar, M.O. An overview of computational methods for chemical equilibrium and kinetic calculations for geochemical and reactive transport modeling. Pure Appl. Chem. 2017, 89, 597-643. [CrossRef]
62. Tsanas, C.; Stenby, E.H.; Yan, W. Calculation of simultaneous chemical and phase equilibrium by the method of Lagrange multipliers. Chem. Eng. Sci. 2017, 174, 112-126. [CrossRef]
63. Zhang, H. A Review on Global Optimization Methods for Phase Equilibrium Modeling and Calculations. Open Thermodyn. J. 2011, 5, 71-92. [CrossRef]
64. Toikka, A.M.; Toikka, M.A. Solubility and critical phenomena in reactive liquid-liquid systems. Pure Appl. Chem. 2009, 81, 1591-1602. [CrossRef]
65. Toikka, M.A.; Toikka, A.M. Peculiarities of phase diagrams of reactive liquid-liquid systems. Pure Appl. Chem. 2012, 85, 277-288. [CrossRef]
66. Toikka, A.M.; Toikka, M.A.; Trofimova, M.A. Chemical equilibrium in a heterogeneous fluid phase system: Thermodynamic regularities and topology of phase diagrams. Russ. Chem. Bull. 2012, 61, 741-751. [CrossRef]
67. Toikka, A.M.; Toikka, M.A.; Pisarenko, Y.A.; Serafimov, L.A. Vapor-liquid equilibria in systems with esterification reaction. Theor. Found. Chem. Eng. 2009, 43, 129-142. [CrossRef]
68. Gromov, D.; Toikka, A. Toward formal analysis of thermodynamic stability: Le Chatelier-Brown principle. Entropy 2020, 22, 1113. [CrossRef]
69. Toikka, A.M.; Jenkins, J.D. Conditions of thermodynamic equilibrium and stability as a basis for the practical calculation of vapour-liquid equilibria. Chem. Eng. J. 2002, 89, 1-27. [CrossRef]
70. Gorovits, B.I.; Toikka, A.M.; Pisarenko, Y.A.; Serafimov, L.A. Thermodynamics of heterogeneous systems with chemical interaction. Theor. Found. Chem. Eng. 2006, 40, 239-244. [CrossRef]
71. Toikka, M.A.; Kuzmenko, P.; Samarov, A.; Trofimova, M. Phase behavior of the oleic acid-methanol-methyl oleate-water mixture as a promising model system for biodiesel production: Brief data review and new results at 303.15 K and atmospheric pressure. Fuel 2022, 319, 123730. [CrossRef]
72. Senina, A.; Samarov, A.; Toikka, M.; Toikka, A. Chemical equilibria in the quaternary reactive mixtures and liquid phase splitting: A system with n-amyl acetate synthesis reaction at 318.15 K and 101.3 kPa . J. Mol. Liq. 2022, 345, 118246. [CrossRef]
73. Toikka, M.A.; Tsvetov, N.S.; Toikka, A.M. Experimental study of chemical equilibrium and vapor-liquid equilibrium calculation for chemical-equilibrium states of the n-propanol-acetic acid-n-propyl acetate-water system. Theor. Found. Chem. Eng. 2013, 47, 554-562. [CrossRef]
74. Samarov, A.; Prikhodko, I.; Shner, N.; Sadowski, G.; Held, C.; Toikka, A. Liquid-Liquid Equilibria for Separation of Alcohols from Esters Using Deep Eutectic Solvents Based on Choline Chloride: Experimental Study and Thermodynamic Modeling. J. Chem. Eng. Data 2019, 64, 6049-6059. [CrossRef]
75. Samarov, A.; Naumkin, P.; Toikka, A. Chemical equilibrium for the reactive system acetic acid+n-butanol+n-butyl acetate+ water at 308.15 K. Fluid Ph. Equilibria 2015, 403, 10-13. [CrossRef]
76. Golikova, A.; Samarov, A.; Trofimova, M.; Rabdano, S.; Toikka, M.; Pervukhin, O.; Toikka, A. Chemical equilibrium for the reacting system acetic acid-ethanol-ethyl acetate-water at $303.15 \mathrm{~K}, 313.15 \mathrm{~K}$ and 323.15 K . J. Solut. Chem. 2017, 46, 374-387. [CrossRef]
77. Grob, S.; Hasse, H. Thermodynamics of phase and chemical equilibrium in a strongly nonideal esterification system. J. Chem. Eng. Data 2005, 50, 92-101. [CrossRef]
78. Riechert, O.; Husham, M.; Sadowski, G.; Zeiner, T. Solvent effects on esterification equilibria. AIChE J. 2015, 61, 3000-3011. [CrossRef]
79. Wangler, A.; Canales, R.; Held, C.; Luong, T.Q.; Winter, R.; Zaitsau, D.H.; Verevkin, S.P.; Sadowski, G. Co-solvent effects on reaction rate and reaction equilibrium of an enzymatic peptide hydrolysis. Phys. Chem. Chem. Phys. 2018, 20, 11317-11326. [CrossRef]
80. Gajardo-Parra, N.; Akrofi-Mantey, H.O.; Ascani, M.; Cea-Klapp, E.; Garrido, J.M.; Sadowski, G.; Held, C. Osmolyte effect on enzymatic stability and reaction equilibrium of formate dehydrogenase. Phys. Chem. Chem. Phys. 2022, 24, 27930-27939. [CrossRef]
81. Wangler, A.; Böttcher, D.; Hüser, A.; Sadowski, G.; Held, C. Prediction and Experimental Validation of Co-Solvent Influence on Michaelis Constants: A Thermodynamic Activity-Based Approach. Chem. A Eur. J. 2018, 24, 16418-16425. [CrossRef]
82. Wangler, A.; Bunse, M.J.; Sadowski, G.; Held, C. Thermodynamic activity-based Michaelis constants. In Kinetics of Enzymatic Synthesis; IntechOpen: London, UK, 2018; pp. 27-49.
83. Jaworek, M.W.; Gajardo-Parra, N.F.; Sadowski, G.; Winter, R.; Held, C. Boosting the kinetic efficiency of formate dehydrogenase by combining the effects of temperature, high pressure and co-solvent mixtures. Colloids Surf. B Biointerfaces 2021, 208, 112127. [CrossRef]
84. Michelsen, M.L. The isothermal flash problem. Part I. Stability. Fluid Phase Equilibria 1982, 9, 1-19. [CrossRef]
85. Michelsen, M.L. The isothermal flash problem. Part II. Phase-split calculation. Fluid Phase Equilibria 1982, 9, 21-40. [CrossRef]
86. Alsaifi, N.M.; Englezos, P. Prediction of multiphase equilibrium using the PC-SAFT equation of state and simultaneous testing of phase stability. Fluid Ph. Equilibria 2011, 302, 169-178. [CrossRef]
87. Boston, J.F.; Britt, H.I. A radically different formulation and solution of the single-stage flash problem. Comput. Chem. Eng. 1978, 2, 109-122. [CrossRef]
88. Xiao, W.; Zhu, K.; Yuan, W.; Chien, H.H. An algorithm for simultaneous chemical and phase equilibrium calculation. AIChE J. 1989, 35, 1813-1820. [CrossRef]
89. Sandler, S.I. Chemical and Engineering Thermodynamics, 3rd ed.; Section 8.5; J. Wiley \& Sons Inc.: Hoboken, NJ, USA, 1999.
90. Smith, W.R.; Missen, R.W. Chemical Reaction Equilibrium Analisis: Theory and Algorithms; Wiley-Interscience: New York, NY, USA, 1982.
91. Ascani, M.; Held, C. Thermodynamics for reactive separations. In Process Intensification by Reactive and Membrane-Assisted Separations, 2nd ed.; De Gruyter: Berlin, Germany; Boston, MA, USA, 2022; ISBN 978-3-11-072045-7.
92. Storonkin, A.V. Thermodynamics of Heterogeneous Systems; Part 1\&2; Publishing House of Leningrad University: Leningrad, Russia, 1967.
93. Ascani, M.; Pabsch, D.; Klinksiek, M.; Gajardo-Parra, N.; Sadowski, G.; Held, C. Prediction of pH in multiphase multicomponent systems with ePC-SAFT advanced. Chem. Commun. 2022, 58, 8436-8439. [CrossRef]
94. Ascani, M.; Sadowski, G.; Held, C. Calculation of Multiphase Equilibria Containing Mixed Solvents and Mixed Electrolytes: General Formulation and Case Studies. J. Chem. Eng. Data 2022, 67, 1972-1984. [CrossRef]
95. Broyden, C.G. A class of methods for solving nonlinear simultaneous equations. Math. Comput. 1965, 19, 577-593. [CrossRef]
96. Yu, W.; Blair, M. DNAD, a simple tool for automatic differentiation of Fortran codes using dual numbers. Comput. Phys. Commun. 2013, 184, 1446-1452. [CrossRef]
97. Schmitt, M. Heterogen Katalysierte Reaktivdestillation: Stoffdaten, Experimente, Simulation und Scale-up am Beispiel der Synthese von Hexylacetat. Ph.D. Thesis, Universität Stuttgart, Stuttgard, Germany, 2006.
98. Schmitt, M.; Hasse, H. Phase equlibria for hexyl acetate reactive distillation. J. Chem. Eng. Data 2005, 50, 1677-1683. [CrossRef]
99. Schmitt, M.; Hasse, H. Chemical equilibrium and reaction kinetics of heterogeneously catalyzed n-hexyl acetate esterification. Ind. Eng. Chem. Res. 2006, 45, 4123-4132. [CrossRef]
100. Stephenson, R.; Stuart, J. Mutual binary solubilities: Water-alcohols and water-esters. J. Chem. Eng. Data 1986, 31, 56-70. [CrossRef]
101. Stephenson, R.; Stuart, J.; Tabak, M. Mutual solubility of water and aliphatic alcohols. J. Chem. Eng. Data 1984, 29, 287-290. [CrossRef]
102. Toikka, M.A.; Vernadskaya, V.; Samarov, A. Solubility, liquid-liquid equilibrium and critical states for quaternary system acetic acid-n-amyl alcohol-n-amyl acetate-water at 303.15 K and atmospheric pressure. Fluid Phase Equilibria 2018, 471, 68-73. [CrossRef]
103. Esquivel, M.M.; Bernardo-Gil, M.G. Liquid-Liquid equilibria for the systems: Water/1-pentanol/acetic acid and water/1hexanol/acetic acid. Fluid Ph. Equilibria 1991, 62, 97-107. [CrossRef]
104. Cameretti, L.F.; Sadowski, G. Modeling of aqueous amino acid and polypeptide solutions with PC-SAFT. Chem. Eng. Process. 2008, 47, 1018-1025. [CrossRef]
105. Gross, J.; Sadowski, G. Application of the perturbed-chain SAFT equation of state to associating systems. Ind. Eng. Chem. Res. 2002, 41, 5510-5515. [CrossRef]
106. Tihic, A.; Kontogeorgis, G.M.; von Solms, N.; Michelsen, M.L. Applications of the simplified perturbed-chain SAFT equation of state using an extended parameter table. Fluid Ph. Equilibria 2006, 248, 29-43. [CrossRef]
107. Pabsch, D.; Lindfeld, J.; Schwalm, J.; Strangmann, A.; Figiel, P.; Sadowski, G.; Held, C. Influence of solvent and salt on kinetics and equilibrium of esterification reactions. Chem. Eng. Sci. 2022, 263, 118046. [CrossRef]
108. Veith, H.; Voges, M.; Held, C.; Albert, J. Measuring and Predicting the Extraction Behavior of Biogenic Formic Acid in Biphasic Aqueous/Organic Reaction Mixtures. ACS Omega 2017, 2, 8982-8989. [CrossRef]
109. Gmehling, J.; Onken, U.; Arlt, W. Vapor-Liquid Equilibrium Data Collection: Organic Hydroxy Compounds: Alcohols and Phenols (Chemistry Data Series, Volume 1, Part 2b); DECHEMA Research Institute: Frankfurt, Germany, 1978.
110. Lee, L.; Liang, S. Phase and reaction equilibria of acetic acid-1-pentanol-water-n-amyl acetate system at 760 mm Hg . Fluid Ph. Equilibria 1998, 149, 57-74. [CrossRef]
111. Lee, M.-J.; Chen, S.-L.; Kang, C.-H.; Lin, H. Simultaneous chemical and phase equilibria for mixtures of acetic acid, amyl alcohol, amyl acetate, and water. Ind. Eng. Chem. Res. 2000, 39, 4383-4391. [CrossRef]
112. Gross, J.; Sadowski, G. Perturbed-chain SAFT: An equation of state based on a perturbation theory for chain molecules. Ind. Eng. Chem. Res. 2001, 40, 1244-1260. [CrossRef]
113. Zwanzig, R.W. High-temperature equation of state by a perturbation method. I. Nonpolar gases. J. Chem. Phys. 1954, 22, 1420-1426. [CrossRef]
114. McQuarrie, D.A. Statistical Mechanics; Sterling Publishing Company: New York, NY, USA, 2000; ISBN 1891389157.
115. Wertheim, M.S. Fluids with highly directional attractive forces. I. Statistical thermodynamics. J. Stat. Phys. 1984, 35, 19-34. [CrossRef]
116. Wertheim, M.S. Fluids with highly directional attractive forces. II. Thermodynamic perturbation theory and integral equations. J. Stat. Phys. 1984, 35, 35-47. [CrossRef]
117. Wertheim, M.S. Fluids with highly directional attractive forces. III. Multiple attraction sites. J. Stat. Phys. 1986, 42, 459-476. [CrossRef]
118. Wertheim, M.S. Fluids with highly directional attractive forces. IV. Equilibrium polymerization. J. Stat. Phys. 1986, 42, 477-492. [CrossRef]

Disclaimer/Publisher's Note: The statements, opinions and data contained in all publications are solely those of the individual author(s) and contributor(s) and not of MDPI and/or the editor(s). MDPI and/or the editor(s) disclaim responsibility for any injury to people or property resulting from any ideas, methods, instructions or products referred to in the content.

