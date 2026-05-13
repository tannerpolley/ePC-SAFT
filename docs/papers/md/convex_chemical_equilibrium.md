This treatment is largely derived from Section 3.5.2 from Rawlings and Ekerdt (2013) with a few modifications. Chemical equilibrium at constant temperature and pressure is achieved when the Gibbs energy is minimized. We have

$$
G\left(T, P, n_{j}\right)=\sum_{j} \mu_{j} n_{j}
$$

and can write $\mu_{j}$ in terms of a reference state and activity:

$$
\mu_{j}=\mu_{j}^{\circ}+R T \log a_{j}
$$


We express the number of moles of $j$ after reaction by

$$
n_{j}=n_{j, 0}+\sum_{i} v_{i j} \varepsilon_{i}
$$


Combining these equations, we have:

$$
\mu_{j} n_{j}=n_{j, 0} \mu_{j}^{0}+\mu_{j}^{0} \sum_{i} v_{i j} \varepsilon_{i}+\left[n_{j, 0}+\sum_{i} v_{i j} \varepsilon_{i}\right] R T \log a_{j}
$$


We define a reduced Gibbs energy for convenience:

$$
\tilde{G}=\frac{G-\sum_{j} n_{j, 0} \mu_{j}^{\circ}}{n_{T, 0} R T}=\sum_{j}\left(\sum_{i} \frac{\mu_{j}^{\circ}}{R T} v_{i j} \varepsilon_{i}^{\prime}+\left[x_{j, 0}+\sum_{i} v_{i j} \varepsilon_{i}^{\prime}\right] \log a_{j}\right)
$$

in which $\varepsilon_{i}^{\prime}=\varepsilon_{i} / n_{T, 0}$. Because $n_{j, 0}, T$, and $\mu_{j}^{\circ}$ are constant in Gibbs energy minimization, a minimizer for $\widetilde{G}$ will be a minimizer for $G$. We have that

$$
\frac{\Delta G_{\mathrm{rxn}, i}}{R T}=\sum_{j} \frac{v_{i j} \mu_{j}^{\circ}}{R T}=-\log K_{i}
$$

and therefore

$$
\tilde{G}=-\sum_{i} \log K_{i} \varepsilon_{i}^{\prime}+\sum_{j}\left[x_{j, 0}+\sum_{i} v_{i j} \varepsilon_{i}^{\prime}\right] \log a_{j}
$$


For the cases of ideal gas mixtures and ideal liquid mixtures, we can express

$$
\begin{equation*}
\log a_{j}=\log \left(y_{j} \frac{P}{P^{\circ}}\right) \quad \text { (1a) } \quad \log a_{j}=\log x_{j} \tag{1b}
\end{equation*}
$$


For now, assume we have an ideal liquid mixture. Define:

$$
\tilde{x}_{j}=x_{j, 0}+\sum_{i} v_{i j} \varepsilon_{i}^{\prime}=\frac{n_{j}}{n_{T, 0}}
$$


Then

$$
\log a_{j}=\log x_{j}=\log \left(\frac{\tilde{x}_{j}}{\tilde{n}_{T}}\right)=\log \tilde{x}_{j}-\log \tilde{n}_{T}
$$

in which $\tilde{n}_{T}=n_{T} / n_{T, 0}=\sum_{j} \tilde{x}_{j}$. We can thus write

$$
\tilde{G}=-\sum_{i} \log K_{i} \varepsilon_{i}^{\prime}+\sum_{j} \tilde{x}_{j} \log \tilde{x}_{j}-\log \tilde{n}_{T} \sum_{j} \tilde{x}_{j}=-\sum_{i} \log K_{i} \varepsilon_{i}^{\prime}+\sum_{j} \tilde{x}_{j} \log \tilde{x}_{j}-\tilde{n}_{\mathrm{T}} \log \tilde{n}_{T}
$$


This final expression is amenable for optimization. The extents of reaction enter the objective function linearly, and $x \log x$ is a strictly convex function. The only problematic term is $\tilde{n}_{\mathrm{T}} \log \tilde{n}_{T}$, which is negative and therefore concave. However, it can be show that if we treat $\tilde{n}_{T}$ as an Expression and apply the chain rule, the resulting Hessian is positive semidefinite (see Rawlings and Ekerdt (2013), Exercise 3.14). Furthermore, if $\tilde{x}_{j}$ is eliminated in favor of $\varepsilon_{i}^{\prime}$, the Hessian is positive definite. However, in formulating this problem for solution in Pyomo and IPOPT, it is advantageous to keep $\tilde{x}_{j}$ a Var because then we can give it a lower bound of 0 . The entropy of mixing term acts effectively as a log barrier term, so strictly speaking that bound is unnecessary, but 1) it is useful to keep them as separate variables in order to discover values of $\varepsilon_{i}^{\prime}$ that lead to $\tilde{x}_{j}>0$ for all $j$ and 2 ) the solver $\log$ will not get flooded with AMPL errors when values of $\varepsilon_{i}^{\prime}$ with some $\tilde{x}_{j} \leq 0$ are used.

So the convex optimization problem solved is then

$$
\begin{aligned}
& \min _{\varepsilon_{i}^{\prime}, \tilde{x}_{j}}-\sum_{i} \log K_{i} \varepsilon_{i}^{\prime}+\sum_{j} \tilde{x}_{j} \log \tilde{x}_{j}-\left(\sum_{j} \tilde{x}_{j}\right) \log \left(\sum_{j} \tilde{x}_{j}\right) \\
& \text { subject to } \quad \tilde{x}_{j}=x_{j, 0}+\sum_{i} v_{i j} \varepsilon_{i}^{\prime} \\
& \quad \tilde{x}_{j} \geq 0
\end{aligned}
$$


It's worth going through the same exercise for an ideal vapor phase. We have

$$
\log a_{j}=\log \left(\boldsymbol{y}_{j} \frac{P}{P^{\circ}}\right)=\log \tilde{y}_{j}-\log \tilde{n}_{j}+\log \frac{P}{P^{\circ}}
$$

and therefore

$$
\begin{aligned}
\tilde{G}=-\sum_{i} \log K_{i} & \varepsilon_{i}^{\prime}+\sum_{j} \tilde{y}_{j} \log \tilde{y}_{j}+\left(-\log \tilde{n}_{T}+\log \frac{P}{P^{\circ}}\right) \sum_{j} \tilde{x}_{j} \\
& =-\sum_{i} \log K_{i} \varepsilon_{i}^{\prime}+\sum_{j} \tilde{x}_{j} \log \tilde{x}_{j}-\tilde{n}_{\mathrm{T}} \log \tilde{n}_{T}+\tilde{n}_{\mathrm{T}} \log \frac{P}{P^{\circ}}
\end{aligned}
$$


The addition of a term linear in $\tilde{n}_{\mathrm{T}}$ does not affect the convexity of the Hessian, so the resulting problem remains convex.