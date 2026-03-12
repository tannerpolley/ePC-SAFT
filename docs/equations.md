# Residual Helmholz Energy

$$\begin{aligned}
    \tilde{a}^{\mathrm{res}}=\tilde{a}^{h c}+\tilde{a}^{\text {disp }}+\tilde{a}^{\text {assoc }}+\tilde{a}^{\text {DH }}+\tilde{a}^{\text {Born }}
\end{aligned}$$

## Hard-Chain Reference Contribution

$$\tilde{a}^{\mathrm{hc}}=\bar{m} \tilde{a}^{\mathrm{hs}}-\sum_i x_i\left(m_i-1\right) \ln \mathrm{g}_{i i}^{\mathrm{hs}}\left(\sigma_{i i}\right)$$

$$\bar{m}=\sum_i x_i m_i$$

$$\begin{aligned}
\tilde{a}^{\mathrm{hs}} =\frac{1}{\zeta_0}\left[\frac{3 \zeta_1 \zeta_2}{\left(1-\zeta_3\right)}+\frac{\zeta_2^3}{\zeta_3\left(1-\zeta_3\right)^2}+\left(\frac{\zeta_2^3}{\zeta_3^2}-\zeta_0\right) \ln \left(1-\zeta_3\right)\right]
\end{aligned}$$

$$\begin{aligned}
& \mathrm{g}_{i j}^{\mathrm{hs}}=\frac{1}{\left(1-\zeta_3\right)}+\left(\frac{d_i d_j}{d_i+d_j}\right) \frac{3 \zeta_2}{\left(1-\zeta_3\right)^2} + \left(\frac{d_i d_j}{d_i+d_j}\right)^2 \frac{2 \xi_2{ }^2}{\left(1-\xi_3\right)^3}
\end{aligned}$$

$$\xi_n=\frac{\pi}{6} \rho \sum_i x_i m_i d_i^n \quad n \in\{0,1,2,3\}$$

$$d_i=\sigma_i\left[1-0.12 \exp \left(-3 \frac{\epsilon_i}{k T}\right)\right]$$

## Dispersion Contribution

$$\tilde{a}^{\mathrm{disp}}=-2 \pi \rho I_1(\eta, \bar{m}) \overline{m^2 \epsilon \sigma^3}-\pi \rho \bar{m} C_1 I_2(\eta, \bar{m}) \overline{m^2 \epsilon^2 \sigma^3}$$

$$C_1 = \left(1+\bar{m} \frac{8 \eta-2 \eta^2}{(1-\eta)^4}+\right.\left.\quad(1-\bar{m}) \frac{20 \eta-27 \eta^2+12 \eta^3-2 \eta^4}{[(1-\eta)(2-\eta)]^2}\right)$$

$$\overline{m^2 \epsilon \sigma^3}=\sum_i \sum_j x_i x_j m_i m_j\left(\frac{\epsilon_{i j}}{k T}\right) \sigma_{i j}^3$$

$$\overline{m^2 \epsilon^2 \sigma^3}=\sum_i \sum_j x_i x_j m_i m_j\left(\frac{\epsilon_{i j}}{k T}\right)^2 \sigma_{i j}{ }^3$$

$$\sigma_{i j}=\frac{1}{2}\left(\sigma_i+\sigma_j\right) \cdot \left(1- l_{ij} \right)$$

$$\epsilon_{ij}=
\begin{cases}
\sqrt{\epsilon_i\,\epsilon_j}\left(1-k_{ij}\right), & z_i z_j \neq 0,\\[6pt]
0, & z_i z_j = 0 .
\end{cases}$$

$$\epsilon_{i j}= 0$$

$$I_1(\eta, \bar{m})=\sum_{i=0}^6 a_i(\bar{m}) \eta^i$$

$$I_2(\eta, \bar{m})=\sum_{i=0}^6 b_i(\bar{m}) \eta^i$$

$$a_i(\bar{m})=a_{0 i}+\frac{\bar{m}-1}{\bar{m}} a_{1 i}+\frac{\bar{m}-1}{\bar{m}} \frac{\bar{m}-2}{\bar{m}} a_{2 i}$$

$$b_i(\bar{m})=b_{0 i}+\frac{\bar{m}-1}{\bar{m}} b_{1 i}+\frac{\bar{m}-1}{\bar{m}} \frac{\bar{m}-2}{\bar{m}} b_{2 i}$$

$$\rho=\frac{6}{\pi} \eta\left(\sum_i x_i m_i d_i^3\right)^{-1}$$

$$\hat{\rho}=\frac{\rho}{N_{\mathrm{AV}}}\left(10^{10} \frac{A^\circ{}}{\mathrm{~m}}\right)^3\left(10^{-3} \frac{\mathrm{kmol}}{\mathrm{~mol}}\right)$$

## Association Contribution

$$\tilde{a}^{\mathrm{assoc}}= \sum_i x_i\sum_{\mathrm{A}_i}\left(\ln X^{\mathrm{A}_i}-\frac{X^{\mathrm{A}_i}}{2}+\frac{1}{2}\right)$$

$$\begin{gathered}
X^{\mathrm{A}_i} =\left[1+ \sum_j \sum_{\mathrm{B}_j} \rho x_jX^{\mathrm{B}_j} \Delta^{\mathrm{A}_i \mathrm{~B}_j}\right]^{-1}  \\ \sum_{\mathrm{B}_j}  \text{over ALL sites on molecule j},  \mathrm{~A}_j, \mathrm{~B}_j, \mathrm{C}_j, \ldots ; \sum_j \text{over all components})
\end{gathered}$$

$$\rho_j=X_j \rho_{\text {mixture }}$$

$$\Delta^{\mathrm{A}, \mathrm{~B}_j}=d_{i j}{ }^3 \mathrm{g}_{i j}^{\mathrm{hs}}\left[\exp \left(\frac{\epsilon^{\mathrm{A}, \mathrm{~B}_j}}{k T}\right)-1\right]$$

$$\varepsilon^{A_{i}B_{j}}=\frac{1}{2}(\varepsilon^{A_{i}B_{i}}+\varepsilon^{A_{j}B_{j}})(1-k_{ij}^{\mathrm{hb}})$$

$$k^{A_iB_j}=\sqrt{k^{A_iB_i}k^{A_jB_j}}\quad\left(\frac{\sqrt{\sigma_i\sigma_j}}{1/2(\sigma_i+\sigma_j)}\right)^3$$

$$d_{i j}=\left(d_{i i}+d_{j j}\right) / 2 .$$

## Debye and Huckel Electrolyte Term Contribution

$$\tilde{a}^{DH}=-\frac{\kappa e^{2}}{12\pi\varepsilon_{0}\varepsilon_{r}k_{B}T}\sum_{i}x_{i}z_{i}^{2}\chi_{i}$$

$$\kappa=\sqrt{\frac{\rho e^{2}}{k_{B}T\varepsilon_{0}\varepsilon_{r}}\sum_{j}x_{j}z_{j}^{2}}$$

$$\chi_i=\frac{3}{\left(\kappa d_i\right)^3}\left[\frac{3}{2}+ln(1+\kappa d_i)-2(1+\kappa d_i)+\frac{1}{2}\left(1+\kappa d_i\right)^2\right]$$

### Ion Diameter

Temperature Independent Version
$$d_i=\sigma_{i} \qquad i\in\mathcal{I}$$

Temperature Dependent Version
$$d_i=\left[1-0.12\right]\sigma_{i} = .88\sigma_{i}  \qquad i\in\mathcal{I}$$

Temperature Dependent Version 2
$$d_i=\left[1-0.12 \exp \left(-3 \frac{\epsilon_{i}}{k T}\right)\right]\sigma_{i}  \qquad i\in\mathcal{I}$$

### Dielectric Constant

Mole-Fraction
$$\varepsilon_{r}=\sum_{j=1}^{N_c}\varepsilon_{r,j}x_j$$

Mass-Fraction
$$\varepsilon_r=\sum_{j=1}^{N_c}\varepsilon_{r,j}\,w_j
=
\frac{\sum_{j=1}^{N_c}\varepsilon_{r,j}\,x_j\,MW_j}{\sum_{j=1}^{N_c}x_j\,MW_j}
=
\frac{\sum_{j=1}^{N_c}\varepsilon_{r,j}\,x_j\,MW_j}{\overline{MW}}$$

Combo
$$\varepsilon_r = x_{sol}\,\varepsilon_{r,sol}^{(w)} + \sum_{j\in\mathcal I}\varepsilon_{r,j}\,x_j$$

$$\varepsilon_{r,sol}^{(w)} \equiv \sum_{j\in\mathcal S}\varepsilon_{r,j}\,w_j^{sol} = \frac{\sum_{j\in\mathcal S}\varepsilon_{r,j}\,x_j\,MW_j}{\sum_{j\in\mathcal S}x_j\,MW_j} = \frac{\sum_{j\in\mathcal S}\varepsilon_{r,j}\,x_j\,MW_j}{\overline{MW}_{sol}}$$

New Mixing Rule (2025)

$$\varepsilon_{r,\mathrm{mix}}=\frac{\varepsilon_{r,\mathrm{solvent,mix}}^{\mathrm{salt-free}}}{1+7.01\cdot x_{\mathrm{ion}}}$$

$$\varepsilon_{r,\mathrm{solvent,mix}}^{\mathrm{salt-free}}=\sum_{\mathrm{solvent}}w_{\mathrm{solvent}}^{\mathrm{salt-free}}\cdot\varepsilon_{r,\mathrm{solvent}}$$

Intermediates

$$\overline{MW}=\sum_{j=1}^{N_c} x_j\,MW_j$$

$$\overline{MW}_{sol}=\sum_{j\in\mathcal S}x_j\,MW_j$$

$$\mathcal S=\{j:\,z_j=0\},\qquad \mathcal I=\{j:\,z_j\neq 0\}$$

$$x_{sol}=\sum_{j\in\mathcal S}x_j$$

### Bjerrum Treatment

$$\tilde{a}^{DH}=-\frac{\kappa e^{2}}{12\pi\varepsilon_{0}\varepsilon_{r}k_{B}T}\sum_{i}\alpha_{i}x_{i}z_{i}^{2}\chi_{i}$$

$$\kappa=\sqrt{\frac{\rho e^{2}}{k_{B}T\varepsilon_{0}\varepsilon_{r}}\sum_{j}\alpha_{j}x_{j}z_{j}^{2}}$$

$$\chi_i=\frac{3}{\left(\kappa R_i\right)^3}\left[\frac{3}{2}+ln(1+\kappa R_i)-2(1+\kappa R_i)+\frac{1}{2}\left(1+\kappa R_i\right)^2\right]$$

#### Bjerrum Length

$$R_i=\begin{array}{c}\{\begin{array}{ccc}d_{ion}&if&d_{ion}>l_B\\l_B&if&d_{ion}<l_B\end{array}\end{array}$$

$$l_B=\frac{\left|z_iz_j\right|e^2}{8\pi\varepsilon_0\varepsilon_rk_BT}$$

#### Dissociation Degree

Solved implicitly

$$\alpha=\frac{-1+\sqrt{1+4x_\pm K_{ip}\frac{\left(\gamma_\pm^*\left(x_f\right)\right)^2}{\gamma_{ip}^*}}}{2x_\pm K_{ip}\frac{\left(\gamma_\pm^*\left(x_f\right)\right)^2}{\gamma_{ip}^*}}$$

$$K_{ip}(T)=\frac{(1-\alpha)\cdot\gamma_{ip}^{*}}{\alpha^{2}\cdot x_{\pm}\cdot\left(\gamma_{\pm}^{*}(x_{f}(\alpha)\right)^{2}} = 4\pi\rho_N\overset{l_B}{\operatorname*{\int_a^{l_B}}}\exp\left(\frac{\left|z_iz_j\right|e^2}{4\pi\varepsilon_0\varepsilon_rk_BT}\cdot\frac{1}{r}\right)r^2dr$$

$$\gamma_{ip}^* \approx 1$$

$$x_\pm = x_f + x_{ip}$$

$$x_\pm = \alpha x_f$$

$$x_\pm = (1 - \alpha) x_{ip}$$

$$x_{\pm}=(x_{c}^{\nu c}\cdot x_{a}^{\nu a})^{\frac{1}{\nu_{c}+\nu_{a}}}$$

$$\mu_{i}^{DH}\left(x_{f}\right)=\left(\frac{\partial A^{DH}}{\partial\rho_{i}\left(x_{f}\right)}\right)_{T,V,N_{j\neq i}}=-\frac{e^{2}z_{i}^{2}\kappa}{24\pi\varepsilon_{0}\varepsilon_{r}}\left[2\chi_{i}+\frac{\sum_{k}x_{k,f}z_{k}^{2}\sigma_{k}}{\sum_{k}x_{k,f}z_{k}^{2}}\right]$$

$$ln\gamma_{i}^{*}\left(x_{f,i}\right)=-\frac{e^{2}z_{i}^{2}\kappa}{24\pi\varepsilon_{0}\varepsilon_{r}}\left[2\chi_{i}+\frac{\sum_{k}x_{k,f}z_{k}^{2}\sigma_{k}}{\sum_{k}x_{k,f}z_{k}^{2}}\right]$$

$$\gamma_{\pm}^{*}=(\gamma_{c}^{*,\nu c}\cdot\gamma_{a}^{*,\nu a})^{\frac{1}{\nu_{c}+\nu_{a}}}$$

## Born Electrolyte Term Contribution

$$\tilde{a}^{\text {Born }}(\varepsilon(x))=-\frac{e^2}{4 \pi \varepsilon_0 k_B T}\left(1-\frac{1}{\varepsilon_{r,\mathrm{bulk}}}\right) \sum_i \frac{x_i z_i^2}{d^{\text{Born}}_{i}}$$

$$d^{\text{Born}}_{i} = d_{i} \qquad i\in\mathcal{I}$$

### Solvation Shell Model + Dielectric Saturation

Version 1a: Equation taken directly from the paper (probably wrong):

$$a_{\mathrm{SSM+DS}}^{\mathrm{Born}}=-\frac{e^{2}}{4\pi\varepsilon_{0}k_{\mathrm{B}}T}\left(1-\frac{1}{\varepsilon_{r,\mathrm{bulk}}}\right)\sum_{i}\frac{x_{i}z_{i}^{2}}{d_{i}^{\mathrm{Born}}+\Delta d_{i}}+\left(1-\frac{1}{\varepsilon_{r,\mathrm{bulk}}}\right)\left(\sum_{i}x_{i}z_{i}^{2}\left(\frac{1}{d_{i}^{\mathrm{Born}}}-\frac{1}{d_{i}^{\mathrm{Born}}+\Delta d_{i}}\right)\right)\$$

Version 1b: Equation from paper with modification to include the first
fraction for the DS term (probably wrong):

$$a_{\mathrm{SSM+DS}}^{\mathrm{Born}}=-\frac{e^{2}}{4\pi\varepsilon_{0}k_{\mathrm{B}}T}\left(1-\frac{1}{\varepsilon_{r,\mathrm{bulk}}}\right)\sum_{i}\frac{x_{i}z_{i}^{2}}{d_{i}^{\mathrm{Born}}+\Delta d_{i}} -\frac{e^{2}}{4\pi\varepsilon_{0}k_{\mathrm{B}}T}\left(1-\frac{1}{\varepsilon_{r,\mathrm{bulk}}}\right)\left(\sum_{i}x_{i}z_{i}^{2}\left(\frac{1}{d_{i}^{\mathrm{Born}}}-\frac{1}{d_{i}^{\mathrm{Born}}+\Delta d_{i}}\right)\right)\$$

Version 2 (my version): potential fix (changing the second
$\varepsilon_{r,\mathrm{solv}}$ to $\varepsilon_{r,\mathrm{ion}}$ and
including the first fraction for the DS term)

$$a_{\mathrm{SSM}}^{\mathrm{Born}}=-\frac{e^{2}}{4\pi\varepsilon_{0}k_{\mathrm{B}}T}\sum_{i}x_{i}z_{i}^{2} \left(1-\frac{1}{\varepsilon_{r,\mathrm{bulk}}}\right)\left(\frac{1}{d_{i}^{\mathrm{Born}}+\Delta d_{i}}\right)$$

$$a_{\mathrm{DS}}^{\mathrm{Born}}=-\frac{e^{2}}{4\pi\varepsilon_{0}k_{\mathrm{B}}T}\sum_{i}x_{i}z_{i}^{2} \left(1-\frac{1}{\varepsilon_{r,\mathrm{ion}}}\right)\left(\frac{1}{d_{i}^{\mathrm{Born}}}-\frac{1}{d_{i}^{\mathrm{Born}}+\Delta d_{i}}\right)$$

$$a_{\mathrm{SSM+DS}}^{\mathrm{Born}}=-\frac{e^{2}}{4\pi\varepsilon_{0}k_{\mathrm{B}}T}\sum_{i}x_{i}z_{i}^{2}\left(1-\frac{1}{\varepsilon_{r,\mathrm{bulk}}}\right)\left(\frac{1}{d_{i}^{\mathrm{Born}}+\Delta d_{i}}\right) -\frac{e^{2}}{4\pi\varepsilon_{0}k_{\mathrm{B}}T}\sum_{i}x_{i}z_{i}^{2} \left(1-\frac{1}{\varepsilon_{r,\mathrm{ion}}}\right)\left(\frac{1}{d_{i}^{\mathrm{Born}}}-\frac{1}{d_{i}^{\mathrm{Born}}+\Delta d_{i}}\right)$$

$$\Delta d_i=\frac{(f_{mix}-1)}{|z_i|}\cdot d_i^{\mathrm{Born}}$$

$$f_{mix}=\sum_{k}x_{k}f_{k}$$

$$d^{\text{Born}}_{i} = d^{\text{Born, fitted}}_{i} \qquad i\in\mathcal{I}$$

# Pressure (Compressibility Factor)

$$P=ZkT\rho\left(10^{10}\frac{\hat{\mathrm{A}}}{\mathrm{m}}\right)^{3}$$

$$Z=1+\eta\left(\frac{\partial\tilde{a}^\mathrm{res}}{\partial\eta}\right)_{T,x_i}$$
$$\eta\left(\frac{\partial\tilde{a}^\mathrm{res}}{\partial\eta}\right)_{T,x_i} = \rho\left(\frac{\partial\tilde{a}^\mathrm{res}}{\partial\rho}\right)_{T,x_i}$$

$$Z=1+\rho\left(\frac{\partial\tilde{a}^\mathrm{res}}{\partial\rho}\right)_{T,x_i}$$

$$Z = 1 + Z^{hc} + Z^{disp} + Z^{assoc} + Z^{DH} + Z^{Born}$$

$$Z^\alpha =\rho\left(\frac{\partial\tilde{a}^\alpha}{\partial\rho}\right)_{T,x}
.$$

## Hard-Chain Reference Contribution

$$Z^{hc} =\rho\left(\frac{\partial\tilde{a}^{hc}}{\partial\rho}\right)_{T,x}
.$$

$$Z^\mathrm{hc}=\bar{m}Z^\mathrm{hs}-\sum_ix_\mathrm{i}(m_i-1)(g_{ii}^\mathrm{hs})^{-1}\rho\frac{\partial g_{ii}^\mathrm{hs}}{\partial\rho}$$

$$Z^{hs} =\rho\left(\frac{\partial\tilde{a}^{hs}}{\partial\rho}\right)_{T,x}
.$$

$$Z^{\mathrm{hs}}=\frac{\zeta_{3}}{(1-\zeta_{3})}+\frac{3\zeta_{1}\zeta_{2}}{\zeta_{0}(1-\zeta_{3})^{2}}+\frac{3\zeta_{2}^{3}-\zeta_{3}\zeta_{2}^{3}}{\zeta_{0}(1-\zeta_{3})^{3}}$$

$$\begin{aligned}
\rho\frac{\partial g_{ij}^{\mathrm{hs}}}{\partial\rho}=\frac{\zeta_{3}}{\left(1-\zeta_{3}\right)^{2}}+\left(\frac{d_{i}d_{j}}{d_{i}+d_{j}}\right)\left(\frac{3\zeta_{2}}{\left(1-\zeta_{3}\right)^{2}}+\frac{6\zeta_{2}\zeta_{3}}{\left(1-\zeta_{3}\right)^{3}}\right)\\+\left(\frac{d_{i}d_{j}}{d_{i}+d_{j}}\right)^{2}\left(\frac{4{\zeta_{2}}^{2}}{\left(1-{\zeta_{3}}\right)^{3}}+\frac{6{\zeta_{2}}^{2}\zeta_{3}}{\left(1-{\zeta_{3}}\right)^{4}}\right)
\end{aligned}$$

## Dispersion Contribution

$$Z^{disp} =\rho\left(\frac{\partial\tilde{a}^{disp}}{\partial\rho}\right)_{T,x}
.$$

$$Z^{\mathrm{disp}}=-2\pi\rho\frac{\partial(\eta I_{1})}{\partial\eta}\overline{m^{2}\epsilon\sigma^{3}}-\\\pi\rho\bar{m}\left[C_{1}\frac{\partial(\eta I_{2})}{\partial\eta}+C_{2}\eta I_{2}\right]\overline{m^{2}\epsilon^{2}\sigma^{3}}$$

$$\frac{\partial(\eta I_1)}{\partial\eta}=\sum_{j=0}^6a_j(\bar{m})(j+1)\eta^j$$

$$\frac{\partial(\eta I_2)}{\partial\eta}=\sum_{i=0}^6b_j(\bar{m})(j+1)\eta^i$$

$$C_{2}=\frac{\partial C_{1}}{\partial\eta}=- C_{1}^{2}\biggl(\bar{m}\frac{-4\eta^{2}+20\eta+8}{\left(1-\eta\right)^{5}}+\\(1-\bar{m})\frac{2\eta^{3}+12\eta^{2}-48\eta+40}{\left[(1-\eta)(2-\eta)\right]^{3}}\biggr)$$

## Association Contribution

$$Z^{assoc} =\rho\left(\frac{\partial\tilde{a}^{assoc}}{\partial\rho}\right)_{T,x}
.$$

$$Z^{assoc}=\rho\sum_{i}x_{i}\sum_{A_{i}}\left(\frac{1}{X^{A_{i}}}-\frac{1}{2}\right)\left(\frac{\partial X^{A_{i}}}{\partial\rho}\right)_{T,x}.$$

$$\left(\frac{\partial X^{A_{i}}}{\partial\rho}\right)_{T,x}=-(X^{A_{i}})^{2}\sum_{j}\sum_{B_{j}}x_{j}\left[X^{B_{j}}\Delta^{A_{i}B_{j}}+\rho\left(\Delta^{A_{i}B_{j}}\left(\frac{\partial X^{B_{j}}}{\partial\rho}\right)_{T,x}+X^{B_{j}}\left(\frac{\partial\Delta^{A_{i}B_{j}}}{\partial\rho}\right)_{T,x}\right)\right]$$

$$\left(\frac{\partial\Delta^{A_iB_j}}{\partial\rho}\right)_{T,x}=d_{ij}^3\kappa^{A_iB_j}\left[\exp(\epsilon^{A_iB_j}/kT)-1\right]\left(\frac{\partial g_{ij}^{hs}}{\partial\rho}\right)_{T,x}$$

## Debye and Huckel Electrolyte Term Contribution

$$Z^{DH} =\rho\left(\frac{\partial\tilde{a}^{DH}}{\partial\rho}\right)_{T,x}
.$$

$$Z^{DH}=-\frac{\kappa e^2}{24\pi kT\epsilon}\sum_{i}x_{i}z_{i}{}^{2}\sigma_{i}$$

$$\sigma_i=\left(\frac{\partial(\kappa\chi_i)}{\partial\kappa}\right)_{T,\mathrm{x}}=-2\chi_i+\frac{3}{1+\kappa a_i}$$

### Bjerrum Treatment

$$Z^{DH}=-\frac{\kappa e^2}{24\pi kT\epsilon}\sum_{i}\alpha _ix_{i}z_{i}{}^{2}\sigma_{i}$$

$$\sigma_i=\left(\frac{\partial(\kappa\chi_i)}{\partial\kappa}\right)_{T,\mathrm{x}}=-2\chi_i+\frac{3}{1+\kappa R_i}$$

## Born Electrolyte Term Contribution

$$Z^{Born} =\rho\left(\frac{\partial\tilde{a}^{Born}}{\partial\rho}\right)_{T,x}
.$$

$$Z^{Born}= 0$$

# Chemical Potential and Fugacity Coefficient

$$\tilde{\mu_{k}}^{\mathrm{res}}=\frac{\mu_k^\mathrm{res}}{kT} = \frac{\hat{\mu}_k^\mathrm{res}}{RT}$$

$$\ln\varphi_k=\frac{\mu_k^\mathrm{res}}{kT}-\ln Z$$

$$\tilde{\mu_{k}}^{\mathrm{res}}=\tilde{a}^{\mathrm{res}}+(\mathrm{Z}-1)+\\\left(\frac{\partial\tilde{a}^{\mathrm{res}}}{\partial x_{k}}\right)_{T,\nu,x_{i\neq k}}-\sum_{j=1}^{N}\left[x_{j}\left(\frac{\partial\tilde{a}^{\mathrm{res}}}{\partial x_{j}}\right)_{T,\nu,x_{i\neq j}}\right]$$

$$\begin{aligned}
\tilde{\mu_{k}}^{\mathrm{res}}=\tilde{\mu_{k}}^{\mathrm{hc}}+\tilde{\mu_{k}}^{\mathrm{disp}}+\tilde{\mu_{k}}^{\mathrm{assoc}}+\tilde{\mu_{k}}^{\mathrm{DH}}+\tilde{\mu_{k}}^{\mathrm{Born}}
\end{aligned}$$

$$\zeta_{n,xk}=\left(\frac{\partial\zeta_{n}}{\partial x_{k}}\right)_{T,\boldsymbol{\rho},x_{j\neq k}}=\frac{\pi}{6}\rho m_{k}(d_{k})^{n}\quad n\in\{0,1,2,3\}$$

## Hard-Chain Reference Contribution

$$\tilde{\mu}^{hc}_k = \tilde{a}^{\mathrm{hc}} + Z^{hc} + \\\left(\frac{\partial\tilde{a}^{\mathrm{hc}}}{\partial x_{k}}\right)_{T,\nu,x_{i\neq k}}-\sum_{j=1}^{N}\left[x_{j}\left(\frac{\partial\tilde{a}^{\mathrm{hc}}}{\partial x_{j}}\right)_{T,\nu,x_{i\neq j}}\right]$$

$$\left(\frac{\partial\tilde{a}^{\mathrm{hc}}}{\partial x_{k}}\right)_{T,\boldsymbol{\rho},x_{j\neq k}}=m_{k}\tilde{a}^{\mathrm{hs}}+\bar{m}\left(\frac{\partial\tilde{a}^{\mathrm{hs}}}{\partial x_{k}}\right)_{T,\boldsymbol{\rho},x_{j\neq k}}-\\\sum_{i}x_{i}(m_{i}-1)(g_{ii}^{\mathrm{hs}})^{-1}\left(\frac{\partial g_{ii}^{\mathrm{hs}}}{\partial x_{k}}\right)_{T,\boldsymbol{\rho},x_{j\neq k}}$$

$$\begin{aligned}\left(\frac{\partial\tilde{a}^{\mathrm{hs}}}{\partial x_{k}}\right)_{T,\boldsymbol{\rho},x_{j\neq k}}=- \frac{\zeta_{0,xk}}{\zeta_{0}}\tilde{a}^{\mathrm{hs}}+\frac{1}{\zeta_{0}}\biggl[\frac{3(\zeta_{1,xk}\zeta_{2}+\zeta_{1}\zeta_{2,xk})}{(1-\zeta_{3})}+\\\frac{3\zeta_{1}\zeta_{2}\zeta_{3,xk}}{\left(1-\zeta_{3}\right)^{2}}+\frac{3\zeta_{2}^{2}\zeta_{2,xk}}{\zeta_{3}(1-\zeta_{3})^{2}}+\frac{\zeta_{2}^{3}\zeta_{3,xk}(3\zeta_{3}-1)}{\zeta_{3}^{2}(1-\zeta_{3})^{3}}+\\\left(\frac{3{\zeta_{2}}^{2}{\zeta_{2,xk}}{\zeta_{3}}-2{\zeta_{2}}^{3}{\zeta_{3,xk}}}{{\zeta_{3}}^{3}}-{\zeta_{0,xk}}\right)\mathrm{ln}(1-{\zeta_{3}})+\\\left(\zeta_{0}-\frac{\zeta_{2}^{3}}{\zeta_{3}^{2}}\right)\frac{\zeta_{3,xk}}{(1-\zeta_{3})}\biggr] \end{aligned}$$

$$\begin{aligned}\left(\frac{\partial g_{ij}^{\mathrm{by}}}{\partial x_{k}}\right)_{T,\boldsymbol{\rho},\boldsymbol{\rho}_{j\neq k}}&=\frac{\zeta_{3,x,k}}{\left(1-\zeta_{3}\right)^{2}}+\\&\left(\frac{d_{i}d_{j}}{d_{i}+d_{j}}\right)\left(\frac{3\zeta_{2,x,k}}{(1-\zeta_{3})^{2}}+\frac{6\zeta_{2}\zeta_{3,x,k}}{(1-\zeta_{3})^{3}}\right)+\\&\left(\frac{d_{i}d_{j}}{d_{i}+d_{j}}\right)^{2}\left(\frac{4\zeta_{2}\zeta_{2,x,k}}{(1-\zeta_{3})^{3}}+\frac{6\zeta_{2}^{2}\zeta_{3,x,k}}{(1-\zeta_{3})^{4}}\right)\end{aligned}$$

## Dispersion Contribution

$$\tilde{\mu}^{disp}_k = \tilde{a}^{\mathrm{disp}} + Z^{disp} + \\\left(\frac{\partial\tilde{a}^{\mathrm{disp}}}{\partial x_{k}}\right)_{T,\nu,x_{i\neq k}}-\sum_{j=1}^{N}\left[x_{j}\left(\frac{\partial\tilde{a}^{\mathrm{disp}}}{\partial x_{j}}\right)_{T,\nu,x_{i\neq j}}\right]$$

$$\overline{(m^{2}\epsilon\sigma^{3}})_{xk}=2m_{k}\sum_{j}x_{j}m_{j}\left(\frac{\epsilon_{kj}}{kT}\right)\sigma_{kj}^{3}$$

$$(\overline{m^2\epsilon^2\sigma^3})_{xk}=2m_k\sum_jx_jm_j\left(\frac{\epsilon_{kj}}{kT}\right)^2\sigma_{kj}^3$$

$$C_{1,xk}=C_{2}\zeta_{3,xk}-\\C_{1}^{2}\left\{m_{k}\frac{8\eta-2\eta^{2}}{\left(1-\eta\right)^{4}}-m_{k}\frac{20\eta-27\eta^{2}+12\eta^{3}-2\eta^{4}}{\left[(1-\eta)(2-\eta)\right]^{2}}\right\}$$

$$I_{1,xk}=\sum_{i=0}^{6}[a_{\mathrm{i}}(\bar{m})i\zeta_{3,xk}\eta^{i-1}+a_{\mathrm{i,x}k}\eta^{i}]$$

$$I_{2,xk}=\sum_{i=0}^{6}[b_{\mathrm{i}}(\bar{m})i\zeta_{3,xk}\eta^{i-1}+b_{\mathrm{i,xk}}\eta^{i}]$$

$$a_{\mathrm{i},xk}=\frac{m_{k}}{\bar{m}^{2}}a_{1i}+\frac{m_{k}}{\bar{m}^{2}}\left(3-\frac{4}{\bar{m}}\right)a_{2i}$$

$$b_{\mathrm{i},xk}=\frac{m_{k}}{\bar{m}^{2}}b_{1i}+\frac{m_{k}}{\bar{m}^{2}}\left(3-\frac{4}{\bar{m}}\right)b_{2i}$$

## Association Contribution

$$\tilde{\mu}^{assoc}_k = \tilde{a}^{\mathrm{assoc}} + Z^{assoc} + \\\left(\frac{\partial\tilde{a}^{\mathrm{assoc}}}{\partial x_{k}}\right)_{T,\nu,x_{i\neq k}}-\sum_{j=1}^{N}\left[x_{j}\left(\frac{\partial\tilde{a}^{\mathrm{assoc}}}{\partial x_{j}}\right)_{T,\nu,x_{i\neq j}}\right]$$

$$\left(\frac{\partial \tilde a^{assoc}}{\partial x_k}\right)_{T,v,x_{i\neq k}}
=
\sum_{A_k}\left(\ln X^{A_k}-\frac{X^{A_k}}{2}+\frac{1}{2}\right)
+
\sum_i x_i\sum_{A_i}
\left(\frac{1}{X^{A_i}}-\frac{1}{2}\right)
\left(\frac{\partial X^{A_i}}{\partial x_k}\right)_{T,v,x_{i\neq k}}$$

$$\left(\frac{\partial X^{A_i}}{\partial x_k}\right)_{T,v,x_{j\neq k}}
=
-(X^{A_i})^{2}\,\rho
\Bigg[
\sum_{B_k}X^{B_k}\Delta^{A_iB_k}
+
\sum_{j}x_j\sum_{B_j}
\left(
\Delta^{A_iB_j}\left(\frac{\partial X^{B_j}}{\partial x_k}\right)_{T,v,x_{j\neq k}}
+
X^{B_j}\left(\frac{\partial \Delta^{A_iB_j}}{\partial x_k}\right)_{T,v,x_{j\neq k}}
\right)
\Bigg]$$

$$\left(\frac{\partial \Delta^{A_iB_j}}{\partial x_k}\right)_{T,v,x_{j\neq k}}
=
d_{ij}^{3}\kappa^{A_iB_j}
\left[\exp\!\left(\frac{\epsilon^{A_iB_j}}{kT}\right)-1\right]
\left(\frac{\partial g_{ij}^{hs}}{\partial x_k}\right)_{T,\rho,x_{j\neq k}}$$

## Debye and Huckel Electrolyte Term Contribution

$$\tilde{\mu}^{DH}_k = \tilde{a}^{\mathrm{DH}} + Z^{DH} + \\\left(\frac{\partial\tilde{a}^{\mathrm{DH}}}{\partial x_{k}}\right)_{T,\nu,x_{i\neq k}}-\sum_{j=1}^{N}\left[x_{j}\left(\frac{\partial\tilde{a}^{\mathrm{DH}}}{\partial x_{j}}\right)_{T,\nu,x_{i\neq j}}\right]$$

The original 2005 version from

$$\begin{aligned}
\tilde{\mu}^{DH}_i= -\frac{q_i^2 \kappa}{24 \pi k T \epsilon}\left[2 \chi_i+\frac{\sum_j x_j q_j^2 \sigma_k}{\sum_j x_j q_j^2}\right]
\end{aligned}$$

$$\sigma_k=\left(\frac{\partial\left(\kappa \chi_k\right)}{\partial \kappa}\right)_{T, \mathrm{~N}}=-2 \chi_k+\frac{3}{1+\kappa d_k}$$

The updated version for concentration dependent dielectric constant from

$$\begin{aligned}
\left(\frac{\partial \tilde{a}^{DH}}{\partial x_{i}}\right)_{T,v,x_{j\neq i}}
=
-\frac{e^{2}}{12\pi\varepsilon_{0}k_{B}T}
\Bigg[
&\frac{1}{\varepsilon_{r}}\left(\frac{\partial\kappa}{\partial x_{i}}\right)_{T,v,x_{j\neq i}}
\sum_{j}x_{j}z_{j}^{2}\chi_{j}
\\
-&\frac{\kappa}{\varepsilon_{r}^{2}}
\left(\frac{\partial\varepsilon_{r}}{\partial x_{i}}\right)_{T,v,x_{j\neq i}}
\sum_{j}x_{j}z_{j}^{2}\chi_{j}
\\
+&\frac{\kappa}{\varepsilon_{r}}
\sum_{j}\chi_{j}z_{j}^{2}
\left(\frac{\partial x_{j}}{\partial x_{i}}\right)_{T,v,x_{k\neq i}}
\\
+&\frac{\kappa}{\varepsilon_{r}}
\sum_{j}x_{j}z_{j}^{2}
\left(\frac{\partial\chi_{j}}{\partial x_{i}}\right)_{T,v,x_{k\neq i}}
\Bigg]
\end{aligned}$$

$$\begin{aligned}
&\left(\frac{\partial\kappa}{\partial x_{i}}\right)_{T,v,x_{j\neq i}}=\frac12\left(\frac{\rho_{N}e^{2}}{k_{B}T\varepsilon_{0}\varepsilon_{r}}\sum_{j}x_{j}z_{j}^{2}\right)^{-\frac12}\\&\left\{-\frac{\rho_Ne^2}{k_BT\left(\varepsilon_0\varepsilon_r\right)^2}\varepsilon_0\frac{\partial\left(\varepsilon_r\right)}{\partial x_i}\sum_jx_jz_j^2+\frac{\rho_Ne^2}{k_BT\varepsilon_0\varepsilon_r}\alpha_iz_i^2\right\}\\&=\left(\frac{\rho_{N}e^{2}}{k_{B}T}\right)^{\frac{1}{2}}\left\{-\frac{1}{2}\left(\varepsilon_{0}\varepsilon_{r}\right)^{-\frac{3}{2}}\varepsilon_{0}\frac{\partial\left(\varepsilon_{r}\right)}{\partial x_{i}}\left[\sum_{j}x_{j}z_{j}^{2}\right]^{\frac{1}{2}}\right.\\&+\frac{1}{2\sqrt{\varepsilon_0\varepsilon_r}}\Bigg[\sum_jx_jz_j^2\Bigg]^{-\frac{1}{2}}\alpha_iz_i^2\Bigg\}
\end{aligned}$$

$$\begin{aligned}
\left(\frac{\partial\chi_{i}}{\partial x_{i}}\right)_{T,v,x_{j\neq i}}&=-\frac{9}{\left(\kappa d_{i}\right)^{4}}\biggl[\frac{3}{2}+ln(1+\kappa d_{i})-2(1+\kappa d_{i})\\&+\frac{1}{2}\left(1+\kappa d_{i}\right)^{2}\bigg]d_{i}\frac{\partial\kappa}{\partial x_{i}}+\frac{3}{\left(\kappa d_{i}\right)^{3}}\bigg[\frac{1}{1+\kappa d_{i}}-1+\kappa d_{i}\bigg]\\d_{i}\frac{\partial\kappa}{\partial x_{i}}&= 3\frac{\partial\kappa}{\partial x_{i}}\left\{-\frac{\chi_{i}}{\kappa}+\frac{d_{i}}{\left(\kappa d_{i}\right)^{3}}\biggl[\frac{1}{1+\kappa d_{i}}-1+\kappa d_{i}\biggr]\right\}
\end{aligned}$$

### Bjerrum Treatment

$$\begin{aligned}
\left(\frac{\partial \tilde{a}^{DH}}{\partial x_{i}}\right)_{T,v,x_{j\neq i}}
=
-\frac{e^{2}}{12\pi\varepsilon_{0}k_{B}T}
\Bigg[
&\frac{1}{\varepsilon_{r}}\left(\frac{\partial\kappa}{\partial x_{i}}\right)_{T,v,x_{j\neq i}}
\sum_{j}\alpha_{j}x_{j}z_{j}^{2}\chi_{j}
\\
-&\frac{\kappa}{\varepsilon_{r}^{2}}
\left(\frac{\partial\varepsilon_{r}}{\partial x_{i}}\right)_{T,v,x_{j\neq i}}
\sum_{j}\alpha_{j}x_{j}z_{j}^{2}\chi_{j}
\\
+&\frac{\kappa}{\varepsilon_{r}}\alpha_{i}z_{i}^{2}\chi_{i}
\\
+&\frac{\kappa}{\varepsilon_{r}}
\sum_{j}\alpha_{j}x_{j}z_{j}^{2}
\left(\frac{\partial\chi_{j}}{\partial x_{i}}\right)_{T,v,x_{k\neq i}}
\Bigg]
\end{aligned}$$

$$\begin{aligned}
&\left(\frac{\partial\kappa}{\partial x_{i}}\right)_{T,v,x_{j\neq i}}=\frac12\left(\frac{\rho_{N}e^{2}}{k_{B}T\varepsilon_{0}\varepsilon_{r}}\sum_{j}\alpha_{j}x_{j}z_{j}^{2}\right)^{-\frac12}\\&\left\{-\frac{\rho_Ne^2}{k_BT\left(\varepsilon_0\varepsilon_r\right)^2}\varepsilon_0\frac{\partial\left(\varepsilon_r\right)}{\partial x_i}\sum_j\alpha_{j}x_jz_j^2+\frac{\rho_Ne^2}{k_BT\varepsilon_0\varepsilon_r}\alpha_iz_i^2\right\}\\&=\left(\frac{\rho_{N}e^{2}}{k_{B}T}\right)^{\frac{1}{2}}\left\{-\frac{1}{2}\left(\varepsilon_{0}\varepsilon_{r}\right)^{-\frac{3}{2}}\varepsilon_{0}\frac{\partial\left(\varepsilon_{r}\right)}{\partial x_{i}}\left[\sum_{j}\alpha_{j}x_{j}z_{j}^{2}\right]^{\frac{1}{2}}\right.\\&+\frac{1}{2\sqrt{\varepsilon_0\varepsilon_r}}\Bigg[\sum_j\alpha_{j}x_jz_j^2\Bigg]^{-\frac{1}{2}}\alpha_iz_i^2\Bigg\}
\end{aligned}$$

$$\begin{aligned}
\left(\frac{\partial\chi_{i}}{\partial x_{i}}\right)_{T,v,x_{j\neq i}}&=-\frac{9}{\left(\kappa R_{i}\right)^{4}}\biggl[\frac{3}{2}+ln(1+\kappa R_{i})-2(1+\kappa R_{i})\\&+\frac{1}{2}\left(1+\kappa R_{i}\right)^{2}\bigg]R_{i}\frac{\partial\kappa}{\partial x_{i}}+\frac{3}{\left(\kappa R_{i}\right)^{3}}\bigg[\frac{1}{1+\kappa R_{i}}-1+\kappa R_{i}\bigg]\\R_{i}\frac{\partial\kappa}{\partial x_{i}}&= 3\frac{\partial\kappa}{\partial x_{i}}\left\{-\frac{\chi_{i}}{\kappa}+\frac{R_{i}}{\left(\kappa R_{i}\right)^{3}}\biggl[\frac{1}{1+\kappa R_{i}}-1+\kappa R_{i}\biggr]\right\}
\end{aligned}$$

### Dielectric Constant

Mole-Fraction
$$\left(\frac{\partial \varepsilon_r}{\partial x_i}\right)_{x_{j\neq i}}
= \varepsilon_{r, i}$$

Mass-Fraction
$$\left(\frac{\partial \varepsilon_r}{\partial x_i}\right)_{x_{j\neq i}}
=
\frac{MW_i}{\overline{MW}}\left(\varepsilon_{r,i}-\varepsilon_r\right)$$

Combo
$$\left(\frac{\partial \varepsilon_r}{\partial x_i}\right)_{x_{j\neq i}}
=
\varepsilon_{r,sol}^{(w)}
+
x_{sol}\,\frac{MW_i}{\overline{MW}_{sol}}
\left(\varepsilon_{r,i}-\varepsilon_{r,sol}^{(w)}\right),
\qquad i\in\mathcal S$$

New Mixing Rule
$$\varepsilon_{r,\mathrm{mix}}(\mathbf{x})
=
\frac{\varepsilon_{sf}(\mathbf{x})}{1+7.01\,x_{\mathrm{ion}}(\mathbf{x})},
\qquad
\varepsilon_{sf}\equiv \varepsilon_{r,\mathrm{solvent,mix}}^{\mathrm{salt\text{-}free}} .$$

$$\left(\frac{\partial \varepsilon_{r,\mathrm{mix}}}{\partial x_i}\right)_{T,v,x_{j\neq i}}
=
\frac{1}{1+7.01\,x_{\mathrm{ion}}}
\left(\frac{\partial \varepsilon_{sf}}{\partial x_i}\right)_{T,v,x_{j\neq i}}
-\frac{7.01\,\varepsilon_{sf}}{\left(1+7.01\,x_{\mathrm{ion}}\right)^2}
\left(\frac{\partial x_{\mathrm{ion}}}{\partial x_i}\right)_{T,v,x_{j\neq i}} .$$

$$\varepsilon_{sf}(\mathbf{x})
=
\sum_{s\in\mathcal{S}} w_{s}^{sf}\,\varepsilon_{r,s},
\qquad
w_{s}^{sf}
=
\frac{x_s M_s}{\displaystyle \sum_{m\in\mathcal{S}} x_m M_m},
\qquad
D\equiv \sum_{m\in\mathcal{S}} x_m M_m .$$

$$\left(\frac{\partial \varepsilon_{sf}}{\partial x_i}\right)_{T,v,x_{j\neq i}}
=
\sum_{s\in\mathcal{S}} \varepsilon_{r,s}
\left(\frac{\partial w_s^{sf}}{\partial x_i}\right)_{T,v,x_{j\neq i}}
=
\begin{cases}
\dfrac{M_i}{D}\left(\varepsilon_{r,i}-\varepsilon_{sf}\right), & i\in\mathcal{S},\\[10pt]
0, & i\notin\mathcal{S}.
\end{cases}$$

$$x_{\mathrm{ion}}(\mathbf{x})=\sum_{m\in\mathcal{I}} x_m,
\qquad
\left(\frac{\partial x_{\mathrm{ion}}}{\partial x_i}\right)_{T,v,x_{j\neq i}}
=
\begin{cases}
1, & i\in\mathcal{I},\\
0, & i\notin\mathcal{I}.
\end{cases}$$

$$\left(\frac{\partial \varepsilon_{r,\mathrm{mix}}}{\partial x_i}\right)_{T,v,x_{j\neq i}}
=
\begin{cases}
\dfrac{1}{1+7.01\,x_{\mathrm{ion}}}\;
\dfrac{M_i}{\displaystyle \sum_{m\in\mathcal{S}} x_m M_m}
\left(\varepsilon_{r,i}-\varepsilon_{sf}\right),
& i\in\mathcal{S},\\[14pt]
-\dfrac{7.01\,\varepsilon_{sf}}{\left(1+7.01\,x_{\mathrm{ion}}\right)^{2}},
& i\in\mathcal{I}.
\end{cases}$$

## Born Electrolyte Term Contribution

$$\tilde{\mu}^{Born}_k = \tilde{a}^{\mathrm{Born}} + Z^{Born} + \\\left(\frac{\partial\tilde{a}^{\mathrm{Born}}}{\partial x_{k}}\right)_{T,\nu,x_{i\neq k}}-\sum_{j=1}^{N}\left[x_{j}\left(\frac{\partial\tilde{a}^{\mathrm{Born}}}{\partial x_{j}}\right)_{T,\nu,x_{i\neq j}}\right]$$

Version 1: How the reference "ePC-SAFT advanced - Part I: Physical
meaning of including a concentration-dependent dielectric constant in
the born term and in the Debye-Hückel theory" has it (equations 3)
$$\begin{aligned}\left(\frac{\partial \tilde{a}^{Born}}{\partial x_{i}}\right)_{T,v_{N},x_{j\neq i}}=-\frac{e^{2}}{4\pi k_{B}T\varepsilon_{0}}\biggl[\left(1-\frac{1}{\varepsilon_{r}}\right)\frac{z_{i}^{2}}{d^{\text{Born}}_{i}} + \left(\frac{1}{\varepsilon_{r}^{2}}\right) \left(\frac{\partial\varepsilon_{r}}{\partial x_{i}}\right)_{T,v,x_{j\neq i}} \biggr]\end{aligned}$$

Version 2: Potentially Fixed Differential from the paper
$$\begin{aligned}\left(\frac{\partial \tilde{a}^{Born}}{\partial x_{i}}\right)_{T,v_{N},x_{j\neq i}}=-\frac{e^{2}}{4\pi k_{B}T\varepsilon_{0}}\biggl[\left(1-\frac{1}{\varepsilon_{r}}\right)\frac{z_{i}^{2}}{d^{\text{Born}}_{i}} + \left(\frac{1}{\varepsilon_{r}^{2}}\right)\sum_j \frac{x_j z_j^2}{d^{\text{Born}}_{j}} \left(\frac{\partial\varepsilon_{r}}{\partial x_{i}}\right)_{T,v,x_{j\neq i}} \biggr]\end{aligned}$$

### Solvation Shell Model + Dielectric Saturation

Version 3

$$\begin{aligned}
\left(\frac{\partial \tilde a_{\mathrm{SSM+DS}}^{\mathrm{Born}}}{\partial x_{k}}\right)_{T,v,x_{j\neq k}}
=
-\frac{e^{2}}{4\pi\varepsilon_{0}k_{\mathrm{B}}T}
\Bigg[
&\left(1-\frac{1}{\varepsilon_{r}}\right)\frac{z_{k}^{2}}{d_{k}^{\mathrm{Born}}+\Delta d_{k}}
\\
&+\left(1-\frac{1}{\varepsilon_{r,\mathrm{ion}}}\right) z_{k}^{2}
\left(\frac{1}{d_{k}^{\mathrm{Born}}}-\frac{1}{d_{k}^{\mathrm{Born}}+\Delta d_{k}}\right)
\\
&+\left(\frac{1}{\varepsilon_{r}^{2}}\right)
\left(\frac{\partial \varepsilon_{r}}{\partial x_{k}}\right)_{T,v,x_{j\neq k}}
\left(\sum_{i} x_{i} z_{i}^{2}\frac{1}{d_{i}^{\mathrm{Born}}+\Delta d_{i}}\right)
\Bigg].
\end{aligned}$$

Version 4 $$\begin{aligned}
\left(\frac{\partial a^{\mathrm{Born}}}{\partial x_i}\right)_{T,v,x_{j\neq i}}
&=
-\frac{e^2}{4\pi\varepsilon_0 k_B T}
\Bigg[
z_i^2\!\left(
\left(1-\frac{1}{\varepsilon_{r,\mathrm{ion}}}\right)\frac{1}{D_i}
+
\left(1-\frac{1}{\varepsilon_{r,\mathrm{solv}}}\right)\left(\frac{1}{d_i^{\mathrm{Born}}}-\frac{1}{D_i}\right)
\right)
\\
&\qquad\qquad
+\sum_j x_j z_j^2
\Bigg(
\left(\frac{1}{\varepsilon_{r,\mathrm{ion}}}-\frac{1}{\varepsilon_{r,\mathrm{solv}}}\right)
\frac{1}{D_j^2}
\frac{\partial \Delta d_m}{\partial x_i}
+
\frac{1}{\varepsilon_{r,\mathrm{solv}}^2}
\frac{\partial \varepsilon_{r,\mathrm{solv}}}{\partial x_i}
\left(\frac{1}{d_j^{\mathrm{Born}}}-\frac{1}{D_j}\right)
\Bigg)
\Bigg],
\end{aligned}$$

$$\frac{\partial \Delta d_m}{\partial x_i}=\frac{d_m^{\text {Born }}}{\left|z_m\right|} \frac{\partial f_{\text {mix }}}{\partial x_i}=\frac{d_m^{\text {Born }}}{\left|z_m\right|} f_i,$$

$$D_m=d_m^{\text {Born }}+\Delta d_m, \Delta d_m=\frac{\left(f_{\text {mix }}-1\right)}{\left|z_m\right|} d_m^{\text {Born }} .$$

# Enthalpy, Entropy, and Gibbs Free Energy

## Entropy

$$\tilde{h}^{\mathrm{res}} = \frac{\hat{h}^{\mathrm{res}}}{RT}=-T\left(\frac{\partial\tilde{a}^{\mathrm{res}}}{\partial T}\right)_{\rho,x_i}+(Z-1)$$

## Entropy

$$\tilde{s}^{\mathrm{res}} = \frac{\hat{s}^{\mathrm{res}}(P,T)}{R}=\frac{\hat{s}^{\mathrm{res}}(\nu,T)}{R}+\ln(Z)$$

$$\tilde{s}^{\mathrm{res}} = \frac{\hat{s}^{\mathrm{res}}(P,T)}{R}=-T\left[\left(\frac{\partial\tilde{a}^{\mathrm{res}}}{\partial T}\right)_{\rho,x_{i}}+\frac{\tilde{a}^{\mathrm{res}}}{T}\right]+\ln(Z)$$

## Gibbs Free Energy

$$\begin{aligned}
\tilde{g}^{\mathrm{res}}=\frac{\hat{g}^{\mathrm{res}}}{RT}=\frac{\hat{h}^{\mathrm{res}}}{RT}-\frac{\hat{s}^{\mathrm{res}}(P,T)}{R}
\end{aligned}$$

$$\tilde{g}^{\mathrm{res}}=\tilde{a}^{\mathrm{res}}+(Z-1)-\ln(Z)$$

Gibbs Energy of Solvation

$$\Delta^{\mathrm{solv},x}G_{i}^{\infty}(T,p,x_{j\neq i},x_{i}\to0)=RT\mathrm{ln}(\varphi_{i}^{\infty}(T,p,x_{j\neq i},x_{i}\to0))\quad(17)$$

Gibbs Energy of Transfer at Infinite Dilution from solvent S1 to solvent
S2
$$\Delta^{\mathrm{tr}}G_{i}^{\infty,\mathrm{S}1\to\mathrm{S}2}(T,p,x_{j\neq i},x_{i}\rightarrow0)=RT\mathrm{ln}\left(\frac{\varphi_{i}^{\infty,\mathrm{S}2}}{\varphi_{i}^{\infty,\mathrm{S}1}}\right)$$

## Temperature Differential Equations

$$\left(\frac{\partial\tilde{a}^\mathrm{res}}{\partial T}\right)_{\rho,x_i} =\left(\frac{\partial\tilde{a}^\mathrm{hc}}{\partial T}\right)_{\rho,x_i} +\left(\frac{\partial\tilde{a}^\mathrm{disp}}{\partial T}\right)_{\rho,x_i}
+\left(\frac{\partial\tilde{a}^\mathrm{assoc}}{\partial T}\right)_{\rho,x_i}
+\left(\frac{\partial\tilde{a}^\mathrm{DH}}{\partial T}\right)_{\rho,x_i}
+\left(\frac{\partial\tilde{a}^\mathrm{Born}}{\partial T}\right)_{\rho,x_i}$$

### Hard-Chain Reference Contribution

$$\left(\frac{\partial\tilde{a}^{\mathrm{hc}}}{\partial T}\right)_{\rho,x_{i}}=\bar{m}\left(\frac{\partial\tilde{a}^{\mathrm{hs}}}{\partial T}\right)_{\rho,x_{i}}-\sum_{i}x_{\mathrm{i}}(m_{\mathrm{i}}-1)(g_{ii}^{\mathrm{hs}})^{-1}\left(\frac{\partial g_{ii}^{\mathrm{hs}}}{\partial T}\right)_{\rho,x_{i}}(\mathrm{A.54})$$

$$\begin{gathered}
\left(\frac{\partial\tilde{a}^{hs}}{\partial T}\right)_{\rho,x_{i}}=\frac{1}{\zeta_{0}}\left[\frac{3(\zeta_{1,T}\zeta_{2}+\zeta_{1}\zeta_{2,T})}{(1-\zeta_{3})}+\frac{3\zeta_{1}\zeta_{2}\zeta_{3,T}}{(1-\zeta_{3})^{2}}+\right.\\\frac{3{\zeta_{2}}^{2}{\zeta_{2,T}}}{\zeta_{3}{(1-\zeta_{3})}^{2}}+\frac{{\zeta_{2}}^{3}{\zeta_{3,T}}(3{\zeta_{3}}-1)}{{\zeta_{3}}^{2}{(1-\zeta_{3})}^{3}}+\\\left(\frac{3{\zeta_{2}}^{2}{\zeta_{2,T}}{\zeta_{3}}-2{\zeta_{2}}^{3}{\zeta_{3,T}}}{{\zeta_{3}}^{3}}\right)\ln(1-{\zeta_{3}})+\left({\zeta_{0}}-\frac{{\zeta_{2}}^{3}}{{\zeta_{3}}^{2}}\right)\frac{\zeta_{3,T}}{(1-\zeta_{3})}
\end{gathered}$$

$$\begin{gathered}
\frac{\partial g_{ii}^{\mathrm{hs}}}{\partial T}=\frac{\zeta_{3,T}}{\left(1-\zeta_{3}\right)^{2}}+\left(\frac{1}{2}d_{i,T}\right)\frac{3\zeta_{2}}{\left(1-\zeta_{3}\right)^{2}}+\\\left(\frac{1}{2}d_{i}\right)\left(\frac{3\zeta_{2,T}}{\left(1-\zeta_{3}\right)^{2}}+\frac{6\zeta_{2}\zeta_{3,T}}{\left(1-\zeta_{3}\right)^{3}}\right)+\left(\frac{1}{2}d_{i}d_{i,T}\right)\frac{2\zeta_{2}^{2}}{\left(1-\zeta_{3}\right)^{3}}+\\\left(\frac{1}{2}d_{i}\right)^{2}\left(\frac{4\xi_{2}\xi_{2,T}}{\left(1-\xi_{3}\right)^{3}}+\frac{6\xi_{2}{}^{2}\xi_{3,T}}{\left(1-\xi_{3}\right)^{4}}\right)
\end{gathered}$$

$$d_{i,T}=\frac{\partial d_i}{\partial T}=\sigma_i\left(3\frac{\epsilon_i}{kT^2}\right)\left[-0.12\exp\left(-3\frac{\epsilon_i}{kT}\right)\right]$$

$$\zeta_{n,T}=\frac{\partial\zeta_{n}}{\partial T}=\frac{\pi}{6}\rho\sum_{i}x_{i}m_{i}nd_{i,T}\left(d_{i}\right)^{n-1}\quad n\in\{1,2,3\}$$

$$\frac{1}{2}d_i=\left(\frac{d_id_i}{d_i+d_i}\right)$$

### Dispersion Contribution

$$\begin{aligned}
\left(\frac{\partial\tilde{a}^{\mathrm{disp}}}{\partial T}\right)_{\rho,x_{i}}&=-2\pi\rho\left(\frac{\partial I_{1}}{\partial T}-\frac{I_{1}}{T}\right)\overline{m^{2}\epsilon\sigma^{3}}-\\&\pi\rho\bar{m}\left[\frac{\partial C_{1}}{\partial T}I_{2}+C_{1}\frac{\partial I_{2}}{\partial T}-2C_{1}\frac{I_{2}}{T}\right]\overline{m^{2}\epsilon^{2}\sigma^{3}}
\end{aligned}$$

### Debye and Huckel Electrolyte Term Contribution

$$\frac{d \tilde{a}^{D H}}{d T}=-\frac{1}{12 \pi k_B \varepsilon_0 \varepsilon_r} \sum_i x_i{\left(z_i e\right)}^{2}\left[\left(-2 \chi_i+\frac{3}{1+\kappa a_i}\right)\left(-\frac{\kappa}{2 T^2}\right)-\frac{\kappa \chi_i}{T^2}\right]$$

### Born Contribution

$$\begin{aligned}
\left(\frac{\partial\tilde{a}^{\mathrm{Born}}}{\partial T}\right)_{\rho,x_{i}}=\frac{e^2}{4 \pi \varepsilon_0 k_B T^2}\left(1-\frac{1}{\varepsilon_r}\right) \sum_i \frac{x_i z_i^2}{a_i}
\end{aligned}$$

# Salt Basis Conversions

$$\gamma_{\pm}^{*,m}=\frac{\gamma_{\pm}^{*,x}}{1+M_{\mathrm{solvent}}\cdot\tilde{m}_{\mathrm{solute}}\cdot\sum_{i}\nu_{i,\mathrm{ion}}}$$
