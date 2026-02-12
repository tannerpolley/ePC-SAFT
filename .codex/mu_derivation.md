# Debye–Hückel (DH) chemical potential derivation (no skipped steps)

This document derives, step-by-step, the closed-form expression for the DH contribution to the *dimensionless* chemical potential
\[
\tilde{\mu}^{DH}_k \equiv \frac{\mu^{DH}_k}{k_B T},
\]
starting from the definitions of \(\tilde a^{DH}\), \(\kappa\), \(\chi_i\), \(Z^{DH}\), and the mixture chemical-potential identity.

---

## 0) Assumptions required for the closed-form result

The final compact expression
\[
\tilde{\mu}^{DH}_k= -\frac{q_k^2 \kappa}{24 \pi k_B T \,\epsilon}\left[2 \chi_k+\frac{\sum_i x_i q_i^2 \sigma_i}{\sum_i x_i q_i^2}\right]
\]
is obtained **only if**:

1. **\(\epsilon=\varepsilon_0\varepsilon_r\) is constant with respect to composition**, i.e.
   \[
   \left(\frac{\partial \epsilon}{\partial x_k}\right)_{T,\nu,x_{i\neq k}} = 0.
   \]
   If \(\varepsilon_r=\varepsilon_r(\mathbf x)\), extra terms appear and the closed form changes.

2. The ion size parameters \(a_i\) (appearing inside \(\chi_i\)) are treated as constants with respect to composition.

---

## 1) Given starting equations

### 1.1 DH reduced Helmholtz energy
\[
\tilde{a}^{DH}=-\frac{\kappa e^{2}}{12\pi\varepsilon_{0}\varepsilon_{r}k_{B}T}\sum_{i}x_{i}z_{i}^{2}\chi_{i}
\tag{1}
\]

Let \(\epsilon\equiv \varepsilon_0\varepsilon_r\). Then
\[
\tilde{a}^{DH}=-\frac{\kappa e^{2}}{12\pi\epsilon k_{B}T}\sum_{i}x_{i}z_{i}^{2}\chi_{i}.
\tag{1'}
\]

### 1.2 Screening parameter
\[
\kappa=\sqrt{\frac{\rho e^{2}}{k_{B}T\varepsilon_{0}\varepsilon_{r}}\sum_{j}x_{j}z_{j}^{2}}
=
\sqrt{\frac{\rho e^{2}}{k_{B}T\epsilon}\sum_{j}x_{j}z_{j}^{2}}.
\tag{2}
\]

### 1.3 \(\chi_i\) function
\[
\chi_i=\frac{3}{\left(\kappa a_i\right)^3}\left[\frac{3}{2}+\ln(1+\kappa a_i)-2(1+\kappa a_i)+\frac{1}{2}\left(1+\kappa a_i\right)^2\right].
\tag{3}
\]

### 1.4 DH compressibility contribution
\[
Z^{DH}=-\frac{\kappa e^2}{24\pi k_B T \epsilon}\sum_{i}x_{i}z_{i}^{2}\sigma_{i}.
\tag{4}
\]

### 1.5 \(\sigma_i\) definition
\[
\sigma_i=\left(\frac{\partial(\kappa\chi_i)}{\partial\kappa}\right)_{T,\mathbf{x}}
=-2\chi_i+\frac{3}{1+\kappa a_i}.
\tag{5}
\]

### 1.6 Mixture chemical-potential identity
\[
\tilde{\mu}^{DH}_k = \tilde{a}^{DH} + Z^{DH}
+ \left(\frac{\partial\tilde{a}^{DH}}{\partial x_{k}}\right)_{T,\nu,x_{i\neq k}}
-\sum_{j=1}^{N}\left[x_{j}\left(\frac{\partial\tilde{a}^{DH}}{\partial x_{j}}\right)_{T,\nu,x_{i\neq j}}\right].
\tag{6}
\]

---

## 2) Notation simplifications (purely algebraic)

Define the constant
\[
C \equiv \frac{e^2}{12\pi k_B T\,\epsilon}.
\tag{7}
\]

Define the three composition-dependent sums:
\[
S \equiv \sum_i x_i z_i^2 \chi_i,
\qquad
T \equiv \sum_i x_i z_i^2 \sigma_i,
\qquad
Q \equiv \sum_i x_i z_i^2.
\tag{8}
\]

Then Eq. (1′) becomes:
\[
\tilde a^{DH} = -C\,\kappa\,S.
\tag{9}
\]

Eq. (4) becomes:
\[
Z^{DH} = -\frac{C}{2}\,\kappa\,T
\quad\text{because}\quad
\frac{e^2}{24\pi k_B T\epsilon} = \frac{1}{2}\frac{e^2}{12\pi k_B T\epsilon}=\frac{C}{2}.
\tag{10}
\]

---

## 3) Derivative of \(\kappa\) with respect to \(x_k\) (no skipped steps)

From Eq. (2):
\[
\kappa = \sqrt{\frac{\rho e^2}{k_B T\epsilon} \, Q}.
\tag{11}
\]

Define
\[
B \equiv \frac{\rho e^2}{k_B T\epsilon},
\quad\Rightarrow\quad
\kappa = \sqrt{B Q}.
\tag{12}
\]

Square both sides:
\[
\kappa^2 = BQ.
\tag{13}
\]

Differentiate w.r.t. \(x_k\) at constant \(T,\nu,x_{i\neq k}\):
\[
\frac{\partial (\kappa^2)}{\partial x_k} = \frac{\partial (BQ)}{\partial x_k}.
\tag{14}
\]

Left-hand side:
\[
\frac{\partial (\kappa^2)}{\partial x_k} = 2\kappa \frac{\partial \kappa}{\partial x_k}.
\tag{15}
\]

Right-hand side: \(B\) is constant (since \(\rho,\epsilon,T\) are constant here), so:
\[
\frac{\partial (BQ)}{\partial x_k} = B\frac{\partial Q}{\partial x_k}.
\tag{16}
\]

Compute \(\partial Q/\partial x_k\) from \(Q=\sum_i x_i z_i^2\):
\[
\frac{\partial Q}{\partial x_k}
=
\frac{\partial}{\partial x_k}\left(\sum_i x_i z_i^2\right)
=
\sum_i z_i^2 \frac{\partial x_i}{\partial x_k}.
\tag{17}
\]

Since \(x_{i\neq k}\) are held fixed:
\[
\frac{\partial x_i}{\partial x_k} =
\begin{cases}
1,& i=k\\
0,& i\neq k
\end{cases}
= \delta_{ik}.
\tag{18}
\]

Therefore:
\[
\frac{\partial Q}{\partial x_k} = \sum_i z_i^2 \delta_{ik} = z_k^2.
\tag{19}
\]

Insert (15), (16), (19) into (14):
\[
2\kappa \frac{\partial \kappa}{\partial x_k} = B z_k^2.
\tag{20}
\]

Solve:
\[
\frac{\partial \kappa}{\partial x_k} = \frac{B z_k^2}{2\kappa}.
\tag{21}
\]

Now use (13) to eliminate \(B\): from \(\kappa^2=BQ\Rightarrow B=\kappa^2/Q\). Substitute:
\[
\frac{\partial \kappa}{\partial x_k} = \frac{(\kappa^2/Q) z_k^2}{2\kappa}
= \frac{\kappa z_k^2}{2Q}.
\tag{22}
\]

**Result:**
\[
\boxed{
\frac{\partial \kappa}{\partial x_k} = \frac{\kappa z_k^2}{2\sum_i x_i z_i^2}
}
\tag{23}
\]

---

## 4) Derivative of \(S=\sum_i x_i z_i^2\chi_i\) with respect to \(x_k\)

Start:
\[
S=\sum_i x_i z_i^2 \chi_i.
\tag{24}
\]

Differentiate:
\[
\frac{\partial S}{\partial x_k}
=
\frac{\partial}{\partial x_k}\left(\sum_i x_i z_i^2 \chi_i\right)
=
\sum_i z_i^2 \frac{\partial}{\partial x_k}(x_i \chi_i).
\tag{25}
\]

Apply product rule to \(x_i\chi_i\):
\[
\frac{\partial}{\partial x_k}(x_i\chi_i)
=
\chi_i\frac{\partial x_i}{\partial x_k} + x_i\frac{\partial \chi_i}{\partial x_k}.
\tag{26}
\]

Therefore:
\[
\frac{\partial S}{\partial x_k}
=
\sum_i z_i^2\left(\chi_i\frac{\partial x_i}{\partial x_k}\right)
+
\sum_i z_i^2\left(x_i\frac{\partial \chi_i}{\partial x_k}\right).
\tag{27}
\]

Use \(\partial x_i/\partial x_k=\delta_{ik}\):
\[
\sum_i z_i^2 \chi_i \delta_{ik} = z_k^2 \chi_k.
\tag{28}
\]

So:
\[
\frac{\partial S}{\partial x_k}
=
z_k^2\chi_k
+
\sum_i x_i z_i^2\frac{\partial \chi_i}{\partial x_k}.
\tag{29}
\]

Now use the chain rule for \(\chi_i(\kappa)\):
\[
\frac{\partial \chi_i}{\partial x_k}
=
\frac{d\chi_i}{d\kappa}\frac{\partial \kappa}{\partial x_k}.
\tag{30}
\]

Insert into (29):
\[
\frac{\partial S}{\partial x_k}
=
z_k^2\chi_k
+
\left(\frac{\partial\kappa}{\partial x_k}\right)\sum_i x_i z_i^2\frac{d\chi_i}{d\kappa}.
\tag{31}
\]

---

## 5) Express \(\frac{d\chi_i}{d\kappa}\) in terms of \(\sigma_i\) (no skipped steps)

Definition:
\[
\sigma_i = \frac{\partial (\kappa\chi_i)}{\partial \kappa}.
\tag{32}
\]

Differentiate \(\kappa\chi_i\) using product rule:
\[
\frac{\partial(\kappa\chi_i)}{\partial\kappa}
=
\chi_i + \kappa\frac{d\chi_i}{d\kappa}.
\tag{33}
\]

Thus:
\[
\sigma_i = \chi_i + \kappa\frac{d\chi_i}{d\kappa}.
\tag{34}
\]

Solve for \(d\chi_i/d\kappa\):
\[
\kappa\frac{d\chi_i}{d\kappa} = \sigma_i - \chi_i
\quad\Rightarrow\quad
\frac{d\chi_i}{d\kappa} = \frac{\sigma_i - \chi_i}{\kappa}.
\tag{35}
\]

Now form the sum in (31):
\[
\sum_i x_i z_i^2\frac{d\chi_i}{d\kappa}
=
\sum_i x_i z_i^2 \frac{\sigma_i-\chi_i}{\kappa}
=
\frac{1}{\kappa}\left(\sum_i x_i z_i^2\sigma_i - \sum_i x_i z_i^2\chi_i\right).
\tag{36}
\]

Recognize the sums:
\[
\sum_i x_i z_i^2\sigma_i = T,
\qquad
\sum_i x_i z_i^2\chi_i = S.
\tag{37}
\]

So:
\[
\sum_i x_i z_i^2\frac{d\chi_i}{d\kappa} = \frac{T-S}{\kappa}.
\tag{38}
\]

Insert into (31):
\[
\frac{\partial S}{\partial x_k}
=
z_k^2\chi_k
+
\left(\frac{\partial\kappa}{\partial x_k}\right)\frac{T-S}{\kappa}.
\tag{39}
\]

Now insert (23): \(\partial\kappa/\partial x_k = \kappa z_k^2/(2Q)\).
Then:
\[
\left(\frac{\partial\kappa}{\partial x_k}\right)\frac{1}{\kappa}
=
\left(\frac{\kappa z_k^2}{2Q}\right)\frac{1}{\kappa}
=
\frac{z_k^2}{2Q}.
\tag{40}
\]

So (39) becomes:
\[
\boxed{
\frac{\partial S}{\partial x_k}
=
z_k^2\chi_k
+
\frac{z_k^2}{2Q}(T-S)
}
\tag{41}
\]

---

## 6) Derivative \(\partial\tilde a^{DH}/\partial x_k\) (no skipped steps)

Recall:
\[
\tilde a^{DH} = -C\kappa S.
\tag{42}
\]

Differentiate w.r.t. \(x_k\):
\[
\frac{\partial \tilde a^{DH}}{\partial x_k}
=
-C \frac{\partial}{\partial x_k}(\kappa S).
\tag{43}
\]

Product rule for \(\kappa S\):
\[
\frac{\partial}{\partial x_k}(\kappa S)
=
\left(\frac{\partial\kappa}{\partial x_k}\right)S
+
\kappa\left(\frac{\partial S}{\partial x_k}\right).
\tag{44}
\]

Thus:
\[
\frac{\partial \tilde a^{DH}}{\partial x_k}
=
-C\left[
\left(\frac{\partial\kappa}{\partial x_k}\right)S
+
\kappa\left(\frac{\partial S}{\partial x_k}\right)
\right].
\tag{45}
\]

Insert (23) and (41).

### 6.1 First term
\[
\left(\frac{\partial\kappa}{\partial x_k}\right)S
=
\left(\frac{\kappa z_k^2}{2Q}\right)S
=
\frac{\kappa z_k^2 S}{2Q}.
\tag{46}
\]

### 6.2 Second term
\[
\kappa\left(\frac{\partial S}{\partial x_k}\right)
=
\kappa\left[z_k^2\chi_k + \frac{z_k^2}{2Q}(T-S)\right]
=
\kappa z_k^2\chi_k + \frac{\kappa z_k^2}{2Q}(T-S).
\tag{47}
\]

### 6.3 Combine inside the bracket
Sum (46) + (47):
\[
\frac{\kappa z_k^2 S}{2Q}
+
\kappa z_k^2\chi_k
+
\frac{\kappa z_k^2}{2Q}(T-S).
\tag{48}
\]

Combine the two \((\cdot)/(2Q)\) pieces:
\[
\frac{\kappa z_k^2 S}{2Q}
+
\frac{\kappa z_k^2}{2Q}(T-S)
=
\frac{\kappa z_k^2}{2Q}\left[S + (T-S)\right]
=
\frac{\kappa z_k^2}{2Q}T.
\tag{49}
\]

So (48) becomes:
\[
\kappa z_k^2\chi_k + \frac{\kappa z_k^2}{2Q}T.
\tag{50}
\]

Therefore:
\[
\boxed{
\frac{\partial \tilde a^{DH}}{\partial x_k}
=
-C\left[\kappa z_k^2\chi_k + \frac{\kappa z_k^2}{2Q}T\right]
}
\tag{51}
\]

Factor common terms:
\[
\boxed{
\frac{\partial \tilde a^{DH}}{\partial x_k}
=
-\,C\,\kappa z_k^2\left[\chi_k + \frac{T}{2Q}\right]
}
\tag{52}
\]

---

## 7) Compute \(\sum_j x_j(\partial\tilde a^{DH}/\partial x_j)\) (no skipped steps)

Start from (52) with \(k\to j\):
\[
\frac{\partial \tilde a^{DH}}{\partial x_j}
=
-\,C\,\kappa z_j^2\left[\chi_j + \frac{T}{2Q}\right].
\tag{53}
\]

Multiply by \(x_j\) and sum over \(j\):
\[
\sum_j x_j\frac{\partial \tilde a^{DH}}{\partial x_j}
=
\sum_j x_j\left(-C\kappa z_j^2\left[\chi_j + \frac{T}{2Q}\right]\right)
=
-\,C\,\kappa\sum_j x_j z_j^2\left[\chi_j + \frac{T}{2Q}\right].
\tag{54}
\]

Split the sum:
\[
\sum_j x_j z_j^2\left[\chi_j + \frac{T}{2Q}\right]
=
\sum_j x_j z_j^2\chi_j + \sum_j x_j z_j^2\left(\frac{T}{2Q}\right).
\tag{55}
\]

The first sum is \(S\):
\[
\sum_j x_j z_j^2\chi_j = S.
\tag{56}
\]

For the second sum, \(\frac{T}{2Q}\) is constant w.r.t. \(j\), so:
\[
\sum_j x_j z_j^2\left(\frac{T}{2Q}\right)
=
\frac{T}{2Q}\sum_j x_j z_j^2
=
\frac{T}{2Q}Q
=
\frac{T}{2}.
\tag{57}
\]

Thus (55) becomes:
\[
\sum_j x_j z_j^2\left[\chi_j + \frac{T}{2Q}\right]
=
S + \frac{T}{2}.
\tag{58}
\]

Insert into (54):
\[
\boxed{
\sum_j x_j\frac{\partial \tilde a^{DH}}{\partial x_j}
=
-\,C\,\kappa\left(S + \frac{T}{2}\right)
}
\tag{59}
\]

---

## 8) Assemble \(\tilde\mu_k^{DH}\) and simplify fully (no skipped steps)

Recall the identity (6):
\[
\tilde{\mu}^{DH}_k = \tilde{a}^{DH} + Z^{DH}
+ \left(\frac{\partial\tilde{a}^{DH}}{\partial x_{k}}\right)
-\sum_{j}\left[x_{j}\left(\frac{\partial\tilde{a}^{DH}}{\partial x_{j}}\right)\right].
\tag{60}
\]

Substitute:
- \(\tilde a^{DH} = -C\kappa S\) from (9)
- \(Z^{DH} = -\frac{C}{2}\kappa T\) from (10)
- \(\partial\tilde a^{DH}/\partial x_k = -C\kappa z_k^2\left[\chi_k + \frac{T}{2Q}\right]\) from (52)
- \(\sum_j x_j(\partial\tilde a^{DH}/\partial x_j) = -C\kappa\left(S + \frac{T}{2}\right)\) from (59)

So:
\[
\tilde\mu_k^{DH}
=
(-C\kappa S)
+
\left(-\frac{C}{2}\kappa T\right)
+
\left(-C\kappa z_k^2\left[\chi_k + \frac{T}{2Q}\right]\right)
-
\left(-C\kappa\left[S+\frac{T}{2}\right]\right).
\tag{61}
\]

Distribute the final minus sign:
\[
\tilde\mu_k^{DH}
=
(-C\kappa S)
+
\left(-\frac{C}{2}\kappa T\right)
+
\left(-C\kappa z_k^2\left[\chi_k + \frac{T}{2Q}\right]\right)
+
C\kappa\left[S+\frac{T}{2}\right].
\tag{62}
\]

Now expand the last term:
\[
C\kappa\left[S+\frac{T}{2}\right] = C\kappa S + \frac{C}{2}\kappa T.
\tag{63}
\]

Insert (63) into (62):
\[
\tilde\mu_k^{DH}
=
(-C\kappa S)
+
\left(-\frac{C}{2}\kappa T\right)
+
\left(-C\kappa z_k^2\left[\chi_k + \frac{T}{2Q}\right]\right)
+
(C\kappa S)
+
\left(\frac{C}{2}\kappa T\right).
\tag{64}
\]

Now cancel like terms:

- \((-C\kappa S) + (C\kappa S) = 0\)
- \(\left(-\frac{C}{2}\kappa T\right) + \left(\frac{C}{2}\kappa T\right) = 0\)

Thus:
\[
\boxed{
\tilde\mu_k^{DH}
=
-\,C\kappa z_k^2\left[\chi_k + \frac{T}{2Q}\right]
}
\tag{65}
\]

Now substitute \(C=\frac{e^2}{12\pi k_B T \epsilon}\):
\[
\tilde\mu_k^{DH}
=
-\frac{e^2}{12\pi k_B T\epsilon}\kappa z_k^2\left[\chi_k + \frac{T}{2Q}\right].
\tag{66}
\]

Rewrite by multiplying inside bracket by 2 and dividing prefactor by 2:
\[
\tilde\mu_k^{DH}
=
-\frac{e^2}{24\pi k_B T\epsilon}\kappa z_k^2\left[2\chi_k + \frac{T}{Q}\right].
\tag{67}
\]

Now insert definitions \(T=\sum_i x_i z_i^2\sigma_i\) and \(Q=\sum_i x_i z_i^2\):
\[
\boxed{
\tilde\mu_k^{DH}
=
-\frac{e^2 z_k^2\,\kappa}{24\pi k_B T\epsilon}
\left[
2\chi_k
+
\frac{\sum_i x_i z_i^2 \sigma_i}{\sum_i x_i z_i^2}
\right]
}
\tag{68}
\]

---

## 9) Convert to the \(q_i\) notation used in the final formula

Define the charge magnitude in Coulombs:
\[
q_i \equiv e z_i \quad\Rightarrow\quad q_i^2 = e^2 z_i^2.
\tag{69}
\]

Then:
- \(e^2 z_k^2 = q_k^2\)
- \(\sum_i x_i z_i^2 = \sum_i x_i \frac{q_i^2}{e^2} = \frac{1}{e^2}\sum_i x_i q_i^2\)
- \(\sum_i x_i z_i^2\sigma_i = \frac{1}{e^2}\sum_i x_i q_i^2\sigma_i\)

Therefore the ratio is unchanged:
\[
\frac{\sum_i x_i z_i^2\sigma_i}{\sum_i x_i z_i^2}
=
\frac{\sum_i x_i q_i^2\sigma_i}{\sum_i x_i q_i^2}.
\tag{70}
\]

So (68) becomes:
\[
\boxed{
\tilde{\mu}^{DH}_k
=
-\frac{q_k^2 \kappa}{24 \pi k_B T\epsilon}
\left[
2\chi_k
+
\frac{\sum_i x_i q_i^2\sigma_i}{\sum_i x_i q_i^2}
\right]
}
\tag{71}
\]

This is the desired closed form.

---

## 10) Derivation of \(\sigma_i = -2\chi_i + \frac{3}{1+\kappa a_i}\) (no skipped steps)

Start from the definition:
\[
\sigma_i \equiv \frac{\partial(\kappa\chi_i)}{\partial\kappa}.
\tag{72}
\]

From the given \(\chi_i(\kappa)\) expression, define:
\[
y \equiv \kappa a_i.
\tag{73}
\]

Let
\[
F(y) \equiv \frac{3}{2}+\ln(1+y)-2(1+y)+\frac{1}{2}(1+y)^2.
\tag{74}
\]

Then:
\[
\chi_i = \frac{3}{y^3}F(y).
\tag{75}
\]

We want:
\[
\sigma_i = \frac{d}{d\kappa}\left(\kappa\chi_i\right).
\tag{76}
\]

Use product rule:
\[
\frac{d}{d\kappa}(\kappa\chi_i) = \chi_i + \kappa\frac{d\chi_i}{d\kappa}.
\tag{77}
\]

So:
\[
\sigma_i = \chi_i + \kappa\frac{d\chi_i}{d\kappa}.
\tag{78}
\]

Now apply chain rule:
\[
\frac{d\chi_i}{d\kappa} = \frac{d\chi_i}{dy}\frac{dy}{d\kappa}.
\tag{79}
\]

Since \(y=\kappa a_i\) and \(a_i\) is constant:
\[
\frac{dy}{d\kappa}=a_i.
\tag{80}
\]

Thus:
\[
\frac{d\chi_i}{d\kappa} = a_i\frac{d\chi_i}{dy}.
\tag{81}
\]

Insert into (78):
\[
\sigma_i = \chi_i + \kappa a_i \frac{d\chi_i}{dy} = \chi_i + y\frac{d\chi_i}{dy}.
\tag{82}
\]

Now compute \(d\chi_i/dy\). From (75):
\[
\chi_i = 3F(y)\,y^{-3}.
\tag{83}
\]

Differentiate w.r.t. \(y\):
\[
\frac{d\chi_i}{dy}
=
3\left(F'(y)y^{-3} + F(y)\frac{d}{dy}(y^{-3})\right).
\tag{84}
\]

Compute:
\[
\frac{d}{dy}(y^{-3}) = -3y^{-4}.
\tag{85}
\]

Insert into (84):
\[
\frac{d\chi_i}{dy}
=
3\left(F'(y)y^{-3} - 3F(y)y^{-4}\right)
=
\frac{3F'(y)}{y^3} - \frac{9F(y)}{y^4}.
\tag{86}
\]

Multiply by \(y\):
\[
y\frac{d\chi_i}{dy}
=
\frac{3F'(y)}{y^2} - \frac{9F(y)}{y^3}.
\tag{87}
\]

Now insert into (82):
\[
\sigma_i
=
\chi_i + \left(\frac{3F'(y)}{y^2} - \frac{9F(y)}{y^3}\right).
\tag{88}
\]

Replace \(\chi_i = \frac{3F(y)}{y^3}\):
\[
\sigma_i
=
\frac{3F(y)}{y^3} + \frac{3F'(y)}{y^2} - \frac{9F(y)}{y^3}
=
\frac{3F'(y)}{y^2} + \left(\frac{3F(y)}{y^3}-\frac{9F(y)}{y^3}\right)
=
\frac{3F'(y)}{y^2} - \frac{6F(y)}{y^3}.
\tag{89}
\]

Note:
\[
\frac{6F(y)}{y^3} = 2\left(\frac{3F(y)}{y^3}\right)=2\chi_i.
\tag{90}
\]

So:
\[
\sigma_i = \frac{3F'(y)}{y^2} - 2\chi_i.
\tag{91}
\]

Now compute \(F'(y)\) from (74):

- \(\frac{d}{dy}\left(\frac32\right)=0\)
- \(\frac{d}{dy}\ln(1+y)=\frac{1}{1+y}\)
- \(\frac{d}{dy}[-2(1+y)]=-2\)
- \(\frac{d}{dy}\left[\frac12(1+y)^2\right]=(1+y)\)

Thus:
\[
F'(y)=\frac{1}{1+y}-2+(1+y).
\tag{92}
\]

Simplify \(-2+(1+y)=y-1\):
\[
F'(y)=\frac{1}{1+y}+y-1.
\tag{93}
\]

Write over the common denominator \(1+y\):
\[
F'(y)=\frac{1}{1+y} + \frac{(y-1)(1+y)}{1+y}
=
\frac{1 + (y-1)(1+y)}{1+y}.
\tag{94}
\]

Expand \((y-1)(1+y)=y(1+y)-(1+y)=y+y^2-1-y=y^2-1\). So:
\[
F'(y)=\frac{1 + (y^2-1)}{1+y}=\frac{y^2}{1+y}.
\tag{95}
\]

Insert into \(\frac{3F'(y)}{y^2}\):
\[
\frac{3F'(y)}{y^2}=\frac{3}{1+y}.
\tag{96}
\]

So (91) becomes:
\[
\sigma_i = -2\chi_i + \frac{3}{1+y}.
\tag{97}
\]

Finally use \(y=\kappa a_i\):
\[
\boxed{
\sigma_i = -2\chi_i + \frac{3}{1+\kappa a_i}
}
\tag{98}
\]

This matches the given expression.

---

## 11) Final “implementation-ready” formulas

### 11.1 Core definitions
\[
Q=\sum_i x_i z_i^2,
\quad
\kappa=\sqrt{\frac{\rho e^2}{k_B T\epsilon}Q}.
\tag{99}
\]

For each ion \(i\):
\[
\chi_i=\frac{3}{(\kappa a_i)^3}\left[\frac{3}{2}+\ln(1+\kappa a_i)-2(1+\kappa a_i)+\frac{1}{2}(1+\kappa a_i)^2\right],
\tag{100}
\]
\[
\sigma_i=-2\chi_i+\frac{3}{1+\kappa a_i}.
\tag{101}
\]

Sums:
\[
S=\sum_i x_i z_i^2\chi_i,
\qquad
T=\sum_i x_i z_i^2\sigma_i.
\tag{102}
\]

### 11.2 Energies
Let \(C=\frac{e^2}{12\pi k_BT\epsilon}\). Then:
\[
\tilde a^{DH} = -C\kappa S,
\qquad
Z^{DH}=-\frac{C}{2}\kappa T.
\tag{103}
\]

### 11.3 Derivative pieces (under constant \(\epsilon\))
\[
\frac{\partial\tilde a^{DH}}{\partial x_k}
=
-\,C\,\kappa z_k^2\left[\chi_k + \frac{T}{2Q}\right],
\tag{104}
\]
\[
\sum_j x_j\frac{\partial\tilde a^{DH}}{\partial x_j}
=
-\,C\,\kappa\left(S+\frac{T}{2}\right).
\tag{105}
\]

### 11.4 Chemical potential (two equivalent ways)

**(A) Composite identity**
\[
\tilde{\mu}^{DH}_k = \tilde a^{DH} + Z^{DH} + \frac{\partial\tilde a^{DH}}{\partial x_k}
-\sum_j x_j\frac{\partial\tilde a^{DH}}{\partial x_j}.
\tag{106}
\]

**(B) Closed form**
\[
\tilde{\mu}^{DH}_k
=
-\frac{e^2 z_k^2\,\kappa}{24\pi k_BT\epsilon}
\left[
2\chi_k + \frac{T}{Q}
\right].
\tag{107}
\]

Or with \(q_i=e z_i\):
\[
\tilde{\mu}^{DH}_k
=
-\frac{q_k^2\,\kappa}{24\pi k_BT\epsilon}
\left[
2\chi_k + \frac{\sum_i x_i q_i^2\sigma_i}{\sum_i x_i q_i^2}
\right].
\tag{108}
\]

---

## 12) Notes for Codex implementation

- Keep a consistent convention for \(\rho\):
  - If your \(\kappa\) formula uses **number density**, set \(\rho = \rho_{\text{molar}}N_A\).
  - If it uses molar density, ensure the prefactor matches the equation set you’re using.
- The equivalence of (106) and (108) holds only under the same assumptions used in the derivation (especially constant \(\epsilon\) vs \(\epsilon(\mathbf x)\)).
- If you later enable \(\epsilon(\mathbf x)\), you must include \(\partial\epsilon/\partial x_k\) contributions in:
  - \(\partial\kappa/\partial x_k\),
  - \(\partial\tilde a^{DH}/\partial x_k\),
  - and the final closed form changes (no longer exactly (108)) unless you derive the modified closed form.

---
