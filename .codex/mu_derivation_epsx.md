# Debye–Hückel (DH) chemical potential — composite derivation with ε0 and εr(x)

We derive the composite mixture DH chemical potential when the dielectric constant is written as
\[
\epsilon(\mathbf{x})=\varepsilon_0\,\varepsilon_r(\mathbf{x}),
\]
where \(\varepsilon_0\) is constant and only \(\varepsilon_r\) depends on composition.

We keep \(d\varepsilon_r/dx_k\) symbolic throughout. If you later set it to zero, you must recover the constant-\(\varepsilon_r\) results.

---

## 0) Starting equations

### 0.1 Reduced Helmholtz energy (dimensionless)
\[
\tilde{a}^{DH}
=
-\frac{\kappa e^2}{12\pi k_BT\,\varepsilon_0\,\varepsilon_r}\sum_i x_i z_i^2\chi_i
\tag{1}
\]

### 0.2 Screening parameter
\[
\kappa
=
\sqrt{\frac{\rho e^2}{k_BT\,\varepsilon_0\,\varepsilon_r}\sum_i x_i z_i^2}
\tag{2}
\]

### 0.3 Compressibility contribution
\[
Z^{DH}
=
-\frac{\kappa e^2}{24\pi k_BT\,\varepsilon_0\,\varepsilon_r}\sum_i x_i z_i^2\sigma_i
\tag{3}
\]

### 0.4 Sigma definition
\[
\sigma_i
\equiv
\left(\frac{\partial(\kappa\chi_i)}{\partial\kappa}\right)_{T,\mathbf{x}}
=
\chi_i+\kappa\frac{d\chi_i}{d\kappa}
\tag{4}
\]

### 0.5 Mixture chemical potential identity
\[
\tilde{\mu}^{DH}_k
=
\tilde a^{DH}
+
Z^{DH}
+
\left(\frac{\partial \tilde a^{DH}}{\partial x_k}\right)_{T,\nu,x_{i\neq k}}
-
\sum_j x_j
\left(\frac{\partial \tilde a^{DH}}{\partial x_j}\right)_{T,\nu,x_{i\neq j}}
\tag{5}
\]

---

## 1) Notation

Define
\[
K_0 \equiv \frac{e^2}{12\pi k_BT},
\qquad
Q \equiv \sum_i x_i z_i^2,
\qquad
S \equiv \sum_i x_i z_i^2\chi_i,
\qquad
T \equiv \sum_i x_i z_i^2\sigma_i.
\tag{6}
\]

Then
\[
\tilde a^{DH} = -K_0\left(\frac{\kappa}{\varepsilon_0\varepsilon_r}\right)S,
\tag{7}
\]
\[
Z^{DH} = -\frac{K_0}{2}\left(\frac{\kappa}{\varepsilon_0\varepsilon_r}\right)T.
\tag{8}
\]

Also define
\[
A \equiv \frac{\rho e^2}{k_BT},
\tag{9}
\]
so Eq. (2) becomes
\[
\kappa^2 = \frac{A}{\varepsilon_0\varepsilon_r}Q.
\tag{10}
\]

---

## 2) Derivative of κ with εr(x)

Differentiate (10) with respect to \(x_k\):

\[
\frac{\partial(\kappa^2)}{\partial x_k}
=
\frac{\partial}{\partial x_k}\left(\frac{A}{\varepsilon_0\varepsilon_r}Q\right).
\tag{11}
\]

Left side:
\[
\frac{\partial(\kappa^2)}{\partial x_k} = 2\kappa\frac{\partial\kappa}{\partial x_k}.
\tag{12}
\]

Right side: \(A\) and \(\varepsilon_0\) are constants, so treat it as \(A\cdot (1/\varepsilon_0)\cdot (Q/\varepsilon_r)\):
\[
\frac{\partial}{\partial x_k}\left(\frac{A}{\varepsilon_0}\frac{Q}{\varepsilon_r}\right)
=
\frac{A}{\varepsilon_0}
\frac{\partial}{\partial x_k}\left(\frac{Q}{\varepsilon_r}\right).
\tag{13}
\]

Product rule on \(Q\cdot (1/\varepsilon_r)\):
\[
\frac{\partial}{\partial x_k}\left(\frac{Q}{\varepsilon_r}\right)
=
\left(\frac{\partial Q}{\partial x_k}\right)\frac{1}{\varepsilon_r}
+
Q\frac{\partial}{\partial x_k}\left(\frac{1}{\varepsilon_r}\right).
\tag{14}
\]

Now:
\[
Q=\sum_i x_i z_i^2
\quad\Rightarrow\quad
\frac{\partial Q}{\partial x_k}=z_k^2
\tag{15}
\]
(holding \(x_{i\neq k}\) fixed).

And
\[
\frac{\partial}{\partial x_k}\left(\frac{1}{\varepsilon_r}\right)
=
-\frac{1}{\varepsilon_r^2}\frac{\partial\varepsilon_r}{\partial x_k}.
\tag{16}
\]

Insert (15) and (16) into (14):
\[
\frac{\partial}{\partial x_k}\left(\frac{Q}{\varepsilon_r}\right)
=
\frac{z_k^2}{\varepsilon_r}
-
\frac{Q}{\varepsilon_r^2}\frac{\partial\varepsilon_r}{\partial x_k}.
\tag{17}
\]

Insert (17) into (13):
\[
\frac{\partial}{\partial x_k}\left(\frac{A}{\varepsilon_0\varepsilon_r}Q\right)
=
\frac{A}{\varepsilon_0}\left[
\frac{z_k^2}{\varepsilon_r}
-
\frac{Q}{\varepsilon_r^2}\frac{\partial\varepsilon_r}{\partial x_k}
\right].
\tag{18}
\]

Equate with (12):
\[
2\kappa\frac{\partial\kappa}{\partial x_k}
=
\frac{A}{\varepsilon_0}\left[
\frac{z_k^2}{\varepsilon_r}
-
\frac{Q}{\varepsilon_r^2}\frac{\partial\varepsilon_r}{\partial x_k}
\right].
\tag{19}
\]

Solve:
\[
\boxed{
\frac{\partial\kappa}{\partial x_k}
=
\frac{A}{2\kappa\varepsilon_0}\left[
\frac{z_k^2}{\varepsilon_r}
-
\frac{Q}{\varepsilon_r^2}\frac{\partial\varepsilon_r}{\partial x_k}
\right]
}
\tag{20}
\]

---

## 3) Derivative of S

Recall:
\[
S=\sum_i x_i z_i^2\chi_i(\kappa).
\tag{21}
\]

Differentiate:
\[
\frac{\partial S}{\partial x_k}
=
\sum_i z_i^2\frac{\partial}{\partial x_k}(x_i\chi_i)
\tag{22}
\]

Product rule:
\[
\frac{\partial}{\partial x_k}(x_i\chi_i)
=
\chi_i\frac{\partial x_i}{\partial x_k} + x_i\frac{\partial\chi_i}{\partial x_k}.
\tag{23}
\]

Using \(\partial x_i/\partial x_k=\delta_{ik}\):
\[
\sum_i z_i^2\chi_i\delta_{ik} = z_k^2\chi_k.
\tag{24}
\]

So:
\[
\frac{\partial S}{\partial x_k}
=
z_k^2\chi_k
+
\sum_i x_i z_i^2\frac{\partial\chi_i}{\partial x_k}.
\tag{25}
\]

Chain rule (only via \(\kappa\)):
\[
\frac{\partial\chi_i}{\partial x_k}
=
\frac{d\chi_i}{d\kappa}\frac{\partial\kappa}{\partial x_k}.
\tag{26}
\]

Thus:
\[
\frac{\partial S}{\partial x_k}
=
z_k^2\chi_k
+
\left(\frac{\partial\kappa}{\partial x_k}\right)
\sum_i x_i z_i^2\frac{d\chi_i}{d\kappa}.
\tag{27}
\]

Use (4):
\[
\sigma_i=\chi_i+\kappa\frac{d\chi_i}{d\kappa}
\quad\Rightarrow\quad
\frac{d\chi_i}{d\kappa}=\frac{\sigma_i-\chi_i}{\kappa}.
\tag{28}
\]

Then:
\[
\sum_i x_i z_i^2\frac{d\chi_i}{d\kappa}
=
\frac{1}{\kappa}\left(
\sum_i x_i z_i^2\sigma_i - \sum_i x_i z_i^2\chi_i
\right)
=
\frac{T-S}{\kappa}.
\tag{29}
\]

Insert into (27):
\[
\boxed{
\frac{\partial S}{\partial x_k}
=
z_k^2\chi_k
+
\left(\frac{\partial\kappa}{\partial x_k}\right)\frac{T-S}{\kappa}
}
\tag{30}
\]

---

## 4) Derivative of a~^DH with εr(x)

From (7):
\[
\tilde a^{DH}=-K_0\left(\frac{\kappa}{\varepsilon_0\varepsilon_r}\right)S.
\tag{31}
\]

Differentiate:
\[
\frac{\partial \tilde a^{DH}}{\partial x_k}
=
-K_0\frac{\partial}{\partial x_k}\left[
\left(\frac{\kappa}{\varepsilon_0\varepsilon_r}\right)S
\right].
\tag{32}
\]

Product rule:
\[
\frac{\partial}{\partial x_k}\left[
\left(\frac{\kappa}{\varepsilon_0\varepsilon_r}\right)S
\right]
=
\left(\frac{\partial}{\partial x_k}\frac{\kappa}{\varepsilon_0\varepsilon_r}\right)S
+
\left(\frac{\kappa}{\varepsilon_0\varepsilon_r}\right)\frac{\partial S}{\partial x_k}.
\tag{33}
\]

Now write:
\[
\frac{\kappa}{\varepsilon_0\varepsilon_r} = \kappa\left(\frac{1}{\varepsilon_0}\right)\left(\frac{1}{\varepsilon_r}\right),
\]
with \(1/\varepsilon_0\) constant, so:
\[
\frac{\partial}{\partial x_k}\left(\frac{\kappa}{\varepsilon_0\varepsilon_r}\right)
=
\frac{1}{\varepsilon_0}\frac{\partial}{\partial x_k}\left(\kappa\frac{1}{\varepsilon_r}\right).
\tag{34}
\]

Product rule:
\[
\frac{\partial}{\partial x_k}\left(\kappa\frac{1}{\varepsilon_r}\right)
=
\left(\frac{\partial\kappa}{\partial x_k}\right)\frac{1}{\varepsilon_r}
+
\kappa\frac{\partial}{\partial x_k}\left(\frac{1}{\varepsilon_r}\right).
\tag{35}
\]

Use (16):
\[
\frac{\partial}{\partial x_k}\left(\frac{1}{\varepsilon_r}\right)
=
-\frac{1}{\varepsilon_r^2}\frac{\partial\varepsilon_r}{\partial x_k}.
\tag{36}
\]

So:
\[
\frac{\partial}{\partial x_k}\left(\frac{\kappa}{\varepsilon_0\varepsilon_r}\right)
=
\frac{1}{\varepsilon_0}\left[
\left(\frac{\partial\kappa}{\partial x_k}\right)\frac{1}{\varepsilon_r}
-
\kappa\frac{1}{\varepsilon_r^2}\frac{\partial\varepsilon_r}{\partial x_k}
\right].
\tag{37}
\]

Insert (37) into (33), then into (32):

\[
\boxed{
\frac{\partial \tilde a^{DH}}{\partial x_k}
=
-K_0\left\{
\left[
\frac{1}{\varepsilon_0}\left(\frac{\partial\kappa}{\partial x_k}\right)\frac{1}{\varepsilon_r}
-
\frac{1}{\varepsilon_0}\kappa\frac{1}{\varepsilon_r^2}\frac{\partial\varepsilon_r}{\partial x_k}
\right]S
+
\left(\frac{\kappa}{\varepsilon_0\varepsilon_r}\right)\frac{\partial S}{\partial x_k}
\right\}
}
\tag{38}
\]

with \(\partial\kappa/\partial x_k\) from (20) and \(\partial S/\partial x_k\) from (30).

---

## 5) Composite chemical potential (final)

\[
\boxed{
\tilde{\mu}^{DH}_k
=
\tilde a^{DH}
+
Z^{DH}
+
\left(\frac{\partial \tilde a^{DH}}{\partial x_k}\right)
-
\sum_j x_j\left(\frac{\partial \tilde a^{DH}}{\partial x_j}\right)
}
\tag{39}
\]

where:
\[
\tilde a^{DH} = -K_0\left(\frac{\kappa}{\varepsilon_0\varepsilon_r}\right)S,
\qquad
Z^{DH} = -\frac{K_0}{2}\left(\frac{\kappa}{\varepsilon_0\varepsilon_r}\right)T,
\tag{40}
\]
and \(\partial\tilde a^{DH}/\partial x_k\) is (38).

---

## 6) Reduction check

If \(\partial\varepsilon_r/\partial x_k \to 0\), then (20) reduces to the constant-\(\varepsilon_r\) derivative, and (38) collapses to the constant-\(\varepsilon_r\) derivative used in the closed-form equivalence proof.

So implementing the εr(x) composite form and evaluating it with \(\partial\varepsilon_r/\partial x_k \approx 0\) must reproduce the same numerical \(\tilde\mu_k^{DH}\) as the constant-\(\varepsilon_r\) methods.
