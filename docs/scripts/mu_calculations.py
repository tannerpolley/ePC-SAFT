#!/usr/bin/env python3
"""
dh_mu_equivalence_demo.py

Proves (numerically) that the Debye–Hückel chemical potential contribution mu~_k^DH
computed by:

(1) The closed-form final expression
    mu~_k^DH = -(q_k^2*kappa)/(24*pi*kB*T*eps) * [ 2*chi_k + (Σ x_i q_i^2 sigma_i)/(Σ x_i q_i^2) ]

matches exactly the value computed by:

(2) The assembled definition
    mu~_k^DH = a~^DH + Z^DH + (∂a~^DH/∂x_k) - Σ_j x_j (∂a~^DH/∂x_j)

for the NaCl(aq) example using your provided parameters.

Notes / assumptions:
- eps_r is treated as constant w.r.t. composition (as required for the closed-form identity).
- Ion "a_i" used inside chi_i is taken from the provided segment diameters s_i (Å -> m).
- kappa is computed using number density: rho_number = rho_molar * N_A.
- mu~ is dimensionless (mu/(kB*T)).
"""

from __future__ import annotations

import numpy as np


def main() -> None:
    # -----------------------------
    # Physical constants (SI)
    # -----------------------------
    e_charge = 1.602176634e-19  # C
    kB = 1.380649e-23  # J/K
    eps0 = 8.8541878128e-12  # F/m
    NA = 6.02214076e23  # 1/mol
    pi = np.pi

    # -----------------------------
    # User-provided parameters
    # 0 = Na+, 1 = Cl-, 2 = H2O
    # -----------------------------
    x = np.asarray([0.0629838206, 0.0629838206, 0.8740323588], dtype=float)
    m = np.asarray([1.0, 1.0, 1.2047], dtype=float)  # not used here, included for completeness
    s = np.asarray([2.8232, 2.7599589, 0.0], dtype=float)  # segment diameters [Å]
    E = np.asarray([230.00, 170.00, 353.9449], dtype=float)  # not used here
    volAB = np.asarray([0.0, 0.0, 0.0451], dtype=float)  # not used here
    eAB = np.asarray([0.0, 0.0, 2425.67], dtype=float)  # not used here
    k_ij = np.asarray([[0.0, 0.317, 0.0], [0.317, 0.0, -0.25], [0.0, -0.25, 0.0]], dtype=float)  # not used here
    z = np.asarray([1.0, -1.0, 0.0], dtype=float)

    ref = 1.116  # not used here
    T = 293.15  # K

    # Temperature-dependent segment diameter for water (Å)
    s[2] = 2.7927 + 10.11 * np.exp(-0.01775 * T) - 1.417 * np.exp(-0.01146 * T)

    # Temperature-dependent k_ij entries (not used here)
    k_ij[0, 2] = -0.007981 * T + 2.37999
    k_ij[2, 0] = -0.007981 * T + 2.37999

    # Dielectric constant of water (user-specified)
    eps_r = 78.09
    eps = eps0 * eps_r

    # Given molar density (mol/m^3) from your earlier message
    rho_molar = 55757.0726  # mol/m^3
    rho_number = rho_molar * NA  # 1/m^3

    # -----------------------------
    # DH setup: ions only
    # -----------------------------
    ions = [0, 1]  # Na+, Cl-
    ion_names = ["Na+", "Cl-"]

    x_ions = x[ions]
    z_ions = z[ions]

    # Q = Σ_ions x_i z_i^2
    Q = float(np.sum(x_ions * z_ions**2))

    # Use s_i (Å) as a_i in chi_i, converted to meters
    a_ion = s[ions] * 1e-10  # m

    # Screening parameter kappa (number-density formulation)
    kappa = float(np.sqrt(rho_number * e_charge**2 / (kB * T * eps) * Q))

    # -----------------------------
    # chi_i and sigma_i
    # -----------------------------
    def chi(kappa_: float, a_i: float) -> float:
        y = kappa_ * a_i
        return float((3.0 / (y**3)) * (1.5 + np.log(1.0 + y) - 2.0 * (1.0 + y) + 0.5 * (1.0 + y) ** 2))

    chi_i = np.asarray([chi(kappa, a_ion[0]), chi(kappa, a_ion[1])], dtype=float)

    # sigma_i = ∂(kappa*chi_i)/∂kappa = -2 chi_i + 3/(1 + kappa a_i)
    sigma_i = -2.0 * chi_i + 3.0 / (1.0 + kappa * a_ion)

    # -----------------------------
    # a~^DH and Z^DH
    # -----------------------------
    C = e_charge**2 / (12.0 * pi * kB * T * eps)  # constant used in a_DH and Z_DH

    # S = Σ x_i z_i^2 chi_i ;  Tsum = Σ x_i z_i^2 sigma_i
    S = float(np.sum(x_ions * z_ions**2 * chi_i))
    Tsum = float(np.sum(x_ions * z_ions**2 * sigma_i))

    # a~^DH = -C * kappa * S
    a_DH = -C * kappa * S

    # Z^DH = -(C/2) * kappa * Tsum
    Z_DH = -(C / 2.0) * kappa * Tsum

    # -----------------------------
    # Method 1: closed-form mu~_k^DH
    # -----------------------------
    def mu_direct(idx_ion: int) -> float:
        # q_i = e*z_i => q_i^2 = (e*z_i)^2
        q2_k = (e_charge * z[ions[idx_ion]]) ** 2

        numer = float(np.sum(x_ions * (e_charge * z_ions) ** 2 * sigma_i))
        denom = float(np.sum(x_ions * (e_charge * z_ions) ** 2))

        return float(-(q2_k * kappa) / (24.0 * pi * kB * T * eps) * (2.0 * chi_i[idx_ion] + numer / denom))

    # -----------------------------
    # Method 2: assembled mu~_k^DH
    # Using the analytical derivative derived earlier:
    #   ∂a~^DH/∂x_k = -C * kappa * z_k^2 * ( chi_k + Tsum/(2Q) )
    #   Σ_j x_j ∂a~^DH/∂x_j = -C * kappa * ( S + Tsum/2 )
    # -----------------------------
    def da_dx(idx_ion: int) -> float:
        return float(-C * kappa * (z_ions[idx_ion] ** 2) * (chi_i[idx_ion] + Tsum / (2.0 * Q)))

    sum_x_dadx = float(-C * kappa * (S + Tsum / 2.0))

    def mu_composite(idx_ion: int) -> float:
        return float(a_DH + Z_DH + da_dx(idx_ion) - sum_x_dadx)

    # -----------------------------
    # Print report
    # -----------------------------
    print("=== Inputs used ===")
    print(f"T [K] = {T}")
    print(f"eps_r = {eps_r}")
    print(f"rho_molar [mol/m^3] = {rho_molar}")
    print(f"x (Na+, Cl-, H2O) = {x}")
    print(f"s (Å) (Na+, Cl-, H2O) = {s}")
    print(f"a_ion (m) (Na+, Cl-) = {a_ion}")

    print("\n=== Intermediate DH quantities ===")
    print(f"Q = sum_ions x_i z_i^2 = {Q}")
    print(f"kappa [1/m] = {kappa}")
    print(f"chi (Na+, Cl-) = {chi_i}")
    print(f"sigma (Na+, Cl-) = {sigma_i}")
    print(f"a~^DH = {a_DH}")
    print(f"Z^DH = {Z_DH}")

    print("\n=== mu~ comparison (dimensionless) ===")
    for idx, name in enumerate(ion_names):
        m1 = mu_direct(idx)
        m2 = mu_composite(idx)
        diff = m2 - m1
        print(f"{name:3s}  direct={m1:.15g}  composite={m2:.15g}  diff={diff:.3g}")

    # Optional: assert near-equality
    for idx in range(2):
        if not np.isclose(mu_direct(idx), mu_composite(idx), rtol=0, atol=1e-12):
            raise AssertionError("Methods do not match within tolerance!")


if __name__ == "__main__":
    main()
