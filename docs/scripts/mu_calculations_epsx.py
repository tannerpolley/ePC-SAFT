#!/usr/bin/env python3
"""
dh_mu_4ways_eps0_epsr.py

Four DH mu~ computation paths:

V1  Closed-form (assumes depsr_dx = 0)
V2  Composite constant-epsr simplified derivative
V3  Composite with epsr(x) terms retained, evaluated with depsr_dx = 0 (should match V1/V2)
V4  Same as V3 but with depsr_dx = 8 (sensitivity test)

All mu~ values are dimensionless: mu/(kB*T).

Important:
- epsilon = eps0 * eps_r
- Only eps_r is treated as composition-dependent. eps0 is constant.
- kappa uses number density: rho_number = rho_molar * N_A
"""

import numpy as np

# --- Constants (SI) ---
e_charge = 1.602176634e-19
kB = 1.380649e-23
eps0 = 8.8541878128e-12
NA = 6.02214076e23
pi = np.pi


def chi(kappa: float, a_i: float) -> float:
    y = kappa * a_i
    return (3.0 / (y**3)) * (1.5 + np.log(1.0 + y) - 2.0 * (1.0 + y) + 0.5 * (1.0 + y) ** 2)


def sigma_from_chi(kappa: float, a_i: float, chi_i: float) -> float:
    return -2.0 * chi_i + 3.0 / (1.0 + kappa * a_i)


def main() -> None:
    # -------------------------
    # User case parameters
    # 0 = Na+, 1 = Cl-, 2 = H2O
    # -------------------------
    x_full = np.asarray([0.0629838206, 0.0629838206, 0.8740323588], dtype=float)
    z_full = np.asarray([1.0, -1.0, 0.0], dtype=float)

    T = 293.15
    eps_r = 78.09  # user-specified dielc_water
    eps = eps0 * eps_r

    rho_molar = 55757.0726  # mol/m^3
    rho_number = rho_molar * NA  # 1/m^3

    # segment diameters [Å]
    s = np.asarray([2.8232, 2.7599589, 0.0], dtype=float)
    s[2] = 2.7927 + 10.11 * np.exp(-0.01775 * T) - 1.417 * np.exp(-0.01146 * T)

    # ions only
    ions = [0, 1]
    x = x_full[ions]
    z = z_full[ions]
    a = s[ions] * 1e-10  # m

    # -------------------------
    # Core DH quantities
    # -------------------------
    Q = np.sum(x * z**2)
    kappa = np.sqrt(rho_number * e_charge**2 / (kB * T * eps) * Q)

    chi_i = np.array([chi(kappa, a[0]), chi(kappa, a[1])], dtype=float)
    sigma_i = np.array(
        [
            sigma_from_chi(kappa, a[0], chi_i[0]),
            sigma_from_chi(kappa, a[1], chi_i[1]),
        ],
        dtype=float,
    )

    S = np.sum(x * z**2 * chi_i)
    Tsum = np.sum(x * z**2 * sigma_i)

    # constants
    C = e_charge**2 / (12.0 * pi * kB * T * eps)  # includes 1/(eps0*epsr)
    K0 = e_charge**2 / (12.0 * pi * kB * T)  # excludes dielectric

    # -------------------------
    # V1: closed-form (valid when depsr_dx = 0)
    # -------------------------
    def mu_closed(idx: int) -> float:
        q2 = (e_charge * z[idx]) ** 2
        numer = np.sum(x * (e_charge * z) ** 2 * sigma_i)
        denom = np.sum(x * (e_charge * z) ** 2)
        return -(q2 * kappa) / (24.0 * pi * kB * T * eps) * (2.0 * chi_i[idx] + numer / denom)

    # -------------------------
    # V2: composite constant-epsr simplified derivative
    # -------------------------
    a_DH = -C * kappa * S
    Z_DH = -(C / 2.0) * kappa * Tsum

    def da_dx_const(idx: int) -> float:
        return -C * kappa * (z[idx] ** 2) * (chi_i[idx] + Tsum / (2.0 * Q))

    sum_x_dadx_const = -C * kappa * (S + Tsum / 2.0)

    def mu_comp_const(idx: int) -> float:
        return a_DH + Z_DH + da_dx_const(idx) - sum_x_dadx_const

    # -------------------------
    # V3/V4: composite with epsr(x) retained (depsr_dx selectable)
    # -------------------------
    # kappa^2 = A * Q / (eps0 * epsr)
    A = rho_number * e_charge**2 / (kB * T)

    # dchi/dkappa from sigma = chi + kappa dchi/dkappa
    dchi_dk = (sigma_i - chi_i) / kappa

    def mu_comp_epsr(idx: int, depsr_dx: np.ndarray) -> float:
        """
        depsr_dx: array([d eps_r / dx_Na, d eps_r / dx_Cl])  (per mole fraction)
        """

        def dkappa_dx(i: int) -> float:
            dQ = z[i] ** 2
            # 2*kappa*dk = A*( dQ/(eps0*epsr) - Q*(eps0*depsr_dx)/(eps0*epsr)^2 )
            return (A * (dQ / (eps0 * eps_r) - Q * (eps0 * depsr_dx[i]) / (eps0 * eps_r) ** 2)) / (2.0 * kappa)

        def d_inv_eps_dx(i: int) -> float:
            # d(1/(eps0*epsr))/dx = -(1/(eps0))* (depsr_dx/epsr^2)
            return -(1.0 / eps0) * (depsr_dx[i] / (eps_r**2))

        def dS_dx(i: int) -> float:
            return z[i] ** 2 * chi_i[i] + dkappa_dx(i) * np.sum(x * z**2 * dchi_dk)

        def a_DH_gen() -> float:
            return -K0 * (kappa / (eps0 * eps_r)) * S

        def Z_DH_gen() -> float:
            return -(K0 / 2.0) * (kappa / (eps0 * eps_r)) * Tsum

        def da_dx_gen(i: int) -> float:
            # a = -K0 * (kappa/(eps0*epsr)) * S
            # d[(kappa/(eps0*epsr))S] = (d(kappa/(eps0*epsr))*S) + (kappa/(eps0*epsr))*dS
            term1 = (((1.0 / (eps0 * eps_r)) * dkappa_dx(i)) + kappa * d_inv_eps_dx(i)) * S
            term2 = (kappa / (eps0 * eps_r)) * dS_dx(i)
            return -K0 * (term1 + term2)

        sum_x_dadx = np.sum(x * np.array([da_dx_gen(0), da_dx_gen(1)]))
        return a_DH_gen() + Z_DH_gen() + da_dx_gen(idx) - sum_x_dadx

    # V3: depsr_dx = 0
    depsr_dx_0 = np.array([0.0, 0.0], dtype=float)

    # V4: depsr_dx = 8 (your test)
    depsr_dx_8 = np.array([8.0, 8.0], dtype=float)

    # -------------------------
    # Print report
    # -------------------------
    print("=== Inputs ===")
    print("T [K] =", T)
    print("eps_r =", eps_r)
    print("rho_molar [mol/m^3] =", rho_molar)
    print("x (Na+, Cl-, H2O) =", x_full)
    print("s (Å) (Na+, Cl-, H2O) =", s)
    print()

    print("=== Intermediate DH quantities ===")
    print("Q =", Q)
    print("kappa [1/m] =", kappa)
    print("chi (Na+, Cl-) =", chi_i)
    print("sigma (Na+, Cl-) =", sigma_i)
    print()

    print("=== mu~ comparison (dimensionless) ===")
    names = ["Na+", "Cl-"]
    for i, nm in enumerate(names):
        v1 = mu_closed(i)
        v2 = mu_comp_const(i)
        v3 = mu_comp_epsr(i, depsr_dx_0)
        v4 = mu_comp_epsr(i, depsr_dx_8)
        print(
            f"{nm:3s}  V1(closed)={v1:.15g}  V2(comp_const)={v2:.15g}  V3(comp_epsr,deps=0)={v3:.15g}  V4(comp_epsr,deps=8)={v4:.15g}"
        )
        print(f"     diffs: V2-V1={v2-v1:.3g}, V3-V1={v3-v1:.3g}, V4-V1={v4-v1:.6g}")

    print("\nNote: depsr_dx has units of 'per mole fraction' (eps_r is dimensionless).")


if __name__ == "__main__":
    main()
