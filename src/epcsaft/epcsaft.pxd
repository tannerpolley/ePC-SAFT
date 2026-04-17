# -*- coding: utf-8 -*-
# setuptools: language=c++
"""
Created on Thu Jul 19 14:23:00 2018

@author: Tanner Polley
"""
from libcpp.vector cimport vector
from libcpp.memory cimport shared_ptr
from libc.stddef cimport size_t

cdef extern from "epcsaft_electrolyte.h":
    double Z_cpp(double t, double rho, vector[double] x, add_args &cppargs)
    vector[double] mures_cpp(double t, double rho, vector[double] x, add_args &cppargs)
    vector[double] lnfug_cpp(double t, double rho, vector[double] x, add_args &cppargs)
    vector[double] fugcoef_cpp(double t, double rho, vector[double] x, add_args &cppargs)
    double den_cpp(double t, double p, vector[double] x, int phase, add_args &cppargs) except +
    double ares_cpp(double t, double rho, vector[double] x, add_args &cppargs)
    double dadt_cpp(double t, double rho, vector[double] x, add_args &cppargs)
    double sres_cpp(double t, double rho, vector[double] x, add_args &cppargs)
    double gres_cpp(double t, double rho, vector[double] x, add_args &cppargs)
    double dielectric_eps_cpp(vector[double] x, add_args &cppargs)
    vector[double] dielectric_diff_cpp(vector[double] x, add_args &cppargs)
    double dielc_eps_cpp(vector[double] x, add_args &cppargs)
    vector[double] dielc_diff_cpp(vector[double] x, add_args &cppargs)

    ctypedef struct add_args:
        vector[double] m
        vector[double] s
        vector[double] e
        vector[double] k_ij
        vector[double] e_assoc
        vector[double] vol_a
        vector[double] z
        vector[double] dielc
        vector[double] mw
        vector[double] mixed_rel_perm_a
        vector[double] mixed_rel_perm_b
        vector[double] mixed_rel_perm_c
        vector[int] mixed_rel_perm_mask
        int mixed_rel_perm_water_index
        int dielc_rule
        int dielc_diff_mode
        int hc_dadx_diff_mode
        int disp_dadx_diff_mode
        int assoc_dadx_diff_mode
        int d_ion_mode
        int mu_DH_diff_mode
        int mu_DH_comp_dep_rel_perm
        int mu_DH_include_sum_term
        int include_born_model
        int d_born_mode
        int born_solvation_shell_model
        int born_dielectric_saturation
        int born_bulk_mode
        int mu_born_diff_mode
        int mu_born_comp_dep_rel_perm
        int mu_born_include_sum_term
        int mu_born_comp_dep_delta_d
        vector[double] d_born
        vector[double] f_solv
        int born_model
        int born_radius_model
        int born_diff_mode
        int born_eps_mode
        int DH_model
        vector[int] assoc_num
        vector[int] assoc_matrix
        vector[double] k_hb
        vector[double] l_ij

    ctypedef struct ActivityCoefficientNative:
        vector[double] component_activity_coefficients
        vector[double] mean_ionic_activity_coefficients_mole_fraction
        vector[double] mean_ionic_activity_coefficients_molality
        vector[double] solvation_free_energy
        vector[double] pair_molality
        vector[double] pair_conversion_factor
        vector[int] cation_indices
        vector[int] anion_indices
        vector[int] solvent_indices
        vector[int] pair_cation_indices
        vector[int] pair_anion_indices
        vector[int] pair_nu_cation
        vector[int] pair_nu_anion
        int solvent_index
        double osmotic_coefficient

    ctypedef struct ScalarContributionTerms:
        double hc
        double disp
        double assoc
        double ion
        double born
        double total

    ctypedef struct CompressibilityFactorResult:
        ScalarContributionTerms raw
        ScalarContributionTerms terms

    ctypedef struct VectorContributionTerms:
        vector[double] hc
        vector[double] disp
        vector[double] assoc
        vector[double] ion
        vector[double] born
        vector[double] total

    ctypedef struct CompositionContributionResult:
        VectorContributionTerms dadx
        ScalarContributionTerms ares
        ScalarContributionTerms sum_x_dadx
        ScalarContributionTerms z_raw
        ScalarContributionTerms z

    ctypedef struct ResidualChemicalPotentialResult:
        VectorContributionTerms mu
        CompositionContributionResult composition

    ctypedef struct FugacityContributionResult:
        VectorContributionTerms mu
        VectorContributionTerms lnfugcoef
        CompositionContributionResult composition

    cdef cppclass ePCSAFTStateNative:
        ePCSAFTStateNative(shared_ptr[ePCSAFTMixtureNative] mixture, double t, vector[double] x,
                          int phase, bint has_p, double p, bint has_rho, double rho) except +
        double temperature() const
        int phase() const
        const vector[double]& composition() const
        double pressure()
        double density()
        double compressibility_factor()
        CompressibilityFactorResult compressibility_factor_result()
        double residual_helmholtz()
        ScalarContributionTerms residual_helmholtz_result()
        double temperature_derivative_residual_helmholtz()
        ScalarContributionTerms temperature_derivative_residual_helmholtz_result()
        double residual_enthalpy()
        double residual_entropy()
        double residual_gibbs()
        vector[double] residual_chemical_potential()
        ResidualChemicalPotentialResult residual_chemical_potential_result()
        CompositionContributionResult composition_derivative_residual_helmholtz_result()
        vector[double] ln_fugacity_coefficient()
        vector[double] fugacity_coefficient()
        FugacityContributionResult fugacity_coefficient_result()
        vector[double] relative_permittivity()
        double osmotic_coefficient()
        vector[double] solvation_free_energy()
        ActivityCoefficientNative activity_coefficient_native(bint include_aux, bint has_solvent_override, int solvent_override_index) except +

    cdef cppclass ePCSAFTMixtureNative:
        ePCSAFTMixtureNative(const add_args& args) except +
        const add_args& args() const
        shared_ptr[ePCSAFTStateNative] state(double t, vector[double] x, int phase,
                                           bool has_p, double p, bool has_rho, double rho)
        size_t ncomp() const
        void clear_runtime_caches()
        void reset_runtime_cache_stats()
        size_t reference_state_cache_hits() const
        size_t reference_state_cache_misses() const
        size_t density_warm_start_hits() const
        size_t density_warm_start_fallbacks() const



