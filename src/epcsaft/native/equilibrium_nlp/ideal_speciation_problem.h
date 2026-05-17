#pragma once

#include "ipopt_adapter.h"

#include <string>
#include <vector>

struct ChemicalEquilibriumOptionsNative;
struct ChemicalEquilibriumResultNative;

namespace epcsaft::native::equilibrium_nlp {

struct IdealSpeciationRequest {
    int species_count = 0;
    int balance_rows = 0;
    std::vector<double> balance_matrix_row_major;
    std::vector<double> total_vector;
    int reaction_rows = 0;
    std::vector<double> reaction_stoichiometry_row_major;
    std::vector<double> log_equilibrium_constants;
    std::vector<double> initial_x;
    std::vector<double> charges;
    double min_mole_fraction = 1.0e-14;
};

struct IdealSpeciationIpoptResult {
    IpoptSolveResult ipopt;
    std::vector<double> amounts;
    std::vector<double> composition;
    std::vector<double> activity_coefficients;
    std::vector<double> mass_balance_residuals;
    double charge_residual = 0.0;
    std::vector<double> reaction_residuals;
    std::vector<double> standard_mu_rt;
};

IdealSpeciationIpoptResult solve_ideal_speciation_ipopt(
    const IdealSpeciationRequest& request,
    const IpoptSolveOptions& options,
    const std::string& derivative_backend
);

ChemicalEquilibriumResultNative solve_ideal_speciation_chemical_equilibrium_ipopt(
    const IdealSpeciationRequest& request,
    const ChemicalEquilibriumOptionsNative& options
);

}  // namespace epcsaft::native::equilibrium_nlp
