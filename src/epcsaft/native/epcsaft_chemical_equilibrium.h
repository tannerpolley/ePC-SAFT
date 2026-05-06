#pragma once

#include <map>
#include <memory>
#include <string>
#include <vector>

class ePCSAFTMixtureNative;

struct ChemicalEquilibriumOptionsNative {
    int max_iterations = 50;
    double tolerance = 1.0e-8;
    double damping = 0.5;
    double min_mole_fraction = 1.0e-14;
    double finite_difference_step = 1.0e-6;
    std::string phase = "liq";
};

struct ChemicalEquilibriumResultNative {
    bool success = false;
    std::string message;
    std::vector<double> composition;
    std::vector<double> activity_coefficients;
    std::vector<double> mass_balance_residuals;
    double charge_residual = 0.0;
    std::vector<double> reaction_residuals;
    std::map<std::string, double> diagnostics_double;
    std::map<std::string, int> diagnostics_int;
    std::map<std::string, bool> diagnostics_bool;
    std::map<std::string, std::string> diagnostics_string;
    std::map<std::string, std::vector<double>> diagnostics_vector;
};

ChemicalEquilibriumResultNative chemical_equilibrium_native(
    const std::shared_ptr<ePCSAFTMixtureNative>& mixture,
    double t,
    double p,
    const std::vector<double>& initial_x,
    const std::vector<double>& balance_matrix_row_major,
    int balance_rows,
    const std::vector<double>& total_vector,
    const std::vector<double>& reaction_stoichiometry_row_major,
    int reaction_rows,
    const std::vector<double>& log_equilibrium_constants,
    const ChemicalEquilibriumOptionsNative& options
);
