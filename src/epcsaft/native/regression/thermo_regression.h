#pragma once

#include "epcsaft_chemical_equilibrium.h"
#include "epcsaft_equilibrium.h"
#include "regression_types.h"

#include <memory>
#include <string>
#include <vector>

class ePCSAFTMixtureNative;

struct NativeThermoRegressionTarget {
    std::string family;
    std::string target;
    int index = -1;
    double observed = 0.0;
    double scale = 1.0;
};

struct NativeThermoRegressionRow {
    std::string row_id;
    std::string row_mode;
    double t = 0.0;
    double p = 1.0e5;
    std::vector<double> initial_x;
    std::vector<double> x_liq;
    std::vector<double> balance_matrix_row_major;
    int balance_rows = 0;
    std::vector<double> total_vector;
    std::vector<double> reaction_stoichiometry_row_major;
    int reaction_rows = 0;
    std::vector<double> log_equilibrium_constants;
    std::vector<int> reaction_standard_states;
    ChemicalEquilibriumOptionsNative speciation_options;
    ElectrolyteBubbleOptionsNative bubble_options;
    std::vector<std::string> targets_species;
    std::vector<std::string> vapor_species;
    std::vector<NativeThermoRegressionTarget> targets;
};

NativeRegressionResidualEvaluation evaluate_native_thermo_regression_rows(
    const std::shared_ptr<ePCSAFTMixtureNative>& mixture,
    const std::vector<std::string>& species,
    const std::vector<NativeThermoRegressionRow>& rows,
    double penalty_residual = 1.0e6
);

NativeRegressionFitResult fit_native_thermo_regression(
    const std::shared_ptr<ePCSAFTMixtureNative>& mixture,
    const std::vector<std::string>& species,
    const std::vector<NativeThermoRegressionRow>& rows,
    const std::vector<NativeRegressionParameterSpec>& parameters,
    const NativeRegressionFitOptions& options
);
