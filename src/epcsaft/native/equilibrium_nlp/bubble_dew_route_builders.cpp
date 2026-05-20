#include "route_builders.h"

#include "eos_phase_block.h"
#include "epcsaft_core_internal.h"
#include "nlp_problem.h"
#include "route_metadata.h"
#include "second_order.h"

#include <algorithm>
#include <cmath>
#include <numeric>
#include <utility>

namespace epcsaft::native::equilibrium_nlp {
namespace {

constexpr double kInitialPressure = 1.0e5;
constexpr double kInitialTemperature = 300.0;
constexpr double kGasConstant = 8.31446261815324;
constexpr double kSeparatedLiquidDensity = 8000.0;
constexpr double kMinimumLiquidVolume = 1.0e-6;
constexpr double kMaximumLiquidVolume = 5.0e-4;
constexpr double kMinimumVaporVolume = 1.0e-3;
constexpr double kMaximumVaporVolume = 1.0e6;
constexpr double kMinimumPhaseVolumeGap = 1.0e-7;

void apply_route_metadata(NeutralTwoPhaseEosNlpContract& out, const RouteMetadata& metadata) {
    out.variable_model = metadata.variable_model;
    out.density_backend = metadata.density_backend;
    out.residual_families = metadata.residual_families;
    out.constraint_families = metadata.constraint_families;
}

void apply_route_metadata(NeutralTwoPhaseEosPostsolve& out, const RouteMetadata& metadata) {
    out.density_backend = metadata.density_backend;
    out.residual_families = metadata.residual_families;
    out.constraint_families = metadata.constraint_families;
}

void apply_route_metadata(NeutralTwoPhaseEosRouteResult& out, const RouteMetadata& metadata) {
    out.variable_model = metadata.variable_model;
    out.density_backend = metadata.density_backend;
    out.residual_families = metadata.residual_families;
    out.constraint_families = metadata.constraint_families;
}

enum class DensitySeedMode {
    PhasePressureRoot,
    SeparatedPhaseRole,
};

void require_size(const std::vector<double>& values, std::size_t expected, const std::string& label) {
    if (values.size() == expected) {
        return;
    }
    throw ValueError(label + " size does not match the fixed-temperature pressure route.");
}

void require_positive_finite(double value, const std::string& label) {
    if (std::isfinite(value) && value > 0.0) {
        return;
    }
    throw ValueError(label + " must be positive and finite.");
}

std::vector<double> normalized_positive_values(const std::vector<double>& values, const std::string& label) {
    if (values.empty()) {
        throw ValueError(label + " requires at least one value.");
    }
    double total = 0.0;
    for (double value : values) {
        require_positive_finite(value, label + " value");
        total += value;
    }
    require_positive_finite(total, label + " total");
    std::vector<double> normalized;
    normalized.reserve(values.size());
    for (double value : values) {
        normalized.push_back(value / total);
    }
    return normalized;
}

std::vector<double> shifted_composition(const std::vector<double>& composition, double shift_sign = 1.0) {
    if (composition.empty()) {
        return {};
    }
    if (composition.size() == 1) {
        return {composition.front()};
    }
    const double triangular_sum = 0.5 * static_cast<double>(composition.size() * (composition.size() + 1));
    std::vector<double> shifted;
    shifted.reserve(composition.size());
    for (std::size_t index = 0; index < composition.size(); ++index) {
        const double triangular = static_cast<double>(index + 1) / triangular_sum;
        shifted.push_back(
            0.8 * composition[index]
            + 0.2 * (shift_sign > 0.0 ? triangular : (1.0 - triangular + 1.0 / triangular_sum))
        );
    }
    return normalized_positive_values(shifted, "shifted composition");
}

std::vector<double> charge_neutral_shifted_composition(
    const std::vector<double>& composition,
    const std::vector<double>& charges,
    const std::string& label,
    double shift_sign = 1.0
) {
    require_size(charges, composition.size(), label + " charge");
    if (composition.size() <= 1) {
        throw ValueError(label + " requires at least two species.");
    }
    double composition_charge = 0.0;
    double charge_square_weight = 0.0;
    for (std::size_t index = 0; index < composition.size(); ++index) {
        if (!std::isfinite(charges[index])) {
            throw ValueError(label + " charge values must be finite.");
        }
        composition_charge += composition[index] * charges[index];
        charge_square_weight += composition[index] * charges[index] * charges[index];
    }
    if (charge_square_weight <= 0.0) {
        throw ValueError(label + " requires at least one charged species.");
    }
    if (std::abs(composition_charge) > 1.0e-10) {
        throw ValueError(label + " fixed composition must be charge neutral.");
    }

    std::vector<double> positions;
    positions.reserve(composition.size());
    const double denominator = static_cast<double>(composition.size() - 1);
    for (std::size_t index = 0; index < composition.size(); ++index) {
        positions.push_back(-1.0 + 2.0 * static_cast<double>(index) / denominator);
    }
    double weighted_position = 0.0;
    for (std::size_t index = 0; index < composition.size(); ++index) {
        weighted_position += composition[index] * positions[index];
    }

    std::vector<double> direction;
    direction.reserve(composition.size());
    double charge_direction = 0.0;
    for (std::size_t index = 0; index < composition.size(); ++index) {
        const double value = positions[index] - weighted_position;
        direction.push_back(value);
        charge_direction += composition[index] * charges[index] * value;
    }
    const double charge_projection = charge_direction / charge_square_weight;
    double max_abs_direction = 0.0;
    for (std::size_t index = 0; index < direction.size(); ++index) {
        direction[index] -= charge_projection * charges[index];
        max_abs_direction = std::max(max_abs_direction, std::abs(direction[index]));
    }
    if (max_abs_direction <= 0.0) {
        throw ValueError(label + " could not construct a charge-neutral initial direction.");
    }

    std::vector<double> shifted;
    shifted.reserve(composition.size());
    double shifted_sum = 0.0;
    for (std::size_t index = 0; index < composition.size(); ++index) {
        const double value =
            composition[index] * (1.0 + 0.2 * shift_sign * direction[index] / max_abs_direction);
        require_positive_finite(value, label + " shifted composition");
        shifted.push_back(value);
        shifted_sum += value;
    }
    require_positive_finite(shifted_sum, label + " shifted composition sum");
    for (double& value : shifted) {
        value /= shifted_sum;
    }
    return shifted;
}

int phase_kind_from_index(int phase, const std::string& problem_name) {
    if (phase == 0) {
        return 0;
    }
    if (phase == 1) {
        return 1;
    }
    throw ValueError(problem_name + " phase index is out of range for density-root seeding.");
}

double phase_density_root_seed(
    const add_args& args,
    double temperature,
    double pressure,
    const std::vector<double>& composition,
    int phase,
    const std::string& problem_name
) {
    const DensitySolveResult density = density_solve_report_cpp(
        temperature,
        pressure,
        composition,
        phase_kind_from_index(phase, problem_name),
        args
    );
    if (!density.valid || !std::isfinite(density.rho) || density.rho <= 0.0) {
        throw ValueError(problem_name + " could not construct a phase-specific pressure-density root seed.");
    }
    return density.rho;
}

double separated_phase_role_density_seed(double temperature, double pressure, int phase, const std::string& problem_name) {
    if (phase == 0) {
        return kSeparatedLiquidDensity;
    }
    if (phase == 1) {
        return std::max(pressure / (kGasConstant * temperature), 1.0e-12);
    }
    throw ValueError(problem_name + " phase index is out of range for separated density seeding.");
}

bool fixed_temperature_pressure_seed_satisfies_volume_bounds(
    const std::vector<double>& variables,
    int species_count
) {
    const int local_variable_count = species_count + 1;
    if (variables.size() != static_cast<std::size_t>(2 * local_variable_count + 1)) {
        return false;
    }
    const double liquid_volume = variables[static_cast<std::size_t>(species_count)];
    const double vapor_volume = variables[static_cast<std::size_t>(local_variable_count + species_count)];
    return liquid_volume >= kMinimumLiquidVolume
        && liquid_volume <= kMaximumLiquidVolume
        && vapor_volume >= kMinimumVaporVolume
        && vapor_volume <= kMaximumVaporVolume
        && vapor_volume - liquid_volume >= kMinimumPhaseVolumeGap;
}

struct NamedInitialVariables {
    std::string seed_name;
    std::vector<double> variables;
};

std::vector<double> build_pressure_route_initial_variables(
    const add_args& args,
    const std::vector<double>& fixed_composition,
    int fixed_phase_index,
    double temperature,
    const std::vector<double>& charges,
    const std::string& problem_name,
    double shift_sign,
    DensitySeedMode density_seed_mode
) {
    const std::vector<double> shifted = charges.empty()
        ? shifted_composition(fixed_composition, shift_sign)
        : charge_neutral_shifted_composition(
            fixed_composition,
            charges,
            problem_name + " shifted composition",
            shift_sign
        );
    std::vector<double> out;
    out.reserve(2 * (fixed_composition.size() + 1) + 1);
    for (int phase = 0; phase < 2; ++phase) {
        const std::vector<double>& composition = phase == fixed_phase_index ? fixed_composition : shifted;
        out.insert(out.end(), composition.begin(), composition.end());
        const double density = density_seed_mode == DensitySeedMode::PhasePressureRoot
            ? phase_density_root_seed(args, temperature, kInitialPressure, composition, phase, problem_name)
            : separated_phase_role_density_seed(temperature, kInitialPressure, phase, problem_name);
        out.push_back(1.0 / density);
    }
    out.push_back(kInitialPressure);
    return out;
}

std::vector<double> build_temperature_route_initial_variables(
    const add_args& args,
    const std::vector<double>& fixed_composition,
    int fixed_phase_index,
    double target_pressure,
    const std::string& problem_name,
    double shift_sign,
    DensitySeedMode density_seed_mode
) {
    const std::vector<double> shifted = shifted_composition(fixed_composition, shift_sign);
    std::vector<double> out;
    out.reserve(2 * (fixed_composition.size() + 1) + 1);
    for (int phase = 0; phase < 2; ++phase) {
        const std::vector<double>& composition = phase == fixed_phase_index ? fixed_composition : shifted;
        out.insert(out.end(), composition.begin(), composition.end());
        const double density = density_seed_mode == DensitySeedMode::PhasePressureRoot
            ? phase_density_root_seed(args, kInitialTemperature, target_pressure, composition, phase, problem_name)
            : separated_phase_role_density_seed(kInitialTemperature, target_pressure, phase, problem_name);
        out.push_back(1.0 / density);
    }
    out.push_back(kInitialTemperature);
    return out;
}

std::vector<NamedInitialVariables> pressure_route_seed_candidates(
    const add_args& args,
    const std::vector<double>& fixed_composition,
    int fixed_phase_index,
    double temperature,
    const std::vector<double>& charges,
    const std::string& problem_name
) {
    std::vector<NamedInitialVariables> out;
    auto append_if_feasible = [&](std::string seed_name, std::vector<double> variables) {
        if (fixed_temperature_pressure_seed_satisfies_volume_bounds(
                variables,
                static_cast<int>(fixed_composition.size())
            )) {
            out.push_back({std::move(seed_name), std::move(variables)});
        }
    };
    if (charges.empty()) {
        append_if_feasible(
            "canonical_phase_density_root",
            build_pressure_route_initial_variables(
                args,
                fixed_composition,
                fixed_phase_index,
                temperature,
                charges,
                problem_name,
                1.0,
                DensitySeedMode::PhasePressureRoot
            )
        );
        append_if_feasible(
            "mirrored_phase_density_root",
            build_pressure_route_initial_variables(
                args,
                fixed_composition,
                fixed_phase_index,
                temperature,
                charges,
                problem_name,
                -1.0,
                DensitySeedMode::PhasePressureRoot
            )
        );
    }
    out.push_back({
        "canonical_shifted_partner_phase",
        build_pressure_route_initial_variables(
            args,
            fixed_composition,
            fixed_phase_index,
            temperature,
            charges,
            problem_name,
            1.0,
            DensitySeedMode::SeparatedPhaseRole
        )
    });
    out.push_back({
        "mirrored_shifted_partner_phase",
        build_pressure_route_initial_variables(
            args,
            fixed_composition,
            fixed_phase_index,
            temperature,
            charges,
            problem_name,
            -1.0,
            DensitySeedMode::SeparatedPhaseRole
        )
    });
    return out;
}

std::vector<NamedInitialVariables> temperature_route_seed_candidates(
    const add_args& args,
    const std::vector<double>& fixed_composition,
    int fixed_phase_index,
    double target_pressure,
    const std::string& problem_name
) {
    std::vector<NamedInitialVariables> out;
    auto append_if_feasible = [&](std::string seed_name, std::vector<double> variables) {
        if (fixed_temperature_pressure_seed_satisfies_volume_bounds(
                variables,
                static_cast<int>(fixed_composition.size())
            )) {
            out.push_back({std::move(seed_name), std::move(variables)});
        }
    };
    append_if_feasible(
        "canonical_phase_density_root",
        build_temperature_route_initial_variables(
            args,
            fixed_composition,
            fixed_phase_index,
            target_pressure,
            problem_name,
            1.0,
            DensitySeedMode::PhasePressureRoot
        )
    );
    append_if_feasible(
        "mirrored_phase_density_root",
        build_temperature_route_initial_variables(
            args,
            fixed_composition,
            fixed_phase_index,
            target_pressure,
            problem_name,
            -1.0,
            DensitySeedMode::PhasePressureRoot
        )
    );
    out.push_back({
        "canonical_shifted_partner_phase",
        build_temperature_route_initial_variables(
            args,
            fixed_composition,
            fixed_phase_index,
            target_pressure,
            problem_name,
            1.0,
            DensitySeedMode::SeparatedPhaseRole
        )
    });
    out.push_back({
        "mirrored_shifted_partner_phase",
        build_temperature_route_initial_variables(
            args,
            fixed_composition,
            fixed_phase_index,
            target_pressure,
            problem_name,
            -1.0,
            DensitySeedMode::SeparatedPhaseRole
        )
    });
    return out;
}

int neutral_route_quality(const NeutralTwoPhaseEosRouteResult& result) {
    if (result.accepted) {
        return 3;
    }
    if (result.solver_accepted) {
        return 2;
    }
    if (result.ran) {
        return 1;
    }
    return 0;
}

RouteSeedAttempt neutral_seed_attempt_from_result(const NeutralTwoPhaseEosRouteResult& result) {
    RouteSeedAttempt out;
    out.seed_name = result.seed_name;
    out.status = result.status;
    out.solver_status = result.solver_status;
    out.application_status = result.application_status;
    out.solver_accepted = result.solver_accepted;
    out.accepted = result.accepted;
    out.iteration_count = result.iteration_count;
    out.objective = result.objective;
    out.phase_distance = result.postsolve.phase_distance;
    out.material_balance_norm = result.postsolve.material_balance_norm;
    out.charge_balance_norm = result.postsolve.charge_balance_norm;
    out.pressure_consistency_norm = result.postsolve.pressure_consistency_norm;
    out.chemical_potential_consistency_norm = result.postsolve.chemical_potential_consistency_norm;
    return out;
}

std::vector<std::vector<double>> pressure_route_phase_amounts(
    const std::vector<double>& variables,
    int species_count
) {
    const int local_variable_count = species_count + 1;
    require_size(
        variables,
        static_cast<std::size_t>(2 * local_variable_count + 1),
        "fixed-temperature pressure route variable"
    );
    std::vector<std::vector<double>> phase_amounts(2, std::vector<double>(static_cast<std::size_t>(species_count)));
    for (int phase = 0; phase < 2; ++phase) {
        const int offset = phase * local_variable_count;
        for (int species = 0; species < species_count; ++species) {
            phase_amounts[static_cast<std::size_t>(phase)][static_cast<std::size_t>(species)] =
                variables[static_cast<std::size_t>(offset + species)];
        }
    }
    return phase_amounts;
}

std::vector<double> pressure_route_phase_volumes(const std::vector<double>& variables, int species_count) {
    const int local_variable_count = species_count + 1;
    require_size(
        variables,
        static_cast<std::size_t>(2 * local_variable_count + 1),
        "fixed-temperature pressure route variable"
    );
    return {
        variables[static_cast<std::size_t>(species_count)],
        variables[static_cast<std::size_t>(local_variable_count + species_count)]
    };
}

std::vector<double> summed_feed_amounts(const std::vector<std::vector<double>>& phase_amounts, int species_count) {
    std::vector<double> feed_amounts(static_cast<std::size_t>(species_count), 0.0);
    for (const auto& phase : phase_amounts) {
        require_size(phase, static_cast<std::size_t>(species_count), "fixed-temperature pressure phase amount");
        for (int species = 0; species < species_count; ++species) {
            feed_amounts[static_cast<std::size_t>(species)] += phase[static_cast<std::size_t>(species)];
        }
    }
    return feed_amounts;
}

double fixed_composition_norm(
    const std::vector<std::vector<double>>& phase_amounts,
    int fixed_phase_index,
    const std::vector<double>& fixed_composition
) {
    const auto& fixed_amounts = phase_amounts[static_cast<std::size_t>(fixed_phase_index)];
    const double total = std::accumulate(fixed_amounts.begin(), fixed_amounts.end(), 0.0);
    require_positive_finite(total, "fixed-temperature pressure fixed-phase amount total");
    double norm = 0.0;
    for (std::size_t species = 0; species < fixed_composition.size(); ++species) {
        norm = std::max(norm, std::abs(fixed_amounts[species] / total - fixed_composition[species]));
    }
    return norm;
}

double phase_total_norm(const std::vector<std::vector<double>>& phase_amounts) {
    double norm = 0.0;
    for (const auto& phase : phase_amounts) {
        const double total = std::accumulate(phase.begin(), phase.end(), 0.0);
        norm = std::max(norm, std::abs(total - 1.0));
    }
    return norm;
}

double phase_charge_norm(
    const std::vector<std::vector<double>>& phase_amounts,
    const std::vector<double>& charges
) {
    if (charges.empty()) {
        return 0.0;
    }
    double norm = 0.0;
    for (const auto& phase : phase_amounts) {
        require_size(phase, charges.size(), "fixed-temperature pressure charge-balance phase amount");
        double phase_charge = 0.0;
        for (std::size_t species = 0; species < charges.size(); ++species) {
            phase_charge += phase[species] * charges[species];
        }
        norm = std::max(norm, std::abs(phase_charge));
    }
    return norm;
}

bool route_has_active_association_sites(const add_args& args) {
    for (int sites : args.assoc_num) {
        if (sites > 0) {
            return true;
        }
    }
    return false;
}

bool route_supports_exact_phase_derivatives(const add_args& args) {
    return !route_has_active_association_sites(args) && (args.z.empty() || args.born_model <= 1);
}

void require_square_block(
    const std::vector<double>& values,
    int variable_count,
    const std::string& label
) {
    require_size(
        values,
        static_cast<std::size_t>(variable_count) * static_cast<std::size_t>(variable_count),
        label
    );
}

void require_third_derivative_tensor(
    const std::vector<double>& values,
    int variable_count,
    const std::string& label
) {
    require_size(
        values,
        static_cast<std::size_t>(variable_count)
            * static_cast<std::size_t>(variable_count)
            * static_cast<std::size_t>(variable_count),
        label
    );
}

void symmetrize_square_block(std::vector<double>& values, int variable_count) {
    const std::size_t n = static_cast<std::size_t>(variable_count);
    for (std::size_t row = 0; row < n; ++row) {
        for (std::size_t col = 0; col < row; ++col) {
            const double value = 0.5 * (values[row * n + col] + values[col * n + row]);
            values[row * n + col] = value;
            values[col * n + row] = value;
        }
    }
}

std::vector<int> fixed_temperature_global_indices(int phase, int local_variable_count) {
    std::vector<int> indices;
    indices.reserve(static_cast<std::size_t>(local_variable_count));
    const int offset = phase * local_variable_count;
    for (int local = 0; local < local_variable_count; ++local) {
        indices.push_back(offset + local);
    }
    return indices;
}

std::vector<int> fixed_pressure_global_indices(
    int phase,
    int local_variable_count,
    int temperature_col
) {
    std::vector<int> indices = fixed_temperature_global_indices(phase, local_variable_count);
    indices.push_back(temperature_col);
    return indices;
}

void add_local_hessian_to_dense(
    std::vector<double>& dense,
    int variable_count,
    const std::vector<int>& global_indices,
    const std::vector<double>& local_hessian,
    double scale,
    const std::string& label
) {
    const int local_variable_count = static_cast<int>(global_indices.size());
    require_square_block(local_hessian, local_variable_count, label);
    const std::size_t n = static_cast<std::size_t>(variable_count);
    require_square_block(dense, variable_count, "dense route Hessian");
    for (int local_row = 0; local_row < local_variable_count; ++local_row) {
        const int global_row = global_indices[static_cast<std::size_t>(local_row)];
        for (int local_col = 0; local_col < local_variable_count; ++local_col) {
            const int global_col = global_indices[static_cast<std::size_t>(local_col)];
            dense[static_cast<std::size_t>(global_row) * n + static_cast<std::size_t>(global_col)] +=
                scale * local_hessian[
                    static_cast<std::size_t>(local_row * local_variable_count + local_col)
                ];
        }
    }
}

void add_local_hessian_to_constraint_tensor(
    std::vector<double>& tensor,
    int constraint_row,
    int variable_count,
    const std::vector<int>& global_indices,
    const std::vector<double>& local_hessian,
    double scale,
    const std::string& label
) {
    const int local_variable_count = static_cast<int>(global_indices.size());
    require_square_block(local_hessian, local_variable_count, label);
    const std::size_t n = static_cast<std::size_t>(variable_count);
    const std::size_t row_offset = static_cast<std::size_t>(constraint_row) * n * n;
    if (tensor.size() < row_offset + n * n) {
        throw ValueError(label + " target tensor is smaller than the requested constraint row.");
    }
    for (int local_row = 0; local_row < local_variable_count; ++local_row) {
        const int global_row = global_indices[static_cast<std::size_t>(local_row)];
        for (int local_col = 0; local_col < local_variable_count; ++local_col) {
            const int global_col = global_indices[static_cast<std::size_t>(local_col)];
            tensor[row_offset + static_cast<std::size_t>(global_row) * n + static_cast<std::size_t>(global_col)] +=
                scale * local_hessian[
                    static_cast<std::size_t>(local_row * local_variable_count + local_col)
                ];
        }
    }
}

std::vector<double> third_derivative_slice(
    const std::vector<double>& tensor,
    int component,
    int variable_count,
    const std::string& label
) {
    require_third_derivative_tensor(tensor, variable_count, label);
    const std::size_t n = static_cast<std::size_t>(variable_count);
    const std::size_t offset = static_cast<std::size_t>(component) * n * n;
    return std::vector<double>(
        tensor.begin() + static_cast<std::ptrdiff_t>(offset),
        tensor.begin() + static_cast<std::ptrdiff_t>(offset + n * n)
    );
}

struct TemperatureRoutePhaseBlock {
    EosPhaseBlockResult block;
    std::vector<double> gradient;
    std::vector<double> objective_hessian_row_major;
    std::vector<double> objective_third_derivative_tensor_row_major;
    std::vector<double> pressure_jacobian_row_major;
    std::vector<double> pressure_hessian_row_major;
};

TemperatureRoutePhaseBlock evaluate_temperature_route_phase_block(
    const add_args& args,
    double temperature,
    double target_pressure,
    const std::vector<double>& amounts,
    double volume
) {
    TemperatureRoutePhaseBlock out;
    out.block = evaluate_eos_phase_block(args, temperature, target_pressure, amounts, volume);
    double cppad_objective = 0.0;
    eos_phase_temperature_variable_derivatives_cpp(
        temperature,
        target_pressure,
        amounts,
        volume,
        args,
        &cppad_objective,
        &out.gradient,
        &out.objective_hessian_row_major,
        &out.objective_third_derivative_tensor_row_major
    );
    const int local_variable_count = static_cast<int>(amounts.size()) + 2;
    if (out.gradient.size() != static_cast<std::size_t>(local_variable_count)
        || out.objective_hessian_row_major.size()
            != static_cast<std::size_t>(local_variable_count * local_variable_count)) {
        throw ValueError("fixed-pressure temperature route CppAD derivative shape did not match variables.");
    }
    require_third_derivative_tensor(
        out.objective_third_derivative_tensor_row_major,
        local_variable_count,
        "fixed-pressure temperature route CppAD third-derivative tensor"
    );
    const double objective_scale = std::max(1.0, std::abs(out.block.objective));
    if (std::abs(cppad_objective - out.block.objective) > 1.0e-8 * objective_scale) {
        throw ValueError("fixed-pressure temperature route CppAD value did not match the EOS block value.");
    }

    const int volume_row = local_variable_count - 2;
    const int temperature_col = local_variable_count - 1;
    const double gas_constant = kb * N_AV;
    out.pressure_jacobian_row_major.reserve(static_cast<std::size_t>(local_variable_count));
    for (int col = 0; col < local_variable_count; ++col) {
        double value = -out.block.gas_constant_temperature
            * out.objective_hessian_row_major[static_cast<std::size_t>(volume_row * local_variable_count + col)];
        if (col == temperature_col) {
            value -= gas_constant * out.gradient[static_cast<std::size_t>(volume_row)];
        }
        out.pressure_jacobian_row_major.push_back(value);
    }
    out.pressure_hessian_row_major.assign(
        static_cast<std::size_t>(local_variable_count * local_variable_count),
        0.0
    );
    for (int row = 0; row < local_variable_count; ++row) {
        for (int col = 0; col < local_variable_count; ++col) {
            double value = -out.block.gas_constant_temperature
                * out.objective_third_derivative_tensor_row_major[
                    static_cast<std::size_t>(
                        volume_row * local_variable_count * local_variable_count
                        + row * local_variable_count
                        + col
                    )
                ];
            if (row == temperature_col) {
                value -= gas_constant * out.objective_hessian_row_major[
                    static_cast<std::size_t>(volume_row * local_variable_count + col)
                ];
            }
            if (col == temperature_col) {
                value -= gas_constant * out.objective_hessian_row_major[
                    static_cast<std::size_t>(volume_row * local_variable_count + row)
                ];
            }
            out.pressure_hessian_row_major[
                static_cast<std::size_t>(row * local_variable_count + col)
            ] = value;
        }
    }
    symmetrize_square_block(out.pressure_hessian_row_major, local_variable_count);
    return out;
}

class NeutralFixedTemperaturePressureProblem final : public NlpProblem {
public:
    NeutralFixedTemperaturePressureProblem(
        add_args args,
        double temperature,
        std::vector<double> fixed_composition,
        int fixed_phase_index,
        std::string problem_name,
        std::vector<double> charges = {},
        double charge_constraint_tolerance = 0.0
    )
        : args_(std::move(args)),
          temperature_(temperature),
          fixed_composition_(normalized_positive_values(fixed_composition, problem_name + " composition")),
          fixed_phase_index_(fixed_phase_index),
          problem_name_(std::move(problem_name)),
          charges_(std::move(charges)),
          charge_constraint_tolerance_(charge_constraint_tolerance) {
        require_positive_finite(temperature_, problem_name_ + " temperature");
        if (!std::isfinite(charge_constraint_tolerance_) || charge_constraint_tolerance_ < 0.0) {
            throw ValueError(problem_name_ + " charge constraint tolerance must be finite and non-negative.");
        }
        if (fixed_phase_index_ < 0 || fixed_phase_index_ >= phase_count()) {
            throw ValueError(problem_name_ + " fixed phase index is out of range.");
        }
        species_count_ = static_cast<int>(fixed_composition_.size());
        if (!charges_.empty()) {
            charge_neutral_shifted_composition(
                fixed_composition_,
                charges_,
                problem_name_ + " fixed composition"
            );
        }
    }

    std::string name() const override {
        return problem_name_;
    }

    int variable_count() const override {
        return phase_count() * local_variable_count() + 1;
    }

    int constraint_count() const override {
        return composition_constraint_count() + 2 * phase_count() + species_count_ + 1
            + charge_constraint_count();
    }

    int jacobian_nonzero_count() const override {
        const int composition_nonzeros = composition_constraint_count() * species_count_;
        const int phase_total_nonzeros = phase_count() * species_count_;
        int charge_nonzeros = 0;
        for (double charge : charges_) {
            if (charge != 0.0) {
                charge_nonzeros += phase_count();
            }
        }
        const int pressure_nonzeros = phase_count() * (local_variable_count() + 1);
        const int chemical_nonzeros = species_count_ * phase_count() * local_variable_count();
        const int gap_nonzeros = 2;
        return composition_nonzeros + phase_total_nonzeros + charge_nonzeros + pressure_nonzeros + chemical_nonzeros
            + gap_nonzeros;
    }

    NlpBounds bounds() const override {
        NlpBounds out;
        out.variable_lower.reserve(static_cast<std::size_t>(variable_count()));
        out.variable_upper.reserve(static_cast<std::size_t>(variable_count()));
        for (int phase = 0; phase < phase_count(); ++phase) {
            for (int species = 0; species < species_count_; ++species) {
                out.variable_lower.push_back(1.0e-14);
                out.variable_upper.push_back(10.0);
            }
            if (phase == liquid_phase_index()) {
                out.variable_lower.push_back(minimum_liquid_volume_);
                out.variable_upper.push_back(maximum_liquid_volume_);
            } else if (phase == vapor_phase_index()) {
                out.variable_lower.push_back(minimum_vapor_volume_);
                out.variable_upper.push_back(maximum_vapor_volume_);
            } else {
                throw ValueError(problem_name_ + " phase role is out of range.");
            }
        }
        out.variable_lower.push_back(1.0);
        out.variable_upper.push_back(1.0e9);
        out.constraint_lower.assign(static_cast<std::size_t>(constraint_count()), 0.0);
        out.constraint_upper.assign(static_cast<std::size_t>(constraint_count()), 0.0);
        if (charge_constraint_count() > 0 && charge_constraint_tolerance_ > 0.0) {
            const int charge_row_start = composition_constraint_count() + phase_count();
            for (int phase = 0; phase < charge_constraint_count(); ++phase) {
                out.constraint_lower[static_cast<std::size_t>(charge_row_start + phase)] =
                    -charge_constraint_tolerance_;
                out.constraint_upper[static_cast<std::size_t>(charge_row_start + phase)] =
                    charge_constraint_tolerance_;
            }
        }
        out.constraint_lower.back() = minimum_phase_volume_gap_;
        out.constraint_upper.back() = 1.0e12;
        return out;
    }

    std::vector<double> initial_point() const override {
        if (charges_.empty()) {
            const std::vector<double> root = build_pressure_route_initial_variables(
                args_,
                fixed_composition_,
                fixed_phase_index_,
                temperature_,
                charges_,
                problem_name_,
                1.0,
                DensitySeedMode::PhasePressureRoot
            );
            if (fixed_temperature_pressure_seed_satisfies_volume_bounds(root, species_count_)) {
                return root;
            }
        }
        return build_pressure_route_initial_variables(
            args_,
            fixed_composition_,
            fixed_phase_index_,
            temperature_,
            charges_,
            problem_name_,
            1.0,
            DensitySeedMode::SeparatedPhaseRole
        );
    }

    double objective(const std::vector<double>& variables) const override {
        double value = 0.0;
        for (const EosPhaseBlockResult& block : phase_blocks(variables)) {
            value += block.objective;
        }
        return value;
    }

    std::vector<double> objective_gradient(const std::vector<double>& variables) const override {
        std::vector<double> out;
        out.reserve(static_cast<std::size_t>(variable_count()));
        double pressure_derivative = 0.0;
        for (const EosPhaseBlockResult& block : phase_blocks(variables)) {
            out.insert(out.end(), block.gradient.begin(), block.gradient.end());
            pressure_derivative += block.volume / block.gas_constant_temperature;
        }
        out.push_back(pressure_derivative);
        return out;
    }

    std::vector<double> constraints(const std::vector<double>& variables) const override {
        const auto amounts = pressure_route_phase_amounts(variables, species_count_);
        std::vector<double> out(static_cast<std::size_t>(constraint_count()), 0.0);
        int row = 0;
        const auto& fixed_amounts = amounts[static_cast<std::size_t>(fixed_phase_index_)];
        const double fixed_total = std::accumulate(fixed_amounts.begin(), fixed_amounts.end(), 0.0);
        for (int species = 0; species < composition_constraint_count(); ++species) {
            out[static_cast<std::size_t>(row++)] =
                fixed_amounts[static_cast<std::size_t>(species)]
                - fixed_composition_[static_cast<std::size_t>(species)] * fixed_total;
        }
        for (int phase = 0; phase < phase_count(); ++phase) {
            const auto& phase_amounts = amounts[static_cast<std::size_t>(phase)];
            out[static_cast<std::size_t>(row++)] =
                std::accumulate(phase_amounts.begin(), phase_amounts.end(), 0.0) - 1.0;
        }
        for (int phase = 0; phase < charge_constraint_count(); ++phase) {
            const auto& phase_amounts = amounts[static_cast<std::size_t>(phase)];
            double phase_charge = 0.0;
            for (int species = 0; species < species_count_; ++species) {
                phase_charge += phase_amounts[static_cast<std::size_t>(species)]
                    * charges_[static_cast<std::size_t>(species)];
            }
            out[static_cast<std::size_t>(row++)] = phase_charge;
        }
        const auto blocks = phase_blocks(variables);
        for (const EosPhaseBlockResult& block : blocks) {
            out[static_cast<std::size_t>(row++)] = block.pressure_consistency_residual;
        }
        for (int species = 0; species < species_count_; ++species) {
            out[static_cast<std::size_t>(row++)] =
                blocks[0].gradient[static_cast<std::size_t>(species)]
                - blocks[1].gradient[static_cast<std::size_t>(species)];
        }
        out[static_cast<std::size_t>(row++)] =
            variables[static_cast<std::size_t>(vapor_volume_col())]
            - variables[static_cast<std::size_t>(liquid_volume_col())];
        return out;
    }

    NlpJacobianStructure jacobian_structure() const override {
        NlpJacobianStructure out;
        out.rows.reserve(static_cast<std::size_t>(jacobian_nonzero_count()));
        out.cols.reserve(static_cast<std::size_t>(jacobian_nonzero_count()));
        int row = 0;
        for (int species = 0; species < composition_constraint_count(); ++species) {
            for (int col = 0; col < species_count_; ++col) {
                out.rows.push_back(row);
                out.cols.push_back(fixed_col(col));
            }
            ++row;
        }
        for (int phase = 0; phase < phase_count(); ++phase) {
            const int offset = phase * local_variable_count();
            for (int species = 0; species < species_count_; ++species) {
                out.rows.push_back(row);
                out.cols.push_back(offset + species);
            }
            ++row;
        }
        for (int phase = 0; phase < charge_constraint_count(); ++phase) {
            const int offset = phase * local_variable_count();
            for (int species = 0; species < species_count_; ++species) {
                if (charges_[static_cast<std::size_t>(species)] == 0.0) {
                    continue;
                }
                out.rows.push_back(row);
                out.cols.push_back(offset + species);
            }
            ++row;
        }
        for (int phase = 0; phase < phase_count(); ++phase) {
            const int offset = phase * local_variable_count();
            for (int col = 0; col < local_variable_count(); ++col) {
                out.rows.push_back(row);
                out.cols.push_back(offset + col);
            }
            out.rows.push_back(row);
            out.cols.push_back(variable_count() - 1);
            ++row;
        }
        for (int species = 0; species < species_count_; ++species) {
            for (int phase = 0; phase < phase_count(); ++phase) {
                const int offset = phase * local_variable_count();
                for (int col = 0; col < local_variable_count(); ++col) {
                    out.rows.push_back(row);
                    out.cols.push_back(offset + col);
                }
            }
            ++row;
        }
        out.rows.push_back(row);
        out.cols.push_back(liquid_volume_col());
        out.rows.push_back(row);
        out.cols.push_back(vapor_volume_col());
        return out;
    }

    std::vector<double> jacobian_values(const std::vector<double>& variables) const override {
        const auto blocks = phase_blocks(variables);
        std::vector<double> out;
        out.reserve(static_cast<std::size_t>(jacobian_nonzero_count()));
        int row = 0;
        for (int species = 0; species < composition_constraint_count(); ++species) {
            for (int col = 0; col < species_count_; ++col) {
                out.push_back((col == species ? 1.0 : 0.0) - fixed_composition_[static_cast<std::size_t>(species)]);
            }
            ++row;
        }
        for (int phase = 0; phase < phase_count(); ++phase) {
            for (int species = 0; species < species_count_; ++species) {
                out.push_back(1.0);
            }
            ++row;
        }
        for (int phase = 0; phase < charge_constraint_count(); ++phase) {
            for (int species = 0; species < species_count_; ++species) {
                const double charge = charges_[static_cast<std::size_t>(species)];
                if (charge == 0.0) {
                    continue;
                }
                out.push_back(charge);
            }
            ++row;
        }
        for (int phase = 0; phase < phase_count(); ++phase) {
            const EosPhaseBlockResult& block = blocks[static_cast<std::size_t>(phase)];
            if (block.pressure_jacobian_row_major.size() != static_cast<std::size_t>(local_variable_count())) {
                throw ValueError(problem_name_ + " pressure Jacobian size did not match variables.");
            }
            for (int col = 0; col < local_variable_count(); ++col) {
                out.push_back(block.pressure_jacobian_row_major[static_cast<std::size_t>(col)]);
            }
            out.push_back(-1.0);
            ++row;
        }
        for (int species = 0; species < species_count_; ++species) {
            for (int phase = 0; phase < phase_count(); ++phase) {
                const EosPhaseBlockResult& block = blocks[static_cast<std::size_t>(phase)];
                if (block.objective_curvature_rows != local_variable_count()
                    || block.objective_curvature_cols != local_variable_count()
                    || block.objective_curvature_row_major.size()
                        != static_cast<std::size_t>(local_variable_count() * local_variable_count())) {
                    throw ValueError(problem_name_ + " chemical-potential Jacobian size did not match variables.");
                }
                for (int col = 0; col < local_variable_count(); ++col) {
                    out.push_back((phase == 0 ? 1.0 : -1.0) * block.objective_curvature_row_major[
                            static_cast<std::size_t>(species * local_variable_count() + col)
                        ]);
                }
            }
            ++row;
        }
        out.push_back(-1.0);
        out.push_back(1.0);
        return out;
    }

    bool has_exact_hessian() const override {
        return route_supports_exact_phase_derivatives(args_);
    }

    int hessian_nonzero_count() const override {
        return LagrangianHessianAssembler(variable_count()).nonzero_count();
    }

    NlpHessianStructure hessian_structure() const override {
        return LagrangianHessianAssembler(variable_count()).structure();
    }

    std::vector<double> hessian_values(
        const std::vector<double>& variables,
        double objective_factor,
        const std::vector<double>& constraint_multipliers
    ) const override {
        if (!has_exact_hessian()) {
            throw ValueError(problem_name_ + " exact Hessian requires direct non-associating phase derivatives.");
        }
        if (constraint_multipliers.size() != static_cast<std::size_t>(constraint_count())) {
            throw ValueError(problem_name_ + " Hessian multiplier vector size does not match the constraint count.");
        }

        const auto blocks = phase_blocks(variables);
        const int n = variable_count();
        ObjectiveSecondOrderData objective;
        objective.variable_count = n;
        objective.value = 0.0;
        objective.gradient.assign(static_cast<std::size_t>(n), 0.0);
        objective.hessian_row_major.assign(static_cast<std::size_t>(n * n), 0.0);
        objective.backend = "cppad_phase_pressure_route";
        const int pressure_col = n - 1;
        for (int phase = 0; phase < phase_count(); ++phase) {
            const EosPhaseBlockResult& block = blocks[static_cast<std::size_t>(phase)];
            require_square_block(
                block.objective_curvature_row_major,
                local_variable_count(),
                problem_name_ + " objective Hessian block"
            );
            require_third_derivative_tensor(
                block.objective_third_derivative_tensor_row_major,
                local_variable_count(),
                problem_name_ + " objective third-derivative tensor"
            );
            const std::vector<int> phase_indices =
                fixed_temperature_global_indices(phase, local_variable_count());
            objective.value += block.objective;
            for (int local = 0; local < local_variable_count(); ++local) {
                objective.gradient[static_cast<std::size_t>(phase_indices[static_cast<std::size_t>(local)])] =
                    block.gradient[static_cast<std::size_t>(local)];
            }
            objective.gradient[static_cast<std::size_t>(pressure_col)] +=
                block.volume / block.gas_constant_temperature;
            add_local_hessian_to_dense(
                objective.hessian_row_major,
                n,
                phase_indices,
                block.objective_curvature_row_major,
                1.0,
                problem_name_ + " objective Hessian block"
            );
            const int phase_volume_col = phase * local_variable_count() + species_count_;
            const double pressure_volume_cross = 1.0 / block.gas_constant_temperature;
            objective.hessian_row_major[
                static_cast<std::size_t>(phase_volume_col * n + pressure_col)
            ] += pressure_volume_cross;
            objective.hessian_row_major[
                static_cast<std::size_t>(pressure_col * n + phase_volume_col)
            ] += pressure_volume_cross;
        }

        ConstraintSecondOrderData constraints_data;
        constraints_data.constraint_count = constraint_count();
        constraints_data.variable_count = n;
        constraints_data.values = this->constraints(variables);
        constraints_data.hessian_tensor_row_major.assign(
            static_cast<std::size_t>(constraint_count() * n * n),
            0.0
        );
        constraints_data.has_hessian.assign(static_cast<std::size_t>(constraint_count()), false);
        constraints_data.backend = "cppad_phase_pressure_route";

        const int pressure_row_start =
            composition_constraint_count() + phase_count() + charge_constraint_count();
        for (int phase = 0; phase < phase_count(); ++phase) {
            const EosPhaseBlockResult& block = blocks[static_cast<std::size_t>(phase)];
            require_square_block(
                block.pressure_hessian_row_major,
                local_variable_count(),
                problem_name_ + " pressure Hessian block"
            );
            const int constraint_row = pressure_row_start + phase;
            constraints_data.has_hessian[static_cast<std::size_t>(constraint_row)] = true;
            add_local_hessian_to_constraint_tensor(
                constraints_data.hessian_tensor_row_major,
                constraint_row,
                n,
                fixed_temperature_global_indices(phase, local_variable_count()),
                block.pressure_hessian_row_major,
                1.0,
                problem_name_ + " pressure Hessian block"
            );
        }

        const int chemical_row_start = pressure_row_start + phase_count();
        for (int species = 0; species < species_count_; ++species) {
            const int constraint_row = chemical_row_start + species;
            constraints_data.has_hessian[static_cast<std::size_t>(constraint_row)] = true;
            for (int phase = 0; phase < phase_count(); ++phase) {
                const EosPhaseBlockResult& block = blocks[static_cast<std::size_t>(phase)];
                const std::vector<double> chemical_hessian = third_derivative_slice(
                    block.objective_third_derivative_tensor_row_major,
                    species,
                    local_variable_count(),
                    problem_name_ + " chemical-potential Hessian block"
                );
                add_local_hessian_to_constraint_tensor(
                    constraints_data.hessian_tensor_row_major,
                    constraint_row,
                    n,
                    fixed_temperature_global_indices(phase, local_variable_count()),
                    chemical_hessian,
                    phase == 0 ? 1.0 : -1.0,
                    problem_name_ + " chemical-potential Hessian block"
                );
            }
        }

        return LagrangianHessianAssembler(n).values(
            objective_factor,
            objective,
            constraints_data,
            constraint_multipliers
        );
    }

    std::string hessian_backend() const override {
        return "cppad_phase_pressure_route";
    }

    NlpScaling scaling() const override {
        NlpScaling out;
        out.objective = 1.0;
        out.variables.assign(static_cast<std::size_t>(variable_count()), 1.0);
        out.variables.back() = 1.0e-5;
        out.constraints.assign(static_cast<std::size_t>(constraint_count()), 1.0);
        return out;
    }

    std::map<std::string, std::string> diagnostics() const override {
        return route_metadata_diagnostics(fixed_temperature_pressure_route_metadata(!charges_.empty()));
    }

    int species_count() const {
        return species_count_;
    }

    int phase_count() const {
        return 2;
    }

private:
    int local_variable_count() const {
        return species_count_ + 1;
    }

    int composition_constraint_count() const {
        return std::max(0, species_count_ - 1);
    }

    int charge_constraint_count() const {
        return charges_.empty() ? 0 : phase_count();
    }

    int fixed_col(int species) const {
        return fixed_phase_index_ * local_variable_count() + species;
    }

    int liquid_phase_index() const {
        return 0;
    }

    int vapor_phase_index() const {
        return 1;
    }

    int liquid_volume_col() const {
        return liquid_phase_index() * local_variable_count() + species_count_;
    }

    int vapor_volume_col() const {
        return vapor_phase_index() * local_variable_count() + species_count_;
    }

    std::vector<EosPhaseBlockResult> phase_blocks(const std::vector<double>& variables) const {
        require_size(variables, static_cast<std::size_t>(variable_count()), problem_name_ + " variable");
        const double pressure = variables.back();
        require_positive_finite(pressure, problem_name_ + " pressure");
        const auto amounts = pressure_route_phase_amounts(variables, species_count_);
        const auto volumes = pressure_route_phase_volumes(variables, species_count_);
        std::vector<EosPhaseBlockResult> blocks;
        blocks.reserve(static_cast<std::size_t>(phase_count()));
        for (int phase = 0; phase < phase_count(); ++phase) {
            blocks.push_back(
                evaluate_eos_phase_block(
                    args_,
                    temperature_,
                    pressure,
                    amounts[static_cast<std::size_t>(phase)],
                    volumes[static_cast<std::size_t>(phase)]
                )
            );
        }
        return blocks;
    }

    add_args args_;
    double temperature_ = 0.0;
    double initial_pressure_ = 1.0e5;
    double minimum_liquid_volume_ = kMinimumLiquidVolume;
    double maximum_liquid_volume_ = kMaximumLiquidVolume;
    double minimum_vapor_volume_ = kMinimumVaporVolume;
    double maximum_vapor_volume_ = kMaximumVaporVolume;
    double minimum_phase_volume_gap_ = kMinimumPhaseVolumeGap;
    std::vector<double> fixed_composition_;
    int fixed_phase_index_ = 0;
    std::string problem_name_;
    std::vector<double> charges_;
    double charge_constraint_tolerance_ = 0.0;
    int species_count_ = 0;
};

class NeutralFixedPressureTemperatureProblem final : public NlpProblem {
public:
    NeutralFixedPressureTemperatureProblem(
        add_args args,
        double target_pressure,
        std::vector<double> fixed_composition,
        int fixed_phase_index,
        std::string problem_name
    )
        : args_(std::move(args)),
          target_pressure_(target_pressure),
          fixed_composition_(normalized_positive_values(fixed_composition, problem_name + " composition")),
          fixed_phase_index_(fixed_phase_index),
          problem_name_(std::move(problem_name)) {
        require_positive_finite(target_pressure_, problem_name_ + " pressure");
        if (fixed_phase_index_ < 0 || fixed_phase_index_ >= phase_count()) {
            throw ValueError(problem_name_ + " fixed phase index is out of range.");
        }
        species_count_ = static_cast<int>(fixed_composition_.size());
    }

    std::string name() const override {
        return problem_name_;
    }

    int variable_count() const override {
        return phase_count() * local_variable_count() + 1;
    }

    int constraint_count() const override {
        return composition_constraint_count() + 2 * phase_count() + species_count_ + 1;
    }

    int jacobian_nonzero_count() const override {
        const int composition_nonzeros = composition_constraint_count() * species_count_;
        const int phase_total_nonzeros = phase_count() * species_count_;
        const int pressure_nonzeros = phase_count() * (local_variable_count() + 1);
        const int chemical_nonzeros = species_count_ * (phase_count() * local_variable_count() + 1);
        const int gap_nonzeros = 2;
        return composition_nonzeros + phase_total_nonzeros + pressure_nonzeros + chemical_nonzeros + gap_nonzeros;
    }

    NlpBounds bounds() const override {
        NlpBounds out;
        out.variable_lower.reserve(static_cast<std::size_t>(variable_count()));
        out.variable_upper.reserve(static_cast<std::size_t>(variable_count()));
        for (int phase = 0; phase < phase_count(); ++phase) {
            for (int species = 0; species < species_count_; ++species) {
                out.variable_lower.push_back(1.0e-14);
                out.variable_upper.push_back(10.0);
            }
            if (phase == liquid_phase_index()) {
                out.variable_lower.push_back(minimum_liquid_volume_);
                out.variable_upper.push_back(maximum_liquid_volume_);
            } else if (phase == vapor_phase_index()) {
                out.variable_lower.push_back(minimum_vapor_volume_);
                out.variable_upper.push_back(maximum_vapor_volume_);
            } else {
                throw ValueError(problem_name_ + " phase role is out of range.");
            }
        }
        out.variable_lower.push_back(minimum_temperature_);
        out.variable_upper.push_back(maximum_temperature_);
        out.constraint_lower.assign(static_cast<std::size_t>(constraint_count()), 0.0);
        out.constraint_upper.assign(static_cast<std::size_t>(constraint_count()), 0.0);
        out.constraint_lower.back() = minimum_phase_volume_gap_;
        out.constraint_upper.back() = 1.0e12;
        return out;
    }

    std::vector<double> initial_point() const override {
        const std::vector<double> root = build_temperature_route_initial_variables(
            args_,
            fixed_composition_,
            fixed_phase_index_,
            target_pressure_,
            problem_name_,
            1.0,
            DensitySeedMode::PhasePressureRoot
        );
        if (fixed_temperature_pressure_seed_satisfies_volume_bounds(root, species_count_)) {
            return root;
        }
        return build_temperature_route_initial_variables(
            args_,
            fixed_composition_,
            fixed_phase_index_,
            target_pressure_,
            problem_name_,
            1.0,
            DensitySeedMode::SeparatedPhaseRole
        );
    }

    double objective(const std::vector<double>& variables) const override {
        double value = 0.0;
        for (const TemperatureRoutePhaseBlock& block : phase_blocks(variables)) {
            value += block.block.objective;
        }
        return value;
    }

    std::vector<double> objective_gradient(const std::vector<double>& variables) const override {
        std::vector<double> out;
        out.reserve(static_cast<std::size_t>(variable_count()));
        double temperature_derivative = 0.0;
        for (const TemperatureRoutePhaseBlock& block : phase_blocks(variables)) {
            out.insert(out.end(), block.gradient.begin(), block.gradient.end() - 1);
            temperature_derivative += block.gradient.back();
        }
        out.push_back(temperature_derivative);
        return out;
    }

    std::vector<double> constraints(const std::vector<double>& variables) const override {
        const auto amounts = pressure_route_phase_amounts(variables, species_count_);
        std::vector<double> out(static_cast<std::size_t>(constraint_count()), 0.0);
        int row = 0;
        const auto& fixed_amounts = amounts[static_cast<std::size_t>(fixed_phase_index_)];
        const double fixed_total = std::accumulate(fixed_amounts.begin(), fixed_amounts.end(), 0.0);
        for (int species = 0; species < composition_constraint_count(); ++species) {
            out[static_cast<std::size_t>(row++)] =
                fixed_amounts[static_cast<std::size_t>(species)]
                - fixed_composition_[static_cast<std::size_t>(species)] * fixed_total;
        }
        for (int phase = 0; phase < phase_count(); ++phase) {
            const auto& phase_amounts = amounts[static_cast<std::size_t>(phase)];
            out[static_cast<std::size_t>(row++)] =
                std::accumulate(phase_amounts.begin(), phase_amounts.end(), 0.0) - 1.0;
        }
        const auto blocks = phase_blocks(variables);
        for (const TemperatureRoutePhaseBlock& block : blocks) {
            out[static_cast<std::size_t>(row++)] = block.block.pressure_consistency_residual;
        }
        for (int species = 0; species < species_count_; ++species) {
            out[static_cast<std::size_t>(row++)] =
                blocks[0].gradient[static_cast<std::size_t>(species)]
                - blocks[1].gradient[static_cast<std::size_t>(species)];
        }
        out[static_cast<std::size_t>(row++)] =
            variables[static_cast<std::size_t>(vapor_volume_col())]
            - variables[static_cast<std::size_t>(liquid_volume_col())];
        return out;
    }

    NlpJacobianStructure jacobian_structure() const override {
        NlpJacobianStructure out;
        out.rows.reserve(static_cast<std::size_t>(jacobian_nonzero_count()));
        out.cols.reserve(static_cast<std::size_t>(jacobian_nonzero_count()));
        int row = 0;
        for (int species = 0; species < composition_constraint_count(); ++species) {
            for (int col = 0; col < species_count_; ++col) {
                out.rows.push_back(row);
                out.cols.push_back(fixed_col(col));
            }
            ++row;
        }
        for (int phase = 0; phase < phase_count(); ++phase) {
            const int offset = phase * local_variable_count();
            for (int species = 0; species < species_count_; ++species) {
                out.rows.push_back(row);
                out.cols.push_back(offset + species);
            }
            ++row;
        }
        for (int phase = 0; phase < phase_count(); ++phase) {
            const int offset = phase * local_variable_count();
            for (int col = 0; col < local_variable_count(); ++col) {
                out.rows.push_back(row);
                out.cols.push_back(offset + col);
            }
            out.rows.push_back(row);
            out.cols.push_back(temperature_col());
            ++row;
        }
        for (int species = 0; species < species_count_; ++species) {
            for (int phase = 0; phase < phase_count(); ++phase) {
                const int offset = phase * local_variable_count();
                for (int col = 0; col < local_variable_count(); ++col) {
                    out.rows.push_back(row);
                    out.cols.push_back(offset + col);
                }
            }
            out.rows.push_back(row);
            out.cols.push_back(temperature_col());
            ++row;
        }
        out.rows.push_back(row);
        out.cols.push_back(liquid_volume_col());
        out.rows.push_back(row);
        out.cols.push_back(vapor_volume_col());
        return out;
    }

    std::vector<double> jacobian_values(const std::vector<double>& variables) const override {
        const auto blocks = phase_blocks(variables);
        std::vector<double> out;
        out.reserve(static_cast<std::size_t>(jacobian_nonzero_count()));
        int row = 0;
        for (int species = 0; species < composition_constraint_count(); ++species) {
            for (int col = 0; col < species_count_; ++col) {
                out.push_back((col == species ? 1.0 : 0.0) - fixed_composition_[static_cast<std::size_t>(species)]);
            }
            ++row;
        }
        for (int phase = 0; phase < phase_count(); ++phase) {
            for (int species = 0; species < species_count_; ++species) {
                out.push_back(1.0);
            }
            ++row;
        }
        for (int phase = 0; phase < phase_count(); ++phase) {
            const TemperatureRoutePhaseBlock& block = blocks[static_cast<std::size_t>(phase)];
            for (int col = 0; col < local_variable_count(); ++col) {
                out.push_back(block.pressure_jacobian_row_major[static_cast<std::size_t>(col)]);
            }
            out.push_back(block.pressure_jacobian_row_major.back());
            ++row;
        }
        for (int species = 0; species < species_count_; ++species) {
            double temperature_value = 0.0;
            for (int phase = 0; phase < phase_count(); ++phase) {
                const TemperatureRoutePhaseBlock& block = blocks[static_cast<std::size_t>(phase)];
                const int block_variable_count = local_variable_count() + 1;
                for (int col = 0; col < local_variable_count(); ++col) {
                    out.push_back((phase == 0 ? 1.0 : -1.0) * block.objective_hessian_row_major[
                            static_cast<std::size_t>(species * block_variable_count + col)
                        ]);
                }
                temperature_value += (phase == 0 ? 1.0 : -1.0)
                    * block.objective_hessian_row_major[
                        static_cast<std::size_t>(species * block_variable_count + block_variable_count - 1)
                    ];
            }
            out.push_back(temperature_value);
            ++row;
        }
        out.push_back(-1.0);
        out.push_back(1.0);
        return out;
    }

    bool has_exact_hessian() const override {
        return route_supports_exact_phase_derivatives(args_);
    }

    int hessian_nonzero_count() const override {
        return LagrangianHessianAssembler(variable_count()).nonzero_count();
    }

    NlpHessianStructure hessian_structure() const override {
        return LagrangianHessianAssembler(variable_count()).structure();
    }

    std::vector<double> hessian_values(
        const std::vector<double>& variables,
        double objective_factor,
        const std::vector<double>& constraint_multipliers
    ) const override {
        if (!has_exact_hessian()) {
            throw ValueError(problem_name_ + " exact Hessian requires direct non-associating phase derivatives.");
        }
        if (constraint_multipliers.size() != static_cast<std::size_t>(constraint_count())) {
            throw ValueError(problem_name_ + " Hessian multiplier vector size does not match the constraint count.");
        }

        const auto blocks = phase_blocks(variables);
        const int n = variable_count();
        const int block_variable_count = local_variable_count() + 1;
        ObjectiveSecondOrderData objective;
        objective.variable_count = n;
        objective.value = 0.0;
        objective.gradient.assign(static_cast<std::size_t>(n), 0.0);
        objective.hessian_row_major.assign(static_cast<std::size_t>(n * n), 0.0);
        objective.backend = "cppad_phase_temperature_route";
        for (int phase = 0; phase < phase_count(); ++phase) {
            const TemperatureRoutePhaseBlock& block = blocks[static_cast<std::size_t>(phase)];
            require_square_block(
                block.objective_hessian_row_major,
                block_variable_count,
                problem_name_ + " objective Hessian block"
            );
            require_third_derivative_tensor(
                block.objective_third_derivative_tensor_row_major,
                block_variable_count,
                problem_name_ + " objective third-derivative tensor"
            );
            const std::vector<int> phase_indices =
                fixed_pressure_global_indices(phase, local_variable_count(), temperature_col());
            objective.value += block.block.objective;
            for (int local = 0; local < block_variable_count; ++local) {
                objective.gradient[static_cast<std::size_t>(phase_indices[static_cast<std::size_t>(local)])] +=
                    block.gradient[static_cast<std::size_t>(local)];
            }
            add_local_hessian_to_dense(
                objective.hessian_row_major,
                n,
                phase_indices,
                block.objective_hessian_row_major,
                1.0,
                problem_name_ + " objective Hessian block"
            );
        }

        ConstraintSecondOrderData constraints_data;
        constraints_data.constraint_count = constraint_count();
        constraints_data.variable_count = n;
        constraints_data.values = this->constraints(variables);
        constraints_data.hessian_tensor_row_major.assign(
            static_cast<std::size_t>(constraint_count() * n * n),
            0.0
        );
        constraints_data.has_hessian.assign(static_cast<std::size_t>(constraint_count()), false);
        constraints_data.backend = "cppad_phase_temperature_route";

        const int pressure_row_start = composition_constraint_count() + phase_count();
        for (int phase = 0; phase < phase_count(); ++phase) {
            const TemperatureRoutePhaseBlock& block = blocks[static_cast<std::size_t>(phase)];
            require_square_block(
                block.pressure_hessian_row_major,
                block_variable_count,
                problem_name_ + " pressure Hessian block"
            );
            const int constraint_row = pressure_row_start + phase;
            constraints_data.has_hessian[static_cast<std::size_t>(constraint_row)] = true;
            add_local_hessian_to_constraint_tensor(
                constraints_data.hessian_tensor_row_major,
                constraint_row,
                n,
                fixed_pressure_global_indices(phase, local_variable_count(), temperature_col()),
                block.pressure_hessian_row_major,
                1.0,
                problem_name_ + " pressure Hessian block"
            );
        }

        const int chemical_row_start = pressure_row_start + phase_count();
        for (int species = 0; species < species_count_; ++species) {
            const int constraint_row = chemical_row_start + species;
            constraints_data.has_hessian[static_cast<std::size_t>(constraint_row)] = true;
            for (int phase = 0; phase < phase_count(); ++phase) {
                const TemperatureRoutePhaseBlock& block = blocks[static_cast<std::size_t>(phase)];
                const std::vector<double> chemical_hessian = third_derivative_slice(
                    block.objective_third_derivative_tensor_row_major,
                    species,
                    block_variable_count,
                    problem_name_ + " chemical-potential Hessian block"
                );
                add_local_hessian_to_constraint_tensor(
                    constraints_data.hessian_tensor_row_major,
                    constraint_row,
                    n,
                    fixed_pressure_global_indices(phase, local_variable_count(), temperature_col()),
                    chemical_hessian,
                    phase == 0 ? 1.0 : -1.0,
                    problem_name_ + " chemical-potential Hessian block"
                );
            }
        }

        return LagrangianHessianAssembler(n).values(
            objective_factor,
            objective,
            constraints_data,
            constraint_multipliers
        );
    }

    std::string hessian_backend() const override {
        return "cppad_phase_temperature_route";
    }

    NlpScaling scaling() const override {
        NlpScaling out;
        out.objective = 1.0;
        out.variables.assign(static_cast<std::size_t>(variable_count()), 1.0);
        out.variables.back() = 1.0e-2;
        out.constraints.assign(static_cast<std::size_t>(constraint_count()), 1.0);
        return out;
    }

    std::map<std::string, std::string> diagnostics() const override {
        return route_metadata_diagnostics(fixed_pressure_temperature_route_metadata());
    }

    int species_count() const {
        return species_count_;
    }

    int phase_count() const {
        return 2;
    }

private:
    int local_variable_count() const {
        return species_count_ + 1;
    }

    int composition_constraint_count() const {
        return std::max(0, species_count_ - 1);
    }

    int fixed_col(int species) const {
        return fixed_phase_index_ * local_variable_count() + species;
    }

    int liquid_phase_index() const {
        return 0;
    }

    int vapor_phase_index() const {
        return 1;
    }

    int liquid_volume_col() const {
        return liquid_phase_index() * local_variable_count() + species_count_;
    }

    int vapor_volume_col() const {
        return vapor_phase_index() * local_variable_count() + species_count_;
    }

    int temperature_col() const {
        return variable_count() - 1;
    }

    std::vector<TemperatureRoutePhaseBlock> phase_blocks(const std::vector<double>& variables) const {
        require_size(variables, static_cast<std::size_t>(variable_count()), problem_name_ + " variable");
        const double temperature = variables.back();
        require_positive_finite(temperature, problem_name_ + " temperature");
        const auto amounts = pressure_route_phase_amounts(variables, species_count_);
        const auto volumes = pressure_route_phase_volumes(variables, species_count_);
        std::vector<TemperatureRoutePhaseBlock> blocks;
        blocks.reserve(static_cast<std::size_t>(phase_count()));
        for (int phase = 0; phase < phase_count(); ++phase) {
            blocks.push_back(
                evaluate_temperature_route_phase_block(
                    args_,
                    temperature,
                    target_pressure_,
                    amounts[static_cast<std::size_t>(phase)],
                    volumes[static_cast<std::size_t>(phase)]
                )
            );
        }
        return blocks;
    }

    add_args args_;
    double target_pressure_ = 0.0;
    double initial_temperature_ = 300.0;
    double minimum_temperature_ = 1.0;
    double maximum_temperature_ = 2000.0;
    double minimum_liquid_volume_ = kMinimumLiquidVolume;
    double maximum_liquid_volume_ = kMaximumLiquidVolume;
    double minimum_vapor_volume_ = kMinimumVaporVolume;
    double maximum_vapor_volume_ = kMaximumVaporVolume;
    double minimum_phase_volume_gap_ = kMinimumPhaseVolumeGap;
    std::vector<double> fixed_composition_;
    int fixed_phase_index_ = 0;
    std::string problem_name_;
    int species_count_ = 0;
};

NeutralTwoPhaseEosNlpContract make_contract(const NlpProblem& problem, int phase_count, int species_count) {
    validate_nlp_problem_shape(problem);

    const std::vector<double> initial = problem.initial_point();
    const NlpBounds bounds = problem.bounds();
    const NlpJacobianStructure structure = problem.jacobian_structure();

    NeutralTwoPhaseEosNlpContract out;
    out.problem_name = problem.name();
    out.derivative_backend = "analytic_cppad";
    apply_route_metadata(out, route_metadata_from_diagnostics(problem.diagnostics()));
    out.phase_count = phase_count;
    out.species_count = species_count;
    out.variable_count = problem.variable_count();
    out.constraint_count = problem.constraint_count();
    out.jacobian_nonzero_count = problem.jacobian_nonzero_count();
    out.initial_point = initial;
    out.variable_lower_bounds = bounds.variable_lower;
    out.variable_upper_bounds = bounds.variable_upper;
    out.constraint_lower_bounds = bounds.constraint_lower;
    out.constraint_upper_bounds = bounds.constraint_upper;
    out.objective_at_initial = problem.objective(initial);
    out.gradient_at_initial = problem.objective_gradient(initial);
    out.constraints_at_initial = problem.constraints(initial);
    out.jacobian_rows = structure.rows;
    out.jacobian_cols = structure.cols;
    out.jacobian_values_at_initial = problem.jacobian_values(initial);
    return out;
}

NeutralTwoPhaseEosPostsolve fixed_temperature_pressure_postsolve(
    const add_args& args,
    double temperature,
    double pressure,
    const std::vector<std::vector<double>>& phase_amounts,
    const std::vector<double>& volumes,
    int fixed_phase_index,
    const std::vector<double>& fixed_composition,
    const std::vector<double>& charges,
    double phase_total_tolerance,
    double pressure_tolerance,
    double charge_tolerance,
    double chemical_potential_tolerance,
    double phase_distance_tolerance
) {
    NeutralTwoPhaseEosPostsolve out = evaluate_neutral_two_phase_eos_postsolve(
        args,
        temperature,
        pressure,
        phase_amounts,
        volumes,
        summed_feed_amounts(phase_amounts, static_cast<int>(fixed_composition.size())),
        phase_total_tolerance,
        pressure_tolerance,
        chemical_potential_tolerance,
        phase_distance_tolerance
    );
    apply_route_metadata(out, fixed_temperature_pressure_route_metadata(!charges.empty()));
    out.fixed_composition_norm = fixed_composition_norm(phase_amounts, fixed_phase_index, fixed_composition);
    out.phase_amount_total_norm = phase_total_norm(phase_amounts);
    out.charge_balance_norm = phase_charge_norm(phase_amounts, charges);
    out.accepted = out.accepted
        && out.fixed_composition_norm <= phase_total_tolerance
        && out.phase_amount_total_norm <= phase_total_tolerance
        && (charges.empty() || out.charge_balance_norm <= charge_tolerance);
    if (!charges.empty() && out.charge_balance_norm > charge_tolerance) {
        out.rejection_reason = "charge_balance";
    } else if (!out.accepted && out.phase_amount_total_norm > phase_total_tolerance) {
        out.rejection_reason = "phase_amount_total";
    } else if (!out.accepted && out.fixed_composition_norm > phase_total_tolerance) {
        out.rejection_reason = "fixed_composition";
    }
    return out;
}

NeutralTwoPhaseEosRouteResult solve_pressure_route(
    const add_args& args,
    double temperature,
    const std::vector<double>& fixed_composition,
    int fixed_phase_index,
    const std::string& problem_name,
    const std::vector<double>& charges,
    const IpoptSolveOptions& options,
    double phase_total_tolerance,
    double pressure_tolerance,
    double charge_tolerance,
    double chemical_potential_tolerance,
    double phase_distance_tolerance
) {
    const IpoptAdapterInfo adapter = native_ipopt_adapter_info();
    NeutralTwoPhaseEosRouteResult best;
    best.compiled = adapter.compiled;
    best.adapter_available = adapter.adapter_available;
    best.adapter_kind = adapter.adapter_kind;
    best.exact_gradient_required = adapter.exact_gradient_required;
    best.exact_jacobian_required = adapter.exact_jacobian_required;
    best.problem_name = problem_name;
    apply_route_metadata(best, fixed_temperature_pressure_route_metadata(!charges.empty()));
    if (!adapter.compiled) {
        best.status = "ipopt_dependency_required";
        return best;
    }

    const std::vector<double> normalized_composition =
        normalized_positive_values(fixed_composition, problem_name + " composition");
    const std::vector<NamedInitialVariables> seeds = pressure_route_seed_candidates(
        args,
        normalized_composition,
        fixed_phase_index,
        temperature,
        charges,
        problem_name
    );
    bool have_best = false;
    std::vector<RouteSeedAttempt> attempts;
    attempts.reserve(seeds.size() + (options.initial_variables.empty() ? 0 : 1));

    auto run_attempt = [&](const std::string& seed_name, const IpoptSolveOptions& attempt_options) {
        NeutralFixedTemperaturePressureProblem problem(
            args,
            temperature,
            normalized_composition,
            fixed_phase_index,
            problem_name,
            charges,
            charge_tolerance
        );
        const IpoptSolveResult solve = solve_ipopt_nlp(problem, attempt_options);
        NeutralTwoPhaseEosRouteResult result;
        result.compiled = adapter.compiled;
        result.adapter_available = adapter.adapter_available;
        result.adapter_kind = adapter.adapter_kind;
        result.exact_gradient_required = adapter.exact_gradient_required;
        result.exact_jacobian_required = adapter.exact_jacobian_required;
        result.problem_name = problem_name;
        result.initial_point_strategy = "deterministic_seed_sweep";
        result.seed_name = seed_name;
        result.ran = solve.solver_ran;
        const bool can_postsolve = solve.accepted || solve.feasible_point;
        result.solver_accepted = can_postsolve;
        result.solver_feasible_point = solve.feasible_point;
        result.solver_status = solve.solver_status;
        result.application_status = solve.application_status;
        apply_ipopt_solve_metadata(result, solve);
        result.objective = solve.objective;
        result.variables = solve.variables;
        result.constraints = solve.constraints;
        if (!can_postsolve) {
            result.status = "solver_rejected";
            attempts.push_back(neutral_seed_attempt_from_result(result));
            if (!have_best || neutral_route_quality(result) > neutral_route_quality(best)) {
                best = result;
                have_best = true;
            }
            return result;
        }

        const int species_count = problem.species_count();
        result.phase_amounts = pressure_route_phase_amounts(solve.variables, species_count);
        result.phase_volumes = pressure_route_phase_volumes(solve.variables, species_count);
        result.postsolve = fixed_temperature_pressure_postsolve(
            args,
            temperature,
            solve.variables.back(),
            result.phase_amounts,
            result.phase_volumes,
            fixed_phase_index,
            normalized_composition,
            charges,
            phase_total_tolerance,
            pressure_tolerance,
            charge_tolerance,
            chemical_potential_tolerance,
            phase_distance_tolerance
        );
        result.accepted = result.postsolve.accepted;
        result.status = result.accepted ? "accepted" : "postsolve_rejected";
        attempts.push_back(neutral_seed_attempt_from_result(result));
        if (!have_best || neutral_route_quality(result) > neutral_route_quality(best)) {
            best = result;
            have_best = true;
        }
        return result;
    };

    if (!options.initial_variables.empty()) {
        const NeutralTwoPhaseEosRouteResult continuation = run_attempt("continuation_state", options);
        if (continuation.accepted) {
            best.seed_attempts = attempts;
            return best;
        }
    }

    for (const auto& seed : seeds) {
        IpoptSolveOptions attempt_options = options;
        attempt_options.initial_variables = seed.variables;
        attempt_options.initial_bound_lower_multipliers.clear();
        attempt_options.initial_bound_upper_multipliers.clear();
        attempt_options.initial_constraint_multipliers.clear();
        const NeutralTwoPhaseEosRouteResult attempt = run_attempt(seed.seed_name, attempt_options);
        if (attempt.accepted) {
            break;
        }
    }

    best.initial_point_strategy = "deterministic_seed_sweep";
    best.seed_attempts = attempts;
    return best;
}

NeutralTwoPhaseEosRouteResult solve_temperature_route(
    const add_args& args,
    double target_pressure,
    const std::vector<double>& fixed_composition,
    int fixed_phase_index,
    const std::string& problem_name,
    const IpoptSolveOptions& options,
    double phase_total_tolerance,
    double pressure_tolerance,
    double chemical_potential_tolerance,
    double phase_distance_tolerance
) {
    const IpoptAdapterInfo adapter = native_ipopt_adapter_info();
    NeutralTwoPhaseEosRouteResult best;
    best.compiled = adapter.compiled;
    best.adapter_available = adapter.adapter_available;
    best.adapter_kind = adapter.adapter_kind;
    best.exact_gradient_required = adapter.exact_gradient_required;
    best.exact_jacobian_required = adapter.exact_jacobian_required;
    best.problem_name = problem_name;
    apply_route_metadata(best, fixed_pressure_temperature_route_metadata());
    if (!adapter.compiled) {
        best.status = "ipopt_dependency_required";
        return best;
    }

    const std::vector<double> normalized_composition =
        normalized_positive_values(fixed_composition, problem_name + " composition");
    const std::vector<NamedInitialVariables> seeds = temperature_route_seed_candidates(
        args,
        normalized_composition,
        fixed_phase_index,
        target_pressure,
        problem_name
    );
    bool have_best = false;
    std::vector<RouteSeedAttempt> attempts;
    attempts.reserve(seeds.size() + (options.initial_variables.empty() ? 0 : 1));

    auto run_attempt = [&](const std::string& seed_name, const IpoptSolveOptions& attempt_options) {
        NeutralFixedPressureTemperatureProblem problem(
            args,
            target_pressure,
            normalized_composition,
            fixed_phase_index,
            problem_name
        );
        const IpoptSolveResult solve = solve_ipopt_nlp(problem, attempt_options);
        NeutralTwoPhaseEosRouteResult result;
        result.compiled = adapter.compiled;
        result.adapter_available = adapter.adapter_available;
        result.adapter_kind = adapter.adapter_kind;
        result.exact_gradient_required = adapter.exact_gradient_required;
        result.exact_jacobian_required = adapter.exact_jacobian_required;
        result.problem_name = problem_name;
        result.initial_point_strategy = "deterministic_seed_sweep";
        result.seed_name = seed_name;
        result.ran = solve.solver_ran;
        const bool can_postsolve = solve.accepted || solve.feasible_point;
        result.solver_accepted = can_postsolve;
        result.solver_feasible_point = solve.feasible_point;
        result.solver_status = solve.solver_status;
        result.application_status = solve.application_status;
        apply_ipopt_solve_metadata(result, solve);
        result.objective = solve.objective;
        result.variables = solve.variables;
        result.constraints = solve.constraints;
        if (!can_postsolve) {
            result.status = "solver_rejected";
            attempts.push_back(neutral_seed_attempt_from_result(result));
            if (!have_best || neutral_route_quality(result) > neutral_route_quality(best)) {
                best = result;
                have_best = true;
            }
            return result;
        }

        const int species_count = problem.species_count();
        const double solved_temperature = solve.variables.back();
        result.phase_amounts = pressure_route_phase_amounts(solve.variables, species_count);
        result.phase_volumes = pressure_route_phase_volumes(solve.variables, species_count);
        result.postsolve = fixed_temperature_pressure_postsolve(
            args,
            solved_temperature,
            target_pressure,
            result.phase_amounts,
            result.phase_volumes,
            fixed_phase_index,
            normalized_composition,
            {},
            phase_total_tolerance,
            pressure_tolerance,
            0.0,
            chemical_potential_tolerance,
            phase_distance_tolerance
        );
        result.accepted = result.postsolve.accepted;
        result.status = result.accepted ? "accepted" : "postsolve_rejected";
        attempts.push_back(neutral_seed_attempt_from_result(result));
        if (!have_best || neutral_route_quality(result) > neutral_route_quality(best)) {
            best = result;
            have_best = true;
        }
        return result;
    };

    if (!options.initial_variables.empty()) {
        const NeutralTwoPhaseEosRouteResult continuation = run_attempt("continuation_state", options);
        if (continuation.accepted) {
            best.seed_attempts = attempts;
            return best;
        }
    }

    for (const auto& seed : seeds) {
        IpoptSolveOptions attempt_options = options;
        attempt_options.initial_variables = seed.variables;
        attempt_options.initial_bound_lower_multipliers.clear();
        attempt_options.initial_bound_upper_multipliers.clear();
        attempt_options.initial_constraint_multipliers.clear();
        const NeutralTwoPhaseEosRouteResult attempt = run_attempt(seed.seed_name, attempt_options);
        if (attempt.accepted) {
            break;
        }
    }

    best.initial_point_strategy = "deterministic_seed_sweep";
    best.seed_attempts = attempts;
    return best;
}

}  // namespace

NeutralTwoPhaseEosNlpContract evaluate_neutral_bubble_p_eos_nlp_contract(
    const add_args& args,
    double temperature,
    const std::vector<double>& liquid_composition
) {
    NeutralFixedTemperaturePressureProblem problem(
        args,
        temperature,
        liquid_composition,
        0,
        "neutral_bubble_p_eos"
    );
    return make_contract(problem, problem.phase_count(), problem.species_count());
}

NeutralTwoPhaseEosNlpContract evaluate_electrolyte_bubble_p_eos_nlp_contract(
    const add_args& args,
    double temperature,
    const std::vector<double>& liquid_composition
) {
    NeutralFixedTemperaturePressureProblem problem(
        args,
        temperature,
        liquid_composition,
        0,
        "electrolyte_bubble_p_eos",
        args.z
    );
    return make_contract(problem, problem.phase_count(), problem.species_count());
}

NeutralTwoPhaseEosNlpContract evaluate_neutral_dew_p_eos_nlp_contract(
    const add_args& args,
    double temperature,
    const std::vector<double>& vapor_composition
) {
    NeutralFixedTemperaturePressureProblem problem(
        args,
        temperature,
        vapor_composition,
        1,
        "neutral_dew_p_eos"
    );
    return make_contract(problem, problem.phase_count(), problem.species_count());
}

NeutralTwoPhaseEosNlpContract evaluate_neutral_bubble_t_eos_nlp_contract(
    const add_args& args,
    double target_pressure,
    const std::vector<double>& liquid_composition
) {
    NeutralFixedPressureTemperatureProblem problem(
        args,
        target_pressure,
        liquid_composition,
        0,
        "neutral_bubble_t_eos"
    );
    return make_contract(problem, problem.phase_count(), problem.species_count());
}

NeutralTwoPhaseEosNlpContract evaluate_neutral_dew_t_eos_nlp_contract(
    const add_args& args,
    double target_pressure,
    const std::vector<double>& vapor_composition
) {
    NeutralFixedPressureTemperatureProblem problem(
        args,
        target_pressure,
        vapor_composition,
        1,
        "neutral_dew_t_eos"
    );
    return make_contract(problem, problem.phase_count(), problem.species_count());
}

NeutralTwoPhaseEosRouteResult solve_neutral_bubble_p_eos_route(
    const add_args& args,
    double temperature,
    const std::vector<double>& liquid_composition,
    const IpoptSolveOptions& options,
    double phase_total_tolerance,
    double pressure_tolerance,
    double chemical_potential_tolerance,
    double phase_distance_tolerance
) {
    return solve_pressure_route(
        args,
        temperature,
        liquid_composition,
        0,
        "neutral_bubble_p_eos",
        {},
        options,
        phase_total_tolerance,
        pressure_tolerance,
        0.0,
        chemical_potential_tolerance,
        phase_distance_tolerance
    );
}

NeutralTwoPhaseEosRouteResult solve_electrolyte_bubble_p_eos_route(
    const add_args& args,
    double temperature,
    const std::vector<double>& liquid_composition,
    const IpoptSolveOptions& options,
    double phase_total_tolerance,
    double pressure_tolerance,
    double charge_tolerance,
    double chemical_potential_tolerance,
    double phase_distance_tolerance
) {
    return solve_pressure_route(
        args,
        temperature,
        liquid_composition,
        0,
        "electrolyte_bubble_p_eos",
        args.z,
        options,
        phase_total_tolerance,
        pressure_tolerance,
        charge_tolerance,
        chemical_potential_tolerance,
        phase_distance_tolerance
    );
}

NeutralTwoPhaseEosRouteResult solve_neutral_dew_p_eos_route(
    const add_args& args,
    double temperature,
    const std::vector<double>& vapor_composition,
    const IpoptSolveOptions& options,
    double phase_total_tolerance,
    double pressure_tolerance,
    double chemical_potential_tolerance,
    double phase_distance_tolerance
) {
    return solve_pressure_route(
        args,
        temperature,
        vapor_composition,
        1,
        "neutral_dew_p_eos",
        {},
        options,
        phase_total_tolerance,
        pressure_tolerance,
        0.0,
        chemical_potential_tolerance,
        phase_distance_tolerance
    );
}

NeutralTwoPhaseEosRouteResult solve_neutral_bubble_t_eos_route(
    const add_args& args,
    double target_pressure,
    const std::vector<double>& liquid_composition,
    const IpoptSolveOptions& options,
    double phase_total_tolerance,
    double pressure_tolerance,
    double chemical_potential_tolerance,
    double phase_distance_tolerance
) {
    return solve_temperature_route(
        args,
        target_pressure,
        liquid_composition,
        0,
        "neutral_bubble_t_eos",
        options,
        phase_total_tolerance,
        pressure_tolerance,
        chemical_potential_tolerance,
        phase_distance_tolerance
    );
}

NeutralTwoPhaseEosRouteResult solve_neutral_dew_t_eos_route(
    const add_args& args,
    double target_pressure,
    const std::vector<double>& vapor_composition,
    const IpoptSolveOptions& options,
    double phase_total_tolerance,
    double pressure_tolerance,
    double chemical_potential_tolerance,
    double phase_distance_tolerance
) {
    return solve_temperature_route(
        args,
        target_pressure,
        vapor_composition,
        1,
        "neutral_dew_t_eos",
        options,
        phase_total_tolerance,
        pressure_tolerance,
        chemical_potential_tolerance,
        phase_distance_tolerance
    );
}

}  // namespace epcsaft::native::equilibrium_nlp
