#pragma once

#include <cstddef>
#include <map>
#include <string>
#include <vector>

namespace epcsaft::native::equilibrium_nlp {

struct RouteMetadata {
    std::string variable_model;
    std::string density_backend;
    std::vector<std::string> residual_families;
    std::vector<std::string> constraint_families;
};

inline std::string join_route_families(const std::vector<std::string>& families) {
    std::string out;
    for (std::size_t index = 0; index < families.size(); ++index) {
        if (index > 0) {
            out += ",";
        }
        out += families[index];
    }
    return out;
}

inline std::vector<std::string> split_route_families(const std::string& text) {
    std::vector<std::string> out;
    std::string current;
    for (char ch : text) {
        if (ch == ',') {
            if (!current.empty()) {
                out.push_back(current);
            }
            current.clear();
        } else if (ch != ' ') {
            current.push_back(ch);
        }
    }
    if (!current.empty()) {
        out.push_back(current);
    }
    return out;
}

inline std::map<std::string, std::string> route_metadata_diagnostics(const RouteMetadata& metadata) {
    return {
        {"variable_model", metadata.variable_model},
        {"density_backend", metadata.density_backend},
        {"residual_blocks", join_route_families(metadata.residual_families)},
        {"constraint_blocks", join_route_families(metadata.constraint_families)},
    };
}

inline RouteMetadata route_metadata_from_diagnostics(
    const std::map<std::string, std::string>& diagnostics
) {
    RouteMetadata out;
    const auto variable_model = diagnostics.find("variable_model");
    if (variable_model != diagnostics.end()) {
        out.variable_model = variable_model->second;
    }
    const auto density_backend = diagnostics.find("density_backend");
    if (density_backend != diagnostics.end()) {
        out.density_backend = density_backend->second;
    }
    const auto residual_blocks = diagnostics.find("residual_blocks");
    if (residual_blocks != diagnostics.end()) {
        out.residual_families = split_route_families(residual_blocks->second);
    }
    const auto constraint_blocks = diagnostics.find("constraint_blocks");
    if (constraint_blocks != diagnostics.end()) {
        out.constraint_families = split_route_families(constraint_blocks->second);
    }
    return out;
}

inline RouteMetadata phase_amount_volume_route_metadata(
    bool has_charge_constraints,
    bool has_phase_distance_constraint
) {
    RouteMetadata out;
    out.variable_model = "phase_species_amounts_plus_phase_volume";
    out.density_backend = "explicit_phase_volume_pressure_constraint";
    out.residual_families = {
        "material_balance",
        "phase_pressure_consistency",
        "phase_equilibrium",
        "phase_distance",
    };
    if (has_charge_constraints) {
        out.residual_families.push_back("phase_charge");
    }
    out.constraint_families = {
        "material_balance",
        "phase_pressure_consistency",
    };
    if (has_charge_constraints) {
        out.constraint_families.push_back("phase_charge");
    }
    if (has_phase_distance_constraint) {
        out.constraint_families.push_back("phase_distance");
    }
    return out;
}

inline RouteMetadata electrolyte_liquid_root_route_metadata(const std::string& variable_model) {
    RouteMetadata out;
    out.variable_model = variable_model;
    out.density_backend = "explicit_log_density_pressure_constraint";
    out.residual_families = {"phase_equilibrium", "material_balance"};
    out.constraint_families = {
        "phase_equilibrium",
        "phase_pressure_consistency",
        "phase_distance",
        "formula_feasibility",
    };
    return out;
}

inline RouteMetadata reactive_liquid_root_route_metadata(
    bool phase_tagged_reaction_constraints_active,
    bool phase_charge_constraints_active
) {
    RouteMetadata out;
    out.variable_model = "log_phase_species_amounts_plus_log_density";
    out.density_backend = "explicit_log_density_pressure_constraint";
    out.residual_families = {
        "conserved_balance",
        "reaction_stationarity",
        "phase_equilibrium",
        "phase_charge",
    };
    out.constraint_families = {"conserved_balance"};
    if (phase_tagged_reaction_constraints_active) {
        out.constraint_families.push_back("reaction_stationarity");
    }
    if (phase_charge_constraints_active) {
        out.constraint_families.push_back("phase_charge");
    }
    out.constraint_families.push_back("phase_pressure_consistency");
    out.constraint_families.push_back("phase_distance");
    return out;
}

inline RouteMetadata reactive_phase_amount_volume_route_metadata(bool has_charge_constraints) {
    RouteMetadata out;
    out.variable_model = "phase_species_amounts_plus_phase_volume";
    out.density_backend = "explicit_phase_volume_pressure_constraint";
    out.residual_families = {
        "conserved_balance",
        "reaction_stationarity",
        "phase_pressure_consistency",
        "phase_equilibrium",
        "phase_distance",
    };
    if (has_charge_constraints) {
        out.residual_families.push_back("phase_charge");
    }
    out.constraint_families = {
        "conserved_balance",
        "phase_pressure_consistency",
    };
    if (has_charge_constraints) {
        out.constraint_families.push_back("phase_charge");
    }
    return out;
}

inline RouteMetadata fixed_temperature_pressure_route_metadata(bool has_charge_constraints) {
    RouteMetadata out;
    out.variable_model = "phase_species_amounts_plus_phase_volume_plus_pressure";
    out.density_backend = "explicit_phase_volume_pressure_constraint";
    out.residual_families = {
        "fixed_composition",
        "phase_amount_total",
        "phase_pressure_consistency",
        "phase_equilibrium",
        "phase_distance",
    };
    if (has_charge_constraints) {
        out.residual_families.push_back("phase_charge");
    }
    out.constraint_families = {
        "fixed_composition",
        "phase_amount_total",
    };
    if (has_charge_constraints) {
        out.constraint_families.push_back("phase_charge");
    }
    out.constraint_families.push_back("phase_pressure_consistency");
    out.constraint_families.push_back("phase_equilibrium");
    out.constraint_families.push_back("phase_volume_gap");
    return out;
}

inline RouteMetadata fixed_pressure_temperature_route_metadata() {
    RouteMetadata out = fixed_temperature_pressure_route_metadata(false);
    out.variable_model = "phase_species_amounts_plus_phase_volume_plus_temperature";
    return out;
}

}  // namespace epcsaft::native::equilibrium_nlp
