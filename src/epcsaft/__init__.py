"""Public package interface for the native ePC-SAFT runtime."""

from __future__ import annotations

import os
from importlib import import_module

_DLL_DIRECTORY_HANDLES: list[object] = []


def _add_runtime_dll_directories() -> None:
    if os.name != "nt" or not hasattr(os, "add_dll_directory"):
        return
    raw_dirs = os.environ.get("EPCSAFT_RUNTIME_DLL_DIRS", "")
    for raw_dir in raw_dirs.split(os.pathsep):
        dll_dir = raw_dir.strip()
        if not dll_dir:
            continue
        _DLL_DIRECTORY_HANDLES.append(os.add_dll_directory(dll_dir))


# Register external native dependency directories before imports that load _core.
_add_runtime_dll_directories()

_EXPORT_GROUPS = {
    "epcsaft._types": ("ActivityCoefficientResult", "InputError", "SolutionError"),
    "epcsaft.dataset_validation": ("validate_dataset_bundle",),
    "epcsaft.electrolyte_bubble": ("ElectrolyteBubbleOptions", "ElectrolyteBubbleResult"),
    "epcsaft.epcsaft": ("ePCSAFTMixture", "ePCSAFTState"),
    "epcsaft.equilibrium": (
        "BubblePoint",
        "DewPoint",
        "ElectrolyteBubblePoint",
        "ElectrolyteLLEProblem",
        "EquilibriumOptions",
        "EquilibriumPhase",
        "EquilibriumProblem",
        "EquilibriumResult",
        "LLEProblem",
        "ReactiveElectrolyteBubbleProblem",
        "ReactivePhaseEquilibriumProblem",
        "ReactiveSpeciationProblem",
        "StabilityAnalysis",
        "StabilityResult",
        "StabilityTrial",
        "TPFlash",
        "bubble_p",
        "bubble_t",
        "dew_p",
        "dew_t",
        "electrolyte_feed_from_molality",
    ),
    "epcsaft.implicit_sensitivity": (
        "ImplicitSolveResult",
        "implicit_backend_for_residual_backend",
        "implicit_sensitivity_from_jacobians",
    ),
    "epcsaft.parameter_schema": (
        "AssociationSite",
        "BinaryRecord",
        "ComponentIdentifier",
        "ParameterSet",
        "PermittivityRecord",
        "PureRecord",
    ),
    "epcsaft.parameter_templates": ("create_parameter_template",),
    "epcsaft.parameters": (
        "DATASET_ROOT",
        "available_datasets",
        "get_prop_dict",
        "molality_to_molefraction",
        "molefraction_to_molality",
    ),
    "epcsaft.properties": ("evaluate_fugacity_coefficients", "evaluate_fugacity_coefficients_batch"),
    "epcsaft.reactive_electrolyte": (
        "ReactiveElectrolyteBubbleOptions",
        "ReactiveElectrolyteBubbleResult",
        "solve_reactive_electrolyte_bubble",
        "solve_reactive_electrolyte_bubble_sweep",
    ),
    "epcsaft.reactive_regression": (
        "ReactiveElectrolyteBatch",
        "ReactiveElectrolyteBatchOptions",
        "ReactiveElectrolyteBatchResult",
        "ReactiveElectrolyteRegressionContext",
        "ReactiveElectrolyteRow",
        "ReactiveElectrolyteRowResult",
        "ReactiveRegressionFitResult",
        "ReactiveRegressionJacobianResult",
        "ReactiveRegressionObjective",
        "ReactiveRegressionObjectiveResult",
        "build_reactive_regression_objective",
        "evaluate_reactive_regression_objective",
        "fit_reactive_electrolyte_parameters",
        "summarize_regression_result",
        "write_regression_parameter_table",
        "write_regression_residual_table",
        "write_regression_row_table",
        "write_regression_summary",
    ),
    "epcsaft.reactive_speciation": (
        "ReactionConstantConvention",
        "ReactionDefinition",
        "ReactiveSpeciationOptions",
        "ReactiveSpeciationResult",
        "solve_reactive_speciation",
        "solve_reactive_speciation_sweep",
    ),
    "epcsaft.reactive_staged": ("ReactiveStagedEquilibriumResult", "solve_reactive_staged_equilibrium"),
    "epcsaft.regression": (
        "BinaryInteraction",
        "FitBounds",
        "FitParameter",
        "FitProblem",
        "FitResult",
        "FitTerm",
        "RelativePermittivityResidual",
        "TargetDataset",
        "TargetRow",
        "evaluate_pure_neutral_derivatives",
        "fit_binary_pair",
        "fit_binary_parameters",
        "fit_liquid_electrolyte_parameters",
        "fit_pure_ion",
        "fit_pure_neutral",
        "fit_pure_parameters",
        "load_regression_records",
        "validate_regression_provenance",
        "write_fit_result",
    ),
    "epcsaft.runtime": ("__git_commit__", "__version__", "capabilities", "runtime_build_info"),
}
_EXPORT_MODULES = {name: module_name for module_name, names in _EXPORT_GROUPS.items() for name in names}
__all__ = sorted(_EXPORT_MODULES)


def __getattr__(name: str) -> object:
    module_name = _EXPORT_MODULES.get(name)
    if module_name is None:
        raise AttributeError(f"module {__name__!r} has no attribute {name!r}")
    value = getattr(import_module(module_name), name)
    globals()[name] = value
    return value


def __dir__() -> list[str]:
    return sorted(set(globals()) | set(__all__))
