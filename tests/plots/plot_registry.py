from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Iterable


@dataclass(frozen=True)
class PlotRecipe:
    source_tests: tuple[str, ...]
    plot_targets: tuple[str, ...]
    output_dirs: tuple[str, ...]


class PlotRecipeLookupError(ValueError):
    pass


PLOT_RECIPES: tuple[PlotRecipe, ...] = (
    PlotRecipe(
        source_tests=(
            "tests/equilibrium/test_api.py",
            "tests/equilibrium/test_vle.py",
        ),
        plot_targets=("tests/plots/test_equilibrium_plot_outputs.py",),
        output_dirs=("equilibrium/vle",),
    ),
    PlotRecipe(
        source_tests=("tests/equilibrium/test_lle.py",),
        plot_targets=("tests/plots/test_equilibrium_plot_outputs.py",),
        output_dirs=("equilibrium/lle",),
    ),
    PlotRecipe(
        source_tests=("tests/equilibrium/test_stability.py",),
        plot_targets=("tests/plots/test_equilibrium_plot_outputs.py",),
        output_dirs=("equilibrium/stability",),
    ),
    PlotRecipe(
        source_tests=(
            "tests/equilibrium/test_electrolyte_lle.py",
            "tests/equilibrium/test_electrolyte_thermo_diagnostics.py",
        ),
        plot_targets=("tests/plots/test_equilibrium_plot_outputs.py",),
        output_dirs=("equilibrium/electrolyte_lle",),
    ),
    PlotRecipe(
        source_tests=(
            "tests/api/test_runtime.py",
            "tests/api/test_package_main.py",
        ),
        plot_targets=("tests/plots/test_api_parity_plot_outputs.py",),
        output_dirs=("api/parity",),
    ),
    PlotRecipe(
        source_tests=(
            "tests/native/test_equation_registry.py",
            "tests/native/test_runtime_contracts.py",
        ),
        plot_targets=("tests/plots/test_native_plot_outputs.py",),
        output_dirs=("native/branches", "native/derivatives"),
    ),
    PlotRecipe(
        source_tests=("tests/regression/test_hydrocarbon.py",),
        plot_targets=("tests/plots/test_regression_plot_outputs.py",),
        output_dirs=("regression/gradients", "regression/hydrocarbon"),
    ),
    PlotRecipe(
        source_tests=("tests/api/test_runtime.py",),
        plot_targets=("tests/plots/test_property_plot_outputs.py",),
        output_dirs=("properties/activity_fugacity", "properties/residual_energy"),
    ),
    PlotRecipe(
        source_tests=("tests/api/test_runtime.py",),
        plot_targets=("tests/plots/test_contribution_plot_outputs.py",),
        output_dirs=("contributions/ionic", "contributions/neutral"),
    ),
)


def normalize_test_target(target: str | Path) -> str:
    raw = str(target).split("::", 1)[0].strip()
    normalized = raw.replace("\\", "/").strip("/")
    if normalized.startswith("./"):
        normalized = normalized[2:]
    return normalized


def all_plot_recipes() -> tuple[PlotRecipe, ...]:
    return PLOT_RECIPES


def resolve_plot_recipes(targets: Iterable[str | Path], *, all_recipes: bool = False) -> tuple[PlotRecipe, ...]:
    if all_recipes:
        return all_plot_recipes()

    normalized_targets = tuple(normalize_test_target(target) for target in targets if str(target).strip())
    matched: list[PlotRecipe] = []
    for recipe in PLOT_RECIPES:
        if any(_target_matches_recipe(target, recipe) for target in normalized_targets):
            matched.append(recipe)

    if not matched:
        raise PlotRecipeLookupError(missing_recipe_message(normalized_targets))
    return tuple(dict.fromkeys(matched))


def producer_targets(recipes: Iterable[PlotRecipe]) -> tuple[str, ...]:
    targets: list[str] = []
    for recipe in recipes:
        targets.extend(recipe.plot_targets)
    return tuple(dict.fromkeys(targets))


def output_dirs(recipes: Iterable[PlotRecipe]) -> tuple[str, ...]:
    dirs: list[str] = []
    for recipe in recipes:
        dirs.extend(recipe.output_dirs)
    return tuple(dict.fromkeys(dirs))


def missing_recipe_message(targets: Iterable[str]) -> str:
    target_list = tuple(targets)
    target = target_list[0] if target_list else "tests/path/test_feature.py"
    return (
        "No plot recipe is registered for the requested test target(s): "
        f"{', '.join(target_list) or '<none>'}.\n"
        "Add an explicit recipe in tests/plots/plot_registry.py, for example:\n\n"
        "    PlotRecipe(\n"
        f"        source_tests=(\"{target}\",),\n"
        "        plot_targets=(\"tests/plots/test_<feature>_plot_outputs.py\",),\n"
        "        output_dirs=(\"<gallery-subfolder>\",),\n"
        "    )\n\n"
        "The workflow does not synthesize placeholder scientific plots."
    )


def _target_matches_recipe(target: str, recipe: PlotRecipe) -> bool:
    if not target:
        return False
    target_path = target.rstrip("/")
    for source in recipe.source_tests:
        if target_path == source:
            return True
        if source.startswith(f"{target_path}/"):
            return True
    return False
