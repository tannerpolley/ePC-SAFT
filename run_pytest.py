import argparse
import os
import shutil
import sys
import uuid
from pathlib import Path

GENERIC_TEST_TARGETS = (
    "tests/api/package/test_package_main.py::test_python_m_epcsaft_reports_package_and_core_status",
    "tests/api/runtime/test_runtime_exports_and_metadata.py::test_package_exports_are_available",
    "tests/api/runtime/test_runtime_neutral_methods.py::test_neutral_scalar_methods_return_expected_values",
    "tests/api/runtime/test_runtime_ionic_methods.py::test_ionic_activity_and_solution_methods_return_expected_values",
    "tests/api/parameters/test_parameter_templates.py::test_runtime_options_accept_cppad_modes_and_preserve_explicit_overrides",
    "tests/api/regression/test_regression_api_native_backends.py::test_public_pure_neutral_regression_is_robust_to_distinct_initial_guesses",
    "tests/regression/core/test_hydrocarbon.py::test_methane_reference_parameters_keep_native_objective_pinned",
    "tests/equilibrium/core/test_vle.py::test_tp_flash_builds_one_native_route_request_before_ipopt_gate",
    "tests/equilibrium/core/test_lle.py::test_lle_flash_builds_one_native_route_request_before_ipopt_gate",
    "tests/equilibrium/core/test_stability.py::test_stability_uses_native_ipopt_route_after_validation",
    (
        "tests/equilibrium/electrolyte/test_electrolyte_lle_smokes.py::"
        "test_electrolyte_lle_builds_native_route_before_ipopt_gate"
    ),
    (
        "tests/workflows/validation/equilibrium_core/test_electrolyte_thermo_diagnostics.py::"
        "test_khudaida_fixture_loads_charge_neutral_explicit_ions"
    ),
    "tests/native/contracts/test_equilibrium_native_contracts.py::test_native_equilibrium_entrypoint_is_exposed",
    "tests/native/runtime/test_runtime_density_closure.py::test_pressure_based_and_density_based_states_match_for_neutral_system",
    "tests/native/contracts/test_equation_registry.py::test_equation_registry_outputs_are_synced",
    "tests/workflows/repo/test_project_structure.py",
    "tests/workflows/repo/test_run_pytest.py::test_list_slices_exits_without_running_pytest",
    "tests/workflows/repo/test_workflow_entrypoints.py::test_docs_make_confidence_suite_the_default_runtime_check",
)
FAST_TEST_TARGETS = GENERIC_TEST_TARGETS
CONFIDENCE_TEST_TARGETS = (
    *GENERIC_TEST_TARGETS,
    "tests/native/runtime/test_runtime_density_closure.py::test_pressure_based_and_density_based_states_match_for_ionic_system",
    "tests/native/runtime/test_runtime_contribution_contracts.py::test_native_residual_helmholtz_and_compressibility_contributions_match_neutral_contract",
    "tests/native/contracts/test_equilibrium_native_contracts.py::test_public_tp_flash_requires_native_ipopt_route",
)
EQUILIBRIUM_CONFIDENCE_TEST_TARGETS = (
    (
        "tests/workflows/validation/equilibrium_core/test_electrolyte_lle_confidence.py::"
        "test_khudaida_benchmark_fixture_loads_charge_neutral_cases"
    ),
    (
        "tests/workflows/validation/equilibrium_core/test_electrolyte_thermo_diagnostics.py::"
        "test_khudaida_package_tieline_fixed_phase_residual_is_internally_consistent"
    ),
)
EQUILIBRIUM_API_TEST_TARGETS = (
    "tests/equilibrium/core/test_vle.py::test_tp_flash_builds_one_native_route_request_before_ipopt_gate",
    "tests/equilibrium/core/test_lle.py::test_lle_flash_builds_one_native_route_request_before_ipopt_gate",
    "tests/equilibrium/core/test_stability.py::test_stability_uses_native_ipopt_route_after_validation",
    (
        "tests/equilibrium/electrolyte/test_electrolyte_lle_smokes.py::"
        "test_electrolyte_lle_builds_native_route_before_ipopt_gate"
    ),
    "tests/api/runtime/test_runtime_exports_and_metadata.py::test_runtime_build_info_and_capabilities_are_json_like",
    (
        "tests/api/reactive/test_reactive_speciation_results.py::"
        "test_solve_reactive_speciation_activity_coupled_state_requires_native_ipopt_route"
    ),
    "tests/api/reactive/test_reactive_speciation_options.py::test_reactive_speciation_options_expose_jacobian_backend_selector",
    "tests/api/reactive/test_reactive_speciation_options.py::test_reactive_speciation_requested_ipopt_routes_ideal_speciation_when_compiled",
    "tests/api/reactive/test_reactive_electrolyte_bubble_setup.py",
    "tests/api/reactive/test_reactive_electrolyte_bubble_results.py",
    (
        "tests/native/equilibrium/test_chemical_equilibrium_native_api.py::"
        "test_native_chemical_equilibrium_residual_evaluator_uses_analytic_jacobian_by_default"
    ),
    (
        "tests/native/equilibrium/test_chemical_equilibrium_native_errors.py::"
        "test_native_chemical_equilibrium_residual_evaluator_rejects_removed_backend"
    ),
)
ALL_TEST_TARGETS = ("tests",)
RUNTIME_TEST_TARGETS = (
    "tests/api/runtime/test_runtime_exports_and_metadata.py",
    "tests/api/runtime/test_runtime_neutral_methods.py",
    "tests/api/runtime/test_runtime_ionic_methods.py",
    "tests/native/runtime/test_runtime_density_closure.py",
    "tests/native/runtime/test_runtime_contribution_contracts.py",
    "tests/native/runtime/test_runtime_cache_contracts.py",
)
API_TEST_TARGETS = (
    "tests/api/runtime/test_runtime_exports_and_metadata.py",
    "tests/api/runtime/test_runtime_neutral_methods.py",
    "tests/api/runtime/test_runtime_ionic_methods.py",
    "tests/api/parameters/test_parameter_templates.py",
    "tests/api/regression/test_regression_api_public_contracts.py",
    "tests/api/regression/test_regression_api_native_backends.py",
    "tests/api/regression/test_regression_api_results_and_errors.py",
)
NATIVE_TEST_TARGETS = (
    "tests/native/runtime/test_runtime_density_closure.py",
    "tests/native/runtime/test_runtime_contribution_contracts.py",
    "tests/native/runtime/test_runtime_cache_contracts.py",
)
PROFILE_TEST_TARGETS = ("tests/profile/test_runtime_profile.py",)
FULL_PROFILE_TEST_TARGETS = (
    "tests/profile/test_runtime_profile.py",
    "tests/profile/test_miac_profile.py",
    "tests/profile/test_regression_profile.py",
)
SLICE_TARGETS = {
    "generic": GENERIC_TEST_TARGETS,
    "all": ALL_TEST_TARGETS,
    "confidence": CONFIDENCE_TEST_TARGETS,
    "equilibrium-confidence": EQUILIBRIUM_CONFIDENCE_TEST_TARGETS,
    "equilibrium-api": EQUILIBRIUM_API_TEST_TARGETS,
    "runtime": RUNTIME_TEST_TARGETS,
    "api": API_TEST_TARGETS,
    "native": NATIVE_TEST_TARGETS,
    "profile": PROFILE_TEST_TARGETS,
    "profile-full": FULL_PROFILE_TEST_TARGETS,
}
FULL_PROFILE_MIN_TIMEOUT_SECONDS = 120
FULL_PROFILE_RUNTIME_NOTE = (
    "--profile-full runs runtime, MIAC, and regression profiling; "
    f"it can take about a minute locally, so allow at least {FULL_PROFILE_MIN_TIMEOUT_SECONDS} seconds."
)
SLICE_SELECTION_NOTE = (
    "Slice flags are mutually exclusive. Developers should normally start with "
    "`uv run python scripts/dev/validate_project.py quick` or `uv run python run_pytest.py -q`. "
    "Use `--all` only when you explicitly need the exhaustive historical suite. "
    "Extra positional pytest targets after a slice are appended and will run in addition to that slice."
)


def _repo_root() -> Path:
    return Path(__file__).resolve().parent


def _pytest_temp(repo_root: Path) -> Path:
    configured_root = os.environ.get("EPCSAFT_PYTEST_TEMP_ROOT")
    if configured_root is not None:
        root = Path(configured_root).expanduser()
        if not root.is_absolute():
            root = (repo_root / root).resolve()
        root = root / "pytest-temp"
    else:
        root = repo_root / "build" / "pytest-temp"

    path = root / f"run-{os.getpid()}-{uuid.uuid4().hex[:8]}"
    path.mkdir(parents=True, exist_ok=True)
    return path


def _pytest_env(pytest_temp: Path, profile: bool = False) -> dict[str, str]:
    env = os.environ.copy()
    env["TMP"] = str(pytest_temp.resolve())
    env["TEMP"] = str(pytest_temp.resolve())
    env["TMPDIR"] = str(pytest_temp.resolve())
    try:
        from scripts.dev.native_runtime_env import apply_native_runtime_env
    except ModuleNotFoundError:
        apply_native_runtime_env = None
    if apply_native_runtime_env is not None:
        apply_native_runtime_env(env)
    if profile:
        env["EPCSAFT_RUN_PERF"] = "1"
        env["ePCSAFT_RUN_PERF"] = "1"
    return env


def _pytest_args(
    pytest_args: list[str],
    pytest_temp: Path,
    generic: bool = False,
    confidence: bool = False,
    equilibrium_confidence: bool = False,
    equilibrium_api: bool = False,
    runtime: bool = False,
    api: bool = False,
    native: bool = False,
    profile: bool = False,
    profile_full: bool = False,
    all_tests: bool = False,
) -> list[str]:
    cmd: list[str] = []
    has_predefined_targets = (
        generic
        or confidence
        or equilibrium_confidence
        or equilibrium_api
        or runtime
        or api
        or native
        or profile
        or profile_full
        or all_tests
    )
    if all_tests:
        cmd.extend(ALL_TEST_TARGETS)
    elif confidence:
        cmd.extend(CONFIDENCE_TEST_TARGETS)
    elif equilibrium_confidence:
        cmd.extend(EQUILIBRIUM_CONFIDENCE_TEST_TARGETS)
    elif equilibrium_api:
        cmd.extend(EQUILIBRIUM_API_TEST_TARGETS)
    elif generic:
        cmd.extend(GENERIC_TEST_TARGETS)
    elif runtime:
        cmd.extend(RUNTIME_TEST_TARGETS)
    elif api:
        cmd.extend(API_TEST_TARGETS)
    elif native:
        cmd.extend(NATIVE_TEST_TARGETS)
    elif profile:
        cmd.extend(PROFILE_TEST_TARGETS)
    elif profile_full:
        cmd.extend(FULL_PROFILE_TEST_TARGETS)

    if has_predefined_targets:
        cmd.extend(pytest_args)
    else:
        has_positional_target = any(not arg.startswith("-") for arg in pytest_args)
        if has_positional_target:
            cmd.extend(pytest_args)
        else:
            cmd.extend(FAST_TEST_TARGETS)
            cmd.extend(pytest_args)

    if not any(arg == "--basetemp" or arg.startswith("--basetemp=") for arg in cmd):
        cmd.extend(["--basetemp", str(pytest_temp)])
    return cmd


def _slice_listing_text() -> str:
    lines = [SLICE_SELECTION_NOTE, "", "Available slices:"]
    for name, targets in SLICE_TARGETS.items():
        lines.append(f"{name}:")
        for target in targets:
            lines.append(f"  {target}")
    return "\n".join(lines)


def _patch_windows_pytest_temp_acl() -> None:
    if os.name != "nt":
        return

    original_mkdir = Path.mkdir

    def sandbox_safe_mkdir(self, mode=0o777, parents=False, exist_ok=False):
        return original_mkdir(self, mode=0o777, parents=parents, exist_ok=exist_ok)

    Path.mkdir = sandbox_safe_mkdir


def _patch_pytest_cleanup() -> None:
    # Pytest's Windows cleanup hook can trip over restricted ACLs after
    # tmp_path tests pass. Keep the per-run basetemp, but skip that hook.
    try:
        import _pytest.pathlib as pytest_pathlib
        import _pytest.tmpdir as pytest_tmpdir

        pytest_pathlib.cleanup_dead_symlinks = lambda root: None
        pytest_tmpdir.cleanup_dead_symlinks = lambda root: None
    except Exception:
        pass


def _failure_message(pytest_temp: Path) -> str:
    cleanup_path = str(pytest_temp.resolve())
    return (
        "Pytest failed; keeping temp directory for triage: "
        f"{cleanup_path}\n"
        "Cleanup with: "
        f"Remove-Item -Recurse -Force '{cleanup_path}' (PowerShell)\n"
        f"or rm -rf {cleanup_path} (POSIX shells)"
    )


def main() -> int:
    parser = argparse.ArgumentParser(description=SLICE_SELECTION_NOTE)
    predefined = parser.add_mutually_exclusive_group()
    predefined.add_argument("--generic", action="store_true", help="Run the core generic test slice")
    predefined.add_argument(
        "--confidence",
        action="store_true",
        help="Run generic fast contracts plus native runtime contract tests",
    )
    predefined.add_argument(
        "--equilibrium-confidence",
        action="store_true",
        help="Run electrolyte equilibrium confidence contract tests; full reports remain env opt-in",
    )
    predefined.add_argument(
        "--equilibrium-api",
        action="store_true",
        help="Run fast equilibrium/speciation API tests for downstream-agent workflows",
    )
    predefined.add_argument("--runtime", action="store_true", help="Run runtime API and native contract tests")
    predefined.add_argument("--api", action="store_true", help="Run public API and regression API tests")
    predefined.add_argument("--native", action="store_true", help="Run native runtime contract tests")
    predefined.add_argument("--profile", action="store_true", help="Run the opt-in runtime-only profile test")
    predefined.add_argument(
        "--profile-full",
        action="store_true",
        help="Run all opt-in runtime, MIAC, and regression profile tests",
    )
    predefined.add_argument(
        "--all",
        dest="all_tests",
        action="store_true",
        help="Run the exhaustive historical test suite under tests/; this is intentionally opt-in",
    )
    parser.add_argument(
        "--list-slices", action="store_true", help="Print named test slices and exit without running pytest"
    )
    args, pytest_args = parser.parse_known_args()

    if args.list_slices:
        print(_slice_listing_text())
        return 0

    repo_root = _repo_root()
    pytest_temp = _pytest_temp(repo_root)
    env = _pytest_env(pytest_temp, profile=args.profile or args.profile_full)
    src_root = repo_root / "src"
    sys.path.insert(0, str(src_root))
    env["PYTHONPATH"] = str(src_root)

    cmd = _pytest_args(
        pytest_args,
        pytest_temp,
        args.generic,
        confidence=args.confidence,
        equilibrium_confidence=args.equilibrium_confidence,
        equilibrium_api=args.equilibrium_api,
        runtime=args.runtime,
        api=args.api,
        native=args.native,
        profile=args.profile,
        profile_full=args.profile_full,
        all_tests=args.all_tests,
    )
    print("Running:", f"{sys.executable} -m pytest", " ".join(cmd), flush=True)
    if args.profile_full:
        print(FULL_PROFILE_RUNTIME_NOTE, flush=True)
    os.environ.update(env)

    _patch_windows_pytest_temp_acl()
    _patch_pytest_cleanup()

    import pytest

    exit_code = int(pytest.main(cmd))
    if exit_code == 0:
        shutil.rmtree(pytest_temp, ignore_errors=True)
    else:
        print(_failure_message(pytest_temp), flush=True)
    return exit_code


if __name__ == "__main__":
    raise SystemExit(main())
