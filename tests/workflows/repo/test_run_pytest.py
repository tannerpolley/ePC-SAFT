import shutil
import subprocess
import sys
from pathlib import Path

import run_pytest
from scripts.dev import doctor, validate_project


def test_confidence_slice_extends_generic_targets_without_changing_generic():
    pytest_temp = Path("build") / "pytest-temp" / "run-test"

    generic_args = run_pytest._pytest_args(["-q"], pytest_temp, generic=True)
    confidence_args = run_pytest._pytest_args(["-q"], pytest_temp, confidence=True)

    assert generic_args[: len(run_pytest.GENERIC_TEST_TARGETS)] == list(run_pytest.GENERIC_TEST_TARGETS)
    assert "tests/native/runtime/test_runtime_density_closure.py" not in generic_args
    assert confidence_args[: len(run_pytest.CONFIDENCE_TEST_TARGETS)] == list(run_pytest.CONFIDENCE_TEST_TARGETS)
    assert (
        "tests/native/runtime/test_runtime_density_closure.py::"
        "test_pressure_based_and_density_based_states_match_for_ionic_system"
    ) in confidence_args
    assert confidence_args[-3:] == ["-q", "--basetemp", str(pytest_temp)]


def test_default_pytest_route_uses_fast_contracts_not_exhaustive_suite():
    pytest_temp = Path("build") / "pytest-temp" / "run-test"

    default_args = run_pytest._pytest_args(["-q"], pytest_temp)

    assert default_args[: len(run_pytest.FAST_TEST_TARGETS)] == list(run_pytest.FAST_TEST_TARGETS)
    assert "tests" not in default_args
    assert not any(target.startswith("tests/plots/") for target in default_args)
    assert "tests/workflows/validation/equilibrium_core/test_electrolyte_lle_confidence.py" not in default_args
    assert default_args[-3:] == ["-q", "--basetemp", str(pytest_temp)]


def test_predefined_pytest_targets_reference_existing_files():
    repo_root = Path(__file__).resolve().parents[3]
    missing: set[str] = set()

    for targets in run_pytest.SLICE_TARGETS.values():
        for target in targets:
            if target == "tests":
                continue
            target_path = target.split("::", 1)[0]
            if not (repo_root / target_path).exists():
                missing.add(target_path)

    assert not sorted(missing)


def test_all_shortcut_is_the_explicit_exhaustive_pytest_route():
    pytest_temp = Path("build") / "pytest-temp" / "run-test"

    all_args = run_pytest._pytest_args(["-q"], pytest_temp, all_tests=True)

    assert all_args[: len(run_pytest.ALL_TEST_TARGETS)] == list(run_pytest.ALL_TEST_TARGETS)
    assert all_args == ["tests", "-q", "--basetemp", str(pytest_temp)]


def test_validate_project_modes_route_to_standard_validation_bundles():
    assert validate_project.CHECK_COMMANDS["quick"] == (
        ("scripts/dev/doctor.py",),
        ("run_pytest.py", "-q"),
    )
    assert all(
        "build_plot_" + "manifest.py" not in command
        for mode in validate_project.CHECK_COMMANDS.values()
        for command in mode
    )
    assert "plots" not in validate_project.CHECK_COMMANDS
    assert validate_project.CHECK_COMMANDS["full"] == (
        ("scripts/dev/doctor.py",),
        ("run_pytest.py", "--all", "-q"),
    )
    assert validate_project.CHECK_COMMANDS["ceres-cppad"] == (
        ("scripts/dev/build_epcsaft.py", "--profile", "full"),
        (
            "run_pytest.py",
            "tests/native/ceres/test_ceres_pure_regression.py",
            "tests/native/ceres/test_ceres_binary_regression.py",
            "-q",
        ),
    )


def test_doctor_tracks_native_symbols_added_by_recent_workflows():
    required = set(doctor.REQUIRED_CORE_SYMBOLS)

    assert "_fit_pure_neutral_native_ceres" in required
    assert "_fit_pure_neutral_native_least_squares" not in required
    assert "_fit_generic_native_least_squares" not in required
    assert "_fit_generic_native_ceres" in required
    assert "_evaluate_generic_native_debug" in required
    assert "_solve_equilibrium_native" in required
    assert "_native_ipopt_smoke" in required


def test_native_regression_source_has_no_eigen_nonlinear_optimizer_route():
    source = Path("src/epcsaft/native/epcsaft_regression.cpp").read_text(encoding="utf-8")

    blocked_terms = (
        "unsupported/Eigen/" + "Levenberg" + "Marquardt",
        "Levenberg" + "Marquardt",
        "Numerical" + "Diff",
    )
    for term in blocked_terms:
        assert term not in source


def test_native_ceres_sources_have_no_numeric_diff_route():
    source = "\n".join(
        path.read_text(encoding="utf-8")
        for path in sorted(Path("src/epcsaft/native").rglob("*"))
        if path.suffix in {".cpp", ".h", ".hpp"}
    )
    source_lower = source.lower()

    blocked_terms = (
        "Numeric" + "Diff",
        "Numeric" + "DiffCostFunction",
        "Dynamic" + "Numeric" + "DiffCostFunction",
        "CERES_" + "NUMERIC" + "_DIFF",
        "numeric" + "_diff",
        "finite" + "_difference",
    )
    for term in blocked_terms:
        assert term.lower() not in source_lower


def test_native_cppad_source_does_not_assemble_removed_backend_status():
    source = "\n".join(
        path.read_text(encoding="utf-8")
        for path in (
            Path("src/epcsaft/native/cppad/cppad_derivative_result.h"),
            Path("src/epcsaft/native/cppad/cppad_scalar.h"),
            Path("src/epcsaft/native/cppad/cppad_smoke_checks.cpp"),
        )
    )

    removed = 'std::string("backend") + "_" + "unavailable"'
    assert removed not in source


def test_doctor_recommends_ninja_migration_for_mingw_build_tree():
    command = doctor._ninja_migration_recommendation("MinGW Makefiles", "C:/tools/ninja.exe")

    assert command == "uv run python scripts\\dev\\build_epcsaft.py --clean --generator ninja"


def test_doctor_does_not_recommend_ninja_migration_when_already_ninja():
    command = doctor._ninja_migration_recommendation("Ninja", "C:/tools/ninja.exe")

    assert command is None


def test_named_shortcuts_expand_to_expected_targets_and_keep_pytest_arg_ordering():
    pytest_temp = Path("build") / "pytest-temp" / "run-test"

    runtime_args = run_pytest._pytest_args(["-q"], pytest_temp, runtime=True)
    api_args = run_pytest._pytest_args(["-q"], pytest_temp, api=True)
    native_args = run_pytest._pytest_args(["-q"], pytest_temp, native=True)
    equilibrium_confidence_args = run_pytest._pytest_args(["-q"], pytest_temp, equilibrium_confidence=True)
    equilibrium_api_args = run_pytest._pytest_args(["-q"], pytest_temp, equilibrium_api=True)
    profile_args = run_pytest._pytest_args(["-q"], pytest_temp, profile=True)
    profile_full_args = run_pytest._pytest_args(["-q"], pytest_temp, profile_full=True)

    assert runtime_args[: len(run_pytest.RUNTIME_TEST_TARGETS)] == list(run_pytest.RUNTIME_TEST_TARGETS)
    assert api_args[: len(run_pytest.API_TEST_TARGETS)] == list(run_pytest.API_TEST_TARGETS)
    assert native_args[: len(run_pytest.NATIVE_TEST_TARGETS)] == list(run_pytest.NATIVE_TEST_TARGETS)
    assert equilibrium_confidence_args[: len(run_pytest.EQUILIBRIUM_CONFIDENCE_TEST_TARGETS)] == list(
        run_pytest.EQUILIBRIUM_CONFIDENCE_TEST_TARGETS
    )
    assert equilibrium_api_args[: len(run_pytest.EQUILIBRIUM_API_TEST_TARGETS)] == list(
        run_pytest.EQUILIBRIUM_API_TEST_TARGETS
    )
    assert profile_args[: len(run_pytest.PROFILE_TEST_TARGETS)] == list(run_pytest.PROFILE_TEST_TARGETS)
    assert profile_full_args[: len(run_pytest.FULL_PROFILE_TEST_TARGETS)] == list(run_pytest.FULL_PROFILE_TEST_TARGETS)
    assert runtime_args[-3:] == ["-q", "--basetemp", str(pytest_temp)]
    assert api_args[-3:] == ["-q", "--basetemp", str(pytest_temp)]
    assert native_args[-3:] == ["-q", "--basetemp", str(pytest_temp)]
    assert equilibrium_confidence_args[-3:] == ["-q", "--basetemp", str(pytest_temp)]
    assert equilibrium_api_args[-3:] == ["-q", "--basetemp", str(pytest_temp)]
    assert profile_args[-3:] == ["-q", "--basetemp", str(pytest_temp)]
    assert profile_full_args[-3:] == ["-q", "--basetemp", str(pytest_temp)]


def test_plot_output_tests_have_no_named_slice():
    assert "plots" not in run_pytest.SLICE_TARGETS
    assert all(
        not target.startswith("tests/plots/") for targets in run_pytest.SLICE_TARGETS.values() for target in targets
    )


def test_slice_targets_use_grouped_test_subpackages():
    all_targets = [
        *run_pytest.GENERIC_TEST_TARGETS,
        *run_pytest.CONFIDENCE_TEST_TARGETS,
        *run_pytest.RUNTIME_TEST_TARGETS,
        *run_pytest.API_TEST_TARGETS,
        *run_pytest.NATIVE_TEST_TARGETS,
        *run_pytest.EQUILIBRIUM_CONFIDENCE_TEST_TARGETS,
        *run_pytest.EQUILIBRIUM_API_TEST_TARGETS,
        *run_pytest.PROFILE_TEST_TARGETS,
        *run_pytest.FULL_PROFILE_TEST_TARGETS,
    ]

    assert all(target.startswith("tests/") for target in all_targets)
    assert all(target.count("/") >= 2 for target in all_targets)
    assert "tests/api/runtime/test_runtime_exports_and_metadata.py" in run_pytest.API_TEST_TARGETS
    assert (
        "tests/equilibrium/core/test_vle.py::"
        "test_ternary_hydrocarbon_basis_tp_flash_closes_material_and_fugacity_balance"
    ) in run_pytest.GENERIC_TEST_TARGETS
    assert "tests/equilibrium/core/test_lle.py::test_methanol_cyclohexane_lle_flash_solves_seeded_phase_split" in (
        run_pytest.GENERIC_TEST_TARGETS
    )
    assert "tests/equilibrium/core/test_lle.py::test_lle_flash_requested_ipopt_requires_native_adapter" in (
        run_pytest.GENERIC_TEST_TARGETS
    )
    assert "tests/native/contracts/test_equation_registry.py::test_equation_registry_outputs_are_synced" in (
        run_pytest.GENERIC_TEST_TARGETS
    )
    assert "tests/regression/core/test_hydrocarbon.py::test_methane_reference_parameters_keep_native_objective_pinned" in (
        run_pytest.GENERIC_TEST_TARGETS
    )
    assert "tests/workflows/repo/test_run_pytest.py" not in run_pytest.GENERIC_TEST_TARGETS


def test_profile_shortcut_sets_perf_environment_flag(monkeypatch):
    pytest_temp = Path("build") / "pytest-temp" / "run-test"
    monkeypatch.delenv("EPCSAFT_RUN_PERF", raising=False)

    env = run_pytest._pytest_env(pytest_temp, profile=True)

    assert env["EPCSAFT_RUN_PERF"] == "1"
    assert env["ePCSAFT_RUN_PERF"] == "1"


def test_equilibrium_confidence_slice_is_listed():
    result = subprocess.run(
        [sys.executable, "run_pytest.py", "--equilibrium-confidence", "--list-slices"],
        capture_output=True,
        text=True,
        check=False,
    )

    assert result.returncode == 0
    assert "equilibrium-confidence:" in result.stdout


def test_equilibrium_api_slice_is_listed():
    result = subprocess.run(
        [sys.executable, "run_pytest.py", "--equilibrium-api", "--list-slices"],
        capture_output=True,
        text=True,
        check=False,
    )

    assert result.returncode == 0
    assert "equilibrium-api:" in result.stdout


def test_equilibrium_confidence_shortcut_keeps_full_report_env_opt_in():
    source = Path(run_pytest.__file__).read_text(encoding="utf-8")

    assert 'env["EPCSAFT_EQUILIBRIUM_CONFIDENCE"] = "1"' not in source


def test_full_profile_shortcut_sets_perf_environment_flag(monkeypatch):
    pytest_temp = Path("build") / "pytest-temp" / "run-test"
    monkeypatch.delenv("EPCSAFT_RUN_PERF", raising=False)

    env = run_pytest._pytest_env(pytest_temp, profile=True)

    assert env["EPCSAFT_RUN_PERF"] == "1"
    assert env["ePCSAFT_RUN_PERF"] == "1"


def test_pytest_temp_root_prefers_configured_root_and_normalizes_relative_paths(monkeypatch, tmp_path):
    repo_root = tmp_path / "repo"
    repo_root.mkdir()
    monkeypatch.setenv("EPCSAFT_PYTEST_TEMP_ROOT", "external-temp")

    pytest_temp = run_pytest._pytest_temp(repo_root)

    assert pytest_temp.is_dir()
    assert pytest_temp.parent == repo_root / "external-temp" / "pytest-temp"
    shutil.rmtree(repo_root / "external-temp", ignore_errors=True)


def test_pytest_env_overrides_perf_flag_only_for_profile_modes(monkeypatch):
    pytest_temp = Path("build") / "pytest-temp" / "run-test"
    monkeypatch.setenv("EPCSAFT_RUN_PERF", "0")
    monkeypatch.setenv("ePCSAFT_RUN_PERF", "0")

    normal_env = run_pytest._pytest_env(pytest_temp, profile=False)
    profile_env = run_pytest._pytest_env(pytest_temp, profile=True)

    assert normal_env["EPCSAFT_RUN_PERF"] == "0"
    assert profile_env["EPCSAFT_RUN_PERF"] == "1"
    assert profile_env["ePCSAFT_RUN_PERF"] == "1"


def test_full_profile_runtime_note_sets_expected_timeout_floor():
    assert run_pytest.FULL_PROFILE_MIN_TIMEOUT_SECONDS >= 120
    assert "about a minute" in run_pytest.FULL_PROFILE_RUNTIME_NOTE
    assert "allow at least 120 seconds" in run_pytest.FULL_PROFILE_RUNTIME_NOTE


def test_slice_listing_text_names_all_targets():
    listing = run_pytest._slice_listing_text()

    assert run_pytest.SLICE_SELECTION_NOTE in listing
    for name, targets in run_pytest.SLICE_TARGETS.items():
        assert f"{name}:" in listing
        for target in targets:
            assert target in listing


def test_list_slices_exits_without_running_pytest():
    result = subprocess.run(
        [sys.executable, "run_pytest.py", "--list-slices"],
        capture_output=True,
        text=True,
        check=False,
    )

    assert result.returncode == 0, result.stdout + result.stderr
    assert "Available slices:" in result.stdout
    assert "Running:" not in result.stdout


def test_slice_selection_note_uses_dev_validate_project_path():
    assert "`uv run python scripts/dev/validate_project.py quick`" in run_pytest.SLICE_SELECTION_NOTE
    assert "`uv run python scripts/validate_project.py quick`" not in run_pytest.SLICE_SELECTION_NOTE


def test_help_mentions_slice_append_semantics():
    result = subprocess.run(
        [sys.executable, "run_pytest.py", "--help"],
        capture_output=True,
        text=True,
        check=False,
    )

    assert result.returncode == 0, result.stdout + result.stderr
    assert "Slice flags are mutually exclusive" in result.stdout
    assert "Extra positional pytest targets" in result.stdout
    assert "exhaustive historical test suite" in result.stdout
