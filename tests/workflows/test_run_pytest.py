import json
import shutil
import subprocess
import sys
from pathlib import Path

import run_pytest
from scripts import build_plot_manifest
from scripts import codex_check


def test_confidence_slice_extends_generic_targets_without_changing_generic():
    pytest_temp = Path("build") / "pytest-temp" / "run-test"

    generic_args = run_pytest._pytest_args(["-q"], pytest_temp, generic=True)
    confidence_args = run_pytest._pytest_args(["-q"], pytest_temp, confidence=True)

    assert generic_args[: len(run_pytest.GENERIC_TEST_TARGETS)] == list(run_pytest.GENERIC_TEST_TARGETS)
    assert "tests/native/test_runtime_contracts.py" not in generic_args
    assert confidence_args[: len(run_pytest.CONFIDENCE_TEST_TARGETS)] == list(run_pytest.CONFIDENCE_TEST_TARGETS)
    assert "tests/native/test_runtime_contracts.py" in confidence_args
    assert confidence_args[-3:] == ["-q", "--basetemp", str(pytest_temp)]


def test_codex_check_modes_route_to_agent_facing_validation_bundles():
    assert codex_check.CHECK_COMMANDS["quick"] == (
        ("scripts/codex_doctor.py",),
        ("run_pytest.py", "--generic", "-q"),
    )
    assert ("scripts/build_plot_manifest.py", "--check") in codex_check.CHECK_COMMANDS["confidence"]
    assert ("scripts/build_plot_manifest.py", "--refresh") in codex_check.CHECK_COMMANDS["plots"]
    assert ("run_pytest.py", "--equilibrium-confidence", "-q") in codex_check.CHECK_COMMANDS["full"]


def test_plot_manifest_validation_rejects_html_and_duplicate_outputs(tmp_path):
    manifest = tmp_path / "manifest.json"
    item = {field: "" for field in build_plot_manifest.REQUIRED_FIELDS}
    item.update(
        {
            "path": "../../scripts/example/out/figure.html",
            "output_path": "scripts/example/out/figure.html",
            "source_folder": "scripts",
        }
    )
    manifest.write_text(
        json.dumps(build_plot_manifest.manifest_payload([item, item])),
        encoding="utf-8",
    )

    errors = build_plot_manifest.validate_manifest(manifest)

    assert any("Duplicate output_path" in error for error in errors)
    assert any("should not reference HTML" in error for error in errors)


def test_named_shortcuts_expand_to_expected_targets_and_keep_pytest_arg_ordering():
    pytest_temp = Path("build") / "pytest-temp" / "run-test"

    runtime_args = run_pytest._pytest_args(["-q"], pytest_temp, runtime=True)
    api_args = run_pytest._pytest_args(["-q"], pytest_temp, api=True)
    native_args = run_pytest._pytest_args(["-q"], pytest_temp, native=True)
    equilibrium_confidence_args = run_pytest._pytest_args(["-q"], pytest_temp, equilibrium_confidence=True)
    profile_args = run_pytest._pytest_args(["-q"], pytest_temp, profile=True)
    profile_full_args = run_pytest._pytest_args(["-q"], pytest_temp, profile_full=True)
    plots_args = run_pytest._pytest_args(["-q"], pytest_temp, plots=True)

    assert runtime_args[: len(run_pytest.RUNTIME_TEST_TARGETS)] == list(run_pytest.RUNTIME_TEST_TARGETS)
    assert api_args[: len(run_pytest.API_TEST_TARGETS)] == list(run_pytest.API_TEST_TARGETS)
    assert native_args[: len(run_pytest.NATIVE_TEST_TARGETS)] == list(run_pytest.NATIVE_TEST_TARGETS)
    assert equilibrium_confidence_args[: len(run_pytest.EQUILIBRIUM_CONFIDENCE_TEST_TARGETS)] == list(
        run_pytest.EQUILIBRIUM_CONFIDENCE_TEST_TARGETS
    )
    assert profile_args[: len(run_pytest.PROFILE_TEST_TARGETS)] == list(run_pytest.PROFILE_TEST_TARGETS)
    assert profile_full_args[: len(run_pytest.FULL_PROFILE_TEST_TARGETS)] == list(run_pytest.FULL_PROFILE_TEST_TARGETS)
    assert plots_args[: len(run_pytest.PLOT_TEST_TARGETS)] == list(run_pytest.PLOT_TEST_TARGETS)
    assert runtime_args[-3:] == ["-q", "--basetemp", str(pytest_temp)]
    assert api_args[-3:] == ["-q", "--basetemp", str(pytest_temp)]
    assert native_args[-3:] == ["-q", "--basetemp", str(pytest_temp)]
    assert equilibrium_confidence_args[-3:] == ["-q", "--basetemp", str(pytest_temp)]
    assert profile_args[-3:] == ["-q", "--basetemp", str(pytest_temp)]
    assert profile_full_args[-3:] == ["-q", "--basetemp", str(pytest_temp)]
    assert plots_args[-3:] == ["-q", "--basetemp", str(pytest_temp)]


def test_plot_slice_stays_out_of_generic_and_confidence_targets():
    assert "tests/plots/test_equilibrium_plot_outputs.py" in run_pytest.PLOT_TEST_TARGETS
    assert "tests/plots/test_property_plot_outputs.py" in run_pytest.PLOT_TEST_TARGETS
    assert "tests/plots/test_contribution_plot_outputs.py" in run_pytest.PLOT_TEST_TARGETS
    assert "tests/plots/test_regression_plot_outputs.py" in run_pytest.PLOT_TEST_TARGETS
    assert "tests/plots/test_native_plot_outputs.py" in run_pytest.PLOT_TEST_TARGETS
    assert "tests/plots/test_api_parity_plot_outputs.py" in run_pytest.PLOT_TEST_TARGETS
    assert ("tests/plots/test_" + "plot" + "ly_backfill.py") not in run_pytest.PLOT_TEST_TARGETS
    assert all(target.startswith("tests/plots/") for target in run_pytest.PLOT_TEST_TARGETS)
    for target in run_pytest.PLOT_TEST_TARGETS:
        assert target not in run_pytest.GENERIC_TEST_TARGETS
        assert target not in run_pytest.CONFIDENCE_TEST_TARGETS


def test_slice_targets_use_grouped_test_subpackages():
    all_targets = [
        *run_pytest.GENERIC_TEST_TARGETS,
        *run_pytest.CONFIDENCE_TEST_TARGETS,
        *run_pytest.RUNTIME_TEST_TARGETS,
        *run_pytest.API_TEST_TARGETS,
        *run_pytest.NATIVE_TEST_TARGETS,
        *run_pytest.EQUILIBRIUM_CONFIDENCE_TEST_TARGETS,
        *run_pytest.PROFILE_TEST_TARGETS,
        *run_pytest.FULL_PROFILE_TEST_TARGETS,
        *run_pytest.PLOT_TEST_TARGETS,
    ]

    assert all(target.startswith("tests/") for target in all_targets)
    assert all(target.count("/") >= 2 for target in all_targets)
    assert "tests/api/test_runtime.py" in run_pytest.API_TEST_TARGETS
    assert "tests/equilibrium/test_vle.py" in run_pytest.GENERIC_TEST_TARGETS
    assert "tests/native/test_equation_registry.py" in run_pytest.GENERIC_TEST_TARGETS
    assert "tests/regression/test_hydrocarbon.py" in run_pytest.GENERIC_TEST_TARGETS
    assert "tests/workflows/test_run_pytest.py" not in run_pytest.GENERIC_TEST_TARGETS


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
    assert "appended and will run in addition" in result.stdout
