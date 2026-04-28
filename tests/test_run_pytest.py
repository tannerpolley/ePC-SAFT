import shutil
from pathlib import Path

import run_pytest


def test_confidence_slice_extends_generic_targets_without_changing_generic():
    pytest_temp = Path("build") / "pytest-temp" / "run-test"

    generic_args = run_pytest._pytest_args(["-q"], pytest_temp, generic=True)
    confidence_args = run_pytest._pytest_args(["-q"], pytest_temp, confidence=True)

    assert generic_args[:len(run_pytest.GENERIC_TEST_TARGETS)] == list(run_pytest.GENERIC_TEST_TARGETS)
    assert "tests/test_native_runtime_contracts.py" not in generic_args
    assert confidence_args[:len(run_pytest.CONFIDENCE_TEST_TARGETS)] == list(run_pytest.CONFIDENCE_TEST_TARGETS)
    assert "tests/test_native_runtime_contracts.py" in confidence_args
    assert confidence_args[-3:] == ["-q", "--basetemp", str(pytest_temp)]


def test_named_shortcuts_expand_to_expected_targets_and_keep_pytest_arg_ordering():
    pytest_temp = Path("build") / "pytest-temp" / "run-test"

    runtime_args = run_pytest._pytest_args(["-q"], pytest_temp, runtime=True)
    api_args = run_pytest._pytest_args(["-q"], pytest_temp, api=True)
    native_args = run_pytest._pytest_args(["-q"], pytest_temp, native=True)
    profile_args = run_pytest._pytest_args(["-q"], pytest_temp, profile=True)
    profile_full_args = run_pytest._pytest_args(["-q"], pytest_temp, profile_full=True)

    assert runtime_args[:len(run_pytest.RUNTIME_TEST_TARGETS)] == list(run_pytest.RUNTIME_TEST_TARGETS)
    assert api_args[:len(run_pytest.API_TEST_TARGETS)] == list(run_pytest.API_TEST_TARGETS)
    assert native_args[:len(run_pytest.NATIVE_TEST_TARGETS)] == list(run_pytest.NATIVE_TEST_TARGETS)
    assert profile_args[:len(run_pytest.PROFILE_TEST_TARGETS)] == list(run_pytest.PROFILE_TEST_TARGETS)
    assert profile_full_args[:len(run_pytest.FULL_PROFILE_TEST_TARGETS)] == list(run_pytest.FULL_PROFILE_TEST_TARGETS)
    assert runtime_args[-3:] == ["-q", "--basetemp", str(pytest_temp)]
    assert api_args[-3:] == ["-q", "--basetemp", str(pytest_temp)]
    assert native_args[-3:] == ["-q", "--basetemp", str(pytest_temp)]
    assert profile_args[-3:] == ["-q", "--basetemp", str(pytest_temp)]
    assert profile_full_args[-3:] == ["-q", "--basetemp", str(pytest_temp)]


def test_profile_shortcut_sets_perf_environment_flag(monkeypatch):
    pytest_temp = Path("build") / "pytest-temp" / "run-test"
    monkeypatch.delenv("EPCSAFT_RUN_PERF", raising=False)

    env = run_pytest._pytest_env(pytest_temp, profile=True)

    assert env["EPCSAFT_RUN_PERF"] == "1"
    assert env["ePCSAFT_RUN_PERF"] == "1"


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
