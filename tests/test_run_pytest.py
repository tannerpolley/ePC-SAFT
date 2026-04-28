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
