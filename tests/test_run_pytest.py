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
