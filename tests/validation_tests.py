"""Compatibility shim for validation plotting.

Canonical scripts:
- validation_2014_repro.py
- validation_version_analysis.py
"""

from validation_version_analysis import run_validation_version_analysis


def run_validation_tests():
    run_validation_version_analysis()


def test_validation_tests_compat():
    run_validation_tests()


if __name__ == "__main__":
    run_validation_tests()
