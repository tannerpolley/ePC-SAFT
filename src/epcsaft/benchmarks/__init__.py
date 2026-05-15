"""Compatibility shims for benchmark helpers.

Benchmark execution lives in :mod:`scripts.benchmarks.helpers`; the runtime
package keeps these imports only for older callers during the cleanup window.
"""

from scripts.benchmarks.helpers.literature import (
    LITERATURE_CASES,
    render_literature_benchmark_table,
    run_literature_benchmarks,
)
from scripts.benchmarks.helpers.neutral_equilibrium import CASE_BUILDERS, run_neutral_equilibrium_benchmarks
from scripts.benchmarks.helpers.reactive_regression import run_reactive_regression_benchmarks

__all__ = [
    "CASE_BUILDERS",
    "LITERATURE_CASES",
    "render_literature_benchmark_table",
    "run_literature_benchmarks",
    "run_neutral_equilibrium_benchmarks",
    "run_reactive_regression_benchmarks",
]
