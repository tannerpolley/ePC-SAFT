from __future__ import annotations

import json
import subprocess
from collections import Counter, OrderedDict
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any

from epcsaft.runtime import __git_commit__, __version__

REPO_ROOT = Path(__file__).resolve().parents[3]

ALLOWED_CLASSIFICATIONS = {
    "implemented",
    "already_supported_with_tests",
    "blocker_requires_followup",
    "out_of_scope_by_roadmap",
}


@dataclass(frozen=True)
class LiteratureBenchmarkEntry:
    case: str
    title: str
    classification: str
    coverage_kind: str
    package_surface: tuple[str, ...]
    validation_paths: tuple[str, ...]
    notes: str

    def to_payload(self) -> dict[str, Any]:
        payload = asdict(self)
        payload["package_surface"] = list(self.package_surface)
        payload["validation_paths"] = list(self.validation_paths)
        return payload


def _git_commit() -> str:
    try:
        completed = subprocess.run(
            ["git", "rev-parse", "HEAD"],
            cwd=REPO_ROOT,
            capture_output=True,
            text=True,
            check=True,
        )
    except (OSError, subprocess.SubprocessError):
        return str(__git_commit__)
    commit = completed.stdout.strip()
    return commit or str(__git_commit__)


LITERATURE_CASES: OrderedDict[str, LiteratureBenchmarkEntry] = OrderedDict(
    (
        (
            "mea_simple_workflow",
            LiteratureBenchmarkEntry(
                case="mea_simple_workflow",
                title="MEA simple workflow benchmark",
                classification="already_supported_with_tests",
                coverage_kind="smoke_regression",
                package_surface=(
                    "tests/regression/literature/test_mea_co2_h2o_pure_parameter_benchmark.py",
                    "tests/regression/literature/test_literature_pure_parameter_regression.py",
                ),
                validation_paths=("tests/regression/literature/test_literature_pure_parameter_regression.py",),
                notes=(
                    "Current coverage exists through a private benchmark-only regression path rather than a "
                    "documented application-specific public API."
                ),
            ),
        ),
        (
            "mdea_epcsaft",
            LiteratureBenchmarkEntry(
                case="mdea_epcsaft",
                title="MDEA ePC-SAFT benchmark",
                classification="blocker_requires_followup",
                coverage_kind="inventory_only",
                package_surface=(),
                validation_paths=(),
                notes=(
                    "No package-owned MDEA literature fixture or workflow test is wired into the current benchmark "
                    "surface."
                ),
            ),
        ),
        (
            "figiel_2025_ssm_ds_born",
            LiteratureBenchmarkEntry(
                case="figiel_2025_ssm_ds_born",
                title="Figiel 2025 SSM+DS Born benchmark",
                classification="already_supported_with_tests",
                coverage_kind="smoke_regression",
                package_surface=(
                    "tests/regression/literature/test_figiel_2025_born_parameter_parity.py",
                    "tests/regression/electrolyte/test_miac_liquid_electrolyte_regression.py",
                    "tests/regression/electrolyte/test_miac_liquid_electrolyte_parity.py",
                ),
                validation_paths=("tests/regression/literature/test_figiel_2025_born_parameter_parity.py",),
                notes="Liquid-electrolyte Born parity and backend-unavailable contracts are already covered by tests.",
            ),
        ),
        (
            "held_2014_revised_epcsaft",
            LiteratureBenchmarkEntry(
                case="held_2014_revised_epcsaft",
                title="Held 2014 revised ePC-SAFT benchmark",
                classification="blocker_requires_followup",
                coverage_kind="inventory_only",
                package_surface=(),
                validation_paths=(),
                notes="The current repo does not expose a package-owned Held 2014 benchmark fixture or dedicated test.",
            ),
        ),
        (
            "non_electrolyte_lle",
            LiteratureBenchmarkEntry(
                case="non_electrolyte_lle",
                title="non-electrolyte LLE benchmark",
                classification="blocker_requires_followup",
                coverage_kind="inventory_only",
                package_surface=("src/epcsaft/benchmarks/neutral_equilibrium.py",),
                validation_paths=(),
                notes=(
                    "A generic seeded neutral LLE benchmark exists, but it is not yet tied to a literature-owned "
                    "fixture or named issue-scope benchmark anchor."
                ),
            ),
        ),
        (
            "ascani_2022_electrolyte_lle",
            LiteratureBenchmarkEntry(
                case="ascani_2022_electrolyte_lle",
                title="Ascani 2022 electrolyte LLE benchmark",
                classification="blocker_requires_followup",
                coverage_kind="inventory_only",
                package_surface=(),
                validation_paths=(),
                notes=(
                    "The package has generic electrolyte LLE capability and Khudaida coverage, but no package-owned "
                    "Ascani 2022 literature benchmark fixture or runnable test entry."
                ),
            ),
        ),
        (
            "ascani_2023_reactive_lle",
            LiteratureBenchmarkEntry(
                case="ascani_2023_reactive_lle",
                title="Ascani 2023 reactive LLE benchmark",
                classification="blocker_requires_followup",
                coverage_kind="inventory_only",
                package_surface=(),
                validation_paths=(),
                notes=(
                    "Reactive LLE benchmark coverage is still blocked on broader generic reactive LLE production "
                    "support and a package-owned literature fixture."
                ),
            ),
        ),
        (
            "khudaida_2026_salting_out_lle",
            LiteratureBenchmarkEntry(
                case="khudaida_2026_salting_out_lle",
                title="Khudaida salting-out LLE benchmark",
                classification="already_supported_with_tests",
                coverage_kind="suite_smoke_and_opt_in",
                package_surface=(
                    "src/epcsaft/equilibrium_core/confidence.py",
                    "tests/workflows/validation/equilibrium_core/test_electrolyte_lle_confidence.py",
                ),
                validation_paths=(
                    "tests/workflows/validation/equilibrium_core/test_electrolyte_lle_confidence.py",
                ),
                notes=(
                    "The package already ships a dedicated Khudaida confidence suite with smoke and opt-in full modes."
                ),
            ),
        ),
        (
            "hubach_yu_lithium_equilibrium",
            LiteratureBenchmarkEntry(
                case="hubach_yu_lithium_equilibrium",
                title="Hubach/Yu lithium-related equilibrium benchmark",
                classification="already_supported_with_tests",
                coverage_kind="fixture_and_opt_in_hard_case",
                package_surface=("tests/equilibrium/electrolyte/test_hubach_electrolyte_lle.py",),
                validation_paths=("tests/equilibrium/electrolyte/test_hubach_electrolyte_lle.py",),
                notes=(
                    "Fixture, seed, and diagnostics coverage exists now; the distinct native split solve remains an "
                    "opt-in hard-case regression."
                ),
            ),
        ),
    )
)


def run_literature_benchmarks(*, case: str | None = None) -> dict[str, Any]:
    if case is not None and case not in LITERATURE_CASES:
        raise ValueError(f"Unknown case {case!r}. Expected one of: {', '.join(LITERATURE_CASES)}")

    selected = [case] if case is not None else list(LITERATURE_CASES)
    entries = [LITERATURE_CASES[name].to_payload() for name in selected]
    classification_counts = Counter(entry["classification"] for entry in entries)
    if set(classification_counts) - ALLOWED_CLASSIFICATIONS:
        unknown = ", ".join(sorted(set(classification_counts) - ALLOWED_CLASSIFICATIONS))
        raise ValueError(f"Unsupported literature benchmark classification(s): {unknown}")
    return {
        "issue": 95,
        "title": "Literature benchmark suite",
        "package_version": str(__version__),
        "git_commit": _git_commit(),
        "selected_cases": selected,
        "classification_counts": dict(classification_counts),
        "cases": entries,
    }


def render_literature_benchmark_table(payload: dict[str, Any]) -> str:
    headers = ["case", "classification", "coverage_kind", "validation_paths"]
    widths = {header: len(header) for header in headers}
    rows: list[list[str]] = []
    for row in payload["cases"]:
        values = [
            str(row["case"]),
            str(row["classification"]),
            str(row["coverage_kind"]),
            str(len(row["validation_paths"])),
        ]
        rows.append(values)
        for header, value in zip(headers, values):
            widths[header] = max(widths[header], len(value))

    def _format(values: list[str]) -> str:
        return "  ".join(value.ljust(widths[header]) for header, value in zip(headers, values))

    header_line = _format(headers)
    divider = "  ".join("-" * widths[header] for header in headers)
    body = "\n".join(_format(values) for values in rows)
    return "\n".join((header_line, divider, body))


def payload_as_json(payload: dict[str, Any]) -> str:
    return json.dumps(payload, indent=2, sort_keys=False)
