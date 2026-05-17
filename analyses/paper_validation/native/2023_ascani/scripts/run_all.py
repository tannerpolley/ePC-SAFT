from __future__ import annotations

import json
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[5]
ANALYSIS_DIR = REPO_ROOT / "analyses" / "paper_validation" / "native" / "2023_ascani"
SOURCE_MD = (
    REPO_ROOT
    / "docs"
    / "papers"
    / "md"
    / "Ascani - 2023 - Simultaneous Predictions of Chemical and Phase Equilibria in Systems with an Esterif.md"
)
RESULTS_DIR = ANALYSIS_DIR / "results" / "reactive_phase_equilibrium"
SUMMARY_JSON = RESULTS_DIR / "summary.json"


def _rel(path: Path) -> str:
    return path.relative_to(REPO_ROOT).as_posix()


def main() -> int:
    source_text = SOURCE_MD.read_text(encoding="utf-8", errors="replace")
    required_markers = ("Table 1. PC-SAFT pure-component parameters", "Table 2. Binary interaction parameters", "Table 4. Obtained")
    present_markers = [marker for marker in required_markers if marker in source_text]
    summary = {
        "status": "blocked",
        "lane": "ascani_2023_reactive_phase_equilibrium",
        "source_records": [_rel(SOURCE_MD)],
        "source_markers_present": present_markers,
        "source_facts": {
            "system_1": {
                "reaction": "acetic acid + 1-pentanol <=> pentyl acetate + water",
                "temperature_K": 318.15,
                "pressure_bar": 1.0,
                "K_a": 43.99,
            },
            "system_2": {
                "reaction": "acetic acid + 1-hexanol <=> hexyl acetate + water",
                "temperature_K": 353.6,
                "pressure_bar": 1.0,
                "K_a": 22.92,
            },
        },
        "blocker": {
            "kind": "missing_source_target_rows",
            "message": (
                "The local paper markdown provides model parameters and equilibrium constants, but this repo "
                "does not yet contain machine-readable source feed and phase-composition target rows for an "
                "Ascani 2023 reactive LLE gate."
            ),
        },
        "not_substituted": [
            "toy reactive LLE fixtures",
            "route construction tests",
            "monkeypatched native payloads",
            "non-literature random benchmarks",
        ],
    }
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    SUMMARY_JSON.write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(summary, indent=2))
    return 1


if __name__ == "__main__":
    raise SystemExit(main())
