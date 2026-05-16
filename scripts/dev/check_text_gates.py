from __future__ import annotations

import subprocess
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]

TEXT_SUFFIXES = {
    ".bib",
    ".cfg",
    ".cmake",
    ".cpp",
    ".h",
    ".hpp",
    ".ini",
    ".md",
    ".py",
    ".pyi",
    ".rst",
    ".sh",
    ".toml",
    ".txt",
    ".yaml",
    ".yml",
    ".ps1",
}

EXCLUDED_PREFIXES = (
    "docs/papers/",
)
EXCLUDED_PATHS: set[str] = set()


def _tracked_files() -> list[Path]:
    result = subprocess.run(
        ["git", "ls-files"],
        cwd=REPO_ROOT,
        text=True,
        capture_output=True,
        check=True,
    )
    paths: list[Path] = []
    for line in result.stdout.splitlines():
        rel = line.strip().replace("\\", "/")
        if not rel:
            continue
        if rel in EXCLUDED_PATHS or any(rel.startswith(prefix) for prefix in EXCLUDED_PREFIXES):
            continue
        path = REPO_ROOT / rel
        if not path.exists():
            continue
        if path.suffix.lower() in TEXT_SUFFIXES:
            paths.append(path)
    return paths


def _blocked_terms() -> tuple[str, ...]:
    first = "backend" + "_" + "unavailable"
    second = "finite" + "_" + "difference"
    third = "finite" + "-" + "difference"
    fourth = "finite" + " " + "difference"
    fifth = "central" + "_" + "perturbation"
    sixth = "source" + "_" + "perturbation"
    seventh = "perturbation" + "_" + "derivative"
    eighth = "numerical" + "diff"
    ninth = "levenberg" + "_" + "marquardt"
    tenth = "manual" + " numeric " + "perturbation"
    eleventh = "numerical" + " perturbation"
    twelfth = "numeric" + " perturbation"
    thirteenth = "perturbation" + " jacobian"
    fourteenth = "perturbation" + " derivative"
    fifteenth = "perturbation" + "-" + "jacobian"
    sixteenth = "perturbation" + "-" + "derivative"
    seventeenth = "perturbation" + "-" + "derived"
    return (
        first,
        second,
        third,
        fourth,
        fifth,
        sixth,
        seventh,
        eighth,
        ninth,
        tenth,
        eleventh,
        twelfth,
        thirteenth,
        fourteenth,
        fifteenth,
        sixteenth,
        seventeenth,
    )


def _source_blocked_terms(rel: str) -> tuple[str, ...]:
    if not rel.startswith(("src/", "tests/", "scripts/")):
        return ()
    return (
        "fall" + "back",
        "return" + "_best" + "_effort",
        "best" + "_effort",
        "_solve" + "_electrolyte" + "_bubble" + "_native",
        "native" + "_log" + "_pressure" + "_bisection",
        "hessian" + "_backend",
        "multi" + "start",
        "line" + "_search",
        "soft" + "_start",
        "damp" + "ing",
        "native" + "_scalar" + "_binary" + "_activity",
        "binary" + "_log" + "_amounts",
    )


def main() -> int:
    terms = _blocked_terms()
    matches: list[str] = []
    for path in _tracked_files():
        rel = path.relative_to(REPO_ROOT).as_posix()
        text = path.read_text(encoding="utf-8", errors="ignore").lower()
        path_terms = terms + _source_blocked_terms(rel)
        for line_number, line in enumerate(text.splitlines(), start=1):
            if any(term in line for term in path_terms):
                matches.append(f"{rel}:{line_number}: blocked solver/derivative gate term")

    if matches:
        print("Strict solver derivative text gate failed:", file=sys.stderr)
        print("\n".join(matches), file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
