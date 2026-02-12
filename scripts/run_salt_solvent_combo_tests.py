#!/usr/bin/env python
"""Run salt/solvent combo tests for available MIAC datasets."""

from __future__ import annotations

import argparse
import os
import re
import subprocess
import sys
from pathlib import Path


SUPPORTED_WATER_SALTS = {"LiCl", "NaCl", "KCl", "LiBr", "NaBr", "KBr"}
METHANOL_SALT_PATTERN = re.compile(r"^(Li|Na|K)(Cl|Br|I)$")


def discover_water_salts(root: Path) -> tuple[list[str], list[str]]:
    data_dirs = [
        root / "data" / "MIAC_m" / "water",
        root / "data" / "MIAC_m" / "csv",
    ]
    salts_set = set()
    for data_dir in data_dirs:
        if data_dir.exists():
            salts_set.update(path.stem for path in data_dir.glob("*.csv"))
    salts = sorted(salts_set)
    supported = [salt for salt in salts if salt in SUPPORTED_WATER_SALTS]
    unsupported = [salt for salt in salts if salt not in SUPPORTED_WATER_SALTS]
    return supported, unsupported


def discover_methanol_salts(root: Path) -> tuple[list[str], list[str]]:
    salts_found = set()
    unsupported_names = set()

    pure_dir = root / "data" / "MIAC_m" / "methanol"
    if pure_dir.exists():
        for path in pure_dir.glob("methanol-*.csv"):
            salt = path.stem.replace("methanol-", "", 1)
            if salt:
                salts_found.add(salt)

    mixed_dir = root / "data" / "MIAC_m" / "water_methanol"
    if mixed_dir.exists():
        for path in mixed_dir.glob("*.csv"):
            stem = path.stem
            if stem.startswith("water-methanol-"):
                salt = stem.replace("water-methanol-", "", 1)
            else:
                salt = stem
            if salt:
                salts_found.add(salt)
            else:
                unsupported_names.add(stem)

    salts = sorted(salts_found)
    supported = [salt for salt in salts if METHANOL_SALT_PATTERN.match(salt)]
    unsupported = sorted(set(salts) - set(supported)) + sorted(unsupported_names)
    return supported, unsupported


def run_pytest_for_combo(root: Path, solvent: str, salt: str, python_exe: str, extra_pytest_args: list[str]) -> int:
    env = os.environ.copy()
    env["MIAC_SOLVENT"] = solvent
    env["MIAC_SALT"] = salt
    cmd = [
        python_exe,
        "-m",
        "pytest",
        str(root / "tests" / "test_DH_born_models.py"),
        "-s",
        "-k",
        "test_DH_born_models and not water_batch",
        *extra_pytest_args,
    ]
    print(f"\n=== Running combo: solvent={solvent}, salt={salt} ===")
    print(" ".join(cmd))
    completed = subprocess.run(cmd, cwd=str(root), env=env)
    return int(completed.returncode)


def run_command(root: Path, cmd: list[str], python_exe: str) -> int:
    cmd = [python_exe if token == "{python}" else token for token in cmd]
    print("\n=== Running command ===")
    print(" ".join(cmd))
    completed = subprocess.run(cmd, cwd=str(root))
    return int(completed.returncode)


def main() -> int:
    parser = argparse.ArgumentParser(description="Run available salt/solvent combo tests.")
    parser.add_argument(
        "--python",
        default=sys.executable,
        help="Python interpreter to use (default: current interpreter).",
    )
    parser.add_argument(
        "--include-sweep-tests",
        action="store_true",
        help="Also run rule-sweep tests (methanol fits and dielectric diff/fit).",
    )
    parser.add_argument(
        "--extra-pytest-arg",
        action="append",
        default=[],
        help="Extra arg passed to each combo pytest invocation (repeatable).",
    )
    args = parser.parse_args()

    root = Path(__file__).resolve().parents[1]
    python_exe = args.python

    water_supported, water_unsupported = discover_water_salts(root)
    methanol_supported, methanol_unsupported = discover_methanol_salts(root)

    print("Discovered water salts:", ", ".join(water_supported) if water_supported else "(none)")
    print("Discovered methanol salts:", ", ".join(methanol_supported) if methanol_supported else "(none)")
    if water_unsupported:
        print("Skipped unsupported water salts:", ", ".join(water_unsupported))
    if methanol_unsupported:
        print("Skipped unsupported methanol salts:", ", ".join(methanol_unsupported))

    failures: list[str] = []

    for salt in water_supported:
        rc = run_pytest_for_combo(root, "water", salt, python_exe, args.extra_pytest_arg)
        if rc != 0:
            failures.append(f"water/{salt}")

    for salt in methanol_supported:
        rc = run_pytest_for_combo(root, "methanol", salt, python_exe, args.extra_pytest_arg)
        if rc != 0:
            failures.append(f"methanol/{salt}")

    if args.include_sweep_tests:
        sweep_cmds = [
            ["{python}", "-m", "pytest", str(root / "tests" / "test_methanol_fits.py"), "-s"],
            ["{python}", "-m", "pytest", str(root / "tests" / "test_water_fits.py"), "-s"],
            ["{python}", "-m", "pytest", str(root / "tests" / "test_dielc_diff.py"), "-s"],
            ["{python}", "-m", "pytest", str(root / "tests" / "test_dielc_fit.py"), "-s"],
        ]
        for cmd in sweep_cmds:
            rc = run_command(root, cmd, python_exe)
            if rc != 0:
                failures.append(" ".join(token for token in cmd if token != "{python}"))

    print("\n=== Summary ===")
    if failures:
        print("Failures:")
        for item in failures:
            print(f"  - {item}")
        return 1

    print("All requested combo tests passed.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
