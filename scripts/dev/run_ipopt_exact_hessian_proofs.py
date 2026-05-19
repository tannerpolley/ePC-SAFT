from __future__ import annotations

import argparse
import json
import os
import subprocess
import sys
from pathlib import Path
from typing import Any

import numpy as np

try:
    from scripts.dev.native_runtime_env import apply_native_runtime_env, apply_to_current_process
except ModuleNotFoundError:
    from native_runtime_env import apply_native_runtime_env, apply_to_current_process


REPO_ROOT = Path(__file__).resolve().parents[2]
DEFAULT_OUTPUT_DIR = REPO_ROOT / "build" / "generated" / "ipopt_exact_hessian_proofs"
SUMMARY_PREFIX = "IPOPT_PROOF_SUMMARY "
CASES = ("vle", "electrolyte_lle")


def _hydrocarbon_workbook_mixture():
    import epcsaft

    params = {
        "m": np.asarray([1.0, 1.6069, 2.0020]),
        "s": np.asarray([3.7039, 3.5206, 3.6184]),
        "e": np.asarray([150.03, 191.42, 208.11]),
        "k_ij": np.asarray(
            [
                [0.0, 3.0e-4, 1.15e-2],
                [3.0e-4, 0.0, 5.10e-3],
                [1.15e-2, 5.10e-3, 0.0],
            ],
            dtype=float,
        ),
    }
    return epcsaft.ePCSAFTMixture.from_params(params, species=["Methane", "Ethane", "Propane"])


def _ascani_electrolyte_mixture():
    import epcsaft

    aq = np.asarray([0.798324680201737, 0.016320352824141723, 0.09267748348706063, 0.09267748348706063])
    org = np.asarray([0.37006036048879404, 0.6214918588210971, 0.004223890345054407, 0.004223890345054407])
    beta_org = 0.613766575013417
    feed = ((1.0 - beta_org) * aq + beta_org * org).tolist()
    mix = epcsaft.ePCSAFTMixture.from_dataset("2022_Ascani", ["H2O", "Butanol", "Na+", "Cl-"], feed, 298.15)
    return mix, feed


def _summary_from_payload(case: str, payload: dict[str, Any]) -> dict[str, Any]:
    return {
        "case": case,
        "status": str(payload.get("status", "")),
        "solver_status": str(payload.get("solver_status", "")),
        "accepted": bool(payload.get("accepted", False)),
        "solver_accepted": bool(payload.get("solver_accepted", False)),
        "hessian_approximation": str(payload.get("hessian_approximation", "")),
        "hessian_backend": str(payload.get("hessian_backend", "")),
        "exact_hessian_available": bool(payload.get("exact_hessian_available", False)),
        "eval_h_calls": int(payload.get("eval_h_calls", 0)),
        "iteration_count": int(payload.get("iteration_count", 0)),
        "iteration_history_size": int(payload.get("iteration_history_size", 0)),
    }


def _run_vle_child() -> dict[str, Any]:
    from epcsaft import _core

    mix = _hydrocarbon_workbook_mixture()
    payload = _core._native_neutral_bubble_p_eos_route_result(
        mix._native,
        233.15,
        [0.1, 0.3, 0.6],
        200,
        1.0e-8,
        0.0,
        "auto",
        8,
        1.0e-8,
        1.0e-3,
        1.0e-8,
        1.0e-8,
        None,
        print_level=5,
    )
    return _summary_from_payload("vle", payload)


def _run_electrolyte_lle_child() -> dict[str, Any]:
    from epcsaft import _core

    mix, feed = _ascani_electrolyte_mixture()
    payload = _core._native_electrolyte_lle_eos_route_result(
        mix._native,
        298.15,
        1.0e5,
        feed,
        500,
        1.0e-8,
        0.0,
        "auto",
        8,
        1.0e-8,
        1.0e-8,
        1.0e-8,
        1.0e-7,
        1.0e-3,
        None,
        print_level=5,
    )
    return _summary_from_payload("electrolyte_lle", payload)


def _run_child(case: str) -> int:
    apply_to_current_process()
    if case == "vle":
        summary = _run_vle_child()
    elif case == "electrolyte_lle":
        summary = _run_electrolyte_lle_child()
    else:
        raise SystemExit(f"Unknown proof case: {case}")
    print(SUMMARY_PREFIX + json.dumps(summary, sort_keys=True))
    return 0


def _parse_summary(log_text: str, case: str) -> dict[str, Any]:
    for line in reversed(log_text.splitlines()):
        if line.startswith(SUMMARY_PREFIX):
            return json.loads(line[len(SUMMARY_PREFIX) :])
    raise RuntimeError(f"{case} proof did not emit a summary line.")


def _assert_proof_summary(summary: dict[str, Any]) -> None:
    case = summary["case"]
    if not summary["accepted"]:
        raise RuntimeError(f"{case} proof was not accepted: {summary}")
    if not summary["solver_accepted"]:
        raise RuntimeError(f"{case} proof solver was not accepted: {summary}")
    if summary["hessian_approximation"] != "exact":
        raise RuntimeError(f"{case} proof did not use exact Hessians: {summary}")
    if summary["hessian_backend"] == "limited-memory":
        raise RuntimeError(f"{case} proof silently used limited-memory Hessians: {summary}")
    if not summary["exact_hessian_available"]:
        raise RuntimeError(f"{case} proof did not report exact Hessian availability: {summary}")
    if summary["eval_h_calls"] <= 0:
        raise RuntimeError(f"{case} proof did not call eval_h: {summary}")


def _assert_log_contains_ipopt_output(case: str, log_text: str) -> None:
    required = ("Ipopt", "Number of Iterations", "EXIT:")
    missing = [token for token in required if token not in log_text]
    if missing:
        raise RuntimeError(f"{case} proof log is missing Ipopt output tokens: {missing}")


def _run_parent(cases: list[str], output_dir: Path) -> int:
    output_dir.mkdir(parents=True, exist_ok=True)
    env = os.environ.copy()
    apply_native_runtime_env(env)
    src_path = str(REPO_ROOT / "src")
    existing_pythonpath = env.get("PYTHONPATH", "")
    env["PYTHONPATH"] = src_path if not existing_pythonpath else os.pathsep.join([src_path, existing_pythonpath])
    summaries: list[dict[str, Any]] = []
    for case in cases:
        command = [sys.executable, str(Path(__file__).resolve()), "--case", case, "--child"]
        result = subprocess.run(command, cwd=REPO_ROOT, env=env, text=True, capture_output=True, check=False)
        log_text = result.stdout
        if result.stderr:
            log_text += "\n[stderr]\n" + result.stderr
        log_path = output_dir / f"{case}.log"
        log_path.write_text(log_text, encoding="utf-8")
        if result.returncode != 0:
            raise RuntimeError(f"{case} proof failed with exit {result.returncode}; see {log_path}")
        summary = _parse_summary(log_text, case)
        _assert_proof_summary(summary)
        _assert_log_contains_ipopt_output(case, log_text)
        summary["log_path"] = str(log_path)
        summaries.append(summary)
        print(f"{case}: accepted exact Hessian proof captured at {log_path}")
    summary_path = output_dir / "summary.json"
    summary_path.write_text(json.dumps({"proofs": summaries}, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(f"summary: {summary_path}")
    return 0


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="Capture representative native Ipopt exact-Hessian proof logs.")
    parser.add_argument("--case", choices=CASES, action="append", dest="cases")
    parser.add_argument("--child", action="store_true", help=argparse.SUPPRESS)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    args = parser.parse_args(argv)
    cases = args.cases or list(CASES)
    if args.child:
        if len(cases) != 1:
            raise SystemExit("--child requires exactly one --case")
        return _run_child(cases[0])
    return _run_parent(cases, args.output_dir)


if __name__ == "__main__":
    raise SystemExit(main())
