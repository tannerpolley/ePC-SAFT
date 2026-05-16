from __future__ import annotations

import argparse
import os
import shutil
import subprocess
import sys
import zipfile
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
DIST_ROOT = REPO_ROOT / "dist"
TEMPFILE_SITE = REPO_ROOT / "scripts" / "sandbox_tempfile_site"


def _run(cmd: list[str], *, env: dict[str, str] | None = None) -> None:
    print("Running:", " ".join(cmd), flush=True)
    subprocess.run(cmd, cwd=str(REPO_ROOT), env=env, check=True)


def _env(parallel: str | None = None) -> dict[str, str]:
    temp_root = REPO_ROOT / "build" / "dist-temp"
    temp_root.mkdir(parents=True, exist_ok=True)
    env = os.environ.copy()
    env["TMP"] = str(temp_root.resolve())
    env["TEMP"] = str(temp_root.resolve())
    env["TMPDIR"] = str(temp_root.resolve())
    env["EPCSAFT_SANDBOX_SAFE_TEMPFILE"] = "1"
    existing_pythonpath = env.get("PYTHONPATH")
    site_path = str(TEMPFILE_SITE.resolve())
    env["PYTHONPATH"] = site_path if not existing_pythonpath else os.pathsep.join([site_path, existing_pythonpath])
    if parallel:
        env["CMAKE_BUILD_PARALLEL_LEVEL"] = str(parallel)
    return env


def _clean_dist() -> None:
    DIST_ROOT.mkdir(exist_ok=True)
    for artifact in DIST_ROOT.glob("epcsaft-*"):
        artifact.unlink()


def _newest_wheel() -> Path:
    wheels = sorted(DIST_ROOT.glob("epcsaft-*.whl"), key=lambda path: path.stat().st_mtime_ns)
    if not wheels:
        raise RuntimeError("uv build completed but no epcsaft wheel was found in dist/")
    return wheels[-1]


def _assert_wheel_has_no_solver_dev_artifacts(wheel: Path) -> None:
    blocked_fragments = (
        "/include/ceres/",
        "/lib/cmake/ceres/",
        "/libceres",
        "numeric" + "_diff",
    )
    with zipfile.ZipFile(wheel) as archive:
        offenders = [
            name
            for name in archive.namelist()
            if any(fragment in f"/{name}".lower() for fragment in blocked_fragments)
        ]
    if offenders:
        sample = "\n".join(f"  - {name}" for name in offenders[:20])
        raise RuntimeError(
            "wheel contains Ceres development artifacts; vendored solver dependencies must not be installed into epcsaft:\n"
            f"{sample}"
        )


def _smoke_wheel(wheel: Path) -> None:
    target = REPO_ROOT / "build" / "wheel-smoke-target"
    shutil.rmtree(target, ignore_errors=True)
    target.mkdir(parents=True, exist_ok=True)
    _run(["uv", "pip", "install", "--target", str(target), str(wheel)])
    code = f"""
import sys
sys.path.insert(0, {str(target)!r})
import numpy as np
import epcsaft
import epcsaft._core
from epcsaft import ePCSAFTMixture
mixture = ePCSAFTMixture.from_params(
    {{"m": np.asarray([1.0]), "s": np.asarray([3.7039]), "e": np.asarray([150.03])}},
    species=["Methane"],
)
state = mixture.state(T=300.0, x=np.asarray([1.0]), rho=100.0)
print("wheel smoke ok", epcsaft.__file__, state.compressibility_factor())
"""
    _run([sys.executable, "-S", "-c", code])


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--skip-smoke", action="store_true")
    parser.add_argument(
        "--parallel",
        default=os.environ.get("EPCSAFT_BUILD_DIST_PARALLEL", "1"),
        help="CMake build parallelism for isolated PEP 517 builds. Defaults to 1 to avoid Windows Ceres OOMs.",
    )
    args = parser.parse_args()

    _clean_dist()
    _run(["uv", "build"], env=_env(args.parallel))
    wheel = _newest_wheel()
    _assert_wheel_has_no_solver_dev_artifacts(wheel)
    if not args.skip_smoke:
        _smoke_wheel(wheel)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
