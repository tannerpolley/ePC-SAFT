from __future__ import annotations

import argparse
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO_ROOT))

import build_config  # noqa: E402

DIST_ROOT = REPO_ROOT / "dist"


def _run_command(cmd: list[str], *, cwd: Path = REPO_ROOT) -> None:
    print("Running:", " ".join(cmd), flush=True)
    subprocess.run(cmd, cwd=str(cwd), check=True)


def _clean_dist_artifacts() -> None:
    DIST_ROOT.mkdir(exist_ok=True)
    for pattern in ("epcsaft-*.whl", "epcsaft-*.tar.gz"):
        for artifact in DIST_ROOT.glob(pattern):
            artifact.unlink()


def _clean_build_artifacts() -> None:
    if not build_config.BUILD_ROOT.exists():
        return
    for child in build_config.BUILD_ROOT.iterdir():
        if child.is_dir() and child.name.startswith(("bdist.", "lib.", "temp.")):
            shutil.rmtree(child, ignore_errors=True)


def _cleanup_generated_cpp() -> None:
    generated_cpp = build_config.generated_cython_cpp()
    if generated_cpp.exists():
        generated_cpp.unlink()


def _newest_wheel() -> Path:
    wheels = sorted(DIST_ROOT.glob("epcsaft-*.whl"), key=lambda path: path.stat().st_mtime_ns)
    if not wheels:
        raise RuntimeError("wheel build completed but no epcsaft wheel was found in dist/")
    return wheels[-1]


def _smoke_check_wheel(wheel: Path) -> None:
    with tempfile.TemporaryDirectory(prefix="epcsaft-wheel-smoke-") as temp_root:
        temp_path = Path(temp_root)
        target = temp_path / "target"
        cwd = temp_path / "cwd"
        cwd.mkdir()
        _run_command([sys.executable, "-m", "pip", "install", "--no-deps", "--target", str(target), str(wheel)])
        smoke_code = f"""
import sys
sys.path.insert(0, {str(target)!r})

import numpy as np
import epcsaft
from epcsaft import ePCSAFTMixture

mixture = ePCSAFTMixture.from_params(
    {{"m": np.asarray([1.0]), "s": np.asarray([3.7039]), "e": np.asarray([150.03])}},
    species=["Methane"],
)
state = mixture.state(T=300.0, x=np.asarray([1.0]), rho=100.0)
print("wheel smoke ok", epcsaft.__file__, state.compressibility_factor())
"""
        _run_command([sys.executable, "-c", smoke_code], cwd=cwd)
        shutil.rmtree(target, ignore_errors=True)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--skip-smoke", action="store_true", help="Build artifacts without installing and importing the wheel.")
    args = parser.parse_args()

    try:
        _clean_build_artifacts()
        _clean_dist_artifacts()
        _run_command([sys.executable, "-m", "build", "--sdist", "--wheel", "--no-isolation"])
        wheel = _newest_wheel()
        print(f"Built wheel: {wheel}")
        if not args.skip_smoke:
            _smoke_check_wheel(wheel)
    finally:
        _cleanup_generated_cpp()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
