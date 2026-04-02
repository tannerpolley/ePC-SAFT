import argparse
import subprocess
import sys
from pathlib import Path


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--skip-build",
        action="store_true",
        help="Run pytest without checking whether the editable install needs a rebuild.",
    )
    parser.add_argument(
        "--force-build",
        action="store_true",
        help="Force an editable reinstall before running pytest.",
    )
    args, pytest_args = parser.parse_known_args()

    repo_root = Path(__file__).resolve().parent
    build_script = repo_root / "scripts" / "build_pcsaft.py"
    if build_script.exists() and not args.skip_build:
        build_cmd = [sys.executable, str(build_script)]
        if args.force_build:
            build_cmd.append("--force")
        print("Running:", " ".join(build_cmd))
        subprocess.run(build_cmd, cwd=str(repo_root), check=True)
    elif not build_script.exists():
        print(f"warning: build script not found at {build_script}")

    cmd = [sys.executable, "-m", "pytest", "tests"]
    cmd.extend(pytest_args)
    print("Running:", " ".join(cmd))
    subprocess.run(cmd, cwd=str(repo_root), check=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
