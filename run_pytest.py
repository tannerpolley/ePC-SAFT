import subprocess
import sys
from pathlib import Path


def main() -> int:
    repo_root = Path(__file__).resolve().parent
    build_script = repo_root / "scripts" / "build_pcsaft.py"
    if build_script.exists():
        build_cmd = [sys.executable, str(build_script)]
        print("Running:", " ".join(build_cmd))
        subprocess.run(build_cmd, cwd=str(repo_root), check=True)
    else:
        print(f"warning: build script not found at {build_script}")

    cmd = [sys.executable, "-m", "pytest", "tests/test_cython.py"]
    if len(sys.argv) > 1:
        cmd.extend(sys.argv[1:])
    print("Running:", " ".join(cmd))
    subprocess.run(cmd, cwd=str(repo_root), check=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
