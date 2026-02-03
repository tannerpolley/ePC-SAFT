import subprocess
import sys
from pathlib import Path


def main() -> int:
    repo_root = Path(__file__).resolve().parents[1]
    setup_py = repo_root / "setup.py"
    if not setup_py.exists():
        print(f"error: expected {setup_py} to exist")
        return 1

    cmd = [sys.executable, "setup.py", "build_ext", "--inplace", "--force"]
    print("Running:", " ".join(cmd))
    subprocess.run(cmd, cwd=str(repo_root), check=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
