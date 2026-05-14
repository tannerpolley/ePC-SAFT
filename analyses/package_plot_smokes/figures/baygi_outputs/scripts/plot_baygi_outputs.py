from __future__ import annotations

import subprocess
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[5]
TEST_PATH = REPO_ROOT / "analyses/package_plot_smokes/tests/plots/test_2015_baygi_outputs.py"


def main() -> None:
    subprocess.run([sys.executable, "run_pytest.py", str(TEST_PATH.relative_to(REPO_ROOT)), "-q"], cwd=REPO_ROOT, check=True)


if __name__ == "__main__":
    main()
