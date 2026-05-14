from __future__ import annotations

import subprocess
import sys


from pathlib import Path
import sys as _bootstrap_sys
from pathlib import Path as _BootstrapPath

for _candidate in _BootstrapPath(__file__).resolve().parents:
    if (_candidate / "scripts" / "plot_outputs.py").is_file():
        if str(_candidate) not in _bootstrap_sys.path:
            _bootstrap_sys.path.insert(0, str(_candidate))
        break
else:
    raise ModuleNotFoundError("Could not locate repo root containing scripts/plot_outputs.py")
from scripts.plot_outputs import REPO_ROOT

TEST_PATH = REPO_ROOT / "analyses/package_validation/package_plot_smokes/tests/plots/test_equilibrium_plot_outputs.py"


def main() -> None:
    subprocess.run([sys.executable, "run_pytest.py", str(TEST_PATH.relative_to(REPO_ROOT)), "-q"], cwd=REPO_ROOT, check=True)


if __name__ == "__main__":
    main()
