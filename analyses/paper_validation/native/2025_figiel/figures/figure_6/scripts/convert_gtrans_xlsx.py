from __future__ import annotations

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
from scripts.plot_outputs import analysis_root
import sys



ROOT = Path(__file__).resolve().parent
ANALYSIS_ROOT = analysis_root(__file__)
if str(ANALYSIS_ROOT) not in sys.path:
    sys.path.insert(0, str(ANALYSIS_ROOT))

from shared import _common as common


def main() -> None:
    targets = [
        Path("data/reference/G_trans/water/methanol/K.xlsx"),
        Path("data/reference/G_trans/water/methanol/Br.xlsx"),
        Path("data/reference/G_trans/water/ethanol/Na.xlsx"),
        Path("data/reference/G_trans/water/ethanol/Cl.xlsx"),
    ]
    for rel in targets:
        xlsx = common.REPO_ROOT / rel
        csv_path = common.write_xlsx_to_csv(xlsx)
        print(csv_path)


if __name__ == "__main__":
    main()
