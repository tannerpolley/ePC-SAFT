from __future__ import annotations

from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parent
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

import _common as common


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
