from __future__ import annotations

from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

import _common as common


def main() -> None:
    common.plot_tables_9_10(Path(__file__).resolve().parent)


if __name__ == "__main__":
    main()
