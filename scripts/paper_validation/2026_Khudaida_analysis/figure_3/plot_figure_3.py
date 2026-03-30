from __future__ import annotations

from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

import _common as common


def main() -> None:
    common.plot_lle_figure(Path(__file__).resolve().parent, 3, 303.15, 0.05)


if __name__ == "__main__":
    main()
