from __future__ import annotations

from pathlib import Path
import sys

ANALYSIS_SCRIPTS = Path(__file__).resolve().parents[3] / "scripts"
if str(ANALYSIS_SCRIPTS) not in sys.path:
    sys.path.insert(0, str(ANALYSIS_SCRIPTS))

import _common as common


def main() -> None:
    common.plot_lle_figure(Path(__file__).resolve().parent, 4, 313.15, 0.05)


if __name__ == "__main__":
    main()

