from __future__ import annotations

import sys
from pathlib import Path

ANALYSIS_ROOT = Path(__file__).resolve().parents[3]
if str(ANALYSIS_ROOT) not in sys.path:
    sys.path.insert(0, str(ANALYSIS_ROOT))

from shared import figure_data


def main() -> None:
    path = figure_data.write_figure("figure_7")
    print(path.relative_to(figure_data.REPO_ROOT))


if __name__ == "__main__":
    main()
