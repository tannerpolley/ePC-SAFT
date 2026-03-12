import sys
from pathlib import Path
sys.path.insert(0, str(Path('scripts/Held_2012_analysis').resolve()))
print('before import', flush=True)
import _common as c
print('after import', flush=True)

