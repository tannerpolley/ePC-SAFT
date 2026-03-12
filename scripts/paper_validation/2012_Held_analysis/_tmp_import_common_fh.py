import faulthandler
faulthandler.dump_traceback_later(15, repeat=False)
import sys
from pathlib import Path
sys.path.insert(0, str(Path('scripts/Held_2012_analysis').resolve()))
print('before import', flush=True)
import _common
print('after import', flush=True)

