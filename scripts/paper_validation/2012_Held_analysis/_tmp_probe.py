import time
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[3]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from scripts._env import require_pcsaft_install

require_pcsaft_install()

sys.path.insert(0, str(Path('scripts/Held_2012_analysis').resolve()))
import _common as c
from pcsaft import pcsaft_den, pcsaft_miac_m

print('loaded', flush=True)
for solvent in ('ethanol', 'methanol'):
    data = c.read_miac_dataset(Path(f'data/MIAC/{solvent}/{solvent}-NaI.csv'), solvent)
    print('solvent', solvent, 'rows', len(data), flush=True)
    mmax = max(float(r['molality']) for r in data)
    print('mmax', mmax, flush=True)
    comp = {solvent: 1.0}
    params = c.build_params('2012_Held', 'NaI', solvent, comp)
    species = c.species_for_combo('NaI', solvent)
    print('params ready', flush=True)
    for m in (0.01, 0.1, 0.5, 1.0, 1.6):
        t0 = time.time()
        x = c.molality_to_species_molefraction(m, 'NaI', solvent, comp)
        rho = pcsaft_den(c.T_REF, c.P_REF, x, params, phase='liq')
        vals = pcsaft_miac_m(c.T_REF, rho, x, params, species=species)
        key = c._resolve_pair_key(vals, 'NaI')
        print('m', m, 'gamma', vals[key], 'dt', round(time.time() - t0, 4), flush=True)

