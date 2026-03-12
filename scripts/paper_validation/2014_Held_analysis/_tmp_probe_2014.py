import time
import sys
from pathlib import Path
sys.path.insert(0, str(Path('scripts/Held_2014_analysis').resolve()))
import _common as c

for salt in ('NaCl','KBr'):
    for strategy in ('2008','2014'):
        print('start',salt,strategy, flush=True)
        for m in (0.05,0.1,0.5,1.0,2.0,3.0,4.0):
            t0=time.time()
            y=c.calc_osmotic_curve(salt,[m],strategy)[0]
            print(' ',m,'->',round(float(y),5),'dt',round(time.time()-t0,3), flush=True)

