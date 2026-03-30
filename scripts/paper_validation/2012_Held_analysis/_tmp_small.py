import time,sys
from pathlib import Path
sys.path.insert(0, str(Path('scripts/Held_2012_analysis').resolve()))
import _common as c
print('start', flush=True)
t0=time.time()
grid,y=c.mean_ionic_activity_curve('2012_Held','NaI','ethanol',{'ethanol':1.0},m_max=1.8,points=5)
print('done',len(grid),'elapsed',time.time()-t0,'first',y[0],'last',y[-1], flush=True)

