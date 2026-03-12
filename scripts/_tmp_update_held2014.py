import csv
from pathlib import Path

root = Path('data/pcsaft_parameters/held_2014')
pure_path = root / 'pure' / 'any_solvent.csv'
kij_path = root / 'mixed' / 'binary_interaction' / 'k_ij.csv'
lij_path = root / 'mixed' / 'binary_interaction' / 'l_ij.csv'
khb_path = root / 'mixed' / 'binary_interaction' / 'k_hb_ij.csv'

# --- update pure/any_solvent.csv ---
with pure_path.open('r', encoding='utf-8-sig', newline='') as f:
    rows = list(csv.DictReader(f))
fieldnames = [
    'component','m','s','e','e_assoc','vol_a','assoc_scheme','dipm','dip_num','z','dielc','d_born','f_solv','MW'
]
by_comp = {r['component']: r for r in rows}

def upsert(component, **vals):
    row = by_comp.get(component, {k: '' for k in fieldnames})
    row['component'] = component
    for k,v in vals.items():
        row[k] = '' if v is None else str(v)
    by_comp[component] = row

# normalize existing water expression
upsert('H2O', s='sigma=2.7927+left(10.11 * e^-0.01775 T / K-1.417 * e^-0.01146 T / Kright)', MW='0.01801528')

# Table 1 non-charged additions (Held 2014)
upsert('Butanol', m='2.7515', s='3.6139', e='259.59', e_assoc='2544.6', vol_a='0.0067', assoc_scheme='2B', z='0', dielc='17.51', MW='0.0741216')
upsert('Benzene', m='2.4653', s='3.6478', e='287.35', assoc_scheme='', z='0', dielc='2.28', MW='0.0781118')
upsert('Toluene', m='2.8149', s='3.7169', e='285.69', assoc_scheme='', z='0', dielc='2.38', MW='0.0921405')
upsert('Glycine', m='4.8507', s='2.3270', e='216.96', e_assoc='2598.1', vol_a='0.0393', assoc_scheme='2B', z='0', dielc='78.0', MW='0.075067')
upsert('Alanine', m='5.4647', s='2.5222', e='287.59', e_assoc='3176.6', vol_a='0.0819', assoc_scheme='2B', z='0', dielc='78.0', MW='0.089094')

# Table 2 strategy-2 ion additions (Held 2014)
ion_rows = {
    'H3O+': {'m':'1.0','s':'3.4654','e':'500.00','z':'1','dielc':'8.0','MW':'0.019023'},
    'Cs+': {'m':'1.0','s':'3.9246','e':'180.00','z':'1','dielc':'8.0','MW':'0.132905'},
    'NH4+': {'m':'1.0','s':'3.5740','e':'230.00','z':'1','dielc':'8.0','MW':'0.018038'},
    'Mg2+': {'m':'1.0','s':'3.1327','e':'1500.00','z':'2','dielc':'8.0','MW':'0.024305'},
    'Ca2+': {'m':'1.0','s':'3.2648','e':'1060.00','z':'2','dielc':'8.0','MW':'0.040078'},
    'Cu2+': {'m':'1.0','s':'3.8379','e':'1610.90','z':'2','dielc':'8.0','MW':'0.063546'},
    'Zn2+': {'m':'1.0','s':'2.9798','e':'1250.00','z':'2','dielc':'8.0','MW':'0.06538'},
    'F-': {'m':'1.0','s':'1.7712','e':'275.00','z':'-1','dielc':'8.0','MW':'0.018998'},
    'NO3-': {'m':'1.0','s':'3.2988','e':'130.00','z':'-1','dielc':'8.0','MW':'0.062005'},
    'OH-': {'m':'1.0','s':'2.0177','e':'650.00','z':'-1','dielc':'8.0','MW':'0.017007'},
    'ClO4-': {'m':'1.0','s':'4.0186','e':'104.16','z':'-1','dielc':'8.0','MW':'0.09945'},
    'HCO3-': {'m':'1.0','s':'2.9296','e':'70.00','z':'-1','dielc':'8.0','MW':'0.061017'},
    'H2PO4-': {'m':'1.0','s':'3.6505','e':'95.00','z':'-1','dielc':'8.0','MW':'0.096987'},
    'Fo-': {'m':'1.0','s':'3.3077','e':'190.00','z':'-1','dielc':'8.0','MW':'0.045017'},
    'Ac-': {'m':'1.0','s':'3.9328','e':'150.00','z':'-1','dielc':'8.0','MW':'0.059044'},
    'SO4^2-': {'m':'1.0','s':'2.6491','e':'80.00','z':'-2','dielc':'8.0','MW':'0.09606'},
    'S2O3^2-': {'m':'1.0','s':'3.1877','e':'80.00','z':'-2','dielc':'8.0','MW':'0.11213'},
    'CO3^2-': {'m':'1.0','s':'2.4422','e':'249.26','z':'-2','dielc':'8.0','MW':'0.06001'},
    'HPO4^2-': {'m':'1.0','s':'2.1621','e':'146.02','z':'-2','dielc':'8.0','MW':'0.09598'},
    'PO4^3-': {'m':'1.0','s':'2.5516','e':'310.00','z':'-3','dielc':'8.0','MW':'0.09497'},
}
for comp, vals in ion_rows.items():
    upsert(comp, **vals)

# preserve stable order: existing first, then appended new in defined sequence
existing_order = [r['component'] for r in rows]
for comp in list(ion_rows.keys()) + ['Butanol','Benzene','Toluene','Glycine','Alanine']:
    if comp not in existing_order:
        existing_order.append(comp)
# butanol etc should appear near solvents
for comp in ['Butanol','Benzene','Toluene','Glycine','Alanine']:
    if comp in existing_order:
        existing_order.remove(comp)
insert_after = existing_order.index('Ethanol') + 1 if 'Ethanol' in existing_order else 0
for idx, comp in enumerate(['Butanol','Benzene','Toluene','Glycine','Alanine']):
    existing_order.insert(insert_after + idx, comp)

with pure_path.open('w', encoding='utf-8', newline='') as f:
    writer = csv.DictWriter(f, fieldnames=fieldnames)
    writer.writeheader()
    for comp in existing_order:
        row = by_comp[comp]
        for key in fieldnames:
            row.setdefault(key, '')
        writer.writerow({k: row.get(k, '') for k in fieldnames})

# --- build expanded interaction matrices ---
components = existing_order
index = {c:i for i,c in enumerate(components)}

def new_matrix():
    return [['0' for _ in components] for _ in components]

kij = new_matrix()
lij = new_matrix()
khb = new_matrix()

# load old kij
with kij_path.open('r', encoding='utf-8-sig', newline='') as f:
    old = list(csv.DictReader(f))
old_cols = [c for c in old[0].keys() if c and c != 'component'] if old else []
for r in old:
    ri = r['component']
    if ri not in index:
        continue
    for cj in old_cols:
        if cj not in index:
            continue
        val = (r.get(cj) or '').strip()
        if val:
            kij[index[ri]][index[cj]] = val


def set_pair(mat, a, b, val):
    if a not in index or b not in index:
        return
    sval = str(val)
    mat[index[a]][index[b]] = sval
    mat[index[b]][index[a]] = sval

# table1 solvent-water entries
set_pair(kij, 'Butanol', 'H2O', '0.000294 T / K-0.102')
set_pair(kij, 'Benzene', 'H2O', '0.000607 T / K-0.155478')
set_pair(kij, 'Toluene', 'H2O', '0.000726 T / K-0.195057')
set_pair(kij, 'Glycine', 'H2O', '-0.0612')
set_pair(kij, 'Alanine', 'H2O', '0.000291 T / K-0.14796')
set_pair(lij, 'Butanol', 'H2O', '-0.0044')
set_pair(khb, 'Butanol', 'H2O', '0.026')

# table2 ion-water entries (strategy 2)
water_ion = {
    'H3O+':'0.25','Li+':'-0.25','Na+':'-0.007981 T / K+2.37999','K+':'-0.004012 T / K+1.3959',
    'Cs+':'0.081','NH4+':'0.064','Mg2+':'-0.25','Ca2+':'0.0041','Cu2+':'0.25','Zn2+':'-0.25',
    'Cl-':'-0.25','Br-':'-0.25','I-':'-0.25','NO3-':'0.098','OH-':'-0.25','ClO4-':'-0.25',
    'HCO3-':'0.00','H2PO4-':'0.25','Fo-':'-0.23','Ac-':'-0.23','SO4^2-':'0.25','S2O3^2-':'0.25',
    'CO3^2-':'-0.25','HPO4^2-':'0.25','PO4^3-':'-0.25','F-':''
}
for ion,val in water_ion.items():
    if val:
        set_pair(kij, ion, 'H2O', val)

# table3 cation-anion kij
cat_cols = ['H3O+','Li+','Na+','K+','Cs+','NH4+','Mg2+','Ca2+','Cu2+','Zn2+']
rows = {
    'F-':      [None,None,'0.665','1.000','1.000',None,None,None,None,None],
    'Cl-':     ['0.654','0.669','0.317','0.064','-0.417','-0.566','0.817','1.000','-0.216','-0.705'],
    'Br-':     ['0.645','0.591','0.290','-0.102','-0.670','-0.639','0.752','0.993',None,None],
    'I-':      ['0.497','0.002','0.018','-0.312','-1.000','-0.787','0.317','0.229',None,None],
    'ClO4-':   ['0.861','0.406','-0.118',None,None,'-1.000','0.122','0.674',None,None],
    'NO3-':    [None,'0.364','-0.300','-0.585','-0.855','-0.419','0.285','-0.101',None,None],
    'H2PO4-':  [None,None,'-0.071','0.018',None,'-1.000',None,None,None,None],
    'HPO4^2-': [None,None,'-1.000','1.000',None,'-0.556',None,None,None,None],
    'PO4^3-':  [None,None,'-1.000','1.000',None,None,None,None,None,None],
    'Ac-':     [None,'-0.998','0.246','1.000','0.785',None,None,None,None,None],
    'Fo-':     [None,None,'-0.370','0.340',None,None,None,None,None,None],
    'OH-':     [None,'-0.566','0.649','1.000','0.564',None,None,None,None,None],
    'SO4^2-':  [None,'-0.952','-1.000','1.000','-1.000','-1.000','-1.000','-0.908','-1.000','-0.446'],
    'CO3^2-':  [None,None,'-1.000','1.000',None,None,None,None,None,None],
    'HCO3-':   [None,None,'-0.514','-0.476',None,None,None,None,None,None],
    'S2O3^2-': [None,None,'0.184',None,None,None,None,None,None,None],
}
for an, vals in rows.items():
    for cat, val in zip(cat_cols, vals):
        if val is not None:
            set_pair(kij, an, cat, val)

# table4 ion-organic kij and lij
set_pair(kij, 'NH4+', 'Butanol', '0.29')
set_pair(kij, 'Cl-', 'Butanol', '0.22')
set_pair(lij, 'NH4+', 'Butanol', '0.140')
set_pair(lij, 'Cl-', 'Butanol', '0.245')
for ion, val in [('Na+','0.35'),('Cl-','0.35'),('Br-','0.15'),('I-','0.07'),('SO4^2-','1.00')]:
    set_pair(kij, ion, 'Benzene', val)
for ion, val in [('Na+','0.35'),('Cl-','0.35'),('Br-','0.15')]:
    set_pair(kij, ion, 'Toluene', val)

# keep diagonal zero
for i in range(len(components)):
    kij[i][i] = lij[i][i] = khb[i][i] = '0'


def write_matrix(path, mat):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open('w', encoding='utf-8', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(['component', *components])
        for comp, row in zip(components, mat):
            writer.writerow([comp, *row])

write_matrix(kij_path, kij)
write_matrix(lij_path, lij)
write_matrix(khb_path, khb)

print('Updated', pure_path)
print('Updated', kij_path)
print('Updated', lij_path)
print('Updated', khb_path)
print('Component count:', len(components))
