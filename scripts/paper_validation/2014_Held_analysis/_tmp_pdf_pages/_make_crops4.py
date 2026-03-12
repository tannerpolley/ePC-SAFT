import matplotlib.pyplot as plt
from pathlib import Path
base=Path('scripts/Held_2014_analysis/_tmp_pdf_pages')
img10=plt.imread(base/'p-10.png')
img11=plt.imread(base/'p-11.png')
img12=plt.imread(base/'p-12.png')
plt.imsave(base/'fig3a_panel.png', img10[180:1060, 240:1050])
plt.imsave(base/'fig3b_panel.png', img10[180:1060, 1030:1760])
plt.imsave(base/'fig4a_panel.png', img11[220:930, 250:870])
plt.imsave(base/'fig4b_panel.png', img11[220:930, 880:1510])
plt.imsave(base/'fig5_panel_clean.png', img11[1550:2270, 900:1570])
plt.imsave(base/'fig6a_panel.png', img12[1600:2330, 210:900])
plt.imsave(base/'fig6b_panel.png', img12[1600:2330, 900:1590])
print('saved4')

