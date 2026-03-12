import matplotlib.pyplot as plt
from pathlib import Path
base = Path('scripts/Held_2014_analysis/_tmp_pdf_pages')
out = base

def save_crop(page, name, x0,x1,y0,y1):
    img = plt.imread(base/page)
    crop = img[y0:y1, x0:x1]
    plt.imsave(out/name, crop)

save_crop('p-10.png','fig3_left_panel.png',320,910,310,1160)
save_crop('p-10.png','fig3_right_panel.png',930,1520,310,1160)
save_crop('p-11.png','fig4_left_panel.png',290,840,300,980)
save_crop('p-11.png','fig4_right_panel.png',860,1350,300,980)
save_crop('p-11.png','fig5_panel.png',930,1510,1710,2410)
save_crop('p-12.png','fig6_left_panel.png',220,860,1710,2450)
save_crop('p-12.png','fig6_right_panel.png',900,1530,1710,2450)
print('saved')

