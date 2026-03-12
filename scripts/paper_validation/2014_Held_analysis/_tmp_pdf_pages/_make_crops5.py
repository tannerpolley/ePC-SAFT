import matplotlib.pyplot as plt
from pathlib import Path
base=Path('scripts/Held_2014_analysis/_tmp_pdf_pages')
img10=plt.imread(base/'p-10.png')
img11=plt.imread(base/'p-11.png')
img12=plt.imread(base/'p-12.png')
plt.imsave(base/'fig3a_plot_only.png', img10[210:1030, 245:925])
plt.imsave(base/'fig3b_plot_only.png', img10[210:1030, 975:1605])
plt.imsave(base/'fig4a_plot_only.png', img11[245:910, 255:850])
plt.imsave(base/'fig4b_plot_only.png', img11[245:910, 885:1495])
plt.imsave(base/'fig5_plot_only.png', img11[1575:2225, 915:1545])
plt.imsave(base/'fig6a_plot_only.png', img12[1620:2295, 225:885])
plt.imsave(base/'fig6b_plot_only.png', img12[1620:2295, 915:1570])
print('saved5')

