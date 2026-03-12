import matplotlib.pyplot as plt
from pathlib import Path
base=Path('scripts/Held_2014_analysis/_tmp_pdf_pages')
img10=plt.imread(base/'p-10.png')
img11=plt.imread(base/'p-11.png')
img12=plt.imread(base/'p-12.png')
plt.imsave(base/'fig3a_axes.png', img10[170:1180, 210:930])
plt.imsave(base/'fig3b_axes.png', img10[170:1180, 920:1640])
plt.imsave(base/'fig4a_axes.png', img11[210:1020, 240:900])
plt.imsave(base/'fig4b_axes.png', img11[210:1020, 860:1520])
plt.imsave(base/'fig5_axes.png', img11[1620:2460, 860:1560])
plt.imsave(base/'fig6a_axes.png', img12[1620:2480, 190:930])
plt.imsave(base/'fig6b_axes.png', img12[1620:2480, 860:1600])
print('saved2')

