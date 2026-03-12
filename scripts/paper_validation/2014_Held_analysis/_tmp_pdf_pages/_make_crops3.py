import matplotlib.pyplot as plt
from pathlib import Path
base=Path('scripts/Held_2014_analysis/_tmp_pdf_pages')
img10=plt.imread(base/'p-10.png')
img11=plt.imread(base/'p-11.png')
img12=plt.imread(base/'p-12.png')
plt.imsave(base/'fig3a_axis_only.png', img10[170:1080, 240:1210])
plt.imsave(base/'fig3b_axis_only.png', img10[170:1080, 900:1600])
plt.imsave(base/'fig4a_axis_only.png', img11[200:940, 250:930])
plt.imsave(base/'fig4b_axis_only.png', img11[200:940, 860:1510])
plt.imsave(base/'fig5_axis_only.png', img11[1600:2340, 880:1530])
plt.imsave(base/'fig6a_axis_only.png', img12[1610:2370, 210:920])
plt.imsave(base/'fig6b_axis_only.png', img12[1610:2370, 880:1580])
print('saved3')

