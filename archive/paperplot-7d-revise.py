#!/usr/bin/env python
import matplotlib
matplotlib.use('Agg')


import numpy as np
import matplotlib.pyplot as plt

groupname = ['ESRRB', 'E2F1', 'TCFCP', 'NMYC', 'NANOG', 'KLF4', 'ZFX']
mtdname   = ['Our Method', 'CENTDIST', 'CORE_TF', 'CEAS']

index = np.arange(7)
width = 0.15
colorlist = ['red', 'yellow', 'blue', 'green']

aucdata = [\
        [0.8892, 0.8828, 0.8437, 0.9432, 0.9148, 0.9109, 0.9006], \
        [0.7869, 0.8761, 0.9072, 0.8889, 0.9699, 0.8550, 0.8758], \
        [0.6627, 0.8202, 0.6889, 0.8052, 0.9399, 0.7075, 0.8353], \
        [0.6203, 0.5789, 0.6333, 0.7255, 0.8510, 0.6883, 0.6327]]

fig = plt.figure(num=None, dpi=1200, facecolor='w')
#ax = plt.subplot(111)
reclist = []
for i in range(4):
    reclist.append(plt.bar(index + 0.2 + i * width, aucdata[i], width, color=colorlist[i]))

plt.ylabel('AUC')
#plt.title('Comparison of four methods\' best result(draft)')
plt.xticks(index + 0.5, groupname)
#plt.ylim([0.5, 0.98])
#yticklabels = plt.gca().get_yticklabels()
#yticklabels[len(yticklabels) - 1].set_color('darkgrey')
plt.legend([x[0] for x in reclist], mtdname, bbox_to_anchor=(0.5, -0.1), loc = 'upper center', frameon = True, ncol = 4)

#plt.show()
plt.savefig('GeneFindingResult.png', format = 'png', dpi = 300, bbox_inches = 'tight')
