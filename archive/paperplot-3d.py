#!/usr/bin/env python
import matplotlib
matplotlib.use('Agg')


import numpy as np
import matplotlib.pyplot as plt

groupname = ['TCFCP', 'E2F1', 'ESRRB']
mtdname   = ['CEAS', 'CORE_TF', 'CENTDIST', 'Our Method']

index = np.arange(3)
width = 0.15
colorlist = ['green', 'blue', 'yellow', 'red']

aucdata = [\
        [0.6333, 0.5789, 0.6203], \
        [0.6889, 0.8202, 0.6627], \
        [0.9072, 0.8761, 0.7869], \
        [0.8437, 0.8828, 0.8892]]

fig = plt.figure(num=None, dpi=1200, facecolor='w')
#ax = plt.subplot(111)
reclist = []
for i in range(4):
    reclist.append(plt.barh(index + 0.2 + i * width, aucdata[i], width, color=colorlist[i]))

plt.xlabel('AUC')
#plt.title('Comparison of four methods\' best result(draft)')
plt.yticks(index + 0.5, groupname)
yticklabels = plt.gca().get_yticklabels()
#yticklabels[len(yticklabels) - 1].set_color('darkgrey')
plt.legend([x[0] for x in reclist[::-1]], mtdname[::-1], bbox_to_anchor=(0.5, -0.101), loc = 'upper center', frameon = True, ncol = 4)

#plt.show()
plt.savefig('GeneFindingResult-a3.eps', format = 'eps', dpi = 1200, bbox_inches = 'tight')
