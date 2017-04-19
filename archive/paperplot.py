#!/usr/bin/env python
import matplotlib
matplotlib.use('Agg')


import numpy as np
import matplotlib.pyplot as plt

groupname = ['CMYC', 'E2F1', 'ESRRB', 'KLF4', 'NANOG', 'NMYC', 'OCT4', 'P300',	'SMAD1', 'SOX2', 'STAT3', 'TCFCP', 'ZFX', 'average']
mtdname   = ['CENTDIST', 'CORE_TF', 'CEAS', 'Our Method']

index = np.arange(14)
width = 0.15
colorlist = ['green', 'blue', 'yellow', 'red']

aucdata = [\
        #[0.9611, 0.8828, 0.8466, 0.8818, 0.8807, \
        # 0.9432, 0.8286, 0.8667, 0.8299, 0.8460, \
        # 0.8017, 0.8207, 0.8494, 0.8646], \
        [0.9957, 0.8761, 0.7869, 0.8550, 0.9699, \
         0.8889, 0.9300, 0.8646, 0.9507, 0.9507, \
         0.9175, 0.9072, 0.8758, 0.9053], \
        [0.9892, 0.8202, 0.6627, 0.7075, 0.9399, \
         0.8052, 0.9067, 0.9397, 0.9430, 0.9145, \
         0.8742, 0.6889, 0.8353, 0.8425], \
        [0.7828, 0.5789, 0.6111, 0.6883, 0.8510, \
         0.7255, 0.8650, 0.7917, 0.8531, 0.8684, \
         0.8067, 0.6333, 0.6327, 0.7369], \
        [0.9556, 0.8920, 0.8892, 0.9109, 0.9148, \
         0.9432, 0.8000, 0.8621, 0.8805, 0.8644, \
         0.8017, 0.8437, 0.9006, 0.8814]]

fig = plt.figure(num=None, dpi=1200, facecolor='w')
#ax = plt.subplot(111)
reclist = []
for i in range(4):
    reclist.append(plt.barh(index + 0.2 + i * width, aucdata[i], width, color=colorlist[i]))

plt.xlabel('AUC')
#plt.title('Comparison of four methods\' best result(draft)')
plt.yticks(index + 0.5, groupname)
yticklabels = plt.gca().get_yticklabels()
yticklabels[len(yticklabels) - 1].set_color('darkgrey')
plt.legend([x[0] for x in reclist[::-1]], mtdname[::-1], bbox_to_anchor=(0.5, -0.101), loc = 'upper center', frameon = True, ncol = 4)

#plt.show()
plt.savefig('GeneFindingDraft.eps', format = 'eps', dpi = 1200, bbox_inches = 'tight')
