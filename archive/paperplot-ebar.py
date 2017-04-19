#!/usr/bin/env python
import matplotlib
matplotlib.use('Agg')


import numpy as np
import matplotlib.pyplot as plt

groupname = ['>20K', '5-20K', '>5k']
mtdname   = ['Our Method', 'CENTDIST', 'CORE_TF', 'CEAS']

index = np.arange(3)
width = 0.1
colorlist = ['red', 'yellow', 'blue', 'green']
aucdata = [\
         [0.8750, 0.9174, 0.8992], [0.8567, 0.8974, 0.8800], [0.7239, 0.8220, 0.7800], [0.6108, 0.7244, 0.6757]]
errdata = [\
         [0.0271, 0.0182, 0.0304], [0.0624, 0.0503, 0.0551], [0.0844, 0.0957, 0.0985], [0.0284, 0.0926, 0.0908]]

fig = plt.figure(num=None, dpi=1200, facecolor='w')
#ax = plt.subplot(111)
reclist = []
for i in range(4):
    reclist.append(plt.bar(index + 0.3 + i * width, aucdata[i], width, color=colorlist[i], yerr = errdata[i]))

plt.xlabel('Size Range')
plt.ylabel('AUC')
#plt.title('Comparison of four methods\' best result(draft)')
plt.xticks(index + 0.5, groupname)
plt.ylim([0.55, 0.98])
#yticklabels = plt.gca().get_yticklabels()
#yticklabels[len(yticklabels) - 1].set_color('darkgrey')
plt.legend([x[0] for x in reclist], mtdname, bbox_to_anchor=(0.5, -0.101), loc = 'upper center', frameon = True, ncol = 4)
#plt.show()
plt.savefig('GeneFindingResult-ebar.png', format = 'png', dpi = 300, bbox_inches = 'tight')
