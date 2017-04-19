#!/usr/bin/env python

import sys
import math
import numpy as np

filename = sys.argv[1]

result = {}
with open(filename, 'r') as msfile:
    for line in msfile:
        score, topk, name, ID, tag, FN, f, v = line.strip('\n').split()
        if not FN in result.keys():
            result[FN] = ([], [], [], tag)
        r = result[FN]
        r[0].append(float(score))
        r[1].append(float(f))
        r[2].append(float(v))
        if r[3] != tag: print "ERROR in ", line
#print result

with open(filename[:-3] + ".mt", 'w') as mtfile, \
     open(filename[:-3] + ".ft", 'w') as ftfile, \
     open(filename[:-3] + ".vt", 'w') as vtfile:
    filelist = [mtfile, ftfile, vtfile]
    for key in result.keys():
        for i in range(len(filelist)):
            f = filelist[i]
            f.write(key + "\t" + result[key][3])
            scores = result[key][i]
            f.write("\t" + "%.5f" % max(scores))
            f.write("\t" + "%.5f" % min(scores))
            f.write("\t" + "%.5f" % np.mean(scores))
            f.write("\t" + "%.5f" % np.std(scores))
            f.write("\t" + "%.5f" % np.median(scores) + "\n")
