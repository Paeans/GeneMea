##
# Author:       G. Pan
# Create Date:  2016 - 03 - 27
# Version:      1.0
#
##
# centdist.py
# program to restore centdist method
# calculate the frequency Z-score and
# the velocity Z-score
# 
# 
##

import math
import time
import operator
import os
import sys

from motif import *
from peakpoint import *

#[[0 for i in range(100)] for j in range(3)]

def getbound(peakname, motifname):
    scoredir = "../GeneData/SCORE"
    chrnamelist = getchrnamelist(peakname)
    maxscore, minscore = -100, 100
    sclistlen = 0
    scorelistdir = scoredir + "/" + peakname + "/" + motifname
    for chrname in chrnamelist:
        scorelistfile = scorelistdir + "/" + chrname + ".sc"
        with open(scorelistfile, 'r') as scorefile:
            for line in scorefile:
                tmp = line.strip('\n').split()
                sclistlen = len(tmp)
                for item in tmp:
                    if not item == "NoN": 
                        maxscore = max(float(item), maxscore)
                        minscore = min(float(item), minscore)
    return sclistlen, maxscore,minscore
                
def genfreqst(peakname, motifnames, costeps = 100, binsize = 20):
    scoredir = "../GeneData/SCORE"
    stdir = "../GeneData/ST"
    logfilename = stdir + "/" + peakname
    os.makedirs(stdir + "/" + peakname)
    chrnamelist = getchrnamelist(peakname)
    for motif in motifnames:
        sclistlen, maxscore, minscore = getbound(peakname, motif)
        
        with open(logfilename + ".log", 'a') as logfile:
            loginfo =  "*********************\n" + \
                "motif name is: " + motif + "\n" + \
                "max score is:  " + str(maxscore) + "\n" + \
                "min score is:  " + str(minscore) + "\n" + \
                "\n"
            logfile.write(loginfo)
                    
        listsize = sclistlen / binsize
        if not sclistlen % binsize == 0:
            listsize += 1
        statlist = [[0 for i in range(listsize)] for j in range(costeps)]
        step = (maxscore - minscore) / (costeps - 1)
        
        scorelistdir = scoredir + "/" + peakname + "/" + motif
        for chrname in chrnamelist:
            if step <= 0:
                continue
            scorelistfile = scorelistdir + "/" + chrname + ".sc"
            with open(scorelistfile, 'r') as scorefile:
                for line in scorefile:
                    counter = 0
                    tmp = line.strip('\n').split()
                    while counter < len(tmp):
                        if not tmp[counter] == "NoN": # and float(tmp[counter]) > cut-off
                            scval = float(tmp[counter])
                            statlist[int((maxscore - scval) / step)][counter / binsize] += 1
                            
                        counter += 1        
        if not listsize == 100:
            print "ERROR the length of " + peakname + \
                " in " + motif + " is not 100 but " + str(len(statlist[0]))
           
        mid = listsize / 2
        if not listsize % 2 == 0:
            mid += 1
        tres = [[0 for i in range(mid)] for j in range(costeps)]
        i = 0
        while i + mid < listsize and mid - i - 1 >= 0:
            for j in range(costeps):
                tres[j][i] += statlist[j][i + mid]
                tres[j][i] += statlist[j][mid - i -1]
            i += 1
        #print tres
        print time.strftime('%X %x %Z') + ": Calc... " + motif + " at " + peakname
        
        peakstname = stdir + "/" + peakname + "/" + motif
        with open(peakstname + ".st", 'w') as stfile:
            #stfile.write(motif)
            for num in tres:
                for item in num:
                    stfile.write(str(item) + '\t')
                stfile.write('\n')
 
def calcfreqsc(stlist, maxd = 25):
    stlistf = [float(x) for x in stlist]
    fscore, d = None, 0
    for tmp in range(1, maxd + 1):
        dt, ds = sum(stlistf), sum(stlistf[0:tmp])
        
        if dt == 0:
            fscore, d = 0, tmp
            continue
        pt = tmp / float(len(stlist))
        zf = (ds - dt * pt) / math.sqrt( dt * pt * (1 - pt) )
        if fscore == None or zf >= fscore:
            fscore, d = zf, tmp
        pass
        
    return fscore, d

def caclvelosc(stlist, maxd = 25, step = 4):
    stlistf = [float(x) for x in stlist]
    fscore, d = calcfreqsc(stlistf, maxd)
    diflist = []
    for i in range(len(stlistf) - step):
        diflist.append(stlistf[i] - stlistf[i + step])
        
    gv, bv = 0, 0
    for item in diflist[0:d]:
        if item > 0:
            gv += item
        else:
            bv -= item
    for item in diflist[d:]:
        if item > 0:
            bv += item
        else:
            gv -= item 
    if (gv + bv) == 0:
        zscore = 0
    else:
        zscore = (gv - bv) / math.sqrt(gv + bv)
    return zscore + fscore, fscore, zscore, d

def getscore(ftdict):
    return ftdict[1][0]
    

def calcmtfreq(peakname, motifname, maxd = 25):
    stdir = "../GeneData/ST"
    ftdir = "../GeneData/FT"
    stfilename = stdir + "/" + peakname + "/" + motifname
    ftfilename = ftdir + "/" + peakname + "/" + motifname
    
    score, fscore, zscore, d = None, None, None, 0
    ftdict = {}
    with open(stfilename + ".st", 'r') as stfile, \
        open(ftfilename + ".ft", 'w') as ftfile:
        stlist = []
        counter = 0
        for line in stfile:
            tmp = line.strip('\n').split()
            stlen = len(stlist)
            if stlen == 0:
                stlist = [float(x) for x in tmp]
            elif not stlen == len(tmp):
                print "Error have unconscious lines in ",
                print motifname + " on " + peakname
                break
            else:
                stlist = [stlist[i] + float(tmp[i]) for i in range(stlen)]
            tmpres = caclvelosc(stlist, maxd)
            #if score == None or tmpres[0] > score:
            #    score, fscore, zscore, d = tmpres
            if fscore == None or tmpres[1] > fscore:
                score, fscore, zscore, d = tmpres
            ftfile.write(str(counter) + " " \
                + str(score) + " "\
                + str(fscore) + " "\
                + str(zscore) + " "\
                + str(d) + '\n')
            counter += 1
        return score, fscore, zscore, d
        
    
def calcfreq(peakname, motifnamelist, maxd = 25):
    ftdict = {}
    ftdir = "../GeneData/FT"
    sffilename = ftdir + "/" + peakname
    os.makedirs(ftdir + "/" + peakname)
    for motif in motifnamelist:
        score, fscore, zscore, d = calcmtfreq(peakname, motif, maxd)
        ftdict[motif] = (score, fscore, zscore, d)

    #sortft = sorted(ftdict.items(), key=operator.itemgetter(1), reverse = True)
    sortft = sorted(ftdict.items(), key=getscore, reverse = True)
    with open(sffilename + ".sf", 'w') as sffile:
        for item in sortft:
            sffile.write(str(item[0]) + "\t" \
                + str(item[1][0]) + "\t" \
                + str(item[1][1]) + "\t" \
                + str(item[1][2]) +"\t" \
                + str(item[1][3]) +"\n")
    pass

    
'''
# old method 
# no cutoff considered

def genfreqst(peakname, motifnames):
    
    scoredir = "../GeneData/SCORE"
    chrnamelist = getchrnamelist(peakname)
    for motif in motifnames:
        sclistlen, maxscore, minscore = getbound(peakname, motifnames)
        statlist = []
        maxscore = -100
        minscore = 100
        scorelistdir = scoredir + "/" + peakname + "/" + motif
        for chrname in chrnamelist:
            scorelistfile = scorelistdir + "/" + chrname + ".sc"
            with open(scorelistfile, 'r') as scorefile:
                for line in scorefile:
                    counter = 0
                    tmp = line.strip('\n').split()
                    while counter < len(tmp):
                        if len(statlist) <= counter:
                            statlist.append(0)
                        if not tmp[counter] == "NoN": # and float(tmp[counter]) > cut-off
                            statlist[counter] += 1
                            maxscore = max(float(tmp[counter]), maxscore)
                            minscore = min(float(tmp[counter]), minscore)
                        counter += 1        
        #print statlist
        
        freqlist = []
        i = 0
        while i < len(statlist):
            if i % 20 == 0:
                freqlist.append(0)
            freqlist[i/20] += statlist[i]
            i += 1    
        #print freqlist
        
        if not len(freqlist) == 100:
            print "ERROR the length of " + peakname + \
                " in " + motif + " is not 100 but " + str(len(freqlist))
        tres = []    
        mid = len(freqlist) / 2
        i = 0
        while i + mid < len(freqlist) and mid - i - 1 >= 0:
            tres.append(0)
            tres[i] += freqlist[i + mid]
            tres[i] += freqlist[mid - i -1]
            i += 1
        #print tres
        print time.strftime('%X %x %Z') + ": Calc... " + motif + " at " + peakname
        peakstname = scoredir + "/" + peakname + "/" + peakname
        with open(peakstname + ".st", 'a') as stfile:
            stfile.write(motif)
            for num in tres:
                stfile.write('\t' + str(num))
            stfile.write('\n')
        with open(peakstname + ".log", 'a') as logfile:
            loginfo =  "*********************\n" + \
                "motif name is: " + motif + "\n" + \
                "max score is:  " + str(maxscore) + "\n" + \
                "min score is:  " + str(minscore) + "\n" + \
                "\n"
            logfile.write(loginfo)
'''    
