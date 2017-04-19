#!/usr/bin/env python

import os
import sys
import time

from scipy.stats import gamma

from fafile import *
from motif  import *

fadir = "../GeneData/FA"
seqdir = "../GeneData/SEQ"
scoredir = "../GeneData/SCORE"
peakfiledir = "../GeneData/GSE11431_RAW"

motifilename = "../GeneData/TRANSFAC-FINAL.txt"

peakfilelist = \
        ["GSM288345_ES_Nanog",  "GSM288346_ES_Oct4",        "GSM288347_ES_Sox2", \
         "GSM288348_ES_Smad1",  "GSM288350_ES_Tcfcp2l1",    "GSM288353_ES_Stat3", \
         "GSM288349_ES_E2f1",   "GSM288351_ES_Ctcf",        "GSM288352_ES_Zfx", \
         "GSM288354_ES_Klf4",   "GSM288355_ES_Esrrb",       "GSM288356_ES_c-Myc", \
         "GSM288357_ES_n-Myc",  "GSM288359_ES_p300",        "GSM288360_ES_Suz12"]

falist = \
        ["chr" + str(i) for i in range(1, 20)] \
            + ["chrM", "chrX", "chrY"]

gammascale = 1
a_shape = 1
gammadist = [gamma.pdf(i, a_shape, loc = 0, scale = gammascale) for i in range(1, 51)]

vscore_param_k = 1

# functions used to extract sequence from fa file
# ***********************************************

def getpeaklist(peakfilename):
    result = []
    with open(peakfiledir + "/" + peakfilename + ".txt", 'r') as peakfile:
        for line in peakfile:
            s = line.split()[0].split(":")
            result.append([s[0], \
                (int(s[1].split("-")[0]) + int(s[1].split("-")[1])) / 2])
    return result

def getseqlist(peakfilename, size = 2000):
    result = []
    peaklist = getpeaklist(peakfilename)
    for faname in falist:
        fa = fafile(fadir + "/" + faname + ".fa")
        for peak in peaklist:            
            if not peak[0] == faname: continue
            subseq = fa.getsequence(peak[1], size)
            if not subseq == None: result.append(subseq)
        fa.close()
    return result

def storeseqlist(peakfilename, size = 2000):
    result = getseqlist(peakfilename)
    with open(seqdir + "/" + peakfilename + ".seq", 'w') as seqfile:
        for line in result:
            seqfile.write(line + "\n")
    pass
    
def restoreseqlist(peakfilename):
    with open(seqdir + "/" + peakfilename + ".seq", 'r') as seqfile:
        if seqfile == None:
            print "No seqfile in seq dir. Exec storeseqlist first"
            return None
        return [line.strip('\n') for line in seqfile]

# ***********************************************

# functions used to calculate pwm score
# ***********************************************


def getmotiflist(pwmfilename, prefix = None):
    motifs = motiflist(pwmfilename)
    result = []
    tmp = motifs.getnext()
    while not tmp == None:
        if prefix == None or \
            tmp.getid().startswith(prefix):
            result.append(tmp)
        tmp = motifs.getnext()
    return result

def calcpwmscore(seq, motifpwm):
    if not len(seq) == len(motifpwm.getpo()):
        print "ERROR: " + seq + " and " + motifpwm.getac() + " not equal size"
        return None
        
    pwmscore = 0
    for i in range(len(seq)):
        b = seq[i]
        pi = motifpwm.getpo()[i]
        index = -1
        if b == 'A' or b == 'a':
            index = 0
        elif b == 'C' or b == 'c':
            index = 1
        elif b == 'G' or b == 'g':
            index = 2
        elif b == 'T' or b == 't':
            index = 3
        elif b == 'N' or b == 'n':
            #print subseq
            index = 4
        else:
            print "ERROR: " + seq + " has unrecog pairwise"
            return None
        if index == 4:
            pwmscore += math.log(1)
        elif pi[index] == 0:
            return None
        else:
            pwmscore += math.log(pi[index] * 4)
    return pwmscore

def calcscmatrix(peakfilename, motifpwm):
    seqlist = restoreseqlist(peakfilename)
    if seqlist == None: return None
    counter = 0
    number = [0 for i in range(2000)]
    result = [[] for i in range(2000)]
    for line in seqlist:
        polen = len(motifpwm.getpo())
        for i in range(0, len(line) - polen + 1):        
            tmp = calcpwmscore(line[i:i + polen], motifpwm)
            if tmp == None or tmp == 0: continue
            #print tmp,
            counter += 1
            number[i] += 1 
            result[i].append(tmp)  
    print counter
    return result, number, counter
    

# ***********************************************

# functions used to operate score file
# ***********************************************

def storescore(scorelist, peakname, motifname):
    if not os.path.exists(scoredir + "/" + peakname):
        os.makedirs(scoredir + "/" + peakname)
    scorename = scoredir + "/" + peakname + \
                        "/" + motifname + ".sc"
    with open(scorename, 'w') as scorefile:
        for res in scorelist:
            for i in res:
                scorefile.write(str(i) + '\t')
            scorefile.write(str(len(res)) + '\n')
    pass
    
def restorescore(peakname, motifname):
    scorelist = [ [] for i in range(2000)]
    statlist = [ 0 for i in range(2000)]
    scorename = scoredir + "/" + peakname + \
                        "/" + motifname + ".sc"
    if not os.path.exists(scorename): return None, None
    index = 0
    with open(scorename, 'r') as scorefile:
        for line in scorefile:
            res = line.strip('\n').split()
            for i in res[0:len(res) - 1]:
                scorelist[index].append(float(i))
            statlist[index] = int(res[len(res) - 1])
            index += 1
    return scorelist, statlist

# ***********************************************

# functions used to operate score file
# ***********************************************

def combinescorematrix(raw_score_matrix, bp = 20):
    bins = len(raw_score_matrix)/2/bp
    result = []
    for i in range(bins):
        tmp = []
        #from 50 + i) * bp to 50 + i + 1) * bp - 1
        #from 49 - i) * bp to 49 - i + 1) * bp - 1
        index1 = (50 + i) * bp
        index2 = (49 - i) * bp
        for j in range(bp):
            tmp += raw_score_matrix[index1 + j]
            tmp += raw_score_matrix[index2 + j]
        result.append(tmp)
    return result

def calcfscore_k(list1, list2):
    if len(list1) != len(list2):
        return None
    if len(list1) == 0: return 0.0
    maxscore = max(max(list1), max(list2))
    result = 0.0
    for i in range(len(list1)):
        result += (list1[i] - list2[i])/maxscore
    return result

def calcfscore_c(list1, list2):
    result = 0.0
    if len(list1) == 0 or len(list2) == 0: return result
    maxscore = max(max(list1), max(list2))
    result += (sum(list1) - sum(list2))/maxscore/len(list2)
    return result

def calcmotifscore_k(binscorematrix, topk = 100, includev = True):
    minlen = topk
    result = 0.0
    subscorelist = []
    for binscorelist in binscorematrix:
        binscorelist.sort(reverse = True)
        #minlen = min(minlen, len(binscorelist))
    for i in range(1, len(binscorematrix) - 1):
        list1 = binscorematrix[i - 1]
        list2 = binscorematrix[i]
        minlen = min(len(list1), len(list2), minlen)
        binscore = calcfscore_k(list1[0:minlen], list2[0:minlen])
        result += binscore * gammadist[i]
        subscorelist.append(binscore)
    if includev == True:
        for i in range(1, len(subscorelist)):
            result += vscore_param_k * math.exp(subscorelist[i - 1] - subscorelist[i])
    #return result, minlen
    return float("{0:.6f}".format(result)), minlen
    
def calcmotifscore_c(binscorematrix, cutoff = 0, includev = True):
    result = 0.0
    cutmatrix = []
    subscorelist = []
    for binscorelist in binscorematrix:
        tmp = []
        for score in binscorelist:
            if score > cutoff:
                tmp.append(score)
        cutmatrix.append(tmp)
    for i in range(1, len(cutmatrix) - 1):
        list1 = binscorematrix[i - 1]
        list2 = binscorematrix[i]
        subscore = calcfscore_c(list1, list2)
        result += subscore * gammadist[i]
        subscorelist.append(subscore)
    if includev == True:
        for i in range(1, len(subscorelist)):
            result += vscore_param_k * math.exp(subscorelist[i - 1] - subscorelist[i])
    #return result, cutoff
    return float("{0:.6f}".format(result)), cutoff

# ***********************************************


if __name__ == "__main__":    
    '''
    for i in range(15):
        print time.strftime('%X %x %Z')
        #result = getseqlist(peakfilelist[i])
        #print peakfilelist[i], len(result)
        #storeseqlist(peakfilelist[i])
        print peakfilelist[i], len(restoreseqlist(peakfilelist[i]))        
    print time.strftime('%X %x %Z')
    '''
    
    '''
    res = getmotiflist(motifilename)
    print res[0].getpo()
    print len(res)
    print len(getmotiflist(motifilename, "V$"))
    subseq = "tCcaatttccaaaacacaacgcatctggcaag\
aaaatcaagaaggaagaccaacgcatggatac\
ttggacaggtgttgcattcctccctagaatag\
ggaataaaatacccatggaaggagttgcagtgaca"
    subseqlen = 12
    for i in range(0, len(subseq) - subseqlen + 1):        
        tmp = calcpwmscore(subseq[i:i + subseqlen], res[0])
        if tmp == None: continue
        print tmp, subseq[i:i + subseqlen]
    
### res = getmotiflist(motifilename, "V$")
    
    startmotif = int(sys.argv[1])
    motifnum = int(sys.argv[2])
    for peakname in peakfilelist:
        for i in range(motifnum):
            if startmotif + i >= len(res): continue
            motifpwm = res[startmotif + i]
            print time.strftime('%X %x %Z'), "calc peak " + peakname + " with " + motifpwm.getid()
            result, statlist, counter = calcscmatrix(peakname, motifpwm)
            storescore(result, peakname, motifpwm.getac())    
    '''
    
    '''
    counter = 0
    for peakname in peakfilelist:
        for motifname in res:
            #print time.strftime('%X %x %Z'), counter
            restore, reslist = restorescore(peakname, motifname.getac())
            if reslist == None: print peakname, motifname.getac(), "Error no score file"
            counter += 1
    print counter
    '''
    
    
    print time.strftime('%X %x %Z')
    for peakname in [peakfilelist[0]]:
        result_k = []
        result_c = []
        for motifpwm in res:
            restore, reslist = restorescore(peakname, motifpwm.getac())
            result = combinescorematrix(restore)
            motifscore_k = calcmotifscore_k(result, includev = False)
            motifscore_c = calcmotifscore_c(result, includev = False)
            result_k.append((motifscore_k[0], motifscore_k[1], motifpwm.getac()))
            result_c.append((motifscore_c[0], motifscore_c[1], motifpwm.getac()))
        print sorted(result_k, key = lambda score: score[0])
        print sorted(result_c, key = lambda score: score[0])
    
    restore, reslist = restorescore(peakfilelist[0], "M00095")
    result = combinescorematrix(restore)
    motifscore_k = calcmotifscore_k(result, includev = False)
    motifscore_c = calcmotifscore_c(result, includev = False)
    print motifscore_k, calcmotifscore_k(result)
    
    print time.strftime('%X %x %Z')
    

