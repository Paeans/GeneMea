#!/usr/bin/env python

import os
import sys
import time
import math

from multiprocessing import Process, Queue

from scipy.stats import gamma
from Bio import motifs
from Bio.Seq import Seq

from fafile import *


pnum         = 20

fadir        = "../GeneData/FA"
seqdir       = "../GeneData/SEQ"
scoredir     = "../GeneData/SCORE"
msdir        = "../GeneData/MS"
peakfiledir  = "../GeneData/GSE11431_RAW"

motifilename = "../GeneData/TRANSFAC-FINAL.txt"

#"GSM288351_ES_Ctcf", "GSM288360_ES_Suz12"
peakfilelist = \
           ["GSM288345_ES_Nanog",  "GSM288346_ES_Oct4",  "GSM288347_ES_Sox2", \
            "GSM288348_ES_Smad1",  "GSM288349_ES_E2f1",  "GSM288350_ES_Tcfcp2l1", \
            "GSM288352_ES_Zfx",    "GSM288353_ES_Stat3", "GSM288354_ES_Klf4", \
            "GSM288355_ES_Esrrb",  "GSM288356_ES_c-Myc", \
            "GSM288357_ES_n-Myc",  "GSM288359_ES_p300"]

bindingfn    = \
           {"GSM288345_ES_Nanog" : ["NANOG", "OCT", "SOX", "ERE", "STAT"], \
            "GSM288346_ES_Oct4"  : ["NANOG", "OCT", "SOX", "ERE", "STAT", "CP2", "E2F", "EBOX"], \
            "GSM288347_ES_Sox2"  : ["NANOG", "OCT", "SOX", "ERE", "STAT", "CP2"], \
            "GSM288348_ES_Smad1" : ["NANOG", "OCT", "SOX", "STAT", "CP2", "ERE"], \
            "GSM288349_ES_E2f1"  : ["OCT", "STAT", "CP2", "E2F", "EBOX", "ZF5"], \
            "GSM288350_ES_Tcfcp2l1":["OCT", "SOX", "STAT", "CP2", "E2F"], \
            "GSM288352_ES_Zfx"   : ["OCT", "CP2", "E2F", "EBOX", "ZF5"], \
            "GSM288353_ES_Stat3" : ["NANOG", "OCT", "SOX", "ERE", "STAT", "CP2", "E2F", "EBOX"], \
            "GSM288354_ES_Klf4"  : ["NANOG", "OCT", "SOX", "ERE", "STAT", "E2F", "EBOX", "ZF5"], \
            "GSM288355_ES_Esrrb" : ["NANOG", "OCT", "SOX", "ERE", "STAT"], \
            "GSM288356_ES_c-Myc" : ["E2F", "EBOX", "ZF5"], \
            "GSM288357_ES_n-Myc" : ["OCT", "STAT", "E2F", "EBOX", "ZF5"], \
            "GSM288359_ES_p300"  : ["NANOG", "OCT", "SOX", "ERE", "STAT", "CP2"]}

falist       = \
           ["chr" + str(i) for i in range(1, 20)] \
            + ["chrM", "chrX", "chrY"]

fglist       = \
   {"V$MYOD_01":"EBOX", "V$E47_01":"EBOX", "V$VMYB_01":"VMYB", "V$CMYB_01":"MYB", \
    "V$AP4_01":"AP4", "V$MEF2_01":"MEF2", "V$ELK1_01":"ETS", "V$SP1_01":"SP1", "V$EVI1_06":"GATA_DIMER", \
    "V$ATF_01":"CREB", "V$HOX13_01":"HOX", "V$E2F_01":"E2F", "V$ELK1_02":"ETS", "V$RSRFC4_01":"MEF2", \
    "V$CETS1P54_01":"ETS", "V$P300_01":"P300", "V$P53_01":"P53", "V$VMAF_01":"VMAF", "V$VJUN_01":"CREB", \
    "V$NFE2_01":"AP1", "V$CREB_01":"CREB", "V$CREBP1_01":"CREB", "V$CREBP1CJUN_01":"CREB", "V$SOX5_01":"SOX", \
    "V$E4BP4_01":"E4BP4", "V$E2F_02":"E2F", "V$NFKAPPAB50_01":"NFKB", "V$NFKAPPAB65_01":"NFKB", \
    "V$CREL_01":"NFKB", "V$NFKAPPAB_01":"NFKB", "V$NMYC_01":"EBOX", "V$MYOGNF1_01":"MYOGNF1", "V$COMP1_01":"COMP1", \
    "V$HEN1_02":"HEN", "V$YY1_01":"CAAT", "V$IRF1_01":"IRF", "V$IRF2_01":"IRF", "V$TAL1BETAE47_01":"EBOX", \
    "V$TAL1ALPHAE47_01":"EBOX", "V$HEN1_01":"HEN", "V$YY1_02":"CAAT", "V$TAL1BETAITF2_01":"EBOX", "V$E47_02":"EBOX", \
    "V$CP2_01":"CP2", "V$DELTAEF1_01":"CACCT", "V$CETS1P54_02":"ETS", "V$GATA1_01":"GATA", "V$GATA2_01":"GATA", \
    "V$GATA3_01":"GATA", "V$EVI1_01":"GATA_DIMER", "V$EVI1_02":"GATA_DIMER", "V$EVI1_03":"GATA_DIMER", \
    "V$EVI1_04":"GATA_DIMER", "V$EVI1_05":"GATA_DIMER", "V$MZF1_01":"MZF1", "V$MZF1_02":"MZF1", "V$ZID_01":"ZID", \
    "V$IK1_01":"IK", "V$IK2_01":"IK", "V$IK3_01":"IK", "V$CDP_01":"ATCGAT", "V$PBX1_01":"PBX", "V$PAX6_01":"PAX", \
    "V$PAX2_01":"PAX", "V$S8_01":"S8", "V$CDXA_01":"TATA", "V$CDXA_02":"TATA", "V$CDP_02":"ATCGAT", "V$CLOX_01":"CLOX", \
    "V$CDPCR1_01":"ATCGAT", "V$CDPCR3_01":"ATCGAT", "V$CDPCR3HD_01":"ATCGAT", "V$E2_01":"E2", "V$NRF2_01":"NRF", \
    "V$CEBPB_01":"CEBP", "V$CREB_02":"CREB", "V$TAXCREB_01":"CREB", "V$TAXCREB_02":"CREB", "V$CEBPA_01":"CEBP", \
    "V$CEBPB_02":"CEBP", "V$MYCMAX_01":"EBOX", "V$MAX_01":"EBOX", "V$USF_01":"EBOX", "V$USF_02":"EBOX", "V$MYCMAX_02":"EBOX", \
    "V$PBX1_02":"PBX", "V$GATA1_02":"GATA", "V$GATA1_03":"GATA", "V$GATA1_04":"GATA", "V$HFH1_01":"FOX", "V$FOXD3_01":"FOX", \
    "V$HNF3B_01":"FOX", "V$HNF1_01":"HNF1", "V$TST1_01":"TST1", "V$HNF4_01":"ERE", "V$OCT1_01":"OCT", "V$OCT1_02":"OCT", \
    "V$OCT1_03":"OCT", "V$OCT1_04":"OCT", "V$AHR_01":"AHR", "V$LYF1_01":"LYF1", "V$PAX5_01":"PAX", "V$PAX5_02":"PAX", \
    "V$BRN2_01":"BRN2", "V$HSF1_01":"HSF", "V$HSF2_01":"HSF", "V$SRY_01":"FOX", "V$BRACH_01":"BRACH", "V$SRF_01":"SRF", \
    "V$ARP1_01":"ARP1", "V$RORA1_01":"ERE", "V$RORA2_01":"ERE", "V$COUP_01":"ERE", "V$CEBP_01":"CEBP", "V$SRY_02":"FOX", \
    "V$OCT1_05":"OCT", "V$OCT1_06":"OCT", "V$AP1FJ_Q2":"AP1", "V$AP1_Q2":"AP1", "V$AP1_Q6":"AP1", "V$AP4_Q5":"AP4", \
    "V$AP4_Q6":"AP4", "V$CREB_Q2":"CREB", "V$CREB_Q4":"CREB", "V$CREBP1_Q2":"CREB", "V$E2_Q6":"E2", "V$MYB_Q6":"MYB", \
    "V$MYOD_Q6":"EBOX", "V$GATA_C":"GATA", "V$GRE_C":"AR", "V$HNF1_C":"HNF1", "V$NFKB_C":"NFKB", "V$NFY_C":"CAAT", \
    "V$OCT_C":"OCT", "V$PADS_C":"TGTGGT", "V$POLY_C":"POLYC", "V$SEF1_C":"SEF1", "V$SRF_C":"SRF", "V$TATA_C":"TATA", \
    "V$USF_C":"EBOX", "V$SREBP1_01":"SREB/EBOX", "V$SREBP1_02":"SREB", "V$HAND1E47_01":"HAND1E47", "V$STAT_01":"STAT", \
    "V$STAT1_01":"STAT", "V$STAT3_01":"STAT", "V$VMYB_02":"VMYB", "V$VBP_01":"EBOX", "V$MEF2_02":"MEF2", "V$MEF2_03":"MEF2", \
    "V$MEF2_04":"MEF2", "V$AHRARNT_01":"AHR", "V$ARNT_01":"EBOX", "V$AHRARNT_02":"AHR", "V$BARBIE_01":"BARBIE", "V$T3R_01":"ERE", \
    "V$NKX25_01":"NKX", "V$NKX25_02":"NKX", "V$PPARA_01":"ERE", "V$EGR1_01":"EGR", "V$NGFIC_01":"EGR", "V$EGR3_01":"EGR", \
    "V$EGR2_01":"EGR", "V$OCT1_07":"OCT", "V$CHOP_01":"CHOP", "V$GFI1_01":"GFI", "V$XBP1_01":"EBOX", "V$TATA_01":"TATA", \
    "V$CAP_01":"CAP", "V$CAAT_01":"CAAT", "V$GC_01":"SP1", "V$NRSF_01":"CREB", "V$RREB1_01":"RREB", "V$ISRE_01":"IRF", \
    "V$HLF_01":"CREB", "V$OLF1_01":"OLF1", "V$STAF_01":"STAF", "V$STAF_02":"STAF", "V$XFD1_01":"FOX", "V$XFD2_01":"FOX", \
    "V$XFD3_01":"FOX", "V$AML1_01":"TGTGGT", "V$P53_02":"P53", "V$R_01":"R", "V$LMO2COM_01":"EBOX", "V$LMO2COM_02":"GATA", \
    "V$MIF1_01":"MIF1", "V$RFX1_01":"RFX", "V$RFX1_02":"RFX", "V$TCF11MAFG_01":"AP1", "V$TCF11_01":"AP1", "V$NFY_01":"CAAT", \
    "V$HFH3_01":"FOX", "V$FREAC2_01":"FOX", "V$FREAC3_01":"FOX", "V$FREAC4_01":"FOX", "V$FREAC7_01":"FOX", "V$HFH8_01":"FOX", \
    "V$NFAT_Q6":"HELIOS", "V$GATA1_05":"GATA", "V$GATA1_06":"GATA", "V$GATA2_02":"GATA", "V$GATA2_03":"GATA", "V$GATA3_02":"GATA", \
    "V$GATA3_03":"GATA", "V$PAX3_01":"PAX", "V$PAX4_01":"PAX", "V$PAX4_02":"PAX", "V$PAX4_03":"PAX", "V$PAX4_04":"PAX", \
    "V$MSX1_01":"MSX1", "V$HOXA3_01":"HOX", "V$EN1_01":"EN", "V$SOX9_B1":"SOX", "V$HNF4_01_B":"ERE", "V$AREB6_01":"CACCT", \
    "V$AREB6_02":"CACCT", "V$AREB6_03":"CACCT", "V$AREB6_04":"CACCT", "V$CART1_01":"CART1", "V$TGIF_01":"TGIF", "V$MEIS1_01":"MEIS1", \
    "V$MEIS1AHOXA9_01":"MEIS1", "V$MEIS1BHOXA9_02":"MEIS1", "V$FOXJ2_02":"FOX", "V$NKX61_01":"NKX", "V$HMX1_01":"HMX1", \
    "V$CHX10_01":"AT_RICH", "V$XVENT1_01":"XVENT1", "V$SPZ1_01":"SPZ", "V$ZIC1_01":"GLI", "V$ZIC2_01":"GLI", "V$ZIC3_01":"GLI", \
    "V$NKX3A_01":"NKX", "V$IRF7_01":"IRF", "V$MRF2_01":"MRF2", "V$FAC1_01":"FAC1", "V$STAT5A_01":"STAT", "V$STAT5B_01":"STAT", \
    "V$STAT5A_02":"STAT", "V$GATA6_01":"GATA", "V$POU3F2_01":"POU", "V$POU3F2_02":"POU", "V$POU6F1_01":"POU", "V$ROAZ_01":"ROAZ", \
    "V$AP2REP_01":"AP2", "V$AP2ALPHA_01":"AP2", "V$AP2GAMMA_01":"AP2", "V$TBP_01":"TATA", "V$FOXO4_01":"FOX", "V$FOXO1_01":"FOX", \
    "V$FOXO1_02":"FOX", "V$FOXO4_02":"FOX", "V$FOXO3_01":"FOX", "V$CDC5_01":"CDC5", "V$LUN1_01":"LUN1", "V$ATF6_01":"CREB", \
    "V$NCX_01":"NCX", "V$NKX22_01":"NKX", "V$PAX2_02":"PAX", "V$BACH2_01":"BACH", "V$MAZR_01":"SP1", "V$BACH1_01":"BACH", \
    "V$STAT1_03":"STAT", "V$STAT3_02":"STAT", "V$STAT4_01":"STAT", "V$STAT5A_04":"STAT", "V$STAT6_02":"STAT", "V$LHX3_01":"AT_RICH", \
    "V$PPARG_01":"ERE", "V$PPARG_02":"ERE", "V$E2F_03":"E2F", "V$AP1_01":"AP1", "V$GCNF_01":"ERE", "V$PPARG_03":"ERE", \
    "V$RP58_01":"RP58", "V$HTF_01":"EBOX", "V$ARNT_02":"EBOX", "V$MYCMAX_03":"EBOX", "V$GR_Q6":"AR", "V$NF1_Q6":"NF1", \
    "V$NFKB_Q6":"NFKB", "V$OCT1_Q6":"OCT", "V$SP1_Q6":"SP1", "V$AP1_C":"AP1", "V$CAAT_C":"CAAT", "V$CEBP_C":"CEBP", \
    "V$USF_Q6":"EBOX", "V$AP1_Q4":"AP1", "V$AP2_Q6":"AP2", "V$CEBP_Q2":"CEBP", "V$ER_Q6":"ERE"}



gammascale   = [1, 3]
gammashape   = [1, 2.5]
gammadist    = [gamma.pdf(i, gammashape[0], loc = 0, scale = gammascale[0]) for i in range(1, 51)]
gammalist    = []

for sc in gammascale:
    for sp in gammashape:
        gammalist.append([gamma.pdf(i, sp, loc = 0, scale = sc) for i in range(1, 51)])

fscore_w     = 1
vscore_w     = 1

fwlist       = [1, 10]
vwlist       = [1, 10]

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


def getmotiflist(pwmfilename, filetype = "TRANSFAC", prefix = None):
    with open(pwmfilename, 'r') as mhandle:
        motiflist = motifs.parse(mhandle, filetype)
    result = []
    for motifpwm in motiflist:
        if prefix == None or \
            motifpwm['ID'].strip().startswith(prefix):
            result.append(motifpwm)
    return result

def calcpssm(motifpwm, pseudocounts, cutoff, seqlist, bpsize, q):
    pssm = motifpwm.counts.normalize(pseudocounts).log_odds()
    alphabeta = motifpwm.alphabet
    pmax, pmin = pssm.max, pssm.min
    tmpnumber = [0 for i in range(bpsize)]
    tmpresult = [[] for i in range(bpsize)]
    for line in seqlist:
        tmp = pssm.calculate(Seq(line, alphabeta))
        for i in range(len(tmp)):
            if math.isinf(tmp[i]) or math.isnan(tmp[i]): continue
            score = tmp[i]
            if not pseudocounts == None:
                score = (tmp[i] - pmin) / (pmax - pmin)
            if score < cutoff: continue
            tmpnumber[i] += 1
            tmpresult[i].append(score)
    q.put(tmpnumber)
    q.put(tmpresult)


def calcscmatrix_mt(peakfilename, motifpwm, pseudocounts = None, cutoff = float("-inf"), bpsize = 2000):
    seqlist = restoreseqlist(peakfilename)
    if seqlist == None: return None
    counter = 0
    number = [0 for i in range(bpsize)]
    result = [[] for i in range(bpsize)]
    qlist = [Queue() for i in range(pnum)]
    
    totaltask = len(seqlist)
    step = totaltask / pnum + 1
    
    pro_list =[ Process(target=calcpssm, \
        args=(motifpwm, pseudocounts, cutoff, \
            seqlist[i * step: i * step + min(step, totaltask - i * step)], \
            bpsize, qlist[i], )) \
            for i in range(pnum) ]
    for prs in pro_list:
        prs.start()
    for q in qlist:
        tmpnum = q.get()
        for i in range(bpsize):
            number[i] += tmpnum[i]
        tmpres = q.get()
        for i in range(bpsize):
            result[i] += tmpres[i]
    for prs in pro_list:
        prs.join()
    return result, number, counter

def calcscmatrix(peakfile, motifpwm, pseudocounts = None, cutoff = float("-inf"), bpsize = 2000):
    #seqlist = restoreseqlist(peakfilename)
    seqlist = peakfile
    if seqlist == None: return None
    pssm = motifpwm.counts.normalize(pseudocounts).log_odds()
    pmax, pmin = pssm.max, pssm.min
    counter = 0
    number = [0 for i in range(bpsize)]
    result = [[] for i in range(bpsize)]
    for line in seqlist:
        tmp = pssm.calculate(Seq(line, motifpwm.alphabet))
        for i in range(len(tmp)):
            if math.isinf(tmp[i]) or math.isnan(tmp[i]): continue
            score = tmp[i]
            if not pseudocounts == None:
                score = (tmp[i] - pmin) / (pmax - pmin)
            if score < cutoff: continue
            counter += 1
            number[i] += 1
            result[i].append(score)
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

def calcmotifscore_k(binscorematrix, topk = 2000, includev = True):
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
        if binscore != 0.0:
            result += binscore * gammadist[i]
            #modify add divide by k(topk)
            #result += binscore * gammadist[i] / minlen
        subscorelist.append(binscore)
    if includev == True:
        for i in range(1, len(subscorelist)):
            result += vscore_w * math.exp(subscorelist[i - 1] - subscorelist[i])
    #return result, minlen
    print minlen, result
    return float("{0:.5f}".format(result)), minlen

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
    #return result, cutoff  "%3.5f" % result
    return float("{0:.5f}".format(result)), cutoff

def calcmotifscore(binscorematrix):
    maxlist = []
    listsize = 0
    for binscorelist in binscorematrix:
        binscorelist.sort(reverse = True)
        maxlist.append(binscorelist[0])
        listsize = max(len(binscorelist), listsize)
    for binscorelist in binscorematrix:
        if len(binscorelist) < listsize:
            binscorelist += [0 for i in range(listsize - len(binscorelist))]
    binsize = len(binscorematrix)
    #listsize = len(binscorematrix[0])
    for i in range(binsize - 1):
        binscorematrix[i][0] -= binscorematrix[i + 1][0]
        for j in range(1, listsize):
            binscorematrix[i][j] -= binscorematrix[i + 1][j]
            binscorematrix[i][j] += binscorematrix[i][j - 1]
    fscorelist = []
    vscorelist = []
    for gdist in gammalist:
        gfscore = []
        for k in range(listsize):
            tmp = 0
            for bindex in range(binsize - 1):
                tmp += gdist[bindex] * binscorematrix[bindex][k] / maxlist[bindex]
            gfscore.append(tmp)
        fscorelist.append(gfscore)
    for j in range(listsize):
        binscorematrix[0][j] = binscorematrix[0][j] / maxlist[0]
    for i in range(1, binsize - 1):
        for j in range(listsize):
            binscorematrix[i][j] = binscorematrix[i][j] / maxlist[i]
            binscorematrix[i - 1][j] -= binscorematrix[i][j]
    for k in range(listsize):
        vscore = 0
        for i in range(binsize - 2):
            vscore += 1/(1 + math.exp(-binscorematrix[i][k]/(k + 1)))
        vscorelist.append(vscore)
    
    result = []
    for fscore in fscorelist:
        maxfvalue = max(fscore)
        maxindex = fscore.index(maxfvalue)
        maxvvalue = vscorelist[maxindex]
        for fw in fwlist:
            for vw in vwlist:
                #result.append((maxindex, maxfvalue, maxvvalue, fw * maxfvalue + vw * maxvvalue))
                result.append((maxindex, float("{0:.5f}".format(maxfvalue)), \
                                float("{0:.5f}".format(maxvvalue)), \
                                float("{0:.5f}".format(fw * maxfvalue + vw * maxvvalue))))
    return result


# ***********************************************

#python ./gammadist.py 0 10 [0-14]
if __name__ == "__main__":
    '''
    every process execute a task to calculate
    motifs from argv[1] to argv[1] + argv[2]
    ChIP-Seq dataset is peakfilelist[argv[3]]
    '''
    
    startmotif = int(sys.argv[1])
    motifnum   = int(sys.argv[2])
    peakindex  = int(sys.argv[3])
    peakname   = peakfilelist[peakindex]
    result     = [[] for i in range(len(gammalist) * len(fwlist) * len(vwlist))]
    
    res = getmotiflist(motifilename, prefix = "V$")
    seqlist = restoreseqlist(peakname)
    for i in range(motifnum):
        if startmotif + i >= len(res):
            continue
        motifpwm = res[startmotif + i]
        
        label = 0
        if fglist[motifpwm['ID'].strip()] in bindingfn[peakname]: label = 1
        
        print time.strftime('%X'), "%04d" % i, "calc peak " + peakname + " with " + motifpwm['ID']
        #result, statlist, counter = calcscmatrix(peakname, motifpwm, pseudocounts = 0.01, cutoff = 0.0)
        scmatrix, statlist, counter = calcscmatrix(seqlist, motifpwm, pseudocounts = 0.01, cutoff = 0.0)
        scorematrix = combinescorematrix(scmatrix)
        '''
        topk = 10
        motifscore_k = calcmotifscore_k(scorematrix, topk = topk, includev = False)
        topk += 5
        while topk < len(seqlist):
            tmp_k = calcmotifscore_k(scorematrix, topk = topk, includev = False)
            if motifscore_k[0] < tmp_k[0]: motifscore_k = tmp_k
            topk += 5
        '''
        motifscore = calcmotifscore(scorematrix)
        for index in range(len(motifscore)):
            ms = motifscore[index]
            result[index].append((ms[3], ms[0], \
                            motifpwm['AC'].strip(), motifpwm['ID'].strip(), \
                            label, fglist[motifpwm['ID'].strip()], \
                            ms[1], ms[2]))
        '''
        result.append((motifscore_k[0], motifscore_k[1], \
                            motifpwm['AC'].strip(), motifpwm['ID'].strip(), \
                            label, fglist[motifpwm['ID'].strip()]))
        '''
        
    print time.strftime('%X'), "end of calc, begin save res"
    if not os.path.exists(msdir):
        os.makedirs(msdir)
    index = 0
    for sc in gammascale: 
        for sp in gammashape:
            for fw in fwlist: 
                for vw in vwlist:
                    msname = msdir + "/" + "%05.2f" % sc + "_" + \
                        "%05.2f" % sp + "_" + \
                        "%05.2f" % fw + "_" + \
                        "%05.2f" % vw + "_" + \
                        peakname[6:9] + "." + \
                        str(startmotif) + ".ms"
                    with open(msname, 'w') as wfile:
                        for rec in result[index]:
                            for re in rec:
                                wfile.write(str(re) + '\t')
                            wfile.write('\n')
                    index += 1
    print time.strftime('%X'), "end of process"
    # "%08.5f" % number # this will display number with total 8 digits and 5 precs
    '''
    msname = msdir + "/" + peakname + \
                   "." + str(startmotif) + ".ms"
    with open(msname, 'w') as wfile:
        for record in result_k:
            for re in record:
                wfile.write(str(re) + '\t')
            wfile.write('\n')
    '''
    

