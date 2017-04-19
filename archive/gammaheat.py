#!/usr/bin/env python
import matplotlib
matplotlib.use('Agg')

import os
import sys
import time
import math

from scipy.stats import gamma
from Bio import motifs
from Bio.Seq import Seq

from fafile import *

import matplotlib.pyplot as plt
import numpy as np

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



#gammascale   = [3.0 + x/10.0 for x in range(10)] #[0.5, 1, 1.5, 2, 2.5, 2.8, 3, 3.5]
#gammashape   = [2.0 + x/10.0 for x in range(10)] #[0.5, 1, 1.5, 2, 2.5, 2.8, 3, 3.5]
gammascale   = [3.5]
gammashape   = [3.6]
#fwlist       = [x/10.0 for x in range(1, 11)] #[1, 10]
#vwlist       = [x/10.0 for x in range(1, 11)] #[1, 10, 50, 100]
fwlist       = [0.1, 1]
vwlist       = [0.1, 100]

gammalist    = []


for sc in gammascale:
    for sp in gammashape:
        gammalist.append([gamma.pdf(i, sp, loc = 0, scale = sc) for i in range(1, 51)])



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

def calcmotifscore(binscorematrix):
    maxlist = []
    listsize = 0
    for binscorelist in binscorematrix:
        if len(binscorelist) == 0: continue
        binscorelist.sort(reverse = True)
        listsize = max(len(binscorelist), listsize)
    if listsize == 0: return None
    for binscorelist in binscorematrix:
        if len(binscorelist) < listsize:
            binscorelist += [0 for i in range(listsize - len(binscorelist))]
        if binscorelist[0] == 0:
            maxlist.append(0.01)
        else:
            maxlist.append(binscorelist[0])
    #print binscorematrix
    print listsize
    binsize = len(binscorematrix)
    for i in range(binsize - 1):
        binscorematrix[i][0] -= binscorematrix[i + 1][0]
        for j in range(1, 100000):
            binscorematrix[i][j] -= binscorematrix[i + 1][j]
    tmpmatrix = []
    nonsortmatrix = []
    for bins in binscorematrix[0:40]:
        nonsortmatrix.append(bins[0:100000])
    maxval = -10
    minval = 10
    for tmp in nonsortmatrix:
        tmp.sort(reverse = True)
        #maxval = max(max(tmp), maxval)
        #minval = min(min(tmp), minval)
    for t in nonsortmatrix:
        tmpmatrix.append([ sum(t[i*1000 : (i+1) * 1000]) for i in range(100)])
    for t in tmpmatrix:
        maxval = max(max(t), maxval)
        minval = min(min(t), minval)
    print tmpmatrix[1]
    print len(tmpmatrix)
    print maxval, minval
    xm = [[i for j in range(100)] for i in range(40)]
    ym = [range(100) for i in range(40)]
    #print xm, ym
    #print tmpmatrix
    plt.pcolormesh(xm, ym, tmpmatrix, vmin=-4, vmax=4)
    plt.title('V$OCT1_07')
    plt.xlabel('bins index')
    plt.ylabel('rows number(K)')
    plt.axis([0, 39, 0, 99])
    plt.colorbar()
    plt.savefig('gammaheat_oct.eps', format='eps', dpi = 1200)
    #for t in tmpmatrix:
    #   print sum(t)


# ***********************************************

#python ./gammadist.py 0 10 [0-14]
if __name__ == "__main__":
    '''
    every process execute a task to calculate
    motifs from argv[1] to argv[1] + argv[2]
    ChIP-Seq dataset is peakfilelist[argv[3]]
    '''
    
    startmotif = 0
    motifnum   = 1
    peakindex  = 10
    cutoff     = 0.0
    peakname   = peakfilelist[peakindex]
    result     = [[] for i in range(len(gammalist) * len(fwlist) * len(vwlist))]
    
    res = getmotiflist(motifilename, prefix = "V$OCT1_07")
    seqlist = restoreseqlist(peakname)
    for i in range(motifnum):
        if startmotif + i >= len(res):
            continue
        motifpwm = res[startmotif + i]
        
        label = 0
        if fglist[motifpwm['ID'].strip()] in bindingfn[peakname]: label = 1
        
        print time.strftime('%X'), "%04d" % i, "calc peak " + peakname + " with " + motifpwm['ID']
        #result, statlist, counter = calcscmatrix(peakname, motifpwm, pseudocounts = 0.01, cutoff = 0.0)
        # no cutoff
        #scmatrix, statlist, counter = calcscmatrix(seqlist, motifpwm, pseudocounts = 0.01, cutoff = 0.0)
        # with cutoff 0.8, use 0 pad each list
        scmatrix, statlist, counter = calcscmatrix(seqlist, motifpwm, pseudocounts = 0.01, cutoff = cutoff)
        scorematrix = combinescorematrix(scmatrix)
        #print scorematrix
        calcmotifscore(scorematrix)

