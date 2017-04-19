
#from genomesplit import *
#from fafile import *
#from motif import *
#from peakpoint import *
#from genscorelist import *
#from centdist import *


from scipy.stats import gamma
from multiprocessing import Process, Queue
import threading
from time import sleep

motifilename = "../GeneData/TRANSFAC-FINAL.txt"

fglist = {"V$MYOD_01":"EBOX", "V$E47_01":"EBOX", "V$VMYB_01":"VMYB", "V$CMYB_01":"MYB", \
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

#motiflist = motiflist(motifilename)
'''
with open("../GeneData/FT/GSM288345_ES_Nanog.sf", "r") as sfile, \
    open("../GeneData/FT/result.sf", "w") as wfile:
    for line in sfile:
        name, res = line.split()
        wfile.write(motiflist.getmotifid(name) + "\t" + line)
'''
#tmp = motiflist.getnext()
#print tmp.getpo()
'''
subseq = "tCcaatttccaaaacacaacgcatctggcaagaaaatcaagaaggaagaccaacgcatggatacttggacaggtgttgcattcctccctagaatagggaataaaatacccatggaaggagttgcagtgaca"
subseqlen = 12
for i in range(0, len(subseq) - subseqlen + 1):
    #print motif.calcpwm(subseq[i:i + subseqlen]),
    res = tmp.calcpwm(subseq[i:i + subseqlen])
    print res
    '''
'''    
for i in range(1, 51):
    print gamma.pdf(i, 1, loc = 0, scale = 1)
    print gamma.pdf(i, 1, loc = 0, scale = 10)
'''
'''
counter = 0
motifpwm = motiflist.getnext()
while motifpwm != None:
    if not motifpwm.getid().startswith("V$"): 
        motifpwm = motiflist.getnext()
        continue
    tmp = []
    for i in motifpwm.getpo():
        tmp.append(sum(i))
    
    tag = tmp[0]
    for t in tmp:
        if tag != t:
            print tmp, motifpwm.getac()
            counter += 1
            break
    
    motifpwm = motiflist.getnext()
print counter


with open(motifilename, 'r') as rfile:
    tag = False
    dotcount = 0
    necount = 0
    ac = ""
    name = ""
    po = []
    result = []
    for line in rfile:
        if line.startswith("AC"): 
            ac = line.strip('\n').split()[1]
            continue
        if line.startswith("ID"): 
            name = line.strip('\n').split()[1]
            continue
        if line.startswith("PO"): 
            tag = True
            po = []
            continue
        if line.startswith("XX"):
            if tag == True:
                dtag = False
                for tmp in po:                    
                    for t in tmp:
                        if "." in t: 
                            print ". in " + ac, name
                            if name.startswith("V$"): dotcount += 1
                            dtag = True
                            break
                    if dtag == True: break
                if dtag == True:
                    c = sum([float(x) for x in po[0]])
                    for p in po:
                        if c != sum([float(x) for x in p]):
                            if name.startswith("V$"): 
                                necount += 1
                                result.append(ac)
                            print "Not equal in " + ac, name
                            
                            #print po
                            break
                else:
                    c = sum([int(x) for x in po[0]])
                    for p in po:
                        if c != sum([int(x) for x in p]):
                            if name.startswith("V$"): 
                                necount += 1
                                result.append(ac)
                            print "Not equal in " + ac, name
                            
                            #print po
                            break
            tag = False
        if tag == True:
            p = line.strip('\n').split()
            po.append([p[1], p[2], p[3], p[4]])
    print "dotcount ", dotcount
    print "necount ", necount
    print len(result)
    print result
'''
#print fglist["V$SOX5_01"]
#print fglist["V$AHR_01"]

def f(i, q):
    q.put([i, 42, None, 'hello'])
    sleep(6)
    q.put([i, 41, None, 'hello'])
    pass

def g(i, q):
    q.put([i, 43, None, 'hello'])
    sleep(12)
    q.put([i, 44, None, 'hello'])
    
def print_number(q):
    while True:
        while q.qsize() != 0:
            m = q.get()
            if m == None: return
            print m
            
if __name__ == '__main__':
    q = Queue()
    p = [Process(target=f, args=(0, q,)), Process(target=g, args=(1, q,))]
    for n in p:
        n.start()
    t = threading.Thread(target=print_number, args=(q,))
    t.start()
    for n in p:
        n.join()
    q.put(None)
    t.join()
    print "FIN"
    

