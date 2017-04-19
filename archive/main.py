#!/usr/bin/env python
##
# Author:       G. Pan
# Create Date:  2016 - 03 - 06
# Version:      1.0
#
##
# main.py
# 
#
##

import os
import sys

from genomesplit import *
from fafile import *
from motif import *
from peakpoint import *
from genscorelist import *
from centdist import *

scoredir = "../GeneData/SCORE"
seqdir = "../GeneData/SEQ"
stdir = "../GeneData/ST"
fadir = "../GeneData/FA"
peakfiledir = "../GeneData/GSE11431_RAW"

motifilename = "../GeneData/TRANSFAC-FINAL.txt"

peakfilelist = \
        ["GSM288345_ES_Nanog"]
        #, "GSM288346_ES_Oct4", "GSM288347_ES_Sox2", \
        #"GSM288348_ES_Smad1", "GSM288350_ES_Tcfcp2l1", "GSM288353_ES_Stat3", \
        #"GSM288349_ES_E2f1", "GSM288351_ES_Ctcf", "GSM288352_ES_Zfx", \
        #"GSM288354_ES_Klf4", "GSM288355_ES_Esrrb", "GSM288356_ES_c-Myc", \
        #"GSM288357_ES_n-Myc", "GSM288359_ES_p300", "GSM288360_ES_Suz12"]

def test():
    
    calcfreq("GSM288345_ES_Nanog", 25)    
    pass

def transpeakfile():
    for filename in peakfilelist:
        with open(filename + ".txt", 'r') as pf, \
            open("result.txt", 'w') as rf:
            for line in pf:
                tmp = line.split()[0].split(":")
                sep = tmp[1].split("-")
                rf.write( tmp[0] + "\t" + sep[0] + "\t" + sep[1] + "\n" )

if __name__ == "__main__":
    #print len(sys.argv)
    #print sys.argv
    '''
    if peak_split(filelist) == 1:
        print "ERROR"
    chrl = chrlist()
    res = chrl.getfile("abc")
   
    test = fafile("datafortest.txt")
    
    print test.indextoseek(1)
    print test.getsequence(20, 1)
    print test.getsequence(20, 50)
    print test.getsequence(20, 51)
    print test.getsequence(20, 97)
   
    print test.seektoindex(7)
    print test.seektoindex(56)
    print test.seektoindex(57)
    print test.seektoindex(58)
    test.close()
    '''
    
    '''
    test = motiflist("motifraw.txt")
    
    tmp = test.getnext()
    print tmp.getpo()
    
    print tmp.calcpwm("ACGTACGTACGT")
    print tmp.calcpwm("GGacaggTGttG")
    print tmp.calcpwm("TTGTACGTACGT")
    print tmp.calcpwm("tCGTACGTACGT")
    print tmp.calcpwm("ACGcATCTGCCA")
    print tmp.calcpwm("ACGcATCTGCCAG")
    print tmp.calcpwm("ACGcDTCTGCCA")
    '''
    
    '''
    while not tmp == None:
        print tmp.getac()
        print tmp.getpo()
        tmp = test.getnext()
    '''
    
    '''
    chrlist = range(1, 20) + ["X", "Y", "M"]
    chrfile = {}
    for chr in chrlist:
        filename = "chr" + str(chr) + ".fa"
        chrfile[chr] = fafile(filename)
    
    print chrfile
    '''    
    
    '''
    for name in peakfilelist:
        print "Current peakfile is: " + name
        ot = sys.argv[1]
        if ot == "e":
            extractsubseq(name)
        elif ot == "c":
            chrnamelist = getchrnamelist(name)
            calcscorelist(motifilename, chrnamelist, name)
        elif ot == "g":
            num = int(sys.argv[2])
            motifnamelist = motiflist(motifilename).getmotifnames()
            for i in range(int(sys.argv[3])):
                genfreqst(peakfilelist[num + i], motifnamelist)
    '''
    
    if len(sys.argv) == 1:
        print "No command parameter is provide, run test"
        print "*****************************************"
        test()
        sys.exit(0)
    
    ot = sys.argv[1]
    if ot == "g" or ot == "m":
        num = int(sys.argv[2])
        motifnamelist = motiflist(motifilename).getmotifnames()
        for i in range(int(sys.argv[3])):
            peakname = peakfilelist[num + i]
            if ot =="g":                
                genfreqst(peakname, motifnamelist)
                print "Finish genfreqst"
            else:
                print "Start calcfreq"
                calcfreq(peakname, motifnamelist, 25)
    else:
        for name in peakfilelist:
            print "Current peakfile is: " + name
            if ot == "e":
                extractsubseq(name)
            elif ot == "c":
                chrnamelist = getchrnamelist(name)
                calcscorelist(motifilename, chrnamelist, name)
            else:
                print "Command line parameter is not right, check command line"
    
    pass


