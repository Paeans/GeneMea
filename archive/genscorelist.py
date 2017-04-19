##
# Author:       G. Pan
# Create Date:  2016 - 03 - 27
# Version:      1.0
#
##
# genscorelist.py
# define method to extract sub gene sequence
# and method to calculate scores from sequence
# output files
# sub gene sequence file is store under:
# GeneData/SEQ/chr?.se
# *****
# peakpoint subsequence
# ... ...
# *****
# default subsequence length is 2000
# +/- 1000 bp region of every peak
#
# scores file is store under:
# GeneData/SCORE/MOTIFNAME/chr?.sc
# *****
# score socre Non ... ...
# ... ...
# *****
# every line is a score list of the motif 
# responds to the subsequence of SEQ files
#
##

import sys
import os

from peakpoint import *
from fafile import *
from motif import *

def extractsubseq(peakname, size = 2000):
    fadir = "../GeneData/FA"
    peakfiledir = "../GeneData/GSE11431_RAW"
    seqdir = "../GeneData/SEQ"
    pseqdir = seqdir + "/" + peakname
    
    peakfilename = peakfiledir + "/" + peakname + ".txt"
    
    fafiles = {}
    seqfiles = {}
    
    os.mkdir(pseqdir)
    peakfile = peakpoint(peakfilename)    
    chrname, peak = peakfile.getnext()
    while not chrname == None:
        if not chrname in fafiles:
            filename = fadir + "/" + chrname + ".fa"
            fafiles[chrname] = fafile(filename)
        
        if not chrname in seqfiles:
            seqfilename = pseqdir + "/" + chrname + ".seq"
            seqfiles[chrname] = open(seqfilename, 'w')
        
        #print "index is: ", peakfile.getindex()
        #print chrname, peak
        subseq = fafiles[chrname].getsequence(peak, size)
        
        if subseq == None:
            #raw_input()
            pass
        else:
            seqfiles[chrname].write(str(peak) + '\t' + subseq + '\n')
        
        chrname, peak = peakfile.getnext()
    
    peakfile.close()
    for fadatafile in fafiles.values():
        fadatafile.close()
    for seqdatafile in seqfiles.values():
        seqdatafile.close()
    return seqfiles.keys()
    
def calcscorelist(motifilename, seqfilesname, peakname):
    scoredir = "../GeneData/SCORE"
    seqdir = "../GeneData/SEQ"
    pscoredir = scoredir + "/" + peakname
    pseqdir = seqdir + "/" + peakname 
    
    counter = 0
    
    motifs = motiflist(motifilename)
    motif = motifs.getnext()
    
    while (not motif == None):    
        if counter < int(sys.argv[2]):        
            motif = motifs.getnext()
            counter += 1
            continue
        if counter >= int(sys.argv[2]) + int(sys.argv[3]):    
            break
            
        subseqlen = motif.getmotiflen()
        os.makedirs(pscoredir + "/" + motif.getac())
        for seqfilename in seqfilesname:
            seqfilenamestr = pseqdir + "/" + seqfilename + ".seq"
            scorenamestr = pscoredir + "/" + motif.getac() + "/" \
                + seqfilename + ".sc"
            
            with open(seqfilenamestr, 'r') as seqfile, \
                open(scorenamestr, 'w') as scorefile:
                for line in seqfile:
                    peak, subseq = line.strip('\n').split()
                    for i in range(0, len(subseq) - subseqlen + 1):
                        #print motif.calcpwm(subseq[i:i + subseqlen]),
                        res = motif.calcpwm(subseq[i:i + subseqlen])
                        if not res == None:                            
                            scorefile.write(str(res) + '\t')
                        else:
                            print "Error on " + subseq[i:i + subseqlen]
                    scorefile.write('\n')
                    #raw_input()
        
        motif = motifs.getnext()
        counter += 1
    pass

