##
# Author:       G. Pan
# Create Date:  2016 - 03 - 14
# Version:      1.0
#
##
# motif.py
# class of motif list
# 
# motif file format is:
# <<< begin from next line >>>
# AC   M00001
# XX
# ID   V$MYOD_01
# XX
# DT   19.10.1992 (created); ewi.
# DT   22.10.1997 (updated); dbo.
# CO   Copyright (C), Biobase GmbH.
# XX
# NA   MyoD
# XX
# DE   myoblast determination gene product
# XX
# BF   T00526 MyoD; Species: mouse, Mus musculus.
# XX
# PO      A      C      G      T
# 01      1      2      2      0      S
# 02      2      1      2      0      R
# ... ...
# XX
# BA   5 functional elements in 3 genes
# XX
# CC   no comment
# XX
# //
# ... ...
# <<< end of motif >>>
# 
# this class used to operate on motif file
# like get each motif matrix
# 
##

import math

class motiflist:
    def __init__(self, filename):
        self.motifs = []
        self.index = 0
        with open(filename, 'r') as motifile:
            for line in motifile:
                if line == None:
                    break
                line = line.lstrip().rstrip()
                if not line.startswith("AC"):
                    if not line == "":
                        print "ERROR in parsing motif file"
                        print line
                    continue
                self.motifs.append(motif(line.split()[1], motifile))
        self.size = len(self.motifs)
        pass
    
    def getnext(self):
        if self.index >= self.size:
            return None
        self.index += 1
        return self.motifs[self.index - 1]
    
    def getmotifsize(self):
        return self.size
    
    def reset(self):
        self.index = 0
        pass
        
    def getmotifnames(self):
        motifnames = []
        for motif in self.motifs:
            motifnames.append(motif.getac())
        return motifnames
        
    def getmotifid(self, motifname):
        for motif in self.motifs:
            if motif.getac() == motifname:
                return motif.getid()
        return None
    


class motif:    
    def __init__(self, AC, motifile):
        self.AC = AC
        tag = False
        for line in motifile:
            if line.startswith("//"):
                break
                
            if line.startswith("ID"):
                self.ID = line.split()[1]
                
            if line.startswith("XX"):
                tag = False
                
            if line.startswith("PO"):
                self.PO = []
                tag = True
                #r = line.split()
                #if r[1] != 'A' or r[2] != 'C' or r[3] !='G' or r[4] != 'T':
                #    print "error"
            elif tag == True:
                tmp = line.split()
                pm = []
                pi = 0
                for i in range(1, 5):                    
                    pi += float(tmp[i])
                for i in range(1, 5):
                    pm.append(float(tmp[i]) / pi)
                self.PO.append(pm)
        pass    
    
    def getac(self):
        return self.AC    
    
    def getid(self):
        return self.ID    
    
    def getpo(self):
        return self.PO
    
    def getmotiflen(self):
        return len(self.PO)
        
    def calcpwm(self, subseq):
        if not len(subseq) == len(self.PO):
            print "ERROR: " + subseq + " and " + self.AC + " not equal size"
            return None
        
        pwmscore = 0
        for i in range(len(subseq)):
            b = subseq[i]
            pi = self.PO[i]
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
                print "ERROR: " + subseq + " has unrecog pairwise"
                return None
            if index == 4:
                pwmscore += math.log(1 * 4)
            elif pi[index] == 0:
                return "NoN"
            else:
                pwmscore += math.log(pi[index] * 4)
        return pwmscore
