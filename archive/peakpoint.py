##
# Author:       G. Pan
# Create Date:  2016 - 03 - 14
# Version:      1.0
#
##
# peakpoint.py
# class of peak point
# 
# file format is:
# <<< begin from next line >>>
# chr*:3053032-3053034	53	1	11.469303
# chr*:3333837-3333843	26	0	11.469303
# chr*:3334422-3334449	12	0	5.522257
# chr*:3473143-3473144	22	0	9.770147
# chr*:3671806-3671822	14	0	6.371835
# chr*:3937230-3937239	12	0	5.522257
# chr*:3985018-3985079	11	0	5.097468
# ... ...
# <<< end of fasta file >>>
# 
# this class used to operate on fasta file
# like extract gene sequence
# 
##

def getchrnamelist(peakname):
    peakfiledir = "../GeneData/GSE11431_RAW"    
    peakfilename = peakfiledir + "/" + peakname + ".txt"
    chrnamelist = []
    peakfile = peakpoint(peakfilename)    
    chrname, peak = peakfile.getnext()
    while not chrname == None:        
        if not chrname in chrnamelist:
            #seqfilename = seqdir + "/" + chrname + ".seq"
            chrnamelist.append(chrname)
        chrname, peak = peakfile.getnext()
    return chrnamelist

class peakpoint:
    def __init__(self, filename):
        self.file = open(filename, 'r')
        self.index = 0
        pass
    
    def getnext(self):
        '''
        return chr?, point
        '''
        self.index += 1
        line = self.file.readline()
        if line == None or line.strip() == "":
            return None, None
        s = line.split()[0].split(":")
        if not len(s) == 2:
            print "Error in line" + str(self.index)
            return None, None
        start, end = int(s[1].split("-")[0]), int(s[1].split("-")[1])
        peakloc = ( start + end ) / 2
        return s[0], peakloc
        
    def getindex(self):
        return self.index
        
    def close(self):
        self.file.close()
        pass
        
          
    
