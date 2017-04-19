##
# Author:       G. Pan
# Create Date:  2016 - 03 - 13
# Version:      1.0
#
##
# fafile.py
# class of fasta file
# 
# fasta file format is:
# <<< begin from next line >>>
# >chr*
# GTTAATGTAGCTTAATAACAAAGCAAAGCACTGAAAATGCTTAGATGGAT
# AATTGTATCCCATAAACACAAAGGTTTGGTCCTGGCCTTATAATTAATTA
# GAGGTAAAATTACACATGCAAACCTCCATAGACCGGTGTAAAATCCCTTA
# AACATTTACTTAAAATTTAAGGAGAGGGTATCAAGCACATTAAAATAGCT
# TAAGACACCTTGCCTAGCCACACCCCCACGGGACTCAGCAGTGATAAATA
# TTAAGCAATAAACGAAAGTTTGACTAAGTTATACCTCTTAGGGTTGGTAA
# ATTTCGTGCCAGCCACCGCGGTCATACGATTAACCCAAACTAATTATCTT
# CGGCGTAAAACGTGTCAACTATAAATAAATAAATAGAATTAAAATCCAAC
# TTATATGTGAAAATTCATTGTTAGGACCTAAACTCAATAACGAAAGTAAT 
# ... ...
# <<< end of fasta file >>>
# 
# this class used to operate on fasta file
# like extract gene sequence
# 
##

class fafile:
    def __init__(self, filename, sizeofline = 50):
        self.filename = filename
        self.file = open(filename, 'r')
        self.lsize = sizeofline
        
        self.header = self.file.readline().rstrip('\n')
        self.hsize = len(self.header) 
        pass
    
    def indextoseek(self, index):
        lines = ( index - 1 ) / 50
        return self.hsize + lines + index
    
    def seektoindex(self, seek):
        '''
        current not used in this program
        '''
        tmp = seek - self.hsize
        if tmp % (self.lsize + 1) == 0:
            return None
        return tmp - tmp / (self.lsize + 1)
    
    def getsequence(self, position = 0, size = 2000):
        seekpos = position - size / 2
        if seekpos < 0:
            return None
        #if seekpos != 0:
        self.file.seek(self.indextoseek(seekpos))
        #else:
        #    self.file.seek(0)
        
        genestr = ""
        while(len(genestr) < size):
            tmpstr = self.file.readline().rstrip('\n')
            if tmpstr == "":
                #print "WARNING: file " + self.filename + " position " + str(position) + " cannot be read anymore"
                return None
            remain = min(size - len(genestr), len(tmpstr))
            genestr = genestr + tmpstr[ 0 : remain ]
        return genestr
    
    def close(self):
        self.file.close()
        pass
    
     
