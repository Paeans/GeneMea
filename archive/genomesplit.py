###
# Author:       G. Pan
# Create Date:  2016 - 03 - 12
# Version:      2.0
#
##
# genomesplit.py
# input:  GSM*
# output: GSM*_chr?_peak_list
# log:    ch? counter; total counter
#
###
# Author:       G. Pan
# Create Date:  2016 - 02 - 02
# Version:      1.0
#
##
# peak_split.py
# input:  GSM*
# output: GSM*_chr?_peak_list
# log:    ch? counter; total counter
#
###

class chrlist:
    def __init__(self):
        self.chrdict = {}
        
    def getfile(self, chrname):
        if not chrname in self.chrdict:
            return None
        return self.chrdict[chrname]
        
    def openfile(self, chrname):
        if not chrname in self.chrdict:
            filename = open(chrname, "w")
            self.chrdict[chrname] = filename
        return self.chrdict[chrname]
        
    def closefiles():
        for key, value in self.chrdict:
            value.close()

def split_file(filename):
    counter = index = maxdis = 0
    for line in filename:
        counter += 1
        s = line.split()[0].split(":")
        if not len(s) == 2:
            return -1
        start = int(s[1].split("-")[0])
        end = int(s[1].split("-")[1])
        if (end - start) > maxdis:
            maxdis = end - start
            index = counter
        peakloc = ( start + end ) / 2
        print s[0],     peakloc
    return maxdis, index, counter

def peak_split(filelist):
    for filename in filelist:
        with open(filename) as gsmfile:
            maxdis, index, counter = split_file(gsmfile)
            #print maxdis, index, counter
    pass

def test():
    print "test in peak_split"
    pass

def getSegment(genomefile, peakpoint, length = 1000):
    '''
    from genome file, extract a gene string
    string length is 2 * 1000
    string is from (peakpoint - 1000) to (peakpoint + 1000)
    
    parameters:
        genomefile: file have gene seq data
        peakpoint:  the position of peakpoint
        length:     size of segment
        
    return:
        the segment extracted, midpoint is peakpoint
    '''
    startpoint = peakpoint - length     #not add the head of file
    genomefile.seek(startpoint)
    genestr = ""
    while(len(genestr) < length * 2):
        tmpstr = genomefile.readline().rstrip('\n')
        genestr = genestr + \
            tmpstr[0 : min(length * 2 - len(genestr), len(tmpstr))]
    return genestr


def saveSegment(filename, peakpoint, segment):
    '''
    current only used to assistent coding process
    save the segment extract from genome file 
    format of data of saved file is:
    every line of the file is a segment information
    start with "peakpoint" and blank space, 
    then is the segment string
    
    parameters:
        filename:   file object used to store the segment
        peakpoint:  midpoint of this segment
        segment:    substring extracted from genome file    
    '''
    pass



def calcScoreList(segment, motif):
    '''
    calculate the scores of the motif in this segment
    from the start position of the segment, 
    calculate scores of every of substring of segment, 
    with length of motif
    '''
    pass


def calcScore(fragment, motif):
    '''
    calculate the score of a fragment on a motif
    compare each point of fragment and the coresponding
    position in motif, then sum the values
    
    parameters:
        fragment:   a substring of segment, same length as motif
        motif:      matrix used to calculate the PWM score
        
    return:
        score value of this fragment, as a evaluation of same or not
    '''
    pass

