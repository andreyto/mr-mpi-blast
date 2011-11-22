#!/usr/bin/env python

import sys
from cStringIO import StringIO
import numpy as npy
 
class FastaReader(object):
    """Class that supports an input iterator protocol for a FASTA file.
    Example that prints an exact copy of the input file:
    for rec in FastaReader(open('seq.fsa','r')).records():
        print rec.header(),
        for line in rec.seqLines():
            print line,
    Instead of rec.seqLines(), you can use the methods which post-process the
    raw sequences lines: seqChunks(), seqArrays(), sequence().
    """

    def __init__(self,infile):
        if not hasattr(infile,"readline"):
            infile = openCompressed(infile,'r')
            self.ownInfile = True
        else:
            self.ownInfile = False
        self.infile = infile
        self.freshHdr = False
        self.maxLineLen = 0
        
    def records(self):
        infile = self.infile
        while True:
            if self.freshHdr:
                self.freshHdr = False
                yield self
                continue
            line = infile.readline()
            
            if not line:
                return
            # skip blank lines
            elif line.isspace():
                continue
            elif line.startswith(">"):
                self.hdr = line
                yield self
    
    def header(self):
        assert self.hdr.startswith('>')
        return self.hdr
        
    def getNCBI_Id(self):
        """Assume that header starts with '>gi|1234567|' and return the string id from second field."""
        return self.hdr.split('|',2)[1]
    
    def getNCBI_GI(self):
        return int(self.getNCBI_Id())

    def getSimpleId(self):
        """Assume that header starts with '>string_no_spaces ' and return that string."""
        return self.hdr.split(None,1)[0][1:]

    def seqLines(self):
        infile = self.infile
        while True:
            line = infile.readline()
            if not line:
                break
            elif line.isspace():
                continue
            elif line.startswith(">"):
                self.hdr = line
                self.freshHdr = True
                return
            self.maxLineLen = max(self.maxLineLen,len(line)-1)
            yield line

    def seqChunks(self,queryLen):
        seq = StringIO()
        for line in self.seqLines():
            seq.write(line.rstrip("\n"))
            if seq.tell() >= queryLen:
                yield seq.getvalue()
                seq.close()
                seq = StringIO()
        if seq.tell() > 0:
            yield seq.getvalue()
        seq.close()

    def seqArrays(self,queryLen):
        for s in self.seqChunks(queryLen):
            yield numpy.fromstring(s,dtype='S1')

    def sequence(self,format='str'):
        seq = StringIO()
        for line in self.seqLines():
            seq.write(line.rstrip("\n"))
        s = seq.getvalue()
        seq.close()
        if format == 'array':
            s = numpy.fromstring(s,dtype='S1')
        return s

    def seqLen(self):
        n = 0
        for line in self.seqLines():
            n += len(line) - 1
            if not line.endswith("\n"):
                n += 1
        return n

    def lineLen(self):
        return self.maxLineLen

    def close(self):
        if self.ownInfile:
            self.infile.close()
 

if __name__ == '__main__':

    if len(sys.argv) != 5:
        print "USAGE: python splitter.py fasta_file queryLen overlapBP outfile\n"
        print "       - fasta_file: input file\n"
        print "       - queryLen: length (bp) of splitted sequences\n"
        print "       - overlapBP: length (bp) of overlap\n"
        print "       - outfile: output file\n"
        sys.exit(0)

    seqFileName = sys.argv[1]       ## SEQ FILE NAME
    queryLen = int(sys.argv[2])    ## CHUNK SIZE
    overlapBP = int(sys.argv[3])   ## OVERLAP BPorigGI
    outFileName = sys.argv[4]
    
    
    seqFile = open(seqFileName, "r")
    outFile = open(outFileName, "w")

    listcount = 0
    cid = 0         ## CHUNK UNIQUE ID
    sid = 0         ## seq id

    for rec in FastaReader(open(seqFileName, 'r')).records():

        sid += 1  
        start = 0
        defLinePart = rec.header().strip().split()[0][1:] 
        currSeq = rec.sequence().strip()
        seqLen = len(currSeq)
            
        ##
        ## 
        ##
        while True:
            
            end = start + queryLen
            cid += 1
            #print "cid=%d, start=%d, end=%d, cnt=%d, seqLen=%d, sid=%d" % (cid, start, end, cnt, seqLen, sid)
            
            if start == 0 and end < seqLen:
                newHeader = ">"+defLinePart+"_"+str(cid)+"_"+"0"+"_"+str(start)+"_"+str(end)+"_"+str(start)+"_"+str(end)
                #print newHeader
                outFile.write(newHeader+"\n")
                outFile.write(currSeq[start:end].upper()+"\n")
                
            elif start == 0 and end >= seqLen:
                newHeader =  ">"+defLinePart+"_"+str(cid)+"_"+"1"+"_"+str(start)+"_"+str(seqLen)+"_"+str(start)+"_"+str(seqLen)
                #print newHeader
                outFile.write(newHeader+"\n")
                outFile.write(currSeq[start:seqLen].upper()+"\n")
                break

            elif start > 0 and end < seqLen:
                newHeader = ">"+defLinePart+"_"+str(cid)+"_"+"2"+"_"+str(start)+"_"+str(end)+"_"+str(start)+"_"+str(end)
                #print newHeader
                outFile.write(newHeader+"\n")
                upper = currSeq[start:end].upper()
                outFile.write(upper+"\n")

            elif start > 0 and end >= seqLen:
                newHeader =  ">"+defLinePart+"_"+str(cid)+"_"+"3"+"_"+str(start)+"_"+str(seqLen)+"_"+str(start)+"_"+str(seqLen)
                #print newHeader
                outFile.write(newHeader+"\n")
                if (seqLen - start > overlapBP):
                    upper = currSeq[start:seqLen].upper()
                    outFile.write(upper+"\n")
                else:
                    outFile.write(currSeq[start:seqLen]+"\n")
                break
                
            else:
                print "What?", seqLen, start, end, cur_record.id
                sys.exit(0)
            
            start = end - overlapBP ################ CAUSTION!
            
        
    outFile.close()
    seqFile.close()
    print "Total number of input sequences = ", sid
    print "Query length = ", queryLen, "bp"
    print "Overlap length = ", overlapBP, "bp"
    print "Total %d queries are generated and saved in %s." % (cid, outFileName)


## EOF
