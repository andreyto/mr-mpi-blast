#!/usr/bin/env python

### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##
#
# See COPYING file distributed along with the MGTAXA package for the
# copyright and license terms.
#
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##
 
import optparse

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

    def seqChunks(self,chunkSize):
        seq = StringIO()
        for line in self.seqLines():
            seq.write(line.rstrip("\n"))
            if seq.tell() >= chunkSize:
                yield seq.getvalue()
                seq.close()
                seq = StringIO()
        if seq.tell() > 0:
            yield seq.getvalue()
        seq.close()

    #def seqArrays(self,chunkSize):
        #for s in self.seqChunks(chunkSize):
            #yield numpy.fromstring(s,dtype='S1')

    def sequence(self,format='str'):
        seq = StringIO()
        for line in self.seqLines():
            seq.write(line.rstrip("\n"))
        s = seq.getvalue()
        seq.close()
        #if format == 'array':
            #s = numpy.fromstring(s,dtype='S1')
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

version = '2.0'
verbose = False
inFileName = ''
outFileName = ''
defFileName = ''
deflineOption = 0
uidOption = 0
seqUID = 0

if __name__ == '__main__':
     
    usage = "python seqindexer.py -i inFile -o outIndexFile -d outDeflineFile -u uidopt -s startNum -b full_or_part_defline"
    parser = optparse.OptionParser(usage=usage)
    parser.add_option("-i", "--input", dest="infilename", 
                  action="store", type="string", help="input fasta file", default=False)
    parser.add_option("-o", "--output", dest="outfilename",  
                  action="store", type="string", help="output index file")
    parser.add_option("-d", "--deflinefilename", dest="deflinefilename",  
                  action="store", type="string", help="defline file name")                  
    parser.add_option("-u", "--uidopt", dest="uidopt",  
                  action="store", type="int", help="uid choice: 0=serial number, 1=gi")
    parser.add_option("-s", "--startno", dest="startno", 
                action="store", type="int", help="uid start number when -u 0")    
    parser.add_option("-b", "--deflineopt", dest="deflineopt",  
                  action="store", type="int", help="defline saving option: 0=part, 1=full")
    (options, args) = parser.parse_args()
 
    if options.infilename:
        inFileName = options.infilename
    else:
        parser.error("Please set the input file name.")
        
    if options.outfilename:
        outFileName = options.outfilename
    else:
        parser.error("Please set the output index file name.")
        
    if options.uidopt is not None:
        uidOption = options.uidopt        
        if uidOption == 0:
            if options.startno != 0:          
                seqUID = options.startno
            print "    Serial number starting from ", options.startno, "is used for unique query ID."
        else:
            print "    GI number is used for unique query ID."

    if options.deflinefilename and options.deflineopt is not None:
        defFileName = options.deflinefilename
        deflineOption = options.deflineopt
    else:
        parser.error("Please set the output def file name.")
        
    numSeq = 0
    currLoc = 0    
    outFile = open(outFileName, "w")
    defFile = open(defFileName, "w")

    for rec in FastaReader(open(inFileName, 'r')).records():
        defline = rec.header()
        
        ### 
        ### Save start location of each query def line and the length
        ###
        loc = currLoc
        currLoc += len(defline)
        seqLen = 0
        for line in rec.seqLines():
            seqLen += len(line.rstrip("\n"))
            currLoc += len(line)
        numSeq += 1
        seqUID += 1
        uid = 0
        
        if uidOption == 1:
            assert defline.rstrip().split("|")[0].find('gi') != -1
            defline.rstrip().split("|")[0]
            uid = defline.rstrip().split("|")[1]
        else:
            uid = seqUID
        outFile.write(str(loc)+"\t"+str(seqLen)+"\t"+str(uid)+"\n")
        deflines = ''
        
        if deflineOption == 0:
            defline2 = defline.rstrip().split(" ")[0]
        else:
            defline2 = defline.rstrip()        
        defFile.write(str(uid)+"\t"+defline2+"\n")
        
    outFile.close()
    defFile.close()
    print "    Total number of sequences = ", numSeq
    
### EOF
