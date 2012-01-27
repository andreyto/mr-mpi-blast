#!/usr/bin/env python

### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##
#
# See COPYING file distributed along with the MGTAXA package for the
# copyright and license terms.
#
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##

import sys
import os
import struct
import tables as t
import string
import optparse

if __name__ == '__main__':

    usage = "python load_hd5.py -b binDir -o outSqlFile"
    parser = optparse.OptionParser(usage=usage)
    parser.add_option("-b", "--bin-dir", dest="directory", 
                  action="store", type="string", help="path to *.bin files")
    parser.add_option("-o", "--output", dest="outHd5FileName",  
                  action="store", type="string", help="output HD5 database file")
    parser.add_option("-i", "--in", dest="inDefFile",  
                  action="store", type="string", help="input *.def file name")
    parser.add_option("-d", "--with-defline", dest="deflineOpt", default=0,
                  action="store", type="int", help="1: add original defline after qid; 0, if not")                                    
    (options, args) = parser.parse_args()
    
    if options.directory and options.outHd5FileName:
        topDir = options.directory
        outHd5FileName = options.outHd5FileName
    else:
        parser.error("Please set the path to *.bin files and output file name.")

    bDefline = 0
    if options.deflineOpt and options.inDefFile:
        bDefline = options.deflineOpt 
        inDefFile = options.inDefFile
    elif options.deflineOpt and not options.inDefFile:
        parser.error("Please set the input defline file.")
 
    ##
    ## Define a user record to characterize some kind of particles
    ##
    ##typedef struct structBlResToSaveHits {
        ##uint64_t    queryId;
        ##char        subjectId[40];
        ##double      identity;
        ##uint32_t    alignLen;
        ##uint32_t    nMismatches;
        ##uint32_t    nGaps;
        ##uint32_t    qStart;
        ##uint32_t    qEnd;
        ##uint32_t    sStart;
        ##uint32_t    sEnd;
        ##double      eValue;
        ##double      bitScore;
    ##} structBlResToSaveHits_t;

    class BlHits(t.IsDescription):
        qId         = t.UInt64Col()     # L          # Unique query id (serial numner)
        qIdDef      = t.StringCol(40)   # 40s        # Each query's original defline from *.def file; set if bDefline = 1
        sId         = t.StringCol(40)   # 40s        # Subject id (either GI or part of defline)
        dIdent      = t.FloatCol()      # d: double  # Percentage identity
        alignLen    = t.UInt32Col()     # I          # Alignment length
        nMismatches = t.UInt32Col()     # I          # Number of mismatches
        nGaps       = t.UInt32Col()     # I          # Number of gaps
        qStart      = t.UInt32Col()     # I          # Start of alignment in query
        qEnd        = t.UInt32Col()     # I          # End of alignment in query
        sStart      = t.UInt32Col()     # I          # Start of alignment in subject
        sEnd        = t.UInt32Col()     # I          # End of alignment in query
        eValue      = t.FloatCol()      # d: double  # Evalue
        bitScore    = t.FloatCol()      # d: double  # Bit score
            
    ##
    ## Load recordData and insert into db table
    ##
    
    ## Open a file in "w"rite mode
    h5file = t.openFile(outHd5FileName+".hd5", mode = "w", title = "BLAST hits")
    ## Create a new group under "/" (root)
    root = h5file.root
    group = h5file.createGroup(root, "blhits", "blhits")
    ## Create one table on it
    table = h5file.createTable(group, "blhitstab", BlHits, "blhitstab")
    ## Fill the table with 10 particles
    BlHits = table.row

    print "Output HD5 database file: ", outHd5FileName+".hd5"
        
    ##
    ## load hit file names from hitfilelist
    ##
    vecHitFileName = []
    subDir = topDir  
    for f in os.listdir(subDir):
        if f.find(".bin") > -1: 
            vecHitFileName.append(f)
    vecHitFileName.sort()
    numHitFiles = len(vecHitFileName)
    
    ##
    ## Read bin files and append to tables
    ##
    if bDefline:
        defFile = open(inDefFile, 'r')
        line = defFile.readline()
        
    totalHits = 0
    structSize = struct.calcsize('L40sdIIIIIIIdd')
    for i in range(numHitFiles):
        subFileName = os.path.join(subDir,vecHitFileName[i])
        hitFile = open(subFileName, "rb")
        recordData = hitFile.read(structSize)
        numHits = 0
        while True:
            try:
                ## Load data from bin file
                s = struct.unpack('L40sdIIIIIIIdd', recordData)
                totalHits += 1
                numHits += 1                   
                BlHits['qId']         = s[0]

                ##
                ## Add the orig defline from .def file after the 'qid' field
                ##
                if bDefline:
                    while int(line.split()[0]) != int(s[0]):
                        line = defFile.readline()
                    BlHits['qIdDef']  = line.split()[1][1:]
                else:
                    BlHits['qIdDef']  = ""
                
                BlHits['sId']         = filter(lambda x: x in string.printable, str(s[1]))
                BlHits['dIdent']      = s[2]
                BlHits['alignLen']    = s[3]
                BlHits['nMismatches'] = s[4]
                BlHits['nGaps']       = s[5]              
                BlHits['qStart']      = s[6]
                BlHits['qEnd']        = s[7]
                BlHits['sStart']      = s[8]
                BlHits['sEnd']        = s[9]
                BlHits['eValue']      = s[10]
                BlHits['bitScore']    = s[11]
                BlHits.append()
                recordData = hitFile.read(structSize)
            except struct.error:
                break
            except:
                print "Unexpected error:", sys.exc_info()[0]
                raise
        hitFile.close()
        print "Number of hits = %d in %s" % (numHits, vecHitFileName[i])
        
    ##
    ## flush recordData
    ##
    table.flush()   # flush recordData in the table
    h5file.flush()  # flush all pending recordData
    print "Total number of hits = ",totalHits
    
    ## Create index
    #table.cols.sId.createIndex()
    #table.cols.qId.createIndex()

  
 
    h5file.close()
    if bDefline:
        defFile.close()
    
"""
    ## HDF5 file info
    print h5file
    
    ## Access columns
    table = h5file.root.blhits.blhitstab
    gi = [ x['qId'] for x in table.iterrows() ]
    sid = [ x['sId'] for x in table.iterrows() ]
    print gi[0], sid[0]
    
    ## Reading rows
    for r in table.iterrows():
        print r['qId'], r['eValue'], r['eValue'], r['bitScore'] 
        #print [r[i] for i in range(0,15)]
        
    h5file.close()
"""


## EOF
