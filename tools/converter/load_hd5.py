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
#print 'tables.__version__',t.__version__
import string
import optparse

if __name__ == '__main__':

    usage = "python load_hd5.py -b binDir -o outSqlFile"
    parser = optparse.OptionParser(usage=usage)
    parser.add_option("-b", "--bin-dir", dest="directory", 
                  action="store", type="string", help="path to *.bin files")
    parser.add_option("-o", "--output", dest="outHd5FileName",  
                  action="store", type="string", help="output HD5 database file")
    (options, args) = parser.parse_args()
    
    if options.directory and options.outHd5FileName:
        topDir = options.directory
        outHd5FileName = options.outHd5FileName
    else:
        parser.error("Please set the path to *.bin files and output file name.")
                    
    ###
    ### Define a user record to characterize some kind of particles
    ###
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
        qId         = t.UInt64Col()     # L
        sId         = t.StringCol(40)   # 40s 
        dIdent      = t.FloatCol()      # d: double
        alignLen    = t.UInt32Col()     # I
        nMismatches = t.UInt32Col()     # I
        nGaps       = t.UInt32Col()     # I
        qStart      = t.UInt32Col()     # I
        qEnd        = t.UInt32Col()     # I
        sStart      = t.UInt32Col()     # I
        sEnd        = t.UInt32Col()     # I
        eValue      = t.FloatCol()      # d: double
        bitScore    = t.FloatCol()      # d: double
            
    ###
    ### Load recordData and insert into db table
    ###
    
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
        
    ###
    ### load hit file names from hitfilelist
    ###
    vecHitFileName = []
    subDir = topDir  
    for f in os.listdir(subDir):
        if f.find(".bin") > -1: 
            vecHitFileName.append(f)
    numHitFiles = len(vecHitFileName)
    
    ###
    ### Read bin files and append to tables
    ###
    totalHits = 0
    structSize = struct.calcsize('L40sdIIIIIIIdd')
    for i in range(numHitFiles):
        subFileName = os.path.join(subDir,vecHitFileName[i])
        hitFile = open(subFileName, "rb")
        recordData = hitFile.read(structSize)
        numHits = 0
        while True:
            try:
                s = struct.unpack('L40sdIIIIIIIdd', recordData)
                totalHits += 1
                numHits += 1                   
                BlHits['qId']         = s[0]
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
        
    ###
    ### flush recordData
    ###
    table.flush()   # flush recordData in the table
    h5file.flush()  # flush all pending recordData
    print "Total number of hits = ",totalHits
    
    ## Create index
    #table.cols.sId.createIndex()
    #table.cols.qId.createIndex()
 
    h5file.close()
    
"""
    ## HDF5 file info
    print h5file
    
    ## Access columns
    table = h5file.root.blhits.blhitstab
    gi = [ x['gi'] for x in table.iterrows() ]
    sid = [ x['sId'] for x in table.iterrows() ]
    print gi[0], sid[0]
    
    ## Reading rows
    for r in table.iterrows():
        print r['gi'], r['eValue'], r['dIdent'], r['dCover'] 
        #print [r[i] for i in range(0,15)]
        
    h5file.close()
"""


### EOF
