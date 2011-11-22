### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##
#
# See COPYING file distributed along with the MGTAXA package for the
# copyright and license terms.
#
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##

#!/usr/bin/env python

import sys
import os
import struct
from numpy import *
import tables as t
print 'tables.__version__',t.__version__
import string

if __name__ == '__main__':

    if len(sys.argv) != 4:
        print "python load_hd5.py topDir out_file_name 0/1_for_saving_csv"
        sys.exit(1)    
    topDir   = sys.argv[1]    
    filename = sys.argv[2]
    bMakeCSV = int(sys.argv[3])

    ###
    ### Define a user record to characterize some kind of particles
    ###
    ##typedef struct structBlResToSaveHits {
        ##uint64_t    queryId;
        ##char        subjectId[80];
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
        sId         = t.StringCol(80)   # 80s 
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
    h5file = t.openFile(filename+".hd5", mode = "w", title = "BLAST hits")
    ## Create a new group under "/" (root)
    root = h5file.root
    group = h5file.createGroup(root, "blhits", "blhits")
    ## Create one table on it
    table = h5file.createTable(group, "blhitstab", BlHits, "blhitstab")
    ## Fill the table with 10 particles
    BlHits = table.row

    if bMakeCSV:
        csvFileName = filename + ".csv"
        csvFile = open(csvFileName, 'w')
    
    totalHits = 0
    
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
    structSize = struct.calcsize('L80sdIIIIIIIdd')
    for i in range(numHitFiles):
        subFileName = os.path.join(subDir,vecHitFileName[i])
        hitFile = open(subFileName, "rb")
        recordData = hitFile.read(structSize)
        numHits = 0
        while True:
            try:
                s = struct.unpack('L80sdIIIIIIIdd', recordData)
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
                
                if bMakeCSV:
                    csvString = str(s[0]) + "," \
                                + filter(lambda x: x in string.printable, str(s[1])) + "," \
                                + str(s[2]) + "," \
                                + str(s[3]) + "," \
                                + str(s[4]) + "," \
                                + str(s[5]) + "," \
                                + str(s[6]) + "," \
                                + str(s[7]) + "," \
                                + str(s[8]) + "," \
                                + str(s[9]) + "," \
                                + str(s[10]) + "," \
                                + str(s[11]).strip() + "\n"
                    csvFile.write(csvString)
                recordData = hitFile.read(structSize)
            except:
                break
        hitFile.close()
        print "Number of hits = %d in %s" % (numHits, vecHitFileName[i])
        
    ###
    ### flush recordData
    ###
    table.flush()   # flush recordData in the table
    h5file.flush()  # flush all pending recordData
    print "Total number of hits = ",totalHits
    if bMakeCSV:
        csvFile.close() 
    
    ## Create index
    #table.cols.sId.createIndex()
    #table.cols.qId.createIndex()
 
    h5file.close()
    
"""
    ## HDF5 file info
    print h5file
    
    ## Access columns
    table = h5file.root.blhits.blhitstab
    gi = [ x['qId'] for x in table.iterrows() ]
    sid = [ x['sId'] for x in table.iterrows() ]
    #names = [ x['name'] for x in table if x['TDCcount'] > 3 and 20 <= x['pressure'] < 50 ]
    #names = [ x['name'] for x in table.where('(TDCcount > 3) & (20 <= pressure) & (pressure < 50)') ]

    print gi[0], sid[0]
    # Reading rows
    #for r in table.iterrows():
        #print r['qId']
        
    h5file.close()
"""

"""
    ## Create column object and read
    gcolumns = h5file.createGroup(root, "columns", "")   
    h5file.createArray(gcolumns, 'gi', array(gi), "gi array")
    h5file.createArray(gcolumns, 'sid', array(sid), "sid array")
    #print h5file    
    gi2 = h5file.root.columns.gi.read()
    #print gi2
    
    ## Reading rows
    #for r in table.iterrows():
        #print r['gi']
    
    ## Access recordData
    print table.cols.gi[0]
    
    ## Create index
    table.cols.gi.createIndex()
    table.cols.sId.createIndex()
    print h5file 
    
"""
    
### EOF
