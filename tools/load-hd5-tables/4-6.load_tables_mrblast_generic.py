#!/usr/bin/env python
#from sqlite3 import *
import sys
import os
import struct
from numpy import *
import tables as t
print 'tables.__version__',t.__version__
#from numexpr import *
import string

def removeNonAscii(s): return "".join(i for i in s if ord(i)<128)


if __name__ == '__main__':

    if len(sys.argv) != 5:
        print "python 4-2.load_tables_mrblast_generic.py topDir out_file_name 0/1_for_saving_csv num_div"
        sys.exit(1)    
    
    topDir   = sys.argv[1]    
    filename = sys.argv[2]
    bMakeCSV = int(sys.argv[3])
    numDir   = int(sys.argv[4])
    
    dirName = [ "div"+str(x) for x in range(numDir) ]

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
    h5file = t.openFile(filename, mode = "w", title = "BLAST hits")
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
    for j in range(numDir):
        ###
        ### load hit file names from hitfilelist
        ###
        vecHitFileName = []
        subDir = topDir + "/" + dirName[j]
        for f in os.listdir(subDir):
            if f.find(".bin") > -1: 
                vecHitFileName.append(f)
        numHitFiles = len(vecHitFileName)
        
        ###
        ### Read bin files and append to tables
        ###
        #structSize = struct.calcsize('IIccIIIIIdfIIffI')
        #structSize = struct.calcsize('QQIIIIdfff')
        structSize = struct.calcsize('L80sdIIIIIIIdd')
        print "struct size = ", structSize
        print "dir name = ", dirName[j]
        
        #numHits = 0
        for i in range(numHitFiles):
            subFileName = subDir+"/"+vecHitFileName[i]
            hitFile = open(subFileName, "rb")
            recordData = hitFile.read(structSize)
            numHits = 0
            while True:
                try:
                    #s = struct.unpack('IIccIIIIIdfIIffI', recordData)
                    #s = struct.unpack('QQIIIIdfff', recordData)
                    s = struct.unpack('L80sdIIIIIIIdd', recordData)
                    #print s[0], s[2], s[3]
                    #BlHits['qId']         = s[0]
                    ##BlHits['readId']     = s[1]
                    ##BlHits['readStrand'] = s[2]
                    ##BlHits['readPairId'] = s[3]                    
                    #BlHits['sId']        = s[1]
                    #BlHits['qStart']     = s[2]
                    #BlHits['qEnd']       = s[3]
                    #BlHits['sStart']     = s[4]
                    #BlHits['sEnd']       = s[5]
                    #BlHits['eValue']     = s[6]
                    #BlHits['bitScore']   = s[7]
                    ##BlHits['upperStart'] = s[11] 
                    ##BlHits['upperEnd']   = s[12]
                    #BlHits['dIdent']     = s[8]
                    #BlHits['dCover']     = s[9]
                    #BlHits.append()
                    totalHits += 1
                    numHits += 1
                    
                    #print len(str(s[1])), removeNonAscii(str(s[1])), len(removeNonAscii(str(s[1])))
                    #print str(s[1]).encode("ascii", "ignore"), len(str(s[1]).encode("ascii", "ignore"))
                    print filter(lambda x: x in string.printable, str(s[1])), len(filter(lambda x: x in string.printable, str(s[1])))

                    
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
                        print csvString
                    recordData = hitFile.read(structSize)
                except:
                    break
            
            hitFile.close()
            print "num hits %d in %s" % (numHits, vecHitFileName[i])
        
        #totalHits += numHits
        
        ###
        ### flush recordData
        ###
        table.flush()   # flush recordData in the table
        h5file.flush()  # flush all pending recordData
        print "%d hits in %s" %(totalHits, dirName[j])
        
    print "total num hits = ",totalHits
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
    table = h5file.root.blhits.blhits
    gi = [ x['gi'] for x in table.iterrows() ]
    sid = [ x['sId'] for x in table.iterrows() ]
    #names = [ x['name'] for x in table if x['TDCcount'] > 3 and 20 <= x['pressure'] < 50 ]
    #names = [ x['name'] for x in table.where('(TDCcount > 3) & (20 <= pressure) & (pressure < 50)') ]

    #print gi, sid
    
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
