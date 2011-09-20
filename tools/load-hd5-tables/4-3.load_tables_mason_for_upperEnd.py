#!/usr/bin/env python
#from sqlite3 import *
import sys
import os
import struct
from numpy import *
import tables as t
print 'tables.__version__',t.__version__
from sqlite3 import *

##
## Read all bin files from all sub directories and save the hits in a HD5 file
## While saving, the cut
##
if __name__ == '__main__':

    if len(sys.argv) != 4:
        print "python 4-3.load_tables_mason_for_upperEnd.py topDir out_file_name 0/1_for_saving_csv"
        print "Example:"
        print "python ./4-3.load_tables_mason_for_upperEnd.py . output.hd5 1"
        sys.exit(1)    
    
    topDir = sys.argv[1]    
    filename = sys.argv[2]
    bMakeCSV = int(sys.argv[3])

    ##
    ## Get the dir list under topDir
    ##
    dirName = []
    dirList = os.listdir(topDir)
    for na in dirList:
		if os.path.isdir(na):
			dirName.append(na)
    print dirName, len(dirName)
    numDir = len(dirName)

    ##
    ## Define a user record to characterize some kind of particles
    ##
    
    ## Orig struct in mrblast.cpp 
    """
    structBlResToSaveHitsMason {
        uint32_t    gi;
        uint32_t    readId;
        char        readStrand;     /// 'f' or 'r'
        char        readPairId;     /// '0' or '1'
        uint32_t    subjectId;
        uint32_t    qStart;
        uint32_t    qEnd;
        uint32_t    sStart;
        uint32_t    sEnd;
        double      eValue;
        float       bitScore;
        uint32_t    upperStart;
        uint32_t    upperEnd;
        float       identity;
        float       coverage;    
    }
    """
    
    class BlHits(t.IsDescription):  
        gi         = t.UInt32Col()
        readId     = t.UInt32Col()
        readStrand = t.StringCol(1)
        readPairId = t.StringCol(1) 
        sId        = t.UInt32Col()
        qStart     = t.UInt32Col()
        qEnd       = t.UInt32Col()
        sStart     = t.UInt32Col()
        sEnd       = t.UInt32Col()        
        eValue     = t.FloatCol()       ## double-precision ==> 'd'
        bitScore   = t.UInt32Col()      ## single-precision ==> 'f'
        upperStart = t.UInt32Col()
        upperEnd   = t.UInt32Col()
        dIdent     = t.Float32Col()     ## single-precision ==> 'f'
        dCover     = t.Float32Col()     ## single-precision ==> 'f'
            
    ##
    ## Load recordData and insert into db table
    ##
    
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
        ##
        ## load hit file names from hitfilelist
        ##
        vecHitFileName = []
        fastaFileName = ""
        subDir = topDir+"/"+dirName[j]
        for f in os.listdir(subDir):
            if f.find(".bin") > -1: 
                vecHitFileName.append(f)
        numHitFiles = len(vecHitFileName)
        
        ##
        ## Connect deflines.hd5 to update cutEnd
        ##
        dbName = topDir+"/"+dirName[j]+"/deflines.hd5"
        deflines_h5file = t.openFile(dbName)
        ## Create a new group under "/" (root)
        deflines_table = deflines_h5file.root.deflines.deflinetab
        
        ##
        ## Read bin files and append to tables
        ##
        structSize = struct.calcsize('IIccIIIIIdfIIffI')
        print "dir name = ", dirName[j]
        for i in range(numHitFiles):
            subFileName = subDir+"/"+vecHitFileName[i]
            hitFile = open(subFileName, "rb")
            recordData = hitFile.read(structSize)
            numHits = 0
            while True:
                try:
                    s = struct.unpack('IIccIIIIIdfIIffI', recordData)
                    
                    ## SQL version
                    #curs.execute('''CREATE TABLE IF NOT EXISTS item
                                  #( gi          integer,
                                    #readId      integer, 
                                    #pairId      integer, 
                                    #strandId    integer,
                                    #readLen     integer,
                                    #cutStart     integer, 
                                    #cutEnd      integer,
                                    #readStart   integer,
                                    #readEnd     integer
                                   #)''')
                    #if int(s[0]) != temp_gi and int(s[1]) != temp_readId and int(s[3]) != temp_pairId:
                        #sqlCmd = "select cutStart, cutEnd from item where gi = " + \
                                #str(s[0]) + " and readId = " +  \
                                #str(s[1]) + " and pairId = " + str(s[3])
                        #curs.execute(sqlCmd)
                        #print curs.fetchone()[0], curs.fetchone()[1]
                        #temp_cutStart = curs.fetchone()[0]
                        #temp_cutEnd   = curs.fetchone()[1] 
                    
                    ## HD5 version
                    mycond = "(gi == " + str(s[0]) + ") & (readId == " + str(s[1]) + ") & (pairId == " + str(s[3]) + ")"
                    new_cutStart = [ x['cutStart'] for x in deflines_table.where(mycond)]
                    new_cutEnd = [ x['cutEnd'] for x in deflines_table.where(mycond)]
                    assert len(new_cutStart) == 1
                    assert len(new_cutEnd) == 1
                    
                    BlHits['gi']         = s[0]
                    BlHits['readId']     = s[1]
                    BlHits['readStrand'] = s[2]
                    BlHits['readPairId'] = s[3]                    
                    BlHits['sId']        = s[4]
                    BlHits['qStart']     = s[5]
                    BlHits['qEnd']       = s[6]
                    BlHits['sStart']     = s[7]
                    BlHits['sEnd']       = s[8]
                    BlHits['eValue']     = s[9]
                    BlHits['bitScore']   = s[10]
                    BlHits['upperStart'] = int(new_cutStart[0])
                    BlHits['upperEnd']   = int(new_cutEnd[0])
                    BlHits['dIdent']     = s[13]
                    BlHits['dCover']     = s[14]
                    BlHits.append()
                    numHits += 1
                    
                    if bMakeCSV:
                        csvString = str(s[0]) + "," \
                            + str(s[1]) + "," \
                            + str(s[2]) + "," \
                            + str(s[3]) + "," \
                            + str(s[4]) + "," \
                            + str(s[5]) + "," \
                            + str(s[6]) + "," \
                            + str(s[7]) + "," \
                            + str(s[8]) + "," \
                            + str(s[9]) + "," \
                            + str(s[10]) + "," \
                            + str(s[11]) + "," \
                            + str(s[12]) + "," \
                            + str(s[13]) + "," \
                            + str(s[14]) + "\n"
                        csvFile.write(csvString)
                    recordData = hitFile.read(structSize)
                except:
                    break
            
            hitFile.close()
            print "num hits %d in %s" % (numHits, vecHitFileName[i])
		
        #curs.close()		
        deflines_h5file.close()
        totalHits += numHits
        
        ##
        ## flush recordData
        ##
        table.flush()   # flush recordData in the table
        h5file.flush()  # flush all pending recordData
        print "%d hits in %s" %(totalHits, dirName[j])
        
    print "total num hits = ",totalHits
    if bMakeCSV:
        csvFile.close() 
    
    ## Create index
    table.cols.gi.createIndex()
    table.cols.sId.createIndex()    
    table.cols.readId.createIndex()
    table.cols.readPairId.createIndex()
    
    h5file.close()
    
""" Read Test

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
    
    h5file.close()
"""
    
## EOF
