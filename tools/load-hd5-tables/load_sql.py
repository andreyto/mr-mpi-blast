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
import string
from sqlite3 import *

if __name__ == '__main__':

    if len(sys.argv) != 4:
        print "python load_hd5.py topDir out_file_name 0/1_for_saving_csv" 
        sys.exit(1)    
        
    topDir   = sys.argv[1]    
    outFileName = sys.argv[2]
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
    
    if bMakeCSV:
        csvFileName = outFileName + ".csv"
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
    print vecHitFileName
    
    ## Create database and table
    dbName = outFileName + ".sqlite"
    print "DB name = ", dbName
    conn = connect(dbName)
    curs = conn.cursor()
    curs.execute('''DROP TABLE if exists hits''')
    conn.commit()
    curs.execute('''CREATE TABLE IF NOT EXISTS hits
        ( 
            qId         BIGINT,      
            sId         VARCHAR(80), 
            dIdent      DOUBLE,      
            alignLen    INT,         
            nMismatches INT,         
            nGaps       INT,         
            qStart      INT,         
            qEnd        INT,         
            sStart      INT,         
            sEnd        INT,         
            eValue      DOUBLE,      
            bitScore    DOUBLE       
        )
        ''')
       
        
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
                cmd = "insert into hits values (" \
                    + str(s[0]) + "," \
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
                    + str(s[11]).strip() + ")" 
                curs.execute(cmd)
                    
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
        conn.commit()  
        hitFile.close()
        print "Number of hits = %d in %s" % (numHits, vecHitFileName[i])
    print "Total number of hits = ",totalHits
    
    if bMakeCSV:
        csvFile.close() 
    curs.close()  
    conn.close()
    
"""        
    ### test
    curs.execute("select count(*) from hits")
    print "num records = ", curs.fetchone()[0]
    curs.close()   
"""


    
### EOF
