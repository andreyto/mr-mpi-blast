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
import string
from sqlite3 import *
import optparse
import linecache

if __name__ == '__main__':
   
    usage = "python load_sql.py -b binDir -o outSqlFile -d deflineOpt -i inDefFile"
    parser = optparse.OptionParser(usage=usage)
    parser.add_option("-b", "--bin-dir", dest="directory", 
                  action="store", type="string", help="path to *.bin files")
    parser.add_option("-o", "--output", dest="outSqlFileName",  
                  action="store", type="string", help="output sqlite database file")
    parser.add_option("-d", "--with-defline", dest="deflineOpt",  default=0,
                  action="store", type="int", help="1: add original defline after qid, 0, if not")                  
    parser.add_option("-i", "--in", dest="inDefFile",  
                  action="store", type="string", help="input def file")
    (options, args) = parser.parse_args()
    
    if options.directory and options.outSqlFileName:
        topDir = options.directory
        outSqlFileName = options.outSqlFileName
    else:
        parser.error("Please set the path to *.bin files and output file name.")
    
    bDefline = 0
    if options.deflineOpt and options.inDefFile:
        bDefline = options.deflineOpt 
        inDefFile = options.inDefFile
    elif options.deflineOpt and not options.inDefFile:
        parser.error("Please set the input defline file")

    ###
    ### If deflineOpt is set, construct a dict for <qid, lineNo> for retirieving
    ### defline when making final tabular output format.
    ###
    qidDict = {}
    lineNum = 1    
    if bDefline:
        with open(inDefFile) as infile:
            for line in infile:
                qidDict[int(line.split()[0])] = lineNum
                lineNum += 1
                        
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
    
    print "Output Sqlite database file: ", outSqlFileName+".sqlite"
         
    ###
    ### load hit file names from hitfilelist
    ###
    vecHitFileName = []
    subDir = topDir  
    for f in os.listdir(subDir):
        if f.find(".bin") > -1: 
            vecHitFileName.append(f)
    numHitFiles = len(vecHitFileName)
    
    ## Create database and table
    dbName = outSqlFileName + ".sqlite"
    conn = connect(dbName)
    curs = conn.cursor()
    curs.execute('''DROP TABLE if exists hits''')
    conn.commit()
    if bDefline:
        curs.execute('''CREATE TABLE IF NOT EXISTS hits
        ( 
            qId         BIGINT,      
            qIdDef      VARCHAR(80), 
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
    else:
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
       
    ###    
    ### Read bin files and append to tables
    ###
    totalHits = 0
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
                    + str(s[0]) + ",\"" \
                    
                ### Add the orig defline from .def file after the 'qid' field
                if bDefline:
                    cmd += linecache.getline(inDefFile, qidDict[int(s[0])]).split()[1][1:] + "\",\""                         
                    
                cmd += filter(lambda x: x in string.printable, str(s[1])) + "\"," \
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
                recordData = hitFile.read(structSize)
            except struct.error:
                break
            except:
                print "Unexpected error:", sys.exc_info()[0]
                raise
        conn.commit()  
        hitFile.close()
        print "Number of hits = %d in %s" % (numHits, vecHitFileName[i])
    print "Total number of hits = ",totalHits
 
    curs.close()  
    conn.close()
    
"""        
    ### test
    curs.execute("select count(*) from hits")
    print "num records = ", curs.fetchone()[0]
    curs.close()   
"""


    
### EOF
