#!/usr/bin/env python

### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##
#
# See COPYING file distributed along with the MGTAXA package for the
# copyright and license terms.
#
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##

##
## This tool is only for converting BLAST hits for classifier project.
## Each hit has two more fields, percentage identity and percentage coverage
## for query. Those values will be extracted and added to the final converted
## output file.
##

import sys
import os
import struct
import string
import optparse
from sqlite3 import *

if __name__ == '__main__':
   
    usage = "python load_sql.py -b binDir -o outSqlFile -d deflineOpt -i inDefFile"
    parser = optparse.OptionParser(usage=usage)
    parser.add_option("-b", "--bin-dir", dest="directory", 
                  action="store", type="string", help="path to *.bin files")
    parser.add_option("-o", "--output", dest="outSqlFileName",  
                  action="store", type="string", help="output sqlite database file")
    parser.add_option("-i", "--in", dest="inDefFile",  
                  action="store", type="string", help="input def file")
    parser.add_option("-d", "--with-defline", dest="deflineOpt",  default=0,
                  action="store", type="int", help="1: add original defline after qid, 0, if not")                  
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
        ##double      percIdent; 
        ##double      percCover;  
    ##} structBlResToSaveHits_t;
    
    print "Output Sqlite database file: ", outSqlFileName+".sqlite"
         
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
    ## Create database and table
    ##
    ##  qId         Unique query id (serial numner)
    ##  qIdDef      Each query's original defline from *.def file; set if bDefline = 1
    ##  sId         Subject id (either GI or part of defline)
    ##  dIdent      Percentage identity
    ##  alignLen    Alignment length
    ##  nMismatches Number of mismatches
    ##  nGaps       Number of gaps
    ##  qStart      Start of alignment in query
    ##  qEnd        End of alignment in query
    ##  sStart      Start of alignment in subject
    ##  sEnd        End of alignment in query
    ##  eValue      Evalue
    ##  bitScore    Bit score
    ##  percIdent   Percentage identity for query (identity count/query length * 100)
    ##  percCover   Percentage coverage for subject ((qend-qstart)/query length * 100)
    ##
    dbName = outSqlFileName + ".sqlite"
    conn = connect(dbName)
    curs = conn.cursor()
    curs.execute('''DROP TABLE if exists hits''')
    conn.commit()
    curs.execute('''CREATE TABLE IF NOT EXISTS hits
                ( 
                    qId         BIGINT,      
                    qIdDef      VARCHAR(40), 
                    sId         VARCHAR(40), 
                    dIdent      DOUBLE,      
                    alignLen    INT,         
                    nMismatches INT,         
                    nGaps       INT,         
                    qStart      INT,         
                    qEnd        INT,         
                    sStart      INT,         
                    sEnd        INT,         
                    eValue      DOUBLE,      
                    bitScore    DOUBLE,
                    percIdent   DOUBLE,
                    percCover   DOUBLE       
                )
                ''')
 
    ##    
    ## Read bin files and append to tables
    ##
    if bDefline:
        defFile = open(inDefFile, 'r')
        line = defFile.readline()
        
    totalHits = 0
    structSize = struct.calcsize('L40sdIIIIIIIdddd')
    for i in range(numHitFiles):
        subFileName = os.path.join(subDir,vecHitFileName[i])
        hitFile = open(subFileName, "rb")
        recordData = hitFile.read(structSize)
        numHits = 0
        while True:
            try:
                ## Load data from bin file
                s = struct.unpack('L40sdIIIIIIIdddd', recordData)
                totalHits += 1
                numHits += 1
                cmd = "insert into hits values (" \
                    + str(s[0]) + ",\"" \
                    
                ##
                ## Add the orig defline from .def file after the 'qid' field
                ##
                if bDefline:
                    while int(line.split()[0]) != int(s[0]):
                        line = defFile.readline()
                    cmd += line.split()[1][1:] + "\",\""
                else:
                    cmd += "\",\""
                    
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
                    + str(s[11]) + "," \
                    + str(s[12]) + "," \
                    + str(s[13]).strip() + ")" 
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
    if bDefline:
        defFile.close()
    
"""        
    ## test
    curs.execute("select count(*) from hits")
    print "num records = ", curs.fetchone()[0]
    curs.close()   
"""


    
## EOF
