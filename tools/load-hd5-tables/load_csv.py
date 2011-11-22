#!/usr/bin/env python

import sys
import os
import struct
from numpy import *
import string

if __name__ == '__main__':

    if len(sys.argv) != 3:
        print "python load_hd5.py topDir out_file_name"
        sys.exit(1)    
    
    topDir   = sys.argv[1]    
    filename = sys.argv[2]

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
        
    print "Total number of hits = ",totalHits
    csvFile.close() 
     
### EOF
