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
import optparse
import pdb

def main():

    usage = "Run 'python load_csv.py --help' for options."
    parser = optparse.OptionParser(usage=usage)
    parser.add_option("-b", "--bin-dir", dest="directory", 
                  action="store", type="string", help="set path to *.bin files")
    parser.add_option("-o", "--output", dest="outfilename",  
                  action="store", type="string", help="set output CSV file name")
    parser.add_option("-i", "--in", dest="inDefFile",  
                  action="store", type="string", help="input *.def file name")
    parser.add_option("-d", "--with-defline", dest="deflineOpt", default=0,
                  action="store", type="int", help="1: add original defline after qid; 0, if not")                                    
    parser.add_option("-c", "--classifier-mode", dest="classifierMode", default=0,
            action="store", type="int", help="1: process output from BLAST tool executed in classifier mode; 0: otherwise")                                    
    (options, args) = parser.parse_args()
    
    if options.directory and options.outfilename:
        topDir = options.directory
        outCsvFileName = options.outfilename
    else:
        parser.error("Please set the path to *.bin files and output file name.")
    
    bDefline = 0
    if options.deflineOpt and options.inDefFile:
        bDefline = options.deflineOpt 
        inDefFile = options.inDefFile
    elif options.deflineOpt and not options.inDefFile:
        parser.error("Please set the input defline file.")
  
    ##
    ## Expected record format
    ##
    ## Unique query id (serial numner)
    ## Each query's original defline from *.def file; set if bDefline = 1
    ## Subject id (either GI or part of defline)
    ## Percentage identity
    ## Alignment length
    ## Number of mismatches
    ## Number of gaps
    ## Start of alignment in query
    ## End of alignment in query
    ## Start of alignment in subject
    ## End of alignment in query
    ## Evalue
    ## Bit score
    ## --- In classifier mode, these fields are also added: ---
    ## Percent identity over entire query 
    ## Percent coverage over entire query
    
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
        ## --- In classifier mode, these fields are also added: ---
        ##double      percIdent;
        ##double      percCover;
    ##} structBlResToSaveHits_t;
 
    csvFileName = outCsvFileName + ".csv"
    csvFile = open(csvFileName, 'w')
    print "Output CSV file: ", csvFileName
    
    ##
    ## load hit file names from hitfilelist
    ##
    vecHitFileName = []
    subDir = topDir  
    for f in os.listdir(subDir):
        if f.endswith(".bin"): 
            vecHitFileName.append(f)
    vecHitFileName.sort()
    numHitFiles = len(vecHitFileName)
    
    ##
    ## Read bin files and append to tables
    ##
    if bDefline:
        defFile = open(inDefFile, 'r')
        line = defFile.readline()
    structDef = 'L40sdIIIIIIIdd'
    if options.classifierMode:
        structDef += 'dd'
    totalHits = 0
    structSize = struct.calcsize(structDef)
    for i in range(numHitFiles):
        subFileName = os.path.join(subDir,vecHitFileName[i])
        hitFile = open(subFileName, "rb")
        numHits = 0
        while True:
            recordData = hitFile.read(structSize)
            lenRead = len(recordData)
            if lenRead == 0:
                break
            elif lenRead != structSize:
                raise ValueError("Partial record read from "+\
                    "binary output file which must be corrupt.")
            ## Load data from bin file
            s = struct.unpack(structDef, recordData)
            totalHits += 1
            numHits += 1
            qIid = int(s[0])
            csvString = str(qIid) + ","

            ##
            ## Add the orig defline from .def file after the 'qid' field
            ##
            if bDefline:
                while int(line.split()[0]) != qIid:
                    line = defFile.readline()
                csvString += line.split()[1][1:] + ","
                
            csvString += s[1].partition(b'\0')[0] + "," \
                    + ','.join((("%s" % (x,)) for x in s[2:])) + "\n"
            csvFile.write(csvString)
        hitFile.close()
        print "Number of hits = %d in %s" % (numHits, vecHitFileName[i])
            
    print "Total number of hits = ",totalHits
    csvFile.close() 
    if bDefline:
        defFile.close()

main()

## EOF
