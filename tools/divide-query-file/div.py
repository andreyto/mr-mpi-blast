#!/usr/bin/env python

import sys
import os
  
if len(sys.argv) != 4:
    print "USAGE: python div.py fasta_file numLinesPerFile totalLine"
    sys.exit(0)
    
seqsubFileName = sys.argv[1]  
numLine     = int(sys.argv[2])   
totalLine   = int(sys.argv[3])
#totalLine   = 49000894

fileNum = 0
#line = 1

numDir = int(totalLine / numLine)
remain = totalLine % numLine
print "numDir, remaining lines = ",numDir, remain

seqFile = open(seqsubFileName, "r")

while fileNum < numDir:
    dirName = "div"+str(fileNum)
    cmd = "mkdir " + dirName
    os.system(cmd)
    subFileName = dirName + "/data.fasta"
    print subFileName
    subFile = open(subFileName, "w")
    li = 0
    while li < numLine:
        defLine = seqFile.readline()
        seqLine = seqFile.readline()
        li += 2        
        subFile.write(defLine)
        subFile.write(seqLine)    
    subFile.close()
    fileNum += 1
    #cmd = "python /work/01471/ssul/work/distros3/ncbi_cxx/ncbi_cxx--Jun_15_2010/src/app/mr-mpi-blast/tools/seqindexer/seqindexer.py  " + subFileName + " " + subFileName + ".idx2"
    cmd = "python ../../../../../src/tools/seqindexer/seqindexer.py " + subFileName + " " + subFileName + ".idx2"
    os.system(cmd)
    cmd = "cp ./mrblast.ini " + dirName
    os.system(cmd)
    
if remain > 0:
    line = 0
    dirName = "div"+str(fileNum)
    cmd = "mkdir " + dirName
    os.system(cmd)
    subFileName = dirName + "/data.fasta"
    print subFileName
    subFile = open(subFileName, "w")
    while line < remain:
        defLine = seqFile.readline()
        seqLine = seqFile.readline()
        line += 2        
        subFile.write(defLine)
        subFile.write(seqLine)    
    subFile.close()
    #cmd = "python /work/01471/ssul/work/distros3/ncbi_cxx/ncbi_cxx--Jun_15_2010/src/app/mr-mpi-blast/tools/seqindexer/seqindexer.py  " + subFileName + " " + subFileName + ".idx2"
    cmd = "python ../../../../../src/tools/seqindexer/seqindexer.py " + subFileName + " " + subFileName + ".idx2"
    os.system(cmd)
    cmd = "cp ./mrblast.ini " + dirName
    os.system(cmd)

seqFile.close()

### EOF
