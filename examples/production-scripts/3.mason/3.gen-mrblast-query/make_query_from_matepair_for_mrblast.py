#!/usr/bin/env python

import HTSeq
import itertools
import sys
#import gzip


if __name__ == '__main__':
     
    if len(sys.argv) != 4:
        print "python xxxx.py fileName-1 filename-2 outFile"
        sys.exit(1)  
    fileName1 = sys.argv[1]
    fileName2 = sys.argv[2]
    outFileName = sys.argv[3]
    
    oFile = open(outFileName, 'w')
    
    ###
    ###
    ###
    qFile = open(fileName1, 'r')  
    #qFile = gzip.open(fileName1, 'rb')  
    #f = gzip.open('/home/joe/file.txt.gz', 'rb')
    #file_content = f.read()
  
    numReads = 0
    for line in qFile:            
        if line[0] == '>':
            numReads += 1
            if numReads % 1000 == 0:
                print numReads
            s = line.split()
                
            #>test_illu.fasta.000000098/0
            #contig=seq2
            #haplotype=0
            #length=100
            #orig_begin=2837
            #orig_end=2937
            #haplotype_infix=TCCTGGTTGTAGCTAACTAACTTCAGAACACCAACTTATACCATAATATATATTTTAAAGGACCAGACCAGCTTTCAAAAAGAAAATTGTTAAAGAGAGC
            #edit_string=MMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM
            #strand=forward
    
            origDefLine = s[1].split("=")[-1]
            #readId = s[0].split(".")[-1]
            pairId = s[0].split("/")[-1]
            readId = (s[0].split("/")[-2]).split(".")[-1]
            readLen = int(s[3].split("=")[-1])
            origBegin = int(s[4].split("=")[-1])
            origEnd = int(s[5].split("=")[-1])
            strand = s[8].split("=")[-1]
            newDef = ">" + origDefLine + "_" + str(int(readId)) + "_" + pairId + "_" + \
                     strand[0] + "_" + str(readLen) + "_" + \
                     str(origBegin) + "_" + str(origEnd) + "_" + \
                     str(origBegin) + "_" + str(origBegin+readLen) + "\n"
            oFile.write(newDef)
        else:
            oFile.write(line)
    qFile.close()
    print "num reads in 1 = ",numReads
    
    ###
    ###
    ###
    qFile = open(fileName2, 'r')    
    for line in qFile:            
        if line[0] == '>':
            numReads += 1
            if numReads % 1000 == 0:
                print numReads
            s = line.split()
            origDefLine = s[1].split("=")[-1]
            pairId = s[0].split("/")[-1]
            readId = (s[0].split("/")[-2]).split(".")[-1]
            readLen = int(s[3].split("=")[-1])
            origBegin = int(s[4].split("=")[-1])
            origEnd = int(s[5].split("=")[-1])
            strand = s[8].split("=")[-1]
            newDef = ">" + origDefLine + "_" + str(int(readId)) + "_" + pairId + "_" + \
                     strand[0] + "_" + str(readLen) + "_" + \
                     str(origBegin) + "_" + str(origEnd) + "_" + \
                     str(origBegin) + "_" + str(origBegin+readLen) + "\n"
            oFile.write(newDef)
        else:
            oFile.write(line)
    qFile.close()
    print "total num reads two files = ",numReads
        
    
    oFile.close()
    


## EOF
