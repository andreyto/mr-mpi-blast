import sys
from Bio import SeqIO
from Bio.Seq import Seq
import os,random
#import numpy as npy

#def chunks(l, n):
    #return [l[i:i+n] for i in range(0, len(l), n)]


if len(sys.argv) != 4:
    print "USAGE: python *.py in_fasta_file chunksize overlapBP\n"
    sys.exit(0)

seqFileName = sys.argv[1]       ## SEQ FILE NAME
chunkSize = int(sys.argv[2])    ## CHUNK SIZE
overlapBP = int(sys.argv[3])   ## OVERLAP BP
seqFile = open(seqFileName, "r")

listcount = 0
cid = 0         ## CHUNK UNIQUE ID
sid = 0         ## seq id
seqLen = 0

for cur_record in SeqIO.parse(seqFile, "fasta") :

    sid += 1  
    origGI = cur_record.id.split()[0][0:] 
    seqLen += len(cur_record.seq)
    start = 0
        
    ##
    ## 
    ##
    while True:
        
        end = start + chunkSize
        cid += 1
        #print "cid=%d, start=%d, end=%d, cnt=%d, seqLen=%d, sid=%d" % (cid, start, end, cnt, seqLen, sid)
        
        if start == 0 and end < seqLen:
            pass
            #newHeader = ">"+origGI+"_"+str(cid)+"_"+"0"+"_"+str(start)+"_"+str(end)+"_"+str(start)+"_"+str(end)
            
        elif start == 0 and end >= seqLen:
            #newHeader =  ">"+origGI+"_"+str(cid)+"_"+"1"+"_"+str(start)+"_"+str(seqLen)+"_"+str(start)+"_"+str(seqLen)
            break

        elif start > 0 and end < seqLen:
            #newHeader = ">"+origGI+"_"+str(cid)+"_"+"2"+"_"+str(start)+"_"+str(end)+"_"+str(start)+"_"+str(end)
            pass

        elif start > 0 and end >= seqLen:
            #newHeader =  ">"+origGI+"_"+str(cid)+"_"+"3"+"_"+str(start)+"_"+str(seqLen)+"_"+str(start)+"_"+str(seqLen)
            break
            
        else:
            print "What?", seqLen, start, end, cur_record.id
            sys.exit(0)
        
        start += overlapBP
        
seqFile.close()
print "# input sequences = ", sid
print "chunk size = ", chunkSize, "bp"
print "overlap = ", overlapBP, "bp"
print "Total %d chunks are saved." % (cid)
print "Total seq len = ", seqLen


