import sys
from Bio import SeqIO
from Bio.Seq import Seq
import os,random
#import numpy as npy

#def chunks(l, n):
    #return [l[i:i+n] for i in range(0, len(l), n)]


if len(sys.argv) != 5:
    print "USAGE: python splitter.py fasta_file chunksize overlapBP outfile"
    print "       - fasta_file: input file"
    print "       - chunksize: length (bp) of splitted sequences"
    print "       - overlapBP: length (bp) of overlap"
    print "       - outfile: output file\n"
    sys.exit(0)

seqFileName = sys.argv[1]       ## SEQ FILE NAME
chunkSize = int(sys.argv[2])    ## CHUNK SIZE
overlapBP = int(sys.argv[3])   ## OVERLAP BP
outFileName = sys.argv[4]
seqFile = open(seqFileName, "r")
outFile = open(outFileName, "w")

listcount = 0
cid = 0         ## CHUNK UNIQUE ID
sid = 0         ## seq id

for cur_record in SeqIO.parse(seqFile, "fasta") :

    sid += 1  
    origGI = cur_record.id.split()[0][0:] 
    #print origGI
    seqLen = len(cur_record.seq)
    start = 0
    cnt = 0;
    
    ##
    ## 
    ##
    while True:
        
        end = start + chunkSize
        cid += 1
        #print "cid=%d, start=%d, end=%d, cnt=%d, seqLen=%d, sid=%d" % (cid, start, end, cnt, seqLen, sid)
        
        if start == 0 and end < seqLen:
            newHeader = ">"+origGI+"_"+str(start)+"_"+str(end)
            print newHeader, sid, cid
            outFile.write(newHeader+"\n")
            upper = cur_record.seq[start:end-overlapBP].upper()
            lower2 = cur_record.seq[end-overlapBP:end].lower()
            outFile.write((upper+lower2).tostring()+"\n")
                        
        elif start == 0 and end >= seqLen:
            newHeader =  ">"+origGI+"_"+str(start)+"_"+str(seqLen) 
            print newHeader, sid, cid
            outFile.write(newHeader+"\n")
            outFile.write(cur_record.seq[start:seqLen].tostring()+"\n")
 
            break

        elif start > 0 and end < seqLen:
            newHeader = ">"+origGI+"_"+str(start)+"_"+str(end) 
            print newHeader, sid, cid
            outFile.write(newHeader+"\n")
            lower1 = cur_record.seq[start:start+overlapBP].lower()
            upper = cur_record.seq[start+overlapBP:end-overlapBP].upper()
            lower2 = cur_record.seq[end-overlapBP:end].lower()
            outFile.write((lower1+upper+lower2).tostring()+"\n")
 
        elif start > 0 and end >= seqLen:
            newHeader =  ">"+origGI+"_"+str(start)+"_"+str(seqLen) 
            print newHeader, sid, cid
            outFile.write(newHeader+"\n")
            if (seqLen - start > overlapBP):
                lower = cur_record.seq[start:start+overlapBP].lower()
                upper = cur_record.seq[start+overlapBP:seqLen].upper()
                outFile.write((lower+upper).tostring()+"\n")
            else:
                outFile.write(cur_record.seq[start:seqLen].tostring()+"\n")
 
            break
            
        else:
            print "What?", seqLen, start, end, cur_record.id
            sys.exit(0)
        
        cnt += 1
        start = end - overlapBP
        
outFile.close()
seqFile.close()
print "# input sequences = ", sid
print "chunk size = ", chunkSize, "bp"
print "overlap = ", overlapBP, "bp"
print "Total %d chunks are saved in %s." % (cid, outFileName)


