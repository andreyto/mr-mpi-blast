import sys
from Bio import SeqIO
from Bio.Seq import Seq
import os,random
#import numpy as npy

#def chunks(l, n):
    #return [l[i:i+n] for i in range(0, len(l), n)]


if len(sys.argv) != 6:
    print "USAGE: python splitter.py fasta_file chunksize lowerPartLen outfile"
    print "       - fasta_file: input file"
    print "       - upperPartLen: length (bp) of uppercase sequences"
    print "       - lowerPartLen: length (bp) of lowercase overlap"
    print "       - upperPartOverlap: length (bp) of uppercase part overlap"
    print "       - outfile: output file\n"
    sys.exit(0)

###
###
###

seqFileName = sys.argv[1]       ## SEQ FILE NAME
upperPartLen = int(sys.argv[2]) ## upperPart length BP
lowerPartLen = int(sys.argv[3])   ## lowerPart length BP
upperPartOverlap = int(sys.argv[4])   ## upperPart overlap length BP
outFileName = sys.argv[5]
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
    cutStart = 0
    upperStart = 0
    cnt = 0;

    while True:
        #end = cutStart + upperPartLen + (lowerPartLen*2)
        cid += 1
        
        if cutStart == 0 and upperStart == 0:
            if seqLen <= upperPartLen + lowerPartLen:
                newHeader = ">"+origGI+"_"+str(cid)+"_"+"0"+"_"+str(cutStart)+"_"+str(seqLen)+"_"+str(upperStart)+"_"+str(seqLen)
                print newHeader
                outFile.write(newHeader+"\n")
                upper = cur_record.seq[cutStart:seqLen].upper()
                outFile.write((upper).tostring()+"\n")
                break
            elif seqLen > upperPartLen + lowerPartLen:
                newHeader = ">"+origGI+"_"+str(cid)+"_"+"1"+"_"+str(cutStart)+"_"+str(upperPartLen + lowerPartLen)+"_"+str(upperStart)+"_"+str(upperPartLen)
                print newHeader
                outFile.write(newHeader+"\n")
                upper = cur_record.seq[cutStart:upperPartLen].upper()
                lower2 = cur_record.seq[upperPartLen:upperPartLen + lowerPartLen].lower()
                outFile.write((upper+lower2).tostring()+"\n")
                cutStart = 0
                upperStart += upperPartOverlap                
     
        elif cutStart == 0 and upperStart == upperPartOverlap:
            if seqLen <= upperStart + upperPartLen + lowerPartLen:
                newHeader = ">"+origGI+"_"+str(cid)+"_"+"2"+"_"+str(cutStart)+"_"+str(seqLen)+"_"+str(upperStart)+"_"+str(seqLen)
                print newHeader
                outFile.write(newHeader+"\n")
                lower1 = cur_record.seq[cutStart:upperStart].lower()
                upper = cur_record.seq[upperStart:seqLen].upper()
                outFile.write((lower1+upper).tostring()+"\n")
                break
            elif seqLen > upperStart + upperPartLen + lowerPartLen:
                newHeader = ">"+origGI+"_"+str(cid)+"_"+"3"+"_"+str(cutStart)+"_"+str(upperStart + upperPartLen + lowerPartLen)+"_"+str(upperStart)+"_"+str(upperStart+upperPartLen)
                print newHeader
                outFile.write(newHeader+"\n")
                lower1 = cur_record.seq[cutStart:upperStart].lower()
                upper = cur_record.seq[upperStart:upperStart+upperPartLen].upper()
                lower2 = cur_record.seq[upperStart+upperPartLen:upperStart+upperPartLen+lowerPartLen].lower()
                outFile.write((lower1+upper+lower2).tostring()+"\n")
                cutStart = 0
                upperStart += upperPartOverlap 
                               
        elif cutStart == 0 and upperStart == upperPartOverlap * 2:
            if seqLen <= upperStart + upperPartLen + lowerPartLen:
                newHeader = ">"+origGI+"_"+str(cid)+"_"+"4"+"_"+str(cutStart)+"_"+str(seqLen)+"_"+str(upperStart)+"_"+str(seqLen)
                print newHeader
                outFile.write(newHeader+"\n")
                lower1 = cur_record.seq[cutStart:upperStart].lower()
                upper = cur_record.seq[upperStart:seqLen].upper()
                outFile.write((lower1+upper).tostring()+"\n")
                break
            elif seqLen > upperStart + upperPartLen + lowerPartLen:
                newHeader = ">"+origGI+"_"+str(cid)+"_"+"5"+"_"+str(cutStart)+"_"+str(upperStart + upperPartLen + lowerPartLen)+"_"+str(upperStart)+"_"+str(upperStart + upperPartLen)
                print newHeader
                outFile.write(newHeader+"\n")
                lower1 = cur_record.seq[cutStart:upperStart].lower()
                upper = cur_record.seq[upperStart:upperStart+upperPartLen].upper()
                lower2 = cur_record.seq[upperStart+upperPartLen:upperStart+upperPartLen+lowerPartLen].lower()
                outFile.write((lower1+upper+lower2).tostring()+"\n")
                cutStart += upperPartOverlap
                upperStart += upperPartOverlap 
            
        elif cutStart > 0 and upperStart > upperPartOverlap:
            if seqLen <= upperStart + upperPartLen + lowerPartLen:
                newHeader = ">"+origGI+"_"+str(cid)+"_"+"6"+"_"+str(cutStart)+"_"+str(seqLen)+"_"+str(upperStart)+"_"+str(seqLen)
                print newHeader
                outFile.write(newHeader+"\n")
                lower1 = cur_record.seq[cutStart:upperStart].lower()
                upper = cur_record.seq[upperStart:seqLen].upper()
                outFile.write((lower1+upper).tostring()+"\n")
                break
            else:
                newHeader = ">"+origGI+"_"+str(cid)+"_"+"7"+"_"+str(cutStart)+"_"+str(upperStart + upperPartLen + lowerPartLen)+"_"+str(upperStart)+"_"+str(upperStart + upperPartLen)
                print newHeader
                outFile.write(newHeader+"\n")
                lower1 = cur_record.seq[cutStart:upperStart].lower()
                upper = cur_record.seq[upperStart:upperStart+upperPartLen].upper()
                lower2 = cur_record.seq[upperStart+upperPartLen:upperStart+upperPartLen+lowerPartLen].lower()
                outFile.write((lower1+upper+lower2).tostring()+"\n")
                cutStart += upperPartOverlap
                upperStart += upperPartOverlap    
                
        #if cutStart == 0:
            #if seqLen <= upperPartLen:
                #newHeader = ">"+origGI+"_"+str(cid)+"_"+"0"+"_"+str(cutStart)+"_"+str(seqLen)+"_"+str(cutStart)+"_"+str(seqLen)
                #print newHeader
                #outFile.write(newHeader+"\n")
                #upper = cur_record.seq[cutStart:seqLen].upper()
                #outFile.write((upper).tostring()+"\n")
                #break
            #elif seqLen > upperPartLen and seqLen <= upperPartLen + lowerPartLen:
                #newHeader = ">"+origGI+"_"+str(cid)+"_"+"1"+"_"+str(cutStart)+"_"+str(seqLen)+"_"+str(cutStart)+"_"+str(upperPartLen)
                #print newHeader
                #outFile.write(newHeader+"\n")
                #upper = cur_record.seq[cutStart:upperPartLen].upper()
                #lower2 = cur_record.seq[upperPartLen:seqLen].lower()
                #outFile.write((upper+lower2).tostring()+"\n")
                #break
            #elif seqLen > upperPartLen + lowerPartLen:
                #newHeader = ">"+origGI+"_"+str(cid)+"_"+"2"+"_"+str(cutStart)+"_"+str(upperPartLen + lowerPartLen)+"_"+str(cutStart)+"_"+str(upperPartLen)
                #print newHeader
                #outFile.write(newHeader+"\n")
                #upper = cur_record.seq[cutStart:upperPartLen].upper()
                #lower2 = cur_record.seq[upperPartLen:upperPartLen + lowerPartLen].lower()
                #outFile.write((upper+lower2).tostring()+"\n")
            #else:
                #print "what?"
                #sys.exit(0)

        #elif cutStart > 0:
            #if seqLen <= cutStart + upperPartLen:
                #newHeader = ">"+origGI+"_"+str(cid)+"_"+"3"+"_"+str(cutStart - lowerPartLen)+"_"+str(seqLen)+"_"+str(cutStart)+"_"+str(seqLen)
                #print newHeader
                #outFile.write(newHeader+"\n")
                #lower1 = cur_record.seq[cutStart-lowerPartLen:cutStart].lower()
                #upper = cur_record.seq[cutStart:seqLen].upper()
                #outFile.write((lower1+upper).tostring()+"\n")
                #break
            #elif seqLen > cutStart + upperPartLen and seqLen <= cutStart + upperPartLen + lowerPartLen:
                #newHeader = ">"+origGI+"_"+str(cid)+"_"+"4"+"_"+str(cutStart - lowerPartLen)+"_"+str(seqLen)+"_"+str(cutStart)+"_"+str(cutStart+upperPartLen)
                #print newHeader
                #outFile.write(newHeader+"\n")
                #lower1 = cur_record.seq[cutStart-lowerPartLen:cutStart].lower()
                #upper = cur_record.seq[cutStart:cutStart+upperPartLen].upper()
                #lower2 = cur_record.seq[cutStart+upperPartLen:seqLen].lower()
                #outFile.write((lower1+upper+lower2).tostring()+"\n")
            #elif seqLen > cutStart + upperPartLen + lowerPartLen:
                #newHeader = ">"+origGI+"_"+str(cid)+"_"+"5"+"_"+str(cutStart - lowerPartLen)+"_"+str(cutStart + upperPartLen + lowerPartLen)+"_"+str(cutStart)+"_"+str(cutStart+upperPartLen)
                #print newHeader
                #outFile.write(newHeader+"\n")
                #lower1 = cur_record.seq[cutStart-lowerPartLen:cutStart].lower()
                #upper = cur_record.seq[cutStart:cutStart+upperPartLen].upper()
                #lower2 = cur_record.seq[cutStart+upperPartLen:cutStart + upperPartLen + lowerPartLen].lower()
                #outFile.write((lower1+upper+lower2).tostring()+"\n")
            #else:
                #print "what??"
                #sys.exit(0)
                
        else:
            print "What???", seqLen, cutStart, end, cur_record.id
            sys.exit(0)
        
        cnt += 1
        #cutStart = cutStart + upperPartOverlap
        #upperStart = 
        
outFile.close()
seqFile.close()
print "# input sequences = ", sid
#print "chunk size = ", chunkSize, "bp"
print "overlap = ", lowerPartLen, "bp"
print "Total %d chunks are saved in %s." % (cid, outFileName)


