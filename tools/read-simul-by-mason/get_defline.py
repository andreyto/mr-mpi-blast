import sys
from Bio import SeqIO
from Bio.Seq import Seq
import os,random
 
if __name__ == '__main__':
    
    if len(sys.argv) != 2:
        sys.exit(0)
     
    seqFileName = sys.argv[1]           ## SEQ FILE NAME
    #outFileName = sys.argv[2]
    seqFile = open(seqFileName, "r")
    
    #outFileName = seqFileName + ".filtered.fasta"
    #outFile = open(outFileName, "w")
    
    total_len = 0
    nrec = 0
    for cur_record in SeqIO.parse(seqFile, "fasta") :
        #origGI = cur_record.id.split()[0][0:] 
        seqLen = len(cur_record.seq)
        #print origGI, segLen
        #print cur_record.description
        #if (cur_record.description.find("plasmid") > -1):
        #if (cur_record.description.find("plasmid") > -1) :
            #pass
        #elif (cur_record.description.find("scaffold") > -1) : 
            #pass
        #elif (cur_record.description.find("|ref|") <= -1):
            #pass 
        #else:          
            #print cur_record.description
            #print cur_record.seq
            #outFile.write(">"+str(cur_record.description)+"\n")
            #outFile.write(str(cur_record.seq)+"\n")
        print ">"+str(cur_record.description), "length=", seqLen
        total_len += seqLen
        nrec+=1
        
    seqFile.close()
    #outFile.close()
    
    print "num gis = ", nrec
    print "total len = ", total_len
    
    
## EOF
