import os,random
import sys
#from random import *
from Bio import SeqIO
#import numpy as npy

if __name__ == "__main__":
    
    inFileName = sys.argv[1] 
    op_case = int(sys.argv[2])
    
    if (op_case == 1):    
        sid = 0
        for cur_record in SeqIO.parse(inFileName, "fasta"):
            print cur_record.id, len(cur_record) 
            sid += 1
        print "OK! Num seqs = %d" % sid
    
    elif (op_case == 2):
        numShortGenome = 0
        for cur_record in SeqIO.parse(inFileName, "fasta"):
            if (len(cur_record) < 3000):
                #print cur_record.id
                #print cur_record.seq
                #print len(cur_record) #len(cur_record.seq) seq length
                header1 = cur_record.id
                seq1 = cur_record.seq.upper()
                print ">"+header1
                print seq1
                numShortGenome += 1
                if (numShortGenome == 10):
                    break
                
    elif (op_case == 3):
        numLongGenome = 0
        for cur_record in SeqIO.parse(inFileName, "fasta"):
            if (len(cur_record) > 5000000):
                #print cur_record.id
                #print repr(cur_record.seq)
                #print cur_record.seq
                #print len(cur_record) #len(cur_record.seq) seq length
                header2 = cur_record.id
                seq2 = cur_record.seq.upper()
                print ">"+header2
                print seq2
                numLongGenome += 1        
                if (numLongGenome == 2):
                    break
                    
    elif (op_case == 4):
        numShortGenome = 0
        for cur_record in SeqIO.parse(inFileName, "fasta"):
            if (len(cur_record) > 10000 and len(cur_record) < 20000):
                #print cur_record.id
                #print cur_record.seq
                #print len(cur_record) #len(cur_record.seq) seq length
                header1 = cur_record.id
                seq1 = cur_record.seq.upper()
                print ">"+header1
                print seq1
                numShortGenome += 1
                if (numShortGenome == 2):
                    break
## EOF        
