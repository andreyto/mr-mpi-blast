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
        for cur_record in SeqIO.parse(inFileName, "fasta"):
            if (len(cur_record) < 1000000):
                #print cur_record.id
                #print cur_record.seq
                #print len(cur_record) #len(cur_record.seq) seq length
                header1 = cur_record.id
                seq1 = cur_record.seq.upper()
            if (len(cur_record) > 5000000):
                #print cur_record.id
                #print repr(cur_record.seq)
                #print cur_record.seq
                #print len(cur_record) #len(cur_record.seq) seq length
                header2 = cur_record.id
                seq2 = cur_record.seq.upper()
                break
        print header1
        print seq1
        print header2
        print seq2
        
        
