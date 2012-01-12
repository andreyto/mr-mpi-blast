#!/usr/bin/env python

#import HTSeq
#import itertools
import sys
import os
import subprocess
import numpy

#path=str(sys.argv[1])
path2=str(sys.argv[1])
#path="./illu_100x_reads_2_3/"
#dirList=os.listdir(path)

#path2=str(sys.argv[2])
#path2="./ranger_illu_100x_reads_2_3_mrblast/"
#path2="./"
dirList2=os.listdir(path2)

#f_taxid_2_3 = open("/home/ssul/work/mrblast-build/data/basic_test/microbial_real/collect-taxid/2.uniq_taxids_from_2_3/COLLECTED_UNIQ_TAXIDS_FROM_2_3.txt", "r")
taxids = {}
#for l in f_taxid_2_3:
    ##taxids.append(int(l))
    #taxids[l.strip()] = 0
#print len(taxids)
#f_taxid_2_3.close()

for fname in dirList2:
    #print  fname
    #taxid_from_file = fname.split(".")[0]
    #print taxid_from_file
    #taxids[taxid_from_file] = 0
    #if taxid_from_file in taxids:
        #print taxid_from_file
        #taxids[taxid_from_file]+=1

    if os.path.exists(path2 + fname +"/output-workitems.txt"):
        if os.path.exists(path2 + fname +"/SUCCESS"): #and os.path.exists(path2 + fname +"/output-hits-.txt.bin"):
            print fname
        else:
            print fname, "no success"
