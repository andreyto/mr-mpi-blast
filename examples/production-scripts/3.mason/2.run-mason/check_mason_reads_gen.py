#!/usr/bin/env python

#import HTSeq
#import itertools
import sys
import os
import subprocess
import numpy 

path=str(sys.argv[1])
dirList=os.listdir(path)

path2=str(sys.argv[2])
dirList2=os.listdir(path2)

#f_taxid_2_3 = open("/home/ssul/work/mrblast-build/data/basic_test/microbial_real/collect-taxid/2.uniq_taxids_from_2_3/COLLECTED_UNIQ_TAXIDS_FROM_2_3.txt", "r")
taxids = {}
#for l in f_taxid_2_3:
    ##taxids.append(int(l))
    #taxids[l.strip()] = 0
#print len(taxids)
#f_taxid_2_3.close()
 
for fname in dirList:
    #print i, fname
    
    taxid_from_file = fname.split(".")[0]
    #print taxid_from_file
    taxids[taxid_from_file] = 0
    #if taxid_from_file in taxids:
        #print taxid_from_file
        #taxids[taxid_from_file]+=1
        
print len(taxids)


for fname2 in dirList2:
    taxid_from_file2 = fname2.split(".")[0]
    if taxid_from_file2 in taxids:
        taxids[taxid_from_file2]+=1
        
#print taxids
#of = open("check_mason_reads_gen.py.out", "w")
for k, v in taxids.iteritems():
    #print k, v
    if v < 2:
        pass
        #print k, v
        #of.write(str(k)+"\n")
    elif v == 2:
        print k, v
        cmd1 = "mv ./100x/" + str(k) + "*.fasta ./100x-for-2-3/"
        print cmd1
        p = subprocess.Popen(cmd1, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
        stdout_value = p.communicate()[0]
        s = stdout_value.split()
        for i in range(len(s)):
            print str(i) + " " + str(s[i])
#of.close()
