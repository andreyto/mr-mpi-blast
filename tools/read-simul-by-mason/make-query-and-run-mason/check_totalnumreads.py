#!/usr/bin/env python

#import HTSeq
#import itertools
import sys
import os
import subprocess
import numpy


run_opt = int(sys.argv[2])
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
dict_taxid_numreads_illu = {}
dict_taxid_numreads_454 = {}
num_reads_check_file = open("/media/disk/SU-estimate2.csv", "r")
for l in num_reads_check_file:
    taxid = (l.split(",")[1]).split(".")[0]
    numreads_illu = l.split(",")[7]
    numreads_454 = l.split(",")[8]
    #print taxid, numreads_illu, numreads_454
    dict_taxid_numreads_illu[int(taxid)] = int(numreads_illu)
    dict_taxid_numreads_454[int(taxid)] = int(numreads_454)
print len(dict_taxid_numreads_illu), len(dict_taxid_numreads_454)


i=0
for fname in dirList2:
    #print  fname, os.path.isdir(fname)

    if run_opt == 0 and os.path.isdir(fname) == True:
        dirList3=os.listdir(path2 + "/" + fname)
        #print dirList3
        for fname2 in dirList3:
            if (fname2.find("TOTAL") > -1):
                #print fname, fname2
                if int(fname) in dict_taxid_numreads_illu:
                    num_generated_reads = int(fname2.split("_")[3])
                    if num_generated_reads == dict_taxid_numreads_illu[int(fname)] or \
                       num_generated_reads - 1 == dict_taxid_numreads_illu[int(fname)]:
                        print "OK", dict_taxid_numreads_illu[int(fname)], "==>", num_generated_reads
                        i+=1
                        #pass
                    else:
                        print "ERR: ", fname, "should be", dict_taxid_numreads_illu[int(fname)], " but I have ", num_generated_reads
                    
                #elif int(fname.split("-")[4]) in dict_taxid_numreads_454:
                    #num_generated_reads = int(fname2.split("_")[3])
                    #if num_generated_reads == dict_taxid_numreads_454[int(fname)] or \
                       #num_generated_reads - 1 == dict_taxid_numreads_454[int(fname)]:
                        ##print "OK", dict_taxid_numreads_454[int(fname)], "==>", num_generated_reads
                        #pass
                    #else:
                        #print "ERR: ", fname, "should be", dict_taxid_numreads_454[int(fname)], " but I have ", num_generated_reads
                    
                        
                else:
                    print "?"
                    
    elif run_opt == 1 and os.path.isdir(fname) == True:
        dirList3=os.listdir(path2 + "/" + fname)
        #print dirList3
        for fname2 in dirList3:
            if (fname2.find("TOTAL") > -1):
                #print fname, fname2
                #if int(fname) in dict_taxid_numreads_illu:
                    #num_generated_reads = int(fname2.split("_")[3])
                    #if num_generated_reads == dict_taxid_numreads_illu[int(fname)] or \
                       #num_generated_reads - 1 == dict_taxid_numreads_illu[int(fname)]:
                        ##print "OK", dict_taxid_numreads_illu[int(fname)], "==>", num_generated_reads
                        #pass
                    #else:
                        #print "ERR: ", fname, "should be", dict_taxid_numreads_illu[int(fname)], " but I have ", num_generated_reads
                    
                ttt = int(fname.split("-")[-1]) 
                #print ttt
                if ttt in dict_taxid_numreads_454:
                    num_generated_reads = int(fname2.split("_")[3])
                    if num_generated_reads == dict_taxid_numreads_454[ttt] or \
                       num_generated_reads - 1 == dict_taxid_numreads_454[ttt] or \
                       num_generated_reads + 1 == dict_taxid_numreads_454[ttt]:
                        print "OK", dict_taxid_numreads_454[ttt], "==>", num_generated_reads
                        i+=1
                        #pass
                    else:
                        print "ERR: ", fname, "should be", dict_taxid_numreads_454[ttt], " but I have ", num_generated_reads
                else:
                    print "?"
                    
    elif run_opt == 2 and os.path.isdir(fname) == True:
		dirList3=os.listdir(path2 + "/" + fname)
        #print dirList3
		for fname2 in dirList3:
			if (fname2.find("TOTAL") > -1):
				ttt = int(fname.split("-")[-1]) 
				#print ttt
				if ttt in dict_taxid_numreads_illu:
					num_generated_reads = int(fname2.split("_")[-1])
					if num_generated_reads == dict_taxid_numreads_illu[ttt] or \
					   num_generated_reads - 1 == dict_taxid_numreads_illu[ttt] or \
					   num_generated_reads + 1 == dict_taxid_numreads_illu[ttt]:
						print "OK", dict_taxid_numreads_illu[ttt], "==>", num_generated_reads
						i+=1
						#pass
					else:
						print "ERR: ", fname, "should be", dict_taxid_numreads_illu[ttt], " but I have ", num_generated_reads
				else:
					print "?"

print "success = ", i
## EOF
