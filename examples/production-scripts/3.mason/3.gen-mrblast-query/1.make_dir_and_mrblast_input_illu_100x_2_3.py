#!/usr/bin/env python

#import HTSeq
#import itertools
import sys
import os
import subprocess
import numpy 

#path=str(sys.argv[1])
path="./illu_100x_reads_2_3/"
dirList=os.listdir(path)

#path2=str(sys.argv[2])
path2="./ranger_illu_100x_reads_2_3_mrblast/"
#dirList2=os.listdir(path2)

taxids = {}
for fname in dirList:
    taxid_from_file = fname.split(".")[0]
    taxids[taxid_from_file] = 0
print len(taxids)

#ERR:  314608 should be 4351589  but I have  3780487
#ERR:  314253 should be 6304713  but I have  3136710
#ERR:  247633 should be 3925629  but I have  3872896
#ERR:  270374 should be 4894744  but I have  7021
#ERR:  232348 should be 2686395  but I have  1347615
#taxids = {}
#taxids["314608"] = 0
#taxids["314253"] = 0
#taxids["247633"] = 0
#taxids["270374"] = 0
#taxids["232348"] = 0
#print len(taxids)

##
## Make each taxid folder in illu_100x_reads_2_3_mrblast
##
for k in taxids.iterkeys():
    cmd1 = "mkdir " + path2 + str(k)
    print cmd1
    p = subprocess.Popen(cmd1, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    stdout_value = p.communicate()[0]
    s = stdout_value.split()
    for i in range(len(s)):
        print str(i) + " " + str(s[i])
        
##
## Make mrblast input
##
make_mrblast_read_exe = "/home/ssul/work/mrblast-build/src/tools/read-simul-by-mason/make_query_from_matepair_for_mrblast.py "
make_index_exe = "/home/ssul/work/mrblast-build/src/tools/seqindexer/seqindexer.py "

for k in taxids.iterkeys():
    input1 = path + str(k) + ".gis.fasta_illu_100x_reads_1.fasta"
    input2 = path + str(k) + ".gis.fasta_illu_100x_reads_2.fasta"
    output = path2 + str(k) + "/" + str(k) + ".gis.fasta_illu_100x_reads_mrblast.fasta"
    
    #if os.path.exists(path2 + str(k) + "/mrblast.ini"):
    if os.path.exists(path2 + str(k) + "/" + str(k) + ".gis.fasta_illu_100x_reads_mrblast.fasta.idx") and \
       os.path.exists(path2 + str(k) + "/" + "mrblast.ini") and \
       os.path.exists(path2 + str(k) + "/" + str(k) + ".gis.fasta_illu_100x_reads_mrblast.fasta"):
        print str(k), "is done"
    else:    
        cmd1 = "python " + make_mrblast_read_exe + " " + input1 + " " + input2 + " " + output
        print cmd1
        p = subprocess.Popen(cmd1, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
        stdout_value = p.communicate()[0]
        s = stdout_value.split()
        print "total num reads = ", s[-1]
        #for i in range(len(s)):
            #print str(i) + " " + str(s[i])
            
        ## make index file
        index_f_name = str(k) + ".gis.fasta_illu_100x_reads_mrblast.fasta.idx"
        cmd1 = "python " + make_index_exe + path2 + str(k) + "/" + str(k) + ".gis.fasta_illu_100x_reads_mrblast.fasta" + " " + path2 + str(k) + "/" + index_f_name
        print cmd1
        p = subprocess.Popen(cmd1, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
        stdout_value = p.communicate()[0]
        s = stdout_value.split()
        print "total num reads = ", s[-1]
        total_num_reads = s[-1]
        #for i in range(len(s)):
            #print str(i) + " " + str(s[i])
            
        ## make mrblast.ini
        #QUERYFILENAME       
        #INDEXFILENAME
        cmd1 = "cp mrblast_part.ini " + path2 + str(k) + "/" + "mrblast.ini"
        print cmd1    
        p = subprocess.Popen(cmd1, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
        stdout_value = p.communicate()[0]
        s = stdout_value.split()
        for i in range(len(s)):
            print str(i) + " " + str(s[i])
        
        if os.path.exists(path2 + str(k) + "/" + str(k) + ".gis.fasta_illu_100x_reads_mrblast.fasta.idx") and \
           os.path.exists(path2 + str(k) + "/" + "mrblast.ini") and \
           os.path.exists(path2 + str(k) + "/" + str(k) + ".gis.fasta_illu_100x_reads_mrblast.fasta"):
           
            f = open(path2 + str(k) + "/" + "mrblast.ini", "a")
            f.write("QUERYFILENAME = " + str(k) + ".gis.fasta_illu_100x_reads_mrblast.fasta"+"\n")
            f.write("INDEXFILENAME = " + index_f_name + "\n")
            f.close()
            
            cmd1 = "touch " + path2 + str(k) + "/TOTAL_NUM_READS_" + total_num_reads
            print cmd1    
            p = subprocess.Popen(cmd1, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
            stdout_value = p.communicate()[0]
            s = stdout_value.split()
            for i in range(len(s)):
                print str(i) + " " + str(s[i])
        else:
            print "ERROR: failed to make index ", str(k)
        
        print "##################################"

#for fname2 in dirList2:
    #taxid_from_file2 = fname2.split(".")[0]
    #if taxid_from_file2 in taxids:
        #taxids[taxid_from_file2]+=1
        
##print taxids
##of = open("check_mason_reads_gen.py.out", "w")
#for k, v in taxids.iteritems():
    ##print k, v
    #if v < 2:
        #pass
        ##print k, v
        ##of.write(str(k)+"\n")
    #elif v == 2:
        #print k, v
        #cmd1 = "mv ./100x/" + str(k) + "*.fasta ./100x-for-2-3/"
        #print cmd1
        #p = subprocess.Popen(cmd1, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
        #stdout_value = p.communicate()[0]
        #s = stdout_value.split()
        #for i in range(len(s)):
            #print str(i) + " " + str(s[i])
##of.close()
