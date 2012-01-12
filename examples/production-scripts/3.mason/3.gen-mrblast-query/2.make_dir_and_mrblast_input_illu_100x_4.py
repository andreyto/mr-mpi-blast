#!/usr/bin/env python

#import HTSeq
#import itertools
import sys
import os
import subprocess
import numpy 

#path=str(sys.argv[1])
path="./illu_100x_reads_4/"
dirList=os.listdir(path)

#path2=str(sys.argv[2])
path2="./ranger_illu_100x_reads_4_mrblast/"
#dirList2=os.listdir(path2)
 

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

##
## Make each taxid folder 
##
for k in taxids.iterkeys():
    outputpath = path2 + "illu-100x-4-" + str(k)
    cmd1 = "mkdir " + outputpath
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
    outputpath = path2 + "illu-100x-4-" + str(k)
    output = outputpath + "/" + str(k) + ".gis.fasta_illu_100x_reads_mrblast.fasta"
    
    if os.path.exists(outputpath + "/" + str(k) + ".gis.fasta_illu_100x_reads_mrblast.fasta.idx") and \
       os.path.exists(outputpath + "/" + "mrblast.ini") and \
       os.path.exists(outputpath + "/" + str(k) + ".gis.fasta_illu_100x_reads_mrblast.fasta") and \
       os.path.getsize(outputpath + "/" + str(k) + ".gis.fasta_illu_100x_reads_mrblast.fasta.idx") > 0 and \
       os.path.getsize(outputpath + "/" + "mrblast.ini") > 0 and \
       os.path.getsize(outputpath + "/" + str(k) + ".gis.fasta_illu_100x_reads_mrblast.fasta") > 0:
        print str(k), "is all done!!!"
    else:    
        cmd1 = "python " + make_mrblast_read_exe + " " + input1 + " " + input2 + " " + output
        print cmd1
        p = subprocess.Popen(cmd1, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
        stdout_value = p.communicate()[0]
        s = stdout_value.split()
        
        for i in range(len(s)):
            print str(i) + " " + str(s[i])
        print "total num reads = ", s[-1]
            
        ## make index file
        index_f_name = str(k) + ".gis.fasta_illu_100x_reads_mrblast.fasta.idx"
        cmd1 = "python " + make_index_exe + outputpath + "/" + str(k) + ".gis.fasta_illu_100x_reads_mrblast.fasta" + " " + outputpath + "/" + index_f_name
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
        cmd1 = "cp mrblast_part.ini " + outputpath + "/" + "mrblast.ini"
        print cmd1    
        p = subprocess.Popen(cmd1, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
        stdout_value = p.communicate()[0]
        s = stdout_value.split()
        for i in range(len(s)):
            print str(i) + " " + str(s[i])
        
        #if os.path.exists(outputpath + "/" + str(k) + ".gis.fasta_illu_100x_reads_mrblast.fasta.idx") and \
           #os.path.exists(outputpath + "/" + "mrblast.ini") and \
           #os.path.exists(outputpath + "/" + str(k) + ".gis.fasta_illu_100x_reads_mrblast.fasta"):
        if os.path.exists(outputpath + "/" + str(k) + ".gis.fasta_illu_100x_reads_mrblast.fasta.idx") and \
           os.path.exists(outputpath + "/" + "mrblast.ini") and \
           os.path.exists(outputpath + "/" + str(k) + ".gis.fasta_illu_100x_reads_mrblast.fasta") and \
           os.path.getsize(outputpath + "/" + str(k) + ".gis.fasta_illu_100x_reads_mrblast.fasta.idx") > 0 and \
           os.path.getsize(outputpath + "/" + "mrblast.ini") > 0 and \
           os.path.getsize(outputpath + "/" + str(k) + ".gis.fasta_illu_100x_reads_mrblast.fasta") > 0:
           
            f = open(outputpath + "/" + "mrblast.ini", "a")
            f.write("QUERYFILENAME = " + str(k) + ".gis.fasta_illu_100x_reads_mrblast.fasta"+"\n")
            f.write("INDEXFILENAME = " + index_f_name + "\n")
            f.close()
            
            cmd1 = "touch " + outputpath + "/TOTAL_NUM_READS_" + total_num_reads
            print cmd1    
            p = subprocess.Popen(cmd1, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
            stdout_value = p.communicate()[0]
            s = stdout_value.split()
            for i in range(len(s)):
                print str(i) + " " + str(s[i])
        else:
            print "ERROR: failed to make mrblast.ini and index for taxid = ", str(k)
        
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
