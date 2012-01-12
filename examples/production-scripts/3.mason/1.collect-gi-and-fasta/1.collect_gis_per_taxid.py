from MGT.Taxa import *
#from MGT.Svm import *
from MGT.FastaIO import *
#from MGT.Shogun.Util import *
#from shogun.Features import *
#import pdb
import sys
from sqlite3 import *
import subprocess
from Bio import SeqIO
from Bio.Seq import Seq


## 
taxid_file = "../FINAL_UNIONED_TAXIDS_2_3_4_5.TXT"
taxids = []
f = open(taxid_file, "r")
for l in f:
    taxids.append(int(l))
print len(taxids)
f.close()

## load gi2taxa
print "loading gi2taxa..."
mgtDbDir = "/home/ssul/work/mgtdata"
taxaDir = os.path.join(mgtDbDir,"taxonomy")
giToTaxa = loadGiTaxBin(os.path.join(taxaDir,"gi_taxid.pkl.gz"))
print "# giToTaxa dict =", len(giToTaxa)

##
#print "indexing refseq fasta..."
#all_fasta_file_name = "/export/blastdb_microbial/all_fasta/refseq_microbial.fasta"
##handle = open(all_fasta_file_name, "rU")
#record_dict = SeqIO.index(all_fasta_file_name, "fasta")
##handle.close()
#print record_dict["gi:"]
 
logf = open("not_found_taxids.txt","w")
i = 0
for taxid in taxids:
    print "###", i, taxid
    all_gis_file_name = "/export/blastdb_microbial/all_fasta/refseq_microbial_all_gi.txt"
    gi_file = open(all_gis_file_name, "r")
    num_gis = 0
    of = open("./gis-files/"+str(taxid)+".gis", "w")
    for gi in gi_file:
        gi = int(gi)        
        taxid2 = giToTaxa[gi]        
        if int(taxid) == int(taxid2):
            of.write(str(gi)+"\n")
            num_gis+=1
    of.close()
    print "%d is saved for %d" % (num_gis,int(taxid))
    
    if num_gis == 0:
        logf.write(str(taxid)+"\n")
    #assert(num_gis > 0)
    i+=1
    gi_file.close()

logf.close()
###
#i = 0
#for taxid in taxids:
    #print "###", i, taxid
    #all_fasta_file_name = "/export/blastdb_microbial/all_fasta/refseq_microbial.fasta"
    #seqFile = open(all_fasta_file_name, "r")    
    #num_gis = 0
    #of = open("./gis-files/"+str(taxid)+".gis", "w")
    #handle = open(all_fasta_file_name, "rU")
    #for cur_record in SeqIO.parse(handle, "fasta") :
        #defline = cur_record.id.split()[0][0:]
        #gi = int(defline.split("|")[1])
        
        #taxid2 = giToTaxa[gi]
        
        #if int(taxid) == int(taxid2):
            #of.write(str(gi)+"\n")
            #num_gis+=1
    #handle.close()
    #of.close()
    #print "%d is saved for %d" % (num_gis,int(taxid))
    #i+=1


## EOF
