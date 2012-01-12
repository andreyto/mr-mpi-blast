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
      

      
if __name__ == '__main__':
     
    if len(sys.argv) != 2:
        print "python *.py src_dir"
        sys.exit(1) 
        
    path=str(sys.argv[1])
    dirList=os.listdir(path)

    
    ## load gi2taxa
    print "loading gi2taxa..."
    mgtDbDir = "/home/ssul/work/mgtdata"
    taxaDir = os.path.join(mgtDbDir,"taxonomy")
    giToTaxa = loadGiTaxBin(os.path.join(taxaDir,"gi_taxid.pkl.gz"))
    print "# giToTaxa dict =", len(giToTaxa)


    i=0
    for fname in dirList:
        print "#######", i, fname
        inFileName = fname
        outFileName = fname + ".fasta"
        
        ##
        ## Final gi2taxa check
        ##
        taxid_list = []
        gis_file = open(inFileName, "r")
        for l in gis_file:
            gi = int(l)
            taxid = giToTaxa[gi]
            taxid_list.append(taxid)
        gis_file.close()
        taxid_list = unique(taxid_list)
        print "num taxids = ", len(taxid_list)
        #assert(len(taxid_list) == 1)
        
        if len(taxid_list) == 1:
            cmd1 = "fastacmd -d /export/blastdb_microbial/refseq_microbial -i " + inFileName + " -t T -o /home/ssul/work/mrblast-build/data/basic_test/microbial_real/collect-taxid/6.collect_gis_from_unioned_taxids/fasta/" + outFileName
            print cmd1
            p = subprocess.Popen(cmd1, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
            stdout_value = p.communicate()[0]
            s = stdout_value.split()
            for i in range(len(s)):
                print str(i) + " " + str(s[i])    
            
            #outFile.write(str(taxid)+","+str(gi)+","+str(seqLen)+"\n")
            #outFile.write(final_seq+"\n")
            #inFile.close()
            #outFile.close()
            i+=1
        else:
            print "###################### 0 #########################"
            print fname


## EOF
