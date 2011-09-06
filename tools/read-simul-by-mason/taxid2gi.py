from MGT.Taxa import *
#from MGT.Svm import *
from MGT.FastaIO import *
#from MGT.Shogun.Util import *
#from shogun.Features import *
#import pdb
import sys
from sqlite3 import *

def unisq(seq):
    return dict(zip(seq, [None,]*len(seq))).keys()
#l = [1, 2, 3, 4, 5, 5, 6, 7, 7, 7, 2]
#print unisq(l)

if __name__ == '__main__':
     
    if len(sys.argv) != 2:
        print "python *.py fileName taxid"
        sys.exit(1)  
    
    taxid = int(sys.argv[1])
    
    ## load gi2taxa
    print "loading gi2taxa..."
    mgtDbDir = "/home/ssul/work/mgtdata"
    #mgtDbDir = "/usr/local/depot/projects/MGTAXA/mgtaxa/mgtdata"
    taxaDir = os.path.join(mgtDbDir,"taxonomy")
    giToTaxa = loadGiTaxBin(os.path.join(taxaDir,"gi_taxid.pkl.gz"))
    print "# giToTaxa dict =", len(giToTaxa)
    
    ##
    ## LOAD ALL GIS FROM REFSEQ
    ##                
    refseqGIs = []
    for line in open("/home/ssul/work/mgtdata/refseq_microbial_all_gi.txt", "r"):
        refseqGIs.append(int(line.strip())) 
    print "# refseg all GIs = ", len(refseqGIs)

    #final_gis = []
    for gi in refseqGIs:
        taxid2 = giToTaxa[int(gi)]
        if taxid == int(taxid2):
            print "taxid, gi = ", taxid, gi
 

## EOF

