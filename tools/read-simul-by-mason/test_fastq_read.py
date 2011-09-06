#!/usr/bin/env python

import HTSeq
import itertools
import sys

if __name__ == '__main__':
     
    if len(sys.argv) != 2:
        print "python xxxx.py fileName"
        sys.exit(1)  
    fileName = sys.argv[1]
    
    ##
    ##
    ##
    #fastq_file = HTSeq.FastqReader( fileName )

    #vecLen = []
    #for read in itertools.islice( fastq_file, 1000000000 ):
        ##print read, len(read)
        #vecLen.append(len(read))
        
    #print "num reads=", len(vecLen)
    #print "min=", min(vecLen)
    #print "max=", max(vecLen)

    
    ##
    ## 
    ##
    for s in HTSeq.FastqReader( fileName ):
        print "Sequence '%s' has length %d." % ( s.name, len(s) )
        newDefline = s.name+"
    


