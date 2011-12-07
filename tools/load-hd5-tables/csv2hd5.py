### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##
#
# See COPYING file distributed along with the MGTAXA package for the
# copyright and license terms.
#
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##

#!/usr/bin/env python

import sys
import os
import struct
from numpy import *
import tables as t
print 'tables.__version__',t.__version__
import string

if __name__ == '__main__':

    if len(sys.argv) != 3:
        print "python csv2hd5.py infile outfile"
        sys.exit(1)    
        
    inFileName = sys.argv[1]    
    outFilename = sys.argv[2]   
    
    ##
    ## Define a user record to characterize some kind of particles
    ##
    
    ## Orig struct in mrblast.cpp 
    """
    structBlResToSaveHitsMason {
        uint32_t    gi;
        uint32_t    readId;
        char        readStrand;     /// 'f' or 'r'
        char        readPairId;     /// '0' or '1'
        uint32_t    subjectId;
        uint32_t    qStart;
        uint32_t    qEnd;
        uint32_t    sStart;
        uint32_t    sEnd;
        double      eValue;
        float       bitScore;
        uint32_t    upperStart;
        uint32_t    upperEnd;
        float       identity;
        float       coverage;    
    }
    """
    #86604733,78067,r,1,87302878,0,400,83655,84053,3.25569269388e-29,136.387893677,1973265,1973665,74.25,100.0
    
    class BlHits(t.IsDescription):  
        gi         = t.UInt32Col()
        readId     = t.UInt32Col()
        readStrand = t.StringCol(1)
        readPairId = t.StringCol(1) 
        sId        = t.UInt32Col()
        qStart     = t.UInt32Col()
        qEnd       = t.UInt32Col()
        sStart     = t.UInt32Col()
        sEnd       = t.UInt32Col()        
        eValue     = t.FloatCol()       ## double-precision ==> 'd'
        bitScore   = t.UInt32Col()      ## single-precision ==> 'f'
        upperStart = t.UInt32Col()
        upperEnd   = t.UInt32Col()
        dIdent     = t.Float32Col()     ## single-precision ==> 'f'
        dCover     = t.Float32Col()     ## single-precision ==> 'f'
        
    ## Open a file in "w"rite mode
    h5file = t.openFile(outFilename, mode = "w", title = "BLAST hits")
    ## Create a new group under "/" (root)
    root = h5file.root
    group = h5file.createGroup(root, "blhits", "blhits")
    ## Create one table on it
    table = h5file.createTable(group, "blhitstab", BlHits, "blhitstab")
    ## Fill the table with 10 particles
    BlHits = table.row
    
    numTotalHits = 0
    with open(inFileName) as infile:
        for line in infile:
            tok = line.strip().split(',')
            assert(len(tok) == 15)
            BlHits['gi']         = int(tok[0])
            BlHits['readId']     = int(tok[1])
            BlHits['readStrand'] = tok[2]
            BlHits['readPairId'] = tok[3]                   
            BlHits['sId']        = int(tok[4])
            BlHits['qStart']     = int(tok[5])
            BlHits['qEnd']       = int(tok[6])
            BlHits['sStart']     = int(tok[7])
            BlHits['sEnd']       = int(tok[8])
            BlHits['eValue']     = float(tok[9])
            BlHits['bitScore']   = float(tok[10])
            BlHits['upperStart'] = int(tok[11])
            BlHits['upperEnd']   = int(tok[12])
            BlHits['dIdent']     = float(tok[13])
            BlHits['dCover']     = float(tok[14])
            BlHits.append()
            numTotalHits += 1
            if numTotalHits % 10000 == 0:
                print numTotalHits
    table.flush()   # flush recordData in the table
    h5file.flush()  # flush all pending recordData
    print "Total number of hits = ", numTotalHits
    
    ## Create index
    #table.cols.gi.createIndex()
    #table.cols.sId.createIndex()
    
    h5file.close()
    infile.close()
    
"""
    ## HDF5 file info
    print h5file
    
    ## Access columns
    table = h5file.root.blhits.blhitstab
    gi = [ x['gi'] for x in table.iterrows() ]
    sid = [ x['sId'] for x in table.iterrows() ]
    #names = [ x['name'] for x in table if x['TDCcount'] > 3 and 20 <= x['pressure'] < 50 ]
    #names = [ x['name'] for x in table.where('(TDCcount > 3) & (20 <= pressure) & (pressure < 50)') ]

    #print gi[0], sid[0]
    # Reading rows
    #for r in table.iterrows():
        #print r['gi'], r['eValue'], r['dIdent'], r['dCover'] 
        #print [r[i] for i in range(0,15)]
        
    h5file.close()
 
"""

"""
    ## Create column object and read
    gcolumns = h5file.createGroup(root, "columns", "")   
    h5file.createArray(gcolumns, 'gi', array(gi), "gi array")
    h5file.createArray(gcolumns, 'sid', array(sid), "sid array")
    #print h5file    
    gi2 = h5file.root.columns.gi.read()
    #print gi2
    
    ## Reading rows
    #for r in table.iterrows():
        #print r['gi']
    
    ## Access recordData
    print table.cols.gi[0]
    
    ## Create index
    table.cols.gi.createIndex()
    table.cols.sId.createIndex()
    print h5file 
    
"""
    
### EOF
