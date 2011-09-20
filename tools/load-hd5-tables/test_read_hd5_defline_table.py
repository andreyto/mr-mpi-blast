#!/usr/bin/env python
#from sqlite3 import *
import sys
import os
import struct
from numpy import *
import tables as t
print 'tables.__version__',t.__version__
#from numexpr import *
#from sqlite3 import *


if __name__ == '__main__':

    if len(sys.argv) != 2:
        print "python make_sqlitedb.py fname"
        sys.exit(1)    
    
    #topDir = sys.argv[1]    
    filename = sys.argv[1]
    #bMakeCSV = int(sys.argv[3])
    #numDir = int(sys.argv[4])
    
    #dirName = [ "div"+str(x) for x in range(numDir) ]
    ## Get the dir list under topDir
    #dirName = []
    #dirList = os.listdir(topDir)
    ##print dirList
    #for na in dirList:
		#if os.path.isdir(na):
			#dirName.append(na)

    #print dirName, len(dirName)
    #numDir = len(dirName)

    ###
    ### Define a user record to characterize some kind of particles
    ###
    #class DefLines(t.IsDescription):
            #gi         = t.UInt32Col()
            #readId     = t.UInt32Col()
            #strandId   = t.StringCol(1)
            #pairId     = t.StringCol(1) 
            #readLen    = t.UInt32Col()
            #cutStart   = t.UInt32Col()
            #cutEnd     = t.UInt32Col()            
            #readStart  = t.UInt32Col()
            #readEnd    = t.UInt32Col()

    ## Open a file in "w"rite mode
    h5file = t.openFile(filename)
    ## Create a new group under "/" (root)
    root = h5file.root
    #group = h5file.createGroup(root, "deflines", "deflines")
    ## Create one table on it
    #table = h5file.createTable(group, "deflinetab", DefLines, "deflinetab")
    ## Fill the table with 10 particles
    #DefLines = table.row
 
    ## HDF5 file info
    print h5file
    
    ## Access columns
    table = h5file.root.deflines.deflinetab
    #gi = [ x['gi'] for x in table.iterrows() ]
    #sid = [ x['sId'] for x in table.iterrows() ]
    #names = [ x['name'] for x in table if x['TDCcount'] > 3 and 20 <= x['pressure'] < 50 ]
    #names = [ x['name'] for x in table.where('(TDCcount > 3) & (20 <= pressure) & (pressure < 50)') ]
    #print gi
    #gi = [row['gi'] for row in ro.where('gi > 10')]
    cutEnd = [ x['cutEnd'] for x in table.where("""(gi > 0)""")]
    print cutEnd
    
    ## Create column object and read
    #gcolumns = h5file.createGroup(root, "columns", "")   
    #h5file.createArray(gcolumns, 'gi', array(gi), "gi array")
    #h5file.createArray(gcolumns, 'sid', array(sid), "sid array")
    ##print h5file    
    #gi2 = h5file.root.columns.gi.read()
    ##print gi2
    
    ## Reading rows
    #for r in table.iterrows():
        #print r['gi']
    
    ## Access recordData
    #print table.cols.gi[0]
    
    ## Create index
    #table.cols.gi.createIndex()
    #table.cols.sId.createIndex()
    #print h5file 
 
    
### EOF
