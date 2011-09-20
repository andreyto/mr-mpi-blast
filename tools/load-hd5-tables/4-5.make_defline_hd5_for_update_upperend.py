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
        print "python 4-5.make_defline_hd5_for_update_upperend.py topDir "
        sys.exit(1)    
    topDir = sys.argv[1]    
    
    ## Get the dir list under topDir
    dirName = []
    dirList = os.listdir(topDir)
    #print dirList
    for na in dirList:
		if os.path.isdir(na): ###################
			dirName.append(na)

    print dirName, len(dirName)
    numDir = len(dirName)
    
	## Read fasta and make sqlite db
    for j in range(numDir):
                   
        fastaFileName = ""
        subDir = topDir+"/"+dirName[j]
        for f in os.listdir(subDir):
            if f.find(".fasta") > -1 and f.find(".idx") <= -1: 
                fastaFileName = f
        print fastaFileName 
        
        ## Create database and table
        dbName = topDir+"/"+dirName[j]+"/deflines.hd5"
        print "DB name = ", dbName
        #conn = connect(dbName)
        #curs = conn.cursor()
        #curs.execute('''DROP TABLE if exists item''')
        #conn.commit()
        #curs.execute('''CREATE TABLE IF NOT EXISTS item
          #( gi          integer,
            #readId      integer, 
            #pairId      integer, 
            #strandId    integer,
            #readLen     integer,
            #cutStart    integer, 
            #cutEnd      integer,
            #readStart   integer,
            #readEnd     integer
           #)''')
           
        ###
        ### Define a user record for HD5 table
        ###
        class DefLines(t.IsDescription):
            gi         = t.UInt32Col()
            readId     = t.UInt32Col()
            pairId     = t.UInt32Col()
            strandId   = t.StringCol(1)            
            readLen    = t.UInt32Col()
            cutStart   = t.UInt32Col()
            cutEnd     = t.UInt32Col()            
            readStart  = t.UInt32Col()
            readEnd    = t.UInt32Col()
            
        ## Open a file in "w"rite mode
        h5file = t.openFile(dbName, mode = "w", title = "deflines")
        ## Create a new group under "/" (root)
        root = h5file.root
        group = h5file.createGroup(root, "deflines", "deflines")
        ## Create one table on it
        table = h5file.createTable(group, "deflinetab", DefLines, "deflinetab")
        ## Fill the table with 10 particles
        DefLines = table.row
        
        ## open fasta and load table
        recNum = 0
        fastaFile = open(topDir+"/"+dirName[j]+"/"+fastaFileName, "r")
        for l in fastaFile:
            if l.find(">") > -1: 
                tok = l.split("_")
                #print tok
                
                ## Simulated read cases #########################################
                #['>gi|9999|', '99', '1', 'f', '100', '982', '1082', '982', '1082\n']
                ## Illumina 100bp
                #>gi|262276709|ref|NZ_GG704932.1|_0_0_f_100_248911_249010_248911_249011
                ## 454 400bp
                #>gi|83952607|ref|NZ_AALY01000004.1|_0_0_f_398_94061_94461_94061_94459
                
                gi = tok[0].split("|")[1]
                
                ## 262276709|82637|0|0|100|559841|559940|559841|559941
                #cmd = "insert into item values (" \
                    #+ str(gi) + "," \
                    #+ str(tok[-8]) + "," \
                    #+ str(tok[-7]) + "," \
                    #+ str(strand) + "," \
                    #+ str(tok[-5]) + "," \
                    #+ str(tok[-4]) + "," \
                    #+ str(tok[-3]) + "," \
                    #+ str(tok[-2]) + "," \
                    #+ str(tok[-1]) + ")"                     
                #curs.execute(cmd)
                
                ##
                ##
                DefLines['gi']        = int(gi)
                DefLines['readId']    = int(tok[-8])
                DefLines['pairId']    = int(tok[-7])
                DefLines['strandId']  = str(tok[-6])
                DefLines['readLen']   = int(tok[-5])                    
                DefLines['cutStart']  = int(tok[-4])
                DefLines['cutEnd']    = int(tok[-3])
                DefLines['readStart'] = int(tok[-2])
                DefLines['readEnd']   = int(tok[-1])
                DefLines.append()
                
                recNum += 1
        
        table.flush()   # flush recordData in the table
        h5file.flush()  # flush all pending recordData
        fastaFile.close()
        #conn.commit()    
        recNum += 1
        
        print recNum
        ### 
        #curs.execute("select count(*) from item")
        #print "num records = ", curs.fetchone()[0]

        #curs.close()
    
        ## Create index
        table.cols.gi.createIndex()
        table.cols.readId.createIndex()
        table.cols.pairId.createIndex()
        
        #print h5file
        h5file.close()


""" Read test

    ## Open a file in "w"rite mode
    h5file = t.openFile(filename)
    ## Create a new group under "/" (root)
    root = h5file.root
 
    ## HDF5 file info
    print h5file
    
    ## Access columns
    table = h5file.root.deflines.deflinetab
    mycond = "(gi > 0)"
    cutEnd = [ x['cutEnd'] for x in table.where(mycond)]
    print cutEnd
    h5file.close()

"""

## EOF
