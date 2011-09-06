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
        print "python make_sqlitedb.py topDir "
        sys.exit(1)    
    
    topDir = sys.argv[1]    
    #filename = sys.argv[2]
    #bMakeCSV = int(sys.argv[3])
    #numDir = int(sys.argv[4])
    
    ## Get the dir list under topDir
    dirName = []
    dirList = os.listdir(topDir)
    #print dirList
    for na in dirList:
		if os.path.isdir(na):
			dirName.append(na)

    print dirName, len(dirName)
    numDir = len(dirName)
	       
	## Read fasta and make sqlite db
    for j in range(numDir):
                   
        fastaFileName = ""
        subDir = topDir+"/"+dirName[j]
        for f in os.listdir(subDir):
            if f.find(".fasta") > -1: 
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
            #strandId   = t.UInt32Col()
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
                t = l.split("_")
                #print t
                
                ## Simulated read case #########################################
                #['>gi|9999|', '99', '1', 'f', '100', '982', '1082', '982', '1082\n']
                ## Illumina 100bp
                #>gi|262276709|ref|NZ_GG704932.1|_0_0_f_100_248911_249010_248911_249011
                ## 454 400bp
                #>gi|83952607|ref|NZ_AALY01000004.1|_0_0_f_398_94061_94461_94061_94459
                
                gi = t[0].split("|")[1]
                #strand = 0
                #if t[-6] == "r":
                    #strand = 1
                
                ## 262276709|82637|0|0|100|559841|559940|559841|559941
                #cmd = "insert into item values (" \
                    #+ str(gi) + "," \
                    #+ str(t[-8]) + "," \
                    #+ str(t[-7]) + "," \
                    #+ str(strand) + "," \
                    #+ str(t[-5]) + "," \
                    #+ str(t[-4]) + "," \
                    #+ str(t[-3]) + "," \
                    #+ str(t[-2]) + "," \
                    #+ str(t[-1]) + ")"                     
                #curs.execute(cmd)
                
                ##
                ##
                DefLines['gi']        = int(gi)
                DefLines['readId']    = int(t[-8])
                DefLines['pairId']    = int(t[-7])
                #DefLines['strandId']  = int(strand)
                DefLines['strandId']  = str(t[-6])
                DefLines['readLen']   = int(t[-5])                    
                DefLines['cutStart']  = int(t[-4])
                DefLines['cutEnd']    = int(t[-3])
                DefLines['readStart'] = int(t[-2])
                DefLines['readEnd']   = int(t[-1])
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
        
        print h5file
        h5file.close()
## EOF
