#!/usr/bin/env python
from sqlite3 import *
import sys
import os
import struct

    
if __name__ == '__main__':

    if len(sys.argv) != 4:
        print "python make_sqlitedb.py hit_file_prefix(+hits-*) database__name 0/1_for_saving_csv"
        sys.exit(1)    
    
    hitFilePrefix = sys.argv[1]    
    dbName = sys.argv[2]
    makeCSV = int(sys.argv[3])
    
    ### load hit file names from hitfilelist
    vecHitFileName = []
    dirname = "./"
    for f in os.listdir(dirname):
        if f.find(hitFilePrefix+"hits-") > -1: 
            vecHitFileName.append(f)
    numHitFiles = len(vecHitFileName)
    print "num of hit files = ", numHitFiles
    
    ###
    ### Create database and table
    ###
    print "DB name = ", dbName
    conn = connect(dbName)
    curs = conn.cursor()
    
    ### Drop table item and itemvendor
    curs.execute('''DROP TABLE if exists item''')
    conn.commit()

    ### Create table
    #curs.execute('''CREATE TABLE IF NOT EXISTS item
      #( recId       integer primary key, 
        #cid         integer,
        #gi          integer,
        #sId         integer,
        #pIdent      double,
        #alignLen    integer,
        #misMatches  integer, 
        #gapOpens    integer,
        #qStart      integer,
        #qEnd        integer,
        #sStart      integer, 
        #sEnd        integer,
        #eValue      double,
        #bitScore    integer,
        #cutStart    integer,
        #cutEnd      integer,
        #dIdent      double, 
        #dCover      double
       #)''')
    curs.execute('''CREATE TABLE IF NOT EXISTS item
      ( recId       integer primary key, 
        gi          integer,
        sId         integer,
        qStart      integer,
        qEnd        integer,
        sStart      integer, 
        sEnd        integer,
        eValue      double,
        bitScore    integer,
        upperStart  integer,
        upperEnd    integer,
        dIdent      double, 
        dCover      double
       )''')

    ###
    ### Load data and insert into db table
    ###
    recNum = 0;
    if makeCSV:
        csvFileName = dbName + ".csv"
        csvFile = open(csvFileName, 'w')
    
    size = struct.calcsize('IIIIIIdfIIffI')

    for i in range(numHitFiles):
        hitFile = open(vecHitFileName[i], "rb")
        data = hitFile.read(size)
        while True:
            try:
                s = struct.unpack('IIIIIIdfIIffI', data)
                cmd = "insert into item values (" \
                    + str(recNum) + "," \
                    + str(s[0]) + "," \
                    + str(s[1]) + "," \
                    + str(s[2]) + "," \
                    + str(s[3]) + "," \
                    + str(s[4]) + "," \
                    + str(s[5]) + "," \
                    + str(s[6]) + "," \
                    + str(s[7]) + "," \
                    + str(s[8]) + "," \
                    + str(s[9]) + "," \
                    + str(s[10]) + "," \
                    + str(s[11]) \
                    + ")" 
                curs.execute(cmd)
                
                if makeCSV:
                    csvString = str(recNum) + "," \
                        + str(s[0]) + "," \
                        + str(s[1]) + "," \
                        + str(s[2]) + "," \
                        + str(s[3]) + "," \
                        + str(s[4]) + "," \
                        + str(s[5]) + "," \
                        + str(s[6]) + "," \
                        + str(s[7]) + "," \
                        + str(s[8]) + "," \
                        + str(s[9]) + "," \
                        + str(s[10]) + "," \
                        + str(s[11]) + "\n"
                    csvFile.write(csvString)
                
                recNum += 1
                data = hitFile.read(size)
            except:
                break
            
        conn.commit()            
        hitFile.close()
        
    if makeCSV:
        csvFile.close()    
    
    ### 
    curs.execute("select count(*) from item")
    print "num records = ", curs.fetchone()[0]

    
    curs.close()
    
### EOF
