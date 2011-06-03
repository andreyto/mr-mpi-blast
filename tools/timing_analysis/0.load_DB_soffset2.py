#!/usr/bin/env python
from sqlite3 import *
import sys
    
if __name__ == '__main__':

    if len(sys.argv) != 4:
        print "python analysis.py dbName logfilenameprefix nproc"
        sys.exit(1)    
    
    dbName = sys.argv[1]
    logfinenameprefix = sys.argv[2]    
    nProc = sys.argv[3]
    
    ###
    ### Create database and table
    print "DB name = ", dbName
    conn = connect(dbName)
    curs = conn.cursor()

    # Drop table item and itemvendor
    curs.execute('''DROP TABLE if exists item''')
    conn.commit()

    ### Create table
    curs.execute('''CREATE TABLE IF NOT EXISTS item
      ( rec____ID integer primary key, 
        rank___ID integer,
        startends integer, 
        mpi_wtime double,
        wallclock double,
        rusage_ut double, 
        rusage_st double,
        db___name char(256),
        call___ID integer, 
        proc_name char(256),
        s__offset double
       )''')
    
    recNum = 1;
    csvFileName = dbName + ".csv"
    csvFile = open(csvFileName, 'w')
    csvFile.write("#################################################################################################################################\n")
    csvFile.write("# UID, rank, OP_CODE, MPI_Wtime (sec), wallclock, rusage user time, rusage system time, db name, call id, proc name, start offset\n")
    csvFile.write("# OP_CODE = 2:query loading start, 3:db loading start, 0:blast call start, 1:blast call end\n")
    csvFile.write("#################################################################################################################################\n")
    for np in range(1, int(nProc)):
        logfilename = logfinenameprefix + str(np) + "-log.txt"
        print logfilename
        logfile = open(logfilename, "r")
        for line in logfile:
            if line.find("Rank:") > -1:
                if line.find("blast") > -1:
                    s = line.split(",") 
                    if (len(s) > 2):                        
                        cmd = "insert into item values (" \
                            + str(recNum) + "," \
                            + str(np) + "," 
                        csvString = str(recNum) + "," + str(np) + "," 
                        if s[0].find("starts") > -1: 
                            cmd += "0,"
                            csvString += "0,"
                        else:
                            cmd += "1,"
                            csvString += "1,"
                            
                        cmd += str(s[1]) + "," \
                            + str(s[2]) + "," \
                            + str(s[3]) + "," \
                            + str(s[4]) + "," \
                            + "'" + str(s[5]) + "'," \
                            + str(s[6]) + "," \
                            + "'" + str(s[7]) + "'," \
                            + str(s[8]) \
                            + ")" 
                            
                        csvString += str(s[1]) + "," \
                            + str(s[2]) + "," \
                            + str(s[3]) + "," \
                            + str(s[4]) + "," \
                            + str(s[5]).strip() + "," \
                            + str(s[6]) + "," \
                            + str(s[7]).strip() + "," \
                            + str(s[8])

                        curs.execute(cmd)
                        recNum += 1
                        csvFile.write(csvString)
                elif line.find("query") > -1:
                    s = line.split(",") 
                    if (len(s) > 2):
                        cmd = "insert into item values (" \
                            + str(recNum) + "," \
                            + str(np) + "," 
                        csvString = str(recNum) + "," + str(np) + "," 
                        if s[0].find("starts") > -1: 
                            cmd += "2,"
                            csvString += "2,"
                        #else:
                            #cmd += "3,"
                            #csvString += "3,"
                            
                        cmd += str(s[1]) + "," \
                            + str(s[2]) + "," \
                            + str(s[3]) + "," \
                            + str(s[4]) + "," \
                            + "'" + str(s[5]) + "'," \
                            + str(s[6]) + "," \
                            + "'" + str(s[7]) + "'," \
                            + str(s[8]) \
                            + ")" 
                            
                        csvString += str(s[1]) + "," \
                            + str(s[2]) + "," \
                            + str(s[3]) + "," \
                            + str(s[4]) + "," \
                            + str(s[5]).strip() + "," \
                            + str(s[6]) + "," \
                            + str(s[7]).strip() + "," \
                            + str(s[8])
                            
                        curs.execute(cmd)
                        recNum += 1
                        csvFile.write(csvString)                   
                elif line.find("db_loading") > -1:
                    s = line.split(",") 
                    if (len(s) > 2):                        
                        cmd = "insert into item values (" \
                            + str(recNum) + "," \
                            + str(np) + "," 
                        csvString = str(recNum) + "," + str(np) + "," 
                        if s[0].find("starts") > -1: 
                            cmd += "3,"
                            csvString += "3,"
                        #else:
                            #cmd += "1,"
                            #csvString += "1,"
                            
                        cmd += str(s[1]) + "," \
                            + str(s[2]) + "," \
                            + str(s[3]) + "," \
                            + str(s[4]) + "," \
                            + "'" + str(s[5]) + "'," \
                            + str(s[6]) + "," \
                            + "'" + str(s[7]) + "'," \
                            + str(s[8]) \
                            + ")" 
                            
                        csvString += str(s[1]) + "," \
                            + str(s[2]) + "," \
                            + str(s[3]) + "," \
                            + str(s[4]) + "," \
                            + str(s[5]).strip() + "," \
                            + str(s[6]) + "," \
                            + str(s[7]).strip() + "," \
                            + str(s[8])

                        curs.execute(cmd)
                        recNum += 1
                        csvFile.write(csvString)
        conn.commit()
    
    csvFile.close()
    
    ### 
    curs.execute("select count(*) from item")
    print curs.fetchone()[0]

    
    curs.close()
    
### EOF

