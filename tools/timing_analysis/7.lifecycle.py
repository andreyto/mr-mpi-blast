#!/usr/bin/env python
from sqlite3 import *

import sys
import numpy as npy
import matplotlib.pyplot as plt
import pylab

### table spec. ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
    #curs.execute('''CREATE TABLE IF NOT EXISTS item
      #( rec____ID integer primary key, 
        #rank___ID integer,
        #startends integer, 
        #mpi_wtime double,
        #wallclock double,
        #rusage_ut double, 
        #rusage_st double,
        #db___name char(256),
        #call___ID integer, 
        #proc_name char(256),
        #s__offset double
       #)''')
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###    
   
if __name__ == '__main__':
     
    if len(sys.argv) != 3:
        print "python lifecycle.py dbName nodenum"
        sys.exit(1)    
    
    dbName = sys.argv[1]
    nd = int(sys.argv[2])
    conn = connect(dbName)
    curs = conn.cursor()
    
    #### 
    curs.execute("select min(wallclock), max(wallclock), min(rusage_ut), max(rusage_ut) from item")
    res = curs.fetchone()
    start = res[0]
    stop =  res[1]
    print start, stop, stop - start
    numSlices = 1000
    inc = (stop - start) / numSlices
    life = npy.zeros(numSlices, dtype=float)
    
    ### Collect node names
    nodenames = []
    curs.execute("select distinct proc_name from item")
    for row in curs:
        nodenames.append(row[0])
        
    print "num node = " , len(nodenames)
    
    #nd = 0
    


    
    #for nd in range(0, len(nodenames)):
    #for nd in range(0, dbfileindex):
    print "nd = " , nd, nodenames[nd]
    
    ####
    dbnames = [] 
    sql = "select distinct db___name from item"
    curs.execute(sql)
    for row in curs:
        dbnames.append(row[0])
    print len(dbnames)
    
    ###
    #dbnames = [] 
    #sql = "select distinct db___name from item where proc_name == '%s' order by wallclock" % (nodenames[nd])
    #curs.execute(sql)
    #for row in curs:
        ##print row[0]          
        #dbnames.append(row[0])
    
    
    ###
    vecStart = [] 
    vecEnd   = []
    vecDbNameIdx = []
    sql = "select startends, wallclock, db___name from item where proc_name == '%s' order by wallclock" % (nodenames[nd])
    curs.execute(sql)
    for row in curs:
        if row[0] == 2: ### query build start time
            vecStart.append(int((row[1] - start) / inc))
        elif row[0] == 1: ### End time
            vecEnd.append(int((row[1] - start) / inc))
            vecDbNameIdx.append(dbnames.index(row[2]) + 1)
    
    print "num blast call = ", len(vecDbNameIdx)
            
         
    curs.close()
    conn.close()
 
    ####
    #### Plotting
    ####
    pylab.figure(1)
    ax = pylab.subplot(111)
    x = range(numSlices)
    ##print x
    #ax.plot(x[vecStart[0]:vecEnd[0]], vecDbNameIdx[0], 'bo', linewidth=1, markersize=6)
    c = ['bx', 'rx', 'gx']
    
    #for nd in range(0, 1):
    for i in range(0, len(vecDbNameIdx)):
        for x in range(vecStart[i], vecEnd[i]+1):
            ax.plot(x, vecDbNameIdx[i], c[1], linewidth=1, markersize=1)
    
    ax.set_ylim(0, 109)
    ax.set_xlim(0, numSlices)

    #ax.set_xticks(x, minor=True)
    #ax.set_yticks(y, minor=True)

    ax.grid(color='k', linestyle='-.', linewidth=0.5)
    ax.yaxis.grid(True, linestyle='-.', which='minor')
    ax.xaxis.grid(True, linestyle='-.', which='minor')

    ax.set_ylabel('DB file (1~109)', fontsize=20)
    ax.set_xlabel('t (%)', fontsize=20)
    
    ax.set_xticklabels(('0', '', '', '', '', '100'))
        
    fontsize=16
    for tick in ax.xaxis.get_major_ticks():
        tick.label1.set_fontsize(fontsize)
    for tick in ax.yaxis.get_major_ticks():
        tick.label1.set_fontsize(fontsize)
        
    #imFileName = dbName + "-life.png"
    #pylab.savefig(imFileName, dpi=(300))
    pylab.show()
    
### EOF




##!/usr/bin/env python
#from sqlite3 import *

#import sys
#import numpy as npy
#import matplotlib.pyplot as plt
#import pylab

#### table spec. ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
    ##curs.execute('''CREATE TABLE IF NOT EXISTS item
      ##( rec____ID integer primary key, 
        ##rank___ID integer,
        ##startends integer, 
        ##mpi_wtime double,
        ##wallclock double,
        ##rusage_ut double, 
        ##rusage_st double,
        ##db___name char(256),
        ##call___ID integer, 
        ##proc_name char(256),
        ##s__offset double
       ##)''')
#### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###    
   
#if __name__ == '__main__':
     
    #if len(sys.argv) != 2:
        #print "python lifecycle.py dbName"
        #sys.exit(1)    
    
    #dbName = sys.argv[1]
    #conn = connect(dbName)
    #curs = conn.cursor()
    
    ##### 
    #curs.execute("select min(wallclock), max(wallclock), min(rusage_ut), max(rusage_ut) from item")
    #res = curs.fetchone()
    #start = res[0]
    #stop =  res[1]
    #print start, stop, stop - start
    #numSlices = 1000
    #inc = (stop - start) / numSlices
    ##life = npy.zeros(numSlices, dtype=float)
    
    #### Collect node names
    #nodenames = []
    #curs.execute("select distinct proc_name from item")
    #for row in curs:
        #nodenames.append(row[0])
    #print "num node = " , len(nodenames)
    
    
    ####
    #dbnames = [] 
    #sql = "select distinct db___name from item"
    #curs.execute(sql)
    #for row in curs:
        #dbnames.append(row[0])
    #print len(dbnames)
        
        
    ##rec = [[] * 3]
    #nd = 0
    ##for nd in range(0, len(nodenames)):
    ##for nd in range(0, 3):
    #print "nd = " , nd, nodenames[nd]
    
    ####
    ##vecStart = [[]] 
    ##vecEnd   = [[]]
    ##vecDbNameIdx = [[]]
    #rec = []
    #sql = "select startends, wallclock, db___name from item where proc_name == '%s' order by wallclock" % (nodenames[nd])
    #curs.execute(sql)
    #for row in curs:
        #if row[0] == 2: ### query build start time
            ##vecStart[nd].append(int((row[1] - start) / inc))
            #temp1  = int((row[1] - start) / inc)
        #elif row[0] == 1: ### End time
            ##vecEnd[nd].append(int((row[1] - start) / inc))
            #temp2 = int((row[1] - start) / inc)
            ##vecDbNameIdx[nd].append(dbnames.index(row[2]) + 1)
            #temp3 = dbnames.index(row[2]) + 1
            #r = temp1, temp2, temp3
            ##rec.append([] * 3)
            #rec.append(r)
    
    #print "num blast call = ", len(rec)
    ##print rec
         
    #curs.close()
    #conn.close()
 
    ##sys.exit()
    #####
    ##### Plotting
    #####
    #pylab.figure(1)
    #ax = pylab.subplot(111)
    #x = range(numSlices)
    ###print x
    ##ax.plot(x[vecStart[0]:vecEnd[0]], vecDbNameIdx[0], 'bo', linewidth=1, markersize=6)
    #c = ['bx', 'rx', 'gx']
    
    ##for nd in range(0, 1):
    ##nd = 0
    #for i in range(0, len(rec)):
        #for x in range(rec[i][0], rec[i][1]):
            #ax.plot(x, rec[i][2], c[0], linewidth=1, markersize=1)
    
    #ax.set_ylim(0, 109)
    #ax.set_xlim(0, numSlices)

    ##ax.set_xticks(x, minor=True)
    ##ax.set_yticks(y, minor=True)

    #ax.grid(color='k', linestyle='-.', linewidth=0.5)
    #ax.yaxis.grid(True, linestyle='-.', which='minor')
    #ax.xaxis.grid(True, linestyle='-.', which='minor')

    #ax.set_ylabel('DB file (1~109)', fontsize=20)
    #ax.set_xlabel('t (%)', fontsize=20)
    
    #ax.set_xticklabels(('0', '', '', '', '', '100'))
        
    #fontsize=16
    #for tick in ax.xaxis.get_major_ticks():
        #tick.label1.set_fontsize(fontsize)
    #for tick in ax.yaxis.get_major_ticks():
        #tick.label1.set_fontsize(fontsize)
        
    ##imFileName = dbName + "-life.png"
    ##pylab.savefig(imFileName, dpi=(300))
    #pylab.show()
    
#### EOF
