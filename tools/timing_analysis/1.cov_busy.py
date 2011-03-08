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
     
    if len(sys.argv) != 2:
        print "python cov_*.py dbName"
        sys.exit(1)    
    
    dbName = sys.argv[1]
    conn = connect(dbName)
    curs = conn.cursor()
    numSlices = 100
    cov_busy = npy.zeros(numSlices,dtype=float)    
 
    vecStart = [] 
    vecEnd   = []
 
    #### 
    curs.execute("select min(wallclock), max(wallclock) from item")
    res = curs.fetchone()
    start = res[0]
    stop =  res[1]
    inc = (stop - start) / numSlices
    print start, stop, stop - start, "step = ",inc
    
    ###
    curs.execute("select mpi_wtime, wallclock, startends from item \
                  order by rec____ID asc")
 
    for row in curs:
        ### If you want to exclude query build time
        #if row[2] == 0: ### blast call start time
        ### If you want to include query build time
        if row[2] == 2: ### query build start time
            vecStart.append(int((row[1] - start) / inc))
        elif row[2] == 1: ### blast call end time
            vecEnd.append(int((row[1] - start) / inc))
 
    print len(vecStart), len(vecEnd)
    
    for i in range(len(vecEnd)):
        cov_busy[vecStart[i]:vecEnd[i]] += 1
   
    
    curs.close()
    conn.close()

    ####
    #### Plotting
    ####
    pylab.figure(1)
    ax=pylab.subplot(111)
    y = range(numSlices)

    ax.plot(y, cov_busy, 'ro', linewidth=1, markersize=5)

    ##ax.set_ylim(0, 800)
    ##ax.set_xlim(16, 2048)

    #ax.set_xticks(x, minor=True)
    #ax.set_yticks(y, minor=True)

    ax.grid(color='k', linestyle='-.', linewidth=0.5)
    ax.yaxis.grid(True, linestyle='-.', which='minor')
    ax.xaxis.grid(True, linestyle='-.', which='minor')

    ax.set_ylabel('Number of busy cores over time', fontsize=20)
    ax.set_xlabel('t (%)', fontsize=20)
        
    fontsize=16
    for tick in ax.xaxis.get_major_ticks():
        tick.label1.set_fontsize(fontsize)
    for tick in ax.yaxis.get_major_ticks():
        tick.label1.set_fontsize(fontsize)
        
    imFileName = dbName + "-cov_busy.png"
    pylab.savefig(imFileName, dpi=(300))
    pylab.show()
    
### EOF

