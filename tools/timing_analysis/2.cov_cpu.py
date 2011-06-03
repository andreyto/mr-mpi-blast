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
        print "python analysis.py dbName"
        sys.exit(1)    
    
    dbName = sys.argv[1]
    conn = connect(dbName)
    curs = conn.cursor()
 
    numSlices = 100
    cov_cpu = npy.zeros(numSlices,dtype=float)
 
    #### 
    curs.execute("select min(wallclock), max(wallclock), min(rusage_ut), \
                  max(rusage_ut) from item")
    res = curs.fetchone()
    start = res[0]
    stop =  res[1]
    rusage_ut_start = res[2]
    rusage_ut_end = res[3]
    
    inc = (stop - start) / numSlices
    print start, stop, stop - start, "step = ", inc
    
    vecStart = [] 
    vecEnd   = []
    vecOrigStart = []
    vecOrigEnd   = []
    vecRusageUStart = []
    vecRusageUEnd   = []
    
    ###
    curs.execute("select mpi_wtime, wallclock, rusage_ut, startends \
                  from item order by rec____ID asc")

    for row in curs:
        ### If you want to exclude query build time
        #if row[3] == 0: ### blast call start time
        ### If you want to include query build time
        if row[3] == 2: ### query build start time
            vecStart.append(int((row[1] - start) / inc))
            vecOrigStart.append(row[1])
            vecRusageUStart.append(row[2])
        elif row[3] == 1: ### End time
            vecEnd.append(int((row[1] - start) / inc))
            vecOrigEnd.append(row[1])
            vecRusageUEnd.append(row[2])

    print len(vecStart), len(vecEnd)
    
    ###
    ### increment by (end_blast_cpu - start_blast_cpu)
    ### /(end_blast_wall - start_blast_wall)
    ### 
    for i in range(len(vecEnd)):
        diff_wallclock = vecOrigEnd[i] - vecOrigStart[i]
        diff_rusageuser = vecRusageUEnd[i] - vecRusageUStart[i]
        cpu_util = diff_rusageuser / diff_wallclock
        cov_cpu[vecStart[i]:vecEnd[i]] += cpu_util
    
    curs.close()
    conn.close()
    
    ####
    #### Plotting
    ####
    pylab.figure(1)
    ax=pylab.subplot(111)
    y = range(numSlices)
    ##print y
    ax.plot(y, cov_cpu, 'bo', linewidth=1, markersize=6)
    
    ##ax.set_ylim(0, 800)
    ##ax.set_xlim(16, 2048)

    #ax.set_xticks(x, minor=True)
    #ax.set_yticks(y, minor=True)

    ax.grid(color='k', linestyle='-.', linewidth=0.5)
    ax.yaxis.grid(True, linestyle='-.', which='minor')
    ax.xaxis.grid(True, linestyle='-.', which='minor')

    ax.set_ylabel('CPU utilization over time', fontsize=20)
    ax.set_xlabel('t (%)', fontsize=20)
        
    fontsize=16
    for tick in ax.xaxis.get_major_ticks():
        tick.label1.set_fontsize(fontsize)
    for tick in ax.yaxis.get_major_ticks():
        tick.label1.set_fontsize(fontsize)
        
    imFileName = dbName + "-cov_cpu.png"
    pylab.savefig(imFileName, dpi=(300))
    pylab.show()
    
### EOF

