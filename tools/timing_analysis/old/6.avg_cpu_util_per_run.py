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

nCores = [64, 128, 256, 512, 1024]
   
if __name__ == '__main__':
     
    if len(sys.argv) != 1:
        print "python cov_*.py"
        sys.exit(1)    
    
    ###
    avgCPUUtil = npy.zeros(len(nCores),dtype=float) 
    
    for nC in range(len(nCores)):
        dbName = str(nCores[nC]) + ".db"
        conn = connect(dbName)
        curs = conn.cursor()
        
        ###
        numSlices = 100 
        curs.execute("select min(wallclock), max(wallclock) from item")
                      
        res = curs.fetchone()
        start = res[0]
        stop =  res[1]
        inc = (stop - start) / numSlices
        print start, stop, stop - start, "step = ", inc
        
        core = nCores[nC]
        print "Num core = ", core
        
        ###
        vecStart = [] 
        vecEnd   = []
        cov_busy = npy.zeros(numSlices,dtype=float)    
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
            
        for i in range(len(vecEnd)):
            cov_busy[vecStart[i]:vecEnd[i]] += 1
        
        cnt = 0
        for elem in reversed(cov_busy):
            #print elem
            if (int(elem) <= core * 0.2):
                cnt += 1
            if (int(elem) > core * 0.2):
                break
                
        print "cnt = ", cnt
        cutoff = stop - (cnt * inc)
        print stop, cutoff
        
        ###
        vecOrigStart = []
        vecOrigEnd   = []
        vecRusageUStart = []
        vecRusageUEnd   = []
        vecCPUUtil = []
        vecDiff_wallclock = [] 
        vecDiff_rusageuser = []
        
        curs.execute("select mpi_wtime, wallclock, rusage_ut, startends \
                      from item where wallclock <= " + str(cutoff) \
                      + "order by rec____ID asc")
                      
        temp1 = ""
        temp2 = ""
        for row in curs:
            ### If you want to exclude query build time
            #if row[3] == 0: ### blast call start time
            ### If you want to include query build time
            if row[3] == 2: ### query build start time
                temp1 = row[1]
                temp2 = row[2]
            elif row[3] == 1: ### End time
                vecOrigStart.append(temp1)
                vecRusageUStart.append(temp2)
                vecOrigEnd.append(row[1])
                vecRusageUEnd.append(row[2])

        print len(vecOrigStart), len(vecOrigEnd), len(vecRusageUStart), len(vecRusageUEnd)
                
        ###
        for i in range(len(vecOrigEnd)):            
            diff_wallclock = 0.0
            diff_rusageuser = 0.0
            diff_wallclock = vecOrigEnd[i] - vecOrigStart[i]
            diff_rusageuser = vecRusageUEnd[i] - vecRusageUStart[i]
            vecDiff_wallclock.append(diff_wallclock)
            vecDiff_rusageuser.append(diff_rusageuser)
                
        ### Avg
        totalDiff_wallclock = 0.0
        totalDiff_rusageuser = 0.0
        for i in range(len(vecOrigEnd)):
            totalDiff_wallclock += vecDiff_wallclock[i]
            totalDiff_rusageuser += vecDiff_rusageuser[i]
        avgCPUUtil[nC] = totalDiff_rusageuser / totalDiff_wallclock
        print avgCPUUtil
        
        ###
        curs.close()
        conn.close()
        
    ####
    #### Plotting
    ####
    pylab.figure(1)
    ax=pylab.subplot(111)
    y = range(0, len(nCores))
    ##print y
    ax.plot(y, avgCPUUtil, 'bo', linewidth=1, markersize=6)
    
    ##ax.set_ylim(0, 800)
    ax.set_xlim(-1, len(nCores))

    #ax.set_xticks(x, minor=True)
    #ax.set_yticks(y, minor=True)

    ax.grid(color='k', linestyle='-.', linewidth=0.5)
    ax.yaxis.grid(True, linestyle='-.', which='minor')
    ax.xaxis.grid(True, linestyle='-.', which='minor')

    ax.set_ylabel('Avg. CPU utilization', fontsize=20)
    ax.set_xlabel('Num cores', fontsize=20)
    
    ax.set_xticklabels(('', '64', '128','256','512','1024', ''))
    
    fontsize=16
    for tick in ax.xaxis.get_major_ticks():
        tick.label1.set_fontsize(fontsize)
    for tick in ax.yaxis.get_major_ticks():
        tick.label1.set_fontsize(fontsize)
        
    imFileName = "avg_cpu_util_over_run.png"
    pylab.savefig(imFileName, dpi=(300))
    pylab.show()
    
### EOF
