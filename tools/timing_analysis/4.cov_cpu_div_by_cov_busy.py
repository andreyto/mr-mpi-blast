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
    print start, stop, stop - start, "step = ", inc
    
    ###
    curs.execute("select mpi_wtime, wallclock, startends from item order by rec____ID asc")
 
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
        
    #############################################
    
    cov_cpu = npy.zeros(numSlices,dtype=float)

    #### 
    curs.execute("select min(wallclock), max(wallclock), min(rusage_ut), max(rusage_ut) from item")
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
    curs.execute("select mpi_wtime, wallclock, rusage_ut, startends from item order by rec____ID asc")

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
    for i in range(len(vecStart)):
        diff_wallclock = vecOrigEnd[i] - vecOrigStart[i]
        diff_rusageuser = vecRusageUEnd[i] - vecRusageUStart[i]
        cpu_util = diff_rusageuser / diff_wallclock
        cov_cpu[vecStart[i]:vecEnd[i]] += cpu_util
        
    cov_cpu /= (cov_busy+0.01)    

    curs.close()
    conn.close()
 
    ####
    #### Plotting
    ####
    pylab.figure(1)
    ax=pylab.subplot(111)
    y = range(numSlices)
    ##print y
    ax.plot(y, cov_cpu, 'bo-', linewidth=1, markersize=6)
    
    
    #res_from_3 = [  9.98315110e-01,   9.98315110e-01,   9.98315110e-01,   9.98315110e-01,
   #9.98315110e-01,   9.98315110e-01,   9.98315110e-01,   9.98315110e-01,
   #9.98315110e-01,   9.98315110e-01,   9.98315110e-01,   9.98315110e-01,
   #9.98315110e-01,   9.98315110e-01,   9.98315110e-01,   9.98315110e-01,
   #9.98315110e-01,   9.98322399e-01,   9.98327431e-01,   9.98329672e-01,
   #9.98345774e-01,   9.98396182e-01,   9.98435782e-01,   9.98459885e-01,
   #9.98538857e-01,   9.98504077e-01,   9.98508818e-01,   9.98508818e-01,
   #9.98508818e-01,   9.98508818e-01,   9.98508818e-01,   9.98508818e-01,
   #9.98508818e-01,   9.98508818e-01,   9.98508818e-01,   9.98508818e-01,
   #9.98508818e-01,   9.98479780e-01,   9.98479780e-01,   9.98451447e-01,
   #9.98471687e-01,   9.98390510e-01,   9.98390510e-01,   9.98390947e-01,
   #9.98401464e-01,   9.98392907e-01,   9.98336014e-01,   9.98333846e-01,
   #9.98332603e-01,   9.98326430e-01,   9.98279970e-01,   9.98279576e-01,
   #9.98279576e-01,   9.98279576e-01,   9.98279576e-01,   9.98279576e-01,
   #9.98295078e-01,   9.98240974e-01,   9.97284057e-01,   9.98245520e-01,
   #9.98278249e-01,   9.98239237e-01,   9.98269611e-01,   9.96235377e-01,
   #9.98108207e-01,   9.98107583e-01,   9.98126227e-01,   9.98101674e-01,
   #9.98135195e-01,   9.98142630e-01,   9.98131385e-01,   9.98131385e-01,
   #9.98131385e-01,   9.98131385e-01,   9.98131385e-01,   9.98154813e-01,
   #9.98181566e-01,   9.98185777e-01,   9.98179912e-01,   9.03537185e-01,
   #8.54748590e-01,   5.93338391e-01,   4.29412836e-01,   3.92322916e-01,
   #3.72808505e-01,   3.26946936e-01,   3.17191272e-01,   2.92788121e-01,
   #2.62541437e-01,   2.19582773e-01,   2.04950880e-01,   1.79572444e-01,
   #1.78596683e-01,   1.78596683e-01,   1.78596683e-01,   1.78596683e-01,
   #1.71773053e-01,   1.21030045e-01,   3.51450614e-02,   9.76368380e-04]

    #ax.plot(y, res_from_3, 'rd-', linewidth=1, markersize=6)
    
    ##ax.set_ylim(0, 800)
    ##ax.set_xlim(16, 2048)

    #ax.set_xticks(x, minor=True)
    #ax.set_yticks(y, minor=True)

    ax.grid(color='k', linestyle='-.', linewidth=0.5)
    ax.yaxis.grid(True, linestyle='-.', which='minor')
    ax.xaxis.grid(True, linestyle='-.', which='minor')

    ax.set_ylabel('"Useful" CPU utilization over time', fontsize=20)
    ax.set_xlabel('Time from start (% of total run-time)', fontsize=20)
        
    ax.set_ylim(0, 1.0)
    #ax.set_xlim(32, 1024)

    fontsize=16
    for tick in ax.xaxis.get_major_ticks():
        tick.label1.set_fontsize(fontsize)
    for tick in ax.yaxis.get_major_ticks():
        tick.label1.set_fontsize(fontsize)
    
    #ax.legend( (
            #'average CPU utilization',
            #'"Userful" CPU utilization'
            #), loc=3)
            
    #leg = plt.gca().get_legend()
    #ltext  = leg.get_texts()
    #plt.setp(ltext, fontsize='14')    # the legend text fontsize

    imFileName = dbName + "-cov_cpu_div_by_cov_busy.png"
    pylab.savefig(imFileName, dpi=(300))
    pylab.show()
    
### EOF
