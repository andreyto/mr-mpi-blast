#!/usr/bin/env python
from sqlite3 import *

import sys
import numpy as npy
import matplotlib.pyplot as plt
import pylab

from matplotlib import pyplot, lines 

#curs.execute('''CREATE TABLE IF NOT EXISTS item
  #( recId       integer primary key, 
    #gi          integer,
    #sId         integer,
    ##pIdent      double,
    ##alignLen    integer,
    ##misMatches  integer, 
    ##gapOpens    integer,
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

if __name__ == '__main__':
     
    if len(sys.argv) != 3:
        print "python 2.identity_plot.py dbName subjectIndex"
        sys.exit(1)  
        
    dbName = sys.argv[1]
    subjectIndex = int(sys.argv[2])
    coverCutoff = 0
    identCutoff = 0
    
    ### DB conn
    conn = connect(dbName)
    curs = conn.cursor()

    ### Get unique subID
    vecSubId = []
    sql = "select sId, count(*), gi from item where dCover >= %s group by sId order by count(*) desc"  % (coverCutoff)
    curs.execute(sql)
    for row in curs:
        vecSubId.append(row[0])
    
    print "num distinct subID = ", len(vecSubId)

    currSubId = vecSubId[subjectIndex]
    
    print "currSubId = ", currSubId
    
    ###
    sql = "select min(sStart), max(sEnd) from item where sId == %s and dIdent >= %s and dCover >= %s" % (currSubId, identCutoff, coverCutoff)
    curs.execute(sql)
    res = curs.fetchone()
    minSStart = res[0]
    maxSEnd   = res[1]
    print "min s.start, max s.end = ", minSStart, maxSEnd
    sql = "select count(*) from item where sId == %s and dIdent >= %s and dCover >= %s" % (currSubId, identCutoff, coverCutoff)
    curs.execute(sql)
    res = curs.fetchone()
    numReads = res[0]
    
    ### Get  X_s, Y_s, dIdent
    vecStart  = npy.zeros((numReads),dtype=int) 
    vecEnd    = npy.zeros((numReads),dtype=int) 
    vecDIdent = npy.zeros((numReads),dtype=float) 
    vecGi     = npy.zeros((numReads),dtype=int) 
    
    sql = "select sStart, sEnd, dIdent, gi, dCover from item where sId == %s and dIdent >= %s and dCover >= %s order by gi" % (currSubId, identCutoff, coverCutoff)
    curs.execute(sql)
    i = 0
    #print "s.start, s.end, identity, coverage, gi, sid"
    for row in curs:
        vecStart[i] = (row[0])
        vecEnd[i] = (row[1])
        vecDIdent[i] = (row[2])
        vecGi[i] = (row[3])
        #print (row[0]), (row[1]), (row[2]), row[4], row[3], currSubId
        i += 1
    print "num hits = ", len(vecDIdent)

    curs.close()
    conn.close()    
    
    ####
    #### Plotting
    ####
    pylab.figure(1)
    ax = pylab.subplot(111)
    x = range(maxSEnd)
    c = ['r', 'b', 'g', 'k', 'm', 'y']
    
    print "gi = ", vecGi[0]
    for i in range(0, len(vecDIdent)):
        xs = vecStart[i]
        xe = vecEnd[i]
        x2 = [xs, xe]
        y2 = [vecDIdent[i], vecDIdent[i]]
        cIndex = 0
        
        if i > 0 and vecGi[i] != vecGi[i-1]:
            cIndex += 1
            if cIndex > len(c):
                cIndex = 0
            print "gi, colorIndex =", vecGi[i], cIndex
        line = lines.Line2D(x2, y2, lw=3., color=c[cIndex], alpha=0.4)
        ax.add_line(line)

    ax.set_ylim(0.0, 100.0)
    ax.set_xlim(minSStart, maxSEnd)
    #ax.set_xticks(x, minor=True)
    #ax.set_yticks(y, minor=True)

    ax.grid(color='k', linestyle='-.', linewidth=0.5)
    ax.yaxis.grid(True, linestyle='-.', which='minor')
    ax.xaxis.grid(True, linestyle='-.', which='minor')

    ax.set_ylabel('identity (%)', fontsize=20)
    ax.set_xlabel('subject genome sequence', fontsize=20)
    
    #ax.set_xticklabels(('0', '', '', '', '', '100'))
        
    fontsize=16
    for tick in ax.xaxis.get_major_ticks():
        tick.label1.set_fontsize(fontsize)
    for tick in ax.yaxis.get_major_ticks():
        tick.label1.set_fontsize(fontsize)
        
    #imFileName = dbName + "-indentity-plot.png"
    #pylab.savefig(imFileName, dpi=(300))
    pylab.show()
                
 
### EOF
    
