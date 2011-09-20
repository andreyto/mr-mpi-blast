#!/usr/bin/env python
#from sqlite3 import *
import sys
import numpy as np
import matplotlib as plt
import pylab
#import random
from matplotlib import lines 
import tables as tb
import colorgradient as gc
from matplotlib.collections import LineCollection

print 'tables.__version__',tb.__version__

def iter5(tbl):
   keys = set(tbl.col('sId'))
   for _key in keys:
      rows = tbl.readWhere('key == _key')
      rows.sort(order = ['value'])
      for row in rows:
         print(row['key'], row['value'])
         
#def random_color():
    #COLOR_RANGE = (50, 255)
    #rgb = list()
    #for c in range(0, 3):
        #rgb.append(random.randrange(COLOR_RANGE[0], COLOR_RANGE[1]))
    #return "".join([hex(c)[2:].upper() for c in rgb])
 
        
#class BlHits(t.IsDescription):
    #gi         = t.UInt32Col()
    #sId        = t.UInt32Col()
    #qStart     = t.UInt32Col()
    #qEnd       = t.UInt32Col()
    #sStart     = t.UInt32Col()
    #sEnd       = t.UInt32Col()        
    #eValue     = t.FloatCol()
    #bitScore   = t.UInt32Col()        
    #upperStart = t.UInt32Col()
    #upperEnd   = t.UInt32Col()
    #dIdent     = t.Float32Col() 
    #dCover     = t.Float32Col() 
        
if __name__ == '__main__':
     
    if len(sys.argv) != 3:
        print "python 5.ident_plot_pytables.py fileName subjectIndex"
        sys.exit(1)  
        
    fileName = sys.argv[1]
    subjectIndex = int(sys.argv[2])
    coverCutoff = 0
    identCutoff = 0
    
    ### 
    ### HDF5 file open
    ### 
    h5file = tb.openFile(fileName, mode = "r")
    #root = h5file.root
    #group = h5file.createGroup(root, "blhits", "blhits")
    #table = h5file.createTable(group, "blhits", BlHits, "blhits")
    table = h5file.root.blhits.blhitstab
    #BlHits = table.row

    ### 
    ### Get unique subID
    ### 
    #vecSubId = []
    #sql = "select sId, count(*), gi from item \
           #where dIdent >= %s and dCover >= %s \
           #group by sId order by count(*) desc" \
           #% (identCutoff, coverCutoff)
    #curs.execute(sql)
    #for row in curs:
        #vecSubId.append(row[0])
        
    #colSid = table.cols.sId
    #uniqSid = dict((sid, set()) for sid in colSid)
    #print 'num unique sid = ', len(uniqSid)
    #currSubId = uniqSid.keys()[subjectIndex]
    #currSubId = table.cols.sId[subjectIndex]
    
    #colSid = table.col('sId')
    #colSid.sort()
    #dif = np.ones(colSid.shape,colSid.dtype)
    #dif[1:] = np.diff(colSid)
    #idx = np.where(dif>0)
    #vals = colSid[idx]
    #count = np.diff(idx)
    #print vals, count

    #freq = [ (a, colSid.count(a)) for a in set(colSid) ]
    #print sorted(freq, key=lambda x: -x[1]) 

    #scol = h5file.getNode('/blhits/_i_blhitstab/sId/sorted')
    #print scol
    #colSid = table.readSorted(sortby, 'sId')
    #print colSid
    
    currGi = 116249766
    currSubId = 116249766
    print "gi selected = ", currGi
    print "subject id selected = ", currSubId
    
    ###
    ### Get min, max and read columns
    ###
    #minSStart = min(table.cols.sStart)    
    #maxSEnd = max(table.cols.sEnd)    
    #print "min s.start, max s.end = ", minSStart, maxSEnd
    
    vecSStart = [ x['sStart'] for x in table.where('(gi == currGi) & (sId == currSubId)') ]
    vecSEnd   = [ x['sEnd']   for x in table.where('(gi == currGi) & (sId == currSubId)') ]
    vecDIdent = [ x['dIdent'] for x in table.where('(gi == currGi) & (sId == currSubId)') ]
    vecGi     = [ x['gi']     for x in table.where('(gi == currGi) & (sId == currSubId)') ]
    assert len(vecSStart) == len(vecSEnd) == len(vecDIdent) == len(vecGi) 
    print "num records = ", len(vecDIdent)
    minSStart = min(vecSStart)
    maxSEnd = max(vecSEnd)
    print "min s.start, max s.end = ", minSStart, maxSEnd
    
    h5file.close()
    
    ####
    #### Plotting
    ####
    fig = pylab.figure(1)
    ax = pylab.subplot(111)
    c = ['r', 'b', 'g', 'k', 'm', 'y']
    #colorCode = random_color()
    #cvalues.append("#"+str(colorCode))
    
    ###
    ### Gen gradient color code
    ###
    #steps = 10
    #generatedColors = gc.grad_colors("#0000FF", "#FF0000", steps)
    #print generatedColors
    
    print "gi = ", vecGi[0]
    cIndex = 0
    for i in range(0, len(vecDIdent), len(vecDIdent)/1000):
        xs = vecSStart[i]
        xe = vecSEnd[i]
        x2 = [xs, xe]
        y2 = [vecDIdent[i], vecDIdent[i]]
        #c = int(vecDIdent[i] / 10)
        line = lines.Line2D(x2, y2, lw=2., color=c[cIndex], alpha=0.5)
        ax.add_line(line)
        #ax.plot(x2, y2, 'r-', linewidth=3, markersize=2, alpha=0.5)

    ax.set_ylim(0.0, 100.0)
    ax.set_xlim(minSStart, maxSEnd)
    ax.set_xlim(0, 5057142)
    #ax.set_xticks(x, minor=True)
    #ax.set_yticks(y, minor=True)

    ax.grid(color='k', linestyle='-.', linewidth=0.5)
    ax.yaxis.grid(True, linestyle='-.', which='minor')
    ax.xaxis.grid(True, linestyle='-.', which='minor')

    #ax.set_ylabel('identity (%)', fontsize=20)
    #ax.set_xlabel('subject genome sequence', fontsize=20)
    
    ax.set_xticklabels(('0', '1000000', '2000000', '3000000', '4000000', '5000000'))        
    
    ### Set lable font size
    fontsize=16
    for tick in ax.xaxis.get_major_ticks():
        tick.label1.set_fontsize(fontsize)
    for tick in ax.yaxis.get_major_ticks():
        tick.label1.set_fontsize(fontsize)
    
    ### Rotate xaxis label
    #for label in ax.get_xticklabels():
        #pylab.setp(label, rotation=45)
    #fig.autofmt_xdate()
    
    ### Save fig
    imgFileName = str(currSubId) + ".png"
    pylab.savefig(imgFileName, dpi=(600))
    
    ### Show fig
    pylab.show()
 
### EOF
    
