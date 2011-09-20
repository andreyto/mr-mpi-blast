#!/usr/bin/env python
#from sqlite3 import *
import sys
import numpy as numpy
import matplotlib as plt
import pylab
import random
from matplotlib import lines 
import tables as tb
import colorgradient as gc
from matplotlib.collections import LineCollection
from matplotlib.colors import colorConverter
import time
from pylab import *
from numpy import *
from matplotlib.collections import LineCollection
import pdb
 
print 'tables.__version__',tb.__version__

if sys.platform == "win32":
     # On Windows, the best timer is time.clock()
     default_timer = time.clock
else:    
    # On most other platforms the best timer is time.time()    
    default_timer = time.time 
        
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
        print "python 8.evalue_plot.py fileName queryLen"
        sys.exit(1)  
        
    fileName = sys.argv[1]
    #subjectIndex = int(sys.argv[2])
    queryLen = int(sys.argv[2])   
    coverCutoff = 0
    identCutoff = 0
    
    ### 
    ### HDF5 file open
    ### 
    h5file = tb.openFile(fileName, mode = "r")
    table = h5file.root.blhits.blhitstab

    ### 
    ### Get unique subID
    ### 
    #colSid = table.cols.sId 
    #currSubId = table.cols.sId[subjectIndex]
    
    currGi = 116249766    
    currSubId = 116249766 #1
    #currSubId = 190889639 #2
    #currSubId = 319779749 #3
    #currSubId = 312112794 #4
    #currSubId = 217976200 #5import matplotlib as plt    
    print "gi selected = ", currGi
    print "subject id selected = ", currSubId
    
    ###
    ### Get min, max and read columns
    ###    
    #vec = [ [x['sStart'],x['sEnd'],x['dIdent'],x['gi'] ] \
    vec = [ [ x['dIdent'],x['bitScore'],x['eValue'] ] \
    #vec = [ [ x['dIdent'],x['bitScore'] ] \
            for x in table.where('(gi == currGi) & (sId == currSubId)') ]
    print "num records = ", len(vec)
    #minSStart = min([r[0] for r in vec])
    #maxSEnd = max([r[1] for r in vec])
    #print "min s.start, max s.end = ", minSStart, maxSEnd
    
    ### 
    ### 
    ###
    print "query len = ", queryLen
    yy = [r[0]*queryLen/100 for r in vec] ## # identical bases
    #yy = [r[0] for r in vec]              ## perc identity
    #xx = [r[1] for r in vec]              ## bitscore
    xx2 = [r[2] for r in vec]             ## evalue
    #print yy
    
    ###
    ###
    ###
    f = figure()
    ax = gca()
    #ax.plot(xx, yy, 'ro', markersize=1)
    ax.plot(xx2, yy, 'ro', markersize=1)
    ax.set_ylim(0, queryLen)
    ax.set_xlim(-0.000002, )
    ax.set_xlabel('evalue', fontsize=16)
    ax.set_ylabel('Number of identical bases', fontsize=16)
    fontsize=16
    for tick in ax.xaxis.get_major_ticks():
        tick.label1.set_fontsize(fontsize)
    for tick in ax.yaxis.get_major_ticks():
        tick.label1.set_fontsize(fontsize)
    ax.grid(color='k', linestyle='-.', linewidth=0.5)
    ax.yaxis.grid(True, linestyle='-.', which='minor')
    ax.xaxis.grid(True, linestyle='-.', which='minor')
    f.autofmt_xdate()
    #savefig(str(queryLen)+"bp-evalue-identbases.png", dpi=(300))                
    show() 
    draw()  

    
    h5file.close()
    
### EOF
    
 

