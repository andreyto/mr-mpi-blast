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
     
    if len(sys.argv) != 2:
        print "python 5.identity_plot.py fileName"
        sys.exit(1)  
        
    fileName = sys.argv[1]
    #subjectIndex = int(sys.argv[2])
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
    #currSubId = 217976200 #5
    
    print "gi selected = ", currGi
    print "subject id selected = ", currSubId
    
    ###
    ### Get min, max and read columns
    ###
    
    vec = [ [x['sStart'],x['sEnd'],x['dIdent'],x['gi'] ] \
        for x in table.where('(gi == currGi) & (sId == currSubId)') ]
    print "num records = ", len(vec)
    minSStart = min([r[0] for r in vec])
    maxSEnd = max([r[1] for r in vec])
    print "min s.start, max s.end = ", minSStart, maxSEnd

    ###
    ### histogram
    ###
    query_len = 100
    num_queries = 101185
    total_query_len = num_queries * query_len
     
    numbins = 6
    mybins = np.linspace(0, 100, numbins)
    htemp, jnk = np.histogram([r[2] for r in vec], mybins)
    #print htemp
    for x in htemp:
        print '%d'%x,
    print ""
    
    hist_query_length = htemp * query_len / float(total_query_len)
    hist_num_query = htemp / float(num_queries)
    #print [float(x) for x in hist_query_length]
    #print hist_num_query
    for x in hist_query_length:
        print '%.3f'%x,
    
    
    h5file.close()
    
### EOF
    
