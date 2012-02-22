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
import time

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
        print "python 5-2.identity_plot.py fileName subjectIndex"
        sys.exit(1)  
        
    fileName = sys.argv[1]
    subjectIndex = int(sys.argv[2])
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
    currGi = 196259856

    #currSubId = 116249766 #1
    #currSubId = 190889639 #2
    #currSubId = 319779749 #3
    #currSubId = 312112794 #4
    #currSubId = 217976200 #5
    currSubId = subjectIndex
    print "gi selected = ", currGi
    print "subject id selected = ", currSubId
    
    ###
    ### Get min, max and read columns
    ###  
    #vec = [ [x['sStart'],x['sEnd'],x['dIdent'],x['gi'] ] \
        #for x in table.where('(sId == currSubId)') ]
    vec = [ [x['sStart'],x['sEnd'],x['dIdent'] ] \
        for x in table.where('(gi == currGi) & (sId == currSubId)') ]
    print "num records = ", len(vec)
    minSStart = min([r[0] for r in vec])
    maxSEnd = max([r[1] for r in vec])
    print "min s.start, max s.end = ", minSStart, maxSEnd
    
    h5file.close()
    
    ####
    #### Plotting
    ####
    fig = pylab.figure(1)
    ax = pylab.subplot(111)
    pylab.interactive(False)
    
    ###
    ### Gen gradient color code
    ###
    t0 = default_timer()
    
    cIndex = 0
    xpairs = []
    ypairs = []
    for i in range(0, len(vec)):
        xs = vec[i][0]
        xe = vec[i][1]
        y  = vec[i][2]
        x2 = [xs, xe]
        y2 = [y, y]
        
        ## 3.
        #xpairs.append(x2)   
        #ypairs.append(y2)
        
        ## 2.
        line = lines.Line2D(x2, y2, lw=2., color='r', alpha=0.5)
        ax.add_line(line)
        
        ## 1.
        #ax.plot(x2, y2, 'r-', linewidth=3, markersize=2, alpha=0.5)
    
        if i % 1000 == 0:
            print i
    
    ## 1.
    #for xends,yends in zip(xpairs,ypairs):
        #ax.plot(xends,yends,'r-',linewidth=3,markersize=2,alpha=0.5)
    
    ## 2. 
    #call_list = []
    #for xends,yends in zip(xpairs,ypairs):
        #call_list.append(xends)
        #call_list.append(yends)
        #call_list.append('r-')
    #ax.plot(*call_list,linewidth=2,alpha=0.5)
    
    ## 3.
    #xlist = []
    #ylist = []
    #for xends,yends in zip(xpairs,ypairs):
        #xlist.extend(xends)
        #xlist.append(None)
        #ylist.extend(yends)
        #ylist.append(None)
    #ax.plot(xlist, ylist, 'r', linewidth=3, markersize=2, alpha=0.5)
        
    ## 4.
    #call_list = []
    #for xends,yends in zip(xpairs,ypairs):
        #call_list.append(xends)
        #call_list.append(yends)
        #call_list.append('r-')
    #ax.plot(*call_list,linewidth=3,markersize=2,alpha=0.5)
    
    
    t1 = default_timer()
    
    print 'Plotting time = %f sec' % (t1-t0)  
    
    ax.set_ylim(0.0, 100.0)
    #ax.set_xlim(minSStart, maxSEnd)
    #ax.set_xlim(0, 5057142)
    ax.set_xlim(0, 6000000)
    #ax.set_xticks(x, minor=True)
    #ax.set_yticks(y, minor=True)

    ax.grid(color='k', linestyle='-.', linewidth=0.5)
    ax.yaxis.grid(True, linestyle='-.', which='minor')
    ax.xaxis.grid(True, linestyle='-.', which='minor')

    #ax.set_ylabel('identity (%)', fontsize=20)
    #ax.set_xlabel('subject genome sequence', fontsize=20)
    #ax.set_xticklabels(('0', '1000000', '2000000', '3000000', '4000000', '5000000'))        
    
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
    #imgFileName = str(currSubId) + ".png"
    #pylab.savefig(imgFileName, dpi=300)
    
    ### Show fig
    pylab.show()
 
### EOF
    
