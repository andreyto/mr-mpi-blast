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
    
def iter5(tbl):
   keys = set(tbl.col('sId'))
   for _key in keys:
      rows = tbl.readWhere('key == _key')
      rows.sort(order = ['value'])
      for row in rows:
         print(row['key'], row['value'])
        
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
        print "python 5.identity_plot.py fileName subjectIndex"
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
    colSid = table.cols.sId 
    currSubId = table.cols.sId[subjectIndex]
 
    
    currGi = 116249766 #1    
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
        #for x in table.where('(sId == currSubId)') ]
    print "num records = ", len(vec)
    minSStart = min([r[0] for r in vec])
    maxSEnd = max([r[1] for r in vec])
    print "min s.start, max s.end = ", minSStart, maxSEnd
    
       
    ####
    #### Plotting
    ####
    #fig = pylab.figure(1)
    #ax = pylab.subplot(111)
    c = ['r', 'b', 'g', 'k', 'm', 'y']
    #colorCode = random_color()
    #cvalues.appebnd("#"+str(colorCode))
    
    ###
    ### Gen gradient color code
    ###
    #steps = 10
    #generatedColors = gc.grad_colors("#0000FF", "#FF0000", steps)
    #print generatedColors
    
    #print "gi = ", vecGi[0]
    #print "gi = ", vec[0]
    cIndex = 0
    f = figure()
    ax = gca()
    #vec = np.random.random((10,3)) * 100
    segs = []
    #pdb.set_trace()

#####################################################################
    #for i in range(0, len(vec)):
        #x1 = vec[i][0] 
        ##x1 = 3000000
        #x2 = vec[i][1] 
        ##x2 = 4000000
        #y = vec[i][2] / 100.0
        #segs.extend( [[(x1,y),(x2,y)]] )
     
    #line_segments = LineCollection(segs, linewidth=3, alpha=0.3, colors ='r', linestyle = 'solid')
    #ax.add_collection(line_segments)

#####################################################################
    xpairs = []
    ypairs = []
    for i in range(0, len(vec)):
        ## 4. 
        xs = float(vec[i][0]) / maxSEnd
        xe = float(vec[i][1]) / maxSEnd
        y  = float(vec[i][2])
        x2 = [xs, xe]
        y2 = [y, y]
        
        ## 3.
        xpairs.append(x2)   
        ypairs.append(y2)
        
        #if i % 1000 == 0:
            #print i
    
    ## 3.
    xlist = []
    ylist = []
    for xends,yends in zip(xpairs,ypairs):
        xlist.extend(xends)
        xlist.append(None)
        ylist.extend(yends)
        ylist.append(None)
    ax.plot(xlist, ylist, 'r', linewidth=2, alpha=0.5)
    #line = lines.Line2D(xlist, ylist, lw=2., color='r', alpha=0.5)
    #line.set_alpha(0.5)    
    #ax.add_line(line)
#####################################################################
    
    #ax.set_xlim(0,1)
    #ax.set_ylim(0,1)
    
    ax.set_ylim(0.0, 100.0)
    #ax.set_xlim(minSStart, maxSEnd)
    #ax.set_xlim(0, float(5057142)/maxSEnd)
    #ax.set_xticks(x, minor=True)
    #ax.set_yticks(y, minor=True)

    ax.grid(color='k', linestyle='-.', linewidth=0.1)
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
    #f.autofmt_xdate()
    
    show() 
    draw()  
        
    h5file.close()
        
    
    '''
    segs = []
 

    for i in range(0, len(vec)):
                
        ## 4. 
        x1 = vec[i][0]
        #x1 = 3000000
        x2 = vec[i][1]
        #x2 = 4000000
        y1 = y2 = vec[i][2]
        #y1 = y2 = 99
        
        segs.extend( [[(x1,y1),(x2,y2)]] )
        
        ## 3.
        #xpairs.append(x2)   
        #ypairs.append(y2)
        
        ## 2.
        #line = lines.Line2D(x2, y2, lw=2., color=c[cIndex], alpha=0.5)
        #ax.add_line(line)
        
        ## 1.
        #ax.plot(x2, y2, 'r-', linewidth=3, markersize=2, alpha=0.5)
    
    t0 = default_timer()
    
    print len(segs)
    print segs[:10]
    
    ## 4.
    line_segments = LineCollection(segs, linewidth=4, alpha=0.3,
                                colors = 'r',
                                linestyle = 'solid')
    ax.add_collection(line_segments)
    #ax.add_collection( LineCollection( [[(1000,91),(5000000,91)]], colors='r', linewidths=4 ) )      
    
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
        ##line = lines.Line2D(xends, yends, lw=2., color=c[cIndex], alpha=0.5)
        ##ax.add_line(line)
        
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
    ax.set_xlim(minSStart, maxSEnd)
    #ax.set_xlim(0, 5057142)
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
    fig.autofmt_xdate()
    
    ### Save fig
    #imgFileName = str(currSubId) + ".png"
    imgFileName = "test.png"
    pylab.savefig(imgFileName, dpi=(600))
    
    ### Show fig
    pylab.show()
    pylab.draw() # force a draw; now it works

'''

### EOF
    
 
