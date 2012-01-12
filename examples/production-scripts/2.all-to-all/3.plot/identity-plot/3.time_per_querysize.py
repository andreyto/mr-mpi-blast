#!/usr/bin/env python
#from sqlite3 import *

import sys
import numpy as npy
import matplotlib.pyplot as plt
import pylab
 
if __name__ == '__main__':
     
    if len(sys.argv) != 1:
        print "python cov_*.py"
        sys.exit(1)    
    
    #dbName = sys.argv[1]
    #conn = connect(dbName)
    #curs = conn.cursor()
    #numSlices = 100
    #cov_busy = npy.zeros(numSlices,dtype=float)    
    
    cores = [1024, 1024, 1024, 2048]
    querysize = [600, 800, 1000, 1200]
    runtime = [54, 79, 95, 64]
    
    x = [1, 2, 3, 4]
    y = [cores[i]*runtime[i]/querysize[i] for i in range(4)]
    
    print y
    
 
 
    ####
    #### Plotting
    ####
    fig = pylab.figure(1)
    #pylab.figure(1)
    ax=pylab.subplot(111)
    #y = range(numSlices)

    ax.plot(x, y, 'ro', linewidth=1, markersize=5)

    ax.set_xlim(0, 5)
    ax.set_ylim(90, 110)

    ##ax.set_xticks(x, minor=True)
    ##ax.set_yticks(y, minor=True)

    ax.grid(color='k', linestyle='-.', linewidth=0.5)
    ax.yaxis.grid(True, linestyle='-.', which='minor')
    ax.xaxis.grid(True, linestyle='-.', which='minor')

    ax.set_ylabel('runtime*cores/querysize', fontsize=20)
    #ax.set_xlabel('t (%)', fontsize=20)
    
    ax.set_xticklabels(('', '600MB\n1024cores', '800MB\n1024cores','1000MB\n1024cores','1200MB\n2048cores', ''))
        
    fontsize=16
    for tick in ax.xaxis.get_major_ticks():
        tick.label1.set_fontsize(fontsize)
    for tick in ax.yaxis.get_major_ticks():
        tick.label1.set_fontsize(fontsize)
    
    #fig.autofmt_xdate()
    
    imFileName = "3.runtime-per-querysize.png"
    pylab.savefig(imFileName, dpi=(300))
    pylab.show()
    
### EOF

