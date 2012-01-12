import pylab
import matplotlib.pyplot as plt

x = [12000, 40000, 80000]

#y1 = [599, 327, 195, 163, 163, 157]
#y2 = [587, 318, 192, 154, 155, 156]
#y3 = [1700, 828, 423, 229, 168, 159]
#y4 = [0,0,0,0,285,270]
#y5 = [4800,1552,756,386,214,160]
#y5 = [4856, 1648, 725, 386, 214, 160]
#z1 = [0,0,0,0,278,212]

y1000s12f = [587, 318, 192, 154, 155, 156]
y1000s40f = [1700, 828, 423, 229, 168, 159]
#y2000s40f = [3278, 1238, 574, 365, 285, 270]
y1000s80f = [4856, 1648, 725, 386, 214, 160]

i=0
y32c = [y1000s12f[i], y1000s40f[i], y1000s80f[i]]
i=1
y64c = [y1000s12f[i], y1000s40f[i], y1000s80f[i]]
i=2
y128c = [y1000s12f[i], y1000s40f[i], y1000s80f[i]]
i=3
y256c = [y1000s12f[i], y1000s40f[i], y1000s80f[i]]
i=4
y512c = [y1000s12f[i], y1000s40f[i], y1000s80f[i]]
i=5
y1024c = [y1000s12f[i], y1000s40f[i], y1000s80f[i]]

pylab.figure(1)
ax=pylab.subplot(111)

#ax.loglog(x, y4, 'md-.', linewidth=1.5, basex=2, basey=10, label='v4')
#ax.loglog(x, zold, 'bx--', linewidth=1.5, basex=2, basey=10, label='v1')
#ax.loglog(x, znew, 'gs-.', linewidth=1.5, basex=2, basey=10, label='v2')
ax.plot(x, y32c, 'md-.', linewidth=2, label='v4')
ax.plot(x, y64c, 'bv--', linewidth=2, label='v4')
ax.plot(x, y128c, 'gs-.', linewidth=2, label='v4')
ax.plot(x, y256c, 'ks:', linewidth=2, label='v4')
ax.plot(x, y512c, 'co-.', linewidth=2, label='v4')
ax.plot(x, y1024c, 'ro-', linewidth=2, label='v4')

#ax.loglog(x, y2000s40f, 'bx--', linewidth=1.5, basex=2, basey=10, label='v4')
#ax.semiogx(x, y1000s40f, 'rd-.', linewidth=2, basex=2, basey=10, label='v5')
#ax.semiogx(x, y1000s12f, 'ko-.', linewidth=2, basex=2, basey=10, label='v3')
#ax.loglog(x, y1, 'bx--', linewidth=1.5, basex=2, basey=10, label='v1')
#ax.loglog(x, y2, 'gs-.', linewidth=1.5, basex=2, basey=10, label='v2')
#ax.loglog(x, z1, '-.', linewidth=1.5, basex=2, basey=10, label='v2')

ax.set_ylim(0, 5000)
ax.set_xlim(10000, 82000)

#ax.set_xticks(x, minor=True)
#ax.set_yticks(y32c, minor=True)
ax.set_xticklabels(('12k', '', '', '40k', '', '', '', '80k'))
ax.yaxis.grid(True, linestyle='-.', which='major')
ax.xaxis.grid(True, linestyle='-.', which='major')
ax.yaxis.grid(True, linestyle='-.', which='minor')
ax.xaxis.grid(True, linestyle='-.', which='minor')

ax.set_xlabel('Dataset size  ', fontsize=20)
ax.set_ylabel('Computation time in minutes ', fontsize=20)
    
fontsize=16
for tick in ax.xaxis.get_major_ticks():
    tick.label1.set_fontsize(fontsize)
for tick in ax.yaxis.get_major_ticks():
    tick.label1.set_fontsize(fontsize)

#for i in range(0,len(x)):
    #ax.annotate(y1000s80f[i], xy=(x[i], y1000s80f[i]),  xycoords='data',
                #xytext=(2, 2), textcoords='offset points', fontsize=10, color='m')
                ##,arrowprops=dict(arrowstyle="->",connectionstyle="arc3,rad=.2")
                ##)
                
#for i in range(0,len(x)):
    #ax.annotate(y2000s40f[i], xy=(x[i], y2000s40f[i]),  xycoords='data',
                #xytext=(-21, -10), textcoords='offset points', fontsize=10, color='b')
                ###,arrowprops=dict(arrowstyle="->",connectionstyle="arc3,rad=.2")
                ###)
                
#for i in range(0,len(x)):
    #ax.annotate(y1000s40f[i], xy=(x[i], y1000s40f[i]),  xycoords='data',
                #xytext=(-21, -10), textcoords='offset points', fontsize=10, color='r')
                ###,arrowprops=dict(arrowstyle="->",connectionstyle="arc3,rad=.2")
                ###)
                
#for i in range(0,len(x)):
    #ax.annotate(y1000s12f[i], xy=(x[i], y1000s12f[i]),  xycoords='data',
                #xytext=(-21, -15), textcoords='offset points', fontsize=10, color='k')
                ###,arrowprops=dict(arrowstyle="->",connectionstyle="arc3,rad=.2")
                ###)
                
ax.legend( (
            '32 cores',
            '64 cores',
            '128 cores',
            '256 cores',
            '512 cores',
            '1024 cores'            
            ), 'upper left' )
            
leg = plt.gca().get_legend()
ltext  = leg.get_texts()
plt.setp(ltext, fontsize='18')    # the legend text fontsize

pylab.savefig('mrblast_plot4.png', dpi=(300))  
pylab.show()
