import pylab
import matplotlib.pyplot as plt

x = [32, 64, 128, 256, 512, 1024]

#y1 = [599, 327, 195, 163, 163, 157]
#y2 = [587, 318, 192, 154, 155, 156]
#y3 = [1700, 828, 423, 229, 168, 159]
#y4 = [0,0,0,0,285,270]
#y5 = [4800,1552,756,386,214,160]
#y5 = [4856, 1648, 725, 386, 214, 160]
#z1 = [0,0,0,0,278,212]
#znew = [0,0,0,0,0,0]
#zold = [0,0,0,0,0,598]

y1000s12f = [587, 318, 192, 154, 155, 156]
y1000s40f = [1700, 828, 423, 229, 168, 159]
y2000s40f = [3278, 1238, 574, 365, 285, 270]
y1000s80f = [4856, 1648, 725, 386, 214, 160]

pylab.figure(1)
ax=pylab.subplot(111)

#ax.loglog(x, y4, 'md-.', linewidth=1.5, basex=2, basey=10, label='v4')
#ax.loglog(x, zold, 'bx--', linewidth=1.5, basex=2, basey=10, label='v1')
#ax.loglog(x, znew, 'gs-.', linewidth=1.5, basex=2, basey=10, label='v2')
ax.loglog(x, y1000s80f, 'md-', linewidth=2, basex=2, basey=10, label='v4')
ax.loglog(x, y2000s40f, 'bv--', linewidth=2, basex=2, basey=10, label='v4')
ax.loglog(x, y1000s40f, 'rd-.', linewidth=2, basex=2, basey=10, label='v5')
ax.loglog(x, y1000s12f, 'ko-.', linewidth=2, basex=2, basey=10, label='v3')
#ax.loglog(x, y1, 'bx--', linewidth=1.5, basex=2, basey=10, label='v1')
#ax.loglog(x, y2, 'gs-.', linewidth=1.5, basex=2, basey=10, label='v2')
#ax.loglog(x, z1, '-.', linewidth=1.5, basex=2, basey=10, label='v2')

ax.set_ylim(1e2, 1e4)
ax.set_xlim(16, 2048)

#ax.set_xticks(x, minor=True)
#ax.set_yticks(y, minor=True)
ax.set_xticklabels(('', '32','64','128','256','512','1024', ''))
ax.yaxis.grid(True, linestyle='-.', which='minor')
ax.xaxis.grid(True, linestyle='-.', which='minor')

ax.set_xlabel('Number of cores ($log_2$)', fontsize=20)
ax.set_ylabel('Computation time in minutes ($log_{10}$)', fontsize=20)
    
fontsize=16
for tick in ax.xaxis.get_major_ticks():
    tick.label1.set_fontsize(fontsize)
for tick in ax.yaxis.get_major_ticks():
    tick.label1.set_fontsize(fontsize)

datalabelfontsize=13
for i in range(0,len(x)):
    ax.annotate(y1000s80f[i], xy=(x[i], y1000s80f[i]),  xycoords='data',
                xytext=(2, 2), textcoords='offset points', fontsize=datalabelfontsize, color='m')
                #,arrowprops=dict(arrowstyle="->",connectionstyle="arc3,rad=.2")
                #)
                
for i in range(0,len(x)):
    ax.annotate(y2000s40f[i], xy=(x[i], y2000s40f[i]),  xycoords='data',
                xytext=(-24, -13), textcoords='offset points', fontsize=datalabelfontsize, color='b')
                ##,arrowprops=dict(arrowstyle="->",connectionstyle="arc3,rad=.2")
                ##)
                
for i in range(0,len(x)):
    ax.annotate(y1000s40f[i], xy=(x[i], y1000s40f[i]),  xycoords='data',
                xytext=(-24, -12), textcoords='offset points', fontsize=datalabelfontsize, color='r')
                ##,arrowprops=dict(arrowstyle="->",connectionstyle="arc3,rad=.2")
                ##)
                
for i in range(0,len(x)):
    ax.annotate(y1000s12f[i], xy=(x[i], y1000s12f[i]),  xycoords='data',
                xytext=(-21, -18), textcoords='offset points', fontsize=datalabelfontsize, color='k')
                ##,arrowprops=dict(arrowstyle="->",connectionstyle="arc3,rad=.2")
                ##)
                
ax.legend( (
            '80,000 sequences in 80 blocks',             
            '80,000 sequences in 40 blocks',             
            '40,000 sequences in 40 blocks',             
            '12,000 sequences in 12 blocks',               
            ) )
            
leg = plt.gca().get_legend()
ltext  = leg.get_texts()
plt.setp(ltext, fontsize='18')    # the legend text fontsize

pylab.savefig('mrblast_plot1.png', dpi=(300))                
pylab.show()
