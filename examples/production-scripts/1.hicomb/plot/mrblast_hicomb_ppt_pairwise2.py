import pylab
import matplotlib.pyplot as plt

#x = [32, 64, 128, 256, 512, 1024]
#x2 = [32, 64, 128, 256, 512]
#x3 = [512]
#x4 = [32, 64, 128, 256, 512, 1024, 2048]
#x4 = [32, 64, 128, 256, 512, 1024]
x = [32, 64, 128, 256, 512, 1024]
 

### for HiComb PPT #############
#s40k_800k_20blocks = [0,0,0,0,0,4.27]
#s40k_400k_40blocks = [133, 60, 30.2, 17.81, 10.37, 6.71] # borrowed from s40k_1000x40_workitemrodering
#s40k_200k_80blocks = [0,0,0,0,0,12.95]
#s40k_100k_160blocks = [0,0,0,0,0,20.9]

### 40k, 120k, 200k seqs
for512cores = [10.37, 24.32, 43.31, 57.26]
for1024cores = [6.28, 16.26, 26.32, 34.26]
################################

pylab.figure(1)
ax=pylab.subplot(111)

lw = 1.5
li = 5
#ax.loglog(x, s40k_1000x40_2nd_scheduler, 'gx-', linewidth=3, basex=2, basey=10, label='v4', markersize=10)
#ax.loglog(x4, s40k_100k_160blocks, 'gv-', linewidth=2, basex=2, basey=10, label='v4', markersize=5)
#ax.loglog(x4, s40k_200k_80blocks, 'md-', linewidth=2, basex=2, basey=10, label='v4', markersize=5)
ax.semilogx(x[1:5], for512cores, 'rd-', linewidth=2,  basex=2, label='v4', markersize=5)
ax.semilogx(x[1:5], for1024cores, 'bo-', linewidth=2, basex=2,   label='v4', markersize=5)


#ax.set_ylim(0, 100)
ax.set_xlim(32, 1024)

ax.set_xticks(x, minor=True)
#ax.set_yticks(y, minor=True)

ax.set_xticklabels(('', '40k', '120k', '200k', '280k', ''))
ax.yaxis.grid(True, linestyle='-.', which='minor')
ax.xaxis.grid(True, linestyle='-.', which='minor')

ax.set_xlabel('Number of sequences', fontsize=20)
ax.set_ylabel('Computation time in minutes', fontsize=20)
    
fontsize=16
for tick in ax.xaxis.get_major_ticks():
    tick.label1.set_fontsize(fontsize)
for tick in ax.yaxis.get_major_ticks():
    tick.label1.set_fontsize(fontsize)

datalabelfontsize=13
#for i in range(1, 5):
    #ax.annotate(for512cores[i], xy=(x[i], for512cores[i]),  xycoords='data',
                #xytext=(2, -12), textcoords='offset points', fontsize=datalabelfontsize, color='r')
                ###,arrowprops=dict(arrowstyle="->",connectionstyle="arc3,rad=.2")
                ###)
#for i in range(1, 5):
    #ax.annotate(for1024cores[i], xy=(x[i], for1024cores[i]),  xycoords='data',
                #xytext=(2, -12), textcoords='offset points', fontsize=datalabelfontsize, color='b')
                ###,arrowprops=dict(arrowstyle="->",connectionstyle="arc3,rad=.2")
                ###)
                
ax.legend( (
            '512 cores',
            '1024 cores'
            ), loc=0)
            
leg = plt.gca().get_legend()
ltext  = leg.get_texts()
plt.setp(ltext, fontsize='14')    # the legend text fontsize

pylab.savefig('mrblast_plot_hicomb_ppt_pairwiseversion2.png', dpi=(300))                
pylab.show()
