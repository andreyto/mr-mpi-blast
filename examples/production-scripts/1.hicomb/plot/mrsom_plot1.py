import pylab
import matplotlib.pyplot as plt

x = [32, 64, 128, 256, 512, 1024]
#y2 = [123.71, 62.28, 33.9, 19.4, 13.14, 11.18]
#y1 = [129.31, 61.37, 32.9, 17.24, 11.23, 6.44]
y3 = [0, 67.24, 32.45, 0, 0, 5.51]
#est = [200.1, 100.698, 50.67, 25.645, 13.136, 6.89]
est = [183, 92.272, 46.47, 23.57, 12.127, 6.4]
#y2 = [587, 318, 192, 154, 155, 156]
#y3 = [1700, 828, 423, 229, 168, 159]
#y4 = [0,0,0,0,285,270]
#y5 = [4856,1648,725,386,214,160]

y1 = [122.45, 61.57, 30.9, 15.55, 8.1, 4.17]
y2 = [123.63, 61.8, 30.85, 15.57, 7.87, 4.02]

pylab.figure(1)
ax=pylab.subplot(111)

#ax.plot(x1700, y, 'ro-', linewidth=1.5)
#ax.set_xscale('log', basex=2, label_minor=True)
#ax.set_yscale('log', basex=10, label_minor=True)


#ax.loglog(x, est, 'md-', linewidth=1.5, basex=2, basey=10, label='est')
#ax.loglog(x, y5, 'rd-', linewidth=1.5, basex=2, basey=10, label='v5')
#ax.loglog(x, y1, 'gs-.', linewidth=2, basex=2, basey=10, label='v2')
ax.loglog(x, y2, 'ro-', linewidth=2, basex=2, basey=10, label='v1')
#ax.loglog(x, y3, 'ko-.', linewidth=1.5, basex=2, basey=10, label='v3')

#ax.set_ylim(1e2, 1e4)
ax.set_xlim(16, 2048)

#ax.set_xticks(x, minor=True)
#ax.set_yticks(y, minor=True)
ax.set_xticklabels(('', '32','64','128','256','512','1024', ''))
ax.yaxis.grid(True, linestyle='-.', which='minor')
ax.xaxis.grid(True, linestyle='-.', which='minor')

ax.set_xlabel('Number of cores ($log_2$)', fontsize=20)
ax.set_ylabel('Running time in minutes ($log_{10}$)', fontsize=20)
    
fontsize=16
for tick in ax.xaxis.get_major_ticks():
    tick.label1.set_fontsize(fontsize)
fontsize=16
for tick in ax.yaxis.get_major_ticks():
    tick.label1.set_fontsize(fontsize)

for i in range(0,len(x)):
    ax.annotate(y2[i], xy=(x[i], y2[i]),  xycoords='data',
                xytext=(2, 2), textcoords='offset points', fontsize=15, color='r')
                ##,arrowprops=dict(arrowstyle="->",connectionstyle="arc3,rad=.2")
                ##)
                
#for i in range(0,len(x)):
    #ax.annotate(y2[i], xy=(x[i], y2[i]),  xycoords='data',
                #xytext=(2, 5), textcoords='offset points', fontsize=10, color='g')
                ##,arrowprops=dict(arrowstyle="->",connectionstyle="arc3,rad=.2")
                ##)
                
#for i in range(0,len(x)):
    #ax.annotate(y5[i], xy=(x[i], y5[i]),  xycoords='data',
                #xytext=(2, 2), textcoords='offset points', fontsize=10, color='r')
                ##,arrowprops=dict(arrowstyle="->",connectionstyle="arc3,rad=.2")
                ##)


#ax.legend( (
            #'81,920 vectors in 1,024 blocks', 
            #'81,920 vectors in 2,048 blocks'
            #'4,096 work items'
            #'1000x12 query files, sort in workers', 
            #'1000x12 query files, sort in master'            
            #) )
            
#leg = plt.gca().get_legend()
#ltext  = leg.get_texts()
#plt.setp(ltext, fontsize='18')    # the legend text fontsize

pylab.savefig('mrsom_plot1.png', dpi=300)            
pylab.show()
