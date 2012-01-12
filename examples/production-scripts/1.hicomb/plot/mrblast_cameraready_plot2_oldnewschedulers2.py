import pylab
import matplotlib.pyplot as plt

x = [32, 64, 128, 256, 512, 1024]
x2 = [32, 64, 128, 256, 512]
x3 = [512]
x4 = [32, 64, 128, 256, 512, 1024, 2048]

#y1000s12f = [587, 318, 192, 154, 155, 156]
#y1000s40f = [1700, 828, 423, 229, 168, 159]
#y2000s40f = [3278, 1238, 574, 365, 285, 270]
#y1000s80f = [4856, 1648, 725, 386, 214, 160]

s6k_new = [120.32, 114.26, 94.3, 64, 34.54]
s10k_new = [200, 77, 32.43, 22, 23]
s20k_new = [220, 196.25, 160.14, 118.57, 66]
s30k_new = [251, 167, 92.32, 47.49, 22]
s40k_new_1000x40 = [255.52, 180.45, 160.49, 125, 66 ]
s40k_new_500x80 = [319.37, 240, 206.53, 129.33, 61 ]
s40k_new_1000x40_2 = [109.65,42.39,22.08,12.26,10.67,11.08]
s40k_new_250x160_2 = [50.76]

##############################
s40k_1000x40_workitemrodering = [133, 60, 30.2, 17.81, 10.37, 6.71, 4.64]
#s40k_1000x40_2nd_scheduler = [0, 0, 31.63, 17.48, 12.6, 8.87]
s40k_1000x40_2nd_scheduler = [0, 0, 29.77, 17.23, 11.15, 7.65]
##############################

s80k_new_1000x80_2 = [23.43]
s160k_new_1000x160_2 = [51.21]


s6k_old = [148.33, 132, 113.58, 92.25, 60]
s10k_old = [209.18, 88.29, 43.21, 27.2, 28.1 ]
s20k_old = [270, 243.2, 185.52, 135.25, 77.22 ]
s30k_old = [286.2, 190, 112.16, 66, 28.1 ]
s40k_old_1000x40 = [319.35, 227.9, 185.3, 162.58, 77.3 ]
s40k_old_500x80 = [382, 282.9, 240, 183.47, 79.54 ]
s40k_old_1000x40_2 = [152.44,65.87,29.34,17.29,12.88,15.2]

pylab.figure(1)
ax=pylab.subplot(111)

y1000s12f = [587, 318, 192, 154, 155, 156]
y1000s40f = [1700, 828, 423, 229, 168, 159]
y2000s40f = [3278, 1238, 574, 365, 285, 270]
y1000s80f = [4856, 1648, 725, 386, 214, 160]

lw = 1.5

#ax.loglog(x, y1000s80f, 'md--', linewidth=lw, basex=2, basey=10, label='v4')
#ax.loglog(x, y2000s40f, 'bv--', linewidth=lw, basex=2, basey=10, label='v4')
#ax.loglog(x, y1000s40f, 'rd--', linewidth=lw, basex=2, basey=10, label='v5')
#ax.loglog(x, y1000s12f, 'ko--', linewidth=lw, basex=2, basey=10, label='v3')

li = 5
#ax.loglog(x2, s40k_old_500x80[:li], 'md-.', linewidth=lw, basex=2, basey=10 )
#ax.loglog(x2, s40k_old_1000x40[:li], 'bv-.', linewidth=lw, basex=2, basey=10 )
##ax.loglog(x2, s30k_old[:li], 'bv-.', linewidth=lw, basex=2, basey=10 )
#ax.loglog(x2, s20k_old[:li], 'rd-.', linewidth=lw, basex=2, basey=10)
##ax.loglog(x2, s10k_old[:li], 'ko-.', linewidth=lw, basex=2, basey=10, label='v3')
#ax.loglog(x2, s6k_old[:li], 'ko-.', linewidth=lw, basex=2, basey=10, label='v3')

#ax.loglog(x2, s40k_new_500x80[:li], 'md-', linewidth=lw, basex=2, basey=10, label='v4')
#ax.loglog(x2, s40k_new_1000x40[:li], 'bv-', linewidth=lw, basex=2, basey=10, label='v4')
##ax.loglog(x2, s30k_new[:li], 'gv-', linewidth=lw, basex=2, basey=10, label='v4')
#ax.loglog(x2, s20k_new[:li], 'rd-', linewidth=lw, basex=2, basey=10, label='v5')
##ax.loglog(x2, s10k_new[:li], 'o-', linewidth=lw, basex=2, basey=10, label='v3')
#ax.loglog(x2, s6k_new[:li], 'ko-', linewidth=lw, basex=2, basey=10, label='v3')

#ax.loglog(x, s40k_old_1000x40_2, 'bv-.', linewidth=lw, basex=2, basey=10 )
#ax.loglog(x, s40k_new_1000x40_2, 'bv-', linewidth=lw, basex=2, basey=10, label='v4')

#ax.loglog(x3, s40k_new_250x160_2, 'rd-', linewidth=lw, basex=2, basey=10, label='v4', markersize=10)
#ax.loglog(x3, s80k_new_1000x80_2, 'go-', linewidth=lw, basex=2, basey=10, label='v4', markersize=10)
#ax.loglog(x3, s160k_new_1000x160_2, 'cv-', linewidth=lw, basex=2, basey=10, label='v4', markersize=10)

ax.loglog(x4, s40k_1000x40_workitemrodering, 'ro-', linewidth=2, basex=2, basey=10, label='v4', markersize=5)
#ax.loglog(x, s40k_1000x40_2nd_scheduler, 'gx-', linewidth=3, basex=2, basey=10, label='v4', markersize=10)



ax.set_ylim(1e0, 1e3)
ax.set_xlim(16, 4096)

#ax.set_xticks(x, minor=True)
#ax.set_yticks(y, minor=True)
ax.set_xticklabels(('', '32','64','128','256','512','1024', '2048', ''))
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
#for i in range(0,len(x)):
    #ax.annotate(y1000s80f[i], xy=(x[i], y1000s80f[i]),  xycoords='data',
                #xytext=(2, 2), textcoords='offset points', fontsize=datalabelfontsize, color='m')
                ##,arrowprops=dict(arrowstyle="->",connectionstyle="arc3,rad=.2")
                ##)
#for i in range(0,len(x2)):
    #ax.annotate(s40k_new_500x80[i], xy=(x2[i], s40k_new_500x80[i]),  xycoords='data',
                #xytext=(2, 2), textcoords='offset points', fontsize=datalabelfontsize, color='m')
                ##,arrowprops=dict(arrowstyle="->",connectionstyle="arc3,rad=.2")
                ##)
                                
#for i in range(0,len(x)):
    #ax.annotate(y2000s40f[i], xy=(x[i], y2000s40f[i]),  xycoords='data',
                #xytext=(-24, -13), textcoords='offset points', fontsize=datalabelfontsize, color='b')
                ###,arrowprops=dict(arrowstyle="->",connectionstyle="arc3,rad=.2")
                ###)
#for i in range(0,len(x2)):
    #ax.annotate(s40k_new_1000x40[i], xy=(x2[i], s40k_new_1000x40[i]),  xycoords='data',
                #xytext=(-24, -13), textcoords='offset points', fontsize=datalabelfontsize, color='b')
                ###,arrowprops=dict(arrowstyle="->",connectionstyle="arc3,rad=.2")
                ###)                
for i in range(0,len(x4)):
    ax.annotate(s40k_1000x40_workitemrodering[i], xy=(x4[i], s40k_1000x40_workitemrodering[i]),  xycoords='data',
                xytext=(2, 2), textcoords='offset points', fontsize=datalabelfontsize, color='r')
                ##,arrowprops=dict(arrowstyle="->",connectionstyle="arc3,rad=.2")
                ##)    
                                
#for i in range(0,len(x)):
    #ax.annotate(y1000s40f[i], xy=(x[i], y1000s40f[i]),  xycoords='data',
                #xytext=(-24, -12), textcoords='offset points', fontsize=datalabelfontsize, color='r')
                ###,arrowprops=dict(arrowstyle="->",connectionstyle="arc3,rad=.2")
                ###)
                
#for i in range(0,len(x)):
    #ax.annotate(y1000s12f[i], xy=(x[i], y1000s12f[i]),  xycoords='data',
                #xytext=(-21, -18), textcoords='offset points', fontsize=datalabelfontsize, color='k')
                ###,arrowprops=dict(arrowstyle="->",connectionstyle="arc3,rad=.2")
                ###)
                
#for i in range(0,len(x)):
    #ax.annotate(s40k_old_1000x40_2[i], xy=(x[i], s40k_old_1000x40_2[i]),  xycoords='data',
                #xytext=(-24, -13), textcoords='offset points', fontsize=datalabelfontsize, color='b')
                ###,arrowprops=dict(arrowstyle="->",connectionstyle="arc3,rad=.2")
                ###)
#for i in range(0,len(x)):
    #ax.annotate(s40k_new_1000x40_2[i], xy=(x[i], s40k_new_1000x40_2[i]),  xycoords='data',
                #xytext=(-24, -13), textcoords='offset points', fontsize=datalabelfontsize, color='b')
                ###,arrowprops=dict(arrowstyle="->",connectionstyle="arc3,rad=.2")
                ###)      
            
                
#ax.legend( (
            ##'40k sequences in 80 blocks (hicomb)',             
            ##'40k sequences in 40 blocks (hicomb)',             
            ##'20k sequences in 40 blocks (hicomb)',             
            ##'6k  sequences in 12 blocks (hicomb)',  
            ##'40k sequences in 80 blocks (old)',             
            ##'40k sequences in 40 blocks (old)',             
            ##'20k sequences in 40 blocks (old)',             
            ##'6k  sequences in 12 blocks (old)',                           
            ##'40k sequences in 80 blocks (new)',             
            ##'40k sequences in 80 blocks (new)',             
            ##'20k sequences in 40 blocks (new)',             
            ##'6k  sequences in 12 blocks (new)', 
            ##'40k sequences in 40 blocks (old2)',                  
            ##'40k sequences in 40 blocks (new2)',
            ##'40k sequences in 160 blocks (new2)',
            ##'80k sequences in 80 blocks (new2)',
            ##'160k sequences in 160 blocks (new2)',
            #'40k sequences in 40 blocks (new3-1)',
            ##'40k sequences in 40 blocks (new3-2)'
            #), loc=1 )
            
#leg = plt.gca().get_legend()
#ltext  = leg.get_texts()
#plt.setp(ltext, fontsize='9')    # the legend text fontsize

pylab.savefig('mrblast_plot2.png', dpi=(300))                
pylab.show()
