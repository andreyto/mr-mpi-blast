import pylab
import matplotlib.pyplot as plt

x = [32, 64, 128, 256, 512, 1024]
x2 = [32, 64, 128, 256, 512]
x3 = [512]
#x4 = [32, 64, 128, 256, 512, 1024, 2048]
x4 = [32, 64, 128, 256, 512, 1024]

#y1000s12f = [587, 318, 192, 154, 155, 156]
#y1000s40f = [1700, 828, 423, 229, 168, 159]
y2000s40f = [3278, 1238, 574, 365, 285, 270]
y1000s80f = [4856, 1648, 725, 386, 214, 160]

s6k_new = [120.32, 114.26, 94.3, 64, 34.54]
s10k_new = [200, 77, 32.43, 22, 23]
s20k_new = [220, 196.25, 160.14, 118.57, 66]
s30k_new = [251, 167, 92.32, 47.49, 22]
s40k_new_1000x40 = [255.52, 180.45, 160.49, 125, 66 ]
s40k_new_500x80 = [319.37, 240, 206.53, 129.33, 61 ]
s40k_new_1000x40_2 = [109.65,42.39,22.08,12.26,10.67,11.08]
s40k_new_250x160_2 = [50.76]

##############################
#s40k_1000x40_workitemrodering = [133, 60, 30.2, 17.81, 10.37, 6.71, 4.64]
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

y1000s12f = [587, 318, 192, 154, 155, 156]
y1000s40f = [1700, 828, 423, 229, 168, 159]
y2000s40f = [3278, 1238, 574, 365, 285, 270]
y1000s80f = [4856, 1648, 725, 386, 214, 160]

### for HiComb PPT #############
s40k_16000k_1blocks = [30, 20.14, 15.16, 15.2, 15.45, 15.2]
s40k_800k_20blocks = [85.2, 42.07, 20.5, 11.53, 7.2, 4.5]
s40k_400k_40blocks = [133, 60.3, 30.2, 18.22, 10.37, 6.28] # borrowed from s40k_1000x40_workitemrodering
s40k_200k_80blocks = [277, 133, 58.1, 28.56, 20.5, 13.48]
s40k_200k_120blocks = [0, 0, 77.18, 50.51, 24.32, 12.25]
s40k_100k_160blocks = [642, 300, 141.2, 70.58, 40.1, 21.19]

#s40k_400k_40blocks_w_m2 = [
s120k_400k_120blocks = [0, 0, 77.18, 41.53, 24.32, 16.26]
s120k_1200k_40blocks = [0, 0, 41.32, 0, 12.18, 7.52]

################################

pylab.figure(1)
ax=pylab.subplot(111)
lw = 1.5
li = 5
 
#ax.loglog(x, s40k_new_1000x40_2, 'bv-', linewidth=lw, basex=2, basey=10, label='v4')

#ax.loglog(x3, s40k_new_250x160_2, 'rd-', linewidth=lw, basex=2, basey=10, label='v4', markersize=10)
#ax.loglog(x3, s80k_new_1000x80_2, 'go-', linewidth=lw, basex=2, basey=10, label='v4', markersize=10)
#ax.loglog(x3, s160k_new_1000x160_2, 'cv-', linewidth=lw, basex=2, basey=10, label='v4', markersize=10)

#ax.loglog(x, s40k_1000x40_2nd_scheduler, 'gx-', linewidth=3, basex=2, basey=10, label='v4', markersize=10)
ax.loglog(x4, s40k_100k_160blocks, 'gv-', linewidth=2, basex=2, basey=10, label='v4', markersize=5)
ax.loglog(x4, s40k_200k_80blocks, 'md-', linewidth=2, basex=2, basey=10, label='v4', markersize=5)
#ax.loglog(x4, s120k_400k_120blocks, 'bv--', linewidth=2, basex=2, basey=10, label='v4', markersize=5)
#ax.loglog(x4, s40k_200k_120blocks, 'rx--', linewidth=2, basex=2, basey=10, label='v4', markersize=5)
#ax.loglog(x4, s120k_1200k_40blocks, 'bo-', linewidth=2, basex=2, basey=10, label='v4', markersize=5)
ax.loglog(x4, s40k_400k_40blocks, 'rx-', linewidth=2, basex=2, basey=10, label='v4', markersize=5)
ax.loglog(x4, s40k_800k_20blocks, 'bo-', linewidth=2, basex=2, basey=10, label='v4', markersize=5)
ax.loglog(x4, s40k_16000k_1blocks, 'gx--', linewidth=2, basex=2, basey=10, label='v4', markersize=5)



ax.set_ylim(1e0, 1e3)
ax.set_xlim(16, 2048)

#ax.set_xticks(x, minor=True)
#ax.set_yticks(y, minor=True)
#ax.set_xticklabels(('', '32','64','128','256','512','1024', '2048', ''))
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
for i in range(0,len(x4)):
    ax.annotate(s40k_100k_160blocks[i], xy=(x4[i], s40k_100k_160blocks[i]),  xycoords='data',
                xytext=(2, 2), textcoords='offset points', fontsize=datalabelfontsize, color='g')
                #,arrowprops=dict(arrowstyle="->",connectionstyle="arc3,rad=.2")
                #)
for i in range(0,len(x4)):
    ax.annotate(s40k_200k_80blocks[i], xy=(x4[i], s40k_200k_80blocks[i]),  xycoords='data',
                xytext=(2, 2), textcoords='offset points', fontsize=datalabelfontsize, color='m')
                #,arrowprops=dict(arrowstyle="->",connectionstyle="arc3,rad=.2")
                #)
                                
#for i in range(0,len(x4)):
    #ax.annotate(s120k_400k_120blocks[i], xy=(x4[i], s120k_400k_120blocks[i]),  xycoords='data',
                #xytext=(2, 2), textcoords='offset points', fontsize=datalabelfontsize, color='b')
                ###,arrowprops=dict(arrowstyle="->",connectionstyle="arc3,rad=.2")
                ###)                
for i in range(0,len(x4)):
    ax.annotate(s40k_400k_40blocks[i], xy=(x4[i], s40k_400k_40blocks[i]),  xycoords='data',
                xytext=(2, 2), textcoords='offset points', fontsize=datalabelfontsize, color='r')
                ##,arrowprops=dict(arrowstyle="->",connectionstyle="arc3,rad=.2")
                ##)    
for i in range(0,len(x4)):
    ax.annotate(s40k_800k_20blocks[i], xy=(x4[i], s40k_800k_20blocks[i]),  xycoords='data',
                xytext=(2, 2), textcoords='offset points', fontsize=datalabelfontsize, color='b')
                #xytext=(-24, -13), textcoords='offset points', fontsize=datalabelfontsize, color='b')
                #,arrowprops=dict(arrowstyle="->",connectionstyle="arc3,rad=.2")
                #)      
            
                
ax.legend( (
            ##'40k sequences in 80 blocks (hicomb)',             
            ##'40k sequences in 40 blocks (hicomb)',             
            ##'20k sequences in 40 blocks (hicomb)',             
            ##'6k  sequences in 12 blocks (hicomb)',  
            ##'40k sequences in 80 blocks',             
            ##'40k sequences in 40 blocks',             
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
            '40k sequences in 160 blocks',
            #'120k sequences in 120 blocks',
            '40k sequences in 80 blocks',
            #'120k sequences in 40 blocks',
            '40k sequences in 40 blocks',
            '40k sequences in 20 blocks',
            '40k sequences in 1 block',
            ), loc=3)
            
leg = plt.gca().get_legend()
ltext  = leg.get_texts()
plt.setp(ltext, fontsize='14')    # the legend text fontsize

pylab.savefig('mrblast_plot_hicomb_ppt_pairwiseversion.png', dpi=(300))                
pylab.show()
