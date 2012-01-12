import pylab
from numpy import *
import string
from collections import defaultdict

datafilename = "m1000_1024-alldb.out.txt"
#datafilename = "test.dat"
datafile = open(datafilename,"r")
count = 0
acc = []
v0 = 0
for line in datafile.readlines():
    parts= line.split(',')
    v1 = parts[0]
    acc.append(int(v1))
datafile.close()

print 'len(acc)=', len(acc)

print '\nget freq\n'
d = defaultdict(int)
for word in acc:
    d[word] += 1

f=open("dict.txt",'w')
for i in d.keys():
    f.write(str(i))
    f.write(' %s'%(d[i]))
    f.write('\n')
f.close()

items = [(v, k) for k, v in d.items()]
items.sort()
print 'len(items)=', len(items)


freq = []
for i in range(0,len(items)):
    freq.append(items[i][0])

d2 = defaultdict(int)
for word in freq:
    d2[word] += 1

f=open("hitcounttable.txt",'w')
for key in sorted(d2.iterkeys()):
    f.write("%s %s\n" % (key, d2[key]))
f.close()


items2 = [(k, v) for k, v in d2.items()]

hit = []
for i in range(0,len(items2)):
    hit.append(items2[i][1])

################################################################################
pylab.figure(1)
ax = pylab.subplot(111)

#bins = range(1, 10)+range(10,100,10)+range(100,1001,100)
#ax.hist(hit, bins=bins)
ax.hist(hit, bins=500)

#ax.set_ylim(1e2, 1e4)
ax.set_xlim(0, 1000)
#ax.yaxis.grid(True, linestyle='-.', which='minor')
#ax.xaxis.grid(True, linestyle='-.', which='minor')

fontsize=16
for tick in ax.xaxis.get_major_ticks():
    tick.label1.set_fontsize(fontsize)
fontsize=16
for tick in ax.yaxis.get_major_ticks():
    tick.label1.set_fontsize(fontsize)
    
ax.set_xlabel('Hit count per query sequence', fontsize=20)
ax.set_ylabel('Number of query sequences per hit', fontsize=20)

pylab.savefig('mrblast_hist.png', dpi=300)            
pylab.show()

 
