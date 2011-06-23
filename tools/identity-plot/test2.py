from matplotlib import pyplot, lines
import numpy

x = numpy.linspace(0,10,100)
y = numpy.sin(x)*(1+x)
fig = pyplot.figure()
ax = pyplot.subplot(111)
ax.plot(x,y,label='a')

# new clear axis overlay with 0-1 limits
#ax2 = pyplot.axes([0,0,1,1], axisbg=(1,1,1,0))
ax2 = pyplot.axes([0,0,1,1])

x,y = numpy.array([[0.05, 0.1, 0.9], [0.05, 0.5, 0.9]])
line = lines.Line2D(x, y, lw=5., color='r', alpha=0.4)
ax2.add_line(line)

pyplot.show()
