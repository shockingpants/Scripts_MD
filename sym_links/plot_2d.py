#!/opt/local/bin/python2.7
from strtup import *
import numpy
import sys
import matplotlib.pyplot as plt

target=sys.argv[1:]
print type(target[0])
print target[0]
data = numpy.genfromtxt(target[0],skip_header=1)
print data
print type(data)
plt.plot(data[:,0], data[:,1])
plt.show()
