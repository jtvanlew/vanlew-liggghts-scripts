# Tool for pulling coordination number/atom out of dump data and performing
# an averaging at each timestep. This saves a file that can be loaded by 
# MATLAB for plotting. Future version should plot directly from Python.
# copyright Jon Van Lew 2013
import numpy as np
import scipy.io
import sys, os
from dump import dump



# Load the first case to find its size and initialize variables
filename = str(sys.argv[1])
forcedata = dump(filename,0)
timestep = forcedata.next()
pebbleIDs = forcedata.vecs(timestep,"id")


index,time,flag = forcedata.iterator(0)

r = forcedata.vecs(time,"radius")
x = forcedata.vecs(time,"x")
y = forcedata.vecs(time,"y")
z = forcedata.vecs(time,"z")

xyzr = x, y, z, r
xyzr = zip(*xyzr)

import csv

data = open('liggghts2comsol.csv','wb')
write = csv.writer(data)

for row in xyzr:
    write.writerow(row)
data.close()


