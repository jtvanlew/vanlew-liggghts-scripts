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
initializeData = dump(filename,0)
timestep = initializeData.next()
pebbleIDs = initializeData.vecs(timestep,"id")

# Initialize variables. they are N x M
# where N is the number of pebbles (actually pebbles+1 to allow appending the snapshot onto the array) 
# and M is the number of snapshots - 1

#coord = np.zeros((np.size(pebbleIDs)+1,len(sys.argv)-1))
temperature = np.zeros((np.size(pebbleIDs)+1,len(sys.argv)-1))
totalSnaps = len(sys.argv)-1
x = np.zeros((np.size(pebbleIDs)+1,len(sys.argv)-1))
y = np.zeros((np.size(pebbleIDs)+1,len(sys.argv)-1))
z = np.zeros((np.size(pebbleIDs)+1,len(sys.argv)-1))


# Loop over all the dump_*.liggghts filenames
timeSnap = 0
for snapshot in sys.argv[1:]:

    filename = str(snapshot)
    inputpath = os.path.abspath(filename)
    inputdir = os.path.split(inputpath)[0] # what does this do?
    
    outputdir = os.path.join(inputdir, 'MATLAB')
    try:
        os.mkdir(outputdir)
    except:
        pass
    forcedata = dump(filename,0)
    timestep = forcedata.next()
    print timestep
    # Here we have a single snapshot loaded, and we're going to pull out the pebble vectors' of data
    while timestep >=0:
        index,time,flag = forcedata.iterator(0)
        #coord[0:, timeSnap] = np.append(timestep, forcedata.vecs(time,"c_coord[0]"))
        temperature[0:, timeSnap] = np.append(timestep, forcedata.vecs(time,"f_Temp[0]"))
        x[0:, timeSnap] = np.append(timestep, forcedata.vecs(time,"x"))
        y[0:, timeSnap] = np.append(timestep, forcedata.vecs(time,"y"))
        z[0:, timeSnap] = np.append(timestep, forcedata.vecs(time,"z"))
        print 'Processing timestep:',str(timestep), ', Step number',str(timeSnap+1),'/',str(totalSnaps)
        timestep = forcedata.next()
        timeSnap += 1



# the snapshots could be loaded in random orders (or just non-sequentially), so here we rearrange all of them
# in order of the snapshot # which is the first entry of each array

x = x.T[x.T[:,0].argsort(),].T
y = y.T[y.T[:,0].argsort(),].T
z = z.T[z.T[:,0].argsort(),].T
temperature = temperature.T[temperature.T[:,0].argsort(),].T
#coord = coord.T[coord.T[:,0].argsort(),].T



print "Now we'll save all the data in the MATLAB subdirectory for any future manipulation"
#scipy.io.savemat('MATLAB/coordinationNumber.mat',mdict={'coordNum':coord})
scipy.io.savemat('MATLAB/temperature.mat',mdict={'temperature':temperature})
scipy.io.savemat('MATLAB/x.mat',mdict={'x':x})
scipy.io.savemat('MATLAB/y.mat',mdict={'y':y})
scipy.io.savemat('MATLAB/z.mat',mdict={'z':z})


