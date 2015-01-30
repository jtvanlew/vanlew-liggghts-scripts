# Tool for pulling coordination number/atom out of dump data and performing
# an averaging at each timestep. This saves a file that can be loaded by 
# MATLAB for plotting. Future version should plot directly from Python.
# copyright Jon Van Lew 2013
import numpy as np
import scipy.io
import sys, os
from dump import dump
import matplotlib.pyplot as plt


# Load the first case to find its size and initialize variables
filename = str(sys.argv[1])
initializeData = dump(filename,0)
timestep = initializeData.next()
pebbleIDs = initializeData.vecs(timestep,"id")

# Initialize variables. they are N x M
# where N is the number of pebbles (actually pebbles+1 to allow appending the snapshot onto the array) 
# and M is the number of snapshots - 1

coord = np.zeros((np.size(pebbleIDs)+1,len(sys.argv)-1))
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
    pebbleIDs = initializeData.vecs(timestep,"id")
    
    coord = np.zeros((np.size(pebbleIDs)+1,len(sys.argv)-1))
    temperature = np.zeros((np.size(pebbleIDs)+1,len(sys.argv)-1))
    totalSnaps = len(sys.argv)-1
    x = np.zeros((np.size(pebbleIDs)+1,len(sys.argv)-1))
    y = np.zeros((np.size(pebbleIDs)+1,len(sys.argv)-1))
    z = np.zeros((np.size(pebbleIDs)+1,len(sys.argv)-1))
    r = np.zeros((np.size(pebbleIDs)+1,len(sys.argv)-1))
    # Here we have a single snapshot loaded, and we're going to pull out the pebble vectors' of data
    while timestep >=0:
        index,time,flag = forcedata.iterator(0)
        coord[0:, timeSnap] = np.append(timestep, forcedata.vecs(time,"c_coord[0]"))
        temperature[0:, timeSnap] = np.append(timestep, forcedata.vecs(time,"f_Temp[0]"))
        r[0:, timeSnap] = np.append(timestep, forcedata.vecs(time,"radius"))
        x[0:, timeSnap] = np.append(timestep, forcedata.vecs(time,"x"))
        y[0:, timeSnap] = np.append(timestep, forcedata.vecs(time,"y"))
        z[0:, timeSnap] = np.append(timestep, forcedata.vecs(time,"z"))
        print 'Processing timestep:',str(timestep), ', Step number',str(timeSnap+1),'/',str(totalSnaps)
        timestep = forcedata.next()
        timeSnap += 1



# the snapshots could be loaded in random orders (or just non-sequentially), so here we rearrange all of them
# in order of the snapshot # which is the first entry of each array
r = r.T[r.T[:,0].argsort(),].T

x = x.T[x.T[:,0].argsort(),].T
y = y.T[y.T[:,0].argsort(),].T
z = z.T[z.T[:,0].argsort(),].T
temperature = temperature.T[temperature.T[:,0].argsort(),].T
coord = coord.T[coord.T[:,0].argsort(),].T



print "Now binning the snapshot", str(x[0,-1]), "to get temperature profiles"
# The variables had the timestep appended onto the first row. This helps for sorting.
# But here we strip that off so we can do our histogram data. We sort along x and then use the same
# sorting to rearrange Temp (like an excel table, for example)


finalData = np.array([x[1:, -1], temperature[1:, -1], coord[1:, -1], r[1:, -1]])
finalData = finalData.T[finalData.T[:,0].argsort(),].T
finalX = finalData[0,:]
finalTemperature = finalData[1,:]
finalCoord = finalData[2,:]
finalRadius = finalData[3,:]
pebbleVolume=0.
for radius in finalRadius:
    pebbleVolume += 4./3. * np.pi * radius**3

totalBins = 50
hist, bins = np.histogram(finalX, bins = totalBins)

count = 0
tempBinMean = np.zeros(totalBins)
xBinMean = np.zeros(totalBins)
coordBinMean = np.zeros(totalBins)
for numInBin in hist:
    
    tempBinMean[count] = np.average(finalTemperature[0:numInBin])
    xBinMean[count] = np.average(finalX[0:numInBin])
    coordBinMean[count] = np.average(finalCoord[0:numInBin])
    
    count += 1
    finalTemperature = np.delete(finalTemperature, np.s_[0:numInBin], None)
    finalX = np.delete(finalX, np.s_[0:numInBin], None)
    finalCoord = np.delete(finalCoord, np.s_[0:numInBin], None)
# Make dimensionless x and temperature
xstar = xBinMean/np.max(finalData[0,:])
tempMax,tempMin = np.max(tempBinMean),np.min(tempBinMean)
deltaT = tempMax - 573
theta = (tempBinMean - tempMin)/deltaT
#scipy.io.savemat("MATLAB/temperature_"+str(timestep), mdict={'tempProfile':tempBinMean})
thetaPredicted = 1-xstar**2

# Find the height of the pebble bed (taken from the highest 50 pebbles)
finalZ = z[1:, -1]
H = np.average(finalZ.T[finalZ.T.argsort(),].T[-50:])


# Find the effective thermal conductivity
q_nuc = 8.e6

L = 15*0.0005
W = .01
N = np.size(pebbleIDs)
V = 4* L * W * H
phi = N * 4/3. * np.pi * 0.0005**3 / ( V )
print "packing fraction: "+str(phi)
print "nuclear power deposited: " + str(phi*q_nuc*V) + " W"
#q_h = 0.9*0.8025697597

print "temperature delta (centerline to wall): "+str(deltaT) + " K"

kcond = q_nuc*phi*L**2/(2*deltaT) 
#kenth = q_h/(H*deltaT)
k = kcond #- kenth

#print "conductivity without correction: "+str(kcond)
#print "'conductivity' of enthalpy: "+str(kenth)
print "effective conductivity: "+str(k)+" W/m-K"


# Plot results along with predicted curve
plt.figure(1)
#plt.plot(xstar, theta, 'bo-', xstar, thetaPredicted, 'rx-')
plt.plot(xstar, tempBinMean, 'bo-')
plt.ylabel('Temperature (K)')
plt.xlabel(r'$x^* = x/L$, Dimensionless Length')
#plt.axis([-1, 1, 0, 1])
#plt.legend((r'DEM, $k_{eff}$ = %.3f'%k, 'Predicted 1D'), loc = 'lower center')
#plt.legend(('DEM simulation', 'Predicted 1D'), loc = 'lower center')
epsFile = 'tempProfile.eps'
epsFile = os.path.join(outputdir, epsFile)
plt.savefig(epsFile)

plt.figure(2)
plt.plot(xstar, coordBinMean, 'bo-')
plt.ylabel(r'$<Z>$, Average coordination number of slice')
plt.xlabel(r'$x^* = x/L$, Dimensionless Length')
plt.axis([-1, 1, 1, 12])
plt.legend((r'DEM, $k_{eff}$ = %.3f'%k, 'Predicted 1D'), loc = 'lower center')

epsFile = 'coordProfile.eps'
epsFile = os.path.join(outputdir, epsFile)
plt.savefig(epsFile)

plt.show()





