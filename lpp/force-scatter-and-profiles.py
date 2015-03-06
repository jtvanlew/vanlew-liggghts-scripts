#!/usr/bin/python
"""
A simple routine to load in a LIGGGHTS hybrid dump file containing
contact and contact force data and convert into a .vtk unstructured
grid which can be used to visualise the force network.

evtk is used to write binary VTK files:
https://bitbucket.org/pauloh/pyevtk

The pizza.py bdump command is used to handle LIGGGHTS dump files and
therefore PYTHONPATH must include the pizza/src location.

NOTE: bdump is NOT included in granular pizza, and should be taken
from the standard LAMMPS pizza package!

NOTE: it is impossible to tell from the bdump header which values
have been requested in the compute, so check that your compute
and dump match the format here - this will be checked in future!

"""

from evtk.vtk import VtkFile, VtkGroup, VtkUnstructuredGrid
from dump import dump
import scipy.io
import numpy as np
import sys, os
import matplotlib.pyplot as plt
import numpy as np

# TODO: use a try/except here to check for missing modules, and fallback to ASCII VTK if evtk not found
# TODO: ask for timestep or timestep range as input (code is NOT efficient and large files = long runtimes!)
# TODO: write celldata for contact area and heat flux (if present)

# Check for command line arguments
#if len(sys.argv) != 2:
#        sys.exit('Usage: dump2forcenetwork.py <filename>, where filename is typically dump.<runname>')
        
#elif len(sys.argv) == 2: # we have one input param, that should be parsed as a filename
#    filename = str(sys.argv[1])
#    if not os.path.isfile(filename):
#        sys.exit('File ' + filename + ' does not exist!')

# The check above wouldn't work for me, so I just assume the load is properly executed...

fig, ax = plt.subplots(1)
alpha = 0.2
plt.xlabel('Dimensionless x')
plt.ylabel(r'$\langle F \rangle ^{1/3}$ ($N^{1/3}$)')
plt.xlim([-1, 1])







timeSnap = 0
for snapshot in sys.argv[1:]:
    filename = str(snapshot)

    splitname = filename.split('.')

    if len(splitname) == 2 and splitname[0].lower() == 'dump':
        fileprefix = splitname[1]
    else:
      fileprefix = splitname[0]

    inputpath = os.path.abspath(filename)
    inputdir = os.path.split(inputpath)[0]

    # create a sub-directory for the output .vtu files
    outputdir = os.path.join(inputdir,fileprefix)
    try:
        os.mkdir(outputdir)
    except:
        pass
    
    # Read in the dump file - since we can have many contacts (i.e. >> nparticles)
    # and many timesteps I will deal with one timestep at a time in memory,
    # write to the appropriate .vtu file for a single timestep, then move on.
    
    forcedata = dump(filename,0)
    
    groupfile = fileprefix
    groupfile = os.path.join(inputdir,groupfile)
    groupfile = VtkGroup(groupfile)
    
    fileindex = 0
    timestep = forcedata.next()

    # check that we have the right number of colums (>11)
    #
    # NOTE: the first timesteps are often blank, and then natoms returns 0, so this doesn't really work...
    #
    #if forcedata.snaps[fileindex].natoms !=0 and len(forcedata.snaps[0].atoms[0]) < 11:
    #    print "Error - dump file requires at least all parameters from a compute pair/gran/local id pos force (12 in total)"
    #    sys.exit()
    
    # loop through available timesteps
    
    while timestep >= 0:
    
        # default data are stored as pos1 (3) pos2 (3) id1 id2 periodic_flag force (3) -> 12 columns
        #
        # if contactArea is enabled, that's one more (13) and heatflux (14)
        #
        # assign names to atom columns (1-N)
        forcedata.map(1,"x1",2,"y1",3,"z1",4,"x2",5,"y2",6,"z2",7,"id1",8,"id2",9,"periodic",10,"fx",11,"fy",12,"fz",13,"contactArea",14,"heatflux")

        # check for contact data (some timesteps may have no particles in contact)
        #
        # NB. if one loads two datasets into ParaView with defined timesteps, but in which
        # one datasets has some missing, data for the previous timestep are still displayed - 
        # this means that it is better here to generate "empty" files for these timesteps.
    

        # number of cells = number of interactions (i.e. entries in the dump file)
        ncells = len(forcedata.snaps[fileindex].atoms)

        # number of periodic interactions
        periodic = np.array(forcedata.snaps[fileindex].atoms[:,forcedata.names["periodic"]],dtype=bool)
        nperiodic = sum(periodic)

        # number of non-periodic interactions (which will be written out)
        nconnex = ncells - nperiodic

        # extract the IDs as an array of integers
        id1 = np.array(forcedata.snaps[fileindex].atoms[:,forcedata.names["id1"]],dtype=long)
        id2 = np.array(forcedata.snaps[fileindex].atoms[:,forcedata.names["id2"]],dtype=long)

        # and convert to lists
        id1 = id1.tolist()
        id2 = id2.tolist()

        # concatenate into a single list
        ids = []
        ids = id1[:]
        ids.extend(id2)

        # convert to a set and back to remove duplicates, then sort
        ids = list(set(ids))
        ids.sort()

        # number of points = number of unique IDs (particles)
        npoints = len(ids)
        f = np.zeros( npoints, dtype=np.float64 )
        # create empty arrays to hold x,y,z data
        x = np.zeros( npoints, dtype=np.float64)
        y = np.zeros( npoints, dtype=np.float64)
        z = np.zeros( npoints, dtype=np.float64)

        print 'Timestep:',str(timestep),'npoints=',str(npoints),'ncells=',str(ncells),'nperiodic=',nperiodic

        # Point data = location of each unique particle
        #
        # The order of this data is important since we use the position of each particle
        # in this list to reference particle connectivity! We will use the order of the 
        # sorted ids array to determine this.

        counter = 0   
        for id in ids:
            if id in id1:
                index = id1.index(id)
                forcetemp = np.sqrt( forcedata.snaps[fileindex].atoms[index,forcedata.names["fx"]]**2 + \
                             forcedata.snaps[fileindex].atoms[index,forcedata.names["fy"]]**2 + \
                             forcedata.snaps[fileindex].atoms[index,forcedata.names["fz"]]**2 )

                xtemp,ytemp,ztemp = forcedata.snaps[fileindex].atoms[index,forcedata.names["x1"]], \
                        forcedata.snaps[fileindex].atoms[index,forcedata.names["y1"]], \
                        forcedata.snaps[fileindex].atoms[index,forcedata.names["z1"]]
            else:
                index = id2.index(id)
                forcetemp = np.sqrt( forcedata.snaps[fileindex].atoms[index,forcedata.names["fx"]]**2 + \
                             forcedata.snaps[fileindex].atoms[index,forcedata.names["fy"]]**2 + \
                             forcedata.snaps[fileindex].atoms[index,forcedata.names["fz"]]**2 )
                xtemp,ytemp,ztemp = forcedata.snaps[fileindex].atoms[index,forcedata.names["x2"]], \
                        forcedata.snaps[fileindex].atoms[index,forcedata.names["y2"]], \
                        forcedata.snaps[fileindex].atoms[index,forcedata.names["z2"]]
            f[counter]=forcetemp**(1/3.)
            x[counter]=xtemp
            y[counter]=ytemp
            z[counter]=ztemp           
            counter += 1
        x = x/0.01
        xspan = np.abs(np.max(x)) + np.abs(np.min(x))
        Nbins = 50
        force_bin = np.zeros([Nbins])
        xbin = np.zeros([Nbins])
        dx = xspan/(Nbins)

        x0 = np.min(x)
        for i in np.arange(0,Nbins):
            x1 = x0 + dx
            if i == 0:
                xbin[i] = x0
            else:
                xbin[i] = x0 + dx
            force_bin[i] = 0
            count = 0
            index = 0
            for xvalue in x:
                if xvalue >= x0 and xvalue < x1:
                    force_bin[i] += f[index]
                    count += 1
                index += 1
            force_bin[i] /= count
            x0 = x1
        force_average = np.mean(force_bin)
        fileindex += 1
        timestep = forcedata.next()

    # end of main loop - close group file

    timeSnap +=1
plt.scatter(x,f,c = 'c', alpha = alpha, label="DEM Data")
plt.ylim([0, np.max(f)])
plt.plot(xbin, force_bin, c = 'k', linewidth = 2, label="Binned average" )
plt.title(r"Average $\langle F \rangle^{1/3}$ = %s $N^{1/3}$"%(np.round(force_average,2)))
plt.legend(loc='best')
epsFile = 'force-profile.eps'
epsFile = os.path.join(outputdir, epsFile)
plt.savefig(epsFile)
pngFile = 'force-profile.png'
pngFile = os.path.join(outputdir, pngFile)
plt.savefig(pngFile)
plt.show()
#scipy.io.savemat('MATLAB/averageContactRadius.mat',mdict={'a':a})
