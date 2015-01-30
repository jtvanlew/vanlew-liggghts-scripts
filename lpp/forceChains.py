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


from dump import dump
import scipy.io
import numpy as np
import sys, os
import matplotlib.pyplot as plt
from matplotlib.ticker import NullFormatter
nullfmt   = NullFormatter()         # no labels


fAve = np.array([0])
figureNumber = 1
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
    outputdir = os.path.join(inputdir,fileprefix+'_'+os.path.split(inputdir)[-1])
    try:
        os.mkdir(outputdir)
    except:
        pass
    
    # Read in the dump file - since we can have many contacts (i.e. >> nparticles)
    # and many timesteps I will deal with one timestep at a time in memory,
    # write to the appropriate .vtu file for a single timestep, then move on.
    
    forcedata = dump(filename,0)
    
    fileindex = 0
    timestep = forcedata.next()
    # loop through available timesteps

    while timestep >= 0:
    
        # default data are stored as pos1 (3) pos2 (3) id1 id2 periodic_flag force (3) -> 12 columns
        #
        # if contactArea is enabled, that's one more (13) and heatflux (14)
        #
        # assign names to atom columns (1-N)
        forcedata.map(1,"x1",2,"y1",3,"z1",4,"x2",5,"y2",6,"z2",7,"id1",8,"id2",9,"periodic",10,"fx",11,"fy",12,"fz")
        
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
        forceArray = np.zeros((npoints,1))
        saveID1s = np.zeros((npoints, 1))
        saveID2s = np.zeros((npoints, 1))
        saveForce1s = np.zeros((npoints, 1))
        saveForce2s = np.zeros((npoints, 1))
        print 'Timestep:',str(timestep),'npoints=',str(npoints),'ncells=',str(ncells),'nperiodic=',nperiodic
        


        Fs = 0




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
                if forcetemp > Fs:
                    saveID1 = forcedata.snaps[fileindex].atoms[index,forcedata.names["id1"]]
                    saveID2 = forcedata.snaps[fileindex].atoms[index,forcedata.names["id2"]]
                    saveForce = forcetemp
                else:
                    saveID1 = 0
                    saveID2 = 0
                    saveForce = 0
            else:
                index = id2.index(id)
                forcetemp = np.sqrt( forcedata.snaps[fileindex].atoms[index,forcedata.names["fx"]]**2 + \
                             forcedata.snaps[fileindex].atoms[index,forcedata.names["fy"]]**2 + \
                             forcedata.snaps[fileindex].atoms[index,forcedata.names["fz"]]**2 )
                if forcetemp > Fs:
                    saveID1 = forcedata.snaps[fileindex].atoms[index,forcedata.names["id1"]]
                    saveID2 = forcedata.snaps[fileindex].atoms[index,forcedata.names["id2"]]
                    saveForce = forcetemp
                else:
                    saveID1 = 0
                    saveID2 = 0
                    saveForce = 0

            forceArray[counter] = forcetemp

            saveID1s[counter] = saveID1
            saveID2s[counter] = saveID2
            saveForce1s[counter] = saveForce
            saveForce2s[counter] = saveForce

            counter += 1


        print np.mean(forceArray)

        saveID1s = saveID1s[np.all(saveID1s != 0, axis=1)]
        saveID2s = saveID2s[np.all(saveID2s != 0, axis=1)]
        saveForce1s = saveForce1s[np.all(saveForce1s != 0, axis=1)]
        saveForce2s = saveForce2s[np.all(saveForce2s != 0, axis=1)]

        # concatenate into a single list
        saveIDs = []
        saveIDs = saveID1s[:,0]
        saveIDs = list(saveIDs)
        saveIDs.extend(saveID2s[:,0])

        saveForces = []
        saveForces = saveForce1s[:,0]
        saveForces = list(saveForces)
        saveForces.extend(saveForce2s[:,0])

        idf = zip(saveIDs, saveForces)

        sortedidf = sorted(idf)

        saveIDs = [ids[0] for ids in sortedidf]
        saveForces = [forces[1] for forces in sortedidf]
        saveForces = saveForces

        print 'percent of pebbles in this chain network = ' +str(len(saveIDs)/8000.*100.)+' %'
        
        scipy.io.savemat('IDs.mat',mdict={'ids':saveIDs})
        scipy.io.savemat('Forces.mat',mdict={'force':saveForces})

        fileindex += 1
        timestep = forcedata.next()
    # end of main loop - close group file

print '\n Finished! \n '
print '\n list of ids in IDs.mat\n '
print '\n list of forces in Forces.mat\n '
