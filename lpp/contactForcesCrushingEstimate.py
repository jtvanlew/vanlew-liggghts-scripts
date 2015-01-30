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
        xs = np.zeros((npoints, 1))
        ys = np.zeros((npoints, 1))
        zs = np.zeros((npoints, 1))
        print 'Timestep:',str(timestep),'npoints=',str(npoints),'ncells=',str(ncells),'nperiodic=',nperiodic
        


        # create a weibull distribution of 'strength' for pebbles. Use the MATLAB file to generate
        # the average and minimum of F for the pebble in ensemble
        
        # E, Rp
        # alpha = 9.8120
        # beta = 0.8284

        # E, <Rp>
        #alpha = 9.7997
        #beta = 0.8263

        # <E>, Rp
        alpha = 9.5855
        beta = 0.6491

        # <E>, <Rp>
        #alpha = 9.8050
        #beta = 0.6347

        Fs = np.random.gamma(alpha, beta, npoints)


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
                if forcetemp > Fs[counter]:
                    x = 0.5*( forcedata.snaps[fileindex].atoms[index,forcedata.names["x1"]] +
                              forcedata.snaps[fileindex].atoms[index,forcedata.names["x2"]] )
                    y = 0.5*( forcedata.snaps[fileindex].atoms[index,forcedata.names["y1"]] +
                              forcedata.snaps[fileindex].atoms[index,forcedata.names["y2"]] )
                    z = 0.5*( forcedata.snaps[fileindex].atoms[index,forcedata.names["z1"]] +
                              forcedata.snaps[fileindex].atoms[index,forcedata.names["z2"]] )
                else:
                    x = 0
                    y = 0
                    z = 0
            else:
                index = id2.index(id)
                forcetemp = np.sqrt( forcedata.snaps[fileindex].atoms[index,forcedata.names["fx"]]**2 + \
                             forcedata.snaps[fileindex].atoms[index,forcedata.names["fy"]]**2 + \
                             forcedata.snaps[fileindex].atoms[index,forcedata.names["fz"]]**2 )
                if forcetemp > Fs[counter]:
                    x = 0.5*( forcedata.snaps[fileindex].atoms[index,forcedata.names["x1"]] +
                              forcedata.snaps[fileindex].atoms[index,forcedata.names["x2"]] )
                    y = 0.5*( forcedata.snaps[fileindex].atoms[index,forcedata.names["y1"]] +
                              forcedata.snaps[fileindex].atoms[index,forcedata.names["y2"]] )
                    z = 0.5*( forcedata.snaps[fileindex].atoms[index,forcedata.names["z1"]] +
                              forcedata.snaps[fileindex].atoms[index,forcedata.names["z2"]] )
                else:
                    x = 0
                    y = 0
                    z = 0

            forceArray[counter] = forcetemp

            xs[counter] = x
            ys[counter] = y
            zs[counter] = z
            counter += 1
        xs = xs[np.all(xs != 0, axis=1)]
        ys = ys[np.all(ys != 0, axis=1)]
        zs = zs[np.all(zs != 0, axis=1)]

        print 'percent of broken pebbles = ' +str(len(xs)/8000.*100.)+' %'
        






        plt.hold(True)
        # definitions for the axes
        left, width = 0.13, 0.63
        bottom, height = 0.1, 0.65
        left_h = left+width+0.02
        bottom_h = bottom+height+0.02

        rect_scatter = [left, bottom, width, height]
        rect_histx = [left, bottom_h, width, 0.2]
        rect_histy = [left_h, bottom, 0.2, height]

        # start with a rectangular Figure
        plt.figure(2, figsize=(8,8))

        axScatter = plt.axes(rect_scatter)
        axHistx = plt.axes(rect_histx)
        axHisty = plt.axes(rect_histy)

        # no labels
        axHistx.xaxis.set_major_formatter(nullfmt)
        axHisty.yaxis.set_major_formatter(nullfmt)
        
        # the scatter plot:
        markerArea = np.pi*5**3
        markerColor = 'c'

        axScatter.scatter(xs, ys, s = markerArea, c = markerColor, alpha = 0.7)

        # now determine nice limits by hand:
        bins = 20
        xlim = 0.005
        ylim = 0.00375
        axScatter.set_xlim( (-xlim, xlim) )
        axScatter.set_ylim( (-ylim, ylim) )
        axScatter.set_xlabel( 'x (m)' )
        axScatter.set_ylabel( 'y (m)' )
        
        axHistx.hist(xs, bins)
        axHisty.hist(ys, bins, orientation='horizontal')

        axHistx.set_xlim( axScatter.get_xlim() )
        axHisty.set_ylim( axScatter.get_ylim() )
        locs, labels = plt.xticks()
        axHisty.set_xticklabels('', rotation=45)
        
        
        epsFile = 'FgtFc_scatter_hist.eps'
        epsFile = os.path.join(outputdir, epsFile)
        plt.savefig(epsFile)
        pngFile = 'FgtFc_scatter_hist.png'
        pngFile = os.path.join(outputdir, pngFile)
        plt.savefig(pngFile)



        plt.figure(3)
        plt.hist(zs, bins)
        plt.xlabel('z coordinate (m)')
        plt.ylabel('count')

        epsFile = 'FgtFc_z_hist.eps'
        epsFile = os.path.join(outputdir, epsFile)
        plt.savefig(epsFile)
        pngFile = 'FgtFc_z_hist.png'
        pngFile = os.path.join(outputdir, pngFile)
        plt.savefig(pngFile)





        fAve = np.append(fAve, np.mean(forceArray))
        n, bins, patches = plt.hist(forceArray, 100, normed=True, histtype = 'step')
        plt.close()
        plt.figure(figureNumber)
        x = np.linspace(0,max(forceArray),100)
        plt.plot(x,n,'o')
        plt.xlabel('Force (N)')
        plt.ylabel('Probability')
        
        
        epsFile = 'forceProfiles_'+str(timestep)+'.eps'
        epsFile = os.path.join(outputdir, epsFile)
        plt.savefig(epsFile)
        pngFile = 'forceProfiles_'+str(timestep)+'.png'
        pngFile = os.path.join(outputdir, pngFile)
        plt.savefig(pngFile)
        #figureNumber += 1
        data = np.zeros([2,len(x)])
        data[0,:] = x
        data[1,:] = n
        scipy.io.savemat('data.mat',mdict={'data':data})

        fileindex += 1
        timestep = forcedata.next()
    # end of main loop - close group file

print '\n Finished! \n '
print '\n see dump directory for images\n '
#plt.show()
#np.savetxt('averageForces.txt', fAve)
#stressarray = np.loadtxt(open('stress.csv'), delimiter = ',', skiprows=1)
#stress = stressarray[0]/stressarray[2]
#fig, ax = plt.subplots(1)
#ax.plot(stress, fAve,label='Average contact force')
#plt.ylabel('Force (N)')
#plt.xlabel('External load (MPa)')
#ax.plot(stress,fc,label='Crush force')
#ax.fill_between(stress, fAve, fc, facecolor = 'yellow', alpha = 0.5)
#plt.axis((0,stress[-1],0,300))
#ax.grid()
#plt.legend()
#plt.show()