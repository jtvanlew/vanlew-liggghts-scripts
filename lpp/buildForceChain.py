from dump import dump
import scipy.io as sio
import numpy as np
import sys, os



ids = sio.loadmat('IDs.mat')['ids']

forces = sio.loadmat('Forces.mat')['force']

dumpfile = sys.argv[1]
d = dump(dumpfile)  

d.sort("id") 
# d.map(1,"id",2,"type",3,"type",4 "x",5,"y",6,"z",7,"ix",8,"iy",9,"iz",\
# 	  10,"vx",11,"vy",12,"vz",13,"fx",14,"fy",15,"fz",\
# 	  16,"omegax",17,"omegay",18,"omegaz",19,"radius",20,"c_coord[0]")

index,timestep,flag = d.iterator(0)


compareString = ''
idset = set()
count = 0
for ID in ids:
	ID = int(ID)
	if ID in idset:
		forces[count] = 0

	compareString+='$id == '+str(int(ID))+' or '
	idset.add(ID)

	count += 1


forces = forces[np.all(forces != 0, axis=1)]


compareString = compareString[:-4]
d.aselect.all()
d.aselect.test(compareString)
d.setv("contactForce",forces)  
d.write("dump_forcechain.liggghts")	 
