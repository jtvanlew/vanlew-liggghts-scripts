""" Filling a box with walls on the sides (x-direction) and periodic y. The top is 'open' unless the user
wants to also have a zwall.

This simulation is an example of using the vl_dem module to create a clean input script file.

We define properties in the 'USER INPUT' section. below that, you shouldn't need to edit anything.

This script will fill a volume to the specified void fraction, run a relaxation (because we allow overlap)
when filling to the void fraction, then allows the ensemble to settle based on the criteria of minimum
necessary kinetic energy. Once the system's KE passes below the threshold defined by the user, a restart
file is saved.

Directories are created automatically for post data and restart files.

Copyright 2015 Jon Van Lew (or, whatever)

"""

import      mpi4py
import      sys, os
from        lammps      import  lammps
import      numpy       as      np
from        vl_dem  	import  *
# invoke LAMMPS session
lmp = lammps()




#-----------------------------------------------------------------------------------------------------------
# USER INPUT

# dt (s) and times for outputs or checks in simulation
dt       		= 5.e-8		# s
dump_time      		= 0.01 		# s
CTE_check_time 		= 0.001 	# s
screen_print_time 	= 0.001 	# s



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Mechanical 
E       = 100e9       # Pa
nu      = 0.24
rho     = 3440        # kg/m3
Rp      = 0.0005      # m
dp      = Rp*2
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# choose only one of the following and use the appropriate insertion scheme
# N       = 10 		# desired number of pebbles in the system
phi 	= 0.62		# desired packing fraction
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# If doing a gaussian distribution on radius, define the average and standard deviation
mu      = Rp
sigma   = Rp/15.
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# geometry limits (these are per side for x and y, and total extent in z)
xlim    = dp*10
ylim    = dp*5
zlim    = dp*25.
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Coefficients of friction and restitution
mu    = 0.2
gamma = 0.1
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Heat transfer
k     = 2.4         # W/m-K
Cp    = 1.          # kJ/kg-K
Ti    = 300         # K
beta  = 15.e-6
Q     = 8.e6
Qp    = Q * (4./3 * 3.1415 * Rp**3)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# some programming variables that can be edited
# manually specifiy processor decomposition of space if desired
processor_layout = '4 2 1'

min_energy  = 1.e-7
#-----------------------------------------------------------------------------------------------------------



#-----------------------------------------------------------------------------------------------------------
# pile up geometric values into this list for putting into functions
geometry = [xlim, ylim, zlim, Rp]

# create the directories for saving files
output_dir  = 'post'
restart_dir = 'restarts'
post_dir    = output_dir+'/filling'
mesh_dir    = post_dir+'/meshes'

make_directory( output_dir, os  )
make_directory( restart_dir, os )
make_directory( post_dir, os    )
make_directory( mesh_dir, os    )
#-----------------------------------------------------------------------------------------------------------




#-----------------------------------------------------------------------------------------------------------
# assign the dt and initialize the lammps commands for the filling script

CTE_check_steps = int(CTE_check_time/dt)
dump_steps     	= int(dump_time/dt)
print_steps    	= int(screen_print_time/dt)

define_dt( dt, lmp)

initialize_filling_liggghts( geometry, lmp)
#-----------------------------------------------------------------------------------------------------------




#-----------------------------------------------------------------------------------------------------------
# Specify material properties
dp = define_properties()
dp.youngs_modulus(E, lmp)
dp.poissons_ratio(nu, lmp)
dp.coefficient_restitution(gamma, lmp)
dp.coefficient_friction(mu, lmp)
dp.thermal_conductivity(k, lmp)
dp.specific_heat(Cp, lmp)
dp.initial_temperature(Ti, lmp)
dp.characteristic_velocity(lmp)

# Initialize the thermal expansion code
thermal_expansion(Rp, Ti, beta, CTE_check_steps, lmp)
#-----------------------------------------------------------------------------------------------------------




#-----------------------------------------------------------------------------------------------------------
# Specify system geometries and create templates for the pebbles
create_horizontal_walls(xlim, lmp)
define_phi_pebbles(geometry, dump_steps, rho, phi, lmp)
#-----------------------------------------------------------------------------------------------------------



#-----------------------------------------------------------------------------------------------------------
# specify how frequently to dump and where to put the files. change screen display to normal custom version
custom_screen_output(print_steps, lmp)
set_dumps(dump_steps, post_dir, lmp)
#-----------------------------------------------------------------------------------------------------------




#-----------------------------------------------------------------------------------------------------------
# the pebbles are inserted with overlap allowed so to quickly fill to the desired void fraction. Now we 
# limit their travel in each integration step and relax toward an equilibirum packing
relax_insertion(dump_steps, Rp, lmp)
#-----------------------------------------------------------------------------------------------------------




#-----------------------------------------------------------------------------------------------------------
# the system is allowed to settle under normal verlet integration until the kinetic energy of the ensemble
# is less than a specified value. 

run_to_kinetic_steady_state(dump_steps, min_energy, lmp)
#-----------------------------------------------------------------------------------------------------------




#-----------------------------------------------------------------------------------------------------------
# save to a restart file for heating or crushing or whatever
lmp.command('write_restart '+restart_dir+'/filled_*.restart')
#-----------------------------------------------------------------------------------------------------------

