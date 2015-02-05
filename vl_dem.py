class define_properties:
	def youngs_modulus(self, E, lmp):
		lmp.command("fix m1 all property/global youngsModulus peratomtype " + str(E))
	def poissons_ratio(self, nu, lmp):
		lmp.command("fix m2 all property/global poissonsRatio peratomtype " + str(nu))
	def coefficient_restitution(self, gamma, lmp):		
		lmp.command("fix m3 all property/global coefficientRestitution peratomtypepair 1 " + str(gamma))
	def coefficient_friction(self, mu, lmp):
		lmp.command("fix m4 all property/global coefficientFriction peratomtypepair 1 " + str(mu))
	def thermal_conductivity(self, k, lmp):		
		lmp.command("fix htc1 all property/global thermalConductivity peratomtype " + str(k))
	def specific_heat(self, Cp, lmp):
		lmp.command("fix htc2 all property/global thermalCapacity peratomtype " + str(Cp))
	def initial_temperature(self, Ti, lmp):		
		lmp.command("fix htc3 all heat/gran/conduction initial_temperature " + str(Ti))
	def characteristic_velocity(self, lmp):
		lmp.command("fix m5 all property/global characteristicVelocity scalar 2.")


def make_directory(name, os):
    try:
        os.mkdir(name)
    except:
        pass

def initialize_filling_z_liggghts(geometry, lmp):
    xlim = geometry[0]
    ylim = geometry[1]
    zlim = geometry[2]
    Rp   = geometry[3]
    lmp.command('atom_style granular')
    #lmp.command('processors '+processor_layout)
    lmp.command('atom_modify map array')
    lmp.command('boundary m p m')
    lmp.command('newton off')
    lmp.command('communicate single vel yes')
    lmp.command('units si')
    lmp.command('region reg block '+str(-xlim)+' '+str(xlim)+' '+str(-ylim)+' '\
        +str(ylim)+' 0 '+str(zlim)+' units box')
    lmp.command('create_box 1 reg')
    lmp.command('neighbor '+str(2*Rp)+' bin')
    lmp.command('neigh_modify delay 0')
    lmp.command('echo both')
    lmp.command('pair_style gran model hertz tangential history')
    lmp.command('pair_coeff * *')
    lmp.command('fix gravi all gravity 9.81 vector 0.0 0.0 -1.0')
    lmp.command('variable energy equal ke')
    lmp.command('variable volume equal vol')

def initialize_filling_x_liggghts(geometry, lmp):
    xlim = geometry[0]
    ylim = geometry[1]
    zlim = geometry[2]
    Rp   = geometry[3]
    lmp.command('atom_style granular')
    #lmp.command('processors '+processor_layout)
    lmp.command('atom_modify map array')
    lmp.command('boundary m p m')
    lmp.command('newton off')
    lmp.command('communicate single vel yes')
    lmp.command('units si')
    lmp.command('region reg block '+str(-xlim)+' '+str(xlim)+' '+str(-ylim)+' '\
        +str(ylim)+' 0 '+str(zlim)+' units box')
    lmp.command('create_box 1 reg')
    lmp.command('neighbor '+str(2*Rp)+' bin')
    lmp.command('neigh_modify delay 0')
    lmp.command('echo both')
    lmp.command('pair_style gran model hertz tangential history')
    lmp.command('pair_coeff * *')
    lmp.command('fix gravi all gravity 9.81 vector 0.0 0.0 -1.0')
    lmp.command('variable energy equal ke')
    lmp.command('variable volume equal vol')

def initialize_restart_liggghts(filename, geometry, lmp):
    xlim = geometry[0]
    ylim = geometry[1]
    zlim = geometry[2]
    Rp   = geometry[3]
    lmp.command('neighbor '+str(2*Rp)+' bin')
    lmp.command('neigh_modify delay 0')
    lmp.command('pair_style gran model hertz tangential history')
    lmp.command('region reg block '+str(-xlim)+' '+str(xlim)+' '+str(-ylim)+' '\
        +str(ylim)+' 0 '+str(zlim)+' units box')
    lmp.command('read_restart '+filename)
    lmp.command('pair_coeff * *')
    lmp.command('newton off')
    lmp.command('communicate single vel yes')
    lmp.command('fix gravi all gravity 9.81 vector 0.0 0.0 -1.0')
    lmp.command('variable energy equal ke')
    lmp.command('variable volume equal vol')
    lmp.command('fix integr all nve/sphere')

    


def define_timestep(dt, lmp):
    lmp.command("timestep " + str(dt))

def create_floor(lmp):
    lmp.command('fix bottom all wall/gran model hertz '\
        'tangential history primitive type 1 zplane 0')

def create_horizontal_walls(xlim, lmp):
    lmp.command('fix xFill1 all wall/gran model hertz '\
                'tangential history primitive type 1 xplane '\
                + str(-xlim))
    lmp.command('fix xFill2 all wall/gran model hertz '\
                'tangential history primitive type 1 xplane ' \
                + str(xlim))
def destroy_filling_walls(lmp):
    lmp.command('unfix xFill1')
    lmp.command('unfix xFill2')

def create_horizontal_walls_hot(xlim, T, lmp):
    lmp.command('fix x_hot_1 all wall/gran model hertz '\
                'tangential history primitive type 1 xplane ' \
                + str(-xlim)+' temperature '+str(T))
    lmp.command('fix x_hot_2 all wall/gran model hertz '\
                'tangential history primitive type 1 xplane ' \
                + str(xlim)+' temperature '+str(T))

def create_vertical_walls(Z, lmp):
    lmp.command('fix xFill1 all wall/gran model hertz '\
        'tangential history primitive type 1 zplane 0')
    lmp.command('fix xFill2 all wall/gran model hertz '\
        'tangential history primitive type 1 zplane ' + str(Z))

def create_vertical_walls_hot(Z, T, lmp):
    lmp.command('fix z_hot_1 all wall/gran model hertz '\
                'tangential history primitive type 1 zplane 0 \
                temperature '+str(T))
    lmp.command('fix z_hot_2 all wall/gran model hertz '\
                'tangential history primitive type 1 zplane ' \
                + str(Z)+' temperature '+str(T))

def define_gaussian_pebbles(rho, mu, sigma, N, lmp):
    lmp.command('fix pts1 all particletemplate/sphere 1 atom_type 1 density constant ' \
                + str(rho) + " radius gaussian number " + str(mu) + ' ' \
                + str(sigma))
    lmp.command('fix pdd1 all particledistribution/discrete 4444 1 pts1 1.0')
    lmp.command('fix ins all insert/pack seed 4910 distributiontemplate pdd1 insert_every \
                100 maxattempt 10000 overlapcheck yes all_in yes vel constant 0. 0. 0. region reg particles_in_region '+str(N))

def insert_N_pebbles(rho, Rp, N, lmp):
    lmp.command('fix pts1 all particletemplate/sphere 1 atom_type 1 density constant '+ str(rho) + " radius constant " + str(Rp))
    lmp.command('fix pdd1 all particledistribution/discrete 4444 1 pts1 1.0')
    lmp.command('fix ins all insert/pack seed 4910 distributiontemplate pdd1 insert_every 1\
         maxattempt 1000 overlapcheck no all_in yes vel constant 0. 0. 0. region reg particles_in_region '+str(N))


def define_phi_pebbles(geometry, dump_steps, rho, phi, lmp):
    Rp   = geometry[3]
    lmp.command('fix pts1 all particletemplate/sphere 1 atom_type 1 density constant '+ str(rho) + " radius constant " + str(Rp))
    lmp.command('fix pdd1 all particledistribution/discrete 4444 1 pts1 1.0')
    lmp.command('fix ins all insert/pack seed 4910 distributiontemplate pdd1 insert_every 1\
         maxattempt 1000 overlapcheck no all_in no vel constant 0. 0. 0. region reg volumefraction_region '+str(phi))
    

def relax_insertion(dump_steps, Rp, lmp):
    lmp.command('fix nve_limit all nve/limit absolute '+str(Rp/50000.))
    lmp.command('run 1')
    lmp.command('unfix ins')
    lmp.command('run ' + str(10*dump_steps))
    lmp.command('unfix nve_limit')
    lmp.command('fix integr all nve/sphere')

def custom_screen_output(print_steps, lmp):
    lmp.command('thermo_style custom step atoms ke f_htc3 vol')
    lmp.command('thermo '+str(print_steps))

def set_all_temp(T, lmp):
    lmp.command('set region reg property/atom Temp '+str(T))


def thermal_expansion(r0, Ti, beta, grow_steps, lmp):
    # warning -- this only works with constant radius pebbles at the moment. 
    # Will need to make a command that pull's out the pebble's original radius then expands it.
	lmp.command('variable r0 equal '+str(r0))
	lmp.command('variable T0 equal '+str(Ti))
	lmp.command('variable beta equal ' + str(beta))
	lmp.command('variable growevery equal '+str(grow_steps))
	lmp.command('variable temp atom f_Temp')
	lmp.command('variable dexpand atom 2*${r0}*(1.+${beta}*(f_Temp-${T0}))')
	lmp.command('fix grow all adapt ${growevery} atom diameter v_dexpand')    

def nuke_all_pebbles(Qp, lmp):
    lmp.command('set group all property/atom heatSource '+str(Qp))



def set_dumps(output_steps, output_directory, lmp):
    lmp.command("compute fc all pair/gran/local pos id force contactArea heatFlux")
    lmp.command("run 0")
    lmp.command("dump dmp1 all custom " + str(output_steps) + ' '+\
    	        output_directory+'/dump_*.liggghts id type type x y z ix iy iz vx vy vz '+\
        	    'fx fy fz omegax omegay omegaz radius f_Temp[0] f_heatSource[0]')
    lmp.command("dump forcechain all local " + str(output_steps) + ' '+\
            	output_directory+'/dump.fc.*.liggghts c_fc[1] c_fc[2] c_fc[3] c_fc[4] '+\
            	'c_fc[5] c_fc[6] c_fc[7] c_fc[8] c_fc[9] c_fc[10] c_fc[11] '+\
            	'c_fc[12] c_fc[13] c_fc[14]')

def run_to_kinetic_steady_state(check_steps, min_energy, lmp):
	lmp.command('run '+str(check_steps))
	ke = lmp.extract_variable('energy','null',0)
	while ke > min_energy:
		lmp.command('run '+str(check_steps))
		ke = lmp.extract_variable('energy','null',0)

def get_height(Rp):
    import glob, os, numpy
    # Nzs is number of pebbles to consider for the height
    last_dump_file = str(max(glob.iglob('post/filling/dump_*.liggghts'), key=os.path.getctime))
    atomdata = numpy.loadtxt(last_dump_file, skiprows=9)
    return max(atomdata[:,5])+Rp

def get_dt(lmp):
    dt = lmp.extract_variable('step','null',0)
    return dt

def get_phi(geometry, lmp):
    xlim = geometry[0]
    ylim = geometry[1]
    zlim = geometry[2]
    Rp   = geometry[3]
    natoms = lmp.get_natoms()
    h = get_height(Rp)
    bed_volume = xlim*2*ylim*2*h
    peb_volume = (4./3) * 3.1415 * Rp**3.
    return peb_volume * natoms / bed_volume

def get_pressure(geometry, lmp):
    xlim = geometry[0]
    ylim = geometry[1]
    zlim = geometry[2]
    Rp   = geometry[3]
    force = lmp.extract_variable('top_force','null',0)
    return (force/(4*xlim*ylim))/10.**6

def hold_still(Rp, lmp):
    lmp.command('fix nve_limit_heating all nve/limit absolute '+str(Rp/10000000.))







def crush_pebbles_x(percent_to_crush, rho, Rp, lmp):
    import numpy as np
    
    natoms = lmp.get_natoms()
    
    # Pull LIGGGHTS data from dumps. Formerly this was done with extract_x but
    # many oddities pushed me into just nabbing dumps.
    import glob, os, numpy
    last_dump_file = str(max(glob.iglob('post/filling/dump_*.liggghts'), key=os.path.getctime))
    last_contact_file = str(max(glob.iglob('post/filling/dump.fc.*.liggghts'), key=os.path.getctime))
    
    forcedata = np.loadtxt(last_contact_file, skiprows=9)   
    atomdata = np.loadtxt(last_dump_file, skiprows=9)
    
 
    # Atom temperature / radius data
    atomIDs   = list(atomdata[:, 0])
    atom_temps = atomdata[:, 19]
    atom_radii = atomdata[:, 18]
    id1  = forcedata[:, 6]
    id2  = forcedata[:, 7]
    x1  = forcedata[:, 0]
    y1  = forcedata[:, 1]
    z1  = forcedata[:, 2]
    x2  = forcedata[:, 3]
    y2  = forcedata[:, 4]
    z2  = forcedata[:, 5]    


    
    import random
    reduced_id_list = list(atomIDs) # the list command makes a copy
    number_to_crush = int(np.round(percent_to_crush * natoms,0))
    atoms_to_crush = np.zeros(number_to_crush)

    for i in np.arange(0, number_to_crush):
        
        # check if the pebble is periodic. right now i have a bug where the periodic pebble
        # can not be replaced since some of the tiny pebbles will be completely outside of the system
        # so, if periodic, take a new pebble.
        periodic_flag = 1
        while periodic_flag == 1:
            atoms_to_crush[i] = random.choice(reduced_id_list)
            index = atomIDs.index(atoms_to_crush[i])
            periodic_flag = forcedata[index, 8]
            if (atomdata[index,4] < -4.5e-3 or atomdata[index,4] > 4.5e-3):
                periodic_flag = 1
            print periodic_flag



        reduced_id_list.remove(atoms_to_crush[i])



    # concatenate into a single list
    id1 = id1.tolist()
    id2 = id2.tolist()
    ids = []
    ids = id1[:]
    ids.extend(id2)
    # convert to a set and back to remove duplicates, then sort
    ids = list(set(ids))
    ids.sort()
    
    # Loop over the IDs from fc, calculate normal forces between neighbors and compare 
    # it to a Weibull distribution of forces (parameters defined above). If it gets
    # flagged as broken, its ID gets added to a group for deletion
    breakIDs = ""
    breakCount = 0
    
    for id in atoms_to_crush:
        #id = int(id)
        # Pull out the temperature & radius of the atom under consideration
        # print id
        # for ids in atomIDs:
        #     if np.abs(ids - id) < 3:
        #         print ids
        #         print atomIDs.index(ids)

        atomIndex = atomIDs.index(id)
        currentPebbleTemperature = atom_temps[atomIndex]
        Rp = atom_radii[atomIndex]
        
        # Analyze the contact data of the atom under consideration
        if id in id1:
            index = id1.index(id)
            xtemp,ytemp,ztemp = x1[index], y1[index], z1[index]
        else:
            index = id2.index(id)
            xtemp,ytemp,ztemp = x2[index], y2[index], z2[index]

        breakCount += 1
        breakIDs += str(np.int(id))+" "
        suffix = str(breakCount)
        
        # Code to insert broken fragments into region where the pebble was
        # REGION
        regionString = 'region brokenPebble' + suffix + ' sphere ' + \
                str(xtemp) + ' ' + str(ytemp) + ' ' + str(ztemp) + ' ' + \
                str(Rp) + ' units box'
        
        # PARTICLE DENSITY, SIZE DEFINITION
        particleTemplateString = 'fix brokeSphereTemplateA' + suffix + \
            ' all particletemplate/sphere 1 atom_type 1 density constant ' \
            + str(rho) + " radius gaussian number " + str(Rp/5.) + ' ' \
            + str(Rp/25.)
        # PARTICLE DISTRIBUTION TEMPLATE
        distributionTemplateString = 'fix bSTd' + suffix + \
            ' all particledistribution/discrete ' + \
            str(np.random.randint(1, 100000)) + '  1  brokeSphereTemplateA' + \
            suffix + ' 1.0'
        insertionString = 'fix bSTins' + suffix + '  all insert/pack seed '\
             + str(np.random.randint(1, 100000)) + ' distributiontemplate bSTd' \
             + suffix + ' vel constant 0. 0. 0. insert_every once ' + \
             'overlapcheck yes all_in yes volumefraction_region 1.0 ' + \
             'region brokenPebble' + suffix

        regionTempString = 'set region brokenPebble' + suffix + \
            ' property/atom Temp ' + str(currentPebbleTemperature)
        
        lmp.command(regionString)
        lmp.command(particleTemplateString)
        lmp.command(distributionTemplateString)
        lmp.command(insertionString)
        lmp.command(regionTempString)

    # Send all IDs of pebbles that broke into a single group to delete them all at once
    breakGroup = "delete_group"
    breakGroupCommand = "group "+breakGroup+" id " + breakIDs
    lmp.command(breakGroupCommand)
    lmp.command("delete_atoms group " + breakGroup + ' compress no')



def crush_pebbles_z(percent_to_crush, rho, Rp, lmp):
    import numpy as np
    
    natoms = lmp.get_natoms()
    
    # Pull LIGGGHTS data from dumps. Formerly this was done with extract_x but
    # many oddities pushed me into just nabbing dumps.
    import glob, os, numpy
    last_dump_file = str(max(glob.iglob('post/filling/dump_*.liggghts'), key=os.path.getctime))
    last_contact_file = str(max(glob.iglob('post/filling/dump.fc.*.liggghts'), key=os.path.getctime))
    
    forcedata = np.loadtxt(last_contact_file, skiprows=9)   
    atomdata = np.loadtxt(last_dump_file, skiprows=9)
    
 
    # Atom temperature / radius data
    atomIDs   = list(atomdata[:, 0])
    atom_temps = atomdata[:, 19]
    atom_radii = atomdata[:, 18]
    id1  = forcedata[:, 6]
    id2  = forcedata[:, 7]
    x1  = forcedata[:, 0]
    y1  = forcedata[:, 1]
    z1  = forcedata[:, 2]
    x2  = forcedata[:, 3]
    y2  = forcedata[:, 4]
    z2  = forcedata[:, 5]    


    
    import random
    reduced_id_list = list(atomIDs) # the list command makes a copy
    number_to_crush = int(np.round(percent_to_crush * natoms,0))
    atoms_to_crush = np.zeros(number_to_crush)

    for i in np.arange(0, number_to_crush):
        
        # check if the pebble is periodic. right now i have a bug where the periodic pebble
        # can not be replaced since some of the tiny pebbles will be completely outside of the system
        # so, if periodic, take a new pebble.
        periodic_flag = 1
        while periodic_flag == 1:
            atoms_to_crush[i] = random.choice(reduced_id_list)
            index = atomIDs.index(atoms_to_crush[i])
            periodic_flag = forcedata[index, 8]
            if (atomdata[index,4] < -7e-3 or atomdata[index,4] > 7e-3 or atomdata[index,5] < -7e-3 or atomdata[index,5] > 7e-3):
                periodic_flag = 1
            print periodic_flag



        reduced_id_list.remove(atoms_to_crush[i])



    # concatenate into a single list
    id1 = id1.tolist()
    id2 = id2.tolist()
    ids = []
    ids = id1[:]
    ids.extend(id2)
    # convert to a set and back to remove duplicates, then sort
    ids = list(set(ids))
    ids.sort()
    
    # Loop over the IDs from fc, calculate normal forces between neighbors and compare 
    # it to a Weibull distribution of forces (parameters defined above). If it gets
    # flagged as broken, its ID gets added to a group for deletion
    breakIDs = ""
    breakCount = 0
    
    for id in atoms_to_crush:
        #id = int(id)
        # Pull out the temperature & radius of the atom under consideration
        # print id
        # for ids in atomIDs:
        #     if np.abs(ids - id) < 3:
        #         print ids
        #         print atomIDs.index(ids)

        atomIndex = atomIDs.index(id)
        currentPebbleTemperature = atom_temps[atomIndex]
        Rp = atom_radii[atomIndex]
        
        # Analyze the contact data of the atom under consideration
        if id in id1:
            index = id1.index(id)
            xtemp,ytemp,ztemp = x1[index], y1[index], z1[index]
        else:
            index = id2.index(id)
            xtemp,ytemp,ztemp = x2[index], y2[index], z2[index]

        breakCount += 1
        breakIDs += str(np.int(id))+" "
        suffix = str(breakCount)
        
        # Code to insert broken fragments into region where the pebble was
        # REGION
        regionString = 'region brokenPebble' + suffix + ' sphere ' + \
                str(xtemp) + ' ' + str(ytemp) + ' ' + str(ztemp) + ' ' + \
                str(Rp) + ' units box'
        
        # PARTICLE DENSITY, SIZE DEFINITION
        particleTemplateString = 'fix brokeSphereTemplateA' + suffix + \
            ' all particletemplate/sphere 1 atom_type 1 density constant ' \
            + str(rho) + " radius gaussian number " + str(Rp/5.) + ' ' \
            + str(Rp/25.)
        # PARTICLE DISTRIBUTION TEMPLATE
        distributionTemplateString = 'fix bSTd' + suffix + \
            ' all particledistribution/discrete ' + \
            str(np.random.randint(1, 100000)) + '  1  brokeSphereTemplateA' + \
            suffix + ' 1.0'
        insertionString = 'fix bSTins' + suffix + '  all insert/pack seed '\
             + str(np.random.randint(1, 100000)) + ' distributiontemplate bSTd' \
             + suffix + ' vel constant 0. 0. 0. insert_every once ' + \
             'overlapcheck yes all_in yes volumefraction_region 1.0 ' + \
             'region brokenPebble' + suffix

        regionTempString = 'set region brokenPebble' + suffix + \
            ' property/atom Temp ' + str(currentPebbleTemperature)
        
        lmp.command(regionString)
        lmp.command(particleTemplateString)
        lmp.command(distributionTemplateString)
        lmp.command(insertionString)
        lmp.command(regionTempString)

    # Send all IDs of pebbles that broke into a single group to delete them all at once
    breakGroup = "delete_group"
    breakGroupCommand = "group "+breakGroup+" id " + breakIDs
    lmp.command(breakGroupCommand)
    lmp.command("delete_atoms group " + breakGroup + ' compress no')
