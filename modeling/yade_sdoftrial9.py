#yadedaily yade_sdoftrial9.py [FILENAME] [LINENUMBER] [batch/timing/interactive] 

#yade has a built in reader utils.readParamsFromTable
#but it only seems to work in batch mode.  seemingly no ability to
#use it interactively.  interactive mode allows us to read command
#line arguments, but batch mode doesnt! therefore build our own quick
#CSV reader AND our own batch vs interactive switch
#first argument is the filename, second argument is the line number
#if no filename is given, uses reasonable defaults.  If linenumber
#is not given by default runs the first line

#if 3rd argument is 'batch', runs in batch mode (i.e. no GUI),
#if 3rd argument is 'timing', runs a few iterations and then gives timing statistics
#otherwise interactive

#if file already exists and checkpoint exists, attempts restart from a checkpoint.  otherwise starts a new file

#system imports
import math 
import numpy as np
import csv
import os.path
import matplotlib.pyplot as pyplot
import shelve
from scipy import signal 

#yade imports
from yade import pack,export,ymport
import yade.timing
import random
import glob

sys.path.append(".") # this is not necessary in a generic python shell, but inside YADE it is for some odd reason
from common_functions import ball_gen_8_output_filename, calc_PLL_gain

#################################################
#checkpoint load/save
#
# checkpoints are automatically created stored (once per hour)
# the simulation time is in the filename
# when you reload a simulation, it creates a new output file so as to not overwrite the first one
# note that as it runs the re-started file will itself generate checkpoint files
# also if you re-run two restart file simulatenously, they will overwrite the output files
# you can change the outputs around on the reload (like enabling/disabling peak force outputs or kinetic energy outputs)
# do not try to change anything else in the spreadsheet, it probably won't work.  

def checkpoint():

        old_checkpoints = glob.glob( os.path.basename(sys.argv[1])+"_checkpoint_" + str(linenumber) + "_*")

        check_prefix = os.path.basename(sys.argv[1])+"_checkpoint_" + str(linenumber) + "_" + "{:.6g}".format(O.time)
        
        O.save( check_prefix + '.xml.bz2') 

        try:
                my_shelf = shelve.open( check_prefix + '.shelve' ,'n') # 'n' for new
        except:
                print('Couldnt open shelve file for some reason. permissions or out of disk space? skipping')
                return

                
        #only the variables in the below list get saved/reloaded in the checkpoint file. "save all" didnt work for some reason
        #but I dont remember why anymore.  If there is something you want saved, add it to this list
        for key in [ 'q1_id',  'wall_clump_list','topWall_list','botWall_list','TopwallF_accum_list', 'wallF_accum_list','f_sweep_start','f_sweep_rate','Modal_Mass','Modal_Damp','Modal_Stiff', 'F_integral','F_amp','theta','extF','base_vel','F_phase','base_displacement_amp','last_output_iter','stopTime','Psi_pocket_list','balls_mass_ratio','last_output_iter', 'ballIds', 'ballsmass_all_pockets', 'wall_listlist','SnapCount', 'Cant_nat_freq_hz', 'zf_X','zf_Y','PhaseErrorAccum','PhaseErrorPrev','PhaseErrorPrevPrev','PhaseError','PLL_Freq_rad','swing_coupled_extF_new_prev_time','PLL_prev_time', 'plotDataVirtPeriod'  ] :

                
                try:
                        my_shelf[key] = globals()[key]
                except (KeyError, TypeError):
                        #
                        # __builtins__, my_shelf, and imported modules can not be shelved.
                        #
                        print('ERROR shelving: {0}'.format(key))
        my_shelf.close()

        file1.flush()

        for f in old_checkpoints:
                os.remove(f)


def reload():
        print("remember: restarting only works if you use the same identical version of yade, python, and this script!")
        #remove both extensions
        check_prefix = os.path.splitext( os.path.splitext( old_checkpoints[0]  )[0] )[0]
        

        print('reloading checkpoint: ')
        print(check_prefix)
        O.load( check_prefix + '.xml.bz2') #fixme, custom name
        
        my_shelf = shelve.open( check_prefix + '.shelve')
        for key in my_shelf:
                #print('unshelving: ' + key)
                globals()[key]=my_shelf[key]
        my_shelf.close()



#################################################
##### process command line arguments
#################################################
        
        
if ( len(sys.argv) <= 3):
        batch = False
        timing = False
else:

        if ( sys.argv[3] == 'batch'):
                batch = True
                timing = False
        elif ( sys.argv[3] == 'timing'):
                batch = True
                timing = True
        elif ( sys.argv[3] == 'interactive'):
                batch = False
                timing = False
        else:
                raise Exception('unknown mode. choices are batch, timing, and interactive')

if ( len(sys.argv) <= 2):
        linenumber = 0
else:
        linenumber = int( sys.argv[2])-1

###############################################3

try:
        out_file_name = os.path.basename( sys.argv[1] ) + "_output_" + str(linenumber) + ".txt"
        file1 = open(out_file_name,"x") # x means create.  will raise an exception if the file already exists. prevents overwritting an existing file (which I've already done so many times!
        restart = False
except FileExistsError:
        file_size = os.path.getsize(out_file_name)
        if (file_size < 2000):
                #no real data here, just overwrite it and keep going
                file1 = open(out_file_name, "w")
                restart = False
        else:        
                #determine if restart file exists
                old_checkpoints = glob.glob( os.path.basename(sys.argv[1])+"_checkpoint_" + str(linenumber) + "_*xml.bz2")

                if (old_checkpoints):
                        restart = True
                
                        restart_number = 1        
                        out_file_name = os.path.basename( sys.argv[1] ) + "_output_" + str(linenumber) + "_restart" + str(restart_number) + ".txt"
                        while (os.path.isfile(out_file_name)):
                                restart_number = restart_number + 1
                                out_file_name = os.path.basename( sys.argv[1] ) + "_output_" + str(linenumber) + "_restart" + str(restart_number) + ".txt"
                        file1 = open(out_file_name,"x")
                else:
                        raise Exception('output file already exists and could not find a checkpoint file')

#################################################        
##### Read data from files
#################################################        

        
if (len(sys.argv) > 1):
        with open(sys.argv[1] ) as csvfile:
                reader = csv.DictReader(csvfile, skipinitialspace=True)
                row_list = list(reader)
                box_x = float( row_list[linenumber]['box_x'])
                box_y = float( row_list[linenumber]['box_y']) #this is the direction of motion
                box_z = float( row_list[linenumber]['box_z'])
                f_sweep_start = float( row_list[linenumber]['f_sweep_start'])
                f_sweep_start_rad = 2 * math.pi * f_sweep_start 
                f_sweep_stop  = float( row_list[linenumber]['f_sweep_stop'])
                f_sweep_rate  = float( row_list[linenumber]['f_sweep_rate'])
                periodic_bc = int(row_list[linenumber]['periodic_bc'])
                particle_dia = float( row_list[linenumber]['particle_dia'])
                F_amp = float( row_list[linenumber]['force_amp'])
                grav_y = float( row_list[linenumber]['grav_y'])
                damp = float( row_list[linenumber]['num_damping'])
                damp_fix = int(row_list[linenumber]['num_damp_fix'])
                damp_gravity = int(row_list[linenumber]['damp_gravity'])
                dt_in = float( row_list[linenumber]['time_step']) #Set -1 for dynamic timestepper
                verlet_dist = float( row_list[linenumber]['verlet_dist'])
                ball_gen_method = int(row_list[linenumber]['ball_gen_method'])
                
                q_factor = float(row_list[linenumber]['q_Factor'])
                
                mode_num = int(row_list[linenumber]['mode_num']) 
                modal_method = row_list[linenumber]['modal_method']
                
                
                beamL = float(row_list[linenumber]['beam_length'])
                beamH = float(row_list[linenumber]['beam_height'])
                pocket_height = float(row_list[linenumber]['pocket_height'])
                pocketLoc1 = float(row_list[linenumber]['pocketLoc1'])
                pocketLoc2 = float(row_list[linenumber]['pocketLoc2'])
                pocketLoc2_Y = float(row_list[linenumber]['pocketLoc2_Y']) #0 for off, 1 for on
                

                solution_method = int(row_list[linenumber]['solution_method'])
                force_method = int(row_list[linenumber]['force_method']) # 1 = force, 2 = base excitation (neglecting base velocity in pocket), 3 = base excitation
                base_excitation = float(row_list[linenumber]['base_excitation']) #in g's
                restitution = float(row_list[linenumber]['restitution_coeff']) #leave zero if you want to use FrictMat physics

                cn = float(row_list[linenumber]['cn']) #leave zero if you want to use FrictMat physics
                                
                wall_frict_angle = float(row_list[linenumber]['wall_frict_angle'])
                sphere_frict_angle = float(row_list[linenumber]['sphere_frict_angle'])
                distribution = row_list[linenumber]['distribution']

                include_vel_scale = int(row_list[linenumber]['include_vel_scale'])

                young = float(row_list[linenumber]['young_mod'])        
                beam_young = float(row_list[linenumber]['beam_young'])        
                kappa1 = float(row_list[linenumber]['beam_kappa'])
                p_ratio = float(row_list[linenumber]['poisson'])
                
                pocket_len =  float(row_list[linenumber]['pocket_len']) #this is the easiest way to incorporate multiple windows per pocket.
                n_windows_per_pocket = int(row_list[linenumber]['n_windows_per_pocket'])

                normalCohesion =  float(row_list[linenumber]['normalCohesion'])        
		#contactPhysics->normalAdhesion = normalCohesion * radius^2

                target_packing_factor = row_list[linenumber]['packing_factor_bg8']
                rollfrict = float(row_list[linenumber]['roll_frict'])
                timestepval = float(row_list[linenumber]['timestep_val'])
                Eroll = float(row_list[linenumber]['eta_roll'])
                angled_model = float(row_list[linenumber]['orientation'])
                
                seedType = int(row_list[linenumber]['seed'])
                sidewall_frict_on = int(row_list[linenumber]['sidewall_frict'])
                SnapCount = 0

                StokesDrag = int(row_list[linenumber]['StokesDrag'])

                DMT_gamma = float(row_list[linenumber]['DMT_gamma'])
                DMT_cutoff = float(row_list[linenumber]['DMT_cutoff'])
                
                BottomWallKFactor = float(row_list[linenumber]['BottomWallKFactor'])
                
                SettledPacking = row_list[linenumber]['settled_packing']
                mindlin = int(row_list[linenumber]['Mindlin_Type'])

                #these two are explicitly not in the restart file, so that we can change them on a restart
                use_PLL= int(row_list[linenumber]['use_PLL'])
                PLL_Qguess = float(row_list[linenumber]['q_guess'])
                PLL_filter_BW_cycles = float(row_list[linenumber]['PLL_filter_BW_cycles'])


#this needs to happen after reading the csv file, because F_amp might get overwritten
if restart:
        reload()


                
#################################################
##### DEFINING VARIABLES AND MATERIALS ######
#################################################
if (not batch):
        from yade import qt, plot
        import matplotlib.pyplot as pyplot
        qt.View() #open the controlling and visualization interfaces


nat_period = 1 /  ((f_sweep_start+f_sweep_stop)/2)

if (not restart):
        sp_frict_angle = atan(sphere_frict_angle) #contact friction during the deviatoric loading
        wl_frict_angle = atan(wall_frict_angle)
        
        rho    = 8230 #generic inco 718, kg/m^3

        assert(pocketLoc1 <= beamL)
        assert(pocketLoc2 <= beamL)

        beamB = 0.0254 #Meter


        rho_c = rho * beamB * beamH # beam linear density kg/m
        beam_I = beamB*beamH**3 / 12
        
    
    
        if (periodic_bc == 2):
                #for the 2D box, the z dimension becomes the particle diameter.                
                box_z = particle_dia

        mn = Vector3(0,      box_y,      0)
        mx = Vector3(box_x,2*box_y,  box_z) #corners of the initial packing
        thick = 2*particle_dia # the thickness of the walls


        if (ball_gen_method == 0):
                plotDataVirtPeriod = nat_period / 20
        elif (ball_gen_method == 2):
                plotDataVirtPeriod = 2 * O.dt # just one particle but we want to know a lot about it
        elif (ball_gen_method == 9):
                plotDataVirtPeriod = nat_period / 40
        else:
                plotDataVirtPeriod = nat_period / 100 #this gives 100 data points per cantilever oscillation cycle, which should be plenty
        
        extFUpdateVirtPeriod = nat_period / 1000 #update theta and external force 1000 times per osc cycle

###########################################################################################
#materials
###########################################################################################
if (not restart):
  #wanted to put zero density on walls to exclude from KE calculation, but apparently if I do that, walls fly off into space (even though they are blocked).
  #so I guess instead make the density finite but small
   if (cn > 0) or (restitution > 0) or (rollfrict > 0):
          if (normalCohesion != 0):
                  raise Exception('error, currently no way to combine viscoelasticity and cohesion')
          elif (cn > 0) and (restitution > 0):
                  raise Exception('error, use cn or en but not both')
          elif (restitution > 0 ):
                  if sidewall_frict_on == 1:
                          raise Exception("sidewall_frict_on not implemented for this case!")
                  O.materials.append(ViscElMat(young=young,poisson=p_ratio,frictionAngle=sp_frict_angle, et = 1.0, en = restitution, density=rho,  label='spheres'))
                  O.materials.append(ViscElMat(young=young,poisson=p_ratio,frictionAngle=wl_frict_angle, et = 1.0, en = restitution, density=1e-16,label='top_walls'))
                  O.materials.append(ViscElMat(young=young*BottomWallKFactor,poisson=p_ratio,frictionAngle=wl_frict_angle, et = 1.0, en = restitution, density=1e-16,label='bot_walls'))
                  O.materials.append(ViscElMat(young=young,poisson=p_ratio,frictionAngle=0,              et = 1.0, en = restitution, density=1e-16,label='sidewalls'))
          elif (cn > 0 ):
                  if sidewall_frict_on == 1:
                          raise Exception("sidewall_frict_on not implemented for this case!")
                  O.materials.append(ViscElMat(young=young,poisson=p_ratio,frictionAngle=sp_frict_angle, cs=cn, cn = cn, density=rho,  label='spheres'))
                  O.materials.append(ViscElMat(young=young,poisson=p_ratio,frictionAngle=wl_frict_angle, cs=cn, cn = cn, density=1e-16,label='top_walls'))
                  O.materials.append(ViscElMat(young=young*BottomWallKFactor,poisson=p_ratio,frictionAngle=wl_frict_angle, cs=cn, cn = cn, density=1e-16,label='bot_walls'))
                  O.materials.append(ViscElMat(young=young,poisson=p_ratio,frictionAngle=0,              cs=cn, cn = cn, density=1e-16,label='sidewalls'))

          elif (rollfrict > 0):
                  if sidewall_frict_on == 1:
                          raise Exception("sidewall_frict_on not implemented for this case!")
                  print('using ViscElMat')
                  O.materials.append(ViscElMat(young=young,poisson=p_ratio,frictionAngle=sp_frict_angle, mR = rollfrict, mRtype = 2, cs=0.0, cn = 0.0, density=rho,  label='spheres'))
                  O.materials.append(ViscElMat(young=young,poisson=p_ratio,frictionAngle=wl_frict_angle, mR = rollfrict, mRtype = 2, cs=0.0, cn = 0.0, density=1e-16,label='top_walls'))
                  O.materials.append(ViscElMat(young=young*BottomWallKFactor,poisson=p_ratio,frictionAngle=wl_frict_angle, mR = rollfrict, mRtype = 2, cs=0.0, cn = 0.0, density=1e-16,label='bot_walls'))
                  O.materials.append(ViscElMat(young=young,poisson=p_ratio,frictionAngle=0, mR = rollfrict, mRtype = 2, cs=0.0, cn = 0.0, density=1e-16,label='sidewalls'))
   
   else:
          if (normalCohesion != 0) or (Eroll !=0):
                  print('using cohfrictmat')
                  O.materials.append(CohFrictMat(young=young,poisson=p_ratio,frictionAngle=sp_frict_angle,density=rho, etaRoll = Eroll, etaTwist = Eroll, isCohesive=True, momentRotationLaw=True,normalCohesion= normalCohesion, shearCohesion=normalCohesion, label='spheres'))
                  O.materials.append(CohFrictMat(young=young,poisson=p_ratio,frictionAngle=wl_frict_angle,density=1e-16, etaRoll = Eroll, etaTwist = Eroll, isCohesive=True, momentRotationLaw=True,normalCohesion= normalCohesion, shearCohesion=normalCohesion, label='top_walls'))
                  O.materials.append(CohFrictMat(young=young*BottomWallKFactor,poisson=p_ratio,frictionAngle=wl_frict_angle,density=1e-16, etaRoll = Eroll, etaTwist = Eroll, isCohesive=True, momentRotationLaw=True,normalCohesion= normalCohesion, shearCohesion=normalCohesion, label='bot_walls'))
                  #including cohesion on the sidewall or not doesn't seem to make much differences (as long as incorporated into swing correctly).  If no real difference,
                  #more physically justifiable to not include it.  real box is ~700x700 particles, so the wall conditions are probably negligible
                  if sidewall_frict_on == 0:
                      O.materials.append(CohFrictMat(young=young,poisson=p_ratio,frictionAngle=0             ,density=1e-16, etaRoll = Eroll, etaTwist = Eroll, isCohesive=False, momentRotationLaw=True, normalCohesion=0, shearCohesion=0, label='sidewalls'))
                  else:
                      O.materials.append(CohFrictMat(young=young,poisson=p_ratio,frictionAngle=wl_frict_angle, density=1e-16, etaRoll = Eroll, etaTwist = Eroll, isCohesive=True, momentRotationLaw=True, normalCohesion=normalCohesion, shearCohesion=normalCohesion, label='sidewalls'))
          else:
                  if sidewall_frict_on == 1:
                          raise Exception("sidewall_frict_on not implemented for this case!")
                  #"poisson" here is a hidden flag.  ACTUALLY the ratio between normal and shear stiffness
                  #young here is a hidden flag.  Actually a contact stiffness.  F = E * D * d
                  #for bc 0, same material is used for all 6 walls including friction. for 2 and 7, friction is only used for top & bottom walls, no friction for others
                  #this is only for backwards compatibility.  never really intended to have friction on side walls
                  print('using FrictMat')
                  O.materials.append(FrictMat(young=young,poisson=p_ratio,frictionAngle=sp_frict_angle,density=rho,label='spheres'))
                  O.materials.append(FrictMat(young=young,poisson=p_ratio,frictionAngle=wl_frict_angle,density=1e-16  ,label='top_walls'))
                  O.materials.append(FrictMat(young=young*BottomWallKFactor,poisson=p_ratio,frictionAngle=wl_frict_angle,density=1e-16  ,label='bot_walls'))
                  O.materials.append(FrictMat(young=young,poisson=p_ratio,frictionAngle=0,             density=1e-16  ,label='sidewalls'))



##################################################################################################
#GENERATING BODIES #####                

#ball gen methods
#0: no particles at all, for testing fully fused or experimental Q
#1: loose packing with yade.pack.SpherePack().  not recommended
#2: only only 1 single ball. only for debugging
#3: yade.pack.randomDensePack doesnt work for normal dist, not really dense enough
#4: gravity deposition
#5: gravity deposition with box shaking
#6: regular hexagonal packing.  doesnt work for norml dist. never really used
#7: gravity deposition like 5, but start with an initial velocity to improve convergence speed. never really used much
#8: overlapping random generation which is then let to settle out under gravity
#9: 1d: one big clump of particles.  essentially generates a single cylinder 
#10 1D stack of particles (not a clump)
#11: like 8, but clump together some of the particles to create non-spherical particles and reduce flowability
#12: unused?
#13: read a pre-packed csv file from a previous simulation
#14: like 8, but with a pressure BC representing the weight of the powder above.  Used for U Tennessee model, not for AFRL model.

def ball_generate(box_off_x,box_off_y):
        
    if (distribution == 'uniform'):
            dia = [particle_dia, particle_dia * 1.001 ]  #these are the points on a distribution
            phi = [0.0,   1.0]
    elif (distribution == 'normal'):
            std_dev = 0.39 * particle_dia #fixme make an input parameter
            #this is 10 points from -2 std dev to 2 std dev, rounded off to 0 and 1 at the ends
            dia = np.array([-2.0000,-1.5556,-1.1111,-0.6667,-0.2222, 0.2222, 0.6667, 1.1111, 1.5556, 2.0000]) * std_dev + particle_dia
            phi =          [0.0, 0.0599, 0.1333, 0.2525, 0.4121, 0.5879, 0.7475, 0.8667, 0.9401, 1.0]
    
            
    box_off = (box_off_x, box_off_y, 0)
    
    ballIds = []
    if (ball_gen_method==1):
        #use a SpherePack object to generate a random **loose** particles packing
        #typical around 25% packing factor, which is way too loose.
        
        sp=yade.pack.SpherePack()
        sp.makeCloud(minCorner=mn+box_off,maxCorner=mx+box_off,periodic=False,psdSizes=dia,psdCumm=phi,distributeMass=False,seed=1)
        # "seed" makes the "random" generation always the same
        ballIds = O.bodies.append([sphere(center,rad,material='spheres',color=(0,0,1)) for center,rad in sp])


        if (distribution == 'uniform'):
                timestepper =  GlobalStiffnessTimeStepper(active=True,timeStepUpdateInterval=100,timestepSafetyCoefficient=0.8)
                #timestepper =  GlobalStiffnessTimeStepper(active=True,timeStepUpdateInterval=100,timestepSafetyCoefficient=0.8, viscEl=True) This method was not working -- viscEl did not have any impact on adjustable timestepper.

        else:
                #time stepper not working good for non-uniform distributions. just set fixed based on smallest particle
                timestepper =  GlobalStiffnessTimeStepper(active=False)
                O.dt = 0.95* yade.utils.SpherePWaveTimeStep( dia[0]/2 , rho, young)


    elif (ball_gen_method == 0):
        #turn off all spheres.  this is only for debugging

        # without any spheres we need to set time step manually
        timestepper =  GlobalStiffnessTimeStepper(active=False,timeStepUpdateInterval=100,timestepSafetyCoefficient=0.8)

        #was doing 100 here. phase calculations were off. think amplitude was okay
        O.dt = nat_period  / 1000

    elif (ball_gen_method == 2):
        #only only 1 single ball.  this is also only for debugging
        #currently broken and I don't know why
        center = ( (mn[0]+mx[0])/2 + box_off_x, (mn[1]+mx[1])/2 + box_off_y, (mn[2]+mx[2])/2  )
        rad = particle_dia/2
        ballIds = [ O.bodies.append(sphere(center,rad,material='spheres',color=(0,0,1))) ]

        timestepper =  GlobalStiffnessTimeStepper(active=False,timeStepUpdateInterval=100,timestepSafetyCoefficient=0.8)
        O.dt = nat_period / 10000
    elif (ball_gen_method == 3):
        #use a random dense pack. originally thought to add some some headroom.  never got there.  think this might be too loose as is
        if (periodic_bc == 2):
                raise Exception('random dense pack is not working in 2D')
            
        #        headroom_fraction = 0.01
        #        headroom = Vector3(0,0, headroom_fraction * box_z)
        #headroom = Vector3(0,0,0)
        
        #predicate = yade._packPredicates.inAlignedBox( mn, mx - headroom)
        predicate = yade._packPredicates.inAlignedBox( mn+box_off, mx+box_off)

        #spheresInCell makes this go *much* quicker. memoizeDb saves results to file for later use
        sp = yade.pack.randomDensePack(predicate, radius=particle_dia/2, material='spheres', color=(1,0,0), seed=1, returnSpherePack=True, spheresInCell=5000, memoizeDb='packing_database' )
        
        #sp.makeCloud(minCorner=mn,maxCorner=mx,periodic=False,psdSizes=dia,psdCumm=phi,distributeMass=True,seed=1)
        # "seed" makes the "random" generation always the same
        
        ballIds = O.bodies.append([sphere(center,rad,material='spheres',color=(0,0,1)) for center,rad in sp])

        timestepper =  GlobalStiffnessTimeStepper(active=True,timeStepUpdateInterval=100,timestepSafetyCoefficient=0.8)
        
    elif (ball_gen_method == 4) or (ball_gen_method == 5)  or (ball_gen_method == 7)  or (ball_gen_method == 8)  or (ball_gen_method == 11) :
        #pre-compute the packing with a separate script and store to a filename
        #filename is automatically computed based on the relavent parameters


        if (ball_gen_method == 8):
            file_name = ball_gen_8_output_filename( target_packing_factor, particle_dia, box_x, box_y, box_z, seedType, periodic_bc, distribution, damp_gravity, std_dev)
        else:            
                if (ball_gen_method == 5):
                        shake = "shake_"
                elif (ball_gen_method == 7):
                        shake = "gen7_"
                elif  (ball_gen_method == 11):
                        shake = "gen8_" + str(target_packing_factor) + "_"
                else:
                        shake = ""
                        
                if (distribution == 'uniform'):
                        file_name = "grav_dep_pack_" + str(np.floor(box_x/1e-6)) + "_" + str(np.floor(box_y/1e-6)) + "_" + str(np.floor(box_z/1e-6)) + "_uniform_" + shake + str(np.floor(particle_dia/1e-6)) + ".csv"
                        clump_name = "clump_" + str(np.floor(box_x/1e-6)) + "_" + str(np.floor(box_y/1e-6)) + "_" + str(np.floor(box_z/1e-6)) + "_uniform_" + shake + str(np.floor(particle_dia/1e-6)) + ".csv" #for ball gen 11 only
                elif (distribution == 'normal') and (seedType == 1):
                        file_name = "grav_dep_pack_" + str(np.floor(box_x/1e-6)) + "_" + str(np.floor(box_y/1e-6)) + "_" + str(np.floor(box_z/1e-6)) + "_normal_" + shake + str(np.floor(particle_dia/1e-6)) + "_" + str(np.floor(std_dev/1e-6)) + ".csv"
                        clump_name = "clump_" + str(np.floor(box_x/1e-6)) + "_" + str(np.floor(box_y/1e-6)) + "_" + str(np.floor(box_z/1e-6)) + "_normal_" + shake + str(np.floor(particle_dia/1e-6)) + "_" + str(np.floor(std_dev/1e-6)) + ".csv"
                elif (distribution == 'normal') and (seedType != 1):
                        file_name = "grav_dep_pack_" + str(np.floor(box_x/1e-6)) + "_" + str(np.floor(box_y/1e-6)) + "_" + str(np.floor(box_z/1e-6)) + "_normal_" + shake + str(np.floor(particle_dia/1e-6)) + "_" + str(np.floor(std_dev/1e-6)) + "_Seed" +str(seedType)+".csv"
                        clump_name = "clump_" + str(np.floor(box_x/1e-6)) + "_" + str(np.floor(box_y/1e-6)) + "_" + str(np.floor(box_z/1e-6)) + "_normal_" + shake + str(np.floor(particle_dia/1e-6)) + "_" + str(np.floor(std_dev/1e-6)) + ".csv"

            
        ballIds = []

        print("opening grav dep file: " + file_name + "\n")
        
        with open(file_name ) as csvfile:
                reader = csv.DictReader(csvfile, skipinitialspace=True)
                row_list = list(reader)
                for row in row_list:
                        x   = float( row['x'])
                        y   = float( row['y'])
                        z   = float( row['z'])
                        rad = float( row['rad'])

                        ball = O.bodies.append([sphere((x+box_off_x,y+box_off_y,z),rad,material='spheres',color=(0,0,1))])
                        ballIds.append(ball[0])

        if (ball_gen_method == 11):
                if (periodic_bc == 2):
                        raise Exception('not tested in 2D.  not sure if need to block the DOF on the clump or not. needs tested')
                
                with open(clump_name ) as csvfile:
                        reader = csv.DictReader(csvfile, skipinitialspace=True)
                        row_list = list(reader)
                        for row in row_list:
                                m1   = int( row['m1'])
                                m2   = int( row['m2'])
                                O.bodies.clump([m1, m2])
                        
        if (distribution == 'uniform'):
                timestepper =  GlobalStiffnessTimeStepper(active=True,timeStepUpdateInterval=100,timestepSafetyCoefficient=0.8)
        else:
            if (rollfrict>0):
                timestepper =  GlobalStiffnessTimeStepper(active=False)
                O.dt = timestepval
#                O.dt = 1e-9
#                timestepper =  GlobalStiffnessTimeStepper(active=True,viscEl=True, timestepSafetyCoefficient=0.1)

            else:
                #time stepper not working good for non-uniform distributions. just set fixed based on smallest particle
                timestepper =  GlobalStiffnessTimeStepper(active=False)
                O.dt = 0.8* yade.utils.SpherePWaveTimeStep( dia[0]/2 , rho, young)
#                timestepper =  GlobalStiffnessTimeStepper(active=True,viscEl=True, timestepSafetyCoefficient=0.1)
                
        if angled_model > 0:
             axis = (0,0,1)
             center = Vector3(0,0,0)
             angl = angled_model*pi/18
             q = Quaternion(axis,angl)
        
             for s in ballIds:
                p = O.bodies[s].state.pos
                newP = center+q*(p-center)
                O.bodies[s].state.pos=O.bodies[s].state.refPos=newP
                O.bodies[s].state.ori = q

                
                
    elif (ball_gen_method == 13):
        ballIds = []
        kw={'material':0}
        temp = O.bodies.append(ymport.text('BallGen_Settled_BallPos_'+SettledPacking,shift=Vector3(box_off_x,box_off_y,0),**kw))
        O.bodies.erase(temp[-1]) #IMPORTANT, REMOVES Q1 GENERATED PARTICLE. Note the order of how it works: Balls get generated, then walls, then clump, then q1
        del temp[-1]
#        print(temp[-1])

        for z in range(len(temp)):
            ballIds.append(temp[z])

        timestepper =  GlobalStiffnessTimeStepper(active=False)
        O.dt = 0.8* yade.utils.SpherePWaveTimeStep( dia[0]/2 , rho, young)
        
        
    elif (ball_gen_method == 6):
        #If 2d do 2d        
        
        #If 3d do this
        #Using the hexagonal packing!
        ballIds = O.bodies.append(pack.regularHexa(pack.inAlignedBox(mn,mx),material = 'spheres', radius=particle_dia/2,gap=0,color=(0,1,0)),)
        if (distribution == 'uniform'):
                timestepper =  GlobalStiffnessTimeStepper(active=True,timeStepUpdateInterval=100,timestepSafetyCoefficient=0.8)
        else:
                raise ValueError('Cannot use a non-uniform packing distribution with hexagonal packing!')                

    elif (ball_gen_method == 9) or (ball_gen_method == 10) :
        #1D just a stack of balls as one big clump
        if (periodic_bc != 8):
                raise Exception("ball gen 9+10 only works with bc 8")

        if (distribution != "uniform"):
                raise Exception("ball gen 9+10 only works with unifom size")
            
        nballs = int( box_y / particle_dia)

        for i in range(0,nballs):
                center = ( (mn[0]+mx[0])/2 + box_off_x, mn[1] + particle_dia/2 + particle_dia * i ,  (mn[2]+mx[2])/2 )
                id = O.bodies.append(sphere(center, particle_dia/2 ,material='spheres',color=(0,0,1)))
                ballIds.append(id)

        if (ball_gen_method == 9):
                #one big clump
                ballClump=O.bodies.clump(ballIds)
                O.bodies[ballClump].state.blockedDOFs = 'xzXYZ'
                timestepper =  GlobalStiffnessTimeStepper(active=False)
                O.dt = 0.8* yade.utils.SpherePWaveTimeStep( particle_dia/2 , rho, young)

        else:
                #one stack
                for i in ballIds:
                        O.bodies[i].state.blockedDOFs = 'xzXYZ'
                timestepper =  GlobalStiffnessTimeStepper(active=False)
                O.dt = 0.1* yade.utils.SpherePWaveTimeStep( particle_dia/2 , rho, young)
        

        #timestepper =  GlobalStiffnessTimeStepper(active=True,timeStepUpdateInterval=100,timestepSafetyCoefficient=0.8)
    else:
            raise Exception('Unknown ball_gen_method')
            
    return ballIds, timestepper

#########################################################################################################################

if (not restart):
        ballIds = []



        for i in range(n_windows_per_pocket):
                ballIds1, timestepper  = ball_generate( i * 2.5 * box_x, i * 2.5 * box_y  )
                ballIds = ballIds + ballIds1


        if (periodic_bc == 2):
                #area packing factor:
                balls_area = 0
                for i in ballIds:
                        balls_area = balls_area + math.pi * O.bodies[i].shape.radius**2

                box_area = box_x * box_y
        
                relative_density = balls_area / box_area
                print("area packing factor: " + str(relative_density))
        
        else:
                #volume packing factor
                ballsmass_1pocket = 0
                for i in ballIds:
                        ballsmass_1pocket = ballsmass_1pocket + O.bodies[i].state.mass
        
                ballsvolume = ballsmass_1pocket / rho
                boxvolume   = box_x * box_y * box_z
                relative_density = ballsvolume / boxvolume
                print("volume packing factor: " + str(relative_density))


        if ( pocketLoc2 > 0 ) and (pocketLoc2_Y == 0):
                for i in range(n_windows_per_pocket):
                        ballIds1, _  = ball_generate( ( i + n_windows_per_pocket) * 2.5 * box_x, 0 )
                        ballIds = ballIds + ballIds1
                        
        if (pocketLoc2_Y > 0):
                for i in range(n_windows_per_pocket):
                        ballIds1, _  = ball_generate(0, ( i + n_windows_per_pocket) * 2.5 * box_y )
                        ballIds = ballIds + ballIds1
  
        if (dt_in > 0):
                timestepper.active=False
                O.dt = dt_in


        ballsmass_all_pockets = 0
        for i in ballIds:
                ballsmass_all_pockets = ballsmass_all_pockets + O.bodies[i].state.mass

        
####################################################################
#galerkin stuff
if (not restart):
  if (modal_method == 'galerkin'):
          if (mode_num > 0):
                  #use Euler-bernoulli beam theory to generate modal parameters for a single mode
                  roots = np.array([  1.8751, 4.6941, 7.8548, 10.9955, 14.1372  ])
                  idx = mode_num-1
                  #                
                  x = np.linspace(0,beamL,1000); 

                  sigma = (math.sin(roots[idx]) + math.sinh(roots[idx])) / (math.cos(roots[idx]) + math.cosh(roots[idx]))
  
                  def Psi(var):
                          return np.sin(roots[idx] * var / beamL) - np.sinh(roots[idx]*var/beamL) - sigma * (np.cos(roots[idx]*var/beamL) - np.cosh(roots[idx]*var/beamL))

                  #since we may now have multiple pockets, normalization to pocket location does not make sense.  always normalize to the tip.
                  Psi_tip = Psi(beamL) 

                  def Psi_n(var):
                          return Psi(var) / Psi_tip

                  Psi_Norm = Psi_n(x)

                  Psi_pocket_list = []

                  for i in range( n_windows_per_pocket):
                          if (n_windows_per_pocket == 1):
                                  loc =  pocketLoc1
                          else:
                                  loc =  pocketLoc1 - pocket_len/2 + i * pocket_len / (n_windows_per_pocket-1)
                        
                          Psi_pocket_list.append( Psi_n(loc))
                          print( 'location: ' + str(loc) + " psi: " + str(Psi_n(loc)))
                
                  if ( pocketLoc2 > 0) and (pocketLoc2_Y == 0):
                          for i in range( n_windows_per_pocket):
                                  if (n_windows_per_pocket == 1):
                                          loc =  pocketLoc2
                                  else:
                                          loc =  pocketLoc2 - pocket_len/2 + i * pocket_len / (n_windows_per_pocket-1)
                        
                                  Psi_pocket_list.append( Psi_n(loc))
                                  print( 'location: ' + str(loc) + " psi: " + str(Psi_n(loc)))
                  elif ( pocketLoc2_Y > 0):
                          for i in range( n_windows_per_pocket):
                                  if (n_windows_per_pocket == 1):
                                          loc =  pocketLoc2
                                  else:
                                          loc =  pocketLoc2 - pocket_len/2 + i * pocket_len / (n_windows_per_pocket-1)
                        
                                  Psi_pocket_list.append( Psi_n(loc))
                                  print( 'location: ' + str(loc) + " psi: " + str(Psi_n(loc)))
                  Psi_Norm_sq = Psi_Norm**2
                  Modal_Mass = float(rho_c * np.trapz(Psi_Norm_sq,x))

                  #second derivative of the mode shape
                  temp = (roots[idx]**2 * (-np.sin(roots[idx]*x/beamL) + sigma*np.cos(roots[idx]*x/beamL) -np.sinh(roots[idx]*x/beamL) + sigma*np.cosh(roots[idx]*x/beamL))/(beamL**2 * Psi_tip))**2


#                  beam_young = 200e9 #note not the same as the ball-ball young's modulus, which might be a scaled value

                  Modal_Stiff = float( beam_young * beam_I * np.trapz(temp,x))

                  F_integral = float(np.trapz(Psi_Norm,x) / np.trapz(Psi_Norm_sq,x)) #base excitation
          elif (mode_num == -1):
                  # use Euler-bernoulli beam theory to generate modal parameters for a multiple modes
                  raise Exception('multi-mode branch was never merged!')
          elif (mode_num == -2):
                  #directly input values computed from ansys
                  #old way: hard coded to thin beam, mode 3, pocket at 0.050
                  if ( pocketLoc2 > 0):
                          raise Exception('only implemented for one pocket')
                  if ( n_windows_per_pocket > 1):
                          raise Exception('only implemented for one window')
                  if (pocket_len != 0.0175):
                          raise Exception('only implemented for thin beam')
                  if (pocketLoc1 != 0.05):
                          raise Exception('only implemented for pocket at 0.05')

                  raise Exception('this method is deprecated. use modal_method=ansys instead')
                  Modal_Mass=0.02557
                  Modal_Stiff=4088889
                  F_integral=0.49702
                  Psi_pocket_list=[0.784]
          else:
                  raise Exception('unknown mode_num choice?')
  elif (modal_method == "ansys"):
          if ( n_windows_per_pocket > 1):
                  raise Exception('only implemented for one window')
        
          ansys_file_name = "ansys_" + str(beamL) + "_" + str(pocketLoc1) + "_" + str(pocketLoc2) + "_" + str(mode_num) + ".txt"
          with open(ansys_file_name ) as csvfile:
                  reader = csv.DictReader(csvfile, skipinitialspace=True)
                  row_list = list(reader)

                  Modal_Mass=float( row_list[0]['Modal_Mass'])
                  Modal_Stiff=float( row_list[0]['Modal_Stiff'])
                  F_integral=float( row_list[0]['F_integral'])
                  Psi_pocket_list=[ float( row_list[0]['Psi_pocket_list_1_1']) ]  #trailing 1 is for future use to implement multiple windows. it will index the windows
                  if (pocketLoc2 > 0):
                          Psi_pocket_list.append(  float( row_list[0]['Psi_pocket_list_2_1']) )

  elif (modal_method == "galerkin2"):
      k1 = 1e10 #linear spring stiffness
#      kappa1 = 7.5e3 #torsional spring stiffness
#      E_mod = 190e9 #Beam bulk young modulus, based on Onome verbal
      kW = k1 * beamL**3 / (beam_young * beam_I)
      kth = kappa1 * beamL / (beam_young * beam_I)
      Cross_Area = beamB*beamH
      if (mode_num > 0):
                  #use Euler-bernoulli beam theory to generate modal parameters for a single mode
                  if (beamH != 0.00635 and beamH != 0.0032):
                      raise Exception('A beam dimension has been chosen that is not a thick or thin beam geometry. This modal method (galerkin2) is not suited for this outside case!')
                  elif beamH == 0.00635:
                      if kappa1 != 7500:
                          raise Exception('Kappa must be 7500 for this beam height')
                      #Also check the kappa value to make sure the spreadsheet is using 7.5e3 and not another value.
                      #Throw an exception here
                      
                      print('Using premade root values for kappa of 7.5e3N-m/rad here')
                      roots = np.array([  1.7524, 4.4453, 7.5026, 10.5716, 13.6594]) #For thick beam 7.5e3 kappa
                  else:
                      if kappa1 != 1250:
                          raise Exception('Kappa must be 1250 for this beam height')
                      #Also check the kappa value to make sure the spreadsheet is using 1.25e3 and not another value.
                      #Throw an exception here
                      print('Using premade root values for kappa of 1.25e3N-m/rad here')
                      roots = np.array([  1.7618, 4.4609, 7.5216, 10.5918, 13.6806 ]) #For thin beam 1.25e3 kappa
                  
                  idx = mode_num-1
                  x = np.linspace(0,beamL,1000); 

                  sigma = ((roots[idx]**3 * (math.cos(roots[idx]) - math.cosh(roots[idx]))) + kW * (math.sin(roots[idx]) + math.sinh(roots[idx]))) / ((roots[idx]**3 * (math.sin(roots[idx]) + math.sinh(roots[idx]))) - kW * (math.cos(roots[idx]) + math.cosh(roots[idx])))
  
                  def Psi(var):
                          return (np.sin(roots[idx] * (1 - (var / beamL))) + np.sinh(roots[idx]*(1 - (var/beamL)))) + sigma * (np.cos(roots[idx]*(1 - (var/beamL))) + np.cosh(roots[idx]*(1 - (var/beamL))))


                  #since we may now have multiple pockets, normalization to pocket location does not make sense.  always normalize to the tip.
                  Psi_tip = Psi(beamL) 

                  def Psi_n(var):
                          return Psi(var) / Psi_tip

                  Psi_Norm = Psi_n(x)

                  Psi_pocket_list = []

                  for i in range( n_windows_per_pocket):
                          if (n_windows_per_pocket == 1):
                                  loc =  pocketLoc1
                          else:
                                  loc =  pocketLoc1 - pocket_len/2 + i * pocket_len / (n_windows_per_pocket-1)
                        
                          Psi_pocket_list.append( Psi_n(loc))
                          print( 'location: ' + str(loc) + " psi: " + str(Psi_n(loc)))
                
                  if ( pocketLoc2 > 0) and (pocketLoc2_Y == 0):
                          for i in range( n_windows_per_pocket):
                                  if (n_windows_per_pocket == 1):
                                          loc =  pocketLoc2
                                  else:
                                          loc =  pocketLoc2 - pocket_len/2 + i * pocket_len / (n_windows_per_pocket-1)
                        
                                  Psi_pocket_list.append( Psi_n(loc))
                                  print( 'location: ' + str(loc) + " psi: " + str(Psi_n(loc)))
                  elif ( pocketLoc2_Y > 0):
                          for i in range( n_windows_per_pocket):
                                  if (n_windows_per_pocket == 1):
                                          loc =  pocketLoc2
                                  else:
                                          loc =  pocketLoc2 - pocket_len/2 + i * pocket_len / (n_windows_per_pocket-1)
                        
                                  Psi_pocket_list.append( Psi_n(loc))
                                  print( 'location: ' + str(loc) + " psi: " + str(Psi_n(loc)))
                  Psi_Norm_sq = Psi_Norm**2
                  Modal_Mass = float(rho_c * np.trapz(Psi_Norm_sq,x))

                  #second derivative of the mode shape
                  temp = (-1*(roots[idx]**2 * (sigma * np.cos(roots[idx]*(1-x/beamL)) - sigma * np.cosh(roots[idx]*(1 - x/beamL)) + np.sin(roots[idx]*(1-x/beamL)) - np.sinh(roots[idx]*(1-x/beamL))) / (beamL**2 * Psi_tip)))**2



                  Modal_Stiff = float( beam_young * beam_I * np.trapz(temp,x))
                  F_integral = float(np.trapz(Psi_Norm,x) / np.trapz(Psi_Norm_sq,x)) #base excitation
      else:
                  raise Exception('multi-mode for torsional boundary condition never created!')
      
      
  Modal_Damp = float(math.sqrt( (Modal_Mass * Modal_Stiff )) / q_factor)

  Cant_nat_freq_rad = math.sqrt( Modal_Stiff / Modal_Mass)
  Cant_nat_freq_hz  = math.sqrt( Modal_Stiff / Modal_Mass) / (2 * math.pi)
#  b = roots[idx]
#  Natfreq = (b**2 / beamL**2 / math.sqrt(rho * Cross_Area / (E_mod*beam_I))) / 2*math.pi

  print("Modal Mass: " + str(Modal_Mass) + " kg\nModal Stiffness: " + str(Modal_Stiff) + " N/m\nModal Damping: " + str(Modal_Damp) + " N*s/m\n")
  print("Cant nat freq = " + str(Cant_nat_freq_hz) + " Hz\n");

        
  galerkin_file_name = os.path.basename( sys.argv[1] ) + "_galerkin_" + str(linenumber) + ".txt"
  file3 = open(galerkin_file_name,"w+") 
  file3.write("ModalMass, ModalStiffness, ModalDamping, f_sweep_start, f_sweep_stop, F_amp, F_integral, base_exc_g\n")
  file3.write(str(Modal_Mass) + "," + str(Modal_Stiff) + "," + str(Modal_Damp) + "," + str(f_sweep_start) + "," + str(f_sweep_stop) + "," + str( np.sqrt( (Modal_Mass)**2 + (Modal_Damp/(f_sweep_start*2*math.pi))**2)* base_excitation*9.81 * F_integral) + ", "+ str(F_integral)+ ", "+ str(base_excitation)+ "\n")
  file3.close()

#########################################################

if (not restart):
  #this allows us to simulate only a portion of the pocket, and then scale up the forces
  #balls_mass_ratio   = 0.0178 * 0.0013 * (0.0178 / n_windows_per_pocket)  / (box_x * box_y * box_z)
  balls_mass_ratio   = pocket_len * pocket_height * (pocket_len / n_windows_per_pocket)  / (box_x * box_y * box_z)

  if (include_vel_scale == 1):
          raise Exception('this method was deleted')

  q1_dummy_volume = (4/3) * math.pi * (particle_dia/2)**3
  q1_dummy_density = Modal_Mass / q1_dummy_volume

  O.materials.append(ElastMat(young=young,poisson=p_ratio, density=q1_dummy_density ,label='dummyq1'))
    
  q1_id = O.bodies.append(sphere( (-box_x, 0, -box_z), particle_dia/2, material='dummyq1',color=(0,1,0)))

  if (solution_method == 3):
          #proscribed displacement.  block (so not affected by gravity) and then apply velocity later
          O.bodies[q1_id].state.blockedDOFs = 'xyzXYZ'
  
#define here so we don't have to de-reference the object every time through swing(), saves ~8 microseconds in a time critical loop
q1_state   = O.bodies[q1_id].state


####################################################################
#generating walls

#boundary condition options
#0: normal 3D space, all 6 walls are the same
#1: 3D with periodic boundary conditions (not sure if it works correctly)
#2: 2D (4 walls)
#7: normal 3D space, top & bottom walls have different material than side walls
#8: 1D space (2 walls)
#9: 3D space, non-smooth top & bottom walls (branch not merged as made things more complicated and didnt improve correlation)

def wall_generate(box_off_x,box_off_y):

    box_off = (box_off_x, box_off_y,0)
    
    if (periodic_bc==1):
        if (box_off_x > 0):
            raise Exception('cant have multiple walls with periodic boundary conditions')
        
        #3D space with periodic boundary conditions. not sure this works right.
        
        #yade's periodic boundaries only really work correctly for particles than are less than half of the cell
        #width. something like a wall that spans the entire length of the cell doesnt work.
        #when the particle crosses the periodic boundary, it's position is not wrapped around.  it just keeps
        #increasing, like there was no boundary. but then the collider applies a modulo operation when it determines
        #interactions. but that modulo is not applied to particles bigger than the cell width.  So, as a workaround,
        #the wall has to be much much much bigger than the cell (and technically speaking this could still fail
        # if a particle happened to wrap around the boundar 1000 times, although that's unlikely)
        
        expand=1000


        # yade.utils.box(  center                             extents= half sizes!
        topWall   = yade.utils.box( (box_x/2,2*box_y+thick/2,box_z/2), (box_x/2*expand,thick/2,box_z/2*expand), material='top_walls')
        topWallId =  O.bodies.append(topWall)
        wallIds = [topWallId]

        botWall   = yade.utils.box( (box_x/2,box_y-thick/2,box_z/2), (box_x/2*expand,thick/2, box_z/2*expand), material='bot_walls')
        botWallId =  O.bodies.append(botWall)
        wallIds.append(botWallId)
        
        O.periodic = True
        
        O.cell.setBox(  Vector3(box_x,3*box_y,box_z) )
    elif (periodic_bc==0):
            raise Exception('periodic_bc=0 removed')        
    elif (periodic_bc==2):
        #2D space.  4 walls
        #special 4 walls case for 2D because of all of the blocked DOF
        expand=1.5
            
                  # yade.utils.box(  center                             extents= half sizes!
        topWall   = yade.utils.box( ( box_off_x+box_x/2,2*box_y+thick/2,box_z/2), (box_x/2*expand,thick/2,box_z/2*expand), material='top_walls')
        topWallId =  O.bodies.append(topWall)
        wallIds = [topWallId]

        botWall   = yade.utils.box( ( box_off_x+box_x/2,box_y-thick/2,box_z/2), (box_x/2*expand,thick/2, box_z/2*expand), material='bot_walls')
        botWallId =  O.bodies.append(botWall)
        wallIds.append(botWallId)
            
            
        rightWall   = yade.utils.box( ( box_off_x+box_x+thick/2,1.5*box_y,box_z/2), (thick/2,box_y/2*expand, box_z/2*expand), material='sidewalls')
        rightWallId =  O.bodies.append(rightWall)
        wallIds.append(rightWallId)

        leftWall   = yade.utils.box(  ( box_off_x+-thick/2,1.5*box_y,box_z/2), (thick/2,box_y/2*expand, box_z/2*expand), material='sidewalls')
        leftWallId =  O.bodies.append(leftWall)
        wallIds.append(leftWallId)

        for i in ballIds:
            O.bodies[i].state.blockedDOFs = 'zXY'
            
            
    elif (periodic_bc==7):
        #normal 3D space, top & bottom walls have different material than side walls
        expand=1.5

        #create walls around the packing
        topWall   = yade.utils.box((box_x/2+box_off_x,2*box_y+box_off_y+thick/2,box_z/2), (box_x/2*expand,thick/2,box_z/2*expand), material='top_walls', wire=True)
        topWallId =  O.bodies.append(topWall)
        wallIds = [topWallId]

        botWall   = yade.utils.box((box_x/2+box_off_x,box_y+box_off_y-thick/2,box_z/2), (box_x/2*expand,thick/2, box_z/2*expand), material='bot_walls', wire=True)
        botWallId =  O.bodies.append(botWall)
        wallIds.append(botWallId)
                        
        rightWall   = yade.utils.box((box_x+thick/2+box_off_x,1.5*box_y+box_off_y,box_z/2), (thick/2,box_y/2*expand, box_z/2*expand), material='sidewalls', wire=True)
        rightWallId =  O.bodies.append(rightWall)
        wallIds.append(rightWallId)

        leftWall   = yade.utils.box((-thick/2+box_off_x,1.5*box_y+box_off_y,box_z/2), (thick/2,box_y/2*expand, box_z/2*expand), material='sidewalls', wire=True)
        leftWallId =  O.bodies.append(leftWall)
        wallIds.append(leftWallId)
        
        frontWall   = yade.utils.box((box_x/2+box_off_x,1.5*box_y+box_off_y,box_z+thick/2), (box_x/2*expand,box_y/2*expand, thick/2), material='sidewalls', wire=True)
        frontWallId =  O.bodies.append(frontWall)
        wallIds.append(frontWallId)
        
        backWall = yade.utils.box((box_x/2+box_off_x,1.5*box_y+box_off_y,-thick/2), (box_x/2*expand,box_y/2*expand, thick/2), material='sidewalls', wire=True)
        backWallId =  O.bodies.append(backWall)
        wallIds.append(backWallId)
        
        
        if angled_model > 0:
             axis = (0,0,1)
             center = Vector3(0,0,0)
             angl = angled_model*pi/180
             q = Quaternion(axis,angl)
        
             for f in wallIds:
                p = O.bodies[f].state.pos
                newP = center+q*(p-center)
                O.bodies[f].state.pos=O.bodies[f].state.refPos=newP
                O.bodies[f].state.ori = q
        
        
    elif (periodic_bc == 8):
        expand=1.5
            
                  # yade.utils.box(  center                             extents= half sizes!
        topWall   = yade.utils.box( ( box_off_x+box_x/2,2*box_y+thick/2,box_z/2), (box_x/2*expand,thick/2,box_z/2*expand), material='top_walls')
        topWallId =  O.bodies.append(topWall)
        wallIds = [topWallId]

        botWall   = yade.utils.box( ( box_off_x+box_x/2,box_y-thick/2,box_z/2), (box_x/2*expand,thick/2, box_z/2*expand), material='bot_walls')
        botWallId =  O.bodies.append(botWall)
        wallIds.append(botWallId)
            
#        for i in ballIds:
#                O.bodies[i].state.blockedDOFs = 'xzXYZ'


    return wallIds, topWallId, botWallId
##################################################################################################



if (not restart):
  topWall_list = []
  botWall_list = []
  wall_clump_list = [] 
  wall_listlist = []


  for i in range(n_windows_per_pocket):
          wallIds, topWallId, botWallId = wall_generate(i * 2.5 * box_x, i*2.5*box_y )
          wall_listlist.append(wallIds)
          topWall_list.append( topWallId)
          botWall_list.append( botWallId)        
          wallClump=O.bodies.clump(wallIds)
          O.bodies[wallClump].state.blockedDOFs = 'xyzXYZ'
          wall_clump_list.append( wallClump)
        

  if (pocketLoc2 > 0) and (pocketLoc2_Y == 0):
          for i in range(n_windows_per_pocket):
                  wallIds, topWallId, botWallId = wall_generate( (i+n_windows_per_pocket) * 2.5 * box_x, 0 )
                  wall_listlist.append(wallIds)
                  topWall_list.append( topWallId)
                  botWall_list.append( botWallId)        
                  wallClump=O.bodies.clump(wallIds)
                  O.bodies[wallClump].state.blockedDOFs = 'xyzXYZ'
                  wall_clump_list.append( wallClump)
                  
  if (pocketLoc2_Y > 0):
          for i in range(n_windows_per_pocket):
                  wallIds, topWallId, botWallId = wall_generate( 0, (i+n_windows_per_pocket) * 2.5 * box_y )
                  wall_listlist.append(wallIds)
                  topWall_list.append( topWallId)
                  botWall_list.append( botWallId)        
                  wallClump=O.bodies.clump(wallIds)
                  O.bodies[wallClump].state.blockedDOFs = 'xyzXYZ'
                  wall_clump_list.append( wallClump)

                  
                  
wall_state_list = [] 
wallF_accum_list = []
TopwallF_accum_list = []
for i in wall_clump_list:
        wallF_accum_list.append(0)
        TopwallF_accum_list.append(0)
        wall_state_list.append( O.bodies[i].state )
        

#####################################################################################################



def base_excitation_set_F_amp():
#note: treating F_amp as real, but its actually complex.  neglecting the phase component. probably okay.
#for mode 2 at the fully fused Q, angle changes by 0.003% over +/-2% of the nat frequency

#this method neglects the rigid body contribution of base velocity on the pocket velocity.  
        
        global F_amp, F_phase

        if (use_PLL):
                omega = PLL_Freq_rad
        else:
                omega = 2*math.pi*(f_sweep_start + f_sweep_rate*O.time)
        
        #F_amp for base excitation in acceleration (g's)
        F_amp = np.sqrt( (Modal_Mass * omega**2)**2 + (Modal_Damp *omega)**2)* base_excitation*9.81/(omega**2) * F_integral
        F_phase = 0
        

def base_excitation_set_F_amp_new():
#this method includes the rigid body contribution of base velocity on the pocket velocity.  
        global F_amp, F_phase, base_velocity_amp, base_displacement_amp

        raise Exception('pll code not updated')
        
        omega = 2*math.pi*(f_sweep_start + f_sweep_rate*O.time)        

        base_displacement_amp = (base_excitation*9.81) / omega**2  
        base_velocity_amp     = (base_excitation*9.81) / omega
        
        F_amp_complex = (-Modal_Mass + 1j * Modal_Damp / omega) * (base_excitation*9.81) * F_integral
                
        F_amp   = np.abs(F_amp_complex)
        F_phase = np.angle(F_amp_complex)
        

                
###########################################################################################
#PLL setup
###########################################################################################

def PLL():
        global zf_X, zf_Y, Amp, Phase, PhaseError, PhaseErrorAccum, PhaseErrorPrev,PhaseErrorPrevPrev, PLL_Freq_rad, PLL_prev_time

        #fixme, delta_time not being perfectly uniform will cause filter bandwidth to change
        delta_time = (O.time - PLL_prev_time)
        PLL_prev_time = O.time
        
        position=q1_state.pos[1] + base_displacement_amp * sin(theta)

        X = position * sin(theta)
        Y = position * cos(theta)

        Xf, zf_X = signal.lfilter( b,a, [X], zi=zf_X)
        Yf, zf_Y = signal.lfilter( b,a, [Y], zi=zf_Y)

        Amp   = 2*sqrt( Xf**2 + Yf**2) #prior to Aug 28th was off by factor of 2
        Phase = atan2( Yf, Xf)

        PhaseError = Phase - ( - math.pi/2)

        #1st order difference
#        PhaseErrorDer = (PhaseError - PhaseErrorPrev) / delta_time
#        PhaseErrorPrev = PhaseError

        #2nd order difference
        PhaseErrorDer = (3*PhaseError - 4*PhaseErrorPrev + PhaseErrorPrevPrev) / (2*delta_time)
        PhaseErrorPrevPrev = PhaseErrorPrev
        PhaseErrorPrev = PhaseError

        
        #give transients time to settle before adjusting freq
        if (O.time > PLL_deadtime):
                PhaseErrorAccum = PhaseErrorAccum + PhaseError * delta_time
                PLL_Freq_rad = f_sweep_start_rad + PhaseError * PLL_Pgain +  PhaseErrorAccum * PLL_Igain + PhaseErrorDer * PLL_Dgain

                if (PLL_Freq_rad < 0):
                        raise Exception('PLL went unstable')
 

if (use_PLL) and (f_sweep_rate != 0):
        raise Exception('PLL cannot be on during frequency sweep')


PLL_deadtime = 1.7 * PLL_Qguess * nat_period


if (use_PLL == 1):
        raise Exception('old buggy version of PLL gain calculation no longer supported')
elif (use_PLL > 1):
        PLL_filterBW = f_sweep_start / PLL_filter_BW_cycles

        Fs = 1 / plotDataVirtPeriod
        b,a = signal.butter(2, PLL_filterBW / (Fs/2) )

        PLL_Pgain, PLL_Igain, PLL_Dgain    = calc_PLL_gain(PLL_Qguess, f_sweep_start_rad, b,a,Fs )

        if (not restart):
                zf_X = numpy.zeros( max(len(a), len(b)) - 1 )
                zf_Y = numpy.zeros( max(len(a), len(b)) - 1 )
                PhaseErrorAccum = 0
                PhaseErrorPrev = 0
                PhaseErrorPrevPrev = 0
                PhaseError = 0
                PLL_Freq_rad = f_sweep_start_rad
                swing_coupled_extF_new_prev_time=0
                PLL_prev_time = 0


#####################################################################################################
#DEFINING ENGINES ####



#move these to the spreadsheet?
ball_ke_output= False #turn these off to reduce file size and slight increase to computation time
peak_force_output= False and (ball_gen_method > 1)
specific_ball_force_output = False
snapshot_enabled = False

if (damp_fix == 1):
    print('Fix for num damping is applied: dummy particle no longer has damping')
    O.bodies[q1_id].state.isDamped = False
#    for i in range(len(wall_listlist[0])):
#        O.bodies[wall_listlist[0][i]].state.isDamped = False
        
else:
    if ( damp > 0):
            print('CURRENTLY NOT USING THE NUMERICAL DAMPING BUG FIG (num_damp_fix = 1). ALL WALLS AND Q1_Dummy ARE UNDER EFFECTS OF NUMERICAL DAMPING')
        

if (snapshot_enabled):
        if (restart):
                snapshot_file_name = os.path.basename( sys.argv[1] ) + "_snapshot_" + str(linenumber) + "_restart.txt"
        else:
                snapshot_file_name = os.path.basename( sys.argv[1] ) + "_snapshot_" + str(linenumber) + ".txt"
                
        file2 = open(snapshot_file_name,"w+") 
        file2.write( str( len(ballIds)) + "," + str(len(botWall_list)) + "\n")
        snapshot_start = 0.6
        snapshot_stop  = 0.7        
        snapshot_interval = nat_period / 20

        if (not restart):
                snapshot = PyRunner(command='snapshot()',virtPeriod=snapshot_interval )
        else:
                #on restart, engines are already defined, so have to modify engine.
                #does depend on the position of the snapshot engine in the list!
                O.engines[-2].dead = False
                O.engines[-2].virtPeriod=snapshot_interval
                
else:
        snapshot = PyRunner(command='snapshot()', dead=True )


if (not restart):
  theta=0
  extF = 0
  base_vel = 0
  F_phase = 0
  base_displacement_amp=0
  last_output_iter = 0
  detect_factor = 1

  #right here, set int_loop based on some variable
  if (restitution > 0) or (cn > 0):
          if (normalCohesion != 0):
                  raise Exception('error, currently no way to combine viscoelasticity and cohesion')
          else:
                  print('using Law2_ScGeom_ViscElPhys_Basic()')
                  int_loop = InteractionLoop(
                          [Ig2_Sphere_Sphere_ScGeom(),Ig2_Box_Sphere_ScGeom()],
                          [Ip2_ViscElMat_ViscElMat_ViscElPhys()],
                          [Law2_ScGeom_ViscElPhys_Basic()] )
  else:
          if (normalCohesion > 0) or (Eroll !=0):
                  if (DMT_gamma > 0):
                          raise Exception('asked for attractive forces but not getting them')
                  
                  print('Using Law2_ScGeom6D_CohFrictPhys_CohesionMoment(useIncrementalForm = True)')
                  int_loop = InteractionLoop(
                          [Ig2_Sphere_Sphere_ScGeom6D(), Ig2_Box_Sphere_ScGeom6D() ],
                          [Ip2_CohFrictMat_CohFrictMat_CohFrictPhys(setCohesionOnNewContacts=True)],
                          [Law2_ScGeom6D_CohFrictPhys_CohesionMoment(useIncrementalForm = True)] )
                  #fixme, may need to think about time step here!
          elif (normalCohesion < 0):
                  print('Law2_ScGeom6D_CohFrictPhys_CohesionMoment_drk()')
                  #hidden flag. negative value for normal cohesion means to use the new model
                  int_loop = InteractionLoop(
                          [Ig2_Sphere_Sphere_ScGeom6D(), Ig2_Box_Sphere_ScGeom6D() ],
                          [Ip2_CohFrictMat_CohFrictMat_CohFrictPhys(setCohesionOnNewContacts=True)],
                          [Law2_ScGeom6D_CohFrictPhys_CohesionMoment_drk()] )
                  
          elif (rollfrict > 0):
                  print('using Law2_ScGeom_ViscElPhys_Basic()')
                  int_loop = InteractionLoop(
                          [Ig2_Sphere_Sphere_ScGeom(),Ig2_Box_Sphere_ScGeom()],
                          [Ip2_ViscElMat_ViscElMat_ViscElPhys()],
                          [Law2_ScGeom_ViscElPhys_Basic()] )  
                  
          elif (mindlin == 2):
              print('Law2_ScGeom_MindlinPhys_MindlinDeresiewitz()')
              int_loop = InteractionLoop(
                          [Ig2_Sphere_Sphere_ScGeom(),Ig2_Box_Sphere_ScGeom()],
                          [Ip2_FrictMat_FrictMat_MindlinPhys()],
                          [Law2_ScGeom_MindlinPhys_MindlinDeresiewitz()] )
              
          elif (mindlin == 1):
              print('Law2_ScGeom_MindlinPhys_Mindlin()')
              int_loop = InteractionLoop(
                          [Ig2_Sphere_Sphere_ScGeom(),Ig2_Box_Sphere_ScGeom()],
                          [Ip2_FrictMat_FrictMat_MindlinPhys(gamma=DMT_gamma)],
                          [Law2_ScGeom_MindlinPhys_Mindlin(includeAdhesion=(DMT_gamma>0))] )
          elif (mindlin == 3):
              print('using Law2_ScGeom_MindlinPhys_Mindlin() with hard coded en=0.5 in Ip2')
              int_loop = InteractionLoop(
                          [Ig2_Sphere_Sphere_ScGeom(),Ig2_Box_Sphere_ScGeom()],
                          [Ip2_FrictMat_FrictMat_MindlinPhys(en=0.5)],
                          [Law2_ScGeom_MindlinPhys_Mindlin()] )
          elif (DMT_gamma > 0):
                  print('using Law2_ScGeom_AdhesiveFrictPhys_CundallStrack()')
                  #bit of a hack, assuming the normal distribution 25 micron mean, 0.39*mean std dev, 2 std dev cutoff
                  #this sets the interaction detection factor so that the smallest particle exactly hits. 
                  min_dia = (1-2*0.39) * particle_dia
                  detect_factor = 1+DMT_cutoff / (min_dia/2)                                    
                  print("detect_factor " + str(detect_factor))
                  int_loop = InteractionLoop(
                          [Ig2_Sphere_Sphere_ScGeom(interactionDetectionFactor=detect_factor),Ig2_Box_Sphere_ScGeom(interactionDetectionFactor=detect_factor)],
                          [Ip2_FrictMat_FrictMat_AdhesiveFrictPhys(gamma=DMT_gamma, cutoff=DMT_cutoff)],
                          [Law2_ScGeom_AdhesiveFrictPhys_CundallStrack()] )                  
          else:
                  print('using Law2_ScGeom_FrictPhys_CundallStrack()')
                  int_loop = InteractionLoop(
                          [Ig2_Sphere_Sphere_ScGeom(),Ig2_Box_Sphere_ScGeom()],
                          [Ip2_FrictMat_FrictMat_FrictPhys()],
                          [Law2_ScGeom_FrictPhys_CundallStrack()] )

                             
                  

if (ball_gen_method == 2):
        file1.write( "iter,time,freq,theta,position, misc\n" )
        file1.write( "0,0," + str(f_sweep_start) + ",0,0,0\n" )
else:
        file1.write( "iter,time,freq,theta,wall_position, q1")
        for i in range(len(wallF_accum_list)):
                file1.write(",wallF" + str(i+1) )                
        for i in range(len(TopwallF_accum_list)):
                file1.write(",TopWallF" + str(i+1) )

                
        file1.write(",ball_cg_corrected, ball_centroid" )
#        file1.write(",top_wall_force,bot_wall_force") #For checking parameters for gap finding

        if (use_PLL):
                file1.write(",PLL_amp, PLL_phase" )
        
        if (ball_ke_output):
                file1.write(",q1_ke,ball_ke" )
        if (peak_force_output):
                temp_walllist = []
                for sublist in wall_listlist:
                        for item in sublist:
                                temp_walllist.append(item)

                file1.write(",mean_force, peak_force" )
                if (specific_ball_force_output):
                    file1.write(",ball1_r, ball2_r, int1_kn, int1_norm, int1_pen, ball3_r, ball4_r, int2_kn, int2_norm, int2_pen" )
                        
        file1.write("\n")


        file1.flush()




###########################################################################################
             
if (not restart):                
        if (force_method == 2):
                #for a sweep, this gets called every ten oscillation cycles, but not until cycle ten, so this sets the
                #initial condition. for fixed frequency, this is the only call ever.  
                base_excitation_set_F_amp()
        elif (force_method == 3):
                #for a sweep, this gets called every ten oscillation cycles, but not until cycle ten, so this sets the
                #initial condition. for fixed frequency, this is the only call ever.  
                base_excitation_set_F_amp_new()
                print("F_amp = " + str(F_amp),", F_phase = " + str(F_phase) )
                
        if ( f_sweep_start == f_sweep_stop):
                stopTime = 12 * q_factor / Cant_nat_freq_hz

                be = PyRunner(command="", dead=True,iterPeriod=1)
                if (f_sweep_rate != 0):
                        raise Exception('you requested f_sweep_start = f_sweep_stop, but f_sweep_rate was not zero.')                                        
        else:
                #sweep
                if ( (f_sweep_stop < f_sweep_stop) and (f_sweep_rate >= 0)):
                        raise Exception('you requested f_sweep_stop < f_sweep_start, but f_sweep_rate is positive.')
                        
                if ( (f_sweep_stop > f_sweep_stop) and (f_sweep_rate <= 0)):
                        raise Exception('you requested f_sweep_stop > f_sweep_start, but f_sweep_rate is negative.')

                        
                stopTime = np.abs(f_sweep_stop - f_sweep_start) / f_sweep_rate
                
                if (force_method== 2):
                        #to save time, we only need to call this about once per cantilever oscillation. More frequently than that makes
                        #no difference as the cantilever can't respond that fast anyway                        
                        #be = PyRunner(command='base_excitation_set_F_amp()',iterPeriod=1)
                        be = PyRunner(command='base_excitation_set_F_amp()',virtPeriod=10/Cant_nat_freq_hz )
                if (force_method == 3):
                        be = PyRunner(command='base_excitation_set_F_amp_new()',virtPeriod=10/Cant_nat_freq_hz )
                else:
                        be = PyRunner(command="", dead=True,iterPeriod=1)

        #stock yade way
                
        if (damp_gravity):
                #this is the default
                print("\nUsing default damp_gravity. This is not using the custom engine, where numerical damping is removed from gravity\n")
                newton=NewtonIntegrator(damping=damp, gravity=(0,grav_y,0))
        else:
                #right now only works with a custom compiled engine.
                print("\nUsing custom engine, where numerical damping is removed from gravity\n")

                newton=NewtonIntegrator(damping=damp, gravity=(0,grav_y,0), dampGravity=False )

        if (solution_method == 1):
                swing_engines = [ PyRunner(command='swing_coupled_extF_new()',virtPeriod=extFUpdateVirtPeriod), PyRunner(command='swing_coupled()',iterPeriod=1)]
        elif (solution_method == 2):
                raise Exception('solution method 2 was deleted. had not been maintained, was not use anymore')
        if (solution_method == 3):
                #all pockets have same motion, an enforced sine wave at a specified acceleration.  
                swing_engines = [ PyRunner(command='swing_decoupled()',iterPeriod=1)]
                

        if (StokesDrag):
                drag = LinearDragEngine( ids=ballIds, nu = 1.8e-5)
        else:
                drag = LinearDragEngine( ids=ballIds, nu = 1.8e-5, dead=True)

        if (use_PLL):
                PLLe = PyRunner(command='PLL()',virtPeriod=plotDataVirtPeriod)
        else:
                PLLe = PyRunner(command='PLL()',virtPeriod=plotDataVirtPeriod, dead=True)
                
        O.engines=[
                ForceResetter(),
                InsertionSortCollider([Bo1_Sphere_Aabb(aabbEnlargeFactor=detect_factor),Bo1_Box_Aabb()], allowBiggerThanPeriod=True, verletDist=verlet_dist ),
                int_loop,
                timestepper,
                be, PLLe] + swing_engines + [drag, newton,                
                PyRunner(command='addPlotData_coupled()',virtPeriod=plotDataVirtPeriod),
                PyRunner(command='check_stop_time()',iterPeriod=1000),
                snapshot ,
                PyRunner(command='checkpoint()', realPeriod=3600),  
        ]


stopTime = 5 # bad hack delete me
        
#this was needed due to bug in O.run, reported on forums and think fixed by now, but haven't recompiled the version on the cluster yet        
def check_stop_time():
        global SnapCount
        if (O.time > 0) and (SnapCount == 0):
            export.text("BallPos_"+str(linenumber)+"_"+str(SnapCount))
            SnapCount+=1
        elif (O.time > 1) and (SnapCount == 1):
            export.text("BallPos_"+str(linenumber)+"_"+str(SnapCount))
            SnapCount+=1
        elif (O.time > 2) and (SnapCount == 2):
            export.text("BallPos_"+str(linenumber)+"_"+str(SnapCount))
            SnapCount+=1
        elif (O.time > 2.6) and (SnapCount == 3):
            export.text("BallPos_"+str(linenumber)+"_2.6")
            SnapCount+=1            
        elif (O.time > 3) and (SnapCount == 4):
            export.text("BallPos_"+str(linenumber)+"_3")
            SnapCount+=1
        elif (O.time > 3.5) and (SnapCount == 5):
            export.text("BallPos_"+str(linenumber)+"_3.5")
            SnapCount+=1
        elif (O.time > 4) and (SnapCount == 6):
            export.text("BallPos_"+str(linenumber)+"_4")
            SnapCount+=1
        elif (O.time > stopTime):
                export.text("BallPos_"+str(linenumber))
                O.pause()

#only call these every 10 steps or so to save time, modal parameters don't change very often
#merging the old routine in here, slight run time penalty for easier management
def swing_coupled_extF_new():
        global theta, extF, base_vel, swing_coupled_extF_new_prev_time

        
        if (use_PLL):
                #rectangular integration, could be more precise, but time steps are small and this was good enough for VEDA
                #just be careful that we might not be called at a perfectly uniform rate!
               theta = theta + PLL_Freq_rad * (O.time - swing_coupled_extF_new_prev_time)
               swing_coupled_extF_new_prev_time = O.time
        elif (f_sweep_rate == 0):
                theta = 2*math.pi * ( f_sweep_start * O.time)
        else:
                theta = 2*math.pi * ( f_sweep_start * O.time + 0.5*f_sweep_rate * O.time**2)

                
        #base displacement is sin(theta), so force is  sin(theta+phase), velocity = cos(theta)
        extF     = F_amp * math.sin(theta + F_phase)

        if (force_method == 2):
                base_vel = 0
        else:
                base_vel = base_velocity_amp *  math.cos(theta )


                
def swing_coupled():
        global wallF_accum_list, TopwallF_accum_list

        #this is right and decently fast
        kF   = -Modal_Stiff *  q1_state.pos[1] 
        cF   = -Modal_Damp  *  q1_state.vel[1]         
        
        #this is the way I'd like to do it, but forces on clumps only get updated INSIDE the NewtonIntegrator
        #by that time it's too late.  So have to do it the long way.  
        #wallF = O.forces.f(wallClump)[1] * balls_mass_ratio

        wallF = 0

#        print(base_vel)
        
        for i in range( len( wall_state_list)):  #loop overs the windows
                wall_state_list[i].vel[1] = (base_vel + q1_state.vel[1] * Psi_pocket_list[i]) 
                if sidewall_frict_on == 0:
                    #no friction with the front/back/left/right walls, therefore there can never be any
                    #vertical force on them. So to save time, don't even check.
                    wallF1 = ( O.forces.f(topWall_list[i] )[1] +  O.forces.f(botWall_list[i] )[1] ) * balls_mass_ratio * Psi_pocket_list[i]
                else:
                    wallF1 = ( O.forces.f(wall_listlist[i][0] )[1] +  O.forces.f(wall_listlist[i][1] )[1]  + O.forces.f(wall_listlist[i][2] )[1] + O.forces.f(wall_listlist[i][3] )[1]  + O.forces.f(wall_listlist[i][4] )[1]  + O.forces.f(wall_listlist[i][5] )[1] ) * balls_mass_ratio * Psi_pocket_list[i]

                wallF = wallF + wallF1

                wallF_accum_list[i]     = wallF_accum_list[i] + wallF1  #average the force over each output period. technically should be a weighted average as dt is not same every time, but this is close enough.  
                TopwallF_accum_list[i] = TopwallF_accum_list[i] + O.forces.f(topWall_list[i] )[1] * balls_mass_ratio #top wall only
                
        O.forces.addF(q1_id , ( 0, extF + kF + wallF + cF - grav_y * Modal_Mass , 0) )




def swing_decoupled():
        global wallF_accum_list, TopwallF_accum_list

        if (use_PLL):
                raise Exception('not implemented')
        elif (f_sweep_rate == 0):
                omega = 2*math.pi * f_sweep_start
                theta = omega * O.time
        else:
                theta = 2*math.pi * ( f_sweep_start * O.time + 0.5*f_sweep_rate * O.time**2)
                omega = 2*math.pi * ( f_sweep_start + O.time * f_sweep_rate)

        wallF = 0
        
        for i in range( len( wall_state_list)):  #loop overs the windows
                wall_state_list[i].vel[1] =  q1_state.vel[1]
                if sidewall_frict_on == 0:
                    #no friction with the front/back/left/right walls, therefore there can never be any
                    #vertical force on them. So to save time, don't even check.
                    wallF1 = ( O.forces.f(topWall_list[i] )[1] +  O.forces.f(botWall_list[i] )[1] ) * balls_mass_ratio
                else:
                    wallF1 = ( O.forces.f(wall_listlist[i][0] )[1] +  O.forces.f(wall_listlist[i][1] )[1]  + O.forces.f(wall_listlist[i][2] )[1] + O.forces.f(wall_listlist[i][3] )[1]  + O.forces.f(wall_listlist[i][4] )[1]  + O.forces.f(wall_listlist[i][5] )[1] ) * balls_mass_ratio

                wallF = wallF + wallF1

                wallF_accum_list[i]     = wallF_accum_list[i] + wallF1  #average the force over each output period. technically should be a weighted average as dt is not same every time, but this is close enough.  
                TopwallF_accum_list[i] = TopwallF_accum_list[i] + O.forces.f(topWall_list[i] )[1] * balls_mass_ratio #top wall only

                
        #fixme, here we impose a velocity on q1.  This could be done quicker with a HarmonicMotionEngine, but CPU time is not the most important thing now. just get it to work.
        q1_state.vel[1] = (base_excitation * 9.81) / ( omega) * cos( theta)
        



def addPlotData_coupled():
        global wallF_accum_list, TopwallF_accum_list, last_output_iter, theta

        #fixme, need to localize excitation frequency information in one place
        if (use_PLL):
                f = PLL_Freq_rad / (2 * np.pi)
        elif (f_sweep_rate != 0):
                f = f_sweep_start + O.time * f_sweep_rate                
        else:
                f = f_sweep_start


        position=q1_state.pos[1] + base_displacement_amp * sin(theta)

        
        #for averaging output force
        delta_iterations = O.iter - last_output_iter       
        last_output_iter = O.iter
                
        if (ball_gen_method == 2):
                position2=O.bodies[ballIds[0]].state.pos[1]
                file1.write("0," +  str(O.time) + "," + str(f) + "," + str(theta) + "," +  str( position) + ","+  str( position2) + "\n" )
        else:
                file1.write( str(O.iter) + "," +  str(O.time) + "," + str(f) + "," + str(theta) + "," +  "{:.7g}".format(position) + "," +  "{:.7g}".format(q1_state.pos[1]) )
                for i in range(len(wallF_accum_list)):
                        file1.write( "," + "{:.5g}".format( wallF_accum_list[i] /delta_iterations))
                        wallF_accum_list[i] =0 
                for i in range(len(TopwallF_accum_list)):
                        file1.write( "," + "{:.5g}".format( TopwallF_accum_list[i] /delta_iterations))
                        TopwallF_accum_list[i] =0 

                        
                ball_cg       = 0
                ball_centroid = 0 #this is what we formerly called c.g.
                for i in ballIds:
                        ball_centroid = ball_centroid + O.bodies[i].state.pos[1] / len(ballIds)
                        ball_cg       = ball_cg       + O.bodies[i].state.pos[1] * O.bodies[i].state.mass / ballsmass_all_pockets

                file1.write("," + "{:.6g}".format(ball_cg- 1.5 * box_y)+ "," + "{:.9g}".format(ball_centroid - 1.5 * box_y) )
#                file1.write(","+ "{:.6g}".format(O.forces.f(wall_listlist[0][0])[1]) + "," + "{:.6g}".format(O.forces.f(wall_listlist[0][1])[1])) #TEMP

                if (use_PLL):
                        file1.write("," + "{:.6g}".format(Amp)+ "," + "{:.6g}".format(Phase) )
                                        
                if (ball_ke_output):
                        q1_ke  = 0.5 * Modal_Mass *( q1_state.vel[1]**2 )
                                                
                        all_ke = yade._utils.kineticEnergy() #this includes q1. and the walls, but walls *should* have trivial density
                        ball_ke = (all_ke - q1_ke) * balls_mass_ratio

                        # ball_ke = 0
                        # for i in ballIds:
                        #         #fixme add rotational
                        #         ball_ke = ball_ke + O.bodies[i].state.mass * ( O.bodies[i].state.vel[0]**2 +O.bodies[i].state.vel[1]**2 +O.bodies[i].state.vel[2] **2 )
                        
                        file1.write( "," + "{:.5g}".format(q1_ke) +  "," + "{:.5g}".format(ball_ke))
                if (peak_force_output):
                        #this is going to alias since we are only outputting every 100 iterations or so,
                        #but doing it otherwise will be too computationally expensive
                        mean_force = 0
                        peak_force = 0

                        part_list = []
                        part_list_other = []
                        force_list = []
                        temp_list = []
                                
                                
#                        cnt = 0
                        for i in (O.interactions):
                                force = i.phys.normalForce.norm()
                                mean_force = mean_force + force
                                peak_force = max( force, peak_force)
                                
                                if (i.id1 not in temp_walllist) and (i.id2 not in temp_walllist) and (i.id1 not in wall_clump_list) and (i.id2 not in wall_clump_list):
                                    part_list.append(i.id1)
                                    part_list_other.append(i.id2)
                                    force_list.append(i.phys.normalForce.norm())
#                                    cnt+=1
                                    
                        
                        
                        temp_list = random.sample(list(enumerate(part_list)),4)
                        
                        mean_force = mean_force / len(O.interactions) #THIS SHOULD BE OUTSIDE LOOP
                        file1.write( "," + "{:.5g}".format(mean_force) +  "," + "{:.5g}".format(peak_force))
                        if (specific_ball_force_output):
#                            file1.write( "," + "{:.4g}".format(force_list[temp_list[0][0]]) +  "," + "{:.4g}".format(O.bodies[temp_list[0][1]].shape.radius) +  "," + "{:.4g}".format(force_list[temp_list[1][0]]) +  "," + "{:.4g}".format(O.bodies[temp_list[1][1]].shape.radius) +  "," + "{:.4g}".format(force_list[temp_list[2][0]]) +  "," + "{:.4g}".format(O.bodies[temp_list[2][1]].shape.radius) +  "," + "{:.4g}".format(force_list[temp_list[3][0]]) +  "," + "{:.4g}".format(O.bodies[temp_list[3][1]].shape.radius))
                            file1.write( "," + "{:.4g}".format(O.bodies[part_list[temp_list[0][0]]].shape.radius) +  \
                                        "," + "{:.4g}".format(O.bodies[part_list_other[temp_list[0][0]]].shape.radius) +  \
                                        "," + "{:.4g}".format(O.interactions[part_list[temp_list[0][0]],part_list_other[temp_list[0][0]]].phys.kn) +  \
                                        "," + "{:.4g}".format(O.interactions[part_list[temp_list[0][0]],part_list_other[temp_list[0][0]]].phys.normalForce.norm()) +  \
                                        "," + "{:.4g}".format(O.interactions[part_list[temp_list[0][0]],part_list_other[temp_list[0][0]]].geom.penetrationDepth) +  \
                                        "," + "{:.4g}".format(O.bodies[part_list[temp_list[1][0]]].shape.radius) +  \
                                        "," + "{:.4g}".format(O.bodies[part_list_other[temp_list[1][0]]].shape.radius) +  \
                                        "," + "{:.4g}".format(O.interactions[part_list[temp_list[1][0]],part_list_other[temp_list[1][0]]].phys.kn) +  \
                                        "," + "{:.4g}".format(O.interactions[part_list[temp_list[1][0]],part_list_other[temp_list[1][0]]].phys.normalForce.norm()) +  \
                                        "," + "{:.4g}".format(O.interactions[part_list[temp_list[1][0]],part_list_other[temp_list[1][0]]].geom.penetrationDepth))

                file1.write("\n")
                        
                

        if (not batch):
                if (ball_gen_method == 2):
                        plot.addData(freq=position2 , position=O.forces.f(ballIds[0])[1] )
                else:
#                        plot.addData(t1=f , PLL_Freq=PLL_Freq_rad / (2*math.pi), t2 =f, pos=position, t3=f, Phase=PhaseError * 180 / math.pi,  t4=f, theta=theta )
                        plot.addData(freq=f , position=np.abs(position))
                        
#                plot.addData(freq=f , wallFR=O.forces.f(wall_listlist[0][2])[1],wallFL=O.forces.f(wall_listlist[0][3])[1],wallFT=O.forces.f(wall_listlist[0][0])[1])



# define how to plot data: 'i' (step number) on the x-axis, raw position and order tracked magnitude on y
# show the plot on the screen, and update while the simulation runs
if ((not batch) and (solution_method == 1)) :
#        plot.plots={'t1':('PLL_Freq'), 't2':('pos'), 't3':('Phase'), 't4':('theta'), }
#        plot.plots={'freq':('wallFR','wallFL','wallFT')}
        plot.plots={'f':('position') }

        plot.plot()

def snapshot():
        if (snapshot_enabled) and (O.time > snapshot_start) and (O.time < snapshot_stop):
                file2.write(str( O.time))
                for i in ballIds:
                        file2.write("," + "{:.5g}".format(O.bodies[i].state.pos[0]) + "," + "{:.5g}".format(O.bodies[i].state.pos[1]) + "," + "{:.5g}".format(O.bodies[i].state.pos[2]))
                for i in botWall_list:
                        file2.write("," + "{:.6g}".format(O.bodies[i].state.pos[1]))                        
                file2.write("\n")
                
                
if (batch):
        if (timing):
                O.timingEnabled = True
                O.run(10000, True)
                yade.timing.stats()
        else:
                O.run(-1, True)  #this will run until the O.pause hits in addPlotData
                
        file1.close()
        if (snapshot_enabled):
                file2.close()
        O.exitNoBacktrace()



        
        
