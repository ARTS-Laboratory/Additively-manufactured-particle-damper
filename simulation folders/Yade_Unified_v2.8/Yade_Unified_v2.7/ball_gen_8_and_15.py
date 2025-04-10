#yadedaily ball_gen_8.py [FILENAME] [LINENUMBER] [INDENT]

#ball_gen_method 8 and 15 for creating particle by random insertion and then letting them settle.
#
#
#  This will take the exact same input format as the other code, and generate the packing
#The final ball locations are stored in a file with the name based on the size of the box and the size
#of the particle.  for example: grav_dep_pack_250.0_1300.0_250.0_uniform_25.0.csv (numbers are in
#micrometers), When yade_sdoftrial9.py runs, if
#ball_gen_method==8. it looks for a suitable CSV file and then loads it.
#If no such file then error message
#

import yade
from yade import pack
import math #import mathematical formulae
import numpy as np
import yade.timing
import csv
import os.path

import sys
sys.path.append(".") # this is not necessary in a generic python shell, but inside YADE it is for some odd reason
from common_functions import wall_generate, ball_gen_8_output_filename

import matplotlib.pyplot as pyplot

#######################################################

#goal,find the packing factor at which the force is zero, but if we added just one more particle, there would be an initial compression


#################################################
##### process command line arguments
#################################################

# SAMUEL MODIFIED CODE
# Removed Batch feature completely in favor of Indent bool as we never use the Batch feature

if ( len(sys.argv) <= 3):
        batch = False
        indent = False
else:
        # batch = True
        batch = False
        indent = True

if ( len(sys.argv) <= 2):
        linenumber = 0
else:
        linenumber = int( sys.argv[2])-1

###############################################3
##### Read data from files (or use default parameters)
#################################################        

if (len(sys.argv) > 1):
        with open(sys.argv[1] ) as csvfile:
                reader = csv.DictReader(csvfile, skipinitialspace=True)
                row_list = list(reader)
                box_x = float( row_list[linenumber]['box_x'])
                box_y = float( row_list[linenumber]['box_y']) #this is the direction of excitation
                box_z = float( row_list[linenumber]['box_z'])
                f_sweep_start = float( row_list[linenumber]['f_sweep_start'])
                f_sweep_stop  = float( row_list[linenumber]['f_sweep_stop'])
                f_sweep_rate  = float( row_list[linenumber]['f_sweep_rate'])

                try:
                        bc = int(row_list[linenumber]['bc'])
                except:
                        bc = int(row_list[linenumber]['periodic_bc']) #backwards compatible


                particle_dia = float( row_list[linenumber]['particle_dia'])
                try:
                        particle_dia2 = float( row_list[linenumber]['particle_dia2'])
                except:
                        #backwards compatibility
                        particle_dia2 = particle_dia
                        
                F_amp = float( row_list[linenumber]['force_amp'])
                try:
                        grav_x = float( row_list[linenumber]['grav_x'])
                except:
                        grav_x = 0
                                
                grav_y = float( row_list[linenumber]['grav_y'])
#                damp = float( row_list[linenumber]['num_damping']) # hard coded for this file

                damp_gravity = int(row_list[linenumber]['damp_gravity'])
                
                dt_in = float( row_list[linenumber]['time_step'])
#                verlet_dist = float( row_list[linenumber]['verlet_dist'])
                ball_gen_method = int(row_list[linenumber]['ball_gen_method'])
                
                q_factor = float(row_list[linenumber]['q_Factor'])
                mode_num = int(row_list[linenumber]['mode_num'])

                young = float(row_list[linenumber]['young_mod'])
                
#                solution_method = int(row_list[linenumber]['solution_method'])
#                force_method = int(row_list[linenumber]['force_method']) # 1 = force, 2 = base excitation
                base_excitation = float(row_list[linenumber]['base_excitation'])

                distribution = row_list[linenumber]['distribution']
                                
                target_packing_factor = float(row_list[linenumber]['packing_factor_bg8'])
                seedType = int(row_list[linenumber]['seed'])

                try:
                        #for U tennesee system
                        NumPocketsForArbitrary = int(row_list[linenumber]['NumPocketsForArbitrary'])
                except:
                        NumPocketsForArbitrary = nan
                        
#################################################
##### DEFINING VARIABLES AND MATERIALS ######
#################################################

# SAMUEL ADDED VARIABLES
# ----------------------
#deleteCount = 0
# ----------------------

if (not batch):
        from yade import qt, plot
        qt.View() #open the controlling and visualization interfaces

if (ball_gen_method != 8) and (ball_gen_method != 15) :
        raise Exception("input file does not ask for ball gen method 8. nothing to do")

if (ball_gen_method == 15):
        #hack in the naming convention
        box_x = box_x * NumPocketsForArbitrary

finalFricDegree = 0   #this lets particles find a home easier


young_sp=young
young_wl=young
        
rho    = 8230 #generic inco 718, kg/m^3


mn = Vector3(0,      box_y,      0)
mx = Vector3(box_x,2*box_y,  box_z) #corners of the initial packing
thick = 2*particle_dia # the thickness of the walls

std_dev_dia = 0.39 * particle_dia #fixme make an input parameter        
###################


out_file_name = ball_gen_8_output_filename( target_packing_factor, particle_dia, box_x, box_y, box_z, seedType, bc, distribution, damp_gravity, std_dev_dia)

###################        
                
status_file_name = os.path.basename( sys.argv[1] ) + "_gravdep_" + str(linenumber) + ".txt"

file1 = open(status_file_name,"w+")
file1.write('iter, time, KE, KE_maxKE, topWallForce, botWallForce\n')

###########################################################################################
#GENERATING BODIES #####
###########################################################################################


maxKE = 0
currentKE = 0
state_change_time = 0
state = 0

        
def checkUnbalanced7():        
        global maxKE, currentKE, prevKE

        currentKE = yade._utils.kineticEnergy()

        if (currentKE > maxKE):
                maxKE = currentKE

        outstr = str(O.iter) + "," + str(O.time) + "," + str(currentKE) + "," + str(currentKE / maxKE) + "," + str(O.forces.f(topWallId)[1] ) + "," + str(O.forces.f(botWallId)[1] )
        file1.write(outstr + "\n")
        print(outstr )

        #after initial movement we can decrease the damping a bit and open up the timestep
        #fixme should be a state machine
        if (O.iter > 10000) and (O.time > 0.005) :
                my_newton.damping = 0.2
                O.dt=.8*PWaveTimeStep()
                        
        if (O.iter > 20000) and (O.time > 0.050)  and (currentKE <= prevKE) and (currentKE < 0.0000002 * maxKE) :
                print('finished')
                checkpoint(ballIds)
                O.pause()
                
        prevKE = currentKE


def checkUnbalancedBig():        
        global maxKE, currentKE, prevKE

        currentKE = yade._utils.kineticEnergy()

        if (currentKE > maxKE):
                maxKE = currentKE

        outstr = str(O.iter) + "," + str(O.time) + "," + str(currentKE) + "," + str(currentKE / maxKE) + "," + str(O.forces.f(topWallId)[1] ) + "," + str(O.forces.f(botWallId)[1] )
        file1.write(outstr + "\n")
        print(outstr )
                
        #after initial movement we can decrease the damping a bit and open up the timestep
        #fixme this needs to be re-written.  last case will never get reached!
        if  (O.time > 0.0005) and (O.time < 0.001) :
                collider.verletDist = -0.02
                my_newton.damping = 0.3
                O.dt=.6*PWaveTimeStep()
        elif  (O.time > 0.001) and (O.time < 0.002) :
                collider.verletDist = -0.04
                my_newton.damping = 0.2
                O.dt=.7*PWaveTimeStep()
        elif  (O.time > 0.002)  and (O.time < 0.005) :
                collider.verletDist = -0.06
                my_newton.damping = 0.1
                O.dt=.8*PWaveTimeStep()
        elif  (O.time > 0.005)  :
                collider.verletDist = -0.1
                my_newton.damping = 0.1
                O.dt=.9*PWaveTimeStep()                        
        elif  (O.time > 0.30)  and (currentKE <= prevKE) and (currentKE < 0.0000002 * maxKE) :
                print('finished')
                checkpoint(ballIds)
                O.pause()
                
        prevKE = currentKE

        
def kill_pos_velocity():
        for b in ballIds:
                if (O.bodies[b].state.vel[0] > 0):
                        O.bodies[b].state.vel[0]=0

def kill_velocity():
        for b in ballIds:                
                O.bodies[b].state.vel[0]=0
                O.bodies[b].state.vel[1]=0
                O.bodies[b].state.vel[2]=0


#for long running simulations, periodically dump out the results so we can monitor (and include velocity so we can potentially restart, velocity will be ignored by main code)
def checkpoint(outList):
        print("starting checkpoint")
        
        file2 = open(out_file_name, 'w+')
        file2.write('rad, x, y, z, vx, vy, vz\n');
        for i in outList:                
                file2.write(str(O.bodies[i].shape.radius ) + "," + str(O.bodies[i].state.pos[0]) + "," +  str(O.bodies[i].state.pos[1]) + "," + str(O.bodies[i].state.pos[2]) + "," + str(O.bodies[i].state.vel[0]) + "," +  str(O.bodies[i].state.vel[1]) + "," + str(O.bodies[i].state.vel[2]) + "\n" )
        file2.close()

                

#gets rid of escaped balls. not called automatically for this code, only manually.
def escaper_check():
        max_wall = [-9e99, -9e99, -9e99];
        min_wall = [ 9e99,  9e99,  9e99];

        for w in wallIds:
                max_wall = np.maximum( max_wall, O.bodies[w].state.pos)
                min_wall = np.minimum( min_wall, O.bodies[w].state.pos)

        for b in ballIds:
                for i in range(3):
                        if ((O.bodies[b].state.pos[i] > max_wall[i]) or (O.bodies[b].state.pos[i] < min_wall[i])):
                                print('Ball ' + str(b) + ' escaped! Deleting it now')
                                print('min wall: ' + str(min_wall) )
                                print('max wall: ' + str(max_wall) )
                                print('state   : ' + str(O.bodies[b].state.pos) )

                                O.bodies.erase(b)
                                ballIds.remove(b)
                                #deleteCount += 1
                                break

# SAMUEL ADDED CODE
# This is to delete balls that generate within the sphere attached to a wall.
# It does so by deleting particles within a box generated by using the sphere's radius.
# Written before understanding escaper_check() fully, not necessary for code as it is redundant.

def escaper_check_sam():
        max_wall = np.array([])
        min_wall = np.array([])
        wallSphereRad = thick*3

        #print('sphere state pos: '+str(O.bodies[wallSphereId].state.pos))
        for i in range(3):
                max_wall = np.append(max_wall, O.bodies[wallSphereId].state.pos[i] + wallSphereRad)
                min_wall = np.append(min_wall, O.bodies[wallSphereId].state.pos[i] - wallSphereRad)

        # print('max_wall: '+str(max_wall))
        # print('min_wall: '+str(min_wall))
        # print('rad: '+str(wallSphereRad))

        for b in ballIds:
                #for i in range(3):
                        if (O.bodies[b].state.pos[0] > (box_x + thick/2)) or (O.bodies[b].state.pos[0] < (0 - thick/2)):
                                print('Samuel\'s escaper check triggered, x')
                                # print('Ball ' + str(b) + ' escaped! Deleting it now')
                                # print('min wall: ' + str(min_wall) )
                                # print('max wall: ' + str(max_wall) )
                                # print('state   : ' + str(O.bodies[b].state.pos))
                                #deleteCount += 1
                                O.bodies.erase(b)
                                ballIds.remove(b)
                        elif (O.bodies[b].state.pos[1] < (box_y - thick/2)) or (O.bodies[b].state.pos[1] > (box_y*2 + thick/2)):
                                print('Samuel\'s escaper check triggered, y')
                                # print('Ball ' + str(b) + ' escaped! Deleting it now')
                                # print('box_y: ' + str(box_y) )
                                # print('state   : ' + str(O.bodies[b].state.pos))
                                #deleteCount += 1
                                O.bodies.erase(b)
                                ballIds.remove(b)
                        elif (O.bodies[b].state.pos[2] < -thick) or (O.bodies[b].state.pos[2] > box_z+thick):
                                print('Samuel\'s escaper check triggered, z')
                                #print('Ball ' + str(b) + ' escaped! Deleting it now')
                                #print('box_y: ' + str(box_z) )
                                #print('state   : ' + str(O.bodies[b].state.pos))
                                #deleteCount += 1
                                O.bodies.erase(b)
                                ballIds.remove(b)

#-----------------

# SAMUEL ADDED CODE
# sets a KE limit to pause the program. Designed to achieve parity between different iterations of the code.
# --------------------
def max_KE():
        if np.log10(currentKE) < -11:
                escaper_check()
                O.pause()
                print('KE threshold reached: '+str(np.log10(currentKE)))
                balls_deleted()
        else:
                print('Current KE: '+str(np.log10(currentKE)))
# --------------------

# SAMUEL ADDED CODE
# prints the particle count, used to keep track of particle deletions
# --------------------
def balls_deleted():
        print('\n')
        print('Current Particle Count: '+str(len(ballIds)))
        print('\n')
# --------------------

random.seed(a=seedType)
        
O.materials.append(FrictMat(young=young  ,poisson=0.3,frictionAngle=finalFricDegree,density=rho,label='spheres'))
        
packing_factor = 0
        
area_accumulator = 0
volu_accumulator = 0
        
ballIds = []        
        
while (packing_factor < target_packing_factor):
        if (distribution == 'uniform'):
                if (particle_dia2 == particle_dia):
                        rad = particle_dia / 2
                else:
                        rad = random.uniform( particle_dia/2,particle_dia2/2)
                                
        elif (distribution == 'normal'):                        
                rad_n = random.normalvariate(0, 1) #std normal
                if (rad_n < -2):
                        rad_n = -2
                if (rad_n > 2):
                        rad_n = 2

                rad = (particle_dia + std_dev_dia * rad_n) / 2

        #make sure it starts inside the box. don't want any escapers
        x = random.uniform( mn[0]+rad, mx[0]-rad)
        y = random.uniform( mn[1]+rad, mx[1]-rad)
        if (bc == 2):
                z = box_z/2
        else:
                z = random.uniform( mn[2]+rad, mx[2]-rad)

        
        id = O.bodies.append(sphere( (x,y,z) ,rad,material='spheres',color=(0,0,1)))
        ballIds.append(id)

        if (bc == 2):
                area_accumulator =  area_accumulator + math.pi * (rad**2)
                packing_factor = area_accumulator  / (box_x * box_y)
        else:
                volu_accumulator =  volu_accumulator + (4/3) * math.pi * (rad**3)                        
                packing_factor =  volu_accumulator  / (box_x * box_y * box_z)

                
if (bc == 2):
        #necessary for 2D to make sense
        for i in ballIds:
                O.bodies[i].state.blockedDOFs = 'zXYZ' #given that we are doing no friction, might as well block torque too

if (bc == 9) or (bc == 10):
        raise Exception('BC 9 or 10 not implemented into 8')

print("num particles: " + str(len(ballIds)) + ", packing factor = " + str(packing_factor) )

# SAMUEL MOVED CODE
# Moved this up here to be used in the sphereGen() function
# shouldn't cause issues as it won't be called unless necessary
O.materials.append(FrictMat(young=young_wl,poisson=0.3,frictionAngle=finalFricDegree, density=8000,  label='temp_walls'))
#ballLoc = 'back'
ballLoc = 'back'

def sphereGen():
        # SAMUEL ADDED CODE
        # This code generates the wall sphere
        # Sphere radius is calculated such that a specific sphere cap volume interacts with the particles using a specified indent depth.
        # Sphere has dense material properties properly work with ballvel(), otherwise sphere would just bounce off the particles due to lower density.
        # ----------------------------------
        global ballLoc
        pocketVolume = box_x * box_y * box_z
        # indentVolume = 0.35 * (1 * (10**(-11)))
        indentVolume = 0.0015 * pocketVolume
        indentDepth = thick
        sphereRad = indentVolume / (math.pi * indentDepth ** 2) + (indentDepth / 3)

        if sphereRad < indentDepth:
                indentDepth = indentDepth/2

        # a = ((2*indentDepth*(3*sphereRad - indentDepth)-(indentDepth**2))/3)**0.5             # old code, no idea what I was smoking when I derived this (from the spherical cap volume formula)

        # Should prevent the sphere from growing to crazy sizes
        print('pocket volume = ' + str(pocketVolume))
        print('indent volume = ' + str(indentVolume))
        while True:
                sphereRad = (indentVolume / (math.pi * indentDepth**2)) + (indentDepth / 3)
                print('sphereRad: ' + str(sphereRad))
                a = ((6*indentVolume/(np.pi*indentDepth)-indentDepth**3)/3)**0.5
                # I would love to use a match case here but if/elif will have to do for compatibility reasons
                if ballLoc in ('back', 'front'):
                        if a > ((0.9 * box_x / 2) or (0.9 * box_y / 2)):
                                indentDepth = indentDepth * 1.5
                                print('box_x: ' + str(box_x) + ', box_y: ' + str(box_y))
                                print('indent depth up, a: ' + str(a))
                                input()
                        #elif a < (0.1 * side / 2):
                        elif (a < min((0.1 * box_x / 2), (0.1 * box_y / 2))) or (sphereRad < (particle_dia)):
                                indentDepth = indentDepth / 2
                                print('indent depth down')
                                #input()
                        else:
                                break
                elif ballLoc in ('top', 'bottom'):
                        print('not implemented yet')
                        if a > ((0.9 * box_x / 2) or (0.9 * box_z / 2)):
                                indentDepth = indentDepth * 1.5
                                print('box_x: ' + str(box_x) + ', box_z: ' + str(box_z))
                                print('indent depth up, a: ' + str(a))
                                input()
                        #elif a < (0.1 * side / 2):
                        elif (a < min((0.1 * box_x / 2), (0.1 * box_z / 2))) or (sphereRad < (particle_dia)):
                                indentDepth = indentDepth / 2
                                print('indent depth down')
                                #input()
                        else:
                                break
                elif ballLoc in ('right', 'left'):
                        print('not implemented yet')
                        if a > ((0.9 * box_z / 2) or (0.9 * box_y / 2)):
                                indentDepth = indentDepth * 1.5
                                print('box_z: ' + str(box_z) + ', box_y: ' + str(box_y))
                                print('indent depth up, a: ' + str(a))
                                input()
                        #elif a < (0.1 * side / 2):
                        elif (a < min((0.1 * box_z / 2), (0.1 * box_y / 2))) or (sphereRad < (particle_dia)):
                                indentDepth = indentDepth / 2
                                print('indent depth down')
                                #input()
                        else:
                                break
                else:
                        print('!!!')
                        print('ballLoc is not defined properly. Choose back, front, right, left, top, or bottom')
                        ballLoc = input('ballLoc = ')
        print('indent volume = '+str((1/3)*math.pi*(indentDepth**2)*(3*sphereRad-indentDepth)))

        print('sphere radius = ' + str(sphereRad))
        # O.materials.append(FrictMat(young=young_wl,poisson=0.3,frictionAngle=finalFricDegree, density=0.1, label='massless'))
        if indent:
                if ballLoc == 'right':
                        wallSphere = yade.utils.sphere(
                                (box_x + thick + sphereRad, 1.5 * box_y + 0, box_z / 2),
                                sphereRad, material='temp_walls', wire=True)
                        # wallSphere = yade.utils.sphere(
                        #         (0 + box_x + 0 + 0.002, 1.5 * box_y + 0, box_z / 2),
                        #         thick * 3, material='temp_walls', wire=True)
                elif ballLoc == 'left':
                        wallSphere = yade.utils.sphere(
                                (-thick - sphereRad, 1.5 * box_y + 0, box_z / 2),
                                sphereRad, material='temp_walls', wire=True)
                elif ballLoc == 'back':
                        wallSphere = yade.utils.sphere(
                                (0 + box_x / 2 + 0, 1.5 * box_y + 0, -thick - sphereRad),
                                sphereRad, material='temp_walls', wire=True)
                        # wallSphere = yade.utils.sphere(
                        #         (0 + box_x / 2 + 0, 1.5 * box_y + 0, box_z / 2 + 0.001),
                        #         thick * 3, material='x_walls', wire=True)
                elif ballLoc == 'front':
                        wallSphere = yade.utils.sphere(
                                (0 + box_x / 2 + 0, 1.5 * box_y + 0, box_z + thick + sphereRad),
                                sphereRad, material='temp_walls', wire=True)
                elif ballLoc == 'top':
                        wallSphere = yade.utils.sphere(
                                (0 + box_x / 2 + 0, 2*box_y + 0 + thick + sphereRad, box_z / 2),
                                sphereRad, material='temp_walls', wire=True)
                elif ballLoc == 'bottom':
                        wallSphere = yade.utils.sphere(
                                (0 + box_x / 2 + 0, box_y -thick - sphereRad, box_z / 2),
                                sphereRad, material='temp_walls', wire=True)

                wallSphereId = O.bodies.append(wallSphere)
        return wallSphereId, indentDepth, sphereRad, indentVolume
        # ----------------------------------

if (bc == 12):
        O.materials.append(FrictMat(young=young_wl,poisson=0.3,frictionAngle=finalFricDegree, density=8000,  label='y_walls'))
        O.materials.append(FrictMat(young=young_wl,poisson=0.3,frictionAngle=finalFricDegree, density=8000,  label='x_walls'))
        O.materials.append(FrictMat(young=young_wl,poisson=0.3,frictionAngle=finalFricDegree, density=8000,  label='z_walls'))

        if (distribution == 'uniform'):
                mean_particle_dia = (particle_dia + particle_dia2)/2
        elif (distribution == 'normal'):
                mean_particle_dia = particle_dia

        wallSphereId, indentDepth, sphereRad, indentVolume = sphereGen()
        
        wallIds, topWallId, botWallId, rightWallId, leftWallId, backWallId =  wall_generate(0, 0 ,0 , bc, mean_particle_dia, box_x, box_y, box_z,mn, mx, thick, False, True, True)

        ClumpId = O.bodies.clump(wallIds)          
        O.bodies[ClumpId].state.blockedDOFs = 'xyzXYZ'
                
else:
        #generate temporary walls                
        #O.materials.append(FrictMat(young=young_wl,poisson=0.3,frictionAngle=finalFricDegree, density=8000,  label='temp_walls'))
        #create walls around the packing. big oversize factor needed b/c fast particles
        walls=utils.aabbWalls([mn,mx],thickness=thick,oversizeFactor=10,material='temp_walls')
                
        wallIds=O.bodies.append(walls)

        topWallId   = wallIds[3];  #this is the in the AFRL system sense. the excitation direction.  
        botWallId   = wallIds[2];
        leftWallId  = wallIds[0];  #x direction (gravity for U tenn)
        rightWallId = wallIds[1];

        wallSphereId, indentDepth, sphereRad, indentVolume = sphereGen()

        ClumpId = O.bodies.clump(wallIds)
        O.bodies[ClumpId].state.blockedDOFs = 'xyzXYZ'

initial_position = O.bodies[wallIds[0]].state.pos

if (len(ballIds)>500000):
        updatePeriod = 15 #plot updates
        #for extremely large models, checkUnbalanced can take more than 2 seconds, so we don't need to call it very frequently
        check_unbal = PyRunner(command='checkUnbalancedBig()', iterPeriod=100 )
        checkpoint_period = 1800        
elif (len(ballIds)>200000):
        updatePeriod = 10
        check_unbal = PyRunner(command='checkUnbalancedBig()', iterPeriod=100 )
        checkpoint_period = 1200        
elif (len(ballIds)>100000):
        updatePeriod = 5
        check_unbal = PyRunner(command='checkUnbalancedBig()',realPeriod=updatePeriod)
        checkpoint_period = 600
elif (len(ballIds)<3000):
        updatePeriod = 1
        check_unbal = PyRunner(command='checkUnbalanced7()',realPeriod=updatePeriod)
        checkpoint_period = 30
else:
        updatePeriod = 1
        check_unbal = PyRunner(command='checkUnbalanced7()',realPeriod=updatePeriod)
        checkpoint_period = 600        
                
                
#higher damping needed b/c fast particles
try:
        my_newton = NewtonIntegrator(gravity=(grav_x,grav_y ,0),damping=0.4, dampGravity=bool(damp_gravity))
except AttributeError :
        #if we get here, we are using still the older version of the code
        if (damp_gravity):
                #this is the default in the older version anyway
                my_newton = NewtonIntegrator(gravity=(grav_x,grav_y ,0),damping=0.4)
        else:
                raise Exception("tried to disable numerical damping of gravity, but using too old a version of the software")
        
        
collider = InsertionSortCollider([Bo1_Sphere_Aabb(),Bo1_Box_Aabb()],verletDist=-0.01 ) #want verlet dist small b/c very fast particles, but totally off is a big performance hit

# SAMUEL ADDED CODE
#
#

lastIter = 0
lastReal = 0
iterSpeed = 0


def simSpeedTracker():
        global lastIter, lastReal, iterSpeed
        currentIter = O.engines[-1].iterLast
        currentReal = O.engines[-1].realLast

        if lastReal > 0:
                iterSpeed = (currentIter - lastIter)/(currentReal - lastReal)
                #print('iter/sec: '+str(iterSpeed))

        lastIter = currentIter
        lastReal = currentReal



# SAMUEL ADDED CODE
# Moves sphere into the damper to simulate the indent
# Done by giving the sphere velocity until it reaches a certain distance inside the damper (dependent on indentDepth)
# Only implemented for 'back' ball location
# --------------------
ballclump = False
showActual = True
initTime = True

def ballvel():
        global initTime, timeStart
        ballVelCond = (indentDepth - sphereRad)

        # ballVelMag is meant to scale based on the simulation speed, as higher simulation speeds returned less accurate
        # indent volumes. I am unsure as to why this is, but th

        #ballVelMag = (3.96*10**-8)*(len(ballIds)**2.69)
        #ballVelMag = (2.52*10**-9)*(len(ballIds)**2.69)
        #ballVelMag = 1000
        if iterSpeed == 0:
                ballVelMag = 0
        else:
                ballVelMag = 50000/iterSpeed
        if initTime:
                timeStart = O.realtime
                initTime = False

        if indent:
                global ballclump
                global showActual
                if ballLoc == "right":

                        if O.bodies[wallSphereId].state.pos[0] > (ballVelCond + box_x):
                                O.bodies[wallSphereId].state.vel = (-5000, 0, 0)
                                #print('ballvel Active')

                        elif ballclump == False:
                                O.bodies[wallSphereId].state.vel = (0, 0, 0)
                                O.bodies.addToClump([wallSphereId], ClumpId)

                                ballclump = True
                                print('ball clumped')
                                #O.bodies.appendClumped(wallSphereId)
                                #print('ballvel Inactive')
                if ballLoc == "left":

                        if O.bodies[wallSphereId].state.pos[0] < (ballVelCond):
                                O.bodies[wallSphereId].state.vel = (5000, 0, 0)
                                #print('ballvel Active')

                        elif ballclump == False:
                                O.bodies[wallSphereId].state.vel = (0, 0, 0)
                                O.bodies.addToClump([wallSphereId], ClumpId)

                                ballclump = True
                                print('ball clumped')
                                #O.bodies.appendClumped(wallSphereId)
                                #print('ballvel Inactive')
                elif ballLoc == 'back':

                        # As of 
                        #ballDist = abs(ballVelCond - O.bodies[wallSphereId].state.pos[2])
                        timeElapsed = (O.realtime - timeStart)
                        ballVelProgress = min(timeElapsed/10, 1)
                        progressSmoothing = 1 - math.cos(ballVelProgress * math.pi / 2)
                        adjustedX = abs((-thick - sphereRad) + progressSmoothing*(ballVelCond - (-thick - sphereRad)))
                        smoothingFactor = 600
                        ballVelMag = adjustedX / (smoothingFactor * O.dt)
                        #print('ballVelMag: '+str(ballVelMag))

                        if O.bodies[wallSphereId].state.pos[2] < ballVelCond:
                                O.bodies[wallSphereId].state.vel = (0, 0, ballVelMag)
                                #print('ballvel Active')

                        elif ballclump == False:
                                O.bodies[wallSphereId].state.vel = (0, 0, 0)
                                O.bodies.addToClump([wallSphereId], ClumpId)

                                ballclump = True
                                print('ball clumped')
                                print('DEBUG')
                                print('')
                                print('wall ball rad: '+str(sphereRad))
                                print('wall ball pos: '+str(O.bodies[wallSphereId].state.pos[2]))
                                #print('back wall pos: '+str(O.bodies[backWallId].state.pos[2]))
                                print('wall thick: '+str(thick))
                                #print('backBallVelCond: '+str(backBallVelCond))
                elif ballLoc == 'front':

                        if O.bodies[wallSphereId].state.pos[2] > (ballVelCond + box_z):
                                O.bodies[wallSphereId].state.vel = (0, 0, -5000)
                                #print('ballvel Active')

                        elif ballclump == False:
                                O.bodies[wallSphereId].state.vel = (0, 0, 0)
                                O.bodies.addToClump([wallSphereId], ClumpId)

                                ballclump = True
                                print('ball clumped')
                elif ballLoc == 'top':

                        if O.bodies[wallSphereId].state.pos[1] > (ballVelCond + 2*box_y):
                                O.bodies[wallSphereId].state.vel = (0, -5000, 0)
                                #print('ballvel Active')

                        elif ballclump == False:
                                O.bodies[wallSphereId].state.vel = (0, 0, 0)
                                O.bodies.addToClump([wallSphereId], ClumpId)

                                ballclump = True
                                print('ball clumped')

                elif ballLoc == 'bottom':

                        if O.bodies[wallSphereId].state.pos[1] < (box_y + ballVelCond):
                                O.bodies[wallSphereId].state.vel = (0, 5000, 0)
                                #print('ballvel Active')

                        elif ballclump == False:
                                O.bodies[wallSphereId].state.vel = (0, 0, 0)
                                O.bodies.addToClump([wallSphereId], ClumpId)

                                ballclump = True
                                print('ball clumped')


                if ballclump:
                        with open('Indent_Loc.csv', 'w') as indentFile:
                                fileWriter = csv.writer(indentFile)
                                data = [['indentX', 'indentY', 'indentZ', 'indentRad', 'indentMat', 'wire'],
                                        [O.bodies[wallSphereId].state.pos[0], O.bodies[wallSphereId].state.pos[1], O.bodies[wallSphereId].state.pos[2],
                                         O.bodies[wallSphereId].shape.radius, O.bodies[wallSphereId].material, O.bodies[wallSphereId].shape.wire]]
                                fileWriter.writerows(data[0:])

                        # SAMUEL ADDED CODE
                        # This is debug code meant to be changed based on what side you want to look at
                        BallLoc = O.bodies[wallSphereId].state.pos[2]
                        Radius = O.bodies[wallSphereId].shape.radius
                        capHeight = BallLoc + Radius
                        indentVolumeActual = (1/3)*math.pi*(capHeight**2)*(3*Radius-capHeight)
                        if showActual:
                                print('DEBUG CODE')
                                print('----------')
                                print('target indent volume'+str(indentVolume))
                                print('actual indent volume'+str(indentVolumeActual))
                                print('----------')
                                showActual = False
# --------------------

#turn on gravity and let it settle
O.engines=[
	ForceResetter(),
        collider,
        InteractionLoop(
                [Ig2_Sphere_Sphere_ScGeom(),Ig2_Box_Sphere_ScGeom()],
                [Ip2_FrictMat_FrictMat_FrictPhys()],
                [Law2_ScGeom_FrictPhys_CundallStrack()] ),        
	my_newton,
        check_unbal,	
        PyRunner(command='addPlotData()',realPeriod=updatePeriod),
        PyRunner(command='file1.flush()', realPeriod=300),
        PyRunner(command='checkpoint(ballIds)',realPeriod=checkpoint_period ),
        PyRunner(command='escaper_check()', iterPeriod=1000, nDo=10,firstIterRun=100  ),
        PyRunner(command='ballvel()', iterPeriod=1),
        #PyRunner(command='escaper_check_sam()', iterPeriod=100),
        PyRunner(command='balls_deleted()', iterPeriod=1000),
        PyRunner(command='max_KE()', iterPeriod=1000),
        PyRunner(command='simSpeedTracker()', iterPeriod=2),
]
#small time step needed b/c fast particles
O.dt=.5*PWaveTimeStep()
                
testing = False

def addPlotData():
        if (not batch):
                if (bc == 12) and (testing == False):
                        plot.addData(time=O.time , top_force_y = np.log10(O.forces.f(topWallId)[1]),  bot_force_y = np.log10(np.abs( O.forces.f(botWallId)[1] ) ), net_force_x = np.log10(np.abs(O.forces.f(ClumpId)[0])),  right_force_x = np.log10(np.abs( O.forces.f(rightWallId)[0] ) ), time2=O.time, KE=np.log10(currentKE) )
                if (bc == 12) and (testing == True):
                        plot.addData(time=O.time , top_force_y = np.log10(O.forces.f(topWallId)[1]),  bot_force_y = np.log10(np.abs( O.forces.f(botWallId)[1] ) ), net_force_x = np.log10(np.abs(O.forces.f(ClumpId)[1])),  right_force_x = np.log10(np.abs( O.forces.f(rightWallId)[0] ) ), time2=O.time, KE=np.log10(currentKE) )
                else:
                        plot.addData(time=O.time , top_force_y = np.log10(O.forces.f(topWallId)[1]),  bot_force_y = np.log10(np.abs( O.forces.f(botWallId)[1] ) ), left_force_x = np.log10(np.abs(O.forces.f(leftWallId)[0])),  right_force_x = np.log10(np.abs( O.forces.f(rightWallId)[0] ) ), time2=O.time, KE=np.log10(currentKE) )
                

# define how to plot data: 'i' (step number) on the x-axis, raw position and order tracked magnitude on y
# show the plot on the screen, and update while the simulation runs

if (not batch):
        if (bc == 12):
                plot.plots={'time':('top_force_y', 'bot_force_y', 'net_force_x', 'right_force_x'), 'time2':('KE') }
        else:
                plot.plots={'time':('top_force_y', 'bot_force_y', 'left_force_x', 'right_force_x'), 'time2':('KE') }
        plot.plot()


if (batch):
#        O.timingEnabled = True
        O.run(-1, True)  #this will run until the O.pause hits in addPlotData
        file1.close()
#        O.run(10000, True)
#        yade.timing.stats()
        O.exitNoBacktrace()




        
