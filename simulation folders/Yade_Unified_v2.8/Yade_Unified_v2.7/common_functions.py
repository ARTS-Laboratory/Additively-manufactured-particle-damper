import yade
import numpy as np
import cmath
from scipy import signal
import inspect


#boundary condition options
#0: normal 3D space, all 6 walls are the same
#2: 2D (4 walls)
#7: normal 3D space, top & bottom walls have different material than side walls

#8: like 7, but one giant column. no intermediate walls in between the pockets.   

#9: 3D space, non-smooth top & bottom (Excitation direction) walls (branch not merged as made things more complicated and didnt improve correlation)
#10: 3D space, non-smooth side walls

#11: like 8, but bug fix in yade_utennessee, for the ordering of the pockets. identical as far as this code is concerned

#12: like 11, but non-smooth rightmost wall (i.e. minus gravity direction for U Tenn)

#13: like 8/11/12, but one big row instead of one big column. 




# #new version, will handle non-smooth side walls (BC 10).  integrated into ball_gen_14, but not yet integrated into yade_utennessee.py
# #output list has changed versus the older version.
# #want to integrate side wall friction first, as that might be sufficient.
# def wall_generate_new(box_off_x,box_off_y, RightWallDelta, periodic_bc, particle_dia, box_x, box_y, box_z,mn, mx):

#     thick = 2*particle_dia # the thickness of the walls
            
#     box_off = (box_off_x, box_off_y,0)
    
#     if (periodic_bc==1):
#         raise Exception('3D space with periodic boundary conditions. Removed')        
#     elif (periodic_bc==0):
#         #normal 3D space, all 6 walls are the same
        
#         osf = 1.5        
            
#         #create walls around the packing
#         walls=utils.aabbWalls([mn+box_off,mx+box_off],thickness=thick,oversizeFactor=osf,material='y_walls')
#         wallIdsAll=O.bodies.append(walls)
        
#         topWallId = wallIds[2];  #fixme these labels are reversed!
#         botWallId = wallIds[3];                
#     elif (periodic_bc==7) or (periodic_bc==10) :
#         #normal 3D space, top & bottom walls have different material than side walls
#         expand=1.5

#         #create walls around the packing
#         topWall   = yade.utils.box((box_x/2+box_off_x,2*box_y+box_off_y+thick/2,box_z/2), (box_x/2*expand,thick/2,box_z/2*expand), material='y_walls', wire=True)
#         topWallId =  O.bodies.append(topWall)
#         wallIdsAll = [topWallId]

#         botWall   = yade.utils.box((box_x/2+box_off_x,box_y+box_off_y-thick/2,box_z/2), (box_x/2*expand,thick/2, box_z/2*expand), material='y_walls', wire=True)
#         botWallId =  O.bodies.append(botWall)
#         wallIdsAll.append(botWallId)
                        
#         rightWall   = yade.utils.box((RightWallDelta+box_x+thick/2+box_off_x,1.5*box_y+box_off_y,box_z/2), (thick/2,box_y/2*expand, box_z/2*expand), material='sidewalls', wire=True)
#         tmp  =  O.bodies.append(rightWall)
#         wallIdsAll.append(tmp)
#         wallIdsRight = [tmp]
        
#         leftWall   = yade.utils.box((-thick/2+box_off_x,1.5*box_y+box_off_y,box_z/2), (thick/2,box_y/2*expand, box_z/2*expand), material='sidewalls', wire=True)
#         leftWallId =  O.bodies.append(leftWall)
#         wallIdsAll.append(leftWallId)
        
#         frontWall   = yade.utils.box((box_x/2+box_off_x,1.5*box_y+box_off_y,box_z+thick/2), (box_x/2*expand,box_y/2*expand, thick/2), material='sidewalls', wire=True)
#         frontWallId =  O.bodies.append(frontWall)
#         wallIdsAll.append(frontWallId)
        
#         backWall = yade.utils.box((box_x/2+box_off_x,1.5*box_y+box_off_y,-thick/2), (box_x/2*expand,box_y/2*expand, thick/2), material='sidewalls', wire=True)
#         backWallId =  O.bodies.append(backWall)
#         wallIdsAll.append(backWallId)

#         if  (periodic_bc == 10):
#                 #additionally add particles so not smooth
        
#                 #note: must make sure that these don't overlap, and that there is enough room for the right wall to move freely between
#                 #the other walls
#                 wall_particles_x = np.arange( mn[0]+box_off[0]+particle_dia,  mx[0]+box_off[0], particle_dia);
#                 wall_particles_y = np.arange( mn[1]+box_off[1]+particle_dia,  mx[1]+box_off[1], particle_dia);

                            
#                 for x in wall_particles_x:
#                         for y in wall_particles_y:
#                                 center = ( x, y, box_z)
#                                 id = O.bodies.append(yade.utils.sphere(center, particle_dia/2 ,material='sidewalls',color=(1,0,0)))                    
#                                 wallIdsAll.append(id)

#                                 center = ( x, y, 0)
#                                 id = O.bodies.append(yade.utils.sphere(center, particle_dia/2 ,material='sidewalls',color=(0,1,0)))                    
#                                 wallIdsAll.append(id)

#                                 center = ( box_x, y,x)
#                                 id = O.bodies.append(yade.utils.sphere(center, particle_dia/2 ,material='sidewalls',color=(0,0,1)))                    
#                                 wallIdsAll.append(id)
#                                 wallIdsRight.append(id)

#                                 center = ( 0, y,x)
#                                 id = O.bodies.append(yade.utils.sphere(center, particle_dia/2 ,material='sidewalls',color=(1,1,0)))                    
#                                 wallIdsAll.append(id)


#                 rightWallClumpId = O.bodies.clump(wallIdsRight)
#                 #rightWallClumpId = []
#         else:
#             rightWallClumpId = []
        
#     else:
#         raise Exception('unknown BC')

#     for i in wallIdsAll:
#             O.bodies[i].state.isDamped = False #this was inadvertently left on in original AFRL code

#     return wallIdsAll, wallIdsRight, topWallId, botWallId, rightWallClumpId, thick


#old version.  does not handle BC 9 non-smooth side walls.  does handle BC 12 non-smooth leftmost wall
def wall_generate(box_off_x,box_off_y, RightWallDelta, bc, particle_dia, box_x, box_y, box_z,mn, mx, thick, backstop=False, leftmost=False, rightmost=False, z_plus_wall=False, z_minus_wall=False, box_off_z=0):
    
    box_off = (box_off_x, box_off_y,box_off_z)
    
    if (bc==1):
        raise Exception('3D space with periodic boundary conditions. Removed')        
    elif (bc==0):
        raise Exception('normal 3D space, all 6 walls are the same. removed')
    elif (bc==2):
        raise Exception('removed')            
    elif (bc==7):
        #normal 3D space, top & bottom walls have different material than side walls
        expand=2.5 #was 1.5, but that was not enough for u tenn system.  some balls were escaping out the right.

        #create walls around the packing
        topWall   = yade.utils.box((box_x/2+box_off_x,2*box_y+box_off_y+thick/2,box_z/2), (box_x/2*expand,thick/2,box_z/2*expand), material='y_walls', wire=True)
        topWallId =  O.bodies.append(topWall)
        wallIds = [topWallId]

        botWall   = yade.utils.box((box_x/2+box_off_x,box_y+box_off_y-thick/2,box_z/2), (box_x/2*expand,thick/2, box_z/2*expand), material='y_walls', wire=True)
        botWallId =  O.bodies.append(botWall)
        wallIds.append(botWallId)
                        
        rightWall   = yade.utils.box((RightWallDelta+box_x+thick/2+box_off_x,1.5*box_y+box_off_y,box_z/2), (thick/2,box_y/2*expand, box_z/2*expand), material='x_walls', wire=True)
        rightWallId =  O.bodies.append(rightWall)
        wallIds.append(rightWallId)

        if (backstop):
            #for utennesee system, to prevent light right walls from escaping
            #for the 212um box, this allows the right wall to move ~113 um, or +50%.
            #not sure if it would scale correctly to other box sizes
            id1 = O.bodies.append(yade.utils.sphere((RightWallDelta+1.7*box_x+box_off_x,1.4*box_y+box_off_y,box_z/2) , particle_dia/2 ,material='x_walls',color=(1,1,1)))
            O.bodies[id1].state.blockedDOFs='xyzXYZ'
            
        leftWall   = yade.utils.box((-thick/2+box_off_x,1.5*box_y+box_off_y,box_z/2), (thick/2,box_y/2*expand, box_z/2*expand), material='x_walls', wire=True)
        leftWallId =  O.bodies.append(leftWall)
        wallIds.append(leftWallId)
        
        frontWall   = yade.utils.box((box_x/2+box_off_x,1.5*box_y+box_off_y,box_z+thick/2), (box_x/2*expand,box_y/2*expand, thick/2), material='z_walls', wire=True)
        frontWallId =  O.bodies.append(frontWall)
        wallIds.append(frontWallId)
        
        backWall = yade.utils.box((box_x/2+box_off_x,1.5*box_y+box_off_y,-thick/2), (box_x/2*expand,box_y/2*expand, thick/2), material='z_walls', wire=True)
        backWallId =  O.bodies.append(backWall)
        wallIds.append(backWallId)
    elif (bc == 8) or (bc == 11) or (bc == 12):
        #one giant column
        # SAMUEL ADDED NOTE
        # Remember to set expansion higher for gravity deposition sims, then back to 1 for yade_sdoftrial9 simulations
        if inspect.currentframe().f_back.f_code.co_filename == 'ball_gen_8_and_15.py':
            expand = 3
        else:
            expand = 1 #make sure one wall doesn't interfere with the next set of particles

        #create walls around the packing
        topWall   = yade.utils.box((box_x/2+box_off_x,2*box_y+box_off_y+thick,box_z/2), (box_x/2*expand,thick,box_z/2*expand), material='y_walls', wire=True)
        topWallId =  O.bodies.append(topWall)
        wallIds = [topWallId]

        botWall   = yade.utils.box((box_x/2+box_off_x,box_y+box_off_y-thick,box_z/2), (box_x/2*expand,thick, box_z/2*expand), material='y_walls', wire=True)
        botWallId =  O.bodies.append(botWall)
        wallIds.append(botWallId)

        if (rightmost):
            rightWall   = yade.utils.box((RightWallDelta+box_x+thick+box_off_x,1.5*box_y+box_off_y,box_z/2), (thick,box_y/2*expand, box_z/2*expand), material='x_walls', wire=True)
            rightWallId =  O.bodies.append(rightWall)
            wallIds.append(rightWallId)
            print('rightmost wall, ' + str(rightWallId) + " " + str(RightWallDelta))                               
        else:
            rightWallId  = None
        
        if (leftmost):
            leftWall   = yade.utils.box((-thick+box_off_x,1.5*box_y+box_off_y,box_z/2), (thick,box_y/2*expand, box_z/2*expand), material='x_walls', wire=True)
            leftWallId =  O.bodies.append(leftWall)
            wallIds.append(leftWallId)

            # SAMUEL ADDED NOTE
            # This generates red balls along one wall. For some reason, not having this active breaks non-indent sims
            if (bc == 12):
                wall_particles_y = np.arange( mn[1]+box_off[1]+particle_dia,  mx[1]+box_off[1], particle_dia);
                wall_particles_z = np.arange( mn[2]+box_off[2]+particle_dia,  mx[2]+box_off[2], particle_dia);

                print(wall_particles_y)
                print(wall_particles_z)

                for z in wall_particles_z:
                        for y in wall_particles_y:
                                center = ( box_off_x, y, z)
                                id = O.bodies.append(yade.utils.sphere(center, particle_dia/2 ,material='x_walls',color=(1,0,0)))
                                wallIds.append(id)
                                
        else:
            leftWallId  = None
            
        frontWall   = yade.utils.box((box_x/2+box_off_x,1.5*box_y+box_off_y,box_z+thick), (box_x/2*expand,box_y/2*expand, thick), material='z_walls', wire=True)
        frontWallId =  O.bodies.append(frontWall)
        wallIds.append(frontWallId)
        
        backWall = yade.utils.box((box_x/2+box_off_x,1.5*box_y+box_off_y,-thick), (box_x/2*expand,box_y/2*expand, thick), material='z_walls', wire=True)
        backWallId =  O.bodies.append(backWall)
        wallIds.append(backWallId)

    elif (bc == 13):
        #one giant row
        expand=1.0 #make sure one wall doesn't interfere with the next set of particles

        topWall   = yade.utils.box((box_x/2+box_off_x,2*box_y+box_off_y+thick/2,box_z/2+box_off_z), (box_x/2*expand,thick/2,box_z/2*expand), material='y_walls', wire=True)
        topWallId =  O.bodies.append(topWall)
        wallIds = [topWallId]

        botWall   = yade.utils.box((box_x/2+box_off_x,box_y+box_off_y-thick/2,box_z/2+box_off_z), (box_x/2*expand,thick/2, box_z/2*expand), material='y_walls', wire=True)
        botWallId =  O.bodies.append(botWall)
        wallIds.append(botWallId)

        
        if (z_plus_wall):
            frontWall   = yade.utils.box((box_x/2+box_off_x,1.5*box_y+box_off_y,box_z+thick/2+box_off_z), (box_x/2*expand,box_y/2*expand, thick/2), material='z_walls', wire=True)
            frontWallId =  O.bodies.append(frontWall)
            wallIds.append(frontWallId)
        else:
            frontWallId = None

        if (z_minus_wall):
            backWall = yade.utils.box((box_x/2+box_off_x,1.5*box_y+box_off_y,-thick/2+box_off_z), (box_x/2*expand,box_y/2*expand, thick/2), material='z_walls', wire=True)
            backWallId =  O.bodies.append(backWall)
            wallIds.append(backWallId)
        else:
            backWallId = None
        
        
        rightWall   = yade.utils.box((RightWallDelta+box_x+thick/2+box_off_x,1.5*box_y+box_off_y,box_z/2+box_off_z), (thick/2,box_y/2*expand, box_z/2*expand), material='x_walls', wire=True)
        rightWallId =  O.bodies.append(rightWall)
        wallIds.append(rightWallId)

        leftWall   = yade.utils.box((-thick/2+box_off_x,1.5*box_y+box_off_y,box_z/2+box_off_z), (thick/2,box_y/2*expand, box_z/2*expand), material='x_walls', wire=True)
        leftWallId =  O.bodies.append(leftWall)
        wallIds.append(leftWallId)
        
    else:
        raise Exception('unknown BC')

    for i in wallIds:
            O.bodies[i].state.isDamped = False #this was inadvertently left on in original AFRL code

    return wallIds, topWallId, botWallId, rightWallId, leftWallId, backWallId

###########################################################################

def ball_gen_8_output_filename( target_packing_factor, particle_dia, box_x, box_y, box_z, seedType, bc, distribution, damp_gravity, std_dev_dia):
        suffix = "gen8_" + str(target_packing_factor) + "_"+ str(np.floor(particle_dia/1e-6))

        prefix = "grav_dep_pack_" + str(np.floor(box_x/1e-6)) + "_" + str(np.floor(box_y/1e-6)) + "_" + str(np.floor(box_z/1e-6))

        if (damp_gravity):
                dgstr = ""
        else:
                dgstr = "_undampGrav"
        
        if (seedType == 1):
                seedstr = ""
        else:
                seedstr = "_Seed" +str(seedType)
                
        if (bc == 12):
                bcstr = "_bc12"
        else:
                bcstr = ""
        
        if (distribution == 'uniform'):
                fn = prefix + "_uniform_" + suffix + bcstr + seedstr + dgstr + ".csv"
        elif (distribution == 'normal'):
                fn = prefix + "_normal_" + suffix + "_" + str(np.floor(std_dev_dia/1e-6)) + bcstr + seedstr + dgstr + ".csv"
        return fn


def ball_gen_14_output_filename( target_packing_factor, particle_dia,particle_dia2, box_x, box_y, box_z, seedType, distribution, std_dev_dia, DistanceToTopCurrentPocket):
        suffix = "gen14_" + str(target_packing_factor) + "_Seed" +str(seedType) + "_" + str(DistanceToTopCurrentPocket) + "_"

        if (distribution == 'uniform'):
                file_name = "grav_dep_pack_" + str(np.floor(box_x/1e-6)) + "_" + str(np.floor(box_y/1e-6)) + "_" + str(np.floor(box_z/1e-6)) + "_uniform_" + suffix + str(np.floor(particle_dia/1e-6)) + "_" + str(np.floor(particle_dia2/1e-6)) + ".csv"        
        elif (distribution == 'normal'):
                file_name = "grav_dep_pack_" + str(np.floor(box_x/1e-6)) + "_" + str(np.floor(box_y/1e-6)) + "_" + str(np.floor(box_z/1e-6)) + "_normal_" + suffix + str(np.floor(particle_dia/1e-6)) + "_" + str(np.floor(std_dev_dia/1e-6)) + ".csv"

        return file_name


def packing_z(ball_gen_method,box_z,NumPocketsForArbitrary,bc11_12_13_numCols ):
    if (ball_gen_method==19):
        return box_z*NumPocketsForArbitrary / bc11_12_13_numCols
    else:
        return box_z
    
def packing_x(ball_gen_method,box_x,NumPocketsForArbitrary,bc11_12_13_numCols):
    if (ball_gen_method == 15):
        return box_x*NumPocketsForArbitrary / bc11_12_13_numCols
    else:
        return box_x

########################################3
# this section mostly copied from VEDA

def ZN_PI( k_ultimate, tau_ultimate):
      Kp = k_ultimate / 2.2
      Ki = 1.2 * Kp / tau_ultimate 
      Kd = 0
      return Kp, Ki, Kd

def ZN_PID( k_ultimate, tau_ultimate):
      Kp = k_ultimate * 0.6
      Ki = 2 * Kp / tau_ultimate 
      Kd = Kp * tau_ultimate / 8
      return Kp, Ki, Kd

#the transfer function for the slow time scale dynamics of the amplitude
def AmpXferFunc(omega_d, Q, omega_n, b,a,Fs):
     w,filtergain = signal.freqz(b,a,  [(omega_d/np.pi/2)], fs=Fs )
     return filtergain[0] * ((omega_n / (2 * Q)) / complex( omega_n / (2 * Q), omega_d))

# second order filter never hits -pi?  am I using something different here than in VEDA?  VEDA comments indicated that 1st order would not work but 2nd order would
 
#    return LockinGain(omega_d, PLL_filterBW) * ((omega_n / (2 * Q)) / complex( omega_n / (2 * Q), omega_d))

#hard coded for second order butterworth
#def LockinGain( omegad, PLL_filterBW ):
#    LockinTC = 1 / PLL_filterBW
#    return 1 / (complex( np.sqrt(2) / 2, np.sqrt(2)/2 +  omegad * LockinTC) * complex( np.sqrt(2)/2, -  np.sqrt(2)/2 + omegad * LockinTC))

def calc_PLL_gain(Q, omega_n, b,a,Fs):
    #VEDA starts with amplitude controller
    #step 1, find the frequency at which the angle of the open loop
    #transfer function is -180 degree.  call it omega_180

    if (Q <= 0):
        raise Exception('Invalid Q guess')
    
    print("Q = " + str(Q) + "omega_n = " + str(omega_n) + " Fs = " + str(Fs) )
    
    #simple bisection root finding
    #essentially, we're looking for the wrap around point
    left  = 0
    right = 10 * omega_n

    while( (right -left) > 0.001):
         guess = (left+right)/2
         aaa = np.angle(AmpXferFunc(guess,Q, omega_n, b,a,Fs))
#         print("guess = " + str(guess) + "angle = " + str(aaa) )
         if ( aaa  > 0 ) :
            right = guess
         else:
            left = guess

    omegad_180 = (left+right)/2

    print("omegad180 = " + str(omegad_180))
    
    #step 2, ultimate proportional gain is that which makes transfer function
    #have unity magnitude at omega_180.
    k_ultimate_amp = 1 / abs( AmpXferFunc(omegad_180,Q, omega_n, b,a,Fs))

    print("k_ultimate_amp = " + str(k_ultimate_amp))
    
    #step 3, oscillation period at ultimate gain. 
    tau_ultimate = (2 * np.pi) / omegad_180 

    print("tau_ultimate = " + str(tau_ultimate))
    
    #step 4, phase dynamics for PLL.  transfer function is actually the same, just scaled differently.
    #!in particular, need the slope of the phase at natural frequency, which is 2 * Q / omega_i
    #this what VEDA had
    k_ultimate_pll = k_ultimate_amp / ( 2 * Q / omega_n )
    #this is what seems to make it work. not sure the discrepancy?
#    k_ultimate_pll = k_ultimate_amp  *  ( 2 * Q / omega_n )

#    k_ultimate_pll =     k_ultimate_pll / (2*np.pi)
    
    print("k_ultimate_pll = " + str(k_ultimate_pll))
    
    return ZN_PID( k_ultimate_pll, tau_ultimate)
