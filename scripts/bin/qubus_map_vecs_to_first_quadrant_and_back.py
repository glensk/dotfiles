#!/usr/bin/env python

import numpy as np
import utils
from scipy.ndimage        import map_coordinates
import sys

def disp_to_longvec(vec0,x,y=False,z=False):
    # in the first case x is the disp vector
    # in the second case x is just the x coordinate of the disp
    if y == False and z == False:
        return vec0-x
    else:
        return vec0-np.array([x,y,z])

def longvec_to_disp(vec0,longvec):
    return vec0-longvec

var = np.load("/Users/glensk/Dropbox/Understand_distributions/displacements_dense/Ir/2x2x2sc_qubus_3x3x3kp/qubus.F.npz")

Fx = var['Fx']
Fy = var['Fy']
Fz = var['Fz']
Ex = var['Ex']
Ey = var['Ey']
Ez = var['Ez']
all_x = var['x']
all_y = var['y']
all_z = var['z']

def get_ef_rest(x,y=False,z=False):
    if type(y) == bool and type(z) == bool:
        x,y,z = x[0],x[1],x[2]
    #print "x,y,z",x,y,z
    #sys.exit()
    scale = (all_x.size-1)/2./all_x.max()
    #print "self.x:",self.x
    #print "self.x.max()",self.x.max(),self.x.size,self.x.size,(self.x.size-1)/2,scale
    def xyz_to_map_coords(x):
        ''' x = -1.3 ---> out = 0
            x = 0    ---> out = 13
            x = 1.3  ---> out = 26 '''
        return (x+all_x.max())*scale

    coords = np.array([[xyz_to_map_coords(x), xyz_to_map_coords(y),xyz_to_map_coords(z)]])
    coords = coords.T
    #print "coords:",coords
    #zi = scipy.ndimage.map_coordinates(q.fall, coords, order=2, mode='nearest')
    ffx = map_coordinates(Fx, coords, order=2, mode='nearest')[0]
    ffy = map_coordinates(Fy, coords, order=2, mode='nearest')[0]
    ffz = map_coordinates(Fz, coords, order=2, mode='nearest')[0]
    eex = map_coordinates(Ex, coords, order=2, mode='nearest')[0]
    eey = map_coordinates(Ey, coords, order=2, mode='nearest')[0]
    eez = map_coordinates(Ez, coords, order=2, mode='nearest')[0]
    if abs(ffx) < 1e-10:ffx = 0.0
    if abs(ffy) < 1e-10:ffy = 0.0
    if abs(ffz) < 1e-10:ffz = 0.0
    if abs(eex) < 1e-10:eex = 0.0
    if abs(eey) < 1e-10:eey = 0.0
    if abs(eez) < 1e-10:eez = 0.0
    ffx = round(ffx,5)
    ffy = round(ffy,5)
    ffz = round(ffz,5)
    eex = round(eex,5)
    eey = round(eey,5)
    eez = round(eez,5)
    return [ffx,ffy,ffz],[eex,eey,eez]

def map_vec_back_to_first_quadrant(vec0_curr,longvec_curr=False,disp_curr=False):
    '''
    ##########################################################################
    # in the end we need:
    # vec0_curr
    # vec_curr  == longvec (correspongs to undisplaced vec0_curr)
    # disp_curr
    #
    # vec0_orig == sollte immer [1.955, 1.955, 0 ] sein
    # disp_orig
    # vec_orig
    ###########################################################################

    vec0_curr: is the vector in the undisplaced reference e.g. [1.955, 1.955, 0.0]
    vecs0 = np.array([ [ 1., 1.,0.], [ -1., 1.,0.], [-1.,-1.,0.], [ 1.,-1.,0.],
                   [ 1., 0.,1.], [ -1., 0.,1.], [-1.,0.,-1.], [ 1.,0.,-1.],
                   [ 0., 1.,1.], [ 0., -1.,1.], [0., -1.,-1.], [ 0., 1.,-1.],
    ])

    Rorig is the original xyz refernze frame with basis [1,0,0]
                                                        [0,1,0]
                                                        [0,0,1]
    R is the reference frme of the current longitudinal vector. Whith this function we
    convert this reference frame to Rorig, get there the energies and the forces, and map
    those back to the current referenze frame
    '''

    Rorig = np.eye(3)
    Rorig[0] = [1.0,0.0,0.0]
    Rorig[1] = [0.0,1.0,0.0]
    Rorig[2] = [0.0,0.0,1.0]
    R = np.zeros((3,3))

    # erste reihe
    R[0,0] = vec0_curr[0]
    if R[0,0] == 0.0:
        R[0,1] = vec0_curr[1]
        R[1,2] = vec0_curr[2]
        R[2] = np.cross(R[0],R[1])

    # zweite reihe
    else:  # R[0,0] ist bereits gefuellt mit einer 1 oder -1
        R[1,1] = vec0_curr[1]
        if R[1,1] != 0.0:
            R[2] = np.cross(R[0],R[1])
        if R[1,1] == 0.0:
            R[1,2] = vec0_curr[2]
            R[2] = np.cross(R[0],R[1])

    ####################### normalize matrix ########
    R = R / np.linalg.norm(R, axis=-1)[:, np.newaxis]

    ######## map vec0_curr to orig #######################
    vec0_orig = np.dot(R,vec0_curr) # sollte immer [1.955,1.955,0.0] sein
    if vec0_orig[2] != 0.0:
        sys.exit("vec0_orig[2] != 0.0")
    if vec0_orig[0] < 0.0:
        sys.exit("vec0_orig[0] < 0.0")
    if vec0_orig[1] < 0.0:
        sys.exit("vec0_orig[1] < 0.0")

    ######## get longvec_curr disp_curr  ###############################
    if type(longvec_curr) == bool and type(disp_curr) == bool:
        sys.exit("either longvec_vec or disp_curr, not both")

    if type(disp_curr) != bool:
        longvec_curr =  disp_to_longvec(vec0_curr,disp_curr)
    else:
        disp_curr = longvec_to_disp(vec0_curr,longvec_curr)

    ######## get longvec_orig disp_orig  ###############################
    longvec_orig = np.dot(R,longvec_curr)
    disp_orig = np.dot(R,disp_curr)

    ####### get forces in original position ############################
    ef_orig = get_ef_rest(disp_orig)
    ef_curr_f = np.dot(R.T,ef_orig[0])
    ef_curr_e = np.dot(R.T,ef_orig[1])


    # bei der energie macht nur die ef_orig[1] sinn, da positiv,
    # NICHT jecho die zurueckgemappte ef_curr_e
    print "{disp,longvec,vec0,disp}_curr",disp_curr,"|",longvec_curr,"|",vec0_curr
    print "{disp,longvec,vec0,disp}_orig",disp_orig,"|",longvec_orig,"|",vec0_orig
    print "ef_{orig,curr}[0]",ef_orig[0],utils.printred(ef_curr_f[0],ef_curr_f[1],ef_curr_f[2])
    print utils.printblue(ef_orig[1]),"--->",np.sum(ef_orig[1])
    print utils.printblue("ENERGY:"),ef_curr_e,"-->",np.sum(ef_curr_e)
    print utils.printyellow("FORCES:"),ef_curr_f
    print ""
    print ""
    return ef_curr_f,ef_curr_e

print "    longvec = vec                   longvec in orig        l{x,y,z}"
print "-------------------------------------------------------------------"

vecs0 = np.array([ [ 1., 1.,0.], [ -1., 1.,0.], [-1.,-1.,0.], [ 1.,-1.,0.],
                   [ 1., 0.,1.], [ -1., 0.,1.], [-1.,0.,-1.], [ 1.,0.,-1.],
                   [ 0., 1.,1.], [ 0., -1.,1.], [0., -1.,-1.], [ 0., 1.,-1.],
    ])
#for i in vecs0:
#    vec0 = i*1.995
#    map_vec_back_to_first_quadrant(vec = np.array([1.855,1.955,0.0]),vec0 = vec0)
#map_vec_back_to_first_quadrant(
#        vec  = np.array([1.855,1.955,0.0]),
#        vec0 = np.array([1.955,1.955,0.0]))
#map_vec_back_to_first_quadrant(
#        vec  = np.array([-0.1,1.955,1.955]),
#        vec0 = np.array([0.0,1.955,1.955]))
xdisp = 0.35355
ydisp = 0.35355
zdisp = 0.00000

xdisp = 0.0
ydisp = 0.1
zdisp = 0.0

map_vec_back_to_first_quadrant(
        vec0_curr       = np.array([0.0,0.0,3.99]),
        disp_curr       = np.array([xdisp,ydisp,zdisp]))
sys.exit()

print utils.printblue("# xy ebene")
map_vec_back_to_first_quadrant(
        vec0_curr       = np.array([1.955,1.955,0.0]),
        disp_curr       = np.array([xdisp,ydisp,zdisp]))
map_vec_back_to_first_quadrant(
        vec0_curr       = np.array([-1.955,1.955,0.0]),
        disp_curr       = np.array([xdisp,ydisp,zdisp]))
map_vec_back_to_first_quadrant(
        vec0_curr       = np.array([-1.955,-1.955,0.0]),
        disp_curr       = np.array([xdisp,ydisp,zdisp]))
map_vec_back_to_first_quadrant(
        vec0_curr       = np.array([1.955,-1.955,0.0]),
        disp_curr       = np.array([xdisp,ydisp,zdisp]))

print utils.printblue("# yz ebene")
map_vec_back_to_first_quadrant(
        vec0_curr       = np.array([0.0,1.955,1.955]),
        disp_curr       = np.array([xdisp,ydisp,zdisp]))
map_vec_back_to_first_quadrant(
        vec0_curr       = np.array([0.0,1.955,-1.955]),
        disp_curr       = np.array([xdisp,ydisp,zdisp]))
map_vec_back_to_first_quadrant(
        vec0_curr       = np.array([0.0,-1.955,1.955]),
        disp_curr       = np.array([xdisp,ydisp,zdisp]))
map_vec_back_to_first_quadrant(
        vec0_curr       = np.array([0.0,-1.955,-1.955]),
        disp_curr       = np.array([xdisp,ydisp,zdisp]))

print utils.printblue("# xz ebene")
map_vec_back_to_first_quadrant(
        vec0_curr       = np.array([1.955,0.0,1.955]),
        disp_curr       = np.array([xdisp,ydisp,zdisp]))
map_vec_back_to_first_quadrant(
        vec0_curr       = np.array([1.955,0.0,-1.955]),
        disp_curr       = np.array([xdisp,ydisp,zdisp]))
map_vec_back_to_first_quadrant(
        vec0_curr       = np.array([-1.955,0.0,1.955]),
        disp_curr       = np.array([xdisp,ydisp,zdisp]))
map_vec_back_to_first_quadrant(
        vec0_curr       = np.array([-1.955,0.0,-1.955]),
        disp_curr       = np.array([xdisp,ydisp,zdisp]))
print '----------------------------------------------------------------------------------'

xdisp = 0.071
ydisp = 0.00000
zdisp = 0.4

#xdisp = 0.5
#ydisp = 0.5
#zdisp = 0.0
# x=0,y=1,z=0: [ 1.658  0.642  0.   ] == ene
# x=1,y=0,z=0: [ 1.272  2.162  0.   ] == ene
#[ 1.272  0.     0.   ]  [ 3.117  6.97   0.   ]
#[ 0.     1.272  0.   ]  [ 6.97   3.117  0.   ]

# for z == 0:   get forces spline for every y
#               get energy spline for every y

#xdisp = 0.7  # [ 1.892  3.076  0.   ]
#xdisp = 0.8  # [ 2.308  4.174  0.   ]
#ydisp = 0.0
map_vec_back_to_first_quadrant(
        vec0_curr       = np.array([1.955,1.955,0.0]),
        disp_curr       = np.array([xdisp,ydisp,zdisp]))
# auslenkung in [0.4,-0.7,1.2] get forces for all 1NN perfect
# /Users/glensk/Dropbox/Understand_distributions/displacements_dense/Ir/2x2x2sc_qubus_3x3x3kp/3.99Ang_0.4_-0.7_1.2
# $cat forces_OUTCAR | head -25 | tail -1 # 0.155477 0.377736 -0.271039
# $cat forces_OUTCAR | head -29 | tail -1 # -0.200123 0.198097 -0.182764
# $cat forces_OUTCAR | head -31 | tail -1 # -0.897345 -0.134890 -0.773204
# $cat forces_OUTCAR | head -27 | tail -1 # 4.007208 -2.937014 -4.512334
# $cat forces_OUTCAR | head -9  | tail -1 # -0.250454 0.992265 0.058652
# $cat forces_OUTCAR | head -10 | tail -1 # -0.021998 0.025864 0.056919
# $cat forces_OUTCAR | head -11 | tail -1 # -17.939098 -52.728393 31.527082
# $cat forces_OUTCAR | head -12 | tail -1 # -0.079064 -0.180217 -0.255576
# $cat forces_OUTCAR | head -17 | tail -1 # 20.805627 10.967347 9.323492
# $cat forces_OUTCAR | head -18 | tail -1 # 0.121554 0.097666 -0.154984
# $cat forces_OUTCAR | head -21 | tail -1 # -2.312723 1.251293 0.225208
# $cat forces_OUTCAR | head -22 | tail -1 # -0.104589 0.063529 -0.028853


# auslenkung in [0.35355,0.35355,0.0] get forces for all 1NN perfect
# /Users/glensk/Dropbox/Understand_distributions/displacements_dense/Ir/2x2x2sc_quer_3x3x3kp/3.99Ang_0.5
# $cat cartesian_coords | head -25 | tail -1 # 1.99500 1.99500 0.00000
#
# $cat cartesian_coords | head -29 | tail -1 # 5.98500 1.99500 0.00000
# $cat forces_OUTCAR    | head -29 | tail -1 # 0.050216 0.001746 0.000000
# $cat cartesian_coords | head -31 | tail -1 # 5.98500 5.98500 0.00000
# $cat forces_OUTCAR    | head -31 | tail -1 # 0.062074 0.062074 0.000000
# $cat cartesian_coords | head -27 | tail -1 # 1.99500 5.98500 0.00000
# $cat forces_OUTCAR    | head -27 | tail -1 # 0.001746 0.050216 0.000000
# $cat cartesian_coords | head -9  | tail -1 # 0.00000 1.99500 1.99500
# $cat forces_OUTCAR    | head -9  | tail -1 # -0.129400 0.578717 0.712410
# $cat cartesian_coords | head -10 | tail -1 # 0.00000 1.99500 5.98500
# $cat forces_OUTCAR    | head -10 | tail -1 # -0.129400 0.578717 -0.712410
# $cat cartesian_coords | head -11 | tail -1 # 0.00000 5.98500 1.99500
# $cat forces_OUTCAR    | head -11 | tail -1 # -0.041858 0.127368 -0.102035
# $cat cartesian_coords | head -12 | tail -1 # 0.00000 5.98500 5.98500
# $cat forces_OUTCAR    | head -12 | tail -1 # -0.041858 0.127368 0.102035
# $cat cartesian_coords | head -17 | tail -1 # 1.99500 0.00000 1.99500
# $cat forces_OUTCAR    | head -17 | tail -1 # 0.578717 -0.129400 0.712410
# $cat cartesian_coords | head -18 | tail -1 # 1.99500 0.00000 5.98500
# $cat forces_OUTCAR    | head -18 | tail -1 # 0.578717 -0.129400 -0.712410
# $cat cartesian_coords | head -21 | tail -1 # 5.98500 0.00000 1.99500
# $cat forces_OUTCAR    | head -21 | tail -1 # 0.127368 -0.041858 -0.102035
# $cat cartesian_coords | head -22 | tail -1 # 5.98500 0.00000 5.98500
# $cat forces_OUTCAR    | head -22 | tail -1 # 0.127368 -0.041858 0.102035
#

