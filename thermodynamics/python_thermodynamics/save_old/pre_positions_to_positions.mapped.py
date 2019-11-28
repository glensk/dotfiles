#!/usr/bin/env python


import hesse
import crystal
import os
import numpy as np
import utils
import sys

posall=np.loadtxt("pre_POSITIONs")[:,0:3]
pos0 = np.loadtxt("pos0")
atoms=32
schritte = posall.shape[0]/atoms


for i in np.arange(schritte):
    pos = posall[i*atoms:i*atoms+32]
    u = pos-pos0

    # mapping back:
    cell=np.loadtxt("cell")
    for ind,j in enumerate(u):
        if  u[ind][0] > cell[0,0]/2:
            u[ind][0] = u[ind][0] - cell[0,0]
        if  u[ind][1] > cell[0,0]/2:
            u[ind][1] = u[ind][1] - cell[0,0]
        if  u[ind][2] > cell[0,0]/2:
            u[ind][2] = u[ind][2] - cell[0,0]

    posall[i*atoms:i*atoms+32] = u + pos0

np.savetxt("pre_POSITIONs_mapped",posall,fmt="%.6f %.6f %.6f")
print "written pre_POSITIONs_mapped"


sys.exit()
#
#
#
#
#
#
#
#
#h = hesse.read_Hessematrix()
#
## pos at zero
#crystal = crystal.crystal()
#
#if os.path.isfile("POSCAR") != True:
#    sys.exit("please provide a POSCAR file, this is necessary to create a cell file")
#if os.path.isfile("EqCoords_direct") != True:
#    sys.exit("please provida a EqCoords_direct file")
#if os.path.isfile("cell") != True:
#    import utils
#    utils.run2("rm -f cell; POSCAR_cell_cartesian.sh > cell") # to create cell
#
##print "---------------- harmonic --------------------"
#crystal.read_rrel_positions( coordfile = "EqCoords_direct", cellname = "cell")
#pos0=crystal.rcar
#pos1=np.loadtxt("cartesian_coords")
#
#
## mapping back:
#cell=np.loadtxt("cell")
#u = pos1-pos0
#for ind,i in enumerate(u):
#    if  u[ind][0] > cell[0,0]/2:
#        u[ind][0] = u[ind][0] - cell[0,0]
#    if  u[ind][1] > cell[0,0]/2:
#        u[ind][1] = u[ind][1] - cell[0,0]
#    if  u[ind][2] > cell[0,0]/2:
#        u[ind][2] = u[ind][2] - cell[0,0]
#
#
#
#
#
#fh = hesse.qh_forces(u,h)
#eh = hesse.qh_energy_cell(u,h)
#np.savetxt("forces",fh)
#np.savetxt("energy",np.array([eh]),fmt="%.6f")
#
#
#eip0, fip0 = hesse.energy_forces_inversepot(positionsfilein = "pos0", cellfilein = "cell", NNshell = 1)
#eip, fip = hesse.energy_forces_inversepot(positionsfilein = "cartesian_coords", cellfilein = "cell", NNshell = 1)
##print "---------------- sum --------------------"
#eip = eip - eip0
#print eh,eip,eh+eip
##print fh+fip
#
#np.savetxt("forces",fh*1.+fip*8.)
#np.savetxt("energy",np.array([eh*1.+eip*8.]))
#
