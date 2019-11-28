#!/usr/bin/env python

import glob
import sys
import os
import numpy as np

def lammps_check_if_trj_lammps_and_trj_lammpsnew_equal(f1,f2,rm=False):
    c1 = os.popen('tail -1 '+f1+" | awk '{print $NF}'").read().rstrip().rstrip()
    c2 = os.popen('tail -1 '+f2+" | awk '{print $NF}'").read().rstrip().rstrip()
    c3 = os.popen('tail -1 '+f1+" | awk '{print $(NF-1)}'").read().rstrip().rstrip()
    c4 = os.popen('tail -1 '+f2+" | awk '{print $(NF-1)}'").read().rstrip().rstrip()
    print "c1,c2:",c1,c2
    print "c3,c4:",c3,c4
    if c1 == c2 and c3==c4:
        print "c1 == c2 and c3 == c4"
        if type(rm) != bool:
            if rm == f1 or rm == f2:
                print "removing:",rm
                os.remove(rm)
        return True

    else:
        print f1,"and",f2,"are not equal!"
        return False


def lammps_split_lammps_trajectory_in_xaa_files_and_app_c_skript(filenamelammps="trj_lammps.out",filenamepos="trj_lammpsnew.out",positionsfilename = "positions.*",infilefilename = "in_file_dynamics.in",linesokformemory = 800000000):
    ''' splits lammps trajectory file in xaa, xab, xac ...'''
    positionsfile = glob.glob(positionsfilename)
    if len(positionsfile) != 1:
        print "positionsfile:",positionsfile
        sys.exit("Not one positionsfile but "+str(len(positionsfile)))

    print "positionsfile:",positionsfile

    infile = glob.glob(infilefilename)
    print "infile:",infile
    if len(infile) != 1:
        print "infile:",infile
        sys.exit("Not one infile but "+str(len(infile)))

    atoms1 = os.popen('wc -l '+positionsfile[0] + "| awk '{print $1}'").read()
    atoms = int(atoms1)-8
    sc = N = round(float((atoms/4.)**(1/3.)),0)
    stepslammps = int(os.popen('grep "^run " '+infile[0] + "| awk '{print $2}'").read().rstrip())
    dump = int(os.popen('grep "^dump dump1" '+infile[0] + "| awk '{print $5}'").read().rstrip())

    print "atoms:",atoms,type(atoms)
    print "sc:",sc
    print "stepslammps:",stepslammps
    print "dump:",dump
    stepswritten = stepslammps/dump+1
    print "stepswritten:",stepswritten
    lineswritten = stepswritten*atoms
    print "lineswritten:",lineswritten
    print "linesokformemory:",linesokformemory
    nsteps = linesokformemory / atoms
    nsteps = int(round(linesokformemory / atoms/10,0)*10)
    print "nsteps:",nsteps
    from string import ascii_lowercase
    stepsremain = stepswritten



    if os.path.isfile(filenamelammps):
        print filenamelammps,"exists!"
        if not os.path.isfile(filenamepos):
            print filenamepos,"does not exist! --> creating it!"
            os.popen('grep "^1" '+filenamelammps + "| awk '{print $2,$3,$4}' > trj_lammpsnew.out").read()

    if os.path.isfile(filenamelammps) and os.path.isfile(filenamepos):
        lammps_check_if_trj_lammps_and_trj_lammpsnew_equal(filenamelammps,filenamepos,rm=filenamelammps)

    print
    print
    splitanzahl=-1
    print "      stepstot","\t","steps  ","\t","stepsremain"
    print "---------------------------------------------------"
    lastfile = False
    lastamoutoflines = False
    # make here an array
    for idx,c in enumerate(ascii_lowercase):
        nstepscurr = nsteps
        if (stepsremain - nsteps) > 0:
            stepsremain = stepsremain - nsteps
            exittrue = False
        else:
            nstepscurr = stepsremain
            stepsremain = 0
            exittrue = True
        splitanzahl = splitanzahl+1
        lines = nstepscurr*atoms
        print idx,"xa"+c,stepswritten,"\t",nstepscurr,"\t",stepsremain,splitanzahl,lines
        lastfile = "xa"+c
        lastamoutoflines = lines
        if exittrue:
            break
    print "splitanzahl:",splitanzahl
    print
    print

    print
    print
    #for i in np.arange(20)+1:
    #    atoms = 4*i**3
    #    print  i,atoms,linesokformemory / atoms, int(round(linesokformemory / atoms/10,0)*10)

    print
    print
    if splitanzahl == 0:
        fileforsaschasskript = filenamepos
    if splitanzahl > 0:
        fileforsaschasskript = "xaa"
        print "splitanzhal > 0, in fact it is:",splitanzahl
        if not os.path.isfile("xaa"):
            print "splitting by:",str(int(nsteps*atoms)),"since xaa does not exist! --> splitting"
            os.popen('split -l '+str(int(nsteps*atoms))+" "+filenamepos).read()
    print "fileforsaschasskript:",fileforsaschasskript

    print "lastfile:",lastfile
    print "lastamoutoflines:",lastamoutoflines
    if splitanzahl > 0 and os.path.isfile("xaa") and os.path.isfile(lastfile) and os.path.isfile(filenamepos):
        lammps_check_if_trj_lammps_and_trj_lammpsnew_equal(filenamepos,lastfile,rm=filenamepos)

    #print "##############################################################################"
    #print "now applying C++ skript"
    #print "##############################################################################"
    #lammps_create_c_SETTINGS_file(infile="xaa",atoms=atoms,steps=)
    return

def lammps_create_c_SETTINGS_file(filename="SETTINGS",infile=False,atoms=False,steps=False):
    ''' infile: xaa'''
    if type(infile) == bool:
        sys.exit("specify infile")
    if type(atoms) == bool:
        sys.exit("specify atoms")
    if not os.path.exists(infile+"_"):
        os.makedirs(infile+"_")
    settingsfile = infile+"_"+"/SETTINGS"
    if os.path.isfile(settingsfile):
        os.remove(settingsfile)

    target = open(settingsfile, 'w')
    target.write("NATOMS = "+str(atoms)+"\n")
    target.write("NSTEPS = "+str(steps)+"\n")
    target.close()


lammps_split_lammps_trajectory_in_xaa_files_and_app_c_skript(filenamelammps="trj_lammps.out",filenamepos="trj_lammpsnew.out",positionsfilename = "positions.*",infilefilename = "in_file_dynamics.in",linesokformemory = 800000000)
