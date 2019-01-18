#!/usr/bin/env python
from __future__ import print_function

import click
import tarfile
import re,os,sys
import subprocess,os
import numpy as np
import time
import scipy.linalg as salg
import scipy.sparse.linalg as spalg
import scipy.integrate as spint
import operator
from iolib import *
from ase.io import read as ase_read
from ase.io import write as ase_write
import myutils as my

CONTEXT_SETTINGS = my.get_click_defaults()
@click.command(context_settings=CONTEXT_SETTINGS)
@click.option('-d1','--db1',required=True,type=str,help="path to DB1 which are preexisting structures")
@click.option('-d2','--db2',required=True,type=str,help="path to DB2 which are current structures")
@click.option('-nsyms','--nsyms',required=False,default={"Mg":64, "Al":64, "Si":64},help="dictionary containing amount of symmetry functions per element")
@click.option('-s','--structures_upperlim',required=False,default=False,type=int,help="upper limit of structures to fps/input.fps.data")


def make_fps(db1,db2,nsyms,structures_upperlim):
    ''' makes the fps considering previous structures

        to get nsyms (and not needing to provide it):
         - grep symfunction_short $dotfiles/scripts/potentials/n2p2_v1ag/input.nn | awk '{print $2}' | uniq | grep -v "<.*" # -> gets the elements
         -  grep symfunction_short $dotfiles/scripts/potentials/n2p2_v1ag/input.nn | awk '{print $2}' | grep "Al" | wc -l # -> gets 64
         -  grep symfunction_short $dotfiles/scripts/potentials/n2p2_v1ag/input.nn | awk '{print $2}' | grep "Mg" | wc -l # -> gets 64
         -  grep symfunction_short $dotfiles/scripts/potentials/n2p2_v1ag/input.nn | awk '{print $2}' | grep "Si" | wc -l # -> gets 64
    '''
    DB1_path = str(os.path.abspath(db1))
    DB2_path = str(os.path.abspath(db2))
    print('DB1_path         ',DB1_path)
    print('DB2_path         ',DB2_path) #,type(DB2_path))
    print('structures_upperlim       ',structures_upperlim)
    print()

    # check if necessary files exist
    print()
    for i in [ DB1_path,DB2_path]:
        if not os.path.isdir(i): sys.exit(i+" does not exist!")
        else: print(i+" exists")
        if not os.path.isfile(i+'/input.data'): sys.exit(i+"/input.data does not exist!")
        else: print(i+"/input.data exists")
        if not os.path.isfile(i+'/function.data'):
            if os.path.isfile(i+'/function.data.tar.bzip2'):
                print('extracting', i+'/function.data')
                tar = tarfile.open(i+'/function.data.tar.bzip2', "r:bz2")
                tar.extractall(i)
                tar.close()
        print()

    # to get function.data for the new file:
    #  - get input.nn from DB1
    #  - get input.data from the current input data
    #  - take the submitskirt and run it. this is fast (< 1min);
    if not os.path.isfile(DB2_path+'/function.data'):
        gfdf = get_function_data_folder = DB2_path+"/get_function_data"
        if os.path.isfile(gfdf+"/function.data"):
            my.cp(gfdf+"/function.data",DB2_path+'/function.data')
        else:
            # setup the job to create function data for the new structures
            if not os.path.isfile(DB1_path+'/input.nn'):
                sys.exit("Need "+DB1_path+'/input.nn')
            scripts = my.scripts()
            submitscript = scripts+"/runner_scripts/submit_scaling_debug.sh"
            if not os.path.isfile(submitscript):
                sys.exit("Need submitscript "+submitscript)
            print()
            print("making "+DB2_path+"/get_function_data")
            my.mkdir(DB2_path+"/get_function_data")
            my.cp(DB1_path+'/input.nn',DB2_path+"/get_function_data/input.nn")
            my.cp(DB2_path+'/input.data',DB2_path+"/get_function_data/input.data")
            my.cp(submitscript,DB2_path+"/get_function_data/submit_scaling_debug.sh")
            my.submitjob(submitdebug=True,jobdir=DB2_path+"/get_function_data",submitskript="submit_scaling_debug.sh")
            my.create_READMEtxt(DB2_path+"/get_function_data/",add="# It it necessary to create the function data for the new structures for fps")
            my.create_READMEtxt(os.getcwd()                   ,add=["#","# Restart this skript once "+DB2_path+"/get_function_data/function.data exists!"])
            sys.exit("submitted job to get function data to debug que; this typically takes only ~50 sec; once finished, restart this job again")

    for i in [ DB1_path,DB2_path]:
        if not os.path.isfile(i+'/function.data'):
            sys.exit(i+"/function.data does not exist! (do nnp-scaling to get it)")
        else: print(i+"/function.data exists")
    print()


    #sys.exit()
    # to create DB1
    if not os.path.isfile(DB1_path+'/function_average.data'):
        print('reading:',DB1_path+'/input.data')
        nat_per_frame_5000 = getnat_per_frame(DB1_path+'/input.data')
        print('read in:',DB1_path+'/input.data',"nat_per_frame:",nat_per_frame_5000)
        print('createing (DB1)',DB1_path+'/function.data')
        DB1 = create_average_SF(DB1_path+'/function.data', nsyms, nat_per_frame_5000,outfile=DB1_path+'/function_average.data')
        print('createing (DB1)',DB1_path+'/function.data DONE!')
    else:
        print('loading ... DB1 from ',DB1_path+'/function_average.data')
        DB1 = np.loadtxt(DB1_path+'/function_average.data')

    print()
    if not os.path.isfile(DB2_path+'/function_average.data'):
        print('reading:',DB2_path+'/input.data')
        nat_per_frame_2509 = getnat_per_frame(DB2_path+'/input.data')
        print(nat_per_frame_2509)
        print('createing (DB2)',DB2_path+'/function.data')
        DB2 = create_average_SF(DB2_path+'/function.data', nsyms, nat_per_frame_2509,outfile=DB2_path+'/function_average.data')
        print('createing (DB2)',DB2_path+'/function.data DONE!')
    else:
        print('loading ... DB2 from ',DB2_path+'/function_average.data')
        DB2 = np.loadtxt(DB2_path+'/function_average.data')

    #sys.exit()
    DB1len = DB1.shape[0]
    DB2len = DB2.shape[0]
    every = 500
    if DB2len < 1000:
        every = 100
    if DB2len < 100:
        every = 10


    if not os.path.isdir("fps"):
        os.mkdir('fps')

    print()
    print()
    # one could in principle check if any of the structures in DB2 are repetitions.
    if not os.path.isfile('fps/dist_vec_from_DB1.dat'):
        print('-------------------------------------------------------------------------------')
        print('making fps/dist_vec.dat to find the structure in DB2 which is furthest from DB1')
        print('-------------------------------------------------------------------------------')

        dist_vec = np.full((DB2len), np.inf)  # 2509
        # this could be easily parallelized ...
        for i in range(DB2len):  # 0 ... 2508
            for j in range(DB1len):  # 0 ... 2508
                dist = salg.norm(DB2[i]-DB1[j])
                if dist < dist_vec[i]:
                    dist_vec[i] = dist
            if i in np.arange(0,DB2len,every):
                print('i',i,'/',DB2.shape[0],dist_vec[i])
        #print('i',i,dist,dist_vec[0])
        #print(dist_vec)
        np.savetxt('fps/dist_vec_from_DB1.dat',dist_vec)  # this is the distance of every structure in DB2 to all structures in DB1
    else:
        print('reading fps/dist_vec_from_DB1.dat')
        dist_vec = np.loadtxt('fps/dist_vec_from_DB1.dat')

    print()
    if not os.path.isfile('fps/dist_vec_fin.dat'):
        print('----------------------------------------------------------')
        print('creating dist_vec_fin.dat which holds the distances of DB2')
        print('----------------------------------------------------------')
        argmax = np.argmax(dist_vec) # this has the largest distance to the previous 5000 struct.
        distmax = dist_vec[argmax]
        # from here on we dont really need dis_vec anymore! Only DB2 distanes will be checked.
        print('sarting loop....')
        print(0,'distmax,argmax :',distmax,argmax) #,np.where(dist_vec == 0)[0])
        #dist_vec[0] = distmax
        length = DB2.shape[0]
        dist_vec_new      = np.zeros((length,3))  # 2509
        dist_vec_new[0,0] = 0
        dist_vec_new[0,1] = distmax
        dist_vec_new[0,2] = argmax

        for j in np.arange(1,DB2len): # geht ueber alle eingraege von DB2, 2509 eintraege, diese sind von interesse.
            for i in range(DB2len): # geht ueber alle eingraege von DB2, 2509 eintraege, diese sind von interesse.
                new_dist = salg.norm(DB2[i]-DB2[argmax])
                if new_dist < dist_vec[i]:
                    dist_vec[i] = new_dist
                    #dist_vec_new[j,1] = new_dist
            #argmax = np.argmax(dist_vec_to2509)
            argmax = np.argmax(dist_vec)
            distmax = dist_vec[argmax]
            #dist_vec[argmax] = 0
            if j in np.arange(0,DB2len,every):
                print(j,'/',DB2len,'distmax,argmax :',distmax,argmax) #,np.where(dist_vec == 0)[0])
            dist_vec_new[j,0] = int(j)
            dist_vec_new[j,1] = distmax
            dist_vec_new[j,2] = argmax
            #np.savetxt('dist_vec_'+str(j),dist_vec)
        np.savetxt('fps/dist_vec_fin.dat',dist_vec_new)
    else:
        print('reading fps/dist_vec_fin.dat')
        dist_vec_new = np.loadtxt('fps/dist_vec_fin.dat')


    print()
    frames = ase_read(DB2_path+'/input.data',':',format='runner')
    fo = 'fps/input.fps'+str(DB2len)+'.data'
    print("----------------------------------------------------")
    print('saving DB2 ',fo)
    print("----------------------------------------------------")
    if structures_upperlim > 0:
        fo = 'fps/input.fps'+str(structures_upperlim)+'.data'
    for idx,i in enumerate(dist_vec_new[:,2].astype(int)):
        #print('i',i)
        #print(frames[int(i)].get_chemical_formula())
        ase_write(fo,frames[i],format='runner',append=True)
        if structures_upperlim > 0 and idx == structures_upperlim:
            break

    my.create_READMEtxt(os.getcwd(),add="Do: xmgrace fps/dist_vec_fin.dat in a log/log plot; Then grep the first structures of interest from fps/input.fps256.data ")
    return


def create_average_SF(datafile, nsyms, nat_per_frame,outfile=False):
    #!TODO nsyms should be discarded, like also nat_per_frame. Both should be read within the function, the only arguments passed
    # should be datafile, logfile, and nlandmarks
    """ This function computes the distance between the frames
    First it reduces the amount of data from function.data in a single array for every frame
    Then it calculates the furthest frames from the ones already selected until the requested number of landmarks is hit.
    datafile --> function.data
    nsyms --> a dictionary that contains the elements and the respective number of SF i.e. nsyms = {"Al":64, "Si":64}
    nat_per_frame --> a dictionary that contains the species and the corresponding number of atoms in each frame,
        i.e. nat_per_frame = {"Al" : np.array([10,20,12]), "Si" : np.array([2,10,12])} --> 3 frames; first contains 10 Al and 2 Si, second 20 Al and 10 Si, third 12 Al and 12 Si
    """

    avg_sf = dict.fromkeys(nsyms.keys())

    print('Preprocessing for landmark selection')
    # reads in symmetry functions from the preprocessed data and computes an "average fingerprint" for the atoms of each specie
    totsym = 0
    natoms = None
    for element in nsyms:
        nsym = nsyms[element]
        totsym += nsym

        nframes = len(nat_per_frame[element])
        #print('nframes',nframes)
        nlandmarks = nframes
        if nlandmarks > nframes:
            print("You requested {:d} landmarks and there are {:d} frames. You can use the whole input.data or re-run CurSel requesting fewer landmarks".format(nlandmarks, nframes))
            return range(nframes)
        if natoms is None:
            natoms = np.array(nat_per_frame[element])
        else: natoms += nat_per_frame[element]
        xmat = ReadFuncdata(datafile, element)
        avg_sf[element] = np.zeros((len(nat_per_frame[element]),nsym))

        iframe = 0
        past_row = 0
        for nrow in nat_per_frame[element]:
            if nrow == 0: # this element is missing in this frame, so we just add a block of zeros
                avg_sf[element][iframe,:]=np.zeros(nsym)
            else:
                data = xmat[past_row:past_row+nrow]
                past_row += nrow
                avg_sf[element][iframe,:]=np.sum(data,axis=0)
            iframe += 1

    # now collates the descriptors for each element, and normalizes properly so we can use the distance between
    # the compounded fingerprints as an indicator of the overall nature of each frame
    frame_sf = np.zeros((nframes, totsym))

    ksym = 0
    for element in nsyms:
        nsym = nsyms[element]
        for k in range(nframes):
            frame_sf[k,ksym:ksym+nsym] = avg_sf[element][k,:] / natoms[k]
        ksym += nsym
    np.savetxt(outfile,frame_sf)
    return frame_sf

def getnat_per_frame(infile):
    frames = ase_read(infile,':',format='runner')
    print('len',len(frames))
    mg = np.zeros(len(frames)).astype(int)
    al = np.zeros(len(frames)).astype(int)
    si = np.zeros(len(frames)).astype(int)
    for idx,i in enumerate(frames):
        #print()
        #print(idx,i.get_number_of_atoms())
        #print(idx,i.get_chemical_formula())

        #chs = i.get_chemical_symbols()
        #elements = set(i.get_chemical_symbols())
        #string = i.get_chemical_formula()

        mg[idx] = int(np.sum(i.numbers == 12))
        #print('mg',mg[idx],type(mg[idx]))
        al[idx] = int(np.sum(i.numbers == 13))
        si[idx] = int(np.sum(i.numbers == 14))
    return {"Mg":mg,"Al":al,"Si":si}




if __name__ == "__main__":
    base = "/scratch/glensk/n2p2_get_function_data_and_make_fps/"
    #DB1_path    = base + "5000_struct"
    #DB2_path    = base + "2509_struct"
    #DB2_path    = base + "5020_struct"
    #nsym = {"Mg":64, "Al":64, "Si":64}
    make_fps()
