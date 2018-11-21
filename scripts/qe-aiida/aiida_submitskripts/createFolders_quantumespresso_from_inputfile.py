#!/usr/bin/env python

import os,sys,copy
import numpy as np
import ase
from shutil import copyfile
import click

# from scripts folder
import myutils as my

# show default values in click
orig_init = click.core.Option.__init__
def new_init(self, *args, **kwargs):
    orig_init(self, *args, **kwargs)
    self.show_default = True
click.core.Option.__init__ = new_init

# get help also with -h
CONTEXT_SETTINGS = dict(help_option_names=['-h', '--help'])
@click.command(context_settings=CONTEXT_SETTINGS)

@click.argument('ase_positionsfiles',required=True,nargs=-1)
@click.option('-calculation','--calculation',default='scf',type=click.Choice(['scf', 'relax']),help="qe calculation type")
@click.option('-kp','--kpoint_density',default=80,type=int,help="k-point density")
@click.option('-forc_conv_thr','--forc_conv_thr',default='2.0d-5',help="qe forc_conv_thr")
@click.option('-espresso_pseudo', envvar='ESPRESSO_PSEUDO',help='environment variable $ESPRESSO_PSEUDO (folder with qe potentials)')
@click.option('-submit/-no-submit', default=False)
@click.option('-submithost','--submithost',default='daint',type=click.Choice(['fidis', 'daint','cosmopc','mac']),help="submithost for the job/jobs")
@click.option('-submitdebug/-no-submitdebug', default=False)
@click.option('-v', '--verbose', count=True)


def adapt_ase_qe_file(
        ase_positionsfiles,
        calculation,
        kpoint_density,
        forc_conv_thr,
        submit,
        submitdebug,
        submithost,
        espresso_pseudo,
        verbose
        ):
    '''
    Create folder from given ase_positionsfiles.
    '''

    if verbose > 0:
        print('verbose',verbose)

    # checks
    submit_file = my.scripts()+"/qe-aiida/aiida_submitskripts/submit-aiida-qe.sh"
    settings_file = my.scripts()+"/qe-aiida/aiida_submitskripts/aiida.in.top"
    my.check_isfile_or_isfiles([submit_file,settings_file],["submit_file","settings_file"])
    for positionsfile_ase in ase_positionsfiles:
        my.check_isfile_or_isfiles([positionsfile_ase],["positionsfile_ase"])
    my.check_isdir_or_isdirs([espresso_pseudo])

    # get the settings (top of the file)
    f = open(settings_file,"r")
    settings = f.readlines()
    f.close()

    ######################################################
    # go over every inputfile
    ######################################################
    for idx_folder,positionsfile_ase in enumerate(ase_positionsfiles):
        # get the structures file
        f = open(positionsfile_ase,"r")
        positions = f.readlines()
        f.close()

        # read in the ase structure
        ase_structure = ase.io.read(positionsfile_ase, format="espresso-in")

        # k-points
        kpoint_density = 80
        kmesh = my.get_kmesh_size_daniel(ase_structure, kpoint_density)
        if verbose:
            print('kmesh',kmesh)

        #########################
        # adapt the positions
        #########################
        listdelete = []
        listdelete_append = True
        adaptelements = False

        idxkpoints = -1
        for idx,line in enumerate(positions):
            if verbose > 1:
                if idx < 30: print('positions idx',idx,'line:'+line.rstrip()+":",listdelete_append,adaptelements)
            if line[:8] == "   ntyp ":
                ntype_line = copy.copy(positions[idx])
            if line[:8] == "   nat  ":
                nat_line = copy.copy(positions[idx])
            if line[:9] == "K_POINTS ":
                idxkpoints = idx+1
                positions[idx] = "K_POINTS automatic\n"
            if idx == idxkpoints:
                positions[idx] = str(int(kmesh[0]))+" "+str(int(kmesh[1]))+" "+str(int(kmesh[2]))+' 0 0 0\n'

            if listdelete_append: listdelete.append(idx)

            # where to start adapting elements
            if line[:14] == "ATOMIC_SPECIES":
                listdelete_append = False
                adaptelements = True
                continue

            # where to stop adapting element lines
            if listdelete_append == False and line.rstrip() == "" and adaptelements == True:
                adaptelements = False

            if verbose > 1:
                if idx < 30: print('positions idy',idx,'line:'+line.rstrip()+":",listdelete_append,adaptelements)

            # adapt the names of the elements
            if adaptelements == True:
                element_str = positions[idx].split()[0]
                potential=my.qe_full_get_path_to_potential(element_str,espresso_pseudo)
                pot = os.path.basename(potential[0])
                if verbose > 1:
                    print()
                    print('potential',potential)

                my_list = positions[idx].split()
                my_list[2] = pot
                positions[idx] = " ".join(my_list)
                positions[idx] = positions[idx]+'\n'


        number_electrons = my.qe_get_numelectrons(ase_structure,espresso_pseudo)
        number_bands = int(max(20, 0.75*number_electrons))

        if verbose > 0:
            print('number_electrons:',number_electrons)
            print('number_bands    :',number_bands)
            print('listdelete:',listdelete)
            print('ntype_line:',ntype_line)
            print('nat_line:',nat_line)

        for i in np.array(listdelete)[::-1]:
            del positions[i]

        #########################
        # adapt the settings
        #########################
        for idx,line in enumerate(settings):
            if line[:15] == "  pseudo_dir = ":
                settings[idx] = "  pseudo_dir = '"+espresso_pseudo+"'\n"
            if line[:15] == "  calculation =":
                settings[idx] = "  calculation = '"+calculation+"'\n"
            if line[:15] == "  forc_conv_thr":
                settings[idx] = "  forc_conv_thr = "+forc_conv_thr+"\n"
            if line[:8] == "   ntyp ":
                settings[idx] = ntype_line
            if line[:8] == "   nat  ":
                settings[idx] = nat_line


        ###############################################
        # create a directory (if more than one file)
        ###############################################
        dirname = "pos_"+str(idx_folder)
        my.mkdir(dirname)
        copyfile(submit_file,dirname+"/submit.sh")

        # write aiida.in
        f = open(dirname+'/aiida.in',"w")
        f.writelines(settings)
        f.writelines(positions)
        f.close()
        print('written '+dirname+'/aiida.in')

        ###############################################
        # submit
        ###############################################
        my.submitjob(submit=submit,submitdebug=submitdebug,jobdir=dirname,submitskript="submit.sh")

        my.create_READMEtxt(os.getcwd())
    return

if __name__ == "__main__":
    adapt_ase_qe_file()
