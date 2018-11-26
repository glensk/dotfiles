#!/usr/bin/env python
from __future__ import print_function
import os,sys
import click
import glob
from shutil import copyfile
import numpy as np
from subprocess import call,check_output
import myutils as my

CONTEXT_SETTINGS = my.get_click_defaults()
@click.command(context_settings=CONTEXT_SETTINGS)

@click.option('-fi','--filename_in',required=False,default="input.data",help="file containing runner frames (positions)")
@click.option('-m','--mode',required=False,default=1,type=int,help="runner mode")
@click.option('-tf','--test_fraction',required=False,default=0,type=float,help="runner test fraction")
@click.option('-t','--test/--no-test', default=False)

def runner_run_mode_1(filename_in,mode,test_fraction,test):
    '''
    This scipt gets some symmetry functions from scratch for making fps:
    runner_run_mode1.py -m 1 -tf 0  (first attempt to get symmetry functions if none are available, uses (runner) mode = 1 and (runner) test_fraction = 0)

    runner mode 1 --- evaluation of the symmetry functions
    runner mode 2 --- creation of the potential

    important runner files (output):
    - function.data (== representation of the stuctrues (~positions)
      [either all or test_fraction] for a chosen set of symmetry functions)
    -

    steps:
    1.1 get symmetry functions (which are derived by CurSel.py)
    1.2 get the correct representation/function.date
        (with CurSel.py derived symmetry functions)
    2.  get the potential
    '''
    ######################################
    # check weather necessary runner input data exists
    ######################################
    scripts = my.scripts()
    runner_exec = my.runner_exec(test=test)
    if not test: runner_exec = my.runner_exec()
    file_input_runner          = scripts + "/runner_scripts/inputAlMgSi.nn"
    my.check_isfile_or_isfiles([filename_in,file_input_runner],["filename_in","inputAlMgSi.nn"])

    ######################################
    # make symfun.output for step 1.1
    ######################################
    if mode == 1 and test_fraction == 0:
        symfun_file = 'symfun.output'
        if os.path.isfile(symfun_file):
            sys.exit(symfun_file+' file does already exist! Exit')

        f = open(symfun_file,"wb")
        # this would be 2250 symmetry functions
        call(["symfun_gen.py","-e","Al,Si,Mg","-c","12","-n","10"],stdout=f) # 729 SF
        call(["symfun_gen.py","-e","Al,Si,Mg","-c","16","-n","10"],stdout=f) # 729 SF
        call(["symfun_gen.py","-e","Al,Si,Mg","-c","20","-n","10"],stdout=f) # 729 SF
        call(["symfun_gen.py","-e","Al,Si,Mg","-c","16","-n","4"],stdout=f) # 297 SF
        f.close()

        #run(["sort "+symfun_file+" | uniq > tmp"],shell=True)
        call(["sort "+symfun_file+" | uniq > tmp"],shell=True)
        call(["mv","tmp",symfun_file])
        print("written",symfun_file)

    ###############################################################
    # if more than 300 structures: select structures
    # however: for 2509 structures it fook less than a day.
    ###############################################################

    ############################################
    # get and adapt input.nn for runner_mode 1
    ############################################
    my.get_inputfile_runner(file_input_runner,"input.nn",symfun_old_delete=True,symfun_file=symfun_file,\
            test_fraction= 0,
            runner_mode= mode,
            number_of_elements= 3,
            elements = "Al Mg Si",
            test_input_data=filename_in)

    ######################################
    # write README
    ######################################
    my.create_READMEtxt(os.getcwd(),add = ["# using Runner executable: "+str(runner_exec),"# running on "+my.hostname()])

    ######################################
    # run runner (whoever not in testmode)
    ######################################
    if test == True:
        sys.exit('Not starting RuNNer in testmode')
    else:
        print("now starting RuNNer ... (in the background)")
        call([runner_exec+" > logfile_mode1.1&"],shell=True)  # call is interactive, run is in the background but since this rediects the output it works as run
    return


if __name__ == "__main__":
    runner_run_mode_1()
