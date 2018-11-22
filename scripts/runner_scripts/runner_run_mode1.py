#!/usr/bin/env python

import os,sys
import click
import glob
from shutil import copyfile
import numpy as np

#from subprocess import call,Popen,check_output,PIPE,run
from subprocess import call,run,check_output


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

@click.option('-fi','--filename_in',required=False,default="input.data",help="file containing runner frames (positions)")
@click.option('-m','--mode',required=False,default=1,type=int,help="runner mode")
@click.option('-tf','--test_fraction',required=False,default=0,type=float,help="runner test fraction")
@click.option('-t','--test/--no-test', default=False)

def runner_run_mode_1(filename_in,mode,test_fraction,test):
    '''
    This scipt gets some symmetry functions from scratch for making fps:
    runner_run_mode1.py -m 1 -tf 0  (first attempt to get symmetry functions if none are available, uses (runner) mode = 1 and (runner) test_fraction = 0)
    '''
    ######################################
    # check weather necessary runner input data exists
    ######################################
    scripts = my.scripts()
    if not test: runner_exec = my.runner_exec()
    file_input_runner          = scripts + "/runner_scripts/inputAlMgSi.nn"
    my.check_isfile_or_isfiles([filename_in,file_input_runner],["filename_in","inputAlMgSi.nn"])


    ######################################
    # make symfun.output
    ######################################
    if mode == 1 and testfraction == 0:
        symfun_file = 'symfun.output'
        if os.path.isfile(symfun_file):
            sys.exit(symfun_file+' file does already exist! Exit')

        f = open(symfun_file,"wb")
        call(["symfun_gen.py","-e","Al,Si,Mg","-c","12","-n","10"],stdout=f)
        call(["symfun_gen.py","-e","Al,Si,Mg","-c","16","-n","10"],stdout=f)
        call(["symfun_gen.py","-e","Al,Si,Mg","-c","20","-n","10"],stdout=f)
        call(["symfun_gen.py","-e","Al,Si,Mg","-c","16","-n","4"],stdout=f)
        f.close()

        run(["sort "+symfun_file+" | uniq > tmp"],shell=True)
        call(["mv","tmp",symfun_file])
        print("written",symfun_file)

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
        call([runner_exec+" > logfile_mode1.1&"],shell=True)  # call is interactive, run is in the background
    return


if __name__ == "__main__":
    runner_run_mode_1()
