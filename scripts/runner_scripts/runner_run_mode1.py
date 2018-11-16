#!/usr/bin/env python

import os,sys
python_version = sys.version_info[0]
if python_version < 3:
    sys.exit('Your python environment uses a python < 3; Exit;')

import click
import glob
from shutil import copyfile
import numpy as np

#from subprocess import call,Popen,check_output,PIPE,run
from subprocess import call,run,check_output


# from scripts folder
import kmc_createjob

# show default values in click
orig_init = click.core.Option.__init__
def new_init(self, *args, **kwargs):
    orig_init(self, *args, **kwargs)
    self.show_default = True
click.core.Option.__init__ = new_init


# get help also with -h
CONTEXT_SETTINGS = dict(help_option_names=['-h', '--help'])
@click.command(context_settings=CONTEXT_SETTINGS)

@click.option('-fi','--filename_in',required=False,default="input.data",help="file containing runner frames")
@click.option('-runner_exec', envvar='runner_exec',help='path to RuNNer executable (or environment variable $runner_exec)')
@click.option('-scripts', envvar='scripts',help='environment variable $scripts (can alse be set here)')
@click.option('-t','--test/--no-test', default=False)

def make_fps(filename_in,runner_exec,scripts,test):
    '''
    This scipt ....
    '''
    ######################################
    # check weather input.data exists
    ######################################
    if not os.path.isfile(filename_in):
        sys.exit('file '+filename_in+" does not exist! Exit")
    file_input_runner          = scripts + "/runner_scripts/inputAlMgSi.nn"

    ######################################
    # check if necessary inputfiles and runner executable is defined
    ######################################
    kmc_createjob.check_isfile([file_input_runner],["inputAlMgSi.nn"])
    if not test:kmc_createjob.check_isfile([runner_exec],["runner_exec"],environment=True)

    ######################################
    # make symfun.output
    ######################################
    symfun = 'symfun.output'
    if os.path.isfile(symfun):
        sys.exit(symfun+' file does already exist! Exit')

    f = open(symfun,"wb")
    call(["symfun_gen.py","-e","Al,Si,Mg","-c","12","-n","10"],stdout=f)
    call(["symfun_gen.py","-e","Al,Si,Mg","-c","16","-n","10"],stdout=f)
    call(["symfun_gen.py","-e","Al,Si,Mg","-c","20","-n","10"],stdout=f)
    call(["symfun_gen.py","-e","Al,Si,Mg","-c","16","-n","4"],stdout=f)
    f.close()

    run(["sort "+symfun+" | uniq > tmp"],shell=True)
    call(["mv","tmp",symfun])
    print("written",symfun)

    ############################################
    # get and adapt input.nn for runner_mode 1
    ############################################
    get_input_runner(file_input_runner,"input.nn",symfun_delete=True,symfun_file=symfun,\
            test_fraction= 0,
            runner_mode= 1,
            number_of_elements= 3,
            elements = "Al Mg Si",
            test_input_data=filename_in)

    ######################################
    # write README
    ######################################
    kmc_createjob.create_READMEtxt(os.getcwd(),add = "# using Runner executable: "+str(runner_exec))

    ######################################
    # run runner (whoever not in testmode)
    ######################################
    if test == True:
        sys.exit('Not starting RuNNer in testmode')
    else:
        print("now starting RuNNer ...")
        run([runner_exec+" > logfile_mode1.1&"],shell=True)
    return


def file_len(fname):
    i = 0
    with open(fname) as f:
        for i, l in enumerate(f,1):
            pass
    return i


def get_input_runner(template,filename,symfun_delete=True,symfun_file=False,
        test_fraction=0,runner_mode=1,number_of_elements=3,elements="Al Mg Si",test_input_data=True):
    if test_input_data:
        len = file_len(test_input_data)
        #print('len',len)
        if len <= 3:sys.exit('file '+test_input_data+' seems too short! Exit;')

    # read in the runner.in template
    f = open(template,"r")
    lines = f.readlines()
    f.close()

    if symfun_delete == True:
        listdelete = []

        # delete the old symmetry functioins
        for idx,line in enumerate(lines):
            #print()
            #print('idx',idx,line)
            #print('idk',idx,line[:18])
            if line[:18] == "symfunction_short ":
                listdelete.append(idx)
            if line[:24] == "# symfunctions for type ":
                listdelete.append(idx)

        for i in np.array(listdelete)[::-1]:
            del lines[i]

        # insert the new symmetry functions
        if symfun_file != False:
            s = open(symfun_file,"r")
            sym = s.readlines()
            s.close()
            for idj,symline in enumerate(np.array(sym)[::-1]):
                #print('sl',symline)
                lines.insert(listdelete[0],symline)

        # set other options
        print('test_fraction        :',test_fraction)
        print('runner_mode          :',runner_mode)
        print('number_of_elements   :',number_of_elements)
        print('elements             :',elements)
        for idx,line in enumerate(lines):
            if line[:14] == "test_fraction ":
                lines[idx] = "test_fraction "+str(test_fraction)+"\n"
            if line[:12] == "runner_mode ":
                lines[idx] = "runner_mode "+str(runner_mode)+"\n"
            if line[:19] == "number_of_elements ":
                lines[idx] = "number_of_elements "+str(number_of_elements)+"\n"
            if line[:9] == "elements ": lines[idx] = "elements "+str(elements)+"\n"


        # write the file
        f = open(filename,"w")
        f.writelines(lines)
        f.close()
        print('written '+filename)
        return


if __name__ == "__main__":
    make_fps()
