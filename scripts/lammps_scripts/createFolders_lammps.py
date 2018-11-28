#!/usr/bin/env python
from __future__ import print_function
import os,sys
from shutil import copyfile
import click
import myutils as my
import glob

CONTEXT_SETTINGS = my.get_click_defaults()
@click.command(context_settings=CONTEXT_SETTINGS)

@click.argument('foldername')
@click.option('--calc_type','-t', type=click.Choice(['static','geopt','ti']),default='static')
@click.option('--pos','-pos', type=str, default="positions.lmp")

@click.option('--pot','-p', type=click.Choice(['eam','runner','n2p2']),default='eam')
@click.option('--pot_id','-pid', type=str, default="Al-LEA.eam.alloy",help="this is the file or the corresponding folder of the potential")
@click.option('--test_static', default=False)

def createFolder_lammps(foldername,pos,calc_type,pot,pot_id,test_static):
    ''' calc_type any of [static|geopt|ti] '''
    if os.path.isdir(foldername):
        sys.exit(foldername+' does already exist; Exit.')
    my.mkdir(foldername)
    print('foldername:',foldername)
    print('positions :',pos)
    print('positions :',os.path.abspath(pos))


    scripts = my.scripts()

    potpath = get_pot(pot,pot_id,verbose=False)

    potlines = potpath_to_lammps_lines(pot,potpath,verbose=True)
    print('potlines',potlines)

    # make in.lmp
    if test_static:
        copyfile(scripts+'/lammps_scripts/positions_dummyfiles/2x2x2sc_fcc_cubic_vac',foldername+'/positions.lmp')
    get_inlmp(foldername,os.path.abspath(pos),potlines)

    my.create_READMEtxt(os.getcwd())
    return


def potpath_to_lammps_lines(pot,potpath,verbose=False):
    if verbose:
        print('pot    :',pot)
        print('potpath:',potpath)
    if pot == 'eam':
        return ["pair_style eam/alloy\n","pair_coeff * * "+potpath+" "+"Al"]
    if pot == 'runner':
        a = 'variable runnerDir string "'+potpath+'"'
        b = "pair_style runner dir ${runnerDir} showewsum 1 showew yes resetew no maxew 1000000"
        c = "variable runnerCutoff equal 10.0"
        d = "pair_coeff * * ${runnerCutoff}"
        return [a+"\n",b+"\n",c+"\n",d+"\n"]
    else:
        return False


def get_pot(type='eam',id="",verbose=False):
    potout = False
    listout_basename = []
    listout_fullpath = []
    scripts = my.scripts()
    pot_folder = scripts + '/potentials/lammps_'+type+'/'
    files = glob.glob(pot_folder+'/*')
    for i in files:
        if verbose:
            print(os.path.basename(i))
        listout_basename.append(os.path.basename(i))
        listout_fullpath.append(i)

    for i in files:
        if os.path.basename(i) == id:
            potout = i

    if verbose:
        print('potout',potout)
    if potout == False:
        print('type',type)
        print('id',id)
        print('pot_folder',pot_folder)
        print("potout",potout)
        sys.exit('potential not found')

    return potout


def get_inlmp(folder,positions_file,potlines,static=True):
    with open(folder+'/in.lmp', "w") as f:
        # stuff which is always added
        f.write("\n")
        f.write("clear\n")
        f.write('variable dump_file string "trj_lammps"\n')
        f.write("\n")
        f.write("units metal\n")
        f.write("boundary  p p p\n")
        f.write("atom_modify sort 0 0.0\n")
        f.write("\n")
        f.write("read_data "+str(positions_file)+"\n")
        f.write("\n")

        # potential
        f.write("## potential\n")
        for i in potlines:
            f.write(i)
        f.write("\n")

        # fixes define how atoms move

        # run
        #f.write("## run\n")
        #f.write("dump dump_all all custom 50 trj_lammps id type x y z fx fy fz\n")
        #f.write("thermo_style custom step cpu pe etotal vol\n")
        #f.write("thermo 10\n")
        f.write("\n")
        f.write("min_style cg\n")
        f.write("#minimize 1.0e-9 1.0e-10 1000 1000\n")
        f.write("\n")
        f.write("run 0\n")
    return



if __name__ == "__main__":
    createFolder_lammps()
