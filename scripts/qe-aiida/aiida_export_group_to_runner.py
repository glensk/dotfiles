#!/usr/bin/env python
from __future__ import print_function
import aiida
aiida.load_dbenv()
from aiida.orm.querybuilder import QueryBuilder
from aiida.orm import Group, WorkCalculation
from aiida.orm.utils import load_node
import click
import myutils,os
import sys
from ase.io import read

# show default values in click
orig_init = click.core.Option.__init__


def new_init(self, *args, **kwargs):
    orig_init(self, *args, **kwargs)
    self.show_default = True


click.core.Option.__init__ = new_init

# get help also with -h
CONTEXT_SETTINGS = dict(help_option_names=['-h', '--help'])


@click.command(context_settings=CONTEXT_SETTINGS)
@click.option('-wg', '--work_group', required=True, prompt=True,
              type=str, help="verid group to be converted to runner format")
@click.option('-f', '--filename', required=False, default=False,
              type=str, help="filename for outputfile, default = work_group+\".input.data\"")

def createjob(work_group,filename):
    ''' e.g.
    ./aiida_export_group_to_runner.py -wg kmc_1000K_4
    ./aiida_export_group_to_runner.py -wg Al6xxxDB_passingsubset works
    '''
    print('work_group:', work_group)
    print('filename  :',filename)
    print()
    qb = QueryBuilder()
    qb.append(Group, filters={'name': work_group}, tag='g')
    qb.append(WorkCalculation, tag='job', member_of='g')
    all_works = [x[0] for x in qb.all()]  # all workchains

    def get_workcalc_runnerdata(worknode):
        ase_structure = worknode.inp.structure.get_ase()
        energy = worknode.out.output_parameters.get_attrs()['energy']  # units?
        forces = worknode.out.output_array.get_array('forces')  # units?
        path = worknode.out.CALL.out.retrieved.get_abs_path()
        return ase_structure, energy, forces, worknode.uuid, path

    angstrom_to_bohrradius = 1.8897261
    eV_to_Hartree = 0.036749325
    eV_per_angstrom_to_hartree_per_bohrradius = 0.019446905

    if filename == False:
        fileOut = open("aiida_export_group_"+work_group+".data", "w")
    else:
        fileOut = open(filename, "w")
    myutils.create_READMEtxt(os.getcwd())

    if False:
        # from verdi calculation list -a -p1
        failed = [102402,102497,102500,102503,102523,102530,102537,102588]
        worked_part = [102202,102217,102230,102261,102272,102298,102313,102323]

        for idx, workchain in enumerate(all_works):
            print()
            #print('idx all:',idx,workchain)
            uuid = workchain.uuid
            pk = workchain.pk

            print('uuid:',uuid) #,workchain.pk)
            if pk in failed:
                print('pk  :',pk,"FAILED") #,workchain.pk)
            else:
                print('pk  :',pk,"WORKED") #,workchain.pk)


            calc = load_node(pk)
            uuidcheck = calc.uuid
            abspath = calc.get_abs_path()
            abspath2=""
            try:
                abspath2= calc.out.retrieved.get_abs_path()
            except AttributeError:
                pass
            #print('path:',abspath) #,workchain.pk)
            #print('path2',abspath2) #,workchain.pk)
            if uuid != uuidcheck:
                sys.exit('uuid prob')

            #try:
            #    calc.inp.CALL.inp.options.get_dict()
            #except AttributeError:
            #    print("Error here")
            #    continue
    if False:
        for idx in failed:
            calc = load_node(idx)
            pathinput= calc.get_abs_path()
            path = calc.out.retrieved.get_abs_path()
            out = myutils.grep(pathinput+"/raw_input/_aiidasubmit.sh","#SBATCH --nodes=")
            atoms= myutils.grep(pathinput+"/raw_input/aiida.in","nat =")
            print("Failed path input:",pathinput)
            #print("Failed           :",path)
            print("out FAILED       :",out,"atoms",atoms)
            dd=calc.inp.CALL.inp.options.get_dict()
            print('dd',dd)
        print()
        print('----------')

        for idx in worked_part:
            calc = load_node(idx)
            pathinput =  calc.get_abs_path()
            path = calc.out.retrieved.get_abs_path()
            out = myutils.grep(pathinput+"/raw_input/_aiidasubmit.sh","#SBATCH --nodes=")
            atoms= myutils.grep(pathinput+"/raw_input/aiida.in","nat =")

            print("worked path input :",pathinput)
            print("worked path output:",path)
            print("out WORKED        :",out,"atoms",atoms)
            dd=calc.inp.CALL.inp.options.get_dict()
            print('dd',dd)

    for idx, workchain in enumerate(all_works):
        try:
            ase_structure, energy, forces, uuid, path = get_workcalc_runnerdata(workchain)
        except AttributeError:
            print('WARNING: idx',idx,"failed for some reason!!! (will be excluded)")
            continue
        print(idx, "ene (eV)", energy, uuid, path)
        forces = forces[0]

        work_uuid = workchain.uuid
        fileOut.write("begin\ncomment uuid: {}\n".format(work_uuid))

        # write the cell
        cell = ase_structure.cell*angstrom_to_bohrradius
        for idx_cell, i in enumerate(cell):
            fileOut.write("lattice %.10f %.10f %.10f\n" %
                          (cell[idx_cell][0], cell[idx_cell][1], cell[idx_cell][2]))

        #  write the positions
        nr_of_atoms = ase_structure.positions.shape[0]
        for idx_pos in range(nr_of_atoms):
            atCor = ase_structure.positions[idx_pos]*angstrom_to_bohrradius
            atFor = forces[idx_pos]*eV_per_angstrom_to_hartree_per_bohrradius
            element = ase_structure.get_chemical_symbols()[idx_pos]
            fileOut.write("atom   %.6f    %.6f   %.6f %s  0.0   0.0  %.10f  %.10f  %.10f\n" %
                          (atCor[0], atCor[1], atCor[2],
                           element,
                           atFor[0], atFor[1], atFor[2]))
        fileOut.write("energy %.15f\n" % (energy*eV_to_Hartree))
        fileOut.write("charge 0\nend\n")

    fileOut.close()
    return


if __name__ == "__main__":
    createjob()
