#!/usr/bin/env python
from __future__ import print_function
import aiida,os
aiida.try_load_dbenv()
from aiida.orm import Node
from aiida.orm.querybuilder import QueryBuilder
from aiida.orm import Calculation, Group, WorkCalculation
from aiida.orm.data.structure import StructureData
from aiida.orm.data.array.trajectory import TrajectoryData
from aiida.orm.utils import load_node, WorkflowFactory
import click
from ase import units
from ase import Atoms
from ase.io import write as ase_write
import aiida_utils
import sys
import numpy as np

#Define unit conversions
ANGSTROM_TO_BOHRRADIUS = 1./units.Bohr
EV_PER_ANGSTROM_TO_HARTREE_PER_BOHRRADIUS = units.Bohr/units.Hartree
EV_TO_HARTREE = 1/units.Hartree
conversion_aiida_ase_forces = 1.0000000908541917
conversion_aiida_ase_energy = 1.0000000945842116
print('exporting taking into account that aiida has its own conversion factors ....')

def get_allnodes_fromgroup(group_name):
    qb = QueryBuilder()
    qb.append(Group, filters={'name': group_name}, tag='g')
    qb.append(Node, tag='job', member_of='g')
    all_nodes = [x[0] for x in qb.all()]
    return all_nodes

def write_runner_commentline(fileout, uuid, extra_comments={}):
    fileout.write("begin\ncomment ")
    fileout.write("uuid: {} ".format(uuid))
    for label in extra_comments:
        fileout.write("{}: {} ".format(label, extra_comments[label]))
    fileout.write("\n")
    return

def write_runner_cell(fileout, cell):
    cell = cell * ANGSTROM_TO_BOHRRADIUS
    for idx_cell, i in enumerate(cell):
        fileout.write("lattice %.10f %.10f %.10f\n" %
                      (cell[idx_cell][0], cell[idx_cell][1], cell[idx_cell][2]))
    return

def write_runner_atomlines(fileout, atomiccoord_array, elements, atomicforce_array=None):
    """
    Assumes input units of eV and angstrom
    """
    if atomicforce_array is None:
        atomicforce_array = np.zeros(atomiccoord_array.shape)

    atomiccoord_array = atomiccoord_array * ANGSTROM_TO_BOHRRADIUS
    atomicforce_array = atomicforce_array * EV_PER_ANGSTROM_TO_HARTREE_PER_BOHRRADIUS

    for i in range(len(atomiccoord_array)):
        atomiccoord = atomiccoord_array[i]
        element = elements[i]
        atomicforce = atomicforce_array[i]

        fileout.write("atom   %.6f    %.6f   %.6f "
                      "%s  0.0   0.0  "
                      "%.10f  %.10f  %.10f\n" %
                      (atomiccoord[0], atomiccoord[1], atomiccoord[2],
                       element,
                       atomicforce[0], atomicforce[1], atomicforce[2]))
    return

def write_runner_finalline(fileout, energy=0, charge=0):
    fileout.write("energy %.15f\n" % (energy*EV_TO_HARTREE))
    fileout.write("charge %.15f\nend\n" % charge)
    return


def write_structure_torunner(fileout, structure_node, extra_comments={},stress=False):
    # get structure path, if applicable
    ase_structure = structure_node.get_ase()

    cell = ase_structure.get_cell()
    print('cell  :',cell)
    print('stress:',stress)
    positions = ase_structure.get_positions()
    elements = ase_structure.get_chemical_symbols()

    write_runner_commentline(fileout, structure_node.uuid, extra_comments=extra_comments)
    write_runner_cell(fileout, cell)
    write_runner_atomlines(fileout, positions, elements)
    write_runner_finalline(fileout)
    return

def get_timesorted_basenodes(relaxworknode):
    q = QueryBuilder()
    q.append(WorkCalculation, filters={"uuid": relaxworknode.uuid}, tag="relaxworknode")
    q.append(WorkCalculation, output_of="relaxworknode",
             project=["id", "ctime",  "*"],  tag="calc")
    q.order_by({"calc": "ctime"})
    timesorted_scf = [x[2] for x in q.all()]
    return timesorted_scf

def get_timesorted_scfs(worknode, relax_worknode=False):
    q = QueryBuilder()
    q.append(WorkCalculation, filters={"uuid": worknode.uuid}, tag="worknode")
    output_tag = "worknode"
    if relax_worknode:
        output_tag = "worknode2"
        q.append(WorkCalculation, tag=output_tag, output_of="worknode")
    q.append(Calculation, output_of=output_tag, project=["id", "ctime",  "*"],  tag="calc")
    q.order_by({"calc": "ctime"})
    timesorted_scf = [x[2] for x in q.all()]
    return timesorted_scf

def write_pwbase_torunner(fileout, pwbasenode, extra_comments={},stress=False):
    print()
    print('using write_pwbase_torunner',pwbasenode.uuid)
    scf_node = get_timesorted_scfs(pwbasenode)[-1]

    ase_structure = scf_node.inp.structure.get_ase()
    cell = ase_structure.get_cell()

    print(scf_node,'cell[0]',cell[0]*ANGSTROM_TO_BOHRRADIUS)
    print('cell  :',cell[0],cell[1],cell[2])
    vvv = ase_structure.get_volume()/ase_structure.get_number_of_atoms()
    print('volume:',ase_structure.get_volume()/ase_structure.get_number_of_atoms())
    print('volume:',ase_structure.get_volume())
    print('stress:',stress)
    print('c44 ST:',stress[2][1]*1000./2.,"GPa")
    e0 = -2149.85053966
    V0 = 16.476413798802568*4.
    strain = 0.2

    positions = ase_structure.get_positions()
    elements = ase_structure.get_chemical_symbols()

    try:
        atomicforce_array = scf_node.out.output_array.get_array('forces')[-1] * conversion_aiida_ase_forces
    except KeyError:
        print('Error forces not obtained, skipping this structure...')
        return
    energy = scf_node.out.output_parameters.get_attr('energy') * conversion_aiida_ase_energy
    vol = ase_structure.get_volume()
    #((-2149.85053966--2149.85051177)*2./((0.2/100.)**2.)/(16.47639732238896*4))/aseunits.GPa
    #((e0--2149.85051177)*2./((0.2/100.)**2.)/(16.47639732238896*4))/aseunits.GPa

    c44 = ((e0/V0-energy/vol)*2./((strain/100.)**2.))/units.GPa
    c44 = ((e0-energy)*2./((0.2/100.)**2.)/(16.47639732238896*4))/units.GPa
    c44 = ((e0-energy)*2./((0.2/100.)**2.)/(vol))/units.GPa
    c44 = ((e0-energy)*2./vol/((0.2/100.)**2.))/units.GPa
    c44 = ((e0/vol-energy/vol)*2./1./((0.2/100.)**2.))/units.GPa
    c44 = ((e0/V0-energy/vol)*2./1./((0.2/100.)**2.))/units.GPa

    print('V0',V0,'vol',vol)
    print('ene   :',energy,"eV")
    print('V0, ene (meV_pa):',vvv,energy/ase_structure.get_number_of_atoms()*1000.)
    print('c44 EN:',c44,"GPa")

    write_runner_commentline(fileout, pwbasenode.uuid, extra_comments=extra_comments)
    write_runner_cell(fileout, cell)
    write_runner_atomlines(fileout, positions, elements, atomicforce_array=atomicforce_array)
    write_runner_finalline(fileout, energy=energy)
    return

def get_timesorted_trajectories(relaxworkcalc):
    q = QueryBuilder()
    q.append(WorkCalculation, filters={"uuid": relaxworkcalc.uuid}, tag="relaxworkcalc")
    q.append(WorkCalculation, tag="baseworkcalc", output_of="relaxworkcalc")
    q.append(Calculation, output_of="baseworkcalc", tag="calc")
    q.append(TrajectoryData, output_of="calc", project=["id", "ctime",  "*"], tag="traj")
    q.order_by({"traj": "ctime"})
    timesorted_trajectories = [x[2] for x in q.all()]
    return timesorted_trajectories

def get_arraysbyname_fromtrajectories(timesorted_trajectories, arrayname):
    timesorted_arrays = [x.get_array(arrayname) for x in timesorted_trajectories]
    return np.concatenate(timesorted_arrays)

def write_pwrelax_torunner(fileout, relax_node, write_only_relaxed, verbose, extra_comments={}):
    print('using write_pwrelax_torunner')
    trajectories = get_timesorted_trajectories(relax_node)

    timesorted_cells = get_arraysbyname_fromtrajectories(trajectories, 'cells')
    timesorted_positions = get_arraysbyname_fromtrajectories(trajectories, 'positions')
    timesorted_forces = get_arraysbyname_fromtrajectories(trajectories, 'forces') * conversion_aiida_ase_forces
    if len(timesorted_forces) == 1:
        energy = relax_node.out.CALL.out.CALL.out.output_parameters.get_attr('energy') * conversion_aiida_ase_energy
        timesorted_energy = [energy]
    else:
        timesorted_energy = get_arraysbyname_fromtrajectories(trajectories, 'energy') * conversion_aiida_ase_energy
    elements = trajectories[0].get_array('symbols') # assume unchangin

    extra_comments={"trajectory_step":None}
    print('relaxation steps:',len(timesorted_cells))

    if write_only_relaxed == False:
        for i in range(len(timesorted_cells)):
            extra_comments["trajectory_step"] = i
            if verbose:
                a = timesorted_cells[i][0]
                b = timesorted_cells[i][1]
                c = timesorted_cells[i][2]
                vol=np.dot(a,np.cross(b,c))
                maxforce = np.abs(timesorted_forces[i]).max()
                print('extra',str(extra_comments).ljust(25," "),
                      "volume",vol,
                      "energy",str(timesorted_energy[i]).ljust(20),
                      "maxforce",maxforce)
            ###################################################################
            # at some point we'll need to update the runner write_runner part
            # since it can not be read by n2p2 in the current form (also it looks ok)
            # the ase_write("outt.runner",frame,format='runner',append=True)
            # works but not in my currently installed aiida version
            ###################################################################
            #frame = Atoms(elements, positions=timesorted_positions[i])
            #frame.set_cell(timesorted_cells[i])
            #ka = extra_comments["trajectory_step"]
            #ase_write("outt.runner",frame,format='runner',append=True) #,comment=str(relax_node.uuid+" trajectory_step "+str(ka)+"\n")) #+" "+extra_comments)
            write_runner_commentline(fileout, relax_node.uuid, extra_comments=extra_comments)
            write_runner_cell(fileout, timesorted_cells[i])
            write_runner_atomlines(fileout,
               timesorted_positions[i], elements, atomicforce_array=timesorted_forces[i])
            write_runner_finalline(fileout, energy=timesorted_energy[i])

    if bool(relax_node.inp.final_scf):
        print('final')
        final_basenode = get_timesorted_basenodes(relax_node)[-1]
        extra_comments["trajectory_step"] = "final_scf"
        extra_comments["parent_uuid"] = relax_node.uuid
        write_pwbase_torunner(fileout, final_basenode, extra_comments=extra_comments)

    return


# show default values in click
orig_init = click.core.Option.__init__
def new_init(self, *args, **kwargs):
    orig_init(self, *args, **kwargs)
    self.show_default = True
click.core.Option.__init__ = new_init
# get help also with -h
CONTEXT_SETTINGS = dict(help_option_names=['-h', '--help'])


@click.command(context_settings=CONTEXT_SETTINGS)
@click.option('-gn', '--group_name', default=None,
              type=str, help="Group to export identified by name")
@click.option('-f', '--filename', required=False, default=False,
         type=str, help="filename for outputfile, default = work_group+\".input.data\"")
@click.option('-wor', '--write_only_relaxed', required=False, default=False,
         is_flag=True, help="only write the final relaxed structure")
@click.option('-sreadme', '--supress_readme', is_flag=True,
         help="supresses the generation of a readme file")
@click.option('-v', '--verbose', is_flag=True,
         type=str, help="Enables verbosity")


def createjob(group_name, filename, write_only_relaxed, supress_readme, verbose):
    ''' e.g.
    ./aiida_export_group_to_runner.py -gn Al6xxxDB_structuregroup
    '''
    all_nodes = get_allnodes_fromgroup(group_name)
    if not supress_readme:
        aiida_utils.create_READMEtxt()

    add_to_filename = "__all_steps"
    if write_only_relaxed == True:
        add_to_filename = "__only_relaxed"

    if filename == False:
        file = "aiida_exported_group_"+group_name+add_to_filename+".input.data"
        fileout = open(file, "w")
    else:
        file = filename
        fileout = open(filename, "w")

    if verbose:
        print('file           :', file)
        print('structure_group:', group_name)

    def get_stress(work):
        return work.out.output_parameters.get_dict()['stress']

    for node in all_nodes:
        stress = get_stress(load_node(node.uuid))

        if verbose:
            print("Writing node: {}".format(node.uuid))

        if isinstance(node, StructureData):
            print('using write_structure_torunner')
            write_structure_torunner(fileout, node,stress=stress)
        elif isinstance(node, WorkCalculation):
            process_label = node.get_attrs()['_process_label']
            if process_label == "PwBaseWorkChain":
                write_pwbase_torunner(fileout, node,stress=stress)
            elif process_label == "PwRelaxWorkChain":
                write_pwrelax_torunner(fileout, node, write_only_relaxed,verbose)
            else:
                print("Could not identify node, skipping")
        else:
            print("Could not identify node, skipping")

    fileout.close()
    return


if __name__ == "__main__":
    createjob()
