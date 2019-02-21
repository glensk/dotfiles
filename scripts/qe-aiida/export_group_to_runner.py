#!/usr/bin/env python
import aiida
aiida.load_dbenv()
from aiida.orm.querybuilder import QueryBuilder
from aiida.orm import Group, WorkCalculation
import click

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
def createjob(work_group):
    ''' e.g.
    ./aiida_export_group_to_runner.py -wg kmc_1000K_4
    ./aiida_export_group_to_runner.py -wg Al6xxxDB_passingsubset works
    '''
    print('work_group:', work_group)
    qb = QueryBuilder()
    qb.append(Group, filters={'name': work_group}, tag='g')
    qb.append(WorkCalculation, tag='job', member_of='g')
    all_works = [x[0] for x in qb.all()]  # all workchains

    def get_workcalc_runnerdata(worknode):
        #TODO enable multi-structure support
        try:
            ase_structure = worknode.out.output_structure.get_ase()
        except Exception:
            ase_structure = worknode.inp.structure.get_ase()

        energy = worknode.out.output_parameters.get_attrs()['energy']  # units?

        #TODO: this section splits for SCF and relax, should fix & merge
        print "worknode: ", worknode
        try:
            # SCF
            forces = worknode.out.output_array.get_array('forces')  # units?
        except Exception:
            # Relax (probably)
            forces = worknode.out.CALL.out.CALL.out.output_trajectory.get_array('forces')

        #TODO: this section splits for SCF and relax, should fix & merge
        try:
            # SCF
            path = worknode.out.CALL.out.retrieved.get_abs_path()
        except Exception:
            # Relax (probably)
            path = "path support only for SCF calcs"
        return ase_structure, energy, forces, worknode.uuid, path

    angstrom_to_bohrradius = 1.8897261
    eV_to_Hartree = 0.036749325
    eV_per_angstrom_to_hartree_per_bohrradius = 0.019446905

    fileOut = open(work_group+".input.data", "w")

    for idx, workchain in enumerate(all_works):
        ase_structure, energy, forces, uuid, path = get_workcalc_runnerdata(workchain)
        print(idx, "ene (eV)", energy, uuid, path)
        forces = forces[-1]

        work_uuid = workchain.uuid
        fileOut.write("begin\ncomment uuid: {}\n".format(work_uuid))

        # write the cell
        cell = ase_structure.cell*angstrom_to_bohrradius
        for idx_cell, i in enumerate(cell):
            fileOut.write("lattice %.10f %.10f %.10f\n" %
                          (cell[idx_cell][0], cell[idx_cell][1], cell[idx_cell][2]))

        #  write the positions
        nr_of_atoms = ase_structure.positions.shape[0]
        print nr_of_atoms, len(forces)
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
