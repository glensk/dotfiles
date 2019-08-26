#!/usr/bin/env python
import click
import time
from aiida.cmdline.params import options, arguments
from aiida.cmdline.params.types import DataParamType
from aiida.cmdline.utils import decorators
import sys
import aiida_utils

import aiida
aiida.load_dbenv()
from aiida.work.workfunctions import workfunction


def retrieve_alluncalculated_structures(structure_group_name,
                                        workchain_group_name=None):
    from aiida.orm.group import Group
    from aiida.orm.data.structure import StructureData
    from aiida.orm.calculation import WorkCalculation
    from aiida.orm.querybuilder import QueryBuilder

    sqb = QueryBuilder()
    sqb.append(Group, filters={'name': structure_group_name}, tag='g')
    sqb.append(StructureData, project='id', tag='s', member_of='g')
    sqb.append(WorkCalculation, tag='job', output_of='s')

    filters = {}
    if workchain_group_name:
        filters = {'name': workchain_group_name}
    sqb.append(Group, group_of='job', filters=filters)

    ids_dealt_with = [_ for _, in sqb.distinct().all()] or [-1]  # prevent empty list

    # # Now the main query:
    qb = QueryBuilder()
    qb.append(Group, filters={'name': structure_group_name}, tag='g')
    qb.append(StructureData, project='*', tag='s', member_of='g',
              filters={'id': {'!in': ids_dealt_with}})  # filter out calculated '!in' for not in

    res = [x[0] for x in qb.all()]
    return res

def retrieve_numactive_calculations():
    from aiida.orm.calculation import JobCalculation
    from aiida.orm.querybuilder import QueryBuilder
    qb = QueryBuilder()
    qb.append(JobCalculation,
              filters={'attributes.process_state':
                       {'!in': ['finished', 'excepted', 'killed']}}
    )
    return len(qb.all())


def get_numelectrons_structure_upffamily(structure, pseudos):

    def parse_numelectrons_upfpath(upfpath):
        for line in open(upfpath):
            if "valence" in line.lower() and "z" in line.lower():
                if len(line.split("=")) == 2:
                    num_e = int(float((line.split("=")[-1].strip().strip('"'))))
                elif len(line.split()) == 3:
                    num_e = int(float(line.split()[0]))
                else:
                    raise Exception("Could not parse {}".format(upfpath))
        return num_e

    def build_upf_numelectrons_dict(structure_ase, pseudos):
        element_nume_dict = {}
        for element in set(structure_ase.get_chemical_symbols()):
            upfpath = pseudos[element].get_file_abs_path()
            element_nume_dict[element] = parse_numelectrons_upfpath(upfpath)
        return element_nume_dict

    structure_ase = structure.get_ase()

    element_nume_dict = build_upf_numelectrons_dict(structure_ase, pseudos)

    num_e = 0
    for element in structure_ase.get_chemical_symbols():
        num_e += element_nume_dict[element]

    return num_e


def get_kmeshfrom_kptper_recipang(aiida_structure, kptper_recipang):
    import numpy as np

    ase_structure = aiida_structure.get_ase()
    reci_cell = ase_structure.get_reciprocal_cell()
    kmesh = [np.ceil(kptper_recipang * np.linalg.norm(reci_cell[i]))
             for i in range(len(reci_cell))]
    return kmesh


def get_nummachines(structure, pseudo_familyname):
    # NOTE: used very adhoc guess for nodes, assuming quadratic scaling
    from aiida.orm.data.upf import get_pseudos_from_structure
    pseudos = get_pseudos_from_structure(structure, pseudo_familyname)
    num_electrons = get_numelectrons_structure_upffamily(structure, pseudos)
    a2 = 1.5*10**-6
    a1 = 5.7*10**-3
    a0 = 2
    numnodes = a2*num_electrons**2+a1*num_electrons+a0
    numnodes = max(round(numnodes/2)*2, 2)  # force even # of nodes
    return numnodes


def get_nk(num_machines, code):
    def nk_nump_evenlydivisible(nk, nump):
        nk = float(nk)
        nump = float(nump)

        if round(nump/nk) == nump/nk:
            return True
        else:
            return False

    nk = str(max(4, int(num_machines/2)))  # adhoc guess

    # check if our choice is valid
    computer = code.get_computer()
    ppm = computer.get_default_mpiprocs_per_machine()
    nump = num_machines * ppm
    if not nk_nump_evenlydivisible(nk, nump):
        raise Exception("Error number processors: {} "
                        "is not divisible by nk: {}".format(nump, nk))

    return nk


@workfunction
def wf_getkpoints(aiida_structure, kptper_recipang):
    from aiida.orm.data.array.kpoints import KpointsData

    def get_kmeshfrom_kptper_recipang(aiida_structure, kptper_recipang):
        import numpy as np
        kptper_recipang = kptper_recipang.value

        ase_structure = aiida_structure.get_ase()
        reci_cell = ase_structure.get_reciprocal_cell()
        kmesh = [np.ceil(kptper_recipang * np.linalg.norm(reci_cell[i]))
                 for i in range(len(reci_cell))]
        return kmesh

    kpoints_mesh = get_kmeshfrom_kptper_recipang(aiida_structure, kptper_recipang)
    kpoints = KpointsData(kpoints_mesh=kpoints_mesh)

    return kpoints


@workfunction
def wf_setupparams(base_parameter, structure,
                   pseudo_familyname, nume2bnd_ratio):
        from aiida.orm.data.upf import get_pseudos_from_structure
        from aiida.orm.data.parameter import ParameterData

        pseudos = get_pseudos_from_structure(structure, pseudo_familyname.value)
        nelec = get_numelectrons_structure_upffamily(structure, pseudos)
        nbnd = nelec * nume2bnd_ratio.value
        nbnd = max(nbnd, 20) # minimum of 20 bands to avoid certain crashes
        parameter_dict = base_parameter.get_dict()
        parameter_dict['SYSTEM']['nbnd'] = nbnd
        parameters = ParameterData(dict=parameter_dict)

        return parameters


@click.command()
@click.option('-c', '--code_node', required=False, default="1",
              help="node of code to use")
@click.option('-sg', '--structure_group_name', required=True,
              help='input group of structures to submit workchains on')
@click.option('-wg', '--workchain_group_name', required=True,
              help='output group of workchains')
@click.option('-bp', '--base_parameter_node', required=True,
              help='node of base ParameterData to setup calculations')
@click.option('-pfn', '--pseudo_familyname', required=True,
              help='name of pseudopotential family to use')
@click.option('-kra', '--kptper_recipang', required=True,
              help='number of kpoints to use per reciprocal angstrom')
@click.option('-ber', '--nume2bnd_ratio', required=True,
              help='band to electron ratio')
@click.option('-cm', '--calc_method', default='scf',
              type=click.Choice(["scf", "relax", "vc-relax"]),
              help='The type of calculation to perform')
@click.option('-mws', '--max_wallclock_seconds', default=8*60*60,
              help='maximum wallclock time per job in seconds')
@click.option('-mac', '--max_active_calculations', default=300,
              help='maximum number of active calculations')
@click.option('-nnd', '--number_of_nodes', default=None,
              help='Force all calculations to use the specified number of nodes')
@click.option('-memgb', '--memory_gb', default=None,
              help='specify the amount of memory for all jobs in GB')
@click.option('-nd', '--ndiag', default=None,
              help='ndiag setting to be passed direct to QE')
@click.option('-nk', '--npools', default=None,
              help='npools setting to be passed direct to QE')
@click.option('-sli', '--sleep_interval', default=10*60,
              help='time to wait (sleep) between calculation submissions')
@click.option('-zmo', '--z_movement_only', is_flag=True,
              help='Restricts movement to the z direction only. For relaxing stacking fault')
@click.option('-kwd', '--keep_workdir', is_flag=True,
              help='Keep the workdir files after running')
@click.option('-dr', '--dryrun', is_flag=True,
              help="Prints inputs but does not launch anything")
@click.option('-rdb', '--run_debug', is_flag=True,
              help='run the script in debug mode. Submits one structure only'
                   ' and does not attach the output to the workchain_group')
def launch(code_node, structure_group_name, workchain_group_name,
           base_parameter_node, pseudo_familyname, kptper_recipang,
           nume2bnd_ratio, calc_method, max_wallclock_seconds, max_active_calculations,
           number_of_nodes, memory_gb, ndiag, npools,
           sleep_interval, z_movement_only, keep_workdir,dryrun, run_debug):
    from aiida.orm.group import Group
    from aiida.orm.utils import load_node, WorkflowFactory
    from aiida.orm.data.base import Bool, Float, Int, Str
    from aiida.orm.data.parameter import ParameterData
    from aiida.work.launch import submit

    valid_calc_methods = ['scf', 'relax', 'vc-relax']
    if calc_method not in valid_calc_methods:
        raise Exception("Invalid calc_method: {}".format(calc_method))

    # setup parameters
    code = load_node(code_node)
    structure_group = Group.get_from_string(structure_group_name)
    workchain_group = Group.get_or_create(name=workchain_group_name)[0]
    base_parameter = load_node(base_parameter_node)
    # announce if running in debug mode
    if run_debug:
        print "Running in debug mode!"

    # Load all the structures in the structure group, not-yet run in workchain_group_name
    uncalculated_structures = retrieve_alluncalculated_structures(
                                structure_group_name,
                                workchain_group_name=workchain_group_name
    )
    if len(uncalculated_structures) == 0:
        print("All structures in {} already have associated workchains in "
              "the group {}".format(structure_group_name, workchain_group_name))
        sys.exit()

    # determine number of calculations to submit
    running_calculations = retrieve_numactive_calculations()
    calcs_to_submit = max_active_calculations - running_calculations


    # submit calculations
    for structure in uncalculated_structures:
        print "Preparing to launch {}".format(structure)

        # ensure no more than the max number of calcs are submitted
        while (calcs_to_submit <= 0):
            running_calculations = retrieve_numactive_calculations()
            calcs_to_submit = max_active_calculations - running_calculations
            if calcs_to_submit <= 0:  # in case jobs finished during submission
                print("{} calcs running,"
                      "max num calcs {} waiting....".format(
                          running_calculations, max_active_calculations))
                time.sleep(sleep_interval)

        # start timer to inspect job submission times
        from timeit import default_timer as timer
        start = timer()

        # determine number of bands & setup the parameters
        parameters = wf_setupparams(base_parameter,
                                    structure,
                                    Str(pseudo_familyname),
                                    Float(nume2bnd_ratio))

        # determine kpoint mesh & setup kpoints
        kpoints = wf_getkpoints(structure, Int(kptper_recipang))

        # determine parallelization & resources (setup the settings & options)
        if number_of_nodes:
            num_machines = int(number_of_nodes)
        else:
            num_machines = get_nummachines(structure, pseudo_familyname)
            max_nodes_to_submit = 20
            if num_machines > max_nodes_to_submit:
                print "{} nodes requested, maximum is {}".format(num_machines, max_nodes_to_submit)
                print "If you wish to launch please choose nodes manually with --number_of_nodes"
                continue
        options_dict = {
            'max_wallclock_seconds': max_wallclock_seconds,
            'resources': {'num_machines': num_machines},
        }
        if memory_gb:
            options_dict['max_memory_kb'] = int(int(memory_gb)*1024*1024)
        if run_debug:
            num_machines = 2
            options_dict['resources']['num_machines'] = num_machines
            options_dict['max_wallclock_seconds'] = int(30*60)
            options_dict['queue_name'] = 'debug'
        workchain_options = ParameterData(dict=options_dict)

        if npools:
            nk = npools
        else:
            nk = get_nk(num_machines, code)
        settings_dict = {
            'cmdline': ['-nk', nk],
            'no_bands': True
            }
        if ndiag:
            settings_dict['cmdline'] += ['-ndiag', ndiag]
        if z_movement_only:
            num_atoms = len(structure.get_ase())
            coordinate_fix = [[True,True,False]]*num_atoms
            settings_dict['fixed_coords'] = coordinate_fix
        settings = ParameterData(dict=settings_dict)

        # setup inputs & submit workchain
        clean_workdir = not keep_workdir
        inputs = {
                  'structure': structure,
                  'settings': settings,
                  'clean_workdir': Bool(clean_workdir)
                  }
        base_inputs = {
                'code': code,
                'pseudo_family': Str(pseudo_familyname),
                'kpoints': kpoints,
                'parameters': parameters,
                'options': workchain_options,
                'settings': settings,
                }
        if calc_method == 'scf':
            PwBaseWorkChain = WorkflowFactory('quantumespresso.pw.base')
            inputs.update(base_inputs)
        elif calc_method == 'relax':
            PwBaseWorkChain = WorkflowFactory('quantumespresso.pw.relax')
            inputs['base'] = base_inputs
            inputs['relaxation_scheme'] = Str('relax')
            inputs['final_scf'] = Bool(False)
            inputs['meta_convergence'] = Bool(False)
        elif calc_method == 'vc-relax':
            PwBaseWorkChain = WorkflowFactory('quantumespresso.pw.relax')
            inputs['base'] = base_inputs
            inputs['relaxation_scheme'] = Str('vc-relax')
            inputs['final_scf'] = Bool(True)
            inputs['meta_convergence'] = Bool(True)
        else:
            raise Exception("Invalid calc_method: {}".format(calc_method))

        def print_timing(start):
            end = timer()
            time_elapsed = end - start
            print "timing: {}s".format(time_elapsed)

        if dryrun:
            print "ase_structure: {}".format(structure.get_ase())
            print "aiida_settings: {}".format(settings.get_dict())
            print "aiida_options: {}".format(workchain_options.get_dict())
            print "aiida_inputs: {}".format(inputs)
            print_timing(start)
        else:
            node = submit(PwBaseWorkChain, **inputs)
            workchain_group.add_nodes([node])
            print "WorkChain: {} submitted".format(node)
            print_timing(start)
        calcs_to_submit -= 1

        if run_debug:
            sys.exit()

    aiida_utils.create_READMEtxt()



if __name__ == "__main__":
    launch()

