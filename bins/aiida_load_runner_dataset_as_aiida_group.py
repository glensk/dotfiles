#!/usr/bin/env python
import aiida
aiida.load_dbenv()
from aiida.orm.group import Group
from aiida.orm.utils import load_node
from aiida.orm.data.structure import StructureData
import ase
import ase.io
import click
import aiida_utils

def add_parentstructure_extras(structurenode, parent_uuid):
    # NOTE: consider adding a check if parent_extras is already assigned
    structure_extras = structurenode.get_extras()
    parent_extras = load_node(parent_uuid).get_extras()
    for key, value in parent_extras.items():
        if key not in structure_extras:
            structurenode.set_extra(key, value)
    structurenode.set_extra('parent_extras', True)
    return


@click.command()
@click.option('-d', '--dataset_path', required=True)
@click.option('-gn', '--group_name', required=True)
@click.option('-gd', '--group_description', default="")
@click.option('-pcp', '--parse_comments_path', is_flag=True)
@click.option('-pcs', '--parse_comments_structure', is_flag=True)
@click.option('-sreadme', '--supress_readme', is_flag=True,
        help="supresses the generation of a readme file")
def launch(dataset_path, group_name, group_description,
           parse_comments_path, parse_comments_structure,supress_readme):
    print "loading dataset: {} to group: {}".format(dataset_path, group_name)

    # Setup/Retrieve the Group
    g = Group.get_or_create(name=group_name, description=group_description)[0]
    # Loop over structures in the dataset_path, storing the nodes then adding them to the group
    i = 0
    while True:
        try:
            ase_structure = ase.io.read(dataset_path, index=i, format="runner")
        except StopIteration:
            break
        # setup and store the ase structure as an aiida StructureData node
        aiida_structure = StructureData()
        aiida_structure.set_ase(ase_structure)
        aiida_structure_stored = aiida_structure.store()

        # add in the dataset_path line if possible
        if parse_comments_path:
            try:
                structure_path = ase_structure.comment.strip().split()[-1][3:]
                aiida_structure_stored.set_extra("structure_path", structure_path)
            except AttributeError:
                print "could not set structure_path on {}".format(ase_structure)
                pass
        # Add in details of parent structure uuid
        if parse_comments_structure:
            try:
                parent_uuid = ase_structure.comment.strip().split()[-1]
                aiida_structure_stored.set_extra("parent_uuid", parent_uuid)
                add_parentstructure_extras(aiida_structure_stored, parent_uuid)
            except AttributeError:
                print "could not set parent_uuid on {}".format(ase_structure)
                pass


        # add in the chemical formula and number of atoms if possible
        try:
            aiida_structure_stored.set_extra("num_atoms",
                                             len(ase_structure))
            aiida_structure_stored.set_extra("chem_formula",
                          ase_structure.get_chemical_formula())
        except AttributeError:
            print "could not set either num_atoms or chemical_formula " \
                  " on {}".format(ase_structure)
            pass

        # add the structure to the group
        g.add_nodes(aiida_structure_stored)
        g.store()
        i += 1

    if not supress_readme:
        aiida_utils.create_READMEtxt()

if __name__ == "__main__":
    try:
        ase.io.read("", format="runner")
    except ValueError:
        raise ValueError("You need a version of ase that can read runner files")
    except IOError:
        pass
    launch()
