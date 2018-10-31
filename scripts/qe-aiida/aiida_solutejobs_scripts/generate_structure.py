import ase
from ase.build import sort
import pandas as pd
import numpy as np
import sys


def supercellsolute_ase(a0, supercell_shape, matrix,
                                   solute):
    a1 = np.array([a0,0.,0.])
    a2 = np.array([0.,a0,0.])
    a3 = np.array([0.,0.,a0])
    t1 = np.array([0.0,0.0,0.0])
    t2 = np.array([0.0,0.5,0.5])
    t3 = np.array([0.5,0.0,0.5])
    t4 = np.array([0.5,0.5,0.0])

    atoms = ase.Atoms('{}4'.format(matrix),
                      cell=[a1,a2,a3], pbc=True,
                      scaled_positions=[t1,t2,t3,t4])
    supercell_atoms = atoms*np.array(supercell_shape, dtype=int)
    supercell_atoms[0].symbol=solute
    return supercell_atoms

def return_nn_distanceAndIndex(supercell):
    supercell_frame = pd.DataFrame(supercell.get_positions(),
                                   columns=['x','y','z'])
    supercell_frame['distance'] = supercell_frame.apply(
                                    np.linalg.norm, axis=1)
    supercell_frame['distance'] = supercell_frame['distance'].round(9)

    supercell_shape = supercell.get_cell()
    lx = supercell_shape[0][0]
    ly = supercell_shape[1][1]
    lz = supercell_shape[2][2]

    supercell_frame = supercell_frame[
                 (supercell_frame['x'] <= lx/2.) &
                 (supercell_frame['y'] <= ly/2.) &
                 (supercell_frame['z'] <= lz/2.) ]

    supercell_frame.sort_values('distance', inplace=True)
    supercell_frame.drop_duplicates(subset='distance',
                                    inplace=True)



    distance_values = supercell_frame['distance'].values
    index_values = supercell_frame.index.values

    return distance_values.tolist(), index_values.tolist()

def main():
    if len(sys.argv) != 7:
        sys.exit("generate_structure.py a0 supercell_shape matrix solute position outputfile")

    a0, supercell_shape, matrix, solute, position, outputfile = sys.argv[1:]

    a0 = float(a0)
    supercell_shape = supercell_shape.split(',')
    if len(supercell_shape) != 3:
        sys.exit("supercell_shape must be of the form Nx,Ny,Nz")

    if position == 'pure':
        structure = supercellsolute_ase(a0, supercell_shape, matrix, matrix)
        structure_label='pure'
    elif position == 'singlesolute':
        structure = supercellsolute_ase(a0, supercell_shape, matrix, solute)
        structure_label='singlesolute'
    else:
        solute_structure = supercellsolute_ase(a0, supercell_shape, matrix, solute)
        solsol_distances, solsol_indexes = return_nn_distanceAndIndex(solute_structure)

        if position == "print_all":
            print "There are {} valid neighbours\nWith distances of {}".format(
                    len(solsol_indexes[1:]), solsol_distances[1:])
            sys.exit()

        position = int(position)
        if position > len(solsol_indexes):
            sys.exit('This structure only has sensible nearest neighbours up to N={}'.format(
                       len(solsol_indexes)))


        structure = solute_structure
        structure[solsol_indexes[position]].symbol = solute
        structure_label = solsol_distances[position]

    print structure_label
    structure = sort(structure)
    structure.write(outputfile, format='espresso-in') 


if __name__ == "__main__":
   main()



#return atomic_distance (or named tag)  & print POSCAR


