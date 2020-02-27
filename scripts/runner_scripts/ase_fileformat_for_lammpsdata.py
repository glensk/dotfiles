import re,copy,sys
import numpy as np
from numpy.linalg import norm


from ase.atoms import Atoms
from ase.parallel import paropen
from ase.calculators.lammpslib import unit_convert
from ase.utils import basestring


def read_lammps_data(fileobj, Z_of_type=None, style='full', sort_by_id=False,
                     units="metal"):
    """Method which reads a LAMMPS data file.

    sort_by_id: Order the particles according to their id. Might be faster to
    switch it off.
    Units are set by default to the style=metal setting in LAMMPS.
    """
    if isinstance(fileobj, basestring):
        f = paropen(fileobj)
    else:
        f = fileobj

    # load everything into memory
    lines = f.readlines()

    # begin read_lammps_data
    comment = None
    N = None
    # N_types = None
    xlo = None
    xhi = None
    ylo = None
    yhi = None
    zlo = None
    zhi = None
    xy = None
    xz = None
    yz = None
    pos_in = {}
    travel_in = {}
    mol_id_in = {}
    mmcharge_in = {}
    mass_in = {}
    vel_in = {}
    bonds_in = []
    angles_in = []
    dihedrals_in = []

    sections = ["Atoms",
                "Velocities",
                "Masses",
                "Charges",
                "Ellipsoids",
                "Lines",
                "Triangles",
                "Bodies",
                "Bonds",
                "Angles",
                "Dihedrals",
                "Impropers",
                "Impropers Pair Coeffs",
                "PairIJ Coeffs",
                "Pair Coeffs",
                "Bond Coeffs",
                "Angle Coeffs",
                "Dihedral Coeffs",
                "Improper Coeffs",
                "BondBond Coeffs",
                "BondAngle Coeffs",
                "MiddleBondTorsion Coeffs",
                "EndBondTorsion Coeffs",
                "AngleTorsion Coeffs",
                "AngleAngleTorsion Coeffs",
                "BondBond13 Coeffs",
                "AngleAngle Coeffs"]
    header_fields = ["atoms",
                     "bonds",
                     "angles",
                     "dihedrals",
                     "impropers",
                     "atom types",
                     "bond types",
                     "angle types",
                     "dihedral types",
                     "improper types",
                     "extra bond per atom",
                     "extra angle per atom",
                     "extra dihedral per atom",
                     "extra improper per atom",
                     "extra special per atom",
                     "ellipsoids",
                     "lines",
                     "triangles",
                     "bodies",
                     "xlo xhi",
                     "ylo yhi",
                     "zlo zhi",
                     "xy xz yz"]
    sections_re = '(' + '|'.join(sections).replace(' ', '\\s+') + ')'
    header_fields_re = '(' + '|'.join(header_fields).replace(' ', '\\s+') + ')'

    section = None
    header = True
    for line in lines:
        if comment is None:
            comment = line.rstrip()
        else:
            line = re.sub("#.*", "", line).rstrip().lstrip()
            if re.match("^\\s*$", line):  # skip blank lines
                continue

        # check for known section names
        m = re.match(sections_re, line)
        if m is not None:
            section = m.group(0).rstrip().lstrip()
            header = False
            continue

        if header:
            field = None
            val = None
            # m = re.match(header_fields_re+"\s+=\s*(.*)", line)
            # if m is not None: # got a header line
            #   field=m.group(1).lstrip().rstrip()
            #   val=m.group(2).lstrip().rstrip()
            # else: # try other format
            #   m = re.match("(.*)\s+"+header_fields_re, line)
            #   if m is not None:
            #       field = m.group(2).lstrip().rstrip()
            #       val = m.group(1).lstrip().rstrip()
            m = re.match("(.*)\\s+" + header_fields_re, line)
            if m is not None:
                field = m.group(2).lstrip().rstrip()
                val = m.group(1).lstrip().rstrip()
            if field is not None and val is not None:
                if field == "atoms":
                    N = int(val)
                # elif field == "atom types":
                #     N_types = int(val)
                elif field == "xlo xhi":
                    (xlo, xhi) = [float(x) for x in val.split()]
                elif field == "ylo yhi":
                    (ylo, yhi) = [float(x) for x in val.split()]
                elif field == "zlo zhi":
                    (zlo, zhi) = [float(x) for x in val.split()]
                elif field == "xy xz yz":
                    (xy, xz, yz) = [float(x) for x in val.split()]

        if section is not None:
            fields = line.split()
            if section == "Atoms":  # id *
                id = int(fields[0])
                if style == 'full' and (len(fields) == 7 or len(fields) == 10):
                    # id mol-id type q x y z [tx ty tz]
                    pos_in[id] = (int(fields[2]), float(fields[4]),
                                  float(fields[5]), float(fields[6]))
                    mol_id_in[id] = int(fields[1])
                    mmcharge_in[id] = float(fields[3])
                    if len(fields) == 10:
                        travel_in[id] = (int(fields[7]),
                                         int(fields[8]),
                                         int(fields[9]))
                elif (style == 'atomic' and
                      (len(fields) == 5 or len(fields) == 8)):
                    # id type x y z [tx ty tz]
                    pos_in[id] = (int(fields[1]), float(fields[2]),
                                  float(fields[3]), float(fields[4]))
                    if len(fields) == 8:
                        travel_in[id] = (int(fields[5]),
                                         int(fields[6]),
                                         int(fields[7]))
                elif ((style == 'angle' or style == 'bond' or
                       style == 'molecular') and
                      (len(fields) == 6 or len(fields) == 9)):
                    # id mol-id type x y z [tx ty tz]
                    pos_in[id] = (int(fields[2]), float(fields[3]),
                                  float(fields[4]), float(fields[5]))
                    mol_id_in[id] = int(fields[1])
                    if len(fields) == 9:
                        travel_in[id] = (int(fields[6]),
                                         int(fields[7]),
                                         int(fields[8]))
                else:
                    raise RuntimeError("Style '{}' not supported or invalid "
                                       "number of fields {}"
                                       "".format(style, len(fields)))
            elif section == "Velocities":  # id vx vy vz
                vel_in[int(fields[0])] = (float(fields[1]),
                                          float(fields[2]),
                                          float(fields[3]))
            elif section == "Masses":
                mass_in[int(fields[0])] = float(fields[1])
            elif section == "Bonds":  # id type atom1 atom2
                bonds_in.append((int(fields[1]),
                                 int(fields[2]),
                                 int(fields[3])))
            elif section == "Angles":  # id type atom1 atom2 atom3
                angles_in.append((int(fields[1]),
                                  int(fields[2]),
                                  int(fields[3]),
                                  int(fields[4])))
            elif section == "Dihedrals": # id type atom1 atom2 atom3 atom4
                dihedrals_in.append((int(fields[1]),
                                     int(fields[2]),
                                     int(fields[3]),
                                     int(fields[4]),
                                     int(fields[5])))

    # set cell
    cell = np.zeros((3, 3))
    cell[0, 0] = xhi - xlo
    cell[1, 1] = yhi - ylo
    cell[2, 2] = zhi - zlo
    if xy is not None:
        cell[1, 0] = xy
    if xz is not None:
        cell[2, 0] = xz
    if yz is not None:
        cell[2, 1] = yz

    # initialize arrays for per-atom quantities
    positions = np.zeros((N, 3))
    numbers = np.zeros((N), int)
    ids = np.zeros((N), int)
    types = np.zeros((N), int)
    if len(vel_in) > 0:
        velocities = np.zeros((N, 3))
    else:
        velocities = None
    if len(mass_in) > 0:
        masses = np.zeros((N))
    else:
        masses = None
    if len(mol_id_in) > 0:
        mol_id = np.zeros((N), int)
    else:
        mol_id = None
    if len(mmcharge_in) > 0:
        mmcharge = np.zeros((N), float)
    else:
        mmcharge = None
    if len(travel_in) > 0:
        travel = np.zeros((N, 3), int)
    else:
        travel = None
    if len(bonds_in) > 0:
        bonds = [""] * N
    else:
        bonds = None
    if len(angles_in) > 0:
        angles = [""] * N
    else:
        angles = None
    if len(dihedrals_in) > 0:
        dihedrals = [""] * N
    else:
        dihedrals = None

    ind_of_id = {}
    # copy per-atom quantities from read-in values
    for (i, id) in enumerate(pos_in.keys()):
        # by id
        ind_of_id[id] = i
        if sort_by_id:
            ind = id-1
        else:
            ind = i
        type = pos_in[id][0]
        positions[ind, :] = [pos_in[id][1], pos_in[id][2], pos_in[id][3]]
        if velocities is not None:
            velocities[ind, :] = [vel_in[id][0], vel_in[id][1], vel_in[id][2]]
        if travel is not None:
            travel[ind] = travel_in[id]
        if mol_id is not None:
            mol_id[i] = mol_id_in[id]
        if mmcharge is not None:
            mmcharge[i] = mmcharge_in[id]
        ids[i] = id
        # by type
        types[ind] = type
        if Z_of_type is None:
            numbers[ind] = type
        else:
            numbers[ind] = Z_of_type[type]
        if masses is not None:
            masses[ind] = mass_in[type]
    # convert units
    positions *= unit_convert("distance", units)
    cell *= unit_convert("distance", units)
    if masses is not None:
        masses *= unit_convert("mass", units)
    if velocities is not None:
        velocities *= unit_convert("velocity", units)

    # create ase.Atoms
    at = Atoms(positions=positions,
               numbers=numbers,
               masses=masses,
               cell=cell,
               pbc=[True, True, True])
    # set velocities (can't do it via constructor)
    if velocities is not None:
        at.set_velocities(velocities)
    at.arrays['id'] = ids
    at.arrays['type'] = types
    if travel is not None:
        at.arrays['travel'] = travel
    if mol_id is not None:
        at.arrays['mol-id'] = mol_id
    if mmcharge is not None:
        at.arrays['mmcharge'] = mmcharge

    if bonds is not None:
        for (type, a1, a2) in bonds_in:
            i_a1 = ind_of_id[a1]
            i_a2 = ind_of_id[a2]
            if len(bonds[i_a1]) > 0:
                bonds[i_a1] += ","
            bonds[i_a1] += "%d(%d)" % (i_a2, type)
        for i in range(len(bonds)):
            if len(bonds[i]) == 0:
                bonds[i] = '_'
        at.arrays['bonds'] = np.array(bonds)

    if angles is not None:
        for (type, a1, a2, a3) in angles_in:
            i_a1 = ind_of_id[a1]
            i_a2 = ind_of_id[a2]
            i_a3 = ind_of_id[a3]
            if len(angles[i_a2]) > 0:
                angles[i_a2] += ","
            angles[i_a2] += "%d-%d(%d)" % (i_a1, i_a3, type)
        for i in range(len(angles)):
            if len(angles[i]) == 0:
                angles[i] = '_'
        at.arrays['angles'] = np.array(angles)

    if dihedrals is not None:
        for (type, a1, a2, a3, a4) in dihedrals_in:
            i_a1 = ind_of_id[a1]
            i_a2 = ind_of_id[a2]
            i_a3 = ind_of_id[a3]
            i_a4 = ind_of_id[a4]
            if len(dihedrals[i_a1]) > 0:
                dihedrals[i_a1] += ","
            dihedrals[i_a1] += "%d-%d-%d(%d)" % (i_a2, i_a3, i_a4, type)
        for i in range(len(dihedrals)):
            if len(dihedrals[i]) == 0:
                dihedrals[i] = '_'
        at.arrays['dihedrals'] = np.array(dihedrals)

    at.info['comment'] = comment

    return at

def is_upper_triangular(arr, atol=1e-8):
    """test for upper triangular matrix based on numpy"""
    # must be (n x n) matrix
    assert len(arr.shape)==2
    assert arr.shape[0] == arr.shape[1]
    return np.allclose(np.tril(arr, k=-1), 0., atol=atol)



def convert_cell(cell,pos):
      """
      Convert a parallelepipedal (forming right hand basis)
      to lower triangular matrix LAMMPS can accept. This
      function transposes cell matrix so the bases are column vectors
      """
      cell = np.matrix.transpose(cell[:])

      if not is_upper_triangular(cell):
          # rotate bases into triangular matrix
          tri_mat = np.zeros((3, 3))
          A = cell[:, 0]
          B = cell[:, 1]
          C = cell[:, 2]
          tri_mat[0, 0] = norm(A)
          Ahat = A / norm(A)
          AxBhat = np.cross(A, B) / norm(np.cross(A, B))
          tri_mat[0, 1] = np.dot(B, Ahat)
          tri_mat[1, 1] = norm(np.cross(Ahat, B))
          tri_mat[0, 2] = np.dot(C, Ahat)
          tri_mat[1, 2] = np.dot(C, np.cross(AxBhat, Ahat))
          tri_mat[2, 2] = norm(np.dot(C, AxBhat))

          # create and save the transformation for coordinates
          volume = np.linalg.det(cell)
          trans = np.array([np.cross(B, C), np.cross(C, A), np.cross(A, B)])
          trans /= volume
          coord_transform = np.dot(tri_mat, trans)
          pos = np.dot(coord_transform, pos.transpose())
          pos = pos.transpose()

          return tri_mat, pos
      else:
          return cell, pos


def write_lammps_data(fileobj, atoms,comment=None,append=False,runner=False):
    ''' atoms is the ase atoms object '''
    frame = atoms
    fout = fileobj
    #print('out',fout)
    #sys.exit()

    if type(frame) != list:
        frame = [frame]
    framework = copy.deepcopy(frame)
    framework = framework[0]
    newcell, newpos = convert_cell(framework.cell, framework.positions)
    framework.set_cell(newcell)
    framework.set_positions(newpos)
    species = framework.get_atomic_numbers()
    unique_species, i_unique_species = np.unique(species, return_index=True)
    masses = framework.get_masses()
    positions = framework.get_positions()
    lammps_indices = species.copy()
    for ii in range(len(unique_species)):
        lammps_indices[species == unique_species[ii]] = ii + 1

    ((hxx, hxy, hxz), (hyx, hyy, hyz) , (hzx, hzy, hzz)) = framework.get_cell()


    #fout = open(outfilename,'w')  # is already opened by lammps before
    fout.write("LAMMPS positionsfile from: "+str(comment)+"\n")
    fout.write("\n")
    fout.write(str(len(species))+" atoms\n")

    if runner is False:  # "normal"  lammps format
        fout.write(str(len(unique_species))+" atom types\n")
    elif runner is True:
        # in principle this can be automatically obtained from the nn potential,
        # for now however I am lazy
        fout.write("3 atom types\n")
    else:
        sys.exit('variable runner is neigher True nor False; Exit.')
    fout.write("\n")
    fout.write("0.0 "+str(hxx)+" xlo xhi\n")
    fout.write("0.0 "+str(hyy)+" ylo yhi\n")
    fout.write("0.0 "+str(hzz)+" zlo zhi\n")
    fout.write(str(hxy)+" "+str(hxz)+" "+str(hyz)+" xy xz yz\n")
    fout.write("\n")
    fout.write("Masses\n")
    fout.write("\n")
    if runner is False:
        for ii in range(len(unique_species)):
            #print(lammps_indices[i_unique_species[ii]],
            #masses[i_unique_species[ii]])
            fout.write(str(lammps_indices[i_unique_species[ii]])+" "+str(masses[i_unique_species[ii]])+"\n")
    elif runner is True:
        # in principle this can be automatically obtained from the nn potential,
        # for now however I am lazy
        fout.write("1 24.305\n")        # Mg
        fout.write("2 26.9815385\n")    # Al
        fout.write("3 28.085\n")        # Si                                                   species_runner = species.copy()
        species_runner[species_runner == 12] = 1
        species_runner[species_runner == 13] = 2
        species_runner[species_runner == 14] = 3
        lammps_indices = species_runner
        #print('sp',species_runner)

    fout.write("\n")
    fout.write("Atoms\n")
    fout.write("\n")

      # take care here for runner / n2p2 if only 2 species (say Mg,Si) for Mg,Al,Si      potential
    for ii in range(len(species)):  # if only Mg and Si get is, there will be only     type 1 and 2 and not 1 and 3!
        fout.write(("%5d %5d %s" % (ii + 1, lammps_indices[ii], "  ".join(map(str,     positions[ii]))) )+"\n")
        #fout.write(("%5d %5d %.10f %.10f %.10f" % (ii + 1, lammps_indices[ii], "  ".  join(map(str, positions[ii]))) )+"\n")
    fout.write("\n")
    fout.close()
    #print('written (4)      : ')
    return
