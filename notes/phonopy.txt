#######################################################
# Setyp phonopy:
#######################################################
is PYTHONPATH set correctly?: echo $PYTHONPATH you have to add (lib64/python and/or lib/python)


#######################################################
# What does not work
#######################################################
kann leider keine fcc4 zelle vervielfachen, kommt bloedsinn raus
schlaegt falsche auslenkungen bei der fcc4 zelle vor (check /direct/cartesian)
- KANN DIREKT ABER NICHT CARTESISCH!!!!!!!!!!!!!!!!!!!!
- KANN NICHT negative volumen verarbeiten!


#######################################################
# run phonopy (al)
#######################################################
phonopy -d --dim="1 1 1" --amplitude .015875352   (=0.03 bohrradius) (needs POSCAR)
phonopy -d --dim="1 1 1" --amplitude .0052917721  (=0.01 bohrradius) (needs POSCAR)
phonopy -d --dim="1 1 1" --amplitude .0052917721 --nodiag  (to create disp only in xyz)
phonopy -d --dim="1 1 1" --amplitude .0052917721 --tolerance 0.0001

phonopy -d --dim="2 2 1" --amplitude 0.0052917721                    (needs POSCAR)
phonopy -f 4.05Ang/vasprun.xml
phonopy --fz job-000/vasprun.xml ../job-0*/vasprun.xml (to substract background forces)
make mesh.conf file:
    ATOM_NAME = Al
    DIM = 1 1 1
    MP = 8 8 8
    FORCE_CONSTANTS = WRITE
phonopy -p mesh.conf  (--tolerance 0.0001) -> erstellt FORCE_CONSTANTS (==HesseMatrix)
phonopy -t mesh.conf

a) create displacements

        - atom:    1
          direction:
            [   1.0000000000000000,  0.0000000000000000,  0.0000000000000000 ]
          displacement:
            [   0.0052917721000000,  0.0000000000000000,  0.0000000000000000 ]


            --> disp.yaml displacement: [vector ] ist das displacement in 
            carteschischen koordinaten, absolute auslenkung des Atoms


c) create FORCE_SETS
            phonopy -f 1-1-pos_shift/vasprun.xml disp-002/vasprun.xml ... this needs disp.yaml!
            --> the displacement is taken directly from disp.yaml!! not from inputfile!


b) erstelle mesh.cof
                        was auch geht: erstelle file settings:
                        DIM = 2 2 2
                        PM = .TRUE.
                        phonopy -d settings

b) do calculation


d) create THERMO
            phonopy -t mesh.conf
                        (needs primitive POSCAR, Fs/glensk/v/pp/al/dynmat_bulk_fcc4/2x2x2sc_450eV-NGXF360-ADDGRID_10x10x10kp-0x0x0shift_0.0eVsmear0.2/fitperatom/getThermodynamics_Mesh/'ORCES_SETS, DOESNT NEED disp.yaml :))) )
                        1. doesn-t need vasprun.xml
                        2. no matter wich POSCAR.PC i put in --> it gives the same Fene -> seems that just FORCES_SETS is read 

INSTALLATION:
1. entpacke phonopy-0.9.5 mit:: tar xfvz phonopy-0.9.5.tar.gz 
2. create lib and lib64 by going to phonopy folder (should contain setup.py):: python setup.py install --home=.
3. in tcsh add: setenv PYTHONPATH /home/glensk/progs/phonopy/phonopy-0.9.4.2/lib64/python
4. in tcsh add: set path = ($HOME/progs/phonopy/phonopy-0.9.4.2/scripts $path)

phonopy -help
--writefc -> gives HessianMatrix

Format of Forces sets:
1. displacement in angstrom  -> the displacement
2. Forces in eV/Angstrom^2 from vasprun
0.01588 == displacement in Angstrom || -0.071773 eV/Angstrom^2 (force)
