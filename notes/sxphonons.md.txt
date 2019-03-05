     module load sphinx/serial/2.0.4
### setting up structure and displacements with sphinx
# generating supercell
     sxstructprint --vasp -i POSCAR > struct.sx
     sxstructrep -r 2x2x2 -i struct.sx -o Sstruct.sx
     sx2poscar -i Sstruct.sx -o SPOSCAR
# generating displacements
     sxuniqdispl -d 0.02 -i Sstruct.sx (displacement of 0.02 bohrradius usually ok)
     sx2poscar -i input-disp-1.sx -o SPOSCAR1
     sx2poscar -i input-disp-2.sx -o SPOSCAR2
     ...
     cp SPOSCAR background_forces/POSCAR (dont forget background forces (sometimes zero!)
# computing stuff

# evaluate
get_sxdynmat.sx

/home/grabowski/sphinx/trunk/src/add-ons/sxdynmat -i sxdynmat.sx
/home/grabowski/sphinx/trunk/src/add-ons/sxdynmat -i sxdynmat.sx -H
                          or  #  sxdynmat -i sxdynmat.sx
                          or  #  sxphonon -s phononSet.sx 

For thermodynamics:
#  sxthermo -T 2000 -dT 2 -p phonon.sxb 

