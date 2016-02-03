#!/bin/sh
echo
echo "########### otool -L /usr/local/bin/lammps ##################"
otool -L /usr/local/bin/lammps
echo
echo
echo
echo
echo "######## cd /usr/local/share/lammps/bench"
cd /usr/local/share/lammps/bench
echo
echo "##### starting ####"
echo "##### if this is not working do (even if you started the seriell version!): brew link open-mpi"
echo
echo
echo "##### which lammps ######"
which lammps
echo
echo
echo
echo
echo
echo "##### mpiexec -n 1 lammps -i in.lj #####"
mpiexec -n 1 lammps -i in.lj
echo "##### done #####"
echo
echo
echo
echo

echo "##### lammps -in in.lj #####"
lammps -in in.lj
echo "##### done #####"
