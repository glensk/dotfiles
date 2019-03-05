cd ~/google_drive/scripts_epfl/lammps_cosmo
git clone https://github.com/cosmo-epfl/lammps.git
cd lammps
git checkout runner-lammps
cd src
make ps
make yes-CLASS2 yes-KSPACE yes-MANYBODY yes-MISC yes-MOLECULE yes-REPLICA yes-RIGID yes-USER-MISC
make yes-USER-RUNNER
module load intel intel-mpi intel-mkl fftw python/2.7.14  # only in case of parrel (fidis) compiilation
make serial  or make daneb
(fidis has 28 and deneb has 16 and 24 depending on partition)


cd /Users/glensk/Dropbox/Albert/scripts/lammps_cosmo/first_input_natascha/test_1
../../lammps/src/lmp_serial -i in.lmp


data.lmp: structure (fcc 4x4x4 supercell with 216 atoms)
scaling.data and weiths.028xxx and input.nn define the NN potential 

-----------------------------------------------------------------

ERROR: neighbor list array too small, increase MAXNEIGH and recompile.
--> change in pair_runner.h the #define MAXNEIGH from 192 to 500 in line 26 in the scr folder of lammps

-----------------------------------------------------------------
then, to make lammps work as a python library (correctly) (miniconda version was not working correctly)

# to make the python libraries:
cd /Users/glensk/Dropbox/Albert/Google_Drive/scripts_epfl/lammps_cosmo/lammps_macbook_20180814/src
make mpi-stubs
make g++_serial mode=shlib

Export environment variables:
export LAMMPS_COMMAND="/usr/local/bin/lmp_serial" (path to binary file)
export LAMMPSPATH="/Users/lopanits/source/lammps" (path to lammps folder)
export PYTHONPATH="$LAMMPSPATH/python:$PYTHONPATH"
export LD_LIBRARY_PATH="$LD_LIBRARY_PATH:$LAMMPSPATH/src"
export GPAW_SETUP_PATH="$GPAW_SETUP_PATH:/Users/lopanits/source/gpaw" (in case of using GPAW calculator, not needed for lammps)
