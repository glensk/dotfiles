
/home/grabowski/projects/implementingTI_dynMat/vasp



Hessematrix_sphinx  hartree/bohr2
Hessematrix		ev/ang2

EqCoords_direct ungestoerte struktur
XXX geth nicht :  EqCoords_cartesian 	cartesiche koordinaten

PRE_EQ_N --> harmonic eq steps (500 oder 1000)

POTIM zeitschritt in fs
IBRION 0	= MD


/usr/mpi/intel/openmpi-1.2.8/bin/mpirun -np 8 /san.backup/grabowski/VASP_src.4.6/parallel_Langevin/vasp


parallel:
/san.backup/grabowski/openmpi-1.4.1_libsAndBins/bin/mpirun -np 8 /san.backup/grabowski/VASP_src.4.6/parallel_Langevin/vasp

seriell:
/san.backup/grabowski/VASP_src.4.6/serial_forExperimenting/vasp



geaenderte dateien
stepver.F langevin
thermodynamicIntegration_reference.F	main.F thermodyn integration
reader.F einleseroutine fuer die Incar
