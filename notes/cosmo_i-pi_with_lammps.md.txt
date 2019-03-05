cd /home/glensk/google_drive/Glensk/KMC-Tutorial/test
../../../scripts_epfl/i-pi-mc/bin/i-pi input-runner.xml&> log &
cat log # -->(need to install numpy with .... probably change to env the i-pi script)
lmp_serial < in.nn  # this uses the wrong lammps (not the one compiled with runner)

study ~/google_drive/scripts_epfl/i-pi-mc/ipi/engine/motion/dynamics.py




cd /home/glensk/google_drive/Glensk/KMC-Tutorial/test_3
i-pi input-runner.xml&> log &
/home/glensk/google_drive/scripts_epfl/lammps_cosmo/lammps_cosmopc_20180821/lammps/src/lmp_serial_cosmopc_runner-lammps2 < in.nn

