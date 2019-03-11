1514  mcd giulio_Al6xxx
 1515  mv ../input.data .
 1517* dot
 1518* cd scripts
 1519* cd runner_scripts
 1520* cp inputAlMgSi.nn /home/glensk/giulio_Al6xxx
 1522* cd examples
 1523* cd nnp-train/H2O_RPBE-D3
 1524* cd ../Cu2S_PBE
 1525* cp input.nn ~/giulio_Al6xxx
 1528* cd ../H2O_RPBE-D3
 1534* cd ..
 1535* cd bin
 1539  vi -O input.nn inputAlMgSi.nn
 1540  grep symfunction_short $dotfiles/scripts/potentials/runner_lammps_nn/v2dg/input.nn
 1541  grep "^symfunction_short" $dotfiles/scripts/potentials/runner_lammps_nn/v2dg/input.nn
 1542  grep "^symfunction_short" $dotfiles/scripts/potentials/runner_lammps_nn/v2dg/input.nn >> input.nn
 1544  /home/glensk/Dropbox/Albert/git/n2p2/bin/nnp-scaling 1
 1545* top
 1546  vi nnpscaling.log
 1548  wc -l nnp*
 1550  cat ~/Dropbox/Albert/scripts/dotfiles/commands/que.commands
 1551  slurm -h
 1552  slurm --help
 1554  rm *histo
 1556  vi nnp-scaling.log.0000
 1558  mpirun --help
 1560  head input.data
 1561  mpirun -np 4 /home/glensk/Dropbox/Albert/git/n2p2/bin/nnp-train
 1562  ml
 1563  module list
 1565  export SET_NUM_THREADS=4
 1567  export OMP_NUM_THREADS=4
 1568  /home/glensk/Dropbox/Albert/git/n2p2/bin/nnp-train
 1569  vi input.nn
 1570  ls
 1571  head scaling.data
 1572  head -20 scaling.data

