!!!!!! rund lammpscheck.sh if problems with lammps
!!!!!! try brew link --overwrite open-mpi  (solved the problem last time) -> makes parallel work but not serial fully ...
!!!!!! try to remove anaconda from path if problems with lammps

HOW TO START LAMMPS:
in /Users/glensk/lammps_tmp/lammpsRuns
either in ipython: run demo.py 


or from command line: 
lammps -i in_file.in



From Homebrew:
%brew install lammps
==> Installing lammps from homebrew/homebrew-science
==> Downloading https://homebrew.bintray.com/bottles-science/lammps-2015.02.10_1.yosemite.bottle.tar.gz
Already downloaded: /Library/Caches/Homebrew/lammps-2015.02.10_1.yosemite.bottle.tar.gz
==> Pouring lammps-2015.02.10_1.yosemite.bottle.tar.gz
==> Caveats
You should run a benchmark test or two. There are plenty available.

  cd /usr/local/share/lammps/bench
  lammps -in in.lj
  # with mpi
  mpiexec -n 2 lammps -in in.lj

The following directories could come in handy

  Documentation:
  /usr/local/share/lammps/doc/Manual.html

  Potential files:
  /usr/local/share/lammps/potentials

  Python examples:
  /usr/local/share/lammps/python-examples

  Additional tools (may require manual installation):
  /usr/local/share/lammps/tools

To use the Python module with Python, you may need to amend your
PYTHONPATH like:
  export PYTHONPATH=/usr/local/lib/python2.7/site-packages:$PYTHONP



# To check if correctly linked
otool -L /usr/local/bin/lammps

otool == ldd (on linux)

env (gives environment vars)


##########################################################
    INPUT FILE
##########################################################

dump dump1 all xyz 10 trj_lammps.out   # usual output
dump_modify dump1 format "%s %4.7f %4.7f %4.7f"   # modifies output to have more digits
