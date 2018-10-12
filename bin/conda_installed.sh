#!/bin/sh
conda install -c conda-forge pyfftw=0.10.4
conda install -c conda-forge lammps                 # lammps for python
conda install -c conda-forge phonopy                # phonopy  (if manyally installed needs to be put into PYTHONPATH manually)
