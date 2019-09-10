#!/usr/bin/env python
# -*- coding: utf-8 -*-
from __future__ import print_function
import os,sys,re,glob
import argparse
import textwrap
import numpy as np
from subprocess import call
from tempfile import mkstemp
from shutil import move
from shutil import copyfile
from os import fdopen, remove
import my_atom
import get_parametrization_for_displacementfolder as get_param_disp
import myutils

def help(p = None):
    string = '''

    - to get u_OUTCAR you can use OUTCAR_ene-potential_energy.sh OUTCAR(.gz) > u_OUTCAR --> better use u_OUTCAR_without_first_substract for catted jobs
      OUTCAR_ene-potential_energy_without_first_substracted.sh > u_OUTCAR

    TODO: --------------------------------------------------------------------------------
    - for checking forces check: /Users/glensk/Dropbox/proj/proj_current/__2017.01_phonon_linewidth_al/__2016.09_phonon_lifetimes_3_nach_elternzeit/michaels_code_simulation_morse_withforces_oldtonew_samesteps_tox0/old_to_new/check_forces_lon_tox2
    - running: check la md for tow different atomic positions (load from POSITIONS and made from skript itself) and check lifetimes with both -> running -> make 500000 steps but write only every 10th
            answer following questions:
            - can a la md without tox reproduce the faltung of DFT?
            - is the anharmonic hamiltonian relevant at all the get faltung right? Is any amount of anharmonicity ok?
    - if deconvolution works for la without tox, check it for all LA potentials for all elements, if its thatn still the same for all elements we can argue that its a universal feature;
    - write mail an michael (see fragen omnifocus)


    - check the harmonic and LA forces to old python sktrips!
    - make ene std work with TOX
    - go quickly through all elements for a given la hamiltonian and get std ...\
    -  .... this has the advantage that std (for energy and forces) can be quickly checked for a certain parametrization or version of LA
    - get the std in energy (about 0.5 meV/atom and forces of the TUTILD Al potential)


    Examples: ----------------------------------------------------------------------------
    export basefolder=/Users/glensk/Dropbox/Albert/Understanding_distributions/
    export testfolder_al_md=$basefolder/ti/Al/30__PTS_dosall_from_2x2x2sc__LON_quer_3x3x3kp_all_morse_LONADD__NO_LON2_NONE_LON2ADD__NO_TOX_NONE_TOY_NONE_TOZ_NONE/lambda1.0/
    export testfolder_al_disp=$basefolder/displacements_/Al/5x5x5sc_quer_2x2x2kp_vasp4

    - argparse_md_code.py -alat_mor 4.13 -a_mor 1.502257 -D_mor 0.238955 -rp $basefolder/ti/ti_al_2x2x2sc/11_newskript_long_tox_redo_WORKS/lambda0.0/tmp/cartesian_coords
    - argparse_md_code.py -rp POSITIONs.6 -alat_mor 4.13 -a_mor 1.502257 -D_mor 0.238955
    - argparse_md_code.py -t1 -0.09 -rp POSITIONs.6
    - argparse_md_code.py -f $testfolder_al_md -v -wa -a_mor 1.30     # this way morse parameters can be changed
      argparse_md_code.py -rp $testfolder_al_md/POSITIONs -a_mor 1.30 -v -wa  # this way morse paramters can be changed.
    - argparse_md_code.py -f $testfolder_al_md -v -wa
          (gives correct std and dudl as python scripts -> $basefolder/ti/summary_al; change to lambda0.0 to see corresponding result)
    - argparse_md_code.py -rp $testfolder_al_disp/POSITIONs -v -D_mor 0.2793 -a_mor 1.432673 -a 4.13 -N 5 -waf
    - argparse_md_code.py  -f $basefolder/ti/ti_original_hesse/Al/lambda0.0_20604/ -v -wa  (for lambda 0.0 and 1.0 gives same ene std harm)

    - argparse_md_code.py -ea --f_md_2x2x2_30 -dbl -vm2   # parametrization from the 2x2x2 supercell displacements (I think)
	==> (Al) ene_std:   8.8 for_std:  0.157  dudl/2:   7.4  ||| alat 4.130 ||| alat_mor 4.130 a_mor 1.43267  D_mor 0.27930 ||(PRL NO )
	==> (Ag) ene_std:  16.6 for_std:  0.447  dudl/2:  26.8  ||| alat 4.310 ||| alat_mor 4.250 a_mor 1.86396  D_mor 0.18942 ||(PRL NO )
	==> (Au) ene_std:  19.4 for_std:  0.299  dudl/2:  24.9  ||| alat 4.250 ||| alat_mor 4.250 a_mor 1.86396  D_mor 0.18942 ||(PRL YES)
	==> (Cu) ene_std:   2.1 for_std:  0.046  dudl/2:   0.5  ||| alat 3.750 ||| alat_mor 3.750 a_mor 1.96836  D_mor 0.13183 ||(PRL YES)
	==> (Ir) ene_std:   9.7 for_std:  0.319  dudl/2:   6.6  ||| alat 3.990 ||| alat_mor 3.990 a_mor 1.80552  D_mor 0.33313 ||(PRL NO )
	==> (Pb) ene_std:   2.9 for_std:  0.060  dudl/2:   1.8  ||| alat 5.130 ||| alat_mor 5.130 a_mor 1.56537  D_mor 0.06893 ||(PRL NO )
	==> (Pd) ene_std:   9.4 for_std:  0.121  dudl/2:   4.0  ||| alat 4.100 ||| alat_mor 4.100 a_mor 1.82081  D_mor 0.17966 ||(PRL NO )
	==> (Pt) ene_std:  31.8 for_std:  0.421  dudl/2:  34.8  ||| alat 4.100 ||| alat_mor 4.100 a_mor 1.78089  D_mor 0.28964 ||(PRL NO )
	==> (Rh) ene_std:   6.1 for_std:  0.201  dudl/2:   3.2  ||| alat 3.980 ||| alat_mor 3.980 a_mor 1.78328  D_mor 0.22483 ||(PRL YES)

    - argparse_md_code.py -ea --f_md_2x2x2_30 -dbl -vm2 -gpf_dispfolder_5x5x5   # parametrization from the 5x5x5 supercells
        ==> (Al) ene_std:   6.5 for_std:  0.110  dudl/2:   5.3  ||| alat 4.130 ||| alat_mor 4.130 a_mor 1.51085  D_mor 0.22675 ||(PRL YES)
        ==> (Ag) ene_std:   2.8 for_std:  0.063  dudl/2:   3.4  ||| alat 4.310 ||| alat_mor 4.310 a_mor 2.02203  D_mor 0.07629 ||(PRL NO )
        ==> (Au) ene_std:  17.8 for_std:  0.272  dudl/2:  22.7  ||| alat 4.250 ||| alat_mor 4.250 a_mor 1.91739  D_mor 0.17024 ||(PRL YES)
        ==> (Cu) ene_std:   1.9 for_std:  0.045  dudl/2:   0.5  ||| alat 3.750 ||| alat_mor 3.750 a_mor 1.99523  D_mor 0.12707 ||(PRL YES)
        ==> (Ir) ene_std:  10.6 for_std:  0.341  dudl/2:   5.9  ||| alat 3.990 ||| alat_mor 3.990 a_mor 1.98876  D_mor 0.22729 ||(PRL YES)
        ==> (Pb) ene_std:   2.1 for_std:  0.053  dudl/2:   1.3  ||| alat 5.130 ||| alat_mor 5.130 a_mor 1.72439  D_mor 0.04996 ||(PRL YES)
        ==> (Pd) ene_std:   8.3 for_std:  0.110  dudl/2:   3.5  ||| alat 4.100 ||| alat_mor 4.100 a_mor 1.85532  D_mor 0.16740 ||(PRL YES)
        ==> (Pt) ene_std:  26.2 for_std:  0.346  dudl/2:  28.4  ||| alat 4.100 ||| alat_mor 4.100 a_mor 1.86887  D_mor 0.23838 ||(PRL NO )
        ==> (Rh) ene_std:   6.8 for_std:  0.217  dudl/2:   3.2  ||| alat 3.980 ||| alat_mor 3.980 a_mor 1.89434  D_mor 0.17953 ||(PRL YES)

    - argparse_md_code.py -ea --f_md_2x2x2_30 -dbl -vm2 -gpf_dispfolder_5x5x5_alat_mor_T0K  # seems a mixed blessing, sometimes works, sometimes not...
  	==> (Al) ene_std:   4.0 for_std:  0.104  dudl/2:   1.9  ||| alat 4.130 ||| alat_mor 4.054 a_mor 1.51085  D_mor 0.22675 ||(PRL YES)+
  	==> (Ag) ene_std:   7.1 for_std:  0.240  dudl/2:   9.9  ||| alat 4.310 ||| alat_mor 4.170 a_mor 2.02203  D_mor 0.07629 ||(PRL NO )--
  	==> (Au) ene_std:  11.3 for_std:  0.185  dudl/2:  11.8  ||| alat 4.250 ||| alat_mor 4.170 a_mor 1.91739  D_mor 0.17024 ||(PRL YES)+
  	==> (Cu) ene_std:   7.0 for_std:  0.260  dudl/2:   0.3  ||| alat 3.750 ||| alat_mor 3.644 a_mor 1.99523  D_mor 0.12707 ||(PRL YES)--
  	==> (Ir) ene_std:  21.7 for_std:  0.708  dudl/2:   4.7  ||| alat 3.990 ||| alat_mor 3.880 a_mor 1.98876  D_mor 0.22729 ||(PRL YES)--
  	==> (Pb) ene_std:   3.3 for_std:  0.086  dudl/2:   0.1  ||| alat 5.130 ||| alat_mor 5.035 a_mor 1.72439  D_mor 0.04996 ||(PRL YES)-
  	==> (Pd) ene_std:   8.0 for_std:  0.312  dudl/2:   1.7  ||| alat 4.100 ||| alat_mor 3.985 a_mor 1.85532  D_mor 0.16740 ||(PRL YES)-
  	==> (Pt) ene_std:  14.3 for_std:  0.358  dudl/2:   8.1  ||| alat 4.100 ||| alat_mor 3.980 a_mor 1.86887  D_mor 0.23838 ||(PRL YES)-
  	==> (Rh) ene_std:  20.0 for_std:  0.569  dudl/2:   3.9  ||| alat 3.980 ||| alat_mor 3.840 a_mor 1.89434  D_mor 0.17953 ||(PRL NO )-

    - argparse_md_code.py -ea --f_md_2x2x2_30 -dbl -vm2 -gpf_dispfolder_5x5x5_alat_mor_current_alat  # consistent lattice with morse
    - argparse_md_code.py -ea --f_md_2x2x2_30 -dbl -vm2 -gpf_dispfolder_5x5x5 -gpf_dispfolder_cf110 999
    ==> (Al) ene_std:   7.1 for_std:  0.121  dudl/2:   5.8  ||| alat 4.130 ||| alat_mor 4.130 a_mor 1.48337  D_mor 0.24190 ||(PRL YES)o
    ==> (Ag) ene_std:   2.5 for_std:  0.058  dudl/2:   2.7  ||| alat 4.310 ||| alat_mor 4.310 a_mor 2.04471  D_mor 0.07274 ||(PRL NO )+
    ==> (Au) ene_std:  16.0 for_std:  0.239  dudl/2:  20.1  ||| alat 4.250 ||| alat_mor 4.250 a_mor 1.96352  D_mor 0.15298 ||(PRL YES)+  (std 6.7 possible)
    ==> (Cu) ene_std:   1.9 for_std:  0.044  dudl/2:   0.5  ||| alat 3.750 ||| alat_mor 3.750 a_mor 2.00776  D_mor 0.12481 ||(PRL YES)+
    ==> (Ir) ene_std:   9.4 for_std:  0.311  dudl/2:   6.2  ||| alat 3.990 ||| alat_mor 3.990 a_mor 1.94468  D_mor 0.25800 ||(PRL YES)+
    ==> (Pb) ene_std:   3.8 for_std:  0.074  dudl/2:   2.2  ||| alat 5.130 ||| alat_mor 5.130 a_mor 1.56174  D_mor 0.07429 ||(PRL NO )-
    ==> (Pd) ene_std:   6.8 for_std:  0.099  dudl/2:   2.6  ||| alat 4.100 ||| alat_mor 4.100 a_mor 1.89152  D_mor 0.15303 ||(PRL YES)+
    ==> (Pt) ene_std:  20.2 for_std:  0.273  dudl/2:  21.2  ||| alat 4.100 ||| alat_mor 4.100 a_mor 1.97067  D_mor 0.18950 ||(PRL YES)+
    ==> (Rh) ene_std:   6.0 for_std:  0.201  dudl/2:   3.1  ||| alat 3.980 ||| alat_mor 3.980 a_mor 1.86984  D_mor 0.19434 ||(PRL YES)+

    NOTES:
    ------
    - argparse_md_code.py -e Au --f_md_2x2x2_30 -dbl -vm2 -gpf_dispfolder_5x5x5 -gpf_dispfolder_cf110 999 -D_mor 0.090
    ==> (Au) ene_std:   6.7 for_std:  0.292  dudl/2:   0.9  ||| alat 4.250 ||| alat_mor 4.250 a_mor 1.96352  D_mor 0.09000 ||(PRL YES)


    POSITIONS @ MG 4.13:  ( -> Hmmm.... seeems to be way better when alat_mor == alat_lattice )
    --------------------
    Reference: Parametrization @ 4.14 -> shifted to @ 4.13
    - argparse_md_code.py -e Al --f_md_2x2x2_30 -dbl -vm2 --gpf_dispfolder_setfolder Al/3x3x3sc_4.14Ang_quer_10x10x10kp_vasp4_ENCUT400 --gpf_dispfolder_shift_parametrization_to_alat 4.13
   ==> (Al) ene_std:   6.3 for_std:  0.106  dudl/2:   5.1  ||| alat 4.130 ||| alat_mor 4.130 a_mor 1.51876  D_mor 0.22239 ||(PRL YES)

    - argparse_md_code.py -e Al --f_md_2x2x2_30 -dbl -vm2 --gpf_dispfolder_setfolder Al/3x3x3sc_4.14Ang_quer_10x10x10kp_vasp4_ENCUT400 --gpf_dispfolder_shift_parametrization_to_alat 4.04
    ==> (Al) ene_std:   8.7 for_std:  0.143  dudl/2:   7.0  ||| alat 4.130 ||| alat_mor 4.040 a_mor 1.44702  D_mor 0.32937 ||(PRL NO )

    Parametrization (original) @ 4.14
    - argparse_md_code.py -e Al --f_md_2x2x2_30 -dbl -vm2 --gpf_dispfolder_setfolder Al/3x3x3sc_4.14Ang_quer_10x10x10kp_vasp4_ENCUT400 --gpf_dispfolder_shift_parametrization_to_alat 4.14
   ==> (Al) ene_std:   6.1 for_std:  0.103  dudl/2:   4.9  ||| alat 4.130 ||| alat_mor 4.140 a_mor 1.52566  D_mor 0.21331 ||(PRL YES)

    Parametrization (original) @ alt = 4.04
    argparse_md_code.py -e Al --f_md_2x2x2_30 -dbl -vm2 --gpf_dispfolder_alat_lattice_T0K
    ==> (Al) ene_std:   7.8 for_std:  0.126  dudl/2:   6.2  ||| alat 4.130 ||| alat_mor 4.040 a_mor 1.46571  D_mor 0.30986 ||(PRL YES)

    Parametrization (original) @ alt = 4.04 -> shifted to @ 4.13
    # THIS IS EXACTLY THE PARAMETRIZATION I LIKE TO BE BEST! (so far)
    - argparse_md_code.py -e Al --f_md_2x2x2_30 -dbl -vm2 --gpf_dispfolder_alat_lattice_T0K --gpf_dispfolder_shift_parametrization_to_alat 4.13
    ==> (Al) ene_std:   5.9 for_std:  0.098  dudl/2:   4.7  ||| alat 4.130 ||| alat_mor 4.130 a_mor 1.52024  D_mor 0.21700 ||(PRL YES)

    Parametrization (original) @ alt = 4.04 -> shifted to 4.05
    - argparse_md_code.py -e Al --f_md_2x2x2_30 -dbl -vm2 --gpf_dispfolder_alat_lattice_T0K --gpf_dispfolder_shift_parametrization_to_alat 4.05
    ==> (Al) ene_std:   7.6 for_std:  0.122  dudl/2:   6.0  ||| alat 4.130 ||| alat_mor 4.050 a_mor 1.47257  D_mor 0.29737 ||(PRL YES)

    POSITIONs @ DISP 4.14;
    ---------------------
    - argparse_md_code.py -e Al -f /Users/glensk/Dropbox/Albert/Understanding_distributions/displacements_/Al/3x3x3sc_4.14Ang_quer_10x10x10kp_vasp4_ENCUT400 -vm2 --gpf_dispfolder_alat_lattice_T0K
    ==> (Al) ene_std:   0.3 for_std:  0.014  dudl/2:   0.0  ||| alat 4.140 ||| alat_mor 4.040 a_mor 1.46571  D_mor 0.30986 ||(PRL -  )

    - argparse_md_code.py -e Al -f /Users/glensk/Dropbox/Albert/Understanding_distributions/displacements_/Al/3x3x3sc_4.14Ang_quer_10x10x10kp_vasp4_ENCUT400 -vm2 --gpf_dispfolder_alat_lattice_T0K --gpf_dispfolder_shift_parametrization_to_alat 4.14
    ==> (Al) ene_std:   0.2 for_std:  0.011  dudl/2:   0.0  ||| alat 4.140 ||| alat_mor 4.140 a_mor 1.52538  D_mor 0.20895 ||(PRL -  )


    NEXT1: CHECK IF THIS IS CONSISTENT WITH FORCES/ENERGIES From the displacements
    NEXT2: for --f_md_2x2x2_30 get morse parametrization with alat_mor 4.05

    NOTES: -------------------------------------------------------------------------------
    - this python script is only a wrapper to start md_lon_tox.c
    - currently supported hamiltonians  : LA / Harmonic
    - currently supported structures    : fcc qubic
    - currently used mass               : Al: 26.982*1.660538e-27

    - md_lon_tox.c can be used in 3 modi (a,b,c):
        a) -md: run MD (molecular dynamics) using a LA/harmonic hamiltonian
                - needs no further input but the executable;
                - if on harmonic potential: needs a Hessematrix_sphinx or a FORCE file from phonopy

        b) -ti: perform Thermodynamic Integration (TI) from harmonic reference to LA
                - needs a Hessematrix_sphinx or a FORCE file from phonopy

        c) -rp: (default) to read a positionsfile and calculate the corresponding forces/energies
                this corresponds to a thermodynamic integration on the corresponding potential where the
                POSITIONs were made
                - needs a positionsfile (POSITIONs) which may additionally have forces
                - can read u_OUTCAR{,_first_substract} to get u_DFT
                - car read u_temp to get temperature
                - can read dUdL to get u_DFT and temperature
                - can read in Hessematrix_sphinx

        analyze files:
            - out_analyze_forces_vs_dft__is__fig3b_prl.dat
              plots the LA forces vs the DFT forces as in Fig 3b of anharmonicity prl
    '''
    p = argparse.ArgumentParser(description=string,
            formatter_class=argparse.RawTextHelpFormatter) #ArgumentDefaultsHelpFormatter)

    p.add_argument('-e',   '--element', required=False,
       help='Element used to set the mass', type=str, default="Al")

    p.add_argument('-ea',   '--element_all' , default=False,required=False,action='store_true',
       help='Go over all available fcc elements')

    p.add_argument('-em',   '--element_mass', required=False,
       help=argparse.SUPPRESS, default=26.9815386000)

    p.add_argument('-sc',   '--supercell', '-N' , required=False,
       help='supercell size (default = 2) -> 2*2*2*4 = 32 atoms for fcc', type=int, default=False)
    p.add_argument('-t',   '--temperature', '-T' , required=False,
       help='starting temperature, approximately 2*real temperature (Ekin+Epot)', type=float, default=0)
    p.add_argument('-dt',   '--timestep' , required=False,
       help='timestep in [fs] of the calculation when performing molecular dynamics (MD) (default = 1 [fs])', type=int, default=1)
    p.add_argument('-s','-steps', '-mdsteps',   '--steps', required=False,
       help='number of ionic steps', type=int, default=False)
    p.add_argument('-a',   '--alat', required=False,
       help='lattie constant in Angstrom (default = False) (Al~4.05-4.14)', type=float, default=False)
    p.add_argument('-l',   '--write_every' , required=False,
       help='output positions every l simulation steps', type=int, default=1)
    p.add_argument('-eh',   '--evolve_harmonic' , required=False,action='store_true',
       help='run/evolve the MD on the harmonic Hessematrix instead of the la potentila (lambda=0.0)')

    ### specify folder to import POSITIONs from
    p.add_argument('-fmd30',   '--f_md_2x2x2_30' , default=False,required=False,action='store_true',
       help='run/evolve the MD on from the displacements/ti folder for set element using 30__PTS_dosall')
    p.add_argument('--f_d_5x5x5','-f_d_5x5x5' , default=False,required=False,action='store_true',
       help='run/evolve the MD on from the displacements folder for set element in the 5x5x5 supercell')
    p.add_argument('-f',   '--folder', required=False,
       help='Define folder to look for POSITIONs or HesseMatrix_sphinx file', type=str, default=False)
    p.add_argument('-fandp',   '--folder_and_parametrization', required=False,
       help='Define folder to look for POSITIONs or HesseMatrix_sphinx file and parametrization file as disp_fit.parameters_morse.dat', type=str, default=False)


    p.add_argument('-rp',   '--read_positions', required=False,
       help='if positions file to read in e.g. -rp POSITIONs ', type=str, default=False)
    p.add_argument('-rp0',   '--read_positions0', required=False,
       help='if positions of the undisplaced structure', type=str, default=False)
    p.add_argument('-rh',   '--read_hesse', required=False,
       help='if HesseMatrix_sphinx file to read in e.g. -rh HesseMatrix_sphinx', type=str, default=False)
    p.add_argument('-rf',   '--read_forces', required=False,
       help=argparse.SUPPRESS, default=False)
    p.add_argument('-ru',   '--read_uoutcar', required=False,
       help=argparse.SUPPRESS, default=False)
    p.add_argument('-cef',   '--calculate_energy_and_forces', required=False,
       help=argparse.SUPPRESS, default=False)

    p.add_argument('-dbl',   '--do_both_lambdas', required=False, action='store_true',
       help=argparse.SUPPRESS, default=False)


    p.add_argument('-sa',          '--sweep_a_mor', required=False, default=False, action='store_true',
       help='sweep a_mor variable and return for_std as funcion of a_mor')
    p.add_argument('-sd',          '--sweep_D_mor', required=False, default=False, action='store_true',
       help='sweep D_mor variable and return for_std as funcion of D_mor')
    p.add_argument('-gpf_dispfolder_setfolder',          '--gpf_dispfolder_setfolder', required=False, default=False, type=str,
       help=argparse.SUPPRESS)
    p.add_argument('-gpf_dispfolder_5x5x5',              '--gpf_dispfolder_5x5x5', required=False, action='store_true',
       help=argparse.SUPPRESS, default=False)
    p.add_argument('-gpf_dispfolder_5x5x5_weighted',     '--gpf_dispfolder_5x5x5_weighted', required=False, action='store_true',
       help=argparse.SUPPRESS, default=False)
    p.add_argument('-gpf_dispfolder_5x5x5_alat_mor_T0K', '--gpf_dispfolder_5x5x5_alat_mor_T0K', required=False, action='store_true',
       help=argparse.SUPPRESS, default=False)
    p.add_argument('-gpf_dispfolder_5x5x5_alat_mor_current_alat', '--gpf_dispfolder_5x5x5_alat_mor_current_alat', required=False, action='store_true',
       help=argparse.SUPPRESS, default=False)
    p.add_argument('-gpf_dispfolder_alat_lattice_T0K',   '--gpf_dispfolder_alat_lattice_T0K', required=False, action='store_true',
            help='get a parametrisation which has been done @T=0K alat', default=False)
    p.add_argument('-gpf_dispfolder_spta', '--gpf_dispfolder_shift_parametrization_to_alat', required=False, type=float, default=False,
             help="at another atom distance, the morse would have a static force which cancels out in equilibrium. Find the static force for the chosed alat_lattice and shift the forces to parametrize correspondingly.")
    p.add_argument('-gpf_dispfolder_cf110', '--gpf_dispfolder_correct_for_110_forces', required=False, type=float, default=False,
             help="correct forces on 05_05_0 atom by the forces on 1_1_0 atom.")



    p.add_argument('-te',   '--test',choices=list(range(1,8))+list(range(20,28))+list(range(40,46)), type=int, default=False,
       help=textwrap.dedent('''\
               Use predefined folder
               1 : MD on POSITIONs (1 structure)  with    tox
               2 : MD on POSITIONs (1 structure)  without tox
               3 : MD on POSITIONs (1 structure)  with    tox
               4 : MD on POSITIONs (1 structure)  without tox 3x3x3 supercell
               5 : MD on POSITIONs (1 structure)  with    tox 3x3x3 supercell
               6 : MD on POSITIONs (1 structure)  with    tox 3x3x3 supercell, other alat
               7 : MD on POSITIONs (1 structure)  with    tox 2x2x2 supercell, other displacement

               20: MD on POSITIONS (2000  steps) lambda=0 (on LA )   without tox
               21: MD on POSITIONS (2001  steps) lambda=1 (on DFT)   without tox
               22: MD on POSITIONS (4000  steps) lambda=0 (on Hesse) without tox
               23: MD on POSITIONS (230   steps) lambda=1 (on DFT)   without tox
               24: MD on POSITIONS (8671  steps) 4x4x4sc DFT job 300K for michael check (WRONG ALAT 4.07)
               25: MD on POSITIONS (79012 steps) 4x4x4sc DFT job 300K 4.14 for michael check
               26: MD on POSITIONS (2001  steps) lambda=1 (900K on DFT) with tox (best std = 1.644)
                                                 -t1 -0.071 -D_mor 0.258 -a_mor 1.437
               27: MD on POSITIONS (2001  steps) lambda=1 (300K on DFT) with tox (best std = 1.644)
                                                 -t1 -0.073 -D_mor 0.300 -a_mor 1.430

               40: evolve MD on HesseMatrix 4x4x4 sc
               41: evolve MD on la 4x4x4 sc by using positions from Hessematrix 4x4x4
               42: evolve MD on la 4x4x4 sc by using positions from c skript
               43: evolve MD on 4x4x4 sc Cu@ 300K alat = 3.65762

               44: evolve MD on Al std best 900K
                   -t1 -0.072 -D_mor 0.258 -a_mor 1.437 -t 1863 -wa -wpf -steps 40000
               45: evolve MD on Al std best 300K
                   -t1 -0.072 -D_mor 0.258 -a_mor 1.437 -t 1863 -wa -wpf -steps 40000


               '''))


    p.add_argument('-wa',   '--write_analyze', required=False, action='store_true',
       help='write analysis files for investigations/correlations/etc. of forces', default=False)
    p.add_argument('-waf',   '--write_analyze_forces', required=False, action='store_true',
       help='add to out_dudl{,av}.dat the info of the forces std, out_forcesdiff.dat, out_forcesdiffav.dat, out_forces_vs_dft.dat', default=False)

    p.add_argument('-wp',   '--write_positions', required=False, action='store_true', help='write out_positions.dat', default=False)
    p.add_argument('-wpr',  '--write_positions_rel', required=False, action='store_true', default=False,
       help='write relative position coordinates instaed of absolute ones')
    p.add_argument('-wpf',   '--write_positions_forces', required=False, action='store_true', help='write out_positions_forces.dat', default=None)


    # morse parameters for 1NN
    p.add_argument('-alat_mor',   '--alat_mor', required=False,
       help='morse parameter for equilibrium lattice constant (default = False)', type=float, default=False)
    p.add_argument('-a_mor',   '--a_mor', required=False,
            help='morse parameter for widths of potential (default = False) typical for al is 1.432673', type=float, default=False)
    p.add_argument('-D_mor',   '--D_mor', required=False,
       help='morse parameter for well depth (default = False); value for al 0.27', type=float, default=False)
    p.add_argument('-t1',   '--t1', required=False,
       help='linear constant in transversal direction out of plane (default = 0.00)', type=float, default=False)
       #help='linear constant in transversal direction out of plane (default = -0.65)', type=float, default=-0.65)

    p.add_argument('-v','--verbose',
            help='verbose', action='count', default=False)
    p.add_argument('-vm1',   '--verbose_m1', required=False, action='store_true',
       help='print only python part to screen (not the part of the c-code)', default=False)
    p.add_argument('-vm2',   '--verbose_m2', required=False, action='store_true',
       help='print only summarizing python part to screen (not the part of the c-code)', default=False)
    p.add_argument('--stdout', required=False, action='store_true',
       help=argparse.SUPPRESS, default=False)
    p.add_argument('--dudl_lambda_0', required=False, action='store_true',
       help=argparse.SUPPRESS, default=False)
    p.add_argument('--dudl_lambda_1', required=False, action='store_true',
       help=argparse.SUPPRESS, default=False)
    p.add_argument('--calculate_energy_and_forces_exists', required=False, action='store_true',
       help=argparse.SUPPRESS, default=0)
    return p




def change_file(file_path, pattern, subst):
    #Create temp file
    fh, abs_path = mkstemp()
    with fdopen(fh,'w') as new_file:
        with open(file_path) as old_file:
            for line in old_file:
                new_file.write(re.sub(pattern,subst,line))
    #Remove original file
    #remove(file_path)

    #Move new file
    move(abs_path, file_path)
    return

def check_if_is_compilefile(compilefile,searchstring):
    out = []; file = open(compilefile, "r")
    for line in file:
         if re.search(searchstring, line): out.append(line)
    if len(out) != 1:
        if len(out) != 2:
            print('len(out)',len(out))
            sys.exit('following string is not in '+compilefile+": "+searchstring)
    return True

def ccc(file):
    from contextlib import closing
    from zipfile import ZipFile
    with closing(ZipFile(file)) as archive:
        count = len(archive.infolist())
    return

def check_for_disp_fit_parameters_morse(agrs):
    if type(args.folder) == bool:
        return False
    file = args.folder+"/disp_fit.parameters_morse.dat"
    if os.path.exists(file):
        a = np.loadtxt(file)
        #print(a)
        if args.D_mor == False or args.a_mor == False:
            args.D_mor = a[0]
            args.a_mor = a[1]
        if args.alat_mor == False:
            args.alat_mor = a[2]*np.sqrt(2)
        #print('D_mor',D_mor)
        #print('a_mor',a_mor)
        #print('alat_mor',alat_mor)
        #sys.exit()
    return

def obtain_parametrization_G_from_calculate_energy_and_forces(args):
    if type(args.folder) == bool:
        return False



    ############################# scheck for energy and forces ############################
    file = args.folder+"/calculate_energy_and_forces"
    if not os.path.exists(file):
        args.calculate_energy_and_forces = False
        if args.verbose:
            print("############ obtain_parametrization_G_from_calculate_energy_and_forces "+printred("(does NOT exist (1)")+ " ... #################")
        return False

    if args.verbose:
        print("############ obtain_parametrization_G_from_calculate_energy_and_forces ... #################################")


    ####### Here we assme now that calculate_energy_and_forces exists!
    args.calculate_energy_and_forces = file
    if args.verbose:
        print("############ check_for_calculate_energy_and_forces "+printgreen("(does exist)")+" ... ###########")
    #print printgreen('calculate_energy_and_forces does exist! taking parametrization from it!')
    #print printgreen(file)
    import imp
    par = imp.load_source('par', file)
    ver = 0
    try:
        par.u1nn_pot
        ver = 1
    except AttributeError:
        pass

    if ver == 0:
        try:
            par = par.parameters()
            ver = 2
            #print 'kk1',par.u1nn_pottype
            #print 'kk2',par.u1nn_pot
            par.u1nn_pot = par.u1nn_pottype
            if par.u1nn_pot == 'morse':
                par.u1nn_pot = 'm'
            #print 'kk3',par.u1nn_pot
        except:
            pass


    #print 'ver',ver
    #sys.exit()
    if ver == 0:
        print('dont understand this',file)
        return


    #print 'parx',par.u1nn_pot

    #print(printgreen(file))
    #print('par.u1nn_pot',par.u1nn_pot)
    #sys.exit()
    if par.u1nn_pot == 'm':
        if False:
            print('args.D_mor 1   ',args.D_mor)
            print('args.a_mor 1   ',args.a_mor)
            print('args.alat_mor 1',args.alat_mor)
        D_mor_tmp, a_mor_tmp, alat_mor_tmp = list(map(float,par.u1nn_potparam.split("_")))
        #print('got in',D_mor_tmp, a_mor_tmp, alat_mor_tmp)
        if type(args.D_mor) == bool:
            args.D_mor = D_mor_tmp
        if type(args.a_mor) == bool:
            args.a_mor = a_mor_tmp
        if type(args.alat_mor) == bool:
            args.alat_mor = alat_mor_tmp*np.sqrt(2.)  # back from NN dist to alat
        if False:
            print('args.D_mor 2   ',args.D_mor)
            print('args.a_mor 2   ',args.a_mor)
            print('args.alat_mor 2',args.alat_mor)
        #sys.exit('112233')
        #print par.u1nn_topx
        if ver == 1:
            if par.u1nn_topx != False:
                if False:
                    print('args.t1 1:',args.t1)
                args.t1 = par.u1nn_topx.split("_")[2]
                if False:
                    print('args.t1 2:',args.t1)
                #sys.exit()
            elif ver == 2:
                args.t1 = 0.
    else:
        print('some olf format')
        sys.exit("do I need to exit here or not?")

    if type(args.t1) == bool:
        args.t1 = 0.
    print_parameters(args,idx="(G)")
    return True

def get_supercell_C_at_least_try_from_eqcoords_or_poscar_if_exist(args):
    if args.verbose:
        print("############ get_supercell_C_at_least_try_from_eqcoords_or_poscar_if_exist ... ########################################")

    #########################################################
    ####################### get supercell size
    #########################################################
    #print 'args.sc 1:',args.supercell
    #print 'args.folder',args.folder
    if type(args.folder) != bool:
        eqcoords = args.folder+'/EqCoords_direct'
        if args.verbose:
            print(eqcoords+" do not exist")

        if os.path.isfile(eqcoords):
            lines = getlines_file(args.folder+'/EqCoords_direct')
            #print('eqcoords lines',lines)
            for i in np.arange(100):
                #print(i,i*i*i*4)
                if i*i*i*4 == lines:
                    args.supercell = int(i)
                    break
            #print 'args.sc 2',args.supercell,type(args.supercell)
    if type(args.supercell) == bool and os.path.isfile(args.folder+'/POSCAR'):
        file = open(args.folder+'/POSCAR','r')
        for idx,i in enumerate(file):
            if idx > 4 and idx < 8:
                #print i.rstrip(),type(i)
                try:
                    atoms = int(i.rstrip())
                    for i in np.arange(100):
                        if i*i*i*4 == atoms:
                            args.supercell = int(i)
                            #print 'atoms',atoms,args.supercell
                            break
                except ValueError:
                    pass
    print_parameters(args,idx="(C)")
    return

def get_alat_mor(args):
    if args.verbose:
        print("############ get_alat_mor ... ########################################")

    ##### try to get sc from cell file
    if type(args.supercell) == bool and args.read_positions != False:
        for i in range(10):
            if str(i)+"x"+str(i)+"x"+str(i)+"sc" in args.read_positions:
                args.supercell = i
                break
    if type(args.supercell) == bool:
        sys.exit('args.supercell NOT found')

    ##### try to get alat
    if type(args.supercell) != bool and type(args.folder) != bool and os.path.isfile(args.folder+'/cell'):
        cell = np.loadtxt(args.folder+'/cell')
        #print cell
        args.alat = cell[0,0]/args.supercell

    if type(args.alat) == bool and args.read_positions != False:
        from itertools import islice
        with open(args.read_positions) as myfile:
                head = list(islice(myfile, 2))
        head0 = head[0].rstrip().split()
        head1 = head[1].rstrip().split()
        #print('head0',head0)
        #print('head1',head1)
        if float(head0[0]) == 0 and float(head0[1]) == 0 and float(head0[2]) == 0 and float(head1[0]) == 0 and float(head1[1]) == 0 and float(head1[2]) != 0:
            args.alat = float(head1[2])
        if type(args.alat) == bool:  # when displacement in xy direction
            if float(head0[0]) == float(head0[1]) and float(head0[2]) == 0 and float(head1[0]) == 0 and float(head1[1]) == 0 and float(head1[2]) != 0:
                args.alat = float(head1[2])
    if type(args.alat) == bool:
        print_parameters(args,idx="(EXIT)")
        sys.exit('args.alat NOT found in get_alat_mor XXXX' )
    #if type(args.alat_mor) == bool:
    #    #print 'args.alat:',args.alat
    #    args.alat_mor = np.copy(args.alat)  # not so sure if I want that
    return

def print_parameters(args,idx="(?)"):
    if args.verbose > 1:
        print(">>>")
        print("args.D_mor ................ "+idx+":",args.D_mor)
        print("args.a_mor ................ "+idx+":",args.a_mor)
        print("args.alat ................. "+idx+":",args.alat)
        print("args.alat_mor ............. "+idx+":",args.alat_mor)
        print("args.t1 ................... "+idx+":",args.t1)
        print("args.element_all .......... "+idx+":",args.element_all)
        print("args.supercell ............ "+idx+":",args.supercell)
        print("args.folder ............... "+idx+":",args.folder)
        print("args.write_positions_forces "+idx+":",args.write_positions_forces)
        print("args.read_positions ....... "+idx+":",args.read_positions)
        print("args.read_hesse ........... "+idx+":",args.read_hesse)
        print("args.read_positions0 ...... "+idx+":",args.read_positions0)
        print("args.read_uoutcar ......... "+idx+":",args.read_uoutcar)
        print("<<<")
        return

def obtain_parametrization_E_from_parametrize_displacements_skript(args): #KKK
    if args.verbose:
        print("############ obtain_parametrization_E_from_parametrize_displacements_skript ... #################################")

    if         args.gpf_dispfolder_5x5x5 \
            or args.gpf_dispfolder_setfolder \
            or args.gpf_dispfolder_5x5x5_weighted \
            or args.gpf_dispfolder_5x5x5_alat_mor_T0K \
            or args.gpf_dispfolder_5x5x5_alat_mor_current_alat \
            or args.gpf_dispfolder_alat_lattice_T0K \
            or args.gpf_dispfolder_shift_parametrization_to_alat \
            or args.gpf_dispfolder_correct_for_110_forces:

        sc = False
        if     args.gpf_dispfolder_5x5x5 \
            or args.gpf_dispfolder_5x5x5_weighted \
            or args.gpf_dispfolder_5x5x5_alat_mor_T0K \
            or args.gpf_dispfolder_5x5x5_alat_mor_current_alat: sc = 5

	if args.gpf_dispfolder_5x5x5_alat_mor_current_alat:
	    args.gpf_dispfolder_shift_parametrization_to_alat = args.alat

        ele = get_param_disp.get_all_disps(
                element=args.element,
                sc = sc,
                dofor = args.gpf_dispfolder_setfolder,
                only_return_parametrization_file = True,
                folder_alat_lattice_T0K = args.gpf_dispfolder_alat_lattice_T0K,
                shift_parametrization_to_alat = args.gpf_dispfolder_shift_parametrization_to_alat,
                correct_for_110_forces = args.gpf_dispfolder_correct_for_110_forces
                )

        paramfile = ele.parametrize_it()
        if args.gpf_dispfolder_5x5x5_weighted == True:
            paramfile = ele.dofor+'/disp_fit.parameters_morse_weighted.dat'
        #print('pf',paramfile)
        #sys.exit()
        param = np.loadtxt(paramfile)

        if len(args.element_all) == 1:
            if type(args.D_mor) == bool: args.D_mor  = param[0]
            if type(args.a_mor) == bool: args.a_mor  = param[1]
            if type(args.alat_mor) == bool: args.alat_mor  = param[2]*np.sqrt(2.)
        elif len(args.element_all) > 1:
            args.D_mor  = param[0]
            args.a_mor  = param[1]
            args.alat_mor = param[2]*np.sqrt(2.)
        #args.alat_mor = 4.04

        if args.gpf_dispfolder_5x5x5_alat_mor_T0K:
            args.alat_mor = my_atom.alatT0K[args.element]

        print_parameters(args,idx="(E)")
        return

def obtain_input_folder_A_for_f_md_2x2x2_30_f_d_5x5x5(args):
    if args.verbose:
        print("############ obtain_input_folder_A_for_f_md_2x2x2_30_f_d_5x5x5 ... #################################")
    if args.f_md_2x2x2_30 == False and args.f_d_5x5x5 == False:
        return False

    if args.f_md_2x2x2_30:
        args.read_positions = False
        args.read_positions0 = False
        args.read_hesse = False
        args.read_uoutcar = False

        def getpath_MD_30(kp="3x3x3",lam="1"):
            base = "/Users/glensk/Dropbox/Albert/Understanding_distributions/ti/"
            path0 = "/30__PTS_dosall_from_2x2x2sc__LON_quer_"
            path1 = "kp_all_morse_LONADD__NO_LON2_NONE_LON2ADD__NO_TOX_NONE_TOY_NONE_TOZ_NONE/lambda"+lam+".0/"
            folder = base + args.element + path0 + kp + path1
            #print('f1',folder)
            if not os.path.isdir(folder):
                kp='2x2x2'
                folder = base + args.element + path0 + kp + path1
                if not os.path.isdir(folder):
                    kp='4x4x4'
                    folder = base + args.element + path0 + kp + path1
                    if not os.path.isdir(folder):
                        sys.exit('folder not found xxx;Exit')
            return folder
        args.folder = getpath_MD_30(kp="3x3x3",lam=args.lambda_)

    print('aaxxx',args.f_d_5x5x5)
    if args.f_d_5x5x5:
        args.folder = "/Users/glensk/Dropbox/Albert/Understanding_distributions/displacements_/"+args.element+"/5x5x5sc_quer_2x2x2kp_vasp4/"
        return "/Users/glensk/Dropbox/Albert/Understanding_distributions/displacements_/"+args.element+"/5x5x5sc_quer_2x2x2kp_vasp4/"

    if args.verbose:
        print('ARGS.FOLDER',args.folder)
    if args.alat == False:
        args.alat = np.loadtxt(args.folder+'/cell')[0,0]/2
    #args.read_hesse = False  # do I need this?
    print_parameters(args,idx="(A)")
    return

def obtain_input_folder_B_and_related_files(args):
    if args.verbose:
        print("############ obtain_input_folder_B_and_related_files ... #################################")

    if args.read_positions:
        if os.path.isfile(args.read_positions) != True:
            sys.exit("positionsfile "+args.read_positions+" does not exist (2)!")
        if os.path.isfile(args.read_positions) == True and not args.folder:
            print('here no args folder')
            args.folder = "/".join(args.read_positions.split("/")[:-1])

    if args.read_positions0:
        if os.path.isfile(args.read_positions0) != True:
            sys.exit("positionsfile for undisplaced structure "+args.read_positions0+" does not exist (3)!")

    if args.read_hesse:
        if os.path.isfile(args.read_hesse) != True:
            sys.exit("Hessematrix "+args.read_hesse+" does not exist (4)!")

    if args.read_uoutcar:
        if os.path.isfile(args.read_uoutcar) != True:
            sys.exit("u_OUTCAR "+args.read_uoutcar+" does not exist (5)!")
    if args.verbose:
        print("args.read_hesse xx",args.read_hesse)

    if args.folder:
        if os.path.isdir(args.folder) != True:
            sys.exit("Folder "+args.folder +" does not exist (6)!")
        else:  # the args.folder exists
            #print('the args.folder exists!')
            #sys.exit()
            if not args.read_positions and os.path.isfile(args.folder+'/POSITIONs'):
                args.read_positions = args.folder+'/POSITIONs'

            if not args.read_positions0 and os.path.isfile(args.folder+'/POSITIONs'):
                args.read_positions0 = args.folder+'/POSITIONs'

            if not args.read_hesse and os.path.isfile(args.folder+'/HesseMatrix_sphinx'):
                args.read_hesse = args.folder+'/HesseMatrix_sphinx'

            if not args.read_uoutcar and os.path.isfile(args.folder+'/u_OUTCAR_no_first_substract'):
                args.read_uoutcar = args.folder+'/u_OUTCAR_no_first_substract'

            if not args.read_uoutcar and os.path.isfile(args.folder+'/u_OUTCAR'):
                args.read_uoutcar = args.folder+'/u_OUTCAR'

            if not args.read_uoutcar:
                if os.path.isfile(args.folder+'/u_OUTCAR'):
                    args.read_uoutcar= args.folder+'/u_OUTCAR'
                elif os.path.isfile(args.folder+"/OUTCAR"):
                    import subprocess
                    int(subprocess.check_output('cd '+args.folder+';OUTCAR_ene-potential_energy.sh > u_OUTCAR', shell=True,stderr=subprocess.STDOUT))
                    #sys.exit('get u_OUTCAR 1@')
                    args.read_uoutcar= args.folder+'/u_OUTCAR'
                elif os.path.isfile(args.folder+"/OUTCAR.gz"):
                    print('args.folder',args.folder)


                    import subprocess
                    try:
                        int(subprocess.check_output('cd '+args.folder+';OUTCAR_ene-potential_energy.sh > u_OUTCAR', shell=True,stderr=subprocess.STDOUT))
                    except ValueError:
                        if os.path.isfile(args.folder+'/u_OUTCAR'):
                            args.read_uoutcar= args.folder+'/u_OUTCAR'
                        else:
                            sys.exit('u_OUTCAR was not created')
                    #sys.exit('get u_OUTCAR 2@')
                    args.read_uoutcar= args.folder+'/u_OUTCAR'
    if args.verbose:
        print("args.read_hesse yy",args.read_hesse)

    if args.read_hesse and not args.read_positions0:
        import myutils as my
        with my.cd(args.folder):
            call(["extractPOSITIONS.sh"])
        print('extracted POSITIONs')
        print('args.read_positions0:',args.read_positions0)

    if args.read_hesse and not args.read_positions0:
        sys.exit("when loading hesse you need the positions0 file -rp0 or --read_positions0 !")

    if args.verbose > 1:
        print("aa args.read_positions   :",args.read_positions)
        print("aa args.read_hesse       :",args.read_hesse)
        print("aa args.read_positions0  :",args.read_positions0)
        print("aa args.read_uoutcar     :",args.read_uoutcar)
        print()
    #sys.exit('33')


    if args.read_positions:
        #print 'args.read_forces 0',args.read_forces
        import subprocess
        columns = int(subprocess.check_output('head -1 '+args.read_positions+' | wc -w', shell=True,\
                stderr=subprocess.STDOUT))
        if columns == 3 and args.read_forces == True:
            print("POSITIONs file:",args.read_positions)
            sys.exit('args.read_forces is set to True but only 3 columns in POSITIONs file')
        elif columns == 3:
            args.read_forces = False
        elif columns == 6:
            args.read_forces = True
        else: sys.exit('Error: File '+args.read_positions+' has neither 3 nor 6 columns')
        #print 'args.read_forces',args.read_forces
        #print "columns:",columns
    print_parameters(args,idx="(B)")
    return



def pvi(var,typee=int,last=False):
    if type(var) == str:
        if last == True:
            return var.split("/")[-1],'1'
        else:
            return var,'1'
    #print 'ok'
    #print 'var',var
    #print 'typee',typee
    #print 'type(var)',type(var)
    #sys.exit()
    #print 'tv',typee(var)
    return var,str(typee(var))

def ov0(var,type=int):
    return pvi(var,type)[0]

def ov1(var,type=int):
    return pvi(var,type)[1]


def prepare_input_file(file_orig,file_comp,args):
    if args.verbose:
        print("############ prepare_input_file ... ###########################################")
    if not os.path.exists(file_orig): sys.exit("ERROR: could not find the c program "+file_orig+"!")
    copyfile(file_orig, file_comp)



    if args.verbose > 3:
        print("agrs:",args)


        #sys.exit('44444444444')
    if args.verbose:
        print('COMPILING (0) with args.supercell        :',args.supercell)
        print('COMPILING (0) with args.element_mass     :',args.element_mass)
        print('COMPILING (0) with args.alat             :',args.alat,type(args.alat))
        print('COMPILING (0) with args.alat_mor         :',args.alat_mor)
        print('COMPILING (0) with args.a_mor            :',args.a_mor)
        print('COMPILING (0) with args.D_mor            :',args.D_mor)
        print('COMPILING (0) with args.t1               :',args.t1)



    change_file(file_comp,'^#define N 2','#define N '+str(args.supercell))
    change_file(file_comp,'^#define mass_element.*','#define mass_element ('+str(args.element_mass)+")")
    change_file(file_comp,'^#define alat_lattice.*','#define alat_lattice ('+str(args.alat)+")")
    change_file(file_comp,'^#define alat_morse.*','#define alat_morse ('+str(args.alat_mor)+")")
    change_file(file_comp,'^#define a_mor.*','#define a_mor ('+str(args.a_mor)+")")
    change_file(file_comp,'^#define D_mor.*','#define D_mor ('+str(args.D_mor)+")")
    if type(args.t1) == bool:
        args.t1 = 0.0
    change_file(file_comp,'^#define ktr_tox.*','#define ktr_tox ('+str(args.t1)+")")

    if args.read_positions:
        change_file(file_comp,'^    const char \*filename_in_positions = .*','    const char *filename_in_positions = "'+str(args.read_positions)+'";')
	if args.verbose:
            print('COMPILING (1) with args.read_positions   :',args.read_positions)
    if args.read_positions0:
        change_file(file_comp,'^    int read_pos0=.*','    int read_pos0=1;')
        change_file(file_comp,'^    const char \*filename_in_positions0 = .*','    const char *filename_in_positions0 = "'+str(args.read_positions0)+'";')
	if args.verbose:
            print('COMPILING (1) with args.read_positions0  :',args.read_positions0)
    if args.write_analyze_forces:
        change_file(file_comp,'^    int write_analyze_forces=.*','    int write_analyze_forces=1;')
        np.savetxt("out_forces_vs_dft.dat.ref",np.array(([-5,-5],[5,5])))
    if args.read_hesse:
        change_file(file_comp,'^    const char \*filename_in_hesse = .*','    const char *filename_in_hesse = "'+str(args.read_hesse)+'";')
	if args.verbose:
            print('COMPILING (1) with args.read_hesse       :',args.read_hesse)
    if args.read_uoutcar:
        change_file(file_comp,'^    const char \*filename_in_uoutcar = .*','    const char *filename_in_uoutcar = "'+str(args.read_uoutcar)+'";')
	if args.verbose:
            print('COMPILING (1) with args.read_uoutcar     :',args.read_uoutcar)

    if type(args.supercell)     not in [int]: sys.exit('\nERROR: args.supercell is wrongly defined')
    if type(args.element_mass)  not in [np.float64,float]: sys.exit('\nERROR: args.element_mass is wrongly defined')
    #print('kka',type(args.alat_mor))
    if type(args.alat)          not in [np.float64,float]: sys.exit('\nERROR: args.alat is wrongly defined !!'+str(type(args.alat)))
    if type(args.alat_mor)      not in [np.float64,float]: sys.exit('\nERROR: args.alat_mor is wrongly defined')
    if type(args.a_mor)         not in [np.float64,float]: sys.exit('\nERROR: args.a_mor is wrongly defined')
    if type(args.D_mor)         not in [np.float64,float]: sys.exit('\nERROR: args.D_mor is wrongly defined')
    #print("t1",args.t1,type(args.t1))
    #print_parameters(args,idx="(tmpttt)")

    columns0 = False
    if args.read_positions0:
        import subprocess
        columns0 = int(subprocess.check_output('head -1 '+args.read_positions0+' | wc -w', shell=True,\
                stderr=subprocess.STDOUT))
        if columns0 != 3 and columns0 != 6:
            sys.exit('Error: File '+args.read_positions0+' has neither 3 nor 6 columns')
        #print "columns0:",columns0

    if columns0:
        change_file(file_comp,'^    int columns0=.*','    int columns0='+str(columns0)+";")
    change_file(file_comp,'^    int read_forces=.*','    int read_forces='+str(int(args.read_forces))+";")
    if args.evolve_harmonic:
        change_file(file_comp,'^    int evolve_md_on_hesse=.*','    int evolve_md_on_hesse=1;')

    if args.write_positions_forces and args.write_positions:
        args.write_positions = False
    if args.write_positions_forces:
        change_file(file_comp,'^    int write_positions_forces=.*','    int write_positions_forces=1;')
    if args.write_positions:
        change_file(file_comp,'^    int write_positions=.*','    int write_positions=1;')


    if args.read_positions:
        check_if_is_compilefile(file_comp,str(args.read_positions))
    if args.read_hesse:
        check_if_is_compilefile(file_comp,str(args.read_hesse))
    if args.read_uoutcar:
        check_if_is_compilefile(file_comp,str(args.read_uoutcar))



    last = False
    add = ""
    if args.folder:
        last = True
        add = "$FLOLDER/"
        if args.verbose:
            print("############ $FOLDER: ",args.folder)
    else:
        if args.verbose:
            print("############ $FOLDER: "+printred("None"))
    return

    def avi(add,var,last=False):
        #print 'add',add
        #print pvi(var,last=last)
        if pvi(var,last=last)[1] == '1':
            return add,pvi(var,last=last)[0],pvi(var,last=last)[1]
        else:
            return pvi(var,last=last)

    if args.verbose:
        print("COMPILING (2) with args.read_positions   :",avi(add,args.read_positions,last=last))
        print("COMPILING (2) with args.read_positions0  :",avi(add,args.read_positions0,last=last))
        print("COMPILING (2) with args.read_hesse       :",avi(add,args.read_hesse,last=last))
        print("COMPILING (2) with args.read_uoutcar     :",avi(add,args.read_uoutcar,last=last))
        print("COMPILING (2) with args.read_forces      :",pvi(args.read_forces))
        print("COMPILING (2) with columns0              :",columns0)
        print("COMPILING (2) with args.evolve_harmonic  :",pvi(args.evolve_harmonic))
        print()
        print("COMPILING (2) with args.write_positions  :",args.write_positions)
        print("COMPILING (2) with args.write_positions_forces:",args.write_positions_forces)
    sys.exit()
    return


def particular_args(args):
    if args.verbose:
        print("args.write_positions_forces 11:",args.write_positions_forces)
        print('sys.argv',sys.argv)
        print("args.test                     :",args.test                 )
        print("args.read_hesse               :",args.read_hesse           )
    if args.test and args.write_positions_forces != None:
        args.write_positions_forces = True
    else:
        return
    if args.verbose:
        print("args.write_positions_forces 12:",args.write_positions_forces)

    understanding = "/Users/glensk/Dropbox/Understanding_distributions/"
    disp          = understanding+"displacements_/"
    tial          = understanding+'ti/Al/'

    #######################################################3
    # run single structures (one snapshot)
    #######################################################3
    if args.test == 1:
        args.folder = understanding+"/ti/Al/30__PTS_dosall_from_2x2x2sc__LON_quer_3x3x3kp_all_morse_LONADD__NO_LON2_NONE_LON2ADD__NO_TOX_NONE_TOY_NONE_TOZ_NONE/lambda1.0/try_with_tox"
    if args.test == 2:
        args.folder = disp+"/Al/2x2x2sc_4.13Ang_xdir_3x3x3kp/4.13Ang_0.3/check_forces_python_without_tox"
    if args.test == 3:
        args.folder = disp+"/Al/2x2x2sc_4.13Ang_xdir_3x3x3kp/4.13Ang_0.3/check_forces_python_with_tox"
    if args.test == 4:
        args.folder = disp+"/Al/3x3x3sc_4.13Ang_xdir_2x2x2kp/4.13Ang_0.3/check_forces_python_without_tox"
    if args.test == 5:
        args.folder = disp+"/Al/3x3x3sc_4.13Ang_xdir_2x2x2kp/4.13Ang_0.3/check_forces_python_with_tox"
    if args.test == 6:
        args.folder = disp+"/Al/3x3x3sc_4.14Ang_xdir_10x10x10kp_vasp4_ENCUT400/4.14Ang_0.3/check_forces_python_with_tox"
    if args.test == 7:
        args.folder = disp+"/Al/2x2x2sc_4.13Ang_xdir_3x3x3kp/4.13Ang_0.5/check_forces_python_with_tox"

    #######################################################3
    # run MD from POSITIONS
    #######################################################3
    if args.test == 20:
        args.folder = tial+ "30__PTS_dosall_from_2x2x2sc__LON_quer_3x3x3kp_all_morse_LONADD__NO_LON2_NONE_LON2ADD__NO_TOX_NONE_TOY_NONE_TOZ_NONE/lambda0.0"
    if args.test == 21 or args.test == 26:
        args.folder = "/Users/glensk/Dropbox/Understanding_distributions/ti/Al/30__PTS_dosall_from_2x2x2sc__LON_quer_3x3x3kp_all_morse_LONADD__NO_LON2_NONE_LON2ADD__NO_TOX_NONE_TOY_NONE_TOZ_NONE/lambda1.0"
    if args.test == 22:
        args.folder = "/Users/glensk/Dropbox/Understanding_distributions/ti/ti_original_hesse/Al/lambda0.0_20604"
    if args.test == 23:
        args.folder = "/Users/glensk/Dropbox/Understanding_distributions/ti/ti_original_hesse/Al/lambda1.0_20620"

    #######################################################
    # no data other but POSITIONS
    #######################################################
    if args.test == 24:
        args.supercell = 4
        args.read_positions = "/Users/glensk/Dropbox/proj/proj_current/__2017.01_phonon_linewidth_al/__2017.01_phonon_lifetimes_4_ab2017/check_42_check_michaels_plot_different_displacements/michaels_POS/POSITIONs_300K_4.07_tats_4.14"
        args.alat = 4.14
        args.write_analyze = True

    if args.test == 25:
        args.supercell = 4
        args.folder = "/Users/glensk/Dropbox/proj/proj_current/__2017.01_phonon_linewidth_al/__2017.01_phonon_lifetimes_4_ab2017/check_42_check_michaels_plot_different_displacements/"
        args.alat = 4.07
        args.write_analyze = True
        args.steps = 5000

    if args.test == 26:
        args.folder = "/Users/glensk/Dropbox/Understanding_distributions/ti/Al/30__PTS_dosall_from_2x2x2sc__LON_quer_3x3x3kp_all_morse_LONADD__NO_LON2_NONE_LON2ADD__NO_TOX_NONE_TOY_NONE_TOZ_NONE/lambda1.0/"
        if type(args.t1) == bool:
            args.t1 = -0.071
        if type(args.D_mor) == bool:
            args.D_mor = 0.258
        if type(args.a_mor) == bool:
            args.a_mor = 1.437


    if args.test == 27:
        args.folder = "/Users/glensk/Dropbox/Albert/v/pp/Al/molecular_dynamics_lifetimes/low_4x4x4sc_250eV_2x2x2kp_EDIFF1E-2__4.07Ang_300K_GGA/SUM_run_1.save/"
        if type(args.t1) == bool:
            args.t1 = -0.073
        if type(args.D_mor) == bool:
            args.D_mor = 0.300
        if type(args.a_mor) == bool:
            args.a_mor = 1.430
        #if type(args.steps) == bool:
        # args.steps = 20000
        args.write_positions_forces = False



    #######################################################
    # evolve MD ... needs data about supercell, alat, ...
    #######################################################
    if args.test == 40:
        args.supercell = 4
        args.read_positions0 = "/Users/glensk/Dropbox/proj/proj_current/__2017.01_phonon_linewidth_al/__2017.01_phonon_lifetimes_4_ab2017/check_36_for_ankit_c_code/Hessematrix_4x4x4/POSITIONs"
        args.read_hesse = "/Users/glensk/Dropbox/proj/proj_current/__2017.01_phonon_linewidth_al/__2017.01_phonon_lifetimes_4_ab2017/check_36_for_ankit_c_code/Hessematrix_4x4x4/hessesym.dat"
        args.temperature = 1820
        args.evolve_harmonic = True

    if args.test == 41:
        args.supercell = 4
        args.temperature = 1814.6
        args.read_positions0 = "/Users/glensk/Dropbox/proj/proj_current/__2017.01_phonon_linewidth_al/__2017.01_phonon_lifetimes_4_ab2017/check_36_for_ankit_c_code/Hessematrix_4x4x4/POSITIONs"
        args.steps = 500000
        args.steps = 500000
        #args.steps = 1000-1
        args.write_every = 1
        #args.write_every = 10

    if args.test == 42:
        args.supercell = 4
        args.temperature = 1814.6
        args.steps = 500000



    if args.test == 43:
        args.supercell = 4
        args.alat = 3.65762
        args.temperature = 300
        args.dt = 1
        args.steps = 10
        args.element = "Cu"
        args.write_positions = True

    if args.test == 44:
        print("args.write_positions_forces 22:",args.write_positions_forces)
        args.temperature = 1852
        args.supercell = 3  # this is working up to 200000 nicely; probably much further
                            # this was however for a temperature of 851K (1852 starting temp)

        args.temperature = 1945
        args.supercell = 3  # this is working up to 200000 nicely; probably much further
                            # this was for a temperature of 900K (1945 starting temp)

        args.temperature = 1860
        args.supercell = 4  # this is working up to 200000; probably much futher
                            # this was in the end a temperauture of 906K

        args.temperature = 1849
        args.supercell = 4  # this is working up to 200000; probably much futher
                            # this was in the end a temperauture of 901K
        #args.steps = 300000 # a jump at 315581
        #args.steps = 300 # a jump at 315581

        ### THIS IS THE PARAMETRIZATION USED FOR THE LA paper !!!
        args.t1 = -0.072
        args.D_mor = 0.258
        args.a_mor = 1.437
        args.alat_mor = 4.13
        if type(args.steps) == bool:
            args.steps = 300000 # a jump at 315581
        args.write_analyze = False
        if args.write_positions_forces != None: args.write_positions_forces = False
        args.write_positions = False
        args.write_positions_forces = True
        print("args.write_positions_forces 33:",args.write_positions_forces)

    if args.test == 45:
        args.temperature = 606.8
        args.supercell = 4  # this is working up to 200000; probably much futher
        args.alat_mor = 4.07
        args.alat = 4.07
        args.t1 = -0.073
        args.D_mor = 0.300
        args.a_mor = 1.430
        args.steps = 800000
        args.write_positions_forces = False
        args.write_positions = True
    return

def printred(text):
    return ('\x1b[4;31;40m' + text + '\x1b[0m')
def printgreen(text):
    return ('\x1b[4;32;40m' + text + '\x1b[0m')

def compare_dudl(args):
    folder = args.folder
    if type(folder) == bool:
        return False
    if not os.path.isfile(folder+"/dUdL"):
        if int(args.steps) > 1:
            #print 'args.steps',args.steps,type(args.steps)
            print(printred('does not have old dUdL file'))
        return
    if not os.path.isfile("out_dudl.dat"):
        print(printred('does not have new out_dudl.dat file'))
        return
    #print "###############################################################################"
    #print "comparing dudl ..."
    #print "###############################################################################"
    dudl_old_whole = np.loadtxt(folder+"/dUdL")
    dudl_new_whole = np.loadtxt("out_dudl.dat")

    dudl_old = dudl_old_whole[:,6]
    dudl_new_la = dudl_new_whole[:,1][1:]
    dudl_new_harm = dudl_new_whole[:,2][1:]
    #print 'dudl_old',dudl_old
    #print 'dudl_new_la',dudl_new_la
    #print 'dudl_new_harm',dudl_new_harm
    dudl_diff = (abs(dudl_new_la - dudl_old)).max()
    dudl_diff_harm = (abs(dudl_new_harm - dudl_old)).max()
    if dudl_diff < 0.03:
        print(printgreen('DUDL   (whole file    ) OK dudl old and new (DFT-LA) are exactly similar, maxdiff '+str(dudl_diff)))
    else:
        if dudl_diff_harm < 0.03:
            print(printgreen('DUDL   (whole file    ) OK dudl old and new (DFT-H) are exactly similar, maxdiff '+str(dudl_diff_harm)))
        else:
            print(printred('dudl old and new (DFT-LA) differ! maxdiff '+str(dudl_diff)))
            print(printred('dudl old and new (DFT-H)  differ! maxdiff '+str(dudl_diff_harm)))
            print('dudl_old   :',folder+"/dUdL")
            print('dudl_new:',"out_dudl.dat")
            return

    if os.path.isfile("out_dudlav.dat"):
        av = np.loadtxt("out_dudlav.dat")
        #print 'av',av[-1]
        print(printgreen("Energy std <(DFT-LA)>       "+str(np.round(av[-1][1],2))))
        print(printgreen("Energy std <(DFT-H )>       "+str(np.round(av[-1][3],2))))

def read_in_results(args):
    # the onces which should always be there
    if os.path.isfile("out_dudlav.dat"):
        a = np.loadtxt("out_dudlav.dat")
    else:
        sys.exit('DONE, not out_dudlav.dat ... assuming MD')

    if len(a.shape) > 1:
        alast = a[-1]
    elif len(a.shape) == 1:
        alast = a
    else:
        sys.exit("len out_dudlav.dat is wired")

    #alast = a[-1]
    #print(a,a.shape)
    #print(len(a.shape),a[-1],'kk',a[False,1])
    #print('alast',alast)
    #if len(a.shape) == 1:


    # Energy std (DFT-LA) (lam 1.0)
    ene_std_dft_min_la          = alast[1]
    ene_dudl_dft_min_la         = alast[2]
    ene_std_dft_min_harmonic    = alast[3]
    ene_dudl_dft_min_harmonic   = alast[4]
    forces_std_la               = alast[10]
    forces_std_harmonic         = alast[11]
    units = "meV/atom"

    #print("ene_std_dft_min_la :",ene_std_dft_min_la)
    #print("ene_dudl_dft_min_la:",ene_dudl_dft_min_la)
    if type(args.calculate_energy_and_forces) != bool and os.path.exists(args.calculate_energy_and_forces):
        element = printgreen(args.element)
        args.calculate_energy_and_forces_exists += 1
    else:
        element = printred(args.element)
    if args.verbose_m1:
        print(element,"")
        print()
        print("Energy std (DFT-H)  :",ene_std_dft_min_harmonic,units)
        print("Energy std (DFT-LA) :",ene_std_dft_min_la,units)
        print()
        print("Energy dudl (DFT-H) :",ene_dudl_dft_min_harmonic,units)
        print("Energy dudl (DFT-LA):",ene_dudl_dft_min_la,units)
        print()
        print("Forces std (DFT-H):",forces_std_harmonic,"eV/angstrom")
        print("Forces std (DFT-LA):",forces_std_la,"eV/angstrom")
    if args.do_both_lambdas and args.lambda_ == "0":
        args.ene_std_lam0 = ene_std_dft_min_la
        args.for_std_lam0 = forces_std_la
        args.dudl_lambda_0 = ene_dudl_dft_min_la
    if args.lambda_ == "1":
        args.ene_std_lam1 = ene_std_dft_min_la
        args.for_std_lam1 = forces_std_la
        args.dudl_lambda_1 = ene_dudl_dft_min_la
    return



def compare_last_forces(folder):
    ''' to create the forces run simply
    /Users/glensk/Thermodynamics/python_thermodynamics/hesse.py al -eneext -vvv
    '''
    #print "###############################################################################"
    #print "comparing last forces ..."
    #print "###############################################################################"
    if type(folder) == bool:
        return False
    ffo = folder+'/forces'
    eeo = folder+'/energy'
    ffn = "out_positions_forces.dat"
    een = "out_dudl.dat"
    forces = False
    ok = True
    if not os.path.isfile(ffo):
        print ('does not have old forces file')
    if not os.path.isfile(eeo):
        print ('does not have old energy file')
    if os.path.isfile(ffo) and os.path.isfile(ffn):
        fo = np.loadtxt(ffo)
        print('lfo',len(fo))
        #print fo
        #print fo.shape,fo.shape[0]
        fna = np.loadtxt(ffn)
        print('lfn',len(fna))
        if len(fna) == 0:
            print(printred(ffn+' is empty'))
        else:
            #print 'fna',fna.shape
            fn = fna[-fo.shape[0]:][:,-3:]
            #print fn
            #print fn.shape
            fd = np.abs(fo-fn).max()
            if fd < 0.0001:
                print(printgreen("FORCES (last structure) OK"))
            else:
                ok=False
                print(printred("FORCES (last structure) PROBLEM"))
            print('force difference max:',np.round(fd,4))
    if os.path.isfile(eeo) and os.path.isfile(een):
        eo = np.loadtxt(eeo)
        en = np.loadtxt(een)
        #print 'leneo',len(eo)
        #print 'lenen',len(en)
        atoms = (args.supercell**3)*4
        eomev = eo*1000/(atoms-1)
        #print 'eomev',eomev
        enmev = 0
        #print 'en',en,en.shape,len(en.shape),en.shape[0]
        #if len(en.shape) == 1 and en.shape[0]==7:
        if en.shape[-1]==7:
            if len(en.shape) == 1:
                enmev = en[5]
            if len(en.shape) == 2:
                enmev = en[-1][5]
            #print 'enmev',enmev
        else:
            print('en',en)
            print('len(en.shape)',len(en.shape),en.shape)
        ed = np.abs(enmev-eomev).max()
        if ed < 0.01:
            print(printgreen("ENERGY (last structure) OK"))
        else:
            ok=False
            print(printred("ENERGY (last structure) PROBLEM"))
        print('energy difference max (meV):',np.round(ed,3))
        print('eomev',np.round(eomev,2))
        print('enmev',enmev)
    #jif ok:
    #j    print printgreen('to understand differences (to to corresponding old folder):')
    #j    print printgreen('/Users/glensk/Thermodynamics/python_thermodynamics/hesse.py al -eneext -vvv')
    #jelse:
    if ok == False:
        print(printred('to understand differences'))
        print(printred('/Users/glensk/Thermodynamics/python_thermodynamics/hesse.py al -eneext -vvv'))
        print('for displacement of 0.3 and foreconsant of ...')
        print('-0.066423494090*0.3^2 = -.005978114  and -.005978114 is the tox energy on one atom on ft2')
        print('the toxenergy contribution per atom is -0.005978*4*2*1000/31/2 == -.771354839')
        print('*1000/31 == faktor_enery_atom (*1/1.602e-19)')
        print('(*1/1.602e-19) == faktor force')
        print('-0.06642349409*0.3*2 == -0.039854 == force on one atom')
        print()
    return forces


def myexit(text):
    print()
    print("#####################################################")
    print(text)
    print("#####################################################")
    sys.exit()

def getlines_file(file):
    import subprocess
    lines = int(subprocess.check_output('wc -l '+file+' | awk \'{print $1}\'', shell=True,stderr=subprocess.STDOUT))
    return lines

def getsteps_POSITIONSfile(file,args):
    lines = getlines_file(file)
    #print 'lines',lines,type(lines)
    #print 'as',args.supercell,type(args.supercell)
    stepsq = lines/((args.supercell**3)*4)
    #print 'stepsq',stepsq,type(stepsq)
    #sys.exit()
    if type(stepsq) == int:
        return stepsq
    else:
        print('lines',lines)
        print('atoms',(args.supercell**3)*4)
        print('stepsq',stepsq)
        sys.exit('somehow lines by atoms is not an integer')

def common_start(*strings):
    """ Returns the longest common substring
        from the beginning of the `strings`
    """
    def _iter():
        for z in zip(*strings):
            if z.count(z[0]) == len(z):  # check all elements in `z` are the same
                yield z[0]
            else:
                return

    return ''.join(_iter())

if __name__ == '__main__':
    p = help()  # this gives the possibility to change some __init__ settings
    args = p.parse_args()
    script_path = os.path.dirname(os.path.realpath(__file__))
    print('script_path',script_path)
    file_orig = script_path+"/md_long_tox.c"
    file_comp = script_path+"/md_long_tox_compile.c"
    file_comp = os.getcwd()+"/md_long_tox_compile.c"
    if args.folder == ".": args.folder = os.getcwd()

    if type(args.folder_and_parametrization) != bool:
        if args.folder_and_parametrization == ".": args.folder_and_parametrization = os.getcwd()
        if type(args.folder) != bool: sys.exit('either args.folder or args.folder_and_parametrization')
        if type(args.gpf_dispfolder_setfolder) != bool: sys.exit('either args.gpf_dispfolder_setfolder or args.folder_and_parametrization')
        args.folder = args.folder_and_parametrization
        args.gpf_dispfolder_setfolder = args.folder_and_parametrization

    if args.verbose:
        print(">>>>>>>>>>>> file_orig:",file_orig,"file_comp:",file_comp)

    if args.test:
        particular_args(args)  # particular folder

    print_parameters(args,idx="(0)")

    if args.element_all == False:
        args.element_all = [args.element]
    if args.element_all == True:
        args.element_all = ["Al","Ag","Au","Cu","Ir","Pb","Pd","Pt","Rh"]

    ### make readme
    myutils.create_READMEtxt(os.getcwd())

    for args.element in args.element_all:
        args.calculate_energy_and_forces_exists = 0
        print_parameters(args,idx="(8)")
        if len(args.element_all) > 1:  # needs to be set, otherwise settings of old element are kept
            args.alat_mor = False
            args.a_mor = False
            args.D_mor = False
            args.t1 = False
            args.alat = False
            args.steps = False
            args.supercell = False

        print_parameters(args,idx="(9)")
        args.element_mass = my_atom.atom([args.element]).mass[0]

        lam01_ = [ "1" ]
        if args.do_both_lambdas:
            lam01_ = ["0","1"]

        args.ene_std_lam0 = "--"
        args.ene_std_lam1 = "--"
        args.for_std_lam0 = "--"
        args.for_std_lam1 = "--"
        args.dudl_lambda_0 = "--"
        args.dudl_lambda_1 = "--"

        def get_sweep_range(args):
            print("vvvvvvvvvvvvvvvvvv")
            D_mor_range = 0.003 # D ~ 0.2
            D_mor_range = 0.02 # D ~ 0.2
            D_mor_points = 5
            if args.sweep_D_mor:
                D_mor_sweep = np.linspace(args.D_mor - D_mor_range/2., args.D_mor + D_mor_range/2., num=D_mor_points)

            else:
                D_mor_sweep = [args.D_mor]
            print("D_mor_sweep (temporary)",D_mor_sweep)

            a_mor_range = 0.1 # a ~ 1.8
            a_mor_points = 10
            if args.sweep_a_mor:
                a_mor_sweep = np.linspace(args.a_mor - a_mor_range/2., args.a_mor + a_mor_range/2., num=a_mor_points)
            else:
                a_mor_sweep = [args.a_mor]
            print("a_mor_sweep (temporary)",a_mor_sweep)

            print()
            sweep = [[False,False,False]]  # first, this is a fake sweep, this will be filled once the
            print()
            print('sweep (1)',sweep)
            print()
            if len(D_mor_sweep) > 1 or len(a_mor_sweep) > 1:
                sweep = []
                for DD in D_mor_sweep:
                    for aa in a_mor_sweep:
                        print('DD, aa:',DD,aa)
                        sweep.append([DD,aa,False])
                print()
                print('sweep (2)',sweep)
                print()
                for ddaa in sweep:
                    print('iii loop:',ddaa)
            else:
                sweep = [[args.D_mor,args.a_mor,False]]
            print("^^^^^^^^^^^^^^^^^^")
            return sweep

        sweep_all = get_sweep_range(args)
        print('sweep_all (1)',sweep_all)



        for sweep_idx,sweep in enumerate(sweep_all):
            for args.lambda_ in lam01_:
                #########################################################################
                # get the parametrization
                #########################################################################
                obtain_input_folder_A_for_f_md_2x2x2_30_f_d_5x5x5(args)  # f_md_2x2x2_30
                obtain_input_folder_B_and_related_files(args)  # once args.folder is defined this is very good
                get_supercell_C_at_least_try_from_eqcoords_or_poscar_if_exist(args)
                obtain_parametrization_E_from_parametrize_displacements_skript(args)
                obtain_parametrization_G_from_calculate_energy_and_forces(args)
                #obtain_parametrization_I_from_
                # get parameters when (gpdisp)  can be either from disp_fit.parameters_morse.dat or from calculate_energy_and_forces.
                #if args.get_parametrization_from_displacements:
                #    parameters = Analyze_Forces_and_make_parametrization.get_all_disps(dofor=args.get_parametrization_from_displacements)
                #has_parametrization   = check_for_disp_fit_parameters_morse(args)
                get_alat_mor(args)
                print_parameters(args,idx="(11)")

                if sweep_idx == 0:
                    print('a!!',args.a_mor)
                    print('D!!',args.D_mor)
                    sweep_all = get_sweep_range(args)
                    print('sweep_all (2)',sweep_all)
                    print('sweep_all (3) x',sweep_all[0])
                    print('sweep_all (4) x',sweep_all[0][0])
                    print('sweep_all (5) y',sweep_all[0][1])

                args.D_mor = sweep_all[sweep_idx][0]
                args.a_mor = sweep_all[sweep_idx][1]

                ##############################################################################
                # now all variables should be defined
                ##############################################################################
                prepare_input_file(file_orig,file_comp,args)
                if args.verbose:
                    print("############ compiling (gcc)... ###############################################")
                if os.path.exists('a.out'):os.remove('a.out')
                if os.path.exists('out_dudl.dat'):os.remove('out_dudl.dat')
                call(["gcc", file_comp])  # on cmmc this shoudl be icc
                if not os.path.exists('a.out'):
                    call(["icc", file_comp])  # on cmmc this shoudl be icc
                if not os.path.exists('a.out'): sys.exit("no file a.out found! compilation problem")
                if args.verbose:
                    print("############ starting with following parameters ... ###########################")
                if args.read_positions and args.read_uoutcar and args.steps == 0:
                    args.steps = getlines_file(args.read_uoutcar)

                if args.read_positions and not args.read_uoutcar and args.steps == 0:
                    args.steps = getsteps_POSITIONSfile(args.read_positions,args)

                args.temperature = str(float(args.temperature))
                timestep = str(int(args.timestep)/1000.)


                if args.verbose:
                    print('START with args.temperature :',pvi(args.temperature,typee=float))
                    print('START with args.timestep    :',args.timestep,'(fs)')
                    print('START with timestep         :',timestep,'(ps)')
                    print('START with args.steps       :',pvi(args.steps))
                    print('START with write_every      :',str(args.write_every))
                    print('START with read_positions   :',ov1(args.read_positions))
                    print('START with verbose          :',ov1(args.verbose))
                    print('START with write_analyze    :',ov1(args.write_analyze))
                    print('START with read_hesse       :',ov1(args.read_hesse))
                    print('START with read_forces      :',ov1(args.read_forces))
                    print('START with read_uoutcar     :',ov1(args.read_uoutcar))
                    #sys.exit('kkkkabel')
                    print()
                if int(args.steps) == 0:
                    myexit("Error: Number of steps is 0")
                if args.verbose:
                    print("############ running md_long_tox_compile.c ... ##################################")
                    #sys.exit('111111111111111')
                    #                , T               ,dt           ,zeitschr  ,l                    ,read_pos, verbose, write_analyze , read_hesse
                    #call(  ["./a.out", args.temperature,timestep,args.steps,str(args.write_every),read_pos,verbose, write_analyze , read_hesse,read_forces,read_uoutcar])
                    print(printred('some input variables ar in COMPILING section and here ..... remove one'))
                    print(printred('forinstance write_analyze'))

                devnull = None   # print to screen
                if args.verbose_m1 or args.verbose_m2:
                    devnull = open(os.devnull, 'w')
                command =   ["./a.out",
                    ov0(args.temperature,float),
                    timestep,
                    ov1(args.steps),
                    str(args.write_every),
                    ov1(args.read_positions),
                    ov1(args.verbose),
                    ov1(args.write_analyze) ,
                    ov1(args.read_hesse),
                    ov1(args.read_forces),
                    ov1(args.read_uoutcar),
                    ov1(args.write_positions_rel)]
                if args.verbose:
                    print('command')
                    print(command)
                    print(" ".join(command))
                call(command,stdout=devnull)
                if args.verbose:
                    print("------- BACK in argparse_my_code.py -------")
                #print('args.stdout 1',args.stdout)
                read_in_results(args)
                #print('args.lambda_',args.lambda_)
                #print("??",args.element,"std:",args.stdout,"l0:",args.dudl_lambda_0,"l1:",args.dudl_lambda_1)
                #print('args.stdout 2',args.stdout)

                if args.test:
                    print("###############################################################################")
                    print("comparing result (energy, forces, dudl ...)")
                    print("###############################################################################")
                    print(args.folder)
                    forces_old = compare_last_forces(args.folder)
                    compare_dudl(args)
                    if args.calculate_energy_and_forces != bool:
                        print(printgreen("has calculate_energy_and_forces"))
                    else:
                        print(printred("does not have old calculate_energy_and_forces"))

                #print('args.dudl_lambda_0',args.dudl_lambda_0)
                #print('args.dudl_lambda_1',args.dudl_lambda_1)
                ### now the loop over lambdas is finished



            if args.do_both_lambdas == False:
	        free_ene_err = 0.
	        PRL_L = "(PRL -  )"
            else:
	        free_ene_err = np.abs((args.dudl_lambda_0 - args.dudl_lambda_1)/2.)
	        PRL_L = "(PRL NO )"
                if free_ene_err <= my_atom.PRL_L[args.element]+0.5:
                    PRL_L = "(PRL YES)"
            #print('args.dudl_lambda_0',args.dudl_lambda_0)
            #print('args.dudl_lambda_1',args.dudl_lambda_1)
            #sys.exit()
	    #import unicodedata
	    #print("I am " + unicodedata.lookup("GREEK SMALL LETTER PI"))
	    #print("I am " + unicodedata.lookup("GREEK SMALL LETTER DELTA"))
	    #DELTA=unicodedata.lookup("GREEK CAPITAL LETTER DELTA")
            print(("==> ("+args.element+") ene_std: {ene_std:5.1f} for_std: {for_std:6.5f}  dudl/2: {free_ene_err:5.1f}  ||| alat {alat:4.3f} ||| alat_mor {alat_mor:4.3f} a_mor {a_mor:7.5f}  D_mor {D_mor:7.5f} ||"+PRL_L).format(ene_std=args.ene_std_lam1, for_std=args.for_std_lam1, free_ene_err=free_ene_err,alat_mor=args.alat_mor,alat=args.alat,delta=1. - args.alat_mor/args.alat,a_mor=args.a_mor,D_mor=args.D_mor))
            sweep_all[sweep_idx][2] = args.for_std_lam1
            #if sweep_idx == 0 and len(sweep_all) > 1:
            if len(sweep_all) > 1:
               print(sweep_all[sweep_idx])

        if len(sweep_all) > 1:
            print()
            print(sweep_all)
            mat = np.asarray(sweep_all)
            np.savetxt('arrout.dat',mat)
            D = np.unique(mat[:,0])
            a = np.unique(mat[:,1])
            print("D",D)
            print("a",a)
            for idx,i in enumerate(D):
                np.savetxt('sweep_'+str(D[idx])+".dat",mat[np.where(mat[:,0] == D[idx])[0]][:,[1,2]])
