#!/bin/sh

echo 'execute this using: phonon_lifeties.py > allout.txt 2>&1 '
cd /cmmc/ptmp/aglen/Understand_phonon_lifetimes/check_3_sc_convergence
~/Thermodynamics/python_thermodynamics/phonon_lifetimes.py -ps -q t2 t2 0 -f "nv*/sc*"  -quiet

cd /cmmc/ptmp/aglen/Understand_phonon_lifetimes/check_3_sc_convergence_300K
~/Thermodynamics/python_thermodynamics/phonon_lifetimes.py -ps -q t2 t2 0 -f "nv*/sc*" -quiet

cd /cmmc/ptmp/aglen/Understand_phonon_lifetimes/check_4_temperatuer_effects
~/Thermodynamics/python_thermodynamics/phonon_lifetimes.py -ps -q t2 t2 0 -f "*00K_*" -quiet

cd /cmmc/ptmp/aglen/Understand_phonon_lifetimes/check_4_temperatuer_effects_LA
~/Thermodynamics/python_thermodynamics/phonon_lifetimes.py -ps -q t2 t2 0 -f "*00K_*" -quiet

cd /cmmc/ptmp/aglen/Understand_phonon_lifetimes/check_6_LA_both_parametrizations_diff_cutoffs
~/Thermodynamics/python_thermodynamics/phonon_lifetimes.py -ps -q t2 t2 0 -f "nv*"  -quiet
