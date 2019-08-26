#!/usr/bin/env python
cd /Users/glensk/Dropbox/Thermodynamics/python_thermodynamics
cp phonon_lifetimes.py phonon_lifetimes_tmp2.py
sed -i 's/\([ ]*print \)\(.*\)\(#\)/\1(\2)/p' phonon_lifetimes_tmp2.py


cd /Users/glensk/Dropbox/Albert/v/pp/Al/molecular_dynamics_lifetimes/low_4x4x4sc_250eV_2x2x2kp_EDIFF1E-2__4.07Ang_300K_GGA_2x2x2/SUM_run_1.save
%~/Thermodynamics/python_thermodynamics/phonon_lifetimes_tmp2.py

