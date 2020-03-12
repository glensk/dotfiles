#!/bin/sh
ff=`ls -1d n2p2_v7ag_neinth20_*`
hier=`pwd`
for f in $ff;do
    cd $hier
    echo $f
    cd $f
        
        #getEnergies_byLammps.py -p . -sys fcc -sys_ele Al -evinet
		getEnergies_byLammps.py -p . --formation_energies beta2 --units eV -ea 
        getEnergies_byLammps.py -p . -sys fcc -sys_ele Al -sys_ncell 1 -ea  # needs to be redone

        #if [ ! -e "summary_formations_eV_peratom_beta2_DFT_T0.dat_per_formula_unit.dat" ]; then
        #    echo need to calculate
		#    getEnergies_byLammps.py -p . --formation_energies beta2 --units eV -ea 
        #fi

		#mkdir interstitials
		#cd interstitials
        #getEnergies_byLammps.py -p .. -i ../../aiida_get_structures_new/aiida_exported_group_Al6xxxDB_structures_calc__all_steps.input.data --pick_atoms_al 33 --units meV_pa --interstitial_form -v 
        #

    cd $hier
done


