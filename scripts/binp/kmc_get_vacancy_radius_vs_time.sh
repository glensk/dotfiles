#!/bin/sh
echo "dotfiles:$dotfiles:"
[ ! -e "vacancy_pos_angstrom.xyz" ] && echo "vacancy_pos_angstrom.xyz does not exist" && exit
[ ! -e "KMC_AL6XXX" ] && echo "KMC_AL6XXX does not exist" && exit
[ ! -e "simulation.pos_0.xyz" ] && echo "simulation.pos_0.xyz does not exist" && exit
[ ! -e "$dotfiles/scripts/i-pi-mc_scripts/vmd_vac_script.tcl" ] && echo "$dotfiles/scripts/i-pi-mc_scripts/vmd_vac_script.tcl does not exist" && exit
[ -e "vacancy_radius_squared.dat" ] && rm vacancy_radius_squared.dat

pbc=`head -2 simulation.pos_0.xyz | tail -1 | awk '{print $3,$4,$5,$6,$7,$8}'`
cp $dotfiles/scripts/i-pi-mc_scripts/vmd_vac_script.tcl .

sed -i 's|mol load xyz.*|mol load xyz "'"`pwd`/vacancy_pos_angstrom.xyz"'"|' vmd_vac_script.tcl 
sed -i 's|pbc set {.*|pbc set {'"$pbc"'} -all|' vmd_vac_script.tcl 
vmd -dispdev text -e vmd_vac_script.tcl
paste vacancy_radius_squared.dat KMC_AL6XXX | awk '{print $2*2.4e-17,$1}' > vacancy_radius_squared_vs_time.dat
[ -e "vacancy_radius_squared.dat" ] && rm vacancy_radius_squared.dat
