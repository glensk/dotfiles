#!/bin/sh
echo "dotfiles:$dotfiles:"
[ ! -e "vacancy_pos.xyz" ] && echo "vacancy_pos.xyz does not exist" && exit
[ ! -e "KMC_AL6XXX" ] && echo "KMC_AL6XXX does not exist" && exit
[ ! -e "simulation.pos_0.xyz" ] && echo "simulation.pos_0.xyz does not exist" && exit
[ ! -e "$dotfiles/scripts/i-pi-mc_scripts/vmd_vac_script.tcl" ] && echo "$dotfiles/scripts/i-pi-mc_scripts/vmd_vac_script.tcl does not exist" && exit
[ ! -e "vacancy_radius.dat" ] && rm vacancy_radius.dat

pbc=`head -2 simulation.pos_0.xyz | tail -1 | awk '{print $3,$4,$5,$6,$7,$8}'`
cp $dotfiles/scripts/i-pi-mc_scripts/vmd_vac_script.tcl .

sed -i 's|mol load xyz.*|mol load xyz "'"`pwd`/vacancy_pos.xyz"'"|' vmd_vac_script.tcl 
sed -i 's|pbc set {.*|pbc set {'"$pbc"'} -all|' vmd_vac_script.tcl 
vmd -dispdev text -e vmd_vac_script.tcl
paste vacancy_radius.dat KMC_AL6XXX | awk '{print $2,$1}' > vacancy_radius_vs_time.dat
