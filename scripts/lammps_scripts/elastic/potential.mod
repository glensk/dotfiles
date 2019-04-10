# NOTE: This script can be modified for different pair styles 
# See in.elastic for more info.

# Choose potential
#pair_style	sw
#pair_coeff * * Si.sw Si
#mass 1 24.305
#mass 2 26.9815385
#mass 3 28.0855

# Runner pot
variable runnerDir string "rv64_1"
variable runnerDir string "/Users/glensk/Dropbox/Albert/scripts/dotfiles/scripts/potentials/runner_v2dg"
pair_style runner dir ${runnerDir} showewsum 1 showew yes resetew no maxew 1000000
pair_coeff * * 7.937658735


# NN POT n2p2
#variable nnpDir string "/home/glensk/Dropbox/Albert/scripts/dotfiles/scripts/potentials/n2p2_v1ag"
#variable nnpDir string "/home/glensk/Dropbox/Albert/scripts/dotfiles/scripts/potentials/n2p2_v2ag"
#pair_style nnp dir ${nnpDir} showew no resetew yes maxew 100000000 cflength 1.8897261328 cfenergy 0.0367493254
#pair_coeff * * 11.0

#pair_style eam/alloy
#pair_coeff * * Al.eam.alloy_cutoff6_seed_127519_std_0.715638 Al
#pair_coeff * * Al99.eam.alloy Al

#pair_style meam/alloy
#pair_style meam
#pair_coeff * * library.meam Al meam.alsimgcufe Al

# Setup neighbor style
neighbor 1.0 nsq
neigh_modify once no every 1 delay 0 check yes

# Setup minimization style
min_style	     cg
min_modify	     dmax ${dmax} line quadratic

# Setup output
thermo		1
thermo_style custom step temp pe press pxx pyy pzz pxy pxz pyz lx ly lz vol
thermo_modify norm no
