#!/bin/sh
#################################################################
# this script returns host dependent environment variables 
# it does know all prior defined environment variables           
# it will not know the here defined environment variables 
# all here defined environment variables are defined by setenv  
# which will (in the end) work for bash / tcsh / zsh                         
#################################################################


##############################################################################
# defaults / fall back options
##############################################################################
setenv mylaptop mac;
setenv GZIP -9;   # an empty string ("") is not possible, better use dummy
setenv myshell zsh; # set this specifically for other cases
setenv myhome_mac /Users/glensk;
setenv myhome_daint /users/aglensk;
setenv myhome_fidis /home/glensk;
setenv myhome_cosmopc /home/glensk;

case $onhost in
mac)
    setenv myprompttime magenta;
    setenv myprompthostuser magenta;
    setenv mypromptpath magenta;
    setenv lmp_exec $scripts/executables/lmp_mac_serial_runner_2018_10_30;
    # conda it seems can not be loaded here
  ;;
cosmopc)
    setenv myprompttime green ;
    setenv myprompthostuser green;
    setenv mypromptpath green;
	setenv SCRATCH /local/scratch/glensk;
	complete -d cd;
  ;;
fidis)
    setenv myprompttime cyan;
    setenv myprompthostuser cyan;
    setenv mypromptpath cyan;
	setenv SCRATCH /scratch/glensk;
    module load git;  # to get git 2.17 (instad of 1.8)
    #module load intel;
    # all of those take too long and are not necessary in interactive work
    #module load intel-mpi;
    #module load intel-mkl;
    #module load gsl;
    #module load eigen;
    #module load quantum-espresso;
    #module load python/2.7.14;  # conflicts with miniconda
    
    ################################# 
    # those might be the most important ones
    # but maybe can be loaded from manually if python necessary
    ################################# 
    #module load intel;          # intel is necessary to load python
    #module load python/2.7.14;  # conflicts with miniconda
  ;;
daint)
    setenv myprompttime blue;
    setenv myprompthostuser blue;
    setenv mypromptpath blue;
	#setenv SCRATCH /scratch/glensk;"  # DONT SET SCRATCH SINCE THIS IS DONE GLOBALY BY SYSADMINS
  ;;
cmpc)
    setenv myprompttime red;
    setenv myprompthostuser red;
    setenv mypromptpath red;
   	module load mathematica1002/10.0.2;
   	limit stacksize unlimited;
   	setenv TMPDIR /scratch/$USER;
   	setenv PRINTER cmcopy1166;
   	# christoph changed in /etc/profile.d/gxhivemgr.sh: if test x"$PS1" = x || echo
  ;;
cmmc)
    setenv myprompttime cyan;
    setenv myprompthostuser cyan;
    setenv mypromptpath cyan;
	module load git intel impi mkl grace mathematica/10.3.1 lammps vasp/5.3.5 anaconda;
	#add="limit stacksize unlimited;"
    #add2="setenv vaspq /u/aglen/vasp/vasp_4.6_lj_morse_alles_copied_FFTWPLANS_works/vasp.v1;"
	echo "setenv vaspq /u/aglen/vasp/vasp_4.6_lj_morse_alles_copied_FFTWPLANS_works/vasp.v1;"
  ;;
helvetios)
    setenv myprompttime yellow;
    setenv myprompthostuser yellow;
    setenv mypromptpath yellow;
	setenv SCRATCH /scratch/glensk;
    module load git; # intel python/2.7.14;  # to get git 2.17 (instad of 1.8)
  ;;
*)
	setenv onhost UNKNOWN;
  ;;
esac

setenv whichalias mkalias;
#[ "$currentshell" = "tcsh" ] && setenv whichalias alias;  # I dont need to support tcsh anymore
