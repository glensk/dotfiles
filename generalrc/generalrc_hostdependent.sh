#!/bin/sh
#################################################################
# this script returns host dependent environment variables 
# it does know all prior defined environment variables           
# it will not know the here defined environment variables 
# all here defined environment variables are defined by setenv  
# which will (in the end) work for bash / tcsh / zsh                         
# now using export instead of setenv
#################################################################


##############################################################################
# defaults / fall back options
##############################################################################
export mylaptop="mac";
export GZIP="-9";   # an empty string ("") is not possible, better use dummy
export myshell="zsh"; # set this specifically for other cases
export myhome_mac="/Users/glensk";
export myhome_daint="/users/aglensk";
export myhome_fidis="/home/glensk";
export myhome_cosmopc="/home/glensk";

#case $onhost in
case $myhost in
mac)
    export myprompthostuser="magenta";
    export mypromptpath="magenta";
    export lmp_execi="$scripts/executables/lmp_mac_serial_runner_2018_10_30";
    # conda it seems can not be loaded here
  ;;
cosmopc)
    export myprompthostuser="green";
    export mypromptpath="green";
	export SCRATCH="/local/scratch/glensk";
	complete -d cd;
  ;;
fidis)
    export myprompthostuser="cyan";
    export mypromptpath="cyan";
	export SCRATCH="/scratch/glensk";
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
    export myprompthostuser="blue";
    export mypromptpath="blue";
	#export SCRATCH /scratch/glensk;"  # DONT SET SCRATCH SINCE THIS IS DONE GLOBALY BY SYSADMINS
    [ "$gettime" = "true" ] && gett=`gt $gett` && echo "general (1) : $gett time before moduel load"
	#module load daint-mc && module switch PrgEnv-cray PrgEnv-intel && module unload cray-libsci && module load GSL/2.5-CrayIntel-18.08 cray-python/2.7.15.1 cray-fftw
    [ "$gettime" = "true" ] && gett=`gt $gett` && echo "general (1) : $gett time after moduel load"
  ;;
cmpc)
    export myprompthostuser="red";
    export mypromptpath="red";
   	module load mathematica1002/10.0.2;
   	limit stacksize unlimited;
   	export TMPDIR="/scratch/$USER";
   	export PRINTER="cmcopy1166";
   	# christoph changed in /etc/profile.d/gxhivemgr.sh: if test x"$PS1" = x || echo
  ;;
cmmc)
    export myprompthostuser="cyan";
    export mypromptpath="cyan";
	module load git intel impi mkl grace mathematica/10.3.1 lammps vasp/5.3.5 anaconda;
	#add="limit stacksize unlimited;"
    #add2="setenv vaspq /u/aglen/vasp/vasp_4.6_lj_morse_alles_copied_FFTWPLANS_works/vasp.v1;"
	echo "export vaspq=\"/u/aglen/vasp/vasp_4.6_lj_morse_alles_copied_FFTWPLANS_works/vasp.v1\";"
  ;;
helvetios)
    export myprompthostuser="yellow";
    export mypromptpath="yellow";
	export SCRATCH="/scratch/glensk";
    module load git; # intel python/2.7.14;  # to get git 2.17 (instad of 1.8)
  ;;
*)
	: # export onhost="UNKNOWN";  : is the bash equivalent for pass in python
  ;;
esac

tab-color $myprompthostuser
#export whichalias mkalias;
#[ "$currentshell" = "tcsh" ] && setenv whichalias alias;  # I dont need to support tcsh anymore
