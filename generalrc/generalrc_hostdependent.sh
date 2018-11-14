#!/bin/bash
#################################################################
# this script returns host dependent environment variables 
# it does know all prior defined environment variables           
# it will not know the here defined environment variables 
# all here defined environment variables are defined by setenv  
# which will (in the end) work for bash / tcsh / zsh                         
#################################################################

#host=`hostname`   # 0.001s
#onhost=`echo $host | sed 's/\..*$//'`
#echo "setenv host $host;
#      setenv onhost $onhost;
#      setenv mylaptop mac;"
echo "setenv mylaptop mac;"
# cmpc25=wolfgang; cmpc07=biswanath;

# Prompt colors: use 
# time:             zsh: magenta;    bash: orange;      tcsh: green
# host& path:       mac: magenta;    cmpc: red;       cmmd: blue;    cmmc:turquoise

##############################################################################
# defaults / fall back options
##############################################################################
echo "setenv GZIP -9;"   # an empty string ("") is not possible, better use dummy
echo "setenv myshell zsh;" # set this specifically for other cases
#echo "setenv ipi_mc $scripts/i-pi-mc/bin/i-pi;"

case $onhost in
mac)
    echo "setenv myprompttime magenta;"
    echo "setenv myprompthostuser magenta;"
    echo "setenv mypromptpath magenta;"
    echo "setenv lmp_exec $scripts/lammps_executables/lmp_serial_mac_runner_2018_10_30;"
  ;;
cosmopc)
    echo "setenv myprompttime green ;"
    echo "setenv myprompthostuser green;"
    echo "setenv mypromptpath green;"
	echo "setenv SCRATCH /local/scratch/glensk;"
	echo "complete -d cd;"
    echo "setenv lmp_exec $scripts/lammps_executables/lmp_serial_cosmopc18_runner_2018_10_30;"
  ;;
fidis)
    echo "setenv myprompttime cyan;"
    echo "setenv myprompthostuser cyan;"
    echo "setenv mypromptpath cyan;"
	echo "setenv SCRATCH /scratch/glensk;"
    echo "setenv lmp_exec $scripts/lammps_executables/lmp_fidis_fidis_runner_2018_10_31;"
  ;;
daint)
    echo "setenv myprompttime blue;"
    echo "setenv myprompthostuser blue;"
    echo "setenv mypromptpath blue;"
	#echo "setenv SCRATCH /scratch/glensk;"  # DONT SET SCRATCH SINCE THIS IS DONE GLOBALY BY SYSADMINS
  ;;
cmpc)
    echo "setenv myprompttime red;"
    echo "setenv myprompthostuser red;"
    echo "setenv mypromptpath red;"
   	echo "module load mathematica1002/10.0.2;"
   	# christoph changed in /etc/profile.d/gxhivemgr.sh: if test x"$PS1" = x || echo
   	echo "limit stacksize unlimited;"
   	echo "setenv TMPDIR /scratch/$USER;"
   	echo "setenv PRINTER cmcopy1166;"
  ;;
cmmc)
    echo "setenv myprompttime cyan;"
    echo "setenv myprompthostuser cyan;"
    echo "setenv mypromptpath cyan;"
	echo "module load git intel impi mkl grace mathematica/10.3.1 lammps vasp/5.3.5 anaconda;"
	#add="limit stacksize unlimited;"
    #add2="setenv vaspq /u/aglen/vasp/vasp_4.6_lj_morse_alles_copied_FFTWPLANS_works/vasp.v1;"
	echo "setenv vaspq /u/aglen/vasp/vasp_4.6_lj_morse_alles_copied_FFTWPLANS_works/vasp.v1;"
  ;;
*)
	echo "setenv onhost UNKNOWN;"
  ;;
esac

echo "setenv whichalias mkalias;"
[ "$currentshell" = "tcsh" ] && echo "setenv whichalias alias;"
