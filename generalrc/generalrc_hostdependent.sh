#!/bin/bash
#################################################################
# this script returns host dependent environment variables 
# it does know all prior defined environment variables           
# it will not know the here defined environment variables 
# all here defined environment variables are defined by setenv  
# which will (in the end) work for bash / tcsh / zsh                         
#################################################################

#################################################################
# determine the current host;
#   on{mac,cmmc,cmmc,cmpc} == true/false
#   myshell == zsh / tcsh / bash
#   whichalias == alias / mkalias  (alias for tcsh) (mkalias for bash/zsh)
#   currenthost == on{mac,cmmc,cmmc,cmpc}
#################################################################

host=`hostname`   # 0.001s
onhost=`echo $host | sed 's/\..*$//'`
echo "setenv host $host;
      setenv onhost $onhost;
      setenv myhost cmpc34;
      setenv myhostmpie cmpc34.mpie.de;
      setenv mylaptop mac;
      setenv submithost cmmc001;
      setenv printloadstat false;"
# cmpc25=wolfgang; cmpc07=biswanath;

# Prompt colors: use 
# time:             zsh: magenta;    bash: orange;      tcsh: green
# host& path:       mac: magenta;    cmpc: red;       cmmd: blue;    cmmc:turquoise

##############################################################################
# fall back options
##############################################################################
onxxx="setenv ka kb"
currenthost="setenv currenthost unknown"
myshell="setenv myshell bash"
myprompttime="setenv myprompttime magenta"
myprompthostuser="setenv myprompthostuser red"
mypromptpath="setenv mypromptpath red"
module="setenv GZIP -9"   # an empty string ("") is not possible, better use dummy
add="setenv ka kb"      # an empty string ("") is not possible, better use dummy

case $onhost in
mac)
    echo "setenv onmac true;"
    echo "setenv currenthost onmac;"
    echo "setenv myshell zsh;"
    echo "setenv myprompttime magenta;"
    echo "setenv myprompthostuser magenta;"
    echo "setenv mypromptpath magenta;"
    echo "setenv GZIP -9;"
  ;;
cosmopc)
	echo "setenv myshell zsh;"
	echo "setenv SCRATCH /local/scratch/glensk;"

  ;;
cmpc)
   	onxxx="setenv oncmpc true"
   	currenthost="setenv currenthost oncmpc"
   	myshell="setenv myshell zsh"
   	myprompttime="setenv myprompttime magenta"
   	myprompthostuser="setenv myprompthostuser red"
   	mypromptpath="setenv mypromptpath red"
   	module="module load mathematica1002/10.0.2"
   	# christoph changed in /etc/profile.d/gxhivemgr.sh: if test x"$PS1" = x || echo
   	add="limit stacksize unlimited;setenv TMPDIR /scratch/$USER;setenv PRINTER cmcopy1166"
  ;;
cmmc)
	onxxx="setenv oncmmc true"
	currenthost="setenv currenthost oncmmc"
	myshell="setenv myshell zsh"
	myprompttime="setenv myprompttime cyan"
	myprompthostuser="setenv myprompthostuser cyan"
	mypromptpath="setenv mypromptpath cyan"
	module="module load git intel impi mkl grace mathematica/10.3.1 lammps vasp/5.3.5 anaconda"
	#add="limit stacksize unlimited;"
    #add2="setenv vaspq /u/aglen/vasp/vasp_4.6_lj_morse_alles_copied_FFTWPLANS_works/vasp.v1;"
    #add3="setenv submithost cmmc002;"
	add="setenv vaspq /u/aglen/vasp/vasp_4.6_lj_morse_alles_copied_FFTWPLANS_works/vasp.v1"
	add2="setenv submithost cmmc002;"
  ;;
daint)
	echo "setenv myshell zsh;"
	#echo "setenv SCRATCH /scratch/glensk;"  # DONT SET SCRATCH SINCE THIS IS DONE GLOBALY BY SYSADMINS
  	currenthost="setenv currenthost ondaint"
  	myprompttime="setenv myprompttime blue"
  	myprompthostuser="setenv myprompthostuser blue"
  	mypromptpath="setenv mypromptpath blue;"
  ;;
fidis)
	echo "setenv SCRATCH /scratch/glensk;"
  ;;
*)
	echo "setenv HOSTUNKNOWN TRUE"
  ;;
esac

echo "setenv whichalias mkalias;"
[ "$currentshell" = "tcsh" ] && whichalias="setenv whichalias alias"
