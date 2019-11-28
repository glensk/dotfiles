#!/bin/bash

#-----set parameters and paths------------------------
Tmax=1000; TstartDef=1; TstepDef=1;
kB=11.6045059554883246
#-----------------------------------------------------


# following 3 lines must always be present
path=`set | grep BASH_SOURCE | sed 's|.*\"\(.*\)/[^/]*\".*|\1|' | sed '/BASH_SOURCE=/s/.*/\./'`
script=`set | grep BASH_SOURCE | sed 's|.*\".*/\([^/]*\)\".*|\1|' | sed '/BASH_SOURCE=/s/.*/\./'`
options=$*; . $path/utilities/functions.include; checkOptions "-h -help -c -e -Tstart -Tstep -Tmax"

# small help for options
# space between option and value must be 1
# space between option and Comment or value and Comment must be >2
h=`getOption -h`
if [ $h == True ]; then
  usage $script
  printNeeded  "-c      concs     arbitrary many concentrations, must add up to one"
  printOptions "-e      elem      element(s) to determine Tmelt and take it as Tmax"     \
               "-Tstart Tstart    start T for Fqh mesh (def: Tstart=$TstartDef)"         \
               "-Tstep  Tstep     T step for Fqh mesh (def: Tstep=$TstepDef)"            \
               "-Tmax   Tmax      max T for Fqh mesh (def: Tmax=$Tmax or melting T)"
  exit
fi

# detailed help
help=`getOption -help`
if [ $help == True ]; then
  details $script
  echo2 "   calculate configurational entropy and free energy for the given" \
        "   concentrations as a function of temperatures"
  echo2 "   if element(s) is(are) given then the max temperature is taken to be Tmelt" \
        "   which for more elements is linearly mixed from the consitutents with the concentrations"
  echo2 "   typical script sequence (for emto):"                                                     \
        "    createFolders_murn_emto.sh -a           (adjust parameters.dat)"                        \
        "    createFolders_murn_emto.sh"                                                             \
        "    submit.sh run.emto.ser.cluster.workdir  (or the local version)"                         \
        "    collectEMTO.sh                          (when finished and inside NxNxNkp_XXeV folder)" \
        "    extractEnergies_emto.sh"                                                                \
        "    getE0KFit.sh -i energies_Tto0K_GGA      (cmpc01)"                                       \
        "    getDebye.sh -e elem/mass -c conc        (creates Debye Fqhs)"                           \
        "    getFqhFit.sh                            (cmpc01)"                                       \
        "    \033[1mgetConfEntropy.sh\033[0m -c conc -e elem           (creates Fconf)"                 \
        "    getThermodynamics.sh                    (collect EVinet,Fqh_fit,Fconf in separate folder)" \
        "                                            (to add Fel check createFolders_elec_emto.sh -help)"
  exit
fi


# we need -c option with concentrations for more than one element
conc=`getOption -c`
if [ $conc != True ]; then error "-c option with concentrations required"; fi
conc=`getValue -c`
nconc=`echo $conc | awk '{print NF}'`
# check if concentrations add up to one
totConc=`echo $conc | awk '{c=0; for (i=1;i<=NF;i++) c=c+$i; if ((c-1)^2>0.0001) print "wrong"}'`
if [ "$totConc" == wrong ]; then error "concentrations in -c option do not add up to one"; fi

el=`getOption -e`
tmaxOp=`getOption -Tmax`
if [ $el == True -a $tmaxOp == True ]; then error "-e option cannot be used together with -Tmax option"; fi
if [ $el == True ]; then
  el=`getValue -e`; if [ -z "$el" ]; then error "value for -e option missing"; fi
  nel=`echo $el | awk '{print NF}'`
  if [ $nconc -ne $nel ]; then error "number of concentrations in -c option is not equal to number of elements in -e option"; fi

  # now we loop over the elements (in $el) and average the masses with the concentrations (in $conc)
  Tmax=0
  for (( i=1; i<=$nel; i++ )) do
    e=`echo $el   | awk '{print $'$i'}'`
    x=`echo $conc | awk '{print $'$i'}'`
    melt=`$path/getMeltingPoint.sh $e`
    Tmax=`echo $Tmax $melt $x | awk '{printf "%.7f",$1+$2*$3}'`
  done
fi

if [ $tmaxOp == True ]; then
  Tmax=`getValue -Tmax`
  if [ -z "$Tmax" ]; then error "value for -Tmax option missing"; fi
fi

# determine T start and step
opT=`getOption -Tstart`; if [ $opT == True ]; then Tstart=`getValue -Tstart`; else Tstart=$TstartDef; fi
if [ -z "$Tstart" ]; then error "no value to -Tstart option"; fi
opT=`getOption -Tstep`; if [ $opT == True ]; then Tstep=`getValue -Tstep`; else Tstep=$TstepDef; fi
if [ -z "$Tstep" ]; then error "no value to -Tstep option"; fi

# now calculate configurational entropy and Fconf
s=`echo $conc $nconc | awk '{s=0; for (i=1;i<=$NF;i++) s=s+$i*log($i); print -s}'`
echo $s $Tstart $Tmax $Tstep | awk '{for (t=$2;t<=$3;t=t+$4) printf "%s %.14f\n",t,-t*$1/'$kB'}' > Fconf
echo; echo " configurational entropy for $nconc elements with concentrations $conc"; echo " $s kB"


