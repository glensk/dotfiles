#!/bin/bash

# following 3 lines must always be present
path=`set | grep BASH_SOURCE | sed 's|.*\"\(.*\)/[^/]*\".*|\1|' | sed '/BASH_SOURCE=/s/.*/\./'`
script=`set | grep BASH_SOURCE | sed 's|.*\".*/\([^/]*\)\".*|\1|' | sed '/BASH_SOURCE=/s/.*/\./'`
options=$*; . $path/utilities/functions.include; checkOptions "-h -help -p -f"

# small help for options
# space between option and value must be 1
# space between option and Comment or value and Comment must be >2
h=`getOption -h`
if [ $h == True ]; then
  usage $script
  printOptions "-p        create parameters.math and exit" \
               "-f        force parameters.math creation even if existing"
  exit
fi

# detailed help
help=`getOption -help`
if [ $help == True ]; then
  details $script
  echo2 "   fit a mesh of Fel(V,T) points"
  echo2 "   parameters.math file must be present (generate with -p)"
  echo2 "   typical script sequence:"                                   \
        "     createFolders_Fel.sh"                                     \
        "     collectOUTCARs.sh          (one folder higher from here)" \
        "     extractFel.sh"                                            \
        "     \033[1mgetFelFit.sh\033[0m"                               \
        "     getThermodynamics.sh       (in separate folder)"
  exit
fi

# check if force overwriting templates or folders
overwrite=`getOption -f`

# if applies create parameters.math
genPar=`getOption -p`
if [ $genPar == True ]; then
  if [ -e parameters.math -a $overwrite != True ]; then error "parameters.math exists; use -f to overwrite"; fi
      
      if [ ! -e "extractFel.log" ]; then
        echo
        echo "  extractFel.log missing (extractFel.sh not used?)"
        echo "  --> lattype will be missing in parameters.dat"
        lattype=XXX
        if [ `ls -1d Fel_*Ang 2> /dev/null | wc -l` == 0 ]; then
          echo
          echo "  Fel_*Ang files also missing"
          echo "  --> alats will be missing in parameters.dat"
          ang=XXX
        else
          ang=`ls -1d Fel_*Ang | sed 's|Fel_||' | sed 's|Ang|,|'`
        fi
      else
        ## get some of the values
        ang=`grep -v "^#" extractFel.log | awk '{print $1","}'`
        lattype=`grep "^# lattype: " extractFel.log | awk '{print $NF}'`
      fi
        aLats=`echo "aLats = {"$ang"};" | sed 's|,}|}|' | sed 's|[ ]*||'`
  echo "
$aLats (* list of lattice constants; use aLats={XXX,XXX,...}  alternatively   *)
maxSigma = \"all\";                     (* keep \"all\" if all temperatures are to be read, else truncate (eV) *)

latType = \"$lattype\";                      (* bcc, fcc, hcp, dfcc (double fcc), omega, or size of supercell (2,3,4,...; needed for defects) *)
                                      (* latType=\"aLatsAreVolumes\" is also possible to assume aLats to be volumes *)
cBya = 1;                             (* 1 for bcc, fcc, supercell, and actual c/a ratio for hcp *)
s = 1;                                (* scaling factor for electronic free energy; typically 1/nAtoms or *)
                                      (* simply 1 for defect supercells *)

Tmin  = 1;                            (* start T (K) for the final T mesh           *)
Tmax  = XXX;                          (* end   T (K); typically melting temperature *)
Tstep = 2;                            (* step  T (K); 2 K has been always fine      *)

Tbasis = {1,T,T^2,T^3,T^4};           (* basis for the DOS expansion along T; typically: {1,T,T^2,T^3,T^4}        *)
Vorder = 3;                           (* basis for the volume parametrization; typically: 3 giving: {1,V,V^2,V^3} *)

" > parameters.math
  echo; echo "parameters.math written  <==  mind: parameters.math needs to be adjusted"
fi

# exit if we created some template
if [ $genPar == True ]; then exit; fi

# mathematica kernel if needed
checkAndSetMath

# check if parameters.math available
check parameters.math

# check if at least Vorder+1 Fel files
ord=`grep Vorder parameters.math | sed 's/.*=.*\([1-9]\);.*/\1/'`
nfiles=`ls -1 Fel_*Ang | wc -l`
if [ $ord -ge $nfiles ]; then error "not enough Fel_*Ang files for Vorder"; fi

# do the fitting
echo; echo " ... fitting ..."; echo;
$math << EOF
<<$path/mathematica/ALL.math;
<<parameters.math;
fitFel[aLats,cBya,latType,s,maxSigma,Tmin,Tmax,Tstep,Tbasis,Vorder]
EOF

