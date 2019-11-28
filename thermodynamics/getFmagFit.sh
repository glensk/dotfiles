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
  printOptions "-p base suffix    create parameters.math (see -help) and exit" \
               "-f                force parameters.math creation even if existing"
  exit
fi

# detailed help
help=`getOption -help`
if [ $help == True ]; then
  details $script
  echo2 "   fit a set of magnetic free energy files as a function of "  \
        "   volume for different T using Log[M + 1]"
  echo2 "   parameters.math file must be present (generate with -p)"
        "   -p option can be given a base and suffix which are used to" \
        "   scan for corresonding files base*suffix from which aLats"   \
        "   in parameters.math are automatically generated"
  echo2 "   typical script sequence:"                                   \
        "     ..."                                                      \
        "     ..."                                                      \
        "     \033[1mgetFmagFit.sh\033[0m"                              \
        "     getThermodynamics.sh       (in separate folder)"
  exit
fi

# check if force overwriting templates or folders
overwrite=`getOption -f`

# if applies create parameters.math
genPar=`getOption -p`
if [ $genPar == True ]; then
  if [ -e parameters.math -a $overwrite != True ]; then error "parameters.math exists; use -f to overwrite"; fi
  base=`getValue -p`
  if [ "$base" != "" ]; then
    suffix=`echo $base | awk '{print $2}'`
    base=`echo $base | awk '{print $1}'`
    aLats=`ls -1 ${base}[0-9.][0-9.]*$suffix | sed -e 's/'$base'//' -e 's/'$suffix'//'`
    aLats="`echo $aLats | awk '{printf("{"); for (i=1;i<=NF-1;i++) printf("%s,",$i); printf("%s}",$NF)}'`"
    base="{\"$base\"}"
    suffix="\"$suffix\""
  else
    aLats="Table[a, {a, XXX, XXX, XXX}]"
    base="{\"Fmag_\"}"
    suffix="\"\""
  fi
  echo " (* !!! mathematica format !!!  *)

aLats = $aLats; (* list of lattice constants; use aLats=Table[a,{a,XXX,XXX}]; or aLats={XXX,XXX,...}; *)

latType = \"fcc\";                      (* bcc, fcc, hcp, or size of supercell (2,3,4,...; needed for defects) *)
cBya = 1;                             (* 1 for bcc, fcc, supercell, and actual c/a ratio for hcp *)
s = 1;                                (* scaling factor for electronic free energy; typically 1/nAtoms or *)
                                      (* simply 1 for defect supercells *)

Vorder = 3;                           (* basis for the volume parametrization; typically: 3 giving: {1,V,V^2,V^3} *)
baseName = "$base";                 (* base names for the free energy files; more base names can be given as list*)
suffix = "$suffix";                  (* sufix after aLat such as e.g. angstrom or Ang or blank i.e. \"\" *)

" > parameters.math
  echo; echo "parameters.math written  <==  mind: parameters.math needs to be adjusted"
fi

# exit if we created some template
if [ $genPar == True ]; then exit; fi

# mathematica kernel if needed
checkAndSetMath

# check if parameters.math available
check parameters.math

# do the fitting
echo
$math << EOF
<<$path/mathematica/ALL.math;
<<parameters.math;
fitF[aLats,cBya,latType,s,Vorder,#,suffix,True]&/@baseName;
EOF

