#!/bin/bash

#-----set parameters and paths-------------------------------------
fileDef=POSCAR; shellDef=13
splitScript=splitPOSCAR.sh; getClustersScript=getNeighborShells.sh
#------------------------------------------------------------------


# following 3 lines must always be present
path=`  set | grep BASH_SOURCE | sed 's|.*\"\(.*\)/[^/]*\".*|\1|' | sed '/BASH_SOURCE=/s/.*/\./'`
script=`set | grep BASH_SOURCE | sed 's|.*\".*/\([^/]*\)\".*|\1|' | sed '/BASH_SOURCE=/s/.*/\./'`
options=$*; . $path/utilities/functions.include; checkOptions "-h -help -i -s -c -S -k"

# small help for options
# space between option and value must be 1
# space between option and Comment or value and Comment must be >2
h=`getOption -h`
if [ $h == True ]; then
  usage $script
  printOptions "-i inpFile  use inpFile as input POSCAR file (default: $fileDef)" \
               "-s shell    get occupations for clusters up to shell (default: $shellDef)" \
               "-c cFile    instead of calculating clusters read from cFile (cartesian coords)" \
               "-S sFile    use sFile as symmetry operations (default: all 48 point symmetries)" \
               "-k          keep clusters file (if created) and occupations file"
  exit
fi

# detailed help
help=`getOption -help`
if [ $help == True ]; then
  details $script
  echo2 "   calculates occupation numbers of pair clusters (figures)"
  echo2 "   the atomic structure (cell+coordinates) are read in from a POSCAR type file"
  echo2 "   clusters are calculated up to given shell ('-s shell' option) or alternatively" \
        "   can be supplied in a clusters file ('-c cFile' option) containing the" \
        "   clusters in cartesian coordinates"
  echo2 "   by default all 48 point symmetries are used as crystal symmetries;" \
        "   alternatively a symOps file can be provided ('-S symOps' option) containing" \
        "   the crystal symmetries in cartesian coordinates, with one symmetry operation" \
        "   per line without commas or brackets, e.g., identity and mirror symmetry:" \
        "     1  0  0    0  1  0    0  0  1" \
        "    -1  0  0    0 -1  0    0  0 -1" \
        "     ..."
  echo2 "   the script is working only for symmorphic symmetries, i.e., point symmetries" \
        "   this means that, e.g., for structures like hcp which contains non-symmorphic" \
        "   symmetries related to the second atom in the unit cell those cannot be treated"
  echo2 "   by default only a log file is written containing clusters and occupations" \
        "   to keep the clusters file (if created) and occupations file use -k option"
  exit
fi

# mathematica kernel if needed
checkAndSetMath

# get and check input files
inpFile=`getOption -i`; sFile=`getOption -S`;
if [ $inpFile == True ]; then inpFile=`getValue -i`; else inpFile=$fileDef; fi
check $inpFile
if [ $sFile == True ]; then sFile=`getValue -S`; check $sFile; else sFile=48; fi
echo; echo " input:  $inpFile  $sFile"

# check if we import clusters (option -c) or calculate
cOp=`getOption -c`;
if [ $cOp == True ]; then
  cFile=`getValue -c`; check $cFile
  echo; echo " reading clusters from file $cFile"
else
  # check maximum shell to calculate clusters (option -s)
  op=`getOption -s`
  if [ $op == True ]; then
    shell=`getValue -s`; c=`checkInteger $shell`
    if [ $c != ok ]; then error "value for shell in -s option empty or wrong"; fi
  else
    shell=$shellDef;
  fi
  echo; echo " calculating clusters ..."
  $path/$getClustersScript -i $inpFile -s $shell -irr $sFile -ns -f > /dev/null
  cat shell_* > clusters
  rm shell_* irrPerShell atomsPerShell radiusPerShell log_getNeighborShells
  cFile=clusters
fi

# split POSCAR into cell cartesian_coordinates and species files
$path/$splitScript -i $inpFile -t cartesian  > /dev/null

# run mathematica to get occupations
rm -f _tmp_math; echo; echo " calculating occupations ..."
$math >> _tmp_math << EOF
  <<$path/mathematica/ALL.math;

  cell=checkImport["cell"];                                       (* cell cartesian_coords species scale files *)
  coords=checkImport["cartesian_coords"];                         (* are coming from splitScript *)
  species=checkImport["species"]//Flatten;
  scale=checkImport["scale"]//Flatten;
  clusters=checkImport["$cFile"];

  If["$sFile"=="48",S=symOps,S=Partition[#,3]&/@checkImport["$sFile"]];            (* symOps contains all 48 point operations *)
  occ=getPairOccupationNumbers[cell,coords,species,clusters,S,0.01];   (* actual work done here; coords, clusters, S need to be cartesian *)

  clusters=1/scale #&/@clusters;                                  (* scale clusters for convenient output *)
  Export["$cFile"<>"_scaled",clusters//N,"Table"];
  Export["occupations",occ//N,"Table"];
EOF
if [ "`grep ERROR _tmp_math`" != "" ]; then error "check _tmp_math"; fi
rm _tmp_math cell cartesian_coords species scale

# write results to stdout
echo > log_getPairOccupation; echo "     cluster        occupation" >> log_getPairOccupation
paste ${cFile}_scaled occupations | awk '{printf(" %5.2f %5.2f %5.2f  %5.2f\n",$1,$2,$3,$4)}' >> log_getPairOccupation
cat log_getPairOccupation

# check if -k option given to keep the occupations file
occ=`getOption -k`
if [ $occ == True ]; then
  # check if we have created a clusters file
  if [ $cOp == True ]; then str="  $cFile"; else str=""; fi
  echo; echo " output:$str  $occ  log_getPairOccupation"
else
  rm occupations
  # if clusters file has been created remove it since no -k option
  if [ $cOp != True ]; then rm $cFile{,_scaled}; fi
  echo; echo " output:  log_getPairOccupation"
fi


