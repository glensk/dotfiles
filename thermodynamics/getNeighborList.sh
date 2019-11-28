#!/bin/bash

#-----set parameters and paths------------------------
fileDef=POSCAR
splitScript=splitPOSCAR.sh
#-----------------------------------------------------


# following 3 lines must always be present
path=`set | grep BASH_SOURCE | sed 's|.*\"\(.*\)/[^/]*\".*|\1|' | sed '/BASH_SOURCE=/s/.*/\./'`
script=`set | grep BASH_SOURCE | sed 's|.*\".*/\([^/]*\)\".*|\1|' | sed '/BASH_SOURCE=/s/.*/\./'`
options=$*; . $path/utilities/functions.include; checkOptions "-h -help -i"

# small help for options
# space between option and value must be 1
# space between option and Comment or value and Comment must be >2
h=`getOption -h`
if [ $h == True ]; then
  usage $script
  printOptions "-i inpFile   use inpFile as input POSCAR file (default: $fileDef)"
  exit
fi

# detailed help
help=`getOption -help`
if [ $help == True ]; then
  details $script
  echo2 "   calculates neighbors list for each atom from a POSCAR file"
  echo2 "   format of the output file 'neighList':" \
        "     atom1  shell1  radius11    shellAtom111  shellAtom112 ... shellAtom11n" \
        "     atom1  shell2  radius12    shellAtom121  shellAtom122 ... shellAtom12m" \
        "      ..." \
        "     atom2  shell1  radius21    shellAtom211  shellAtom212 ... shellAtom21o" \
        "      ..." \
        "      ..." \
        "     atomN  shellM  radiusNM    shellAtomNM1  shellAtomNM2 ... shellAtomNMp"
  echo2 "   atom numbers refer to the index in the POSCAR file"
  exit
fi

# mathematica kernel if needed
checkAndSetMath

# if no type given get and check input files
inpFile=`getOption -i`
if [ $inpFile == True ]; then inpFile=`getValue -i`; else inpFile=$fileDef; fi
check $inpFile

# split POSCAR into cell cartesian_coordinates and species files
echo; echo " splitting $inpFile"; $path/$splitScript -i $inpFile -t cartesian > /dev/null

# remove previous output and write info to stdout
rm -f _tmp_math shell_* species_*;
echo; echo " calculating neighbors list ..."

# run mathematica to get shells
$math >> _tmp_math << EOF
<<$path/mathematica/ALL.math;
  cell=checkImport["cell"];
  coords=checkImport["cartesian_coords"];

  allNeighs={};
  Do[
    neighs={};
    dist={};
    Do[
      If[i==j,Continue[]];
      d=periodicDistance[coords[[i]],coords[[j]],cell];
      If[Length[dist]==0,
        AppendTo[dist,d];
        AppendTo[neighs,{j}];
        Continue[];
      ];
      Do[
        If[d+GLOBALEPS<dist[[k]],
          dist=Insert[dist,d,k];
          neighs=Insert[neighs,{j},k];
          Break[];
        ];
        If[Norm[d-dist[[k]]]<GLOBALEPS,
          AppendTo[neighs[[k]],j];
          Break[];
        ];
        If[Length[dist]==k,
          AppendTo[dist,d];
          AppendTo[neighs,{j}];
        ];
      ,{k,Length[dist]}];
    ,{j,coords//Length}];
    Do[AppendTo[allNeighs,{i,j,SetPrecision[dist[[j]],4],neighs[[j]]}//Flatten],{j,dist//Length}];
  ,{i,coords//Length}];

  Export["neighList",allNeighs,"Table"];
EOF
# rm unneeded files
rm -f _tmp_math _tmp_scale cell cartesian_coords species scale
echo; echo " ouput file: neighList"

