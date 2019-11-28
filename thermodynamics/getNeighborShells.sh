#!/bin/bash

#-----set parameters and paths------------------------
fileDef=POSCAR; nSdef=5; atDef=1
supported="fcc bcc hcp sc"
splitScript=splitPOSCAR.sh
#-----------------------------------------------------


# following 3 lines must always be present
path=`set | grep BASH_SOURCE | sed 's|.*\"\(.*\)/[^/]*\".*|\1|' | sed '/BASH_SOURCE=/s/.*/\./'`
script=`set | grep BASH_SOURCE | sed 's|.*\".*/\([^/]*\)\".*|\1|' | sed '/BASH_SOURCE=/s/.*/\./'`
options=$*; . $path/utilities/functions.include; checkOptions "-h -help -i -t -s -r -c -a -A -R -S -nc -fc -irr -ns -f -l"

# small help for options
# space between option and value must be 1
# space between option and Comment or value and Comment must be >2
h=`getOption -h`
if [ $h == True ]; then
  usage $script
  printOptions "-i inpFile   use inpFile as input POSCAR file (default: $fileDef)" \
               "-t type      alternatively calculate shells for type (supported: $supported)" \
               "-s shell     get neighbor shells up to shell (default: $nSdef)" \
               "-r radius    use radius (Ang) instead of maximum shell for obtaining shells" \
               "-a atom      use the atom'th atom from POSCAR as center of shells (default: $atDef)" \
               "-A           do it for all atoms in POSCAR" \
               "-R           return shells in reduced coordinates (default: cartesian)" \
               "-S           return cartesian coordinates but scaled with the lattice constant from POSCAR" \
               "-nc          no centering of the shells (default: no centering if -A, otherwise centering)" \
               "-fc          force centering of the shells if -A" \
               "-irr symOps  print only irreducible atoms per shell using symOps (see -help)" \
               "-ns          no species info when calculating irreducible atoms (for cluster calculation)" \
               "-f           force overwritting old ouput files" \
               "-l           print just log file"
  exit
fi

# detailed help
help=`getOption -help`
if [ $help == True ]; then
  details $script
  echo2 "   calculates neighboring shells around an atom up to a maximum shell or radius"
  echo2 "   the atomic structure (cell+coordinates) are read in from a POSCAR type file"
  echo2 "   alternatively to the POSCAR input, the '-t type' option can be used to" \
        "   calculate shells for given type (supported: $supported)"
  echo2 "   all atoms (or just the irreducible with -irr) of each shell are written to" \
        "   to a separate file (shell_1, shell_2, ...)"
  echo2 "   to each shell file corresponds a mapping file (mapping_1, mapping_2, ...)" \
        "   if not -irr option is given; the mapping files contain the index of each"  \
        "   atom onto the POSCAR file"
  echo2 "   to each shell file corresponds a species file (species_1, species_2, ...)" \
        "   which contains the number of the species for each atom of the shell" \
        "   the number corresponds to the order of the species in the POSCAR file"
  echo2 "   by default the shells are centered around the origin after they are obtained" \
        "   this can be circumvented by using the -nc option"
  echo2 "   the '-irr symOps' option can be used to extract only the irreducible atoms for each" \
        "   shell; symOps must be either the name of the file containing the symmetry operations" \
        "   of the crystal in cartesian coords (use $path/getCrystalSymmetries.sh)" \
        "   or for symOps=48 all 48 point symmetries are applied (valid e.g. for fcc, bcc, sc)"
  echo2 "   when calculating the irreducible atoms, the -ns option can be used to prevent" \
        "   distinguishing between different species, i.e, all species are effectively the" \
        "   same this is useful for instance when calculating pair cluster figures for a" \
        "   cluster expansion"
  exit
fi

# mathematica kernel if needed
checkAndSetMath

# all atoms option
allOp=`getOption -A`

# check if previous ouput should be overwritten
force=`getOption -f`
if [ $allOp == True ]; then
  previous=`ls -1d atom_[0-9]* 2> /dev/null | awk 'END{print NR}'`
else
  previous=`ls atomsPerShell radiusPerShell shell_[0-9]* species_[0-9]* 2> /dev/null | awk 'END{print NR}'`
fi
if [ $previous != 0 -a $force != True ]; then error "previous ouput files exist; use -f to force overwritting"; fi

# check if -t type option given
tOp=`getOption -t`
if [ $tOp == True ]; then
  strType=`getValue -t`
  c=`echo $supported | xargs -n1 | grep "\<$strType\>"`
  if [ "$c" != $strType ]; then error "type not supported (supported: $supported)"; fi
else
  # if no type given get and check input files
  inpFile=`getOption -i`
  if [ $inpFile == True ]; then inpFile=`getValue -i`; else inpFile=$fileDef; fi
  check $inpFile

  # split POSCAR into cell cartesian_coordinates and species files
  echo; echo " splitting $inpFile"; $path/$splitScript -i $inpFile -t cartesian > /dev/null

  # check if atom option given
  aOp=`getOption -a`
  if [ $aOp == True ]; then
    atom=`getValue -a`; c=`checkInteger $atom`
    if [ $c != ok ]; then error "value for atom in -a option empty or wrong"; fi
  else
    atom=$atDef;
  fi

  # check if scaled coordinates option given
  scaled=`getOption -S`
  if [ $scaled == True ]; then
    if [ $reduced == True ]; then error "choose only one of the -R and -S option"; fi
    awk 'NR==2{print}' $inpFile > _tmp_scale
  fi
fi

# check if reduced coordinates, no centering, or irreducible atom option given
reduced=`getOption -R`; nocenter=`getOption -nc`; irrOp=`getOption -irr`
if [ $irrOp == True ]; then
  irr=`getValue -irr`;
  if [ -z "$irr" ]; then error "no value given to -irr option"; fi;
  if [ $irr != 48 ]; then check $irr; fi;
  ns=`getOption -ns`
fi

# compute clusters up to given shell or within given radius
op=`getOption -r`;
if [ $op == True ]; then
  value=`getValue -r`; c=`checkReal $value`
  if [ $c != ok ]; then error "value for radius in -r option empty or wrong"; fi
  type="radius"
else
  op=`getOption -s`
  if [ $op == True ]; then
    value=`getValue -s`; c=`checkInteger $value`
    if [ $c != ok ]; then error "value for nShell in -s option empty or wrong"; fi
  else
    value=$nSdef;
  fi
  type="shell"
fi

# -A option cannot be run with some other options
if [ $allOp == True -a $irrOp == True ]; then error "-A option only without -irr option possible"; fi
if [ $allOp == True -a   $tOp == True ]; then error "-A option only without '-t type' option possible"; fi
if [ $allOp == True -a "$aOp" == True ]; then error "-A option only without '-a atom' option possible"; fi

# for -A option we do not center and we need to create folders for each atom
if [ $allOp == True ]; then
  nocenter=True;
  nAtoms=`awk 'END{print NR}' cartesian_coords`
  for (( i=1; i<=$nAtoms; i++)); do rm -fr atom_$i; mkdir atom_$i; done
  dir=`pwd`
fi

# force centering if -fc option given
fcOp=`getOption -fc`
if [ $fcOp == True ]; then nocenter=False; fi

# remove previous output and write info to stdout
rm -f _tmp_math shell_* species_*;
echo; echo " calculating shells up to $type $value"; echo; echo " ..."

# run mathematica to get shells
$math >> _tmp_math << EOF
<<$path/mathematica/ALL.math;

(* -------- exportShells module start -------- *)
exportShells:=Module[{},                           (* small module for exporting shells *)
  If["$nocenter"=="True",                          (* if -nc option is given we need to shift back because we are centered at this point  *)
    shells=#-coords[[atom]]&/@#&/@shells
  ];
  If["$reduced"=="True",                           (* transform back to reduced if -R option given *)
    shells=toReducedCoords[#,cell]&/@shells;
  ];
  If["$scaled"=="True",                            (* -S option --> division of all shells by scale *)
    shells=1/scale*#&/@#&/@shells;
  ];
  Do[                                              (* export shell_* and species_* files and exit *)
    str=ToString[PaddedForm[i,                     (* to get nice file names we pad with 0s *)
      StringLength[ToString[Length[shells]]],      
      NumberPadding->"0",NumberSigns->{"",""}]];   (* turn off number sign because mathematica leaves otherwise space for the sign *)
    Export["shell_"<>str,
      Append[shells[[i]],{}],"Table"];
    If["$ns"!="True",                              (* export species only if we are not in -ns mode which does not distinguish species *)
      Export["species_"<>str,
        Append[species[[i]],{}],"Table"]
    ];
    If["$irrOp"!="True",
      Export["mapping_"<>str,
        Append[mapping[[i]],{}],"Table"];
    ];
  ,{i,shells//Length}];
];
(* -------- exportShells module end ---------- *)

If["$tOp"=="True",
  coordsIn={{0,0,0}};
  speciesIn={1};
  atom=1;
  Switch["$strType",                                 (* if type given, choose out of the available *)
    "fcc",cell=fccCell[],
    "bcc",cell=bccCell[],
    "sc",cell=scCell[],
    "hcp",cell=hcpCell[];
          coordsIn=hcpCoordsRel;
          speciesIn={1,1},
    _,error["cell type not known"]
  ];
,
  cell=checkImport["cell"];
  coordsIn=checkImport["cartesian_coords"];
  speciesIn=checkImport["species"]//Flatten;
  coordsIn=toReducedCoords[coordsIn,cell];             (* input coords to getNeighborShells need to be in reduced coordinates *)
  atom=1 $atom;
];

cont=False;
While[cont==False,
  {dist,n,shells,species,mapping}=
    getNeighborShells[cell,coordsIn,speciesIn,           (* here the actual work is done *)
                      atom,$value,"$type"];

  If["$allOp"=="True",
    dir="$dir/atom_"<>ToString[atom];
    SetDirectory[dir];
  ];

  Export["radiusPerShell",dist,"Table"];             (* dist (radius of shells) and n (number of atoms per shell) are not modified *)
  Export["atomsPerShell",n,"Table"];                 (* below so we can export already *)

  shells=toCartesianCoords[#,cell]&/@#&/@shells;     (* we work in cartesian coordinates from now on *)
  coords=toCartesianCoords[#,cell]&/@coordsIn;

  If["$scaled"=="True",                              (* this is for -S option; scaling with lattice constant from POSCAR *)
    scale=Import["_tmp_scale","List"]//Flatten;      (* _tmp_scale was prepared above *)
    If[Length[scale]==1,scale=scale[[1]]];           (* scale can have either 1 or 3 values *)
    If[scale<0,scale=(-scale)^(1/3)];                (* if negative we have volume *)
  ];

  shells=#-coords[[atom]]&/@#&/@shells;              (* centering around origin for symmetrization *)
  exportShells;

  If["$allOp"!="True",cont=True];
  If[atom==Length[coords],cont=True,atom+=1];
];
If["$irrOp"!="True",Exit[]];                         (* if we are not reducing to irreducible not much modfication is needed --> quick exit *)

(* ----------- END here if not irrOp -------------- *)


If["$irr"=="48",S=symOps,S=checkImport["$irr"]];   (* from here, we are in irreducible mode; if irr==48 we take all 48 *)
                                                   (* point symmetry operations otherwise we read them from file contained in irr *)
If["$ns"=="True",species*=0];                      (* if -ns option we do not distinguish species and set therefore all to 0 *)

newShells=Table[{},{shells//Length}];
newSpecies=Table[{},{shells//Length}];

Do[                                                (* loop over all shells *)
   union=Union[species[[i]]];                      (* we will symmetrize only among the same species *)
   subshells=Table[{},{union//Length}];            (* we therefore subdivide into subshells containing each the same species *)
   Do[ Do[
     If[species[[i,j]]==union[[k]],
       AppendTo[subshells[[k]],shells[[i,j]]];
       Break[];
     ];
   ,{k,union//Length}],{j,species[[i]]//Length}];

   Do[                                             (* do symmetrization for each subshell *)
     CELL=10 cell;        (*  IMPORTANT  *)        (* we need to increase the cell because getIrreducibleAtoms will otherwise detect the atoms *)
                                                   (* to be periodically equal; we want however only pure point symmetries to be in effect *)
     irrIndMat=
       getIrreducibleAtoms[subshells[[j]],         (* actual symmetrization done here; shells and symmetries S need to be cartesian *)
       CELL,S,"irreducibleStars"][[1]];            (* we need only the first element returned which are the indices grouped in reducible groups *)
     Do[ Do[                                       (* try to find for each reducible star irreducible representative with positive coords *)
       n=irrIndMat[[k,l]];
       If[l==Length[irrIndMat[[k]]]||              (* we append if positive *)
         Abs[subshells[[j,n]]]==subshells[[j,n]],  (* or if we are at last element *)
         AppendTo[newShells[[i]],subshells[[j,n]]];
         AppendTo[newSpecies[[i]],union[[j]]];
         Break[];
       ];
     ,{l,irrIndMat[[k]]//Length}],{k,irrIndMat//Length}];
   ,{j,subshells//Length}];
,{i,shells//Length}];

Export["irrPerShell",Length/@newSpecies,"Table"];
shells=newShells; species=newSpecies; exportShells
EOF

# rm unneeded files
rm -f _tmp_math _tmp_scale cell cartesian_coords species scale

# no log file if -A option; quit then
if [ $allOp == True ]; then 
  echo " output to folders written:"; ls -d atom_[0-9]*
  exit
fi

# additional output if irreducible option
if [ $irrOp == True ]; then str="irr"; else str=""; fi

# write relevant ouput to stdout
if [ $tOp == True ]; then unit=aLat; else unit=Ang; fi
echo > log_getNeighborShells; echo " nrShell radius($unit) atoms $str" >> log_getNeighborShells
if [ $irrOp == True ]; then
  paste radiusPerShell atomsPerShell irrPerShell | awk '{printf("%5d  %9.3f  %6d  %3d\n",NR,$1,$2,$3,$4)}' >> log_getNeighborShells
else
  paste radiusPerShell atomsPerShell | awk '{printf("%5d  %9.3f  %6d\n",NR,$1,$2,$3)}' >> log_getNeighborShells
fi
cat log_getNeighborShells

# if -l option we keep only log file
log=`getOption -l`
if [ $log == True ]; then rm -f atomsPerShell radiusPerShell irrPerShell mapping_[0-9]* shell_[0-9]* species_[0-9]*; fi

# write which files have been produced
echo; echo " files written:"; ls atomsPerShell radiusPerShell irrPerShell log_getNeighborShells 2> /dev/null;
ls shell_[0-9]* 2> /dev/null; ls mapping_[0-9]* 2> /dev/null; ls species_[0-9]* 2> /dev/null


