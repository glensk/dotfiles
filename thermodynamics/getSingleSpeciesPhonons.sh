#!/bin/bash

#-----Paths---------------------------------------------
kpMeshFile=irreducible16x16x16kPmesh.dat
kpMeshWeightFile=irreducible16x16x16kPmesh_Weights.dat
#-------------------------------------------------------


path=`set | grep BASH_SOURCE | sed 's|.*\"\(.*\)/[^/]*\".*|\1|' | sed '/BASH_SOURCE=/s/.*/\./'`
script=`set | grep BASH_SOURCE | sed 's|.*\".*/\([^/]*\)\".*|\1|' | sed '/BASH_SOURCE=/s/.*/\./'`
options=$*; . $path/utilities/functions.include; checkOptions "-h -help -D -H -fe -fm -fa -f -Fe -Fm -F -P -A -a -p -c -afm -C -o -v -vv -vac -disp -ge -gc -gs -ga -gt -gn -gm -gbd -folder -makefits"



h=`getOption -h`
if [ $h = True ]; then
  usage $script
  printOptions "-folder [dir]   evaluate not in parent working directory (pwd) but in specified folder" \
               "-D              DynMatR" \
               "-H              HesseMatrix" \
               "-fe             ExactFreqs" \
               "-fm             MeshFreqs" \
               "-fa             MeanFreqs (average)" \
               "-f              ExactFreqs+MeshFreqs" \
               "-Fe             Fqh_fromExactFreqs" \
               "-Fm             Fqh_fromMeshFreqs" \
               "-F              Fqh_fromExactFreqs+Fqh_fromMeshFreqs" \
               "-P              PhononDispersion" \
               "-A              all above" \
               "-a aLat         override aLats from parameters.dat with aLat" \
               "-p              create parameters.dat template and exit" \
               "-c              switch to classical free energy (default: qm)" \
               "-afm            switch to afm case with forcesX!=forcesZ (see -help)" \
               "-C coords       use coords file instead of the default coordinates" \
               "-o              overwrite parameters.dat if it exists" \
               " " \
               "-vac            get HesseMatrix and Fvib for vacancy calculation" \
               "-disp [value]   defines displacement in bohrradius of vacancy calculation" \
               " " \
               "-ge             echo element and exit" \
               "-gc             echo cellType and exit" \
               "-gs             echo supercell and exit" \
               "-ga             echo alats and exit" \
               "-gt             echo tmelt and exit" \
               "-gn             echo number of atoms and exit" \
               "-gm             echo mass and exit" \
               "-gbd            echo bulk (if bulkcell) or def (for defect cell)" \
               " " \
               "-makefits       create fitperatom and fitsupercell Folder with corresponding fit " \
               " " \
               "-v              be verbose" \
               "-vv             be even more verbose"
  echo "Note:    at least one option (other than -a, -c, -afm) must be given" 1>&2
  echo "Example: getSingleSpeciesPhonons.sh -H -fe -Fe" 1>&2
  exit
fi

help=`getOption -help`
if [ $help = True ]; then
  details $script
  echo2 "   calculates harmonic phonons and free energies from finite differences" \
        "   of input forces for single species materials of fcc or bcc structure"
  echo2 "   parameters.dat file must exist containing (mathematica format, but # comments):" \
        "      element=\"XX\"          # e.g., XX=Ca" \
        "      cellType=\"XX\"         # XX=fcc or bcc" \
        "      supercell=XX          # XX=2,3,4,..." \
        "      aLats={a1,a2,a3,...}  # with a1 ... lattice constants in Ang" \
        "      TRange={T1,T2,dT}     # optionally: start/end/step temperature T1/T2/dT in K"
  echo2 "   for each lattice constant {a1,a2,...} a forces.a1 and disp.a1 file must exist"
  echo2 "   a forces file contains forces in eV/Ang; each force (x,y,z) on a single line" \
        "   the order of the forces must be consistent with the coordinates_SCxSCxSCsc file" \
        "   which is placed in $path/utilities/cellType/" \
        "   (SC=supercell, cellType=fcc or bcc)"
  echo2 "   the '-C coords' option can be used to provide a different coordinates file" \
        "   the supplied coordinates file must be consistent with the given supercell and" \
        "   cellType, but it can contain the atomic coordinates in a different order" \
        "   the displaced atom must be however on the first position"
  echo2 "   a disp file contains a single number which gives the displacement in Ang"
  echo2 "   to extract forces and disp files with this format from corresponding OUTCARs use:" \
        "   $path/extractForces.sh"
  echo2 "   to create vasp input for runs generating the necessary OUTCARs use:" \
        "   $path/createFolders_dynMat.sh"
  echo2 "   for -fe and -Fe an exactRedKPoints_SCxSCxSCsc_noGamma.dat must be available in" \
        "   $path/utilities/cellType/"
  echo2 "   for -fm and -Fm an $kpMeshFile and $kpMeshWeightFile" \
        "   must be available in $path/utilities/cellType/"
  echo2 "   to get quatities at a lattice constant aLatNew which lies in between the calculated" \
        "   ones (e.g. interpolated phonon dispersion at aLatNew) make sure that all displacements" \
        "   in the disp.* files are equal and use the following:" \
        "   (disp.aLat is thereby any one of the available disp files)"
  echo2 "       $path/fitForces.sh -B forces. -a aLatNew" \
        "       cp disp.aLat disp.aLatNew" \
        "       $script -P -a aLatNew"
  echo2 "   the -afm option allows to evaluate the quantities for an afm spin type configuration" \
        "   effectively this means that instead of a forces.a1 one needs a forcesX.a1 and forcesZ.a1" \
        "   where the first file contains the forces when the atom is displaced within the plane" \
        "   of same spin and where the second file contains forces when the atom is displaced out" \
        "   of the plane; forcesY.a1 is obtained by symmetry from forcesX.a1"
  echo2 "   typical script sequence:"                                   \
        "     createFolders_dynMat.sh -a ..."                           \
        "     createFolders_dynMat.sh"                                  \
        "     collectOUTCARs.sh"                                        \
        "     extractForces.sh"                                         \
        "     \033[1mgetSingleSpeciesPhonons.sh -A\033[0m                  (cmpc01)"  \
        "     getFqhFit.sh                                   (cmpc01)"  \
        "     getThermodynamics.sh               (in separate folder)"
  exit
fi
checkAndSetMath
# check if at least one option
if [ "$#" = "0" ]; then error "at least one option needed (use -h for help)"; fi

####################################################################
# -folder [dir] : go to folder if specified
####################################################################
if [ "`getOption -folder`" = "True" ];then
    goto=`getValue -folder`
    [ ! -e "$goto" ] && error "folder $goto does not exist" && exit
    cd $goto
fi
[ "`getOption -v`" = "True" ] && echo hier: `pwd`

#########################################################################################
# -p -gt -gn : get necessary information so we dont have to fill parameters.dat manually
#########################################################################################
genPar=`getOption -p`
overwrite=`getOption -o`

#[ "`getOption -gt`" = "True" ] || [ "`getOption -gn`" = "True" ] || [ "`getOption -gbd`" = "True" ] || [ "`getOption -makefits`" = "True" ] && outputinformation="True"
[ "`getOption -gn`" = "True" ] || [ "`getOption -gbd`" = "True" ] || [ "`getOption -makefits`" = "True" ] && outputinformation="True"
if [ "$genPar" = "True" ] || [ "$outputinformation" = "True" ] ; then
        # get element
        [ "`getOption -v`" = "True" ] && echo in
        
        # get element and atoms
        outcar_tmp=`find -L . -maxdepth 3 -mindepth 1 -type f -name "OUTCAR*" -print -quit`   # depth 3 necessary for defects sometimes
        [ "`getOption -v`" = "True" ] && echo outcar_tmp: $outcar_tmp
        element=XX;atoms=XX;element_tmp="";atoms_tmp="";
        [ -e "$outcar_tmp" ] && element=`OUTCAR_elements.sh $outcar_tmp` && atoms=`OUTCAR_number_of_atoms.sh $outcar_tmp`;

        [ "`getOption -v`" = "True" ] && echo atoms:$atoms:
        [ "`getOption -gn`" = "True" ] && [ "`isnumber.sh $atoms`" = "yes" ] && echo $atoms && exit
        [ "`getOption -ge`" = "True" ] && [ "$element" != "" ] && echo $element && exit
        [ ! -e "$outcar_tmp" ] && [ -e "POTCAR" ] && element=`POTCAR_element.sh POTCAR`
        [ "`getOption -ge`" = "True" ] && [ "$element" != "" ] && echo $element && exit

        tmelt=1000;tmelt_tmp=`getMeltingPoint.sh $element -r`;[ "`isnumber.sh $tmelt_tmp`" = "yes" ] && tmelt=$tmelt_tmp
        [ "`getOption -v`" = "True" ] && echo tmelt:$tmelt
        [ "`getOption -gt`" = "True" ] && echo $tmelt && exit
        
        # alats
        alats="a1,a2,...";
        outcar_alats=`find -L . -maxdepth 2 -mindepth 1 -type f -name "OUTCAR*" | sed 's|\.gz||' | grep -o "[0-9]*[.][0-9]*Ang\|OUTCAR.[0-9]*[.][0-9]*" | sed 's|Ang||g' | sed 's|OUTCAR.||g' | sort -n | uniq | xargs`
        [ "`getOption -v`" = "True" ] && echo outcar_alats:$outcar_alats:
        folder_alats=`find . -maxdepth 1 -mindepth 0 -type d -name "*Ang" | sed 's|./||g' | sed 's|Ang||g' | xargs`
        [ "`getOption -v`" = "True" ] && echo folder_alats:$folder_alats:

        [ "`isnumber.sh $outcar_alats`" != "yes" ] && outcar_alats=""
        [ "`isnumber.sh $folder_alats`" != "yes" ] && folder_alats=""
        all_alats=`echo "$outcar_alats $folder_alats" | xargs -n1 | sort -n | uniq | xargs`
        [ "`getOption -v`" = "True" ] && echo all_alats:$all_alats:
        #[ "$outcar_alats" = "" ] && outcar_alats=`echo "$outcar_all" | sed 's|Ang.*||' | xargs -n1 | sed 's|./||' | xargs`
        


        alats=`echo "$all_alats" | sed 's| |,|g'`
        tmelt=1000;tmelt_tmp=`getMeltingPoint.sh $element -r`;[ "`isnumber.sh $tmelt_tmp`" = "yes" ] && tmelt=$tmelt_tmp
        #echo na:$atoms_tmp:
        ct=XX;sc=XX;bulkdef=XX
        atoms_tmp=$atoms
        if [ "`isnumber.sh $atoms_tmp`" = "yes" ];then
            fccatoms=(4 32 108 256 500)
            fccatomsvac=(3 31 107 255 499)
            fccsc=(1 2 3 4 5)
            bccatoms=(2 16 54 128 249)
            bccatomsvac=(1 15 53 127 248)
            bccsc=(1 2 3 4 5)
            sequence=`awk 'BEGIN {x=-1; while(++x<='"4"'){print x; }; exit}'`
            for s in `echo $sequence`;do
                f=`echo ${fccatoms[$s]}`
                fv=`echo ${fccatomsvac[$s]}`
                b=`echo ${bccatoms[$s]}`
                bv=`echo ${bccatomsvac[$s]}`
                [ "$f" = "$atoms_tmp" ] && atoms=$atoms_tmp && ct=fcc && sc=`echo ${fccsc[$s]}` && bulkdef=bulk
                [ "$fv" = "$atoms_tmp" ] && atoms=$atoms_tmp && ct=fcc && sc=`echo ${fccsc[$s]}` && bulkdef=def
                [ "$b" = "$atoms_tmp" ] && atoms=$atoms_tmp && ct=bcc && sc=`echo ${bccsc[$s]}` && bulkdef=bulk
                [ "$bv" = "$atoms_tmp" ] && atoms=$atoms_tmp && ct=bcc && sc=`echo ${bccsc[$s]}` && bulkdef=def
                #echo s:$s f:$f atoms:$atoms ct:$ct  sc:$sc
            done
            fi
            
            # this is necessary in case parameters.dat does not exist
            [ "`getOption -gt`" = "True" ] && echo $tmelt && exit
            [ "`getOption -gbd`" = "True" ] && echo $bulkdef && exit
            [ "`getOption -gn`" = "True" ] && echo $atoms && exit
        fi

#####################################################################
# write parameters.dat 
#####################################################################
if [ "$genPar" = "True" ] ; then
if [ -e parameters.dat -a $overwrite != True ]; then error "parameters.dat existing; use -o to overwrite"; fi
echo "
element=\"$element\"          # e.g., XX=Ca
cellType=\"$ct\"         # XX=fcc or bcc
supercell=$sc          # XX=2,3,4,...
aLats={$alats}     # with a1 ... lattice constants in Ang
TRange={1,"$tmelt",1}     # optionally: start/end/step temperature T1/T2/dT in K
" > parameters.dat
 echo; echo "parameters.dat written"; 
 [ "`getOptionNr`" = "1" ] && exit # but actually shouldnt leave if other options specified
fi


#####################################################################
# read from parameters.dat :  element, cellType, and supercell size
#####################################################################
[ ! -e "parameters.dat" ] && echored "ERROR: please create parameters.dat in `pwd` first using `basename $0` -p" && exit
check parameters.dat
element=`get element`; cellType=`get cellType`; supercell=`get supercell`;sc=`get sc`;
[ "$supercell" = "" ] && [ "$sc" != "" ] && supercell=$sc
[ "$element" = "" ] && error "please define element in parameters.dat (e.g. element=\"Al\")"
[ "$supercell" = "" ] && error "please define supercell in parameters.dat (e.g. sc=3) for a 3x3x3 supercell"
[ "$cellType" = "" ] && error "please define cellType in parameters.dat (e.g. cellType=\"fcc\")"
[ "`echo $cellType | grep -o '"' | wc -w | sed 's|[ ]*||g'`" = "0" ] && cellType=\"$cellType\"   # make "bcc" out of bcc; necessary when reading input of createFolder_dynmat_vacancy.sh
[ "`echo $element | grep -o '"' | wc -w | sed 's|[ ]*||g'`" = "0" ] && element=`echo "$element" | python -c "print raw_input().capitalize()"` && element=\"$element\"
if [ "`getOption -v`" = "True" ];then
    echo cellType:$cellType
    echo element:$element
    echo supercell:$supercell
fi
checkInput "$element" "$cellType" "$supercell"
# get atomic mass
element=`echo $element | sed 's/"//g'`

# get tmelt (has to be after element got rid of "")
check $path/getAtomicMass.sh; mass=`$path/getAtomicMass.sh $element`
tmelt=1000;tmelt_tmp=`getMeltingPoint.sh $element -r`;[ "`isnumber.sh $tmelt_tmp`" = "yes" ] && tmelt=$tmelt_tmp
[ "`getOption -v`" = "True" ] && echo tmelt:$tmelt tmelt_tmp:$tmelt_tmp
[ "`getOption -gt`" = "True" ] && echo $tmelt && exit

######################################################
# get aLats , element, atomicmass
######################################################
# check if -a option given
aa=`getOption -a`
if [ $aa = True ]; then
  # override aLats from parameters.dat
  aa=`getValue -a`
  aLats=`echo $aa | awk '{printf("{"); for (i=1;i<NF;i++) printf("%s,",$i); printf("%s}",$NF)}'`
else
  # get aLats from parameters.dat otherwise
  aLats=`get aLats`; checkInput "$aLats"
fi
[ "`echo $aLats | grep -o '{' | wc -w | sed 's|[ ]* ||g'`" = "0" ] && aLats=`echo "{\`echo $aLats | sed 's| |,|g'\`}"`


##########################################################################################
# -fa (MeanFreqs)
##########################################################################################
if [ "`getOption -fa`" = "True" ];then
        for a in `echo $aLats | sed 's|{||' | sed 's|,| |g' | sed 's|}||'`;do
            [ ! -e "MeshFreqs_$a" ] && options="$options -fm" && meshFreqs=`getOption -fm`
            [ ! -e "ExactFreqs_$a" ] && options="$options -fe" && exactFreqs=`getOption -fe`
        done
    fi
    
##########################################################################################
# -makefits (checks)
##########################################################################################
if [ "`getOption -makefits`" = "True" ];then
    ## make checks
        for a in `echo $aLats | sed 's|{||' | sed 's|,| |g' | sed 's|}||' | sed 's|;||'`;do             
            #[ -e "Fqh_fromExactFreqs_$a" ] && echo exists  one could think to exclude existg files from recalculation ...
            [ ! -e "Fqh_fromExactFreqs_$a" ] && [ "`getOption -Fe`" != "True" ] && options="$options -a $a -Fe" && echored "redo Fqh_fromExactFreqs_$a"
            [ ! -e "Fqh_fromMeshFreqs_$a" ] && options="$options -a $a -Fm"  && echored "redo Fqh_fromMeshFreqs_$a"
        done
fi
[ "`getOption -v`" = "True" ] && echo options:$options:

# read temperature range from parameters.dat
[ "`getOption -v`" = "True" ] && echo FqhExact:$FqhExact: FqhMesh:$FqhMesh:
#echo :ex:$FqhExact:  mesh:$FqhMesh:

TRange="{1,2000,2}";[ "`isnumber.sh $tmelt`" = "yes" ] && TRange="{1,$tmelt,2}"
TRange_param=`get TRange`;[ ! -z "$TRange_param" ] && TRange=$TRange_param
#if [ "`getOption -Fe`" = "True" -o "`getOption -Fm`" = "True" -o "`getOption -A`" ]; then
#  TRange=`get TRange`;if [ -z "$TRange" ]; then error "TRange missing"; fi
#else
#  TRange="{0,0,0}"
#fi



# output information
[ "`getOption -ge`" = "True" ] && echo $element | sed 's|"||g' && exit
[ "`getOption -gc`" = "True" ] && echo $cellType | sed 's|"||g' && exit
[ "`getOption -gs`" = "True" ] && echo $supercell && exit
[ "`getOption -ga`" = "True" ] && echo $aLats && exit
[ "`getOption -gm`" = "True" ] && echo $mass && exit

[ "`getOption -gt`" = "True" ] && echo $tmelt && exit
[ "`getOption -gn`" = "True" ] && echo $atoms && exit


#######################################################################################
# get Options for mathematica (necessary before reading and checking parameters.dat 
#######################################################################################
        dynMatR=`getOption -D`; HesseMatrix=`getOption -H` exactFreqs=`getOption -fe`;
        meshFreqs=`getOption -fm`; freqs=`getOption -f`; FqhExact=`getOption -Fe`;
        FqhMesh=`getOption -Fm`; Fqh=`getOption -F`; phonons=`getOption -P`; all=`getOption -A`;
        if [ $all = True ]; then dynMatR=True; HesseMatrix=True; freqs=True; Fqh=True; phonons=True; fi
        if [ $freqs = True ]; then exactFreqs=True; meshFreqs=True; fi
        if [ $Fqh = True ]; then FqhExact=True; FqhMesh=True; fi
        flags="{$dynMatR,$HesseMatrix,$exactFreqs,$meshFreqs,$FqhExact,$FqhMesh,$phonons}"
        #echo fl$flags
        [ "`getOption -v`" = "True" ] && echo flags:$flags

####################################################################
# prepare coordinates file
####################################################################
cstr=`echo $cellType | sed 's/"//g'`; utilityPath=$path/utilities/$cstr
f=`getOption -C`
if [ $f = True ]; then
  f=`getValue -C`
  if [ "$f" = "" ]; then error "no value given to -C option"; fi
else
  f=$utilityPath/coordinates_$supercell\x$supercell\x$supercell\sc
fi
check $f; cp $f coordinates

# prepare kpMeshes
if [ $exactFreqs = True -o $FqhExact = True ]; then
  f=$utilityPath/exactRedKPoints_$supercell\x$supercell\x$supercell\sc_noGamma.dat
  check $f; cp $f exactKP
fi
if [ $meshFreqs = True -o $FqhMesh = True ]; then
  f=$utilityPath/$kpMeshFile; check $f; cp $f meshKP;
  f=$utilityPath/$kpMeshWeightFile; check $f; cp $f meshWeights;
fi

# get atomic mass
element=`echo $element | sed 's/"//g'`
check $path/getAtomicMass.sh; mass=`$path/getAtomicMass.sh $element`

# check if we are calculating classical free energy
cOp=`getOption -c`
if [ $cOp = True ]; then type="classical"; else type="quantum"; fi

# check if -afm option for running an afm spint type configuration is given
afmOp=`getOption -afm`
if [ $afmOp = True ]; then spinType=afm; else spinType=nm; fi

# otuput
if [ "`getOption -v`" = "True" ];then
    echo disp:$disp:
    echo mass:$mass:
    echo cellType:$cellType:
    echo supercell:$supercell:
    echo aLats:$aLats:
    echo TRange:$TRange:
    echo flags:$flags:
    echo type:$type:
    echo spinType:$spinType:
fi

########################################################################################
# checks before run mathematica
########################################################################################
# mathematica kernel if needed
checkAndSetMath
runmath="True";[ `echo $flags | grep True | wc -l | sed 's|[ ]*||'` = 0 ] && runmath="no"
runvak="no";[ "`getOption -vac`" = "True" ] && runvak="True"
#[ "`getOption -vac`" != "True" ] && if [ `echo $flags | grep True | wc -l | sed 's|[ ]*||'` = 0 ]; then error "unknown flag"; fi


########################################################################################
# run mathematica bulk
########################################################################################
#if [ "`getOption -vac`" != "True" ];then
if [ "$runmath" = "True" ];then
    # create forces.$a and disp.$a if not existing
    aLatsin=$aLats
    for a in `echo $aLats | sed 's|{||' | sed 's|,| |g' | sed 's|}||' | sed 's|;||'`;do           
        [ ! -e forces.$a ] && [ "`getOption -afm`" != "True" ] && echo "WARNING: forces.$a file does not exist; I will run extractForces.sh first" && extractForces.sh

        forceswords="";[ "`getOption -afm`" != "True" ] && forceswords=`wc -w forces.$a | awk '{print $1}'` 
        #echo forceswords:$forceswords:
        if [ "$forceswords" = "0" ];then
            echored "ERROR: forces.$a are empty, please check if your run finished correctly"
            [ "`echo $aLats | grep -o ",$a"`" = ",$a" ] && aLats="{`echo "$alats" | sed 's|,'"$a"'||'`}"
            [ "`echo $aLats | grep -o "$a,"`" = "$a," ] && aLats="{`echo "$alats" | sed 's|'"$a"',||'`}"
        fi
        #echo "aLatsin:$aLatsin:"
        #echo "  aLats:$aLats:"
        # check number of lines in forces
    done
    echo "aLats for Mathematica:$aLats:"
    # run mathematica
$math << EOF
<<$path/mathematica/ALL.math;
Timing[singleSpeciesQH[$mass,$cellType,$supercell,$aLats,$TRange,$flags,"$type","$spinType"]]
EOF
fi


########################################################################################
# run mathematica vacancy
########################################################################################
if [ "$runvak" = "True" ];then
for a in `echo $aLats | sed 's|{||' | sed 's|,| |g' | sed 's|}||'`;do
    [ ! -e $a\Ang ] && error "$a\Ang folder does not exist" && exit
    [ ! -e $a\Ang/allForces.dat ] && error "$a\Ang/allForces.dat file does not exist run createFolders_dynMat_vacancy.sh -g" && exit
    atoms=`cat relaxedCoords_$a | wc -l | sed 's|[ ]*||'`
    # ToDo: include check of other files which are imported further down in mathematica
done
disp=`get disp`;[ "`isnumber.sh $disp`" != "yes" ] && echo disp is missing in parameters.dat && exit
tmelt=`getMeltingPoint.sh $element -r`; [ "`isnumber.sh $tmelt`" != "yes" ] && echo Melting point unknown && exit
[ "`isnumber.sh $atoms`" != "yes" ] && echo could not determine number of atoms && exit
#echo "atoms: $atoms  <---- THIS has to be the number of atoms in the defect cell (31 for fcc vak in 2x2x2sc or 107 for fcc vac in 3x3x3sc)"
#echo disp: $disp
#echo tmelt: $tmelt
$math << EOF
<<$path/mathematica/ALL.math;
SetDirectory["`pwd`"]
aLatList = $aLats;
Print["aLats: ",aLatList];
atoms = $atoms;  (* das aendert nur FvibSupercell(_perAtom) *)
Print["ATOMS: ",atoms, " e.g. 31 for 2x2x2 fcc vak or 107 for 3x3x3 fcc vak"];
sc = $supercell; (* supercell *)
disp = $disp;
mass = $mass; (*63.546 for Cu and 26.9815386 for Al*)
type = $cellType;
TStart = 1;
TEnd = $tmelt;
TStep = 2;
Print["supercell: ",sc];
Print["disp: ",disp];
Print["mass: ",mass];
Print["type: ",type];
(* keep this: Note: this is just the usual scaling for N-1 atoms due to degrees of freedom*)
(* oder anders: Es handelt sich und die 3nuller Frequenzen *)
scale=1. atoms/(atoms-1);
Print["scale: ",scale];
eVbyAngToharbyBohr=0.0194469054353735189;
Print["get in Do loop ..."];
Do[
  aLat = aLatList[[i]];
  Print["aLat:",aLat];
(*  Print["str:  ", 
   "~/Thermodynamics/utilities/" <> type <> "/coordinates_vacancy_" <>
     ToString[sc] <> "x" <> ToString[sc] <> "x" <> ToString[sc] <> 
    "sc"];
    *)
  str = Import[
    "~/Thermodynamics/utilities/" <> type <> "/coordinates_vacancy_" <>
      ToString[sc] <> "x" <> ToString[sc] <> "x" <> ToString[sc] <> 
     "sc", "Table"];
(*Print["str: ",str];*)
  eqList = 
   Import["~/Thermodynamics/utilities/" <> type <> 
     "/eqList_vacancy_" <> ToString[sc] <> "x" <> ToString[sc] <> 
     "x" <> ToString[sc] <> "sc", "Table"];
     (*Print["eqList:",eqList];*)
  sList = 
   Partition[Partition[#, 3], 3] & /@ 
    Import["~/Thermodynamics/utilities/" <> type <> 
      "/sList_vacancy_" <> ToString[sc] <> "x" <> ToString[sc] <> 
      "x" <> ToString[sc] <> "sc", "Table"];
(*Print["sList: ",sList];*)
Print["---"];
  forces = 
   eVbyAngToharbyBohr Partition[#, 3] & /@ 
    Import[ToString[aLat] <> "Ang/allForces.dat","Table"];
(*Print["---",forces];*)
If[forces==\$Failed,{Print["ERROR: COULD NOT READ ",ToString[aLat]<>"Ang/allForces.dat"];Exit[]}];
  supercell = scCell[1. sc];
  
d=getFullDynamicalMatrix[str,supercell,forces,eqList,sList];
  Print["getFullDynamicalMatrix .... done      disp: ",disp,"        mass: ",mass];
  dynMat=d/disp/mass;
  hesse=-d/disp;
  
   Export["DynMatRSupercell_"<>ToString[aLat],dynMat,"Table",FieldSeparators->" "];

   Export["HesseMatrix_"<>ToString[aLat],hesse,"Table",FieldSeparators->" "];
   
   Print["dynmat//Dimensions ",dynMat//Dimensions];
   (*Export["hallo.dat",dynMat,"Table"];*)
   ev=Eigenvalues[-dynMat];
   Print["eigenvalues//Length ",ev//Length];
   ka=Drop[Eigenvalues[-dynMat],-3];
   freqs = dynMatToFreq*Sqrt[Drop[Eigenvalues[-dynMat],-3]];
   If[freqs // Length < 10, {Print["problems with eigenvalues"];Exit[]}];

  Export["FreqsSupercell_" <> ToString[aLat], freqs, "Table"];

  Print["FreqsSupercell_ ... exported"];
   F = getFvib[TStart, TEnd, TStep, freqs];
   Fs=Table[{F[[i]][[1]],1. F[[i]][[2]]*scale},{i,1,Length[F]}];
   Fsa=Table[{F[[i]][[1]],1. F[[i]][[2]]*scale/atoms},{i,1,Length[F]}];
   Export["FvibSupercell_" <> ToString[aLat], Append[Fs, {Null, Null}],"Table"];
   Export["FvibSupercell_perAtom_" <> ToString[aLat], Append[Fsa, {Null, Null}],"Table"];
  , {i, 1, aLatList // Length}];
EOF
fi


########################################################################################
# -fa
########################################################################################
if [ "`getOption -fa`" = "True" -o "`getOption -A`" = "True" ];then
#if [ -e mean_freqs -a $overwrite != True ]; then error "MeanFreqs existing; use -o to overwrite"; fi
rm -f mean_freqs
for a in `echo $aLats | sed 's|{||' | sed 's|,| |g' | sed 's|}||' | sed 's|;||g' `;do
    #echo a:$a
    [ ! -e "ExactFreqs_$a" ] && echored "ExactFreqs_$a does not exist, create it first using `basename $0` -fe"
    mean=`cat ExactFreqs_$a | xargs -n1 | awk '{sum1+=$1} {printf "%.10f\n", sum1/NR}' | tail -1`
    echo $a $mean >> mean_freqs
done
fi



##########################################################################################
# -makefits (checks)
##########################################################################################
if [ "`getOption -makefits`" = "True" ];then
    # mathematica kernel if needed
    checkAndSetMath
    mkfolders="fitperatom fitsupercell"
    hier=`pwd`
    for mkfolder in $mkfolders;do
#         echo "celltype:$cellType:"
         if [ "$mkfolder" = "fitperatom" ]; then
         scaling=1 && sctake=1
        [ "$cellType" = '"bcc"' ] && struct=2                 ## 4 for fcc 2 for bcc                 
        [ "$cellType" = '"fcc"' ] && struct=4
          fi
        echo "struct"$struct
        [ "$mkfolder" = "fitsupercell" ] && scaling=$atoms && sctake=$supercell && struct=1
        echo "#######################################################################"
        echo "# --> $mkfolder" # , scaling:$scaling atoms:$atoms"
        echo "#######################################################################"
#        if [ -e $mkfolder -a $overwrite != True ]; then error "$mkfolder existing; use -o to overwrite"; fi
        rm -rf $mkfolder
        mkdir $mkfolder
        cd $hier
        cd $mkfolder
        rm -f parameters.math
        echo "(* !!! mathematica format !!!  *)" > parameters.math
        echo "aLats = $aLats (* list of lattice constants; use aLats=Table[a,{a,XXX,XXX}]; or aLats={XXX,XXX,...}; *)" >> parameters.math
        #echo "latType = "$cellType";                (* bcc, fcc, hcp, or size of supercell (2,3,4,...; needed for defects) *)" >> parameters.math
        #echo "cBya = 1;                             (* 1 for bcc, fcc, supercell, and actual c/a ratio for hcp *)" >> parameters.math
        #echo "s = $scaling;                                (* scaling factor for electronic free energy; typically 1/nAtoms or *)" >> parameters.math
        #echo "                                      (* simply 1 for defect supercells *)" >> parameters.math
        #echo "Vorder = 3;                           (* basis for the volume parametrization; typically: 3 giving: {1,V,V^2,V^3} *)" >> parameters.math
        #echo "baseName = {\"Fqh_fromExactFreqs_\", \"Fqh_fromMeshFreqs_\"};    (* base names for the free energy files; followed by aLat WITHOUT UNIT (e.g. no angstrom) *)" >> parameters.math
        #echo "suffix = \"\";                          (* sufix after aLat such as e.g. angstrom or Ang or blank i.e. "" *)" >> parameters.math


        echo "structureFactor = $struct;                                          (*  4: fcc(bulk)  2: bcc(bulk)  1: supercells(defect) *)" >> parameters.math
        echo "sc = $sctake;                                                       (*  supercell if defect; 1 for bulk *)" >> parameters.math
        echo "s = $scaling;                                                 (*  scaling factor for Fqh bulk:1  defect:[number of atoms in sc]*)" >> parameters.math
        echo "fitOrder = 2;                                                 (*  typically: 2 giving: {1,V,V^2} *)" >> parameters.math
        echo "baseNames = {\"Fqh_fromExactFreqs_\",\"Fqh_fromMeshFreqs_\"}; (*  Fvib_fromExactFreqs_  *)" >> parameters.math
        echo "<<\"$path/mathematica/fitFqh.MATH\"" >> parameters.math
        echo "parameters.math created in "`pwd`

        
        ## copy baseName Folder

#        echo $aLats
        for a in `echo $aLats | sed 's|{||' | sed 's|,| |g' | sed 's|}||' | sed 's|;||'`;do      
        echo "a: $a"
        if [ -e "../flipped_Fqh/Fqh_fromExactFreqs_$a" ];then                                  #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! if schleife

        cp ../flipped_Fqh/Fqh_fromExactFreqs_$a .
            else
            cp ../Fqh_fromExactFreqs_$a .
            cp ../Fqh_fromMeshFreqs_$a .
            fi
        done
        ls
        $math < parameters.math   ## fitting the surface
        

        ############################################
        ## getThermodynamcis only for the case of peratom
        ############################################
        [ "$mkfolder" = "fitsupercell" ] && continue
        therm="Exact Mesh"
        hierm=`pwd`
        for th in $therm;do
            
            echo "################################"
            echo "# --> $mkfolder --> $th "
            echo "################################"
            cd $hierm
            rm -rf getThermodynamics_$th
            mkdir getThermodynamics_$th
            cd getThermodynamics_$th
            cp ../Fqh_from$th\Freqs_fit_order2 Fqh
            if [ -e "$hier/EVinet" ];then
                cp $hier/EVinet .
                getThermodynamics.sh
            else
                echored "DID NOT FIND EVinet and can not make getThermodynamics.sh; please add EVinet to $hier"
            fi
            cd $hierm
        done
        cd $hier
    done
    exit
fi
