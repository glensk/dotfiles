#!/bin/bash
out=no #yes #(print additional info for debugging when yes)
path=`set | grep BASH_SOURCE | sed 's|.*\"\(.*\)/[^/]*\".*|\1|' | sed '/BASH_SOURCE=/s/.*/\./'`;[ "$out" = "yes" ] && echo path: $path
script=`set | grep BASH_SOURCE | sed 's|.*\".*/\([^/]*\)\".*|\1|' | sed '/BASH_SOURCE=/s/.*/\./'`;[ "$out" = "yes" ] && echo script: $script
options=$*; . $path/../utilities/functions.include; checkOptions "-h -help -v -f -o";[ "$out" = "yes" ] && echo options: $optionsA

if [ `getOption -h` = True ] || [ `getOption -help` = True ]; then
  usage "`basename $0` can be run without parameters; this skript has to be run in the folder \"auswertung_lowplushighFahbest/Fah_fre_best_changedHesseref\""
  printOptions " " \
               "-f          create fit.input" \
               "-o          overwrite if file exists" \ 
               "-v          verbose mode"
  exit
fi

# check mathematica
checkAndSetMath

[ "`pwd | sed 's|.*/||'`" = "auswertung_lowplushighFahbest" ] && [ -e "Fah_fre_best_changedHesseref" ] && cd Fah_fre_best_changedHesseref
[ "`pwd | sed 's|.*/||'`" != "Fah_fre_best_changedHesseref" ] && echo 'you have to run this skript in the folder: "auswertung_lowplushighFahbest/Fah_fre_best_changedHesseref"' && exit


# Fah_from_fit_tangens
[ ! -e "Fah_from_fit_tangens" ] && echo 4 > Fah_from_fit_tangens

# Fah_surface
[ ! -e "Fah_surface" ] && echo Fah_surface does not exist && exit

# mean.freqs
lowpath=`grep "lowpath=" ../../parameters.dat | awk '{print $1}' | sed 's|lowpath=||'`
if [ ! -e "$lowpath" ];then
   hier=`pwd`
  cd ../..
  lowpath=`ti_high_0_create_Folders_vasp.sh -l`
  cd $hier
  [ ! -e "$lowpath" ] && echo lowpath: $lowpath : not found1 && exit  
fi
[ "`getOption -v`" = "True" ] && echo lowpath: $lowpath

[ ! -e "$lowpath/parameters.dat" ] && echo $lowpath/parameters.dat not found2 && exit
hessepath=`grep "hessepath=" $lowpath/parameters.dat | awk '{print $1}' | sed 's|hessepath=||'`
if [ ! -e "$hessepath" ];then
    hier=`pwd`
    cd $lowpath
    hessepath=`ti_low_0_create_Folders_vasp.sh -e`
    cd $hier
    [ ! -e "$hessepath" ] && echo hessepath: $hessepath : not found2 && exit
fi

[ "`getOption -v`" = "True" ] && echo hessepath: $hessepath

alats=`cat Fah_surface | awk '{print $2}' | sort | uniq | xargs`
[ "`getOption -v`" = "True" ] && echo alats: $alats

if [ ! -e "mean_freqs" ];then
for a in $alats;do
    [ "`getOption -v`" = "True" ] && echo a $a
    hesse=$hessepath/HesseMatrix_$a
    [ ! -e "$hesse" ] && echo $hesse does not exisst && exit
    freq=`dynmat.py -htmf $hesse`
    [ "`isnumber.sh $freq`" != "yes" ] && echo couldnt find freq from $hesse :$freq: && exit
    [ "`getOption -v`" = "True" ] && echo a $a   freq:$freq
    echo $a $freq >> mean_freq_tmp
done
mv mean_freq_tmp mean_freqs
fi

# fit.surface
amin=`cat Fah_surface | awk '{print $2}' | sort -n | head -1`
amax=`cat Fah_surface | awk '{print $2}' | sort -n | tail -1`
[ "`getOption -v`" = "True" ] && echo amin:$amin amax:$amax

[ ! -e "$hessepath/POSCAR_$amax" ] && echo $hessepath/POSCAR_$amax not found3 && exit
numatoms=`POSCAR_numatoms.sh $hessepath/POSCAR_$amax`
[ "`getOption -v`" = "True" ] && echo numatoms:$numatoms:

sc=""
[ "$numatoms" = "31" ] && sc=2
[ "$numatoms" = "107" ] && sc=3
[ "$sc" = "" ] && echo add data for sc && exit


element=""; [ -e "$hessepath/POTCAR" ] && element=`POTCAR_element.sh $hessepath/POTCAR`
[ "$element" = "" ] && tmelt=2000
[ "$element" != "" ] && tmelt=`getMeltingPoint.sh $element -r`
[ "`getOption -v`" = "True" ] && echo tmelt:$tmelt


if [ "`getOption -f`" = "True" ];then
    overwrite=`getOption -o`
    if [ -e fit.input -a $overwrite != True ]; then error "parameters.dat existing; use -o to overwrite"; fi
    rm -f fit.input
echo "(* adjustable parameters start *)" > fit.input
echo "" >> fit.input
echo "FsurfFile = \"Fah_surface\";" >> fit.input
echo "type = 1;                                                  (*  1: T(K)  aLat(Ang/at)  F(meV/at)   *)" >> fit.input
echo "                                                           (*  2: T(K)  V(Ang/at)     F(meV/at)   *)" >> fit.input
echo "                                                           (*  3: T(K)  V(Ang/cell)   F(meV/cell) *)" >> fit.input
echo "" >> fit.input
echo "min = $amin;                                                (*  aLat or volume range (same format as Fsurf) for the 2. fit *)" >> fit.input
echo "max = $amax;                                                (*  typically: Vmin=Veq(T=0K) and Vmax=Veq(Tmelt) *)" >> fit.input
echo "mesh = 100;                                                (* 100 is good and should be kept *)" >> fit.input
echo "" >> fit.input
echo "structureFactor = 1;                                       (*  4: fcc  2: bcc  1: supercells *)" >> fit.input
echo "sc = $sc;                                                    (*  supercell, 1 for bulk *)" >> fit.input
echo "nAtoms = $numatoms;                                               (*  1 for bulk *)" >> fit.input
echo "" >> fit.input
echo "fitType = \"Fvib\";                                          (*  "Fvib"  or  "poly"  fit type for 1. fit; take "Fvib" for Fah or Fel *)" >> fit.input
echo "basis[V_, T_] := {1,T, V }                                 (*  for Fah typically: "Fvib" and {1,T,V} *)" >> fit.input
echo "                                                           (*  for Fel typically: "Fvib" and {1,T, V,T V,T^2,V^2,V^3} *)" >> fit.input
echo "basis2[V_]:={1, V, V^2, V^3}                               (*  should be more than sufficient: {1,V,V^2,V^3} *)" >> fit.input
echo "" >> fit.input
echo "minT = 1;                                                  (* typically 1     *)" >> fit.input
echo "maxT = $tmelt;                                               (*           Tmelt *)" >> fit.input
echo "stepT = 1;                                                 (*           2     *)" >> fit.input
echo "" >> fit.input
echo "useMeanFreqs=True;                                         (* if True "mean_freqs" file must be available; format as Fsurf, e.g. aLat(Ang) meanFreq(meV) *)" >> fit.input
echo "                                                           (* meanFreqs are then used in the fit formula (check fitSurface.math) *)" >> fit.input
echo "(* adjustable parameters end *)" >> fit.input
echo "" >> fit.input
echo "<<"$path/../mathematica/fitSurface.MATH" >> fit.input" >> fit.input

fi    


# run mathematica
math < fit.input
