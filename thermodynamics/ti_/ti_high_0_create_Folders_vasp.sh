#!/bin/bash

# up to creating UPTILD folder everything is mac ready; but after that .... who knows .... not done yet
debug=no # yes (print additional info for debugging when yes)
out=no #yes #(print additional info for debugging when yes)
path=`set | grep BASH_SOURCE | sed 's|.*\"\(.*\)/[^/]*\".*|\1|' | sed '/BASH_SOURCE=/s/.*/\./'`;[ "$out" = "yes" ] && echo path: $path
script=`set | grep BASH_SOURCE | sed 's|.*\".*/\([^/]*\)\".*|\1|' | sed '/BASH_SOURCE=/s/.*/\./'`;[ "$out" = "yes" ] && echo script: $script
options=$*; . $path/../utilities/functions.include; checkOptions "-h -help -p -i -k -o -a -t -c -r -d -l -e -ua -ut -ul -un -Hp -He -Hc -Hs -Hn -Hbd";[ "$out" = "yes" ] && echo options: $options
lom=`getlinuxormac`;[ "$lom" = "Linux" ] && add="";[ "$lom" = "mac" ] && add="'' -e"


if [ `getOption -h` = True ] || [ `getOption -help` = True ] || [ $# = 0 ]; then
  usage `basename $0`
  printOptions " " \
               "-p          create parameters.dat template and exit" \
               "-i          create INCAR template and exit" \
               "-k          create KPOINTS template and exit" \
               "-t          get POTCAR from lowrun if possible and exit" \
               "-o          force overwrite of file if it exists" \
               "-a          all above and exit (use this with -o to overwrite existing files)" \
               "-ua [Real]  use this lattice constant instead of the one from parameters.dat" \
               "-ut [Real]  use this temperature      instead of the one from parameters.dat" \
               "-ul [Real]  use this lambda           instead of the one from parameters.dat" \
               "-un [INT]   use this n to force additional job creation (`basenaem $0` * -m 0.05)" \
               "-c          create folder" \
               " " \
               "-l          echo lowpath and exit" \
               "-Hp         echo hessepath and exit" \
               "-He         echo element (from hessepath) and exit " \
               "-Hc         echo cellType (from hessepath) and exit " \
               "-Hs         echo supercell (from hessepath) and exit " \
               "-Hn         echo number of atoms (from hessepath) and exit" \
               "-Hbd        echo bulk or def (bulk/defectcell) and exit" \
               "-r          echo refpath ( e.g. lowpath/ref_high_2x2x2sc )" \
               "-e          echo refpathjob ( e.g. lowpath/ref_high_2x2x2sc/400eV_8x8x8kp )" \
               "-d          debug mode for plotting inbetween information"
  exit
fi

# check if at least one option
[ `getOption -d` = True ] && debug=yes

#################################################
# xxx get data from pfad
#################################################
## stuff always to compare to parameters.dat
[ "$debug" = "yes" ] && echo xxx get data from pfad
hier=`pwd -L`
cutoff_pfad=`echo $hier | sed 's|.*[-_]\([0-9]*\)eV.*|\1|'`; [ "`isnumber $cutoff_pfad`" != "yes" ] && cutoff_pfad="???"
kpoints_pfad=`echo $hier | grep -o "[0-9]*[xmg][0-9]*[xmg][0-9]*kp\|kp[0-9]*[xmg][0-9]*[xmg][0-9]*" | tail -1 | sed 's|[xmkp]| |g' | sed 's|^[0]*||' | sed 's| 0| |g' | sed -e 's/^[ \t]*//' | sed 's/^[ \t]*//;s/[ \t]*$//'`;[ "`isnumber $kpoints_pfad`" != "yes" ] && kpoints_pfad="???"


# stuff just to get if parameters.dat is written
#genPar=`getOption -p`
#allop=`getOption -a`
#overwrite=`getOption -o`
#if [ $genPar = True ] || [ $allop = True ]; then
hier=`pwd -L`
lowpath_pfad=`pwd -L | sed 's|\(.*\)\([__]*\)high_.*|\1|' | sed 's|__$||'`; [ ! -e "$lowpath_pfad" ] && lowpath_pfad="???"
alats_folder="???"; [ -e "$lowpath_pfad" ] && alats_folder=`ls -1d $lowpath_pfad/[0-9.]*Ang_[0-9.]*K | sed 's|'"$lowpath_pfad"'/\([0-9.]*\)Ang_\([0-9.]*\)K|\1|' | sort -n | uniq | xargs`
temps_folder="???"; [ -e "$lowpath_pfad" ] && temps_folder=`ls -1d $lowpath_pfad/[0-9.]*Ang_[0-9.]*K | sed 's|'"$lowpath_pfad"'/\([0-9.]*\)Ang_\([0-9.]*\)K|\2|' | sort -n | uniq | xargs`
nbands_pfad="???";
ngxf_pfad="none";
ngyf_pfad="none";
ngzf_pfad="none";

[ "$debug" = "yes" ] && echo xx
if [ -e "$lowpath_pfad" ];then
  cd $lowpath_pfad

  [ "$debug" = "yes" ] && echo pwd:`pwd`
  [ ! -e "parameters.dat" ] && echo parameters.dat not existing in `pwd` please make sure it does && exit
  nbands_pfad=`get NBANDS`
  ngxf_pfad=`get NGXF`
  ngyf_pfad=`get NGYF`
  ngzf_pfad=`get NGZF`
  
  [ "$debug" = "yes" ] && echo nbands_pfad:$nbands_pfad
  if [ "$nbands_pfad" = "" ];then
        incar=`find -L . -mindepth 3 -maxdepth 3 -name INCAR -print -quit`
        [ "$debug" = "yes" ] && echo incar: $incar
        [ "$incar" != "" ] && nbands_pfad=`INCAR_NBANDS.sh $incar`
        fi
  cd $hier
  fi
[ "$ngxf_pfad" = "" ] && ngxf_pfad="none";
[ "$ngyf_pfad" = "" ] && ngyf_pfad="none";
[ "$ngzf_pfad" = "" ] && ngzf_pfad="none";

[ "$debug" = "yes" ] && echo YY
#    && INCARS=`find -L $lowpath_pfad/*K/lambda*/ -maxdepth 2 -mindepth 1 -type f -name INCAR`
#nbands_pfad="???"; [ -e "$lowpath_pfad" ] && nbands_pfad=`sed -n '/NBANDS/p' $INCARS | awk '{print $3}' | sort | uniq | xargs` && [ "`isnumber $nbands_pfad`" != "yes" ] && nbands_pfad="???"
#ngxf_pfad="none"; [ -e "$lowpath_pfad" ] && ngxf_pfad=`sed -n '/NGXF/p' $INCARS | awk '{print $3}' | sort | uniq | xargs` && [ "`isnumber $ngxf_pfad`" != "yes" ] && ngxf_pfad="none"
#ngyf_pfad="none"; [ -e "$lowpath_pfad" ] && ngyf_pfad=`sed -n '/NGYF/p' $INCARS | awk '{print $3}' | sort | uniq | xargs` && [ "`isnumber $ngyf_pfad`" != "yes" ] && ngyf_pfad="none"
#ngzf_pfad="none"; [ -e "$lowpath_pfad" ] && ngzf_pfad=`sed -n '/NGZF/p' $INCARS | awk '{print $3}' | sort | uniq | xargs` && [ "`isnumber $ngzf_pfad`" != "yes" ] && ngzf_pfad="none"
#lambdas_pfad="???"; [ -e "$lowpath_pfad" ] && lambdas_pfad=`ls -1d $lowpath_pfad/*K/lambda* | sed 's|.*lambda\([0-9.]*\)_.*|\1|' | sort | uniq | xargs`
lambdas_pfad="0.0 0.15 0.5 0.85 1.0" # everything else takes to longe ... possibly take also from parameters.dat.. 

### HESSEMAT_pfad and refpath_pfad
[ "$debug" = "yes" ] && echo lowpath_pfad: $lowpath_pfad
subproject_tmp=""; [ -e "$lowpath_pfad" ] && subproject_tmp=`echo $lowpath_pfad | sed 's|\(.*\)/.*|\1|'`   ## only necessary to find -L refpath this is the
[ "$debug" = "yes" ] && echo subproject_tmp:$subproject_tmp
sc_string_tmp=`pwd | grep -o "[0-9]*x[0-9]*x[0-9]*sc"`
#hessefolder="$subproject_tmp/Hessematrix_$sc_string_tmp"


hiertmp=`pwd`
hessefolder=`cd $lowpath_pfad; grep "hessepath.*=" parameters.dat | sed 's|hessepath[^=]*=||' | sed 's| ||g'`

[ ! -e "$hessefolder" ] && hessefolder=???
refpath_pfad=$subproject_tmp/ref_high_$sc_string_tmp

        if [ "$debug" = "yes" ];then
            echo "#########################################"
            echo "## get data from path "
            echo "########################################"
            echo lowpath_pfad:$lowpath_pfad
            echo cut:$cutoff_pfad
            echo kp:$kpoints_pfad
            echo alats_folder:$alats_folder
            echo temps_folder:$temps_folder
            echo nbands_pfad:$nbands_pfad
            echo lam:$lambdas_pfad
            echo ngxf:$ngxf_pfad
            echo ngyf:$ngyf_pfad
            echo ngzf:$ngzf_pfad
            echo hessefolder:$hessefolder
            fi

[ "`getOption -Hp`" = "True" ] && echo $hessefolder && exit
# echo lowpath
if [ "`getOption -Hp`" = "True" ] || [ "`getOption -Hc`" = "True" ] || [ "`getOption -He`" = "True" ] || [ "`getOption -Hs`" = "True" ] || [ "`getOption -Hn`" = "True" ] || [ "`getOption -Hbd`" = "True" ];then
    #lowpath_out=`get lowpath`
    #echo ccc

    [ "$debug" = "yes" ] && echo in
    #lowpath_parameters=$lowpath_out
    #echo aa
    #[ ! -e "$lowpath_parameters" ] && lowpath_out=$lowpath_pfad
    #[ ! -e "$lowpath_pfad" ] && echo lowpath not found : $lowpath_parameters : $lowpath_out && exit
    #[ "`getOption -l`" = "True" ] && echo $lowpath_out && exit

    # if we have -Hp -Hc -He -Hs we need HessePath
    #Hp=`cd $lowpath_out && ti_low_0_create_Folders_vasp.sh -e`
    Hp=$hessefolder
    [ "$debug" = "yes" ] && echo Hp:$hessefolder

    He=`getSingleSpeciesPhonons.sh -folder $Hp -ge`
    [ "`getOption -He`" = "True" ] && echo $He && exit

    [ "$debug" = "yes" ] && echo he:$He
    Hc=`getSingleSpeciesPhonons.sh -folder $Hp -gc`
    [ "`getOption -Hc`" = "True" ] && echo $Hc && exit

    [ "$debug" = "yes" ] && echo hc:$Hc
    Hs=`getSingleSpeciesPhonons.sh -folder $Hp -gs`
    [ "`getOption -Hs`" = "True" ] && echo $Hs && exit
    [ "$debug" = "yes" ] && echo hs:$Hs

    Hn=`getSingleSpeciesPhonons.sh -folder $Hp -gn`
    [ "`getOption -Hn`" = "True" ] && echo $Hn && exit

    [ "$debug" = "yes" ] && echo hn:$Hn
    Hbd=`getSingleSpeciesPhonons.sh -folder $Hp -gbd`
    [ "`getOption -Hbd`" = "True" ] && echo $Hbd && exit
    [ "$debug" = "yes" ] && echo hbd:$Hbd
fi

if [ "`getOption -e`" = "True" ];then
    refpath_out=`get refpath`
    refpath_parameters=$refpath_out
    [ ! -e "$refpath_parameters" ] && refpath=$refpath_pfad
    [ -e "$refpath_parameters" ] && refpath=$refpath_parameters
        kpoints=`get kp | sed 's/^[ \t]*//;s/[ \t]*$//'`
        [ "`isnumber $kpoints`" != "yes" ] && kpoints=`echo $hier | grep -o "[0-9]*[xmg][0-9]*[xmg][0-9]*kp\|kp[0-9]*[xmg][0-9]*[xmg][0-9]*" | tail -1 | sed 's|[xmkp]| |g' | sed 's|^[0]*||' | sed 's| 0| |g' | sed -e 's/^[ \t]*//' | sed 's/^[ \t]*//;s/[ \t]*$//'`;[ "`isnumber $kpoints_pfad`" != "yes" ] && kpoints_pfad="???"
        [ "`isnumber.sh $kpoints`" != "yes" ] && echo could not find refpath error1 && exit
        kpstring=`echo $kpoints | sed 's| |x|g'`
        kp1=`echo $kpoints | awk '{print $1}'`
        kp2=`echo $kpoints | awk '{print $2}'`
        kp3=`echo $kpoints | awk '{print $3}'`
        # cutoff
        cutoff=`get cutoff`
        [ "`isnumber $cutoff`" != "yes" ] && cutoff=`echo $hier | sed 's|.*[-_]\([0-9]*\)eV.*|\1|'`
        [ "`isnumber $cutoff`" != "yes" ] && echo could not find ferfpath error2 && exit
        [ "`getOption -d`" = "True" ] && echo ----------------- $kpstring ----------------cutoff:$cutoff------------------refpath:$refpath-----------
        refpathjobseV=`find -L $refpath/ -maxdepth 1 -mindepth 1 -type d -name '*'$cutoff\eV'*'`
        
        [ "$debug" = "yes" ] && echo refpathjobseV:$refpathjobseV
        refpathjobs=`echo "$refpathjobseV" | grep '.*[0]*'"$kp1"'x[0]*'"$kp2"'x[0]*'"$kp3"'[kK][pP].*'`
        [ "$debug" = "yes" ] && echo refpathjobs:$refpathjobs
## h    ere we create the refpathjob if it did not exist
        [ "`echo $refpathjobs | wc -w | sed 's|[ ]*||g'`" = "0" ] || [ "$refpathjobs" = "" ] && refpathjobs=$refpath/$cutoff\eV_$kpstring\kp && echored "refpathjob did not exist and was created : $refpathjobs"
        if [ "`echo $refpathjobs | wc -w | sed 's|[ ]*||g'`" -ge "2" ];then
            echo more than one refpathjobs found:
            echo $refpathjobs | xargs -n1
            echo ""
            echo not implemented yet
            exit
        fi
        [ "$debug" = "yes" ] && echo refpathjobs:$refpathjobs
         
        ## get refpathjob
        refpathjob=""; [ "`echo $refpathjobs | wc -w | sed 's|[ ]*||g'`" = "1" ] && refpathjob=$refpathjobs
        [ "$debug" = "yes" ] && echo refpathjob:$refpathjob and ls && ls $refpathjob

        [ ! -e "$refpathjob" ] && echored "refpathjob problem here is should be known $refpathjob"  && exit
        echo $refpathjob 
        exit
fi

##########################################################
# xxx create parameters.dat template and exit
##########################################################
if [ "$debug" = "yes" ];then
    echo "#######################"
    echo wring parameters.dat
    echo "#######################"
fi
genPar=`getOption -p`
allop=`getOption -a`
overwrite=`getOption -o`
if [ $genPar = True ] || [ $allop = True ]; then
    if [ -e parameters.dat -a $overwrite != True ]; then error "parameters.dat existing; use -o to overwrite"; fi
[ "`echo $nbands_pfad | wc -w | sed 's|[ ]*||g'`" != "1" ] && nbandstmp=$nbands_pfad && nbands_pfad=`echo $nbands_pfad | xargs -n1 | sort | tail -1` && echored "- TAKE CARE: NBANDS where in low_runs: $nbandstmp ; set to $nbands_pfad"
[ "???" = "$nbands_pfad" ] && incarjob_tmp=`find -L *Ang_*K/lambda*/ -name INCAR -print -quit` && bands_tmp=`grep NBANDS $incarjob_tmp | awk '{print $3+5}'` && [ "`isnumber $bands_tmp`" = "yes" ] && nbands_pfad=$bands_tmp
sepa=`echo $alats_folder | wc -c | sed 's|[ ]*||g'`
sept=`echo $temps_folder | wc -c | sed 's|[ ]*||g'`
sepl=`echo $lambdas_pfad | wc -c | sed 's|[ ]*||g'`
sepc=`echo $cutoff_pfad | wc -c | sed 's|[ ]*||g'`;
sepk=`echo $kpoints_pfad | wc -c | sed 's|[ ]*||g'`
sepb=`echo $nbands_pfad | wc -c | sed 's|[ ]*||g'`
konst1=`echo "$sepa $sept $sepl $sepc $sepk $sepb" | xargs -n1 | sort -n | tail -1`
konst=`expr $konst1 + 11 `
    str="                                                                                                   "
    sepa=`echo $alats_folder | wc -c | sed 's|[ ]*||g'`;sepa=`expr $konst - $sepa - 6 `;sepa=`echo "$str" | cut -c 1-$sepa`;
    sept=`echo $temps_folder | wc -c | sed 's|[ ]*||g'`;sept=`expr $konst - $sept - 6 `;sept=`echo "$str" | cut -c 1-$sept`
    sepl=`echo $lambdas_pfad | wc -c | sed 's|[ ]*||g'`;sepl=`expr $konst - $sepl - 8`;sepl=`echo "$str" | cut -c 1-$sepl`
    sepc=`echo $cutoff_pfad | wc -c | sed 's|[ ]*||g'`;sepc=`expr $konst - $sepc - 7`;sepc=`echo "$str" | cut -c 1-$sepc`
    sepk=`echo $kpoints_pfad | wc -c | sed 's|[ ]*||g'`;sepk=`expr $konst - $sepk - 3`;sepk=`echo "$str" | cut -c 1-$sepk`
    sepn=`echo 10 | wc -c | sed 's|[ ]*||g'`;sepn=`expr $konst - $sepn - 2`;sepn=`echo "$str" | cut -c 1-$sepn`
    sepb=`echo $nbands_pfad | wc -c | sed 's|[ ]*||g'`;sepb=`expr $konst - $sepb - 7`;sepb=`echo "$str" | cut -c 1-$sepb`
    if [ "$debug" = "yes" ];then
        echo a:"$sepa"$alats_folder:
        echo t:"$sept"$temps_folder:
        echo l:"$sepl"$lambdas_pfad:
        echo c:"$sepc"$cutoff_pfad:
        echo k:"$sepk"$kpoints_pfad:
        fi
    echo "aLats=$alats_folder `echo "$sepa"`# lattice constants: e.g. aLats=3.74 3.62 3.68"   > parameters.dat
    echo "temps=$temps_folder `echo "$sept"`# temperatures: e.g. temps=450 250 800 1100 1360" >> parameters.dat 
    echo "lambdas=$lambdas_pfad `echo "$sepl"`# lambdas: e.g. lambdas=0.0 0.5 1.0" >> parameters.dat  
    echo "cutoff=$cutoff_pfad `echo "$sepc"`# cutoff in eV: e.g. 500" >> parameters.dat
    echo "kp=$kpoints_pfad `echo "$sepk"`# kpoints: e.g. 4 4 4" >> parameters.dat
    echo "n=10 `echo "$sepn"`# number of snapshots for UPTILD per lambda;either one value for all temps or" >> parameters.dat
    echo "     `echo "$sepn"`# one per temp; 1st,2nd,3rd n is mapped to 1st,2nd,3rd temp " >> parameters.dat
    echo "NBANDS=$nbands_pfad `echo "$sepb"`# NBANDS for structure" >> parameters.dat
    echo "NGXF=$ngxf_pfad" >> parameters.dat
    echo "NGYF=$ngyf_pfad" >> parameters.dat
    echo "NGZF=$ngzf_pfad" >> parameters.dat
    echo "lowpath=$lowpath_pfad  # path to low run folder" >> parameters.dat
    echo "refpath=$refpath_pfad  # path to reference folder for High run" >> parameters.dat
    echored "- TAKE CARE: n (n=numbers of molecular dynamics snapshots for the UPTILD procedure) was set to 10; pleas change n to the amount of spes which are necessary to converge your UPTILD; 10 should usually be good."
    echo; echo "parameters.dat written";
    [ "`getOption -a`" != "True" ] && exit
    fi

    
##################################################################
# xxx read element, cellType, and supercell size from parameters.dat
##################################################################
if [ "$debug" = "yes" ];then
    echo "#######################"
    echo read from parameters.dat
    echo "#######################"
fi
check parameters.dat
[ "$debug" = "yes" ] && echo ok112
aLats=`get aLats`; temps=`get temps`; lambdas=`get lambdas`;cutoff=`get cutoff`;nbands=`get NBANDS`
n=`get n`;kpoints=`get kp | sed 's/^[ \t]*//;s/[ \t]*$//'`;lowpath=`get lowpath`; NGXF=`get NGXF`; NGYF=`get NGYF`; NGZF=`get NGZF`
refpath=`get refpath`

##################################################################
# xxx chenge values if -ua -ut -ul -un given
##################################################################
[ "`getOption -ua`" = "True" ] && aLats=`getValue -ua`
[ "`getOption -ut`" = "True" ] && temps=`getValue -ut`
[ "`getOption -ul`" = "True" ] && lambdas=`getValue -ul`
[ "`getOption -un`" = "True" ] && n=`getValue -un`

for i in "$aLats" "$temps" "$lambdas" "$cutoff" "$nbands" "$n" "$kpoints";do
    [ -z "$i" ] && echored parameters.dat needs to contain: aLats, temps, cutoff, nbands, n but something is missing && exit
    done
checkInput "$aLats" "$temps" "$lambdas" "$cutoff" "$nbands" "$n" "$kpoints"

if [ "$debug" = "yes" ];then
        echo ""
        echo READING FROM parameters.dat
        echo ""
        echo aLats:$aLats
        echo temps:$temps
        echo lambdas:$lambdas
        echo cutoff:$cutoff
        echo kpoints:$kpoints
        echo n:$n
        echo lowpath: $lowpath
        echo nbands:$nbands
        echo NGXF:$NGXF
        echo NGYF:$NGYF
        echo NGZF:$NGZF
        echo refpath: $refpath
        kkk=($temps)
        xxx=($n)
        len=`echo ${#kkk[*]}`
        for i in `awk 'BEGIN {x=-1; while(++x<='"$len"'-1){print x; }; exit}'`;do
            echo "$i "Temp:"${kkk[$i]}" K,        n:"${xxx[$i]}"
        done
        fi
[ "`isnumber $aLats`" != "yes" ] && echored "ERROR: aLats:$aLats: in parameters.dat is no number" && exit
[ "`isnumber $temps`" != "yes" ] && echored ERROR: temps:$temps: in parameters.dat is no number && exit
[ "`isnumber $lambdas`" != "yes" ] && echored ERROR: lambdas:$lambdas: in parameters.dat is no number && exit
[ "`isnumber $cutoff`" != "yes" ] && echored ERROR: cutoff:$cutoff: in parameters.dat is no number && exit
[ "`isnumber $kpoints`" != "yes" ] && echored ERROR: kpoints:$kpoints: in parameters.dat is no number && exit
[ "`echo $kpoints | wc -w | sed 's|[ ]*||g'`" != "3" ] && echored ERROR: kpoints:$kpoints: in parameters.dat needs to have 3 INTEGERS && exit
[ "`isnumber $n`" != "yes" ] && echored ERROR: n:$n: in parameters.dat is no number && exit
[ "`isnumber $nbands`" != "yes" ] && echored ERROR: nbands:$nbands: in parameters.dat is no number && exit
[ "$NGXF" != "none" ] && [ "`isnumber $NGXF`" != "yes" ] && echored ERROR: NGXF:$NGXF: in parameters.dat is no number && exit
[ "$NGYF" != "none" ] && [ "`isnumber $NGYF`" != "yes" ] && echored ERROR: NGYF:$NGYF: in parameters.dat is no number && exit
[ "$NGZF" != "none" ] && [ "`isnumber $NGZF`" != "yes" ] && echored ERROR: NGZF:$NGZF: in parameters.dat is no number && exit
[ ! -e "$lowpath" ] && echored ERROR: lowpath: $lowpath : in parameters.dat does not exist a && exit
[ ! -w "$lowpath" ] && echored ERROR: lowpath: $lowpath : in parameters.dat now write permissions && exit
[ `getOption -l` = True ] && echo $lowpath && exit
[ `getOption -r` = True ] && echo $refpath && exit
#n has to be either one number (fall all temperatures) or there have to be as many temperatures as n values
number_temps=`echo $temps | wc -w | sed 's|[ ]*||'`
number_ns=`echo $n | wc -w | sed 's|[ ]*||'`
if [ "$number_ns" != "1" ];then
    [ "$number_ns" != "$number_temps" ] && echo "(n:$n: number_ns:$number_ns: number_temps:$number_temps:) n has to be either one number (fall all temperatures) or there have to be as many temperatures as n values" && exit
fi
#hessefolder=/nas/glensk/v/pp/binary_si-mg/MgXXSi1-hcp-PAW_PBE_questMurn/Mg35Si-hcp3x3x2/ti_dilute_si1mg35_hcp3x3x2/HesseMatrix_3x3x2
[ ! -e "$hessefolder" ] && echo "no Hessematrix_...: $hessefolder : folder found: $hessefolder THIS IS CRUCIAL FOR THIS SCRIPT and has to exist" && exit
for a in $aLats;do
    [ ! -e "$hessefolder/POSCAR_$a" ] && echo $hessefolder/POSCAR_$a does not exist b && exit
done
[ "$debug" = "yes" ] && echo read in parameters.dat

    #[ ! -e "$refpath" ] && echo refpath does not exit try: `basename $0` -c to create it && exit #mkdir -p $refpath && echo $refpath created
    [ ! -e "$refpath" ] && echored "refpath created since it did not exist before!: $refpath check the refpath is ok, otherwise change folder location in parameters.dat and delete or move created folder" &&  mkdir -p $refpath 

## get refpathjobs (or make if does not existing) 
    kpstring=`echo $kpoints | sed 's| |x|g'`
    kp1=`echo $kpoints | awk '{print $1}'`
    kp2=`echo $kpoints | awk '{print $2}'`
    kp3=`echo $kpoints | awk '{print $3}'`
    
    [ "$debug" = "yes" ] && echo kpstring:$kpstring:
    [ "$debug" = "yes" ] && echo kp1:$kp1: kp2:$kp2: kp3:$kp3:

    refpathjobseV=`find -L $refpath/ -maxdepth 1 -mindepth 1 -type d -name '*'$cutoff\eV'*'`
    
    [ "$debug" = "yes" ] && echo refpathjobseV:$refpathjobseV
    refpathjobs=`echo "$refpathjobseV" | grep '.*[0]*'"$kp1"'x[0]*'"$kp2"'x[0]*'"$kp3"'[kK][pP].*'`
    [ "$debug" = "yes" ] && echo refpathjobs:$refpathjobs
## here we create the refpathjob if it did not exist
    [ "`echo $refpathjobs | wc -w | sed 's|[ ]*||g'`" = "0" ] || [ "$refpathjobs" = "" ] && mkdir -p $refpath/$cutoff\eV_$kpstring\kp && refpathjobs=$refpath/$cutoff\eV_$kpstring\kp && echored "refpathjob did not exist and was created : $refpathjobs"
    if [ "`echo $refpathjobs | wc -w | sed 's|[ ]*||g'`" -ge "2" ];then
        echo more than one refpathjobs found:
        echo $refpathjobs | xargs -n1
        echo ""
        echo not implemented yet
        exit
    fi
    [ "$debug" = "yes" ] && echo refpathjobs:$refpathjobs
     
    ## get refpathjob
    refpathjob=""; [ "`echo $refpathjobs | wc -w | sed 's|[ ]*||g'`" = "1" ] && refpathjob=$refpathjobs
    [ "$debug" = "yes" ] && echo refpathjob:$refpathjob and ls && ls $refpathjob

    [ ! -e "$refpathjob" ] && echored "refpathjob problem here is should be known $refpathjob"  && exit
    [ `getOption -e` = True ] && echo $refpathjob && exit
    [ "$debug" = "yes" ] && echo rrr:$refpathjob

#  check for consistency
[ "$cutoff_pfad" != "???" ] && [ "$cutoff_pfad" != "$cutoff" ] && echored "cutoff from path: $cutoff_pfad ;; paramerts.dat: $cutoff" && exit
[ "$kpoints_pfad" != "???" ] && [ "$kpoints_pfad" != "$kpoints" ] && echored "kpoints from path:$kpoints_pfad:  from parameters.dat:$kpoints:" && exit




[ "$debug" = "yes" ] && echo abc

########################################################################
# xxx create INCAR template and exit
########################################################################
genPar=`getOption -i`
overwrite=`getOption -o`
if [ $genPar = True ] || [ $allop = True ]; then


    if [ -e INCAR -a $overwrite != True ]; then error "INCAR existing; use -o to overwrite"; fi
    echo "NPAR=1" > INCAR
    echo "" >> INCAR
    echo "ENCUT  =    xxxCUTOFFxxx" >> INCAR
    echo "ISMEAR =     -1" >> INCAR
    echo "SIGMA  =    0.1" >> INCAR
    echo "EDIFF  =   1E-3" >> INCAR
    echo "" >> INCAR
    echo "ISYM   =      0" >> INCAR
    echo "ADDGRID=    TRUE" >> INCAR
    echo "PREC   =  Accurate" >> INCAR
    echo "LREAL =.FALSE." >> INCAR
    echo "NELMDL = -5" >> INCAR
    echo "NBANDS =    xxxNBANDSxxx" >> INCAR
    echo "ALGO   =    FAST" >> INCAR
    echo "" >> INCAR
    echo "LWAVE  =      F" >> INCAR
    echo "LCHARG =      F" >> INCAR
    echo "LVTOT = F" >> INCAR
    echo "LELF = F" >> INCAR
    echo "NWRITE = 0" >> INCAR

    [ "$NGXF" != "none" ] && echo "NGXF=$NGXF" >> INCAR
    [ "$NGYF" != "none" ] && echo "NGYF=$NGYF" >> INCAR
    [ "$NGZF" != "none" ] && echo "NGZF=$NGZF" >> INCAR
    echo "INCAR written"; 
    [ "`getOption -a`" != "True" ] && exit
    fi

########################################################################
# xxx create KPOINTS template and exit
########################################################################
genPar=`getOption -k`
overwrite=`getOption -o`
if [ "$genPar" = "True" ] || [ "$allop" = "True" ]; then
    if [ -e KPOINTS -a $overwrite != True ]; then error "KPOINTS existing; use -o to overwrite"; fi
    path=`find -L $lowpath -type f -name KPOINTS -print -quit`
    monkhorst_gamma=`sed -n '3p' $path`
    if [ "`echo "$monkhorst_gamma" | cut -c 1`" != "M" ];then
        if [ "`echo "$monkhorst_gamma" | cut -c 1`" != "m" ];then
            if [ "`echo "$monkhorst_gamma" | cut -c 1`" != "G" ];then
                if [ "`echo "$monkhorst_gamma" | cut -c 1`" != "g" ];then
                    echo "KPOINTS file: $path : neither Gamma nor Monkhorst"
                    exit
        fi
        fi
        fi
        fi
    echo "K-Points" > KPOINTS 
    echo " 0" >> KPOINTS
    echo "$monkhorst_gamma" >> KPOINTS
    echo "xxxkp1xxx xxxkp2xxx xxxkp3xxx" >> KPOINTS
    echo "0 0 0" >> KPOINTS
    echo "KPOINTS written"; 
    [ "`getOption -a`" != "True" ] && exit
    fi



########################################################################
# xxx get POTCAR from lowrun if possible 
########################################################################
genPar=`getOption -t`
overwrite=`getOption -o`
if [ "$genPar" = "True" ] || [ "$allop" = "True" ]; then
    if [ -e POTCAR -a $overwrite != True ]; then error "POTCAR existing; use -o to overwrite"; fi
    #echo low:$lowpath
    path=`find -L $lowpath -type f -name POTCAR -print -quit`
    [ "$path" != "" ] && [ -e "$path" ] && cp $path .
    echo "POTCAR found and copied"; 
    [ "`getOption -a`" != "True" ] && exit
    fi

##################################################################
# xxx checks (if all input files available, if in folder high)"
##################################################################
[ "$debug" = "yes" ] && echo 123
check parameters.dat
check INCAR
check KPOINTS
check POTCAR


# exit when an option but -c or -a
[ "`getOption -a`" = "True" ] && exit
[ "`getOption -p`" = True ] || [ "`getOption -i`" = True ] || [ "`getOption -k`" = True ] || [ "`getOption -a`" = True ] || [ "`getOption -t`" = True ] && exit  


if [ `getOption -c` = True ] || [ `getOption -e` = True ];then

###############################
#echo "# creating refpath folder ... (WAVECAR will be written! this speeds up UPTILDS)"
###############################
rm -f jobList
[ "$lom" != "Linux" ] && echo dont create jobs on a mac.... this couses trouble && exit
[ ! -e "$refpath" ] && mkdir -p $refpath 
[ ! -w "$refpath" ] && echored ERROR: refpath: $refpath : no write permissions && exit
## get existing refpathjobang
    for a in $aLats;do
        echo a: $a
        folder=$refpathjob/$a\Ang
        if [ -e "$folder" ];then
            echo $folder exists
            continue
        else
           echo creating refpath folder: $folder
           mkdir -p $folder
           cp KPOINTS $folder/KPOINTS
           cp INCAR $folder/INCAR
           sed -i $add 's|.*xxxkp1xxx.*|'"$kpoints"'|' $folder/KPOINTS
           sed -i $add 's|.*NBANDS.*|NBANDS = '"$nbands"'|' $folder/INCAR
           sed -i $add 's|.*ENCUT.*|ENCUT = '"$cutoff"'|' $folder/INCAR
           sed -i $add 's|.*EDIFF.*|EDIFF = 1e-6|' $folder/INCAR
#           sed -i $add 's|.*LWAVE.*|LWAVE=T|' $folder/INCAR
           cp POTCAR $folder/POTCAR
           echo $hessefolder/POSCAR_$a
           cp $hessefolder/POSCAR_$a $folder/POSCAR
           touch jobList
           echo $folder >> jobList
        fi
    done

[ "$debug" = "yes" ] && echo 888
###############################
echo ""
echo "# creating UPTILD folder ... (LINK TO WAVECAR will not be created)"
###############################

      for a in $aLats; do for t in $temps; do for l in $lambdas; do ## schleife ueber a=aLats, t=temps, l=lambdas, n=anzahljobs
      nthistemp="no"
      [ "$number_ns" = "1" ] && nthistemp=$n
      if [ "$number_ns" != "1" ];then
            [ "$number_ns" != "$number_temps" ] && echo "n has to be either one number (fall all temperatures) or there have to be as many temperatures as n values" && exit
            kkk=($temps)
            xxx=($n)
            len=`echo ${#kkk[*]}`
            for i in `awk 'BEGIN {x=-1; while(++x<='"$len"'-1){print x; }; exit}'`;do
                #echo "$i "Temp:"${kkk[$i]}" K,        n:"${xxx[$i]}"
                [ "$t" = "${kkk[$i]}" ] && nthistemp=${xxx[$i]}
            done
      fi
      [ "$nthistemp" = "no" ] && echo n not found for this temperature && exit
      done
      done
      done



for a in $aLats; do 
    outc=`find -L $lowpath/${a}Ang* -name "OUTCAR*" -print -quit`
    cell_=`OUTCAR_cell-last-cartesian-ARRAY.sh $outc`
    [ "`echo "$cell_" | wc -w`" != "9" ] && echo "$cell" .. cell does not have 9 words && exit

    
  for t in $temps; do for l in $lambdas; do ## schleife ueber a=aLats, t=temps, l=lambdas, n=anzahljobs
  folder=$a\Ang_$t\K
  folder_high_lambda=$folder/lambda$l
  #echo atl: $a $t $l
  #echo -en "\r$folder_high_lambda | "

  ## 1. checke in jedem $folder ob schon high_jobs da sind; wieviele fehlen noch??
  ## 1. checke in jedem $folder ob schon high_jobs da sind; wieviele fehlen noch??
  jobs_avail=""
  [ -e "$folder_high_lambda" ] && jobs_avail=`find -L $folder_high_lambda -maxdepth 1 -type d -name "*_*" | sed 's|.*/||'`
  jobs_avail_anz=`echo $jobs_avail | wc -w | sed 's|[ ]*||g'`

  ## now: change n depending in which temperature we are in!!!
  nthistemp="no"
  [ "$number_ns" = "1" ] && nthistemp=$n
  if [ "$number_ns" != "1" ];then
        [ "$number_ns" != "$number_temps" ] && echo "n has to be either one number (fall all temperatures) or there have to be as many temperatures as n values" && exit
        kkk=($temps)
        xxx=($n)
        len=`echo ${#kkk[*]}`
        for i in `awk 'BEGIN {x=-1; while(++x<='"$len"'-1){print x; }; exit}'`;do
            #echo "$i "Temp:"${kkk[$i]}" K,        n:"${xxx[$i]}"
            [ "$t" = "${kkk[$i]}" ] && nthistemp=${xxx[$i]}
        done
  fi
  [ "$nthistemp" = "no" ] && echo n not found for this temperature && exit

  jobs_missg_anz=` expr $nthistemp - $jobs_avail_anz`; [ "$jobs_missg_anz" -lt "0" ] && jobs_missg_anz=0
  [ "`getOption -un`" = "True" ] && jobs_missg_anz=`getValue -un`
  jobs_created=0
  #echo 0
  #echo -e "\r$folder_high_lambda | existing jobs: $jobs_avail_anz"
  [ "$jobs_missg_anz" -le "0" ] && echo -e "\r$folder_high_lambda | existing jobs: $jobs_avail_anz | missing: $jobs_missg_anz | created: $jobs_created | OK" && continue
  ## 2. get all possible (uncorrelated) jobs to create; look for structures_uncor_lambda$l
  ## 2. get all possible (uncorrelated) jobs to create; look for structures_uncor_lambda$l
  uncor_path=$lowpath/$folder/structures_uncor_lambda$l
  #echo path to uncorrelated structures: $uncor_path
  [ ! -e "$uncor_path" ] && ti_low_3_get_structures_uncor.sh $a $t $l > /dev/null
  [ ! -e "$uncor_path" ] && echo -e "\r$folder_high_lambda | --> NO LOW JOB <-- did you run ti_auswerten.sh in the low_ run? -> this is necessary" && continue   # uncor structures $uncor_path cant be created" && continue
  all_uncor=`cat $uncor_path`

  ## 3. get jobs to create but remove existing high jobs from the list
  ## 3. get jobs to create but remove existing high jobs from the list
  status=""
          for job in $all_uncor;do  ## im "job" steht nur 21342_31 oder so
          createfolder=""
          #echo jobsav:$jobs_avail
          [ "`echo $jobs_avail | grep -o $job | wc -w | sed 's|[ ]*||g'`" != "0" ] && status="$status $job\_already_exists" && continue
          createfolder=`pwd -L`/$folder_high_lambda/$job; [ -e "$createfolder" ] && status="$status $job\_high_exists" && continue

          seed=`echo $job | sed 's|\(.*\)_\(.*\)|\1|'`
          struct=`echo $job | sed 's|\(.*\)_\(.*\)|\2|'`
          #echo lp: $lowpath
          low_folder=$lowpath/$folder_high_lambda\_$seed; 
          #echo low_folder: $low_folder
          [ ! -e "$low_folder" ] && $lowpath/$folder_high_lambda\_$seed/workDir
          [ ! -e "$low_folder" ] && status="$status $low_folder not found" && continue

          # change Blazej 28.07.2014
          # low POSCAR was only needed for NIONS; changing this now to OUTCAR.gz, since POSCAR is removed by VASP_zipOUTPUT.sh
          #low_folder_POSCAR=$low_folder/POSCAR; [ ! -e "$low_folder_POSCAR" ] && status="$status $job\_lowpath_POSCAR_not_found" && continue
          low_folder_OUTCAR=$low_folder/OUTCAR.gz; [ ! -e "$low_folder_OUTCAR" ] && status="$status $job\_lowpath_OUTCAR.gz_not_found" && continue

          #echo ll: $low_folder
          # low_folder_structures=/nas/glensk/v/pp/al/ti_bulk_fcc4/low_2x2x2sc_250eV_02x02x02kp_EDIFF1E-1/4.13Ang_934K/lambda0.0_21057/structures_vasprun.gz
          low_folder_structures=`find -L $low_folder -maxdepth 1 -mindepth 1 -type f -name "structures*" -print -quit`; 
          echo struc: $low_folder_structures

          #exit
          [ ! -e "$low_folder_structures" ] && low_folder_structures=`find -L $low_folder/workDir -maxdepth 1 -mindepth 1 -type f -name "structures*" | head -1`; 



          [ ! -e "$low_folder_structures" ] && status="$status $job\_lowpath_structures_not_found" && continue
          [ "$low_folder_structures" = "" ] && status="$status $job\_lowpath_structures_not_found" && continue

          ## 4. create folder
          ## 4. create folder
#echo  4   
           mkdir -p $createfolder
           echo manual > $createfolder/POSCAR
           echo 1.0 >> $createfolder/POSCAR
           echo "$cell_" >> $createfolder/POSCAR
           
           # change Blazej 28.07.2014
           # see few lines above
           #head -6 $low_folder/POSCAR | tail -1 >> $createfolder/POSCAR
	 
	   #changed by sascha 2.7.2015
           nAtoms=`OUTCAR_species.sh $low_folder/OUTCAR.gz`
           echo $nAtoms >> $createfolder/POSCAR
           
           echo Cartesian >> $createfolder/POSCAR
           ## check if scaling factor is 1 !!!
           scaling=`sed -n '2p' $createfolder/POSCAR`
           #echo scaling:$scaling:
           first=`echo $scaling | cut -c 1`
           #echo first:$first:
           
           #if [ "`echo $scaling | awk '{print $1-1}'`" != "0" ];then
           #if [ "$first" != "-" ];then
           #     echo changing vecotor of POSCAR ... 
           #     sed -i $add '2 s|.*|1.0|' $createfolder/POSCAR
           #     vec1=`sed -n '3p' $createfolder/POSCAR | awk '{a='"$scaling"';printf "%.10f %.10f %.10f\n",$1*a,$2*a,$3*a}'`;sed -i $add '3 s|.*|'"$vec1"'|' $createfolder/POSCAR
           #     vec2=`sed -n '4p' $createfolder/POSCAR | awk '{a='"$scaling"';printf "%.10f %.10f %.10f\n",$1*a,$2*a,$3*a}'`;sed -i $add '4 s|.*|'"$vec2"'|' $createfolder/POSCAR
           #     vec3=`sed -n '5p' $createfolder/POSCAR | awk '{a='"$scaling"';printf "%.10f %.10f %.10f\n",$1*a,$2*a,$3*a}'`;sed -i $add '5 s|.*|'"$vec3"'|' $createfolder/POSCAR
           #     #volold=`POSCAR_volume.sh $low_folder/POSCAR`
           #     #volnew=`POSCAR_volume.sh $createfolder/POSCAR`

           #      #echo ok
           #      #echo check volume!! with low!
           #     fi
           #     fi
#echo 5

           zgrep -a --text "^$struct " $low_folder_structures | awk '{print $2,$3,$4}' >> $createfolder/POSCAR
           cp KPOINTS $createfolder/KPOINTS
           sed -i $add 's|.*xxxkp1xxx.*|'"$kpoints"'|' $createfolder/KPOINTS
           cp INCAR $createfolder/INCAR
           sed -i $add 's|.*NBANDS.*|NBANDS = '"$nbands"'|' $createfolder/INCAR
           sed -i $add 's|.*ENCUT.*|ENCUT = '"$cutoff"'|' $createfolder/INCAR
           cp POTCAR $createfolder/POTCAR
           sed -i $add '7 s|.*|Cartesian|' $createfolder/POSCAR
           #ln -s $refpathjob/$a\Ang/WAVECAR $createfolder/WAVECAR   # $refpath could be a /nas/$USER/ file (and link)

           echo $createfolder >> jobList
           
           jobs_created=` expr $jobs_created + 1 `; [ "$jobs_created" = "$jobs_missg_anz" ] && break
           done


  ## 4. check how many jobs created / still missing
  ## 4. check how many jobs created / still missing
  jobs_still_missing=` expr $jobs_missg_anz - $jobs_created `
  [ "$jobs_still_missing" = "0" ]  && echo -e "\r$folder_high_lambda | existing jobs: $jobs_avail_anz | missing: $jobs_missg_anz | created: $jobs_created | $cutoff eV | $kpoints kp | OK"
  [ "$jobs_still_missing" != "0" ] && echo -e "\r$folder_high_lambda | existing jobs: $jobs_avail_anz | missing: $jobs_missg_anz | created: $jobs_created | $cutoff eV | $kpoints kp | STILL MISSING: $jobs_still_missing"
  if [ "$jobs_still_missing" != "0" ];then
  echo "  status of all structures_uncor_lambda$l:" 
  for stat in $status; do 
  echo "  --> $stat"
  done 
  echo;
  fi
done; done; done; ## schleife ueber a=aLats, t=temps, l=lambdas, n=anzahljobs
fi
echo 
