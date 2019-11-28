#!/bin/bash

hier=`pwd`
out=no #yes #(print additional info for debugging when yes)
path=`set | grep BASH_SOURCE | sed 's|.*\"\(.*\)/[^/]*\".*|\1|' | sed '/BASH_SOURCE=/s/.*/\./'`;[ "$out" = "yes" ] && echo path: $path
script=`set | grep BASH_SOURCE | sed 's|.*\".*/\([^/]*\)\".*|\1|' | sed '/BASH_SOURCE=/s/.*/\./'`;[ "$out" = "yes" ] && echo script: $script
options=$*; . $path/../utilities/functions.include; checkOptions "-h -help -hi -hit -p -i -k -o -a -d -c -e -gn";[ "$out" = "yes" ] && echo options: $options


help=`getOption -help`
if [ $help == True ]; then
    details $script
echo "1) mkdir low_2x2x2sc_300eV_3x3x3kp_EDIFF1E-3    "
echo "   create the low folder and use following foldername: low_2x2x2sc_300eV_3x3x3kp_EDIFF1E-3"
echo "   (if you use this naming for the folder, this scripts will fill out the corresponding variables automatically)"
echo "   this name just tells you:"
echo "   2x2x2 = you use a 2x2x2 supercell"
echo "   300eV = cutoff"
echo "   3x3x3kp = kpoints"
echo "   EDIFF1E-3 = accuracy for electronic convergence"
echo ""
echo "2) ln -s path_to_quasiharmonic_referece Hessematrix_2x2x2sc   # e.g. ls -s ../PBE_dynmat_fcc/2x2x2_220eV..._0.0eV"
echo "   if you do ls now you should see 2 folders:"
echo "   Hessematrix_2x2x2sc"
echo "   low_2x2x2sc_220eV_2x2x2kp_EDIFF1E-2"
echo "   alternatively you can create the folder Hessematrix_2x2x2sc and copy necessary HesseMatrix_sphinx_xx," 
echo "   POSCAR_xx and so on there  (mkdir Hessematrix_2x2x2sc) "
echo "3) cd low_2x2x2sc_220eV_2x2x2kp_EDIFF1E-2"
echo "4) ti_low_0_create_Folders_vasp.sh (you might get a message saying that NBANDS was set to ???)"
echo "5) vi parameters.dat # change now the NBANDS to the amount of NBANDS you need."
echo "6) change values in parameters.dat to the values you want to use for your Molecular dynaics run;"
echo "   the values for alats have to also in the reference folder Hessematrix_2x2x2sc "
echo "   e.g.:"
echo "   aLats= 4.01 4.03  # put in the alats you want"
echo "   temps=1000 1500 1830 "
echo "   lambdas=0.0 0.15 0.5 0.85 1.0"
echo "   nSeeds=3"
echo "   hessepath=/home/sahara/Pd/PBE_anharmonic_fcc/Hessematrix_2x2x2sc"
echo ""
echo "   preequilibration=2000"
echo "   reftype=hesse"
echo "   gamma=0.01"
echo "   timesetp=10  # molecular dynamic timestep in fs"
echo "   ionicsteps=5000"
echo ""
echo "   cutoff=220"
echo "   kp=2 2 2"
echo "   NBANDS=180"
echo "   EDIFF=1E-2"
echo "   NGXF=NONE"
echo "   NGYF=NONE"
echo "   NGZF=NONE"
echo ""
echo "7) ti_low_0_create_Folders_vasp.sh -k"
echo "8) ti_low_0_create_Folders_vasp.sh -i"
echo "9) cp /home/sahara/Pd/PBE_anharmonic_fcc/Hessematrix_2x2x2sc/4.01Ang/POSCAR /home/sahara/Pd/PBE_anharmonic_fcc/Hessematrix_2x2x2sc/POSCAR_4.01"
echo "   cp /home/sahara/Pd/PBE_anharmonic_fcc/Hessematrix_2x2x2sc/4.03Ang/POSCAR /home/sahara/Pd/PBE_anharmonic_fcc/Hessematrix_2x2x2sc/POSCAR_4.03"
echo "   # basically copy all POSCARS which you want to calculate (and which are in parameters.dat) to the current folder and name them POSCAR_xx where xx ist the alat you specified in parameters.dat"
echo "  " 
echo "10) cp ~/Thermodynamics/utilities/fcc/EqCoords_direct_fcc_2x2x2sc /home/sahara/Pd/PBE_anharmonic_fcc/Hessematrix_2x2x2sc/EqCoords_direct_4.01"
echo "11) cp ../../PBE_dynmat_fcc_POTCAR  .  # copy your POTCAR in this folder"
echo "12) ti_low_0_create_Folders_vasp.sh    # if it says done at the end it means that your input files are OK! congrats :)"
echo "13) ti_low_0_create_Folders_vasp.sh -c"
echo "14) submit.sh <TAB> run.vasp5-tdi.par.cluster.cmmd"
exit
fi


h=`getOption -h`
if [ "$h" = "True" ];then
    usage "`basename $0` [OPTION] (if no option is specified: check of inputfiles)"
  printOptions " " \
  	       "-p       create parameters.dat template and exit" \
               "-o       overwrite file if it exists" \
	       " " \
               "-k       create KPOINTS template and exit" \
               "-i       create INCAR template and exit" \
               "-a       both above and exit (use this with -o to overwrite existing files)" \
	           " " \
	           " " \
               "-d       debug mode for plotting inbetween information" \
               "-hi      echo all corresponding high folder and exit" \
               "-hit     echo corresponding high TAKE folder or only if one highfolder exists and exit" \
               "-e       echo hessepath and exit" \
               "-gn      echo nbands and exit" \
               " " \
               "-c       create folders" 
  exit
fi

[ `getOption -d` = True ] && debug=yes
#####################################################################
# -hi  -hit (list all high folder go to high TAKE folder)
#####################################################################
if [ "`getOption -hi`" = "True" ] || [ "`getOption -hit`" = "True" ];then
    hier=`pwd`
    sub=`pwd | sed 's|.*/||'`
    #echo sub:$sub:
    cd ../
    hier1=`pwd`
        found=`find $sub* -maxdepth 0 -mindepth 0 -type d -name "*__high_*"`
        allout=""
        highTAKEout=""
        for i in $found;do
            cd $hier1
            cd $i
            [ "`pwd`" != "$hier" ] && allout="$allout `pwd`" && [ "`getOption -hi`" = "True" ] && echo `pwd`
        done
    cd $hier

    if [ "`getOption -hit`" = "True" ];then
        allout=`echo $allout | xargs -n1`
        [ "`echo $allout | wc -w | sed 's|[ ]*||'`" = "1" ] && echo $allout && exit
        if [ "`echo $allout | wc -w | sed 's|[ ]*||'`" != "1" ];then
            hit=`echo "$allout" | grep 'TAKE$'`
            [ "`echo $hit | wc -w | sed 's|[ ]*||'`" = "1" ] && echo $hit && exit
        fi
    fi
    exit
fi

if [ "`getOption -e`" = "True" ];then
    hessepath=`grep hessepath= parameters.dat | sed 's|hessepath=||'`
    [ -e "$hessepath" ] && echo $hessepath && exit

    sc_string_tmp=`pwd | grep -o "[0-9]*x[0-9]*x[0-9]*sc"`
    hessefolder="xxx"
    hessefoldercheck=../Hessematrix_$sc_string_tmp
    [ "`getOption -d`" = "True" ] && echo hesefoldercheck: $hessefoldercheck
    [ -e "$hessefoldercheck" ] && hier=`pwd` && cd $hessefoldercheck && hessefolder=`pwd` && cd $hier
    [ "`getOption -d`" = "True" ] && echo hesefolder: $hessefolder
    [ ! -e "$hessefolder" ] && echored "hessefolder does not exist !!"
    echo $hessefolder
    exit
fi
########################################################################
# xxx create KPOINTS template and exit
########################################################################
genPar=`getOption -k`
overwrite=`getOption -o`
aop=`getOption -a`
if [ "$genPar" = "True" ] || [ "$aop" = "True" ]; then
if [ -e KPOINTS -a $overwrite != True ]; then error "KPOINTS existing; use -o to overwrite"; fi
echo "K-Points" > KPOINTS
echo " 0" >> KPOINTS
echo "Monkhorst Pack" >> KPOINTS
echo "xxkpxx xxkpxx xxkpxx" >> KPOINTS
echo "0 0 0" >> KPOINTS
echo; echo "KPOINTS written";
[ "`getOption -k`" = "True" ] && exit
fi


########################################################################
# xxx create INCAR template and exit
########################################################################
genPar=`getOption -i`
overwrite=`getOption -o`
aop=`getOption -a`
if [ "$genPar" = "True" ] || [ "$aop" = "True" ]; then
if [ -e INCAR -a $overwrite != True ]; then error "INCAR existing; use -o to overwrite"; fi
echo "LAMBDA       = xxxLAMBDAxxx" > INCAR
echo "REF_TYPE     = xxxREF_TYPExxx" >> INCAR
echo "PRE_EQ_N     = xxxPREEQUILIBRATIONxxx" >> INCAR
echo "GAMMA_LD     = xxxGAMMA_LDxxx" >> INCAR
echo "SEED         = xxxSEEDxxx" >> INCAR
echo "POTIM  = xxxPOTIMxxx" >> INCAR
echo "TEBEG  = xxxTEMPxxx" >> INCAR
echo "NSW    = xxxNSWxxx" >> INCAR
echo " " >> INCAR
echo " " >> INCAR
echo "NPAR=1" >> INCAR
echo "" >> INCAR
echo "ENCUT  =    xxxCUTOFFxxx" >> INCAR
echo "ISMEAR =     -1" >> INCAR
echo "SIGMA  =    0.1" >> INCAR
echo "EDIFF  =    xxxEDIFFxxx" >> INCAR
echo "" >> INCAR
echo "ISYM   =      0" >> INCAR
echo "ADDGRID=    TRUE" >> INCAR
echo "PREC   =    Accurate" >> INCAR
echo "LREAL =     .FALSE." >> INCAR
echo "NBANDS =    xxxNBANDSxxx" >> INCAR
echo "ALGO   =    FAST" >> INCAR
echo "" >> INCAR
echo "LWAVE  =      F" >> INCAR
echo "LCHARG =      F" >> INCAR
echo "LVTOT = F" >> INCAR
echo "LELF = F" >> INCAR
echo "NWRITE = 0" >> INCAR

#[ "$NGXF" != "none" ] && echo "NGXF=$NGXF" >> INCAR
#[ "$NGYF" != "none" ] && echo "NGYF=$NGYF" >> INCAR
#[ "$NGZF" != "none" ] && echo "NGZF=$NGZF" >> INCAR

echo; echo "INCAR written";
[ "`getOption -i`" = "True" ] && exit
fi


##########################################################
# xxx create parameters.dat template and exit
##########################################################
if [ `getOption -p` = True ] || [ "`getOption -gn`" = "True" ];then
    if [ "$debug" = "yes" ];then
        echo "#######################"
        echo wring parameters.dat
        echo "#######################"
    fi  
    overwrite=`getOption -o`
    sc_string_tmp=`pwd | grep -o "[0-9]*x[0-9]*x[0-9]*sc"`
    hessefolder="xxx"
    hessefoldercheck=../Hessematrix_$sc_string_tmp
    [ "`getOption -d`" = "True" ] && echo hesefoldercheck: $hessefoldercheck
    [ -e "$hessefoldercheck" ] && hier=`pwd` && cd $hessefoldercheck && hessefolder=`pwd` && cd $hier
    [ "`getOption -d`" = "True" ] && echo hesefolder: $hessefolder
    [ ! -e "$hessefolder" ] && echored "hessefolder does not exist !!"
    parameterssc=../parameters_$sc_string_tmp
    gamma="xxx";timestep="xxx";ionicsteps="xxx";reftype="xxx";
    alats="xxx";[ -e "$hessefolder" ] && alats_tmp=`find $hessefolder -name 'HesseMatrix_*' | sed 's|.*HesseMatrix_||' | xargs` && [ "`isnumber.sh $alats_tmp`" = "yes" ] && alats=$alats_tmp; [ "$alats" = "xxx" ] && alats_tmp=`find . -maxdepth 1 -mindepth 1 -type d -name "*Ang_*K" | sed 's|./||g' | sed 's|Ang_.*||' | sort | uniq | xargs`
    temp="xxx"
    temps_tmp=`find . -maxdepth 1 -mindepth 1 -type d -name "*Ang_*K" | sed 's|./||g' | sed 's|.*Ang_||' | sed 's|K||g' | sort | uniq | xargs`
    [ "`isnumber.sh $temps_tmp`" = "yes" ] && temp=$temps_tmp

    [ "`getOption -d`" = "True" ] && echo alats_tmp:$alats_tmp
    [ "`isnumber.sh $alats_tmp`" = "yes" ] && [ "$alats" = "xxx" ] && alats=$alats_tmp
    [ -e "$hessefolder" ] && outcar_tmp=`find -L $hessefolder -name 'OUTCAR*' -print -quit` 
    [ "`getOption -d`" = "True" ] && echo outcar_tmp:$outcar_tmp
    [ -e "$outcar_tmp" ] && element_tmp=`OUTCAR_elements.sh $outcar_tmp` 
    [ "`getOption -d`" = "True" ] && echo ele: $element_tmp
    [ "`echo $element_tmp | grep -o error`" = "error" ] && element_tmp=`get.py -fpse`
    [ "`getOption -d`" = "True" ] && echo elepy: $element_tmp :$elements_tmp
	incarjob_tmp=`find -L *Ang_*K/lambda*_*/ -name INCAR -print -quit`
	[ -e "$incarjob_tmp" ] && timestep_tmp=`grep POTIM $incarjob_tmp | awk '{print $3}'`
	[ -e "$incarjob_tmp" ] && gamma_tmp=`grep GAMMA_LD $incarjob_tmp | awk '{print $3}'`
	[ -e "$incarjob_tmp" ] && reftype_tmp=`grep REF_TYPE $incarjob_tmp | awk '{print $3}'`
	[ -e "$incarjob_tmp" ] && nsw_tmp=`grep NSW $incarjob_tmp | awk '{print $3}'`
	[ -e "$incarjob_tmp" ] && bands_tmp=`grep NBANDS $incarjob_tmp | awk '{print $3+5}'`
    [ "`isnumber $timestep_tmp`" = "yes" ] && timestep=$timestep_tmp
    [ "`isnumber $gamma_tmp`" = "yes" ] && gamma=$gamma_tmp
    [ "`isnumber $reftype_tmp`" = "yes" ] && reftype=$reftype_tmp
    [ "`isnumber $nsw_tmp`" = "yes" ] && ionicsteps=$nsw_tmp
    [ "`isnumber $bands_tmp`" = "yes" ] && nbands=$bands_tmp

    [ "$temp" = "xxx" ] && [ "`echo $element_tmp | wc -w | sed 's|[ ]*||'`" -eq "1" ] && temp_check=`getMeltingPoint.sh $element_tmp -r` 

    [ "`getOption -d`" = "True" ] && echo kkkkkkkkkkk
    [ "$temp" = "xxx" ] && [ "`isnumber.sh $temp_check`" = "yes" ] && temp=$temp_check

    [ "`getOption -d`" = "True" ] && echo ----------------------
    [ "`getOption -d`" = "True" ] && echo temp_check: $temp_check
    [ -z "$nbands" ] && nbands="???"; [ -e "$hessefolder/INCAR" ] && nbands_tmp=`INCAR_NBANDS.sh $hessefolder/INCAR` && [ "`isnumber.sh $nbands_tmp`" = "yes" ] && nbands=$nbands_tmp

    [ "`getOption -d`" = "True" ] && echo nbands_tmp:$nbands_tmp nbands:$nbands


    ## echo stuff
    [ "`getOption -gn`" = "True" ] && echo $nbands && exit

    ## write stuff
    if [ -e parameters.dat -a $overwrite != True ]; then error "parameters.dat existing; use -o to overwrite"; fi
    echo "# GENERAL part" > parameters.dat
    echo "aLats=$alats" >> parameters.dat
    echo "temps=$temp" >> parameters.dat
    echo "lambdas=0.0 0.15 0.5 0.85 1.0" >> parameters.dat
    echo "nSeeds=3" >> parameters.dat
    echo "hessepath=$hessefolder" >> parameters.dat
    echo "" >> parameters.dat
    echo "" >> parameters.dat
    echo "# INCAR and KPOINTS part" >> parameters.dat
    [ -e "$parameterssc" ] && [ -z "$gamma" ] && gamma=`grep gamma= $parameterssc | sed 's|gamma=||' | awk '{print $1}'` 
    [ -e "$parameterssc" ] && [ -z "$reftype" ] && reftype=`grep reftype= $parameterssc | sed 's|reftype=||' | awk '{print $1}'` 
    [ -e "$parameterssc" ] && [ -z "$timestep" ] && timestep=`grep timestep= $parameterssc | sed 's|timestep=||' | awk '{print $1}'` 
    [ -e "$parameterssc" ] && [ -z "$nbands" ] && ionicsteps=`grep ionicsteps= $parameterssc | sed 's|ionicsteps=||' | awk '{print $1}'` 
    [ -e "$parameterssc" ] && [ "$nbands" = "???" ] && nbands=`grep NBANDS= $parameterssc | sed 's|NBANDS=||' | awk '{print $1}'` 

#	echo timestep:$timestep:
    echo "preequilibration=1000        # time steps for preequilibration, this is fast and can easily be 10000" >> parameters.dat
    echo "reftype=$reftype               # used reference: hesse or extern" >> parameters.dat
    echo "gamma=$gamma                  # friction parameter of Langevin Thermostat (e.g. 0.01)" >> parameters.dat
    echo "timestep=$timestep                # timestep of ionic motion in fs (e.g. 5 10)" >> parameters.dat
    echo "ionicsteps=$ionicsteps            # how many ionic speps should be performed in md (e.g. 4000)" >> parameters.dat
    echo "" >> parameters.dat
    cutoff_pfad=`echo $hier | sed 's|.*[-_]\([0-9]*\)eV.*|\1|'`; [ "`isnumber $cutoff_pfad`" != "yes" ] && cutoff_pfad="???"
    kpoints_pfad=`echo $hier | grep -o "[0-9]*[xmg][0-9]*[xmg][0-9]*kp\|kp[0-9]*[xmg][0-9]*[xmg][0-9]*" | tail -1 | sed 's|[xmkp]| |g' | sed 's|^[0]*||' | sed 's| 0| |g' | sed -e 's/^[ \t]*//' | sed 's/^[ \t]*//;s/[ \t]*$//'`;[ "`isnumber $kpoints_pfad`" != "yes" ] && kpoints_pfad="???"
    echo "cutoff=$cutoff_pfad                    # cutoff in eV: e.g. 500" >> parameters.dat
    echo "kp=$kpoints_pfad                      # kpoints: e.g. 4 4 4" >> parameters.dat
    echo "NBANDS=$nbands                    # INCAR NBANDS" >> parameters.dat
    echored "TAKE CARE, NBANDS in parameters.dat was set to $nbands which is taken from $hessefolder/INCAR; please check this carefully "
    ediff="xxx"; ediff_tmp=`pwd | grep -o "EDIFF.E.." | sed 's|EDIFF||'`; [ "$ediff_tmp" != "" ] && ediff=$ediff_tmp
    echo "EDIFF=$ediff                    # INCAR EDIFF e.g. 1E-3" >> parameters.dat
    ngxf="none";ngxf_tmp=`pwd | grep -o "[0-9]*NGXF" | sed 's|NGXF||'`;[ "`isnumber.sh $ngxf_tmp`" = "yes" ] && ngxf=$ngxf_tmp
    echo "NGXF=$ngxf                     # INCAR NGXF; leafe none if no specification for NGXF" >> parameters.dat
    echo "NGYF=$ngxf                     # INCAR NGYF; leafe none if no specification for NGXF" >> parameters.dat
    echo "NGZF=$ngxf                     # INCAR NGZF; leafe none if no specification for NGXF" >> parameters.dat
    exit
    fi

##################################################################
# xxx read parameters.dat
##################################################################
[ ! -e "INCAR" ] && echo "please create INCAR template using `basename $0` -i" && exit
[ ! -e "KPOINTS" ] && echo "please create KPOINTS template using `basename $0` -k" && exit
echo checking INCAR and KPOINTS template
kpointsstring=`grep -o "[x]*KP*[x]* [x]*KP[x]* [x]*KP[x]*" KPOINTS`
[ "$kpointsstring" = "" ] && kpointsstring=`grep -o "[x]*kp*[x]* [x]*kp[x]* [x]*kp[x]*" KPOINTS`

[ "`echo $kpointsstring | wc -w | sed 's|[ ]*||g'`" != "3" ] && echo "not found in KPOINTS: xxkpxx xxkpxx xxkpxx but found: $kpointsstring" && exit

## checking INCAR with parameters.dat values
[ ! -e "INCAR" ] && echo INCAR missing && exit
[ "`grep -o "^NBANDS =    xxxNBANDSxxx" INCAR | wc -w | sed 's|[ ]*||g'`" != "3" ] && echo "not found in INCAR: NBANDS =    xxxNBANDSxxx" && exit
[ "`grep -o "^NSW    = xxxNSWxxx" INCAR | wc -w | sed 's|[ ]*||g'`" != "3" ] && echo "not found in INCAR: NSW = xxxNSWxxx" && exit
[ "`grep -o "^GAMMA_LD     = xxxGAMMA_LDxxx" INCAR | wc -w | sed 's|[ ]*||g'`" != "3" ] && echo "not found in INCAR: GGAMMA_LD     = xxxGAMMA_LDxxx" && exit
[ "`grep -o "^POTIM  = xxxPOTIMxxx" INCAR | wc -w | sed 's|[ ]*||g'`" != "3" ] && echo "not found in INCAR: POTIM  = xxxPOTIMxxx" && exit
[ "`grep -o "^EDIFF  =    xxxEDIFFxxx" INCAR | wc -w | sed 's|[ ]*||g'`" != "3" ] && echo "not found in INCAR: EDIFF  =    xxxEDIFFxxx" && exit
[ "`grep -o "^REF_TYPE     = xxxREF_TYPExxx" INCAR | wc -w | sed 's|[ ]*||g'`" != "3" ] && echo not found in INCAR: "REF_TYPE     = xxxREF_TYPExxx" && exit
[ "`grep -o "^PRE_EQ_N     = xxxPREEQUILIBRATIONxxx" INCAR | wc -w | sed 's|[ ]*||g'`" != "3" ] && echo "please add to INCAR:\"PRE_EQ_N = xxxPREEQUILIBRATIONxxx\"  I suggest to use `basename $0` -i" && exit
echo checking INCAR and KPOINTS template .... done

echo checking parameters.dat
if [ "$debug" = "yes" ];then
    echo "#######################"
    echo reading from parameters.dat
    echo "#######################"
fi
check parameters.dat
[ "$debug" = "yes" ] && echo parameter.at checking ...
aLats=`get aLats`; temps=`get temps`; lambdas=`get lambdas`;cutoff=`get cutoff`;nbands=`get NBANDS`;reftype=`get reftype`;gamma=`get gamma`;timestep=`get timestep`;ionicsteps=`get ionicsteps`;ediff=`get EDIFF`;preequilibration=`get preequilibration`
nSeeds=`get nSeeds`;kpoints=`get kp | sed 's/^[ \t]*//;s/[ \t]*$//'`;hessepath=`get hessepath`; NGXF=`get NGXF`; NGYF=`get NGYF`; NGZF=`get NGZF`
for i in "$aLats" "$temps" "$lambdas" "$cutoff" "$nbands" "$nSeeds" "$kpoints" "$gamma" "$reftype" "$timestep" "$ionicsteps" "$ediff";do

    if [ -z "$i" ];then
        echo aLats:"$aLats" 
        echo temps:"$temps" 
        echo lambdas:"$lambdas"
        echo cutoff:"$cutoff" 
        echo nbands:"$nbands" 
        echo nSeeds:"$nSeeds" 
        echo kpoints:"$kpoints" 
        echo gamma:"$gamma" 
        echo reftype:"$reftype"
        echo timestep:"$timestep" 
        echo ionicstep:"$ionicsteps" 
        echo ediff:"$ediff"
        echo preequilibration:"$preequilibration"
        echored parameters.dat needs to contain: aLats, temps, cutoff, reftype, nbands, nSeeds but something is missing && exit
        fi
    done
checkInput "$aLats" "$temps" "$lambdas" "$cutoff" "$nbands" "$nSeeds" "$kpoints" "$reftype" "$gamma" "$timestep" "$ionicsteps" "$ediff" "$preequilibration"

if [ "$debug" = "yes" ];then
echo ""
echo READ FROM parameters.dat:
echo aLats:$aLats
echo temps:$temps
echo lambdas:$lambdas
echo cutoff:$cutoff
echo
echo reftype:$reftype
echo gamma:$gamma
echo timestep:$timestep
echo ionicsteps:$ionicsteps
echo
echo kpoints:$kpoints
echo nSeeds:$nSeeds
echo hessepath:$hessepath
echo nbands:$nbands
echo ediff:$ediff
echo NGXF:$NGXF
echo NGYF:$NGYF
echo NGZF:$NGZF
echo ""
fi
[ "`isnumber $aLats`" != "yes" ] && echored "ERROR: aLats:$aLats: in parameters.dat is no number. Create parameters.dat using `basename $0` -p -o" && exit
[ "`isnumber $temps`" != "yes" ] && echored "ERROR: temps:$temps: in parameters.dat is no number. Create parameters.dat using `basename $0` -p -o" && exit
[ "`isnumber $preequilibration`" != "yes" ] && echored "ERROR: preequilibration:$preequilibration: in parameters.dat is no number. Create parameters.dat using `basename $0` -p -o" && exit
[ "`isnumber $lambdas`" != "yes" ] && echored "ERROR: lambdas:$lambdas: in parameters.dat is no number. Create parameters.dat using `basename $0` -p -o" && exit
[ "`isnumber $gamma`" != "yes" ] && echored "ERROR: gamma:$gamma: in parameters.dat is no number. Create parameters.dat using `basename $0` -p -o" && exit
[ "$reftype" != "hesse" ] && [ "$reftype" != "extern"  ] && [ "$reftype" != "pot" ] && echored "ERROR: reftype:$reftype: in parameters.dat is neither hesse nor extern nor pot. Create parameters.dat using `basename $0` -p -o and adjust it" && exit
[ "`isnumber $timestep`" != "yes" ] && echored "ERROR: timestep:$timestep: in parameters.dat is no number. Create parameters.dat using `basename $0` -p -o" && exit
[ "`isnumber $ionicsteps`" != "yes" ] && echored "ERROR: ionicsteps:$ionicsteps: in parameters.dat is no number. Create parameters.dat using `basename $0` -p -o" && exit
[ "`isnumber $cutoff`" != "yes" ] && echored "ERROR: cutoff:$cutoff: in parameters.dat is no number. Create parameters.dat using `basename $0` -p -o" && exit
[ "`isnumber $kpoints`" != "yes" ] && echored "ERROR: kpoints:$kpoints: in parameters.dat is no number. Create parameters.dat using `basename $0` -p -o" && exit
[ "`echo $kpoints | wc -w | sed 's|[ ]*||g'`" != "3" ] && echored "ERROR: kpoints:$kpoints: in parameters.dat needs to have 3 INTEGERS" && exit
[ "`isnumber $nSeeds`" != "yes" ] && echored "ERROR: nSeeds:$nSeeds: in parameters.dat is no number. Create parameters.dat using `basename $0` -p -o" && exit
[ "`isnumber $nbands`" != "yes" ] && echored "ERROR: nbands:$nbands: in parameters.dat is no number. Create parameters.dat using `basename $0` -p -o" && exit
[ "$NGXF" != "none" ] && [ "`isnumber $NGXF`" != "yes" ] && echored "ERROR: NGXF:$NGXF: in parameters.dat is no number. Create parameters.dat using `basename $0` -p -o" && exit
[ "$NGYF" != "none" ] && [ "`isnumber $NGYF`" != "yes" ] && echored "ERROR: NGYF:$NGYF: in parameters.dat is no number. Create parameters.dat using `basename $0` -p -o" && exit
[ "$NGZF" != "none" ] && [ "`isnumber $NGZF`" != "yes" ] && echored "ERROR: NGZF:$NGZF: in parameters.dat is no number. Create parameters.dat using `basename $0` -p -o" && exit

[ ! -e "$hessepath" ] && echored "ERROR: hessepath: $hessepath : in parameters.dat does not exist" && exit
ediffcheck=`echo $ediff | awk '{print $1*1}'`
[ "`isnumber.sh $ediffcheck`" != "yes" ] && echored "ERROR: EDIFF:$ediff: in parameters.dat is not a number" && exit


for a in $aLats;do
    if [ ! -e "$hessepath/POSCAR_$a" ];then
    ptmp=`find -L $a\Ang_*K -mindepth 2 -maxdepth 2 -name "POSCAR" -print -quit`   
    [ -e "$ptmp" ] && cp $ptmp $hessepath/POSCAR_$a
    fi
    [ ! -e "$hessepath/POSCAR_$a" ] && echored "$hessepath/POSCAR_$a does not exist" && exit

    # check whether POSCAR contains Selective Dynamics, if so quit, because it can destroy whole MD
    # Selective Dynamics switch can be in the 7th or 8th line and only first letter matters (s or S)
    selective=`sed -n '7,8{/^[sS]/p}' $hessepath/POSCAR_$a`
    if [ "$selective" != "" ]; then
      error "$hessepath/POSCAR_$a\n       contains 'Selective Dynamics' which can destroy MD, remove (possibly in other POSCARs) and run again"
    fi

    if [ ! -e "$hessepath/HesseMatrix_$a" ];then
        ptmp=`find -L $a\Ang_*K -mindepth 2 -maxdepth 2 -name "HesseMatrix_sphinx" -print -quit`   
        [ -e "$ptmp" ] && cp $ptmp $hessepath/HesseMatrix_$a
    fi
    [ $reftype = "hesse"  ] && [ ! -e "$hessepath/HesseMatrix_$a" ] && echored "$hessepath/HesseMatrix_$a does not exist" && exit
    if [ ! -e "$hessepath/EqCoords_direct_$a" ];then
        ptmp=`find -L $a\Ang_*K -mindepth 2 -maxdepth 2 -name "EqCoord_direct" -print -quit`   
        [ -e "$ptmp" ] && cp $ptmp $hessepath/EqCoords_direct_$a
    fi
    [ ! -e "$hessepath/EqCoords_direct_$a" ] && echo $hessepath/EqCoords_direct_$a does not exist && exit
done
echo checking parameters.dat ... done 
[ "`getOption -a`" = "True" ] && exit



########################################################################
# xxx create folder
########################################################################
[ ! -e "KPOINTS" ] && echo KPOINTS missing && exit
[ ! -e "INCAR" ] && echo INCAR missing && exit
if [ "$NGXF" != "none" ] && [ "`isnumber $NGXF`" = "yes" ];then
# checke ob in INcar eine Zeile mit ngxf drin steht
    if [ "`grep -o NGXF INCAR`" = "" ];then
        echo NGXF= >> INCAR
        sed -i 's|.*NGXF.*|NGXF='"$NGXF"'|' INCAR
    fi
fi
if [ "$NGYF" != "none" ] && [ "`isnumber $NGYF`" = "yes" ];then
# checke ob in INcar eine Zeile mit ngyf drin steht
    if [ "`grep -o NGYF INCAR`" = "" ];then
        echo NGYF= >> INCAR
        sed -i 's|.*NGYF.*|NGYF='"$NGYF"'|' INCAR
    fi
fi
if [ "$NGZF" != "none" ] && [ "`isnumber $NGZF`" = "yes" ];then
# checke ob in INcar eine Zeile mit ngzf drin steht
    if [ "`grep -o NGZF INCAR`" = "" ];then
        echo NGZF= >> INCAR
        sed -i 's|.*NGZF.*|NGZF='"$NGZF"'|' INCAR
    fi
fi

#kpcheck=`KPOINTS_cat-kp.sh | sed 's|x| |g'`
#[ "$kpcheck" != "$kpoints" ] && echo kpoints in parameters.dat:$kpoints in KPOINTS:$kpcheck && exit
#echo "checking KPOINTS        ... done"

[ ! -e "POTCAR" ] && echo POTCAR missing && exit
if [ `getOption -c` = "True" ];then
rm -f jobList
#seed=`tail -n1 usedSeeds`
for a in $aLats; do for t in $temps; do for l in $lambdas; do for seedi in `seq 1 $nSeeds`;do
    EqCoords_direct=$hessepath/EqCoords_direct_$a
    POSCAR_pfad=$hessepath/POSCAR_$a
    [ ! -e "$POSCAR_pfad" ] && error $POSCAR_pfad does not exist 44 && exit
    [ ! -e "$EqCoords_direct" ] && error $EqCoords_direct does not exist 55 && exit

    seed=`date +%s%N | cut -c 14-20`
    [ "$seed" = "" ] && seed=`perl -MTime::HiRes=gettimeofday -e 'print int(1000*gettimeofday()).qq(\n);' | cut -c 8-13` # necessary on mac
    [ "$seed" = "" ] && echo no seed! && exit ## necessary to check at mac
    #echo "seed: $seed"
    [ "`isnumber.sh $cutoff $seed $l $t $a`" != "yes" ] && echo no number in cutoff:$cutoff seed:$seed l:$l t:$t a:$a && exit

    f=$a\Ang_$t\K/lambda$l\_$seed
    #rm -fr $f;

    [ -e $f ] && echo $f already exists && continue
    mkdir -p $f
    #echo ok1
    #echo cutoff:$cutoff
    #echo seed:$seed
    #echo lam:$l
    echo Temp:$t ALAT:$a LAM:$l SEED:$seed CUTOFF:$cutoff kp:$kpoints
    sed -e 's/xxxCUTOFFxxx/'$cutoff'/' \
        -e 's/xxxSEEDxxx/'$seed'/' \
        -e 's/xxxPREEQUILIBRATIONxxx/'$preequilibration'/' \
        -e 's/xxxLAMBDAxxx/'$l'/' \
        -e 's/xxxNBANDSxxx/'$nbands'/' \
        -e 's/xxxNSWxxx/'$ionicsteps'/' \
        -e 's/xxxGAMMA_LDxxx/'$gamma'/' \
        -e 's/xxxREF_TYPExxx/'$reftype'/' \
        -e 's/xxxPOTIMxxx/'$timestep'/' \
        -e 's/xxxEDIFFxxx/'$ediff'/' \
        -e 's/xxxTEMPxxx/'$t'/' INCAR > $f/INCAR
    sed -e 's/[x]*[KkpP]*xx.*/'"$kpoints"'/' KPOINTS > $f/KPOINTS



    #kpchk1=`sed -n '4p' $f/KPOINTS | awk '{print $1}'`;[ "$kpchk1" != "$kp1" ] && error "could not substitue KPOINS in $f/KPOINTS correctly" && exit

    #kpchk2=`sed -n '4p' $f/KPOINTS | awk '{print $2}'`;[ "$kpchk2" != "$kp2" ] && error "could not substitue KPOINS in $f/KPOINTS correctly" && exit
    #kpchk3=`sed -n '4p' $f/KPOINTS | awk '{print $3}'`;[ "$kpchk3" != "$kp3" ] && error "could not substitue KPOINS in $f/KPOINTS correctly" && exit
    #echo "kp1:$kp1: kp2:$kp2: kp3:$kp3:"
    #kpcheck=`grep -o "'"$kp1"'*'"$kp2"'*'"$kp3"'" $f/KPOINTS | xargs`
    #[ "`echo "$kpcheck" | wc -w | sed 's|[ ]*||g'`" != "3" ] && error "could not substitue KPOINS in $f/KPOINTS correctly kpcheck:$kpcheck: " && exit
    #cp KPOINTS $f/KPOINTS
    #sed -e 's/xxxALATxxx/'$a'/' POSCAR > $f/POSCAR
    cp $POSCAR_pfad $f/POSCAR
    cp POTCAR $f/POTCAR
    #echo ee:$EqCoords_direct 
    #echo ff:$f/EqCoords_direct

    cp $EqCoords_direct $f/EqCoords_direct
    if [ $reftype = "hesse" ]
      then cp $hessepath/HesseMatrix_$a $f/HesseMatrix_sphinx
    fi
#    cp $hessepath/HesseMatrix_$a $f/HesseMatrix_sphinx

    echo `pwd`/$f >> jobList
  #done
done; done; done; done
fi
echo done
