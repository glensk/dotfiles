#!/bin/bash
out=no
path=`set | grep BASH_SOURCE | sed 's|.*\"\(.*\)/[^/]*\".*|\1|' | sed '/BASH_SOURCE=/s/.*/\./'`;[ "$out" = "yes" ] && echo path: $path
script=`set | grep BASH_SOURCE | sed 's|.*\".*/\([^/]*\)\".*|\1|' | sed '/BASH_SOURCE=/s/.*/\./'`;[ "$out" = "yes" ] && echo script: $script
options=$*; . $path/utilities/functions.include; checkOptions "-h -help -i -k -o -p -a -alat -c -g -d -vac -fit -fitref -t -t1 -t2";[ "$out" = "yes" ] && echo options: $options

if [ `getOption -h` = True ] || [ `getOption -help` = True ] ; then
    usage "`basename $0` [OPTION]  (without an option: check inputfiles)"
    echo "" 
    echo ""
    echo "         This script needs a POTCAR and your relaxedCoords_xxx files as a starting point"

  printOptions " " \
               "-p           create parameters.dat template and exit" \
               "-i           create INCAR template and exit" \
               "-k           create KPOINTS template and exit" \
               "-o           force overwrite of file if it exists" \
               "-a           all above and exit (use this with -o to overwrite existing files)" \
               "-alat VALUE  do this only for a certaion lattice constant (e.g. -alat 3.72)" \
               "-c           create folder (this also creates the corresponding reference structure)" \
               " " \
               "-g           1. create allForces.dat files when run is finished" \
			   "-vac         2. create allForces.dat.x folder (getSingleSpediesPhonons.sh -vac for HesseMatrix and Fvib of vacancy)" \
               "-fitref      3. create fitFqh.input and fit surface for reference bulk system" \
               "-fit         4. create fitFqh.input and fit surface" \
               " " \
               "-t1          4. get EVinet(T=0K) formation energy" \
               "-t2          5. get Fqh(T) formation energy" \
               " " \
               "-d           debug mode for plotting inbetween information" 
  exit
fi

toAng=0.529177208278835
########################################################################
# xxx create parameters.dat template and exit
########################################################################
genPar=`getOption -p`
overwrite=`getOption -o`
hier=`pwd`
if [ "$genPar" = "True" ] || [ "`getOption -a`" = "True" ]; then
    if [ -e parameters.dat -a $overwrite != True ]; then error "parameters.dat existing; use -o to overwrite"; fi
 # get element
        if [ ! -e "POTCAR" ];then
            potcar=`find . -maxdepth 2 -mindepth 1 -name "POTCAR" -print -quit` 
            #echo POTCAR:$potcar
            [ -e "$potcar" ] && cp $potcar .
        fi
        
        numatomspos=""
        if [ ! -e "POSCAR" ];then
            poscar=`find . -maxdepth 2 -mindepth 1 -name "POSCAR*" -print -quit` 
            #echo POSCAR:$potcar
            [ -e "$poscar" ] && numatomspos=`POSCAR_numatoms.sh $poscar`
        fi
        
        if [ "$numatomspos" = "" ];then
        relcoords=`find . -name "relaxedCoords_*" | head -1`
        if [ "$relcoords" != "" ];then
            numatomspos=`cat $relcoords | awk '{print $1}' | wc -w`
        else
            error "please copy your relaxedCoords_xxx to this folder (e.g. relaxedCoords_3.7 ... relaxedCoords_3.8)" && exit
            
        fi
        fi

        outcar_all=`find . -maxdepth 2 -mindepth 1 -name "OUTCAR*"`
        [ "`getOption -d`" = "True" ] && echo outcar_all: $outcar_all
        outcar_tmp=`echo "$outcar_all" | head -1`
        [ "`getOption -d`" = "True" ] && echo outcar_tmp: $outcar_tmp
        outcar_alats=`echo "$outcar_all" | sed 's|.*OUTCAR||' | sed 's|^.||' | sed 's|\([0-9.]*\).*|\1|' | xargs`
        [ "`getOption -d`" = "True" ] && echo outcar_alats: $outcar_alats
        [ "$outcar_alats" = "" ] && outcar_alats=`echo "$outcar_all" | sed 's|Ang.*||' | xargs -n1 | sed 's|./||' | xargs`
        [ "`getOption -d`" = "True" ] && echo outcar_alats: "$outcar_alats"
        alats="???";[ "`isnumber.sh $outcar_alats`" = "yes" ] && alats=`echo "$outcar_alats" | sed 's| |,|g'`
        [ "$alats" = "???" ] && alats_tmp=`ls -1d [0-9.]*Ang | sed 's|Ang.*||' | xargs` && [ "$alats_tmp" != "" ] && alats=$alats_tmp

        [ "$alats" = "???" ] && relalats=`find . -maxdepth 1 -mindepth 1 -name "relaxedCoords_*" | sed 's|relaxedCoords_||g' | sed 's|./||g' | xargs` && [ "$relalats" != "" ] && alats=$relalats
        [ "`getOption -d`" = "True" ] && echo alats_tmp:$alats_tmp

        alats=`echo "$alats" | xargs -n1 | sort -n | sed 's|0*$||g'`
        [ -e "$outcar_tmp" ] && element_tmp=`OUTCAR_elements.sh $outcar_tmp`;element=XX; atoms_tmp=`OUTCAR_number_of_atoms.sh $outcar_tmp`;atoms=XX
        [ ! -e "$outcar_tmp" ] && [ -e "POTCAR" ] && element=`POTCAR_element.sh POTCAR`
        tmelt=1000;tmelt_tmp=`getMeltingPoint.sh $element_tmp -r`;[ "`isnumber.sh $tmelt_tmp`" = "yes" ] && tmelt=$tmelt_tmp && element=$element_tmp
        cutoff_pfad=`echo $hier | sed 's|.*[-_]\([0-9]*\)eV.*|\1|'`;
        
        [ "`getOption -d`" = "True" ] && echo cutoff_pfad:$cutoff_pfad
        [ "`isnumber $cutoff_pfad`" != "yes" ] && cutoff_pfad="???"
        [ "$cutoff_pfad" = "???" ] && [ -e "$outcar_tmp" ] && cutoff_pfad=`OUTCAR_ENCUT.sh $outcar_tmp`
        hier=`pwd`
        kpoints_pfad=`echo $hier | grep -o "[0-9]*[xmg][0-9]*[xmg][0-9]*[KPkp]\|[KPkp][0-9]*[xmg][0-9]*[xmg][0-9]*" | tail -1 | sed 's|[xmkpKP]| |g' | sed 's|^[0]*||' | sed 's| 0| |g' | sed -e 's/^[ \t]*//' | sed 's/^[ \t]*//;s/[ \t]*$//'`;[ "`isnumber $kpoints_pfad`" != "yes" ] && echo kkk:$kpoints_pfad: && kpoints_pfad="???"

        #echo na:$atoms_tmpk:
        
        ct=XX;sc=XX;
        [ "`isnumber.sh $atoms_tmp`" != "yes" ] && atoms_tmp=$numatomspos
        [ "`getOption -d`" = "True" ] && echo atoms_tmp:$atoms_tmp
        if [ "`isnumber.sh $atoms_tmp`" = "yes" ];then
            fccatoms=(4 32 108 256 500 31 107 255 499)
            fccsc=(1 2 3 4 5 2 3 4 5)
            bccatoms=(2 16 54 128 249 15 53 127 248)
            bccsc=(1 2 3 4 5 2 3 4 5)
            sequence=`awk 'BEGIN {x=-1; while(++x<='"8"'){print x; }; exit}'`
            for s in `echo $sequence`;do
                f=`echo ${fccatoms[$s]}`
                b=`echo ${bccatoms[$s]}`
                [ "$f" = "$atoms_tmp" ] && atoms=$atoms_tmp && ct=fcc && sc=`echo ${fccsc[$s]}`
                [ "$b" = "$atoms_tmp" ] && atoms=$atoms_tmp && ct=bcc && sc=`echo ${bccsc[$s]}`
                #echo s:$s f:$f atoms:$atoms ct:$ct  sc:$sc
            done
            fi
        [ "`getOption -d`" = "True" ] && echo ct:$ct sc:$sc

        # vakrefFolder
        vakrefFolder=XX
        base=`pwd | sed 's|.*/||'`; basecheck=`find .. -type d -name "$base*ref*"`;hier=`pwd`;baseadd=""
        [ "`getOption -d`" = "True" ] && echo base: $base
        for i in $basecheck;do
        cd $hier;cd $i; [ "`pwd`" != "$hier" ] && baseadd="$baseadd `pwd`";cd $hier;
        done
        [ "`echo $baseadd | wc -w | sed 's|[ ]*||'`" = "1" ] && baseadd=`echo $baseadd | sed 's|[ ]*||'` && vakrefFolder=$baseadd 
        [ "`getOption -d`" = "True" ] && echo baseadd:$baseadd
        [ "`getOption -d`" = "True" ] && echo alats:$alats
        
        ## ngxf
        ngxf=xx
        ngxf_path=`pwd | sed -n 's|.*_\([0-9]*\)NGXF\([0-9]*\).*|\1\2|p'`;[ "`getOption -d`" = "True" ] && echo ngxf_path: $ngxf_path
        [ "`isnumber.sh $ngxf_path`" = "yes" ] && ngxf=$ngxf_path

        ## nbands
        nbands=xx
        INCARS=`find . -name INCAR`;
        if [ "$INCARS" != "" ];then
            for i in $INCARS;do
                nbands_tmp=`INCAR_NBANDS.sh $i`
                #echo $i
                [ "`isnumber.sh $nbands_tmp`" = "yes" ] && nbands=$nbands_tmp && break
            done
        fi
        
        
echo "aLats=`echo $alats | sed 's|,| |g'`     # Lattice constant in Angstrom" > parameters.dat
echo "disp=0.03      # displacement in bohrradius" >> parameters.dat
echo "cutoff=$cutoff_pfad     # in eV" >> parameters.dat
echo "" >> parameters.dat
echo "element=$element  # Al, al, Cu, Ti, ...(ignores case)" >> parameters.dat
echo "cellType=$ct   # fcc or bcc " >> parameters.dat
echo "sc=$sc           # 2,3,4,.... dimensions of supercell" >> parameters.dat
echo "nbands=$nbands     # 150  480/560  1200 2400   MP/Fermi" >> parameters.dat
echo "ngxf=$ngxf       # 120  180  240  300" >> parameters.dat
echo "kp=$kpoints_pfad           # 6    4    3    2" >> parameters.dat
echo "kpshift=0      # 0    0    .5   0" >> parameters.dat
echo "ediff=1E-8     # 1E-8 1E-8 1E-8 1E-8" >> parameters.dat
echo "" >> parameters.dat
echo "sigma=0.1      # smearing in eV; 0.1 sigma is 1160.4506 Kelvin" >> parameters.dat
echo "smear=-1       # smearing sceme (fermi/Methfessel-Paxton/...)" >> parameters.dat
echo "vakrefFolder=$vakrefFolder            # folder where the reference will be calculated using the same convergence parameters" >> parameters.dat
echo -e "\033[0m\033[1;4;31mparameters.dat wirtten .... please check input carefully!\033[0m"
echo
    fi

##########################################################
# xxx create INCAR template and exit
##########################################################
if [ "$debug" = "yes" ];then
    echo "#######################"
    echo wring INCAR
    echo "#######################"
fi
genPar=`getOption -i`
overwrite=`getOption -o`
if [ "$genPar" = "True" ] || [ "`getOption -a`" = "True" ] || [ ! -e "INCAR" ]; then
    if [ -e INCAR -a $overwrite != True ]; then error "INCAR existing; use -o to overwrite"; fi

echo "NPAR    =   1" > INCAR
echo "" >> INCAR
echo "NGXF    =   xxxNGXFxxx" >> INCAR
echo "NGYF    =   xxxNGXFxxx" >> INCAR
echo "NGZF    =   xxxNGXFxxx" >> INCAR
echo "ADDGRID =   .TRUE." >> INCAR
echo "" >> INCAR
echo "PREC    =   Accurate" >> INCAR
echo "ISTART  =   0             !WAVECAR:0-new,1-cont,2-samecut" >> INCAR
echo "ICHARG  =   2             !charge:1-file,2-overlapping_atom,10-const" >> INCAR
echo "" >> INCAR
echo "ENCUT   =   xxxCUTOFFxxx" >> INCAR
echo "NBANDS  =   xxxNBANDSxxx" >> INCAR
echo "LREAL   =   .FALSE." >> INCAR
echo "NELMDL  =   -5" >> INCAR
echo "EDIFF   =   xxxEDIFFxxx" >> INCAR
echo "ALGO    =   FAST" >> INCAR
echo "ISMEAR  =   xxxSMEARxxx" >> INCAR
echo "SIGMA   =   xxxSIGMAxxx" >> INCAR
echo "NELM    =   80" >> INCAR
echo "" >> INCAR
echo "LWAVE   =   F" >> INCAR
echo "LCHARG  =   F" >> INCAR
echo "LVTOT   =   F" >> INCAR
echo "LELF    =   F" >> INCAR
echo INCAR written
fi


########################################################################
# xxx create KPOINTS template and exit
########################################################################
genPar=`getOption -k`
overwrite=`getOption -o`
if [ "$genPar" = "True" ] || [ "`getOption -a`" = "True" ] || [ ! -e "KPOINTS" ]; then
    if [ -e KPOINTS -a $overwrite != True ]; then error "KPOINTS existing; use -o to overwrite"; fi
    echo "K-Points" > KPOINTS 
    echo " 0" >> KPOINTS
    echo "Monkhorst Pack" >> KPOINTS
    echo "xxxKPxxx xxxKPxxx xxxKPxxx" >> KPOINTS
    echo "0 0 0" >> KPOINTS
    echo "KPOINTS written"; 
    fi


########################################################################
# xxx read parameters.dat and other checks
########################################################################
echo checking syntax of parameters.dat ...
[ ! -e "parameters.dat" ] && error "please create first parameters.dat using `basename $0` -p" && exit
aLats=`get aLats`;
aLats=`echo "$aLats" | xargs -n1 | sort -n | sed 's|0*$||g'`  ## otherwise mathematica cant read 4.00/allForces.dat
disp=`get disp`;
cutoff=`get cutoff`;    ngxf=`get ngxf`; cellType=`get cellType`
element=`get element`;
kp=`get kp`;            nbands=`get nbands`;
sc=`get sc`; sigma=`get sigma`;      smear=`get smear`;
ediff=`get ediff`
disp=`echo $disp | awk '{printf("%.8f",$1*'$toAng')}'`
vakrefFolder=`get vakrefFolder`
if [ "`getOption -d`" = "True" ];then
    echo aLats:$aLats
    echo disp:$disp
    echo cutoff:$cutoff
    echo ngxf:$ngxf
    echo cellType:$cellType
    echo element:$element
    echo vakrefFolder:$vakrefFolder
    echo -----------------------------------------------
fi
check() {
  if [ ! -e $1 ]; then echo "" 1>&2; echo -e "\033[31m\033[1mERROR\033[0m: file $1 not existing" 1>&2; exit 1; fi
}



## it depends on the task which value we need to check ... but this we can implement later
[ "`isnumber $aLats`" != "yes" ] && error "aLats: $aLats : is not a number" && exit
[ "`isnumber $cutoff`" != "yes" ] && error "cutoff: $cutoff : is not a number" && exit
[ "`isnumber $kp`" != "yes" ] && error "kp: $kp : is not a number" && exit
[ "`isnumber $sc`" != "yes" ] && error "sc: $sc : is not a number" && exit
[ "`isnumber $disp`" != "yes" ] && error "disp: $disp : is not a number" && exit
[ "`isnumber $ngxf`" != "yes" ] && error "ngxf: $ngxf : is not a number" && exit
[ "`isnumber $nbands`" != "yes" ] && error "nbands: $nbands : is not a number" && exit
[ "$cellType" != "fcc" ] && [ "$cellType" != "bcc" ] && error "cellType is neither fcc nor bcc but it should be one of both; add e.g. in parameters.dat: cellType=fcc" && exit
if [ ! -e eqList_vacancy_$sc\x$sc\x$sc\sc ];then
    file=`find $path/utilities/$cellType -name eqList_vacancy_$sc\x$sc\x$sc\sc`
    [ "`getOption -d`" = "True" ] && echo ff: $file
    if [ -e $file ];then
        [ "`getOption -d`" = "True" ] && echo eqList exists
        cp $path/utilities/$cellType/eqList_vacancy_$sc\x$sc\x$sc\sc .
    fi
    fi
check eqList_vacancy_$sc\x$sc\x$sc\sc
d=`awk '{print $1}' eqList_vacancy_$sc\x$sc\x$sc\sc`
[ "`getOption -d`" = "True" ] && echo d:$d
[ ! -e "INCAR" ] && error "INCAR does not exist" && exit
[ ! -e "KPOINTS" ] && error "KPOINTS does not exist" && exit
[ ! -e "POTCAR" ] && error " POTCAR does not exist" && exit
numatoms=""
for aLat in $aLats; do
    [ ! -e relaxedCoords_$aLat ] && error "relaxedCoords_$aLat does not exist; either copy the relaxedCoords_xxx to this folder or define relaxedCoordFolder in parameters.dat" && exit
    add=`cat relaxedCoords_$aLat | wc -l | sed 's|[ ]*||'`
    
    [ "`getOption -d`" = "True" ] && echo add:$add:
    numatoms="$numatoms $add"
done
[ "`getOption -d`" = "True" ] && echo numatoms:$numatoms:
numatoms=`echo "$numatoms" | xargs -n1 | sort | uniq | sed 's|[ ]*||'`
[ "`getOption -d`" = "True" ] && echo numatoms:$numatoms:
[ "`echo $numatoms | wc -w | sed 's|[ ]*||'`" != "1" ] && echored "numatoms:$numatoms: not found,check your relaedCoords_xxx files"  && exit
refstr=""
[ "$numatoms" = "107" ] && refstr=$path/utilities/fcc/coordinates_3x3x3sc
[ "$numatoms" = "31" ] && refstr=$path/utilities/fcc/coordinates_2x2x2sc
[ "`getOption -d`" = "True" ] && echo refstr:$refstr:
[ ! -f $refstr ] && error "refstr: $refstr not found" && exit
refatoms=` expr $numatoms + 1 | sed 's|[ ]*||'`
[ "`getOption -d`" = "True" ] && echo refatoms:$refatoms:
tmelt=2000;tmelt_tmp=`getMeltingPoint.sh $element -r`;[ "`isnumber.sh $tmelt_tmp`" = "yes" ] && tmelt=$tmelt_tmp
[ "`getOption -d`" = "True" ] && echo tmelt:$tmelt
[ "$vakrefFolder" = "" ] && echo vakrefFolder does not exist, please put it in parameters.dat && exit
[ ! -e "$vakrefFolder"  ] && echo "vakrefFolder $vakrefFolder does not exist, please create it maually (mkdir $vakreffolder)" && exit
echo "                                  ... OK"
echo
checkAndSetMath


########################################################################
# xxx check if -alat defined
########################################################################
[ "`getOption -alat`" = "True" ] && aLats=`getValue -alat`
[ "`getOption -d`" = "True" ] && echo aLats:$aLats

########################################################################
# xxx create folders
########################################################################
if [ "`getOption -c`" = "True" ];then
echo creating folders ...
directory=`pwd`
rm -f jobList
touch jobList
echo "alat all: -> $aLats"
echo 

lom=`getlinuxormac`;[ "$lom" = "Linux" ] && add="";[ "$lom" = "mac" ] && add="'' -e"



# for every lattice constant
for aLat in $aLats; do
    cd $directory
    echogreen "alat -> $aLat"
    i=relaxedCoords_$aLat
    check $i   # check if relaxedCoords_$alat exists
    [ ! -e $aLat\Ang ] && mkdir $aLat\Ang  # heir nicht "$aLat\Ang" sondern $aLat\Ang !
    cd $aLat\Ang; 
    rm -f mapping.dat   
    scALat=`echo $aLat | awk '{printf("%.11f",$1*'$sc')}'`
    
    sed -e 's/xxxCUTOFFxxx/'$cutoff'/' \
        -e 's/xxxNBANDSxxx/'$nbands'/' \
        -e 's/xxxEDIFFxxx/'$ediff'/' \
        -e 's/xxxSIGMAxxx/'$sigma'/' \
        -e 's/xxxSMEARxxx/'$smear'/' \
        -e 's/xxxNGXFxxx/'$ngxf'/' $directory/INCAR > INCAR
    kp1=`echo $kp | awk '{print $1}'`
    kp2=`echo $kp | awk '{print $2}'`
    kp3=`echo $kp | awk '{print $3}'`
    sed -e 's/xxxKPxxx.*/'$kp1' '$kp2' '$kp3'/g' $directory/KPOINTS > KPOINTS

    echo "generic" > POSCAR
    echo "1" >> POSCAR
    echo " $scALat 0 0" >> POSCAR
    echo " 0 $scALat 0" >> POSCAR
    echo " 0 0 $scALat" >> POSCAR
    echo `cat $directory/$i | wc -l | sed 's|[ ]*||'` >> POSCAR
    echo "Cartesian" >> POSCAR
    
    
    ###################################################
    # background_forces folder
    ###################################################
    back=background_forces
    if [ ! -e $back ];then
        mkdir $back 
        cp POSCAR $back/
        cat $directory/$i >> $back/POSCAR
        cp $directory/POTCAR INCAR KPOINTS $back/
        sed -i $add 's|.*LWAVE.*|LWAVE=T|' $back/INCAR
        sed -i $add 's|.*LCHARG.*|LCHARG=T|' $back/INCAR
        echo `pwd`/background_forces >> $directory/jobList
    else
        echo $back exists
        fi
    
    ###################################################
    # reference_structure folder (braucht natuerlich auch eine auslenkung)
    ###################################################
    ref=$vakrefFolder/$aLat\Ang
    if [ ! -e $ref ];then
        echo CREATED vakrefFolder: $ref
        mkdir -p $ref
        cp $directory/POTCAR INCAR KPOINTS $ref
        #sed -i $add 's|.*LWAVE.*|LWAVE=T|' $ref/INCAR
        #sed -i $add 's|.*LCHARG.*|LCHARG=T|' $ref/INCAR
        echo "generic" > $ref/POSCAR
        echo "1" >> $ref/POSCAR
        echo " $scALat 0 0" >> $ref/POSCAR
        echo " 0 $scALat 0" >> $ref/POSCAR
        echo " 0 0 $scALat" >> $ref/POSCAR
        echo ` expr $numatoms + 1 ` >> $ref/POSCAR
        echo "Cartesian" >> $ref/POSCAR
        echo "$disp 0.0 0.0" >> $ref/POSCAR
        cat $refstr | awk '{printf(" %25.15f %25.15f %25.15f\n",$1*'"$scALat"'/'"$sc"',$2*'"$scALat"'/'"$sc"',$3*'"$scALat"'/'"$sc"')}' | tail -n+2 >> $ref/POSCAR
        ## auslenkung
        echo $ref >> $directory/jobList   
    else
        echo vakrefFolder $ref exists
    fi

    # for every displacement
    for j in $d; do
            at=`echo $j | awk '{printf("%d",($1-1)/3+1)}'`
            dir=`echo $j | awk '{printf("%d",(($1-1)%3)+1)}'`

            if [ ! -e $at\_$dir ];then
                echo $at\_$dir creating
                mkdir $at\_$dir
                cp POSCAR $at\_$dir/
                awk 'NR=='$at'&&'$dir'==1{printf(" %.15f  %s %s\n",$1+'$disp',$2,$3)} \
                     NR=='$at'&&'$dir'==2{printf(" %s  %.15f %s\n",$1,$2+'$disp',$3)} \
                     NR=='$at'&&'$dir'==3{printf(" %s  %s %.15f\n",$1,$2,$3+'$disp')} \
                     NR!='$at'{printf(" %s  %s %s\n",$1,$2,$3)}' $directory/$i >> $at\_$dir/POSCAR
                    cp $directory/POTCAR INCAR KPOINTS $at\_$dir/    # does this syntax work?
                    #ln -s ../$ref/WAVECAR $at\_$dir/WAVECAR
                echo `pwd`/$at\_$dir/ >> $directory/jobList
                else
                echo `pwd`/$at\_$dir exists
            fi


            if [ ! -e $at\_-$dir ];then
             echo $at\_-$dir creating
             mkdir $at\_-$dir
             cp POSCAR $at\_-$dir/
             awk 'NR=='$at'&&'$dir'==1{printf(" %.15f  %s %s\n",$1-'$disp',$2,$3)} \
                  NR=='$at'&&'$dir'==2{printf(" %s  %.15f %s\n",$1,$2-'$disp',$3)} \
                  NR=='$at'&&'$dir'==3{printf(" %s  %s %.15f\n",$1,$2,$3-'$disp')} \
                  NR!='$at'{printf(" %s  %s %s\n",$1,$2,$3)}' $directory/$i >> $at\_-$dir/POSCAR
             echo $j $at $dir >> mapping.dat
            cp $directory/POTCAR INCAR KPOINTS $at\_-$dir/    # does this syntax work?
                    #ln -s ../$ref/WAVECAR $at\_-$dir/WAVECAR
            echo `pwd`/$at\_-$dir/ >> $directory/jobList
        else
            echo `pwd`/$at\_-$dir exists
        fi
    done

    cd $directory
done

    echo creating folders ... done
fi


########################################################################
# xxx get Forces
########################################################################
if [ "`getOption -g`" = "True" ];then
		ineq=`awk '{print $1}' eqList_vacancy_$sc\x$sc\x$sc\sc`  # 1 3 7 8 19 22 23 34 35 (in 2x2x2 fcc sc)
directory=`pwd`
echo "get allForces.dat .... (also checking for necessary files,"
echo "if those are finished, if the amplitude of the displacement"
echo "is correct, ...."

l=`ls -1d *Ang`

for i in $l; do     # go over all lattice constants ..... but just take the one from parameters.dat (see 2 lines below, alatscheck)
    grepfor=`echo $i | sed 's|Ang| |'`
    alatscheck=`echo $aLats | xargs | sed 's|\(.*\)|\1 |'`
    checkifinpar=`echo "$alatscheck" | grep -o "$grepfor"`
    #echo alatscheck:"$alatscheck":
    #echo "checkifinpar:$checkifinpar:"
    #echo i:$i:check:$checkifinpar:
    [ "$checkifinpar" = "" ] && continue   ## only do files specified in parameters.dat
   echo folder: $i
   cd $directory
   cd $i
   lattconst=`echo $i | sed 's|Ang||'`
   rm -f allForces.dat allForces.dat.p allForces.dat.pn allForces.dat.pb allForces.dat.pnb
    


   ## check if background forces exist
   bgforces=""
   if [ -e "background_forces" ];then
        if [ "`OUTCAR_finished.sh background_forces`" = "yes" ];then
            bgforces=`OUTCAR_forces-last-ARRAY.sh background_forces`
            rm -f background_forces/forces
            echo "$bgforces" > background_forces/forces
        else
            echo;echo "!!!! background forces not finished !!!";echo
        fi
   else
       echo;echo ... no background forces found;echo
   fi
    [ "`getOption -d`" = "True" ] && echo bgforces: "$bgforces":
    [ "`getOption -d`" = "True" ] && echo ineq:$ineq:

	for j in $ineq; do   # ineq = 1 3 7 8 19 22 23 34 35 ( in fcc 2x2x2sc)
       ######################
       # get mapping dat
       ######################
        if [ ! -e "mapping.dat" ];then
        map=$path/utilities/$cellType/vacancy_mapping_$sc\x$sc\x$sc\sc
        [ -e $map ] && cp $map mapping.dat   # $map needs to be without ""
        fi
        [ ! -e "mapping.dat" ] && error "mapping.dat does not exist in `pwd`" && exit

      #####################################################
      # get positive displacements
      #####################################################
      at=`grep "^$j " mapping.dat | awk '{print $2}'`
      dir=`grep "^$j " mapping.dat | awk '{print $3}'`
      rm -f $at\_$dir/forces $at\_$dir/forcesminbg
      outcar=$directory/$i/$at\_$dir/OUTCAR.gz
      [ ! -e "$outcar" ] && gzip $directory/$i/$at\_$dir/OUTCAR
      #echo "    folder: "$at\_$dir""
      #echo "outcar:  $outcar"
      #echo "moutcar: $moutcar"
      [ ! -e "$outcar" ] && echo "   $directory/$i/$at\_$dir/OUTCAR.gz does not exist ... moving to next lattice constant" && break
      [ "`OUTCAR_finished.sh $outcar`" != "yes" ] && error "OUTCAR not finished in $at\_$dir/OUTCAR.gz" && exit
      at=`grep "^$j " mapping.dat | awk '{print $2}'`
      dir=`grep "^$j " mapping.dat | awk '{print $3}'`
      nions=`zgrep NIONS $at\_$dir/OUTCAR.gz | awk '{print $12}'`
      a=`zgrep NIONS $at\_$dir/OUTCAR.gz | awk '{print $12+1}'`
      zgrep -A$a POSITION $at\_$dir/OUTCAR.gz  | tail -n$nions | awk '{print $4,$5,$6}' > $at\_$dir/forces
      zgrep -A$a POSITION $at\_$dir/OUTCAR.gz  | tail -n$nions | awk '{print $1,$2,$3}' > $at\_$dir/position

      #########################################################
      # create forcesminbg if background exists (for positive displacements)
      #########################################################
      ## substract background forces
      take=forces
      [ "$bgforces" != "" ] && take=forcesminbg && paste $at\_$dir/forces background_forces/forces | awk '{printf "%.6f  %.6f  %.6f\n", $1-$4,$2-$5,$3-$6}' > $at\_$dir/forcesminbg 

      [ "`getOption -d`" = "True" ] && echo bgforces:$bgforces:
      [ "`getOption -d`" = "True" ] && echo take:$take

      ##########################################################
      # check if displacement is as big as stated in parameters.dat (for positive displacements)
      ##########################################################
      dispang=`paste $at\_$dir/position $directory/relaxedCoords_$lattconst | awk '{print sqrt(($1-$4)^2),sqrt(($2-$5)^2),sqrt(($3-$6)^2)}' | xargs -n1 | sort | uniq | awk 'max==""|| $1>max {max=$1}END{ print max}'`
      [ "`getOption -d`" = "True" ] && echo dispang: $dispang $disp
      checkdisp=`echo $dispang $disp | awk '{print sqrt(($1-$2)^2)}' | awk '$1 > 0.0001 {print "FALSE"}'`
      [ "$checkdisp" != "" ] && echo dispang in parameters.dat: $disp   disp in $at\_$dir: $dispang please check carefully && exit
      
      echo -e "$i/$at\e__$dir                  OK"
      
      
      #############################################################
      # get negative displacement
      #############################################################
      moutcar=$directory/$i/$at\_-$dir/OUTCAR.gz
      [ ! -e "$moutcar" ] && [ -e $directory/$i/$at\_-$dir/OUTCAR ] && gzip $directory/$i/$at\_-$dir/OUTCAR
      [ ! -e "$moutcar" ] && [ "$take" = "forces" ] && error "you need either to have background forces or done all displacements or have +- displacements" && exit
         if [ -e "$moutcar" ];then
              #echo " ...... negative disp found ... "
              # && echo "   $directory/$i/$at\_-$dir/OUTCAR.gz does not exist ... moving to next lattice constant" && break
              [ "`OUTCAR_finished.sh $moutcar`" != "yes" ] && error "OUTCAR not finished in $at\_$dir/OUTCAR.gz" && exit

              #echo "           "$at\_-$dir""
              rm -f $at\_-$dir/forces $at\_-$dir/forcesminbg
              zgrep -A$a POSITION $at\_-$dir/OUTCAR.gz | tail -n$nions | awk '{print $4,$5,$6}' > $at\_-$dir/forces

              #########################################################
              # create forcesminbg if background exists (for negative displacements)
              #########################################################
              [ "$bgforces" != "" ] && paste $at\_-$dir/forces background_forces/forces | awk '{printf "%.6f  %.6f  %.6f\n", $1-$4,$2-$5,$3-$6}' > $at\_-$dir/forcesminbg 
              zgrep -A$a POSITION $at\_-$dir/OUTCAR.gz | tail -n$nions | awk '{print $1,$2,$3}' > $at\_-$dir/position


              ##########################################################
              # check if displacement is as big as stated in parameters.dat (for negative displacements)
              ##########################################################
              dispang=`paste $at\_-$dir/position $directory/relaxedCoords_$lattconst | awk '{print sqrt(($1-$4)^2),sqrt(($2-$5)^2),sqrt(($3-$6)^2)}' | xargs -n1 | sort | uniq | awk 'max==""|| $1>max {max=$1}END{ print max}'`
              #[ "`echo $dispang | awk '$1>
              [ "`getOption -d`" = "True" ] && echo dispang: $dispang $disp
              checkdisp=`echo $dispang $disp | awk '{print sqrt(($1-$2)^2)}' | awk '$1 > 0.0001 {print "FALSE"}'`
              if [ "$checkdisp" != "" ];then   ## map positions back in supercell
                  scalat=`OUTCAR_cell-last-cartesian-ARRAY.sh $at\_-$dir/OUTCAR.gz | xargs -n1 | sort | uniq | awk 'max==""|| $1>max {max=$1}END{ print max}'`
                  [ "`getOption -d`" = "True" ] && echo dispang: $dispang disp:$disp scalat:$scalat
                  dispang=`echo 1 | awk '{print sqrt(('"$dispang"'-'"$scalat"')^2)}'`
                  checkdisp=`echo $dispang $disp | awk '{print sqrt(($1-$2)^2)}' | awk '$1 > 0.0001 {print "FALSE"}'`
                  [ "`getOption -d`" = "True" ] && echo dispppp: $dispang
              fi
              [ "$checkdisp" != "" ] && error "dispang in parameters.dat: $disp   disp in $at\_-$dir: $dispang please check carefully" && exit
              fi
              echo -e "$i/$at\e__-$dir                 OK"
        
      #######################################################################
      ## only positive forces (with and without background)
      #######################################################################
      if [ -e $at\_$dir/forces ];then
      # file is removed when cd $i
      cat $at\_$dir/forces | awk '{printf("%.7f %.7f %.7f\n",$1,$2,$3)}' | \
        xargs >> allForces.dat.p
      fi
      if [ -e $at\_$dir/forcesminbg ];then
      cat $at\_$dir/forcesminbg | awk '{printf("%.7f %.7f %.7f\n",$1,$2,$3)}' | \
        xargs >> allForces.dat.pb
      fi

      #######################################################################
      ## pos and neg forces (with and without background)
      #######################################################################
      if [ -e $at\_$dir/forces ] && [ -e $at\_-$dir/forces ];then
      # file is removed when cd $i
      paste $at\_{,-}$dir/forces | awk '{printf("%.7f %.7f %.7f\n",($1-$4)/2,($2-$5)/2,($3-$6)/2)}' | \
                xargs >> allForces.dat.pn
      fi
      if [ -e $at\_$dir/forcesminbg ] && [ -e $at\_-$dir/forcesminbg ];then
      paste $at\_{,-}$dir/forcesminbg | awk '{printf("%.7f %.7f %.7f\n",($1-$4)/2,($2-$5)/2,($3-$6)/2)}' | \
                xargs >> allForces.dat.pnb
      fi
        
    done
      files="allForces.dat.p allForces.dat.pb allForces.dat.pn allForces.dat.pnb"
      #echo celltype:$cellType:
      #echo numatoms:$numatoms:
      for ii in $files;do
          [ ! -e "$ii" ] && continue
          words=`cat $ii | wc -w | sed 's|[ ]*||'`
          [ "$cellType" = "fcc" ] && [ "$numatoms" = "107" ] && [ "$words" = "6099" ] && echo "                      $ii    OK" && continue  ## this is good!
          [ "$cellType" = "fcc" ] && [ "$numatoms" = "107" ] && [ "$words" != "6099" ] && echored "expected :6099: words in $ii and found only :$words: ... move $ii $ii.wrong" && mv $ii $ii.wrong && continue
          [ "$cellType" = "fcc" ] && [ "$numatoms" = "31" ] && [ "$words" = "837" ] && echo "                          $ii    OK" && continue  ## this is good!
          [ "$cellType" = "fcc" ] && [ "$numatoms" = "31" ] && [ "$words" != "837" ] && echored "expected :837: words in $ii and found only :$words: ... move $ii $ii.wrong" && mv $ii $ii.wrong && continue
            [ "$ii" = "allForces.dat.p" ] && echo "OK .... Background forces substracted from $ii"
            [ "$ii" = "allForces.dat.pb" ] && echo "OK .... Background forces substracted from $ii"
            [ "$ii" = "allForces.dat.pn" ] && echo "OK .... Background forces substracted from $ii"
            [ "$ii" = "allForces.dat.pnb" ] && echo "OK .... Background forces substracted from $ii"
        done
    cd $directory
done
fi


########################################################################
# xxx getSingleSpeciesPhonons.sh creates allForces.dat
########################################################################
files="allForces.dat.p allForces.dat.pb allForces.dat.pn allForces.dat.pnb"
if [ "`getOption -vac`" = "True" ];then
for i in $files;do  # files="allForces.dat.p allForces.dat.pb allForces.dat.pn allForces.dat.pnb"
    echo "############################################"
    echo "###  $i"
    echo "############################################"
  check=`find . -name $i | sed 's|\./||g'`
  [ "$check" = "" ] && continue
  folder=`echo "$check" | sed 's|'"$i"'||g'`
  for j in $folder;do
      echo cp $j/$i $j/allForces.dat
  done
  for j in $folder;do
      cp $j/$i $j/allForces.dat
  done
  #echo "$folder"
  #echo "$check"
  echo
    echo
    echo STARTING getSingleSpeciesPhonons.sh -vac -a $aLats
    echo
    getSingleSpeciesPhonons.sh -vac -a $aLats
    mkdir -p $i
    mv FreqsSupercell_* $i
    mv DynMatRSupercell_* $i
    mv FvibSupercell_* $i
    #mv FvibSupercell_perAtom_* $i
    mv HesseMatrix_* $i
    echo 
    #echo `pwd`
    echo
done
fi

########################################################################
# xxx fit surface vacancy ( here we already expect allForces.dat.xx folder to be available
########################################################################
folder="allForces.dat.p allForces.dat.pb allForces.dat.pn allForces.dat.pnb"
hier=`pwd`
if [ "`getOption -fit`" = "True" ];then
    checkAndSetMath
    overwrite=`getOption -o`
    echo numatoms:$numatoms
    echo refatoms:$refatoms
    evinet_d=EVinet_d_$numatoms
    evinet_b=EVinet_b_$refatoms
    #fqh_d=Fqh_d_$numatoms
    #fqh_b=$vakrefFolder/Fqh_b_$refatoms
    [ "`getOption -d`" = "True" ] && echo evinet_d:$evinet_d  evinet_b:$evinet_b
    [ ! -e "$evinet_d" ] || [ ! -e "$evinet_b" ] && error "please copy $evinet_d and $evinet_b to this folder" && continue
    #[ ! -e "$fqh_d" ] && error "please create $fqh_d first by invoking `basename $0` -fit" && exit
    #[ ! -e "$fqh_b" ] && error "please create $fqh_b first by invoking `basename $0` -fitref" && exit
    for i in $folder;do
        [ ! -e "$i" ] && continue
        cd $hier
        cd $i
          echo "############################################"
          echo "###  $i"
          echo "############################################"
          if [ -e fitFqh.input -a $overwrite != True ]; then error "fitFqh.input existing; use -o to overwrite"; fi
          allalats=`ls -1d FvibSupercell_[0-9.]* | sed 's|FvibSupercell_||' | sort -n | xargs | sed 's| |,|g'`
          allalats_=`ls -1d FvibSupercell_[0-9.]* | sed 's|FvibSupercell_||' | sort -n | xargs | sed 's| |_|g'`
          rm -f fitFqh.input
          echo
          echo "creating fitFqh.input ...."
          echo "                      .... OK"
          echo
          echo "aLats = {$allalats};" > fitFqh.input
          echo "structureFactor = 1;                                          (*  4: fcc  2: bcc  1: supercells *)" >> fitFqh.input
          echo "sc = $sc;                                                     (*  supercell, 1 for bulk *)" >> fitFqh.input
          echo "s = 1;                                                        (*  scaling factor for Fvib *)" >> fitFqh.input
          echo "fitOrder = 2;                                                 (*  typically: 2 giving: {1,V,V^2} *)" >> fitFqh.input
          #echo "baseNames = {\"Fvib_fromExactFreqs_"$numatoms"_\"};                   (*  Fvib_fromExactFreqs_  *)" >> fitFqh.input
          echo "baseNames = {\"FvibSupercell_\"};                             (*  Fvib_fromExactFreqs_  *)" >> fitFqh.input
          echo "<<\"$path/mathematica/fitFqh.MATH\"" >> fitFqh.input
          echo 
          echo "starting mathematica ........................."
          echo
          math < fitFqh.input
          echo "                     ......................... OK"
          echo ""
          mv FvibSupercell_fit_order2 FvibSupercell_fit_order2_$allalats_
          [ -e "gform" ] && echo "....removign gform old ..." && rm -rf gform
          mkdir gform
          cd gform
                echo get evinet_d: $evinet_d
                cp ../../$evinet_d .
                echo get evinet_b: $evinet_b
                cp ../../$evinet_b .
                echo get ../FvibSupercell_fit_order2_$allalats_
                cp ../FvibSupercell_fit_order2_$allalats_ Fqh_d_$numatoms  #this has to be the Fvib which was used fot the ti run

                ## general
                echo get Fqh_b_$refatoms
                FQH_B_REF=$vakrefFolder/Fqh_b_$refatoms
               [ ! -e "$FQH_B_REF" ] && echored "run this skript first for"
               [ ! -e "$FQH_B_REF" ] && echored "the reference, `basename $0`"
               [ -e "$FQH_B_REF" ] && cp $FQH_B_REF Fqh_b_$refatoms #this has to be the Fvib which was used fot the ti run



                ## Cu 3x3x3sc
#                cp /nas/glensk/v/pp/cu/dynmat_bulk_fcc4/3x3x3sc_290eV-NGXF288_4x4x4kP_own_TAKE/Fqh_whole_cell/Fvib_fromExactFreqs_108_fit_order2 Fqh_b_$refatoms
#                cp /nas/glensk/v/pp/cu/ti_bulk_fcc4/low_3x3x3sc_260eV_2x2x2kp_NGXF160_TAKE__high_400eV_4x4x4kpm0_NGXF120-TAKE/auswertung_lowplushighFahbest/Fah_fre_best_changedHesseref/Fah_surface_vakref_fit Fah_b_$refatoms
#                cp /nas/glensk/v/pp/cu/ti_vak_fcc4/low_3x3x3sc_260eV_2x2x2kp_240NGXF_TAKE__high_400eV_6x6x6kp_240NGXF_Vol-10_TAKE/auswertung_lowplushighFahbest/Fah_fre_best_changedHesseref/Fah_surface_fit Fah_d_$numatoms


                ## Al 3x3x3sc

                #FQH_B_REF="/nas/glensk/v/pp/al/ti_bulk_fcc4/Hessematrix_3x3x3sc/wholeCell/Fqh_fromExactFreqs_fit_order2"
                #[ -e "$FQH_B_REF" ] && cp $FQH_B_REF Fqh_b_$refatoms #this has to be the Fvib which was used fot the ti run
                #FAH_B_REF="/nas/glensk/v/pp/al/ti_bulk_fcc4/low_3x3x3sc_250eV_3x3x3kp_EDIFF1E-2__high_450eV_4x4x4kp__TAKE/auswertung_lowplushighFahbest/Fah_fre_best_changedHesseref/Fah_surface_fit"
               #[ -e "$FAH_B_REF" ] && cp $FAH_B_REF Fah_b_$refatoms
                #FAH_D="/nas/glensk/v/pp/al/ti_vak_fcc4/low_3x3x3sc_250eV_3x3x3kp_EDIFF1E-2_TAKE__high_450eV_4x4x4kp/auswertung_lowplushighFahbest/Fah_fre_best_changedHesseref/Fah_surface_fit"
               #[ -e "$FAH_D" ] && cp $FAH_D Fah_d_$numatoms 

                ## getGibbsEnergyOfFomation
                makegibbs="yes"
                [ "$makegibbs" = "yes" ] && [ ! -e "Fqh_b_$refatoms" ] && makegibbs="missing Fqh_b_$refatoms"
                [ "$makegibbs" = "yes" ] && [ ! -e "Fah_b_$refatoms" ] && makegibbs="missing Fah_b_$refatoms"
                [ "$makegibbs" = "yes" ] && [ ! -e "Fah_d_$numatoms" ] && makegibbs="missing Fah_d_$numatoms"
                if [ "$makegibbs" = "yes" ]
                    then getGibbsEnergyOfFormationLink.sh
                    else echored "CANNOT MAKE getGibbsEnergyOfFormation.sh
                        $makegibbs; did you make the fit for the vacancy
                        reference?"
                    fi

          cd $hier
    done
fi

########################################################################
# xxx fit surface reference vacancy
########################################################################
if [ "`getOption -fitref`" = "True" ];then
    [ ! -e "$vakrefFolder" ] && error "vakrefFolder: $vakrefFolder does not exist, please create reference jobs" && exit
    hier=`pwd`
    echo;echo "--> --> going to vakrefFolder: $vakrefFolder";echo
    cd $vakrefFolder
    overwrite=`getOption -o`
    if [ -e fitFqh.input -a $overwrite != True ]; then error "fitFqh.input existing in vakrefFolder `pwd`; use -o to overwrite"; fi
    for i in `find . -maxdepth 1 -mindepth 1 -name "Fqh_fromExactFreqs_[0-9.]*"`;do
      #echo $i
      [ "`cat $i | wc -l | sed 's|[ ]*||'`" = "0" ] && rm $i
    done  
    allref=`find . -maxdepth 1 -mindepth 1 -name "Fqh_fromExactFreqs_[0-9.]*" | sed 's|./||g' | sed 's|Fqh_fromExactFreqs_||g' | sort -n` 
    [ "`getOption -d`" = "True" ] && echo allref1: $allref
    [ ! -e "parameters.dat" ] && echo && echo "CREATING parameters.dat in `pwd` since it dit not exist" && echo && getSingleSpeciesPhonons.sh -p
    [ "$allref" = "" ] && getSingleSpeciesPhonons.sh -Fe -Fm
    allref=`find . -maxdepth 1 -mindepth 1 -name "Fqh_fromExactFreqs_[0-9.]*" | sed 's|./||g' | sed 's|Fqh_fromExactFreqs_||g' | sort -n` 
    [ "`getOption -d`" = "True" ] && echo allref2: $allref
    [ "$allref" = "" ] && error "did not find Fqh_fromExactfreqs_xxx in `pwd`; please run getSingleSpeciesPhonons.sh -Fe here manually" && eixt
    allalats=`echo "$allref" | xargs | sed 's| |,|g'`
    [ "`getOption -d`" = "True" ] && echo allalats: $allalats
    allalats_=`echo "$allref" | xargs | sed 's| |_|g'`
    [ "`getOption -d`" = "True" ] && echo allalats in reference folder: $allalats
    [ "$allalats" = "" ] && error "no Fqh_fromExactFreqs_xxx files found; please run first getSingelespeciesDynmat... in `pwd`" && exit
    rm -f fitFqh.input
    echo;echo "--> --> creating fitFqh.input ....";echo
    echo "aLats = {$allalats};" > fitFqh.input
    echo "structureFactor = 1;                                          (*  4: fcc  2: bcc  1: supercells *)" >> fitFqh.input
    echo "sc = $sc;                                                     (*  supercell, 1 for bulk *)" >> fitFqh.input
    echo "s = $refatoms;                                                        (*  scaling factor for Fvib *)" >> fitFqh.input
    echo "fitOrder = 2;                                                 (*  typically: 2 giving: {1,V,V^2} *)" >> fitFqh.input
    #echo "baseNames = {\"Fvib_fromExactFreqs_"$numatoms"_\"};                   (*  Fvib_fromExactFreqs_  *)" >> fitFqh.input
    echo "baseNames = {\"Fqh_fromExactFreqs_\"};                             (*  Fvib_fromExactFreqs_  *)" >> fitFqh.input
    echo "<<\"$path/mathematica/fitFqh.MATH\"" >> fitFqh.input
    math < fitFqh.input
    echo;echo " .... done";echo
    mv Fqh_fromExactFreqs_fit_order2 Fqh_fromExactFreqs_fit_order2_$allalats_
    [ -e "Fqh_b_$refatoms" ] && unlink Fqh_b_$refatoms
    [ "`getOption -d`" = "True" ] && echo allalats_: $allalats_  refatoms: $refatoms
    ln -s Fqh_fromExactFreqs_fit_order2_$allalats_ Fqh_b_$refatoms
    echo;echo done ...;echo
    cd $hier
fi

########################################################################
# xxx -t1 T-0K formation energy
########################################################################
if [ "`getOption -t1`" = "True" ];then
    evinet_d=EVinet_d_$numatoms
    evinet_b=EVinet_b_$refatoms
    [ "`getOption -d`" = "True" ] && echo evinet_d:$evinet_d  evinet_b:$evinet_b
    [ ! -e "$evinet_d" ] || [ ! -e "$evinet_b" ] && error "please copy $evinet_d and $evinet_b to this folder" && exit
    folder=nurEVinet_TAKE
    if [ ! -e $folder ];then
        mkdir $folder
        cp $evinet_d $folder
        cp $evinet_b $folder
        [ ! -e "$folder/temperature_range" ] && echo "1 $tmelt 2" > $folder/temperature_range
        hier=`pwd`
        cd $folder
        getGibbsEnergyOfFormationLink.sh
        cd $hier
    fi
fi

########################################################################
# xxx -t2 Fqh
########################################################################
if [ "`getOption -t2`" = "True" ];then
    evinet_d=EVinet_d_$numatoms
    evinet_b=EVinet_b_$refatoms
    fqh_d=Fqh_d_$numatoms
    fqh_b=$vakrefFolder/Fqh_b_$refatoms
    [ "`getOption -d`" = "True" ] && echo evinet_d:$evinet_d  evinet_b:$evinet_b
    [ ! -e "$evinet_d" ] || [ ! -e "$evinet_b" ] && error "please copy $evinet_d and $evinet_b to this folder" && exit
    [ ! -e "$fqh_d" ] && error "please create $fqh_d first by invoking `basename $0` -fit" && exit
    [ ! -e "$fqh_b" ] && error "please create $fqh_b first by invoking `basename $0` -fitref" && exit
    start=1
    while [ -e "quasiharmonic_$start" ];do
        start=` expr $start + 1 `
    done
    [ "`getOption -d`" = "True" ] && echo start: $start
    folder=quasiharmonic_$start
        mkdir $folder
        cp $evinet_d $folder
        cp $evinet_b $folder
        cp $fqh_d $folder
        cp $fqh_b $folder
        hier=`pwd`
        cd $folder
        getGibbsEnergyOfFormationLink.sh
        cd $hier
fi
