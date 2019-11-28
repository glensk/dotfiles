#!/bin/bash

# following 3 lines must always be present
path=`set | grep BASH_SOURCE | sed 's|.*\"\(.*\)/[^/]*\".*|\1|' | sed '/BASH_SOURCE=/s/.*/\./'`
script=`set | grep BASH_SOURCE | sed 's|.*\".*/\([^/]*\)\".*|\1|' | sed '/BASH_SOURCE=/s/.*/\./'`
options=$*; . $path/utilities/functions.include; checkOptions "-h -help -p -I -s -sp -si -K -a -c -ki -kpu -ff -f -e -efcc -epm -eh"

# small help for options
# space between option and value must be 1
# space between option and Comment or value and Comment must be >2
h=`getOption -h`
if [ $h == True -o $# == 0 ]; then
  usage $script
  printOptions "-p               create parameters.dat and exit" \
               "-I               create INCAR template and exit" \
               "-K               create KPOINTS template and exit" \
               "-a               create create all of the above and exit" \
               "-c               create Folders and exit" \
               "-ki              create Folders and keep INCAR as it is for all jobs" \
               "-kpu             create Folders and keep POSCAR with undisplaced atoms (for reference structrue) for all jobs" \
               "-ff              force running even if EDIFF>1e-8" \
               "-f               force template or folder creation even if existing" \
               "-e               go in every dynmat folder and run VASP_sphinx_get_sxdynmatsx.sh" \
               "-efcc            go in every dynmat folder and run VASP_sphinx_get_sxdynmatsx.sh -fcc" \
               "-eh              go in every dynmat folder and run sxdynmat -i sxdynmat.sx & -H" \
               "-epm             go in every dynmat folder and create HesseMatrix (use masses from POTCAR instead of OUTCAR)" \
               "-s [murnfolder]  create simply from all POSCAR_xxx jobs; murnfolder can be specified to get corresponding (relaxed) CONTCARS" \
               "-si              increase the range of taken CONTCARS by 0.6% (this may include one more structures close to the minimum ot the cold curve)" \
               "-sp              to be used with -s option. From the murnfolder given instead of the CONTCAR the POSCAR is taken"
  exit
fi

# detailed help
help=`getOption -help`
if [ $help == True ]; then
  details $script
  echo2 "   creates input folders for a dynamical matrix calculation"
  echo2 "   parameters.dat file must be present (generate with -p)"
  echo2 "   INCAR KPOINTS templates must be present" \
        "   (generate with -I -K or -a (including -p))"
  echo2 "   typical script sequence:"       \
        "     \033[1mcreateFolders_dynMat.sh\033[0m -a" \
        "     \033[1mcreateFolders_dynMat.sh\033[0m"    \
        "     collectOUTCARs.sh"               \
        "     extractForces.sh"                \
        "     getSingleSpeciesPhonons.sh -A                  (cmpc01)"  \
        "     getFqhFit.sh                                   (cmpc01)"  \
        "     getThermodynamics.sh               (in separate folder)"
  exit
fi

if [[ "`getOption -e`" = "True" || "`getOption -epm`" = "True" || "`getOption -efcc`" = "True" ]];then
  folder=`ls -1d [0-9]*`
  [ "`echo $folder | wc -w`" = "0" ] && folder=`ls -1d Vol[0-9]*`
  hier=`pwd`
  for i in $folder;do
    cd $hier
    cd $i
    echo 
    echo 
    echo "##################################################################################"
    echo `pwd`
    echo "##################################################################################"
    [ "`getOption -epm`" = "True" ] && VASP_sphinx_get_sxdynmatsx.sh -pm
    [ "`getOption -efcc`" = "True" ] && VASP_sphinx_get_sxdynmatsx.sh -fcc
    [ "`getOption -e`" = "True" ] && VASP_sphinx_get_sxdynmatsx.sh

    cd $hier
  done
fi

if [ "`getOption -eh`" = "True" ];then
  folder=`ls -1d [0-9]*`
  hier=`pwd`
  for i in $folder;do
    cd $hier
    echo $i
    cd $i
    sxdynmat -i sxdynmat.sx
    sxdynmat -i sxdynmat.sx -H
    vol=`POSCAR_volume.sh`
    echo vol:$vol

    ## get elements
    elements=`POTCAR_element.sh`
    elements_anz=`echo $elements | wc -w`
    echo $elements,$elements_anz
    
    atoms_frompos=`POSCAR_numatomsline.sh`
    atoms_frompos_anz=`echo $atoms_frompos | wc -w`
    [ "$elements_anz" != "$atoms_frompos_anz" ] && echo aa: $elements bb: $atoms_frompos && exit
    echo $atoms_frompos
    
    [ "$elements_anz" = "1" ] && forhesse=$elements
    [ "$elements_anz" = "2" ] && forhesse=`echo $elements $atoms_frompos | awk '{print $1,$3,$2,$4}'`
    echo forhesse: $forhesse

    hesse.py $forhesse -we -wf
    file=`ls -1d Fqh_fromExactFreqs_*`
    echo $file && cp $file ../$file\_$vol
    cp HesseMatrix_sphinx ../HesseMatrix_$vol

    cd $hier
  done
fi

[ "`getOption -sp`" = "True" ] && [ "`getOption -s`" != "True" ] && echo -sp option has to be used with -s option && exit
if [ "`getOption -s`" = "True" ];then
    echo "getOption -s:"`getOption -s`
    echo "getValue  -s:"`getValue -s`
    if [ "`getValue -s`" != "" ];then
        murnpath=`getValue -s`
        [ ! -e "$murnpath" ] && echo murnpath $murnpath does not exist && exit
        echo murnpath: $murnpath

        contcars=`find -L $murnpath -name "CONTCAR*"`
        [ "`getOption -sp`" = "True" ] && contcars=`find -L $murnpath -name "POSCAR"`
        atoms=0
        minvol=0
        vol=0
        for i in $contcars;do
            echo ""
            echo ""
            echo i:$i
            contcarpath=`echo $i | sed 's|CONTCAR.*||'`
            contcarpathname=`echo $i | sed 's|.*/CONTCAR|CONTCAR|'`
            echo cpn: $contcarpathname
            echo contcarpath: $contcarpath
            [ "$atoms" = "0" ] && [ "`getOption -sp`" != "True" ] && atoms=`OUTCAR_number_of_atoms.sh $contcarpath`
            [ "$atoms" = "0" ] && [ "`getOption -sp`" = "True" ] && atoms=`POSCAR_numatoms.sh $contcarpath`
            [ "`isnumber.sh $atoms`" != "yes" ] && atoms=`OUTCAR_number_of_atoms.sh $contcarpath`
            echo atoms1: $atoms
            [ "`isnumber.sh $atoms`" != "yes" ] && atoms=`OUTCAR_number_of_atoms.sh $i`
            echo atoms2: $atoms
            [ "`isnumber.sh $atoms`" != "yes" ] && atoms=`POSCAR_numatoms.sh $contcarpath`

            echo atoms3: $atoms
            volnecessary="OK"
            [ -e "$murnpath/EVinet_$atoms" ] && minvol=`head -1 $murnpath/EVinet_$atoms | awk '{print $2}'`
            [ -e "$murnpath/EVinet_$atoms" ] && [ "`getOption -si`" = "True" ] && minvol=`head -1 $murnpath/EVinet_$atoms | awk '{print $2-$2*0.006}'`
            [ "`getOption -s`" = "True" ] && [ "`getOption -sp`" != "True" ] && vol=`OUTCAR_volume-lastexact.sh $contcarpath`
            [ "`getOption -s`" = "True" ] && [ "`getOption -sp`" =  "True" ] && vol=`POSCAR_volume.sh $contcarpath`
           
            echo minvol:$minvol:
            echo vol:$vol:
            char1=`echo $vol | cut -c 1`
            charcheck=`echo $vol | grep "[a-zA-Z]"`
            [ "$vol" = "" ] && echo contcarpath $contcarpath has no volume && volnecessary="NO" && continue
            [ "$char1" = "$charcheck" ] && echo contcarpath $contcarpath has no volume && volnecessary="NO" && continue
            [ "`echo "$vol" | wc -w`" != "1" ] && volnecessary="NO"
            

            [ "$minvol" != "0" ] && volnecessary=`echo $vol $minvol | awk '$1 > $2 {print "OK"}'`
            #echo contcarpath:$contcarpath $vol:$vol volnecessary:$volnecessary i:$i
            if [ "$volnecessary" = "OK" ];then
                if [ "$contcarpathname" = "CONTCAR.gz" ];then
                    cp $i POSCAR_$vol.gz
                    gunzip POSCAR_$vol.gz
                else
                    cp $i POSCAR_$vol
                fi
                echo $vol necessary
            else
                echo $vol not NECESSARY
            fi
        done
    echo $murnpath >> ANMERKUNG
    fi
    vols=`ls -1d POSCAR_* | sed 's|POSCAR_||' | xargs`
    hier=`pwd`
    rm -f jobList
    for i in $vols;do
        cd $hier
        [ -e "$i" ] && echo $i already exists && continue
        [ ! -e "$i" ] && mkdir $i 
        cd $i
        echo $i
        gvinp ../
        VASP_sphinx_create_all_dynmat_jobs.sh
        cat jobList >> ../jobList

        cd $hier
    done
    echo "CHECK NBANDS in your INCAR!!"
    exit
fi

# if -a option we create all templates and exit
all=`getOption -a`

# check if force overwriting templates or folders
overwrite=`getOption -f`

# if applies create parameters.dat
genPar=`getOption -p`
if [ $genPar == True -o $all == True ]; then
  if [ -e parameters.dat -a $overwrite != True ]; then error "parameters.dat existing; use -f to overwrite"; fi
  echo "
aLats= 4.               # Angstrom
cutoff= 220             # eV

type= fcc               # fcc/bcc/POSCAR; for POSCAR the given POSCAR template is used
                        # supercell sensitive values, typical values:
sc= 2                   # 2sc      3sc      4sc       (note: these are only guide values)
ngxf= 120               # 120      180      240       (might need higher settings for very high cutoffs >400eV??)
kp= 4                   # 12       8        6         for a well converged bcc and very well converged fcc calculation
nbands= 100             # 150/180  480/580  1100/xxx  MP/Fermi with 0.1eV for fcc Ca with 8 valence electrons

                        # supercell independent values:
disp=0.03               # displacement in Bohr (typically 0.03 is ok)
smear=1                 # 1=Methfessel, -1=Fermi
sigma=0.1               # smearing in eV
kpshift=0               # typically 0 or 0.5

addString=              # if addINCAR is used it is useful to add some string to name; leave blank otherwise
addINCAR=               # add additional flags to INCAR, e.g., addINCAR=ISPIN=2; leave blank if not needed
                        # more flags can be added by separating them with ;
" > parameters.dat
  echo; echo "parameters.dat written  <==  mind: parameters.dat needs to be adjusted"
fi

# if applies create INCAR
genINC=`getOption -I`
if [ $genINC == True -o $all == True ]; then
  if [ -e INCAR -a $overwrite != True ]; then error "INCAR existing; use -f to overwrite"; fi
  if [[ `hostname` == cmmc* ]]; then NPAR=4; else NPAR=1; fi
  echo "
 NPAR = $NPAR    ! 1 for cmmd and 4 for cmmc

 NGXF    =   xxxNGXFxxx
 NGYF    =   xxxNGXFxxx
 NGZF    =   xxxNGXFxxx
 ADDGRID =   .TRUE.

 ENCUT  =    xxxCUTOFFxxx
 ISMEAR =    xxxSMEARxxx
 SIGMA  =    xxxSIGMAxxx
 NBANDS =    xxxNBANDSxxx

 PREC   =    Accurate
 LREAL  =    .FALSE.
 ALGO   =    NORMAL
 EDIFF  =    1E-8
 NELM   =    200

 LWAVE  =    F   
 LCHARG =    F   
" > INCAR
  echo; echo "INCAR written"
fi

# if applies create KPOINTS
genKPO=`getOption -K`
if [ $genKPO == True -o $all == True ]; then
  if [ -e KPOINTS -a $overwrite != True ]; then error "KPOINTS existing; use -f to overwrite"; fi
  echo "K-Points
 0
Monkhorst Pack
 xxxKPxxx xxxKPxxx xxxKPxxx
 xxxKPSHIFTxxx xxxKPSHIFTxxx xxxKPSHIFTxxx
" > KPOINTS
  echo; echo "KPOINTS written"
fi

# exit if we created some template
if [ $genPar == True -o $genINC == True -o $genKPO == True -o $all == True ]; then exit; fi

# create folders
if [ "`getOption -c`" = "True" ] || [ "`getOption -ki`" = "True" ];then
    # check if all input files available
    input="parameters.dat INCAR KPOINTS POTCAR"
    for i in $input; do check $i; done
    
    # read in all parameters from parameters.dat
    aLats=`get aLats`; cutoff=`get cutoff`; type=`get type`; sc=`get sc` kp=`get kp`; kpshift=`get kpshift`; addString=`get addString`
    smear=`get smear`; sigma=`get sigma`; nbands=`get nbands`; addINCAR=`get addINCAR`; ngxf=`get ngxf`; disp=`get disp`
    checkInput "$aLats" "$cutoff" "$type" "$sc" "$kp" "$kpshift" "$smear" "$sigma" "$nbands" "$ngxf" "$disp"
    
    if [ $smear == -1 ]; then s=$sigma; else s=0.0; fi
    
    # check if structure type known
    if [ "$type" != "bcc" -a "$type" != "fcc" -a "$type" != "POSCAR" ]; then
      error "structure type in parameters.dat not known"
    fi

    # if type is not POSCAR we need to create a fcc or bcc POSCAR according to supercell size sc in parameters.dat
    if [ "$type" != "POSCAR" ]; then
      c=`checkInteger $sc`
      if [ $c != ok ]; then error "supercell size non integer"; fi 
      if [ "$sc" -lt 1 -o "$sc" -gt 6 ]; then error "supercell size not supported"; fi
      if [ -e POSCAR -a $overwrite != True ]; then error "POSCAR existing; use -f to overwrite"; fi

      # generate POSCAR
      natoms=`awk 'BEGIN{c=0}; NF==3{c=c+1}; END{print c}' $path/utilities/$type/coordinates_$sc\x$sc\x$sc\sc`
      echo "DynMat $type ${sc}x${sc}x${sc}sc
xxxALATxxx
$sc 0 0
0 $sc 0
0 0 $sc
$natoms
Cartesian
xxxDISPxxx   0.0   0.0" > POSCAR
      awk 'NR>1{print}' $path/utilities/$type/coordinates_${sc}x${sc}x${sc}sc >> POSCAR
    fi
    
    # check if templates ok
    if [ "`getOption -ki`" != "True" ];then
      cc=`grep -e xxxCUTOFFxxx -e xxxNGXFxxx -e xxxSMEARxxx -e xxxSIGMAxxx -e xxxNBANDSxxx INCAR | wc -l`
      if [ $cc != 7 ]; then error "INCAR template wrong"; fi
      cc=`grep -e xxxKPxxx -e xxxKPSHIFTxxx KPOINTS | xargs -n1 | grep -e xxxKPxxx -e xxxKPSHIFTxxx | wc -l`
      if [ $cc != 6 ]; then error "KPOINTS template wrong"; fi
      cc=`grep -e xxxALATxxx -e xxxDISPxxx POSCAR | wc -l`
      if [ $cc != 2 ]; then error "POSCAR template wrong"; fi
    fi
    
    # check EDIFF
    ediffForce=`getOption -ff`
    cc=`sed -nr 's/^[ ]*EDIFF[ ]*=([^;!]*).*/\1/p' INCAR | awk '{if ($1>1e-8) print "larger"}'`
    cc_i=`sed -nr 's/^[ ]*EDIFF[ ]*=([^;!]*).*/\1/p' INCAR`
    if [ "$cc" == larger -a $ediffForce != True ]; then error "EDIFF>1e-8 ($cc, $cc_i); run with -ff option to override"; fi
    
    # bohrradius to angstrom
    toAng=0.529177208278835
    
    # check if folder exist and should be overwritten
    for a in $aLats; do
      f=$sc\x$sc\x$sc\sc_$cutoff\eV-NGXF$ngxf\-ADDGRID_$kp\x$kp\x$kp\kp-$kpshift\x$kpshift\x$kpshift\shift_$s\eV$addString/$a\Ang
      if [ -e "$f" -a "$overwrite" != "True" ]; then
        error "folder $f exists; run with -f option to force overwriting";
      fi
    done
    
    # create the input folders
    echo; rm -f jobList
    for a in $aLats; do
      f=$sc\x$sc\x$sc\sc_$cutoff\eV-NGXF$ngxf\-ADDGRID_$kp\x$kp\x$kp\kp-$kpshift\x$kpshift\x$kpshift\shift_$s\eV$addString/$a\Ang
      echo $f; rm -fr $f; mkdir -p $f
    
      # transform displacement to angstrom
      d=`echo "$disp * $toAng / $a" | bc -l`
  
      ##################################
      # INCAR
      ##################################
      if [ "`getOption -ki`" = "True" ];then
          ENCUT=`INCAR_ENCUT.sh INCAR`
          ENCUT_check=`echo $cutoff $ENCUT | awk '{print $1-$2}'`
        #echo kk `getOption -ki`
        # if keep incar dont change ENCUT
          #[ "`getOption -ki`" = "False" ] && [ "$ENCUT_check" != "0" ] && echo "changing INCAR ENCUT to parameters.dat value" && INCAR_change.sh INCAR ENCUT $cutoff
          #echo "cutoff parameters: $cutoff INCAR: $ENCUT check:$ENCUT_check" && exit
          cp INCAR $f/INCAR
      else
          sed -e 's/xxxCUTOFFxxx/'$cutoff'/' \
              -e 's/xxxNGXFxxx/'$ngxf'/'     \
              -e 's/xxxSMEARxxx/'$smear'/'   \
              -e 's/xxxSIGMAxxx/'$sigma'/'   \
              -e 's/xxxNBANDSxxx/'$nbands'/' INCAR > $f/INCAR
          if [ -n "$addINCAR" ]; then
            echo "$addINCAR" | xargs -n1 -d \; >> $f/INCAR;
          fi
      fi
      ##################################
      # KPOINTS
      ##################################
      sed -e 's/xxxKPxxx/'$kp'/g' \
          -e 's/xxxKPSHIFTxxx/'$kpshift'/g' KPOINTS > $f/KPOINTS

      ##################################
      # POSCAR
      ##################################
      if [ "`getOption -kpu`" = "True" ];then
          sed -e 's/xxxALATxxx/'$a'/' \
              -e 's/xxxDISPxxx/0.0/' POSCAR > $f/POSCAR
      else
          sed -e 's/xxxALATxxx/'$a'/' \
              -e 's/xxxDISPxxx/'$d'/' POSCAR > $f/POSCAR
      fi
      cp POTCAR $f/POTCAR
      echo `pwd`/$f >> jobList
    done
    
    echo; echo " jobList file written"
fi

