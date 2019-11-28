#!/bin/bash

# following 3 lines must always be present
path=`set | grep BASH_SOURCE | sed 's|.*\"\(.*\)/[^/]*\".*|\1|' | sed '/BASH_SOURCE=/s/.*/\./'`
script=`set | grep BASH_SOURCE | sed 's|.*\".*/\([^/]*\)\".*|\1|' | sed '/BASH_SOURCE=/s/.*/\./'`
options=$*; . $path/utilities/functions.include; checkOptions "-h -help -p -I -K -P -a -k -ff -f -l -c -cp"

# small help for options
# space between option and value must be 1
# space between option and Comment or value and Comment must be >2
h=`getOption -h`
if [ $h == True ]; then
  usage $script
  printOptions "-p         create parameters.dat and exit" \
               "-I         create INCAR template and exit" \
               "-K         create KPOINTS template and exit" \
               "-P type    create POSCAR template of type (fcc/bcc/hcp/dhcp/sc) and exit" \
               "-a type    create create all of the above and exit" \
               "-k kpc     override automatic setting of kpoint along c by kpc (hcp only)" \
               "-ff        force running even if EDIFF>1e-5" \
               "-f         force template or folder creation even if existing" \
               "-c [path]  find all CONTCARS(.gz) in path and copy here to POSCAR_volume" \
               "-cp        just loop over POSCAR_volume folder and create folder" \
               "-l         just loop over aLats= from parameter.dat and create folder"
  exit
fi

# detailed help
help=`getOption -help`
if [ $help == True ]; then
  details $script
  echo2 "   creates input folders for a cold curve (T=0K) run"
  echo2 "   parameters.dat file must be present" \
        "   (generate with -p)"
  echo2 "   INCAR KPOINTS POSCAR templates must be present" \
        "   (generate with -I -K -P bcc or fcc or hcp or dhcp or sc)"
  echo2 "   typical script sequence:"       \
        "     \033[1mcreateFolders_murn.sh\033[0m -a ..." \
        "     \033[1mcreateFolders_murn.sh\033[0m"        \
        "     collectOUTCARs.sh"               \
        "     extractEnergies.sh"              \
        "     getE0KFit.sh           (cmpc01)" \
        "     getThermodynamics.sh"
  exit
fi

# get CONTCARS
if [ "`getOption -c`" = "True" ];then
    murnpath=`getValue -c`
    [ ! -e "$murnpath" ] && echo murnpath $murnpath does not exist && exit
    contcars=`find -L $murnpath -name "CONTCAR*"`
     for i in $contcars;do
         echo "i: $i"
         contcarpath=`echo $i | sed 's|CONTCAR.*||'`
         contcarpathname=`echo $i | sed 's|.*/CONTCAR|CONTCAR|'`
         echo cp: $contcarpath
         #vol=`POSCAR_volume.sh $contcarpath`
         vol=`POSCAR_volume.sh $i`
         echo i:$i vol:$vol
         [ "`echo $vol | wc -w`" != "1" ] && echo $contcarpath no volume recognized .. continue && continue
         echo "CONTCAR_$vol -> POSCAR_$vol"
         if [ "$contcarpathname" = "CONTCAR.gz" ];then
            cp $i POSCAR_$vol.gz
            gunzip POSCAR_$vol.gz
         else
            cp $i POSCAR_$vol
         fi
         sed -i 's|^  0.00000000E+00  0.00000000E+00  0.00000000E+00||g' POSCAR_$vol
         #sed -i '6,6d' POSCAR_$vol
         sed -i '/^ *$/d' POSCAR_$vol
     done
fi

## get POTCARS
#if [ "`getOption -c`" = "True" ];then
#    murnpath=`getValue -c`
#    [ ! -e "$murnpath" ] && echo murnpath $murnpath does not exist && exit
#    contcars=`find -L $murnpath -name "CONTCAR*"`
#     for i in $contcars;do
#         contcarpath=`echo $i | sed 's|CONTCAR.*||'`
#         vol=`POSCAR_volume.sh $contcarpath`
#         [ "`echo $vol | wc -w`" != "1" ] && echo $contcarpath no volume recognized .. continue && continue
#         echo "CONTCAR_$vol -> POSCAR_$vol"
#         cp $i POSCAR_$vol
#     done

# loop through POSCARS
if [ "`getOption -cp`" = "True" ];then
    [ ! -e "POTCAR" ] && echo POTCAR does not exist && exit
    [ ! -e "KPOINTS" ] && echo KPOINTS does not exist && exit
    [ ! -e "INCAR" ] && echo INCAR does not exist && exit
    rm -f jobList
    for a in `ls -1d POSCAR_*`;do
        vol=`echo $a | sed 's|POSCAR_||'`
        echo $a $vol
        [ -e "$vol" ] && echo $vol exists && continue
        mkdir $vol
        cp POTCAR $vol
        cp $a $vol/POSCAR
        cp INCAR $vol
        cp KPOINTS $vol
        echo `pwd`/$vol >> jobList
        done
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
aLats=3.22 3.24 3.26    # Angstrom, example bcc Ti
cutoff=200              # eV; more values possible (separated with space) then a loop is performed
kp=8                    # gives NxNxN for fcc/bcc and NxNxInt(N/cBya); also here more values are possible for looping
smear=1                 # 1=Methfessel, -1=Fermi
sigma=0.1               # smearing in eV
nbands=AUTO             # if AUTO then NBANDS is left out of INCAR so vasp determines it automatically (from NELEC)
cBya=1                  # 1 for bcc/fcc/sc and c/a ratio for hcp or dhcp (ideal=1.632993161855452)

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
  echo "
 NPAR = 1

 ENCUT  = xxxCUTOFFxxx
 ISMEAR = xxxSMEARxxx
 SIGMA  = xxxSIGMAxxx
 EDIFF  =   1E-5

 ADDGRID=    TRUE
 PREC   =    Accurate
 NBANDS =    xxxNBANDSxxx
 ALGO   =    NORMAL

 LWAVE  =      F   
 LCHARG =      F   
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
 xxxKPxxx xxxKPxxx xxxKPCxxx
 0 0 0
" > KPOINTS
  echo; echo "KPOINTS written"
fi

# copy POSCAR
genPOS=`getOption -P`
if [ $genPOS == True -o $all == True ]; then
  if [ -e POSCAR -a $overwrite != True ]; then error "POSCAR existing; use -f to overwrite"; fi
  if [ $genPOS == True ]; then type=`getValue -P`; fi
  if [ $all == True ]; then type=`getValue -a`; fi
  if [ -z "$type" ]; then error "type missing"; fi

  case "$type" in
    sc) echo "Murn sc
xxxALATxxx
1 0 0
0 1 0
0 0 1
1
Cartesian
0 0 0" > POSCAR;;

    bcc) echo "Murn bcc
xxxALATxxx
1 0 0
0 1 0
0 0 1
2
Cartesian
0 0 0
.5 .5 .5" > POSCAR;;

   fcc) echo "Murn fcc
xxxALATxxx
1 0 0
0 1 0
0 0 1
4
Cartesian
0 0 0
0 .5 .5
.5 0 .5
.5 .5 0" > POSCAR;;

    hcp) sin=0.866025403784438597 # sin(pi/3) where pi/3 is angle between hcp lattice vectors
         echo "Murn hcp
xxxALATxxx
1.0  0.0  0.0
0.5  $sin  0.0
0.0  0.0  xxxCBYAxxx
2
Direct
0.000000000000 0.000000000000 0.0
0.333333333333 0.333333333333 0.5
    " > POSCAR;;

    dhcp) sin=0.866025403784438597 # sin(pi/3) where pi/3 is angle between hcp lattice vectors
         echo "Murn dhcp
xxxALATxxx
1.0  0.0  0.0
0.5  $sin  0.0
0.0  0.0  xxxCBYAxxx
4
Direct
0.000000000000 0.000000000000 0.0
0.333333333333 0.333333333333 0.25
0.000000000000 0.000000000000 0.5
0.666666666667 0.666666666667 0.75
    " > POSCAR;;

    *) "structure type $type not known";;
  esac



  echo; echo "POSCAR written"
fi

# exit if we created some template
if [ $genPar == True -o $genINC == True -o $genKPO == True -o $genPOS == True -o $all == True ]; then exit; fi

# check if all input files available
input="parameters.dat INCAR KPOINTS POTCAR POSCAR"
for i in $input; do check $i; done


# if just loop:
if [ "`getOption -l`" = "True" ];then
    aLats=`get aLats`
    rm -f jobList
    for a in $aLats;do
        echo $a
        mkdir $a
        cp POTCAR $a
        cp POSCAR $a
        cp INCAR $a
        cp KPOINTS $a
        #fillin=`echo "$a" | awk '{printf "%.16f", ($1*64.)^(1./3.)}'`
        #fillin=`echo "$a" | awk '{printf "%.16f", $1*3}'`
        fillin=`echo "$a" | awk '{printf "%.16f", $1}'`
        sed -i 's|xxxALATxxx|'"$fillin"'|' $a/POSCAR
        echo `pwd`/$a >> jobList
        done
    exit
fi




# read in all parameters from parameters.dat
aLats=`get aLats`; cutoffs=`get cutoff`; kps=`get kp`; addString=`get addString`; cBya=`get cBya`
smear=`get smear`; sigma=`get sigma`; nbands=`get nbands`; addINCAR=`get addINCAR`
checkInput "$aLats" "$cutoffs" "$kps" "$smear" "$sigma" "$nbands" "$cBya"

if [ $smear == -1 ]; then s=$sigma; else s=0.0; fi

# check if structure type known
type=`awk 'NR==1{print $2}' POSCAR`
echo type $type
if [ "$type" != "sc" -a "$type" != "bcc" -a "$type" != "fcc" -a "$type" != "hcp" -a "$type" != "dhcp" ]; then
  error "structure type in POSCAR not known";
fi

# if dhcp, scale c axis (which now equals cBya) times 2
if [ "$type" == "dhcp" ]; then
  cBya=`echo $cBya | awk '{printf "%.10f",$cBya*2}'`
fi

# loop over k-points and cutoff
for kp in $kps; do
for cutoff in $cutoffs; do

  # scale KPC if hcp or dhcp
  if [ $type == hcp -o $type == dhcp ]; then
    kpc=`echo $kp $cBya | awk '{printf("%d",$1/$2)}'`
    if [ "`getOption -k`" = True ]; then
      kpc=`getValue -k`
      if [ "$kpc" == "" ]; then error "no value to -k option provided"; fi
      c=`checkInteger $kpc`
      if [ $c != ok ]; then error "value give to -k is not an integer"; fi
    fi
  else
    kpc=$kp
  fi

  # if nbands=AUTO leave it out to let vasp determine nbands automatically
  if [ $nbands == "AUTO" ]; then nbandsStr=""; else nbandsStr="NBANDS = $nbands"; fi

  # check if folder exist and should be overwritten
  for a in $aLats; do
    f=$cutoff\eV_$kp\x$kp\x$kpc\kp_$s\eV$addString/$a\Ang
    if [ -e $f -a $overwrite != True ]; then
      error "folder $f exists; run with -f option to force overwriting";
    fi
  done

  # check if templates ok
  cc=`grep -e xxxCUTOFFxxx -e xxxSMEARxxx -e xxxSIGMAxxx -e xxxNBANDSxxx INCAR | wc -l`
  if [ $cc != 4 ]; then error "INCAR template wrong"; fi
  cc=`grep -e xxxKPxxx -e xxxKPCxxx KPOINTS | xargs -n1 | grep -e xxxKPxxx -e xxxKPCxxx | wc -l`
  if [ $cc != 3 ]; then error "KPOINTS template wrong"; fi
  cc=`grep -e xxxALATxxx POSCAR | wc -w`
  if [ $cc != 1 ]; then error "POSCAR template wrong"; fi
  cc=`grep -e xxxCBYAxxx POSCAR | wc -l`
  if [[ ( $type == "hcp" || $type == "dhcp" ) && $cc != 1 ]]; then error "type hcp or dhcp but no xxxCBYAxxx in POSCAR"; fi

  # check EDIFF
  ediffForce=`getOption -ff`
  cc=`sed -nr 's/^[ ]*EDIFF[ ]*=([^;!]*).*/\1/p' INCAR | awk '{if ($1>1e-5) print "larger"}'`
  if [ "$cc" == larger -a $ediffForce != True ]; then error "EDIFF>1e-5; run with -ff option to override"; fi

  # create the input folders
  echo; rm -f jobList
  for a in $aLats; do
    f=$cutoff\eV_$kp\x$kp\x$kpc\kp_$s\eV$addString/$a\Ang
    echo $f; rm -fr $f; mkdir -p $f
    sed -e 's/xxxCUTOFFxxx/'$cutoff'/' \
        -e 's/xxxSMEARxxx/'$smear'/' \
        -e 's/xxxSIGMAxxx/'$sigma'/' \
        -e 's/NBANDS.*xxxNBANDSxxx/'"$nbandsStr"'/' INCAR > $f/INCAR
    if [ -n "$addINCAR" ]; then
      echo "$addINCAR" | xargs -n1 -d \; >> $f/INCAR;
    fi
    sed -e 's/xxxKPxxx/'$kp'/g' \
        -e 's/xxxKPCxxx/'$kpc'/' KPOINTS > $f/KPOINTS
    sed -e 's/xxxALATxxx/'$a'/' \
        -e 's/xxxCBYAxxx/'$cBya'/' POSCAR > $f/POSCAR
    cp POTCAR $f/POTCAR
    echo `pwd`/$f >> jobList
  done

done # end loop over k-points
done # end loop over cutoffs

echo; echo " jobList file written"

