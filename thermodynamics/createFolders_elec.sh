#!/bin/bash

# following 3 lines must always be present
path=`set | grep BASH_SOURCE | sed 's|.*\"\(.*\)/[^/]*\".*|\1|' | sed '/BASH_SOURCE=/s/.*/\./'`
script=`set | grep BASH_SOURCE | sed 's|.*\".*/\([^/]*\)\".*|\1|' | sed '/BASH_SOURCE=/s/.*/\./'`
options=$*; . $path/utilities/functions.include; checkOptions "-h -help -p -I -K -P -a -ff -f"

# small help for options
# space between option and value must be 1
# space between option and Comment or value and Comment must be >2
h=`getOption -h`
if [ $h == True ]; then
  usage $script
  printOptions "-p        create parameters.dat and exit" \
               "-I        create INCAR template and exit" \
               "-K        create KPOINTS template and exit" \
               "-P type   create POSCAR template of type (fcc/bcc/hcp/dhcp) and exit" \
               "-a type   create create all of the above and exit" \
               "-ff       force running even if EDIFF>1e-5" \
               "-f        force template or folder creation even if existing"
  exit
fi

# detailed help
help=`getOption -help`
if [ $help == True ]; then
  details $script
  echo2 "   creates VASP input folders for electronic free energy runs"
  echo2 "   parameters.dat file must be present (generate with -p)"
  echo2 "   INCAR KPOINTS POSCAR templates must be present" \
        "   (generate with -I -K -P bcc or fcc or hcp or dhcp)"
  echo2 "   typical script sequence:"                                   \
        "     \033[1mcreateFolders_Fel.sh\033[0m"                       \
        "     collectOUTCARs.sh          (one folder higher from here)" \
        "     extractFel.sh"                                            \
        "     getFelFit.sh"                                             \
        "     getThermodynamics.sh       (in separate folder)"
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
aLats=a1 a2 a3 a4 ...   # Angstrom (typically from Veq at T=0K up to Veq at T^melt)
sigmas=s1 s2 s3 s4 ...  # smearing in eV (typically 0.01 up to T^melt in eV)

cutoff=                 # eV
kp=N                    # gives NxNxN for fcc/bcc and NxNxInt(N/cBya)
nbands=
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
 ISMEAR =    -1
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

# read in all parameters from parameters.dat
aLats=`get aLats`; cutoff=`get cutoff`; kp=`get kp`; addString=`get addString`
sigmas=`get sigmas`; nbands=`get nbands`; addINCAR=`get addINCAR`; cBya=`get cBya`
checkInput "$aLats" "$cutoff" "$kp" "$sigmas" "$nbands" "$cBya"

# check if folder exist and should be overwritten
for a in $aLats; do for s in $sigmas; do
  f=$cutoff\eV_$kp\x$kp\x$kpc\kp$addString/$a\_$s
  if [ -e $f -a $overwrite != True ]; then
    error "folder $f exists; run with -f option to force overwriting";
  fi
done; done

# check if structure type known
type=`awk 'NR==1{print $2}' POSCAR`
echo type $type
if [ "$type" != "bcc" -a "$type" != "fcc" -a "$type" != "hcp" -a "$type" != "dhcp" ]; then
  error "structure type in POSCAR not known";
fi

# if dhcp, scale c axis (which now equals cBya) times 2
if [ "$type" == "dhcp" ]; then
  cBya=`echo $cBya | awk '{printf "%.10f",$cBya*2}'`
fi

# scale KPC if hcp or dhcp
if [ $type == hcp -o $type == dhcp ]; then
  kpc=`echo $kp $cBya | awk '{printf("%d",$1/$2)}'`
else
  kpc=$kp
fi

# check if templates ok
cc=`grep -e xxxCUTOFFxxx -e xxxSIGMAxxx -e xxxNBANDSxxx INCAR | wc -l`
if [ $cc != 3 ]; then error "INCAR template wrong"; fi
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
for a in $aLats; do for s in $sigmas; do
  f=$cutoff\eV_$kp\x$kp\x$kpc\kp$addString/$a\_$s
  echo $f; 
  [ -d $f ] && echo $f exists! remove it first by hand && exit
  rm -fr $f; mkdir -p $f
  sed -e 's/xxxCUTOFFxxx/'$cutoff'/' \
      -e 's/xxxSIGMAxxx/'$s'/' \
      -e 's/xxxNBANDSxxx/'$nbands'/' INCAR > $f/INCAR
  if [ -n "$addINCAR" ]; then
    echo "$addINCAR" | xargs -n1 -d \; >> $f/INCAR;
  fi
  sed -e 's/xxxKPxxx/'$kp'/g' \
      -e 's/xxxKPCxxx/'$kpc'/' KPOINTS > $f/KPOINTS
  sed -e 's/xxxALATxxx/'$a'/' \
      -e 's/xxxCBYAxxx/'$cBya'/' POSCAR > $f/POSCAR
  cp POTCAR $f/POTCAR
  echo `pwd`/$f >> jobList
done; done

echo; echo " jobList file written"

