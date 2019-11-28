#!/bin/sh

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
  printOptions "-p                 create parameters.dat and exit" \
               "-I                 create INCAR template and exit" \
               "-K                 create KPOINTS template and exit" \
               "-P scType scSize   create POSCAR template of scType (fcc/bcc) and scSize (1-5) and exit" \
               "-a scType scSize   create create all of the above and exit" \
               "-ff                force running even if EDIFF>1e-8" \
               "-f                 force template or folder creation even if existing"
  exit
fi

# detailed help
help=`getOption -help`
if [ $help == True ]; then
  details $script
  echo2 "   creates input folders for a thermodynamic integration calculation"
  echo2 "   parameters.dat file must be present" \
        "   (generate with -p)"
  echo2 "   INCAR KPOINTS templates must be present" \
        "   (generate with -I -K -P bcc or fcc)" \
        "   ../Hessematrix_?x?x?sc folder must be presnet, is not generated automatically at the moment." \
        "   ../Hessematrix_?x?x?sc folder hast to include: POSCAR_{latticeconstant}, HesseMatrix_{latticeconstant}, EqCoords_direct_{latticeconstant}" \
  echo2 "   typical script sequence:"       \
        "     \033[1mcreateFolders_ti.sh\033[0m -a ..." \
        "     \033[1mcreateFolders_dynMat.sh\033[0m"        \
        "     collectOUTCARs.sh"               \
        "     extractForces.sh"                \
        "     getSingleSpeciesPhonons.sh -A                  (cmpc01)"  \
        "     getFqhFit.sh                                   (cmpc01)"  \
        "     getThermodynamics.sh               (in separate folder)"
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
aLats= 4.               # lattice constants in Angstrom; suggested to used volume at the melting point for convergence checks
temps= 3000             # temperatures in K; suggested to use melting temperature for convergence checks
lambdas=1.0             # recommended:0.0 0.15 0.5 0.85 1.0; suggested to use 1.0 for convergence checks
nSeeds=4                # amount of runs for every lambda value; 4 are recommended; 3 are also OK. Depending on ionic steps you do per run
cutoff= 220             # eV

                        # supercell sensitive values, typical values:
                        # 2sc      3sc      4sc       (note: these are only guide values)
ngxf= 120               # 120      180      240       (might need higher settings for very high cutoffs >400eV??)
kp= 4                   # 12       8        6         for a well converged bcc and very well converged fcc calculation
nbands= 100             # 150/180  480/580  1100/xxx  MP/Fermi with 0.1eV for fcc Ca with 8 valence electrons

                        # supercell independent values:
ediff=-2                # electronic convergence in 10^{ediff} eV
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
  echo "
 NPAR = 1

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
 NELM   =    100

 LWAVE  =    F   
 LCHARG =    F   
" > INCAR
  echo; echo "INCAR written"
fi

# cif applies reate KPOINTS
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

# generate POSCAR
genPOS=`getOption -P`
if [ $genPOS == True -o $all == True ]; then
  if [ -e POSCAR -a $overwrite != True ]; then error "POSCAR existing; use -f to overwrite"; fi
  if [ $genPOS == True ]; then type=`getValue -P`; fi
  if [ $all == True ]; then type=`getValue -a`; fi
  if [ -z "$type" ]; then error "supercell type and size type missing"; fi
  sc=`echo $type | sed -e 's/bcc//' -e 's/fcc//' -e 's/ *//g'`
  if [ "$sc" == "" ]; then error "supercell size missing"; fi
  if [ "$sc" -lt 1 -o "$sc" -gt 5 ]; then error "supercell size not supported"; fi
  type=`echo $type | sed -e 's/'"$sc"'//' -e 's/ *//g'`
  if [ $type != bcc -a $type != fcc ]; then error "supercell type not known"; fi
  natoms=`wc -l $path/utilities/$type/coordinates_$sc\x$sc\x$sc\sc | awk '{print $1}'`
  echo "DynMat $type ${sc}x${sc}x${sc}sc
xxxALATxxx
$sc 0 0
0 $sc 0
0 0 $sc
$natoms
Cartesian
xxxDISPxxx   0.0   0.0" > POSCAR
  awk 'NR>1{print}' $path/utilities/$type/coordinates_${sc}x${sc}x${sc}sc >> POSCAR
  echo; echo "POSCAR written"
fi

# exit if we created some template
if [ $genPar == True -o $genINC == True -o $genKPO == True -o $genPOS == True -o $all == True ]; then exit; fi

# check if all input files available
input="parameters.dat INCAR KPOINTS POTCAR POSCAR"
for i in $input; do check $i; done

# read in all parameters from parameters.dat
aLats=`get aLats`; cutoff=`get cutoff`; kp=`get kp`; kpshift=`get kpshift`; addString=`get addString`
smear=`get smear`; sigma=`get sigma`; nbands=`get nbands`; addINCAR=`get addINCAR`; ngxf=`get ngxf`; disp=`get disp`
checkInput "$aLats" "$cutoff" "$kp" "$kpshift" "$smear" "$sigma" "$nbands" "$ngxf" "$disp"

if [ $smear == -1 ]; then s=$sigma; else s=0.0; fi

# check if structure type known
type=`awk 'NR==1{print $2}' POSCAR`
sc=`awk 'NR==1{print $3}' POSCAR | sed 's/^\(.*\)x.*x.*sc/\1/'`
echo; echo " supercell $type ${sc}x${sc}x${sc}sc"
if [ "$type" != "bcc" -a "$type" != "fcc" ]; then
  error "supercell type in POSCAR not known (generate with -P option; see -help)";
fi
if [ "$sc" == "" ]; then error "supercell size in POSCAR empty"; fi 
c=`checkInteger $sc`
if [ $c != ok ]; then error "supercell size non integer"; fi 

# check if templates ok
cc=`grep -e xxxCUTOFFxxx -e xxxNGXFxxx -e xxxSMEARxxx -e xxxSIGMAxxx -e xxxNBANDSxxx INCAR | wc -l | sed 's|[ ]*||g'`
if [ $cc != 7 ]; then error "INCAR template wrong"; fi
cc=`grep -e xxxKPxxx -e xxxKPSHIFTxxx KPOINTS | xargs -n1 | grep -e xxxKPxxx -e xxxKPSHIFTxxx | wc -l | sed 's|[ ]*||g'`
if [ $cc != 6 ]; then error "KPOINTS template wrong"; fi
cc=`grep -e xxxALATxxx -e xxxDISPxxx POSCAR | wc -l | sed 's|[ ]*||g'`
if [ $cc != 2 ]; then error "POSCAR template wrong"; fi

# check EDIFF
ediffForce=`getOption -ff`
cc=`sed -nr 's/^[ ]*EDIFF[ ]*=([^;!]*).*/\1/p' INCAR | awk '{if ($1>1e-8) print "larger"}'`
if [ "$cc" == larger -a $ediffForce != True ]; then error "EDIFF>1e-8; run with -ff option to override"; fi

# bohrradius to angstrom
toAng=0.529177208278835

# check if folder exist and should be overwritten
for a in $aLats; do
  f=$sc\x$sc\x$sc\sc_$cutoff\eV-NGXF$ngxf\-ADDGRID_$kp\x$kp\x$kp\kp-$kpshift\x$kpshift\x$kpshift\shift_$s\eV$addString/$a\Ang
  if [ -e $f -a $overwrite != True ]; then
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

  sed -e 's/xxxCUTOFFxxx/'$cutoff'/' \
      -e 's/xxxNGXFxxx/'$ngxf'/'     \
      -e 's/xxxSMEARxxx/'$smear'/'   \
      -e 's/xxxSIGMAxxx/'$sigma'/'   \
      -e 's/xxxNBANDSxxx/'$nbands'/' INCAR > $f/INCAR
  if [ -n "$addINCAR" ]; then
    echo "$addINCAR" | xargs -n1 -d \; >> $f/INCAR;
  fi
  sed -e 's/xxxKPxxx/'$kp'/g' \
      -e 's/xxxKPSHIFTxxx/'$kpshift'/g' KPOINTS > $f/KPOINTS
  sed -e 's/xxxALATxxx/'$a'/' \
      -e 's/xxxDISPxxx/'$d'/' POSCAR > $f/POSCAR
  cp POTCAR $f/POTCAR
  echo `pwd`/$f >> jobList
done

echo; echo " jobList file written"

