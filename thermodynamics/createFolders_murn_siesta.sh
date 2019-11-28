#!/bin/bash

# following 3 lines must always be present
path=`set | grep BASH_SOURCE | sed 's|.*\"\(.*\)/[^/]*\".*|\1|' | sed '/BASH_SOURCE=/s/.*/\./'`
script=`set | grep BASH_SOURCE | sed 's|.*\".*/\([^/]*\)\".*|\1|' | sed '/BASH_SOURCE=/s/.*/\./'`
options=$*; . $path/utilities/functions.include; checkOptions "-h -help -p -i -a -f"

# small help for options
# space between option and value must be 1
# space between option and Comment or value and Comment must be >2
h=`getOption -h`
if [ $h == True ]; then
  usage $script
  printOptions "-p    create parameters.dat and exit" \
               "-i    create input.fdf template and exit" \
               "-a    create both of the above and exit" \
               "-f    force template or folder creation even if existing"
  exit
fi

# detailed help
help=`getOption -help`
if [ $help == True ]; then
  details $script
  echo2 "   creates input folders for a cold curve (T=0K) run with SIESTA"
  echo2 "   parameters.dat file must be present (generate with -p)"
  echo2 "   input.fdf template must be present (generate with -i bcc or fcc)"
  echo2 "   typical script sequence:"                     \
        "     \033[1mcreateFolders_murn.sh\033[0m -a ..." \
        "     getALatRangeMurn.sh -e ..."                 \
        "     \033[1mcreateFolders_murn.sh\033[0m"        \
        "     extractEnergies_siesta.sh"                  \
        "     getE0KFit.sh           (cmpc01)"            \
        "     getThermodynamics.sh   (needs Fqh additionally)"
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
aLats=3.22 3.2O 3.26    # Angstrom                example bcc Ti 4 valence electrons

elem=Ti
type=bcc                # bcc or fcc
spin=NM                 # NM or FM
xc=PBE                  # LDA or PBE

basis=DZP               # basis size: SZ, DZ, SZP, DZP
kp=8                    # gives NxNxN
smear=1                 # 1=Methfessel, -1=Fermi
sigma=0.1               # smearing in eV

addString=              # if addInput is used it is useful to add some string to name; leave blank otherwise
addInput=               # add additional flags to input.fdf, e.g., addInput=WriteKpoints; leave blank if not needed
                        # more flags can be added by separating them with ;
" > parameters.dat
  echo; echo "parameters.dat written  <==  mind: parameters.dat needs to be adjusted"
fi

# if applies create input.fdf
genINC=`getOption -i`
if [ $genINC == True -o $all == True ]; then
  if [ -e input.fdf -a $overwrite != True ]; then error "input.fdf existing; use -f to overwrite"; fi
  echo "
NumberOfSpecies        1       

%block ChemicalSpeciesLabel
  1  xxxZxxx  xxxELEMxxx
%endblock ChemicalSpeciesLabel

PAO.BasisSize         xxxBASISxxx

LatticeConstant       xxxALATxxx Ang

%block LatticeVectors
 1 0 0
 0 1 0
 0 0 1
%endblock LatticeVectors

AtomicCoordinatesFormat     Fractional
%include fractional_coords.fdf

%block kgrid_Monkhorst_Pack
   xxxKPxxx   0    0    0.5
   0   xxxKPxxx    0    0.5
   0    0   xxxKPxxx    0.5
%endblock kgrid_Monkhorst_Pack

xc.functional         xxxXC1xxx
xc.authors            xxxXC2xxx

MaxSCFIterations              100
DM.Require.Energy.Convergence true

DM.UseSaveDM           true    
DM.NumberPulay         3

Diag.DivideAndConquer .false.
SolutionMethod        diagon  
OccupationFunction     xxxSMEARxxx
ElectronicTemperature  xxxSIGMAxxx eV
" > input.fdf
  echo; echo "input.fdf written"
fi

# exit if we created some template
if [ $genPar == True -o $genINC == True -o $all == True ]; then exit; fi

# check if all input files available
input="parameters.dat input.fdf"
for i in $input; do check $i; done

# read in all parameters from parameters.dat
aLats=`get aLats`; basis=`get basis`; kp=`get kp`; addString=`get addString`
smear=`get smear`; sigma=`get sigma`; addInput=`get addInput`
xc=`get xc`; type=`get type`; spin=`get spin`; elem=`get elem`
checkInput "$aLats" "$basis" "$kp" "$smear" "$sigma" "$xc" "$type" "$spin" "$elem"

if [ $smear == -1 ]; then s=$sigma; else s=0.0; fi

# check if folder exist and should be overwritten
for a in $aLats; do
  f=${basis}_$kp\x$kp\x$kp\kp_$s\eV$addString/$a\Ang
  if [ -e $f -a $overwrite != True ]; then
    error "folder $f exists; run with -f option to force overwriting";
  fi
done

# check if templates ok
cc=`grep -e xxxZxxx -e xxxBASISxxx -e xxxALATxxx -e xxxKPxxx -e xxxXC1xxx -e xxxXC2xxx -e xxxSMEARxxx -e xxxSIGMAxxx input.fdf | wc -l`
if [ $cc != 10 ]; then error "input.fdf template wrong"; fi

# get atomic number from database
file=$path/utilities/atomic_weights_and_ionic_compositions_NIST
check $file
Z=`awk 'BEGIN{Z=0};
        /Atomic Number =/{Ztmp=$NF};
        /Atomic Symbol = '"$elem"'$/{if (Z!=0&&Z!=Ztmp) Z=-1; else Z=Ztmp};
        END{print Z}' $file`
if [ $Z == 0 ];  then error "element $elem not existing in database"; fi
if [ $Z == -1 ]; then error "isotopes have different Atomic Number"; fi

# check if pseudopotential file present
check $elem.psf

# check if given strucuture type ok (fcc or bcc so far)
case "$type" in
  bcc) 
echo "
NumberOfAtoms 2
%block AtomicCoordinatesAndAtomicSpecies
.0 .0 .0 1
.5 .5 .5 1
%endblock AtomicCoordinatesAndAtomicSpecies" > fractional_coords.fdf;;
  fcc)
echo "
NumberOfAtoms 4
%block AtomicCoordinatesAndAtomicSpecies
.0 .0 .0 1
.0 .5 .5 1
.5 .0 .5 1
.5 .5 .0 1
%endblock AtomicCoordinatesAndAtomicSpecies" > fractional_coords.fdf;;
  *) error "structure type $type not supported"
esac

# check if basis ok
if [ $basis != SZ -a $basis != SZP -a $basis != DZ -a $basis != DZP ]; then
 error "basis $basis not supported"
fi

# check if XC ok
case "$xc" in
  LDA) xc1=LDA; xc2=CA;;
  PBE) xc1=GGA; xc2=PBE;;
  *) error "xc type $xc not supported"
esac

# check if smear ok
case "$smear" in
  1) smear=MP;;
  2) smear=FD;;
  *) error "smear type $smear not supported"
esac

# create the input folders
echo; rm -f jobList
for a in $aLats; do
  f=${basis}_$kp\x$kp\x$kp\kp_$s\eV$addString/$a\Ang
  echo $f; rm -fr $f; mkdir -p $f
  sed -e 's/xxxZxxx/'$Z'/'         \
      -e 's/xxxELEMxxx/'$elem'/'   \
      -e 's/xxxBASISxxx/'$basis'/' \
      -e 's/xxxALATxxx/'$a'/'      \
      -e 's/xxxKPxxx/'$kp'/g'      \
      -e 's/xxxXC1xxx/'$xc1'/'     \
      -e 's/xxxXC2xxx/'$xc2'/'     \
      -e 's/xxxSIGMAxxx/'$sigma'/' \
      -e 's/xxxSMEARxxx/'$smear'/' input.fdf > $f/input.fdf
  if [ -n "$addInput" ]; then
    echo "$addInput" | xargs -n1 -d \; >> $f/input.fdf;
  fi
  cp fractional_coords.fdf $elem.psf $f/
  echo `pwd`/$f >> jobList
done

echo; echo " jobList file written"

