#!/bin/bash

#-----set parameters and paths------------------------------------
tDef=1000
eDef="EIGENVALUES"
wDef="weights"
damDef=1
nepsDef="1e-6"
#-----------------------------------------------------------------


path=`set | grep BASH_SOURCE | sed 's|.*\"\(.*\)/[^/]*\".*|\1|' | sed '/BASH_SOURCE=/s/.*/\./'`
script=`set | grep BASH_SOURCE | sed 's|.*\".*/\([^/]*\)\".*|\1|' | sed '/BASH_SOURCE=/s/.*/\./'`
options=$*; . $path/utilities/functions.include; checkOptions "-h -help -t -e -w -d -n -f"

# small help for options
h=`getOption -h`
if [ $h == True ]; then
  usage $script
  printOptions "-t TEMPERATURE  temperature to use for Fel (default: $tDef K)"             \
               "-t T1 T2 Td     calculate on mesh instead (from T1 to T2 in Td steps)"     \
               "-e EIGENVALUES  eigenvalues file(s) to use for Fel (default: $eDef)"       \
               "-w weights      weights file (default: $wDef)"                             \
               "-d damping      damping factor for search of Fermi level shift (default: $damDef)" \
               "-n neps         eps in nElec for stopping the Fermi level search (default: $nepsDef)" \
               "-f              overwrite previous output"
  exit
fi

# detailed help
help=`getOption -help`
if [ $help == True ]; then
  details $script
  echo2 "   script calculates the electronic free energy from a set of eigenvalues"
  echo2 "   the eigenvalues are expected in an EIGENVALUES file (change with -e)"     \
        "   one eigenvalue per row and in eV; use extractEIGENVALUES.sh for creation"
  echo2 "   two EIGENVALUES files can be given to the -e option corresponding to"     \
        "   spin up and spin down (e.g. EIGENVALUES_up and EIGENVALUES_down)"
  echo2 "   for the '-t T1 T2 Td' option results are written to files (not on screen)"
  echo2 "   with '-d damping' the damping of the Fermi level search magnitude can"    \
        "   be changed; this can be helpful for systems with stronger varying"        \
        "   DOS that do not converge with the default setting (try e.g. -d 0.1)"
  echo2 "   with '-n neps' the eps in the number of electrons for stopping the Fermi" \
        "   level search can be changed; this can be helpful for large systems with"  \
        "   many electrons where the default value is difficult to achieve"           \
        "   (e.g. for the HfScTiZr HEA with 2000 electrons '-n 1e-5' worked)"
  echo2 "   WARNING:" \
        "   to obtain free energies per atom, one needs to explicitly scale by the"   \
        "   number of atoms used in the corresponding eigenvalue calculation"
  exit
fi
 
# get input files names
eOp=`getOption -e`
wOp=`getOption -w`
if [ $eOp == True ]; then evFiles=`getValue -e`; else evFiles=$eDef; fi
if [ $wOp == True ]; then weFile=`getValue -w`; else weFile=$wDef; fi
if [ "$evFiles" == "" ]; then error "no value to -e option given"; fi
if [ "$weFile" == "" ]; then error "no value to -w option given"; fi
if [ "`echo $evFiles | awk '{if (NF!=1&&NF!=2) print "error"}'`" == error ]; then error "only up to two files (spin up and down) can be given to -e option"; fi
if [ "`echo $weFile  | awk '{if (NF!=1) print "error"}'`" == error ]; then error "only a single weights file  can be given to -w option"; fi
check $evFiles $weFile

# get number of electrons
check number_of_electrons_NELECT
nelec=`cat number_of_electrons_NELECT`

# get damping factor
dOp=`getOption -d`
if [ $dOp == True ]; then
  damping=`getValue -d`
  if [ "$damping" == "" ]; then error "-d option requires value"; fi
else
  damping=$damDef
fi

# get neps factor
nOp=`getOption -n`
if [ $nOp == True ]; then
  neps=`getValue -n`
  if [ "$neps" == "" ]; then error "-n option requires value"; fi
else
  neps=$nepsDef
fi

# check if single temperature or mesh given
tOp=`getOption -t`
if [ $tOp == True ]; then
  nt=`getValue -t | awk '{print NF}'`
  if [ $nt != 1 -a $nt != 3 ]; then error "1 or 3 values must be given to -t option"; fi
  T1=`getValue -t | awk '{print $1}'`
  c1=`checkReal $T1`
  if [ "$c1" != ok ]; then error "given temperature wrong; must be real value"; fi

  if [ $nt == 3 ]; then
    T2=`getValue -t | awk '{print $2}'`
    Td=`getValue -t | awk '{print $3}'`
    c2=`checkReal $T2`
    c3=`checkReal $Td`
    if [ "$c2" != ok -o "$c3" != ok ]; then error "given temperature mesh wrong; must be real values"; fi

    # check whether previous output exist
    if [ -e Electronic_entropy_kB -o -e Electronic_free_energy_meV -o -e Electronic_internal_energy_meV -o -e Fermi_level_shift_meV ]; then
      if [ `getOption -f` != True ]; then error "previous output exists; use -f to overwrite"; fi
    fi
  else
    T2=$T1
    Td=1
  fi
else
  nt=1
  T1=$tDef
  T2=$tDef
  Td=1
fi

spin=0
rm -fr _tmp_EIGENVALUES _tmp_weights
for i in $evFiles; do
  cat $i >> _tmp_EIGENVALUES
  cat $weFile >> _tmp_weights
  spin=`expr $spin + 1`
done
$path/fortran/getFel.x $T1 $T2 $Td $nelec $spin $damping $neps
rm _tmp_EIGENVALUES _tmp_weights


