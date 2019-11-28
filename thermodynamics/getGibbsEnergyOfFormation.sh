#!/bin/bash

#-----set parameters and paths---------------------------------
out=no #yes #(print additional info for debugging when yes)
nPar=4 # number of needed coefficients in each fit
fits="EVinet EMurn EBirch ECubic"
#--------------------------------------------------------------


path=`set | grep BASH_SOURCE | sed 's|.*\"\(.*\)/[^/]*\".*|\1|' | sed '/BASH_SOURCE=/s/.*/\./'`;[ "$out" = "yes" ] && echo path: $path
script=`set | grep BASH_SOURCE | sed 's|.*\".*/\([^/]*\)\".*|\1|' | sed '/BASH_SOURCE=/s/.*/\./'`;[ "$out" = "yes" ] && echo script: $script
options=$*; . $path/utilities/functions.include; checkOptions "-h -help -v";[ "$out" = "yes" ] && echo options: $options

lom=`getlinuxormac`;[ "$lom" = "Linux" ] && add="";[ "$lom" = "mac" ] && add="'' -e"

h=`getOption -h`
if [ $h = True ]; then
  usage $script
  printOptions
  exit
fi

# detailed help
help=`getOption -help`
if [ $help = True ]; then
  details $script
  echo2 "   the script calculates the Gibbs energy of defect formation from a bulk and defect free energy surface"
  echo2 "   the output consists of a Gibbs energy surface as a function of pressure and temperature that can be"                                         \
        "   directly used for determination of thermodynamic properties employing the getThermodynamics.sh script"                                       \
        "   additionally Gibbs energy files at given pressure, volume, and temperature are written for plotting"                                         \
        "   and checking purposes; the default values for these files are P=0,V=Veq,T=0K and T=Tmelt"
  echo2 "   an optional \"param\" file can be used to set different P,V,T values"                                                                        \
        "   in \"param\" a line starts with the letter P or V or T standing for pressure, volume, temperature "                                          \
        "   the letter is followed by one/more values separated by blanks at which FForm should be calculated (GPa,Ang^3,K)"                             \
        "   the volume refers to the number of atoms in the T=0K energy fit"                                                                             \
        "   to calculate at V=Veq put -1"                                                                                                                
  echo2 "   an optional \"volume_range\" or \"volume_range_atom\" file can exist with \"Vmin Vmax nV\" in the first line"                                \
        "   Vmin and Vmax in Ang^3 and refering to the # of atoms in T=0K energy fit (volume_range) or per atom (volume_range_atom)"                     \
        "   they give the range and mesh (nV) at which FForm at const T is calculated"                                                                   \
        "   and also where it is fitted (for further processing in getThermodynamics.sh)"                                                                \
        "   if the file \"volume_range\" does not exist, the defaults are: Vmin=0.97*Veq, Vmax=1.12*Veq, nV=100"                                         
  echo2 "   if only T=0K fits are supplied an additional \"temperature_range\" file must exists"                                                         \
        "   it needs to contain \"Tmin Tmax Tstep\" (all in K) in the first line"                                                                        \
        "   they give the temperature range at which FForm is printed"
  echo2 "   input files containing the (free) energies have the following naming convention:"                                                            \
        "   defect cell: \033[1mE\033[0mxxx\033[1m_d_\033[0mNNN   \033[1mF\033[0myyy\033[1m_d_\033[0mnn1   \033[1mF\033[0myyy\033[1m_d_\033[0mnn2   ..." \
        "   bulk cell:   \033[1mE\033[0mxxx\033[1m_b_\033[0mMMM   \033[1mF\033[0myyy\033[1m_b_\033[0mmm1   \033[1mF\033[0myyy\033[1m_b_\033[0mmm2   ..." \
        "     where xxx stands for one of the T=0K parametrizations: Vinet, Murn, Birch"                                                                 \
        "           yyy are arbitrary names (e.g., qh, el, ah, mag)"                                                                                     \
        "           NNN is the number of atoms in the T=0K parametrization of the defect cell"                                                           \
        "           MMM is the number of atoms in the T=0K parametrization of the bulk cell"                                                             \
        "           nni is the number of atoms in the i'th free energy parametrization of the defect cell"                                               \
        "           mmi is the number of atoms in the i'th free energy parametrization of the bulk cell"                                                 \
        "     it must hold that NNN+-1=MMM and nni+-1=mmi (+ for vacancy, - for interstitial)"                                                           \
        "     further NNN>=nni and MMM>=mmi"
  echo2 "   Example (3x3x3/2x2x2 vacancy fcc cell):  EVinet_d_107  Fqh_d_31  Fel_d_107  Fah_d_31  Fmag_d_31"                                             \
        "                                            EVinet_b_108  Fqh_b_32  Fel_b_108  Fah_b_32  Fmag_b_32"
  echo2 "   the format of the T=0K fits is:  E0 Veq BM BMder       (EVinet,EMurn,EBirch)"                                                                \
        "                                    E0 Veq BM 3rdOrdCoef  (ECubic)             "
  echo2 "   E0 and V0 in first and second column are important, since this is assumed later in the fortran program"
  exit
fi

echo; echo -n "  checking input ...  "

rm -f _bul_F* _bul_T _bul_input
rm -f _def_F* _def_T _def_input
rm -fr output/
string=""

# check if input files are existing
nd=`ls -1 *_d_* 2> /dev/null | wc -l | sed 's|[ ]*||'`
nb=`ls -1 *_b_* 2> /dev/null | wc -l | sed 's|[ ]*||'`
if [ "$nd" == 0  -o "$nd" != "$nb" ]; then error "number of input files wrong"; fi

if [ $nd == 1 ]; then
  if [ ! -e temperature_range ]; then error "no free energy contributions and no temperature_range file"; fi
  string="$string   temperature_range"
  cat temperature_range | awk '{for (i=$1;i<=$2;i=i+$3) print i}' > _def_T; cp _def_T _bul_T
  cat temperature_range | awk '{for (i=$1;i<=$2;i=i+$3) print 0}' > _def_F__1
  cat temperature_range | awk '{for (i=$1;i<=$2;i=i+$3) print 0}' > _bul_F__1
fi

# check if every defect contribution has a corresponding bulk contribution
def=`ls -1 *_d_*`
for i in $def; do
  bul=`echo $i | sed 's|\(.*\)_d_.*|\1|'`
  nAtv=` echo $i | sed 's|.*_d_\(.*\)|\1|' | awk '{print $1+1}'`
  nAtdv=`echo $i | sed 's|.*_d_\(.*\)|\1|' | awk '{print $1+2}'`
  nAti=` echo $i | sed 's|.*_d_\(.*\)|\1|' | awk '{print $1-1}'`
  if [ ! -e $bul\_b_$nAtv -a ! -e $bul\_b_$nAtdv -a ! -e $bul\_b_$nAti ]; then error "no matching bulk contribution to $i"; fi
done

# we do not support yet magnetization files
M=`ls -1 M*_{b,d}_* 2> /dev/null | wc -l | sed 's|[ ]*||'`
if [ $M != 0 ]; then error "magnetization files (M*_{b,d}_*) are currently not supported"; fi

# check T=0K parametrization
dfits=`echo $fits | xargs -n1 | awk '{print $1"_d_*"}'`
nE=`ls -1 $dfits 2> /dev/null | wc -l | sed 's|[ ]*||'`
if [ "$nE" == 0  ]; then error "no supported Efit existing (supported: $fits)"; fi
if [ "$nE" != 1  ]; then error "too many Efits exisiting (`ls $fits | xargs`)"; fi

# check number of coefficients in defect T=0K parametrization
Edfit=`ls -1 $dfits 2> /dev/null`
EdefAtoms=`echo $Edfit | sed 's|E.*_d_\(.*\)|\1|'`
nP=`head -n1 $Edfit | xargs -n1 | wc -l | sed 's|[ ]*||'`
if [ "$nP" !=  "$nPar" ]; then error "wrong # of coefficients in $Edfit"; fi

# check number of coefficients in bulk T=0K parametrization
bfits=`echo $fits | xargs -n1 | awk '{print $1"_b_*"}'`
Ebfit=`ls -1 $bfits 2> /dev/null`
EbulkAtoms=`echo $Ebfit | sed 's|E.*_b_\(.*\)|\1|'`
nP=`head -n1 $Ebfit | xargs -n1 | wc -l | sed 's|[ ]*||'`
if [ "$nP" !=  "$nPar" ]; then error "wrong # of coefficients in $Ebfit"; fi

# check if temperatures from various free energy contributions match and if temperature step is ok
if [ $nd != 1 ]; then
  Ffit=`ls F*_d_* F*_b_*`
  for f in $Ffit; do wc -l $f | awk '{print $1}' > _$f;  awk '{print $1}' $f >> _$f; done
  ok=`paste _F* | awk 'BEGIN{ok="true";step="ok"};
                           NR==2{t=$1}
                           {for (i=2;i<=NF;i++) if ($i!=$(i-1)) ok="false"}
                           NR>2{if ($1-t>2||$1-t<=0) step="wrong"; t=$1}
                           END{if (ok=="false") print ok; else print step}'`
  if [ "$ok" = false ]; then error "mismatch between the Ffits (number of T points or Tstep)"; fi
  if [ "$ok" = wrong ]; then error "error in Tstep (larger than 2 K, zero or negative)"; fi
  mv _$f _def_T; sed -i $add '1d' _def_T; rm _F*; cp _def_T _bul_T
fi
nT=`wc -l _def_T | awk '{print $1}'`;

# check and prepare parameters file
P="0 0.0001"; V=-1; T="`head -n1 _def_T` `tail -n1 _def_T`"
if [ -e param ]; then
  string="$string   param"
  P=`awk 'BEGIN{p=0};$1=="P"{p=1;for (i=2;i<=NF;i++) printf("%f ",$i)};END{if (p==0) print '"$P"'}' param`
  V=`awk 'BEGIN{p=0};$1=="V"{p=1;for (i=2;i<=NF;i++) printf("%f ",$i)};END{if (p==0) print '"$V"'}' param`
  T=`awk 'BEGIN{p=0};$1=="T"{p=1;for (i=2;i<=NF;i++) printf("%f ",$i)};END{if (p==0) print '"$T"'}' param`
fi
nnP=`echo $P | xargs -n1 | wc -l | sed 's|[ ]*||'`
nnV=`echo $V | xargs -n1 | wc -l | sed 's|[ ]*||'`
nnT=`echo $T | xargs -n1 | wc -l | sed 's|[ ]*||'`
echo $nnP $nnV $nnT > _additional_input;
echo $P >> _additional_input; echo $V >> _additional_input; echo $T >> _additional_input;

# first part of _def_input file
if [ $nd == 1 ]; then nF=1; else nF=`ls -1 F*_d_* 2> /dev/null | wc -l | sed 's|[ ]*||'`; fi
Etype=`echo $Edfit | sed 's|\(.*\)_d_.*|\1|'`
#
# second last 0 is for number of defect Gibbs energy files which only applies to perfect bulk calculation
# last 0 is is for nr of magnetization files (currently not supported)
echo defect 0 $nT $Etype $nPar $nF 0 0 > _def_input
cp _def_input _bul_input # first line the same
cat $Edfit >> _def_input
cat $Ebfit >> _bul_input

if [ $nd != 1 ]; then
  # check for consistency in free energy contributions (number of columns, i.e., coefficients)
  # copy them to _def_F_* files without temperatures for fortran program input
  # add the number of columns, i.e., of coefficients, to _def_input
  Fdfit=`ls -1 F*_d_*`; c=1
  rm -f defAtoms
  for f in $Fdfit; do
    awk 'NR==1{n=NF;error=0};
         {for (i=2;i<=NF;i++) printf("%s ",$i); printf("\n"); if (NF!=n) error=1;};
         END{if (error==1) print "ERROR"}' $f > _def_F__$c;
    if [ "`tail -n1 _def_F__$c`" == "ERROR" ]; then error "inconsistency in columns of $f"; fi
    nV=`head -n1 _def_F__$c | awk '{print NF}'`
    echo $nV >> _def_input
    nAt=`echo $f | sed 's|F.*_d_\(.*\)|\1|'`
    echo $nAt >> defAtoms
    c=`expr $c + 1`;
  done

  # the same as above but for bulk
  rm -f bulkAtoms
  Fbfit=`ls -1 F*_b_*`; c=1
  for f in $Fbfit; do
    awk 'NR==1{n=NF;error=0};
         {for (i=2;i<=NF;i++) printf("%s ",$i); printf("\n"); if (NF!=n) error=1;};
         END{if (error==1) print "ERROR"}' $f > _bul_F__$c;
    if [ "`tail -n1 _bul_F__$c`" == "ERROR" ]; then error "inconsistency in columns of $f"; fi
    nV=`head -n1 _bul_F__$c | awk '{print NF}'`
    echo $nV >> _bul_input
    nAt=`echo $f | sed 's|F.*_b_\(.*\)|\1|'`
    echo $nAt >> bulkAtoms
    c=`expr $c + 1`;
  done
else
  echo 1 >> _bul_input
  echo 1 >> _def_input
fi

if [ -e volume_range -a -e volume_range_atom ]; then error "delete one of volume_range or volume_range_atom"; fi
if [ -e volume_range ]; then
  string="$string   volume_range"
  awk '{for (i=1;i<=NF;i++) printf("%s ",$i); for (i=NF+1;i<=4;i++) printf("-1 ")}; END{printf("\n")}' \
    volume_range >> _additional_input
else
  if [ -e volume_range_atom ]; then
    string="$string   volume_range_atom"
      awk '{for (i=1;i<=NF;i++) if(i<3) printf("%s ",'$EbulkAtoms'*$i); else printf("%s ",$i); for (i=NF+1;i<=4;i++) printf("-1 ")};END{printf("\n")}' \
          volume_range_atom >> _additional_input
  else
    echo -1 -1 -1 -1 >> _additional_input
  fi
fi

echo $EbulkAtoms $EdefAtoms >> _additional_input
if [ $nd != 1 ]; then
  # check if EbulkAtoms > bulkAtoms and EdefAtoms > defAtoms
  ok=`awk 'BEGIN{s="ok"};$1>'$EbulkAtoms'{s="error"};END{print s}' bulkAtoms`
  if [ $ok == error ]; then error "atoms of some free energy bulk paremtrization are larger than for T=0K fit"; fi
  ok=`awk 'BEGIN{s="ok"};$1>'$EdefAtoms'{s="error"};END{print s}' defAtoms`
  if [ $ok == error ]; then error "atoms of some free energy defect paremtrization are larger than for T=0K fit"; fi
  cat bulkAtoms | xargs >> _additional_input
  cat defAtoms  | xargs >> _additional_input
  rm bulkAtoms defAtoms
else
  echo $EbulkAtoms >> _additional_input
  echo $EdefAtoms >> _additional_input
fi

echo "input ok"; echo
if [ -n "$string" ]; then echo -e " \033[1m$string\033[0m   file(s) read in"; echo; fi
echo "  contributions read in:"
echo "    defect: `ls E*_d_* F*_d_* 2> /dev/null | xargs`"
echo "    bulk:   `ls E*_b_* F*_b_* 2> /dev/null | xargs`"
echo
mkdir output

$path/fortran/getGibbsEnergyOfFormation.x

