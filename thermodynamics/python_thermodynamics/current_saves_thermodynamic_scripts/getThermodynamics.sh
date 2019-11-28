#!/bin/bash

#-----set parameters and paths---------------------------------
Pdef=0.0001 # default pressure in GPa
fits="EVinet EMurn EBirch ECubic"
nPar=4 # number of needed coefficients in each fit
strTypeList="bcc fcc hcp none" # structure types supported
eps=0.0500; # pressure eps for Bprime fitting
# the eps values for the other derivatives are given in
# $path/fortran/constants.f90
#--------------------------------------------------------------


# following 3 lines must always be present
path=`set | grep BASH_SOURCE | sed 's|.*\"\(.*\)/[^/]*\".*|\1|' | sed '/BASH_SOURCE=/s/.*/\./'`
script=`set | grep BASH_SOURCE | sed 's|.*\".*/\([^/]*\)\".*|\1|' | sed '/BASH_SOURCE=/s/.*/\./'`
options=$*; . $path/utilities/functions.include; checkOptions "-h -help -helpin -helpout -P -s -f -e -o -v"


# small help for options
# space between option and value must be 1
# space between option and Comment or value and Comment must be >2
h=`getOption -h`
if [ $h = True ]; then
  usage $script
  printOptions "-P pressure   pressure in GPa to use for the thermodynamics (default: $Pdef GPa)"               \
               "-s strType    use strType ($strTypeList) as structure type to calculate lattice expansions" \
               "-f            force execution even if temperature step is too large (>2K)" \
               "-e            make Fqh Fah Fel ... files have equal length (will remove Temperatures which are not present in both)" \
               "-o            override strType from path with strType from -s or file"
               "-v            be more verbose"
  echo "Note:    as an alternative to the '-s strType' option you can either place an empty file with the name" 1>&2
  echo "         'strType' into the work folder, e.g., touch fcc, or put the strType into the path of the folder" 1>&2
  exit
fi

# detailed help
help=`getOption -help`
if [ $help = True ]; then
  details $script
  echo2 "   calculates the temperature dependence of various thermodynamic quantities at constant"    \
        "   pressure from the sum of the various supplied free energy surfaces"
  echo2 "   the input must consist of a single T=0K parametrization and at least one free energy"     \
        "   contribution; several free energy contributions are possible"                             
  echo2 "   the T=0K parametrization can be one of: EVinet, EMurn, EBirch, ECubic; use '-helpin'"     \
        "   for details on these files"                                                               
  echo2 "   the free energy contributions must have file names that start with 'F' but the 'F'"       \
        "   can be followed by an arbitrary string; for convenience it is suggested to name the"      \
        "   free energy files according to the physical contribution they describe, e.g., for a"      \
        "   quasiharmonic, anharmonic, electronic, and magnetic calculation:"                         \
        "      Fqh  Fah  Fel  Fmag"
  echo2 "   each free energy file must contain a parametrized free energy surface as a function"      \
        "   of volume and temperature; use '-helpin' for further details"
  echo2 "   in addition to the free energy files, files containing Gibbs energies of formation of"    \
        "   various defects can be supplied; the file names must start with 'Gform_' and can have"    \
        "   an arbitrary addition, e.g., Gform_vac or Gform_int"
  echo2 "   correctly formatted Gform files (check '-helpin' for details) can be obtained with:"      \
        "   $path/getGibbsEnergyOfFormation.sh"
  echo2 "   the final Gibbs energy surface used to calculate the thermodynamics is constructed as:"   \
        "     G(P,T) = sum_i F_i(V_P,T) + P*V_P(T) - kB*T sum_j exp(-Gform_j/kB*T)"                   \
        "   where i runs over the various free energy files and j over the Gform files; further P"    \
        "   is the given pressure and V_P is the volume adjusted such that -dF/dV (with F=sum_i F_i)" \
        "   corresponds to P"
  echo2 "   the generated output is put into a separate folder with name 'output_xxGPa', where"       \
        "   xx is the given pressure (default: $Pdef GPa; change with '-P pressure' option)"
  echo2 "   details and units of the generated files and quantities can be obtained with '-helpout'"
  echo2 "   by default only the volume expansion is generated; if however the '-s strType' option is" \
        "   supplied, a lattice expansion is generated according to strType (currently supported:"    \
        "   $strTypeList); instead of using the '-s strType' option, the structure type can be"       \
        "   supplied by placing an empty file into the work folder with name strType or by putting"   \
        "   strType into the path of the work folder; if strType==hcp an additional cBya file"        \
        "   containing the c/a ratio is needed (check '-helpin')"
  echo2 "   for calculation of SFEs within the ANNNI model an additional feature is provided:"        \
        "   a volume_expansion file from a previous $script run can be provided as input"             \
        "   in order to generate a free_energy_at_volume_expansion file; this file can be then used"  \
        "   in conjunction with a previously created Gibbs_energy file to determine the SFE by"       \
        "   calculating the difference of the two; typically, the first $script run will"             \
        "   correspond to an fcc structure and the second to an hcp structure"
  echo2 "   WARNING: since the script works with certain file names, avoid keeping files in the"      \
        "   folder other than the input files in order to prevent errors"
  exit
fi

# more details on the format of the input files
helpin=`getOption -helpin`
if [ $helpin = True ]; then
  echo; echo  -e "\033[1m\033[31mINPUT FORMAT DETAILS\033[30m: $script\033[0m"
  echo2 "   \033[1mT=0K file:\033[0m"                                                                \
        "   the T=0K parametrization can be a Vinet, Murnaghan, Birch-Murnaghan, or a third"         \
        "   order fit and the corresponding files must be named EVinet, EMurn, EBirch, or ECubic"    \
        "   the first three files must consist of a single line containing the equilibrium"          \
        "   properties at T=0K (energy, volume, bulk modulus, pressure derivative of bulk modulus):" \
        "     E0(meV)  V0(Ang^3)  B0(GPa)  B0'"                                                      \
        "   whereas the format for ECubic is:"                                                       \
        "     E0(meV)  V0(Ang^3)  B0(GPa)  thirdOrderCoeff(meV/Ang^6)"
  echo2 "   \033[1mF* files:\033[0m"                                                                 \
        "   the free energy surface in each free energy file must be discretely parametrized along"  \
        "   the temperature axis due to the typically different analytical behavior of the various"  \
        "   physical contributions as function of T; in contrast the volume dependence needs to be"  \
        "   parametrized with polynomials since they are a well converging basis along this axis"
  echo2 "   the format of a free energy file is therefore:"                                          \
        "     T_1(K)  a_1(meV)  b_1(meV/Ang^3)  c_1(meV/Ang^6)  ..."                                 \
        "     T_2(K)  a_2(meV)  b_2(meV/Ang^3)  c_2(meV/Ang^6)  ..."                                 \
        "      .       .         .               .               . "                                 \
        "      .       .         .               .               . "                                 \
        "     T_N(K)  a_N(meV)  b_N(meV/Ang^3)  c_N(meV/Ang^6)  ..."                                 \
        "   the index runs here over the temperature mesh and a,b,c,... are the 0., 1., 2., ..."     \
        "   order coefficients of the volume parametrization"
  echo2 "   a single free energy file must have the same number of coefficients in each line;"       \
        "   the number of volume coefficients may vary among different free energy files, however"   \
        "   the number and the values of the temperature mesh points must be exactly the same"
  echo2 "   the delta for the temperature mesh needs to be smaller or equal to 2 K and it is"        \
        "   suggested to start at T_1 = 1 K; T_N is typically the melting temperature"
  echo2 "   \033[1mGform_* files:\033[0m"                                                            \
        "   the Gibbs energies of formation are parametrized in a similar format as the free"        \
        "   energy files; they must also contain the same number of temperature points, but for"     \
        "   each temperature a pressure parametrization is expected:"                                \
        "     T_1(K)  a_1(meV)  b_1(meV/GPa)  c_1(meV/GPa^2)  ..."                                   \
        "     T_2(K)  a_2(meV)  b_2(meV/GPa)  c_2(meV/GPa^2)  ..."                                   \
        "      .       .         .             .               . "                                   \
        "      .       .         .             .               . "                                   \
        "     T_N(K)  a_N(meV)  b_N(meV/GPa)  c_N(meV/GPa^2)  ..."
  echo2 "   \033[1mcBya file:\033[0m"                                                                \
        "   if strType=hcp, a file 'cBya' containing in the first line the T=0K volume dependence"   \
        "   of the aLat/cLat ratio needs to be provided; format:   a  b(1/Ang^3)  c(1/Ang^6)  ..."   \
        "   with the 0., 1., 2., ... order coefficients given by a, b, c, ...; if the volume"        \
        "   dependence is not known just provide the 'a' coefficient, i.e., a single number in cBya"
  exit
fi

# more help on the format of the output files
helpout=`getOption -helpout`
if [ $helpout = True ]; then
  echo; echo  -e "\033[1m\033[31mOUTPUT FORMAT DETAILS\033[30m: $script\033[0m"
  echo2 "   the script generates:"                                                                   \
        "     \033[1mfileName                       unit          comment\033[0m"                    \
        "     volume_expansion               Ang^3         V(T)"                                     \
        "     volume_expansion_relative      %             v(T)=[V(T)-V1]/V1, where V1 corresponds"  \
        "                                                  to the 1. supplied temperature, i.e.,"    \
        "                                                  in a typical case: V1=V(T=1K)"            \
        "     volume_expansion_coefficient   10^-5 K^-1    b(T)=[V(T)-V(T-dT)]/dT/V(T), where dT is" \
        "                                                  the supplied temperature step (1..2 K)"   \
        "     Gibbs_energy                   meV           G(T)"                                     \
        "     enthalpy                       meV           H(T)"                                     \
        "     PV_term                        meV           P*V(T)"                                   \
        "     entropy                        kB            S(T)"                                     \
        "     heat_capacity_isobaric         kB            CP(T)"                                    \
        "     heat_capacity_isochoric        kB            CV(T)"                                    \
        "     compressibility_isothermal     GPa^-1        k(T)"                                     \
        "     bulk_modulus_isothermal        GPa           BT(T)=1/k(T)"                             \
        "     bulk_modulus_adiabatic         GPa           BA(T)=BT(T)*CP(T)/CV(T)"                  \
        "     Bprime_isothermal              unitless      pressure derivative of BT(T)"             \
        "     Bprime_adiabatic               unitless      pressure derivative of BA(T)"             \
        "     TVbetaSqrBT                    kB            =T*V(T)*b(T)^2*BT(T)"                     \
        "     Grueneisen_parameter           unitless      =V(T)*b(T)*BT(T)*CV(T)"
  echo2 "   the generated temperature dependent thermodynamic quantities have the same temperature"  \
        "   mesh as the input free energy files (unit: K)"
  echo2 "   additionally the file free_energy_noET0K_vsV_atTmax is produced which contains the free" \
        "   energy as a function of the volume without the T=0K contribution; the starting volume"   \
        "   is V1, the equilibrium volume at the first temperature, and the end volume is the"       \
        "   equilibrium volume at the largest supplied temperature"
  echo2 "   if a structure type is provided (either by -s option or by strType file or by path)"     \
        "   lattice expansions are calculated in addition to the volume expansion:"                  \
        "     aLat_expansion                 Ang           a(T)=f(V(T)), where f() is a function"    \
        "                                                  that depends on the structure type"       \
        "     aLat_expansion_relative        %             e(T)=[a(T)-a1]/a1"                        \
        "     aLat_expansion_coefficient     10^-5 K^-1    b(T)=[a(T)-a(T-dT)]/dT/a(T)"
  echo2 "   the function f is, e.g., (4.*V(T))^(1./3.) for fcc, (2.*V(T))^(1./3.) for bcc, or"       \
        "   ((V(T)/sin(pi/3.))/cBya(T))^(1./3.) for hcp"
  echo2 "   for hcp cLat_expansion files are additionally written having the same format as the"     \
        "   aLat_expansion files and also a cBya file containing the temperature dependence of the"  \
        "   aLat/cLat ratio cBya(T); the latter is obtained implicitly from the supplied volume"     \
        "   dependence of the cBya ratio at T=0K (cBya file; see '-helpin') as cBya(V(T))"
  echo2 "   if a volume_expansion file is provided (see comment on SFEs in -help) then additionally" \
        "   a free_energy_at_volume_expansion file is generated containing the free_energy at the"   \
        "   provided volume_expansion in meV; note that the P*V term is not included"
  exit
fi

echo; echo -n "  checking input ...  "
rm -f _F* _G* _T _input _M* _volume_expansion

# make Fx files have equal length
makeequal() {
python << END

import glob
import numpy as np
import utils as u
import sys

possible_filenames = [ "Fqh*", "Fah*", "Fel*"]
filenames = []
f = []

for i in possible_filenames:
    filename = glob.glob(i)
    if len(filename) == 1:
        filenames.append(filename[0])

#print "filenames:",filenames
tempslist = []
for i in filenames:
    data = np.loadtxt(i)
    tempslist.append(data[:,0])
    f.append(data)
temps = u.schnittmenge_several_lists(tempslist, verbose = False)
#print temps




# getting equally lenghy fqh's (only intersecting temperatures)"
for ind,i in enumerate(f):
    o = np.empty((len(temps),i.shape[1]))
    #print ind,i
    for tind, t in enumerate(temps):
        # get index of corresponding data
        iind = np.nonzero( t == i[:,0])[0][0]
        #print tind,t,iind
        o[tind] = i[iind]
    np.savetxt(filenames[ind],o)

END
#lom=`getlinuxormac`
#gjoin="----------";
#[ "$lom" = "mac" ] && gjoin=`which gjoin`
#join=`which join`
#[ -e "$gjoin" ] && join=$gjoin
#filesdef=`ls -1d F{el,qh,ah}\_[db]_[0-9]* 2> /dev/null`
#filesbulk=`ls -1d F{el,qh,ah} 2> /dev/null`
#files=`echo "$filesdef" "$filesbulk"`
#echo ... files: $files
#for file1 in $files;do
#    for file2 in $files;do
#        [ "$file1" = "$file2" ] && continue
#        wcf1=`wc -l $file1 | awk '{print $1}'`
#        wcf2=`wc -l $file2 | awk '{print $1}'`
#        cf1=`cat $file1 | awk '{print NF}' | sort | uniq`
#        cf2=`cat $file2 | awk '{print NF}' | sort | uniq`
#        #echo $cf1 | wc -w | sed 's|[ ]*||g'
#        #echo $cf2 | wc -w | sed 's|[ ]*||g'
#        if [ "`echo $cf1 | wc -w | sed 's|[ ]*||g'`" != "1" ] || [ "`echo $cf2 | wc -w | sed 's|[ ]*||g'`" != "1" ];then
#            echo "it seems file1:$file1: contains different rows:cf1:$cf1:"
#            echo "it seems file2:$file2: contains different rows:cf2:$cf2:"
#            echo "file1: $file1    wcf1:$wcf1   ch1:$cf1:"
#            echo "file2: $file2    wcf2:$wcf2   ch2:$cf2:"
#            exit
#            fi
#        [ "`echo $cf2 | wc -w | sed 's|[ ]*||g'`" != "1" ] && echo "it seems $file2 contains different rows:$cf2:" && exit 
#        rm -f $file1\_tmp  $file2\_tmp
#        #nf1=`cat $file1 | awk '{print NF}' | sort | uniq`
#        #nf2=`cat $file2 | awk '{print NF}' | sort | uniq`
#        #[ "`echo $nf1 | wc -w`" != "1" ] && echo lines unequal in $file1
#        #[ "`echo $nf2 | wc -w`" != "1" ] && echo lines unequal in $file2
#        #
#
#        #nf=`awk 'FNR==NR{a[$1]=$2 $3;next}{ print $0, a[$1]}' $file1 $file2 | awk '{print NF}' | sort -n | uniq | tail -1`
#        #if [ "$nf" = "6" ];then
#        #    awk 'FNR==NR{a[$1]=$2 $3;next}{ print $0, a[$1]}' Fqh Fah | awk 'NF==6{print $0}' | awk '{print $1,$2,$3,$4}' > $file1\_tmp;
#        #    awk 'FNR==NR{a[$1]=$2 $3;next}{ print $0, a[$1]}' Fqh Fah | awk 'NF==6{print $0}' | awk '{print $1,$5,$6,$4}' > $file1\_tmp;
#
#        $join --nocheck-order -1 1 -2 1 $file1 $file2 | cut -d" " -f1-$cf1 > $file1\_tmp;
#        $join --nocheck-order -1 1 -2 1 $file2 $file1 | cut -d" " -f1-$cf2 > $file2\_tmp; 
#         
#        mv $file1\_tmp $file1
#        mv $file2\_tmp $file2      
#done
#done
}
[ "`getOption -e`" = "True" ] && makeequal

# check input pressure
P=`getOption -P`; if [ $P = True ]; then P=`getValue -P`; else P=$Pdef; fi
c=`checkReal $P`; if [ "$c" != ok ]; then error "value to -P option empty or wrong"; fi

# error if supplied pressure has more than four digits after the dot because this is the number
# of digits used for the output folder and we want to prevent ambigous output
# such a precision should be sufficient, if not change simply the format (one and three lines below)
c=`echo $P | awk '{printf("%.4f",$1)}' | awk '{if ($1=='$P') print "ok"}'`
if [ "$c" != ok ]; then error "provided pressure should have no more than four digits after the dot"; fi
P=`echo $P | awk '{printf("%.4f",$1)}'`

# this variable keeps track of where the strType comes from
from=""

# check if input structure type given by -s option
strOp=`getOption -s`;
if [ $strOp = True ]; then
  str=`getValue -s | sed 's| ||g'`
  c=`echo $strTypeList | xargs -n1 | awk '$1=="'"$str"'"{print "ok"}'`
  if [ "$c" != ok ]; then error "value given to -s option empty or wrong (supported strTypes: $strTypeList) str:$str: c:$c:"; fi
  from="(from -s option)"; strType=$str
fi



# check if at most one strType file exists
[ "`getOption -v`" = "True" ] && echo "# check if at most one strType file exists"
nStrFi=`ls -1 $strTypeList 2> /dev/null | wc -l`;
strFi=`ls $strTypeList 2> /dev/null | xargs`;
if [ $nStrFi -gt 1 ]; then error "more than 1 structure type file existing: $strFi"; fi

# error if strType exists and -s option given and both are unequal
[ "`getOption -v`" = "True" ] && echo "# error if strType exists and -s option given and both are unequal" 
if [ $nStrFi = 1 -a $strOp = True -a "$str" != "$strFi" ]; then
  error "strType from -s option ($str) in conflict with strType file ($strFi)";
fi

# take strType from file is exists and no -s option
if [ $nStrFi = 1 -a $strOp != True ]; then from="(from strType file)"; strType=$strFi; fi

# check for strType in path
[ "`getOption -v`" = "True" ] && echo "# check for strType in path" 
for i in $strTypeList; do
  strPath=`pwd | sed 's/.*\('"$i"'\).*/\1/'`
  if [ "$strPath" = $i ]; then
    oOp=`getOption -o`
    # if existing, check if strType from path in conflict with the one from -s option or from the file
    if [ $strOp = True -a "$str" != $strPath -a $oOp != True ]; then
      error "strType from -s option ($str) in conflict with strType extracted from path ($strPath); use -o to override path with -s option";
    fi
    if [ $nStrFi = 1 -a "$strFi" != $strPath -a $oOp != True ]; then
      error "strType file ($strFi) in conflict with strType extracted from path ($strPath); use -o to override path with strType from file";
    fi
    if [ $strOp != True -a $nStrFi != 1 ]; then from="(from path)"; strType=$strPath; fi
  fi
done

# if $from is still empty at this place we have no strType provided
[ "`getOption -v`" = "True" ] && echo "# if $from is still empty at this place we have no strType provided" 
if [ -z "$from" ]; then strType=none; fi

# if strType=hcp check if cBya file provided
if [ $strType = hcp -a ! -e cBya ]; then error "strType = hcp but no cBya file provided"; fi

# check T=0K parametrization
[ "`getOption -v`" = "True" ] && echo "# check T=0K parametrization" 
nE=`ls -1 $fits 2> /dev/null | wc -l | sed 's/^[ \t]*//;s/[ \t]*$//'`
[ "`getOption -v`" = "True" ] && echo "# check T=0K parametrization :nE:$nE"
if [ "$nE" = 0  ]; then error "no supported Efit existing (supported: $fits)"; fi
if [ "$nE" != 1  ]; then error "too many Efits exisiting (`ls $fits | xargs`)"; fi

# check number of coefficients in T=0K parametrization
[ "`getOption -v`" = "True" ] && echo "# check check number of coefficients in T=0K parametrization" 
Efit=`ls -1 $fits 2> /dev/null`
[ "`getOption -v`" = "True" ] && echo "# check check number of coefficients in T=0K parametrization :Efit:$Efit"
nP=`head -n1 $Efit | xargs -n1 | wc -l | sed 's/^[ \t]*//;s/[ \t]*$//'`
[ "`getOption -v`" = "True" ] && echo "# check check number of coefficients in T=0K parametrization :nP:$nP"
if [ "$nP" !=  "$nPar" ]; then error "wrong # of coefficients in $Efit"; fi

# check if enough free energy contributions given
nF=`ls -1 F* 2> /dev/null | wc -l | sed 's/^[ \t]*//;s/[ \t]*$//'`
if [ "$nF" = 0 ]; then error "no F* files available (must be other than Gform_* and M*)"; fi

# prepend number of lines to each free energy file and extract only the temperature into same file with _ prepended
Ffit=`ls F* Gform_* M* volume_expansion 2> /dev/null`
for f in $Ffit; do awk 'END{print NR}' $f  > _$f;  awk '{print $1}' $f >> _$f; done

# check if temperatures from various free energy contributions match (stored now in the _F*, _Gform_*, _M*, and _volume_expansion files)
[ "`getOption -v`" = "True" ] && echo "# check if temperatures from various free energy contributions match"
list=`ls _F* _Gform_* _M* _volume_expansion 2> /dev/null`
ok=`paste $list | \
             awk 'BEGIN{ok="true";step="ok"};
                  NR==2{t=$1}
                  {for (i=2;i<=NF;i++) if ($i!=$(i-1)) ok="false"}
                  NR>2{if ($1-t>2||$1-t<=0) step="wrong"; t=$1}
                  END{if (ok=="false") print ok; else print step}'`
if [ "$ok" = false ]; then
  error "mismatch between the F* and/or Gform_* and/or M* and/or volume_expansion files (number of T points or Tstep)"
fi

# if temperature step too large, i.e., $ok==wrong, run only if -f option is given
[ "`getOption -v`" = "True" ] && echo "# check if temperature step too large"
force=`getOption -f`;
if [ "$ok" = wrong -a $force != True ]; then
  error "error in Tstep (larger than 2 K, zero or negative); run with -f to force calculation if Tstep>2K";
fi

# get the temperature file and number of temperature points from the last of the free energy files
# $f stores the last value from the above for loop

# lom see ~/Thermodynamics/utilities/compatibility_issues
lom=`getlinuxormac`;
[ "`getOption -v`" = "True" ] && echo "# lom :lom:$lom"
[ "$lom" = "Linux" ] && add="";[ "$lom" = "mac" ] && add="'' -e"
[ "`getOption -v`" = "True" ] && echo "# add :add:$add"
#echo lom:$lom
#echo add:$add
mv _$f _T; sed -i $add '1d' _T; nT=`wc -l _T | awk '{print $1}'`;
rm -f _F* _Gform_* _M* _volume_expansion

# first part of _input file
[ "`getOption -v`" = "True" ] && echo "# first part of _input file"
nF=`ls -1 F* | wc -l | sed 's/^[ \t]*//;s/[ \t]*$//'`
nG=`ls -1 Gform_* 2> /dev/null | wc -l | sed 's/^[ \t]*//;s/[ \t]*$//'`
nM=`ls -1 M* 2> /dev/null | wc -l | sed 's/^[ \t]*//;s/[ \t]*$//'`
echo $strType $P $nT $Efit $nPar $nF $nG $nM > _input
awk 'NR==1{print $0}' $Efit >> _input

# if hcp prepare cBya input
if [ $strType = hcp ]; then
  nV=`awk 'BEGIN{nV=0};NR==1{nV=NF};END{print nV}' cBya`
  if [ $nV -lt 1 ]; then error "cBya file corrupt"; fi
  echo $nV >> _input
  head -n1 cBya >> _input
fi

# check for consistency in free energy contributions (number of columns, i.e., coefficients)
# copy them to _F_* files without temperatures for fortran program input
# add the number of columns, i.e., of coefficients, to _input
Ffit=`ls -1 F*`; c=1
for f in $Ffit; do
  awk 'NR==1{n=NF;error=0};
       {for (i=2;i<=NF;i++) printf("%s ",$i); printf("\n"); if (NF!=n) error=1;};
       END{if (error==1) print "ERROR"}' $f > _F__$c;
       #This command produces an 'end-line' symbol at the end of each line, which fortran subsequently gets confused by (on some systems)
       #This can be fixed by adding the following two lines after the above command
       #^M has to be entered Ctrl+V Ctrl+M
       tr -d '' < _F__$c > newxx
       mv newxx _F__$c
  if [ "`tail -n1 _F__$c`" = "ERROR" ]; then error "inconsistency in columns of $f"; fi
  nV=`head -n1 _F__$c | awk '{print NF}'`
  echo $nV >> _input
  c=`expr $c + 1`;
done

# the same as above but now for Gibbs energy of formation contributions
# additionally add also their name to _input
Gfit=`ls -1 Gform_* 2> /dev/null`; c=1
for f in $Gfit; do
  awk 'NR==1{n=NF;error=0};
       {for (i=2;i<=NF;i++) printf("%s ",$i); printf("\n"); if (NF!=n) error=1;};
       END{if (error==1) print "ERROR"}' $f > _G__$c;
       #This command produces an 'end-line' symbol at the end of each line, which fortran subsequently gets confused by (on some systems)
       #This can be fixed by adding the following two lines after the above command
       #^M has to be entered Ctrl+V Ctrl+M
       tr -d '' < _G__$c > newxx
       mv newxx _G__$c
  if [ "`tail -n1 _G__$c`" = "ERROR" ]; then error "inconsistency in columns of $f"; fi
  nV=`head -n1 _G__$c | awk '{print NF}'`
  name=`echo $f | sed 's|Gform_\(.*\)|\1|'`
  echo $nV $name >> _input
  c=`expr $c + 1`;
done

# the same now for the M* files containing magnetic moments
Mfit=`ls -1 M* 2> /dev/null`; c=1
for m in $Mfit; do
  awk 'NR==1{n=NF;error=0};
       {for (i=2;i<=NF;i++) printf("%s ",$i); printf("\n"); if (NF!=n) error=1;};
       END{if (error==1) print "ERROR"}' $m > _M__$c;
       #This command produces an 'end-line' symbol at the end of each line, which fortran subsequently gets confused by (on some systems)
       #This can be fixed by adding the following two lines after the above command
       #^M has to be entered Ctrl+V Ctrl+M
       tr -d '' < _M__$c > newxx
       mv newxx _M__$c
  if [ "`tail -n1 _M__$c`" == "ERROR" ]; then error "inconsistency in columns of $m"; fi
  nV=`head -n1 _M__$c | awk '{print NF}'`
  echo $nV >> _input
  c=`expr $c + 1`;
done

# now for the volume_expansion file if it exists (must have two columns)
if [ -e volume_expansion ]; then
  awk 'NR==1{error=0};
       {print $2; if (NF!=2) error=1;};
       END{if (error==1) print "ERROR"}' volume_expansion > _volume_expansion;
  if [ "`tail -n1 _volume_expansion`" == "ERROR" ]; then error "inconsistency in columns of volume_expansion"; fi
fi

echo "  input ok"; echo
# print important parameters
if [ "$lom" = "mac" ];then
        echo -e "  structure type:    \033[31m\033[1m$strType\033[0m $from"
        echo -e "  contributions:     \033[1m`ls E* F* M* Gform_* 2> /dev/null | xargs`\033[0m"
    else
        /bin/echo -e "  structure type:    \033[31m\033[1m$strType\033[0m $from"
        /bin/echo -e "  contributions:     \033[1m`ls E* F* M* Gform_* 2> /dev/null | xargs`\033[0m"
    fi
echo "  pressure:          $P GPa"
T1=`awk 'NR==1{print $1}' F*`
T2=`awk 'END{print $1}' F*`
Tstep=`head -n2 F* | awk 'NR==1{t1=$1}; NR==2{t2=$1}; END{print t2-t1}'`
echo "  temperature range: ${T1} ... ${T2} in ${Tstep} K steps"; echo

# prepare temp output folder needed in the fortran program
rm -fr _output; mkdir _output
rm -f tmpxyz

# save inputfiles
folder=inputfiles_original
[ ! -e "$folder" ] && mkdir $folder
cp -a _* $folder

# check and run fortran program where the actual work is done
check $path/fortran/getThermodynamics.x
$path/fortran/getThermodynamics.x | tee -a tmpxyz

# check if converged and rerun to max T if not
check=`grep -o "running ...  cannot converge within" tmpxyz`
if [ "$check" = "running ...  cannot converge within" ];then
    TT=`grep "running ...  cannot converge within" tmpxyz | sed 's|.*T=||' | awk '{print $1-1}'`
    echo "#########################################"
    echo "# Trying again until Temperature: $TT K #"
    echo "#########################################"
    #rm -fr _F__* _G__* _M__* _input _T _output _volume_expansion inputfiles_original tmpxyz 
    #rm -r _T*
    files=`ls -l F* | grep ^- | awk '{print $NF}' | xargs`
    files2=`cd $folder;ls -l _F* | grep ^- | awk '{print $NF}' | xargs`
    files=`ls -l F* | grep ^- | awk '{print $NF}' | xargs`
    for i in $files;do
        mv $i $folder/$i
        cat $folder/$i | awk '$1<='"$TT"' {print $0}' > $i
        done
    cp $folder/_input .
    lines_orig=`wc -l $folder/_T | awk '{print $1}'`
    cat $folder/_T | awk '$1<='"$TT"' {print $0}' > _T
    lines=`wc -l _T | awk '{print $1}'`
    sed -i $add 's| '"$lines_orig"' | '"$lines"' |' _input
    #echo ::$lines::
    #echo ::$lines_orig::
    for i in $files2;do
        head -n+$lines $folder/$i > $i
        #echo $i
        done
    check $path/fortran/getThermodynamics.x
    $path/fortran/getThermodynamics.x
    for i in $files;do
        cp $folder/$i .
        done
    fi
rm -rf tmpxyz $folder

# check if convergence problems occured in which case no files would have been produced
if [ ! -e _output/volume_expansion ]; then
  rm -fr _F__* _G__* _M__* _input _T _output _volume_expansion
  error "aborting due to convergence problems";
fi

# move ouput files to final ouput folder
outputFolder="output_${P}GPa"; rm -fr $outputFolder;  mv _output $outputFolder

# calculate now everything at a pressure +eps and -eps for Bprime input
peps=`echo $P $eps | awk '{print $1+$2}'`
sed -i $add '1s/^\([^ ]\+\) \+[^ ]\+ \(.*\)/\1 '$peps' \2/' _input
mkdir _output; $path/fortran/getThermodynamics.x > /dev/null; mv _output _output_peps
meps=`echo $P $eps | awk '{print $1-$2}'`
sed -i $add '1s/^\([^ ]\+\) \+[^ ]\+ \(.*\)/\1 '$meps' \2/' _input
mkdir _output; $path/fortran/getThermodynamics.x > /dev/null; mv _output _output_meps

# check if convergence problems occured in which case no files would have been produced
if [ ! -e _output_peps/volume_expansion -o ! -e _output_meps/volume_expansion ]; then
  rm -fr _F__* _G__* _M__* _input _T _output_peps _output_meps _volume_expansion
  error "convergence problems in calculating Bprime; other results might still be useful";
fi

# calculate Bprime_isothermal and Bprime_adiabatic as pressure derivatives
for i in isothermal adiabatic; do
    paste _output_meps/bulk_modulus_$i $outputFolder/bulk_modulus_$i _output_peps/bulk_modulus_$i | awk 'BEGIN{eps='"$eps"'}; \
        {m1=($6-$4)/eps; m2=($4-$2)/eps; m=(m2+m1)/2; print m}' > $outputFolder/Bprime_$i
done

# remove temp input and output files
rm -fr _F__* _G__* _M__* _input _T _output_peps _output_meps _volume_expansion
rm -f _T\'\' _input\'\'
echo "Successful!"


