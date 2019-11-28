#!/bin/bash

#-----set parameters and paths------------------------
min=-1; max=1; # in % around equilibrium cBya
d="0.010"; # delta aLat in Ang
#-----------------------------------------------------


# following 3 lines must always be present
path=`set | grep BASH_SOURCE | sed 's|.*\"\(.*\)/[^/]*\".*|\1|' | sed '/BASH_SOURCE=/s/.*/\./'`
script=`set | grep BASH_SOURCE | sed 's|.*\".*/\([^/]*\)\".*|\1|' | sed '/BASH_SOURCE=/s/.*/\./'`
options=$*; . $path/utilities/functions.include; checkOptions "-h -help -e -s -r -d -c"

# small help for options
# space between option and value must be 1
# space between option and Comment or value and Comment must be >2
h=`getOption -h`
if [ $h == True ]; then
  usage $script
  printNeeded  "-e elem       element(s)"
  printOptions "-r min max    range in % of cBya (default: min=$min% max=$max%)"    \
               "-d delta      delta cBya between cBya points (default: $d Ang)"       \
               "-c concs      if more elements then concentrations need to be given here"
  exit
fi

# detailed help
help=`getOption -help`
if [ $help == True ]; then
  details $script
  echo2 "   create a number of cBya points around the" \
        "   equilibrium cBya for an element specified with -e"
  echo2 "   the spacing between the cBya values can" \
        "   be controlled by the -d option"
  echo2 "   the range can be controlled by the -r option" \
        "   which expects min and max in % of the equilibrium" \
        "   cBya, i.e., min = (cByaMin-cByaEq)/cByaEq * 100 and" \
        "   max = (cByMax-cByaEq)/cByaEq * 100" 
  echo2 "   the information for cByaEq is taken from:" \
        "   $path/utilities/db_cBya"
  echo2 "   more elements can be supplied to the -e option for which"   \
        "   a linear interpolation is used to determine the cByaEq for" \
        "   a certain concentration which needs to be supplied to the"  \
        "   -c option; e.g.:" \
        "   $script -e Ti Fe -c 0.5 0.5"
  echo2 "   values are for GGA-PBE"
  exit
fi

# equilibrium volume or lattice constant?
eq=`getOption -e`; if [ $eq != True ]; then error "elem missing (-e elem)"; fi
nel=1
eq=`getValue -e`; if [ -z "$eq" ]; then error "value for -e option missing"; fi
c=`echo $eq | awk '/^[0-9.]+$/{print "aLat";exit};{print "elem"}'`
check $path/utilities/db_cBya
# check if more elements for linear mixing
nel=`echo $eq | awk '{print NF}'`
if [ $nel -gt 1 ]; then
  # for multiple elements we run some checks concerning -e and -c option
  # we need -c option with concentrations for more than one element
  conc=`getOption -c`
  if [ $conc != True ]; then error "more elements in -e option provided but no -c option with concentrations given"; fi
  conc=`getValue -c`
  # check if number of concentrations equal to number of elements and if they add up to one
  nconc=`echo $conc | awk '{print NF}'`
  if [ $nconc -ne $nel ]; then error "number of concentrations in -c option is not equal to number of elements in -e option"; fi
  totConc=`echo $conc | awk '{c=0; for (i=1;i<=NF;i++) c=c+$i; if ((c-1)^2>0.0001) print "wrong"}'`
  if [ "$totConc" == wrong ]; then error "concentrations in -c option do not add up to one"; fi
else
  # for one element we set the concentration simply to one
  conc="1"
fi

# now we loop over the elements (in $eq) and average with the concentrations (in $conc)
avg=0
for (( i=1; i<=$nel; i++ )) do
  e=`echo $eq   | awk '{print $'$i'}'`
  x=`echo $conc | awk '{print $'$i'}'`
  info=`grep -e "^$e " $path/utilities/db_cBya | awk '{print $2}'`
  if [ -z "$info" ]; then error "element $e not in $path/utilities/db_cBya"; fi 
  avg=`echo "$info" | awk 'NR==1{print $1*'$x'+'$avg'}'` # here we average bcc
done

r=`getOption -r`
if [ $r == True ]; then
  min=`getValue -r | awk '{print $1}'`
  max=`getValue -r | awk '{print $2}'`
fi

dOp=`getOption -d`
if [ $dOp == True ]; then d=`getValue -d`; fi

mincBya=`echo $avg | awk '{printf("%5.3f",(1+('$min'/100))*$1)}'`
maxcBya=`echo $avg | awk '{printf("%5.3f",(1+('$max'/100))*$1)}'`
range=`echo $avg | awk '{aStart=$1; while ((aStart-'$mincBya')^2>((aStart-'$d')-'$mincBya')^2) aStart=aStart-'$d';
  while ((aStart-'$maxcBya')^2>((aStart+'$d')-'$maxcBya')^2) {printf("%s ",aStart); aStart=aStart+'$d'}};END{printf("%s",aStart)}'`

mincByaNew=`echo $range | awk '{printf "%.3f", $1}'`
maxcByaNew=`echo $range | awk '{printf "%.3f", $NF}'`
mincByaNewP=`echo $range | awk '{printf "%.2f", ($1-'$avg')/'$avg'*100}'`
maxcByaNewP=`echo $range | awk '{printf "%.2f", ($NF-'$avg')/'$avg'*100}'`
dP=`echo $d $avg | awk '{printf("%.2f",$1/$2*100)}'`
n=`echo $range | awk '{print NF}'`

echo
echo -e "cBya eq:  \033[1m\033[31m--->\033[0m $avg \033[1m\033[31m<---\033[0m"
echo
echo "cBya min:    $mincByaNew ($mincByaNewP%)"
echo "cBya max:    $maxcByaNew ($maxcByaNewP%)"
echo "cBya step:   $d ($dP%)"
echo
echo "cBya range:  $range  (d=$d, n=$n)"


