#!/bin/bash

#-----set parameters and paths------------------------
min=-8; max=12; # in % around equilibrium volume
d=0.02; # delta aLat in Ang
cBya=1.632993161855452
pi=3.141592653589793
hcpFactor=`echo $cBya $pi | awk '{printf "%.10f",2/($1*sin($2/3))}'`
#-----------------------------------------------------


# following 3 lines must always be present
path=`set | grep BASH_SOURCE | sed 's|.*\"\(.*\)/[^/]*\".*|\1|' | sed '/BASH_SOURCE=/s/.*/\./'`
script=`set | grep BASH_SOURCE | sed 's|.*\".*/\([^/]*\)\".*|\1|' | sed '/BASH_SOURCE=/s/.*/\./'`
options=$*; . $path/utilities/functions.include; checkOptions "-h -help -e -V -s -r -d -c"

# small help for options
# space between option and value must be 1
# space between option and Comment or value and Comment must be >2
h=`getOption -h`
if [ $h == True ]; then
  usage $script
  printNeeded  "-e elem/aLat  element(s) OR equilibrium lattice constant in Ang"
  printOptions "-V vol        equilibrium volume can be used instead of -e aLat" \
               "-s strType    structure type (bcc/fcc/hcp); needs to be given if -e aLat" \
               "-r min max    range in % of volume (default: min=$min% max=$max%)"    \
               "-d delta      delta aLat between aLat points (default: $d Ang)"       \
               "-c concs      if more elements then concentrations need to be given here"
  exit
fi

# detailed help
help=`getOption -help`
if [ $help == True ]; then
  details $script
  echo2 "   create a number of aLat points around the" \
        "   equilibrium lattice constant for a cold curve run"
  echo2 "   the spacing between the lattice constants can" \
        "   be controlled by the -d option"
  echo2 "   the range can be controlled by the -r option" \
        "   which expects min and max in % of the equilibrium" \
        "   volume, i.e., min = (Vmin-Veq)/Veq * 100 and" \
        "   max = (Vmax-Veq)/Veq * 100" 
  echo2 "   instead of specifying the equilibrium lattice constant" \
        "   an element name can be specified with the same option" \
        "   e.g., -e Al (case insensitive, i.e., -e al works as well)," \
        "   and the corresponding aLatEq is taken from:" \
        "   $path/utilities/volumes.dat"
  echo2 "   more elements can be supplied to the -e option for which"   \
        "   a linear interpolation is used to determine the aLatEq for" \
        "   a certain concentration which needs to be supplied to the"  \
        "   -c option; e.g.:" \
        "   $script -e Ti Fe -c 0.5 0.5"
  echo2 "   values are for GGA-PBE"
  exit
fi

# equilibrium volume or lattice constant?
vOp=`getOption -V`;
eq=`getOption -e`; if [ $eq != True -a $vOp != True ]; then error "eqALat or eqVol or elem missing (-e eqALat or -e elem or -V eqVol)"; fi
if [ $eq == True -a $vOp == True ]; then error "-e and -V option cannot be used together"; fi
strType=`getOption -s`
nel=1
if [ $eq == True ]; then
  eq=`getValue -e`; if [ -z "$eq" ]; then error "value for -e option missing"; fi
  c=`echo $eq | awk '/^[0-9.]+$/{print "aLat";exit};{print "elem"}'`
  if [ $c == elem ]; then
    check $path/utilities/volumes.dat
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
    eqBcc=0; eqFcc=0
    for (( i=1; i<=$nel; i++ )) do
      e=`echo $eq   | awk '{s=$'$i'; s1=toupper(substr(s,1,1)); s2=tolower(substr(s,2)); print s1 s2}'`
      x=`echo $conc | awk '{print $'$i'}'`
      info=`grep -e ''$e'$' -e "${e}_" $path/utilities/volumes.dat | \
        awk 'BEGIN{f=0};$NF=="'$e'"{f=1;print $0,-1};END{if(f==0&&NR>0) print}'`
      if [ -z "$info" ]; then error "element $e not in $path/utilities/volumes.dat"; fi 
      eqBcc=`echo "$info" | awk 'NR==1{print $4*'$x'+'$eqBcc'}'` # here we average bcc
      eqFcc=`echo "$info" | awk 'NR==1{print $5*'$x'+'$eqFcc'}'` # here we average fcc
      orig=`echo "$info" | awk 'NR==1{print $NF}'`
      if [ "$orig" != "-1" ]; then echo; echo; echo " WARN: Using $orig instead of $e !"; fi
    done

    # transform fcc aLat to hcp aLat
    eqHcp=`echo $eqFcc | awk '{printf "%.2f", ('$hcpFactor'*($1^3/4))^(1/3) }'`
  else
    if [ $strType != True ]; then error "-e aLat but NO -s option (strType)"; fi
    eqBcc=$eq; eqFcc=$eq; eqHcp=$eq
  fi
else
  c=vol
  eq=`getValue -V`; if [ -z "$eq" ]; then error "value for -V option missing"; fi
  eqBcc=`echo $eq | awk '{printf "%.10f", (2*$1)^(1/3)}'`
  eqFcc=`echo $eq | awk '{printf "%.10f", (4*$1)^(1/3)}'`
  eqHcp=`echo $eq | awk '{printf "%.10f", ('$hcpFactor'*$1)^(1/3)}'`
fi

if [ $strType == True ]; then 
  printAll=False
  strType=`getValue -s`
  if [ "$strType" != bcc -a "$strType" != fcc -a "$strType" != hcp ]; then error "wrong or no value for -s option"; fi
else
  printAll=True
fi

r=`getOption -r`
if [ $r == True ]; then
  min=`getValue -r | awk '{print $1}'`
  max=`getValue -r | awk '{print $2}'`
fi
dOp=`getOption -d`
if [ $dOp == True ]; then d=`getValue -d`; fi

toBohr=1.88972613289636593

printOut() { s=$1; a=$2; fac=$3
eqBohr=`echo $a | awk '{printf("%.2f",$1*'$toBohr')}'`
eqVol=`echo $a | awk '{printf("%.3f",$1^3/'$fac')}'`
eqVolBohr=`echo $eqBohr | awk '{printf("%.3f",$1^3/'$fac')}'`
minVol=`echo $eqVol | awk '{printf("%.3f",(1+('$min'/100))*$1)}'`
maxVol=`echo $eqVol | awk '{printf("%.3f",(1+('$max'/100))*$1)}'`
range=`echo $a | awk '{aStart=$1; while ((aStart^3/'$fac'-'$minVol')^2>((aStart-'$d')^3/'$fac'-'$minVol')^2) aStart=aStart-'$d';
  while ((aStart^3/'$fac'-'$maxVol')^2>((aStart+'$d')^3/'$fac'-'$maxVol')^2) {printf("%s ",aStart); aStart=aStart+'$d'}};END{printf("%s",aStart)}'`
dVol=`echo $range | awk '{printf("%.3f",$2^3/'$fac'-$1^3/'$fac')}'`
dVolP=`echo $dVol | awk '{printf("%.1f",100*$1/'$eqVol')}'`
dVolMax=`echo $range | awk '{printf("%.3f",$NF^3/'$fac'-$(NF-1)^3/'$fac')}'`
dVolPMax=`echo $dVolMax | awk '{printf("%.1f",100*$1/'$eqVol')}'`
n=`echo $range | awk '{print NF}'`
echo -e "\033[1m\033[31m$s $a\033[0m"
echo -e "${s}_eqAlat:     $a Ang       $eqBohr Bohr"
echo "${s}_eqVol:    $eqVol Ang3    $eqVolBohr Bohr3"; echo;
echo "${s}_VolMin:    $min % --> $minVol Ang3"
echo "${s}_VolMax:    $max % --> $maxVol Ang3"
echo "${s}_StepMin:  $dVolP % -->  $dVol Ang3"
echo "${s}_StepMax:  $dVolPMax % -->  $dVolMax Ang3"
echo; echo "${s}_aLatRange Ang: $range  (d=$d, n=$n)"
}

# no concentrations output for one element; adjust for output if more elements
if [ $nel -eq 1 ]; then conc=""; else conc=" $conc "; fi

if [ "$strType" == bcc -o $printAll == True ]; then
  echo
  if [ $c == elem ]; then echo -n -e "\033[31m\033[1m$eq$conc \033[0m"; fi
  printOut bcc $eqBcc 2;
  echo
fi

if [ "$strType" == fcc -o $printAll == True ]; then
  echo
  if [ $c == elem ]; then echo -n -e "\033[31m\033[1m$eq$conc \033[0m"; fi
  printOut fcc $eqFcc 4;
  echo
fi

if [ "$strType" == hcp -o $printAll == True ]; then
  echo
  if [ $c == elem ]; then echo -n -e "\033[31m\033[1m$eq$conc \033[0m"; fi
  printOut hcp $eqHcp $hcpFactor;
  echo
fi


