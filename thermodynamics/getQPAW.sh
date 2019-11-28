#!/bin/bash

#-----set parameters and paths------------------------
wps=WPS_; wae=WAE_
#-----------------------------------------------------


# following 3 lines must always be present
path=`set | grep BASH_SOURCE | sed 's|.*\"\(.*\)/[^/]*\".*|\1|' | sed '/BASH_SOURCE=/s/.*/\./'`
script=`set | grep BASH_SOURCE | sed 's|.*\".*/\([^/]*\)\".*|\1|' | sed '/BASH_SOURCE=/s/.*/\./'`
options=$*; . $path/utilities/functions.include; checkOptions "-h -help -f"

# small help for options
# space between option and value must be 1
# space between option and Comment or value and Comment must be >2
h=`getOption -h`
if [ $h == True ]; then
  usage $script
  printOptions "-f   overwrite QPAW_old if existing"
  exit
fi

# detailed help
help=`getOption -help`
if [ $help == True ]; then
  details $script
  echo2 "   uses files $wae* and $wps* as obtained from splitPOTCAR" \
        "   to calculate QPAW (moments of augmentation charge)"
  echo2 "   QPAW(l,l') = \int (WAE(l,n,r)^2  -  WPS(l',n,r)^2) dr"
  echo2 "   if existing, old QPAW file is saved to QPAW_old" \
        "   and QPAW and QPAW_old are compared"
  exit
fi

n=`ls -1 WAE_[1-9]* | awk 'END{print NR}'`
if [ $n == 0 ]; then error "no WAE_* files"; fi
for (( i=1; i<=$n; i++ )) do
  check WPS_$i; check=`paste WAE_$i WPS_$i | awk 'END{if(NF!=4)print "error"}'`
  if [ "$check" == error ]; then error "WAE_$i and WPS_$i inconsistent"; fi
done

if [ -e QPAW ]; then
  overwrite=`getOption -f`
  if [ -e QPAW_old -a $overwrite != True ]; then error "QPAW and QPAW_old existing; use -f to overwrite QPAW_old"; fi
  mv QPAW QPAW_old; echo; echo old QPAW saved to QPAW_old;
fi

for (( i=1; i<=$n; i++ )) do
  for (( j=1; j<=$n; j++ )) do
    paste WAE_$i WAE_$j WPS_$i WPS_$j | \
      awk 'BEGIN{s=0};
           NR>1{d=($1-x);yy=($2*$4-$6*$8);
           f=y+(yy-y)/2;s=s+f*d};
           {x=$1;y=$2*$4-$6*$8};
           END{printf("%.10f\n",s)}' >> QPAW
  done
done


