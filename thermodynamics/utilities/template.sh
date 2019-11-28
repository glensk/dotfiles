#!/bin/bash

#-----set parameters and paths------------------------
path1=...
parameter1=...
...
#-----------------------------------------------------


# following 3 lines must always be present
path=`set | grep BASH_SOURCE | sed 's|.*\"\(.*\)/[^/]*\".*|\1|' | sed '/BASH_SOURCE=/s/.*/\./'`
script=`set | grep BASH_SOURCE | sed 's|.*\".*/\([^/]*\)\".*|\1|' | sed '/BASH_SOURCE=/s/.*/\./'`
options=$*; . $path/utilities/functions.include; checkOptions "-h -help -o -o1 ..."

# small help for options
# space between option and value must be 1
# space between option and Comment or value and Comment must be >2
h=`getOption -h`
if [ $h == True ]; then
  usage $script
  printNeeded  "-O         Comment" \
               "-O1        Comment"    # flags that must be given
  printOptions "-o         Comment" \
               "-o1        Comment" \
               "-o2 nn     Comment"    # option with value nn
  echo "Note:    some note ..." 1>&2 
  echo "Example: $script ..." 1>&2   
  exit
fi

# detailed help
help=`getOption -help`
if [ $help == True ]; then
  details $script
  echo2 "   some details ..." \
        "   some more details ..."
  echo2 "   yet more ..."
  echo2 "   and so on..." \
        "   ..."
  exit
fi

# mathematica kernel if needed
checkAndSetMath

# create parameters.dat template and exit
genPar=`getOption -p`
if [ $genPar == True ]; then
  echo "
param1=...   # Comment
param2=...   # Comment
...
" > parameters.dat
  echo; echo "parameters.dat written"; exit
fi

# check if all input files available
input="parameters.dat INCAR KPOINTS POTCAR POSCAR"
for i in $input; do check $i; done

# get options
option=`getOption -o`; option1=`getOption -o1`; option2=`getOption -o2`;
if [ $option  == True ]; then ...; fi
if [ $option1 == True ]; then ...; fi
if [ $option2 == True ]; then option2val=`getValue -o2` fi  # getValue for options with value

# read in all parameters from parameters.dat
check parameters.dat
param1=`get param1`; param2=`get param2`; ...
checkOptions "$param1" "$param2" ...

# files preparation or other operations
...

# run mathematica
$math << EOF
<<$path/mathematica/ALL.math;
MODULE[...]
EOF

