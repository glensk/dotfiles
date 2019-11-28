#!/bin/bash

# following 3 lines must always be present
path=`set | grep BASH_SOURCE | sed 's|.*\"\(.*\)/[^/]*\".*|\1|' | sed '/BASH_SOURCE=/s/.*/\./'`
script=`set | grep BASH_SOURCE | sed 's|.*\".*/\([^/]*\)\".*|\1|' | sed '/BASH_SOURCE=/s/.*/\./'`
options=$*; . $path/utilities/functions.include; checkOptions "-h -help"

# small help for options
# space between option and value must be 1
# space between option and Comment or value and Comment must be >2
h=`getOption -h`
if [ $h == True ]; then
  usage $script
  printNeeded  "mathScript     mathScript for evaluation; must be given as FIRST argument"
  exit
fi

# detailed help
help=`getOption -help`
if [ $help == True ]; then
  details $script
  echo2 "   checks for variables in math module that are not defined in the header of module"
  echo2 "   the 'mathScript' for evaluation must be provided as the FIRST argument with"
  exit
fi

# check run script
mathscript=$1
if [ "$mathscript" == "" ]; then error "mathScript not given"; fi
check $mathscript
mm=`echo $mathscript | sed 's|.*/||' | sed 's|\.math||'`

# get input variables of the module
l=`cat $mathscript | xargs | sed -e 's/.*'$mm'[a-zA-Z]*\[//' -e 's/\].*//' -e 's/_[^,]*,/,/g' -e 's/_.*//'`

checkAndSetMath
$math << EOF
<<$mathscript
<<$path/mathematica/checkForUndefined.math
n=checkForUndefined;
n=ToExpression[n];

inp={$l};
Print["input variables:"]; Print[inp]; Print[""];

Print["notdefined variables:"];
notdefined=Complement[n,inp,{$mm,allModulesLoaded}];
If[Length[notdefined]==0,Print["NONE"],Print[notdefined]];
EOF


