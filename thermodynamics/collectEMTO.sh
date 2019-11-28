#!/bin/bash

# following 3 lines must always be present
path=`set | grep BASH_SOURCE | sed 's|.*\"\(.*\)/[^/]*\".*|\1|' | sed '/BASH_SOURCE=/s/.*/\./'`
script=`set | grep BASH_SOURCE | sed 's|.*\".*/\([^/]*\)\".*|\1|' | sed '/BASH_SOURCE=/s/.*/\./'`
options=$*; . $path/utilities/functions.include; checkOptions "-h -help -r -f -c"

# small help for options
h=`getOption -h`
if [ $h == True ]; then
  usage $script
  printOptions "-r    remove folders after moving OUT-EMTO files" \
               "-c    copy the OUT-EMTO files instead of moving them"
  exit
fi

# detailed help
help=`getOption -help`
if [ $help == True ]; then
  details $script
  echo2 "   the script moves all OUT-EMTO filesfrom the 1. subfolder level to the present"  \
        "   directory appending the subfolder name to the OUT-EMTO file"
  echo2 "   any alphabetic charaters before the numeric ones are deleted and the"       \
        "   suffix 'Ang' is as well deleted; a typical application looks like:"         \
        "     3.16Ang/OUT-EMTO.gz  -->  OUT-EMTO.3.16.gz"                               \
        "     3.18Ang/OUT-EMTO.gz  -->  OUT-EMTO.3.18.gz"                               \
        "   or"                                                                         \
        "     lattice.3.16/OUT-EMTO.gz  -->  OUT-EMTO.3.16.gz"                          \
        "     lattice.3.18/OUT-EMTO.gz  -->  OUT-EMTO.3.16.gz"                      
  echo2 "   the script checks whether the corresponding calculations are finished"      \
        "   giving a message if not and also not copying the corresponding log file"
  echo2 "   the -c option keeps the log files in the subfolders, i.e., it just copies"  \
        "   them into the present directory instead of moving (doing still the name change)"
  exit
fi

# check if we are removing the folders
remove=`getOption -r`

# check if we only copy the OUTCARs
copy=`getOption -c`

# -r and -c option do not work together
if [ $remove == True -a $copy == True ]; then error "-r and -c option cannot be used together"; fi

# we work through all OUTCARs one folder level up
l=`ls -1d */OUT-EMTO.gz 2> /dev/null`

for i in $l; do
  # check whether job finished
  c=`zgrep " KFCD: OK  Finished at:" $i`
  if [ "$c" = "" ]; then
    echo "$i unfinished";
    continue
  fi

  # deleting Ang is for consistency reasons with previous versions
  folder=`echo $i | sed 's|/OUT-EMTO\.gz||'`
  s=`echo $folder | sed -e 's|Ang||' -e 's|^[a-zA-Z_\.=]*||'`
  echo "$i  -->  OUT-EMTO.${s}.gz"
  if [ $copy == True ]; then cp $i OUT-EMTO.${s}.gz; else mv $i OUT-EMTO.${s}.gz; fi
  if [ $remove == True ]; then rm -fr $folder ; fi
done

if [ -z "$l" ]; then error " no */OUT-EMTO.gz files; nothing done"; fi

