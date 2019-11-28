#!/bin/bash

# following 3 lines must always be present
path=`set | grep BASH_SOURCE | sed 's|.*\"\(.*\)/[^/]*\".*|\1|' | sed '/BASH_SOURCE=/s/.*/\./'`
script=`set | grep BASH_SOURCE | sed 's|.*\".*/\([^/]*\)\".*|\1|' | sed '/BASH_SOURCE=/s/.*/\./'`
options=$*; . $path/utilities/functions.include; checkOptions "-h -help -r -f -c -F -ff"

if [ `getOption -ff` != True ]; then
  echo
  echo "# =========================================================================="
  echo; echo "    this script is obsolete; use VASP_zipOUTPUT.sh instead"; echo
  echo "# =========================================================================="
  exit;
fi

# small help for options
h=`getOption -h`
if [ $h == True ]; then
  usage $script
  printOptions "-r           remove folders after moving OUTCARs"         \
               "-f           force running also through unzipped OUTCARs" \
               "-c           copy the OUTCARs instead of moving them"     \
               "-F fileName  do the operations on fileName instead of OUTCAR"
  exit
fi

# detailed help
help=`getOption -help`
if [ $help == True ]; then
  details $script
  echo2 "   the script moves all OUTCARs from the 1. subfolder level to the present"  \
        "   directory appending the subfolder name to the OUTCARs"
  echo2 "   any alphabetic charaters before the numeric ones are deleted and the"     \
        "   suffix 'Ang' is as well deleted; a typical application looks like:"       \
        "     5.53Ang/OUTCAR.gz  -->  OUTCAR.5.53.gz"                                 \
        "     5.55Ang/OUTCAR.gz  -->  OUTCAR.5.55.gz"                                 \
        "   or"                                                                       \
        "     lattice.2.86/OUTCAR.gz  -->  OUTCAR.2.86.gz"                            \
        "     lattice.2.90/OUTCAR.gz  -->  OUTCAR.2.90.gz"
  echo2 "   by default script works only for zipped OUTCARs, but the -f option can"   \
        "   be used to override this (in which case the OUTCARs are zipped)"
  echo2 "   the -c option keeps the OUTCARs in the subfolders, i.e., it just copies"  \
        "   them into the present directory instead of moving (doing still the name change)"
  echo2 "   the '-F fileName' option can be used to perform the operations (i.e.,"    \
        "   copying or moving) on fileName instead of on OUTCARs"
  exit
fi

# check if we use a different file name
Fop=`getOption -F`
if [ $Fop == True ]; then
  fileName=`getValue -F`
  if [ "$fileName" == "" ]; then error "no value to -F option"; fi
else
  fileName=OUTCAR
fi

# check if we are removing the folders
remove=`getOption -r`

# check if we only copy the OUTCARs
copy=`getOption -c`

# -r and -c option do not work together
if [ $remove == True -a $copy == True ]; then error "-r and -c option cannot be used together"; fi

# we work through all OUTCARs one folder level up
l=`ls -1d */$fileName.gz 2> /dev/null`

for i in $l; do
  # deleting Ang is for consistency reasons with previous versions
  folder=`echo $i | sed 's|/'$fileName'\.gz||'`
  s=`echo $folder | sed -e 's|Ang||' -e 's|^[a-zA-Z_\.=]*||'`
  echo "$i  -->  $fileName.$s.gz"
  if [ $copy == True ]; then cp $i $fileName.$s.gz; else mv $i $fileName.$s.gz; fi
  if [ $remove == True ]; then rm -fr $folder ; fi
done

# check if we need to do the same for unzipped OUTCARs
force=`getOption -f`
if [ $force == True ]; then
  l2=`ls -1d */$fileName 2> /dev/null`
  for i in $l2; do
    # deleting Ang is for consistency reasons with previous versions
    folder=`echo $i | sed 's|/'$fileName'||'`
    s=`echo $folder | sed 's|Ang||'`
    echo "$i  -->  $fileName.$s.gz"
    if [ $copy == True ]; then cp $i $fileName.$s; else mv $i $fileName.$s; fi
    gzip $fileName.$s
    if [ $remove == True ]; then rm -fr $folder ; fi
  done
fi

if [ -z "$l" -a "$force" != True ]; then error " no */$fileName.gz files; nothing done"; fi
if [ -z "$l" -a -z "$l2" ]; then error " no */$fileName and */$fileName.gz files; nothing done"; fi

