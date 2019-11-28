#!/bin/bash

#-----set parameters and paths------------------------
toKeep=".in0 .in1 .in1_orig .in1_before_rerun .in2 .inc .inm .int .inst .struct .klist
   .dos1 .dos1ev .scf .scf1 .scf2 .clmsum .error .dos1 .dos1ev .qtl"
#-----------------------------------------------------

# following 3 lines must always be present
path=`set | grep BASH_SOURCE | sed 's|.*\"\(.*\)/[^/]*\".*|\1|' | sed '/BASH_SOURCE=/s/.*/\./'`
script=`set | grep BASH_SOURCE | sed 's|.*\".*/\([^/]*\)\".*|\1|' | sed '/BASH_SOURCE=/s/.*/\./'`
options=$*; . $path/utilities/functions.include; checkOptions "-h -help -d -force"

# small help for options
# space between option and value must be 1
# space between option and Comment or value and Comment must be >2
h=`getOption -h`
if [ $h == True ]; then
  usage $script
  printOptions "-d dirs   directories to work through (default present)" \
               "-force    force cleaning (see -help)"
  exit
fi

# detailed help
help=`getOption -help`
if [ $help == True ]; then
  details $script
  echo2 "   script deletes all typically unnecessary files from a wien2k run"
  echo2 "   the following files are NOT deleted (error files only if non-empty):" \
        "   $toKeep"
  echo2 "   if '-d dirs' option is not given then the script works through" \
        "   the present directory otherwise thorugh 'dirs'"
  exit
fi

# check if directories are given (-d option)
dOp=`getOption -d`
if [ $dOp == True ]; then dirs=`getValue -d`; else dirs="."; fi

# check if force cleaning option is given
force=`getOption -force`

# prepare $toKeep variable
new=""
for i in $toKeep; do
  new="$new *$i"
done
toKeep=$new

# work through all directories
dir=`pwd`
sumbef=0; sumaft=0
for d in $dirs; do
  if [ ! -d $d ]; then continue; fi
  cd $d
  bef=`du -b . | awk '{print $1}'`
  sumbef=`du -b . | awk '{print $1+'"$sumbef"'}'`
  if [ "$dirs" != "." ]; then echo -n "$d   "; fi
  # check if we are in a wien2k working folder, i.e., if there are correct input files
  case=`pwd | sed 's|.*/\([^/]*\)|\1|'`
  if [ ! -e $case.struct -o ! -e $case.in0 -o ! -e $case.in1 -o ! -e $case.in2 ]; then
    if [ $force != True ]; then
      error "directory $d seems not to be a wien2k work folder (no input files); to force cleaning use -force option"
    fi
  fi
  mkdir _tmp_dir; mv $toKeep _tmp_dir/ 2> /dev/null; rm -f * 2> /dev/null; mv _tmp_dir/* .; rm -fr _tmp_dir

  # delete empty error files
  l=`ls *.error 2> /dev/null`
  for i in $l; do if [ ! -s $i ]; then rm $i; fi; done

  # print the memory reduction
  aft=`du -b . | awk '{print $1}'`
  sumaft=`du -b . | awk '{print $1+'"$sumaft"'}'`
  bef=`echo $bef | awk '$1>=1024*1024{printf("%dM",$1/1024/1024)}; $1<1024*1024{printf("%dK",$1/1024)}'`
  aft=`echo $aft | awk '$1>=1024*1024{printf("%dM",$1/1024/1024)}; $1<1024*1024{printf("%dK",$1/1024)}'`
  echo "${bef}  -->  ${aft}"

  cd $dir
done

sumbef=`echo $sumbef | awk '$1>=1024*1024{printf("%dM",$1/1024/1024)}; $1<1024*1024{printf("%dK",$1/1024)}'`
sumaft=`echo $sumaft | awk '$1>=1024*1024{printf("%dM",$1/1024/1024)}; $1<1024*1024{printf("%dK",$1/1024)}'`
echo; echo " TOTAL: ${sumbef}  -->  ${sumaft}"

