#!/bin/sh

# following 3 lines must always be present
path=`set | grep BASH_SOURCE | sed 's|.*\"\(.*\)/[^/]*\".*|\1|' | sed '/BASH_SOURCE=/s/.*/\./'`
path="$path/../../"
script=`set | grep BASH_SOURCE | sed 's|.*\".*/\([^/]*\)\".*|\1|' | sed '/BASH_SOURCE=/s/.*/\./'`
options=$*; . $path/utilities/functions.include; checkOptions

l=`ls -1d lambda*/ | wc -l`; if [ "$l" == 0 ]; then error "no lambda* folders available"; fi
l=`ls -1d lambda*/`; dir=`pwd`
echo; echo "# lambda  jDist jRadius"
echo "# lambda  jDist jRadius" > jumps
for i in $l; do
  lambda=`echo $i | sed 's/lambda\(.*\)_.*/\1/'`
  cd $i
  $path/extractPOSITIONS.sh; $path/fortran/getJump.x; rm POSITIONs cell atoms_volume_steps
  j=`awk 'BEGIN{s=0};{s=s+$1};END{print s}' jumps_dist`
  jr=`awk 'BEGIN{s=0};{s=s+$1};END{print s}' jumps_radius`
  # rm jumps_radius
  cd $dir
  echo $lambda $j $jr >> jumps
  echo "    $lambda     $j     $jr"
done

awk 'BEGIN{l=-1};l==$1{s=s+$2;s2=s2+$3};l!=$1&&l!=-1{print l,s,s2};l!=$1{l=$1;s=$2;s2=$3};END{print l,s,s2}' jumps > jumps_avg

