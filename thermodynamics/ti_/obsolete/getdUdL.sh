#!/bin/sh

corSteps=15

math=/opt/mathematica/8.0.4/math
if [ ! -e $math ]; then
  echo "" 1>&2
  echo -e "\033[31m\033[1mERROR\033[0m: $math not available" 1>&2
  exit 1
fi

l=`ls -1d lambda*/ | wc -l`
if [ "$l" == 0 ]; then
  echo "" 1>&2
  echo -e "\033[31m\033[1mERROR\033[0m: no lambda* folders available" 1>&2
  exit 1
fi

if [ $# -lt 1 ]; then
  echo "" 1>&2
  echo -e "\033[31m\033[1mUSAGE\033[0m:" 1>&2
  echo -e " getdUdL.sh lambda1 offset1 [lambda2 offset2 ...]    (given lambdas with given offsets)" 1>&2
  echo -e " getdUdL.sh fileWithValues                           (format:  1.row: lambda1 offset1 )" 1>&2
  echo -e "                                                     (        [2.row: lambda2 offset2])" 1>&2
  echo -e "                                                     (        [...]                   )" 1>&2
  exit 1
fi

if [ $# -gt 1 ]; then
  s=`echo $* | xargs -n2 | awk '{print $1"_"$2}'`
  folders=""
  for i in $s; do
    lambda=`echo $i | sed 's/\(.*\)_\(.*\)/\1/'`
    offset=`echo $i | sed 's/\(.*\)_\(.*\)/\2/'`
    str=`ls -1d lambda$lambda\_* | awk '{print '$offset'"offset"$1}'`
    folders="$folders $str"
  done
else
  if [ ! -e $1 ]; then
    echo "" 1>&2
    echo -e "\033[31m\033[1mERROR\033[0m: file $1 with lambda and offset values missing" 1>&2
    exit 1
  fi
  s=`cat $1 | awk '{print $1"_"$2}'`
  folders=""
  for i in $s; do
    lambda=`echo $i | sed 's/\(.*\)_\(.*\)/\1/'`
    offset=`echo $i | sed 's/\(.*\)_\(.*\)/\2/'`
    str=`ls -1d lambda$lambda\_* 2> /dev/null | awk 'NF==1{print '$offset'"offset"$1}'`
    folders="$folders $str"
  done
fi

echo
echo -e "\033[31m\033[1mWARNING\033[0m: default correlation length of \033[31m\033[1m$corSteps\033[0m steps used (for error only)"
echo

getAvgAndDev() {
    rm -f tmp2
    $math >> tmp2 << EOF
in=Import["tmp","Table"];
time=Transpose[in][[1]]; temp=Transpose[in][[2]]; dUdL=Transpose[in][[3]];
time=time[[-1]]-time[[1]];
devTemp=StandardDeviation[temp]; devdUdL=StandardDeviation[dUdL];
errdUdL=devdUdL/Sqrt[dUdL//Length];
uncorTemp=temp[[1;;-1;;$corSteps]]; uncordUdL=dUdL[[1;;-1;;$corSteps]];
devTemp2=StandardDeviation[uncorTemp]; devdUdL2=StandardDeviation[uncordUdL];
errdUdL2=devdUdL2/Sqrt[uncordUdL//Length];
Print["dUdLCor ",Mean[dUdL], " ",devdUdL," ",errdUdL," ",errdUdL/2," ",Mean[temp]," ",devTemp," ",Length[dUdL]," ",time/1000]
Print["dUdLUnc ",Mean[uncordUdL]," ",devdUdL2," ",errdUdL2," ",errdUdL2/2," ",Mean[uncorTemp]," ",devTemp2," ",Length[uncordUdL]]
Print["dUdL ",   Mean[dUdL]," ",errdUdL2/2," ",devdUdL," ",Mean[temp]," ",devTemp," ",Length[dUdL]," ",time/1000]
EOF
}

rm -f tmp_lambda*_offset* tmp_time_lambda*_offset*
if [ -e dUdL ]; then mv dUdL dUdL_old; string="; saving old dUdL to dUdL_old"; else string=""; fi
echo "creating dUdL (check dUdL_all for details$string)"; echo
echo "# correlation length $corSteps" > dUdL
echo "# correlation length $corSteps" > dUdL_all
echo "# all energies per atom"
echo "# all energies per atom" >> dUdL
echo "# all energies per atom" >> dUdL_all
echo "# lambda  <dUdL>(meV) errUnc/2(meV)     stdDev(meV)  <T>(K)  stdDev(K)  steps   time(ps) offset"
echo "# lambda  <dUdL>(meV) errUnc/2(meV)     stdDev(meV)  <T>(K)  stdDev(K)  steps   time(ps) offset" >> dUdL
echo -n "# lambda seed offset         <dUdL>(meV) stdDev(meV) err(meV) err/2(meV)   <T>(K) stdDev(K)   steps  time(ps)      " >> dUdL_all
echo "Uncor: <dUdL>(meV) stdDev(meV) err(meV) err/2(meV)   <T>(K) stdDev(K)   steps" >> dUdL_all
for ff in $folders; do
  offset=`echo $ff | sed 's/\(.*\)offset\(.*\)/\1/' | awk '{printf("%3d",$1)}'`
  f=`echo $ff | sed 's/\(.*\)offset\(.*\)/\2/'`
  if [ -e $f/dUdL ]; then
     file=$f/dUdL
  else
    if [ -e $f/workDir/dUdL ]; then
      file=$f/workDir/dUdL
    else
      file=error
    fi
  fi
  if [ $file == "error" ]; then
    echo $f no dUdL
  else
    l=`echo $f | sed 's/lambda\(.*\)_\(.*\)/\1/'`
    s=`echo $f | sed 's/lambda\(.*\)_\(.*\)/\2/' | awk '{printf("%5d",$1)}'`
    awk 'NR>'"$offset"'+1{print $2,$3,$7}' $file > tmp
    str=`echo $offset | awk '{print $1}'`
    cat tmp >> tmp_lambda$l\_offset$str
    nUncor=""
    while [ "$nUncor" == "" ]; do
      getAvgAndDev
      nUncor=`grep dUdLUnc tmp2 | awk '{printf("%d",$7)}'`
    done
    if [ "$nUncor" -gt 2 ]; then
      dUdLall1=`grep dUdLCor tmp2 | awk '{printf( "%8.2f %9.2f %9.2f %9.2f %12.1f %8.1f %8d %9.3f",$3,$4,$5,$6,$7,$8,$9,$10)}'`
      dUdLall2=`grep dUdLUnc tmp2 | awk '{printf( "%22.2f %9.2f %9.2f %9.2f %12.1f %8.1f %8d",$3,$4,$5,$6,$7,$8,$9)}'`
      dUdLall="$dUdLall1$dUdLall2"
      dUdL=`grep "dUdL " tmp2 | awk '{printf("%6.2f %9.2f %19.2f %10.1f %9.1f %8d %9.3f", $3,$4,$5,$6,$7,$8,$9)}'`
      grep "dUdL " tmp2 | awk '{printf("%.5f\n",$9)}' >> tmp_time_lambda$l\_offset$str
      echo "   $l  $s  $offset           $dUdLall" >> dUdL_all
      echo "    $l    $dUdL    $offset" >> dUdL
      echo "    $l    $dUdL    $offset"
    else
      echo "#   $l  seed: $s  too few steps: all `awk 'END{print NR-1}' $file`  off: `cat tmp | wc -l`  unc: $nUncor" >> dUdL_all
      echo "    $l  seed: $s  too few steps: all `awk 'END{print NR-1}' $file`  off: `cat tmp | wc -l`  unc: $nUncor"
    fi
  fi
done

if [ -e avg_dUdL ]; then mv avg_dUdL avg_dUdL_old; string=" (saving old avg_dUdL to avg_dUdL_old)"; else string=""; fi
echo; echo creating avg_dUdL"$string"; echo
echo "# correlation length $corSteps" > avg_dUdL
echo "# all energies per atom" >> avg_dUdL
echo "# lambda  <dUdL>(meV) errUnc/2(meV)     stdDev(meV)  <T>(K)  stdDev(K)  steps   time(ps) offset"
echo "# lambda  <dUdL>(meV) errUnc/2(meV)     stdDev(meV)  <T>(K)  stdDev(K)  steps   time(ps) offset" >> avg_dUdL
files=`ls tmp_lambda*_offset* 2> /dev/null`
for f in $files; do
  l=`echo $f | sed 's/tmp_lambda\(.*\)_offset\(.*\)/\1/'`
  offset=`echo $f | sed 's/tmp_lambda\(.*\)_offset\(.*\)/\2/' | awk '{printf("%3d",$1)}'`
  mv $f tmp
  nUncor=""
  while [ "$nUncor" == "" ]; do
    getAvgAndDev
    nUncor=`grep dUdLUnc tmp2 | awk '{printf("%d",$7)}'`
  done
  dUdL=`grep "dUdL " tmp2 | awk '{printf("%6.2f %9.2f %19.2f %10.1f %9.1f %8d", $3,$4,$5,$6,$7,$8)}'`
  str=`echo $offset | awk '{print $1}'`
  time=`awk 'BEGIN {s=0}; NF==1{s=s+$1}; END{printf("%9.3f",s)}' tmp_time_lambda$l\_offset$str`
  echo "    $l    $dUdL $time    $offset" >> avg_dUdL
  echo "    $l    $dUdL $time    $offset"
done

rm -f tmp tmp2 tmp_time_lambda*_offset*
