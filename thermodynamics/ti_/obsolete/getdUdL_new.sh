#!/bin/bash

corSteps=15

l=`ls -1d lambda*`
echo "# lambda  <dUdL>(meV) errUnc(meV)     stdDev(meV)  <T>(K)  stdDev(K)  steps   time(ps)" > avg_dUdL

for i in $l; do
  last=`tail -n1 $i/dUdL`
  t=`echo $last | awk '{printf "%8.1f", $4}'`            # mean temperature
  dUdL=`echo $last | awk '{printf "%6.2f", $9}'`         # mean dUdL value
  steps=`echo $last | awk '{printf "%7d", $1}'`          # total number of steps
  timeTot=`echo $last | awk '{printf "%7.3f", $2/1000}'` # total time

  # standard deviation dUdL and T
  std=`awk 'BEGIN {s=0;sT=0;c=0}
            NR>1{s=s+($7-'"$dUdL"')^2; sT=sT+($3-'"$t"')^2; c=c+1}
            END {printf "%6.2f %6.2f", sqrt(s/(c-1)), sqrt(sT/(c-1))}' $i/dUdL`

  stddUdL=`echo $std | awk '{printf "%6.2f", $1}'`
  stdT=`echo $std | awk '{printf "%6.2f", $2}'`

  err=`echo $stddUdL | awk '{printf "%6.2f", $1/sqrt('"$steps"'/'$corSteps')}'`

  lam=`echo $i | sed 's|lambda\(.*\)/*|\1|' | awk '{printf "%5.2f",$1}'`
  echo "   $lam   $dUdL   $err             $stddUdL   $t   $stdT   $steps   $timeTot" >> avg_dUdL
done

