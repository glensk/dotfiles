#!/usr/bin/env bash
set -u

RESULT_FILE=results.log
if [ -f $RESULT_FILE ]; then
  echo "Result file $RESULT_FILE already exists!"
  exit 1
fi
RESULT_FILE_HEADER="DISTANCE,Energy,Entropy,PXX,PYY,PZZ,PXY,PYZ,PZX"

rybohr3_2_kbar () {
  RYBOHR3_2_KBAR="1.4710516 * 10 ^ 4"
  echo "$1 * $RYBOHR3_2_KBAR" | bc -l
}

echo $RESULT_FILE_HEADER >> $RESULT_FILE

for WORKDIR in ./calculations/*;
  do
  DISTANCE=$(echo $(basename $WORKDIR) | cut -d _ -f 2)

  ERROR_FLAG=$(grep -i "ERROR" $WORKDIR/_scheduler-stdout.txt)
  if [ ! -z "$ERROR_FLAG" ]; then
    echo $WORKDIR has encourntered an error!
    continue
  fi
  FINISHED=$(grep "JOB DONE" $WORKDIR/aiida.out)
  if [ -z "$FINISHED" ] && [ ! $CHECK_FINISHED -eq 0 ]; then
    echo $WORKDIR has not yet converged!
    continue
  fi

  ENERGY_RY=$(grep -iE "! +total +energy" $WORKDIR/aiida.out | tail -n 1 | awk {'print $5'})
  RY2EV=13.605698066
  ENERGY=$(echo "$ENERGY_RY * $RY2EV"| bc)
  ENTROPY_RY=$(grep "smearing contrib" $WORKDIR/aiida.out | awk {'print $5'})
  ENTROPY=$(echo "$ENTROPY_RY * $RY2EV"| bc)
  VOLUME=$(grep -i "new unit-cell volume" $WORKDIR/aiida.out | tail -n 1 | awk {'print $8'})
  STRESS_ROW1=$(grep -iE -A4 "total +stress" $WORKDIR/aiida.out | head -n2 | tail -n1)
  STRESS_ROW2=$(grep -iE -A4 "total +stress" $WORKDIR/aiida.out | head -n3 | tail -n1)
  STRESS_ROW3=$(grep -iE -A4 "total +stress" $WORKDIR/aiida.out | head -n4 | tail -n1)
  PXX=$(echo $STRESS_ROW1 | awk {'print $1'})
  PYY=$(echo $STRESS_ROW2 | awk {'print $2'})
  PZZ=$(echo $STRESS_ROW3 | awk {'print $3'})
  PXY=$(echo $STRESS_ROW1 | awk {'print $2'})
  PZX=$(echo $STRESS_ROW1 | awk {'print $3'})
  PYZ=$(echo $STRESS_ROW2 | awk {'print $3'})
  PXX=$(rybohr3_2_kbar $PXX)
  PYY=$(rybohr3_2_kbar $PYY)
  PZZ=$(rybohr3_2_kbar $PZZ)
  PXY=$(rybohr3_2_kbar $PXY)
  PZX=$(rybohr3_2_kbar $PZX)
  PYZ=$(rybohr3_2_kbar $PYZ)
  TIME=$(grep -i time $WORKDIR/aiida.out | tail -n 1 | awk {'print $9'})

  echo "$DISTANCE,$ENERGY,$ENTROPY,$PXX,$PYY,$PZZ,$PXY,$PYZ,$PZX" >> $RESULT_FILE
done
