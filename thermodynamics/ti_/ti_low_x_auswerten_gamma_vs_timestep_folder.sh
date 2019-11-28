#!/bin/sh
jobs=`ls -1d *t[0-9]*_gamma*`
job1=`ls -1d *t[0-9]*_gamma* | head -1`;
temp=`ls -1d $job1/*Ang_*K  | sed 's|K/|K|' | sed 's|.*/||' | sed 's|.*_||' | sed 's|K||'`
echo temp:$temp
#in=`echo $temp k 3/2`

## You have: 1100 k K 3/2
## You want: meV
##         * 142.18615
##         / 0.0070330339
## You have: 1 k K 3/2
## You want: meV
##         * 0.12926013
##         / 7.7363373

refene=`echo 1 | awk '{print '"$temp"'*0.12926013}'`
echo refene:$refene



rm -f jobListREST
touch jobListREST
for job in $jobs;do
#echo $job
mean="---"
thisjobs=`find -L $job -name pre_equilibration | wc -l`
firstpart=`ls -1d $job/*Ang_*K | sed 's|K/|K|' | sed 's|.*/||'`
#echo job: $job tj:$thisjob fp:$firstpart
[ "$thisjobs" = "0" ] && echo $job && ls -1d $job/*Ang_*K/* >> jobListREST && continue
mean=`tail -n+2000 $job/$firstpart/lambda0.0*/pre_equilibration | awk '{print $6}' |  sed 's|=.*||g' | grep -o "[-0-9.]*" | awk '{sum+=$NF} END{print sum/NR}' `
t=`tail -n+2000 $job/$firstpart/lambda0.0*/pre_equilibration | awk '{print $4}' | sed 's|=.*||g' | grep -o "[-0-9.]*" | awk '{sum+=$NF} END{print sum/NR}' `


echo "$job || $mean ($refene meV) || $t ($temp K)" | tee -a out

done

[ "`cat jobListREST | wc -l`" = "0" ] && rm -f jobListREST

################################
### prepare output
################################
ts=`grep -o "_t[0-9.]*_gamma" out | sed 's|_t||' | sed 's|_gamma.*||' | sort | uniq`

for t in $ts;do
gamma=`grep -o "_t[0-9.]*_gamma[0-9.]*" out | sed 's|.*gamma||' | sort | uniq`
rm -f dUdL_Uref_$t
rm -f dUdL_Uref_ref
cat out | grep "_t$t\_gamma" | sed 's|.*gamma||' | sort -n | awk '{print $1,$3}' > dUdL_Uref_$t
cat out | grep "_t$t\_gamma" | sed 's|.*gamma||' | sort -n | awk '{print $1,'"$refene"'}' > dUdL_Uref_ref
rm -f dUdL_dt_$t
rm -f dUdL_dt_ref
cat out | grep "_t$t\_gamma" | sed 's|.*gamma||' | sort -n | awk '{print $1,$7}' > dUdL_dt_$t
cat out | grep "_t$t\_gamma" | sed 's|.*gamma||' | sort -n | awk '{print $1,'"$temp"'}' > dUdL_dt_ref


done
