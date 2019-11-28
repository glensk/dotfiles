#!/bin/sh

# this skript creates *Ang_*K/dUdL_all (and *Ang_*K/lambda*/avg_dUdL.dat)
# als input brauch dieses skript das dUdL_allinfonojums file
#oldest_file_input=`ls -la lambda*/dUdLallInfo_noJumps -rtl | head -1 | awk '{print $NF}'`
#newest_file_output=`ls -la dUdL_all lambda*/avg_dUdL.dat -rtl | tail -1 | awk '{print $NF}'`
#
#u_oldest_file_input=`stat -c %Y $oldest_file_input`
#u_newest_file_output=`stat -c %Y $newest_file_output`
#echo old:$u_oldest_file_input
#echo new:$u_newest_file_output
#diff=`expr $u_newest_file_output - $u_oldest_file_input`
#echo diff:$diff
#exit



#############################################################################################
## define default values
#############################################################################################
corSteps=15
minsteps=70
offset=80                   ## der offset wird von der dUdLallInfo_noJumps weggeschnitten
                            ## bestimmt auch wieviele schritte mindestens da sein muessen damit der lambdax_XX lauf gewertet wird
offset=30
getdUdLs=getdUdLs             ## foler where all this suff will be saved wird direkt im angK folder gemacht!

#############################################################################################
## check if mathematica folder exists
#############################################################################################
#math=/opt/mathematica/8.0/bin/math
#if [ ! -e $math ]; then
#  echo "" 1>&2
#  echo -e "\033[31m\033[1mERROR\033[0m: $math not available" 1>&2
#  exit 1
#fi

#############################################################################################
## check which lambda folder exist
#############################################################################################
l=`ls -1d lambda*/ | wc -l | sed 's|[ ]*||g'`
[ "$l" = 0 ] && echo "" 1>&2 && echo -e "\033[31m\033[1mERROR\033[0m: no lambda* folders available" 1>&2 && exit
folders=`ls -1d lambda*_[0-9]*`
[ "`echo $folders | wc -w | sed 's|[ ]*||g'`" = 0 ] && echo "" 1>&2 && echo -e "\033[31m\033[1mERROR\033[0m: no lambda* folders available" 1>&2 && exit
dUdLscheck=`ls -1d lambda*_[0-9]*/$dUdLs`
[ "`echo $dUdLscheck | wc -w | sed 's|[ ]*||g'`" = 0 ] && echo "" 1>&2 && echo -e "\033[31m\033[1mERROR\033[0m: no lambda*/$dUdLs folders available" 1>&2 && exit

#echo
#echo -e "\033[31m\033[1mWARNING\033[0m: default correlation length of \033[31m\033[1m$corSteps\033[0m steps used (for error only)"
#echo


#############################################################################################
## difine necessary mathematica scripts
#############################################################################################
## gets tmp2; list like: 
## line1:310.0 264.1 -1.30
## line2:320.0 269.8 -1.38


## mean_dUdL=`awk '{sum+=$3} {print sum/NR}' $file | tail -1`
## mean_dUdL_uncor=`awk 'NR%$corSteps==1' $file | awk '{sum+=$3} {print sum/NR}' | tail -1`

## mean_temp=`awk '{sum+=$2} {print sum/NR}' $file | tail -1`
## mean_temp_uncor=`awk 'NR%$corSteps==1' $file | awk '{sum+=$2} {print sum/NR}' | tail -1`

## dev_dUdL=`awk '{d3=$3-avg3;avg3+=d3/NR;m3+=d3*($3-avg3)}{print sqrt(m3/NR)}' $file | tail -1`
## dev_dUdL_uncor=`awk 'NR%$corSteps==1' $file | awk '{d3=$3-avg3;avg3+=d3/NR;m3+=d3*($3-avg3)}{print sqrt(m3/NR)}' | tail -1`

## error_dUdL=`awk '{d3=$3-avg3;avg3+=d3/NR;m3+=d3*($3-avg3)}{print sqrt(m3/NR)/sqrt(NR)}' $file | tail -1`
## error_dUdL_uncor=`awk 'NR%$corSteps==1' $file | awk '{d3=$3-avg3;avg3+=d3/NR;m3+=d3*($3-avg3)}{print sqrt(m3/NR)/sqrt(NR)}' | tail -1`



getAvgAndDev() {
rm -f tmp2
## hier wird alles aus der file tmp gelesen
dt=`head -2 tmp | awk '{print $1}' | xargs | awk '{print $2-$1}'`
lines_dUdL=`awk 'END {print NR}' tmp`
lines_dUdL_uncor=`awk 'NR%'"$corSteps"'==1' tmp | awk 'END {print NR}'`

mean_dUdL=` awk '{sum+=$3} {print sum/NR}' tmp | tail -1`
mean_dUdL_uncor=`awk 'NR%'"$corSteps"'==1' tmp | awk '{sum+=$3} {print sum/NR}' | tail -1`

mean_temp=` awk '{sum+=$2} {print sum/NR}' tmp | tail -1`
mean_temp_uncor=`awk 'NR%'"$corSteps"'==1' tmp | awk '{sum+=$2} {print sum/NR}' | tail -1`

dev_dUdL=`awk '{d3=$3-avg3;avg3+=d3/NR;m3+=d3*($3-avg3)}{print sqrt(m3/NR)}' tmp | tail -1`
dev_dUdL_uncor=`awk 'NR%'"$corSteps"'==1' tmp | awk '{d3=$3-avg3;avg3+=d3/NR;m3+=d3*($3-avg3)}{print sqrt(m3/NR)}' | tail -1`

dev_temp=`awk '{d2=$2-avg2;avg2+=d2/NR;m2+=d2*($2-avg2)}{print sqrt(m2/NR)}' tmp | tail -1`
dev_temp_uncor=`awk 'NR%'"$corSteps"'==1' tmp | awk '{d2=$2-avg2;avg2+=d2/NR;m2+=d2*($2-avg2)}{print sqrt(m2/NR)}' | tail -1`

error_dUdL=`awk '{d3=$3-avg3;avg3+=d3/NR;m3+=d3*($3-avg3)}{print sqrt(m3/NR)/sqrt(NR)}' tmp | tail -1`
error_dUdL_uncor=`awk 'NR%'"$corSteps"'==1' tmp | awk '{d3=$3-avg3;avg3+=d3/NR;m3+=d3*($3-avg3)}{print sqrt(m3/NR)/sqrt(NR)}' | tail -1`
###  error/dUdL/2 is inserted -----------------------------------V
echo dUdLCor $mean_dUdL       $dev_dUdL         $error_dUdL        $mean_temp       $dev_temp       | awk '{printf "%s %.2f %.2f %.2f %.2f %.2f %.2f %.0f  %.2f\n"   ,$1,$2,$3,$4,$4/2,$5,$6,'"$lines_dUdL"','"$dt"'*'"$lines_dUdL"'/1000}' >> tmp2
echo dUdLUnc $mean_dUdL_uncor $dev_dUdL_uncor   $error_dUdL_uncor  $mean_temp_uncor $dev_temp_uncor | awk '{printf "%s %.2f %.2f %.2f %.2f %.2f %.2f %.0f \n"        ,$1,$2,$3,$4,$4/2,$5,$6,'"$lines_dUdL_uncor"'}' >> tmp2
echo dUdL    $mean_dUdL       $dev_dUdL         $error_dUdL_uncor  $mean_temp       $dev_temp       | awk '{printf "%s    %.2f %.2f %.2f %.2f %.2f %.2f %.0f  %.2f\n",$1,$2,$3,$4,$4/2,$5,$6,'"$lines_dUdL"','"$dt"'*'"$lines_dUdL"'/1000}' >> tmp2
#    $math >> tmp2 << EOF
#in=Import["tmp","Table"];
#time=Transpose[in][[1]]; temp=Transpose[in][[2]]; dUdL=Transpose[in][[3]];
#time=time[[-1]]-time[[1]];
#devTemp=StandardDeviation[temp]; devdUdL=StandardDeviation[dUdL];
#errdUdL=devdUdL/Sqrt[dUdL//Length];
#uncorTemp=temp[[1;;-1;;$corSteps]]; uncordUdL=dUdL[[1;;-1;;$corSteps]];
#devTemp2=StandardDeviation[uncorTemp]; devdUdL2=StandardDeviation[uncordUdL];
#errdUdL2=devdUdL2/Sqrt[uncordUdL//Length];
#Print["dUdLCor ",Mean[dUdL], " ",devdUdL," ",errdUdL," ",errdUdL/2," ",Mean[temp]," ",devTemp," ",Length[dUdL]," ",time/1000]
#Print["dUdLUnc ",Mean[uncordUdL]," ",devdUdL2," ",errdUdL2," ",errdUdL2/2," ",Mean[uncorTemp]," ",devTemp2," ",Length[uncordUdL]]
#Print["dUdL ",   Mean[dUdL]," ",errdUdL2/2," ",devdUdL," ",Mean[temp]," ",devTemp," ",Length[dUdL]," ",time/1000]
#EOF
}


#############################################################################################
## create im lambda folder: dUdL_all (und avg_dUdL.dat)
#############################################################################################
echo
echo "#############################################################################################"
echo "## create dUdL_all"
echo "#############################################################################################"
rm -f tmp_lambda*_offset* tmp_time_lambda*_offset*
rm -rf $getdUdLs  ## notwendig da alles immer neu gemacht werden soll
[ ! -e "$getdUdLs" ] && mkdir $getdUdLs
[ -e "$getdUdLs/dUdL_all" ] && mv $getdUdLs/dUdL_all $getdUdLs/dUdL_all_old
#rm -f $getdUdLs/dUdL $getdUdLs/dUdL_all $getdUdLs/dUdL_all_old $getdUdLs/tmp_time_lambda* $getdUdLs/dUdL_lambda*_all_*_old
#[ "`find -L getdUdLs -maxdepth 1 -mindepth 1 -type f -name "dUdL_lambda*_all_fre" | wc -w | sed 's|[ ]*||g'`" != "0" ] && mmv "$getdUdLs/dUdL_lambda*_all_fre" "$getdUdLs/dUdL_lambda#1_all_fre_old"  ##2&>/dev/null
#[ "`find -L getdUdLs -maxdepth 1 -mindepth 1 -type f -name "dUdL_lambda*_all_ene" | wc -w | sed 's|[ ]*||g'`" != "0" ] && mmv "$getdUdLs/dUdL_lambda*_all_ene" "$getdUdLs/dUdL_lambda#1_all_ene_old"  ##2&>/dev/null
#[ "`find -L getdUdLs -maxdepth 1 -mindepth 1 -type f -name "dUdL_lambda*_all_eS0" | wc -w | sed 's|[ ]*||g'`" != "0" ] && mmv "$getdUdLs/dUdL_lambda*_all_eS0" "$getdUdLs/dUdL_lambda#1_all_eS0_old"  ##2&>/dev/null
echo "# correlation length $corSteps" | tee -a $getdUdLs/dUdL | tee -a $getdUdLs/dUdL_all
echo "# all energies per atom" | tee -a $getdUdLs/dUdL | tee -a $getdUdLs/dUdL_all
echo "# lambda  <dUdL>(meV) errUnc/2(meV)     stdDev(meV)  <T>(K)  stdDev(K)  steps   time(ps) offset" | tee -a $getdUdLs/dUdL
echo -n "# lambda seed offset         <dUdL>(meV) stdDev(meV) err(meV) err/2(meV)   <T>(K) stdDev(K)   steps  time(ps)      " >> $getdUdLs/dUdL_all
echo "Uncor: <dUdL>(meV) stdDev(meV) err(meV) err/2(meV)   <T>(K) stdDev(K)   steps     job_ID" >>  $getdUdLs/dUdL_all


for f in $folders; do  ## f={lambda0.0_13872,lambda0.0_13873,..., jeder einzelne lmabda folder}
    l=`echo $f | sed 's/lambda\(.*\)_\(.*\)/\1/'`
    s=`echo $f | sed 's/lambda\(.*\)_\(.*\)/\2/' | awk '{printf("%5d",$1)}'`
    id=`OUTCAR_ID.sh $f`
   
    ## hier bekommt er die dUdLallInfo_noJumps als inptu (==file)
    ## hier bekommt er die dUdLallInfo_noJumps als inptu (==file)
    file=error
    dUdLs=dUdLallInfo_noJumps            ## input file
    #[ "`frompath_elestoich.sh`" = "Si" ] && dUdLs=dUdLallInfo
    [ -e "$f/$dUdLs" ] && file=$f/$dUdLs ## file={lambda0.0_13872/dUdLallInfo_noJumps,lambda0.0_13873/dUdLallInfo_noJumps, ...}
    [ -e $f/workDir/$dUdLs ] && file=$f/workDir/$dUdLs
    [ "$file" = "error" ] && [ "`echo $* | grep -o "\-ns"`" != "-ns" ] && echo $f no $dUdLs file
    [ "$file" = "error" ] && continue
    steps=`wc -l $file | awk '{print $1-1}'`
    [ "$steps" -le "$offset" ] && echo "   $l  $s $file does not have enough steps (steps=$steps)" | tee -a $getdUdLs/dUdL_all && continue
    [ "$steps" -le "$minsteps" ] && echo "   $l  $s $file does not have enough steps (steps=$steps)" | tee -a $getdUdLs/dUdL_all && continue

    uUncor=0
 
    ## create f/getAvgAndDev_{fre,ene,eS0}_input file=lambda0.5_15355/dUdLallInfo_noJumps f={lambda0.0_13872,lambda0.0_13873,...}
    ## create f/getAvgAndDev_{fre,ene,eS0}_input file=lambda0.5_15355/dUdLallInfo_noJumps
    ## create f/getAvgAndDev_{fre,ene,eS0}_input file=lambda0.5_15355/dUdLallInfo_noJumps

    list="fre ene eS0"
    for var in $list; do
            rm -f tmpfile
            ## hinzufuegen (tee -a) zu getdUdL folder ---------------------------------------->>>       v v v v v v v v v v v v v v v v v 
            [ "$var" = "fre" ] && awk 'NR>'"$offset"'+1{print $2,$3,$13}' $file | tee -a tmpfile | tee -a $getdUdLs/dUdL_lambda$l\_all_$var  > /dev/null ## =? f/getAvgAndDev_input_fre
            [ "$var" = "ene" ] && awk 'NR>'"$offset"'+1{print $2,$3,$14}' $file | tee -a tmpfile | tee -a $getdUdLs/dUdL_lambda$l\_all_$var  > /dev/null ## =? f/getAvgAndDev_input_ene
            [ "$var" = "eS0" ] && awk 'NR>'"$offset"'+1{print $2,$3,$15}' $file | tee -a tmpfile | tee -a $getdUdLs/dUdL_lambda$l\_all_$var  > /dev/null ## =? f/getAvgAndDev_input_eS0
            #[ -e $f/getAvgAndDev_input_$var ] && [ "`diff tmpfile $f/getAvgAndDev_input_$var`" = "" ] && rm -f tmpfile && continue

            ## ueberschreiben lambda0.5_12345 folder --->>> v v v v v v v v v
            mv tmpfile $f/getAvgAndDev_input_$var
            cp $f/getAvgAndDev_input_$var tmp;
            #echo xxxxxxx
            #pwd
            #echo xxxxxxx


            ##################################################################
            ##################################################################
            getAvgAndDev
            ##################################################################
            ##################################################################
            #exit
            cp tmp2 $f/getAvgAndDev_$var
            rm -f tmpfile tmp2
            done

    ## create lambda0.5_12345/avg_dUdL.dat
    rm -f $f/avg_dUdL.dat
    cat $f/getAvgAndDev_input_fre | awk '{print $3}' | awk '{sum+=$NF} {print sum/NR}' > $f/avg_dUdL.dat


    ### add to dUdL_all (dUdL all just contains free energy, not ewe and eS0)
    ### add to dUdL_all (dUdL all just contains free energy, not ewe and eS0)
    nUncor=`grep dUdLUnc $f/getAvgAndDev_fre | awk '{printf("%d",$7)}'`  ### das macht keinen sinn!!!!!
    nUncor=`grep dUdLUnc $f/getAvgAndDev_fre | awk '{printf("%d",$8)}'`  ### das macht keinen sinn!!!!!

    #echo nUncor: $nUncor
    #echo $f/getAvgAndDev_fre
    if [ "$nUncor" -gt 2 ]; then
      #dUdLall1=`grep dUdLCor $f/getAvgAndDev_fre | awk '{printf( "%8.2f %9.2f %9.2f %9.2f %12.1f %8.1f %8d %9.3f",$3,$4,$5,$6,$7,$8,$9,$10)}'` ### das macht auch keinen sinn!
      dUdLall1=`grep dUdLCor $f/getAvgAndDev_fre | awk '{printf( "%8.2f %9.2f %9.2f %9.2f %12.1f %8.1f %8d %9.3f",$2,$3,$4,$5,$6,$7,$8,$9)}'`
      #echo dUdLall1: $dUdLall1
      #dUdLall2=`grep dUdLUnc $f/getAvgAndDev_fre | awk '{printf( "%22.2f %9.2f %9.2f %9.2f %12.1f %8.1f %8d",$3,$4,$5,$6,$7,$8,$9)}'` ### macht auch keinen sinn
      dUdLall2=`grep dUdLUnc $f/getAvgAndDev_fre | awk '{printf( "%22.2f %9.2f %9.2f %9.2f %12.1f %8.1f %8d",$2,$3,$4,$5,$6,$7,$8)}'` 
      #echo dUdLall2: $dUdLall2
      dUdLall="$dUdLall1$dUdLall2"
      #dUdL=`grep "dUdL " $f/getAvgAndDev_fre | awk '{printf("%6.2f %9.2f %19.2f %10.1f %9.1f %8d %9.3f", $3,$4,$5,$6,$7,$8,$9)}'` ### macht keinen sinn
      dUdL=`grep "dUdL " $f/getAvgAndDev_fre | awk '{printf("%6.2f %9.2f %19.2f %10.1f %9.1f %8d %9.3f", $2,$5,$3,$6,$7,$8,$9)}'`
      #echo dUdL: $dUdL
      grep "dUdL " $f/getAvgAndDev_fre | awk '{printf("%.5f\n",$9)}' >> $getdUdLs/tmp_time_lambda$l
      ## $s = seed = 17551
      ## dUdLall = 8.11 0.15 0.07 936.05 137.1 2969.0 14 0.000 8.14 0.58 0.29 939.89 140.6 198.0 0
      echo "   $l  $s  $offset           $dUdLall  $id" >> $getdUdLs/dUdL_all        ### das ist jetzt ok
      ## $l = lambda = 0.0
      ## dUdL = 0.29 8.11 4.05 936.0 137.1 2969 14.850
      ## offset = 30
      echo "    $l    $dUdL    $offset" | tee -a $getdUdLs/dUdL
    else
      echo "#   $l  seed: $s  too few steps: all `awk 'END{print NR-1}' $file`  off: `cat $f/getAvgAndDev_input_fre | wc -l | sed 's|[ ]*||g'`  unc: $nUncor" | tee -a $getdUdLs/dUdL_all
    fi
    #echo jojojojojojo
    #echo $l $dUdL $offset
    #exit
done
cp $getdUdLs/dUdL_all dUdL_all






#############################################################################################"
## create avg_dUdL
#############################################################################################"
echo
echo "#############################################################################################"
echo "## create avg_dUdL"
echo "#############################################################################################"
for end in $list;do ## list="fre ene eS0"
    [ -e "$getdUdLs/avg_dUdL_$end" ] && mv $getdUdLs/avg_dUdL_$end $getdUdLs/avg_dUdL_$end\_old; echo;echo;
    echo "# correlation length $corSteps" > $getdUdLs/avg_dUdL_$end; 
    echo "# all energies per atom" >>  $getdUdLs/avg_dUdL_$end
    echo "# lambda  <dUdL>(meV) errUnc/2(meV)     stdDev(meV)  <T>(K)  stdDev(K)  steps   time(ps) offset" | tee -a  $getdUdLs/avg_dUdL_$end
    files=`ls $getdUdLs/dUdL_lambda[.0-9]*_all_$end 2> /dev/null` ## schleife ueber lambdas {getdUdLs/dUdL_lambda0.0_all_$end getdUdLs/dUdL_lambda0.2_all_$end ...}
    
    for f in $files; do ## f={getdUdLs/dUdL_lambda1.0_all_{fre,ene,eS0},getdUdLs/dUdL_lambda0.9_all_{fre,ene,eS0}}
    lines=`wc -l $f | awk '{print $1}'`;[ "$lines" = "0" ] && continue   ### necessary not to write 0 in avg_dUdL when nothin was calculated
    #echo f: $f  
      cp $f tmp
      getAvgAndDev
      #echo $f
      cp tmp2 $f\_getAvgAndDev
      nUncor=`grep dUdLUnc tmp2 | awk '{printf("%d",$8)}'`

      l=`echo $f | sed 's/.*dUdL_lambda\(.*\)_all\(.*\)/\1/'` ## lambda: 0.4 oder 1.0 ...
      #dUdL=`grep "dUdL " $f\_getAvgAndDev | awk '{printf("%6.2f %9.2f %19.2f %10.1f %9.1f %8d", $2,$3,$4,$5,$6,$7)}'`  ## so nicht
      dUdL=`grep "dUdL " $f\_getAvgAndDev | awk '{printf("%6.2f %9.2f %19.2f %10.1f %9.1f %8d", $2,$5,$3,$6,$7,$8,$9)}'` 
      time="   000"; [ -e "$getdUdLs/tmp_time_lambda$l" ] && time=`awk 'BEGIN {s=0}; NF==1{s=s+$1}; END{printf("%9.3f",s)}' $getdUdLs/tmp_time_lambda$l`
      echo "    $l    $dUdL $time    $offset" | tee -a $getdUdLs/avg_dUdL_$end
    done

cp $getdUdLs/avg_dUdL_$end avg_dUdL_$end
done
rm -f tmp tmp2 tmp_time_lambda*_offset*
