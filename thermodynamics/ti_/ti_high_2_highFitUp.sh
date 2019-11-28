#!/bin/sh
        #################################################################################################
        ## create Fah_high 
        #################################################################################################
        #################################################################################################
        ## create Fah_lowplushigh
        #################################################################################################
        #################################################################################################
        ## create Fah_lowplushighmean
        #################################################################################################

echo "######################################################"
echo "creating AvgEn..."
echo "######################################################"
out=no #yes #(print additional info for debugging when yes)
path=`set | grep BASH_SOURCE | sed 's|.*\"\(.*\)/[^/]*\".*|\1|' | sed '/BASH_SOURCE=/s/.*/\./'`;[ "$out" = "yes" ] && echo path: $path
script=`set | grep BASH_SOURCE | sed 's|.*\".*/\([^/]*\)\".*|\1|' | sed '/BASH_SOURCE=/s/.*/\./'`;[ "$out" = "yes" ] && echo script: $script
options=$*; . $path/../utilities/functions.include; checkOptions "-h -help -p -i -k -o -a -d -c -v";[ "$out" = "yes" ] && echo options: $optionsA

if [ `getOption -h` = True ] || [ `getOption -help` = True ] || [ $# = 0 ]; then
  usage "`basename $0` *Ang*K folder (to rund on all *Ang_*K folders)"
  printOptions " " \
               "-a                (to run on all *Ang_*K folders)" \
               "-c                (to run in current *Ang_*K folder)" \
               "-v                be more verbose" \
               " "
  exit
fi

#echo math:$math
checkAndSetMath

########################################################################
## get folder
########################################################################
l="$*"
[ "`getOption -a`" = "True" ] && l=`ls -1d *Ang_*K | sed 's|\([0-9.]*\)Ang_\([0-9]*\)K.*|\2 \1Ang_\2K|' | sort -n | awk '{print $2}'`
if [ "`getOption -c`" = "True" ];then
       l=`echo $hier | sed 's|.*/||'`
       cd ..
       [ ! -e "parameters.dat" ] && echo "!!!!!!!!!!!! parameters.dat does not exist !!!!!!!!!!!!!!!!"
       #[ -e "parameters.dat" ] && echo parameters.dat does exist
        [ ! -e "$l" ] && echo "PROBLEM: not found $l in $hier" && exit
fi







echo "IF THIS IS NOT WORKING: CHECK LICENCE!"
contributions="eS0 ene fre"
########################################################################
## get low folder with avg_dUdL
########################################################################
#low_folder=`pwd | sed 's|__high.*||'`
low_folder=`ti_high_0_create_Folders_vasp.sh -l`

########################################################################
## check necessary scripts
########################################################################
#fit=fah_low_fit_avg_dUdL.sh
fit=ti_low_5_fit_avg_dUdL.sh
[ ! -e "`which $fit`" ] && echo PROBLEM: $fit command not available && exit -1


########################################################################
## mathematica module linearAndQuadFit 
########################################################################
linearAndQuadFit() {
$math >>/dev/null << EOF 
  avg=Import["tmp","Table"];
  lin=Fit[avg,{1,x},x];
  qua=Fit[avg,{1,x,x^2},x];
  {linfit,quafit}=Append[Table[{i,#/.x->i},{i,0,1,0.01}],{Null,Null}]&/@{lin,qua};
  Export["fit_high_linear",linfit,"Table"]
  Export["fit_high_quad",quafit,"Table"]
  {lindev,quadev}=Sum[Abs[(#/.x->avg[[i,1]])-avg[[i,2]]],{i,Length[avg]}]/Length[avg]&/@{lin,qua};
  {linerr,quaerr}={lindev,quadev}/Sqrt[Length[avg]];
  {lin,qua}=Integrate[#,{x,0,1}]&/@{lin,qua};
  s=Round[100{lin,lindev}]/100//N;
  sl=Round[100{lin,linerr,lindev}]/100//N;
  Export["AvgEn_linear",ToString[sl[[1]]]<>" "<>ToString[sl[[2]]/2]<>" "<>ToString[sl[[3]]]<>" "<>ToString[sl[[2]]],"String"]
  s=Round[100{qua,quadev}]/100//N;
  sq=Round[100{qua,quaerr,quadev}]/100//N;
  Export["AvgEn_quad",ToString[sq[[1]]]<>" "<>ToString[sq[[2]]/2]<>" "<>ToString[sq[[3]]]<>" "<>ToString[sq[[2]]],"String"]
  If[(avg[[1 ;; -1, 1]] // Union // Length) > 2, sb = sq, sb = sl];
  Export["AvgEn_best",ToString[sb[[1]]]<>" "<>ToString[sb[[2]]/2]<>" "<>ToString[sb[[3]]]<>" "<>ToString[sb[[2]]],"String"]
  x1 = Select[avg[[1 ;; -1, 1]] // Union, # <= 0.5 &] // Max;
  x2 = Select[avg[[1 ;; -1, 1]] // Union, # >= 0.5 &] // Min;
  y1 = (Select[avg[[1 ;; -1]], #[[1]] == x1 &][[1 ;; -1, 2]]) // Mean;
  y2 = (Select[avg[[1 ;; -1]], #[[1]] == x2 &][[1 ;; -1, 2]]) // Mean;
If[x1 == x2,
  {yout = 
    Round[(Select[avg[[1 ;; -1]], #[[1]] == 0.5 &][[1 ;; -1, 2]]) // 
      Mean, 0.01]},
  {yout = 
     Round[((y2 - y1)/(x2 - x1) (0.5 - x1)) + y1, 0.01] // ToString;}
  ];
  (*yout = Round[((y2 - y1)/(x2 - x1) (0.5 - x1)) + y1, 0.01] // ToString;*)
  Export["AvgEn_mitte",ToString[yout]<>" "<>ToString[0.0]<>" "<>ToString[0.0]<>" "<>ToString[0.0],"String"]
EOF
}

dir=`pwd`
#for i in $l;do
#echo $i
#done
anz=`echo $l | wc -w | sed 's|[ ]*||g'`
count=1
for i in $l; do               ## i sind jeder einzelne *Ang_*K folder
  cd $dir
 if [ -d $i ]; then
  echo "####################"; 
  echo "## ($count/$anz) $i"; 
  echo "####################" ;
  cd $i        ## jetzt im *Ang*K folder
  count=` expr $count + 1 `

  rm -f Fah Fah_high Fah_lowplushigh Fah_lowplushighmean fit_* Fah_lowplushighFahbest Fah_lowplushighFahmitte
  rm -f Fah Fah_diff Fah_diff_plus_low Fah Fah_lowplus_Fah_high


  AllEn=`ls -1d AllEn`
  [ "`echo $AllEn | wc -w | sed 's|[ ]*||g'`" != "1" ] && echo "PROBLEM:                    AllEn not found in `pwd`" && continue
  echo " fit      Fah_high(meV/at) err/2(meV/at) dev(meV/at) err(meV/at)" > Fah_high #_$contrib #\__$ext

  for contrib in $contributions;do  ## contrib is: eS0 ene fre
        #################################################################################################
        ## create Fah_high (is always new, removed before)
        #################################################################################################
        rm -f tmp
        [ "$contrib" = "fre" ] && awk '{print $1,$2}' $AllEn > tmp
        [ "$contrib" = "ene" ] && awk '{print $1,$3}' $AllEn > tmp
        [ "$contrib" = "eS0" ] && awk '{print $1,$4}' $AllEn > tmp

        rm -f AvgEn_linear AvgEn_quad AvgEn_best
        
        
        linearAndQuadFit  ## dies ist nur fuer high up, da nur lineares/quadratisches fitting
      
        rm -f tmp
        [ ! -e "AvgEn_linear" ] && echo PROBLEM: AvgEn_linear not created && continue
        [ ! -e "AvgEn_quad" ] && echo PROBLEM: AvgEn_quad not created && continue
        [ ! -e "AvgEn_best" ] && echo PROBLEM: AvgEn_best not created && continue
        error=`echo "$contrib _best       \`cat AvgEn_best\`"   | awk '{print $3}'` 
        error_dudl_max=`cat avg_dUdL_low_fre | grep -v "#" | awk '{print $3}' | sort -n | tail -1`
        echo error: $error :  error_max : $error_max

        echo "$contrib _linear     `cat AvgEn_linear`" | awk '{printf("%-3s%-8s      %.2f        %.2f            %.2f      %.2f\n",$1,$2,$3,$4,$5,$6)}' | tee -a Fah_high #_$contrib #\__$ext
        echo "$contrib _quad       `cat AvgEn_quad`"   | awk '{printf("%-3s%-8s      %.2f        %.2f            %.2f      %.2f\n",$1,$2,$3,$4,$5,$6)}' | tee -a Fah_high #_$contrib #\__$ext
        echo "$contrib _best       `cat AvgEn_best`"   | awk '{printf("%-3s%-8s      %.2f        %.2f            %.2f      %.2f\n",$1,$2,$3,$4,$5,$6)}' | tee -a Fah_high #_$contrib #\__$ext
        echo "$contrib _mitte      `cat AvgEn_mitte`"  | awk '{printf("%-3s%-8s      %.2f        %.2f            %.2f      %.2f\n",$1,$2,$3,$4,$5,$6)}' | tee -a Fah_high
        rm -f AvgEn_linear AvgEn_quad AvgEn_best AvgEn_mitte


        #################################################################################################
        ## create Fah_lowplushigh (is always new, removed before)
        #################################################################################################
        ## avg_dUdL (output file)
        echo check for files
        avg_dUdL=avg_dUdL_lowplushigh_$contrib
        cp avg_dUdL_lowplushigh_$contrib avg_dUdL_$contrib  # necessary for fitting skript
        echo vor $fit skript
        $fit > /dev/null  ## uses fah_low_fit_avg_dUdL.sh (same file is used in fah_low_getdUdLs.sh) -> low_fit_avg_dUdL.txt (
        echo nach $fit skript
        cat Fah >> Fah_lowplushigh
        rm -f fit_cubic_ene fit_cubic_eS0 fit_tangens_eS0 fit_tangens_ene
        [ -e "fit_cubic_fre" ] && mv fit_cubic_fre fit_lowplushigh_cubic_fre
        [ -e "fit_tangens_fre" ] && mv fit_tangens_fre fit_lowplushigh_tangens_fre
        rm -f avg_dUdL_$contrib


        #################################################################################################
        ## create Fah_lowplushighFahbest (is always new, removed before)
        #################################################################################################
        # Unterschied zw. Fah_lowplushigh und Fah_lowplushighFah$ll: 
        # one first adds energies from low+high and integrates thereafter
        # the other does integrate seperately and adds energies from low and high afterwards
        list="best mitte"
        for ll in $list;do
        diff_line=`cat Fah_high | grep "$contrib\_$ll"`
        ene=`echo "$diff_line" | awk '{print $2}'`
        error=`echo "$diff_line" | awk '{print $3}'`
           # echo ene:$ene
           # echo error: $error
        [ -e "Fah_low" ] && cat Fah_low | grep "$contrib" | awk '{printf "%-20s %7.2f %7.2f %7.2f\n",$1,$2+'"$ene"',$3+'"$error"',$4}' >> Fah_lowplushighFah$ll
        done


  done
  rm -f Fah
  cd $dir
 else
  echo; echo "no directory $i"; echo
 fi
 echo 
done

