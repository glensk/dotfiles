#!/bin/sh

out=no #yes #(print additional info for debugging when yes)
path=`set | grep BASH_SOURCE | sed 's|.*\"\(.*\)/[^/]*\".*|\1|' | sed '/BASH_SOURCE=/s/.*/\./'`;[ "$out" = "yes" ] && echo path: $path
script=`set | grep BASH_SOURCE | sed 's|.*\".*/\([^/]*\)\".*|\1|' | sed '/BASH_SOURCE=/s/.*/\./'`;[ "$out" = "yes" ] && echo script: $script
options=$*; . $path/../utilities/functions.include; checkOptions "-h -help -p -i -m -k -o -a -d -c -v";[ "$out" = "yes" ] && echo options: $options

if [ `getOption -h` = True ] || [ `getOption -help` = True ] || [ $# = 0 ]; then
  usage "`basename $0` *Ang*K folder (to rund on all *Ang_*K folders)"
  printOptions " " \
               "-a                to run on all *Ang_*K folders" \
               "-c                to run in current *Ang_*K folder" \
               "-m [REAL]         make/create new jobs if ErrorUncor/2 > INT (e.g.: `basename $0` * -m 0.05)" \
               "-v                be more verbose" \
               " "
  exit
fi


## run script in high folder (not in *Ang*K folder)
## OUTCARS koennen in folgender Ordnerstruktur hinterlegt sien: a) lambda0.0/2345_1200/6x6x6kp_400eV/OUTCAR{.gz} b) lambda0.0/2345_1200/OUTCAR{.gz} 
## altes skript saved in ~/scripts/fah/old.../highUp.sh
## --> braucht KEIN mathematica !!!

## das skript erstellt 
##                     lambda0.0/avg
##                     lambda0.0/dif__400eV_6x6x6kpm0_120NGXF_ED1E-3__low 
##                     lambda0.0/ene__400eV_6x6x6kpm0_120NGXF_ED1E-3 
##                     lambda0.0/ene__low
##                     lambda0.0/std__400eV_6x6x6kpm0_120NGXF_ED1E-3__low 

##                     avg_diff_free__400eV_6x6x6kpm0_120NGXF_ED1E-3__low 
##                     avg_diff_enwe__400eV_6x6x6kpm0_120NGXF_ED1E-3__low 
##                     avg_diff_enw0__400eV_6x6x6kpm0_120NGXF_ED1E-3__low 

###########################################################################################
## file to create: avg_dUdL_high_{fre,ene,eS0}, avg_dUdL_lowplushigh_{fre,ene,eS0} AllEn
## file to create: avg_dUdL_high_{fre,ene,eS0}, avg_dUdL_lowplushigh_{fre,ene,eS0} AllEn
###########################################################################################
vars="fre ene eS0"
std=Fah_up 
avg_dUdL_high=avg_dUdL_high                  # --> immer neu erstellt {avg_dUdL_high_{fre,ene,eS0}
file_AllEn=AllEn                             # --> immer neu erstellt
avg_dUdL_lowplushigh=avg_dUdL_lowplushigh    # --> immer neu erstellt {avg_dUdL_lowplushigh_{fre,ene,eS0}
summary=highUp_summary                       # --> immer neu erstellt
summary_lambdas=highUp_summary_lambdas       # --> immer neu erstellt

######################################################################################################################
## to set by user
######################################################################################################################
all=yes   ## defines if all subfolder are updated (yes -> takes longer, no -> just AllEn is being evaluated)

######################################################################################################################
## check input
######################################################################################################################


hier=`pwd`;
rm -f $hier/tmpout

######################################################################################################################
## 1. get paths (reference_high, lowfolder, ...)
######################################################################################################################

##################################### 
### angK : get angK folder to calculate (and check) angK={3.74Ang_1360K, 3.66Ang_800K, ...}
#####################################
if [ "`getOption -c`" = "True" ];then
       angKs=`echo $hier | sed 's|.*/||'`
       cd ..
        [ ! -e "$angKs" ] && echo "PROBLEM: not found $angKs in $hier" && exit
else 
[ "`find -L . -maxdepth 1 -mindepth 1 -type d -name "*Ang_*K"`" = "" ] && echo 'PROBLEM: no *Ang*K folder found in '`pwd` && exit -1
fi
[ "`getOption -v`" = "True" ] && echo jo: `pwd`
[ "`getOption -v`" = "True" ] && echo ja: $angKs
[ "`getOption -m`" = "True" ] && errorsoll=`getValue -m` && [ "`isnumber.sh $errorsoll`" != "yes" ] && echored "-m expects number but got $errorsoll" && exit
[ "`getOption -a`" != "True" ] && [ "`getOption -c`" != "True" ] && angKs=`echo $* | xargs -n1 | grep "Ang_" | xargs`
[ "`getOption -a`" = "True" ] && angKs=`ls -1d *Ang_*K | sed 's|\([0-9.]*\)Ang_\([0-9]*\)K.*|\2 \1Ang_\2K|' | sort -n | awk '{print $2}'`
[ "`getOption -v`" = "True" ] && echo ja: $angKs
allang=`echo $angKs | xargs -n1 | sed 's|Ang.*||' | sort | uniq`
#####################################
## get low_path && ext low 
##################################### 
ext_low=low
#low_path=`pwd | sed 's|\(.*\)\([__]*\)high_.*|\1|' | sed 's|__$||'`
low_path=`ti_high_0_create_Folders_vasp.sh -l`
## check low path
[ ! -e "$low_path" ] && error "PROBLEM2: low_path: $low_path does not exist" && exit -1
[ "$low_path" = "`pwd`" ] && error "PROBLEM3: low_path = pwd" && exit -1
for i in $angKs;do [ ! -e "$low_path/$i" ] && error "PROBLEM4: low_path/Ang_K: $low_path/$i does not exist" && exit -1; done


##################################### 
## get refPath_high
#####################################
refPath_high=`ti_high_0_create_Folders_vasp.sh -e`
    #refPath_high=`ti.py --get_path_ref_high_job`
    #refPath_high=/data/korbmacher/titanium/anharmonic/hcp/ref_high_3x3x2sc/500eV_06x06x06kp
echo "||| refPath_high: $refPath_high "
[ ! -e "$refPath_high" ] && echo "PROBLEM: refPath_high: $refPath_high does not exist (low_path: $low_path) (refPath_sc: $refPath_sc) (kp1:$kp1 kp2:$kp2 kp3:$kp3)" && exit -1

###########################################################################################
# checks ref_path (check if path exists and if it contains correspondign *Ang folder)
###########################################################################################
for i in $allang;do [ "`ls $refPath_high | grep "$i\Ang" | wc -w | sed 's|[ ]*||g'`" != "1" ] && echo PROBLEM: ref_path: $refPath_high/$i\Ang does not exist && exit -1; done
#[ -e "parameters.dat" ] && [ ! -e "`grep "refPath=" parameters.dat | sed 's|refPath=||'`" ] && echo "refPath=$refPath_high" >> parameters.dat 
[ -e "parameters.dat" ] && [ "`grep "refPath=" parameters.dat | sed 's|refPath=||'`" = "../" ] && sed -i 's|refPath=.*|refPath='"$refPath_high"'|' parameters.dat
[ ! -e "parameters.dat" ] && echo refPath=$refPath_high > parameters.dat

##################################### 
# get high_path und ext_high
#####################################
ext_high=`pwd | sed 's|.*high_||'`
high_path=`pwd`

#echo "######################################################################################################################"
#echo "######################################################################################################################"
#echo "######################################################################################################################"
#echo ""
#echo "$angKs"
#echo ""
#echo "ext_high:         $ext_high"
#echo "low_path:         $low_path"
#echo "refPath_high:     $refPath_high"
#echo ""
#echo "######################################################################################################################"
#echo "######################################################################################################################"
#echo "######################################################################################################################"


######################################################################################################################
## 2. skript
######################################################################################################################

dir=`pwd`
anz_folder=`echo $angKs | wc -w | sed 's|[ ]*||g'`
nr_folder=0
rm -f jobList_additional
for angK in $angKs; do                             ## i=3.74Ang_1300K
  cd $dir
  nr_folder=` expr $nr_folder + 1 `
  [ ! -d "$angK" ] && echo && echo no directory $angK && echo && continue
  cd $angK                   ##                        jetzt bin ich im *Ang_*K folder
  #echo `pwd`
  rm -f $summary $summary_lambdas 
  rm -f avg* A* f* F* lambda*/{asd*,avg*,dif*,Fah*,ion_energ*}
  a=`echo $angK | sed 's|\(.*\)Ang_.*|\1|'`    ## a=3.74
  t=`echo $angK | sed 's|\(.*\)Ang_\(.*\)K.*|\2|'`    ## a=3.74
  all_lambdas=`ls -1d lambda* | sed 's|.*lambda||'`
  echo "######################################################################################################"
  echo "#### in folder ($nr_folder/$anz_folder): $angK    lambdas:" $all_lambdas
  echo "######################################################################################################"

 
  #############################################################################################################  
  ############################################################################################################# 
  # get reference high energies (reference OUTCAR)
  ############################################################################################################# 
  ############################################################################################################# 
               [ ! -e "$refPath_high" ] && echo PROBLEM: refPath_high:$refPath_high does not exist && exit -1
               rph=$refPath_high/$a\Ang
               [ ! -e "$rph" ] && echo PROBLEM: refPath_high/$a\Ang: $refPath_high/$a\Ang does not exist && exit -1
               refOUTCAR=`find -L $refPath_high/$a\Ang -maxdepth 2 -mindepth 1 -name "OUTCAR*"`  # it can also be in workDir
               #echo ref:$refOUTCAR
               if [ ! -e "$refOUTCAR" ]; then
                 echo 1>&2; echo -e "\033[31m\033[1mERROR\033[0m: PROBLEM: reference OUTCAR missing $refOUTCAR" 1>&2; exit
               fi
               nions=`zgrep -a --text "NIONS" $refOUTCAR | awk '{print $12}'`
               nionsm1=` expr $nions - 1 `
               ref_free=`zgrep -a --text "free  en" $refOUTCAR | tail -1 | awk '{print $5}'`  ## wichtung 31/32 oder 30/31 spaeter
               ref_we=`zgrep -a --text "energy  w" $refOUTCAR | tail -1 | awk '{print $4}'`   ## wichtung 31/32 oder 30/31 spaeter
               ref_we0=`zgrep -a --text "energy  w" $refOUTCAR | tail -1 |awk '{print $7}'`  ## wichtung 31/32 oder 30/31 spaeter
			   if [ "$ref_free" = "" ];then
					   #echo in1 ref:$refOUTCAR
					   conv1=`checkConvergence.sh -i $refOUTCAR` 
					   #echo conv1:"$conv1"
					   conv2=`echo "$conv1" | grep -v "^#" | grep -v '^$' | awk '{print $7}'`
					   #echo conv2:$conv2
					   if [ "$conv2" = "conv" -o "$conv2" = "ok" ];then
							   #echo in2
							   ref_free=`OUTCAR_ene-free-last.sh $refOUTCAR`
							   ref_we=`OUTCAR_ene-inner-last.sh $refOUTCAR`
							   ref_we0=`OUTCAR_ene-sigma0-last.sh $refOUTCAR`
					   fi

			   fi
               [ "$ref_free" = "" ] && echored "WARNING: could not get reference energy from $refOUTCAR; nions:$nions:$nionsm1:___ene:$ref_free:$ref_we:$ref_we0:either is is not finished or not calculated" && continue
               [ "`isnumber.sh $ref_free`" != "yes" ] && echored "WARNING: reference energy from $refOUTCAR is not a number but:$ref_free: either run is is not finished or not calculated" && continue
                [ "`getOption -v`" = "True" ] && echo ref_free: $ref_free 
  #echo ref_free: $ref_free 
  
  #############################################################################################################  
  ############################################################################################################# 
  # get high OUTCARS path
  ############################################################################################################# 
  ############################################################################################################# 
    allene=""
    allerror123=""
    maxdiflamall=""
    minfinOUTCARS_anz=100000000
    for lambda in $all_lambdas;do                     ### Schleife ueber lambdas
            #echo lambda:$lambda:
            file_energies_high=lambda$lambda/ion_energies           #immer neu erstellt
            file_energies_low=lambda$lambda/ion_energies__$ext_low; #immer neu erstellt
            file_diff=lambda$lambda/dif                             #immer neu erstellt
            file_avg=lambda$lambda/avg                              #immer neu erstellt
            file_std=lambda$lambda/$std;                            #immer neu erstellt
            rm -f $file_energies_high $file_energies_low $file_diff $file_avg $file_std
            touch $file_energies_high $file_energies_low $file_diff $file_avg $file_std
            head=head
            rm -f $file_avg$head
            echo "#n   avg{fre,we,we0}     ? seed     StdDev{fre,we,we0}    ERROR{fre,we,we0}   lambda" > $file_avg$head
            ## get all OUTCARS
            ## get all OUTCARS
            #echo vor find `pwd` lambda:$lambda:
            allOUTCARS=`find -L lambda$lambda -maxdepth 2 -mindepth 1 -name "OUTCAR*"` # this is slower but will also find -L OUTCARS in workDir
            #echo nnn find:$allOUTCARS 
		allOUTCARS_tmp=""
		for iii in $allOUTCARS;do
                        #echo iii:$iii
			pfad=`echo $iii | sed 's|/OUTCAR.*||g'`
			#echo pfad:$pfad
                        #echo allOUTCARS_tmp:$allOUTCARS_tmp
                        checkgrep=`echo $allOUTCARS_tmp | grep $pfad | wc -w`
                        #echo checkgrep:$checkgrep
			[ "$checkgrep" = "0" ] && allOUTCARS_tmp="$allOUTCARS_tmp $iii"
		done
		
	allOUTCARS=$allOUTCARS_tmp
            allOUTCARS_anz=`echo $allOUTCARS | wc -w | sed 's|[ ]*||g'`
            finOUTCARS_anz=0

  
            ## check if all OUTCARS are finished
            ## check if all OUTCARS are finished
            ## check if all OUTCARS are finished
            notfin=0;finOUTCARS=""
            count=0
            #echo fin: $finOUTCARS            
			#echo allOUTCARS:
			#echo $allOUTCARS
            for i in $allOUTCARS;do
		#echo $i
            	count=` expr $count + 1 `


		conv1c=`checkConvergence.sh -i $i` 
		#echo conv1c:"$conv1c"
		conv2c=`echo "$conv1c" | grep -v "^#" | grep -v '^$' | awk '{print $7}'`
		#echo conv2c:$conv2c
            	## if non packed OUTCAR
            	## if non packed OUTCAR

	        [ "$conv2" = "conv" -o "$conv2c" = "ok" ] && finOUTCARS="$finOUTCARS $i" && finOUTCARS_anz=` expr $finOUTCARS_anz + 1 ` && continue
				
            	notfin=` expr $notfin + 1 `
            	#echo "`pwd`/$i" >> ~/VASPnotfin
            done
            [ "$finOUTCARS_anz" -le "$minfinOUTCARS_anz" ] && minfinOUTCARS_anz=$finOUTCARS_anz



            ll_anz=`echo $finOUTCARS | wc -w | sed 's|[ ]*||g'`

            [ "$finOUTCARS" = "" ] && echo "   l=$lambda   --> problem: no finished OUTCARS in `pwd`" && continue ##(continue geht in die naechste iteration der lambdaschleife)
  
            number=0  ## um cooler plotten zu koennen
            #echo fin: $finOUTCARS
            
            
            ## XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
            ## OUTCAR Schleife
            ## XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
            for ii in $finOUTCARS; do  ## ii runs over every OUTCAR lambda0.0/12355_2000/OUTCAR lambda0.8/12356d_1896/OUTCAR
                    #continue
                    #echo ii: $ii
                    seed=`echo $ii | sed 's|lambda\([^/]*\)/\([^/]*\)_\([^/]*\)/.*|\2|'`        #12355 14675d solll das d hier mit rein?
                    structure=`echo $ii | sed 's|lambda\([^/]*\)/\([^/]*\)_\([^/]*\)/.*|\3|'`   #2000
  
                    #############################################################################################################  
                    ############################################################################################################# 
                    # get high energies (OUTCAR)
                    ############################################################################################################# 
                    ############################################################################################################# 
                    [ "`getOption -v`" = "True" ] && echo ref_free:$ref_free: nionsm1:$nionsm1: ii:$ii:: 
#zgrep -a --text "energy  w" $ii | tail -1                    
#echo ref_free:$ref_free
                    #echo nionsm1:$nionsm1
                    high_ene_free=`zgrep -a --text "free  en" $ii | tail -1`
                [ "`getOption -v`" = "True" ] && echo ii           :::$ii
                [ "`getOption -v`" = "True" ] && echo high_ene_free:::$high_ene_free:::
                [ "`getOption -v`" = "True" ] && echo ref_free     :::$ref_free:::
                [ "`getOption -v`" = "True" ] && echo nionsm1      :::$nionsm1:::
                [ "`getOption -v`" = "True" ] && echo "1000*(high_ene_free-ref_free)/nionsm1))"
                    high_ene_free=`zgrep -a --text "free  en" $ii | tail -1 | awk '{printf("%.2f",1000*($5-('$ref_free'))/'$nionsm1')}'`
                    #echo AAAAAAA: $high_ene_free
                [ "`getOption -v`" = "True" ] && echo high_ene_free:::$high_ene_free:::
                    high_ene_we=`zgrep -a --text "energy  w" $ii | tail -1 | awk '{printf("%.2f",1000*($4-('$ref_we'))/'$nionsm1')}'`
                    high_ene_we0=`zgrep -a --text "energy  w" $ii | tail -1 | awk '{printf("%.2f",1000*($7-('$ref_we0'))/'$nionsm1')}'`
                   # echo high_ene_free:$high_ene_free 
                    [ "$high_ene_free" = "" ] && high_ene_free=`zgrep -a --text "free en" $ii | tail -1 | awk '{printf("%.2f",1000*($5-('$ref_free'))/'$nionsm1')}'`
                    [ "$high_ene_we" = "" ] && high_ene_we=`zgrep -a --text "energy w" $ii | tail -1 | awk '{printf("%.2f",1000*($4-('$ref_we'))/'$nionsm1')}'`
                    [ "$high_ene_we0" = "" ] && high_ene_we0=`zgrep -a --text "energy w" $ii | tail -1 | awk '{printf("%.2f",1000*($7-('$ref_we0'))/'$nionsm1')}'`
          #echo "++>>>>> `zgrep -a --text "free en" $ii | tail -1 `"
          #echo "-->>>>> ii: $ii high: $high_ene_free ref: $ref_free nions: $nionsm1"
                [ "$high_ene_free" = "" ] && continue    
                [ "`getOption -v`" = "True" ] && echo high_ene_free:::$high_ene_free:::
                    #############################################################################################################  
                    ############################################################################################################# 
                    # get low energies (OUTCAR) (from ion_energies, this is the file which was created either form vasprun or from OUTCAR)
                    ############################################################################################################# 
                    #############################################################################################################  
                    low_enes=""
  # auskommen      tieren wenn es schnell gehen soll
  #                 ## 1. Moegleichkeit: find -Let energie im high folder (in $file_energies_low)
  #                 ## 2. Moegleichkeit: geht direkt zum low folder
  #                 [ -e "$file_energies_low" ] && low_enes=`awk '$5=='$seed'&&$6=='$structure'{print $2,$3,$4}' $file_energies_low | head -1`
  #                 # echo "low_enes:$low_enes: <-- file_energies_low: $file_energies_low"
  
                    ## wenn er hier ist hat er hier keine low energies gefunden -> 2te. moeglichkeit
                    #echo "lambda:$lambda seed:$seed ii:$ii"
                    
                    #if [ "`echo $low_enes | wc -w | sed 's|[ ]*||g'`" != "3" ];then
                      #ion_energies_low=`ls -1d $low_path/$angK/lambda$lambda\_$seed*/ion_energies`  #kein sternchen vor lambda, sonst koennte er was anderes find -Len
                      #echo ::: $low_path/$angK/lambda$lambda\_$seed
                      checkpathionene=$low_path/$angK/lambda$lambda\_$seed
                      #echo lll:$checkpathionene
                      ion_energies_low=`find -L $checkpathionene -maxdepth 2 -name ion_energies`  #kein sternchen vor lambda, sonst koennte er was anderes find -Len
                      #echo iii: $ion_energies_low
                      #exit
                      seedgrep=`echo $seed | grep -o "[0-9]*"`
                    [ "`getOption -v`" = "True" ] && echo "low_enes     > $ion_energies_low"
                    [ "`getOption -v`" = "True" ] && echo "seedgrep     > $seedgrep"
                    [ "`getOption -v`" = "True" ] && echo "structure    > $structure"
                      [ -e "$ion_energies_low" ] && low_enes=`awk '$5=='$seedgrep'&&$6=='$structure'{print $2,$3,$4}' $ion_energies_low | head -1`
                      
                      #[ ! -e "$ion_energies_low" ] && ion_energies_low=`ls -1d $low_path/$angK/'*'lambda$lambda\_$seed*/ion_energies`
                    [ "`getOption -v`" = "True" ] && echo "seedgrep     > $seedgrep"
                    [ "`getOption -v`" = "True" ] && echo "seed         > $seed"

                      [ ! -e "$ion_energies_low" ] && echo PROBLEM: ion_energies_low does not exist: $ion_energies_low && continue
                #echo -e "\nlambda:$lambda seed:$seed low_enes:$low_enes: <-- ion_energies_low: $ion_energies_low"
                #echo -e "\nlambda:$lambda seed:$seed high_ene:$high_ene_free $high_ene_we $high_ene_we0: <-- ion_energies_hih: $ii nionsm1: $nionsm1"
                    #fi
                    #[ "$low_enes" != "" ] && [ "`echo $low_enes | wc -w | sed 's|[ ]*||g'`" != "3" ] && echo low enes in $file_energies_low are not 3 numbers $low_enes && exit -1
                    ###################
                    #echo le: $low_enes ion_energies_low: $ion_energies_low fel: $file_energies_low seed: $seed structure: $structure ext_low: $ext_low
                    ###################

  #                 # low: ion_energies_$ext does not exist: && ref low is not low folders of TDI
                    [ "$low_enes" = "" ] && [ "$ext_low" != "low" ] && echo PROBLEM: low reference ion_energies not found: $file_energies_low && continue  ## continue sollte auch gehen statt exit, oder
  #                 # low: check
                    [ "`echo $low_enes | wc -w | sed 's|[ ]*||g'`" != "3" ] && echo "problem:(try not to cut away last step in ion_energies) low_enes: $low_enes" try: awk '$5=='$seedgrep'&&$6=='$structure'{print $2,$3,$4}' $ion_energies_low | head -1  && continue
                    [ "`echo $low_enes | wc -w | sed 's|[ ]*||g'`" != "3" ] && echo PROBLEM: low enes are not 3 numbers:"$low_enes": from ion_energies_low: $ion_energies_low : && continue  ## continue sollte auch gehen statt exit, oder
                    
  #                 # diff 
                    diff=`echo $low_enes | awk '{print '"$high_ene_free"'-$1,'"$high_ene_we"'-$2,'"$high_ene_we0"'-$3}'`
                    #echo low: $low_enes
                    #echo high: $high_ene_free
                    #echo diff: $diff
                    #echo low_enes: $low_enes 
                    #echo high_ene_free:$high_ene_free 
                    #echo high_ene_we:$high_ene_we 
                    #echo high_ene_we0:$high_ene_we0
                    #echo diff: $diff
                    
                [ "`getOption -v`" = "True" ] && echo "low_enes     > $low_enes"
                [ "`getOption -v`" = "True" ] && echo "high_ene_free> $high_ene_free <"
                [ "`getOption -v`" = "True" ] && echo "high_ene_diff: $diff" && echo
                #[ "`echo $diff | awk '{print $1}'`" = "17.77" ] && exit
                #[ "$l" = "0.85" ] && exit
                
                    ## | awk '{printf "%.0f   %.2f %.2f %.2f    %.0f %.0f \n", $1,$2,$3,$4,$5,$6}'
                    echo $number $low_enes $seed $structure | awk '{printf "%.0f   %.2f %.2f %.2f    %.0f %.0f \n", $1,$2,$3,$4,$5,$6}' >> $file_energies_low
                    rm -f tmp
                    cat $file_energies_low | sort -n | uniq > tmp; mv tmp $file_energies_low 
                    echo $number $high_ene_free $high_ene_we $high_ene_we0 $seed $structure | awk '{printf "%.0f   %.2f %.2f %.2f    %.0f %.0f \n", $1,$2,$3,$4,$5,$6}' >> $file_energies_high
                    #echo $number $diff $seed $structure | awk '{printf "%.0f   %.2f %.2f %.2f    %.0f %.0f %.1f\n", $1,$2,$3,$4,$5,$6,'"$lambda"'}' 
                    echo $number $diff $seed $structure | awk '{printf "%.0f   %.2f %.2f %.2f    %.0f %.0f %.1f\n", $1,$2,$3,$4,$5,$6,'"$lambda"'}' >> $file_diff

                    cat $file_diff | awk '{d1=$2-avg1;d2=$3-avg2;d3=$4-avg3;avg1+=d1/NR;avg2+=d2/NR;avg3+=d3/NR;m1+=d1*($2-avg1);m2+=d2*($3-avg2);m3+=d3*($4-avg3)}\
                    {printf "%.0f   %.2f %.2f %.2f    %.0f %.0f      %.2f %.2f %.2f       %.2f %.2f %.2f  %-4s\n",\
                    $1,avg1,avg2,avg3,$4,$5,sqrt(m1/NR),sqrt(m2/NR),sqrt(m3/NR),sqrt(m1/NR)/sqrt(NR),sqrt(m2/NR)/sqrt(NR),sqrt(m3/NR)/sqrt(NR),'"$lambda"'}' > $file_avg

                    #echo $number $diff $seed $structure
                    
                    ##################################################################################################################################
                    #echo -en "\rl=$lambda ($number ($finOUTCARS_anz/$allOUTCARS_anz) OUTCARS) $ii"
                    #echo -en "\r$ii"
                    ##################################################################################################################################
          number=`expr $number + 1 `
          ##################################################################################################################################
          #echo -en "\rl=$lambda ($number/$allOUTCARS_anz)"
          echo -en "\r$number/$allOUTCARS_anz"
          ##################################################################################################################################


  done   ## runs over all outcar
  
 ## if just one value in avg file, delete it, std error only meaningfull if we have 2 values at least
 for i in `find lambda* -name avg`;do
     #echo ">> $i"
     linecheckavg=`cat $i | wc -l`
     #echo lc:$linecheckavg:
     [ "$linecheckavg" = "1" ] && rm -f $i && touch $i

 done
          ## XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
          ## OUTCAR Schleife finished
          ## XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
            
          maxdiflam=`awk '{print $2}' $file_diff | awk '{if(min==""){min=max=$1}; if($1>max) {max=$1}; if($1<min) {min=$1}; total+=$1; count+=1} END {printf "%.2f",  max-min}'`
          ##################################################################################################################################
          #echo -en "\r($number / $allOUTCARS_anz)"
          ##################################################################################################################################

          ##################################################################################################################################
          file_avg_print=`awk '{print $10}' $file_avg | tail -1`
          fre_print=`awk '{print $2}' $file_avg | tail -1`
          echo -e "\rl=$lambda ($number / $ll_anz / $allOUTCARS_anz) fre:$fre_print  +/- $file_avg_print  ( $maxdiflam )" | tee -a $summary_lambdas




          ### decide if new jobs to create
          if [ "`getOption -m`" = "True" ];then
              error_is=$file_avg_print
              error_soll=`getValue -m`
              jobs_is=$number 
              faktor=`echo $error_is $error_soll | awk '{printf "%.1f",($1/$2)^2}'`  # Standard Error = sigma / Sqrt(steps)
              jobs_soll=`echo $jobs_is $faktor | awk '{print $1*$2}'`
              newjobs=`echo $jobs_is $jobs_soll | awk '{printf "%.f", $2-$1}'`

              # if newjobs < 0 : continue ... no new jobs necessary
              check=`echo $newjobs | awk '{if ($1<=0) print "done";else print "calc"}'`
              if [ "$check" != "done" ] ;then
                
                # newjobs should be maximal amount of jobs_is, not more (we could change this later)
                newjobs=`echo $newjobs $jobs_is | xargs -n1 | sort -n | head -1`  ## take the smaller number
                if [ "`getOption -v`" = "True" ];then
                echo error_is: $error_is
                echo error_soll: $error_soll
                echo faktor:$faktor
                echo jobs_is:$jobs_is
                echo jobs_soll:$jobs_soll
                echo newjobs:$newjobs
                echo angK: $angK
                echo alat: $a
                echo temp: $t
                echo lambda: $lambda
                echo pwd1: `pwd`
                fi
                # change parameters.dat and submitt stuff
                vor_goto_high=`pwd`
                cd ..
                #echo pwd2:`pwd`
                ti_high_0_create_Folders_vasp.sh -ua $a -ut $t -ul $lambda -un $newjobs -c
                touch jobList_additional
                cat jobList >> jobList_additional
                cd $vor_goto_high
              fi
          fi
          #echo ""
          ##################################################################################################################################
          #allene="$allene `echo $mean123 | awk '{print $1}'`"
          #allerror123="$allerror123 `echo $error123 | awk '{print $1}'`"
          maxdiflamall="$maxdiflamall $maxdiflam"
          
  done    ## runs over all lambdas
  ## XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
  ## lambda Schleife finished
  ## XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX  

  ## maximale differenz zwichen den einzeln berechneten high werten 
  enediff=`ls -1d lambda*/avg | sed 's|\(.*\)|tail -1 \1|' | sh | awk '{print $2}' | xargs -n1 | awk 'min=="";max==""|| $1<min {min=$1};$1>max {max=$1};{ print $1,max-min}' | awk 'END{print $2}'`
  maxerror=`ls -1d lambda*/avg | sed 's|\(.*\)|tail -1 \1|' | sh | awk '{print $10}' | xargs -n1 | awk 'max==""|| $1>max {max=$1}END{ print max}'`
  mm=`echo $maxdiflamall | xargs -n1 | awk '{if(min==""){min=max=$1}; if($1>max) {max=$1}; if($1<min) {min=$1}; total+=$1; count+=1} END {printf "%.2f", max}'`
  echo -e "`pwd` || max diff: fre \033[31m\033[1m $enediff \033[0m +/-\033[31m\033[1m $maxerror \033[0m (\033[31m\033[1m $mm \033[0m) (min $minfinOUTCARS_anz OUTCARS)" | tee -a $summary

  echo "" 
  ##################################################################################################################################
  # create AllEn
  ################################################################################################################################## 
  rm -f $file_AllEn $avg_dUdL_high* $avg_dUdL_lowplushigh*
  cat lambda*/dif | awk '{print $7,$2,$3,$4,$5,$6}' >> $file_AllEn; [ "`cat $file_AllEn`" = "" ] && rm -f $file_AllEn && echored "PROBLEM:                 nothing calculated in `pwd`" && cd .. && continue
  
  for var in $vars;do
        echo "### create avg_dUdL_high_$var"
        ##################################################################################################################################
        # create avg_dUdL_high_{fre,ene,eS0}
        ################################################################################################################################## 
        echo "#" | tee -a $avg_dUdL_high\_$var > /dev/null
        echo "# all energies per atom #" | tee -a $avg_dUdL_high\_$var > /dev/null
        echo "# lambda avg(meV/at) err/2(meV/at)    stdDev(meV/at) err(meV/at)" | tee -a $avg_dUdL_high\_$var > /dev/null
        [ "$var" = "fre" ] && ls -1d lambda*/avg | sed 's|\(.*\)|tail -1 \1|' | sh | awk '{printf "%-4s      %.2f        %.2f             %.2f           %.2f\n", $13,$2,$10/2,$7,$10}' >> $avg_dUdL_high\_$var
        [ "$var" = "ene" ] && ls -1d lambda*/avg | sed 's|\(.*\)|tail -1 \1|' | sh | awk '{printf "%-4s      %.2f        %.2f             %.2f           %.2f\n", $13,$3,$11/2,$8,$11}' >> $avg_dUdL_high\_$var
        [ "$var" = "eS0" ] && ls -1d lambda*/avg | sed 's|\(.*\)|tail -1 \1|' | sh | awk '{printf "%-4s      %.2f        %.2f             %.2f           %.2f\n", $13,$4,$12/2,$9,$12}' >> $avg_dUdL_high\_$var
        [ "`cat $avg_dUdL_high\_$var | wc -l | sed 's|[ ]*||g'`" = "3" ] && echo PROBLEM: $avg_dUdL_high\_$var not created in `pwd` && rm -f $avg_dUdL_high\_$var && break ## wenn nur 3 zeilen hat nix reingeschrieben
        
        ##################################################################################################################################
        # get avg_dUdL_low_{fre,ene,eS0}
        # get Fah (low)
        ################################################################################################################################## 
        avg_dUdL_low=`pwd | sed 's|__high.*/|/|' | sed 's|\(.*\)|\1/avg_dUdL_'"$var"'|'`
        [ -e "$avg_dUdL_low" ] && cp $avg_dUdL_low avg_dUdL_low_$var
        Fah_low=`pwd | sed 's|__high.*/|/|' | sed 's|\(.*\)|\1/Fah|'`
        [ -e "$Fah_low" ] && cp $Fah_low Fah_low

        [ ! -e "$avg_dUdL_low" ] && echo "PROBLEM (with avg_dUdL_low): $avg_dUdL_low does not EXIST" && continue
        [ ! -e "$Fah_low" ] && echo "PROBLEM (with Fah low): $Fah_low does not EXIST"

        ##################################################################################################################################
        # create avg_dUdL_lowplushigh_{fre,ene,eS0}
        #echo "# create avg_dUdL_lowplushigh_{fre,ene,eS0}"
        #echo "lowplushigh:      $avg_dUdL_lowplushigh\_$var"
        #echo "low:              avg_dUdL_low_$var"
        #echo "high:             $avg_dUdL_high\_$var"
        #echo "###########"
        #[ -e "avg_dUdL_low_$var" ] && echo doesEXIST
        #rm -f ka_$var
        #paste avg_dUdL_low_$var $avg_dUdL_high\_$var | tee -a ka_$var
        
        ################################################################################################################################## 
        ## avg_dUdL_lowplushigh_$var    
        ################################################################################################################################## 
        echo "#"  | tee -a $avg_dUdL_lowplushigh\_$var > /dev/null
        echo "# all energies per atom #" | tee -a $avg_dUdL_lowplushigh\_$var > /dev/null
        echo "# lambda avg(meV/at) err/2(meV/at)    stdDev(meV/at) err(meV/at)" | tee -a $avg_dUdL_lowplushigh\_$var > /dev/null
        lllamd=`awk '{print $1}' avg_dUdL_low_$var | tail -n+4`
        [ "`paste avg_dUdL_low_$var $avg_dUdL_high\_$var | tail -n+4 | awk '(NF==14){print $0}' | wc -l | sed 's|[ ]*||g'`" = "0" ] && printerror=yes && rm -f $avg_dUdL_lowplushigh\_$var && echo prob 33 && continue
        for lll in $lllamd;do
            #echo lll:$lll
            check44=`awk '$1=='"$lll"'{print $0}' avg_dUdL_low_$var | wc -l | sed 's|[ ]*||g'`
            check55=`awk '$1=='"$lll"'{print $0}' $avg_dUdL_high\_$var | wc -l | sed 's|[ ]*||g'`
            [ "$check44" = "0" ] && echo prob 44:$check44:lll:$lll && continue
            [ "$check55" = "0" ] && echo "for lambda $lll not uptilds run! | check55:$check55: | `pwd` $lll $avg_dUdL_high\_$var" && continue
        awk '$1=='"$lll"'{print $0}' avg_dUdL_low_$var $avg_dUdL_high\_$var | awk '{sum1+=$1 ; sum2+=$2 ; sum3+=$3; sum4+=$4} {printf "%.2f %.2f %.2f %2f\n", $1,sum2,sum3,sum4}' | tail -1 >> $avg_dUdL_lowplushigh\_$var

        done
        #echo done
        #[ "`paste avg_dUdL_low_$var $avg_dUdL_high\_$var | tail -n+4 | awk '(NF==14){print $0}' | wc -l | sed 's|[ ]*||g'`" = "0" ] && printerror=yes && rm -f $avg_dUdL_lowplushigh\_$var && continue
        #   paste avg_dUdL_low_$var $avg_dUdL_high\_$var | tail -n+4 | awk '(NF==14){print $0}' | awk '{printf "%-5s  %-7s   %-6s   %-5s \n", $1,$2+$11,$3+$12,$4+$13}' >> $avg_dUdL_lowplushigh\_$var
  done  ## dieses done ist fuer var: fre, ene, eS0

  cd ..
  ##################################################################################################################################
  echo -en ""
  ################################################################################################################################## 
done            ## laeuft ueber alle *Ang_*K folder angK=3.74Ang_1300K)


echo
echo
#cat $hier/tmpout | tee -a $hier/highUp_summary
rm -f $hier/tmpout
