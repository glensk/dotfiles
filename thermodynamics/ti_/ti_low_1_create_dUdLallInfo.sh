#!/bin/sh

###### check where you are (how deep in folder structure)
gotofolder=`pwd`
checkfolder=`pwd`
gotolam_folder=`find -L $checkfolder -maxdepth 2 -mindepth 1 -type d -name lambda'[0-9.]*_[0-9]*'`  ## keep maxdepth 2 in case you are in *Ang_*K folder
#echo lam:$gotolam_folder
gotolam_anz=`echo $gotolam_folder | wc -w | sed 's|[ ]*||g'`
[ ! -e "OUTCAR" ] && [ ! -e "OUTCAR.gz" ] && [ ! -e "vasprun.xml" ] && [ ! -e "vasprun.xml.gz" ] && [ -d "workDir" ] && gotofolder=`pwd`/workDir
[ ! -e "OUTCAR" ] && [ ! -e "OUTCAR.gz" ] && [ ! -e "vasprun.xml" ] && [ ! -e "vasprun.xml.gz" ] && [ "$gotolam_anz" != "0" ] && gotofolder=$gotolam_folder



####### start createdUdLallInfo.sh
for i in $gotofolder;do ## schleife ueber lambdafolder (oder den lambdafolder in welchem man gerade drin ist
    cd $i 
    [ ! -e "OUTCAR" ] && [ ! -e "OUTCAR.gz" ] && [ ! -e "vasprun.xml" ] && [ ! -e "vasprun.xml.gz" ] && [ -d "workDir" ] && cd workDir

    echo pwd: `pwd`
    pfad=OUTCAR   ## pfad is link to OUTCAR
    [ ! -e "$pfad" ] && pfad=OUTCAR.gz
    [ -e "$pfad" ] && nions=`zgrep -a --text " number of ions     NIONS = " $pfad | awk '{print $12}'`

if [ ! -e "$pfad" ];then
    pfad=vasprun.xml
    [ ! -e "$pfad" ] && pfad=vasprun.xml.gz
    [ -e "$pfad" ] && nions=`zgrep -a --text "<atoms>" $pfad | sed 's|<atoms>||' | sed 's|</atoms>||' | sed 's|^[ \t]*||;s|[ \t]*$||'`
    [ ! -e "$pfad" ] && echo -e "`pwd` did not find neither OUTCAR(.gz) nor vasprun.xml(.gz) ... this run did not start yet" && continue
    fi
# echo nions: $nions
#echo pfad: $pfad

#f1_check=no
#f2_check=no
#f3_check=no
#[ -f "ion_energies" ] && f1_check=ion_energies
#[ -f "ion_energies.gz" ] && f1_check=ion_energies.gz
#[ -f "ion_energies_vasprun.xml" ] && f1_check=ion_energies_vasprun.xml
#[ -f "ion_energies_vasprun.xml.gz" ] && f1_check=ion_energies_vasprun.xml.gz
#[ -f "dUdLallInfo" ] && f2_check=dUdLallInfo
#[ -f "structures" ] && f3_check=structures
#[ -f "structures.gz" ] && f3_check=structures.gz
#[ -f "structures_vasprun.gz" ] && f3_check=structures_vasprun.gz
#echo in `pwd`
#echo f1:$f1_check
#echo f2:$f2_check
#echo f3:$f3_check


# if next line is fulfilled we know the the necessary files exist, but do we know if they are up to datae really?
#[ -f "$f1_check" ] && [ -f "$f2_check" ] && [ -f "$f3_check" ] && [ "$f1_check" -nt "$pfad" ] && [ "$f2_check" -nt "$pfad" ] && [ "$f3_check" -nt "$pfad" ] && echo "`pwd` up to date" && continue
#echo `pwd` danach
#echo f1:$f1_check
#echo f2:$f2_check
#echo f3:$f3_check
#exit



# the removing of the stuff should come only if files are not up to date

#
#      rm -f ion_energies_vasprun.xml ion_energies
#      rm -f dUdLallInfo dUdLallInfo_noJumps   ## wenn mal dUdLallInfo sicherlich so richtig ist und es so bleibt, muss es nicht mehr geloescht werden, sondern wird evtl. von dUdLallInfo ueberschrieben
#      rm -f structures.gz structures_vasprun.gz structures*
#      rm -f getAvgAndDev*
#      rm -f avg_dUdL.dat





## DIESES SKRIPT WIRD DIREKT IM ENTSPRECHENDEM ORDNER (lambdafolder) AUSGEFUEHRT in dem sich das dUdL file befind -Let
## 
## ES ERSTELLT dUdL_allinfo / ion_energies / structures(_vasprun).gz !!!!!!!!!!!!!!!!!
#
## 0) a) check if dUdL exists b) check if dUdLallInfo is already ok
## 1) module to create dUdLallInfo
## 2) create dUdLallInfo from a) vasprun or b) OUTCAR
## 3) create structures from a) vasprun or b) OUTCAR  ( vasprun -> structures_vasprun.xm.lgz       ; OUTCAR -> structures.gz )



#################################################################################################
## convert path to OUTCAR 
#################################################################################################


#################################################################################################
## 0) a) dUdLallInfo should be created by dUdL file and ion_energies file (check if ion energies has 12 rows)
#################################################################################################
[ ! -e "dUdL" ] && echo -e "!!! PROBLEM !!! dUdL does not exist in `pwd` (mayby this job just started)" && continue   

#################################################################################################
## 0) b) maybe dUdLallInof is already ok like it is!? -> if yes, continue
#################################################################################################

      #[ "`tail -1 dUdLallInfo | wc -w | sed 's|[ ]*||g'`" = "15" ] && [ "`tail -1 dUdL | awk '{print $1}'`" = "`tail -1 dUdLallInfo | awk '{print $1}'`" ] && echo `pwd` : dUdLallinfo up to date && continue


#echo 1
#################################################################################################
## 2) module to create dUdLallInfo , needs ion_energies, so make ion_energies first, the thin is that ion_energies can be calculated from a OUTCAR, vasprun, ...
#################################################################################################
createdudlallinfo() {
# dont echo enything otherwishe this module will not work
schritte_ion_energies=`tail -1 ion_energies | awk '{print $1}'`        ## gibt z.b. 2000 bei 2000 schritten
schritte_dUdL=`tail -1 dUdL | awk '{print $1+1}'`
    #echo str:$schrittestr dUdL:$schrittedUdL
[ "$schritte_ion_energies" != "$schritte_dUdL" ] && OUTCAR_create_ion_energies.sh #&& echo kkk
[ ! -e "ion_energies" ] && echo `pwd` : ion_energies couldnt be created && continue
############################
## 1a) check if dUdL is old or newstyle -> if old, rescale U, Uref, dUdL, (average ist egal)
############################
style=new
U=`head -2 ion_energies | tail -1 | awk '{print $2}'`
UdUdL=`head -2 dUdL | tail -1 | awk '{print $5}'`
[ "`echo $U $UdUdL | awk 'sqrt(($1-$2)^2)<=0.02{print $0}'`" = "" ] && style=old
#echo nions: $nions
[ "$style" = "old" ] && nions=$nions && nionsm1=` expr $nions - 1 `
[ "$style" = "old" ] && UdUdLold=`head -2 dUdL | tail -1 | awk '{printf "%.2f\n", $5*'$nions'/'$nionsm1'}'`
   #echo style:$style U:$U UdUdL:$UdUdL nions:$nions nionsm1:$nionsm1 UdUdLold:$UdUdLold
[ "$style" = "old" ] && rm -f dUdLnew && awk '{printf "%7.0f %9.1f %8.1f %8.1f %13.2f %9.2f %13.2f %9.2f %9.2f\n",$1,$2,$3,$4,$5*'"$nions"'/'"$nionsm1"',$6*'"$nions"'/'"$nionsm1"',$7*'"$nions"'/'"$nionsm1"',0,0}' dUdL > dUdLnew
dUdL=dUdL; [ "$style" = "old" ] && dUdL=dUdLnew
###########################
## 2) append energies to dUdLallInfo and check if dUdLenergies in dUdL and dUdLallInfo match
###########################

    # do all this onlly if we have a newly created ion_energies, otherwise keep dUdLallInfo
    match=yes
    doithere="no"
    [ ! -e "dUdLallInfo" ] && doithere="doit since dUdLallInfo does not exist" #&& echo doesnotexist
    if [ -e "dUdL" ];then
       if [ -e "dUdLallInfo" ];then
        c1=`tail -1 dUdLallInfo | awk '{print $1}'`
        c2=`tail -1 dUdL | awk '{print $1}'`
        [ "$c1" != "$c2" ] && doithere="yes since not up to date" 
    fi 
    fi
   
    if [ "$doithere" != "no" ];then
        #echo ck1:$check1:exit
        #exit
        rm -f t0 t1 t2 tt0 tt1 tt2 dUdLallInfo
        awk '{print $1-1,$2}' ion_energies | sed 's|0 0.00|# U(meV/at)|' > t0     ## 2 spalten
        awk '{print $1-1,$3}' ion_energies | sed 's|0 0.00|# Uwe(meV/at)|' > t1   ## 2 spalten
        awk '{print $1-1,$4}' ion_energies | sed 's|0 0.00|# Us0(meV/at)|' > t2   ## 2 spalten
        
        awk 'FNR==NR{a[$1]=$2 $3;next}{ printf "%s %9.2f\n",$0,a[$1]}' t0 $dUdL > tt0
        awk 'FNR==NR{a[$1]=$2 $3;next}{ printf "%s %10.2f\n",$0,a[$1]}' t1 tt0 > tt1
        awk 'FNR==NR{a[$1]=$2 $3;next}{ printf "%s %12.2f\n",$0,a[$1]}' t2 tt1 > tt2
        sed -i '1  s|.*|#  step   time(fs)  temp(K) average       U(meV/at)    Uref          dUdL   average    offset|' tt2
        awk '{printf "%s %10.2f %9.2f %9.2f\n",$0,$10-$6,$11-$6,$12-$6}' tt2 | sed 's|offset.*|offset    U(meV/at)  Uwe(meV/at)  Us0(meV/at)  dUdL    dUdLwe    dUdLs0|' > dUdLallInfo
        rm -f tt0 tt1 tt2 t0 t1 t2 dUdLnew

        # checks if dUdLenergies match
        # checks if dUdLenergies match
        linesdUdL=`tail -n+2 dUdL | wc -l | sed 's|[ ]*||g'`
        linesdUdLallInfo=`tail -n+2 dUdLallInfo | awk 'sqrt(($5-$10)^2)<=0.02{print $0}' | wc -l | sed 's|[ ]*||g'` ## check if energies match!!!
        match=yes; [ "$linesdUdL" != "$linesdUdLallInfo" ] && match=linesdUdL_is_$linesdUdL\_linesdUdLallInfo_is_$linesdUdLallInfo;
        
        ### if macht is no, koennte man checken ob jede iteration vorhanden ist in der OUTCAR, wenn nicht, qdel this job
        if [ "$match" != "yes" ];then
            echo in mach schleife
        OUTCAR_iteration-check.sh OUTCAR
        fi
    fi
    echo $match
}


#################################################################################################
## 1) create ion_energies   
## a) from vasprun (is usually more reliable + has higher accuracy in forces)
## b) from OUTCAR
## c) check if ion_energies was created (if not, -> continue)
#################################################################################################

## a) from vasprun
  #echo 2
modcreateionenergies() {
  [ ! -e "ion_energies" ] && vasprun_create_ion_energies.sh && mv ion_energies_vasprun.xml ion_energies
  [ -e "ion_energies" ] && [ "`head -1 ion_energies | wc -w | sed 's|[ ]*||g'`" != "12" ] && rm -f ion_energies && vasprun_create_ion_energies.sh && mv ion_energies_vasprun.xml ion_energies
  ###########################################################################
  #### change this when fully without OUTCAR
  ## delete last line if not corresponding, this is not necessary when dUdL is full created from vasprun or OUTCAR (&& Hessematrix)
  ###########################################################################
  [ -e "dUdL" ] && [ -e "ion_energies" ] && [ "`tail -1 ion_energies | awk '{print $1}'`" != "`tail -1 dUdL | awk '{print $1+1}'`" ] && sed -i '$d' ion_energies
  #echo lines111a: `tail -1 ion_energies | awk '{print $1}'`
  #echo lines222a: `tail -1 dUdL | awk '{print $1+1}'`
  [ -e "dUdL" ] && [ -e "ion_energies" ] && [ "`tail -1 ion_energies | awk '{print $1}'`" != "`tail -1 dUdL | awk '{print $1+1}'`" ] && rm -f ion_energies &&  vasprun_create_ion_energies.sh && mv ion_energies_vasprun.xml ion_energies
  #echo lines111x: `tail -1 ion_energies | awk '{print $1}'`
  #echo lines222x: `tail -1 dUdL | awk '{print $1+1}'`
  yesno_energies=`createdudlallinfo`  ## if correct energeis -> try structures from OUTCAR otherwise make first energeis from vasprun and structures from vasprun
  #echo JJJ: $yesno_energies
  [ "$yesno_energies" != "yes" ] && rm -f ion_energies

## b) from OUTCAR
  #echo 3
  if [ "$yesno_energies" != "yes" ];then
  echo "small PROBLEM: ion_energies from vasprun did not work out yesno:$yesno_energies"
  rm -f ion_energies
  [ ! -e "ion_energies" ] && OUTCAR_create_ion_energies.sh
  [ -e "ion_energies" ] && [ "`head -1 ion_energies | wc -w | sed 's|[ ]*||g'`" != "12" ] && rm -f ion_energies && OUTCAR_create_ion_energies.sh
  [ -e "dUdL" ] && [ -e "ion_energies" ] && [ "`tail -1 ion_energies | awk '{print $1}'`" != "`tail -1 dUdL | awk '{print $1+1}'`" ] && rm -f ion_energies &&  OUTCAR_create_ion_energies.sh 
  #createdudlallinfo
  yesno_energies=`createdudlallinfo`  ## if correct energeis -> try structures from OUTCAR otherwise make first energeis from vasprun and structures from vasprun
  #echo JJJ: $yesno_energies
  [ "$yesno_energies" != "yes" ] && rm -f ion_energies
  fi
}
## c) check if ion_energies was created (if not, -> continue)
modcreateionenergies
c1=`tail -1 dUdLallInfo | awk '{print $1}'`
c2=`tail -1 dUdL | awk '{print $1}'`
[ "$c1" != "$c2" ] && echo "c1:$c1: c2:$c2: -> make new dUdLallinfo, and structures.gz file" && createdudlallinfo #&> /dev/null
[ ! -e "ion_energies" ] && echo -e "\033[31m\033[1mPROBLEM\033[0m PROBLEM: ion_energies couldnt be created; neither from vasprun nor from OUTCAR;mach:$yesno_energies ( `pwd` )" && continue
#################################################################################################
## 3) create structures, this is the part which takes long, make it quicker
##    thisone needs definitely to be done when ion_energies is longer / shorter than dUdL
##    thisone needs only to be done if we have new ion_energies file or a new dUdLallinfo file -> check if dUdLallinfo or ion_energies is older than 30 seconds
##    dUdLallInfo is created every time, stop that
##    checke NUR ion_energies da dUdLallInfo jedes mal neu erstellt wird
## a) from vasprun
## b) from OUTCAR
#################################################################################################

## from vasprun
doit="no"
#check1=`find dUdLallInfo -mmin -1`
check2=`find ion_energies -mmin -2`  # wenn check2 leer -> ion energies schon lange her erstellt, keine bedarf an 
[ "$check2" != "" ] && doit="ion energies just now created -> new structures necessary"
check3=`find -name "structures*"`
[ "$check3" = "" ] && doit="since no structures fiels found -> create it"  ## if neither structures.gz nor structures_vasprun.gz exists .... than definitely do it

[ "$c1" != "$c2" ] && doit="dUdLallInfo has:$c1: lines and dUdL only :$c2: -> make new structures file"
    if [ "$doit" != "no" ];then 
        #echo doit: "$doit"
        #echo "ion_energies just created"
        structures=structures_vasprun.gz
        outstr=`vasprun_create_structures.sh`
        [ "$outstr" = "" ] && yesno_structures=yes   ## naja das heisst schonmal das das file ueberhaupt evtl. erstellt worden ist :)
        [ ! -e "$structures" ] && yesno_structures=no
        lines_ion_energies=`tail -1 ion_energies | awk '{print $1}'`
        [ -e "$structures" ] && [ "`zgrep -a --text -A0 "^$lines_ion_energies" $structures | tail -1 | awk '{print $1}'`" != "$lines_ion_energies" ] && yesno_structures=no
        [ "$yesno_structures" != "yes" ] && rm -f $structures 
        #echo 6
        [ "$yesno_structures" = "yes" ] && continue
        #echo 7
        
        ## from OUTCAR
        #echo "`pwd` small PROBLEM: structures_vasprun.xml failed"
        structures=structures.gz
        outstr=`OUTCAR_create_structures.sh`
        [ "$outstr" = "" ] && yesno_structures=yes
        [ ! -e "$structures" ] && yesno_structures=no
        [ -e "$structures" ] && [ "`zgrep -a --text -A0 "^$lines_ion_energies" $structures | tail -1 | awk '{print $1}'`" != "$lines_ion_energies" ] && yesno_structures=no
        [ "$yesno_structures" != "yes" ] && rm -f $structures # && echo OUTCAR_create_ion_energies.sh ok but not OUTCAR_create_structures.sh not in `pwd` hopefully vasprun is ok
        #echo 5d
        [ "$yesno_structures" = "yes" ] && continue
        
        
        ## if it comes to this point -> error
        if [ "$yesno_structures" != "yes" ];then
            echo -e "`pwd` \033[31m\033[1m PROBLEM \033[0m PROBLEM OUTCAR/vasprun_create_ion_energies.sh ok but OUTCAR/vasprun_create_structures.sh not"
            rm -f ion_energies
            fi
        #echo 6d
        fi
done
