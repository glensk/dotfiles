#!/bin/bash

out=noyes #(print additional info for debugging when yes)
path=`set | grep BASH_SOURCE | sed 's|.*\"\(.*\)/[^/]*\".*|\1|' | sed '/BASH_SOURCE=/s/.*/\./'`;[ "$out" = "yes" ] && echo path: $path
script=`set | grep BASH_SOURCE | sed 's|.*\".*/\([^/]*\)\".*|\1|' | sed '/BASH_SOURCE=/s/.*/\./'`;[ "$out" = "yes" ] && echo script: $script
options=$*; . $path/../utilities/functions.include; checkOptions "-h -help -v -d -s";[ "$out" = "yes" ] && echo options: $options

if [ `getOption -h` = True ] || [ `getOption -help` = True ] ; then
  usage `basename $0`
  printOptions " " \
               "-d          dont create folder, just make fit(s)" \
               "-s          show fit data and exit" \
               "-v          verbose mode " \ 
  exit
fi
#echo "checkAndSetMath:"
checkAndSetMath
#echo "checkAndSetMath:"`which math`
[ ! -e "`which math`" ] && echo mathematica not set && exit

[ "`pwd | grep -o high`" != "high" ] && auswertung=low  && files="Fah"  
[ "`pwd | grep -o high`" = "high" ]  && auswertung=high && files="Fah_high Fah_lowplushigh Fah_lowplushighFahbest Fah_lowplushighFahmitte"
[ "`pwd | grep -o high`" = "high" ]  && auswertung=high && files="Fah_high Fah_lowplushighFahbest Fah_lowplushighFahmitte"
[ "`pwd | grep -o high`" = "high" ]  && auswertung=high && files="Fah_high Fah_lowplushighFahbest"  # just now to make it quicker

l=`ls -1d [0-9.]*Ang_*K`

createfoldernew="yes"; [ "`getOption -d`" = "True" ] && createfoldernew="no"
hier=`pwd`
#[ "`getOption -s`" = "True" ] && options="$options -v" 
#[ "`getOption -v`" = "True" ] && echo options:$options:

#############################################################################################################################
# get some information for fitting
#############################################################################################################################
if [ "$auswertung" = "high" ];then
    makefit="yes"

    # .... getting data for fit ... 1. hessepath"
    #echo -en "\r.... getting data for fit ... [1/9]. hessepath"
    echo  ".... getting data for fit ... [1/9]. hessepath"
    Hp=`ti_high_0_create_Folders_vasp.sh -Hp`
    [ ! -e "$Hp" ] && makefit="no hessepath" && echored "COULD NOT FIND Hessepath:$Hp: (for fit surface)"
    [ "`getOption -v`" = "True" ] && echo "Hp:$Hp                  "

     
    # .... getting data for fit ... 2. supercell"
    if [ "$makefit" = "yes" ];then
    #echo -en "\r.... getting data for fit ... [2/9]. supercell"
    echo  ".... getting data for fit ... [2/9]. supercell"
    Hs=`ti_high_0_create_Folders_vasp.sh -Hs`
    [ "`checkInteger $Hs`" != "ok" ] && makefit="no Hs" && echored "COULD NOT FIND supercell:$Hs: (for fit surface) in $Hp/parameters.dat"
    [ "`getOption -v`" = "True" ] && echo "Hs:$Hs             "
    fi


    # .... getting data for fit ... 3. cellType"
    if [ "$makefit" = "yes" ];then
    #echo -en "\r.... getting data for fit ... [3/9]. cellType"
    echo  ".... getting data for fit ... [3/9]. cellType"
    Hc=`ti_high_0_create_Folders_vasp.sh -Hc`
    [ "$Hc" != "fcc" ] && [ "$Hc" != "bcc" ] && makefit="no Hc" && echored "COULD NOT FIND structure:$Hc: (fcc/bcc) (for fit surface) in $Hp/parameters.dat"
    structnr="no";[ "$Hc" = "fcc" ] && structnr=4; [ "$Hc" = "bcc" ] && structnr=2; 
    [ "`getOption -v`" = "True" ] && echo "Hc:$Hc               "
    fi


    # .... getting data for fit ... 4. number of atoms"
    if [ "$makefit" = "yes" ];then
    #echo -en "\r.... getting data for fit ... [4/9]. number of atoms"
    echo  ".... getting data for fit ... [4/9]. number of atoms"
    Hn=`ti_high_0_create_Folders_vasp.sh -Hn`
    [ "`checkInteger $Hn`" != "ok" ] && makefit="no Hn" && echo && echored "COULD NOT FIND number of atoms:$Hn: (for fit surface) in $Hp/parameters.dat" && echo
    [ "`getOption -v`" = "True" ] && echo "Hn:$Hn                     "
    fi

    
    # .... getting data for fit ... 5. element"
    if [ "$makefit" = "yes" ];then
    #echo -en "\r.... getting data for fit ... [5/9]. element"
    echo ".... getting data for fit ... [5/9]. element"
    He=`ti_high_0_create_Folders_vasp.sh -He`
    [ "`getOption -v`" = "True" ] && echo "He:$He            "
    fi


    # .... getting data for fit ... 6. bulk or defect"
    if [ "$makefit" = "yes" ];then
    #echo -en "\r.... getting data for fit ... [6/9]. bulk or defect"
    echo ".... getting data for fit ... [6/9]. bulk or defect"
    Hbd=`ti_high_0_create_Folders_vasp.sh -Hbd`
    [ "`getOption -v`" = "True" ] && echo "Hbd:$Hbd                 "
    fi


    # .... getting data for fit ... 7. bulk or Temperature at Melting"
    if [ "$makefit" = "yes" ];then
    #echo -en "\r.... getting data for fit ... [7/9]. bulk or Temperature at Melting"
    echo  ".... getting data for fit ... [7/9]. bulk or Temperature at Melting"
    tm=`getMeltingPoint.sh $He -r`;[ "`isnumber.sh $tm`" != "yes" ] && tm=2000
    [ "`getOption -v`" = "True" ] && echo "tm:$tm                 "
    fi

    [ "`getOption -v`" = "True" ] && echo makefit66:$makefit
    [ "$makefit" = "no" ] && options="$options -v" 


    # .... getting data for fit ... 8. bulk or range alats"
    # .... getting data for fit ... 9 .mean_freqs"
    if [ "$makefit" = "yes" ];then
        #echo -en "\r.... getting data for fit ... [8/9]. bulk or range alats"
        echo ".... getting data for fit ... [8/9]. bulk or range alats"
        rangemin=`ls -1d [0-9.]*Ang_*K | sed 's|Ang.*||' | sort -n | uniq | head -1`
        rangemax=`ls -1d [0-9.]*Ang_*K | sed 's|Ang.*||' | sort -n | uniq | tail -1`
        [ "`checkReal $rangemin`" != "ok" ] && makefit="no" && echored "COULD NOT FIND minimal lattice constant (for fit surface) min: $rangemin"
        [ "`checkReal $rangemax`" != "ok" ] && makefit="no" && echored "COULD NOT FIND maximal lattice constant (for fit surface) max: $rangemax"
        [ "`getOption -v`" = "True" ] && echo rangemin:$rangemin rangemax:$rangemax
        path_mean_freqs=$Hp/mean_freqs
        #echo -en "\r.... getting data for fit ... [9/9]. mean_freqs"
        echo ".... getting data for fit ... [9/9]. mean_freqs"
        if [ ! -e "$path_mean_freqs" ];then
            [ "$Hbd" = "bulk" ] && echo running  getSingleSpeciesPhonons.sh -folder $Hp -fa
            [ "$Hbd" = "bulk" ] && getSingleSpeciesPhonons.sh -folder $Hp -fa
            [ "$Hbd" = "def" ] && echored "CREATE mean_freqs for defect in $Hp"
        fi
        [ ! -e "$path_mean_freqs" ] && makefit="no" && echored "COULD NOT FIND mean_freqs file in $Hp (for fit surface)"
        [ "`getOption -v`" = "True" ] && echo makefit:$makefit
    fi
    [ "`getOption -v`" = "True" ] && echo makefit77:$makefit
    [ "`getOption -v`" = "True" ] && echo Hbd:$Hbd:
    
    makethermodynamics="no"
    ###################################################
    ## for bulk Thermodynamics
    ###################################################
    if [ "$Hbd" = "bulk" ];then  # for bulk we always need the fitperatom auswertung
        echo .... getting data for getThermodynamics.sh ...
        evinet=$Hp/EVinet;[ ! -e "$evinet" ] && echored "WARNING: evinet does not exist in hessepath: $evinet : you might consider to copy it there to get Thermodynamic relations" && evinet=""
        fqhmesh=$Hp/fitperatom/Fqh_fromMeshFreqs_fit_order2;
        if [ ! -e "$fqhmesh" ];then
        [ "`getOption -v`" = "True" ] && echo "getSingleSpeciesPhonons.sh -folder $Hp -makefits"
        getSingleSpeciesPhonons.sh -folder $Hp -makefits
        [ ! -e "$fqhmesh" ] &&  echored "WARNING: fqhmesh does not exist: $fqhmesh :" && fqhmesh=""
        fi
        fqhexact=$Hp/fitperatom/Fqh_fromExactFreqs_fit_order2;[ ! -e "$fqhexact" ] && echored "WARNING: fqhexact does not exist: $fqhexact :" && fqhexact=""
        [ "`getOption -v`" = "True" ] && echo "evinet  : $evinet"
        [ "`getOption -v`" = "True" ] && echo "fqhmesh : $fqhmesh "
        [ "`getOption -v`" = "True" ] && echo "fqhexact: $fqhexact"
        [ -e "$fqhmesh" ] && [ -e "$fqhexact" ] && [ -e "$evinet" ] && makethermodynamics="yes"
        [ "`getOption -v`" = "True" ] && echo maketh88:$makethermodynamics:
        fi
fi

[ "`getOption -v`" = "True" ] && echored "makefit:$makefit"
[ "`getOption -v`" = "True" ] && echored "makethermodynamics:$makethermodynamics"
[ "`getOption -v`" = "True" ] && echored "makeformation_NOchangehesse:$makeformation_NOchangehesse"
#echo -en "\r.... getting data for fit ... DONE                                                         "
echo ".... getting data for fit ... DONE                                                         "














#############################################################################################################################
# create auswertung_low auwertung_high auswertung_lowplushigh
#############################################################################################################################
[ "`getOption -s`" = "True" ] && exit
[ "$createfoldernew" = "yes" ] && rm -rf auswertung_low auswertung_low_plus_high_new auswertung_high auswertung_low_plus_high auswertung_lowplushighFahbest auswertung_lowplushigh auswertung_low_plus_high_old
#[ "$createfoldernew" = "yes" ] && rm -rf auswertung_*    no!! this deletes all "privat" auswertungen, (like auswertung_June_10th)

###############################################################################################################################
### folder = auswertung_low 
### folder = auswertung_lowplushigh
### folder = auswertung_high
###############################################################################################################################
for file in $files;do  # schleife ueber file="Fah Fah_lowplushigh Fah_high Fah_lowplushighFahbest Fah_lowplushighFahmitte"
      [ "$auswertung" = "low" ]  && [ "$file" = "Fah" ]                     && folder=auswertung_low                 
      [ "$auswertung" = "high" ] && [ "$file" = "Fah_lowplushigh" ]         && folder=auswertung_lowplushigh         
      [ "$auswertung" = "high" ] && [ "$file" = "Fah_high" ]                && folder=auswertung_high                
      [ "$auswertung" = "high" ] && [ "$file" = "Fah_lowplushighFahbest" ]  && folder=auswertung_lowplushighFahbest  
      [ "$auswertung" = "high" ] && [ "$file" = "Fah_lowplushighFahmitte" ]  && folder=auswertung_lowplushighFahmitte
      echo
      echo "#####################################################################################################################################################"
      echo "$folder (=folder) in pwd: `pwd` file: $file"
      echo "#####################################################################################################################################################"
      [ "$createfoldernew" = "yes" ] && rm -rf $folder
      mkdir -p $folder

      #echo "aLat(Ang)   T(K)  Fah(meV/at) err(meV/at)"
      #########################
            
      [ "`getOption -v`" = "True" ] && echo vor createfoldernew
      if [ "$createfoldernew" = "yes" ];then
            for i in $l; do   ## schleife ueber Ang_*K folder l=3.75Ang_1360K 3.66Ang_800K  
                              ##(in jedem folder find -Len sich files: "Fah Fah_lowplushigh Fah_high Fah_lowplushighFahbest Fah_lowplushighFahmitte"
            echo -en "\ri: $i |"
            #########################
            find=`find -L . -maxdepth 2 -mindepth 2 -type f -name $file`   # file: Fah_lowplushighFahbest Fah_high

            # continue if no Fah_xxx files found in *Ang_*K folder
            [ "$find" = "" ] && echo "$file does not exist (Fah_lowplushighFahbest will never exist is just lambda 1.0) " && continue

            fits=`cat *Ang_*K/$file | awk '{print $1}' | sort | uniq | grep "_"` ## fits:{fre,ene,eS0}_{best,cubic,l0.5fromcub,l0.5pointsnext,linear,tangens,tangenssym}
            #contributions=`echo "$fits" | sed 's|_.*||' | sort | uniq`
                  #########################
                  for fit in $fits;do  ## schleife ueber fits:fre_best fre_cubic fre_l0.5fromcub fre_l0.5pointsnext fre_linear fre_tangens fre_tangenssym
                  #echo fit: $fit
                  #########################
                  #echo -en "\rusing fit \033[31m\033[1m$fit\033[0m";echo -en "\r"
                  
                        [ ! -e "$i/$fah" ] && k=k && continue #echo no $fah in $i && continue
                        
                        [ ! -e "$i/$file" ] && continue
                        fah=`ls $i/$file |  sed 's|.*/||'`
                        is=`grep "$fit " $i/$fah | wc -l | sed 's|[ ]*||g'`
                        #echo This will not work if you just got lambda 1.0
                        ### HIER WIRD GEGREPPT: grep "fre_best "  3.75Ang_1360K/Fah_lowplushighFahbest
                        ##[ "$is" != 1 ] && echo is:$is: problem in :$fit: fah: $fah i::: $i/$file && continue   
                        ## ausgabe zu normalen auswertungszwecken unnoetig: {fre,ene,eS0}_l0.5points{next} nicht immer vorhanden
                        [ "$is" != 1 ] && continue ##echo is:$is: problem in :$fit: fah: $fah i::: $i/$file && continue
                        
                        a=`echo $i | sed 's/\(.*\)Ang_\(.*\)K[\/]*/\1/'`
                        #v=`echo $a | awk '{printf("%.2f",$1^3/'$structureFactor')}'`
                        t=`echo $i | sed 's/\(.*\)Ang_\(.*\)K[\/]*/\2/'`
                        f=`grep "$fit " $i/$fah | awk '{print $2,"   ",$3}'`          #### HIER WIRD auch GEGREPPT: grep "fre_best "  3.75Ang_1360K/Fah_lowplushighFahbest
                        
                        echo $t $f >> $folder/Fah_$fit\_$a\Ang
                        sort -g  $folder/Fah_$fit\_$a\Ang > $folder/tmp; mv  $folder/tmp  $folder/Fah_$fit\_$a\Ang
                        echo $a $f >>  $folder/Fah_$fit\_$t\K
                       
                        #echo $v $f >>  $folder/$file\_$t\K_vol
                        sort -g  $folder/Fah_$fit\_$t\K >  $folder/tmp; mv  $folder/tmp  $folder/Fah_$fit\_$t\K
                        #sort -g  $folder/$file\_$t\K_vol >  $folder/tmp; mv  $folder/tmp  $folder/$file\_$t\K_vol
                        echo $t $a $f >>  $folder/Fah_$fit\_surface
                        #echo $t $v $f >>  $folder/$file\_surface_vol
                        #echo "   $a     $t     $f"
                        done ## schleife ueber fits:fre_best fre_cubic fre_l0.5fromcub fre_l0.5pointsnext fre_linear fre_tangens fre_tangenssym
                        done ## schleife ueber Ang_*K folder
            fi

            [ "`getOption -v`" = "True" ] && echo nach createfoldernew
            echo

    #######################
    # cp files to Fah_{ene,fre,eS0}_{best,l0.5pointsnext}_changehesseref for vak and for bulk
    # create HesseChanged surface --> is only needed when there is a difference between the old and the new Fqh (not in vacancy case)
    #######################
    [ "$auswertung" = "low" ] && continue  ## in this case we dont need the fit
    [ "$folder" = "auswertung_high" ] && continue  ## in this case we dont need the fit

    # remove folder if it is empty 
    [ -e $folder ] && [ "`ls $folder | wc -w | sed 's|[ ]*||g'`" = "0" ] && rm -rf $folder && echo "fond empty folder $folder ... continue" && continue
    
    # in case of lowplus high
    # a) mkdir changedHesseref
    [ "`echo $l | wc -w | awk '{print $1}'`" = "1" ] && echo only one folder here && exit
    ###################### makefors ###################################################
    # supercell 
    # per_atom 
    # supercell_correct_vfluct 
    # per_atom_correct_vfluct  
    ###################################################################################
    if [ "`echo $folder | grep lowplushigh | wc -w | sed 's|[ ]*||g'`" = "1" ];then
    contribs="fre_best"
    makefors="supercell per_atom supercell_correct_vfluct per_atom_correct_vfluct"
    [ "$Hbd" = "def" ] && makefors="supercell"
    [ "$structnr" = "no" ] && makefors="supercell"




    #####################################
    # schleifen
    #####################################
    vorschleife=`pwd`
        for contrib in $contribs;do   # fre_best
           for makefor in $makefors;do # supercell per_atom supercell_correct_vfluct per_atom_correct_vfluct
                cd $vorschleife
                [ "$makefor" = "supercell" ] || [ "$makefor" = "supercell_correct_vfluct" ] && sctake=$Hs && atomstake=$Hn && structtake=1
                [ "$makefor" = "per_atom" ] || [ "$makefor" = "per_atom_correct_vfluct" ] && sctake=1 && atomstake=1 && structtake=$structnr

                ############################### subfolder ###############################################
                ## subfolder = auswertung_lowplushigh/fit_surface_fre_best_nochangedHesseref_per_atom
                ## subfolder
                ## subfolder
                ## subfolder
                subfolder=$folder/fit_surface_$contrib\_nochangedHesseref_$makefor
                #echo " -->  subfolder: $subfolder "

                rm -rf $subfolder
                allfiles=`ls $folder/Fah_$contrib\_*`

                echogreen "------------------------------------------------------------------------------------------------"
                echogreen "MAKEFOLDER $subfolder ... done"
                echogreen "pwd: `pwd`"
                echogreen "allfiles:$allfiles"
                echogreen "------------------------------------------------------------------------------------------------"


                #if [ "$Hbd" = "bulk" ];then
                    ### check if can be hessechanged: check it $fqh/fromfit exists
                #[ "$makefor" = "supercell_correct_vfluct" ] || [ "$makefor" = "per_atom_correct_vfluct" ] && echo "CHECK if can be Hessechanged"

                mkdir $subfolder

                ## copy all files in $subfolder 
                #################################### 
                for ii in $allfiles;do   ## $folder/Fah_{fre,ene,eS0}_{best,l0.5pointsnext}_{3.63,3.65,...}Ang {1100,1250,...}K
                      #[ "$contrib" != "fre" ] && continue
                      #echo ii: $ii
                      ##############################
                      ff=`echo $ii | sed 's|'"$folder"'/Fah_'"$contrib"'_||'`
                      echo ff:$ff subfolder:$subfolder
                      cp $ii $subfolder/Fah_$ff
                      done
                echo 4 > $subfolder/Fah_from_fit_tangens


          #######################################################################################
          #######################################################################################
          # makefit
          #######################################################################################
          #######################################################################################

            [ "`getOption -v`" = "True" ] && echo vor makefit:$makefit
            [ "$makefit" != "yes" ] && echored "     -> COULD NOT MAKE FIT .... SOME INFROMATIION MISSING ... SEE ABOVE" 
            if [ "$makefit" = "yes" ];then  # if_makefit

                echo "------------------------------------------------------------------------------------------------"
                echo "MAKEFIT (Surface) in $subfolder "
                echo "------------------------------------------------------------------------------------------------"
                fi=$subfolder/fit.input
                cp $path_mean_freqs $subfolder
                rm -f $fi
                    echo "(* adjustable parameters start *)" > $fi
                    echo "FsurfFile = \"Fah_surface\";" >> $fi
                    echo "type = 1;                                                  (*  1: T(K)  aLat(Ang/at)  F(meV/at)   *)" >> $fi
                    echo "                                                           (*  2: T(K)  V(Ang/at)     F(meV/at)   *)" >> $fi
                    echo "                                                           (*  3: T(K)  V(Ang/cell)   F(meV/cell) *)" >> $fi
                    echo "" >> $fi
                    echo "min = $rangemin;                                                (*  aLat or volume range (same format as Fsurf) for the 2. fit *)" >> $fi
                    echo "max = $rangemax;                                                (*  typically: Vmin=Veq(T=0K) and Vmax=Veq(Tmelt) *)" >> $fi
                    echo "mesh = 100;                                                (* 100 is good and should be kept *)" >> $fi
                    echo "" >> $fi
                    echo "structureFactor = $structtake;                                       (*  4: fcc  2: bcc  1: supercells *)" >> $fi
                    echo "sc = $sctake;                                                    (*  supercell, 1 for bulk *)" >> $fi
                    echo "nAtoms = $atomstake;                                               (*  1 for bulk *)" >> $fi
                    echo "" >> $fi
                    echo "fitType = \"Fvib\";                                          (*  "Fvib"  or  "poly"  fit type for 1. fit; take "Fvib" for Fah or Fel *)" >> $fi
                    echo "basis[V_, T_] := {1,T, V }                                 (*  for Fah typically: "Fvib" and {1,T,V} *)" >> $fi
                    echo "                                                           (*  for Fel typically: "Fvib" and {1,T, V,T V,T^2,V^2,V^3} *)" >> $fi
                    echo "basis2[V_]:={1, V, V^2, V^3}                               (*  should be more than sufficient: {1,V,V^2,V^3} *)" >> $fi
                    echo "" >> $fi
                    echo "minT = 1;                                                  (* typically 1     *)" >> $fi
                    echo "maxT = $tm;                                               (*           Tmelt *)" >> $fi
                    echo "stepT = 1;                                                 (*           2     *)" >> $fi
                    echo "" >> $fi
                    echo "useMeanFreqs=True;                                         (* if True "mean_freqs" file must be available; format as Fsurf, e.g. aLat(Ang) meanFreq(meV) *)" >> $fi
                    echo "                                                           (* meanFreqs are then used in the fit formula (check fitSurface.math) *)" >> $fi
                    echo "(* adjustable parameters end *)" >> $fi
                    echo "" >> $fi
                    echo '(*<<"~/scripts/Thermodynamics/fitSurface.math"*)' >> $fi
                    echo "<<\"$path/../mathematica/fitSurface.MATH\"" >> $fi


                ## make the fit
                hier=`pwd`
                cd $subfolder   # auswertung_lowplushighFahbest/fit_surface_fre_best_nochangedHesseref_per_atom
                math < fit.input
                get_fromSurface_oneVolume.sh
                echo "------------------------------------------------------------------------------------------------"
                echo ""
                    ## make getThermodynamis for bulk systems
                    echo "makethermodynamics:$makethermodynamics: makefor:$makefor:"
                    if [ "$makethermodynamics" = "yes" ] && [ "$makefor" = "per_atom" ] || [ "$makefor" = "per_atom_correct_vfluct" ];then  # only for bulk and if evinet and Fqh exists
                        echo "--------------------------------------------------------------------------------------------"
                        echo "getThermodynamics in $subfolder "
                        echo "--------------------------------------------------------------------------------------------"
                        [ ! -e "Fah_surface_fit" ] && echored "ERROR: Fah_surface_fit not created in `pwd`" && continue
                        exme="FromExactFreqs FromMeshFreqs"
                        con="fqh fqhfah"
                        before=`pwd`
                        for em in $exme;do  # "FromExactFreqs FromMeshFreqs"
                            for co in $con;do
                            [ "$em" = "FromExactFreqs" ] && fqhtake=$fqhexact
                            [ "$em" = "FromMeshFreqs" ] && fqhtake=$fqhmesh
                            cd $before
                            folderth=getThermodynamics_$em\_$co
                            echo "--> ## $folderth ## <--"
                            rm -f $folderth
                            mkdir $folderth
                            cd $folderth
                            cp $evinet .
                            cp $fqhtake Fqh  ####  Fqh / Fah
                            [ "$co" = "fqhfah" ] && cp ../Fah_surface_fit Fah
                            getThermodynamics.sh -e
                            cd $before
                            done
                        done

                    fi 

                    ## make getGibbsEnergyOfFormation.sh for defect systems
                

                cd $hier

            fi # if_makefit
        done  # for_contrib
    done # for_makefor
    fi

done ##  schleife ueber Fah_lowplushigh Fah_high Fah_lowplushighFahbest Fah_lowplushighFahmitte






#####################################################################
### getGibbsEnergyOfFormation.sh NOCHANGEHESSE
#####################################################################
###################################################
# for vak getGibbsEnergyOfFormation.sh ... 
###################################################
if [ "$Hbd" = "def" ];then
echo ""
echo "###################################################"
echo "# ti_low_6_auswertung_vacancy_NO-YES-CHANGEHESSE.sh"
echo "###################################################"
    ti_low_6_auswertung_vacancy_NO-YES-CHANGEHESSE.sh
fi
