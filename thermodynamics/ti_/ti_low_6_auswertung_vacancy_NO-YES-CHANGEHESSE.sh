#!/bin/sh

out=noyes #(print additional info for debugging when yes)
path=`set | grep BASH_SOURCE | sed 's|.*\"\(.*\)/[^/]*\".*|\1|' | sed '/BASH_SOURCE=/s/.*/\./'`;[ "$out" = "yes" ] && echo path: $path
script=`set | grep BASH_SOURCE | sed 's|.*\".*/\([^/]*\)\".*|\1|' | sed '/BASH_SOURCE=/s/.*/\./'`;[ "$out" = "yes" ] && echo script: $script
options=$*; . $path/../utilities/functions.include; checkOptions "-h -help -v -d -s -n -cb -cd -cbd -nvb -nvv";[ "$out" = "yes" ] && echo options: $options

if [ `getOption -h` = True ] || [ `getOption -help` = True ] ; then
  usage `basename $0`
  echo "This script reads the files auswertung_NOCHANGEHESSE_input CHANGEHESSEinput_bulk CHANGEHESSEinput_defect and creates corresponding folder";
  printOptions " " \
               "-n          make only auswertung_NOCHANGEHESSE" \
               "-nvb        make only auswertung_NOCHANGEHESSE_vcorrbulk" \
               "-nvv        make only auswertung_NOCHANGEHESSE_vcorrvak" \
               "-cb         make only auswertung_CHANGEHESSE_bulk"\
               "-cd         make only auswertung_CHANGEHESSE_defect"\
               "-cbd        make only auswertung_CHANGEHESSE_bulk_defect"\
               "-s          show fit data and exit" \
               "-v          verbose mode" \ 
  exit
fi
checkAndSetMath
options_all="-n -cb -cd -bd -nvb -nvv $options"
[ "`getOption -n`" != "True" ] && [ "`getOption -cb`" != "True" ] && [ "`getOption -cd`" != "True" ] && [ "`getOption -cbd`" != "True" ] && [ "`getOption -nvb`" != "True" ] && [ "`getOption -nvv`" != "True" ] && options="$options_all"
hier=`pwd`

##################################################################################
## check input files NOCHANGEHESSEinput
##################################################################################
if [ ! -f "auswertung_NOCHANGEHESSE_input" ];then   
    echored "ERROR: file auswertung_NOCHANGEHESSE_input not found! This is prerequisit for auswertung_NOCHANGEHESSE and auswertung_CHANGEHESSE"
    echored "CONSIDER creating auswertung_NOCHANGEHESSE_input file"
    echo "  EVinet_b_$refatoms path_to_file"
    echo "  EVinet_d_$Hn path_to_file"
    echo ""
    echo "  Fqh_b_$refatoms    path_to_file"
    echo "                  (Fqh_b_$refatoms: path to HesseMatrix folder of the bulk thermodynamic Integration (TI) ....."
    echo "                   ..... it is important that this is the original HesseMatrix bulk TI folder.)"
    echo "  Fah_b_$refatoms    path_to_file"
    echo "                  (path to the auswertung_lowplushighFahbest/fit_surface_fre_best_nochangedHesseref_supercell/Fah_surface_fit)"
    echo ""
    echo "  Fqh_d_$Hn path_to_file"
    echo "  Fah_d_$Hn (not necessary, is known)"
    exit 
else
    echogreen "auswertung_NOCHANGEHESSE_input found .... reading"
    evinet_b=`grep EVinet_b auswertung_NOCHANGEHESSE_input | awk '{print $2}'`
    evinet_d=`grep EVinet_d auswertung_NOCHANGEHESSE_input | awk '{print $2}'`
    fqh_d=`grep Fqh_d auswertung_NOCHANGEHESSE_input | awk '{print $2}'`
    fqh_b=`grep Fqh_b auswertung_NOCHANGEHESSE_input | awk '{print $2}'`
    fah_b=`grep Fah_b auswertung_NOCHANGEHESSE_input | awk '{print $2}'`
    fah_d=`grep Fah_d auswertung_NOCHANGEHESSE_input | awk '{print $2}'`
    if [ "`getOption -v`" = "True" ];then
    echo " "
    echo "evinet_b:$evinet_b"
    echo "evinet_d:$evinet_d"
    echo " "
    echo "   fqh_b:$fqh_b"
    echo "   fah_b:$fah_b"
    echo ""
    echo "   fqh_d:$fqh_d"
    echo "   fah_d:$fah_d"
    echo ""
    fi
    ## NOTES: Fqh_d and Fah_b SHOUD NOT BE SEARCHED AUTOMATICALLY; YOU WANT TO HAVE CONSISTENT VOLUME RANGES FOR FQH_B and FQH_D
    [ ! -f "$evinet_b" ] && echored "  --> ERROR: evinet_b in auswertung_NOCHANGEHESSE: $evinet_b does not exist;" && rmOption -n 
    [ ! -f "$evinet_d" ] && echored "  --> ERROR: evinet_d in auswertung_NOCHANGEHESSE: $evinet_d does not exist;" && rmOption -n
    [ ! -f "$fqh_b" ] &&    echored "  --> ERROR: fqh_b in auswertung_NOCHANGEHESSE: $fqh_b does not exist;" && rmOption -n
    [ ! -f "$fah_b/Fah_surface_fit" ] &&    echored "  --> ERROR: fah_b in auswertung_NOCHANGEHESSE: $fah_b does not exist;" && rmOption -n
    [ ! -f "$fqh_d" ] &&    echored "  --> ERROR: fqh_d in auswertung_NOCHANGEHESSE: $fqh_d does not exist;" && rmOption -n
    [ ! -f "$fah_d/Fah_surface_fit" ] &&    echored "  --> ERROR: fah_d in auswertung_NOCHANGEHESSE: $fah_d does not exist;" && rmOption -n
fi
#[ "`getOption -n`" != "True" ] && echored "SOME PROBLEM WITH auswertung_NOCHANGEHESSE .... cant go on" && exit

#############################################################################################################################
# get some information for fitting
#############################################################################################################################
    makefit="yes"

    # .... getting data for fit ... 1. hessepath"
    echo -en "\r.... getting data for fit ... [1/3]. hessepath"
    Hp=`ti_high_0_create_Folders_vasp.sh -Hp`
    [ ! -e "$Hp" ] && makefit="no hessepath" && echored "COULD NOT FIND Hessepath:$Hp: (for fit surface)"
    [ "`getOption -v`" = "True" ] && echo "Hp:$Hp                  "

    if [ "$makefit" = "yes" ];then
    echo -en "\r.... getting data for fit ... [2/3]. number of atoms"
    Hn=`ti_high_0_create_Folders_vasp.sh -Hn`
    [ "`checkInteger $Hn`" != "ok" ] && makefit="no Hn" && echo && echored "COULD NOT FIND number of atoms:$Hn: (for fit surface) in $Hp/parameters.dat" && echo
    [ "`getOption -v`" = "True" ] && echo "Hn:$Hn                     "
    fi

    # .... getting data for fit ... 6. bulk or defect"
    if [ "$makefit" = "yes" ];then
    echo -en "\r.... getting data for fit ... [3/3]. bulk or defect"
    Hbd=`ti_high_0_create_Folders_vasp.sh -Hbd`
    [ "`getOption -v`" = "True" ] && echo "Hbd:$Hbd                 "
    fi
    [ "$Hbd" != "def" ] && echored "ERROR: This does not seem to be a defect structure! Hbd=$Hbd;" && exit

[ "`getOption -v`" = "True" ] && echored "makefit:$makefit"
[ "`getOption -v`" = "True" ] && echored "makethermodynamics:$makethermodynamics"
[ "`getOption -v`" = "True" ] && echored "makeformation_NOchangehesse:$makeformation_NOchangehesse"
echo -en "\r.... getting data for fit ... DONE     "

#####################################################################
### auswertung_NOCHANGEHESSE
#####################################################################
refatoms=` expr $Hn + 1 `
cd $hier
if [ "`getOption -n`" = "True" ];then
                fnoch=auswertung_NOCHANGEHESSE
                echo ""
                echogreen "------------------------------------------------------------------------------------------------"
                echogreen ".... $fnoch makefolder"
                echogreen "------------------------------------------------------------------------------------------------"
                rm -rf $fnoch
                mkdir $fnoch
                cp $evinet_b $fnoch/EVinet_b_$refatoms
                cp $evinet_d $fnoch/EVinet_d_$Hn
                cp $fqh_b $fnoch/Fqh_b_$refatoms
                cp $fah_b/Fah_surface_fit $fnoch/Fah_b_$refatoms
                cp $fqh_d $fnoch/Fqh_d_$Hn
                cp $fah_d/Fah_surface_fit $fnoch/Fah_d_$Hn


                echo EVinet_b_$refatoms $evinet_b > $fnoch/FOLDER
                echo EVinet_d_$Hn $evinet_d >> $fnoch/FOLDER
                echo Fqh_b_$refatoms $fqh_b >> $fnoch/FOLDER
                echo Fah_b_$refatoms $fah_b/Fah_surface_fit >> $fnoch/FOLDER
                echo Fqh_d_$Hn $fqh_d >> $fnoch/FOLDER
                echo Fah_d_$Hn $fah_d/Fah_surface_fit >> $fnoch/FOLDER
                zuvor=`pwd`
                cd $fnoch
                echo "------------------------------------------------------------------------------------------------"
                echo .... $fnoch getGibbsEnergyOfFormation
                echo "------------------------------------------------------------------------------------------------"
                getGibbsEnergyOfFormationLink.sh
                cd $zuvor
fi


##################################################################################
## check input files CHANGEHESSEinput_{bulk,vak}
##################################################################################
[ "`getOption -v`" = "True" ] && echo 0:$options:
[ "`getOption -cb`" = "True" ] && [ ! -e "auswertung_CHANGEHESSE_bulk_input" ] && echored "auswertung_CHANGEHESSE_bulk_input not found found" && rmOption -cb
[ "`getOption -v`" = "True" ] && echo 1:$options:
if [ "`getOption -cb`" = "True" ];then
        [ "`getOption -v`" = "True" ] && echo inschleife CHANGEHESSEinput_bulk options:$options:
        echogreen "auswertung_CHANGEHESSE_bulk_input found .... reading"
        oldb=`grep oldFqh_b auswertung_CHANGEHESSE_bulk_input | awk '{print $2}'`
        newb=`grep newFqhfolder_b auswertung_CHANGEHESSE_bulk_input | awk '{print $2}'`
        fqh_b_new=`grep newFqhfit_b auswertung_CHANGEHESSE_bulk_input | awk '{print $2}'`
        [ ! -d "$oldb" ] && echored "   --> WARNING: oldFqh_b (directory with Fqh_classical) does not esist in auswertung_CHANGEHESSE_bulk_input file: $oldb :" && rmOption -cb 
        [ ! -d "$newb" ] && echored "   --> WARNING: newFqh_b (directory with Fqh_classical files) does not esist in auswertung_CHANGEHESSE_bulk_input file: $newb :" && rmOption -cb 
        [ ! -f "$fqh_b_new" ] && echored "   --> WARNING: newFqh_b (file new quasiharmonic surface) does not esist in auswertung_CHANGEHESSE_bulk_input file: $fqh_b_new :" && rmOption -cb 
        [ ! -f "$fah_b/fit.input" ] && echored "   --> WARNING: $ahb/fit.input does not exist" && rmOption -cb
        [ -f "$fah_b/HesseRef_changed" ] && echored "   --> WARNING: $ahb/HesseRef_changed does already exist" && rmOption -cb 
fi


[ "`getOption -v`" = "True" ] && echo 2:options:
[ "`getOption -cd`" = "True" ] && [ ! -e "auswertung_CHANGEHESSEinput_defect" ] && echored "auswertung_CHANGEHESSEinput_defect not found found" && rmOption -cd
[ "`getOption -v`" = "True" ] && echo 3:options:
if [ "`getOption -cd`" = "True" ];then
        [ "`getOption -v`" = "True" ] && echo inschleife CHANGEHESSEinput_defect options:$options:
        echogreen "auswertung_CHANGEHESSEinput_defect found .... reading"
        oldd=`grep oldFqh_d CHANGEHESSEinput_defect | awk '{print $2}'`
        newd=`grep newFqh_d CHANGEHESSEinput_defect | awk '{print $2}'`
        ahd=`grep Fah_d CHANGEHESSEinput_defect | awk '{print $2}'`
        [ ! -e "$oldd" ] && echored "WARNING: oldFqh_d does not esist in CHANGEHESSEinput_bulk file: $oldd :" && rmOption -cd 
        [ ! -e "$newd" ] && echored "WARNING: newFqh_d does not esist in CHANGEHESSEinput_bulk file: $newd :" && rmOption -cd 
        [ ! -e "$ahd" ] && echored "WARNING: Fah_d does not esist in CHANGEHESSEinput_bulk file: $ahd :"      && rmOption -cd 
        [ ! -d "$ahd" ] && echored "WARNING: Fah_d in CHANGEHESSEinput_bulk file is not a directory: $ahd :"  && rmOption -cd 
fi

hier=`pwd`
#[ "`getOption -s`" = "True" ] && options="$options -v" 
#[ "`getOption -v`" = "True" ] && echo options:$options:


#####################################################################
### auswertung_NOCHANGEHESSE_vcorrbulk
#####################################################################
cd $hier

if [ ! -e "auswertung_NOCHANGEHESSE_vcorrbulk_input" ];then
    echored "auswertung_NOCHANGEHESSE_vcorrbulk_input not found"
else
    ob=`grep old auswertung_NOCHANGEHESSE_vcorrbulk_input | awk '{print $2}'`
    nb=`grep new auswertung_NOCHANGEHESSE_vcorrbulk_input | awk '{print $2}'`
if [ "`getOption -nvb`" = "True" ] ;then
                fnoch=auswertung_NOCHANGEHESSE_vcorrbulk
                echo ""
                echogreen "------------------------------------------------------------------------------------------------"
                echogreen ".... $fnoch makefolder"
                echogreen "------------------------------------------------------------------------------------------------"
                [ ! -e "$ob" ] && echored "$ob does not exist"
                [ ! -e "$nb" ] && echored "$nb does not exist"
                if [ -e "$ob" ] && [ -e "$nb" ];then
                rm -rf $fnoch
                mkdir $fnoch
                cd $fnoch
                zuvor=`pwd`
                mkdir make_changehesse_bulk
                cd make_changehesse_bulk
                echo fah_b:$fah_b
                cp $fah_b/* . >& /dev/null
                echo cp done
                rm -f Fah_surface_fit
                echo ob: $ob
                echo nb: $nb
                changeHesseReference.sh $ob $nb -qm
                getitfrom=`pwd`/Fah_surface_fit
                if [ ! -e "HesseRef_changed" ];then
                    echored "HesseRef_changed was not created"
                else
                    math < fit.input
                    get_fromSurface_oneVolume.sh
                fi

                cd $hier
                
                if [ -f "$fnoch/make_changehesse_bulk/Fah_surface_fit" ];then
                cp $evinet_b $fnoch/EVinet_b_$refatoms
                cp $evinet_d $fnoch/EVinet_d_$Hn
                cp $fqh_b $fnoch/Fqh_b_$refatoms
                cp $getitfrom $fnoch/Fah_b_$refatoms
                cp $fqh_d $fnoch/Fqh_d_$Hn
                cp $fah_d/Fah_surface_fit $fnoch/Fah_d_$Hn


                echo EVinet_b_$refatoms $evinet_b > $fnoch/FOLDER
                echo EVinet_d_$Hn $evinet_d >> $fnoch/FOLDER
                echo Fqh_b_$refatoms $fqh_b >> $fnoch/FOLDER
                echo Fah_b_$refatoms $getitfrom >> $fnoch/FOLDER
                echo Fqh_d_$Hn $fqh_d >> $fnoch/FOLDER
                echo Fah_d_$Hn $fah_d/Fah_surface_fit >> $fnoch/FOLDER
                else
                    echored "$fnoch/make_changehesse_bulk/Fah_surface_fit does not exist"
                fi

                echo "------------------------------------------------------------------------------------------------"
                echo .... $fnoch getGibbsEnergyOfFormation
                echo "------------------------------------------------------------------------------------------------"
                cd $fnoch
                getGibbsEnergyOfFormationLink.sh
                cd $zuvor
            fi
fi
fi


#####################################################################
### auswertung_NOCHANGEHESSE_vcorrvak
#####################################################################
cd $hier

if [ ! -e "auswertung_NOCHANGEHESSE_vcorrvak_input" ];then
    echored "auswertung_NOCHANGEHESSE_vcorrvak_input not found"
else
    ob=`grep old auswertung_NOCHANGEHESSE_vcorrvak_input | awk '{print $2}'`
    nb=`grep new auswertung_NOCHANGEHESSE_vcorrvak_input | awk '{print $2}'`
if [ "`getOption -nvv`" = "True" ] ;then
                fnoch=auswertung_NOCHANGEHESSE_vcorrvak
                echo ""
                echogreen "------------------------------------------------------------------------------------------------"
                echogreen ".... $fnoch makefolder"
                echogreen "------------------------------------------------------------------------------------------------"
                [ ! -e "$ob" ] && echored "$ob does not exist"
                [ ! -e "$nb" ] && echored "$nb does not exist"
                if [ -e "$ob" ] && [ -e "$nb" ];then
                rm -rf $fnoch
                mkdir $fnoch
                cd $fnoch
                zuvor=`pwd`
                mkdir make_changehesse_bulk
                cd make_changehesse_bulk
                echo fah_d:$fah_d
                cp $fah_d/* . >& /dev/null
                echo -------------
                ls
                echo -------------
                echo cp done
                rm -f Fah_surface_fit
                echo ob: $ob
                echo nb: $nb
                changeHesseReference.sh $ob $nb -qm
                getitfrom=`pwd`/Fah_surface_fit
                if [ ! -e "HesseRef_changed" ];then
                    echored "HesseRef_changed was not created"
                else
                    math < fit.input
                    get_fromSurface_oneVolume.sh
                fi

                cd $hier
                
                if [ -f "$fnoch/make_changehesse_bulk/Fah_surface_fit" ];then
                cp $evinet_b $fnoch/EVinet_b_$refatoms
                cp $evinet_d $fnoch/EVinet_d_$Hn
                cp $fqh_b $fnoch/Fqh_b_$refatoms
                cp $fah_b/Fah_surface_fit $fnoch/Fah_b_$refatoms
                cp $fqh_d $fnoch/Fqh_d_$Hn
                cp $getitfrom $fnoch/Fah_d_$Hn


                echo EVinet_b_$refatoms $evinet_b > $fnoch/FOLDER
                echo EVinet_d_$Hn $evinet_d >> $fnoch/FOLDER
                echo Fqh_b_$refatoms $fqh_b >> $fnoch/FOLDER
                echo Fah_b_$refatoms $fah_b/Fah_surface_fit >> $fnoch/FOLDER
                echo Fqh_d_$Hn $fqh_d >> $fnoch/FOLDER
                echo Fah_d_$Hn $getitfrom >> $fnoch/FOLDER
                else
                    echored "$fnoch/make_changehesse_bulk/Fah_surface_fit does not exist"
                fi

                echo "------------------------------------------------------------------------------------------------"
                echo .... $fnoch getGibbsEnergyOfFormation
                echo "------------------------------------------------------------------------------------------------"
                cd $fnoch
                getGibbsEnergyOfFormationLink.sh
                cd $zuvor
            fi
fi
fi

#####################################################################
### auswertung_CHANGEHESSE_bulk
#####################################################################
cd $hier
if [ "`getOption -cb`" = "True" ] ;then
                fnoch=auswertung_CHANGEHESSE_bulk
                echo ""
                echogreen "------------------------------------------------------------------------------------------------"
                echogreen ".... $fnoch makefolder"
                echogreen "------------------------------------------------------------------------------------------------"
                rm -rf $fnoch
                mkdir $fnoch
                cd $fnoch
                zuvor=`pwd`
				echo 1 $zuvor
                mkdir make_changehesse_bulk
                cd make_changehesse_bulk
				echo 2 `pwd`
                cp $fah_b/* . >& /dev/null
                rm -f Fah_surface_fit
				echo 3 ka
                changeHesseReference.sh $oldb $newb
				echo 4 kb
                getitfrom=`pwd`/Fah_surface_fit
                if [ ! -e "HesseRef_changed" ];then
                    echored "HesseRef_changed was not created"
                else
                    math < fit.input
                    get_fromSurface_oneVolume.sh
                fi

                cd $hier
                
                if [ -f "$fnoch/make_changehesse_bulk/Fah_surface_fit" ];then
                    cp $evinet_b $fnoch/EVinet_b_$refatoms
                    cp $evinet_d $fnoch/EVinet_d_$Hn
                    cp $fqh_b_new $fnoch/Fqh_b_$refatoms
                    cp $getitfrom $fnoch/Fah_b_$refatoms
                    cp $fqh_d $fnoch/Fqh_d_$Hn
                    cp $fah_d/Fah_surface_fit $fnoch/Fah_d_$Hn
                else
                    echored "$fnoch/make_changehesse_bulk/Fah_surface_fit does not exist"
                fi

                echo "------------------------------------------------------------------------------------------------"
                echo .... $fnoch getGibbsEnergyOfFormation
                echo "------------------------------------------------------------------------------------------------"
                cd $fnoch
                getGibbsEnergyOfFormationLink.sh
                cd $zuvor
fi
