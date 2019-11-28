#!/bin/bash

### was man hier noch verbessern muss  ist die automatische suche ....
### ... nach allen EVinet EVinet_b_32 EVinet_b_108 EVinet_d_xx Fqh Fqh_xx uws.

####################################################
## check necessary scripts
####################################################
out=no #yes #(print additional info for debugging when yes)
path=`set | grep BASH_SOURCE | sed 's|.*\"\(.*\)/[^/]*\".*|\1|' | sed '/BASH_SOURCE=/s/.*/\./'`;[ "$out" = "yes" ] && echo path: $path
script=`set | grep BASH_SOURCE | sed 's|.*\".*/\([^/]*\)\".*|\1|' | sed '/BASH_SOURCE=/s/.*/\./'`;[ "$out" = "yes" ] && echo script: $script
options=$*; . $path/utilities/functions.include; checkOptions "-h -help -v -i";[ "$out" = "yes" ] && echo options: $options

if [ `getOption -h` = True ] || [ `getOption -help` = True ]; then
  usage "`basename $0` calculate the formation energy of a defect"
  echo ""
  	  echo " - it has to be started in the defect folder"
	  echo " - it recognizes the references folder if it has the ending __reference"
      echo "    Example folder names:"
      echo "    4x4x4sc_carling_PWps_06x06x06kp_vakanz"
      echo "    4x4x4sc_carling_PWps_06x06x06kp_vakanz__reference"
  printOptions " " \
			   "-i          ignore if forces ar not converged"\
               "-v          verbose mode"
  exit
fi


fitVinet=$path/mathematica/fitToVinet.MATH
OUTCAR_A=OUTCAR_A.sh
OUTCAR_volumelast=OUTCAR_volume-lastexact.sh
OUTCAR_ene_lastsigma0=OUTCAR_ene-sigma0-last.sh
OUTCAR_number_of_atoms=OUTCAR_number_of_atoms.sh
OUTCAR_forces_last=OUTCAR_forces-last-max-ARRAY.sh
getThermodynamics=getThermodynamics.sh

checkAndSetMath

[ ! -e "$fitVinet" ] && echo $fitVinet does not exist && exit
[ ! -e "`which $OUTCAR_A`" ] && echo $OUTCAR_A does not exist && exit
[ ! -e "`which $OUTCAR_volumelast`" ] && echo $OUTCAR_volumelast does not exist && exit
[ ! -e "`which $OUTCAR_ene_lastsigma0`" ] && echo $OUTCAR_ene_lastsigma0 does not exist && exit
[ ! -e "`which $OUTCAR_number_of_atoms`" ] && echo $OUTCAR_number_of_atoms does not exist && exit
[ ! -e "`which $OUTCAR_forces_last`" ] && echo $OUTCAR_forces_last does not exist && exit
[ ! -e "`which getGibbsEnergyOfFormationLink.sh`" ] && echo getGibbsEnergyOfFormationLink.sh does not exist && exit
[ ! -e "`which $getThermodynamics`" ] && echo $getThermodynamics does not exist && exit

echon() {
    echo ""
    echo ""
    echo "####################################################################################"
    echo "$1" `pwd`
    echo "####################################################################################"
}

####################################################
## how to run mathematica
####################################################
mathjob() {
      math < $fitVinet
}

calc() {
calc=defect
[ "`pwd | grep -o "/murn_"`" != "/murn_" ] && echo couldnt find murn folder && exit
[ "`pwd | grep -o "/murn_bulk"`" = "/murn_bulk" ] && calc=bulk && echo bulk
[ "`pwd | grep -i "__REF"`" = "`pwd`" ] && calc=ref && echo ref
[ "$1" = "bulk" ] && calc=ref && echo ref
[ "$calc" = "defect" ] && echo defect
}

mmmurn() {
echo "####################################################"
echo "## check if bulk or defect or reference calculation"
echo "####################################################"
calc=`calc`
#calc=defect
#[ "`pwd | grep -o "/murn_"`" != "/murn_" ] && echo couldnt find murn folder && exit
#[ "`pwd | grep -o "/murn_bulk"`" = "/murn_bulk" ] && calc=bulk && echon BULK
#[ "`pwd | grep -i "__REF"`" = "`pwd`" ] && calc=ref && echon REFERENCE
#[ "$1" = "bulk" ] && calc=ref && echon REFERENCE
#[ "$calc" = "defect" ] && echon DEFECT

###################################################
## create NBANDS file (for further/ following calculations) it does not matter if this is not working
####################################################
#get.py -cn

echo "####################################################"
echo "## create energy.dat"
echo "####################################################"
  ### remove input files
  rm -f murn.in EMurn* EVinet* Murnaghan* Vinet* EVinet tmp EMurn Murn* Vinet* fitToVinet.out
  ## remove energy files
  rm -f energy.dat energy_1.dat energy_supercell.dat energy_*.dat
  ## results
  rm -rf results
  ## create energy.dat --->>> evtl. from vasprun!!!

  $OUTCAR_A $OUTCAR_volumelast $OUTCAR_ene_lastsigma0 $OUTCAR_number_of_atoms | awk '{print $2,$3,$4}' | head -n-1 | tee -a energy.dat
  [ ! -e "energy.dat" ] && echo PROBLEM: no energy.dat created && exit -1

echo "####################################################"
echo "## create forces_last.dat"
echo "####################################################"
rm -f forces_last.dat
$OUTCAR_A $OUTCAR_forces_last | awk '{print $2}' | head -n-1 | tee -a forces_last.dat
for i in `cat forces_last.dat`;do
    echo "->" $i `isnumber.sh $i`
    [ "`isnumber.sh $i`" != "yes" ] && echo no Forces in one of the OUTCARS! && exit
    checkforce=`echo $i | awk '$1>0.001{print "TOHIGH"}'`
    [ "$checkforce" = "TOHIGH" ] && [ "`getOption -i`" != "True" ] && echo to high Force! not converged! && exit
    #echo "HH" $i
done
echo "####################################################"
echo "## check energy.dat"
echo "####################################################"
[ "`cat energy.dat | wc -w`" = "0" ] && echo PROBLEM: energy.dat empty && exit -1
numatoms=`cat energy.dat | awk '{print $3}' | sort | uniq`
[ "`echo $numatoms | wc -w`" != "1" ] && echo numatoms: $numatoms differ! && exit -1
numrows=`cat energy.dat | awk '{print NF}' | sort | uniq`
[ "`echo $numrows`" != "3" ] && echo && echo energy.dat is wierd: "`cat energy.dat`" && exit -1

echo "####################################################"
echo "## run Mathematica"
echo "####################################################"
mathjob
[ -e "VinetFit.input" ] && echo " " >> VinetFit.input
[ -e "EVinet" ] && sed -i 's|[\t]| |g' EVinet
[ -e "EVinet" ] && echo " " >> EVinet
[ ! -e "EVinet" ] && echo PROBLEM: EVinet was not created && exit -1

echo "####################################################"
echo "## umbenennen EVinet :calc:$calc: numatoms:$numatoms:"
echo "###################################################"
[ "$calc" = "bulk" ]   && echo EVinet does not need to be changed in name when bulk
[ "$calc" = "defect" ] && mv EVinet EVinet_d_$numatoms
[ "$calc" = "ref" ]    && mv EVinet EVinet_b_$numatoms
### create EVinet_1
#if [ "$numatoms" != "1" ];then


echo "#####################################################"
echo "## split Vinet Parameters (an easy konvergence can be done later)"
echo "######################################################"
cutoff=`pwd | grep -o "[0-9]*eV" | tail -1 | sed 's|eV||'`
kpatom=`pwd | grep -o "[0-9]*x[0-9]*x[0-9]*[kK][pP]" | tail -1 | sed 's|x| |g' | sed 's|[Kk]||' | sed 's|[Pp]||' | awk '{print $1*$2*$3*'"$numatoms"'}'`
echo KP:$kpatom
B0=`grep B0GPa VinetParameters.dat | awk '{print $2}'`
B0der=`grep B0der VinetParameters.dat | awk '{print $2}'`
V0=`grep V0InAng3 VinetParameters.dat | awk '{print $2}'`
aLat=`grep aLatInAng VinetParameters.dat | awk '{print $2}'`
E0=`grep E0InmeV VinetParameters.dat | awk '{print $2}'`
echo "B0:       $B0";rm -f B0
echo "B0der:    $B0der";rm -f B0der
echo "V0:       $V0";rm -f V0
echo "aLat:     $aLat";rm -f aLat
echo "E0:       $E0";rm -f E0

rm -f cutoff_B0; echo $cutoff $B0 > cutoff_B0
rm -f cutoff_B0der; echo $cutoff $B0der > cutoff_B0der
rm -f cutoff_V0; echo $cutoff $V0 > cutoff_V0
rm -f cutoff_aLat; echo $cutoff $aLat > cutoff_aLat
rm -f cutoff_E0; echo $cutoff $E0 > cutoff_E0
rm -f kpatom_B0; echo $kpatom $B0 > kpatom_B0
rm -f kpatom_B0der; echo $kpatom $B0der > kpatom_B0der
rm -f kpatom_V0; echo $kpatom $V0 > kpatom_V0
rm -f kpatom_aLat; echo $kpatom $aLat > kpatom_aLat
rm -f kpatom_E0; echo $kpatom $E0 > kpatom_E0
rm -f atoms_E0
}

echo -------------------------------------------
mmmurn
echo -------------------------------------------
calc=`calc`
echo -------------------------------------------
echo calc:$calc:
echo -------------------------------------------

echo "######################################################"
echo "## if bulk  :   1)  make results  folder"
echo "##              2)  getThermodynamics.sh"
echo "######################################################"
if [ "$calc" = "bulk" ];then
    echo "calc is bulk"
    [ ! -e "../results__references" ] && echo ../results__references not found && exit
    ######### 1)
    rm -rf results
    mkdir results
    cp EVinet       results    ### das ist OK im /murn_bulk folder ABER NICHT im /murn_vakdivak (=ref_folder)
    cp ../results__references/* results/    ### das ist OK im /murn_bulk folder ABER NICHT im /murn_vakdivak (=ref_folder)
    ######### 2)
    echo "################### $getThermodynamics (Bulk) ################"
    cd results;$getThermodynamics;cd ..
fi

if [ "$calc" = "ref" ];then
    echo "######################################################"
    echo "## is  ref  :   1)  make results  folder"
    echo "##              2)  getThermodynamics.sh"
    echo "calc is ref"
    [ ! -e "../results__references/Fqh" ] && echo ../results__references/Fqh not found && exit
    ######### 1)
    rm -rf results
    mkdir results
    awk '{printf "%.10f %.10f %.10f %.10f\n", $1/'"$numatoms"',$2/'"$numatoms"',$3,$4}' EVinet_b_$numatoms > results/EVinet
    cp ../results__references/{Fqh,Fel,Fqh} results/
    ######### 2)
    echo "################### $getThermodynamics (REFERENCE) ################"
    cd results;$getThermodynamics;cd ..
fi



if [ "$calc" = "defect" ];then
    echo "######################################################"
    echo "## is defect:   1) make murn for reference"
    echo "##              2) make results  folder"
    echo "##              3) getGibbsEnergyOfFormation.sh"
    echo "calc is defect"
    ####### 1)
    def_path=`pwd`
    ref_path=`ls -1d $def_path\__[rR][eE][fF]*`
    [ ! -e "$ref_path" ] && echo "" && echo NO REFPATH FOUND && echo && exit
        cd $ref_path
        echo; echo REFERENCE:;echo `pwd`; echo
        #murn_.sh ref ### nicht so, koennte in endlosschleife enden
        mmmurn ref ### nicht so, koennte in endlosschleife enden
        cd $def_path
        echo --------------------------------
        echo ref_path:$ref_path
        echo --------------------------------
        echo ref_path:$ref_path
        EVinet_refbulk=`ls -1d $ref_path/EVinet_b_*`  ### dass ist automatisch die richtige referenz: EVinet_b_32/EVinet_b_108
    [ "`echo $EVinet_refbulk | wc -w`" != "1" ] && echo "didnt find exactly one reference EVinet_b_* in $ref_path" && exit
    [ ! -e "$EVinet_refbulk" ] && echo "didnt find exactly one reference EVinet_b_* in $ref_path" && exit
    refatoms=`echo "$EVinet_refbulk" |  sed 's|.*EVinet_b_\(.*\).*|\1|'`

    ###### 2)
    rm -rf results
    mkdir results
    echo "---> pwd: `pwd`"
    echo "---> EVinet_refbulk: $EVinet_refbulk"
    refatoms=`echo $EVinet_refbulk | sed 's|.*/EVinet_b_||'`
    #echo refatoms:$refatoms
    numatoms=`echo $refatoms | awk '{print $1-1}'`
    #echo "---> nua:$numatoms"

    cp $EVinet_refbulk results
    cp EVinet_d_$numatoms  results

B0=`grep B0GPa VinetParameters.dat | awk '{print $2}'`
B0der=`grep B0der VinetParameters.dat | awk '{print $2}'`
V0=`grep V0InAng3 VinetParameters.dat | awk '{print $2}'`
aLat=`grep aLatInAng VinetParameters.dat | awk '{print $2}'`
E0=`grep E0InmeV VinetParameters.dat | awk '{print $2}'`
echo "B0:       $B0";rm -f B0
echo "B0der:    $B0der";rm -f B0der
echo "V0:       $V0";rm -f V0
echo "aLat:     $aLat";rm -f aLat
echo "E0:       $E0";rm -f E0
    echo "1 2500 2" > results/temperature_range

    [ -e "$def_path/../results__references/" ] && cp $def_path/../results__references/*_[db]_{$numatoms,$refatoms} results
    rm -f results/{Fqh,Fel,Fah}        ### die werden nicht fuer defekte gebraucht
    echo "E0:$E0:"
    E0ref=`cat $EVinet_refbulk | awk '{print $1}'`
    Eform=`echo $E0 $E0ref | awk '{print $1-('"$numatoms"'/'"$refatoms"')*$2}'`
    echo :E0:$E0: E0ref:$E0ref :numatoms:$numatoms:   refatoms:$refatoms:
    echo "Eform0:$Eform"
    rm -f results/kpatom_Eform0 results/cutoff_Eform0 results/atoms_Eform0
    echo "$kpatom $Eform" > results/kpatom_Eform0
    echo "$cutoff $Eform" > results/cutoff_Eform0
    echo "$numatoms $Eform" > results/atoms_EForm0



    ###### 2)
    echo "################### getGibbsEnergyOfFormationLink.sh (Defect) ################"
    cd results;getGibbsEnergyOfFormationLink.sh;cd ..



fi
