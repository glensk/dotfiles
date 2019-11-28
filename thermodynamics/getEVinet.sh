#!/bin/bash

### was man hier noch verbessern muss  ist die automatische suche ....
### ... nach allen EVinet EVinet_b_32 EVinet_b_108 EVinet_d_xx Fqh Fqh_xx uws.

####################################################
## check necessary scripts
####################################################
path=`echo \`readlink -f "$0"\` | sed 's|'\`basename $0\`'||'`
fitVinet=$path/mathematica/fitToVinet.MATH
OUTCAR_A=OUTCAR_A.sh
OUTCAR_volumelast=OUTCAR_volume-lastexact.sh
OUTCAR_ene_lastsigma0=OUTCAR_ene-sigma0-last.sh
OUTCAR_number_of_atoms=OUTCAR_number_of_atoms.sh
OUTCAR_forces_last=OUTCAR_forces-last-max-ARRAY.sh
[ ! -e "$fitVinet" ] && echo $fitVinet does not exist && exit
[ ! -e "`which $OUTCAR_A`" ] && echo $OUTCAR_A does not exist && exit
[ ! -e "`which $OUTCAR_volumelast`" ] && echo $OUTCAR_volumelast does not exist && exit
[ ! -e "`which $OUTCAR_ene_lastsigma0`" ] && echo $OUTCAR_ene_lastsigma0 does not exist && exit
[ ! -e "`which $OUTCAR_number_of_atoms`" ] && echo $OUTCAR_number_of_atoms does not exist && exit
[ ! -e "`which $OUTCAR_forces_last`" ] && echo $OUTCAR_forces_last does not exist && exit
#[ ! -e "`which getGibbsEnergyOfFormationLink.sh`" ] && echo getGibbsEnergyOfFormationLink.sh does not exist && exit
#[ ! -e "`which getThermodynamicsLink.sh`" ] && echo getThermodynamicsLink.sh does not exist && exit
[ ! -e "`which getGibbsEnergyOfFormation.sh`" ] && echo getGibbsEnergyOfFormationLink.sh does not exist && exit
[ ! -e "`which getThermodynamics.sh`" ] && echo getThermodynamicsLink.sh does not exist && exit
echon() {
    echo ""
    echo ""
    echo "####################################################################################"
    echo "$1" `pwd`
    echo "####################################################################################"
}

####################################################
## check if bulk or defect or reference calculation
####################################################
calc=defect
[ "`pwd | grep -o "/murn_"`" != "/murn_" ] && echo couldnt find murn folder && exit
[ "`pwd | grep -o "/murn_bulk"`" = "/murn_bulk" ] && calc=bulk && echon BULK
[ "`pwd | grep -i "__REF"`" = "`pwd`" ] && calc=ref && echon REFERENCE
[ "$1" = "bulk" ] && calc=ref && echon REFERENCE
[ "$calc" = "defect" ] && echon DEFECT


####################################################
## how to run mathematica
####################################################
mathjob() {
if [ "`hostname`" = "$myhost" ];then
      math < $fitVinet
else
    echo TODO: so far mathematica scripts can only be run at $myhost
     # mathjobs=/home/glensk/jobs.cmpc
     # mathjobs_out=/home/glensk/jobs.cmpc.out
     # echo "`pwd` math < $fitVinet" >> $mathjobs
     # while [ "`wc -l $mathjobs | awk '{print $1}'`" != "0" ];do sleep 5;done;cat $mathjobs_out;rm -f $mathjobs_out 
fi
}

####################################################
## create NBANDS file (for further/ following calculations) it does not matter if this is not working
####################################################
#get.py -cn

####################################################
## create energy.dat
####################################################
  ### remove input files
  rm -f murn.in EMurn* EVinet* Murnaghan* Vinet* EVinet tmp EMurn Murn* Vinet* fitToVinet.out 
  ## remove energy files
  rm -f energy.dat energy_1.dat energy_supercell.dat energy_*.dat
  ## results
  rm -rf results
  ## create energy.dat --->>> evtl. from vasprun!!!
  $OUTCAR_A $OUTCAR_volumelast $OUTCAR_ene_lastsigma0 $OUTCAR_number_of_atoms list | tee -a energy.dat

####################################################
## create forces_last.dat
####################################################
rm -f forces_last.dat
$OUTCAR_A $OUTCAR_forces_last list | tee -a forces_last.dat
for i in `cat forces_last.dat`;do
    echo "->" $i `isnumber.sh $i`
    [ "`isnumber.sh $i`" != "yes" ] && echo no Forces in in one of the OUTCARS! && exit
    checkforce=`echo $i | awk '$1>0.001{print "TOHIGH"}'`
    [ "$checkforce" = "TOHIGH" ] && echo to high Force! not converged! && exit 
    #echo "HH" $i
done
####################################################
## check energy.dat
####################################################
[ ! -e "energy.dat" ] && echo PROBLEM: no energy.dat created && exit -1
[ "`cat energy.dat | wc -w`" = "0" ] && echo PROBLEM: energy.dat empty && exit -1
numatoms=`cat energy.dat | awk '{print $3}' | sort | uniq`
[ "`echo $numatoms | wc -w`" != "1" ] && echo numatoms: $numatoms differ! && exit -1
numrows=`cat energy.dat | awk '{print NF}' | sort | uniq`
[ "`echo $numrows`" != "3" ] && echo && echo energy.dat is wierd: "`cat energy.dat`" && exit -1

####################################################
## run Mathematica
####################################################
mathjob
[ -e "VinetFit.input" ] && echo " " >> VinetFit.input
[ -e "EVinet" ] && sed -i 's|[\t]| |g' EVinet
[ -e "EVinet" ] && echo " " >> EVinet
[ ! -e "EVinet" ] && echo PROBLEM: EVinet was not created && exit -1

####################################################
## umbenennen EVinet
###################################################
[ "$calc" = "bulk" ]   && echo EVinet does not need to be changed in name when bulk
[ "$calc" = "defect" ] && mv EVinet EVinet_d_$numatoms  
[ "$calc" = "ref" ]    && mv EVinet EVinet_b_$numatoms  
### create EVinet_1
#if [ "$numatoms" != "1" ];then
#awk '{printf "%.10f %.10f %.10f %.10f\n", $1/'"$numatoms"',$2/'"$numatoms"',$3,$4}' EVinet > EVinet_1
#fi
#cp EVinet_1 EVinet

######################################################
## split Vinet Parameters (in this way an easy konvergence can be done later)
######################################################
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

######################################################
## if bulk  :   1)  make results  folder
##              2)  getThermodynamics.sh
######################################################
if [ "$calc" = "bulk" ];then
    [ ! -e "../results__references" ] && echo ../results__references not found && exit
    ######### 1)
    rm -rf results
    mkdir results
    cp EVinet       results    ### das ist OK im /murn_bulk folder ABER NICHT im /murn_vakdivak (=ref_folder)
    cp ../results__references/* results/    ### das ist OK im /murn_bulk folder ABER NICHT im /murn_vakdivak (=ref_folder)
    ######### 2)
    echo "################### getThermodynamicsLink.sh (Bulk) ################"
    cd results;getThermodynamicsLink.sh;cd ..
fi

######################################################
## if ref  :   1)  make results  folder
##             2)  getThermodynamics.sh
######################################################
if [ "$calc" = "ref" ];then
    [ ! -e "../results__references/Fqh" ] && echo ../results__references/Fqh not found && exit
    ######### 1)
    rm -rf results
    mkdir results
    awk '{printf "%.10f %.10f %.10f %.10f\n", $1/'"$numatoms"',$2/'"$numatoms"',$3,$4}' EVinet_b_$numatoms > results/EVinet
    cp ../results__references/{Fqh,Fel,Fqh} results/    
    ######### 2)
    echo "################### getThermodynamicsLink.sh (REFERENCE) ################"
    cd results;getThermodynamicsLink.sh;cd ..
fi
######################################################
## if defect:   1) make murn for reference 
##              2) make results  folder
##              3) getGibbsEnergyOfFormation.sh
######################################################
if [ "$calc" = "defect" ];then
    ####### 1)
    def_path=`pwd`
    ref_path=`ls -1d $def_path\__[rR][eE][fF]*`
    [ ! -e "$ref_path" ] && echo "" && echo NO REFPATH FOUND && echo && exit
        cd $ref_path
        echo; echo REFERENCE:;echo `pwd`; echo
        murn_.sh ref ### nicht so, koennte in endlosschleife enden
        cd $def_path
        EVinet_refbulk=`ls -1d $ref_path/EVinet_b_*`  ### dass ist automatisch die richtige referenz: EVinet_b_32/EVinet_b_108
    [ "`echo $EVinet_refbulk | wc -w`" != "1" ] && echo "didnt find exactly one reference EVinet_b_* in $ref_path" && exit 
    [ ! -e "$EVinet_refbulk" ] && echo "didnt find exactly one reference EVinet_b_* in $ref_path" && exit
    refatoms=`echo "$EVinet_refbulk" |  sed 's|.*EVinet_b_\(.*\).*|\1|'`
    
    ###### 2) 
    rm -rf results
    mkdir results
    cp $EVinet_refbulk results
    cp EVinet_d_$numatoms  results
    echo "1 2500 2" > results/temperature_range

    [ -e "$def_path/../results__references/" ] && cp $def_path/../results__references/*_[db]_{$numatoms,$refatoms} results
    rm -f results/{Fqh,Fel,Fah}        ### die werden nicht fuer defekte gebraucht
    echo "E0:$E0:"
    E0ref=`cat $EVinet_refbulk | awk '{print $1}'`
    Eform=`echo $E0 $E0ref | awk '{print $1-('"$numatoms"'/'"$refatoms"')*$2}'`
    echo "Eform0:$Eform"
    rm -f results/kpatom_Eform0 results/cutoff_Eform0 results/atoms_Eform0
    echo "$kpatom $Eform" > results/kpatom_Eform0
    echo "$cutoff $Eform" > results/cutoff_Eform0
    echo "$numatoms $Eform" > results/atoms_EForm0


    
    ###### 2) 
    echo "################### getGibbsEnergyOfFormationLink.sh (Defect) ################"
    cd results;getGibbsEnergyOfFormationLink.sh;cd ..
    

fi
