#!/bin/sh
echo 0
find . -name "MEMCPU" -exec rm -rf {} \;
filename="vasprun.xml"
echo "now testing for both"
for j in `echo "vasprun.xml OUTCAR CONTCAR WAVECAR CHG DOSCAR EIGENVAL IBZKPT OSZICAR PCDAT"`;do
filename=$j
echo $filename
	files=`find . -name "$filename.gz"`
	for i in $files;do
	    folder=`echo $i | sed 's|'"$filename"'.gz||'`
		echo $folder
	    if [ -e "$folder/$filename" ];then # then both exist
	        echo $folder has both && ls $folder && rm $folder/$filename.gz
        else    # only *.gz exists
            gunzip $i
            zpaq a $folder/$filename.zpaq $folder/$filename 
            rm $folder/$filename
        fi

	done
done
echo "vasprun.xml"
find . -name "vasprun.xml" -exec zpaq a {}.zpaq {} \; -exec rm {} \; # rm is only run when first command successfully
echo "OUTCAR"
find . -name "OUTCAR" -exec zpaq a {}.zpaq {} \; -exec rm {} \; # rm is only run when first command successfully

echo 1
#exit
find . -name "tmp" -exec rm -rf {} \;
echo 2
find . -name "vasp.cmmd*" -exec rm -rf {} \;
find . -name "tmp2" -exec rm -rf {} \;
find . -name "CHGCAR.gz" -exec rm {} \;
find . -name "CHG.gz" -exec rm {} \;
find . -name "EIGENVAL.gz" -exec rm {} \;
find . -name "IBZKPT.gz" -exec rm {} \;
echo 1.1
find . -name "atoms_volume_steps" -exec rm {} \;
find . -name "avg_dUdL.dat" -exec gzip {} \;
find . -name "cell" -exec rm {} \;

echo 2
find . -name "dUdLallInfo" -exec gzip {} \;
find . -name "dUdLallInfo_noJumps" -exec gzip {} \;
find . -name "dUdL" -exec gzip {} \;
find . -name "getAvgAndDev_input_ene" -exec gzip {} \;
find . -name "getAvgAndDev_input_eS0" -exec gzip {} \;
find . -name "getAvgAndDev_input_fre" -exec gzip {} \;
find . -name "HesseMatrix_sphinx" -exec gzip {} \;
echo 3
find . -name "ion_energies" -exec gzip {} \;
find . -name "jumps_dist" -exec gzip {} \;
find . -name "jumps_radius" -exec gzip {} \;

find . -name "equilibrium" -exec rm {} \;

echo 4
find . -name "da_0.0" -exec rm {} \;
find . -name "da_1.0" -exec rm {} \;
find . -name "dan_0.0" -exec rm {} \;
find . -name "dan_1.0" -exec rm {} \;
find . -name "di_0.0" -exec rm {} \;
find . -name "di_1.0" -exec rm {} \;
find . -name "din_0.0" -exec rm {} \;
find . -name "din_1.0" -exec rm {} \;
find . -name "do_0.0" -exec rm {} \;
find . -name "do_1.0" -exec rm {} \;
find . -name "don_0.0" -exec rm {} \;
find . -name "don_1.0" -exec rm {} \;

echo "dfn_0.15"
find . -name "dfn_0.0" -exec rm {} \;
find . -name "dfn_0.15" -exec rm {} \;
find . -name "dfn_0.5" -exec rm {} \;
find . -name "dfn_0.85" -exec rm {} \;
find . -name "dfn_1.0" -exec rm {} \;

echo "pos_0.0"
find . -name "pos_0.0" -exec rm {} \;
find . -name "pos_0.15" -exec rm {} \;
find . -name "pos_0.5" -exec rm {} \;
find . -name "pos_0.85" -exec rm {} \;
find . -name "pos_1.0" -exec rm {} \;

echo "dln_0.15"
find . -name "dln_0.0" -exec rm {} \;
find . -name "dln_0.15" -exec rm {} \;
find . -name "dln_0.5" -exec rm {} \;
find . -name "dln_0.85" -exec rm {} \;
find . -name "dln_1.0" -exec rm {} \;

echo ".cmbackup1"
find . -name ".cmbackup*" -exec rm -rf {} \;
echo ".cmbackup2"
find . -name "cmbackup*" -exec rm -rf {} \;

echo "CHG"
find . -name "CHG" -exec rm {} \;
find . -name "CHG.gz" -exec rm {} \;
echo "DOSCAR"
find . -name "DOSCAR" -exec rm {} \;
find . -name "DOSCAR.gz" -exec rm {} \;
echo "EIGENVAL"
find . -name "EIGENVAL" -exec rm {} \;
find . -name "EIGENVAL.gz" -exec rm {} \;
echo "IBZKPT"
find . -name "IBZKPT" -exec rm {} \;
find . -name "IBZKPT.gz" -exec rm {} \;
echo "PCDAT"
find . -name "PCDAT" -exec rm {} \;
find . -name "PCDAT.gz" -exec rm {} \;
echo "OSZICAR"
find . -name "OSZICAR" -exec rm {} \;
echo "OSZICAR.gz"
find . -name "OSZICAR.gz" -exec rm {} \;
echo "WAVECAR"
find . -name "WAVECAR" -exec rm {} \;
echo "WAVECAR.gz"
find . -name "WAVECAR.gz" -exec rm {} \;
echo "meinjob*"
find . -name "meinjob*.o*" -exec rm {} \;
find . -name "meinjob*.po*" -exec rm {} \;




# delete stuff; this is easy
echo "POSITIONs"
find . -name "POSITIONs" 
find . -name "POSITIONs" -exec rm {} \;
find . -name "POSITIONs.gz" -exec rm {} \;

filename="vasprun.xml"
echo "now testing for both"
for j in `echo "vasprun.xml OUTCAR CONTCAR WAVECAR CHG DOSCAR EIGENVAL IBZKPT OSZICAR PCDAT"`;do
filename=$j
echo $filename
	files=`find . -name "$filename.gz"`
	for i in $files;do
	    folder=`echo $i | sed 's|'"$filename"'.gz||'`
		echo $folder
	    [ -e "$folder/$filename" ] && echo $folder has both && rm $folder/$filename.gz 
	done
done

echo "pre_equilibration.gz"
find . -name "pre_equilibration.gz" 
find . -name "pre_equilibration.gz" -exec rm {} \;

echo "pre_equilibration"
find . -name "pre_equilibration" 
find . -name "pre_equilibration" -exec rm {} \;

echo "structures_vasprun.gz"
find . -name "structures_vasprun.gz"
find . -name "structures_vasprun.gz" -exec rm {} \;

echo "structures_vasprun"
find . -name "structures_vasprun" 
find . -name "structures_vasprun" -exec rm {} \;

find . -name "POTCAR" -exec zpaq a {}.zpaq {} \; -exec rm {} \; # rm is only run when first command successfully
find . -name "CONTCAR" -exec zpaq a {}.zpaq {} \; -exec rm {} \; # rm is only run when first command successfully



####################################################################
# hier sollte man nicht erst alles auspacken und dann wieder packen
####################################################################


#files=`find . -name "vasprun.xml.gz"`
#for i in $files;do
#    filegz=$i
#    file=`echo $i | sed 's|.gz$||'`
#    echo $filegz
#    echo $file
#    exit
#done

