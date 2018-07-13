#!/bin/sh
filenameall="OUTCAR vasprun.xml XDATCAR POTCAR OSZICAR ion_energies dUdLallInfo_noJumps dUdLallInfo"
for filename in $filenameall;do
        echo ""
        echo "############## $filename ############"

        files=`find . -name "$filename.gz"`
        hier=`pwd`
        for i in $files;do
            echo $i
            folder=`echo $i | sed 's|'"$filename"'.gz$||'`
            cd $hier
            cd $folder
            gunzip $filename.gz
            lbzip2 $filename
            cd $hier
        done

        files=`find . -name "$filename"`
        hier=`pwd`
        for i in $files;do
            echo $i
            folder=`echo $i | sed 's|'"$filename"'$||'`
            cd $hier
            cd $folder
            #gunzip $filename.gz
            lbzip2 $filename
            cd $hier
        done



done


