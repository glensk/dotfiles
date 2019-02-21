#!/bin/sh
hier=`pwd`
tarfiles="KMC_QCACHE KMC_ECACHE KMC_AL6XXX"
for tar in $tarfiles;do
    echo
    echo tar $tar
    files=`find . -name $tar`
    for i in $files;do
        cd $hier
        folder=`echo $i | sed 's|/'"$tar"'$||g'`
        cd $hier
        cd $folder
        pwd
        [ "$tar" == "KMC_AL6XXX" ] && [ -e "KMC_AL6XXX" ] && [ -e "KMC_AL6XXX_out" ] && rm KMC_AL6XXX_out   # just if previous KMC_AL6xxx_out existed 
        [ "$tar" == "KMC_AL6XXX" ] && [ -e "KMC_AL6XXX" ] && kmc_get_vacancy_radius_vs_time.sh 
        [ "$tar" == "KMC_AL6XXX" ] && awk '{print $1*2.4e-17,$3}' KMC_AL6XXX > KMC_AL6XXX_out
        rm -f .sendx_*

        #du -sh $tar 
        tar.bzip2_file_or_folder_and_remove_originals.sh $tar
        cd $hier
    done
done

# all log files
files=`find . -name log.i-pi`
for i in $files;do
    cd $hier
    folder=`echo $i | sed 's|/log.i-pi$||g'`
    cd $hier
    cd $folder
    pwd
    tar.bzip2_file_or_folder_and_remove_originals.sh log*
    cd $hier
done

