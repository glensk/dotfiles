#!/bin/sh

conda_ase_io_folder1=`ls -1d $HOME/miniconda3/lib/python3.*/site-packages/ase/io/`
echo "conda_ase_io_folder2:$conda_ase_io_folder1:"
conda_ase_io_folder2=`ls -1d $HOME/miniconda3/envs/*/lib/python3.*/site-packages/ase/io/`
echo "conda_ase_io_folder2:$conda_ase_io_folder2:"
echo
echo
conda_ase_io_folderall="$conda_ase_io_folder1 $conda_ase_io_folder2"
#echo "conda_ase_io_folderall:$conda_ase_io_folderall:"
for conda_ase_io_folder in $conda_ase_io_folderall;do
echo $conda_ase_io_folder
done
echo
echo
for conda_ase_io_folder in $conda_ase_io_folderall;do
echo $conda_ase_io_folder

    [ "`echo "$conda_ase_io_folder" | wc -c`" -le "55" ] && echo "ERROR wrong conda folder" && exit
    [ "`echo "$conda_ase_io_folder" | wc -w`" != "1" ] && echo "ERROR found several conda ase folder" && exit
    [ "$conda_ase_io_folder" = "" ] && echo "ERROR conda folder not found" && exit
    file=$scripts/runner_scripts/ase_runner_fileformats.py
    [ ! -e "$file" ] && echo runner file $file does not exist && exit
    
    
    cp $file $conda_ase_io_folder/runner.py
    
    zeile=`grep -n "'abinit': ('ABINIT input file', '1F')," $conda_ase_io_folder/formats.py | sed 's|abinit.*||' | sed 's|:.*||'`
    #echo zeile $zeile
    #sed -i ''"$zeile"'i      'runner': (text)' $conda_ase_io_folder/formats.py 
    #sed -i '49i'"'"'mytext'"'"': 16/16' $conda_ase_io_folder/formats.py
    testit=`grep -c "'runner': ('Runner input file', '+F')," $conda_ase_io_folder/formats.py`
    [ "$testit" == "1" ] && echo "successfull! since runner is already added." && continue
    
    echo adapting $conda_ase_io_folder/formats.py
    cp $conda_ase_io_folder/formats.py $conda_ase_io_folder/formats.py.save
    sed -i '49i'"'"'runner'"'"': ('"'"'Runner input file'"'"', '"'"'+F'"'"'),' $conda_ase_io_folder/formats.py
    isone=`grep -c "'runner': ('Runner input file', '+F')," $conda_ase_io_folder/formats.py`
    [ "$isone" == "1" ] && echo "successfull!"
    [ "$isone" != "1" ] && echo THERE WAS A PROBLEM, isone: $isone, compare to original $conda_ase_io_folder/formats.py.save 

echo
done
