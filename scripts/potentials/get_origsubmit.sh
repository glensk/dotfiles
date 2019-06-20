#!/bin/sh
ff=`ls -1d n2p2_*`
hier=`pwd`
for f in $ff;do
    #echo
    cd $hier
    cd $f
    file="train.data"
    #[ -e "$file" ] && echo "YES $f; $file exists" && continue
    [ -e "$file" ] && continue #echo "YES $f; $file exists" && continue
    
    ################################################
    # get the original folder (original_folder)
    ################################################
    ll=`find . -maxdepth 1 -name "README_*" | wc -l`
    la=`find . -maxdepth 1 -name "README_*"`
    [ "$ll" == "0" ] && echo "NO  $f; NO README" && continue
    #echo "ll:$ll:"
    #echo "$la"
    original_folder=`cat README_* | grep "# pwd:" | awk '{print $3}' | sed 's|/potential$||' | sort | uniq`
    #echo "OF $f; original_folder: $original_folder"
    
    # check original_folder
    [ ! -d "$original_folder" ] && echo "NFF $f; NO original_folder: $original_folder :" && continue

    # copy the file
    [ -e "$original_folder/$file" ] && cp $original_folder/$file $file
    [ ! -e "$original_folder/$file" ] && echo "NOF $f; NO original_folder/file: $original_folder/$file :" && continue
    
    [ "$file" = "submit_n2p2_train.sh" ] && [ -e "$original_folder/submit_training.sh" ] && cp $original_folder/submit_training.sh $file
   
    
    #echo "-->" $f :$ll:$ll: $original_folder
    [ ! -e "$file" ] && echo "N?? $f"
    cd $hier
done


