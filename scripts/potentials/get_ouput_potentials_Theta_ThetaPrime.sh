#!/bin/sh
#cd /scratch/glensk/daniel_comparen2p2_pots_Theta
#foldern="Theta"

cd /scratch/glensk/daniel_comparen2p2_pots_ThetaPrime
foldern="ThetaPrime"

ff=`ls -1d n2p2_alcu_v2dm_*`
hier=`pwd`
for pot in $ff;do
    cd $hier
    to="/home/glensk/Dropbox/Albert/scripts/dotfiles/scripts/potentials/$pot"
    if [ -d "$to" ];then
        #from="$pot/evinet"
        from="$pot/fqh/thermo_3rd/"
        #from="$pot/fah/Fah_surfacefit_upto1000/"
        #from1="$pot/fah/Fah_surface"
        #from2=`ls $pot/README_2020*`
        ##from3=`ls $pot/POSCAR_ThetaADPRelax`
        #from3=`ls $pot/POSCAR_ThetaPrime`
        if [ -d "$from" ];then
            echo YES $pot 
            echo YES from $from 
            echo YES to $to/$foldern/fqh
            mkdir -p $to/$foldern/fqh
            #cp -r $from $to/$foldern/fah
            #cp $from1 $to/$foldern/fah
            #cp $from2 $to/$foldern
            #cp $from3 $to/$foldern
            #cp -r $from $to/$foldern
            cp -r $from $to/$foldern/fqh
        else
            echo NO $pot
        fi
    
    else
        echo NO $pot
    fi
    echo 
    cd $hier
done


