#!/bin/sh
[ "$1" = "" ] && echo 'need $1 to the the foldername' && exit
[ -e "$1.tar.gz" ] && echo $1.tar.gz does already exist && exit
tar -zcvf "$1".tar.gz "$1"

# tar from other folder (here /media/glensk@lammm/My\ Passport\ White/v/pp/) the PROJECT_eqalats_fcc_bcc_lda folder
#tar -cjf PROJECT_eqalats_fcc_bcc_lda.tar.gz -C /media/glensk@lammm/My\ Passport\ White/v/pp/ PROJECT_eqalats_fcc_bcc_lda
