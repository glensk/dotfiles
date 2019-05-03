#!/bin/sh
[ "$1" = "" ] && echo 'need $1 to the the foldername' && exit
saveas=""
[ "`echo $* | wc -w`" = "1" ] && saveas=$1
[ "`echo $* | wc -w`" != "1" ] && saveas=`echo $* | xargs -n1 | sed -e '1{h;d;}' -e 'G;s,\(.*\).*\n\1.*,\1,;h;$!d' | sed 's/\.$//g'`
[ "$saveas" = "" ] && saveas="archive"
[ "`echo $saveas | wc -c`" = "0" ] && saveas="archive"
[ -e "$saveas.tar.gz" ] && echo $saveas.tar.gz does already exist && exit
echo saveas $saveas
#exit
tar -zcvf $saveas.tar.gz $*

# tar from other folder (here /media/glensk@lammm/My\ Passport\ White/v/pp/) the PROJECT_eqalats_fcc_bcc_lda folder
#tar -cjf PROJECT_eqalats_fcc_bcc_lda.tar.gz -C /media/glensk@lammm/My\ Passport\ White/v/pp/ PROJECT_eqalats_fcc_bcc_lda
