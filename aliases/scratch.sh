#!/bin/sh

onhost=`echo $host | sed 's/\..*$//'`
[ "`echo $onhost | grep -o cosmopc`" = "cosmopc" ] && echo /local/scratch/glensk
[ "`echo $onhost | grep -o cosmopc`" = "fidis" ] && echo /scratch/glensk 
[ "`echo $onhost | grep -o cosmopc`" = "daint" ] && echo /scratch/snx3000/aglensk
