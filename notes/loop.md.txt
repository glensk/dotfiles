to loop over lines:
IFSbefore=$IFS
IFS='
'
for x in `ls -l $1`; do echo $x; done
IFS=$IFSbefore
