******** LINKS ************

ln -s path_to_realfile/realfolder LINKname

create link in project folder after deleted job
ln -s /home/glensk/PHInaX/jobs/001/job.823/ job.001-823

rm job.001-823/ does not work
rm job.001-823 does work but removes the content of the parent "real folder" as well!!!!!!!!!!!

unlink job.001-823/ does not work
unlink job.001-823  does the trick

readling link shows target of link:
readlink -f "$0"


ln -s -f target linkname   # to change target of a link without deleting the link

check if variable is defined:   [ ! -z "$variable" ] is true if variable is defined
                                [ -z "$variable" ]   is true if variable is NOT defined   
-e file exists
-f file is regular file
-d file is directory
-h file is a symbolic link
-L file is a symbolic link
-r file has read permissions
-w file has write permissions
-x file has execute permissions
f1 -nt f2 f1 is newer than f2
fq -ot f2 f1 is older than f2
! not (reverse stuff)
