comvert pic.jpg pic.png
yum search pyyaml


********************     PHONONS  *****************************************
/home/grabowski/SFHIngX/structures			structures (hcp)

**********************shows .jpg, .tiff, .bmp files
xv file			shows .jpg, .tiff, .bmp files

*************** scp ************************
scp albert@cecampc22:fullerene/g1.xyz g1.xyz   (in own shell) (copy file via ssh)
scp (-r) user@host:directory/SourceFile TargetFile    (vom Fremden pc auf den gerade eingeloggten -r für ganzen ordner)

************* run command in different shell *****************
/bin/bash -c "echo -e fg31 '\E['31';'01'm SS64'; tput sgr0"

************* uppercase lowercase ****************************
$ echo $VAR_NAME | tr '[:upper:]' '[:lower:]'
$ echo $VAR_NAME | tr '[:lower:]' '[:upper:]'
tr '[A-Z]' '[a-z]'

************ wlan passwort ***********************************
wlan passwort		:		texturtexturx

*********************************#check mem during run ************
ssh knotenname top -b -n1

******************************************************
*************** OTHER ********* OTHER  ***************
******************************************************

sshfs -o idmap=user glensk@cecampc4:/home/glensk ~/mnt/cecampc4  (mounts cecampc4)
fusermount -u ~/mnt/cecampc4                                    (umounts cecampc4)
http://voku-online.de/comment.php?comment.news.130    (hat darueber die Infos)

sshfs -o idmap=user albert@cecampc22:/home/albert ~/mnt/cecampc22

ls -1d *
ls -1d */
ls -1d */*/

cate    <ordner>     					suche <ordnername>     
ldd ~/scripts/vasp-4.6-par_mpich   		Zeigt geladene module an
ldd cpmd.x      						zeigt an welche libs er laden kann in einer dynamical executable
eog     eye of Gnome Bildbetrachtung

zgrep text file

rpm -i file 							install rpm package automatically

cat inputfile | tr '4' '3'         		ersetzt alle 4er durch 3er
paste file1 file2 file3    				fuegt files in eine file
source .bashrc    oder   . .bashrc    	laed die bashrc

cat /proc/cpuinfo | grep processor		#how many cpu cpu's processors

ssh cmdft015 top -b -n1 | head -20


print next to each other
pr -s -t -m file1 file2

seq 1 24				range list of numbers

df -h		schreibe speicherverbrauch von /home/ /raid/ usw

gvimdiff	diff like kompare 

/home/glensk/.snapshot/... sicherung alles ohne links!

nslookup cmmd002 gibt die ip adresse aus

netstat zeigt welche ports offen sind

/etc/ssh/sshd_config  ==> UseDNS no hinzufuegen stellt reverselookup der ip per ssh ab


module load vasp/parallel/5.2.11
mpirun -np 24 vasp


##############33
find . -maxdepth 1 -mindepth 1 -type d -exec sh -c 't=${0%/*}/$(printf %s "${0##*/}" | tr "[:upper:]" "[:lower:]");[ "$t" = "$0" ] || mv -i "$0" "$t"' {} \;
 --> make all subfolder lowercase

 l=`ls exactFreqs_3.*`; for i in $l; do a=`echo $i | sed 's|exactFreqs_||'`; echo $a `cat $i | xargs`>> exact_vs_aLat; done


$cat /proc/version   # to get version of system
uname -m
Linux version 3.9.6-200.fc18.x86_64 (mockbuild@bkernel02) (gcc version 4.7.2 20121109 (Red Hat 4.7.2-8) (GCC) ) #1 SMP Thu Jun 13 18:56:55 UTC 2013

