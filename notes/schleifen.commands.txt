[ "$var" = "find" ] && echo found  
		INSTADS of 
		if [ "$var" = "find" ]; then
		echo "found"
		fi
-gt greaterthen (zahlen)
-lt lessthen	(zahlen)
-ne	notequal	(zahlen)
!	not			(ausdruck)
eq	equal		(zahlen)
=	equal		(strings)
-n	ist wahr wenn zeichenkette nicht leer ist

Ausdruck	Beispiel	Erklärung


-d verzeichnis	[ -d /tmp ]		Ist wahr, wenn die Datei existiert und ein Verzeichnis ist.
-f datei	[ -f txt.txt ]		Ist wahr, wenn die Datei existiert und eine normale Datei ist.


-w datei	[ -w text.txt ]		Ist wahr, wenn die Datei existiert und den Schreibzugriff erlaubt.
-x datei	[ -x script.sh ]	Ist wahr, wenn die Datei existiert und die Ausführung erlaubt.
-n string	[ -n "$name" ]		Ist wahr, wenn die übergebene Zeichenkette nicht leer ist.
str1 = str2	[ "$1" = "Hallo" ]	Ist wahr, wenn beide Zeichenketten identisch sind.
z1 -eq z2	[ 1 -eq $summe ]	Ist wahr, wenn beide Zahlen gleich groß sind 
								(in Bedingungen wird zwischen Zahlen und Zeichenketten unterschieden).
z1 -lt z2	[ 17 -lt $zahl ]	Ist wahr, wenn die erste Zahl kleiner als die zweite Zahl ist (lt = lower then).
z1 -gt z2	[ 28 -gt $tag ]		Ist wahr, wenn die erste Zahl größer als die zweite Zahl ist.
z1 -ne z2	[ $zahl -ne 7 ]		Ist wahr, wenn beide Zahlen ungleich sind.
! ausdruck	[ ! 1 -eq $zahl ]	Ist wahr, wenn der Ausdruck falsch ist (also eine Negierung).

while [[ "$bohrang" = "none"  || "$disp" = "none" ]];do
	echo wie
	break



break		verlaesst schleife
continue	faeht mit der naechsten iteration der schleife fort



########## get whole line as variable in loop
grep "^$VOL.*" POSCARLASTlist | while read LINE ; do	N=$((N+1))
		echo "Line $N = $LINE"
		forcemax=`echo "$LINE" | awk '{print $2}'`
		kp=`echo "$LINE" | awk '{print $3}'`
		kptyp=`echo "$LINE" | awk '{print $4}'`
		ENCUT=`echo "$LINE" | awk '{print $5}'`
		ENAUG=`echo "$LINE" | awk '{print $6}'`

		echo VOL:$VOL
		echo forcemax:$forcemax
		echo kp:$kp
		echo kptyp:$kptyp
		echo ENCUT:$ENCUT
		echo ENAUG:$ENAUG

		echo ""	
	done


##########################################
h=1
k=1
i=0
[ "$h" = "1" ] && echo h1
[ "$k" = "1" ] && echo k1
[[ "$k" = "1" && "$h" = "1" ]] && echo hk1
[[ "$k" = "1" && "$h" = "1" && "$i" = "0" ]] && echo hk1 i0
[[ "$k" = "1" && "$h" = "1" && "$i" = "1" ]] && echo hk1 i1
[[ "$k" = "1" || "$i" = "1" ]] && echo hi1 ---------------------> oder!


###############################################
Sometimes its convienient to have a file
read into memory to work on it.  The
form that you take to accomplish this is
an array data structure.  In ksh88 the
maximum is 1024 elements, however, on 
some of the more modern versions you can 
go much higher.

To do this the following can be done:

#/usr/bin/ksh
typeset -i cnt=0

while read line
do
  myarray[$cnt]=$line
  ((cnt = cnt + 1))
done < myfile
# end of file---------

Now, if I want to access any line of that
file, I simply use:

${<arrayname>[<subscript>]}

echo ${myarray[4]}
