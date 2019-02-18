sed is case sensitive
sed 's|/| |g' | awk '{print $1}' ## to print first second third... folder
sed -n 's/.*"\(.*\)"/\1/p'	#everything between " and " and after that
sed -n 's/.*\(.*\)"/\1/p'	#write everything from last " !OK
sed -n 's|.*kkk||p'			#everything from kkk
sed -n 's/kkx.*kky //p'		#deletes everyting between kkx kky
sed -n 's/kkk.*//p'			#deltees everyting from kkk
sed -n 's/.*://p'			#deletes everything up to :
sed '/parameters.dat/d'   	#deletes parameters.dat 
sed 's|.*"\(.*\)".*|\1|'	#grept alles zwischen erstem " und letztem "
sed -n 's|\(.*\).*|\1H|p'	#grep everything anhd append H
sed -n 's|\(.*\)/.*|\1|p'	#grep everything until last /
sed -n 's|.*/\(.*\)|\1|p'	#grep everything until last /
sed -n '6p' POSCAR			#print line 6
sed -n '6,$p' POSCAR			#print line 6 up to last line
sed -n '6,9p' POSCAR			#print line 6 up to line 9
sed -i '$d' file                        # delete last line from file 


sed 's|.*/\(.*_...-shift/POSCAR\).*|\1|'	#deletes everything up to expression
sed 1d						#deletes first line
sed 's|\(.*\)_disp:\(.*\)\/.*|\1_disp:\2|'   # ...._disp..../ searches for
_disp and after that / and prints everything

sed -n '/mach/p'				#sed acts like grep
sed -e 's|^.*\(.\)$|\1|'		#echo last character of a string
sed -e 's/^[ \t]*//'			#deletes leading whitespace
sed 's/[ \t]*$//'				#deletes trailing whitespace
sed 's|^[ \t]*||;s|[ \t]*$||'	#delete leading & trailing whitespace
sed 's/^[ \t]*//;s/[ \t]*$//'	#delete leading & trailing whitespace
sed -i 7d						#writes imediately to file
sed -n '$='						#counts lines = wc -l
sed -n '/word/='				#linenumber of word
sed -e :a -e '/^\n*$/{$d;N;ba' -e '}'	#delete all trailing blank lines at end of file
sed '2,4d' 					#deletes line 2 to 4
sed -n 's|\(^[ ]*[-+0-9.]*[ ]*[-+0-9.]*[ ]*[-+0-9.]*\).*|\1|p'	#plot first 3 numbers
sed -i '15  s|\(^[ ]*[-+0-9.]*[ ]*[-+0-9.]*[ ]*[-+0-9.]*\).*|\1|' POSCAR	#takes just care of line 15
sed -n 's|\([ ]*[+-]*[0-9]*[.]*[0-9]*\)\([ ]*[+-]*[0-9]*[.]*[0-9]*\)\([ ]*[+-]*[0-9]*[.]*[0-9]*\).*|a\1b\2c\3d|p'
sed 's|\t| |g'  ## change all tabs to whitespace

sw='\([ ]*[\t]*[+-]*[0-9]*[.]*[0-9]*\)'		#one coordinate in SPOSCAR, POSCAR..
sed -n '7 s|'"$sw$sw$sw"'.*|\1|p' SPOSCAR	#prints first coordinate in line7

########################
find . -name DISP.phon | sed 's|\(.*\)DISP.phon|mv & \1DISP|' | sh
sed 's/.*\(.$\)/\1/'			#letztes zeichen

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
sed '/^$/d'						#delete empty lines
sed '/^ *$/d'               #delete empty lines better


sed -i '1itexttext' file.dat	#adds texttext in line1 in file.dat
sed -i '$atexttext' file.dat	#append a line in file.dat with texttext
sed -e '1~2d' Fqh				# entfernt jede zweite zeile aus Fqh

sed 's| ||g'   entfernt alle leerzeichen


sed -n '1p' tmp   ## prints a certain line
sed -e's/  */ /g'   ### multiple blanks/spaces to one

--> lom=`getlinuxormac`;[ "$lom" = "Linux" ] && add="";[ "$lom" = "mac" ] && add="'' -e"
--> sed -i $add     instead of sed -i


--> try sed -i 's|sed -i|sed -i $add|g' file
--> try sed -i "s|wc -\([wlc]\)|wc -\1 \| sed 's\|[ ]\*\|\|g'|g" file

sed -i 's/\r//g' md_long_tox.c     #to remove ^M