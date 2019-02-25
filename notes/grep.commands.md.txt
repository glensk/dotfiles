grep -r -i ENAUG ./   				# grept in allen unterliegenden foldern, auch in gelinkten foldern
grep ENAUG INCAR parameters.dat # grept in ENAUG und INCAR
grep -o INCAR					# prints just INCAR as often as found
grep -o "/.x.x./"				# only prints everything matched . is every character
grep "energy  w\|energy w" file 		#greps for both terms
echo 90 | grep -Eo '^[0-9]+$'			#greps only integers	
grep -Eo '^[0-9]+$' | wc -l			#is number Integer?
grep -o "bohr\|ang"				#greps for bohr or ang
grep -m 3 POSITON $pfad				# stops after 3rd occurance of POSITION in $pfad

grep -P '^Al[\t ]' MeltingPoint*                #grept nur nach dem wort
grep -A 1000 "job_number:[ ]*384740" | grep -m 1 -B 1000 "==="   grep nach 384740 und bis zum naechsten ===

grep -n "hallo" file        # prints the line number(s) where hallo is found


grep 'hallo$'       will work
grep "hallo$"       wont work
