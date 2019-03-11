http://en.wikipedia.org/wiki/Sort_%28Unix%29

sort -g         # sorts alseo with exponential number 1.45342e-12
sort +1 -2		#sortiert nach der 2ten spalte (+1) und hoert vor der dirtten spalte wieder auf (-2) um dann nach der ersten ....spalte zu sortieren
#######
!!!!!!!!!!!!!!!!

sort +3 -4   sortiert nach der vierten !!!!!!!!! also das - sorgt fuer das sortieren...??


sort +1 -2  +3 -4 	#sortiert erst nach der 2ten und dann nach der 4ten spalte

sort -n +2 -3 +0 -2     # sortiert nach der 3ten spalte  und dann nach der 1ten spalte numerisch

sort:
text-1
text-333
text-1111
text-22

cat sort | sort -t - -k2 -n 		-->> defines how separator looks -t

text-1
text-22
text-333
text-1111

sort -n -k4 g0.xyz
ls -1d $lowfolder/*Ang_*K/ | sed 's|/$||' | sed 's|.*/||' | sed 's|\([0-9.]*\)Ang_\([0-9]*\)K.*|\2 \1Ang_\2K|' | sort -n | awk '{print $2}'
