## if else:
cat mapping | awk '{for (i=1;i<=12;i++) printf "%s ",$i; printf "\n"}'  # schleife ueber $i

echo $number | awk '{printf "%.0f\n", $1}'  ### round number 9.45 -> 9; 9.7 -> 10
cat POSITIONs | awk '{if ($1<4) print $1; else print $1-7.5}'

echo "3.14592653589" | awk '{printf "%.2f\n", $1}' ## -> 3.15 == runden mit awk
echo 2.1 2 | awk '($1<=$2){print "yes"};($1>$2){print "no"}'    ## smaller greater with awk
cat file | awk '($3<=0.20){print $0}
cat file | awk '($3<=0.20 && $7>100){print $0}'

echo text | awk '{if(min==""){min=max=$1}; if($1>max) {max=$1}; if($1<min) {min=$1}; total+=$1; count+=1} END {print total/count, max,  min}'

NR                                              ## zeilenNummer
FNR                                             ## dateinummer wenn mehr als eine datei an awk uebergeben wird
print NF                                        ## anzahl spalten
print $NF                                       ## letzte spalte
awk 'END{print NR}'                             ## Das END ist ein Muster, das besagt, dass die Aktion erst ausgeführt werden soll, wenn alle Zeilen bearbeitet wurden.
awk '{sum+=$2} END{print sum}' datei.txt        ## summiert die zweite spalte
awk '{sum+=$NF} {print $1,sum}' free_energy     ## running sum
awk '{sum+=$NF} {print $1,sum/NR}' free_energy  ## running average

awk '{sum+=$3} {print sum/NR}' free_energy      ## running average einer bestimmten spalte
awk '{d3=$3-avg3;avg3+=d3/NR;m3+=d3*($3-avg3)}{print sqrt(m3/NR)}' $file | tail -1  ## standard deviation

awk 'NR%3==1' file ## print every 3rd line


paste file1 file2 can be used to merge files with diffrent length of rows

awk 'BEGIN{size=30} {mod=NR%size; if(NR<=size){count++}else{sum-=array[mod]};sum+=$1;array[mod]=$1;print sum/count}'  ## rolling average nimmt immer die 30 zeilen oder so


awk '{sum1+=$1;sum2+=$2;sum3+=$3} {printf "%.2f %.2f %.2f \n", sum1/NR,sum2/NR,sum3/NR}'    ## running average von 3 spalten
awk '$1 >= 200 {print $1,$2}' free_energy       ## print free energy ab 200K
awk '{print exp($1)}' test
int(x)
sqrt(x)
exp(x)
log(x)
sin(x)
cos(x)
atan2(y,x)
rand()

          

          ##FORMAT
          ##FORMAT
          %-4s is cool  schreibt (falls nicht am . ausgereichtete werte, kann so platz geschaffen werden dass alle spalten gleich gorss sind)

1    hallo
2500 hallo
 awk '{printf(" %25.15f %25.15f %25.15f\n",$1*'"$scALat"'/'"$sc"',$2*'"$scALat"'/'"$sc"',$3*'"$scALat"'/'"$sc"')}'   Buendies am . plotten von spalten
          awk '{ printf "%-8s %-8s %-8s %-8s %-8s %-8s\n", $1, $2, $3,$4,$5,$6,$7}'
          awk '{printf "%.6f %.6f %.6f\n", $1-$4, $2-$5, $3-$6}'


          min=`echo $allene | xargs -n1 | awk 'min==""|| $1<min {min=$1}END{ print min}'`
          max=`echo $allene | xargs -n1 | awk 'max==""|| $1>max {max=$1}END{ print max}'`

arr[5]=7; # setzt den 6. Wert des Arrays auf 7
arr[5,3]="hallo"; # setzt den 4. Wert des 6. Arrays auf "hallo". (Genauer gesagt setzt es "5 SUPSEP 3" auf "hallo")
arr["first"]=8; # assoziatives Array

function name(arg1, arg2){

...
}
awk 'program' inputfile1 inputfile2 ...
awk -f programfile inputfile1 inputfile2 ...




#convert to bohr
awk '{a=0.529177;print $1,$2*a,$3*a,$4*a}' g0 > g0.xyz

# prints out coords from pdb centering in a box of 13 A
awk '{s=6.5;printf "%10.5f %10.5f %10.5f\n",$6+s,$7+s,$8+s}' fullerene.pdb

# sum a list of values
grep "LOOP:  VPU " OUTCAR | awk '{print $7}'| awk 'BEGIN{S=0}{F+=$0}END{print F}'

#awk format print
awk '{printf "%20.9f  %20.9f\n",$1,$2}' Cp_Paper_Z.K.Liu_2010.dat  

awk 'BEGIN{s=0};{s=s+$1};END{print s}' forces.4.05

awk '{a=0.529177;print $1,$2*a,$3*a,$4*a}' g0 > g0.xyz
awk '{a=print $1," \t" $2," \t" $3}' ../2.cutoff_50Ry/GEOMETRY
awk '$1 >1 {print $1,$2; }'            (print nur wenn $1 >1 ist
awk '$1 ~ /0.81/ {print $1,$2; }'            (print nur wenn $1 hat bestimmtes Muster hier 0.81 wobei zwischen den zeichen alles moeglich noch sein darf
awk 'BEGIN {sum=0;} {sum=sum+$2;} END  {print sum,NR,sum/NR; }' file  (Summiert alle Zahlen in Spalte 2)
awk '  '{print $1}'   file                        (Print column)
awk '  '{print $2}'   file                       (like: cat file)
awk '  '{print $1" text  "$3}'   file        (Print column1 text column3)
awk /STEP/ file                                    like: grep STEP file
awk /7...2/ file                                    gibt alles aus das 7xxx2 findet
ak /imi|usi/ file                gibt alle imi und usi dateien aus
awk /\./ file                     sucht nach punkten
awk '/PLOTTED COLUMNS ARE :/, $NF ~ /BLOCKEND/ ' file.dat 	#plots every line (including) PLOTT... and BLOCKEND

echo -e "pfad \t $n1 \t $n2 \t $n3 \t $n4 \t $n5 \t $n6" | awk '{ printf "%-'"$la"'s %-8s %-8s %-8s %-8s %-8s %-8s\n", $1, $2, $3,$4,$5,$6,$7}'

awk '{for(i=1;i<=NF;i++) t+=$i; print t; t=0}'		#summiert eine zeile
awk '{for(i=2;i<=NF;i++) t+=$i; print t; t=0}'		#summiert zweite bis letzte spalte


## further formatted printing
paste forces_pos/forces_1-1 forces_pos/forces-background | awk '{printf "%.6f %.6f %.6f\n", $1-$4, $2-$5, $3-$6}'  #prints number instead of 1e-06


cat file | awk 'ORS=(FNR%1000)?FS:RS' 			prints everything in one line, seperated with " "

multsoll=`echo "$cellvolume_soll/$cellvolume_ist" | bc -l | awk '{printf "%10.15f",$1^(1/3)}'`   ## root

	     awk 'BEGIN{sum=0;c=0;}
	          NF==8{sum=sum+$8; c=c+1;
	          printf("%.11f\n",sum/c)}' > avgdUdL_l$i

      >>>>>>>>>> BETTER:
                  
	     awk 'BEGIN{sum=0;c=0;}
	          {sum=sum+$8; c=c+1;
	          printf("%.11f\n",sum/c)}' > avgdUdL_l$i

awk 'NR==1; END{print}' output_0.000GPa/free_energy  ## prints first and last line of file


### total cool!!! joins lines which belong together from file 1 file 2
awk 'FNR==NR{a[$1]=$2 $3;next}{ print $0, a[$1]}' f1 f2 | awk 'NF==3{print $0}'
besser ist noch::: join --nocheck-order -1 1 -2 1 $file1 $file2    (evtl. mit | cut -d" " -f1-3    wobei f1-3 die anzahl der spalten sind)


### was auch geht umd zusammengehoerige lines zusammenzubringen:
awk '$3=='$seed'&&$4=='$structure'{print $1}'


## krasses konstrukt
zgrep -a --text -B $zeilen "^  <energy>" $pfad | awk '{vecanteil1[NR]=$2;vecanteil2[NR]=$3;vecanteil3[NR]=$4};BEGIN{coff=0;nions='$nions';nionsm1='$nionsm1';\
a1=99999;a2=99999;a3=99999;\
b1=99999;b2=99999;b3=99999;\
c1=99999;c2=99999;c3=99999;\
d1=99999;d2=99999;d3=99999;\
e1=99999;e2=99999;e3=99999;\
f1=99999;f2=99999;f3=99999;\
};
         /<i name="e_fr_energy">/{coff++;{linefr=NR;free=$3;next}};
         NR==linefr+1 {ewe=$3};
         NR==linefr+2 {eS0=$3};
         NR==linefr+8 {a1=$2;a2=$3;a3=$4};
         NR==linefr+9 {b1=$2;b2=$3;b3=$4};
         NR==linefr+10 {c1=$2;c2=$3;c3=$4};
         NR==linefr+14 {d1=$2;d2=$3;d3=$4};
         NR==linefr+15 {e1=$2;e2=$3;e3=$4};
         NR==linefr+16 {f1=$2;f2=$3;f3=$4};

        coff>0 && NR>linefr+19+'"$nions"'+3 && NR<=linefr+19+'"$nions"'+3+'"$nions"' {print \
        coff,NR,linefr,\
        "|||realcoords>",\
        vecanteil1[NR-3-'"$nions"']*a1+vecanteil2[NR-3-'"$nions"']*b1+vecanteil3[NR-3-'"$nions"']*c1,\
        vecanteil1[NR-3-'"$nions"']*a2+vecanteil2[NR-3-'"$nions"']*b2+vecanteil3[NR-3-'"$nions"']*c2,\
        vecanteil1[NR-3-'"$nions"']*a3+vecanteil2[NR-3-'"$nions"']*b3+vecanteil3[NR-3-'"$nions"']*c3,\
        "|||recipcoord>",\
        vecanteil1[NR-3-'"$nions"'],\
        vecanteil2[NR-3-'"$nions"'],\
        vecanteil3[NR-3-'"$nions"'],\
        "|||forces>",$2,$3,$4,"|||realcell>",a1,a2,a3,b1,b2,b3,c1,c2,c3,"||recipcell>",d1,d2,d3,e1,e2,e3,f1,f2,f3,"|||free>",free,"|||ewe>",ewe,"|||eS0>",eS0}' | sed 's|</i>||g'


awk 'BEGIN {x=0; while(++x<=10){print x; }; exit}'  ## print a sequence with awk

awk '{for (i=2;i<=NF;i++)print $i}'
awk 'NR%2==1 {print }' Fqh ( to get every second line )
