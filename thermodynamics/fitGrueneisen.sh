#!/bin/bash

#-----set parameters and paths--------------------------------------------
inBaseDef="MeshFreqs_"; outBaseDef="Grueneisen"; strDef=fcc
bin=3; deltaBin=0.1; # size of the bin (bin) used for calculating
                     # the mean gamma and delta (deltaBin) used for moving
                     # the bin; both in meV if input phonons are in meV
#-------------------------------------------------------------------------


# following 3 lines must always be present
path=`set | grep BASH_SOURCE | sed 's|.*\"\(.*\)/[^/]*\".*|\1|' | sed '/BASH_SOURCE=/s/.*/\./'`
script=`set | grep BASH_SOURCE | sed 's|.*\".*/\([^/]*\)\".*|\1|' | sed '/BASH_SOURCE=/s/.*/\./'`
options=$*; . $path/utilities/functions.include; checkOptions "-h -help -i -b -B -s -g"

# small help for options
# space between option and value must be 1
# space between option and Comment or value and Comment must be >2
h=`getOption -h`
if [ $h == True ]; then
  usage $script
  printOptions "-i inpFiles  input files to use (default: \${inBase}[0-9.]*)" \
               "-b inBase    input base of file name (default: $inBaseDef)" \
               "-B outBase   output base name (default: $outBaseDef)"  \
               "-s str       structure (bcc, fcc, sc, vol; see -help)" \
               "-g           produce Grueneisen fitted phonons at input aLats"
  exit
fi

# detailed help
help=`getOption -help`
if [ $help == True ]; then
  details $script
  echo2 "   fits the Grueneisen function to a set of phonon frequencies"
  echo2 "   the input files are expected to have a consistent number of columns and rows"
  echo2 "   each input file name must contain the corresponding lattice constant (aLat)" \
        "      format: BASEaLatSUFFIX    with default BASE=$inBaseDef" \
        "                                with aLat format [0-9]*\.*[0-9]\+" \
        "                                and arbitrary SUFFIX (e.g. angstrom)"
  echo2 "   each aLat is converted to a volume according to structure which can be" \
        "   adjusted with '-s str' option (implemented: bcc,fcc,sc,vol; default: $strDef)" \
        "   for str=vol no transformtion is performed, i.e., aLat is assumed to be a volume"
  echo2 "   by default the Grueneisen fitted phonons are not written; this can be achieved" \
        "   by using the -g option (output file names: \${inBase}Grueneisen_aLat)"
  exit
fi

# mathematica kernel if needed
checkAndSetMath

# get in and out base names
inBase=`getOption -b`; if [ $inBase == True ]; then inBase=`getValue -b`; else inBase=$inBaseDef; fi
outBase=`getOption -B`; if [ $outBase == True ]; then outBase=`getValue -B`; else outBase=$outBaseDef; fi

# get input files
inp=`getOption -i`;
if [ $inp == True ]; then
  inpFiles=`getValue -i`
  for i in $inpFiles; do check $i; done
else
  inpFiles=`ls ${inBase}[0-9.]* 2> /dev/null`
  if [ "$inpFiles" == "" ]; then error "no input files: ${inBase}[0-9.]*"; fi
fi

# check if sufficient nr of input files
n=`echo $inpFiles | xargs -n1 | awk 'END{print NR}'`
if [ $n -le 1 ]; then error "insufficient nr of input files"; fi

# check if all input files have same nr of columns and lines
nprev=""
for i in $inpFiles; do
  n=`awk 'BEGIN{err=""};NR==1{nf=NF};nf!=NF{err="error"};END{if (err=="error") print err; else print NF,NR}' $i`
  if [ "$n" == "error" ]; then error "inconsistency in file $i; wrong nr of columns somewhere"; fi
  if [ "$n" != "$nprev" -a "$nprev" != "" ]; then error "files $i and $iprev have inconsistent nr of lines or columns"; fi
  nprev=$n; iprev=$i
done

# get structure and check
str=`getOption -s`
if [ $str == True ]; then str=`getValue -s`; else str=$strDef; fi
if [ "$str" != fcc -a "$str" != bcc -a "$str" != sc -a "$str" != vol ]; then error "str in -s option unknown"; fi

# prepare the _tmp_allPhonons with all forces and corresponding _tmp_aLats
rm -f _tmp_aLats _tmp_allPhonons
for i in $inpFiles; do
  aLat=`echo $i | sed 's/'$inBase'\([0-9]*\.*[0-9]\+\).*/\1/'`
  echo $aLat >> _tmp_aLats
  cat $i | xargs >> _tmp_allPhonons
done

# get -g option for writing out fitted phonons
gOp=`getOption -g`
if [ $gOp == True ]; then nCol=`echo $inpFiles | xargs -n1 | awk '{print $1}' | xargs awk 'END{print NF}'`; fi

# print output to log and stdout
log=log_fitGrueneisen
echo  > $log; echo " input files: "$inpFiles >> $log
echo >> $log; echo " aLats: `cat _tmp_aLats | xargs`" >> $log
echo >> $log; echo " structure:  $str" >> $log;
echo >> $log; cat $log
echo " ++++ fitting ++++"

# run mathematica and do the fitting; fitted forces exported directly in correct format
rm -f _tmp_math
$math >> _tmp_math << EOF
aLats=Import["_tmp_aLats","Table"]//Flatten;
Switch["$str",
  "fcc",vols=aLats^3/4,
  "bcc",vols=aLats^3/2,
  "sc",vols=aLats^3,
  "vol",vols=aLats,
  _,error["unknown str"];
];
phonons=Import["_tmp_allPhonons","Table"]//Transpose;
m=Mean/@phonons;
g=FindFit[{vols,#}//Transpose, a V^-gamma, {{a, 10000}, {gamma, 2.5}}, V] & /@ phonons // Quiet;
gList=Table[{m[[i]],gamma/.g[[i,2]]},{i,m//Length}];
gList=Sort[gList];
Export["$outBase",gList,"Table"];

If["$gOp"=="True",
  fitted=Table[(a vols[[i]]^-gamma/.g[[j]]),{i,vols//Length},{j,g//Length}];
  Do[Export["${inBase}Grueneisen_"<>ToString[aLats[[i]]],Partition[fitted[[i]],$nCol],"Table"],{i,fitted//Length}];
];

deltaList=Table[{m[[i]],Max[Abs[ phonons[[i]] - ((a #^-gamma/.g[[i]]) &/@vols) ]]},{i,m//Length}];
Export["${outBase}_delta",deltaList,"Table"];

min=Min[Transpose[gList][[1]]];
max=Max[Transpose[gList][[1]]];
m={};
Do[
  bin = Select[gList,(x <= #[[1]] <= x+$bin)&];
  mean = Mean[Transpose[bin][[2]]];
  AppendTo[m,{x,mean}];
,{x,min,max,$deltaBin}];
Export["${outBase}_mean",m,"Table"];

Print["Max delta:  ", Max[deltaList//Transpose//Last]];
Print["Mean delta: ", Mean[deltaList//Transpose//Last]];
Print["Mean gamma: ", Mean[gList//Transpose//Last]];

EOF

# print out produced output files
outFiles="log_fitGrueneisen ${outBase} ${outBase}_mean ${outBase}_delta"
echo " output files: $outFiles" >> $log
echo; tail -n 3 $log

# print out error and mean gamma
maxdelta=`grep "Max delta:" _tmp_math | awk '{printf "%.2f",$NF}'`
meandelta=`grep "Mean delta:" _tmp_math | awk '{printf "%.2f",$NF}'`
meangamma=`grep "Mean gamma:" _tmp_math | awk '{printf "%.2f",$NF}'`
echo  >> $log;
echo " Max  delta:      $maxdelta"  >> $log
echo " Mean delta:      $meandelta" >> $log
echo " Mean Grueneisen: $meangamma" >> $log
tail -n 4 $log

echo; echo >> $log
rm -f _tmp_math _tmp_allPhonons _tmp_aLats

