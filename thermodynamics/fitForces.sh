#!/bin/bash

#-----set parameters and paths------------------------
inBaseDef="forces."; orderDef="3"
#-----------------------------------------------------


# following 3 lines must always be present
path=`set | grep BASH_SOURCE | sed 's|.*\"\(.*\)/[^/]*\".*|\1|' | sed '/BASH_SOURCE=/s/.*/\./'`
script=`set | grep BASH_SOURCE | sed 's|.*\".*/\([^/]*\)\".*|\1|' | sed '/BASH_SOURCE=/s/.*/\./'`
options=$*; . $path/utilities/functions.include; checkOptions "-h -help -i -b -B -o -a -d -f"

# small help for options
# space between option and value must be 1
# space between option and Comment or value and Comment must be >2
h=`getOption -h`
if [ $h == True ]; then
  usage $script
  printOptions "-i inpFiles  input files to use (default: \${inBase}[0-9.]*)" \
               "-b inBase    input base of file name (default: $inBaseDef)" \
               "-B outBase   output base of file name (default: \${inBase}fitted.)" \
               "-o order     order of the fit (default: $orderDef)" \
               "-a aLat(s)   produce fitted forces at aLat(s) (default: all input aLat)" \
               "-d           export delta of fit for first atom in x (not compatible with -a)" \
               "-f           force fitting if other format"
  exit
fi

# detailed help
help=`getOption -help`
if [ $help == True ]; then
  details $script
  echo2 "   fits a set of forces (or other quantities) as a function of the lattice constant"
  echo2 "   the input files are expected to have the format (without the comments in brackets):" \
        "      F1x  F1y  F1z   (forces on first atom)" \
        "      F2x  F2y  F2z   (forces on second atom)" \
        "      ..."
  echo2 "   to fit other quantities with other formats use -f option"
  echo2 "   each input file name must contain the corresponding lattice constant (aLat)" \
        "      format: base.aLatSUFFIX    with default base=$inBaseDef" \
        "                                 with aLat format [0-9]*\.*[0-9]\+" \
        "                                 and arbitrary SUFFIX (e.g. angstrom)"
  echo2 "   by default the fitted forces are produced at the same lattice constants as" \
        "   contained in the input files; this can be changed with the -a option"
  echo2 "   maximum deviation of the fitted forces from the input forces is calculated" \
        "   if the input aLats are equal to the output aLats"
  exit
fi

# mathematica kernel if needed
checkAndSetMath

# get in and out base names
inBase=`getOption -b`; if [ $inBase == True ]; then inBase=`getValue -b`; else inBase=$inBaseDef; fi
ouBase=`getOption -B`; if [ $ouBase == True ]; then ouBase=`getValue -B`; else ouBase=${inBase}fitted.; fi

# get input files
inp=`getOption -i`;
if [ $inp == True ]; then
  inpFiles=`getValue -i`
  for i in $inpFiles; do check $i; done
else
  inpFiles=`ls ${inBase}[0-9.]*`
fi

# get fitting order and check if positive integer
order=`getOption -o`; if [ $order == True ]; then order=`getValue -o`; else order=$orderDef; fi
order=`echo "$order" | awk '{if ((int($1)^2)^(1/2)!=$1) print "error"; else print int($1)}'`
if [ "$order" == error ]; then error "given order $order empty or wrong (must be positive integer)"; fi

# check if sufficient nr of input files for given order
n=`echo $inpFiles | xargs -n1 | awk 'END{print NR}'`
if [ $n -le $order ]; then error "insufficient nr of input files ($n) for given order ($order)"; fi

# check if all input files have same nr of columns and lines
nprev=""
for i in $inpFiles; do
  n=`awk 'BEGIN{err=""};NR==1{nf=NF};nf!=NF{err="error"};END{if (err=="error") print err; else print NF,NR}' $i`
  if [ "$n" == "error" ]; then error "inconsistency in file $i; wrong nr of columns somewhere"; fi
  if [ "$n" != "$nprev" -a "$nprev" != "" ]; then error "files $i and $iprev have inconsistent nr of lines or columns"; fi
  nprev=$n; iprev=$i
done

# get number of columns and proceed only if ==3 or -f option given
ncolumn=`echo $n | awk '{print $1}'`
force=`getOption -f`
if [ $ncolumn != 3 -a $force != True ]; then error "nr columns in input files != 3; use -f option to force"; fi

# prepare the _tmp_allForces with all forces and corresponding _tmp_aLats
rm -f _tmp_aLats _tmp_allForces
for i in $inpFiles; do
  aLat=`echo $i | sed 's/'$inBase'\([0-9]*\.*[0-9]\+\).*/\1/'`
  echo $aLat >> _tmp_aLats
  cat $i | xargs >> _tmp_allForces
done

# check if -a option given
aLatOut=`getOption -a`
if [ $aLatOut == True ]; then
  aLatOut=`getValue -a`
  echo $aLatOut | xargs -n1 > _tmp_aLatsOut
else
  cp _tmp_aLats _tmp_aLatsOut
fi

# print output to log and stdout
log=log_fitForces
echo  > $log; echo " input files:" >> $log; echo $inpFiles >> $log
echo >> $log; echo " input aLats:" >> $log; cat _tmp_aLats | xargs >> $log
echo >> $log; echo " output aLats:" >> $log; cat _tmp_aLatsOut | xargs >> $log
echo >> $log; cat $log
echo -e " fitting order: $order" >> $log
echo -e " \033[1mfitting order: $order\033[0m"; echo; echo " ++++ fitting ++++"

# export delta fit?
delta=`getOption -d`

# run mathematica and do the fitting; fitted forces exported directly in correct format
rm -f _tmp_math
$math >> _tmp_math << EOF
aLats=Import["_tmp_aLats","Table"]//Flatten;
aLatsOut=Import["_tmp_aLatsOut","Table"]//Flatten;
forces=Import["_tmp_allForces","Table"]//Transpose;
fitPointsList={}; maxdiff=0;
Do[
  fit=Fit[{aLats,forces[[i]]}//Transpose,Table[x^j,{j,0,$order}],x];
  fitPoints=Table[fit/.x->aLatsOut[[j]],{j,aLatsOut//Length}];
  AppendTo[fitPointsList,fitPoints];
  If[aLats==aLatsOut,
    maxdiff=Max[Abs[forces[[i]]-fitPoints],maxdiff];
    fit = {aLats,forces[[i]]-fitPoints}//Transpose;
    (* remove i==1 if you need all deltas  *)
    If[i==1&&"$delta"=="True",Export["delta_fit_forces_"<>ToString[i],Append[fit,{Null,Null}],"Table"]];
  ];
,{i,forces//Length}]
If[aLats==aLatsOut,Print["maxdiff ",maxdiff]];
fitPointsList=Transpose[fitPointsList];
Do[
  forExport=Partition[fitPointsList[[j]],$ncolumn];
  Print[""];
  Export["$ouBase"<>ToString[aLatsOut[[j]]],forExport,"Table"];
  Print["exported: $ouBase"<>ToString[aLatsOut[[j]]]];
,{j,fitPointsList//Length}]
EOF

# print out produced output files
ouFiles=`grep "^exported: " _tmp_math | awk '{print $NF}'`
echo  >> $log; echo " output files:" >> $log; echo $ouFiles >> $log
tail -n 3 $log

# if aLats == aLatsOut then give maximum deviation of fitted values
maxdiff=`grep maxdiff _tmp_math | awk '{print $NF}'`
if [ -n "$maxdiff" ]; then
  echo; echo -e " max dev: \033[1m\033[31m$maxdiff\033[0m"
  echo >> $log; echo -e " max dev: $maxdiff" >> $log
fi
echo; echo >> $log
rm -f _tmp_math _tmp_allForces _tmp_aLats _tmp_aLatsOut

