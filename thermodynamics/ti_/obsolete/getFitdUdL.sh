#!/bin/sh

math=`which math`
if [ ! -e $math ]; then
		echo "" 1>&2
		echo -e "\033[31m\033[1mERROR\033[0m: $math not available" 1>&2
		exit 1
fi

if [ -e high_avg_dUdL ]; then
  echo; echo -e "\033[31m\033[1mWARNING\033[0m: high_avg_dUdL found, will be used for fitting!"
  avgFile=high_avg_dUdL
else
  if [ ! -e avg_dUdL ]; then
    echo "" 1>&2
    echo -e "\033[31m\033[1mERROR\033[0m: no avg_dUdL file" 1>&2
    exit 1
  fi
  avgFile=avg_dUdL
fi

if [ "$1" != "-p" ]; then
  if [ "$1" != "" ]; then
    file=$1
  else
    file=~/Thermodynamics/ti_/obsolete/parameters_fit_actinides
  fi
  if [ ! -e $file ]; then
    echo "" 1>&2
    echo -e "\033[31m\033[1mERROR\033[0m: no $file file" 1>&2
    exit 1
  fi
  ta0=`awk '{print $1}' $file`; ta1=`awk '{print $2}' $file`
  ta2=`awk '{print $3}' $file`; ta3=`awk '{print $4}' $file`
fi

linearAndCubicFit() {
$math >> tmp2 << EOF 
  avg=Import["tmp","Table"];
  lin=Fit[avg,{1,x},x];
  cub=Fit[avg,{1,x,x^2,x^3},x];
  {linfit,cubfit}=Append[Table[{i,#/.x->i},{i,0,1,0.01}],{Null,Null}]&/@{lin,cub};
  Export["fit_linear",linfit,"Table"]
  Export["fit_cubic",cubfit,"Table"]
  {lindev,cubdev}=Sum[Abs[(#/.x->avg[[i,1]])-avg[[i,2]]],{i,Length[avg]}]/Length[avg]&/@{lin,cub};
  {lin,cub}=Integrate[#,{x,0,1}]&/@{lin,cub};
  s=Round[100{lin,lindev}]/100//N;
  Export["Fah_linear",ToString[s[[1]]]<>" "<>ToString[s[[2]]],"String"]
  s=Round[100{cub,cubdev}]/100//N;
  Export["Fah_cubic",ToString[s[[1]]]<>" "<>ToString[s[[2]]],"String"]
EOF
}

tanFit() {
$math >> tmp2 << EOF 
  avg=Import["tmp","Table"];
  tanFunction[x_,a0_,a1_,a2_,a3_]:=a0*Cot[Pi (a1 x + a2)]+a3;
  (* old tan fit:  tanFunction[x_,a0_,a1_,a2_,a3_]:=-a0*Tan[Pi ((1-a1)x+a2+0.5)]+a3; *)
  tan=FindFit[avg,{tanFunction[x,a0,a1,a2,a3],
        0 < a0 < 200, 0. < a1 < 1, 0 < a2 < a1},
        {{a0,$ta0}, {a1,$ta1}, {a2,$ta2}, {a3,$ta3}},x];
  {b0,b1,b2,b3}={a0,a1,a2,a3}/.tan;
  Print["tangens ",b0," ",b1," ",b2," ",b3]
  tanfit=Append[Table[{i,tanFunction[i,a0,a1,a2,a3]/.tan},{i,0,1,0.01}],{Null,Null}];
  Export["fit_tangens",tanfit,"Table"];
  tandev=Sum[Abs[(tanFunction[avg[[i,1]],a0,a1,a2,a3]/.tan)-avg[[i,2]]],{i,Length[avg]}]/Length[avg];
  tanFah=Chop[NIntegrate[tanFunction[x,a0,a1,a2,a3]/.tan,{x,0,1}],10^(-3)];
  s=Round[100{tanFah,tandev}]/100//N;
  Export["Fah_tangens",ToString[s[[1]]]<>" "<>ToString[s[[2]]],"String"]
EOF
}

tanFit2() {
$math >> tmp2 << EOF 
  avg=Import["tmp","Table"];
  tanFunction[x_,a0_,a1_,a2_,a3_,a4_]:=a0*Cot[Pi (a4 x^2 + a1 x + a2)]+a3;
  (* old tan fit:  tanFunction[x_,a0_,a1_,a2_,a3_]:=-a0*Tan[Pi ((1-a1)x+a2+0.5)]+a3; *)
  tan=FindFit[avg,{tanFunction[x,a0,a1,a2,a3,a4],
        0 < a0 < 200, 0. < a1 < 1, 0 < a2 < a1, -1 < a4 <1},
        {{a0,$ta0}, {a1,$ta1}, {a2,$ta2}, {a3,$ta3}, {a4,0}},x];
  {b0,b1,b2,b3,b4}={a0,a1,a2,a3,a4}/.tan;
  Print["tangens ",b0," ",b1," ",b2," ",b3," ",b4]
  tanfit=Append[Table[{i,tanFunction[i,a0,a1,a2,a3,a4]/.tan},{i,0,1,0.01}],{Null,Null}];
  Export["fit_tangens2",tanfit,"Table"];
  tandev=Sum[Abs[(tanFunction[avg[[i,1]],a0,a1,a2,a3,a4]/.tan)-avg[[i,2]]],{i,Length[avg]}]/Length[avg];
  tanFah=Chop[NIntegrate[tanFunction[x,a0,a1,a2,a3,a4]/.tan,{x,0,1}],10^(-3)];
  s=Round[100{tanFah,tandev}]/100//N;
  Export["Fah_tangens",ToString[s[[1]]]<>" "<>ToString[s[[2]]],"String"]
EOF
}

# fermiFit() {
# $math >> tmp2 << EOF 
#   avg=Import["tmp","Table"];
#   fermiFunction[x_,a0_,a1_,a2_,a3_]:= If[x==0,
#                         a0*Log[(1+a1+a2)/a1-1]+a3,
#                         a0*Log[(1+a1+a2)/(a1+x)-1]+a3];
#   fermi=FindFit[avg,{fermiFunction[x,a0,a1,a2,a3],
#         0.0001<=a1,0.0001<=a2},
#         {{a0,$fa0}, {a1,$fa1}, {a2,$fa2}, {a3,$fa3}},x]
#   {b0,b1,b2,b3}={a0,a1,a2,a3}/.fermi;
#   Print["fermi ",b0," ",b1," ",b2," ",b3]
#   fermifit=Append[Table[{i,fermiFunction[i,a0,a1,a2,a3]/.fermi},{i,0,1,0.01}],{Null,Null}];
#   Export["fit_fermi",fermifit,"Table"];
#   fermiFah=Integrate[fermiFunction[x,a0,a1,a2,a3]/.fermi,{x,0,1}];
#   Export["Fah_fermi",Round[100fermiFah]/100//N,"Table"];
# EOF
# }

rm -f tmp2
sed '/^[ ]*#/d' $avgFile | awk '{print $1,$2+$3}' > tmp
rm -f Fah_linear
while [ ! -e Fah_linear ] ; do
  linearAndCubicFit
  errlin=`awk '{print $1}' Fah_linear`; errcub=`awk '{print $1}' Fah_cubic`
done
rm Fah_linear
sed '/^[ ]*#/d' $avgFile | awk '{print $1,$2}' > tmp
while [ ! -e Fah_linear ] ; do
  linearAndCubicFit
  errlin=`awk '{printf("%.2f        %.2f       %.2f",$1,(($1-1*('"$errlin"'))^2)^(1/2),$2)}' Fah_linear > tmp`; mv tmp Fah_linear
  errcub=`awk '{printf("%.2f        %.2f       %.2f",$1,(($1-1*('"$errcub"'))^2)^(1/2),$2)}' Fah_cubic > tmp`; mv tmp Fah_cubic
done

if [ -e Fah ]; then mv Fah Fah_old; string1=", Fah_old"; else string1=""; fi
echo " fit      Fah(meV/at) err/2(meV/at) dev(meV/at)" > Fah
echo " linear      `cat Fah_linear`" >> Fah
echo " cubic       `cat Fah_cubic`" >> Fah
rm Fah_linear Fah_cubic; echo; cat Fah
if [ "$1" != "-p" ]; then
  sed '/^[ ]*#/d' $avgFile | awk '{print $1,$2+$3}' > tmp
  rm -f Fah_tangens
  while [ ! -e Fah_tangens ]; do
    tanFit
    errtan=`awk '{print $1}' Fah_tangens`
  done
  sed '/^[ ]*#/d' $avgFile | awk '{print $1,$2}' > tmp
  rm Fah_tangens 
  while [ ! -e Fah_tangens ]; do
    tanFit
    errtan=`awk '{printf("%.2f        %.2f       %.2f",$1,(($1-1*('"$errtan"'))^2)^(1/2),$2)}' Fah_tangens > tmp`; mv tmp Fah_tangens
  done
  sed -n -e 's/.*tangens \(.*\)/\1/p' tmp2 > parameters_fit
  echo " tangens     `cat Fah_tangens`" >> Fah; tail -n1 Fah
  string2=", fit_tangens, parameters_fit"
  rm Fah_tangens

  # run now tanFit with one more, quadratic fitting parameter (tanFit2)
  # take as start values the ones from the previous fit (and a4=0)
  ta0=`awk '{print $1}' parameters_fit`; ta1=`awk '{print $2}' parameters_fit`
  ta2=`awk '{print $3}' parameters_fit`; ta3=`awk '{print $4}' parameters_fit`

  sed '/^[ ]*#/d' $avgFile | awk '{print $1,$2+$3}' > tmp
  rm -f Fah_tangens
  while [ ! -e Fah_tangens ]; do
    tanFit2
    errtan=`awk '{print $1}' Fah_tangens`
  done
  sed '/^[ ]*#/d' $avgFile | awk '{print $1,$2}' > tmp
  rm Fah_tangens 
  while [ ! -e Fah_tangens ]; do
    tanFit2
    errtan=`awk '{printf("%.2f        %.2f       %.2f",$1,(($1-1*('"$errtan"'))^2)^(1/2),$2)}' Fah_tangens > tmp`; mv tmp Fah_tangens
  done
  sed -n -e 's/.*tangens \(.*\)/\1/p' tmp2 > parameters_fit
  echo " tangens2    `cat Fah_tangens`" >> Fah; tail -n1 Fah
  string2=", fit_tangens, parameters_fit"
  rm Fah_tangens
fi

echo >> Fah; tail -n1 Fah
echo " l=0.5(cub)  `awk '$1==0.5{printf("%.2f",$2)}' fit_cubic`" >> Fah; tail -n1 Fah
if [ "$1" != "-p" ]; then
  echo " l=0.5(tan)  `awk '$1==0.5{printf("%.2f",$2)}' fit_tangens`" >> Fah; tail -n1 Fah
fi
echo; echo "Files written: Fah"$string1", fit_linear, fit_cubic$string2"
rm tmp2


