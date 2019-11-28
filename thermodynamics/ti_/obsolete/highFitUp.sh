#!/bin/sh

math=/opt/mathematica/8.0.4/math
if [ ! -e $math ]; then
  echo "" 1>&2
  echo -e "\033[31m\033[1mERROR\033[0m: $math not available" 1>&2
  exit 1
fi

if [ $# -lt 1 ]; then
  echo "" 1>&2
  echo -e "\033[31m\033[1mUSAGE\033[0m:" 1>&2
  echo -e " getUpFit.sh -a                (to run on all *Ang_*K folders)" 1>&2
  echo -e " getUpFit.sh *Ang*K_folders    (to run on given folders)" 1>&2
  exit 1
fi
if [ "$1" == "-a" ]; then l=`ls -1d *Ang_*K`; else l="$*"; fi

linearAndQuadFit() {
$math >>/dev/null << EOF 
  avg=Import["tmp","Table"];
  con=Fit[avg,{1},x];
  lin=Fit[avg,{1,x},x];
  qua=Fit[avg,{1,x,x^2},x];
  {confit,linfit,quafit}=Append[Table[{i,#/.x->i},{i,0,1,0.01}],{Null,Null}]&/@{con,lin,qua};
  Export["fit_const",confit,"Table"]
  Export["fit_linear",linfit,"Table"]
  Export["fit_quad",quafit,"Table"]
  {condev,lindev,quadev}=Sum[Abs[(#/.x->avg[[i,1]])-avg[[i,2]]],{i,Length[avg]}]/Length[avg]&/@{con,lin,qua};
  {con,lin,qua}=Integrate[#,{x,0,1}]&/@{con,lin,qua};
  s=Round[100{con,condev}]/100//N;
  Export["AvgEn_const",ToString[s[[1]]]<>" "<>ToString[s[[2]]],"String"]
  s=Round[100{lin,lindev}]/100//N;
  Export["AvgEn_linear",ToString[s[[1]]]<>" "<>ToString[s[[2]]],"String"]
  s=Round[100{qua,quadev}]/100//N;
  Export["AvgEn_quad",ToString[s[[1]]]<>" "<>ToString[s[[2]]],"String"]
EOF
}

dir=`pwd`
for i in $l; do
 if [ -d $i ]; then
  echo; echo $i; cd $i
  ll=`ls avg_diff__*`
  for ii in $ll; do
    ext=`echo $ii | sed 's|avg_diff__\(.*\)|\1|'`
    sed '/^[ ]*#/d' $ii | awk '{printf("%5.2f %10.2f\n", $2,$3)}' > tmp
    if [ `cat tmp | wc -l` == 1 ]; then
      awk '{for (i=0;i<1.1;i=i+0.01) print i,$1}' tmp > fit_const__$ext
      echo "# only one lambda value given" > AvgEn__$ext
      echo " fit      AvgEn(meV/at) err/2(meV/at)" >> AvgEn__$ext
      echo " const        `cat tmp`" >> AvgEn__$ext
      echo; echo -e "  $ext   \033[31m\033[1mWARNING\033[0m only 1 lambda value found --> const fit only"
      head -n2 AvgEn__$ext | tail -n1; echo -e " const        \033[1m`cat tmp`\033[0m"
    else
      sed '/^[ ]*#/d' $ii | awk '{print $1,$2+$3}' > tmp
      linearAndQuadFit
      conlin=`awk '{print $1}' AvgEn_const`; errlin=`awk '{print $1}' AvgEn_linear`; errcub=`awk '{print $1}' AvgEn_quad`
      sed '/^[ ]*#/d' $ii | awk '{print $1,$2}' > tmp
      linearAndQuadFit
      mv fit_const fit_const__$ext
      mv fit_linear fit_linear__$ext
      mv fit_quad fit_quad__$ext
      conlin=`awk '{printf("%.2f         %.2f        %.2f",$1,(($1-1*('"$conlin"'))^2)^(1/2),$2)}' AvgEn_const > tmp`; mv tmp AvgEn_const
      errlin=`awk '{printf("%.2f         %.2f        %.2f",$1,(($1-1*('"$errlin"'))^2)^(1/2),$2)}' AvgEn_linear > tmp`; mv tmp AvgEn_linear
      errcub=`awk '{printf("%.2f         %.2f        %.2f",$1,(($1-1*('"$errcub"'))^2)^(1/2),$2)}' AvgEn_quad > tmp`; mv tmp AvgEn_quad

      echo " fit      AvgEn(meV/at) err/2(meV/at) dev(meV/at)" > AvgEn__$ext
      echo " const        `cat AvgEn_const`" >> AvgEn__$ext
      echo " linear       `cat AvgEn_linear`" >> AvgEn__$ext
      echo " quad         `cat AvgEn_quad`" >> AvgEn__$ext
      echo; echo "  $ext"; head -n1 AvgEn__$ext
      echo -e " const        \033[1m`awk '{print $1}' AvgEn_const`\033[0m    `awk '{printf("%10.2f %10.2f",$2,$3)}' AvgEn_const`"
      echo -e " linear       \033[1m`awk '{print $1}' AvgEn_linear`\033[0m    `awk '{printf("%10.2f %10.2f",$2,$3)}' AvgEn_linear`"
      echo -e " quad         \033[1m`awk '{print $1}' AvgEn_quad  `\033[0m    `awk '{printf("%10.2f %10.2f",$2,$3)}' AvgEn_quad`"
      rm AvgEn_const AvgEn_linear AvgEn_quad
    fi
  done
  cd $dir
 else
  echo; echo no directory $i; echo
 fi
done

