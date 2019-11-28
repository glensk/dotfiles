#!/bin/bash

# following 3 lines must always be present
path=`set | grep BASH_SOURCE | sed 's|.*\"\(.*\)/[^/]*\".*|\1|' | sed '/BASH_SOURCE=/s/.*/\./'`
script=`set | grep BASH_SOURCE | sed 's|.*\".*/\([^/]*\)\".*|\1|' | sed '/BASH_SOURCE=/s/.*/\./'`
options=$*; . $path/utilities/functions.include; checkOptions "-h -help -i -d -l"

# small help for options
h=`getOption -h`
if [ $h == True ]; then
  usage $script
  printOptions "-i inputFile(s)  default: find -name 'OUTCAR*'" \
               "-d dirs          run find command only on dirs" \
               "-l               follow symbolic links, i.e., find -L -name ..."
  exit
fi

# detailed help
help=`getOption -help`
if [ $help == True ]; then
  details $script
  echo2 "   prints relevant information about convergence from OUTCAR"                     \
        "   giving warnings when convergence has not been reached"
  echo2 "   the output is:"                                                                \
        "     it-max:        nr of electronic iteration steps in last ionic step"          \
        "     NELM:          NELM value from OUTCAR (max. nr of allowed electronic steps)" \
        "     change-last:   energy change between the last and last but one electronic"   \
        "                    step within the last ionic step (energy-change tag);"         \
        "                    both the total energy and the eps (band energies) are"        \
        "                    considered, the larger one is taken similarly as in VASP"     \
        "     EDIFF:         EDIFF value from OUTCAR (energy convergence criterion)"       \
        "     status:        noIt:        if no Iteration tag can be found, i.e., job"     \
        "                                 is still running or dead" \
        "                    unfinished:  if Iteration tag can be found but not the"       \
        "                                 \"Total CPU time used (sec):\" tag"              \
        "                    FINISHED:    if \"Total CPU time used (sec):\" tag is found"  \
        "     convergence:   notConv:     if unfinished and energy-change>EDIFF"           \
        "                    conv:        if unfinished and energy-change<EDIFF"           \
        "                    ok:          if FINISHED   and energy-change<EDIFF"           \
        "                    UNCONV:      if FINISHED   and energy-change>EDIFF"           \
        "     last|LOOP:     time in hr since last write on OUTCAR and average time in"    \
        "                    hr for electronic iteration step (LOOP tag); these two times" \
        "                    allow to quickly check whether job is running or dead"
  echo2 "   a summary is written if the output is longer than 30 OUTCARs"
  echo2 "   IMPORTANT NOTE:  only the last ionic step is evaluated"                        \
        "                    for MD runs it might be necessary to check all ionic steps"
  exit
fi

# get files to process
inp=`getOption -i`
if [ $inp == True ]; then
  l=`getValue -i`;
  if [ -z "$l" ]; then error "-i option has no value"; fi
  check $l
else
  dOp=`getOption -d`
  if [ $dOp == True ]; then
    dd=`getValue -d`; if [ "$dd" == "" ]; then error "no value to -d option"; fi
  else
    dd="."
  fi
  link=`getOption -l`
  if [ $link == True ]; then link="-L"; else link=""; fi
  l=`find $link $dd -type f -name "OUTCAR*" -print`
  if [ -z "$l" ]; then error "no OUTCARs in this folder nor in the subfolders"; fi
fi

# get the longest path for nice printout
max=0
for i in $l; do
  max=`echo $i | wc -c | awk '{if ($1>'$max') print $1; else print '$max'}'`
done

# work through all OUTCARs
kpinfo=`getOption -k`; pathInfo=`getOption -p`; echo
fo=0; fu=0; uc=0; un=0; ni=0;
filename=`echo file | awk '{printf "%-'$max's",$1}'`
echo "# itMax NELM  changeLast  EDIFF      last|LOOP       status          $filename"
for i in $l; do
  # main work done here
  result=`zgrep ".*" $i | awk '
    BEGIN{it=-100;NELM=-100;change=-100;ediff=-100;status="unfinished ";s=0;c=0}
    /LOOP:/{s=s+$NF; c=c+1}
    /NELM   =/{split($3,a,";");NELM=a[1]}; /EDIFF  =/{EDIFF=$3};
    /Iteration/{split($0,a,"(");split(a[2],b,")");it=b[1]}
    /energy-change/{n=split($0,a,"(");split(a[n],b,")");change1=b[1];change1=(change1^2)^(1/2);
                    n=split($0,a,":");split(a[n],b," ");change2=b[1];change2=(change2^2)^(1/2);
                    if (change1>change2) change=change1; else change=change2}
    /Total CPU time used \(sec\):/{status="FINISHED"}
    END{if(c>0) printf "|%.1f ",s/c/3600; else printf "|---  ";
        print it,NELM,change,EDIFF,status,LOOP}'`
  
  # put data from results variable into nice format
  LOOP=`  echo $result | awk '{printf("%s",$1)}'`
  it=`    echo $result | awk '{printf("%4d",$2)}'`
  nelm=`  echo $result | awk '{printf("%4d",$3)}'`
  change=`echo $result | awk '{printf("%8.1e",$4)}'`
  ediff=` echo $result | awk '{printf("%8.1e",$5)}'`
  status=`echo $result | awk '{printf("%12s",$6)}'`

  # check if warning needs to be printed
  convergence=`echo $change $ediff $status | 
          awk '{if($1^2> $2^2&&$3=="unfinished")   printf "unConv";
                if($1^2<=$2^2&&$3=="unfinished")   printf "conv  ";
                if($1^2<=$2^2&&$3=="FINISHED")     printf "ok    ";
                if($1^2> $2^2&&$3=="FINISHED")     printf "UNCONV"}'`

  # increase counters
  if [ "$convergence" = "UNCONV" ]; then fu=`expr $fu + 1`; fi
  if [ "$convergence" = "ok    " ]; then fo=`expr $fo + 1`; fi
  if [ "$convergence" = "unConv" ]; then un=`expr $un + 1`; fi
  if [ "$convergence" = "conv  " ]; then uc=`expr $uc + 1`; fi

  # use color and bold if unconverged
  if [ "$convergence" = "UNCONV" ]; then convergence="\033[31m\033[1mUNCONV\033[0m"; fi

  # check if OUTCAR has no Iterations
  if [ $it == -100 ]; then
    status="        noIt";
    convergence="      ";
    it="----";
    change=" -------"
    ni=`expr $ni + 1`

    # reduce unConv/notConv because it was increased above also for OUTCARs with no iterations (and we want to treat both separately)
    un=`expr $un - 1`
  fi
  if [ $nelm == -100 ]; then nelm=" ---"; ediff=" -------"; fi

  # bring outar path to a consistent length
  outcar=`echo $i | awk '{printf "%s",$1}'`

  # get the time of last writing of OUTCAR which together with LOOP obtained above allows to have a quick check whether job still runs
  if [ $status == unfinished -o $status == noIt ]; then
    present=`date +%s`
    delta=`ls -tl --time-style=+%s $i | awk '{d=('$present'-$6)/3600; if (d<9999) printf "%6.1f",d; else printf "******"}'`
  else
    delta="      "
  fi

  echo -e "  $it $nelm    $change $ediff   $delta$LOOP  $status  $convergence  $outcar"
done

ll=`echo $l | awk '{print NF}'`
if [ $ll -gt 10 ]; then
  echo "# itMax NELM  changeLast  EDIFF      last|LOOP       status          $filename"
fi

# print summary if long output (more than 30 OUTCARs)
n=`echo $l | xargs -n1 | wc -l`
echo "#"
if [ $n -gt 30 ]; then
  echo "# ===========>>        noIt              $ni"
  echo "# ===========>>  unfinished notConv      $un"
  echo "# ===========>>  unfinished conv         $uc"
  echo "# ===========>>    FINISHED ok           $fo"
  if [ $fu != 0 ]; then
    echo -e "\033[1m\033[31m# ===========>>    FINISHED UNCONVERGED  $fu\033[0m"
  else
    echo -e "# ===========>>    FINISHED UNCONVERGED  $fu"
  fi
  echo "#"
fi
echo "# Note: Only the last ionic step of each OUTCAR is considered"


