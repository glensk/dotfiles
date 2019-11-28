#!/bin/bash

#-----set parameters and paths--------------------------
VERYEMPTY=20    # when to give a note about empty bands
VERYFULL=0.001  # when to give stronger warning
#-------------------------------------------------------


# following 3 lines must always be present
path=`set | grep BASH_SOURCE | sed 's|.*\"\(.*\)/[^/]*\".*|\1|' | sed '/BASH_SOURCE=/s/.*/\./'`
script=`set | grep BASH_SOURCE | sed 's|.*\".*/\([^/]*\)\".*|\1|' | sed '/BASH_SOURCE=/s/.*/\./'`
options=$*; . $path/utilities/functions.include; checkOptions "-h -help -i -d"

# small help for options
# space between option and value must be 1
# space between option and Comment or value and Comment must be >2
h=`getOption -h`
if [ $h == True ]; then
  usage $script
  printOptions "-i inputFile(s)  default: find -name 'OUTCAR*'" \
               "-d               print more details (see -help)"
  exit
fi

# detailed help
help=`getOption -help`
if [ $help == True ]; then
  details $script
  echo2 "   prints relevant information about band occupation from OUTCAR" \
        "   giving warnings when the highest band is occupied"
  echo2 "   a stronger warning is given when highest band is higher than $VERYFULL"
  echo2 "   a note is also issued when there are more than $VERYEMPTY bands free"
  echo2 "   in normal mode (no -k option) the output is:" \
        "     irr-kp:        nr of irreducible k-points" \
        "     elec:          nr of electrons" \
        "     ionStepsOfAll: number of ionic steps which contain k-point information (1. column)" \
        "                    and number of all ionic steps (2. column)" \
        "     lowOcc:        nr of lowest occupied band" \
        "     highOcc:       nr of highest occupied band" \
        "     nbands:        nr of bands" \
        "     fullOfAll:     nr of full k-points (1. column)" \
        "                    and nr of all checked k-points (2. column)" \
        "     ERROR:         if applies a warning with the highest occupancy" \
        "     NOTE:          if applies a note about the number of empty bands"
  echo2 "   note that lowOcc, highOcc, fullOfAll, ERROR, and NOTE refer only to the" \
        "   k-points which contain explicit band occupancy information in the OUTCAR" \
        "   (the output of the latter is controlled with the NWRITE flag in vasp)"
  echo2 "   with the -k option, for lowOcc and highOcc the number of the respective" \
        "   k-point is additionally printed; if information about more ionic steps" \
        "   is available in the OUTCAR, the number of the respective ionic step is" \
        "   also printed"
  echo2 "   note if all k-points differ from irr-kp this is due to several ionic steps" \
        "   included in the counting of all k-points"
  exit
fi

# get files to process
inp=`getOption -i`
if [ $inp == True ]; then
  l=`getValue -i`;
  if [ -a "$i" ]; then error "-i option has no value"; fi
  check $i
else
  l=`find -L . -type f -name "OUTCAR*" -print`
  if [ -z "$l" ]; then error "no OUTCARs in this folder nor in the subfolders"; fi
fi

# get the longest path for nice printout
max=0
for i in $l; do
  max=`echo $i | wc -c | awk '{if ($1>'$max') print $1; else print '$max'}'`
done

# work through all OUTCARs
ok=0; empty=0; full=0; very=0; nokp=0; c=0; details=`getOption -d`; pathInfo=`getOption -p`; echo
filename=`echo file | awk '{printf "%'$max's",$1}'`
if [ $details == True ]; then
  echo "#$filename  irr-kp  elec  nband  ionStepsOfAll    lowOcc at kp        highOcc at kp      full of all   status"
else
  echo "#$filename  irr-kp  elec      lowOc highOc nband     status"
fi
for i in $l; do
  # main work done here
  result=`zgrep ".*" $i | awk '
    BEGIN{f1=0;f2=0;fullkp=0;highestValue=0;allkp=0;highOcc=0;lowOcc=10e5;nmd=0;nkp=0;n=0};  # f1=0(1): out of(in) k-point list, f2=0(1): out of(in) single k-point
    /irreducible/{irrkp=$2}; /NELECT/{elec=$3}; /NBANDS=/{nbands=$NF}; /energy  w/{n=n+1};   # irreducible k-points, nr of electrons, nr of bands, nr of all md steps
    /-----------/{f1=0};                                                                     # pattern signifies end of k-point list
    /potential at core/{f1=1;nmd=nmd+1;nkp=0;                                                # pattern signifies beginning of k-point list
                        if("'"$details"'"!="True") {highestValue=0;highOcc=0;lowOcc=10e5}};  # if running without details we consider only kp info from last ionic step
    /band No/{if (f1==1) {f2=1;nkp=nkp+1;next}}                                              # pattern signifies start of a single k-point
    NF!=3&&f2==1{f2=0; allkp=allkp+1;                                                        # here we are at end of single k-point (NF!=3)
      if (highest<=lowOcc) { lowOcc=highest;  lowOccKP=nkp;  lowOccKP2=nmd};                 # if new highest is lower than all previous store nr of band/kp/MD step
      if (highest>=highOcc){highOcc=highest; highOccKP=nkp; highOccKP2=nmd};                 # if new highest is higher than all previous store nr of band/kp/MD step
      if (highest==nbands) {fullkp=fullkp+1; if(value^2>highestValue^2) highestValue=value}} # if completely full store occupation number if larger than all previous
    f1==1&&f2==1{if ($3!=0) {highest=$1; value=$3}}                                          #    square is needed because Methfessel-Paxton gives negative occupations
    END{print irrkp,elec,nbands,lowOcc,lowOccKP,lowOccKP2,highOcc,                           # both flags on so we collect nr of band and occupation if >0
              highOccKP,highOccKP2,fullkp,allkp,highestValue,nmd,n,"ss'$details'ss"}'`

  # bring outar path to a consistent length
  outcar=`echo $i | awk '{printf "%'$max's",$1}'`

  # put data from results variable into nice format
  irrkp=`  echo $result | awk '{printf("%5d",$1)}'`
  nelec=`  echo $result | awk '{printf("%5d",$2)}'`
  nbands=` echo $result | awk '{printf("%4d",$3)}'`
  lowOcc=`   echo $result | awk '{printf("%5d",$4)}'`
  lowOccKP=` echo $result | awk '{printf("%4d",$5)}'`
  lowOccKP2=`echo $result | awk '$1!=$11{printf("%4d",$6)}; $1==$11{print "    "}'`
  highOcc=`    echo $result | awk '{printf("%5d",$7)}'`
  highOccKP=`  echo $result | awk '{printf("%4d",$8)}'`
  highOccKP2=` echo $result | awk '$1!=$11{printf("%4d",$9)}; $1==$11{print "    "}'`
  fullkp=` echo $result | awk '{printf("%5d",$10)}'`
  allkp=`  echo $result | awk '{printf("%5d",$11)}'`
  highest=`echo $result | awk '{printf("%8.5f",$12)}'`
  nmd=`echo $result | awk '{printf("%4d",$13)}'`
  nmdall=`echo $result | awk '{printf("%4d",$14)}'`

  # check if warning needs to be printed
  status=`echo $highest | awk '{x=($1^2)^(1/2)}; x==0{print "ok";exit}; x>0&&x<'$VERYFULL'{print "full"}; x>='$VERYFULL'{print "very"}'`
  case "$status" in
    ok)   err="";;
    full) err="\033[1mFULL\033[0m  $highest"; c=`expr $c + 1`; full=`expr $full + 1`;;
    very) err="\033[1mFULL\033[0m  $highest  \033[31m\033[1mVERY!!!\033[0m"; c=`expr $c + 1`; very=`expr $very + 1`;;
  esac

  # check if there are very empty bands
  veryempty=`expr $nbands - $highOcc`
  if [ $veryempty -gt $VERYEMPTY ]; then veryempty="\033[1m\033[32mEMPTY\033[0m \033[1m$veryempty\033[0m"; empty=`expr $empty + 1`; else veryempty=""; fi

  # if no k-point information in OUTCAR
  if [ $nmd == 0 ]; then echo " $outcar $irrkp  $nelec      ++++ no kp info ++++"; nokp=`expr $nokp + 1`;continue; fi

  # if no error and no empty then ok
  if [ "$err" == "" -a "$veryempty" == "" ]; then err=ok; ok=`expr $ok + 1`; fi

  # two outputs possible more (-k option) or less details
  if [ $details == True ]; then
    echo -e " $outcar $irrkp  $nelec   $nbands  | $nmd  $nmdall  | $lowOcc $lowOccKP  $lowOccKP2  | $highOcc $highOccKP  $highOccKP2  | $fullkp $allkp    $err$veryempty"
  else
    echo -e " $outcar $irrkp  $nelec      $lowOcc $highOcc  $nbands       $err$veryempty"
  fi
done

# add additional description at bottom if long ouput
ll=`echo $l | awk '{print NF}'`
if [ $ll -gt 10 ]; then
if [ $details == True ]; then
  echo "#$filename  irr-kp  elec  nband  ionStepsOfAll    lowOcc at kp        highOcc at kp      full of all   status"
else
  echo "#$filename  irr-kp  elec      lowOc highOc nband     status"
fi
fi

# print summary
echo
echo    "# ===========>>  nokpInfo  $nokp"
if [ $empty != 0 ]; then
  echo -e "\033[32m\033[1m# ===========>>  empty     $empty\033[0m"
else
  echo    "# ===========>>  empty     $empty"
fi
echo    "# ===========>>  ok        $ok"
if [ $full != 0 ]; then
  echo -e "\033[1m# ===========>>  full      $full\033[0m"
else
  echo -e "# ===========>>  full      $full"
fi
if [ $very != 0 ]; then
  echo -e "\033[31m\033[1m# ===========>>  veryFull  $very\033[0m"
else
  echo -e "# ===========>>  veryFull  $very"
fi
echo


