#!/bin/bash

#-----set parameters and paths------------------------
# parameters should be ok in most cases
fft="default"    # fft mesh for xc potential; "default" takes the value provided by dstart; should be ok in most cases
lmax=12          # maximum l value for partial waves used inside atomic spheres
nslmax=6         # maximum l value for partial waves used in the computation of non- muffin-tin matrix elements
R0="0.00000500"  # step in radial mesh; from wien2k/SRC_structeditor/SRC_structgen/module.f
NPT=781          # number of points for radial mesh; from wien2k/SRC_structeditor/SRC_structgen/module.f
emax=4.0         # (Ry) the maximum energy for the window in which eigenvalues are calculated; we increase with respect
                 # to the standard 2 Ry to have a better description in the case of band structure calculations
                 # this increases CPU time but timing is not really an issue for 'standard' Murnaghan runs
if [ -e /data/grabowski/ ]; then wienPath=/data/grabowski/wien2k_src/; else wienPath=/home/grabowski/wien2k_src/; fi
#-----------------------------------------------------


toRy=0.0734986508052261761 # eV to Rydberg
toBohr=1.88972613289636593 # Angstrom to Bohr

# following 3 lines must always be present
path=`set | grep BASH_SOURCE | sed 's|.*\"\(.*\)/[^/]*\".*|\1|' | sed '/BASH_SOURCE=/s/.*/\./'`
script=`set | grep BASH_SOURCE | sed 's|.*\".*/\([^/]*\)\".*|\1|' | sed '/BASH_SOURCE=/s/.*/\./'`
options=$*; . $path/utilities/functions.include; checkOptions "-h -help -l -rerun -o -f -p -r -R"

# small help for options
h=`getOption -h`
if [ $h == True ]; then
  usage $script
  printOptions "-l nr=en                use \"en\" at \"nr\" as energy in in.in1"    \
               "-rerun l old new1 new2  rerun with optimized linearization energies" \
               "-o optdir               run with optimized case.in1 from optdir"     \
               "-R R=l1 a=aLat1 aLat2   run Rmt convergence for specified values"    \
               "-f                      force folder creation even if existing"      \
               "-p                      create parameters.dat template and exit"     \
               "-r                      automatically determine rmt and exit"
  exit
fi

# detailed help
help=`getOption -help`
if [ $help == True ]; then
  details $script
  echo2 "   creates input folders/files (in.in0,in.in1,...) for a Murnaghan calculation with wien2k and stores the"       \
        "   corresponding folders in jobList file"
  echo2 "   parameters.dat must be present (generate with -p and adjust)"
  echo2 "   check header of this script ($script) for some defaults (fft, lmax, nslmax, R0, NPT)"
  echo2 "   virtual crystal approximation (VCA) runs are possible by adjusting addCharge in parameters.dat"
  echo2 "   the -r option determines the Rmt using the nn program from wien2k (no parameters.dat needed); for accuracy"   \
        "   purposes it is however better to determine an optimum Rmt using the '-R R=l1 RK=l2 a=aLat' discussed below"
  echo2 "   if '-l nr=en' option is given then the energy of the nr'th exception in in.in1 is changed to the value of en"
  echo2 "   the '-l nr=en' option is helpful when a QTL-B Error occurs and the calculation stops adjusting the starting"  \
        "   values of the linearization energies in in.in1 can then help"
  echo2 "   several exceptions can be given, e.g.: $script -l 2=1 3=.5 6=.6"
  echo2 "   in some cases the '-l nr=en' option will not be sufficient to fully remove all QTL-B warnings (e.g., fcc Cu)" \
        "   then partial DOS calculations are necessary to determine optimal energies and an additional orbital"          \
        "   (typically for d channel); run such an optimization for a single volume of the curve in some separate"        \
        "   folder; instructions can be found in: utilities/wien2k_info"
  echo2 "   having determined the two optimal linearization energies for the problematic channel run"                     \
        "   '$script -rerun l old new1 new2' where l is the channel, old its linearization energy"                        \
        "   before optimization and new1 and new2 the two linearization energies obtained from the optimization"          \
        "   procedure; the energies are correctly scaled for the different volumes"
  echo2 "   NOTE: the rerun calculations must be started without the -in1new flag"
  echo2 "   if warnings still persist for some volumes, you can repeat the rerun calculation with readjusted new1/new2"
  echo2 "   after optimization, the optimized case.in1 files can be reused for instance for the case of convergence"      \
        "   calculations; use then the '-o optdir' option where optdir is the folder which contains the aLat folders of"  \
        "   the optimized case.in1 files"
  echo2 "   in several cases the Rmt value suggested by wien2k will not be the optimal one; to determine the optimal"     \
        "   value use the '-R R=l1 RK=l2 a=aLat' option, where l1 is the list of Rmt values to be calculated, l2 the"     \
        "   list of RmtKmax to be scanned and aLat the lattice constant; the minimum with respect to Rmt and RmtKmax"     \
        "   gives the optimum Rmt and RmtKmax; see also utilities/wien2k_info"
  echo2 "   typical script sequence:"        \
        "     \033[1mcreateFolders_murn_wien.sh\033[0m -p" \
        "     \033[1mcreateFolders_murn_wien.sh\033[0m -r" \
        "     \033[1mcreateFolders_murn_wien.sh\033[0m"    \
        "     submit.sh -r run.wien..."      \
        "     extractEnergies_wien.sh"       \
        "     getE0KFit.sh"                  \
        "     getThermodynamics.sh"
  exit
fi

# check if force overwriting templates or folders
overwrite=`getOption -f`

# check if we only generate a parameters.dat file; if so generate and exit
genPar=`getOption -p`
if [ $genPar == True ]; then
  if [ -e parameters.dat -a $overwrite != True ]; then error "parameters.dat existing; use -f to overwrite"; fi
  echo "
aLats=a1 a2 a3 a4 ...   # Angstrom

elem=XX         # Al, Ca, Cr, etc.
spin=XX         # FM or NM
addCharge=0.00  # add valence and nuclei charge (0..1 and %.2f) to do VCA (virtual crystal approx.)
xc=XX           # LDA or PBE
type=XX         # fcc or bcc or sc
semi=XX         # Rydberg (separation energy to core states; typically -6 or -7)
                
rmt=XX          # Bohr (muffin tin radius, typically 2.2 to 2.5, or use -r option for auto detect)
rmtkmax=XX      # 9..12
gmax=XX         # 12..18
                
kp=XX           
sh=XX           # true or false; if true: kp shift=(0.5 0.5 0.5)
                
smear=XX        # Fermi (native wien2k settings also possible (e.g., TETRA) see manual)
sigma=XX        # eV
" > parameters.dat
  echo; echo "parameters.dat written"; exit
fi

writeLog() {
  if [ ! -e log_createFolders_murn ]; then
    echo "" | awk '{printf("%16s %45s %72s\n","date","options","createdFolder")}' > log_createFolders_murn;
  fi
  if [ $rmtauto == True ]; then
    folder="---";
  else
    if [ $Ropt == True ]; then
      folder="RmtConvergence_${rmtkmax}RmtKmax_Gmax=$gmax${fftString}_${kp}x${kp}x${kp}kp${str}_$smear${sigma}eV"
    else
      folder=${rmtkmax}RmtKmax_Rmt=${rmt}_Gmax=$gmax${fftString}_${kp}x${kp}x${kp}kp${str}_$smear${sigma}eV
    fi
  fi
  options=`echo $* | awk '{printf("           %-60s           ",$0)}'`
  line="`date` $options $folder"
  echo "$line" >> log_createFolders_murn 
}

# get first part of parameters from parameters.dat needed even if we do rmt detection and exit
check parameters.dat
elem=`get elem`; aLats=`get aLats`; type=`get type`
checkInput "$elem" "$aLats" "$type" EOF

# get atomic number from database
file=$path/utilities/atomic_weights_and_ionic_compositions_NIST
check $file
Z=`awk 'BEGIN{Z=0};
        /Atomic Number =/{Ztmp=$NF};
        /Atomic Symbol = '"$elem"'$/{if (Z!=0&&Z!=Ztmp) Z=-1; else Z=Ztmp};
        END{print Z}' $file`
if [ $Z == 0 ];  then error "element $elem not existing in database"; fi
if [ $Z == -1 ]; then error "isotopes have different Atomic Number"; fi

# check if given strucuture type ok (fcc or bcc or sc so far)
case "$type" in
  fcc) type="F  ";;
  bcc) type="B  ";;
  sc) type="P  ";;
  *) error "structure type $type not known";;
esac

# prepare formatted input for in.struct
elemStr=`echo $elem | awk '{printf("%-3s",$1)}'`
Zstr=`echo $Z | awk '{printf("%-7.2f",'$addCharge'+$1)}'`
nptStr=`echo $NPT | awk '{printf("%4d",$1)}'`
R0str=`echo $R0 | awk '{printf("%10.8f",$1)}'`

# check if automatic determination of rmt; if so determine and exit
rmtauto=`getOption -r`
if [ $rmtauto == True ]; then
  dir=`pwd`; check $wienPath/setrmt_lapw;
  export PATH=$PATH:$wienPath # wien path needed inside setrmt_lapw
  echo; echo -e "\033[31m\033[1mAutomatic determination of RMT\033[0m"
  rm -fr _tmp_in; mkdir -p _tmp_/in; cd _tmp_/in

  # get smallest aLat for rmt determination
  amin=`echo $aLats | xargs -n1 | awk 'BEGIN{amin=10e5};$1<amin{amin=$1};END{print amin}'`
  astr=`echo $amin | awk '{printf("%9.6f",$1*'$toBohr')}'`
  echo "Murn
$type LATTICE,NONEQUIV.ATOMS:  1
MODE OF CALC=RELA unit=bohr
 $astr $astr $astr 90.000000 90.000000 90.000000
ATOM   1: X=0.00000000 Y=0.00000000 Z=0.00000000
          MULT= 1          ISPLIT= 2
$elemStr        NPT= $nptStr  R0=$R0str RMT=   2.00     Z: $Zstr
LOCAL ROT MATRIX:    1.0000000 0.0000000 0.0000000
                     0.0000000 1.0000000 0.0000000
                     0.0000000 0.0000000 1.0000000
                     0" > in.struct                     
  $wienPath/setrmt_lapw in
  rmt=`awk 'NR==7{printf("%.2f",$(NF-2))}' in.struct_setrmt`
  cd $dir; rm -fr _tmp_; echo

  # print out determined rmt and exit
  echo; echo;
  echo "----------------------------------------------------------------------------------"
  echo " RMT = $rmt Bohr determined by setrmt using smallest lattice constant ($amin Ang)"
  echo "----------------------------------------------------------------------------------"
  writeLog $*
  exit
fi

# read in remaining parameters from parameters.dat
rmtkmax=`get rmtkmax`; kp=`get kp`; sh=`get sh`; xc=`get xc`; addCharge=`get addCharge`
gmax=`get gmax`; smear=`get smear`; sigma=`get sigma`; rmt=`get rmt`; semi=`get semi`; spin=`get spin`
checkInput "$rmtkmax" "$kp" "$sh" "$xc" "$gmax" "$smear" "$sigma" "$rmt" "$semi" "$addCharge" "$spin" EOF
c=`checkReal $rmt`; if [ "$c" != ok ]; then error "rmt value in parameters.dat wrong"; fi

# check if given xc ok
case "$xc" in
  LDA) xc=5;;
  PBE) xc=13;;
  *)   error "xc $xc not known";;
esac

# k-points and shift
kpAll=`echo $kp | awk '{print $1^3}'`
if [ $sh == true ]; then str="-shift"; sh=1; else sh=0; fi

# check smearing method
if [ "$smear" != TETRA -a "$smear" != GAUSS -a "$smear" != TEMP -a "$smear" != TEMPS -a "$smear" != Fermi ]; then
  error "smearing method $smear not supported";
fi

# check fft parameter
if [ $fft == default ]; then fftString=""; else fftString="_fft=$fft"; fi

# check additional charge for VCA
c=`echo $addCharge | awk '$1>=-1&&$1<=1&&int($1*100)-$1*100==0 {print "ok"}'`
if [ "$c" != ok ]; then error "addCharge in parameters.dat is incorrect (should be 0<=addCharge<=1 and %.2f)"; fi

# check spin parameter
if [ "$spin" != FM -a "$spin" != NM ]; then
  error "spin type not supported (use FM or NM)"
fi

# check if we are doing a -rerun calculation
rerun=`getOption -rerun`
if [ $rerun == True ]; then
  # get rerun parameters
  rerun=`getValue -rerun`
  channel=`echo $rerun | awk '{print $1}'`; old=`echo $rerun | awk '{print $2}'`
  new1=`echo $rerun | awk '{print $3}'`; new2=`echo $rerun | awk '{print $4}'`
  if [ -z "$channel" -o -z "$old" -o -z "$new1" -o -z "$new2" ]; then
    error "parameters for -rerun option missing; check -help"; exit
  fi

  # deltas will be added on top of old energies; this will ensure appropriate scaling with volume
  delta1=`echo $old $new1 | awk '{printf("%.3f",$2-1*($1))}'`
  delta2=`echo $old $new2 | awk '{printf("%.3f",$2-1*($1))}'`

  # run over all lattice constants
  dir=`pwd`; rm -f jobList; echo ""
  for a in $aLats; do
    f="$rmtkmax"\RmtKmax_Rmt=$rmt\_Gmax=$gmax${fftString}_$kp\x$kp\x$kp\kp$str\_$smear${sigma}eV/$a\Ang/in
    echo $f
    if [ ! -d $f ]; then error "-rerun calculation but no folder $f"; fi
    if [ ! -e $f/in.in1 ]; then error "-rerun calculation but no previous $f/in.in1 file"; fi

    # check if -rerun is run for second time
    n=`awk 'BEGIN{s=0}; NR>3&&$1=='$channel'{s=s+1}; END{print s}' $f/in.in1`

    cd $f
    if [ $n != 1 ]; then 
      if [ -e in.in1_before_rerun ]; then
        echo "WARN: using in.in1_before_rerun"; echo
        mv in.in1_before_rerun in.in1
      else
        error "nr of orbitals for channel $channel different from 1 and no in.in1_before_rerun";
      fi
    fi
    awk 'NR<3{print}
         NR==3{printf("%7.5f %3d %3d       global e-param with N other choices, napw\n",$1,$2+1,$3)}
         NR>3&&$1!='$channel'{print}
         NR>3&&$1=='$channel'{printf("%2d %8.3f  %8.3f %s %d \n",$1,$2+1*('$delta1'),$3,$4,$5);
                              printf("%2d %8.3f  %8.3f %s %d \n",$1,$2+1*('$delta2'),$3,$4,$5)}' in.in1 > tmp
    mv in.in1 in.in1_before_rerun; mv tmp in.in1
    rm -f *broy*
    cd $dir; echo `pwd`/$f >> jobList
  done
  echo; echo "  RERUN calculation created (jobList)"

  writeLog $*
  exit
fi

# check if -R option given to run an Rmt optimization instead of murn
Ropt=`getOption -R`
if [ $Ropt == True ]; then
  Rvalues=`getValue -R`
  
  # option has several values: '-R R=l1 RK=l2 a=aLat'; split them here
  aValues=`echo $Rvalues | sed 's/R=\(.*\)a=\(.*\)/\2/'`
  Rvalues=`echo $Rvalues | sed 's/R=\(.*\)a=\(.*\)/\1/'`

  # now merge the Rmt and RKvalues into the single 1d array 'aLats' for easier implementation below
  aLats=""; for i in $Rvalues; do for j in $aValues; do aLats="$aLats Rmt=${i}_${j}Ang"; done; done
fi

# check if folders exist and should be overwritten
for a in $aLats; do
  if [ $Ropt == True ]; then
    # in this case aLats actually contains Rmt and RKvalues (see above)
    f="RmtConvergence_${rmtkmax}RmtKmax_Gmax=$gmax${fftString}_${kp}x${kp}x${kp}kp${str}_$smear${sigma}eV/$a/in"
  else
    f="$rmtkmax"\RmtKmax_Rmt=$rmt\_Gmax=$gmax${fftString}_$kp\x$kp\x$kp\kp$str\_$smear${sigma}eV/$a\Ang/in
  fi
  if [ -e $f -a $overwrite != True ]; then
    error "folder '$f' exists
       run './createFolders.sh -f' to force overwriting
       or with -rerun option for optimized calculation";
  fi
done

# create struct file template; aLat is left undefined here and set in later loop
elemStr=`echo $elem | awk '{printf("%-3s",$1)}'`
Zstr=`echo $Z | awk '{printf("%-7.2f",'$addCharge'+$1)}'`
rmtStr=`echo $rmt | awk '{printf("%5.2f",$1)}'`
nptStr=`echo $NPT | awk '{printf("%4d",$1)}'`
R0str=`echo $R0 | awk '{printf("%10.8f",$1)}'`
echo "Murn
$type LATTICE,NONEQUIV.ATOMS:  1
MODE OF CALC=RELA unit=bohr
 XXXXXXXXX XXXXXXXXX XXXXXXXXX 90.000000 90.000000 90.000000
ATOM   1: X=0.00000000 Y=0.00000000 Z=0.00000000
          MULT= 1          ISPLIT= 2
$elemStr        NPT= $nptStr  R0=$R0str RMT=  $rmtStr      Z: $Zstr
LOCAL ROT MATRIX:    1.0000000 0.0000000 0.0000000
                     0.0000000 1.0000000 0.0000000
                     0.0000000 0.0000000 1.0000000
                     0" > template.struct

# check if we are using optimized in.in1
oop=`getOption -o`;
if [ $oop == True ]; then
  if [ $Ropt == True ]; then error "-R and -o option not compatible"; fi
  optdir=`getValue -o`; if [ -z "$optdir" ]; then error "no value to -o option"; fi
fi

# now create the folders looping over the lattice constant
dir=`pwd`; rm -f jobList; check $wienPath/x; check $wienPath/instgen_lapw
for a in $aLats; do

  # folder name
  if [ $Ropt == True ]; then
    # in this case aLats actually contains Rmt and RKvalues (see above)
    f="RmtConvergence_${rmtkmax}RmtKmax_Gmax=$gmax${fftString}_${kp}x${kp}x${kp}kp${str}_$smear${sigma}eV/$a/in"
  else
    f="$rmtkmax"\RmtKmax_Rmt=$rmt\_Gmax=$gmax${fftString}_$kp\x$kp\x$kp\kp$str\_$smear${sigma}eV/$a\Ang/in
  fi
  echo; echo -e "\033[1m\033[31m$f\033[0m"; rm -fr $f; mkdir -p $f; cd $f

  # lattice constant to in.struct
  if [ $Ropt == True ]; then
    # if we have -R option aLat contains the lattice parameter and a the Rmt value which needs to be changed here
    astr=`echo $a | sed 's/Rmt=\(.*\)_\(.*\)Ang/\2/'`;
    c=`checkReal $astr`; if [ "$c" != ok ]; then error "parameters in '-R R=l1 a=aLat1 aLat2' option not correct (a=... probably)"; fi
    astr=`echo $astr | awk '{printf("%9.6f",$1*'$toBohr')}'`
    rmtStr=`echo $a | sed 's/Rmt=\(.*\)_\(.*\)Ang/\1/'`;
    c=`checkReal $astr`; if [ "$c" != ok ]; then error "parameters in '-R R=l1 a=aLat1 aLat2' option not correct (R=... probably)"; fi
    rmtStr=`echo $rmtStr | awk '{printf("%5.2f",$1)}'`
    sed -e 's/XXXXXXXXX/'"$astr"'/g' \
        -e 's/RMT=.*Z:/RMT=  '"$rmtStr"'      Z:/' $dir/template.struct > in.struct
  else
    astr=`echo $a | awk '{printf("%9.6f",$1*'$toBohr')}'`
    sed 's/XXXXXXXXX/'"$astr"'/g' $dir/template.struct > in.struct
  fi

  # create in.inst file
  $wienPath/instgen_lapw

  # check rmt
  $wienPath/x nn > tmp_log_nn << EOF
2
EOF
  if [ "`grep ERROR tmp_log_nn`" == "   ERROR !!!!!!!!!!!!!!!" ]; then
    error "RMT is too large!";
  fi

  # create symmetries for in.struct
  $wienPath/x sgroup
  mv in.struct_sgroup in.struct

  # renew Zstr because sgroup cuts digits after comma
  sed -i '7s/\(^.\{56\}\).*/\1'"$Zstr"'/' in.struct

  # run symmetry and lstart and prepare input files
  $wienPath/x symmetry
  $wienPath/x lstart > tmp_log_lstart 2> tmp_error << EOF
$xc
$semi
EOF
  if [ "`cat tmp_error`" != "LSTART ENDS" ]; then error "error in lstart (check semi parameter)"; fi
  cp in.in0_st in.in0; cp in.inc_st in.inc; cp in.inm_st in.inm
  cat in.in2_ls in.in2_sy > in.in2

  # change rmtkmax, lmax, nslmax, emax to provided values in in.in1
  in1Str=`echo $rmtkmax $lmax $nslmax | awk '{printf("%6.2f       %2d    %1d ",$1,$2,$3)}'`
  sed -e 's/.*\((R-MT\*K-MAX;.*\)/'"$in1Str"'\1/' \
      -e 's/\(K-VECTORS FROM UNIT:[0-9]\+[ ]\+[-0-9.]\+[ ]\+\)[-0-9.]\+\(.*\)/\1'"$emax"'\2/' in.in1_st > in.in1

  # check if manual setting of linearisation energies; if so change corresponding line in in.in1 for all given values
  lop=`getOption -l`;
  if [ $lop == True ]; then
    le=`getValue -l`
    for i in $le; do
      line=`echo $i | sed 's/\(.*\)=\(.*\)/\1/'`
      energy=`echo $i | sed 's/\(.*\)=\(.*\)/\2/'`
      awk 'NR!='$line'+3{print};
           NR=='$line'+3{printf("%2d  %7.3f     %5.3f %s %1d \n",$1,'$energy',$3,$4,$5)}' in.in1 > tmp
      mv tmp in.in1
    done
  fi

  # if -o option, get optimized energies from in.in1
  if [ $oop == True ]; then
    check $dir/$optdir/$a\Ang/in/in.in1
    awk 'NR<3{print}' in.in1 > _tmp
    awk 'NR>=3&&$1!="K-VECTORS"{print}' $dir/$optdir/$a\Ang/in/in.in1 >> _tmp
    awk 'END{print}' in.in1 >> _tmp
    mv _tmp in.in1
  fi
  
  # change smear, sigma, gmax in in.in2 to provided values
  if [ $smear == Fermi ]; then
    smearStr=`echo  TEMPS $sigma | awk '{printf("%-6s   %8.6f       ",$1,$2*'$toRy')}'`
  else
    smearStr=`echo $smear $sigma | awk '{printf("%-6s   %8.6f       ",$1,$2*'$toRy')}'`
  fi
  gmaxStr=`echo $gmax | awk '{printf("%6.2f          ",$1)}'`
  valStr=`awk 'NR==2{printf("%5.2f",$2+'$addCharge')}' in.in2`
  sed -e '2s/^\(.\{15\}\).\{5\}\(.\)/\1'"$valStr"'\2/' \
      -e 's/.*\((GAUSS,ROOT.*\)/'"$smearStr"'\1/' \
      -e 's/.*\(GMAX\)/'"$gmaxStr"'\1/' in.in2 > tmp; mv tmp in.in2

  # create k-points
  $wienPath/x kgen > /dev/null << EOF
$kpAll
0
$sh
EOF

  # generate starting density and check if corresponding FFT mesh needs to be changed to given values
  $wienPath/x dstart

  # check if FM calculation
  if [ $spin == FM ]; then
    $wienPath/x dstart -up
    $wienPath/x dstart -dn
  fi

  if [ $fft == default ]; then
    mv in.in0_std in.in0
  else
    fftStr=`echo $fft | awk '{printf(" %3d %3d %3d",$1,$1,$1)}'`
    sed '/IFFT-parameters/s/^[ ]*[0-9]*[ ]*[0-9]*[ ]*[0-9]*\(.*\)/'"$fftStr"'\1/' in.in0_std > in.in0
  fi

  cd $dir; echo `pwd`/$f >> jobList
done

echo; rm template.struct

# print note if additional charge added
echo $addCharge | awk '$1>0{print ""; for (i=1;i<70;i++) printf("-"); print "";
  print " NOTE: Additional charge has been added ",$1," electrons"; for (i=1;i<70;i++) printf("-"); print ""}'

writeLog $*

