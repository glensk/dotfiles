#!/bin/bash

# following 3 lines must always be present
path=`set | grep BASH_SOURCE | sed 's|.*\"\(.*\)/[^/]*\".*|\1|' | sed '/BASH_SOURCE=/s/.*/\./'`
script=`set | grep BASH_SOURCE | sed 's|.*\".*/\([^/]*\)\".*|\1|' | sed '/BASH_SOURCE=/s/.*/\./'`
options=$*; . $path/utilities/functions.include; checkOptions "-h -help -p -I -K -P -a -k -ff -f"

# small help for options
# space between option and value must be 1
# space between option and Comment or value and Comment must be >2
h=`getOption -h`
if [ $h == True ]; then
  usage $script
  printOptions "-p         create parameters.dat and exit" \
               "-I         create INCAR_{01_DFT,02_EXX,03_diag,04_ACFDT,05_HFc} templates and exit" \
               "-K         create KPOINTS template and exit" \
               "-P type    create POSCAR template of type (fcc/bcc/hcp/dhcp/sc) and exit" \
               "-a type    create create all of the above and exit" \
               "-k kpc     override automatic setting of kpoint along c by kpc (hcp only)" \
               "-ff        force running even if EDIFF>1e-5" \
               "-f         force template or folder creation even if existing"
  exit
fi

# detailed help
help=`getOption -help`
if [ $help == True ]; then
  details $script
  echo2 "   creates input folders for a cold curve (T=0K) run"
  echo2 "   parameters.dat file must be present" \
        "   (generate with -p)"
  echo2 "   INCAR KPOINTS POSCAR templates must be present" \
        "   (generate with -I -K -P bcc or fcc or hcp or dhcp or sc)"
  echo2 "   typical script sequence:"       \
        "     \033[1mcreateFolders_murn.sh\033[0m -a ..." \
        "     \033[1mcreateFolders_murn.sh\033[0m"        \
        "     collectOUTCARs.sh"               \
        "     extractEnergies.sh"              \
        "     getE0KFit.sh           (cmpc01)" \
        "     getThermodynamics.sh"
  exit
fi

# if -a option we create all templates and exit
all=`getOption -a`

# check if force overwriting templates or folders
overwrite=`getOption -f`

# if applies create parameters.dat
genPar=`getOption -p`
if [ $genPar == True -o $all == True ]; then
  if [ -e parameters.dat -a $overwrite != True ]; then error "parameters.dat existing; use -f to overwrite"; fi
  echo "
aLats=3.22 3.24 3.26    # Angstrom                example bcc Ti 4 valence electrons
cutoff=300              # eV
kp=8                    # gives NxNxN for fcc/bcc and NxNxInt(N/cBya)
cBya=1                  # 1 for bcc/fcc/sc and c/a ratio for hcp or hcp (ideal=1.632993161855452)

# IMPORTANT: cutoffGW should be 2/3 of cutoff for the extrapolation of the ACFDT energy to work
# 2/3 cutoff is obtained by setting DEFAULT
cutoffGW=DEFAULT        # cutoff for response function
nbandsACFDT=all         # all takes all possible determined by cutoff
nomega=16               # number of frequency integration points; 16 seems to be good
omegatl=800             # max. frequency for integration (eV); 800 seems good

addString=              # if addINCAR is used it is useful to add some string to name; leave blank otherwise
addINCAR=               # add additional flags to INCAR, e.g., addINCAR=ISPIN=2; leave blank if not needed
                        # more flags can be added by separating them with ;
                        # the flags are added to each INCAR of every step (01_DFT,02_EXX,03_diag,04_ACFDT,05_HFc)
" > parameters.dat
  echo; echo "parameters.dat written  <==  mind: parameters.dat needs to be adjusted"
fi

# if applies create INCAR
genINC=`getOption -I`
if [ $genINC == True -o $all == True ]; then
  if [ \( -e INCAR_00_noXC -o -e INCAR_01_DFT -o -e INCAR_02_EXX -o -e INCAT_03_diag -o -e INCAR_04_ACFDT -o -e INCAR_05_HFc \) -a $overwrite != True ]; then
    error "one or more of the INCARs existing; use -f to overwrite"
  fi

  echo "
 ENCUT  = xxxCUTOFFxxx
 ISMEAR =    1
 SIGMA  =   0.2

 ADDGRID=    TRUE
 PREC   =    Accurate
 PRECFOCK =  Accurate
 LMAXFOCKAE = 2
 NMAXFOCKAE = 2

 LWAVE  =      F   
           
 ALGO = EIGENVAL ; NELM = 1
 LHFCALC = .TRUE. ; AEXX = 0.0 ; ALDAC = 0.0 ; AGGAC = 0.0
 ALDAX = 0.0 ; AGGAX = 0.0

" > INCAR_00_noXC
  echo; echo "INCAR_00_noXC written"

  echo "
 NPAR = 1

 ENCUT  = xxxCUTOFFxxx
 ISMEAR =    1
 SIGMA  =   0.2
 EDIFF  =   1E-8

 ADDGRID=    TRUE
 PREC   =    Accurate
 ALGO   =    NORMAL

 LWAVE  =      T   
 LCHARG =      F   
" > INCAR_01_DFT
  echo; echo "INCAR_01_DFT written"

  echo "
 ENCUT  = xxxCUTOFFxxx
 ISMEAR =    1
 SIGMA  =   0.2

 ADDGRID=    TRUE
 PREC   =    Accurate
 PRECFOCK =  Accurate
 LMAXFOCKAE = 2
 NMAXFOCKAE = 2

 LWAVE  =      F   
           
 ALGO = EIGENVAL ; NELM = 1
 LHFCALC = .TRUE. ; AEXX = 1.0 ; ALDAC = 0.0 ; AGGAC = 0.0

" > INCAR_02_EXX
  echo; echo "INCAR_02_EXX written"

  echo "
 ENCUT  = xxxCUTOFFxxx
 ISMEAR =    1
 SIGMA  =   0.2
 EDIFF  =   1E-8
 PREC   =  Accurate

 NBANDS = xxxNBANDSACFDTxxx
 ALGO   = Exact
 NELM   = 1
 ! LOPTICS = .TRUE. ! not needed for metals
" > INCAR_03_diag
  echo; echo "INCAR_03_diag written"

  echo "
 ENCUT   = xxxCUTOFFxxx
 ENCUTGW = xxxCUTOFFGWxxx
 ISMEAR  =    1
 SIGMA   =   0.2
 ADDGRID = TRUE
 PREC    = Accurate

 ALGO    = ACFDT
 NBANDS  = xxxNBANDSACFDTxxx
 NOMEGA  = xxxNOMEGAxxx
 OMEGATL = xxxOMEGATLxxx
" > INCAR_04_ACFDT
  echo; echo "INCAR_04_ACFDT written"

  echo "
 ENCUT   = xxxCUTOFFxxx
 ISMEAR  =    1
 SIGMA   =   0.2
 ADDGRID =   TRUE
 PREC       =  Accurate
 PRECFOCK   =  Accurate
 LMAXFOCKAE =  2
 NMAXFOCKAE =  2

 ALGO    = HFc
 NOMEGA  = xxxNOMEGAxxx
 OMEGATL = xxxOMEGATLxxx
" > INCAR_05_HFc
  echo; echo "INCAR_05_HFc written"
fi

# if applies create KPOINTS
genKPO=`getOption -K`
if [ $genKPO == True -o $all == True ]; then
  if [ -e KPOINTS -a $overwrite != True ]; then error "KPOINTS existing; use -f to overwrite"; fi
  echo "K-Points
 0
Gamma
 xxxKPxxx xxxKPxxx xxxKPCxxx
 0 0 0
" > KPOINTS
  echo; echo "KPOINTS written"
fi

# copy POSCAR
genPOS=`getOption -P`
if [ $genPOS == True -o $all == True ]; then
  if [ -e POSCAR -a $overwrite != True ]; then error "POSCAR existing; use -f to overwrite"; fi
  if [ $genPOS == True ]; then type=`getValue -P`; fi
  if [ $all == True ]; then type=`getValue -a`; fi
  if [ -z "$type" ]; then error "type missing"; fi

  case "$type" in
    sc) echo "Murn sc
xxxALATxxx
1 0 0
0 1 0
0 0 1
1
Cartesian
0 0 0" > POSCAR;;

    bcc) echo "Murn bcc
xxxALATxxx
-.5 .5 .5
.5 -.5 .5
.5 .5 -.5
1
Cartesian
0 0 0" > POSCAR;;

   fcc) echo "Murn fcc
xxxALATxxx
0 .5 .5
.5 0 .5
.5 .5 0
1
Cartesian
0 0 0" > POSCAR;;

    hcp) sin=0.866025403784438597 # sin(pi/3) where pi/3 is angle between hcp lattice vectors
         echo "Murn hcp
xxxALATxxx
1.0  0.0  0.0
0.5  $sin  0.0
0.0  0.0  xxxCBYAxxx
2
Direct
0.000000000000 0.000000000000 0.0
0.333333333333 0.333333333333 0.5
    " > POSCAR;;

    dhcp) sin=0.866025403784438597 # sin(pi/3) where pi/3 is angle between hcp lattice vectors
         echo "Murn dhcp
xxxALATxxx
1.0  0.0  0.0
0.5  $sin  0.0
0.0  0.0  xxxCBYAxxx
4
Direct
0.000000000000 0.000000000000 0.0
0.333333333333 0.333333333333 0.25
0.000000000000 0.000000000000 0.5
0.666666666667 0.666666666667 0.75
    " > POSCAR;;

    *) "structure type $type not known";;
  esac

  echo; echo "POSCAR written"
fi

# exit if we created some template
if [ $genPar == True -o $genINC == True -o $genKPO == True -o $genPOS == True -o $all == True ]; then exit; fi


# from here on folder creation
# check if all input files available
input="parameters.dat INCAR_00_noXC INCAR_01_DFT INCAR_02_EXX INCAR_03_diag INCAR_04_ACFDT INCAR_05_HFc KPOINTS POSCAR"
for i in $input; do check $i; done
if [ ! -e POTCAR ]; then error "POTCAR missing; remember to take GW POTCAR"; fi

# read in all parameters from parameters.dat
aLats=`get aLats`; cutoff=`get cutoff`; kp=`get kp`; addString=`get addString`; cBya=`get cBya`
cutoffGW=`get cutoffGW`; nbandsACFDT=`get nbandsACFDT`; addINCAR=`get addINCAR`; nomega=`get nomega`; omegatl=`get omegatl`
checkInput "$aLats" "$cutoff" "$kp" "$cBya" "$cutoffGW" "$nbandsACFDT" "$nomega" "$omegatl"

# check if structure type known
type=`awk 'NR==1{print $2}' POSCAR`
echo type $type
if [ "$type" != "sc" -a "$type" != "bcc" -a "$type" != "fcc" -a "$type" != "hcp" -a "$type" != "dhcp" ]; then
  error "structure type in POSCAR not known";
fi

# if dhcp, scale c axis (which now equals cBya) times 2
if [ "$type" == "dhcp" ]; then
  cBya=`echo $cBya | awk '{printf "%.10f",$cBya*2}'`
fi

# scale KPC if hcp or dhcp
if [ $type == hcp -o $type == dhcp ]; then
  kpc=`echo $kp $cBya | awk '{printf("%d",$1/$2)}'`
  if [ "`getOption -k`" = True ]; then
    kpc=`getValue -k`
    if [ "$kpc" == "" ]; then error "no value to -k option provided"; fi
    c=`checkInteger $kpc`
    if [ $c != ok ]; then error "value give to -k is not an integer"; fi
  fi
else
  kpc=$kp
fi

# check if templates ok
cc=`grep -e xxxCUTOFFxxx -e "ALGO.*=.*EIGENVAL" INCAR_00_noXC | wc -l`
if [ $cc != 2 ]; then error "INCAR_00_noXC template wrong"; fi
cc=`grep -e xxxCUTOFFxxx INCAR_01_DFT | wc -l`
if [ $cc != 1 ]; then error "INCAR_01_DFT template wrong"; fi
cc=`grep -e xxxCUTOFFxxx -e "ALGO.*=.*EIGENVAL" INCAR_02_EXX | wc -l`
if [ $cc != 2 ]; then error "INCAR_02_EXX template wrong"; fi
cc=`grep -e xxxCUTOFFxxx -e xxxNBANDSACFDTxxx -e "ALGO.*=.*Exact" INCAR_03_diag | wc -l`
if [ $cc != 3 ]; then error "INCAR_03_diag template wrong"; fi
cc=`grep -e xxxCUTOFFxxx -e xxxNBANDSACFDTxxx -e xxxNOMEGAxxx -e xxxOMEGATLxxx -e xxxCUTOFFGWxxx -e "ALGO.*=.*ACFDT" INCAR_04_ACFDT | wc -l`
if [ $cc != 6 ]; then error "INCAR_04_ACFDT template wrong"; fi
cc=`grep -e xxxCUTOFFxxx -e xxxNOMEGAxxx -e xxxOMEGATLxxx -e "ALGO.*=.*HFc" INCAR_05_HFc | wc -l`
if [ $cc != 4 ]; then error "INCAR_05_HFc template wrong"; fi
cc=`grep -e xxxKPxxx -e xxxKPCxxx KPOINTS | xargs -n1 | grep -e xxxKPxxx -e xxxKPCxxx | wc -l`
if [ $cc != 3 ]; then error "KPOINTS template wrong"; fi
cc=`grep -e xxxALATxxx POSCAR | wc -w`
if [ $cc != 1 ]; then error "POSCAR template wrong"; fi
cc=`grep -e xxxCBYAxxx POSCAR | wc -l`
if [[ ( $type == "hcp" || $type == "dhcp" ) && $cc != 1 ]]; then error "type hcp or dhcp but no xxxCBYAxxx in POSCAR"; fi

# check EDIFF
ediffForce=`getOption -ff`
cc=`sed -nr 's/^[ ]*EDIFF[ ]*=([^;!]*).*/\1/p' INCAR_01_DFT | awk '{if ($1>1e-8) print "larger"}'`
if [ "$cc" == larger -a $ediffForce != True ]; then error "EDIFF>1e-8; run with -ff option to override"; fi

# if nbandsACFDT==all we have to determine the actual number from a pseudo vasp run
if [ "$nbandsACFDT" == all ]; then
  if [ "`hostname`" != cmmc001 -a "`hostname`" != cmmc002 ]; then
    error "nbandsACFDT=all is supported only on cmmc"
  fi
  echo "running quick vasp calculation to get nbandsACFDT"
  mkdir _tmp_vasp_run
  sed -e 's/xxxCUTOFFxxx/'$cutoff'/' -e 's/.*NELM.*/NELM = 1/' \
      -e 's/.*LWAVE.*/LWAVE = FALSE/' INCAR_01_DFT               > _tmp_vasp_run/INCAR
  echo 'NELM = 1' >> _tmp_vasp_run/INCAR
  sed -e 's/xxxKPxxx/1/g' -e 's/xxxKPCxxx/1/' KPOINTS            > _tmp_vasp_run/KPOINTS
  a=`echo $aLats | awk '{print $1}'`
  sed -e 's/xxxALATxxx/'$a'/' -e 's/xxxCBYAxxx/'$cBya'/' POSCAR  > _tmp_vasp_run/POSCAR
  cp POTCAR _tmp_vasp_run/POTCAR
  
  cd _tmp_vasp_run
  /bin/tcsh $path/run_scripts/run.vasp5.ser.local.cmmc
  nbandsACFDT=`grep "maximum number of plane-waves" OUTCAR | awk '{print $NF}'`
  echo; echo "nbandsACFDT = $nbandsACFDT"

  # vasp will change NBANDS in parallel mode, we do the change ourselves to avoid a later change
  # we assume that we will be running on 20 cores
  cores=20
  nbandsACFDT=`echo $nbandsACFDT $cores | awk '{print int(($1+$2-1)/$2)*$2}'`
  echo "adjustment to parallel mode assuming $cores cores"
  echo "nbandsACFDT = $nbandsACFDT"

  cd ..
  rm -fr _tmp_vasp_run
fi

# check if folder exist and should be overwritten
for a in $aLats; do
  f=$cutoff\eV_$kp\x$kp\x$kpc\kp__${nbandsACFDT}_${cutoffGW}eV_${nomega}_${omegatl}eV$addString/$a\Ang
  if [ -e $f -a $overwrite != True ]; then
    error "folder $f exists; run with -f option to force overwriting";
  fi
done

# create the input folders
echo; rm -f jobList
for a in $aLats; do
  f=$cutoff\eV_$kp\x$kp\x$kpc\kp__${nbandsACFDT}_${cutoffGW}eV_${nomega}_${omegatl}eV$addString/$a\Ang
  echo $f; rm -fr $f; mkdir -p $f
  mkdir -p $f/00_noXC/ $f/01_DFT/ $f/02_EXX/ $f/03_diag/ $f/04_ACFDT/ $f/05_HFc/
  sed -e 's/xxxCUTOFFxxx/'$cutoff'/'           INCAR_00_noXC > $f/00_noXC/INCAR
  sed -e 's/xxxCUTOFFxxx/'$cutoff'/'           INCAR_01_DFT > $f/01_DFT/INCAR
  sed -e 's/xxxCUTOFFxxx/'$cutoff'/'           INCAR_02_EXX > $f/02_EXX/INCAR
  sed -e 's/xxxCUTOFFxxx/'$cutoff'/'     \
      -e 's/xxxNBANDSACFDTxxx/'$nbandsACFDT'/' INCAR_03_diag > $f/03_diag/INCAR
  sed -e 's/xxxCUTOFFxxx/'$cutoff'/'     \
      -e 's/xxxCUTOFFGWxxx/'$cutoffGW'/' \
      -e 's/xxxNOMEGAxxx/'$nomega'/'     \
      -e 's/xxxOMEGATLxxx/'$omegatl'/'   \
      -e 's/xxxNBANDSACFDTxxx/'$nbandsACFDT'/' INCAR_04_ACFDT | \
      sed '/ENCUTGW *= *DEFAULT/d' > $f/04_ACFDT/INCAR   # if cutoffGW==DEFAULT remove flag in INCAR to let vasp determine default
  sed -e 's/xxxCUTOFFxxx/'$cutoff'/'     \
      -e 's/xxxNOMEGAxxx/'$nomega'/'     \
      -e 's/xxxOMEGATLxxx/'$omegatl'/'         INCAR_05_HFc > $f/05_HFc/INCAR

  if [ -n "$addINCAR" ]; then
    echo "$addINCAR" | xargs -n1 -d \; >> $f/00_noXC/INCAR;
    echo "$addINCAR" | xargs -n1 -d \; >> $f/01_DFT/INCAR;
    echo "$addINCAR" | xargs -n1 -d \; >> $f/02_EXX/INCAR;
    echo "$addINCAR" | xargs -n1 -d \; >> $f/03_diag/INCAR;
    echo "$addINCAR" | xargs -n1 -d \; >> $f/04_ACFDT/INCAR;
    echo "$addINCAR" | xargs -n1 -d \; >> $f/05_HFc/INCAR;
  fi
  sed -e 's/xxxKPxxx/'$kp'/g' \
      -e 's/xxxKPCxxx/'$kpc'/' KPOINTS > $f/KPOINTS
  sed -e 's/xxxALATxxx/'$a'/' \
      -e 's/xxxCBYAxxx/'$cBya'/' POSCAR > $f/POSCAR
  cp POTCAR $f/POTCAR

  cp $f/{KPOINTS,POSCAR,POTCAR} $f/00_noXC/
  cp $f/{KPOINTS,POSCAR,POTCAR} $f/01_DFT/
  cp $f/{KPOINTS,POSCAR,POTCAR} $f/02_EXX/
  cp $f/{KPOINTS,POSCAR,POTCAR} $f/03_diag/
  cp $f/{KPOINTS,POSCAR,POTCAR} $f/04_ACFDT/
  cp $f/{KPOINTS,POSCAR,POTCAR} $f/05_HFc/

  echo `pwd`/$f >> jobList
done

echo; echo " jobList file written"

