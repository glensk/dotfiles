#!/bin/bash

#-----set parameters and paths------------------------------------
f="OUTCAR OUTCAR.gz"
#-----------------------------------------------------------------


path=`set | grep BASH_SOURCE | sed 's|.*\"\(.*\)/[^/]*\".*|\1|' | sed '/BASH_SOURCE=/s/.*/\./'`
script=`set | grep BASH_SOURCE | sed 's|.*\".*/\([^/]*\)\".*|\1|' | sed '/BASH_SOURCE=/s/.*/\./'`
options=$*; . $path/utilities/functions.include; checkOptions "-h -help -i -o"

# small help for options
h=`getOption -h`
if [ $h == True ]; then
  usage $script
  printOptions "-i inpFile  input file to use (default: check for $f)" \
               "-o          extract occupancies additionally"
  exit
fi

# detailed help
help=`getOption -help`
if [ $help == True ]; then
  details $script
  echo2 "   $script extracts the Kohn-Sham eigenvalues from an OUTCAR-like file"
  echo2 "   for ISPIN=1 the eigenvalues are written to an EIGENVALUES file"             \
        "   for ISPIN=2 two files are produced corresponding to spin up and down"
  echo2 "   in addition a weights file is written containing the weights of the"        \
        "   irreducible k-points; the weights file is needed for constructing the"      \
        "   DOS and for calculating the electronic free energy"
  echo2 "   the electronic DOS can be obtained subsequently using (depending on spin):" \
        "     getDOS.sh -w EIGENVALUES"                                                 \
        "   or"                                                                         \
        "     getDOS.sh -w EIGENVALUES_up; mv DOS DOS_up"                               \
        "     getDOS.sh -w EIGENVALUES_down; mv DOS DOS_down"
  echo2 "   to obtain a reasonable DOS the k-point mesh needs to be dense (e.g. 20^3)"  \
        "   it might be also necessary to change the smearing to resolve the DOS nicely"\
        "   a good value for smearing should be 0.1 eV, i.e. use:"                      \
        "     getDOS.sh -w  EIGENVALUES 0.1"
  echo2 "   the electronic free energy can be obtained using:"                          \
        "     getFelFromEIGENVALUES.sh"                                                 \
        "   or"                                                                         \
        "     getFelFromEIGENVALUES.sh -e EIGENVALUES_up"                               \
        "     getFelFromEIGENVALUES.sh -e EIGENVALUES_down"                             \
        "   add Fel from up and down and divide by two"
  echo2 "   $script can optionally (-o option) extract the occupancies"                 \
        "   from an OUTCAR into a file OCCUPANCIES (OCCUPANCIES_{up,down} for ISPIN=2)" \
        "   which then contains the eigenvalues and occupancies (first/second column)"
  echo2 "   the eigenvalues and occupancies in the OCCUPANCIES files are in the same"   \
        "   order as in the EIGENVALUES/weights files; additionally (for -o option)"    \
        "   OCCUPANCIES_sorted_for_plotting files are created which contain sorted"     \
        "   eigenvalues which are useful for plotting"
  echo2 "   the OCCUPANCIES file (not the sorted one) can be used for e.g. recovering"  \
        "   the number of electrons NELECT as:"                                         \
        "     paste OCCUPANCIES weights | awk 'BEGIN{s=0} {s=s+\$2*\$3} END{print s}'"
  echo2 "   if the OUTCAR contains several ionic steps (e.g. relaxation) each with"     \
        "   eigenvalues then the last step is used for extraction"

  exit
fi

# get and check input file
input=`getOption -i`
if [ $input == True ]; then
  file=`getValue -i`; check $file;
else
  file=0; for i in $f; do if [ -e $i ]; then file=$i; break; fi; done
  if [ $file == 0 ]; then error "no input file existing; checked for $f"; fi
fi
echo; echo " extracting EIGENVALUES from $file ..."

# check whether several ionic steps (relaxation) with eigenvalues, then take only last one
nsteps=`zgrep "potential at core" $file | awk 'END{print NR}'`

# pre-extract the needed information to save time (typically long OUTCAR)
zgrep ".*" $file | awk 'BEGIN{write=0;wsum=0; check=0; nstepsCount=1}
                        /ISPIN  =/{spin=$3}
                        /NELECT =/{nelec=$3}
                        /irreducible k-points/{irr=$2;l=NR}
                        NR>l+3&&NR<=l+3+irr{print $NF;wsum=wsum+$NF}
                        /number of bands    NBANDS/{nbands=$NF}
                        /potential at core/{if (nstepsCount=='$nsteps') write=1; else nstepsCount=nstepsCount+1}
                        write==2{check=check+1; printf("%.4f %.5f\n",$2-efermi,$3); if ($1==nbands) write=1}
                        /E-fermi :/{if (write==1) efermi=$3}
                        /band No.  band energies/{if (write==1) write=2}
                        END{if (irr*nbands*spin!=check) error="True"; else error="False";
                            print irr,nbands,spin,wsum,error,check,nelec}' > _tmp_EIGENVALUES

irr=`   tail -n1 _tmp_EIGENVALUES | awk '{print $1}'`
nbands=`tail -n1 _tmp_EIGENVALUES | awk '{print $2}'`
spin=`  tail -n1 _tmp_EIGENVALUES | awk '{print $3}'`
wsum=`  tail -n1 _tmp_EIGENVALUES | awk '{print $4}'`
error=` tail -n1 _tmp_EIGENVALUES | awk '{print $5}'`
nelec=` tail -n1 _tmp_EIGENVALUES | awk '{print $7}'`

if [ $error == True ]; then error "input file $file is corrupt (nbands*irr*spin != eigenvalues)"; fi

# only spin 1 or 2 supported
if [ $spin != 1 -a $spin != 2 ]; then rm _tmp_EIGENVALUES; error "ISPIN=$spin not supported"; fi

# write EIGENVALUES and weights and number_of_electrons_NELECT
awk 'NR<='$irr' {for (i=1;i<='$nbands';i++) printf("%.15f\n", $1/'$wsum')}' _tmp_EIGENVALUES > weights
awk 'NR>'$irr'&&NR<='$irr'+'$irr'*'$nbands' {print $1}' _tmp_EIGENVALUES > EIGENVALUES
echo $nelec > number_of_electrons_NELECT


# check for -o option controlling OCCUPANCIES output
oOp=`getOption -o`
if [ $oOp == True ]; then 
  awk 'NR>'$irr'&&NR<='$irr'+'$irr'*'$nbands' {print $1,$2}' _tmp_EIGENVALUES > OCCUPANCIES
  cat OCCUPANCIES | sort -n > OCCUPANCIES_sorted_for_plotting
  occ=" OCCUPANCIES  OCCUPANCIES_sorted_for_plotting "
fi

if [ $spin == 1 ]; then
  echo " files written: EIGENVALUES $occ weights  number_of_electrons_NELECT"
else
  # write files for both spin channels if spin==2
  mv EIGENVALUES EIGENVALUES_up
  awk 'NR>'$irr'+'$irr'*'$nbands'&&NR<='$irr'+'$irr'*'$nbands'*'$spin' {print $1}' _tmp_EIGENVALUES > EIGENVALUES_down
  if [ $oOp == True ]; then
    mv OCCUPANCIES OCCUPANCIES_up
    mv OCCUPANCIES_sorted_for_plotting OCCUPANCIES_up_sorted_for_plotting
    awk 'NR>'$irr'+'$irr'*'$nbands'&&NR<='$irr'+'$irr'*'$nbands'*'$spin' {print $1,$2}' _tmp_EIGENVALUES > OCCUPANCIES_down
    cat OCCUPANCIES_down | sort -n > OCCUPANCIES_down_sorted_for_plotting
    occ=" OCCUPANCIES_{up,down}  OCCUPANCIES_{up,down}_sorted_for_plotting "
  fi
  echo " files written: EIGENVALUES_{up,down} $occ weights  number_of_electrons_NELECT"
fi
rm _tmp_EIGENVALUES

