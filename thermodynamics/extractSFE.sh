#!/bin/bash

#-----set parameters and paths------------------------
# scale1=0.1666666666666667 # 1/6 not used anymore
scale2=16.021765          # meV/angstrom^2 --> mJ/m^2
outFolder=SFE_results
#-----------------------------------------------------

# following 3 lines must always be present
path=`set | grep BASH_SOURCE | sed 's|.*\"\(.*\)/[^/]*\".*|\1|' | sed '/BASH_SOURCE=/s/.*/\./'`
script=`set | grep BASH_SOURCE | sed 's|.*\".*/\([^/]*\)\".*|\1|' | sed '/BASH_SOURCE=/s/.*/\./'`
options=$*; . $path/utilities/functions.include; checkOptions "-h -help -f -n -t -H -P -d -a"

# small help for options
# space between option and value must be 1
# space between option and Comment or value and Comment must be >2
h=`getOption -h`
if [ $h == True ]; then
  usage $script
  printNeeded  "-f scZ SFEFolder refFolder    scZ supercell along Z, refFolder = fcc bulk (order of arguments must be preserved)"
  printOptions "-n      do not rerun hcp/dhcp" \
               "-t      T=0K SFE"              \
               "-H      harmonic SFE"          \
               "-P      const P SFE"           \
               "-a      = -t -H -P"
  echo "Example: $script -f 2 PBE_sfe_q PBE_ref_q  # when SF calculated in a 1x1x2 supercell, see -help" 1>&2   
  exit
fi

# detailed help
help=`getOption -help`
if [ $help == True ]; then
  details "$script -f fccFolder hcpFolder [dhcpFolder]"
  echo2 "   calculate explicit stacking fault energy from a set of" \
        "   getThermodynamics folders given to the -f option"
  echo2 "   by default the sfe folder is recalculated on the volume expansion" \
        "   of fcc to produce the correct const V SFE; the explicit rerunning" \
        "   can be disabled with -n"
  echo2 "   sfe is in the units mJ/m^2"
  echo2 "   the scZ argument corresponds to the number of repeats perpendicular" \
        "   to the stacking fault used in the calculation, e.g., 2 for 1x1x2 or" \
        "   for 2x3x2, 4 for 1x1x4, etc.; this is needed to get the scaling factor"
  echo2 "   scaleSFE.sh can be used to bring different free energy contributions" \
        "   to the same scZ value"
  exit
fi


# rerun option
if [ "`getOption -n`" == True ]; then rerunvol=False; else rerunvol=True; fi


# -a option ?
if [ `getOption -a` == True ]; then options="$options -t -H -P"; fi


# get scZ and folders from -f option
if [ "`getOption -f`" != True ]; then error "-f option with scZ and SFE/ref folders mandatory (see -h)"; fi
folders=`getValue -f`
if [ "`echo $folders | awk 'NF!=3{print "error"}'`" == error ]; then
  error "number of arguments given to -f option incorrect"
fi
scZ=`   echo $folders | awk '{print $1}'`
sfeFolder=`echo $folders | awk '{print $2}'`
fccFolder=`echo $folders | awk '{print $3}'`
if [ `checkInteger $scZ` != ok ]; then error "first argument of -a option not an integer (scZ)"; fi
echo
echo "  scZ: $scZ"
echo
echo "  sfeFolder: $sfeFolder"
echo "  fccFolder: $fccFolder"


# ouput string
out=`echo $sfeFolder $fccFolder | sed 's|/||g' | sed 's|thermo||g' | sed 's|fit_order|ord|g' | awk '{printf "%s_%s",$1,$2}'`


# fcc volume expansion
fccVol="$fccFolder/output_0.0001GPa/volume_expansion"
check $fccVol


# check whether sfe and fcc have same mesh
check $sfeFolder/Fqh
c=`paste $fccVol $sfeFolder/Fqh | awk 'BEGIN{f="ok"} {if ($1!=$3) f="error"} END{print f}'`
if [ "$c" == "error" ]; then error "ERROR: mismatch in the T mesh between sfe and fcc ref; rerun either one on the mesh of the other"; fi


# check if output folder exists, create otherwise
if [ ! -e $outFolder ]; then mkdir $outFolder; fi


# ------ T0K values   START -----------------------------------------------------------------------------------
if [ `getOption -t` == True ]; then
  echo; echo; echo "  === T0K ==="
  check $fccFolder/EVinet $sfeFolder/EVinet
  T0Kvol=`awk '{print $2}' $fccFolder/EVinet`
  T0Kfcc=`awk '{print $1}' $fccFolder/EVinet`
  T0Ksfe=`awk '{print $1}' $sfeFolder/EVinet`

  cat $fccVol | \
    awk '{T    = $1
          vol  = '$T0Kvol'
          aLat = (4*vol)^(1/3)
          A    = sqrt(3)/2 * aLat^2
          fcc  = '$T0Kfcc'
          sfe  = '$T0Ksfe'
          SFE  = '$scale2' * '$scZ'*6  * (sfe - 1*(fcc))/A;  # 6 is to the number of atoms per primitive cell

          printf "%.1f  %8.3f\n", T, SFE}' > $outFolder/SFE_____${out}__0K
fi
# ------ T0K values   END -------------------------------------------------------------------------------------


# ====== harmonic values?   START =============================================================================
if [ `getOption -H` == True ]; then
  echo; echo; echo "  === harm ==="

  if [ "$rerunvol" == True ]; then
    dir=`pwd`
    # prepare a file that contains T=0K fcc volume for all T
    awk 'NR==1{v=$2} {print $1,v}' $fccVol > $fccFolder/volume_expansion
    cd $fccFolder
    $path/getThermodynamics.sh
    cp output_0.0001GPa/free_energy_at_volume_expansion free_energy_at_volume_expansion_harmonic
    cd $dir
    
    # now for sfe
    cp $fccFolder/volume_expansion $sfeFolder/
    cd $sfeFolder
    $path/getThermodynamics.sh
    cp output_0.0001GPa/free_energy_at_volume_expansion free_energy_at_volume_expansion_harmonic
    cd $dir
  fi

  fcc=$fccFolder/free_energy_at_volume_expansion_harmonic
  sfe=$sfeFolder/free_energy_at_volume_expansion_harmonic
  check $fcc $sfe

  paste $fccVol $fcc $sfe | \
    awk '{T    = $1
          vol  = $2
          aLat = (4*vol)^(1/3)
          A    = sqrt(3)/2 * aLat^2
          fcc  = $4
          sfe  = $6
          SFE  = '$scale2' * '$scZ'*6  * (sfe - 1*(fcc))/A;  # 6 is to the number of atoms per primitive cell

          printf "%.1f  %8.3f\n", T, SFE}' > $outFolder/SFE_____${out}__harm
fi
# ====== harmonic values?   END ===============================================================================



# ------ full values   START ----------------------------------------------------------------------------------
echo; echo; echo "  === full ==="
if [ "$rerunvol" == "True" ]; then
  dir=`pwd`
  cp $fccVol $sfeFolder/
  cd $sfeFolder
  $path/getThermodynamics.sh
  cd $dir
fi

fcc="$fccFolder/output_0.0001GPa/Gibbs_energy"
sfe="$sfeFolder/output_0.0001GPa/free_energy_at_volume_expansion"
check $fcc $sfe

paste $fccVol $fcc $sfe | \
  awk '{T    = $1
        vol  = $2
        aLat = (4*vol)^(1/3)
        A    = sqrt(3)/2 * aLat^2
        fcc  = $4
        sfe  = $6
        SFE  = '$scale2' * '$scZ'*6  * (sfe - 1*(fcc))/A;  # 6 is to the number of atoms per primitive cell

        printf "%.1f  %8.3f\n", T, SFE}' > $outFolder/SFE_____${out}

# ------ full values   END ------------------------------------------------------------------------------------



# ====== const P?   START =====================================================================================
# needs to be here at end to have rerun in full values before
if [ `getOption -P` == True ]; then
  echo; echo; echo "  === const P ==="

  sfe="$sfeFolder/output_0.0001GPa/Gibbs_energy"
  check $sfe

  paste $fccVol $fcc $sfe | \
    awk '{T    = $1
          vol  = $2
          aLat = (4*vol)^(1/3)
          A    = sqrt(3)/2 * aLat^2
          fcc  = $4
          sfe  = $6
          SFE  = '$scale2' * '$scZ'*6  * (sfe - 1*(fcc))/A;  # 6 is to the number of atoms per primitive cell

          printf "%.1f  %8.3f\n", T, SFE}' > $outFolder/SFE_____${out}__constP
fi
# ====== const P?   END =======================================================================================

echo

