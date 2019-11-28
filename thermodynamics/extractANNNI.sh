#!/bin/bash

#-----set parameters and paths------------------------
scale=16.021765        # meV/angstrom^2 --> mJ/m^2
outFolder=SFE_results
meVByAngstrom3ToGPa=0.1602176462
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
  printNeeded  "-f fccFolder hcpFolder [dhcpFolder]   for ANNNI [or DANNNI] (order of folders must be preserved)"
  printOptions "-n      do not rerun hcp/dhcp" \
               "-t      T=0K ANNNI/DANNNI"     \
               "-H      harmonic ANNNI/DANNNI" \
               "-P      const P ANNNI/DANNNI"  \
               "-d      delta hcp-fcc"         \
               "-a      = -t -H -P -d"
  echo "Example: $script PBE_fcc_q108  PBE_hcp_q288  PBE_dhcp_q128" 1>&2   
  exit
fi

# detailed help
help=`getOption -help`
if [ $help == True ]; then
  details "$script -f fccFolder hcpFolder [dhcpFolder]"
  echo2 "   calculate stacking fault energy in the ANNNI or DANNNI approximation" \
        "   from a set of getThermodynamics folders given to the -f option"
  echo2 "   by default the hcp/dhcp folders are recalculated on the volume" \
        "   expansion of fcc to produce the correct const V SFE; the explicit rerunning" \
        "   can be disabled with -n"
  echo2 "   sfe is in the units mJ/m^2"
  exit
fi


# rerun option
if [ "`getOption -n`" == True ]; then rerunvol=False; else rerunvol=True; fi


# -a option ?
if [ `getOption -a` == True ]; then options="$options -t -H -P -d"; fi


# get folders from -f option
if [ "`getOption -f`" != True ]; then error "-f option with folders mandatory (see -h)"; fi
folders=`getValue -f`
if [ "`echo $folders | awk 'NF!=2&&NF!=3{print "error"}'`" == error ]; then error "number of folders given to -f option incorrect"; fi
fccFolder=`echo $folders | awk '{print $1}'`
hcpFolder=`echo $folders | awk '{print $2}'`
dhcpFolder=`echo $folders | awk '{print $3}'`
echo;
echo "  fccFolder:  $fccFolder"
echo "  hcpFolder:  $hcpFolder"
if [ -n "$dhcpFolder" ]; then echo "  dhcpFolder: $dhcpFolder"; fi


# ouput string
outannni=`echo $fccFolder $hcpFolder | sed 's|/||g' | awk '{print $1"___"$2}'`
outdannni=`echo $fccFolder $hcpFolder $dhcpFolder | sed 's|/||g' | awk '{print $1"___"$2"___"$3}'`


# check if output folder exists, create otherwise
if [ ! -e $outFolder ]; then mkdir $outFolder; fi

# check whether only EVinet files available, then produce only T=0K SFE
c1=`ls -1 $fccFolder/F* 2> /dev/null | awk 'END{print NR}'`
c2=`ls -1 $hcpFolder/F* 2> /dev/null | awk 'END{print NR}'`

opT=`getOption -t`
if [ $c1 == 0 ]; then echo; echo ' No free energy files in fcc folder, trying T=0K'; opT=True; fi
if [ $c2 == 0 ]; then echo; echo ' No free energy files in hcp folder, trying T=0K'; opT=True; fi


# fcc volume expansion
  fccVol="$fccFolder/output_0.0001GPa/volume_expansion"
if [ $c1 == 0 -o $c2 == 0 ]; then
  mkdir -p $fccFolder/output_0.0001GPa
  check $fccFolder/EVinet
  awk '{print 1,$2}' $fccFolder/EVinet > $fccVol
else
  check $fccVol
fi


# ------ T0K values   START -----------------------------------------------------------------------------------
if [ $opT == True ]; then
  echo; echo; echo "  === T0K ==="
  check $fccFolder/EVinet $hcpFolder/EVinet
  T0Kvol=`awk '{print $2}' $fccFolder/EVinet`
  T0Kfcc=`awk '{print $1}' $fccFolder/EVinet`
  T0Khcp=`awk '{V='$T0Kvol'; E0=$1; V0=$2; BM=$3/'$meVByAngstrom3ToGPa'; Bder=$4;
                Vhcp = E0 + (4*BM*V0)/(BMder - 1)^2 - 2*V0*BM*(BMder - 1)^(-2)*(5 + 3*BMder*((V/V0)^(1/3) - 1) - 3*(V/V0)^(1/3))*exp(-(3/2)*(BMder - 1)*((V/V0)^(1/3) - 1));
                printf "%8.15f",Vhcp}' $hcpFolder/EVinet`

  cat $fccVol | \
    awk '{T    = $1
          vol  = '$T0Kvol'
          aLat = (4*vol)^(1/3)
          A    = sqrt(3)/4 * aLat^2
          fcc  = '$T0Kfcc'
          hcp  = '$T0Khcp'
          SFE  = '$scale' * 2 * (hcp - 1*(fcc))/A;

          printf "%.1f  %8.3f\n", T, SFE}' > $outFolder/ANNNI___${outannni}__0K

  if [ -n "$dhcpFolder" ]; then
    check $dhcpFolder/EVinet
    T0Kdhcp=`awk '{print $1}' $dhcpFolder/EVinet`
    cat $fccVol | \
     awk '{T    = $1
           vol  = '$T0Kvol'
           aLat = (4*vol)^(1/3)
           A    = sqrt(3)/4 * aLat^2
           fcc  = '$T0Kfcc'
           hcp  = '$T0Khcp'
           dhcp = '$T0Kdhcp'
           SFE  = '$scale' * (hcp + 2*(dhcp) - 3*(fcc))/A;
   
           printf "%.1f  %8.3f\n", T, SFE}' > $outFolder/DANNNI__${outdannni}__0K
  fi
fi
# ------ T0K values   END -------------------------------------------------------------------------------------


# if no free energy files either for fcc or hcp folder then stop here
if [ $c1 == 0 -o $c2 == 0 ]; then exit; fi


# check whether fcc and hcp have same mesh
check $hcpFolder/Fqh
c=`paste $fccVol $hcpFolder/Fqh | awk 'BEGIN{f="ok"} {if ($1!=$3) f="error"} END{print f}'`
if [ "$c" == "error" ]; then error "ERROR: mismatch in the T mesh between fcc and hcp; rerun either one on the mesh of the other"; fi


if [ -n "$dhcpFolder" ]; then
  # check whether fcc and dhcp have same mesh
  check $dhcpFolder/Fqh
  c=`paste $fccVol $dhcpFolder/Fqh | awk 'BEGIN{f="ok"} {if ($1!=$3) f="error"} END{print f}'`
  if [ "$c" == "error" ]; then error "ERROR: mismatch in the T mesh between fcc and dhcp; rerun either one on the mesh of the other"; fi
fi



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
    
    # now for hcp
    cp $fccFolder/volume_expansion $hcpFolder/
    cd $hcpFolder
    $path/getThermodynamics.sh
    cp output_0.0001GPa/free_energy_at_volume_expansion free_energy_at_volume_expansion_harmonic
    cd $dir
    
    # and for dhcp if needed
    if [ -n "$dhcpFolder" ]; then
      cp $fccFolder/volume_expansion $dhcpFolder/
      cd $dhcpFolder
      $path/getThermodynamics.sh
      cp output_0.0001GPa/free_energy_at_volume_expansion free_energy_at_volume_expansion_harmonic
      cd $dir
    fi
  fi

  fcc=$fccFolder/free_energy_at_volume_expansion_harmonic
  hcp=$hcpFolder/free_energy_at_volume_expansion_harmonic
  check $fcc $hcp

  paste $fccVol $fcc $hcp | \
    awk '{T    = $1
          vol  = $2
          aLat = (4*vol)^(1/3)
          A    = sqrt(3)/4 * aLat^2
          fcc  = $4
          hcp  = $6
          SFE  = '$scale' * 2 * (hcp - 1*(fcc))/A;
          printf "%.1f  %8.3f\n", T, SFE}' > $outFolder/ANNNI___${outannni}__harm


  if [ -n "$dhcpFolder" ]; then
    dhcp=$dhcpFolder/free_energy_at_volume_expansion_harmonic
    check $dhcp

    paste $fccVol $fcc $hcp $dhcp | \
      awk '{T    = $1;
            vol  = $2;
            aLat = (4*vol)^(1/3);
            A    = sqrt(3)/4 * aLat^2;
            fcc  = $4;
            hcp  = $6;
            dhcp = $8;
            SFE  = '$scale' * (hcp + 2*(dhcp) - 3*(fcc))/A;
            printf "%.1f  %8.3f\n", T, SFE}' > $outFolder/DANNNI__${outdannni}__harm
  fi
fi
# ====== harmonic values?   END ===============================================================================



# ------ full values   START ----------------------------------------------------------------------------------
echo; echo; echo "  === full ==="
if [ "$rerunvol" == "True" ]; then
  dir=`pwd`
  cp $fccVol $hcpFolder/
  cd $hcpFolder
  $path/getThermodynamics.sh
  cd $dir

  if [ -n "$dhcpFolder" ]; then
    cp $fccVol $dhcpFolder/
    cd $dhcpFolder
    $path/getThermodynamics.sh
    cd $dir
  fi
fi

fcc="$fccFolder/output_0.0001GPa/Gibbs_energy"
hcp="$hcpFolder/output_0.0001GPa/free_energy_at_volume_expansion"
check $fcc $hcp

paste $fccVol $fcc $hcp | \
  awk '{T    = $1
        vol  = $2
        aLat = (4*vol)^(1/3)
        A    = sqrt(3)/4 * aLat^2
        fcc  = $4
        hcp  = $6
        SFE  = '$scale' * 2 * (hcp - 1*(fcc))/A;
        printf "%.1f  %8.3f\n", T, SFE}' > $outFolder/ANNNI___${outannni}

if [ -n "$dhcpFolder" ]; then
  dhcp="$dhcpFolder/output_0.0001GPa/free_energy_at_volume_expansion"
  check $dhcp
  paste $fccVol $fcc $hcp $dhcp | \
    awk '{T    = $1;
          vol  = $2;
          aLat = (4*vol)^(1/3);
          A    = sqrt(3)/4 * aLat^2;
          fcc  = $4;
          hcp  = $6;
          dhcp = $8;
          SFE  = '$scale' * (hcp + 2*(dhcp) - 3*(fcc))/A;
          printf "%.1f  %8.3f\n", T, SFE}' > $outFolder/DANNNI__${outdannni}
fi
# ------ full values   END ------------------------------------------------------------------------------------



# ------ delta hcp - fcc START --------------------------------------------------------------------------------
if [ `getOption -d` == True ]; then
  echo; echo; echo "  === delta hcp - fcc ==="
  hcpG="$hcpFolder/output_0.0001GPa/Gibbs_energy"
  check $hcpG
  paste $hcpG $fcc | awk '{printf "%.1f  %8.3f\n", $1, $2-$4}' > $outFolder/deltaG__${outannni}
fi
# ------ delta hcp - fcc END ----------------------------------------------------------------------------------



# ====== const P?   START =====================================================================================
# needs to be here at end to have rerun in full values before
if [ `getOption -P` == True ]; then
  echo; echo; echo "  === const P ==="

  hcp="$hcpFolder/output_0.0001GPa/Gibbs_energy"
  check $hcp

  paste $fccVol $fcc $hcp | \
    awk '{T    = $1
          vol  = $2
          aLat = (4*vol)^(1/3)
          A    = sqrt(3)/4 * aLat^2
          fcc  = $4
          hcp  = $6
          SFE  = '$scale' * 2 * (hcp - 1*(fcc))/A;
          printf "%.1f  %8.3f\n", T, SFE}' > $outFolder/ANNNI___${outannni}__constP

  if [ -n "$dhcpFolder" ]; then
    dhcp="$dhcpFolder/output_0.0001GPa/Gibbs_energy"
    check $dhcp
    paste $fccVol $fcc $hcp $dhcp | \
      awk '{T    = $1;
            vol  = $2;
            aLat = (4*vol)^(1/3);
            A    = sqrt(3)/4 * aLat^2;
            fcc  = $4;
            hcp  = $6;
            dhcp = $8;
            SFE  = '$scale' * (hcp + 2*(dhcp) - 3*(fcc))/A;
            printf "%.1f  %8.3f\n", T, SFE}' > $outFolder/DANNNI__${outdannni}__constP
  fi
fi
# ====== const P?   END =======================================================================================

echo

