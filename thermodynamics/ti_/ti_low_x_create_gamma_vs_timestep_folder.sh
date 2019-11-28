#!/bin/sh
#sc=3                                                #for 3x3x3sc
#cutoff=300                                          #cutoff

#echo create a folder with the named like low_3x3x3sc_300eV_gamma_test
ti_low_0_create_Folders_vasp.sh -o -i
ti_low_0_create_Folders_vasp.sh -o -k
#[ -e "parameters.dat" ] && echo parameters.dat already exists, please remove it by hand && exit

#ti_low_0_create_Folders_vasp.sh -o -p


#echo we want to probe tmelt
sed -i 's|lambdas.*|lambdas=0.0|' parameters.dat
sed -i 's|seeds.*|seeds=0.0|' parameters.dat
#echo in der incar kan EDIFF=1E-1


EDIFF="1E-2"
NBANDS="300"
temps="1358"
time="10 8 4 3 5 2"                                         # POTIM in INCAR
#time="0.5 1 2 3 4 5 10"                                         # POTIM in INCAR
#gamma="0.0001 0.0005 0.05 0.07 0.09 0.15 0.2 0.25 0.5 1 2"  # GAMMA_LD in INCAR
gamma="0.1 0.04 0.02 0.01 0.008 0.006 0.004 0.001"  # GAMMA_LD in INCAR


#echo "this will create jobs with named:  low_2x2x2sc_280eV_02x02x02kp_EDIFF1E-3_tTiMESTEP\_gammaGAMMA"
#echo "KPOINTS: one kpoint (1 1 1 monkhorst pack) is enough since we just run on lambda 0.0 and an harmonic potential only"


[ ! -e "parameters.dat" ] && echo parameters.dat missing && exit -1
[ ! -e "KPOINTS" ] && echo KPOINTS missing && exit -1
[ ! -e "POTCAR" ] && echo POTCAR missing && exit -1
[ ! -e "INCAR" ] && echo INCAR missing && exit -1



rm -f jobList
for temp in $temps;do
for t in $time;do
echo $t
  for g in $gamma;do
    #[ "$t" = "10" ] && [ "$g" = "0.01" ] && echo had t10 g0.01!!!!!!!! && continue
    echo $g
   job=low_gammacheck_t$t\_gamma$g
   mkdir $job
   cd $job
   cp ../POTCAR .
   cp ../KPOINTS .
   cp ../INCAR .
   cp ../parameters.dat .
   sed -i 's|gamma=.*|gamma='"$g"'|' parameters.dat         # friction parameter of Langevin Thermostat (e.g. 0.01)
   sed -i 's|timestep=.*|timestep='"$t"'|' parameters.dat         # timestep of ionic motion in fs (e.g. 5 10)
   sed -i 's|ionicsteps=.*|ionicsteps=1|' parameters.dat         # how many ionic speps should be performed in md (e.g. 4000)
   sed -i 's|preequilibration=.*|preequilibration=100000|' parameters.dat         
   sed -i 's|temps=.*|temps='"$temp"'|' parameters.dat
   sed -i 's|kp=.*|kp=1 1 1|' parameters.dat
   sed -i 's|EDIFF=.*|EDIFF='"$EDIFF"'|' parameters.dat
   sed -i 's|NBANDS=.*|NBANDS='"$NBANDS"'|' parameters.dat

   ti_low_0_create_Folders_vasp.sh -c


    echo here we create the ti job
    cat jobList >> ../jobList
    cd ../
   done

echo
done
done
