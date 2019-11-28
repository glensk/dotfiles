#!/bin/bash


cC="0.005"      # C concentration in atomic fraction
rRepulsive=1.1 # radius (in units of a) around a C atom where no other C atom should be placed
sc=10          # supercell size
nAvg=100       # number of configurations to average over
nAvgStart=1
a=2.86         # lattice constant, only as start condition, will be relaxed

asc=`echo $sc $a | awk '{print $1*$2}'` # supercell lattice constant

dir=`pwd`
rm -f jobList
seed=`date +%s`
for c in $cC; do
  nFe=`echo 2 $sc | awk '{print $1*$2^3}'`
  nC=`echo $nFe $c | awk '{printf "%d",$1*$2/(1-$2)}'`  # nC = int(nFe * c/(1-c))
  nC=`echo $nC | awk '{print $1+3-($1%3)}'`            # we need nC to be dividable by 3 for an equal nr of C atoms on each sublattice
  cC=`echo $nC $nFe | awk '{print $1/($1+$2)}'`         # cC can change due to discrete nr of atoms
  nAt=`echo $nC $nFe | awk '{print $1+$2}'`

  mkdir -p Cconc_$cC
  cd Cconc_$cC
  echo Cconc_$cC

  # generate basis bcc lattice
  echo "" > coords
  echo "$nAt atoms" >> coords
  echo "2 atom types" >> coords
  echo "" >> coords
  echo "0 $asc xlo xhi" >> coords
  echo "0 $asc ylo yhi" >> coords
  echo "0 $asc zlo zhi" >> coords
  echo "0 0 0 xy xz yz" >> coords
  echo "" >> coords
  echo "Atoms" >> coords
  echo "" >> coords

  $dir/generateLattice.x $a $sc

  dir2=`pwd`
  for (( i=$nAvgStart; i<=`echo $nAvgStart $nAvg | awk '{print $1+$2}'`; i++ )); do
    mkdir -p dis_$i
    echo " dis_$i"
    cp coords dis_$i
    cp $dir/in.file_iso dis_$i/in.file
    cd dis_$i
    echo $seed > seed
    sleep 1 # to force new seed
    $dir/generateRandomC_on_sublattices.x $nC $rRepulsive $a $sc $seed 3
    seed=`expr $seed + 1`
    pwd >> $dir/jobList
    cd $dir2

    mkdir -p ord_$i
    echo " ord_$i"
    cp coords ord_$i
    cp $dir/in.file_aniso ord_$i/in.file
    cd ord_$i
    echo $seed > seed
    sleep 1 # to force new seed
    $dir/generateRandomC_on_sublattices.x $nC $rRepulsive $a $sc $seed 1
    seed=`expr $seed + 1`
    pwd >> $dir/jobList
    cd $dir2
  done

  cd $dir
done

