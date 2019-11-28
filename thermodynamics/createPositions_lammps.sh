#!/bin/sh

# following 3 lines must always be present
path=$(set | grep BASH_SOURCE | sed 's|.*\"\(.*\)/[^/]*\".*|\1|' | sed '/BASH_SOURCE=/s/.*/\./')
script=$(set | grep BASH_SOURCE | sed 's|.*\".*/\([^/]*\)\".*|\1|' | sed '/BASH_SOURCE=/s/.*/\./')
options=$*; . $path/utilities/functions.include; checkOptions "-h -o"

# small help for options
# space between option and value must be 1
# space between option and Comment or value and Comment must be >2
h=$(getOption -h)
if [ "$h" = "True" ]; then
  #usage $script STRUCTURE SC ALAT 
  echo 
  echo "How to use this script:"
  echo
  echo "$script STRUCTURE SC ALAT  "
  echo
  echo "e.g. $script fccqubic 3 4.13"
  echo
  echo "STRUCTURE: fccqubic"
  echo "           bccqubic"
  echo "           fccprimitive (currently not implemented)"
  echo
  echo "SC: supercell use integer (e.g. 2 for a 2x2x2sc)"
  echo 
  echo "ALAT: lattice constant in Angstrom (e.g. 4.13)"
  echo
  printOptions "-p         option not in use currently" \
               "-o         outputfilename; default = positions.dat"
  exit
fi


# number of atoms
struct=none
[ "$1" = "fccqubic" ] && struct=fccqubic
[ "$1" = "bccqubic" ] && struct=bccqubic
[ "$struct" = "none" ] && echo "struct ($1) not found" && exit

# sc size
sc=none
[ "$2" = "" ] && echo 'need $2 to be the sc size; I exit here' && exit
sc=$2

posfile=none
if [ "$struct" = "fccqubic" ];then
    posfile=$path/utilities/fcc/EqCoords_direct_fcc_$sc\x$sc\x$sc\sc
    if [ ! -e "$posfile" ];then
        posfilec=$path/utilities/fcc/coordinates_$sc\x$sc\x$sc\sc
        [ -e "$posfilec" ] && cat $posfilec | awk '{printf "%.10f   %.10f   %.10f\n", $1/'"$sc"',$2/'"$sc"',$3/'"$sc"'}' > $posfile
    fi
fi
if [ "$struct" = "bccqubic" ];then
    posfile=$path/utilities/bcc/EqCoords_direct_bcc_$sc\x$sc\x$sc\sc
    if [ ! -e "$posfile" ];then
        posfilec=$path/utilities/bcc/coordinates_$sc\x$sc\x$sc\sc
        [ -e "$posfilec" ] && cat $posfilec | awk '{printf "%.10f   %.10f   %.10f\n", $1/'"$sc"',$2/'"$sc"',$3/'"$sc"'}' > $posfile
    fi
fi



[ ! -e "$posfile" ] && echo "could not find posfile: $posfile" && exit

# alat
alat=none
[ "$3" = "" ] && echo 'need $3 to be alat; I exit here' && exit
alat=$3


# scalat
scalat=$(echo $alat $sc | awk '{print $1*$2}')
echo "scalat: $scalat"

# number of atoms
atoms=none
[ "$struct" = "fccqubic" ] && atoms=$(awk 'END { print NR }' $posfile)
[ "$struct" = "bccqubic" ] && atoms=$(awk 'END { print NR }' $posfile)
[ "$atoms" = "none" ] && echo "number of atoms not counted" && exit
echo "atoms: $atoms"

# filename
filename=positions.$alat.dat
[ "$(getOption -o)" = "True" ] && filename=$(getValue -o)
[ -e "$filename" ] && echo "filename $filename exists; I therefore exit here!" && exit
echo "filename: $filename"


# create File
touch "$filename"
echo "title placeholder" >> "$filename"
echo "     $atoms atoms" >> "$filename"
echo "       1 atom types" >> "$filename"
if [[ "$struct" = "fccqubic" || "$struct" = "bccqubic" ]];then
echo "    0.00000   $scalat xlo xhi" >> "$filename"
echo "    0.00000   $scalat ylo yhi" >> "$filename"
echo "    0.00000   $scalat zlo zhi" >> "$filename"
fi
echo "Atoms   " >> "$filename"
echo "   " >> "$filename"
awk '{print NR,1,$1*'"$alat"'*'"$sc"',$2*'"$alat"'*'"$sc"',$3*'"$alat"'*'"$sc"'}' "$posfile" >> "$filename"
