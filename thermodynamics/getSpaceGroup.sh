#!/bin/bash

#-----set parameters and paths------------------------
inpFile=POSCAR; tolDef=0.0001
#-----------------------------------------------------


# following 3 lines must always be present
path=`set | grep BASH_SOURCE | sed 's|.*\"\(.*\)/[^/]*\".*|\1|' | sed '/BASH_SOURCE=/s/.*/\./'`
script=`set | grep BASH_SOURCE | sed 's|.*\".*/\([^/]*\)\".*|\1|' | sed '/BASH_SOURCE=/s/.*/\./'`
options=$*; . $path/utilities/functions.include; checkOptions "-h -help -P -i -c -f -t"

# small help for options
# space between option and value must be 1
# space between option and Comment or value and Comment must be >2
h=`getOption -h`
if [ $h == True ]; then
  usage $script
  printOptions "-P POS   use file POS as input POSCAR file (default: $inpFile)" \
               "-i POS   same as -P (for compatibility with other scripts " \
               "-c       use cell, cartesian_coords or reduced_coords, and species files as input" \
               "-f       full output" \
               "-t TOL   use TOL as tolerance (default: $tolDef)"
  exit
fi

# detailed help
help=`getOption -help`
if [ $help == True ]; then
  details $script
  echo2 "   script is a wrapper for the sgroup program determining the space group"   \
        "   of an atomic structure in a POSCAR format"
  echo2 "   a cell, cartesian_coords or reduced_coords, and species files can be"     \
        "   used as an alternative input (-c option); species contains as many lines" \
        "   as the coordinates file (cartesian_coords or reduced_coords) and on each" \
        "   the name of the corresponding atom (numbers are also possible)"
  echo2 "   the wrapper produces an appropriate input for sgroup which determines"    \
        "   the primitive cell, the Bravais lattice, space group, irreducible atoms"  \
        "   and symmetry operations"
  echo2 "   the wrapper additionally tries to determine the prototype and the"        \
        "   Strukturbericht of the structure by comparing to the list in"             \
        "   $path/utilities/SpaceGroupDefinitions.dat"                                \
        "   obtained from cst-www.nrl.navy.mil"
  echo2 "   by default the wrapper produces a short output with the most imortant"    \
        "   information (log_space); with the -f option the full output is kept:"     \
        "     \033[1msgroup_output\033[0m: original sgroup output file"               \
        "     \033[1mirreducible_mapping\033[0m: mapping of the symmetry equivalent atoms onto" \
        "       each other, each line containing the indices of a symmetry"                     \
        "       equivalent set of atoms with the indices corresponding to the order"            \
        "       in the POSCAR file (or the coordinates files)"                                  \
        "     \033[1msymOps\033[0m: all symmetry matrices of the atomic structure in Cartesian" \
        "       coordinates; a matrix is given on three lines in three columns"                 \
        "     \033[1mshiftVecs\033[0m: shift vectors corresponding to the symOps file in"       \
        "       in Cartesian coordinates; one vector per line"
  exit
fi

# check if we are using cell and coordinates files as input (-c option)
cellinp=`getOption -c`
if [ $cellinp == True ]; then
  check cell
  if [ -e cartesian_coords ]; then
    if [ -e reduced_coords ]; then error "provide only either cartesian_coords OR reduced_coords"; fi
    $path/fortran/toReduced.x > /dev/null; rm cartesian_coords
  else
    if [ ! -e reduced_coords ]; then error "-c option used but no cartesian_coords nor reduced_coords available"; fi
  fi
  if [ ! -e species ]; then
    # if no species file available assume only one species
    sed 's/.*/1/' reduced_coords > species
  fi
else
  pos=`getOption -P`
  inp=`getOption -i`
  if [ $pos == True -a $inp == True ]; then error "-P and -i cannot be used together"; fi
  if [ $pos == True ]; then
    pos=`getValue -P`
    if [ "$pos" == "" ]; then error "no value to -P option provided"; fi
  else
    pos=$inpFile
  fi
  if [ $inp == True ]; then
    pos=`getValue -i`
    if [ "$pos" == "" ]; then error "no value to -i option provided"; fi
  else
    pos=$inpFile
  fi

  check $pos $path/splitPOSCAR.sh
  $path/splitPOSCAR.sh -i $pos -t reduced > /dev/null
fi

# start preparing input file for sgroup
# first: Bravais lattice (always primitive), cell vector lengths and angles
echo P > _tmp_input
awk 'function norm(x,y,z) { return (x^2+y^2+z^2)^(1/2) }
     function dot(x,y,z,X,Y,Z) { return x*X+y*Y+z*Z }
     function acos(x) { return atan2(sqrt(1-x*x), x) }
     BEGIN{deg=180/2/acos(0)}
     NR==1{a1=$1;a2=$2;a3=$3; na=norm(a1,a2,a3)}
     NR==2{b1=$1;b2=$2;b3=$3; nb=norm(b1,b2,b3)}
     NR==3{c1=$1;c2=$2;c3=$3; nc=norm(c1,c2,c3)}
     END{ alpha=deg*acos(dot(b1,b2,b3,c1,c2,c3)/norm(b1,b2,b3)/norm(c1,c2,c3));
          beta =deg*acos(dot(a1,a2,a3,c1,c2,c3)/norm(a1,a2,a3)/norm(c1,c2,c3));
          gamma=deg*acos(dot(a1,a2,a3,b1,b2,b3)/norm(a1,a2,a3)/norm(b1,b2,b3));
          printf(" %.10f %.10f %.10f  %.10f %.10f %.10f\n",na,nb,nc,alpha,beta,gamma)}' cell > _tmp_lengths_angles
cat _tmp_lengths_angles >> _tmp_input

# now nr of atoms
nAt=`awk 'END{print NR}' species`
echo >> _tmp_input; echo $nAt>> _tmp_input; echo >> _tmp_input

# add coordinates and element names
paste reduced_coords species | awk '{print $1,$2,$3; print $4}' >> _tmp_input

# get tolerance
if [ `getOption -t` == True ]; then
  tol=`getValue -t`
  c=`checkReal $tol`; if [ "$c" != ok ]; then error "tolerance value from -t option wrong"; fi
else
  tol=$tolDef
fi

# run sgroup
check $path/sgroup/sgroup
$path/sgroup/sgroup -set-TOL=$tol -prim _tmp_input > sgroup_output

# print info for original cell
echo > log_space; echo "Original cell:" >> log_space
awk '{printf(" Vector lengths:   %.2f  %.2f  %.2f\n",$1,$2,$3); printf(" Vector angles:    %.2f %.2f %.2f\n",$4,$5,$6)}' _tmp_lengths_angles >> log_space
ns=`cat species | sort -u | wc -l`
echo " Nr of species:    $ns" >> log_space
echo " Nr of all atoms:  $nAt" >> log_space
stoichio=`awk 'BEGIN{n=0};{ins=0; for (i=1;i<=n;i++) {if (a[i]==$1) {ins=1;c[i]=c[i]+1}} if (ins==0) {n=n+1;a[n]=$1;c[n]=1}};
              END{for (i=1;i<=n;i++) print c[i]}' species | sort -g | xargs | sed 's/ /:/g'`
echo " Stoichiometry:    $stoichio" >> log_space

# grep relevant info for log
echo >> log_space; echo "Primitive cell:" >> log_space
sed -n -e '/Bravais lattice/s/.*:\(.*\)/ Bravais lattice: \1/p' \
       -e '/     a    /{n;s/^ *\([^ ]\+\...\)[^ ]\+ \+\([^ ]\+\...\)[^ ]\+ \+\([^ ]\+\...\)[^ ]\+ */ Vector lengths:   \1   \2   \3/p}' \
       -e '/     alpha/{n;s/^ *\([^ ]\+\...\)[^ ]\+ \+\([^ ]\+\...\)[^ ]\+ \+\([^ ]\+\...\)[^ ]\+ */ Vector angles:    \1  \2  \3/p}' \
       -e '/Number of atoms in cell:/s/==== .*cell:\(.*\)/ Nr of all atoms: \1/p' \
       -e '/Number of nonequivalent sorts/s/.*sorts:\(.*\)/ Nr of irr atoms: \1/p' \
       -e '/Number and name of space group/s/.*group:\(.*\)/ Space group:     \1/p' \
       -e '/names of point group/{n;s/^ .*  \+\(.*\)/ Point group(Sf):  \1/p}' sgroup_output >> log_space

# get prototype and Strukturbericht for space group
nSp=`grep "Space group:" log_space | awk '{print $3}'`
nIrr=`grep "Nr of irr atoms" log_space | awk '{print $NF}'`
check $path/utilities/SpaceGroupDefinitions.dat
awk 'NR>2&&$7=='$nSp'{if('$ns'!=$2||'$nIrr'!=$(NF-1)) next;
                      n1=split($3,a1,":"); n2=split("'$stoichio'",a2,":");
                      if (n1!=n2) next;
                      r=a2[1]/a1[1];
                      for (i=2;i<=n1;i++) if (r!=a2[i]/a1[i]) next;
                      print}' $path/utilities/SpaceGroupDefinitions.dat > _tmp_prototype
echo " Prototype:        `awk '{printf("%s  ",$1)};END{printf("\n")}' _tmp_prototype`" >> log_space
echo " Strukturbericht:  `awk '{printf("%s  ",$5)};END{printf("\n")}' _tmp_prototype`" >> log_space
echo " Notes:            `awk '{printf("%s  ",$NF)};END{printf("\n")}' _tmp_prototype`" >> log_space

# extract nr of symOps, also nr of symmorphic and non-symmorphic operations
awk '/Number of symmetry operations/{n=$NF; nn=0; printf(" Nr of symOps:     %s\n",n)}
     /^Operation:/{c=0;next}
     c==0&&NF==4&&$4!=0{c=1;nn=nn+1}
     END{printf(" Symmorphic:       %d\n",n-nn); printf(" Non-symmorph:     %d\n",nn)}' sgroup_output  >> log_space

cat log_space

rm cell reduced_coords scale species _tmp_input _tmp_lengths_angles _tmp_prototype

# we are done if -f option not given
full=`getOption -f`
if [ $full != True ]; then rm sgroup_output; exit; fi



