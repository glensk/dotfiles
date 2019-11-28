#!/bin/bash

error() { echo; echo -e "\033[31m\033[1mERROR:\033[0m $1"; exit
}

# check input
if [ $# != 1 ]; then error "provide # of steps for new memory"; fi
nSteps=$1

# check files
if [ ! -e MEMORYFILE ]; then error "no MEMORYFILE"; fi
if [ ! -e freq_history ]; then error "no freq_history"; fi

# check if enough steps
n=`awk 'END{print NR}' freq_history`
if [ $n -lt $nSteps ]; then error "not sufficient steps in freq_history"; fi

# check consistency between MEMORYFILE and freq_history
n1=`awk 'BEGIN{f=0};f==0&&NF==2{print NR-2;f=1}' MEMORYFILE`
n2=`tail -n1 freq_history | xargs -n1 | awk 'END{print NR-1}'`
if [ $n1 != $n2 ]; then error "\# of frequencies inconsistent"; fi

i1=`head -n1 MEMORYFILE`
i2=`awk 'END{print $1+1}' freq_history`
if [ $i1 != $i2 ]; then error "iterations inconsistent"; fi

# construct new memory
echo $i1 | awk '{printf("%4d\n",$1)}' > tmp_MEM
tail -n$nSteps freq_history | \
  awk 'BEGIN{for (i=2;i<=NF;i++) a[i]=0}
            {for (i=2;i<=NF;i++) a[i]=a[i]+$i}
         END{for (i=2;i<=NF;i++) printf("%30.20f\n",a[i]/'$nSteps')}' >> tmp_MEM
awk "NR>$n1+1{print}" MEMORYFILE >> tmp_MEM
mv MEMORYFILE MEMORYFILE_old
mv tmp_MEM MEMORYFILE

