#!/usr/bin/env sh
SETTING_FILE=$1
POSITION_FILE=$2
INPUT_FILE=$3

if [ "$#" -ne 3 ]; then
    USAGE="$0 SETTING_FILE POSITION_FILE INPUT_FILE"
    echo $USAGE
    exit
fi
if [ ! -e $SETTING_FILE ]; then
   echo "No setting file found at $SETTING_FILE"
   exit 1
fi
if [ ! -e $POSITION_FILE ]; then
   echo "No position file found at $POSITION_FILE"
   exit 1
fi
if [ -e $INPUT_FILE ]; then
   echo "$INPUT_FILE already exists!"
   exit 1
fi

#testing to make sure setting file is only settings
ISTHERE_CELLATOMIC_INFO=$(grep -iE "CELL_PARAMETERS|ATOMIC_POSITIONS" $SETTING_FILE)
if [ ! -z "$ISTHERE_CELLATOMIC_INFO" ]; then
   echo "$SETTING_FILE must not contain cell parameters or atomic positions!"
   exit 1
fi


cp $SETTING_FILE $INPUT_FILE

NATOM=$(grep nat $POSITION_FILE | awk {'print $3'})
grep CELL_PARAMETERS -A3 $POSITION_FILE >> $INPUT_FILE
grep ATOMIC_POSITIONS -A$NATOM $POSITION_FILE >> $INPUT_FILE

NATLINE=$(grep nat $POSITION_FILE)
sed -i "s/nat.*/$NATLINE/g" $INPUT_FILE

NTYPLINE=$(grep ntyp $POSITION_FILE)
sed -i "s/ntyp.*/$NTYPLINE/g" $INPUT_FILE
