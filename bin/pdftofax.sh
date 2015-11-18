#!/bin/sh -e

OPTS=`getopt -o hbl -n pdftofax -- "$@"`
if [ $? != 0 ] ; then exit 1 ; fi
eval set -- "$OPTS"

pbmopts=""
rotate="cat"
res="200x98"
size="1728x1145"
echo 1
while true; do
    case "$1" in
        -h) cat <<EOF
Usage: pdftofax [-b] [-l] [-h] <[pdf-file] >[tiff-file]

Reads a PDF file on standard input and generates a multi-page tiff
file in Fax-G3 encoding on standard output with the correct
(asymetric) resolution.

    -h  Show this help message

    -b  Generate a non-dithered binary tiff. Use this for scanned
        documents

    -l  Landscape mode. The pages will be rotated by 90 degrees
        counter-clock-wise

EOF
            exit 0
            ;;
        -b)
            pbmopts="-threshold -value 0.9"
            shift
            ;;
        -l)
            rotate="pnmflip -ccw"
            res="98x200"
            size="1145x1728"
            shift
            ;;
        --)
            shift ; break ;;
        *)
            echo "Internal error"; exit 1 ;;
    esac
done

echo 2
temp=/tmp/pdftofax-$$
mkdir "$temp"
trap 'cd /; rm -rf "$temp"' 0 1 2 15
cd "$temp"

echo 3
echo 3 res $res
echo 3 size $size
which gs
gs -sDEVICE=pgmraw -r${res} -g${size} -sOutputFile="%d.pgm" -dSAFER -dNOPAUSE -dBATCH -dNOINTERPOLATE -q -


echo 4
for page in *.pgm; do
    pgmtopbm $pbmopts <$page \
        | ${rotate} \
        | pnmtotiff -g3 - >${page%.pgm}.tiff
done

tiffcp *.tiff fax.tiff
cat fax.tiff
