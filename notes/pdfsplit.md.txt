# how to execute it:
# cmd + space to open alfred
# 'filename
# cmd + y
# pdfsplit{1,2}
#
# or
# select file in finder
# cmd + y
# pdfsplit{1,2}
 

folder=`echo "{query}" | sed -n 's|\(.*\)/.*|\1|p'`
file=`basename "{query}"`
fileopdf=`echo $file | sed 's|.pdf||'`
time=`date +%s`
cd $folder
#/usr/local/bin/pdftk "{query}" burst
#$HOME/Dropbox/scripts/dotfiles/cpdf-binaries/OSX-Intel/cpdf "{query}" -split
#-chunk 1 -o $fileopdf\_$time\_%%%.pdf



folder=`echo "{query}" | sed -n 's|\(.*\)/.*|\1|p'`
file=`basename "{query}"`
fileopdf=`echo $file | sed 's|.pdf||'`
time=`date +%s`
cd $folder
#/usr/local/bin/pdftk "{query}" burst
#$HOME/Dropbox/scripts/dotfiles/cpdf-binaries/OSX-Intel/cpdf "{query}" -split
#-chunk 2 -o $fileopdf\_$time\_%%%.pdf
