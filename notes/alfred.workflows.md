folder=`echo "{query}" | sed -n 's|\(.*\)/.*|\1|p'`
time=`date +%s`
cd $folder
#/usr/local/bin/pdftk "{query}" burst
$HOME/Dropbox/scripts/dotfiles/bin/cpdf-binaries/OSX-Intel/cpdf "{query}" -split -chunk 2 -o $time\_out%%%.pdf
