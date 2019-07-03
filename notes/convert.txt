convert distribution_al-si_allpoints.eps distribution_al-si_allpoints.png
convert -geometry 110% distribution_al-si_allpoints.eps -flatten distribution_al-si_allpoints.jpg
gs -dNOPAUSE -r300 -sDEVICE=jpeg -sOutputFile=output.jpg distribution_al-si_allpoints.eps
gs -dEPSCrop -dNOPAUSE -r300 -sDEVICE=jpeg -sOutputFile=output.jpg distribution_al-si_allpoints.eps
gs -sDEVICE=pngalpha -dEPSCrop -dNOPAUSE -r300 -sOutputFile=output.png distribution_al-si_allpoints.eps
