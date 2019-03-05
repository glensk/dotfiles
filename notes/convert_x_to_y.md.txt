convert -density 500 file.pdf file.jpg

convert distribution_al-si_allpoints.eps distribution_al-si_allpoints.png
convert -geometry 110% distribution_al-si_allpoints.eps -flatten distribution_al-si_allpoints.jpg
gs -dNOPAUSE -r300 -sDEVICE=jpeg -sOutputFile=output.jpg distribution_al-si_allpoints.eps
gs -dEPSCrop -dNOPAUSE -r300 -sDEVICE=jpeg -sOutputFile=output.jpg distribution_al-si_allpoints.eps
gs -sDEVICE=pngalpha -dEPSCrop -dNOPAUSE -r300 -sOutputFile=output.png distribution_al-si_allpoints.eps

convert animation.gif target.png  (== to split animated gif in several png files)
convert -coalesce animation.gif target.png  (if gif has transparent areas)
convert image.png image.gif  (convert every png to gif)
gifsicle -U --disposal=previous --transparent="#ffffff" -O2 target1.gif > target01.gif   (remove background)
ffmpeg -i %03d.png output.gif   (join png files to single gif)

ffmpeg -i man2.mov -vcodec h264 -acodec aac -strict -2 man2.mp4
ffmpeg -i man2.mov -pix_fmt rgb24 output.gif     ## reduce size of *.giv:  convert -layers Optimize output.gif output_optimized.gif    
qlmanage -t -s 1000 -o . Thermal_conductivity.svg    # svg to png
