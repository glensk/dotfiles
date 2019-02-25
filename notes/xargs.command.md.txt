find -type d | xargs chmod +r
echo hallo wie gehts | xargs -n 1 echo  		# prints 1 word per line
echo hallo wie gehts | xargs -n 6 echo  		# prints 6 word per line

ls -1d *Ang | xargs -n 1 cp INCAR
find . -mindepth 2 -name INCAR | xargs -n 1 cp INCAR
