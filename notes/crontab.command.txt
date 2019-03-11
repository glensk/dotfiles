crontab albert.cron		# start crontab albert.cron
crontab -r				#stop all crontabs

crontab -e			#edit the corntab, after closig it will be started
crontab -l 			#lists all cronjobs

to uncomment job in crontab:		just use # at the beginning of the line
after adding line to crontab:		stop crontab and start again to make line work

#*/1 * * * * /home/glensk/scripts/general/priority.check.sh  #repeats script every minute
