date +%s    gives unixtime (only seconds)
date +%s%N  gives unixtime (nanoseconds) 
stat -c %Y file gives timestempt (unixtime) of last modification of the file

date +"%m-%d-%y"                --> 10-26-16
date +"%Y"                      --> 2016
date +"%y.%m.%d__%H:%M:%S__"    --> 16.12.22__18:02:19__
date +"%Y_%m_%d"                --> 2018_10_30