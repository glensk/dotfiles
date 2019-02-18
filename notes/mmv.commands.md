mmv -r "*.0529177_1-1-1" VOL=#1  				#Benennt alle um. #1 wird durch * ersetzt
mmv "*page*.htm*" "#1page_#2.htm#3" # ©2007 dsplabs.com.au"	#same with more wildcards
mmv "./foo_*/bar/*.txt" "./foo_#1/#2.txt" # ©2007 dsplabs.com.au"	#include entire path names
mmv -r "foo*/*.txt" "~#2.txt" # ©2007 dsplabs.com.au"		#works for subdirectories
mmv "./getdUdLs/dUdL_lambda*_all_*" "./getdUdLs/dUdL_lambda#1_all_#2_old"     # rename files in subdirectories (remove -r)

mmv -r "*" #u1		#make all uppercase
mmv -r "*" #l1		#make all lowercase

mmv [-m|r|c|o|a|l|s] [-h] [-v|n] [from to]

       -m : move source file to target name
       -r : rename source file or directory to target name
       -c : copy source file to target name
       -o : overwrite target name with source file
       -a : append contents of source file to target name
       -l : link  target name to source file
       -s : same as -l, but use symbolic links instead of hard links

       -v : verbose mode
       -n : no-execute mode

       -h : help information

mmv -r "?ALL_*" 0#1ALL_#2  	#just changes the ones with 1ALL_, 2ALL_,...
							#generally ?, *, [] are allowed
							#?? are just 2 wildcarts
mmv -r "?*.txt" "#1u1#2.txt 	makes first letter uppercase

mmv -r ";*VOL*" "#2Vol#3"		renames all subdirectories  #1 is here assigned to the first directory (;) of the subdirectories

mmv -r ";400eV_6x6x6kpm0_NGXF120_ED1E-3" "400eV_6x6x6kpm0_120NGXF_ED1E-3"       renames all subsubdirectories

mmv -r ";2x2x2sc*_vasp4" "2x2x2sc#2"  remames all subdirectories #1 is here assigned for the first directory (;) of the subdirectories
mmv -r ";3x3x3sc*_vasp4" "3x3x3sc#2"

mmv -r ";*dudl" "#2dUdL"

##################################################################################################
# At mac
##################################################################################################
mmv \*.new \#1   # rename all *.new files to *
mmv -r sum_all_new_\?\?\?.npy sum_all_new_\#1_\#2_\#3.npy    # first, second, thirc number/ character
mmv ps10_0_0_1200000.dat\* save.ps10_0_0_1200000.dat\#10_0_0_1200000.dat\* save.ps10_0_0_1200000.dat\#1