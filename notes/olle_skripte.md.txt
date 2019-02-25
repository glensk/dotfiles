
#################################################################
# olle stuff
#################################################################
[ "$printinfo" = "True" ] && echo "--> .bash_alias : olle stuff"
alias olle_1_process_outcar='$olletrunk/bin/process_outcar_4.6.py'
alias olle_2_workflow='cat $olletrunkshared/WORKFLOW'
alias olle_2_workflowvi='vi $olletrunkshared/WORKFLOW'
alias olle_3_getforceconstants100='getforceconstants.skript infile.ucposcar infile.ssposcar 100'
alias olle_3_getforceconstants3='getforceconstants.skript infile.ucposcar infile.ssposcar 3'
alias olle_4_extract_forceconstants='./extract_forceconstants'
alias olle_5_link='ln -s {out,in}file.forceconstant'
alias olle_6_hessematrix='$olletrunk/remap_forceconstant/remap_forceconstant'
alias process_outcar_4.6.py='/Users/glensk/Dropbox/scripts/codes_to_try/olle_trunk_savefiles/process_outcar_4.6.py'

