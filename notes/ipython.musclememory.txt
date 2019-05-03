  _i: stores previous input.
  _ii: next previous.
  _iii: next-next previous.
  _ih : a list of all input _ih[n] is the input from line n.
  _i<n> == _ih[<n>]

  _ (one underscore): previous output.
  __ (two underscores): next previous.
  ___ (three underscores): next-next previous.
 _<n> output n

  _i: stores previous input.
  _ii: next previous.
  _iii: next-next previous.
  _ih : a list of all input _ih[n] is the input from line n.
  _oh : list of all linex which created output


  _dh : show all visited paths

  edit -x eos.einet_   opens vi at function vinet_

!shellfunction			#to run a shell function
!!shellfunction			#to run a function and get output in list
or:
In [1]: %%bash
. ~/.bashrc
<my_fancy_bash_function> <function_argument>



who			shows all declared vars
whos		shows all declared vars and types
psearch v*  shows all variables starting with v
psearch v* int  shows all variables starting with v and are integer
store var   stores var for comming sessions
store var > ~/textfile sotres var to textfile
store	show all stored vars
reset		deletes all variables
store -r	resets stored vars
store -z	deletes storage
logstart	starts logging a session
lonon logoff halt/resume logging
lsmagic		shows magic command
magic		shows help about magic commands

import sys
print sys.path		in this paths python looks for modules
page sys.path		(list) in this paths python looks for modules 

function  4, 5	    == function(4, 5)
,function 4, 5		== function("4","5")
help('for') help(sys) shows help about for or sys command
help() then try topics

pdef re.match		shows what arguments re.match needs
pfile re.match		show the whole moduel code where re.match is defined
edit -x re.macht	open vi where re.match is defined
edit 1-4			opens vi and puts input lines 1-4 in it
edit xyz			directly opens vi to change xyz function (even if only defined in ipyther terminal)

psource xyz zeigt direkt 

cpaste				paset source from web directly into ipyton

bookmark name		bookmarks current path to name
bookmark -l			shows all avail bookmarks
cd name				brings one to the directory desired if it is in the bookmarks (tab completition works)

dhist				shows history of paths
lf					list files (only)
ldir				list directories (only)
alias ka echo you said :%s   makes an alias which takes argument

x = !ls
x			whoses the list
x.s			shows joined / connectes list
x.l			shows x as list (list with newlines)
x.n			whos seperated with \n

run				run script
run -d			run in debugger

_i1			what was the first command executed
In[3]		what was the third command executed
exec _i1	execute command nuber 1
exec In[1:4] execute commands 1 to 4

hist		shows the history
hist 60 07  shows corresponding history items

time xyz()	gives the time function xyz needs to execute
timeit xyz()	gives the (averaged) time function xyz needs to execute (runs serveral times)
prun xyz()	gives the time function xyz needs to execute of all steps in the function (list)

run -p file.py same as above for a file


########################################################
# importing modules
########################################################
module load utils
dir(utils)   shows all definitions in module utils
utils.__file__     shows the location of the .pyc file which was loaded on: module load utils


to masure time:

python -m cProfile -o prof /Users/glensk/Thermodynamics/python_thermodynamics/pot_energy_forces.py -dudlref -s 20
pstats.Stats('prof').strip_dirs().sort_stats("cumulative").print_stats(30)'''
