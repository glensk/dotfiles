:xx,yy s/^/# /            vi uncomment line number xx bis yy
kommentieren / unkommentieren im visual mode: http://zinformatik.de/tipps-tricks/vim-tipps/schnelles-ein-und-auskommentieren-mit-vim/

:1,$s/Vol/VOL/g			substitutes every Vol to VOL
:%s/Vol/VOL/c           c == convirm
                        g == global




############################################
working with tabs
############################################

vim -p file1 file2 opens all file in vim with tabs

:n  go to next file
:prev go to previous file
:args show which files are open

:tabn nexttab
:tabp previoustab
:tabnew file to open additionally
:tabedit file (opens file in new tab)
gt when in normal mode to go to next tab
gT wechselt zum vorigen tab
:qa – schließt alle Tabs (ohne zu speichern) und vim
in ~/.vimrc: set showtabline=2 to show tabs always
nmap <C-H> :tabprev<CR>
nmap <C-L> :tabnext<CR>
:tabm moves tab to last spot
:tabm 2 moves tab to 3rd place (vim starts counting with 0)
:tabdo %s/foo/bar/g   That will run through each open tab and run the search and replace command (%s/foo/bar/g) in each one.

addint a tab: a) :tabnew b)go to new tab in vi c) vi --remote filenew

    "let name = fnamemodify(name,":t")   
    ":t nur filename
    ":p ganzer pfad
    ":r relativer pfad

############################################
makros aufnehmen
############################################
1.    Start recording by pressing q, followed by a lower case character to name the macro
2.    Perform any typical editing, actions inside Vim editor, which will be recorded
3.    Stop recording by pressing q
4.    Play the recorded macro by pressing @ followed by the macro name
5.    To repeat macros multiple times, press : NN @ macro name. NN is a number

############################################
open everything just in one vim instance
############################################

gvim --servername server1 ovim                              (open window)
gvim --servername server1 --remote-tab file1.txt file2.txt  (add files to window)       
gvim --serverlist

gvim --servername gvim --remote-tab-silent file1 file2  (better than above)



/\t              # suche nach tabs
/                # einleer nach / suche nach leerstellen

u               Undo
Ctrl-R          Redo
vim --version
VIM - Vi IMproved 7.3 (2010 Aug 15, compiled May  2 2013 15:16:05), Included patches: 1-244, 246-762

## open 2 files vertically splitted
vi -o /Users/glensk/proj/literature/play.bib /Users/glensk/proj/literature/play.ris
vi -o /Users/glensk/proj/literature/play.bib /Users/glensk/proj/literature/play.ris

            ahlloap



############################################
great shortcuts
############################################
in command mode:
0	beginning of line
$	end of line
shift i		first non whitespace char
shift a		append to end of line
o			create a new line below
shift o		create a new line above
shift c		change entire line
x			deletes
shift x		deletes back
.			repeates the last command


shift h		high
shift m		middle
shift l		low

shift j		fuege nexte zeile an die aktuelle
