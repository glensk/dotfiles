" vim:fdm=marker
"
" ================ syntax for mathematica/html/python/... ================== {{{
au! BufRead,BufNewFile *.math       setfiletype mma
au! BufRead,BufNewFile *.MATH       setfiletype mma
au! BufRead,BufNewFile .bash_alias call SetFileTypeSH("bash")
au! BufRead,BufNewFile bash_set call SetFileTypeSH("bash")
"au! BufRead,BufNewFile *alias* call SetFileTypeSH("bash")
au! BufRead,BufNewFile generalrc call SetFileTypeSH("bash")
au! BufRead,BufNewFile *tcshrc call SetFileTypeSH("bash")
au! BufRead,BufNewFile .tcshrc call SetFileTypeSH("bash")
au BufNewFile,BufRead *.F set filetype=fortran95
autocmd BufRead *.htm,*.html set syntax=html
autocmd BufNewFile *.htm,*.html set syntax=html
"au BufNewFile,BufRead *.md,*.markdown,*.mdown,*.mkd,*.mkdn,README.md  setf markdown
"hi def htmlItalic              term=italic cterm=italic gui=italic

" }}}

" ================ TMPLATES (.sh, .py) ================== {{{
au BufNewFile *.py 0r ~/.vim/py.template
au BufNewFile *.sh 0r ~/.vim/sh.template
au BufWritePost,BufFilePost *.py call system("chmod +x ".expand("%"))
au BufWritePost,BufFilePost *.sh call system("chmod +x ".expand("%"))
" }}}


"================ LOOK ADN FEEL ================== {{{
set cursorline  " highlights the current line
set ruler       " line at the bottom of vi showing line number and char number

" Set 4 lines to the cursor - when moving vertically using j/k
set scrolloff=0  "Start scrolling when we're 4 lines away from margins before =4
set sidescrolloff=15
set sidescroll=1

if !&scrolloff
  set scrolloff=1
endif
if !&sidescrolloff
  set sidescrolloff=7
endif


" If you lose track of the current mode, you get a convenient --INSERT-- indicator
"set showmode
"Show partial/incomplete command in the last line of the screen when you type them
set showcmd


" switches in Insertmode between INSERT and INSERT paste
set pastetoggle=<F2>

" syntax highlighing: execute the command:source $VIMRUNTIME/syntax/syntax.vim
" also resets folding and hi Folded !!!!!!!!!!
"syntax on       " syntax highlighting dont do this here, breakes colorescheme


"" This currently makes not so nice colors in python comments, --> comment in!
"if filereadable($HOME.'/.vim/syntax/syntax.vim')
"    "echo "SpecificFile exists"
"    so ${HOME}/.vim/syntax/syntax.vim
"endif

" try to detect filetypes
filetype on         " try to detect filetypes
" enable loading indent file for filetype / necessary for file completition
filetype plugin indent on       " enable loading indent file for filetype/ ncesssary for file comletion


" to set number is always save, relativenumber breaks wraps in .tex files
set number 
if version >= 703
    " puts reltive numbers on left side
    " we will try this for higher version numbers currently it makes vim slow
    " down
    "set relativenumber  " BREAKS in .tex files the wrapping!
    set number    " shows the current line number instead of 0 
    " colorcolumn sets the marks the end of the document by a line
    set colorcolumn=91
else
    " puts linenumber on left sinde
    set number
endif
set relativenumber   " try now with neovim
" }}}
" necessary to handle correctly tex files
"set nolist wrap linebreak "breakat

