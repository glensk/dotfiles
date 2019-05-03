iabbrev pritn print
iabbrev pirnt print
iabbrev prit print
iabbrev paht path
iabbrev sefl self
iabbrev retrun return
iabbrev Flase False
iabbrev Treu True 
" turn spell checking on
"set spell spelllang=de,en_us
set spell spelllang=en_us
autocmd BufRead,BufNewFile *.tex,*.txt setlocal spell spelllang=en_us

" to turn off  :set nospell

highlight clear SpellBad
highlight SpellBad term=standout ctermfg=4 term=underline cterm=underline
hi        SpellBad cterm=undercurl ctermfg=4
highlight clear SpellCap
highlight SpellCap term=underline cterm=underline
highlight clear SpellRare
highlight SpellRare term=underline cterm=underline
highlight clear SpellLocal
highlight SpellLocal term=underline cterm=underline


" \(\<\w\+\>\)\_s*\1  findes duplicate words with search
" 
" syn match TmlDoubleWords /\c\<\(\S\+\)\_s\+\1\ze\_s/    correct syntax but not working
" syn match TmlDoubleWords /\c\<\(\S\+\)\_s\+\1\ze\_s/
" syn match TmlDoubleWords /\<\(\w\+\)\s\+\1\> does not work
" syn match TmlDoubleWords 
" hi def link TmlDoubleWords ToDo 

"syn match TmlDoubleWords \(\<\w\+\>\)\_s*\1 
"hi def link TmlDoubleWords ToDo


" ANLEITUNG:
" ]s go to next     misspelled word
" [s go to previous misspelled word
" zg to add a word to dictionary: zg
" zug   Remove word from regular dictionary
" for words with umlaut: go to visual mode, mark word, zg
" z=, --> Once the cursor is on the word, use z=, and Vim will suggest alternatives 
" go to http://vimdoc.sourceforge.net/htmldoc/spell.html for more information
