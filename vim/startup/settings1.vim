" vim:fdm=marker
" ================ wrapping of lines (and words) {{{
set wrap                    " The next section makes Vim handle long lines correctly:
" set showbreak=...         " to make wrappled lines indent nicly; DONT WANT THIS
set textwidth=83            " Wrap text after a certain number of characters python: 79
set linebreak               " By default setting wrap breaks words in the middle, linebreak avoids this 
"set nolist                  " list disables linebreak
set nolist wrap linebreak   " necessary to handle correctly tex files; excluded breakat here, was not wroking
" }}}


if !isdirectory($HOME."/.vimviewdir")
    call mkdir($HOME."/.vimviewdir", "p")
endif

" DESTROYS WRAPS ON TEX FILES
set viewdir=$HOME/.vimviewdir
"au BufWinLeave * silent! mkview
"au BufWinEnter * silent! loadview
" THIS DOES NOT DESTROYS WRAP on yy! so we can savely use it
"au BufWinLeave *.tex mkview  " DO NOT USE DESTROYS WRAPS
"au VimEnter *.tex loadview   " DO NOT USE DESTROYS WRAPS

"au BufWinLeave * mkview  " DO NOT USE DESTROYS WRAPS
"au VimEnter * loadview   " DO NOT USE DESTROYS WRAPS

" ================ Persistent Undo across sessions ================== {{{
" Starting from vim 7.3 undo can be persisted across sessions
" http://www.reddit.com/r/vim/comments/kz84u/what_are_some_simple_yet_mindblowing_tweaks_to/c2onmqe
if has("persistent_undo")
    set undodir=~/.vimundodir
    set undofile   "Maintain undo history between sessions
    if !isdirectory($HOME."/.vimundodir")
        call mkdir($HOME."/.vimundodir", "p")
    endif
endif
" }}}
