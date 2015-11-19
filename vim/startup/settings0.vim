" vim:fdm=marker
" ================ GENERAL set's ================

" This makes vim act like all other editors, buffers
" can exist in the background without being in a window.
" http://items.sjbach.com/319/configuring-vim-right
" you can have unwritten changes to a file and open a
" new file using :e, without being forced to write or undo your changes first.
set hidden


" Use the OS clipboard by default (on versions compiled with `+clipboard`)
" on cmpc cmmd vi/vim have -clipboard but gvim -v can be used which has +clipboard
" to make this work on mac (vim) and cmpc (gvim) and cmmc (gvim) have to use
" have to use unnamed,unnamedplus
set clipboard=unnamed,unnamedplus

" Enhance command-line completion / Better? completion on command line
set wildmenu

set esckeys   " use arrows in insertmode (for Blazej)

" ================ source vimrc directly:    :so after changing .vimrc {{{
autocmd! bufwritepost $HOME/.vimrc source %
" }}}

" ================ automatically remove all trailing whitespace {{{
" autocmd BufWritePre *.py :%s/\s\+$//e
" automaticaly remove all trailing whitespace and keep cursor position
fun! <SID>StripTrailingWhitespaces()
    let l = line(".")
    let c = col(".")
    %s/\s\+$//e
    call cursor(l, c)
endfun
autocmd FileType c,cpp,java,php,ruby,python autocmd BufWritePre <buffer> :call <SID>StripTrailingWhitespaces()
" }}}
