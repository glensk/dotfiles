" expand tab        (Makes tab exansion possible, supertab) {{{
Plug 'ervandew/supertab'
" NEVER SET: let g:SuperTabDefaultCompletionType = "context"  --> this destroys
" tabcompletion for class variables
"
" the following 3 lines we only need if we dont use jedi-vim
filetype plugin on
set omnifunc=syntaxcomplete#Complete
autocmd FileType python set omnifunc=pythoncomplete#Complete

set completeopt=menuone,longest,preview
" NEVER SET: let g:SuperTabDefaultCompletionType = "context"
"
" au FileType python set omnifunc=pythoncomplete#Complete
" let g:SuperTabDefaultCompletionType = "context"
" let g:SuperTabContextDefaultCompletionType = "<c-n>"
