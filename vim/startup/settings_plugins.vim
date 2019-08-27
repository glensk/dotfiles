""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
" solarized
" Plug 'iCyMind/NeoSolarized'                     " color scheme
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
colorscheme NeoSolarized   	" set colorscheme, can only be done after plug#end
set termguicolors 		    " for colorscheme

""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
" tagbar            
" Plugin 'https://github.com/majutsushi/tagbar'
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
nnoremap <silent><leader>t :TagbarToggle<cr><C-w>l
nnoremap <leader>t :TagbarToggle<cr>
"nnoremap <silent><leader>t :TagbarOpenAutoClose<cr>
" autofocus on Tagbar open, autoclose on selection, show tag if folded, compact
let g:tagbar_autofocus = 1
"let g:tagbar_autoclose = 1
let g:tagbar_autoshowtag = 1
let g:tagbar_compact = 1
"let g:tagbar_sort = 0   " for 1 it sorts by name on 0 it sorts as in file


""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
" supertab 
" Plugin 'ervandew/supertab'
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
" expand tab        (Makes tab exansion possible, supertab) {{{
" NEVER SET: let g:SuperTabDefaultCompletionType = "context"  --> this destroys
" tabcompletion for class variables
" the following 3 lines we only need if we dont use jedi-vim
filetype plugin on
set omnifunc=syntaxcomplete#Complete
autocmd FileType python set omnifunc=pythoncomplete#Complete
set completeopt=menuone,longest,preview
" NEVER SET: let g:SuperTabDefaultCompletionType = "context"
" au FileType python set omnifunc=pythoncomplete#Complete
" let g:SuperTabDefaultCompletionType = "context"
" let g:SuperTabContextDefaultCompletionType = "<c-n>"


""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
" switch 
" Plugin 
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
nnoremap - :Switch<cr>


""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
" latex
" Plugin 'jcf/vim-latex.git'                          " latex
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
let g:Tex_IgnoredWarnings = 
    \'Underfull'."\n".
    \'Label(s) may have changed.'."\n".
    \'Overfull'."\n".
    \'specifier changed to'."\n".
    \'You have requested'."\n".
    \'Missing number, treated as zero.'."\n".
    \'A float is stuck'."\n".
    \'There were undefined references'."\n".
    \'Citation %.%# undefined'."\n".
    \'Float too large for'."\n".
    \'Float too large for page by'."\n".
    \'Double space found.'."\n"
let g:Tex_GotoError=0  " dont go to error but stay where you started compiling
let g:tex_flavor='latex'
let g:Tex_DefaultTargetFormat='pdf'
let g:Tex_ViewRule_pdf='Skim'
let g:Tex_SmartKeyBS = 0
let g:Tex_SmartKeyQuote = 0
let g:Tex_SmartKeyDot = 0
let g:Imap_UsePlaceHolders = 0 " dont add <++> which is a time saver and is removed with ctrl+J == ctrl+shift+j ... once ready this brings you you out of $ formula $ and deletes <++>
"let g:Tex_MultipleCompileFormats='dvi,pdf'
let g:Imap_FreezeImap=1


""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
" easymotin
" Plug 'easymotion/vim-easymotion'                    " <Leader><Leader>{w,b}  :hi show all colors
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
hi link EasyMotionTarget2First EasyMotionTargetDefault
hi EasyMotionTarget2Second ctermbg=none ctermfg=red guifg=#dc322f guibg=#002b36 gui=bold

""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
" easyalign
" Plug 'junegunn/vim-easy-align' ... use in a line gaip= to allign all next lines on = <Paste>
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
nmap ga <Plug>(EasyAlign)
xmap ga <Plug>(EasyAlign)
