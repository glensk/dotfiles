" tagbar            (Class/module browser Tagbar: <leader>t) {{{
Plugin 'https://github.com/majutsushi/tagbar'

nnoremap <silent><leader>t :TagbarToggle<cr><C-w>l
nnoremap <leader>t :TagbarToggle<cr>
"nnoremap <silent><leader>t :TagbarOpenAutoClose<cr>
" autofocus on Tagbar open, autoclose on selection, show tag if folded, compact
let g:tagbar_autofocus = 1
let g:tagbar_autoclose = 1
let g:tagbar_autoshowtag = 1
let g:tagbar_compact = 1
let g:tagbar_sort = 0   " for 1 it sorts by name on 0 it sorts as in file
