Plugin 'bling/vim-airline'
let g:airline#extensions#tabline#enabled = 0   " to show buffers or not
let g:airline#extensions#tabline#left_sep = ' '
let g:airline#extensions#tabline#left_alt_sep = '|'
let g:airline_powerline_fonts = 1
"let g:airline_theme             = 'powerlineish'
"let g:airline_theme             = 'badwolf'
let g:Powerline_symbols = 'fancy'
" let g:airline_enable_branch     = 1
let g:airline#extensions#branch#enabled = 1
let g:airline#extensions#syntastic#enabled = 1
" old let g:airline_enable_syntastic  = 1

" vim-powerline symbols
"let g:airline_powerline_fonts = 1
let g:airline_left_sep          = ''
let g:airline_left_alt_sep      = ''
let g:airline_right_sep         = ''
let g:airline_right_alt_sep     = ''
"""let g:airline_branch_prefix     = ''
"""let g:airline_readonly_symbol   = ''
"let g:airline_linecolumn_prefix = ''
""" let g:airline_symbols.linenr    = '-'
"let g:airline_section_c = '%F'   " shows full path
let g:airline_section_c = '%f'
let g:airline_section_x = ''   " get rid of filetype (python, vim, ...)
let g:airline_section_y = ''   " to get rid of utf-8[unix] == fileencoding
let g:airline_section_warning = '' " to get rid of trailing[51]
