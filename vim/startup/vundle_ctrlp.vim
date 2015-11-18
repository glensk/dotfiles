Plugin 'kien/ctrlp.vim'
let g:ctrlp_map = '<c-p>'
let g:ctrlp_cmd = 'CtrlP'
let g:ctrlp_max_height = 30
" ignores when opening files / ctrlp(plugin) ================ {{{
"set wildignore+=*/tmp/*,*/undodir/*,*.so,*.swp,*.zip,*.exe,*.pdf,*.ps     " MacOSX/Linux
"set wildignore+=*/Documents/*
"set wildignore+=*.pdf
"set wildignore+=*.ps
"set wildignore+=*.pyc
set wildignore+=*/tmp/*,*.so,*.swp,*.zip,*.pdf,*.ps,*.nb,*.djvu,*.cdf,*.docx,*.dmg,*.tar.gz,*.png     " MacOSX/Linux
"let g:ctrlp_custom_ignore = {'dir':'\.git$\|\.hg$\|\.svn$\|\.yardoc\|public\/images\|public\/system\|data\|log\|tmp$', 'file': '\.exe$\|\.so$\|\.dat$'}
let g:ctrlp_custom_ignore = '\.git$\|\.hg$\|\.svn$\|\.yardoc\|public\/images\|public\/system\|data\|Library\|anaconda\|/undodir\|Movies\|/view\|log\|tmp$'
