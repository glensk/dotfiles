"Plugin 'jcf/vim-latex.git'                          " latex
Plug 'jcf/vim-latex'                          " latex
"Plugin 'https://github.com/vim-latex/vim-latex.git'  " latex
" Plugin 'gerw/vim-latex-suite.git'                   " seems to be newer
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

