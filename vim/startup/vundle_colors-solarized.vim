" colors-soarized           (gives vim nice colors) {{{
Bundle 'altercation/vim-colors-solarized.git'
syntax enable

" works mac
"colors solarized
"colorscheme solarized
"set background=dark
"set t_Co=256


"" works for solarized light
"let g:solarized_termcolors=256 "this is what fixed it for me
"colorscheme solarized  "needs to be defined after let g:solarized_termcolors=256
"colors solarized       "if this is at the bottom, it initializes solarized-light colors with terminator 


" works for solarized dark but a bit wierd
"let g:solarized_termcolors=256 "this is what fixed it for me
"colorscheme solarized  "needs to be defined after let g:solarized_termcolors=256
"colors solarized       "if this is at the bottom, it initializes solarized-light colors with terminator 
"set background=dark

" works for solarized dark perfectly
set background=dark
"let g:solarized_termtrans = 1 " This gets rid of the grey background
colorscheme solarized
"let g:solarized_contrast="high"
"let g:solarized_visibility="high"


" }}}
