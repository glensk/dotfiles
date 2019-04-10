" python-mode /jedi-vim          ( {{{
Bundle 'davidhalter/jedi-vim'
" I don't want the docstring window to popup during completion
autocmd FileType python setlocal completeopt-=preview
let g:jedi#popup_on_dot = 0
let g:jedi#popup_select_first = 0
"let g:jedi#auto_vim_configuration = 0
""Finally, if you don't want completion, but all the other features, use:
"let g:jedi#completions_enabled = 0
"
""    Completion <C-Space>
""    Goto assignments <leader>g (typical goto function)
""    Goto definitions <leader>d (follow identifier as far as possible, includes imports and statements)
""    Show Documentation/Pydoc K (shows a popup with assignments)
""    Renaming <leader>r
""    Usages <leader>n (shows all the usages of a name)
""    Open module, e.g. :Pyimport os (opens the os module)
