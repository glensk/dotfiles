" vim:fdm=marker
" to use this :set foldmethod=marker
" ctags (installed with homebrew)           {{{
" to update ctags file in the curret folder: in iterm: 
" /usr/local/Cellar/ctags/5.8/bin/ctags -R .   -> tags file
"
" }}}

" vundle                    (plug in manager) {{{
" set nocompatible is a must for vundle, solarized, ... and a lot of vi annoyances
set nocompatible              " be iMproved, required
filetype off                  " required

" set the runtime path to include Vundle and initialize
set rtp+=~/.vim/bundle/Vundle.vim
call vundle#begin()
Bundle 'gmarik/vundle'

    " Brief help
    " https://github.com/gmarik/Vundle.vim
    " :BundleList          - list configured bundles
    " :BundleInstall(!)    - install (update) bundles
    " :BundleSearch(!) foo - search (or refresh cache first) for foo
    " :BundleClean(!)      - confirm (or auto-approve) removal of unused bundles
    "
    " see :h vundle for more details or wiki for FAQ
    " NOTE: comments after Bundle commands are not allowed.
" }}}

" colors-soarized           (gives vim nice colors) {{{
Bundle 'altercation/vim-colors-solarized.git'
"
" }}}

" restore_view.vim          (restores view on files) {{{
Bundle 'vim-scripts/restore_view.vim.git'
"
" }}}

" syntax checker            (python syntax checker = pyflakes-vim) {{{
Bundle 'scrooloose/syntastic.git'
"
" }}}

" Style checker {{{
Bundle 'DamienCassou/textlint'
"Bundle 'https://github.com/vim-scripts/LanguageTool.git'
"let g:languagetool_jar='$HOME/Dropbox/scripts/dotfiles/bin/LanguageTool-2.6/languagetool-commandline.jar'
"
" }}}

syntax enable
colors solarized
colorscheme solarized
set background=dark
set t_Co=256

" surround (https://github.com/tpope/vim-surround.git) {{{
Bundle 'https://github.com/tpope/vim-surround.git'
"
" }}}

" mru
Bundle 'https://github.com/yegappan/mru.git'

" comment out
Bundle 'https://github.com/tomtom/tcomment_vim.git'

" nerdcommenter     (comment in/out lines: <leader>ci to toggel comand line) {{{
"Bundle 'scrooloose/nerdcommenter.git'
    " [count]<leader>cc |NERDComComment|
    " Comment out the current line or text selected in visual mode.
    " 
    " [count]<leader>cn |NERDComNestedComment|
    " Same as <leader>cc but forces nesting.
    " 
    " [count]<leader>c |NERDComToggleComment|
    " Toggles the comment state of the selected line(s). If the topmost selected line is commented, all selected lines are uncommented and vice versa.
    " 
    " [count]<leader>cm |NERDComMinimalComment|
    " Comments the given lines using only one set of multipart delimiters.
    " 
    " [count]<leader>ci |NERDComInvertComment|
    " Toggles the comment state of the selected line(s) individually.
    " 
    " [count]<leader>cs |NERDComSexyComment|
    " Comments out the selected lines ``sexily''
    " 
    " [count]<leader>cy |NERDComYankComment|
    " Same as <leader>cc except that the commented line(s) are yanked first.
    " 
    " <leader>c$ |NERDComEOLComment|
    " Comments the current line from the cursor to the end of line.
    " 
    " <leader>cA |NERDComAppendComment|
    " Adds comment delimiters to the end of line and goes into insert mode between them.
    " 
    " |NERDComInsertComment|
    " Adds comment delimiters at the current cursor position and inserts between. Disabled by default.
    " 
    " <leader>ca |NERDComAltDelim|
    " Switches to the alternative set of delimiters.
    " 
    " [count]<leader>cl
    " [count]<leader>cb |NERDComAlignedComment|
    " Same as |NERDComComment| except that the delimiters are aligned down the left side (<leader>cl) or both sides (<leader>cb).
    " 
    " [count]<leader>cu |NERDComUncommentLine|
    " Uncomments the selected line(s). 
" }}}

" ctrlP             (find files and code quickly files fuzzy finder; USE: ctrl + P) {{{
Bundle 'kien/ctrlp.vim'


let g:ctrlp_map = '<c-p>'
let g:ctrlp_cmd = 'CtrlP'
" ================ ignores when opening files / ctrlp(plugin) ================ {{{
set wildignore+=*/tmp/*,*/undodir/*,*.so,*.swp,*.zip,*.exe,*.pdf,*.ps     " MacOSX/Linux
set wildignore+=*/Documents/*

"let g:ctrlp_custom_ignore = '\v[\/]\.(git|hg|svn)$'
" Sane Ignore For ctrlp
"let g:ctrlp_custom_ignore = {
"  \ 'dir':  '\.git$\|\.hg$\|\.svn$\|\.yardoc\|public\/images\|public\/system\|data\|log\|tmp$',
"  \ 'file': '\.exe$\|\.so$\|\.dat$'
"  \ }
"let g:ctrlp_custom_ignore = '\v\.(exe|so|dll)$'
"let g:ctrlp_custom_ignore = {'dir':'\.git$\|\.hg$\|\.svn$\|\.yardoc\|public\/images\|public\/system\|data\|log\|tmp$', 'file': '\.exe$\|\.so$\|\.dat$'}
let g:ctrlp_custom_ignore = '\.git$\|\.hg$\|\.svn$\|\.yardoc\|public\/images\|public\/system\|data\|Library\|anaconda\|/undodir\|Movies\|/view\|log\|tmp$'
" }}}
" }}}

" switch True/False (Toggle True / False USE: - to toggle) {{{
Bundle 'https://github.com/AndrewRadev/switch.vim'
nnoremap - :Switch<cr>
" }}}

" vim-latex         (use <leader>ll to compie in vim, <leader>ls to jump to skim) {{{
Bundle 'jcf/vim-latex.git'
" REQUIRED. This makes vim invoke Latex-Suite when you open a tex file.
" filetype plugin on
"
" IMPORTANT: win32 users will need to have 'shellslash' set so that latex
" can be called correctly.
set shellslash
"
" IMPORTANT: grep will sometimes skip displaying the file name if you
" search in a singe file. This will confuse Latex-Suite. Set your grep
" program to always generate a file-name.
set grepprg=grep\ -nH\ $*
"
" OPTIONAL: Starting with Vim 7, the filetype of empty .tex files defaults to
" 'plaintex' instead of 'tex', which results in vim-latex not being loaded.
" The following changes the default filetype back to 'tex':
"set g:tex_flavor='latex'
"set g:Tex_DefaultTargetFormat='pdf'
"set g:Tex_ViewRule_pdf='Skim'
"set g:Tex_MultipleCompileFormats='dvi,pdf'

" }}}

" set marks         (m {a-z}: create mark '{a-z}: jump to mark: d{a-z} to delete mark) {{{
Bundle 'kshenoy/vim-signature.git'
"  m[a-zA-Z]    : Toggle mark
"  m<Space>     : Delete all marks
"  m,           : Place the next available mark
"  ]`           : Jump to next mark
"  [`           : Jump to prev mark
"  ]'           : Jump to start of next line containing a mark
"  ['           : Jump to start of prev line containing a mark
"  `]           : Jump by alphabetical order to next mark
"  `[           : Jump by alphabetical order to prev mark
"  ']           : Jump by alphabetical order to start of next line containing a mark
"  '[           : Jump by alphabetical order to start of prev line containing a mark
"
"  m[0-9]       : Toggle the corresponding marker !@#$%^&*()
"  m<S-[0-9]>   : Remove all markers of the same type
"  ]-           : Jump to next line having same marker
"  [-           : Jump to prev line having same marker
"  m<BackSpace> : Remove all markers
" }}}

" easymotion        (to navigate quickly using <leader><leader>w or b or f) {{{
Bundle 'Lokaltog/vim-easymotion'
" change the default EasyMotion shading to something more readable with Solarized
"hi link EasyMotionTarget ErrorMsg
"hi link EasyMotionShade  Comment
"hi link EasyMotionShade  ErrorMsg
" EasyMotionTarget2First
" EasyMotionTarget2Second
"let g:EasyMotion_do_shade = 0
" possible colors: Brown, red, white, Magenta, DarkRed
" show all colors with :hi
hi link EasyMotionTarget2First EasyMotionTargetDefault
"hi EasyMotionTarget2Second ctermbg=none ctermfg=red guifg=#000000 guibg=#cccccc gui=bold
hi EasyMotionTarget2Second ctermbg=none ctermfg=red guifg=#dc322f guibg=#002b36 gui=bold
" }}}

" tagbar            (Class/module browser Tagbar: <leader>t) {{{
Bundle 'https://github.com/majutsushi/tagbar'

nnoremap <silent><leader>t :TagbarToggle<cr><C-w>l
nnoremap <leader>t :TagbarToggle<cr>
"nnoremap <silent><leader>t :TagbarOpenAutoClose<cr>
" autofocus on Tagbar open, autoclose on selection, show tag if folded, compact
let g:tagbar_autofocus = 1
let g:tagbar_autoclose = 1
let g:tagbar_autoshowtag = 1
let g:tagbar_compact = 1
let g:tagbar_sort = 0   " for 1 it sorts by name on 0 it sorts as in file
" }}}

" results counter   (index-search) {{{
Bundle 'https://github.com/henrik/vim-indexed-search'
"
" }}}

" clipboard         (Provide pseudo clipboard register for non-GUI version of Vim) {{{
Bundle 'https://github.com/kana/vim-fakeclip'
"
" }}}

" expand tab        (Makes tab exansion possible, supertab) {{{
Bundle 'ervandew/supertab'
" NEVER SET: let g:SuperTabDefaultCompletionType = "context"  --> this destroys
" tabcompletion for class variables
"
" the following 3 lines we only need if we dont use jedi-vim
"filetype plugin on
"set omnifunc=syntaxcomplete#Complete
"autocmd FileType python set omnifunc=pythoncomplete#Complete

set completeopt=menuone,longest,preview
" NEVER SET: let g:SuperTabDefaultCompletionType = "context"
" }}}

" you complete me               {{{
"Bundle 'Valloric/YouCompleteMe'
""
" }}}

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


"}}}
"Plugin 'bling/vim-airline'
"set laststatus=2
"let g:airline_powerline_fonts = 1
"
"
" ################################### DEPRACTED PLUGINS #########################
" nerdtree          (very nice file browser) never used it PtclP is great {{{
"Bundle 'scrooloose/nerdtree.git'
"let NERDTreeChDirMode=2
"map <localleader>m :NERDTreeToggle %<CR>    " NerdtreeToggle
"map <localleader>n :NERDTreeToggle<CR>      " NerdtreeToggle
"map <leader>n :NERDTreeToggle<CR>
"let g:NERDTreeDirArrows=0                   " crucial to make it work
"nmap <esc>:nt<cr> :NERDTree<cr>B
"autocmd bufenter * if (winnr("$") == 1 && exists("b:NERDTreeType") && b:NERDTreeType == "primary") | q | endif "close vim if the only NERDTree is open
" }}}

" Efficient-python-folding          (DEPRACTED, destroys nesting) {{{
" nicely folds python stuff ... but does only to level folding (no nesting) but
" I want nesting
" Bundle 'vim-scripts/Efficient-python-folding.git'
" }}}

" vim-autoclose                     (DEPRACTED, destroys iabbrev) {{{
" Autoclose brackets emacs style (when writng and you write ) vim wont put it there
" This breaks iabbrev autocorrection
" Bundle 'Townk/vim-autoclose' " This breaks iabbrev autocorrection
" }}}

" Airline Lean & mean status/tabline (DEPRACTED) {{{
" Bundle 'bling/vim-airline'
" }}}

"" Ultimate auto-completion system for Vim. (Why dont I use it?){{{
"Bundle 'https://github.com/Shougo/neocomplcache.vim'
"" Note: This option must set it in .vimrc(_vimrc).  NOT IN .gvimrc(_gvimrc)!
"" Disable AutoComplPop.
"let g:acp_enableAtStartup = 0
"" Use neocomplcache.
"let g:neocomplcache_enable_at_startup = 1
"" " Use smartcase.
"let g:neocomplcache_enable_smart_case = 1
"" " Set minimum syntax keyword length.
"let g:neocomplcache_min_syntax_length = 3
"let g:neocomplcache_lock_buffer_name_pattern = '\*ku\*'
"
"
"" Enable heavy features.
"" Use camel case completion.
""let g:neocomplcache_enable_camel_case_completion = 1
"" Use underbar completion.
""let g:neocomplcache_enable_underbar_completion = 1
"
"" Define dictionary.
"let g:neocomplcache_dictionary_filetype_lists = {
"    \ 'default' : '',
"    \ 'vimshell' : $HOME.'/.vimshell_hist',
"    \ 'scheme' : $HOME.'/.gosh_completions'
"        \ }
"
"" Define keyword.
"if !exists('g:neocomplcache_keyword_patterns')
"    let g:neocomplcache_keyword_patterns = {}
"endif
"let g:neocomplcache_keyword_patterns['default'] = '\h\w*'
"
"" Plugin key-mappings.
"inoremap <expr><C-g>     neocomplcache#undo_completion()
"inoremap <expr><C-l>     neocomplcache#complete_common_string()
"
"" Recommended key-mappings.
"" <CR>: close popup and save indent.
"inoremap <silent> <CR> <C-r>=<SID>my_cr_function()<CR>
"function! s:my_cr_function()
"  return neocomplcache#smart_close_popup() . "\<CR>"
"  " For no inserting <CR> key.
"  "return pumvisible() ? neocomplcache#close_popup() : "\<CR>"
"endfunction
"" <TAB>: completion.
"inoremap <expr><TAB>  pumvisible() ? "\<C-n>" : "\<TAB>"
"" <C-h>, <BS>: close popup and delete backword char.
"inoremap <expr><C-h> neocomplcache#smart_close_popup()."\<C-h>"
"inoremap <expr><BS> neocomplcache#smart_close_popup()."\<C-h>"
"inoremap <expr><C-y>  neocomplcache#close_popup()
"inoremap <expr><C-e>  neocomplcache#cancel_popup()
"" Close popup by <Space>.
""inoremap <expr><Space> pumvisible() ? neocomplcache#close_popup() : "\<Space>"
"
"" Enable omni completion.
"autocmd FileType css setlocal omnifunc=csscomplete#CompleteCSS
"autocmd FileType html,markdown setlocal omnifunc=htmlcomplete#CompleteTags
"autocmd FileType javascript setlocal omnifunc=javascriptcomplete#CompleteJS
"autocmd FileType python setlocal omnifunc=pythoncomplete#Complete
"autocmd FileType xml setlocal omnifunc=xmlcomplete#CompleteTags



" }}}

" }}}
call vundle#end()            " required
filetype plugin indent on    " required

