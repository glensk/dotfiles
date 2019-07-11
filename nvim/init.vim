" get plug if it does not exist
if empty(glob('$dotfiles/nvim/autoload/plug.vim'))
  silent !curl -fLo $dotfiles/nvim/autoload/plug.vim --create-dirs
    \ https://raw.githubusercontent.com/junegunn/vim-plug/master/plug.vim
  autocmd VimEnter * PlugInstall --sync | source $MYVIMRC
endif

"""""""""""""""""""""""""""""""""""""""""""""""""""""""
" load Plugins (dont place comments after Plug command)
"""""""""""""""""""""""""""""""""""""""""""""""""""""""
call plug#begin('$dotfiles/nvim/plugged')
Plug 'tpope/vim-surround'
Plug 'Yggdroot/indentLine'
Plug 'iCyMind/NeoSolarized'                     " color scheme
source ~/.vim/startup/vundle_tagbar.vim         " show all functions in current file: <leader>t
source ~/.vim/startup/vundle_supertab.vim       " tab completion in vim
Plug 'AndrewRadev/switch.vim'                   "  True -> False (simply press minus)
nnoremap - :Switch<cr>
source ~/.vim/startup/vundle_latex.vim          " ,ll and ,ls in latex
Plug 'kshenoy/vim-signature'                    " set markers m[a-zA-Z]; delall marks :m<Space>
source ~/.vim/startup/vundle_easymotion.vim     " <Leader><Leader>{w,b} :hi show all colors
Plug 'henrik/vim-indexed-search'                " results counter
"Plug 'scrooloose/nerdtree'                      " I dont use it
source ~/.vim/startup/cmdline-complete.vim      " use ctrl+p in command mode (now tab)
                                                " never change this to tab since this breaks
                                                " completion in :Plugin<TAB> command line mode
Plug 'junegunn/vim-easy-align'                  " use in a line gaip= to allign all next lines on =
nmap ga <Plug>(EasyAlign)
xmap ga <Plug>(EasyAlign)
"" Comments
"Plug 'tpope/vim-commentary'
call plug#end()

""""""""""""""""""""""""""""""""""""""""""
" clipboard
""""""""""""""""""""""""""""""""""""""""""
" to copy out of nvim: use ,yy
" to copy out of nvim: select with mouse and press y 
" to compy into nvim: cmd+v (in inster mode)
""""""""""""""""""""""""""""""""""""""""""
" settings      
""""""""""""""""""""""""""""""""""""""""""
" default value is "normal", Setting this option to "high" or "low" does use the 
" same Solarized palette but simply shifts some values up or down in order to 
" expand or compress the tonal range displayed.
let g:neosolarized_contrast = "normal"

" Special characters such as trailing whitespace, tabs, newlines, when displayed 
" using ":set list" can be set to one of three levels depending on your needs. 
" Default value is "normal". Provide "high" and "low" options.
let g:neosolarized_visibility = "normal"

" I make vertSplitBar a transparent background color. If you like the origin solarized vertSplitBar
" style more, set this value to 0.
let g:neosolarized_vertSplitBgTrans = 1

" If you wish to enable/disable NeoSolarized from displaying bold, underlined or italicized 
" typefaces, simply assign 1 or 0 to the appropriate variable. Default values:  
let g:neosolarized_bold = 1
let g:neosolarized_underline = 1
let g:neosolarized_italic = 0

colorscheme NeoSolarized   	" set colorscheme, can only be done after plug#end
set termguicolors 		" for colorscheme
source ~/.vim/startup/settings1.vim   	" memorize folding, wrapping of lines, persistent undo
source ~/.vim/startup/keybindings.vim 	" keybindings for: navigation, folding, paste
source ~/.vim/startup/settings0.vim   	" general sets, remove trailing whitespaces, currently breaks yanking on mac!
source ~/.vim/startup/settings2.vim     " syntax! highlighting, templates, look-and-feel, this changes color of comments
source ~/.vim/startup/settings3.vim   " Folding, Search, Tab-behavior, indentation
set thesaurus+=$HOME/scripts/dotfiles/vim/thesaurus/mthesaur.txt  "Ctrl x + Ctrl t in insert mode
""source ~/.vim/startup/vundle_neocomplete.vim          " breaks spell correction
"""inoremap <expr><Tab>  pumvisible() ? "\<C-n>" : neocomplete#start_manual_complete()
source ~/.vim/startup/spellcorrection.vim  " Rechtschreibung,

set clipboard=unnamed           " to send copied stuff to system clipboard ( needs to set in ierm: keyboard shortcut: "cmd+c" -> Action: "send escape sequence" Esc+: "y"
set mouse=a                     " get scrolling inside nvim using mouse
" config iTerm2 keys: Esc+Ac, Esc+As, Esc+Aa
"vnoremap <M-A>c "+y
"nnoremap <M-A>s :up<CR>
"inoremap <M-A>s <C-o>:up<CR>
"nnoremap <M-A>a ggVG
