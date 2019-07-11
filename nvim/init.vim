"""""""""""""""""""""""""""""""""""""""""""""""""""""""
" get/install all plugins if some/all do not exist
"""""""""""""""""""""""""""""""""""""""""""""""""""""""
if empty(glob('$HOME/sources/nvim/autoload/plug.vim'))
  silent !curl -fLo $HOME/sources/nvim/autoload/plug.vim --create-dirs
    \ https://raw.githubusercontent.com/junegunn/vim-plug/master/plug.vim
  autocmd VimEnter * PlugInstall --sync | source $MYVIMRC
endif

"if empty(glob('$HOME/.local/share/nvim/site/autoload/plug.vim'))
"  silent !curl -fLo $HOME/.local/share/nvim/site/autoload/plug.vim --create-dirs https://raw.githubusercontent.com/junegunn/vim-plug/master/plug.vim
"  autocmd VimEnter * PlugInstall --sync | source $MYVIMRC
"endif

"function! PlugLoaded(name)
"    return (
"        \ has_key(g:plugs, a:name) &&
"        \ isdirectory(g:plugs[a:name].dir) &&
"        \ stridx(&rtp, g:plugs[a:name].dir) >= 0)
"endfunction

"""""""""""""""""""""""""""""""""""""""""""""""""""""""
" vim cheatsheat 
"""""""""""""""""""""""""""""""""""""""""""""""""""""""
" ~/.vim/vim_cheatsheet     " move, lowercase, search, indent, copy-paste, 

"""""""""""""""""""""""""""""""""""""""""""""""""""""""
" load Plugins (dont place comments after Plug command)
"""""""""""""""""""""""""""""""""""""""""""""""""""""""
call plug#begin('$HOME/sources/nvim/plugged')
" Plug 'tpope/vim-surround'                     " I dont seem to use it
Plug 'Yggdroot/indentLine'                      " shows ['|', '¦', '┆', '┊'] to indicate indents
Plug 'iCyMind/NeoSolarized'                     " color scheme
Plug 'majutsushi/tagbar'                        " show functions in current file: <leader>t
Plug 'ervandew/supertab'                        " tab completion for vim
Plug 'tomtom/tcomment_vim'                      " use gcc to un/comment line
Plug 'AndrewRadev/switch.vim'                   "  True -> False (simply press minus)
Plug 'jcf/vim-latex'                            " ,ll and ,ls in latex to compie stuff and open the pdf
Plug 'kshenoy/vim-signature'                    " set markers m[a-zA-Z]; delall marks :m<Space>
Plug 'easymotion/vim-easymotion'                " <Leader><Leader>{w,b}  :hi show all colors
Plug 'henrik/vim-indexed-search'                " results counter
"Plug 'scrooloose/nerdtree'                      " I dont use it
Plug 'junegunn/vim-easy-align'                  " use in a line gaip= to allign all next lines on =
"Plug 'tpope/vim-commentary'                    " Comments
"source ~/.vim/startup/vundle_neocomplete.vim   " breaks spell correction
call plug#end()

""""""""""""""""""""""""""""""""""""""""""
" settings      
""""""""""""""""""""""""""""""""""""""""""
source ~/.vim/startup/cmdline-complete.vim    " use ctrl+p in command mode (now tab)
source ~/.vim/startup/settings_plugins.vim   	" settings related to plugins
source ~/.vim/startup/keybindings.vim 	" keybindings for: navigation, folding, paste
source ~/.vim/startup/settings1.vim   	" memorize folding, wrapping of lines, persistent undo
source ~/.vim/startup/settings0.vim   	" general sets, remove trailing whitespaces, currently breaks yanking on mac!
source ~/.vim/startup/settings2.vim     " syntax! highlighting, templates, look-and-feel, this changes color of comments
source ~/.vim/startup/settings3.vim   " Folding, Search, Tab-behavior, indentation
set thesaurus+=$HOME/scripts/dotfiles/vim/thesaurus/mthesaur.txt  "Ctrl x + Ctrl t in insert mode
"""inoremap <expr><Tab>  pumvisible() ? "\<C-n>" : neocomplete#start_manual_complete()
source ~/.vim/startup/spellcorrection.vim  " Rechtschreibung,

set clipboard=unnamed           " to send copied stuff to system clipboard ( needs to set in ierm: keyboard shortcut: "cmd+c" -> Action: "send escape sequence" Esc+: "y"
set mouse=a                     " get scrolling inside nvim using mouse
" hi Comment  guifg=#80a0ff ctermfg=darkred           " Color for comments (red) " I like the actual gray better
