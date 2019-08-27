""""""""""""""""""""""""""""""""""""""""""""""""""""""""
" get/install all plugins if some/all do not exist
"""""""""""""""""""""""""""""""""""""""""""""""""""""""
if empty(glob('$HOME/sources/nvim/autoload/plug.vim'))
  silent !curl -fLo $HOME/sources/nvim/autoload/plug.vim --create-dirs
    \ https://raw.githubusercontent.com/junegunn/vim-plug/master/plug.vim
  autocmd VimEnter * PlugInstall --sync | source $MYVIMRC
endif

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
"Plug 'farmergreg/vim-lastplace'                 " saves the exact position where file was opened last time
Plug 'vim-scripts/restore_view.vim'             " saves exact position and folds
Plug 'rickhowe/diffchar.vim'                    " show diffs (vim -d) for words and not for lines
"Plug 'tpope/vim-commentary'                    " Comments
"source ~/.vim/startup/vundle_neocomplete.vim   " breaks spell correction
call plug#end()

""""""""""""""""""""""""""""""""""""""""""
" update plug every week
""""""""""""""""""""""""""""""""""""""""""
function! OnVimEnter() abort
  " Run PlugUpdate every week automatically when entering Vim.
  if exists('g:plug_home')
    let l:filename = printf('%s/.vim_plug_update', g:plug_home)
    if filereadable(l:filename) == 0
      call writefile([], l:filename)
    endif

    let l:this_week = strftime('%Y_%V')
    let l:contents = readfile(l:filename)
    if index(l:contents, l:this_week) < 0
      call execute('PlugUpdate')
      call writefile([l:this_week], l:filename, 'a')
    endif
  endif
endfunction
autocmd VimEnter * call OnVimEnter()

""""""""""""""""""""""""""""""""""""""""""
" settings      
""""""""""""""""""""""""""""""""""""""""""
"let g:DiffUnit="Word1" " for diffchar.Vim
"let g:DiffUnit="Char" " for diffchar.Vim
source ~/.vim/startup/cmdline-complete.vim    " use ctrl+p in command mode (now tab)
source ~/.vim/startup/settings_plugins.vim   	" settings related to plugins
source ~/.vim/startup/keybindings.vim 	" keybindings for: navigation, folding, paste
source ~/.vim/startup/settings1.vim   	" memorize folding, wrapping of lines, persistent undo
source ~/.vim/startup/settings0.vim   	" general sets, remove trailing whitespaces, currently breaks yanking on mac!
source ~/.vim/startup/settings2.vim     " syntax! highlighting, templates, look-and-feel, this changes color of comments
source ~/.vim/startup/settings3.vim     " Folding, Search, Tab-behavior, indentation, mouse
set thesaurus+=$HOME/scripts/dotfiles/vim/thesaurus/mthesaur.txt  "Ctrl x + Ctrl t in insert mode
"""inoremap <expr><Tab>  pumvisible() ? "\<C-n>" : neocomplete#start_manual_complete()
source ~/.vim/startup/spellcorrection.vim  " Rechtschreibung,

" as it currently is (on mac) is best since can copy stuff in vim without the
" cursor changing (have removed the "Send escape sequence" Esc+: "y" in iterm key
" settings. set clipboard=unnamed is enabled. Scrolling works fine with the iterm
" command used: defaults write com.googlecode.iterm2 AlternateMouseScroll -bool true
" copy works (on mac) just by selecting the text (so cmd+c) can be pressed as well.
set clipboard=unnamed           " to send copied stuff to system clipboard ( needs to set in ierm: keyboard shortcut: "cmd+c" -> Action: "Send escape sequence" Esc+: "y"
" executed in iterm2: defaults write com.googlecode.iterm2 AlternateMouseScroll -bool true to get mouse scrolling
" hi Comment  guifg=#80a0ff ctermfg=darkred           " Color for comments (red) " I like the actual gray better

" Uncomment the following to have Vim jump to the last position when reopening
" a file (only line, not exact position) 
"if has("autocmd")
"  au BufReadPost * if line("'\"") > 1 && line("'\"") <= line("$") | exe "normal! g'\"" | endif
"endif
