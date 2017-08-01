" vim:fdm=marker

" ================ Foldind (best without plugins) / Remember FOLDS & VIEW DEFINE AFTER SYNTAX ON ! (LOOK AND FEEL) ================== {{{
" to remember folds (next 2 lines)
autocmd BufWinLeave *.* mkview
autocmd BufWinEnter *.* silent loadview
" suggested by restore_view.vim plugin
set viewoptions=cursor,folds,slash,unix
"set foldmethod=indent   " foldmethod=marker   
set foldmethod=syntax    " now trying for c code
set foldmarker={,}
set foldnestmax=2
set foldlevel=99
" only fold if more than 3 ines are indented (2 are not enoug)
set foldminlines=3
" open folds when search or jump to specific line (:22)
"set foldopen=all
set foldopen=block,hor,insert,jump,mark,percent,quickfix,search,tag,undo
"set foldopen=jump
" foldcolumn will display a sidebar that indicates what will be folded.
" I stopped liking it, not useful to me
" set foldcolumn=2
" no markers like -------------------- to be in folds (default)
set fillchars="fold: " 
" ctermbg=0 (for iterm2) makes the folded line appear slightly brighter compared to normal background
" :hi Folded            returns            term=bold,underline cterm=bold,underline ctermfg=12 ctermbg=0 guifg=Cyan guibg=DarkGrey
":hi Folded term=NONE cterm=NONE gui=NONE ctermbg=None ctermfg=12 ctermbg=0 guifg=Cyan guibg=DarkGrey
":hi Folded term=NONE cterm=NONE gui=NONE ctermbg=None ctermbg=0 guibg=DarkGrey
hi Folded term=NONE cterm=NONE gui=NONE ctermbg=0
" often the text of the first line is not usefull but in .vimrc it is and can be
" sometimes in other files, therefore comment out next line
" "set foldtext="" 
" }}}

" ================ SEARCH ================== {{{
set smartcase  " if a search string contains UPPER letters than ignorecase is switched off
set hlsearch  " highlight all search matches  | to turn off set nohlsearch
" Press Space to turn off highlighting and clear any message already displayed.
nnoremap <silent> <Space> :nohlsearch<Bar>:echo<CR>
set ignorecase " If you dont care for case-sensitive searches and substitutions, turn it off completely
set incsearch  " Start searching before pressing enter.
" set showmatch " Show matching [] and {} by jumping to first one on closing .... not cool 
setlocal iskeyword-=:  " when searching for abc also find abc in abc:  from http://superuser.com/questions/552781/prevent-whole-word-search-from-matching-colon
" }}}

" ================ TAB behaviour | INDENTATION ================== {{{
" The next 3 options are crucial for python (to be in line with Pep8)
" Each tab that you type is converted to an equivalent number of spaces (Pep8)
set expandtab
" To control the number of space characters that will be inserted when the
" tab key is pressed set the 'tabstop' option: == a tab is four spaces
set tabstop=4
" To change the number of space characters inserted for indentation, use the 'shiftwidth' option:  
" If you prefer to work with spaces, then it is preferable to ensure that softtabstop == shiftwidth. 
" This way, you can expect the same number of spaces to be inserted whether you press the tab key in 
" insert mode, or use the indentation commands in normal/visual modes.
set shiftwidth=4
" If you prefer to work with tab characters then it is a good idea to ensure that tabstop == softtabstop. 
" This makes it less likely that youll end up with a mixture of tabs and spaces for indentation.
set softtabstop=4

" always set autoindenting on (if tab in oneline tab will be inserted in next line)
set autoindent  
" when press tab cursor moves to the appropriate location on the current
" line (not sorter)
set smarttab
" does the right thing (mostly) in programs
set smartindent 
set cindent     
" copy the previous indentation on autoindenting
set copyindent
" use multiple of shiftwidth when indenting with '<' and '>' shiftwidth, not
" tabstop
set shiftround
" }}}

" ================ further from sensible.vim (stuff most people use with vim) ================== {{{
set complete-=i
set nrformats-=octal
set ttimeout
set ttimeoutlen=1
set display+=lastline
set autoread
" }}}

" ================= vim annoyances  http://blog.sanctum.geek.nz/vim-annoyances/ ================== {{{


set nobackup 
set nowritebackup
" Backup files are a nuisance {{{
"   If you’re developing with a version control system, you might find
"   the in-place backups Vim keeps for saved files with the ~ suffix more
"   annoying than useful. You can turn them off completely with nobackup:
" }}}

set noswapfile
" Swap files are a nuisance {{{
"   Swap files can be similarly annoying, and unnecessary on systems with a
"   lot of memory. If you dont need them, you can turn them off
"   completely:
"   }}}

" backspace behaves normal (:help 'bs') bs=2 : same as :set backspace=indent,eol,start"
"   If youre in insert mode, by default you cant use backspace to
"   delete text past the start of the insert operation; that is, you
"   cant backspace over where you first pressed insert. This is old vi
"   behavior, and if you dont like it, you can make backspace work
"   everywhere instead:
set backspace=indent,eol,start

" Cant move into blank space in visual block mode
"   If you need to define a block in visual block mode with bounds outside the 
"   actual text (that is, past the end of lines), you can allow this with:
"   This will let you move around the screen freely while in visual block mode 
"   to define your selections. As an example, this can make selecting the contents 
"   of the last column in a formatted table much easier.
set virtualedit=block

" Filename completion on command line is confusing
"   If youe used to the behavior of shell autocomplete functions for completing 
"   filenames, you can emulate it in Vim’s command mode:
"   With this set, the first Tab press (or whatever your wildchar is set to) will 
"   expand a filename or command in command mode to the 
"   longest common string it can, and a second press will display a list of all 
"   possible completions above the command line.
" set wildmode=longest,list
set wildmode=longest,list,full

"Automatically cd into the directory that the file
set autochdir

" Disable error bells
set noerrorbells
" somehow set visualbell turnes visual bell off which was on by default
set visualbell 
set t_vb=

" Don’t reset cursor to start of line when moving around. (like gg or shift+G)
set nostartofline
" }}}

" Enable mouse in all modes  # I DONT WANT THIS ... WITHOUT mouse=a EVERYTHING {{{
" CAN BE COPIED :) ... but mouse=r or mouse=v might also work (better!?)
" it seems that everything can be compied also with mouse=a
if has("mouse")
set mouse=a
endif

" Use UTF-8 without BOM
set encoding=utf-8 nobomb
set termencoding=utf-8
set complete-=i
set ttyfast
set laststatus=2
set formatoptions=qrn1

" Show “invisible” characters
" set list turnes on list mode (but words are broken in the middle when line break)
" set list !! I DONT WATN THIS
"set listchars=tab:▸\ ,eol:¬
set listchars=nbsp:¬,tab:»·,trail:·


"hi User1 ctermfg=196 guifg=#eea040 guibg=#222222
"hi User2 ctermfg=75 guifg=#dd3333 guibg=#222222
"hi User3 guifg=#ff66ff guibg=#222222
"hi User4 ctermfg=239 guifg=#a0ee40 guibg=#222222
"hi User5 guifg=#eeee40 guibg=#222222

"Invisible character colors
"You can customise the syntax highlighting colours of invisible characters
"with the NonText and SpecialKey keywords.
"highlight NonText guifg=#4a4a59
"highlight SpecialKey guifg=#4a4a59
"highlight NonText guibg=#060606
"highlight Folded guibg=#0A0A0A guifg=#9090D0

set history=100         " remember more commands and search history
set undolevels=100      " use many muchos levels of undo
if has("wildmenu")
    set wildignore+=*.swp,*.bak,*.pyc,*.class
    set wildignore+=*.a,*.o
    set wildignore+=.DS_Store,.git,.hg,.svn
    set title                " change the terminal's title
    set novisualbell         " don't beep
    set noerrorbells         " don't bee
endif
"autocmd filetype python set expandtab

""""" folgende einstellungen sollen gut sein fuer python
"set modeline     "

" tabcompletion for vim
" Omni completion provides smart autocompletion for programs. When invoked, the
" text before the cursor is inspected to guess what might follow. A popup menu
" offers word completion choices that may include struct and class members,
" system functions, and more. A similar feature in Microsoft Visual Studio is
" known as IntelliSense


" }}}

" change cursorline in insert/command mode
autocmd InsertEnter * set cul
autocmd InsertLeave * set nocul

"hi textit term=NONE cterm=NONE gui=NONE ctermbg=None ctermfg=12 ctermbg=0 
"match textit /\\underline{\zs.\{-}\ze}/
"hi VisualNOS term=NONE cterm=NONE gui=NONE ctermbg=None ctermfg=12 ctermbg=0
"hi StatusLine term=NONE cterm=NONE gui=NONE ctermbg=None ctermfg=12 ctermbg=0
"nmap <silent> dsa ds}dF\
