" vim:fdm=marker
let mapleader = ","
let maplocalleader = ","

" Open new split panes to right and bottom, which feels more natural than Vim's
" default:
set splitbelow
set splitright

" Ctrl+[ works in general; 
" other way to avoid esc
imap jj <Esc>
inoremap jj <Esc>


" NAVIGATION moving in vim (general / splits / tabs /{{{


" general: go to next row instead of go to next line
nnoremap j gj
nnoremap k gk

" general: want usually the whole word
nmap w W
nmap b B

" scroll with control -> does not work!
"nnoremap <C-J> <C-E>
"map <C-J> <C-E>
"map <C-K> <C-Y>


" splits, move between splits (control||capslock + hjkl) ###############
nnoremap <C-h> <C-w>h
nnoremap <C-j> <C-w>j
nnoremap <C-k> <C-w>k
nnoremap <C-l> <C-w>l

"inoremap <C-j> <Left>
"inoremap <C-k> <Right>

" make hjkl movements accessible from insert mode via the <Alt> modifier key
inoremap <C-k> <Right>
inoremap <A-h> <C-o>h
inoremap <A-j> <C-o>j
inoremap <A-k> <C-o>k
inoremap <A-l> <C-o>l


"inoremap <M-h> <C-o>h   
"inoremap <D-h> <Left>
"inoremap '^[[D' <Left>
"inoremap <D-h> :echo 'yo'
"inoremap <M-h> <Left>
"
"
"inoremap <D-l> <C-o>l
"inoremap <D-l> <Right>
"inoremap <D-j> <C-o>j
" tabs
" open new tab, close tab, got tab to right/left
nmap tt <esc>:tabnew<cr>
" nmap ct <esc>:tabclose<cr>   NO chatge to  like ct) deletes everything to )
nmap tp <esc>:tabp<cr>
nmap t, <esc>:tabp<cr>
nmap th <esc>:tabp<cr>
nmap tn <esc>:tabn<cr>
nmap >  <esc>:tabn<cr>
nmap <  <esc>:tabp<cr>
nmap t. <esc>:tabn<cr>
nmap tl <esc>:tabn<cr>

" tabs: let label=v:lnum   " put number in tab
map t1 :tabn 1<CR>
map t2 :tabn 2<CR>
map t3 :tabn 3<CR>
map t4 :tabn 4<CR>
map t5 :tabn 5<CR>
map t6 :tabn 6<CR>
map t7 :tabn 7<CR>
map t8 :tabn 8<CR>
map t9 :tabn 9<CR>

" move fast
nnoremap <C-n> 10j
nnoremap <C-m> 10k

" }}}




" use ,ev to directly open .vimrc
"nnoremap <leader>ev <C-w><C-v><C-l>:e $MYVIMRC<cr>
nnoremap <leader>ev :sp $MYVIMRC<cr>

" use ,es to source vim directly
nnoremap <leader><leader>ev :so $MYVIMRC<cr>
nnoremap <leader>es :so $MYVIMRC<cr>

" folding
nmap FF zM
nmap ff za

nnoremap ; :
"nnoremap : ;
"nmap - :bd<CR>

" when deleting a line go to the biginning of line since it most time does not make
" sence to stat at the same spot
nmap dd dd^
" Fast saving
"nmap <leader>w :w!<cr>

nmap <silent> <Leader>u :nohlsearch<CR>
" make intend not loose highlighing
vnoremap > >gv
noremap < <gv


" swich between True / False
" nnoremap <leader>t :Switch<cr> 
nnoremap <leader>t :TagbarToggle<cr>

imap <D-v> ^O:set paste<Enter>^R+^O:set nopaste<Enter>
imap <D-V> ^O"+p
" <D-v> should mean cmd-v. The ^O and ^R are literal control-O and control-R, which
" you can type with ^V^O (control-v control-o) and ^V^R (control-v control-r).
" Control-O in insert mode allows you to execute one command then return to insert
" mode; here you can use it to put from the clipboard register.


" yank whole word -0.05 instead of just -
nnoremap yw yW

nnoremap <leader>u :redo<CR>

" Remap VIM 0 to first non-blank character
map 0 ^

" Move a line of text using ALT+[jk] or Comamnd+[jk] on mac
" nmap <M-j> mz:m+<cr>`z
" nmap <M-k> mz:m-2<cr>`z
" vmap <M-j> :m'>+<cr>`<my`>mzgv`yo`z
" vmap <M-k> :m'<-2<cr>`>my`<mzgv`yo`z

" if has("mac") || has("macunix")
"   nmap <D-j> <M-j>
"   nmap <D-k> <M-k>
"   vmap <D-j> <M-j>
"   vmap <D-k> <M-k>
" endif

" Toggle paste mode on and off
map <leader>pp :setlocal paste!<cr>

" Indent lines with cmd+[ and cmd+]
nmap <D-]> >>
nmap <D-[> <<
vmap <D-[> <gv
vmap <D-]> >gv
