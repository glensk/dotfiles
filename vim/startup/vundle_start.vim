" vundle                    (plug in manager) {{{
" set nocompatible is a must for vundle, solarized, ... and a lot of vi annoyances
set nocompatible              " be iMproved, required,  get rid of Vi compatibility mode. SET FIRST!
filetype off                  " required for vundle, will be enabled later on

" set the runtime path to include Vundle and initialize
set rtp+=~/.vim/bundle/Vundle.vim   
"set rtp+=~/.vim/bundle/vundle
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
