#!/usr/bin/env python
import os,sys
from subprocess import call

dotfiles = os.environ['dotfiles']

def gitclone(what):
    call(["git","clone","--depth","1",what])
    return

def git_getfoldername(link):
    out = link.split('/')[-1]
    out = out.split('.git')[0]
    #print(out)
    return out

def file_to_folder(filepath):
    folder = "/".join(filepath.split('/')[:-1])
    return folder


bundle = dotfiles+'/vim/bundle/'
if not os.path.exists(bundle):
    os.makedirs(bundle)

def checkit(package):
    out = git_getfoldername(package)
    checkfolder  = bundle+"/"+out
    if not os.path.isdir(checkfolder):
        print("dne ",checkfolder)
        os.chdir(bundle)
        gitclone(package)
    else:
        print("OK ",checkfolder)
    return

vundle =  "https://github.com/gmarik/Vundle.vim.git"
vimsignature = "https://github.com/kshenoy/vim-signature.git"
#solarized = "git://github.com/altercation/vim-colors-solarized.git" # works also without
#restoreview = "https://github.com/vim-scripts/restore_view.vim.git" # works also without
#pyfoldering = "git://github.com/vim-scripts/Efficient-python-folding.git" # works also without
systastic = "https://github.com/scrooloose/syntastic.git"
#nerdtree = "git://github.com/scrooloose/nerdtree.git"
vimlatex = "https://github.com/jcf/vim-latex.git"


#pymode = "git://github.com/klen/python-mode.git"
#taglist = "git://github.com/vim-scripts/taglist.vim.git"
#commenter = "git://github.com/scrooloose/nerdcommenter.git"
#git submodule add https://github.com/jcf/vim-latex.git dotfiles/vim/bundle/vim-latex
tabcompletion = "https://github.com/ervandew/supertab.git"
toggletruefalse = "https://github.com/AndrewRadev/switch.vim.git"
checkit(vundle)         # necessary
checkit(vimsignature)   # necessary  set a mark (prss m c in normal mode)
checkit(toggletruefalse) # necessary switsches True to False and so non
checkit(tabcompletion)   # necesary comletes words when <TAB>

#checkit()
