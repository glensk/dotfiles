#!/bin/sh

# to install homebrew
ruby -e "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/master/install)"

echoit() {
echo ""
echo ""
echo ""
echo "#####################################################################"
echo "#####################################################################"
echo "#####################################################################"
echo "$*"
echo "#####################################################################"
echo "#####################################################################"
echo "#####################################################################"
echo ""
echo ""
echo ""
}




echoit "creating /nas /data /u"
owner=$USER
user=glensk;sys=nas;
sudo mkdir -p /$sys/$user;sudo chown $owner /$sys;sudo chown $owner /$sys/$user; echo /$sys/$user
user=glensk;sys=data;
sudo mkdir -p /$sys/$user;sudo chown $owner /$sys;sudo chown $owner /$sys/$user; echo /$sys/$user
user=grabowski;sys=nas;
sudo mkdir -p /$sys/$user;sudo chown $owner /$sys;sudo chown $owner /$sys/$user; echo /$sys/$user
user=grabowski;sys=data;
sudo mkdir -p /$sys/$user;sudo chown $owner /$sys;sudo chown $owner /$sys/$user; echo /$sys/$user
user=korbmacher;sys=nas;
sudo mkdir -p /$sys/$user;sudo chown $owner /$sys;sudo chown $owner /$sys/$user; echo /$sys/$user
user=aglen;sys=u;
sudo mkdir -p /$sys/$user;sudo chown $owner /$sys;sudo chown $owner /$sys/$user; echo /$sys/$user
user=aglen;sys=cmmc/u;
sudo mkdir -p /$sys/$user;sudo chown $owner /$sys;sudo chown $owner /$sys/$user; echo /$sys/$user
sudo mkdir -p /cmmc/ptmp/aglen;sudo chown $owner /cmmc/ptmp; sudo chown $owner /cmmc/ptmp/aglen
read -p "Did everything go smoothly? If not ctrc+c; otherwise just press enter" 
ln -s /nas/glensk/v $HOME/v





echoit "brew intall all kind of stuff"
# notes from http://www.moncefbelyamani.com/how-to-install-xcode-homebrew-git-rvm-ruby-on-mac/

brew doctor  # necesssary after first install
read -p "#### Did everything go smoothly with brew doctor? If not ctrc+c; otherwise just press enter" 
brew update 
read -p "#### Did everything go smoothly with brew updat? If not ctrc+c; otherwise just press enter" 


brew tap 'caskroom/cask'
brew install 'brew-cask'
# check at http://caskroom.io/search if your package is available
brew cask install --appdir="/Applications" iterm2
brew cask install --appdir="/Applications" google-chrome
brew cask install --appdir="/Applications" dropbox
brew cask install --appdir="/Applications" alfred
brew cask install --appdir="/Applications" vlc 
brew cask install --appdir="/Applications" pdftk
#brew cask install --appdir="/Applications" evernote
brew cask install --appdir="/Applications" mendeley-desktop
brew cask install --appdir="/Applications" gimp
brew cask install --appdir="/Applications" anaconda
brew cask install --appdir="/Applications" aquaterm
brew cask install --appdir="/Applications" macvim
brew cask install --appdir="/Applications" djview
brew cask install --appdir="/Applications" skim
brew cask install --appdir="/Applications" skype
brew cask install --appdir="/Applications" utorrent
brew cask install --appdir="/Applications" omnifocus
brew cask install --appdir="/Applications" omnifocus-clip-o-tron
brew cask install --appdir="/Applications" paintbrush
brew cask install --appdir="/Applications" xquartz
brew cask install --appdir="/Applications" seil
brew cask install --appdir="/Applications" cuda   # necessary for vmd (make movies from VASP MD)
brew cask install --appdir="/Applications" flux 
brew cask install --appdir="/Applications" bartender
brew cask install --appdir="/Applications" bettertouchtool
brew cask install --appdir="/Applications" textexpander
brew cask install --appdir="/Applications" spyder
brew cask install --appdir="/Applications" calibre  # to get ebook-convert file.pdf file.epub --enable-heuristics



# install lammps (from http://lammps.sandia.gov/download.html#git)
brew tap homebrew/science
#brew install lammps              # serial version
brew install lammps --with-mpi   # mpi support 
brew install rmtrash # 
brew install fondu  # for matplotlib fonts
brew install lbzip2 # to parallel zip files
brew install boost  # c libraries necessary fot fftw (for phonon lifetimes project), also for saschas c++ skrpit 


brew install octave        # for phonon-lifetime code of michael leitner aus muenchen
brew install ssh-copy-id   # to copy public key to remote host
brew install macvim        # this will be used as macvim -v == gvim -v for clipboard support between cmmc/cmpc gvim sessions
brew install git
git config --global credential.helper osxkeychain

brew install coreutils findutils gnu-tar gnu-sed gawk gnutls gnu-indent gnu-getopt --with-default-names   # --default-names ist wichtig!
# brew uninstall coreutils findutils gnu-tar gnu-sed gawk gnutls gnu-indent gnu-getopt  # just in case wrongly installed
brew install fftw --with-fortran
brew install cpulimit   # to limit dropbox usage
brew install unrar      # to unrar *.rar files
brew install colordiff
brew install gcc        # which now also installs gfortran gfortran-4.9
brew install mmv        # --> requires XQUARTZ
brew install grace      # xmgrace  (needs openmotif)
brew install youtube-dl
brew install imagemagick  # to change pdf to jpg
#brew install gnuplot --with-x 
brew install gnuplot --with-aquaterm
brew install p7zip      # necessary to unzip michael leitners zipped files
brew install popper     # pdftops
brew install tesseract  # get orc for pdfs (make pdfs readable)
brew install gs         # also necessay for making pdfs readable  
brew install qpdf       # to compresss pdfs
brew install shellcheck
brew install ffmpeg     # to convert mp4 in mp3
#brew install qscintilla2 # for joergs work bench: from PyQt4 import QtGui,QtCore, QtSql;PyQt4.Qsci import QsciScintilla 
brew install bash       # to get latest version of bash; osx bash is too old (no colors of files on tab)
echo /usr/local/bin/bash|sudo tee -a /etc/shells;chsh -s /usr/local/bin/bash

chsh -s /usr/local/bin/bash $USER   # change schell to brew shell
                                    # additional note in commands/bash
# it can be that you also need to add this ot /etc/shells on fist place
brew install wget --with-iri
brew install svn
#brew install vim   # to get +clipboard on vim --version
brew install vim --with-client-server --with-lua --with-features=huge --with-xterm_clipboard
brew install sshfs  # to make possible to mount nas / data  # also installes osxfuse!
                    # which is necessary to mout nas / data
                    # is this is not working (yosemite) then go to 
                    # http://osxfuse.github.io and install fuse and sshfs package
brew install ack    # like grep but cooler partly
brew install bash-completion
brew install ctags  # for tags in vim

#brew install macvim --with-cscope --with-lua --HEAD
#brew install macvim --custom-icons --override-system-vim --with-lua --with-luajit
brew install macvim --with-cscope --custom-icons --override-system-vim --with-lua --with-luajit

#    CMake (also for macvim and youcompleteme)
#            in /Users/glensk/ycm_build:
#            /usr/local/Cellar/cmake/2.8.12.2/bin/cmake -G "Unix Makefiles" . ~/.vim/bundle/YouCompleteMe/third_party/ycmd/cpp

#brew install https://github.com/downloads/zolrath/wemux/wemux.rb
           

#brew install mp3info  # to change mp3 tags massively
#mp3info  -> does not??? support id3v2
#id3tool  -> does not support id3v2
#id3lib  ----> installs id3tag
brew install homebrew/science/armadillo  # for phonon linewidths



echo "Force the Dock to only show running applications System"
defaults write com.apple.dock static-only -bool TRUE

echo "Disable the switching Spaces animation in the terminal by typing:"
defaults write com.apple.dock workspaces-swoosh-animation-off -bool YES && killall Dock

echo "to display full path in finder"
defaults write com.apple.finder _FXShowPosixPathInTitle -bool YES


echo ""
echo ""
echo run: brew info osxfuse and do the stuff which needs to be done as sudo
echo ""
echo ""
echo then:
echo ""
echo "#################### reboot system once :) ######################"
echo "#################### reboot system once :) ######################"
echo "#################### reboot system once :) ######################"
echo "#################### reboot system once :) ######################"
echo "#################### reboot system once :) ######################"

# latex (necessary for general latex commands)
hier=`pwd`
echo " get MacTex package which should be installed then manually if the MacTex.pkg was not exected (should be .... since open)"
echo " go to $HOME/Downloads and install the MacTeX.pkg (just klick on it)"
cd $HOME/Downloads
echo
echo BETTER DOWNLOAD WITH SAFARI SINCE MUCH QUICKER THAN CHROME
echo
wget -r http://mirror.ctan.org/systems/mac/mactex/MacTeX.pkg
cd $hier
open $HOME/Downloads/mirror.ctan.org/systems/mac/mactex/MacTex.pkg

echo " also try texstudio, i did like it better working at the cmpc"


#in texShop preferences: engine
#pdflatex: /usr/texbin changed to /Library/TeX/texbin
