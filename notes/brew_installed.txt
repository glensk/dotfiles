brew uninstall coreutils findutils gnu-tar gnu-sed gawk gnutls gnu-indent gnu-getopt
brew install coreutils findutils gnu-tar gnu-sed gawk gnutls gnu-indent gnu-getopt --default-names   # --default-names ist wichtig!

brew install
    wget
    ctags
    svn
    vim   # to get +clipboard on vim --version
    sshfs # to make possible to mount nas / data
    brew install mp3info  # to change mp3 tags massively

    mp3info  -> does not??? support id3v2
    id3tool  -> does not support id3v2
    id3lib  ----> installs id3tag
    ack
    bash-completion

    macvim with:
    brew install macvim --with-cscope --with-lua --HEAD

    CMake (also for macvim and youcompleteme)
            in /Users/glensk/ycm_build:
            /usr/local/Cellar/cmake/2.8.12.2/bin/cmake -G "Unix Makefiles" . ~/.vim/bundle/YouCompleteMe/third_party/ycmd/cpp

brew install https://github.com/downloads/zolrath/wemux/wemux.rb
            
