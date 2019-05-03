sudo port selfupdate
sudo port upgrade outdated
sudo port -v install vim +python +x11
vim --version
VIM - Vi IMproved 7.3 (2010 Aug 15, compiled May  2 2013 15:16:05), Included patches: 1-244, 246-762
port select --list python
sudo port select --set python python27


to uninstall port:
sudo port -fp uninstall installed
sudo rm -rf \
    /opt/local \
    /Applications/DarwinPorts \
    /Applications/MacPorts \
    /Library/LaunchDaemons/org.macports.* \
    /Library/Receipts/DarwinPorts*.pkg \
    /Library/Receipts/MacPorts*.pkg \
    /Library/StartupItems/DarwinPortsStartup \
    /Library/Tcl/darwinports1.0 \
    /Library/Tcl/macports1.0 \
    ~/.macports
