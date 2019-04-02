spaetestens bei make braucht man su rechte!   su manke install

1) yum install nautilus-devel
2) yum install python-docutils
3) tar xjf ./nautilus-dropbox-1.4.0.tar.bz2
4) cd ./nautilus-dropbox-1.4.0.tar.bz2; ./configure; make; (su) make install;


[824] 11:26 glensk@cmpc08  ~/progs/nautilus-dropbox-1.4.0                                                                                                                                                            $make install
Making install in data
make[1]: Entering directory `/data/glensk/progs/nautilus-dropbox-1.4.0/data'
Making install in icons
make[2]: Entering directory `/data/glensk/progs/nautilus-dropbox-1.4.0/data/icons'
Making install in hicolor
make[3]: Entering directory `/data/glensk/progs/nautilus-dropbox-1.4.0/data/icons/hicolor'
Making install in 16x16
make[4]: Entering directory `/data/glensk/progs/nautilus-dropbox-1.4.0/data/icons/hicolor/16x16'
Making install in apps
make[5]: Entering directory `/data/glensk/progs/nautilus-dropbox-1.4.0/data/icons/hicolor/16x16/apps'
make[6]: Entering directory `/data/glensk/progs/nautilus-dropbox-1.4.0/data/icons/hicolor/16x16/apps'
make[6]: Nothing to be done for `install-exec-am'.
test -z "/usr/share/icons/hicolor/16x16/apps" || /bin/mkdir -p "/usr/share/icons/hicolor/16x16/apps"
 /usr/bin/install -c -m 644 dropbox.png '/usr/share/icons/hicolor/16x16/apps'
/usr/bin/install: cannot create regular file `/usr/share/icons/hicolor/16x16/apps/dropbox.png': Permission denied
make[6]: *** [install-iconsDATA] Error 1
make[6]: Leaving directory `/data/glensk/progs/nautilus-dropbox-1.4.0/data/icons/hicolor/16x16/apps'
make[5]: *** [install-am] Error 2
make[5]: Leaving directory `/data/glensk/progs/nautilus-dropbox-1.4.0/data/icons/hicolor/16x16/apps'
make[4]: *** [install-recursive] Error 1
make[4]: Leaving directory `/data/glensk/progs/nautilus-dropbox-1.4.0/data/icons/hicolor/16x16'
make[3]: *** [install-recursive] Error 1
make[3]: Leaving directory `/data/glensk/progs/nautilus-dropbox-1.4.0/data/icons/hicolor'
make[2]: *** [install-recursive] Error 1
make[2]: Leaving directory `/data/glensk/progs/nautilus-dropbox-1.4.0/data/icons'
make[1]: *** [install-recursive] Error 1
make[1]: Leaving directory `/data/glensk/progs/nautilus-dropbox-1.4.0/data'
make: *** [install-recursive] Error 1

as su:
apt-get install libnautilus-extension-dev

http://ogfomk.com/wordpress/2011/06/29/dropbox-for-nautilus-debian/
