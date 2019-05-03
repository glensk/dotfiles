      vasprun.xml                                   92 MB
gzip  vasprun.xml                                   26 MB
tarfolder.sh vasprun.xml                            25 MB
7z a vasprun.xml.7z vasprun.xml                     16 MB
zpaq a archive.zpaq vasprun.xml -m5                 14 MB (16:20:45 - 16:25:30)
zpaq a archive.zpaq vasprun.xml -m5 -t2             14 MB (16:26:38 - 16:31:15
tar -c --lzma -f vasprun.xml.tar.lzma vasprun.xml   18 MB (16:49:13 - 16:5

zpaq x archive.zpaq              # is enough to unpack stuff

