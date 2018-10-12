Entpacken:

tar xfv filename.tar 
tar xfvz filename.tar.gz 
tar -xvf file.tar.gz
tar xfvi filename.tar.bz2
gzip -d file.gz
tar -xvf file.tgz
uncompress POTCAR.Z
zcat POTCAR.Z > POTCAR
unzip file.zip (auf cmpc08, not cmmd010)

Packen:
tar -zcvf archive-name.tar.gz directory-name  ## whole directory
gzip file --> file.gz
