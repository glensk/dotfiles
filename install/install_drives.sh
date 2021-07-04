#!/bin/sh


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
read -p "Did everything go smoothly? If not ctrc+c; otherwise just press enter" 
ln -s /nas/glensk/v $HOME/v


