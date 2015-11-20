#!/bin/sh

#Host cmmc*
#    HostName %h.bc.rzg.mpg.de
#    User aglen
#    StrictHostKeyChecking no
#    UserKnownHostsFile /dev/null
#
#\ssh aglen@cmmc002.bc.rzg.mpg.de  -> ssh cmmc002

echo on mac do
ssh-keygen -t rsa
cat ~/.ssh/id_rsa.pub | ssh glensk@cmpc34 'cat >> .ssh/authorized_keys'
cat ~/.ssh/id_rsa.pub | ssh glensk@cmpc34 'cat >> .ssh/authorized_keys2'

###### NOT DONE YET:
cat ~/.ssh/id_rsa.pub | ssh glensk@cmpc34.mpie.de 'cat >> .ssh/authorized_keys'
cat ~/.ssh/id_rsa.pub | ssh glensk@cmpc34.mpie.de 'cat >> .ssh/authorized_keys2'
chmod 700 ~/.ssh;chmod 640  ~/.ssh/authorized_keys2  ~/.ssh/authorized_keys
#Put the public key in .ssh/authorized_keys2
#Change the permissions of .ssh to 700
#Change the permissions of .ssh/authorized_keys2 to 640
ssh-copy-id -i ~/.ssh/id_rsa.pub glensk@cmpc34


# http://www.softpanorama.org/Net/Application_layer/SSH/passwordless_ssh_login.shtml
# cd:cd: key has expired: /home/glensk/.ssh
# cd: key has expired: settitlepath.sh
# try manual
# \ssh glensk@cmpc34
#
# ls: cannot open directory .: Key has expired
# try now with cmpc where login is set to 16000 and not 1600 anymore
