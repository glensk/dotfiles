#!/bin/sh

echo "## this works :) (for mac) however, the /nas/glensk/ptmp mount on cmpc disappears after a certain time in contrast to /nas/glensk/cmmc"
echo "## on cmpc ptmp was mounted by: .... christophs command, try now with same command below on cmpc"
sshfs -o idmap=user,uid=10010,gid=10010 aglen@cmmc001.bc.rzg.mpg.de:/cmmc/ptmp/aglen /nas/glensk/ptmp/ -o reconnect -C -o workaround=all,Ciphers="blowfish-cbc",transform_symlinks,BatchMode=yes
#echo "THIS DOES NOT WORK ON mac", sshfs aglen@cmmc001.bc.rzg.mpg.de:/cmmc/ptmp/aglen /nas/glensk/ptmp -o reconnect -C -o workaround=all,transform_symlinks

echo to umount: fusermount -u /nas/glensk/ptmp

