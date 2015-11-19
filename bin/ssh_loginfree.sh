#!/bin/sh

echo on mac do
ssh-keygen -t rsa
cat ~/.ssh/id_rsa.pub | ssh glensk@cmpc34 'cat >> .ssh/authorized_keys'
cat ~/.ssh/id_rsa.pub | ssh glensk@cmpc34 'cat >> .ssh/authorized_keys2'

###### NOT DONE YET:
cat ~/.ssh/id_rsa.pub | ssh glensk@cmpc34.mpie.de 'cat >> .ssh/authorized_keys'
cat ~/.ssh/id_rsa.pub | ssh glensk@cmpc34.mpie.de 'cat >> .ssh/authorized_keys2'
#Put the public key in .ssh/authorized_keys2
#Change the permissions of .ssh to 700
#Change the permissions of .ssh/authorized_keys2 to 640

