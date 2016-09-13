#!/bin/sh

echo Add the following 4 lines to your macs $HOME/.ssh/config:
echo Host gate
echo     User aglen
echo     Hostname gate.rzg.mpg.de 
echo     IdentityFile /Users/glensk/.ssh/gate


cd $HOME/.ssh
ssh-keygen -t rsa -f gate
scp gate.pub gate:.ssh/

echo "now log in to gate (again: use normal password, skip key passphrase)"
echo "ssh gate"
echo "cd .ssh/"
echo
echo "(if ./authorized_keys does not exist, use"
echo "touch ./authorized_keys"
echo "chmod 644 ./authorized_keys"
echo ")"
echo
echo "cat ./garching.pub >> ./authorized_keys"
echo
echo "logout from gate:"
echo "exit"
echo
echo "Now, you may test the private key"
echo "ssh gate"
