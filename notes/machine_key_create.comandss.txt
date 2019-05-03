****************Create a machine key : **********************
> ssh-keygen -t dsa
(always type 'Enter')

Copy the public key tothe remote machine
> scp ~/.ssh/id_dsa.pub glensk@<remote>:~/ssh_key.pub

SSH Remote Machine then :
> cat ~/ssh_key.pub >> ~/.ssh/authorized_keys

or

ssh-keygen -t dsa  (immer Enter Pw AUSSUCHEN)
touch .ssh/authorized_keys2
cat .ssh/id_dsa.pub >> .ssh/authorized_keys2
scp -r .ssh .ssh2 cmpc72:                 where cmpc72 is any remotepc
scp ~/.ssh/id_dsa.pub cmpc72:
