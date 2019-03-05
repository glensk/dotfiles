ssh cmmd010 "module list;module load vasp/parallel/4.6-tdi;module list;ls;cd scripts;ls"    ## es gehen auch mehrere befehle auf einmal

#############################################################
## change permissions to login to nodes directly
## you want to add cmmd010 to cmpc08
#############################################################

quick:  @ cmpc08:    cd /home/$USER  
        @ cmpc08:    ssh-keygen -t rsa  (and just type return; in case you did this before or /home/$USER/.ssh/id_rsa.pub existst, skipt this)
        @ cmpc08:    cat .ssh/id_rsa.pub | ssh `echo $USER`@cmmd002 'cat >> .ssh/authorized_keys'         ### type password last time
        chmod 700  .ssh  (on either machine)
        chmod 640  .ssh/authorized_keys2  .ssh/authorized_keys  ## (on either machine)


long: 
How to do it

First log in on A as user a and generate a pair of authentication keys. Do not enter a passphrase:

a@A:~> ssh-keygen -t rsa
Generating public/private rsa key pair.
Enter file in which to save the key (/home/a/.ssh/id_rsa): 
Created directory '/home/a/.ssh'.
Enter passphrase (empty for no passphrase): 
Enter same passphrase again: 
Your identification has been saved in /home/a/.ssh/id_rsa.
Your public key has been saved in /home/a/.ssh/id_rsa.pub.
The key fingerprint is:
3e:4f:05:79:3a:9f:96:7c:3b:ad:e9:58:37:bc:37:e4 a@A

Now use ssh to create a directory ~/.ssh as user b on B. (The directory may already exist, which is fine):

a@A:~> ssh b@B mkdir -p .ssh
b@B's password: 

Finally append a's new public key to b@B:.ssh/authorized_keys and enter b's password one last time:

a@A:~> cat .ssh/id_rsa.pub | ssh `echo $user`@cmmd010 'cat >> .ssh/authorized_keys'
b@B's password: 

From now on you can log into B as b from A as a without password:

a@A:~> ssh b@B hostname
B


#############################################################################################################################################

quick:  @ cmpc08:    cd /home/$USER  
        @ cmpc08:    ssh-keygen -t rsa  (and just type return; in case you did this before or /home/$USER/.ssh/id_rsa.pub existst, skipt this)
        @ cmpc08:    cat .ssh/id_rsa.pub | ssh `echo $USER`@cmmd002 'cat >> .ssh/authorized_keys'         ### type password last time
        chmod 700 .ssh;chmod 640 .ssh/authorized_keys2  .ssh/authorized_keys  ## (on either machine)
