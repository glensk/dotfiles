###################
## STRATO
###################
sshfs glensk@sftp.hidrive.strato.com:/users/glensk /data/glensk/strato -o reconnect -C -o workaround=all,Ciphers="blowfish-cbc" (/data/glensk/strato needs to exist)
# sshfs -o IdentityFile=/path/to/some/ssh_private_key user@scp.hidrive.strato.com:/ /some/local/mount_point 

###################
## NOTBOOK at MPIE
###################
sshfs glensk@cmpc08:/nas/glensk/ /nas/glensk -o reconnect -C -o workaround=all,Ciphers="blowfish-cbc"   (/nas/glensk needs to exist)
sshfs glensk@cmpc08:/data/glensk/ /data/glensk -o reconnect -C -o workaround=all,Ciphers="blowfish-cbc" (/data/glensk needs to exist)
sshfs glensk@cmpc08:/home/grabowski/ /home/grabowski -o reconnect -C -o workaround=all,Ciphers="blowfish-cbc" (/home/grabowski needs to exist)

############
## GRZ
############
sshfs aglen@rzgate:/afs/ipp-garching.mpg.de/home/a/aglen /home/glensk/rzg -o reconnect -C -o workaround=all,Ciphers="blowfish-cbc"
sshfs glensk@cmpc08:/nas/glensk/v /home/glensk/v -o reconnect -C -o workaround=all,Ciphers="blowfish-cbc"


###############
## ??
###############
ssh als SOCKS-Proxy: ssh -D 1234 user@box_B          ## Surfen ueber ein gate wenn Zensur
        --> siehe        http://www.daniel-ritter.de/blog/6-nutzliche-dinge-die-man-mit-ssh-tun-kann


change rights of folder within /
sudo chown X /var/www -R
Where X is the username.

##################################################
#found:
sshfs#username@myserver:/storage/Data    /home/myuser/Documents/Data        fuse comment=sshfs,rw,exec,uid=1000,allow_other,reconnect,transform_symlinks,BatchMode=yes,nonempty,noauto 0 0

################################################## 
## fidis
################################################## 
sshfs glensk@fidis.epfl.ch:/scratch/glensk/ /Users/glensk/fidis -o reconnect -C


######################################
### unmount
######################################
linux: fusermount -u /data/glensk/strato
mac: umount -f /nas/glensk

