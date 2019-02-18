sshfs glensk@sftp.hidrive.strato.com:/users/glensk /data/glensk/strato -o reconnect -C -o workaround=all,Ciphers="blowfish-cbc"
sshfs glensk@sftp.hidrive.strato.com:/users/glensk /data/glensk/strato -o reconnect,Ciphers="blowfish-cbc"
fusermount -u /data/glensk/strato

sshfs glensk@cmpc08:/nas/glensk/ /home/glensk/n -o reconnect -C -o workaround=all,Ciphers="blowfish-cbc"




# sshfs -o IdentityFile=/path/to/some/ssh_private_key user@scp.hidrive.strato.com:/ /some/local/mount_point 

ssh als SOCKS-Proxy: ssh -D 1234 user@box_B          ## Surfen ueber ein gate wenn Zensur
        --> siehe        http://www.daniel-ritter.de/blog/6-nutzliche-dinge-die-man-mit-ssh-tun-kann

daten abgleichen mit unison
