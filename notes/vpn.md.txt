reaction:~ user$ ps -ef | grep racoon
    0   265     1   0   0:00.22 ??         0:00.34 /usr/sbin/racoon
      501   339   335   0   0:00.00 ttys001    0:00.00 grep racoon
      reaction:~ user$ sudo kill -9 265
      Password:
reaction:~ user$ ps -ef | grep racoon
  501   345   335   0   0:00.00 ttys001    0:00.00 grep racoon
  reaction:~ user$ sudo /usr/sbin/racoon
  reaction:~ user$ ps -ef | grep racoon
      0   347     1   0   0:00.00 ??         0:00.01 /usr/sbin/racoon -x
        501   349   335   0   0:00.00 ttys001    0:00.00 grep racoon



oder:
sudo launchctl stop com.apple.racoon
sudo launchctl start com.apple.racoon
