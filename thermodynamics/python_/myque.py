#!/usr/bin/python

def qls():
    import socket
    import os
    import my

    host = socket.gethostname()
    print "host:",host

    if host == "cmmd010":
        print "qls is at cmmd010 merely an output and cant be used further with python"
        qls = os.system("qls.origin -u")
    else:
        qls = my.run("qls")
    return qls

#def qstat():
#    import socket
#    import os
#    import my
#
#    host = socket.gethostname()
#
#    if host == "cmmd010":
#        print "qstat is at cmmd010 merely an output and cant be used further with python"
#        qstat = os.system("/opt/sge/bin/lx26-amd64//qstat")
#
#    else:
#        qstat = my.run("ssh cmmd010 qstat")
#    return qstat
