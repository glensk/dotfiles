#!/usr/bin/env python
# -*- coding: utf-8 -*-
from __future__ import print_function
from socket import gethostname
import os
#hostname = gethostname()   # fidis, helvetios, daint105, h332, g037, f221


def check_for_known_hosts(exit=False):
    hostname = gethostname()   # fidis, helvetios, daint105, h332, g037, f221
    known_hosts =  ['fidis','helvetios', "daint", "mac", "cmmd", 'cosmopc' ]
    for i in known_hosts:  # also works with daint102, helvetios101, ...
        if i in hostname: return i

    if hostname[0] == 'h' and is_int(hostname[1:]) == True:
        return 'helvetios'
    if hostname[0] in ['g','f'] and is_int(hostname[1:]) == True:
        return 'fidis'

    # if not known host
    if exit == True:
        print("known hosts:",known_hosts)
        sys.exit(hostname+" is not in the list of known hosts!")
    else:
        return False

if __name__ == '__main__':
    hostknown = check_for_known_hosts(exit=False)
    print(hostknown)
