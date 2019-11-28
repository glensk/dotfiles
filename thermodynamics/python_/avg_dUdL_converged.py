import os
import sys
import numpy             as np


if len(sys.argv) is not 2:
    print "ERROR: you need to give the convergence criterion"
    sys.exit(0)
else:
    #print "--",sys.argv[1]
    #print "--",type(sys.argv[1])
    #print "00000"
    #try:
    conv = float(sys.argv[1])
    #print "00000"
    #print "conv:",conv,type(conv)
d = np.loadtxt("avg_dUdL_fre")
l = np.loadtxt("avg_dUdL_fre",dtype=str)
di = d[:,[0,2]]
li = l[:,[0,2]]
if di.size == 2:
    if di[1] <= conv:
	print li[0]
else:
#print "0000:"
#print len(di)
#print "x",di.size
#print di[:,1]
#print li
#print di[:,1]
#print "0000:"
#dl = di[di[:,1] <= conv][:,0]
    ll = li[di[:,1] <= conv][:,0]
    #print p
    #print np.array_str(ll)
    #print np.array_repr(ll)
    #print np.ndarray.tolist(ll)
    #print "---"
    #print ll
    #print [''.join(row) for row in ll]
    for row in ll:
        print str(row)
    #print str(dl)
