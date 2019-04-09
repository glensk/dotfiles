
import numpy as np
from scipy.optimize import minimize_scalar
import sympy
import sys
import warnings
import feos

V = sympy.Symbol('V')

def surf_eos(V):
    return feos.vinet(V, 1.2,1.2,1.2,1.2) #, convert_to_sympy=sympy_)

def f2(V):
    return (V - 2)*V*(V + 2)**2

def f2der(V):
    return sympy.diff(f2(V),V)

def f3(V):
    return 3*V**3 - 8* V - 7* V**2
def f4(V):
    return 3*V**3 - 8* V - 7* V**2+2*V**4
def f4der(V):
    return sympy.diff(f4(V),V)

f4der_ = 3*V**3 - 8* V - 7* V**2+2*V**4

print feos.vinet(V, 1.2,1.2,1.2,1.2, convert_to_sympy=True)
print minimize_scalar(surf_eos).x
print minimize_scalar(f2).x
print "<<"
####################################################################
# minimize_scalar
#########################################################################
print minimize_scalar(f3, bounds=(-4, -1), method='bounded').x
print "vvv bei -1.79 und 1.15 vvvvvvvvv"
print minimize_scalar(f4, bounds=(-3, -1), method='bounded').x
print minimize_scalar(f4, bounds=(0, 1000), method='bounded').x

#########################################################################
# optimize.fmin_bfgs
#########################################################################
from scipy import optimize
print optimize.fmin_bfgs(f4, 7)
#print optimize.fmin_cg(f4, [2, 2])
print "vvv DER: ca 0.45 vvvvvvvvv"
#print sympy.diff(f4(V), V)
#print minimize_scalar(sympy.diff(f4(V),V))


#        def surf_eos_der(V):
#            #return feos.vinet(V, *self.feos.parameters) #, convert_to_sympy=sympy_)
