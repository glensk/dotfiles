#!/usr/bin/env python


import numpy as np
from scipy.optimize import minimize_scalar
import sympy
from sympy import init_printing
init_printing(use_unicode=False, wrap_line=False, no_global=True)
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
    a= 3*V**3 - 8* V - 7* V**2+2*V**4
    print "aaaaaaaaaaaaaaaa|",a,"||",type(a)
    return a

def f4der(V):
    return 8*V**3 + 9*V**2 - 14*V - 8

f4_ = 3*V**3 - 8* V - 7* V**2+2*V**4
f4der_ = 8*V**3 + 9*V**2 - 14*V - 8

print feos.vinet(V, 1.2,1.2,1.2,1.2, convert_to_sympy=True)
print minimize_scalar(surf_eos).x
print minimize_scalar(f2).x
print "<<"
print minimize_scalar(f3, bounds=(-4, -1), method='bounded').x
print "vvv bei -1.79 und 1.15 vvvvvvvvv"
print minimize_scalar(f4, bounds=(-3, -1), method='bounded').x
print minimize_scalar(f4, bounds=(0, 1000), method='bounded').x

print "DER bei 0.47 und xxx vvvvvvvvv"
print minimize_scalar(f4der, bounds=(-3, 10), method='bounded').x

print "bfgs:vvv"
from scipy import optimize
print optimize.fmin_bfgs(f4, 7)
#print optimize.fmin_cg(f4, [2, 2])
print "vvv DER: ca 0.45 vvvvvvvvv"
d4 = sympy.diff(f4(V), V)
d44 = f4(V).diff(V)
def e(V):
    #return eval(f4(V).diff(V))
    #return f4(V).diff(V)
    #return sympy.diff(f4(V), V)
    return eval(sympy.diff(f4(V), V))

print "d44:",d44,type(d44)
print "d4 :",d4,type(d44)
print 'eeeeeeeeeeeeeeeeeeeeeeeeee',e,type(e)
#print minimize_scalar(e)
#print "f4:", sympy.solve(f4(V), V)
#print "f4der:", sympy.solve(f4der(V), V)
#from scipy.optimize import fsolve
#print "f4:", fsolve(f4,0)
###print minimize_scalar(f4(V).diff(V))
##print d4,type(d4),type(f4)
###print minimize_scalar(d4)
##print sympy.solve(f4der_)
##print sympy.solve(f4_)
#print "-----------"*10
#from sympy.utilities.lambdify import lambdify, implemented_function
#f = lambdify([V], 3*V**3 - 8* V - 7* V**2+2*V**4, 'numpy')
#
#x_vec = np.arange(0, 10, 0.1)
#y_vec = f(x_vec)
#print y_vec

#        def surf_eos_der(V):
#            #return feos.vinet(V, *self.feos.parameters) #, convert_to_sympy=sympy_)
#            return sympy.diff(feos.vinet(V,*self.feos.parameters,convert_to_sympy=True),V)
