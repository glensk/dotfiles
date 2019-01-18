#!/usr/bin/python
import argparse
import re
import subprocess,os
import numpy as np
import time
import scipy.linalg as salg
import scipy.sparse.linalg as spalg
import sys
import scipy.integrate as spint

###########################################################
## I WILL TEMPORARILY DEFINE HERE THE SYMMETRY FUNCTIONS ##
## IN THE FUTURE THIS SHOULD BE DONE BETTER              ##
###########################################################

def nvfc(x):
    if x>1: return 0
    return np.tanh(1.0-x)**3
fc = np.vectorize(nvfc)

def g2(r, rc, eta, rs=0.0):
    return np.exp(-eta*(r-rs)**2)*fc(r/rc)

def g3(rij, rik, theta, rc, eta, zeta, l):
    costh = np.cos(theta)
    rjk = np.sqrt(rij**2+rik**2-2*rij*rik*costh)
    return np.exp(-eta*(rij*rij+rik*rik+rjk*rjk))*fc(rij/rc)*fc(rik/rc)*fc(rjk/rc)*2**(1.0-zeta)*(1+l*costh)**zeta

def g9(rij, rik, theta, rc, eta, zeta, l):
    costh = np.cos(theta)
    return np.exp(-eta*(rij*rij+rik*rik) ) * fc(rij/rc)*fc(rik/rc) * 2**(1.0-zeta)*(1.0+l*costh)**zeta

###########################################################

def SF_integrate(sym, rho={}):
    """Calculate the integral of the SF for an homogeneous gas"""
    nsym = len(sym)
    SF_int = np.zeros(nsym)
    for idx in xrange(nsym):
        if sym[idx][1] == "2":  # Two body
            rc = float(sym[idx][-1])
            rs = float(sym[idx][-2])
            eta = float(sym[idx][-3])
            rho1 = rho[sym[idx][2]]
            #print('aa 2')
            SF_int[idx] = rho1*np.pi*4*spint.quad(lambda x: x**2*g2(x,rc,eta,rs),0,rc)[0]
        elif sym[idx][1] == "3":  # Three body
            rc = float(sym[idx][-1])
            zeta = float(sym[idx][-2])
            l = float(sym[idx][-3])
            eta = float(sym[idx][-4])
            rho1 = rho[sym[idx][2]]
            rho2 = rho[sym[idx][3]]
            #print('aa 3')
            SF_int[idx] = rho1*rho2*np.pi**2*8*spint.tplquad(lambda x, y, th: (x**2*y**2*np.sin(th)*g3(x, y, th,rc, eta, zeta, l)), 0,np.pi,lambda a:0,lambda a:rc,lambda a,b:0,lambda a,b:rc )[0]
        elif sym[idx][1] == "9":
            rc = float(sym[idx][-1])
            zeta = float(sym[idx][-2])
            l = float(sym[idx][-3])
            eta = float(sym[idx][-4])
            g9itheta = lambda rij, rik: (np.exp(-eta*(rij*rij+rik*rik))*fc(rij/rc)*fc(rik/rc)*2**(1.0-zeta)*((1.0+l)**(1+zeta)-(1.0-l)**(1+zeta))/(l*(1.0+zeta)) )
            rho1 = rho[sym[idx][2]]
            rho2 = rho[sym[idx][3]]
            SF_int[idx] = rho1*rho2*np.pi**2*8*spint.dblquad(lambda x, y: (x**2*y**2*g9itheta(x,y)), 0,rc,lambda a:0,lambda a:rc )[0]
        else:
            raise ImplementationError('Not yet implemented with the symmetry function choice you selected')
        print idx,"out of",nsym,"--",sym[idx], SF_int[idx]

    return SF_int

def GetCosts(sym, rho):
    # estimates the number of pairs/triplets to be computed for each SF
    costs = np.zeros(len(sym))
    i = 0
    ftpi = np.pi * 4.0/3.0
    for s in sym:
        radius = float(s[-1])
        # radial syms have a cost scaling with the cube of radius,
        # angular use triplets so it's the sixth power
        if s[1] == '2':
            costs[i] = ftpi*radius**3 * rho[s[0]] * rho[s[2]]
        else:
            costs[i] = ftpi**2*radius**6 * rho[s[0]] * rho[s[2]] * rho[s[3]]
        i += 1

    return costs
