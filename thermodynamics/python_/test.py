#!/usr/bin/python
import get, my, argparse, sys, os, db
print "ahllo"
my.exit("joexit")
#def f1(V,E0,V0,BM,BMder):
#    return E0+(4*BM*V0)/(BMder-1)**2-2*V0*BM*(BMder-1)**(-2)*(5+3*BMder*((V/V0)**(1./3.)-1)-3*(V/V0)**(1./3.))*np.exp(-(3./2.)*(BMder-1)*((V/V0)**(1/3.)-1))
#

import numpy as np
def f(x):
    return x**2 + 10*np.sin(x)
xdata = np.linspace(-10,10,num=20)
ydata = f(xdata) + np.random.randn(xdata.size)
#print "xdata",xdata
#print "ydata",ydata
def f22(x,a,b):
    return a*x**2 + b*np.sin(x)
guess = [2,2]
from scipy import optimize
params, params_covariance = optimize.curve_fit(f22,xdata,ydata,guess)
#print params

#print "OK"

#if path == None:
#    path = my.pwd()
energyfile = my.readfile("energy.dat")

atoms = None
xdata = []
ydata = []
for line in energyfile:
    words = line.split()
    if len(words) != 3:
        continue
    if atoms == None:
        atoms = int(words[2])
    volume = float(words[0])
    ene = float(words[1])
    xdata.append(volume)
    ydata.append(ene)
    #print atoms,volume,ene
#print "xdata",xdata
#print "ydata",ydata
def Vinet(V,E0,V0,BM,BMder):
    return E0+(4*BM*V0)/(BMder-1)**2-2*V0*BM*(BMder-1)**(-2)*(5+3*BMder*((V/V0)**(1./3.)-1)-3*(V/V0)**(1./3.))*np.exp(-(3./2.)*(BMder-1)*((V/V0)**(1/3.)-1))

guess = [-3.7,14,150,5]   ### ene//mean, vol//mean, 100,5
params, params_covariance = optimize.curve_fit(Vinet,xdata,ydata,guess)
paramsout = params
#print paramsout,len(paramsout),paramsout[2]*160.2176462
#print "OK"
paramsout[2] = paramsout[2]*160.2176462
#paramsout[0] = paramsout[0]/1000
print paramsout
#return paramsout
