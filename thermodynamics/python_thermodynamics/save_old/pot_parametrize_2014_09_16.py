#!/usr/bin/env python


import math
import sys
import os
import copy
import glob
import numpy as np
import pandas as pd
import itertools
import pickle
import time

from scipy                import optimize
from scipy.interpolate    import UnivariateSpline


import hesse as pot_energy
#reload(pot_energy)
import utils
import my_atom
reload(my_atom)
import crystal_generator
reload(crystal_generator)

np.set_printoptions(threshold=np.nan)  # print the whole array
np.set_printoptions(linewidth=240)    # print only 6 digist after .
np.set_printoptions(precision=3)    # print only 6 digist after .

pd.set_option('display.max_columns', 500)
pd.set_option('display.width', 200)
pd.set_option('display.max_columns', 25)
pd.set_option('display.max_colwidth', 6)


class fit_to_func( object ):
    def __init__(self, data = False, function = False, fixzeroat = False, weights = False ):
        ''' data can be either a 2d numpy array or a string containing the path to the file to import '''
        self.verbose    = False
        self.data       = data

        self.fixzeroat  = fixzeroat # nearest neighbor distance
        self.function   = function
        self.function_known = [ 'l', 'm', 'morse', 'mc1', 'poly' ]
        if self.function not in self.function_known:
            sys.exit("function "+function+" not in "+str(self.function_known))
        self.weights    = weights

        self.verbose = False
        self.parameters = False
        self.deltas = False
        self.deltasshifted = False
        self.polyfit = find_best_polyfit()

        self.get_data()
        self.datax = self.data[:,0]
        self.datay = self.data[:,1]
        #self.fitx = np.linspace(self.datax.min()-1.0,self.datax.max()+1.0, num=100)
        self.fit = False
        self.fitx = self.datax

        if type(self.weights) != bool:
            # bring weights on same grid as datax, datay
            self.weights_samegrid = utils.return_function2_on_grid_of_f1(self.data,self.weights)

        if self.function == 'poly':
            self.polyfit = find_best_polyfit(
                data = self.data, fixzeroat = self.fixzeroat)
            return
        else:
            self.fit_to_function()

    def get_data(self):
        if type(self.data) == str:
            self.data = np.loadtxt(self.data)
        return

    def fit_to_function(self):
        from scipy import optimize

        if self.function == 'l':
            self.parameters, self.function_covariance = \
                optimize.curve_fit(lambda r, eps: pot_energy.LJ_derivative(r, eps, self.fixzeroat), self.datax, self.datay, maxfev=1000)
            self.fity = pot_energy.LJ_derivative(self.fitx, 0.0342863, self.fixzeroat)

        if self.function == 'm' or self.function == 'morse':
            self.parameters, self.function_covariance = \
                    optimize.curve_fit(lambda r, De, aa: pot_energy.Morse_derivative(r, De, aa, self.fixzeroat), self.datax, self.datay, p0=[0.02, 1.96], maxfev=100000)
                #optimize.curve_fit(lambda r, De, aa: pot_energy.Morse_derivative(r, De, aa, self.fixzeroat), datax, datay,p0=[0.02, 1.96],sigma = weights[:,1], absolute_sigma = False, maxfev=10000000)
                #optimize.curve_fit(Morse_derivativeNN, data[:,0], data[:,1], maxfev=1000)
            self.parameters = np.array([self.parameters[0],self.parameters[1],self.fixzeroat])
            self.fity = pot_energy.Morse_derivative(self.fitx, *self.parameters)
            fitdelta = pot_energy.Morse_derivative(self.datax, *self.parameters) - self.datay

        if self.function == 'mc1':
            try:
                self.parameters, self.function_covariance = \
                    optimize.curve_fit(lambda r, De, aa, A, B: pot_energy.mc1_derivative(r, De, aa, self.fixzeroat, A, B), self.datax, self.datay, maxfev=100000)
                #optimize.curve_fit(lambda r, De, aa, A, B: pot_energy.mc1_derivative(r, De, aa, NN, A, B), datax, datay, sigma = weights[:,1], absolute_sigma = False, maxfev=1000000)
            except RuntimeError:
                return
            self.parameters = np.array([self.parameters[0],self.parameters[1],self.fixzeroat,self.parameters[2],self.parameters[3]])

            self.fity = pot_energy.mc1_derivative(self.fitx, *self.parameters)
            fitdelta = pot_energy.mc1_derivative(self.datax, *self.parameters) - self.datay

        self.fit    = np.transpose([self.datax,   self.fity    ])
        self.deltas = np.transpose([self.datax, - fitdelta])
        self.deltasshifted = np.transpose([self.datax-self.fixzeroat, -fitdelta])

        # now fit deltas with polyfit
        self.polyfit = find_best_polyfit(
                data = self.deltas, fixzeroat = self.fixzeroat)

        #polycoefsrealene,polycoefsrealforces,polydeltas = polyfit(
        #foldername=deltasfolderfit,
        #filename=filenameorig+"_add_",
        #        zeroat=NN,
        #        data=deltas,
        #        verbose2=False)


        #polyfit(
        #foldername=deltasfolderfitshifted,
        #filename=filenameorig+"_add_",
        #        zeroat=0.0,
        #        data=deltasshifted,
        #        verbose2=False)
        return


class find_best_polyfit( object ):
    def __init__(self, data = False, fixzeroat = False ):
        ''' data can be either a 2d numpy array or a string containing the path to the file to import '''
        self.data       = data

        self.fixzeroat  = fixzeroat # nearest neighbor distance

        self.verbose = False

        self.ordermax = False
        self.coefsrealene = False
        self.coefsrealforces = False
        self.deltas = False
        self.deltasshifted = False

        if type(self.data) != bool:
            self.get_data()
            self.datax = self.data[:,0]
            self.datay = self.data[:,1]
            self.fitx = np.linspace(self.datax.min()-1.0,self.datax.max()+1.0, num=100)
            self.polyfit()

    def get_data(self):
        if type(self.data) == str:
            self.data = np.loadtxt(self.data)
        return

    def polyfit(self):
        if self.ordermax == False:
            self.ordermax = self.datax.shape[0]/2.
            if self.ordermax > 9:
                self.ordermax = 9
        else:
            self.ordermax = int(self.ordermax)

        orderstry = np.arange(1,self.ordermax+1)  # 7 goes to only to order 6

        def findbestpolyfit(orderstry = self.ordermax, dataxminfit = False ,dataxmaxfit = False):
            self.ordermaxdelta = np.zeros((len(orderstry),18))
            self.ordermaxdelta[:] = np.nan
            for orderidx,order in enumerate(orderstry):
                self.ordermaxdelta[orderidx,0] = order
                order = int(order)
                #print ""
                if self.verbose:
                    print "order    :",order,type(int(order))
                    print "self.fixzeroat   :",self.fixzeroat,type(float(self.fixzeroat))
                coefs = polyfit_with_fixed_points(int(order),self.datax,self.datay,np.array([float(self.fixzeroat)]),np.array([0.0]))  # coefs for forces
                e0 = np.polyint(np.poly1d(coefs[::-1]))(float(self.fixzeroat))
                #coefs[0] = coefs[0] + e0
                #print "coefs",coefs


                # coefs:    coefs for forces
                # coefse:   coefs for energy
                coefse = np.polyint(coefs[::-1])[::-1]  # turn coefs for foces into coefs for energy
                coefse[0] = -e0
                e = np.poly1d(coefs[::-1])(float(self.fixzeroat))
                f = np.polyder(np.poly1d(coefse[::-1]))(float(self.fixzeroat))

                # get coefs as those would be read in (==used)
                # We want to do this to see how much of accuracy we loose due to exporting a
                # finite number of digits

                self.coefsstringforces = utils.list_to_string(coefs)
                self.coefsrealforces = utils.string_to_list(self.coefsstringforces)

                self.coefsstringene = utils.list_to_string(coefse)
                self.coefsrealene = utils.string_to_list(self.coefsstringene)

                if self.verbose:
                    print "coefs:",coefs
                er = np.poly1d(self.coefsrealene[::-1])(float(self.fixzeroat))
                fr = np.polyder(np.poly1d(self.coefsrealene[::-1]))(float(self.fixzeroat))
                function = np.polyder(np.poly1d(self.coefsrealene[::-1]))
                if self.verbose:
                    print "function:",function

                # calculate new x's and y's
                if type(dataxminfit) == bool:
                    dataxminfit = self.datax.min()-1
                if type(dataxmaxfit) == bool:
                    dataxmaxfit = self.datax.max()+1
                self.fitx = np.linspace(dataxminfit,dataxmaxfit, num=101)
                self.fity = function(self.fitx)

                self.fitdelta = function(self.datax) - self.datay
                #print "self.fitdeltamax.max():",'%f' % fitdeltamax
                #ft = np.poly1d(coefs[::-1])(float(self.fixzeroat))
                #et = np.polyint(np.poly1d(self.coefsrealene[::-1]))(float(self.fixzeroat))
                #f = np.poly1d(coefs[::-1])(float(self.fixzeroat))
                #e = np.polyint(np.poly1d(self.coefsrealene[::-1]))(float(self.fixzeroat))
                #coefse = np.polyint(np.poly1d(self.coefsrealene[::-1]))
                #print ""
                self.ordermaxdelta[orderidx,1] = self.fitdelta.max()
                self.ordermaxdelta[orderidx,2] = er - e
                self.ordermaxdelta[orderidx,3] = fr - f
                if False:
                    print "order:",order, '%f' % self.fitdelta.max(),"f:",f,"e:",e,"er:",er,"fr:",fr #,"coefs:",coefs
                #print "coefs:",coefs
                #print "coefse:",coefse
                #print "coefse:",coefse
            return
            #return ordermaxdelta,self.fitdelta,self.coefsstringene,self.coefsrealene,self.coefsstringforces,self.coefsrealforces,order,self.fitx,self.fity

        ###############################################################################################################################
        # try fit with several orders 1 ... 9  (orderstry)
        ###############################################################################################################################
        #self.ordermaxdelta,self.fitdelta,self.coefsstringene,self.coefsrealene,self.coefsstringforces,self.coefsrealforces,order,self.fitx,self.fity =   findbestpolyfit(self.datax,self.datay,self.fixzeroat,orderstry,dataxminfit,dataxmaxfit)
        findbestpolyfit(orderstry = orderstry)
        if self.verbose:
            print "(order)      (dforce)  (er - e)  (fr - f)"
            np.set_printoptions(linewidth=240)    # print only 6 digist after .
            print self.ordermaxdelta
            print self.ordermaxdelta[:,1].argmin()
            print self.ordermaxdelta.argmin()
        #orderstry = np.arange(1,ordermax+1)  # 7 goes to only to order 6
        ###############################################################################################################################
        # now get the best fitted function (np.arange(ordermaxdelta[:,1].argmin()+1,ordermaxdelta[:,1].argmin()+2))
        ###############################################################################################################################
        #self.ordermaxdelta,self.fitdelta,self.coefsstringene,self.coefsrealene,self.coefsstringforces,self.coefsrealforces,order,self.fitx,self.fity = findbestpolyfit(self.datax,self.datay,self.fixzeroat,
        #        np.arange(self.ordermaxdelta[:,1].argmin()+1,seff.poly_ordermaxdelta[:,1].argmin()+2),
        #        self.dataxminfit,dataxmaxfit)
        findbestpolyfit(    orderstry = np.arange(self.ordermaxdelta[:,1].argmin()+1,self.ordermaxdelta[:,1].argmin()+2))

        #                    dataxminfit = self.dataxminfit,
        #                    dataxmaxfit = dataxmaxfit)


        if self.verbose:
            print self.ordermaxdelta
            print "-------"
        self.deltas = np.transpose([self.datax, - self.fitdelta])
        self.deltasshifted = np.transpose([self.datax-self.fixzeroat, - self.fitdelta])

        ##############################################################################################################################
        # Save file to disk
        ##############################################################################################################################
        if False: # part for saving files
            if type(filename) != bool:
                deltasfolder = foldername+"/deltas/"
                deltasfoldershifted = foldername+"/deltasshifted/"
                if len(deltasfolder.split("/")) > 16 or len(deltasfolder.split("/")) < 8:
                    print "in polyfit: len:",len(deltasfolder.split("/"))
                    print "in polyfit: foldername:",foldername
                    print "in polyfit: deltasfolder:",deltasfolder
                    sys.exit("in polyfit: deltasfolder is wrong")

                deltasfolderfilename = deltasfolder+filenameorig+"_poly"+str(int(ordermax))+"_"+str(int(order))+"_delta___"+self.coefsstringene
                deltasfolderfilenameshifted = deltasfoldershifted+filenameorig+"_poly"+str(int(ordermax))+"_"+str(int(order))+"_delta___"+self.coefsstringene
                if len(deltasfolderfilename.split("/")) > 16:
                    print "in polyfit: deltasfolderfilename:",deltasfolderfilename
                    sys.exit("in polyfit: deltasfolder is wrong")

                if os.path.isdir(deltasfolder) != True: os.makedirs(deltasfolder)
                if os.path.isdir(deltasfoldershifted) != True: os.makedirs(deltasfoldershifted)
                np.savetxt(
                        filename+"_poly"+str(int(ordermax))+"_"+str(int(order))+"_fit___"+self.coefsstringene,
                        np.transpose([self.fitx, self.fity]),
                        fmt='%.18f',
                        delimiter='  ')   # X is an array
                np.savetxt(
                        # filenameorig seems always to be
                        #deltasfolder+filenameorig+"_poly"+str(int(ordermax))+"_"+str(int(order))+"_delta___"+self.coefsstringene,
                        deltasfolderfilename,
                        deltas,
                        fmt='%.18f',
                        delimiter='  ')   # X is an array
                np.savetxt(
                        #deltasfolder+filenameorig+"_poly"+str(int(ordermax))+"_"+str(int(order))+"_delta___"+self.coefsstringene,
                        deltasfolderfilenameshifted,
                        self.deltasshifted,
                        fmt='%.18f',
                        delimiter='  ')   # X is an array

        #if return_deltasfilename:
        #    return deltasfilename
        #else:
        #    return self.coefsrealene,self.coefsrealforces,self.poly_deltas
        return




def get_fit_forces_to_pot(filename = False, NN = False, pot = False, data = False, weights = False, foldername = "", verbose = True):
    ''' fits forces to morse, lj, mc1 '''
    filenameorig = copy.copy(filename)
    filename = foldername+filename
    if type(data) == bool:
        if type(filename) == bool:
            sys.exit("please provide a filename")
        if os.path.isfile(filename) != True:
            sys.exit("file "+filename+" does not exist")
        data=np.loadtxt(filename)
    else:
        data = data
    #alat =  float(filename.split("_")[-1])
    #NN = alat/np.sqrt(2.)
    if type(NN) == bool:
        if pot == '135' or pot == '1357' or pot == '13579':
            pass
        else:
            sys.exit("please provide NN")
    if type(NN) == str:
        NN= float(NN)
    if type(pot) == False:
        sys.exit("please define a pot")

    datax=data[:,0]
    datay=data[:,1]

    if type(weights) != bool:
        # bring weights on same grid as datax, datay
        weights = utils.return_function2_on_grid_of_f1(data,weights)

    from scipy import optimize
    if pot == 'l':
        parameters, function_covariance = \
            optimize.curve_fit(lambda r, eps: pot_energy.LJ_derivative(r, eps, NN), datax, datay, maxfev=1000)
    if pot == 'm' or pot == 'morse':
        #print 'datx:'
        #print datax.shape
        #print datay.shape
        parameters, function_covariance = \
                optimize.curve_fit(lambda r, De, aa: pot_energy.Morse_derivative(r, De, aa, NN), datax, datay, p0=[0.02, 1.96], maxfev=100000)
            #optimize.curve_fit(lambda r, De, aa: pot_energy.Morse_derivative(r, De, aa, NN), datax, datay,p0=[0.02, 1.96],sigma = weights[:,1], absolute_sigma = False, maxfev=10000000)
            #optimize.curve_fit(Morse_derivativeNN, data[:,0], data[:,1], maxfev=1000)
        parameters = np.array([parameters[0],parameters[1],NN])
    if pot == 'mc1':
        try:
            parameters, function_covariance = \
                optimize.curve_fit(lambda r, De, aa, A, B: pot_energy.mc1_derivative(r, De, aa, NN, A, B), datax, datay, maxfev=100000)
            #optimize.curve_fit(lambda r, De, aa, A, B: pot_energy.mc1_derivative(r, De, aa, NN, A, B), datax, datay, sigma = weights[:,1], absolute_sigma = False, maxfev=1000000)
        except RuntimeError:
            return False, False, False, False, False, False
        parameters = np.array([parameters[0],parameters[1],NN,parameters[2],parameters[3]])

    fitx = np.linspace(datax.min()-1.0,datax.max()+1.0, num=100)
    if pot == 'l':
        #fity = LJ_derivative(fitx, parameters[0], NN)
        fity = pot_energy.LJ_derivative(fitx, 0.0342863, NN)
    if pot == 'm' or pot == 'morse':
        fity = pot_energy.Morse_derivative(fitx, parameters[0], parameters[1], NN)
        fitdelta = pot_energy.Morse_derivative(datax, parameters[0], parameters[1], NN) - datay
        #print "fd:",fitdelta
    if pot == 'mc1':
        fity = pot_energy.mc1_derivative(fitx, parameters[0], parameters[1], NN, \
                parameters[2], parameters[3])
        fitdelta = pot_energy.mc1_derivative(datax, parameters[0], parameters[1], NN, \
                parameters[2], parameters[3]) - datay


    #print "saving..."
    deltas = np.transpose([datax, -fitdelta])
    deltasshifted = np.transpose([datax-NN, -fitdelta])

    deltasfolder = foldername+"/deltas/"
    deltasfoldershifted = foldername+"/deltasshifted/"

    datashiftedfolder = foldername+"/datashifted"
    datashifted = np.transpose([datax-NN, datay])

    deltasfolderfit = foldername+"/deltasfit/"
    deltasfolderfitshifted = foldername+"/deltasfitshifted/"

    if len(deltasfolder.split("/")) > 14:
        print "in get_fit_forces_to_pot: deltasfolder in :",deltasfolder,len(deltasfolder)
        sys.exit("in get_fit_forces_to_pot: deltasfolder is wrong")
    if os.path.isdir(deltasfolder) != True: os.makedirs(deltasfolder)
    if os.path.isdir(deltasfoldershifted) != True: os.makedirs(deltasfoldershifted)
    if os.path.isdir(deltasfolderfit) != True: os.makedirs(deltasfolderfit)
    if os.path.isdir(datashiftedfolder) != True: os.makedirs(datashiftedfolder)

    deltasfolderfilename = deltasfolder+filenameorig+"_fit_delta___"+str(parameters[0])[:8]
    if len(deltasfolderfilename.split("/")) > 14:
        print "in get_fit_forces_to_pot: deltasfolderfilename:",deltasfolderfilename
        sys.exit("in get_fit_forces_to_pot: deltasfolder is wrong")

    polycoefsrealene,polycoefsrealforces,polydeltas = polyfit(
    foldername=deltasfolderfit,
    filename=filenameorig+"_add_",
            zeroat=NN,
            data=deltas,
            verbose2=False)

    polyfit(
    foldername=deltasfolderfitshifted,
    filename=filenameorig+"_add_",
            zeroat=0.0,
            data=deltasshifted,
            verbose2=False)

    if pot == 'l':
        np.savetxt(
            filename+"_fit_lj_"+str(parameters[0])[:8]
            +"_"+str(NN)[:8],
                np.transpose([fitx, fity]),
                fmt='%.18f',
                delimiter='  ')   # X is an array

    if pot == 'm' or pot == 'morse':
        if type(filename) != bool:
            np.savetxt(
                filename+"_fit___"+str(parameters[0])[:8]
                +"_"+str(parameters[1])[:8]
                +"_"+str(NN)[:8],
                    np.transpose([fitx, fity]),
                    fmt='%.18f',
                    delimiter='  ')   # X is an array
            np.savetxt(
                datashiftedfolder+filenameorig+"_fit___"+str(parameters[0])[:8]
                +"_"+str(parameters[1])[:8]
                +"_"+str(NN)[:8],
                    datashifted,
                    fmt='%.18f',
                    delimiter='  ')   # X is an array

            np.savetxt(
                deltasfolder+filenameorig+"_fit_delta___"+str(parameters[0])[:8]
                +"_"+str(parameters[1])[:8]
                +"_"+str(NN)[:8],
                    deltas,
                    fmt='%.18f',
                    delimiter='  ')   # X is an array
            np.savetxt(
                deltasfoldershifted+filenameorig+"_fit_delta___"+str(parameters[0])[:8]
                +"_"+str(parameters[1])[:8]
                +"_"+str(NN)[:8],
                    deltasshifted,
                    fmt='%.18f',
                    delimiter='  ')   # X is an array

    if pot == 'mc1':
        #print "parss mc1:",parameters
        #parameters = np.array([parameters[0],parameters[1],NN,parameters[2],parameters[3]])
        #print "parss mc1:",parameters
        #print NN
        #sys.exit()
        if type(filename) != bool:
            np.savetxt(
                filename+"_fit___"+str(parameters[0])[:8]
                +"_"+str(parameters[1])[:8]
                +"_"+str(NN)[:8]
                +"_"+str(parameters[2])[:8]
                +"_"+str(parameters[3])[:8],
                    np.transpose([fitx, fity]),
                    fmt='%.18f',
                    delimiter='  ')   # X is an array
            #print "filenameorig:",filenameorig

            np.savetxt(
                deltasfolder+filenameorig+"_fit_delta___"+str(parameters[0])[:8]
                +"_"+str(parameters[1])[:8]
                +"_"+str(NN)[:8]
                +"_"+str(parameters[2])[:8]
                +"_"+str(parameters[3])[:8],
                    deltas,
                    fmt='%.18f',
                    delimiter='  ')   # X is an array
            np.savetxt(
                deltasfoldershifted+filenameorig+"_fit_delta___"+str(parameters[0])[:8]
                +"_"+str(parameters[1])[:8]
                +"_"+str(NN)[:8]
                +"_"+str(parameters[2])[:8]
                +"_"+str(parameters[3])[:8],
                    deltasshifted,
                    fmt='%.18f',
                    delimiter='  ')   # X is an array
    return parameters, deltas, deltasshifted, polycoefsrealene,polycoefsrealforces,polydeltas



def polyfit_with_fixed_points(n, x, y, xf, yf) :
    #print "n:",n
    #print "x:",x
    #print "y:",y
    #print "xf:",xf,type(xf)
    #print "yf:",yf,type(yf)
    #print n,len(xf),type(n),type(len(xf))
    mat = np.empty((n + 1 + len(xf),) * 2)
    vec = np.empty((n + 1 + len(xf),))
    x_n = x**np.arange(2 * n + 1)[:, None]
    yx_n = np.sum(x_n[:n + 1] * y, axis=1)
    x_n = np.sum(x_n, axis=1)
    idx = np.arange(n + 1) + np.arange(n + 1)[:, None]
    mat[:n + 1, :n + 1] = np.take(x_n, idx)
    xf_n = xf**np.arange(n + 1)[:, None]
    mat[:n + 1, n + 1:] = xf_n / 2
    mat[n + 1:, :n + 1] = xf_n.T
    mat[n + 1:, n + 1:] = 0
    vec[:n + 1] = yx_n
    vec[n + 1:] = yf
    #print "mat:",mat
    #print "vec:",vec
    params = np.linalg.solve(mat, vec)
    return params[:n + 1]

def polyfit(filename = False, zeroat = False, ordermax = False, data = False, dataxminfit= False, dataxmaxfit = False, verbose = False, verbose2=False, foldername = "", return_deltasfilename = False ):
    ''' fits data in filename to polynomial of order x '''
    filenameorig = copy.copy(filename)
    #print "foldername>>>:",foldername
    #print "filename>>>:",filename
    if type(foldername) == bool and type(filename) == bool:
        filename = foldername
    else:
        filename = foldername+filename
    if type(data) == bool:
        if type(filename) == bool:
            sys.exit("please provide a filename or data")
        if os.path.isfile(filename) != True:
            sys.exit("file "+filename+" does not exist")
        data=np.loadtxt(filename)
    else:
        data = data

    if type(zeroat) == bool:
        sys.exit("arguments 1: valid Filename; 2: zeroat;")
    datax=data[:,0]
    datay=data[:,1]
    #print "#################: datax:",zeroat
    #print datax
    #print "#################: datay:",zeroat
    #print datay

    if ordermax == False:
        ordermax = datax.shape[0]/2.
        if ordermax > 9:
            ordermax = 9
    else:
        ordermax = int(ordermax)
    #print "ordermax:",ordermax

    orderstry = np.arange(1,ordermax+1)  # 7 goes to only to order 6

    def findbestpolyfit(datax,datay,zeroat,orderstry,dataxminfit,dataxmaxfit):
        #ordermaxdelta = np.zeros((ordermax+1,20))
        ordermaxdelta = np.zeros((len(orderstry),18))
        ordermaxdelta[:] = np.nan
        for orderidx,order in enumerate(orderstry):
            ordermaxdelta[orderidx,0] = order
            order = int(order)
            #print ""
            if verbose:
                print "order    :",order,type(int(order))
                print "zeroat   :",zeroat,type(float(zeroat))
            coefs = polyfit_with_fixed_points(int(order),datax,datay,np.array([float(zeroat)]),np.array([0.0]))  # coefs for forces
            e0 = np.polyint(np.poly1d(coefs[::-1]))(float(zeroat))
            #coefs[0] = coefs[0] + e0
            #print "coefs",coefs


            # coefs:    coefs for forces
            # coefse:   coefs for energy
            coefse = np.polyint(coefs[::-1])[::-1]  # turn coefs for foces into coefs for energy
            coefse[0] = -e0
            e = np.poly1d(coefs[::-1])(float(zeroat))
            f = np.polyder(np.poly1d(coefse[::-1]))(float(zeroat))

            # get coefs as those would be read in (==used)
            # We want to do this to see how much of accuracy we loose due to exporting a
            # finite number of digits

            coefsstringforces = utils.list_to_string(coefs)
            coefsrealforces = utils.string_to_list(coefsstringforces)

            coefsstringene = utils.list_to_string(coefse)
            coefsrealene = utils.string_to_list(coefsstringene)

            if verbose:
                print "coefs:",coefs
            er = np.poly1d(coefsrealene[::-1])(float(zeroat))
            fr = np.polyder(np.poly1d(coefsrealene[::-1]))(float(zeroat))
            function = np.polyder(np.poly1d(coefsrealene[::-1]))
            if verbose:
                print "function:",function

            # calculate new x's and y's
            if type(dataxminfit) == bool:
                dataxminfit = datax.min()-1
            if type(dataxmaxfit) == bool:
                dataxmaxfit = datax.max()+1
            fitx = np.linspace(dataxminfit,dataxmaxfit, num=101)
            fity = function(fitx)

            fitdelta = function(datax) - datay
            fitdeltamax = abs(fitdelta).max()
            #print "fitdeltamax:",'%f' % fitdeltamax
            #ft = np.poly1d(coefs[::-1])(float(zeroat))
            #et = np.polyint(np.poly1d(coefsrealene[::-1]))(float(zeroat))
            #f = np.poly1d(coefs[::-1])(float(zeroat))
            #e = np.polyint(np.poly1d(coefsrealene[::-1]))(float(zeroat))
            #coefse = np.polyint(np.poly1d(coefsrealene[::-1]))
            #print ""
            ordermaxdelta[orderidx,1] = fitdeltamax
            ordermaxdelta[orderidx,2] = er - e
            ordermaxdelta[orderidx,3] = fr - f
            if False:
                print "order:",order, '%f' % fitdeltamax,"f:",f,"e:",e,"er:",er,"fr:",fr #,"coefs:",coefs
            #print "coefs:",coefs
            #print "coefse:",coefse
            #print "coefse:",coefse
        return ordermaxdelta,fitdelta,coefsstringene,coefsrealene,coefsstringforces,coefsrealforces,order,fitx,fity




    #print "orderstry:",orderstry
    ###############################################################################################################################
    # try fit with several orders 1 ... 9  (orderstry)
    ###############################################################################################################################
    ordermaxdelta,fitdelta,coefsstringene,coefsrealene,coefsstringforces,coefsrealforces,order,fitx,fity = findbestpolyfit(datax,datay,zeroat,orderstry,dataxminfit,dataxmaxfit)
    if verbose2:
        print "(order)      (dforce)  (er - e)  (fr - f)"
        np.set_printoptions(linewidth=240)    # print only 6 digist after .
        print ordermaxdelta
        print "xxxxxxyy"
        print ordermaxdelta[:,1].argmin()
        print ordermaxdelta.argmin()
        print "-------"
    #orderstry = np.arange(1,ordermax+1)  # 7 goes to only to order 6
    ###############################################################################################################################
    # now get the best fitted function (np.arange(ordermaxdelta[:,1].argmin()+1,ordermaxdelta[:,1].argmin()+2))
    ###############################################################################################################################
    ordermaxdelta,fitdelta,coefsstringene,coefsrealene,coefsstringforces,coefsrealforces,order,fitx,fity = findbestpolyfit(datax,datay,zeroat,
            np.arange(ordermaxdelta[:,1].argmin()+1,ordermaxdelta[:,1].argmin()+2),
            dataxminfit,dataxmaxfit)
    if verbose2:
        print ordermaxdelta
        print "-------"
    deltas = np.transpose([datax, -fitdelta])
    deltasshifted = np.transpose([datax-zeroat, -fitdelta])
    deltamax = fitdelta



    if type(filename) != bool:
        deltasfolder = foldername+"/deltas/"
        deltasfoldershifted = foldername+"/deltasshifted/"
        if len(deltasfolder.split("/")) > 16 or len(deltasfolder.split("/")) < 8:
            print "in polyfit: len:",len(deltasfolder.split("/"))
            print "in polyfit: foldername:",foldername
            print "in polyfit: deltasfolder:",deltasfolder
            sys.exit("in polyfit: deltasfolder is wrong")

        deltasfolderfilename = deltasfolder+filenameorig+"_poly"+str(int(ordermax))+"_"+str(int(order))+"_delta___"+coefsstringene
        deltasfolderfilenameshifted = deltasfoldershifted+filenameorig+"_poly"+str(int(ordermax))+"_"+str(int(order))+"_delta___"+coefsstringene
        if len(deltasfolderfilename.split("/")) > 16:
            print "in polyfit: deltasfolderfilename:",deltasfolderfilename
            sys.exit("in polyfit: deltasfolder is wrong")

        if os.path.isdir(deltasfolder) != True: os.makedirs(deltasfolder)
        if os.path.isdir(deltasfoldershifted) != True: os.makedirs(deltasfoldershifted)
        np.savetxt(
                filename+"_poly"+str(int(ordermax))+"_"+str(int(order))+"_fit___"+coefsstringene,
                np.transpose([fitx, fity]),
                fmt='%.18f',
                delimiter='  ')   # X is an array
        np.savetxt(
                # filenameorig seems always to be
                #deltasfolder+filenameorig+"_poly"+str(int(ordermax))+"_"+str(int(order))+"_delta___"+coefsstringene,
                deltasfolderfilename,
                deltas,
                fmt='%.18f',
                delimiter='  ')   # X is an array
        np.savetxt(
                #deltasfolder+filenameorig+"_poly"+str(int(ordermax))+"_"+str(int(order))+"_delta___"+coefsstringene,
                deltasfolderfilenameshifted,
                deltasshifted,
                fmt='%.18f',
                delimiter='  ')   # X is an array
    #return fitdelta
    #print "cfd:",coefsrealene
    #print "2;",coefsrealforces
    #print "3;",deltas
    if return_deltasfilename:
        return deltasfilename
    else:
        return coefsrealene,coefsrealforces,deltas




def getforcepot(longvec,pot,potparam,pot2 = False):
    usepot=['nan','nan','nan','nan','nan','nan','nan']
    if pot == 'm': usepot[0] = "Morse"
    if pot == 'l': usepot[0] = "LJ"
    if pot == 'i': usepot[0] = "inversepot"
    if pot == 'mc1': usepot[0] = "mc1"
    for idx,p in enumerate(potparam):
        usepot[idx+1] = float(p)

    longvecnorm=np.linalg.norm(longvec)
    getfrompot2 = False
    if type(pot2) != bool:
        if usepot[0] == "Morse" or usepot[0] == "mc1":
            NNdist = float(pot2[3])
            #print "NNdist:",NNdist,longvecnorm
            if longvecnorm > NNdist:
                getfrompot2 = True
        else:
            sys.exit("pot2 only m or mc1")

    if usepot[0] == "inversepot":
        elong = pot_energy.inversepot(longvecnorm,usepot[1],usepot[2],usepot[3])
        flongnorm = pot_energy.inversepot_derivative(longvecnorm,usepot[1],usepot[2],usepot[3])

    if usepot[0] == "Morse":
        if getfrompot2 == False:
            elong = pot_energy.Morse(longvecnorm,usepot[1],usepot[2],usepot[3])
            flongnorm = pot_energy.Morse_derivative(longvecnorm,usepot[1],usepot[2],usepot[3])
        else:
            if pot2[0] != 'm':
                sys.exit("not m")
            elong =                Morse(longvecnorm,float(pot2[1]),float(pot2[2]),float(pot2[3]))
            flongnorm = Morse_derivative(longvecnorm,float(pot2[1]),float(pot2[2]),float(pot2[3]))


    if usepot[0] == "mc1" or usepot[0] == 'mc1':
        #print "usepot:",usepot
        if getfrompot2 == False:
            elong = pot_energy.mc1(longvecnorm,usepot[1],usepot[2],usepot[3],usepot[4],usepot[5])
            flongnorm = pot_energy.mc1_derivative(longvecnorm,usepot[1],usepot[2],usepot[3],usepot[4],usepot[5])
        else:
            if pot2[0] != 'mc1':
                print "pot2:",pot2
                sys.exit("not mc1")
            elong =                pot_energy.mc1(longvecnorm,float(pot2[1]),float(pot2[2]),float(pot2[3]),float(pot2[4]),float(pot2[5]))
            flongnorm = pot_energy.mc1_derivative(longvecnorm,float(pot2[1]),float(pot2[2]),float(pot2[3]),float(pot2[4]),float(pot2[5]))


    if usepot[0] == "LJ":
        elong = pot_energy.LJ(longvecnorm,usepot[1],usepot[2])
        flongnorm = pot_energy.LJ_derivative(longvecnorm,usepot[1],usepot[2])

    # flong ist die longitudinale Kraft
    flong = longvec/longvecnorm*flongnorm   # der volle vektor == voller longitudinale vektor
    if flong[0] == -0.0:
        flong[0] = 0.0
    if flong[1] == -0.0:
        flong[1] = 0.0
    if flong[2] == -0.0:
        flong[2] = 0.0
    #print np.around(flong[0],decimals=7),np.around(flong[1],decimals=7),np.around(flong[2],decimals=7)
    return flong


def dos_to_pot(dos,tmelt):
    ''' gets dos (2d numpy array) and tmelt in Kelvin and returns the potential '''
    pot = np.copy(dos)
    pot[:,1] = -np.log(dos[:,1])*tmelt*0.086173423
    pot[:,1] = pot[:,1]+abs(pot[:,1].min())
    #pot[:,1] = pot_
    return pot

class forcesneighbors( object ):
    '''
    for fcc:
    form quer we want to get    long 1NN (  this appears to be pure compared to the
                                            x-direction which contains also an inplane
                                            part)

    from xdir we want to get    long 1NN (this we can get also quer disp might be better)
                                long 2NN
                                tox 1NN
                                tox 2NN

    analyzes the forces of the different neighbor shells
    this class is designed for a certain alat which has to be defined beforehand
    f[1]    : all forces of first alat WRONG NOW
    f[1,1]  : all forces of first alat and first displacement WRONG NOW
    f[1,1,1]: all forces fo first alat and first displacement in x,y,z direction WRONG NOW
    '''
    def __init__(self):
        self.a       = False      # alat
        self.disp    = False      # defines if the current displacement is a disp in 'xdir' or 'quer'
                                 # quer direction
        self.sc      = False
        self.folder  = False      # list of folders
        self.dnorm   = False      # norm of displacement (from origin) (should be of len a)
        self.dvec    = False      # vec of displacement (from origin) (should be of len a)
        self.rcar    = False
        self.f       = False
        self.p       = False      # positions (should be of lenght a)
        self.struct  = False      # fcc, bcc ...
        self.element = False
        self.nnforcesaddname = False    # Fase means to add nothing to nnforces
        self.initialize_dataframe()
        self.pkl = False
        self.fit_file = False
        self.verbose = False
        self.parametrize_npz         = False
        self.parametrize_data_pkl    = False
        self.parametrize_fml_npz = False
        self.parametrize_fml_npz_path              = False

        self.fml                    = False
        self.fml_u1nn_pottype       = False
        self.fml_u1nn_potparam      = False
        self.fml_u1nn_potadd        = False
        self.fml_u1nn_potaddparam   = False
        self.FML                    = False
        self.angles                 = False
        self.fprojlongold           = False

    def initialize_dataframe(self):
        ''' initializes an empty dataframe which will be filled afterwards '''

        self.datacolumns = [ 'shells', 'contr', 'disp' ,'range' ]

        self._contr  = [ 'lon', 'tox', 'toy','toz', 'tix', 'tiy', 'tiz' ] ## -> dies ist die column coordinate
        self._contr  = [ 'lon', 'tox' ] ## -> dies ist die column coordinate
        self._diff   = self._contr
        self._shells = [ 1, 2, 3 ]
        self._func   = [ 'morse', 'mc1' , 'poly']
        self._range  = [ 'alle', 'dos', 'doscut' ] # dont take all, gets confused with pandas
        self._disp   = [ 'xdir', 'quer', 'midd', '3nnd', '4nnd' ]


        self.datarows= [ 'points', 'force', 'file', 'dos', 'pot', 'pointsshifted', 'dosshifted', 'potshifted' ] + self._func + [ i+'fit' for i in self._func]

        #self.datarows= [ 'points', 'morse' ]
        #values = [ [ 'lon',
        #          lon
        #          lonadd
        #          lon2
        #          lon2add
        #          to{x,y,z}
        #          ti{x,y,z}
        #                      {1,2,3,4...}
        #                                  'm'
        #                                  'mc1'
        #                                  'poly'
        #                                           'alle'
        #                                           'dos'
        #                                           'doscut'
        #                                                   'xdir'
        #                                                   'quer'

        #########################################################################################
        # self.data
        a = [ eval("self._"+i) for i in self.datacolumns ]
        #print "a:",a
        b = list(itertools.product(*a))
        #print "b:",b
        #print "self.datacolumns:",self.datacolumns
        index = pd.MultiIndex.from_tuples(b, names=self.datacolumns)    # hirarchie / how many layers
        #print "index:",index

        self.data = pd.DataFrame(columns = index)
        for i,iname in enumerate(self.datarows):
            self.data.loc[iname] = [ False for i in range(len(b)) ]
        self.data = self.data.astype(object)
        #print self.data
        #return


        #########################################################################################
        # self.datafit
        self.datacolumnsfit = self.datacolumns + [ 'func' ]
        self.datarowsfit= [ 'diff', 'file', 'diffshifted', 'polybeste', 'polybestf', 'polybestdiff', 'polybestdiffshifted',  'splinedata', 'splinediff', 'splineshifteddata', 'splineshifteddiff' ]
        a = [ eval("self._"+i) for i in self.datacolumnsfit ]
        b = list(itertools.product(*a))
        index = pd.MultiIndex.from_tuples(b, names=self.datacolumnsfit)    # hirarchie / how many layers
        self.fit = pd.DataFrame(columns = index)
        for i,iname in enumerate(self.datarowsfit):
            self.fit.loc[iname] = [ False for i in range(len(b)) ]
        self.fit = self.fit.astype(object)

        ###########################################################################################################
        # Example Multidim
        ###########################################################################################################
        #from itertools import product

        #np.random.seed(1)

        #team_names = ['Yankees', 'Mets', 'Dodgers']
        #jersey_numbers = [35, 71, 84]
        #game_numbers = [1, 2]


        #observer_names = ['Bill', 'John', 'Ralph']
        #observer_names = self._contr
        #observation_types = ['Speed', 'Strength']
        #observation_types = self._shells
        #disp = self._disp
        #rang = self._range
        #func = self._func

        #row_indices = list(product(team_names, jersey_numbers, game_numbers, observer_names, observation_types, disp, rang, func))
        #observation_values = np.random.randn(len(row_indices))

        #tns, jns, gns, ons, ots, disp, rang, func = zip(*row_indices)

        #data = pd.DataFrame({'team': tns, 'jersey': jns, 'game': gns, 'diff': ons, 'shells': ots, 'disp': disp, 'rang': rang, 'func': func, 'value': observation_values})

        #data = data.set_index(['team', 'jersey', 'game', 'diff', 'shells', 'disp', 'rang', 'func'])

        ## set the columns
        #data = data.unstack(['diff', 'shells', 'disp', 'rang', 'func' ])
        #data.columns = data.columns.droplevel(0)

        #print data
        return self.data

    def loadforces(self, folder = False, a = False):
        ''' loads data from folder (has to be present) '''
        if type(folder) == bool:
            sys.exit("pot_parametrize.forcesneighbors.loadforces: please specify the folder where to load data")
        if os.path.isdir(folder) != True:
            print "folder:",folder
            sys.exit("pot_parametrize.forcesneighbors.loadforces: specifies folder does not exist!")
        self.folderparam = folder
        if self.parametrize_npz == False:
            if os.path.isfile(self.folderparam+"/nnforces.npz") == True:
                self.loadData(self.folderparam,"nnforces")
                return

        self.eqcoords = np.loadtxt(self.folderparam+"/EqCoords_direct")
        self.atoms = self.eqcoords.shape[0]
        if self.atoms == 54: self.struct = 'bcc'
        if self.atoms == 32: self.struct = 'fcc'
        if self.atoms == 108: self.struct = 'fcc'
        if type(self.struct) == bool:
            sys.exit("self.struct not known")

        #self.folderparam = os.getcwd()
        self.element = my_atom.get_element_from_path(self.folderparam)
        self.atom = my_atom.atom([self.element])
        if type(self.a) == bool:
            self.a = a
        if type(self.a) == bool:
            sys.exit("when loading loadforces you have to specify \"a\"")
        self.folder = utils.lsn(self.folderparam+"/"+str(self.a)+"Ang_"+"*")
        if len(self.folder) == 0:
            print "self.a",self.a
            sys.exit("Error in loadforces: no folder found with name: "+str(self.a)+"Ang_"+"*")

        self.f = np.zeros((np.array(self.folder).shape[0],self.atoms,3))
        self.p = np.zeros((np.array(self.folder).shape[0],self.atoms,3))
        self.pr = np.zeros((np.array(self.folder).shape[0],self.atoms,3))  # removed mapping
        self.e = np.zeros(np.array(self.folder).shape[0])
        self.dstring = range(len(self.folder))
        self.dvec = np.zeros((np.array(self.folder).shape[0],3))
        self.dnorm = np.zeros((np.array(self.folder).shape[0]))


        for folderidx,folder in enumerate(self.folder):
            #print self.f[folderidx]
            self.dstring[folderidx] = utils.string_to_num_list(folder)[-1]
            self.dvec[folderidx] = self.p[folderidx][0]



            if os.path.isfile(folder+"/forces_OUTCAR") != True:
                utils.run2("OUTCAR_forces-last-ARRAY.sh "+folder+" > "+folder+"/forces_OUTCAR")
            self.f[folderidx] = np.loadtxt(folder+"/forces_OUTCAR")
            print str(folderidx)+"/"+str(len(self.folder)-1),folder,self.dstring[folderidx],self.f[folderidx][0]

            if os.path.isfile(folder+"/cartesian_coords") != True:
                utils.run2("OUTCAR_positions-last-ARRAY.sh "+folder+" > "+folder+"/cartesian_coords")
            self.p[folderidx] = np.loadtxt(folder+"/cartesian_coords")

            if os.path.isfile(folder+"/ene_free_last") != True:
                utils.run2("OUTCAR_ene-free-last.sh "+folder+" > "+folder+'/ene_free_last')
            #self.e[folderidx] = np.loadtxt(folder+"/ene_free_last")  # convert to meV per atom
            self.e[folderidx] = np.loadtxt(folder+"/ene_free_last")*1000./(self.atoms-1)  # convert to meV per atom
            #self.energymevtopy = np.sum(self.energytrantopy)/2.*1000/(self.numberofatoms-1)

            if os.path.isfile(folder+"/cell") != True:
                utils.run2("OUTCAR_cell-last-cartesian-ARRAY.sh "+folder+" > "+folder+"/cell")
            self.sc = np.loadtxt(folder+"/cell")

            # create also positions with removedmapping into original cell (needs crysal0)
            if folderidx == 0:
                crystal0 = crystal_generator.crystal()
                crystal0.load_positions_cell(cell = self.sc, coord_rrel = self.eqcoords)
                copy_of_crystal0_instance = copy.deepcopy(crystal0)
                self.nndist = copy_of_crystal0_instance.get_NNlist(0, 1,
                cell = copy_of_crystal0_instance.cellvec,
                coord_rrel = copy_of_crystal0_instance.rrel,
                return_NNdist = True)
            crystal1 = crystal_generator.crystal()
            crystal1.load_positions_cell(cell = self.sc, coord_cart = self.p[folderidx])
            crystalremovedmapping = crystal_generator.remove_mapping_into_originalcell()
            crystalremovedmapping.remove_mapping(crystal1,crystal0)
            self.pr[folderidx] = crystalremovedmapping.rcar
            self.dvec[folderidx] = self.pr[folderidx][0]
            self.dnorm[folderidx] = np.linalg.norm(self.pr[folderidx][0])


        # get other variables
        allposforces = np.zeros((self.p.shape[0]*self.p.shape[1]+self.atoms,6))  # + self.atoms um die erste (unausgelenkte struktur) hinzuzufuegen
        allpos = self.p.flatten().reshape(self.p.shape[0]*self.p.shape[1],3)
        allforces = self.f.flatten().reshape(self.f.shape[0]*self.f.shape[1],3)
        allposforces[self.atoms:,0:3] = allpos
        allposforces[self.atoms:,3:6] = allforces


        # ene_vs_disp
        self.positions = allposforces.reshape(self.p.shape[0]*self.p.shape[1]+self.atoms,6)
        self.ene_vs_disp = np.array((self.dstring,self.e)).transpose()

        # POSITIONs
        minidx = np.nonzero(self.ene_vs_disp[:,0] == 0.0)[0][0]
        minene = self.ene_vs_disp[minidx][1]
        minposforces = self.positions[minidx*self.atoms+self.atoms:minidx*self.atoms+self.atoms+self.atoms]
        self.positions[:self.atoms] = minposforces

        # udftmev_vs_disp
        #udftmev = (self.ene_vs_disp[:,1] - minene)*1000/(self.atoms)
        udftmev = (self.ene_vs_disp[:,1] - minene)
        self.udftmev_vs_disp = np.array((self.dstring,udftmev)).transpose()

        # dUdL
        self.dudl = np.zeros((np.array(self.folder).shape[0],9))
        self.dudl[:,4] = self.udftmev_vs_disp[:,1]
        self.dudl[:,0] = self.udftmev_vs_disp[:,0]

        np.savetxt(self.folderparam+"/ene_vs_disp",self.ene_vs_disp)
        np.savetxt(self.folderparam+"/POSITIONs",self.positions,fmt="%.8f")
        np.savetxt(self.folderparam+"/udftmev_vs_disp",self.udftmev_vs_disp)
        np.savetxt(self.folderparam+"/dUdL",self.dudl,fmt="%7.3f%10.1f%9.1f%9.1f%14.2f%10.2f%14.2f%10.2f%10.2f",header=" step   time(fs)  temp(K) average       U(meV/at)    Uref          dUdL   average    offset")

        takeevery = 1
        #print "lines:",self.dudl.shape[0]
        if self.dudl.shape[0] > 100:  # 100 - 200
            takeevery = 6
        if self.dudl.shape[0] > 200:
            takeevery = 8
        if self.dudl.shape[0] > 300:
            takeevery = 10
        positions_short = self.positions[self.atoms:].flatten().reshape(self.dudl.shape[0],self.atoms*6)[::takeevery]
        nr1 = positions_short.shape[0]
        nr2 = positions_short.shape[1]
        tmppositions_short = positions_short.flatten().reshape(nr1*nr2/6,6)
        #print "kkk:",tmppositions_short.shape
        self.positions_short = np.zeros((tmppositions_short.shape[0]+self.atoms,6))
        self.positions_short[:self.atoms] = self.positions[:self.atoms]
        self.positions_short[self.atoms:] = tmppositions_short

        self.dudl_short = self.dudl[::takeevery]
        self.ene_vs_disp_short = self.ene_vs_disp[::takeevery]
        self.udftmev_vs_disp_short = self.udftmev_vs_disp[::takeevery]
        print "self.udftmev_vs_disp_short:"
        print self.udftmev_vs_disp_short

        np.savetxt(self.folderparam+"/ene_vs_disp_short",self.ene_vs_disp_short)
        np.savetxt(self.folderparam+"/POSITIONs_short",self.positions_short,fmt="%.8f")
        np.savetxt(self.folderparam+"/udftmev_vs_disp_short",self.udftmev_vs_disp_short)
        np.savetxt(self.folderparam+"/dUdL_short",self.dudl_short,fmt="%7.3f%10.1f%9.1f%9.1f%14.2f%10.2f%14.2f%10.2f%10.2f",header=" step   time(fs)  temp(K) average       U(meV/at)    Uref          dUdL   average    offset")
        self.repeat_supercell_positions_forces()
        self.saveData(self.folderparam,"nnforces")
        return

    def repeat_supercell_positions_forces(self):
        ''' before the whole cell can be repeated, the mapping to the original cell
        has to be remove (in other words: the atoms have to be close to their origial
        undisplaced positions, otherwise crystal_generator.center_atoms_around_atom
        does not work properly --> Therefore we use self.pr rather then self.p'''
        nsc = 2  # for now on we will just double the supercell, in general maybe factor
        self.P = np.zeros((self.pr.shape[0],self.pr.shape[1]*nsc**3,3))
        self.F = np.zeros((self.pr.shape[0],self.pr.shape[1]*nsc**3,3))
        self.Pr = np.zeros((self.pr.shape[0],self.pr.shape[1]*nsc**3,3)) #  (centered around 0 atom)
        for idx,i in enumerate(self.dvec):
            # coords
            coord_cart = np.copy(self.pr[idx])
            self.crystal = crystal_generator.crystal()
            self.crystal.load_positions_cell(coord_cart = coord_cart, cell = self.sc)


            # coords big (repeat supercell)
            #print "coords:"
            self.crystalbig = crystal_generator.supercell()
            self.crystalbig.create_supercell(  self.crystal, nsc, nsc, nsc, newsorting = True )
            self.P[idx] = self.crystalbig.rcar
            self.crystalbigcentered = copy.deepcopy(self.crystalbig)
            self.crystalbigcentered.center_atoms_around_atom(0,coord_cart=self.crystalbigcentered.rcar,cell=self.crystalbigcentered.cellvec)
            self.Pr[idx] = self.crystalbigcentered.rcar


            # forces
            #print "forces:"
            forces = np.copy(self.f[idx])
            self.crystalforces = crystal_generator.crystal()
            self.crystalforces.load_positions_cell(coord_cart = forces , cell = self.sc)

            # forces big (repeat supercell)
            self.crystalforcesbig = crystal_generator.supercell()
            self.crystalforcesbig.create_supercell(  self.crystalforces, nsc, nsc, nsc, newsorting = True, shiftforces = True )
            self.F[idx] = self.crystalforcesbig.rcar
            #posvec = self.crystalposbig.rcar[atpos]
            #negvec = self.crystalnegbig.rcar[atneg]
            #if 'xdir' in self.folderparam.split("/")[-1]:
            #    self.disp = 'xdir'
            #if 'quer' in self.folderparam.split("/")[-1]:
            #    self.disp = 'quer'
            self.disp = self.folderparam.split("/")[-1].split("_")[1]

        #################################################################################
        # also repeat eqcoords to get crystal0
        #################################################################################
        self.crystal0 = crystal_generator.crystal()
        self.crystal0.load_positions_cell(coord_rrel = self.eqcoords, cell = self.sc)
        self.crystal0big = crystal_generator.supercell()
        self.crystal0big.create_supercell(  self.crystal0, nsc, nsc, nsc, newsorting = True )
        self.P0 = self.crystal0big.rcar
        self.crystal0bigcentered = copy.deepcopy(self.crystal0big)
        self.crystal0bigcentered.center_atoms_around_atom(0,coord_cart=self.crystal0bigcentered.rcar,cell=self.crystal0bigcentered.cellvec)
        self.Pr0 = self.crystal0bigcentered.rcar
        self.F0 = np.zeros((self.Pr0.shape[0],self.Pr0.shape[1]))
        return

    def saveData( self, outdir, outfname ):
        '''
        '''
        print utils.printred("      saving npz : "+str(outdir+'/'+outfname+".npz"))
        self.npz = outdir + "/" + outfname + ".npz"
        np.savez_compressed( outdir + "/" + outfname,
        npz             = self.npz,
        element         = self.element,
        folderparam     = self.folderparam,
        nndist          = self.nndist,
        f               = self.f,
        p               = self.p,
        pr              = self.pr,
        F               = self.F,
        P               = self.P,
        Pr              = self.Pr,
        P0              = self.P0,
        F0              = self.F0,
        Pr0             = self.Pr0,
        dstring         = self.dstring,
        dvec            = self.dvec,
        dnorm           = self.dnorm,
        atoms           = self.atoms,
        eqcoords        = self.eqcoords,
        folder          = self.folder,
        sc              = self.sc,
        disp            = self.disp,
        a               = self.a,
        positions       = self.positions,
        ene_vs_disp     = self.ene_vs_disp,
        dudl            = self.dudl,
        udftmev_vs_disp = self.udftmev_vs_disp,
        struct          = self.struct,
        e               = self.e
        )
        # class instances (as self.atom) seem to be a problem for np.load, therefore are
        # not saved but class has to be reloaded every time (this is ok for small and
        # quick loads (like my_atom) but not in general
        self.parametrize_lon()
        return

    def loadData( self, indir, infname ):
        '''
        '''
        print utils.printgreen("     loading npz : "+str(indir+'/'+infname+".npz"))
        var = np.load( indir + "/" + infname + ".npz" )
        self.npz                = str(var['npz'])
        self.element            = str(var['element'])               # this is true for xdisp and quer
        self.folderparam        = str(var['folderparam'])
        self.nndist             = var['nndist']
        self.f                  = var['f']
        self.p                  = var['p']
        self.pr                 = var['pr']
        self.F                  = var['F']
        self.P                  = var['P']
        self.Pr                 = var['Pr']
        self.P0                 = var['P0']
        self.F0                 = var['F0']
        self.Pr0                = var['Pr0']
        self.dstring            = var['dstring']
        self.dvec               = var['dvec']
        self.dnorm              = var['dnorm']
        self.atoms              = var['atoms']                      # this is true for xdisp and quer
        self.eqcoords           = var['eqcoords']                   # this is true for xdisp and quer
        self.folder             = var['folder']
        self.sc                 = var['sc']                         # this is true for xdisp and quer
        self.disp               = str(var['disp'])
        self.struct             = str(var['struct'])
        self.a                  = var['a']
        self.positions          = var['positions']
        self.ene_vs_disp        = var['ene_vs_disp']
        self.dudl               = var['dudl']
        self.udftmev_vs_disp    = var['udftmev_vs_disp']
        self.e                  = var['e']

        # class instances seem to be a problem for np.load, therefore:
        self.atom = my_atom.atom([self.element])
        #print "type:",type(self.fml)
        #if type(self.fml) != bool:
            #print "self.fml.shape",self.fml.shape
            #print self.fml
            #print self.fml[3,3]
        if os.path.isdir(self.folderparam) != True:
            sys.exit(self.folderparam+" does not exist, rm nnforces.npz and rerun parametrization")
        self.parametrize_lon()
        return

    def saveDatapd(self, filename):
        ''' saves pandas dataframe '''
        self.pkl = filename
        print utils.printred("      saving pkl : "+str(filename))
        with open(filename, 'w') as f:
            #pickle.dump(self.list_save_load, f)
            pickle.dump([self.data,self.fit], f)
        return

    def loadDatapd(self, filename):
        ''' loads pandas dataframe '''
        if type(filename) != str:
            print "filename:",filename
            sys.exit("loadDatapd  needs as input a picklefile (str) but got "+str(filename)+" of type "+str(type(filename)))
        if os.path.isfile(filename) != True:
            sys.exit("file "+filename+" does not exist!")
        #self.data = pd.read_pickle(filename)
        print utils.printgreen("     loading pkl : "+str(filename))
        with open(filename) as f:
            self.data, self.fit = pickle.load(f)
            #self.data, self.eqcoords, self.element = pickle.load(f)
            #*self.list_save_load() = pickle.load(f)
            #for i in enumerate(pickle.load(f)):
            #    self.list_save_load[i] = i
        return

    def saveFML( self, outdir, outfname ):
        print utils.printpink("      saving npz : "+str(outdir+'/'+outfname+".npz"))
        self.parametrize_fml_npz_path = outdir + "/" + outfname + ".npz"
        np.savez_compressed( outdir + "/" + outfname,
        parametrize_fml_npz_path           = self.parametrize_fml_npz_path,
        fml                  = self.fml,
        FML                  = self.FML,
        fml_u1nn_pottype     = self.fml_u1nn_pottype,
        fml_u1nn_potparam    = self.fml_u1nn_potparam,
        fml_u1nn_potadd      = self.fml_u1nn_potadd,
        fml_u1nn_potaddparam = self.fml_u1nn_potaddparam,
        angles               = self.angles,
        fprojlongold         = self.fprojlongold
        )
        return

    def loadFML( self, indir, infname ):
        '''
        '''
        print utils.printgreen("     loading npz : "+str(indir+'/'+infname+".npz"))
        var = np.load( indir + "/" + infname + ".npz" )
        self.parametrize_fml_npz_path            = str(var['parametrize_fml_npz_path'])
        self.fml                   = var['fml']
        self.FML                   = var['FML']
        self.fml_u1nn_pottype      = str(var['fml_u1nn_pottype'])
        self.fml_u1nn_potparam     = var['fml_u1nn_potparam']
        self.fml_u1nn_potadd       = str(var['fml_u1nn_potadd'])
        self.fml_u1nn_potaddparam  = var['fml_u1nn_potaddparam']
        print "--->",self.FML.shape
        return

    def infodos(self, verbose = False):
        ''' get DOS, DOScut for corresponding element;
            needs to be called after loadforces or loadData in order to have the self.element information'''
        if verbose:
            print utils.printgreen("     infodos ...")

        # get DOS / DOScut
        self.jobvorlage_all = '/Users/glensk/Dropbox/Understand_distributions/jobvorlage_all/'
        self.jobvorlage_search = self.jobvorlage_all+'2x2x2sc_'+self.atom.symbol[0]+"_*"
        self.jobvorlage = glob.glob(self.jobvorlage_search)
        if len(self.jobvorlage) != 1:
            print "self.jobvorlage_search:",self.jobvorlage_search
            print "self.jobvorlage:",self.jobvorlage
            sys.exit("self.jobvorlage not found (or not just one path)")
        else:
            self.jobvorlage = self.jobvorlage[0]

        ##################################################################################
        ##################################################################################
        ##################################################################################
        # load data from jobvorlage (at least the path)
        ##################################################################################
        ##################################################################################
        ##################################################################################
        self.jobvorlage_avg_dudl_low_fre_file = self.jobvorlage+"/avg_dUdL_low_fre"
        self.jobvorlage_avg_dudl_low_fre = np.loadtxt(self.jobvorlage_avg_dudl_low_fre_file)


        ######################################################################
        # tveclonall.dat
        self.jobvorlage_vecnormlon_file = self.jobvorlage+"/tveclonall.dat"   # somewhat shorter
        self.jobvorlage_vecnormlon_file = self.jobvorlage+"/atoms_1nn_all_f/dfn_1.0"
        if os.path.isfile(self.jobvorlage_vecnormlon_file) != True:
            sys.exit(self.jobvorlage_vecnormlon_file + " self.jobvorlage_vecnormlon_file does not exist 3 !")
        if verbose:
            print "     self.jobvorlage_vecnormlon_file   :",self.jobvorlage_vecnormlon_file

        ######################################################################
        # DOSlon  DOSDOSDOS
        ######################################################################
        shells = 5
        self.DOSlon_fileorig    = [ False for shell in range(shells+1) ]
        self.DOSlon_file        = [ False for shell in range(shells+1) ]
        self.DOSlon_filecut     = [ False for shell in range(shells+1) ]
        self.DOSlon             = [ False for shell in range(shells+1) ]
        self.DOSloncut          = [ False for shell in range(shells+1) ]
        for shell in range(shells): self.DOSlon_fileorig[shell+1] = self.jobvorlage+"/longvecnorm_"+str(shell+1)
        for shell in range(shells): self.DOSlon_file[shell+1] = self.jobvorlage+"/longvecnorm_"+str(shell+1)+"_DOS"
        for shell in range(shells): self.DOSlon_filecut[shell+1] = self.jobvorlage+"/longvecnorm_"+str(shell+1)+"_DOScut"
        for shell in range(shells):
            if os.path.isfile(self.DOSlon_file[shell+1]) == True:
                self.DOSlon[shell+1] = np.loadtxt(self.DOSlon_file[shell+1])

        # create self.DOSlon_filecut in case it does not exist
        for shell in range(shells):
            if os.path.isfile(self.DOSlon_filecut[shell+1]) != True:
                if os.path.isfile(self.DOSlon_fileorig[shell+1]) == True:
                    print "creating DOScut"
                    dos,doscut = utils.DOS_cut(
                        longvecnormfile = self.DOSlon_fileorig[shell+1],
                        dos = self.DOSlon[shell+1],
                        nndist = self.nndist[shell+1],
                        filenamedos = self.DOSlon_file[shell+1])
        for shell in range(shells):
            if os.path.isfile(self.DOSlon_filecut[shell+1]) == True:
                self.DOSloncut[shell+1] = np.loadtxt(self.DOSlon_filecut[shell+1])

        #sys.exit()
        ##self.jobvorlage_DOSlon_file = self.jobvorlage+"/DOS_dfn_1.0_py3.0"  # gebraucht in der joberstellung
        #self.jobvorlage_DOSlon_filecut = self.jobvorlage+"/DOS_dfn_1.0_py3.0cut"  # gebraucht in der joberstellung
        #self.jobvorlage_DOSlon_filecutshiftednndist = self.jobvorlage+"/DOS_dfn_1.0_py3.0cutshiftednndist"

        #self.jobvorlage_DOSlon = np.loadtxt(self.jobvorlage_DOSlon_file)   # gebraucht fuer die parametrisierung
        #self.jobvorlage_DOSloncut = np.loadtxt(self.jobvorlage_DOSlon_filecut)  # gebraucht fuer die parametrisierung
        #self.jobvorlage_DOSlon_max = self.jobvorlage_DOSlon[:,0].max()
        #self.jobvorlage_DOSlon_min = self.jobvorlage_DOSlon[:,0].min()
        #self.jobvorlage_DOSloncut_max = self.jobvorlage_DOSloncut[:,0].max()
        #self.jobvorlage_DOSloncut_min = self.jobvorlage_DOSloncut[:,0].min()
        ##if type(nndist) != False:
        #if False:
        #    # currently we just do this for the 1NN self.nndist[1]
        #    self.jobvorlage_DOSloncutshiftednndist = np.array([self.jobvorlage_DOSloncut[:,0]-self.nndist[1],self.jobvorlage_DOSloncut[:,1]]).transpose()
        #    np.savetxt(self.jobvorlage_DOSlon_filecutshiftednndist,self.jobvorlage_DOSloncutshiftednndist)
        if verbose:
            print "     self.DOSlon: OK"
            print "     self.DOSloncut: OK"
            #print "     self.DOSlon_file   :",self.jobvorlage_DOSlon_file,self.jobvorlage_DOSlon_min,self.jobvorlage_DOSlon_max
            #print "     self.DOSlon_filecut:",self.jobvorlage_DOSlon_filecut,self.jobvorlage_DOSloncut_min,self.jobvorlage_DOSloncut_max
        return

    def parametrize_whatatom(self, shell = False, contrib = False):
        ''' defines which parametrization we are interested in '''
        if shell == False:
            sys.exit("please define the shell for which you want to get the parametrization")
        if contrib not in [ 'lon', 'to', 'ti' ]:
            sys.exit("please specify a contribution for the parametrization lon/to/ti")
        #                   V structure
        #                        V shell
        #                             V   V ... atoms to check
        listcheck = [
        [ 'fcc', 1, 'lon',32, 24, 126 ],   # [ 1.995,  1.995,  0.   ], [-1.995, -1.995,   0.   ]
        [ 'fcc', 2, 'lon',32,  4,  36 ],   # [ 3.99,  0.  ,  0.     ], [-3.99,  0.  ,   0.     ]
        [ 'fcc', 3, 'lon',32, 12, 239 ],   # [ 3.99 ,  1.995,  1.995], [-3.99 ,  -1.995, -1.995]
        [ 'fcc', 4, 'lon',32,  6, 102 ],   # [ 3.99,  3.99,  0.  ]  , [-3.99,  -3.99,  0.  ]

        [ 'fcc', 1, 'tox',32, 8, 203 ],   # [ 0.   ,  2.065,  2.065],[ 0.   , -2.065, -2.065]
        [ 'fcc', 2, 'tox',32, 2,  66 ],   # [ 0.  ,  4.13,  0.00], [0,  -4.13  ,   0.     ]
        [ 'fcc', 3, 'tox',32, 0,  0  ],   # [ 3.99 ,  1.995,  1.995], [-3.99 ,  -1.995, -1.995]
        [ 'fcc', 4, 'tox',32, 0,  0  ],   # [ 3.99,  3.99,  0.  ]  , [-3.99,  -3.99,  0.  ]

        [ 'fcc', 1, 'tix',32, 0,  0  ],   # [ 0.   ,  2.065,  2.065],[ 0.   , -2.065, -2.065]
        [ 'fcc', 2, 'tix',32, 0,  0  ],   # [ 0.  ,  4.13,  0.00], [0,  -4.13  ,   0.     ]
        [ 'fcc', 3, 'tix',32, 0,  0  ],   # [ 3.99 ,  1.995,  1.995], [-3.99 ,  -1.995, -1.995]
        [ 'fcc', 4, 'tix',32, 0,  0  ],   # [ 3.99,  3.99,  0.  ]  , [-3.99,  -3.99,  0.  ]

        [ 'fcc', 1, 'lon',108, 81, 429 ],   # [ 1.995,  1.995,  0.   ], [-1.995, -1.995,  0.   ]
        [ 'fcc', 2, 'lon',108,  9, 126 ],   # [ 3.99,  0.  ,  0.     ], [-3.99,  0.  ,  0.     ]
        [ 'fcc', 3, 'lon',108, 36, 809 ],   # [ 3.99 ,  1.995,  1.995], [-3.99 ,  -1.995,-1.995]
        [ 'fcc', 4, 'lon',108, 12, 348 ],   # [ 3.99,  3.99,  0.  ]   , [-3.99,  -3.99,  0.  ]

        [ 'fcc', 1, 'lon',108, 0,  0 ],   # [ 1.995,  1.995,  0.   ], [-1.995, -1.995,  0.   ]
        [ 'fcc', 2, 'lon',108, 0,  0 ],   # [ 3.99,  0.  ,  0.     ], [-3.99,  0.  ,  0.     ]
        [ 'fcc', 3, 'lon',108, 0,  0 ],   # [ 3.99 ,  1.995,  1.995], [-3.99 ,  -1.995,-1.995]
        [ 'fcc', 4, 'lon',108, 0,  0 ],   # [ 3.99,  3.99,  0.  ]   , [-3.99,  -3.99,  0.  ]

        ]

        # for 2x2x2sc:
        # 1     12     2.82136 [ 8  9 10 11 16 17 20 21 24 26 28 30]
        # 2     3         3.99 [1 2 4]
        # 3     12     4.88673 [12 13 14 15 18 19 22 23 25 27 29 31]
        # 4     3      5.64271 [3 5 6]
        # 5     1      6.91088 [7]
        #
        # NN     atoms       distance atomindex
        # ----------------------------------------------------------------------------
        # 1     12     2.82136 [  8  16  24  52  60  74  90 126 137 145 181 203]
        # 2     6         3.99 [  1   2   4  36  66 129]
        # 3     24     4.88673 [ 12  18  25  44  54  61  78  82  91 110 118 127 141 147 153]  ...
        # 4     12     5.64271 [  3   5   6  37  38  67  70 102 131 133 165 195]
        # 5     24     6.30874 [  9  10  17  20  26  28  48  53  56  62  72  75  88  94 122]  ...
        # ----------------------------------------------------------------------------
        # array([[ 0.   ,  1.995,  1.995],  8
        #        [ 1.995,  0.   ,  1.995],  16
        #        [ 1.995,  1.995,  0.   ],  24
        #        [-1.995,  0.   ,  1.995],  52
        #        [-1.995,  1.995,  0.   ],  60
        #        [ 0.   , -1.995,  1.995],  74
        #        [ 1.995, -1.995,  0.   ],  90
        #        [-1.995, -1.995,  0.   ],  126
        #        [ 0.   ,  1.995, -1.995],  137
        #        [ 1.995,  0.   , -1.995],  145
        #        [-1.995,  0.   , -1.995],  181
        #        [ 0.   , -1.995, -1.995]]) 203


        self.param_atoms = False
        for i in listcheck:
            if self.struct == i[0]:
                if shell == i[1]:
                    if contrib == i[2]:
                        if self.atoms == i[3]:
                            self.param_atoms = np.array(i[4:])
        if type(self.param_atoms) == bool:
            sys.exit("self.param_atoms for parametrize_lon_whatatom not found!")
        return self.param_atoms

    def plot_lon(self, shell = 1, out = 'd', nice = False):
        ''' plot the lon datapoints for a certain shell
            out:{\'a\' == all,\'d\' == dos, \'dc\' == doscut '''
        import matplotlib.pyplot as plt
        if nice == True:
            from matplotlib import rc
            plt.rc('text', usetex=True)
            plt.rc('font', family='serif')
        plt.clf()
        if out == 'a':
            plotfunc = self.lonall
        if out == 'd':
            plotfunc = self.lonallDOS
        plt.plot(plotfunc[shell][:,0],plotfunc[shell][:,1],'b.')
        plt.plot(np.array([self.nndist[shell]]),np.array([0.0]),'ro')

        xliml,xlimr,ylimmin,ylimmax = utils.plot_find_ylim_from_xlim(
            plotfunc[shell],
            xliml = self.nndist[shell]-1.0,
            xlimr = self.nndist[shell]+1.0)

        plt.xlim(xliml,xlimr)
        plt.ylim(ylimmin,ylimmax)
        plt.xlabel('Internuclear distance (\AA)', fontsize=14)
        plt.ylabel('Internuclear force (eV/\AA)')
        plt.title(str(shell)+' nn')
        plt.grid(True, which='both')
        return

    def parametrize_lon(self):
        ''' creates parametrization for long vectors
            make sure loadforces() was called before '''
        # Example: Ir, 2x2x2sc 1nn 0.5 quer displacement
        # Folder: /Users/glensk/Dropbox/understand_distributions/displacements_dense/Ir/2x2x2sc_quer_3x3x3kp/3.99Ang_0.5
        # positions VASP: 0.35355 0.35355 0.00000 = OUTCAR_positions-last-ARRAY.sh | head -1
        # np.linalg.norm([0.35355 ,0.35355, 0.00000]) = 0.5
        # OUTCAR_positions-last-ARRAY.sh | head -25 | tail -1 = 3.220865 3.220865 0.000000
        # np.linalg.norm([3.220865, 3.220865, 0.000000]) = 4.5549909655728191
        #
        # self.u1nn_pottype          = 'morse'
        # self.u1nn_potparam = '0.15282579184613_2.23016908333589_2.8213560574986'
        # lon: 313.66 313.66
        # VASP: 322.1 322.1
        # diff: -8.4 -8.4 0.0
        # pot.param.data[1].lon.quer.dos.morse
        # array([ 0.25044761787896,  1.97031460132089,  2.82135605693432])
        #
        # pot.param.data[1].lon.quer.dos.points =
        # [ 2.32136085195732, -4.55499096557282]
        #
        # self.u1nn_potaddparam = [-1844.659573665, 4050.2563377033, -3616.049646782, 1631.8553962654, -362.0183968268, 30.058737768156, -5.9095020203, 5.3612411576369, -1.669404234424, 0.2239125568498, -0.011411984973]
        # diff: -0.2 -0.2 0.0
        #
        print utils.printblue("     pot_parametrize.parametrize_lon from T=0K displacements ...")

        # if pkl file does not exist
        if type(self.pkl) != bool:
            if os.path.isfile(self.pkl) != True:
                print "os.path.isfile(self.pkl):", os.path.isfile(self.pkl)
                print "os.path.isfile(self.pkl) is not true .... going to create it"
                self.parametrize_data_pkl = True

        # in case there is no pkl we have to parametrize pkl
        #if os.path.isfile(self.folderparam+'/data.pkl') != True:
        #    self.parametrize_data_pkl = True

        # in case we parametrize npz we also want to parametrize pkl
        if self.parametrize_npz == True:
            self.parametrize_data_pkl = True

        # if we dont have to parametrize pkl and the pkl file exists
        if self.verbose == True:
            print "self.parametrize_data_pkl:",self.parametrize_data_pkl
        if self.parametrize_data_pkl == False:
            #if os.path.isfile(self.folderparam+"/data.pkl") == True:
            #print utils.printred("      loading "+self.folderparam+'/data.pkl')
            #self.loadDatapd(self.folderparam+'/data.pkl')
            #print "self.pkllllllllllllll:",self.pkl
            self.loadDatapd(self.pkl)
            #self.loadFitpd(self.folderparam+'/fit.pkl')
            # xdir has usually all the data whil quer has only quer data
            return

        # in case pkl file exists, we want to load it and add to/replace items of this file

        #########################################################################
        # create self.func_lon  (wrong for 2nd shell 2x2x2sc xdir)
        #########################################################################
        self.func_lon = np.zeros((self.Pr.shape[1],self.Pr.shape[0],2))
        if self.verbose:
            print "self.Pr.shape:",self.Pr.shape
        for i in np.arange(self.Pr.shape[1]):   # 0 ... 255
            #sign_whole = -np.sum(self.P[:,0],axis=1)[i]
            self.xvalues = np.linalg.norm(self.Pr[:,i],axis=1)
            self.yvalues = np.linalg.norm(self.F[:,i],axis=1)
            self.ysign = np.sum(-self.F[:,i],axis=1)
            if self.element == 'Pb':
                if self.disp == '3nnd':
                    if i == 36 or i == 809:
                        # for the Pb 3NN the atom at [3.99, 1.99, 1.99] moves in -x direction towards the the displaced atom (displacement e.g. [0.1, 0.05, 0.05]) contrary to the usual behaviour (which would be the positive x direction) , while the y and z direction to the usual thing and behave repulsive.
                        pass
                        #self.ysign = np.ones(self.ysign.shape[0]) * -1
            # print pot.param.xvalues
            # print pot.param.yvalues
            self.yvalues_and_sign = np.copysign(self.yvalues,self.ysign)
            self.func_lon[i] = np.array([self.xvalues,self.yvalues_and_sign]).transpose()

            # for the 2x2x2sc and the 2nd shell we dont have to do anything with the sign!

            #if self.disp == 'xdir':
            #    #self.ysign = np.sum(-self.F[:,i],axis=1) # does not work for Ir xdir
            #    self.ysign = -self.F[:,i][:,1] # if we displace in xdir and atom wants to go gown towards the xaxis, this is probably attractive
            #if self.disp == 'quer':
            #    #self.ysign = -self.F[:,i][:,0]
            #    self.ysign = np.sum(-self.F[:,i],axis=1) # does not work for Ir xdir .... mayby it does work better in the end?
            #self.func_lon[i] = np.array([self.xvalues,self.yvalues]).transpose()
            # generell:
            #       - t{o,i}{x,y,z} wird in komponenten aufgesplittet, also kein problem mit vorzeichen
            #       - lon: hier koennten die x,y,z komponenten moeglicherweise gemischt sein
            #       - lon: positive Kraft kann etweder attraktiv oder repulsiv sein, je nachdem wohin wir auslenken
            #       - lon: am einfachsten wird es sein wenn man das vorzeichen (bestimmt attraktiv, repulsiv)
            #               einfach an die auslenkung koppelt. Bei positiver auslenkun bleibt alles, bei negativer
            #               auslenkung wird das vorzeichen umgedreht.

            ##################################################################
            if np.sum(self.Pr[:,i]) < 0.0:
                self.func_lon[i][:,1] = -1.*self.func_lon[i][:,1]



        #########################################################################
        # make corrections if jump in forces which appear obviously wrong (only if forces close to 0)
        #########################################################################
        # 1 NN quer 2x2x2 all behave well
        # 1 NN quer 3x3x3 all behave well
        # 2 NN xdir 2x2x2 all behave well
        # 2 NN xdir 3x3x3 all behave well
        # 3 NN xdir 3x3x3 Al, Rh, Ir well behaving functions but partly wrong sign
        # 3 NN xdir 2x2x2 all BIG JUMPS, unnecessary
        #
        # 2 NN behave well in 2x2x2sc and 3x3x3sc, for Cu in quer disp we have to ensure that the forces does not jump close to 0
        if self.verbose:
            print "self._shells:",self._shells
        for self.shell in self._shells:
            self.param_atoms = self.parametrize_whatatom(shell = self.shell, contrib = 'lon')
            if self.shell == 2 and self.element == 'Cu' and self.disp == 'quer':
                print "NN:",self.nndist[2]
                print np.linalg.norm(self.Pr[:,i],axis=1)
                changesign = np.nonzero(np.linalg.norm(self.Pr[:,i],axis=1) < self.nndist[2])
                print changesign
                pass


        #########################################################################
        # sort forces only here, since from here on those are not sorted as self.P self.Pr self.F anymore
        #########################################################################
        for i in np.arange(self.Pr.shape[1]):   # 0 ... 255
            self.func_lon[i] = self.func_lon[i][self.func_lon[i][:,0].argsort()]

            #########################################################################
            # forces gained from OUTCAR's seems to be ok
            #########################################################################
            #print "np.linalg.norm(pot.parametrize_xdir.Pr[:,36][50])",np.linalg.norm(self.Pr[:,36][50])
            #print "np.linalg.norm(pot.parametrize_xdir.F[:,36][50])",np.linalg.norm(self.F[:,36][50])
            #print "jup, (4/36) sind die atome in der 2x2x2sc, (9/126) sind die atome in der 3x3x3sc, checke beides:"

            #print "3x3x3sc: (repulsive)"
            #print "pot.parametrize_xdir.Pr[:,9][50]: == array([ 3.49,  0.  ,  0.  ])",pot.parametrize_xdir.Pr[:,9][50]
            #print "pot.parametrize_xdir.P[:,9][50]: == array([ 3.99,  0.  ,  0.  ])",pot.parametrize_xdir.P[:,9][50]
            #print "pot.parametrize_xdir.F[:,9][50]: == array([ 0.480286,  0.      ,  0.      ])",pot.parametrize_xdir.F[:,9][50]
            #print "OUTCAR_positions-last-ARRAY.sh  3.99Ang_0.5/ | head -10 | tail -1 --> 3.99000 0.00000 0.00000"
            #print "OUTCAR_forces-last-ARRAY.sh  3.99Ang_0.5/ | head -10 | tail -1 --> 0.480286 0.000000 0.000000"

            #print "pot.parametrize_xdir.Pr[:,126][50]: == array([-4.49,  0.  ,  0.  ])",pot.parametrize_xdir.Pr[:,126][50]
            #print "pot.parametrize_xdir.P[:,126][50]: == array([ 19.95,   0.  ,   0.  ])",pot.parametrize_xdir.P[:,126][50]
            #print "pot.parametrize_xdir.F[:,126][50]: == array([ 0.18607,  0.     ,  0.     ])",pot.parametrize_xdir.F[:,126][50]
            #print "OUTCAR_positions-last-ARRAY.sh 3.99Ang_0.5/ | head -19 | tail -1 --> 7.98000 0.00000 0.00000"
            #print "OUTCAR_forces-last-ARRAY.sh  3.99Ang_0.5/ | head -19 | tail -1 --> 0.186070 0.000000 0.000000"


        self.fitfolder = self.folderparam+'/nnforces/'
        if os.path.isdir(self.fitfolder) != True:
            os.makedirs(self.fitfolder)
        if os.path.isdir(self.fitfolder+'/datashifted/') != True:
            os.makedirs(self.fitfolder+'/datashifted/')
        if os.path.isdir(self.fitfolder+'/diffshifted/') != True:
            os.makedirs(self.fitfolder+'/diffshifted/')
        if os.path.isdir(self.fitfolder+'/diffshiftedpolydiff/') != True:
            os.makedirs(self.fitfolder+'/diffshiftedpolydiff/')


        ##############################################################################
        # define lon_{1,2,3,4}nn_{alle,dos,doscut}
        # - points
        # - dos
        ##############################################################################
        self.infodos()
        for self.shell in self._shells:
            shellstr = str(self.shell)
            self.param_atoms = self.parametrize_whatatom(shell = self.shell, contrib = 'lon')
            print "     shell:",self.shell,"self.param_atoms:",self.param_atoms,"nndist[shell]:",self.nndist[self.shell],"self.disp:",self.disp


            # concat pos / negside and append 0.0 point and sort
            # pot.param.parametrize_whatatom(2,'lon') array([ 24, 126])
            # pot.param.parametrize_whatatom(2,'lon') array([ 4, 36])
            # pot.param.parametrize_whatatom(3,'lon') array([ 12, 239])
            # pot.param.parametrize_whatatom(4,'lon') array([ 6, 102])
            lonall = np.concatenate((self.func_lon[self.param_atoms[0]],self.func_lon[self.param_atoms[1]]))
            lonall = np.concatenate((lonall,np.array([[self.nndist[self.shell],0.0]])))

            # remove duplicates
            lonall = utils.remove_duplicates_of_2d_array_within_tolerance(lonall,1e-6, 1e-4)
            lonall = lonall[lonall[:,0].argsort()]
            self.data[self.shell,'lon',self.disp,'alle'].ix['points'] = lonall
            self.data[self.shell,'lon',self.disp,'alle'].ix['pointsshifted'] = np.transpose([lonall[:,0]-self.nndist[self.shell],lonall[:,1]])

            # at this place we have to make sure we loaded DOSlon stuff (DOSinfo)
            if type(self.DOSlon[self.shell]) != bool:
                lonallDOS      = utils.cut_function_at_DOS(lonall,self.DOSlon[self.shell])
                lonallDOScut   = utils.cut_function_at_DOS(lonall,self.DOSloncut[self.shell])

                self.data[self.shell,'lon',self.disp,'dos'].ix['points'] = lonallDOS
                self.data[self.shell,'lon',self.disp,'doscut'].ix['points'] = lonallDOScut
                self.data[self.shell,'lon',self.disp,'dos'].ix['pointsshifted'] = np.transpose([lonallDOS[:,0]-self.nndist[self.shell],lonallDOS[:,1]])
                self.data[self.shell,'lon',self.disp,'doscut'].ix['pointsshifted'] = np.transpose([lonallDOScut[:,0]-self.nndist[self.shell],lonallDOScut[:,1]])
            else:
                self.data[self.shell,'lon',self.disp,'dos'].ix['points'] = False
                self.data[self.shell,'lon',self.disp,'doscut'].ix['points'] = False
                self.data[self.shell,'lon',self.disp,'dos'].ix['pointsshifted'] = False
                self.data[self.shell,'lon',self.disp,'doscut'].ix['pointsshifted'] = False
                print utils.printred("DOS not cut for shell "+str(self.shell))



            self.data[self.shell,'lon',self.disp,'alle'].ix['dos']    = self.DOSlon[self.shell]
            self.data[self.shell,'lon',self.disp,'dos'].ix['dos']    = self.DOSlon[self.shell]
            self.data[self.shell,'lon',self.disp,'doscut'].ix['dos'] = self.DOSloncut[self.shell]

            if type(self.DOSlon[self.shell]) != bool:
                self.data[self.shell,'lon',self.disp,'alle'].ix['pot']    =    dos_to_pot(self.DOSlon[self.shell],self.atom.melting_rounded[0])
                self.data[self.shell,'lon',self.disp,'dos'].ix['pot']    =    dos_to_pot(self.DOSlon[self.shell],self.atom.melting_rounded[0])
                self.data[self.shell,'lon',self.disp,'alle'].ix['dosshifted']    = np.transpose([self.DOSlon[self.shell][:,0]-self.nndist[self.shell],self.DOSlon[self.shell][:,1]])
                self.data[self.shell,'lon',self.disp,'dos'].ix['dosshifted']    = np.transpose([self.DOSlon[self.shell][:,0]-self.nndist[self.shell],self.DOSlon[self.shell][:,1]])
                self.data[self.shell,'lon',self.disp,'alle'].ix['potshifted']    =    dos_to_pot(self.data[self.shell,'lon',self.disp,'alle'].ix['dosshifted'],self.atom.melting_rounded[0])
                self.data[self.shell,'lon',self.disp,'dos'].ix['potshifted']    =    dos_to_pot(self.data[self.shell,'lon',self.disp,'dos'].ix['dosshifted'],self.atom.melting_rounded[0])
            if type(self.DOSloncut[self.shell]) != bool:
                self.data[self.shell,'lon',self.disp,'doscut'].ix['pot'] = dos_to_pot(self.DOSloncut[self.shell],self.atom.melting_rounded[0])
                self.data[self.shell,'lon',self.disp,'doscut'].ix['dosshifted'] = np.transpose([self.DOSloncut[self.shell][:,0]-self.nndist[self.shell],self.DOSloncut[self.shell][:,1]])
                self.data[self.shell,'lon',self.disp,'doscut'].ix['potshifted'] =    dos_to_pot(self.data[self.shell,'lon',self.disp,'doscut'].ix['dosshifted'],self.atom.melting_rounded[0])

            ##############################################################################
            # save              lon_{1,2,3,4}nn_{alle,dos,doscut} if not already existing
            # save  datashifted/lon_{1,2,3,4}nn_{alle,dos,doscut} if not already existing
            ##############################################################################
            for rang in self._range: # ['alle', 'dos', 'doscut']
                filename = self.fitfolder+'lon'+"_"+shellstr+"nn_"+rang
                filenameshifted = self.fitfolder+'/datashifted/lon'+"_"+shellstr+"nn_"+rang
                self.data[self.shell,'lon',self.disp,rang].ix['file'] = filename
                if os.path.isfile(filename) != True:
                    data = self.data[self.shell,'lon',self.disp,rang].ix['points']
                    if type(data) != bool:
                        print "saving:",filename
                        #print data[:3],type(data),data.shape
                        np.savetxt(filename,self.data[self.shell,'lon',self.disp,rang].ix['points'],fmt="%.7f")
                if os.path.isfile(filenameshifted) != True:
                    data = self.data[self.shell,'lon',self.disp,rang].ix['pointsshifted']
                    if type(data) != bool:
                        print "saving shifted data:",filenameshifted
                        np.savetxt(filenameshifted,self.data[self.shell,'lon',self.disp,rang].ix['pointsshifted'],fmt="%.7f")


        def check_if_exists(foldername, filename):
            ''' checks if corresponding fit exists '''
            searchfor = foldername+"/"+filename+"_fit___*"
            #print "searchfor:",searchfor
            found = glob.glob(searchfor)
            #print "found:",found
            if len(found) == 1:
                return True
            else:
                return False

        ##############################################################################
        # PARAMETRIZATION (of a certain displacement (xdir/quer/3nnd)
        ##############################################################################
        for self.shell in self._shells:
            shellstr = str(self.shell)

            for rang in self._range: # ['alle', 'dos', 'doscut']
                for func in self._func: # ['morse', 'mc1']

                    if func in [ 'morse', 'mc1'] and type(self.data[self.shell,'lon',self.disp,rang].ix['points']) != bool:
                        filenamecheck = 'lon_'+shellstr+"nn_"+rang+"_"+func+"_"

                        # NEW WAY (in some cases it might happen that no fit is found
                        #print "-----> shell:",self.shell,"self.disp:",self.disp,"rang:",rang,"func:",func,"0@t:",self.nndist[self.shell]
                        #print "-----> dat:",self.data[self.shell,'lon',self.disp,rang].ix['points'].shape
                        self.fitting = fit_to_func( fixzeroat = self.nndist[self.shell],function = func,
                                        data=self.data[self.shell,'lon',self.disp,rang].ix['points'],
                                        weights = False)

                        #print "-----> sfd:",self.fitting.deltas[-1]
                        # OLD WAY
                        #params, deltas, deltasshifted, polycoefsrealene,polycoefsrealforces,polydeltas = get_fit_forces_to_pot(
                        #        foldername=self.fitfolder,
                        #        filename='lon'+"_"+shellstr+"nn_"+rang+"_"+func,
                        #        NN=self.nndist[self.shell],pot = func,
                        #        data=self.data[self.shell,'lon',self.disp,rang].ix['points'],
                        #        weights = False)
                        if self.verbose > 0:
                            print "         ",self.shell,rang,func,self.fitting.parameters

                        ##################################################################
                        # write to dataframe
                        ##################################################################
                        self.data[self.shell,'lon',self.disp,rang].ix[func] = self.fitting.parameters
                        self.data[self.shell,'lon',self.disp,rang].ix[func+'fit'] = self.fitting.fit
                        self.fit[ self.shell,'lon',self.disp,rang,func].ix['diff'] = self.fitting.deltas
                        self.fit[ self.shell,'lon',self.disp,rang,func].ix['diffshifted'] = self.fitting.deltasshifted
                        self.fit[ self.shell,'lon',self.disp,rang,func].ix['polybeste'] = self.fitting.polyfit.coefsrealene
                        self.fit[ self.shell,'lon',self.disp,rang,func].ix['polybestf'] = self.fitting.polyfit.coefsrealforces
                        self.fit[ self.shell,'lon',self.disp,rang,func].ix['polybestdiff'] = self.fitting.polyfit.deltas
                        self.fit[ self.shell,'lon',self.disp,rang,func].ix['polybestdiffshifted'] = self.fitting.polyfit.deltasshifted

                        ##################################################################
                        # spline fit
                        ##################################################################
                        if self.verbose > 0:
                            print "         ","spline FIT"
                        if type(self.fitting.deltas) != bool:  # for cases where no fitting was found/done
                            ii    = 1
                            spl   = UnivariateSpline( self.fitting.deltas[:,0][::ii], self.fitting.deltas[:,1][::ii], k=3, s=0.00000000000*self.nndist[self.shell]  )
                            splshifted   = UnivariateSpline( self.fitting.deltasshifted[:,0][::ii], self.fitting.deltasshifted[:,1][::ii], k=3, s=0.00000000000*self.nndist[self.shell]  )
                            self.fit[ self.shell,'lon',self.disp,rang,func].ix['splinedata'] = spl
                            self.fit[ self.shell,'lon',self.disp,rang,func].ix['splineshifteddata'] = splshifted
                            spldiff = np.transpose([self.fitting.deltas[:,0], self.fitting.deltas[:,1] - spl(self.fitting.deltas[:,0])])
                            splshifteddiff = np.transpose([self.fitting.deltasshifted[:,0], self.fitting.deltasshifted[:,1] - splshifted(self.fitting.deltasshifted[:,0])])
                            self.fit[ self.shell,'lon',self.disp,rang,func].ix['splinediff'] = spldiff
                            self.fit[ self.shell,'lon',self.disp,rang,func].ix['splineshifteddiff'] = splshifteddiff

                        ##################################################################
                        # write to hard disk
                        ##################################################################
                        filenameshifteddiff = self.fitfolder+'/diffshifted/lon'+"_"+shellstr+"nn_"+rang+"_"+func
                        filenameshifteddiffpolydiff = self.fitfolder+'/diffshiftedpolydiff/lon'+"_"+shellstr+"nn_"+rang+"_"+func
                        filenameshifteddiffsplinediff = self.fitfolder+'/diffshiftedpolydiff/lon'+"_"+shellstr+"nn_"+rang+"_"+func+"_splinediff"
                        if os.path.isfile(filenameshifteddiff) != True:
                            if type(self.fit[ self.shell,'lon',self.disp,rang,func].ix['polybestdiff']) != bool:
                                print "saving:",filenameshifteddiff
                                np.savetxt(filenameshifteddiff,self.fit[ self.shell,'lon',self.disp,rang,func].ix['diffshifted'],fmt='%.7f')
                                np.savetxt(filenameshifteddiffpolydiff,self.fit[ self.shell,'lon',self.disp,rang,func].ix['polybestdiffshifted'],fmt='%.7f')
                        np.savetxt(filenameshifteddiffsplinediff,splshifteddiff,fmt='%.7f')

                    if type(self.fitting.parameters) == bool:
                        print "         ",self.shell,rang,func,utils.printred(self.fitting.parameters)
                    else:
                        print "         ",self.shell,rang,func,self.fitting.parameters

                            #self.data[self.shell,'lon',self.disp,rang].ix[func] = params
                    #if func in [ 'poly' ]:
                    #    for orderfit in range(1,10):
                    #        filenamecheck = 'lon'+"_poly"+str(orderfit)+"_[0-9]"
                    #        if check_if_exists(self.fitfolder,filenamecheck) != True:
                    #            print "fitting:",filenamecheck
                    #            polyfit(
                    #            foldername=self.fitfolder,
                    #            filename=i[2],
                    #            zeroat=self.nndist[self.shell],
                    #            data=i[3],
                    #            ordermax=orderfit,
                    #            verbose2=False
                    #            )
        #################################################################################
        # SAVING makes sens of both , xdir and quer, are loaded in!
        # do it therefore in externally
        #self.pkl = self.folderparam+'/data.pkl'
        #self.saveDatapd(self.pkl)
        #################################################################################
        print utils.printblue("     pot_parametrize_function complete ...")
        print ""
        return

    def parametrize_tox(self):
        '''
        KEEP IN MIND:
        neither morse nor mc1 are perfect (especially the attractive part)
        deviations from the strong repulsive part:
        with a d of 0.2: deviations of 0.0006 are to be expected in forces for fcc Al
        with a d of 0.3: deviations of 0.005  are to be expected in forces for fcc Al
        THIS ALSO MEANS THAT THIS ERROR WILL PROPAGATE in quantities where we substract
        the longitudinal part:
                tox (will be only "accurate" to this error, only to this displacement)
                tix ( -- "" --  )

                (-> try EVinet,EBrich, other parametrizations)

        For small displacements tox is much bigger than toy toz (factor 3);
        at x displacements of 0.5 tox is similar in magnitude to toy toz; BUT: we have to
        check the statement above once the longitudinal contribution is substractet!
        --> tox including long
        --> tox pure  ( here we should extract the quer very well fitted long part)
        '''
        print utils.printblue("     pot_parametrize.parametrize_tox ...")
        pass


class forcesneighbors_all_disp( object ):
    def __init__(self):
        self.verbose                            = False
        self.folder_displacement_direction_all  = False
        self.pkl                                = False
        self.npz                                = False
        self.alatsstringlist                    = False  # [3.99]
        self.parametrize_npz                    = False
        self.parametrize_data_pkl               = False
        self.parametrize_fml_npz            = False

        self.alle = forcesneighbors()
        self.update_what_to_calculate_and_write_self_xdir_quer()
        self.xdir = False
        self.quer = False
        self.midd = False
        self.dnnd = False

        # folder_displacement_direction_all  ( necessary for loading forces )
        # [ /Users/glensk/Dropbox/Understand_distributions//displacements_dense/Ir/2x2x2sc_quer_3x3x3kp,
        #   /Users/glensk/Dropbox/Understand_distributions//displacements_dense/Ir/2x2x2sc_xdir_3x3x3kp]
        #
        # pkl : /Users/glensk/Dropbox/Understand_distributions//displacements_dense/Ir//2x2x2sc_data_3x3x3kp.pkl
        #
        # alatsstringlist = [ 3.99 ]
        # parametrize = True or False

    def update_what_to_calculate_and_write_self_xdir_quer(self, initialize = False, onlyprint = False):
        if onlyprint == True:
            print "1:npz:",self.alle.parametrize_npz,"self.alle.parametrize_npz"
            print "2:pkl:",self.alle.parametrize_data_pkl,"self.alle.parametrize_data_pkl"
            print "3:anz:",self.alle.parametrize_fml_npz,"self.alle.parametrize_fml_npz"
            return

        self.alle.parametrize_npz                    = self.parametrize_npz
        self.alle.parametrize_data_pkl               = self.parametrize_data_pkl
        self.alle.parametrize_fml_npz            = self.parametrize_fml_npz
        self.alle.pkl = self.pkl
        if self.alle.parametrize_npz == True:
            self.alle.parametrize_data_pkl = True
            self.alle.parametrize_fml_npz = True

        if self.alle.parametrize_data_pkl == True:
            self.alle.parametrize_fml_npz = True


        if self.alle.disp == 'xdir':
            if initialize == True:
                self.xdir = forcesneighbors()
            self.xdir = copy.deepcopy(self.alle)

        if self.alle.disp == 'quer':
            if initialize == True:
                self.quer = forcesneighbors()
            self.quer = copy.deepcopy(self.alle)

        if self.alle.disp == '3nnd':
            if initialize == True:
                self.dnnd = forcesneighbors()
            self.dnnd = copy.deepcopy(self.alle)

        if self.alle.disp == 'midd':
            if initialize == True:
                self.midd = forcesneighbors()
            self.midd = copy.deepcopy(self.alle)
        return


    def load_all_pkl_npz(self):
        print utils.printblue("     pot_parametrize.load_all_pkl_npz ...")
        ###################################################################################
        # load forces
        ###################################################################################
        #print "_________||||||||",self.pkl
        self.update_what_to_calculate_and_write_self_xdir_quer()
        self.update_what_to_calculate_and_write_self_xdir_quer(onlyprint=True)
        for folder in self.folder_displacement_direction_all:  # quer,xdir,midd,3nnd, ...
            print "}}}",folder
            if os.path.isfile(folder+"/EqCoords_direct") != True:

                sys.exit(folder+"/EqCoords_direct does not exist; need routing to get Eqcoords")
                #shutil.copyfile(self.eqcoordsfile_getpot,folder+"/EqCoords_direct")

            # this needs to happen in the corresponding class to have all variables there for xdir and quer
            ###################################################################################
            # load/or parametrize forces for every displacement
            ###################################################################################
            for a in self.alatsstringlist:
                # self.alle.loadforces macht zum schluss
                # --> self.parametrize_lon()  which makes
                # -->       --> infodos
                self.alle.loadforces(folder = folder, a = a)    # macht zum schluss    # self.parametrize_lon()  which makes infodos
                #self.alle.infodos()
                #print 'self.pkl:',self.pkl  # == /Users/glensk/Dropbox/Understand_distributions//displacements_dense/Ir//2x2x2sc_data_3x3x3kp.pkl

                #if self.parametrize_data_pkl == True or os.path.isfile(self.pkl) != True:
                #    self.alle.parametrize_lon()


                # TODO: 1. merge data und fit into one thing: evtl: auch all variablen! die sachen
                # TODO: 2. uebergeben lon_tix den pfad zu den gemergten sachen die es dann einladen muss
                self.update_what_to_calculate_and_write_self_xdir_quer(initialize = True)
                #if self.alle.disp == 'xdir':
                #    self.xdir = forcesneighbors()
                #if self.alle.disp == 'quer':
                #    self.quer = forcesneighbors()
                #if self.alle.disp == '3nnd':
                #    self.dnnd = forcesneighbors()
                #if self.alle.disp == 'midd':
                #    self.midd = forcesneighbors()


                #if os.path.isfile(self.alle.parametrize_fml_npz_path) != True or self.parametrize_fml_npz == True:
                #    needtoparametrize = True
                #if os.path.isfile(self.alle.parametrize_fml_npz_path) == True and self.parametrize_fml_npz != True:
                #print "||||",self.alle.parametrize_fml_npz_path
                #print "||||",self.alle.parametrize_fml_npz
                if self.alle.parametrize_fml_npz == False and os.path.isfile(folder+"/FML.npz") == True:
                    #print 'lonading ',folder,'FML'
                    self.alle.loadFML(folder,"FML")
                self.update_what_to_calculate_and_write_self_xdir_quer()

        ###################################################################################
        # load (and evtl save) data.pkl fits if no parametrization necessary
        ###################################################################################
        if self.parametrize_data_pkl != True and os.path.isfile(self.pkl) == True:
            self.alle.loadDatapd(self.pkl)
        elif self.parametrize_data_pkl == True or os.path.isfile(self.pkl) != True:   # (xdir, quer, 3nnd, 4nnd, False = * = all)
            if len(self.folder_displacement_direction_all) > 1: # we want to have xdir and quer
                self.alle.saveDatapd(self.pkl)
                self.parametrize_data_pkl = False
                self.alle.parametrize_data_pkl = False
            else:
                sys.exit(self.pkl+" was not saved since you would need to load at least xdir and quer but you have only: "+str(self.folder_displacement_direction_all[0]))
        else:
            sys.exit("pdk not loaded ... check why!")
        ###################################################################################
        # load (and evtl. save) FML
        ###################################################################################
        #if self.parametrize_fml_npz != True and os.path.isfile(self.npz) == True:
        #    self.loadFML(self.npz)
        #elif self.parametrize_fml_npz == True or os.path.isfile(self.npz) != True:

        #print 'quer;',self.quer.data[1,'lon','quer','dos'].ix['morse']
        #print 'xdir;',self.xdir.data[1,'lon','quer','dos'].ix['morse']
        #print 'alle;',self.alle.data[1,'lon','quer','dos'].ix['morse']   -- << -- only here

        # das naechste muss eigentlich nur gemacht werden wenn FML noch nicht vorhanden ist und irgendeine parametrisierung angestanden hat (dann sollte man updaten)
        self.update_what_to_calculate_and_write_self_xdir_quer(onlyprint = True)
        #print "self.xdir.FML.shape:",self.xdir.FML.shape
        #print "self.quer.FML.shape:",self.quer.FML.shape
        #print "self.midd.FML.shape:",self.midd.FML.shape
        #sys.exit()
        if self.alle.parametrize_fml_npz == True or os.path.isfile(folder+"/FML.npz") == False:
            self.get_forces_minus_longforce()
        print utils.printblue("     pot_parametrize.load_all_pkl_npz complete ...")
        return


    def get_forces_minus_longforce(self):
        ''' this can only be performed once the .pkl file exists

        creates parametrization for tix vectors
        by substracting quer displacement from xdir

        concrete: go through every displacement in xdirection and substract morse + correction from force on , the ramaining forces are the ti{x,y,z} contributions
        the error will be as high as in the diffshiftedpolydiff/ == self.fit[1,'lon','xdir','dos','morse'].ix['polybesteiff']
        im haupt *.pkl file: /Users/glensk/Dropbox/understand_distributions/displacements_dense/Ir/2x2x2sc_data_3x3x3kp.pkl
             '''
        print ""
        print utils.printblue("     get_forces_minus_longforce (FML) ...")
        # xdir has usually all the data whil quer has only quer data
        import pot_energy_forces
        #print 'pad;',self.pot_energy_forces_params.u1nn_potaddparam

        ##################################################################################
        # go through every displacement and get long forces in xyz direction
        ##################################################################################
        self.update_what_to_calculate_and_write_self_xdir_quer(onlyprint = True)
        self.alle.loadDatapd(self.pkl)
        for folderidx,folder in enumerate(self.folder_displacement_direction_all):
            for a in self.alatsstringlist:
                ##########################################################################
                # get parametrization current long
                ##########################################################################
                self.pot_energy_forces = pot_energy_forces.pot_energy_forces_class()
                #self.pot_energy_forces.params.u1nn_pottype = 'morse'
                #self.pot_energy_forces.params.u1nn_potparam = self.alle.data[1,'lon','quer','dos'].ix['morse']
                #self.pot_energy_forces.params.u1nn_potadd = 'poly'
                #self.pot_energy_forces.params.u1nn_potaddparam = self.alle.fit[1,'lon','quer','dos','morse'].ix['polybeste']
                self.pot_energy_forces.params.u1nn_pottype = 'morse'
                self.pot_energy_forces.params.u1nn_potparam = self.alle.data[1,'lon','quer','alle'].ix['morse']
                self.pot_energy_forces.params.u1nn_potadd = 'spline'
                #self.pot_energy_forces.params.u1nn_potaddparam = self.alle.fit[1,'lon','quer','alle','morse'].ix['splinedata']
                self.pot_energy_forces.params.u1nn_potaddparam = 'lon.quer.alle.morse.splinedata'
                self.pot_energy_forces.pot_parametrize = copy.deepcopy(self.alle)

                def check_if_parametrizations_long_same():
                    if self.alle.fml_u1nn_pottype == self.pot_energy_forces.params.u1nn_pottype and \
                       np.array_equal(self.alle.fml_u1nn_potparam, self.pot_energy_forces.params.u1nn_potparam) == True and \
                       self.alle.fml_u1nn_potadd == self.pot_energy_forces.params.u1nn_potadd and \
                       np.array_equal(self.alle.fml_u1nn_potaddparam, self.pot_energy_forces.params.u1nn_potaddparam) == True and \
                       type(self.alle.fml) != bool and \
                       type(self.alle.FML) != bool:
                       return True
                    else:
                       return False

                ##########################################################################
                # load everything from npz
                # this is important since later we write to this object again
                ##########################################################################
                self.alle.parametrize_fml_npz_path = folder+"/FML.npz"
                self.alle.loadforces(folder = folder, a = a)  # loads npz (forces)
                self.alle.infodos()
                needtoparametrize = None
                if os.path.isfile(self.alle.parametrize_fml_npz_path) != True or self.parametrize_fml_npz == True:
                    needtoparametrize = True
                if os.path.isfile(self.alle.parametrize_fml_npz_path) == True and self.parametrize_fml_npz != True:
                    self.alle.loadFML(folder,"FML")
                    if type(self.alle.FML) != bool:
                        if check_if_parametrizations_long_same() == True:
                            needtoparametrize = False
                        else:
                            needtoparametrize = True
                    else:
                        needtoparametrize = True
                        print utils.printred("self.alle.FML is of type bool ! , need to redo ...")

                ##########################################################################
                # in case it is neccesary to parametrize
                ##########################################################################
                if needtoparametrize == True:
                    self.pot_energy_forces.load_parameters_from_file = False
                    self.pot_energy_forces.cell = self.alle.sc
                    self.pot_energy_forces.coord0_rrel = self.alle.eqcoords
                    self.alle.fml = np.zeros((self.alle.p.shape[0],self.alle.p.shape[1],3))

                    self.alle.fml_u1nn_pottype = 'morse'
                    self.alle.fml_u1nn_potparam = self.alle.data[1,'lon','quer','alle'].ix['morse']
                    #self.alle.fml_u1nn_potadd = 'poly'
                    #self.alle.fml_u1nn_potaddparam = self.alle.fit[1,'lon','quer','dos','morse'].ix['polybeste']
                    self.alle.fml_u1nn_potadd = 'spline'
                    #self.alle.fml_u1nn_potaddparam = self.alle.fit[1,'lon','quer','alle','morse'].ix['splinedata']
                    self.alle.fml_u1nn_potaddparam = 'lon.quer.alle.morse.splinedata'

                    #self.tox = np.zeros((self.xdir.p.shape[0],3))
                    #self.tix = np.zeros((self.xdir.p.shape[0],3))
                    #[ 'fcc', 1, 'lon',32, 24, 126 ],   # [ 1.995,  1.995,  0.   ], [-1.995, -1.995,   0.   ]
                    #[ 'fcc', 2, 'lon',32,  4,  36 ],   # [ 3.99,  0.  ,  0.     ], [-3.99,  0.  ,   0.     ]
                    #[ 'fcc', 3, 'lon',32, 12, 239 ],   # [ 3.99 ,  1.995,  1.995], [-3.99 ,  -1.995, -1.995]
                    #[ 'fcc', 4, 'lon',32,  6, 102 ],   # [ 3.99,  3.99,  0.  ]  , [-3.99,  -3.99,  0.  ]

                    #[ 'fcc', 1, 'tox',32, 8, 203 ],   # [ 0.   ,  2.065,  2.065],[ 0.   , -2.065, -2.065]
                    #[ 'fcc', 2, 'tox',32, 2,  66 ],   # [ 0.  ,  4.13,  0.00], [0,  -4.13  ,   0.     ]
                    #[ 'fcc', 3, 'tox',32, 0,  0  ],   # [ 3.99 ,  1.995,  1.995], [-3.99 ,  -1.995, -1.995]
                    #[ 'fcc', 4, 'tox',32, 0,  0  ],   # [ 3.99,  3.99,  0.  ]  , [-3.99,  -3.99,  0.  ]
                    nsc = 2  # for now on we will just double the supercell, in general maybe factor
                    self.alle.FML = np.zeros((self.alle.p.shape[0],self.alle.p.shape[1]*nsc**3,3))
                    for idx,i in enumerate(self.alle.p):
                        #self.pot_energy_forces.verbose = 2
                        self.pot_energy_forces.coord_cart = i
                        self.pot_energy_forces.pot_run()
                        print idx,"/",self.alle.dstring.shape[0]
                        self.alle.fml[idx] = self.alle.f[idx]-self.pot_energy_forces.forces

                        # for FML:
                        forces = np.copy(self.alle.fml[idx])
                        self.crystalforces = crystal_generator.crystal()
                        self.crystalforces.load_positions_cell(coord_cart = forces , cell = self.alle.sc)

                        # forces big (repeat supercell)
                        self.crystalforcesbig = crystal_generator.supercell()
                        self.crystalforcesbig.create_supercell(  self.crystalforces, nsc, nsc, nsc, newsorting = True, shiftforces = True )
                        self.alle.FML[idx] = self.crystalforcesbig.rcar


                self.update_what_to_calculate_and_write_self_xdir_quer()
                #print 'ntp:',needtoparametrize

                if needtoparametrize == True:
                    self.alle.saveFML(folder,"FML")
        print utils.printblue("         get_forces_minus_longforce complete ...")
        print ""
        return

    def analyze_restforces(self):
        ''' at this point we have everything loaded: self.xdir and self.quer '''
        print ""
        print ""
        print utils.printpink("     analyze_restforces ...")
        for folderidx,folder in enumerate(self.folder_displacement_direction_all):
            for a in self.alatsstringlist:
                self.alle.loadData(     folder, 'nnforces' )  # loads npz (forces)
                self.alle.loadFML(  folder, 'FML'  )
                print "     direction:",self.alle.disp
                self.alle.angles = np.zeros((self.alle.Pr.shape[0],self.alle.Pr.shape[1]))
                self.alle.fprojlongold = np.zeros((self.alle.Pr.shape[0],self.alle.Pr.shape[1],3))
                self.alle.fprojlongnew = np.zeros((self.alle.Pr.shape[0],self.alle.Pr.shape[1],3))
                for dispidx,i in enumerate(self.alle.p): # geht ueber jedes displacement
                    for atomidx in np.arange(self.alle.Pr.shape[1]): # geht ueber jedes atom der vergroesserten superzelle
                        # project vector onto old longvec which is for all 1nndisp the quer direction
                        # np.array([ utils.coord_transform_to_quer(i) for i in pot.disp.xdir.FML[:,60]])

                        # project vector onto a referece frame which is spanned by the new longitudinal vecotr
                        # it would be best to save this vector as a function of the winkel
                        # a) get winkel between old longvec and current longvec for every displacement for every atom of interest
                        # b) project force onto old longvec
                        # c) project force onto new longvec the new reference frame which is turned by the winkel
                        # d) mache das gleiche fuer die tox
                        # e) vergleiche das ganze zur 3x3x3 sc
                        # f) compare the forces as function of winkel for diffrent atoms
                        # newvec = pot.disp.xdir.Pr[1,60]

                        # a)
                        # pot.disp.xdir.angles[:,24] ist der winkel auf das [2.065, 2.065, 0 ] atom wenn wir in x richtug auslenken
                        self.alle.angles[dispidx,atomidx] = utils.anglevec(self.alle.Pr0[atomidx],self.alle.Pr[dispidx,atomidx])

                        # b) + c)
                        if abs(self.alle.FML[dispidx,atomidx][2]) > 1e-9:  # wir gehen hier ueber jedes atom
                            self.alle.fprojlongold[dispidx,atomidx] = np.array([0.0,0.0,0.0])
                            self.alle.fprojlongnew[dispidx,atomidx] = np.array([0.0,0.0,0.0])
                            #print "check  :",self.alle.FML[dispidx,atomidx][2]
                            #print "atomidx:",atomidx
                            #print "dispidx:",dispidx
                        else:
                            self.alle.fprojlongold[dispidx,atomidx] = utils.coord_transform_to_quer(self.alle.FML[dispidx,atomidx])
                            self.alle.fprojlongnew[dispidx,atomidx] = utils.coord_transform_to_quer(self.alle.FML[dispidx,atomidx],self.alle.angles[dispidx,atomidx])

                self.update_what_to_calculate_and_write_self_xdir_quer()
                self.alle.saveFML(folder,"FML")

                ##########################################################################
                # analysis tix
                ##########################################################################
                disp = [ 'xdir', 'quer' ]
                atomslon = [ 24, 60 ]

                #for d in disp:
                #    for at in atomslon:
                #        if
                #        forces = eval('self.'+d+'.FML[:,'+str(atomslon)+'][:,'+xy+']')
                def save_analyze(disp,atom,force_in = False, xfunc = False):
                    if xfunc not in [ 'winkel', 'longvec' ]:
                        sys.exit("xfunc !! really not in winkel / longvec")

                    if disp not in [ 'xdir', 'quer', 'midd' ]:
                        sys.exit("disp !! really not in xdir quer")
                    if force_in not in [ 'xyz', 'quo', 'qun' ]:
                        sys.exit("force_in !! really not in xyz, quo, qun")
                    if xfunc == 'winkel':
                        if disp == 'xdir': xfunctake = self.xdir.angles[:,atom]
                        if disp == 'quer': xfunctake = self.quer.angles[:,atom]
                        if disp == 'midd': xfunctake = self.midd.angles[:,atom]
                    if xfunc == 'longvec':
                        if disp == 'xdir': xfunctake = self.xdir.dnorm
                        if disp == 'quer': xfunctake = self.quer.dnorm
                        if disp == 'midd': xfunctake = self.midd.dnorm

                    if force_in == 'xyz':
                        getfrom = 'FML'
                    if force_in == 'quo':
                        getfrom = 'fprojlongold'
                    if force_in == 'qun':
                        getfrom = 'fprojlongnew'

                    print "disp:",disp
                    print "getfrom:",getfrom
                    print "atom:",atom
                    forcesx = eval('self.'+disp+'.'+getfrom+'[:,'+str(atom)+'][:,0]')
                    forcesy = eval('self.'+disp+'.'+getfrom+'[:,'+str(atom)+'][:,1]')


                    folder = '/Users/glensk/Dropbox/Understand_distributions/analyze/'+str(self.alle.element)+'/'
                    if os.path.isdir(folder) != True:
                        os.makedirs(folder)
                    filenamex = folder+force_in+"_"+xfunc+"_"+disp+"_"+str(atom)+"_x.dat"
                    filenamey = folder+force_in+"_"+xfunc+"_"+disp+"_"+str(atom)+"_y.dat"
                    print filenamex
                    print filenamey
                    print xfunctake.shape,forcesx.shape,forcesy.shape
                    if xfunctake.shape != forcesx.shape:
                        print "xfunctake:"
                        print xfunctake
                        print "forcesx:"
                        print forcesx
                        sys.exit("nope")
                    if os.path.isfile(filenamex) != True:
                        np.savetxt(filenamex,np.transpose([xfunctake,forcesx]))
                    if os.path.isfile(filenamey) != True:
                        np.savetxt(filenamey,np.transpose([xfunctake,forcesy]))
                    return

                if self.alle.disp == 'xdir':
                    save_analyze('xdir',8 ,'xyz','winkel')  # tox
                    save_analyze('xdir',8 ,'quo','winkel')  # tox
                    save_analyze('xdir',24,'xyz','winkel')
                    save_analyze('xdir',60,'xyz','winkel')
                    save_analyze('xdir',24,'quo','winkel')
                    save_analyze('xdir',60,'quo','winkel')
                    save_analyze('xdir',24,'qun','winkel')
                    save_analyze('xdir',60,'qun','winkel')

                    save_analyze('xdir',8 ,'xyz','longvec')  # tox
                    save_analyze('xdir',8 ,'quo','longvec')  # tox
                    save_analyze('xdir',24,'xyz','longvec')
                    save_analyze('xdir',60,'xyz','longvec')
                    save_analyze('xdir',24,'quo','longvec')
                    save_analyze('xdir',60,'quo','longvec')
                    save_analyze('xdir',24,'qun','longvec')
                    save_analyze('xdir',60,'qun','longvec')

                if self.alle.disp == 'midd':
                    save_analyze('midd',8 ,'xyz','longvec')  # tox
                    save_analyze('midd',8 ,'quo','longvec')  # tox
                    save_analyze('midd',24,'xyz','longvec')
                    save_analyze('midd',60,'xyz','longvec')
                    save_analyze('midd',24,'quo','longvec')
                    save_analyze('midd',60,'quo','longvec')
                    save_analyze('midd',24,'qun','longvec')
                    save_analyze('midd',60,'qun','longvec')

                if self.alle.disp == 'quer':
                    save_analyze('quer',8 ,'xyz','winkel') # tox
                    save_analyze('quer',24,'xyz','winkel')
                    save_analyze('quer',60,'xyz','winkel')
                    save_analyze('quer',60,'quo','winkel')
                    save_analyze('quer',60,'qun','winkel')

                    save_analyze('quer',8 ,'xyz','longvec') # tox
                    save_analyze('quer',24,'xyz','longvec')
                    save_analyze('quer',60,'xyz','longvec')
                    save_analyze('quer',60,'quo','longvec')
                    save_analyze('quer',60,'qun','longvec')

        return


    ######################################################################################

    def getforces(self):
        '''
        KEEP IN MIND:
        neither morse nor mc1 are perfect (especially the attractive part)
        deviations from the strong repulsive part:
        with a d of 0.2: deviations of 0.0006 are to be expected in forces for fcc Al
        with a d of 0.3: deviations of 0.005  are to be expected in forces for fcc Al
        THIS ALSO MEANS THAT THIS ERROR WILL PROPAGATE in quantities where we substract
        the longitudinal part:
                tox (will be only "accurate" to this error, only to this displacement)
                tix ( -- "" --  )

                (-> try EVinet,EBrich, other parametrizations)

        For small displacements tox is much bigger than toy toz (factor 3);
        at x displacements of 0.5 tox is similar in magnitude to toy toz; BUT: we have to
        check the statement above once the longitudinal contribution is substractet!
        --> tox including long
        --> tox pure  ( here we should extract the quer very well fitted long part)
        '''
        #################################################
        # get xdir or quer
        #################################################
        xdir = False
        quer = False
        q111 = False
        xdirfolder = False
        longreffolder = False
        querxdir = os.path.basename(os.getcwd())

        if 'xdir' in querxdir: xdir = True
        if 'quer' in querxdir: quer = True
        if 'q111' in querxdir: q111 = True
        if '3NNdir' in querxdir: xdir = True
        if xdir == True and quer == True:
            sys.exit("xdir and quer found in pwd, can just understand one")
        if xdir != True and quer != True and q111 != True:
            sys.exit("neither xdir, quer nor q111 found in pwd")

        #################################################
        # if xdir get corresponding quer folder
        #################################################
        def replace_right(source, target, replacement, replacements=None):
            return replacement.join(source.rsplit(target, replacements))
        # a) get NN atoms of (0/0/0) atom by using eqcoords
        # b) take the atom on (2.065,2.065,0) and opposite (-2.065,-2.065,0)
        # this here is a list for our current positions in fcc
        if self.atoms == 108:  # for fcc
            atfuncpos = 81      # [ 2.065,  2.065,  0.   ]
            atfuncneg = 105     # [ 10.325,  10.325,   0.   ]
            atlong2nnpos = 9       # [ 4.13,  0.  ,  0.  ]
            atlong2nnneg = 18      # [ 8.26,  0.  ,  0.  ]

            # works
            attox1nnpos = 27       # [ 0.   ,  2.065,  2.065]
            attox1nnneg = 27       # [ 0.   ,  2.065,  2.065]
            attox2nnpos = 3        # [ 0.  ,  4.13,  0.  ] or atom 1: [ 0.  ,  0.  ,  4.13])
            attox2nnneg = 3        # [ 0.  ,  4.13,  0.  ] or atom 1: [ 0.  ,  0.  ,  4.13])


            attix1nnpos = 87       # [2.07   10.32    0.  ]
            attix1nnneg = 99       # [2.07   10.32    0.  ]
            attiy1nnpos = 87       # [2.07   10.32    0.  ]
            attiy1nnneg = 99       # [2.07   10.32    0.  ]

            attix1nnpos = 81       # [2.07   2.07    0.  ]  # this atoms correspond to x/y basis for tix and not basis parallel and senkr to vec0
            attix1nnneg = 81       # [10.32  2.07    0.  ]
            attiy1nnpos = 81       # [2.07   2.07    0.  ]
            attiy1nnneg = 81       # [10.32  2.07    0.  ]

            attix1nnpos = 87       # [2.07   2.07    0.  ]  # this atoms correspond to x/y basis for tix and not basis parallel and senkr to vec0
            attix1nnneg = 87       # [10.32  2.07    0.  ]
            attiy1nnpos = 87       # [2.07   2.07    0.  ]
            attiy1nnneg = 87       # [10.32  2.07    0.  ]

            atlong3nnpos = 36       # array([ 4.13 ,  2.065,  2.065])
            atlong3nnneg = 53       # array([  8.26 ,  10.325,  10.325])


        if self.atoms == 32:    # for fcc 2x2x2sc
            atfuncpos = 24      # [ 2.065,  2.065,  0.   ]
            atfuncneg = 30      # [ 6.195,  6.195,  0.   ]
            atlong2nnpos = 4             # [ 4.13,  0.  ,  0.  ]
            atlong2nnneg = 36             # [ -4.13,  0.  ,  0.  ]

            # not yet
            attox1nnpos  =  8       # [ 0.     ,  2.065  ,  2.065  ]
            attox1nnneg  =  8       # [ 0.     ,  2.065  ,  2.065  ]
            attox2nnpos  =  1       # [ 0.   ,  0.   ,  4.13 ] same as [ 0.   ,  4.13 ,  0.   ],
            attox2nnneg  =  1       # [ 0.   ,  0.   ,  4.13 ] same as [ 0.   ,  4.13 ,  0.   ],
            attiatom1nn  =  16     # [ 2.065  ,  0.     ,  2.065  ]  (does not have a tox part)

            attix1nnpos = 26       # [2.07   6.07    0.  ]  # this atoms correspond to x/y basis for tix and not basis parallel and senkr to vec0
            attix1nnneg = 26       # [2.07   6.07    0.  ]
            attiy1nnpos = 26       # [2.07   6.07    0.  ]
            attiy1nnneg = 26       # [2.07   6.07    0.  ]


        if self.atoms == 54:  # bcc 3x3x3sc
            atfuncpos = 27
            atfuncneg = 53
            atlong2nnpos = 9
            atlong2nnneg = 18

            # tox
            attox1nnpos = 45        # in bcc: das 27 atom hat dann gemischt 45:  [ 7.725  ,  1.545  ,  1.545  ],
            attox1nnneg = 45
            attox2nnpos = 3        #
            attox2nnneg = 3        #


        # here we should do a loop over the corresponding neighbor
        def printatom(atom,contr,funcalld):
            print ""
            print utils.printblue("###############################################################"*2)
            print "SHELL:",atom,"    CONTRIBUTION:",contr,"   d:",funcalld
            print utils.printblue("###############################################################"*2)



        self.fitfolder = self.folderparam+'/nnforces/'
        if type(self.DOScut) != bool: #self.DOScut is a numpy array
            if type(self.nnforcesaddname) == bool:
                self.fitfolder = self.folderparam+'/nnforcesDOScut/'
            else:
                self.fitfolder = self.folderparam+'/nnforces'+self.nnforcesaddname+'/'


        if os.path.isdir(self.fitfolder) != True:
            os.makedirs(self.fitfolder)
        ############################################################################
        # start schleife
        ############################################################################
        to = [ 'tox', 'toy', 'toz' ]
        ti = [ 'tix', 'tiy', 'tiz' ]
        toti = [ 'tox', 'toy', 'toz', 'tix', 'tiy', 'tiz' ]
        # 0.25 and 0.35 need to be included AS FIRST d's !!! since thos are currently reference for tox
        funcalldall = [ 0.25, 0.35, 0.15, 0.2, 0.25, 0.3, 0.35 ]
        #funcalldall = [ 0.25, 0.35 ] #, 0.2, 0.25, 0.3, 0.35 ]
        funcalldall = [ 0.15, 0.2 ] #, 0.3 ]
        funcalldall = [ 0.25, 0.35 ]
        shellall       = [ 1,2,3]
        shellall       = [ 1]
        # hier ist entweder xdir == True oder quer == True or q111 == True
        for x in shellall:  # schell 1,2,3  ( == first NN, second NN, third NN ...)
            shell = str(int(x))
            ############################################################################
            # set parameters for schleifen
            ############################################################################
            # GETFORCESVEC1
            # LONTOXTIX
            #shellall       = [ 3]
            #shellall       = [ 1]
            if self.struct == 'fcc' and quer == True:   shellall = [ 1]
            if self.struct == 'fcc' and xdir == True:   shellall = [ 1]

            if xdir == True and x == 1:                 contrib = [ 'lon' , 'tox' ,'toy','toz' ,'tix', 'tiy' ]  # hier lon fuer 2NN
            if xdir == True and x == 1:                 contrib = [ 'lon' ,'tox' ] # , 'tox' ,'toy','toz' ] # ,'tix', 'tiy' ]  # hier lon fuer 2NN
            if xdir == True and x == 1:                 contrib = [ 'lon' ] # , 'tox' ,'toy','toz' ] # ,'tix', 'tiy' ]  # hier lon fuer 2NN
            #if xdir == True and x == 1:                 contrib = [ 'tix' ] # ,'tix', 'tiy' ]  # hier lon fuer 2NN
            if xdir == True and x == 2:                 contrib = [ 'lon' , 'tox' ] #, 'toy', 'toz', 'tix', 'tiy']  # hier lon fu
            if xdir == True and x == 2:                 contrib = [ 'lon' ]

            if quer == True:                            contrib = [ 'lon' , 'tox', 'toy', 'toz' , 'tix' , 'tiy' ]
            if quer == True:                            contrib = [ 'lon' , 'tix' , 'tiy' ]
            #if quer == True:                            contrib = [ 'lon' , 'tox', 'toy', 'toz' ]
            if quer == True:                            contrib = [ 'lon' ]
            #if quer == True:                            contrib = [ 'tix' ] #,'tix', 'tiy' ]
            if q111 == True:                            contrib = [ 'lon' , 'tox' ]

            if type(self.DOScut) != bool:
                # currently unclear why but I get always errors with tox when DOScut
                if xdir == True and x == 1:                 contrib = [ 'lon' ] # , 'tox' ,'toy','toz' ] # ,'tix', 'tiy' ]  # hier lon fuer 2NN
                if quer == True:                            contrib = [ 'lon' ]

            for contr in contrib:    # 'lon', 'tox', 'toy', 'toz', 'tix', 'tiy', 'tiz'
                if contr != 'lon': funcalldall = [ 1.0 ]


                printatom(x,contr,'-')

                if x == 1 and contr == 'lon' and self.struct == 'fcc' and xdir == True:
                    print "THIS WILL BE DONE IN QUER"
                    pass


                ### x == 1
                if x == 1 and contr == 'lon':
                    atpos = atfuncpos
                    atneg = atfuncneg
                if x == 1 and contr == 'tox':
                    atpos = attox1nnpos
                    atneg = attox1nnneg
                if x == 1 and contr == 'toy':
                    atpos = attox1nnpos
                    atneg = attox1nnneg
                if x == 1 and contr == 'toz':
                    atpos = attox1nnpos
                    atneg = attox1nnneg
                if x == 1 and contr == 'tix':
                    atpos = attix1nnpos
                    atneg = attix1nnneg
                if x == 1 and contr == 'tiy':
                    atpos = attix1nnpos
                    atneg = attix1nnneg

                ### x == 2
                if x == 2 and contr == 'lon':
                    atpos = atlong2nnpos
                    atneg = atlong2nnneg
                if x == 2 and contr == 'tox':
                    atpos = attox2nnpos
                    atneg = attox2nnneg

                ### x == 2
                if x == 3 and contr == 'lon':
                    atpos = atlong3nnpos
                    atneg = atlong3nnneg


                #print "XXX: x",x,"CONTR:",contr,"atpos:",atpos,"self.p[0,atpos]:",self.p[0,atpos]
                self.nnxdist = np.linalg.norm(self.p[0,atpos] - np.array([0.0,0.0,0.0]))


                if len(self.p[0]) == 108:
                    self.nnxdist = np.linalg.norm(self.p[0][27])
                if len(self.p[0]) == 32:
                    self.nnxdist = np.linalg.norm(self.p[0][16])

                if x == 2:
                    self.nnxdist = self.p[0,1][2]

                if x == 3:
                    self.nnxdist = np.linalg.norm(n.p[0][atlong3nnpos])

                print "self.nnxdist:",self.nnxdist
                if x == 1:
                    if self.nnxdist > 4.0:
                        sys.exit("self.nnxdist > 4.0")
                    if self.nnxdist < 2.0:
                        sys.exit("self.nnxdist < 2.0")


                # pos side + left side + 0.0 -> 1+2.*
                self.funcpos = np.zeros((np.array(self.folder).shape[0],2))
                self.funcneg = np.zeros((np.array(self.folder).shape[0],2))

                ###########################################################
                # get all posvec and all negvec
                ###########################################################
                posvecall = np.zeros((len(self.dvec),3))
                negvecall = np.zeros((len(self.dvec),3))

                ###########################################################
                # get vectors (go through every dvec)
                ###########################################################
                for idx,i in enumerate(self.dvec):
                    #################################################################
                    # pos displacement
                    #################################################################
                    # small
                    coord_cart_pos = np.copy(self.p[idx])
                    self.crystalpossmall = crystal_generator.crystal()
                    self.crystalpossmall.load_positions_cell(coord_cart = coord_cart_pos, cell = self.sc)

                    # big (repeat supercell)
                    self.crystalposbig = crystal_generator.supercell()
                    nsc = 2  # for now on we will just double the supercell, in general maybe factor
                    self.crystalposbig.create_supercell(  self.crystalpossmall, nsc, nsc, nsc, newsorting = True )
                    self.lookat = copy.deepcopy(self.crystalposbig)
                    self.crystalpossmall.center_atoms_around_atom(0,coord_cart=self.crystalposbig.rcar,cell=self.crystalposbig.cellvec)



                    #################################################################
                    # neg displacement
                    #################################################################
                    # small
                    coord_cart_neg = np.copy(self.p[idx])
                    self.crystalnegsmall = crystal_generator.crystal()
                    self.crystalnegsmall.load_positions_cell(coord_cart = coord_cart_neg, cell = self.sc)

                    # big (repeat supercell)
                    self.crystalnegbig = crystal_generator.supercell()
                    nsc = 2  # for now on we will just double the supercell, in general maybe factor
                    self.crystalnegbig.create_supercell(  self.crystalnegsmall, nsc, nsc, nsc, newsorting = True )
                    self.crystalnegbig.center_atoms_around_atom(0,coord_cart=self.crystalnegbig.rcar,cell=self.crystalnegbig.cellvec)


                    posvec = self.crystalposbig.rcar[atpos]
                    negvec = self.crystalnegbig.rcar[atneg]
                    posvecall[idx] = posvec
                    negvecall[idx] = negvec

                for idx,i in enumerate(self.dvec):
                    #################################################################
                    # pos displacement
                    #################################################################
                    # small
                    coord_cart_pos = np.copy(self.p[idx])
                    self.crystalpossmall = crystal_generator.crystal()
                    self.crystalpossmall.load_positions_cell(coord_cart = coord_cart_pos, cell = self.sc)

                    # big (repeat supercell)
                    self.crystalposbig = crystal_generator.supercell()
                    nsc = 2  # for now on we will just double the supercell, in general maybe factor
                    self.crystalposbig.create_supercell(  self.crystalpossmall, nsc, nsc, nsc, newsorting = True )
                    self.lookat = copy.deepcopy(self.crystalposbig)
                    self.crystalpossmall.center_atoms_around_atom(0,coord_cart=self.crystalposbig.rcar,cell=self.crystalposbig.cellvec)



                    #################################################################
                    # neg displacement
                    #################################################################
                    # small
                    coord_cart_neg = np.copy(self.p[idx])
                    self.crystalnegsmall = crystal_generator.crystal()
                    self.crystalnegsmall.load_positions_cell(coord_cart = coord_cart_neg, cell = self.sc)

                    # big (repeat supercell)
                    self.crystalnegbig = crystal_generator.supercell()
                    nsc = 2  # for now on we will just double the supercell, in general maybe factor
                    self.crystalnegbig.create_supercell(  self.crystalnegsmall, nsc, nsc, nsc, newsorting = True )
                    self.crystalnegbig.center_atoms_around_atom(0,coord_cart=self.crystalnegbig.rcar,cell=self.crystalnegbig.cellvec)


                    posvec = self.crystalposbig.rcar[atpos]
                    negvec = self.crystalnegbig.rcar[atneg]

                    ###########################################################
                    # BEDINGUNGSHOW
                    ###########################################################
                    bedingungshow = False
                    if idx == 0 or idx == 1 or idx == 2 or np.linalg.norm(self.p[idx,0]) == 0.3 or np.linalg.norm(self.p[idx,0]) == 1.0 or idx == len(self.dvec)-1:
                        bedingungshow = True
                    if idx == 0 or idx == len(self.dvec)-1:
                        bedingungshow = True
                    bedingungshow = False


                    if bedingungshow:
                        stratpos = str(atpos)
                        stratneg = str(atneg)
                        if len(stratpos) == 2: stratpos = " "+stratpos
                        if len(stratneg) == 2: stratneg = " "+stratneg
                        strpixpos = self.p[idx,atpos]
                        strpixneg = self.p[idx,atneg]
                        if len(strpixpos) > len(strpixneg):
                            strpixneg = strpixneg+" "*len(strpixpos)-len(strpixneg)
                        if len(strpixneg) > len(strpixpos):
                            strpixpos = strpixpos+" "*len(strpixneg)-len(strpixpos)

                        print utils.printred("DISP:"+str("")+"  LONGPART:"+str(self.dvec[idx]))
                        print "idx:",idx,"self.p[idx,0]:",self.p[idx,0],"atpos:",stratpos,"self.p[idx,atpos:",strpixpos,"posvec:",posvec\
                                ,"|posvec|:",np.linalg.norm(posvec),"force"+stratpos+":",self.f[idx][atpos]
                        print "idx:",idx,"self.p[idx,0]:",self.p[idx,0],"atneg:",stratneg,"self.p[idx,atneg:",strpixpos,"negvec:",negvec\
                                ,"|negvec|:",np.linalg.norm(negvec),"force"+stratneg+":",self.f[idx][atneg]
                    if self.verbose:
                        print "neg:",self.p[idx,atneg],abziehen,self.p[idx,0]
                        print self.dstring[idx],"pos:",posvec,np.linalg.norm(posvec),\
                                "f:",self.f[idx,atpos]
                        print self.dstring[idx],"neg:",negvec,np.linalg.norm(negvec),\
                                "f:",self.f[idx,atneg]

                    self.funcneg[idx,0] = np.linalg.norm(posvecall[idx])
                    self.funcpos[idx,0] = np.linalg.norm(negvecall[idx])
                    print "idx:",idx,"self.funcneg[idx],self.funcpos[idx]",self.funcneg[idx],self.funcpos[idx]
                    #print "self.funcneg[idx,0]:",self.funcneg[idx,0]
                    #print "self.funcpos[idx,0]:",self.funcpos[idx,0]

                    #print self.dvec
                    if self.funcpos[idx,0] == self.funcneg[idx,0]:
                        if self.dnorm[idx] == 0.0:
                            pass
                        else:
                            print "x(==shell):",x
                            print "idx:",idx
                            print "self.funcpos[idx]",self.funcpos[idx]
                            print "self.funcneg[idx]",self.funcneg[idx]
                            print "self.funcpos[idx,0]",self.funcpos[idx,0]
                            print "self.funcneg[idx,0]",self.funcneg[idx,0]
                            print "self.dnorm[idx]:",self.dnorm[idx]

                            print "this happens if atoms are at boundary or outside of supercell, then the mapping does not work....well, kindoff"
                            print "it would be best to increase the size of the supercell from the beginnin on, then there would not be no mapping back to shorter place"
                            print "dnorm",idx, self.dnorm[idx]
                            sys.exit("both equal")

                    if x == 1 or x == 2:
                        self.funcneg[idx,1] = -np.linalg.norm(self.f[idx,atpos])  # - zeichen da repulsive kraft
                        self.funcpos[idx,1] = np.linalg.norm(self.f[idx,atneg])   # + bleibt da attraktive kraft
                    if x == 3:
                        self.funcneg[idx,1] = np.linalg.norm(self.f[idx,atpos])  # - zeichen da repulsive kraft
                        self.funcpos[idx,1] = -np.linalg.norm(self.f[idx,atneg])   # + bleibt da attraktive kraft

                    if self.dvec[idx,0] < 0:
                        if contr == 'lon':
                            self.funcneg[idx,1] =  np.linalg.norm(self.f[idx,atpos])  # - zeichen da repulsive kraft
                            self.funcpos[idx,1] = -np.linalg.norm(self.f[idx,atneg])   # + bleibt da attraktive kraft
                        if contr in toti:
                            self.funcneg[idx,0] = -np.linalg.norm(posvecall[idx])
                            self.funcpos[idx,0] = np.linalg.norm(negvecall[idx])
                    #        self.funcneg[idx,1] = -np.linalg.norm(self.f[idx,atneg])  # - zeichen da repulsive kraft
                    #        self.funcpos[idx,1] = np.linalg.norm(self.f[idx,atpos])   # + bleibt da attraktive kraft


                    if bedingungshow:
                        print "    f1(kurze seite):self.f[idx,atpos]:",self.f[idx,atpos],"|self.f[idx,atpos]|:",np.linalg.norm(self.f[idx,atpos]),"--> self.funcneg[idx]:",self.funcneg[idx]
                        print "    f2(lange seite):self.f[idx,atneg]:",self.f[idx,atneg],"|self.f[idx,atneg]|:",np.linalg.norm(self.f[idx,atneg]),"--> self.funcneg[idx]:",self.funcpos[idx]
                        print " "

                    ########################################################################################
                    ########################################################################################
                    ########################################################################################
                    ## CONTRIBUTINON TO TI
                    ########################################################################################
                    ########################################################################################
                    ########################################################################################
                    if contr in toti:
                        # for to{x,y,z}/ti{x,y,z} we want now to have here the
                        #   lon_{1,2}nn_neg_rlv_0.3_mc1_fit_mc1_xxxCOEFSxxxx
                        #   lon_{1,2}nn_pos_rlv_0.3_mc1_fit_mc1_xxxCOEFSxxxx (NOT)
                        #   FROM THE QUER FOLDER
                        #def longreffolder(
                        if x == 1:
                        #if xdir == True and x == 1:
                            if self.struct == 'fcc':
                                longreffolder = replace_right(os.path.abspath(os.getcwd()),"xdir","quer",1)
                                #longreffolder = replace_right(os.path.abspath(os.getcwd()),"xdir","xdir",1)
                            if self.struct == 'bcc':
                                longreffolder = replace_right(os.path.abspath(os.getcwd()),"xdir","q111",1)
                        #if xdir == True and x == 2:
                        if x == 2:
                            longreffolder = os.getcwd()

                        if type(longreffolder) == bool:
                            sys.exit("longreffolder not found, TYPE bool")
                        if os.path.isdir(longreffolder) != True:
                            sys.exit(longreffolder+" not found")

                        #print "longreffolder:",longreffolder

                        ##################################################################
                        # here we need a module which looks for the longitudinal
                        # reference ( for tox_1nn and for tox_2nn we need the longitudinal
                        # reference), gets the parametrization (morse,mc1,polynomial)
                        # and substracts this force from the current atom
                        ##################################################################
                        # wir sollten hier nicht nur nn_neg sondern ganz nn nehmen evtl.
                        # naja evtl. ist 0.3 doch besser
                        def getlongforce(longreffolder,x,posvec,tito=False,dispdir=False,verbose=False,nnxdist=False,posvecall=False):
                            '''
                            - longreffolder is the folder where it is searchd for
                            the long refence
                            - x is the shell it is looked for (e.g. lon_2nn) x=2
                            '''
                            # we will usually take the accurately parametrized part
                            # in case of 1NN tox berechnung:
                            # in x displacement the longvec will always just be
                            # attractive (never repulsive) --> therefore pos
                            # since vector is longer than equilibrium vector
                            if type(dispdir) == bool:
                                sys.exit("dispdir is bool")
                            if len(dispdir) != 3:
                                sys.exit("len dispdir != 3")
                            xdir, quer, q111 = dispdir

                            to = [ 'tox', 'toy', 'toz' ]
                            ti = [ 'tix', 'tiy', 'tiz' ]
                            if type(tito) == bool:
                                sys.exit("tito is bool")
                            if tito in ti or tito in to:
                                pass
                            else:
                                sys.exit("tito not known")
                            if tito in to:  # here we are only interested in the
                                            # attractive part since d(long) > d1NN
                                sided = '0.35'
                            if tito in ti:  # here we are only interested in the
                                            # attractive part since d(long) > d1NN
                                sided = '0.25' # reicht locker bei jetzigen auslenkungen
                                # lon_1nn_pos_mc1_fit_mc1_0.442778_1.039120_2.920351_0.801952_-0.23830

                            if xdir == True and tito in ti:
                                sided = '0.35'  # macht bei hohen d's "weniger" feherl
                            if type(nnxdist) == False:
                                sys.exit("nnxdist missing")
                            allverh =  np.array([ np.linalg.norm(j)/nnxdist for j in posvecall ])
                            allverhgetmax = abs(allverh.max()-1.0)
                            allverhgetmin = abs(allverh.max()-1.0)
                            allverhget = allverhgetmax
                            if allverhgetmax > allverhgetmax:
                                allverhget = allverhgetmin

                            #print "nnx:",nnxdist,"posvec:",posvec,np.linalg.norm(posvec),"||| VERH:",np.linalg.norm(posvec)/nnxdist,"|| VERHGET:",allverhget
                            if np.linalg.norm(posvec) >= nnxdist:
                                side = "pos"
                                sided = '_rlv_0.35'
                            else:
                                side = "neg"
                                sided = '_rlv_0.25'
                            #print "allverhget:",allverhget
                            #if allverhget <= 0.35:
                            #    sided = '_rlv_0.35'
                            #if allverhget < 0.25:
                            #    sided = '_rlv_0.25'
                            #if allverhget > 0.35:
                            #    sided = ''  # this will give us all the available points and the corresponding fit




                            parameters_long_neg_search = longreffolder+\
                                "/nnforces/lon_"+str(x)+\
                                "nn_"+side+sided+"_mc1_fit___[0-9-]*"
                            #if tito in ti:
                            #    parameters_long_neg_search = longreffolder+\
                            #        "/nnforces/lon_"+str(x)+\
                            #        "nn_"+side+"_mc1_fit_mc1_[0-9-]*"
                            #in case of 2NN tox berechnung:
                            parameters_long_neg_searchalt = longreffolder+\
                                "/nnforces/lon_"+shell+\
                                "nn_"+side+"_poly_fit___*"
                            parameters_long_neg_file = glob.glob(parameters_long_neg_search)
                            parameters_long_neg_filealt =\
                            glob.glob(parameters_long_neg_searchalt)
                            len1 = len(parameters_long_neg_file)
                            len2 = len(parameters_long_neg_filealt)
                            parameters_long_neg_filetake = False
                            found = False
                            if len1 == 0 and len2 == 0:
                                print "parameters_long_neg_search",parameters_long_neg_search
                                print "parameters_long_neg_searchalt",parameters_long_neg_searchalt
                                print "parameters_long_neg_file",parameters_long_neg_file
                                print "parameters_long_neg_filealt",parameters_long_neg_filealt
                                sys.exit("reference not found")
                            if len1 == 1 and len2 == 1:
                                if "_poly_" in parameters_long_neg_filealt[0] and "_poly_" not in parameters_long_neg_file[0]:
                                    parameters_long_neg_filetake =\
                                        parameters_long_neg_filealt[0]
                                    #print "yo1",type(parameters_long_neg_filealt)
                                if "_poly_" not in parameters_long_neg_filealt[0] and "_poly_" in parameters_long_neg_file[0]:
                                    parameters_long_neg_filetake =\
                                        parameters_long_neg_file[0]
                                    #print "yo2",type(parameters_long_neg_file)

                            if len1 != 0 and len2 != 0 and parameters_long_neg_filetake == False:
                                print parameters_long_neg_file
                                print parameters_long_neg_filealt
                                if "_poly_" in parameters_long_neg_filealt:
                                    print "yo1",type(parameters_long_neg_file)
                                print "1"
                                if "_poly_" in parameters_long_neg_file:
                                    print "yo2",type(parameters_long_neg_file)
                                print "2"
                                sys.exit("two references found!?")
                            if len1 != 0 and len2 != 0 and parameters_long_neg_filetake == False:
                                print parameters_long_neg_file
                                print parameters_long_neg_filealt
                                sys.exit("two references found!?")
                            if len1 == 0 and len2 == 0 and parameters_long_neg_filetake == False:
                                if len2 != 1: sys.exit("len2 is not 1")
                                parameters_long_neg_filetake =\
                                        parameters_long_neg_filealt[0]
                            if len1 != 0 and len2 == 0 and parameters_long_neg_filetake == False:
                                #print(parameters_long_neg_file)
                                #print(parameters_long_neg_filealt)
                                if len1 != 1: sys.exit("len1 is not 1")
                                parameters_long_neg_filetake = parameters_long_neg_file[0]
                            if verbose:
                                print "file:",parameters_long_neg_filetake
                            if parameters_long_neg_filetake == False:
                                sys.exit( \
                                        "parameters_long_neg_filetak not found")
                            parameters_long_neg = False

                            if "mc1" in os.path.basename(parameters_long_neg_filetake):
                                    parameters_long_neg_str = \
                                    parameters_long_neg_filetake.\
                                    split("nn_"+side+sided+\
                                    "_mc1_fit___")[1]
                                    parameters_long_neg =\
                                            utils.string_to_list(parameters_long_neg_str)

                            if "poly" in os.path.basename(parameters_long_neg_filetake):
                                parameters_long_neg_str = \
                                parameters_long_neg_filetake.\
                                split("order_")[1]
                                parameters_long_neg =\
                                        utils.string_to_list(parameters_long_neg_str)
                            #print 'parameters:',parameters_long_neg
                            if parameters_long_neg == False:
                                sys.exit("parameters_long_neg not found")

                            forcelong = False
                            #print "||",parameters_long_neg
                            #print "|||",np.linalg.norm(-posvec)
                            if "mc1" in os.path.basename(parameters_long_neg_filetake):
                                forcelong = getforcepot(-posvec,'mc1',parameters_long_neg)
                            if "poly" in os.path.basename(parameters_long_neg_filetake):
                                #sys.exit("POLY!")
                                #forcelong = getforcepot(-posvec,'mc1',parameters_long_neg)
                                e,forcelong = pot_energy.getefvec(-posvec,parameters_long_neg,pot = False)
                            if type(forcelong) == bool:
                                sys.exit("forcelong not found")
                            return forcelong



                        dispdir = [ xdir, quer, q111 ]
                        if bedingungshow:
                            print utils.printyellow("TO - TI - PART:")
                            print "dispdir:",dispdir

                        # for tox posvec will be enough
                        forcelong = getlongforce(longreffolder,x,posvec,contr,dispdir,verbose=bedingungshow,nnxdist=self.nnxdist,posvecall=posvecall)
                        forcelongneg = getlongforce(longreffolder,x,negvec,contr,dispdir,verbose=bedingungshow,nnxdist=self.nnxdist,posvecall=posvecall)


                        if bedingungshow:
                            print "posvec    ",posvec,"norm:",np.linalg.norm(posvecall[idx])
                            print "negvec    ",negvec,"norm:",np.linalg.norm(negvecall[idx])
                            print "force     ",self.f[idx,atneg]
                            print "forcelong   >",forcelong
                            print "forcelongneg>",forcelongneg


                        #print "pos>>:",posvec
                        #print "neg>>:",negvec
                        # der posvec/negvec sind immer da um die kraft fuer die longitudinale komponente zu bekommen, von daher muessen diese immer zwischen 1.5 udn 3.5 sein
                        if x == 1 and np.linalg.norm(posvec) > 5.6:    # posvec has always to be smaller then equilibrium lattice constant
                            sys.exit("posvec > 5.6:"+str(posvec)+" "+str(np.linalg.norm(posvec)))

                        if x == 1 and np.linalg.norm(posvec) < 1.5:
                            sys.exit("posvec < 1.5:"+str(posvec)+" "+str(np.linalg.norm(posvec)))
                        if x == 1 and np.linalg.norm(negvec) < 1.5:
                            sys.exit("negvec < 1.5:"+str(negvec)+" "+str(np.linalg.norm(negvec)))
                        if x == 1 and np.linalg.norm(negvec) > 5.6:
                            sys.exit("negvec > 5.6:"+str(negvec)+" "+str(np.linalg.norm(negvec)))


                            # in case of tox (2NN)  we need to get the forcelong from
                            # the polyfit! (since lon2NN was parametrized with a polynom)

                        #print "f:",self.f[idx,atneg]
                        #print 'fcl:',forcelong
                        toxrest = self.f[idx,atneg] - forcelong
                        if bedingungshow:
                            print "toxrest:  ",toxrest

                        if bedingungshow:
                            print "toxrestNEW",toxrest


                        # hier fuer tox erstmal nur die x komponente
                        #print "kkk:",self.funcpos[idx],self.funcneg[idx]
                        if contr in to or contr == 'lon':
                            if contr == 'tox':takeidx = 0
                            if contr == 'toy':takeidx = 1
                            if contr == 'toz':takeidx = 2
                            self.funcneg[idx,0] =  -posvec[0]
                            self.funcpos[idx,0] =   negvec[0]
                            self.funcneg[idx,1] =  toxrest[takeidx]
                            self.funcpos[idx,1] = -toxrest[takeidx]
                        if contr in ti:
                            if contr == 'tix':takeidx = 0
                            if contr == 'tiy':takeidx = 1
                            if contr == 'tiz':takeidx = 2

                        if xdir == True and contr == 'toy' or contr == 'toz':
                            if self.dvec[idx,0] >= 0:
                                self.funcneg[idx,0] =  posvec[0]
                                self.funcpos[idx,0] =  -negvec[0]
                                self.funcneg[idx,1] =  -toxrest[takeidx]
                                self.funcpos[idx,1] =  toxrest[takeidx]
                            if self.dvec[idx,0] < 0:
                                self.funcneg[idx,0] =  -posvec[0]
                                self.funcpos[idx,0] =  negvec[0]
                                self.funcneg[idx,1] =  -toxrest[takeidx]
                                self.funcpos[idx,1] =  toxrest[takeidx]

                        if xdir == True and contr == 'tix' or contr == 'tiy':
                            if self.dvec[idx,0] >= 0:
                                self.funcneg[idx,0] =  self.dvec[idx][0]
                                self.funcpos[idx,0] =  self.dvec[idx][0]
                                self.funcneg[idx,1] =  toxrest[takeidx]
                                self.funcpos[idx,1] =  toxrest[takeidx]
                            if self.dvec[idx,0] < 0:
                                self.funcneg[idx,0] =  self.dvec[idx][0]
                                self.funcpos[idx,0] =  self.dvec[idx][0]
                                self.funcneg[idx,1] =  toxrest[takeidx]
                                self.funcpos[idx,1] =  toxrest[takeidx]

                        if quer == True and contr == 'tix' or contr == 'tiy':
                            if self.dvec[idx,0] >= 0:
                                #self.funcneg[idx,0] =  self.dvec[idx][0]
                                #self.funcpos[idx,0] =  self.dvec[idx][0]
                                self.funcneg[idx,0] = math.copysign(np.linalg.norm([self.dvec[idx][0],self.dvec[idx][1],0.0]),self.dvec[idx][0])
                                self.funcpos[idx,0] = math.copysign(np.linalg.norm([self.dvec[idx][0],self.dvec[idx][1],0.0]),self.dvec[idx][0])
                                self.funcneg[idx,1] =  toxrest[takeidx]
                                self.funcpos[idx,1] =  toxrest[takeidx]

                                #self.funcnegquermagnitude[idx,0] = np.linalg.norm(self.dvec[idx])
                                #self.funcposquermagnitude[idx,0] = np.linalg.norm(self.dvec[idx])
                                #self.funcnegquermagnitude[idx,1] = self.funcneg[idx,1]
                                #self.funcposquermagnitude[idx,1] = self.funcpos[idx,1]



                            if self.dvec[idx,0] < 0:
                                #self.funcneg[idx,0] =  self.dvec[idx][0]
                                #self.funcpos[idx,0] =  self.dvec[idx][0]
                                self.funcneg[idx,0] = math.copysign(np.linalg.norm([self.dvec[idx][0],self.dvec[idx][1],0.0]),self.dvec[idx][0])
                                self.funcpos[idx,0] = math.copysign(np.linalg.norm([self.dvec[idx][0],self.dvec[idx][1],0.0]),self.dvec[idx][0])
                                self.funcneg[idx,1] =  toxrest[takeidx]
                                self.funcpos[idx,1] =  toxrest[takeidx]




                        if bedingungshow:
                            print "funcpos[idx]>>:",self.funcpos[idx]
                            print "funcneg[idx]>>:",self.funcneg[idx]
                            print "-----"
                            print ""
                    # schleife ueber toti hoert hier auf

                self.funcneg = utils.remove_duplicates_of_2d_array(self.funcneg)
                self.funcpos = utils.remove_duplicates_of_2d_array(self.funcpos)
                print "self.funcneg:"
                print self.funcneg
                print "self.funcpos:"
                print self.funcpos
                if contr == 'lon':
                    self.funcpos = self.funcneg
                sys.exit()
                # here: in case we do have dvec:
                # bring vectors in right order
                if type(self.DOScut) != bool: #self.DOScut is a numpy array
                    self.funcneg = utils.cut_function_at_DOS(self.funcneg,self.DOScut)
                    self.funcpos = utils.cut_function_at_DOS(self.funcpos,self.DOScut)

                self.funcall = np.concatenate((self.funcneg,self.funcpos))
                if contr == 'lon':
                    zeroat = self.nnxdist
                else:  # 'to{x,y,z}','ti{x,y,z}'
                    zeroat = 0.0
                self.funcall = np.concatenate((self.funcall,np.array([[zeroat,0.0]])))
                #print "sfa: ----------------------"
                #print self.funcall
                #print "sfa: ----------------------"
                #sys.exit()

                self.funcall = self.funcall[self.funcall[:,0].argsort()]
                self.funcneg = np.concatenate((self.funcneg,np.array([[zeroat,0.0]])))
                self.funcpos = np.concatenate((self.funcpos,np.array([[zeroat,0.0]])))
                self.funcneg = self.funcneg[self.funcneg[:,0].argsort()]
                self.funcpos = self.funcpos[self.funcpos[:,0].argsort()]

                if contr == 'lon':
                    self.funcpos = self.funcpos[self.funcpos[:,0] >= self.nnxdist]
                    self.funcneg = self.funcneg[self.funcneg[:,0] <= self.nnxdist]
                #print "self.funcneg:"
                #print self.funcneg
                #print "self.funcpos:"
                #print self.funcpos



                #forlinfit = np.array([self.funcpos[1],self.funcpos[0],self.funcneg[-2]])
                if self.verbose:
                    print "--"
                    print self.funcpos[0]
                    print self.funcpos[1]
                    print self.funcneg[-2]

                print utils.printyellow("################# "+contr+": ##################")
                print self.funcall[:4]
                print self.funcall[-4:]
                print ""


                ###################################################
                ## fit all forces
                ###################################################
                for funcalld in funcalldall:   #     0.2 oder 0.3

                    print utils.printgreen("################# funcalld: "+str(funcalld)+" ########################")
                    self.funcalld = funcalld
                    dstr = str(self.funcalld)
                    # erst an dieser stelle wird das d interessant
                    # get min and max d
                    self.funcallmax = self.nnxdist+self.funcalld*self.nnxdist
                    self.funcallmin = self.nnxdist-self.funcalld*self.nnxdist

                    # get rlv region
                    self.funcnegrlv = \
                            self.funcneg[self.funcneg[:,0]>self.funcallmin]
                    self.funcposrlv = \
                            self.funcpos[self.funcpos[:,0]<self.funcallmax]
                    self.funcallrlv = \
                            np.concatenate((self.funcnegrlv,self.funcposrlv))
                    self.funcallrlv = \
                            self.funcallrlv[self.funcallrlv[:,0].argsort()]

                    self.funcallrlvallpos = \
                            np.concatenate((self.funcnegrlv,self.funcpos))
                    self.funcallrlvallpos = \
                            self.funcallrlvallpos[self.funcallrlvallpos[:,0].argsort()]

                    # save vectors
                    #if self.verbose:
                    print "saving:"
                    print  self.fitfolder+contr+"_"+shell+"nn_all"
                    print  self.fitfolder+contr+"_"+shell+"nn_pos"
                    print  self.fitfolder+contr+"_"+shell+"nn_neg"
                    if type(self.DOScut) == bool:
                        print  self.fitfolder+contr+"_"+shell+"nn_pos_rlv_"+dstr
                        print  self.fitfolder+contr+"_"+shell+"nn_neg_rlv_"+dstr

                    np.savetxt(self.fitfolder+contr+"_"+shell+"nn_all",self.funcall,fmt="%.6f")
                    np.savetxt(self.fitfolder+contr+"_"+shell+"nn_all_shiftedto0",np.transpose([self.funcall[:,0]-self.nnxdist,self.funcall[:,1]]),fmt="%.6f")

                    if quer == True and contr in 'tiy':  # at this point we usually already do have tix
                        tix = np.loadtxt(self.fitfolder+"tix"+"_"+shell+"nn_all")
                        tiy = np.loadtxt(self.fitfolder+"tiy"+"_"+shell+"nn_all")
                        #print "tix:"
                        #print tix
                        #print "tiy:"
                        #print tiy
                        out = np.zeros((tix.shape[0],3))
                        #print "coordinate transformation for plotting quer function"
                        for outidx,outline in enumerate(out):
                            out[outidx][0] = tix[outidx][0]
                            #print  tix[outidx],tiy[outidx],utils.coord_transform_to_quer(np.array([tix[outidx][1],tiy[outidx][1],0.0]))
                            coordtransformed = utils.coord_transform_to_quer(np.array([tix[outidx][1],tiy[outidx][1],0.0]))
                            #print "coordtransformed:",coordtransformed
                            out[outidx][1] = coordtransformed[0]
                            out[outidx][2] = coordtransformed[1]
                        quermagnitude_coordstransformed_x = out[:,:2]
                        quermagnitude_coordstransformed_y = out[:,[0,2]]

                        np.savetxt(self.fitfolder+'tix'+"_"+shell+"nn_all_coordstransformed_quer",quermagnitude_coordstransformed_x,fmt="%.6f")
                        np.savetxt(self.fitfolder+'tiy'+"_"+shell+"nn_all_coordstransformed_quer",quermagnitude_coordstransformed_y,fmt="%.6f")
                        # let's now parametrize this functions
                        polyfit(
                            foldername=self.fitfolder,
                            filename='tix'+"_"+shell+"nn_all_coordstransformed_quer_poly",
                                    zeroat=0.0,
                                    data=quermagnitude_coordstransformed_x)
                        polyfit(
                            foldername=self.fitfolder,
                            filename='tiy'+"_"+shell+"nn_all_coordstransformed_quer_poly",
                                    zeroat=0.0,
                                    data=quermagnitude_coordstransformed_y)



                        ################################################################################################################################
                        # Now we want to get the difference between the xdisplacement and the querdisplacemetn to determine the ti parallel part
                        ################################################################################################################################
                        # go thround every xdir-displacement and get corresponding quer force
                        #for i in self.dvec:  #his is ment to be the dvec in x-dir which goes up to 1.6
                        #for outidx,outline in enumerate(out):
                        #    vec0 = np.array([ 2.065, -2.065,  0.   ])
                        #    vecs = np.array([1.0,1.0,0.0])
                        #    longvec = np.array([-2.065+i[0],2.065,0.0])
                        #    print i,longvec,utils.anyvec_to_inplane_senkr_parallel(longvec,vec0,vecs)
                        #    quervec = utils.anyvec_to_inplane_senkr_parallel()


                    if contr == 'lon':
                        np.savetxt(self.fitfolder+contr+"_"+shell+"nn_pos",self.funcpos,fmt="%.6f")
                        np.savetxt(self.fitfolder+contr+"_"+shell+"nn_neg",self.funcneg,fmt="%.6f")
                        if type(self.DOScut) == bool:
                            np.savetxt(self.fitfolder+contr+"_"+shell+"nn_pos_rlv_"+dstr,self.funcposrlv,fmt="%.6f")
                            np.savetxt(self.fitfolder+contr+"_"+shell+"nn_neg_rlv_"+dstr,self.funcnegrlv,fmt="%.6f")
                            np.savetxt(self.fitfolder+contr+"_"+shell+"nn_all_rlv_"+dstr,self.funcallrlv,fmt="%.6f")
                            np.savetxt(self.fitfolder+contr+"_"+shell+"nn_all_rlv_"+dstr+"_allpos",self.funcallrlvallpos,fmt="%.6f")

                    #####################################################################################################
                    # make fit (it would be best to do this seperately for the left and for the right side
                    #####################################################################################################
                    # fit long forces
                    #####################################################################################################
                    if x == 1 and contr == 'lon':

                        if type(self.DOScut) == bool:
                            # fit neg rlv  (m, mc1)
                            self.funcnegmorse = get_fit_forces_to_pot(
                            foldername=self.fitfolder,
                            filename=contr+"_"+shell+"nn_neg_rlv_"+dstr+"_mc1",
                                    NN=self.nnxdist,pot = 'mc1',
                                    data=self.funcnegrlv)

                            self.funcnegmorse = get_fit_forces_to_pot(
                            foldername=self.fitfolder,
                            filename=contr+"_"+shell+"nn_neg_rlv_"+dstr+"_morse",
                                    NN=self.nnxdist,pot = 'm',
                                    data=self.funcnegrlv)

                            # fit pos rlv (m, mc1)
                            self.funcposmorse = get_fit_forces_to_pot(
                            foldername=self.fitfolder,
                            filename=contr+"_"+shell+"nn_pos_rlv_"+dstr+"_mc1",
                                    NN=self.nnxdist,pot = 'mc1',
                                    data=self.funcposrlv)


                            # fit all rlv (links und rechts den gleichen bereich)
                            self.funcallmorse = get_fit_forces_to_pot(
                            foldername=self.fitfolder,
                            filename=contr+"_"+shell+"nn_all_rlv_"+dstr+"_morse",
                                    NN=self.nnxdist, pot = 'm',
                                    data=self.funcallrlv)

                            self.funcallmorse = get_fit_forces_to_pot(
                            foldername=self.fitfolder,
                            filename=contr+"_"+shell+"nn_all_rlv_"+dstr+"_mc1",
                                    NN=self.nnxdist, pot = 'mc1',
                                    data=self.funcallrlv)


                            # fit allrvl + all pos on the right side
                            self.funcallmorse = get_fit_forces_to_pot(
                            foldername=self.fitfolder,
                            filename=contr+"_"+shell+"nn_all_rlv_"+dstr+"_allpos_morse",
                                    NN=self.nnxdist, pot = 'm',
                                    data=self.funcallrlvallpos)

                            self.funcallmorse = get_fit_forces_to_pot(
                            foldername=self.fitfolder,
                            filename=contr+"_"+shell+"nn_all_rlv_"+dstr+"_allpos_mc1",
                                    NN=self.nnxdist, pot = 'mc1',
                                    data=self.funcallrlvallpos)


                        # fit pos all
                        self.funcallmorse = get_fit_forces_to_pot(
                        foldername=self.fitfolder,
                        filename=contr+"_"+shell+"nn_pos_mc1",
                                NN=self.nnxdist, pot = 'mc1',
                                data=self.funcpos)

                        self.funcallmorse = get_fit_forces_to_pot(
                        foldername=self.fitfolder,
                        filename=contr+"_"+shell+"nn_pos_morse",
                                NN=self.nnxdist, pot = 'm',
                                data=self.funcpos)

                        # fit neg all
                        self.funcallmorse = get_fit_forces_to_pot(
                        foldername=self.fitfolder,
                        filename=contr+"_"+shell+"nn_neg_mc1",
                                NN=self.nnxdist, pot = 'mc1',
                                data=self.funcneg)



                        # plyfits for the pos part
                        for orderfit in np.arange(1,10):
                            # makes something like lon_1nn_pos_poly9_6_fit___157.32847059439_-253.3149803360_172.60763615111_-64.55258342696_14.303612909803_-1.874831186120_0.1343007378748_-0.004043636776 BUT NOT A CORRESPONDING DELTAFOLDER
                            # here the deltas are saved in deltas but no fit is made afterwards
                            polyfit(
                                foldername=self.fitfolder,
                                filename=contr+"_"+shell+"nn_pos",
                                zeroat=zeroat,
                                data=self.funcpos,
                                ordermax=orderfit,
                                verbose2=False
                                )

                        #######################################################
                        # fit all above but for DOScut bereich
                        #######################################################
                        if type(self.DOScut) != bool:
                            # fit neg rlv  (m, mc1)
                            self.funcnegmorse = get_fit_forces_to_pot(
                            foldername=self.fitfolder,
                            filename=contr+"_"+shell+"nn_all_mc1",
                                    NN=self.nnxdist,pot = 'mc1',
                                    data=self.funcall)

                            self.funcnegmorse = get_fit_forces_to_pot(
                            foldername=self.fitfolder,
                            filename=contr+"_"+shell+"nn_all_morse",
                                    NN=self.nnxdist,pot = 'm',
                                    data=self.funcall)


                    else:
                        #####################################################################################################
                        # fit tox forces (da else) ( == immer wenn contribution nicht lon ist )
                        # verstehen nicht wiso hier wider auf lon angespielt wird, sollte hier keinen sinn machen
                        #####################################################################################################
                        if contr == 'lon':  # im fall von x = 2 also second NN
                            #sys.exit("hier sollte NIE noch lon sein, ansonsten verstehen ich was nicht; scheint auch der fall zu sein")
                            pass

                        if contr == 'lon':
                            zeroat = self.nnxdist
                        if contr == 'tox' or contr == 'toy' or contr == 'toz' or contr == 'tix' or contr == 'tiy' or contr == 'tiz':
                            zeroat = 0.0
                        print "zeroat:",zeroat,"self.nnxdist:",self.nnxdist

                        for orderfit in np.arange(1,14):
                            # makes something like lon_1nn_pos_poly9_6_fit___157.32847059439_-253.3149803360_172.60763615111_-64.55258342696_14.303612909803_-1.874831186120_0.1343007378748_-0.004043636776 BUT NOT A CORRESPONDING DELTAFOLDER
                            # here the deltas are saved in deltas but no fit is made afterwards
                            polyfit(
                                foldername=self.fitfolder,
                                filename=contr+"_"+shell+"nn_all",
                                zeroat=zeroat,
                                data=self.funcall,
                                ordermax=orderfit,
                                verbose2=False
                                )


                        #if contr == 'lon':
                        #    # will keinen polyfit ueber die ganze range, macht keinen sinn, hat keinen exponenten wie morse
                        #    polyfit(
                        #    foldername=self.fitfolder,
                        #    filename=contr+"_"+shell+"nn_pos",
                        #            zeroat=zeroat,
                        #            data=self.funcpos)
                    #################################################################
                    # tix
                    #################################################################
                    if xdir == True and contr in ti and x == 1:  # first shell
                        print ""
                        print "Building ti force component in parallel direction ..."
                        print ""
                        # For tox it is easy since there is "only" one tox direction, therefore 2NN should work well when included 2nn_lon and 2nn_tox
                        #
                        # a1) DONE; get the ture ti forces in x-direction
                        # a2) DONE; get parametrization ti of forces in x-direction (just done above)
                        # b1) DONE; get the true ti forces in quer-direction  (ti{x,y}_1nn_all)
                        # b2) DONE; get parametrization ti of forces in quer-direction ( from quer folder)
                        # c)  DONE; get a function which gives for any longvec (rather tvec) the correct projection on senkr/=quer and parallel (see utils)
                        # d)  DONE; when in quer auslenkung, get forces in x and y direction (hopefully this is ok, but should be only a matter of coordinate transformation)
                        #           correctly one should plot this as a function of magnitude in quer direction; do we save this as a function of x coordinate os as a function of magnitued of disp?
                        #           it is easier to save this as a function of x coordinate then no transformation is necessary
                        # e)  DONE; make a seperate quer kraft for plotting which is as function of magnitude
                        # f)  DONE; make a seperate quer kraft for plotting which is as function of magnitude and the force shows in quer and par direction! (par since the force will have
                        #           a y component)
                        # g)  DONE; get the minvalues and maxvalues to which the parallel plot still does have to be extrapolated
                        # h)  DONE; get the difference x-dir ti{x,y} (a2) (=full force on 87)  minus d) for the correct mesh in e) == parkraft
                        #           example: if xdir auslenkung = 1.271/0/0 we want to substract the force form quer from (0.9192343/0.9192343/0) (make module longvec_quer_par_force)
                        #           therefore x = (1.271/0/0) would be the maximum displacement
                        #           therefore: go thround every xdir-displacement and get corresponding quer force
                        #           |xdir| = 1.271  --> |quer| = 1.3
                        #           |xdir| = 2.065  --> |quer| = 2.92 (==nndist)
                        #           |xdir| = 0.5    --> |quer| = 0.40225220472457851
                        #           --> |xdir| and |quer| are very similar but for large vectors
                        #           !!! IN PRINCIPLE ONE COULD ALTHO THINK TO DEFINE THE FORCE PARAMETRIZED as a function of |xdir| or |quer|; make those equal and project those on par/quer?
                        #           !!! NO THIS DOES NOT MAKE SENSE; Just think of quer auslenkung! There is no parallel part!
                        # i) once we do have the x and y coordinate of the parkraft: transform h) (parkraft) to parbasis (needs coordinate transformation) (once on quer basis once on par basis)
                        #
                        # ?) problem?: the parkraft will only be defined to 0.6 but it should be in principle be defined as far as |quer|. Here there has to be a certain symmetry!? At the parts were we dont have the quer direction but just the pardirection.
                        # parallel part of ti woks full in xdir; going from xdir to quer or par direction reduces forces of par
                        # !! Wenn wir in einer naeherung bleiben wollen welche aehnlich der Harmonischen naeherung ist (ohne winkelabhaengigkeit) dann koennen wir "nur" eine ti-mode
                        # samplen welche senkrecht zur longitudinalen mode ist


                        #
                        # [351]utils.anyvec_to_inplane_senkr_parallel(np.array([(-2.065+1.271),  2.065+1.271,  0.0   ]),np.array([2.065, -2.065,  0.]),np.array([2.0,2.0,0.0]))
                        # Out[351]: (array([ 1.271,  1.271,  0.   ]), array([ 0.,  0.,  0.]))
                        #
                        # [352]utils.anyvec_to_inplane_senkr_parallel(np.array([(-2.065+1.271),  2.065+0.0,  0.0   ]),np.array([2.065, -2.065,  0.]),np.array([2.0,2.0,0.0]))
                        # Out[352]: (array([ 0.918019,  0.918019,  0.      ]), array([ 0.6355, -0.6355,  0.    ]))
                        #
                        # [350]utils.anyvec_to_inplane_senkr_parallel(np.array([(-2.065+1.271),  2.065-1.271,  0.0   ]),np.array([2.065, -2.065,  0.]),np.array([2.0,2.0,0.0]))
                        # Out[350]: (array([ 0.,  0.,  0.]), array([ 1.271, -1.271,  0.   ]))
                        #
                        #
                        # WHEN this is done, the inplane forces can then be calculated
                        #   a) check the quer auslenkungen where there is no parkraft
                        #   b) chekc the xdisp auslenkung where we have senkkraft and parkraft, all the x auslenkungen should now be perfect (for the 1NN).
                        #   c) the quer auslenkungen: check that symmetries are correctly implemented; how well forces are described ..... we will see; probably/hopefully better then without ti forces

                        ############################
                        # a1)
                        ############################
                        tidataxdirxy = np.copy(self.funcall)
                        #print "tidataxy:"
                        #print tidataxy

                        ############################
                        # a2)
                        ############################
                        print "zeroat:",zeroat
                        ti_xdir_coefsrealene,ti_xdir_coefsrealforces,deltas = polyfit(
                                foldername=False,
                                filename=False,
                                zeroat=zeroat,
                                data=tidataxdirxy,
                                verbose2=False)
                        print "ti_xdir_coefsrealene,ti_xdir_coefsrealforces:",ti_xdir_coefsrealene,ti_xdir_coefsrealforces

                        ############################
                        # b1)
                        ############################
                        print "longreffolder:",longreffolder
                        if contr == 'tix':
                            #tidataquerxy = np.loadtxt(longreffolder+self.fitfolder+'/tix_1nn_all')
                            tidataquerxy = np.loadtxt(longreffolder+'/nnforces/tix_1nn_all')
                        if contr == 'tiy':
                            #tidataquerxy = np.loadtxt(longreffolder+self.fitfolder+'/tiy_1nn_all')
                            tidataquerxy = np.loadtxt(longreffolder+'/nnforces/tiy_1nn_all')

                        ############################
                        # b2)
                        ############################
                        print "zeroat:",zeroat
                        ti_quer_coefsrealene,ti_quer_coefsrealforces,deltas = polyfit(
                                foldername=False,
                                filename=False,
                                zeroat=zeroat,
                                data=tidataquerxy,
                                verbose2=False)
                        print "ti_quer_coefsrealene,ti_quer_coefsrealforces:",ti_quer_coefsrealene,ti_quer_coefsrealforces

                        ############################
                        # c) see utils
                        ############################

                        ############################
                        # f) search for: quermagnitude_coordstransformed_y
                        ############################


                        ############################
                        # g) mesh
                        ############################
                        meshmin = tidataxdirxy[:,0].min()
                        meshmax = tidataxdirxy[:,0].max()
                        if tidataquerxy[:,0].min() > meshmin:
                            meshmin = tidataquerxy[:,0].min()
                        if tidataquerxy[:,0].max() < meshmax:
                            meshmax = tidataquerxy[:,0].max()

                        print "min:",meshmin
                        print "max:",meshmax


                        ################################################################################################################################
                        # h) Now we want to get the difference between the xdisplacement and the querdisplacemetn to determine the ti parallel part
                        ################################################################################################################################
                        # go thround every xdir-displacement and get corresponding quer force
                        def longvec_quer_par_force():
                            pass


                        print ""
                        print "for x-dir auslenkung"

                        tiparxy= np.zeros((len(self.dvec),2))
                        for dvecidx,i in enumerate(self.dvec):
                            vec0 = np.array([ 2.065, -2.065,  0.   ])
                            vecs = np.array([1.0,1.0,0.0])
                            longvec = np.array([-2.065+i[0],2.065,0.0])
                            senkr, parallel = utils.anyvec_to_inplane_senkr_parallel(longvec,vec0,vecs)

                            vec0par = np.array([ 2.065, -2.065,  0.   ])
                            vecspar = np.array([1.0,1.0,0.0])
                            longvecpar = np.array([-2.065+i[0],2.065+i[0],0.0])
                            senkrpar, parallelpar = utils.anyvec_to_inplane_senkr_parallel(longvecpar,vec0par,vecspar)

                            # dies wird bei 'tix' fuer die x komponente und bei und bei 'tiy' fuer die y komponente gemacht
                            # f ist nun die kraft aus quer
                            equer,fquer = pot_energy.getef(senkr[0],ti_quer_coefsrealene)
                            exdir,fxdir = pot_energy.getef(i[0],ti_xdir_coefsrealene)
                            #print "i:",i,"longvec:",longvec,"senkr:",senkr,"parallel:",parallel
                            #print "i[0]:",i[0],"senkr[0]:",senkr[0],"fquer",fquer,"i[0]:",i[0],"fxir:",fxdir
                            #[ 1.3  0.   0. ] [-0.765  2.065  0.   ] senkr: [ 0.948587  0.948587  0.      ] parallel: [ 0.65 -0.65  0.  ] senkr[0]: 0.948586572438 querfunktion(senkr[0]) 0.00992614445689
                            # the corresponding force we now need is force of the x displacment
                            dpar = parallel[0]
                            fpar =  fxdir - fquer
                            tiparxy[dvecidx][0] = dpar
                            tiparxy[dvecidx][1] = fpar
                            #print "this should only be plotted up to 1.3, not to 1.6 (this would extrapolate) the quer regin too far which is only defined up to 0.919 \
                            #        (which is 1.3 absolute in quer direction)"
                            #print "senkr[0]",senkr[0],"meshmax:",meshmax,"meshmin:",meshmin
                            if senkr[0] > meshmax+0.05:
                                tiparxy[dvecidx][0] = 0.0
                                tiparxy[dvecidx][1] = 0.0
                            if senkr[0] < meshmin-0.05:
                                tiparxy[dvecidx][0] = 0.0
                                tiparxy[dvecidx][1] = 0.0

                        tiparxy = np.concatenate((tiparxy,np.array([[zeroat,0.0]])))
                        tiparxy = tiparxy[tiparxy[:,0].argsort()]
                        np.savetxt(self.fitfolder+contr+"_"+shell+"nn_all_parallel_f_in_xy",tiparxy,fmt="%.6f")


                        if contr == 'tiy': #ensured that tix does already exist
                            tix = np.loadtxt(self.fitfolder+"tix"+"_"+shell+"nn_all_parallel_f_in_xy")
                            tiy = np.loadtxt(self.fitfolder+"tiy"+"_"+shell+"nn_all_parallel_f_in_xy")
                            #print "tix:"
                            #print tix
                            #print "tiy:"
                            #print tiy
                            out = np.zeros((tix.shape[0],3))
                            #print "coordinate transformation for plotting quer function"
                            for outidx,outline in enumerate(out):
                                out[outidx][0] = tix[outidx][0]
                                #print  tix[outidx],tiy[outidx],utils.coord_transform_to_quer(np.array([tix[outidx][1],tiy[outidx][1],0.0]))
                                coordtransformed = utils.coord_transform_to_quer(np.array([tix[outidx][1],tiy[outidx][1],0.0]))
                                print "tix:",tix[outidx],"coordtransformed:",coordtransformed
                                out[outidx][1] = coordtransformed[0]
                                out[outidx][2] = coordtransformed[1]
                            quermagnitude_coordstransformed_x = out[:,:2]
                            quermagnitude_coordstransformed_y = out[:,[0,2]]

                            np.savetxt(self.fitfolder+'tix'+"_"+shell+"nn_all_parallel_f_coordstransformed_x",quermagnitude_coordstransformed_x,fmt="%.6f")
                            np.savetxt(self.fitfolder+'tiy'+"_"+shell+"nn_all_parallel_f_coordstransformed_y",quermagnitude_coordstransformed_y,fmt="%.6f")
                            # let's now parametrize this functions
                            polyfit(
                                foldername=self.fitfolder,
                                filename='tix'+"_"+shell+"nn_all_parallel_f_coordstransformed_x_poly",
                                        zeroat=0.0,
                                        data=quermagnitude_coordstransformed_x)
                            polyfit(
                                foldername=self.fitfolder,
                                filename='tiy'+"_"+shell+"nn_all_parallel_f_coordstransformed_y_poly",
                                        zeroat=0.0,
                                        data=quermagnitude_coordstransformed_y)




                #print self.funcallmin
                #print self.funcallmax
                #print "0.2:",self.nnxdist-0.2*self.nnxdist,self.nnxdist,self.nnxdist+0.2*self.nnxdist
                #print "0.3:",self.nnxdist-0.3*self.nnxdist,self.nnxdist,self.nnxdist+0.3*self.nnxdist
                # b) save the fit
                # c) in xdir get corresponding fit in quer auslenkung
                # d) calculate tox by substraction of (quer) long part from tox atoms
                # e) calculate tix by substraction of (quer) long part from tix atoms

        return
