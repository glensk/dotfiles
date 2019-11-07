#!/usr/bin/env python

# for starting something use pot_info_startjob.py !!! not pot_parametrize.

#########################################################################################
# see pot.ANMERKUNG for information regarding this and the other pot scripts
#########################################################################################



#########################################################################################
# all following packages are included in anaconda
#########################################################################################
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
import inspect


from scipy                import optimize
from scipy.interpolate    import UnivariateSpline
from scipy.ndimage        import map_coordinates

import matplotlib.pyplot as plt


import hesse as pot_energy
#reload(pot_energy)
import utils_rename as utils
import my_atom
import imp
imp.reload(my_atom)
import crystal_generator
imp.reload(crystal_generator)
#import pot_energy_forces
#reload(pot_energy_forces)

#np.set_printoptions(threshold=np.nan)  # print the whole array
np.set_printoptions(linewidth=240)    # print only 6 digist after .
np.set_printoptions(precision=3)    # print only 6 digist after .

pd.set_option('display.max_columns', 500)
pd.set_option('display.width', 200)
pd.set_option('display.max_columns', 25)
pd.set_option('display.max_colwidth', 6)

def get_new_local_basis(longvec,longvec0,outdirection0):
    ''' longvec0 und outdirection0 spannnen die basis auf
        longvec ist der atuelle longvec
        outdirection ist die richtng in welche die kraft wirkt/projiziert wird
    '''
    l  = longvec
    l0 = longvec0
    t1 = outdirection0
    t2 = np.cross(l0,t1)

    t11_unnorm = np.cross(l,t2)
    t22_unnorm = np.cross(l,t11_unnorm)

    # die * -1. braucht man fuer die to sachen damit die vorzeichen der vektoren stimmen
    if np.linalg.norm(t11_unnorm) == 0.0:
        t11 = np.array([0.0,0.0,0.0])
    else:
        t11 = t11_unnorm/np.linalg.norm(t11_unnorm) * -1.

    if np.linalg.norm(t22_unnorm) == 0.0:
        t22 = np.array([0.0,0.0,0.0])
    else:
        t22 = t22_unnorm/np.linalg.norm(t22_unnorm) * -1.

    return t11,t22


class fit_to_func( object ):
    def __init__(self, data = False, function = False, fixzeroat = False, weights = False ):
        ''' data can be either a 2d numpy array or a string containing the path to the file to import;
            What is the difference to get_fit_forces_to_pot?
                - fit_to_func is apparently the "NEW" way and get_fit_forces_to_pot the "OLD" way
                - get_fit_forces_to_pot is not used anymore!
        '''
        self.verbose    = False
        self.data       = data

        self.fixzeroat  = fixzeroat # nearest neighbor distance
        self.function   = function
        self.function_known = [ 'l', 'm', 'morse', 'mc1', 'poly', 'ma' ]
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
        self.fitx_interpolated = False
        self.fity_interpolated = False

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
        if self.function == 'l':
            self.parameters, self.function_covariance = \
                optimize.curve_fit(lambda r, eps: pot_energy.LJ_derivative(r, eps, self.fixzeroat), self.datax, self.datay, maxfev=1000)
            self.fity = pot_energy.LJ_derivative(self.fitx, 0.0342863, self.fixzeroat)

        if self.function == 'm' or self.function == 'morse':
            self.model = pot_energy.Morse_derivative
            #print('morse roxx')
            self.parameters, self.function_covariance = \
                    optimize.curve_fit(lambda r, De, aa: pot_energy.Morse_derivative(r, De, aa, self.fixzeroat), self.datax, self.datay, p0=[0.02, 1.96], maxfev=100000)
                #optimize.curve_fit(lambda r, De, aa: pot_energy.Morse_derivative(r, De, aa, self.fixzeroat), datax, datay,p0=[0.02, 1.96],sigma = weights[:,1], absolute_sigma = False, maxfev=10000000)
                #optimize.curve_fit(Morse_derivativeNN, data[:,0], data[:,1], maxfev=1000)
            self.parameters = np.array([self.parameters[0],self.parameters[1],self.fixzeroat])
            self.fity = pot_energy.Morse_derivative(self.fitx, *self.parameters)
            fitdelta = pot_energy.Morse_derivative(self.datax, *self.parameters) - self.datay
            self.fitx_interpolated = np.arange(self.fitx.min(),self.fitx.max(),0.001)
            self.fity_interpolated = pot_energy.Morse_derivative(self.fitx_interpolated, *self.parameters)

        if self.function == 'ma':
            self.model = pot_energy.Morse_repulsive_derivative_to_normaldata
            #print('morse roxx')
            self.parameters, self.function_covariance = \
                    optimize.curve_fit(lambda r, De, aa: pot_energy.Morse_repulsive_derivative_to_normaldata(r, De, aa, self.fixzeroat), self.datax, self.datay, p0=[0.02, 1.96], maxfev=100000)
                #optimize.curve_fit(lambda r, De, aa: pot_energy.Morse_derivative(r, De, aa, self.fixzeroat), datax, datay,p0=[0.02, 1.96],sigma = weights[:,1], absolute_sigma = False, maxfev=10000000)
                #optimize.curve_fit(Morse_derivativeNN, data[:,0], data[:,1], maxfev=1000)
            self.parameters = np.array([self.parameters[0],self.parameters[1],self.fixzeroat])
            self.fity = self.model(self.fitx, *self.parameters)
            fitdelta = self.model(self.datax, *self.parameters) - self.datay
            self.fitx_interpolated = np.arange(self.fitx.min(),self.fitx.max(),0.001)
            self.fity_interpolated = self.model(self.fitx_interpolated, *self.parameters)



        if self.function == 'mc1':
            self.model = pot_energy.mc1_derivative
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
        self.fit_interpolated  = np.transpose([self.fitx_interpolated,   self.fity_interpolated    ])
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
                    print("order    :",order,type(int(order)))
                    print("self.fixzeroat   :",self.fixzeroat,type(float(self.fixzeroat)))
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
                    print("coefs:",coefs)
                er = np.poly1d(self.coefsrealene[::-1])(float(self.fixzeroat))
                fr = np.polyder(np.poly1d(self.coefsrealene[::-1]))(float(self.fixzeroat))
                function = np.polyder(np.poly1d(self.coefsrealene[::-1]))
                if self.verbose:
                    print("function:",function)

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
                    print("order:",order, '%f' % self.fitdelta.max(),"f:",f,"e:",e,"er:",er,"fr:",fr) #,"coefs:",coefs
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
            print("(order)      (dforce)  (er - e)  (fr - f)")
            np.set_printoptions(linewidth=240)    # print only 6 digist after .
            print(self.ordermaxdelta)
            print(self.ordermaxdelta[:,1].argmin())
            print(self.ordermaxdelta.argmin())
        #orderstry = np.arange(1,ordermax+1)  # 7 goes to only to order 6
        ###############################################################################################################################
        # now get the best fitted function (np.arange(ordermaxdelta[:,1].argmin()+1,ordermaxdelta[:,1].argmin()+2))
        ###############################################################################################################################
        #self.ordermaxdelta,self.fitdelta,self.coefsstringene,self.coefsrealene,self.coefsstringforces,self.coefsrealforces,order,self.fitx,self.fity = findbestpolyfit(self.datax,self.datay,self.fixzeroat,
        #        np.arange(self.ordermaxdelta[:,1].argmin()+1,seff.poly_ordermaxdelta[:,1].argmin()+2),
        #        self.dataxminfit,dataxmaxfit)
        #print("findbestpolyfit(    orderstry = np.arange(self.ordermaxdelta[:,1].argmin()+1,self.ordermaxdelta[:,1].argmin()+2))")
        findbestpolyfit(    orderstry = np.arange(self.ordermaxdelta[:,1].argmin()+1,self.ordermaxdelta[:,1].argmin()+2))
        #print("done")
        #                    dataxminfit = self.dataxminfit,
        #                    dataxmaxfit = dataxmaxfit)


        if self.verbose:
            print(self.ordermaxdelta)
            print("-------kkk")
        self.deltas = np.transpose([self.datax, - self.fitdelta])
        self.deltasshifted = np.transpose([self.datax-self.fixzeroat, - self.fitdelta])

        ##############################################################################################################################
        # Save file to disk
        ##############################################################################################################################
        if False: # part for saving files
            ## FILENAME is not defined here!
            print("type(filename)",type(filename),"FILENAME is not defined here!")
            sys.exit()
            if type(filename) != bool:
                deltasfolder = foldername+"/deltas/"
                deltasfoldershifted = foldername+"/deltasshifted/"
                if len(deltasfolder.split("/")) > 16 or len(deltasfolder.split("/")) < 8:
                    print("in polyfit: len:",len(deltasfolder.split("/")))
                    print("in polyfit: foldername:",foldername)
                    print("in polyfit: deltasfolder:",deltasfolder)
                    sys.exit("in polyfit: deltasfolder is wrong")

                deltasfolderfilename = deltasfolder+filenameorig+"_poly"+str(int(ordermax))+"_"+str(int(order))+"_delta___"+self.coefsstringene
                deltasfolderfilenameshifted = deltasfoldershifted+filenameorig+"_poly"+str(int(ordermax))+"_"+str(int(order))+"_delta___"+self.coefsstringene
                if len(deltasfolderfilename.split("/")) > 16:
                    print("in polyfit: deltasfolderfilename:",deltasfolderfilename)
                    sys.exit("in polyfit: deltasfolder is wrong")

                print("name:",filename+"_poly"+str(int(ordermax))+"_"+str(int(order))+"_fit___"+self.coefsstringene)
                print("deltasfolder:",deltasfolder)
                print("deltasfoldershifted:",deltasfoldershifted)
                sys.exit()
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


#def get_fit_forces_to_pot(filename = False, NN = False, pot = False, data = False, weights = False, foldername = "", verbose = True):
#    ''' fits forces to morse, lj, mc1; what is the difference to fit_to_func???
#           --> this function is not used anymore, use fit_to_func class
#           --> this is a simple function, not a class'''
#
#    filenameorig = copy.copy(filename)
#    filename = foldername+filename
#    if type(data) == bool:
#        if type(filename) == bool:
#            sys.exit("please provide a filename")
#        if os.path.isfile(filename) != True:
#            sys.exit("file "+filename+" does not exist")
#        data=np.loadtxt(filename)
#    else:
#        data = data
#    #alat =  float(filename.split("_")[-1])
#    #NN = alat/np.sqrt(2.)
#    if type(NN) == bool:
#        if pot == '135' or pot == '1357' or pot == '13579':
#            pass
#        else:
#            sys.exit("please provide NN")
#    if type(NN) == str:
#        NN= float(NN)
#    if type(pot) == False:
#        sys.exit("please define a pot")
#
#    datax=data[:,0]
#    datay=data[:,1]
#    print "get_fit_forces_to_pot datax:",datax
#
#    if type(weights) != bool:
#        # bring weights on same grid as datax, datay
#        weights = utils.return_function2_on_grid_of_f1(data,weights)
#
#    from scipy import optimize
#    if pot == 'l':
#        parameters, function_covariance = \
#            optimize.curve_fit(lambda r, eps: pot_energy.LJ_derivative(r, eps, NN), datax, datay, maxfev=1000)
#    if pot == 'm' or pot == 'morse':
#        #print 'datx:'
#        #print datax.shape
#        #print datay.shape
#        parameters, function_covariance = \
#                optimize.curve_fit(lambda r, De, aa: pot_energy.Morse_derivative(r, De, aa, NN), datax, datay, p0=[0.02, 1.96], maxfev=100000)
#            #optimize.curve_fit(lambda r, De, aa: pot_energy.Morse_derivative(r, De, aa, NN), datax, datay,p0=[0.02, 1.96],sigma = weights[:,1], absolute_sigma = False, maxfev=10000000)
#            #optimize.curve_fit(Morse_derivativeNN, data[:,0], data[:,1], maxfev=1000)
#        parameters = np.array([parameters[0],parameters[1],NN])
#    if pot == 'mc1':
#        try:
#            parameters, function_covariance = \
#                optimize.curve_fit(lambda r, De, aa, A, B: pot_energy.mc1_derivative(r, De, aa, NN, A, B), datax, datay, maxfev=100000)
#            #optimize.curve_fit(lambda r, De, aa, A, B: pot_energy.mc1_derivative(r, De, aa, NN, A, B), datax, datay, sigma = weights[:,1], absolute_sigma = False, maxfev=1000000)
#        except RuntimeError:
#            return False, False, False, False, False, False
#        parameters = np.array([parameters[0],parameters[1],NN,parameters[2],parameters[3]])
#
#    fitx = np.linspace(datax.min()-1.0,datax.max()+1.0, num=100)
#    if pot == 'l':
#        #fity = LJ_derivative(fitx, parameters[0], NN)
#        fity = pot_energy.LJ_derivative(fitx, 0.0342863, NN)
#    if pot == 'm' or pot == 'morse':
#        fity = pot_energy.Morse_derivative(fitx, parameters[0], parameters[1], NN)
#        fitdelta = pot_energy.Morse_derivative(datax, parameters[0], parameters[1], NN) - datay
#        #print "fd:",fitdelta
#    if pot == 'mc1':
#        fity = pot_energy.mc1_derivative(fitx, parameters[0], parameters[1], NN, \
#                parameters[2], parameters[3])
#        fitdelta = pot_energy.mc1_derivative(datax, parameters[0], parameters[1], NN, \
#                parameters[2], parameters[3]) - datay
#
#
#    #print "saving..."
#    deltas = np.transpose([datax, -fitdelta])
#    deltasshifted = np.transpose([datax-NN, -fitdelta])
#
#    deltasfolder = foldername+"/deltas/"
#    deltasfoldershifted = foldername+"/deltasshifted/"
#
#    datashiftedfolder = foldername+"/datashifted"
#    datashifted = np.transpose([datax-NN, datay])
#
#    deltasfolderfit = foldername+"/deltasfit/"
#    deltasfolderfitshifted = foldername+"/deltasfitshifted/"
#
#    if len(deltasfolder.split("/")) > 14:
#        print "in get_fit_forces_to_pot: deltasfolder in :",deltasfolder,len(deltasfolder)
#        sys.exit("in get_fit_forces_to_pot: deltasfolder is wrong")
#    if os.path.isdir(deltasfolder) != True: os.makedirs(deltasfolder)
#    if os.path.isdir(deltasfoldershifted) != True: os.makedirs(deltasfoldershifted)
#    if os.path.isdir(deltasfolderfit) != True: os.makedirs(deltasfolderfit)
#    if os.path.isdir(datashiftedfolder) != True: os.makedirs(datashiftedfolder)
#
#    deltasfolderfilename = deltasfolder+filenameorig+"_fit_delta___"+str(parameters[0])[:8]
#    if len(deltasfolderfilename.split("/")) > 14:
#        print "in get_fit_forces_to_pot: deltasfolderfilename:",deltasfolderfilename
#        sys.exit("in get_fit_forces_to_pot: deltasfolder is wrong")
#
#    polycoefsrealene,polycoefsrealforces,polydeltas = polyfit(
#    foldername=deltasfolderfit,
#    filename=filenameorig+"_add_",
#            zeroat=NN,
#            data=deltas,
#            verbose2=False)
#
#    polyfit(
#    foldername=deltasfolderfitshifted,
#    filename=filenameorig+"_add_",
#            zeroat=0.0,
#            data=deltasshifted,
#            verbose2=False)
#
#    if pot == 'l':
#        np.savetxt(
#            filename+"_fit_lj_"+str(parameters[0])[:8]
#            +"_"+str(NN)[:8],
#                np.transpose([fitx, fity]),
#                fmt='%.18f',
#                delimiter='  ')   # X is an array
#
#    if pot == 'm' or pot == 'morse':
#        if type(filename) != bool:
#            np.savetxt(
#                filename+"_fit___"+str(parameters[0])[:8]
#                +"_"+str(parameters[1])[:8]
#                +"_"+str(NN)[:8],
#                    np.transpose([fitx, fity]),
#                    fmt='%.18f',
#                    delimiter='  ')   # X is an array
#            np.savetxt(
#                datashiftedfolder+filenameorig+"_fit___"+str(parameters[0])[:8]
#                +"_"+str(parameters[1])[:8]
#                +"_"+str(NN)[:8],
#                    datashifted,
#                    fmt='%.18f',
#                    delimiter='  ')   # X is an array
#
#            np.savetxt(
#                deltasfolder+filenameorig+"_fit_delta___"+str(parameters[0])[:8]
#                +"_"+str(parameters[1])[:8]
#                +"_"+str(NN)[:8],
#                    deltas,
#                    fmt='%.18f',
#                    delimiter='  ')   # X is an array
#            np.savetxt(
#                deltasfoldershifted+filenameorig+"_fit_delta___"+str(parameters[0])[:8]
#                +"_"+str(parameters[1])[:8]
#                +"_"+str(NN)[:8],
#                    deltasshifted,
#                    fmt='%.18f',
#                    delimiter='  ')   # X is an array
#
#    if pot == 'mc1':
#        #print "parss mc1:",parameters
#        #parameters = np.array([parameters[0],parameters[1],NN,parameters[2],parameters[3]])
#        #print "parss mc1:",parameters
#        #print NN
#        #sys.exit()
#        if type(filename) != bool:
#            np.savetxt(
#                filename+"_fit___"+str(parameters[0])[:8]
#                +"_"+str(parameters[1])[:8]
#                +"_"+str(NN)[:8]
#                +"_"+str(parameters[2])[:8]
#                +"_"+str(parameters[3])[:8],
#                    np.transpose([fitx, fity]),
#                    fmt='%.18f',
#                    delimiter='  ')   # X is an array
#            #print "filenameorig:",filenameorig
#
#            np.savetxt(
#                deltasfolder+filenameorig+"_fit_delta___"+str(parameters[0])[:8]
#                +"_"+str(parameters[1])[:8]
#                +"_"+str(NN)[:8]
#                +"_"+str(parameters[2])[:8]
#                +"_"+str(parameters[3])[:8],
#                    deltas,
#                    fmt='%.18f',
#                    delimiter='  ')   # X is an array
#            np.savetxt(
#                deltasfoldershifted+filenameorig+"_fit_delta___"+str(parameters[0])[:8]
#                +"_"+str(parameters[1])[:8]
#                +"_"+str(NN)[:8]
#                +"_"+str(parameters[2])[:8]
#                +"_"+str(parameters[3])[:8],
#                    deltasshifted,
#                    fmt='%.18f',
#                    delimiter='  ')   # X is an array
#    return parameters, deltas, deltasshifted, polycoefsrealene,polycoefsrealforces,polydeltas
#
def polyfit_with_fixed_points(n, x, y, xf, yf) :
    if False:
        print("n:",n)
        print("x:",x)
        print("y:",y)
        print("xf:",xf,type(xf))
        print("yf:",yf,type(yf))
        print(n,len(xf),type(n),type(len(xf)))
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
    #print("mat:",mat)
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

    print("polyfit datax:",datax)
    #print "#################: datax:",zeroat
    #print datax
    #print "#################: datay:",zeroat
    #print datay

    if ordermax == False:
        ordermax = datax.shape[0]/2.
        #if ordermax > 9:
        #    ordermax = 9
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
            #print "verbose:",verbose

            if verbose:
                print("order    :",order,type(int(order)))
                print("zeroat   :",zeroat,type(float(zeroat)))
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
                print("coefs:",coefs)
            er = np.poly1d(coefsrealene[::-1])(float(zeroat))
            fr = np.polyder(np.poly1d(coefsrealene[::-1]))(float(zeroat))
            function = np.polyder(np.poly1d(coefsrealene[::-1]))
            if verbose:
                print("function:",function)

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
                print("order:",order, '%f' % fitdeltamax,"f:",f,"e:",e,"er:",er,"fr:",fr) #,"coefs:",coefs
            #print "coefs:",coefs
            #print "coefse:",coefse
            #print "coefse:",coefse
        return ordermaxdelta,fitdelta,coefsstringene,coefsrealene,coefsstringforces,coefsrealforces,order,fitx,fity




    #print "orderstry:",orderstry
    ###############################################################################################################################
    # try fit with several orders 1 ... 9  (orderstry)
    ###############################################################################################################################
    print("datax:",datax)
    print("datay:",datay)
    ordermaxdelta,fitdelta,coefsstringene,coefsrealene,coefsstringforces,coefsrealforces,order,fitx,fity = findbestpolyfit(datax,datay,zeroat,orderstry,dataxminfit,dataxmaxfit)
    if verbose2:
        print("(order)      (dforce)  (er - e)  (fr - f)")
        np.set_printoptions(linewidth=240)    # print only 6 digist after .
        print(ordermaxdelta)
        print("xxxxxxyy")
        print(ordermaxdelta[:,1].argmin())
        print(ordermaxdelta.argmin())
        print("-------")
    #orderstry = np.arange(1,ordermax+1)  # 7 goes to only to order 6
    ###############################################################################################################################
    # now get the best fitted function (np.arange(ordermaxdelta[:,1].argmin()+1,ordermaxdelta[:,1].argmin()+2))
    ###############################################################################################################################
    ordermaxdelta,fitdelta,coefsstringene,coefsrealene,coefsstringforces,coefsrealforces,order,fitx,fity = findbestpolyfit(datax,datay,zeroat,
            np.arange(ordermaxdelta[:,1].argmin()+1,ordermaxdelta[:,1].argmin()+2),
            dataxminfit,dataxmaxfit)
    if verbose2:
        print("ordermaxdelta;",ordermaxdelta)
        print("-------")
    deltas = np.transpose([datax, -fitdelta])
    deltasshifted = np.transpose([datax-zeroat, -fitdelta])
    deltamax = fitdelta


    #print "type(filename):",type(filename)

    if type(filename) != bool:
        deltasfolder = foldername+"/deltas/"
        deltasfoldershifted = foldername+"/deltasshifted/"
        if len(deltasfolder.split("/")) > 16 or len(deltasfolder.split("/")) < 8:
            print("in polyfit: len:",len(deltasfolder.split("/")))
            print("in polyfit: foldername:",foldername)
            print("in polyfit: deltasfolder:",deltasfolder)
            sys.exit("in polyfit: deltasfolder is wrong")

        deltasfolderfilename = deltasfolder+filenameorig+"_poly"+str(int(ordermax))+"_"+str(int(order))+"_delta___"+coefsstringene
        deltasfolderfilenameshifted = deltasfoldershifted+filenameorig+"_poly"+str(int(ordermax))+"_"+str(int(order))+"_delta___"+coefsstringene
        if len(deltasfolderfilename.split("/")) > 16:
            print("in polyfit: deltasfolderfilename:",deltasfolderfilename)
            sys.exit("in polyfit: deltasfolder is wrong")

        #print "deltasfolder:",deltasfolder
        #print "deltasfoldershifted:",deltasfoldershifted
        #print "name:",filename+"_poly"+str(int(ordermax))+"_"+str(int(order))+"_fit___"+coefsstringene

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

def dos_to_pot(dos,tmelt):
    ''' gets dos (2d numpy array) and tmelt in Kelvin and returns the potential '''
    pot = np.copy(dos)
    pot[:,1] = -np.log(dos[:,1])*tmelt*0.086173423
    pot[:,1] = pot[:,1]+abs(pot[:,1].min())
    #pot[:,1] = pot_
    return pot




class forcesneighbors( object ):
    '''
    analyzes the forces of the different neighbor shells
    this class is designed for a certain alat which has to be defined beforehand
    f[1]    : all forces of first alat WRONG NOW
    f[1,1]  : all forces of first alat and first displacement WRONG NOW
    f[1,1,1]: all forces fo first alat and first displacement in x,y,z direction WRONG NOW

    for fcc:
    form quer we want to get    long 1NN (  this appears to be pure compared to the
                                            x-direction which contains also an inplane
                                            part)

    from xdir we want to get    long 1NN (this we can get also quer disp might be better)
                                long 2NN
                                tox 1NN
                                tox 2NN

    '''
    def __init__(self):
        self.a                          = False      # alat
        self.disp                       = False      # defines if the current displacement is a disp in 'xdir' or 'quer'
                                                     # quer direction
        self.sc                         = False
        self.folder                     = False      # list of folders
        self.dnorm                      = False      # norm of displacement (from origin) (should be of len a)
        self.dvec                       = False      # vec of displacement (from origin) (should be of len a)
        self.rcar                       = False
        self.f                          = False
        self.p                          = False      # positions (should be of lenght a)
        self.struct                     = False      # fcc, bcc ...
        self.element                    = False
        self.nnforcesaddname            = False      # Fase means to add nothing to nnforces
        self.disp_all                   = [ 'xdir', 'quer', 'midd', 'qubus', '3nnd', '4nnd' ]
        self.pkl                        = False
        self.fit_file                   = False
        self.verbose                    = False
        self.parametrize_nnforcesnpz            = False
        self.parametrize_data_pkl       = False
        self.parametrize_fml_npz        = False
        self.parametrize_fml_npz_path   = False

        self.fml                        = False
        self.fml_u1nn_pottype           = False
        self.fml_u1nn_potparam          = False
        self.fml_u1nn_potadd            = False
        self.fml_u1nn_potaddparam       = False
        self.FML                        = False

        self.fmlmti                     = False
        self.FMLmti                     = False

        self.fmlmto                     = False
        self.FMLmto                     = False

        self.fml_u1nn_tox               = False
        self.fml_u1nn_tix               = False

        self.angles                     = False
        self.fprojlongold               = False

        self.folder_base                = "/Users/glensk/Dropbox/Understanding_distributions/"
        if os.path.isdir(self.folder_base) != True:
            sys.exit("ERROR: folder_base (in pot_parametrize.py) does not exist:"+str(self.folder_base))

        self.initialize_dataframe()

    def initialize_dataframe(self):
        ''' initializes an empty dataframe which will be filled afterwards '''

        self.datacolumns = [ 'shells', 'contr', 'disp' ,'range' ]

        self._contr  = [ 'lon', 'to', 'ti','tox', 'toy','toz', 'tix', 'tiy', 'tiz' ] ## -> dies ist die column coordinate
        self._contr  = [ 'lon', 'to', 'ti' ] ## -> dies ist die column coordinate
        self._diff   = self._contr
        self._shells = [ 1, 2, 3 ]
        self._func   = [ 'morse', 'mc1' , 'poly']
        self._range  = [ 'alle', 'dos', 'doscut', 'dosneg', 'dospos' ] # dont take all, gets confused with pandas
        self._disp   = self.disp_all


        self.datarows= [ 'points', 'force', 'file', 'dos', 'pot', 'pointsshifted', 'dosshifted', 'potshifted' ] + self._func + [ i+'fit' for i in self._func]

        # self.datacolumns = [ 'x', 'y', 'z' ]
        # x = [ -1.3 .... 1.3 ]
        # y = [ -1.3 .... 1.3 ]
        # z = [ -1.3 .... 1.3 ]

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

    def creatennforcesnpz(self, folder = False, a = False,energyauswertung=False):
        ''' loads nnforces.npz from folder if present
            if not present: create it
            energyauswertung creates:
                  - dUdL*
                  - ene_vs_disp*
                  - udftmev_vs_disp
                  - POSITIONs
        '''
        if type(folder) == bool:
            sys.exit("pot_parametrize.forcesneighbors.creatennforcesnpz: please specify the folder where to load data")
        if os.path.isdir(folder) != True:
            print("folder:",folder)
            sys.exit("pot_parametrize.forcesneighbors.creatennforcesnpz: specifies folder does not exist!")
        self.folderparam = folder
        if self.parametrize_nnforcesnpz == False:
            if os.path.isfile(self.folderparam+"/nnforces.npz") == True:
                self.loadnnforcesnpz(self.folderparam,"nnforces")
                return

        self.eqcoords = np.loadtxt(self.folderparam+"/EqCoords_direct")
        #self.eqcoords_cart = pot.disp.alle.eqcoords*pot.disp.alle.a*int(pot.forparametrization_sc)


        self.atoms = self.eqcoords.shape[0]
        if type(self.struct) == bool:
            if self.atoms == 54: self.struct = 'bcc'
            if self.atoms == 32: self.struct = 'fcc'
            if self.atoms == 108: self.struct = 'fcc'
        if type(self.struct) == bool:
            sys.exit("self.struct not known, ABER BEKANNNT: pot.structure = 'fcc'  und pot.forparametrization_sc = '5', wenn sollte das uebergeben werden; number of atoms kann auch schon hier berechnet sein.")

        #self.folderparam = os.getcwd()
        self.element = my_atom.get_element_from_path(self.folderparam)
        self.atom = my_atom.atom([self.element])
        if type(self.a) == bool:
            self.a = a
        if type(self.a) == bool:
            sys.exit("when loading creatennforcesnpz you have to specify \"a\"")

        if self.struct == 'fcc':
            self.scc = round(float((self.atoms/4.)**(1/3.)),0)
            self.sc = np.array([[self.scc*self.a,0.,0.],[0.,self.scc*self.a,0.],[0.,0.,self.scc*self.a]])
            self.eqcoords_cart = self.eqcoords*self.a*self.scc




        self.folder = utils.lsn(self.folderparam+"/"+str(self.a)+"Ang_"+"*")
        #self.folder = utils.lsn(self.folderparam+"/"+str(self.a)+"Ang_"+"*")[:20]
        # hier sollte man nun erstmal durch den folder gehen und schauen welche folder
        # drin bleiben (fertig) und welche nicht (nicht fertig), dann sollte self.folder neu
        # definiert werden
        if len(self.folder) == 0:
            print("self.a",self.a)
            sys.exit("Error in creatennforcesnpz: no folder found with name: "+str(self.a)+"Ang_"+"*")
        #print "||||||||||||||||",len(self.folder)
        #print "||||||||||||||||",self.folderparam
        #print "||||||||||||||||",self.atoms
        ####################################################################################
        ####################################################################################
        # das wird hier geladen
        ####################################################################################
        ####################################################################################
        self.f = np.zeros((np.array(self.folder).shape[0],self.atoms,3)) # forces orig cell
        self.p = np.zeros((np.array(self.folder).shape[0],self.atoms,3)) # positions orig cell
        self.pr = np.zeros((np.array(self.folder).shape[0],self.atoms,3))  # removed mapping
        self.e = np.zeros(np.array(self.folder).shape[0])
        self.dstring = list(range(len(self.folder)))
        self.dstring1 = list(range(len(self.folder)))
        self.dstring2 = list(range(len(self.folder)))
        self.dstring3 = list(range(len(self.folder)))
        self.dvec = np.zeros((np.array(self.folder).shape[0],3))
        self.dvec1 = np.zeros((np.array(self.folder).shape[0],3))
        self.dvec2 = np.zeros((np.array(self.folder).shape[0],3))
        self.dvec3 = np.zeros((np.array(self.folder).shape[0],3))
        self.dnorm = np.zeros((np.array(self.folder).shape[0]))
        self.dnorm1 = np.zeros((np.array(self.folder).shape[0]))
        self.dnorm2 = np.zeros((np.array(self.folder).shape[0]))
        self.dnorm3 = np.zeros((np.array(self.folder).shape[0]))

        nsc = 2  # for now on we will just double the supercell, in general maybe factor
        self.P = np.zeros((self.pr.shape[0],self.pr.shape[1]*nsc**3,3))
        self.F = np.zeros((self.pr.shape[0],self.pr.shape[1]*nsc**3,3))
        self.Pr = np.zeros((self.pr.shape[0],self.pr.shape[1]*nsc**3,3)) #  (centered around 0 atom)
        self.P0 = False
        self.F0 = False
        self.Pr0 = False
        self.positions = False
        self.ene_vs_disp = False
        self.dudl = False
        self.udftmev_vs_disp = False

        if os.path.isfile(self.folderparam+"/nnforces.alldata.npz") == True:
            print("-------->>>",self.folderparam+"/nnforces.alldata.npz EXISTS, therefore I dont need to go through inputfolder")
            self.loadnnforcesnpz(self.folderparam,"nnforces.alldata", parametrize_data_pkl=False)
        else:
            for folderidx,folder in enumerate(self.folder):
                print("___>>>",folderidx,folder)
                #print "dstring:",utils.string_to_num_list(folder)
                # hier haengt es nun davon ab wie viele atome ausgelenkt worden sind.
                self.dstring1[folderidx] = utils.string_to_num_list(folder)[-3]
                self.dstring2[folderidx] = utils.string_to_num_list(folder)[-2]
                self.dstring3[folderidx] = utils.string_to_num_list(folder)[-1]
                self.dstring[folderidx] = utils.string_to_num_list(folder)[-1]
                #print "disp:",self.dstring1[folderidx],self.dstring2[folderidx],self.dstring3[folderidx]
                #self.dvec[folderidx] = self.p[folderidx][0]
                #self.dvec1[folderidx] = self.p[folderidx][0]
                #self.dvec2[folderidx] = self.p[folderidx][375]
                #self.dvec3[folderidx] = self.p[folderidx][495]


                if os.path.isfile(folder+"/forces_OUTCAR") != True:
                    if os.path.isfile(folder+"/forces") != True:
                        if os.path.isfile(folder+"/OUTCAR") == True or os.path.isfile(folder+"/OUTCAR.gz") == True:
                            utils.run2("OUTCAR_forces-last-ARRAY.sh "+folder+" > "+folder+"/forces_OUTCAR")

                # The next 4 lines are important when OUTCAR is not finished! dont remove them
                if os.path.isfile(folder+"/forces_OUTCAR") != True and os.path.isfile(folder+"/forces") != True:
                    try:
                        self.f[folderidx] = np.loadtxt(folder+"/forces_OUTCAR")
                    except ValueError:
                        print((folder+" redo forces_OUTCAR"))
                        sys.exit()
                #print "sf:",self.f[folderidx].shape
                if os.path.isfile(folder+"/forces_OUTCAR") == True:
                    self.f[folderidx] = np.loadtxt(folder+"/forces_OUTCAR")
                elif os.path.isfile(folder+"/forces") == True:
                    self.f[folderidx] = np.loadtxt(folder+"/forces")
                else:
                    sys.exit("could not find forces in "+folder)

                print(str(folderidx)+"/"+str(len(self.folder)-1),folder,"force[0]:",self.dstring[folderidx],self.f[folderidx][0])


                ################################################
                # positions
                ################################################
                pos = False
                if os.path.isfile(folder+"/pos") == True:
                    self.p[folderidx] = np.loadtxt(folder+"/pos")
                    pos = True
                if os.path.isfile(folder+"/cartesian_coords") == True:
                    self.p[folderidx] = np.loadtxt(folder+"/cartesian_coords")
                    pos = True
                elif os.path.isfile(folder+"/cartesian_coords") != True:
                    utils.run2("OUTCAR_positions-last-ARRAY.sh "+folder+" > "+folder+"/cartesian_coords")
                    self.p[folderidx] = np.loadtxt(folder+"/cartesian_coords")
                    pos = True

                if pos == False:
                    sys.exit("could not load positions@!")

                if os.path.isfile(folder+"/ene_free_last") != True:
                    if os.path.isfile(folder+"/OUTCAR") == True or os.path.isfile(folder+"/OUTCAR.gz") == True:
                        utils.run2("OUTCAR_ene-free-last.sh "+folder+" > "+folder+'/ene_free_last')
                        #self.e[folderidx] = np.loadtxt(folder+"/ene_free_last")  # convert to meV per atom
                        self.e[folderidx] = np.loadtxt(folder+"/ene_free_last")*1000./(self.atoms-1)  # convert to meV per atom
                #self.energymevtopy = np.sum(self.energytrantopy)/2.*1000/(self.numberofatoms-1)

                if os.path.isfile(folder+"/cell") != True:
                    if os.path.isfile(folder+"/OUTCAR") == True or os.path.isfile(folder+"/OUTCAR.gz") == True:
                        utils.run2("OUTCAR_cell-last-cartesian-ARRAY.sh "+folder+" > "+folder+"/cell")
                        self.sc = np.loadtxt(folder+"/cell")

                if type(self.sc) == bool and self.struct == 'fcc':
                    self.scc = float((self.atoms/4.)**(1/3.))
                    self.sc = np.array([[scc*self.a,0.,0.],[0.,scc*self.a,0.],[0.,0.,scc*self.a]])

                #print self.p[folderidx]
                #print "sc:|",self.sc
                #print "alat:",self.a
                #print "str:",self.struct
                #print "atoms:",self.atoms
                # run  /Users/glensk/Thermodynamics/python_thermodynamics/pot_info_startjob.py  -e al -a 4.14 -sp 5 -kp 2x2x2kp -dispdirection quer -pnnforcesnpz -v -stradd "_vasp4_ENCUT250_GGA_2_atoms_displaced"
                #

                #print "sc:|",self.sc
                #sys.exit()
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
                # this is only necessary when having 2 displacements
                #print "ka",self.atoms
                if self.atoms > 495:
                    self.dvec1[folderidx] = self.pr[folderidx][0]-self.eqcoords_cart[0]
                    self.dvec2[folderidx] = self.pr[folderidx][375]-self.eqcoords_cart[375]
                    self.dvec3[folderidx] = self.pr[folderidx][495]-self.eqcoords_cart[495]
                    # pot.disp.alle.p[0,0]   --> array([ 0.   ,  0.   ,  0.  ])
                    # pot.disp.alle.p[0,375] --> array([ 2.07 ,  2.07 ,  0.  ])
                    # pot.disp.alle.p[7,495] --> array([ 18.63,  18.63,  0.  ])
                #print "kkkk:",self.pr[folderidx][0]
                self.dnorm[folderidx] = np.linalg.norm(self.pr[folderidx][0])
                if self.atoms > 495:
                    self.dnorm1[folderidx] = round(np.sign(self.dvec1[folderidx][0])*np.linalg.norm(self.dvec1[folderidx]),4)
                    self.dnorm2[folderidx] = round(np.sign(self.dvec2[folderidx][0])*np.linalg.norm(self.dvec2[folderidx]),4)
                    self.dnorm3[folderidx] = round(np.sign(self.dvec3[folderidx][0])*np.linalg.norm(self.dvec3[folderidx]),4)
                    if self.verbose:
                        print("dvec{1,2,3} :",self.dvec1[folderidx],self.dvec2[folderidx],self.dvec3[folderidx])
                        print("dnorm{1,2,3}:",self.dnorm1[folderidx],self.dnorm2[folderidx],self.dnorm3[folderidx])

                    #print "kk::",self.p[folderidx][0]
                    #print "kk::",self.p[folderidx][1]
                    #print "kk::",self.p[folderidx][2]



            ####################################################################################
            # an diseer stelle sollten die daten gesichert werden, nicht spaeter
            # spaeter koennte ein fehler auftreten
            # dis hier ist direkt nach der forschleife ueber alle displacements
            self.savennforcesnpz(self.folderparam,"nnforces.alldata", parametrize_lon=False)
            ####################################################################################
        print("-----------")
        print(self.dnorm)
        print("-----------")
        #sys.exit()
        ###################################################################################
        # get other variables (energy evaluation)
        ###################################################################################
        if energyauswertung:
            allposforces = np.zeros((self.p.shape[0]*self.p.shape[1]+self.atoms,6))  # + self.atoms um die erste (unausgelenkte struktur) hinzuzufuegen
            allpos = self.p.flatten().reshape(self.p.shape[0]*self.p.shape[1],3)
            allforces = self.f.flatten().reshape(self.f.shape[0]*self.f.shape[1],3)
            allposforces[self.atoms:,0:3] = allpos
            allposforces[self.atoms:,3:6] = allforces


            # ene_vs_disp
            self.positions = allposforces.reshape(self.p.shape[0]*self.p.shape[1]+self.atoms,6)
            self.ene_vs_disp = np.array((self.dstring,self.e)).transpose()

            # POSITIONs
            print("kkkkkkkkkkkk",self.ene_vs_disp.shape)  # macht probleme wenn nicht gerechnet die perfekte struktur ohne auslenkung
            print(self.ene_vs_disp)
            minidx = np.nonzero(self.ene_vs_disp[:,0] == 0.0)[0][0]  # dies ist nicht immer berechnet
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
            print("self.udftmev_vs_disp_short:")
            print(self.udftmev_vs_disp_short)

            np.savetxt(self.folderparam+"/ene_vs_disp_short",self.ene_vs_disp_short)
            np.savetxt(self.folderparam+"/POSITIONs_short",self.positions_short,fmt="%.8f")
            np.savetxt(self.folderparam+"/udftmev_vs_disp_short",self.udftmev_vs_disp_short)
            np.savetxt(self.folderparam+"/dUdL_short",self.dudl_short,fmt="%7.3f%10.1f%9.1f%9.1f%14.2f%10.2f%14.2f%10.2f%10.2f",header=" step   time(fs)  temp(K) average       U(meV/at)    Uref          dUdL   average    offset")
        self.disp = self.folderparam.split("/")[-1].split("_")[1]
        if self.disp not in self.disp_all:
            for i in self.disp_all:
                if i in self.folderparam.split("/")[-1].split("_"):
                    self.disp = i
            print(utils.printred("self.disp_all:",self.disp_all))
        print(utils.printred("self.disp:",self.disp))
        #sys.exit()
        # here self.P / self.Pr and self.F are created
        self.repeat_supercell_positions_forces()
        print("skk:",self.verbose)
        self.savennforcesnpz(self.folderparam,"nnforces")
        return

    def get_forces_from_3_displaced_atoms(self,one,two,three):
        ''' one, two, three stand for the displaced atoms;

        one, two, three can be:
            - displacement  like 0.1, -0.6, 0.0 for the corresponding atoms
            - :             like :  ,  0.0, 0.0 for seeing all forces when only first atom is displaced


        to run this:

        run  /Users/glensk/Thermodynamics/python_thermodynamics/pot_info_startjob.py  -e al -a 4.14 -sp 5 -kp 2x2x2kp -dispdirection quer -pnnforcesnpz -stradd "_vasp4_ENCUT250_GGA_2_atoms_displaced"  -v

        pot.disp.alle.get_forces_from_3_displaced_atoms('all',0.0,0.0)
        '''
        #if one not in np.unique(pot.disp.alle.dnorm1):
        if one != "all":
            if one not in np.unique(self.dnorm1): # or one != 'all':
                print("one asked for :",one)
                print("one may be    :",np.unique(self.dnorm1))
                sys.exit("Error: one has to be one of above but is not!")
        if two != "all":
            if two not in np.unique(self.dnorm2): # or two != 'all':
                print("two asked for :",two)
                print("two may be    :",np.unique(self.dnorm2))
                sys.exit("Error: two has to be one of above but is not!")
        if three != "all":
            if three not in np.unique(self.dnorm3): # or three != 'all':
                print("three asked for :",three)
                print("three may be    :",np.unique(self.dnorm3))
                sys.exit("Error: three has to be one of above but is not!")

        #check1 = np.where(pot.disp.alle.dnorm1 == one  )[0]
        #check2 = np.where(pot.disp.alle.dnorm2 == two  )[0]
        #check3 = np.where(pot.disp.alle.dnorm3 == three)[0]
        if one == 'all':
            check1 = np.where(self.dnorm1 < 9999999999999  )[0]
        else:
            check1 = np.where(self.dnorm1 == one  )[0]

        if two == 'all':
            check2 = np.where(self.dnorm2 < 9999999999999  )[0]
        else:
            check2 = np.where(self.dnorm2 == two  )[0]


        if three == 'all':
            check3 = np.where(self.dnorm3 < 9999999999999  )[0]
        else:
            check3 = np.where(self.dnorm3 == three)[0]

        # schnittmenge
        schnittmenge12 = [val for val in check1 if val in check2 ]
        schnittmenge123= [val for val in schnittmenge12 if val in check3 ]
        out=schnittmenge123
        #print "check1:",check1
        #print "check2:",check2
        #print "check3:",check3
        #print "schnittmenge12:",schnittmenge12
        #print "schnittmenge123:",schnittmenge123

        # pot.disp.alle.p[pot.disp.alle.get_forces_from_3_displaced_atoms('all',0.0,0.0),0]
        # pot.disp.alle.f[pot.disp.alle.get_forces_from_3_displaced_atoms('all',0.0,0.0),0]
        # am schoensten waere es wenn man als letzten index (also hier die 0)
        # 1 fuer 1NN angibt
        # -1  fuer den 1NN auf der anderen seite
        # 2 fuer den NN x = [ 4.14 4.14 0 ]
        # - 2 fuer den NN x = [ -4.14 -4.14 0 ]
        # usw.
        #
        #
        # TODO:
        # 1) den normalen morse anzeigen fuer one, two, three = 'all', 0.0, 0.0                                         DONE
        # 2) dann die forces auf atome bei [ 4.14 4.14 0 ] und [ -4.14 -4.14 0 ] abziehen von den forces bei 1          DONE
        # 3) auch interessant ist wie sich die Kraft auf das erste atom aendert wenn sich das zweite bewegt worden ist
        #
        #   pot.disp.alle.p[pot.disp.alle.get_forces_from_3_displaced_atoms('all',0.0,0.0),0]
        #   pot.disp.alle.p[pot.disp.alle.get_forces_from_3_displaced_atoms('all',0.0,0.0),375]
        #   pot.disp.alle.p[pot.disp.alle.get_forces_from_3_displaced_atoms('all',0.0,0.0),495]
        #
        #   np.linalg.norm(pot.disp.alle.f[pot.disp.alle.get_forces_from_3_displaced_atoms('all',0.0,0.0),375],axis=1)
        #   np.linalg.norm(pot.disp.alle.f[pot.disp.alle.get_forces_from_3_displaced_atoms('all',0.0,0.0),495],axis=1)
        #
        # FORCES ON 1NN:
        # frep = np.linalg.norm(pot.disp.alle.f[pot.disp.alle.get_forces_from_3_displaced_atoms('all',0.0,0.0),375],axis=1)
        # fatt = np.linalg.norm(pot.disp.alle.f[pot.disp.alle.get_forces_from_3_displaced_atoms('all',0.0,0.0),495],axis=1)
        # DISP ON 1NN:
        #  pot.disp.alle.dnorm[pot.disp.alle.get_forces_from_3_displaced_atoms('all',0.0,0.0)]
        #  pot.disp.alles.sc        ==  [ 20.7 ....]

        #  pot.disp.alle.p[320,375] ==  [ 2.07,  2.07,  0.   ] rep  (vasp line 376)
        #  pot.disp.alle.p[1,495] ==    [ 18.63,  18.63,   0.   ] attr (vasp line 496)

        #  pot.disp.alle.p[1,30] ==     [ 4.14,  4.14,  0.  ] rep
        #  pot.disp.alle.p[320,120] ==  [ 16.56,  16.56,   0.  ] attr
        #
        #  pot.disp.alle.p[320,405] ==  [ 6.21,  6.21,  0.  ] rep (vasp line 406)
        #  pot.disp.alle.p[320,465] ==  [ 14.49,  14.49,   0.  ] attr (vasp line 466)
        #
        #  pot.disp.alle.p[320,60]==    [ 8.28,  8.28,  0.  ] rep (vasp line 61)
        #  pot.disp.alle.p[320,90] ==   [ 12.42,  12.42,   0.  ] attr (vasp line 91)
        #
        #  out sind einfach die kraefte fuer die displacements die wir wollen ('all, 0.0, 0.0)
        #  und 357 oder 495 oder 30 oder wie auch immer ist das atom dass wir anschauen
        frep = np.linalg.norm(self.f[out,375],axis=1)
        fatt = np.linalg.norm(self.f[out,495],axis=1) # vasp 496

        frep2 = np.linalg.norm(self.f[out,30],axis=1)
        fatt2 = np.linalg.norm(self.f[out,120],axis=1) # vasp 121

        frep3 = np.linalg.norm(self.f[out,405],axis=1)
        fatt3 = np.linalg.norm(self.f[out,465],axis=1)

        frep4 = np.linalg.norm(self.f[out,60],axis=1)
        fatt4 = np.linalg.norm(self.f[out,90],axis=1)

        dd = self.dnorm[out]

        print("frep:",frep,frep.shape)
        print("fatt:",fatt,fatt.shape)
        print("dd  :",dd,dd.shape)
        print("kkk")
        np.savetxt("frep.dat",np.transpose((-dd,-frep)))
        np.savetxt("fatt.dat",np.transpose((dd,fatt)))

        np.savetxt("frep2.dat",np.transpose((-dd,-frep2)))
        np.savetxt("fatt2.dat",np.transpose((dd,fatt2)))

        np.savetxt("frep3.dat",np.transpose((-dd,-frep3)))
        np.savetxt("fatt3.dat",np.transpose((dd,fatt3)))

        np.savetxt("frep4.dat",np.transpose((-dd,-frep4)))
        np.savetxt("fatt4.dat",np.transpose((dd,fatt4)))


        # TODO (next spteps:)
        # 1) parametrize only the 1NN morse and see if equivalent to prior assessments.         --> DONE
        #       fitting.parameters: array([ 0.212,  1.532,  2.927])

        # 2) new parametrization with substracted 4NN                                           --> DONE
        #       fittingcorr.parameters: array([ 0.23 ,  1.495,  2.927])
        #
        # 3) run MD and look at:    a) std energy  (auswerten mit michaels skript)
        #                           b) std forces (and pic forces)
        # 4) rerun MD and look at: std energy, std forces (and pic forces)
        #       erwartung: davor "schiefes" forces bild im vergleich zu DFT, danach nicht mehr.
        # 5) use fast MD skript and check other (sweep) morse values


        #params, deltas, deltasshifted, polycoefsrealene,polycoefsrealforces,polydeltas = get_fit_forces_to_pot(
        #        foldername=self.fitfolder,
        #        filename='lon'+"_"+shellstr+"nn_"+rang+"_"+func,
        #        NN=self.nndist[self.shell],pot = func,
        #        data=self.data[self.shell,'lon',self.disp,rang].ix['points'],
        #        weights = False)

        allx = np.append(-dd[::-1],dd[1:])+np.sqrt(2*(2.07**2.))
        ally = np.append(-frep[::-1],fatt[1:])
        allycorrectedby4NN = np.append(-frep[::-1]-frep2[::-1],fatt[1:]+fatt2[1:])
        data = np.array([allx,ally])
        datacorrectedby4nn = np.array([allx,allycorrectedby4NN])

        print('data:',np.transpose(data))
        print()
        print('datacorrectedby4nn:',np.transpose(datacorrectedby4nn))
        print()
        fitting = fit_to_func( fixzeroat = np.sqrt(2*(2.07**2.)),function = 'morse',data=np.transpose(data),weights = False)
        fittingcorr = fit_to_func( fixzeroat = np.sqrt(2*(2.07**2.)),function = 'morse',data=np.transpose(datacorrectedby4nn),weights = False)

        #return schnittmenge123,frep, fatt, dd,fitting
        # fitting.parameters
        # fitting.data
        # fitting.fit
        # fitting.deltas
        # fitting.deltasshifted
        # fitting.polyfit.coefsrealene
        # fitting.polyfit.coefsrealforces
        # fitting.polyfit.deltas
        # fitting.polyfit.deltasshifted
        #
        # run  /Users/glensk/Thermodynamics/python_thermodynamics/pot_info_startjob.py  -e al -a 4.14 -sp 5 -kp 2x2x2kp -dispdirection quer -pnnforcesnpz -stradd "_vasp4_ENCUT250_GGA_2_atoms_displaced"  -v
        #
        # data,fitting = pot.disp.alle.get_forces_from_3_displaced_atoms('all',0.0,0.0)
        #
        # In [381]: fitting.parameters
        # Out[381]: array([ 0.212,  1.532,  2.927])
        #
        # In [382]: fittingcorr.parameters
        # Out[382]: array([ 0.23 ,  1.495,  2.927])
        #
        return data,fitting,fittingcorr

    def repeat_supercell_positions_forces(self):
        ''' before the whole cell can be repeated, the mapping to the original cell
        has to be remove (in other words: the atoms have to be close to their origial
        undisplaced positions, otherwise crystal_generator.center_atoms_around_atom
        does not work properly --> Therefore we use self.pr rather then self.p'''
        for idx,i in enumerate(self.dvec):
            # coords
            coord_cart = np.copy(self.pr[idx])
            self.crystal = crystal_generator.crystal()
            self.crystal.load_positions_cell(coord_cart = coord_cart, cell = self.sc)


            # coords big (repeat supercell)
            #print "coords:"
            self.crystalbig = crystal_generator.supercell()
            nsc = 2
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

    def savennforcesnpz( self, outdir, outfname, parametrize_lon = True ):
        '''
        '''
        print(utils.printred("      savennforcesnpz : "+str(outdir+'/'+outfname+".npz")))
        self.npz = outdir + "/" + outfname + ".npz"
        np.savez_compressed( outdir + "/" + outfname,
        #npz             = self.npz,  # better not, it changes rather frequently
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
        dvec1           = self.dvec1,
        dvec2           = self.dvec2,
        dvec3           = self.dvec3,
        dnorm           = self.dnorm,
        dnorm1          = self.dnorm1,
        dnorm2          = self.dnorm2,
        dnorm3          = self.dnorm3,
        atoms           = self.atoms,
        eqcoords        = self.eqcoords,
        folder          = self.folder,
        sc              = self.sc,
        scc             = self.scc,
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
        if parametrize_lon:
            self.parametrize_lon()
        return

    def loadnnforcesnpz( self, indir, infname, parametrize_data_pkl = True ):
        '''
        '''
        print(utils.printgreen("     loadnnforcesnpz : "+str(indir+'/'+infname+".npz")))
        var = np.load( indir + "/" + infname + ".npz" )
        #self.npz                = str(var['npz'])
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
        try:
            self.dvec1              = var['dvec1']
            self.dvec2              = var['dvec2']
            self.dvec3              = var['dvec3']
            self.dnorm1             = var['dnorm1']
            self.dnorm2             = var['dnorm2']
            self.dnorm3             = var['dnorm3']
        except KeyError:
            self.dvec1          = False
            self.dvec2          = False
            self.dvec3          = False
            self.dnorm1         = False
            self.dnorm2         = False
            self.dnorm3         = False
        try:
            self.scc            = var['scc']
        except KeyError:
            self.scc            = False
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
        if os.path.isdir(self.folderparam) != True:
            print("os.path.isdir(self.folderparam):",os.path.isdir(self.folderparam))
            sys.exit(self.folderparam+" does not exist, rm nnforces.npz and rerun parametrization")
        self.npz = indir + "/" + infname + ".npz"
        #print "npz:",self.parametrize_nnforcesnpz
        #print "lon:",self.parametrize_lon
        if parametrize_data_pkl:
            if self.parametrize_data_pkl == True:
                self.parametrize_lon()
        return

    def savedatapkl(self, filename):
        ''' saves pandas dataframe '''
        self.pkl = filename
        print(utils.printred("      savedatapkl : "+str(filename)))
        with open(filename, 'w') as f:
            #pickle.dump(self.list_save_load, f)
            pickle.dump([self.data,self.fit], f)
        return

    def loaddatapkl(self, filename):
        ''' loads pandas dataframe '''
        if type(filename) != str:
            print("filename:",filename)
            sys.exit("loaddatapkl  needs as input a picklefile (str) but got "+str(filename)+" of type "+str(type(filename)))
        if os.path.isfile(filename) != True:
            sys.exit("file "+filename+" does not exist! (you could redo pot_parametrize.py -e ELEMENT -p)")
        #self.data = pd.read_pickle(filename)
        print(utils.printgreen("     loaddatapkl : "+str(filename)))
        if os.path.isfile(filename) != True:
            sys.exit("Error 77 "+str(filename)+" does not exist!")
        with open(filename) as f:
            self.data, self.fit = pickle.load(f)
            #self.data, self.eqcoords, self.element = pickle.load(f)
            #*self.list_save_load() = pickle.load(f)
            #for i in enumerate(pickle.load(f)):
            #    self.list_save_load[i] = i
        self.pkl = filename
        return

    def saveTOTInpz( self, outdir, outfname ):
        print(utils.printpink("      saveTOTInpz : "+str(outdir+'/'+outfname+".npz")))
        np.savez_compressed( outdir + "/" + outfname,
        forcelo              = self.forcelo,
        forcet1              = self.forcet1,
        forcet2              = self.forcet2,
        forcelo_norm         = self.forcelo_norm,
        forcet1_norm         = self.forcet1_norm,
        forcet2_norm         = self.forcet2_norm,
        fmlforcelo           = self.fmlforcelo,
        fmlforcet1           = self.fmlforcet1,
        fmlforcet2           = self.fmlforcet2,
        fmlforcelo_norm      = self.fmlforcelo_norm,
        fmlforcet1_norm      = self.fmlforcet1_norm,
        fmlforcet2_norm      = self.fmlforcet2_norm,
        dnorm_p_to           = self.dnorm_p_to,
        dnorm_p_ti           = self.dnorm_p_ti,
        coefs_ene_to         = self.coefs_ene_to,
        coefs_ene_ti         = self.coefs_ene_ti,
        angles               = self.angles
        )
        return

    def loadTOTInpz( self, indir, infname ):
        '''
        '''
        loadfile = indir+'/'+infname+".npz"
        print(utils.printgreen("     loadTOTInpz : "+str(loadfile)))
        if os.path.isfile(loadfile) != True:
            print(loadfile,"seems not to exist ...")
        var = np.load(loadfile)
        self.forcelo               = var['forcelo']
        self.forcet1               = var['forcet1']
        self.forcet2               = var['forcet2']
        self.forcelo_norm          = var['forcelo_norm']
        self.forcet1_norm          = var['forcet1_norm']
        self.forcet2_norm          = var['forcet2_norm']
        self.fmlforcelo               = var['fmlforcelo']
        self.fmlforcet1               = var['fmlforcet1']
        self.fmlforcet2               = var['fmlforcet2']
        self.fmlforcelo_norm          = var['fmlforcelo_norm']
        self.fmlforcet1_norm          = var['fmlforcet1_norm']
        self.fmlforcet2_norm          = var['fmlforcet2_norm']
        self.dnorm_p_to = var['dnorm_p_to']
        self.dnorm_p_ti = var['dnorm_p_ti']
        self.coefs_ene_to = var['coefs_ene_to']
        self.coefs_ene_ti = var['coefs_ene_ti']
        if self.coefs_ene_to.any() == False:
            self.coefs_ene_to = False
        if self.coefs_ene_ti.any() == False:
            self.coefs_ene_ti = False
        if type(self.coefs_ene_to) == bool:
            self.coefs_ene_to = False
        if type(self.coefs_ene_ti) == bool:
            self.coefs_ene_ti = False
        #print "|||||||:",self.coefs_ene_to,type(self.coefs_ene_to)
        #print "|||||||:",self.coefs_ene_ti,type(self.coefs_ene_ti)
        self.angles = var['angles']
        #print "--->",self.FML.shape
        return

    def saveFMLnpz( self, outdir, outfname ):
        print(utils.printpink("      saveFMLnpz : "+str(outdir+'/'+outfname+".npz")))
        self.parametrize_fml_npz_path = outdir + "/" + outfname + ".npz"
        np.savez_compressed( outdir + "/" + outfname,
        parametrize_fml_npz_path           = self.parametrize_fml_npz_path,
        fml                  = self.fml,
        FML                  = self.FML,
        fml_u1nn_pottype     = self.fml_u1nn_pottype,
        fml_u1nn_potparam    = self.fml_u1nn_potparam,
        fml_u1nn_potadd      = self.fml_u1nn_potadd,
        fml_u1nn_potaddparam = self.fml_u1nn_potaddparam,
        fmlmti               = self.fmlmti,
        FMLmti               = self.FMLmti,
        fmlmto               = self.fmlmto,
        FMLmto               = self.FMLmto,
        fml_u1nn_tox         = self.fml_u1nn_tox,
        fml_u1nn_tix         = self.fml_u1nn_tix,
        angles               = self.angles,
        fprojlongold         = self.fprojlongold
        )
        return

    def loadFMLnpz( self, indir, infname ):
        '''
        '''
        loadfile = indir+'/'+infname+".npz"
        print(utils.printgreen("     loadFMLnpz : "+str(loadfile)))
        if os.path.isfile(loadfile) != True:
            print(loadfile,"seems not to exist ...")
        var = np.load(loadfile)
        #self.parametrize_fml_npz_path            = str(var['parametrize_fml_npz_path'])
        self.parametrize_fml_npz_path = indir + "/" + infname + ".npz"
        self.fml                   = var['fml']
        self.FML                   = var['FML']
        self.fml_u1nn_pottype      = str(var['fml_u1nn_pottype'])
        self.fml_u1nn_potparam     = var['fml_u1nn_potparam']
        self.fml_u1nn_potadd       = str(var['fml_u1nn_potadd'])
        self.fml_u1nn_potaddparam  = var['fml_u1nn_potaddparam']
        self.fmlmti                = var['fmlmti']
        self.FMLmti                = var['FMLmti']
        self.fmlmto                = var['fmlmto']
        self.FMLmto                = var['FMLmto']
        self.fml_u1nn_tox          = str(var['fml_u1nn_tox'])
        self.fml_u1nn_tix          = str(var['fml_u1nn_tix'])
        print("--->",self.FML.shape)
        return

    def infodos(self, verbose = False):
        ''' get DOS, DOScut for corresponding element and lattice constant;
            needs to be called after creatennforcesnpz or loadnnforcesnpz in order to have the self.element information'''
        if verbose:
            print(utils.printgreen("     infodos ..."))

        # define variables
        shells = 5
        #print "xx:",self._shells
        shells = np.array(self._shells).max()
        self.DOSlon_fileorig    = [ False for shell in range(shells+1) ]
        self.DOSlon_file        = [ False for shell in range(shells+1) ]
        self.DOSlon_filecut     = [ False for shell in range(shells+1) ]
        self.DOSlon             = [ False for shell in range(shells+1) ]
        self.DOSloncut          = [ False for shell in range(shells+1) ]

        # get DOS / DOScut
        self.jobvorlage_all = self.folder_base+'/jobvorlage_all/'
        self.jobvorlage_search = self.jobvorlage_all+'2x2x2sc_'+self.atom.symbol[0]+"_"+str(self.a)
        self.jobvorlage = glob.glob(self.jobvorlage_search)
        if len(self.jobvorlage) != 1:
            print("self.jobvorlage_search:",self.jobvorlage_search)
            print("self.jobvorlage:",self.jobvorlage)
            self.jobvorlage = False
            #sys.exit("self.jobvorlage not found (or not just one path)")
        else:
            self.jobvorlage = self.jobvorlage[0]

        ##################################################################################
        ##################################################################################
        ##################################################################################
        # load data from jobvorlage (at least the path)
        ##################################################################################
        ##################################################################################
        ##################################################################################
        # this is not a problem if it does not exist
        if self.jobvorlage:
            self.jobvorlage_avg_dudl_low_fre_file = self.jobvorlage+"/avg_dUdL_low_fre"
            if os.path.isfile(self.jobvorlage_avg_dudl_low_fre_file) == True:
                self.jobvorlage_avg_dudl_low_fre = np.loadtxt(self.jobvorlage_avg_dudl_low_fre_file)


            ######################################################################
            # tveclonall.dat
            self.jobvorlage_vecnormlon_file = self.jobvorlage+"/tveclonall.dat"   # somewhat shorter
            self.jobvorlage_vecnormlon_file = self.jobvorlage+"/atoms_1nn_all_f/dfn_1.0"
            if os.path.isfile(self.jobvorlage_vecnormlon_file) != True:
                print(self.jobvorlage_vecnormlon_file + " self.jobvorlage_vecnormlon_file does not exist 3 !")
                print("first: go to the lowfolder with all your lambda_* folders and run DOS_POSITIONS_auswerten.py -c")
                print("first: go to the lowfolder with all your lambda_* folders and run DOS_POSITIONS_auswerten.py -corb -element al")

                sys.exit(self.jobvorlage_vecnormlon_file + " self.jobvorlage_vecnormlon_file does not exist 3 ! (this has to be done with DOS_POSITIONS_auswerten.py -corb -dosmax -inte")
            if verbose:
                print("     self.jobvorlage_vecnormlon_file   :",self.jobvorlage_vecnormlon_file)

            ######################################################################
            # DOSlon  DOSDOSDOS
            ######################################################################
            for shell in range(shells): self.DOSlon_fileorig[shell+1] = self.jobvorlage+"/longvecnorm_"+str(shell+1)
            for shell in range(shells): self.DOSlon_file[shell+1] = self.jobvorlage+"/longvecnorm_"+str(shell+1)+"_DOS"
            for shell in range(shells): self.DOSlon_filecut[shell+1] = self.jobvorlage+"/longvecnorm_"+str(shell+1)+"_DOScut"
            for shell in range(shells):
                if os.path.isfile(self.DOSlon_file[shell+1]) == True:
                    self.DOSlon[shell+1] = np.loadtxt(self.DOSlon_file[shell+1])
                else:
                    print(utils.printred(self.DOSlon_file[shell+1]+" not found!"))

            # create self.DOSlon_filecut in case it does not exist
            for shell in range(shells):
                if os.path.isfile(self.DOSlon_filecut[shell+1]) != True:
                    if os.path.isfile(self.DOSlon_fileorig[shell+1]) == True:
                        print("creating DOScut")
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
                print("     self.DOSlon: OK")
                print("     self.DOSloncut: OK")
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
        #                           V   V ... atoms to check
        #      pot.disp.alle.sc  ==   array([[ 7.24,  0.  ,  0.  ], ----------- 14.48^3 ist vergroesserte sc -------------
        #                                    [ 0.  ,  7.24,  0.  ], ------------------------------------------VVVV
        #                                    [ 0.  ,  0.  ,  7.24]])
        #                                                                fuer xx = 3.99                   fuer Ni = 3.62
        listcheck = [
        [ 'fcc', 1, 'lon',32, 24, 126 ],   # [ 1.995,  1.995,  0.   ], ([-1.995, -1.995,   0.   ] == [ 12.67,  12.67,   0.  ])
        [ 'fcc', 2, 'lon',32,  4,  36 ],   # [ 3.99,  0.  ,  0.     ], ([-3.99,  0.  ,   0.     ] == [ 10.86,   0.  ,   0.  ])
        [ 'fcc', 3, 'lon',32, 12, 239 ],   # [ 3.99 ,  1.995,  1.995], ([-3.99 ,  -1.995, -1.995] == [ 10.86,  12.67,  12.67])
        [ 'fcc', 4, 'lon',32,  6, 102 ],   # [ 3.99,  3.99,  0.  ]  ,  ([-3.99,  -3.99,  0.  ]    == [ 10.86,  10.86,   0.  ])

        [ 'fcc', 1, 'tox',32, 8, 203 ],   # [ 0.   ,  2.065,  2.065],[ 0.   , -2.065, -2.065]
        [ 'fcc', 2, 'tox',32, 2,  66 ],   # [ 0.  ,  4.13,  0.00], [0,  -4.13  ,   0.     ]
        [ 'fcc', 3, 'tox',32, 0,  0  ],   # [ 3.99 ,  1.995,  1.995], [-3.99 ,  -1.995, -1.995]
        [ 'fcc', 4, 'tox',32, 0,  0  ],   # [ 3.99,  3.99,  0.  ]  , [-3.99,  -3.99,  0.  ]

        [ 'fcc', 1, 'tix',32, 60,  0  ],   # [ -2.065   ,  -2.065, 0.0],[ 0.   , -2.065, -2.065]
        [ 'fcc', 2, 'tix',32, 0,  0  ],   # [ 0.  ,  4.13,  0.00], [0,  -4.13  ,   0.     ]
        [ 'fcc', 3, 'tix',32, 0,  0  ],   # [ 3.99 ,  1.995,  1.995], [-3.99 ,  -1.995, -1.995]
        [ 'fcc', 4, 'tix',32, 0,  0  ],   # [ 3.99,  3.99,  0.  ]  , [-3.99,  -3.99,  0.  ]

        [ 'fcc', 1, 'lon',108, 81, 429 ],   # [ 1.995,  1.995,  0.   ], [-1.995, -1.995,  0.   ]
        [ 'fcc', 2, 'lon',108,  9, 126 ],   # [ 3.99,  0.  ,  0.     ], [-3.99,  0.  ,  0.     ]
        [ 'fcc', 3, 'lon',108, 36, 809 ],   # [ 3.99 ,  1.995,  1.995], [-3.99 ,  -1.995,-1.995]
        [ 'fcc', 4, 'lon',108, 12, 348 ],   # [ 3.99,  3.99,  0.  ]   , [-3.99,  -3.99,  0.  ]

        [ 'fcc', 1, 'lon',108, 81,  429 ],   # [ 1.995,  1.995,  0.   ], [-1.995, -1.995,  0.   ]   Ni: sc = 10.86^3 grosse sc = 21.72^3  == [ 19.91,  19.91,   0.  ])
        [ 'fcc', 2, 'lon',108, 9 ,  126 ],   # [ 3.99,  0.  ,  0.     ], [-3.99,  0.  ,  0.     ]                                            [ 18.1   0.    0. ]
        [ 'fcc', 3, 'lon',108, 36,  809 ],   # [ 3.99 ,  1.995,  1.995], [-3.99 ,  -1.995,-1.995]                                            [ 18.1   19.91  19.91]
        [ 'fcc', 4, 'lon',108, 12,  348 ],   # [ 3.99,  3.99,  0.  ]   , [-3.99,  -3.99,  0.  ]                                              [ 18.1  18.1   0. ]
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

    def parametrize_lon(self,verbose=False):
        ''' creates parametrization for long vectors
            make sure creatennforcesnpz() was called before '''
        # Example: Ir, 2x2x2sc 1nn 0.5 quer displacement
        # Folder: /Users/glensk/Understand_distributions/displacements_dense/Ir/2x2x2sc_quer_3x3x3kp/3.99Ang_0.5
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
        print(utils.printblue("     pot_parametrize.parametrize_lon from T=0K displacements ... 0"))

        # if pkl file does not exist
        if type(self.pkl) != bool:
            if os.path.isfile(self.pkl) != True:
                print("os.path.isfile(self.pkl):", os.path.isfile(self.pkl))
                print("os.path.isfile(self.pkl) is not true .... going to create it")
                self.parametrize_data_pkl = True

        # in case there is no pkl we have to parametrize pkl
        #if os.path.isfile(self.folderparam+'/data.pkl') != True:
        #    self.parametrize_data_pkl = True

        # in case we parametrize npz we also want to parametrize pkl
        if self.parametrize_nnforcesnpz == True:
            self.parametrize_data_pkl = True

        # if we dont have to parametrize pkl and the pkl file exists
        if self.verbose == True:
            print("self.parametrize_data_pkl:",self.parametrize_data_pkl)
        if self.parametrize_data_pkl == False:
            #if os.path.isfile(self.folderparam+"/data.pkl") == True:
            #print utils.printred("      loading "+self.folderparam+'/data.pkl')
            #self.loaddatapkl(self.folderparam+'/data.pkl')
            #print "self.pkllllllllllllll:",self.pkl
            self.loaddatapkl(self.pkl)
            #self.loadFitpd(self.folderparam+'/fit.pkl')
            # xdir has usually all the data whil quer has only quer data
            return

        # in case pkl file exists, we want to load it and add to/replace items of this file

        #########################################################################
        # create self.func_lon  (wrong for 2nd shell 2x2x2sc xdir)
        #########################################################################
        print(utils.printblue("     pot_parametrize.parametrize_lon from T=0K displacements ... 1 (pot.disp.quer.Pr)"))
        self.func_lon = np.zeros((self.Pr.shape[1],self.Pr.shape[0],2))
        self.func_to = np.zeros((self.Pr.shape[1],self.Pr.shape[0],2))
        self.func_ti = np.zeros((self.Pr.shape[1],self.Pr.shape[0],2))

        if self.verbose:
            print("self.Pr.shape:",self.Pr.shape)

        for i in np.arange(self.Pr.shape[1]):   # 0 ... 255
            ############################################################
            # Pr.shape: (200, 256, 3)
            # --> 200 displacements
            # --> 256 atoms in repeated supercell from 32 atomcell
            # --> 3   for xyz
            # pot.disp.quer.P[:,24 ] will always give the same entry of the position of the first NN [ 1.81,  1.81,  0.  ]
            # pot.disp.quer.P[:,126] will always give the same entry of the position of the first NN [ 12.67,  12.67,   0.  ]
            # pot.disp.quer.Pr[:,24] will give the vector from the displaced atom to the 1NN so from [1.81 , 1.81 ,  0.   ] to  [ 0.403,  0.403,  0.   ]
            # pot.disp.quer.Pr[:,126] will give the vector from the displaced atom to the 1NN so from [-1.81 , -1.81 ,  0.   ] to [-3.217, -3.217,  0.   ]
            # np.linalg.norm(pot.disp.quer.Pr[95,24])-3.62/np.sqrt(2.) = -0.9499979
            # np.linalg.norm(pot.disp.quer.Pr[96,24])-3.62/np.sqrt(2.) = -0.9599979
            # np.linalg.norm(pot.disp.quer.Pr[97,24])-3.62/np.sqrt(2.) = -0.9699979
            # np.linalg.norm(pot.disp.quer.F[95,24]) = 20.669026784717321   [ 14.615,  14.615,   0.   ] == cat 3.62Ang_0.95/forces_OUTCAR | head -25 | tail -1
            # np.linalg.norm(pot.disp.quer.F[96,24]) = 23.904260925961502   [ 16.903,  16.903,   0.   ] == cat 3.62Ang_0.96/forces_OUTCAR | head -25 | tail -1
            # np.linalg.norm(pot.disp.quer.F[97,24]) = 22.231901062553515   [ 15.72,  15.72,   0.  ]    == cat 3.62Ang_0.97/forces_OUTCAR | head -25 | tail -1
            # in Ni waren die laeufe (quer auslenkung) 3.62_0.96 und 3.62_1.12 kaputt; habe diese geloescht
            #
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
        print(utils.printblue("     pot_parametrize.parametrize_lon from T=0K displacements ... 2 (pot.disp.quer.Pr)"))
        if self.verbose:
            print("self._shells:",self._shells)
        for self.shell in self._shells:
            self.param_atoms = self.parametrize_whatatom(shell = self.shell, contrib = 'lon')
            if self.shell == 2 and self.element == 'Cu' and self.disp == 'quer':
                print("NN:",self.nndist[2])
                print(np.linalg.norm(self.Pr[:,i],axis=1))
                changesign = np.nonzero(np.linalg.norm(self.Pr[:,i],axis=1) < self.nndist[2])
                print(changesign)
                pass


        #########################################################################
        # sort forces only here, since from here on those are not sorted as self.P self.Pr self.F anymore
        #########################################################################
        print(utils.printblue("     pot_parametrize.parametrize_lon from T=0K displacements ... 3 (pot.disp.quer.Pr)"))
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


        print(utils.printblue("     pot_parametrize.parametrize_lon from T=0K displacements ... 4 (pot.disp.quer.Pr)"))
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
        print(utils.printblue("     pot_parametrize.parametrize_lon from T=0K displacements ... 5 (pot.disp.quer.Pr)"))
        for self.shell in self._shells:
            shellstr = str(self.shell)
            self.param_atoms = self.parametrize_whatatom(shell = self.shell, contrib = 'lon')
            print("     shell:",self.shell,"self.param_atoms:",self.param_atoms,"nndist[shell]:",self.nndist[self.shell],"self.disp:",self.disp)


            # concat pos / negside and append 0.0 point and sort
            # pot.param.parametrize_whatatom(2,'lon') array([ 24, 126])
            # pot.param.parametrize_whatatom(2,'lon') array([ 4, 36])
            # pot.param.parametrize_whatatom(3,'lon') array([ 12, 239])
            # pot.param.parametrize_whatatom(4,'lon') array([ 6, 102])
            lonall = np.concatenate((self.func_lon[self.param_atoms[0]],self.func_lon[self.param_atoms[1]]))
            # add (0,0) point
            lonall = np.concatenate((lonall,np.array([[self.nndist[self.shell],0.0]])))

            # remove duplicates
            lonall = utils.remove_duplicates_of_2d_array_within_tolerance(lonall,1e-6, 1e-4)
            lonall = lonall[lonall[:,0].argsort()]
            self.data[self.shell,'lon',self.disp,'alle'].ix['points'] = lonall
            self.data[self.shell,'lon',self.disp,'alle'].ix['pointsshifted'] = np.transpose([lonall[:,0]-self.nndist[self.shell],lonall[:,1]])

            # at this place we have to make sure we loaded DOSlon stuff (DOSinfo)
            print("#############################")
            print(lonall)
            print(self.DOSlon[1])
            print("#############################")
            xNNdist = lonall[np.where(lonall[:,1] == 0)[0]][0,0]
            print(xNNdist)
            print("self.shell:",self.shell,self.nndist[self.shell])
            print("&&&&&&&&&&&&:",type(self.DOSlon[self.shell]))
            print("self.DOSlon[self.shell]",self.DOSlon[self.shell])
            # 1NN: doscut is always -0.6 bis + 1.15
            # 1NN: dos is always -0.8 bis + 1.5   (das sollte noch fuer alle funktionieren, ist aber nicht so sicher wie die cut werte)
            # 2NN: doscut is always -1 bis 1
            # 3NN: doscut is always -1 bis 1
            if type(self.DOSlon[self.shell]) == bool:
                nndist = self.nndist[self.shell]
                if self.shell == 1:
                    self.DOSlon[self.shell] = np.array([[nndist-0.8,0.],[nndist,1.],[nndist+1.5,0]])
                    self.DOSloncut[self.shell] = np.array([[nndist-0.6,0.],[nndist,1.],[nndist+1.15,0]])
                if self.shell == 2 or self.shell == 3:
                    self.DOSlon[self.shell] = np.array([[nndist-1.,0.],[nndist,1.],[nndist+1.,0]])
                    self.DOSloncut[self.shell] = np.array([[nndist-1.,0.],[nndist,1.],[nndist+1.,0]])

                print(self.DOSlon[self.shell])

            if type(self.DOSlon[self.shell]) != bool:
                lonallDOS      = utils.cut_function_at_DOS(lonall,self.DOSlon[self.shell])
                lonallDOScut   = utils.cut_function_at_DOS(lonall,self.DOSloncut[self.shell])

                self.data[self.shell,'lon',self.disp,'dos'].ix['points'] = lonallDOS
                self.data[self.shell,'lon',self.disp,'doscut'].ix['points'] = lonallDOScut
                self.data[self.shell,'lon',self.disp,'dosneg'].ix['points'] = lonallDOS[np.nonzero(lonallDOS[:,0] <= self.nndist[self.shell]+0.1)]  # we want the +-0.1 to add some data ... otherwise weights necessary
                self.data[self.shell,'lon',self.disp,'dospos'].ix['points'] = lonallDOS[np.nonzero(lonallDOS[:,0] >= self.nndist[self.shell]-0.1)]  # we watn the +-0.1 to add some data ... otherwies weights necessary
                print(utils.printyellow("done",self.shell))


                self.data[self.shell,'lon',self.disp,'dos'].ix['pointsshifted'] = np.transpose([lonallDOS[:,0]-self.nndist[self.shell],lonallDOS[:,1]])
                self.data[self.shell,'lon',self.disp,'doscut'].ix['pointsshifted'] = np.transpose([lonallDOScut[:,0]-self.nndist[self.shell],lonallDOScut[:,1]])
            else:
                self.data[self.shell,'lon',self.disp,'dos'].ix['points'] = False
                self.data[self.shell,'lon',self.disp,'doscut'].ix['points'] = False
                self.data[self.shell,'lon',self.disp,'dos'].ix['pointsshifted'] = False
                self.data[self.shell,'lon',self.disp,'doscut'].ix['pointsshifted'] = False
                print(utils.printred("DOS not cut for shell "+str(self.shell)+"since self.DOSlon[self.shell] not defined"))



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
                        print("saving:",filename)
                        #print data[:3],type(data),data.shape
                        np.savetxt(filename,self.data[self.shell,'lon',self.disp,rang].ix['points'],fmt="%.7f")
                if os.path.isfile(filenameshifted) != True:
                    data = self.data[self.shell,'lon',self.disp,rang].ix['pointsshifted']
                    if type(data) != bool:
                        print("saving shifted data:",filenameshifted)
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
        print(utils.printblue("     pot_parametrize.parametrize_lon from T=0K displacements ... 6 (pot.disp.quer.Pr)"))
        for self.shell in self._shells:
            shellstr = str(self.shell)

            for rang in self._range: # ['alle', 'dos', 'doscut', 'dosneg' ]
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
                            print("         ",self.shell,rang,func,self.fitting.parameters)

                        ##################################################################
                        # write to dataframe
                        # --> pot.disp.alle.data[1,'lon','quer','dos']
                        # --> pot.disp.alle.fit[1,'lon','quer','dos','morse']
                        #
                        # pot.disp.alle.fit[1,'lon','quer','dos']
                        # pot.param.data[1,'lon','quer','dos']
                        # pot.param.fit[1,'lon','quer','alle','morse']
                        #
                        ##################################################################
                        print("pot.param.data[",self.shell,'lon',self.disp,rang,func)
                        #print "self.fitting.parameters:",self.fitting.parameters
                        #print "self.fitting.deltas:",self.fitting.deltas
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
                            print("         ","spline FIT")
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
                                print("saving:",filenameshifteddiff)
                                np.savetxt(filenameshifteddiff,self.fit[ self.shell,'lon',self.disp,rang,func].ix['diffshifted'],fmt='%.7f')
                                np.savetxt(filenameshifteddiffpolydiff,self.fit[ self.shell,'lon',self.disp,rang,func].ix['polybestdiffshifted'],fmt='%.7f')
                        np.savetxt(filenameshifteddiffsplinediff,splshifteddiff,fmt='%.7f')

                    if type(self.fitting.parameters) == bool:
                        print("         ",self.shell,rang,func,utils.printred(self.fitting.parameters))
                    else:
                        print("         ",self.shell,rang,func,self.fitting.parameters)

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
        #self.savedatapkl(self.pkl)
        #################################################################################
        print(utils.printblue("     pot_parametrize_function complete ... :) "))
        print("")
        return



class forcesneighbors_all_disp( object ):
    def __init__(self):
        self.verbose                            = False
        self.folder_displacement_direction_all  = False
        self.pkl                                = False
        self.npz                                = False
        self.alatsstringlist                    = False  # [3.99]
        self.parametrize_nnforcesnpz            = False
        self.parametrize_data_pkl               = False
        self.parametrize_fml_npz                = False
        self.create_qubus_npz                   = False

        self.folder_base                = "/Users/glensk/Understand_distributions/"

        self.alle = forcesneighbors()
        self.disp_all = self.alle.disp_all
        self.update_what_to_calculate_and_write_self_xdir_quer()
        self.xdir = False
        self.quer = False
        self.midd = False
        self.dnnd = False
        self.qubus = False

        # folder_displacement_direction_all  ( necessary for loading forces )
        # [ /Users/glensk/Understand_distributions//displacements_dense/Ir/2x2x2sc_quer_3x3x3kp,
        #   /Users/glensk/Understand_distributions//displacements_dense/Ir/2x2x2sc_xdir_3x3x3kp]
        #
        # pkl : /Users/glensk/Understand_distributions//displacements_dense/Ir//2x2x2sc_data_3x3x3kp.pkl
        #
        # alatsstringlist = [ 3.99 ]
        # parametrize = True or False

    def update_what_to_calculate_and_write_self_xdir_quer(self, initialize = False, onlyprint = False):
        if onlyprint == True:
            print("1:npz:(nnforces.npz):",self.alle.parametrize_nnforcesnpz,"self.alle.parametrize_nnforcesnpz")
            print("2:pkl:",self.alle.parametrize_data_pkl,"self.alle.parametrize_data_pkl")
            print("3:anz:",self.alle.parametrize_fml_npz,"self.alle.parametrize_fml_npz")
            return

        self.alle.parametrize_nnforcesnpz            = self.parametrize_nnforcesnpz
        self.alle.parametrize_data_pkl               = self.parametrize_data_pkl
        self.alle.parametrize_fml_npz                = self.parametrize_fml_npz
        self.alle.pkl = self.pkl

        #self.alle.angles        = self.
        #self.alle.forcelo[dispidx,atomidx] = forcelo
        #self.alle.forcet1[dispidx,atomidx] = forcet1
        #self.alle.forcet2[dispidx,atomidx] = forcet2
        #self.alle.forcelo_norm[dispidx,atomidx] = forcelo_norm
        #self.alle.forcet1_norm[dispidx,atomidx] = forcet1_norm
        #self.alle.forcet2_norm[dispidx,atomidx] = forcet2_norm



        if self.alle.parametrize_nnforcesnpz == True:
            self.alle.parametrize_data_pkl = True
            self.alle.parametrize_fml_npz = True

        if self.alle.parametrize_data_pkl == True:
            self.alle.parametrize_fml_npz = True


        if self.alle.disp == 'xdir':
            if initialize == True:
                self.xdir = forcesneighbors()
            self.xdir           = copy.deepcopy(self.alle)
            #self.xdir.angles    = copy.deepcopy(self.alle.angles)
            #self.xdir.forceslo  = copy.deepcopy(self.alle.forceslo)

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

        if self.alle.disp == 'qubus':
            if initialize == True:
                self.qubus = forcesneighbors()
            self.qubus = copy.deepcopy(self.alle)
        return


    def load_all_pkl_npz(self):
        print(utils.printblue("     pot_parametrize.load_all_pkl_npz ..."))
        ###################################################################################
        # 1. load nnforces.npz for every displacement
        # 2. load parametrization_long == xxx_data_xxx.pkl
        # 3. load FML.npz == forces (for corresponding direction) - longforce in quer direction
        ###################################################################################
        self.update_what_to_calculate_and_write_self_xdir_quer()
        self.update_what_to_calculate_and_write_self_xdir_quer(onlyprint=True)

        ###################################################################################
        ###################################################################################
        # in case we do have everything :
        # hier sollte eine schnelloption hin welche
        # 1. alle npz files einlaedt ( sowohl nnforces.npz als auch FML.npz
        # 2. das data_pkl einlaedt
        ###################################################################################
        ###################################################################################
        coefs_ene_ti = False
        coefs_ene_to = False
        if self.alle.parametrize_nnforcesnpz == False and self.alle.parametrize_data_pkl == False and self.alle.parametrize_fml_npz == False:
            for folder in self.folder_displacement_direction_all:
                #print "iii:",folder
                if os.path.isfile(folder+"/EqCoords_direct") != True:
                    sys.exit(folder+"/EqCoords_direct does not exist; need routing to get Eqcoords")
                for a in self.alatsstringlist:
                    #print "aaa:",a
                    if os.path.isfile(folder+"/nnforces.npz") != True:
                        print(utils.printred("     "+folder+"/nnforces.npz does not exist! have to create it 1!"))
                        print(utils.printred("      run first pot_info_startjob.py -e element -a alat -dispdirection {xdir,quer,...} -pnnforcesnpz"))
                        sys.exit()
                    if os.path.isfile(folder+"/FML.npz") != True:
                        print(utils.printred("     "+folder+"/FML.npz does not exist! have to create it 2!"))
                        print(utils.printred("      run first pot_info_startjob.py -e element -a alat -dispdirection {xdir,quer,...} -ppkldata"))
                        sys.exit()
                    self.alle.loadnnforcesnpz(folder,"nnforces")
                    #self.alle.creatennforcesnpz(folder = folder, a = a)    # macht zum schluss    # self.parametrize_lon()  which makes infodos
                    self.alle.loadFMLnpz(folder,"FML")

                    if os.path.isfile(folder+"/TOTI.npz") == True:
                        self.alle.loadTOTInpz(folder,"TOTI")
                    else:
                        print("currently self.parametrize_tox is not run! next ten lines grayed out")
                        #self.parametrize_tox(go_through_one_folder=folder)

                    #self.update_what_to_calculate_and_write_self_xdir_quer()
                    #if self.alle.disp == 'quer':
                    #    coefs_ene_ti = self.quer.coefs_ene_ti
                    #if self.alle.disp == 'xdir':
                    #    coefs_ene_to = self.xdir.coefs_ene_to


            self.alle.loaddatapkl(self.pkl)

        if self.alle.parametrize_nnforcesnpz == False and self.alle.parametrize_data_pkl == False and self.alle.parametrize_fml_npz == True:
            for folder in self.folder_displacement_direction_all:
                for a in self.alatsstringlist:
                    self.alle.loadnnforcesnpz(folder,"nnforces")
                    self.alle.loadTOTInpz(folder,"TOTI")
                    self.update_what_to_calculate_and_write_self_xdir_quer()
                    if self.alle.disp == 'quer':
                        coefs_ene_ti = self.quer.coefs_ene_ti
                    if self.alle.disp == 'xdir':
                        coefs_ene_to = self.xdir.coefs_ene_to

            print("coefs_ene_ti:",coefs_ene_ti)
            print("coefs_ene_to:",coefs_ene_to)
            if type(coefs_ene_ti) == bool or type(coefs_ene_to) == bool:
                sys.exit('coefs_ene_ti or coefs_ene_to is bool')
            self.get_FMLnpz(shell = 1, to = coefs_ene_to, ti = coefs_ene_ti)
            self.alle.parametrize_nnforcesnpz = False
            return

        ###################################################################################
        # 1. in case we neeed to create nnforces.npz --> a) xxx_data_xxx.pkl will be recreated automatically ---> then we also need to recreated FML.npz
        # after this disable creating
        ###################################################################################
        if self.alle.parametrize_nnforcesnpz == True: # or self.alle.parametrize_fml_npz == True or self.alle.parametrize_data_pkl == True:
            for folder in self.folder_displacement_direction_all:  # quer,xdir,midd,3nnd, ...
                print("creating for:",folder)
                if os.path.isfile(folder+"/EqCoords_direct") != True:
                    sys.exit(folder+"/EqCoords_direct does not exist; need routing to get Eqcoords")
                    #shutil.copyfile(self.eqcoordsfile_getpot,folder+"/EqCoords_direct")

                # this needs to happen in the corresponding class to have all variables there for xdir and quer
                ###################################################################################
                # load/or parametrize forces for every displacement
                ###################################################################################
                for a in self.alatsstringlist:
                    # self.alle.creatennforcesnpz macht zum schluss
                    # --> self.parametrize_lon()  which makes
                    # -->       --> infodos
                    self.alle.creatennforcesnpz(folder = folder, a = a)    # macht zum schluss    # self.parametrize_lon()  which makes infodos
                    #self.alle.infodos()
                    #print 'self.pkl:',self.pkl  # == /Users/glensk/Understand_distributions//displacements_dense/Ir//2x2x2sc_data_3x3x3kp.pkl

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
                    #if self.alle.parametrize_fml_npz == False and os.path.isfile(folder+"/FML.npz") == True:
                        #print 'lonading ',folder,'FML'
                    #    self.alle.loadFMLnpz(folder,"FML")
                    self.update_what_to_calculate_and_write_self_xdir_quer()

            ##########################################################################
            # since we do have new nnforces.npz we also need new FML.npz
            # this can only be done once all FML.npz exist (for all quer/xdir/... folder)
            # apparently to create FML.npz we first need the pkl file! (which makes ...
            # ... sense since pkl contains the longitudinal + correction fits)
            # but this is ok since new xxx_data_xxx.pkl file is created with every new
            #
            ##########################################################################
        ##################################################################################
        # 1. at this stage we do have for sure the nnforces.npz loaded
        #    -> NECESSARY to save xxx_data_xxx pkl since otherwies the parametrize_long
        #    which created the xxx_data_xxx is not saved to xxx_data_xxx.pkl
        #    -> NECESSARY to make new FML.npz ... curretnly this will be done for every folder xdir/quer/midd ...
        #    --> set self.alle.parametrize_fml_npz = True
        ##################################################################################
        if self.alle.parametrize_data_pkl == True or os.path.isfile(self.pkl) == False:
            for folder in self.folder_displacement_direction_all:  # quer,xdir,midd,3nnd, ...
                for a in self.alatsstringlist:
                    self.alle.creatennforcesnpz(folder,a)
                    self.update_what_to_calculate_and_write_self_xdir_quer()
                    self.alle.savedatapkl(self.pkl)
            self.alle.parametrize_nnforcesnpz   = False
            self.alle.parametrize_data_pkl      = False
            self.parametrize_data_pkl           = False

        #sys.exit('out here')
        #if self.alle.parametrize_nnforcesnpz == True or self.alle.parametrize_data_pkl == True:
        #    self.alle.savedatapkl(self.pkl)
        #    self.alle.parametrize_fml_npz = True

        #self.alle.parametrize_nnforcesnpz = False



        ###################################################################################
        # 2. load (and evtl save) data.pkl fits if no parametrization necessary
        ###################################################################################
        if self.parametrize_data_pkl == False and os.path.isfile(self.pkl) == True:
            self.alle.loaddatapkl(self.pkl)
        elif self.parametrize_data_pkl == True or os.path.isfile(self.pkl) == False:   # (xdir, quer, 3nnd, 4nnd, False = * = all)
            if len(self.folder_displacement_direction_all) > 1: # we want to have xdir and quer
                self.alle.savedatapkl(self.pkl)
                self.parametrize_data_pkl = False
                self.alle.parametrize_data_pkl = False
            else:
                sys.exit(self.pkl+" was not saved since you would need to load at least xdir and quer but you have only: "+str(self.folder_displacement_direction_all[0]))
        else:
            sys.exit("pdk not loaded ... check why!")

        print(utils.printredbold('jetzt sollte xxx_data_xxx.pkl geladen sein'))


        ###################################################################################
        # load (and evtl. save) FML
        ###################################################################################
        #if self.parametrize_fml_npz != True and os.path.isfile(self.npz) == True:
        #    self.loadFMLnpz(self.npz)
        #elif self.parametrize_fml_npz == True or os.path.isfile(self.npz) != True:

        #print 'quer;',self.quer.data[1,'lon','quer','dos'].ix['morse']
        #print 'xdir;',self.xdir.data[1,'lon','quer','dos'].ix['morse']
        #print 'alle;',self.alle.data[1,'lon','quer','dos'].ix['morse']   -- << -- only here

        # das naechste muss eigentlich nur gemacht werden wenn FML noch nicht vorhanden ist und irgendeine parametrisierung angestanden hat (dann sollte man updaten)
        self.update_what_to_calculate_and_write_self_xdir_quer(onlyprint = True)

        if self.alle.parametrize_fml_npz == True or os.path.isfile(folder+"/FML.npz") == False:
            self.get_FMLnpz()
        print(utils.printblue("     pot_parametrize.load_all_pkl_npz complete ..."))
        return

    def parametrize_tox(self, go_through_one_folder = False):
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
        check the statement above once the longitudinal contribution is subtracted!
        --> tox including long
        --> tox pure  ( here we should extract the quer very well fitted long part)
        '''
        print(utils.printblue("     pot_parametrize.parametrize_tox ..."))
        self.alle.coefs_ene_to = False
        self.alle.coefs_ene_ti = False
        #@for folder in self.folder_displacement_direction_all:  # quer,xdir,midd,3nnd, ...
        #@    print "creating for:",folder
        #@    if os.path.isfile(folder+"/EqCoords_direct") != True:
        #@        sys.exit(folder+"/EqCoords_direct does not exist; need routing to get Eqcoords")
        #@        #shutil.copyfile(self.eqcoordsfile_getpot,folder+"/EqCoords_direct")

        #@    # this needs to happen in the corresponding class to have all variables there for xdir and quer
        #@    ###################################################################################
        #@    # load/or parametrize forces for every displacement
        #@    ###################################################################################
        #@    for a in self.alatsstringlist:
        #@        # self.alle.creatennforcesnpz macht zum schluss
        #@        # --> self.parametrize_lon()  which makes
        #@        # -->       --> infodos
        #@        for atom in np.arange(self.Pr.shape[1]):   # 0 ... 255
        #@            ############################################################
        #@            # Pr.shape: (200, 256, 3)
        #@            # --> 200 displacements
        #@            # --> 256 atoms in repeated supercell from 32 atomcell
        #@            # --> 3   for xyz
        #@            # pot.disp.quer.P[:,24 ] will always give the same entry of the position of the first NN [ 1.81,  1.81,  0.  ]
        #@            # pot.disp.quer.P[:,126] will always give the same entry of the position of the first NN [ 12.67,  12.67,   0.  ]
        #@            # pot.disp.quer.Pr[:,24] will give the vector from the displaced atom to the 1NN so from [1.81 , 1.81 ,  0.   ] to  [ 0.403,  0.403,  0.   ]
        #@            # pot.disp.quer.Pr[:,126] will give the vector from the displaced atom to the 1NN so from [-1.81 , -1.81 ,  0.   ] to [-3.217, -3.217,  0.   ]
        #@            # np.linalg.norm(pot.disp.quer.Pr[95,24])-3.62/np.sqrt(2.) = -0.9499979
        #@            # np.linalg.norm(pot.disp.quer.Pr[96,24])-3.62/np.sqrt(2.) = -0.9599979
        #@            # np.linalg.norm(pot.disp.quer.Pr[97,24])-3.62/np.sqrt(2.) = -0.9699979
        #@            # np.linalg.norm(pot.disp.quer.F[95,24]) = 20.669026784717321   [ 14.615,  14.615,   0.   ] == cat 3.62Ang_0.95/forces_OUTCAR | head -25 | tail -1
        #@            # np.linalg.norm(pot.disp.quer.F[96,24]) = 23.904260925961502   [ 16.903,  16.903,   0.   ] == cat 3.62Ang_0.96/forces_OUTCAR | head -25 | tail -1
        #@            # np.linalg.norm(pot.disp.quer.F[97,24]) = 22.231901062553515   [ 15.72,  15.72,   0.  ]    == cat 3.62Ang_0.97/forces_OUTCAR | head -25 | tail -1
        #@            # in Ni waren die laeufe (quer auslenkung) 3.62_0.96 und 3.62_1.12 kaputt; habe diese geloescht
        #@            for disp in np.arange(self.Pr.shape[0]):   # 0 ... 255
        #@                force = self.F[disp,i]
        #@                forcelo  = utils.project_vector(force,longvec)
        #@                forcet1 = utils.project_vector(force,t1)
        #@                forcet2 = utils.project_vector(force,t2)
        #@                print atom,disp,"force:",force,"forcelo:",forcelo,"ft1:",forcet1,"ft2:",forcet2,"CHECK:",(force-forcet1)-forcet2
        #print "gg1;",type(go_through_one_folder)
        #print "gg;",go_through_one_folder,type(go_through_one_folder)
        if type(go_through_one_folder) == bool:
            go_through_folder = self.folder_displacement_direction_all
        elif type(go_through_one_folder) == str:
            go_through_folder = [go_through_one_folder]
        elif type(go_through_one_folder) == list:
            go_through_folder = go_through_one_folder


        #for folderidx,folder in enumerate(self.folder_displacement_direction_all):
        for folderidx,folder in enumerate(go_through_folder):
            print("folder:",folder)
            #sys.exit()
            for a in self.alatsstringlist:
                if os.path.isfile(folder+'/TOTI.npz') == True:
                    self.alle.loadTOTInpz(folder, 'TOTI')
                else:
                    print(utils.printred("TOTI.npz does not exist, a:",a,"##################:",self.alle.disp))

                    self.alle.loadnnforcesnpz(     folder, 'nnforces' )  # loads npz (forces)
                    self.alle.loadFMLnpz(  folder, 'FML'  )
                    print("     direction:",self.alle.disp)
                    self.alle.angles        = np.zeros((self.alle.Pr.shape[0],self.alle.Pr.shape[1]))

                    self.alle.dvec_p_to0    = np.zeros((self.alle.Pr.shape[0],self.alle.Pr.shape[1],3))
                    self.alle.dvec_p_ti0    = np.zeros((self.alle.Pr.shape[0],self.alle.Pr.shape[1],3))
                    self.alle.dnorm_p_to    = np.zeros((self.alle.Pr.shape[0],self.alle.Pr.shape[1]))
                    self.alle.dnorm_p_ti    = np.zeros((self.alle.Pr.shape[0],self.alle.Pr.shape[1]))

                    self.alle.forcelo_norm  = np.zeros((self.alle.Pr.shape[0],self.alle.Pr.shape[1]))
                    self.alle.forcet1_norm  = np.zeros((self.alle.Pr.shape[0],self.alle.Pr.shape[1]))
                    self.alle.forcet2_norm  = np.zeros((self.alle.Pr.shape[0],self.alle.Pr.shape[1]))
                    self.alle.forcelo       = np.zeros((self.alle.Pr.shape[0],self.alle.Pr.shape[1],3))
                    self.alle.forcet1       = np.zeros((self.alle.Pr.shape[0],self.alle.Pr.shape[1],3))
                    self.alle.forcet2       = np.zeros((self.alle.Pr.shape[0],self.alle.Pr.shape[1],3))

                    self.alle.fmlforcelo_norm  = np.zeros((self.alle.Pr.shape[0],self.alle.Pr.shape[1]))
                    self.alle.fmlforcet1_norm  = np.zeros((self.alle.Pr.shape[0],self.alle.Pr.shape[1]))
                    self.alle.fmlforcet2_norm  = np.zeros((self.alle.Pr.shape[0],self.alle.Pr.shape[1]))
                    self.alle.fmlforcelo       = np.zeros((self.alle.Pr.shape[0],self.alle.Pr.shape[1],3))
                    self.alle.fmlforcet1       = np.zeros((self.alle.Pr.shape[0],self.alle.Pr.shape[1],3))
                    self.alle.fmlforcet2       = np.zeros((self.alle.Pr.shape[0],self.alle.Pr.shape[1],3))

                    self.alle.fmlmtiforcelo_norm  = np.zeros((self.alle.Pr.shape[0],self.alle.Pr.shape[1]))
                    self.alle.fmlmtiforcet1_norm  = np.zeros((self.alle.Pr.shape[0],self.alle.Pr.shape[1]))
                    self.alle.fmlmtiforcet2_norm  = np.zeros((self.alle.Pr.shape[0],self.alle.Pr.shape[1]))
                    self.alle.fmlmtiforcelo       = np.zeros((self.alle.Pr.shape[0],self.alle.Pr.shape[1],3))
                    self.alle.fmlmtiforcet1       = np.zeros((self.alle.Pr.shape[0],self.alle.Pr.shape[1],3))
                    self.alle.fmlmtiforcet2       = np.zeros((self.alle.Pr.shape[0],self.alle.Pr.shape[1],3))
                    ##############################################################################
                    ##############################################################################
                    ##############################################################################
                    ##############################################################################
                    ## hier muessen wir zusehen das die werte in den jeweiligen auslenkungen abgespeichert werden
                    ##############################################################################
                    ##############################################################################
                    ##############################################################################

                    nndistances_all = np.sort(utils.unique_to_decimals(np.linalg.norm(self.alle.Pr0,axis=1),6))
                    print("nnd: (parametrize tox)",nndistances_all)
                    for dispidx,i in enumerate(self.alle.p): # geht ueber jedes displacement (xdir,quer,3nnd,...)
                        #print "direction:",self.alle.disp,"dispidx:",dispidx
                        for atomidx in np.arange(self.alle.Pr.shape[1]): # geht ueber jedes atom der vergroesserten superzelle
                            self.alle.angles[dispidx,atomidx] = utils.anglevec(self.alle.Pr0[atomidx],self.alle.Pr[dispidx,atomidx])
                            longvec0 = self.alle.Pr0[atomidx]


                            ################################################################################
                            # get shell
                            ################################################################################
                            sig_fig=5
                            a = np.linalg.norm(longvec0)
                            b = nndistances_all[1]
                            self.shell = 0
                            if a == b or int(a*10**sig_fig) == int(b*10**sig_fig):
                                self.shell = 1
                            b = nndistances_all[2]
                            if a == b or int(a*10**sig_fig) == int(b*10**sig_fig):
                                self.shell = 2
                            b = nndistances_all[3]
                            if a == b or int(a*10**sig_fig) == int(b*10**sig_fig):
                                self.shell = 3

                            if self.shell == 0:
                                continue

                            ################################################################################
                            # def get_to_ti_direction0_fcc(self):
                            ################################################################################
                            self.vec0 = longvec0
                            self.todirection0  = np.array([0.0,0.0,0.0])
                            self.tidirection0  = np.array([0.0,0.0,0.0])
                            #################################################################################
                            # shell 1
                            #################################################################################
                            if self.shell == 1:
                                if self.vec0[0] == 0.0:  # dies gilt nur fuer die ersten Nachbarn, fuer die 2NN
                                    self.todirection0 = np.array([1.0,0.0,0.0])  # minus xdir  will also be correct
                                if self.vec0[1] == 0.0:  # dies gilt nur fuer die ersten Nachbarn, fuer die 2NN
                                    self.todirection0 = np.array([0.0,1.0,0.0])  # minus xdir  will also be correct
                                if self.vec0[2] == 0.0:  # dies gilt nur fuer die ersten Nachbarn, fuer die 2NN
                                    self.todirection0 = np.array([0.0,0.0,1.0])  # minus xdir  will also be correct

                                self.tidirection0 = np.cross(self.vec0,self.todirection0)


                            #################################################################################
                            # shell 2
                            #################################################################################
                            if self.shell == 2:
                                if self.vec0[0] != 0.0:
                                    self.todirection0 = np.array([0.0,1.0,0.0])
                                    self.tidirection0 = np.array([0.0,0.0,1.0])
                                if self.vec0[1] != 0.0:
                                    self.todirection0 = np.array([1.0,0.0,0.0])
                                    self.tidirection0 = np.array([0.0,0.0,1.0])
                                if self.vec0[2] != 0.0:
                                    self.todirection0 = np.array([1.0,0.0,0.0])
                                    self.tidirection0 = np.array([0.0,1.0,0.0])

                            #################################################################################
                            # shell 3
                            #################################################################################
                            if self.shell == 3:
                                if self.vec0[0] != self.vec0[1]:
                                    self.todirection0 = np.cross(self.vec0,np.array([-self.vec0[0],-self.vec0[1],0.0]))
                                if self.vec0[0] != self.vec0[2]:
                                    self.todirection0 = np.cross(self.vec0,np.array([-self.vec0[0],0.0,-self.vec0[2]]))
                                if self.vec0[1] != self.vec0[2]:
                                    self.todirection0 = np.cross(self.vec0,np.array([0.0,-self.vec0[1],-self.vec0[2]]))

                                self.tidirection0 = np.cross(self.vec0,self.todirection0)

                            self.todirection0 = (self.todirection0/np.linalg.norm(self.todirection0))*np.linalg.norm(self.vec0)
                            self.tidirection0 = (self.tidirection0/np.linalg.norm(self.tidirection0))*np.linalg.norm(self.vec0)


                            longvec = self.alle.Pr[dispidx,atomidx]
                            self.longvec = longvec

                            self.dvec_p_to0 = utils.project_vector(self.longvec,self.todirection0) # hat aber die richtige richtung
                            self.dvec_p_ti0 = utils.project_vector(self.longvec,self.tidirection0) # hat aber die richtige richtung

                            # leave for later
                            #self.dnorm_r_to = utils.reject_vector(self.longvec,self.todirection0) # hat aber die richtige richtung (
                            #self.dnorm_r_ti = utils.reject_vector(self.longvec,self.tidirection0) # hat aber die richtige richtung (


                            self.dnorm_p_to = np.linalg.norm(self.dvec_p_to0)
                            self.dnorm_p_ti = np.linalg.norm(self.dvec_p_ti0)
                            if self.dnorm_p_to == 0.0:
                                self.dvec_p_to0_normalized = np.array([0.0,0.0,0.0])
                            else:
                                self.dvec_p_to0_normalized = self.dvec_p_to0/self.dnorm_p_to
                            if self.dnorm_p_ti == 0.0:
                                self.dvec_p_ti0_normalized = np.array([0.0,0.0,0.0])
                            else:
                                self.dvec_p_ti0_normalized = self.dvec_p_ti0/self.dnorm_p_ti


                            t11_norm, t22_norm = get_new_local_basis(\
                                    longvec = self.longvec,
                                    longvec0 = self.vec0,
                                    outdirection0 = self.todirection0)

                            ###############################################################
                            #### force VASP ###############################################
                            ###############################################################
                            #print "################",self.alle.disp,dispidx,atomidx
                            # self.alle ist entweder quer oder xdir.... schleife laeuft ueber beides
                            force   = self.alle.F[dispidx,atomidx]
                            forcelo = utils.project_vector(force,longvec)  # dies sollte die long kraft sein
                            forcet1 = utils.project_vector(force,t11_norm)
                            forcet2 = utils.project_vector(force,t22_norm)


                            # hier muesste man eigentlich auch noch auf die Vorzeichen achten;
                            # fuer den longitudinalen kann man einfach checken ob dieser parallel oder antiparallel zu longvec ist;
                            # fuer die transversalen muss man eigenlicht auch einfach nur schauen ob diese parallel zu dnorm oder antiparalle zu dnorm wirken;
                            # in principle one can define a shpere which two half spheres seperated by lonvec0; half spehre which is parallel to dnorm: posivive; other = neg
                            veclon = utils.anglevec(forcelo,longvec)
                            veclonvz = 1.0
                            if abs(veclon) > 1.57:
                                veclonvz = -1.0

                            vect1n = utils.anglevec(forcet1,t11_norm)
                            vect1nvz = 1.0
                            if abs(vect1n) > 1.57:
                                vect1nvz = -1.0

                            vect2n = utils.anglevec(forcet2,t22_norm)
                            vect2nvz = 1.0
                            if abs(vect2n) > 1.57:
                                vect2nvz = -1.0

                            forcelo_norm = veclonvz * np.linalg.norm(forcelo)
                            forcet1_norm = vect1nvz * np.linalg.norm(forcet1)
                            forcet2_norm = vect2nvz * np.linalg.norm(forcet2)

                            ###############################################################
                            #### force fml ################################################
                            ###############################################################
                            fmlforce        = self.alle.FML[dispidx,atomidx]
                            fmlforcelo      = utils.project_vector(fmlforce,longvec)
                            fmlforcet1      = utils.project_vector(fmlforce,t11_norm)
                            fmlforcet2      = utils.project_vector(fmlforce,t22_norm)

                            fmlveclon = utils.anglevec(fmlforcelo,longvec)
                            fmlveclonvz = 1.0
                            if abs(fmlveclon) > 1.57:
                                fmlveclonvz = -1.0

                            fmlvect1n = utils.anglevec(fmlforcet1,t11_norm)
                            fmlvect1nvz = 1.0
                            if abs(fmlvect1n) > 1.57:
                                fmlvect1nvz = -1.0

                            fmlvect2n = utils.anglevec(fmlforcet2,t22_norm)
                            fmlvect2nvz = 1.0
                            if abs(fmlvect2n) > 1.57:
                                fmlvect2nvz = -1.0

                            fmlforcelo_norm = fmlveclonvz * np.linalg.norm(fmlforcelo)
                            fmlforcet1_norm = fmlvect1nvz * np.linalg.norm(fmlforcet1)
                            fmlforcet2_norm = fmlvect2nvz * np.linalg.norm(fmlforcet2)

                            ###############################################################
                            #### force fmlmti #############################################
                            ###############################################################
                            #print 'kkk;',self.alle.fmlmti.shape
                            fmlmtiforce        = self.alle.FMLmti[dispidx,atomidx]
                            fmlmtiforcelo      = utils.project_vector(fmlforce,longvec)
                            fmlmtiforcet1      = utils.project_vector(fmlforce,t11_norm)
                            fmlmtiforcet2      = utils.project_vector(fmlforce,t22_norm)

                            fmlmtiveclon = utils.anglevec(fmlmtiforcelo,longvec)
                            fmlmtiveclonvz = 1.0
                            if abs(fmlmtiveclon) > 1.57:
                                fmlmtiveclonvz = -1.0

                            fmlmtivect1n = utils.anglevec(fmlmtiforcet1,t11_norm)
                            fmlmtivect1nvz = 1.0
                            if abs(fmlmtivect1n) > 1.57:
                                fmlmtivect1nvz = -1.0

                            fmlmtivect2n = utils.anglevec(fmlmtiforcet2,t22_norm)
                            fmlmtivect2nvz = 1.0
                            if abs(fmlmtivect2n) > 1.57:
                                fmlmtivect2nvz = -1.0

                            fmlmtiforcelo_norm = fmlmtiveclonvz * np.linalg.norm(fmlmtiforcelo)
                            fmlmtiforcet1_norm = fmlmtivect1nvz * np.linalg.norm(fmlmtiforcet1)
                            fmlmtiforcet2_norm = fmlmtivect2nvz * np.linalg.norm(fmlmtiforcet2)


                            if self.shell == 1:
                                #print atomidx,dispidx,"longvec0:",self.vec0,longvec,"forceVASP:",force,"forcelo:",forcelo,"ft1:",forcet1,"ft2:",forcet2,"CHECK:",((force-forcet1)-forcet2)-forcelo
                                #print atomidx,dispidx,np.linalg.norm(longvec),"forceVASP:",force,"forcelo:",forcelo,"ft1:",forcet1,"ft2:",forcet2,"CHECK:",((force-forcet1)-forcet2)-forcelo
                                pass
                            ##############################################################################
                            # always make the check:
                            ##############################################################################
                            def printinfo1():
                                print("displacement:",self.alle.P[dispidx,0])
                                print("self.dvec_p_to0_normalized:",self.dvec_p_to0_normalized)
                                print("self.dvec_p_ti0_normalized:",self.dvec_p_ti0_normalized)
                                print("")
                                print("------")
                                print("force(VASP==total):",force)
                                print("longforce         :",force - fmlforce,"norm:", np.linalg.norm(force-fmlforce))
                                print("fmlforce = f - l  :",fmlforce)
                                print("longvec:",longvec,"t11_norm:",t11_norm,"t22_norm:",t22_norm)
                                print("longvecnorm",np.linalg.norm(longvec),"0:",nndistances_all[1],"diff:",np.linalg.norm(longvec)-nndistances_all[1])
                                print("------")
                                print("")
                                print("forcet1:",forcet1)
                                print("forcet2:",forcet2)
                                print("forcelo:",forcelo)
                                print("fmlforcet1:",fmlforcet1)
                                print("fmlforcet2:",fmlforcet2)
                                print("fmlforcelo:",fmlforcelo)
                                print()
                                print("forcet1n:",forcet1_norm)
                                print("forcet2n:",forcet2_norm)
                                print("forcelon:",forcelo_norm)
                                print("fmlforcet1n:",fmlforcet1_norm)
                                print("fmlforcet2n:",fmlforcet2_norm)
                                print("fmlforcelon:",fmlforcelo_norm)
                                print()
                                print("self.shell:",self.shell)
                                print(self.alle.disp)
                                print("   CHECK:",((force-forcet1)-forcet2)-forcelo)
                                print("FMLCHECK:",((fmlforce-fmlforcet1)-fmlforcet2)-fmlforcelo)
                                sys.exit()

                            check = ((force-forcet1)-forcet2)-forcelo
                            if abs(np.linalg.norm(check)) > 0.0001:
                                sys.exit("check")

                            if False and atomidx == 60 and dispidx == 37 and self.alle.disp == 'quer':
                                printinfo1()
                            if False and atomidx == 60 and dispidx == 37 and self.alle.disp == 'xdir':
                                printinfo1()
                            if False and atomidx == 60 and dispidx == 53 and self.alle.disp == 'quer':
                                printinfo1()
                            if False and atomidx == 60 and dispidx == 53 and self.alle.disp == 'xdir':
                                printinfo1()
                            if False and atomidx == 60 and dispidx == 92 and self.alle.disp == 'xdir':
                                printinfo1()
                            if False and atomidx == 60 and dispidx == 96 and self.alle.disp == 'quer': # [0.96, 0, 0]
                                printinfo1()
                            if False and atomidx == 8 and dispidx == 160 and self.alle.disp == 'xdir':
                                printinfo1()

                            # atomidx = 60 oder atomidx = 28 == atom bei [ -1, -1, 0 ]
                            # dispidx = 27 == displacement in quer of [0.26,0.26,0]
                            self.alle.angles[dispidx,atomidx] = utils.anglevec(self.alle.Pr0[atomidx],self.alle.Pr[dispidx,atomidx])

                            self.alle.dnorm_p_to[dispidx,atomidx] = self.dnorm_p_to
                            self.alle.dnorm_p_ti[dispidx,atomidx] = self.dnorm_p_ti
                            self.alle.dvec_p_to0[dispidx,atomidx] = self.dvec_p_to0
                            self.alle.dvec_p_ti0[dispidx,atomidx] = self.dvec_p_ti0

                            self.alle.forcelo[dispidx,atomidx] = forcelo
                            self.alle.forcet1[dispidx,atomidx] = forcet1
                            self.alle.forcet2[dispidx,atomidx] = forcet2
                            self.alle.forcelo_norm[dispidx,atomidx] = forcelo_norm
                            self.alle.forcet1_norm[dispidx,atomidx] = forcet1_norm
                            self.alle.forcet2_norm[dispidx,atomidx] = forcet2_norm

                            self.alle.fmlforcelo[dispidx,atomidx]         = fmlforcelo
                            self.alle.fmlforcet1[dispidx,atomidx]         = fmlforcet1
                            self.alle.fmlforcet2[dispidx,atomidx]         = fmlforcet2
                            self.alle.fmlforcelo_norm[dispidx,atomidx]    = fmlforcelo_norm
                            self.alle.fmlforcet1_norm[dispidx,atomidx]    = fmlforcet1_norm
                            self.alle.fmlforcet2_norm[dispidx,atomidx]    = fmlforcet2_norm

                            self.alle.fmlmtiforcelo[dispidx,atomidx]         = fmlmtiforcelo
                            self.alle.fmlmtiforcet1[dispidx,atomidx]         = fmlmtiforcet1
                            self.alle.fmlmtiforcet2[dispidx,atomidx]         = fmlmtiforcet2
                            self.alle.fmlmtiforcelo_norm[dispidx,atomidx]    = fmlmtiforcelo_norm
                            self.alle.fmlmtiforcet1_norm[dispidx,atomidx]    = fmlmtiforcet1_norm
                            self.alle.fmlmtiforcet2_norm[dispidx,atomidx]    = fmlmtiforcet2_norm


                # if shifted one tab to the right will not be parametrized every time
                # new tox:
                # plt.plot(pot.disp.xdir.dnorm_p_to[:,8],pot.disp.xdir.forcet1_norm[:,8])
                # plt.plot(pot.disp.xdir.dnorm_p_to[:,8],pot.disp.xdir.fmlforcet1_norm[:,8])
                # plt.plot(pot.disp.xdir.dnorm_p_to[:,8],pot.disp.xdir.fmlforcelo_norm[:,8])
                #
                # new tix:
                # plt.plot(pot.disp.quer.dnorm_p_ti[:,60],pot.disp.quer.forcet2_norm[:,60])
                # plt.plot(pot.disp.quer.dnorm_p_ti[:,60],pot.disp.quer.forcelo_norm[:,60])
                # plt.plot(pot.disp.quer.dnorm_p_ti[:,60],pot.disp.quer.fmlforcelo_norm[:,60])
                ###########################################################################
                # dies hier wird gemacht immer wenn parametrize_tox gestartet wurde
                ###########################################################################
        #for folderidx,folder in enumerate(go_through_folder):
        #    print "folder:",folder
        #    #sys.exit()
        #    for a in self.alatsstringlist:

                if self.alle.disp == 'xdir':
                    print(utils.printyellow("folder:",folder))
                    to  = np.transpose([self.alle.dnorm_p_to[:,8],self.alle.fmlforcet1_norm[:,8]])
                    tolc = np.transpose([self.alle.dnorm_p_to[:,8],self.alle.fmlforcelo_norm[:,8]])

                    to_ = np.transpose([-self.alle.dnorm_p_to[:,8],-self.alle.fmlforcet1_norm[:,8]])
                    to = np.concatenate((to,to_))

                    tolc_ = np.transpose([-self.alle.dnorm_p_to[:,8],-self.alle.fmlforcelo_norm[:,8]])
                    tolc = np.concatenate((tolc,tolc_))

                    tody = np.abs(to[:,1]) - np.abs(tolc[:,1])
                    tod = np.transpose([to[:,0],tody])
                    tod = utils.remove_duplicates_of_2d_array_within_tolerance(tod)
                    tod = tod[tod[:,0].argsort()]
                    np.savetxt(folder+"/tod.dat",tod,fmt="%.7f %.7f")

                    to = utils.remove_duplicates_of_2d_array_within_tolerance(to)
                    to = to[to[:,0].argsort()]
                    tolc = utils.remove_duplicates_of_2d_array_within_tolerance(tolc)
                    tolc = tolc[tolc[:,0].argsort()]

                    np.savetxt(folder+"/to.dat",to,fmt="%.7f %.7f")
                    np.savetxt(folder+"/tolcorr.dat",tolc,fmt="%.7f %.7f")



                    ti = np.transpose([self.alle.dnorm_p_ti[:,60],self.alle.fmlforcet1_norm[:,60]])
                    tilc = np.transpose([self.alle.dnorm_p_ti[:,60],self.alle.fmlforcelo_norm[:,60]])
                    ti_ = np.transpose([-self.alle.dnorm_p_ti[:,60],-self.alle.fmlforcet1_norm[:,60]])
                    ti = np.concatenate((ti,ti_))
                    tilc_ = np.transpose([-self.alle.dnorm_p_ti[:,60],-self.alle.fmlforcelo_norm[:,60]])
                    tilc = np.concatenate((tilc,tilc_))
                    tidy = np.abs(ti[:,1]) - np.abs(tilc[:,1])
                    tid = np.transpose([ti[:,0],tidy])
                    tid = utils.remove_duplicates_of_2d_array_within_tolerance(tid)
                    tid = tid[tid[:,0].argsort()]
                    np.savetxt(folder+"/tid.dat",tid,fmt="%.7f %.7f")

                    ti = utils.remove_duplicates_of_2d_array(ti)
                    ti = ti[ti[:,0].argsort()]
                    tilc = utils.remove_duplicates_of_2d_array_within_tolerance(tilc)
                    tilc = tilc[tilc[:,0].argsort()]

                    np.savetxt(folder+"/ti60.dat",ti,fmt="%.7f %.7f")
                    np.savetxt(folder+"/ti60lcorr.dat",tilc,fmt="%.7f %.7f")



                if self.alle.disp == 'quer':
                    print("self.alle.atoms:",self.alle.atoms)
                    if self.alle.atoms == 32:
                        queratom = 60
                    if self.alle.atoms == 108:
                        queratom = 99
                    print(utils.printyellow("self.alle.disp quer, folder:",folder))

                    ti = np.transpose([self.alle.dnorm_p_ti[:,queratom],self.alle.fmlforcet2_norm[:,queratom]])
                    print("ti:!!!: ",ti)
                    print("rest wird ausgegraut da gerade nicht funktionierend in /Users/glensk/Understand_distributions//displacements_dense/Ir//3x3x3sc_3.99Ang_quer_2x2x2kp/ folder")
                    #tilc = np.transpose([self.alle.dnorm_p_ti[:,queratom],self.alle.fmlforcelo_norm[:,queratom]])

                    #ti_ = np.transpose([-self.alle.dnorm_p_ti[:,queratom],-self.alle.fmlforcet2_norm[:,queratom]])
                    #ti = np.concatenate((ti,ti_))

                    #tilc_ = np.transpose([-self.alle.dnorm_p_ti[:,queratom],-self.alle.fmlforcelo_norm[:,queratom]])
                    #tilc = np.concatenate((tilc,tilc_))

                    ## newparam mti
                    ##timm = np.transpose([self.alle.dnorm_p_ti[:,60],self.alle.fmlmtiforcet2_norm[:,60]])
                    ##timmlc = np.transpose([self.alle.dnorm_p_ti[:,60],self.alle.fmlforcelo_norm[:,60]])
                    ##tid = np.transpose([self.alle.dnorm_p_ti[:,60],np.abs(ti[:,1]) - np.abs(tilc[:,1])])
                    ##np.savetxt(folder+"/tid.dat",tid,fmt="%.7f %.7f")

                    #tidy = np.abs(ti[:,1]) - np.abs(tilc[:,1])
                    #tid = np.transpose([ti[:,0],tidy])
                    #tid = utils.remove_duplicates_of_2d_array_within_tolerance(tid)
                    #tid = tid[tid[:,0].argsort()]
                    #np.savetxt(folder+"/tid.dat",tid,fmt="%.7f %.7f")

                    #ti = utils.remove_duplicates_of_2d_array(ti)
                    #ti = ti[ti[:,0].argsort()]
                    #tilc = utils.remove_duplicates_of_2d_array_within_tolerance(tilc)
                    #tilc = tilc[tilc[:,0].argsort()]

                    #print utils.printgreen("save ",folder,"ti.dat")
                    #np.savetxt(folder+"/ti.dat",ti,fmt="%.7f %.7f")
                    #np.savetxt(folder+"/tilcorr.dat",tilc,fmt="%.7f %.7f")


                ##############################################################
                # now parametrize
                ##############################################################
                if self.alle.disp == 'xdir':
                    print("parametrize to")
                    for orderfit in np.arange(1,10):
                        # makes something like lon_1nn_pos_poly9_6_fit___157.32847059439_-253.3149803360_172.60763615111_-64.55258342696_14.303612909803_-1.874831186120_0.1343007378748_-0.004043636776 BUT NOT A CORRESPONDING DELTAFOLDER
                        # here the deltas are saved in deltas but no fit is made afterwards
                        coefsrealene,coefsrealforces,deltas = polyfit(
                            foldername=folder,
                            filename="to_1nn_all",
                            zeroat=0.0,
                            data=to,
                            ordermax=9,
                            verbose2=False
                            )
                        self.alle.coefs_ene_to = coefsrealene
                    print("coefsrealene to:",coefsrealene)

                if self.alle.disp == 'quer':
                    print("parametrize ti")
                    for orderfit in np.arange(1,10):
                        # makes something like lon_1nn_pos_poly9_6_fit___157.32847059439_-253.3149803360_172.60763615111_-64.55258342696_14.303612909803_-1.874831186120_0.1343007378748_-0.004043636776 BUT NOT A CORRESPONDING DELTAFOLDER
                        # here the deltas are saved in deltas but no fit is made afterwards
                        print("data ti xxx:",ti)
                        coefsrealene,coefsrealforces,deltas = polyfit(
                            foldername=folder,
                            filename="ti_1nn_all",
                            zeroat=0.0,
                            data=ti,
                            ordermax=9,
                            verbose2=False
                            )
                        self.alle.coefs_ene_ti = coefsrealene
                    print("coefsrealene ti:",coefsrealene)
                        #print "coefsrealforces:",coefsrealforces
                    #print "parametrize ti"

                #print '-->typeti',self.alle.coefs_ene_ti,type(self.alle.coefs_ene_ti)
                #print '-->typeto',self.alle.coefs_ene_to,type(self.alle.coefs_ene_to)
                # puts self.alle to self.xdir or self.quer ...
                # afterwrs self.xdir or self.quer should be defined
                self.update_what_to_calculate_and_write_self_xdir_quer()



                ## after loop save to TOTI.npz file
                #self.alle.saveTOTInpz(folder,"TOTI")


                ## puts self.alle to self.xdir or self.quer ...
                #self.update_what_to_calculate_and_write_self_xdir_quer()
                #print '+++-->typeti',self.alle.coefs_ene_ti,type(self.alle.coefs_ene_ti)
                #print '+++-->typeto',self.alle.coefs_ene_to,type(self.alle.coefs_ene_to)
                # after loop save to TOTI.npz file
                self.alle.saveTOTInpz(folder,"TOTI")

        print(utils.printblue("     pot_parametrize.parametrize_tox complete ..."))
        return

    def get_FMLnpz(self, shell = False, to = False, ti = False):
        ''' this can only be performed once the .pkl file exists

        subtracts from the forces on every atom the parametrized longforces

        by substracting quer displacement from xdir
        (creates parametrization for tix vectors --> NO, it does not do that)

        concrete: go through every displacement in xdirection and substract morse + correction from force on , the ramaining forces are the ti{x,y,z} contributions
        the error will be as high as in the diffshiftedpolydiff/ == self.fit[1,'lon','xdir','dos','morse'].ix['polybesteiff']
        im haupt *.pkl file: /Users/glensk/Understand_distributions/displacements_dense/Ir/2x2x2sc_data_3x3x3kp.pkl
             '''
        print("")
        print(utils.printblue("     get_FMLnpz (FML) ..."))
        # xdir has usually all the data whil quer has only quer data
        from . import pot_energy_forces

        ##################################################################################
        # go through every displacement and get long forces in xyz direction
        ##################################################################################
        self.update_what_to_calculate_and_write_self_xdir_quer(onlyprint = True)
        self.alle.loaddatapkl(self.pkl)
        for folderidx,folder in enumerate(self.folder_displacement_direction_all):
            for a in self.alatsstringlist:
                print("xxxfolder:",folder)
        #sys.exit('so weit')
        for folderidx,folder in enumerate(self.folder_displacement_direction_all):
            for a in self.alatsstringlist:
                #print "disp:",self.alle.disp
                #print "direction:",direction
                #print "coefsto:",paramsto
                #print "coefsti:",paramsti
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

                self.pot_energy_forces.params.u1nn_tox = to
                self.pot_energy_forces.params.u1nn_tix = ti

                print("self.pot_energy_forces.params.u1nn_tox:",self.pot_energy_forces.params.u1nn_tox)
                print("self.pot_energy_forces.params.u1nn_tix:",self.pot_energy_forces.params.u1nn_tix)

                self.pot_energy_forces.pot_parametrize = copy.deepcopy(self.alle)

                def check_if_parametrizations_long_same():
                    if self.alle.fml_u1nn_pottype == self.pot_energy_forces.params.u1nn_pottype and \
                       np.array_equal(self.alle.fml_u1nn_potparam, self.pot_energy_forces.params.u1nn_potparam) == True and \
                       self.alle.fml_u1nn_potadd == self.pot_energy_forces.params.u1nn_potadd and \
                       np.array_equal(self.alle.fml_u1nn_potaddparam, self.pot_energy_forces.params.u1nn_potaddparam) == True and \
                       np.array_equal(self.alle.fml_u1nn_tox, self.pot_energy_forces.params.u1nn_tox) == True and \
                       np.array_equal(self.alle.fml_u1nn_tix, self.pot_energy_forces.params.u1nn_tix) == True and \
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
                self.alle.creatennforcesnpz(folder = folder, a = a)  # loads npz (forces)
                self.alle.infodos()
                needtoparametrize = None  # see further up
                print("needtoparametrize1:",needtoparametrize)
                if os.path.isfile(self.alle.parametrize_fml_npz_path) != True or self.parametrize_fml_npz == True:
                    print(utils.printred("1--------------------------------------------------------"))
                    print("os.path.isfile(self.alle.parametrize_fml_npz_path)",os.path.isfile(self.alle.parametrize_fml_npz_path))
                    print("self.parametrize_fml_npz:",self.parametrize_fml_npz)
                    needtoparametrize = True
                    print("needtoparametrize2:",needtoparametrize)
                    print(utils.printred("1--------------------------------------------------------"))
                print("needtoparametrize2:",needtoparametrize)
                # if FML file exists we, in any case, want to check if FML has to be redone
                if needtoparametrize == False:
                    if os.path.isfile(self.alle.parametrize_fml_npz_path) == True: # and self.parametrize_fml_npz == False:
                        self.alle.loadFMLnpz(folder,"FML")
                        if type(self.alle.FML) != bool:
                            if check_if_parametrizations_long_same() == True:
                                needtoparametrize = False
                                print(utils.printgreen("    old parameters from FML.npz and new ones are the same; no need for parametrization!; loaded FML.npz"))
                            else:
                                needtoparametrize = True
                                print(utils.printred("----------- compare parametrizations; they seem to differ -----------"))
                                print(self.alle.fml_u1nn_pottype,"-----",self.pot_energy_forces.params.u1nn_pottype)
                                print(self.alle.fml_u1nn_potparam,"-----",self.pot_energy_forces.params.u1nn_potparam)
                                print(self.alle.fml_u1nn_potadd,"-----",self.pot_energy_forces.params.u1nn_potadd)
                                print(self.alle.fml_u1nn_potaddparam,"-----",self.pot_energy_forces.params.u1nn_potaddparam)
                                print(self.alle.fml_u1nn_tox,"-----",self.pot_energy_forces.params.u1nn_tox)
                                print(self.alle.fml_u1nn_tix,"-----",self.pot_energy_forces.params.u1nn_tix)
                        else:
                            needtoparametrize = True
                            print(utils.printred("self.alle.FML is of type bool ! , need to redo ..."))
                print("needtoparametrize3:",needtoparametrize)
                ##########################################################################
                # in case it is neccesary to parametrize
                ##########################################################################
                if needtoparametrize == False:
                    print(utils.printblue("         get_FMLnpz complete ..."))
                    return

                #if needtoparametrize == True:
                self.pot_energy_forces.load_parameters_from_file = False
                self.pot_energy_forces.cell = self.alle.sc
                self.pot_energy_forces.coord0_rrel = self.alle.eqcoords
                self.alle.fml    = np.zeros((self.alle.p.shape[0],self.alle.p.shape[1],3))
                self.alle.fmlmti = np.zeros((self.alle.p.shape[0],self.alle.p.shape[1],3))
                self.alle.fmlmto = np.zeros((self.alle.p.shape[0],self.alle.p.shape[1],3))


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
                self.alle.FML    = np.zeros((self.alle.p.shape[0],self.alle.p.shape[1]*nsc**3,3))
                self.alle.FMLmti = np.zeros((self.alle.p.shape[0],self.alle.p.shape[1]*nsc**3,3))
                self.alle.FMLmto = np.zeros((self.alle.p.shape[0],self.alle.p.shape[1]*nsc**3,3))


                ##########################################################################
                # parametrization
                ##########################################################################
                self.alle.fml_u1nn_pottype = 'morse'
                self.alle.fml_u1nn_potparam = self.alle.data[1,'lon','quer','dos'].ix['morse']
                #self.alle.fml_u1nn_potadd = 'poly'
                #self.alle.fml_u1nn_potaddparam = self.alle.fit[1,'lon','quer','dos','morse'].ix['polybeste']
                self.alle.fml_u1nn_potadd = 'spline'
                #self.alle.fml_u1nn_potaddparam = self.alle.fit[1,'lon','quer','alle','morse'].ix['splinedata']
                #self.alle.fml_u1nn_potaddparam = 'lon.quer.alle.morse.splinedata'
                self.alle.fml_u1nn_potaddparam = 'lon.quer.dos.morse.splinedata'

                print("self.alle.fml_u1nn_tox        :",to)
                print("self.alle.fml_u1nn_tix        :",ti)
                print("---------")
                print("self.pot_energy_forces.params.u1nn_tix:",self.pot_energy_forces.params.u1nn_tix)
                print("self.pot_energy_forces.params.u1nn_tox:",self.pot_energy_forces.params.u1nn_tox)

                substrall = [ 'lon', 'lonto', 'lonti' ]
                if to == False and ti == False:
                    substrall = [ 'lon' ]

                for substr in substrall:
                    if substr == 'lon':
                        self.alle.fml_u1nn_tox = False
                        self.alle.fml_u1nn_tix = False
                        self.pot_energy_forces.params.u1nn_tipx = False
                        self.pot_energy_forces.params.u1nn_topx = False
                    elif substr == 'lonti':
                        self.alle.fml_u1nn_tix = ti
                        self.alle.fml_u1nn_tox = False
                        self.pot_energy_forces.params.u1nn_tipx = ti
                        self.pot_energy_forces.params.u1nn_topx = False
                    elif substr == 'lonto':
                        self.alle.fml_u1nn_tix = False
                        self.alle.fml_u1nn_tox = to
                        self.pot_energy_forces.params.u1nn_tipx = False
                        self.pot_energy_forces.params.u1nn_topx = to
                    else:
                        sys.exit("substr not known")
                    print("###############################################################")
                    print(folder)
                    print(a)
                    print("self.alle.fml_u1nn_pottype    :",self.alle.fml_u1nn_pottype)
                    print("self.alle.fml_u1nn_potparam   :",self.alle.fml_u1nn_potparam)
                    print("self.alle.fml_u1nn_potadd     :",self.alle.fml_u1nn_potadd)
                    print("self.alle.fml_u1nn_potaddparam:",self.alle.fml_u1nn_potaddparam)
                    print("self.alle.fml_u1nn_tox        :",self.alle.fml_u1nn_tox)
                    print("self.alle.fml_u1nn_tix        :",self.alle.fml_u1nn_tix)
                    print("self.pot_energy_forces.params.u1nn_tix:",self.pot_energy_forces.params.u1nn_tipx)
                    print("self.pot_energy_forces.params.u1nn_tox:",self.pot_energy_forces.params.u1nn_topx)
                    print("###############################################################")

                    ###############################################################################
                    ################# we go through all positions from the npz file ###############
                    ###############################################################################
                    if True:
                        for idx,i in enumerate(self.alle.p):                                    # we go through all displacements  quer,xdir,...
                            #self.pot_energy_forces.verbose = 2
                            self.pot_energy_forces.coord_cart = i                               # we need the cartesian coordinates
                            self.pot_energy_forces.pot_run()

                            # HERE WE DO THE SUBTRACTION
                            if substr == 'lon':
                                self.alle.fml[idx] = self.alle.f[idx]-self.pot_energy_forces.forces # we need the vasp forces
                                forces = np.copy(self.alle.fml[idx])
                            elif substr == 'lonti':
                                #print "kkk:|,", self.alle.f[idx]-self.pot_energy_forces.forces # we need the vasp forces
                                self.alle.fmlmti[idx] = self.alle.f[idx]-self.pot_energy_forces.forces # we need the vasp forces
                                forces = np.copy(self.alle.fmlmti[idx])
                            elif substr == 'lonto':
                                self.alle.fmlmto[idx] = self.alle.f[idx]-self.pot_energy_forces.forces # we need the vasp forces
                                forces = np.copy(self.alle.fmlmto[idx])
                            else:
                                sys.exit("substr not known")

                            # for FML:
                            self.crystalforces = crystal_generator.crystal()
                            self.crystalforces.load_positions_cell(coord_cart = forces , cell = self.alle.sc)

                            # forces big (repeat supercell)
                            self.crystalforcesbig = crystal_generator.supercell()
                            self.crystalforcesbig.create_supercell(  self.crystalforces, nsc, nsc, nsc, newsorting = True, shiftforces = True )
                            if substr == 'lon':
                                self.alle.FML[idx]    = self.crystalforcesbig.rcar
                            if substr == 'lonti':
                                self.alle.FMLmti[idx] = self.crystalforcesbig.rcar
                            if substr == 'lonto':
                                self.alle.FMLmto[idx] = self.crystalforcesbig.rcar



                            ################################################
                            # print to screen 183/199
                            ################################################
                            #if type(self.alle.dstring) == list:
                            #    print idx,"/",len(self.alle.dstring),self.alle.FML[37,60],self.alle.FMLmti[37,60],self.alle.FMLmto[37,60]
                            #else:
                            #    print idx,"/",self.alle.dstring.shape[0],self.alle.FML[37,60],self.alle.FMLmti[37,60],self.alle.FMLmto[37,60]
                            #if idx == 37:
                            #    print "tix:",self.pot_energy_forces.params.u1nn_tipx
                            #    print "tox:",self.pot_energy_forces.params.u1nn_topx
                            #    print "substr:",substr
                            #    print "---- forces --- ()"
                            #    print forces
                            #    print '--- frceslong ---'
                            #    print self.pot_energy_forces.forceslong
                            #    print '--- frcestrantopx ---'
                            #    print self.pot_energy_forces.forcestrantopx
                            #    print '--- frcestrantipx ---'
                            #    print self.pot_energy_forces.forcestrantipx


                    # hoehe schleife 'lon', 'lonto', 'lonti'
                    # puts self.alle to self.xdir or self.quer ...
                    self.update_what_to_calculate_and_write_self_xdir_quer()

                    #if needtoparametrize == True:
                    self.alle.saveFMLnpz(folder,"FML")


                # hoehe schleife 'lon', 'lonto', 'lonti'
                # puts self.alle to self.xdir or self.quer ...
                self.update_what_to_calculate_and_write_self_xdir_quer()

                #if needtoparametrize == True:
                self.alle.saveFMLnpz(folder,"FML")




                ###############################################################################
                ################# we go through a POSITIONs file  #############################
                ###############################################################################
                if False:
                    if False:
                        out = np.array([[999,999,999]])
                        def getpos_from_POSITIONs(schritte,pos,atoms):
                            ''' returnes from POITIONs file (pos) the only coordinates of a certain given step '''
                            if type(schritte) == np.ndarray:
                                for i in schritte:
                                    print(pos[atoms*schritte:atoms*schritte+atoms])
                            return pos[atoms*schritte:atoms*schritte+atoms]

                        POSITIONs = np.loadtxt(self.folder_base+"/displacements_dense/Ir/abc_2x2x2sc_qubus_3x3x3kp/POSITIONs")
                        print("POSITIONs.shape[0]:",POSITIONs.shape[0])
                        schritte  = POSITIONs.shape[0]/32
                        print("schritte:",schritte)
                        print()
                        out = np.loadtxt("FML")
                        for i in np.arange(schritte):
                            print(i)
                            if i < 5011:
                                continue
                            posforceaktuell = getpos_from_POSITIONs(i,POSITIONs,self.alle.atoms)
                            self.pot_energy_forces.coord_cart = posforceaktuell[:,:3]
                            fDFT = posforceaktuell[:,3:]

                            self.pot_energy_forces.pot_run()
                            fml = fDFT - self.pot_energy_forces.forces
                            out = np.append(out,fml)
                            out = out.reshape(out.shape[0]/3,3)

                            if i in np.arange(5000,6000,100):
                                print("PWD:",os.getcwd())
                                np.savetxt("FML",out,fmt="%.7f %.7f %.7f")

                        #posforceaktuell = self.getpos_from_POSITIONs(i,pos,atoms)
                        #posaktuell = posforceaktuell[:,:3]  # those are the positions generated by the run
                        #fDFT = posforceaktuell[:,3:]

                        out = out[1:]
                        np.savetxt("FML",out,fmt="%.7f %.7f %.7f")


        print(utils.printblue("         get_FMLnpz complete ..."))
        print("")
        return

    def analyze_restforces(self):
        ''' at this point we have everything loaded: self.xdir and self.quer '''
        print("")
        print("")
        print(utils.printpink("     analyze_restforces ..."))
        # bis hierher wird nix geschrieben wenn nicht -pnpz -ppkldata -pnpzfml


        for folderidx,folder in enumerate(self.folder_displacement_direction_all):
            print("folder:",folder)
            for a in self.alatsstringlist:
                print("a:",a)
                if self.alle.parametrize_fml_npz == False and os.path.isfile(folder+"/FML.npz") == True:
                    self.alle.loadFMLnpz(folder,"FML")
                else:
                    self.alle.loadnnforcesnpz(     folder, 'nnforces' )  # loads npz (forces)
                    self.alle.loadFMLnpz(  folder, 'FML'  )
                    print("     direction:",self.alle.disp)
                    self.alle.angles = np.zeros((self.alle.Pr.shape[0],self.alle.Pr.shape[1]))
                    self.alle.fprojlongold = np.zeros((self.alle.Pr.shape[0],self.alle.Pr.shape[1],3))
                    self.alle.fprojlongnew = np.zeros((self.alle.Pr.shape[0],self.alle.Pr.shape[1],3))
                    for dispidx,i in enumerate(self.alle.p): # geht ueber jedes displacement
                        print("dispidx:",dispidx)
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

                    sys.exit("after that saves stuff")
                    self.update_what_to_calculate_and_write_self_xdir_quer()

                    # here we save fprojlongold (in case of qubus this is 15 MB)
                    self.alle.saveFMLnpz(folder,"FML")
        if True:
            if True:
                if True:
                    print('yes')
                    ##########################################################################
                    # analysis tix
                    ##########################################################################
                    disp = [ 'xdir', 'quer' ]
                    atomslon = [ 24, 60 ]  # only for 2x2x2 cell?

                    def save_analyze(disp,atom,force_in = False, xfunc = False):
                        ''' save_analyze '''
                        if xfunc not in [ 'winkel', 'longvec', 'dnorm' ]:
                            sys.exit("xfunc !! really not in winkel / longvec")

                        if disp not in self.disp_all:
                            sys.exit("disp !! really not in xdir quer")
                        if force_in not in [ 'xyz', 'quo', 'qun', 'xyzorig' ]:
                            print("force_in:",force_in)
                            sys.exit("force_in !! really not in xyz, quo, qun")
                        if xfunc == 'winkel':
                            try:
                                if disp == 'xdir': xfunctake = self.xdir.angles[:,atom]
                                if disp == 'quer': xfunctake = self.quer.angles[:,atom]
                                if disp == 'midd': xfunctake = self.midd.angles[:,atom]
                                if disp == 'qubus': xfunctake = self.qubus.angles[:,atom]
                            except TypeError:
                                print("ERROR 1:",force_in,xfunc,disp,atom)
                                return
                        if xfunc == 'dnorm':
                            if disp == 'xdir': xfunctake = self.xdir.dnorm
                            if disp == 'quer': xfunctake = self.quer.dnorm
                            if disp == 'midd': xfunctake = self.midd.dnorm
                            if disp == 'qubus': xfunctake = self.qubus.dnorm
                        if xfunc == 'longvec':
                            if disp == 'xdir': xfunctake = np.linalg.norm(self.xdir.Pr[:,atom],axis=1)
                            if disp == 'quer': xfunctake = np.linalg.norm(self.quer.Pr[:,atom],axis=1)
                            if disp == 'midd': xfunctake = np.linalg.norm(self.midd.Pr[:,atom],axis=1)
                            if disp == 'qubus': xfunctake = np.linalg.norm(self.qubus.Pr[:,atom],axis=1)
                            pass

                        if force_in == 'xyz':
                            getfrom = 'FML'
                        if force_in == 'xyzorig':
                            getfrom = 'F'
                        if force_in == 'quo':
                            getfrom = 'fprojlongold'
                        if force_in == 'qun':
                            getfrom = 'fprojlongnew'

                        if True:
                            print("disp:",disp)
                            print("getfrom:",getfrom)
                            print("atom:",atom)

                        folder = self.folder_base+'/analyze/'+str(self.alle.element)+'/'
                        if os.path.isdir(folder) != True:
                            os.makedirs(folder)

                        if type(eval('self.'+disp+'.'+getfrom)) == bool:
                            print("disp:",disp)
                            print("getfrom:",getfrom)
                            sys.exit('this is bool')
                        forcesx = eval('self.'+disp+'.'+getfrom+'[:,'+str(atom)+'][:,0]')
                        forcesy = eval('self.'+disp+'.'+getfrom+'[:,'+str(atom)+'][:,1]')
                        forcesz = eval('self.'+disp+'.'+getfrom+'[:,'+str(atom)+'][:,2]')

                        filenamex = folder+force_in+"_"+xfunc+"_"+disp+"_"+str(atom)+"_x.dat"
                        filenamey = folder+force_in+"_"+xfunc+"_"+disp+"_"+str(atom)+"_y.dat"
                        filenamez = folder+force_in+"_"+xfunc+"_"+disp+"_"+str(atom)+"_z.dat"

                        #print filenamex
                        #print filenamey
                        #print xfunctake.shape,forcesx.shape,forcesy.shape
                        #if xfunctake.shape != forcesx.shape:
                        #    print "xfunctake:"
                        #    print xfunctake
                        #    print "forcesx:"
                        #    print forcesx
                        #    sys.exit("nope")
                        okx = oky = okz = True
                        writeanyhow = False
                        if okx == True:
                            if os.path.isfile(filenamex) != True or writeanyhow == True:
                                #print 'write'
                                np.savetxt(filenamex,np.transpose([xfunctake,forcesx]))
                        if oky == True:
                            if os.path.isfile(filenamey) != True or writeanyhow == True:
                                #print 'write'
                                np.savetxt(filenamey,np.transpose([xfunctake,forcesy]))
                        if okz == True:
                            if os.path.isfile(filenamez) != True or writeanyhow == True:
                                #print 'write'
                                np.savetxt(filenamez,np.transpose([xfunctake,forcesz]))
                        return

                    # 8:    array([ 0.  ,  2.08,  2.08])
                    # 24:   array([ 2.08,  2.08,  0.  ])
                    # 60:   array([-2.08,  2.08,  0.  ])
                    # 126:  array([-2.08, -2.08,  0.  ])

                    #for a in [ 8,24,60,126 ]:
                    #    for b in [ 'xyz', 'quo', 'qun' ]:
                    #        for c in [ 'winkel', 'dnorm', 'longvec' ]:
                    #            #save_analyze('xdir', a, b, c)  # tox
                    #            save_analyze(self.alle.disp, a, b, c)
                    for disp in [ 'xdir', 'quer' ]:
                        for a in [ 8, 60, 126 ]:
                            for b in [ 'xyz', 'xyzorig' ]:
                                for c in [ 'dnorm' ]:
                                    #save_analyze('xdir', a, b, c)  # tox
                                    save_analyze(disp, a, b, c)
        return

    def analyze_qubus(self):
        ''' analysisi of all displacements
            in the pandas pkl file: FML; the pkl is automatically
        '''
        def saveQubuspd(filename):
            ''' saves pandas dataframe '''
            self.pkl = filename
            print(utils.printred("      saving pkl : "+str(filename)))
            with open(filename, 'w') as f:
                #pickle.dump(self.list_save_load, f)
                pickle.dump(self.qubus.data, f)
            return

        def loadQubuspd(filename):
            ''' loads pandas dataframe '''
            if type(filename) != str:
                print("filename:",filename)
                sys.exit("loaddatapkl  needs as input a picklefile (str) but got "+str(filename)+" of type "+str(type(filename)))
            if os.path.isfile(filename) != True:
                sys.exit("file "+filename+" does not exist!")
            #self.data = pd.read_pickle(filename)
            print(utils.printgreen("     loadQubuspd : "+str(filename)))
            with open(filename) as f:
                self.alle.data = pickle.load(f)
                self.update_what_to_calculate_and_write_self_xdir_quer()
            return

        def create_qubus_pkl():
            '''
            entfalte halben (viertel) gerechneten qubus
            DIESE SYMMETRIEN (DIE UNTEN angewendet werden) GELDEN NUR FUER ATOM 24 !!!!!!!!!!
            (und nicht fuer andere Nachbarn) -> KONSEQUENZ: in der MD muessen alle atom zurueck
            gemapped werden in den ersten quadranten!
            ######### erstelle erstmal array und gehe davon aus dass alle dimensionen / werte
            ######### gleich zu z
            '''
            print(utils.printred("hier werden die symmetrien angewandt! und die qubus.pkl erstellt"))
            self.qubus.w = np.array(np.unique(self.qubus.p[:,0][:,2]))
            self.qubus.x = np.unique(np.ndarray.flatten(np.array([self.qubus.w,-self.qubus.w])))
            self.qubus.y = np.unique(np.ndarray.flatten(np.array([self.qubus.w,-self.qubus.w])))
            self.qubus.z = np.unique(np.ndarray.flatten(np.array([self.qubus.w,-self.qubus.w])))

            self.qubus.forces_Fx = False
            self.qubus.forces_Fy = False
            self.qubus.forces_Fz = False
            self.qubus.fml_Fx = False
            self.qubus.fml_Fy = False
            self.qubus.fml_Fz = False
            self.qubus.data = pd.DataFrame(columns = [ \
                    'x', 'y', 'z', \
                    'Fx', 'Fy', 'Fz', \
                    'fmlx', 'fmly', 'fmlz' \
                    ])
            #self.qubus.data = self.qubus.data.astype(object)
            self.qubus.data = self.qubus.data.astype(float)

            for zidx,z in enumerate(self.qubus.z):
                longvec_idx = np.nonzero(self.qubus.p[:,0][:,2] == z)[0]
                print("z:",z) #,longvec_idx


                for idx in longvec_idx:
                    Fx = self.qubus.F[idx,24][0]  # this is the true (full) vasp force
                    Fy = self.qubus.F[idx,24][1]  # this is the true (full) vasp force
                    Fz = self.qubus.F[idx,24][2]  # this is the true (full) vasp force
                    fmlx = self.qubus.fml[idx,24][0]  # this is the true (full) vasp force
                    fmly = self.qubus.fml[idx,24][1]  # this is the true (full) vasp force
                    fmlz = self.qubus.fml[idx,24][2]  # this is the true (full) vasp force
                    #lx = self.qubus.Pr[idx,24][0] - self.qubus.P[idx,24][0]
                    #ly = self.qubus.Pr[idx,24][1] - self.qubus.P[idx,24][1]
                    #lz = self.qubus.Pr[idx,24][2] - self.qubus.P[idx,24][2]
                    #lx = self.qubus.p[idx,0][0]
                    #ly = self.qubus.p[idx,0][1]
                    #lz = self.qubus.p[idx,0][2]
                    lx = round(self.qubus.P[idx,24][0] - self.qubus.Pr[idx,24][0],6)
                    ly = round(self.qubus.P[idx,24][1] - self.qubus.Pr[idx,24][1],6)
                    lz = round(self.qubus.P[idx,24][2] - self.qubus.Pr[idx,24][2],6)
                    # P[idx,24] ist immer [ 1.995,  1.995,  0.   ] und aendert sich nicht
                    # Pr[idx,24] ist jeweils der longvec (aendert sich mit jedem displacement)
                    #print "idx(1456):",idx,self.qubus.p[idx,0],"l:",lx,ly,lz,"F:",Fx,Fy,Fz


                    #########################################################################
                    # Flaeche
                    #########################################################################
                    self.qubus.forces_Fx = utils.append_row_to_2d_array(
                            inarray = self.qubus.forces_Fx, addrow=[ lx, ly, lz, Fx ])
                    self.qubus.forces_Fy = utils.append_row_to_2d_array(
                            inarray = self.qubus.forces_Fy, addrow=[ lx, ly, lz, Fy])
                    self.qubus.forces_Fz = utils.append_row_to_2d_array(
                            inarray = self.qubus.forces_Fz, addrow=[ lx, ly, lz, Fz])

                    self.qubus.fml_Fx = utils.append_row_to_2d_array(
                            inarray = self.qubus.fml_Fx, addrow=[ lx, ly, lz, fmlx ])
                    self.qubus.fml_Fy = utils.append_row_to_2d_array(
                            inarray = self.qubus.fml_Fy, addrow=[ lx, ly, lz, fmly])
                    self.qubus.fml_Fz = utils.append_row_to_2d_array(
                            inarray = self.qubus.fml_Fz, addrow=[ lx, ly, lz, fmlz])

                    self.qubus.data.loc[self.qubus.data.shape[0]] = \
                                [ lx, ly, lz, Fx, Fy, Fz, fmlx, fmly, fmlz ]


                    #########################################################################
                    # add mirroring x/y (for Fx/Fy/Fz force)
                    #########################################################################
                    #if lx != ly:
                    #    self.qubus.forces_Fx = utils.append_row_to_2d_array(
                    #            inarray = self.qubus.forces_Fx, addrow=[ ly, lx, lz, Fy])
                    #    self.qubus.forces_Fy = utils.append_row_to_2d_array(
                    #            inarray = self.qubus.forces_Fy, addrow=[ ly, lx, lz, Fx])
                    #    self.qubus.forces_Fz = utils.append_row_to_2d_array(
                    #            inarray = self.qubus.forces_Fz, addrow=[ ly, lx, lz, Fz])

                    #    self.qubus.fml_Fx = utils.append_row_to_2d_array(
                    #            inarray = self.qubus.fml_Fx, addrow=[ ly, lx, lz, fmly])
                    #    self.qubus.fml_Fy = utils.append_row_to_2d_array(
                    #            inarray = self.qubus.fml_Fy, addrow=[ ly, lx, lz, fmlx])
                    #    self.qubus.fml_Fz = utils.append_row_to_2d_array(
                    #            inarray = self.qubus.fml_Fz, addrow=[ ly, lx, lz, fmlz])

                    #    self.qubus.data.loc[self.qubus.data.shape[0]] = \
                    #            [ ly, lx, lz, Fy, Fx, Fz, fmly, fmlx, fmlz ]

                    #########################################################################
                    # add mirroring in z (for Fx/Fy/Fz force)
                    #########################################################################
                    if lz != 0.0:
                        self.qubus.forces_Fx = utils.append_row_to_2d_array(
                                inarray = self.qubus.forces_Fx, addrow=[ lx, ly, -lz, Fx ])
                        self.qubus.forces_Fy = utils.append_row_to_2d_array(
                                inarray = self.qubus.forces_Fy, addrow=[ lx, ly, -lz, Fy])
                        self.qubus.forces_Fz = utils.append_row_to_2d_array(
                                inarray = self.qubus.forces_Fz, addrow=[ lx, ly, -lz, -Fz])

                        self.qubus.fml_Fx = utils.append_row_to_2d_array(
                                inarray = self.qubus.fml_Fx, addrow=[ lx, ly, -lz, fmlx ])
                        self.qubus.fml_Fy = utils.append_row_to_2d_array(
                                inarray = self.qubus.fml_Fy, addrow=[ lx, ly, -lz, fmly])
                        self.qubus.fml_Fz = utils.append_row_to_2d_array(
                                inarray = self.qubus.fml_Fz, addrow=[ lx, ly, -lz, -fmlz])

                        self.qubus.data.loc[self.qubus.data.shape[0]] = \
                                    [ lx, ly, -lz, Fx, Fy, -Fz, fmlx, fmly, -fmlz ]

                    if lx != ly and lz != 0.0:
                        self.qubus.forces_Fx = utils.append_row_to_2d_array(
                                inarray = self.qubus.forces_Fx, addrow=[ ly, lx, -lz, Fy])
                        self.qubus.forces_Fy = utils.append_row_to_2d_array(
                                inarray = self.qubus.forces_Fy, addrow=[ ly, lx, -lz, Fx])
                        self.qubus.forces_Fz = utils.append_row_to_2d_array(
                                inarray = self.qubus.forces_Fz, addrow=[ ly, lx, -lz, -Fz])

                        self.qubus.fml_Fx = utils.append_row_to_2d_array(
                                inarray = self.qubus.fml_Fx, addrow=[ ly, lx, -lz, fmly])
                        self.qubus.fml_Fy = utils.append_row_to_2d_array(
                                inarray = self.qubus.fml_Fy, addrow=[ ly, lx, -lz, fmlx])
                        self.qubus.fml_Fz = utils.append_row_to_2d_array(
                                inarray = self.qubus.fml_Fz, addrow=[ ly, lx, -lz, -fmlz])

                        self.qubus.data.loc[self.qubus.data.shape[0]] = \
                                [ ly, lx, -lz, Fy, Fx, -Fz, fmly, fmlx, -fmlz ]

            saveQubuspd(folderqubus+"/qubus.pkl")
            return




        def create_qubus_npz(folder):
            from . import qubus_interpolate
            q = qubus_interpolate.qubus_splines()
            q.loadQubuspd(folder+'/qubus.pkl')
            q.ef_in_whole_qubus(folder) # creates qubus.npz
            return

        def load_qubus_npz(folder):
            print(utils.printred("loading "+folder+'/qubus.fml.npz'))
            var = np.load(folder+'/qubus.fml.npz')
            self.fmlx = var['fmlx']
            self.fmly = var['fmly']
            self.fmlz = var['fmlz']
            self.emlx = var['emlx']
            self.emly = var['emly']
            self.emlz = var['emlz']
            self.x = var['x']
            self.y = var['y']
            self.z = var['z']

            print(utils.printred("loading "+folder+'/qubsu.F.npz'))
            var = np.load(folder+'/qubus.F.npz')
            self.Fx = var['Fx']
            self.Fy = var['Fy']
            self.Fz = var['Fz']
            self.Ex = var['Ex']
            self.Ey = var['Ey']
            self.Ez = var['Ez']
            self.x = var['x']
            self.y = var['y']
            self.z = var['z']
            return




        ###################################################################################
        # load in data or whole qubus
        ###################################################################################
        for folder in self.folder_displacement_direction_all:
            if 'qubus' in folder:
                folderqubus = folder
                if os.path.isfile(folder+'/qubus.pkl'):
                    self.alle.disp = 'qubus'
                    loadQubuspd(folderqubus+"/qubus.pkl")
                else:
                    self.alle.loadnnforcesnpz(folder,"nnforces")
                    self.alle.loadFMLnpz(folder,"FML")
                    self.update_what_to_calculate_and_write_self_xdir_quer()
                    self.alle.loaddatapkl(self.pkl)

                if os.path.isfile(folder+'/qubus.F.npz') != True or \
                   os.path.isfile(folder+'/qubus.fml.npz') != True or \
                   self.create_qubus_npz == True:
                    create_qubus_npz(folder)
                else:
                    load_qubus_npz(folder)


        ###################################################################################
        # if pkl does not exist then create it
        ###################################################################################
        if os.path.isfile(folder+'/qubus.pkl') != True:
            print("folder"+'/qubus.pkl does not exist ... --> creating qubus')
            create_qubus_pkl()

        def plotdata(z = False, plot_to = False):
            ''' F: {Fx,Fy,Fz,fmlx,fmly,fmlz}
                z: {-1.3,-1.2, ..., 0.0, ..., 1.2, 1.3} '''

            plot_to = 1.0

            ########### xy
            self.qubus.xy   = self.qubus.data[(self.qubus.data.z == z) & \
                    (self.qubus.data.x <=  plot_to) & \
                    (self.qubus.data.x >= -plot_to) & \
                    (self.qubus.data.y <=  plot_to) & \
                    (self.qubus.data.y >= -plot_to) \
                    ][['x','y']].values

            ########### z
            values_forces_x = self.qubus.data[(self.qubus.data.z == z) & \
                    (self.qubus.data.x <=  plot_to) & \
                    (self.qubus.data.x >= -plot_to) & \
                    (self.qubus.data.y <=  plot_to) & \
                    (self.qubus.data.y >= -plot_to) \
                    ]['Fx'].values
            values_forces_y = self.qubus.data[(self.qubus.data.z == z) & \
                    (self.qubus.data.x <=  plot_to) & \
                    (self.qubus.data.x >= -plot_to) & \
                    (self.qubus.data.y <=  plot_to) & \
                    (self.qubus.data.y >= -plot_to) \
                    ]['Fy'].values
            values_forces_z = self.qubus.data[(self.qubus.data.z == z) & \
                    (self.qubus.data.x <=  plot_to) & \
                    (self.qubus.data.x >= -plot_to) & \
                    (self.qubus.data.y <=  plot_to) & \
                    (self.qubus.data.y >= -plot_to) \
                    ]['Fz'].values

            values_fml_x = self.qubus.data[(self.qubus.data.z == z) & \
                    (self.qubus.data.x <=  plot_to) & \
                    (self.qubus.data.x >= -plot_to) & \
                    (self.qubus.data.y <=  plot_to) & \
                    (self.qubus.data.y >= -plot_to) \
                    ]['fmlx'].values
            values_fml_y = self.qubus.data[(self.qubus.data.z == z) & \
                    (self.qubus.data.x <=  plot_to) & \
                    (self.qubus.data.x >= -plot_to) & \
                    (self.qubus.data.y <=  plot_to) & \
                    (self.qubus.data.y >= -plot_to) \
                    ]['fmly'].values
            values_fml_z = self.qubus.data[(self.qubus.data.z == z) & \
                    (self.qubus.data.x <=  plot_to) & \
                    (self.qubus.data.x >= -plot_to) & \
                    (self.qubus.data.y <=  plot_to) & \
                    (self.qubus.data.y >= -plot_to) \
                    ]['fmlz'].values



            # TODO: - make slice for every z
            #       - save all in npz file
            #       - only load npz file and play with that

            grid_x, grid_y = np.mgrid[-1.3:1.3:27j,-1.3:1.3:27j]
            grid_x, grid_y = np.mgrid[-1.0:1.0:21j,-1.0:1.0:21j]

            from scipy.interpolate import griddata
            forces_x2 = griddata(self.qubus.xy, values_forces_x, (grid_x, grid_y), method='cubic')
            forces_y2 = griddata(self.qubus.xy, values_forces_y, (grid_x, grid_y), method='cubic')
            forces_z2 = griddata(self.qubus.xy, values_forces_z, (grid_x, grid_y), method='cubic')

            fml_x = griddata(self.qubus.xy, values_fml_x, (grid_x, grid_y), method='cubic')
            fml_y = griddata(self.qubus.xy, values_fml_y, (grid_x, grid_y), method='cubic')
            fml_z = griddata(self.qubus.xy, values_fml_z, (grid_x, grid_y), method='cubic')
            return forces_x2,forces_y2,forces_z2,fml_x,fml_y,fml_z


        forces_x2,forces_y2,forces_z2,fml_x,fml_y,fml_z = plotdata(z=0.0)

        def plot_fml(xyz=False,z=False):
            forces_x2,forces_y2,forces_z2,fml_x,fml_y,fml_z = plotdata(z=z)
            if xyz == 'x':return fml_x
            if xyz == 'y':return fml_y
            if xyz == 'z':return fml_z
        def plot_forces(xyz=False,z=False):
            forces_x,forces_y,forces_z,fml_x,fml_y,fml_z = plotdata(z=z)
            if xyz == 'x':return forces_x
            if xyz == 'y':return forces_y
            if xyz == 'z':return forces_z


        def plot_row(   where = False,
                        what = False,
                        title = False):
            '''     where=[331,332,333],
                    what=[plot_fml('x',0.0),plot_fml('x',0.0),plot_fml('x',0.0)]
            '''


            for idx,i in enumerate(where):
                #plt.figure()
                plt.subplot(i)
                #print "dim:",what[idx].T.shape
                # xi,yi,zi = what[idx].T # too many values to unpack
                #plt.imshow(what[idx].T, extent=(-1.3,1.3,-1.3,1.3), cmap=plt.cm.RdBu, origin='lower')
                # fuer fml
                fmax = 2.4
                gr = 1.3
                gr = 1.0
                contdist = 0.2

                # fuer forces
                fmax = 70.0
                gr = 1.0
                contdist = 1.5

                #print what[idx].shape,what[idx].max(),what[idx].min(),"clindist:",contdist
                print(what[idx].max(),what[idx].min())
                x = y = np.linspace(-gr,gr,gr*10*2+1)

                cp1 = plt.contourf(x,y,what[idx].T,51,cmap=plt.cm.RdBu_r,\
                        vmin=-fmax, vmax=fmax)
                #cp1 = plt.contourf(x,y,what[idx].T,51,\
                #        vmin=-fmax, vmax=fmax)
                mylvls = np.arange( -fmax, fmax, contdist )
                cp2 = plt.contour(x,y,what[idx].T,levels=mylvls,colors='k')

                # hinzufuegen 0- lines
                fp1 = plt.plot( [-gr,gr], np.zeros(2), '-r' )
                fp1 = plt.plot( np.zeros(2), [-gr,gr], '-r' )
                ka = np.loadtxt(self.folder_base+"/ti/Ir/31__PTS_dosall_from_2x2x2sc__LON_quer_3x3x3kp_all_morse_LONADD_YES_LON2_NONE_LON2ADD__NO_TOX_NONE_TOY_NONE_TOZ_NONE/lambda1.0/tests/longecxyz/DOS2d_0.0_0.1__0.005_0.006")
                fp1 = plt.plot( ka[:,0], ka[:,1], 'or' )
                #plt.contour(what[idx].T,15, cmap=plt.cm.rainbow)

                #plt.contourf(xi, yi, zi, 15, cmap=plt.cm.rainbow,
                #    norm=plt.Normalize(vmax=abs(zi).max(), vmin=-abs(zi).max()))
                #plt.clim([-fmax, fmax])
                plt.colorbar(cp1)
                #print "-->",str(what[idx])
                #print "|| ",inspect.stack()
                #for line in inspect.stack():
                #    print "--->",line
                #    #sources.write(re.sub(stringsearch, stringreplace, line))
                #print "||0",inspect.stack()[0]
                #print "||1",inspect.stack()[1]
                print("")
                if title == False:
                    plt.title("")
                else:
                    plt.title(title[idx])
                plt.clabel(cp2, inline=1, fontsize=8)
                plt.axis('tight')
                #plt.axis('equal')

        plt.clf()
        #plot_row(where  =[331,332,333],
        #         what   =[plot_forces('x',0.0),plot_forces('y',0.0),plot_forces('z',0.0)])

        if False:
            z = 0.0
            plot_row(where  =[331,332,333],
                     what   =[plot_fml('x',z),plot_fml('y',z),plot_fml('z',z)],
                     title  =['fmlx z='+str(z),'fmly z='+str(z), 'fmlz z='+str(z)])
            z = 0.2
            plot_row(where  =[334,335,336],
                     what   =[plot_fml('x',z),plot_fml('y',z),plot_fml('z',z)],
                     title  =['fmlx z='+str(z),'fmly z='+str(z), 'fmlz z='+str(z)])

            z = 0.4
            plot_row(where  =[337,338,339],
                     what   =[plot_fml('x',z),plot_fml('y',z),plot_fml('z',z)],
                     title  =['fmlx z='+str(z),'fmly z='+str(z), 'fmlz z='+str(z)])
        if True:
            z = 0.0
            plot_row(where  =[331,332,333],
                     what   =[plot_forces('x',z),plot_forces('y',z),plot_forces('z',z)],
                     title  =['fmlx z='+str(z),'fmly z='+str(z), 'fmlz z='+str(z)])
            z = 0.2
            plot_row(where  =[334,335,336],
                     what   =[plot_forces('x',z),plot_forces('y',z),plot_forces('z',z)],
                     title  =['fmlx z='+str(z),'fmly z='+str(z), 'fmlz z='+str(z)])

            z = 0.4
            plot_row(where  =[337,338,339],
                     what   =[plot_forces('x',z),plot_forces('y',z),plot_forces('z',z)],
                     title  =['fmlx z='+str(z),'fmly z='+str(z), 'fmlz z='+str(z)])
        return

    def get_ef_rest(self,x,y,z):
        scale = (self.x.size-1)/2./self.x.max()
        #print "self.x:",self.x
        #print "self.x.max()",self.x.max(),self.x.size,self.x.size,(self.x.size-1)/2,scale
        def xyz_to_map_coords(x):
            ''' x = -1.3 ---> out = 0
                x = 0    ---> out = 13
                x = 1.3  ---> out = 26 '''
            return (x+self.x.max())*scale

        def map_coords_to_xyz(map_coords):
            return map_coords*0.1-1.3

        coords = np.array([[xyz_to_map_coords(x), xyz_to_map_coords(y),xyz_to_map_coords(z)]])
        coords = coords.T
        #zi = scipy.ndimage.map_coordinates(q.fall, coords, order=2, mode='nearest')
        fx = map_coordinates(self.fmlx, coords, order=2, mode='nearest')[0]
        fy = map_coordinates(self.fmly, coords, order=2, mode='nearest')[0]
        fz = map_coordinates(self.fmlz, coords, order=2, mode='nearest')[0]
        ex = map_coordinates(self.emlx, coords, order=2, mode='nearest')[0]
        ey = map_coordinates(self.emly, coords, order=2, mode='nearest')[0]
        ez = map_coordinates(self.emlz, coords, order=2, mode='nearest')[0]
        return [fx,fy,fz],[ex,ey,ez]

