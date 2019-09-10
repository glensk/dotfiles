#!/usr/bin/env python

import numpy as np
import sys
import os
import glob
#import pylab
import argparse   # vasp 2.6.6 would be too old
import filecmp
import myutils

np.set_printoptions(suppress=True)   # display arrays withou 000001
np.set_printoptions(precision=6)    # print only 6 digist after .

# functions with _name are not seen when imported
def _printred(var):
    ENDC = '\033[0m'
    red = '\033[31m'
    print(red + str(var) + ENDC)

def _printgreen(var):
    ENDC = '\033[0m'
    red = '\033[32m'
    print(red + str(var) + ENDC)

def is_int(x):
    try:
        a = float(x)
        b = int(a)
    except ValueError:
        return False
    else:
        return a == b

def vinet(V, e0=None, v0=None, b0=None, b0der=None, convert_to_sympy = False):
    """
    vinet Equation of state (EOS) as a function of:
    V [Angstrom^3/atom], e0[meV/atom], v0[Angsrom^3], b0[GPa], b0der[arb.]
    output: energy [meV/atom]

    vinet equation from PRB 70, 224107

    convert_to_sympy makes this function use with the sympy package
    """
    if convert_to_sympy:
        from sympy import exp
    else:
        from numpy import exp
    #print "type exp:",type(exp)
    meVByAngstromToGPa = .1602176462
    b0 = b0 / meVByAngstromToGPa
    return e0 + (4. * b0 * v0) / (b0der - 1.) ** 2. \
        - 2. * v0 * b0 * (b0der - 1.) ** (-2.) \
        * (5. + 3. * b0der * ((V / v0) ** (1. / 3.) - 1.) - 3 * (V / v0) ** (1. / 3.)) \
        * exp(-(3. / 2.) * (b0der - 1.) * ((V / v0) ** (1. / 3.) - 1.))

def vinet_diff_1(V, e0, v0, b0, b0der, convert_to_sympy = False):
    ''' first derivative of vinet eos'''
    if convert_to_sympy:
        from sympy import exp
    else:
        from numpy import exp
    meVByAngstromToGPa = .1602176462
    b0 = b0 / meVByAngstromToGPa
    return -12.4830194890231*b0*v0*(b0der - 1.0)**(-2.0)*(1.0*b0der*(V/v0)**(1./3.)/V - 1.0*(V/v0)**(1./3.)/V)*exp((-1.5*b0der + 1.5)*((V/v0)**(1./3.) - 1.0)) - 4.16100649634102*b0*v0*(V/v0)**(1./3.)*(-1.5*b0der + 1.5)*(b0der - 1.0)**(-2.0)*(3.0*b0der*((V/v0)**(1./3.) - 1.0) - 3*(V/v0)**(1./3.) + 5.0)*exp((-1.5*b0der + 1.5)*((V/v0)**(1./3.) - 1.0))/V

def lennardjones(r, eps=None, rm=None, convert_to_sympy = False):
    """
    Lennard-Jones Equation of state (EOS) as a function of:
    r [Angstrom] (distance between particles),
    rm [Anstrom] (minimum position) (Nearest neighbor distance),
    eps[eV] (depth of potential well),
    output: energy [meV/atom]

    vinet equation from PRB 70, 224107

    convert_to_sympy makes this function use with the sympy package
    """
    if convert_to_sympy:
        from sympy import exp
    else:
        from numpy import exp
    # *1000 since we now expect everything in meV
    return eps*( (rm/r)**12. - 2.*(rm/r)**6. )*1000

def lennardjones_diff_1():
    pass

def murn(V, e0=None, v0=None, b0=None, b0der=None):
    """
    Murnaghan Equation of state (EOS) as a function of:
    V [Angstrom^3], e0[eV], v0[Angsrom^3], b0[GPa], b0der[arb.]
    output: energy [meV]

    From PRB 28,5480 (1983)
    """
    meVByAngstromToGPa = .1602176462
    b0 = b0 / meVByAngstromToGPa
    return e0 + (b0 * V)/(b0der * (b0der - 1.))*(b0der*(1. - v0/V) + (v0/V)**b0der - 1.)

def murn_diff_1():
    pass

def birch(V, e0=None, v0=None, b0=None, b0der=None):
    """ Birch Equation of state (EOS) as a function of:
    V [Angstrom^3], e0[eV], v0[Angsrom^3], b0[GPa], b0der[arb.]
    output: energy [meV] """
    meVByAngstromToGPa = .1602176462
    b0 = b0 / meVByAngstromToGPa
    return e0 + (9.*v0*b0)/16.*(((v0/V)**(2./3.)-1.)**3.*b0der + ((v0/V)**(2./3.) - 1.)**2.*(6. - 4. *(v0/V)**(2./3.)))

def run2(command=None, dont_raise_exceptino = False):
    """
    constantly prints output, not just at the end
    """
    import subprocess
    process = subprocess.Popen(
        command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    output = ''

    # Poll process for new output until finished
    for line in iter(process.stdout.readline, ""):
        print(line) #, end=' ')
        output += line
    process.wait()
    exitCode = process.returncode

    if dont_raise_exceptino == True:
        return output

    if (exitCode == 0):
        return output
    else:
        raise Exception(command, exitCode, output)

def run(befehl):
    """
    z.B. run('ls -l')
    z.B. run('frompath_stoich.sh')
    """
    import subprocess
    return subprocess.check_output(befehl, shell=True, stderr=subprocess.STDOUT)

def create_energy_dat_from_underlying_OUTCARS():
    _printgreen("creating energy.dat")
    filename = "energy.dat"
    if os.path.isfile(filename):
        filename = "energy.new.dat"
    run2("OUTCAR_A.sh OUTCAR_volume-lastexact.sh OUTCAR_ene-sigma0-last.sh OUTCAR_number_of_atoms.sh list > "+str(filename), dont_raise_exceptino = True)
    # old run2("OUTCAR_A.sh OUTCAR_volume-last.sh OUTCAR_ene-sigma0-last.sh OUTCAR_number_of_atoms.sh list > "+str(filename), dont_raise_exceptino = True)
    if os.path.isfile("energy.new.dat") == True:
        if os.path.isfile("energy.dat") == True:
            if filecmp.cmp('energy.new.dat', 'energy.dat') == True:
                filename = "energy.dat"
                os.remove("energy.new.dat")
    return filename


class eos(object):
    def __init__(self, args = None):
        # __init__ make all importat definitions and sys.exit() if sth necessary does not exist!!
        # __init__ options can be changed
        # a) after initializing eos()
        # b) Interactively using  help
        self.equation = "vinet"
        self.inputfile = "energy.dat"
        self.inputfile_parameters = None
        #specific to a certaion eos {vinet, birch, murn, ...}
        #self.equation = args.equation #"vinet"
        self.v0 = None
        self.e0 = None
        self.b0 = None
        self.b0der = None
        self.parameters = None # [e0, v0, b0, b0der]
        self.fitvolumes = None
        self.fitenergies = None
        self.fitdeltas = None
        self.data = False

        #specifit for importing
        self._verbose = False
        self._scale_output  = 1
        self._input_volume_divide   = 1.
        self._input_energy_divide   = 1.
        self._alat_fcc = False
        self._alat_bcc = False
        self._atoms = 1
        self._considerpoints = False
        self._units_energy = 'eV'       # per atom
        self._units_volume = 'Ang^3'    # per atom
        self._outfile = "E"+self.equation[:1].upper()+self.equation[1:].lower()+"_"+str(self._scale_output)
        self._write_energydat = False
        self._write_energyperatomdat = False
        self._inputfiles = None
        self._c = None # create energy.dat

        if args:
            self.equation = args.equation
            self._verbose = args.verbose
            self._scale_output  = args.atoms
            self._input_volume_divide   = args.input_volume_divide
            self._input_energy_divide   = args.input_energy_divide
            self._alat_fcc = args.fcc
            self._alat_bcc = args.bcc
            self._considerpoints = args.considerpoints
            self._units_energy = args.ue
            self._units_volume = args.uv
            self._outfile = "E"+self.equation[:1].upper()+self.equation[1:].lower()+"_"+str(self._scale_output)
            self._write_energydat = args.w
            self._write_energyperatomdat = args.wa
            self.inputfile = args.inputfile
            self._inputfiles = [ os.path.abspath(folder + '/' + args.inputfile) for folder in args.folder]
            self._c = args.c

        #help (may generally only change __init__ optons when calling from system shell)
        self.help()

    def help(self, p = None):
        string = '''
        Fit of energies vs volume curve to equation-of-state (EOS).
        EOS types implemented for fitting: vinet Murn Birch.
        The energy vs volume curve is read from a inputfile (default: energy.dat).
        Expected format of inputfile:
            - volume: (1st column) [Ang^3/per atom]
            - energy: (2nd column) [eV/per atom]
        OUTPUT:
            - Evinet_1 containing E0, V0, B0, B0der in [meV/atom, Ang^3/atom, GPa, no units]
        '''
        if p == None:
            from argparse import ArgumentDefaultsHelpFormatter
            #from argparse import RawDescriptionHelpFormatter
            #from argparse import RawTextHelpFormatter
            #p = argparse.ArgumentParser(description=string, formatter_class=RawTextHelpFormatter)
            p = argparse.ArgumentParser(description=string,
                    formatter_class=ArgumentDefaultsHelpFormatter)
                    #formatter_class=RawTextHelpFormatter) # No
                    #formatter_class=RawDescriptionHelpFormatter)  # No
        p.add_argument('-c',
            help='create energy.dat file from underlying OUTCARS',
            action='store_true', default=False)
        p.add_argument('-e',   '--equation', choices=['vinet', 'birch', 'murn', 'lennardjones'],
           help='specify equation of state', default=self.equation)
        p.add_argument('-i',   '--inputfile',
           help='specify name of inputfile',
           type=str, default=self.inputfile)
        p.add_argument('-ue', choices=['eV', 'meV', 'Hartree'],
            help='specify units of inputfile energy', default=self._units_energy)
        p.add_argument('-uv', choices=['Ang^3', 'Bohr^3'],
            help='specify units of inputfile volume', default=self._units_volume)
        p.add_argument('-ivd', '--input_volume_divide',
            help='scales the 1st column (=volume) of the inputfile by 1/float',
            type=float, default=self._input_volume_divide)
        p.add_argument('-ied', '--input_energy_divide',
            help='scales the 2nd column (=energy) of inputfile by 1/float',
            type=float, default=self._input_energy_divide)
        p.add_argument('-fcc',
            help='1st column (=volume) of inputfile (energy.dat) is fcc alat e.g. 4.05-4.16 for Aluminum',
            action='store_true', default=self._alat_fcc)
        p.add_argument('-bcc',
            help='1st column (=volume) of inputfile (energy.dat) is bcc alat',
            action='store_true', default=self._alat_bcc)
        p.add_argument('-n',   '--considerpoints',
            help='take only the n lowest energies (and correstpoinding volumes) for fitting',
        type=int, default=self._considerpoints)

        #p.add_argument('-o',   '--outputfile',
        #    help='change name of outputfile', type=str, default=self._outfile)
        p.add_argument('-w',
            help='write energy.dat file', action='store_true', default=False)
        p.add_argument('-wa',
            help='write energyperatom.dat file', action='store_true', default=False)
        p.add_argument('-a',  '--atoms',
        help='create Evinet_ATOMS  and not EVinet_1. This scales e0 and v0 to number of ATOMS',
        type=int, default=1)
        p.add_argument('-f','--folder',nargs='+',
                help='Do evaluation in every given folder; You can use wildcards (*.{}...)', default=[os.getcwd()])
        p.add_argument('-v', '--verbose',
                help='verbose', action='store_true', default=self._verbose)
        return p

    def import_energy_vs_volume_data(self,data=False):
        """ Imports the first tow columns from a file containing:
        1st column: Volumes per atom [Angstrom^3]
        2nd column: Energies per atom [eV] """
        if type(data) == bool:
            if os.path.isfile(self.inputfile) is not True:
                sys.exit("Necessary inputfile \"" + str(self.inputfile) + "\" not found.")
            else:
                self.inputfile = os.path.abspath(self.inputfile)

            #_printgreen("importing "+str(os.path.abspath(self.inputfile))+" ...")
            _printgreen("importing "+str(os.path.relpath(self.inputfile))+" ...")

            self.folder = os.path.split(os.path.realpath(self.inputfile))[0]

            #data = pylab.loadtxt(self.inputfile)
            data = np.loadtxt(self.inputfile)
        else:
            data = data
        self.data = data[data[:, 0].argsort()]
        self.data_in = self.data
        #print('skkkkkkk',self.data)
        #sys.exit()

        # if we have a third column in energy.dat file
        if self.data.shape[1] == 3:
            # if all the numbers in 3rd column are unique, those are the number of atoms
            if len(np.unique(self.data[:,2])) == 1:
                check_atoms = np.unique(self.data[:,2])[0]
                if is_int(check_atoms) != True:
                    _printred("3rd column of "+str(os.path.relpath(self.inputfile))+" has float values")
                    sys.exit()
                if self._input_energy_divide == 1.:
                    _printred("3rd column of "+str(os.path.relpath(self.inputfile))+": "+str(check_atoms)+" atoms.")
                    self._atoms = int(check_atoms)
                    pass

                    #check_atoms = np.unique(self.data[:,2])[0]
                    #if is_int(check_atoms):
                    #    self._input_energy_divide = int(np.unique(self.data[:,2])[0])
                    #    self._input_volume_divide = int(np.unique(self.data[:,2])[0])
                    #    _printgreen(">>> scaling energy and volume to "+\
                    #            str(self._input_volume_divide))
            else:
                _printred("2rd row of "+str(self.inputfile)+" has diffrent number of atoms")
                sys.exit()

        self.data = np.copy(self.data[:,0:2])

        # scale data
        dataxall = self.data[:, 0]/self._input_volume_divide
        if self._alat_fcc is True:
            dataxall = self.data[:, 0]**3/4
        elif self._alat_bcc is True:
            dataxall = self.data[:, 0]**3/2

        if self._units_volume == "Ang^3":
            pass
        elif self._units_volume == "Bohr^3":
            dataxall = dataxall * 6.7483346


        datayall = self.data[:, 1]*1000/self._input_energy_divide  # *1000 since murn and so on need energy in meV
        if self._units_energy == 'eV':
            pass
        elif self._units_energy == 'meV':
            datayall = datayall / 1000.
        elif self._units_energy == 'Hartree':
            datayall = datayall * 27.211384

        datax = dataxall
        datay = datayall

        # just for the case if we would only like to fit to 10 lowest y values
        if self._considerpoints:
            mask = datayall.argsort()[:self._considerpoints]
            datax = dataxall[mask]
            datay = datayall[mask]

        self.data[:,0] = datax
        self.data[:,1] = datay
        self.datax = datax
        self.datay = datay

        return self.data

    def import_parameters_data(self, filename = None):
        possible_files = [ "EVinet*" , "Evinet*", "EMurn*", "Emurn*", "Ebirch*", "EBirch*" ]

        if filename == None:
            for i in possible_files:
                filecheck = glob.glob(i)
                if filecheck:
                    if len(filecheck) == 1:  ## if only one file found
                        filename = filecheck[0]

        if filename == None:
            sys.exit("EOS file (e.g. EVinet) not found!")
        if os.path.isfile(filename) is not True:
            sys.exit("Necessary inputfile \"" + str(filename) + "\" not found; exit")

        #_printgreen("importing "+filename+" ..............................................")
        self.inputfile_parameters = filename
        #parameters = pylab.loadtxt(filename)
        parameters = np.loadtxt(filename)
        self._outfile = filename
        self.e0 = parameters[0]
        self.v0 = parameters[1]
        self.b0 = parameters[2]
        self.b0der = parameters[3]
        self.parameters = [self.e0, self.v0, self.b0, self.b0der]
        #_printgreen("importing "+filename+" .............................................. DONE")
        #print ""
        return self.parameters

    def import_surface_expr(self, filename = None, temperatures = None, convert_to_sympy = True):
        ''' '''
        if temperatures == None:
            sys.exit("pleas provide temperatures as input")
        possible_filename = [ "EVinet*", "EMurn*", "EBirch*" ]
        if filename == None:
            for i in possible_filename:
                filecheck = glob.glob(i)
                if filecheck:
                    if len(filecheck) == 1:
                        self.parameters = self.import_parameters_data(filecheck[0])
                        if "inet" in filecheck[0]:
                            self.equation = "vinet"
                        elif "urn" in filecheck[0]:
                            self.equation = "murn"
                        elif "irch" in filecheck[0]:
                            self.equation = "birch"
                        if self.equation == None:
                            sys.exit("which eos is this?")
                        self.parameters = np.loadtxt(filecheck[0])
        if self.parameters == None:
            print("Did not find EOS file (EViet*, EMurn*, EBirch*)")
            return None
        import sympy
        V = sympy.Symbol('V')
        add = vinet(V, *self.parameters, convert_to_sympy=convert_to_sympy)
        self.surface_expr = np.array([])
        for i in temperatures:
            self.surface_expr = np.append(self.surface_expr, add)
        return self.surface_expr


    def fit_to_energy_vs_volume_data(self, datax=False, datay=False, points=100):
        """
        defines:
        self.v0         : equilibrium/minimum volume per atom [Angstrom^3]
        self.e0         : energy at v0 in [meV]
        self.b0         : Bulk modulus in [GPa]
        self.b0der      : Bulk modulus derivative [arb. units]
        self.parameters: [e0[meV], v0[Angsrom^3], b0[GPa], b0der[arb.]]

        input:
            volumes_per_atom : [ v1, v2, ... ] in Angstrom^3
            energies_per_atom: [ e1, e2, ... ] in eV

        optional input:
            points : defines the number of volume points for the fit
        """
        #_printgreen("fitting energy vs volume data to: >>> "+str(self.equation))
        eos = eval(self.equation) #self.vinet

        v0smallererror = 3
        v0biggererror = 50
        v0smallerwarnin = 8
        v0biggerwarnin = 35



        if type(datax) == bool and type(datay) == bool:
            if type(self.data) == bool:
                sys.exit('no data given!!!')
            else:
                datax = self.data[:,0]
                datay = self.data[:,1]

        if type(datax) != bool and type(datay) != bool:
            data = np.transpose([datax,datay])
            self.import_energy_vs_volume_data(data=data)
            datax = self.datax
            datay = self.datay

        self.datax = datax
        self.datay = datay
        if self._verbose:
            print("datax (in Angstrom^3/atom):",datax)
            print("datay (in meV/atom)       :",datay)

        if len(datax) is not len(datay):
            sys.exit("volumes_per_atom have different lenght than energies_per_atom")

        # guess for output parameters
        b0guess = 111.
        b0derguess = 4.
        enemin = np.min(datay)
        ind = np.nonzero(datay == np.min(datay))[0]
        volume_enemin = datax[ind][0]
        guess = [enemin, volume_enemin, b0guess, b0derguess]
        if self.equation == "lennardjones":
            #guess = [ eps, rm ]
            guess = [ 16. , volume_enemin ]
        if self._verbose:
            print("guess:",guess)

        # fit to data
        from scipy import optimize
        parameters, evinet_covariance = \
            optimize.curve_fit(eos, datax, datay, guess, maxfev=1000)
        if self._verbose:
            _printgreen(str(eos))
            _printgreen(str(parameters))

        # define self.parameters
        self.v0 = parameters[1]
        if v0smallerwarnin*self._atoms < self.v0 < v0biggerwarnin*self._atoms:
            pass # everything outside of warning range
        else:
            _printred("Your minimum volume ( v0 = " + str(self.v0) + ") seems not to be scaled per atom." )
            _printred("It is ok if you dont't want to scale per atom. This is just a note for you to make a sanity check.")
            #_printred('seems souspiciously small/large. Did you really scale per atom?')
            _printred('Warning range: '+str(v0smallerwarnin)+' > v0 > '+str(v0biggerwarnin))
            if v0smallererror < self.v0 < v0biggererror:
                pass # everything outsind of error range
            else:
                _printred('Exit range   : '+str(v0smallererror)+' > v0 > '+str(v0biggererror))
                #sys.exit("I will exit here")
        if self.equation == "lennardjones":
            self.parameters = [ parameters[0],parameters[1] ]
        else:
            self.e0 = parameters[0] #* 1000
            self.b0 = parameters[2] #* meVByAngstromToGPa * 1000
            self.b0der = parameters[3]
            self.parameters = [self.e0, self.v0, self.b0, self.b0der]
        if self._verbose:
            _printgreen(">>> "+str(self.equation)+" <<< ",)
            print(("pars: "+str(self.parameters)))

        # fit (it would be more nice to solve maximum energy and plot to corr. volume)
        dx = max([abs(datax.min() - self.v0), abs(datax.max() - self.v0)])*1.03
        if self._verbose:
            print("datax",datax)
            print("datax.min()",datax.min())
            print("self.v0",self.v0)
            print("datax.max",datax.max())
            print("dx",dx)
        self.fitvolumes = np.linspace(self.v0 - dx, self.v0 + dx, num=points)
        if self._verbose:
            print(("fm:",self.fitvolumes,"aaa",parameters))
        self.fitenergies = eos(self.fitvolumes, *parameters)


        # in case energies were scaled, shift fitenergies accordingly
        self.fitvolumes_in  = self.fitvolumes*self._input_volume_divide
        if self._alat_fcc is True:
            self.fitvolumes_in = (self.fitvolumes*4.)**(1./3.)
        if self._alat_bcc is True:
            self.fitvolumes_in = (self.fitvolumes*2.)**(1./3.)
        if self._units_volume == "Bohr^3":
            self.fitvolumes_in = self.fitvolumes / 6.7483346


        self.fitenergies_in = self.fitenergies*self._input_energy_divide/1000  # /1000 since we expected everything in eV
        if self._units_energy == 'meV':
            self.fitenergies_in = self.fitenergies * 1000.
        if self._units_energy == 'Hartree':
            self.fitenergies_in = self.fitenergies / 27.211384


        # delta to fit
        self.fitdeltas = eos(datax, *parameters) - datay
        if self._verbose:
            print("datay:",datay)
            print("datan:", eos(datax, *parameters))
            print("fitde:",self.fitdeltas)

        # check if delta to fit has no hop (e.g. due to change in NGX(F))
        #d = np.diff(self.fitdeltas)/np.diff(self.datax)
        self._fitdeltas_grad = np.gradient(self.fitdeltas)/np.gradient(self.datax)

        # print a WARNING in case of a jump with 2 standard deviations
        std2 = self._fitdeltas_grad.std()*3.0
        self._fitdeltas_grad_greater = [i for i in abs(self._fitdeltas_grad) if i >= std2]
        if len(self._fitdeltas_grad_greater) != 0:
            m1="WARNING: The difference: (fitdata - datain) seem to have a jump; "
            m0="std: "+str(std2)
            m2="Have a look at eos.fitdeltas; "
            m3="A possible cause is if NGX is changed between volumes. "
            m4="VASP might do this automatically!"
            print("abs gradient fitdeltas:",abs(self._fitdeltas_grad))
            print("abs greater:",self._fitdeltas_grad_greater)
            _printred(m1+m2+m3+m4+m0)
            #np.savetxt('datain.dat',np.transpose([datax,datay]))
            #np.savetxt('datafit.dat',np.transpose([datax,eos(datax, *parameters)]))


        #return [self.e0, self.v0, self.b0, self.b0der]
        return self.fitvolumes_in, self.fitenergies_in

    #def get_surface(self, volumes, temperatures, example_surface):
    #    eos = eval(self.equation) #self.vinet
    #    energies0 = eos(volumes, *self.parameters)*1000
    #    e0 = np.copy(example_surface)
    #    for ind, temp in enumerate(temperatures):
    #        e0[:,ind] = energies0
    #    self.surface = e0
    #    return self.surface


    def write_data(self,file1="energy.dat.fit",file2="energy.fitdeltas"):
        """ Save files:
                - energy.dat.fit
                - Evinet_ATOMS
                - energy.fitdelta """
        folder_ = os.path.split(os.path.relpath(self.inputfile))[0]
        folder = os.path.split(os.path.realpath(self.inputfile))[0]

        # energy.dat.fit
        print("writing: "+file1)
        np.savetxt(
                folder + "/" + file1,
                np.transpose([self.fitvolumes_in, self.fitenergies_in]),
                fmt='%.18f',
                delimiter='  ')   # X is an array

        # Evinet.parameters get number of atoms
        add = "_"+str(self._atoms)
        if self._scale_output != 1:
            add = "_"+str(self._scale_output)
        if self._verbose:
            print("self._atoms:",self._atoms)
            print("self._scale_output:",self._scale_output)
            print("add:",add)

        # get filename if different fome (energy)_b_32(.dat)
        print("os.path.basename(self.inputfile)[:6]",os.path.basename(self.inputfile)[:6])
        add_check = False
        if os.path.basename(self.inputfile)[:6] == "energy":
            if os.path.basename(self.inputfile)[-4:] == ".dat":
                add_check = os.path.basename(self.inputfile)[6:-4]
                if add_check == '':
                    pass
                else:
                    add = add+add_check
        if self._verbose:
            print("add_check:",add_check)
            print("add:",add)



        #print "_scale_output:",self._scale_output
        #print "atoms:",self._atoms
        self._outfile = "E"+self.equation[:1].upper()+self.equation[1:].lower()+str(add)
        print("writing",self._outfile)
        print("e0:",self.e0)
        print("v0:",self.v0)
        print("b0:",self.b0)
        print("b0der:",self.b0der)
        if self.e0 != None and self.b0 != None and self.b0der != None:
            np.savetxt(
            folder + "/" + self._outfile,
            [[self.e0*self._scale_output, self.v0*self._scale_output, self.b0, self.b0der]],
            fmt='%.12f')
        else:
            np.savetxt(
            folder + "/" + self._outfile,
            [self.parameters],fmt='%.12f')


        if False:
            np.savetxt(
                    file2,
                    np.transpose([self.data[:, 0], self.fitdeltas]),
                    fmt='%.18f',
                    delimiter='  ')   # X is an array

        ###################################################
        # write energy.dat
        ###################################################
        if self._write_energydat:
            folder = os.path.split(os.path.realpath(self.inputfile))[0]
            filename="energy.dat"
            np.savetxt(folder + "/" + filename,
                    np.transpose([self.data[:,0],self.data[:,1]]),
                    fmt='%.18f',
                    delimiter='  ')   # X is an array
        ###################################################
        # write energy_peratom.dat and energy_peratom.dat.fit
        ###################################################
        if self._write_energyperatomdat:
            folder = os.path.split(os.path.realpath(self.inputfile))[0]
            filename="energyperatom.dat"
            np.savetxt(folder + "/" + filename,
                    np.transpose([self.data[:,0],self.data[:,1]]),
                    fmt='%.18f',
                    delimiter='  ')   # X is an array
            filename="energyperatom.dat.fit"
            np.savetxt(
                    folder + "/" + filename,
                    np.transpose([self.fitvolumes, self.fitenergies]),
                    fmt='%.18f',
                    delimiter='  ')   # X is an array

        #print('sd',self.data)
        #print('fd',self.fitdeltas)
        plotdiff = np.transpose([self.data[:,0],self.fitdeltas])
        np.savetxt("delta_evinet_fit-data.dat",plotdiff)
        return

    def plot(self, plot=2):
        import matplotlib
        import matplotlib.pyplot as plt

        #plt.ion()
        #plt.ion()
        #plt.cfg()
        plt.rc('text' , usetex=True)
        plt.rc('font' , family='serif')

        if plot is 1: # only energyes vs voluem plot
            plt.grid(True)  # Trunes on grid
            plt.plot(self.data[:, 0], self.data[:, 1], 'ro', self.fitvolumes, self.fitenergies)
        if plot is 2: # energes vs volume and deltaenergies vs volume plot
            # this defines the matplotlib backend
            matplotlib.use('gtkagg')



            # create all axes we need
            ax0 = plt.subplot(211)
            plt.ylabel('Energy [meV]')
            ax2 = plt.subplot(212)
            plt.ylabel('delta to fit [meV]')
            plt.xlabel('Volume per atom [\AA$^3$]')
            #plt.xlabel('Volume per atom [\AA$^3$]')
            #ax.set_xlabel('$a_0$ (\AA)')
            #plt.set_xlabel('$a_{0}$ (\AA)')
            # share the secondary axes
            ax0.get_shared_x_axes().join(ax0, ax2)

            ax0.plot(self.data[:, 0], self.data[:, 1], 'r.', self.fitvolumes, self.fitenergies)
            ax2.plot(self.data[:, 0], self.fitdeltas , 'o-')

            ax0.grid(True)  # Trunes on grid
            ax2.grid(True)  # Trunes on grid
            plt.show()
        return


class qht(object):
    def __init__(self):
        self.ka = 1

    def help(self, p = None):
        string = '''
        new
        '''
        if p == None:
            from argparse import ArgumentDefaultsHelpFormatter
            #from argparse import RawDescriptionHelpFormatter
            #from argparse import RawTextHelpFormatter
            #p = argparse.ArgumentParser(description=string, formatter_class=RawTextHelpFormatter)
            p = argparse.ArgumentParser(description=string,
                    formatter_class=ArgumentDefaultsHelpFormatter)
                    #formatter_class=RawTextHelpFormatter) # No
                    #formatter_class=RawDescriptionHelpFormatter)  # No
        p.add_argument('-r',   '--rrr', action='store_true',
           help='specify equation of state', default=False)
        #p.add_argument('-v', '--verbose',
        #        help='verbose', action='store_true', default=False)
        return p

if __name__ == '__main__':
    p = eos().help()  # this gives the possibility to change some __init__ settings
    p = qht().help(p)  # this gives the possibility to change some __init__ settings
    args = p.parse_args()

    eos = eos(args)

    #import data and export (just in case several _inputfiles are defined)
    for inputfile in eos._inputfiles:
        eos.inputfile = inputfile
        if os.path.isfile(inputfile) != True:
            args.c = True
        if args.c:
            eos.inputfile = create_energy_dat_from_underlying_OUTCARS()
        eos.import_energy_vs_volume_data()
        eos.fit_to_energy_vs_volume_data(eos.datax, eos.datay)
        eos.write_data()
        _printgreen("DONE!")
        print("")
    #eos.plot()
    myutils.create_READMEtxt()
