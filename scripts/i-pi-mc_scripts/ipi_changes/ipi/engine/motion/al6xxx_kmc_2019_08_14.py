"""Contains the classes that deal with the different dynamics required in
different types of ensembles.

Holds the algorithms required for normal mode propagators, and the objects to
do the constant temperature and pressure algorithms. Also calculates the
appropriate conserved energy quantity for the ensemble of choice.
"""

# This file is part of i-PI.
# i-PI Copyright (C) 2014-2015 i-PI developers
# See the "licenses" directory for full license information.


from __future__ import print_function
from ase.build import bulk as ase_build_bulk
import time,os,sys
from ase.io import write as ase_write
from collections import defaultdict
import pickle
import threading
#import logging
import numpy as np
from ipi.utils.softexit import softexit
from ase.calculators.lammpslib import LAMMPSlib
from ase.optimize import LBFGS
from ase.optimize import GPMin
import timeit,time

from ipi.engine.motion import Motion, GeopMotion
from ipi.utils.depend import dstrip, depend_value, dobject, dpipe, dd
from ipi.engine.thermostats import Thermostat
from ipi.engine.cell import Cell
from ipi.engine.normalmodes import NormalModes
from ipi.engine.barostats import Barostat
from ipi.utils.units import Constants
import ipi.utils.io as io


def get_new_state(state,svac,sneigh):
        nstate = state.copy() # ['S' 'S' 'S' 'M' 'M' 'M' 'A' 'A' 'A' 'A' ....'V' 'A' 'A' ]
        nstate[svac], nstate[sneigh] = state[sneigh], state[svac]
        nstr = "".join(nstate) # this is the string that corresponds to the new state
        return nstate,nstr


def ase_minimize(atomsc_in=False,minimizer="LBFGS"):
    ''' minimizer: LBFGS,GPmin
        currently: expects positions in angstrom, returns positions in bohrradius
    '''
    if minimizer.lower() not in ["lbfgs", "gpmin"]:
        sys.exit("minimizer not known")
    #print('ASE -1',atomsc_in.positions.shape)
    #print('ASE -2',len(atomsc_in.positions))
    atomsc = atomsc_in.copy()
    if False:
        for idx,i in enumerate(atomsc.positions):
            if idx < 8:
                print('ASE pos:',idx,atomsc.get_chemical_symbols()[idx],atomsc.positions[idx])

    del atomsc[-1]  # delete the V(acancy)
    # get calculator
    atom_types = {'Mg':1,'Al':2,'Si':3}
    lmpcmd = ['mass 1 24.305', 'mass 2 26.9815385', 'mass 3 28.0855', 'variable nnpDir string "/Users/glensk/Dropbox/Albert/scripts/dotfiles/scripts/potentials/runner_v3ag_5000_46489_2"', 'pair_style runner dir ${nnpDir} showewsum 1 showew yes resetew no maxew 1000000', 'pair_coeff * * 14.937658735']
    asecalcLAMMPS = LAMMPSlib(lmpcmds=lmpcmd, atom_types=atom_types,keep_alive=True)
    atomsc.set_calculator(asecalcLAMMPS)
    #ene = atomsc.get_potential_energy()
    #print('ene in (eV)',ene,"in hartree:",ene*0.036749322)
    #print('ene in  (ASE)',ene*0.036749322,"(hartree)")

    #starttime = timeit.timeit()
    #starttime = time.time()
    #fmax=0.03
    fmax=0.09 # this seems to be closeset to the ipi geop settings
    logfile="-" # output to screen
    #logfile="XXLMPLOG" # output to file and not to screen
    #logfile=False
    logfile=None  # no output at all
    if minimizer.lower() == "lbfgs":
        opt = LBFGS(atomsc,logfile=logfile)
    elif minimizer.lower() == "gpmin":
        opt = GPMin(atomsc,logfile=logfile)
    opt.run(fmax=fmax)
    #print('pos-out')
    #print(atomsc.positions[:6])
    ene = atomsc.get_potential_energy()
    #print('ene out (eV)',ene,"in hartree:",ene*0.036749322)
    #endtime = timeit.timeit()
    #endtime = time.time()
    #print('ene out (ASE)',str(ene*0.036749322).ljust(25),"(hartree) in",endtime - starttime,"seconds")
    #positions_out = atomsc.positions
    #print('positions_out ASE')
    #print(positions_out)
    #print('ASE 1',positions_out.shape)
    #print('ASE 2',len(positions_out))
    positions_out = np.append(atomsc.positions, [atomsc_in.positions[-1]], axis=0)
    if True:
        atomsc_out = atomsc_in.copy()
        atomsc_out.set_positions(positions_out)
    #print('ASE 4',positions_out.shape)
    #print('ASE 5',len(positions_out))
    return ene*0.036749322,positions_out.flatten()/0.52917721,atomsc_out

def gs(NN1,search="S"):
    ''' sums elements in array and return frequency '''
    return (NN1==search).sum()

def str_to_symbols(string):
    sym = []
    for i in string:
        if i == 'A':
            sym.append("Al")
        elif i == 'S':
            sym.append("Si")
        elif i == 'M':
            sym.append("Mg")
        else:
            sys.exit('not A S or M')
    return sym

def get_cliques(pairs):
    ''' Build a graph using the pairs '''
    nodes = defaultdict(lambda: [])
    for a, b in pairs:
        if b is not None:
            nodes[a].append((b, nodes[b]))
            nodes[b].append((a, nodes[a]))
        else:
            nodes[a]  # empty list

    # Add all neighbors to the same group
    visited = set()
    def _build_group(key, group):
        if key in visited:
            return
        visited.add(key)
        group.add(key)
        for key, _ in nodes[key]:
            _build_group(key, group)
    groups = []
    for key in nodes.keys():
        if key in visited: continue
        groups.append(set())
        _build_group(key, groups[-1])
    return groups


def load_cache_file(ecache_file):
    if not os.path.isfile(ecache_file):
        print(ecache_file,"does not exist! ")
        return {}

    print('>> (7) loading',ecache_file,"ECACHE/QCACHE ~10/80sec for 0.2GB/9GB file")
    f = open(ecache_file, "rb")
    start_time = time.time()
    x = pickle.load(f)
    #z = x.copy()
    #z = x
    #print('loaded IN :',x)
    #lst = []
    while 1:
        try:
            y = pickle.load(f)
            x.update(y)
            #lst.append(y)
            #print('loaded in :',y)
        except EOFError:
            #f.close()
            break
    #print('loaded len   :',len(lst),lst)
    #print('loaded FIN[:]:',lst[:])
    #out = lst[:len(lst)]
    #print('loaded FIN[0]:',out)
    #print('x',x)
    print('>> (7) loading',ecache_file,'done. In',time.time() - start_time,'seconds, containing',len(x),'structures')
    return x

class AlKMC(Motion):
    """Stepper for a KMC for Al-6xxx alloys.

    Gives the standard methods and attributes needed in all the
    dynamics classes.

    Attributes:
        beads: A beads object giving the atoms positions.
        cell: A cell object giving the system box.
        forces: A forces object giving the virial and the forces acting on
            each bead.
        prng: A random number generator object.
        nm: An object which does the normal modes transformation.

    Depend objects:
        econs: The conserved energy quantity appropriate to the given
            ensemble. Depends on the various energy terms which make it up,
            which are different depending on the ensemble.he
        temp: The system temperature.
        dt: The timestep for the algorithms.
        ntemp: The simulation temperature. Will be nbeads times higher than
            the system temperature as PIMD calculations are done at this
            effective classical temperature.
    """

    def __init__( self, mode, geop, nstep, a0, ncell, nvac, nsi, nmg,
                  neval, diffusion_barrier_al, diffusion_prefactor_al,
                  diffusion_barrier_mg, diffusion_prefactor_mg,
                  diffusion_barrier_si, diffusion_prefactor_si,
                  idx=[], tottime=0, ecache_file="",
                  qcache_file="", thermostat=None, barostat=None,
                  fixcom=False, fixatoms=None, nmts=None):
        """Initialises a "dynamics" motion object.

        Args:
            dt: The timestep of the simulation algorithms.
            fixcom: An optional boolean which decides whether the centre of mass
                motion will be constrained or not. Defaults to False.
        """


        # This will generate a lattice model based on a primitive FCC cell. the lattice is represented in three ways:
        # 1. as a string in which each lattice site is identified by a letter
        # 2. by a list of the lattice sites in 3D space, in 1-1 mapping with the letters
        # 3. by a list of the atoms, whose lattice position is indicated by an integer


        self.nstep = nstep
        self.ncell = ncell
        self.nvac = nvac
        self.nsi = nsi
        self.nmg = nmg
        self.nsites = self.ncell**3

        self.ase            = True  # ase yes? ase_yes_no
        self.normal         = True  # use normal evaluation of energy
        self.filled         = True  # use evaalutation of energy where beyond cutoff only Al
        self.conventional   = False  # use qubic supercell?
        self.parallel       = True  # run with threading
        self.write_pos      = True   # wrte in_pos_ && out_pos
        self.write_strings  = False  # write AASMAAAAVAASAAA ....
        self.step_isel      = False #[ 1,10,6,8,2,7,4,3,9,5,11 ] # or False
        self.cutoff_filled  = 1.1  # 1.1 == 1NN +2NN = 18 atoms
        self.cutoff_filled  = 1.5  # 1.5 == 34 atoms  (??? in 7x7x7) (54 in 5x5x5)
        self.cutoff_filled  = 2.0  # 1.5 == 52 atoms  (140 in 7x7x7) (98 in 5x5x5)
        self.diff_mev_normal = []
        self.diff_mev_filled = []
        self.diff_diff_abs_mev_filled = []
        self.ediff_vs_solutes = []

        if self.ase == False and self.filled == True:
            sys.exit("if you want self.filled you need to ues ase!")

        if self.conventional:
            self.ncell = 2
            self.nsites = 4*self.nsites


        self.natoms = self.nsites - self.nvac
        self.neval = neval
        self.diffusion_barrier_al = diffusion_barrier_al
        self.diffusion_prefactor_al = diffusion_prefactor_al
        if diffusion_barrier_mg > 0:
            self.diffusion_barrier_mg = diffusion_barrier_mg
        else:
            self.diffusion_barrier_mg = diffusion_barrier_al
        if diffusion_barrier_si > 0:
            self.diffusion_barrier_si = diffusion_barrier_si
        else:
            self.diffusion_barrier_si = diffusion_barrier_al
        if diffusion_prefactor_mg > 0:
            self.diffusion_prefactor_mg = diffusion_prefactor_mg
        else:
            self.diffusion_prefactor_mg = diffusion_prefactor_al
        if diffusion_prefactor_si > 0:
            self.diffusion_prefactor_si = diffusion_prefactor_si
        else:
            self.diffusion_prefactor_si = diffusion_prefactor_al
        self.barriers = { "A": self.diffusion_barrier_al,
                          "M": self.diffusion_barrier_mg,
                          "S": self.diffusion_barrier_si }
        self.prefactors = { "A": self.diffusion_prefactor_al,
                          "M": self.diffusion_prefactor_mg,
                          "S": self.diffusion_prefactor_si }

        self.a0 = a0
        cell=np.zeros((3,3))
        cell[0]=[0.7071067811865475, 0.35355339059327373, 0.35355339059327373]
        cell[1]=[0.,0.6123724356957945, 0.20412414523193154]
        cell[2]=[0.,0.,0.5773502691896258]

        if self.conventional:
            cell[0] = [1,0,0]
            cell[1] = [0,1,0]
            cell[2] = [0,0,1]

        self.scell=self.a0*cell
        self.dcell = Cell()
        self.dcell.h = self.scell*self.ncell


        print("LATTICE PARAM (self.a0)", self.a0)
        # this is the list of lattice sites, in 3D coordinates
        ix,iy,iz = np.meshgrid(range(self.ncell), range(self.ncell), range(self.ncell), indexing='ij')
        self.sites = np.dot(np.asarray([ix.flatten(),iy.flatten(),iz.flatten()]).T, self.scell.T)

        if self.conventional:
            fcc2 = "/Users/glensk/Thermodynamics/utilities/fcc/EqCoords_direct_fcc_2x2x2sc"
            self.sites = np.loadtxt(fcc2)*(self.a0*0.52917721*ncell)/0.52917721



        self.neigh  = self.determing_connectivity_of_the_fcc_lattice(cutoff = self.a0*8./9.,verbose=False)
        print()
        self.neigh2 = self.determing_connectivity_of_the_fcc_lattice(cutoff = self.a0) # ,mindist=self.a0/2.)
        print()
        #self.neighc = self.determing_connectivity_of_the_fcc_lattice(cutoff = self.a0*1.5)
        #self.neighc = self.determing_connectivity_of_the_fcc_lattice(cutoff = self.a0*2.5)
        #self.neighc = self.determing_connectivity_of_the_fcc_lattice(cutoff = self.a0*0.9,verbose=False)  # 1NN    12
        #self.neighc = self.determing_connectivity_of_the_fcc_lattice(cutoff = self.a0*1.1,verbose=False)  # (1+2)NN 15
        self.neighc = self.determing_connectivity_of_the_fcc_lattice(cutoff = self.a0*self.cutoff_filled,verbose=False)  # (1+2)NN 15
        print('sn1',self.neigh.shape,self.neigh.shape[1])
        print('sn2',self.neigh2.shape)
        print('snc',self.neighc.shape)

        self.atomsc         = self.set_up_ase_structure()
        self.atomsc_filled  = self.set_up_ase_structure(alloy=False)
        self.atomsc_neigh   = self.set_up_ase_structure(atoms=self.neigh.shape[1])
        self.atomsc_neigh2  = self.set_up_ase_structure(atoms=self.neigh2.shape[1])
        self.atomsc_neighc  = self.set_up_ase_structure(atoms=self.neighc.shape[1])

        print('self.neigh[-1]',self.neigh[-1],'self.neigh.shape',self.neigh.shape)
        #print('self.neigh[63]',self.neigh[63])
        #print('self.neigh[64]',self.neigh[64])

        #sys.exit()
        #self.atomsc_neigh.set_positions(self.sites[self.NN1idx]*0.52917721)


        #print len(self.sites), self.nsites, "###"
        ## now we build list of nearest neighbors (fcc-lattice hardcoded!)
        #self.neigh = np.zeros((self.nsites,12),int)
        #nneigh = np.zeros(self.nsites,int)
        ## could be done in a more analytic way but whatever, I'm too lazy
        #a02 = 1.01*0.5*self.a0**2                                        # perhaps 1.01 it is not enough, must check!
        #for i in xrange(self.nsites): # determines the connectivity of the lattice
        #    rij = self.sites.copy().flatten()
        #    for j in xrange(self.nsites):
        #        rij[3*j:3*j+3] -= self.sites[i]
        #    self.dcell.array_pbc(rij)
        #    rij.shape = (self.nsites,3)
        #    for j in xrange(i):
        #        if np.dot(rij[j],rij[j]) < a02: # found nearest neighbor
        #            self.neigh[i,nneigh[i]] = j
        #            self.neigh[j,nneigh[j]] = i
        #            nneigh[i]+=1
        #            nneigh[j]+=1

        self.idx = idx
        print('self.idx',self.idx)


        # the KMC step is variable and so it cannot be stored as proper timing
        dd(self).dt = depend_value(name="dt", value = 0.0)
        self.fixatoms = np.asarray([])
        self.fixcom = True
        self.geop = [None] * self.neval
        # geop should not trigger exit if there is early convergence, but just carry on.
        # we hard-code this option to avoid early-termination that would be hard to debug for a user
        geop["exit_on_convergence"] = False
        for i in xrange(self.neval):
            # geometry optimizer should not have *any* hystory dependence
            self.geop[i] = GeopMotion(fixcom=fixcom, fixatoms=fixatoms,**geop) #mode="cg", ls_options={"tolerance": 1, "iter": 20,  "step": 1e-3, "adaptive": 0.0}, tolerances={"energy": 1e-7, "force": 1e-2, "position": 1e-4}, ) #!TODO: set the geop parameters properly

        # dictionary of previous energy evaluations.
        self.ecache_file = ecache_file
        self.qcache_file = qcache_file

        #self.try_load_ECACHE_QCACHE_files()
        self.ecache = load_cache_file(ecache_file)
        self.qcache = load_cache_file(qcache_file)
        self.ecache_n = {}
        self.qcache_n = {}

        #try:
        #    ff = open(self.ecache_file, "rb")
        #    self.ecache = pickle.load(ff)
        #    ff.close()
        #    ff = open(self.qcache_file, "rb")
        #    self.qcache = pickle.load(ff)
        #    ff.close()
        #    print "Loaded %d cached energies" % (len(self.ecache))
        #except:
        #    print "Couldn't load cache files "+self.ecache_file+","+self.qcache_file+" - resetting"
        #    self.ecache = {}
        #    self.qcache = {}

        self.ncache = len(self.ecache)
        self.ncache_stored = self.ncache


        # no TS evaluation implemented yet
        self.tscache = {}

        # todo make these optional and restarted
        self.kmcfile = open("KMC_AL6XXX","w+")
        self.tottime = tottime

    def print_pos_ipi(self,text="",first=None):
        pos_tmp = self.dbeads[0].q[0,:]
        pos_tmp.shape = (self.nsites,3)
        pos_tmp = pos_tmp*0.52917721
        for idx,i in enumerate(pos_tmp[:first]):
            print(text,idx,pos_tmp[idx])
        return

    def set_up_ase_structure(self,atoms=False,alloy=True,cubic=False):
        atom = ase_build_bulk("Al",crystalstructure='fcc',a=self.a0,cubic=self.conventional)
        atomsc = atom.repeat(self.ncell)
        if alloy == True:
            for i in np.arange(self.nsi):  atomsc[i].symbol = 'Si'
            for i in np.arange(self.nsi,self.nmg+self.nsi): atomsc[i].symbol = 'Mg'
            for i in np.arange(1,self.nvac+1)*-1:    atomsc[i].symbol = 'V'

        atomsc.set_cell(self.dcell.h.T*0.52917721)
        atomsc.set_positions(self.sites*0.52917721)
        print('numat',atomsc.get_number_of_atoms())
        if atoms != False:
            while atomsc.get_number_of_atoms() > atoms:
                del atomsc[-1]
            #for idx,i in enumerate(np.arange(atomsc.get_number_of_atoms())[::-1]):
            #    print('idx',idx,i,atomsc.positions[idx])
            #    if i > atoms-1:
            #        print('atoms',atoms,'deleting atom ',i)
            #        del atomsc[-1]
            #        print('remaining',atomsc.get_number_of_atoms())
            print()
            for idx,i in enumerate(np.arange(atomsc.get_number_of_atoms())):
                print('idx',idx,i,atomsc.positions[idx])
            print()
        return atomsc

        def get_dbeads_names(self,ieval,state,filled = False,state_filled=False):
            ''' for unfilled states
                state = self.state (for initial configurations)
                state = nstate (for swiched cofigurations)

                if filled == True also define:
                state_filled = state_filled (for initial configuration) or
                state_filled = nstate_filled '''
            unique_idx       = self.unique_idx(state)
            self.dbeads[ieval].q[0,:] = self.sites[unique_idx].flatten() # CORRECT
            self.dbeads[ieval].names = self.dbeads_names
            if filled == True:
                for solute_idx in range(self.nsi + self.nmg): # is already the index
                    if state_filled[unique_idx][solute_idx] == 'A':
                        (self.dbeads[ieval].names)[solute_idx] = "Al"
            return


    def showpos_around(self,step,svac):
        ''' svac is the index around which the neighbors will be shown '''
        # a) first do the reference (fill up around every vacancy)
        print('step',step,'svac',svac,"bohr    :",self.sites[svac])
        print('step',step,'svac',svac,"angstrom:",self.sites[svac]*0.52917721)

        NNcidx = self.neighc[svac]                              # atoms in the sphere
        print('neighbors')
        print('step',step,'svac',svac)
        for idx,i in enumerate(NNcidx):
            print(self.state[NNcidx[idx]],self.sites[NNcidx[idx]]*0.52917721)
        tmp_neigh = (self.sites[NNcidx]-self.sites[svac])
        print()
        print('neighbors distances')
        for tnidx,tmp_neigh_v in enumerate(tmp_neigh):
            for tnxyz,xyz in enumerate(tmp_neigh_v):
                if xyz > 1.05*(self.a0*self.ncell/2):
                    tmp_neigh[tnidx,tnxyz] = xyz-self.a0*self.ncell
                if xyz < -1.05*(self.a0*self.ncell/2):
                    tmp_neigh[tnidx,tnxyz] = xyz+self.a0*self.ncell
        print('ORIG STRUCT: step',step,'svac',svac,'diff')
        for idx,i in enumerate(tmp_neigh):
            print(self.state[NNcidx[idx]],tmp_neigh[idx]*0.52917721)
        return


    def fill_state_with_Al(self,idx1,idx2,state,svac_idx):
        ''' idx1: has to be idx1 (a vacancy)
            state: can be self.state or nstate or another state
        '''
        NNcidx_v    = self.neighc[idx1]                              # atoms in the sphere
        NNcidx_n    = self.neighc[idx2]                              # atoms in the sphere
        NNcidx      = np.union1d(NNcidx_n,NNcidx_v)
        NNcidx_not  = np.delete(np.arange(len(self.idx)),NNcidx) # remaining atoms
        state_filled = state.copy()
        state_filled[NNcidx_not] = 'A'
        state_filled[svac_idx] = 'V'
        ostr_filled = "".join(state_filled)
        ## up to here it is enough when in cache
        dbead_names_filled = self.dbeads_names.copy()
        for solute_idx in range(self.nsi + self.nmg): # is already the index
            if state_filled[self.unique_idx(state)][solute_idx] == 'A':
                dbead_names_filled[solute_idx] = "Al"
        return state_filled,ostr_filled,dbead_names_filled



    def enumerate_possible_reaction_events(self,ostr,step=False,filled=False,filled_ref=False,verbose=False):
        '''
        one particular nevent = [svac, sneigh, self.ecache[nstr], self.qcache[nstr], 0.0]
        one particular nevent = [svac, sneigh, self.ecache[nstr], self.qcache[nstr], ostr, nstr, 0.0]
        '''

        #format = "%(asctime)s: %(message)s"
        #logging.basicConfig(format=format, level=logging.INFO, datefmt="%H:%M:%S")
        #indexlog = 0
        # enumerates possible reactive events (vacancy swaps)
        levents = []
        self.in_cache = 0
        self.not_in_cache = 0
        self.in_cache_filled = 0
        self.not_in_cache_filled = 0
        ethreads = [None] * self.neval  # defined in input.nn (posible parallel events)
        # loops over the vacancy
        #print('maxdiff+-:',0.95*self.a0*self.ncell*0.52917721)
        for ivac in xrange(self.natoms, self.natoms + self.nvac):
            svac = self.idx[ivac] # lattice site associated with this vacancy

            #print("######## ORIGINAL POS",'svac',svac,"#######################")
            #self.showpos_around(step,svac)
            #print("######## ORIGINAL POS done ###################")
            # a) first do the reference (fill up around every vacancy)
            # b) do swap from original state
            # first do the swap, then fill up!  --> but, what is with the reference??!!

            if self.state[svac] != "V":
                raise IndexError("Something got screwed and a vacancy state is not a vacancy anymore!")

            # loops over the neighbors of the selected vacancy


            for isel_sidxx,sneigh in enumerate(self.neigh[svac]):
                #print('isel_sidxx',isel_sidxx,'svac',svac,'sneigh',sneigh)
                # if the neighbor is a vacancy, move on. does not make sense to swap two vacancies!
                if self.state[sneigh] == "V" : continue

                # a) swop state (normal)
                nstate,nstr = get_new_state(self.state,svac,sneigh)  # is necessary in any case
                nstate_ref = nstate.copy()  # this is necessary to set later on the swaped positions in dbeads


                if filled == True or filled_ref == True:
                    ## this loop needs also to be parallel
                    #######################################################################
                    # creates nstate, nstr for filled state
                    # reference structure (filled), and ---contrary the the original
                    # structure--- needs to be calculated for every swap
                    # but, only the dbead_names need to be chagned, not the positions
                    #######################################################################
                    state_filled, ostr_filled,dbead_names_filled = self.fill_state_with_Al(svac,sneigh,self.state,svac)
                    # this is the reference
                    #nstate, ostr_filled,dbead_names_filled = self.fill_state_with_Al(svac,sneigh,self.state,svac)
                    ostr = ostr_filled
                    nstr = ostr_filled

                    # This needs to be done in a separate loop of "enumerate_possible_reaction_events(...)"
                    if filled == True:
                        # check if referece is available
                        if not ostr_filled in self.ecache:
                            softexit.trigger("Error: not calculated reference ...! Exit.")

                        #if False:
                        #    # make sure to calculate this only when filled_ref == True
                        #    if not ostr_filled in self.ecache:
                        #        self.not_in_cache_filled += 1
                        #        # self.state is the unswapped state
                        #        # nstate bzw. nstate_ref are the swapped state
                        #        self.dbeads[ieval].q[0,:] = self.sites[self.unique_idx(self.state)].flatten()
                        #        #self.dbeads[ieval].names = dbead_names_filled
                        #        rv = [0, 0, 0, 0, 0]
                        #        savename = "pos_step_"+str(step)+"_filled_"+str(svac)+"_"+str(sneigh)+"_normal"
                        #        # nope, calculate this further down!
                        #        #self.geop_thread(ieval, ostr_filled, rv, None,savename=savename)
                        #    else:
                        #        self.in_cache_filled += 1



                    #######################################################################
                    ## creates (filled) swapped sate
                    #######################################################################
                    if filled == True:
                        nstate, nstr, dbead_names_filled = self.fill_state_with_Al(svac,sneigh,nstate,sneigh)



                if not nstr in self.ecache:
                    print('calc not in cache')
                    if filled == True:
                        self.not_in_cache_filled += 1
                    else:
                        self.not_in_cache += 1

                    if  self.ridx[svac]!= ivac or self.idx[ivac]!=svac:
                        raise  IndexError("Something got screwed and the index does not correspond anymore to site occupancy")

                    ###################################
                    ## stuff for geop
                    ###################################
                    ieval = self.find_eval(ethreads)

                    if filled == False and filled_ref == False: # normal case
                        self.dbeads[ieval].q[0,:] = self.sites[self.unique_idx(nstate_ref)].flatten()
                        #self.dbeads[ieval].q[0,:] = self.sites[self.unique_idx(nstate)].flatten()
                        self.dbeads[ieval].names = self.dbeads_names

                    # filled case
                    if filled == True or filled_ref == True:
                        if filled == True:  # swapped pos
                            # here, the positions are the same as in nstate (no need to change dbeads.q
                            # and only the atom names change
                            self.dbeads[ieval].q[0,:] = self.sites[self.unique_idx(nstate_ref)].flatten()
                            self.dbeads[ieval].names = dbead_names_filled
                        elif filled_ref == True:
                            self.dbeads[ieval].q[0,:] = self.sites[self.unique_idx(self.state)].flatten()
                            self.dbeads[ieval].names = dbead_names_filled
                            ostr = ostr_filled
                            nstr = ostr_filled
                            print('fill ref ...')


                    # nevent = [svac, sneigh, ene, pos, 0.0]
                    # nevent = [svac, sneigh, 0.0, 0.0, 0.0]
                    nevent   = [svac, sneigh, 0.0, 0.0, ostr, nstr, 0.0]

                    # runs a geometry optimization
                    savename = "pos_step_"+str(step)+"_swap_"+str(svac)+"_"+str(sneigh)+"_normal"

                    if filled == True:
                        savename = "pos_step_"+str(step)+"_filled_"+str(svac)+"_"+str(sneigh)+"_swap"
                    if filled_ref == True:
                        savename = "pos_step_"+str(step)+"_filled_"+str(svac)+"_"+str(sneigh)+"_normal"

                    if self.parallel != True:
                        self.geop_thread(ieval=ieval, nstr=nstr, nevent=nevent,ostr=None,savename=savename)
                    else:
                        #logging.info("Main    : create and start thread %d.", indexlog)
                        #indexlog+=1
                        if step == 1 and sneigh == 20:
                            print('praallel nstr', nstr)
                        st = threading.Thread(target=self.geop_thread, name=str(ieval), kwargs={"ieval":ieval, "nstr":nstr, "nevent" : nevent, "ostr": ostr, "savename": savename})
                        #print('starting thread',sneigh)
                        st.daemon = True
                        st.start()
                        ethreads[ieval] = st
                    ###################################
                    ## stuff for geop stop
                    ###################################


                else:
                    if filled == True:
                        self.in_cache_filled += 1
                    else:
                        self.in_cache += 1
                    #print('nstr (swap) normal/filled :)')
                    #print "Found state ", nstr, " retrieving cached energy ", self.ecache[nstr]

                    # fetch energy from previous calculation
                    #nevent = [svac, sneigh, self.ecache[nstr], self.qcache[nstr], 0.0]
                    nevent = [svac, sneigh, self.ecache[nstr], self.qcache[nstr], ostr, nstr, 0.0]

                    # EVALUATION OF TS ENERGY IS DISABLED FOR THE MOMENT...
                    # we might still need to compute the TS energy!
                    # if not nstr in self.tscache[ostr]:
                    #   print "Computing TS"
                    #   ieval = self.find_eval(ethreads)
                    #    st = threading.Thread(target=self.ts_thread, name=str(ieval), kwargs={"ieval":ieval, "ostr": ostr, "nstr":nstr, "nevent" : nevent})
                    #    st.daemon = True
                    #    st.start()
                    #    ethreads[ieval] = st
                    #else:
                    #    print "Found TS"
                    #    nevent[3] = self.tscache[ostr][nstr]
                nevent[-1] = self.state[sneigh]  # this is just an "A" or an "S" or an "M"
                levents.append( nevent )
        # wait for all evaluators to finish
        if self.parallel == True:
            print("parallel: wait for all evaluators to finish")
            for st in ethreads:
                #logging.info("Main    : before joining thread %d.", indexlog)
                while not st is None and st.isAlive():
                    st.join(2)
                    #logging.info("Main    : thread %d done", indexlog)
        print("all evaluators done")
        return levents


    def print_pos_ase(self,text="",first=None):
        for idx,i in enumerate(self.atomsc.positions[:first]):
            print(text,idx,self.atomsc.positions[idx]*0.52917721)
        return

    def determing_connectivity_of_the_fcc_lattice(self,cutoff,mindist=0.,verbose=False):
        ''' cutoff is the cutoff like self.a0 (which would include 2NN in a conventional fcc
        unit cell
        '''
        if verbose:
            print(">> (3) determines the connectivity of the lattice ...")
        maxatoms = 999999
        selfneigh  = np.ones((self.nsites,maxatoms),int)*-1
        nneigh  = np.zeros(self.nsites,int)
        print('cutoff/self.a0',(cutoff/self.a0))
        cutoff_ = 1.01*(cutoff/self.a0)*self.a0**2                                        # perhaps 1.01 it is not enough, must check!
        cutoff_ = 1.01*((cutoff/self.a0)*self.a0)**2                                        # perhaps 1.01 it is not enough, must check!
        print('cutoff_',cutoff_)
        mindist_ = 1.01*((mindist/self.a0)*self.a0)**2                                        # perhaps 1.01 it is not enough, must check!
        print('mindist_',mindist)
        #a02     = 1.01*(0.5)*self.a0**2                                        # perhaps 1.01 it is not enough, must check!
        #print("1NNs",selfneigh[0])
        #print(nneigh)
        for i in xrange(self.nsites): # determines the connectivity of the lattice
            neighline = []
            #print('i',i)
            rij = self.sites.copy().flatten()
            for j in xrange(self.nsites):
                rij[3*j:3*j+3] -= self.sites[i]
            self.dcell.array_pbc(rij)
            rij.shape = (self.nsites,3)
            for j in xrange(i):
                #print('i',i,'j',j)
                dist = np.dot(rij[j],rij[j])
                if verbose: # and i == 31:
                    print('i',i,'j',j,dist,self.sites[i])
                #if dist < a02: #cutoff_: # found nearest neighbor
                if dist > mindist_ and dist < cutoff_: #a02: #cutoff_: # found nearest neighbor
                    #print(i,'j',j,np.dot(rij[j],rij[j]))
                    selfneigh[i,nneigh[i]] = j
                    selfneigh[j,nneigh[j]] = i
                    nneigh[i]+=1
                    nneigh[j]+=1
                    neighline.append(j)
            #if verbose:
            #    sys.exit('asdf;ff')
        if verbose:
            print(">> (4) 1NNs (of atom at idx=0):",selfneigh[0])
            print('>> (4) nneigh',nneigh)
        amount_neighbors_max = np.where(selfneigh[0] > 0)[0].max()+1
        return selfneigh[:,:amount_neighbors_max]

    def barriers_energies_analysis(self,step,levents,rates):
        if step is not None:  # only than we want to write
            if not hasattr(self, 'file_KMC_barriers'):
                self.file_KMC_barriers = open("KMC_barriers","a+")
            if not hasattr(self, 'file_KMC_rates'):
                self.file_KMC_rates = open("KMC_rates","a+")
            if not hasattr(self, 'file_KMC_energies'):
                self.file_KMC_energies = open("KMC_energies","a+")
            if not hasattr(self, 'file_KMC_denergies'):
                self.file_KMC_denergies = open("KMC_denergies","a+")

            self.file_KMC_energies.write("%18.11e"% (self.ecurr))
            for i in xrange(len(levents)):
                self.file_KMC_energies.write("%20.11e"% (levents[i][2]))
                self.file_KMC_barriers.write("%20.11e"% (self.barriers[levents[i][-1]]))
                self.file_KMC_rates.write("%20.11e"% rates[i])
                self.file_KMC_denergies.write("%20.11e"% (-0.5*(self.ecurr - levents[i][2])))

            for i in [self.file_KMC_energies,self.file_KMC_barriers,self.file_KMC_rates,self.file_KMC_denergies]:
                i.write("\n")
                i.flush()
        return

    def vacancy_neighborhood_analysis(self):
        ''' can be performed directy ins step without any further calculations '''
        # loops over the vacancies
        for ivac in xrange(self.natoms, self.natoms + self.nvac):
            svac = self.idx[ivac] # lattice site associated with this vacancy

            # this are the indices of the vacancy neighbors
            self.NN1idx = self.neigh[svac]
            self.NN2idx = self.neigh2[svac]
            self.NNcidx = self.neighc[svac]

            #print('svac',svac,'self.NN1idx',self.NN1idx,type(self.NN1idx))
            #print('svac',svac,'self.NNcidx',self.NNcidx,type(self.NNcidx))
            #print('svac',svac,'self.idx.shape',self.idx.shape)
            #print(np.arange(len(self.idx)))

            self.NNcidx_not = np.delete(np.arange(len(self.idx)),self.NNcidx)
            #print('self.state[svac]',self.state[svac])

            # this is just an "A" or an "S" or an "M"
            self.NN1 = self.state[self.NN1idx]
            self.NN2 = self.state[self.NN2idx]
            self.NNc = self.state[self.NNcidx]

            if False:
                self.NN1pos = self.sites[self.NN1idx]
                self.NN2pos = self.sites[self.NN2idx]
                self.NNcpos = self.sites[self.NNcidx]
            #print('self.NN1',self.NN1,type(self.NN1))
            #print('self.NN2',self.NN2)
            #print('self.NNc',self.NNc)
            #print('self.NN1 count S',gs(self.NN1,"S")) #(self.NN1=="S").sum())
            #print('self.NN1 count M',gs(self.NN1,"M")) #(self.NN1=="M").sum())
            #print('self.NN1 count A',gs(self.NN1,"A")) #(self.NN1=="A").sum())

            #print('self.NN1idx',self.NN1idx,type(self.NN1idx))
        return


    def vacancy_neighborhood_analysis_write(self,step,dt,cdf,nrand2):
        if step is not None:  # only than we want to write
            # this need to be separated from step==0 for RESTART's

            if step == 0:
                if os.path.isfile("KMC_vacancyneigh"): os.remove("KMC_vacancyneigh")
                if not hasattr(self, 'file_KMC_vacancyneigh'):
                    self.file_KMC_vacancyneigh = open("KMC_vacancyneigh","a+")
                self.file_KMC_vacancyneigh.write("# column   1     --> step\n")
                self.file_KMC_vacancyneigh.write("# column   2     --> time{picosecond} : The elapsed simulation time.\n")
                self.file_KMC_vacancyneigh.write("# column   3     --> time{picosecond} : The time of the current timestep.\n")
                self.file_KMC_vacancyneigh.write("# column   4     --> time{picosecond} : The time of the timestep (using fixed random number).\n")
                self.file_KMC_vacancyneigh.write("# column   5     --> ecurr (meV)\n")
                self.file_KMC_vacancyneigh.write("# column   6     --> cdf (sum)\n")
                self.file_KMC_vacancyneigh.write("# column   7     --> random number time\n")
                for ivac in range(self.nvac):
                    #print('ivac',ivac)
                    self.file_KMC_vacancyneigh.write("# column   8     --> vac:"+str(ivac+1)+" 1st NN Al\n")
                    self.file_KMC_vacancyneigh.write("# column   9     --> vac:"+str(ivac+1)+" 1st NN Mg\n")
                    self.file_KMC_vacancyneigh.write("# column   10    --> vac:"+str(ivac+1)+" 1st NN Si\n")
                    self.file_KMC_vacancyneigh.write("# column   11    --> vac:"+str(ivac+1)+" 2nd NN Al\n")
                    self.file_KMC_vacancyneigh.write("# column   12    --> vac:"+str(ivac+1)+" 2nd NN Mg\n")
                    self.file_KMC_vacancyneigh.write("# column   13    --> vac:"+str(ivac+1)+" 2nd NN Si\n")

            dtr = -1.0/cdf*np.log(1.0 - 0.5)
            self.file_KMC_vacancyneigh.write("%12.0f  %12.4e  %12.4e  %12.4e  %18.11e %12.4e  %8.5f "% (step, self.tottime*2.418884326509e-17, dt*2.418884326509e-17, dtr*2.418884326509e-17, self.ecurr, cdf, nrand2))
            #,self.NN1.count("A"),self.NN1.count("M"),self.NN1.count("S"),self.NN2.count("A"),self.NN2.count("M"),self.NN2.count("S")) )

                #self.file_KMC_vacancyneigh.write(" %4.0f %4.0f %4.0f %4.0f %4.0f %4.0f "% (self.NN1.count("A"),self.NN1.count("M"),self.NN1.count("S"),self.NN2.count("A"),self.NN2.count("M"),self.NN2.count("S")) )
            for ivac in xrange(self.natoms, self.natoms + self.nvac):
                self.file_KMC_vacancyneigh.write(" %4.0f %4.0f %4.0f %4.0f %4.0f %4.0f "% (gs(self.NN1,"A"),gs(self.NN1,"M"),gs(self.NN1,"S"),gs(self.NN2,"A"),gs(self.NN2,"M"),gs(self.NN2,"S")) )
            self.file_KMC_vacancyneigh.write("\n")
            self.file_KMC_vacancyneigh.flush()
        return

    def KMC_cluster_analysis(self,step=None,verbose=False):
        ''' determine size of clusters and mg/si content
            KMC_cluster_analysis(nsi=self.nsi,nmg=self.nmg,nvac=self.nvac,idx=self.idx,state=self.state,ridx=self.ridx,neigh=self.neigh)
        '''
        #all_solute_idx = np.zeros(nsi+self.nmg)
        #print('>> step',step) #,"step 1813 makes problems")
        allpairs = np.empty(((self.nsi+self.nmg+self.nvac)*12,2))
        allpairs[:] = np.nan
        running_idx = 0
        for solute in range(self.nsi + self.nmg): # is already the index
            sidx = self.idx[solute] # lattice site associated with this solute
            #all_solute_idx[solute] = sidx
            for sneigh in self.neigh[sidx]:
                running_idx += 1
                #print('running_idx',running_idx,self.state[sneigh])
                if self.state[sneigh] != 'A':
                    #print(">>++",self.state[sidx],'solute:',solute,'(out of',self.nsi+self.nmg,') now sidx',sidx,'neighbor:',self.state[sneigh],'sneigh',sneigh,'ridx',self.ridx[sneigh],'running_idx',running_idx)
                    allpairs[running_idx] = [sidx,sneigh]

        allpairs = allpairs[~np.isnan(allpairs[:,0])].astype(int)  # get rid of nans
        #print('allpairs rid of nans')
        #print(allpairs)

        # sort (and make unique) allpairs in case there are several.
        if len(allpairs) > 0:
            allpairs.sort(axis=1)  # sort so it will be possible to get the uniqe pairs
            #print('allpairs sorted')
            #print(allpairs)
            allpairs = np.unique(allpairs,axis=0) # get unique pairs
            #print('allpairs unique')
            #print(allpairs)
            #print("-- all unique pairs (1) --")
            #print(allpairs)
        clusters = get_cliques(allpairs)
        #clusters_max = self.nmg+self.nsi+self.nvac)/2
        if False:
            print('clusters',clusters)
            print("how many cluster?:",len(clusters))
            #print("theoretical clusters max:",clusters_max)
        cluster_sizes = np.zeros((self.nmg+self.nsi+self.nvac)/2).astype(int)
        cluster_mg = np.zeros((self.nmg+self.nsi+self.nvac)/2).astype(int)
        cluster_si = np.zeros((self.nmg+self.nsi+self.nvac)/2).astype(int)
        cluster_vac = np.zeros((self.nmg+self.nsi+self.nvac)/2).astype(int)
        for cluster_idx,cluster in enumerate(clusters):
            cluster_sizes[cluster_idx] = len(cluster)
            #print('cluster:',cluster)
            for index in cluster:
                #print('index',index,'type',self.state[index])
                if self.state[index] == "M": cluster_mg[cluster_idx] += 1
                if self.state[index] == "S": cluster_si[cluster_idx] += 1
                if self.state[index] == "V": cluster_vac[cluster_idx] += 1
                #solute_type = self.state(index)
                #print("i",solute_type)
        self.cluster_sizes = cluster_sizes
        self.cluster_mg    = cluster_mg
        self.cluster_si    = cluster_si
        self.cluster_vac   = cluster_vac
        if verbose:
            print(step,"cluster_sizes",self.cluster_sizes)
            print(step,"cluster_mg   ",self.cluster_mg)
            print(step,"cluster_si   ",self.cluster_si)
            print(step,"cluster_vac  ",self.cluster_vac)

        # save the stuff to file
        if step is not None:  # only than we want to write
            if step == 0:
                f = "KMC_cluster_"
                for i in [ "sizes", "mg", "si", "vac" ]:
                    if os.path.isfile(f+i): os.remove(f+i)
                    #print('removed',f+i)


            # this need to be separated from step==0 for RESTART's
            if not hasattr(self, 'file_KMC_cluster_sizes'):
                self.file_KMC_cluster_sizes = open("KMC_cluster_sizes","a+")
                self.file_KMC_cluster_mg = open("KMC_cluster_mg","a+")
                self.file_KMC_cluster_si = open("KMC_cluster_si","a+")
                self.file_KMC_cluster_vac = open("KMC_cluster_vac","a+")


            self.file_KMC_cluster_sizes.write("%12.0f "% (step))
            for i in range(len(cluster_sizes)):
                self.file_KMC_cluster_sizes.write("%5.0f "% (cluster_sizes[i]))
            self.file_KMC_cluster_sizes.write("\n")
            self.file_KMC_cluster_sizes.flush()

            self.file_KMC_cluster_mg.write("%12.0f "% (step))
            for i in range(len(cluster_mg)):
                self.file_KMC_cluster_mg.write("%5.0f "% (cluster_mg[i]))
            self.file_KMC_cluster_mg.write("\n")
            self.file_KMC_cluster_mg.flush()

            self.file_KMC_cluster_si.write("%12.0f "% (step))
            for i in range(len(cluster_si)):
                self.file_KMC_cluster_si.write("%5.0f "% (cluster_si[i]))
            self.file_KMC_cluster_si.write("\n")
            self.file_KMC_cluster_si.flush()

            self.file_KMC_cluster_vac.write("%12.0f "% (step))
            for i in range(len(cluster_vac)):
                self.file_KMC_cluster_vac.write("%5.0f "% (cluster_vac[i]))
            self.file_KMC_cluster_vac.write("\n")
            self.file_KMC_cluster_vac.flush()
            #print(step,"saving")
        return

    def try_load_ECACHE_QCACHE_files(self):
        print(">> (7) loading cache ...")
        if os.path.isfile(self.ecache_file):
            print('>>     found',self.ecache_file,'file')
        else:
            print('>>     did not find',self.ecache_file)
        if os.path.isfile(self.qcache_file):
            print('>>     found',self.qcache_file,'file')
        else:
            print('>>     did not find',self.qcache_file)

        if not os.path.isfile(self.qcache_file) or not os.path.isfile(self.ecache_file):
           print(">>      Couldn't load cache files ",self.qcache_file,self.ecache_file," - resetting")
           self.ecache = {}
           self.qcache = {}
           return

        if os.path.isfile(self.qcache_file) and os.path.isfile(self.ecache_file):
            # first get ecache file
            try:
                print()
                print('>>     pickle.load(ecache_file) [ca 10 sec for 200MB file]',self.ecache_file)
                self.ecache = load_cache_file(ecache_file)
                #start_time = time.time()
                #ff = open(self.ecache_file, "rb")
                #self.ecache = pickle.load(ff)
                #ff.close()
                print('>>     pickle.load(ecache_file) done.')
                print(">>     in --- %s seconds ---" % (time.time() - start_time))
                print(">>     Loaded %d cached energies" % (len(self.ecache)))
                print()
            except:
                print(">>     TRYING TO REPAIR DAMAGED ECACHE FILE")
                self.repair_cache_file(file=self.ecache_file)
                try:
                    start_time = time.time()
                    ff = open(self.ecache_file, "rb")
                    self.ecache = pickle.load(ff)
                    ff.close()
                    print('>>     pickle.load(ecache_file) done.')
                    print(">>     in --- %s seconds ---" % (time.time() - start_time))
                    print(">>     Loaded %d cached energies" % (len(self.ecache)))
                except:
                    print(">>      Couldn't load cache file "+self.ecache_file+" - resetting")
                    self.ecache = {}
                    self.qcache = {}
                    return
            key0 = self.ecache.keys()[0]
            print('key0',key0,len(key0))
            print('self.nsites',self.nsites)
            print('self.natoms',self.natoms)

            try:
                print()
                print('>>     pickle.load(qcache_file) [ca 80 sec for 9GB file]',self.qcache_file)
                start_time = time.time()
                ff = open(self.qcache_file, "rb")
                self.qcache = pickle.load(ff)
                ff.close()
                print('>>     pickle.load(qcache_file) done.')
                print(">>     in --- %s seconds ---" % (time.time() - start_time))
                print(">>     Loaded %d cached positions" % (len(self.qcache)))
                print()
            except:
                print(">>     TRYING TO REPAIR DAMAGED QCACHE FILE")
                self.repair_cache_file(file=self.qcache_file)
                try:
                    start_time = time.time()
                    ff = open(self.qcache_file, "rb")
                    self.qcache = pickle.load(ff)
                    ff.close()
                    print('>>     pickle.load(qcache_file) done.')
                    print(">>     in --- %s seconds ---" % (time.time() - start_time))
                    print(">>     Loaded %d cached positions" % (len(self.qcache)))
                    print()
                except:
                    print(">>      Couldn't load cache file "+self.qcache_file+" - resetting")
                    self.ecache = {}
                    self.qcache = {}
                    return


        if len(self.ecache.keys()) != len(self.qcache.keys()):
            print(">>     ECACHE and QCAHCE files of different lenght; removing unmached pairs ...")
            for key in self.qcache.keys():
                if key not in self.ecache.keys():
                    #print('deleting q',key)
                    del self.qcache[key]
            for key in self.ecache.keys():
                if key not in self.qcache.keys():
                    #print('deleting e',key)
                    del self.ecache[key]

        if len(self.ecache.keys()) != len(self.qcache.keys()):
            sys.exit(">>     ECACHE and QCAHCE files of different lenght; Exit")

        print(">>     Available %d cached energies." % (len(self.ecache)))
        print(">>     Available %d cached positions." % (len(self.qcache)))
        return

    def repair_cache_file(self,file=None):
        struct=0
        lastbs=0
        # read in the file
        line_list = [line.rstrip('\n') for line in open(file)]
        for idx,i in enumerate(line_list):
            #try:
            #print('i',i,i[:2])
            if i[:2] == 'bs':
                struct+=1
                lastbs=idx+1

        # repair the cache file and write to tmp file
        with open(file+"_tmp", 'wb') as text_file:
            for i in line_list[:lastbs-1]:
                text_file.write(i+"\n")
        with open(file+"_tmp", 'a') as text_file:
            text_file.write("bs.")

        os.rename(file+'_tmp',file)
        return

    def bind(self, ens, beads, nm, cell, bforce, prng):
        """Binds ensemble beads, cell, bforce, and prng to the dynamics.

        This takes a beads object, a cell object, a forcefield object and a
        random number generator object and makes them members of the ensemble.
        It also then creates the objects that will hold the data needed in the
        ensemble algorithms and the dependency network. Note that the conserved
        quantity is defined in the init, but as each ensemble has a different
        conserved quantity the dependencies are defined in bind.

        Args:
            beads: The beads object from whcih the bead positions are taken.
            nm: A normal modes object used to do the normal modes transformation.
            cell: The cell object from which the system box is taken.
            bforce: The forcefield object from which the force and virial are
                taken.
            prng: The random number generator object which controls random number
                generation.
        """


        self.prng = prng
        self.beads = beads
        self.cell = cell
        self.ens = ens
        self.forces = bforce

        # this is the index for the atoms, self.idx[i] indicates the lattice site of atom i.
        # atoms are in the order Si1 Si2 ... Mg1 Mg2 .... Al1 Al2 .... Vac1 Vac2 ...
        f_restart = True
        if self.idx is None or len(self.idx)==0:
            f_restart = False
            idx =np.asarray(range(self.ncell**3), int) # initialize random distribution of atoms
            if self.conventional:
                idx =np.asarray(range(4*self.ncell**3), int) # initialize random distribution of atoms
            self.prng.rng.shuffle(idx)
            self.idx = idx


        # initialize state based on the index
        # this is a string indicating the state of the system. each character corresponds to a lattice site
        # and can be S=Si, M=Mg, A=Al, V=Vac.
        # create a random string
        names = np.asarray(list("S" * self.nsi + "M" * self.nmg + (self.natoms - self.nsi - self.nmg) * "A" +  "V" * self.nvac))
        state = names.copy()

        # this maps the string to random sites
        state[self.idx] = names  #backshuffle!
        state = "".join(state)

        # reverse lookup index [i.e. ridx[i] gives the index of the atoms at site i]
        self.ridx = np.zeros(self.nsites,int)
        self.ridx[self.idx] = range(self.nsites)

        self.state = np.asarray(list(state))
        print ("".join(self.state))


        print("CHECKING INITIAL ASSIGNMENTS")
        for i in xrange(self.nsites):
            if self.ridx[i]<self.natoms:
                if self.state[i] in ['M','S']:
                    print(i,self.beads.names[self.ridx[i]], self.state[i])
            else:
                print(i,"V", self.state[i])
            if self.idx[self.ridx[i]] != i:
                print("inconsistent site string for atom ",i, " and site ", self.ridx[i])

        if not f_restart:
            self.beads.q[0,:] = self.sites[self.idx].flatten() # also sets global q so we can use it further down

        if True:
            for ivac in xrange(self.natoms, self.natoms + self.nvac):
                svac = self.idx[ivac] # lattice site associated with this vacancy
                print('ivac',ivac,'svac',svac)
                NN1 = self.neigh[svac]
                print('NN1 (indexes for NN1)',NN1)
                self.NN1 = self.state[NN1]
                print('self.NN1 = self.state[NN1]',self.NN1)
                self.NN1idx = self.idx[NN1]
                print('self.NN1idx',self.NN1idx,type(self.NN1idx))
        #sys.exit()

        self.dbeads = [None] * self.neval
        self.dforces = [None] * self.neval
        self.dnm = [None] * self.neval
        self.dens = [None] * self.neval
        self.dbias = [None] * self.neval
        for i in xrange(self.neval):
            self.dbeads[i] = beads.copy()
            self.dforces[i] = bforce.copy(self.dbeads[i], self.dcell)
            self.dnm[i] = nm.copy()
            self.dens[i] = ens.copy()
            self.dnm[i].bind(self.dens[i], self, beads=self.dbeads[i], forces=self.dforces[i])
            self.dbias[i] = ens.bias.copy(self.dbeads[i], self.dcell)
            self.dens[i].bind(self.dbeads[i], self.dnm[i], self.dcell, self.dforces[i], self.dbias[i])
            self.geop[i].bind(self.dens[i], self.dbeads[i], self.dnm[i], self.dcell, self.dforces[i], prng)

        #self.ebeads = self.dbeads.copy()
        #self.eforces = self.dforces.copy()
        #self.enm = self.dnm.copy()
        #self.eens= self.dens.copy()
        #self.ebias = self.dbias.copy()

        self.feval = np.ones(self.neval,int)
        self._threadlock = threading.Lock()
        self.dbeads_names = self.dbeads[0].names.copy()
        print('sdn',self.dbeads[0].names)
        print('sdn',self.dbeads_names)
        #sys.exit()

    # threaded geometry optimization
    def geop_thread(self, ieval, nstr, nevent, ostr=None,savename=None):
        if type(savename) == str and self.write_pos == True:
            # in positions
            p3 = np.reshape(self.dbeads[ieval].q[0,:],(-1,3))
            atomsxx = self.atomsc.copy()
            atomsxx.set_positions(p3*0.52917721)
            ase_write("in_new_"+savename+".extxyz",atomsxx,format="extxyz",append=False)
        #if nevent[0] == 27 and nevent[1] == 11:
        #    print('this this is isel 1 in step 0')
        #    print('the starting positions are:')
        #    print('self.dbeads[ieval].q[0,:][7:]')
        #    print('now dbeads',ieval,self.dbeads[ieval].q[0,:])
        if self.ase == True:
            p3 = np.reshape(self.dbeads[ieval].q[0,:],(-1,3))
            atomsc = self.atomsc.copy()
            atomsc.set_positions(p3*0.52917721)
            atomsc.set_chemical_symbols(self.dbeads[ieval].names)
            #if nevent[0] == 27 and nevent[1] == 8:
                #print('s{vac,neigh}',nevent[0],nevent[1],'idx',idx,atomsc.get_chemical_symbols()[idx],atomsc.positions[idx])
            #if type(savename) == str and self.write_pos == True:
            #    ase_write("in_"+savename+".extxyz",atomsc,format="extxyz",append=False)

            newpot,newq,atomsc = ase_minimize(atomsc_in=atomsc,minimizer="gpmin") #lbfgs")
            #if type(savename) == str and self.write_pos == True:
            #    ase_write("out_"+savename+".extxyz",atomsc,format="extxyz",append=False)
            self.dbeads[ieval].q[0] = newq

        else:
            self.geop[ieval].reset()
            ipot = self.dforces[ieval].pot
            #print('geop ipot',ipot)

            #print("geop step", 'A','ene', self.dforces[ieval].pot)
            for i in xrange(self.nstep):
                #print("geop step", i,'ieval',ieval,'nstr',nstr) #,ostr,'forces', self.dforces[ieval].pot)
                #print("geop step", i,'ieval',ieval,'ostr',ostr)
                #print("geop step", i,'ene', self.dforces[ieval].pot)
                self.geop[ieval].step(i)
                #if self.geop[ieval].converged[0]: break
            #if nevent[0] == 27 and nevent[1] == 11:
            #    print('IPI final positions == newq',self.dforces[ieval].pot)
            #    print('IPI dbeads',ieval,self.dbeads[ieval].q[0,:])
            newq = dstrip(self.dbeads[ieval].q[0]).copy()  # final positoins
            newpot =  self.dforces[ieval].pot

        # print "geop ", self.nstep, self.dforces[ieval].pot
        with self._threadlock:
            #print('save gope nstr',nstr)
            self.ecache[nstr] = newpot
            self.qcache[nstr] = newq
            self.ncache += 1
            nevent[2] = newpot # self.ecache[nstr]
            nevent[3] = newq #self.qcache[nstr]


            if type(savename) == str and self.write_pos == True:
                # out positions
                p2 = np.reshape(newq,(-1,3))
                atomsxx.set_positions(p2*0.52917721)
                ase_write("out_new_"+savename+".extxyz",atomsxx,format="extxyz",append=False)

        # launches TS calculation
        #if not ostr is None:
        #    self.ts_thread(ieval, ostr, nstr, nevent)
        #else:
        #    with self._threadlock:
        #        self.feval[ieval] = 1
        with self._threadlock:
            self.feval[ieval] = 1

        with self._threadlock:
            if False:
                print("Finished ", nstr)
                print("Energy, initial - TS - final: ", ipot, nevent[-1], newpot)

    # threaded ts evaluation
    def ts_thread(self, ieval, ostr, nstr, nevent, setfev=1):
        # computes TS energy by linearly interpolating initial & final state
        # interpolates between two structures considering PBCs
        qstart = self.qcache[ostr]
        qend = self.qcache[nstr].copy()
        # finds atom matching assuming initial and final states differ only by a vacancy swap and are based on unique_idx lists
        midx = self.unique_idx_match(list(ostr), list(nstr))[0:self.natoms]
        qend.shape = (self.natoms, 3)
        qend[:] = qend[midx]
        qend.shape = (3*self.natoms)

        qts = qend - qstart
        self.dcell.array_pbc(qts) # use pbc!
        qts *= 0.5
        qts += qstart
        self.dbeads[ieval].q[0,:] = qts

        tspot = self.dforces[ieval].pot
        io.print_file("xyz", self.beads[0], self.dcell, self.tslist, title=("START  "  ), key="positions", dimension="length", units="angstrom", cell_units="angstrom" )
        io.print_file("xyz", self.dbeads[ieval][0], self.dcell, self.tslist, title=("TS: %s %s  Energy: %15.8e  %15.8e " % (ostr, nstr, tspot, tspot-self.forces.pot) ), key="positions", dimension="length", units="angstrom", cell_units="angstrom" )
        self.dbeads[ieval].q[0] = qend
        io.print_file("xyz", self.dbeads[ieval][0], self.dcell, self.tslist, title=("END  "  ), key="positions", dimension="length", units="angstrom", cell_units="angstrom" )
        self.tslist.flush()

        with self._threadlock:
            # sets the tscache for both FW and BW transition (uses a dictionary to be super-safe and lazy, although of course this could be just a list)
            self.tscache[ostr][nstr] = tspot
            if not nstr in self.tscache: self.tscache[nstr]  = {}
            self.tscache[nstr][ostr] = tspot
            self.feval[ieval] = 1
            nevent[4] = tspot

    def find_eval(self, ethreads):
        # finds first free evaluator
        while self.feval.sum() == 0: # if all evaluators are busy, wait for one to get free
            for st in ethreads:
                st.join(1e-2)
                if st is None or not st.isAlive():
                    break
        with self._threadlock:
            # finds free evaluator
            for e in xrange(self.neval):
                if self.feval[e] == 1:
                    ieval = e
                    self.feval[ieval] = 0
                    break
        return ieval

    def unique_idx(self, state):
        # generates a starting lattice configuration that corresponds to a given state vector (a string of atomic types)
        # makes sure that the same occupation string corresponds to the same atom positions,
        # even though the actual atom indices might be different basically, makes the configurations
        # independent of same-atom permutations
        ksi = 0
        kmg = self.nsi
        kal = self.nsi+self.nmg
        kvac = self.natoms
        k = 0
        idx = np.zeros(self.nsites, int)
        for s in state:
            if s == "S":
                idx[ksi] = k
                ksi += 1
            elif s == "M":
                idx[kmg] = k
                kmg += 1
            elif s == "A":
                idx[kal] = k
                kal += 1
            elif s == "V":
                idx[kvac] = k
                kvac += 1
            k += 1
        return idx

    def unique_idx_match(self, state_1, state_2):
        # finds the best matching between the atoms in states state_1 and state_2,
        # assuming there is only one site swap
        # (should raise an error otherwise but it's not trivial without doing too many extra checks)
        # this basically says what is the motion of atoms from state 1 to state 2. This is useful
        # to interpolate between atomic coordinates in the two states, e.g. to get the TS geometry

        # gets the unique mapping of atoms from the state to the site ids
        uid_1 = self.unique_idx(state_1)
        ru1 = np.zeros(self.nsites,int)
        # this is the reverse map. what is the atom index that sits in a given site?
        ru1[uid_1] = np.asarray(range(self.nsites),int)

        uid_2 = self.unique_idx(state_2)
        ru2 = np.zeros(self.nsites,int)
        ru2[uid_2] = np.asarray(range(self.nsites),int)  # this says which atom is in a given site in u2

        iu12 = ru2[uid_1]
        iu21 = ru1[uid_2]
        for i in xrange(self.natoms):
            if iu12[i] >= self.natoms:
                #print "found vacancy swap 1->2", i, u1[i], u2[iu12[i]], iu12[i]
                i1vac2 = i
            if iu21[i] >= self.natoms:
                i2vac1 = i
                #print "found vacancy swap 2->1", i, u2[i], u1[iu21[i]], iu21[i]
        iu12[i1vac2] = i2vac1

        return iu12

    def get_list_of_rates(self,step,levents,kT,filled="UNKNOWN"):
        # get list of rates
        if type(self.step_isel) == bool:
            isel_loc = -1
        else:
            isel_loc = self.step_isel[step]
        rates = np.zeros(len(levents), float)
        rates_new = np.zeros(len(levents), float)
        crates = np.zeros(len(levents), float)
        cdf = 0.0
        print()
        print('------------------------------------------------------------------------------------------------')
        print('------------------------------------------------------------------------------------------------')
        #rint("Event       E_curr (self.ecurr) >> highest (ets)       >>  E_fin = levents[i][2] || E_fin-E_curr")
        print("Event       E_curr (self.ecurr) >>                     >>  E_fin = levents[i][2] || E_fin-E_curr")
        print('------------------------------------------------------------------------------------------------')

        def kaa(text,l=14):
            return str(text).ljust(l)
        def kab(text,l=14):
            a = kaa(text,l=l)
            if a[0] == "-":
                return a
            else:
                return " "+a

        def kac(text,l1=11,l2=8):
            ka = "{:"+str(l1)+"."+str(l2)+"f}"
            #print('ka',ka)
            return (ka.format(text)).rjust(l1+1)

        hartree_to_mev = 27211.386
        #diff_mev = np.zeros((12,2))

        for i in xrange(len(levents)):
            svac    = levents[i][0]
            sneigh  = levents[i][1]
            ene_fin = levents[i][2]   # == same as self.ecache[nstr]
            #pos = levents[i][3]
            ostr    = levents[i][4]
            nstr    = levents[i][5]
            ostr_solutes = ostr.count('M')+ostr.count('S')
            nstr_solutes = nstr.count('M')+nstr.count('S')
            if ostr_solutes != nstr_solutes:
                softexit.trigger("Error: ostr_solutes != nstr_solutes. Exit")
            print('ostr_solutes',ostr_solutes)
            print('nstr_solutes',nstr_solutes)
            ene_ref = self.ecache[ostr]
            ene_fin2 = self.ecache[nstr]
            #print('ostr',ostr,'ene_ref',ene_ref)
            #print('nstr',nstr,'ene_fin',ene_fin,ene_fin2)
            ene_base = self.ecurr
            #print('ene_base',ene_base)
            diff_to_base = ene_ref - ene_base
            diff_to_base_meV = diff_to_base*hartree_to_mev
            #print('ene_base',ene_base,'ene_ref',ene_ref,'diff_to_base',diff_to_base,'mev',diff_to_base_meV)



            #print ("Barrier, naive: %f, static: %f" % (0.5*(ecurr + levents[i][2]) + self.diffusion_barrier_al, levents[i][4]))
            A_B = (str(svac)+"/"+str(sneigh)).ljust(6)  # 53/42   = svac/sneigh
            pp = "++"
            if filled == "normal":
                pp = "+N"
            if filled == "filled":
                pp = "+F"
            if i == isel_loc:
                pp = "**"
            #diffh = levents[i][2]-self.ecurr
            ene_diff = ene_fin - ene_ref # this also works for the filled structures
            # if ene_diff is messed up: ostr or nstr are wrongly saved.
            ene_diff_mev = (ene_diff)*hartree_to_mev #/self.neigh.shape[0]

            ets  = 0.5*(ene_ref + levents[i][2]) + self.barriers[levents[i][-1]]  #naive heuristic for the ts energy
            etsnew = 0.5*ene_diff+ene_ref        + self.barriers[levents[i][-1]]

            #ekold = ets - ene_ref
            #eknew = 0.5*diffh
            #print('ets==etsnew?',ets==etsnew,ets,etsnew,ene_diff_mev,"mev/def")
            #print('eko==ekn   ?',ekold==eknew,ekold,eknew)
            if True:
                etsp  = 0.5*(ene_ref + ene_fin) #*hartree_to_mev
                print(pp,step,str(i).ljust(3),
                        kaa(levents[i][-1],3),
                        A_B,
                        #kaa(ene_ref+diff_to_base),
                        kac(ene_ref-diff_to_base),
                        ">>",
                        #("{:10.8f}".format(ene_fin)).rjust(11),
                        kac(ene_fin-diff_to_base),
                        "||",
                        kac(ene_diff),
                        '== meV/def:',
                        kac(ene_diff_mev,7,2)) #,etsp)
            if np.abs(ene_diff_mev) > 10000:
                softexit.trigger("Error: ene diff too high. Exit")
                #print("Error: ene diff too high. Exit")
            rates[i] = self.prefactors[levents[i][-1]] * np.exp(-(ets-ene_ref)/kT)
            rates_new[i] = self.prefactors[levents[i][-1]] * np.exp(-(etsnew-ene_ref)/kT)


            cdf += rates[i]
            crates[i] = cdf
            print('ets',ets,'rates[i]    :',rates[i])
            print('ets',ets,'rates[i]_new:',rates_new[i])

            #diff_mev[i][0] = step*12+i
            #diff_mev[i][0] = step+(i/len(levents))
            #diff_mev[i][1] = ene_diff_mev

            x = step+((i*1.0)/len(levents))
            #print('i',i,'x',x)
            y = ene_diff_mev


            if filled == "filled":
                self.diff_mev_filled.append([x,y])
                ya = np.abs(y-self.diff_mev_normal_thisstep[i][1])
                self.diff_diff_abs_mev_filled.append([x,ya])

                solutes_diff = (self.nsi+self.nmg) - ostr_solutes
                self.ediff_vs_solutes.append([solutes_diff,ya])
            elif filled == "normal":
                self.diff_mev_normal.append([x,y])
                self.diff_mev_normal_thisstep.append([x,y])

        if filled == "filled":
            np.savetxt("nnrates_diff_filled_"+str(self.cutoff_filled),self.diff_mev_filled)
            np.savetxt("nnrates_diff_filledabs_"+str(self.cutoff_filled),self.diff_diff_abs_mev_filled)
            np.savetxt("nnrates_solutes_vs_energy_diff_"+str(self.cutoff_filled),self.ediff_vs_solutes)
            kk = self.diff_mev_normal
            kk = np.array([self.diff_mev_normal,self.diff_mev_filled])
            #print('kk')
            #print(kk)
            kkx = kk[0][:,1]
            kky = kk[1][:,1]
            kk = np.array([kkx,kky]).T
            #print('kk out')
            #print(kk)

            np.savetxt("nnrates__diff_"+str(self.cutoff_filled),kk)
            #if step == 1:
            #    softexit.trigger("Error: kk. Exit.")
        elif filled == "normal":
            np.savetxt("nnrates_diff_normal",self.diff_mev_normal)

        #np.savetxt("rates_ORIG_step_"+str(step)+"_"+filled,rates)
        #np.savetxt("rates_NEW_step_"+str(step)+"_"+filled,rates_new)
        if False:
            out = np.zeros(len(levents)+1)
            out[0] = self.ecurr
            #print(kaa(self.ecurr))
            for i in xrange(len(levents)):
                out[i+1] = levents[i][2]
            if os.path.isfile("step_"+str(step)+".dat"):
                os.remove("step_"+str(step)+".dat")
            np.savetxt("step_"+str(step)+".dat",out)
        print('------------------------------------------------------------------------------------------------')
        print('------------------------------------------------------------------------------------------------')
        print()
        return rates,crates,cdf



    def step(self, step=None):
        print()
        print(">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>")
        print('>> step',step)
        if False:
            self.atomsc.set_positions(self.sites[self.idx]*0.52917721)
            ase_write("simulation.pos_0.sitesidx.extxyz",self.atomsc,format="extxyz",append=True)

            self.atomsc.set_positions(self.sites[self.unique_idx(self.state.copy())]*0.52917721)
            ase_write("simulation.pos_0.sitesidxunique.extxyz",self.atomsc,format="extxyz",append=True)


        #### neighborhood sphere
        print(">> step",step,'vacancy_neighborhood_analysis()')
        self.vacancy_neighborhood_analysis()

        if False:
            self.atomsc_neigh.set_positions(self.sites[self.NN1idx]*0.52917721)
            self.atomsc_neigh.set_chemical_symbols(str_to_symbols(self.NN1))
            ase_write("simulation.pos_0.sphere1.extxyz",self.atomsc_neigh,format="extxyz",append=True)

            self.atomsc_neigh2.set_positions(self.sites[self.NN2idx]*0.52917721)
            self.atomsc_neigh2.set_chemical_symbols(str_to_symbols(self.NN2))
            ase_write("simulation.pos_0.sphere2.extxyz",self.atomsc_neigh2,format="extxyz",append=True)


            self.atomsc_neighc.set_positions(self.sites[self.NNcidx]*0.52917721)
            self.atomsc_neighc.set_chemical_symbols(str_to_symbols(self.NNc))
            ase_write("simulation.pos_0.spherec.extxyz",self.atomsc_neighc,format="extxyz",append=True)

        if False:
            self.atomsc_filled.set_chemical_symbols("Al")
            p3 = np.reshape(self.beads.q[0,:],(-1,3))*0.52917721; 	#-> all atoms are relaxed, not mapped into original cell.
            self.atomsc_filled.set_positions(p3)


        if False:
            p2 = np.reshape(self.dbeads[0].q[0,:],(-1,3))
            if True:
                self.atomsc.set_positions(p2*0.52917721)
                ase_write("simulation.pos_0.dbeads.extxyz",self.atomsc,format="extxyz",append=True)
                ene,positions,atomsc = ase_minimize(atomsc_in=self.atomsc,minimizer="gpmin") #lbfgs")


        #print('p3',p3.shape)
        #np.savetxt("pos_beads_"+str(step)  ,p3*0.52917721)
        if True:
            p3 = np.reshape(self.beads.q[0,:],(-1,3))
            self.atomsc.set_positions(p3*0.52917721)
            ase_write("simulation.pos_0.beads.extxyz",self.atomsc,format="extxyz",append=True)

        if False:
            p4 = np.reshape(dstrip(self.beads.q[0,:]),(-1,3))
            #print('p4',p4.shape)
            #np.savetxt("pos_beadsds_"+str(step)  ,p4*0.52917721)
            if True:
                self.atomsc.set_positions(p4*0.52917721)
                ase_write("simulation.pos_0.beadsds.extxyz",self.atomsc,format="extxyz",append=True)

        print(">> step",step,'self.KMC_cluster_analysis(step=step)')
        self.KMC_cluster_analysis(step=step)

        kT = Constants.kb  * self.ens.temp

        ##############################################################
        # computes current energy, self.ecurr (if not already stored)
        ##############################################################
        ostr = "".join(self.state)  # this is a unique id string that charactrizes the current state
        self.tscache[ostr] = {}
        print('>> step',step,"calc ostr")
        if not ostr in self.ecache:
            ieval = 0
            #self.get_dbeads_names(ieval,self.state,filled = False,nstate_filled=False)
            self.dbeads[ieval].q[0,:] = self.sites[self.unique_idx(self.state)].flatten()
            self.dbeads[ieval].names = self.dbeads_names
            print('self.dbeads[0].q[0,:] in (this is the structure of the very first step 0)')
            print(self.dbeads[0].q[0,:][:7])
            rv = [0, 0, 0, 0, 0]
            savename = "pos_step_"+str(step)+"_normal"
            self.geop_thread(ieval, ostr, rv, None,savename)
        self.ecurr = self.ecache[ostr]
        print('>> step',step,"calc ostr done")

        ##############################################################
        # enumerates possible reactive events (vacancy swaps)
        ##############################################################
        if self.normal:
            starttime = time.time()
            levents        = self.enumerate_possible_reaction_events(ostr,step=step,filled=False)
            endtime = time.time()
            self.diff_mev_normal_thisstep = []
            self.solutes_thisstep = []
            rates,crates,cdf = self.get_list_of_rates(step,levents,kT,filled="normal")
            print('++++++++++++++++++++++++++++++++++++++++++++++++++++++++++',endtime-starttime)
            #print('serial   6 steps: 11:08:02 to 11:13:20 == 5 min 18 sec')
            # 0   1   2   3   4   5
            # 19, 11, 24, 37, 28, 32,  -> ser
            # 17, 11, 25, 31, 27, 31,  -> par, not much faster indeed!
            #print('parallel 6 steps: 11:14:57 to 11:18:04 == 3 min 7 sec')
            #print('par: 11:52:00 11:53:41 == 1 min 41 sec  # with cache
            #print('time serial   step 0: 22.1343359947')
            #print('time serial   step 0: 18.6')
            #print('time serial   step 0: 18.5')
            #print('time parallel step 0: 12.2087500095')
            #print('time parallel step 0: 13.9')
            #print('time parallel step 0: 17.2')

        if self.filled == True:
            print('calculaeting filled reference ..............')
            levents_filled = self.enumerate_possible_reaction_events(ostr,step=step,filled=False,filled_ref=True)
            print('calculaeted filled reference .............. done')
            print('++++++++++++++++++++++++++++++++++++++++++++++++++++++++++')
            print()
            print('calculaeting filled structures ..............')


            levents_filled = self.enumerate_possible_reaction_events(ostr,step=step,filled=True)
            #for i in levents_filled:
            #    print('old',i[4],'new',i[5])
            xrates,xcrates,xcdf = self.get_list_of_rates(step,levents_filled,kT,filled="filled")

        # KMC selection rule based on the rate
        fpick = self.prng.u*cdf
        isel = 0
        while fpick > crates[isel]:
            isel += 1
        nrand2 = self.prng.u
        dt = -1.0/cdf*np.log(1.0-nrand2)

        self.vacancy_neighborhood_analysis_write(step,dt,cdf,nrand2)

        #print ("Time spent %12.5e at %s nrg %12.5e" % (dt, ostr,ecurr))
        if type(self.step_isel) == bool:
            pass
        else:
            isel = self.step_isel[step]
        print(">> step:",step,"Selected event ", isel, " with ene",levents[isel][2]) #rate ", rates[isel], " / ", cdf)

        iev = levents[isel] # levents[self.prng.randint(len(levents))]
        svac, sneigh = iev[0], iev[1]
        ivac, ineigh = self.ridx[svac], self.ridx[sneigh]

        # does the swap (never reject, for the moment)
        self.state[svac], self.state[sneigh] = self.state[sneigh], self.state[svac] # [ 'A', 'S', ...]
        self.ridx[svac], self.ridx[sneigh] = self.ridx[sneigh], self.ridx[svac]
        self.idx[ivac], self.idx[ineigh] = self.idx[ineigh], self.idx[ivac]


        # we got a new configuration but the residence time is linked to the previous configuration so we output that
        self.kmcfile.write("%12.5e  %12.5e  %18.11e  %s\n"% (self.tottime, dt, self.ecurr, ostr) )
        self.kmcfile.flush()
        self.tottime += dt
        self.ens.time += dt  # updates time counter
        #print  "Finishing step at ", "".join(self.state)

        print('>> step:',step,"computed:",str(self.not_in_cache).ljust(3),"+",str(self.not_in_cache_filled).ljust(3),'from cache:',str(self.in_cache).ljust(3),"+",str(self.in_cache_filled).ljust(3),"clusters_sizes:",self.cluster_sizes,"events: ", len(levents), " possible reactions. Cache len:",len(self.ecache),'self.ecurr',self.ecurr,'isel',isel,'energy[isel]',levents[isel][2]) #,'ene ipi??',self.dforces[ieval].pot) #,"new cache len:",len(self.ecache_n))
        # updates the positions
        self.cell.h = self.dcell.h

        uidx = self.unique_idx(self.state)
        ruidx = np.zeros(self.nsites,int)
        ruidx[uidx] = range(self.nsites)

        self.sites[self.unique_idx(self.state)]
        oldq = dstrip(self.beads.q[0]).copy()

        #print('oldq (beads)')
        #print(oldq[:7])

        newq = np.zeros(self.nsites*3, float)
        # we want continuity (modulo PBC jumps, that we'll take care of later...)
        #print('iev[0] = svac, iev[1]=sneigh:',iev[0],iev[1])
        #print('iev[3] at the end')
        #print(iev[3][:7])
        #print()
        for i in xrange(self.nsites):
            # in which site sits atom i?
            #isite = self.idx[i]
            # which atom sits in this site in the unique-mapped structure?
            iuid = ruidx[self.idx[i]]
            #print('iuid',iuid.shape)
            #print('iev[2].shape',iev[2].shape)
            #print('iev[3]',iev[3])
            #if i == 0:
            #    print('iev[3]',iev[3])
            #    print('iev[3]',iev[3].shape)
            newq[3*i:3*i+3] = iev[3][3*iuid:3*(iuid+1)]
        ### UP TO HERE NO EXTRAPOLATION PROBLEMS
        if False:
            p6 = np.reshape(oldq,(-1,3))
            self.atomsc.set_positions(p6*0.52917721)
            ase_write("simulation.pos_0.oldq.extxyz",self.atomsc,format="extxyz",append=True)
            p5 = np.reshape(newq,(-1,3))
            self.atomsc.set_positions(p5*0.52917721)
            ase_write("simulation.pos_0.newq.extxyz",self.atomsc,format="extxyz",append=True)

        #print('newqbefore changeing')
        #print(newq[:7])
        newq-=oldq
        if False:
            p5 = np.reshape(newq,(-1,3))
            self.atomsc.set_positions(p5*0.52917721)
            ase_write("simulation.pos_0.newqq.extxyz",self.atomsc,format="extxyz",append=True)
        #print('newq after newq-=oldq (==difference will be added to beads)')
        #print(newq[:7])
        ### UP TO HERE NO EXTRAPOLATION PROBLEMS
        self.cell.array_pbc(newq)
        self.beads.q[0]+=newq
        #print('final beads')
        #print(self.beads.q[0][7:])
        ### HERE EXTRAPLATION PROBS`
        #print('1234aaa')
        #sys.exit('1234aaa')


