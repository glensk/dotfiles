"""Contains the classes that deal with the different dynamics required in
different types of ensembles.

Holds the algorithms required for normal mode propagators, and the objects to
do the constant temperature and pressure algorithms. Also calculates the
appropriate conserved energy quantity for the ensemble of choice.
"""

from __future__ import print_function
# This file is part of i-PI.
# i-PI Copyright (C) 2014-2015 i-PI developers
# See the "licenses" directory for full license information.

import sys,os
import time
import pickle
import threading
import numpy as np
from collections import defaultdict
#from ase import units as aseunits  # dont import ase, since it will not start ipi correctly on cluster

from ipi.engine.motion import Motion, GeopMotion
from ipi.utils.depend import dstrip, depend_value, dobject, dpipe, dd
from ipi.engine.thermostats import Thermostat
from ipi.engine.cell import Cell
from ipi.engine.normalmodes import NormalModes
from ipi.engine.barostats import Barostat
from ipi.utils.units import Constants
import ipi.utils.io as io

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

    print('----------- loading',ecache_file,"ECACHE/QCACHE ~10/80sec for 0.2GB/9GB file")
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
    print('----------- loading',ecache_file,'done. In',time.time() - start_time,'seconds, containing',len(x),'structures')
    #print('----------- loading ----- done',ecache_file)
    #print()
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
        self.scell=self.a0*cell
        self.dcell = Cell()
        self.dcell.h = self.scell*self.ncell

        print(">> (1) LATTICE PARAM ", self.a0,"(bohrradius) == ",self.a0*0.52917721,"(Angstrom)")
        # this is the list of lattice sites, in 3D coordinates
        ix,iy,iz = np.meshgrid(range(self.ncell), range(self.ncell), range(self.ncell), indexing='ij')
        self.sites = np.dot(np.asarray([ix.flatten(),iy.flatten(),iz.flatten()]).T, self.scell.T)
        print(">> (2) self.sites, self.nsites",len(self.sites), self.nsites, "###")
        # now we build list of nearest neighbors (fcc-lattice hardcoded!)
        self.neigh  = np.zeros((self.nsites,12),int)
        self.neigh2 = np.zeros((self.nsites,6),int)
        nneigh  = np.zeros(self.nsites,int)
        nneigh2 = np.zeros(self.nsites,int)
        # could be done in a more analytic way but whatever, I'm too lazy
        a02 = 1.01*0.5*self.a0**2                                        # perhaps 1.01 it is not enough, must check!
        a022 = 1.01*self.a0**2                                        # perhaps 1.01 it is not enough, must check!
        #print('a02',a02)

        print(">> (3) determines the connectivity of the lattice ...")
        #print("1NNs",self.neigh[0])
        #print(nneigh)
        for i in xrange(self.nsites): # determines the connectivity of the lattice
            #print('i',i)
            rij = self.sites.copy().flatten()
            for j in xrange(self.nsites):
                rij[3*j:3*j+3] -= self.sites[i]
            self.dcell.array_pbc(rij)
            rij.shape = (self.nsites,3)
            for j in xrange(i):
                if np.dot(rij[j],rij[j]) < a02: # found nearest neighbor
                    #print(i,'j',j,np.dot(rij[j],rij[j]))
                    self.neigh[i,nneigh[i]] = j
                    self.neigh[j,nneigh[j]] = i
                    nneigh[i]+=1
                    nneigh[j]+=1
        print(">> (4) 1NNs (of atom at idx=0):",self.neigh[0])
        #print(nneigh)

        print(">> (4) determines the connectivity of the lattice (for 2NN) ...")
        #print("2NNs",self.neigh2[0])
        #print(nneigh2)
        for i in xrange(self.nsites): # determines the connectivity of the lattice
            rij = self.sites.copy().flatten()
            for j in xrange(self.nsites):
                rij[3*j:3*j+3] -= self.sites[i]
            self.dcell.array_pbc(rij)
            rij.shape = (self.nsites,3)
            for j in xrange(i):
                dist = np.dot(rij[j],rij[j])
                if dist > a02 and dist < a022: # found second nearest neighbor
                    self.neigh2[i,nneigh2[i]] = j
                    self.neigh2[j,nneigh2[j]] = i
                    nneigh2[i]+=1
                    nneigh2[j]+=1
        print(">> (5) 2NNs (of atom at idx=0):",self.neigh2[0])
        #print(nneigh2)

        self.idx = idx

        # the KMC step is variable and so it cannot be stored as proper timing
        dd(self).dt = depend_value(name="dt", value = 0.0)
        self.fixatoms = np.asarray([])
        self.fixcom = True
        self.geop = [None] * self.neval
        # geop should not trigger exit if there is early convergence, but just carry on.
        # we hard-code this option to avoid early-termination that would be hard to debug for a user
        geop["exit_on_convergence"] = False
        print(">> (6) init geop ..")
        for i in xrange(self.neval):
            # geometry optimizer should not have *any* hystory dependence
            self.geop[i] = GeopMotion(fixcom=fixcom, fixatoms=fixatoms,**geop) #mode="cg", ls_options={"tolerance": 1, "iter": 20,  "step": 1e-3, "adaptive": 0.0}, tolerances={"energy": 1e-7, "force": 1e-2, "position": 1e-4}, ) #!TODO: set the geop parameters properly

        #sys.exit('test2')
        #print('still not out')
        #raise IndexError('test1')
        #print('still not out')

        # dictionary of previous energy evaluations.
        self.ecache_file = ecache_file
        self.qcache_file = qcache_file

        #self.try_load_ECACHE_QCACHE_files()
        self.ecache = load_cache_file(ecache_file)
        self.qcache = load_cache_file(qcache_file)
        self.ecache_n = {}
        self.qcache_n = {}

        self.ncache = len(self.ecache)
        self.ncache_stored = self.ncache


        # no TS evaluation implemented yet
        self.tscache = {}

        # todo make these optional and restarted
        print('>> (8) writing KMC_files ...')
        self.kmcfile = open("KMC_AL6XXX","w+")
        #self.kmcfile_energies = open("KMC_energies","w+")
        self.kmcfile_denergies= open("KMC_denergies","w+")
        #self.kmcfile_rates= open("KMC_rates","w+")
        #self.kmcfile_barriers = open("KMC_barriers","w+")
        self.tottime = tottime

    def barriers_energies_analysis(self,step,ecurr,levents,rates):
        if step is not None:  # only than we want to write
            if not hasattr(self, 'file_KMC_barriers'):
                self.file_KMC_barriers = open("KMC_barriers","a+")
            if not hasattr(self, 'file_KMC_rates'):
                self.file_KMC_rates = open("KMC_rates","a+")
            if not hasattr(self, 'file_KMC_energies'):
                self.file_KMC_energies = open("KMC_energies","a+")
            if not hasattr(self, 'file_KMC_denergies'):
                self.file_KMC_denergies = open("KMC_denergies","a+")

            self.file_KMC_energies.write("%18.11e"% (ecurr))
            for i in xrange(len(levents)):
                self.file_KMC_energies.write("%20.11e"% (levents[i][2]))
                self.file_KMC_barriers.write("%20.11e"% (self.barriers[levents[i][-1]]))
                self.file_KMC_rates.write("%20.11e"% rates[i])
                self.file_KMC_denergies.write("%20.11e"% (-0.5*(ecurr - levents[i][2])))

            for i in [self.file_KMC_energies,self.file_KMC_barriers,self.file_KMC_rates,self.file_KMC_denergies]:
                i.write("\n")
                i.flush()
        return

    def vacancy_neighborhood_analysis(self,step,dt,ecurr,cdf,nrand2):
        if step is not None:  # only than we want to write
            # this need to be separated from step==0 for RESTART's
            if not hasattr(self, 'file_KMC_vacancyneigh'):
                self.file_KMC_vacancyneigh = open("KMC_vacancyneigh","a+")

            if step == 0:
                self.file_KMC_vacancyneigh.write("# column   1     --> step\n")
                self.file_KMC_vacancyneigh.write("# column   2     --> time{picosecond} : The elapsed simulation time.\n")
                self.file_KMC_vacancyneigh.write("# column   3     --> time{picosecond} : The time of the current timestep.\n")
                self.file_KMC_vacancyneigh.write("# column   4     --> time{picosecond} : The time of the timestep (using fixed random number).\n")
                self.file_KMC_vacancyneigh.write("# column   5     --> ecurr (meV)\n")
                self.file_KMC_vacancyneigh.write("# column   6     --> cdf (sum)\n")
                self.file_KMC_vacancyneigh.write("# column   7     --> random number time\n")
                for ivac in range(self.nvac):
                    print('ivac',ivac)
                    self.file_KMC_vacancyneigh.write("# column   8     --> vac:"+str(ivac+1)+" 1st NN Al\n")
                    self.file_KMC_vacancyneigh.write("# column   9     --> vac:"+str(ivac+1)+" 1st NN Mg\n")
                    self.file_KMC_vacancyneigh.write("# column   10    --> vac:"+str(ivac+1)+" 1st NN Si\n")
                    self.file_KMC_vacancyneigh.write("# column   11    --> vac:"+str(ivac+1)+" 2nd NN Al\n")
                    self.file_KMC_vacancyneigh.write("# column   12    --> vac:"+str(ivac+1)+" 2nd NN Mg\n")
                    self.file_KMC_vacancyneigh.write("# column   13    --> vac:"+str(ivac+1)+" 2nd NN Si\n")

            dtr = -1.0/cdf*np.log(1.0 - 0.5)
            self.file_KMC_vacancyneigh.write("%12.0f  %12.4e  %12.4e  %12.4e  %18.11e %12.4e  %8.5f "% (step, self.tottime*2.418884326509e-17, dt*2.418884326509e-17, dtr*2.418884326509e-17, ecurr, cdf, nrand2))
            #,NN1.count("A"),NN1.count("M"),NN1.count("S"),NN2.count("A"),NN2.count("M"),NN2.count("S")) )

            # loops over the vacancy
            for ivac in xrange(self.natoms, self.natoms + self.nvac):
                NN1 = []
                NN2 = []
                svac = self.idx[ivac] # lattice site associated with this vacancy
                for sneigh in self.neigh2[svac]:
                    NN2.append(self.state[sneigh]) # this is just an "A" or an "S" or an "M"
                for sneigh in self.neigh[svac]:
                    NN1.append(self.state[sneigh]) # this is just an "A" or an "S" or an "M"


                #self.file_KMC_vacancyneigh.write("%12.0f  %12.4e  %12.4e  %12.4e  %18.11e %12.4e  %8.5f %4.0f %4.0f %4.0f %4.0f %4.0f %4.0f\n"% (step, self.tottime*2.418884326509e-17, dt*2.418884326509e-17, dtr*2.418884326509e-17, ecurr, cdf, nrand2,NN1.count("A"),NN1.count("M"),NN1.count("S"),NN2.count("A"),NN2.count("M"),NN2.count("S")) )
                self.file_KMC_vacancyneigh.write(" %4.0f %4.0f %4.0f %4.0f %4.0f %4.0f "% (NN1.count("A"),NN1.count("M"),NN1.count("S"),NN2.count("A"),NN2.count("M"),NN2.count("S")) )
            self.file_KMC_vacancyneigh.write("\n")
            self.file_KMC_vacancyneigh.flush()
        return

    def KMC_cluster_analysis(self,step=None,verbose=False):
        ''' determine size of clusters and mg/si content
            KMC_cluster_analysis(nsi=self.nsi,nmg=self.nmg,nvac=self.nvac,idx=self.idx,state=self.state,ridx=self.ridx,neigh=self.neigh)
        '''
        #all_solute_idx = np.zeros(nsi+self.nmg)
        #print('>> step',step) #,"step 1813 makes problems")
        allpairs = np.empty(((self.nsi+self.nmg)*12,2))
        allpairs[:] = np.nan
        running_idx = 0
        for solute in range(self.nsi + self.nmg): # is already the index
            sidx = self.idx[solute] # lattice site associated with this solute
            #all_solute_idx[solute] = sidx
            for sneigh in self.neigh[sidx]:
                running_idx += 1
                if self.state[sneigh] != 'A':
                    #print(">>++",self.state[sidx],'solute:',solute,'now sidx',sidx,'neighbor:',self.state[sneigh],'sneigh',sneigh,'ridx',self.ridx[sneigh])
                    allpairs[running_idx] = [sidx,sneigh]
                    #print('running_idx',running_idx,[sidx,sneigh])

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
            # this need to be separated from step==0 for RESTART's
            if not hasattr(self, 'file_KMC_cluster_sizes'):
                self.file_KMC_cluster_sizes = open("KMC_cluster_sizes","a+")
                self.file_KMC_cluster_mg = open("KMC_cluster_mg","a+")
                self.file_KMC_cluster_si = open("KMC_cluster_si","a+")
                self.file_KMC_cluster_vac = open("KMC_cluster_vac","a+")

            if step == 0:
                if os.path.isfile(file_KMC_cluster_sizes): os.remove(file_KMC_cluster_sizes)
                if os.path.isfile(file_KMC_cluster_mg): os.remove(file_KMC_cluster_mg)
                if os.path.isfile(file_KMC_cluster_si): os.remove(file_KMC_cluster_si)
                if os.path.isfile(file_KMC_cluster_vac): os.remove(file_KMC_cluster_vac)
                self.file_KMC_cluster_sizes.write("# column   1     --> step\n")
                self.file_KMC_cluster_sizes.write("# column   2 ... --> cluster sizes\n")

                self.file_KMC_cluster_mg.write("# column   1     --> step\n")
                self.file_KMC_cluster_mg.write("# column   2 ... --> amount Mg atoms in cluster \n")

                self.file_KMC_cluster_si.write("# column   1     --> step\n")
                self.file_KMC_cluster_si.write("# column   2 ... --> amount Si atoms in cluster \n")

                self.file_KMC_cluster_vac.write("# column   1     --> step\n")
                self.file_KMC_cluster_vac.write("# column   2 ... --> amount Vacancies in cluster \n")

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
            self.prng.rng.shuffle(idx)
            self.idx = idx


        # initialize state based on the index
        # this is a string indicating the state of the system. each character corresponds to a lattice site
        # and can be S=Si, M=Mg, A=Al, V=Vac.
        # create a random string
        names = np.asarray(list("S" * self.nsi + "M" * self.nmg + (self.natoms - self.nsi - self.nmg) * "A" +  "V" * self.nvac))
        state = names.copy()
        #print('names',names)
        #print('self.idx',self.idx)
        #print('state (0)',state)
        #print('state[7] (0)',state[7])
        #print('state[25] (0)',state[25])
        #print('state[48] (0)',state[48])

        # this maps the string to random sites
        state[self.idx] = names  #backshuffle!
        #print('state (1)',state)
        state = "".join(state)
        #print('state (1)',state)

        # reverse lookup index [i.e. ridx[i] gives the index of the atoms at site i]
        self.ridx = np.zeros(self.nsites,int)
        self.ridx[self.idx] = range(self.nsites)

        self.state = np.asarray(list(state))
        print('>>     ',"".join(self.state))
        print('>>     ','this is how one can get the positions of the solutes')
        print('>>     ','self.idx[0]',self.idx[0]  ," self.state[self.idx[0]] == self.state[7] :",self.state[self.idx[0]],'self.ridx[7] ',self.ridx[self.idx[0]])
        print('>>     ','self.idx[1]',self.idx[1]  ,"self.state[self.idx[1]] == self.state[25]:",self.state[self.idx[1]],'self.ridx[25]',self.ridx[self.idx[1]])
        print('>>     ','self.idx[2]',self.idx[2]  ,"self.state[self.idx[2]] == self.state[53]:",self.state[self.idx[2]],'self.ridx[53]',self.ridx[self.idx[2]])
        print('>>     ','self.idx[3]',self.idx[3]  ,"self.state[self.idx[3]] == self.state[40]:",self.state[self.idx[3]],'self.ridx[40]',self.ridx[self.idx[3]])
        print('>>     ','self.idx[4]',self.idx[4]  ,"self.state[self.idx[4]] == self.state[22]:",self.state[self.idx[4]],'self.ridx[22]',self.ridx[self.idx[4]])
        print('>>     ','self.idx[5]',self.idx[5]  ,"self.state[self.idx[5]] == self.state[33]:",self.state[self.idx[5]],'self.ridx[33]',self.ridx[self.idx[5]])


        self.KMC_cluster_analysis()

        print(">>      CHECKING INITIAL ASSIGNMENTS (this is already randomized)")
        for i in xrange(self.nsites):
            if self.ridx[i]<self.natoms:
                #if i < 4 or i > self.nsites-4 or self.state[i] in ['M','S']:
                if self.state[i] in ['M','S']:
                    print(i,self.beads.names[self.ridx[i]], self.state[i])
            else:
                print(i,"V", self.state[i])
            if self.idx[self.ridx[i]] != i:
                print("inconsistent site string for atom ",i, " and site ", self.ridx[i])

        if not f_restart:
            self.beads.q[0,:] = self.sites[self.idx].flatten() # also sets global q so we can use it further down

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
        self.feval = np.ones(self.neval,int)
        self._threadlock = threading.Lock()

    # threaded geometry optimization
    def geop_thread(self, ieval, nstr, nevent, ostr=None):
        self.geop[ieval].reset()
        ipot = self.dforces[ieval].pot

        #print('ieval',ieval)
        #print('nstr',nstr)
        #print('nevent',nevent)
        #print('ostr',ostr)
        for i in xrange(self.nstep):
            # print "geop ", i, self.dforces[ieval].pot
            self.geop[ieval].step(i)

            #print('--',i)
            #if self.geop[ieval].converged[0]: break

        newq = dstrip(self.dbeads[ieval].q[0]).copy()
        newpot =  self.dforces[ieval].pot

        # print "geop ", self.nstep, self.dforces[ieval].pot
        with self._threadlock:
            self.ecache[nstr] = newpot
            self.qcache[nstr] = newq
            self.ecache_n[nstr] = newpot
            self.qcache_n[nstr] = newq
            #self.newcache[nstr] = [newpot,newq]
            self.ncache += 1
            nevent[2] = self.ecache[nstr]
            nevent[3] = self.qcache[nstr]

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
        '''
        # computes TS energy by linearly interpolating initial & final state
        # interpolates between two structures considering PBCs
        '''
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
        ''' # finds first free evaluator '''
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
        '''
            # generates a starting lattice configuration that corresponds to a given state vector (a string of atomic types)
            # makes sure that the same occupation string corresponds to the same atom positions,
            # even though the actual atom indices might be different basically, makes the configurations
            # independent of same-atom permutations
        '''
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
        ''' # finds the best matching between the atoms in states state_1 and state_2,
            # assuming there is only one site swap
            # (should raise an error otherwise but it's not trivial without doing too many extra checks)
            # this basically says what is the motion of atoms from state 1 to state 2. This is useful
            # to interpolate between atomic coordinates in the two states, e.g. to get the TS geometry
        '''
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


    def step(self, step=None):
        kT = Constants.kb  * self.ens.temp
        # computes current energy (if not already stored)
        ostr = "".join(self.state)  # this is a unique id string that charactrizes the current state
        #print('step',step,ostr)
        self.tscache[ostr] = {}
        #print('step',step,self.state)
        self.KMC_cluster_analysis(step=step)

        if not ostr in self.ecache:
            self.dbeads[0].q[0,:] = self.sites[self.unique_idx(self.state)].flatten()
            rv = [0, 0, 0, 0, 0]
            self.geop_thread(0, ostr, rv)
            #self.beads.q[0,:] = self.dbeads[0].q[0,:] # also updates current position
            # self.forces.transfer_forces(self.dforces[0]) # forces have already been computed here...

        ecurr = self.ecache[ostr]

        # enumerates possible reactive events (vacancy swaps)
        levents = []
        in_cache = 0
        not_in_cache = 0
        ethreads = [None] * self.neval
        # loops over the vacancy
        for ivac in xrange(self.natoms, self.natoms + self.nvac):
            svac = self.idx[ivac] # lattice site associated with this vacancy
            #print('svac',svac)
            if self.state[svac] != "V":
                raise IndexError("Something got screwed and a vacancy state is not a vacancy anymore!")
            # loops over the neighbors of the selected vacancy
            for sneigh in self.neigh[svac]:
                # if the neighbor is a vacancy, move on. does not make sense to swap two vacancies!
                if self.state[sneigh] == "V" : continue

                # creates a new state vector with swapped atoms-vacancy and the associated label string
                nstate = self.state.copy()
                nstate[svac], nstate[sneigh] = self.state[sneigh], self.state[svac]
                nstr = "".join(nstate) # this is the string that corresponds to the new state
                if not nstr in self.ecache:
                    not_in_cache += 1
                    # new state, must compute!
                    # creates a swapped index
                    nidx = self.idx.copy()
                    if  self.ridx[svac]!= ivac or self.idx[ivac]!=svac:
                        raise  IndexError("Something got screwed and the index does not correspond anymore to site occupancy")

                    #ivac = self.ridx[svac]
                    ineigh = self.ridx[sneigh]   # gets index of atom associated with the neighboring site
                    nidx[ivac], nidx[ineigh] = self.idx[ineigh], self.idx[ivac]

                    ieval = self.find_eval(ethreads)
                    # launches evaluator
                    self.dbeads[ieval].q[0,:] = self.sites[self.unique_idx(nstate)].flatten()

                    nevent = [svac, sneigh, 0.0, 0.0, 0.0]

                    # runs a geometry optimization
                    #self.geop_thread(ieval=ieval, nstr=nstr, nevent=nevent)
                    st = threading.Thread(target=self.geop_thread, name=str(ieval), kwargs={"ieval":ieval, "nstr":nstr, "nevent" : nevent, "ostr": ostr})
                    st.daemon = True
                    st.start()
                    ethreads[ieval] = st
                else:
                    in_cache += 1
                    #print "Found state ", nstr, " retrieving cached energy ", self.ecache[nstr]

                    # fetch energy from previous calculation
                    nevent = [svac, sneigh, self.ecache[nstr], self.qcache[nstr], 0.0]

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
        for st in ethreads:
            while not st is None and st.isAlive():
                st.join(2)

        #print("Computed ", len(levents), " possible reactions. Cache len ", len(self.ecache))

        # get list of rates
        rates = np.zeros(len(levents), float)
        crates = np.zeros(len(levents), float)
        cdf = 0.0
        #self.kmcfile_energies.write("%18.11e"% (ecurr))

        for i in xrange(len(levents)):
            #print ("Barrier, naive: %f, static: %f" % (0.5*(ecurr + levents[i][2]) + self.diffusion_barrier_al, levents[i][4]))
            #print('levents[i][-1]',levents[i][-1],'barrier',self.barriers[levents[i][-1]],'prefact',self.prefactors[levents[i][-1]])  # currently we have only one barrier/prefactor, so final energy decides everything. (barrier: 0.01910964952, prefactor: 0.0025229174680000003)
            #print 'eurr',ecurr,'levents[i][2]',levents[i][2],'SUM',(ecurr + levents[i][2])
            #print '@@0.5dE',(0.5*(ecurr - levents[i][2])),'barrier',self.barriers[levents[i][-1]]

            #@@print '@@ b/dE',self.barriers[levents[i][-1]]/(0.5*(ecurr - levents[i][2]))
            ets = 0.5*(ecurr + levents[i][2]) + self.barriers[levents[i][-1]]      #naive heuristic for the ts energy
            #@@print "Event ", i, levents[i][-1], ecurr, ">>", ets, ">>", levents[i][2]
            rates[i] = self.prefactors[levents[i][-1]] * np.exp(-(ets-ecurr)/kT)   # this is k_x = k[i] in Voters paper   (rate constant k_{ij}
            cdf += rates[i]                                                        # this is k_{TOT} in Voters paper
            crates[i] = cdf                                                        # this is k_{TOT} however subdivided in smaller boxes k_x

        self.barriers_energies_analysis(step,ecurr,levents,rates)

        # KMC selection rule based on the rate (needs: crates,
        #                                       wobei cdf = crates[-1]
        fpick = self.prng.u*cdf
        isel = 0
        while fpick > crates[isel]:
            isel += 1
        nrand2 = self.prng.u
        dt = -1.0/cdf*np.log(1.0 - nrand2)

        self.vacancy_neighborhood_analysis(step,dt,ecurr,cdf,nrand2)
        #print ("Time spent %12.5e at %s nrg %12.5e" % (dt, ostr,ecurr))
        #print "Selected event ", isel, " with rate ", rates[isel], " / ", cdf

        if False:
            print('pos!!')
            print(self.dbeads[0].q[0,:])  # this is in bohrradius
            time.sleep(10)


        ###################################################################
        # create new state
        ###################################################################
        iev = levents[isel] # levents[self.prng.randint(len(levents))]
        svac, sneigh = iev[0], iev[1]
        ivac, ineigh = self.ridx[svac], self.ridx[sneigh]

        # does the swap (never reject, for the moment)
        self.state[svac], self.state[sneigh] = self.state[sneigh], self.state[svac]
        self.ridx[svac], self.ridx[sneigh] = self.ridx[sneigh], self.ridx[svac]
        self.idx[ivac], self.idx[ineigh] = self.idx[ineigh], self.idx[ivac]


        # we got a new configuration but the residence time is linked to the previous configuration so we output that
        self.kmcfile.write("%12.5e  %12.5e  %18.11e  %s\n"% (self.tottime, dt, ecurr, ostr) )
        self.kmcfile.flush()
        if False:
            print('self.tottime*2.418884326509e-17:',self.tottime*2.418884326509e-17)
            print('dt*2.418884326509e-17          :',dt*2.418884326509e-17)
            print('cdf                            :',cdf)
            print('ecurr                          :',ecurr)
            print('nrand2                         :',nrand2)
        print('>> step:',step,"computed:",not_in_cache,'from cache:',in_cache,"clusters_sizes:",self.cluster_sizes,"events: ", len(levents), " possible reactions. Cache len ", len(self.ecache))
        #print('ivac',ivac,'svac = self.idx[ivac]',self.idx[ivac])


        self.tottime += dt
        self.ens.time += dt  # updates time counter
        #print  "Finishing step at ", "".join(self.state)

        # updates the positions
        self.cell.h = self.dcell.h

        uidx = self.unique_idx(self.state)
        ruidx = np.zeros(self.nsites,int)
        ruidx[uidx] = range(self.nsites)

        self.sites[self.unique_idx(self.state)]
        oldq = dstrip(self.beads.q[0]).copy()

        newq = np.zeros(self.nsites*3, float)
        # we want continuity (modulo PBC jumps, that we'll take care of later...)
        for i in xrange(self.nsites):
            # in which site sits atom i?
            isite = self.idx[i]
            # which atom sits in this site in the unique-mapped structure?
            iuid = ruidx[self.idx[i]]
            newq[3*i:3*i+3] = iev[3][3*iuid:3*(iuid+1)]
        newq-=oldq
        self.cell.array_pbc(newq)
        self.beads.q[0]+=newq


