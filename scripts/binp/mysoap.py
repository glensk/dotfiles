#!/usr/bin/env python
from __future__ import print_function
import numpy as np
import os,sys
import quippy
import myutils as my
try:
    from tqdm import tqdm_notebook, tnrange
except ImportError:
    pass

##################################################################################
## SOAP funcions
##################################################################################
def get_rawsoap(at, gs=0.3, co=3.0, cotw=0.5, nmax=9, lmax=6, cw=1.0, nrm=True, zlist=[], central_z=0, x_cmd=""):
    """
    e.g.: zliat, soaps = get_rawsoap(at, nmax=nmax, lmax=lmax, co=co, gs=gs, cotw=cotw, nrm=nrm, cw=1, zlist=zlist, central_z=central_z,
                                   x_cmd="cutoff_dexp=2 cutoff_scale=3.0 cutoff_rate=2.0")

    gs gaussian widths
    co cutoff radius [Angstrom]
    cotw you don't care
    n radial basis functions
    l spherical harmonics basis
    zlist list of atomic numbers the SOAPS will see
    central_z list of atomic numbers you will center on : IT DEFINES THE NUMBER OF ROWS of the rawsoaps

    """
    #np.savetxt("/Users/glensk/tmppos.dat",at.positions)
    #pos = np.loadtxt("/Users/glensk/tmppos.dat")
    #maxdiff = np.abs(pos-at.positions).max()
    #print(at.symbols)
    #print(at.numbers)
    #print(at.positions[:2])
    #print('nmax:',nmax,'lmax:',lmax,'co',co,'gs',gs,'cotw',cotw,'nrm',nrm,'cw',cw,'zlist',zlist,'central_z',central_z) #,'x_cmd',x_cmd)
    #print('nmax:',nmax,'lmax:',lmax,'co',co,'gs',gs,'cotw',cotw,'nrm',nrm,'cw',cw,'zlist','x_cmd',x_cmd)
    #print(gs, co, cotw, nmax, lmax, cw, nrm, zlist, central_z, x_cmd,at.numbers,at.symbols,at.cell) #,maxdiff)
    if len(zlist)==0:
        zlist = np.unique(np.asarray(at.z))
    nz = len(zlist)
    lspecies = "n_species="+str(nz)+" species_Z={"+(" ".join(np.asarray(zlist,dtype=str)))+"}"
    soapstr=("soap central_reference_all_species=F normalise="+("T" if nrm else "F")+" central_weight="+str(cw)+
               "  covariance_sigma0=0.0 atom_sigma="+str(gs)+" cutoff="+str(co)+" cutoff_transition_width="+str(cotw)+
               " n_max="+str(nmax)+" l_max="+str(lmax)+' '+lspecies+' Z='+str(central_z)+' '+x_cmd)
    #print(':'+soapstr+":")
    desc = quippy.descriptors.Descriptor(soapstr )
    #print(co,soapstr)
    #print(quippy.__file__)
    #print(quippy.descriptors.__file__)
    at.set_cutoff(co);
    at.calc_connect();
    return zlist, desc.calc(at)["descriptor"]

def get_soap2_vec(rawsoap, zlist, zglobal=None,nmax=None,lmax=None):
    if nmax is None:
        sys.exit('Please provide nmax!')
    if lmax is None:
        sys.exit('Please provide lmax!')
    isoap = 0
    isqrttwo = 1.0/np.sqrt(2.0)
    #njsoap = {}
    #ipair = {}

    #print('zglobal:',zglobal)
    #print('zlist  :',zlist,'ka')
    if zglobal is None:
        zglobal = zlist
    #print('len(zglobal)',len(zglobal),nmax,len(zglobal),lmax)
    njsoap = np.zeros((len(zglobal),nmax,len(zglobal),nmax,lmax+1))
    zmap = np.zeros(len(zlist),int)
    for s in xrange(len(zlist)):
        zmap[s] = zglobal.index(zlist[s])
    for s1 in xrange(len(zlist)):
        z1 = zmap[s1]
        for n1 in xrange(nmax):
            for s2 in xrange(s1+1):
                z2 = zmap[s2]
                for n2 in xrange(nmax if s2<s1 else n1+1):
                    soap_element = rawsoap[isoap:isoap+lmax+1]
                    if (s1 != s2 or n1 != n2):
                        soap_element *= isqrttwo  # undo the normalization since we will actually sum over all pairs in all directions!
                    njsoap[z1,n1,z2,n2,:] = soap_element
                    njsoap[z2,n2,z1,n1,:] = soap_element
                    isoap+=lmax+1
    return njsoap

def do_fps(x, d=0):
    if d == 0 : d = len(x)
    n = len(x)
    iy = np.zeros(d, int)
    # faster evaluation of Euclidean distance
    n2 = np.sum(x**2,axis=1)
    iy[0] = np.random.randint(0, n)
    dl = n2 + n2[iy[0]] - 2* np.dot(x, x[iy[0]])
    lmin = np.zeros(d)
    for i in range(1,d):
        iy[i] = np.argmax(dl)
        lmin[i-1] = dl[iy[i]]
        if i%1000 ==0: print("max min dist %f" %( dl[iy[i]]))
        nd = n2 + n2[iy[i]] - 2*np.dot(x,x[iy[i]])
        dl = np.minimum(dl, nd)
    return iy, lmin

def get_soaps(kmcxyz = False, nmax = 8, lmax = 6, co = 4, gs = 0.5, zlist = [12,13,14], central_z=23,verbose=False,showtdqm=True,nrm=True ):
    ''' co          = cutoff [Angstrom]
        gs gaussian widths
    '''
    ###################################
    # first get rsoap
    ###################################
    rsoap = [] # is done by rsoap
    #nmax = 8   # 12
    #lmax = 6  # 9
    #co = 4  # cutoff radius, 3.0  ## is this in angstrom?
    #gs = 0.5 # gaussian width, 0.3
    cotw = gs
    #det = "c"+str(co)+"-g"+str(gs)+"-nrm-cw"
    #zlist = [12,13,14]  # why is this working, even without specifying the vacancy/Vanadium ....? [12,13,14,23]
    ztot = len(zlist)
    nsoap = ztot**2 * nmax**2 * (1+lmax)
    ntot = len(kmcxyz)
    if verbose:
        print('rsoap ...')
    if showtdqm == True:
        gothrough = tqdm_notebook(kmcxyz)
    else:
        gothrough = kmcxyz

    #print('nmax:',nmax,'lmax:',lmax,'co',co,'gs',gs,'cotw',cotw,'nrm',nrm,'cw',1,'zlist',zlist,'central_z',central_z)
    for at in gothrough:
        #progress(i,len(kmcxyz))
        # andrea used a zentral_z=0
        # the next line could be potentially done parallel...
        #print('type',type(at))  # type <class 'quippy.atoms.Atoms'>
        zliat, soaps = get_rawsoap(at, nmax=nmax, lmax=lmax, co=co, gs=gs, cotw=cotw, nrm=nrm, cw=1, zlist=zlist, central_z=central_z, x_cmd="cutoff_dexp=2 cutoff_scale=3.0 cutoff_rate=2.0")
        rsoap.append((zliat, soaps))
    if verbose:
        print('rsoap done ...')
    print('rsoap in get_soaps (1)',rsoap)

    ###################################
    # now get esoap2 (from rsoap)
    ###################################
    soap2 = []  # is done by soap2vec
    if verbose:
        print('soap2 ...')

    if showtdqm == True:
        gothrough = tnrange(ntot)
    else:
        gothrough = range(ntot)

    for i in gothrough:
        #progress(i,ntot)
        izl, soaps = rsoap[i]
        print('izl',izl,'soaps',soaps)
        lenv = np.zeros((len(soaps), ztot, nmax, ztot, nmax, 1+lmax))
        for s in range(len(soaps)):
            lenv[s] = get_soap2_vec(soaps[s], izl, zlist,nmax=nmax,lmax=lmax)
        soap2.append(lenv)
    if verbose:
        print('soap2 done ...')
        print('esoap2 ...')
    print('rsoap in get_soaps (2)',rsoap)
    esoap2 = []
    for s in soap2:
        #progress(i,len(soap2))
        for e in s:
            esoap2.append(e.flatten())
    if verbose:
        print('esoap2 done ...')
    esoap2 = np.asarray(esoap2)
    return rsoap, esoap2


def comparesoap(rsoap,rsoap_new,tags=False):
    ''' rsoap: are all the previously known rsoaps '''
    print("idx   close   exact   diff:  xxdiff to rsoap_new")
    for idx in range(len(rsoap)): ## idx is the number of the soap
        maxdiff = 0.
        for i in range(len(rsoap[0][1][0])): ## i is the index of every soap element
            diff = rsoap[idx][1][0][i]-  rsoap_new[1][0][i]
            if np.abs(diff) > maxdiff:
                maxdiff = np.abs(diff)
                #print('maxdiff',maxdiff)
        isclose = np.alltrue(np.isclose(rsoap[idx][1][0],rsoap_new[1][0]))
        isexact = np.alltrue(rsoap[idx][1][0]==rsoap_new[1][0])
        #isclose = "??"
        #print str(idx).ljust(5),str(isclose).ljust(7),"diff:", str(np.round(maxdiff,7)).ljust(10),str(np.round(maxdiff2,7)).ljust(10),"T:",tags[idx].ljust(25),"M",str(rsoap[idx][1][0].max()).ljust(15),str(rsoap_new[1][0].max()).ljust(15)
        if type(tags) == bool: tags_ = "--"
        else: tags_ = str(tags[idx])

        #print('ii',idx,isclose,maxdiff,maxdiff2,'tag',str(tags[idx]))
        #print(str(idx).ljust(5),str(isclose).ljust(7),"diff:", str(np.round(maxdiff,7)).ljust(10),str(np.round(maxdiff2,7)).ljust(10),"T:",str(tags_).ljust(25),"M",rsoap[idx][1][0].max()-rsoap_new[1][0].max())
        if isclose == True:
            print(my.printgreen(str(idx).ljust(5)+str(isclose).ljust(7)+str(isexact).ljust(7)+tags_.ljust(30)+"maxdiff:"+str(maxdiff)))
        else:
            print(my.printred(str(idx).ljust(5)+str(isclose).ljust(7)+str(isexact).ljust(7)+tags_.ljust(30)+"maxdiff:"+str(maxdiff)))

    return


def checksoap():
    inp1 = quippy.AtomsList("/Users/glensk/tmp/dataxx.quippy.xyz")  # cubic lattice with 1st and 2NN around the vacancy at [0,0,0]
    inp1_0 = quippy.AtomsList("/Users/glensk/tmp/dataxx0.quippy.xyz")

    filename = fileName= "/Users/glensk/tmp/dataxx.name"
    tags = [line.rstrip('\n') for line in open(filename)]

    co              = 10.6 # cutoff
    co              = 10.6 # cutoff
    gs              = 1.5 # gaussian_width = 0.5
    nmax            = 2 # quadratisch 12
    lmax            = 1 # linear 9,4
    showtdqm        = False
    verbose         = False
    showMatrix      = False
    oldscatterplot  = False
    zlist           = [12,13,14]
    #zlist          = [13]
    central_z       = 23
    x_cmd           = "cutoff_dexp=2 cutoff_scale=3.0 cutoff_rate=2.0"

    rsoap_a = []
    for i in range(len(inp1)):
        rsoap_a.append(get_rawsoap(inp1[i], gs=gs, co=co, cotw=gs, nmax=nmax, lmax=lmax, cw=1, nrm=True, zlist=zlist, central_z=central_z, x_cmd=x_cmd))
    rsoap_e0         = get_rawsoap(inp1[0], gs=gs, co=co, cotw=gs, nmax=nmax, lmax=lmax, cw=1, nrm=True, zlist=zlist, central_z=central_z, x_cmd=x_cmd)
    rsoap_e1         = get_rawsoap(inp1[1], gs=gs, co=co, cotw=gs, nmax=nmax, lmax=lmax, cw=1, nrm=True, zlist=zlist, central_z=central_z, x_cmd=x_cmd)
    rsoap_e2         = get_rawsoap(inp1[2], gs=gs, co=co, cotw=gs, nmax=nmax, lmax=lmax, cw=1, nrm=True, zlist=zlist, central_z=central_z, x_cmd=x_cmd)

    print()
    print('tt',type(rsoap_a))
    print()
    print('how many structrues known:',len(rsoap_a),"correct: 0?")
    if False: sys.exit('99')

    #print(rsoap_a[0][1][0][:20])
    comparesoap(rsoap_a,rsoap_new=rsoap_a[0],tags=tags)
    print('a',len(rsoap_a),"correct: 1? (and therefore also 3)")
    #print(rsoap_a[1][1][0][:20])
    comparesoap(rsoap_a,rsoap_new=rsoap_a[1],tags=tags)
    print('e0',len(rsoap_a),"correct: 0?")
    comparesoap(rsoap_a,rsoap_new=rsoap_e0,tags=tags)

    print('e1',len(rsoap_a),"correct: 1? (and therefore also 3)")
    comparesoap(rsoap_a,rsoap_new=rsoap_e1,tags=tags)
    print('e2',len(rsoap_a),"correct: 2? (and therefore also 4)")
    comparesoap(rsoap_a,rsoap_new=rsoap_e2,tags=tags)
    print('a[5]',len(rsoap_a),"correct: 2? (and therefore also 4)")
    comparesoap(rsoap_a,rsoap_new=rsoap_a[14],tags=tags)

    ##########################################################
    # now create structure from ase (or otherwise, from ipi)
    ##########################################################
    # a) get ase structure
    print('##########')
    print('get_sphere')
    print('##########')
    atomsc_sphere_ase = my.create_al_sphere(a0=4.05,matrix_element="Al",cubic=True,ncell=4,nvac=1,cutoff=4.05,vacidx=0)
    frame1_ase = atomsc_sphere_ase.copy()
    frame1_ase[1].symbol = "Si"
    frame2_ase = atomsc_sphere_ase.copy()
    frame2_ase[3].symbol = "Mg"

    # b) convert ase structure to quippy
    from quippy.atoms import Atoms as QuippyAtoms
    q0 = QuippyAtoms(atomsc_sphere_ase)
    q1 = QuippyAtoms(frame1_ase)
    q2 = QuippyAtoms(frame2_ase)
    rsoap_q0         = get_rawsoap(q0, gs=gs, co=co, cotw=gs, nmax=nmax, lmax=lmax, cw=1, nrm=True, zlist=zlist, central_z=central_z, x_cmd=x_cmd)
    rsoap_q1         = get_rawsoap(q1, gs=gs, co=co, cotw=gs, nmax=nmax, lmax=lmax, cw=1, nrm=True, zlist=zlist, central_z=central_z, x_cmd=x_cmd)
    rsoap_q2         = get_rawsoap(q2, gs=gs, co=co, cotw=gs, nmax=nmax, lmax=lmax, cw=1, nrm=True, zlist=zlist, central_z=central_z, x_cmd=x_cmd)
    print('rsoap_q0:',rsoap_q0)
    print('rsoap_q1:',rsoap_q1)
    print('rsoap_q2:',rsoap_q2)
    print('q0',len(rsoap_a),"correct: 2? (and therefore also 4)")
    comparesoap(rsoap_a,rsoap_new=rsoap_q0,tags=False)
    print()
    comparesoap(rsoap_a,rsoap_new=rsoap_q1,tags=False)
    print()
    comparesoap(rsoap_a,rsoap_new=rsoap_q2,tags=False)
    print(rsoap_q2)
    sys.exit()

    def get_descriptor(at):
        soapstr = "soap central_reference_all_species=F normalise=T central_weight=1  covariance_sigma0=0.0 atom_sigma=0.5 cutoff=10.6 cutoff_transition_width=0.5 n_max=5 l_max=4 n_species=3 species_Z={12 13 14} Z=23 cutoff_dexp=2 cutoff_scale=3.0 cutoff_rate=2.0"
        soapstr = "soap central_reference_all_species=F normalise=T central_weight=1  covariance_sigma0=0.0 atom_sigma=0.5 cutoff=10.6 cutoff_transition_width=0.5 n_max=2 l_max=1 n_species=3 species_Z={12 13 14} Z=23 cutoff_dexp=2 cutoff_scale=3.0 cutoff_rate=2.0"
        desc = quippy.descriptors.Descriptor(soapstr )
        co = 10.6
        at.set_cutoff(co);
        at.calc_connect();
        descript = desc.calc(at)["descriptor"]
        return descript

    desc1 = get_descriptor(inp1[0])
    print('from get_descriptor(inp1[0])')
    print(desc1[0])

    print('e1',len(rsoap_e1),"0?")
    comparesoap(rsoap_a,rsoap_new=rsoap_e1[0],tags=False)
    print('e2',len(rsoap_e2),"0?")
    comparesoap(rsoap_a,rsoap_new=rsoap_e2,tags=False)
    return

if __name__ == '__main__':
    checksoap()
