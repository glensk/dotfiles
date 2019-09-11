import os
import sys
import numpy as np
import copy

def ListToStr(varlist):
    mystr = varlist[0]
    for ii, vstr in enumerate(varlist):
        if not(ii==0):
            mystr = mystr + " " + vstr

    return mystr

class crystal( object ):
    '''
        General crystal class
    '''
    def __init__( self ):
        '''
        '''
        # title of the supercell
        self.title = 'title placeholder'
        self.crystaltype = []

        # supercell lattice vectors
        self.scale   = 1.0               # scale of cellvec
        self.cellvec = np.zeros([3,3])   # defines the supercell box

        # cartesian coordinates of atoms
        self.xcar = []
        self.ycar = []
        self.zcar = []
        self.rcar = []

        # relative coordinates of atoms
        self.xrel = []
        self.yrel = []
        self.zrel = []
        self.rrel = []

        # atom types and ids
        self.itype  = []                 # interger type
        self.stype  = []                 # string type
        self.nid    = []
        self.natoms = []
        self.ntypes = []

    def rotate_atoms( self, Rmatrix ):
        '''
            This function rotates the crystal by Rmatrix
        '''

        # rotate atoms
        rcar_rotated    = np.dot(self.rcar   , Rmatrix)
        cellvec_rotated = np.dot(self.cellvec, Rmatrix)
        #cellvec_rotated = cellvec_rotated.T

        # export variables
        self.xcar  = rcar_rotated[:,0]
        self.ycar  = rcar_rotated[:,1]
        self.zcar  = rcar_rotated[:,2]
        self.cellvec = cellvec_rotated

        # update position variabls
        self.update_rcar_from_xyzcar()
        self.update_xyzrel_from_xyzcar()
        self.update_rrel_from_rcar()

    def reflect_atoms( self ):
        '''
            Reflects atoms outside the cell back into the cell
        '''
        idref = np.nonzero( self.xrel < 0 )[0]
        self.xrel[idref] = self.xrel[idref] + 1
        idref = np.nonzero( self.xrel > 1 )[0]
        self.xrel[idref] = self.xrel[idref] - 1

        idref = np.nonzero( self.yrel < 0 )[0]
        self.yrel[idref] = self.yrel[idref] + 1
        idref = np.nonzero( self.yrel > 1 )[0]
        self.yrel[idref] = self.yrel[idref] - 1

        idref = np.nonzero( self.zrel < 0 )[0]
        self.zrel[idref] = self.zrel[idref] + 1
        idref = np.nonzero( self.zrel > 1 )[0]
        self.zrel[idref] = self.zrel[idref] - 1

        # update atoms coordinates
        self.update_rrel_from_xyzrel
        self.update_xyzcar_from_xyzrel()
        self.update_rcar_from_xyzcar()

    def translate_atoms_cart( self, rcar_trans ):
        '''
            Translates atoms in the cell by cartesian translation vector rcar_trans
        '''
        self.rcar = self.rcar - rcar_trans
        self.update_xyzcar_from_rcar()
        self.update_rrel_from_rcar()
        self.update_xyzrel_from_xyzcar()

    def center_atoms_rel( self ):
        '''
            Centers the supercell about zero
        '''
        idshift = np.nonzero( self.xrel >  0.5 )[0]
        self.xrel[idshift]  = self.xrel[idshift] - 1
        idshift = np.nonzero( self.yrel >  0.5 )[0]
        self.yrel[idshift]  = self.yrel[idshift] - 1
        idshift = np.nonzero( self.zrel >  0.5 )[0]
        self.zrel[idshift]  = self.zrel[idshift] - 1
        idshift = np.nonzero( self.xrel < -0.5 )[0]
        self.xrel[idshift]  = self.xrel[idshift] + 1
        idshift = np.nonzero( self.yrel < -0.5 )[0]
        self.yrel[idshift]  = self.yrel[idshift] + 1
        idshift = np.nonzero( self.zrel < -0.5 )[0]
        self.zrel[idshift]  = self.zrel[idshift] + 1

        # update other coordinates
        self.update_rrel_from_xyzrel()
        self.update_xyzcar_from_xyzrel()
        self.update_rcar_from_rrel()

    def update_rcar_from_xyzcar( self ):
        '''
            Updates r-cartesian based on x-,y-,z-cartesian coordinates
        '''
        self.rcar = np.zeros([ len(self.xcar), 3 ])
        self.rcar[:,0] = self.xcar
        self.rcar[:,1] = self.ycar
        self.rcar[:,2] = self.zcar

    def update_xyzcar_from_rcar( self ):
        '''
            Updates x-,y-,z-cartesian based on r-cartesian coordinates
        '''
        self.xcar = np.zeros([ np.shape(self.rcar)[0], 1 ])
        self.ycar = np.zeros([ np.shape(self.rcar)[0], 1 ])
        self.zcar = np.zeros([ np.shape(self.rcar)[0], 1 ])
        self.xcar = self.rcar[:,0]
        self.ycar = self.rcar[:,1]
        self.zcar = self.rcar[:,2]

    def update_rrel_from_xyzrel( self ):
        '''
            Updates r-relative based on x-,y-,z-relative coordinates
        '''
        self.rrel = np.zeros([ len(self.xrel), 3 ])
        self.rrel[:,0] = self.xrel
        self.rrel[:,1] = self.yrel
        self.rrel[:,2] = self.zrel

    def update_xyzrel_from_rrel( self ):
        '''
            Updates x-,y-,z-relative based on r-relative coordinates
        '''
        self.xrel = np.zeros([ np.shape(self.rrel)[0], 1 ])
        self.yrel = np.zeros([ np.shape(self.rrel)[0], 1 ])
        self.zrel = np.zeros([ np.shape(self.rrel)[0], 1 ])
        self.xrel = self.rrel[:,0]
        self.yrel = self.rrel[:,1]
        self.zrel = self.rrel[:,2]

    def update_rrel_from_rcar( self ):
        '''
            Updates r-relative based on r-cartesian coordinates
        '''
        cellvec_inv = np.linalg.inv( self.cellvec )
        self.rrel   = np.dot( self.rcar, cellvec_inv )

    def update_rcar_from_rrel( self ):
        '''
            Updates r-cartesian based on r-relative coordinates
        '''
        self.rcar   = np.dot( self.rrel, self.cellvec )

    def update_xyzrel_from_xyzcar( self ):
        '''
            Updates xyz-relative based on xyz-cartesian coordinates
        '''
        self.update_rcar_from_xyzcar()
        self.update_rrel_from_rcar()
        self.update_xyzrel_from_rrel()

    def update_xyzcar_from_xyzrel( self ):
        '''
            Updates xyz-cartesian based on xyz-relative coordinates
        '''
        self.update_rrel_from_xyzrel()
        self.update_rcar_from_rrel()
        self.update_xyzcar_from_rcar()

    def updateCellDetails( self ):
        '''
            Updates crystal details (natoms, ntypes, nid)
        '''
        self.natoms = int( len( self.xcar ) )
        self.nid    = np.array( list(range( self.natoms)) ) + 1
        self.ntypes = int( np.max( self.itype ) )

    def deleteAtoms( self, idnull ):
        '''
            Deletes elements with indices contained in idnull
        '''
        self.xcar = np.delete( self.xcar, idnull )
        self.ycar = np.delete( self.ycar, idnull )
        self.zcar = np.delete( self.zcar, idnull )
        self.xrel = np.delete( self.xrel, idnull )
        self.yrel = np.delete( self.yrel, idnull )
        self.zrel = np.delete( self.zrel, idnull )
        self.rcar = np.delete( self.rcar, idnull, 0 )
        self.rrel = np.delete( self.rrel, idnull, 0 )

        self.itype = np.delete( self.itype, idnull )
        self.stype = np.delete( self.stype, idnull )
        self.nid   = np.delete( self.nid  , idnull )

        self.updateCellDetails()

    def exportLammps( self, directory, filename = None ):
        '''
            Exports data to lammps format
            NOTE: only works for orthagonal cell vectors
        '''
        # initialize
        self.xlo = 0.0
        self.ylo = 0.0
        self.zlo = 0.0
        self.xhi = self.cellvec[0,0]
        self.yhi = self.cellvec[1,1]
        self.zhi = self.cellvec[2,2]

        # open file
        fh = open( directory + "/" + filename, 'w' )
        #
        # write title
        fh.write( self.title + "\n" )
        fh.write( "%8i atoms \n"             %( self.natoms ) )
        fh.write( "%8i atom types \n"        %( self.ntypes ) )
        fh.write( "%10.5f %10.5f xlo xhi \n" %( self.xlo, self.xhi ) )
        fh.write( "%10.5f %10.5f ylo yhi \n" %( self.ylo, self.yhi ) )
        fh.write( "%10.5f %10.5f zlo zhi \n" %( self.zlo, self.zhi ) )
        fh.write( "Atoms \n" )
        fh.write( "\n" )
        for ii in range(self.natoms):
            fh.write("%6i %3i %20.15f %20.15f %20.15f \n" % \
                    (self.nid[ii] , self.itype[ii], self.xcar[ii], \
                     self.ycar[ii], self.zcar[ii]) )
        # close file
        fh.close()

        return

    def exportLammpsTriclinic( self, directory, filename = None ):
        '''
            Exports data to lammps format
            NOTE: only works for orthagonal cell vectors
        '''

        ''' Create rotation matrix '''
        a1 = self.cellvec[0,:]
        a2 = self.cellvec[1,:]
        a3 = self.cellvec[2,:]
        ep1 = a1 / np.linalg.norm(a1)
        ep2 = ( a2 - (np.dot(a2,ep1))*ep1 ); ep2 = ep2/np.linalg.norm(ep2)
        ep3 = np.cross( ep1, ep2 )

        Rinv = np.array([ep1, ep2, ep3])
        Rmat = Rinv.T

        ''' Rotate unit cell '''
        self.rotate_atoms( Rmat )

        ''' initialize output '''
        xlo = 0.0
        ylo = 0.0
        zlo = 0.0
        xhi = self.cellvec[0,0]
        yhi = self.cellvec[1,1]
        zhi = self.cellvec[2,2]
        xy  = self.cellvec[1,0]
        xz  = self.cellvec[2,0]
        yz  = self.cellvec[2,1]

        # open file
        fh = open( directory + "/" + filename, 'w' )
        #
        # write title
        fh.write( self.title + "\n" )
        fh.write( "%8i atoms \n"             %( self.natoms ) )
        fh.write( "%8i atom types \n"        %( self.ntypes ) )
        fh.write( "%15.10f %15.10f     xlo xhi \n" %( xlo, xhi ) )
        fh.write( "%15.10f %15.10f     ylo yhi \n" %( ylo, yhi ) )
        fh.write( "%15.10f %15.10f     zlo zhi \n" %( zlo, zhi ) )
        fh.write( "%15.10f %15.10f %15.10f xy xz yz \n" %( xy, xz, yz ) )
        fh.write( "Atoms \n" )
        fh.write( "\n" )
        for ii in range(self.natoms):
            fh.write("%6i %3i %20.15f %20.15f %20.15f \n" % \
                    (self.nid[ii], self.itype[ii], \
                     self.xcar[ii], self.ycar[ii], self.zcar[ ii]) )
        # close file
        fh.close()

        return

    def read_lammps_input_file( self, directory, filename = None ):
        '''
            Reads existing Lammps input file into the class
            NOTE: only works for orthagonal cells
        '''
        fh    = open( directory + "/" + filename, 'r' )
        lines = fh.readlines()
        fh.close()

        # assign data to variables
        for ii, line in enumerate(lines):
            split = line.split()
            if   ( ii == 0 ):
                self.title = ListToStr(split)
            elif ( ii == 1 ):
                self.natoms = int( split[0] )
            elif ( ii == 2 ):
                self.ntypes = int( split[0] )
            elif ( ii == 3 ):
                self.xlo = float( split[0] )
                self.xhi = float( split[1] )
            elif ( ii == 4 ):
                self.ylo = float( split[0] )
                self.yhi = float( split[1] )
            elif ( ii == 5 ):
                self.zlo = float( split[0] )
                self.zhi = float( split[1] )
            elif ( ii in range(8,8+self.natoms) ):
                self.nid.append( int( split[0] ) )
                self.itype.append( int( split[1] ) )
                self.xcar.append( float(split[2]) )
                self.ycar.append( float(split[3]) )
                self.zcar.append( float(split[4]) )
                self.rcar.append( [ float(split[2]), float(split[3]), float(split[4]) ] )

        # convert to numpy array
        self.nid   = np.array( self.nid   )
        self.itype = np.array( self.itype )
        self.xcar  = np.array( self.xcar  )
        self.ycar  = np.array( self.ycar  )
        self.zcar  = np.array( self.zcar  )
        self.rcar  = np.array( self.rcar  )
        self.cellvec[0,0] = self.xhi - self.xlo
        self.cellvec[1,1] = self.yhi - self.ylo
        self.cellvec[2,2] = self.zhi - self.zlo

    def read_POSCAR( self, directory, filename = None ):
        '''
            Reads poscar file
        '''
        fh    = open( directory + "/" + filename, 'r' )
        lines = fh.readlines()
        fh.close()

        cellvec = []
        self.dof = []
        self.vel = []

        # assign data to list
        for ii, line in enumerate(lines):
            if ii == 0:
                ll = line.split()
                self.title = ListToStr(ll)
            elif ii == 1:
                self.scale = float(line.split()[0])
            elif ii in [2,3,4]:
                ll = line.split()
                cellvec.append([ float(ll[0]), float(ll[1]), float(ll[2]) ])
            elif ii in [5]:
                self.atype = line.split()
            elif ii in [6]:
                self.ntypes = list(map(int,line.split()))
            elif ii in [7]:
                self.method = ListToStr(line.split())
            elif ii in [8]:
                self.coord = line.split()
            elif ii in range(9,9+np.sum(self.ntypes)):
                ll = line.split()
                self.rrel.append([ float(ll[0]), float(ll[1]), float(ll[2]) ])
                self.dof.append(ll[3:6])
            elif ii in range(10+np.sum(self.ntypes),10+2*np.sum(self.ntypes)):
                ll = line.split()
                self.vel.append([ float(ll[0]), float(ll[1]), float(ll[2]) ])

        # convert data into numpy matrix
        self.cellvec = np.array( cellvec )
        self.scale  = np.array( self.scale  )
        self.ntypes = np.array( self.ntypes )
        self.rrel    = np.array( self.rrel    )
        self.update_xyzrel_from_rrel()
        self.update_xyzcar_from_xyzrel()
        self.update_rcar_from_rrel()

    def read_rcar_positions( self, coordfile = "cartesian_coords", cellname = "cell"):
        '''
        read in cartesian coords
        '''
        self.cellvec = np.loadtxt(cellname)
        self.rcar = np.loadtxt(coordfile)
        self.update_rrel_from_rcar()
        self.update_xyzcar_from_rcar()
        self.update_xyzrel_from_xyzcar()
        return self.rcar

    def read_rrel_positions( self, coordfile = "EqCoords_direct", cellname = "cell"):
        '''
        read in relative coords
        '''
        self.cellvec = np.loadtxt(cellname)
        self.rrel = np.loadtxt(coordfile)
        self.update_rcar_from_rrel()
        self.update_xyzrel_from_rrel()
        self.update_xyzcar_from_xyzrel()
        return self.rcar

    def load_positions_cell(self, coordfile_cart = False, coord_cart = False, cellfile = False, cell = False, coordfile_rrel = False, coord_rrel = False):
        ####################################
        # cell first
        ####################################
        cellread = False
        if type(cell) != bool and type(cellfile) != bool:
            sys.exit("You need to provide either a cellfile or numpy vectors of the cell")
        if type(cell) != bool:
            self.cellvec = cell
            cellread = True
        if type(cellfile) != bool:
            self.cellvec = np.loadtxt(cellfile)
            cellread = True
        if cellread != True:
            sys.exit("could not load a cell")

        ####################################
        # positions after cell
        ####################################
        posread = False
        if type(coord_cart) != bool:
            self.rcar = coord_cart
            self.update_rrel_from_rcar()
            self.update_xyzcar_from_rcar()
            self.update_xyzrel_from_xyzcar()
            posread = True
            self.itype = np.ones( self.xcar.shape )
            self.stype = ['Al' for i in np.arange(self.xcar.size)]
        #print "11:",posread

        if type(coord_rrel) != bool:
            if posread == True:
                sys.exit("you got several pos")
            self.rrel = coord_rrel
            self.update_rcar_from_rrel()
            self.update_xyzrel_from_rrel()
            self.update_xyzcar_from_xyzrel()
            self.itype = np.ones( self.xcar.shape )
            self.stype = ['Al' for i in np.arange(self.xcar.size)]

            posread = True

        if type(coordfile_cart) != bool:
            if posread == True:
                sys.exit("you got several pos")
            if type(coordfile_cart) != str:
                sys.exit("coordfile_cart has to be a string")
            self.rcar = np.loadtxt(coordfile_cart)
            self.update_rrel_from_rcar()
            self.update_xyzcar_from_rcar()
            self.update_xyzrel_from_xyzcar()
            self.itype = np.ones( self.xcar.shape )
            self.stype = ['Al' for i in np.arange(self.xcar.size)]
            posread = True

        if type(coordfile_rrel) != bool:
            if posread == True:
                sys.exit("you got several pos")
            if type(coordfile_rrel) != str:
                sys.exit("coordfile_rrel has to be a string")
            self.rrel = np.loadtxt(coordfile_rrel)
            self.update_rcar_from_rrel()
            self.update_xyzrel_from_rrel()
            self.update_xyzcar_from_xyzrel()
            self.itype = np.ones( self.xcar.shape )
            self.stype = ['Al' for i in np.arange(self.xcar.size)]
            posread = True
        if posread == False:
            sys.exit("no positions read in")
        return


    def center_atoms_around_atom(self, atomnr, coordfile_cart = False, coord_cart = False, cellfile = False, cell = False, coordfile_rrel = False, coord_rrel = False):
        self.load_positions_cell(coordfile_cart = coordfile_cart, coord_cart = coord_cart, cellfile = cellfile, cell = cell, coordfile_rrel = coordfile_rrel, coord_rrel = coord_rrel)
        #print "---------------------"
        #print atomnr
        #print self.rcar
        #print "---------------------"
        self.translate_atoms_cart(self.rcar[atomnr])  # one example
        self.center_atoms_rel()
        return self.rcar

    def read_Hessematrix( self , hessematrixfile = "HesseMatrix_sphinx" ):
        self.hessematrix = np.loadtxt(hessematrixfile)*97.173617
        #You have:hartree/bohrradius^2
        #You want:eV/angstrom^2
        #    *97.173617
        #    /0.010290859
        return self.hessematrix

    def qh_forces( self, dpos = None, h = None ):
        '''
        dpos = deltas of positions in cartesian coords
        h = hessematrix
        '''
        f = np.dot(-h, np.ndarray.flatten(dpos))
        return f.reshape(dpos.shape)

    def qh_energy_cell( self, dpos = None, h = None ):
        f = self.qh_forces( dpos = dpos, h = h )
        e = -np.dot(np.ndarray.flatten(dpos),np.ndarray.flatten(f))/2.
        if e < 1e-8 and e > -1e-8:
            e = 0
        return e

    def qh_energy_atom( self, dpos = None, h = None ):
        ecell = self.qh_energy_cell( dpos = dpos, h = h )
        atoms = len(np.ndarray.flatten(dpos))/3.-1.
        e = e/atoms*1000.
        if e < 1e-8 and e > -1e-8:
            e = 0
        return e

    def get_NNlist(self, atom, NNschale, coordfile_cart = False, coord_cart = False, cellfile = False, cell = False, coordfile_rrel = False, coord_rrel = False, decimals = 6, verbose = False, return_NNdist = False, return_result_and_NNdist = False, return_d_NNidx=False):
        ''' get NN list of atom: atom with all NN included in NNschale
            e.g. atom = 1; [integer] center all atoms around 1st atom
            e.g. NNschale = 1; [integer] get all indizes of nearest neighbors
            e.g. NNschale = 2; [integer] get all indizes up to second nearest neighbor
        '''
        #def get_NNlist(self, atom, NNschale, coordfile_cart = False, cellfile = "cell", coordfile_rrel = False, decimals = 4, verbose = False,return_NNdist = False):
        #print "cc22:",coordfile_cart
        #print "cr22:",coordfile_rrel
        if False:
            print("------coordfile_cart------")
            print(coordfile_cart)
            print("------coord_cart------")
            print(coord_cart)
            print("------cellfile------")
            print(cellfile)
            print("------cell------")
            print(cell)
            print("------coordfile_rrel------")
            print(coordfile_rrel)
            print("------coord_rrel------")
            print(coord_rrel)

        self.load_positions_cell(coordfile_cart = coordfile_cart, coord_cart = coord_cart, cellfile = cellfile, cell = cell, coordfile_rrel = coordfile_rrel, coord_rrel = coord_rrel)
        self.center_atoms_around_atom(atom, coord_cart=self.rcar, cell=self.cellvec)
        #print "--"
        #print self.rcar
        #print "--"
        #d = np.zeros(self.rcar.shape[0])
        d_NNidx = np.zeros(self.rcar.shape[0])

        neighborlist = np.empty(self.rcar.shape[0])
        neighborlist[:] = np.nan

        #for ind,i in enumerate(self.rcar):
        #    d[ind] = np.linalg.norm(i)
        #    #print i, np.linalg.norm(i)
        d = np.linalg.norm(self.rcar,axis=1)

        def unique_to_decimals(liste,decimals):
            d = liste
            NNdist = np.unique(np.sort(np.around(d,decimals = decimals))) #NEVER accuracy is necessary
            NNdist4 = np.unique(np.sort(np.around(d,decimals = 4))) #NEVER accuracy is necessary
            decimalstake = 4
            for kk in np.arange(4,18):
                NNdistcheck = np.unique(np.sort(np.around(d,decimals = kk))) #NEVER accuracy is necessary
                #print kk,len(NNdistcheck),NNdistcheck
                if len(NNdistcheck) == len(NNdist4):
                    NNdist = NNdistcheck
                    decimalstake = kk
            if decimalstake < decimals:
                sys.exit("not enough digits to get precise next nearest neighbor distance")
            #print "decimalstake:",decimalstake
            return NNdist

        NNdist = unique_to_decimals(d,6)

        #NNdist = np.unique(np.sort(np.around(d,decimals = decimals))) #NEVER accuracy is necessary

        #NNdist4 = np.unique(np.sort(np.around(d,decimals = 4))) #NEVER accuracy is necessary
        #decimalstake = 4
        #for kk in np.arange(4,18):
        #    NNdistcheck = np.unique(np.sort(np.around(d,decimals = kk))) #NEVER accuracy is necessary
        #    #print kk,len(NNdistcheck),NNdistcheck
        #    if len(NNdistcheck) == len(NNdist4):
        #        NNdist = NNdistcheck
        #        decimalstake = kk
        #if decimalstake < 6:
        #    sys.exit("not enough digits to get precise next nearest neighbor distance")

        #print "-----------"
        #print "decimalstake:",decimalstake,NNdist
        #print "-----------"
        #print d
        #########################################
        #print "NNdist:",NNdist
        #NNdist: [  0.         2.821356   3.99       4.886732   5.642712   6.308744
        #           6.910883   7.464606   7.98       8.464068   8.921911   9.357379
        #              9.773464  10.172544  11.285424  11.632749  11.97      13.821765]

        #########################################
        if return_NNdist == True:
            return NNdist

        if verbose == True:
            print("NNdist:(vorhanden)",NNdist)

        def find_nearest_value(array,value):
            idx = (np.abs(array-value)).argmin()
            return array[idx]

        def find_nearest_index(array,value):
            idx = (np.abs(array-value)).argmin()
            return idx

        for ind,i in enumerate(d):
             d_NNidx[ind] = find_nearest_index(NNdist,i)

        if verbose == True:
            print("d_NNidx:",d_NNidx)
        #########################################
        #print "d_NNidx:",d_NNidx
        #########################################
        if return_d_NNidx == True:
            for idx,i in enumerate(d_NNidx):
                d_NNidx[idx] = int(i)
            return d_NNidx
        #result_and_0_all = np.nonzero(d_NNidx <= NNschale)[0]   # including smaller shells
        #print "alll",result_and_0_all
        #if type(NNschale) ==
        if type(NNschale) == int:
            result_and_0 = np.nonzero(d_NNidx == NNschale)[0]    # only certain shell
            #print "part",result_and_0
            all_no_0 = np.nonzero(d_NNidx > 0)[0]

            if return_result_and_NNdist == True:
                return np.intersect1d(result_and_0,all_no_0),NNdist
            else:
                return np.intersect1d(result_and_0,all_no_0)
        if type(NNschale) == list:
            out = []
            for i in NNschale:
                result_and_0 = np.nonzero(d_NNidx == i)[0]    # only certain shell
                #print "part",result_and_0
                all_no_0 = np.nonzero(d_NNidx > 0)[0]
                add = np.intersect1d(result_and_0,all_no_0)
                #print "add:",add
                out.append(add)
                #print "out:",out

            if return_result_and_NNdist == True:
                return out,NNdist
            else:
                return out

class unit_cell( crystal ):
    '''
        Creates a unit cell in the crystal format.
        Predefined basis:
            cubic bcc
            cubic bcc with octahedral sites
            cubic bcc with octahedral sites (unique labels)
            cubic fcc
            cubic fcc with octahedral sites
    '''
    def __init__( self ):
        super(unit_cell,self).__init__()

    def load_sc_cell( self, aParam ):
        '''
            generates a conventional sc sell with lattice parameter aParam
        '''
        self.crystaltype = 'sc'
        # clear variabls
        self.rrel = []
        #
        # define lattice vectors
        self.scale        = 1.0
        self.cellvec[0,:] = aParam * np.array([1,0,0])
        self.cellvec[1,:] = aParam * np.array([0,1,0])
        self.cellvec[2,:] = aParam * np.array([0,0,1])
        #
        # define basis in relative coordinates
        np.append( self.rrel, [0.0,0.0,0.0] )
        self.rrel.append([0.0,0.0,0.0])
        self.itype.append(1)
        self.stype.append('Atype')
        #
        # convert lists into numpy arrays
        self.rrel = np.array( self.rrel )
        self.itype = np.array( self.itype )
        #
        # updates other variabls in the class
        self.update_xyzrel_from_rrel()
        self.update_xyzcar_from_xyzrel()
        self.updateCellDetails()

    def load_bcc_cell( self, aParam ):
        '''
            generates a conventional bcc sell with lattice parameter aParam
        '''
        self.crystaltype = 'bcc'
        # clear variabls
        self.rrel = []
        #
        # define lattice vectors
        self.scale        = 1.0
        self.cellvec[0,:] = aParam * np.array([1,0,0])
        self.cellvec[1,:] = aParam * np.array([0,1,0])
        self.cellvec[2,:] = aParam * np.array([0,0,1])
        #
        # define basis in relative coordinates
        np.append( self.rrel, [0.0,0.0,0.0] )
        self.rrel.append([0.0,0.0,0.0])
        self.itype.append(1)
        self.stype.append('Atype')
        self.rrel.append([0.5,0.5,0.5])
        self.itype.append(1)
        self.stype.append('Atype')
        #
        # convert lists into numpy arrays
        self.rrel = np.array( self.rrel )
        self.itype = np.array( self.itype )
        #
        # updates other variabls in the class
        self.update_xyzrel_from_rrel()
        self.update_xyzcar_from_xyzrel()
        self.updateCellDetails()

    def load_bcc_cell_woct( self, aParam ):
        '''
            Generates a conventional bcc sell with octahedtral sites and
            lattice parameter aParam.
        '''
        self.crystaltype = 'bcc'
        # clear variabls
        self.rrel = []
        #
        # define lattice vectors
        self.cellvec[0,:] = aParam * np.array([1,0,0])
        self.cellvec[1,:] = aParam * np.array([0,1,0])
        self.cellvec[2,:] = aParam * np.array([0,0,1])
        #
        # define basis in relative coordinates
        self.rrel.append([0.0,0.0,0.0])
        self.itype.append(1)
        self.stype.append('Atype')
        #
        self.rrel.append([0.5,0.5,0.5])
        self.itype.append(1)
        self.stype.append('Atype')
        #
        self.rrel.append([0.5,0.0,0.0])
        self.itype.append(4)
        self.stype.append('oct1')
        #
        self.rrel.append([0.0,0.5,0.5])
        self.itype.append(4)
        self.stype.append('oct1')
        #
        self.rrel.append([0.0,0.5,0.0])
        self.itype.append(4)
        self.stype.append('oct2')
        #
        self.rrel.append([0.5,0.0,0.5])
        self.itype.append(4)
        self.stype.append('oct2')
        #
        self.rrel.append([0.0,0.0,0.5])
        self.itype.append(4)
        self.stype.append('oct3')
        #
        self.rrel.append([0.5,0.5,0.0])
        self.itype.append(4)
        self.stype.append('oct3')
        #
        # convert lists into numpy arrays
        self.rrel = np.array(  self.rrel  )
        self.itype = np.array( self.itype )
        #
        # updates other variabls in the class
        self.update_xyzrel_from_rrel()
        self.update_xyzcar_from_xyzrel()
        self.updateCellDetails()

    def load_bcc_cell_woct2( self, aParam ):
        '''
            Generates a conventional bcc sell with octahedtral sites and
            lattice parameter aParam. Octahedral sites are labeled as follows:
            4 - x-aligned sites
            5 - y-aligned sites
            6 - z-aligned sites
        '''
        self.crystaltype = 'bcc'
        # clear variabls
        self.rrel = []
        #
        # define lattice vectors
        self.cellvec[0,:] = aParam * np.array([1,0,0])
        self.cellvec[1,:] = aParam * np.array([0,1,0])
        self.cellvec[2,:] = aParam * np.array([0,0,1])
        #
        # define basis in relative coordinates
        self.rrel.append([0.0,0.0,0.0])
        self.itype.append(1)
        self.stype.append('Atype')
        #
        self.rrel.append([0.5,0.5,0.5])
        self.itype.append(1)
        self.stype.append('Atype')
        #
        self.rrel.append([0.5,0.0,0.0])
        self.itype.append(4)
        self.stype.append('oct1')
        #
        self.rrel.append([0.0,0.5,0.5])
        self.itype.append(4)
        self.stype.append('oct1')
        #
        self.rrel.append([0.0,0.5,0.0])
        self.itype.append(5)
        self.stype.append('oct2')
        #
        self.rrel.append([0.5,0.0,0.5])
        self.itype.append(5)
        self.stype.append('oct2')
        #
        self.rrel.append([0.0,0.0,0.5])
        self.itype.append(6)
        self.stype.append('oct3')
        #
        self.rrel.append([0.5,0.5,0.0])
        self.itype.append(6)
        self.stype.append('oct3')
        #
        # convert lists into numpy arrays
        self.rrel = np.array(  self.rrel  )
        self.itype = np.array( self.itype )
        #
        # updates other variabls in the class
        self.update_xyzrel_from_rrel()
        self.update_xyzcar_from_xyzrel()
        self.updateCellDetails()

    def load_bcc_cell_wtet( self, aParam ):
        '''
            Generates a conventional bcc sell with tetrahedral sites and
            lattice parameter aParam. Octahedral sites are labeled as follows:
            4 - x-aligned sites
            5 - y-aligned sites
            6 - z-aligned sites
        '''
        self.crystaltype = 'bcc'
        # clear variabls
        self.rrel   = []
        self.itype2 = []
        #
        # define lattice vectors
        self.cellvec[0,:] = aParam * np.array([1,0,0])
        self.cellvec[1,:] = aParam * np.array([0,1,0])
        self.cellvec[2,:] = aParam * np.array([0,0,1])
        #
        # define basis in relative coordinates
        self.rrel.append([0.00,0.00,0.00])
        self.itype2.append(1)
        self.itype.append(1)
        self.stype.append('Atype')
        self.rrel.append([0.50,0.50,0.50])
        self.itype2.append(1)
        self.itype.append(1)
        self.stype.append('Atype')
        #
        self.rrel.append([0.00,0.50,0.25])
        self.itype.append(4)
        self.itype2.append(4)
        self.stype.append('tet1')
        self.rrel.append([0.00,0.50,0.75])
        self.itype.append(4)
        self.itype2.append(4)
        self.stype.append('tet1')
        self.rrel.append([0.00,0.25,0.50])
        self.itype.append(4)
        self.itype2.append(4)
        self.stype.append('tet1')
        self.rrel.append([0.00,0.75,0.50])
        self.itype.append(4)
        self.itype2.append(4)
        self.stype.append('tet1')
        #
        self.rrel.append([0.50,0.00,0.25])
        self.itype2.append(5)
        self.itype.append(4)
        self.stype.append('tet2')
        self.rrel.append([0.50,0.00,0.75])
        self.itype2.append(5)
        self.itype.append(4)
        self.stype.append('tet2')
        self.rrel.append([0.25,0.00,0.50])
        self.itype2.append(5)
        self.itype.append(4)
        self.stype.append('tet2')
        self.rrel.append([0.75,0.00,0.50])
        self.itype2.append(5)
        self.itype.append(4)
        self.stype.append('tet2')
        #
        self.rrel.append([0.50,0.25,0.00])
        self.itype2.append(6)
        self.itype.append(4)
        self.stype.append('tet3')
        self.rrel.append([0.50,0.75,0.00])
        self.itype2.append(6)
        self.itype.append(4)
        self.stype.append('tet3')
        self.rrel.append([0.25,0.50,0.00])
        self.itype2.append(6)
        self.itype.append(4)
        self.stype.append('tet3')
        self.rrel.append([0.75,0.50,0.00])
        self.itype2.append(6)
        self.itype.append(4)
        self.stype.append('tet3')
        #
        # convert lists into numpy arrays
        self.rrel   = np.array( self.rrel   )
        self.itype  = np.array( self.itype  )
        self.itype2 = np.array( self.itype2 )
        #
        # updates other variabls in the class
        self.update_xyzrel_from_rrel()
        self.update_xyzcar_from_xyzrel()
        self.updateCellDetails()

    def load_fcc_cell( self, aParam ):
        '''
            Loads a conventional (cubic) fcc cell with
            lattice parameter aParam
        '''
        self.crystaltype = 'fcc'
        # clear variabls
        self.rrel = []
        #
        # define lattice vectors
        self.cellvec[0,:] = aParam * np.array([1,0,0])
        self.cellvec[1,:] = aParam * np.array([0,1,0])
        self.cellvec[2,:] = aParam * np.array([0,0,1])
        #
        # define basis in relative coordinates
        self.rrel.append([0.0,0.0,0.0])
        self.itype.append(1)
        self.stype.append('Atype')
        self.rrel.append([0.0,0.5,0.5])
        self.itype.append(1)
        self.stype.append('Atype')
        self.rrel.append([0.5,0.0,0.5])
        self.itype.append(1)
        self.stype.append('Atype')
        self.rrel.append([0.5,0.5,0.0])
        self.itype.append(1)
        self.stype.append('Atype')
        #
        # convert lists into numpy arrays
        self.rrel = np.array( self.rrel )
        self.itype = np.array( self.itype )
        #
        # updates other variabls in the class
        self.update_xyzrel_from_rrel()
        self.update_xyzcar_from_xyzrel()
        self.updateCellDetails()

    def load_fcc_cell_woct( self, aParam ):
        '''
            Loads a conventional (cubic) fcc cell with
            lattice parameter aParam
        '''
        self.crystaltype = 'fcc'
        # clear variabls
        self.rrel = []
        #
        # define lattice vectors
        self.scale        = 1.0
        self.cellvec[0,:] = aParam * np.array([1,0,0])
        self.cellvec[1,:] = aParam * np.array([0,1,0])
        self.cellvec[2,:] = aParam * np.array([0,0,1])
        #
        # define basis in relative coordinates
        self.rrel.append([0.0,0.0,0.0])
        self.itype.append(1)
        self.stype.append('Atype')
        self.rrel.append([0.0,0.5,0.5])
        self.itype.append(1)
        self.stype.append('Atype')
        self.rrel.append([0.5,0.0,0.5])
        self.itype.append(1)
        self.stype.append('Atype')
        self.rrel.append([0.5,0.5,0.0])
        self.itype.append(1)
        self.stype.append('Atype')
        self.rrel.append([0.0,0.0,0.5])
        self.itype.append(4)
        self.stype.append('oct')
        self.rrel.append([0.0,0.5,0.0])
        self.itype.append(4)
        self.stype.append('oct')
        self.rrel.append([0.5,0.0,0.0])
        self.itype.append(4)
        self.stype.append('oct')
        self.rrel.append([0.5,0.5,0.5])
        self.itype.append(4)
        self.stype.append('oct')#
        # convert lists into numpy arrays
        self.rrel = np.array( self.rrel )
        self.itype = np.array( self.itype )
        #
        # updates other variabls in the class
        self.update_xyzrel_from_rrel()
        self.update_xyzcar_from_xyzrel()
        self.updateCellDetails()

    def loadHcpOrtho( self, aParam, cParam ):
        '''
            Loads an othagonal hcp cell with
            lattice parameters aParam and cParam
        '''
        self.crystaltype = 'hcp'
        # clear variabls
        self.rrel = []
        #
        # define lattice vectors
        self.scale        = 1.0
        self.cellvec[0,:] =              aParam * np.array([1,0,0])
        self.cellvec[1,:] = np.sqrt(3) * aParam * np.array([0,1,0])
        self.cellvec[2,:] =              cParam * np.array([0,0,1])
        #
        # define basis in relative coordinates
        self.rrel.append([0.0,0.0,0.0])
        self.itype.append(1)
        self.stype.append('Atype')
        self.rrel.append([0.5,0.5,0.0])
        self.itype.append(1)
        self.stype.append('Atype')
        self.rrel.append([0.0,1./3.,0.5])
        self.itype.append(1)
        self.stype.append('Atype')
        self.rrel.append([0.5,5./6.,0.5])
        self.itype.append(1)
        self.stype.append('Atype')
        #
        # convert lists into numpy arrays
        self.rrel = np.array( self.rrel )
        self.itype = np.array( self.itype )
        #
        # updates other variabls in the class
        self.update_xyzrel_from_rrel()
        self.update_xyzcar_from_xyzrel()
        self.updateCellDetails()

    def loadHcpOrtho_wtet( self, aParam, cParam ):
        '''
            Loads an othagonal hcp cell with
            lattice parameters aParam and cParam
        '''
        self.crystaltype = 'hcp'
        # clear variabls
        self.rrel = []
        #
        # define lattice vectors
        self.scale        = 1.0
        self.cellvec[0,:] =              aParam * np.array([1,0,0])
        self.cellvec[1,:] = np.sqrt(3) * aParam * np.array([0,1,0])
        self.cellvec[2,:] =              cParam * np.array([0,0,1])
        #
        # define basis in relative coordinates
        self.rrel.append([ 0./2., 0./2., 0./2.])
        self.itype.append(1)
        self.stype.append('host')
        self.rrel.append([ 1./2., 1./2., 0./2.])
        self.itype.append(1)
        self.stype.append('host')
        self.rrel.append([ 0./3., 1./3., 1./2.])
        self.itype.append(1)
        self.stype.append('host')
        self.rrel.append([ 1./2., 5./6., 1./2.])
        self.itype.append(1)
        self.stype.append('host')

        self.rrel.append([ 0./2., 1./3., 1./6.])
        self.itype.append(3)
        self.stype.append('tetra1')
        self.rrel.append([ 1./2., 5./6., 1./6.])
        self.itype.append(3)
        self.stype.append('tetra1')
        self.rrel.append([ 0./2., 1./3., 5./6.])
        self.itype.append(3)
        self.stype.append('tetra1')
        self.rrel.append([ 1./2., 5./6., 5./6.])
        self.itype.append(3)
        self.stype.append('tetra1')

        self.rrel.append([ 0./2., 0./2., 2./6.])
        self.itype.append(4)
        self.stype.append('tetra2')
        self.rrel.append([ 1./2., 1./2., 2./6.])
        self.itype.append(4)
        self.stype.append('tetra2')
        self.rrel.append([ 0./2., 0./2., 4./6.])
        self.itype.append(4)
        self.stype.append('tetra2')
        self.rrel.append([ 1./2., 1./2., 4./6.])
        self.itype.append(4)
        self.stype.append('tetra2')
        #
        # convert lists into numpy arrays
        self.rrel = np.array( self.rrel )
        self.itype = np.array( self.itype )
        #
        # updates other variabls in the class
        self.update_xyzrel_from_rrel()
        self.update_xyzcar_from_xyzrel()
        self.updateCellDetails()

class remove_mapping_into_originalcell( crystal ):
    ''' before the whole cell can be repeated, the mapping to the original cell
        has to be remove (in other words: the atoms have to be close to their origial
        undisplaced positions, otherwise crystal_generator.center_atoms_around_atom
        does not work properly '''
    def __init__( self ):
        super(remove_mapping_into_originalcell,self).__init__()

    def remove_mapping( self, cryst, cryst0 ):
        '''
        '''

        cryst1 = copy.deepcopy(cryst)  # doing this would not chenage the instance of cryst input
        u = cryst1.rcar-cryst0.rcar
        for ind,i in enumerate(u):
            if  u[ind][0] > cryst1.cellvec[0,0]/2.:
                u[ind][0] = u[ind][0] - cryst1.cellvec[0,0]
            if  u[ind][1] > cryst1.cellvec[0,0]/2.:
                u[ind][1] = u[ind][1] - cryst1.cellvec[0,0]
            if  u[ind][2] > cryst1.cellvec[0,0]/2.:
                u[ind][2] = u[ind][2] - cryst1.cellvec[0,0]
        self.rcar = u + cryst0.rcar
        self.cellvec = cryst.cellvec
        self.update_rrel_from_rcar()
        self.update_xyzcar_from_rcar()
        self.update_xyzrel_from_xyzcar()
        return


class supercell( crystal ):
    '''
        Creates a supercell by tranlating a unit cell
    '''
    def __init__( self ):
        super(supercell,self).__init__()

    def create_supercell( self, unit, N1, N2, N3, newsorting = False, shiftforces = False ):
        '''
            This function creates a supercell by repeating a unit cell
            by Nx times in the x-direction, Ny times in the y-direction
            and Nz times in the z-direction
        '''
        self.crystaltype = unit.crystaltype
        # define the supercell vector
        self.cellvec[0,:] = unit.cellvec[0,:] * N1
        self.cellvec[1,:] = unit.cellvec[1,:] * N2
        self.cellvec[2,:] = unit.cellvec[2,:] * N3
        #
        # create translation matrix
        var1 = np.mgrid[ 0:N1, 0:N2, 0:N3 ]
        ntrans1 = var1[0,:,:].flatten(1)
        ntrans2 = var1[1,:,:].flatten(1)
        ntrans3 = var1[2,:,:].flatten(1)
        transmat = np.zeros( [ len(ntrans1),3] )
        transmat[:,0] = ntrans1
        transmat[:,1] = ntrans2
        transmat[:,2] = ntrans3
        #
        # calculate supercell lattice positions
        self.sclattpos = np.dot( transmat, unit.cellvec )

        if shiftforces != False and shiftforces != True:
            sys.exit("shiftforces has to be False or True")
        if shiftforces == True:
            #print "in---------"
            self.sclattpos = np.zeros(self.sclattpos.shape)
        #print "shiftforces:",shiftforces
        #print ">>self.sclattpos:",self.sclattpos

        #
        # intialize self.rcar and self.itype variable
        self.rcar = []

        # calculate atomic positions for basis
        #print sclattpos
        #print "----------"
        #for ij, sclatt in enumerate(sclattpos):
        #    for ii, current_ibasis in enumerate(unit.xcar):
        #        print sclatt[ij]+unit.rcar[ii]

        #        #if len(self.rcar)==0:
        #        #    self.rcar  = sclatt[ij] + unit.rcar[ii]
        #        #    self.itype = np.ones(len(sclattpos) ) * unit.itype[ii]
        #        #    self.stype = np.shape( sclattpos )[0] * [unit.stype[ii]]
        #        #else:
        #        #    self.rcar  = np.append( self.rcar, sclatt[ij] + unit.rcar[ii], axis=0 )
        #        #    self.itype = np.append( self.itype, np.ones(len(sclattpos) ) * \
        #        #            unit.itype[ii], axis=0 )
        #        #    self.stype = self.stype + ( np.shape( sclattpos )[0] * [unit.stype[ii]] )

        ##############################################################
        # gerards sorting
        for ii, current_ibasis in enumerate(unit.xcar):
            if len(self.rcar)==0:
                self.rcar  = self.sclattpos + unit.rcar[ii,:]
                self.itype = np.ones(len(self.sclattpos) ) * unit.itype[ii]
                self.stype = np.shape( self.sclattpos )[0] * [unit.stype[ii]]
            else:
                self.rcar  = np.append( self.rcar, self.sclattpos + unit.rcar[ii,:], axis=0 )
                self.itype = np.append( self.itype, np.ones(len(self.sclattpos) ) * \
                        unit.itype[ii], axis=0 )
                self.stype = self.stype + ( np.shape( self.sclattpos )[0] * [unit.stype[ii]] )


        # update other coordinates
        self.update_xyzcar_from_rcar()
        self.update_xyzrel_from_xyzcar()
        self.update_rrel_from_rcar()
        self.updateCellDetails()

        # new sorting if wished () ()
        if newsorting != False and newsorting != True:
            sys.exit("newsorting has to be False or True")
        if newsorting == True:
            if N1 != N2 or N2 != N3:
                sys.exit("N1 == N2 == N3 necessary for creating supercess with !!newsorting!!")
            self.crystalsc=np.zeros(self.rcar.shape)
            for ii in range( N1**3 ):  # ii: 0 ... 7  wenn nsc == 2; 0 .. 26 wenn nsc == 3
                for jj in range( unit.rcar.shape[0] ):  # jj:  1 .. 32
                    self.crystalsc[unit.rcar.shape[0]*ii+jj] = self.rcar[jj*(N1**3)+ii]
            #print "unit.rcar[0]",unit.rcar[0]
            #print "self.crystalsc[32]:",self.crystalsc[32]
            #print ""
            #print "self.sclattpos:",self.sclattpos
            #print "unit.rcar[1]",unit.rcar[1]
            #print "self.crystalsc[33]:",self.crystalsc[33]
            #print "unit.f[0][1]",unit.rcar[0][1]
            #print "unit.P[0][1]",unit.rcar[0][1]
            #print "unit.P[0][1]",unit.rcar[0][1]

            self.rcar = self.crystalsc
            self.update_xyzcar_from_rcar()
            self.update_xyzrel_from_xyzcar()
            self.update_rrel_from_rcar()
            self.updateCellDetails()

class dislocation( supercell ):
    '''
        Generates a dislocation core structure with line direction t, slip direction m
        and glide plane normal n
    '''
    def __init__( self ):
        super( dislocation ,self ).__init__()
        #
        # supercell orientation
        self.m = np.array([1,0,0]) # glide direction
        self.n = np.array([0,1,0]) # glide plane normal
        self.t = np.array([0,0,1]) # line direction
        #
        # dislocation positons and burgers vector
        self.rcarDisl = []
        self.bDisl    = []
        #
        # material parameters
        self.aParam = []
        self.cParam = []
        self.nu     = []
        #
        # Rotation matrix based on m, n and t
        self.R = []

    def getRotationMatrix( self ):
        '''
            Generates a bcc Fe dislocation
        '''
        Rinv = np.zeros([3,3])
        Rinv[0,:] = self.m/np.linalg.norm(self.m)
        Rinv[1,:] = self.n/np.linalg.norm(self.n)
        Rinv[2,:] = self.t/np.linalg.norm(self.t)
        self.R = np.linalg.inv( Rinv )

    def determine_tperiodicity( self, unit ):
        '''
            Determines the periodicity in the t-direction.
        '''
        unit = unit_cell()
        if   ( self.crystaltype == 'bcc' ):
            unit.load_sc_cell( self.aParam )
        elif ( self.crystaltype == 'fcc' ):
            unit.load_sc_cell( self.aParam )

        sctemp = supercell()
        sctemp.create_supercell( unit, 10, 10, 10 )

    def loadBccFeEdge( self ):
        '''
            Loads a Fe bcc edge dislocation
        '''
        # define material parameters
        self.aParam = 2.855315
        self.nu     = 0.29

        # define system coordinates
        self.m = np.array([  1.0,  1.0,  1.0 ])
        self.n = np.array([  1.0, -1.0,  0.0 ])
        self.t = np.cross( self.m, self.n )
        self.getRotationMatrix()

    # define unit cell
    # TODO: this routine is unfinished!


if __name__ == '__main__':
    import sys
    crystal = crystal()
    print(crystal.get_NNlist(0,1,coordfile_cart = False, coordfile_rrel = "EqCoords_direct", cellfile = "cell"))  # atomindex = 15, up_to_schale = 2
    print("-----------")
    print(crystal.get_NNlist(0,1,coordfile_cart = "pos0", coordfile_rrel = False, cellfile = "cell"))  # atomindex = 15, up_to_schale = 2
    sys.exit()

    aParam = 2.855315
    bccunit = unit_cell()
    bccunit.load_bcc_cell(aParam)
    bccsc = supercell()
    bccsc.create_supercell( bccunit, 11, 11, 11 )
    bccsc.translate_atoms_cart( aParam*np.array([ 5, 5, 5]) )
    m = np.array([ 1, 1, 1 ])
    n = np.array([ 1,-1, 0 ])
    m = m/np.linalg.norm(m)
    n = n/np.linalg.norm(n)
    t = np.cross( m, n )
    Rinv = np.zeros([3,3])
    Rinv[0,:] = m
    Rinv[1,:] = n
    Rinv[2,:] = t
    R = np.linalg.inv(Rinv)
    bccsc.rotate_atoms(R)
    bccsc.updateCellDetails()
