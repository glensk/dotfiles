import os
import sys
import numpy             as np
import matplotlib.pyplot as plt
import matplotlib        as mpl
import scipy             as sp
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.ticker    import LinearLocator, FormatStrFormatter
from scipy                import optimize
from scipy.interpolate    import UnivariateSpline
from scipy.interpolate    import InterpolatedUnivariateSpline

class splineStuff( object ):
    def __init__( self ):
        self.x       = None
        self.y       = None
        self.indir   = 'data'
        self.infname = None

        return

    def loadData( self ):
        z5 = np.loadtxt("/Users/glensk/Dropbox/Understand_distributions/z_0.5")
        z5 = np.loadtxt("/Users/glensk/Dropbox/Understand_distributions/z_0.5")
        z5 = np.loadtxt("/Users/glensk/Dropbox/Understand_distributions/z_0.0")
        #data = np.loadtxt( self.indir + "/" + self.infname )
        data = z5[np.nonzero(z5[:,1] == 0.0)][:,[0,6]]
        print data
        u, idu = np.unique(data[:,0], return_index=True)
        self.x = data[idu,0]
        self.y = data[idu,1]

        return

class splineStuff_visual( splineStuff ):
    def __init__( self ):
        self.x       = None
        self.y       = None
        self.indir   = 'data'
        self.infanme = None

        self.fignum = 1
        self.figsize = [5.5,5.]

        # matplotlib default settings
        self.usetex = True
        self.family = 'serif'
        self.fc     = 'c'
        self.lw     = 1.0
        self.mew    = 1.0
        self.ms     =  5
        self.axeslw = 1.5

        return

    def initializeMpl( self ):
        '''
            Initializes matplotlib with latex, etc
        '''
        # generate plots
        plt.rc('text' , usetex     = self.usetex )
        plt.rc('font' , family     = self.family )
        mpl.rc('patch', facecolor  = self.fc     )
        mpl.rc('lines', linewidth  = self.lw     )
        mpl.rc('lines', mew        = self.mew    )
        mpl.rc('lines', markersize = self.ms     )
        mpl.rc('axes' , linewidth  = self.axeslw )
        plt.ion();

        return

    def plotForces( self ):

        a0 = 3.99
        b  = a0 / np.sqrt(2.)
        delta = 1.0

        ''' cubic spline stuff  '''
        ii    = 1
        #spl   = UnivariateSpline( self.x[::ii], self.y[::ii], k=2, s=0.00000000000*b  )
        print "y->",type(self.x),self.x
        print "x->",type(self.y),self.y
        spl   = UnivariateSpline( self.x, self.y, k=2, s=0.0  )
        print "cc:",spl.get_coeffs()
        ispl  = spl.antiderivative()
        spl2  = ispl.derivative()

        plt.ion();
        fig = plt.figure( self.fignum )
        fig.set_size_inches( self.figsize )
        fig.clf()
        ax = fig.add_subplot(211)
        ax.plot( self.x, self.y, 'kx' )
        ax.plot( self.x[::ii], self.y[::ii], 'rx' )
        ax.plot( [self.x.min()-0.3,self.x.max()+0.3], np.zeros(2), '--r' )
        #ax.set_xlim([b-delta, b+delta])
        #ax.set_ylim([-0.01,0.01])
        ax.axis('tight')

        xx      = np.linspace( self.x.min(), self.x.max(), 1001 )
        self.x2 = xx
        self.y2 = spl2( xx )
        ax.plot( self.x2, self.y2, '-b' )

        xx1      = np.linspace( self.x.min()-0.3, self.x.max()+0.3, 1001 )
        self.xx2 = xx1
        self.yy2 = spl2( xx1 )
        ax.plot( self.xx2, self.yy2, '--r' )
        ymin, ymax = spl(self.x[::ii]+0.3).min(),spl(self.x[::ii]+0.3).max()
        print "ymin,ymax",ymin,ymax
        ax.set_ylim( [ymin-.3, ymax+0.3] )

        ax = fig.add_subplot(212)
        #ax.plot( self.x2, ispl(self.x2), '-b' )
        ax.plot( self.x[::ii], self.y[::ii] - spl(self.x[::ii]), 'kx')
        #ax.set_xlim([b-delta, b+delta])
        ax.axis('tight')

        return spl


if __name__ == '__main__':
    data = splineStuff_visual()
    data.indir   = 'data'
    data.infname = 'data2.dat'
    data.loadData()

    data.initializeMpl()
    spl = data.plotForces()
