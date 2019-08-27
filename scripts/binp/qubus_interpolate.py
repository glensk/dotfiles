#!/usr/bin/env python

##############################################################################
# for the fitting it is best to start from z=0 plane and find the integration
# constant for the Energy fit
##############################################################################

# np.savetxt("z_0.0",pot.disp.qubus.data[(pot.disp.qubus.data.z == 0.0)][['x','y','z','Fx','Fy','Fz','fmlx','fmly','fmlz']].values)
import os
import sys
import utils
import pickle
import numpy as np

from scipy                import optimize
from scipy.interpolate    import UnivariateSpline
from scipy.interpolate    import RectBivariateSpline
import matplotlib.pyplot as plt

class qubus_splines(object):
    ''' all splines are as functioni of x
    TODO:
            1. get xy spline for all 26 planes (seperately for Fx, Fy, and Fz, and E)
            2. fix x,y point and get Fx(z) for every of 26 planes
            3. make a spline of 2 to get Fx Fy and Fz

            1a. read in qubus pkl and save all fml data sp
    '''
    def __init__(self):
        self.verbose = False

        ################# for plotting
        self.fignum = 1
        self.figsize = [5.5,5.]

    def loadQubuspd(self,filename):
        ''' loads pandas dataframe '''
        if type(filename) != str:
            print "filename:",filename
            sys.exit("loadDatapd  needs as input a picklefile (str) but got "+str(filename)+" of type "+str(type(filename)))
        if os.path.isfile(filename) != True:
            sys.exit("file "+filename+" does not exist!")
        #self.data = pd.read_pickle(filename)
        print utils.printgreen("     loading pkl : "+str(filename))
        with open(filename) as f:
            self.data = pickle.load(f)
            self.get_qubus_info()
        return

    def get_qubus_info(self):
        self.xyplane = self.data[(self.data.z == 0.0)]\
                [['x','y','z','fmlx','fmly','fmlz']].values
        self.x = np.unique(self.xyplane[:,0])
        self.y = np.unique(self.xyplane[:,1])
        self.z = self.x
        return


    def spline_for_z_y_plot(self, z = False, y = False):
        self.spline_for_z_y(z = z, y = y)

        plt.ion();
        fig = plt.figure( self.fignum )
        fig.set_size_inches( self.figsize )
        fig.clf()

        # measured points
        ax = fig.add_subplot(211)
        ax.plot( self.x_all, self.fmlx_all, 'kx' )
        ax.axis('tight')

        # zero line
        ax.plot( [self.x_all.min()-0.3,self.x_all.max()+0.3], np.zeros(2), '--r' )

        # spline line
        self._xs = np.linspace( self.x_all.min()-0.3, self.x_all.max()+0.3, 1001 )
        self._ys = self.splfmlx_f( self._xs )
        ax.plot( self._xs, self._ys, '-b' )

        # xmin,xmax,ymin,ymax
        ymin = self._ys.min()
        ymax = self._ys.max()
        xmin = self._xs.min()
        xmax = self._xs.max()
        print "xmin,xmax",xmin,xmax
        print "ymin,ymax",ymin,ymax
        ax.set_xlim( [xmin, xmax] )
        ax.set_ylim( [ymin, ymax] )

        #self._      = np.linspace( self.x.min()-0.3, self.x.max()+0.3, 1001 )
        #self.xx2 = xx1
        #self.yy2 = spl2( xx1 )
        #ax.plot( self.xx2, self.yy2, '--r' )

        ax = fig.add_subplot(212)
        #ax.plot( self.x2, ispl(self.x2), '-b' )
        ax.plot( self.x_all, self.fmlx_all - self.splfmlx_f(self.x_all), 'kx')
        #ax.set_xlim([b-delta, b+delta])
        ax.axis('tight')

        return


    def spline_for_z_y(self, z = 0.0, y = 0.0):
        ''' '''
        self.plane = self.data[(self.data.z == z) & (self.data.y == y)]\
                [['x','y','z','Fx','Fy','Fz','fmlx','fmly','fmlz']].values

        self.x_all = self.plane[:,0]

        self.Fx_all = self.plane[:,3]
        self.Fy_all = self.plane[:,4]
        self.Fz_all = self.plane[:,5]

        self.fmlx_all = self.plane[:,6]
        self.fmly_all = self.plane[:,7]
        self.fmlz_all = self.plane[:,8]

        # this is only necessary in case .pkl file was created as object type
        #self.x_all = np.array(self.plane[:,0], dtype='float')
        #self.fmlx_all = np.array(self.plane[:,3],dtype='float')
        #self.fmly_all = np.array(self.plane[:,4],dtype='float')
        #self.fmlz_all = np.array(self.plane[:,5],dtype='float')

        # forces parametrization
        self.splfmlx_f   = UnivariateSpline( self.x_all, self.fmlx_all, k=2, s=0.0)
        self.splfmly_f   = UnivariateSpline( self.x_all, self.fmly_all, k=2, s=0.0)
        self.splfmlz_f   = UnivariateSpline( self.x_all, self.fmlz_all, k=2, s=0.0)
        self.splFx          = UnivariateSpline( self.x_all, self.Fx_all, k=2, s=0.0)
        self.splFy          = UnivariateSpline( self.x_all, self.Fy_all, k=2, s=0.0)
        self.splFz          = UnivariateSpline( self.x_all, self.Fz_all, k=2, s=0.0)

        #self.splfmlx_e = - self.splfmlx_achsenabschnitt + self.splfmlx_f.antiderivative()
        return

    def get_achsenabschnitt_at_0(self):
        self.spline_for_z_y(z = 0.0, y = 0.0)
        # for energy parametrization
        self.splfmlx_achsenabschnitt  = self.splfmlx_f.antiderivative()(0)
        self.splfmly_achsenabschnitt  = self.splfmly_f.antiderivative()(0)
        self.splfmlz_achsenabschnitt  = self.splfmlz_f.antiderivative()(0)
        self.Fx_achsenabschnitt  = self.splFx.antiderivative()(0)
        self.Fy_achsenabschnitt  = self.splFy.antiderivative()(0)
        self.Fz_achsenabschnitt  = self.splFz.antiderivative()(0)

    def splfmlx_e(self, x):
        return self.splfmlx_f.antiderivative()(x) - self.splfmlx_achsenabschnitt
    def splfmly_e(self, x):
        return self.splfmly_f.antiderivative()(x) - self.splfmly_achsenabschnitt
    def splfmlz_e(self, x):
        return self.splfmlz_f.antiderivative()(x) - self.splfmlz_achsenabschnitt
    def splEx(self, x):
        return self.splFx.antiderivative()(x) - self.Fx_achsenabschnitt
    def splEy(self, x):
        return self.splFy.antiderivative()(x) - self.Fy_achsenabschnitt
    def splEz(self, x):
        return self.splFz.antiderivative()(x) - self.Fz_achsenabschnitt

    def ka(self):
        '''
        - die integration der energien sollte erst gemacht werden wenn der ganze F{x,y,z}
        wuerfel erstellt ist
        - zuerst sollte fuer die 3 hauptaxen (x,y,z) die integration gemacht werden
          dies entspricht: Fx[:,13,13], Fx[13,:,13], Fx[13,13,:]
          dies natuerlich fuer Fy und Fz wobei spaeter nachgeprueft werden sollte
          ob die integration das die gleiche energie liefert fuer jeden wuerfel.
        - die energien an den hauptachsen koennen jetzt als aufpunkte fuer die energie
          integration angesehen werden.
        - zum schluss muss immer fuer die energien alle symmetrien gelten:
          Ex[x,y,z=0] == Ey[y,x,z=0]
        - super pragmatisch waere einfach die VASP kraft fuer die entsprechende auslenkung
          zu nehmen. (ah, diese gilt natuerlich nur fuer die ganze zelle, also lieber nicht)
        - um dies ganz zu verstehen sollte man evtl. erstmal
          a) eine lineare kette studieren und von kraeften die energien integrieren
          b) eine flaeche studieren und von kraeften die energien integrieren
          c) dies dann auf den qubus uebertragen
        '''
        # check die sache mit den aufpunkten erstmal fuer die z = 0 flaeche und
        # a) mache fuer Fx :
        #           Ex bei y=0 sowie Ey bei x = 0, dies gibt die aufpunkte
        # b)        Ex bei y=1 und setze auf den aufpunkt bei (0/1)
        #           Ey bei x=1 und setze auf den aufpunkt bei (1/0)
        #           nun sollten die Energien symmetrisch sein sowie ben gleichen wert bei (1/1)
        #           --> somit auch bei (0.5/0.5)
        # b) gegenchecken wenn die aufpunkte entlang der quer richtung sind.
        # c) mach das gleiche fuer Fy sowie Fz:
        #           nun vergleichen ob (x=y=z, 0.5,0.5,0.5) die gleiche energie gibt fuer F{x,y,z}
        # a)
        self.fx0x = self.Fx[:,13,13]
        self.fx0y = self.Fx[13,:,13]
        self.fx0z = self.Fx[13,13,:]
        self.ex0x_aa = UnivariateSpline(self.x_all,self.fx0x,k=2,s=0.0).antiderivative()(0)
        self.ex0y_aa = UnivariateSpline(self.x_all,self.fx0y,k=2,s=0.0).antiderivative()(0)
        self.ex0z_aa = UnivariateSpline(self.x_all,self.fx0z,k=2,s=0.0).antiderivative()(0)
        self.ex0x = UnivariateSpline(self.x_all,self.fx0x,k=2,s=0.0).antiderivative()(self.x_all)
        self.ex0y = UnivariateSpline(self.x_all,self.fx0y,k=2,s=0.0).antiderivative()(self.x_all)
        self.ex0z = UnivariateSpline(self.x_all,self.fx0z,k=2,s=0.0).antiderivative()(self.x_all)
        self.ex0x = self.ex0x - self.ex0x_aa
        self.ex0y = self.ex0y - self.ex0y_aa
        self.ex0z = self.ex0z - self.ex0z_aa

        # b)  1 ist in diesemfall 0.5
        ap = 18 # 18 == 0.5
        self.fx18x =self.Fx[:,ap,13]
        self.fx18y =self.Fx[ap,:,13]
        self.fx18z =self.Fx[13,ap,:]   # ist ehe 0
        self.ex18x_aa = UnivariateSpline(self.x_all,self.fx18x,k=2,s=0.0).antiderivative()(0)
        self.ex18y_aa = UnivariateSpline(self.x_all,self.fx18y,k=2,s=0.0).antiderivative()(0)
        self.ex18z_aa = UnivariateSpline(self.x_all,self.fx18z,k=2,s=0.0).antiderivative()(0)
        self.ex18x=UnivariateSpline(self.x_all,self.fx18x,k=2,s=0.0).antiderivative()(self.x_all)
        self.ex18y=UnivariateSpline(self.x_all,self.fx18y,k=2,s=0.0).antiderivative()(self.x_all)
        self.ex18z=UnivariateSpline(self.x_all,self.fx18z,k=2,s=0.0).antiderivative()(self.x_all)
        self.ex18x = self.ex18x - self.ex18x_aa
        self.ex18y = self.ex18y - self.ex18y_aa
        self.ex18z = self.ex18z - self.ex18z_aa

        self.ex18x = self.ex18x + self.ex0y[ap]
        self.ex18y = self.ex18y + self.ex0x[ap]
        self.ex18z = self.ex18z + self.ex0z[ap]
        print "done:"
        tofolder="/Users/glensk/Dropbox/Understand_distributions/displacements_dense/Ir/2x2x2sc_qubus_3x3x3kp/3.99Ang_0.5_0.0_0.4/mit_spline_correction_long/tmp_check_energyies"
        np.savetxt(tofolder+"/ex0x",np.transpose([self.x,self.ex0x]))
        np.savetxt(tofolder+"/ex0y",np.transpose([self.x,self.ex0y]))
        np.savetxt(tofolder+"/ex0z",np.transpose([self.x,self.ex0z]))
        np.savetxt(tofolder+"/ex18x",np.transpose([self.x,self.ex18x]))
        np.savetxt(tofolder+"/ex18y",np.transpose([self.x,self.ex18y]))
        return

    def kb(self):
        xa = self.x_all
        self.Exn = np.zeros((self.x.size,self.y.size,self.z.size))
        self.Eyn = np.zeros((self.x.size,self.y.size,self.z.size))
        self.Ezn = np.zeros((self.x.size,self.y.size,self.z.size))


        for idz,zv in enumerate(self.z):
            for idy,yv in enumerate(self.y):
                self.ex_aa = UnivariateSpline(xa,self.Fx[:,idy,idz],k=2,s=0.0).antiderivative()(0)
                self.ex = UnivariateSpline(xa,self.Fx[:,idy,idz],k=2,s=0.0).antiderivative()(xa)
                self.ex = self.ex - self.ex_aa
                self.Exn[:,idy,idz] = self.ex


        for idz,zv in enumerate(self.z):
            for idx,xv in enumerate(self.y):
                self.ey_aa = UnivariateSpline(xa,self.Fy[idx,:,idz],k=2,s=0.0).antiderivative()(0)
                self.ey = UnivariateSpline(xa,self.Fy[idx,:,idz],k=2,s=0.0).antiderivative()(xa)
                self.ey = self.ey - self.ey_aa
                self.Eyn[idx,:,idz] = self.ey

        for idy,zv in enumerate(self.z):
            for idx,xv in enumerate(self.y):
                self.ez_aa = UnivariateSpline(xa,self.Fz[idx,idy,:],k=2,s=0.0).antiderivative()(0)
                self.ez = UnivariateSpline(xa,self.Fz[idx,idy,:],k=2,s=0.0).antiderivative()(xa)
                self.ez = self.ez - self.ez_aa
                self.Ezn[idx,idy,:] = self.ez
        self.Ex = self.Exn
        self.Ey = self.Eyn
        self.Ez = self.Ezn
        return


    def ef_in_whole_qubus(self,folder):
        self.spline_for_z_y(z = 0.0, y = 0.0)
        self.get_achsenabschnitt_at_0()
        print "self.splfmlx_achsenabschnitt:",self.splfmlx_achsenabschnitt
        print "self.splfmly_achsenabschnitt:",self.splfmly_achsenabschnitt
        print "self.splfmlz_achsenabschnitt:",self.splfmlz_achsenabschnitt
        print "self.Fx_achsenabschnitt:",self.Fx_achsenabschnitt
        print "self.Fy_achsenabschnitt:",self.Fy_achsenabschnitt
        print "self.Fz_achsenabschnitt:",self.Fz_achsenabschnitt
        print self.x
        print self.y
        print self.z
        self.fmlx = np.zeros((self.x.size,self.y.size,self.z.size))
        self.fmly = np.zeros((self.x.size,self.y.size,self.z.size))
        self.fmlz = np.zeros((self.x.size,self.y.size,self.z.size))
        self.emlx = np.zeros((self.x.size,self.y.size,self.z.size))
        self.emly = np.zeros((self.x.size,self.y.size,self.z.size))
        self.emlz = np.zeros((self.x.size,self.y.size,self.z.size))
        self.Fx = np.zeros((self.x.size,self.y.size,self.z.size))
        self.Fy = np.zeros((self.x.size,self.y.size,self.z.size))
        self.Fz = np.zeros((self.x.size,self.y.size,self.z.size))
        self.Ex = np.zeros((self.x.size,self.y.size,self.z.size))
        self.Ey = np.zeros((self.x.size,self.y.size,self.z.size))
        self.Ez = np.zeros((self.x.size,self.y.size,self.z.size))


        # es muss erst die ganze F{x,y,z} oberflaeche aufgebaut werden bevor mann E{x,y,z}
        # integriert, da sonst zwar fuer Ex alle daten vorliegen, jedoch nicht fuer Ey,Ez
        TESTRUN = False
        if TESTRUN == True:
            zv = 0; idz = 13
            for idy,yv in enumerate(self.y):
                print "z:",zv,"y:",yv,'idy:',idy
                for idx,xv in enumerate(self.x):
                    #print "z:",zv,"y:",yv,"x:",xv
                    self.spline_for_z_y(z = zv, y = yv)
                    self.Fx[idx,idy,idz] = self.splFx(xv)
                    self.Fy[idx,idy,idz] = self.splFy(xv)
                    self.Fz[idx,idy,idz] = self.splFz(xv)

            #zv = 0; idz = 13
            #for idy,yv in enumerate(self.y):
            #    print "z:",zv,"y:",yv,'idy:',idy
            #    for idx,xv in enumerate(self.x):
                    self.Ex[idx,idy,idz] = self.splEx(xv)
                    self.Ey[idx,idy,idz] = self.splEy(xv)
                    self.Ez[idx,idy,idz] = self.splEz(xv)
            return
            sys.exit('yo this was just the testrun')
            print "saving stuff:"
            for idy,yv in enumerate(self.y):
                print "z:",zv,"y:",yv,'idy:',idy,"-->"
                tofolder = "/Users/glensk/Dropbox/Understand_distributions/displacements_dense/Ir/2x2x2sc_qubus_3x3x3kp_analyze"
                np.savetxt(tofolder+"/Fx_"+str(yv),np.transpose([self.x,self.Fx[:,idy,13]]))
                np.savetxt(tofolder+"/Ex_"+str(yv),np.transpose([self.x,self.Ex[:,idy,13]]))

                np.savetxt(tofolder+"/Fxt_"+str(yv),np.transpose([self.x,self.Fy[:,idy,13]]))
                np.savetxt(tofolder+"/Ext_"+str(yv),np.transpose([self.x,self.Ey[:,idy,13]]))

                np.savetxt(tofolder+"/Fy_"+str(yv),np.transpose([self.x,self.Fy[idy,:,13]]))
                np.savetxt(tofolder+"/Ey_"+str(yv),np.transpose([self.x,self.Ey[idy,:,13]]))

        for idx,xv in enumerate(self.x): # [-1.3 -1.2 -1.1 -1.  -0.9 -0.8 -0.7 -0.6 -0.5 -0.4 -0.3 -0.2 -0.1  0.   0.1  0.2  0.3  0.4  0.5  0.6  0.7  0.8  0.9  1.   1.1  1.2  1.3]
            print "x:",xv,"___"
            for idy,yv in enumerate(self.y): #
                for idz,zv in enumerate(self.z):
                    #print "x:",xv,"y:",yv,"z:",zv
                    # a) fix y This is always equal!
                    self.spline_for_z_y(z = zv, y = yv)

                    # b) get at certain (defined) x
                    self.fmlx[idx,idy,idz] = self.splfmlx_f(xv)
                    self.fmly[idx,idy,idz] = self.splfmly_f(xv)
                    self.fmlz[idx,idy,idz] = self.splfmlz_f(xv)

                    self.emlx[idx,idy,idz] = self.splfmlx_e(xv)
                    self.emly[idx,idy,idz] = self.splfmly_e(xv)
                    self.emlz[idx,idy,idz] = self.splfmlz_e(xv)

                    self.Fx[idx,idy,idz] = self.splFx(xv)
                    self.Fy[idx,idy,idz] = self.splFy(xv)
                    self.Fz[idx,idy,idz] = self.splFz(xv)

                    self.Ex[idx,idy,idz] = self.splEx(xv)
                    self.Ey[idx,idy,idz] = self.splEy(xv)
                    self.Ez[idx,idy,idz] = self.splEz(xv)

                    #print "in:",x,yv,z,fx,ex
                    if self.verbose and zv == 0.0:
                        print "x:",xv,"y:",yv,"z:",zv , '\t',"fx:",round(fx,9),'\t',"ex:",round(ex,9),"||"
        self.save_npz_fml_F(folder)
        return

    def save_npz_fml_F(self,folder):
        print utils.printred("saving "+folder+'/qubus.fml.npz')
        np.savez_compressed(folder+"/qubus.fml.npz",
                x    = self.x,
                y    = self.y,
                z    = self.z,
                fmlx = self.fmlx,
                fmly = self.fmly,
                fmlz = self.fmlz,
                emlx = self.emlx,
                emly = self.emly,
                emlz = self.emlz
                #splfmlx_achsenabschnitt = self.splfmlx_achsenabschnitt,
                #splfmly_achsenabschnitt = self.splfmly_achsenabschnitt,
                #splfmlz_achsenabschnitt = self.splfmlz_achsenabschnitt
                )


        #self.splfmlx_achsenabschnitt  = self.splfmlx_f.antiderivative()(0)
        #self.splfmly_achsenabschnitt  = self.splfmly_f.antiderivative()(0)
        #self.splfmlz_achsenabschnitt  = self.splfmlz_f.antiderivative()(0)
        #self.Fx_achsenabschnitt  = self.splFx.antiderivative()(0)
        #self.Fy_achsenabschnitt  = self.splFy.antiderivative()(0)
        #self.Fz_achsenabschnitt  = self.splFz.antiderivative()(0)

        print utils.printred("saving "+folder+'/qubus.F.npz')
        np.savez_compressed(folder+"/qubus.F.npz",
                x    = self.x,
                y    = self.y,
                z    = self.z,
                Fx   = self.Fx,
                Fy   = self.Fy,
                Fz   = self.Fz,
                Ex   = self.Ex,
                Ey   = self.Ey,
                Ez   = self.Ez
                #Fx_achsenabschnitt = self.Fx_achsenabschnitt,
                #Fy_achsenabschnitt = self.Fy_achsenabschnitt,
                #Fz_achsenabschnitt = self.Fz_achsenabschnitt
                )
        return

    def ef_for_z(self, x = False, y = False, z = False):
        ''' z has to be in the pkl
            y and x can be chosen freely
        for every y in the plane do the fitting at x '''
        self.spline_for_z_y(z = 0.0, y = 0.0)
        self.get_achsenabschnitt_at_0()
        if self.verbose:
            print "x \t y \t z \t f \t e"
            print "z=0, y=0:",self.splfmlx_achsenabschnitt
        y_f_e = False

        #################################################################################
        # THE FOLLOWING LOOP GOES THROUTH PREDEFINED Y AND Z
        #################################################################################
        # here we go through all y's [-1.3, ...1.3]  for fixed z [0.5]
        # at this point we could introduce a lookup table with the coeffitients
        # for every y we have ca 27 coefs and there are 27 lines per z plane
        # make a matrix and explicitly save the whole splinemethod
        for yv in self.y:
            # a) fix y This is always equal!
            self.spline_for_z_y(z = z, y = yv)


            # b) get at certain (defined) x
            fx = self.splfmlx_f(x)
            fy = self.splfmly_f(x)
            fz = self.splfmlz_f(x)
            ex = self.splfmlx_e(x)
            ey = self.splfmly_e(x)
            ez = self.splfmlz_e(x)
            #print "in:",x,yv,z,fx,ex
            if self.verbose:
                print "x:",x,"y:",yv,"z:",z , '\t',"fx:",round(fx,9),'\t',"ex:",round(ex,9),"||"

            y_f_e = utils.append_row_to_2d_array(inarray = y_f_e, addrow = [yv, fx,fy,fz, ex,ey,ez])
        #print "schleife done"
        #print y_f_e
        fx  = UnivariateSpline( y_f_e[:,0], y_f_e[:,1], k=2, s=0.0)(y)
        fy  = UnivariateSpline( y_f_e[:,0], y_f_e[:,2], k=2, s=0.0)(y)
        fz  = UnivariateSpline( y_f_e[:,0], y_f_e[:,3], k=2, s=0.0)(y)

        ex  = UnivariateSpline( y_f_e[:,0], y_f_e[:,4], k=2, s=0.0)(y)
        ey  = UnivariateSpline( y_f_e[:,0], y_f_e[:,5], k=2, s=0.0)(y)
        ez  = UnivariateSpline( y_f_e[:,0], y_f_e[:,6], k=2, s=0.0)(y)
        if self.verbose:
            print utils.printred("x: "+str(x)+" y: "+str(y)+" z: "+str(z)+" fx: "+str(fx)+" ex: "+str(ex))
        #print "DONE"
        #print ""
        return fx,fy,fz,ex,ey,ez

    def ef(self, x = False, y = False, z = False):
        # 1) get correct achsenabschnitt
        self.spline_for_z_y(z = 0.0, y = 0.0)
        self.get_achsenabschnitt_at_0()
        z_f_e = False
        #print "yo",self.z
        for zv in self.z:
            fx,fy,fz,ex,ey,ez = self.ef_for_z(x = x, y = y, z = zv)
            #print "fxK",fx
            #print "eyK",ey
            z_f_e = utils.append_row_to_2d_array(inarray = z_f_e, \
                    addrow = [zv, fx,fy,fz, ex,ey,ez])
            if self.verbose:
                print "x:",x,"y:",y,"z:",zv , '\t',"fx:",round(fx,9),'\t',"ex:",round(ex,9),"|"

        fx  = UnivariateSpline( z_f_e[:,0], z_f_e[:,1], k=2, s=0.0)(z)
        fy  = UnivariateSpline( z_f_e[:,0], z_f_e[:,2], k=2, s=0.0)(z)
        fz  = UnivariateSpline( z_f_e[:,0], z_f_e[:,3], k=2, s=0.0)(z)

        ex  = UnivariateSpline( z_f_e[:,0], z_f_e[:,4], k=2, s=0.0)(z)
        ey  = UnivariateSpline( z_f_e[:,0], z_f_e[:,5], k=2, s=0.0)(z)
        ez  = UnivariateSpline( z_f_e[:,0], z_f_e[:,6], k=2, s=0.0)(z)
        if self.verbose:
            print utils.printred("x: "+str(x)+" y: "+str(y)+" z: "+str(z)+" fx: "+str(fx)+" ex: "+str(ex))
        self.f = np.array([fx,fy,fz])
        self.e = np.array([ex,ey,ez])
        return

    def efnew(self, x = False, y = False, z = False):
        #self.xx, self.yy = np.meshgrid(self.x,self.y)
        #self.xxx, self.yyy = np.mgrid[-1.3:1.3:27j, -1.3:1.3:27j]

        self.fmlx = np.zeros((self.x.size,self.y.size,self.z.size))
        self.fmly = np.zeros((self.x.size,self.y.size,self.z.size))
        self.fmlz = np.zeros((self.x.size,self.y.size,self.z.size))
        for idx, x in enumerate(self.x):
            print x
            for idy, y in enumerate(self.y):
                for idz, z in enumerate(self.z):
                    self.fmlx[idx,idy,idz] = self.data[(self.data.x == x) & \
                            (self.data.y == y) & (self.data.z == z)]['fmlx'].values
                    self.fmly[idx,idy,idz] = self.data[(self.data.x == x) & \
                            (self.data.y == y) & (self.data.z == z)]['fmly'].values
                    self.fmlz[idx,idy,idz] = self.data[(self.data.x == x) & \
                            (self.data.y == y) & (self.data.z == z)]['fmlz'].values
        np.savez_compressed("/Users/glensk/Dropbox/Understand_distributions/fall.npz",
                fmlx = self.fmlx,
                fmly = self.fmly,
                fmlz = self.fmlz,
                x    = self.x,
                y    = self.y,
                z    = self.z,
                )
        var = np.load("/Users/glensk/Dropbox/Understand_distributions/fall.npz")
        self.fall = var['fall']

        xidx = np.nonzero(self.x == x)[0][0]
        yidx = np.nonzero(self.y == y)[0][0]
        zidx = np.nonzero(self.z == z)[0][0]
        print 'x:',x,'xidx:',xidx
        print 'y:',y,'yidx:',yidx
        print 'z:',z,'zidx:',zidx
        print 'q.fall[',xidx,yidx,zidx,']',self.fall[xidx,yidx,zidx]
        print RectBivariateSpline(self.x, self.y, self.fall[zidx],kx=2,ky=2,s=0.0)(x,y)

if __name__ == '__main__':
    q = qubus_splines()
    folder = "/Users/glensk/Dropbox/Understand_distributions/displacements_dense/Ir/2x2x2sc_qubus_3x3x3kp"
    #q.loadQubuspd("/Users/glensk/Dropbox/Understand_distributions/displacements_dense/Ir/2x2x2sc_qubus_3x3x3kp/qubus.pkl")
    q.loadQubuspd(folder+"/qubus.pkl")

    q.ef_in_whole_qubus(folder)
    q.kb()
    q.save_npz_fml_F(folder)
    sys.exit()

    #print "spline fit:"
    x = -0.73
    y = 0.01
    z = 0.932
    #q.ef(x = x, y = y, z = z)

    #print "f:",q.f
    #print "e:",q.e

    def load_qubus_npz(filename="/Users/glensk/Dropbox/Understand_distributions/displacements_dense/Ir/2x2x2sc_qubus_3x3x3kp/qubus.npz"):
        print utils.printred("loading "+filename)
        var = np.load(filename)
        q.fmlx = var['fmlx']
        q.fmly = var['fmly']
        q.fmlz = var['fmlz']
        q.x = var['x']
        q.y = var['y']
        q.z = var['z']
        return

    def load_qubus_npz2(filename="/Users/glensk/Dropbox/Understand_distributions/feall.npz"):
        print utils.printred("loading "+filename)
        var = np.load(filename)
        q.ofmlx = var['fmlx']
        q.ofmly = var['fmly']
        q.ofmlz = var['fmlz']
        q.oemlx = var['emlx']
        q.oemly = var['emly']
        q.oemlz = var['emlz']
        q.ox = var['x']
        q.oy = var['y']
        q.oz = var['z']

    load_qubus_npz()
    load_qubus_npz2()
    sys.exit()
    #print "new:"
    #q.efnew(x = 0.7, y=0.8, z=0.9)
    print "new:"
    var = np.load("/Users/glensk/Dropbox/Understand_distributions/fall.npz")
    q.fall = var['fall']

    #q.fall[:,:,13] # z = 0.0
    #q.fall[:,:,17] # z = 0.4
    #q.fall[:,:,23] # z = 1.0

    # x,y,z == [ -1.3 .. 1.3 ] --> 0, 26
    # wobei man hier den fit nehmen sollte von -1.6 bis 1.6 ne doch nicht
    #0  == -1.3 |  1.3 -1.3
    #26 ==  1.3 |  2.6 -1.3
    def xyz_to_map_coords(x):
        ''' x = -1.3 ---> out = 0
            x = 0    ---> out = 13
            x = 1.3  ---> out = 26 '''
        return (x+1.3)*10

    def map_coords_to_xyz(map_coords):
        return map_coords*0.1-1.3

    #for i in np.arange(q.x.size):
    #for i in q.x:
    #    print i,xyz_to_map_coords(i),map_coords_to_xyz(xyz_to_map_coords(i))
    import scipy

    #idxz = np.nonzero(q.x == z)[0][0]
    #data = q.fall[:,:,idxz]
    #def gc(x,y):
    #    coords = np.array([[xyz_to_map_coords(x), xyz_to_map_coords(y)]])
    #    coords = coords.T
    #    zi = scipy.ndimage.map_coordinates(data, coords, order=3, mode='nearest')
    #    return zi
    #print gc(x,y)

    def gcc(x,y,z):
        coords = np.array([[xyz_to_map_coords(x), xyz_to_map_coords(y),xyz_to_map_coords(z)]])
        coords = coords.T
        zi = scipy.ndimage.map_coordinates(q.fall, coords, order=2, mode='nearest')
        return zi
    print gcc(x,y,z)
