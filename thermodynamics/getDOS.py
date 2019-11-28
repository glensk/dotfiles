#!/usr/bin/env python

import sys
import os
import numpy as np



def get_dos(data, points = 250):
    ''' gets a 1d numpy array and creates a DOS out of that '''

    # get x grid
    dist_space = np.linspace( min(data), max(data), points )

    # define smearing
    fact = 3.0
    fact = 1.0
    #fact = 0.5
    #fact = 0.1
    #fact = 0.01
    def my_kde_bandwidth(obj, fac=fact):
        """We use Scott's Rule, multiplied by a constant factor."""
        return np.power(obj.n, -1./(obj.d+4)) * fac
    data = data.astype(np.double)
    #print data[:30]
    #print "data:>>",type(data),data.shape,type(data[0])
    #print ""
    ##data = np.random.rand(len(data))
    #print data[:30]
    #print "data:>>",type(data),data.shape,type(data[0])

    #print data
    from scipy.stats.kde import gaussian_kde
    kde = gaussian_kde( data, bw_method=my_kde_bandwidth)
    print "getting dos ... this may take a while"
    dos = np.array([dist_space, kde(dist_space)]).transpose()
    print "dos done ..."
    return dos

def get_and_save_dos(filename):
    print "loading "+filename+" ..."
    if os.path.isfile(filename) != True:
        sys.exit(filename+" does not exist; exit;")
    data = np.loadtxt(filename)
    dos = get_dos(data)
    np.savetxt(filename+"_DOS",dos)
    return dos

if __name__ == '__main__':
    #print sys.argv
    #print "len:",len(sys.argv)#,"third arg:",sys.argv[2]
    #temp = sys.argv[2]
    #print "temp:",temp,type(temp),float(temp)
    filename = sys.argv[1]
    dos = get_and_save_dos(filename)

    if False:
        #if len(sys.argv) > 2:
        #    temp = float(sys.argv[2])
        #    print "getting pot ..."
        #    import pot_parametrize
        #    pot = pot_parametrize.dos_to_pot(dos, temp)
        #    np.savetxt(filename+"_POT_"+str(temp),pot)
        temp = int(os.getcwd().split("/")[-1].split("_")[-1].split("K")[0])
        dos0 = np.loadtxt("DOS_dfn_0.0")
        dos0max = dos0[np.nonzero(dos0[:,1] == dos0[:,1].max())[0][0]][0]
        print "dos0max:",dos0max
        np.savetxt("harm",np.array([[temp,dos0max]]))
        dos0 = np.loadtxt("DOS_dfn_1.0")
        dos0max = dos0[np.nonzero(dos0[:,1] == dos0[:,1].max())[0][0]][0]
        print "dos0max:",dos0max
        np.savetxt("anharm",np.array([[temp,dos0max]]))
        sys.exit()

        import pot_parametrize
        temp = 1235
        dos0 = np.loadtxt("dfn_0.0_DOS")
        pot0 = pot_parametrize.dos_to_pot(dos0, temp)
        dos1 = np.loadtxt("dfn_1.0_DOS")
        pot1 = pot_parametrize.dos_to_pot(dos1, temp)
        dos0max = np.nonzero(dos0[:,1] == dos0[:,1].max())[0][0]
        dos1max = np.nonzero(dos1[:,1] == dos1[:,1].max())[0][0]
        dos0shifted = np.array([dos0[:,0] - dos0[dos0max,0],dos0[:,1]]).transpose()
        dos1shifted = np.array([dos1[:,0] - dos1[dos1max,0],dos1[:,1]]).transpose()

        pot0shifted = np.array([pot0[:,0] - pot0[dos0max,0],pot0[:,1]]).transpose()
        pot1shifted = np.array([pot1[:,0] - pot1[dos1max,0],pot1[:,1]]).transpose()

        np.savetxt("dfn_0.0_DOS_shifted",dos0shifted)
        np.savetxt("dfn_1.0_DOS_shifted",dos1shifted)
        np.savetxt("dfn_0.0_POT_1235_shiftedorig",pot0shifted)

        print pot0shifted[:4]
        print pot0shifted[-4:]
        print ""
        print pot1shifted[:4]
        print pot1shifted[-4:]
        print ""
        import utils
        #pot1shifted2 = utils.return_function2_on_grid_of_f1(pot0shifted,pot1shifted)  # tn
        pot0shifted2,pot1shifted2 = utils.return_fuctions_on_same_grid(pot0shifted,pot1shifted)
        print pot1shifted2[:4]
        print pot1shifted2[-4:]

        #np.savetxt("dfn_0.0_POT_1235_shifted",pot0shifted)
        #np.savetxt("dfn_1.0_POT_1235_shifted",pot1shifted)
        np.savetxt("dfn_0.0_POT_1235_shifted",pot0shifted2)
        np.savetxt("dfn_1.0_POT_1235_shifted",pot1shifted2)

        diff = np.array([pot0shifted2[:,0],pot1shifted2[:,1]-pot0shifted2[:,1]]).transpose()
        np.savetxt("dfn_1.0_POT_1235_shifted2-0.0_shifted2",diff)
