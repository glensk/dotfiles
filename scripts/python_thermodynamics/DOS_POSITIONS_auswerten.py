#!/usr/bin/env python

import numpy as np
import sys
import os
import glob
import utils
reload(utils)
import argparse
import inspect, re
import hesse
reload(hesse)
import pot_info_startjob
#reload(pot_info_startjob)

#atoms=32
deletefirst=30  # deletes first 30 structures of trajectory

np.set_printoptions(suppress=True)   # display arrays withou 000000
np.set_printoptions(precision=14)    # print only 6 digist after .




def getmapping(atoms,pos0,eqposfile='eqpos',cellfile='cell'):
    import crystal_generator
    crystal = crystal_generator.crystal()
    nnlist = np.zeros((atoms,atoms))
    nnlist[:] = np.nan
    for idxi,pos0i in enumerate(pos0):
        nnlisti = crystal.get_NNlist(idxi,1,coordfile_cart=eqposfile,cellfile=cellfile)
        #print idxi,nnlisti
        lennnlist = len(nnlisti)
        for idxj,val in enumerate(nnlisti):
            nnlist[idxi,idxj] = val+1
        #nnlist[idxi] = crystal.get_NNlist(idxi,1,coordfile_cart=eqposfile,cellfile=cellfile)

    #print nnlist
    # now remove nans
    nnlistout = np.zeros((atoms,lennnlist))
    nnlistout[:] = np.nan
    for idxi,i in enumerate(nnlistout):
        for idxj,j in enumerate(nnlistout[idxi]):
            nnlistout[idxi,idxj] = nnlist[idxi,idxj]
    return nnlistout

def makedos(foldername, lambd, dn, tmelt, d = None, dnstr = None):
    """
    foldername (str): where to put the dos (is created if it does not exist
    lambda     (str): lambdastring
    dn         (1d numpy array): containing d (distance, distance_long, ...)
    temelt     (integer       ): melting temperature
    d          (numpy array of vectors): corresponding matrix of vectors"""
    if os.path.isdir(foldername) != True:
        os.makedirs(foldername)
    def namestr(obj, namespace):
        return [name for name in namespace if namespace[name] is obj]
    if d != None:
        dstr = namestr(d, globals())[0]
        np.savetxt(foldername+"/"+dstr+"_"+lambd,d,fmt="%.6f %.6f %.6f")

    if dnstr == None:
        dnstr = namestr(dn, globals())[0]
    print "1/3 saving: "+dnstr+"_"+lambd
    np.savetxt(foldername+"/"+dnstr+"_"+lambd,dn,fmt="%.6f")


    hier=os.getcwd()
    os.chdir(foldername)
    print "2/3 getDOS.py "+foldername+"/"+dnstr+"_"+lambd
    utils.run2("getDOS.py "+dnstr+"_"+lambd)

    # push DOS slightly up so there are no negative values
    dos_read = np.loadtxt(dnstr+"_"+lambd+"_DOS")
    dos_read_min = dos_read[:,1].min()
    dos_read[:,1] = dos_read[:,1]+abs(dos_read_min)+0.000000000001
    np.savetxt("DOS",dos_read)

    utils.run2("mv DOS DOS_"+dnstr+"_"+lambd)
    os.chdir(hier)
    print "3/3 mv DOS "+foldername+"/DOS_"+dnstr+"_"+lambd
    dos = np.loadtxt(foldername+"/DOS_"+dnstr+"_"+lambd)

    pot_ = -np.log(dos[:,1])*tmelt*0.086173423
    pot_ = pot_+abs(pot_.min())
    pot = np.copy(dos)
    pot[:,1] = pot_
    np.savetxt(foldername+"/pot_"+dnstr+"_"+lambd,pot)
    potref = np.copy(dos)
    potref[:,1][:] = tmelt*0.086173423
    np.savetxt(foldername+"/potref_"+dnstr+"_"+lambd,potref)


    com = np.array([[np.mean(dn),0],[np.mean(dn),dos[:,1].max()]])
    np.savetxt(foldername+"/DOS_"+dnstr+"_"+lambd+"_center_mass",com,fmt="%.6f %.6f")

def varname(p):
    for line in inspect.getframeinfo(inspect.currentframe().f_back)[3]:
        m = re.search(r'\bvarname\s*\(\s*([A-Za-z_][A-Za-z0-9_]*)\s*\)', line)
        if m:
            return m.group(1)

def nearly_equal(a,b,sig_fig=5):
    return ( a==b or
    int(a*10**sig_fig) == int(b*10**sig_fig))

def nearly_equal(a,b,sig_fig=5):
    return ( a==b or
             int(a*10**sig_fig) == int(b*10**sig_fig)
           )

def nearly_equal(a,b,sig_fig=5):
    return ( a==b or
             int(a*10**sig_fig) == int(b*10**sig_fig)
           )

def getpos(schritte,pos,atoms):
    if type(schritte) == np.ndarray:
        for i in schritte:
            print pos[atoms*schritte:atoms*schritte+atoms]
    return pos[atoms*schritte:atoms*schritte+atoms]


def std():
    #print "args.folder:",args.folder,len(args.folder)
    startfolder = os.getcwd()
    element = startfolder.split("/ti/")[-1]
    element = element.split("/")[0]
    print utils.printred("ELEMENT:"+element)
    insgesamt =160
    print " "*(insgesamt+3)+utils.printred("step")+"  "+utils.printgreen("std   (sup)")+" "+utils.printblue("dudl +UP  ")
    print " "*(insgesamt+3)+utils.printred("----")+"  "+utils.printgreen("-----------")+" "+utils.printblue("----------")
    listone = args.folder
    listtwo = utils.list_sorted(args.folder)


    dataarray = False # initialize with False; None will not work
    #for foldername in listone:
    for foldername in listtwo:
        #os.chdir(startfolder)
        #os.chdir(foldername)
        #print "###############################################"*2
        #print "#",foldername
        #print "###############################################"*2

        filename=args.std[0]
        import warnings
        warnings.filterwarnings("ignore")
        if os.path.isfile(foldername+"/"+filename) != True:
            print utils.printred(foldername+" "+filename+ ' does not exist')
            continue

        try:
            dudl=np.loadtxt(foldername+"/"+filename)
        except ValueError:
            print utils.printred(foldername+"/"+filename+ ' something went wrong with loading dudl')
            continue
        except UserWarning:
            print utils.printred(foldername+"/"+filename+ ' something went wrong with loading dudl (dudl empty)')
            continue
        #print "shpe:",dudl.shape
        if dudl.shape[0] == 0:
            print utils.printred(foldername+"/"+filename+ ' is empty')
            continue
        if len(dudl.shape) == 1:
            print utils.printred(foldername+"/"+filename+ ' is has just one line')
            continue

        if dudl.shape[0] <= 2:
            print utils.printred(foldername+"/"+filename+ ' is has just two line')
            continue

        np.set_printoptions(precision=6)    # print only 6 digist after .
        np.set_printoptions(linewidth=200)
        #print dudl
        #print "dudl.shape[1]:",dudl.shape[1]
        std=np.zeros((dudl.shape[0],2))
        dudlout=np.zeros((dudl.shape[0],2))
        for idx,i in enumerate(dudl):
            std[idx,0] = idx
            dudlout[idx,0] = idx
            if idx == 0:
                std[idx,1] = 0.0
                #if filename == "dUdL" or "dUdL" in filename or "dudl" in filename or dudl.shape[1] == 9:
                if dudl.shape[1] == 9:
                    #print "filename:",filename
                    #print "dudl.shape",dudl.shape
                    #print "len(dudl.shape)",len(dudl.shape)
                    #print "dudl.shape[0]",dudl.shape[0]
                    #print "dudl[:,1]:",dudl[:,1]
                    #print "dudl[:,4]:",dudl[:,4]
                    #print "dudl[:,5]:",dudl[:,5]

                    dudlout[idx,1] = (dudl[:,4]-dudl[:,5])[0]  # for dUdL file
                else:
                    dudlout[idx,1] = (dudl[:,0]-dudl[:,1])[0]
            else:
                #if filename == "dUdL" or "dUdL" in filename or "dudl" in filename or dudl.shape[1] == 9:
                #print "__> dudl.shape[1]:",dudl.shape[1]
                if dudl.shape[1] == 9:
                    u = dudl[:,[0,4]]
                    uref = dudl[:,[0,5]]
                    std[idx,1] = (dudl[:,4]-dudl[:,5])[:idx].std()
                    dudlout[idx,1] = (dudl[:,4]-dudl[:,5])[:idx].mean()
                    #print "__> dudl.shape[1]:",dudl.shape[1],u,uref,std[idx,1]
                else:
                    std[idx,1] = (dudl[:,0]-dudl[:,1])[:idx].std()
                    dudlout[idx,1] = (dudl[:,0]-dudl[:,1])[:idx].mean()

        #print "std:",std[:10]
        #sys.exit()
        #lambd9 = foldername.split("/lambda")[-1]
        #lambd8 = lambd9.split('/')[0]
        #lambd  = float(lambd8)
        #print 'lambd:',lambd

        if std[-1][1] > 999:
            std[-1][1] = 999

        stepsoutnr = int(std[-1][0])
        stepsout = utils.number_to_string_of_certain_length(stepsoutnr,0,4)

        stdoutnr = std[-1][1]
        stdout = utils.number_to_string_of_certain_length(stdoutnr,1,5)

        #####################################################################################
        ## get dUdL from hesse reference
        #####################################################################################
        low_std = False
        lowplushigh_dudl = False
        uptildout = "--"
        printit = False
        lambd = foldername.split("/lambda")[-1]
        if printit:
            print "foldername:",type(foldername),foldername.split("/lambda")
            print "000:",lambd
        lambd = foldername.split("/lambda")[-1].split("/")
        if printit:
            print "1init",lambd
        lambd = foldername.split("/lambda")[-1].split("/")[0]
        if printit:
            print foldername.split("/lambda")[-1].split("/")
            print "2",lambd
        lambd = foldername.split("/lambda")[-1].split("/")[0].split("lambda")
        if printit:
            print "3",lambd
        lambd = foldername.split("/lambda")[-1].split("/")[0].split("lambda")[-1]
        if printit:
            print "4",lambd
        lambd = lambd.split("_")[0]
        if printit:
            print "5",lambd[0]
        #print "lambd:",lambd
        lambd = 1.0
        lambd = float(lambd)
        #lambd = float(foldername.split("/lambda")[-1].split("/")[0].split("lambda")[-1])
        #print "lambd:",lambd

        TOX = foldername.split("_TOX_")[-1].split("_TOY_")[0]
        if TOX == 'NONE':
            TOX = False
        else:
            TOX = True
        #print "tox:",TOX,type(TOX)
        #element = ""
        #POTCAR = foldername+"/POTCAR"
        #if os.path.isfile(POTCAR) == True:
            #element = utils.read_only_certain_lines_of_file(POTCAR,0,0)[0].split()[1]
        lowfile = glob.glob("/Users/glensk/Dropbox/Understand_distributions/jobvorlage_all/2x2x2sc_"+element+"_*/avg_dUdL_low_fre")
        lowhighfile = glob.glob("/Users/glensk/Dropbox/Understand_distributions/jobvorlage_all/2x2x2sc_"+element+"_*/avg_dUdL_lowplushigh_fre")
        #print "len:",len(lowfile),lowfile
        if len(lowfile) == 1:
            low = np.loadtxt(lowfile[0])
            #print "lambd:",lambd
            #print "ow:",low
            for i in low:
                #print i
                if i[0] == lambd:
                    low_dudl = i[1]
                    low_std = i[3]
                    #print ">> std:",low_std
        #print "###########"
        if len(lowhighfile) == 1:
            low = np.loadtxt(lowhighfile[0])
            for i in low:
                #print i
                if i[0] == lambd:
                    lowplushigh_dudl = i[1]

        #print element,"low_std:",low_std,"lowplushigh_dudl:",lowplushigh_dudl
        speedupstr = str(0)
        speedup = 0
        #print "low_std:",low_std
        if type(low_std) != bool:
            speedup = low_std**2./stdoutnr**2.
            #print "low_std:",low_std,"stdoutnr:",stdoutnr,"speedup:",speedup
            speedupstr = str(int(speedup))
            speedupadd = 3 - len(speedupstr)
            speedupstr = speedupstr + " "*speedupadd
            if type(lowplushigh_dudl) != bool:
                uptild = lowplushigh_dudl - low_dudl
                uptildout = utils.number_to_string_of_certain_length(uptild,1,4)

        #print element,"low_std:",low_std,"lowplushigh_dudl:",lowplushigh_dudl,"uptild:",uptildout
        #print ""


        dudloutnr = dudlout[-1][1]
        if dudlout[-1][1] > 999:
            dudlout[-1][1] = 999
        if dudlout[-1][1] == 999:
            if len(dudlout) > 1500:
                if dudlout[1500][1] < 999:
                    dudloutnr = dudlout[1500][1]
        dudlcurrent = np.array([u[:,0],u[:,1]-uref[:,1]]).transpose()
        dudlcurrentuncor = dudlcurrent[::10]
        #dudlcurrentuncor_avg_std = [ std(dudlcurrentuncor[:,1][:i]) for i in dudlcurrentuncor ]
        #for jidx,j in enumerate(dudlcurrentuncor):
        dudlcurrentuncor_avg_std = [ dudlcurrentuncor[:jidx,1].std() for jidx in range(1,len(dudlcurrentuncor)) ]
        #dudlcurrentuncor_avg_std = np.array([j,dudlcurrentuncor_avg_std[j])
        #print dudlcurrentuncor_avg_std

        dudloutstr= utils.number_to_string_of_certain_length(dudloutnr,1,6)

        foldernameout = foldername[:insgesamt]
        foldernameoutadd = insgesamt - len(foldernameout)
        foldernameout = foldernameout+ " "*foldernameoutadd
        zeichen = len(foldername)+len(filename)
        leer = insgesamt - zeichen


        # dataarray
        #print  " np.array([[lambd, stepsoutnr, stdoutnr, speedup, dudloutnr]])"
        addrow = [lambd, stepsoutnr, stdoutnr, speedup, dudloutnr, TOX]
        #print "addrow:",addrow
        #print "dataarray:",dataarray
        dataarray = utils.append_row_to_2d_array(inarray = dataarray, addrow=addrow)
        if len(foldername) > insgesamt:
            print foldername,filename
            print " "*(insgesamt+2),utils.printred(stepsout),utils.printgreen(stdout+" ("+speedupstr+")"),utils.printblue(dudloutstr+" "+uptildout)
        if len(foldername) <= insgesamt:
            print foldername,filename," "*leer,utils.printred(stepsout),utils.printgreen(stdout+" ("+speedupstr+")"),utils.printblue(dudloutstr+" "+uptildout)

        if dudl.shape[0] == 0:
            print foldername,filename,"dUdL is empty!"
            if os.path.isfile(foldername+"/avg_std_"+filename+"_schritte") == True:
                os.remove(foldername+"/avg_std_"+filename+"_schritte")
            if os.path.isfile(foldername+"/avg_dudl_"+filename+"_schritte") == True:
                os.remove(foldername+"/avg_dudl_"+filename+"_schritte")
        else:
            #print "avg_std_"+filename+"_schritte written .."
            #print "avg_dudl_"+filename+"_schritte written .."
            try:
                np.savetxt(foldername+"/avg_std_"+filename+"_schritte",std)
                np.savetxt(foldername+"/avg_dudl_"+filename+"_schritte",dudlout)
            except IOError:
                pass

    #try:
    #    np.savetxt("dataarray_",dataarray)
    #    #np.savetxt("dataarray_dudl0.0_vs_std",dataarray[np.nonzero((dataarray[:,0] == 0.) & (dataarray[:,1] > 400)  & (dataarray[:,2] < 200.) )[0]][:,[2,4]])
    #    #np.savetxt("dataarray_dudl1.0_vs_std",dataarray[np.nonzero((dataarray[:,0] == 1.) & (dataarray[:,1] > 400)  & (dataarray[:,2] < 200.) )[0]][:,[2,4]])
# & #    (dataarray[:,2] >= 999)
    #    np.savetxt("dataarray_dudl0.0_vs_std_with_tox",dataarray[np.nonzero((dataarray[:,0] == 0.) & (dataarray[:,5] == 1.) & (dataarray[:,1] > 400)  & (dataarray[:,2] < 200.) )[0]][:,[2,4]])
    #    np.savetxt("dataarray_dudl0.0_vs_std_wout_tox",dataarray[np.nonzero((dataarray[:,0] == 0.) & (dataarray[:,5] == 0.) & (dataarray[:,1] > 400)  & (dataarray[:,2] < 200.) )[0]][:,[2,4]])
    #    np.savetxt("dataarray_dudl1.0_vs_std_with_tox",dataarray[np.nonzero((dataarray[:,0] == 1.) & (dataarray[:,5] == 1.) & (dataarray[:,1] > 400)  & (dataarray[:,2] < 200.) )[0]][:,[2,4]])
    #    np.savetxt("dataarray_dudl1.0_vs_std_wout_tox",dataarray[np.nonzero((dataarray[:,0] == 1.) & (dataarray[:,5] == 0.) & (dataarray[:,1] > 400)  & (dataarray[:,2] < 200.) )[0]][:,[2,4]])


    #    #np.savetxt("dataarray_dudl1.0_vs_std",dataarray[np.nonzero(dataarray[:,0] == 1.)[0]][:,[2,4]])
    #    #for foldername in listone:
    #    #    print foldername
    #    #    continue
    #    pot = pot_info_startjob.pot()
    #    pot.element    = element             # Rh, Al, ...
    #    pot.folder_displacement_from = "displacements"
    #    pot.folder_displacement_from = "displacements_dense"
    #    pot.sc_forparams = 2            # parameters for parametrization are taken from here (2x2x2sc disp quer xdir)
    #    pot.kpstring = pot.get_element_2x2x2sc_info(pot.element)[1]   # gets 3x3x3kp string
    #    pot.sc_create = 2               # create a 2x2x2sc or a 3x3x3sc ? (is alredy used in pot.init_variables)
    #    #pot.init_variables( kpstringoptional = pot.kpstring)  # only to get avg_dUdL_low_fre (wholefile)
    #    print "done:"
    #    np.savetxt("dataarray_dudl_vs_std_harmonic",pot.jobvorlage_avg_dudl_low_fre[:,[3,1]])
    #except IOError:
    #    pass
    return dataarray

def dudlvseps():
    print args.dudlvseps
    #enesall= glob.glob("ene_vs_eps_all_lj_2.9203510063")
    enesall=glob.glob("ene_vs_eps_all_morse_aa_1.544")
    #enesall= glob.glob("ene_vs_eps_all_lj_2.86378246381")
    #enesall=glob.glob("ene_vs_eps_all")
    for i in enesall:
        #enes = np.loadtxt("ene_vs_eps_all")
        #eps = np.loadtxt("ene_vs_eps_all_eps")
        enes = np.loadtxt(i)
        eps = np.loadtxt(i+"eps")
        outstd = np.zeros((enes.shape[1]-1,2))
        outmean = np.zeros((enes.shape[1]-1,2))
        for epsidx,eneidx in enumerate(np.arange(1,enes.shape[1])):
            dudl = enes[:,0]-enes[:,eneidx]
            outstd[epsidx,0 ] = eps[epsidx]
            outstd[epsidx,1 ] = dudl.std()
            outmean[epsidx,0 ] = eps[epsidx]
            outmean[epsidx,1 ] = dudl.mean()
        print "-----------",i
        np.savetxt(i+"sigma_vs_eps",outstd)
        np.savetxt(i+"dUdL_vs_eps",outmean)



def dudlposc(dudlposcm = False, sb = False, se = False, verbose = False):
    if os.path.isfile("POSITIONs") != True:
        sys.exit("need POSITIONs")
    if os.path.isfile("dUdL") != True:
        sys.exit("need dUdL")
    if os.path.isfile("cell") != True:
        sys.exit("need cell")
    cell = np.loadtxt("cell")
    if os.path.isfile("EqCoords_direct") != True:
        sys.exit("need EqCoords_direct")
    pos0rel = np.loadtxt("EqCoords_direct")
    pos = np.loadtxt("POSITIONs")
    listdudl = np.loadtxt("dUdL")
    eqpos_rrel = np.loadtxt("EqCoords_direct")

    FEXPLODE = False #True
    h = False
    if FEXPLODE:
        h = hesse.read_Hessematrix("HesseMatrix_sphinx")


    ########################################################################
    # from here we dont have to load anything
    ########################################################################
    listdudlnew = np.copy(listdudl)
    listdudlnew[:,5] = np.nan  # Uref
    listdudlnew[:,6] = np.nan  # dUdL
    schritte_dudl = listdudl.shape[0]+1

    atoms = eqpos_rrel.shape[0]
    print atoms
    schritte = pos.shape[0]/atoms
    print "schritte:",schritte
    print "atoms:",atoms
    out = np.zeros((schritte,atoms))

    ########################################################################
    # in case we also want the harmonic forces
    if dudlposcm:
        for i in np.arange(schritte):
            if i > 0 and type(sb) != bool and i < sb:
                continue
            if i > 0 and type(se) != bool and i > se:
                continue
            posforceaktuell = getpos(i,pos,atoms)
            posaktuell = posforceaktuell[:,:3]
            if i == 0:
                pos0 = np.copy(posaktuell)
            print "i:",i,"/",schritte
            diffpos = posaktuell - pos0
            du = hesse.get_energy_forces(pot="h",hessefile="HesseMatrix_sphinx",coord_cart=posaktuell,coord0_cart=pos0,returndu=True)
            #print "du:",du
            for a in np.arange(len(du)):
                #print "i,a:",i,a
                out[i,a] = np.linalg.norm(du[a])
        #print out
        for idx in np.arange(atoms):
            print idx
            np.savetxt("atom_"+str(idx),out[:,idx])
        sys.exit()

    fhar_vs_steps = False
    fdft_vs_steps = False
    fref_vs_steps = False
    eref_vs_steps = False
    elon_vs_steps = False
    etox_vs_steps = False


    for i in np.arange(schritte):
        if i > 0 and type(sb) != bool and i < sb:
            continue
        if i > 0 and type(se) != bool and i > se:
            continue
        posforceaktuell = getpos(i,pos,atoms)
        posaktuell = posforceaktuell[:,:3]
        fDFT = posforceaktuell[:,3:]
        print "CURRENT:",os.getcwd()
        energy,forces,energymevlong,energymevtopx = hesse.get_energy_forces_pairpot(usepot=['extern'],coord_cart = posaktuell, coord0_rrel = pos0rel, cell = cell, verbose=True)
        #print "|||||| energy,forces:",energy,forces,energymevlong,energymevtopx


        energy = energy
        fREF = forces
        if FEXPLODE:
            fHAR,eHAR = hesse.get_energy_forces(
                pot='h',
                pot2=False,
                #potparam=args.potparam,
                #hessefile = "HesseMatrix_sphinx",
                h = h,
                #hessefile1nn = args.inputfile1nn,
                #coordfile_cart = "cartesian_coords",
                coord_cart = posaktuell,
                only_return_forces_harmonic = True,
                verbose=verbose
                )
            fhar_vs_steps = utils.append_row_to_2d_array(fhar_vs_steps, np.insert(fHAR.flatten(),0,i).flatten())
        energymev = energy*1000/(atoms-1)
        print "listdudl.shape:",listdudl.shape,"i:",i

        eDFT  = listdudl[i-1][4]
        eREF = listdudl[i-1][5]
        print utils.printred(">>>>>>"*4+" DOS_POSITIONS_auswerten.py "+"<<<<<<"*4)
        print utils.printred(">>>>>>>",i,"energymev:",energymev,"eDFT:",eDFT,"eREF:",eREF)
        print utils.printred(">>>>>>"*4+" DOS_POSITIONS_auswerten.py "+"<<<<<<"*4)
        listdudlnew[i-1][5] = energymev
        listdudlnew[i-1][6] = eDFT - energymev
        listdudlsave = listdudlnew[:i]
        np.savetxt("dUdLnew",listdudlsave,fmt="%7.2f%10.1f%9.1f%9.1f%14.2f%16.2f%14.2f%10.2f%10.2f",header=" step   time(fs)  temp(K) average       U(meV/at)    Uref          dUdL   average    offset")

        eref_vs_steps = utils.append_row_to_2d_array(eref_vs_steps, np.array([i,energymevlong+energymevtopx]))
        elon_vs_steps = utils.append_row_to_2d_array(elon_vs_steps, np.array([i,energymevlong]))
        etox_vs_steps = utils.append_row_to_2d_array(etox_vs_steps, np.array([i,energymevtopx]))
        fdft_vs_steps = utils.append_row_to_2d_array(fdft_vs_steps, np.insert(fDFT.flatten(),0,i).flatten())
        fref_vs_steps = utils.append_row_to_2d_array(fref_vs_steps, np.insert(fREF.flatten(),0,i).flatten())



        if i in np.arange(schritte)[::5]:
            np.savetxt("elon_vs_steps",elon_vs_steps,fmt="%.0f %.6f")
            np.savetxt("etox_vs_steps",etox_vs_steps,fmt="%.0f %.6f")
            np.savetxt("fdft_vs_steps",fdft_vs_steps,fmt="%.6f")
            np.savetxt("fref_vs_steps",fref_vs_steps,fmt="%.6f")
            ka = fdft_vs_steps[:,1:].flatten()
            print "shape.ka:",ka.shape
            kb = fref_vs_steps[:,1:].flatten()
            #print "shape.kb:",kb.shape
            fref_vs_fdft = np.zeros((ka.shape[0],2))
            print "shape.kb:",kb.shape
            fref_vs_fdft[:,0] = ka
            fref_vs_fdft[:,1] = kb
            np.savetxt("fref_vs_fdft.dat",fref_vs_fdft,fmt="%.6f %.6f")

            minimum = ka.min();
            if kb.min() < minimum:
                minimum = kb.min()
            maximum = ka.max();
            if kb.max() > maximum:
                maximum = kb.max()
            np.savetxt("fdft_vs_fdft.dat",np.array([[minimum*1.2,minimum*1.2],[maximum*1.2,maximum*1.2]]),fmt="%.6f %.6f")
            if FEXPLODE:
                kc = fhar_vs_steps[:,1:].flatten()
                fhar_vs_fdft = np.zeros((kc.shape[0],2))
                fhar_vs_fdft[:,0] = ka
                fhar_vs_fdft[:,1] = kc
                np.savetxt("fhar_vs_fdft.dat",fhar_vs_fdft,fmt="%.6f %.6f")


    np.savetxt("elon_vs_steps",elon_vs_steps,fmt="%.6f")
    np.savetxt("etox_vs_steps",etox_vs_steps,fmt="%.6f")
    np.savetxt("fdft_vs_steps",fdft_vs_steps,fmt="%.6f")
    np.savetxt("fref_vs_steps",fref_vs_steps,fmt="%.6f")


def dudl():
    #################################################################################
    # start editing here
    #################################################################################
    atoms=108
    pot="m"

    # LJ
    if pot == "l":
        alleps = np.arange(0.06,0.15,0.01)

    # Morse
    if pot == 'm':
        alleps = np.arange(0.01,0.36,0.02)  # De  # bis 0-4 checken
        alleps = np.arange(0.05,0.2,0.01)   # De  # bis 0-4 checken
        allaa = np.arange(1.5,2.2,0.1)    # aa  # bis 2 checken
        alleps = np.array([0.17,0.18,0.19,0.2,0.21,0.22])   # De  # bis 0-4 checken
        #0.211, 1.544, 2.92
        alleps = np.array([111])   # De  # bis 0-4 checken
        allaa = np.array([111])  # 0.02 0.04 0.06 0.08 0.1
        print "allaa",allaa

    NN = 4.13/np.sqrt(2.)
    print "NN:",NN
    #NN = 2.852
    #################################################################################
    # stop editing here
    #################################################################################



    #folders = glob.glob("lambda1.0_"+"*")   # do all lambda folder
    folders=[os.getcwd()] # do only current folder


    import hesse
    hierbegin=os.getcwd()
    for makefolder in folders:
        print "makefolder:",makefolder
        os.chdir(hierbegin)
        os.chdir(makefolder)
        ene_vs_eps_all = False
        if os.path.isfile("structures_vasprun.gz") == True:
            utils.run2("gunzip structures_vasprun.gz")
        if os.path.isfile("structures_vasprun") != True:
            print "structures_vasprun dne"
            continue
        if os.path.isfile("structures_vasprun_dudl") == True:
            os.remove("structures_vasprun_dudl")
        utils.run2("awk \'{print $1,$2,$3,$4,$34}\' structures_vasprun > structures_vasprun_topy")
        if os.path.isfile("structures_vasprun_topy") != True or os.path.isfile("dUdL") != True:
            print "structures_vasprun_topy or dUdL dne"
            continue
        if os.path.isfile("POSITIONs") != True:
            print "POSITIONs dne"
            continue
        listposfre = np.loadtxt("structures_vasprun_topy")
        listdudl = np.loadtxt("dUdL")
        positions = np.loadtxt("POSITIONs")
        schritte_listposfre=np.shape(listposfre)[0]/atoms
        schritte_positions = np.shape(positions)[0]/atoms
        schritte_dudl = listdudl.shape[0]+1

        if schritte_listposfre != schritte_positions:
            sys.exit("schritte_listposfre:"+str(schritte_listposfre)+"schritte_positions:"+str(schritte_positions))
        if schritte_listposfre != schritte_dudl:
            sys.exit("schritte_listposfre:"+str(schritte_listposfre)+"schritte_dudl:"+str(schritte_dudl))
        print "schritte_listposfre:",schritte_listposfre
        print "schritte_positions :",schritte_positions
        print "listdudl.shape:     ",schritte_dudl

        if os.path.isfile("equilibrium") != True:
            sys.exit("equilibrium missing")
        if os.path.isfile("HesseMatrix_sphinx") != True:
            sys.exit("HesseMatrix_sphinx missing")
        if os.path.isfile("cell") != True:
            sys.exit("cell missing")
        cell = np.loadtxt("cell")
        equilibrium = np.loadtxt("equilibrium")

        if pot == 'l':
            allaa = np.array([2.2])  # this number does not matter for lennarj johnes
        for idxaa,aa in enumerate(allaa):

            evc0 = False
            ev0 = False
            for i in np.arange(schritte_listposfre):
                writeif = np.arange(0,schritte_listposfre,1)
                pos = getpos(i,positions,atoms)
                posandene = getpos(i,listposfre,atoms)
                evc = posandene[-1,-1]
                if i == 0:
                    evc0 = posandene[-1,-1]
                    ev0 = posandene[-1,-1]*1000./(atoms-1.)
                ev  = evc*1000./(atoms-1.)
                evd  = listdudl[i-1][4]
                evdh = listdudl[i-1][5]
                if i in writeif:
                    if i == 0:
                        continue
                    eh, ehc, fh = hesse.get_energy_forces(pot='h', coord_cart=pos, cell=cell, printresult=False)
                    checkh = np.sqrt((eh - evdh)**2)
                    check = np.sqrt((ev - ev0 - evd)**2)
                    checkdudl = (ev-ev0)*1000./(atoms-1.)
                    if check > 0.01:
                        sys.exit("check "+str(check))
                    if checkh > 0.01:
                        sys.exit("check "+str(check))
                print i,evd #eh,"||",evdh, evd,"||",ev-ev0,evd
                if i == 0:
                    continue


                #for idxaa,aa in enumerate(allaa):
                #print "aa:",aa
                ene_vs_eps = np.zeros(alleps.shape[0]+1)
                ene_vs_eps[:] = np.nan
                ene_vs_eps[0] = evd

                for idxeps,eps in enumerate(alleps):
                    if pot == "l":
                        el, elc, fl = hesse.get_energy_forces(pot='l', potparam=[ eps,NN], coord_cart=pos, cell=cell, printresult=False)
                    if pot == "m":
                        #el, elc, fl = hesse.get_energy_forces(pot='m', potparam=[ eps,aa, NN], coord_cart=pos, cell=cell, printresult=False, trans_1nn_from_h=True)

                        # hesse.py al -p mc1 -pp 0.221583 1.520757 2.920351 8.541780 0.196864 -ene -p2 mc1 0.447778 1.051045 2.920351 0.684528 -0.24431
                        #el, elc, fl = hesse.get_energy_forces(pot='mc1', potparam=[ 0.221583, 1.520757, 2.920351, 8.541780, 0.196864], coord_cart=pos, cell=cell, printresult=False, trans_1nn_from_h=True)

                        #el,elc,fl = hesse.get_energy_forces(
                        #                pot=args.pot,
                        #                pot2=args.pot2,
                        #                potparam=args.potparam,
                        #                #hessefile = args.inputfile,
                        #                #hessefile1nn = args.inputfile1nn,
                        #                coord_cart = pos,
                        #                cell = cell
                        #                #coordfile_cart = "cartesian_coords",
                        #                #add_transverse = args.tm,
                        #                #hrest1 = args.hrest1
                        #                )
                        elc = np.loadtxt("energy")

                        #print evd,eps,el
                    ene_vs_eps[idxeps+1] = el
                if type(ene_vs_eps_all) == bool:
                    ene_vs_eps_all = ene_vs_eps
                else:
                    ene_vs_eps_all = np.vstack((ene_vs_eps_all,ene_vs_eps))


                #if i in writeif:
                if pot == "m":
                    np.savetxt("ene_vs_eps_all_morse_aa_"+str(aa),ene_vs_eps_all,fmt='%.4f')
                    np.savetxt("ene_vs_eps_all_morse_aa_"+str(aa)+"eps",alleps,fmt='%.4f')
                if pot == "l":
                    np.savetxt("ene_vs_eps_all_lj_"+str(NN),ene_vs_eps_all,fmt='%.4f')
                    np.savetxt("ene_vs_eps_all_lj_"+str(NN)+"_eps",alleps,fmt='%.4f')
            # write ene_vs_es_all for every lambda
            #np.savetxt("ene_vs_eps_all_morse_aa_"+str(aa),ene_vs_eps_all,fmt='%.4f')
            #np.savetxt("ene_vs_eps_all_morse_aa_"+str(aa)+"eps",alleps,fmt='%.4f')







def c():
    makefolders = [ "lambda1.0_", "lambda0.85_", "lambda0.5_", "lambda0.15_", "lambda0.0_", "lambda0.4_", "lambda0.6_"]
    if args.l == '1.0':
        makefolders = [ 'lambda1.0_']
    if args.l == '0.0':
        makefolders = [ 'lambda0.0_']
    if args.l == '0.15':
        makefolders = [ 'lambda0.15_']
    if args.l == '0.85':
        makefolders = [ 'lambda0.85_']
    if args.l == '0.5':
        makefolders = [ 'lambda0.5_']

    print "makefoders:",makefolders
    hierbegin=os.getcwd()
    for makefolder in makefolders:
        os.chdir(hierbegin)
        folder_one = glob.glob(makefolder+"*")
        print "folder_one:",folder_one
        if len(folder_one) == 0:
            print "no folder with ",makefolder
            continue
        print ""
        print ""
        print "folder_one:",folder_one
        hier=os.getcwd()

        pos = None
        for i in folder_one:
            print ">> ",i," << extractpositions.sh -eq -p -r"
            if os.path.isfile(i+"/OUTCAR.gz") != True and os.path.isfile(i+"/outcar") != True:
                print "continue since no outcar in ",i
                os.chdir(hierbegin)
                continue
            os.chdir(hier)
            os.chdir(i)
            print "=== i:",i,"extractpositioins.sh"
            utils.run2("extractPOSITIONS.sh -eq -p -r",dont_raise_exceptino=True)
            print utils.printgreen("---- xtractPOSITIONS.sh -eq -p -r DONE")
            if os.path.isfile("POSITIONs") != True:
                print "POSITIONs was not created in",i
                os.chdir(hierbegin)
                continue
            utils.run2("wc -l POSITIONs")
            if os.path.isfile("equilibrium") == True:
                pos0=np.loadtxt("equilibrium")
                atoms = np.loadtxt("equilibrium").shape[0]
            else:
                print "no equilibrium file found, extractpositions.sh did not work"
                continue
            utils.run2("wc -l POSITIONs")
            if os.path.isfile(hier+"/eqpos") != True:
                utils.run2("cp "+hier+"/"+i+"/equilibrium "+str(hier)+"/eqpos")
            if os.path.isfile(hier+"/cell") != True:
                utils.run2("cp "+hier+"/"+i+"/cell "+str(hier)+"/cell")
            utils.run2("wc -l POSITIONs")
            print "------------- 1 -----------------"
            utils.run2("rm -f ka")
            utils.run2("tail -n+"+str(atoms*deletefirst+1)+" POSITIONs > ka")
            utils.run2("wc -l POSITIONs")
            postmp=np.loadtxt("ka")
            utils.run2("rm -f ka")
            print "postmp.spape:",postmp.shape[0]
            print "postmp.spape:",len(postmp.shape)
            if len(postmp.shape) == 1:
                if postmp.shape[0] == 0:
                    print "------------- 2 -----------------"
                    postmp=np.loadtxt("POSITIONs")

            schritte=np.shape(postmp)[0]/atoms
            # make sure to read the wohle file in case of just few schritten
            if schritte < 100:
                print "------------- 3 -----------------"
                postmp=np.loadtxt("POSITIONs")
                schritte=np.shape(postmp)[0]/atoms
                if schritte < 100:
                    if os.path.isfile("pre_positions_nojumps"):
                        postmp=np.loadtxt("pre_positions_nojumps")


            if postmp.shape[0] == 0:
                print "positions in ",i," is empty, .... continue"
                continue
            print "pos:",pos
            if pos == None:
                pos = postmp
            else:
                pos = np.concatenate((pos,postmp))
            print "np.shape(pow)",np.shape(pos)
            os.chdir(hier)


        print "--- pos  -------"
        print pos
        #print "--- pos0 -------"
        #print pos0
        print "--- len(np.shape(pos)):",len(np.shape(pos))
        print "--- np.shape(pos),atoms:",np.shape(pos),atoms
        schritte=np.shape(pos)[0]/atoms
        print "schritte:",schritte

        os.chdir(hierbegin)
        if makefolder == "lambda1.0_":
            np.savetxt("pos_1.0",pos,fmt="%.6f %.6f %.6f")
            print("written pos_1.0")
        if makefolder == "lambda0.0_":
            np.savetxt("pos_0.0",pos,fmt="%.6f %.6f %.6f")
            print("written pos_0.0")
        if makefolder == "lambda0.15_":
            np.savetxt("pos_0.15",pos,fmt="%.6f %.6f %.6f")
            print("written pos_0.15")
        if makefolder == "lambda0.85_":
            np.savetxt("pos_0.85",pos,fmt="%.6f %.6f %.6f")
            print("written pos_0.85")
        if makefolder == "lambda0.5_":
            np.savetxt("pos_0.5",pos,fmt="%.6f %.6f %.6f")
            print("written pos_0.5")
        if makefolder == "lambda0.4_":
            np.savetxt("pos_0.4",pos,fmt="%.6f %.6f %.6f")
            print("written pos_0.4")
        if makefolder == "lambda0.6_":
            np.savetxt("pos_0.6",pos,fmt="%.6f %.6f %.6f")
            print("written pos_0.6")


        posmapped=np.copy(pos)
        for i in np.arange(schritte):
            a = getpos(i,pos,atoms)-pos0
            #print "aa",a
            posmapped[atoms*i:atoms*i+atoms] = getpos(i,pos,atoms)-pos0



        if args.ci2d:
            print "getposmapped"
            print posmapped
            a = posmapped
            print "apply fcc symmmetry ... to get pos_{0,1}.0xy"
            xy   = np.transpose(np.array([a[:,0],a[:,1]]))
            xy_  = np.transpose(np.array([a[:,0],-a[:,1]]))
            xz   = np.transpose(np.array([a[:,0],a[:,2]]))
            xz_  = np.transpose(np.array([a[:,0],-a[:,2]]))
            yz   = np.transpose(np.array([a[:,1],a[:,2]]))
            yz_  = np.transpose(np.array([a[:,1],-a[:,2]]))

            x_y   = np.transpose(np.array([-a[:,0],a[:,1]]))
            x_y_  = np.transpose(np.array([-a[:,0],-a[:,1]]))
            x_z   = np.transpose(np.array([-a[:,0],a[:,2]]))
            x_z_  = np.transpose(np.array([-a[:,0],-a[:,2]]))
            y_z   = np.transpose(np.array([-a[:,1],a[:,2]]))
            y_z_  = np.transpose(np.array([-a[:,1],-a[:,2]]))

            posmappedsym = np.concatenate((xy,xy_,xz,xz_,yz,yz_,x_y,x_y_,x_z,x_z_,y_z,y_z_))
            #posmappedsym = np.concatenate((xy,xy_))
            if makefolder == "lambda1.0_":
                np.savetxt("pos_1.0xy",posmappedsym,fmt="%.5f %.5f")
                if args.ci2d:
                    utils.run2("getDOS2d.sh pos_1.0xy")
                    utils.run2("mv DOS2d DOS_from_pos_1.0xy")

            if makefolder == "lambda0.0_":
                np.savetxt("pos_0.0xy",posmappedsym,fmt="%.5f %.5f")
                if args.ci2d:
                    utils.run2("getDOS2d.sh pos_0.0xy")
                    utils.run2("mv DOS2d DOS_from_pos_0.0xy")


def corb(element=False):
    ''' element e.g.: al,ni,...'''
    import utils
    import my_atom as atom
    import shutil
    import os
    hier=os.getcwd()
    #########################################
    # get element
    #########################################
    #element_try = os.getcwd().split("/")
    ##print ">>>",element_try,len(element_try)
    #if len(element_try) >= 5:
    #    element = element_try[5]
    #else:
    #    element = 'al'
    #    print "ELEMENT",element," UNKNOWN! taking al instead for tmelt!"
    #if element == "anharmonic":
    #    element = 'al'
    #element = 'al'
    print utils.printred("element: "+element)
    ka = atom.atom([element])
    tmelt = ka.melting_rounded[0]
    print utils.printred( "tmelt: "+str(tmelt)+" K")



    #########################################
    # get pos0 and a, alat, sc
    #########################################
    if os.path.isfile("eqpos") == True:
        pos0=np.loadtxt("eqpos")
        atoms = pos0.shape[0]
    else:
        print utils.printred( "no eqpos found; run this script with -c option")
        sys.exit()
    if os.path.isfile("cell") == True:
        cell=np.loadtxt("cell")
    else:
        print utils.printred( "no file with name \"cell\" found; run this script with -c option")
        sys.exit()

    #########################################
    # lambdas
    #########################################
    lambdall = args.l
    if args.l == '1.0':
        lambdall = [ '1.0']
    if args.l == '0.0':
        lambdall = [ '0.0']

    print os.getcwd()
    f = os.getcwd()+'/test_loesch_'+str(element)
    if os.path.isdir(f) != True:
        os.makedirs(f)

    if os.path.isfile(f+"/POSITIONs"):
        os.remove(f+"/POSITIONs")
    if os.path.isfile(f+"/cell"):
        os.remove(f+"/cell")
    if os.path.isfile(f+"/equilibrium"):
        os.remove(f+"/equilibrium")
    if os.path.isfile(f+"/distances"):
        os.remove(f+"/distances")
    if os.path.isfile(f+"/distances_long"):
        os.remove(f+"/distances_long")
    if os.path.isfile(f+"/mapping"):
        os.remove(f+"/mapping")
    print "---empty--"

    for lambd in lambdall:
        os.chdir(hier)
        mapping = getmapping(atoms,pos0,eqposfile='eqpos',cellfile='cell')
        np.savetxt(f+"/mapping",mapping,fmt="%.0f")
        #shutil.copyfile("/data/glensk/mapping",f+"/mapping")
        if os.path.isdir(f) != True:
            os.makedirs(f)

        if os.path.isfile("pos_"+lambd) != True:
            print "no pos_"+lambd,"file"
            continue
        shutil.copyfile("eqpos",f+"/equilibrium")
        shutil.copyfile("cell",f+"/cell")
        shutil.copyfile("pos_"+lambd,f+"/POSITIONs")
        #print "---ssh-- at cmmd002 get_d_from_POSITIONs"
        #utils.run2("ssh cmmd002 \"cd /home/glensk/test_loesch_"+str(element)+";/home/grabowski/Thermodynamics/fortran/get_d_from_POSITIONs.x\"")
        utils.run2("cd "+f+";/Users/glensk/Thermodynamics/fortran/get_d_from_POSITIONs.x")


        print "---distdone--"
        print "dfn:",f+'/distances'
        print "dln:",f+'/distances_long'
        dfn = np.loadtxt(f+"/distances")
        dln = np.loadtxt(f+"/distances_long")
        makedos("atoms_1nn_all_f", lambd, dfn, tmelt, None, "dfn")
        makedos("atoms_1nn_all_f", lambd, dln, tmelt, None, "dln")

        if os.path.isfile(f+"/POSITIONs"):
            os.remove(f+"/POSITIONs")
        if os.path.isfile(f+"/cell"):
            os.remove(f+"/cell")
        if os.path.isfile(f+"/equilibrium"):
            os.remove(f+"/equilibrium")
        if os.path.isfile(f+"/distances"):
            os.remove(f+"/distances")
        if os.path.isfile(f+"/distances_long"):
            os.remove(f+"/distances_long")
        if os.path.isfile(f+"/mapping"):
            os.remove(f+"/mapping")




def corpyslow():
    neighbor = 1   # for 1st neares neighbor
    neighbor = 2   # for 2nd neares neighbor
    #########################################
    # get element
    #########################################
    import utils
    import my_atom as atom
    #element_try = os.getcwd().split("/")
    ##print ">>>",element_try,len(element_try)
    #if len(element_try) >= 5:
    #    element = element_try[5]
    #else:
    #    element = 'al'
    #    print "ELEMENT",element," UNKNOWN! taking al instead for tmelt!"
    #utils.printred( "element: "+element)
    #ka = atom.atom([element])
    #tmelt = ka.melting_rounded[0]
    #element = 'al'
    tmelt = '934'
    print utils.printred( "tmelt: "+str(tmelt)+" K")



    #########################################
    # get pos0 and a, alat, sc
    #########################################
    if os.path.isfile("eqpos") == True:
        pos0=np.loadtxt("eqpos")
        atoms = pos0.shape[0]
    else:
        print utils.printred( "no eqpos found; run this script with -c option")
        sys.exit()
    if os.path.isfile("cell") == True:
        cell=np.loadtxt("cell")
    else:
        print utils.printred( "no file with name \"cell\" found; run this script with -c option")
        sys.exit()



    #########################################
    # lambdas
    #########################################
    lambdall = args.l
    if args.l == '1.0':
        lambdall = [ '1.0']
    if args.l == '0.0':
        lambdall = [ '0.0']



    #########################################
    # try new way by going over every atom
    #########################################
    import crystal_generator
    import time
    reload(crystal_generator)
    crystal = crystal_generator.crystal()
    crystal0 = crystal_generator.crystal()
    crystal.cellvec = cell
    crystal0.cellvec = cell
    start=os.getcwd()

    nnlist = np.zeros((atoms,atoms))
    nnlist[:] = np.nan
    for idxi,pos0i in enumerate(pos0):
        nn1idx_doublecount = crystal.get_NNlist(idxi,neighbor,coordfile_cart='eqpos',cellfile='cell')
        #print nn1idx_doublecount
        nn1idx = nn1idx_doublecount[np.nonzero(nn1idx_doublecount > idxi)[0]]
        for jidx,j in enumerate(nn1idx):
            nnlist[idxi,jidx] = nn1idx[jidx]
        #print idxi,pos0i,nn1idx




    for lambd in lambdall:
        os.chdir(start)
        print "lambda:",lambd," loading positions"
        pos = np.loadtxt("pos_"+lambd)
        print "lambda:",lambd," positions loaded"
        schritte = pos.shape[0]/atoms
        poscurrent0 = pos0
        crystal0.rcar = poscurrent0
        df_array = np.array([])
        dl_array = np.array([])
        dt_array = np.array([])
        dfn_array = np.array([])
        dln_array = np.array([])
        dtn_array = np.array([])

        for i in np.arange(schritte):
            if i == 200:
                print 200
            if i == 450:
                print 450
            if i == 800:
                print 800
            for anz in np.arange(0,schritte,schritte/20):
                if i == anz: print i,"of",schritte

            poscurrent = getpos(i,pos,atoms)
            crystal.rcar = poscurrent
            #print "-- 0.poscurrent: --"
            #print crystal.rcar
            crystal.update_rrel_from_rcar()
            crystal.update_xyzcar_from_rcar()
            crystal.update_xyzrel_from_xyzcar()
            crystal0.update_rrel_from_rcar()
            crystal0.update_xyzcar_from_rcar()
            crystal0.update_xyzrel_from_xyzcar()
            #print "zuvor"
            # go through every atom
            for idxi in np.arange(atoms):
                nn1idx = nnlist[idxi][~np.isnan(nnlist[idxi])]
                #print idxi,nn1idx
                start = time.time()
                crystal.translate_atoms_cart(crystal.rcar[idxi])
                crystal.center_atoms_rel()
                crystal0.translate_atoms_cart(crystal0.rcar[idxi])
                crystal0.center_atoms_rel()
                end = time.time()
                #print "crystal:",end - start
                start = time.time()
                for idxj in nn1idx:

                    df = crystal.rcar[idxj]
                    dfn = np.linalg.norm(df)

                    d0 = crystal0.rcar[idxj]
                    d0n = np.linalg.norm(d0)

                    dl = utils.project_vector(df,d0)
                    dln = np.linalg.norm(dl)

                    dt = df - dl
                    dtn = np.linalg.norm(dt)

                    ###########################################
                    # concatenate values
                    ###########################################
                    if len(dfn_array) == 0:
                        dfn_array = np.array([dfn])
                    else:
                        dfn_array = np.vstack((dfn_array,np.array([dfn])))

                    if len(dln_array) == 0:
                        dln_array = np.array([dln])
                    else:
                        dln_array = np.vstack((dln_array,np.array([dln])))

                    if len(dtn_array) == 0:
                        dtn_array = np.array([dtn])
                    else:
                        dtn_array = np.vstack((dtn_array,np.array([dtn])))


                    if len(df_array) == 0:
                        df_array = df
                    else:
                        df_array = np.vstack((df_array,df))

                    if len(dl_array) == 0:
                        dl_array = dl
                    else:
                        dl_array = np.vstack((dl_array,dl))

                    if len(dt_array) == 0:
                        dt_array = dt
                    else:
                        dt_array = np.vstack((dt_array,dt))
                    #print df,d0,dfn,d0n
            #print "danach"
        makedos("atoms_"+str(neighbor)+"nn_all", lambd, dfn_array, tmelt, None, lambd)
        makedos("atoms_"+str(neighbor)+"nn_all", lambd, dln_array, tmelt, None, lambd)
        makedos("atoms_"+str(neighbor)+"nn_all", lambd, dtn_array, tmelt, None, lambd)

    sys.exit()




def atom_8_16_24():
    #for i in pos0:
    #    print i
    #print ""
    a = pos0[1,2]
    alat = a
    if pos0.size == 96:
        sc = np.array([[2*a,0,0],[0,2*a,0],[0,0,2*a]])
        eckatome = np.array([
       [ 0.  ,  0.  ,  0.  ],
       [ 0.  ,  0.  ,  a ],
       [ 0.  ,  a ,  0.  ],
       [ 0.  ,  a ,  a ],
       [ a ,  0.  ,  0.  ],
       [ a ,  0.  ,  a ],
       [ a ,  a ,  0.  ],
       [ a ,  a ,  a ]])
    else:
        sys.exit("change code if you got cell > 2x2x2 fcc")


    #########################################
    # get NN of ecka00toms
    #########################################
    Eidx = np.array([])
    Eidxy = np.array([])
    for idx,coord in enumerate(pos0):
        for e in eckatome:
            if e[0] == coord[0] and e[1] == coord[1] and e[2] == coord[2]:
                ide = idx
                Eidx = np.append(Eidx,idx)
                print ""
                print "!!! yes",coord,"ide:",ide
                #print "Eidx:",Eidx,Eidx.size,type(Eidx.size),type(4)
                check = 0
                for idy,j in enumerate(pos0):
                    ideNN = False
                    #print "-->",j[0],j[1],j[2],"||",e[0]+a/2.,e[1]+a/2.,e[2]+a/2.
                    tol = 1 # 5te nachkommastelle
                    if      nearly_equal(j[0],e[0]+a/2.,tol) and \
                            nearly_equal(j[1],e[1]+a/2.,tol) and \
                            nearly_equal(j[2],e[2],tol):
                                check = check+1
                                ideNN = idy
                                #print "NN:",j,"ideNN:",idy
                    if      nearly_equal(j[0],e[0]+a/2.,tol) and \
                            nearly_equal(j[1],e[1],tol) and \
                            nearly_equal(j[2],e[2]+a/2.,tol):
                                check = check+1
                                ideNN = idy
                                #print "NN:",j,"ideNN:",idy
                    if      nearly_equal(j[0],e[0],tol) and \
                            nearly_equal(j[1],e[1]+a/2.,tol) and \
                            nearly_equal(j[2],e[2]+a/2.,tol):
                                check = check+1
                                ideNN = idy
                                #print "NN:",j,"ideNN:",idy
                    ########################################################
                    # here and False and True
                    ########################################################
                    if ideNN != False: # and False:
                        ide = idx
                        idn = ideNN
                        d0 = pos0[idn]-pos0[ide]
                        d0_ = np.copy(d0)
                        for index,kkk in enumerate(d0):
                            #print "kkk:",kkk
                            if d0[index] == 0.0:
                                d0_[index] = d0.max()
                            else:
                                d0_[index] = 0.0
                            #print d0[index],d0_[index]


                        print "ide:",ide,"ideNN:",idy,j,d0,d0_
                        Eidxy = np.append(Eidxy,[ide,idy])

                if check != 3:
                    sys.exit("3 NN not found")
            if Eidx.size == int(args.e):
                break

    print "Eidx:",Eidx
    Eidxy = np.reshape(Eidxy,(Eidxy.size/2,-1))


    lambdall = args.l
    if args.l == '1.0':
        lambdall = [ '1.0']
    if args.l == '0.0':
        lambdall = [ '0.0']
    NNidxs = np.array([8])  # only 8 == [ 0.    2.05  2.05] atom
    NNidxs = np.array([8,24])  # only 8 == [ 0.    2.05  2.05] atom
    NNidxs = np.array([8,16,24])  # only 8 == [ 0.    2.05  2.05] atom

    ######################################
    # foldername
    ######################################
    foldername = "atoms"
    for i in Eidxy[:,1]:
        foldername = foldername+"_"+str(int(i))
    print "foldername:",foldername

    #######################################################
    # LOOP new
    # 0. mg und si ins bild holen fuer diss
    # 0. schauen was fuer mg-si1 noch nachgerechnet werden muss
    # 0. pt vac formation energy
    # 0. cu and pt and ag excess volume
    # 0. correct pd harmonic for diss
    # 1. make xmgrace for every element
    # 2. make 2nd NN
    # 3. mache xy picture fuer 2Ddos
    # 4. mache pot DOS 2NN fuer Si
    # 5. auswerten al-si dilute PW91
    # 6. auswerten mg-si dilute limit (and check effect chemical potential mg2si / si)
    #
    # plan:
    #   - checke if verschiebung nur durch den 1.5**2 und 0.5**2 kommt obwohl um 1 vibriert
    #       - np.linalg.norm(np.array([1.5,1.5,1.5])) = 2.598076211353316
    #       - np.linalg.norm(np.array([0.5,0.5,0.5])) = 0.8660254037844386
    #       - np.linalg.norm(np.array([1.,1.,1.])) = 1.7320508075688772
    #
    #       - (2.598076211353316+0.8660254037844386)/2 = 1.7320508075688772 is not 1 !!!
    #
    #   - second NN (can this be approximated harmonic?)
    #   - get correlation
    #   - get picture 2dDOS projected
    #   - show in which direction the atom goes towards:
    #           - inplane:      more towards 2NN
    #           - out of plane: more
    #   - make harmonic model with morse abhaengigkeit and start harmonic to morse on cmmd001
    #   - do MD with from Harmonic to Harmonic with morse 1NN
    #######################################################
    for lambd in lambdall:
        print "lambda:",lambd
        pos = np.loadtxt("pos_"+lambd)
        schritte = pos.shape[0]/atoms
        write_information = np.arange(0,schritte,schritte/10)
        if type(schritte) != int:
            sys.exit("schritte is not int")
        da = np.zeros((schritte*Eidxy.size/2,3))
        dan = np.zeros(schritte*Eidxy.size/2)

        di = np.zeros((schritte*Eidxy.size/2,3))
        din = np.zeros(schritte*Eidxy.size/2)

        do = np.zeros((schritte*Eidxy.size/2,3))
        don = np.zeros(schritte*Eidxy.size/2)

        #da = np.zeros((schritte*Eidxy.size/2,3))
        #dan = np.zeros(schritte*Eidxy.size/2)


        dc = np.zeros((schritte*Eidxy.size/2,3))
        senkr_vec = np.zeros((schritte*Eidxy.size/2,3))
        dcn = np.zeros(schritte*Eidxy.size/2)

        #dproj = np.zeros((schritte*Eidxy.size/2,3))
        #dprojn = np.zeros(schritte*Eidxy.size/2)
        #bn = np.zeros(schritte*Eidxy.size/2)
        #dcorr = np.empty((schritte,2))
        ##dcorr = np.array([])
        d0_length = np.linalg.norm(pos0[0]-pos0[8])  # das ist nur die laenge in einer ebene
        for n,eidxy in enumerate(Eidxy):
            ide = eidxy[0]
            idn = eidxy[1]
            d0 = pos0[idn]-pos0[ide]
            d0_ = np.copy(d0)
            for index,kkk in enumerate(d0):
                #print "kkk:",kk/k
                if d0[index] == 0.0:
                    d0_[index] = d0.max()
                else:
                    d0_[index] = 0.0
                #print d0[index],d0_[index]

            print "n,ide,idn:",n,ide,idn,"|| pos0[ide],po0[idn]:",pos0[ide],pos0[idn],d0,"schritte:",schritte
            for i in np.arange(schritte):
                da[schritte*n+i] = getpos(i,pos,atoms)[idn] - getpos(i,pos,atoms)[ide] # current vec but not mapped to first quadrant
                nnvec = da[schritte*n+i]
                dan[schritte*n+i] = np.linalg.norm(da[schritte*n+i]) #-d0_length
                if dan[schritte*n+i] == 0.0 and i == 0:
                    # in the pre_equilibration it may often happen thath the 1st struct is the equilibrium struct
                    continue
                if dan[schritte*n+i] == 0.0:
                    sys.exit("i!!!: "+str(i))

                di[schritte*n+i] = utils.project_vector(da[schritte*n+i],d0)
                din[schritte*n+i] = np.linalg.norm(di[schritte*n+i]) #-d0_length
                #print "din:",din[schritte*n+i]
                # eine projektion auf die ebene welche angeschaut wird
                #da[schritte*n+i] = utils.project_vector(da[schritte*n+i],d0_)
                #dan[schritte*n+i] = np.linalg.norm(da[schritte*n+i])

                #ct_vector(da[schritte*n+i],d0)
                don[schritte*n+i] = np.sqrt(abs(dan[schritte*n+i]+d0_length)**2.-abs(din[schritte*n+i]+d0_length)**2.)
                # is equal to:
                do[schritte*n+i] = da[schritte*n+i]-di[schritte*n+i]
                don[schritte*n+i] = np.linalg.norm(da[schritte*n+i]) #-di[schritte*n+i])
                #print "don:",don[schritte*n+i]

                # dcn equals to don
                #dc[schritte*n+i] = utils.project_vector(da[schritte*n+i],do[schritte*n+i])
                #print "dc:",dc[schritte*n+i]
                #dcn[schritte*n+i] = np.linalg.norm(dc[schritte*n+i])
                #print "dcn:",dcn[schritte*n+i]



def dosmax():
    import glob
    import utils
    import sys
    import crystal_generator
    import fah
    lam = [ '0.0', '1.0']
    mode = [ 'dfn', 'dln']

    #####################################################
    # high energies
    #####################################################
    fre_high = None
    lowpath = os.getcwd()
    print lowpath
    highpath = glob.glob(lowpath+"__high_*TAKE")
    if len(highpath) == 1:
        highpath = highpath[0]
    else:
        print("high path not found")
        highpath = "----------- does not exist or not found ------------"
        #sys.exit("high path not found")


    hierlow=os.getcwd()
    f=hierlow+"/auswertung_nn/"
    if os.path.isdir(f) != True:
        os.makedirs(f)
    dat = False
    for angkfolder in utils.lsn("*Ang_*K"):
        #print angkfolder
        angkfolder_high = highpath+"/"+angkfolder
        alat, temp = utils.string_to_num_list(angkfolder)
        NNdist = 0
        os.chdir(hierlow)
        os.chdir(angkfolder)
        angkfolder_low = os.getcwd()


        print "#########",alat, temp, angkfolder
        # get Fah_fre_low_best:
        fre_low = np.nan
        fre_high = np.nan
        fre_lowplushigh = None
        if os.path.isfile("Fah") == True:
            Fah = utils.readfile("Fah",split = True)
            for i in Fah:
                #print i
                for j in i:
                    #print j
                    if j == 'fre_best':
                        fre_low = float(i[1])

        if os.path.isfile(angkfolder_high+"/Fah_lowplushigh") == True:
            Fah_lowplushigh = utils.readfile(angkfolder_high+"/Fah_lowplushigh",split = True)
            for i in Fah_lowplushigh:
                #print i
                for j in i:
                    #print j
                    if j == 'fre_best':
                        fre_lowplushigh = float(i[1])
        else:
            fre_lowplushigh = np.nan



        dosfolder=glob.glob("atoms_1nn_all_f")
        if len(dosfolder) != 1:
            print("len dosfolder is not 1 in "+str(angkfolder_low))
            continue


        #get eqpos
        if os.path.isfile("eqpos") != True:
            sys.exit("eqpos missing in "+angkfolder_low)
        if os.path.isfile("cell") != True:
            sys.exit("cell missing in "+angkfolder_low)
        eqpos=np.loadtxt("eqpos")
        if NNdist == 0:
            crystal = crystal_generator.crystal()
            NNdist = crystal.get_NNlist(1,1,coordfile_cart="eqpos",cellfile="cell",return_NNdist = True)[1]
            #print alat,NNdist

        dosfolder=dosfolder[0]

        ##########################################
        # get avg_dudl
        ##########################################
        if os.path.isfile(dosfolder+"/avg_dUdL_d"):
            avg_dudl = np.loadtxt(dosfolder+"/avg_dUdL_d")
            fah_d = fah.get_dudlmeanfit(avg_dudl[:,0],avg_dudl[:,1])
        else:
            avg_dudl = np.nan
            fah_d = np.nan
        #print "fah_d:",fah_d
        #print "##############################################"




        os.chdir(dosfolder)
        #   von hier an sind wir im dosfolder! "atoms_8_16_24_9_17_25_10_18_26"

        for l in lam:
            for m in mode:
                if os.path.isfile("DOS_"+m+"_"+l):
                    dos=np.loadtxt("DOS_"+m+"_"+l)
                    dosmax = dos[:,1].max()
                    alat = os.getcwd().split("/")
                    for a in alat:
                        split = a.split("Ang_")
                        if len(split) == 2:
                            alat = float(split[0])
                            #print len(split)
                    dosxnr = np.nonzero(dosmax == dos[:,1])[0]
                    dosmax_x = dos[dosxnr][0,0]

                    if m == 'dfn':m_=1
                    if m == 'dln':m_=2
                    if m == 'dtn':m_=3
                    if l == '0.0':l_=0
                    if l == '1.0':l_=1
                    #print "alat,temp,l_:",alat,temp,l_,m_,dosmax_x,fre_low,NNdist
                    if type(dat) == bool:
                        dat =                np.array([alat,temp,l_,m_,    dosmax_x,fre_low,NNdist,fre_lowplushigh,fah_d])
                    else:                           #   0    1   2  3         4        5      6       7
                        dat = np.vstack((dat,np.array([alat,temp,l_,m_,    dosmax_x,fre_low,NNdist,fre_lowplushigh,fah_d])))
                #print dat
    os.chdir(hierlow)
    if False:
        print "--------------------------"
        print dat
        print "--------------------------"
    datdmin = np.concatenate((dat[:,4],dat[:,6])).min()
    datdmax = np.concatenate((dat[:,4],dat[:,6])).max()
    as_ = np.unique(dat[:,0])
    ts = np.unique(dat[:,1])
    ls = np.unique(dat[:,2])
    ms = np.unique(dat[:,3])
    for t in ts:  #temperatures
        for l in ls: #lambdas
            for m in ms:
                idxt = np.nonzero(dat[:,1] == t)[0]
                idxl = np.nonzero(dat[:,2] == l)[0]
                idxm = np.nonzero(dat[:,3] == m)[0]
                #print "idxt:",idxt
                #print "idxl:",idxl
                #print "idxm:",idxm
                idxu = np.intersect1d(idxt,idxl)
                #print "idxu:",idxu
                idxu = np.intersect1d(idxu,idxm)
                #print "idxu:",idxu

                datsave = dat[idxu][:,[6,4]] #6=NNdist 4=dosmax_x
                if m == 1:
                    np.savetxt(f+"/dfn_dosmax_vs_nn_"+str(l)+"_"+str(t),datsave)
                if m == 2:
                    np.savetxt(f+"/dln_dosmax_vs_nn_"+str(l)+"_"+str(t),datsave)
                if l == ls [-1]:  # halbierende
                    datsave_h = np.array([[datdmin,datdmin],[datdmax,datdmax]])
                    np.savetxt(f+"/dhh_dosmax_vs_nn_"+str(l)+"_",datsave_h)



    print "####################### Fah_d_diff"
    print "####################### Fah_d_diff"
    for a in as_:
        f_vs_d_m1 = False
        f_vs_d_m2 = False
        for t in ts:
            for m in ms:
                if False:
                    print a,t,m
                idxl1 = np.nonzero(dat[:,2] == 1.)[0]
                idxl0 = np.nonzero(dat[:,2] == 0.)[0]
                idxt = np.nonzero(dat[:,1] == t)[0]
                idxa = np.nonzero(dat[:,0] == a)[0]
                idxm = np.nonzero(dat[:,3] == m)[0] # dfn => n == 1
                idxu = np.intersect1d(idxt,idxa)
                idxu = np.intersect1d(idxu,idxm)
                #print idxu
                #print dat[idxu]
                #print "len:",len(dat[idxu])

                if len(dat[idxu]) == 0 or len(dat[idxu]) == 1:
                    continue
                #print dat[idxu]
                dfnmax0 = dat[idxu][0,4]
                dfnmax1 = dat[idxu][1,4]
                nndist = dat[idxu][1,6]
                d0_m_nn = dfnmax0-nndist
                d1_m_nn = dfnmax1-nndist
                dudl_nn = -(d1_m_nn + d0_m_nn)
                dd = dfnmax1 - dfnmax0
                flowplushigh = dat[idxu][0,7]
                flowplushigh = dat[idxu][0,5]
                if True:
                    if m == 1:
                        print a,nndist,t,d0_m_nn,d1_m_nn,dudl_nn #,-dd,flowplushigh

                if m == 1:  #dfn
                    if type(f_vs_d_m1) == bool:
                        f_vs_d_m1 =                np.array([dudl_nn,flowplushigh,a,t])
                    else:                           #   0    1   2  3         4        5      6       7
                        f_vs_d_m1 = np.vstack((f_vs_d_m1,np.array([dudl_nn,flowplushigh,a,t])))
                if m == 2:  #dfn
                    if type(f_vs_d_m2) == bool:
                        f_vs_d_m2 =                np.array([dudl_nn,flowplushigh,a,t])
                    else:                           #   0    1   2  3         4        5      6       7
                        f_vs_d_m2 = np.vstack((f_vs_d_m2,np.array([dudl_nn,flowplushigh,a,t])))
        print "------------"
        print f_vs_d_m1
        np.savetxt(f+"/f_vs_d_dfn_"+str(a)+"Ang",f_vs_d_m1)
        np.savetxt(f+"/f_vs_d_dln_"+str(a)+"Ang",f_vs_d_m2)
            #if m == 1:
            #    np.savetxt(f+"/f_vs_d_dfn_"+str(a)+"_"+str(t),np.array([[-d,fe]]))
            #if m == 2:
            #    np.savetxt(f+"/f_vs_d_dln_"+str(a)+"_"+str(t),np.array([[-d,fe]]))
    print "####################### Fah_d_lattice_constant"
    print "####################### Fah_d_lattice_constant"
    for a in as_:
        fah_d_ang = False
        for t in ts:
            idxt = np.nonzero(dat[:,1] == t)[0]
            idxa = np.nonzero(dat[:,0] == a)[0]
            idxm = np.nonzero(dat[:,3] == 1)[0]
            idxu = np.intersect1d(idxt,idxa)
            idxu = np.intersect1d(idxu,idxm)
            #print "dat[idxu]:",dat[idxu]
            if len(dat[idxu]) == 0:
                continue
            fah_d = dat[idxu][0][8]
            if type(fah_d_ang) == bool:
                fah_d_ang =                np.array([t,fah_d])
            else:                           #   0    1   2  3         4        5      6       7
                fah_d_ang = np.vstack((fah_d_ang,np.array([t,fah_d])))
        print t,"K  ",fah_d_ang
        np.savetxt(f+"/Fah_d_"+str(a)+"Ang",fah_d_ang)

    print "####################### Fah_d_temperatrue"
    print "####################### Fah_d_temperatrue"
    for t in ts:
        fah_d_k = False
        for a in as_:
            idxt = np.nonzero(dat[:,1] == t)[0]
            idxa = np.nonzero(dat[:,0] == a)[0]
            idxm = np.nonzero(dat[:,3] == 1)[0]
            idxu = np.intersect1d(idxt,idxa)
            idxu = np.intersect1d(idxu,idxm)
            #print idxu
            #print dat[idxu]
            if len(dat[idxu]) == 0:
                continue
            fah_d = dat[idxu][0][8]
            if type(fah_d_k) == bool:
                fah_d_k =                np.array([a,fah_d])
            else:                           #   0    1   2  3         4        5      6       7
                fah_d_k = np.vstack((fah_d_k,np.array([a,fah_d])))
        print a,"Ang  ",fah_d_k
        np.savetxt(f+"/Fah_d_"+str(t)+"K",fah_d_k)







def inte():
    import glob
    import os
    import sys
    from scipy.interpolate import interp1d


    folder = "atoms_1nn_all_f"
    if os.path.isdir(folder) != True:
        sys.exit(folder+" does not exist")
    mode = [ 'dfn', 'dln' ]
    m = 'dfn'
    ls = [ '0.0', '0.15', '0.4', '0.5', '0.6','0.85', '1.0' ]
    pot0 = np.loadtxt(folder+"/pot_"+m+"_0.0")
    pot1 = np.loadtxt(folder+"/pot_"+m+"_1.0")
    pot0xmin = pot0[:,0].min()
    pot1xmin = pot1[:,0].min()
    potmin = np.array([pot0xmin,pot1xmin]).max()
    pot0xmax = pot0[:,0].max()
    pot1xmax = pot1[:,0].max()
    potmax = np.array([pot0xmax,pot1xmax]).min()
    x = pot1[:,0]
    y = pot1[:,1]
    f = interp1d(x, y)
    f2 = interp1d(x, y, kind='cubic')

    xvaluesidx = np.nonzero( pot0[:,0] <= potmax)[0]
    xvaluesidx = np.nonzero( pot0[:,0][xvaluesidx] >= potmin)[0]
    xvalues = pot0[:,0][xvaluesidx]
    pot0y = pot0[:,1][xvaluesidx]


    xnew = xvalues
    import matplotlib.pyplot as plt
    #plt.plot(x,y,'o',xnew,f(xnew),'-', xnew, f2(xnew),'--')
    #plt.legend(['data', 'linear', 'cubic'], loc='best')
    #plt.show()
    ################### dudl ############################
    dudl_ar = np.zeros((len(xvalues),3))
    for ind,x in enumerate(xvalues):
        dudl_ar[ind][0] = xvalues[ind]
        dudl_ar[ind][1] = pot0y[ind]
        dudl_ar[ind][2] = f2(x)

    dudly = dudl_ar[:,2] - dudl_ar[:,1]
    dudl = np.transpose((dudl_ar[:,0],dudly))
    np.savetxt(folder+"/dudl_is_pot0.0_minus_pot1.0",dudl)

    ################### dos ############################
    dudl_avg = np.zeros((len(ls),2))
    dudl_avg[:,1] = np.nan
    for indl,l in enumerate(ls):
        if os.path.isfile(folder+'/DOS_'+m+"_"+l) != True:
            #print folder+'/DOS_'+m+"_"+l, "does not exist, continue"
            continue
        dos = np.loadtxt(folder+'/DOS_'+m+"_"+l)
        #print l,"dos::::",dos[-1][0],dos[-1][1]
        #sys.exit()
        x = dos[:,0]
        y = dos[:,1]
        f2 = interp1d(x, y,bounds_error=False, fill_value = 0.0,  kind='cubic')
        dos_ar = np.zeros((len(xvalues),2))
        dudl_mal_dos_ar = np.zeros((len(xvalues),2))
        for ind,x in enumerate(xvalues):
            #print 'ind,x:',ind,x
            dos_ar[ind][0] = xvalues[ind]
            dos_ar[ind][1] = f2(x)
            dudl_mal_dos_ar[ind][0] = xvalues[ind]
            dudl_mal_dos_ar[ind][1] = dudl[ind][1]*dos_ar[ind][1]
        np.savetxt(folder+"/dos_"+m+"_"+l,dos_ar)
        np.savetxt(folder+"/dudl_mal_dos_"+m+"_"+l,dudl_mal_dos_ar)
        dudl_avg[indl][0] = float(l)
        dudl_avg[indl][1] = np.sum(dudl_mal_dos_ar[:,1])
        print "dudl_avg[indl]:",dudl_avg[indl]
    #print ":::",dudl_avg
    a=dudl_avg
    #print a[~np.isnan(a).any(axis=1)]
    np.savetxt(folder+"/avg_dUdL_d",a)



if __name__ == '__main__':
    from argparse import RawTextHelpFormatter
    p = argparse.ArgumentParser(description='''help string''', formatter_class=RawTextHelpFormatter)
    p.add_argument('-element', default=False, nargs='+',
            help="element")
    p.add_argument('-c', default=False, action='store_true',
            help='1) in folder *Ang_*K\n (this should rund on every system (also directly on mac))'
            'run this to create pos_0.0 pos_1.0 files from all lambda 0.0 and 1.0 folders')
    p.add_argument('-ci2d', default=False, action='store_true',
        help='make 2dDOS for pos_{0,1}.0xy (use this with the -c option)')
    p.add_argument('-corb', default=False, action='store_true',
        help='2) in same folder as 1) make atom_xxx folder including dos of 1st NN\n'
        'this needs get_d_from_POSITIONs.x currently with ssh therefore do at cmpc')
    p.add_argument('-atom_8_16_24', default=False, action='store_true',
        help='make atom_xxx folder including dos of 1st NN')
    p.add_argument('-corpyslow', default=False, action='store_true',
        help='make atom_xxx folder including dos of 1st NN')
    p.add_argument('-inte', default=False, action='store_true',
        help='integrate over d; has to be started in jobvorlage/2x2x2_folder and needs atoms_1nn_all_f folder')
    p.add_argument('-dudl', default=False, action='store_true',
        help='go over positions and get dudl')
    p.add_argument('-dudlposc', default=False, action='store_true',
        help='go over POSITIONS and get energy with calculate_energy_and_forces')
    p.add_argument('-sb', nargs='?', type=int, default=False,
        help='start with steps given')
    p.add_argument('-se', nargs='?', type=int, default=False,
        help='end with steps given')
    p.add_argument('-dudlposcm', default=False, action='store_true',
        help='go over POSITIONS and get longvec from undisplaced structure for every atom; to be used with -dudlposc')
    p.add_argument('-std', default=False, nargs='+',
        help='get std of dUdL file (e.g. DOS_POSITIONS_auswerten.py -std dUdL -f */lam*/)')

    p.add_argument('-dudlvseps', default=False, type=str,
        help='get dudl vs eps')
    p.add_argument('-f','--folder',nargs='+',
        help='Do evaluation in every given folder; You can use wildcards (*.{}...)', default=[os.getcwd()])



    p.add_argument('-p', '--pot',choices=['h', 'm', 'l', 'mc1', 'i'], default=False, help='define potential (h=hesse (needs Hessematrix with -i option), m=Morse, l=Lennard-Jones, mc1=Morsec1, i=inversepotential (Dario))')
    p.add_argument('-p2', '--pot2', default=False, nargs='+', help='define potential 2 for positive (attractive) longitudinal part e.g.: mc1 0.447778 1.051045 2.920351 0.684528 -0.24431 (parameters defined in potparam go to repulsive part')
    p.add_argument('-pp', '--potparam',nargs='+', default=False, help='define parameter for potential (is to be used with -p option)')




    p.add_argument('-l', choices=['1.0', '0.0', '0.15', '0.5', '0.85'],
       help='define which lambda to take', default=[ '0.0', '0.15', '0.5', '0.85', '1.0', '0.4', '0.6'])
    p.add_argument('-e', choices=['1','2','3','4','5','6','7','8'],
       help='go through how many eckatome?', default='3')
    p.add_argument('-dosmax', default=False, action='store_true',
        help='get dos maximum')
    p.add_argument('-v', '--verbose',
                help='verbose', action='store_true', default=False)
    args = p.parse_args()


    if args.std:
        dataarray = std()


    if args.dudlposc:
        dudlposc(dudlposcm = args.dudlposcm, sb = args.sb, se = args.se, verbose = args.verbose)

    if args.dudlvseps:
        dudlvseps()

    if args.dudl:
        dudl()

    if args.c:
        c()

    if args.corb:
        #print type(args.element),args.element
        corb(element=args.element[0])

    if args.corpyslow:
        corpyslow()

    if args.atom_8_16_24:
        atom_8_16_24()


    if args.dosmax:
        dosmax()

    if args.inte:
        inte()
