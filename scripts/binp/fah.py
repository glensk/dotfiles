#!/usr/bin/env python
# -*- coding: utf-8 -*-
from __future__ import print_function

from subprocess import check_output,call
import argparse
import os,shutil,random,time
from ase.io import read as ase_read
from ase.io import write as ase_write
import sys
import glob
import utils_rename as utils
import shutil
import numpy as np
import pylab
import myutils as my
from scipy.integrate import quad
from scipy.optimize import curve_fit
import hesse as h

np.set_printoptions(suppress=True)   # display arrays withou 000000
np.set_printoptions(precision=6)    # print only 6 digist after .

def _printred(var):
    ENDC = '\033[0m'
    red = '\033[31m'
    print(red + str(var) + ENDC)

def _printgreen(var):
    ENDC = '\033[0m'
    red = '\033[32m'
    print(red + str(var) + ENDC)

def find_surface_filename(sysexit = True):
    ''' '''
    possible_filename = [ "Fah", "Fah_b_*", "Fah_d_*", "Fah_*" ]
    for i in possible_filename:
        filecheck = glob.glob(i)
        if filecheck:
            if len(filecheck) == 1:
                surface_filename = filecheck[0]
                return surface_filename

    # if we came to this point we have to exit
    if sysexit:
        sys.exit("Did not found any Fah surface file in current folder")
    else:
        #print("Did not found any Fah surface file in current folder")
        return None

def get_dudlmeanfit(lambdas_array, energies_array, return_fit = False, return_parameters = False,verbose=False,savefit=False):
    ''' Fits tangence function to dudl vs lambda data.

        input:
        v: the volume you want to fit
        t: the temperature you want to fit
    How can we in general give n parameters to any function without having to type parameters[0], parameters[1], parameters[2], .....?
    '''
    lambdas_array = np.array(lambdas_array)
    def tanFunction(x, a0, a1, a2, a3):
        return -a0*np.tan(np.pi*((1-a1)*x+a2+0.5))+a3
    def tanFunctionguess(energies_array):
        return [ energies_array[0],0.5,0.3,-0.6]

    #lambdas_array = self.l
    #energies_array = self.dudlmean[v, t, :]

    if np.isnan(lambdas_array).any() == True or np.isnan(lambdas_array).any() == True:
        a = np.array([lambdas_array, energies_array])  # this still has nans
        print("a||:---",a)
        if np.isnan(a).any() == True:
            print("WHY:")
            lambdas_array, energies_array = np.ma.compress_rows(np.ma.fix_invalid(a.T)).T  #no nas
    #print "--> compressed:",lambdas_array, energies_array
    if len(energies_array) is 0: return np.nan

    parameters, evinet_covariance = \
        curve_fit(tanFunction, lambdas_array, energies_array, tanFunctionguess(energies_array), maxfev=30000)
    if verbose:
        print('energies_array:',energies_array)
        print('lambdas_array :',lambdas_array)
        print('parameters    :',parameters)

    if parameters[1] < 0:
        print(energies_array)
        _printred("parameter 2 is sammer 0 we might have a probleme in fitting")
        _printred("parameter 2 is sammer 0 we might have a probleme in fitting")
        _printred("parameter 2 is sammer 0 we might have a probleme in fitting")
        _printred("parameter 2 is sammer 0 we might have a probleme in fitting")
    if parameters[2] < 0:
        print(energies_array)
        _printred("parameter 3 is sammer 0 we might have a probleme in fitting")
        _printred("parameter 3 is sammer 0 we might have a probleme in fitting")
        _printred("parameter 3 is sammer 0 we might have a probleme in fitting")
        _printred("parameter 3 is sammer 0 we might have a probleme in fitting")

    # fit
    points = 101
    fit_x_min = 0
    fit_x_max = 1

    fit_x_values = np.linspace(fit_x_min, fit_x_max, points)
    fit_y_values = tanFunction(fit_x_values, *parameters)
    fitdat = np.transpose([fit_x_values,fit_y_values])
    if savefit:
        np.savetxt(savefit,fitdat)

    # delta to fit
    #print('lambdas_array',lambdas_array,type(lambdas_array))
    fitdeltas = tanFunction(lambdas_array, *parameters) - energies_array
    #print('fitdeltas1',fitdeltas)

    # delta to fit max
    fitdeltasmax = abs(fitdeltas).max()

    fah = quad(tanFunction, 0, 1, args=(parameters[0],parameters[1],parameters[2],parameters[3]))[0]
    if return_fit == True:
        return fah, [fit_x_values,fit_y_values]
    if return_parameters == True:
        return fah, [fit_x_values,fit_y_values], parameters
    return fah #,fitdeltasmax

def write_Fah_surface(a, temps, dudlmeanfit, filename = 'Fah_surface'):
    print("a    :",a)
    print("temps:",temps)
    anz = np.sum(np.isnan(dudlmeanfit))
    out  = np.empty([a.size * temps.size - anz, 3])

    out[:]      = np.NAN
    i = 0
    for ia,a in enumerate(a):
        for it,t in enumerate(temps):
            #print ">",t,a,dudlmeanfit[ia,it]
            if np.isnan(dudlmeanfit[ia,it]) == True:
                continue
            out[i][0] = t
            out[i][1] = a
            out[i][2] = np.around(dudlmeanfit[ia,it], decimals=2)
            i = i+1
    print("filename:",filename)
    print(out)
    np.savetxt(filename,out,fmt="%.1f %.7f %.3f")
    return


def get_dudl_per_atom_min1_from_energies_eV_cell(
        energies_lambda_0_eV_cell=False,
        energies_lambda_1_eV_cell=False,
        number_of_atoms=0,
        align_lambda_0_first_to_lambda_1_yes_no=True,
        verbose=False):
    '''
    dudl =  my.get_dudl_from_energies_eV_cell(
            energies_lambda_0_eV_cell=ene_har_eV_cell,
            energies_lambda_1_eV_cell=ene_har_eV_cell,
            number_of_atoms=nat)
    '''
    if verbose:
        print('energies_lambda_1_eV_cell')
        print(energies_lambda_1_eV_cell)
        print('energies_lambda_1_eV_cell')
        print(energies_lambda_1_eV_cell)
        print('number_of_atoms')
        print(number_of_atoms)
        print('energies_lambda_0_eV_cell[0]')
        print(energies_lambda_0_eV_cell[0])
        print('energies_lambda_1_eV_cell[0]')
        print(energies_lambda_1_eV_cell[0])

    if align_lambda_0_first_to_lambda_1_yes_no == True:
        d = energies_lambda_1_eV_cell[0] - energies_lambda_0_eV_cell[0]
        if verbose:
            print('d',d)
        energies_lambda_0_eV_cell = energies_lambda_0_eV_cell + d
        if verbose:
            print('AE energies_lambda_0_eV_cell[0]')
            print(energies_lambda_0_eV_cell[0])
            print('AE energies_lambda_1_eV_cell[0]')
            print(energies_lambda_1_eV_cell[0])
    #return (a[:60,2]-a[:60,1])*27.211386/31*1000
    dudl = (energies_lambda_1_eV_cell - energies_lambda_0_eV_cell)/(number_of_atoms-1)*1000.
    dudl = dudl[~np.isnan(dudl)]
    if verbose:
        print('dudl',type(dudl))
        print(dudl)
    return dudl

def get_dudl_or_other_average_from_array(dudl):
    dudlav = np.copy(dudl)
    for idx,i in enumerate(dudlav):
        dudlav[idx] = (dudl[:idx+1]).mean()
    #print('dudlav')
    #print(dudlav)
    return dudlav

###################################
######### job creation ############
###################################
def fah_create_jobs_from_fqh(ace):
    if not os.path.isdir('fqh'):
        sys.exit('fqh does not exist')
    if os.path.isdir('fah'):
        sys.exit('fah does exist already')
    hesse_vol_pos = my.load_hessefiles_volumes_positionsfiles_from_fqh_folder()
    print("###################################")
    print("# NOW creating anharmonic folders #")
    print("###################################")
    fqhfolder = get_fqh_folder()
    print('fqhfolder',fqhfolder)
    Tmax = get_Tmax_from_fqh_thermo_2nd(verbose=False)
    print('Tmax',Tmax)
    temperatures1 = False
    temperatures2 = False
    for i in np.arange(100,1000,100): # [ 100, 200, ..., 1000 ]
        temp_test = np.arange(200,Tmax+i,i)
        if len(temp_test) >= 5 and len(temp_test) <=10:
            temperatures = temp_test
            break

    # extend temperature range
    temperatures = np.concatenate((np.array([25, 50,100,150]),temperatures))

    #print('len(i)',i,temperatures,len(temperatures))
    #temperatures = [200, 400, 600, 800, 1000 ]
    lambdas = [ 0.0, 0.15, 0.5, 0.85, 1.0 ]
    steps = 40000
    steps = 9000
    steps = 90000
    print("temperatrues",temperatures)
    #temperatures = [400 ]
    for idx_vol,i in enumerate(hesse_vol_pos):
        for idx_temp,temperature in enumerate(temperatures):
            hessefile = i[0]
            volume  = i[1]
            posfile = i[2]
            print("->",hessefile,'vol',volume,'T:',temperature) #, posfile',posfile)
            ipi_thermodynamic_integration_from_fqh(ace,volume,temperature,hessefile,posfile,lambdas,steps)
    return

def ipi_thermodynamic_integration_from_fqh(ace,volume,temperature,hessefile,posfile,lambdas,steps):
    seeds=2
    for l in lambdas:
        for rn in np.arange(seeds):
            rand_nr = random.randint(1,99999)
            #rand_nr = '1234567'
            folder = os.getcwd()+"/fah/"+str(volume)+"_"+str(temperature)+"K/lambda"+str(l)+"_"+str(rand_nr)
            if os.path.isdir(folder):
                sys.exit(folder+" does already exist!")
            os.makedirs(folder)
            ipi_inp = os.environ['HOME']+"/Dropbox/Albert/scripts/dotfiles/scripts/i-pi-mc_scripts/ipi_input_thermodynamic_integration_template.xml"

            with my.cd(folder):
                print('hessefile',hessefile)
                print('to folder',folder)
                print()
                hessefile_basename = os.path.basename(hessefile)
                shutil.copy2(hessefile, folder)
                print('hfbn',hessefile_basename)
                ipi_inp_basename = 'input.xml'
                shutil.copy2(ipi_inp, folder+'/'+ipi_inp_basename)
                print('ipi_inp_basename',ipi_inp_basename)
                print('posfile',posfile)
                print('to folder',folder)
                print()
                pos_basename = os.path.basename(posfile)
                frame = ase_read(posfile)
                if False:   # this stuff is not used
                    ene = ace.ene(frame.copy())  # needs a copy here, otherwise the DFT energy is evaluated and not the NN energy
                    ene_hartree = ene*0.036749322
                    print('ene',ene,"eV")
                    print('ene_hartree',ene_hartree,"hartree")
                ase_write(folder+"/pos.ipi.xyz",frame,format='ipi')
                #ase_write(folder+"/pos.lmp",frame,format='lammps-data')

                ##
                if ace.pot.pottype in [ 'runner' , 'n2p2' ]:
                    #frame.write(folder+'/pos.lmp',format='lammps-runner')
                    # here I would need an ase object
                    ase_write(folder+'/pos.lmp',frame,format='lammps-runner',pot=ace.pot)
                elif ace.pot.pottype in [ 'eam', 'eam-alloy' ]:
#                    frame.write(folder+'/pos.lmp',format='lammps-data')
                    ase_write(folder+'/pos.lmp',frame,format='lammps-data')
                else:
                    sys.exit('44321 Error: dont know how to write this lammps file')

                lammps_in_file = 'in.lmp'
                ace.ipi = True

                #print('ace.ipi',ace.
                ace.pot_get_and_ase_lmp_cmd(kmc=False,temp=False,nsteps=2000000000,ffsocket='inet',address=False)
                my.lammps_write_inputfile(folder=folder,filename=lammps_in_file,positions='pos.lmp',ace=ace)


                #print('fp',frame.positions)
                #print('fp',frame.positions.flatten())
                ang_to_bohr = 1.8897261
                np.savetxt(folder+"/x_reference.data",frame.positions.flatten()*ang_to_bohr,newline=" ",fmt="%3.15f")
                nat = frame.get_number_of_atoms()
                my.sed(folder+'/'+ipi_inp_basename,'xxx123',str(steps))  # steps
                my.sed(folder+'/'+ipi_inp_basename,'hessian.data',hessefile_basename)
                my.sed(folder+'/'+ipi_inp_basename,'init.xyz','pos.ipi.xyz')
                my.sed(folder+'/'+ipi_inp_basename,'xxx600',str(temperature))
                my.sed(folder+'/'+ipi_inp_basename,'xxx1.0',str(l))
                my.sed(folder+'/'+ipi_inp_basename,'xxx0.0',str(1.-l))
                my.sed(folder+'/'+ipi_inp_basename,'96,96',str(nat*3)+","+str(nat*3))
                my.sed(folder+'/'+ipi_inp_basename,'32342',str(rand_nr))
                my.sed(folder+'/'+ipi_inp_basename,'xxxene',str(0))
                timestamp = 'md_ff_'+str(int(round(time.time() * 100000)))
                my.sed(folder+'/'+ipi_inp_basename,'md_ff',timestamp)
                my.sed(folder+'/'+lammps_in_file,'md_ff',timestamp)
                print(folder)
                print('next thing to do: put the socket in the in.lmp!')
                print('executing this with:')
                print('python $HOME/sources/ipi/bin/i-pi input.xml')
                print('~/Dropbox/Albert/scripts/dotfiles/scripts/executables/lmp_mac < in.lmp')
                print()
                print('gives simulation.ti which has in 3rd column the total energy for the nn; independent of how lambdas are defined.')
                print('3rd column in simulation.ti is also independent of v_reference> -2.89829817e+00')
    return


def go_through_all_fah_jobs_and_create_joblist():
    ''' can create joblist if some jobs are already finished '''
    get_into_fah_folder(verbose=False)
    hier = os.getcwd()
    fvta = glob.glob(os.getcwd()+"/*_*K")
    print('fvta',fvta)
    sj=0
    for fvt in fvta:
        os.chdir(fvt)
        print('os.getcwd() fvt',os.getcwd())
        flsa = glob.glob(os.getcwd()+"/lambda*_*")
        for fls in flsa:
            os.chdir(fls)
            print('os.getcwd() fls',os.getcwd())


            #########################################
            # define what to do
            #########################################
            submitjob = True

            if os.path.isfile('simulation.ti'):
                submitjob = False
                sim = np.loadtxt('simulation.ti')
                print('len simulation.ti:',len(sim))
                if len(sim) == 0:
                    submitjob = True
                #if len(sim) != 5001:
                #    submitjob = True
                #    print('len',len(sim))
            if submitjob:
                sj += 1
                print('submit!!',sj,"myutils.py -ef ipi_sart_job")
                #call(["myutils.py -ef ipi_sart_job"],shell=True)

                f = open(hier+"/joblist.dat", "a")
                f.write(os.getcwd()+"\n")
                f.close()

    submitfile = os.environ['dotfiles']+'/scripts/i-pi-mc_scripts/submit_ipi_ti_joblist.sh'
    if os.path.isfile(submitfile):
        print('submitfile',submitfile,'does exist')
        get_into_fah_folder(verbose=False)
        shutil.copy2(submitfile,os.getcwd())
        print('copied submitfile to fah folder; now start the job by sbatch...')
        if os.environ['myhost'] in ['helvetios']:
            print()
            print()
            print()
            print('now submitting jobs')
            call(["sbatch submit_ipi_ti_joblist.sh"],shell=True)
            print()
            print()
            print()
        else:
            print('submit job by sbatch submit_ipi_ti_joblist.sh')
    else:
        print('submitfile',submitfile)
        print('submitfile does not exist!')
        print('for 32 atoms and only 150 jobs do: sbatch -p debug -t 00:59:00 submit_ipi_ti_joblist.sh')
    return

def fah_submit_job_if_on_cluster():
    get_into_fah_folder(verbose=False)
    hier = os.getcwd()

###################################
######### job analysis ############
###################################
def get_fqh_folder(verbose=False):
    hier = os.getcwd()
    get_into_fah_folder(verbose=False)
    ahfolder = os.getcwd()
    os.chdir(hier)
    fqhfolder = ahfolder[:-4]+'/fqh'
    #print('fqhfolder',fqhfolder)
    return fqhfolder

def get_evinet_folder(verbose=False):
    hier = os.getcwd()
    get_into_fah_folder(verbose=False)
    ahfolder = os.getcwd()
    os.chdir(hier)
    evinetfolder = ahfolder[:-4]+'/evinet'
    #print('evinetfolder',evinetfolder)
    return evinetfolder

def get_Tmax_from_fqh_thermo_2nd(verbose=False):
    fqhfolder = get_fqh_folder(verbose=False)
    fqhsurface = fqhfolder+'/thermo_2nd/Fqh'
    a = np.loadtxt(fqhsurface)
    #print(a)
    #print(int(a[-1][0]))
    return int(a[-1][0])

def get_into_fah_folder(verbose=False):
    hier = os.getcwd()
    if verbose:
        print('xx1 os.getcwd()',os.getcwd())
        print('xx2 os.getcwd()',hier[-4:])
        print('xx3 os.getcwd()',hier[:-4])

    # in case already in fah folder
    if hier[-4:] == '/fah':
        #os.chdir(hier[:-4])
        return

    # in case one can directly enter fah folder
    if os.path.isdir(hier+'/fah'):
        os.chdir(hier+"/fah")
        return

    # in case further down in fah folder
    if '/fah/' in hier:
        if verbose:
            print('kk','/fah' in hier)
            print('kk',hier.split('/fah/')[0]+'/fah')
        os.chdir(hier.split('/fah/')[0]+'/fah')
        return

    hier2 = os.getcwd()
    #print('hier2',hier2)
    if not os.path.isdir(hier2+'/fah'):
        if os.path.isdir(hier2+'/fqh'):
            os.makedirs('fah')
            os.chdir(os.getcwd()+'/fah')
            return
    hier2 = os.getcwd()
    if not os.path.isdir(hier2+'/fah'):
        sys.exit('no fah folder found')
    return

def fah_go_through_all_angK_folder_and_exec_function(function=False):
    get_into_fah_folder(verbose=False)
    fvta = glob.glob(os.getcwd()+"/*_*K")
    for fvt in fvta:
        os.chdir(fvt)
        print('os.getcwd() fvt',os.getcwd())
        #########################################
        # define what to do
        #########################################
        #get_Fah_and_avg_dudl_in_one_ang_K_folder_from_ipi_job()
        function()
    return

def get_Fah_and_avg_dudl_in_one_ang_K_folder_from_ipi_job(savefit="avg_dudl_fit.dat",verbose=True):
    ''' savefit = False or filename like avg_dudl_fit.dat '''
    f = get_sorted_lambda_folder(verbose=False)
    hier = os.getcwd()
    # das hier sollte schon das average ueber verschiedene seeds sein ...
    l_ = []
    dudl_mean_ = []
    err_uncor_ = []
    dudl_std_ = []
    for i in f:
        os.chdir(hier)
        os.chdir(i)
        ## here should be the loop over all seeds
        try:
            lam, dudl_mean, std_err_mean_uncorr, dudl_std =   \
                        get_dudl_from_ipi_job(verbose=False)
        except IOError:
            print('nope')
            return # returns nothing and nothing is written if not all lambdas calculated
        l_.append(lam)
        dudl_mean_.append(dudl_mean)
        dudl_std_.append(dudl_std)
        err_uncor_.append(std_err_mean_uncorr)
        print(str(lam).ljust(7),str(dudl_mean).ljust(20),str(std_err_mean_uncorr).ljust(20),str(dudl_std).ljust(20))

    if verbose:
        print('tt l',l_)
        print('tt dudl_mean_:',dudl_mean_)
        print('tt err_uncor_:',err_uncor_)
        print('tt dudl_std_ :',dudl_std_)
    avg_dudl = np.zeros((len(l_),4))
    avg_dudl[:,0] = l_
    avg_dudl[:,1] = dudl_mean_
    avg_dudl[:,2] = err_uncor_
    avg_dudl[:,3] = dudl_std_
    os.chdir(hier)

    fit = get_dudlmeanfit(l_, dudl_mean_, return_fit = True,savefit=savefit)
    fah = fit[0]
    print('fah',fah)
    np.savetxt("avg_dudl",avg_dudl,fmt='   %.2f       %.3f          %.3f             %.3f',header="lambda dudl_mean(meV/at)  err_uncor(meV/at)  stdDev(meV/at)\n# averages to Fah (using tangens):"+str(fah)+" (meV/atom)")
    #print('fit[1]',fit[1])
    #self.dudlmeanfit[inda,indt] = fit[0]
    #self.dudlmeanfit_tanfit[inda,indt] = np.transpose(fit[1])
    if os.path.isfile('Fah'):
        os.remove('Fah')
    f = open('Fah', "a")
    f.write("tangens "+str(fah)+" "+str(np.array(err_uncor_).max())+" 0.0\n")
    f.close()
    return

def get_sorted_lambda_folder(verbose=False):
    f = glob.glob("lambda*_*")
    hier = os.getcwd()
    if verbose:
        print('hier',hier)
        for i in f:
            print(i)
    f = utils.list_sorted(f)
    return f

def fah_get_Fah_surface():
    get_into_fah_folder(verbose=False)
    fvta = glob.glob(os.getcwd()+"/*_*K")

    if os.path.isfile('Fah_surface'):
        os.remove('Fah_surface')
    for fvt in fvta:
        v_str = fvt.split('fah/')[-1].split("_")[0]
        t_str = fvt.split('fah/')[-1].split("_")[1][:-1]
        #dudl  = np.loadtxt(fvt+'/Fah')
        fah_content = my.grep(fvt+'/Fah','tangens')[0].split()
        fah_str = fah_content[1]
        err_str = fah_content[2]
        print('fvt',fvt,t_str,v_str,fah_str,err_str)

        f = open('Fah_surface', "a")
        f.write(t_str+" "+v_str+" "+fah_str+" "+err_str+"\n")
        f.close()
    return


def fit_fah_surface_lmfit2d(x=None,y=None,z=None,verbose=False):
    if x is None and y is None and z is None and os.path.isfile("Fah_surface"):
        print('reading Fah_surface')
        data = np.loadtxt('Fah_surface')[:,[0,1,2]]
        x = data[:,0];
        y = data[:,1];
        z = data[:,2];
    if verbose:
        print('x',x)
        print('y',y)
        print('z',z)
    if x.ndim != 1:
        raise ValueError("x must be 1-dim.")
    if y.ndim != 1:
        raise ValueError("y must be 1-dim.")
    if z.ndim != 1:
        raise ValueError("z must be 1-dim.")
    if x.size != y.size or x.size != z.size:
        print('x.size',x.size)
        print('y.size',y.size)
        print('z.size',z.size)
        raise ValueError("x, y, and z must have the same size.")
    data = np.zeros((len(x),3))
    data[:,0] = x
    data[:,1] = y
    data[:,2] = z

    # add points at T=0K having fah = 0 mev/atom
    if True:
        all1v = np.unique(data[:,1])
        for v in all1v:
            data = np.append(data, [[0, v, 0]], axis=0)
    #print('data')
    #print(data)
    #sys.exit()

    #def exp1(x,y):
    #    return x*0+1, x,       y,       x**2,       x**2*y,       x**2*y**2,       y**2,       x*y**2,       x*y

    #def exp1sw(x,y):
    #    return x*0+1, y,       x**2,       y**2,       x**2*y,       x**2*y**2,       x,       x*y**2,       x*y

    #def exp2(x,y):
    #    return x,  y,  x**2,       x**2*y,       x**2*y**2,       y**2,       x*y**2,       x*y

    def exp17coef(xy,c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12,c13,c14,c15,c16,c17):
        x = xy[:,0]
        y = xy[:,1]
        return c1 *x**0*y**1 + \
               c2 *x**0*y**2 + \
               c3 *x**1*y**0 + \
               c4 *x**1*y**1 + \
               c5 *x**1*y**2 + \
               c6 *x**2*y**0 + \
               c7 *x**2*y**1 + \
               c8 *x**2*y**2 + \
               c9 *x**3*y**0 + \
               c10*x**3*y**1 + \
               c11*x**3*y**2 + \
               c12*x**3*y**3 + \
               c13*x**4*y**0 + \
               c14*x**4*y**1 + \
               c15*x**4*y**2 + \
               c16*x**4*y**3 + \
               c17*x**4*y**4 + \
               0


    #################################
    #####  fit the surface ##########
    #################################
    from lmfit import Model,minimize
    func = exp17coef
    model = Model(func)
    parameter_names = model.param_names
    independent_variable = model.independent_vars
    #print('parameter_names',parameter_names)
    #print('independent_variable',independent_variable)
    #print('parameter_names',type(parameter_names))
    for i in parameter_names: model.set_param_hint(i ,value=1)
    result = model.fit(data[:, 2], xy=data[:, 0:2])
    print(result.fit_report())
    #print('vgl (obtained by fit)',result.best_values)
    coef_lmfit = []
    for i in parameter_names: coef_lmfit.append(result.best_values.get(i))
    print('coef_lmfit',coef_lmfit)
    #print()
    #print(result.best_fit)
    diffmax = np.abs(result.best_fit-data[:,2]).max()
    print('diffmax',diffmax)
    #print('888',data[0,0:2])

    #################################
    #####  print comparision curves #
    #################################

    #######################
    # a) save original data
    #######################
    folder = "Fah_surface_analyze"
    if os.path.isdir(folder):
        shutil.rmtree(folder)
    os.mkdir(folder)
    def save_v_lines_from2darray(usearray=None,string=None,folder=None):
        if usearray is None:
            sys.exit('provide usearray1 and string')
        if string is None:
            sys.exit('provide usearray2 and string')
        if folder is None:
            sys.exit('provide folder')
        if type(string) != str:
            sys.exit('has tobe a sring')
        all1t = np.unique(usearray[:,0])
        all1v = np.unique(usearray[:,1])
        for v in all1v:
            ka = usearray[np.where(usearray[:,1]==v)[0]][:,[0,2]]
            arr = ka[ka[:,0].argsort()]
            np.savetxt(folder+"/out_ang3_vs_temp_"+string+"_"+str(v)+'.dat',arr)
        for t in all1t:
            ka = usearray[np.where(usearray[:,0]==t)[0]][:,[1,2]]
            arr = ka[ka[:,0].argsort()]
            #np.savetxt(folder+"/out_"+string+"_"+str(int(t))+'_K_vs_vol.dat',arr)
            np.savetxt(folder+"/out_temp_vs_ang3_"+string+"_"+str(int(t))+'.dat',arr)
        return
    save_v_lines_from2darray(usearray=data,string='data',folder=folder)

    #######################
    # b) save fitted data
    #######################
    denset = np.arange(0,x.max()+1,1)
    densev = np.arange(y.min(),y.max(),0.1)
    all1v = np.unique(data[:,1])
    all1t = np.unique(data[:,0])
    string="fit"
    for v in all1v:
        densev = np.repeat(v,len(denset))
        dd = np.array([denset,densev]).T
        arr = func(dd,*coef_lmfit)
        np.savetxt(folder+"/out_ang3_vs_temp_"+string+"_"+str(v)+'.dat',arr)
    for t in all1t:
        denset = np.repeat(t,len(densev))
        dd = np.array([densev,denset]).T
        arr = func(dd,*coef_lmfit)
        np.savetxt(folder+"/out_temp_vs_ang3_"+string+"_"+str(int(t))+'.dat',arr)
    return

def fah_get_thermo():
    get_into_fah_folder(verbose=False)
    if not os.path.isfile('Fah_surface'):
        sys.exit('no Fah_surface')
    Fah_surface = np.loadtxt('Fah_surface')[:,[0,1,2]]  # temp, vol, fah

    x = Fah_surface[:,0]
    y = Fah_surface[:,1]
    z = Fah_surface[:,2]
    fit_fah_surface_lmfit2d(x=x,y=y,z=z)

    #Tmax = Fah_surface[:,0].max()
    Tmax = get_Tmax_from_fqh_thermo_2nd(verbose=False)
    Vmin = Fah_surface[:,1].min()
    Vmax = Fah_surface[:,1].max()
    print('Tmax',Tmax)
    print('Vmin',Vmin)
    print('Vmax',Vmax)
    fqh = get_fqh_folder(verbose=False)
    evinet = get_evinet_folder(verbose=False)
    print('fqh',fqh)
    print('evinet',evinet)

    # interpolate the Fah_surface
    foldername='thermo'
    if os.path.isdir(foldername):
        shutil.rmtree(foldername)
    os.mkdir(foldername)
    os.chdir(foldername)
    fah_write_fit_input(Vmin,Vmax,Tmax,for_atoms=1,filename='fit.input')
    call(["cp ../Fah_surface Fah_surface"],shell=True)
    call(['[ "`command -v wolframscript`" = "" ] && module load mathematica; wolframscript -file fit.input'],shell=True)
    with open('Fah_surface_fit') as f:
        first_line = f.readline()
    print('first line')
    print(first_line)

    f = open('Fah_surface_fit','r')
    temp = f.read()
    f.close()
    os.remove("Fah_surface_fit")
    f = open('Fah_surface_fit', 'w')
    f.write("0"+first_line[1:])
    f.write(temp)
    f.close()

    #call(["wolframscript -file fit.input"],shell=True) # execute mathematica
    # helvetios: module load mathematica


    hier = os.getcwd()
    for i in ["1st","2nd","3rd"]:
        os.chdir(hier)
        os.mkdir("thermo_"+str(i))
        with my.cd("thermo_"+str(i)):
            print('now termo_'+str(i),os.getcwd())
            call(["cp "+evinet+"/EVinet_1 EVinet"],shell=True)
            call(["cp "+fqh+"/Fqh_*at_cell_per_atom_Surface_"+i+"_order__* Fqh"],shell=True)
            call(["cp ../Fah_surface_fit Fah"],shell=True)
            print('now gethermodynamcs.sh',os.getcwd())
            call(["$dotfiles/thermodynamics/getThermodynamics.sh"],shell=True)

    return

def fah_write_fit_input(Vmin,Vmax,Tmax,for_atoms=1,filename='fit.input'):
    Fahsurfaceminalat = Vmin
    Fahsurfacemaxalat = Vmax
    structureFactor = 1
    sc1 = 1
    tmelt = Tmax
    with open(filename, "w") as f:
        f.write('(* adjustable parameters start *)'+"\n")
        f.write("\n")
        f.write('FsurfFile = "Fah_surface";'+"\n")
        f.write('type = 2                           (*  1: T(K)  aLat(Ang/at)  F(meV/at)   *)'+"\n")
        f.write('                                   (*  2: T(K)  V(Ang/at)     F(meV/at)   *)'+"\n")
        f.write('                                   (*  3: T(K)  V(Ang/cell)   F(meV/cell) *)'+"\n")
        f.write("\n")
        f.write("\n")
        f.write('min = '+str(Fahsurfaceminalat)+";  (*  aLat or volume range (same format as Fsurf) for the 2. fit *)"+"\n")
        f.write('max = '+str(Fahsurfacemaxalat)+";  (*  typically: Vmin=Veq(T=0K) and Vmax=Veq(Tmelt) *)"+"\n")
        f.write("mesh = 100;                        (* 100 is good and should be kept *)"+"\n")
        f.write("\n");
        f.write("structureFactor = "+str(structureFactor)+";  (*  4: fcc  2: bcc  1: supercells *)"+"\n")
        f.write("sc = "+str(sc1)+";                           (*  supercell, 1 for bulk *)"+"\n")
        f.write("nAtoms = "+str(for_atoms)+";                     (*  1 for bulk *)"+"\n")
        f.write("\n")
        f.write("fitType = \"Fvib\";     (*  \"Fvib\"  or  \"poly\"  fit type for 1. fit; take \"Fvib\" for Fah or Fel *)"+"\n")
        #f.write("basis[V_, T_] := {1,T, V }      (*  for Fah typically: \"Fvib\" and {1,T,V} *)"+"\n")
        f.write("basis[V_, T_] := {1,T,T^2,T^3,V,V^2,V^3}      (*  for Fah typically: \"Fvib\" and {1,T,V} *)"+"\n")
        f.write("\n")
        f.write("                                (*  for Fel typically: \"Fvib\" and {1,T, V,T V,T^2,V^2,V^3} *)"+"\n")
        f.write("basis2[V_]:={1, V, V^2, V^3}    (*  should be more than sufficient: {1,V,V^2,V^3} *)"+"\n")
        f.write("\n")
        f.write("minT = 1;                       (* typically 1     *)"+"\n")
        f.write("maxT = "+str(tmelt)+";          (*           Tmelt *)"+"\n")
        f.write("stepT = 1;                      (*           2     *)"+"\n")
        f.write("\n")
        f.write("useMeanFreqs=False;              (* if True \"mean_freqs\" file must be available; format as Fsurf, e.g. aLat(Ang) meanFreq(meV) *)"+"\n")
        f.write("                                (* meanFreqs are then used in the fit formula (check fitSurface.math) *)"+"\n")
        f.write("(* adjustable parameters end *)"+"\n")
        f.write("\n")
        f.write("(*<<\"~/scripts/Thermodynamics/fitSurface.math\"*)"+"\n")
        thermodynamics = os.environ['thermodynamics'];
        print('thermodynamics:',thermodynamics)
        f.write("<<\""+thermodynamics+"/mathematica/fitSurface.MATH\""+"\n")
        f.close()
    return

def get_dudl_from_ipi_job(verbose=True):
    ''' help '''
    #################################################
    # get the temperature / equilibration time
    #################################################
    desired_temp = int(my.grep('input.xml',"temperature units")[0].split()[2])
    lam = float(my.grep('input.xml',"force forcefield=.*weig")[1].split("'")[3])
    if verbose:
        print('lambda :',lam)
        print('desired:',desired_temp)

    out = np.loadtxt('simulation.out')
    t = out[:,3]
    #if verbose:
    #    print("<T>    K:",t.mean())
    #    print("<Tstd> K:",t.std())

    tall = np.where(np.logical_and(t>=desired_temp-t.std()/2, t<=desired_temp+t.std()/2))
    equilibration = tall[0][0]
    if verbose:
        print('equilibration',equilibration)
    t = t[equilibration:]
    if verbose:
        print("<T>    K:",t.mean())
        print("<Tstd> K:",t.std())

    #################################################
    # get dudl
    #################################################
    dudl = get_dudl_from_file_with_energies_lambda_0_1(filepath=False,number_of_atoms=False,verbose=False)
    dudl = dudl[equilibration:]
    #if False:
    #corr = np.correlate(dudl,dudl,mode='full')
    #print('corr',corr)
    #np.savetxt('corr.dat',corr)
    #print(os.getcwd())
    #sys.exit()
    std_err_mean        = dudl.std()/np.sqrt(len(dudl))
    std_err_mean_uncorr = dudl[::10].std()/np.sqrt(len(dudl[::10]))
    if verbose:
        print("<dudl>             (meV):",dudl.mean())
        print("stdDev             (meV):",dudl.std())
        print("std_err_mean       (meV):",std_err_mean)
        print("std_err_mean_uncorr(meV):",std_err_mean_uncorr,"needs to be improved")
    tempav = get_dudl_or_other_average_from_array(t)

    #################################################
    # get dudl std as function of time
    #################################################
    dudlstd = np.copy(dudl)
    for idx,i in enumerate(dudl):
        #print('idx',idx,dudl[:idx+1].std()/np.sqrt(idx+1))
        dudlstd[idx] = dudl[:idx+1].std()/np.sqrt(idx+1)
    #print('len(dudl)',len(dudl))
    #np.savetxt('dudlstd.dat',dudlstd)
    #print('t1',tempav)
    #np.savetxt('tempav.dat',tempav)
    #np.savetxt('temp.dat',t)
    return lam, dudl.mean(), std_err_mean_uncorr, dudl.std()

def get_dudl_from_file_with_energies_lambda_0_1(filepath=False,number_of_atoms=False,verbose=False):
    if filepath == False:
        if os.path.isfile('simulation.ti'): filepath = 'simulation.ti'
    if filepath == False:
        sys.exit('please specify the filepath!')

    #print('filepath',filepath)
    if number_of_atoms == False and filepath == 'simulation.ti':
        f = glob.glob("Hessematrix_*")
        if len(f) == 1:
            h = np.loadtxt(f[0])
            #print(h.shape)
            number_of_atoms = h.shape[0]/3
            #print(number_of_atoms)
        #print('f',f)
    if number_of_atoms == False:
        sys.exit('please specify the number of atoms')


    l01 = np.loadtxt(filepath)
    if verbose:
        print(l01.shape)
    l0_column = 0
    l0 = l01[:,l0_column]
    if verbose:
        print(l0)
        print(np.diff(l0))
    if len(np.unique(np.diff(l0))) == 1:
        l0_column = 1
        l0 = l01[:,l0_column]

    l1_column = l0_column+1
    l1 = l01[:,l1_column]

    #print('l0')
    #print(l0)
    #print('l1')
    #print(l1)
    l0 = my.convert_energy(l0,"hartree","eV",number_of_atoms,verbose=False)
    l1 = my.convert_energy(l1,"hartree","eV",number_of_atoms,verbose=False)
    #print('l0o')
    #print(l0)

    dudl = get_dudl_per_atom_min1_from_energies_eV_cell(
            energies_lambda_0_eV_cell=l0,
            energies_lambda_1_eV_cell=l1,
            number_of_atoms=number_of_atoms,
            align_lambda_0_first_to_lambda_1_yes_no=True)
    #dudlav = get_dudl_or_other_average_from_array(dudl)
    #print('dudl')
    #print(dudl)
    #print('dudlav')
    #print(dudlav)
    #np.savetxt('dudlav.dat',dudlav)
    #sys.exit()
    return dudl

def save_dudl_to_file(folder,dudl,filename="dudl"):
    if os.path.isifle(folder+"/"+filename):
        print(folder+"/dudl does already exist!;exit")
        sys.exit()
    #  step   time(fs)  temp(K) average       U(meV/at)    Uref          dUdL   average    offset
    dudl = np.zeros((len(dudl)-1,3))
    dudl[:,0] = np.arange(len(dudl))
    dudl[:,1] = dudl
    np.savetxt(folder+"/"+dudl)
    return


######### old fah stuff ############
class ah():
    ''' anharmonic contribution to Gibbs free energy '''
    def __init__(self):
        self._verbose       = False
        self._filestring    = None
        self._pwd           = os.getcwd()
        self._dudl_cutoff_first_structures      = 50
        self._dudl_correlation_length           = 10

        self.files          = None
        self.surface        = None
        self.surface_expr   = None
        self.surface_filename = None
        self.a              = None      # array([ 4.05, 4.09, 4.13])
        self.t              = None      # array([ 250.,  500.,  700.,  934.])
        self.l              = None      # array([ 0.  ,  0.15,  0.5 ,  0.85,  1.  ])
        self.s              = None      #
        self.e              = None      #
        self.dudl           = None      # --> same as dUdL(_noJumps) file
        self.dudlmean       = None      # --> see below (mean dUdL for every lambda and temp)
        self.dudlmeanfit    = None      # --> see below same as avg_dUdL_fre
        self.dudlste        = None      # mean dUdl ste for every lambda
        self.dudlstd        = None      # mean dUdl std for every lambda
        self.data           = None      #
        self.ah             = None      # anharmonic corrections at certain temp and vol
        self.ahste          = None      # anharmonic ste         at certain temp and vol
        self.ahstd          = None      # anharmonic std         at certain temp and vol

        # [158]f.dudlmean (f.dudlstd, f.dudlste)
        # Out[158]:
        # array([[[  0.51,   0.27,   0.01,  -0.53,  -0.64],     <- 250  <- 4.05
        #         [  2.88,   1.66,  -0.16,  -1.57,  -2.45],     <- 500  <- 4.05
        #         [  6.09,   3.81,   0.12,  -2.74,  -4.43],     <- 700  <- 4.05
        #         [ 10.79,   6.62,   0.58,  -5.46,  -8.57]],    <- 934  <- 4.05
        #
        #        [[  0.73,   0.6 ,  -0.09,  -0.5 ,  -0.92],     <- 250  <- 4.05
        #         [  3.76,   2.34,   0.02,  -1.92,  -2.71],     <- 500  <- 4.05
        #         [  7.06,   4.63,   0.62,  -3.43,  -5.23],     <- 700  <- 4.05
        #         [ 13.32,   8.45,   1.25,  -5.57, -10.17]]     <- 934  <- 4.05
        #
        #        [[  1.28,   0.61,   0.04,  -0.57,  -0.83],     <- 250  <- 4.05
        #         [  4.88,   3.73,   0.36,  -2.06,  -3.07],     <- 500  <- 4.05
        #         [  9.27,   5.84,   1.09,  -3.36,  -5.84],     <- 700  <- 4.05
        #         [ 17.44,  11.05,   1.96,  -5.83, -11.71]]])   <- 934  <- 4.05

        # [163]f.dudlmeanfit
        # Out[163]:
        # array([[-0.064615, -0.02395 ,  0.390583,  0.621287],          <-4.05
        #        [-0.019864,  0.158969,  0.634918,  1.375083],          <-4.05
        #        [ 0.04227 ,  0.647443,  1.216228,  2.375211]])         <-4.05


        return

    def initialize_files_data(self):
        """ defines:    self.files
                        self.a
                        self.l
                        self.s
                        self.data filled with NaNs """
        print("initialize_files_data ...")
        self.files = utils.lsn(self._filestring)
        self.filesinfo = []

        atls = np.array([ utils.string_to_num_list(string) for string in self.files ])
        if self._verbose:
            print("atls:",atls)
            print("atls.shape:",atls.shape)
            print("atls.shape[1]:",atls.shape[1])

        # wir setzen a(alts) und t(temperatures) voraus, alles andere optional
        # oder kann spaeter hinzugefuegt werden
        self.a      = np.unique(atls[:,0])
        self.t      = np.unique(atls[:,1])
        if atls.shape[1] <= 2:
            return
        if atls.shape[1] > 2:
            self.l      = np.unique(atls[:,2])

        # determine maximun amount of seeds (seeds are not distinguishable)
        nr_s = 0
        for a in self.a:
            array_a = atls[np.nonzero(atls[:,0] == a)[0]]
            #print "a:",array_a
            for t in self.t:
                array_t = array_a[np.nonzero(array_a[:,1] == t)[0]]
                if len(array_t) is 0: continue # in case a temp/vol missing
                #print "t:",array_t,"len_t:",len(array_t)
                for l in self.l:
                    array_l = array_t[np.nonzero(array_t[:,2] == l)[0]]
                    if len(array_l) > nr_s: nr_s= len(array_l)

        self.s      = np.arange(nr_s)
        self.data   = np.empty((
            len(self.a),
            len(self.t),
            len(self.l),
            len(self.s)), dtype=object)  # necessary for arrays of unequal length
        self.data[:] = np.NAN
        return

    def import_dudl_files_data(self):
        """ Imports all files to self.data; Working with this data (e.g.
        extracting dudls) should be separated; Doing so this can be used for
        sphinx and vasp """
        def get_inds(a,t,l):
            for s in self.s:
                #print "a,t,l,s:",a,t,l,s
                if np.isnan(np.any(self.data[a,t,l,s])):
                    return s

        for file in self.files:
            dudl = pylab.loadtxt(file)
            atls = utils.string_to_num_list(file)
            a = atls[0]
            t = atls[1]
            l = atls[2]
            s = atls[3]

            inda = np.nonzero(self.a == a)[0][0]  # index
            indt = np.nonzero(self.t == t)[0][0]  # index
            indl = np.nonzero(self.l == l)[0][0]  # index
            inds = get_inds(inda,indt,indl)
            print("file:",file,a,t,l,s,"(ind:",inda,indt,indl,inds,")")
            self.filesinfo.append([file,"ind:",[inda,indt,indl,inds]])
            self.data[inda, indt, indl, inds] = dudl
        return

    def import_avg_dudl_data(self, avg_dudl_filenames = "*Ang_*K/avg_dUdL_lowplushigh_fre", verbose = False):
        '''
        can import:
         - "*Ang_*K/avg_dUdL_lowplushigh_fre"   (in high folder)
         - "*Ang_*K/avg_dUdL_fre"               (in low  folder) ????

        '''
        #################################################################################
        # get alat, temperatures, lambdas
        #################################################################################
        avg_dudl_files = utils.lsn(avg_dudl_filenames)
        self.atls = np.array([ utils.string_to_num_list(string) for string in avg_dudl_files ])
        self.a      = np.unique(self.atls[:,0])
        self.t      = np.unique(self.atls[:,1])

        self.l = np.array([])  # here we want to load all data! not exclude all lambda 0.5 if one lambda 0.5 is missing
        for i in avg_dudl_files:
            #print "i:",i
            avg_dudl = np.loadtxt(i)
            #print "avg_dudl:",avg_dudl
            #print "avg_dudl.shape:",avg_dudl.shape
            if len(avg_dudl.shape) == 1:
                continue
            #print "avg_dudl.shape:",len(avg_dudl.shape)
            l_tmp = avg_dudl[:,0]
            #print "l_tmp:",l_tmp
            for ll in l_tmp:
                self.l = np.append(self.l,ll)
            self.l      = np.unique(self.l)
            #print avg_dudl
            #print avg_dudl.shape
            #print "l:",l
        if verbose:
            print("a:",self.a)
            print("t:",self.t)
            print("l:",self.l)
        #if len(self.l) == 9:
        #    print type(self.l)
        #    self.l = np.array([0.0, 0.15, 0.5, 0.85, 1.0])

        self.dudlmean     = np.empty([self.a.size, self.t.size, self.l.size])
        self.dudllen      = np.empty([self.a.size, self.t.size, self.l.size])
        self.dudlste      = np.empty([self.a.size, self.t.size, self.l.size])
        self.dudlstd      = np.empty([self.a.size, self.t.size, self.l.size])
        self.dudlmeanlen  = np.empty([self.a.size, self.t.size])
        self.dudlmeanfit  = np.empty([self.a.size, self.t.size])
        self.dudlmeanfit_tanfit  = np.empty([self.a.size, self.t.size,101,2])

        self.dudlmean[:]     = np.NAN
        self.dudllen[:]      = np.NAN
        self.dudlste[:]      = np.NAN
        self.dudlstd[:]      = np.NAN
        self.dudlmeanlen[:]  = np.NAN
        self.dudlmeanfit[:]  = np.NAN

        #self.dudlavg = np.copy(self.dudl)

        np.set_printoptions(suppress=True)   # display arrays withou 000000
        np.set_printoptions(precision=6)    # print only 6 digist after .

        #################################################################################
        # got trhount avg_dUdL files and fit
        #################################################################################
        for i in avg_dudl_files:
            if verbose:
                print("------------------>i:",i)
            atls = utils.string_to_num_list(i)
            a = atls[0]
            t = atls[1]
            inda = np.nonzero(a == self.a)[0]
            indt = np.nonzero(t == self.t)[0]
            avg_dudl = np.loadtxt(i)
            #print alll
            #print avg_dudl,"|",len(avg_dudl),"|",avg_dudl.shape,len(avg_dudl.shape)
            if len(avg_dudl) == 0:
                continue # empty
            if len(avg_dudl.shape) == 1:
                avg_dudl = np.array([avg_dudl])

            for lline in avg_dudl:
                if verbose:
                    print("lline:",lline)
                l = lline[0]
                indl = np.nonzero(l == self.l)[0]
                if verbose:
                    print("a,inda:",a,inda)
                    print("t,indt:",t,indt)
                    print("l,indl:",l,indl)
                self.dudlmean[inda,indt,indl] = lline[1]
                self.dudlste[inda,indt,indl] = lline[2] #*2   # so we now have ste and not ste/2
                self.dudlstd[inda,indt,indl] = lline[3]

            ## get energy vs temperatrues data
            lambdas_array = avg_dudl[:,0]
            #energies_array = self.dudlmean[inda, indt, :]
            energies_array = avg_dudl[:,1]

            ## only make a fit if we have > 3 lambda values
            if len(lambdas_array) <= 3:
                if verbose:
                    print("only following lambdas:",lambdas_array," therefore no fit!")
                if len(lambdas_array) ==1:
                    if lambdas_array[0] == 0.5:
                        #print "yes!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
                        self.dudlmeanfit[inda,indt] = energies_array[0]
                ## but if only lambda 0.5 one could take this value instead...
                continue

            if verbose:
                print("La;",lambdas_array)
                print("LA;",energies_array)
            if self._verbose:
                print("i:",i,lambdas_array,type(lambdas_array),energies_array[0],type(energies_array))
            fit = get_dudlmeanfit(lambdas_array, energies_array, return_fit = True)
            self.dudlmeanfit[inda,indt] = fit[0]
            self.dudlmeanfit_tanfit[inda,indt] = np.transpose(fit[1])
        return


    def save_data(self):
        print("saving data ...")
        np.savez_compressed( "dudl",                \
                files       = self.files,           \
                data        = self.data,            \
                dudl        = self.dudl,            \
                a           = self.a,               \
                t           = self.t,               \
                l           = self.l,               \
                s           = self.s,               \
                dudlmean    = self.dudlmean,        \
                dudlmeanfit = self.dudlmeanfit,     \
                dudlste     = self.dudlste,         \
                dudlstd     = self.dudlstd,         \
                ah          = self.ah,              \
                ahste       = self.ahste,           \
                ahstd       = self.ahstd,           \
                _filestring = self._filestring,     \
                _pwd        = self._pwd             \
                )
        return


    def load_data(self):
        ''' '''
        print("loading data ...")
        var = np.load( "dudl" + ".npz" )

        self.files          = var['files']
        self.data           = var['data']
        self.a              = var['a']
        self.t              = var['t']
        self.l              = var['l']
        self.s              = var['s']
        self.dudlmean       = var['dudlmean']
        self.dudlmeanfit    = var['dudlmeanfit']
        self.dudlste        = var['dudlste']
        self.dudlstd        = var['dudlstd']
        self.ah             = var['ah']
        self.ahste          = var['ahste']
        self.ahstd          = var['ahstd']
        return


    def get_dudls_and_averages_from_self_data(self):
        ''' 1st index: volume / lattice constant    (self.a)
            2nd index: temperature                  (self.t)
            3rd index: lambda                       (self.s)'''
        print("building dudls ...")
        self.dudl = np.empty((  self.a.size, \
                                self.t.size, \
                                self.l.size, \
                                self.s.size), dtype=object)
        self.dudl[:] = np.NAN

        self.dudlmean     = np.empty([self.a.size, self.t.size, self.l.size])
        self.dudllen      = np.empty([self.a.size, self.t.size, self.l.size])
        self.dudlste      = np.empty([self.a.size, self.t.size, self.l.size])
        self.dudlstd      = np.empty([self.a.size, self.t.size, self.l.size])
        self.dudlmeanlen  = np.empty([self.a.size, self.t.size])
        self.dudlmeanfit  = np.empty([self.a.size, self.t.size])

        self.dudlmean[:]     = np.NAN
        self.dudllen[:]      = np.NAN
        self.dudlste[:]      = np.NAN
        self.dudlstd[:]      = np.NAN
        self.dudlmeanlen[:]  = np.NAN
        self.dudlmeanfit[:]  = np.NAN

        #self.dudlavg = np.copy(self.dudl)

        np.set_printoptions(suppress=True)   # display arrays withou 000000
        np.set_printoptions(precision=6)    # print only 6 digist after .

        for a in range(self.a.size):
            for t in range(self.t.size):
                for l in range(self.l.size):
                    data = np.array([])
                    for s in range(self.s.size):
                        if np.isnan(np.any(f.data[a,t,l,s])):
                            # isnan -> no information
                            #print a,t,l,s,"NONE"
                            pass
                        else:
                            # This is the only input file specifit line
                            # sphinx
                            self.dudl[a,t,l,s] = self.data[a,t,l,s][:,7]

                            # remove first steps of equilibration
                            cut = self._dudl_cutoff_first_structures
                            self.dudl[a,t,l,s] = self.dudl[a,t,l,s][cut:]

                            # concat all dudls
                            data = np.concatenate((data,self.dudl[a,t,l,s]))
                            #print a,t,l,s,"--> DATA"

                            # get dudlavg
                            #self.dudlavg[a,t,l,s] = \
                            #[self.dudl[a,t,l,s][0:i].mean() for i in \
                            #np.arange(len(self.dudl[a,t,l,s]))+1] # sphinx

                    #print a,t,l,s,"||", self.a[a], self.t[t], self.l[l]
                    #print "data:",data
                    corl = self._dudl_correlation_length
                    if data is not np.array([]):
                        self.dudlmean[a,t,l]      = data.mean()
                        self.dudlste[a,t,l]   = \
                            data[0::corl].std()/np.sqrt(len(data[0::corl]))
                        self.dudlstd[a,t,l]   = data.std()
                        self.dudllen[a,t,l]   = len(data)
                #print "a:",a,"t:",t
                # a and t are just an index
                # building dudls ...
                # a: 0 t: 0
                # a: 0 t: 1
                # a: 0 t: 2
                # a: 0 t: 3
                # a: 0 t: 4
        self.dudlmeanfit = duelmean_to_dudlmeanfit(self.dudlmean, self.l)

                #lambdas_array = self.l
                #energies_array = self.dudlmean[a, t, :]
                ##self.dudlmeanfit[a,t] = self.get_dudlmeanfit(a, t)
                #print "ii:",lambdas_array,type(lambdas_array),energies_array,type(energies_array)
                #self.dudlmeanfit[a,t] = get_dudlmeanfit(lambdas_array, energies_array)
                ##print "sdmf:",self.dudlmeanfit[a,t]
        return


    def dudlmean_to_dudlmeanfit(dudlmean = None, lambdas = None):
        if dudlmean == None:
            sys.exit("expecting dudlmean to be an array[a,t,l] of energies")
        if type(dudlmean) != np.ndarray:
            sys.exit("expecting dudlmean to be an array[a,t,l] of energies")
        if lambdas == None:
            sys.exit("expecting lambdas to be an array; e.g. [0.0, 0.15, 0.5, 0.85, 1.0]")
        if len(lambdas) != dudlmean.shape[-1]:
            sys.exit("len(lambdas) != dudlmean.shpae[-1] which it has to be!")

        dudlmeanfit  = np.empty([dudlmean.shape[0], dudlmean.shape[1]])
        dudlmeanfit[:]  = np.NAN
        for inda in np.arange(dudlmean.shape[0]):
            for indt in np.arange(dudlmean.shape[1]):
                energies = dudlmean[inda, indt, :]
                dudlmeanfit[a,t] = get_dudlmeanfit(lambdas, energies)
        return dudlmeanfit

    def get_dudl_vs_temp_folder(self, verbose = False):
        ''' needs:
            - self.dudlmean
            - self.a
            - self.t
            - self.dudldste
            '''
        def quadFunction(x, a0):
            #return a0*x**2.718
            return a0*x**2.
            #return -1.+np.exp(a0*x)
            #return -1.+x*np.exp(a0)
            #return -1.+a0*np.exp(a0*x)

        folder = "auswertung_dudl_vs_temp"

        # remove folder if existing and create new one
        if os.path.isdir(folder):
           shutil.rmtree(folder)
        os.makedirs(folder)

        self.dudlmean_from_quadfit = np.copy(self.dudlmean)

        #for inda,a in np.arange(len(self.a)):
        for inda,a in enumerate(self.a):
            #for indl,l in np.arange(len(self.l)):
            for indl,l in enumerate(self.l):
                if verbose:
                    print("")
                    print("a:",a," l:",l)
                #print "a:",self.a[a],"l:",self.l[l]
                x = self.t
                y = self.dudlmean[inda,:,indl]
                dy = self.dudlste[inda,:,indl]
                dellist = np.array([])
                for iiind,ii in enumerate(y):
                    #print "Y:",iiind,y[iiind]
                    if np.isnan(y[iiind]) == True:
                        dellist = np.append(dellist,iiind)
                        #print "K:",iiind
                #print "dellst:",dellist,len(dellist)
                #print ">x:",x
                #print ">y:",y
                #print "?1:",y[~np.isnan(y)]   #~ reverts np.isnan
                #print "?2:",~np.isnan(y)   #~ reverts np.isnan
                #print "?x:",x
                #print "?y:",y
                y = y[~np.isnan(y)]   #~ reverts np.isnan
                # NO THIS IS HERE WRONG: !x = x[~np.isnan(y)]   #~ reverts np.isnan
                if len(dellist) > 0:
                    x = np.delete(x,dellist)
                    dy = np.delete(dy,dellist)
                    ## NO THIS IS HERE WRONG: y = np.delete(x,dellist)
                #print "!x:",x
                #print "!y:",y

                if len(x) != len(y) != len(dy):
                    sys.exit("len(x):"+str(len(x))+" is not len y:"+str(len(y)))
                if verbose:
                    print("x:",x)
                    print("y:",y)
                if len(x) != len(y):
                    print("len x != len y !!!!!!!!! ERROR")
                    continue
                if len(x) == 0:
                    continue
                #print "dy:",dy
                #print ""
                filebase = folder+"/"+str(a)+"_"+str(l)
                np.savetxt(filebase+"_data",np.transpose([x,y,dy]), fmt="%.0f %.3f %.3f")

                from scipy.optimize import curve_fit
                parameters, covariance = \
                    curve_fit(quadFunction, x, y, maxfev=33000)

                #print "pars:",parameters
                #print "pars:",parameters[0]
                # fit
                points = 101
                fit_x_min = 0
                fit_x_max = self.t.max()

                fit_x_values = np.linspace(fit_x_min, fit_x_max, points)
                fit_y_values = quadFunction(fit_x_values, *parameters)
                np.savetxt(filebase+"_fit",np.transpose([fit_x_values,fit_y_values]), fmt="%.0f %.6f")

                for indt,t in enumerate(self.t):
                    self.dudlmean_from_quadfit[inda,indt,indl] = quadFunction(t,*parameters)

        # make Fqh surface from quad fit
        self.dudlmean_from_quadfit_samenans = np.copy(self.dudlmean_from_quadfit)
        self.dudlmean_from_quadfit_samenansfit = np.copy(self.dudlmeanfit)
        self.dudlmean_from_quadfit_samenansfit[:] = np.NAN
        for a in np.arange(self.dudlmean.shape[0]):
            for t in np.arange(self.dudlmean.shape[1]):
                for e in np.arange(self.dudlmean.shape[2]):
                    print("a:",a," t:",t,self.dudlmean[a,t,e])
                    if np.isnan(self.dudlmean[a,t,e]) == True:
                        self.dudlmean_from_quadfit_samenans[a,t,e] = np.nan
                #print "???",self.l, self.dudlmean_from_quadfit_samenans[a,t]
                #print "???", get_dudlmeanfit(self.l, self.dudlmean_from_quadfit_samenans[a,t])
                self.dudlmean_from_quadfit_samenansfit[a,t] = get_dudlmeanfit(self.l, self.dudlmean_from_quadfit_samenans[a,t])
        for a in np.arange(self.dudlmeanfit.shape[0]):
            for t in np.arange(self.dudlmeanfit.shape[1]):
                if np.isnan(self.dudlmeanfit[a,t]) == True:
                    self.dudlmean_from_quadfit_samenansfit[a,t] = np.nan

        # make Fqh surface from quad fit only highest temperature
        self.dudlmean_from_quadfit_samenans_highestT = np.copy(self.dudlmean)
        self.dudlmean_from_quadfit_samenans_highestTfit = np.copy(self.dudlmeanfit)
        self.dudlmean_from_quadfit_samenans_highestTfit[:] = np.NAN
        for a in np.arange(self.dudlmeanfit.shape[0]):
            for t in np.arange(self.dudlmean.shape[1]):
                for e in np.arange(self.dudlmean.shape[2]):
                    if t != self.dudlmean.shape[1]-1:
                        self.dudlmean_from_quadfit_samenans_highestT[a,t,e] = np.nan
                print("??? a:",a,"t:",t,self.l, self.dudlmean_from_quadfit_samenans_highestT[a,t])
                print("??? a:",a,"t:",t, get_dudlmeanfit(self.l, self.dudlmean_from_quadfit_samenans_highestT[a,t]))
                print("!!!", np.isnan(self.dudlmean_from_quadfit_samenans_highestT[a,t]))
                if np.isnan(get_dudlmeanfit(self.l, self.dudlmean_from_quadfit_samenans_highestT[a,t])):
                    print("isnan!!!!!!!!!!!!!")
                self.dudlmean_from_quadfit_samenans_highestTfit[a,t] = get_dudlmeanfit(self.l, self.dudlmean_from_quadfit_samenans_highestT[a,t])

        return


    def plot_dudl_vs_lambda(self, a, t, l):
        ''' TODO: NOT READY YET
            1st index: a
            2nd index: t
            3rd index: l '''
        energies_array = None
        if utils.is_int(a) is True:
            a = int(a)
        elif a is ":":
            a = eval(":")
            #energies_array = self.dudlmean[:]
        else: sys.exit("index a:"+str(a)+" not correct")

        energies_array = self.dudlmean[a]

        #if utils.is_int(t) is True:
        #    energies_array = energies_array[int(t)]
        #elif t is ":":
        #    energies_array = energies_array[:]

        print("energies_array:",energies_array)
        return


    def import_fah_data(self, filename = None):
        ''' '''
        self.surface = None
        possible_filename = [ "Fah_surface" ]
        if filename == None:
            for i in possible_filename:
                filecheck = glob.glob(i)
                if filecheck:
                    if len(filecheck) == 1:
                        _printgreen("importing fah data : "+filecheck[0]+" "+"."*20)
                        self.surface = np.loadtxt(filecheck[0])
                        break

        if self.surface == None:
            sys.exit("fah surface not found!")
        _printgreen("importing fah data ............................... DONE")
        print("")
        return self.surface


    def import_surface(self, surface_filename = None, interpolate = False, sysexit = True):
        ''' '''
        surface = None
        if type(surface_filename) == bool or surface_filename == None:
            surface_filename = find_surface_filename(sysexit = sysexit)

        # check surface_filename
        #print "000:",surface_filename
        if type(surface_filename) == bool or surface_filename == None:
            if sysexit:
                sys.exit("Did not found any Fah surface file in current folder")
            else:
                return False

        if os.path.isfile(surface_filename) != True:
            if sysexit:
                sys.exit("Did not found any Fah surface file in current folder")
            else:
                return False

        self.surface_filename = surface_filename
        surface = np.loadtxt(surface_filename)
        if interpolate:
            surface = h.fqh_surface_interpolate(surface)
        return surface


    def import_surface_expr(self, surface_filename = None, surface = None, sysexit = True):
        if surface == None:
            surface = self.import_surface(surface_filename = surface_filename, sysexit = sysexit)
        try:
            if surface == False:
                return False
        except ValueError:
            pass # in this case we have an array and everything is correct
        self.surface_expr = np.array([])
        import sympy
        V = sympy.Symbol('V')
        for ind, temp in enumerate(np.arange(surface.shape[0])):
            add =    surface[ind][1] \
                    +surface[ind][2]*V \
                    +surface[ind][3]*V*V
            if surface.shape[1] == 4:
                self.surface_expr = np.append(self.surface_expr, add)
            elif surface.shape[1] == 5:
                self.surface_expr = np.append(self.surface_expr, add \
                    +surface[ind][4]*V*V*V)
        return self.surface_expr

if __name__ == '__main__':
    #f = ah()
    string = '''
    Examples for using this script:
    -------------------------------
    # creating jobs:
    # -----------------
    % fah.py -ef fah_create_jobs_from_fqh       (step 0)
    % fah.py -ef go_through_all_fah_jobs_and_create_joblist       (step 1)

    # analyzing jobs:
    # -----------------
    % fah.py -ef fah_go_through_all_angK_folder_and_exec_function get_Fah_and_avg_dudl_in_one_ang_K_folder_from_ipi_job   (step 2)
    % fah.py -ef fah_get_Fah_surface (step 3)
    % fah.py -ef fah_get_thermo (step 4)
    '''
    p = argparse.ArgumentParser(description=string,
            formatter_class=argparse.RawTextHelpFormatter)
    p.add_argument('-ws', default=False, action='store_true',
        help='Write Fqh_surface (Fah_surface and Fah_surface_x2) start this in high folder')
    p.add_argument('-ef','--execute_function', required=False, type=str,nargs='+',default='', help="function to run from this file.")
    p.add_argument('-lla', default=False, action='store_true',
        help='load low  avg_dUdL_fre files'),
    p.add_argument('-v', '--verbose',
                help='verbose', action='store_true', default=False)
    args = p.parse_args()

    if args.verbose:
        my.print_args(args)

    if args.execute_function:
        print('args.execute_function',args.execute_function)
        #sys.exit()
        if len(args.execute_function) == 1:
            function = eval(args.execute_function[0])
            print('function',function)
            function()
        elif len(args.execute_function) == 2:
            function = eval(args.execute_function[0])
            argument = eval(args.execute_function[1])
            print('function(argument)')
            print('function',function)
            print('argument',argument)
            function(argument)
        sys.exit()



    if args.ws:
        if args.verbose:
            print("###################### 1")
        f.import_avg_dudl_data(avg_dudl_filenames = "*Ang_*K/avg_dUdL_lowplushigh_fre", verbose = args.verbose)

        if args.verbose:
            print("###################### 2")
        f.get_dudl_vs_temp_folder(verbose = args.verbose)
        if args.verbose:
            print("###################### 3")
        write_Fah_surface(a = f.a, temps = f.t, dudlmeanfit = f.dudlmeanfit)
        write_Fah_surface(a = f.a, temps = f.t, dudlmeanfit = f.dudlmean_from_quadfit_samenansfit, filename = 'Fah_surface_x2')


    if args.lla:
        f.import_avg_dudl_data(avg_dudl_filenames = "*Ang_*K/avg_dUdL_fre", verbose = False)

    #read_from_files = True
    #if read_from_files:
    #    f = ah()
    #    f._filestring             = "4.16Ang_*K/lambda*/moldyn.dat*"
    #    f.initialize_files_data()
    #    f.import_dudl_files_data()
    #    f.save_data()
    #    f.get_dudls_and_averages_from_self_data()
    #    f.save_data()
    #else:
    #    pass
    #    #f.load_data()


    #f.get_dudls_and_averages_from_self_data()


    #a=1
    #plot(f.l,f.dudlmean[a,0],f.l,f.dudlmean[a,1],f.l,f.dudlmean[a,2],f.l,f.dudlmean[a,3],f.l,f.dudlmean[a,4])


    ## Fah vs volume
    #a=0
    #plot(f.a,f.dudlmeanfit[a])
    #

    ## Fah vs temperature
    #t=0
    #plot(f.t,f.dudlmeanfit[:,t]
