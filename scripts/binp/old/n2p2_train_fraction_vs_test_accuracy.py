#!/usr/bin/env python

from __future__ import print_function
import numpy as np
import glob,sys,os,argparse,subprocess,time
import myutils as my
from myutils import ase_calculate_ene
#print('imported all ...')
start_time = time.time()

def help(p = None):
    string = ''' helptext '''
    p = argparse.ArgumentParser(description=string,
            formatter_class=argparse.RawTextHelpFormatter)
    p.add_argument('-sbtrain','--sort_by_trainfraction', action='count',default=False, help='sort output by trainfractino')
    p.add_argument('-s'  ,'--from_subfolder'     , action='count',      default=True,  help='from_subfolder')
    p.add_argument('-o'  ,'--only_best_converged', action='count',      default=True,  help='show only best converged job when several restarts')
    p.add_argument('-q'  ,'--from_que'           , action='count',      default=False, help='from_que')
    p.add_argument('-p'  ,'--potential'          , action='count',      default=False, help='potential is in "potential" subfolder')
    p.add_argument('-c44','--getc44'             , action='count',      default=False, help='get c44 in potential folder')
    p.add_argument('-v'  ,'--verbose'            , action='count',      default=False, help='verbose')
    p.add_argument('-b'  ,'--both'               , action='store_true', default=False, help='evaluate for runner & n2p2 (instad of only n2p2)')
    p.add_argument('-e'  ,'--execute'            , action='store_true', default=False, help='execute part')
    p.add_argument('-f'  ,'--find'               , default=-1.,type=str, nargs='*',required=False, help="find (& restrict selectiion to) pattern in foldername")
    return p

p = help()
args = p.parse_args()
if args.from_que == True: args.from_subfolder = False
verbose = args.verbose

if args.verbose:
    print('args.from_subfolder',args.from_subfolder)
    print('args.from_que      ',args.from_que)
    print('args.verbose       ',args.verbose)

if args.verbose:
    print('begin (1) ...')

fo=sorted(glob.glob(os.getcwd()+"/*/learning-curve.out"))

if verbose:
    print('----- currently considered pathes -----')
    for i in fo:
        print(i)
    print()

if verbose:
    print('----- currently in the que -----')
id,stat,path = my.q()
if verbose:
    for idx,i in enumerate(id):
        print(idx,id[idx],stat[idx],path[idx])
    print('----- que done ------')

checkpath_run = []
checkpath_que = []
subfolderbase = "/scratch/glensk/"
subfolder = "n2p2_24_24/"
subfolder_len = len(subfolderbase+subfolder)

if args.verbose:
    print('gc ...')
for idx,i in enumerate(path):
    for g in fo:
        gc = g.split("/learning-curve.out")[0]
        #print('gc',gc)
        if gc == i and stat[idx] == 'R':
        #if gc in i and stat[idx] == 'R':
        #if gc in i.split("/") and stat[idx] == 'R':
            #print('--> gc',gc)
            #print('--> i ',i)
            checkpath_run.append(gc)
        if gc in i.split("/") and stat[idx] == 'Q':
            checkpath_que.append(gc)

if verbose:
    print('----- given the status R -----')
    for i in checkpath_run:
        print(i)
    print()
    print('----- given the status Q -----')
    for i in checkpath_que:
        print(i)
#sys.exit()
def pp():
    a="-----------"
    print(12*a)
pp()
print(' - n2p2 has spikes for t_t_l')
pp()
pp()
print("#           ||    ENERGY (RMSE)    ||          || FORCES (RMSE)|| C44  || KMC  || [stes]")
print("#train      || test / train        ||          || train/ test  || C44  || std  || [total]")
print("#frac       || min  /       (@step)||          || min  /       || C44  ||      ||   path")
pp()


##########################################################
#for c in ["*"]:  # jst use the subfolder
subfolder = ["*"]
if args.from_que == True:
    subfolder = path
allfolder = []
drawn1 = False
drawn2 = False
if args.verbose:
    print('c ...')


def foldername_search_restartname(folder):
    iteration = folder.split("_")[-1] # 1 2 tmp
    preiteration = "_".join(folder.split("_")[:-1]) # 1 2 tmp
    #print('i       ',i)              # runner_4998_21_3/log.fit
    #print('folder  ',folder)         # runner_4998_21_3
    #print('iteration',iteration)      # 1 2 tmp
    #print('preiteration',preiteration)      # 1 2 tmp
    try:
        search = preiteration+"_"+str(int(iteration)+1)
        #print('search',search)
    except ValueError:
        return False
    return search

def get_try_get_kmcstd_from_kmc_epoch_xx(kmc,folder,index=False):
    file = folder+"/kmc/ene_std_epoch_"+str(index)+".dat"
    if not os.path.isfile(file):
        return kmc
    else:
        kmc = np.loadtxt(file)
        kmc = kmc[-1]
        return kmc
    return

def get_try_get_c44_from_elastic_all(c44,elastic_all,index=False):
    if not os.path.isfile(elastic_all):
        #print('is not')
        return c44
    else:
        #print('y1')
        c44_try = np.loadtxt(elastic_all)
        #print(c44_try)
        takeline = np.where(c44_try[:,0]==index)[0][0]
        try:
            c44out = c44_try[takeline][1]
        except IndexError:
            c44out = c44
        #print('y2',takeline,c44out)
        return c44out

    return c44

###################################################################
# from here on search all the learning curve files
###################################################################
for c in subfolder:    # from the ones in the que
    fm1 = ""
    if args.potential:
        fm1 = "/potential"
    if verbose > 3:
        print('lookfor',c+fm1+"/learning-curve.out")
    ##########################
    ## search for n2p2 jobs
    ##########################
    all_learning_curve_files = fn =sorted(glob.glob(c+fm1+"/learning-curve.out"))
    if verbose > 3:
        print('fn n2p2 jobs',fn)
    ##########################
    ## search for n2p2 jobs
    ##########################
    if args.both:
        ru=sorted(glob.glob(c+fm1+"/log.fit"))
        all_learning_curve_files = fn+ru

    out2=[]
    if args.verbose:
        print("for every fo")
    ##########################################
    # go over all all_learning_curve_files
    ##########################################
    for i in all_learning_curve_files:
        if verbose > 0:
            print('learning_curve_file == i:',i)
        # i        : runner_v0_64/log.fit
        # basename : log.fit
        # folder   : runner_v0_64
        basename = os.path.basename(i)
        folder = i.replace(basename, '')[:-1]
        allfolder.append(folder)
        #if verbose:
        #################################################################
        # in case there is something to do in n2p2/runner folder
        #################################################################
        if args.execute == True:
            print('learning_curve_file folder:',folder.ljust(40),os.getcwd())
            hier = os.getcwd()
            os.chdir(folder)
            subprocess.call(["n2p2_get_potential_folder_from_nr.py"],shell=True)
            #subprocess.call(["getEnergies_byLammps.py -p . -ctrain"],shell=True)
            #subprocess.call(["getEnergies_byLammps.py -p . -ctest"],shell=True)
            #subprocess.call(["getEnergies_byLammps.py -p . -ckmc"],shell=True)
            #subprocess.call(["getEnergies_byLammps.py -p . -cinput"],shell=True)
            os.chdir(hier)

        if args.getc44:
            with my.cd(folder):
                if verbose:
                    print("--------->>",os.getcwd())
                subprocess.call("getEnergies_byLammps.py -p . -e",shell=True)
        #print(os.getcwd())
        #sys.exit()

        foldern = foldername_search_restartname(folder)
        if verbose > 3:
            print('folder  ',folder)         # runner_4998_21_3
            print('foldern ',foldern)
        inputnn     =i.replace(basename, 'input.nn')
        inputdata   =i.replace(basename, 'input.data')

        input_structures = my.inputdata_get_nuber_of_structures(inputnn)
        kmcstdfile  =i.replace(basename, 'kmc57/ene_std.npy')
        elastic     =i.replace(basename, 'elastic.dat')
        elastic_ene =i.replace(basename, 'elastic_ene.dat')
        elastic_all =i.replace(basename, 'elastic_c44_all.dat')
        runner_n2p2 = my.inputnn_runner_or_n2p2(inputnn)
        if verbose > 3:
            print('inputnn          :',inputnn)
            print('inputdata        :',inputdata)
            print('input_structures :',input_structures)
            print('runner_n2p2      :',runner_n2p2)
        #pot_epoch = my.inputnn_get_potential_number_from_weightsfile(inputnn)
        if True:
            ace = ase_calculate_ene(pot=False,
                    potpath=folder,
                    use_different_epoch=False,
                    units="meV_pa",
                    verbose=args.verbose)

        #print(inputnn,'runner_n2p2',runner_n2p2)
        testf       = my.inputnn_get_testfraction(inputnn)
        random_seed = my.inputnn_get_random_seed(inputnn)
        nodes_short = my.inputnn_get_nodes_short(inputnn,as_string=True)
        activation_short = my.inputnn_get_activation_short(inputnn)
        nn = nodes_short+"__"+activation_short
        if os.path.isdir(folder+"/kmc"):
            nn = nodes_short+"**"+activation_short

        #print('nodes_short',nodes_short)
        #print('activation_short',activation_short)
        #print("nn",nn)
        #pot_elements, pot_atom_energy = my.inputnn_get_atomic_symbols_and_atom_energy_dict(inputnn)
        pot_elements, [al,mg,si]= my.inputnn_get_atomic_symbols_and_atom_energy_list(inputnn)
        al = int(al)
        mg = int(mg)
        si = int(si)

        #print('pot_atom_energy Mg',mg,si,al,inputnn)
        if os.path.isfile(kmcstdfile):
            kmcstd_all = np.loadtxt(kmcstdfile)
            kmcstd = kmcstd_all[-1]
        else:
            kmcstd = '-'

        if os.path.isfile(elastic):
            c44 = np.loadtxt(elastic)
        else:
            c44 = 0

        if os.path.isfile(elastic_ene):
            c44e = np.loadtxt(elastic_ene)
        else:
            c44e = 0

        #print('c44',int(c44),int(c44e),folder)
        train_fraction =1.-testf

        learning_curve = lc = my.n2p2_runner_get_learning_curve(inputnn)
        #print(lc)
        #print()
        #print(lc[:,[0,2]])
        #print('folder',folder)
        if False:
            np.savetxt(folder+"/learning-curve.e.train.out",lc[:,[0,1]])
            np.savetxt(folder+"/learning-curve.e.test.out",lc[:,[0,2]])
            np.savetxt(folder+"/learning-curve.f.train.out",lc[:,[0,3]])
            np.savetxt(folder+"/learning-curve.f.test.out",lc[:,[0,4]])
        #sys.exit()
        if False:
            print('i',i)
            print('len',len(lc),lc.shape)
            print('lc')
            print(lc)
            sys.exit()

        #if len(lc) == 1:
        #    print('---')
        #    print(lc)
        #    print('---')
        #print('lc',lc.shape)
        epochs_ = len(lc[:,1])
        #print('len',len(lc),epochs_)
        #print()
        def get_try_get_kmcstd_from_kmc_epoch_xx(kmcstd,folder,index=False):
            file = folder+"/kmc/ene_std_epoch_"+str(index)+".dat"
            if not os.path.isfile(file):
                return kmcstd
            else:
                #print('loading',file)
                kmcstd_all = np.loadtxt(file)
                kmcstd = kmcstd_all[-1]
                return kmcstd
            return

        def get_try_get_c44_from_elastic_all(c44,elastic_all,index=False):
            if not os.path.isfile(elastic_all):
                #print('is not')
                return c44
            else:
                #print('y1')
                c44_try = np.loadtxt(elastic_all)
                #print(c44_try)
                takeline = np.where(c44_try[:,0]==index)[0][0]
                try:
                    c44out = c44_try[takeline][1]
                except IndexError:
                    c44out = c44
                #print('y2',takeline,c44out)
                return c44out

            return c44

        for bl in [ 'best','last']:
            if bl == 'best':
                test             = lc[:,2].min()                      # best test RMSE
                testmin_idx      = np.where(test==lc[:,2])[0][0]   # best test index
            elif bl == 'last':
                testmin_idx = len(lc)-1
                #print('lc',lc)
                #print(inputnn)
                #sys.exit()
            #sys.exit()
            # #print()
            #print('c44 in',c44)
            c44 = get_try_get_c44_from_elastic_all(c44,elastic_all,index=testmin_idx)
            kmcstd = get_try_get_kmcstd_from_kmc_epoch_xx(kmcstd,folder,index=testmin_idx)
            #print('c44 ou',c44)
            #print()

            test     = lc[:,2][testmin_idx]                      # best test RMSE
            train    = lc[:,1][testmin_idx]               # train RMSE @ test index

            trainminf_at_testmin    = lc[:,3][testmin_idx]
            testminf_at_testmin     = lc[:,4][testmin_idx]

            path__ = i.replace(os.getcwd()+'/',"")

            out2.append([
                round(train_fraction,2),                            # j[0]

                round(test,2),                                   # j[1]
                round(train,2),                      # j[2]
                testmin_idx,                                        # j[3]

                round(al,2),                                       # j[4]
                round(mg,2),                                       # j[5]
                round(si,2),                                       # j[6]

                round(trainminf_at_testmin,2),                      # j[7]
                round(testminf_at_testmin,2),                       # j[8]
                testmin_idx,                                        # j[9]

                bl,
                input_structures,
                random_seed,                                        # j[10]
                runner_n2p2,                                        # j[1]
                foldern,                                            # j[1] path_
                nn,                                                 # j[1]
                kmcstd,                                             # j[1]
                epochs_-1,                                            # j[1] epochs_
                c44,                                                # j[1] c44_
                c44e,                                               # j[1] c44_
                folder                                              # j[1] path_
                ])

    np.set_printoptions(precision=2)
    np.set_printoptions(suppress=True)
    if args.sort_by_trainfraction:
        out2=sorted(out2,key=lambda x: x[0])

    #out2=sorted(out2,key=lambda x: x[2])               # das ist nach dem train ergebnis
    #out2=sorted(out2,key=lambda x: x[len(out2[0])-3])  # das ist nach dem train ergebnis
    if args.verbose:
        print('++++++++++++++++++++++++++++++')

    if args.both: n2p2_or_runner_loop_ = ["n2p2","runner"]
    else: n2p2_or_runner_loop_ = ["n2p2"]
    for n2p2_or_runner_loop in n2p2_or_runner_loop_:
        if args.verbose > 3:
            print('AA',n2p2_or_runner_loop)
        if True:
        #for exec_conv_unconv_b in ["orange","green","white",'red','blue']:
        #    exec_conv_unconv = [n2p2_or_runner_loop,exec_conv_unconv_b]
            for idj,j in enumerate(out2): # for every line
                if args.verbose > 3:
                    print('CC',j)
                run = "    "
                NJC = "    "
                que = "    "
                bl = best_last = j[len(j) - 11]   # best & last should be printd right one after another
                ist = input_structures = j[len(j) - 10]
                random_seed = rnd = j[len(j) - 9]
                runner_n2p2 = j[len(j) - 8]
                foldern     = j[len(j) - 7]
                if foldern in allfolder:
                    if args.verbose > 3:
                        print('DD is in --> continue',foldern)
                    continue
                nn              = j[len(j) - 6]
                kmc             = j[len(j) - 5]
                epochs          = j[len(j) - 4]
                c44             = j[len(j) - 3]
                c44e_unused     = j[len(j) - 2]
                path = folder   =  j[len(j) - 1]
                print('pathx',path)
                sys.exit()

                if True:
                    pot = my.mypot(False,folder,use_different_epoch=False,verbose=args.verbose)
                    pot.get(exit=False)
                    pot.get_my_assessments(get_outliers=True)
                    #epoch_best = pot.potepoch_bestteste
                    #print(folder,pot.potepoch_all,epoch_best)
                    if len(pot.potepoch_all) > 0:
                        epoch_last = pot.potepoch_all[-1]
                    else:
                        epoch_last = False
                    ab_test  = folder+"/assess_test_"+str(pot.potepoch_bestteste)
                    al_test  = folder+"/assess_test_"+str(epoch_last)
                    ab_train = folder+"/assess_train_"+str(pot.potepoch_bestteste)
                    al_train = folder+"/assess_train_"+str(epoch_last)
                    # first, check if still running;
                    if os.path.isfile(pot.potpath+'/elastic_c44_all.dat'):
                        elastic_c44_all = np.loadtxt(pot.potpath+'/elastic_c44_all.dat')
                        print(len(elastic_c44_all),'potpath',pot.potpath)
                        if len(elastic_c44_all) > 300:
                            pass

                    agree="?"
                    if bl == 'best':
                        kmcbl = pot.kmc57_b
                        if os.path.isfile(ab_test+'/ene_std.npy') and os.path.isfile(ab_train+'/ene_std.npy'):
                            agree = "-"
                            ab_test_ = float(my.read_lastline(ab_test+'/ene_std.npy'))
                            ab_train_ = float(my.read_lastline(ab_train+'/ene_std.npy'))
                            #agree=[ab_test_,ab_train_]
                            if np.abs(ab_test_  - j[1]) < 0.1 and np.abs(ab_train_ - j[2]) < 0.1: agree = "|"

                    if bl == 'last':
                        kmcbl = pot.kmc57_l
                        #if folder == "n2p2_v3ag_4998_new":
                        #    print(al_test+'/ene_std.npy')
                        #    print(al_train+'/ene_std.npy')
                        if os.path.isfile(al_test+'/ene_std.npy') and os.path.isfile(al_train+'/ene_std.npy'):
                            agree = "-"
                            al_test_ = float(my.read_lastline(al_test+'/ene_std.npy'))
                            al_train_ = float(my.read_lastline(al_train+'/ene_std.npy'))
                            #agree=[al_test_,al_train_]
                            if np.abs(al_test_  - j[1]) < 0.1 and np.abs(al_train_ - j[2]) < 0.1: agree = "|"

                            #print('diff',np.abs(al_test_  - j[1]))
                            #print('diff',np.abs(al_train_ - j[2]))
                    if False:
                        if bl == 'best':
                            print('best                         ',pot.potepoch_bestteste,agree,j[1],j[2],folder)
                        if bl == 'last':
                            print('last                         ',epoch_last,agree,j[1],j[2],folder)
                        if folder == "n2p2_v3ag_5000":sys.exit('33')


                path_before_learningcurve = path.split("/learning-curve.out")[0]  # random_seed_2234125
                #print('1path',path)
                path_before_learningcurve = os.getcwd()+"/"+path_before_learningcurve
                #print('2path',path_before_learningcurve)
                if path_before_learningcurve in checkpath_run:
                    run = "(R) "
                if path_before_learningcurve in checkpath_que:
                    run = "(Q) "
                if epochs - 250 < j[3]:
                    NJC = "NJC "
                if len(path.split("vorlage_parallel_mode2")) == 2:
                    if args.verbose > 3:
                        print('exit reason 1')
                    continue   # no difference just makes fitting faster/slower
                if len(path.split("vorlage_parallel_mode3")) == 2:
                    if args.verbose > 3:
                        print('exit reason 2')
                    continue   # no difference just makes fitting faster/slower
                if len(path.split("ff0.02115_kalman1")) == 2:
                    if args.verbose > 3:
                        print('exit reason 3')
                    continue   # does not run correctly
                if len(path.split("ff0.02115_vorlage_kalman1")) == 2:
                    if args.verbose > 3:
                        print('exit reason 4')
                    continue   # does not run correctly
                if len(path.split("ff0.02115_normalize_nodes")) == 2:
                    if args.verbose > 3:
                        print('exit reason 4')
                    continue # too slow
                if len(path.split("vorlage_force_weight_12")) == 2:
                    if args.verbose > 3:
                        print('exit reason 5')
                    continue

                #if len(path.split("tf_0.2_job_21maa_fr_ff0.01215")) == 2:
                #    print('j[epochs_]',j[epochs_],'j[3]',j[3])
                #    print('test1',j[epochs_] - 250)
                #print('j j j',j[1])
                if j[1] > 999: j[1] = 999.9
                if j[2] > 999: j[2] = 999.9
                if j[4] > 999: j[4] = 999.9
                if j[5] > 999: j[5] = 999.9
                if j[7] > 999: j[7] = 999.9
                if j[8] > 999: j[8] = 999.9
                if kmc > 99: kmc = 99

                #print('c44a',c44)
                if type(c44) != str:
                    #c44=str(int(np.round(c44,0))).ljust(3)
                    c44=str(round(c44,1)).ljust(4)
                #print('c44b',c44)

                if args.verbose > 3:
                    print('DD')

                e1 = str(int(j[4])).ljust(2)
                e2 = str(int(j[5])).ljust(2)
                e3 = str(int(j[6])).ljust(2)
                if bl == 'last':
                    e1 = ''.ljust(2)
                    e2 = ''.ljust(2)
                    e3 = ''.ljust(2)
                    path = ''.ljust(2)
                #kmc = str(round(kmc,1)).ljust(4)
                kmcbl = str(round(kmcbl,1)).ljust(4)

                stringout = run+NJC+"%0.1f |"+agree+"%5.1f /%5.1f  (%4.0f) || %s %s %s || %5.1f /%5.1f |C %s |K %s |n %4.0f || [%4.0f] | %s | %8.0f | %s"
                elementout = (        j[0]          ,  j[1],  j[2],   j[3],    e1,e2,e3,   j[7],  j[8],  c44, kmcbl,   ist ,    epochs,   nn,  rnd,   path)

                conv_unconv = "unconv"
                ## erstmak die komischen aussortieren
                if (j[1]+j[2])/2. < 0.1 or j[4] == 0.0:
                    takecolor = "blue"  # wiered
                elif c44 < 10:
                    takecolor = "blue"  # wiered
                elif j[4] < 3:  # old
                    takecolor = 'orange'
                elif (j[1]+j[2])/2. > 4. or (j[7]+j[8])/2. >= 35.:
                    takecolor = "red"
                elif (j[1]+j[2])/2. < 4.0 and (j[7]+j[8])/2. < 35.:
                    takecolor = "green"
                    conv_unconv = "conv"
                else:
                    takecolor = "white"

                if drawn1 == False and n2p2_or_runner_loop == 'runner':
                    print('------------------------------------------------------------------'*2)
                    drawn1 = True
                if args.verbose > 3:
                    print('EE',n2p2_or_runner_loop,runner_n2p2)
                    #print('FF',exec_conv_unconv[1],takecolor)
                if n2p2_or_runner_loop == runner_n2p2: # and exec_conv_unconv[1] == takecolor:
                    if takecolor == "orange":
                        print(my.printorange(stringout)%elementout)
                    if takecolor == "green":
                        print(my.printgreen(stringout)%elementout)
                    if takecolor == "blue":
                        print(my.printblue(stringout)%elementout)
                    if takecolor == "red":
                        print(my.printred(stringout)%elementout)
                    if takecolor == "white":
                        print(stringout%elementout)
#print('aall',allfolder)
end_time = time.time()
print("TIME:",end_time - start_time)
