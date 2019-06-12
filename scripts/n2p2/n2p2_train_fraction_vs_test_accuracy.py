#!/usr/bin/env python

from __future__ import print_function
import numpy as np
import glob,sys,os,argparse,subprocess
import myutils as my
#print('imported all ...')

def help(p = None):
    string = ''' helptext '''
    p = argparse.ArgumentParser(description=string,
            formatter_class=argparse.RawTextHelpFormatter)
    p.add_argument('-sbtrain','--sort_by_trainfraction'      , help='sort output by trainfractino'      , action='count', default=False)
    p.add_argument('-s','--from_subfolder', help='from_subfolder', action='count', default=True)
    p.add_argument('-o','--only_best_converged', help='show only best converged job when several restarts', action='count', default=True)
    p.add_argument('-q','--from_que'      , help='from_que'      , action='count', default=False)
    p.add_argument('-p','--potential'     , help='potential is in "potential" subfolder'      , action='count', default=False)
    p.add_argument('-c44','--getc44'     , help='get c44 in potential folder'      , action='count', default=False)
    p.add_argument('-v','--verbose', help='verbose', action='count', default=False)
    return p

p = help()
args = p.parse_args()
if args.from_que == True: args.from_subfolder = False
verbose = args.verbose

if args.verbose:
    print('args.from_subfolder',args.from_subfolder)
    print('args.from_que      ',args.from_que)
    print('args.verbose       ',args.verbose)

#ru=sorted(glob.glob(os.getcwd()+"/*/log.fit"))
#for i in ru:
#    print(i)
#sys.exit()
if args.verbose:
    print('fo ...')
fo=sorted(glob.glob(os.getcwd()+"/*/learning-curve.out"))
if verbose:
    print('----- currently considered pathes -----')
    for i in fo:
        print(i)
    print()
id,stat,path = my.q()
if verbose:
    print('----- currently in the que -----')
    for idx,i in enumerate(id):
        print(idx,id[idx],stat[idx],path[idx])
    print()

checkpath_run = []
checkpath_que = []
subfolderbase = "/scratch/glensk/"
subfolder = "n2p2_24_24/"
subfolder_len = len(subfolderbase+subfolder)

if args.verbose:
    print('gc ...')
for idx,i in enumerate(path):
    #print('i',i)
    #print(i.split("/scratch/glensk/n2p2_jobs"))
    #print('j',i.split("/"))
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

for c in subfolder:    # from the ones in the que
    fm1 = ""
    if args.potential:
        fm1 = "/potential"
    if verbose > 3:
        print('lookfor',c+fm1+"/learning-curve.out")
    fn=sorted(glob.glob(c+fm1+"/learning-curve.out"))
    if verbose > 3:
        print('fn',fn)
    #ru=sorted(glob.glob(os.getcwd()+"/*/log.fit"))
    ru=sorted(glob.glob(c+fm1+"/log.fit"))
    all_learning_curve_files = fn+ru
    out2=[]
    out_runner_conv = []
    out_runner_unconv = []
    out_n2p2_conv = []
    out_n2p2_unconv = []
    out_conv = []
    out_unconv = []
    if args.verbose:
        print("for every fo")
    for i in all_learning_curve_files:
        if verbose > 0:
            print('learning_curve_file == i:',i)
        # i        : runner_v0_64/log.fit
        # basename : log.fit
        # folder   : runner_v0_64
        basename = os.path.basename(i)
        folder = i.replace(basename, '')[:-1]
        allfolder.append(folder)
        if verbose:
            print('learning_curve_file folder',folder)
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
        input_structures = my.inputdata_get_nuber_of_structures(inputdata)
        kmcstdfile  =i.replace(basename, 'kmc57/ene_std.npy')
        elastic     =i.replace(basename, 'elastic.dat')
        elastic_ene =i.replace(basename, 'elastic_ene.dat')
        runner_n2p2 = my.inputnn_runner_or_n2p2(inputnn)
        if verbose > 3:
            print('inputnn          :',inputnn)
            print('inputdata        :',inputdata)
            print('input_structures :',input_structures)
            print('runner_n2p2      :',runner_n2p2)
        #pot_epoch = my.inputnn_get_potential_number_from_weightsfile(inputnn)

        #print(inputnn,'runner_n2p2',runner_n2p2)
        testf       = my.inputnn_get_testfraction(inputnn)
        random_seed = my.inputnn_get_random_seed(inputnn)
        nodes_short = my.inputnn_get_nodes_short(inputnn,as_string=True)
        activation_short = my.inputnn_get_activation_short(inputnn)
        nn = nodes_short+"__"+activation_short

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
        train_fraction=1.-testf
        train_fraction=1.-testf

        learning_curve = lc = my.n2p2_runner_get_learning_curve(inputnn)
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
        epochs_ = len(lc[:,1])
        #print('len',len(lc),epochs_)
        #print()

        for bl in [ 'best','last']:
            if bl == 'best':
                test             = lc[:,2].min()                      # best test RMSE
                testmin_idx      = np.where(test==lc[:,2])[0][0]   # best test index
                #if pot_epoch != False:
                #    testmin_idx         = pot_epoch
            elif bl == 'last':
                testmin_idx = len(lc)-1
                #print('lc',lc)
                #print(inputnn)
                #sys.exit()
            #sys.exit()

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
    #print('oo',len(out2))
    #print('oo',len(out2[0]))
    #print('oo',len(out2[0])-3)
    #sys.exit()
    #for i in out2:
    #    print(i[15])
    if args.sort_by_trainfraction:
        out2=sorted(out2,key=lambda x: x[0])

    #out2=sorted(out2,key=lambda x: x[2])               # das ist nach dem train ergebnis
    #out2=sorted(out2,key=lambda x: x[len(out2[0])-3])  # das ist nach dem train ergebnis

    #for exec_conv_unconv in [ ["n2p2","conv"],[ "n2p2", "unconv"], ['runner','conv'],['runner','unconv']]:
    #for exec_conv_unconv in [ ["n2p2","green"],[ "n2p2", "white"], ['n2p2','red'],['n2p2','blue'],["runner","green"],[ "runner", "white"], ['runner','red'],['runner','blue']]:
    for exec_conv_unconv_a in ["n2p2","runner"]:
        if args.verbose > 3:
            print('AA',exec_conv_unconv_a)
        for exec_conv_unconv_b in ["orange","green","white",'red','blue']:
            if args.verbose > 3:
                print('BB',exec_conv_unconv_b)
            exec_conv_unconv = [exec_conv_unconv_a,exec_conv_unconv_b]
            for idj,j in enumerate(out2): # for every line
                if args.verbose > 3:
                    print('CC',j)
                run = "    "
                NJC = "    "
                que = "    "
                bl = best_last = j[len(j) - 11]   # best & last should be printd right one after another
                ist = input_structures = j[len(j) - 10]
                random_seed =len(j) - 9
                runner_n2p2 =len(j) - 8
                foldern =len(j) - 7
                if j[foldern] in allfolder:
                    if args.verbose > 3:
                        print('DD is in --> continue',j[foldern])
                    continue
                nn     =len(j) - 6
                kmc    = j[len(j) - 5]
                epochs =j[len(j) - 4]
                c44    =j[len(j) - 3]
                c44e_  =len(j) - 2
                path_  = folder = len(j) - 1
                #print('-->',j[nn],j[foldern])
                #print('a',path_,j,'--->',j[11])
                #sys.exit()
                #print('kk',j[path_].split("/learning-curve.out"))
                #print('0j[path_]',j[path_])
                path_before_learningcurve = j[path_].split("/learning-curve.out")[0]  # random_seed_2234125
                #print('1j[path_]',j[path_])
                path_before_learningcurve = os.getcwd()+"/"+path_before_learningcurve
                #print('2j[path_]',path_before_learningcurve)
                if path_before_learningcurve in checkpath_run:
                    run = "(R) "
                if path_before_learningcurve in checkpath_que:
                    run = "(Q) "
                if epochs - 250 < j[3]:
                    NJC = "NJC "
                if len(j[path_].split("vorlage_parallel_mode2")) == 2:
                    if args.verbose > 3:
                        print('exit reason 1')
                    continue   # no difference just makes fitting faster/slower
                if len(j[path_].split("vorlage_parallel_mode3")) == 2:
                    if args.verbose > 3:
                        print('exit reason 2')
                    continue   # no difference just makes fitting faster/slower
                if len(j[path_].split("ff0.02115_kalman1")) == 2:
                    if args.verbose > 3:
                        print('exit reason 3')
                    continue   # does not run correctly
                if len(j[path_].split("ff0.02115_vorlage_kalman1")) == 2:
                    if args.verbose > 3:
                        print('exit reason 4')
                    continue   # does not run correctly
                if len(j[path_].split("ff0.02115_normalize_nodes")) == 2:
                    if args.verbose > 3:
                        print('exit reason 4')
                    continue # too slow
                if len(j[path_].split("vorlage_force_weight_12")) == 2:
                    if args.verbose > 3:
                        print('exit reason 5')
                    continue

                #if len(j[path_].split("tf_0.2_job_21maa_fr_ff0.01215")) == 2:
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
                    c44=str(int(np.round(c44,0))).ljust(3)
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

                stringout = run+NJC+"%0.1f ||%5.1f /%5.1f  (%4.0f) || %s %s %s || %5.1f /%5.1f |C %s |M %2.0f |S %4.0f  || [%4.0f] | %s | %8.0f | %s"
                elementout = (        j[0] ,  j[1],  j[2],   j[3],    e1,e2,e3,   j[7],  j[8], c44,     kmc,     ist ,     epochs, j[nn],j[random_seed],j[path_])

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

                if drawn1 == False and exec_conv_unconv[0] == 'runner':
                    print('------------------------------------------------------------------'*2)
                    drawn1 = True
                #if exec_conv_unconv[0] == j[runner_n2p2] and exec_conv_unconv[1] == conv_unconv:
                if args.verbose > 3:
                    print('EE',exec_conv_unconv[0],j[runner_n2p2])
                    print('FF',exec_conv_unconv[1],takecolor)
                if exec_conv_unconv[0] == j[runner_n2p2] and exec_conv_unconv[1] == takecolor:
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
