#!/usr/bin/env python

from __future__ import print_function
import numpy as np
import glob,sys,os,argparse
import myutils as my

def help(p = None):
    string = ''' helptext '''
    p = argparse.ArgumentParser(description=string,
            formatter_class=argparse.RawTextHelpFormatter)
    p.add_argument('-s','--from_subfolder', help='from_subfolder', action='count', default=True)
    p.add_argument('-q','--from_que'      , help='from_que'      , action='count', default=False)
    p.add_argument('-v','--verbose', help='verbose', action='count', default=False)
    return p

p = help()
args = p.parse_args()
if args.from_que == True: args.from_subfolder = False
verbose = args.verbose

print('args.from_subfolder',args.from_subfolder)
print('args.from_que      ',args.from_que)
print('args.verbose       ',args.verbose)

#ru=sorted(glob.glob(os.getcwd()+"/*/log.fit"))
#for i in ru:
#    print(i)
#sys.exit()
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
pp()
print("#             || rmse  /  rmse        || rmse / rmse         || [stes]")
print("#train        || test  /  train       || train/ test         || [total]")
print("#frac         || min   /      (@step) || min  /      (@step) ||         path")
pp()

##########################################################
#for c in ["*"]:  # jst use the subfolder
subfolder = ["*"]
if args.from_que == True:
    subfolder = path
for c in subfolder:    # from the ones in the que
    #fo=glob.glob("tf_*_"+c+"_cores*/learning-curve.out")
    fn=sorted(glob.glob(c+"/learning-curve.out"))
    ru=sorted(glob.glob(os.getcwd()+"/*/log.fit"))
    fo = fn+ru
    #for i in fo:
    #    print('fo',i,os.path.basename(i))

    ##if len(fo) == 0:
    #    fo=glob.glob("tf_*_"+c+"*learning-curve.out")
    out=[]
    out2=[]
    for i in fo:
        #print('test_que:',i)
        basename = os.path.basename(i)
        inputnn=i.replace(basename, 'input.nn')
        elastic=i.replace(basename, 'elastic.dat')
        elastic_ene=i.replace(basename, 'elastic_ene.dat')
        testf = my.inputnn_get_testfraction(inputnn)
        pot_elements, pot_atom_energy = my.inputnn_get_atomic_symbols_and_atom_energy(inputnn)
        #print('pe',pot_atom_energy,inputnn)
        try:
            mg = int(pot_atom_energy["Mg"])*-1
        except KeyError:
            mg = 0
        try:
            al = int(pot_atom_energy["Al"])*-1
        except KeyError:
            al = 0
        try:
            si = int(pot_atom_energy["Si"])*-1
        except KeyError:
            si = 0

        #print('pot_atom_energy Mg',mg,si,al,inputnn)

        if os.path.isfile(elastic):
            c44 = np.loadtxt(elastic)
        else:
            c44 = 0

        if os.path.isfile(elastic_ene):
            c44e = np.loadtxt(elastic_ene)
        else:
            c44e = 0
        train_fraction=1.-testf
        train_fraction=1.-testf

        learning_curve = lc = my.n2p2_runner_get_learning_curve(i)

        len_ = len(lc[:,1])
        trainmin                = lc[:,1].min()                      # best train RMSE
        trainmin_idx            = np.where(trainmin==lc[:,1])[0][0]  # best train index
        testrmse_at_trainmin    = lc[:,2][trainmin_idx]              # test RMSE @ train index

        testmin                 = lc[:,2].min()                      # best test RMSE
        testmin_idx             = np.where(testmin==lc[:,2])[0][0]   # best test index
        trainrmse_at_testmin    = lc[:,1][testmin_idx]               # train RMSE @ test index

        trainminf_at_testmin    = lc[:,3][testmin_idx]
        testminf_at_testmin     = lc[:,4][testmin_idx]

        path__ = i.replace(os.getcwd()+'/',"")

        trainmin = al
        testrmse_at_trainmin = mg
        trainmin_idx = si

        out2.append([
            round(train_fraction,2),                # j[0]

            #round(testmin*1000.*27.211384,2),       # j[1]
            round(testmin,2),       # j[1]
            #round(trainrmse_at_testmin*1000.*27.211384,2),     # j[2]
            round(trainrmse_at_testmin,2),     # j[2]
            testmin_idx,                               # j[3]

            #round(trainmin*1000.*27.211384,2),      # j[4]
            round(trainmin,2),      # j[4]
            #round(testrmse_at_trainmin*1000.*27.211384,2),    # j[5]
            round(testrmse_at_trainmin,2),    # j[5]
            trainmin_idx,                              # j[6]

            #round(trainminf_at_testmin*51.422063*1000,2),      # j[7]
            #round(testminf_at_testmin*51.422063*1000,2),    # j[8]
            round(trainminf_at_testmin,2),      # j[7]
            round(testminf_at_testmin,2),    # j[8]
            testmin_idx,                              # j[9]

            len_,                                   # epochs_
            c44,                                    # c44_
            c44e,                                    # c44_
            path__                                       # path_
            ])

    np.set_printoptions(precision=2)
    np.set_printoptions(suppress=True)
    out2=sorted(out2,key=lambda x: x[0])
    for idj,j in enumerate(out2): # for every line
        run = "    "
        NJC = "    "
        que = "    "
        epochs_=len(j) - 4
        c44e_=len(j) - 3
        c44_=len(j) - 2
        path_=len(j) - 1
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
        if j[epochs_] - 250 < j[3]:
            NJC = "NJC "
        if len(j[path_].split("vorlage_parallel_mode2")) == 2:
            continue   # no difference just makes fitting faster/slower
        if len(j[path_].split("vorlage_parallel_mode3")) == 2:
            continue   # no difference just makes fitting faster/slower
        if len(j[path_].split("ff0.02115_kalman1")) == 2:
            continue   # does not run correctly
        if len(j[path_].split("ff0.02115_vorlage_kalman1")) == 2:
            continue   # does not run correctly
        if len(j[path_].split("ff0.02115_normalize_nodes")) == 2:
            continue # too slow
        if len(j[path_].split("vorlage_force_weight_12")) == 2:
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
        print(run+NJC+"%0.2f  || %5.1f /%5.1f  (%4.0f) ||%5.1f /%5.1f (%4.0f)  || %5.1f /%5.1f (%4.0f) || c44 %3.1f %3.1f || [%4.0f] %s"%(j[0],j[1],j[2],j[3],j[4],j[5],j[6],   j[7],j[8],j[9],     j[c44_], j[c44e_] ,j[epochs_],j[path_]))



