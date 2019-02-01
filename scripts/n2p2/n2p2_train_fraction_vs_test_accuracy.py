#!/usr/bin/env python

from __future__ import print_function
import numpy as np
import glob,sys,os
from myutils import q

#id,stat,path = q()
path=""
checkpath_run = []
checkpath_que = []
subfolderbase = "/scratch/glensk/"
subfolder = "n2p2_24_24/"
subfolder_len = len(subfolderbase+subfolder)
#for idx,i in enumerate(path):
#    #print('i',i)
#    #print(i.split("/scratch/glensk/n2p2_jobs"))
#    if i[:subfolder_len] != subfolderbase+subfolder:
#        continue
#    if type(i.split(subfolderbase+subfolder)) == list:
#        print('i >>>',i)
#        print("i -->",i.split(subfolderbase+subfolder))
#        #print('i >>>',i[:3])
#        if i[:3] != "/sc":
#            continue
#        if stat[idx] == 'R':
#            print(i[:3])
#            checkpath_run.append(i.split(subfolderbase+subfolder)[1])
#        else:
#            checkpath_que.append(i.split(subfolderbase+subfolder)[1])
#

def pp():
    a="-----------"
    print(12*a)


for c in ["*"]:
    #fo=glob.glob("tf_*_"+c+"_cores*/learning-curve.out")
    fo=sorted(glob.glob(c+"/learning-curve.out"))
    print('fo::',fo)
    #if len(fo) == 0:
    #    fo=glob.glob("tf_*_"+c+"*learning-curve.out")
    out=[]
    out2=[]
    for i in fo:
        #print('fo',i)
        #testf=float(i.split("_")[1])
        testf=0.2
        train_fraction=1.-testf
        lc = np.loadtxt(i) #+'/learning-curve.out')
        len_ = len(lc[:,1])
        #print(lc[:,1])
        #print(len(lc[:,1]))
        #sys.exit()

        trainmin                = lc[:,1].min()
        trainmin_idx            = np.where(trainmin==lc[:,1])[0][0]
        testrmse_at_trainmin    = lc[:,2][trainmin_idx]

        testmin       = lc[:,2].min()
        testmin_idx   = np.where(testmin==lc[:,2])[0][0]
        trainrmse_at_testmin  = lc[:,1][testmin_idx]

        #trainminf   = lc[:,3].min()
        #trainminf_  = np.where(trainminf==lc[:,3])[0][0]
        trainminf_at_testmin = lc[:,3][testmin_idx]

        #testminf    = lc[:,4].min()
        #testminf_   = np.where(testminf==lc[:,4])[0][0]
        testminf_at_testmin  = lc[:,4][testmin_idx]

        #print('train_fraction',train_fraction,'trainmin',np.where(trainmin==lc[:,1])[0][0],trainmin*1000.,'testmin',np.where(testmin==lc[:,2])[0][0],testmin*1000.)
        #out.append([
        #    train_fraction,                         # j[0]
        #    testmin*1000.*27.211384,                # j[1]
        #    trainmin*1000.*27.211384,               # j[2]
        #    i,                                      # j[3] path
        #    testmin_idx,                               # j[4]
        #    trainmin_idx,                              # j[5]
        #    len_                                    # j[6]
        #    ])
        out2.append([
            round(train_fraction,2),                # j[0]

            round(testmin*1000.*27.211384,2),       # j[1]
            round(trainrmse_at_testmin*1000.*27.211384,2),     # j[2]
            testmin_idx,                               # j[3]

            round(trainmin*1000.*27.211384,2),      # j[4]
            round(testrmse_at_trainmin*1000.*27.211384,2),    # j[5]
            trainmin_idx,                              # j[6]

            round(trainminf_at_testmin*51.422063*1000,2),      # j[4]
            round(testminf_at_testmin*51.422063*1000,2),    # j[5]
            testmin_idx,                              # j[6]

            len_,                                   # j[7] len
            i                                       # j[8] path
            ])

    np.set_printoptions(precision=2)
    np.set_printoptions(suppress=True)
    #print("#            (eV)               (eV)")
    #print("#train_frac testmin (step)    trainmin (step) path")
    #print(out)
    #print(out[out[:,0].argsort()])
    out=sorted(out,key=lambda x: x[0])
    out2=sorted(out2,key=lambda x: x[0])
    #print(out2)
    #print()
    #for j in out2:
    #    print(j)
    #print('-----')
    #l=out2
    #print('\n'.join(['%i: %s' % (n, l[n]) for n in xrange(len(l))]))
    #print('My list:', *out2, sep='\n- ')
    #for idj,j in enumerate(out): # for every line
    #    print("%0.2f  %6.1f  (%4.0f) %5.1f  (%4.0f) [%4.0f]   %s"%(j[0],j[1],j[4],j[2],j[5],j[6],j[3]))
    print()
    print()
    pp()
    pp()
    print("#             || rmse  /  rmse        || rmse / rmse         || [stes]")
    print("#train        || test  /  train       || train/ test         || [total]")
    print("#frac         || min   /      (@step) || min  /      (@step) ||         path")
    pp()
    for idj,j in enumerate(out2): # for every line
        run = "    "
        NJC = "    "
        que = "    "
        path_=len(j) - 1
        epochs_=len(j) - 2
        #print('a',path_,j,'--->',j[11])
        #sys.exit()
        if j[path_].split("/learning-curve.out")[0] in checkpath_run:
            run = "(R) "
        if j[path_].split("/learning-curve.out")[0] in checkpath_que:
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

        print(run+NJC+"%0.2f  || %5.1f /%5.1f  (%4.0f) ||%5.1f /%5.1f (%4.0f)  || %5.1f /%5.1f (%4.0f) || [%4.0f]  %s"%(j[0],j[1],j[2],j[3],j[4],j[5],j[6],   j[7],j[8],j[9]      ,j[epochs_],j[path_]))

