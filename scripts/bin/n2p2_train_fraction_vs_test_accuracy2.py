#!/usr/bin/env python

from __future__ import print_function
import numpy as np
import glob,sys,os,argparse,subprocess,time
import myutils as my
from myutils import ase_calculate_ene
from filecmp import cmp

start_time = time.time()

def help(p = None):
    string = '''
    n2p2_train_fraction_vs_test_accuracy2.py -f v4ag ppl
    '''
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
    p.add_argument('-ex_c44'    ,'--ex_c44'            , action='store_true', default=False, help='execute part c44')
    p.add_argument('-ex_kmc57'  ,'--ex_kmc57'            , action='store_true', default=False, help='execute part kmc57')
    p.add_argument('-ex_test'  ,'--ex_test'            , action='store_true', default=False, help='execute part test')
    p.add_argument('-ex_train'  ,'--ex_train'            , action='store_true', default=False, help='execute part train')
    p.add_argument('-f'  ,'--find'               , default=-1.,type=str, nargs='*',required=False, help="find (& restrict selectiion to) pattern in foldername")
    p.add_argument('-pnn'  ,'--pick_nn'          , default=-1.,type=str, nargs='*',required=False, help="find (& restrict selectiion to) pattern in nn as 24_24__t_t_l")
    p.add_argument('-pas'  , '--pick_amount_sructures'   , default=-1.,type=str, nargs='*',required=False, help="find (& restrict selectiion to) potentials with certain number of structures")

    return p

def gettesttrain(test,ljust_=4,greater=3.0,round_=1):
    #print('test',test)
    if type(test) == str:
        return test.ljust(ljust_)
    if test is None:
        return "-".ljust(ljust_)
    test_ = round(test,round_)
    test_ = str(test_).ljust(ljust_)
    if test > greater: test_ = my.printred(test_)
    else: test_ = my.printnormal(test_)
    return test_

def getf1f2(f1,greater=55,ljust_=5):
    f1_ = str(f1).ljust(ljust_)
    if f1 > greater: f1_ = my.printred(f1_)
    elif 35.0 < f1 < greater: f1_ = my.printorange(f1_)
    else: f1_ = my.printnormal(f1_)
    return f1_

p = help()
args = p.parse_args()
if args.verbose:
    my.print_args(args)
if args.from_que == True: args.from_subfolder = False

##########################################################################################
print('>> (1) searching learning-curve.out files....')
##########################################################################################
maxdepth = 2
if args.potential: maxdepth = 3
all_learning_curve_files = my.find_files(os.getcwd(),'learning-curve.out',maxdepth=maxdepth)
if args.potential:
    for i in all_learning_curve_files:
        if not "potential/learning-curve.out" in i: all_learning_curve_files.remove(i)

if args.verbose:
    print()
    print('>> (1) ################################## all_learning_curve_files ########')
    for i in all_learning_curve_files:
        print('>> (1) i',i)
    print('>> (1) ################################## all_learning_curve_files ########')
    print()

from utils_rename import list_sorted
all_learning_curve_files = list_sorted(all_learning_curve_files)
if args.find != -1:
    all_learning_curve_files_ = []
    for i in all_learning_curve_files:
        #print()
        use = True
        for j in args.find:
            #print(i,j in i)
            if j not in i:
                use = False
                #print('cont')
                continue
        #print(">>",i,use)
        if use == True:
            all_learning_curve_files_.append(i)
    all_learning_curve_files = all_learning_curve_files_
    # for i in all_learning_curve_files:
    #     print('i',i)
    # sys.exit()


print('>> (3) my.q() ....')
que_id_all,que_stat_all,que_folder_all= my.q(args.verbose)


##########################################################################################
print('>> (4) go over all learning curve_files ....')
##########################################################################################
print("#         ||    ENERGY (RMSE)    ||          || FORCES (RMSE)|| C44  || KMC  || [stes]")
print("#train    || test / train (@step)||          || train/ test  || C44  || std  || [total]")
print("#frac     ||      /       (@step)||          || min  /       || C44  ||      ||   path")

for nnidx,i in enumerate(all_learning_curve_files):
    learning_curve_filename = os.path.basename(i)             # 'learning-curve.out'
    folder = i.replace(learning_curve_filename, '')[:-1]


    #################################################################
    # check if still  in que
    #################################################################
    que_status = " "
    que_id = False
    que_path = False
    for que_idx,que_folder in enumerate(que_folder_all):
        if que_folder in folder:
            que_status = que_stat_all[que_idx]
            que_id     = que_id_all[que_idx]
            que_folder = que_folder_all[que_idx]
            break

    #################################################################
    # get the potential
    #################################################################
    pot = my.mypot(False,folder,use_different_epoch=False,verbose=args.verbose)
    pot.get(exit=False,showerrors=False)
    pot.get_my_assessments()  # gets kmc57_{b,l}, train_{b,l}, test_{b,l}

    nodes_short         = my.inputnn_get_nodes_short(pot.inputnn,as_string=True)
    activation_short    = my.inputnn_get_activation_short(pot.inputnn)
    if type(args.pick_nn) == list and nodes_short+"__"+activation_short not in args.pick_nn:
        continue
    nn                  = str(nodes_short+"__"+activation_short).ljust(12)
    if os.path.isdir(folder+"/kmc"):
        nn = str(nodes_short+"**"+activation_short)
    input_structures    = str(my.inputdata_get_nuber_of_structures(pot.inputnn))
    if type(args.pick_amount_sructures) == list and str(input_structures) not in args.pick_amount_sructures:
        continue
    input_structures    = str(my.inputdata_get_nuber_of_structures(pot.inputnn)).ljust(6)

    if len(pot.potepoch_all) > 0:
        epoch_last = pot.potepoch_all[-1]
    else:
        epoch_last = False
    ab_test  = folder+"/assess_test_"+str(pot.potepoch_bestteste)
    al_test  = folder+"/assess_test_"+str(epoch_last)
    ab_train = folder+"/assess_train_"+str(pot.potepoch_bestteste)
    al_train = folder+"/assess_train_"+str(epoch_last)


    #################################################################
    # check if to kill the job
    #################################################################
    if que_status == "R" and que_id != False and os.path.isfile(pot.potpath+'/elastic_c44_all.dat'):
        elastic_c44_all = np.loadtxt(pot.potpath+'/elastic_c44_all.dat')
        #print(len(elastic_c44_all),'potpath',pot.potpath)
        qel = False
        if len(elastic_c44_all) > 500:
            if elastic_c44_all[300][1] > 40 and elastic_c44_all[500][1] > 39.9:
                qel = "qdel1"
                print('qdel1',elastic_c44_all[300][1],elastic_c44_all[500][1])
        if len(elastic_c44_all) > 600 and elastic_c44_all[600][1] > 39.65:
            print('qdel2',elastic_c44_all[600][1])
            qel = "qdel2"
        if len(elastic_c44_all) > 1000 and elastic_c44_all[1000][1] > 38.5:
            print('qdel3',elastic_c44_all[1000][1],folder)
            qel = "qdel3"
        if qel:
            #print(">300",elastic_c44_all[300])
            print('killing id ',que_id,que_status,que_folder)
            print()
            kill_or_not = my.get_from_prompt_True_or_False("elastic constant at step 300 is still > 40 GPa, ("+str(elastic_c44_all[300][1])+") should I kill the job?")
            if kill_or_not == True:
                print('!! scancel',que_id," #",que_status,que_folder)
                subprocess.call(["scancel "+str(que_id)],shell=True)
            print()

    #################################################################
    # print to screen
    #################################################################
    path                = folder.replace(os.getcwd()+'/',"")
    rnd                 = str(my.inputnn_get_random_seed(pot.inputnn)).ljust(8)
    testf               = round(my.inputnn_get_testfraction(pot.inputnn),2)
    train_fraction      = round(1.-testf,2)
    pot_elements, [al,mg,si]= my.inputnn_get_atomic_symbols_and_atom_energy_list(pot.inputnn)
    al                  = int(al)
    mg                  = int(mg)
    si                  = int(si)


    inputnn = pot.inputnn
    if args.potential:
        folder_up = folder.split("/")[:-1]
        inputnn_up = "/".join(folder_up)+"/input.nn"
        same = cmp(pot.inputnn, inputnn_up)
        if same: inputnn = inputnn_up
    lc                  = my.n2p2_runner_get_learning_curve(inputnn)

    if args.verbose > 1:
        print('---lc---')
        print(lc)
    epochs_max          = len(lc[:,1])-1
    agree               = "?"
    ol                  = ""

    if args.execute or args.ex_c44 or args.ex_kmc57 or args.ex_test or args.ex_train:
        print('learning_curve_file folder:',folder.ljust(40),os.getcwd())
        hier = os.getcwd()
        os.chdir(folder)
        add = ""
        if args.execute: args.ex_c44 = args.ex_kmc57 = args.ex_test = args.ex_train = True
        if args.ex_c44: add   = add + ' -ex_c44 '
        if args.ex_kmc57: add = add + ' -ex_kmc57 '
        if args.ex_test: add  = add + ' -ex_test '
        if args.ex_train: add = add + ' -ex_train '
        print('add',add)
        #sys.exit()
        subprocess.call(["n2p2_get_potential_folder_from_nr.py "+add],shell=True)
        os.chdir(hier)
        #print('done')
        #sys.exit()


    print()
    go_through = pot.assessed_epochs
    #print('go_through (0)','len(lc)-1             ',len(lc)-1)
    #print('go_through (0)','pot.potepoch_bestteste',pot.potepoch_bestteste)
    #print('go_through (1)',go_through)
    if len(go_through) == 0:
        go_through = [pot.potepoch_bestteste]
    #print('go_through (2)',go_through)
    if go_through[-1] != len(lc)-1:
        go_through = np.append(go_through,len(lc)-1)
    #print('go_through (3)',go_through)


    ###########################################################################
    # check epochs to go trough
    ###########################################################################
    for eidx, epoch in enumerate(go_through):
        if epoch in pot.assessed_epochs:
            #test    = round(pot.assessed_test[eidx],2)
            #train   = round(pot.assessed_train[eidx],2)
            kmcbl   = pot.assessed_kmc57[eidx]
            c44     = pot.assessed_c44[eidx]
        else:
            kmcbl   = '-'
            c44     = '-'

        has_outliers,has_outliers_, outliers_epochs, outliers_idx,outliers_diff = pot.get_my_assessments_check_outliers(specific_epoch=epoch,verbose=args.verbose)
        if has_outliers and has_outliers_ == "!!O":
            #print(has_outliers,has_outliers_,outliers_diff,type(outliers_diff))
            #print(outliers_diff.max())
            #diffint(np.max(outliers_diff))
            has_outliers_ = has_outliers_+"_"+str(int(outliers_diff.max()))

        ##################################################################################
        # calculate the "??O" folder
        ##################################################################################
        if has_outliers_ == "??O" and False:
            hier = os.getcwd()
            os.chdir(folder)
            print('oe',outliers_epochs,folder)
            for doepoch in outliers_epochs:
                print('doepoch',doepoch)
                #subprocess.call(["n2p2_get_potential_folder_from_nr.py -ex_test"+add],shell=True)
                subprocess.call(["getEnergies_byLammps.py -p . -ctest -pe "+str(doepoch)],shell=True)
            os.chdir(hier)


        test    = round(lc[:,2][epoch],1)
        train   = round(lc[:,1][epoch],1)
        f1 = trainminf_at_testmin    = round(lc[:,3][epoch],1)
        f2 = testminf_at_testmin     = round(lc[:,4][epoch],1)


        ####################################################
        # format output
        ####################################################
        kmcbl_ = gettesttrain(kmcbl,ljust_=4,greater=3.0)
        test_ = gettesttrain(test,ljust_=4,greater=3.0)
        train_ = gettesttrain(train,ljust_=4,greater=3.0)
        c44_ = gettesttrain(c44,ljust_=4,greater=37)
        f1_ = getf1f2(f1)
        f2_ = getf1f2(f2)
        fstr = (str(que_status)+" "+str(has_outliers_)).ljust(9)
        epoch = "("+(str(epoch).ljust(4))+") ||"
        epochs_max_ = "["+(str(epochs_max).ljust(4))+"] |"
        if eidx == 0: path_ = path
        else: path_ = ""

        ####################################################
        # print output
        ####################################################
        print(fstr, train_fraction,"|",test_,"/", train_,  epoch ,   al,mg,si,"||",   f1_,"/",f2_,"|C",    c44_,"|K",  kmcbl_,"|n",   input_structures ,     epochs_max_,nn,"|",  rnd,"|",   path_, nnidx)

end_time = time.time()
print("TIME:",end_time - start_time)

