#!/usr/bin/env python

from __future__ import print_function
import numpy as np
import glob,sys,os,argparse,subprocess,time
import myutils as my
from myutils import ase_calculate_ene

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
    p.add_argument('-ex_c44'    ,'--ex_c44'            , action='store_true', default=False, help='execute part c44')
    p.add_argument('-ex_kmc57'  ,'--ex_kmc57'            , action='store_true', default=False, help='execute part kmc57')
    p.add_argument('-ex_test'  ,'--ex_test'            , action='store_true', default=False, help='execute part test')
    p.add_argument('-ex_train'  ,'--ex_train'            , action='store_true', default=False, help='execute part train')
    return p

p = help()
args = p.parse_args()
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

print('>> (3) my.q() ....')
que_id_all,que_stat_all,que_folder_all= my.q(args.verbose)


##########################################################################################
print('>> (4) go over all learning curve_files ....')
##########################################################################################
for i in all_learning_curve_files:
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
    pot.get(exit=False)
    pot.get_my_assessments()  # gets kmc57_{b,l}, train_{b,l}, test_{b,l}
    has_outliers, outliers_epochs, outliers_diff = pot.get_my_assessments_check_outliers(verbose=args.verbose)
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
        ask = False
        if len(elastic_c44_all) > 500:
            if elastic_c44_all[300][1] > 40 and elastic_c44_all[500][1] > 39.9:
                ask = True
                print('ask1',elastic_c44_all[300][1],elastic_c44_all[500][1])
        if len(elastic_c44_all) > 600 and elastic_c44_all[600][1] > 39.65:
            print('ask2',elastic_c44_all[600][1])
            aks = True
        if ask == True:
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
    rnd                 = my.inputnn_get_random_seed(pot.inputnn)
    testf               = round(my.inputnn_get_testfraction(pot.inputnn),2)
    train_fraction      = round(1.-testf,2)
    nodes_short         = my.inputnn_get_nodes_short(pot.inputnn,as_string=True)
    activation_short    = my.inputnn_get_activation_short(pot.inputnn)
    nn                  = nodes_short+"__"+activation_short
    if os.path.isdir(folder+"/kmc"):
        nn = nodes_short+"**"+activation_short
    pot_elements, [al,mg,si]= my.inputnn_get_atomic_symbols_and_atom_energy_list(pot.inputnn)
    al                  = round(int(al),2)
    mg                  = round(int(mg),2)
    si                  = round(int(si),2)
    lc                  = my.n2p2_runner_get_learning_curve(pot.inputnn)
    if args.verbose > 1:
        print('---lc---')
        print(lc)
    epochs_max          = len(lc[:,1])-1
    agree               = "?"
    ol                  = ""
    if has_outliers == True:
        ol = "! OUTL ! "

    if args.execute or args.ex_c44 or args.ex_kmc57 or args.ex_test or args.ex_train:
        print('learning_curve_file folder:',folder.ljust(40),os.getcwd())
        hier = os.getcwd()
        os.chdir(folder)
        add = ""
        if args.execute: args.ex_c44 = args.ex_kmc57 = args.ex_test = args.ex_train = True
        if args.ex_c44: add = add + ' -ex_c44 '
        if args.ex_kmc57: add = add + ' -ex_kmc57 '
        if args.ex_test: add = add + ' -ex_test '
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

    for eidx, epoch in enumerate(go_through):
        if epoch in pot.assessed_epochs:
            test    = round(pot.assessed_test[eidx],2)
            train   = round(pot.assessed_train[eidx],2)
            kmcbl   = round(pot.assessed_kmc57[eidx],2)
            c44     = round(pot.assessed_c44[eidx],1)
        else:
            test    = round(lc[:,2][epoch],2)
            train   = round(lc[:,1][epoch],2)
            kmcbl   = '-'
            c44     = '-'

        f1 = trainminf_at_testmin    = round(lc[:,3][epoch],2)
        f2 = testminf_at_testmin     = round(lc[:,4][epoch],2)

        #print('kmcbl',kmcbl,'epoch',epoch,'epochs_max',epochs_max)
        #if kmcbl == '-':
        #    pot.print_variables_mypot(text=">> (4)",print_nontheless=True)
        if type(kmcbl) != str:
            #kmcbl = str(round(kmcbl,1)).ljust(4)
            kmcbl = str(round(kmcbl,1))
        stringout = str(ol)+str(que_status)+" %0.1f |"+agree+"%5.1f /%5.1f  (%4.0f) || %s %s %s || %5.1f /%5.1f |C %4s |K %3s |n %4.0f || [%4.0f] | %s | %8.0f | %s"
        #elementout = (    train_fraction     ,        test, train, epoch,      al,mg,si,   f1,     f2,    c44, kmcbl,   ist ,    epochs_max,   nn,  rnd,   path)
        elementout = (    train_fraction     ,          test, train, epoch ,    al,mg,si,   f1,     f2,    c44, kmcbl,   0 ,    epochs_max,   nn,  rnd,   path)

        ##################################
        # get right color
        ##################################
        conv_unconv = "unconv"
        ## erstmak die komischen aussortieren
        if (test+train)/2. < 0.1 or al == 0.0:
            takecolor = "blue"  # wiered
        #elif c44 < 10:
        #    takecolor = "blue"  # wiered
        elif al < 3:  # old
            takecolor = 'orange'
        elif (test+train)/2. > 4. or (f1+f2)/2. >= 35.:
            takecolor = "red"
        elif (test+train)/2. < 4.0 and (f1+f2)/2. < 35.:
            takecolor = "green"
            conv_unconv = "conv"
        else:
            takecolor = "white"

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



        #print(stringout%elementout)
end_time = time.time()
print("TIME:",end_time - start_time)

