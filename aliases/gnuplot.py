#!/usr/bin/env python
 # -*- coding: utf-8 -*-
from __future__ import print_function
import numpy as np
import os,sys,argparse,re
from subprocess import call

def help(p = None):
    string = ''' helptext '''
    p = argparse.ArgumentParser(description=string,
            formatter_class=argparse.RawTextHelpFormatter)
    p.add_argument("inputfile",nargs='+') #,help"name of the inputfile(s)")
    p.add_argument('-mc','--max_columns',required=False, type=int,default=False, help="plot maximally first x columns")
    p.add_argument('-ll', '--log_log', action='store_true', default=False,help='make a x and y logarithmic (log log plot).')
    p.add_argument('-lx', '--log_x', action='store_true', default=False,help='make x axis logarithmic')
    p.add_argument('-ly', '--log_y', action='store_true', default=False,help='make y axis logarithmic')
    p.add_argument('-v','--verbose', help='verbose', action='count', default=False)

    p.add_argument('-xmin','--xmin' , required=False, action='append', type=float,
            help=argparse.SUPPRESS, default=False)
    p.add_argument('-xmax','--xmax' , required=False, action='append', type=float,
            help=argparse.SUPPRESS, default=False)
    p.add_argument('-ymin','--ymin' , required=False, action='append', type=float,
            help=argparse.SUPPRESS, default=False)
    p.add_argument('-ymax','--ymax' , required=False, action='append', type=float,
            help=argparse.SUPPRESS, default=False)
    p.add_argument('-n','--noplot', help=argparse.SUPPRESS, action='count', default=False)
    return p

def set_args_defaults(args,inputfile):
    basename = os.path.basename(inputfile)
    if args.verbose:
        print("######################## set_args_defaults #######################")
        print('basename',basename)
        print('args.inputfile       ',args.inputfile)
        print('args.verbose         ',args.verbose)
        print('args.log_log     (in)',args.log_log)
        print('args.max_columns (in)',args.max_columns)
    if basename == "learning-curve.out":
        args.log_log = True
        args.max_columns = 2
    if args.verbose:
        print('args.log_log     (in)',args.log_log)
        print('args.max_columns (in)',args.max_columns)
        print("######################## set_args_defaults #######################")
        print()
    return

def ca(text,newline = True,verbose=False):
    verbosity_level = 2
    global c
    if verbose > verbosity_level:
        print('c inside1:',c)
        print('c inside2:',text)
    c = c + text
    if newline == True:
        c = c + "\n"
    #print('c inside2:',c)
    return c

def gnuplot_plotline(inputfile,using=False,columns_tot=1):
    verbosity_level = 2
    if args.verbose > verbosity_level:
        print('--> inputfile        :',inputfile)
        print('--> using            :',using)
    global pl
    if pl == "": pl = "plot "

    if args.verbose > verbosity_level:
        print('--> xx using          ',using)
        print('--> xx pl (in)        ',pl)
    pladd = "\""+inputfile+"\" using "+using+" with linespoints"

    # legend
    if columns_tot == 1:
        #pladd = pladd + " notitle,"
        #pladd = pladd + " title \""+inputfile+"\","
        pladd = pladd + " title \""+inputfile+" "+str(using)+"\","
    else:
        pladd = pladd + " title \""+inputfile+" "+str(using)+"\","
    pl = pl + pladd
    if args.verbose > verbosity_level:
        print('--> xx pladd          ',pladd)
        print('--> xx pl (out)       ',pl)
    return pl

def gnuplot_defaults(args):
    ''' general settings '''
    global c
    c = "gnuplot --persist << EOF\n"

    if True:
        ca("set macros")
        # 1) change the default colors to more pleasant ones and make the lines
        # a little bit thicker
        #ca("set style line 1 lc rgb '#8b1a0e' pt 1 ps 1 lt 1 lw 2 # --- red")
        ca("set style line 1 lc rgb '#dd181f' lt 1 lw 2 pt 7   # red")
        ca("set style line 2 lc rgb '#5e9c36' pt 6 ps 1 lt 1 lw 2 # --- green")
        ca("set style line 3 lc rgb '#0060ad' lt 1 lw 2 pt 5   # blue")

        # 2) put the border more to the background by applying it only on the left
        # and bottom part and put it and the tics in gray
        ca("set style line 11 lc rgb '#808080' lt 1")
        ca("set border 3 back ls 11")
        ca("set tics nomirror")

        # 3) add a slight grid to make it easier to follow the exact position
        # of the curves
        ca("set style line 12 lc rgb '#808080' lt 0 lw 1")
        ca("set grid back ls 12")

        #set mxtics
        #set mytics
        #set style line 12 lc rgb '#ddccdd' lt 1 lw 1.5
        #set style line 13 lc rgb '#ddccdd' lt 1 lw 0.5
        #set grid xtics mxtics ytics mytics back ls 12, ls 13
    ca("#set yrange [0:10]")
    ca("#set xrange [0:10]")
    ca("set pointsize 2")
    if args.log_log:
        ca("set logscale xy")
    if args.log_x:
        ca("set logscale x")
    if args.log_y:
        ca("set logscale y")
    return

def y_min_max(args,column):
    cmin = column.min()
    #print('cmin',cmin)
    if type(args.ymin) == bool:
        args.ymin = cmin
    else:
        if cmin < args.ymin:
            args.ymin = cmin

    cmax = column.max()
    #print('cmax',cmax)
    if type(args.ymax) == bool:
        args.ymax = cmax
    else:
        if cmax > args.ymax:
            args.ymax = cmax
    return

def x_min_max(args,column):
    cmin = column.min()
    #print('cmin',cmin)
    if type(args.xmin) == bool:
        args.xmin = cmin
    else:
        if cmin < args.xmin:
            args.xmin = cmin

    cmax = column.max()
    #print('cmax',cmax)
    if type(args.xmax) == bool:
        args.xmax = cmax
    else:
        if cmax > args.xmax:
            args.xmax = cmax
    return

def gnuplot_plot(args):
    ################################################
    # make the plot
    ################################################
    verbosity_level = 1
    if args.verbose > verbosity_level:
        print("######################## gnuplot_plot      #######################")
    global pl
    pl = ""
    for idx,inputfile in enumerate(args.inputfile):
        if idx == 0:
            set_args_defaults(args,inputfile)
            gnuplot_defaults(args)
        if args.verbose > verbosity_level:
            print('##########################################')
            print('inputfile',inputfile)
        input = np.loadtxt(inputfile)
        if args.verbose > verbosity_level:
            print('input',input)
            print('input.shape',input.shape)
        if len(input.shape) == 1:
            using = "1"
            if args.verbose > verbosity_level:
                print('input',input,input.min(),input.max())
            y_min_max(args,column=input)
            x_min_max(args,column=input)
            text = gnuplot_plotline(inputfile,using = using,columns_tot = 1)
            ca(text)
        elif len(input.shape) == 2:
            x_min_max(args,column=input[:,0])
            columns = input.shape[1] - 1
            if args.verbose > verbosity_level:
                print('args.max_columns',args.max_columns)
            if type(args.max_columns) != bool and args.max_columns < columns:
                columns = args.max_columns
            for i in np.arange(columns)+2:
                if args.verbose > verbosity_level:
                    print('ixx',i)
                using = "1:"+str(i)
                #using = "1:(\$"+str(i)+"*27211.386)" # has to be a default
                if args.verbose > verbosity_level:
                    print('input column:',input[:,i-1])
                y_min_max(args,column=input[:,i-1])
                if args.verbose > verbosity_level:
                    print('args.ymin',args.ymin)
                    print('args.ymax',args.ymax)
                text =  gnuplot_plotline(inputfile,using = using,\
                        columns_tot = columns)
                ca(text)
        if args.verbose > verbosity_level:
            print('##########################################')
            print()
            print()
    #    ca("plot \"plot.gnu\" using 1:2 with linespoints notitle,    \"plot.gnu\" using 1:3 with linespoints notitle")
    ca("\nEOF")

    if args.verbose > 2:
        print("######################## show command BEFORE #####################")
        print(c)

    global c
    c = re.sub(r"#set yrange .*", "set yrange ["+str(args.ymin)+":"+str(args.ymax)+"]", c)

    c = re.sub(r"#set xrange .*", "set xrange ["+str(args.xmin)+":"+str(args.xmax)+"]", c)

    #c_new = c
    if args.verbose:
        print('args.ymin',args.ymin)
        print('args.ymax',args.ymax)
    if args.verbose:
        print("######################## show command      #######################")
        print(c)
        print("######################## show command      #######################")

    #"set style line 1 lc rgb '#8b1a0e' pt 1 ps 1 lt 1 lw 2 # --- red"
    #command = ["gnuplot --persist << EOF\nplot sin(x)\nEOF"]
    if args.noplot == False:
        call(c,shell=True)
    return

if __name__ == '__main__':
    p = help()
    args = p.parse_args()
    gnuplot_plot(args)


