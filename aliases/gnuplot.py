#!/usr/bin/env python
 # -*- coding: utf-8 -*-
from __future__ import print_function
import numpy as np
import os,sys,argparse,re,socket
from subprocess import call

def help(p = None):
    string = '''
    gnuplot.py KMC_analyze -c 7 8 9
    gnuplot.py KMC_analyze -c 2
    '''
    p = argparse.ArgumentParser(description=string,
            formatter_class=argparse.RawTextHelpFormatter)
    p.add_argument("inputfile",nargs='+') #,help"name of the inputfile(s)")
    p.add_argument('-mc','--max_columns',required=False, type=int,default=False, help="plot maximally first x columns")
    p.add_argument('-ps','--pointsize',required=False, type=float,default=1., help="pointsize")
    p.add_argument('-c','--columns',required=False, action='append',nargs='+', type=str, help="which columns to plot")
    p.add_argument('-ll', '--log_log', action='store_true', default=False,help='make a x and y logarithmic (log log plot).')
    p.add_argument('-lx', '--log_x', action='store_true', default=False,help='make x axis logarithmic')
    p.add_argument('-ly', '--log_y', action='store_true', default=False,help='make y axis logarithmic')
    p.add_argument('-a','--all',  help='print all columns as a function of first column', action='count', default=False)
    p.add_argument('-v','--verbose', help='verbose', action='count', default=False)
    p.add_argument('-s','--scale', help='scale axis', action='count', default=False)
    p.add_argument('-x11','--x11', help='ser term to x11', action='count', default=False)

    p.add_argument('-xlabel','--xlabel' , required=False, action='append', type=float,
            help=argparse.SUPPRESS, default=False)
    p.add_argument('-ylabel','--ylabel' , required=False, action='append', type=float,
            help=argparse.SUPPRESS, default=False)
    p.add_argument('-xmin','--xmin' , required=False, action='append', type=float,
            help=argparse.SUPPRESS, default=False)
    p.add_argument('-xmax','--xmax' , required=False, action='append', type=float,
            help=argparse.SUPPRESS, default=False)
    p.add_argument('-ymin','--ymin' , required=False, action='append', type=float,
            help=argparse.SUPPRESS, default=False)
    p.add_argument('-ymax','--ymax' , required=False, action='append', type=float,
            help=argparse.SUPPRESS, default=False)
    p.add_argument('-scale_y','--scale_y' ,required=False, action='append', type=float,
            help=argparse.SUPPRESS, default=1.0)
    p.add_argument('-n','--noplot', help=argparse.SUPPRESS, action='count', default=False)

    p.add_argument('-ls','--line_style', choices=["lines","points","linespoints"], default="line", help='plot lines as')
    return p

def is_int(x):
    try:
        a = float(x)
        b = int(a)
    except ValueError:
        return False
    else:
        return b==a

def get_default_y_label(args,inputfile,column):
    basename = os.path.basename(inputfile)
    print('bn',basename)
    if basename == "KMC_analyze":
        if column == 2: args.ylabel = "time (ps) : The time psed simulation time."
        if column == 3: args.ylabel = "time (ps) : The time of the current timestep"
        if column == 4: args.ylabel = "ecurr (meV)"
        if column == 5: args.ylabel = "cdf"
        if column == 6: args.ylabel = "random number 2"
    if basename in [ "KMC_analyze_np_analyze.dat", "KMC_analyze_akmc_ipi.dat" ]:
        if column == 2: args.ylabel = "residency time <400 MD steps>"
        if column == 3: args.ylabel = "Si atoms in 1NN"
        if column == 4: args.ylabel = "Mg atoms in 1NN"
        if column == 5: args.ylabel = "Si+Mg atoms in 1NN"
        args.xlabel = "MD step"
    print('kk',args.ylabel)
    #sys.exit()
    return


def set_args_defaults(args,inputfile):
    basename = os.path.basename(inputfile)
    if args.verbose:
        print("######################## set_args_defaults #######################")
        print('basename',basename)
        print('args.inputfile       ',args.inputfile)
        print('args.verbose         ',args.verbose)
        print('args.log_log     (in)',args.log_log)
        print('args.max_columns (in)',args.max_columns)
    if basename in ["learning-curve.out", 'learning-curve-runner.out']:
        args.log_log = True
        args.max_columns = 2
        args.columns = [2,3,4,5]
        args.line_style = "linespoints"
        #args.scale = True
        args.xlabel = "epochs"
        args.ylabel = "RMSE (meV/at)"
    if basename == "learning-curve.out":
        args.scale_y = [27211.386,27211.386,51422.063,51422.063]   # hartree to meV , here, as many as columns
    if basename == "KMC_analyze":
        args.xlabel = "step"
        pass
    if socket.gethostname() == "mac":
        args.x11 = True  # since aquaterm does open figurs in the background (but the first one)
    if args.verbose:
        print('args.log_log     (out)',args.log_log)
        print('args.max_columns (out)',args.max_columns)
        print('args.scale       (out)',args.scale)
        print("######################## set_args_defaults #######################")
        print()
    return

def args_show(args):
    if args.verbose:
        print()
        print("######################## args_show(args) #########################")
        print('args.inputfile   ',args.inputfile)
        print('args.max_columns ',args.max_columns)
        print('args.all         ',args.all)

        print('args.log_log     ',args.log_log)
        print('args.log_x       ',args.log_x)
        print('args.log_y       ',args.log_y)
        print('args.max_columns ',args.max_columns)
        print('args.scale       ',args.scale)
        print('args.verbose     ',args.verbose)
        print('args.x11         ',args.x11)
        print('args.noplot      ',args.noplot)
        print('args.scale_y     ',args.scale_y)
        print()
        print('args.xmin        ',args.xmin)
        print('args.xmax        ',args.xmax)
        print('args.ymin        ',args.ymin)
        print('args.ymax        ',args.ymax)
        print()
        print('args.line_style  ',args.line_style)
        print('args.columns     ',args.columns,"type:",type(args.columns))
        print("######################## args_show(args) #########################")
        print()
        return

def gnuplot_defaults_labels(args):
    if args.xlabel:
        ca("set xlabel \""+args.xlabel+"\"")
    if args.ylabel:
        ca("set ylabel \""+args.ylabel+"\"")
    return

def gnuplot_defaults(args):
    ''' general settings '''
    global c
    c = "gnuplot --persist << EOF\n"

    ca("set macros")
    if True:
        if False:
            ca("set macros")
        # 1) change the default colors to more pleasant ones and make the lines
        # a little bit thicker
        #ca("set style line 1 lc rgb '#8b1a0e' pt 1 ps 1 lt 1 lw 2 # --- red")

        # pt 7 --> circle
        # pt 5 --> square
        # pt 9 --> triangle
        #ca("set style line 1 lc rgb '#0060ad' lt 1 lw 2 pt 5   # blue")
        if False:
            ca("set style line 1 lc rgb '#dd181f' lt 1 lw 2 pt 7   # red")
            #ca("set style line 2 lc rgb '#5e9c36' lt 1 lw 2 pt 6 ps 1 # --- green")
            ca("set style line 2 lc rgb '#5e9c36' lt 1 lw 2 pt 5   # green")
            ca("set style line 3 lc rgb '#0060ad' lt 1 lw 2 pt 9   # blue")

        # 2) put the border more to the background by applying it only on the left
        # and bottom part and put it and the tics in gray
        if False:
            ca("set style line 11 lc rgb '#808080' lt 1")
            #ca("set border 3 back ls 11")
            ca("set tics nomirror")

        # 3) add a slight grid to make it easier to follow the exact position
        # of the curves
        if True:
            ca("")
            ca("## linestyle for the grid")
            ca("set style line 12 lc rgb '#808080' lt 0 lw .1")
            ca("set grid back ls 12")

        ca("")
        ca("## font size for legend")
        #ca("set terminal enhanced font \'Verdana,10\'")
        #ca('set term x11 font "sans,12"')
        ca('set term x11 font "Times-New-Roman,12"')
        #ca('set term aqua font "Times-New-Roman,12"')
        #ca("set terminal x11 font \'Verdana,10\'")
        #ca("set key font \",10\"")
        #set mxtics
        #set mytics
        #set style line 12 lc rgb '#ddccdd' lt 1 lw 1.5
        #set style line 13 lc rgb '#ddccdd' lt 1 lw 0.5
        #set grid xtics mxtics ytics mytics back ls 12, ls 13
    ca("")
    ca("#set yrange [0:10]")
    ca("#set xrange [0:10]")
    ca("")
    ca("set pointsize "+str(args.pointsize))  # how large should the symbols be?

    if args.x11:
        #ca("set terminal x11")
        #ca("set terminal aqua")
        pass
    if args.log_log:
        ca("set logscale xy")
    if args.log_x:
        ca("set logscale x")
    if args.log_y:
        ca("set logscale y")
    return

def y_min_max(args,column,inputfile=False):
    cmin = column.min()
    #print('cmin',cmin)
    if type(args.ymin) == bool:
        args.ymin = cmin
    else:
        if cmin < args.ymin:
            args.ymin = cmin

    cmax = column.max()
    if inputfile == "learning-curve.out":
        cmax = column[1:].max()

    #print('cmax',cmax)
    if type(args.ymax) == bool:
        args.ymax = cmax
    else:
        if cmax > args.ymax:
            args.ymax = cmax

    return

def x_min_max(args,column,inputfile=False):
    cmin = column.min()
    #print('cmin',cmin)
    if type(args.xmin) == bool:
        args.xmin = cmin
    else:
        if cmin < args.xmin:
            args.xmin = cmin

    cmax = column.max()
    if inputfile == "learning-curve.out":
        cmax = cmax*1.4
    #print('cmax',cmax)
    if type(args.xmax) == bool:
        args.xmax = cmax
    else:
        if cmax > args.xmax:
            args.xmax = cmax
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
    if args.line_style == "line":
        pladd = "\""+inputfile+"\" using "+using+" with line"
    if args.line_style == "linespoints":
        pladd = "\""+inputfile+"\" using "+using+" with linespoints"
    if args.line_style == "points":
        #pladd = "\""+inputfile+"\" using "+using+" with points"
        pladd = "\""+inputfile+"\" using "+using+" w p pt 7" # ls 2"

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

def column_scale_y(column_index,gocolumns,args):
    if args.verbose:
        print("## column_scale_y: column_index:",column_index)
        print("## column_scale_y: gocolumns   :",gocolumns,type(gocolumns))
        print("## column_scale_y: args.scale_y:",args.scale_y,type(args.scale_y))
    if isinstance(args.scale_y, (int, long, float, complex)):
        return args.scale_y
    elif type(args.scale_y) == list:
        return args.scale_y[column_index]
    else:
        print("## column_scale_y: column_index:",column_index,type(column_index))
        print("## column_scale_y: gocolumns   :",gocolumns,type(gocolumns))
        print("## column_scale_y: args.scale_y:",args.scale_y,type(args.scale_y))
        sys.exit('type of args.scale_y not known')
    return

def string_to_using__columns(i,verbose=False):
    if args.verbose:
        print('&& string_to_using__columns: && i',i)
    if is_int(i):
        return "1:"+str(i)
    # first check if there is an ":" if it is, number before is the x achsis
    #print('yes')
    #print('yes',int(i))

    ba = i.split(":")
    if args.verbose:
        print('&& string_to_using__columns: && ba',ba,len(ba))

    ## split x and y columns
    if len(ba) == 2:
        x = ba[0]
        y = ba[1]
    else:
        x = "1"
        y = ba[0]
    print('&& string_to_using__columns: && x',x)
    print('&& string_to_using__columns: && y',y)
    bb = x+":("+y+")"
    if is_int(y):
        bb = x+":($"+y+")"
    print('&& string_to_using__columns: && bb',bb)

    ## make $ to \$
    c = bb.split("$")
    print('&& string_to_using__columns: && c',c)
    using = d = "\$".join(c)
    print('&& string_to_using__columns: && d',d)
    # works:
    #   -c '1:$2'
    #
    # does not work:
    #   -c '1:2'
    #sys.exit()
    return using



def gnuplot_plot(args):
    global c
    ################################################
    # make the plot
    ################################################
    verbosity_level = 1
    global pl
    pl = ""
    ###############################################
    # go over all inputfiles
    ###############################################
    for idx,inputfile in enumerate(args.inputfile):
        if idx == 0:
            set_args_defaults(args,inputfile)
            gnuplot_defaults(args)
        if args.verbose > verbosity_level:
            print("######################## gnuplot_plot begin ######################")
            print('%% inputfile',inputfile)

        #input = np.loadtxt(inputfile)
        from numpy import genfromtxt
        input = genfromtxt(inputfile,filling_values=99)
        print('input',input)
        if args.verbose > verbosity_level+1:
            print('input',input)
        if args.verbose > verbosity_level:
            print('%% input.shape',input.shape)
            print('%% len(input.shape)',len(input.shape))
        if len(input.shape) == 1:
            ###############################################
            # in case only one colume, print it as function of [1,2,3, ...]
            ###############################################
            using = "1"
            if args.scale_y != 1.0:
                using = "$1*"+str(args.scale_y)+")" # has to be a default
            #if False:
            y_min_max(args,column=input,inputfile=inputfile)
            x_min_max(args,column=input,inputfile=inputfile)
            if args.verbose > verbosity_level:
                print("## @@@@@@@@@@@@@@@@@ one column (begin) @@@@@@@@@@@@@@@@@@@@@@")
                print('## args.xmin',args.xmin)
                print('## args.xmax',args.xmax)
                print('## args.ymin',args.ymin)
                print('## args.ymax',args.ymax)
            text = gnuplot_plotline(inputfile,using = using,columns_tot = 1)
            ca(text)
            if args.verbose > verbosity_level:
                print("## @@@@@@@@@@@@@@@@@ one column (end) @@@@@@@@@@@@@@@@@@@@@@@@")
        elif len(input.shape) == 2:
            ###############################################
            # go over all columns,
            # take column 0 as x-column
            ###############################################
            x_min_max(args,column=input[:,0],inputfile=inputfile)
            columns = input.shape[1] - 1
            if args.verbose > verbosity_level:
                print('%% args.max_columns',args.max_columns)
            if type(args.max_columns) != bool and args.max_columns < columns:
                columns = args.max_columns
            if args.all == False:
                print("%% args.all == False --> columns = 1")
                columns = 1
            gocolumns = np.arange(columns)+2
            if args.columns:
                gocolumns = np.array(args.columns)
            if args.verbose > verbosity_level:
                print('%% columns  ',columns)
                print('%% gocolumns',gocolumns)
            for idxc,i in enumerate(gocolumns): # gocolumns [2, 3, 4, 5]
                ###############################################
                # go over all columns,
                # take first column as x-column
                ###############################################
                if args.verbose > verbosity_level:
                    print("## @@@@@@@@@@@@@@@@@ TWO columns (begin) @@@@@@@@@@@@@@@@@@@@@")
                    print('## i',i)
                if type(args.columns) == list and i not in args.columns:
                    print('## --> skipping column',i,type(i),'since not in args.columns',args.columns)
                    continue
                else:
                    print('## --> adding   column',i,type(i),'since in args.columns',args.columns)

                using = string_to_using__columns(i,verbose=args.verbose)
                if args.verbose > verbosity_level:
                    print('## && using result (1):',using)

                #using = "1:($"+str(i)+")"
                #using = "1:"+str(i)
                #using = str(i)
                scale_y = column_scale_y(idxc,gocolumns,args)
                if scale_y != 1.0:
                    using = "1:(\$"+str(i)+"*"+str(scale_y)+")" # has to be a default

                if args.verbose > verbosity_level:
                    print('## using (2):',using)

                if args.verbose > verbosity_level+1:
                    print('## input column:',input[:,i-1])
                if False:
                    y_min_max(args,column=input[:,i-1],inputfile=inputfile)
                if args.verbose > verbosity_level:
                    print('## args.xmin',args.xmin)
                    print('## args.xmax',args.xmax)
                    print('## args.ymin',args.ymin)
                    print('## args.ymax',args.ymax)
                    print('## columns  ',columns)
                get_default_y_label(args,inputfile,i)
                if args.verbose > verbosity_level:
                    print('##xx args.ylabel',args.ylabel)
                    print('##xx args.xlabel',args.xlabel)
                gnuplot_defaults_labels(args)
                text =  gnuplot_plotline(inputfile,using = using,\
                        columns_tot = columns)
                if args.verbose > verbosity_level:
                    print('## text:',text)
                ca(text)
                if args.verbose > verbosity_level:
                    print("##",text)
                    print("## @@@@@@@@@@@@@@@@@ TWO columns (end) @@@@@@@@@@@@@@@@@@@@@@@")
        if args.verbose > verbosity_level:
            print("######################## gnuplot_plot end  #######################")
            print()
    #    ca("plot \"plot.gnu\" using 1:2 with linespoints notitle,    \"plot.gnu\" using 1:3 with linespoints notitle")
    ca("\nEOF")

    if args.verbose > 2:
        print("######################## show command BEFORE #####################")
        print(c)
    if args.log_log or args.log_x:
        if args.scale: c = re.sub(r"#set yrange .*", "set yrange ["+str(args.ymin*args.scale_y)+":"+str(args.ymax*args.scale_y)+"]", c)

    if args.log_log or args.log_x:
        if args.xmin <= 0.0:
            args.xmin = 0.1
        if args.scale: c = re.sub(r"#set xrange .*", "set xrange ["+str(args.xmin)+":"+str(args.xmax)+"]", c)
    #c = re.sub(r"#set xrange .*", "set xrange ["+str(4000)+":"+str(7000)+"]", c)

    #c_new = c
    if args.verbose:
        print()
        print("######################## show command (begin) ####################")
        print(c)
        print("######################## show command (end) ######################")
        print()

    #"set style line 1 lc rgb '#8b1a0e' pt 1 ps 1 lt 1 lw 2 # --- red"
    #command = ["gnuplot --persist << EOF\nplot sin(x)\nEOF"]
    if args.noplot == False:
        call(c,shell=True)
    return

if __name__ == '__main__':
    p = help()
    args = p.parse_args()
    if type(args.columns) == list:
        columns_new = []
        for i in args.columns[0]:
            try:
                columns_new.append(int(i))
            except:
                columns_new.append(i)
        #print('cn',columns_new)
        args.columns = columns_new
    #print('ag',args.columns)
    #sys.exit()
    args_show(args)
    #sys.exit()
    gnuplot_plot(args)
    args_show(args)
