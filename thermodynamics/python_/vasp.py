#!/usr/bin/python
import get, my, argparse, sys, os, db
###############################################################################
## definitions start
## V V V V V V V V V
###############################################################################
all_builtin =  dir() ## necessat to define all defs in this script
#noinspection PyBroadException
def testall():
    """
    -ta test all of the following optional arguments
    """
    out = alldefs()
    #    print "out:",out
    out.remove('alldefs')
    out.remove('testall')
    #    print "out;",out
    #    print ""
    for i in out:
        try:
            ha = str(out.index(i)+1)+"/"+str(len(out)).ljust(7)  +   i.ljust(45)
            hb = str(eval(i)()).ljust(70)
            print "OK " + ha + hb
        except: # catch *all* exceptions
            e = sys.exc_info()[1]
            f = str(e)
            #print hbe = sys.exc_info()[0]
            #            print "e:",e,type(e),str(e)
            #noinspection PyUnboundLocalVariable
            my.print_red_bold("-> " + ha + f.ljust(45))

    my.print_black_bold("-------------------------------> testall <-- finished")
    return ""
def alldefs(script=None):
    """
    -gdcs get all definitions of this script
    """
    return my.get_definitions_of_current_script(script)

#######################
## path definitions
#######################
def get_path_to_file():
    """
    -gp
    """
    pass


def potcar_electrons(path=None):
    """
    -pote
    """
    if path == None:
        path = get_path_to_file()

    electrons = my.readfile(path=path,line=2)
    #print electrons
    if my.is_number(electrons) != True:
        my.exit(error="how many electrons does "+str(path)+" have?")
    return float(electrons)




###############################################################################
## ^ ^ ^  ^ ^  ^ ^
## definitions stop
###############################################################################

############################################################################
## main/data/glensk/progs/vasp/TEST_1h_at_12_cores_VORLAGE/cmmd010@16cores@vasp.5.2.12_m3_comp-intelcompiler_run-intelcompiler_2
############################################################################
if __name__ == "__main__":

    p = argparse.ArgumentParser("TODO: all ERRORS to my.exit()")
    g = p.add_mutually_exclusive_group(required=True)  ## besser hinzufuegen da ansonsten so sachen wie -one 8 jj -two erlaubt sind
#    print  my.get_definitions_of_current_script_additionalinfo()
    for deff in my.get_definitions_of_current_script_additionalinfo():
#        print "--> deff:  ",deff,"1:",deff[0]
        definition = deff[0]
        short = deff[1]
        docstring = deff[2]
        arguments = deff[3]
#        print "--> definition:",definition
#        print "     --> short:",short
#        print "     --> help:",docstring
#        print "     --> type:",arguments
#        print ""

        ### add to parser
        if short != "":
            g.add_argument('--'+definition,short, dest='def_to_run', action='store_const',const=definition, help=docstring+arguments)
        else:
            g.add_argument('--'+definition, dest='def_to_run', action='store_const',const=definition, help=docstring+arguments)

    p.add_argument('arguments', nargs='*')
    args = p.parse_args()
#    print "def_to_run        :",args.def_to_run
#    print "arguments         :",args.arguments
#    print "len(arguments)    :",len(args.arguments)
#    print "args.arguments[:] :",args.arguments[:]
#    import inspect
#    print "args of this def  :",inspect.getargspec(eval(args.def_to_run))[0]

    if len(args.arguments) == 0:
    #        print "0:::",args.def_to_run
        print eval(args.def_to_run)()
    elif len(args.arguments) == 1:
#        print "1::::",args.arguments
#        print "2::::",args.arguments[0]
#        print ""
#        print eval(args.def_to_run)(",".join(args.arguments))
#        print eval(args.def_to_run)(outfile="hallo")  ## funktioniert
        print eval(args.def_to_run)(*args.arguments)
    elif len(args.arguments) > 1:
#        print "3:::::::::::",args.arguments
        print eval(args.def_to_run)(*args.arguments)
