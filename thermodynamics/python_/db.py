#!/usr/bin/python
import get, my, argparse, sys, os
all_builtin =  dir() ## necessat to define all defs in this script

###############################################################################
## definitions start
## V V V V V V V V V
###############################################################################
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
def read_database(path_to_db=None,element_or_string=None,path=None):
    """
    -rd
    """
    if path_to_db == None: my.exit(error=" no path to database defined")
    my.checkfile(path_to_db)

    if element_or_string == None and path == None:
        my.exit(error=" please define either element or path")

    if path != None and element_or_string != None:
        my.exit(error=" please define either element or path")

    if element_or_string == None and path != None:
        element_or_string = get.from_path_string_element(path=path)
        element_or_string = str(element_or_string)

    if path == None and element_or_string != None:
        element_or_string = str(element_or_string)


    #    print "element: ",element_or_string

    file = my.readfile(path_to_db)
#    print "file: ",file

    out = my.grep("^"+element_or_string,file,options="-i")
#    quit()
#    import re
#    out=[]
#    for line in file:
#        print "line: ",line,len(line)
#        add=re.search(element_or_string,line)
#        if add != None:
#            add = re.search(element_or_string,line).group()
#            out.insert(1,line)
#    print "GR:",re.search(element_or_string,file[1])
#    print "OUT: ",out
#    print "len(OUT): ",len(out)
    if out == []: return None

    ### did we find one line or several?
    if len(out) != 1:
        return my.exit(error=" more than one line: "+str(out))
    return out[0].split()[1]


def nbands_peratom(element=None,path=None): ## = "/home/glensk/db/nbands_occupied.dat",element=None):
    """
    -nb
    """
    if path == None:
        path = my.pwd()

    path_to_db = get.set_path_db()+"/nbands_occupied.dat"
    print "path_to_db: ",path_to_db
    return read_database(path_to_db=path_to_db,element_or_string=element,path=path)

def AtomicWeight(element=None):
    """
    -aw
    """
    if element == None:
        element = get.element_real()
    else:
        element = get.element_real(element)

    #print "el:",element
    w = my.run("getAtomicMassStandard.sh "+str(element))
    #w = "a"
    if my.is_number(w):
        weight = w
    #if element == "al":
    #    weight = 26.981539
    #elif element == "cu":
    #    weight = 63.546
    #elif element == "ti":
    #    weight = 47.867
    #elif element == "pt":
    #    weight = 195.080
    else:
        #print "3"
        weight = my.exit(error=" Still to do db.AtomcWeight")  ## TODO
    #print "ja"
    if my.is_number(weight) != True:
        my.exit(error=" Atomic weight of "+element+" wrong!")
    #print "jojo"
    return float(weight)

def MeltingPoint(element=None,path=None):
    """
    -mp
    """
    if element == None:
        element = get.element_real()
    else:
        element = get.element_real(element)

    #print "ele:",element
    w = my.run("getMeltingPoint.sh "+str(element)+" -r")
    if my.is_number:
        Tm = w
    #if element == "al":
    #    Tm = 934
    elif element == "si mg":  ## dass sollte erstmal nur fuer Mg2si temporaer sein
        Tm = 1346
    #elif element == "cu":
    #    Tm = 1360
    #elif element == "ti":
    #    Tm = 1942
    #elif element == "pt":
    #    Tm = 2042
    else:
        Tm = my.exit(error=" Still to do db.MeltingPoint")  ## TODO
    if my.is_number(Tm) != True:
        my.exit(error=" MeltingPoint of "+element+" is not a number!")
    return Tm

###############################################################################
## ^ ^ ^  ^ ^  ^ ^
## definitions stop
###############################################################################

############################################################################
## main
############################################################################
if __name__ == "__main__":

    p = argparse.ArgumentParser("TODO: all ERRORS to my.exit()")
    g = p.add_mutually_exclusive_group(required=True)  ## besser hinzufuegen da ansonsten so sachen wie -one 8 jj -two erlaubt sind

    for deff in my.get_definitions_of_current_script_additionalinfo():
    #    print "deff",deff,"1:",deff[0]
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
    #    print "args of this def  :",inspect.getargspec(eval(args.def_to_run))

    if len(args.arguments) == 0:
    #        print "0:::",args.def_to_run
        print eval(args.def_to_run)()
    elif len(args.arguments) == 1:
    #        print "1::::"
        print eval(args.def_to_run)(",".join(args.arguments))
    elif len(args.arguments) > 1:
    #        print "3:::::::::::",args.arguments
        print eval(args.def_to_run)(*args.arguments)
