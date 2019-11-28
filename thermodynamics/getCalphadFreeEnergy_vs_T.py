#!/usr/bin/python
import sys
import numpy as np
#from sympy import *
import os



try:  
       th = os.environ["thermodynamics"]
except KeyError: 
       print "Please set the environment variable \"thermodynamics\" which point to your Thermodynamics folder"
       sys.exit(1)

calphad = str(th) + "/utilities/db_CALPHAD_PURE5_SGTE.TDB"

def is_number(s):
    try:
        float(s)
        return True
    except (ValueError, TypeError):
        return False


def is_int(s):
    try:
        int(s)
        return True
    except ValueError:
        return False


def is_float(s):
    try:
        float(s)
        return True
    except ValueError:
        return False


def run(befehl):
    """
    z.B. run('ls -l')
    z.B. run('frompath_stoich.sh')
    """
    import subprocess
    return subprocess.check_output(befehl, shell=True, stderr=subprocess.STDOUT)

def iroundup(x):
    x = x + 0.5
    """iround(number) -> integer
    Round a number to the next integer. 190.0 -> 190; 191.6 -> 192"""
    return int(round(x))


def convert_to_number(string):
    if is_int(string) == True:
        value = int(string)
    elif is_float(string) == True:
        value = float(string)
    else:
        import sys
        sys.exit("seems value", string, "is not a number")
    return value

############################################################################
## 0. all possible Elements
############################################################################
allelementslist = "Ac Ag Al Am Ar As At Au Bb Ba Be Bi Bb Cc Ca Cd Ce Cf Cl "\
                  "Cm Co Cr Cs Cu D1 Dt D2 Dy Er Es Eu Ff Fe Fm Fr Ga Gd Ge "\
                  "Hh He Hf Hg Ho Ii In Ir Kk Kr La Li Lu Mg Mn Mo Nn Na Nb "\
                  "Nd Ne Ni Np Oo Os Pp Pa Pb Pd Pm Po Pr Pt Pu Ra Rb Re Rh "\
                  "Rn Ru Ss Sb Sc Se Si Sm Sn Sr T1 T2 Ta Tb Tc Te Th Ti Tl "\
                  "Tm Uu Vv Ww Xe Yy Yb Zn Zr"
allelements = allelementslist.split()
#print allelements

############################################################################
## 1. Get the element
############################################################################


def help():
    print ""
    print("       This skript greps the unary energy functions defined int he PURE5 SGTE database")
    print("       returnes their energy in [meV] as a function of temperature")
    print ""
    print("Usage: ")
    print(str(sys.argv[0]).split("/")[-1] + " [ELEMENT] [OPTIONS] (the element input is not case sensitive)")
    print ""
    print("[Elements]")
    print allelementslist
    print ""
    print("[OPTIONS]:")
    print "-h           print this help"
    print "-help        print this help"
    print "-j           output in [j/(mol K)] instead of [meV]"
    print "-a           only print algebraic expressions"

    quit()

if "-h" in str(sys.argv[:]) or "-help" in str(sys.argv[:]):
    help()

if len(sys.argv[:]) <= 1:
    help()

#for element in allelements:
element = str(sys.argv[1]).upper()
string = str(run("sed -n '/FUNCT GHSER" + str(element) + \
            "/,/\!/p' " + str(calphad)))


## 2. Get the temperature ranges (getrennt durch ;)
import re
freeene = []
alltemp = []
for num, i in enumerate(string.split(";")):
    #print "-"*40

    oneline = " ".join(i.splitlines()).replace("FUNCT", "")
    oneline = oneline.replace("LN", "np.log")
    #print oneline
    onelinetemp = re.sub('\s+', ' ', oneline).strip().split()[0:2]
    onelineterm = "".join(re.sub('\s+', ' ', oneline).strip().split()[2:])

    if is_number(onelinetemp[0]) == True:
        if is_number(onelinetemp[1]) == False:
            temp = convert_to_number(onelinetemp[0])
    if is_number(onelinetemp[1]) == True:
        if is_number(onelinetemp[0]) == False:
            temp = convert_to_number(onelinetemp[1])
    #print "num:",num
    #print "temp:",temp
    #print "term:",onelineterm
    freeene.append([temp, onelineterm])
    alltemp.append(temp)

mintemp = freeene[0][0]
maxtemp = freeene[-1][0]
tmelt = freeene[-2][0]
#print "alltemp:", alltemp
#print "mintemp:", mintemp
#print "maxtemp:", maxtemp
#print "tmelt:", tmelt
#print element, len(freeene)
all = list(enumerate(freeene))
kjpermol_to_mev = 10.364269
jpermol_to_mev = 0.010364269
mevperk_in_kb = 11.604506
for index in range(len(all)-1):
    temp_start = all[index][1][0]
    temp_end = all[index+1][1][0]
    ene_expression = all[index][1][1]

    #T = Symbol('T')
    #y = v
    #print "ex:", ene_expression
    #quit()
    #ene_der1 = ka.diff(T)
    #print ene_der1
    if "-a" in str(sys.argv[:]):
        print index, temp_start, temp_end , ene_expression
        continue
    def free():
        #return lambda T: eval(ene_expression)
        return lambda T: eval(ene_expression)

    for T in range(iroundup(temp_start),iroundup(temp_end)):
        if "-j" in str(sys.argv[:]):
            print T, free()(T)
        else:
            print T, free()(T)*jpermol_to_mev
