#!/usr/bin/env python
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


def string_to_num_list(string):
    """ Turn a string into a list of int/float chunks.
        "4.Ang_300K_.2m-9m=.77k" -> [4.0, 300.0, 0.2, 9.0, 0.77] """

    def tryfloat(stringtest):
        try:
            return float(stringtest)
        except:
            return stringtest
    import re
    return [ tryfloat(c) for c in re.findall(r"[-+]?\d*\.\d+|\d+", string)]

def lsn(searchstring):
    """ Returns files sorted numerically
        searchstring: e.g. "*Ang_*K"            """
    import glob
    files = glob.glob(searchstring)
    return sorted((files), key=string_to_num_list)


class bcolors:
    ENDC = '\033[0m'
    red = '\033[31m'
    green = '\033[32m'
    yellow = '\033[33m'
    blue = '\033[34m'
    okblue = '\033[94m'
    pink = '\033[35m'
    blackbold = '\033[38m\033[1m'
    redbold = '\033[31m\033[1m'


def printred(var):
    print bcolors.red + str(var) + bcolors.ENDC


def printgreen(var):
    print bcolors.green + str(var) + bcolors.ENDC


def printyellow(var):
    print bcolors.yellow + str(var) + bcolors.ENDC


def printblue(var):
    print bcolors.blue + str(var) + bcolors.ENDC


def printokblue(var):
    print bcolors.okblue + str(var) + bcolors.ENDC


def printpink(var):
    print bcolors.pink + str(var) + bcolors.ENDC


def printblackbold(var):
    print bcolors.blackbold + str(var) + bcolors.ENDC


def printredbold(var):
    print bcolors.redbold + str(var) + bcolors.ENDC


def printredboldd(var):
    bcolors.redbold + str(var) + bcolors.ENDC


def printredbolddd(var):
    bcolors.blackbold + str(var) + bcolors.ENDC
