#!/usr/bin/python
import sys
import numpy as np
import re
import os


def printred(var):
    print '\033[31m' + str(var) + '\033[0m'


def printhelp():
    printred("Usage: " + os.path.basename(sys.argv[0]) + " [inputfiles] [OPTIONS]")
    printred("  e.g. " + os.path.basename(sys.argv[0]) + " Fqh_ExactFreqs_1_13.82658525 Fqh_ExactFreqs_1_14.94411775")
    printred("  e.g. " + os.path.basename(sys.argv[0]) + " 3.*/Fqh_MeshFreqs_1_*")
    printred("  e.g. " + os.path.basename(sys.argv[0]) + " Fqh_fromExactFreqs_3_63.{22,55,88}")
    print ""
    print "INPUTFILES:"
    print("         - inputfiles should have following format:")
    print("                1st column: Temperature [K]")
    print("                2nd column: Helmholz Free Energy [meV]")
    print ""
    print("         - inputfiles naming convention: anything_NUMATOMS_VOLUMEperATOM in [Angstrom^3]")
    print("                         e.g. Fvib_fromExactFreqs_1_21.345234565")
    print("                         e.g. ANYNAME_32_16.0023")
    print ""
    print "OPTIONS:"
    print "   -h          print this help"
    print "   -d [2,3,4]  Degree of the fitting polynomial, default is 2"
    quit()


if len(sys.argv) < 2:
    printred("ERROR inputfiles missing!")
    print ""
    printhelp()

if sys.argv[1] == "--help" or sys.argv[1] == "-h":
    printhelp()

volumefolder = sys.argv[1:]
deg = 2
arguments = volumefolder
#print "V1", volumefolder
if "-d" in arguments:
    d_index = arguments.index("-d")
    deg_index = d_index + 1
    #print "deg_place", deg_place
    deg = arguments[deg_index]
    #print "-d found", deg
    if int(deg) == 1 or int(deg) == 2 or int(deg) == 3:
        deg = int(deg)
        del volumefolder[d_index]
        del volumefolder[d_index]
    else:
        sys.exit("-d option can only be 1 or 2 or 3; INPUT -d " + str(deg))


#print "VL", volumefolder
temps = np.array(np.loadtxt(volumefolder[1]))[:, 0]
nr_volumes = len(volumefolder)
ene = np.empty(shape=(len(temps), nr_volumes), dtype=float)


volumes = []
for i, j in enumerate(volumefolder):
    inputarray = np.loadtxt(j)

    # check if temperatures all all the same
    temps_check = inputarray[:,0]
    if all(temps_check!=temps):
        sys.exit("temperatures different in " + str(sys.argv[1]) + " and " + str(j))

    # check if volume is ok
    #print i, j, re.findall("[0-9.]*$", j)
    if len(re.findall("[0-9.]*$", j)) != 2:
        sys.exit("Filename wrong format: " + str(j))
    volume = re.findall("[0-9.]*$", j)[0]
    output_filename = j.split(volume)[0]
    volumes.append(float(volume))
    print('volumes:',volumes)

    # append array
    new_col = inputarray[:, 1]
    ene[:, i] = new_col
    #print new_col

# print ene
print("VL", volumes)
volumes=np.array(volumes)**3./4.
print("VL", volumes)
out = np.empty(shape=(len(temps), deg + 2), dtype=float)
for index, temp in enumerate(temps):
    # print "vol:",volumes
    # print "ene:",ene[index]
    coefs_absteigend = np.polyfit(volumes, ene[index], deg=deg)
    coefs = coefs_absteigend[::-1]
    out[index, 0] = temp
    #print "-->",out[index]
    out[index, 1:] = coefs

np.savetxt(output_filename + 'surface', out, fmt='%1.15f')   # x,y,z equal sized 1D arrays
