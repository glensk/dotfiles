#!/usr/bin/python
import sys
import argparse
import numpy as np

"""Small script to generate symmetry functions for n elements by combining them as intended, without repetitions. The symfunctions used are those generated previously.

It should be noted that the generated symmetry functions are always the same, unless the values inside the arrays are changed. For the future it is better to introduce an input that lets the user choose the dimension of the smallest radial symfunction to probe different environments. """


def main(elements, cutoff, N, rmin=0, optcutoff=False):
    el_list = []

    index = np.arange(N+1, dtype=float)
    shift_array = cutoff*(1./N)**(index/(len(index)-1))
    eta_array = 1./shift_array**2.

    for i in range(len(elements)):
        if elements[i] in ["H","He","Li","Be","B","C","N","Ni","O","F","Ne","Na","Mg","Al","Si","P","S","Cl","Ar","K","Ca","Sc","Ti","V","Cr","Mn","Fe","As","Ga"]:
            el_list.append(elements[i])
        else:
            raise ValueError('The element',elements[i],' is not in the list (you can add it in the program) or does not exist')
    for fel in elements:
        for sel in elements:
            print "# symfunctions for type %s 2 %s" %(fel, sel)
            for eta in eta_array:
                if optcutoff:
                    cutoff = min(np.sqrt(6.908/eta), 12.0) ### Variable cutoff depending on the gaussian
                if 3*np.sqrt(1/eta) > rmin:
                    print "symfunction_short %s 2 %s %.4f 0.000 %.3f" %(fel, sel, eta, cutoff)
            for i in xrange(len(shift_array)-1):
                if optcutoff: #!TODO To test
                    cutoff = min(np.sqrt(6.908/eta) + shift_array[i], 12.0) ### Variable cutoff depending on the gaussian
                eta = 1./((shift_array[N-i] - shift_array[N-i-1])**2)
                if shift_array[i] + 3*np.sqrt(1/eta) > rmin:
                    print "symfunction_short %s 2 %s %.4f %.3f %.3f" %(fel, sel, eta, shift_array[N-i], cutoff)

    eta_array = np.logspace(-3,0,N/2.)
    zeta_array = [1.000, 4.000, 16.000]
    for fel in elements:
        ang_elements = list(elements)
        for sel in elements:
            for tel in ang_elements:
                print "# symfunctions for type %s 3 %s %s" %(fel, sel, tel)
                for eta in eta_array:
                    #cutoff = min(np.sqrt(3.454/eta),12.000) ### Variable cutoff depending on the gaussian
                    for zeta in zeta_array:
                        if 3*np.sqrt(1/eta) > rmin:
                            print "symfunction_short %s 3 %s %s %.4f  1.000 %.3f %.3f" %(fel, sel, tel, eta, zeta, cutoff)
                            print "symfunction_short %s 3 %s %s %.4f -1.000 %.3f %.3f" %(fel, sel, tel, eta, zeta, cutoff)
            ang_elements.pop(0)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=None)
    parser.add_argument('-e', '--element', type = str, default = 'H', help = 'The elements that will be used to generate the symmetry functions, separated by commas. Any number of elements is accepted. Default is H.')
    parser.add_argument('-c', '--cutoff', type = float, default = 12.0, help = 'The desired cutoff for the symmetry functions. Default is 12.0.')
    parser.add_argument('--optcutoff', action = 'store_true', help = 'If called it optimizes the cutoff to reduce computational load. Reduces (slightly) the accuracy of the fit')
    parser.add_argument('-n', '--ntot', type = int, default = 20, help = 'The number of intervals in which the space is divided. It impacts how many symmetry functions will be generated. Higher numbers indicate that the grid will be thicker and more precise. Default is 20 ')
    parser.add_argument('-r', '--rmin', type = float, default = 0.0, help = 'Distance in [bohr] to the first nearest neighbor. Eliminates the symmetry functions that investigate the space between 0 and rmin. Default is 0.0')
    args = parser.parse_args()
    elements = args.element.split(",")
    sys.exit(main(elements, args.cutoff, args.ntot, args.rmin, args.optcutoff,))
