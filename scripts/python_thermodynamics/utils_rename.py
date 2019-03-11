#!/usr/bin/env python

from __future__ import division
import os
import numpy as np
import sys
import math
import time
import datetime
import inspect

def timeend(measure, ts, funcname):
    '''
        e.g.:
        t1 = datetime.datetime.now()
        self.get_or_load_parameters()
        utils.timeend(measure, t1, "get_or_load_parameters")

    '''
    if measure != True:
        return
    if measure == True:
        t2 = datetime.datetime.now()
        #print datetime.timedelta
        toutms = round(float(str(t2-ts).split("0:00:0")[1])*1000,2)
        print toutms,"ms",printblue(funcname)

def measure_time(process, measure = False):
    if measure == True:
        #ts0 = time.clock()
        #ts1 = time.time()
        ts2 = datetime.datetime.now()
    process
    if measure == True:
        funcname = inspect.stack()[1][4][0].split("(")[1]
        #print round((time.time()-t1)*1000,2), "milli seconds wall time",printblue(funcname)
        print datetime.datetime.now()-ts2, "datetime time",printblue(funcname)
        #print "  ",time.clock() - ts0, "seconds process time"
        #print "  ",time.time() - ts1, "seconds wall time"
    return

def splinefit(xlist,ylist, xstart, xstop, xstep):
    from scipy.interpolate import interp1d
    f = interp1d(xlist,ylist)
    f2 = interp1d(xlist,ylist, kind='cubic')
    xnew = np.arange(xstart,xstop,xstep)

def project_vector(vec, on_vec):
    """  returnes the projected vector;
    vec: current vector
    on_vec: vec is projected on on_vec"""
    verbose = False
    if verbose:
        print ""
    u = vec
    s = on_vec
    if type(u) != np.ndarray:
        sys.exit("u has to be numpy array")
    if type(s) != np.ndarray:
        sys.exit("s has to be numpy array")
    for i in np.arange(u.size):
        u[i] = float(u[i])
    for i in np.arange(s.size):
        s[i] = float(s[i])
    if verbose:
        print "u            :",u
        print "s            :",s
        print "np.dot(u,s)  :",np.dot(u,s)
        print "np.dot(s,s)  :",np.dot(s,s)
    tol = 1e-13
    #if abs(u[0]) <= tol and abs(u[1]) <= tol and abs(u[2]) <= tol and abs(s[0]) <= tol and abs(s[1]) <= tol and abs(s[2]) <= tol:
    #    return 0.0
    if abs(u[0]) <= tol and abs(u[1]) <= tol and abs(u[2]) <= tol:
        return np.array([0.0, 0.0, 0.0])
    if abs(s[0]) <= tol and abs(s[1]) <= tol and abs(s[2]) <= tol:
        return np.array([0.0, 0.0, 0.0])

    #import warnings
    #warnings.filterwarnings('error')
    #try:
    #    np.dot((np.dot(u,s)/abs(np.dot(s,s))),s)
    #except RuntimeWarning:
    #    print "u            :",u,u[0],u[1],u[2]
    #    print "s            :",s
    #    print "np.dot(u,s)  :",np.dot(u,s)
    #    print "np.dot(s,s)  :",np.dot(s,s)
    #    if u[0] == 0.0: # and u[1] == 0.0 and u[2] == 0.0 and s[0] == 0.0 and s[1] == 0.0 and s[2] == 0.0:
    #        print "in1:"
    #    if abs(u[1]) <= tol: # and u[1] == 0.0 and u[2] == 0.0 and s[0] == 0.0 and s[1] == 0.0 and s[2] == 0.0:
    #        print "in2:"
    #    if u[2] == 0.0: # and u[1] == 0.0 and u[2] == 0.0 and s[0] == 0.0 and s[1] == 0.0 and s[2] == 0.0:
    #        print "in3:"
    #    if u[0] == 0.0 and u[1] == 0.0 and u[2] == 0.0 and s[0] == 0.0 and s[1] == 0.0 and s[2] == 0.0:
    #        print "in:"
    #    sys.exit("oops")
    return np.dot((np.dot(u,s)/abs(np.dot(s,s))),s)

def reject_vector(vec, on_vec):
    """  have a look on wiki: http://en.wikipedia.org/wiki/Vector_projection#Scalar_projection
    a: vec     (from this one we only have the direction, not the absolute correct length
    a1: on_vec

    returns: a2 (the rejected vector)
    vec: current vector (longitudinal vec; only the direction is of interest, not the length)
    on_vec: vec is projected on on_vec"""
    b = a1 = on_vec
    adir = lvec = vec
    if False:
        print "b             :",b
        print "adir          :",adir
        print "np(linalg.norm(b)*np.linalg.norm(adir):",(np.linalg.norm(b)*np.linalg.norm(adir))
        #print "np.dot(b,adir):",np.dot(b,np.transpose(adir))
        print "np.dot(b,adir):",np.dot(adir,np.transpose(b))
    cosalpha = np.dot(b,adir)/(np.linalg.norm(b)*np.linalg.norm(adir))
    anormforlength = a1/cosalpha
    alength = np.linalg.norm(anormforlength)
    a = vec*alength/np.linalg.norm(vec)
    a2 = a - a1
    a2length = np.linalg.norm(a2)
    if False:
        print "b = a1 = on_vec  :",b
        print "adir = lvec = vec:",vec
        print "cosalpha         :", cosalpha
        print "anormforlength   :",anormforlength
        print "alength          :",alength
        print "a                :",a
        print "a2               :",a2
        print "a2length         :",a2length
        print "alength*cosalpha :",alength*cosalpha
    return a2

def anglevec(v1,v2):
    """ Returns the angle in radians between vectors 'v1' and 'v2'    """
    cosang = np.dot(v1, v2)
    sinang = np.linalg.norm(np.cross(v1, v2))
    return np.arctan2(sinang, cosang)

def R_2vect(vector_orig, vector_fin, fixed_rotation_axis = False):
    """ from: http://svn.gna.org/svn/relax/tags/1.3.4/maths_fns/frame_order_matrix_ops.py
    with slight changes for calculating the angle

    Calculate the rotation matrix required to rotate from one vector to another.

    For the rotation of one vector to another, there are an infinit series of rotation matrices
    possible.  Due to axial symmetry, the rotation axis can be any vector lying in the symmetry
    plane (should this be "just" a vector?) between the two vectors.  Hence the axis-angle convention will be used to construct the
    matrix with the rotation axis defined as the cross product of the two vectors.  The rotation
    angle is the arccosine of the dot product of the two unit vectors.

    Given a unit vector parallel to the rotation axis, w = [x, y, z] and the rotation angle a,
    the rotation matrix R is::

              |  1 + (1-cos(a))*(x*x-1)   -z*sin(a)+(1-cos(a))*x*y   y*sin(a)+(1-cos(a))*x*z |
        R  =  |  z*sin(a)+(1-cos(a))*x*y   1 + (1-cos(a))*(y*y-1)   -x*sin(a)+(1-cos(a))*y*z |
              | -y*sin(a)+(1-cos(a))*x*z   x*sin(a)+(1-cos(a))*y*z   1 + (1-cos(a))*(z*z-1)  |


    @param R:           The 3x3 rotation matrix to update.
    @type R:            3x3 numpy array
    @param vector_orig: The unrotated vector defined in the reference frame.
    @type vector_orig:  numpy array, len 3
    @param vector_fin:  The rotated vector defined in the reference frame.
    @type vector_fin:   numpy array, len 3
    """
    verbose = False
    if verbose:
        print "vector_orig:",vector_orig
        print "vector_fin:",vector_fin
    R = np.zeros((3,3))
    R[:] = np.nan
    if verbose:
        print "R:",R

    # Convert the vectors to unit vectors.

    import warnings
    warnings.filterwarnings('error')
    try:
        vector_orig = vector_orig / np.linalg.norm(vector_orig)
    except RuntimeWarning:
        print "vector_orig:",vector_orig
        print "np.linalg.norm(vector_orig):",np.linalg.norm(vector_orig)
        print "R:"
        print R

    vector_orig = vector_orig / np.linalg.norm(vector_orig)
    vector_fin = vector_fin / np.linalg.norm(vector_fin)
    if verbose:
        print "vector_orig:",vector_orig
        print "vector_fin:",vector_fin

    # The rotation axis (normalised).
    axis = np.cross(vector_orig, vector_fin)
    axis_len = np.linalg.norm(axis)
    if axis_len != 0.0:
        axis = axis / axis_len

    if type(fixed_rotation_axis) != bool:
        axis = fixed_rotation_axis / np.linalg.norm(fixed_rotation_axis)


    ## Alias the axis coordinates.
    x = axis[0]
    y = axis[1]
    z = axis[2]

    # The rotation angle.
    #print "vector_orig:",vector_orig
    #print "vector_fin:",vector_fin
    #print "np.dot(vector_orig, vector_fin):",np.dot(vector_orig, vector_fin)
    #angle = math.acos(np.dot(vector_orig, vector_fin))
    angle =  anglevec(vector_orig, vector_fin)
    if verbose:
        print "angle:",angle

    ## Trig functions (only need to do this maths once!).
    ca = math.cos(angle)
    sa = math.sin(angle)

    # Calculate the rotation matrix elements.
    R[0,0] = 1.0 + (1.0 - ca)*(x**2 - 1.0)
    R[0,1] = -z*sa + (1.0 - ca)*x*y
    R[0,2] = y*sa + (1.0 - ca)*x*z
    R[1,0] = z*sa+(1.0 - ca)*x*y
    R[1,1] = 1.0 + (1.0 - ca)*(y**2 - 1.0)
    R[1,2] = -x*sa+(1.0 - ca)*y*z
    R[2,0] = -y*sa+(1.0 - ca)*x*z
    R[2,1] = x*sa+(1.0 - ca)*y*z
    R[2,2] = 1.0 + (1.0 - ca)*(z**2 - 1.0)
    return R

def anyvec_to_inplane_senkr_parallel(longvec,vec0,vecs):
    '''
        get longvec and return its inplane perpendicular and parallel parts
        longvec: is the current vector between two atoms (the corresponding vector of the undisplaces structure is vec0)
        vec0:       is the original vecotor of the undisplaced structure
        vecs:       is the vertor which is perpendicular to vec0; only the direction matters, not the magnitude

        - This function will not work well when an atom leavs its Wigner Seitz cell
        - the

        Examples:
        ---------
        vec0:   array([ 2.065, -2.065,  0.   ])
        vecs:   np.array([2.0,2.0,0.0])
                np.array([1.0,1.0,0.0])
                np.array([-2.0,-2.0,0.0])       # all equivalent

        [209]utils.project_vector(utils.reject_vector(np.array([(-2.065+1.271),  2.065,  0.0   ]), vec0),np.array([2.0,2.0,0.0]))
        Out[209]: array([ 0.918019,  0.918019,  0.      ])

        [210]utils.project_vector(utils.reject_vector(np.array([(-2.065-1.271),  2.065,  0.0   ]), vec0),np.array([2.0,2.0,0.0]))
        Out[210]: array([-0.48595, -0.48595,  0.     ])


        utils.project_vector(utils.reject_vector(np.array([(-2.065+0.176),  2.065,  0.0   ]), vec0),np.array([2.0,2.0,0.0]))
        utils.project_vector(utils.reject_vector(np.array([(-2.065+0.176),  2.065,  0.0   ]), vec0),np.array([1.0,1.0,0.0]))
        utils.project_vector(utils.reject_vector(np.array([(-2.065+0.176),  2.065,  0.2   ]), vec0),np.array([-2.0,-2.0,0.0]))
        --> array([ 0.091917,  0.091917,  0.      ])
        --> array([ 0.091917,  0.091917,  0.      ])
        --> array([ 0.091917,  0.091917,  0.      ])


        '''
    dot = np.dot(np.array(vec0),np.array(vecs))
    if dot != 0.0:
        print "vec0:",vec0
        print "vecs:",vecs
        print "dot:",dot
        sys.exit("vectors vec0 and vecs are not perpedicular")

    perpendicular = project_vector(reject_vector(longvec, vec0),vecs)
    projection = project_vector(longvec,vec0)
    parallel = projection + vec0
    return perpendicular,parallel

def coord_transform_to_quer(vec, addwinkel = 0.0):
    ''' in the end rotates vec by + 45 degrees ...
        (np.pi/4 = 0.7853981633974483)
            map vec from a reference system in x,y,z to a refernece system wich has
            x direction in (1/1/0) and y direction in (-1/1/0)
            utils.coord_transform_to_quer([-2,1.0,0]) --> array([-0.70710678118655,  2.12132034355964,  0.              ])
            utils.coord_transform_to_quer([1,1.0,0])  --> array([ 1.41421356237309,  0.              ,  0.              ])
            (see in Dropbox/Screenshots/coordinate_transformation_1.png (or other ending))
    '''
    x = vec[0]
    y = vec[1]
    if abs(vec[2]) > 1e-9:
        #return np.array([0.0, 0.0, 0.0])
        sys.exit("vec[2] has to be 0.0")

    xout = x*np.cos(np.pi/4.+addwinkel)+y*np.sin(np.pi/4.+addwinkel)
    yout = -x*np.sin(np.pi/4.+addwinkel)+y*np.cos(np.pi/4.+addwinkel)
    return np.array([xout,yout,0.0])
    #toxrest = coord_transform_to_quer(toxrest)

def coord_transform_to_par(vec):
    ''' in the end rotate vec by -45 degrees
    utils.coord_transform_to_par([1.0,1.0,0]) --> array([ 0.              ,  1.41421356237309,  0.              ])
    '''
    x = vec[0]
    y = vec[1]
    if abs(vec[2]) != 0.0:
        sys.exit("vec[2] != 0.0")
    xout = x*np.cos(-np.pi/4.)+y*np.sin(-np.pi/4.)
    yout = -x*np.sin(-np.pi/4.)+y*np.cos(-np.pi/4.)
    return np.array([xout,yout,0.0])
    #toxrest = coord_transform_to_quer(toxrest)


def get_index_of_numpyarray(array_all_data=None,sublist=None):
    ''' this is merely a reminder; array indexof are lists'''
    out = np.array([])
    for t in sublist:
        #print "t:",t,"len:",len(np.nonzero( t == tempall)),"||",np.nonzero( t == tempall)[0]
        #print "t:",t,"len:",len(np.nonzero( t == tempall)[0]),"||",np.nonzero( t == tempall)[0]
        if len(np.nonzero( t == array_all_data)[0]) == 0:
            continue
        #print np.nonzero( t == tempall)[0]
        out = np.append(out, np.nonzero( t == array_all_data)[0][0])
    return out
    #return [ np.nonzero( t == array_all_data)[0] for t in sublist]

def remove_duplicates_of_2d_array(inarray = False):
    ''' example:
    inarray =
    array([[1, 1],
           [2, 3],
           [1, 1],
           [5, 4],
           [2, 3]])
    results in:
    array([[1, 1],
           [2, 3],
           [5, 4]])
    '''
    a = inarray[inarray[:,1].argsort()]

    diff = np.diff(a, axis=0)
    ui = np.ones(len(a), 'bool')
    ui[1:] = (diff != 0).any(axis=1)
    #out = a[ui]
    return a[ui]

def remove_duplicates_of_2d_array_within_tolerance(inarray = False, deltaxmax = 1e-6, deltaymax = 1e-4):
    # first three lines sort everything by second column
    # for forces from vasp: deltaxmax = 1e-6, deltaymax = 1e-4
    # it would be nice to  have the possibility to average the canceled values

    a = inarray[inarray[:,1].argsort()]
    #print a
    diff = np.diff(a, axis=0)
    #print "inarray.shape",a.shape
    #print "diff.shape",diff.shape
    xidx = np.nonzero(np.fabs(diff)[:,0] <= np.abs(deltaxmax))[0]
    yidx = np.nonzero(np.fabs(diff)[:,1] <= np.abs(deltaymax))[0]
    remove =  intersect(xidx,yidx)
    #print remove
    #print "------",len(remove)
    new = np.delete(a,remove,0)
    return new[new[:,0].argsort()]

def unique_to_decimals(liste,decimals):
    ''' takes a nupy array (liste) and print uniqe elements but cuts away decials
    e.g. unique_to_decimals(liste,6)
    gets uniqe elements wich are unique to 6th digit
    '''
    d = liste
    NNdist = np.unique(np.sort(np.around(d,decimals = decimals))) #NEVER accuracy is necessary
    NNdist4 = np.unique(np.sort(np.around(d,decimals = 4))) #NEVER accuracy is necessary
    decimalstake = 4
    for kk in np.arange(4,18):
        NNdistcheck = np.unique(np.sort(np.around(d,decimals = kk))) #NEVER accuracy is necessary
        #print kk,len(NNdistcheck),NNdistcheck
        if len(NNdistcheck) == len(NNdist4):
            NNdist = NNdistcheck
            decimalstake = kk
    if decimalstake < decimals:
        sys.exit("not enough digits to get precise next nearest neighbor distance")
    #print "decimalstake:",decimalstake
    return NNdist


def append_row_to_2d_array(inarray = False, addrow = False):
    ''' example:
        dataarray = False
        for i in loop:
            addrow = [lambd, stepsoutnr, stdoutnr, speedup, dudloutnr]
            dataarray = utils.append_row_to_2d_array(inarray = dataarray, addrow=addrow)
    '''
    #print "--> inarray:",inarray
    #print "-->  addrow:",addrow
    if type(inarray) == bool: # or type(inarray) == NoneType:
        #inarray = np.array([[addrow[0], addrow[1]]])
        inarray = np.array([addrow])
    else:
        #newrow = [addrow[0],addrow[1]]
        try:
            inarray = np.vstack([inarray, addrow])
        except ValueError:
            print "inarray:",inarray
            print printred("inarray.shape",inarray.shape)
            print "----------------------"
            print "addrow:",addrow,type(addrow)
            #print printred(" addrow.shape",addrow.shape)
            inarray = np.vstack([inarray, addrow])
    #print ">>> out1   :",inarray,inarray.shape,len(inarray.shape)
    if len(inarray.shape) == 3:
        o1 = inarray.shape[0]
        o2 = inarray.shape[1]
        o3 = inarray.shape[2]
        o12  = o1*o2
        o123 = o1*o2*o3
        #print printred("changing shape:")
        inarray = np.reshape(inarray,(o12,o3))
    #print ">>> out2   :",inarray,inarray.shape,len(inarray.shape)
    return inarray

def write_inputdata(filename = None, text = None):
    ''' append text to filename '''
    if filename == None:
        print "No output written since no filename"
        return
    if text == None:
        text = ""
    if os.path.isfile != "True":
        with open(filename, "a") as myfile:
                myfile.write(text+"\n")
    return

def sort_list_to_sortlist(liste = False, sortlist = False):
    ''' e.g.
    liste = [ 'hallo3nnd', 'halloxdir', 'hallomidd', 'halloquer' ]
    sortlist = [ 'quer' , 'xdir']
    --> out: [ 'halloquer', 'halloxdir', 'hallo3nnd', 'hallomidd' ]
    '''

    outlist = []
    for sort in sortlist:
        for folderidx,folder in enumerate(liste):
            #print folderidx,folder
            if sort in folder:
                outlist.append(folder)
    for folder in liste:
        if folder not in outlist:
            outlist.append(folder)
    if len(liste) != len(outlist):
        for i in liste:
            print i
        print ''
        for i in outlist:
            print i
        sys.exit('len liste is not len outlist')
    return outlist

def get_index_of_list(listall=None, elementslist=None):
    ''' get the indizes of elementslist in list '''
    [listall.index(a) for a in listall]

getVar = lambda searchList, ind: [searchList[i] for i in ind]

def sed(filename,stringsearch,stringreplace):
    import re
    with open(filename, "r") as sources:
        lines = sources.readlines()
    with open(filename, "w") as sources:
        for line in lines:
            sources.write(re.sub(stringsearch, stringreplace, line))
    return

def is_int(x):
    try:
        a = float(x)
        b = int(a)
    except ValueError:
        return False
    else:
        return a == b

def is_float(s):   # "88.0" -> False  ; "88" -> True
    if is_int(s) == True:
        return False
    else:
        try:
            float(s)
            return True
        except ValueError:
            return False

def string_to_list(string):
    ''' gets: 0.224830_1.783284_2.814284 and returnes [ 0.225  1.783  2.814] (with all the digits)
        gets: '1,2,3' and returnes [1, 2, 3]'''
    #print "string:",string
    if type(string) == str:
        #print "string:",string,type(string)
        string = string.replace(",", " ")
        #print "string:",string
        stringout = string.split()
        #print "stringout:",stringout
        if len(stringout) == 1:
            # in case it is seperated with "_"
            stringout1 =  [ float(i) for i in string.split("_") ]
            #print "stringout1:",stringout1
            return stringout1
        else:
            stringout2 = [ float(i) for i in stringout ]
            #print "stringout2:",stringout2
            return stringout2
    if type(string) == bool:
        return string
    if type(string) == list:
        return string
    if type(string) == np.ndarray:
        return list(string)

def list_to_string(listin, digits=14):
    if type(listin) == bool:
        return listin
    stringout = ""
    for i in listin:  # go through every number
        #print "i:",i
        accuracy = '%.'+str(digits)+'f'
        if abs(i) <= math.pow(10,-digits):
            accuracy = '%.1f'

        stringadd = str(accuracy % i)[:digits+1]
        if stringadd == '-0.0': stringadd = '0.0'
        a = [x for x in stringadd if x != '-' ]
        b = [x for x in a if x != '.' ]
        #print "j:",stringadd,list(stringadd),"||",b,">>",set(b),len(set(b)),list(set(b))[0]
        if len(set(b)) == 1:
            #print "yes1:"
            if list(set(b))[0] == '0':
                #print "yes2:"
                stringadd = '0.0'


        if stringout == "":
            stringout = stringout+stringadd
        else:
            stringout = stringout+"_"+stringadd
        #print ">",stringadd
        #print ""
    #print "stringout:",stringout
    return stringout

def string_to_num_list(string, tostring = False):
    """ Turn a string into a list of int/float chunks.
        "4.Ang_300K_.2m-9m=.77k" -> [4.0, 300.0, 0.2, 9.0, 0.77] """

    def tryfloat(stringtest):
        try:
            return float(stringtest)
        except:
            return stringtest
    import re
    out = [ tryfloat(c) for c in re.findall(r"[-+]?\d*\.\d+|\d+", string)]
    if tostring:
        return [ str(i) for i in out]
    return out

def number_to_string_of_certain_length(number, digits, chars):
    ''' number: the number to print
        digits: the number of digits to show e.g. 2: 0.34
        chars: the amount of characters to use for this output
        '''
    stringout = str(round(number,digits))
    if digits == 0:
        stringout = str(int(number))
    stringoutadd = chars - len(stringout)  # keep 4
    aaa = ' '*stringoutadd+stringout
    return aaa

def string_add_spaces(string,chars,left = True):
    charsadd = chars - len(str(string))
    if left == True:
        return " "*charsadd+str(string)
    else:
        return str(string)+" "*charsadd

def file_len(fname):
    with open(fname) as f:
        for i, l in enumerate(f):
            pass
    return i + 1

def lsn(searchstring):
    """ Returns files sorted numerically
        searchstring: e.g. "*Ang_*K"            """
    if type(searchstring) == str:
        import glob
        files = glob.glob(searchstring)
    elif type(searchstring) == list:
        files = searchstring
    else:
        sys.exit("need string or list for lsn")
    return sorted((files), key=string_to_num_list)

def list_sorted(liste):
    """ Returns list sorted numerically
        searchstring: e.g. "*Ang_*K"            """
    return sorted((liste), key=string_to_num_list)

def string_to_mathematical_expression_does_not_work_properly(x):
    def parse(x):
        operators = set('+-*/')
        op_out = []    #This holds the operators that are found in the string (left to right)
        num_out = []   #this holds the non-operators that are found in the string (left to right)
        buff = []
        for c in x:  #examine 1 character at a time
            print "c:",c
            if c in operators:
                #found an operator.  Everything we've accumulated in `buff` is
                #a single "number". Join it together and put it in `num_out`.
                num_out.append(''.join(buff))
                buff = []
                op_out.append(c)
            else:
                #not an operator.  Just accumulate this character in buff.
                buff.append(c)
                #buff.append(float(c))
        num_out.append(''.join(buff))
        #num_out = filter(None, num_out)
        #op_out = filter(None, op_out)
        return num_out,op_out

    import operator
    def my_eval(nums,ops):

        nums = list(nums)
        ops = list(ops)
        operator_order = ('*/','+-')  #precedence from left to right.  operators at same index have same precendece.
                                      #map operators to functions.
        op_dict = {'*':operator.mul,
                   '/':operator.div,
                   '+':operator.add,
                   '-':operator.sub}
        Value = None
        for op in operator_order:                   #Loop over precedence levels
            while any(o in ops for o in op):        #Operator with this precedence level exists
                idx,oo = next((i,o) for i,o in enumerate(ops) if o in op) #Next operator with this precedence
                ops.pop(idx)                        #remove this operator from the operator list
                print "ind:",ind,"idx:",idx #"nnn:",nums[ind:idx+2]
                values = map(float,nums[idx:idx+2]) #here I just assume float for everything
                value = op_dict[oo](*values)
                nums[idx:idx+2] = [value]           #clear out those indices

        return nums[0]
    print "parse(x):",parse(x)
    return my_eval(*parse(x))



### hesse.py al -cif cartesian_coords DOES NOT WORK ON CMMD002 if following 16 lines are commented in
#from pyparsing import (Literal,CaselessLiteral,Word,Combine,Group,Optional,
#                       ZeroOrMore,Forward,nums,alphas,oneOf)
#import math
#import operator
#
#__author__='Paul McGuire'
#__version__ = '$Revision: 0.0 $'
#__date__ = '$Date: 2009-03-20 $'
#__source__='''http://pyparsing.wikispaces.com/file/view/fourFn.py
#http://pyparsing.wikispaces.com/message/view/home/15549426
#'''
#__note__='''
#All I've done is rewrap Paul McGuire's fourFn.py as a class, so I can use it
#more easily in other places.
#'''

class NumericStringParser(object):
    '''
    Most of this code comes from the fourFn.py pyparsing example

    '''
    def pushFirst(self, strg, loc, toks ):
        self.exprStack.append( toks[0] )
    def pushUMinus(self, strg, loc, toks ):
        if toks and toks[0]=='-':
            self.exprStack.append( 'unary -' )
    def __init__(self):
        """
        expop   :: '^'
        multop  :: '*' | '/'
        addop   :: '+' | '-'
        integer :: ['+' | '-'] '0'..'9'+
        atom    :: PI | E | real | fn '(' expr ')' | '(' expr ')'
        factor  :: atom [ expop factor ]*
        term    :: factor [ multop factor ]*
        expr    :: term [ addop term ]*
        """
        point = Literal( "." )
        e     = CaselessLiteral( "E" )
        fnumber = Combine( Word( "+-"+nums, nums ) +
                           Optional( point + Optional( Word( nums ) ) ) +
                           Optional( e + Word( "+-"+nums, nums ) ) )
        ident = Word(alphas, alphas+nums+"_$")
        plus  = Literal( "+" )
        minus = Literal( "-" )
        mult  = Literal( "*" )
        div   = Literal( "/" )
        lpar  = Literal( "(" ).suppress()
        rpar  = Literal( ")" ).suppress()
        addop  = plus | minus
        multop = mult | div
        expop = Literal( "^" )
        pi    = CaselessLiteral( "PI" )
        expr = Forward()
        atom = ((Optional(oneOf("- +")) +
                 (pi|e|fnumber|ident+lpar+expr+rpar).setParseAction(self.pushFirst))
                | Optional(oneOf("- +")) + Group(lpar+expr+rpar)
                ).setParseAction(self.pushUMinus)
        # by defining exponentiation as "atom [ ^ factor ]..." instead of
        # "atom [ ^ atom ]...", we get right-to-left exponents, instead of left-to-right
        # that is, 2^3^2 = 2^(3^2), not (2^3)^2.
        factor = Forward()
        factor << atom + ZeroOrMore( ( expop + factor ).setParseAction( self.pushFirst ) )
        term = factor + ZeroOrMore( ( multop + factor ).setParseAction( self.pushFirst ) )
        expr << term + ZeroOrMore( ( addop + term ).setParseAction( self.pushFirst ) )
        # addop_term = ( addop + term ).setParseAction( self.pushFirst )
        # general_term = term + ZeroOrMore( addop_term ) | OneOrMore( addop_term)
        # expr <<  general_term
        self.bnf = expr
        # map operator symbols to corresponding arithmetic operations
        epsilon = 1e-12
        self.opn = { "+" : operator.add,
                "-" : operator.sub,
                "*" : operator.mul,
                "/" : operator.truediv,
                "^" : operator.pow }
        self.fn  = { "sin" : math.sin,
                "cos" : math.cos,
                "tan" : math.tan,
                "abs" : abs,
                "trunc" : lambda a: int(a),
                "round" : round,
                "sgn" : lambda a: abs(a)>epsilon and cmp(a,0) or 0}
    def evaluateStack(self, s ):
        op = s.pop()
        if op == 'unary -':
            return -self.evaluateStack( s )
        if op in "+-*/^":
            op2 = self.evaluateStack( s )
            op1 = self.evaluateStack( s )
            return self.opn[op]( op1, op2 )
        elif op == "PI":
            return math.pi # 3.1415926535
        elif op == "E":
            return math.e  # 2.718281828
        elif op in self.fn:
            return self.fn[op]( self.evaluateStack( s ) )
        elif op[0].isalpha():
            return 0
        else:
            return float( op )
    def eval(self,num_string,parseAll=True):
        self.exprStack=[]
        results=self.bnf.parseString(num_string,parseAll)
        val=self.evaluateStack( self.exprStack[:] )
        return val


def string_to_mathematical_expression(x):
    nsp=NumericStringParser()
    return nsp.eval(x)


def common_prefix(strings, common_suffix = False):
    """ Find the longest string that is a prefix (or suffix) of all the strings of the input list.
    """
    if not strings:
        return ''

    if common_suffix == True:
        strings = [i[::-1] for i in strings]

    prefix = strings[0]
    for s in strings:
        if len(s) < len(prefix):
            prefix = prefix[:len(s)]
        if not prefix:
            return ''
        for i in range(len(prefix)):
            if prefix[i] != s[i]:
                prefix = prefix[:i]
                break
    if common_suffix == True:
        return prefix[::-1]
    else:
        return prefix


def schnittmenge(list1, list2):
    """
    nur elemente die sowolh in list1 als auch in list2 sine
    """
    ## end of http://code.activestate.com/recipes/576694/ }}}
    #import collections
    #class OrderedSetn(collections.OrderedDict, collections.MutableSet):
    #    def update(self, *args, **kwargs):
    #        if kwargs:
    #            raise TypeError("update() takes no keyword arguments")
    #        for s in args:
    #            for e in s:
    #                self.add(e)
    #    def add(self, elem):
    #        self[elem] = None
    #    def discard(self, elem):
    #        self.pop(elem, None)
    #    def __le__(self, other):
    #        return all(e in other for e in self)
    #    def __lt__(self, other):
    #        return self <= other and self != other
    #    def __ge__(self, other):
    #        return all(e in self for e in other)
    #    def __gt__(self, other):
    #        return self >= other and self != other
    #    def __repr__(self):
    #        return 'OrderedSet([%s])' % (', '.join(map(repr, self.keys())))
    #    def __str__(self):
    #        return '{%s}' % (', '.join(map(repr, self.keys())))
    #    difference = property(lambda self: self.__sub__)
    #    difference_update = property(lambda self: self.__isub__)
    #    intersection = property(lambda self: self.__and__)
    #    intersection_update = property(lambda self: self.__iand__)
    #    issubset = property(lambda self: self.__le__)
    #    issuperset = property(lambda self: self.__ge__)
    #    symmetric_difference = property(lambda self: self.__xor__)
    #    symmetric_difference_update = property(lambda self: self.__ixor__)
    #    union = property(lambda self: self.__or__)

    #import myclasses
    #sl1 = myclasses.OrderedSetn(list1)
    #sl2 = myclasses.OrderedSetn(list2)
    #sl1 = OrderedSetn(list1)
    #sl2 = OrderedSetn(list2)
    #    print sl1
    #    print sl2
    return [val for val in list1 if val in list2]
    #return list(sl1 & sl2)

def schnittmenge_several_lists(lists, verbose = False):
    """
    nur elemente die in allen listen sind
    """
    list1 = lists[0]
    schnittmenge = lists[0]
    for ind, i in enumerate(lists):
        #print "i :",i
        l1 = schnittmenge
        l2 = i
        if verbose:
            print "i:",ind,"l1:",l1
            print "i:",ind,"l2:",l2
        schnittmenge = [val for val in l1 if val in l2]
        if verbose:
            print "i:",ind,"schn:",schnittmenge
            print "------------"*10
        #schnittmenge = [val for val in li] if val in schnittmenge]
    return schnittmenge

""" NOTES:
      - requires Python 2.4 or greater
      - elements of the lists must be hashable
      - order of the original lists is not preserved
"""
def unique(a):
    """ return the list with duplicate elements removed """
    return list(set(a))

def intersect(a, b):
    """ return the intersection of two lists """
    return list(set(a) & set(b))

def union(a, b):
    """ return the union of two lists """
    return list(set(list(a)) | set(list(b)))

#if __name__ == "__main__":
#    a = [0,1,2,0,1,2,3,4,5,6,7,8,9]
#    b = [5,6,7,8,9,10,11,12,13,14]

#    print unique(a)
#    print intersect(a, b)
#    print union(a, b)

#       Results:
#
#       [0, 1, 2, 3, 4, 5, 6, 7, 8, 9]
#       [8, 9, 5, 6, 7]
#       [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14]


def get_commom_part_of_filenames(liste = None):
    """
    """
    ka = u.common_prefix(liste)
    import re
    ka = re.sub('\.$','',ka)
    ka = re.sub('[0-9]*$','',ka)
    ka = re.sub('\.$','',ka)
    return ka

def find_closest_sorted(A, target):
    #A must be sorted
    idx = A.searchsorted(target)
    idx = np.clip(idx, 1, len(A)-1)
    left = A[idx-1]
    right = A[idx]
    idx -= target - left < right - target
    return idx

def return_function2_on_grid_of_f1(f1,f2):
    ''' grid is taken from f1; f1, f2 are both 2d arrays;
    the returned f2 is on the xgrid of f1 , ypoints are interpolated
    it might occur that the output has nans! '''
    import scipy
    from scipy.interpolate import interp1d
    xlist = f2[:,0]
    ylist = f2[:,1]
    func2 = interp1d(f2[:,0],f2[:,1],kind='linear', axis=-1, copy=True, bounds_error=False, fill_value=np.nan)
    return np.array([f1[:,0],func2(f1[:,0])]).transpose()
    sys.exit()


    xlistb = functionsout[1][:,0]
    ylistb = functionsout[1][:,1]
    #f2 = interp1d(xlist,ylist,kind='cubic' , axis=-1, copy=True, bounds_error=False, fill_value=np.nan)

    f1b = interp1d(xlistb,ylistb,kind='linear', axis=-1, copy=True, bounds_error=False, fill_value=np.nan)

    arrayout=np.zeros((len(xliste[2:]),2))
    arrayoutb=np.zeros((len(xliste[2:]),2))
    for idx,x in enumerate(xliste[2:]):
        #print x,f1(x),f2(x)
        arrayout[idx,0] = x
        arrayout[idx,1] = f1(x)
        arrayoutb[idx,0] = x
        arrayoutb[idx,1] = f1b(x)


    #np.savetxt("/Users/glensk/tmppo",np.array([xlist,ylist]).transpose())
    #np.savetxt("/Users/glensk/tmpp1",arrayout)
    #np.savetxt("/Users/glensk/tmpp1b",arrayoutb)
    return arrayout,arrayoutb

def return_functions_on_same_grid(f1,f2,verbose=True):
    ''' f1, f2 are both 2d arrays, the grid created is a new grid neither from f1 nor from f2'''
    mmin = mmax = False
    lenmax = 0
    functions = [ f1,f2]
    functionsout = [ f1,f2]

    # get minimum interval
    for f in functions:
        mmintmp = f[:,0].min()
        mmaxtmp = f[:,0].max()
        lencheck = f[:,0].shape[0]
        if lencheck > lenmax:
            lenmax = lencheck
        if verbose == True:
            print "mmintmp,mmaxtmp:",mmintmp,mmaxtmp
        if type(mmin) == bool:
            mmin = mmintmp
        if type(mmax) == bool:
            mmax = mmaxtmp

        if mmintmp > mmin:
            mmin = mmintmp
        if mmaxtmp < mmax:
            mmax = mmaxtmp
    if verbose == True:
        print "mmin,mmax:",mmin,mmax
        print "lenmax:",lenmax
    #liste = np.arange(len(*functions))
    #print "liste:",liste
    xvalues = np.linspace(mmin,mmax,lenmax)
    from scipy.interpolate import UnivariateSpline
    b = 2.9/np.sqrt(2.)
    s1 = UnivariateSpline(f1[:,0], f1[:,1],k=3,s=0.00000000000*b)
    s2 = UnivariateSpline(f2[:,0], f2[:,1],k=3,s=0.00000000000*b)
    arrayout =np.zeros((len(xvalues),2))
    arrayoutb=np.zeros((len(xvalues),2))
    arrayout[:,0] = arrayoutb[:,0] = xvalues
    arrayout[:,1]  = s1(xvalues)
    arrayoutb[:,1] = s2(xvalues)
    return arrayout,arrayoutb

    #for fidx,f in enumerate(functions):
    #    b1 = np.nonzero((f[:,0] >= mmin))[0]
    #    b2 = np.nonzero(f[:,0] <= mmax)[0]
    #    bintersect = intersect(b1,b2)
    #    functionsout[fidx] = f[bintersect]
    #    if fidx == 0:
    #        xliste = functionsout[fidx][:,0]
    #    else:
    #        if functionsout[fidx].shape[0] < functionsout[0].shape[0]:
    #            xliste = functionsout[fidx][:,0]
    #print xliste
    #print "zzzzzz"
    #print functionsout[0]
    #print ""
    #print functionsout[1]
    #print "zzzzzz"

    #def splinefit(xlist,ylist, xstart, xstop, xstep):
    #    f = interp1d(xlist,ylist)
    #    f2 = interp1d(xlist,ylist, kind='cubic')
    #    xnew = np.arange(xstart,xstop,xstep)


    import scipy
    from scipy.interpolate import interp1d
    xlist = functionsout[0][:,0]
    ylist = functionsout[0][:,1]
    xlistb = functionsout[1][:,0]
    ylistb = functionsout[1][:,1]
    print "xlist:"
    print type(xlist)
    print "y;"
    print type(ylist)
    print type(f1[:,0])
    print type(f1[:,1])
    from scipy.interpolate import UnivariateSpline
    s = UnivariateSpline(f1[:,0], f1[:,1], s=1)
    print "s:",s
    ys = s(xs)


    f3 = interp1d(f1[:,0],f1[:,1],kind='linear', axis=-1, copy=True, bounds_error=False, fill_value=np.nan)


    f1 = interp1d(xlist,ylist,kind='linear', axis=-1, copy=True, bounds_error=False, fill_value=np.nan)
    #f2 = interp1d(xlist,ylist,kind='cubic' , axis=-1, copy=True, bounds_error=False, fill_value=np.nan)
    f3 = interp1d(f1[:,0],f1[:,1],kind='linear', axis=-1, copy=True, bounds_error=False, fill_value=np.nan)

    f1b = interp1d(xlistb,ylistb,kind='linear', axis=-1, copy=True, bounds_error=False, fill_value=np.nan)
    #f2b = interp1d(xlistb,ylistb,kind='cubic' , axis=-1, copy=True, bounds_error=False, fill_value=np.nan)
    f3b = interp1d(f2[:,0],f2[:,1],kind='linear', axis=-1, copy=True, bounds_error=False, fill_value=np.nan)

    #print ">>>>>>"
    #print xliste
    #print ylist
    xliste=xvalues
    for idx,x in enumerate(xliste):
        #print x,f1(x),f2(x)
        arrayout[idx,0] = x
        arrayout[idx,1] = f1(x)
        arrayoutb[idx,0] = x
        arrayoutb[idx,1] = f1b(x)

    #arrayout=np.zeros((len(xliste[2:]),2))
    #arrayoutb=np.zeros((len(xliste[2:]),2))
    #for idx,x in enumerate(xliste[2:]):
    #    #print x,f1(x),f2(x)
    #    arrayout[idx,0] = x
    #    arrayout[idx,1] = f1(x)
    #    arrayoutb[idx,0] = x
    #    arrayoutb[idx,1] = f1b(x)


    #np.savetxt("/Users/glensk/tmppo",np.array([xlist,ylist]).transpose())
    #np.savetxt("/Users/glensk/tmpp1",arrayout)
    #np.savetxt("/Users/glensk/tmpp1b",arrayoutb)
    return arrayout,arrayoutb

def cut_function_at_DOS(function,DOS):
    xmin = DOS[:,0].min()
    xmax = DOS[:,0].max()
    b1 = np.nonzero((function[:,0] >= xmin))[0]
    b2 = np.nonzero((function[:,0] <= xmax))[0]
    bintersect = intersect(b1,b2)
    return function[bintersect]

def getDOS_of_1d_array(data, normalize_to_1 = False):
    #data = np.loadtxt(self.jobvorlage_vecnormlon_file)[2:]
    #data = np.loadtxt("/Users/glensk/Dropbox/Understand_distributions/jobvorlage_all/2x2x2sc_Pt_4.1/tests/tveclonall.dat")
    dist_space = np.linspace( min(data), max(data), 100 )
    fact = 3.0
    def my_kde_bandwidth(obj, fac=fact):
        """We use Scott's Rule, multiplied by a constant factor."""
        return np.power(obj.n, -1./(obj.d+4)) * fac

    from scipy.stats.kde import gaussian_kde
    kde = gaussian_kde( data ,bw_method=my_kde_bandwidth)
    #kde = gaussian_kde( data )
    #np.savetxt(self.jobvorlage+"/DOSlonpy"+str(fact),np.array([dist_space, kde(dist_space)]).transpose())
    #np.savetxt("DOS",np.array([dist_space, kde(dist_space)]).transpose())
    dos = np.array([dist_space, kde(dist_space)]).transpose()
    if normalize_to_1 == True:
        dosymax = dos[:,1].max()
        dosy = dos[:,1]/dosymax
        dos[:,1] = dosy
    return dos

def DOS_cut(longvecnormfile, dos, nndist = False, filenamedos = False):
    ''' longvecnormfile is the file which was used to actually make the dos '''
    # takes time to load in, just load in if necessary
    longvecnorm = np.loadtxt(longvecnormfile)
    lonmin = round(longvecnorm.min()-0.01,2)
    lonmax = round(longvecnorm.max()+0.01,2)
    #print "lonmin/max:",lonmin,lonmax
    #print self.jobvorlage_DOSlon
    #print "smaller----"
    #print "doss:"
    #print dos
    smaller = np.nonzero(dos[:,0]<=lonmax)[0]
    greater = np.nonzero(dos[:,0]>=lonmin)[0]
    intersect = schnittmenge(smaller,greater)

    greaterconservative = np.nonzero(dos[:,1]>=0.003)[0]
    #unionhere = union(greaterconservative,intersect)

    DOSloncut = np.array([dos[:,0][greaterconservative],dos[:,1][greaterconservative]]).transpose()
    DOSloncutshifted = False
    if nndist != False:
        DOSloncutshifted = np.array([DOSloncut[:,0]-nndist,DOSloncut[:,1]]).transpose()
    if filenamedos != False:
        np.savetxt(filenamedos+"cut",DOSloncut)
        np.savetxt(filenamedos+"cutshifted",DOSloncutshifted)
    return DOSloncut,DOSloncutshifted

def MatrixPlot(M = None, temp = 0.05):
    ''' makes a MatrixPlot similar to mathematica '''
    import matplotlib
    import matplotlib.pyplot as plt
    plt.clf()
    plt.imshow(M, interpolation='nearest', cmap='RdBu',vmin=M.min(), vmax=M.max())
    plt.grid(True)
    #plt.imshow(M, interpolation='nearest')
    plt.colorbar()   # makes the bar showing the vales
    print "a colormap would be great ...."
    #plt.clim([-temp,temp])  # this is not good for matrix np.array(([1,2],[3,4]))

def read_only_certain_lines_of_file(filename,start,stop):
    fp = open(filename)
    out = False
    for i, line in enumerate(fp):
        line = line.rstrip()
        if i == start:
            if type(out) == bool:
                out = [ line ]
            # 26th line
        elif i < stop:
            out.extend([line])
            #print line
            # 30th line
        elif i >= stop:
            break
    fp.close()
    return out


def run2(command=None, dont_raise_exceptino = False):
    """
    constantly prints output, not just at the end
    """
    import subprocess
    process = subprocess.Popen(
        command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    output = ''

    # Poll process for new output until finished
    for line in iter(process.stdout.readline, ""):
        if line != "" or line != " ":
            print line,
        output += line
    process.wait()
    exitCode = process.returncode

    if dont_raise_exceptino == True:
        return output

    if (exitCode == 0):
        return output
    else:
        raise Exception(command, exitCode, output)

def run(befehl):
    """
    z.B. run('ls -l')
    z.B. run('frompath_stoich.sh')
    """
    import subprocess
    return subprocess.check_output(befehl, shell=True, stderr=subprocess.STDOUT)

def readfile(path=None, line=None, split = False):
    if path is None:
        print "i was called from: " + str(which_function_calledme())
        exit(error="Please specify path")
    import os
    if os.path.isfile(path) != True:
        import sys
        sys.exit(path+" does not exist!")

    f = open(path, "r")
    try:
        out = f.read()
    finally:
        f.close()
        if line is None:
            returnto = out.splitlines()
        else:
            returnto = out.splitlines()[line - 1]
    if split == True:
        returnto = [ i.split() for i in returnto]
    return returnto

def isfile(filename):
    if os.path.isfile(filename) != True:
        sys.exit(filename+ "  does not exist")
    return filename

def printoutcolor(red,var,ENDC):
    if len(var) == 1:
        return red + str(var[0]) + ENDC
    else:
        return red + str(var) + ENDC

def printnormal(*var):
    ENDC = '\033[0m'
    #return ENDC + str(var) + ENDC
    return printoutcolor(red,var,ENDC)
def printred(*var):
    #print "lenred:",len(var)
    red = '\033[31m'
    ENDC = '\033[0m'
    return printoutcolor(red,var,ENDC)
def printgreen(*var):
    red = '\033[32m'
    ENDC = '\033[0m'
    #return red + str(var) + ENDC
    return printoutcolor(red,var,ENDC)
def printyellow(*var):
    #print "lenvar:",len(var)
    red = '\033[33m'
    ENDC = '\033[0m'
    #return red + str(var) + ENDC
    return printoutcolor(red,var,ENDC)
def printblue(*var):
    red = '\033[34m'
    ENDC = '\033[0m'
    #return red + str(var) + ENDC
    return printoutcolor(red,var,ENDC)
def printokblue(*var):
    red = '\033[94m'
    ENDC = '\033[0m'
    #return red + str(var) + ENDC
    return printoutcolor(red,var,ENDC)
def printpink(*var):
    red = '\033[35m'
    ENDC = '\033[0m'
    #return red + str(var) + ENDC
    return printoutcolor(red,var,ENDC)
def printblackbold(*var):
    red = '\033[38m\033[1m'
    ENDC = '\033[0m'
    #return red + str(var) + ENDC
    return printoutcolor(red,var,ENDC)
def printredbold(*var):
    red = '\033[31m\033[1m'
    ENDC = '\033[0m'
    #return red + str(var) + ENDC
    return printoutcolor(red,var,ENDC)


def printarray(array,color=printnormal,decimals=2,prepend=False):
    def printline(line,decimals=decimals,color=color,prepend=prepend):
        strout = ""
        for idx,i in enumerate(line):
            if type(decimals) == list:
                decimalshier = decimals[idx]
                #print "tt:",type(decimals),"decimalshier:",decimalshier,"i:",i,np.around(i,decimals=decimalshier)
            else:
                decimalshier = decimals
            if decimalshier == 0:
                string = str(int(np.around(i, decimals=decimalshier)))
            else:
                string = str(np.around(i, decimals=decimalshier))
            # convert -0 to proper format
            if string == '-0':
                i = 0.0
                if decimalshier == 0:
                    string = str(int(np.around(i, decimals=decimalshier)))
                else:
                    string = str(np.around(i, decimals=decimalshier))

            if string[0] != "-":
               string = " "+string

            prependhier=""
            if type(prepend) == list:
                #print 'idx:',idx,"len:",len(prepend)
                if len(prepend) > idx:
                    prependhier = str(prepend[idx])
                else:
                    prependhier = ''


            strout = strout+prependhier+string+'\t'
        #if color == "green":
        #print "color:",color
        #col = utils.eval(str(color))
        #print "col:",col
        return color(strout)

    for i in array:
        print printline(i,color=color)


def plot_find_ylim_from_xlim(xyfunc, xliml,xlimr, addmargin = True):
    ''' gets a 2d nupy array with x and y data and finds for given
    xliml und ximr the corresponding ylims '''
    if xyfunc[:,0].min() > xliml:
        xliml = xyfunc[:,0].min()
    if xyfunc[:,0].max() < xlimr:
        xlimr = xyfunc[:,0].max()
    import matplotlib.pyplot as plt
    plt.xlim(xliml,xlimr)
    def find_nearest(array,value):
        idx = (np.abs(array-value)).argmin()
        return idx
    ylimlidx = find_nearest(xyfunc[:,0],xliml)
    ylimridx = find_nearest(xyfunc[:,0],xlimr)
    ylimmin = xyfunc[ylimlidx:ylimridx,1].min()
    ylimmax = xyfunc[ylimlidx:ylimridx,1].max()
    if addmargin == True:
        d = np.linalg.norm(ylimmax - ylimmin)
        #print "d:",d
        part = 1/20.
    ylimmin = ylimmin-d*part
    ylimmax = ylimmax+d*part
    #print "ylimin,ylimmax;",ylimmin,ylimmax
    return xliml,xlimr,ylimmin,ylimmax

def file2dat(filename):
    with open(filename, 'r') as file: return file.readlines()

def dat2file(filename, data):
    with open(filename, 'w') as file: file.writelines(data)
