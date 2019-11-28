#!/usr/bin/python
import get, my, argparse, sys, os, db, murn
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
def from_path_string_displacement(path=None,exit=True):
    """
    -fpsd
    """
    if path == None:
        path = my.pwd()
    string = path.split("/")
    if len(string) <= 9:
        if exit == False:
            return None
        else:
            my.exit("Not deep enoug to bi in displacement")
    check = string[9]
    check2 = my.grep("[0-9]*_[0-9]*",check)
    if len(check2) != 1:
        my.exit(error="displacement is wierd")
    if check != check2[0]:
        my.exit(error="displacement is wierd2")
    return check

def from_path_string_displacedatom(path=None,exit=True):
    """
    -fpsda
    """
    if path == None:
        path = my.pwd()
    return int(from_path_string_displacement(path=path).split("_")[0])

def from_path_string_displacedxyzcoord(path=None,exit=True):
    """
    -fpsdxyz
    """
    if path == None:
        path = my.pwd()
    return int(from_path_string_displacement(path=path).split("_")[1])



def get_path_displacement(path=None,Error_if_notexist=True):
    """
    -gpd
    """
    if path == None:
        path = my.pwd()
    from_path_string_displacement(path=path,exit=Error_if_notexist)
    out = "/".join(path.split("/")[0:10])
    return out

def HesseMatrix_to_ExactFreqs_wo_zeroFreqs(infile="HesseMatrix",outfile=None,output_as_nparray=False):
    """
    -chef
    """
    ################################
    ### check input file
    ################################
    my.checkfile(infile)
    infile = os.path.abspath(infile)

    debug = "yes"
    debug = "no"
    #print "infile:",infile
    #print "outfile",outfile

    ## You have: hbar(hartree/(bohrradius^2 u))^(1/2)
    ## You want: meV
    ##    * 637.33912
    dynMatToFreq = 637.339113585464
    infile = str(infile)

    ##########################################
    ### import and extract eigenvalues
    ##########################################
    import numpy as np
    from numpy import linalg
    h = np.loadtxt(infile)
    #### check wether we have just a unary or an alloy (imortant for masses)
    elements = get.element_real(path = infile)
    elements_amount = len(elements.split(" "))
    #print "elements:", elements,len(elements.split(" "))
    if debug == "yes":
        print "---h---"
        print h
        print "---h---",elements_amount
    if elements_amount == 1:
        #print "iii"
        mass = db.AtomicWeight()
        mass = 26.982
        d = h/mass
        #print "innn"
    else:
        #print "infile",infile
        file_w = infile+"_weights_all"
        #print "ww:", file_w
        my.checkfile(file_w)
        massesall = np.loadtxt(file_w)
        d = h/np.sqrt(massesall)
    if debug == "yes":
        print "---d---"
        print d
        print "---d---"

    #quit()
#    d_sym = ((np.conjugate(d) - d)/2 + d)/mass
#    print "dc:",dc
#    print d
#    print ""
    #print "yes"
    eigenvalues,   eigenvectors = linalg.eig(d)
    eigenvalues_realpart = np.sort(eigenvalues.real)
#    print eigenvalues_realpart


#    eigenvalues_sym, eigenvectors_sym = linalg.eig(d_sym)
#    eigenvalues_realpart_sym = np.sort(eigenvalues_sym.real)

    if debug == "yes":
        print eigenvalues
        print "--------------"
        print eigenvalues_realpart
        print "--------------"

    #### Chop vales smaller 1e-8 (this is little)
    tol=1e-8
    from numpy import ma
    eigenvalues_out = ma.masked_inside(eigenvalues_realpart,-tol, tol)
    ev =  eigenvalues_out.filled(0)


    if debug == "yes":
        print "############################### ev ##################"
        print ev

    #ev[0] = 0
    #ev[1] = 0
    #ev[2] = 0

    ### check if we have negative parts
    if ev[0] < 0:
        print "Negative Eigenvalues in "+infile
        print ev[0]
        my.exit("ERROR: Negative Eigenvalues in "+infile)

    ### some conversions
    a = np.sqrt(ev)*dynMatToFreq
    if debug == "yes":
        print "+++++++++++++++++++++++++++++++++++++++ a "
        print a

    ### check if first 3 eigenvalues are 0
    if a[0] != 0 : my.exit("a[0] "+str(a[0])+" is not 0")
    if a[1] != 0 : my.exit("a[1] "+str(a[1])+" is not 0")
    if a[2] != 0 : my.exit("a[2] "+str(a[2])+" is not 0")

    out = np.nan_to_num(a[3:])
    out = out[::-1]

    if output_as_nparray == True: return out

    ## from here: print to screen or print to file
    if outfile != None:
#        from numpy import savetxt
        np.savetxt(str(outfile),out, fmt='%.10f')
        return str(outfile)+" written"
    else:
#        print list(out)
        for i in out:
            print i
    return " "

def HesseMatrix_to_MeanFreqs(infile="HesseMatrix",outfile=None,output_as_nparray=True):
    """
    -htmf
    """
    my.checkfile(infile)
    infile = os.path.abspath(infile)
    exact_freqs = HesseMatrix_to_ExactFreqs_wo_zeroFreqs(infile=infile,outfile=outfile,output_as_nparray=output_as_nparray)
    import numpy as np
    return np.mean(exact_freqs)

def HesseMatrix_to_Fqh_fromExactFreqs(infile="HesseMatrix",Tmax=None,outfile=None,forSupercell=None,output_as_nparray=False):
    """
    -htfqhe gibes Fqh per atom in [eV]
    """
    debug="no"
    my.checkfile(infile)
    infile = os.path.abspath(infile)

    if debug == "yes":
        print "infile: ",infile
        print "forSupercell: ",forSupercell
    if forSupercell == True:
        numatoms = get.from_path_get_numatoms(infile)
    #print "numatoms:", numatoms
    #quit()
    if debug == "yes":
        print "yo3"
    freqs = HesseMatrix_to_ExactFreqs_wo_zeroFreqs(infile=infile,output_as_nparray=True)

    if debug == "yes":
        print "ff ",freqs
    if Tmax == None: Tmax = db.MeltingPoint()
    #print "yo2",Tmax
#    for i in freqs:
#        print i
    #print "yo",Tmax
    out = []
    for i in range(int(Tmax)+1):
        out.append([i,ExactFreqs_to_Fqh_fromExactFreqs(freqs, i )])
#        print "y:",type(ExactFreqs_to_Fqh_fromExactFreqs(freqs, i ))
#        print "%.f %.16f" % (i, ExactFreqs_to_Fqh_fromExactFreqs(freqs, i ))
#    print ""
#    print list
#    if outfile == None
#    import numpy as np
#
#    np.savetxt("ka.da",list,fmt='%.0f  %20.20f')
#    return " "
    #kprint 'yes'
    if forSupercell == True:
        out = out*numatoms

    if output_as_nparray == True:
        return out

    ## from here: print to screen or print to file
    if outfile != None:
        import numpy as np
        from numpy import savetxt
        np.savetxt(str(outfile),out,fmt='%.0f  %20.20f')
        return str(outfile)+" written"
    else:
    #        print list(out)
        for i in out:
            print i[0],i[1]
    return " "

def ExactFreqs_to_Fqh_fromExactFreqs(ExactFreqs=None,T=None,includes_zero_Freqs=False):
    """
    -ceffe
    """
    ############################
    ### T
    ############################
    if T == None:
        my.exit(error="please provide Temperature")
        #Tmax = db.MeltingPoint()

    if my.is_number(T) != True:
        my.exit(error="Temp: "+str(T)+" is not a number")
    T = float(T)

    ############################
    ### Exact Freqs
    ############################
    if ExactFreqs == None: my.exit("ERROR: please specify ExactFreqs, can be either a file or a numpy.array")
    import numpy as np
    if T == 0: T = 1  ## just to avoid durch 0 teilen
    if type(ExactFreqs) == np.ndarray:
        ExactFreqs = ExactFreqs
    elif my.checkfile(str(ExactFreqs)):
        ExactFreqs = np.loadtxt(str(ExactFreqs))
    else:
        my.exit("ERROR: What is "+str(ExactFreqs)+" can be either a file or a numpy.array")
    kB=0.086173423

    #############################
    ### check if exact Freqs have 0 Frequencies included
    #############################


    ### skaliertung auf
    #print "SUM:",np.sum(ExactFreqs/2)
    #print (1-np.exp(-ExactFreqs/(kB * T)))
    #quit()
    #print ExactFreqs/2+(kB*T*np.log(1-np.exp(-ExactFreqs/(kB * T))))
    #print "T:",T

    Fqh_whole_cell = np.sum(ExactFreqs/2+(kB*T*np.log(1-np.exp(-ExactFreqs/(kB * T)))) )
    #print "Fqh_whole_cell:",Fqh_whole_cell
    ### skaliertung auf Fqh_peratom
    Fqh_per_atom = Fqh_whole_cell/(len(ExactFreqs)/3)
    return Fqh_per_atom

def create_fitFqhInput_Blazej_skripte(atoms=None,test=False,allfolder=None):
    """
    -gfi
    """
    if path == None:
        path = my.pwd()

    if atoms == None: my.exit("ERROR: for bulk or supercell? specify atoms to fit Fvib_fromExactFreqs_1 or Fvib_fromExactFreqs_32 ...")
    get.from_path_string_job_is(path=path,job="dynmat")
    get.from_path_string_details_real()  ## just to chek if details are available

    if allfolder == None:
        allfolder = ','.join([get.from_path_get_ang(my.pwd()+"/"+x) for x in os.walk('.').next()[1]])
    else:
        allfolder = allfolder

    print "allfolder:",allfolder," type:",  type(allfolder)

    job = get.from_path_string_structure(path=path)
    if "vak" or 'divak' is job:
        structureFactor = 1
    elif job is "bulk":
        my.exit("BLAZEJ: nochmal besprechen wie das hier mit dem input file gemacht wird oder abeer neues eigenes Fittigskript allgemein!")
    sc1 = get.from_path_string_details_supercell(path=subproj,sc123="1")
    sc2 = get.from_path_string_details_supercell(path=subproj,sc123="2")
    sc3 = get.from_path_string_details_supercell(path=subproj,sc123="3")
    if sc1 != sc2 !=sc3:
        my.exit(error="Supercell has different dimensions: sc1xsc2xsc3"+str(sc1)+" "+str(sc2)+" "+sc2)

    s = 1
    fitOrder = 2
    if test == True: return True

    with open("fitFqh.input", "w") as f:
        f.write("aLats = {"+allfolder+"};"+"\n")
        f.write("structureFactor = "+str(structureFactor)+";                    (*  4: fcc  2: bcc  1: supercells *) \n")
        f.write("sc = "+str(sc1)+";                                              (*  supercell, 1 for bulk *)\n")
        f.write("s = "+str(s)+";                                                (*  scaling factor for Fvib *)\n")
        f.write("fitOrder = 2;                                                  (*  typically: 2 giving: {1,V,V^2} *) \n")
        f.write("baseNames = {\"Fqh_fromExactFreqs_"+str(atoms)+"_\"};                      (*  Fqh_fromExactFreqs_  *)\n")
        f.write('<<"/home/glensk/db/utils/job__dynmat/EXE_fitFqh.math"')
    return True

def create_file_mean_freqs_from_all_HesseMatrix_(path=None):
    """
    -cmf
    """
    if path == None:
        path = my.pwd()

    import glob
    hessemats = sorted(glob.glob(str(path)+"/HesseMatrix_[0-9.]*"))
    #hessemats = sorted(glob.glob(str(path)+"/HesseMatrix_!(*[_]*)"))
    file = str(path)+"/mean_freqs"
    print "file:",file
    for i in hessemats:
        print i

    open(file,'w').close()

    alatmeanfreq = []
    for i in hessemats:
        alatstring = i.split("/")[-1].split("_")[-1]
        meanfreq = HesseMatrix_to_MeanFreqs(infile=i)
        alatmeanfreq.append(str(alatstring)+" "+str(meanfreq))
        print alatstring,meanfreq
    with open(file,'w') as f:
        for line in alatmeanfreq:
            f.write(str(line)+"\n")
    return True

def create_parametersdat_Blazej_skript(alat=None,destinationfolder=None,test=False):
    """
    -cpb
    """
    ## 0 checks: dynmat, fcc4,
    ## 0 checks
    get.from_path_string_job_is(job="dynmat")
    get.from_path_string_job_cell_is(cell="fcc4")
    get.from_path_string_structure_is(structure="bulk")

    if test == False and destinationfolder == None:
        destinationfolder = my.pwd()
    my.checkdir(destinationfolder)

    if alat == None:
        my.exit("ERROR: pleas spezify alat(s)")

    ### for parameters.dat get input
    element = get.from_path_string_element()
    sc = get.from_path_string_details_supercell(sc123=1)
    tmelt = int(db.MeltingPoint())+1

    if test == True:
        return True

    ### write file
    parametersdat = destinationfolder+"/parameters.dat"
    with open(parametersdat, "w") as j:
        j.write("element=\""+str(element)+"\"   # e.g., XX=Ca\n")
        j.write("cellType=\""+"fcc"+"\"   # XX=fcc or bcc\n")
        j.write("supercell="+str(sc)+"   # XX=2,3,4\n")
        j.write("aLats={"+str(alat)+"}     #  linst in angstrom\n")
        j.write("TRange={1,"+str(tmelt)+",1}      # optionally: start/end/step temperature T1/T2/dT in K\n")


def auswerten_bulk_fcc4_Blazejs_skript():
    """
    -ab
    """
    ## 0 checks: dynmat, fcc4,
    ## 0 checks
    get.from_path_string_job_is(job="dynmat")
    get.from_path_string_job_cell_is(cell="fcc4")
    get.from_path_string_structure_is(structure="bulk")
    get.from_path_string_details_real()  ### sicherstellen das diese definiert sind da ansonsten auswertung ueber alle 2x2x2sc_300eV... 2x2x2sc_400eV ...
    ### wenn man wirklicih alle subfolder auswerten moechte, dann obrige zeile einfach rauss

    ## get OUTCARS + create single Fvib
    ## get OUTCARS + create single Fvib
    allOUTCARS = my.run('find -L `pwd` -type f -name "OUTCAR*" | sort -n').split("\n")
    allOUTCARS = filter(None,allOUTCARS)
#    print "ok"
    alats=""
    for idx,OUTCAR in enumerate(allOUTCARS):
        #print "OUTCAR:",OUTCAR
        alat = get.from_path_get_ang(path=OUTCAR)


        #        print "al:",alat
        stringzus=","
        if idx == 0: stringzus=""
        alats=alats+stringzus+str(alat)


        forcesfile = "/".join(OUTCAR.split("/")[0:-1])+"/forces."+str(alat)
        my.run("OUTCAR_forces-last.sh "+str(OUTCAR)+" > "+str(forcesfile))
        displacement = my.run("OUTCAR_positions_last.sh "+str(OUTCAR)).split("\n")[0].split()[0]
        if displacement == "0.00529":
            pass
        else:
            my.exit("ERROR: displacement is nto 0.00529 in "+str(OUTCAR))
        dipfile = "/".join(OUTCAR.split("/")[0:-1])+"/disp."+str(alat)
        my.run("rm -f "+str(dipfile)+" ;echo 0.0052917721 > "+str(dipfile))
        destinationfolder = "/".join(OUTCAR.split("/")[0:-1])
        startsktript = "/home/grabowski/Thermodynamics/getSingleSpeciesPhonons.sh -A"
        create_parametersdat_Blazej_skript(alat=alat,destinationfolder=destinationfolder)
        print ""
        print "starting mathematica in "+str(destinationfolder)
        print my.run('hier=`pwd`;cd '+str(destinationfolder)+';pwd;'+str(startsktript)+';cd $hier; cp '+str(destinationfolder)+"/Fqh_from*Freqs* .")
        print "done ..."
        #FqhExact_file="/".join(OUTCAR.split("/")[0:-1])+"/Fqh_fromExactFreqs_"+str(alat)

        print " "
#        print alat,displacement

    ### now, fit the surface
    ### now, fit the surface
    print my.run('mmv -r "Fqh_from*Freqs_*" "Fqh_from#1Freqs_1_#2" ')
    print "alats:",alats
    create_fitFqhInput_Blazej_skripte(atoms=1,allfolder=str(alats))
    return True


def auswerten_defect_fcc4_from_HesseMatrix():
    """
    -ad
    """
    ## 0 checks: dynmat, fcc4,
    ## 0 checks
    get.from_path_string_job_is(job="dynmat")
    get.from_path_string_job_cell_is(cell="fcc4")
    get.from_path_string_structure_is(structureisnot="bulk")
    get.from_path_string_details_real()  ### sicherstellen das diese definiert sind da ansonsten auswertung ueber alle 2x2x2sc_300eV... 2x2x2sc_400eV ...
    ### wenn man wirklicih alle subfolder auswerten moechte, dann obrige zeile einfach rauss



    ##################################
    ## 1 get subfolder
    ##################################
    allfolder = [my.pwd()+"/"+x for x in os.walk('.').next()[1]]
    allHessefiles = [x+"/HesseMatrix" for x in allfolder]
#    print allfolder
#    print allHessefiles
#    print "ok   "

    ##################################
    ## 2 check if everything there
    ##################################
    for Hessefile in allHessefiles:
        if my.isfile(Hessefile) != True:
            continue
#        print "just skip if Hessefile doees not exist", Hessefile
#        print "check1"
        my.checkfile(Hessefile)
#        print "cehck2"
        get.from_path_get_ang(path=Hessefile)
        atoms = get.from_path_get_numatoms(path=Hessefile)
#        print "cehck3",atoms
        create_fitFqhInput_Blazej_skripte(test=True, atoms=atoms)
#    print 'sowei'
    ##################################
    ## 3 get Fvib
    ##################################
    alats=""
    for idx,Hessefile in enumerate(allHessefiles):
        if my.isfile(Hessefile) != True:
            print "HesseMatrix "+str(Hessefile)+" does not exist"
            continue
        alat = get.from_path_get_ang(path=Hessefile)
#        print "al:",alat
        stringzus=","
        if idx == 0: stringzus=""
        alats=alats+stringzus+str(alat)
        HesseMatrix_to_Fqh_fromExactFreqs(infile=Hessefile,outfile="Fqh_fromExactFreqs_1_"+str(alat))
        atoms = get.from_path_get_numatoms(path=Hessefile)
        HesseMatrix_to_Fqh_fromExactFreqs(infile=Hessefile,outfile="Fqh_fromExactFreqs_"+str(atoms)+"_"+str(alat))


    ##################################
    ## 4 fit Fvibs
    ##################################
    print "alats:",alats
    for atoms in [1,atoms]:
        create_fitFqhInput_Blazej_skripte(atoms=atoms,allfolder=str(alats))
        print "___________________________"*3
        #        print out.split("\n")
        print "running mathematica ..."
        out = my.run("math < fitFqh.input")
        delta=my.readfile("Fqh_fromExactFreqs_"+str(atoms)+"_fit_order2_delta")
        print delta
#        Exponent_diff = my.grep("x  delta in meV E*",out.split("\n"))
#        Exponent_diff = Exponent_diff[0].split()
#        print "EEEEEEEE:",Exponent_diff
#        Exponent_diff = float(Exponent_diff)
#        print "max diff: Exponent:"+str(Exponent_diff)+"                  (-1 is less 0.1 meV; -2 is less 0.01 meV)"
#        print "___________________________"*3
#        Exponent_diff = float(my.grep("x  delta in meV E*",out.split("\n"))[0].split()[-1])
#        if Exponent_diff < -1:   ### -1 is less 0.1 meV; -2 is less 0.01 meV
#            pass
#        else:
#            return my.exit("ERROR: diff zu hoch beim Fitten! Exponent:"+str(Exponent_diff))
    return True

def displacements_necessary(path=None):
    """
    -dn
    """
    if path == None:
        path = my.pwd()
    ### 0. ensure we are in dynmat
    get.from_path_string_job_is(job="dynmat",path=path)
    subproj = get.get_path_subproject(path=path)

    ## if between project path and ang path: -> check all necessary displacements
    ## if in displacement path *(e.g. 1_3) than on this displacement necessary

    disp = from_path_string_displacement(path=path,exit=False)
    if disp != None:
        return [disp]
    else:
        structure = get.from_path_string_structure(path=path)

        if structure != "bulk":
            ## defects -> displace first of all all atoms
            #print structure
            atoms = get.from_path_get_numatoms(path=path)
            #print [[[disp,atom] for disp in range(1,atom)] for atom in range(1,atoms)]
            #print [atom for atom in range(1,atoms)]
            #print ""
            return [str(atom)+"_"+str(disp) for atom in range(1,atoms+1) for disp in range(1,3+1)]

        else:
            my.exit(error="TODO heer")

def create_parametersdat(path=None):
    """
    -cp
    """
    if path == None:
        path = my.pwd()
    return True

def create_job_vakanz(path=None):
    """
    -cjv
    """
    if path == None:
        path = my.pwd()

    ## 1. make sure this job is a fcc4 vak job
    ## 2. get parameters dat infos
    ## 3. check check eqList_vacancy_$sc\x$sc\x$sc\sc
    ## 4. check and get referencepath to relaxed coords
    ## 5. create parameters.dat
    ## 6. get the KPOINTS with xxxKPxxx and INCAR with xxx???xxx


    ## 1. make sure this job is a fcc4 vak job
    get.from_path_string_job_is(path=path, job = "dynmat")  ## just to ensure we are in dynmat
    get.from_path_string_structure_is(path = path, structure = "vak")   ## just to ensure we are in vak job
    get.from_path_string_cell_is(path = path, cell = "fcc4")   ## so far we just know fcc4 vak symmetries
    get.check_path_is_subprojectpath(path = path, exit = True) ## you have to be in subproject so far...

    ## 2. get parameters.dat infos
    ALATS = get.from_path_get_ALATS(path=path)
    disp = 0.03
    cutoff = get.from_path_string_details_cutoff(path = path)
    get.from_path_string_details_supercell_ensure_sc123_isequal(path = path )  ## just to ensure 2x2x2 or 3x3x3 but not 3x4x5
    sc = get.from_path_string_details_supercell(sc123 = "1")
    nbands = get.from_path_get_nbands_thisstruct(path = path)
    ngxf = None
    if sc == "2": ngxf = "120"
    if sc == "3": ngxf = "180"
    if sc == "4": ngxf = "240"
    if sc == "5": ngxf = "300"
    ngxffrompath = get.from_path_string_details_ngxf(path = path)
    if ngxffrompath != None:
        if ngxffrompath != ngxf:
            my.exit(error = "ngxf in parameters: "+str(ngxf)+ "ngxf frompath: " + str(ngxffrompath))
    if ngxf == None: my.exit(error = "which ngxf?")
    kp = get.from_path_string_details_kpoints(path = path)
    if kp[0] != kp[1]:
        my.exit(error = "invest some time for kpoints 0 1")
    if kp[0] != kp[2]:
        my.exit(error = "invest some time for kpoints 0 2")
    kpoints = str(kp[0])

    print "ALATS:", ALATS
    print "cutoff:", cutoff
    print "sc:", sc
    print "nbands:", nbands
    print "ngxf:", ngxf
    print "kp:", kpoints

    ## 3. check check eqList_vacancy_$sc\x$sc\x$sc\sc
    eqList = str(get.set_path_db_utils())+"/cell__fcc4/eqList_vacancy_"+sc+"x"+sc+"x"+sc+"sc"
    #target = str(path)+"/eqList_vacancy_"+sc+"x"+sc+"x"+sc+"sc"
    my.checkfile(eqList)
    import shutil
    shutil.copy(eqList, path)
    print "eqList:", eqList

    ## 4. check and get referencepath to relaxed coords
    folder = get.from_path_get_ALATS_VOLS_folderlist_to_create_for_subproj_or_angvoltemp(path = path,rise_error_if_exists = True,silent=True)
    for i in folder:
        #print "I:",i
        out = murn.get_path_to_relaxedcoords_cartes_angvoltempfolder(thisjob = i , silent = True)
        #print "yo",i
        #print "out:",out
        #shutil.copy(out,)
        #quit()
        ang = get.from_path_get_ang(path = i)
        #print "out:",out,path+"/relaxedCoords_"+str(ang)
        shutil.copy(out,path+"/relaxedCoords_"+str(ang))
        #print "ang:",ang, "folder:",i
    #print "folder:", folder


    ## 5. create parameters.dat
    with open(path+"/parameters.dat",'w') as f:
        f.write("aLats="+str(" ".join(ALATS))+"\n")
        f.write("disp="+str(disp)+"\n")
        f.write("cutoff="+str(cutoff)+"\n")
        f.write("\n")
        f.write("\n")
        f.write("sc="+str(sc)+"           # 2sc  3sc      4sc  5sc\n")
        f.write("nbands="+str(nbands)+"     # 150  480/560  1200 2400   MP/Fermi\n")
        f.write("ngxf="+str(ngxf)+"       # 120  180  240  300\n")
        f.write("kp="+str(kpoints)+"          # 6    4    3    2\n")
        f.write("kpshift=0      # 0    0    .5   0\n")
        f.write("ediff=1E-8     # 1E-8 1E-8 1E-8 1E-8\n")
        f.write("sigma=0.1\n")
        f.write("smear=-1\n")

    ## 6. get KPOINTS / INCAR
    get.create_POTCAR(path = path)
    shutil.copy(str(get.set_path_db_utils()+"/job__dynmat/KPOINTS_vorlage_vak"),"KPOINTS")
    shutil.copy(str(get.set_path_db_utils()+"/job__dynmat/INCAR_vorlage_vak"),"INCAR")
    shutil.copy(str(get.set_path_db_utils()+"/job__dynmat/createFolders_vak.sh"),"createFolders_vak.sh")
    my.run("./createFolders_vak.sh")



def create_jobs(path=None):
    """
    -cj
    """
    if path == None:
        path = my.pwd()
    ### 0. ensure we are in dynmat
    ### 1. get subproj (later we can change this to start from proj)
    ### 2. get Alats (if in Alat folder take this alat, else from paht get ALATS)
    ### 3. get corresponding Murn subproj
    ### 4. check if corresponding relaxedcoords are available
    ### 5. get necessary displacements
    ### 6. check if jobs can be created
    ### 7. create corresponding jobs && the reference if necessary

    ### 0.) 1.)
    get.from_path_string_job_is(job="dynmat",path=path)
    subproj = get.get_path_subproject(path=path,Error_if_notexist=False)

    ### 2)
    #print "2",path
    #print ""
    ang_given = get.from_path_get_ang(path=path,exit=False)
    if ang_given == None:
        alats = get.from_path_get_ALATS(path=path,string=True)
    else:
        alats = [ang_given]
    #print "3",alats
    ### 3) (if defect -> consider same kpoints,cutoff; if bulk: hierfuer sollten wir nicht alle auslenkungen rechnen)
    structure = get.from_path_string_structure(path=path)
    if structure == "bulk": my.exit(error="for bulk: hierfuer sollten wir nicht alle auslenkungen rechnen")
    reference_murn_project = murn.get_reference_murn_project(path=path)
    current_cutoff = get.from_path_string_details_cutoff(path=path)
    current_kpoints = get.from_path_string_details_kpoints(path=path)
    reference = get.get_path_subproject_giving_projectpath_cutoff_kpoints(projectpath=reference_murn_project,cutoff=current_cutoff,kpoints=current_kpoints,exclude_referencefolder=True)
    #print "ref:",reference

    ### 4)
    import glob
    for run in ["test","notest"]:
        for alat in alats:
            #print "a",alat

            checkpath = glob.glob(reference+"/"+str(alat)+"*[Aa][Nn][Gg]")
            if len(checkpath) != 1: my.exit(error="no or more than one checkpath")
            relaxed_coords_path = checkpath[0]
            my.checkdir(relaxed_coords_path) ## path has to exist
            #/relaxed_coords_direct")
            #print "relaxed_coords_path:",relaxed_coords_path
            if my.isfile(relaxed_coords_path+"/relaxed_coords_direct") == True:
                relaxed_coords_file_red = relaxed_coords_path+"/relaxed_coords_direct"
            else:
                my.run("CONTCAR_positions_direct_write_to_path_whenfinished_and_converged.sh "+str(relaxed_coords_path))
                if my.isfile(relaxed_coords_path+"/relaxed_coords_direct") == True:
                    relaxed_coords_file_red = relaxed_coords_path+"/relaxed_coords_direct"
                else:
                    my.exit(error="NOT FOUND "+relaxed_coords_path+"/relaxed_coords_direct")
            print ">> relaxed_coords_file_red:",relaxed_coords_file_red

        ### 5) displacement necessary depends on path:
        disps = displacements_necessary(path=path)
        #print "disps:",disps

        ### 6)
        #print "%%%%%%%%%%%%%%%%%%%%%%"*8
        if len(disps) > 1: folderlist = [subproj+"/"+str(i)+"Ang"+"/"+str(disp) for i in alats for disp in disps]
        #print "1",folderlist,disps,type
        #print ""
        #print ""
        if len(disps) == 1: folderlist = [subproj+"/"+str(i)+"Ang"+"/"+disps[0] for i in alats]
        for folder in folderlist:
            print ">>",folder
            get.check_jobfolder_to_create(onepath=folder,writeinfo=True,direct_position_file=relaxed_coords_file_red)



def EVinet_read(path=None):
    """
    -er
    """
    if path == None:
        path = my.pwd()

    EVinet = my.readfile(path=path+"/EVinet")[0].split()
    if len(EVinet) != 4:
        my.exit(error="len EVinet is not 4")
    E0=float(EVinet[0])
    V0=float(EVinet[1])
    BM=float(EVinet[2])
    BMder=float(EVinet[3])
    return E0,V0,BM,BMder

def fit_EVinet(path=None):
    """
    -fe
    """
    pass

def fit_surface(path=None):
    """
    -fs das ist erstmal nur spielerei
    """
    if path == None:
        path = my.pwd()

    E0, V0, BM, BMder = EVinet_read(path=path)
    print murn.Vinet(V=V0,E0=E0,V0=V0,BM=BM,BMder=BMder)



    return E0, V0, BM, BMder



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
