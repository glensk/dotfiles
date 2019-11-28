#!/usr/bin/python
import get, my, argparse, sys, os, db, dynmat
###############################################################################
## definitions start
## V V V V V V V V V
###############################################################################
all_builtin =  dir() ## necessat to define all defs in this script
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
def __POTENTIAL__o__ELEMENT__o___PROJECT__________o__SUBPROJECT_o__ANGVOLTEMPFOLDER():
    pass
def __POTENTIAL__o__ELEMENT__o___job_structure_cell__o__DETAILS______o__ANGVOLTEMPFOLDER():
    pass
def __POTENTIAL__o__ELEMENT__o___ti_divak_fcc4_______o__DETAILS______o__ANGVOLTEMPFOLDER():
    pass
def __POTENTIAL__o__ELEMENT__o___ti_divak_fcc4_______o__ref_high_____o__ref_high_job__o():
    pass

def check_if_lowpath(path=None,exit=True):
    """
    -cil
    """
    if path == None:
        path = my.pwd()
    get.from_path_string_job_is(path=path,job="ti")
    str = get.from_path_string_details().split("_")[0]
    if str == "low":
        pass
    else:
        return my.exit("ERROR: You are not in a PROJECTpath (here: not in low path)")

    #### check if is in Highfolder!
    if "high" in get.from_path_string_details().split("_"):
        return my.exit("ERROR: You are HIGH folder")
    else:
        return True

def check_if_highpath(path=None):
    """
    -cih
    """
    if path == None:
        path = my.pwd()

    get.from_path_string_job_is(path=path,job="ti")
    #    print "KK",get.from_path_string_details().split("_")
    if "high" in get.from_path_string_details().split("_"):
        return True
    else:
        return my.exit("ERROR: You are not in high folder")


def get_path_to_HesseMatrix_folder(path=None):
    """
    -gph
    """
    if path == None:
        path = my.pwd()
    get.from_path_string_job_is(path=path,job="ti")
    sc = get.from_path_string_details_supercell(path=path)
    proj = get.get_path_project(path=path)
    hessefolder = proj+"/Hessematrix_"+str(sc)
    return my.checkdir(hessefolder)

def get_path_mean_freqs(path=None):
    """
    -gpmf
    """
    if path == None:
        path = my.pwd()
    hessematfolder = get_path_to_HesseMatrix_folder(path=path)
    file = hessematfolder+"/mean_freqs"

    if my.isfile(file) != True:
        dynmat.create_file_mean_freqs_from_all_HesseMatrix_(path=hessematfolder)
    return my.checkfile(hessematfolder+"/mean_freqs")

def get_path_Fah_file(path=None):
    """
    -gpfah
    """



def low_create_job():
    """
    -lcj
    """
    if get.create_POTCAR() != True: return sys.exit("get.create_POTCAR() problem")
    if get.create_KPOINTS() != True: return sys.exit("get.create_KPOINTS() problem")
    if get.create_POSCAR() != True: return sys.exit("get.create_POSCAR() problem")

    return sys.exit("TODO") ##TODO

def low_create_parameters_file():
    return sys.exit("TODO") ##TODO

def low_check_parameters():
    """
    -lcp
    """
    if check_if_lowpath() != True:
    #        print "ERROR: to check low parameters you have to be in low path"
        sys.exit("ERROR: to check low parameters you have to be in low path")
    par = my.readfile(my.pwd()+"/parameters.dat")
    return par

def get_path_ref_high(path=None):
    """
    -hrp  e.g. ../ref_high_2x2x2sc or ../ref_high_3x3x3sc
    """
    if path == None:
        path = my.pwd()

    #print "2",path
    ### check ob wir in ti sind
    get.from_path_string_job_is(path=path,job="ti")

    #print "3        "
    ## get sc for refpath
    sc = get.from_path_string_details_supercell(path=path)


    ## get /home/glensk/v/PAW_PBE/Al/ti_divak_fcc4 oder so, der pfad in welchem
    pathout = get.get_path_job_type_cell(path=path)

    #print "4",path
    ref = my.checkdir(pathout+"/ref_high_"+sc,create=True)
    return ref

def get_path_ref_high_job(path=None):
    """
    -prhj
    """
    if path == None:
        path = my.pwd()
    get.from_path_string_job_is(path=path,job="ti")

    import re
    ref = get_path_ref_high(path=path)
    cut = get.from_path_string_details_cutoff(path=path)
    kp = get.from_path_string_details_kpoints(path=path)
    availall = get_path_ref_high_jobs_avail(path=path)
    #print "len:",len(avail)
    #    print "cut:",cut,"kp:",kp
    if not len(availall):
        avail = my.checkdir(ref+"/"+cut+'eV_'+kp[0]+"x"+kp[1]+"x"+kp[2]+"kp", create=True)
    #        print "no ref path in:",ref, "cut:",cut,"kp:",kp
    #        my.exit()
    out =[]
    #print "avail:",avail
    for avail in availall:
    #    print ":::",avail
    #        print
        ## cutoff check
        if re.search(cut+"eV", avail) != None:
            #kp chekck
            if re.search("[0]*"+kp[0]+"x"+"[0]*"+kp[1]+"x"+"[0]*"+kp[2]+"[kK][pP]",avail) != None:
                #print avail
                #add = re.search(cut+"eV", avail).group()
                out.insert(1,avail)
    #print "out:",out
    if len(out) > 1:
        ngxf = get.from_path_string_details_ngxf(path=path)
        outtmp = []
        for avail in out:
            #print avail
            jo1=None
            jo2=None
            if re.search(str(ngxf)+"NGXF*",avail) != None:
                jo1 = avail
            if re.search("NGXF"+str(ngxf),avail) != None:
                jo2 = avail
            if jo1 == None and jo2 == None: continue
            if jo1 != None and jo2 != None: continue
            outtmp.insert(1,avail)
        out = outtmp
    #print "outtmp:",outtmp

    if len(out) == 0:
        my.exit("ERROR: no ref paths job, ref path out: "+ref)

    if len(out) != 1:
        my.exit("found several ref paths: "+str(out))
    pathout = my.checkdir(ref+'/'+out[0])
    return pathout

def get_path_ref_high_job_ang(path=None):
    """
    -hrpja
    """
    if path == None:
        path = my.pwd()
    get.from_path_string_job_is(path=path,job="ti")

    pathref = get_path_ref_high_job(path=path)
    return my.exit("ERROR: TODO get Ang"+pathref)

def get_path_ref_high_jobs_avail(path=None):
    """
    -hra
    """
    if path == None:
        path = my.pwd()
    thedir = get_path_ref_high(path=path)
    return [ name for name in os.listdir(thedir) if os.path.isdir(os.path.join(thedir, name)) ]


def high_get_avail_temperatures(path=None):
    """
    -gat
    """
    if path == None:
        path = my.pwd()

    check_if_highpath(path=path)
    #print "1"
    highpaths = get.get_path_all_angvoltempfolder(path=path)
    #print "2"
    angs = []
    temps = []
    for highpath in highpaths:

        str = highpath.split("/")[-1].split("Ang_")
        alat = str[0]
        T = str[1].split("K")[0]
        if my.is_number(alat) != True:
            print "alat problem: "+str(alat)
        if my.is_number(T) != True:
            print "Temperatrue problem: "+str(T)
        angs.append(alat)
        temps.append(T)
    #print sorted(set(angs))
    return sorted(set(temps))

def high_get_avail_angs(path=None):
    """
    -gaa
    """
    if path == None:
        path = my.pwd()

    check_if_highpath(path=path)
    #print "1"
    highpaths = get.get_path_all_angvoltempfolder(path=path)
    #print "2"
    angs = []
    temps = []
    for highpath in highpaths:

        str = highpath.split("/")[-1].split("Ang_")
        alat = str[0]
        T = str[1].split("K")[0]
        if my.is_number(alat) != True:
            print "alat problem: "+str(alat)
        if my.is_number(T) != True:
            print "Temperatrue problem: "+str(T)
        angs.append(alat)
        temps.append(T)
    return sorted(set(angs))
    #return sorted(set(temps))

def high_create_auswertung_high(path=None):
    """
    -hcah
    """
    if path == None:
        path = my.pwd()

    ## 1. ensure to be in ti und high
    get.from_path_string_job_is(job="ti",path=path)
    get.get_path_subproject(path=path) ## just to ensure we have a subproj

    ## 2. get all angK folder
    angKfolder = get.get_path_all_angvoltempfolder(path=path)
    Fah_high_files = [folder+"/Fah_high" for folder in angKfolder]

    enes = ['eS0','ene','fre']
    fits = ['linear','quad','best','mitte','linear']

    for i in Fah_high_files:
        print i

def high_create_changeHesseref_folder(oldfolder=None,newfolder=None,path=None,test=False):
    """
    -hch
    """
    if path == None:
        path = my.pwd()

    if oldfolder == None or newfolder == None:
        my.exit(error="you have to provide oldfolder and newfolder")
    if get.get_path_element(path=oldfolder) != get.get_path_element(path=newfolder):
        my.exit(error="path to element oldfolder "+str(newfolder)+" is not path to element newfolder "+str(newfolder))
    if get.from_path_get_numatoms(oldfolder) != get.from_path_get_numatoms(newfolder):
        my.exit(error="numatoms oldfolder "+str(newfolder)+" is not numatoms newfolder "+str(newfolder))
    if get.from_path_string_cell(oldfolder) != get.from_path_string_cell(newfolder):
        my.exit(error="cell oldfolder "+str(newfolder)+" is not cell newfolder "+str(newfolder))


    Fahsurface = my.checkfile(str(path)+"/Fah_fre_best_surface")

    ########################################################################################
    ## 1. get alats Fah_surface
    Fahsurfaceread = my.readfile(Fahsurface)
    alatlist = []
    temperatures = []
    for line in Fahsurfaceread:
        words = line.split()
        if len(words) < 2: continue
        temp = words[0]
        alat = words[1]
        if my.is_number(alat) != True: continue
        if my.is_number(temp) != True: continue
        alatlist.append(float(alat))
        temperatures.append(int(temp))
        #print temp,alat,my.is_number(alat) #words,len(words)
    alatlist = list(set(alatlist))   ## remove duplicates
    temperatures = list(set(temperatures))   ## remove duplicates
    if test != False:
        print "Fah_surface alats       :",alatlist
        print "Fah_surface temperatures:",temperatures
        print ""

    ########################################################################################
    ## 2. check (old) folder

    ########################################################################################
    ## 3. check hessematfolder: Fqh_fromExactFreqs_[alats] exist? with corresponding temperatures?
    ## if not:
    ##              /home/grabowski/Thermodynamics/collectOUTCARs.sh  !!IT REMOVES STUFF!!
    ##              /home/grabowski/Thermodynamics/extractForces.sh
    ##              /home/grabowski/Thermodynamics/getSingleSpeciesPhonons.sh

    for folder in [oldfolder,newfolder]:
        if folder == None: my.exit(error="you have to provide the folder to check")
        print my.error_write(errortext="THSI IS CURRENTLY ONLY DONE FOR EXACT FREQS (for vak/divak)")
        if test != False: print "CHECKING folder:",folder
        for alat in alatlist:  ### check every /nas/glensk/v/pp/al/ti_bulk_fcc4/Hessematrix_2x2x2sc/Fqh_fromExactFreqs_x.xx file
            fqhexact = str(folder)+"/Fqh_fromExactFreqs_"+str(alat)
            my.checkfile(fqhexact)
            fqhexact_read = my.readfile(fqhexact)
            tempsavail = []
            for line in fqhexact_read:
                if len(line.split()) != 2: continue
                temp = line.split()[0]
                if my.is_number(temp) != True: continue
                tempsavail.append(int(line.split()[0]))
                #print line.split()
            for t in temperatures:
                if t not in tempsavail:
                    my.exit(error="Temperature "+str(t)+" is not in "+str(fqhexact))
        if test != False:print "CHECKING folder:",folder," is OK to changeHesse"
        if test != False:print ""

    proj = get.from_path_string_project(path=newfolder)
    subproj = get.from_path_string_subproject(path=newfolder)
    if test == True: return True

    ########################################################################################
    ## 4. get + create new changeHesseFolder
    createfolder = path+"/changeHesseref___"+proj+"___"+subproj
    import glob,shutil
    if my.isdir(createfolder) == True:
        shutil.rmtree(createfolder)
    os.makedirs(createfolder)

    if test != False:print "createfolder:",createfolder
    if test != False:print "path:",path

    for file in glob.glob(path+"/Fah_*"):
        last = file.split("_")
        if len(last) <= 1: continue
        last = last[-1]
        shutil.copyfile(file,createfolder+'/Fah_'+str(last))
        #print last

    os.chdir(createfolder)
    my.run2(command='rm -f Fah_from_fit_tangens; echo 4 > Fah_from_fit_tangens')
    my.run2(command='/home/grabowski/Thermodynamics/changeHesseReference.sh '+str(oldfolder)+" "+str(newfolder))
    return True




def high_oooooooooooooo_create_Fah_surface_changedHesseref(path=None,test=False):
    """
    -hcch which changehesseref do exist?
    ##from            ->  to
    ##this bulk    -> new_current_bulk
    ##this bulk    -> vak__reference     e.g. --dynmat_vak_fcc4---2x2x2sc_300eV_NGXF120_ADDGRID_16x16x16kP_TAKE__reference
    ##this bulk    -> divac__reference   e.g. --dynmat_divak_fcc4---
    ##this vak     -> new_current_vak
    ##this vak     -> bulk__reference
    ##this divak   -> new_current_divak
    ##this divak   -> bulk__reference
    ---> old ist always current Hessematrix_xxxx folder!
    """
    if path == None:
        path = my.pwd()

    ########################################################################################
    ## 0. check for necessary auswertung_lowplushighFahbest/
    check_if_highpath(path=path)  ## just to ensure we deal with a high folder
    subproj = get.get_path_subproject(path=path)  ## get path to (high) subproject

    ########################################################################################
    ## 0. check for necessary auswertung_lowplushighFahbest/(files) exist
    Fahsurface = my.checkfile(str(path)+"/Fah_fre_best_surface")
    print "Fahsurface:",Fahsurface

    ########################################################################################
    ## 1. get alats Fah_surface
    Fahsurfaceread = my.readfile(Fahsurface)
    alatlist = []
    temperatures = []
    for line in Fahsurfaceread:
        words = line.split()
        if len(words) < 2: continue
        temp = words[0]
        alat = words[1]
        if my.is_number(alat) != True: continue
        if my.is_number(temp) != True: continue
        alatlist.append(float(alat))
        temperatures.append(int(temp))
        #print temp,alat,my.is_number(alat) #words,len(words)
    alatlist = list(set(alatlist))   ## remove duplicates
    temperatures = list(set(temperatures))   ## remove duplicates
    print "alats:",alatlist
    print "temperatures:",temperatures

    ########################################################################################
    ## 2. get this (old) folder
    hessematfolder = get_path_to_HesseMatrix_folder(path=path)
    print "Hessematrix folder (HesseMatrix from) : ",hessematfolder

    ########################################################################################
    ## 3. check hessematfolder: Fqh_fromExactFreqs_[alats] exist? with corresponding temperatures?
    ## if not:
    ##              /home/grabowski/Thermodynamics/collectOUTCARs.sh  !!IT REMOVES STUFF!!
    ##              /home/grabowski/Thermodynamics/extractForces.sh
    ##              /home/grabowski/Thermodynamics/getSingleSpeciesPhonons.sh
    print ""
    for alat in alatlist:  ### check every /nas/glensk/v/pp/al/ti_bulk_fcc4/Hessematrix_2x2x2sc/Fqh_fromExactFreqs_x.xx file
        fqhexact = str(hessematfolder)+"/Fqh_fromExactFreqs_"+str(alat)
        my.checkfile(fqhexact)
        fqhexact_read = my.readfile(fqhexact)
        tempsavail = []
        for line in fqhexact_read:
            if len(line.split()) != 2: continue
            temp = line.split()[0]
            if my.is_number(temp) != True: continue
            tempsavail.append(int(line.split()[0]))
            #print line.split()
        for t in temperatures:
            if t not in tempsavail:
                my.exit(error="Temperature "+str(t)+" is not in "+str(fqhexact))
            #else:
            #    print alat,t,"ok"

    ########################################################################################
    ## 3. get NEW folder for changeHesseref
    ##    a) NEW folder muss gleiche anzahl an atomen haben
    ##    b) mussen gleichen element path haben e.g. ~/v/pp/al
    ##    c) muss ein dynmat job sein
    ## if bulk -> make vak (if exists) : get sc and get dynmat_fcc4/samesc* TAKE*__reference
    ## if bulk -> make divak (if exists)

    structure = get.from_path_string_structure(path=path)
    elementpath = get.get_path_element(path=path)
    cell = get.from_path_string_cell(path=path)
    print "str:",structure
    if structure == "bulk":
        for structure in ['vak','divak']:
            print "changeHessetoTAKE",structure,"-->    TODO ",elementpath+"/dynmat_"+structure+"_"+str(cell)

    newfolder =  "/home/glensk/v/pp/al/dynmat_vak_fcc4/2x2x2sc_300eV_NGXF120_ADDGRID_16x16x16kP__reference"



    return "TODO",hessematfolder

def high_fit_Fah_surface(path=None,test=False):
    """
    -hfs
    """
    if path == None:
        path = my.pwd()

    ########################################################################################
    ## 1. check for necessary auswertung_lowplushighFahbest/
    check_if_highpath(path=path)  ## just to ensure we deal with a high folder
    subproj = get.get_path_subproject(path=path)  ## e.g. /nas/glensk/v/pp/al/ti_bulk_fcc4/low_2x2x2sc_250eV_03x03x03kp_EDIFF1E-2_TAKE__high_400eV_6x6x6kp_TAKE
    #print "subproj:",subproj

    ########################################################################################
    ## 2. check for necessary auswertung_lowplushighFahbest/(files) exist
    Fahsurface = my.checkfile(str(path)+"/Fah_surface")

    ########################################################################################
    ## 3. get all necessary information for Fah_surface.input
    atoms = get.from_path_get_numatoms(path=subproj)
    tmelt = db.MeltingPoint(path=subproj)

    job = get.from_path_string_structure(path=subproj)
    print "job:",job," if job is bulk we should do both: vakrefsurface and bulk surface: here: vak(defect) reference"
    if job is "bulk":
        filename = "Fah_b_"+str(atoms)
        filename = "Fah" ## and Fah_surface_fit

    if "vak" or 'divak' is job:
        structureFactor = 1
        filename = "Fah_d_"+str(atoms)

    elif job is "bulk":
        my.exit("BLAZEJ: nochmal besprechen wie das hier mit dem input file gemacht wird oder abeer neues eigenes Fittigskript allgemein!")
    sc1 = get.from_path_string_details_supercell(path=subproj,sc123="1")
    sc2 = get.from_path_string_details_supercell(path=subproj,sc123="2")
    sc3 = get.from_path_string_details_supercell(path=subproj,sc123="3")
    if sc1 != sc2 !=sc3:
        my.exit(error="Supercell has different dimensions: sc1xsc2xsc3"+str(sc1)+" "+str(sc2)+" "+sc2)

    Fahsurfaceread = my.readfile(Fahsurface)
    alatlist = []
    for line in Fahsurfaceread:
        words = line.split()
        if len(words) < 2: continue
        temp = words[0]
        alat = words[1]
        if my.is_number(alat) != True: continue
        if my.is_number(temp) != True: continue
        alatlist.append(float(alat))
        #print temp,alat,my.is_number(alat) #words,len(words)
    alatlist = list(set(alatlist))   ## remove duplicates
    Fahsurfacemaxalat = max(alatlist)
    Fahsurfaceminalat = min(alatlist)


    if test == True:
        return True  ## we already have all information necessary

    ########################################################################################
    ## 4. ceck if mean freqs exist and contain correct alats
    meanfreqfile = get_path_mean_freqs(path=subproj)
    meanfreqfilelist = my.readfile(meanfreqfile)
    alats_in_meanfreqfilelist = []
    for line in meanfreqfilelist:
        words = line.split()
        if len(words) <= 1: continue
        alatcheck = words[0]
        if my.is_number(alatcheck) != True: continue
        alats_in_meanfreqfilelist.append(alatcheck)
    #print "alats_in_meanfreqfilelist",alats_in_meanfreqfilelist


    for alatcalculated in alatlist:
        #print "calculatedalat",alatcalculated
        if str(alatcalculated) in alats_in_meanfreqfilelist:
            pass
        else:
            my.exit(error="Alat "+str(alatcalculated)+" not in "+str(meanfreqfile))

    import shutil
    #print meanfreqfile
    #print path
    shutil.copyfile(meanfreqfile,path+"/mean_freqs")

    ########################################################################################
    ## 5. create Fah_surface.input
    file = str(path)+"/Fah_surface_input"

    if my.isfile(file) == True:
        os.remove(file)

    with open(file, "w") as f:
        f.write('(* adjustable parameters start *)'+"\n")
        f.write("\n")
        f.write('FsurfFile = "Fah_surface";'+"\n")
        f.write('type = 1                                               (*  1: T(K)  aLat(Ang/at)  F(meV/at)   *)'+"\n")
        f.write('                                                       (*  2: T(K)  V(Ang/at)     F(meV/at)   *)'+"\n")
        f.write('                                                       (*  3: T(K)  V(Ang/cell)   F(meV/cell) *)'+"\n")
        f.write("\n")
        f.write("\n")
        f.write('min = '+str(Fahsurfaceminalat)+";                      (*  aLat or volume range (same format as Fsurf) for the 2. fit *)"+"\n")
        f.write('max = '+str(Fahsurfacemaxalat)+";                      (*  typically: Vmin=Veq(T=0K) and Vmax=Veq(Tmelt) *)"+"\n")
        f.write("mesh = 100;                                            (* 100 is good and should be kept *)"+"\n")
        f.write("\n")
        f.write("structureFactor = "+str(structureFactor)+";                                       (*  4: fcc  2: bcc  1: supercells *)"+"\n")
        f.write("sc = "+str(sc1)+";                                                    (*  supercell, 1 for bulk *)"+"\n")
        f.write("nAtoms = "+str(atoms)+";                                               (*  1 for bulk *)"+"\n")
        f.write("\n")
        f.write("fitType = \"Fvib\";                                          (*  \"Fvib\"  or  \"poly\"  fit type for 1. fit; take \"Fvib\" for Fah or Fel *)"+"\n")
        f.write("basis[V_, T_] := {1,T, V }                                 (*  for Fah typically: \"Fvib\" and {1,T,V} *)"+"\n")
        f.write("\n")
        f.write("                                                           (*  for Fel typically: \"Fvib\" and {1,T, V,T V,T^2,V^2,V^3} *)"+"\n")
        f.write("basis2[V_]:={1, V, V^2, V^3}                               (*  should be more than sufficient: {1,V,V^2,V^3} *)"+"\n")
        f.write("\n")
        f.write("minT = 1;                                                  (* typically 1     *)"+"\n")
        f.write("maxT = "+str(tmelt)+";                                               (*           Tmelt *)"+"\n")
        f.write("stepT = 2;                                                 (*           2     *)"+"\n")
        f.write("\n")
        f.write("useMeanFreqs=True;                                         (* if True \"mean_freqs\" file must be available; format as Fsurf, e.g. aLat(Ang) meanFreq(meV) *)"+"\n")
        f.write("                                                           (* meanFreqs are then used in the fit formula (check fitSurface.math) *)"+"\n")
        f.write("(* adjustable parameters end *)"+"\n")
        f.write("\n")
        f.write("(*<<\"~/scripts/Thermodynamics/fitSurface.math\"*)"+"\n")
        f.write("<<\"/home/grabowski/Thermodynamics/mathematica/fitSurface.MATH\""+"\n")

    if os.path.exists(path+str("/Fah_surface_fit")) == True:
        os.remove(path+str("/Fah_surface_fit"))

    print my.run2("math < Fah_surface_input")
    my.checkfile("Fah_surface_fit")
    import shutil
    shutil.copyfile(path+"/Fah_surface_fit",path+"/"+str(filename))
    return True




def high_create_results(path=None):
    """
    -hcr
    """
    if path == None:
        path = my.pwd()

    ########################################################################################
    ## 1. check for necessary auswertung_lowplushighFahbest/
    check_if_highpath(path=path)  ## just to ensure we deal with a high folder
    subproj = get.get_path_subproject(path=path)  ## get path to (high) subproject
    ## subproj e.g. /nas/glensk/v/pp/al/ti_bulk_fcc4/low_2x2x2sc_250eV_03x03x03kp_EDIFF1E-2_TAKE__high_400eV_6x6x6kp_TAKE
    subproj_lowplushigh = my.checkdir(subproj+"/auswertung_lowplushighFahbest")

    surface = my.checkfile(subproj_lowplushigh+"/Fah_fre_best_surface")
    import glob
    Fah_ang_files = glob.glob(subproj_lowplushigh+"/Fah_fre_best_*Ang")
    Fah_K_files = glob.glob(subproj_lowplushigh+"/Fah_fre_best_*K")

    ########################################################################################
    ## 1. check which results exist and get nextresults folder
    results = glob.glob(subproj+"/results*")
    if len(results) == 0:
        nextresults = subproj+"/results1"
        print "next:",nextresults
    else:
        i=0
        while i < 20:
            i = i+1
            if i == 19:
                my.exit(error="Hier gibt es schon zu viele results folder")
            nextresults = subproj+"/results"+str(i)

            #print "check",nextresults
            if my.isdir(nextresults) == True:
                continue
            else:
                break
    print "nextresults --> ",nextresults
    #print ""


    ########################################################################################
    ## 1. create results folder
    os.makedirs(nextresults)
    import shutil
    for file in glob.glob(subproj_lowplushigh+"/Fah_fre_best_*[Ang,K,surface]"):
        #print "KKKKKKK",file.split("_"),len(file.split("_"))
        if len(file.split("_")) <= 1: continue
        last = file.split("_")[-1]
        if last == "changedHesseref": continue
        shutil.copy(file,nextresults)
        #print last

    return True
    #return subproj_lowplushigh,Fah_ang_files


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

