#!/usr/bin/python


import my
import db
import os
import sys
import argparse
import murn

    ###############################################################################
## definitions start
## V V V V V V V V V
###############################################################################
all_builtin = dir()  # necessat to define all defs in this script


def testall():
    """
    -ta test all of the following optional arguments
    """
    out = alldefs()
    #    print "out:",out
    out.remove('alldefs')
    out.remove('testall')
    #   print "out;",out
    #   print ""
    for i in out:
        try:
            ha = str(
                out.index(i) + 1) + "/" + str(len(out)).ljust(7) + i.ljust(45)
            hb = str(eval(i)()).ljust(70)
            print "OK " + ha + hb
        except:  # catch *all* exceptions
            e = sys.exc_info()[1]
            f = str(e)
            #print hbe = sys.exc_info()[0]
    #            print "e:",e,type(e),str(e)
            #noinspection PyUnboundLocalVariable
            my.print_red_bold("-> " + ha + f.ljust(45))

    my.print_black_bold(
        "-------------------------------> testall <-- finished")
    return ""


def alldefs(script=None):
    """
    -gdcs get all definitions of this script
    """
    return my.get_definitions_of_current_script(script)

#######################
## path definitions
#######################


def __POTENTIAL__o__ELEMENT__o___PROJECT_____________o__SUBPROJECT_o__ANGVOLTEMPFOLDER():
    pass


def __POTENTIAL__o__ELEMENT__o___job_structure_cell____o__DETAILS______o__ANGVOLTEMPFOLDER():
    pass


def __POTENTIAL__o__ELEMENT__o___murn_divak_fcc4_____o__DETAILS______o__ANGVOLTEMPFOLDER():
    pass


def set_path_vasp_potentials(path=os.environ['HOME'] + "/Thermodynamics/vasp_potentials"):
    return my.checkdir(path)


def get_path_vasp_potentials_name(pot=None):
    """
    -ssss make you path recognize the potentials in your set_path_vasp_potentials folder
            in FUNCTIONALNAME of WORKING FOLDER --> out FUNCTIONALNAME of set_path_vasp_potentials
            in pp                               --> out pp
            in pp                               --> out PAW_PBE
            in PAW_PBE                          --> out pp
    """
    if pot is None:
        my.exit("please specify a pot e.g. pp, pl, PAW_PBE, ...")
    #print "pot: ", pot
    list = [["pp", "PAW-GGA-PBE"], ["pl", "PAW-LDA"]]
    #print list

    pots_avail = avail_pots()
    for avail in pots_avail:
        #print "list:  ", list
        #print "avail: ", avail
        for i in list:
            #print "i0:",i[0],"i1:",i[1],"pot:",pot,":",type(i[0]),type(pot)
            if str(pot) == i[0]:
                return i[1]
            if str(pot) == i[1]:
                return i[1]
            #print "i:  ", i[0]
            #print "i:  ", i[1]
       # print ""
    return "Could not find the path to potential: " + str(pot)


def set_path_db(path=os.environ['HOME'] + "/db"):
    return my.checkdir(path)


def set_path_db_utils(path=os.environ['HOME'] + "/db/utils"):
    """-pdu"""
    return my.checkdir(path)


##########################
##  avail
##########################


def avail_v(path=None):
    """
    -av
    """
    if path is None:
        path = str(my.pwd())

    string = str(path).split("/")
    if string[3] != "v":
        my.exit(error=" you are not in ~/v")
    return True


def avail_pots():
    """
    -ap lists all folder in set_path_vasp_potentials
        z.B. ['pl', 'US_GGA', 'PAW_GGA', 'US_LDA', 'pp']
        only folder are listed
    """
    #[ name for name in  os.listdir(set_path_vasp_potentials()) if os.path.isdir()
    out = []
    for name in os.listdir(set_path_vasp_potentials()):
        path = set_path_vasp_potentials() + "/" + name
        if os.path.isdir(path):
            out.append(name)
    return out
    #return filter(os.path.isdir, os.listdir(set_path_vasp_potentials()))


def avail_pots_elements(path=None):
    """
    -ape z.B. ['cu','cu_pv',... for current potential]
    """
    if path is None:
        path = my.pwd()

    pot = from_path_string_potential(path=path)
    #print "pot: ",pot
    path = set_path_vasp_potentials() + '/' + pot
    my.checkdir(path)
    import glob
    availpots = glob.glob(path + "/*")
    outavailpots = []
    for i in availpots:
        #print "i",i
        if os.path.isfile(i):
            continue
        if len(i.split("/")) <= 1:
            continue
        add = i.split("/")[-1]
        #print "add",add,type(add)
        outavailpots.append(add)
    #print outavailpots
    return outavailpots


def avail_jobs():
    """
    -aj z.b. ['elcp', 'eqvol', 'ti', 'murn', 'dynmat']
    """
    import glob
    list = glob.glob(set_path_db_utils() + "/job__*")
    #print list
    strip_list = [x.replace(set_path_db_utils() + "/job__", "") for x in list]
    out = strip_list
    #newgrades = {'eqcol': 'A', 'murn':'B'}
    #return sorted(out, key=newgrades.__getitem__)
    #    print out
    #    #if "murn" in out: print "yo"
    #    a, b = out.index('eqvol'), out.index('murn')
    #    print a,b
    #    print out
    #    print a,b
    #    print out
    #print out[b], out[a] = out[a], out[b]
    return out


def avail_cells():
    """
    -ac
    """
    import glob
    list = glob.glob(set_path_db_utils() + "/cell__*")
    #print list
    strip_list1 = [x.replace(set_path_db_utils() + "/cell__", "")
                   for x in list]
    #    print strip_list1
    strip_list = [x.split("=")[0] for x in strip_list1]
    #    for x in strip_list1:
    #        print x.split("=")[0]
    #    print strip_list
    return strip_list


# /nas/glensk/v/PAW_PBE/Al/murn_divak_fcc4/2x2x2sc_250eV-NGXF120-ADDGRID_6x6x6kP
# /nas/glensk/v/__potential_string__/__element_string__/__job_string__/__job_specification__/__job_specification__

#########################
## from_path_string
#########################
def from_path_string_potential(path=None):
    """
    -fpsp #gets the potential: paw_pbe, us_lda, ... from the path
    """
    #print "yo"
    if path is None:
        path = my.pwd()
    else:
        path = str(path)
    #print "||||||||", str(path)

    string = path.split("/")
    avail_v(path=path)

    if len(string) <= 4:
        my.exit("Which potential ?")
    return string[4]


def from_path_string_element(path=None):
    """
    -fpse #gets the element string: Al, Al_h, ... from the path
    """
    if path is None:
        path = my.pwd()
    else:
        path = str(path)
    #    print "||||||||",str(path)

    string = path.split("/")
    avail_v(path=path)

    if len(string) <= 5:
        my.exit(error="Which element ?")
    return string[5]


def from_path_string_element_list(path=None):
    """
    -fpsel
    """
    if path is None:
        path = my.pwd()
    else:
        path = str(path)


def from_path_string_project(path=None, exit=True):
    """
    -fpsjtc
    """
    ### if reference calculation:

    if path is None:
        path = my.pwd()
    else:
        path = str(path)
    #    print "||||||||",str(path)

    string = path.split("/")
    if len(string) <= 6:
        if exit is False:
            return None
        else:
            my.exit("Which project ? (job_type_cell)")
    #print "STRING:",string
    proj = string[6]  # murn_vak_fcc4
    #print "proj: ",proj,"        split: ",proj.split("_"),"        len:",len(proj.split("_"))
    if len(proj.split("_")) != 3:
        my.exit(error="project: \"" + str(proj) +
                "\" does not contain 3 words like murn_vak_fcc4")
    return proj


def from_path_string_subproject(path=None, exit=True):
    """
    -fpssp
    """
    return from_path_string_details_real(path=path, exit=exit)


def from_path_string_angvoltemp(path=None, exit=True):
    """
    -fpsa
    """
    if path is None:
        path = my.pwd()
    string = path.split("/")
    if len(string) <= 8:
        if exit is False:
            return False
        else:
            my.exit("Not deep enoug to bi in Ang/Vol/Temp")
    return string[8]


def from_path_string_job(path=None):
    """
    -fpsj   (z.B. murn, dynmat, ti, elcp, eqvol)
    """
    if path is None:
        path = my.pwd()
    return from_path_string_project(path=path).split("_")[0]


def from_path_string_job_is(job=None, path=None):
    """
    -fpsjit
    """
    if path is None:
        path = my.pwd()
    if job is None:
        my.exit(error=" specify jobtype to test: e.g. ti,dynmat, ...")
    check = str(from_path_string_job(path=path))
    if str(job) == check:
        return True
    else:
        return my.exit(error=" You are not in " + job + " folder")


def from_path_string_structure(path=None, exit=True, exit_if_not=None):
    """
    -fpss
    """
    if path is None:
        path = my.pwd()
    #print "in fps_jobtype: path:",path
    ### reference jobs are bulk always!
    jobtypeget = from_path_string_project(
        path=path, exit=exit)  # dynmat_vak_fcc4
    if jobtypeget is None and exit is False:
        return None
    if jobtypeget is None and exit is True:
        my.exit(error="No project found")

    #print "2",jobtypeget
    jobtype = jobtypeget.split("_")[1]
    #print "TT",jobtype
    try:
        alldetails = from_path_string_details(path=path)
    #        print "One:",alldetails
    except:
        alldetails = None
    #        print "two"
    #        print "twosss",alldetails

    #    print "jetzt:",alldetails

    if alldetails is None:
        out = jobtype
        if exit_if_not is None:
            return out
        else:
            if exit_if_not != out:
                my.exit(error="structure is not " + str(exit_if_not))

    if alldetails is not None:
        ref = from_path_get_isreferencejob(path=path)
        if ref is True:
    #        print "out:bulk"
            return "bulk"
        else:
    #        print "out: ",jobtype
            return jobtype


def from_path_string_structure_is(path=None, structure=None, structureisnot=None):
    """
    -fpsjsi
    """
    if path is None:
        path = my.pwd()
    if structure is None and structureisnot is None:
        my.exit(error=" specify structure to test")
    if structure is not None and structureisnot is not None:
        my.exit(error=" just define structure or structureisnot ; not both")

    if structure is not None and structureisnot is None:
        if str(from_path_string_structure(path=path)) == str(structure):
            return True
        else:
            return my.exit(error=" structure is not " + str(structure))

    if structure is None and structureisnot is not None:
        if str(from_path_string_structure(path=path)) != str(structureisnot):
            return True
        else:
            return my.exit(error=" structure acutally is " + str(structureisnot))


def from_path_string_cell(path=None):
    """
    -fpsc
    """
    if path is None:
        path = my.pwd()
    #print "from_path_string_cell: path: ",path
    proj = from_path_string_project(path=path)
    #print "from_path_string_cell: proj: ",proj
    return proj.split("_")[2]


def from_path_string_cell_is(path=None, cell=None):
    """
    -fpsjci
    """
    if path is None:
        path = my.pwd()
    if cell is None:
        my.exit(error="specify jobcell to test")
    if from_path_string_cell(path=path) == str(cell):
        return True
    else:
        return my.exit(error="Job Cell is not " + str(cell))


def from_path_string_details(path=None, justtest=False):
    """
    -fpsd
    """
    if path is None:
        path = my.pwd()
    else:
        path = str(path)

    string = path.split("/")
    #    print "string;",string
    avail_v(path=path)

    if len(string) <= 7:
        if justtest is True:
            return False
        else:
            my.exit(error="Your ar not deep enought for jobdetails e.g. 2x2x2sc_300eV_NGXF120_ADDGRID_16x16x16kP")
    outstring = string[7:]
    if len(outstring) == 1:
        return outstring[0]
    else:
        return "_".join(string[7:])


def from_path_string_details_real(path=None, exit=True):
    """
    -fpsdr
    """
    #print "path: ",path
    #print my.pwd()
    #quit()
    if path is None:
        path = my.pwd()
    else:
        path = str(path)
    #print "||||||||",str(path)
    #return "yo"
    string = path.split("/")
    avail_v(path=path)

    #print "STRING:",string
    #print "inpath:",path
    if len(string) <= 7:
        if exit is True:
            my.exit(error="Cutoff ? Kpoints ? Supercell ? --> Jobdetails missing!")
        if exit is False:
            return None
    return string[7]
    #return "/".join(string[7:])


def from_path_string_details_supercell(sc123=None, path=None):
    #def from_path_string_details_supercell(sc123=None,ka=1,kb=2):
    """
    -fpsds z.B. 2x3x5 sc123=1 --> 2; sc123=3 --> 5
    """
    #    print "BEGIN sc123:",sc123," type:",type(sc123)
    #    print "BEGIN sc123:",sc123," type:",type(sc123)," ka:",ka," kb:",kb
    from_path_string_details_real(path=path)  # just to ensure we have details
    if path is None:
        path_tocheck = my.pwd()
    else:
        #path_tocheck = str(from_path_string_details(path=path))
        path_tocheck = str(path)
    #        print "||||||||",str(path_tocheck)

    if sc123 is None:
        pass
    if sc123 is not None:
        if not my.is_int(sc123):
            my.exit("sc123 can only be 1 2 or 3; you defined sc123: " +
                    str(sc123))
        sc123 = int(sc123)
        if sc123 == 1 or 2 or 3:
            pass
    #        elif sc123 == 2:
    #            pass
    #        elif sc123 == 3:
    #            pass
        else:
            my.exit("sc123 can only be 1 2 or 3; you defined sc123: " +
                    str(sc123))

    get = "[0-9]*x[0-9]*x[0-9]*sc"
    import re
    out = []
    for word in path_tocheck.split("_")[:]:
        if re.search(get, word) is not None:
            add = re.search(get, word).group()
            out.insert(1, add)
    if len(out) is 1:
        if sc123 is None:
            return out[0]
        if sc123 is 1:
            return out[0].split("sc")[0].split("x")[0]
        if sc123 is 2:
            return out[0].split("sc")[0].split("x")[1]
        if sc123 is 3:
            return out[0].split("sc")[0].split("x")[2]
    else:
        my.exit('got ' + str(out) + " " + str(len(out)) +
                " number of sc from_path, I need exactlt one")


def from_path_string_details_supercell_ensure_sc123_isequal(path=None):
    """
    -fpsdsee
    """
    if path is None:
        path = my.pwd()
    else:
        path = str(path)

    sc1 = from_path_string_details_supercell(sc123="1", path=path)
    sc2 = from_path_string_details_supercell(sc123="2", path=path)
    sc3 = from_path_string_details_supercell(sc123="3", path=path)

    if sc1 == sc2:
        if sc2 == sc3:
            return True
    else:
        my.exit(error="supercell is not x^3 but has different dimensions")


def from_path_string_details_kpoints(path=None):
    """
    -fpsdk
    """
    if path is None:
        path = my.pwd()
    else:
        path = str(path)

    from_path_string_details_real(path=path)  # just to ensure we have details

    if path is None:
        path_tocheck = my.pwd()
    else:
        #path_tocheck = str(from_path_string_details())
        path_tocheck = str(path)
    #    print "||||||||",str(path)

    get = "[0-9]*x[0-9]*x[0-9]*[k][p]"  # see later ignorecase
    import re
    out = []
    for word in path_tocheck.split("_")[:]:
    #        print "word:",word
        string = re.search(get, word, re.IGNORECASE)

        if string is not None:  # er hat irgendwas mit ...kp gefunden
    #            print "get:",get
    #            print "word:",word
    #            print "str;",str
            add = string.group()
    #            print "add:",add
            out.insert(1, add)
    #        print "out;",out
    #
    #    print path.split("_")
    if len(out) is 0:
        my.exit(error=" couldnt find number of kpoints from_path")

    kp = out[-1]
    kp = kp.lower().strip("kp").split("x")

    if len(kp) != 3:
        my.exit("from_path kapoints " + kp + " not length of 3")
    kp[0] = re.sub("^0+", "", kp[0])
    kp[1] = re.sub("^0+", "", kp[1])
    kp[2] = re.sub("^0+", "", kp[2])
    #    print kp[0],kpn
    return kp


def from_path_string_details_cutoff(path=None):
    """
    -fpsdc
    """
    from_path_string_details_real(path=path)  # just to ensure we have details
    if path is None:
        path_tocheck = my.pwd()
    else:
        #path_tocheck = str(from_path_string_details(path=path))
        path_tocheck = str(path)
    #    print "||||||||",str(path)
    get = "[0-9]*eV"
    import re
    out = []
    for word in path_tocheck.split("_")[:]:
        if re.search(get, word) is not None:
            add = re.search(get, word).group()
            out.insert(1, add)
    if len(out) != 0:
        cutoff = out[-1].replace('eV', '')
    else:
        my.exit("cutoff not in path")
        #my.exit()
    if my.is_number(cutoff) is not True:
        sys.exit("cutoff " + str(cutoff) + ' is no number')
    else:
        return cutoff


def from_path_string_details_ediff(path=None):
    """
    -fpsde
    """
    from_path_string_details_real(path=path)  # just to ensure we have details
    ediff = None

    ###########################
    ###########################
    ## if EDIFF is defined in path
    ############################
    ############################
    if path is None:
        path_tocheck = my.pwd()
    else:
        #path_tocheck = str(from_path_string_details(path=path))
        path_tocheck = str(path)
    #    print "||||||||",str(path)
    get = "ediFF[0-9]*E-[0-9]*"
    import re
    out = []
    for word in path_tocheck.split("_")[:]:
        if re.search(get, word, re.IGNORECASE) is not None:
            add = re.search(get, word, re.IGNORECASE).group()
            out.insert(1, add)
    #    print "ED",out
    if len(out) != 0:
        ediff = out[-1].lower().replace('ediff', '').upper()
        ediff = str(ediff)

    if ediff is not None:
        return ediff

    ###########################
    ###########################
    ## else standardwerte
    ############################
    ############################
    ## if dynmat: 1E-8
    job = from_path_string_job(path=path)
    #    print "job;",job,type(job)
    if job == "dynmat":
        return "1E-8"

    ## if murn: 1E-6
    #    print "job;",job,type(job)
    if job == "murn":
        return "1E-6"

    ## if reference:
    ## ...

    ### in case nothing got
    return my.exit(error=" EDIFF not evalueated")


def from_path_string_details_ismear(path=None):
    """
    -fpsdi
    """
    from_path_string_details_real(path=path)  # just to ensure we have details
    ismear = None

    ###########################
    ###########################
    ## if ISMEAR is defined in path
    ############################
    ############################
    if path is None:
        path_tocheck = my.pwd()
    else:
        #path_tocheck = str(from_path_string_details(path=path))
        path_tocheck = str(path)
    #    print "||||||||",str(path)
    get = "IsMeAr[-0-9]*"
    import re
    out = []
    for word in path_tocheck.split("_")[:]:
        if re.search(get, word, re.IGNORECASE) is not None:
            add = re.search(get, word, re.IGNORECASE).group()
            out.insert(1, add)
        #    print "ED",out
    if len(out) != 0:
        ismear = out[-1].lower().replace('ismear', '').upper()
        ismear = str(ismear)

    if ismear is not None:
        return ismear

    ###########################
    ###########################
    ## else standardwerte
    ############################
    ############################
    ## 1,2 MethPaxton,-1 Fermi
    job = from_path_string_job(path=path)

    if job == "dynmat":
        return "-1"
    if job == "murn":
        return "1"  # but not for semiconcuctors TODO: Check if semiconductor or insulator
    if job == "ti":
        return "-1"  # this may be changed TODO: check ISMEAR for TI
    if job == "elcp":
        return "-1"
    ## if reference:
    ## ...

    ### in case nothing got
    return my.exit(error=" ISMEAR not evalueated")


def from_path_string_details_ngxf(path=None):
    """
    -fpsdn
    """
    ngxf = None

    ###########################
    ###########################
    ## if ngxf is defined in path
    ############################
    ############################
    if path is None:
        path = my.pwd()
    else:
        #path = str(from_path_string_details())
        path = str(path)
    #    print "||||||||",str(path)
    get = "ngXF[0-9]*"
    import re
    out = []
    for word in path.split("_")[:]:
        if re.search(get, word, re.IGNORECASE) is not None:
            add = re.search(get, word, re.IGNORECASE).group()
            out.insert(1, add)
        #    print "ED",out
    if len(out) != 0:
        ngxf = out[-1].lower().replace('ngxf', '').upper()
        ngxf = str(ngxf)

    if ngxf is not None:
        return ngxf
    else:
        return None


def from_path_get_ang(path=None, exit=True):
    """
    -fpga
    """
    if path is None:
        path = my.pwd()
    else:
        path = str(path)
    #    print "||||||||",str(path)
    angs = my.grep("[0-9]*.[0-9]*[Aa][Nn][Gg]", str(path), "-v")
    angs2 = my.grep("OUTCAR.[0-9]*.[0-9]*", str(path), "-v")
    #    print "angs:",angs
    #    print "angs2;",angs2,angs2[0]
    if len(angs) == 1:
        return my.grep("[0-9]*.[0-9]*", angs[0], "-v")[0]
    elif len(angs2) == 1:
    #        print ""
        alat = angs2[0][7:]
        return alat

    elif len(angs) is 0:
        if exit is True:
            my.exit(error=" it rather seems you are not in ang folder _____; path: " +
                    str(path))
        else:
            return None

    else:
        #return "Several Ang folder found:",angs
        if exit is True:
            my.exit(error=" Several Ang folder found:" + str(angs))
        else:
            return None
    #    #if re.search("[0-9]*ang",my.pwd()) != None:


##########################
##  from_path_get
##########################
def from_path_get_isreferencejob(path=None):
    """
    -fpsiref
    """

    if path is None:
        path = my.pwd()
    else:
        path = str(path)
    details = from_path_string_details_real(path=path)

    one = len(my.grep("REF[_]", str(details), options='-i'))
    two = len(my.grep("[_]REF", str(details), options='-i'))
    #    print "o",one,two
    if one or two is 1:
    #        print "out:bulk"
        return True

    #print "weiter"
    lookforref = len(my.grep("referenCe", details, "-i"))
    #print "L",lookforref
    if lookforref >= 1:
        return True

    ## if no "reference" found print False
    return False


def from_path_get_isrelaxjob(path=None):
    """
    -fpsirel
    """
    if path is None:
        path = my.pwd()
    else:
        path = str(path)

    from_path_string_details_real(path=path)  # just to ensure we have details
    job = from_path_string_job(path=path)
    if job == "ti":
        return False
    if job == "dynmat":
        return False
    if job == "elcp":
        return False

    if job is "murn":  # grundsaetzlich ja, ausser es ist eine Reference --> get reference
        ref = from_path_get_isreferencejob(path=path)
        if ref is True:
            return False

        if ref is False:
            return True

    return my.exit(error=" is this a relaxjob?")


def from_path_get_numatoms(path=None):
    """
    -fpgn
    """
    #print "pathin:",path

    if path is None:
        path = get_path_positions_direct(path=my.pwd())
    else:
        path = get_path_positions_direct(path=path)
    #print "pathdanach:",path
    file = my.readfile(path)
    while '' in file:
        file.remove('')
    #    print file
    #    print "one",file[1],"len:",len(file[1])
    return len(file)


def from_path_get_nbands_thisstruct(path=None):
    """
    -fpgbands
    """
    if path is None:
        path = my.pwd()
    else:
        path = path

    atoms_now = from_path_get_numatoms(path=path)
    #print "atoms_now:",atoms_now

    ################### searching for NBANDS files ##############
    NBANDSfile = False
    ################### SUCHVORGANG #############################
    ## 0. look for NBANDS: in current folder
    ## 1. look for NBANDS: within subproject (you have to be at leas within subproject or deeper)
    ## 2. look for NBANDS: within project --> get all NBANDS in all subpfojectsfolder
    ## 3. look for NBANDS: within equal job (murn_divak, murn_vak, murn_bulk, ...)

    ## a) get NBANDS list [[30, 50], [30, 52], [30, 53]]
    ## b) get atoms list [30, 32, 106, 108]

    NBANDSfile = False
    #print "YO"
    ## 0. look for NBANDS: in current folder
    if my.isfile(str(path) + "/NBANDS") is True:
        NBANDSfile = [my.checkfile(str(path) + "/NBANDS")]
        #print "-> NBANDS file(s) within (currentfolder) found! "+str(path)

    else:
        #print "YAST"
    ## 1. look for NBANDS: within subproject (you have to be at leas within subproject or deeper)
        #print "NO NBANDS file(s) within (currentfolder): "+str(path)

        subprojpath = get_path_subproject(path=path, Error_if_notexist=False)  # what happens here if reference folder is given as folder which does
        ## not  exisst? --> egal erstmal, subprojpath
        #print "YY",subprojpath
    #        quit()
        if my.isfile(str(subprojpath) + "/NBANDS") is True:
            NBANDSfile = [my.checkfile(str(subprojpath) + "/NBANDS")]
            #print "-> NBANDS file(s) within ( subproject  ) found! "+str(subprojpath)

        else:
        ## 2. look for NBANDS: within project --> get all NBANDS in all subpfojectsfolder
            #print "NO NBANDS file(s) within ( subproject  ): "+str(subprojpath)
            projpath = get_path_project(path=path)
            find = my.run("find " + str(projpath) + " -maxdepth 2 -mindepth 2 -type f -name NBANDS | xargs").split()
            if len(find) > 0:
                NBANDSfile = find
                #print "-> NBANDS file(s) within (   project   ) found! : "+str(projpath)

            else:
            ## 3. look for NBANDS: within equal job (murn_divak, murn_vak, murn_bulk, ...)
                #print "NO NBANDS file(s) within (   project   ) "+str(projpath)
                print "--> TODO: here we could thin of trying other projects first:"
                electrons = potcar_electrons(
                    path=path)  # path = None: go to standard POTCAR
                print "!!!!!!!!! ESTIMATING NBANDS: electrons per Atom in POSCAR: ", electrons, " atoms: ", atoms_now, " !!!!!!!!"
                return my.iroundup(electrons * atoms_now * 1.37 * 0.5)  # *0.5 every band is occupied by two electrons
    #print "KKK"
    ## a) get NBANDS list [[30, 50], [30, 52], [30, 53]]
    allNBANDSfiles = " ".join(NBANDSfile)
    cat = my.run("cat " + allNBANDSfiles + " | sort -n | uniq").split("\n")
    list = []
    for i in cat:
        if len(i) == 0:
            continue
        line = i.split()
        if len(line) != 2:
            continue
        atoms = line[0]
        bands = line[1]
        if my.is_int(atoms) is not True:
            continue
        if my.is_int(bands) is not True:
            continue
        list.append([int(atoms), int(bands)])
        #print "i",atoms,bands,my.is_int(atoms)
    #print "list:    ",list

    ## b) get atoms list [30, 32, 106, 108]
    import numpy as np
    foo = np.array(list)

    atomlist = sorted(set(foo[:, 0]))

    #print "atomlist:",atomlist
    #print "atoms_now:",atoms_now
    smaller_structures_avail = False
    smaller_structures_list = []
    for i in atomlist:
        if i < atoms_now:
            smaller_structures_avail = True
            smaller_structures_list.append(i)
    #print np.array(atomlist),np.array(atomlist)-np.array([atoms_now])

    if atoms_now in atomlist:  # exact number of atoms
        atom_cat = atoms_now
        sicherheit = 0.01
    else:  # second option are smaller structures with less atoms
        if smaller_structures_avail is True:
            atom_cat = max(smaller_structures_list)
            sicherheit = 0.01
        else:
            atom_cat = min(atomlist)
            sicherheit = 0.04

    #print "atom_cat:",atom_cat
    faktor_bands = float(atoms_now) / float(atom_cat)
    #print "faktor:",faktor_bands
    #bandlist = sorted(set(foo[:,1]))
    bands_out = max([y for x, y in list if x == atom_cat])
    sicherheit_out = my.iroundup(bands_out * sicherheit * faktor_bands)
    #print "bands_out:",bands_out
    bands_min = my.iroundup(float(bands_out) * float(faktor_bands))
    #print "bands_out_mal_faktor",
    #print "sicherheit:",sicherheit_out
    return bands_min + sicherheit_out

    ######### WENN SONST NIX GEFUNDNE SUCHE IN /home/glensk/db/nbands_occupied.dat
    bands_peratom = db.nbands_peratom(path=path)
    #print "bands_peratom:",bands_peratom
    if my.is_number(bands_peratom) is not True:
        my.exit(error=" NBANDS not found")
    #    print "jnum:",my.is_number(bands_peratom)
    #print "bands_peratom:",bands_peratom
    bands = int(atoms) * float(bands_peratom)
    sicherheit = bands * 0.1
    import math
    retint = int(math.ceil((bands + sicherheit)))
    if my.is_int(retint) is not True:
        my.exit(error=" nbands problem in calculation, this error is easy to solve!")
    return retint


def from_path_get_kpoints_product(path=None):
    """
    -fpgkp
    """
    if path is None:
        path = my.pwd()
    else:
        path = str(path)
    kp = from_path_string_details_kpoints(path=path)
    return int(kp[0]) * int(kp[1]) * int(kp[2])


def from_path_get_kpoints_mal_atoms(path=None):
    """
    -fpgkma
    """
    if path is None:
        path = my.pwd()
    else:
        path = str(path)
    atoms = from_path_get_numatoms(path=path)
    kp = from_path_string_details_kpoints(path=path)
    return atoms * int(kp[0]) * int(kp[1]) * int(kp[2])


def from_path_get_ALATS(path=None, string=True):
    """
    -fpgA
    """
    if path is None:
        path = my.pwd()
    else:
        path = str(path)
    #print "path:",path
    pathalat = None

    try:
        #print "10"
        pathalat = my.checkfile(str(path) + "/ALATS")
        #print "11",pathalat,my.checkfile(pathalat)
    except:
        try:
            #print "20"
            pathalat = my.checkfile(str(path) + "/../ALATS")
            #print "21"
        except:
            try:
                #print "30"
                #from_path_get_ang(path=path)   ### just to chekc if we have an ang, in that case we can go further down
                #print "31"
                pathalat = my.checkfile(str(path) + "/../../ALATS")
            except:
                #print "30"
                #from_path_get_ang(path=path)   ### just to chekc if we have an ang, in that case we can go further down
                #print "31"
                pathalat = my.checkfile(str(path) + "/../../../ALATS")

    alats = my.readfile(pathalat)[0].split()
    if string == True:
        return alats
    else:
        floats = map(float, alats)
        return floats


def from_path_get_ALATS_VOLS_folderlist_to_create_for_subproj_or_angvoltemp(path=None, rise_error_if_exists=False, silent=True):
    """
    -gf looks in current folder and subfolder for ALATS file and prints out jobfolder to be created
    """
    #print "path:ganz:",path
    if path is None:
        path = my.pwd()
    else:
        path = str(path)

    from_path_string_details(path=path)  # just to ensure we are in details
    #    print "danach:"
    #    print ""
    try:
        ang = from_path_get_ang(path=path)  # if we are in ang folder we dont want to create all! angfolders but just this
    except:  # catch *all* exceptions
        ang = None

        #print "ang:",ang
    #    print ""
    if ang is None:
        alats = from_path_get_ALATS(string=True)
    else:
        alats = [str(ang)]
        #print "alats:",alats
    #    print "riseerror:?",rise_error_if_exists

    alatsFolder = [path + "/" + x + "Ang" for x in alats]
    #print "alatsFolder:",alatsFolder,len(alats),str(path.split("/")[-1]),str(alats[0])+"Ang"
    ########################################
    ### check if path has already x.xxAng/
    #print "YYY",path.split("/")[-1],len(alats)

    if len(alats) == 1 and str(path.split("/")[-1]) == str(alats[0]) + "Ang":
        alatsFolder = [path for x in alats]
        #print "alatsFolder:",alatsFolder
    #print "yes"
    outfolder = []
    for folder in alatsFolder:
    #        print "folder:",folder
        if my.isdir(folder) == True and my.pwd() != folder:  # WENN er existier, nicht hinzufuegen, ausser wenn man gerade in dem folder ist!
            if silent == False:
                print folder, " already done"
            if rise_error_if_exists == True:
                folderout = "/".join(folder.split("/")[-2:])
                #               print "F>>>>>>",folderout
                return my.exit(error=" Folder already exists1: .../" + folderout)
        else:  # ANSONSTEN hinzufuegen
            if silent == False:
                print folder
            outfolder.append(folder)
            #    p
            # rint "OUT:",outfolder
        #print "yo"
    return outfolder

###########################
## from_path getpath
###########################


def get_path_potential(path=None):
    """
    -gppot   prints part of path up to e.g. ~/v/pp
    """
    if path is None:
        path = my.pwd()
    else:
        path = str(path)

    avail_v(path=path)  # just to ensure we are in ~/v
    #print "path:",len(path.split("/"))
    if len(path.split("/")) <= 4:
        my.exit(error="You are not deep enough to be in a potential")
    potpath = "/".join(path.split("/")[0:5])
    ## ensure you are in path w
    #print elementpath
    return my.checkdir(potpath)


def get_path_element(path=None):
    """
    -gpe
    """
    if path is None:
        path = my.pwd()
    else:
        path = str(path)
    avail_v(path=path)  # just to ensure we are in ~/v
    #print "path:",len(path.split("/"))
    if len(path.split("/")) <= 5:
        my.exit(error="You are not deep enough to be in an element")
    elementpath = "/".join(path.split("/")[0:6])
    ## ensure you are in path w
    #print elementpath
    return my.checkdir(elementpath)


def get_path_project(path=None, just_tell_if_exists=False):
    """
    -gpp
    """
    if path is None:
        path = my.pwd()
    else:
        path = str(path)

    words = len(path.split("/"))
    #print "words:",words, "words:",path.split("/")
    if just_tell_if_exists == True:
        if words < 7:
            return False
        else:
            return True
    projectpath = "/".join(path.split("/")[0:7])
    from_path_string_project(path=path)  # just to ensure we are deep enough
    return my.checkdir(projectpath)


def get_path_subproject(path=None, Error_if_notexist=True):
    """
    -gps  ### will ich hier check ob der existiert oder nicht?
    """
    if path is None:
        path = my.pwd()
    else:
        path = str(path)

    job = from_path_string_job(path=path)
    if job == "eqvol":
        return get_path_project(path=path)
    #print "job:",job

    #print "aaa,bbb"
    from_path_string_details(
        path=path)  # just to ensure we have details(=subproject)
    #print  "ccc"
    projectpath = "/".join(path.split("/")[0:8])
    #print "ddd"
    #print "PP:",projectpath
    #print my.checkdir(projectpath)
    #print "Ka"
    if Error_if_notexist == False:
        return projectpath
    #        try:
    #            if my.isdir(projectpath) == False:
    #                return False
    #            my.checkdir(projectpath)
    #        except:
    #            return False

    return my.checkdir(projectpath)


def get_path_angvoltempfolder(path=None, Error_if_notexist=True):
    """
    -gpa
    """
    if path is None:
        path = my.pwd()
    from_path_string_angvoltemp(path=path, exit=Error_if_notexist)
    out = "/".join(path.split("/")[0:9])
    return out


def get_path_all_elements(list=True, path=None):
    """
    -gpae
    """
    if path is None:
        path = str(my.pwd())

    avail_v(path=path)  # just to ensure we have path
    from_path_string_potential(path=path)  # ust to ensure we have a potential
    potpath = get_path_potential(path=path)

    return potpath


def get_path_all_projects(list=True, path=None):
    """
    -gpap
    """
    #print "path:",path
    #print "list:",list
    if path is None:
        path = str(my.pwd())

    avail_v(path=path)  # just to ensure we have path
    from_path_string_potential(path=path)  # ust to ensure we have a potential
    from_path_string_element(path=path)  # ust to ensure we have an element

    av_jobs = avail_jobs()
    elepath = get_path_element(path=path)
    allproj = []
    import glob
    for job in av_jobs:
        #print "job:",job
        string = elepath + "/" + str(job) + "_*"
        #print "string:",string
    #    print "string:",string
        add = glob.glob(string)
        if len(add) == 0:
            continue
        if len(add) != 0:
            allproj.extend(add)
        #print "add",add
        #print ""
    #print 'vor',list
    if str(list) == "True":
        #print "in1"
        return ("".join([i + "\n" for i in allproj]))[0:-1]
    else:
        #print "in2"
        return allproj


def get_path_all_subprojects(path=None):
    """
    -gpas
    """
    if path is None:
        path = str(my.pwd())
    #print "path:",path

    if check_path_is_projectpath(path=path) != True:
        my.exit(error="to get all subprojects you have to bin in project....erstmal")

    if from_path_string_job(path=path) == "eqvol":
        return [path]
    #print ":",path
    #import os
    #import os.path

    all_outfolder = []
    for x in os.listdir(path):
        if os.path.isfile(x) == True:
            continue
        else:
            #if len(my.grep("x[0-9]*sc",x,options="-i")) <= 0: continue
            if len(my.grep("x[0-9]*kp", x, options="-i")) <= 0:
                continue
            if len(my.grep("[0-9]*ev", x, options="-i")) <= 0:
                continue
            #print "x -> ",x,my.grep("sc",x)
            subproj = str(path) + "/" + str(x)
            if os.path.isfile(subproj) == True:
                continue
            all_outfolder.append(subproj)
            #print subproj
    #for xx in all_outfolder:
    #    print "gpa",xx,os.path.isfile(xx)
    #print ""
    #for i in all_outfolder:
    #    print i                        ### some definitions need the [1, 2, 3]

    return all_outfolder


def get_path_all_angvoltempfolder(path=None):
    """
    -gpaavt
    """
    if path is None:
        path = str(my.pwd())

    subprojpath = get_path_subproject(path=path)
    folder = my.ListFolders(subprojpath)
    #print "yo",subprojpath
    folderout = []
    for i in folder:
        check = i.split("/")[-1]
        #print "check"
        c1 = check, my.grep("AnG", check, options="-i")
        c2 = check, my.grep("K", check, options="-i")
        c3 = check, my.grep("VOl", check, options="-i")
        if c1[-1] == [] and c2[-1] == [] and c3[-1] == []:
            continue  # none of Ang, K Vol grept
        folderout.append(i)
        #print check,"                          00",c1,len(c1)
        #print len(c1),len(c2),len(c3)
        #print len(c1),len(c2),len(c3)
    return folderout


def get_path_subproject_giving_projectpath_cutoff_kpoints(projectpath=None, cutoff=None, kpoints=None, exclude_referencefolder=False):
    """
    -gpsf
    """
    if projectpath is None or cutoff is None or kpoints is None:
        my.exit(error="please provied projectpath, cutoff and kpoints")
    import glob
    cutoff_list = glob.glob(projectpath + "/*" + str(cutoff) + "[eEVv]*")
    kpoints_list = glob.glob(projectpath + "/*" + str(kpoints[0]) + "[xX]*" + "*" + str(kpoints[0]) + "[xX]*" + "*" + str(kpoints[0]) + "[kKpP]*")
    #for i in cutoff_list:
    #    print i
    #print "_)"*40
    #for i in kpoints_list:
    #    print i
    #print "||||||"*40
    schnittmenge = my.schnittmenge(cutoff_list, kpoints_list)
    #for i in schnittmenge:
    #    print i
    #print "++"*30
    reference_list = []
    if exclude_referencefolder == True:
        reference_list = glob.glob(projectpath + "/*" + "_[rR][eE][fF]*")
        #for i in reference_list:
        #    print i
    #print "--"*30
    out = my.differenzmenge(schnittmenge, reference_list)
    if len(out) != 1:
        print ""
        for i in out:
            print "--> ", i

        my.exit(error="more thann one subproject with cutoff and kpoints in projectpath found")
    return out[0]


def get_path_db_utils_cell(path=None):
    """
    -pduc
    """
    #    from_pathcell = from_path_string_jobcell()
    #    availcell = avail_cells()
    if path is None:
        path = my.pwd()
    else:
        path = str(path)

    import glob
    string = set_path_db_utils(
    ) + "/cell__" + from_path_string_cell(path=path) + "*"
    #    print "string:",string
    pathdb = glob.glob(string)
    #    print "path:",path
    if len(pathdb) != 1:
        my.exit("couldnd find path to cell: " + string)
    pathdb = pathdb[0]
    #    print "yo"
    my.checkdir(pathdb)
    #    print "len:",len(path)
    #    print "list:",list
    #    print "fp:",from_pathcell
    #    print "av:",availcell
    return pathdb


def get_path_db_utils_cell_cellfile(path=None):
    """
    -pducc
    """
    if path is None:
        path = my.pwd()
    pathout = my.checkfile(get_path_db_utils_cell(path=path) + "/CELL")
    return pathout


def get_path_vasp_potential(path=None):
    """
    -gpvp
    """
    if path is None:
        path = my.pwd()

    pot = from_path_string_potential(path=path)
    #print "pot: ", pot
    potnew = get_path_vasp_potentials_name(pot)
    #print "potnew: ", potnew
    element = from_path_string_element(path=path)
    #print "element:", element
    pathcheck = set_path_vasp_potentials() + '/' + potnew + '/'
    my.checkdir(pathcheck)
    pathout = set_path_vasp_potentials() + '/' + potnew + '/' + element

    # if we got the correct element "string"
    if os.path.isdir(pathout):
        return pathout
    thedir = pathcheck

    # if not, check for lowercae / uppercase
    elements = [name for name in os.listdir(thedir) if os.path.isdir(os.path.join(thedir, name))]
    #print elements
    #print "-------------"
    #print "el",element,"|", type(element)
    #if "al" == str(element):   # to check for equality of string, == is necessary and "is" will not work
    #    print "YES"
    for i in elements:
        elementavail = i
        elementlower = i.lower()
        if str(element) == str(elementlower):
            #print "done",element, elementavail
            pathout = set_path_vasp_potentials() + '/' + potnew + '/' + elementavail
            my.checkdir(pathout)
            return pathout
    my.exit("could not find element" + str(element) + " in " + str(pathcheck))


def get_path_vasp_potential_POTCAR(path=None):
    """
    -pp
    """
    if path is None:
        path = my.pwd()

    pathelem = get_path_vasp_potential(path=path)
    #print "path:",path
    my.checkdir(pathelem)
    pathpot = pathelem + '/' + 'POTCAR'
    #print "pathpot:",pathpot
    my.checkfile(pathpot)
    return pathpot


def potcar_electrons(path=None):
    """
    -pote
    """
    if path is None:
        path = get_path_vasp_potential_POTCAR(path=path)
    else:
        check = path.split("/")[-1]
        if check != "POTCAR":
            path = get_path_vasp_potential_POTCAR(path=path)
        else:
            path = path

    #print "ok"
    electrons = my.readfile(path=path, line=2)
    #print electrons
    if my.is_number(electrons) != True:
        my.exit(error="how many electrons does " + str(path) + " have?")
    return float(electrons)


def get_path_positions_direct(path=None):
    """
    -gppd
    """
    if path is None:
        path = my.pwd()
    #    print "pathinn:",path
    cell = get_path_db_utils_cell(path=path)
    #    print "ok1"
    sc = from_path_string_details_supercell(path=path)
    #    print "ok2"
    job = from_path_string_structure(path=path)
    #    print "path_db_utils_cell:",cell, "sc:",sc,"job:",job
    file = cell + "/" + job + "__" + sc + "__EqCoords_direct"
    #    print "file:",file
    pathout = my.checkfile(file)
     #+"/"+from_path_string_structure()+"__*"+from_path_string_details_supercell()+"__EqCoords_direct"
    #print 'string:',string
    #import glob
    #path = glob.glob(get_path_db_utils_cell()+"/"+from_path_string_structure()+"__*"+from_path_string_details_supercell()+"__EqCoords_direct")
    #    print "path:",path
    #    if len(path) != 1:
    #        print "couldnd find path to cellLLL: ",path
    #        print "1",get_path_db_utils_cell()
    #        print "2","/"+from_path_string_structure()
    #        print "3",from_path_string_details_supercell()
    #        my.exit()
    #    path = path[0]
    #    my.checkfile(path)
    return pathout

#def get_path_positions_cartesian_dynmat_ausgelenkt(coords=None,CELL=None):
#    """
#    -gppcd
#    """
#    return True


def get_path_job_type_cell(path=None):
    """
    -gpj
    """
    if path is None:
        path = my.pwd()
    else:
        path = str(path)

    pot = from_path_string_potential(path=path)
    ele = from_path_string_element(path=path)
    jobtypecell = from_path_string_project(path=path)
    #sc = from_path_string_details_supercell()
    pathout = my.checkdir(
        "/home/glensk/v/" + pot + "/" + ele + "/" + jobtypecell)
    return pathout


###########################
## getfile
###########################
def create_POTCAR(test=False, path=None):
    """
    -cpot
    """
    if path is None:
        path = my.pwd()
    else:
        path = str(path)
    pot = get_path_vasp_potential_POTCAR(path=path)
    import shutil
    if str(test) != "False":
        return True

    ### WRITE FILE
    shutil.copy(pot, "POTCAR")
    return True


def create_KPOINTS(test=False, path=None):
    """
    -ck
    """
    if path is None:
        path = my.pwd()
    else:
        path = str(path)

    #    print "ys",path
    kp = from_path_string_details_kpoints(
        path=path)  # .lower().strip("kp").split("x")
    #    print "kp:",kp
    #    if len(kp) != 3:
    #        print "from_path kapoints ",kp,"not length of 3" #.replace("x"," ")
    file = "KPOINTS"
    #    print "kp>",kp
    for ka in kp:
    #        print "ka",ka
    #        print "isint:",my.is_int2(ka)
        if my.is_int(ka) != True:
            my.exit("kpgot: " + kp + " not integer")
            #print "ka->",ka,"ty:",type(ka),"is?",my.is_int(ka)
    if os.path.isfile(file) == True:
        os.remove(file)
    if str(test) != "False":
        return True

    ### WRITE FILE
    filew = open(file, "w")
    filew.write("K-Points\n")
    filew.write(" 0\n")
    filew.write("Monkhorst Pack\n")
    filew.write(kp[0] + " " + kp[1] + " " + kp[2] + "\n")
    filew.write("0" + " " + "0" + " " + "0" + "\n")
    filew.close()

    return True


def create_POSCAR(path=None, test=False, fromrelaxedcoords=False, displace_atom=None, displace_direction=None, displacement_angstrom=None, direct_position_file=None, printoutput=False):
    """
    -cpos
    """
    if path is None:
        path = my.pwd()
    else:
        path = str(path)
    file = "POSCAR"
    #    print "in,path",path
    numatoms = from_path_get_numatoms(path=path)
    #    print "numat",numatoms
    #    print "in getfilePOS:",path
    cell = get_cell(path=path)
    #print "cell:",cell ##  [[8.24, 0, 0], [0, 8.24, 0], [0, 0, 8.24]]
    if printoutput:
        print "in getfilePOS:", cell

    ####################################################################################################
    ### what jobtype? ##                                  (z.B. murn, dynmat, ti, elcp, eqvol)
    ####################################################################################################
    #####################################################################################################
    ### What structure: ### bulk, vak, divak, subst, interst, ...
    #####################################################################################################
    jobtype = from_path_string_job(
        path=path)  # (z.B. murn, dynmat, ti, elcp, eqvol)
    structure = str(from_path_string_structure(
        path=path))  # bulk, vak, divak, subst, interst, ...
    if printoutput:
        print "jobtype::", jobtype         # murn
    if printoutput:
        print "structure::", structure     # vak

    ##########################
    ### dynmat --> with displacement fcc4
    ##########################
    if jobtype == "dynmat":
        if structure == "bulk":
            #print "xxxxxx"
            cartes_pos = get_positions_cartesian(path=path, displace_1_1=True)
            #print 'yyyyyyyy'
        else:
            print "get ref coords"
            if direct_position_file is not None:
                dir_pos = my.readfile(path=direct_position_file)
                cartes_pos = positions_reduced_to_cartesian(
                    path_coords=direct_position_file)
                print cartes_pos
    #        print "cart_pos:",cartes_pos
        string = "Cartes"
        pos = cartes_pos

    ##########################
    ### only for murn: (due to extrapolation of atomic positions for murn curve)
    ##########################
    ##if jobtype == "murn" or "vak" or "divak":  ## funktioniert so!
    if jobtype == "murn":
        #print "am in list"
        #quit()
        if structure == "bulk":
            dir_pos = get_positions_direct(path=path)
            #    print "dir_pos:",dir_pos
            string = "Direct"
            pos = dir_pos

        elif structure == "vak" or "divak":
                ### check if this is firststructure to relax -> take bulk coords
            grep = my.grep("relxfrombulk", my.pwd(), "-v")
            if len(grep) == 1:
                dir_pos = get_positions_direct(
                    path=path)  # falsch! gibt die ausgangasstruktur der vak
                                                            #### nicht die von bulk!
                string = "Direct"
                pos = dir_pos

            else:
                ### else take relaxedcoords of refpath ->
                #print "hall1"
                dir_pos_path = murn.get_path_to_relaxedcoords_direct_angvoltempfolder(thisjob=path)
                #print "hall2"
                dir_pos = my.readfile(dir_pos_path)
                #print "relaxed coords from:",dir_pos_path
                string = "Direct"
                pos = dir_pos
        else:
            print "need to get relaxed coords"
            my.exit(error=" Positions only known it this is bulk, not for divak,vak,...GOT:" + str(from_path_string_structure(path=path)))

    #######################################
    ## check all positions
    #######################################
    #print "JJJJJJJJJ"
    #print "DD",dir_pos
    #print "len(pos)::",len(pos)
    if numatoms != len(pos):
        my.exit(error="This structure has: " + str(numatoms) + " atoms; From the positions given I have only " + str(len(dir_pos)))
    for atom_pos in pos:
        #print "atom_pos:",atom_pos
        if len(atom_pos.split()) != 3:
        #if len(atom_pos) != 3:
            my.exit(error="position: " + str(atom_pos) + " does not have 3 coordinates in x,y,z but" + str(len(atom_pos)))

    if str(test) != "False":
        return True

    ###################################
    ### WRITE NEW FILE
    ###################################
    if os.path.isfile(file) == True:
        os.remove(file)
    filew = open(file, "w")
    filew.write("created by pythonscript\n")
    filew.write("1.0\n")
    filew.write(str(
        cell[0][0]) + "  " + str(cell[0][1]) + "  " + str(cell[0][2]) + "\n")
    filew.write(str(
        cell[1][0]) + "  " + str(cell[1][1]) + "  " + str(cell[1][2]) + "\n")
    filew.write(str(
        cell[2][0]) + "  " + str(cell[2][1]) + "  " + str(cell[2][2]) + "\n")
    filew.write(
        str(numatoms) + "\n")  # TODO: create POSCAR for more than one species
    filew.write(string + "\n")

    #################################
    ## here we go over the positions (pos)
    #################################
    for atom_pos in pos:
    #        print atom_direct
        filew.write(str(atom_pos) + "\n")
    filew.close()
    return True


def create_INCAR(path=None, test=False, printoutput=False):
    """
    -ci
    """
    if path is None:
        path = my.pwd()
    else:
        path = str(path)
    cutoff = from_path_string_details_cutoff(path=path)
    if printoutput:
        print "cut: ", cutoff
    nbands = str(from_path_get_nbands_thisstruct(path=path))
    if printoutput:
        print "nbands: ", nbands
    ediff = from_path_string_details_ediff(path=path)
    if printoutput:
        print "ediff: ", ediff
    ngxf = from_path_string_details_ngxf(path=path)
    if printoutput:
        print "ngxf: ", ngxf
    ismear = from_path_string_details_ismear(path=path)
    if printoutput:
        print "ismear: ", ismear
    #    print "INCAR     :"
    #    print "    CUT   :",cutoff
    #    print "    NBANDS:",nbands
    #    print "    EDIFF :",ediff
    #    print "    NGXF  :",ngxf
    #    print "    one issue with reference is that one should look in the non-reference structure and calc ..."
    #    print "    ... corresponding EDIFF and ...."
    job = from_path_string_job(path=path)
    if printoutput:
        print "job:", job
    #####################
    ## murn
    #####################
    murnadd = "!"
    relx = from_path_get_isrelaxjob(path=path)
    if printoutput:
        print "relx: ", relx
    if relx == True:
        murnadd = ""

    always = [
        "  PREC   =Accurate",
        "  ISTART =0             !WAVECAR:0-new,1-cont,2-samecut",
        "  ICHARG =2             !charge:1-file,2-overlapping_atom,10-const",
        "  ",
        "  ENCUT  =" + cutoff,
        "  NBANDS =" + nbands,
        "  LREAL  =.FALSE.",
        "  NELMDL =-5",
        "  EDIFF  =" + ediff,
        "  ALGO   =FAST",
        "  ISMEAR =" + ismear,  # -1 Fermi; 0 Gaussian; 1..N Methfessel-Paxton; -4-5 tetrahedron w/wo Bloechel corrections
        "  SIGMA  =0.1",
        "  ",
        "!Murn",
        "  " + murnadd + "IBRION = 1",
        "  " + murnadd + "NSW    = 40",
        "  " + murnadd + "EDIFFG = -0.0001",
        "  " + murnadd + "ISIF =2        !2ion,7volume",
        "  ",
        "  LWAVE  =F",
        "  LCHARF =F",
        "  LVTOT  =F",
        "  LELF   =F",
        "  NWRITE =0",
        "  ",
        "NPAR = 1",
        "ADDGRID=.TRUE."
    ]

    if ngxf is not None:
        always.append("NGXF =" + str(ngxf))
        always.append("NGYF =" + str(ngxf))
        always.append("NGZF =" + str(ngxf))

    ########################################
    #### write stuff out
    ########################################
    file = "INCAR"
    if test == True:
        return True

    filew = open(file, "w")
    for i in always:
        filew.write(i + "\n")
    filew.close()
    return True


def create_vaspinput(path=None, test=False, printoutput=False):
    """
    -cv
    """
    if path is None:
        path = my.pwd()
    if printoutput:
        print "POTCAR---" * 10
    create_POTCAR(test=True, path=path)
    if printoutput:
        print "KPOINTS--" * 10
    create_KPOINTS(test=True, path=path)
    if printoutput:
        print "POSCAR---" * 10
    create_POSCAR(path=path, test=True, printoutput=printoutput)
    if printoutput:
        print "INCAR--" * 10
    create_INCAR(path=path, test=True, printoutput=printoutput)

    if test == True:
        return True

    #print "((((((((((((((((((((((((((((((((((((((("
    create_POTCAR(path=path)
    create_KPOINTS(path=path)
    create_POSCAR(path=path, printoutput=False)
    create_INCAR(path=path, printoutput=False)


def find_TAKE_folder(structure=None, path=None):
    """
    -ftf
    """
    if path is None:
        path = my.pwd()
    else:
        path = str(path)

    ###################################################################################
    ### structure   :   bulk? vak? divak?
    ### Which EVinet:   bulk?   -> EVinet
    ###                 vak?    -> EVinet_b_32+EVinet_d_31 oder EVinet_b_108+EVinet_d_107
    ###                 divak?  -> EVinet_b_32+EVinet_d_30 oder EVinet_b_108+EVinet_d_106
    #quit()
    if structure is None:
        structure = from_path_string_structure(path=path, exit=False)
        if structure is None:
            #structure = "all"
            iam = my.which_function_ami()
            out1 = eval(iam)(structure="bulk", path=path)
            out2 = eval(iam)(structure="vak", path=path)
            out3 = eval(iam)(structure="divak", path=path)
            #print type(out1),out1
            #print type(out2),out2
            #print type(out3),out3
            return out1 + "\n" + out2 + "\n" + out3
            #print "structure:",structure
    #quit()
    #print structure,type(structure)
    if structure != "bulk" and structure != "vak" and structure != "divak":
        my.exit(error="Which structure: bulk, vak or divak? ")
    #print structure,type(structure)

    ###################################################################################
    ### Which EVinet:   fcc?    -> murn_xxxx_fcc
    ###                 bcc?    -> murn_xxxx_bcc
    try:
        cell = get.from_path_string_cell(path=path)
    except:
        cell = "*"
    if cell == "fcc4" or cell == "fcc1":
        cell = "fcc"
        #my.exit(error="YOU need to be in an folder which displays the cell: fcc, bcc, ... , EVinet for bcc or fcc?")
    #print "cell",cell

    ###################################################################################
    ### search for TAKE
    ###
    elementpath = get_path_element(path=path)
    import glob

    #paths_string = str(elementpath)+"/dynmat_"+structure+"_"+str(cell)+"*"+"/*TAKE*/Fqh*"
    paths_string = str(elementpath) + "/murn_" + structure + "_" + str(
        cell) + "*" + "/*TAKE*/EVinet*"
    #paths_string = str(elementpath)+"/dynmat_"+structure+"_"+str(cell)+"*"+"/*TAKE*/Fqh*"
    #print "string:",paths_string
    checkpath = glob.glob(paths_string)

    if len(checkpath) == 0:
        my.exit(error="Counldnt find TAKE in " + str(paths_string))

    ###################################################################################
    ### bulk
    if structure == "bulk" and len(checkpath) == 1:
        #print "in"
        file = my.checkfile(checkpath[0])
        return file  # +"\n"
    if structure == "bulk" and len(checkpath) > 1:
        my.exit(error="Found more than one bulk EVinet " + str(checkpath))

        ###################################################################################
        ### vak & divak
    ### 1) sicherstellen dass zu jedem EVinet_d_xy auch ein passendes EVinet_b_xy da ist
    ### 2) sicherstellen dass nix doppelt ist  e.g. 2mal EVinet_b_107
    ### 3) sicherstellen dass gleich vilele EVinet_d_xxx wie EVinet_b_xxx
    ### ANMERKUNG) hier werden erstmal alle Zellengroessen returned

    ### 1) sicherstellen dass zu jedem EVinet_d_xy auch ein passendes EVinet_b_xy da ist
    out = []
    for i in checkpath:
        EVineti = "EVinet" + i.split("EVinet")[1]
        EVinetpathi = i.split("EVinet")[0]
        bdi = EVineti.split("EVinet_")[1].split("_")[0]
        if bdi == "d":
            EVinetpathi = EVinetpathi[0:-1]
            #print ""
            #print "YYYi:",i,"||",EVinetpathi
            out.append(i)
            bdjlist = []
            for j in checkpath:
                EVinetj = "EVinet" + j.split("EVinet")[1]
                EVinetpathj = j.split("EVinet")[0]
                bdj = EVinetj.split("EVinet_")[1].split("_")[0]
                if bdj == "b":
                    bdjlist.append(j)
                    #print "    YYY:",j,"||",EVinetpathj
                    #print "passend:",bdjlist
            outcheck = my.grep(EVinetpathi, bdjlist)
            if len(outcheck) == 1:
                out.append(outcheck[0])

    ### 2) sicherstellen dass nix doppelt ist  e.g. 2mal EVinet_b_107
    EVinetlist = ["EVinet" + i.split("EVinet")[1] for i in out]
    if my.get_duplicate_items(EVinetlist) != []:
        for i in checkpath:
            print i
        my.exit(error="Doppelter EVinet: " + str(
            my.get_duplicate_items(EVinetlist)))

    ### 3) sicherstellen dass gleich vilele EVinet_d_xxx wie EVinet_b_xxx
    EVinetlistb = []
    EVinetlistd = []
    for i in EVinetlist:
        #print "I:",i
        bd = i.split("EVinet_")[1].split("_")[0]
        atoms = i.split("EVinet_")[1].split("_")[1]
        if bd == "d":
            EVinetlistd.append(i)
        if bd == "b":
            EVinetlistb.append(i)
        #print "Eb:",EVinetlistb
    #print "Ed:",EVinetlistd
    if len(EVinetlistd) != len(EVinetlistb):
        my.exit(error="Unequal length of EVinet_d_xxx and EVinet_b_xxx _b_: " +
                str(EVinetlistb) + "  _d_: " + str(EVinetlistd))

    ###################################################################################
    ### OUTPUT
    string = ("".join([i + "\n" for i in out]))
    #print string
    #print "--"
    return string[0:-1]

    #def get_jobfolderlist_to_create_from_onejobfolder(path=None,rise_error_if_exists=True,silent=False):
    #    """
    #    -gjtc looks in current folder and subfolder for ALATS file and prints out jobfolder to be created
    #    """
    #    #print "path:ganz:",path
    #    if path == None:
    #        path = my.pwd()
    #    else:
    #        path = str(path)
    #
    #    from_path_string_details(path=path)  ## just to ensure we are in details
    ##    print "danach:"
    ##    print ""
    #    try:
    #        ang = from_path_get_ang(path=path)      ## if we are in ang folder we dont want to create all! angfolders but just this
    #    except: # catch *all* exceptions
    #        ang = None
    #
    #    #print "ang:",ang
    ##    print ""
    #    if ang == None:
    #        alats = from_path_get_ALATS(string=True)
    #    else:
    #        alats = [str(ang)]
    #    #print "alats:",alats
    ##    print "riseerror:?",rise_error_if_exists
    #
    #    alatsFolder = [path+"/"+x+"Ang" for x in alats]
    #    #print "alatsFolder:",alatsFolder,len(alats),str(path.split("/")[-1]),str(alats[0])+"Ang"
    #    ########################################
    #    ### check if path has already x.xxAng/
    #    #print "YYY",path.split("/")[-1],len(alats)
    #
    #    if len(alats) == 1 and str(path.split("/")[-1]) == str(alats[0])+"Ang":
    #        alatsFolder = [path for x in alats]
    #    #print "alatsFolder:",alatsFolder
    #
    #    outfolder = []
    #    for folder in alatsFolder:
    ##        print "folder:",folder
    #        if my.isdir(folder) == True and my.pwd() != folder:    ## WENN er existier, nicht hinzufuegen, ausser wenn man gerade in dem folder ist!
    #            if silent == False: print folder," already done"
    #            if rise_error_if_exists == True:
    #                folderout="/".join(folder.split("/")[-2:])
    #                #               print "F>>>>>>",folderout
    #                return my.exit(error=" Folder already exists1: .../"+folderout)
    #        else:                           ## ANSONSTEN hinzufuegen
    #            if silent == False: print folder
    #            outfolder.append(folder)
    ##    print "OUT:",outfolder
    #    return outfolder

    #def check_jobfolder_to_create(allfolderlist=None,onepath=None,rise_error_if_exists=True,silent=False,writeinfo=False,direct_position_file=None):
    #    """
    #    -cjtc checks if POSCAR, KPOINTS, INCAR, POTCAR can be created in given jobfolder.
    #    """
    #
    #    if allfolderlist == None and onepath == None:
    #        onepath = my.pwd()
    #    if allfolderlist != None and onepath != None:
    #        my.exit(er2x2x2sc_16x16x16kp_400eV/ror="either or or, not both")
    #    if onepath != None:
    #        allfolderlist = [onepath]
    #
    ##    print "allfolderlist1:",allfolderlist
    #    if allfolderlist == None:
    #        allfolderlist = get_jobfolderlist_to_create_from_onejobfolder(path=onepath,rise_error_if_exists=rise_error_if_exists,silent=silent)
    #
    #
    #
    #    if writeinfo == True: print "info allfolderlist:",allfolderlist
    #
    #    outfolder = []
    #    for foldername in allfolderlist:
    #        #### check if folder does exist
    #        print "  -->> foldername:",foldername
    #        if my.isdir(foldername) == True and my.pwd() != foldername:
    #            if dont_create_if_exists == False:  ### das ist nicht definiert
    #                foldernameout = foldername.split("/")
    #                print "FA------>",foldernameout
    #                my.exit(error=" Folder already exists: "+str(foldernameout[-1]))
    #        outfolder.append(foldername)
    #
    #        if writeinfo == True: print "info POTCAR ..."
    #        create_POTCAR(test=True,path=onepath)
    #        if writeinfo == True: print "info POTCAR OK"
    #
    #        if writeinfo == True: print "info KPOINTS ..."
    #        create_KPOINTS(test=True,path=onepath)
    #        if writeinfo == True: print "info KPOINTS OK"
    #
    #        if writeinfo == True: print "info POSCAR ..."
    #        create_POSCAR(test=True,path=foldername,direct_position_file=direct_position_file)
    #        if writeinfo == True: print "info POSCAR OK"
    #
    #        if writeinfo == True: print "info INCAR ..."
    #        create_INCAR(path=foldername,test=True)
    #        if writeinfo == True: print "info INCAR OK"
    #    return outfolder


def create_jobs(infolder=None, joblist=None, rise_error_if_exists=False, silent=True):
    """
    -cj
    """
    if infolder is not None and joblist is not None:
        my.exit(error="please specify just onefolder or folderlist, not both")

    if infolder is None and joblist is None:
        infolder = my.pwd()

    ## jobList:
    jobList = my.pwd() + "/jobList"

    if joblist is None:
        joblist = get_jobfolderlist_to_create_from_onejobfolder(path=infolder, rise_error_if_exists=rise_error_if_exists, silent=silent)

    for i in joblist:
        print i

    if len(joblist) == 0:
        return "No jobs"
    #print joblist,len(joblist)

    #print ""
    #print ""
    print "jobList:", jobList
    #print ""
    #print ""

    import os
    if my.isfile(jobList):
        os.remove(jobList)
    #print "___>>>>>",joblist

    for job in joblist:
        if os.path.isdir(job) != True:
            os.makedirs(job)
        else:
            continue
        print "----> creating job:", job
        create_POTCAR(path=job)
        create_KPOINTS(path=job)
        create_POSCAR(path=job)
        create_INCAR(path=job)

        if my.pwd() != job:
            import shutil
            shutil.move("POTCAR", job)
            shutil.move("POSCAR", job)
            shutil.move("KPOINTS", job)
            shutil.move("INCAR", job)

    with open(jobList, "w") as jobListSaveWrite:
        for job in joblist:
            jobListSaveWrite.write(job + "\n")

    return True


    #noinspection PyUnboundLocalVariable,PyUnboundLocalVariable
def get_cell(path=None):
    """
    -gc
    """
    if path is None:
        path = my.pwd()
    #    print "in get_cell: path: ",path
    lattice_const = from_path_get_ang(path=path)
    #print "l",lattice_const
    ### Bei jeder Zelle muss entweder die lattice constand, volumen, bekannt sein!!!!

    cellpath = get_path_db_utils_cell_cellfile(path=path)
    #print "cellpath:",cellpath
    cell = my.readfile(cellpath)
    #    print "cell:",cell
    outcell = []

    ########## VECTOR
    for vec in cell:
        outvec = []
        coordsvec = vec.split()
    #        print ""
    #        print "::coordsvec:",coordsvec
        ############ COORDINATE
        for single_coordinate in coordsvec:
    #            print "-->", single_coordinate
            coordout = []
            single_coordinate_string = single_coordinate.split("*")
            for string in single_coordinate_string:
    #                print "string:",string
    #                print "   --> coordin:",coordin
    #
                if my.is_number(string) == True:
                    value = my.convert_to_number(string)
                    coordout.insert(1, value)
    #
                if my.is_number(string) != True:
                    function_in = string.split("(")
                    function_in_arguments = string.split("(")
    #                    print "(ff",function_in,")"
    #                    print "(fa",function_in_arguments,")"
    #                    print "(",len(function_in),")"
                    if len(function_in) == 1:
                        function = function_in[0]
                        function_arguments = []
                    if len(function_in) >= 2:
                        function = function_in[0]
                        function_arguments = [function_in_arguments[
                            1].split(")")[0]][0].split(",")
                        function_arguments = ','.join(function_arguments)
    #                    print "function:",function
    #                    print "function_arguments",function_arguments

                    ## check if function exists
                    if function not in alldefs(script=my.which_scriptpath_isit_whew_i_am()):
                        my.exit(error=" function " + function +
                                " is not in alldefs: " + str(alldefs()))

                    ## get value
    #                    print "function_arguments:",function_arguments
    #                    print "function          :",function
                    #noinspection PyUnboundLocalVariable,PySimplifyBooleanCheck
                    if path is None:
    #                        print "firstif"
                        if function_arguments == []:
                            value = eval(function)()
                        #noinspection PySimplifyBooleanCheck
                        if function_arguments != []:
                            value = eval(function)(function_arguments)
    #                    print "value:",value
                    if path is not None:
    #                        print "hier Drin"
                        if function_arguments == []:
                            value = eval(function)(path=path)
                        #noinspection PySimplifyBooleanCheck
                        if function_arguments != []:
                            value = eval(
                                function)(function_arguments, path=path)
                        #                    print "value:",value

                    ## check if value is number
    #                    print "value:",value
                    if my.is_number(value) != True:
    #                        return my.error_write(value)
                        return my.error_write("evaluated " + value + " is not a number ")
                            #" evaluated: (",value,') and not a number from function: (',function, ")with args:(",function_arguments,")")

                    ## convert value to number
                    #noinspection PyUnboundLocalVariable
                    value = my.convert_to_number(value)
    #                    print "value:",value

                    ## add value to coordout
                    coordout.append(value)

    #            print "   ---> coordout:",coordout#,out#,result
    #            print ""
            ### multiply all coords
            out = 1
            for i in coordout:
                out = out * i
                #print "i,",i,"out,",out
            coordout = out
    #            print "CCC:",out
            outvec.append(out)
        #print "outvec:",outvec
        outcell.append(outvec)
    return outcell


def get_positions_direct(path=None):
    """
    -gpd
    """
    if path is None:
        path = my.pwd()
    dir_pos = my.readfile(get_path_positions_direct(path=path))
    return dir_pos


def positions_cartesian_to_reduced(path_coords=None, path_cell_cartesian=None, path=None, numpyarray_coords=None, array_cell=None, return_numpyarray=False):
    """
    -pcr
    """
    if path is None:
        path = my.pwd()

    coords_nrs = None
    cell_cartesian = None
    if path_coords is None and numpyarray_coords is None:
        my.exit(
            error="provide coords either as numpyarray or path to files(s)")
    if path_coords is not None and numpyarray_coords is not None:
        my.exit(
            error="provide coords either as numpyarray or path to files(s)")
    if path_cell_cartesian is not None and array_cell is not None:
        my.exit(error="provide cell as file or as numpyarray")

    ############ coords ###########################
    import numpy as np
    if path_coords is not None:
        my.checkfile(path_coords)
        coords_nrs = np.genfromtxt(path_coords)

    if coords_nrs is None:
        my.exit(error="coords?")
    #print coords_nrs

    ############# cell #############################
    if path_cell_cartesian is None:
        #print "path:",path
        cell_cartesian = get_cell(path=path)

    if cell_cartesian is None:
        my.exit(error="cell?")

    A_cart = np.array(cell_cartesian)

    ############# inverse transpose #############################
    import numpy.linalg as nl
    import scipy as sp
    #print A
    CI = nl.inv(A_cart)
    #CIT = sp.transpose(CI)
    CIT = CI.T
    A_red = CIT
    out = sp.dot(coords_nrs, A_red)

    if return_numpyarray == True:
        return out
    return my.numpyarray_to_screenoutput(numpyarray=out)


def positions_reduced_to_cartesian(path_coords=None, path_cell_cartesian=None, path=None, numpyarray_coords=None, array_cell=None, return_numpyarray=False, displace_atom=None, displace_xyz=None, displace_axis=None, displacement_ang=0.0052917721):
    """
    -prc
    """
    if path is None:
        path = my.pwd()

    coords_nrs = None
    cell_cartesian = None
    if path_coords is None and numpyarray_coords is None:
        my.exit(
            error="provide coords either as numpyarray or path to files(s)")
    if path_coords is not None and numpyarray_coords is not None:
        my.exit(
            error="provide coords either as numpyarray or path to files(s)")
    if path_cell_cartesian is not None and array_cell is not None:
        my.exit(error="provide cell as file or as numpyarray")

    ############ coords ###########################
    import numpy as np
    if path_coords is not None:
        my.checkfile(path_coords)
        coords_nrs = np.genfromtxt(path_coords)

    if coords_nrs is None:
        my.exit(error="coords?")
    #print coords_nrs

    ############# cell #############################
    #print "--cell--"
    if array_cell is not None:
        cell_cartesian = array_cell
    else:
        if path_cell_cartesian is None:
            #print "path:",path
            cell_cartesian = get_cell(path=path)
    #print "cell2"

    if cell_cartesian is None:
        my.exit(error="cell?")

    A_cart = np.array(cell_cartesian)
    #print "cell:"
    ############# inverse transpose #############################
    import numpy.linalg as nl
    import scipy as sp
    out = sp.dot(coords_nrs, A_cart)

    #    print "zuvor",displace_atom
    #    if displace_atom != None:
    #        if displace_xyz == None and displace_axis == None:
    #            my.exit(error="either axis or xyz")
    #        if displace_xyz != None and displace_axis != None:
    #            my.exit(error="either axis or xyz")
    #        if displace_axis != None:
    #            my.exit(error="todo yet")
    #        if displace_xyz != None:
    #            print out[int(atom)-1]
    #    quit()
    #print ""
    #out[6] = [1,2,3]
    #print out[6]
    #print ""
    #print out
    #print

    if return_numpyarray == True:
        return out
    return my.numpyarray_to_screenoutput(numpyarray=out)


def get_positions_cartesian(path=None, displace_1_1=None, displacement_1_1_angstrom=0.0052917721):
    """
    -gpc
    """
    if path is None:
        path = my.pwd()
    dir_pos = get_positions_direct(path=path)
    cell = get_cell(path=path)
    if from_path_string_cell(path=path) != "fcc4":
        my.exit(error=" chekc ob die cartesichen coordinaten auch fuer andere strukturen als fcc4 corret sind")
    #    print "cell:",cell
    #    print ""
    #print "pos:",dir_pos,type(dir_pos)
    #    print dir_pos[0],type(dir_pos[0])
    #    print ""

    cartes_all = []
    for idx, pos in enumerate(dir_pos):

    #        print pos
        x = float(pos.split()[0])
        y = float(pos.split()[1])
        z = float(pos.split()[2])  # [0],pos[1],pos[2]
        avec = cell[0]
        bvec = cell[1]
        cvec = cell[2]
        ax = float(avec[0])
        ay = float(avec[1])
        az = float(avec[2])
        bx = float(bvec[0])
        by = float(bvec[1])
        bz = float(bvec[2])
        cx = float(cvec[0])
        cy = float(cvec[1])
        cz = float(cvec[2])
        cartes_pos = x * ax + x * bx + x * cx, y * ay + y * by + \
            y * cy, z * az + z * bz + z * cz
    #        print "idx:",idx,"dis:",displace_1_1
        ### idx 0 ist das erste atom....
        if idx is 0 and displace_1_1 is True:
    #            print "in",cartes_pos
            print "DISPLACED ATOMD 1 in x!"
            cartes_pos = x * ax + x * bx + x * cx + displacement_1_1_angstrom, y * ay + y * by + y * cy, z * az + z * bz + z * cz
    #            print "out",cartes_pos
        cartes_pos_string = " ".join(map(str, cartes_pos))
        cartes_all.append(cartes_pos_string)
    #        print x,y,z," ||| ",cartes_pos
        #print x*ax #+x*bvec[0]+x*cvec[0]
        #out = " ".join(map(str, cartes_all[0])) #(list(cartes_all[0]))
    return cartes_all


def element_real(element=None, path=None):
    """
    -er ## this does not care about correct sorting of the elements
    """
    if path is None:
        path = my.pwd()
    else:
        path = str(path)

    if element is None and path is None:
        my.exit(error=" either element or path")
    if element is not None and path is not None:
        my.exit(error=" either element or path")

    if element is None:
        element = from_path_string_element(path=path)
    element = str.replace(element, "binary_", "")
    element = str.replace(element, "-", " ")
    return element  # TODO STILL not finished


def db_NBANDS_OCCUPIED(path=None):
    if path is None:
        path = my.pwd()
    atoms = from_path_get_numatoms(path=path)
    return sys, exit("TODO")  # TODO

###########################
## fromfile
###########################


def KPOINTS_kpoints():
    """
    -Kk
    """
    file = "KPOINTS"
    out = my.readfile(file)
    kp = out[3].split()
    if len(kp) != 3:
        print "len kp nicht 3 kp:", kp
        my.exit()
    for i in kp:
        if my.is_int(i) != True:
            print "is not int :", kp
            my.exit()
    return kp


def KPOINTS_kpointsproduct():
    """
    -Kkp
    """
    kp = KPOINTS_kpoints()
    return int(kp[0]) * int(kp[1]) * int(kp[2])


def check_path_is_potentialpath(path=None, exit=False):
    """
    -cppot
    """
    if path is None:
        path = str(my.pwd())
    avail_v(path=path)  # just to ensure we are in v
    #print "path:",my.pwdsplit()
    if len(my.pathsplit(path=path)) == 5:
        return True
    else:
        if exit == True:
            my.exit(error="Your are not in POTENTIALpath")
        else:
            return False


def check_path_is_elementpath(path=None, exit=False):
    """
    -cpe
    """
    if path is None:
        path = str(my.pwd())
    avail_v(path=path)  # just to ensure we are in v
    #print "path:",my.pwdsplit()
    if len(my.pathsplit(path=path)) == 6:
        return True
    else:
        if exit == True:
            my.exit(error="Your are not in ELEMENTpath")
        else:
            return False


def check_path_is_projectpath(path=None, exit=False):
    """
    -cpp
    """
    if path is None:
        path = str(my.pwd())
    #print "pp",path

    avail_v(path=path)  # just to ensure we are in v
    if len(my.pathsplit(path=path)) == 7:
        return True
    else:
        if exit == True:
            my.exit(error="Your are not in PROJECTpath")
        else:
            return False


def check_path_is_subprojectpath(path=None, exit=False):
    """
    -cps
    """
    if path is None:
        path = str(my.pwd())

    avail_v(path=path)  # just to ensure we are in v
    if len(my.pathsplit(path=path)) == 8:
        return True
    else:
        if exit == True:
            my.exit(error="Your are not in SUBPROJECTpath")
        else:
            return False


def create_NBANDS(path=None):
    """
    -cn   goes through PROJECT and writes/addes/changes NBANDS list: atoms,NBANDSOCCUPIED
    """
    ## soll duch alle subprojects gehen und NBANDS file schreiben!
    if path is None:
        path = my.pwd()

    avail_v(path=path)  # ensure you are at least in v

    ### find out: are you in projectpath or in subprojectpath
    current_in_potential = check_path_is_potentialpath(path=path)
    current_in_element = check_path_is_elementpath(path=path)
    current_in_proj = check_path_is_projectpath(path=path)
    current_in_subproj = check_path_is_subprojectpath(path=path)

    check_if_only_one_TRUE = [current_in_potential, current_in_element,
                              current_in_proj, current_in_subproj]
    #print "CHECK:",check_if_only_one_TRUE
    if check_if_only_one_TRUE.count(True) != 1:
        print "OUT:", check_if_only_one_TRUE
        my.exit(error="Are you in Elementpath, Projectpath or subprojectpath?" + str(check_if_only_one_TRUE))

    ## get all subprojects: either just the one you are in
    if current_in_element == True:
        all_subproj = []
        allproj = get_path_all_projects(path=path, list=False)
        for i in allproj:
            all_subproj_add = get_path_all_subprojects(path=i)
            all_subproj.extend(all_subproj_add)

    if current_in_potential == True:
        all_subproj = []
        my.exit(error="currently not implemented TODO")
    #        allelements =
    #
    #        allproj = get_path_all_projects(path=path)
    #        for i in allproj:
    #            all_subproj_add = get_path_all_subprojects(path=i)
    #            all_subproj.extend(all_subproj_add)

    if current_in_proj == True:
        #print "current projectpath:",path
        all_subproj = get_path_all_subprojects(path=path)
        #print "all_subproj:", all_subproj

    if current_in_subproj == True:
        all_subproj = [path]

    #print "2",all_subproj
    #for jj in all_subproj:
    #    print jj
    #print ""
    #print ""
    #quit()
    file = "NBANDS"
    #print "all:", all_subproj
    if not all_subproj:
        my.exit(error="which all_subproj?")
    for subproj in all_subproj:
        ##########################################
        ### wenn NBANDS file besteht -> fertig, ansonsten NBANDS file schreiben
        ##########################################
        ## sowas wie cd path
        #print ""
        #print "hallo,",os.path.exists(subproj),str(subproj)
        if os.path.exists(subproj) != True:
            continue

        os.chdir(subproj)
        NBANDS_file = str(subproj) + "/" + str(file)
        isNBANDS = my.isfile(NBANDS_file)
        from time import gmtime, strftime
        if isNBANDS == True:
            #print "--> ",subproj," NBANDS exist"
            print strftime("%Y-%m-%d %H:%M:%S", gmtime()), "exist ", subproj, " NBANDS exist"
            continue
        #from datetime import datetime

        print strftime("%Y-%m-%d %H:%M:%S", gmtime()), "----> ", subproj, my.pwd()

        out = my.run("OUTCAR_A.sh OUTCAR_number_of_atoms.sh OUTCAR_KPOINTS_OCCUPIED-max-occ-band.sh")
        open(NBANDS_file, 'w').close()
             ### besser hier schon erstellen, wenn die OUTCAR nix greppen kann
                                            ### dann muss in diesem Folder nicht wieder OUTCAR_A.sh ... laufen
        atoms_nbands = []
        for line in out.split("\n"):
            if len(line.split()) != 3:
                continue  # geth zur nexten line
            if str(line.split()[0]) == "pfad":
                continue
            atoms = line.split()[1]
            NBANDS = line.split()[2]
            if my.is_int(atoms) != True:
                continue
            if my.is_int(NBANDS) != True:
                continue
            #print atoms,NBANDS,type(atoms),type(NBANDS),my.is_int(atoms),my.is_int(NBANDS)
            #print "line:",line," ||| ",len(line.split())," +++ ",line.split()[0]
            #with open(NBANDS_file,'w') as f:
            #    f.write("hallo")
            #f.closed
            atoms_nbands.append(str(atoms) + " " + str(NBANDS))

            ########################
            ### WRITE TO FILE
            ########################
            open(NBANDS_file, 'w').close()
            with open(NBANDS_file, 'w') as f:
                for line in atoms_nbands:
                    f.write(str(line) + "\n")
            #print atoms_nbands
            #print f.closed
            #print subproj
    return True


def dos(infile=None, outfile=None):
    """
    -dos create density of states
    """
    #infile = "positionsxyz_aufnull"
    if infile is None:
        exit("Please define an inputfile")
    if outfile is None:
        outfile = "dos.dat"

    input = my.readfile(infile)
    #print input
    #print len(input)
    #quit()
    inputp = [float(i) for i in input]

    inputm = [-float(i) for i in input]
    inputf = []
    inputf.extend(inputp)
    inputf.extend(inputm)

    inputf = input

    import matplotlib.pyplot as plt
    import numpy as np
    import matplotlib.mlab as mlab
    from math import sqrt

    #    mean = 0
    #    variance = 1
    #    sigma = sqrt(variance)
    #    x = np.linspace(-3,3,100)
    #    plt.plot(x,mlab.normpdf(x,mean,sigma))
    #    plt.show()
    #    return True

    from numpy.random import standard_normal
    #   data = standard_normal(100)
    import statistics
    import numpy as np
    y, x = statistics.pdf(inputf)
    import matplotlib.pyplot as plt
    #print "x,y:",x,y,type(x)
    plt.plot(x, y)
    #np.savetxt('test.txt', [x,y])
    out = []
    for i, a in enumerate(x):
        out.append([a, y[i]])
        print a, y[i]
    #print out
    np.savetxt(outfile, out, fmt='%10.10f  %20.20f')
    return True
    #plt.show()

    #    import matplotlib.pyplot as plt
    #    import numpy as np
    #    import matplotlib.mlab as mlab
    #    #from scipy.stats import norm
    #    import scipy
    #    import scipy.stats
    #
    #    import numpy as np
    #    import matplotlib.pyplot as plt
    #    from scipy.stats import norm
    #
    #    # Plot between -10 and 10 with .001 steps.
    #    range = np.arange(-10, 10, 0.001)
    #    # Mean = 0, SD = 2.
    #    plt.plot(range, norm.pdf(range,0,2))
    #    y = pdf(inputf)
    #    #plt.show()
    return max(inputf), min(inputf)


def Fqhfit_to_Fqh_vs_V():
    """
    -c
    """
    atoms = ""
    alat1 = 3.6
    alat2 = 3.75
    sc = 2
    V1 = my.iroundup((alat1 * sc) ** 3)
    V2 = my.iroundup((alat2 * sc) ** 3)
    #print "a",a1,a2
    v = range(V1, V2, 1)
    #print "V",v
    Fqh = my.readfile(
        "Fvib_fromExactFreqs" + str(atoms) + "_fit_order2")  # Fqh_surface_fit")
    for line in Fqh:  # every temperature
        words = line.split()
        if len(words) < 2:
            continue
        t, a, b, c = words[0], words[1], words[2], words[3]
        if my.is_number(t) != True:
            continue
        t, a, b, c = float(
            words[0]), float(words[1]), float(words[2]), float(words[3])
        if t != 1101:
            continue

        #print t,a,b,c
        #print ""
        import math
        for V in v:  # go through volume
            #print V, (a+b*V+c*V*V)
            #a=V**(1/3)
            #print a, (a+b*V+c*V*V)
            a = (V ** (1 / 3.)) / sc
            Vpa = a ** 3
            #print a, (a+b*Vpa+c*Vpa*Vpa)
            print a, (a + b * V + c * V * V)

    ###############################################################################
    ## ^ ^ ^  ^ ^  ^ ^
    ## definitions stop
    ###############################################################################
    ############################################################################
    ## main
    ############################################################################
if __name__ == "__main__":
    #print "yes"
    p = argparse.ArgumentParser("TODO: all ERRORS to my.exit()")
    g = p.add_mutually_exclusive_group(required=True)  # besser hinzufuegen da ansonsten so sachen wie -one 8 jj -two erlaubt sind
    #    print my.get_definitions_of_current_script_additionalinfo()
    #print "yo"
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
            g.add_argument('--' + definition, short, dest='def_to_run', action='store_const', const=definition, help=docstring + arguments)
        else:
            g.add_argument('--' + definition, dest='def_to_run', action='store_const', const=definition, help=docstring + arguments)

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
