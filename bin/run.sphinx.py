#!/usr/bin/env python

# Parameters
testalats = False
prog = "sphinx"
scs = [ 2, 3 ]
scs = [ 2, 3, 4, 5, 6, 7 ]
scs = [ 3 ]
string = ""





import socket
import sys
import glob
import os
import numpy as np
import shutil
import subprocess

def getinfo():
    SPHINX = "/home/glensk/sphinx_bin/sphinx.sh"
    if os.path.isfile(SPHINX) != True:
        sys.exit("SPHINX executable "+str(SPHINX)+" not found")
    hostname = None
    hostname = socket.gethostname()
    rooth="/home"
    if hostname == "mac": rooth = "/data"
    if hostname == "cmpc03.mpie.de": rooth = "/data"
    if hostname == "cmmmd002.mpie.de": rooth = "/home"
    print "hostname:", hostname
    if hostname == None: sys.exit("no hostname found")
    root = rooth+"/glensk/v/pp/cu/ti_bulk_fcc4/__FIRST_TESTS__low_2x2x2sc_230eV_kp02m02m02_sphinx"
    if os.path.isdir(root) != True:
        sys.exit("could not find root folder")
    print "root:",root
    home = os.environ['HOME']
    print "home:",home
    path = os.getcwd()
    print "path:",path
    try:
        pot = path.split("sphinx/2013_10_11_")[1].split("/")[0]
    except IndexError:
        sys.exit("you have to go to murn_bulk_fcc4 or ti_vak_fcc4 folder or similar")
    if pot == "imd_ercolessi":
        pot = "ercolessi"
    print "pot:",pot
    try:
        element = path.split("sphinx/2013_10_11_")[1].split("/")[1]
    except IndexError:
        sys.exit("you have to go to murn_bulk_fcc4 or ti_vak_fcc4 folder or similar")
    print "element:",element
    try:
        proj = path.split("sphinx/2013_10_11_")[1].split("/")[2]
    except IndexError:
        sys.exit("you have to go to murn_bulk_fcc4 or ti_vak_fcc4 folder or similar")
    print "proj:",proj
    job = proj.split("_")[0]
    print "job:",job
    mb = proj.split("_")[1]
    print "mb:",mb
    # checks
    if element not in [ "al", "cu" ]:
        sys.exit("element not al or cu")
    if pot not in ["mishin","ercolessi","meidavenport","meidavenport-largerCutoff"]:
        sys.exit("pot ("+str(pot)+") not found")
    if job not in [ "murn", "dynmat", "ti" ]:
        sys.exit("job not found")
    if mb not in [ "bulk", "vak" ]:
        sys.exit("bulk or vak not found")
    if element == "al" and job == "dynmat":
        alats = np.arange(3.94,4.16,0.01)
        if testalats == True:
            #alats = np.arange(3.94,3.96,0.02)
            alats = [3.96]
            #alats = np.arange(3.94,3.96,0.02)
            #alats = np.arange(4.04,4.16,0.02)
            #alats = np.arange(4.04,4.16,0.02)
        #if job == "murn":
        #   alats = np.arange(3.94,4.16,0.02)  # for bulk and vak
    if element == "al":
        mass = 26.982
        elname = "Al"
        elnamel = "Aluminum"

    if element == "al" and job == "dynmat":
        alats = np.arange(4.14,4.16,0.02)
    if element == "al" and job == "murn":
    	alats = np.arange(4.02, 4.06, 0.005)
    if element == "cu":
        alats = np.arange(3.5,3.75,0.01)
        mass = 63.546
        elname = "Cu"
        elnamel = "Copper"
    inputsxt = root+str("/input_template.sx")
    print "element:",element
    print "mass:",mass
    print "alats:",alats
    print "inputsxt:",inputsxt
    if os.path.isfile(inputsxt) != True:
        sys.exit("could not find input_template.sx")
    details = path.split(root)[1].split("/")
    rootd = root+"/"+details[1]+"/"+details[2]+"/"
    print "rootd:",rootd
    return path, pot, element, proj, job, mb, alats, inputsxt, mass, elname, elnamel, root, SPHINX, rootd


def density_embeded_pair(pot,root):
    density = None
    potpot = root+"/SphinxEmpiricalPotentials/"
    sphinx_extractForces = potpot+"/sphinx_extractForces.sh"
    if os.path.isfile(sphinx_extractForces) != True:
        sys.exit("sphinx_extractForces.sh")
    if pot == "ercolessi" and element == "al":
        density = potpot+"al_ercolessi.den.dat"
        embeded = potpot+"al_ercolessi.embeded.dat"
        pair = potpot+"al_ercolessi.pair.dat"
    if pot == "mishin" and element == "al":
        density = potpot+"al_mishin.den.dat"
        embeded = potpot+"al_mishin.embeded.dat"
        pair = potpot+"al_mishin.pair.dat"
    if pot == "meidavenport" and element == "al":
        density = potpot+"al_meidavenport-largerCutoff.den.dat"
        embeded = potpot+"al_meidavenport.embedded.dat"
        pair = potpot+"al_meidavenport-largerCutoff.pair.dat"
    if density == None:
        sys.exit("density missing in density_embeded_pair() function")
    if os.path.isfile(density) != True:
        sys.exit("density missing")
    if os.path.isfile(embeded) != True:
        sys.exit("embeded missing")
    if os.path.isfile(pair) != True:
        sys.exit("pair missing")
    print "density:",density
    print "embeded:",embeded
    print "pair:",pair
    return density, embeded, pair, potpot, sphinx_extractForces


def get_inputsx(folder,a,atomsx):
    inputsx = folder+'/input.sx'
    #print "inputsxt:",inputsxt
    #print "inputsx :",inputsx
    shutil.copyfile(inputsxt, inputsx)
    cmd = "sed -i 's/%s/%s/' %s"%("xxxmassxxx", str(mass), inputsx)
    os.system(cmd)
    cmd = "sed -i 's/%s/%s/' %s"%("xxxelnamelxxx", str(elnamel), inputsx)
    os.system(cmd)
    cmd = "sed -i 's/%s/%s/' %s"%("xxxelnamexxx", str(elname), inputsx)
    os.system(cmd)
    cmd = "sed -i 's/%s/%s/' %s"%("xxxALATxxx", str(a), inputsx)
    os.system(cmd)
    cmd = "sed -i 's/%s/%s/' %s"%("xxxscxxx", str(sc), inputsx)
    os.system(cmd)
    cmd = "sed -i 's|%s|%s|' %s"%("density", str(density), inputsx)
    os.system(cmd)
    cmd = "sed -i 's|%s|%s|' %s"%("embeded", str(embeded), inputsx)
    os.system(cmd)
    cmd = "sed -i 's|%s|%s|' %s"%("pair", str(pair), inputsx)
    os.system(cmd)
    cmd = "sed -i 's|%s|%s|' %s"%("xxxatom.sx", str(atomsx), inputsx)
    os.system(cmd)
    #print "a:",a
    #print "embeded:",embeded
    #print "pair:",pair
    #print "at:",atoms
    return


def get_xxatomsx(folder,filename):
    xxatomsx = potpot+str(filename)
    # print "folder:",folder
    # print "filename:",filename
    # print "potpot:",potpot
    # print "xxatomsx:",xxatomsx
    if os.path.isfile(xxatomsx) != True:
        sys.exit(str(xxatomsx)+" does not exist")
    shutil.copyfile(xxatomsx, str(folder)+'/'+str(filename))
    return


def get_lammps_structure(alat):
    #print "alat:",alat
    from ase.lattice.cubic import FaceCenteredCubic
    bulk = FaceCenteredCubic(size=[sc,sc,sc], symbol='Al', directions=((1,0,0), (0,1,0),(0,0,1)), latticeconstant = alat, pbc = 1)
    print "####### mb:",mb
    if mb == "vak":
        del bulk[1]   # this is when we have a vacancy
    #from ase.calculators.lammps import write_lammps_data  # depends on versin
    from ase.calculators.lammpsrun import write_lammps_data
    write_lammps_data('positions'+str(alat),bulk)
    return

def get_lammps_alpot(folder):
    filename='Al-LEA.eam.alloy'
    alpot="/home/glensk/lammps/"+filename
    if os.path.isfile(alpot) != True:
        sys.exit("al pot lammps does not exist")
    shutil.copyfile(alpot, folder+"/"+filename)
    return


def get_lammps_infile(folder,alat):
    infile="/home/glensk/lammps/in_file"
    if os.path.isfile(infile) != True:
        sys.exit("infile lammps does not exist")
    shutil.copyfile(infile, folder+"/in_file")
    cmd = "sed -i 's|%s|%s|' %s"%("read_data.*", "read_data positions"+str(alat), folder+"/in_file")
    os.system(cmd)
    return
    

def sphinx_run():
    subprocess.call([SPHINX],shell=True)
    if os.path.isfile("sphinx.log") != True:
        sys.exit("SPHINX did not create sphinx.log file")

def sphinx_ene(filename):
    hartreetoev=27.211384
    if os.path.isfile(filename) != True:
        sys.exit("did not find "+str(filename))
    enecmd = "grep TOTAL "+str(filename)+" | awk '{print $4}'"
    proc = subprocess.Popen([enecmd], stdout=subprocess.PIPE, shell=True)
    (ene, err) = proc.communicate()
    return float(ene)*hartreetoev

def lammps_ene(filename):
    if os.path.isfile(filename) != True:
        sys.exit("did not find "+str(filename))
    enecmd = "cat "+str(filename)+" | grep Loop -B 1 | head -1 | awk '{print $3}'" 
    proc = subprocess.Popen([enecmd], stdout=subprocess.PIPE, shell=True)
    (ene, err) = proc.communicate()
    return float(ene)
    

def sphinx_forces(filename):
    if os.path.isfile(filename) != True:
        sys.exit("did not find "+str(filename))
    forces = 0
    return forces


def create_murnjobs_sphinx():
    print "in"
    for alat in alats:
        volume = alat**3/4
        print "volume:",volume
        atomsx = str(atoms)+"atom.sx"
        print "alat:",alat,"sc:",sc,"atoms:",atoms,"atomsv:",atomsv,"atom.sx:",atomsx
        get_inputsx(os.getcwd(),alat,atomsx)
        get_xxatomsx(os.getcwd(),atomsx)

        sphinx_run()

        # auswerten
        filename = "sphinx.log."+str(alat)
        os.rename("sphinx.log",filename)
        eneperatom = sphinx_ene(filename)/atoms
        print "volume,ene:",volume,eneperatom
        open('energy.dat', 'a').close()
        with open("energy.dat", "a") as myfile:
                myfile.write(str(volume)+"  "+str(eneperatom)+"\n")
    return


def create_murnjobs_lammps():
    get_lammps_alpot(os.getcwd())
    for alat in alats:
        volume = alat**3/4
        get_lammps_structure(alat)
        get_lammps_infile(os.getcwd(),alat)
        
        #start lammps
        exe = "/home/glensk/lammps/lmp_cmmd < in_file" 
        proc = subprocess.Popen([exe], stdout=subprocess.PIPE, shell=True)
        (exe, err) = proc.communicate()

        
        filename="log.lammps"+str(alat)
        os.rename("log.lammps",filename)
        os.rename("trj_lammps","trj_lammps"+str(alat))
        eneperatom = lammps_ene(filename)/atoms
        print "volume,ene:",volume,eneperatom
        open('energy.dat', 'a').close()
        with open("energy.dat", "a") as myfile:
                myfile.write(str(volume)+"  "+str(eneperatom)+"\n")
    return


def create_murnjobs_imd():
    for alat in alats:
        volume = alat**3/4
        get_coords_imd(alat)
        #get_lammps_structure(alat)
        #get_lammps_infile(os.getcwd(),alat)
        
        #start lammps
        exe = "/home/glensk/lammps/lmp_cmmd < in_file" 
        proc = subprocess.Popen([exe], stdout=subprocess.PIPE, shell=True)
        (exe, err) = proc.communicate()

        
        filename="log.lammps"+str(alat)
        os.rename("log.lammps",filename)
        os.rename("trj_lammps","trj_lammps"+str(alat))
        eneperatom = lammps_ene(filename)/atoms
        print "volume,ene:",volume,eneperatom
        open('energy.dat', 'a').close()
        with open("energy.dat", "a") as myfile:
                myfile.write(str(volume)+"  "+str(eneperatom)+"\n")
    return

def get_coords_imd(alat):
    if mb == "bulk":
        add=""
    elif mb == "vak":
        add=" 1"
    else:
        sys.exit("specify bulk or vak")
    cmd="/home/glensk/v/pp/cu/ti_bulk_fcc4/__FIRST_TESTS__low_2x2x2sc_230eV_kp02m02m02_sphinx/imd_potentials/createIMDcoords.sh "+str(sc)+" "+str(alat)+str(add)
    os.system(cmd)
    shutil.copyfile("coords", "coords"+str(alat))
    print "mb:",mb
    return

def get_struct_sphinx():
    pass

def create_dynmatjobs(prog=None):
    comeback = os.getcwd()
    print "atoms:",atoms
    print "atomsv:",atomsv
    print "comeback:",comeback
    if prog == None:
        sys.exit("please specify sphinx or imd")
    if prog == "sphinx" or prog == "imd":
        pass
    else:
        sys.exit("please specify sphinx or imd")

    if mb == "bulk":
        iterate = atomsb
        posneg = [ 'forces+', 'forces-' ]
        if os.path.isdir(posneg[0]) == True:
            sys.exit("path "+posneg[0]+" exists")
        if os.path.isdir(posneg[1]) == True:
            sys.exit("path "+posneg[1]+" exists")

        for pn in posneg:
            os.chdir(comeback)
            os.makedirs(pn)
            os.chdir(pn)

            if prog == "sphinx":
                for alat in alats:
                    print "alat:",alat,"pn:",pn
                    if pn == posneg[0]: atomsx = str(atoms)+"atom+.sx"
                    if pn == posneg[1]: atomsx = str(atoms)+"atom-.sx"
                    get_inputsx(os.getcwd(),alat,atomsx)
                    get_xxatomsx(os.getcwd(),atomsx)
                    print "sphinx_run()"
                    sphinx_run()

                    # auswerten alle alats
                    if os.path.isfile("sphinx.log") != True:
                        sys.exit("shinx.log was not created")
                    filename = "sphinx.log."+str(alat)
                    os.rename("sphinx.log",filename)
                
                # auswerten alle extract forces zum schluss
                print "extr_forces",sphinx_extractForces
                out = subprocess.call([sphinx_extractForces])   
                if out != 0:
                    sys.exit("sphinx_extractForces did not work")

            if prog == "imd":
                for alat in alats:
                    print "alat:",alat,"pn:",pn
                    get_coords_imd(alat)


        os.chdir(comeback)
        # cp averageForces.math .
        # ausfUehren getAvgForces.sh PER SSH! (mathematica script)

    
    
    if mb == "vak":
        iterate = atomsv
        for alat in alats:
            os.chdir(comeback)
            stralat = str(alat)
            if alat == 4.0:
                stralat = "4"
            print "alat:",alat
            
            #folderrelcoords = "/home/glensk/v/pp/cu/ti_bulk_fcc4/__FIRST_TESTS__low_2x2x2sc_230eV_kp02m02m02_sphinx/2013_10_11_imd_ercolessi/al/murn_vak_fcc4/"+str(sc)+"x"+str(sc)+"x"+str(sc)+"sc_de-14/"
            murn_vak_coords = str(rootd)+"/murn_vak_fcc4/"+str(sc)+"x"+str(sc)+"x"+str(sc)+"sc/"+"/Al.final.chkpt."+stralat
            if os.path.isfile(murn_vak_coords) != True:
                sys.exit("ref pos do not exist:"+str(murn_vak_coords))
            print "GOT STRUCTURE FROM: ",murn_vak_coords," and written to relaxedCoords"
            cmd = "cat "+str(murn_vak_coords)+" | tail -n +9 | awk '{print $4,$5,$6}' > relaxedCoords_"+str(alat)
            #poseq = pylab.loadtxt("relaxedCoords_"+str(alat))
            os.system(cmd)
            folderalat = comeback+"/"+str(alat)
            if os.path.isdir(folderalat) == True:
                print "folderalat : "+folderalat+" already exists"
            else:
                os.makedirs(folderalat)
            os.chdir(folderalat)
            for atom in range(1,atomsv+1):
                for disp in range(1,3+1):
                    folderdisp=folderalat+"/"+str(atom)+"_"+str(disp)

                    if os.path.isdir(folderdisp) == True:
                        print "folderdisp: "+folderdisp+" already exists"
                    else:
                        os.makedirs(folderdisp)
                    os.chdir(folderdisp)
                    print "alat:",alat,", atom, dips: ",atom,disp
                    cmd=str(root)+"/SphinxEmpiricalPotentials/create_position_displaced_struct.sh "+str(rootd)+"/dynmat_vak_fcc4/"+str(sc)+"x"+str(sc)+"x"+str(sc)+"sc/relaxedCoords_"+str(alat)+" "+str(atom)+" "+str(disp)+" "+str(atomsv)+"atom.sx"
                    os.system(cmd)
                    atomsx = str(atomsv)+"atom.sx"
                    get_inputsx(os.getcwd(),alat,atomsx)
                    sphinx_run()


                    os.chdir(folderalat)
            os.chdir(comeback)
    return

def get_positions_to_atomsx(pos, disp):
    print pos


path, pot, element, proj, job, mb, alats, inputsxt, mass, elname, elnamel, root, SPHINX, rootd = getinfo()
density, embeded, pair, potpot, sphinx_extractForces = density_embeded_pair(pot,root)


# create folder
for sc in scs:
    os.chdir(path)
    print "--------"
    print "sc:",sc
    print "--------"
    atoms = sc**3*4
    atomsb = atoms
    atomsv= sc**3*4-1
    folder = path+"/"+str(sc)+"x"+str(sc)+"x"+str(sc)+"sc"+str(string)
    if os.path.isdir(folder) == True:
        #shutil.rmtree(folder)
        #sys.exit("folder "+str(folder)+" exists")
        print "folder already exists"
    else:
        os.makedirs(folder)
    os.chdir(folder)
    print "-"*(len(folder)+9)
    print "folder:",folder
    print "-"*(len(folder)+9)
    #if job == "murn": create_murnjobs_imd()
    #if job == "murn": create_murnjobs_lammps()
    #details = path.split(root)[1].split("/")
    #print root+"/"+details[1]+"/"+details[2]+"/"
    print "job:",job,type(job)
    print "porg:",prog,type(prog)
    if job == "murn" and prog == "lammps": create_murnjobs_lammps()
    if job == "murn" and prog == "sphinx": create_murnjobs_sphinx()
    if job == "murn" and prog == "imd": create_murnjobs_imd()
    if job == "dynmat": create_dynmatjobs(prog)
    os.chdir(path)
