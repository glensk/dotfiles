///////////////////////////////////////////////////////////////////////////////
//
// RuNNer pair style for LAMMPS
// ----------------------------
// 
// author: Andreas Singraber
// date  : 2015-03-10
// email : andreas.singraber@univie.ac.at
//
// Copyright (c) 2013 - 2015 Andreas Singraber
//
////////////////////////////////////////////////////////////////////////////////

#include <unistd.h>
#include "math.h"
#include "stdlib.h"
#include "string.h"
#include "pair_runner.h"
#include "atom.h"
#include "force.h"
#include "domain.h"
#include "comm.h"
#include "neighbor.h"
#include "neigh_list.h"
#include "neigh_request.h"
#include "memory.h"
#include "error.h"
#include "update.h"
//#include "assert.h"

using namespace LAMMPS_NS;

////////////////////////////////////////////////////////////////////////////////

PairRuNNer::PairRuNNer(LAMMPS *lmp) : Pair(lmp) {}

////////////////////////////////////////////////////////////////////////////////

PairRuNNer::~PairRuNNer()
{
    if(allocated) {
        memory->destroy(setflag);
        memory->destroy(cutsq);
    }

    if(setup_completed == 1) {

        // L
        for(int i=0; i<ne; i++) {
            for(int j=0; j<nl+1; j++) {
                delete L[i][j];
            }
        }
        for(int i=0; i<ne; i++) {
            delete[] L[i];
        }
        delete[] L;

        // G
        for(int i=0; i<ne; i++) {
            for(int j=0; j<nsym[i]; j++) {
                delete G[i][j];
            }
        }
        for(int i=0; i<ne; i++) {
            delete[] G[i];
        }
        delete[] G;

        // GG
        for(int i=0; i<ne; i++) {
            for(int j=0; j<nsymg[i]; j++) {
                delete GG[i][j];
            }
        }
        for(int i=0; i<ne; i++) {
            delete[] GG[i];
        }
        delete[] GG;

        // a
        for(int i=0; i<na; i++) {
            delete a[i];
        }
        delete[] a;

        // e
        for(int i=0; i<ne; i++) {
            delete e[i];
        }
        delete[] e;

        // numnode, actfunc
        for(int i=0; i<ne; i++) delete[] numnode[i];
        delete[] numnode;
        for(int i=0; i<ne; i++) delete[] actfunc[i];
        delete[] actfunc;

        // nsym, nsymg
        delete[] nsym;
        delete[] nsymg;

    }

}

////////////////////////////////////////////////////////////////////////////////

void PairRuNNer::compute(int eflag, int vflag) {

    if(eflag || vflag) ev_setup(eflag, vflag);
    else evflag = vflag_fdotr = eflag_global = eflag_atom = 0;

    // allocate and set per-timestep arrays
    pt_allocset();

    // calculate neighbor list
    calc_neighbor_list();
    //calc_neighbor_list_lammps();

    // allocate space for force calculation
    init_forces();

    // calculate symmetry functions for all local atoms
    //calc_all_symfunc();
    calc_all_symfunc_group();

    // perform neural network calculation for all local atoms
    calc_nn();

    // calculate force contributions to all atoms
    calc_forces();

    // add free atom energies
    for(int i=0; i<nla; i++) {
        a[la[i]]->en += e[a[la[i]]->e]->atomic_energy;
    }

    // unit conversion for energies and forces
    for(int i=0; i<nla; i++) {
        a[la[i]]->en /= CONST_EV2HA;
    }
    for(int i=0; i<nla+atom->nghost; i++) {
        atom->f[i][0] *= CONST_A2B / CONST_EV2HA;
        atom->f[i][1] *= CONST_A2B / CONST_EV2HA;
        atom->f[i][2] *= CONST_A2B / CONST_EV2HA;
    }

    // sum up total energy
    entot = 0.0;
    for(int i=0; i<nla; i++) {
        entot += a[la[i]]->en;
    }

    // add total potential of this processor to global potential energy
    if(eflag_global) {
        ev_tally(0, 0, nla, 1, entot, 0.0, 0.0, 0.0, 0.0, 0.0);
    }

    if(eflag_atom) {
        for(int i=0; i<nla; i++) {
            eatom[i] = a[la[i]]->en;
        }
    }
    
    // if virial needed calculate via F dot r
    if(vflag_fdotr) {
        virial_fdotr_compute();
    }

    // print extrapolation warning summary if requested
    if(showewsum > 0 && update->ntimestep % showewsum == 0) {
        mpicollect_ew();
        if(comm->me == 0) {
            if(screen) {
                fprintf(screen, "# RuNNer EW summary: %ld %ld %f\n", update->ntimestep, ewcountsum, double(ewcountsum) / showewsum);
            }
            if(logfile) {
                fprintf(logfile, "# RuNNer EW summary: %ld %ld %f\n", update->ntimestep, ewcountsum, double(ewcountsum) / showewsum);
            }
        }
        ewcountsum = 0;
    }

    // delete per-timestep arrays
    pt_delete();

}

////////////////////////////////////////////////////////////////////////////////
// allocate all arrays 
////////////////////////////////////////////////////////////////////////////////

void PairRuNNer::allocate() {

    int n = atom->ntypes;
    
    allocated = 1;
    memory->create(setflag, n+1, n+1, "pair:setflag");
    for(int i=1; i<=n; i++) {
        for(int j=i; j<=n; j++) {
            setflag[i][j] = 0;
        }
    }
    
    memory->create(cutsq, n+1, n+1, "pair:cutsq");

    setup_completed = 0;

}

////////////////////////////////////////////////////////////////////////////////
// global settings 
////////////////////////////////////////////////////////////////////////////////

void PairRuNNer::settings(int narg, char **arg) {

    int iarg = 0;

    if(narg == 0) error->all(FLERR, "Illegal pair_style command");

    // default settings
    strcpy(dirRuNNer, "RuNNer");
    showew = 1;
    showewsum = 0;
    resetew = 0;
    maxew = 0;
    ntransformations = 0; //HELLSTROM


    while(iarg<narg) {
        // set RuNNer directory
        if(strcmp(arg[iarg], "dir") == 0) {
            if(iarg + 2 > narg) {
                error->all(FLERR, "Illegal pair_style command");
            }
            strcpy(dirRuNNer, arg[iarg+1]);
            iarg += 2;
        }
        // show extrapolation warnings
        else if(strcmp(arg[iarg], "showew") == 0) {
            if(iarg + 2 > narg) {
                error->all(FLERR, "Illegal pair_style command");
            }
            if(strcmp(arg[iarg+1], "yes") == 0) {
                showew = 1;
            }
            else if(strcmp(arg[iarg+1], "no") == 0) {
                showew = 0;
            }
            else {
                error->all(FLERR, "Illegal pair_style command");
            }
            iarg += 2;
        }
        // show extrapolation warning summary
        else if(strcmp(arg[iarg], "showewsum") == 0) {
            if(iarg + 2 > narg) {
                error->all(FLERR, "Illegal pair_style command");
            }
            sscanf(arg[iarg+1], "%d", &showewsum);
            iarg += 2;
        }
        // reset extrapolation warning counter
        else if(strcmp(arg[iarg], "resetew") == 0) {
            if(iarg + 2 > narg) {
                error->all(FLERR, "Illegal pair_style command");
            }
            if(strcmp(arg[iarg+1], "yes") == 0) {
                resetew = 1;
            }
            else if(strcmp(arg[iarg+1], "no") == 0) {
                resetew = 0;
            }
            else {
                error->all(FLERR, "Illegal pair_style command");
            }
            iarg += 2;
        }
        // maximum allowed extrapolation warnings
        else if(strcmp(arg[iarg], "maxew") == 0) {
            if(iarg + 2 > narg) {
                error->all(FLERR, "Illegal pair_style command");
            }
            sscanf(arg[iarg+1], "%d", &maxew);
            iarg += 2;
        }
        else if(strcmp(arg[iarg], "transform") == 0) { //HELLSTROM
            // transform 4 1 ===> atoms of type 4 will be "chemically" treated as atoms of type 1 (symmetry functions etc. are calculated as for atom type 1)
            // the benefit is that you can put more atom types in the structure file if you want several isotopes of the same element
            // it is necessary to use a different atom type because LAMMPS only allows masses per atom type (not per individual atom, other than for some specific strange atom styles)
            // the user does not need to change anything in input.nn, weights.xxx.data, scaling.data
            if(iarg + 3 > narg) {
                error->all(FLERR, "Illegal pair_style (transform) command");
            }
            ntransformations += 1;
            if (ntransformations > MAXTRANSFORMATIONS) error->all(FLERR, "Too many transformations!");
            sscanf(arg[iarg+1], "%d", &transformfrom[ntransformations-1]);
            sscanf(arg[iarg+2], "%d", &transformto[ntransformations-1]);
            iarg += 3;
        }
        else {
            error->all(FLERR, "Illegal pair_style command");
        }
    }
    sprintf(dirRuNNer, "%s/", dirRuNNer);

    // show settings
    if(comm->me == 0) {
        if(screen) {
            fprintf(screen, "RuNNer pair style directory           : %s\n", dirRuNNer);
            if     (showew == 1) fprintf(screen, "Show extrapolation warnings           : %s\n", "yes");
            else if(showew == 0) fprintf(screen, "Show extrapolation warnings           : %s\n", "no");
            if     (showewsum > 0) fprintf(screen, "Show extrapolation warning summary    : every %d timesteps\n", showewsum);
            else                   fprintf(screen, "Show extrapolation warning summary    : %s\n", "no");
            if     (resetew == 1) fprintf(screen, "Reset extrapolation warning counter   : %s\n", "yes");
            else if(resetew == 0) fprintf(screen, "Reset extrapolation warning counter   : %s\n", "no");
            fprintf(screen, "Maximum allowed extrapolation warnings: %d\n", maxew);
        }
        if(logfile) {
            fprintf(logfile, "RuNNer pair style directory           : %s\n", dirRuNNer);
            if     (showew == 1) fprintf(logfile, "Show extrapolation warnings           : %s\n", "yes");
            else if(showew == 0) fprintf(logfile, "Show extrapolation warnings           : %s\n", "no");
            if     (showewsum > 0) fprintf(logfile, "Show extrapolation warning summary    : every %d timesteps\n", showewsum);
            else                   fprintf(logfile, "Show extrapolation warning summary    : %s\n", "no");
            if     (resetew == 1) fprintf(logfile, "Reset extrapolation warning counter   : %s\n", "yes");
            else if(resetew == 0) fprintf(logfile, "Reset extrapolation warning counter   : %s\n", "no");
            fprintf(logfile, "Maximum allowed extrapolation warnings: %d\n", maxew);
        }
    }
    
}

////////////////////////////////////////////////////////////////////////////////
// set coeffs for one or more type pairs
////////////////////////////////////////////////////////////////////////////////

void PairRuNNer::coeff(int narg, char **arg) {

    if(narg > 3) error->all(FLERR, "Incorrect args for pair coefficients");
    if(!allocated) allocate();

    int ilo,ihi,jlo,jhi;
    force->bounds(FLERR, arg[0], atom->ntypes, ilo, ihi);
    force->bounds(FLERR, arg[1], atom->ntypes, jlo, jhi);

    cut_global = force->numeric(FLERR, arg[2]);

    int count = 0;
    for(int i=ilo; i<=ihi; i++) {
        for(int j=MAX(jlo,i); j<=jhi; j++) {
            setflag[i][j] = 1;
            count++;
        }
    }

    if(count == 0) error->all(FLERR, "Incorrect args for pair coefficients");

}

////////////////////////////////////////////////////////////////////////////////
// initialization of pair style
////////////////////////////////////////////////////////////////////////////////

void PairRuNNer::init_style() {

    int irequest = neighbor->request((void *) this);
    neighbor->requests[irequest]->pair = 1;
    neighbor->requests[irequest]->half = 0;
    neighbor->requests[irequest]->full = 1;

    if(setup_completed == 0) {
        if(comm->me == 0) {
            if(screen) {
                fprintf(screen, "-----------------------------------------------------------------------\n");
                fprintf(screen, "Starting RuNNer pair style setup...\n");
                fprintf(screen, "-----------------------------------------------------------------------\n");
            }
            if(logfile) {
                fprintf(logfile, "-----------------------------------------------------------------------\n");
                fprintf(logfile, "Starting RuNNer pair style setup...\n");
                fprintf(logfile, "-----------------------------------------------------------------------\n");
            }
        }

        // INITIALIZE SOME VARIABLES
        init_vars();

        // GET INFORMATION FROM input.nn
        if(comm->me == 0) {
            analyze_input_nn();
        }

        // MPI: DISTRIBUTE NN SETUP
        mpidistribute_nn();

        // ALLOCATE SPACE FOR NETWORK LAYERS
        init_layers();

        // READ NEURAL NETWORK WEIGHTS FROM weights.???.data
        if(comm->me == 0) {
            read_weights();
        }

        // MPI: DISTRIBUTE NEURAL NETWORK WEIGHTS
        mpidistribute_nnweights();

        // SORT SYMMETRY FUNCTIONS
        if(comm->me == 0) {
            sort_symfuncs();
        }

        // READ IN SCALING VALUES FROM scaling.data
        if(comm->me == 0) {
            read_scaling();
        }

        // MPI: DISTRIBUTE SYMMETRY FUNCTIONS
        mpidistribute_symfunc();

        // FIND MAXIMUM GLOBAL CUTOFF
        rcmax = rc_max();

        // CREATE SF GROUPS
        create_symfunc_groups();

        // ALLOCATE RuNNer_atom ARRAY
        na = atom->natoms;
        a = new RuNNer_atom*[na];   
        for(int i=0; i<na; i++) a[i] = new RuNNer_atom(*this);

        setup_completed = 1;

        if(comm->me == 0) {
            if(screen) {
                fprintf(screen, "RuNNer pair style setup completed.\n");
                fprintf(screen, "-----------------------------------------------------------------------\n");
            }
            if(logfile) {
                fprintf(logfile, "RuNNer pair style setup completed.\n");
                fprintf(logfile, "-----------------------------------------------------------------------\n");
            }
        } 
    }

}

////////////////////////////////////////////////////////////////////////////////
// init for one type pair i,j and corresponding j,i
////////////////////////////////////////////////////////////////////////////////

double PairRuNNer::init_one(int i, int j) {

    if(setflag[i][j] == 0) {
    }
    
    if(offset_flag) {
    } 
    else {
    }
    
    return cut_global;

}

////////////////////////////////////////////////////////////////////////////////
// proc 0 writes to restart file
////////////////////////////////////////////////////////////////////////////////

void PairRuNNer::write_restart(FILE *fp) {

    write_restart_settings(fp);
    
    int i,j;
    for(i=1; i<=atom->ntypes; i++)
        for(j=i; j<=atom->ntypes; j++) {
            fwrite(&setflag[i][j], sizeof(int), 1, fp);
            if(setflag[i][j]) {
        }
    }

}

////////////////////////////////////////////////////////////////////////////////
// proc 0 reads from restart file, bcasts
////////////////////////////////////////////////////////////////////////////////

void PairRuNNer::read_restart(FILE *fp) {

    read_restart_settings(fp);

    allocate();

    int i, j;
    int me = comm->me;
    for(i=1; i<=atom->ntypes; i++) {
        for(j=i; j<=atom->ntypes; j++) {
            if(me == 0) fread(&setflag[i][j], sizeof(int), 1, fp);
            MPI_Bcast(&setflag[i][j], 1, MPI_INT, 0, world);
            if(setflag[i][j]) {
                if(me == 0) {
                }
            }
        }
    }

}

////////////////////////////////////////////////////////////////////////////////
// proc 0 writes to restart file
////////////////////////////////////////////////////////////////////////////////

void PairRuNNer::write_restart_settings(FILE *fp) {

    fwrite(&offset_flag, sizeof(int), 1, fp);
    fwrite(&mix_flag, sizeof(int), 1, fp);

}

////////////////////////////////////////////////////////////////////////////////
// proc 0 reads from restart file, bcasts
////////////////////////////////////////////////////////////////////////////////

void PairRuNNer::read_restart_settings(FILE *fp) {

    if(comm->me == 0) {
        fread(&offset_flag, sizeof(int), 1, fp);
        fread(&mix_flag, sizeof(int), 1, fp);
    }
    MPI_Bcast(&offset_flag, 1, MPI_INT, 0, world);
    MPI_Bcast(&mix_flag, 1, MPI_INT, 0, world);

}

////////////////////////////////////////////////////////////////////////////////

double PairRuNNer::single(int i, int j, int itype, int jtype, double rsq,
    double factor_coul, double factor_lj, double &fforce) {

    error->all(FLERR, "PairRuNNer::single not possible for pair_style runner");

    return 0.0;

}


////////////////////////////////////////////////////////////////////////////////
// PairRuNNer member functions
////////////////////////////////////////////////////////////////////////////////

// initialize some variables
void PairRuNNer::init_vars() {

    setup_completed = 0;
    na = 0;
    nla = 0; 
    ne = 0;
    nl = 0;
    scalesym = 0;
    centersym = 0;
    normnodes = 0;
    ewcount = 0;
    ewcountsum = 0;
    rcmax = 0.0;
    entot = 0.0;
    scalemin = 0.0;
    scalemax = 0.0;

}

int PairRuNNer::transform_atom_type(int original_atom_type) const { //HELLSTROM
    for (int j=0; j<ntransformations; j++) {
        if (transformfrom[j] == original_atom_type) {
            //fprintf(stderr, "Transforming atom from %d to %d\n", original_atom_type, transformto[j]);
            return transformto[j];
        }
    }
    return original_atom_type;
}

// allocate and set per-timestep arrays
void PairRuNNer::pt_allocset() {

    int acount = 0;
    int index = 0;    

    // create list of local atoms
    nla = atom->nlocal;
    la = new int[nla];
    nlape = new int[ne];
    for(int i=0; i<ne; i++) nlape[i] = 0;
    for(int i=0; i<nla; i++) {
        la[i] = atom->tag[i] - 1; 
        //nlape[atom->type[i]-1]++;
        nlape[transform_atom_type(atom->type[i])-1]++; //HELLSTROM
    }

    // allocate RuNNer_neighbor array
    N = new RuNNer_neighbor*[nla];
    for(int i=0; i<nla; i++) N[i] = new RuNNer_neighbor(*this);

    // transfer local atom posisitons to RuNNer_atom array
    for(int i=0; i<nla; i++) {
        index = la[i];    
        a[index]->x = CONST_A2B * atom->x[i][0];
        a[index]->y = CONST_A2B * atom->x[i][1];
        a[index]->z = CONST_A2B * atom->x[i][2];
        //a[index]->e = atom->type[i]-1;
        a[index]->e = transform_atom_type(atom->type[i])-1; //HELLSTROM
    }

    // allocate SF value array
    laipe = new int*[ne];
    G_svalue = new double**[ne];
    for(int i=0; i<ne; i++) {
        laipe[i] = new int[nlape[i]];
        G_svalue[i] = new double*[nlape[i]];
        for(int j=0; j<nlape[i]; j++) {
            laipe[i][j] = 0;
            G_svalue[i][j] = new double[nsym[i]];
            for(int k=0; k<nsym[i]; k++) {
                G_svalue[i][j][k] = 0.0;
            }
        }
    }

    // set per-element atom index array
    for(int i=0; i<ne; i++) {
        acount = 0;
        for(int j=0; j<nla; j++) {
            if(a[la[j]]->e == i) {
                laipe[i][acount] = j;
                acount++;
            }
        }
    }

    // reset extrapolation warning counter if requested
    if(resetew == 1) {
        ewcount = 0;
    }

}

// delete per-timestep arrays
void PairRuNNer::pt_delete() {

    // free space
    for(int i=0; i<ne; i++) {
        for(int j=0; j<nlape[i]; j++) {
            delete[] G_svalue[i][j];
        }
        delete[] G_svalue[i];
        delete[] laipe[i];
    }
    delete[] G_svalue;
    delete[] laipe;
    delete[] la;
    delete[] nlape;
    for(int i=0; i<nla; i++) {
        delete N[i];
        delete F[i];
    }
    delete[] N;
    delete[] F;

}

// sort elements according to nuclear charge
void PairRuNNer::sort_elements() {

    int i;
    RuNNer_element temp(*this);

    i = 0;
    do {
        if(e[i]->nucch > e[i+1]->nucch) {
            temp = *e[i];
            *e[i] = *e[i+1];
            *e[i+1] = temp;           
            i = -1;
        }
        i++;
    } while(i<ne-1);

    return;

}

// find element number for a given element string
int PairRuNNer::find_element(char *tmpname) {

    int i;

    for(i=0; i<ne; i++) {
        if(strcmp(e[i]->name, tmpname) == 0) return i;
    }
    fprintf(stderr, "ERROR: element %s not specified!\n", tmpname);
    exit(EXIT_FAILURE);

}

// initialize force classes
void PairRuNNer::init_forces() {

    F = new RuNNer_forces*[nla];
    for(int i=0; i<nla; i++) {
        F[i] = new RuNNer_forces(*this, nl, numnode[N[i]->eatom], N[i]->eatom, i, N[i]->num);   
    }

    return;

}

// calculate forces on local and ghost atoms from saved data
void PairRuNNer::calc_forces() {
 
    int i = 0;
    int j = 0;
    int ca = 0;
    int ica = 0;
    int neigh = 0;
    int eneigh = 0;
    int ineighcaind = 0;

    // loop over all local + ghost atoms
    for(ica=0; ica<nla+atom->nghost; ica++) {
        // loop over all possible neighbors = local atoms
        for(i=0; i<nla; i++) {
            neigh = la[i];
            eneigh = a[neigh]->e;
            ineighcaind = N[i]->find_nindex(ica);
            if(ineighcaind >= 0) {
//            if(N[i]->ilal2.find(ica)!=N[i]->ilal2.end()) {
//                ineighcaind = N[i]->ilal2.at(ica);
//                //assert(ineighcaind==N[i]->find_nindex(ica));
                for(j=0; j<nsym[eneigh]; j++) {
                    atom->f[ica][0] -= F[i]->dEdG[j] * F[i]->dGdx[j][ineighcaind] * G[eneigh][j]->sf;
                    atom->f[ica][1] -= F[i]->dEdG[j] * F[i]->dGdy[j][ineighcaind] * G[eneigh][j]->sf;
                    atom->f[ica][2] -= F[i]->dEdG[j] * F[i]->dGdz[j][ineighcaind] * G[eneigh][j]->sf;
                }   
            }
        }
    }
    // add force contribution of central atom itself
    for(ica=0; ica<nla; ica++) {
        ca = la[ica];
        for(j=0; j<nsym[a[ca]->e]; j++) {
            atom->f[ica][0] -= F[ica]->dEdG[j] * F[ica]->dGdx[j][N[ica]->num] * G[a[ca]->e][j]->sf;
            atom->f[ica][1] -= F[ica]->dEdG[j] * F[ica]->dGdy[j][N[ica]->num] * G[a[ca]->e][j]->sf;
            atom->f[ica][2] -= F[ica]->dEdG[j] * F[ica]->dGdz[j][N[ica]->num] * G[a[ca]->e][j]->sf;
        }   
    }

    return;

}

// allocate space for layers
void PairRuNNer::init_layers() {

    L = new RuNNer_layer**[ne];
    for(int i=0; i<ne; i++) {
        L[i] = new RuNNer_layer*[nl+1];
    }
    for(int i=0; i<ne; i++) {
        for(int j=0; j<nl+1; j++) {
            L[i][j] = new RuNNer_layer(*this, numnode[i][j], numnode[i][j+1], actfunc[i][j+1]);
            if(comm->me==0) {
                if(screen) {
                    fprintf(screen, "Element %2d (%2s): layer %2d allocated, dim: (%3d,%3d), act. func.: %2d\n",
                        i, e[i]->name, j, numnode[i][j], numnode[i][j+1], actfunc[i][j+1]);
                }
                if(logfile) {
                    fprintf(logfile, "Element %2d (%2s): layer %2d allocated, dim: (%3d,%3d), act. func.: %2d\n",
                        i, e[i]->name, j, numnode[i][j], numnode[i][j+1], actfunc[i][j+1]);
                }
            }
        }
    }
    if(comm->me==0) {
        if(screen) fprintf(screen, "-----------------------------------------------------------------------\n");
        if(logfile) fprintf(logfile, "-----------------------------------------------------------------------\n");
    }

    return;

}

// open weights.???.data file and call get_weights
void PairRuNNer::read_weights() {

    FILE *fp;
    char fname[LLINE];
    int nweight = 0;
    int nbias = 0;

    for(int i=0; i<ne; i++) {
        sprintf(fname, "%sweights.%03d.data", dirRuNNer, e[i]->nucch);
        if(screen) fprintf(screen, "Element %2d (%2s): Reading %s...\n", i, e[i]->name, fname);
        if(logfile) fprintf(logfile, "Element %2d (%2s): Reading %s...\n", i, e[i]->name, fname);
        fp = fopen(fname, "r");
        if( fp == NULL ) {
            fprintf(stderr, "ERROR: Could not open %s.\n", fname);
            exit(EXIT_FAILURE);
        }
        else {
            nweight = 0;
            nbias = 0;
            for(int j=0; j<nl+1; j++) {
                L[i][j]->get_weights(fp);
                nweight += numnode[i][j] * numnode[i][j+1];
                nbias += numnode[i][j+1];
            }
            fclose(fp);
        }
        if(screen) {
            fprintf(screen, "Element %2d (%2s): weights: %6d, bias: %4d, total: %6d\n",
                i, e[i]->name, nweight, nbias, nweight+nbias);
        }
        if(logfile) {
            fprintf(logfile, "Element %2d (%2s): weights: %6d, bias: %4d, total: %6d\n",
                i, e[i]->name, nweight, nbias, nweight+nbias);
        }
    }
    if(screen) fprintf(screen, "-----------------------------------------------------------------------\n");
    if(logfile) fprintf(logfile, "-----------------------------------------------------------------------\n");

    return;

}

// calculate NN for all atoms
void PairRuNNer::calc_nn() {

    int ca = 0;

    // loop over all elements 
    for(int i=0; i<ne; i++) {
        // loop over all atoms of element i
        for(int ica=0; ica<nlape[i]; ica++) {
            ca = laipe[i][ica];
            // set symmetry functions values as input layer nodes
            for(int isym=0; isym<nsym[i]; isym++) {
                L[i][0]->input[isym] = G_svalue[i][ica][isym];   
            }
            // calculate all nodes of NN
            L[i][0]->calc_all_nodes();
            for(int k=1;k<nl+1;k++) {
                L[i][k]->get_input(L[i][k-1]);
                L[i][k]->calc_all_nodes();
            }
            // save data needed for force calculation
            F[ca]->get_dfdx(L[i]);
            F[ca]->calc_dEdG(L[i]);
            // atomic energy contribution is value of output node
            a[la[ca]]->en = L[i][nl]->node[0];
        }
    }

    return;

}

// calculate neighbor list for all atoms
void PairRuNNer::calc_neighbor_list() {

    int i = 0;

    // loop over all atoms
//#ifdef _OPENMP
//    #pragma omp parallel for private(i)
//#endif
    for(i=0; i<nla; i++) {
        N[i]->calc(a, la[i]);
    }

    return;

}

// calculate neighbor list for all atoms
void PairRuNNer::calc_neighbor_list_lammps() {

    int i = 0;
    int inum = list->inum;
    int* ilist = list->ilist;

    // loop over all atoms
    for(int ii=0; ii<inum; ii++) {
        i = ilist[ii];
        N[i]->calcl(a, ii);
    }

    return;

}

// Read necessary information from input.nn
void PairRuNNer::analyze_input_nn() {

    FILE *fp;
    int i, j;
    int *isf;
    int tmpe;
    int tmptype;
    char fname[LLINE];
    char line[LLINE];
    char gsfline[LLINE];
    char keyword[LLINE];
    char tmpename[3];
    int ierr;
    double tmpd;

    ne = 0;
    nl = 0;
    cutoff_type = 1;
    sprintf(fname, "%sinput.nn", dirRuNNer);
    fp = fopen(fname, "r");
    if(screen) {
        fprintf(screen, "Analyzing input.nn...\n");
        fprintf(screen, "-----------------------------------------------------------------------\n");
    }
    if(logfile) {
        fprintf(logfile, "Analyzing input.nn...\n");
        fprintf(logfile, "-----------------------------------------------------------------------\n");
    }
    if(fp == NULL) {
        fprintf(stderr, "ERROR: Could not open input.nn.\n");
        exit(EXIT_FAILURE);
    }
    else {
        while(fgets(line, LLINE, fp) != NULL) {
            keyword[0]='\0';
            ierr = sscanf(line, "%s", keyword);
            // detect number of elements
            if( (ne == 0) && (strcmp(keyword, "number_of_elements") == 0) ) {
                ierr = sscanf(line, "%*s%d", &ne);
                e = new RuNNer_element*[ne];    
                nsym = new int[ne];
                for(i=0; i<ne; i++) {
                    e[i] = new RuNNer_element(*this);
                    nsym[i] = 0;
                }
            }
            // read names of elements
            if( (ne != 0) && (strcmp(keyword, "elements") == 0) ) {
                ierr = sscanf(line, "%*s %[^\n]", line);
                for(i=0; i<ne; i++) {
                    ierr = sscanf(line, "%s %[^\n]", e[i]->name, line);
                    e[i]->getnucch();
                }
                if( ne > 1) sort_elements();
                if(screen) {
                    fprintf(screen, "Number of elements: %2d\n", ne);
                    for(i=0; i<ne; i++) fprintf(screen, "Element %2d: %2s (%2d)\n", i, e[i]->name, e[i]->nucch);
                    fprintf(screen, "-----------------------------------------------------------------------\n");
                }
                if(logfile) {
                    fprintf(logfile, "Number of elements: %2d\n", ne);
                    for(i=0; i<ne; i++) fprintf(logfile, "Element %2d: %2s (%2d)\n", i, e[i]->name, e[i]->nucch);
                    fprintf(logfile, "-----------------------------------------------------------------------\n");
                }
            }
            // read single symmetry functions
            if( (ne != 0) && (strcmp(keyword, "symfunction_short") == 0) ) {
                strcpy(tmpename, "?");
                ierr = sscanf(line, "%*s%s", tmpename);
                nsym[find_element(tmpename)]++;
            }
            // read element symmetry functions
            if( (ne != 0) && (strcmp(keyword, "element_symfunction_short") == 0) ) {
                tmpe = 0;
                strcpy(tmpename, "?");
                ierr = sscanf(line, "%*s%s%d", tmpename, &tmptype);
                tmpe = find_element(tmpename);
                if(tmptype == 2) {
                    nsym[tmpe] += ne;
                }
                else if(tmptype == 3) {
                    nsym[tmpe] += ne * (ne + 1) / 2;
                }
                else {
                    fprintf(stderr, "ERROR: symmetry function type %d not implemented.\n", tmptype);
                    exit(EXIT_FAILURE);
                }
            }
            // read global symmetry functions
            if( (ne != 0) && (strcmp(keyword, "global_symfunction_short") == 0) ) {
                ierr = sscanf(line, "%*s%d", &tmptype);
                if(tmptype == 2) {
                    for(int i=0; i<ne; i++) {
                        nsym[i] += ne;
                    }
                }
                else if(tmptype == 3) {
                    for(int i=0; i<ne; i++) {
                        nsym[i] += ne * (ne + 1) / 2;
                    }
                }
                else {
                    fprintf(stderr, "ERROR: symmetry function type %d not implemented.\n", tmptype);
                    exit(EXIT_FAILURE);
                }
            }
            // detect number of layers
            if( (ne != 0) && (strcmp(keyword, "global_hidden_layers_short") == 0) ) {
                ierr = sscanf(line, "%*s%d", &nl);
                if(screen) {
                    fprintf(screen, "Number of hidden layers: %2d\n", nl);
                    fprintf(screen, "-----------------------------------------------------------------------\n");
                }
                if(logfile) {
                    fprintf(logfile, "Number of hidden layers: %2d\n", nl);
                    fprintf(logfile, "-----------------------------------------------------------------------\n");
                }
                numnode = new int*[ne];
                actfunc = new int*[ne];
                for(i=0; i<ne; i++) {
                    numnode[i] = new int[nl+2];
                    actfunc[i] = new int[nl+2];
                    for(j=0; j<nl+2; j++) {
                        numnode[i][j]=0;
                        actfunc[i][j]=0;
                    }
                    numnode[i][nl+1]=1;
                    actfunc[i][0]=0;
                }
            }
            // detect activation function for each layer
            if( (nl != 0) && (strcmp(keyword, "global_nodes_short") == 0) ) {
                ierr = sscanf(line, "%*s %[^\n]", line);
                for(i=1; i<nl; i++) {
                    ierr = sscanf(line, "%d %[^\n]", &numnode[0][i], line);
                    for(j=1; j<ne; j++) numnode[j][i] = numnode[0][i];
                }
                ierr = sscanf(line, "%d", &numnode[0][nl]);
                for(j=1; j<ne; j++) numnode[j][nl] = numnode[0][nl];
                for(i=0; i<ne; i++) {
                    for(j=1; j<nl+2; j++) {
                        if(screen) fprintf(screen, "Element %2d (%2s): layer %2d: %3d nodes\n", i, e[i]->name, j-1, numnode[i][j]);
                        if(logfile) fprintf(logfile, "Element %2d (%2s): layer %2d: %3d nodes\n", i, e[i]->name, j-1, numnode[i][j]);
                    }
                }
            }
            // detect number of nodes of each layer
            if( (nl != 0) && (strcmp(keyword, "global_activation_short") == 0) ) {
                char tmpact[2];
                ierr = sscanf(line, "%*s %[^\n]", line);
                for(i=1; i<nl+1; i++) {
                    ierr = sscanf(line, "%s %[^\n]", tmpact, line);
                    if(strcmp(tmpact, "l") == 0) actfunc[0][i]=1;
                    if(strcmp(tmpact, "t") == 0) actfunc[0][i]=2;
                    for(j=1; j<ne; j++) actfunc[j][i] = actfunc[0][i];
                }
                ierr = sscanf(line, "%s", tmpact);
                if(strcmp(tmpact, "l") == 0) actfunc[0][nl+1] = 1;
                if(strcmp(tmpact, "t") == 0) actfunc[0][nl+1] = 2;
                for(j=1; j<ne; j++) actfunc[j][nl+1] = actfunc[0][nl+1];
                for(i=0; i<ne; i++) {
                    for(j=1; j<nl+2; j++) {
                        if(screen) {
                            fprintf(screen, "Element %2d (%2s): layer %2d: activation function: %2d\n",
                                i, e[i]->name, j-1, actfunc[i][j]);
                        }
                        if(logfile) {
                            fprintf(logfile, "Element %2d (%2s): layer %2d: activation function: %2d\n",
                                i, e[i]->name, j-1, actfunc[i][j]);
                        }
                    }
                }
                if(screen) fprintf(screen, "-----------------------------------------------------------------------\n");
                if(logfile) fprintf(logfile, "-----------------------------------------------------------------------\n");
            }
            // detect scale_symmetry_functions
            if( strcmp(keyword, "scale_symmetry_functions") == 0 ) {
                scalesym = 1;
                if(screen) {
                    fprintf(screen, "scale_symmetry_functions detected.\n");
                    fprintf(screen, "-----------------------------------------------------------------------\n");
                }
                if(logfile) {
                    fprintf(logfile, "scale_symmetry_functions detected.\n");
                    fprintf(logfile, "-----------------------------------------------------------------------\n");
                }
            }
            // detect center_symmetry_functions
            if( strcmp(keyword, "center_symmetry_functions") == 0 ) {
                centersym = 1;
                if(screen) {
                    fprintf(screen, "center_symmetry_functions detected.\n");
                    fprintf(screen, "-----------------------------------------------------------------------\n");
                }
                if(logfile) {
                    fprintf(logfile, "center_symmetry_functions detected.\n");
                    fprintf(logfile, "-----------------------------------------------------------------------\n");
                }
            }
            // set scalemin value
            if( strcmp(keyword, "scale_min_short") == 0 ) {
                ierr = sscanf(line, "%*s%lf", &scalemin);
                if(screen) {
                    fprintf(screen, "scalemin = %f\n", scalemin);
                    fprintf(screen, "-----------------------------------------------------------------------\n");
                }
                if(logfile) {
                    fprintf(logfile, "scalemin = %f\n", scalemin);
                    fprintf(logfile, "-----------------------------------------------------------------------\n");
                }
            }
            // set scalemax value
            if( strcmp(keyword, "scale_max_short") == 0 ) {
                ierr = sscanf(line, "%*s%lf", &scalemax);
                if(screen) {
                    fprintf(screen, "scalemax = %f\n", scalemax);
                    fprintf(screen, "-----------------------------------------------------------------------\n");
                }
                if(logfile) {
                    fprintf(logfile, "scalemax = %f\n", scalemax);
                    fprintf(logfile, "-----------------------------------------------------------------------\n");
                }
            }
            // detect normalize_nodes
            if( strcmp(keyword, "normalize_nodes") == 0 ) {
                normnodes = 1;
                if(screen) {
                    fprintf(screen, "normalize_nodes detected.\n");
                    fprintf(screen, "-----------------------------------------------------------------------\n");
                }
                if(logfile) {
                    fprintf(logfile, "normalize_nodes detected.\n");
                    fprintf(logfile, "-----------------------------------------------------------------------\n");
                }
            }
            // detect cutoff type
            if( (strcmp(keyword, "cutoff_type") == 0) ) {
                ierr = sscanf(line, "%*s%d", &cutoff_type);
                if(screen) {
                    fprintf(screen, "cutoff type = %d\n", cutoff_type);
                    fprintf(screen, "-----------------------------------------------------------------------\n");
                }
                if(logfile) {
                    fprintf(logfile, "cutoff type = %d\n", cutoff_type);
                    fprintf(logfile, "-----------------------------------------------------------------------\n");
                }
            }
            // read atomic energies
            if( (strcmp(keyword, "atom_energy") == 0) ) {
                ierr = sscanf(line, "%*s%s%lf", tmpename, &tmpd);
                tmpe = find_element(tmpename);
                e[tmpe]->atomic_energy = tmpd;
                if(screen) {
                    fprintf(screen, "Free atom energy of element %2d (%2s) = %16.8E\n", tmpe, e[tmpe]->name, e[tmpe]->atomic_energy);
                    fprintf(screen, "-----------------------------------------------------------------------\n");
                }
                if(logfile) {
                    fprintf(logfile, "Free atom energy of element %2d (%2s) = %16.8E\n", tmpe, e[tmpe]->name, e[tmpe]->atomic_energy);
                    fprintf(logfile, "-----------------------------------------------------------------------\n");
                }
            }
        } 
        if(screen) {
            for(i=0; i<ne; i++) fprintf(screen, "Element %2d (%2s): %d symmetry functions\n", i, e[i]->name, nsym[i]);
            fprintf(screen, "-----------------------------------------------------------------------\n");
        }
        if(logfile) {
            for(i=0; i<ne; i++) fprintf(logfile, "Element %2d (%2s): %d symmetry functions\n", i, e[i]->name, nsym[i]);
            fprintf(logfile, "-----------------------------------------------------------------------\n");
        }
        G = new RuNNer_symfunc**[ne];   
        for(i=0; i<ne; i++) {
            G[i] = new RuNNer_symfunc*[nsym[i]];    
            for(j=0; j<nsym[i]; j++) G[i][j] = new RuNNer_symfunc(*this);
            numnode[i][0] = nsym[i];
        }
        // read input.nn again and set symmetry function parameters
        rewind(fp);
        isf = new int[ne];
        for(i=0; i<ne; i++) isf[i] = 0;
        while(fgets(line, LLINE, fp) != NULL) {
            keyword[0] = '\0';
            ierr = sscanf(line, "%s", keyword);
            if( (ne != 0) && (strcmp(keyword, "symfunction_short") == 0) ) {
                tmpe = 0;
                strcpy(tmpename, "?");
                ierr = sscanf(line, "%*s%s", tmpename);
                tmpe = find_element(tmpename);
                ierr = sscanf(line, "%*s %[^\n]", line);
                G[tmpe][isf[tmpe]]->setup(line);
                isf[tmpe]++;
            }
            if( (ne != 0) && (strcmp(keyword, "element_symfunction_short") == 0) ) {
                tmpe = 0;
                strcpy(tmpename, "?");
                ierr = sscanf(line, "%*s%s%d", tmpename, &tmptype);
                tmpe = find_element(tmpename);
                ierr = sscanf(line, "%*s %*s %*d %[^\n]", line);
                if(tmptype == 2) {
                    for(int i=0; i<ne; i++) {
                        sprintf(gsfline, "%s %d %s %s", tmpename, tmptype, e[i]->name, line);
                        G[tmpe][isf[tmpe]]->setup(gsfline);
                        isf[tmpe]++;
                    }
                }
                else if(tmptype == 3) {
                    for(int i=0; i<ne; i++) {
                        for(int j=i; j<ne; j++) {
                            sprintf(gsfline, "%s %d %s %s %s", tmpename, tmptype, e[i]->name, e[j]->name, line);
                            G[tmpe][isf[tmpe]]->setup(gsfline);
                            isf[tmpe]++;
                        } 
                    }
                }
            }
            if( (ne != 0) && (strcmp(keyword, "global_symfunction_short") == 0) ) {
                tmpe = 0;
                ierr = sscanf(line, "%*s%d", &tmptype);
                ierr = sscanf(line, "%*s %*d %[^\n]", line);
                if(tmptype == 2) {
                    for(int i=0; i<ne; i++) {
                        for(int j=0; j<ne; j++) {
                            sprintf(gsfline, "%s %d %s %s", e[i]->name, tmptype, e[j]->name, line);
                            G[i][isf[i]]->setup(gsfline);
                            isf[i]++;
                        }
                    }
                }
                else if(tmptype == 3) {
                    for(int i=0; i<ne; i++) {
                        for(int j=0; j<ne; j++) {
                            for(int k=j; k<ne; k++) {
                                sprintf(gsfline, "%s %d %s %s %s", e[i]->name, tmptype, e[j]->name, e[k]->name, line);
                                G[i][isf[i]]->setup(gsfline);
                                isf[i]++;
                            } 
                        }
                    }
                }
            }
        } 
        delete[] isf;
        fclose(fp);
    }
    if(screen) {
        fprintf(screen, "Finished analyzing input.nn.\n");
        fprintf(screen, "-----------------------------------------------------------------------\n");
    }
    if(logfile) {
        fprintf(logfile, "Finished analyzing input.nn.\n");
        fprintf(logfile, "-----------------------------------------------------------------------\n");
    }

    return;

}

// sort symmetry functions
void PairRuNNer::sort_symfuncs() {

    int i, j;
    RuNNer_symfunc temp(*this);

    if(screen) fprintf(screen, "Sorting according to   type.\n");
    if(logfile) fprintf(logfile, "Sorting according to   type.\n");
    for(i=0; i<ne; i++) {
        if(nsym[i]>1) {
            j = 0;
            do {
                if(G[i][j]->type > G[i][j+1]->type) {
                    temp = *G[i][j];
                    *G[i][j] = *G[i][j+1];
                    *G[i][j+1] = temp;
                    j = -1;
                }
                j++;
            } while(j<nsym[i]-1);
        }
    }

    if(screen) fprintf(screen, "Sorting according to     rc.\n");
    if(logfile) fprintf(logfile, "Sorting according to     rc.\n");
    for(i=0; i<ne; i++) {
        if(nsym[i]>1) {
            j = 0;
            do {
                if( (G[i][j]->type == G[i][j+1]->type) &&
                    (G[i][j]->rc > G[i][j+1]->rc) ) {
                    temp = *G[i][j];
                    *G[i][j] = *G[i][j+1];
                    *G[i][j+1] = temp;
                    j = -1;
                }
                j++;
            } while(j<nsym[i]-1);
        }
    }

    if(screen) fprintf(screen, "Sorting according to    eta.\n");
    if(logfile) fprintf(logfile, "Sorting according to    eta.\n");
    for(i=0; i<ne; i++) {
        if(nsym[i]>1) {
            j = 0;
            do {
                if( (G[i][j]->type == G[i][j+1]->type) &&
                    (G[i][j]->rc == G[i][j+1]->rc) &&
                    (G[i][j]->eta > G[i][j+1]->eta) ) {
                    temp = *G[i][j];
                    *G[i][j] = *G[i][j+1];
                    *G[i][j+1] = temp;
                    j = -1;
                }
                j++;
            } while(j<nsym[i]-1);
        }
    }

    if(screen) fprintf(screen, "Sorting according to   zeta.\n");
    if(logfile) fprintf(logfile, "Sorting according to   zeta.\n");
    for(i=0; i<ne; i++) {
        if(nsym[i]>1) {
            j = 0;
            do {
                if( (G[i][j]->type == G[i][j+1]->type) &&
                    (G[i][j]->rc == G[i][j+1]->rc) &&
                    (G[i][j]->eta == G[i][j+1]->eta) && 
                    (G[i][j]->zeta > G[i][j+1]->zeta) ) {
                    temp = *G[i][j];
                    *G[i][j] = *G[i][j+1];
                    *G[i][j+1] = temp;
                    j = -1;
                }
                j++;
            } while(j<nsym[i]-1);
        }
    }

    if(screen) fprintf(screen, "Sorting according to lambda.\n");
    if(logfile) fprintf(logfile, "Sorting according to lambda.\n");
    for(i=0; i<ne; i++) {
        if(nsym[i]>1) {
            j = 0;
            do {
                if( (G[i][j]->type == G[i][j+1]->type) &&
                    (G[i][j]->rc == G[i][j+1]->rc) &&
                    (G[i][j]->eta == G[i][j+1]->eta) && 
                    (G[i][j]->zeta == G[i][j+1]->zeta) &&
                    (G[i][j]->lambda > G[i][j+1]->lambda) ) {
                    temp = *G[i][j];
                    *G[i][j] = *G[i][j+1];
                    *G[i][j+1] = temp;
                    j = -1;
                }
                j++;
            } while(j<nsym[i]-1);
        }
    }

    if(screen) fprintf(screen, "Sorting according to     rs.\n");
    if(logfile) fprintf(logfile, "Sorting according to     rs.\n");
    for(i=0; i<ne; i++) {
        if(nsym[i]>1) {
            j = 0;
            do {
                if( (G[i][j]->type == G[i][j+1]->type) &&
                    (G[i][j]->rc == G[i][j+1]->rc) &&
                    (G[i][j]->eta == G[i][j+1]->eta) && 
                    (G[i][j]->zeta == G[i][j+1]->zeta) &&
                    (G[i][j]->lambda == G[i][j+1]->lambda) &&
                    (G[i][j]->rs > G[i][j+1]->rs) ) {
                    temp = *G[i][j];
                    *G[i][j] = *G[i][j+1];
                    *G[i][j+1] = temp;
                    j = -1;
                }
                j++;
            } while(j<nsym[i]-1);
        }
    }

    if(screen) fprintf(screen, "Sorting according to     e1.\n");
    if(logfile) fprintf(logfile, "Sorting according to     e1.\n");
    for(i=0; i<ne; i++) {
        if(nsym[i]>1) {
            j = 0;
            do {
                if( (G[i][j]->type == G[i][j+1]->type) &&
                    (G[i][j]->rc == G[i][j+1]->rc) &&
                    (G[i][j]->eta == G[i][j+1]->eta) && 
                    (G[i][j]->zeta == G[i][j+1]->zeta) &&
                    (G[i][j]->lambda == G[i][j+1]->lambda) &&
                    (G[i][j]->rs == G[i][j+1]->rs) &&
                    (G[i][j]->e1 > G[i][j+1]->e1) ) {
                    temp = *G[i][j];
                    *G[i][j] = *G[i][j+1];
                    *G[i][j+1] = temp;
                    j = -1;
                }
                j++;
            } while(j<nsym[i]-1);
        }
    }

    if(screen) fprintf(screen, "Sorting according to     e2.\n");
    if(logfile) fprintf(logfile, "Sorting according to     e2.\n");
    for(i=0; i<ne; i++) {
        if(nsym[i]>1) {
            j = 0;
            do {
                if( (G[i][j]->type == G[i][j+1]->type) &&
                    (G[i][j]->rc == G[i][j+1]->rc) &&
                    (G[i][j]->eta == G[i][j+1]->eta) && 
                    (G[i][j]->zeta == G[i][j+1]->zeta) &&
                    (G[i][j]->lambda == G[i][j+1]->lambda) &&
                    (G[i][j]->rs == G[i][j+1]->rs) &&
                    (G[i][j]->e1 == G[i][j+1]->e1) &&
                    (G[i][j]->e2 > G[i][j+1]->e2) ) {
                    temp = *G[i][j];
                    *G[i][j] = *G[i][j+1];
                    *G[i][j+1] = temp;
                    j = -1;
                }
                j++;
            } while(j<nsym[i]-1);
        }
    }

    if(screen) {
        fprintf(screen, "Printing symmetry functions...\n");
        fprintf(screen, "-----------------------------------------------------------------------\n");
    }
    if(logfile) {
        fprintf(logfile, "Printing symmetry functions...\n");
        fprintf(logfile, "-----------------------------------------------------------------------\n");
    }
    for(i=0; i<ne; i++) {
        if(screen) {
            fprintf(screen, "short range atomic symmetry functions element %2s\n", e[i]->name);
            fprintf(screen, "-----------------------------------------------------------------------\n");
        }
        if(logfile) {
            fprintf(logfile, "short range atomic symmetry functions element %2s\n", e[i]->name);
            fprintf(logfile, "-----------------------------------------------------------------------\n");
        }
        for(j=0; j<nsym[i]; j++) {
            G[i][j]->num = j;
            G[i][j]->info();
        }
        if(screen) fprintf(screen, "-----------------------------------------------------------------------\n");
        if(logfile) fprintf(logfile, "-----------------------------------------------------------------------\n");
    }

    return;

}

// read scaling.data and extract minimum, maximum and average SF values
void PairRuNNer::read_scaling() {

    FILE *fp;
    char fname[LLINE];
    char line[LLINE];
    char *cerr;
    int ierr;
    int tmpi, tmpj;
    
    sprintf(fname, "%sscaling.data", dirRuNNer);
    fp = fopen(fname, "r");
    if(screen) fprintf(screen, "Reading scaling.data...\n");
    if(logfile) fprintf(logfile, "Reading scaling.data...\n");
    if(fp == NULL) {
        fprintf(stderr, "ERROR: Could not open scaling.data.\n");
        exit(EXIT_FAILURE);
    }
    else {
        for(int i=0; i<ne; i++) {
            for(int j=0; j<nsym[i]; j++) {
                cerr = fgets(line, LLINE, fp);
                ierr = sscanf(line, "%d %d %[^\n]", &tmpi, &tmpj, line);
                tmpi--;
                tmpj--;
                // save the scaling factor for later usage
                if(tmpi==i && tmpj==j) {
                    ierr = sscanf(line, "%lf%lf%lf", &G[i][j]->min, &G[i][j]->max, &G[i][j]->ave);
                    if     (scalesym == 1 && centersym == 0) G[i][j]->sf = (scalemax - scalemin) / (G[i][j]->max - G[i][j]->min);
                    else if(scalesym == 1 && centersym == 1) G[i][j]->sf = 1.0 / (G[i][j]->max - G[i][j]->min);
                    else G[i][j]->sf = 1.0;
                }
                else {
                    fprintf(stderr, "ERROR: scaling.data inconsistent with input.nn.\n");
                    exit(EXIT_FAILURE);
                }
            }
        }
        fclose(fp);
    }
    if(screen) {
        fprintf(screen, "Finished reading scaling.data.\n");
        fprintf(screen, "-----------------------------------------------------------------------\n");
    }
    if(logfile) {
        fprintf(logfile, "Finished reading scaling.data.\n");
        fprintf(logfile, "-----------------------------------------------------------------------\n");
    }

    return;

}

// create symmetry function groups
void PairRuNNer::create_symfunc_groups() {

    int     kstart    = 0;
    int**   nsympg    = NULL;
    int*    tmptype   = NULL;
    int*    tmpe1     = NULL;
    int*    tmpe2     = NULL;
    double* tmprc     = NULL;
    double* tmpeta    = NULL;
    double* tmpzeta   = NULL;
    double* tmplambda = NULL;
    double* tmprs     = NULL;

    if(comm->me == 0) {
        if(screen) {
            fprintf(screen, "Grouping symmetry functions...\n");
        }
        if(logfile) {
            fprintf(logfile, "Grouping symmetry functions...\n");
        }
    }

    nsymg = new int[ne];
    nsympg = new int*[ne];
    GG = new RuNNer_symfuncGroup**[ne];
    for(int i=0; i<ne; i++) {
        nsymg[i] = 0;
        nsympg[i] = new int[nsym[i]];
        for(int j=0; j<nsym[i]; j++) {
            nsympg[i][j] = 0;
        }
    }

    // detect number of symmetry function groups (assume symmetry functions are already sorted!)
    for(int i=0; i<ne; i++) {

        tmptype   = new int[nsym[i]];
        tmprc     = new double[nsym[i]];
        tmpeta    = new double[nsym[i]];
        tmpzeta   = new double[nsym[i]];
        tmplambda = new double[nsym[i]];
        tmprs     = new double[nsym[i]];
        tmpe1     = new int[nsym[i]];
        tmpe2     = new int[nsym[i]];
        for(int j=0; j<nsym[i]; j++) {
            tmptype[j]   = 0;
            tmprc[j]     = 0.0;
            tmpeta[j]    = 0.0;
            tmpzeta[j]   = 0.0;
            tmplambda[j] = 0.0;
            tmprs[j]     = 0.0;
            tmpe1[j]     = 0;
            tmpe2[j]     = 0;
        }


        for(int j=0; j<nsym[i]; j++) {
            if(G[i][j]->type == 2) {
                if(j == 0) {
                    tmptype[0] = G[i][j]->type;
                    tmprc  [0] = G[i][j]->rc;
                    tmpe1  [0] = G[i][j]->e1;
                    nsymg[i]++;
                    nsympg[i][0]++;
                    continue;
                }
                for(int k=0; k<nsymg[i]; k++) {
                    if((G[i][j]->type == tmptype[k]) && 
                       (G[i][j]->rc   == tmprc  [k]) && 
                       (G[i][j]->e1   == tmpe1  [k]) ) {
                        nsympg[i][k]++;
                        break;
                    }
                    if(k == nsymg[i] - 1) {
                        tmptype[nsymg[i]] = G[i][j]->type;
                        tmprc  [nsymg[i]] = G[i][j]->rc  ;
                        tmpe1  [nsymg[i]] = G[i][j]->e1  ;
                        nsympg[i][nsymg[i]]++;
                        nsymg[i]++;
                        break;
                    }
                }
            }
            else if(G[i][j]->type == 3) {
                if(tmptype[nsymg[i] - 1] == 2) {
                    kstart = nsymg[i];
                    tmptype[kstart] = G[i][j]->type;
                    tmprc  [kstart] = G[i][j]->rc;
                    tmpe1  [kstart] = G[i][j]->e1;
                    tmpe2  [kstart] = G[i][j]->e2;
                    nsymg[i]++;
                    nsympg[i][kstart]++;
                    continue;
                }
                for(int k=kstart; k<nsymg[i]; k++) {
                    if((G[i][j]->type == tmptype[k]) && 
                       (G[i][j]->rc   == tmprc  [k]) && 
                       (G[i][j]->e1   == tmpe1  [k]) &&
                       (G[i][j]->e2   == tmpe2  [k]) ) {
                        nsympg[i][k]++;
                        break;
                    }
                    if(k == nsymg[i] - 1) {
                        tmptype[nsymg[i]] = G[i][j]->type;
                        tmprc  [nsymg[i]] = G[i][j]->rc  ;
                        tmpe1  [nsymg[i]] = G[i][j]->e1  ;
                        tmpe2  [nsymg[i]] = G[i][j]->e2  ;
                        nsympg[i][nsymg[i]]++;
                        nsymg[i]++;
                        break;
                    }
                }
            }
        }

        // allocate symmetry function group objects
        GG[i] = new RuNNer_symfuncGroup*[nsymg[i]];
        for(int j=0; j<nsymg[i]; j++) GG[i][j] = new RuNNer_symfuncGroup(*this, j, nsympg[i][j]);

        // loop over symmetry functions again and setup group objects
        for(int j=0; j<nsym[i]; j++) {
            tmptype[j]   = 0;
            tmprc[j]     = 0.0;
            tmpeta[j]    = 0.0;
            tmpzeta[j]   = 0.0;
            tmplambda[j] = 0.0;
            tmprs[j]     = 0.0;
            tmpe1[j]     = 0;
            tmpe2[j]     = 0;
        }

        nsymg[i] = 0;
        delete[] nsympg[i];
        nsympg[i] = new int[nsym[i]];
        for(int j=0; j<nsym[i]; j++) {
            nsympg[i][j] = 0;
        }


        for(int j=0; j<nsym[i]; j++) {
            if(G[i][j]->type == 2) {
                if(j == 0) {
                    tmptype[0] = G[i][j]->type;
                    tmprc  [0] = G[i][j]->rc;
                    tmpe1  [0] = G[i][j]->e1;
                    GG[i][0]->setup(tmptype[0], i, tmpe1[0], tmpe2[0], tmprc[0]);
                    GG[i][0]->setsc(0, G[i][j]->num, G[i][j]->min, G[i][j]->max, G[i][j]->ave, G[i][j]->sf);
                    GG[i][0]->settype2(0, G[i][j]->eta, G[i][j]->rs);
                    nsymg[i]++;
                    nsympg[i][0]++;
                    continue;
                }
                for(int k=0; k<nsymg[i]; k++) {
                    if((G[i][j]->type == tmptype[k]) && 
                       (G[i][j]->rc   == tmprc  [k]) && 
                       (G[i][j]->e1   == tmpe1  [k]) ) {
                        GG[i][k]->setsc(nsympg[i][k], G[i][j]->num, G[i][j]->min, G[i][j]->max, G[i][j]->ave, G[i][j]->sf);
                        GG[i][k]->settype2(nsympg[i][k], G[i][j]->eta, G[i][j]->rs);
                        nsympg[i][k]++;
                        break;
                    }
                    if(k == nsymg[i] - 1) {
                        tmptype[nsymg[i]] = G[i][j]->type;
                        tmprc  [nsymg[i]] = G[i][j]->rc  ;
                        tmpe1  [nsymg[i]] = G[i][j]->e1  ;
                        GG[i][nsymg[i]]->setup(tmptype[nsymg[i]], i, tmpe1[nsymg[i]], tmpe2[nsymg[i]], tmprc[nsymg[i]]);
                        GG[i][nsymg[i]]->setsc(nsympg[i][nsymg[i]], G[i][j]->num, G[i][j]->min, G[i][j]->max, G[i][j]->ave, G[i][j]->sf);
                        GG[i][nsymg[i]]->settype2(nsympg[i][nsymg[i]], G[i][j]->eta, G[i][j]->rs);
                        nsympg[i][nsymg[i]]++;
                        nsymg[i]++;
                        break;
                    }
                }
            }
            else if(G[i][j]->type == 3) {
                if(tmptype[nsymg[i] - 1] == 2) {
                    kstart = nsymg[i];
                    tmptype[kstart] = G[i][j]->type;
                    tmprc  [kstart] = G[i][j]->rc;
                    tmpe1  [kstart] = G[i][j]->e1;
                    tmpe2  [kstart] = G[i][j]->e2;
                    GG[i][kstart]->setup(tmptype[kstart], i, tmpe1[kstart], tmpe2[kstart], tmprc[kstart]);
                    GG[i][kstart]->setsc(nsympg[i][kstart], G[i][j]->num, G[i][j]->min, G[i][j]->max, G[i][j]->ave, G[i][j]->sf);
                    GG[i][kstart]->settype3(nsympg[i][kstart], G[i][j]->eta, G[i][j]->zeta, G[i][j]->lambda, G[i][j]->optpow);
                    nsymg[i]++;
                    nsympg[i][kstart]++;
                    continue;
                }
                for(int k=kstart; k<nsymg[i]; k++) {
                    if((G[i][j]->type == tmptype[k]) && 
                       (G[i][j]->rc   == tmprc  [k]) && 
                       (G[i][j]->e1   == tmpe1  [k]) &&
                       (G[i][j]->e2   == tmpe2  [k]) ) {
                        GG[i][k]->setsc(nsympg[i][k], G[i][j]->num, G[i][j]->min, G[i][j]->max, G[i][j]->ave, G[i][j]->sf);
                        GG[i][k]->settype3(nsympg[i][k], G[i][j]->eta, G[i][j]->zeta, G[i][j]->lambda, G[i][j]->optpow);
                        nsympg[i][k]++;
                        break;
                    }
                    if(k == nsymg[i] - 1) {
                        tmptype[nsymg[i]] = G[i][j]->type;
                        tmprc  [nsymg[i]] = G[i][j]->rc  ;
                        tmpe1  [nsymg[i]] = G[i][j]->e1  ;
                        tmpe2  [nsymg[i]] = G[i][j]->e2  ;
                        GG[i][nsymg[i]]->setup(tmptype[nsymg[i]], i, tmpe1[nsymg[i]], tmpe2[nsymg[i]], tmprc[nsymg[i]]);
                        GG[i][nsymg[i]]->setsc(nsympg[i][nsymg[i]], G[i][j]->num, G[i][j]->min, G[i][j]->max, G[i][j]->ave, G[i][j]->sf);
                        GG[i][nsymg[i]]->settype3(nsympg[i][nsymg[i]], G[i][j]->eta, G[i][j]->zeta, G[i][j]->lambda, G[i][j]->optpow);
                        nsympg[i][nsymg[i]]++;
                        nsymg[i]++;
                        break;
                    }
                }
            }
        }

        for(int j=0; j<nsymg[i]; j++) GG[i][j]->postset();

        delete[] tmptype  ;
        delete[] tmprc    ;
        delete[] tmpeta   ;
        delete[] tmpzeta  ;
        delete[] tmplambda;
        delete[] tmprs    ;
        delete[] tmpe1    ;
        delete[] tmpe2    ;

    }

    if(comm->me == 0) {
        if(screen) {
            fprintf(screen, "Printing symmetry function groups...\n");
            fprintf(screen, "-----------------------------------------------------------------------\n");
        }
        if(logfile) {
            fprintf(logfile, "Printing symmetry function groups...\n");
            fprintf(logfile, "-----------------------------------------------------------------------\n");
        }

        for(int i=0; i<ne; i++) {
            if(screen) {
                fprintf(screen, "short range atomic symmetry functions groups of element %2s\n", e[i]->name);
                fprintf(screen, "-----------------------------------------------------------------------\n");
            }
            if(logfile) {
                fprintf(logfile, "short range atomic symmetry functions groups of element %2s\n", e[i]->name);
                fprintf(logfile, "-----------------------------------------------------------------------\n");
            }
            for(int j=0; j<nsymg[i]; j++) {
                GG[i][j]->info();
            }
            if(screen) {
                fprintf(screen, "-----------------------------------------------------------------------\n");
            }
            if(logfile) {
                fprintf(logfile, "-----------------------------------------------------------------------\n");
            }           
        }
    }


    for(int i=0; i<ne; i++) delete[] nsympg[i];
    delete[] nsympg;

    return;

}

// calculate the maximum global cutoff for all symmetry functions
double PairRuNNer::rc_max() {
    
    double rcmax = 0.0;

    for(int i=0; i<ne; i++) {
        for(int j=0; j<nsym[i]; j++) {
            rcmax = MAX(rcmax, G[i][j]->rc);
        }
    }
    if(comm->me==0) {
        if(screen) {
            fprintf(screen, "rcmax = %f\n", rcmax);
            fprintf(screen, "-----------------------------------------------------------------------\n");
        }
        if(logfile) {
            fprintf(logfile, "rcmax = %f\n", rcmax);
            fprintf(logfile, "-----------------------------------------------------------------------\n");
        }
    }

    return rcmax;

}

// soft cutoff function
inline double PairRuNNer::fc(double r, double rcinv) {

    if(cutoff_type == 1) {
        return 0.5 * (cos(M_PI * r * rcinv) + 1.0);
    }
    else if(cutoff_type == 2) {
        double temp = tanh(1.0 - r * rcinv);
        return temp * temp * temp;
    }

    return -1;

}

// soft cutoff function derivative
inline double PairRuNNer::dfc(double r, double rcinv) {

    if(cutoff_type == 1) {
        return -0.5 * M_PI * sin(M_PI * r * rcinv) * rcinv;
    }
    else if(cutoff_type == 2) {
        double temp = tanh(1.0 - r * rcinv);
        temp *= temp;
        return 3.0 * temp * (temp - 1.0) * rcinv;
    }

    return -1;

}

// optimized power function for integers (gsl_pow_int)
inline double PairRuNNer::powint(double x, unsigned int n) {

    double value = 1.0;

    /* repeated squaring method 
     * returns 0.0^0 = 1.0, so continuous in x
     */
    do {
        if(n & 1) value *= x;  /* for n odd */
        n >>= 1;
        x *= x;
    } while(n);
  
    return value;
}


// calculate one symmetry function
double PairRuNNer::calc_G_dG(int ca, RuNNer_forces *F, RuNNer_symfunc *G, RuNNer_neighbor *N) {

    int num    = G->num ;
    int type   = G->type;
    //int ec     = G->ec  ;
    int e1     = G->e1  ;
    int e2     = G->e2  ;
    int optpow = G->optpow;
    double value  = 0.0     ;
    double svalue = 0.0     ;
    double rc     = G->rc    ;
    double eta    = G->eta   ;
    double zeta   = G->zeta  ;
    double lambda = G->lambda;
    double rs     = G->rs    ;
    double min    = G->min   ;
    double max    = G->max   ;
    double ave    = G->ave   ;
    double sf     = G->sf    ;

    int Nej, Nek;
    int Nnum;
    double rcinv = 1.0 / rc;
    double rc2 = rc * rc;
    double rij, rik, rjk, rinvijik;
    double r2ij, r2ik;
    double dxij, dyij, dzij;
    double dxik, dyik, dzik;
    double dxjk, dyjk, dzjk;
    double costijk;
    double pexp, plambda, pfc, pdfc, pnorm, pzl;
    double pfcij, pfcik, pfcjk;
    double pdfcij;
    double p1, p2, p3;
    double p2etaplambda;
    double fg;
    double *xptr, *yptr, *zptr;

    // pointer to force class
    xptr = F->dGdx[num];
    yptr = F->dGdy[num];
    zptr = F->dGdz[num];

    // type 2 symmetry function
    if(type==2) {
        Nnum = N->num;
        for(int j=0; j<Nnum; j++) {
            Nej = N->e[j];
            rij = N->d[j];
            if(rij<=rc && e1 == Nej) {  
                // energy calculation
                dxij = N->dx[j];
                dyij = N->dy[j];
                dzij = N->dz[j];
                pexp = exp(-eta * (rij - rs) * (rij - rs));
                pfc = fc(rij, rcinv);
                value += pexp * pfc;
                // force calculation
                pdfc = dfc(rij, rcinv);
                p1 = (pdfc - 2.0 * eta * (rij - rs) * pfc) * pexp / rij;
                dxij = p1 * dxij;
                dyij = p1 * dyij;
                dzij = p1 * dzij;
                // save force contributions in force class
                xptr[Nnum] += dxij;
                yptr[Nnum] += dyij;
                zptr[Nnum] += dzij;
                xptr[j] -= dxij;
                yptr[j] -= dyij;
                zptr[j] -= dzij;
            }
        }
    }

    // type 3 symmetry function
    if(type==3) {
        Nnum = N->num;
        pnorm = pow(2.0, 1.0 - zeta);
        pzl = zeta * lambda;
        for(int j=0; j<Nnum-1; j++) {
            Nej = N->e[j];
            rij = N->d[j];
            if( rij<=rc && (e1==Nej || e2==Nej) ) {
                r2ij = rij * rij;
                pfcij = fc(rij, rcinv);
                pdfcij = dfc(rij, rcinv);
                for(int k=j+1; k<Nnum; k++) {
                    Nek = N->e[k];
                    if( (e1==Nej && e2==Nek) ||
                        (e2==Nej && e1==Nek) ) {
                        // energy calculation
                        rik = N->d[k];
                        if(rik <= rc) {
                            dxjk = N->dx[k] - N->dx[j];
                            dyjk = N->dy[k] - N->dy[j];
                            dzjk = N->dz[k] - N->dz[j];
                            rjk = dxjk * dxjk + dyjk * dyjk + dzjk * dzjk;
                            if(rjk <= rc2) {
                                rjk = sqrt(rjk);
                                rinvijik = 1.0 / rij / rik;
                                r2ik = rik * rik;
                                dxij = N->dx[j];
                                dyij = N->dy[j];
                                dzij = N->dz[j];
                                dxik = N->dx[k];
                                dyik = N->dy[k];
                                dzik = N->dz[k];
                                pfcik = fc(rik, rcinv);
                                pfcjk = fc(rjk, rcinv);
                                costijk = dxij * dxik + dyij * dyik + dzij * dzik; 
                                costijk *= rinvijik;
                                pfc = pfcij * pfcik * pfcjk; 
                                pexp = exp(-eta * (r2ij + r2ik + rjk * rjk));
                                plambda = 1.0 + lambda * costijk;
                                if(plambda <= 0.0) fg = 0.0;
                                else {
                                    if(optpow) fg = powint(plambda, optpow - 1) * pexp;
                                    else fg = pow(plambda, zeta - 1.0) * pexp;
                                }
                                value += fg * plambda * pfc;
                                // force calculation
                                fg *= pnorm;
                                rinvijik *= pzl;
                                costijk *= pzl;
                                p2etaplambda = 2.0 * eta * plambda;
                                p1 = pfc * (rinvijik - costijk / r2ij - p2etaplambda) + pfcik * pfcjk * pdfcij * plambda / rij;
                                p2 = pfc * (rinvijik - costijk / r2ik - p2etaplambda) + pfcij * pfcjk * dfc(rik, rcinv) * plambda / rik;
                                p3 = pfc * (rinvijik                  + p2etaplambda) - pfcij * pfcik * dfc(rjk, rcinv) * plambda / rjk;
                                p1 *= fg;
                                p2 *= fg;
                                p3 *= fg;
                                dxij *= p1;
                                dyij *= p1;
                                dzij *= p1;
                                dxik *= p2;
                                dyik *= p2;
                                dzik *= p2;
                                dxjk *= p3;
                                dyjk *= p3;
                                dzjk *= p3;
                                // save force contributions in force class
                                xptr[Nnum] += dxij + dxik;
                                yptr[Nnum] += dyij + dyik;
                                zptr[Nnum] += dzij + dzik;
                                xptr[j] -= dxij + dxjk;
                                yptr[j] -= dyij + dyjk;
                                zptr[j] -= dzij + dzjk;
                                xptr[k] -= dxik - dxjk;
                                yptr[k] -= dyik - dyjk;
                                zptr[k] -= dzik - dzjk;
                            } // rjk <= rc
                        } // rik <= rc
                    } // elem
                } // k
            } // rij <= rc
        } // j
        value *= pnorm;
    } // type 3

    // scale and/or center symmetry functions
    if     (scalesym == 1 && centersym == 0) svalue = (value - min) * sf + scalemin;
    else if(scalesym == 0 && centersym == 1) svalue = value - ave;
    else if(scalesym == 1 && centersym == 1) svalue = (value - ave) / (max - min);
    else svalue = value;

    // if SF value is out of boundaries give a warning
    if(value < min || value > max) {
        if(showew == 1) {
            if(screen) {
                fprintf(screen, "PROC %d: ### EXTRAPOLATION WARNING ### TS: %ld ATOM: %d SYMFUNC: %d TYPE: %d VALUE: %f MIN: %f MAX: %f\n",
                    comm->me, update->ntimestep, ca, num, type, value, min, max);
            }
            if(logfile) {
                fprintf(logfile, "PROC %d: ### EXTRAPOLATION WARNING ### TS: %ld ATOM: %d SYMFUNC: %d TYPE: %d VALUE: %f MIN: %f MAX: %f\n",
                    comm->me, update->ntimestep, ca, num, type, value, min, max);
            }
        }
//#ifdef _OPENMP
//        #pragma omp atomic
        ewcount++;
        ewcountsum++;
//#endif
    }

    return svalue;

}

// calculate one symmetry function group
void PairRuNNer::calc_G_dG_group(int ca, RuNNer_forces *F, RuNNer_symfuncGroup *GG, RuNNer_neighbor *N, double*& l_Gsv) {

    int numSF       = GG->numSF ;
    int type        = GG->type  ;
    int e1          = GG->e1    ;
    int e2          = GG->e2    ;
    int*& index     = GG->index ;
    double rc       = GG->rc    ;
    double*& min    = GG->min   ;
    double*& max    = GG->max   ;
    double*& ave    = GG->ave   ;
    double*& sf     = GG->sf    ;

    int*& optpow    = GG->optpow;
    int*& etaind    = GG->etaind;
    double*& eta    = GG->eta   ;
    double*& zeta   = GG->zeta  ;
    double*& lambda = GG->lambda;    
    double*& zetalambda  = GG->zetalambda  ;
    double*& etaizl = GG->etaizl;    
    double*& rs     = GG->rs    ;

    int Nej, Nek;
    int Nnum;
    double rcinv = 1.0 / rc;
    double rc2 = rc * rc;
    double rij, rik, rjk, rinvijik;
    double r2ij, r2ik;
    double pr1, pr2, pr3;
    double dxij, dyij, dzij;
    double dxik, dyik, dzik;
    double dxjk, dyjk, dzjk;
    double costijk;
    double pexp, plambda, pfc, pdfc;
    double pfcij, pfcik, pfcjk;
    double pdfcij, pdfcik, pdfcjk;
    double pfczl; 
    double p1, p2, p3, pdd;
    double q1, q2, q3;    
    double p1dxij, p1dyij, p1dzij;
    double p2dxik, p2dyik, p2dzik;
    double p3dxjk, p3dyjk, p3dzjk;    
    double p2etaplambda;
    double fg;

    double* value  = NULL;
    double* svalue = NULL;
    double* pnorm = NULL;
    double r2sum = 0.0;
    double vexp = 0.0;

    double* xptr = NULL;
    double* yptr = NULL;
    double* zptr = NULL;

    value = new double[numSF];
    svalue = new double[numSF];
    for(int i=0; i<numSF; i++) {
        value[i] = 0.0;
        svalue[i] = 0.0;
    }

    // type 2 symmetry function
    if(type == 2) {
        Nnum = N->num;
        for(int j=0; j<Nnum; j++) {
            Nej = N->e[j];
            rij = N->d[j];
            if(rij<=rc && e1 == Nej) {  
                dxij = N->dx[j];
                dyij = N->dy[j];
                dzij = N->dz[j];
                pfc = fc(rij, rcinv);
                pdfc = dfc(rij, rcinv);
                for(int isym=0; isym<numSF; isym++) {
                    // energy calculation
                    pexp = exp(-eta[isym] * (rij - rs[isym]) * (rij - rs[isym]));
                    //pexp = fastexp(-eta[isym] * (rij - rs[isym]) * (rij - rs[isym]));
                    value[isym] += pexp * pfc;
                    // force calculation
                    p1 = (pdfc - 2.0 * eta[isym] * (rij - rs[isym]) * pfc) * pexp / rij;
                    // save force contributions in force class
                    xptr = F->dGdx[index[isym]];
                    yptr = F->dGdy[index[isym]];
                    zptr = F->dGdz[index[isym]];
                    xptr[Nnum] += p1 * dxij;
                    yptr[Nnum] += p1 * dyij;
                    zptr[Nnum] += p1 * dzij;
                    xptr[j] -= p1 * dxij;
                    yptr[j] -= p1 * dyij;
                    zptr[j] -= p1 * dzij;
                }
            }
        }
    }

    // type 3 symmetry function
    if(type == 3) {
        pnorm = new double[numSF];
        for(int isym=0; isym<numSF; isym++) {
            pnorm[isym] = pow(2.0, 1.0 - zeta[isym]);
        }
        Nnum = N->num;
        for(int j=0; j<Nnum-1; j++) {
            Nej = N->e[j];
            rij = N->d[j];
            if( rij<=rc && (e1==Nej || e2==Nej) ) {
                r2ij = rij * rij;
                pfcij = fc(rij, rcinv);
                pdfcij = dfc(rij, rcinv);
                dxij = N->dx[j];
                dyij = N->dy[j];
                dzij = N->dz[j];
                for(int k=j+1; k<Nnum; k++) {
                    Nek = N->e[k];
                    if( (e1==Nej && e2==Nek) ||
                        (e2==Nej && e1==Nek) ) {
                        // energy calculation
                        rik = N->d[k];
                        if(rik <= rc) {
                            dxjk = N->dx[k] - dxij;
                            dyjk = N->dy[k] - dyij;
                            dzjk = N->dz[k] - dzij;
                            rjk = dxjk * dxjk + dyjk * dyjk + dzjk * dzjk;
                            if(rjk <= rc2) {                                
                                r2ik = rik * rik;
                                r2sum = r2ij + r2ik + rjk;
                                rjk = sqrt(rjk);
                                rinvijik = 1.0 / (rij * rik);
                                dxik = N->dx[k];
                                dyik = N->dy[k];
                                dzik = N->dz[k];
                                pfcik = fc(rik, rcinv);
                                pfcjk = fc(rjk, rcinv);
                                costijk = dxij * dxik + dyij * dyik + dzij * dzik; 
                                costijk *= rinvijik;
                                
                                pfc = pfcij * pfcik * pfcjk; 
                                pdfcik = dfc(rik, rcinv);
                                pdfcjk = dfc(rjk, rcinv);
                                
                                pr1 = pfcik * pfcjk * pdfcij / rij;
                                pr2 = pfcij * pfcjk * pdfcik / rik;
                                pr3 = pfcij * pfcik * pdfcjk / rjk;
                                q1 = rinvijik - costijk / r2ij;
                                q2 = rinvijik - costijk / r2ik;
                                q3 = rinvijik;
                                for(int isym=0; isym<numSF; isym++) {
                                    plambda = 1.0 + lambda[isym] * costijk;
                                    if(etaind[isym] == isym) {
                                        vexp = exp(-eta[isym] * r2sum);
                                        //vexp = fastexp(-eta[isym] * r2sum);
                                    }
                                    if(plambda <= 0.0) fg = 0.0;
                                    else {
                                        if(optpow[isym]) fg = powint(plambda, optpow[isym] - 1) * vexp;
                                        else fg = pow(plambda, zeta[isym] - 1.0) * vexp;
                                    }
                                    
                                    //! TODO one could save a heap by screening the function for value, i.e. not including contributions that are smaller than a tiny fraction of the current value
                                    pfczl = fg * plambda * pfc; 
                                    //if (pfczl < 1e-6*value[isym]) continue;
                                    value[isym] += pfczl;
                                    
                                    // force calculation
                                    fg *= pnorm[isym];
                                    p2etaplambda = plambda * etaizl[isym];
                                    /*
                                    pfczl = pfc*zl; 
                                    p1 = fg * (pfczl * (rinvijik - costijk / r2ij - p2etaplambda) + pr1 * plambda);
                                    p2 = fg * (pfczl * (rinvijik - costijk / r2ik - p2etaplambda) + pr2 * plambda);
                                    p3 = fg * (pfczl * (rinvijik                  + p2etaplambda) - pr3 * plambda);
                                    */
                                    pfczl = pfc*fg*zetalambda[isym];
                                    plambda *= fg; 
                                    p1 = (pfczl * (q1 - p2etaplambda) + pr1 * plambda);
                                    p2 = (pfczl * (q2 - p2etaplambda) + pr2 * plambda);
                                    p3 = (pfczl * (q3 + p2etaplambda) - pr3 * plambda);
                                    
                                    // save force contributions in force class
                                    xptr = F->dGdx[index[isym]];
                                    yptr = F->dGdy[index[isym]];
                                    zptr = F->dGdz[index[isym]];
                                                                          
                                    p1dxij=p1*dxij; p1dyij=p1*dyij; p1dzij=p1*dzij;
                                    p2dxik=p2*dxik; p2dyik=p2*dyik; p2dzik=p2*dzik;
                                    xptr[Nnum] += p1dxij + p2dxik;
                                    yptr[Nnum] += p1dyij + p2dyik;
                                    zptr[Nnum] += p1dzij + p2dzik;                   
                                    p3dxjk=p3*dxjk; p3dyjk=p3*dyjk; p3dzjk=p3*dzjk;                                     
                                    xptr[k] -= p2dxik - p3dxjk;
                                    yptr[k] -= p2dyik - p3dyjk;
                                    zptr[k] -= p2dzik - p3dzjk;
                                    xptr[j] -= p1dxij + p3dxjk;
                                    yptr[j] -= p1dyij + p3dyjk;
                                    zptr[j] -= p1dzij + p3dzjk;
                                } // isym
                            } // rjk <= rc
                        } // rik <= rc
                    } // elem
                } // k
            } // rij <= rc
        } // j
        for(int isym=0; isym<numSF; isym++) {
            value[isym] *= pnorm[isym];
        }
        delete[] pnorm;
    } // type 3

    for(int isym=0; isym<numSF; isym++) {
        // scale and/or center symmetry functions
        if     (scalesym == 1 && centersym == 0) svalue[isym] = (value[isym] - min[isym]) * sf[isym] + scalemin;
        else if(scalesym == 0 && centersym == 1) svalue[isym] = value[isym] - ave[isym];
        else if(scalesym == 1 && centersym == 1) svalue[isym] = (value[isym] - ave[isym]) / (max[isym] - min[isym]);
        else svalue[isym] = value[isym];

        // if SF value is out of boundaries give a warning
        if(value[isym] < min[isym] || value[isym] > max[isym]) {
            if(showew == 1) {
                if(screen) {
                    fprintf(screen, "PROC %d: ### EXTRAPOLATION WARNING ### TS: %ld ATOM: %d SYMFUNC: %d TYPE: %d VALUE: %f MIN: %f MAX: %f\n",
                        comm->me, update->ntimestep, ca, index[isym], type, value[isym], min[isym], max[isym]);
                }
                if(logfile) {
                    fprintf(logfile, "PROC %d: ### EXTRAPOLATION WARNING ### TS: %ld ATOM: %d SYMFUNC: %d TYPE: %d VALUE: %f MIN: %f MAX: %f\n",
                        comm->me, update->ntimestep, ca, index[isym], type, value[isym], min[isym], max[isym]);
                }
            }
            ewcount++;
            ewcountsum++;
        }
        l_Gsv[index[isym]] = svalue[isym];
    }

    delete[] value;
    delete[] svalue;

    return;

}

// calculate all symmetry functions for all atoms
void PairRuNNer::calc_all_symfunc() {

    int ca = 0;
    int ica = 0;
    int isym = 0;

    // loop over all elements
    for(int i=0; i<ne; i++) {
        // loop over all atoms of element i
//#ifdef _OPENMP
//        #pragma omp parallel for private(ica, ca, isym)
//#endif
        for(ica=0; ica<nlape[i]; ica++) {
            ca = laipe[i][ica];
            // calculate all SF of atom ca
            for(isym=0; isym<nsym[i]; isym++) {
                G_svalue[i][ica][isym] = calc_G_dG(la[ca], F[ca], G[i][isym], N[ca]);
            }
        }
    }

    if(ewcount > maxew) {
        if(screen)  fflush(screen);
        if(logfile) fflush(logfile);
        fflush(stderr);
        usleep(1000000);
        fprintf(stderr, "ERROR: Too many extrapolation warnings.\n");
        exit(EXIT_FAILURE);
    }

    return;

}

// calculate all symmetry function groups for all atoms
void PairRuNNer::calc_all_symfunc_group() {

    int ca = 0;
    int ica = 0;
    int isymg = 0;

    // loop over all elements
    for(int i=0; i<ne; i++) {
        // loop over all atoms of element i
        for(ica=0; ica<nlape[i]; ica++) {
            ca = laipe[i][ica];
            // calculate all SF groups of atom ca
            for(isymg=0; isymg<nsymg[i]; isymg++) {
                calc_G_dG_group(la[ca], F[ca], GG[i][isymg], N[ca], G_svalue[i][ica]);
            }
        }
    }

    if(ewcount > maxew) {
        if(screen)  fflush(screen);
        if(logfile) fflush(logfile);
        fflush(stderr);
        usleep(1000000);
        fprintf(stderr, "ERROR: Too many extrapolation warnings.\n");
        exit(EXIT_FAILURE);
    }

    return;

}

// send neural network parameters to all processors
void PairRuNNer::mpidistribute_nn() {

    // ALL
    MPI_Bcast(&ne,          1, MPI_INT,    0, world);
    MPI_Bcast(&nl,          1, MPI_INT,    0, world);
    MPI_Bcast(&scalesym,    1, MPI_INT,    0, world);
    MPI_Bcast(&centersym,   1, MPI_INT,    0, world);
    MPI_Bcast(&normnodes,   1, MPI_INT,    0, world);
    MPI_Bcast(&cutoff_type, 1, MPI_INT,    0, world);
    MPI_Bcast(&scalemin,    1, MPI_DOUBLE, 0, world);
    MPI_Bcast(&scalemax,    1, MPI_DOUBLE, 0, world);

    // PROC >0
    if(comm->me>0) {
        e = new RuNNer_element*[ne];    
        nsym = new int[ne];
        for(int i=0; i<ne; i++) {
            e[i] = new RuNNer_element(*this);
            nsym[i] = 0;
        }
        numnode = new int*[ne];
        actfunc = new int*[ne];
        for(int i=0; i<ne; i++) {
            numnode[i] = new int[nl+2];
            actfunc[i] = new int[nl+2];
            for(int j=0; j<nl+2; j++) {
                numnode[i][j] = 0;
                actfunc[i][j] = 0;
            }
            numnode[i][nl+1] = 1;
            actfunc[i][0] = 0;
        }
    }

    // ALL
    for(int i=0; i<ne; i++) {
        MPI_Bcast(&e[i]->nucch        , 1, MPI_INT   , 0, world);
        MPI_Bcast(&e[i]->atomic_energy, 1, MPI_DOUBLE, 0, world);
        MPI_Bcast(&e[i]->name         , 3, MPI_CHAR  , 0, world);
    }
    MPI_Bcast(nsym, ne, MPI_INT, 0, world);
    for(int i=0; i<ne; i++) {
        MPI_Bcast(numnode[i], nl + 2, MPI_INT, 0, world);
        MPI_Bcast(actfunc[i], nl + 2, MPI_INT, 0, world);
    }

    return;

}

// send symmetry function parameters to all processors
void PairRuNNer::mpidistribute_symfunc() {

    // PROC >0
    if(comm->me>0) {
        G = new RuNNer_symfunc**[ne];   
        for(int i=0; i<ne; i++) {
            G[i] = new RuNNer_symfunc*[nsym[i]];    
            for(int j=0; j<nsym[i]; j++) G[i][j] = new RuNNer_symfunc(*this);
        }
    }

    // ALL
    for(int i=0; i<ne; i++) {
        for(int j=0; j<nsym[i]; j++) {
            MPI_Bcast(&G[i][j]->num,    1, MPI_INT,    0, world);
            MPI_Bcast(&G[i][j]->type,   1, MPI_INT,    0, world);
            MPI_Bcast(&G[i][j]->ec,     1, MPI_INT,    0, world);
            MPI_Bcast(&G[i][j]->e1,     1, MPI_INT,    0, world);
            MPI_Bcast(&G[i][j]->e2,     1, MPI_INT,    0, world);
            MPI_Bcast(&G[i][j]->optpow, 1, MPI_INT,    0, world);
            MPI_Bcast(&G[i][j]->value,  1, MPI_DOUBLE, 0, world);
            MPI_Bcast(&G[i][j]->svalue, 1, MPI_DOUBLE, 0, world);
            MPI_Bcast(&G[i][j]->rc,     1, MPI_DOUBLE, 0, world);
            MPI_Bcast(&G[i][j]->eta,    1, MPI_DOUBLE, 0, world);
            MPI_Bcast(&G[i][j]->zeta,   1, MPI_DOUBLE, 0, world);
            MPI_Bcast(&G[i][j]->lambda, 1, MPI_DOUBLE, 0, world);
            MPI_Bcast(&G[i][j]->rs,     1, MPI_DOUBLE, 0, world);
            MPI_Bcast(&G[i][j]->min,    1, MPI_DOUBLE, 0, world);
            MPI_Bcast(&G[i][j]->max,    1, MPI_DOUBLE, 0, world);
            MPI_Bcast(&G[i][j]->ave,    1, MPI_DOUBLE, 0, world);
            MPI_Bcast(&G[i][j]->sf,     1, MPI_DOUBLE, 0, world);
        }
    }

    return;

}

// send neural network weights to all processors
void PairRuNNer::mpidistribute_nnweights() {

    // ALL
    for(int i=0; i<ne; i++) {
        for(int j=0; j<nl+1; j++) {
            for(int k=0; k<L[i][j]->nprevnode; k++) {
                MPI_Bcast(L[i][j]->weight[k], L[i][j]->nnode, MPI_DOUBLE, 0, world);
            }
            MPI_Bcast(L[i][j]->bias, L[i][j]->nnode, MPI_DOUBLE, 0, world);
        }
    }

    return;

}

// sum up extrapolation warning counters of all processors
void PairRuNNer::mpicollect_ew() {

    long ewcountsumglobal = 0;

    MPI_Allreduce(&ewcountsum, &ewcountsumglobal, 1, MPI_LONG, MPI_SUM, world);

    ewcountsum = ewcountsumglobal;

    return;

}

////////////////////////////////////////////////////////////////////////////////
// PairRuNNer::RuNNer_atom member functions
////////////////////////////////////////////////////////////////////////////////

PairRuNNer::RuNNer_atom::RuNNer_atom(PairRuNNer& tmp): par(tmp) {
        
    e  = 0;
    x  = 0.0;
    y  = 0.0;
    z  = 0.0;
    en = 0.0;
    fx = 0.0;
    fy = 0.0;
    fz = 0.0;

}

PairRuNNer::RuNNer_atom::~RuNNer_atom() {}

////////////////////////////////////////////////////////////////////////////////
// PairRuNNer::RuNNer_element member functions
////////////////////////////////////////////////////////////////////////////////

PairRuNNer::RuNNer_element::RuNNer_element(PairRuNNer& x): par(x) {
    
    nucch = 0;
    strcpy(name, "??");
    atomic_energy = 0.0;

}

PairRuNNer::RuNNer_element::~RuNNer_element() {}

PairRuNNer::RuNNer_element& PairRuNNer::RuNNer_element::operator=(const RuNNer_element& rhs) {

    if(this != &rhs) {
        nucch = rhs.nucch;
        strcpy(name, rhs.name);
    }

    return *this;

}

// get nuclear charge from name string
int PairRuNNer::RuNNer_element::getnucch() {

    nucch = 0;

    if     (strcmp(name, "H" ) == 0) nucch = 1;
    else if(strcmp(name, "He") == 0) nucch = 2;
    else if(strcmp(name, "Li") == 0) nucch = 3;
    else if(strcmp(name, "Be") == 0) nucch = 4;
    else if(strcmp(name, "B" ) == 0) nucch = 5;
    else if(strcmp(name, "C" ) == 0) nucch = 6;
    else if(strcmp(name, "N" ) == 0) nucch = 7;
    else if(strcmp(name, "O" ) == 0) nucch = 8;
    else if(strcmp(name, "F" ) == 0) nucch = 9;
    else if(strcmp(name, "Ne") == 0) nucch = 10;
    else if(strcmp(name, "Na") == 0) nucch = 11;
    else if(strcmp(name, "Mg") == 0) nucch = 12;
    else if(strcmp(name, "Al") == 0) nucch = 13;
    else if(strcmp(name, "Si") == 0) nucch = 14;
    else if(strcmp(name, "P" ) == 0) nucch = 15;
    else if(strcmp(name, "S" ) == 0) nucch = 16;
    else if(strcmp(name, "Cl") == 0) nucch = 17;
    else if(strcmp(name, "Ar") == 0) nucch = 18;
    else if(strcmp(name, "K" ) == 0) nucch = 19;
    else if(strcmp(name, "Ca") == 0) nucch = 20;
    else if(strcmp(name, "Sc") == 0) nucch = 21;
    else if(strcmp(name, "Ti") == 0) nucch = 22;
    else if(strcmp(name, "V" ) == 0) nucch = 23;
    else if(strcmp(name, "Cr") == 0) nucch = 24;
    else if(strcmp(name, "Mn") == 0) nucch = 25;
    else if(strcmp(name, "Fe") == 0) nucch = 26;
    else if(strcmp(name, "Co") == 0) nucch = 27;
    else if(strcmp(name, "Ni") == 0) nucch = 28;
    else if(strcmp(name, "Cu") == 0) nucch = 29;
    else if(strcmp(name, "Zn") == 0) nucch = 30;
    else if(strcmp(name, "Ga") == 0) nucch = 31;
    else if(strcmp(name, "Ge") == 0) nucch = 32;
    else if(strcmp(name, "As") == 0) nucch = 33;
    else if(strcmp(name, "Se") == 0) nucch = 34;
    else if(strcmp(name, "Br") == 0) nucch = 35;
    else if(strcmp(name, "Kr") == 0) nucch = 36;
    else if(strcmp(name, "Rb") == 0) nucch = 37;
    else if(strcmp(name, "Sr") == 0) nucch = 38;
    else if(strcmp(name, "Y" ) == 0) nucch = 39;
    else if(strcmp(name, "Zr") == 0) nucch = 40;
    else if(strcmp(name, "Nb") == 0) nucch = 41;
    else if(strcmp(name, "Mo") == 0) nucch = 42;
    else if(strcmp(name, "Tc") == 0) nucch = 43;
    else if(strcmp(name, "Ru") == 0) nucch = 44;
    else if(strcmp(name, "Rh") == 0) nucch = 45;
    else if(strcmp(name, "Pd") == 0) nucch = 46;
    else if(strcmp(name, "Ag") == 0) nucch = 47;
    else if(strcmp(name, "Cd") == 0) nucch = 48;
    else if(strcmp(name, "In") == 0) nucch = 49;
    else if(strcmp(name, "Sn") == 0) nucch = 50;
    else if(strcmp(name, "Sb") == 0) nucch = 51;
    else if(strcmp(name, "Te") == 0) nucch = 52;
    else if(strcmp(name, "I" ) == 0) nucch = 53;
    else if(strcmp(name, "Xe") == 0) nucch = 54;
    else if(strcmp(name, "Cs") == 0) nucch = 55;
    else if(strcmp(name, "Ba") == 0) nucch = 56;
    else if(strcmp(name, "La") == 0) nucch = 57;
    else if(strcmp(name, "Ce") == 0) nucch = 58;
    else if(strcmp(name, "Pr") == 0) nucch = 59;
    else if(strcmp(name, "Nd") == 0) nucch = 60;
    else if(strcmp(name, "Pm") == 0) nucch = 61;
    else if(strcmp(name, "Sm") == 0) nucch = 62;
    else if(strcmp(name, "Eu") == 0) nucch = 63;
    else if(strcmp(name, "Gd") == 0) nucch = 64;
    else if(strcmp(name, "Tb") == 0) nucch = 65;
    else if(strcmp(name, "Dy") == 0) nucch = 66;
    else if(strcmp(name, "Ho") == 0) nucch = 67;
    else if(strcmp(name, "Er") == 0) nucch = 68;
    else if(strcmp(name, "Tm") == 0) nucch = 69;
    else if(strcmp(name, "Yb") == 0) nucch = 70;
    else if(strcmp(name, "Lu") == 0) nucch = 71;
    else if(strcmp(name, "Hf") == 0) nucch = 72;
    else if(strcmp(name, "Ta") == 0) nucch = 73;
    else if(strcmp(name, "W" ) == 0) nucch = 74;
    else if(strcmp(name, "Re") == 0) nucch = 75;
    else if(strcmp(name, "Os") == 0) nucch = 76;
    else if(strcmp(name, "Ir") == 0) nucch = 77;
    else if(strcmp(name, "Pt") == 0) nucch = 78;
    else if(strcmp(name, "Au") == 0) nucch = 79;
    else if(strcmp(name, "Hg") == 0) nucch = 80;
    else if(strcmp(name, "Tl") == 0) nucch = 81;
    else if(strcmp(name, "Pb") == 0) nucch = 82;
    else if(strcmp(name, "Bi") == 0) nucch = 83;
    else if(strcmp(name, "Po") == 0) nucch = 84;
    else if(strcmp(name, "At") == 0) nucch = 85;
    else if(strcmp(name, "Rn") == 0) nucch = 86;
    else if(strcmp(name, "Fr") == 0) nucch = 87;
    else if(strcmp(name, "Ra") == 0) nucch = 88;
    else if(strcmp(name, "Ac") == 0) nucch = 89;
    else if(strcmp(name, "Th") == 0) nucch = 90;
    else if(strcmp(name, "Pa") == 0) nucch = 91;
    else if(strcmp(name, "U" ) == 0) nucch = 92;
    else if(strcmp(name, "Np") == 0) nucch = 93;
    else if(strcmp(name, "Pu") == 0) nucch = 94;
    else if(strcmp(name, "Am") == 0) nucch = 95;
    else if(strcmp(name, "Cm") == 0) nucch = 96;
    else if(strcmp(name, "Bk") == 0) nucch = 97;
    else if(strcmp(name, "Cf") == 0) nucch = 98;
    else if(strcmp(name, "Es") == 0) nucch = 99;
    else if(strcmp(name, "Fm") == 0) nucch = 100;
    else if(strcmp(name, "Md") == 0) nucch = 101;
    else if(strcmp(name, "No") == 0) nucch = 102;
    else {
        fprintf(stderr, "ERROR: Unknown element %s.\n", name);
        exit(EXIT_FAILURE);
    }

    return nucch;

}

////////////////////////////////////////////////////////////////////////////////
// PairRuNNer::RuNNer_forces member functions
////////////////////////////////////////////////////////////////////////////////

PairRuNNer::RuNNer_forces::RuNNer_forces(PairRuNNer& x, int nlayer, int *nnlayer, int element, int ai, int nn): par(x) { 

    e = element;
    nlf = nlayer + 1;
    nsym = nnlayer[0];
    indatom = ai;
    num = nn;
    nnode = new int[nlf];
    dEdG = new double[nsym];
    dfdx = new double*[nlf];
    for(int i=0; i<nlf; i++) {
        nnode[i] = nnlayer[i+1];
        dfdx[i] = new double[nnlayer[i+1]];
        for(int j=0; j<nnlayer[i+1]; j++) {
            dfdx[i][j] = 0.0;
        }
    }
    dGdx = new double*[nsym];
    dGdy = new double*[nsym];
    dGdz = new double*[nsym];
    for(int i=0; i<nsym; i++) {
        dEdG[i] = 0.0;
        dGdx[i] = new double[num+1]; 
        dGdy[i] = new double[num+1]; 
        dGdz[i] = new double[num+1]; 
        for(int j=0; j<num+1; j++) {
            dGdx[i][j] = 0.0;
            dGdy[i][j] = 0.0;
            dGdz[i][j] = 0.0;
        }
    }

}

PairRuNNer::RuNNer_forces::~RuNNer_forces() {

    delete[] nnode;
    for(int i=0; i<nlf; i++) {
        delete[] dfdx[i];
    }
    delete[] dfdx;
    for(int i=0; i<nsym; i++) {
        delete[] dGdx[i];
        delete[] dGdy[i];
        delete[] dGdz[i];
    }
    delete[] dGdx;
    delete[] dGdy;
    delete[] dGdz;
    delete[] dEdG;

}

// get saved df/dx from layer class
void PairRuNNer::RuNNer_forces::get_dfdx(RuNNer_layer **L) {

    for(int i=0; i<nlf; i++) {
        for(int j=0; j<nnode[i]; j++) {
            dfdx[i][j] = L[i]->dfdx[j];
        }
    }

    return;

}

// calculate dE/dG from saved values in layer class
void PairRuNNer::RuNNer_forces::calc_dEdG(RuNNer_layer **L) {

    double **tmpinner, **tmpouter;

    tmpinner = new double*[nlf-1];
    tmpouter = new double*[nlf-1];
    for(int i=0; i<nlf-1; i++) {
        tmpinner[i] = new double[nnode[i]];
        tmpouter[i] = new double[nnode[i+1]];
    }

    for(int k=0; k<nsym; k++) {
        for(int i=0; i<nnode[0]; i++) {
            tmpinner[0][i] = dfdx[0][i] * L[0]->weight[k][i];
            if(par.normnodes) tmpinner[0][i] /= L[0]->nprevnode;
        }
        for(int l=0; l<nlf-1; l++) {
            for(int i=0; i<nnode[l+1]; i++) {
                tmpouter[l][i] = 0.0;
                for(int j=0; j<nnode[l]; j++) {
                    tmpouter[l][i] += L[l+1]->weight[j][i] * tmpinner[l][j];
                }
                tmpouter[l][i] *= dfdx[l+1][i];
                if(par.normnodes) tmpouter[l][i] /= L[l+1]->nprevnode;
                if(l<nlf-2) tmpinner[l+1][i] = tmpouter[l][i];
            }
        }
        dEdG[k] = tmpouter[nlf-2][0];
    }

    for(int i=0; i<nlf-1; i++) {
        delete[] tmpinner[i];
        delete[] tmpouter[i];
    }
    delete[] tmpinner;
    delete[] tmpouter;

    return;

}

////////////////////////////////////////////////////////////////////////////////
// PairRuNNer::RuNNer_layer member functions
////////////////////////////////////////////////////////////////////////////////

PairRuNNer::RuNNer_layer::RuNNer_layer(PairRuNNer& x, int a, int b, int c): par(x) {

    nprevnode = a;
    nnode = b;
    af = c;
    input = new double[a];
    node = new double[b];
    dfdx = new double[b];
    bias = new double[b];
    weight = new double*[a];
    for(int i=0; i<a; i++) {
        input[i] = 0.0;
        weight[i] = new double[b];
        for(int j=0; j<b; j++) weight[i][j] = 0.0;
    }
    for(int i=0; i<b; i++) {
        node[i] = 0.0;
        dfdx[i] = 0.0;
        bias[i] = 0.0;
    }

}

PairRuNNer::RuNNer_layer::~RuNNer_layer() {

    delete[] input;
    delete[] node;
    delete[] dfdx;
    delete[] bias;
    for(int i=0; i<nprevnode; i++) {
        delete[] weight[i];
    }
    delete[] weight;

}

// get NN bias and weights from weights.???.data
void PairRuNNer::RuNNer_layer::get_weights(FILE *fp) {

    char line[LLINE];
    int ierr;
    char *cerr;

    // read weights
    for(int w1=0; w1<nprevnode; w1++) {
        for(int w2=0; w2<nnode; w2++) {
            cerr = fgets(line, LLINE, fp);
            ierr = sscanf(line, "%lf", &weight[w1][w2]);
        }
    }
    // read bias
    for(int b=0; b<nnode; b++) {
        cerr = fgets(line, LLINE, fp);
        ierr = sscanf(line, "%lf", &bias[b]);
    }

    return;

}

// calc all nodes of layer
void PairRuNNer::RuNNer_layer::calc_all_nodes() {

    for(int i=0; i<nnode; i++) {
        calc_one_node(i);
    } 

    return;

}

// calc one node of layer
void PairRuNNer::RuNNer_layer::calc_one_node(int inode) {

    double nv = 0.0;
    
    for(int j=0; j<nprevnode; j++) {
        nv += weight[j][inode] * input[j];
    }
    nv += bias[inode];
    if(par.normnodes) nv /= nprevnode;

    // save node values and derivatives for later usage
    if(af==1) {
        dfdx[inode] = 1.0;
        node[inode] = nv;
    }
    else if(af==2) {
        dfdx[inode] = 1.0 - pow(tanh(nv), 2.0);
        node[inode] = tanh(nv);
    }
    else {
        fprintf(stderr, "ERROR: Unknown activation function.\n");
        exit(EXIT_FAILURE);
    }   
    
    return;

}

// get input node values from previous layer
void PairRuNNer::RuNNer_layer::get_input(RuNNer_layer *lprev) {

    for(int i=0; i<nprevnode; i++) {
        input[i] = lprev->node[i];
    }

    return;

}

////////////////////////////////////////////////////////////////////////////////
// PairRuNNer::RuNNer_neighbor member functions
////////////////////////////////////////////////////////////////////////////////

PairRuNNer::RuNNer_neighbor::RuNNer_neighbor(PairRuNNer& x): par(x) {

    iatom = 0;
    eatom = 0;
    num   = 0;
    for(int i=0; i<MAXNEIGH; i++) {
        ilal[i] = 0;
        e[i]    = 0;
        d[i]    = 0.0;
        dx[i]   = 0.0;
        dy[i]   = 0.0;
        dz[i]   = 0.0;
    }

}

PairRuNNer::RuNNer_neighbor::~RuNNer_neighbor() {}

// calculate neighbor list using LAMMPS neighbor list
int PairRuNNer::RuNNer_neighbor::calcl(RuNNer_atom **a, int ii) {

    int i = par.list->ilist[ii];
    int ila = par.la[i];
    int j = 0;
    double tmpd;
    double tmpdx, tmpdy, tmpdz;
    int numneigh = par.list->numneigh[i];
    int* firstneigh = par.list->firstneigh[i];

    eatom = a[ila]->e;

    for(int jj=0; jj<numneigh; jj++) {
        j = firstneigh[jj];
        tmpdx = a[ila]->x - CONST_A2B * par.atom->x[j][0];
        tmpdy = a[ila]->y - CONST_A2B * par.atom->x[j][1];
        tmpdz = a[ila]->z - CONST_A2B * par.atom->x[j][2];
        tmpd  = sqrt(tmpdx * tmpdx + tmpdy * tmpdy + tmpdz * tmpdz);
        // if distance <= maximal cutoff add atom to neighbor list
        if(tmpd <= par.rcmax) {
            if(num >= MAXNEIGH) {
                fprintf(stderr, "ERROR: neighbor list array too small, increase MAXNEIGH and recompile.\n");
                exit(EXIT_FAILURE);
            }
            ilal[num] = j;
            //ilal2[j]   = num;
            e[num]    = par.transform_atom_type(par.atom->type[j])-1; //HELLSTROM
            d[num]    = tmpd;
            dx[num]   = tmpdx;
            dy[num]   = tmpdy;
            dz[num]   = tmpdz;
            num++;
        }
    }

    return num;

}

// calculate neighbor list
int PairRuNNer::RuNNer_neighbor::calc(RuNNer_atom **a, int i) {

    int j;
    double tmpd;
    double tmpdx, tmpdy, tmpdz;

    eatom = a[i]->e;

    for(j=0; j<par.atom->nlocal+par.atom->nghost; j++) {
        if( !( i==par.atom->tag[j]-1 && j<par.nla ) ) { 
            tmpdx = a[i]->x - CONST_A2B * par.atom->x[j][0];
            tmpdy = a[i]->y - CONST_A2B * par.atom->x[j][1];
            tmpdz = a[i]->z - CONST_A2B * par.atom->x[j][2];
            tmpd  = sqrt(tmpdx * tmpdx + tmpdy * tmpdy + tmpdz * tmpdz);
            // if distance <= maximal cutoff add atom to neighbor list
            if(tmpd <= par.rcmax) {
                if(num >= MAXNEIGH) {
                    fprintf(stderr, "ERROR: neighbor list array too small, increase MAXNEIGH and recompile.\n");
                    exit(EXIT_FAILURE);
                }
                ilal[num] = j;
                //ilal2[j]   = num;
                e[num]    = par.transform_atom_type(par.atom->type[j])-1; //HELLSTROM
                d[num]    = tmpd;
                dx[num]   = tmpdx;
                dy[num]   = tmpdy;
                dz[num]   = tmpdz;
                num++;
            }
        }
    }

    return num;

}

// find LAMMPS atom list index in RuNNer neighbor list
int PairRuNNer::RuNNer_neighbor::find_nindex(int n) {

    for(int i=0; i<num; i++) {
        if(ilal[i]==n) return i;
    }

    return -1;

}

////////////////////////////////////////////////////////////////////////////////
// PairRuNNer::RuNNer_symfunc member functions
////////////////////////////////////////////////////////////////////////////////

PairRuNNer::RuNNer_symfunc::RuNNer_symfunc(PairRuNNer& x): par(x) {

    num    = 0;
    type   = 0;
    ec     = 0;
    e1     = 0;
    e2     = 0;
    optpow = 0;
    value  = 0.0;
    svalue = 0.0;
    rc     = 0.0;
    eta    = 0.0;
    zeta   = 0.0;
    lambda = 0.0;
    rs     = 0.0;
    min    = 0.0;
    max    = 0.0;
    ave    = 0.0;
    sf     = 0.0;

}

PairRuNNer::RuNNer_symfunc::~RuNNer_symfunc() {}

PairRuNNer::RuNNer_symfunc& PairRuNNer::RuNNer_symfunc::operator=(const RuNNer_symfunc& rhs) {

    if(this != &rhs) {
        num    = rhs.num;
        type   = rhs.type;
        ec     = rhs.ec;
        e1     = rhs.e1;
        e2     = rhs.e2;
        optpow = rhs.optpow;
        value  = rhs.value;
        svalue = rhs.svalue;
        rc     = rhs.rc;
        eta    = rhs.eta;
        zeta   = rhs.zeta;
        lambda = rhs.lambda;
        rs     = rhs.rs;      
        min    = rhs.min;
        max    = rhs.max;
        ave    = rhs.ave;
        sf     = rhs.sf;
    }

    return *this;

}

// set up one symmetry function 
void PairRuNNer::RuNNer_symfunc::setup(char *line) {

    int tmp;
    char tmpename[3];
    char tmpe1name[3];
    char tmpe2name[3];

    // analyze char string passed from analyze_input_nn()
    sscanf(line, "%s %d", tmpename, &type);
    ec = par.find_element(tmpename);
    if(type==2) {
        sscanf(line, "%*s %*s %s %lf %lf %lf", tmpe1name, &eta, &rs, &rc);   
        e1 = par.find_element(tmpe1name);
    }
    else if(type==3) {
        sscanf(line, "%*s %*s %s %s %lf %lf %lf %lf", tmpe1name, tmpe2name, &eta, &lambda, &zeta, &rc);    
        e1 = par.find_element(tmpe1name);
        e2 = par.find_element(tmpe2name);
        if(e1>e2) {
            tmp = e1;
            e1 = e2;
            e2 = tmp;
        }
        if(fabs(zeta - round(zeta)) < SMALL_DOUBLE) {
            optpow = round(zeta);
        }
        //optpow = round(zeta+1);
    }

    return;

}

// info output
void PairRuNNer::RuNNer_symfunc::info() {

    if(type==2) {
        if(par.screen) fprintf(par.screen, "%4d %2s %2d %2s            %7.3f %7.3f %7.3f\n", num + 1, par.e[ec]->name, type, par.e[e1]->name, eta, rs, rc);
        if(par.logfile) fprintf(par.logfile, "%4d %2s %2d %2s            %7.3f %7.3f %7.3f\n", num + 1, par.e[ec]->name, type, par.e[e1]->name, eta, rs, rc);
    }
    else if(type==3) {
        if(par.screen) fprintf(par.screen, "%4d %2s %2d %2s %2s %7.3f %7.3f %7.3f %7.3f\n", num + 1, par.e[ec]->name, type, par.e[e1]->name, par.e[e2]->name, eta, lambda, zeta, rc);
        if(par.logfile) fprintf(par.logfile, "%4d %2s %2d %2s %2s %7.3f %7.3f %7.3f %7.3f\n", num + 1, par.e[ec]->name, type, par.e[e1]->name, par.e[e2]->name, eta, lambda, zeta, rc);
    }

    return;

}

////////////////////////////////////////////////////////////////////////////////
// PairRuNNer::RuNNer_symfuncGroup member functions
////////////////////////////////////////////////////////////////////////////////

PairRuNNer::RuNNer_symfuncGroup::RuNNer_symfuncGroup(PairRuNNer& x, int l_num, int l_numSF): par(x) {

    num    = l_num;
    numSF  = l_numSF;
    type   = 0;
    ec     = 0;
    e1     = 0;
    e2     = 0;
    index  = NULL;
    optpow = NULL;
    etaind = NULL;
    rc     = 0.0;
    eta    = NULL;
    zeta   = NULL;
    lambda = NULL;
    zetalambda = NULL;
    etaizl = NULL;
    rs     = NULL;
    min    = NULL;
    max    = NULL;
    ave    = NULL;
    sf     = NULL;

}

PairRuNNer::RuNNer_symfuncGroup::~RuNNer_symfuncGroup() {

    if(index  != NULL) delete[] index ;
    if(optpow != NULL) delete[] optpow;
    if(etaind != NULL) delete[] etaind;
    if(eta    != NULL) delete[] eta   ;
    if(zeta   != NULL) delete[] zeta  ;
    if(lambda != NULL) delete[] lambda;
    if(zetalambda   != NULL) delete[] zetalambda  ;
    if(etaizl != NULL) delete[] etaizl  ;    
    if(rs     != NULL) delete[] rs    ;
    if(min    != NULL) delete[] min   ;
    if(max    != NULL) delete[] max   ;
    if(ave    != NULL) delete[] ave   ;
    if(sf     != NULL) delete[] sf    ;

}

void PairRuNNer::RuNNer_symfuncGroup::setup(int l_type, int l_ec, int l_e1, int l_e2, double l_rc) {

    // SF group parameters
    type = l_type;
    ec   = l_ec;
    e1   = l_e1;
    e2   = l_e2;
    rc   = l_rc;

    // allocate type-specific arrays
    if(type == 2) {
        eta    = new double[numSF]; 
        rs     = new double[numSF]; 
        for(int i=0; i<numSF; i++) {
            eta   [i] = 0.0;
            rs    [i] = 0.0;
        }
    }
    else if(type == 3) {
        optpow = new int   [numSF]; 
        etaind = new int   [numSF]; 
        eta    = new double[numSF]; 
        zeta   = new double[numSF]; 
        lambda = new double[numSF]; 
        zetalambda = new double[numSF];
        etaizl = new double[numSF];
        for(int i=0; i<numSF; i++) {
            optpow[i] = 0  ;
            etaind[i] = 0  ;
            eta   [i] = 0.0;
            zeta  [i] = 0.0;
            lambda[i] = 0.0;
            zetalambda[i] = 0.0;
            etaizl[i] = 0.0;
        }

    }

    // allocate general SF arrays
    index  = new int   [numSF]; 
    min    = new double[numSF]; 
    max    = new double[numSF]; 
    ave    = new double[numSF]; 
    sf     = new double[numSF]; 
    for(int i=0; i<numSF; i++) {
        index [i] = 0  ;
        min   [i] = 0.0;
        max   [i] = 0.0;
        ave   [i] = 0.0;
        sf    [i] = 0.0;
    }

    return;

}

void PairRuNNer::RuNNer_symfuncGroup::info() {

    if(type == 2) {
        if(par.screen) {
            fprintf(par.screen, "%4d %2s %2d %2s                  *       * %7.3f %4d\n", num + 1, par.e[ec]->name, type, par.e[e1]->name, rc, numSF);
            for(int i=0; i<numSF; i++) {
                fprintf(par.screen, "                         %7.3f %7.3f              %4d %4d\n", eta[i], rs[i], i + 1, index[i] + 1);
            }
        }
        if(par.logfile) {
            fprintf(par.logfile, "%4d %2s %2d %2s                  *       * %7.3f %4d\n", num + 1, par.e[ec]->name, type, par.e[e1]->name, rc, numSF);
            for(int i=0; i<numSF; i++) {
                fprintf(par.logfile, "                         %7.3f %7.3f              %4d %4d\n", eta[i], rs[i], i + 1, index[i] + 1);
            }
        }
    }
    else if(type == 3) {
        if(par.screen) {
            fprintf(par.screen, "%4d %2s %2d %2s %2s       *       *       * %7.3f %4d\n", num + 1, par.e[ec]->name, type, par.e[e1]->name, par.e[e2]->name, rc, numSF);
            for(int i=0; i<numSF; i++) {
                fprintf(par.screen, "                 %7.3f %7.3f %7.3f              %4d %4d %4d\n", eta[i], lambda[i], zeta[i], i + 1, index[i] + 1, etaind[i] + 1);
            }
        }
        if(par.logfile) {
            fprintf(par.logfile, "%4d %2s %2d %2s %2s       *       *       * %7.3f %4d\n", num + 1, par.e[ec]->name, type, par.e[e1]->name, par.e[e2]->name, rc, numSF);
            for(int i=0; i<numSF; i++) {
                fprintf(par.logfile, "                 %7.3f %7.3f %7.3f              %4d %4d %4d\n", eta[i], lambda[i], zeta[i], i + 1, index[i] + 1, etaind[i] + 1);
            }
        }
    }

    return;

}

void PairRuNNer::RuNNer_symfuncGroup::setsc(int i, int l_index, double l_min, double l_max, double l_ave, double l_sf) {

    index [i] = l_index ;
    min   [i] = l_min   ; 
    max   [i] = l_max   ; 
    ave   [i] = l_ave   ; 
    sf    [i] = l_sf    ; 

    return;

}

void PairRuNNer::RuNNer_symfuncGroup::settype2(int i, double l_eta, double l_rs) {

    eta   [i] = l_eta   ; 
    rs    [i] = l_rs    ; 

    return;

}

void PairRuNNer::RuNNer_symfuncGroup::settype3(int i, double l_eta, double l_zeta, double l_lambda, int l_optpow) {

    eta   [i] = l_eta   ; 
    zeta  [i] = l_zeta  ;     
    lambda[i] = l_lambda; 
    zetalambda[i] = l_zeta*l_lambda; // precomputes some products
    etaizl[i] = 2.0*l_eta/(l_zeta*l_lambda);
    optpow[i] = l_optpow; 

    return;

}

void PairRuNNer::RuNNer_symfuncGroup::postset() {

    if(type == 3) {

        int index = 0;
        double etatmp = 0.0;

        index = 0;
        etaind[0] = index;
        etatmp = eta[0];
        for(int i=1; i<numSF; i++) {
            if(etatmp == eta[i]) {
                etaind[i] = index;
            } 
            else {
                index = i;
                etaind[i] = index;
                etatmp = eta[i];
            }
        }
    }

    return;

}
