////////////////////////////////////////////////////////////////////////////////
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

#ifdef PAIR_CLASS

PairStyle(runner, PairRuNNer)

#else

#ifndef LMP_PAIR_RUNNER_H
#define LMP_PAIR_RUNNER_H

#include "pair.h"
//#include <unordered_map>

#define MAXNEIGH 500
#define LLINE 256                                               // maximum number of characters per line
#define MPITAG_STD 0                                            // macro for standard MPI message tag
#define CONST_A2B 1.8897261328                                  // conversion factor Angstrom to Bohr
#define CONST_EV2HA 0.0367493254                                // conversion factor eV to Hartree
#define SMALL_DOUBLE 1.0E-15                                    // conversion factor eV to Hartree
#define MAXTRANSFORMATIONS 8                                    // HELLSTROM: maximum number of atom type transformations (e.g. from atom type "4" to "1", where "4" is deuterium and "1" is H and you want the same chemistry for both)

namespace LAMMPS_NS {

class PairRuNNer : public Pair {

    public:
    PairRuNNer(class LAMMPS *);
    virtual ~PairRuNNer();
    virtual void compute(int, int);
    void allocate();
    void settings(int, char **);
    void coeff(int, char **);
    void init_style();
    virtual double init_one(int, int);
    void write_restart(FILE *);
    void read_restart(FILE *);
    void write_restart_settings(FILE *);
    void read_restart_settings(FILE *);
    virtual double single(int, int, int, int, double, double, double, double &);

    class RuNNer_atom;                                          // contains atom positions and forces
    class RuNNer_element;                                       // contains element string and nuclear charge
    class RuNNer_forces;                                        // contains data for forces saved during SF and NN calculation
    class RuNNer_layer;                                         // contains one NN layer with nodes, weights and bias values
    class RuNNer_neighbor;                                      // contains neighbor list of one atom
    class RuNNer_symfunc;                                       // contains SF parameters


    ////////////////////////////////////////////////////////////////////////////////
    // class RuNNer_atom
    ////////////////////////////////////////////////////////////////////////////////

    class RuNNer_atom {
    
        public:
        int e;                                                  // element integer
        double x;                                               // x coordinate
        double y;                                               // y coordinate
        double z;                                               // z coordinate
        double en;                                              // atom energy
        double fx;                                              // force x component
        double fy;                                              // force y component
        double fz;                                              // force z component
    
        RuNNer_atom(PairRuNNer&);
        ~RuNNer_atom();
        
        private:
        PairRuNNer& par;
    
    };


    ////////////////////////////////////////////////////////////////////////////////
    // class RuNNer_element
    ////////////////////////////////////////////////////////////////////////////////

    class RuNNer_element {

        public:
        int nucch;                                              // nuclear charge of atom
        double atomic_energy;                                   // free atom energy
        char name[3];                                           // name string of element

        RuNNer_element(PairRuNNer&);
        ~RuNNer_element();
        RuNNer_element& operator=(const RuNNer_element& rhs);
        int getnucch();                                         // get nuclear charge from name
    
        private:
        PairRuNNer& par;

    };


    ////////////////////////////////////////////////////////////////////////////////
    // class RuNNer_forces
    ////////////////////////////////////////////////////////////////////////////////

    class RuNNer_forces {
    
        public:
        int e;                                                  // element integer
        int nlf;                                                // number of hidden layers + 1
        int nsym;                                               // number of SF
        int indatom;                                            // atom index
        int num;                                                // number of neighbors
        int *nnode;                                             // number of nodes of each layer
        double *dEdG;                                           // dE/dG for each SF
        double **dfdx;                                          // df/dx (where x=argument of act. func.) for each layer and node
        double **dGdx;                                          // dG/dx for all SFs and each neighbor atom
        double **dGdy;                                          // dG/dy for all SFs and each neighbor atom
        double **dGdz;                                          // dG/dz for all SFs and each neighbor atom
    
        RuNNer_forces(PairRuNNer&, int, int *, int, int, int);
        ~RuNNer_forces();
        void get_dfdx(RuNNer_layer **);                         // get saved df/dx from layer class
        void calc_dEdG(RuNNer_layer **);                        // calculate dE/dG from saved node values in layer class
    
        private:
        PairRuNNer& par;
    
    };

    ////////////////////////////////////////////////////////////////////////////////
    // class RuNNer_layer
    ////////////////////////////////////////////////////////////////////////////////

    class RuNNer_layer {

        public:
        int nnode;                                              // number of nodes of this layer
        int nprevnode;                                          // number of nodes of the previous layer
        int af;                                                 // activation function number
        double *input;                                          // node values of previous layer = input nodes
        double *node;                                           // node values of this layer
        double *dfdx;                                           // derivatives of act. func. with respect to its argument
        double *bias;                                           // bias values
        double **weight;                                        // weights

        RuNNer_layer(PairRuNNer&, int, int, int);
        ~RuNNer_layer();
        void get_weights(FILE *);                               // get weights and bias from file
        void calc_all_nodes();                                  // calculate all node values
        void calc_one_node(int);                                // calculate one node value
        void get_input(RuNNer_layer *);                         // get input node values from previous layer
        void info();                                            // print layer class info
    
        private:
        PairRuNNer& par;

    };

    ////////////////////////////////////////////////////////////////////////////////
    // class RuNNer_neighbor
    ////////////////////////////////////////////////////////////////////////////////

    class RuNNer_neighbor {

        public:
        int iatom;                                              // atom index
        int eatom;                                              // atom element number
        int num;                                                // number of neighbors (with box duplicates)
        int ilal[MAXNEIGH];                                     // index in LAMMPS atom list
        //std::unordered_map<int,int> ilal2;                      // index in LAMMPS atom list
        int e[MAXNEIGH];                                        // neighbor atom elements 
        double d[MAXNEIGH];                                     // distance to neighbors
        double dx[MAXNEIGH];                                    // x component of distance vector
        double dy[MAXNEIGH];                                    // y component of distance vector
        double dz[MAXNEIGH];                                    // z component of distance vector

        RuNNer_neighbor(PairRuNNer&);
        ~RuNNer_neighbor();
        int calc(RuNNer_atom **, int);                          // calculate neighbors list
        int calcl(RuNNer_atom **, int);                         // calculate neighbors list using LAMMPS neighbor list
        void info();                                            // print neighbor class info
        int find_nindex(int);                                   // give array index for given neighbor atom number
    
        private:
        PairRuNNer& par;

    };


    ////////////////////////////////////////////////////////////////////////////////
    // class RuNNer_symfunc
    ////////////////////////////////////////////////////////////////////////////////

    class RuNNer_symfunc {
    
        public:
        int num;                                                // SF number
        int type;                                               // SF type
        int ec;                                                 // SF element of central atom
        int e1;                                                 // SF element of neigbor 1
        int e2;                                                 // SF element of neigbor 2
        int optpow;                                             // wether optimized power function is used (zeta is an integer)
        double value;                                           // SF value (obsolete, replaced by G_svalue)
        double svalue;                                          // SF scaled value (obsolete, replaced by G_svalue)
        double rc;                                              // SF cutoff radius
        double eta;                                             // SF eta
        double zeta;                                            // SF zeta
        double lambda;                                          // SF lambda
        double rs;                                              // SF shift radius
        double min;                                             // SF minimum value (scaling.data)
        double max;                                             // SF maximum value (scaling.data)
        double ave;                                             // SF average value (scaling.data)
        double sf;                                              // SF scaling factor
    
        RuNNer_symfunc(PairRuNNer&);
        ~RuNNer_symfunc();
        RuNNer_symfunc& operator=(const RuNNer_symfunc& rhs);
        void setup(char *);                                     // set SF parameters from given text line
        void info();                                            // SF screen output after sorting
        
        private:
        PairRuNNer& par;
    
    };


    ////////////////////////////////////////////////////////////////////////////////
    // class RuNNer_symfuncGroup
    ////////////////////////////////////////////////////////////////////////////////

    class RuNNer_symfuncGroup {
    
        public:
        int num;                                                // SF group number
        int numSF;                                              // number of SF in this group
        int type;                                               // SF type
        int ec;                                                 // SF element of central atom
        int e1;                                                 // SF element of neigbor 1
        int e2;                                                 // SF element of neigbor 2
        int* index;                                             // SF index
        int* optpow;                                            // wether optimized power function is used (zeta is an integer)
        int* etaind;                                            // reuse of exponential function value
        double rc;                                              // SF cutoff radius
        double* eta;                                            // SF eta
        double* zeta;                                           // SF zeta
        double* lambda;                                         // SF lambda
        double* zetalambda;                                     // SF zeta*lambda
        double* etaizl;                                         // SF 2 eta/(zeta lambda)
        double* rs;                                             // SF shift radius
        double* min;                                            // SF minimum value (scaling.data)
        double* max;                                            // SF maximum value (scaling.data)
        double* ave;                                            // SF average value (scaling.data)
        double* sf;                                             // SF scaling factor

        RuNNer_symfuncGroup(PairRuNNer&, int, int);
        ~RuNNer_symfuncGroup();
        void setup(int, int, int, int, double);                 // set SF group parameters (type, ec, e1, e2, rc)
        void info();                                            // SF group screen output
        void setsc(int, int, double, double, double, double);   // add data for SF scaling
        void settype2(int, double, double);                     // add data for type 2 SF
        void settype3(int, double, double, double, int);        // add data for type 3 SF
        void postset();                                         // additional grouping
        
        private:
        PairRuNNer& par;
    
    };


    ////////////////////////////////////////////////////////////////////////////////
    // MEMBER VARIABLES of class PairRuNNer
    ////////////////////////////////////////////////////////////////////////////////

    char dirRuNNer[LLINE];                                      // directory with RuNNer files (input.nn, scaling.data, weights.???.data)
    int setup_completed;                                        // wether init_style is already completed
    int na;                                                     // total number of atoms
    int nla;                                                    // number of local atoms
    int *nlape;                                                 // number of local atoms per element
    int *la;                                                    // list of indices of local atoms
    int **laipe;                                                // list of indices of la array per element
    int ne;                                                     // number of elements
    int nl;                                                     // number of hidden layers
    int scalesym;                                               // scale symmetry functions
    int centersym;                                              // center symmetry functions
    int normnodes;                                              // normalize nodes
    int cutoff_type;                                            // cutoff type
    int ewcount;                                                // counter for extrapolation warnings
    int showew;                                                 // print extrapolation warnings
    int showewsum;                                              // print extrapolation warning summary
    int resetew;                                                // reset extrapolation warning counter every timestep
    int maxew;                                                  // maximum number of extrapolation warnings allowed
    int *nsym;                                                  // number of symmetry functions of each element
    int *nsymg;                                                 // number of symmetry function groups of each element
    int **numnode;                                              // number of nodes for each element and layer
    int **actfunc;                                              // activation function for each element and layer
    long ewcountsum;                                            // counter for extrapolation warnings for summary
    double rcmax;                                               // global maximum cutoff radius
    double entot;                                               // total energy
    double lat[3][3];                                           // box vectors
    double scalemin;                                            // min value for symmetry function scaling
    double scalemax;                                            // max value for symmetry function scaling
    double ***G_svalue;                                         // array for SF values
    
    RuNNer_atom **a;                                            // atoms
    RuNNer_element **e;                                         // element information
    RuNNer_forces **F;                                          // data needed for force calculation for each atom
    RuNNer_layer ***L;                                          // neural network layers
    RuNNer_neighbor **N;                                        // neighbor list
    RuNNer_symfunc ***G;                                        // symmetry functions
    RuNNer_symfuncGroup ***GG;                                  // symmetry function groups

    ////////////////////////////////////////////////////////////////////////////////
    // MEMBER FUNCTIONS of class PairRuNNer
    ////////////////////////////////////////////////////////////////////////////////

    // initizalization and memory allocation
    void init_vars();
    void pt_allocset();
    void pt_delete();

    // elements.h
    void sort_elements();
    int find_element(char *);
    
    // forces.h
    void init_forces();
    void calc_forces();

    // layer.h
    void init_layers();
    void read_weights();
    void calc_nn();

    // neighbor.h
    void calc_neighbor_list();
    void calc_neighbor_list_lammps();
    
    // readfiles.h
    void analyze_input_nn();

    // symfunc.h
    void sort_symfuncs();
    void read_scaling();
    void create_symfunc_groups();
    double rc_max();
    double fc(double, double);
    double dfc(double, double);
    double powint(double, unsigned int);
    double calc_G_dG(int, RuNNer_forces *, RuNNer_symfunc *, RuNNer_neighbor *);
    void calc_G_dG_group(int, RuNNer_forces *, RuNNer_symfuncGroup *, RuNNer_neighbor *, double *&);
    void calc_all_symfunc();
    void calc_all_symfunc_group();

    // MPI communication
    void mpidistribute_nn();
    void mpidistribute_nnweights();
    void mpidistribute_symfunc();
    void mpicollect_ew();

    protected:
    double cut_global;
    char **elemname;

    int ntransformations;  ///// HELLSTROM
    int transformfrom[MAXTRANSFORMATIONS]; //// HELLSTROM
    int transformto[MAXTRANSFORMATIONS]; //// HELLSTROM
    int transform_atom_type(int) const; //HELLSTROM

};

}

inline double fastexp(double x) {
    // fast exponential by scaling and squaring
    x /= 1024.; // 2**10
    // truncated Taylor expansion
    x = 1.0 + x * (1.0 + 0.5 * x * (1.0 + x * 0.3333333333333333 * (1.0 + 0.25 * x)));
    x *= x; x *= x; x *= x; x *= x;  
    x *= x; x *= x; x *= x; x *= x; 
    x *= x; x *= x;  
    return x;
}
#endif
#endif
