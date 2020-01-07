#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <string.h>
#include <unistd.h>

// start 4000 steps with
// gcc md_long_tox.c ;./a.out 1775.5 .001 160000 40 > tmp
// ./a.out 1787 .001 160000 40 > tmp   # 4000 steps; cpu-time: 670.552619 sec # T 897.85 with tox -0.05
// python -c 'import numpy as np;a=np.loadtxt("tmp");print a.mean()'

//////////////////////////// start definitions //////////////////////////////
#define N 2			//N*N*N*4 atome, in letzter dimension ineinandergeschoben
#define mass_element (26.982)	    // is fixed! masse Aluminum kg

// 900K GGA
#define alat_lattice  (4.13)     // if r0_mor is equal to alat
#define alat_morse    (4.13)     // if alat = 4.07 but morse is
#define a_mor         (1.432673)	   // m^-1 a  Morese parameter
#define D_mor         (0.279300)	   // J    D  Morse parameter
#define ktr_tox       (-0.00)	   //N/m	transversale kraftkonstante

#define michael_poly_yes_no 1
// for polynomial: a + b*r**(-1) + c*r**(-2) + d*r**(-3)
// for rcut = 0.88
//params_abcd   [-21.769  54.463 -45.499  12.692]
//params_abcd__ [-19.557990008508305, 49.93759327856454, -42.453041782507114, 12.015227100764047]   pc: 0.99935
//#define aa            (-19.557990008508305)
//#define bb            (49.93759327856454)
//#define cc            (-42.453041782507114)
//#define dd            (12.015227100764047)

// for rcut = 0.84
//params_abcd   [-21.769  54.463 -45.499  12.692]                                                   pc: 0.99935
//params_abcd__ [-20.452444362904057, 51.786965651529364, -43.70836806717401, 12.296391796377739]/  pc: 0.99936
//#define aa            (-21.769)
//#define bb            (54.463)
//#define cc            (-45.499)
//#define dd            (12.692)
#define poly_aa            (-20.452444362904057)
#define poly_bb            (51.786965651529364)
#define poly_cc            (-43.70836806717401)
#define poly_dd            (12.296391796377739)
#define poly_ee            (0.0)
#define poly_rcut          (0.88)

//////////////////////////// end definitions //////////////////////////////



// Mache zelle groesser/kleiner als morse constante
//#define d0 (3037000500./N)              // N = 1 --> d0 = 3.037e+09
//#define d0 (3037000500./N*1.001)        // N = 1 --> d0 = 3.04004e+09
//#define d0 (3037000500./N)              // N = 2 --> d0 = 1.5185e+09
//#define d0 (3037000500./N*1.001)        // N = 2 --> d0 = 1.52002e+09
//a0 = 2.14748e+09                      // N = 2 --> stays fiexed
//a0 = 4.29497e+09                      // N = 1 --> stays fiexed
//a0 describes the positions of the atoms initially
//d0 just describes the equilibrium distance and is just used in the morse
//to calculate the forces (and velocities)
//r0_eq_morse is always 2.92742e-10    (2.92742 = 4.14/sqrt(2))

//#define d0 (3037000500./N*1.0001)       // d0=2.920351 = 4.13/np.sqrt(2) -> d0*1.0001
#define d0 (3037000500./N*alat_morse/alat_lattice)              // is fixed! (d0=2.920351 = 4.13/np.sqrt(2))
#define a0 (4294967296./N)              // == 2^32 ; is fixed! (a0=4.13); 2^32 = 4,294,967,296 bytes; 4,294,967,296 / (1,024 x 1,024) = 4,096 MB = 4GB
#define a0_alat_lattice (a0/alat_lattice)
#define da0_alat_lattice (1./a0_alat_lattice)
#define req ((alat_lattice/sqrt(2.))/da0_alat_lattice)
//#define poly_0 (-(poly_dd/(2.*poly_rcut*poly_rcut)) - poly_cc/poly_rcut + poly_aa*poly_rcut + poly_bb*log(poly_rcut))
#define a_par (a_mor*1e10)	            // m^-1 a  Morese parameter
#define D_par (D_mor*1.602176e-19)	    // J    D  Morse parameter
#define r0_mor (alat_morse/sqrt(2))	    // m    r0 Morse parameter for a = 4.14/sqrt(2)
#define r0_eq_morse (r0_mor*1e-10)	    // m    r0 Morse parameter for a = 4.14
#define one_over_r0_eq_mor (1/d0*r0_eq_morse)
#define ktr_par (ktr_tox*2*alat_lattice*1.602176e-9)	// N/m	transversale kraftkonstante --> das alat_lattice kuerzt sich dann weg  -> da das x an sich mit alat_lattice skaliert --> x = alat_lattice/a0 (das tox ist unabhaengig von der gitterkonstante)
#define a0ktr_par (ktr_par/a0)     // This saves 2% of total time;
#define m_element (mass_element*1.660538e-27)
#define k_B 1.380648e-23	            // is fixed! J/K
#define atoms (N*N*N*4)                 // number of atoms
#define first_neighbors (12)
#define NELEMS(x)  (sizeof(x) / sizeof((x)[0]))
#define stepsmax     (99940000)
#define stepsdudlmax (200000)
#define faktor_temperature (m_element/k_B/(12.*N*N*N))

double faktor_rel_pos=1./(a0*N); // show positions in relative coordinates
double faktor=alat_lattice/a0;   // show positions in absolute (angstrom)
//double faktor=1./a0;           // show positions in internal coords

double faktor_hesse=-1000./(atoms-1.)/2.;
//double faktor_hesse2=faktor*faktor*faktor_hesse;
double faktor_hesse2=(alat_lattice/a0)*(alat_lattice/a0)*(-1000./(atoms-1.)/2.);
double faktor_energy_cell=(1/1.602176e-19);
double faktor_energy_atom=(1/1.602176e-19)*1000./(atoms-1.);
double faktor_per_atom=1000./(atoms-1.);
double faktor_force=1./1.602176e-9;
double faktor_harm_vel=alat_lattice*1.602176e-9/a0;
//printf("hallokk %d\n",faktor);
double rmax=0;
double rmin=a0;
double projx=0;
double projy=0;
double projz=0;

// File in
FILE *file_in_positions;
FILE *file_in_hesse;

FILE *file_out_positions;
FILE *file_out_temp;
FILE *file_out_temp_av;
FILE *file_out_dudl;
FILE *file_out_dudl_av;
FILE *file_forces;
FILE *file_forces_av;
FILE *file_forces_vs_forces_dft;
FILE *file_out_check_dist_xyz;
FILE *file_out_check_dist_r;
FILE *file_out_check_dist_nn;
FILE *file_out_check_dist_nn_proj;
FILE *file_out_prl15_2a;
FILE *file_out_prl15_2au;

FILE *file_out_new1;
FILE *file_out_new2;

unsigned int MEIN_RAND_STATE=4;

unsigned int mein_rand(void){
	MEIN_RAND_STATE=1664525L*MEIN_RAND_STATE+1013904223L;
	return MEIN_RAND_STATE;
}

struct Pos  {unsigned int x,y,z;} * pos;
struct Pos0 {unsigned int x,y,z;} * pos0;
struct Du   {signed int x,y,z;} * du; // difference between pos and pos0
struct SqrtDuDuMean   {double x,y,z;} * sqrtdudumean; // sqrt((pos-pos0)^2).mean() for every atom
struct Vel  {double x,y,z;} * vel;
struct Velff  {double x,y,z;} * velff;
struct L1nn    {unsigned int ind1,ind2,i1,i2,i3,j1,j2,j3,k1,k2,k3,x0,y0,z0,x1,y1,z1;} * l1nn;
struct NN1     {unsigned int at1,at2,j1,j2,j3;} * nn1;

struct Proj {double x,y,z;} proj;
struct Forcetmp {double x,y,z;} * forcetmp;
struct Forcepoly {double x,y,z;} * forcepoly;
struct Force {double x,y,z;} * force;
struct Forcedft {double x,y,z;} * forcedft;
struct Forceharm {double x,y,z;} * forceharm;

struct Dudl {double u_dft,u_la,u_harm,\
    u_dft_sum,u_la_sum,u_harm_sum,\
    u_dft_av,u_la_av,u_harm_av,\
    dudl_dft_la,dudl_dft_harm,\
    dudl_dft_la_sum,dudl_dft_harm_sum,\
    dudl_dft_la_av,dudl_dft_harm_av,\
    dudl_dft_la_std,dudl_dft_harm_std
    ;} * dudl;



double u_la=0;
double u_pl=0;
double u_la_av=0;

double u_la_tox=0;
double u_la_tox_av=0;

double u_harm=0;
double u_harm_av=0;

double u_dft=0;
double u_dft_av=0;

double temperature=0;
double temperature_av=0;
double temperature_index=0;

double dudl_dft_la=0;
double dudl_dft_la_av=0;

double dudl_dft_harm=0;
double dudl_dft_harm_av=0;

double fdiff_dft_la_av=0;
double fdiff_dft_harm_av=0;
double fdudl_dft_harm_av=0;

double fstd_dft=0;
double fstd_la=0;
double fcov_dft_la=0;

double forces_diffmax=0;

double getaverage(double a, int T, double val)
{ //definition
    //return a+b;
    //printf("---> a %4.5f T %d val %4.5f\n",a,T,val);
    //printf("---> a*T+val %4.5f\n",a*T+val);
    //printf("---> (a*T+val) %4.5f\n",(a*T+val)/(T+1.));
    //double out;
    //out = ((a*T+val)/(T+1.));
    //printf("---> out %4.5f\n",out);
    return((a*T+val)/(T+1.));
    //return out;
}

////////////// print to screen
void init_l1nn() {
	int i,i1,i2,i3,j1,j2,j3,steps,ind1,ind2,a,b,c;
	i =0;
	for (i3=0;i3<N*2;i3++) for (i2=0;i2<N*2;i2++) for (i1=(i3+i2)%2;i1<N*2;i1+=2) {
	    // i{0,1,2,3} sind die absoluten koordinaten dat atome von 0 bis 2*N-1 (0-3 in 2x2x2sc)
	    // i{0,1,2,3} sind einfach die koordinaten
	    // {0,0,0} --> xyz = {0,0,0}
	    // {1,1,0} --> xyz = {2.065,2.065,0}
	    // {2,0,0} --> xyz = {4.13,0,0}
	    // {3,1,0} --> xyz = {6.195,2.065,0}
	    // --> i3 ist die xy ebene --> die atomstruktur wird in z richtung xy-ebene fuer-xy-ebene aufgebaut
	    //                            also erst alle z=0, dann alle z=2.065 usw.
	    //                            i3 geht von 0 bis 2*N-1 (0-2 in 2x2x2sc)
	    // --> i2 ist die xz koordinate --> i2 geht von 0 bis 2*N-1 (0-2 in 2x2x2sc)
	    // --> i1 ist die yz koordinate und geht im prinzip auch ueber alle indizes
	    // von 0 bis 2*N-1 (0-3 in 2x2x2sc) wobei jedoch nie ein atom bei {1,0,0}
	    // sitzt: {0,0,0} ist das aufatom und {2,0,0} ist der zweit naechste
	    // nachbar und {1,1,0} ist der naechste nachbar; von daher das modulo
	    // (i3+i2)%2 das dafuer sorgt dass man die richtige fcc struktur hat.
		ind1=((i3/2)*2*N+i2)*2*N+i1;
		j1=-1;  // das sind die relativen coordinaten der naechsten nachbarn
		j2=1;   // das sind die relativen coordinaten der naechsten nachbarn
		j3=0;   // das sind die relativen coordinaten der naechsten nachbarn

        /// Analysis
        // not all of the atoms pos[ind1].{x,y,z} exist
        //printf("ind1: %d pos[ind1].{x,y,z} %3.3f %3.3f %3.3f\n",ind1,pos[ind1].x*faktor,pos[ind1].z*faktor,pos[ind1].z*faktor);
		goto JUMPIN;
		for (j3=-1;j3<=1;j3++) for (j2=-1;j2<=1;j2++) for (j1=-1+(j3+j2+3)%2;j1<=1;j1+=2) {
		    // j{1,2,3} sind die relativen coordinaten der naechsten nachbarn
		    // z.b. {-1,1,0},{1,1,0},{0,-1,1},{-1,0,1},{1,0,1},{0,1,1}
		    // j{2,3} gehen ueber {-1,0,1} (also die moeglichen positionen der
		    // Naechsten nachbarn
		    // j1 im prinzip auch
JUMPIN:
			ind2=(((i3+j3+2*N)%(2*N)/2)*2*N+(i2+j2+2*N)%(2*N))*2*N+(i1+j1+2*N)%(2*N);
            //printf("ind1: %d pos[ind1].{x,y,z} %3.3f %3.3f %3.3f ind2: %d pos[ind2] %3.3f %3.3f %3.3f\n",ind1,pos[ind1].x*faktor,pos[ind1].y*faktor,pos[ind1].z*faktor,ind2,pos[ind2].x*faktor,pos[ind2].y*faktor,pos[ind2].z*faktor);

            l1nn[i].ind1=ind1;
            l1nn[i].ind2=ind2;
            l1nn[i].i1=i1;
            l1nn[i].i2=i2;
            l1nn[i].i3=i3;

            l1nn[i].j1=j1;
            l1nn[i].j2=j2;
            l1nn[i].j3=j3;

            a=l1nn[i].i1+l1nn[i].j1;
            b=l1nn[i].i2+l1nn[i].j2;
            c=l1nn[i].i3+l1nn[i].j3;
            a=a%(2*N);
            b=b%(2*N);
            c=c%(2*N);
            if (a==-1) a=2*N-1;
            if (b==-1) b=2*N-1;
            if (c==-1) c=2*N-1;
            l1nn[i].k1=a;
            l1nn[i].k2=b;
            l1nn[i].k3=c;

            l1nn[i].x0=pos[ind1].x;
            l1nn[i].y0=pos[ind1].y;
            l1nn[i].z0=pos[ind1].z;
            l1nn[i].x1=pos[ind2].x;
            l1nn[i].y1=pos[ind2].y;
            l1nn[i].z1=pos[ind2].z;


            //printf("i: %d\n",i);
            i++;
            }
        }



    ///////////////////////////////////////////////////////////////////////////
    // just for checking
    ///////////////////////////////////////////////////////////////////////////
    //i=0;
    //for (i=0;i<atoms*first_neighbors/2;i++) {
	//	printf("%d %d %d   %-2d %-2d %-2d   ind1: %-3d ind2: %-3d at1:%3.3f %3.3f %3.3f    at2:%3.3f %3.3f %3.3f --> i:%-3d\n",l1nn[i].i1,l1nn[i].i2,l1nn[i].i3,l1nn[i].j1,l1nn[i].j2,l1nn[i].j3,l1nn[i].ind1,l1nn[i].ind2,pos[l1nn[i].ind1].x*faktor,pos[l1nn[i].ind1].y*faktor,pos[l1nn[i].ind1].z*faktor,pos[l1nn[i].ind2].x*faktor,pos[l1nn[i].ind2].y*faktor,pos[l1nn[i].ind2].z*faktor,i);
    //}
    //printf("arr   : %lu\n",sizeof(l1nn));
    //printf("arr[0]: %lu\n",sizeof(l1nn[0]));
    //printf("arr[1]: %lu\n",sizeof(l1nn[1]));
    //exit(1);
}

void init(double T) {
	int i1,i2,i3,ind;
	double vx=0.,vy=0.,vz=0.;

	pos=(struct Pos *)malloc(sizeof(struct Pos)*atoms);
	pos0=(struct Pos0 *)malloc(sizeof(struct Pos0)*atoms);
	du=(struct Du *)malloc(sizeof(struct Du)*atoms);
	sqrtdudumean=(struct SqrtDuDuMean *)malloc(sizeof(struct SqrtDuDuMean)*atoms);
	vel=(struct Vel *)malloc(sizeof(struct Vel)*atoms);
	velff=(struct Velff *)malloc(sizeof(struct Velff)*atoms);
	l1nn=(struct L1nn *)malloc(sizeof(struct L1nn)*atoms*first_neighbors/2); // 12 neighbors, to exclude double counting /2 (every atoms interacts just ones with its nearest neighbor)
	nn1=(struct NN1 *)malloc(sizeof(struct NN1)*atoms*40);

    forcetmp=(struct Forcetmp *)malloc(sizeof(struct Forcetmp)*atoms);
    forcepoly=(struct Forcepoly *)malloc(sizeof(struct Forcepoly)*atoms);
	force=(struct Force *)malloc(sizeof(struct Force)*atoms);
	forcedft=(struct Forcedft *)malloc(sizeof(struct Forcedft)*atoms);
	forceharm=(struct Forceharm *)malloc(sizeof(struct Forceharm)*atoms);
	dudl=(struct Dudl *)malloc(sizeof(struct Dudl)*stepsdudlmax);  // this keeps
	//increasing the memory

	for (i3=0;i3<N*2;i3++) for (i2=0;i2<N*2;i2++) for (i1=(i3+i2)%2;i1<N*2;i1+=2) {
	    //printf("kkk %d %d %d\n",i1,i2,i3);
		ind=((i3/2)*2*N+i2)*2*N+i1;
		vx+=(vel[ind].x=(((signed)mein_rand())/2147483648.*sqrt(3.*k_B*T/m_element)));
		vy+=(vel[ind].y=(((signed)mein_rand())/2147483648.*sqrt(3.*k_B*T/m_element)));
		vz+=(vel[ind].z=(((signed)mein_rand())/2147483648.*sqrt(3.*k_B*T/m_element)));
		pos[ind].x=(unsigned int)(i1*a0/2.);
		pos[ind].y=(unsigned int)(i2*a0/2.);
		pos[ind].z=(unsigned int)(i3*a0/2.);
		pos0[ind].x=(unsigned int)(i1*a0/2.);
		pos0[ind].y=(unsigned int)(i2*a0/2.);
		pos0[ind].z=(unsigned int)(i3*a0/2.);
	    //printf("%d\t%.9f\t%.9f\t%.9f\n",ind,pos0[ind].x/a0/N,pos0[ind].y/a0/N,pos0[ind].z/a0/N);
	}
	vx/=atoms;
	vy/=atoms;
	vz/=atoms;

	for (ind=0;ind<atoms;ind++) {
		vel[ind].x-=vx;
		vel[ind].y-=vy;
		vel[ind].z-=vz;
		velff[ind].x=vel[ind].x;
		velff[ind].y=vel[ind].y;
		velff[ind].z=vel[ind].z;
	}

   init_l1nn();
}

void printsep(char stringin[]) {
    int i,maxlength,addchar;
    maxlength=130;
    maxlength=100;
    addchar = maxlength-strlen(stringin);
	//printf("%d\n",strlen(stringin));
	//printf("%s ",stringin);
	printf("%s",stringin);
	for (i=0;i<addchar;i++) {
        printf("=");
    //printf("%.*s", 25, "=");
    }
    //printf("%.*s", 25, "=================");
    printf("\n");
}

void print_l1nn_to_screen(int i,int j,char stringin[]) {
    int ind1,ind2,c1,c2,c3,c4,c5,c6;
    double x,y,z,r,faktrr;
	faktrr=faktor/alat_lattice*2;
	printsep(stringin);
    for (i=i;i<j;i++) {
        ind1=l1nn[i].ind1;
        ind2=l1nn[i].ind2;
		x=(signed)(pos0[ind1].x-pos0[ind2].x);
		y=(signed)(pos0[ind1].y-pos0[ind2].y);
		z=(signed)(pos0[ind1].z-pos0[ind2].z);
		r=sqrt(x*x+y*y+z*z);
        c1=pos0[ind1].x*faktrr-l1nn[i].i1;
        c2=pos0[ind1].y*faktrr-l1nn[i].i2;
        c3=pos0[ind1].z*faktrr-l1nn[i].i3;
        c4=pos0[ind2].x*faktrr-l1nn[i].k1;
        c5=pos0[ind2].y*faktrr-l1nn[i].k2;
        c6=pos0[ind2].z*faktrr-l1nn[i].k3;


	    printf("i:%-3d  %d %d %d   %-2d %-2d %-2d  %d %d %d  ind1/2: %-3d %-3d at1(init):%3.2f %3.2f %3.2f at2(init):%-5.2f %-5.2f %-5.2f at1:%-5.2f %-5.2f %-5.2f   at2:%-5.2f %-5.2f %-5.2f (%d %d %d, %d %d %d) r:%3.2f xyz: %3.2f %3.2f %3.2f\n",i,\
	    l1nn[i].i1,l1nn[i].i2,l1nn[i].i3,\
	    l1nn[i].j1,l1nn[i].j2,l1nn[i].j3,\
	    l1nn[i].k1,l1nn[i].k2,l1nn[i].k3,\
	    l1nn[i].ind1,l1nn[i].ind2,\
	    l1nn[i].x0*faktor,l1nn[i].y0*faktor,l1nn[i].z0*faktor,\
	    l1nn[i].x1*faktor,l1nn[i].y1*faktor,l1nn[i].z1*faktor,\
	    pos0[ind1].x*faktor,pos0[ind1].y*faktor,pos0[ind1].z*faktor,\
	    pos0[ind2].x*faktor,pos0[ind2].y*faktor,pos0[ind2].z*faktor,\
	    c1,c2,c3,c4,c5,c6,\
	    r*faktor,x*faktor,y*faktor,z*faktor);
	}
	printsep(stringin);
}

void projection(double x, double y, double z, double x0,double y0,double z0)
{
    double v0v,nv0,fakt;
    nv0 = sqrt(x0*x0+y0*y0+z0*z0); // Norm[v0]
    v0v = x*x0+y*y0+z*z0; // v.v0
    // v0v/nv0 ---> scalar
	fakt=v0v/nv0/nv0;
	//printf("in projection %.10f %.10f %.10f %.10f %.10f %.10f\n",x,y,z,x0,y0,z0);
	//printf("in projection %.10f %.10f %.10f %.10f %.10f %.10f\n",x,y,z,x0*fakt,y0*fakt,z0*fakt);
	//return(x0*fakt,y0*fakt,z0*fakt);
	projx = x0*fakt;
	projy = y0*fakt;
	projz = z0*fakt;
}


void print_pos_0_to_screen(int idxmax,char stringin[]) {
    int j;
    if (idxmax > atoms)
        idxmax = atoms;
	printsep(stringin);
	for (j=0;j<idxmax;j++) printf("%-5d %16.6f %16.6f %16.6f\n",j,pos0[j].x*faktor,pos0[j].y*faktor,pos0[j].z*faktor);
	printsep(stringin);
}

void get_1NN(int idxmax,char stringin[]) {
    int j,at1,at2,add,ii;
    int i=0;
    double dx,dy,dz,r,dxdx,dydy,dzdz;
    if (idxmax > atoms)
        idxmax = atoms;
	printsep(stringin);
	//for (j=0;j<atoms;j++) printf("%-5d %16.6f %16.6f %16.6f\n",j,pos0[j].x*faktor/alat_lattice,pos0[j].y*faktor/alat_lattice,pos0[j].z*faktor/alat_lattice);

	for (at1=0;at1<atoms;at1++) {
	    //printf("%-5d %16.6f %16.6f %16.6f\n",at1,pos0[at1].x*faktor/alat_lattice,pos0[at1].y*faktor/alat_lattice,pos0[at1].z*faktor/alat_lattice);
	    for (at2=0;at2<atoms;at2++) {
	        dx = (pos0[at2].x-pos0[at1].x)*faktor/alat_lattice;
	        dy = (pos0[at2].y-pos0[at1].y)*faktor/alat_lattice;
	        dz = (pos0[at2].z-pos0[at1].z)*faktor/alat_lattice;
            if (dx > N/2.) dx = dx - N;
            if (dy > N/2.) dy = dy - N;
            if (dz > N/2.) dz = dz - N;
		    r=sqrt(dx*dx+dy*dy+dz*dz);
		    //if r < 0.9 an r > 0.:
		    //    list = [at1,at2]

	        //dxm = dx%N;
	        //if ((at2 == 24) || (at2 == 26) || (at2 == 28)) {
	        if ((r < 0.9) && (r > 0.1) && (at1 < 2)) {
	        //printf("diff at1 %-5d at2 %-5d %16.6f %16.6f %16.6f r %16.6f i %-5d\n",at1,at2,dx,dy,dz,r,i);
            }

	        if ((r < 0.9) && (r > 0.1)) {
	            add = 1;
	            // check if it needs to be added
	            for (ii=0;ii<atoms*20;ii++) { // max 40 NN
                    if ((nn1[ii].at2 == at1) && (nn1[ii].at1 == at2)) {add = 0;};
	            }
	        if (add == 1) {
            nn1[i].at1=at1;
            nn1[i].at2=at2;

	        //printf("!add 1  at1 %-5d at2 %-5d %16.6f %16.6f %16.6f r %16.6f i %-5d\n",at1,at2,dx,dy,dz,r,i);
            // for tox stuff fcc
            //nn1[i].j1 = dx;
            //nn1[i].j2 = dy;
            //nn1[i].j3 = dz;
            if (dx*dx < 0.01) {nn1[i].j1=0;nn1[i].j2=1;nn1[i].j3=1;};
            if (dy*dy < 0.01) {nn1[i].j1=1;nn1[i].j2=0;nn1[i].j3=1;};
            if (dz*dz < 0.01) {nn1[i].j1=1;nn1[i].j2=1;nn1[i].j3=0;};
	        //printf("!add 2  at1 %-5d at2 %-5d %-5d %-5d %-5d r %16.6f i %-5d\n",at1,at2,nn1[i].j1,nn1[i].j2,nn1[i].j3,r,i);
            i++;
            }
	        };
	        };
	    //if (at1 == 1) {exit(1);};
    }
    //for (i=0;i<atoms*6+10;i++) {
	//        printf("final i %-3d at1 %-3d at2 %-3d dx %-3d dy %-3d dz %-3d\n",i,nn1[i].at1,nn1[i].at2,nn1[i].j1,nn1[i].j2,nn1[i].j3);}
    //exit(1);

	printsep(stringin);
}

void print_pos_to_screen(int idxmax,char stringin[],int step) {
    int j;
    //double faktor=1./(a0*N);            // show in relative
    //double faktor=alat_lattice/a0;   // show in absolute
    //double faktor=1./a0;            // show in internal coords

    //char str3[11] = "Call home!";
    //char str1[] = strcat(str3,"joooooooooooooo\n");
    if (idxmax > atoms)
        idxmax = atoms;
    //printf("%s",str3);
    //printf(str1);
    //printf(strcat(stringin, str3));
	//printf(strcat(stringin,"-------------------------------- pos 0 start ---------------------------\n"));
	//printf(stringin,"-------------------------------- pos 0 start ---------------------------\n");
	//printf(stringin);
	//printf("%s\n",stringin);
	printf("vvvvvvvv STEP %d ",step);
	//printf("vvvvvvvv abs coords vvvvvvvvvvvv ");
	//printf(stringin);
	printf("(%s)",stringin);
	printf(" vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv\n");
	for (j=0;j<idxmax;j++) printf("%-5d %16.6f %16.6f %16.6f\n",j,pos[j].x*faktor,pos[j].y*faktor,pos[j].z*faktor);
	//printf("%s\n",stringin);
	//printf("^^^^^^^^ abs coords ^^^^^^^^^^^^ ");
	printf("^^^^^^^^ STEP %d ",step);
	printf("%s",stringin);
	printf(" ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^\n\n");
	//printf(stringin,"-------------------------------- pos 0 end ---------------------------\n");
}

void print_pos_forces_to_screen(int idxmax,char stringin[],double faktor,int step) {
    int j;
    if (idxmax > atoms)
        idxmax = atoms;
	printf("STEP %d ",step);
	////printf("vvvvvvvv abs coords vvvvvvvvvvvv ");
	//printf("\"%s\"",stringin);
	//printf(" vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv\n");
	printsep(stringin);
	for (j=0;j<idxmax;j++) printf("%d\t%.9f\t%.9f\t%.9f\t%.9f\t%.9f\t%.9f\n",j,pos[j].x*faktor,pos[j].y*faktor,pos[j].z*faktor,force[j].x*faktor_force,force[j].y*faktor_force,force[j].z*faktor_force);
	printf("STEP %d ",step);
	//printf("vvvvvvvv STEP %d ",step);
	//printf("^^^^^^^^ abs coords ^^^^^^^^^^^^ ");
	//printf("%s",stringin);
	//printf("\"%s\"",stringin);
	printsep(stringin);
	//printf(" ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^..............................................\n\n");
}

void print_pos_forcestmp_to_screen(int idxmax,char stringin[],double faktor,int step) {
    int j;
    if (idxmax > atoms)
        idxmax = atoms;
	printf("STEP %d ",step);
	////printf("vvvvvvvv abs coords vvvvvvvvvvvv ");
	//printf("\"%s\"",stringin);
	//printf(" vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv\n");
	printsep(stringin);
	for (j=0;j<idxmax;j++) printf("%d\t%.9f\t%.9f\t%.9f\t%.9f\t%.9f\t%.9f\n",j,pos[j].x*faktor,pos[j].y*faktor,pos[j].z*faktor,forcetmp[j].x*faktor_force,forcetmp[j].y*faktor_force,forcetmp[j].z*faktor_force);
	printf("STEP %d ",step);
	//printf("vvvvvvvv STEP %d ",step);
	//printf("^^^^^^^^ abs coords ^^^^^^^^^^^^ ");
	//printf("%s",stringin);
	//printf("\"%s\"",stringin);
	printsep(stringin);
	//printf(" ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^..............................................\n\n");
}

void print_vel_velff_to_screen(int idxmax,char stringin[],int step) {
    int j;
    if (idxmax > atoms)
        idxmax = atoms;
	printf("STEP %d ",step);
	////printf("vvvvvvvv abs coords vvvvvvvvvvvv ");
	//printf("\"%s\"",stringin);
	//printf(" vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv\n");
	printsep(stringin);
	for (j=0;j<idxmax;j++) printf("%d\t%.9f\t%.9f\t%.9f\t%.9f\t%.9f\t%.9f\n",j,vel[j].x,vel[j].y,vel[j].z,velff[j].x,velff[j].y,velff[j].z);
	printf("STEP %d ",step);
	//printf("vvvvvvvv STEP %d ",step);
	//printf("^^^^^^^^ abs coords ^^^^^^^^^^^^ ");
	//printf("%s",stringin);
	//printf("\"%s\"",stringin);
	printsep(stringin);
	//printf(" ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^..............................................\n\n");
}

void print_vel_to_screen(int idxmax,char stringin[],int step) {
    int j;
    if (idxmax > atoms)
        idxmax = atoms;
	printf("STEP %d ",step);
	////printf("vvvvvvvv abs coords vvvvvvvvvvvv ");
	//printf("\"%s\"",stringin);
	//printf(" vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv\n");
	printsep(stringin);
	for (j=0;j<idxmax;j++) printf("%d\t%.9f\t%.9f\t%.9f\n",j,vel[j].x,vel[j].y,vel[j].z);
	printf("STEP %d ",step);
	//printf("vvvvvvvv STEP %d ",step);
	//printf("^^^^^^^^ abs coords ^^^^^^^^^^^^ ");
	//printf("%s",stringin);
	//printf("\"%s\"",stringin);
	printsep(stringin);
	//printf(" ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^..............................................\n\n");
}

void print_forces_to_screen(int idxmax,char stringin[],double faktor,int step) {
    int j;
    if (idxmax > atoms)
        idxmax = atoms;
	printf("STEP %d ",step);
	printsep(stringin);
	for (j=0;j<idxmax;j++) printf("%d\t%.9f\t%.9f\t%.9f\n",j,force[j].x*faktor_force,force[j].y*faktor_force,force[j].z*faktor_force);
	printf("STEP %d ",step);
	printsep(stringin);
}

void print_velocities_to_screen(int idxmax,char stringin[],double faktor,int step) {
    int j;
    if (idxmax > atoms)
        idxmax = atoms;
	printf("STEP %d ",step);
	printsep(stringin);
	for (j=0;j<idxmax;j++) printf("%d\t%.9f\t%.9f\t%.9f\n",j,vel[j].x,vel[j].y,vel[j].z);
	printf("STEP %d ",step);
	printsep(stringin);
}

void print_forcesharm_to_screen(int idxmax,char stringin[],double faktor,int step) {
    int j;
    if (idxmax > atoms)
        idxmax = atoms;
	printf("STEP %d ",step);
	////printf("vvvvvvvv abs coords vvvvvvvvvvvv ");
	//printf("\"%s\"",stringin);
	//printf(" vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv\n");
	printsep(stringin);
	//for (j=0;j<idxmax;j++) printf("%d\t%.9f\t%.9f\t%.9f\t%.9f\t%.9f\t%.9f\n",j,pos[j].x*faktor,pos[j].y*faktor,pos[j].z*faktor,force[j].x*faktor_force,force[j].y*faktor_force,force[j].z*faktor_force);
	for (j=0;j<idxmax;j++) printf("%d\t%.9f\t%.9f\t%.9f\n",j,forceharm[j].x*faktor,forceharm[j].y*faktor,forceharm[j].z*faktor);
	printf("STEP %d ",step);
	//printf("vvvvvvvv STEP %d ",step);
	//printf("^^^^^^^^ abs coords ^^^^^^^^^^^^ ");
	//printf("%s",stringin);
	//printf("\"%s\"",stringin);
	printsep(stringin);
	//printf(" ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^..............................................\n\n");
}

void print_pos_forcesdft_to_screen(int idxmax,char stringin[],double faktor,int step) {
    int j;
    if (idxmax > atoms)
        idxmax = atoms;
	printf("STEP %d ",step);
	////printf("vvvvvvvv abs coords vvvvvvvvvvvv ");
	//printf("\"%s\"",stringin);
	//printf(" vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv\n");
	printsep(stringin);
	for (j=0;j<idxmax;j++) printf("%d\t%.9f\t%.9f\t%.9f\t%.9f\t%.9f\t%.9f\n",j,pos[j].x*faktor,pos[j].y*faktor,pos[j].z*faktor,forcedft[j].x,forcedft[j].y,forcedft[j].z);
	printf("STEP %d ",step);
	//printf("vvvvvvvv STEP %d ",step);
	//printf("^^^^^^^^ abs coords ^^^^^^^^^^^^ ");
	//printf("%s",stringin);
	//printf("\"%s\"",stringin);
	printsep(stringin);
	//printf(" ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^..............................................\n\n");
}

void print_forcesdft_to_screen(int idxmax,char stringin[],double faktor,int step) {
    int j;
    if (idxmax > atoms)
        idxmax = atoms;
	printf("STEP %d ",step);
	////printf("vvvvvvvv abs coords vvvvvvvvvvvv ");
	//printf("\"%s\"",stringin);
	//printf(" vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv\n");
	printsep(stringin);
	//for (j=0;j<idxmax;j++) printf("%d\t%.9f\t%.9f\t%.9f\t%.9f\t%.9f\t%.9f\n",j,pos[j].x*faktor,pos[j].y*faktor,pos[j].z*faktor,force[j].x*faktor_force,force[j].y*faktor_force,force[j].z*faktor_force);
	//for (j=0;j<idxmax;j++) printf("%d\t%.9f\t%.9f\t%.9f\n",j,forcedft[j].x*faktor_force,forcedft[j].y*faktor_force,forcedft[j].z*faktor_force);  WRONG FOR DFT FORCES WE DONT NEED THE FACTOR
	for (j=0;j<idxmax;j++) printf("%d\t%.9f\t%.9f\t%.9f\n",j,forcedft[j].x,forcedft[j].y,forcedft[j].z);
	printf("STEP %d ",step);
	//printf("vvvvvvvv STEP %d ",step);
	//printf("^^^^^^^^ abs coords ^^^^^^^^^^^^ ");
	//printf("%s",stringin);
	//printf("\"%s\"",stringin);
	printsep(stringin);
	//printf(" ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^..............................................\n\n");
}


////////////// write to file
void write_out_info(double T,double dt,int zeitschritte, int l,int read_pos,int read_pos0,int read_uoutcar,int verbose,int write_analyze,const char *filename_in_positions,const char *filename_in_positions0, const char *filename_in_hesse,int read_hesse,int read_forces,int columns, int columns0, const char *filename_in_uoutcar,int evolve_md_on_hesse,int write_positions_forces,int write_positions, int write_positions_rel) {
	FILE *anmerkung;
	int dtorig=dt*1e12;
	//@@char buff[20];
	//@@struct tm *sTm;
	//@@time_t now = time (0);
	//@@//sTm = gmtime (&now);
	//@@sTm = localtime (&now);
	//@@strftime (buff, sizeof(buff), "%Y-%m-%d %H:%M:%S", sTm);
	//@@printf ("PROGRAM START  : %s\n", buff);
	anmerkung=(FILE *)fopen("out_info.txt","wb");
	fprintf(anmerkung,"./a.out %g %g %d %d %d %d %d\n", T, dt,zeitschritte,l,read_pos,verbose,write_analyze);
	fprintf(anmerkung,"\n");
	//@@fprintf(anmerkung,"PROGRAM START  : %s\n", buff);
	fprintf(anmerkung,"N              : %d\n",N);
	fprintf(anmerkung,"alat_lattice   : %-20.10f (~4.14) \n",alat_lattice);
	fprintf(anmerkung,"alat_morse     : %-20.10f (~4.14) (~4.14)  \n",alat_morse);
	fprintf(anmerkung,"d0             : %-10.10f\n",d0);
	fprintf(anmerkung,"a0             : %g (4294967296./N)\n",a0);
	fprintf(anmerkung,"r0_mor         : %g\n",r0_mor);
	fprintf(anmerkung,"r0_eq_morse    : %g\n",r0_eq_morse);
	fprintf(anmerkung,"a_mor          : %-20.10f (~1.43)\n",a_mor);
	fprintf(anmerkung,"D_mor          : %-20.10f (~0.26)\n",D_mor);
	fprintf(anmerkung,"ktr_tox        : %-20.10f\n",ktr_tox);
	fprintf(anmerkung,"m_element      : %-20.10f\n",m_element/1.660538e-27);

	fprintf(anmerkung,"T              : %g\n",T);
	fprintf(anmerkung,"dt             : %g\n",dt);
	fprintf(anmerkung,"zeitschritte   : %d\n",zeitschritte);
	fprintf(anmerkung,"l              : %d\n",l);

	fprintf(anmerkung,"\n");
	if (read_pos==1){fprintf(anmerkung,"read_pos       : %d (True)\n",read_pos);}
    if (read_pos!=1){fprintf(anmerkung,"read_pos       : %d (False)\n",read_pos);}
	if (read_pos0==1){fprintf(anmerkung,"read_pos0      : %d (True)\n",read_pos0);}
    if (read_pos0!=1){fprintf(anmerkung,"read_pos0      : %d (False)\n",read_pos0);}
	if (read_uoutcar==1){fprintf(anmerkung,"read_uoutcar   : %d (True)\n" ,read_uoutcar);}
    if (read_uoutcar!=1){fprintf(anmerkung,"read_uoutcar   : %d (False)\n",read_uoutcar);}
	fprintf(anmerkung,"filename_in_positions : %s\n",filename_in_positions);
	fprintf(anmerkung,"filename_in_positions0: %s\n",filename_in_positions0);
	fprintf(anmerkung,"filename_in_hesse     : %s\n",filename_in_hesse);
	fprintf(anmerkung,"filename_in_uoutcar   : %s\n",filename_in_uoutcar);
	if (verbose>=1) {fprintf(anmerkung,"verbose        : %d (True)\n",verbose);}
	if (verbose<=0) {fprintf(anmerkung,"verbose        : %d (False)\n",verbose);}
	if (write_analyze==1){fprintf(anmerkung,"write_analyze  : %d (True)\n",write_analyze);}
    if (write_analyze!=1){fprintf(anmerkung,"write_analyze  : %d (False)\n",write_analyze);}
	if (read_hesse==1){fprintf(anmerkung,"read_hesse     : %d (True)\n",read_hesse);}
    if (read_hesse!=1){fprintf(anmerkung,"read_hesse     : %d (False)\n",read_hesse);}
	if (read_forces==1){fprintf(anmerkung,"read_forces    : %d (True)\n",read_forces);}
    if (read_forces!=1){fprintf(anmerkung,"read_forces    : %d (False)\n",read_forces);}
	fprintf(anmerkung,"columns        : %d (read in from positionsfile [3 or 6])\n",columns);
	fprintf(anmerkung,"columns0       : %d (read in from positions0file [3 or 6])\n",columns0);
	if (evolve_md_on_hesse==1){fprintf(anmerkung,"evolve_md_on_hesse: %d (True)\n" ,evolve_md_on_hesse);}
    if (evolve_md_on_hesse!=1){fprintf(anmerkung,"evolve_md_on_hesse: %d (False)\n",evolve_md_on_hesse);}
	if (write_positions_forces==1){fprintf(anmerkung,"write_positions_forces: %d (True)\n" ,write_positions_forces);}
    if (write_positions_forces!=1){fprintf(anmerkung,"write_positions_forces: %d (False)\n",write_positions_forces);}
	if (write_positions==1){fprintf(anmerkung,"write_positions: %d (True)\n" ,write_positions);}
    if (write_positions!=1){fprintf(anmerkung,"write_positions: %d (False)\n",write_positions);}

	if (write_positions_rel==1){fprintf(anmerkung,"write_positions_rel: %d (True)\n" ,write_positions_rel);}
    if (write_positions_rel!=1){fprintf(anmerkung,"write_positions_rel: %d (False)\n",write_positions_rel);}
	fclose(anmerkung);


	// print to screen
	printsep("");
	printf("N              : %d\n",N);
	printf("alat_lattice   : %-20.10f (~4.14)\n",alat_lattice);
	printf("alat_morse     : %-20.10f (~4.14)\n",alat_morse);
	printf("d0             : %g (3037000500./N*alat_morse/alat_lattice) \n",d0);
	printf("a0             : %g (4294967296./N)==bins per a0\n",a0);
	printf("r0_mor         : %g\n",r0_mor);
	printf("r0_eq_morse    : %g\n",r0_eq_morse);
	printf("a_mor          : %-20.10f (~1.43)\n",a_mor);
	printf("D_mor          : %-20.10f (~0.26)\n",D_mor);
	printf("ktr_tox        : %-20.10f\n",ktr_tox);
	printf("m_element      : %-20.10f\n",m_element/1.660538e-27);
	printf("\n");
	printf("T              : %g\n",T);
	printf("dt             : %g\n",dt);
	printf("zeitschritte   : %d\n",zeitschritte);
	printf("l              : %d\n",l);
	printf("\n");
	printf("filename_in_positions : %s\n",filename_in_positions);
	printf("filename_in_positions0: %s\n",filename_in_positions0);
	printf("filename_in_hesse     : %s\n",filename_in_hesse);
	printf("filename_in_uoutcar   : %s\n",filename_in_uoutcar);
	if (read_pos==1){printf("read_pos       : %d (True)\n",read_pos);}
	if (read_pos!=1){printf("read_pos       : %d (False)\n",read_pos);}
	if (read_pos0==1){printf("read_pos0      : %d (True)\n",read_pos0);}
	if (read_pos0!=1){printf("read_pos0      : %d (False)\n",read_pos0);}
	if (read_uoutcar==1){printf("read_uoutcar   : %d (True)\n" ,read_uoutcar);}
	if (read_uoutcar!=1){printf("read_uoutcar   : %d (False)\n",read_uoutcar);}
	if (verbose>=1) {printf("verbose        : %d (True)\n",verbose);}
	if (verbose<=0) {printf("verbose        : %d (False)\n",verbose);}
	if (write_analyze==1) {printf("write_analyze  : %d (True)\n",write_analyze);}
	if (write_analyze!=1) {printf("write_analyze  : %d (False)\n",write_analyze);}
	if (read_hesse==1){printf("read_hesse     : %d (True)\n",read_hesse);}
    if (read_hesse!=1){printf("read_hesse     : %d (False)\n",read_hesse);}
	if (read_forces==1){printf("read_forces    : %d (True)\n",read_forces);}
    if (read_forces!=1){printf("read_forces    : %d (False)\n",read_forces);}
	printf("columns        : %d (read in from positionsfile [3 or 6])\n",columns);
	printf("columns0       : %d (read in from position0sfile [3 or 6])\n",columns0);
	if (evolve_md_on_hesse==1) {printf("evolve_md_on_hesse: %d (True)\n" ,evolve_md_on_hesse);}
	if (evolve_md_on_hesse!=1) {printf("evolve_md_on_hesse: %d (False)\n",evolve_md_on_hesse);}
	if (write_positions_forces==1) {printf("write_positions_forces: %d (True)\n" ,write_positions_forces);}
	if (write_positions_forces!=1) {printf("write_positions_forces: %d (False)\n",write_positions_forces);}
	if (write_positions==1) {printf("write_positions: %d (True)\n" ,write_positions);}
	if (write_positions!=1) {printf("write_positions: %d (False)\n",write_positions);}
	if (write_positions_rel==1) {printf("write_positions_rel: %d (True)\n" ,write_positions_rel);}
	if (write_positions_rel!=1) {printf("write_positions_rel: %d (False)\n",write_positions_rel);}
	printsep("");
	//printf("\n");
	//printf("\n");
	//printf("\n");
}

void write_out_info_add_timing(double start, double end,char stringin[]) {
	FILE *anmerkung;
	char buff[20];
	struct tm *sTm;
	time_t now = time (0);
	sTm = localtime (&now);
	strftime (buff, sizeof(buff), "%Y-%m-%d %H:%M:%S", sTm);
	printf ("PROGRAM STOP   : %s\n", buff);
	anmerkung=(FILE *)fopen("out_info.txt","a");
	fprintf(anmerkung,"PROGRAM %s  : %s\n",stringin, buff);

	fprintf(stderr,"#cpu-time %s     : %f sec\n",stringin,((double) (end - start)) / CLOCKS_PER_SEC);
	fprintf(anmerkung,"#cpu-time      : %f sec\n",((double) (end - start)) / CLOCKS_PER_SEC);
	fclose(anmerkung);
}

void write_out_pos0_to_out_EqCoords_direct() {
    int j;
	FILE *currentfile;
	currentfile=(FILE *)fopen("out_EqCoords_direct","wb");
	for (j=0;j<atoms;j++) fprintf(currentfile,"%g %g %g\n",
	    pos0[j].x/a0/N,pos0[j].y/a0/N,pos0[j].z/a0/N);
	fclose(currentfile);
}

void write_out_analyze_positions_distances_from_equilibrium_xyzmean() {
    int j;
	FILE *currentfile;
	currentfile=(FILE *)fopen("out_analyze_positions_distances_from_equilibrium_x_average.dat","wb");
	//for (j=0;j<atoms;j++) fprintf(currentfile,"%d %2.4f\n",j,sqrtdudumean[j].x);
	for (j=0;j<atoms;j++) fprintf(currentfile,"%2.4f\n",sqrtdudumean[j].x);
	fclose(currentfile);
	currentfile=(FILE *)fopen("out_analyze_positions_distances_from_equilibrium_y_average.dat","wb");
	//for (j=0;j<atoms;j++) fprintf(currentfile,"%d %2.4f\n",j,sqrtdudumean[j].y);
	for (j=0;j<atoms;j++) fprintf(currentfile,"%2.4f\n",sqrtdudumean[j].y);
	fclose(currentfile);
	currentfile=(FILE *)fopen("out_analyze_positions_distances_from_equilibrium_z_average.dat","wb");
	//for (j=0;j<atoms;j++) fprintf(currentfile,"%d %2.4f\n",j,sqrtdudumean[j].z);
	for (j=0;j<atoms;j++) fprintf(currentfile,"%2.4f\n",sqrtdudumean[j].z);
	fclose(currentfile);
}

void write_out_cell() {
	FILE *currentfile;
	currentfile=(FILE *)fopen("cell","wb");
	//fprintf(currentfile," %g %g %g\n", pos0[0].x/a0/N,pos0[1].y/a0/N,pos0[2].z/a0/N);
	fprintf(currentfile,"%g %g %g\n", alat_lattice*N,0.0,0.0);
	fprintf(currentfile,"%g %g %g\n", 0.0,alat_lattice*N,0.0);
	fprintf(currentfile,"%g %g %g\n", 0.0,0.0,alat_lattice*N);
	fclose(currentfile);
}

void write_positions_forces_tofile(FILE *file_out_positions) {
    int j=0;
    for (j=0;j<atoms;j++) fprintf(file_out_positions,"%14.10f %14.10f %14.10f %14.10f %14.10f %14.10f\n",pos[j].x*faktor,pos[j].y*faktor,pos[j].z*faktor,force[j].x/1.602176e-9,force[j].y/1.602176e-9,force[j].z/1.602176e-9);
}

void write_positions_tofile(FILE *file_out_positions, int write_positions_rel) {
    int j=0;
	if (write_positions_rel==1){
        for (j=0;j<atoms;j++) fprintf(file_out_positions,"%14.10f %14.10f %14.10f\n",pos[j].x*faktor_rel_pos,pos[j].y*faktor_rel_pos,pos[j].z*faktor_rel_pos);}
    else {
        for (j=0;j<atoms;j++) fprintf(file_out_positions,"%14.10f %14.10f %14.10f\n",pos[j].x*faktor,pos[j].y*faktor,pos[j].z*faktor);}
}

void write_dudl_head(FILE *file_out_dudl,FILE *file_out_dudl_av) {
       fprintf(file_out_dudl,"#step          DFT-LA     DFT-H    temp(K)        u_DFT         u_la        u_harm\n");
    fprintf(file_out_dudl_av,"#step  std<DFT-LA> dudl<DFT-LA>  std<DFT-H> dudl<DFT-H>  <temp(K)>  <u_DFT>  <u_la> <u_harm> <u_h_perf> <for DFT-LA> <for DFT-H>\n");
}

void write_forces_head(FILE *file_forces,FILE *file_forces_av) {
    fprintf(file_forces,"#step*atoms*3  forces(DFT-LA) std<forces_diffmax>\n");
    fprintf(file_forces_av,"#step*atoms*3  forces<DFT-LA> \n");
}


void write_analyze_forces_vs_dft_head(FILE *file) {
       fprintf(file,"#f(DFT) f(LA) f(DFT-LA) <f(DFT-LA)>\n");
}

void write_dudl(int step, FILE *file_out_dudl,FILE *file_out_dudl_av,int read_hesse,int read_uoutcar) {
    // step, u,
    //printf("p:a\n");
    //printf("p:a %d\n",step);
    //double temperature_av=temperature_sum/step;
    //if (step==0){temperature_av=temperature_sum;};
    double urefclassical_av=1.5*0.086173423*temperature_av*(atoms/(atoms-1.)); // 3/2 kB T   * (32/31)
    //printf("step: %-10d in  u_la_av: %8.4f\n",step,u_la_av);
    u_la_av   = getaverage(u_la_av,step,u_la);
    u_dft_av  = getaverage(u_dft_av,step,u_dft);
    //printf("step: %-10d out u_la_av: %8.4f\n",step,u_la_av);

    dudl_dft_la        = u_dft - u_la;
	if (read_uoutcar==0){dudl_dft_la=0;}

    dudl[step].dudl_dft_la = dudl_dft_la;

    dudl_dft_la_av   = getaverage(dudl_dft_la_av,step,dudl_dft_la);
    //printf("step: %-10d dudl_dft_la_av: %8.4f\n",step,dudl_dft_la_av);

    // calculation of the variance and standard deviation
    double sum_var_la=0;
    double sum_var_harm=0;
    double diff;
    int j;
    for (j=0;j<=step;j++) {
        diff=(dudl[j].dudl_dft_la - dudl_dft_la_av);
        sum_var_la += diff*diff;
        printf("step--:%-10d j:%-10d dudl[j]:%8.4f diff:%8.4f diff**2:%8.4f sqrt(sum_var_la):%8.4f\n",step,j,dudl[j].dudl_dft_la,diff,diff*diff,sqrt(sum_var_la/step));
        if (read_hesse==1) {
            diff=(dudl[j].dudl_dft_harm - dudl_dft_harm_av);
            sum_var_harm += diff*diff;
        } else {
            sum_var_harm = 0;
        }
    }
    //printf("step: %-10d sum_var_la: %8.4f\n",step,sum_var_la);
    dudl[step].dudl_dft_la_std      = sqrt(sum_var_la/(step)); //-0.9999));    is equivalent to michaels code
    dudl[step].dudl_dft_harm_std    = sqrt(sum_var_harm/(step)); //-0.9999));

    if (read_hesse==1) {
        u_harm_av = getaverage(u_harm_av,step,u_harm);
        dudl_dft_harm      = u_dft - u_harm;
        dudl[step].dudl_dft_harm = dudl_dft_harm;
        dudl_dft_harm_av = getaverage(dudl_dft_harm_av,step,dudl_dft_harm);
    } else {
        u_harm_av = getaverage(u_harm_av,step,u_harm);
        dudl_dft_harm      = 0;
        dudl[step].dudl_dft_harm = 0;
        dudl_dft_harm_av = 0;
    }

    fprintf(file_out_dudl,   "%-10d  %8.4f   %8.4f  %8.1f     %8.4f      %8.4f     %8.2f\n",\
            step,\
            dudl_dft_la,\
            dudl_dft_harm,\
            temperature,\
            u_dft,\
            u_la,\
            u_harm
            );

    //printf("p:step\n");
    //
    //leave the %8.5f %8.5 for the forces! (was necessary in LA case once)
    fprintf(file_out_dudl_av,"%-8d %8.4f   %8.4f     %8.2f  %8.2f     %8.1f  %8.2f %8.2f %8.2f %8.2f %8.5f %8.5f\n",\
            step,\
            dudl[step].dudl_dft_la_std,\
            dudl_dft_la_av,\
            dudl[step].dudl_dft_harm_std,\
            dudl_dft_harm_av,\
            temperature_av,\
            u_dft_av,\
            u_la_av,\
            u_harm_av,\
            urefclassical_av,\
            sqrt(fdiff_dft_la_av),\
            sqrt(fdiff_dft_harm_av)
           );

    //fprintf(file_forces_av,"%-8d    %8.6f        %8.6f\n",\
    //        step,\
    //        forces_diffmax,\
    //        fdiff_dft_la_av
    //        );
    //printf("p:j\n");
}

void analyze_forces(int i,FILE *file_forces,FILE *file_forces_av,FILE *file_forces_vs_forces_dft, int write_analyze_forces, int read_hesse) {
    // - amount of lines: number_of_atoms*3*steps_total
    // - first  line: force in x: DFT,LA
    // - second line: force in y
    //
    // fuer eine std brauchen wir noch:
    // sum1 += sqrt(dftforces^2) und davon die summe ueber alle eintraege
    // sum2 +=sqrt(laforces^2) und davon die summe ueber alle eintraege
    // sum1/eintraege

    //for (i=0;i<atoms*first_neighbors/2;i++) {   // 12 * atoms_supercell / 2
    //    ind1=l1nn[i].ind1;
    //    ind2=l1nn[i].ind2;
    //    j1=l1nn[i].j1;
    //    j2=l1nn[i].j2;
    //    j3=l1nn[i].j3;

    //    // abstaende nn
	//	x=(signed)(pos[ind1].x-pos[ind2].x);
	//	y=(signed)(pos[ind1].y-pos[ind2].y);
	//	z=(signed)(pos[ind1].z-pos[ind2].z);
	//	r=sqrt(x*x+y*y+z*z);
    //};

    int j;
    int iat3;
    double forcelax,forcelay,forcelaz,dx,dy,dz,dx_abs,dy_abs,dz_abs,ddx,ddy,ddz,ddx_abs,ddy_abs,ddz_abs,diffmax;
    //int anz=i*atoms*3;
	//printf("%.10f\n",2.0);
	forces_diffmax = 0;

	for (j=0;j<atoms;j++) { // goes through all atoms
        forcelax = force[j].x*faktor_force; // force la
        forcelay = force[j].y*faktor_force; // force la
        forcelaz = force[j].z*faktor_force; // force la
	    //printf("--> FDFT, FLA, FPL: %5.4f  %5.4f %5.4f\n",forcedft[j].x,forcelax,forcepoly[j].x*faktor_force);
        //dx = fabs(forcedft[j].x-x);  // difference in forces, needs actually be saved for every step
        //dy = fabs(forcedft[j].y-y);
        //dz = fabs(forcedft[j].z-z);
        //ddx = fabs(forcedft[j].x-forceharm[j].x*faktor);
        //ddy = fabs(forcedft[j].y-forceharm[j].y*faktor);
        //ddz = fabs(forcedft[j].z-forceharm[j].z*faktor);

        dx  = forcedft[j].x-forcelax;  // difference in forces, needs actually be saved for every step
        dy  = forcedft[j].y-forcelay;
        dz  = forcedft[j].z-forcelaz;

        // get forces_diffmax
        dx_abs = fabs(dx);
        dy_abs = fabs(dy);
        dz_abs = fabs(dz);
        if (dx_abs > forces_diffmax) {forces_diffmax = dx_abs;};
        if (dy_abs > forces_diffmax) {forces_diffmax = dy_abs;};
        if (dz_abs > forces_diffmax) {forces_diffmax = dz_abs;};
        //dx = dx*dx;
        //dy = dy*dy;
        //dz = dz*dz;

	    if (read_hesse==1) {
            ddx = forcedft[j].x-forceharm[j].x*faktor;
            ddy = forcedft[j].y-forceharm[j].y*faktor;
            ddz = forcedft[j].z-forceharm[j].z*faktor;
            ddx_abs = fabs(ddx);
            ddy_abs = fabs(ddy);
            ddz_abs = fabs(ddz);
            //ddx = ddx*ddx;
            //ddy = ddy*ddy;
            //ddz = ddz*ddz;
            // currently unused
            //fdudl_dft_harm_av = getaverage(fdudl_dft_harm_av,i,ddx);
            //fdudl_dft_harm_av = getaverage(fdudl_dft_harm_av,i,ddy);
            //fdudl_dft_harm_av = getaverage(fdudl_dft_harm_av,i,ddz);
	    } else {
	        ddx = 0;
	        ddy = 0;
	        ddz = 0;
	        ddx_abs = 0;
	        ddy_abs = 0;
	        ddz_abs = 0;
            fdiff_dft_harm_av   = 0;
	    }

        iat3 = i*atoms*3;
	    //printf("diffmax %.10f dx \%.10f\n",forces_diffmax,dx);
        //fprintf(fab,"%5.3f %5.3f %5.3f\n",dx,dy,dz);
        // write to file if necessary
        fdiff_dft_la_av    = getaverage(fdiff_dft_la_av,  iat3+j*3+0,dx*dx); // since <x> = 0 : (x-<x>) = x
        fdiff_dft_la_av    = getaverage(fdiff_dft_la_av,  iat3+j*3+1,dy*dy);
        fdiff_dft_la_av    = getaverage(fdiff_dft_la_av,  iat3+j*3+2,dz*dz);

        fdiff_dft_harm_av  = getaverage(fdiff_dft_harm_av,iat3+j*3,ddx*ddx); // since <x> = 0 : (x-<x>) = x
        fdiff_dft_harm_av  = getaverage(fdiff_dft_harm_av,iat3+j*3,ddy*ddy); // since <x> = 0 : (x-<x>) = x
        fdiff_dft_harm_av  = getaverage(fdiff_dft_harm_av,iat3+j*3,ddz*ddz); // since <x> = 0 : (x-<x>) = x


        fstd_dft       = getaverage(fstd_dft,     iat3+j*3,forcedft[j].x*forcedft[j].x);
        fstd_dft       = getaverage(fstd_dft,     iat3+j*3,forcedft[j].y*forcedft[j].y);
        fstd_dft       = getaverage(fstd_dft,     iat3+j*3,forcedft[j].z*forcedft[j].z);
        fstd_la        = getaverage(fstd_la ,     iat3+j*3,forcelax*forcelax);
        fstd_la        = getaverage(fstd_la ,     iat3+j*3,forcelay*forcelay);
        fstd_la        = getaverage(fstd_la ,     iat3+j*3,forcelaz*forcelaz);
        fcov_dft_la    = getaverage(fcov_dft_la,  iat3+j*3,forcedft[j].x*forcelax);
        fcov_dft_la    = getaverage(fcov_dft_la,  iat3+j*3,forcedft[j].y*forcelay);
        fcov_dft_la    = getaverage(fcov_dft_la,  iat3+j*3,forcedft[j].z*forcelaz);

        //printf("x i %d dx %.10f fdiff_dft_la_av %.10f
        //\n",iat3+j*3,dx,fdiff_dft_la_av);
	    if (write_analyze_forces==1) {
	        fprintf(file_forces_vs_forces_dft,"%5.4f  %5.4f\n",forcedft[j].x,forcelax);
	        fprintf(file_forces_vs_forces_dft,"%5.4f  %5.4f\n",forcedft[j].y,forcelay);
	        fprintf(file_forces_vs_forces_dft,"%5.4f  %5.4f\n",forcedft[j].z,forcelaz);

	        fprintf(file_forces   ,"%d  %5.5f\n",iat3+j*3+0,dx);
	        fprintf(file_forces   ,"%d  %5.5f\n",iat3+j*3+1,dy);
	        fprintf(file_forces   ,"%d  %5.5f\n",iat3+j*3+2,dz);

	        fprintf(file_forces_av,"%d  %5.5f\n",iat3+j*3+0,sqrt(fdiff_dft_la_av));
	        fprintf(file_forces_av,"%d  %5.5f\n",iat3+j*3+1,sqrt(fdiff_dft_la_av));
	        fprintf(file_forces_av,"%d  %5.5f\n",iat3+j*3+2,sqrt(fdiff_dft_la_av));
	    };

        //printf("y i %d dy %.10f fdiff_dft_la_av %.10f
        //\n",iat3+j*3+1,dy,fdiff_dft_la_av);
	    //if (write_analyze_forces==1) {
	    //};

        //printf("z i %d dz %.10f fdiff_dft_la_av %.10f
        //\n",iat3+j*3+2,dz,fdiff_dft_la_av);
	    //if (write_analyze_forces==1) {
	    //};
	}
	//exit(1);
}

void write_analyze_positions(int i) { //,FILE *file_out_check_dist_xyz,FILE *file_out_check_dist_r,FILE *file_out_check_dist_nn,FILE *file_out_check_dist_nn_proj,FILE *file_out_prl15_2a,FILE *file_out_prl15_2au) {
    int j=0;
    //double dx,dy,dz,rrx,rry,rrz,rr;
    double dx,dy,dz,rr,p1,p2,p3,projlen,dtrans;
    double x,y,z,r,x0,y0,z0,xn,yn,zn;
	int k,j1,j2,j3,ind1,ind2;
	////////////////////////////////
	// schleife ueber alle atome
	////////////////////////////////
		for (j=0;j<atoms;j++) {  // schleife ueber alle atome
		    dx = (pos[j].x-pos0[j].x)*da0_alat_lattice;  // sollte in angstrom sein
		    if (dx > N*alat_lattice/2)
		        dx = dx - N*alat_lattice;
		    dy = (pos[j].y-pos0[j].y)*da0_alat_lattice;
		    if (dy > N*alat_lattice/2)
		        dy = dy - N*alat_lattice;
		    dz = (pos[j].z-pos0[j].z)*a0_alat_lattice;
		    if (dz > N*alat_lattice/2)
		        dz = dz - N*alat_lattice;
		    // this would be the distances from equilibrium
            //printf("[dx,dy,dz] of atom %d: %.6f %.6f %.6f\n",j,dx,dy,dz);
            //printf("a (rel cell) pos  of atom %d: %.6f %.6f %.6f\n",j, pos[j].x/a0,              pos[j].y/a0,             pos[j].z/a0);
            //printf("b (rel cell) pos0 of atom %d: %.6f %.6f %.6f\n",j,pos0[j].x/a0,             pos0[j].y/a0,            pos0[j].z/a0);
            //printf("b (angstrom) pos? of atom %d: %.6f %.6f %.6f\n",j,(pos0[j].x-pos[j].x)/a0,(pos0[j].y-pos[j].y)/a0,(pos0[j].z-pos[j].z)/a0);
            //printf("c (angstrom) pos  of atom %d: %.6f %.6f %.6f\n",j, pos[j].x/a0*alat_lattice, pos[j].y/a0*alat_lattice,pos[j].z/a0*alat_lattice);
            //printf("d (angstrom) pos0 of atom %d: %.6f %.6f %.6f\n",j,pos0[j].x/a0*alat_lattice,pos0[j].y/a0*alat_lattice,pos0[j].z/a0*alat_lattice);
            //printf("e (angstrom) posd of atom %d: %.6f %.6f %.6f\n",j,(pos0[j].x-pos[j].x)/a0*alat_lattice,(pos0[j].y-pos[j].y)/a0*alat_lattice,(pos0[j].z-pos[j].z)/a0*alat_lattice);
		    fprintf(file_out_check_dist_xyz,"%.10f %.10f %.10f\n",dx,dy,dz);
            //rrx=sqrt(dx*dx);
            //rry=sqrt(dy*dy);
            //rrz=sqrt(dz*dz); // hiervon waere der mean fuer jedes atom interessant

            //if (j==0) {printf("a %2.5f %2.5f\n",sqrtdudumean[j].x,rrx);};
	        sqrtdudumean[j].x = getaverage(sqrtdudumean[j].x,i,fabs(dx));
            //if (j==0) {printf("b %2.5f\n",sqrtdudumean[j].x);};
	        sqrtdudumean[j].y = getaverage(sqrtdudumean[j].y,i,fabs(dy));
	        sqrtdudumean[j].z = getaverage(sqrtdudumean[j].z,i,fabs(dz));
			rr=sqrt(dx*dx+dy*dy+dz*dz);  // hiervon waere auch der mean interessant fuer jedes atom in x y und z
			//if (rr > alat_lattice/2) {  // is ca
			if (rr > alat_lattice/2) {
		        printf("warning step ?? atom %d atomdistance to orig: %.10f\n",j,rr);
                printf("pos  of atom %d: %.6f %.6f %.6f\n",j,pos[j].x/a0,pos[j].y/a0,pos[j].z/a0);
                printf("pos0 of atom %d: %.6f %.6f %.6f\n",j,pos0[j].x/a0,pos0[j].y/a0,pos0[j].z/a0);
                exit(1);
            }
		    fprintf(file_out_check_dist_r,"%.10f\n",rr);
        }

	////////////////////////////////
	// schleife ueber naechste nachbarn
	////////////////////////////////
    for (k=0;k<atoms*first_neighbors/2;k++) {   // 12 * atoms_supercell / 2
        ind1=l1nn[k].ind1;
        ind2=l1nn[k].ind2;
        j1=l1nn[k].j1;
        j2=l1nn[k].j2;
        j3=l1nn[k].j3;

        // abstaende nn
		x=(signed)(pos[ind1].x-pos[ind2].x)*da0_alat_lattice;
		y=(signed)(pos[ind1].y-pos[ind2].y)*da0_alat_lattice;
		z=(signed)(pos[ind1].z-pos[ind2].z)*da0_alat_lattice;
		xn=(signed)(pos[ind1].x-pos[ind2].x)*da0_alat_lattice;
		yn=(signed)(pos[ind1].y-pos[ind2].y)*da0_alat_lattice;
		zn=(signed)(pos[ind1].z-pos[ind2].z)*da0_alat_lattice;
		xn = sqrt(xn*xn);
		yn = sqrt(yn*yn);
		zn = sqrt(zn*zn);
		x0=(signed)(pos0[ind1].x-pos0[ind2].x)*da0_alat_lattice;
		y0=(signed)(pos0[ind1].y-pos0[ind2].y)*da0_alat_lattice;
		z0=(signed)(pos0[ind1].z-pos0[ind2].z)*da0_alat_lattice;
        //printf("nn0 (angstrom) pos? of atom %d:%d:%d %.6f %.6f %.6f\n",k,ind1,ind2,x0,y0,z0);
        //printf("nn! (angstrom) pos? of atom %d:%d:%d %.6f %.6f %.6f\n",k,ind1,ind2,x,y,z);
        //printf("nn1 (angstrom) pos? of atom %d:%d:%d %.6f %.6f %.6f\n",k,ind1,ind2,xn,yn,zn);
        //printf("nn2 (angstrom) pos? of atom %d:%d:%d %.6f %.6f %.6f\n",k,ind1,ind2,xn,yn,zn);
		r=sqrt(x*x+y*y+z*z);
		fprintf(file_out_check_dist_nn,"%.10f\n",r);

        // get projection of x,y,z on x0,y0,z0
		projection(x,y,z,x0,y0,z0);
		projlen = sqrt(projx*projx+projy*projy+projz*projz);
		fprintf(file_out_check_dist_nn_proj,"%.10f\n",projlen);  // v = x,y,z in getDistance.math

		dtrans=sqrt((x-projx)*(x-projx)+(y-projy)*(y-projy)+(z-projz)*(z-projz));
		fprintf(file_out_new1,"%.10f\n",dtrans);  // v = x,y,z in getDistance.math
		fprintf(file_out_new2,"%.10f %.10f\n",projlen,dtrans);  // v = x,y,z in getDistance.math

        // make the 2d plot for the prl
        if (x0==0) {
		    fprintf(file_out_prl15_2a,"%.10f %.10f\n",yn,zn);
		    fprintf(file_out_prl15_2a,"%.10f %.10f\n",zn,yn);
		    //fprintf(file_out_prl15_2au,"%.10f %.10f\n",z,y);
		    //printf("x0=0 step %d:%d:%d %.10f %.10f %.10f\n",k,ind1,ind2,x,y,z);
	        //printf("xxxx %.10f %.10f %.10f %.10f %.10f %.10f\n",x,y,z,x0,y0,z0);
		    //printf("xxxy %.10f %.10f %.10f\n",projx,projy,projz);
		    //exit(1);
		    //fprintf(file_out_new1,"%.10f %.10f\n",projlen,xn);  // v = x,y,z in getDistance.math
		    //fprintf(file_out_new2,"%.10f %.10f\n",r,xn);  // v = x,y,z in getDistance.math
        }
        if (y0==0) {
		    fprintf(file_out_prl15_2a,"%.10f %.10f\n",xn,zn);
		    fprintf(file_out_prl15_2a,"%.10f %.10f\n",zn,xn);
		    //fprintf(file_out_prl15_2au,"%.10f %.10f\n",z,x);
		    //fprintf(file_out_prl15_2au,"%.10f %.10f\n",x,z);
		    //fprintf(file_out_new1,"%.10f %.10f\n",projlen,yn);  // v = x,y,z in getDistance.math
		    //fprintf(file_out_new2,"%.10f %.10f\n",r,yn);  // v = x,y,z in getDistance.math
        }
        if (z0==0) {
		    fprintf(file_out_prl15_2a,"%.10f %.10f\n",xn,yn);
		    fprintf(file_out_prl15_2a,"%.10f %.10f\n",yn,xn);
		    //fprintf(file_out_new1,"%.10f %.10f\n",projlen,zn);  // v = x,y,z in getDistance.math
		    //fprintf(file_out_new2,"%.10f %.10f\n",r,zn);  // v = x,y,z in getDistance.math
		    //fprintf(file_out_prl15_2au,"%.10f %.10f\n",xn,zn);
		    //fprintf(file_out_prl15_2au,"%.10f\n",z);
		    //printf("z0=0 %d:%d:%d %.10f %.10f %.10f\n",k,ind1,ind2,xn,yn,z);
		    //fprintf(file_out_prl15_2au,"%.10f %.10f\n",y,x);
        }
    }
}

void check_pos_vs_pos0_if_same_positions_and_check_if_correct_alat() {
	int j=0;
	int check=0;
	int wrongatom=9999999;
	//printf("checking now...\n");
	for (j=0;j<atoms;j++) {
	    //printf("aaa %d %d\n",j,check);
	    //printf("%d pos0.x pos.x %2.4f %2.4f %2.4f\n",j,pos0[j].x*faktor,pos[j].x*faktor,pos0[j].x-pos[j].x);
	    if (pos0[j].x-pos[j].x != 0) {printf("not zerox\n");exit(1);};
	    if (pos0[j].y-pos[j].y != 0) {printf("not zeroy\n");exit(1);};
	    if (pos0[j].z-pos[j].z != 0) {printf("not zeroz\n");exit(1);};
	    //printf("%d pos0.x alat_lattice %2.4f %2.4f\n",j,pos0[j].x*faktor,alat_lattice);
	    //if (pos0[j].x*faktor==alat_lattice) {check=1;wrongatom=j;};  // i dont get this line
	    //printf("kkk %d %d %2.4f %2.4f\n",j,check,pos0[j].x*faktor,alat_lattice);
    }
	//printf("check %d \n",check);
    if (check==1) {
        printf("pos and pos0 seem to have different positions ...atom %d\n",wrongatom);
        printf("pos and pos0 seem to have different positions ...atom %d\n",wrongatom);
        printf("pos and pos0 seem to have different positions ...atom %d\n",wrongatom);

        print_pos_0_to_screen(12, "POS_0 nach init");
        print_pos_to_screen(12, "POS",0);
    };
    if (check==1) {printf("check not 1!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n");exit(1);};
}

void check_pos0_for_duplicates(double faktor,const char *filename_in_positions) {
	int i,j=0;
	double a,b,c;
    //printf("check_pos0_for_duplicates\n");
	for (j=0;j<atoms-1;j++) {
	    a = pos0[j].x;
	    b = pos0[j].y;
	    c = pos0[j].z;
	    //printf("j: %d\ta: %.6f b: %.6f c: %.6f\n",j,a*faktor,b*faktor,c*faktor);
	    for (i=j+1;i<atoms;i++) {
	        //printf("j: %d a: %.6f b: %.6f c: %.6f\n",j,a,b,c);
	        //printf("%d %d\n",j,i);
	        if (pos0[i].x == pos0[j].x) {
	        if (pos0[i].y == pos0[j].y) {
	        if (pos0[i].z == pos0[j].z) {
	            printf("\n");
	            printf("j: %d\ta: %.6f b: %.6f c: %.6f\n",j,pos0[j].x*faktor,pos0[j].y*faktor,pos0[j].z*faktor);
	            printf("i: %d\ta: %.6f b: %.6f c: %.6f\n",i,pos0[i].x*faktor,pos0[i].y*faktor,pos0[i].z*faktor);
	            printf("\n");
	            printf("ERROR: Duplicates found! Exit! Either ... \n");
	            printf("       - the structure of the first pos is really bad or\n");
	            printf("       - somehow the positions have forces which was not recognized here\n");
	            printf("       - you scpcified a supercell size (N=%d) which is not consisten with the POSTIOINs file (\"%s\").\n",N,filename_in_positions);
                print_pos_0_to_screen(5, "A) pos_0 999");
	            printf("\n");
	            exit(1);
            }
            }
            }
        }
    }
}

void check_arguments(int i,char *argv[]) {
	if (i!=12) {
	    // 11 is the total number of input variables including filepath
	    printf("i is %d and not 12!\n",i);
	    //                   1   2      3  4   5   6  7  8  9  10  11 12
		printf("  usage: ./a.out T     dt  L   l   r  v  w  h  f   u  wpr   \n");
		printf("         ./a.out 1775 .001 120 20  0  0  0  0  0   0  0     \n\n");
        printf("    where:\n");
        printf("\n");
        printf("    T:  %-7s temperature (in K) corresponding to total energy per degree of freedom\n",argv[1]);
        printf("    dt: %-7s time step (in units of ps)\n",argv[2]);
        printf("    L:  %-7s number of time steps to simulate\n",argv[3]);
        printf("    l:  %-7s output positions every l simulation steps\n",argv[4]);
        printf("    r:  %-7s read positions (1) or not (any number but 1)\n",argv[5]);
        printf("    v:  %-7s verbose; 0==no;1=verbose;2==very_verbose\n",argv[6]);
        printf("    w:  %-7s write files for MD analyis; 1=True;else==False\n",argv[7]);
        printf("    h:  %-7s read hessematrix_sphinx   ; 1=True;else==False\n",argv[8]);
        printf("    f:  %-7s read forces               ; 1=True;else==False\n",argv[9]);
        printf("    u:  %-7s read u_OUTCAR             ; 1=True;else==False\n",argv[10]);
        printf("    u:  %-7s read u_OUTCAR             ; 1=True;else==False\n",argv[10]);
        printf("    wrp:%-7s write positions relative  ; 1=True;else==False\n",argv[11]);
        exit(1);
	}
}

void set_forces_energy_zero() {
	//force=(struct Force *)malloc(sizeof(struct Force)*atoms);
	int j=0;
	for (j=0;j<atoms;j++) {
	    force[j].x=0;
	    force[j].y=0;
	    force[j].z=0;
    }
    u_la=0;
    u_pl=0;
    u_la_tox=0;
}

void set_ref_forces_energy_zero() {
	//force=(struct Force *)malloc(sizeof(struct Force)*atoms);
	int j=0;
	for (j=0;j<atoms;j++) {
	    forceharm[j].x=0;
	    forceharm[j].y=0;
	    forceharm[j].z=0;
    }
    u_harm=0;
}


// read in
void get_pos_kfrom_externalfile(FILE *tmp) {



    int system(const char *command);

	double dumread[3];
	int j=0;
	for (j=0;j<atoms;j++) {
	    fscanf(tmp," %lf %lf %lf\n",dumread,dumread+1,dumread+2);
	    //pos[j].x=a0*(dumread[0]/alat_lattice); // RICHTIG  (Absolute -> internal)
	    //pos[j].y=a0*(dumread[1]/alat_lattice); // RICHTIG
	    //pos[j].z=a0*(dumread[2]/alat_lattice); // RICHTIG
	    pos[j].x=a0_alat_lattice*dumread[0]; // RICHTIG  (Absolute -> internal)
	    pos[j].y=a0_alat_lattice*dumread[1]; // RICHTIG
	    pos[j].z=a0_alat_lattice*dumread[2]; // RICHTIG
	    printf("\n");
	    printf("pp1:%2.6f %2.6f %2.6f\n",pos[j].x*faktor, pos[j].y*faktor,pos[j].z*faktor);
	    if (j==3) exit(1);
    }
}

//char* cmd_system(const char* command)
//{
//    char* result = "";
//    FILE *fpRead;
//    fpRead = popen(command, "r");
//    char buf[1024];
//    memset(buf,'\0',sizeof(buf));
//    while(fgets(buf,1024-1,fpRead)!=NULL)
//    {
//        result = buf;
//    }
//    if(fpRead!=NULL)
//       pclose(fpRead);
//    return result;
//}

int run_shellcommand_and_recieve_output(const char* command) {
    // e.g. call by         //columns = run_shellcommand_and_recieve_output("head -1 POSITIONs | wc -w");
    int output;
    FILE *fp = popen(command, "r");

    fscanf(fp, "%d", &output);
    pclose(fp);

    //printf("-------------kkkk\n");
    //printf("%d\n",output);
    //printf("-------------kkkk\n");
    return(output);
}

void get_pos_from_externalfile(FILE *tmp,int columns) {
    double dumread[columns];
    // !!!!!!!!! Version where MD is in correct order which probably is more suitable to adapt to read positions
    // /Users/glensk/Dropbox/proj/proj_current/__2017.01_phonon_linewidth_al/__2016.09_phonon_lifetimes_3_nach_elternzeit/michaels_code_adapted/sim_fcc_morse_with_forces.c
    //
    // !!!!!!!!!!!!! Version where read in positions works !!!!!!!!!!!!!!!!!!!!
    // /Users/glensk/Dropbox/proj/proj_current/__2017.01_phonon_linewidth_al/__2016.09_phonon_lifetimes_3_nach_elternzeit/michaels_code_simulation_morse_withforces/michaels_code_original_tox_is_exactly_as_in_python_samepositions_as_michael_random3_WORKS/forces_fcc_morse.c
    //
    //printf("col: %d\n",columns);
    //exit(1);
	int j=0;
	for (j=0;j<atoms;j++) {
	    if (columns==6) {
	        fscanf(tmp," %lf %lf %lf %lf %lf %lf\n",dumread,dumread+1,dumread+2,dumread+3,dumread+4,dumread+5);
	        //pos[j].x=a0*(dumread[0]/alat_lattice); // RICHTIG  (Absolute -> internal)
	        //pos[j].y=a0*(dumread[1]/alat_lattice); // RICHTIG
	        //pos[j].z=a0*(dumread[2]/alat_lattice); // RICHTIG
	        pos[j].x=a0_alat_lattice*dumread[0]; // RICHTIG  (Absolute -> internal)
	        pos[j].y=a0_alat_lattice*dumread[1]; // RICHTIG
	        pos[j].z=a0_alat_lattice*dumread[2]; // RICHTIG
	        forcedft[j].x=dumread[3];
	        forcedft[j].y=dumread[4];
	        forcedft[j].z=dumread[5];
            }
        else {
	        fscanf(tmp," %lf %lf %lf\n",dumread,dumread+1,dumread+2);
	        pos[j].x=a0_alat_lattice*dumread[0]; // RICHTIG  (Absolute -> internal)
	        pos[j].y=a0_alat_lattice*dumread[1]; // RICHTIG
	        pos[j].z=a0_alat_lattice*dumread[2]; // RICHTIG
        }

	    //printf("\n");
	    //printf("pp1:%2.6f %2.6f %2.6f\n",pos[j].x*faktor, pos[j].y*faktor,pos[j].z*faktor);
	    //printf("ff1:%2.6f %2.6f %2.6f\n",forcedft[j].x, forcedft[j].y,forcedft[j].z);
	    //if (j==3) exit(1);
    }
}

void get_pos0_from_externalfile_if_readpos(int read_pos,const char *filename_in_positions, int columns) {
    // needs: the inputs above (read_pos, filename_in_positions, columns)
    // needs: globally defined: atoms, a0, alat_lattice, N
	int j=0;
	FILE *file_in_positions;
	double a,b,c,d,e,f,dumread[columns];
	//printf("kkkkkk%s\n",filename_in_positions);
    //double *faktor2 = 4.13/8492;;   // show in absolute
    if (read_pos == 1) {
	    //file_in_positions=(FILE *)fopen("POSITIONs_1_steps_1128","rb");
	    file_in_positions=(FILE *)fopen(filename_in_positions,"rb");
        //printf("read_pos from file: yes (read_pos == %d)\n",read_pos);
        // get pos_0 from first structure
		for (j=0;j<atoms;j++) {
	        //pos[j].x=0;
	        //pos[j].y=0;
	        //pos[j].z=0;
	    if (columns==3) {
	        fscanf(file_in_positions," %lf %lf %lf\n",dumread,dumread+1,dumread+2);
        }
        else {
	        fscanf(file_in_positions," %lf %lf %lf %lf %lf %lf\n",dumread,dumread+1,dumread+2,dumread+3,dumread+4,dumread+5);
        }
            //if (j<=2) {
            //    printf("dumread %d     : %.6f %.6f %.6f\n",j,dumread[0],dumread[1],dumread[2]);
            //    }
            // dumread ist z.B. in 4x4x4 superzelle mit a = 4.14 Ang:
            //                             undisplace =
            //                             Angstrom
            // 0.01713 0.00467 16.54989  = 0 0 16.56  = 0 0 4
            // 16.53580 16.52625 4.10595 = 0 0 4.14   =
            // 16.52267 0.02688 8.30128  = 0 0 8.28   =
            // 0.02215 16.52375 12.42981 = 0 0 12.42  =
	        //pos[j].x=a0*(dumread[0]/alat_lattice/N);
	        //pos[j].y=a0*(dumread[1]/alat_lattice/N);
	        //pos[j].z=a0*(dumread[2]/alat_lattice/N);
	        d=a0*(dumread[0]/alat_lattice/N);
	        e=a0*(dumread[1]/alat_lattice/N);
	        f=a0*(dumread[2]/alat_lattice/N);
	        // debugging
            //if (j<=2) {
            //printf("pos0 of atom %d: %.6f %.6f %.6f\n",j,pos[j].x/a0,pos[j].y/a0,pos[j].z/a0);
            //}

            // deriving
            // a) get max value; if <= 1 --> *N; if <= N --> OK; if > N: --> /alat_lattice
            //p0rel=np.round([(p0*2*N+0.5)%(2*N)-0.5])/2.
            //a=((pos[j].x/a0)*2.*N+0.5); //%(2*N);
            //a=round(fmod(a,2*N)-0.5)/2.;
            a=round(fmod(((d/a0)*2.*N+0.5),2*N)-0.5)/2./N;
            b=round(fmod(((e/a0)*2.*N+0.5),2*N)-0.5)/2./N;
            c=round(fmod(((f/a0)*2.*N+0.5),2*N)-0.5)/2./N;
            if (a == 0) a = 0;
            if (b == 0) b = 0;
            if (c == 0) c = 0;
            pos0[j].x = a*a0*N;
            pos0[j].y = b*a0*N;
            pos0[j].z = c*a0*N;
            // change pos correspondingly
            pos[j].x = pos0[j].x;
            pos[j].y = pos0[j].y;
            pos[j].z = pos0[j].z;
            //printf("pos0 xxxxxxx %d: %.6f \n",j,a);
            //printf("pos0 xxxxxxx %d: %.6f \n",j,7 % 3);
            //printf("pos0 xxxxxxx %d: %.6f \n",j,(fmod(a,2*N)-0.5));
            //printf("pos0 xxxxxxx %d: %.6f \n",j,round(fmod(a,2*N)-0.5));
            //printf("pos0 xxxxxxx %d: %.6f \n",j,a);
            //printf("pos0 xxxxxxx %d: %.3f %.3f %.3f\n",j,a,b,c);
            //if (j<=2) {
	        //    printf("uuu         %d\t%.9f\t%.9f\t%.9f\n",j,pos0[j].x/a0/N,pos0[j].y/a0/N,pos0[j].z/a0/N);
            //    }
        }
	    fclose(file_in_positions);
	    //printf("hallo\n");
	    //printf("halloxxx %2.16f\n",faktor2);
        //check_pos0_for_duplicates(faktor2);
	    //exit(1);
    }
    //else {
    //    printf("NO read_pos\n");
    //}

    //for (j=0;j<atoms;j++) {
    //    printf("pos0 yyyyyyy %d: %.3f %.3f %.3f\n",j,pos0[j].x,pos0[j].y,pos0[j].z);
    //    if (j==2) {exit(1);}
    //}

}

void check_if_file_exists(const char *filename) {
    if( access( filename, F_OK ) != -1 ) {
        // file exists
      //printf ("It exists\n");
    } else {
        // file doesn't exist
      printf ("ERROR: File %s DOES NOT exist! Exit!\n",filename);
      printf ("Yout can run: OUTCAR_ene-dudl.sh > u_OUTCAR\n");
      exit(1);
    }
}

int get_lines_of_file(const char *filename) {
    check_if_file_exists(filename);
    FILE* myfile = fopen(filename, "r");
    int ch, number_of_lines = 0;

    do
    {
        ch = fgetc(myfile);
        if(ch == '\n')
            number_of_lines++;
    } while (ch != EOF);

    // last line doesn't end with a new line!
    // but there has to be a line at least before the last line
    if(ch != '\n' && number_of_lines != 0)
        number_of_lines++;

    fclose(myfile);
    printf("filename %s\n",filename);
    printf("number of lines in filename = %d\n", number_of_lines);
    return(number_of_lines-1);
}

void read_u_OUTCAR(const char *filename,double* arr_u_dft,int number_of_lines) {
//void read_u_OUTCAR(const char *filename,int number_of_lines) {
    // OUTCAR_ene-dudl.sh > u_OUTCAR

    //int number_of_lines;
    //number_of_lines = get_lines_of_file(filename);
    //printf("number of lines in %s = %d\n",filename, number_of_lines);
    double dum;
     FILE *myFile;
        myFile = fopen(filename, "r");
        //myFile = fopen(filename, "r");

        //read file into array
        //double arr_u_dft[number_of_lines];
        int i;

        for (i = 0; i < number_of_lines; i++)
        {
            //fscanf(myFile, "%lf\n", dum);
            //dudl[i].u_dft=dum;
            fscanf(myFile, "%lf\n", &arr_u_dft[i]);
            //fscanf(myFile, "%lf\n", &dudl[i].u_dft);
        }

        fclose(myFile);

        //printf("0 Number is: %2.4f\n", arr_u_dft[0]);
        //printf("1 Number is: %2.4f\n", arr_u_dft[1]);
        //printf("last Number is: %2.4f\n", arr_u_dft[number_of_lines-1]);
        //return(arr_u_dft);
}

// related to hessematrix
void get_du() {
	int j;
    //for (j=0;j<6;j++) printf("iii%-5d %-15d %-12.6f|| %-15d %-12.6f||%-15d %-12.6f\n",j,pos0[j].x,pos0[j].x*faktor,pos[j].x,pos[j].x*faktor,pos[j].x-pos0[j].x,(pos[j].x-pos0[j].x)*faktor);
	for (j=0;j<atoms;j++) {
	du[j].x = pos[j].x-pos0[j].x;
	du[j].y = pos[j].y-pos0[j].y;
	du[j].z = pos[j].z-pos0[j].z;
    }
    //for (j=0;j<6;j++) printf("-->%-5d %-15d %-12.6f|| %-15d %-12.6f||%-15d %-12.6f\n",j,pos0[j].x,pos0[j].x*faktor,pos[j].x,pos[j].x*faktor,du[j].x,du[j].x*faktor);
	//printf("\n");
}

void calculate_forces_energy_hesse(double** hessemat) {
    int j,k,k3,k3p1,k3p2,j3,j3p1,j3p2;
    set_ref_forces_energy_zero();
    get_du();
    for (j=0;j<atoms;j++) { //spalte
        j3=j*3;
        j3p1=j*3+1;
        j3p2=j*3+2;
    for (k=0;k<atoms;k++) { //zeile
        k3=k*3;
        k3p1=k*3+1;
        k3p2=k*3+2;
        //j=1;
        //printf("pos.x %9.16f\n",pos[0].x*faktor);
        //printf("du0.x %9.16f\n",du[0].x*faktor);
        //printf("h00 %9.16f\n",hessemat[0][0]);
        //printf("h00*du0.x %9.16f\n",hessemat[0][k*3+0]*du[k].x*faktor); // erste zeile == 1x
        //printf("h00*du0.y %9.16f\n",hessemat[0][k*3+1]*du[k].y*faktor);
        //printf("h00*du0.z %9.16f\n",hessemat[0][k*3+2]*du[k].z*faktor);
        //printf("h10*du0.x %9.16f\n",hessemat[1][k*3+0]*du[k].x*faktor); // zweite == 1y
        //printf("h10*du0.y %9.16f\n",hessemat[1][k*3+1]*du[k].y*faktor);
        //printf("h10*du0.z %9.16f\n",hessemat[1][k*3+2]*du[k].z*faktor);
        //printf("h20*du0.x %9.16f\n",hessemat[2][k*3+0]*du[k].x*faktor); // dritte == 1z
        //printf("h20*du0.y %9.16f\n",hessemat[2][k*3+1]*du[k].y*faktor);
        //printf("h20*du0.z %9.16f\n",hessemat[2][k*3+2]*du[k].z*faktor);

        //printf("h30*du0.x %9.16f\n",hessemat[3][k*3+0]*du[k].x*faktor); // vierte == 1x
        //printf("h30*du0.y %9.16f\n",hessemat[3][k*3+1]*du[k].y*faktor);
        //printf("h30*du0.z %9.16f\n",hessemat[3][k*3+2]*du[k].z*faktor);

        forceharm[j].x += hessemat[j3][k3]*du[k].x+hessemat[j3][k3p1]*du[k].y+hessemat[j3][k3p2]*du[k].z; // erste
        forceharm[j].y += hessemat[j3p1][k3]*du[k].x+hessemat[j3p1][k3p1]*du[k].y+hessemat[j3p1][k3p2]*du[k].z; // zweite
        forceharm[j].z += hessemat[j3p2][k3]*du[k].x+hessemat[j3p2][k3p1]*du[k].y+hessemat[j3p2][k3p2]*du[k].z;
    }
    //printf("j %dx %4.16f %4.16f %4.16f\n",j,forceharm[j].x*faktor,du[j].x*faktor,forceharm[j].x*du[j].x*faktor*faktor);
    //printf("j %dy %4.16f %4.16f %4.16f\n",j,forceharm[j].y*faktor,du[j].y*faktor,forceharm[j].y*du[j].y*faktor*faktor);
    //printf("j %dz %4.16f %4.16f %4.16f\n",j,forceharm[j].z*faktor,du[j].z*faktor,forceharm[j].z*du[j].z*faktor*faktor);


    u_harm+= (forceharm[j].x*du[j].x+forceharm[j].y*du[j].y+forceharm[j].z*du[j].z);
    }
    u_harm= u_harm*faktor_hesse2;
    if (u_harm>800) {

        printf("u_harm...... something is wrong, u_harm is too high\n");
        printf("u_harm...... %2.3f\n",u_harm);
        exit(1);
        }
}

void read_hessematrix(double** hessemat,int read_hesse,const char *filename_in_hesse,int verbose) {

    int i;
    int j;
    //printf("aa1\n");
    if (read_hesse==1) {

        /*matrix*/
        /*Use double , you have floating numbers not int*/

        //double** hessemat=malloc(3*atoms*sizeof(double*));
        for(i=0;i<3*atoms;++i)
        hessemat[i]=malloc(3*atoms*sizeof(double));


        FILE *file;
        file=fopen(filename_in_hesse, "r");
        if (verbose>1) {
	        printf("read hessematrix ...  : %s\n",filename_in_hesse);
            }
        //printf("read hessematrix ...\n");
        for(i = 0; i < 3*atoms; i++) {
              for(j = 0; j < 3*atoms; j++) {
                 //Use lf format specifier, %c is for character
                 if (!fscanf(file, "%lf", &hessemat[i][j]))
                    break;
                 hessemat[i][j]=hessemat[i][j]*-97.173617;
                // hessemat[i][j] -= '0';
                //printf("%d %d %lf\n",i,j,hessemat[i][j]); //Use lf format specifier, \n is for new line
              }
        }
        fclose(file);
        if (verbose>1) {
            printf("read hessematrix ... done\n");
        }
        //printf("%9.16f\n",hessemat[0][0]);
        //printf("%9.16f\n",hessemat[95][2]);
        //printf("%9.16f\n",hessemat[95][92]);
    }
    else {
         //printf("kk1\n");
         for(i = 0; i < 3*atoms; i++) {
            hessemat[i]=malloc(3*atoms*sizeof(double));
            //printf("kk2\n");
              for(j = 0; j < 3*atoms; j++) {
                  //printf("kk3\n");
                  hessemat[i][j] = 0;
                  //printf("kk4\n");
              }

        }
 }

    if (read_hesse==1) {
        // check read hesse last line if we have values different from 0
        //i=3*atoms-1; // last line
        //for(j = 0; j < 3*atoms; j++) {
        //        printf("%d %d %lf\n",i,j,hessemat[i][j]); //Use lf format specifier, \n is for new line
        //}
        //exit(1);

        // check if matrix is symmetric
        double check;
        double checkmax=0;
        for(i = 0; i < 3*atoms; i++) {
              for(j = 0; j < 3*atoms; j++) {
                  check=fabs(hessemat[i][j]-hessemat[j][i]);
                  if (check>checkmax) checkmax=check;
                  //printf("checkmax %d %d %2.4f %2.6f %2.6f\n",i,j,hessemat[i][j]/-97.173617,hessemat[j][i]/-97.173617,checkmax); //hessemat[i][j]-hessemat[j][i]);
                  if (checkmax>0.03) break;
                  //exit(1);

        }}
        if (checkmax>0.03) {
        printf("ERROR: hessematrix is not symmetric or mayby for other (smaller/larger?) cell;EXIT\n");
        exit(1);
        }
    }
}

//////////////////////////////
// related to la
//////////////////////////////
struct twodouble{
    double v1;
    double v2;
};

struct twodouble michaelpoly3_forces_energy(double r) {
    struct twodouble f_e;
    double rang,rrel,force,fevang,fxpoly,fypoly,fzpoly,eev1,eev2,a1,a2,a3,a4;
    ///////////////////////////
    // polynomial forces
    ///////////////////////////
    rang = r*da0_alat_lattice;
    rrel = rang/alat_lattice;
    //fevang= aa + bb*rang**(-1) + cc*rang**(-2) + dd*rang**(-3);
    // this force is positive defined for positive r (as is the poly)
    fevang= -(poly_aa + poly_bb/rrel + poly_cc/(rrel*rrel) + poly_dd/(rrel*rrel*rrel));  // seems that the - sign is necessary here to be consistent with michaels definition
    force = fevang/faktor_force;
    f_e.v1 = force;   // this is the force which is given further
    // this can be commented in to get he saparate poly forces
    //fxpoly = (x/r)*force;  // doing this saves 20% of computational time
    //fypoly = (y/r)*force;  // doing this saves 20% of computational time
    //fzpoly = (z/r)*force;  // doing this saves 20% of computational time
    //forcepoly[ind1].x+=fxpoly;	//in N
	//forcepoly[ind2].x-=fxpoly;
	//forcepoly[ind1].y+=fypoly;
	//forcepoly[ind2].y-=fypoly;
	//forcepoly[ind1].z+=fzpoly;
	//forcepoly[ind2].z-=fzpoly;
    //eev_eq = -(dd/(2.*rrel*rrel)) - cc/rrel + aa*rrel + bb*log(rrel); //- 17.1065183970574175;
    //a1 = -(0.5*poly_dd/(rrel*rrel));
    //a2 = - poly_cc/rrel;
    //a3 = poly_aa*rrel;
    //a4 = poly_bb*log(rrel);
	//printf("\n");
	//printf("----------\n");
    //printf("### aa) a1,a2,a3,a4 %2.5f %2.5f %2.5f %2.5f SUM %2.5f\n",a3,a2,a1,a4,a1+a2+a3+a4);
    //a1 = a1*faktor_energy_atom*-alat_lattice/faktor_energy_cell;
    //a2 = a2*faktor_energy_atom*-alat_lattice/faktor_energy_cell;
    //a3 = a3*faktor_energy_atom*-alat_lattice/faktor_energy_cell;
    //a4 = a4*faktor_energy_atom*-alat_lattice/faktor_energy_cell;
    //printf("### bb) a1,a2,a3,a4 %2.5f %2.5f %2.5f %2.5f SUM %2.5f\n",a3,a2,a1,a4,a1+a2+a3+a4);
    eev1 = -(0.5*poly_dd/(rrel*rrel)) - poly_cc/rrel + poly_aa*rrel + poly_bb*log(rrel); // - poly_0;
    //eev2 = eev1*-alat_lattice/faktor_energy_cell;
    eev2 = eev1*alat_lattice/faktor_energy_cell;
	//printf("--> eev1 (the real eev1 in eV/angstrom) %2.16f \n",eev1);
    // written to out_positions_forces.dat: forces: force[j].z/1.602176e-9
    // with forces[j] = fx = (x/r)*f;  # where f = eev2
    f_e.v2 = eev2;
	//printf("--> eev2 %2.16f \n",eev2);
	//printf("--> eev2*-alat_lattice! %2.16f \n",eev2*-alat_lattice);
	//if (D_par*dumm1*dumm1*faktor_energy_atom > 0.3) {exit(1);};
	//
	// can be commented in
	//u_pl+=eev2; // alle drei definitionen gleichwertig
	//
	//
		    // reassign
		    //f = force;
		    //u_la+=eev2;
		    // SOME ANALYSIS (can be commented in)
		    // SOME ANALYSIS (can be commented in)
		    // SOME ANALYSIS (can be commented in)
            //double faktor_analysis;
	        ////faktor_analysis = 4.8;
		    ////if (f*f*faktor_force*faktor_force > faktor_analysis) {
            //    //req = (alat_lattice/sqrt(2.))/da0_alat_lattice;
            //    //if (rang < 2.8991378) {
            //    if (rang < 5.3) {
            //    printf("--> rangxx         : %5.10f (in    angstrom)\n",rang);
            //    printf("--> r*da0_alat...  : %5.10f (in    angstrom)\n",r*da0_alat_lattice);
            //    printf("--> req*da0_alat...: %5.10f (in    angstrom)\n",req*da0_alat_lattice);
            //    printf("--> (r-req)*da0_...: %5.10f (r-req angstrom)\n",(r-req)*da0_alat_lattice);
            //    printf("--> rrel           : %5.10f (in rel coords)\n",rrel);
            //    double f_ev_ang;
            //    f_ev_ang = force/1.602e-9;
            //    printf("--> force/1.602e-9 : %5.10f (in eV/angsrom)\n",f_ev_ang);
            //    printf("--> force/1.602e-9 in x-dir in case of displacement along (110): %5.10f (in eV/angsrom)\n",sqrt(f_ev_ang*f_ev_ang/2));
		    //    printf("--> eev1 %2.16f \n",eev1);
		    //    printf("--> eev2 %2.16f \n",eev2);
		    //    printf("--> eev1*-alat_lattice! %2.16f \n",eev1*-alat_lattice);
		    //    printf("--> eev2*-alat_lattice! %2.16f \n",eev2*-alat_lattice);
		    //    printf("--> eev2*faktor_energy_atom %2.16f \n",eev2*faktor_energy_atom);
		    //    printf("--> a1,a2,a3,a4 %2.5f %2.5f %2.5f %2.5f \n",a1,a2,a3,a4);
		    //    printf("--> ep: %2.16f (in eV for this bond = Poly energy)\n",eev2*faktor_energy_cell);
            //    };
                //printf("--> fp: %5.10f (in eV/angstrom)\n",f*faktor_force);
            //    printf("\n");
		    //    //printf("--> ep: %2.16f (in eV for this bond = Morse energy)\n",D_par*dumm1*dumm1*faktor_energy_cell);
		    //    //
		    //    // eev1*-alat_lattice = 0.894 != eev2*faktor_energy_cell
            //}
    return f_e;
}

struct twodouble morse_forces_energy(double r) {
    double dum, dumm1,f,e;
    struct twodouble f_e;
    ///////////////////////////
    // Morse forces
    ///////////////////////////
    //
        //r=1518500250;  // == 1.5185e+09 == d0 --> r=r0 --> r*faktor:2.920351
        //r=1403923935;  // r = 2.7000
        // morse force long
        // Morse  == energ == De*(1.-np.exp(-aa*(r-re)))**2
        //                    (a − b)^2 = (a − b)(a − b) = a^2 − 2ab + b^2
        // Morse' == force == 2.*aa*De*np.exp(-aa*(r-re))*(1.-np.exp(-aa*(r-re)))
        //           --> subst: dum == np.exp(-aa*(r-re))
        // Morse' == force == 2.*De*aa*dum*(1.-dum)

    //
    //
    //printf("i1 %16.5f\n",ret.v1*faktor_force);
	//dum=exp(-a_par*(r/d0*r0_eq_morse-r0_eq_morse));  // part = e^(-a(r-re)) == dum
	//dumm1=dum-1.;
	//// this force is positive defined (for atoms that are close together)
	//f=(2.*D_par*a_par)*dumm1*dum;               // das ist die erste ableitung (also die Kraft)
    //printf("i2 %16.5f\n",f*faktor_force);

	// SOME ANALYSIS (can be commented in)
	// SOME ANALYSIS (can be commented in)
	// SOME ANALYSIS (can be commented in)
	//faktor_analysis = 4.8;
	//if (f*f*faktor_force*faktor_force > faktor_analysis) {
    //    //req = (alat_lattice/sqrt(2.))/da0_alat_lattice;
	//    printf("\n");
	//    printf("**********************************************\n");
    //    printf("--> r  : %5.10f (req   angstrom)\n",r*da0_alat_lattice);
    //    printf("--> req: %5.10f (req   angstrom)\n",req*da0_alat_lattice);
    //    printf("--> drm: %5.10f (r-req angstrom)\n",(r-req)*da0_alat_lattice);
    //    printf("-->  rm: %5.10f (in    angstrom)\n",r*da0_alat_lattice);
    //    printf("-->  fm: %5.10f (in eV/angstrom)\n",f*faktor_force);
	//    printf("-->  em: %2.16f (in eV for this bond = Morse energy)\n",D_par*dumm1*dumm1*faktor_energy_cell);
	//    //printf("--> em: %2.16f (in meV/(atom-1) for this bond = Morse energy)\n",D_par*dumm1*dumm1*faktor_energy_atom);
	//    // @ displacement of 0.6angstrom in Al, force is ~-3.6 eV/angstrom; == f*faktor_force
	//    // and energy is there 7.3 eV
	//    //
    //    }

    // In [2]: def Morseder(r,De,aa,re):
    //    ...:     return 2.*aa*De*np.exp(-aa*(r-re))*(1.-np.exp(-aa*(r-re)))
    // In [1]: def Morse(r,De,aa,re):
    //    ...:     return De*(1.-np.exp(-aa*(r-re)))**2
    //
    // In [12]: Morse(2.7,0.2793,1.432673,4.13/np.sqrt(2))
    // Out[12]: 0.038485917929653904
    // --> r:1403923935 r*faktor:2.700000 f*faktor:0.407348662 u_la:0.038485918
    //
    //
    // In [13]: Morseder(2.7,0.2793,1.432673,4.13/np.sqrt(2))
    // Out[13]: -0.40734866425804167
    // --> r:1403923935 r*faktor:2.700000 f*faktor:0.407348662 u_la:0.038485918
    //
    //
	//u_la = D_par*(1.-dum)*(1.-dum); // alle drei definitionen gleichwertig
	//dum=exp(-a_par*(r/d0*r0_eq_morse-r0_eq_morse));  // part = e^(-a(r-re)) == dum
	dum=exp(-a_par*(r*one_over_r0_eq_mor-r0_eq_morse));  // part = e^(-a(r-re)) == dum
	dumm1=dum-1.;
	// this force is positive defined (for atoms that are close together)
	f_e.v1=(2.*D_par*a_par)*dumm1*dum;               // das ist die erste ableitung (also die Kraft)

    ///////////////////////////
    // Morse energy
    ///////////////////////////
	//u_la+=D_par*dumm1*dumm1; // alle drei definitionen gleichwertig
	//u_la+=e; // alle drei definitionen gleichwertig
	//u_la = D_par*(1.-2.*dum+dum*dum); // alle drei definitionen gleichwertig
	// faktor in u_la: 1.60217598e-19
	//printf("r:%2.0f r*faktor:%2.6f f*faktor:%2.9f u_la: %2.9f\n",r,r*faktor,f*faktor_force,u_la*faktor_energy_cell);
    //print_pos_forces_to_screen(7,"",faktor,-99);
    // this is the acceleration on particle ind1 and and opposit
    // acceleration on ind2
    // dv = dt a


    f_e.v2=D_par*dumm1*dumm1; // alle drei definitionen gleichwertig
    return f_e;
}

void calculate_forces_tox(int j1,int j2, int j3, int ind1, int ind2, double x, double y, double z) {
		if ( j1 == 0 ) {  // alle Naechsten nachbarn in der yz-ebene --> dx == 0 fuer pos0
		    // x koordinate zum NN ist 0, also alle in yz ebene
            //printf("Hello, Worldxx! %.10f \n",dux/a0*ktr_par);
		    //dux=(signed)(pos0[ind1].x-pos0[ind2].x)-x;
		    //printf("%2.4f   %2.4f  cur:%2.4f dux * faktor: %2.6f\n",pos0[ind1].x*faktor,pos0[ind2].x*faktor,x*faktor,dux*faktor);
		    //printf("%2.6f\n",dux*a0ktr_par);
            //forcetmp[ind1].x +=-x*a0ktr_par;
            //forcetmp[ind2].x -=-x*a0ktr_par;
            force[ind1].x +=-x*a0ktr_par;
            force[ind2].x -=-x*a0ktr_par;
            u_la_tox+=x*faktor*x*faktor*ktr_tox;
            //printf("x %d %2.5f %2.3f %2.8f %2.16f\n",ind1,x,(x*faktor*x*faktor),x*a0ktr_par*faktor_force,a0ktr_par*faktor_force);
            //printf("x %d %2.3f %2.4f %2.7f\n",ind1,(x*faktor*x*faktor),ktr_tox,x*faktor*x*faktor*ktr_tox);
		    //@@vel[ind1].x   +=-x*dtm*a0ktr_par;	//in m/s
		    //@@vel[ind2].x   -=-x*dtm*a0ktr_par;
        }
		if ( j2 == 0 ) {  // alle NN in der xz-ebene --> dy==0 fuer delta pos0
            //printf("Hello, Worldy! \n");
		    // y koordinate zum NN ist 0, also alle in xz ebene
		    //duy=(signed)(pos0[ind1].y-pos0[ind2].y)-y;
            //forcetmp[ind1].y+=-y*a0ktr_par;
            //forcetmp[ind2].y-=-y*a0ktr_par;
            force[ind1].y+=-y*a0ktr_par;
            force[ind2].y-=-y*a0ktr_par;
            u_la_tox+=y*faktor*y*faktor*ktr_tox;
            //printf("y %2.4f\n",y);
		    //@@vel[ind2].y  -=-y*dtm*a0ktr_par;
        }
		if ( j3 == 0 ) {  // alle NN in der xy-ebene --> dz==0 fuer delta pos0
            //printf("z %2.4f\n",z);
            //printf("Hello, Worldz! \n");
		    // z koordinate zum NN ist 0, also alle in xy ebene
		    //duz=(signed)(pos0[ind1].z-pos0[ind2].z)-z;
            forcetmp[ind1].z+=-z*a0ktr_par;
            forcetmp[ind2].z-=-z*a0ktr_par;
            force[ind1].z+=-z*a0ktr_par;
            force[ind2].z-=-z*a0ktr_par;
            u_la_tox+=z*faktor*z*faktor*ktr_tox;
		    //@@vel[ind1].z  +=-z*dtm*a0ktr_par;	//in m/s
		    //@@vel[ind2].z  -=-z*dtm*a0ktr_par;
        }
        //printf("tox %2.8f\n",u_la_tox);

}

void calculate_forces_energy_la_from_simplelist(double dt,int verbose) {
    int i,ind1,ind2,j1,j2,j3;
    double x,y,z,r,f,e,fx,fy,fz;
    struct twodouble f_e;
    set_forces_energy_zero();
    //printf("\n");
    //printf("new step\n");
    for (i=0;i<atoms*20;i++) {  // Max 40 neighbors
        ind1=nn1[i].at1;
        ind2=nn1[i].at2;
        j1=nn1[i].j1;
        j2=nn1[i].j2;
        j3=nn1[i].j3;
        //printf("i %5d ind1 %5d ind2 %5d\n",i,ind1,ind2);
        if ((ind1==0) && (ind2==0)) {break;};  // this is important to not add 0 to the forces

        // abstaende nn  == Dx Dy Dz (distvec = NN1_distances[s,a1,a2])
		x=(signed)(pos[ind1].x-pos[ind2].x);
		y=(signed)(pos[ind1].y-pos[ind2].y);
		z=(signed)(pos[ind1].z-pos[ind2].z);
		r=sqrt(x*x+y*y+z*z);
		if ( r > rmax ) {rmax = r;};
		if ( r < rmin ) {rmin = r;};
        //printf("--> r!!: %5.10f --> %5.10f (in angsrom) aa: %5.10f\n",r,r/a0*alat_lattice,aa);
        //printf("--> r?1: %5.10f --> %5.10f (in angsrom) aa: %5.10f\n",r,r/a0_alat_lattice,aa);
        //printf("--> r?2: %5.10f --> %5.10f (in angsrom) aa: %5.10f\n",r,r*da0_alat_lattice,aa);

        if ( michael_poly_yes_no == 0 ) {
            f_e=morse_forces_energy(r);
            f=f_e.v1;
            e=f_e.v2;
		    u_la+=e;
        }
        else
        {   // if michael_poly_yes_no == 1
            f_e=michaelpoly3_forces_energy(r);
            f=f_e.v1;
            e=f_e.v2;
		    u_la+=e;
        }
        fx = (x/r)*f;  // doing this saves 20% of computational time
        fy = (y/r)*f;  // doing this saves 20% of computational time
        fz = (z/r)*f;  // doing this saves 20% of computational time
        //printf("--fx: %5.10f",fx/1.602176e-9);

        ////// calculate forces explicitly for writing out
        force[ind1].x+=fx;	//in N
		force[ind2].x-=fx;
		force[ind1].y+=fy;
		force[ind2].y-=fy;
		force[ind1].z+=fz;
		force[ind2].z-=fz;

        // calculate tox forces if necessary
        if (ktr_tox != 0) {
        calculate_forces_tox(j1,j2,j3,ind1,ind2,x,y,z);
        };
    }
    u_la=u_la*faktor_energy_atom;
    //u_pl=u_pl*faktor_energy_atom;
	//printf("--> EM (out1: %2.16f, SUM: %2.16f\n",u_la,u_la);
    //printf("tox FINAL %2.8f\n",u_la_tox);


    ////////////////////////////
    // comment in if tox
    ////////////////////////////
    u_la_tox=u_la_tox*faktor_per_atom;
    //printf("tox FINAL MEV %2.8f\n",u_la_tox);
    u_la=u_la+u_la_tox;
}

void calculate_forces_energy_la_NOW_USE_FROM_SIMPLELIST(double dt,int verbose) {
    // dies funktioniert zur zeit nur wenn interne coordinated von der init(T)
    // genutzt werden und nicht bei readpos generell (was es aber soll!)
	int i,j,i1,i2,i3,j1,j2,j3,steps,ind1,ind2;
	double x,y,z,r,rang,rrel,faktor_analysis,fevang,eev1,eev2,fpoly,fxpoly,fypoly,fzpoly,f,e,fx,fy,fz; //,energy;
    struct twodouble f_e;
    //dtm = dt/m_element;              // This saves 8% of computational time;
    set_forces_energy_zero();

	//printf("--> EM (in): %2.16f, SUM: %2.16f\n",u_la*faktor_energy_atom,u_la*faktor_energy_atom);
	//printf("--> EPL(in): %2.16f, SUM: %2.16f\n",u_pl*faktor_energy_atom,u_pl*faktor_energy_atom);
    for (i=0;i<atoms*first_neighbors/2;i++) {   // 12 * atoms_supercell / 2
        // atoms == number of atoms == 32 in 2x2x2 fcc supercell
        // fist_neighbors = 12
        // i = 0 .. 192
        ind1=l1nn[i].ind1;
        ind2=l1nn[i].ind2;
        //printf("i: %3d xxind1: %3d pos[ind1].{x,y,z} %3.3f %3.3f %3.3f ind2: %3d pos[ind2] %3.3f %3.3f %3.3f dxdydz %3.3f %3.3f %3.3f\n",i,ind1,pos[ind1].x*faktor,pos[ind1].y*faktor,pos[ind1].z*faktor,ind2,pos[ind2].x*faktor,pos[ind2].y*faktor,pos[ind2].z*faktor,(pos[ind2].x-pos[ind1].x)*faktor,(pos[ind2].y-pos[ind1].y)*faktor,(pos[ind2].z-pos[ind1].z)*faktor);

        j1=l1nn[i].j1;
        j2=l1nn[i].j2;
        j3=l1nn[i].j3;

        // abstaende nn  == Dx Dy Dz (distvec = NN1_distances[s,a1,a2])
		x=(signed)(pos[ind1].x-pos[ind2].x);
		y=(signed)(pos[ind1].y-pos[ind2].y);
		z=(signed)(pos[ind1].z-pos[ind2].z);
		r=sqrt(x*x+y*y+z*z);

        // Analysis
        //printf("--> r  : %5.5f (req   angstrom) xyz: %3.5f %3.5f %3.5f\n",r*da0_alat_lattice,x*da0_alat_lattice,y*da0_alat_lattice,z*da0_alat_lattice);
		if ( r > rmax ) {rmax = r;};
		if ( r < rmin ) {rmin = r;};
        //printf("--> r!!: %5.10f --> %5.10f (in angsrom) aa: %5.10f\n",r,r/a0*alat_lattice,aa);
        //printf("--> r?1: %5.10f --> %5.10f (in angsrom) aa: %5.10f\n",r,r/a0_alat_lattice,aa);
        //printf("--> r?2: %5.10f --> %5.10f (in angsrom) aa: %5.10f\n",r,r*da0_alat_lattice,aa);

		//fprintf(file_out_check_dist_r,"%.10f\n",rr);

        if ( michael_poly_yes_no == 0 ) {
            f_e=morse_forces_energy(r);
            f=f_e.v1;
            e=f_e.v2;
		    u_la+=e;
        }
        else
        {   // if michael_poly_yes_no == 1
            f_e=michaelpoly3_forces_energy(r);
            f=f_e.v1;
            e=f_e.v2;
		    u_la+=e;
        }

        fx = (x/r)*f;  // doing this saves 20% of computational time
        fy = (y/r)*f;  // doing this saves 20% of computational time
        fz = (z/r)*f;  // doing this saves 20% of computational time


		//printf("r:%3.13f f:%3.18f\n",r,f*10000000000);
		//printf("I: %d %d %d   %d %d %d\n",i1,i2,i3,j1,j2,j3);
		//printf("%d %d %d   %-2d %-2d %-2d   ind1: %-3d ind2: %-3d at1:%3.3f %3.3f %3.3f    at2:%3.3f %3.3f %3.3f --> r:%3.3f\n",l1nn[i].i1,l1nn[i].i2,l1nn[i].i3,l1nn[i].j1,l1nn[i].j2,l1nn[i].j3,ind1,ind2,pos[l1nn[i].ind1].x*faktor,pos[l1nn[i].ind1].y*faktor,pos[l1nn[i].ind1].z*faktor,pos[l1nn[i].ind2].x*faktor,pos[l1nn[i].ind2].y*faktor,pos[l1nn[i].ind2].z*faktor,r*faktor);



        ////// calculate forces explicitly for writing out
        force[ind1].x+=fx;	//in N
		force[ind2].x-=fx;
		force[ind1].y+=fy;
		force[ind2].y-=fy;
		force[ind1].z+=fz;
		force[ind2].z-=fz;

        //@@// das koennte man doch auch ganz zum schluss einfach umrechnen , oder?

        //////////////////////////////////////////////////////////////
        // TODOXXX:
        // tox forces (those are like in hessematrix and therefore noot so
        // good)
        // --> need to be redifined: for tox between atom 0 and atom 1:
        // at0 @ [0,0,0] has 11 other neighbors with certain distances, think which ones
        //   are relevant [1a,2a,3a,4a,...,11a]
        // at1 @ [0.5, 0.5, 0] has also 11 other neighbors, think which ones are relevant
        //   [1b,2b,3b,4b,...,11b]
        // ... probably only the atoms which are in the intersection of the
        // fist and second list!
        // this would be the 4 atoms on at3 == [0.5 0 0.5] at4 == [ 0.5 0 -0.5] at5 == [ 0 0.5 0.5] and at6 == [0 0.5 -0.5]
        // - wie wuerde eine parametrisierung aussehen:
        // r == dist between at0 and at1
        // for every pair, loop over [ a3, a4, a5, a6 ]
        //           atx-aty
        //  Forces[r,at0-a3 ] makes force on at0 and at3
        //  Forces[r,at0-a4 ] makes forces on atx and aty
        //  Forces[r,at0-a5 ]
        //  Forces[r,at0-a6 ]
        //  Forces[r,at1-a3 ]
        //  Forces[r,at1-a4 ]
        //  Forces[r,at1-a5 ]
        //  Forces[r,at1-a6 ]
        //  ---> remaining force on at3 has to be split between at0 and at1
        //  ---> a) get remaining force on at3
        //  ---> b) get vectors between [at0 and at3] and [at1 and at3]
        //
        // - wie wuerde eine parametrisierung aussehen: Force[r,..] =
        //
        //
        // Was man sich zuerst anschauen sollte:
        // -------------------------------------
        // 01> DIE CURVED tox einbauen, so dass die MD nicht mit tox laenger laeuft.
        //
        //             a{x}) == ANALYSE
        //              {x}) == Parametrisierungen
        //
        // aa) ANALYSE: schaue fuer alle elemente ob (nur der morse) der plot
        //     out_analyze_forces_vs_dft.dat  gerade oder gedreht? wie kriegt man
        //     das richtig hin? (kann man das aus den T=0K displacements richtig
        //     hinbekommen? sieht man da schon den shift? diesen plot mal mit
        //     den T=0K displacements machen (nur die die mit den richtigen
        //     abstaenden von der T=0K sturuktur; a.b. bei Al war das bis max
        //     0.7Ang oder so, nicht mehr, hat man nie in der MD gesehen.)
        // ab) ANALYSE: mache eine analyse: von jedem atom den average NN abstand zu allen
        //     NN (==<NN>) und dann als funktion von <NN> das delta force between la und DFT
        //
        // a) wenn quer verschiebung:
        //      - ist kraft auf atom 4.13 0.00 0.00 linear? (diese kraft ist eine
        //      korrektor fuer den morse!)
        //      - ist kraft auf atom 4.13 4.13 0.00 linear?
        //
        // b) die allerserte naeherung waere die repulsive/attraktive kraft auf
        // alle! ersten nachbarn gegeben den anderen abstaenden: also wenn das atom
        // auf 0.0 0.0 0.0 in richtung des atoms bei 2.06 2.06 0.0 bewegt wird gibt
        // es (hier eine repulsion) aller atom um das atom bei 2.06 2.06 0.0; also
        // erfaehrt z.b. das atob bei 4.13 0 0 sowie das atom bei 0 4.14 0 sowie
        // das atom auf 4.14 4.13 eine repulsive kraft!
        //
        // c) die zweite (oder erste?) naeherung waere die 2NN zu parametriseren,
        // jeoch schwaehcer also die tatsaechlichen morse, da hier die die 3body WW
        // genannt in b) wirken!
        //
        //////////////////////////////////////////////////////////////




        // WISO IST IMMER pos0[ind1].x=pos0[ind2].x?
        // WAS IST NOCHMAL j1?
        // in python:
        // ktr_tox = -0.06642349409
        // ktr_tox * 0.3 *2 == -0.039854096453999996 == force on one atom in python and c skript
        //
        // ktr_par = ktr_tox*1.32339686e-8
        // faktor_force = 1./1.602176e-9
        // N = 2
        // a0=4294967296./N
        // a0ktr_par=ktr_par/a0
        // a0ktr_par*faktor_force  == -2.5548871944591834e-10  (genau so in ipython und im c skript
        // x = 155991548
        // a0ktr_par*faktor_force*155991548.00000 == -0.0398540808429065 == force on one atom in c skript
        //
        //
        // a0ktr_par*faktor_force  == -2.5548871944591834e-10  (genau so in ipython und im c skript
        // == ktr_par/a0 * 1./1.602176e-9
        // == ktr_tox*1.32339686e-8 / 4294967296./N * 1./1.602176e-9
        // -2.5548871944591834e-10/ktr_tox = 3.846360733443891e-09
        // 3.846360733443891e-09 = 1.32339686e-8 / 4294967296./N * 1./1.602176e-9

		if ( j1 == 0 ) {  // alle Naechsten nachbarn in der yz-ebene --> dx == 0 fuer pos0
		    // x koordinate zum NN ist 0, also alle in yz ebene
            //printf("Hello, Worldxx! %.10f \n",dux/a0*ktr_par);
		    //dux=(signed)(pos0[ind1].x-pos0[ind2].x)-x;
		    //printf("%2.4f   %2.4f  cur:%2.4f dux * faktor: %2.6f\n",pos0[ind1].x*faktor,pos0[ind2].x*faktor,x*faktor,dux*faktor);
		    //printf("%2.6f\n",dux*a0ktr_par);
            forcetmp[ind1].x +=-x*a0ktr_par;
            forcetmp[ind2].x -=-x*a0ktr_par;
            force[ind1].x +=-x*a0ktr_par;
            force[ind2].x -=-x*a0ktr_par;
            u_la_tox+=x*faktor*x*faktor*ktr_tox;
            //printf("x %d %2.5f %2.3f %2.8f %2.16f\n",ind1,x,(x*faktor*x*faktor),x*a0ktr_par*faktor_force,a0ktr_par*faktor_force);
            //printf("x %d %2.3f %2.4f %2.7f\n",ind1,(x*faktor*x*faktor),ktr_tox,x*faktor*x*faktor*ktr_tox);
		    //@@vel[ind1].x   +=-x*dtm*a0ktr_par;	//in m/s
		    //@@vel[ind2].x   -=-x*dtm*a0ktr_par;
        }
		if ( j2 == 0 ) {  // alle NN in der xz-ebene --> dy==0 fuer delta pos0
            //printf("Hello, Worldy! \n");
		    // y koordinate zum NN ist 0, also alle in xz ebene
		    //duy=(signed)(pos0[ind1].y-pos0[ind2].y)-y;
            forcetmp[ind1].y+=-y*a0ktr_par;
            forcetmp[ind2].y-=-y*a0ktr_par;
            force[ind1].y+=-y*a0ktr_par;
            force[ind2].y-=-y*a0ktr_par;
            u_la_tox+=y*faktor*y*faktor*ktr_tox;
            //printf("y %2.4f\n",y);
		    //@@vel[ind2].y  -=-y*dtm*a0ktr_par;
        }
		if ( j3 == 0 ) {  // alle NN in der xy-ebene --> dz==0 fuer delta pos0
            //printf("z %2.4f\n",z);
            //printf("Hello, Worldz! \n");
		    // z koordinate zum NN ist 0, also alle in xy ebene
		    //duz=(signed)(pos0[ind1].z-pos0[ind2].z)-z;
            forcetmp[ind1].z+=-z*a0ktr_par;
            forcetmp[ind2].z-=-z*a0ktr_par;
            force[ind1].z+=-z*a0ktr_par;
            force[ind2].z-=-z*a0ktr_par;
            u_la_tox+=z*faktor*z*faktor*ktr_tox;
		    //@@vel[ind1].z  +=-z*dtm*a0ktr_par;	//in m/s
		    //@@vel[ind2].z  -=-z*dtm*a0ktr_par;
        }
        //printf("tox %2.8f\n",u_la_tox);
		//} // innere schleife ueber j
	} // aeusere schleife ueber i

    // afer u_la*faktor_energy_atom the energy (u_la) is in meV/atom
	//printf("--> EM (out0: %2.16f, SUM: %2.16f\n",u_la*faktor_energy_atom,u_la*faktor_energy_atom);
    u_la=u_la*faktor_energy_atom;
    //u_pl=u_pl*faktor_energy_atom;
	//printf("--> EM (out1: %2.16f, SUM: %2.16f\n",u_la,u_la);
    //printf("tox FINAL %2.8f\n",u_la_tox);
    u_la_tox=u_la_tox*faktor_per_atom;
    //printf("tox FINAL MEV %2.8f\n",u_la_tox);
    u_la=u_la+u_la_tox;

    //forcelax = force[j].x*faktor_force; // force la
    //forcelay = force[j].y*faktor_force; // force la
    //forcelaz = force[j].z*faktor_force; // force la
    //
    // Analysis
	// printf("--> EM (out2: %2.16f, SUM: %2.16f\n",u_la,u_la);
	// printf("--> EP (outX: %2.16f, SUM: %2.16f\n",u_pl,u_pl);
	// 0.0 out2: 0.0000000000000000, outx: 113.7627341525623450,
	// 0.1 out2: 0.2001199123119097, outx: 113.9200635744249297,
	// 0.2 out2: 0.8121808717818103, outx: 114.4057222777140339,
	// 0.3 out2: 1.8726345409233012, outx: 115.2622177763509228,
	// 0.4 out2: 3.4468617213209667, outx: 116.5655730451549772,
	// 0.5 out2: 5.6365938075227620, outx: 118.4337801760384394,
	// 0.6 out2: 8.5912804506724108, outx: 121.0400939144639949,
	// 0.7 out2: 12.52443834075111 , outx: 124.6328463556048689,
	// delta     3.94                      3.59
}

void forces_to_velocities(double dt) {
    //printf("evolve on la\n");
    double dtm;
    int j=0;
    dtm = dt/m_element;              // This saves 8% of computational time;
    for (j=0;j<atoms;j++) {
        vel[j].x+=force[j].x*dtm;
        vel[j].y+=force[j].y*dtm;
        vel[j].z+=force[j].z*dtm;
        //velff[j].x+=force[j].x*dtm;
        //velff[j].y+=force[j].y*dtm;
        //velff[j].z+=force[j].z*dtm;

        //printf("%d aaa (vel sum of all previous) %14.10f\n",j,vel[j].x);
        //printf("%d bbb (this is added to vel   ) %14.10f\n",j,force[j].x*dtm);
        //printf("%d ccc (velff                  ) %14.10f\n",j,velff[j].x);
        //fprintf(file_out_positions,"%14.10f %14.10f %14.10f %14.10f %14.10f %14.10f\n",pos[j].x*faktor,pos[j].y*faktor,pos[j].z*faktor,force[j].x/1.602176e-9,force[j].y/1.602176e-9,force[j].z/1.602176e-9);
    }
}

void forcesharm_to_velocities(double dt) {
    //printf("evolve on hesse\n");
    double dtm;
    int j=0;
    dtm = dt/m_element;              // This saves 8% of computational time;
    for (j=0;j<atoms;j++) {
        vel[j].x+=forceharm[j].x*faktor_harm_vel*dtm;
        vel[j].y+=forceharm[j].y*faktor_harm_vel*dtm;
        vel[j].z+=forceharm[j].z*faktor_harm_vel*dtm;
        //velff[j].x+=force[j].x*dtm;
        //velff[j].y+=force[j].y*dtm;
        //velff[j].z+=force[j].z*dtm;

        //printf("%d aaa (vel sum of all previous) %14.10f\n",j,vel[j].x);
        //printf("%d bbb (this is added to vel   ) %14.10f\n",j,force[j].x*dtm);
        //printf("%d ccc (velff                  ) %14.10f\n",j,velff[j].x);
        //fprintf(file_out_positions,"%14.10f %14.10f %14.10f %14.10f %14.10f %14.10f\n",pos[j].x*faktor,pos[j].y*faktor,pos[j].z*faktor,force[j].x/1.602176e-9,force[j].y/1.602176e-9,force[j].z/1.602176e-9);
    }
}

void calculate_forces_energy_la_old(double dt,int verbose) {
    // dies funktioniert zur zeit nur wenn interne coordinated von der init(T)
    // genutzt werden und nicht bei readpos generell (was es aber soll!)
	int j,i1,i2,i3,j1,j2,j3,steps,ind1,ind2;
	double x,y,z,r,dum,forcex,forcey,forcez,dux,duy,duz,dtm;
    dtm = dt/m_element;              // This saves 8% of computational time;
    set_forces_energy_zero();
	for (i3=0;i3<N*2;i3++) for (i2=0;i2<N*2;i2++) for (i1=(i3+i2)%2;i1<N*2;i1+=2) {
	    // i{0,1,2,3} sind die absoluten koordinaten dat atome von 0 bis 2*N-1 (0-3 in 2x2x2sc)
	    // i{0,1,2,3} sind einfach die koordinaten
	    // {0,0,0} --> xyz = {0,0,0}
	    // {1,1,0} --> xyz = {2.065,2.065,0}
	    // {2,0,0} --> xyz = {4.13,0,0}
	    // {3,1,0} --> xyz = {6.195,2.065,0}
	    // --> i3 ist die z ebene --> die atomstruktur wird z-ebene fuer-z ebene aufgebaut
	    //                            also erst alle z=0, dann alle z=2.065 usw.
	    //                            i3 geht von 0 bis 2*N-1 (0-2 in 2x2x2sc)
	    // --> i2 ist die y koordinate --> i2 geht von 0 bis 2*N-1 (0-2 in 2x2x2sc)
	    // --> i1 ist die x koordinate und geht im prinzip auch ueber alle indizes
	    // von 0 bis 2*N-1 (0-3 in 2x2x2sc) wobei jedoch nie ein atom bei {1,0,0}
	    // sitzt: {0,0,0} ist das aufatom und {2,0,0} ist der zweit naechste
	    // nachbar und {1,1,0} ist der naechste nachbar; von daher das modulo
	    // (i3+i2)%2 das dafuer sorgt dass man die richtige fcc struktur hat.
		ind1=((i3/2)*2*N+i2)*2*N+i1;
		j1=-1;  // das sind die relativen coordinaten der naechsten nachbarn
		j2=1;   // das sind die relativen coordinaten der naechsten nachbarn
		j3=0;   // das sind die relativen coordinaten der naechsten nachbarn
		goto JUMPIN;
		for (j3=-1;j3<=1;j3++) for (j2=-1;j2<=1;j2++) for (j1=-1+(j3+j2+3)%2;j1<=1;j1+=2) {
		    // j{1,2,3} sind die relativen coordinaten der naechsten nachbarn
		    // z.b. {-1,1,0},{1,1,0},{0,-1,1},{-1,0,1},{1,0,1},{0,1,1}
		    // j{2,3} gehen ueber {-1,0,1} (also die moeglichen positionen der
		    // Naechsten nachbarn
		    // j1 im prinzip auch
JUMPIN:
			ind2=(((i3+j3+2*N)%(2*N)/2)*2*N+(i2+j2+2*N)%(2*N))*2*N+(i1+j1+2*N)%(2*N);

            // wenn man ab hier kiene i{0,1,2,3} oder j{0,1,2,3} braeucht, dann
            // waere list [ind1,ind2] genug, alternative liste:
            // [i1,i2,i3,j1,j2,j3,ind1,ind2]

            // abstaende nn
			x=(signed)(pos[ind1].x-pos[ind2].x);
			y=(signed)(pos[ind1].y-pos[ind2].y);
			z=(signed)(pos[ind1].z-pos[ind2].z);
			r=sqrt(x*x+y*y+z*z);

            // morse force long
			dum=exp(-a_par*(r/d0*r0_eq_morse-r0_eq_morse));  // dum = e^(-a(r-re))
			dum=(2.*D_par*a_par)*(dum-1.)*dum;     // dass soll der andere teil vom morse sein, check it!
            // this is the acceleration on particle ind1 and and opposit
            // acceleration on ind2
            // dv = dt a
            forcex = (x/r)*dum;  // doing this saves 20% of computational time
            forcey = (y/r)*dum;  // doing this saves 20% of computational time
            forcez = (z/r)*dum;  // doing this saves 20% of computational time


			//printf("r:%3.13f dum:%3.18f\n",r,dum*10000000000);
			//printf("I: %d %d %d   %d %d %d\n",i1,i2,i3,j1,j2,j3);
			printf("%d %d %d   %-2d %-2d %-2d   ind1: %-3d ind2: %-3d at1:%3.3f %3.3f %3.3f    at2:%3.3f %3.3f %3.3f --> r:%3.3f\n",i1,i2,i3,j1,j2,j3,ind1,ind2,pos[ind1].x*faktor,pos[ind1].y*faktor,pos[ind1].z*faktor,pos[ind2].x*faktor,pos[ind2].y*faktor,pos[ind2].z*faktor,r*faktor);



            ////// calculate forces explicitly for writing out
            /// force long  ... tox above
            force[ind1].x+=forcex;	//in N
		    force[ind2].x-=forcex;
		    force[ind1].y+=forcey;
		    force[ind2].y-=forcey;
		    force[ind1].z+=forcez;
		    force[ind2].z-=forcez;

			vel[ind1].x+=dtm*forcex;	//in m/s
			vel[ind2].x-=dtm*forcex;
			vel[ind1].y+=dtm*forcey;
			vel[ind2].y-=dtm*forcey;
			vel[ind1].z+=dtm*forcez;
			vel[ind2].z-=dtm*forcez;

            //////////////////////////////////////////////////////////////
            // tox forces (those are like in hessematrix and therefore noot so
            // good)
            // --> need to be redifined: for tox between atom 0 and atom 1:
            // 0 has 11 other neighbors with certain distances, think which ones
            //   are relevant [1a,2a,3a,4a,...,11a]
            // 1 has also 11 other neighbors, think which ones are relevant
            //   [1b,2b,3b,4b,...,11b]
            // ... probably only the atoms which are in the intersection of the
            // fist and second list!
            //////////////////////////////////////////////////////////////
		    if ( j1 == 0 ) {
		        // x koordinate zum NN ist 0, also alle in yz ebene
                //printf("Hello, Worldxx! %.10f \n",dux/a0*ktr_par);
		        dux=(signed)(pos0[ind1].x-pos0[ind2].x)-x;
                force[ind1].x +=dux*a0ktr_par;
                force[ind2].x -=dux*a0ktr_par;
			    vel[ind1].x   +=dtm*dux*a0ktr_par;	//in m/s
			    vel[ind2].x   -=dtm*dux*a0ktr_par;
            }
		    if ( j2 == 0 ) {
                //printf("Hello, Worldy! \n");
		        // y koordinate zum NN ist 0, also alle in xz ebene
		        duy=(signed)(pos0[ind1].y-pos0[ind2].y)-y;
                force[ind1].y+=duy*a0ktr_par;
                force[ind2].y-=duy*a0ktr_par;
			    vel[ind1].y+=dtm*duy*a0ktr_par;	//in m/s
			    vel[ind2].y-=dtm*duy*a0ktr_par;
            }
		    if ( j3 == 0 ) {
                //printf("Hello, Worldz! \n");
		        // z koordinate zum NN ist 0, also alle in xy ebene
		        duz=(signed)(pos0[ind1].z-pos0[ind2].z)-z;
                force[ind1].z+=duz*a0ktr_par;
                force[ind2].z-=duz*a0ktr_par;
			    vel[ind1].z+=dtm*duz*a0ktr_par;	//in m/s
			    vel[ind2].z-=dtm*duz*a0ktr_par;
            }
		} // innere schleife ueber j
	} // aeusere schleife ueber i
}

// mathematical functions
//double sum_array(double a[], int num_elements) {
//       int i;
//       double sum=0;
//       for (i=0; i<num_elements; i++) {
//            sum += a[i];
//
//       //printf("kkk %5.4f %5.4f\n",a[i],sum);
//       }
//       return(sum);
//}

//double get_average(double average,int step ,double val) {

//struct Bar{
//    double sum;
//    double mean;
//    double variance;
//    double std;
//};
//
//
//
//
//struct Bar funct();
//struct Bar funct(double a[], int num_elements){
//    struct Bar result;
//    //printf("kkk %5.4f %5.4f\n",a[i],sum);
//    int i;
//    result.sum=0;
//    result.mean=0;
//    result.variance=0;
//    result.std=0;
//    //result.sum
//    for (i=0; i<num_elements; i++) {
//        result.sum+= a[i];
//    }
//
//    //result.mean
//    result.mean = result.sum/num_elements;
//
//    //result.variance
//    for (i=0; i<num_elements; i++) {
//        //printf("%2.4f %2.4f %2.4f\n",a[i],result.mean,pow(a[i]-result.mean,2));
//        result.variance+= pow(a[i]-result.mean,2);
//        //printf("%2.4f %2.4f %2.4f\n",a[i],result.mean,result.variance);
//    }
//    result.variance = result.variance/num_elements;
//
//    //result.std
//    result.std = sqrt(result.variance);
//    return result;
//}

// for indexing the atoms according to nn list
int check_l1nn_for_pos0(int exitonerr, int verbose) {
    int i,j,ind1,ind2,a,b,c,needs_resorting,c1,c2,c3,c4,c5,c6,c123,c456;
	double x,y,z,r,r0,faktrr,rtest,r0test;
	faktrr=faktor/alat_lattice*2;
    i=0;
    c123=0; // 0 ist falsch (False)
    c456=0; // 0 ist falsch (False)
    needs_resorting=0;
    //printf("needs_resorting (0): %d\n",needs_resorting);
    r0=0;
    for (i=0;i<atoms*first_neighbors/2;i++) {
        ind1=l1nn[i].ind1;
        ind2=l1nn[i].ind2;
		x=(signed)(pos0[ind1].x-pos0[ind2].x);
		y=(signed)(pos0[ind1].y-pos0[ind2].y);
		z=(signed)(pos0[ind1].z-pos0[ind2].z);
		r=sqrt(x*x+y*y+z*z);
		rtest=fabs(x)+fabs(y)+fabs(z);
        if (i==0) {r0=r;r0test=rtest;};
        c1=pos0[l1nn[i].ind1].x*faktrr-l1nn[i].i1;
        c2=pos0[l1nn[i].ind1].y*faktrr-l1nn[i].i2;
        c3=pos0[l1nn[i].ind1].z*faktrr-l1nn[i].i3;
        c4=pos0[l1nn[i].ind2].x*faktrr-l1nn[i].k1;
        c5=pos0[l1nn[i].ind2].y*faktrr-l1nn[i].k2;
        c6=pos0[l1nn[i].ind2].z*faktrr-l1nn[i].k3;
        if (c1==c2 && c2==0) {c123=1;} else {c123=0;}
        if (c4==c5 && c5==0) {c456=1;} else {c456=0;}



        //if (verbose>=1) {
        //    print_l1nn_to_screen(i,i,"in check_l11_for_pos0");
        //    exit(1);};
		//if ((r!=r0) || (c123!=1) || (c456!=1)) {
		if ((fabs(rtest-r0test)>1) || (c123!=1) || (c456!=1)) {
		    needs_resorting=1;
            //printf("needs_resorting here: %d\n",needs_resorting);
		    if ((verbose==1) || (exitonerr==1))
		    {

		    //if (r!=r0) {
		    if (fabs(rtest-r0test)>1) {
		        printf("ERROR: i: %d The l1nn list has not the correct 1nn r!=r0  r: %4.9f   r0:%4.9f\n",i,r*faktor,r0*faktor);
		        printf("ERROR: i: %d The l1nn list has not the correct 1nn r!=r0  r: %2.0f   r0:%2.0f diff:%2.0f\n",i,rtest,r0test,fabs(rtest-r0test));
		        ind1=0;ind2=11;
		        printf("ind:%d r:%d %d %d\n",ind1,pos0[ind1].x,pos0[ind1].y,pos0[ind1].z);
		        printf("ind:%d r:%d %d %d\n",ind2,pos0[ind2].x,pos0[ind2].y,pos0[ind2].z);
		        ind1=0;ind2=7;
		        printf("ind:%d r:%d %d %d\n",ind1,pos0[ind1].x,pos0[ind1].y,pos0[ind1].z);
		        printf("ind:%d r:%d %d %d\n",ind2,pos0[ind2].x,pos0[ind2].y,pos0[ind2].z);
		        //printf("%d\n",pos0[7].x+pos0[11].x);
		        printf("xyz: %2.0f %2.0f %2.0f\n",x,y,z);
            }

		    if ((c123!=1) || (c456!=1)) {
		    printf("ERROR: i: %d   The (x x x, x x x) are not all 0! Exit!\n",i);
            print_l1nn_to_screen(0,i+1,"");
            }
		    }


		    if (exitonerr==1) {exit(1);}
            else {break;};
        }




    }
    return(needs_resorting);
}

int get_index_atom_from_pos0(unsigned int x,unsigned int y, unsigned int z) {
    int j;
    //printf("xyz:%d %d %d   %3.4f %3.4f %3.4f\n",x,y,z,x*faktor,y*faktor,z*faktor);
    for (j=0;j<atoms;j++) {
        if ((pos0[j].x == x) && (pos0[j].y == y) && (pos0[j].z == z)) return j;
    }
    printf("ERROR: index not found!\n");
    exit(1);
}

void resorting_l1nn(int verbose) {
    int i; //,idx,idx1,idx2;
    for (i=0;i<atoms*first_neighbors/2;i++) {
        l1nn[i].ind1 = get_index_atom_from_pos0(l1nn[i].x0,l1nn[i].y0,l1nn[i].z0);
        l1nn[i].ind2 = get_index_atom_from_pos0(l1nn[i].x1,l1nn[i].y1,l1nn[i].z1);
        //printf("1indold: %d indnew: %d\n",l1nn[i].ind1,idx1);
        //printf("2indold: %d indnew: %d\n",l1nn[i].ind2,idx2);
        //exit(1);
    }
}

void check_l1nn_for_pos0_and_resort_1nn() {
    int needs_resorting;
    // fist check
    needs_resorting = check_l1nn_for_pos0(0,0);
    printf("l1nn needs_resorting (a) (1==True): %d\n",needs_resorting);
    //exit(1);
    // resort if necessary and recheck
    if (needs_resorting==1) {
        //print_l1nn_to_screen(0,12,"");
        resorting_l1nn(0);
        //check indices of l1nn list
        //print_l1nn_to_screen(0,12,"");
        printf("oiu\n");
        needs_resorting = check_l1nn_for_pos0(1,1);
        printf("l1nn needs_resorting (b) (1==True): %d\n",needs_resorting);
        if (needs_resorting==1) {
        printf("needs_resorting did not work out: %d\n",needs_resorting);
        exit(1);
        }
    }
}

void get_temperature_from_velocities() {
    int ind1;
    double dum;
	dum=0;
	for (ind1=0;ind1<atoms;ind1++) {
		dum+=vel[ind1].x*vel[ind1].x+vel[ind1].y*vel[ind1].y+vel[ind1].z*vel[ind1].z;
		//dum+=velff[ind1].x*velff[ind1].x+velff[ind1].y*velff[ind1].y+velff[ind1].z*velff[ind1].z;
	}
	//temperatue = dum*(m_element/k_B/(12.*N*N*N));
	temperature = dum*faktor_temperature;

    //temperature_sum+=temperature;
    //temperature_av=temperature_sum/temperature_index;
    temperature_av = getaverage(temperature_av,temperature_index,temperature);
    temperature_index+=1;
}

void write_temperature(FILE *file_out_temp,int step) {
        //printf("----> j: %4.10f\n",temperature);
		fprintf(file_out_temp,"%.2f\n",temperature);
		if (temperature > 2000.) {
		    printf("Temperature exceeds 2000 in step %d something in the MD went wrong?\n",step);
		    exit(1);
        };
}
void write_temperature_average(FILE *file_out_temp) {
		fprintf(file_out_temp,"%.2f\n",temperature_av);
}

// related to md
// most easy: standard verlet algorithm
// better: velocity verlet
// best: SHAKE method
void verletstep_from_velocities(double dt,FILE *file_out_temp,double faktor_verlet) {
        int ind1,ind2;
        double dum; //,tempout; //j,ka,kb,kc;
		dum=0;

        //double dtm;
        //int j=0;
        //dtm = dt/m_element;              // This saves 8% of computational time;
        //for (j=0;j<atoms;j++) {
        //    velff[j].x+=force[j].x*dtm;
        //    velff[j].y+=force[j].y*dtm;
        //    velff[j].z+=force[j].z*dtm;
        //    //printf("%d aaa (vel sum of all previous) %14.10f\n",j,vel[j].x);
        //    //printf("%d bbb (this is added to vel   ) %14.10f\n",j,force[j].x*dtm);
        //    //printf("%d ccc (velff                  ) %14.10f\n",j,velff[j].x);
        //}

        //print_vel_velff_to_screen(4,"======> VEL VELFF ",0);

		for (ind1=0;ind1<atoms;ind1++) {
			pos[ind1].x+=(unsigned)(signed int)round(faktor_verlet*vel[ind1].x);
			pos[ind1].y+=(unsigned)(signed int)round(faktor_verlet*vel[ind1].y);
			pos[ind1].z+=(unsigned)(signed int)round(faktor_verlet*vel[ind1].z);
			//dum+=vel[ind1].x*vel[ind1].x+vel[ind1].y*vel[ind1].y+vel[ind1].z*vel[ind1].z;
		}


		//tempout = dum*(m_element/k_B/(12.*N*N*N));
        //get_temperature_from_velocities();
		//printf("temp1 %.6f %.6f\n",tempout,temperature);		// prints out temperatur to screen
		//fprintf(file_out_temp,"%.6f %.6f\n",tempout,temperature);		// prints out temperatur to screen

	    //if ( tempout > 1900 ) {
	    //printf("The temperature reached > 2.000 Kelvin ...\n");
	    //printf("... this probably means that the structure exploded. I exit here.\n");
	    //exit(1);
        //}
		//printf("inposold2:%.10f\n",posold[0].x/a0);
}

void preequilibration_if_no_readpos(double dt, int read_pos,FILE *file_out_temp,double** hessemat,double faktor_verlet,int read_hesse, int evolve_md_on_hesse) {
    int i;
    printf("dt      :%2.16f\n",dt);
    //double dt_preeq = 0.0000000000000010;
    //if (dt < dt_preeq) {dt_preeq=dt;};  # somehow does not work

    double dt_preeq = dt;
	clock_t start,end;
    printf("dt_preeq:%2.16f\n",dt_preeq);
    if (read_pos!=1) {
        int preeq=3000000/atoms;
	    start = clock();
        printf("... preequilibration start (depending on cellsize...) using %d steps on timestep of 1fs or below: %2.16f\n",preeq,dt_preeq);
        for (i=0;i<preeq;i++) {

            calculate_forces_energy_la_from_simplelist(dt_preeq,0);
            if (read_hesse==1) {calculate_forces_energy_hesse(hessemat);}
            if (evolve_md_on_hesse==1) {forcesharm_to_velocities(dt_preeq);}
            else {forces_to_velocities(dt_preeq);}

            //calculate_forces_energy_hesse(hessemat);   // 1.6-1.7 sec
            //if (i%l==l-1) {write_positions_forces(file_out_positions,faktor);}; // 4.3 sec
            verletstep_from_velocities(dt_preeq,file_out_temp,faktor_verlet);
        }
        printf("... preequilibration finished\n");
	    end = clock();
        write_out_info_add_timing(start,end,"preeq");
        }
}




int main(int argc,char *argv[]){
    printf("--------------------------------------  from here in c code ----------------------------------------\n");
    printf("--------------------------------------  from here in c code ----------------------------------------\n");
    printf("--------------------------------------  from here in c code ----------------------------------------\n");
    int i=0;        // i zaehlt die inputvariablen
    int j=0;
    int k=0;
    int ii=0;
    int zeitschritte,l;
    clock_t start,end;

    //FILE *file_out_check_dist_fcheck;
    double T,dt,a,b,c,rr;
    double dudl_sum=0;
    double dumread[3];
    double** hessemat=malloc(3*atoms*sizeof(double*));
    int read_pos=0;         // 0: dont read POSITIONs; 1: read POSITIONs
    int read_pos0=0;        // 0: dont read POSITIONs; 1: read POSITIONs
    int read_uoutcar=0;     // 0: dont read u_OUTCAR; 1: read u_OUTCAR
    int read_forces=0;      // 0: dont read forces; 1: read forces
    int verbose=0;          // 0: not verbose        ; 1: verbose
    int write_analyze=0;    // 0: dont write analyse ; 1: write the analyse
    int write_analyze_forces=0;    // 0: dont write analyse ; 1: write the analyse
    int write_positions_forces=0;    // 0: dont write analyse ; 1: write the analyse
    int write_positions_rel=0;
    int write_positions=0;    // 0: dont write analyse ; 1: write the analyse
    int read_hesse=0;       // 0: dont read the hesse
    int evolve_md_on_hesse=0;  // perform MD on hessematrix (lambda=0.0)
    const char *filename_in_positions = NULL;   //"POSITIONs";
    const char *filename_in_positions0 = NULL;   //"POSITIONs";
    const char *filename_in_hesse = NULL;       //"HesseMatrix_sphin";
    const char *filename_in_uoutcar = NULL;     //"u_OUTCAR";
    int columns0=3;  // columns when reading in pos0 for hessematrix separately

    //printf("columns xx    %d\n",columns);
    printf("argc    %d\n",argc);
    printf("i (in)  %d\n",i);
    printf("T (in)  %lf\n",T);
    //for (i=0;i<(argc);i++) {
    //    if (i==0) {printf("(a) Step %d   \n",i);}
    //}
    if (argc > 1)  {sscanf(argv[1],"%lf",&T);};
    if (argc > 2)  {sscanf(argv[2],"%lf",&dt);};
    if (argc > 3)  {sscanf(argv[3],"%d",&zeitschritte);};
    if (argc > 4)  {sscanf(argv[4],"%d",&l);};
    if (argc > 5)  {sscanf(argv[5],"%d",&read_pos);};
    if (argc > 6)  {sscanf(argv[6],"%d",&verbose);};
    if (argc > 7)  {sscanf(argv[7],"%d",&write_analyze);};
    if (argc > 8)  {sscanf(argv[8],"%d",&read_hesse);};
    if (argc > 9)  {sscanf(argv[9],"%d",&read_forces);};
    if (argc > 10) {sscanf(argv[10],"%d",&read_uoutcar);};
    if (argc > 11) {sscanf(argv[11],"%d",&write_positions_rel);};
	//if (argc==11) {
	//	i+=sscanf(argv[1],"%lf",&T);
	//	i+=sscanf(argv[2],"%lf",&dt);
	//	i+=sscanf(argv[3],"%d",&zeitschritte);
	//	i+=sscanf(argv[4],"%d",&l);
	//	i+=sscanf(argv[5],"%d",&read_pos);
	//	i+=sscanf(argv[6],"%d",&verbose);
	//	i+=sscanf(argv[7],"%d",&write_analyze);
	//	i+=sscanf(argv[8],"%d",&read_hesse);
	//	i+=sscanf(argv[9],"%d",&read_forces);
	//	i+=sscanf(argv[10],"%d",&read_uoutcar);
	//}
    printf("i (out) %d\n",i);
    printf("T (out) %lf\n",T);


    check_arguments(argc,argv);

    if (verbose>2) {printf("POSxx 0  \n");};
    int columns=0;          // 3 columns in POSITIONs when only positions; 6 when positions and forces;
	if (read_pos==1) {columns=3;};  // after read in of arguments
	if (read_forces==1) {columns=6;};  // after read in of arguments
    if ((columns !=0) && (columns!=6) && (columns!=3)) {printf("neither 0 nor 3 nor 6! exit\n");exit(1);}

    write_out_info(T,dt,zeitschritte,l,read_pos,read_pos0,read_uoutcar,verbose,write_analyze,filename_in_positions,filename_in_positions0,filename_in_hesse,read_hesse, read_forces,columns,columns0,filename_in_uoutcar,evolve_md_on_hesse,write_positions_forces,write_positions,write_positions_rel);

	dt*=1e-12;	                    //in s
    double faktor_verlet=(d0*dt/r0_eq_morse);
    // ok : \rm a.out;gcc md_long_tox.c;./a.out 1787 .001 940000 1 0 0 0 1   (ca. 7.6MB memory)
    // SEG: \rm a.out;gcc md_long_tox.c;./a.out 1787 .001 9940000 1 0 0 0 1

    //if (read_hesse ==0) {filename_in_hesse=NULL;};
	//if (read_pos   ==0) {filename_in_positions=NULL;};
    //printf("read_forces1 %d\n",read_forces);
    //printf("columns xx    %d\n",columns);

    //// arr_u_{dft,la,harm}
    double* arr_u_dft      =malloc(stepsdudlmax*sizeof(double*));


    if (verbose>2) {printf("POSxx 3  \n");};
    if (read_uoutcar==1) {read_u_OUTCAR(filename_in_uoutcar,arr_u_dft,zeitschritte);};
    if (verbose>2) {printf("POSxx 10, %d zeitschritte \n",zeitschritte);};

    ////////////////////////////////////////////////////////////
    // checks
    ////////////////////////////////////////////////////////////
    if (verbose>2) {printf("POS 1, init(T)  \n");};
    if (zeitschritte>stepsmax) {
        printf("zeitschritte (too high for dudlarray): %d stepsmax: %d\n",zeitschritte,stepsmax);
        exit(1);
    }


	init(T);                        // initialize pos0
    //print_pos_0_to_screen(12, "POS_0 nach init");

    //exit(1);

    if (verbose>2) {printf("POS 3  \n");};

    if (write_positions==1) {
	    file_out_positions          =(FILE *)fopen("out_positions.dat","w");};
    if (write_positions_forces==1) {
	    file_out_positions          =(FILE *)fopen("out_positions_forces.dat","w");};


    print_l1nn_to_screen(0,13,"A) l1nn  nach init(T)");
    //write_out_pos0_to_out_EqCoords_direct(); // only for checking forces to python scripts
    //write_out_cell();       // only for checking forces to python scripts
    //exit(1);

    if (write_analyze==1) {
	    file_out_check_dist_xyz =(FILE *)fopen("out_analyze_positions_distances_from_equilibrium_xyz.dat","w");
	    file_out_check_dist_r   =(FILE *)fopen("out_analyze_positions_distances_from_equilibrium_r.dat","w");
	    file_out_check_dist_nn  =(FILE *)fopen("out_analyze_positions_distances_nn.dat","w");
	    file_out_check_dist_nn_proj  =(FILE *)fopen("out_analyze_positions_distances_nn_proj.dat","w");
	    file_out_prl15_2a       =(FILE *)fopen("out_analyze_positions_prl15_2a.dat","w");
	    file_out_prl15_2au      =(FILE *)fopen("out_analyze_positions_prl15_2au.dat","w");
	    file_out_new1            =(FILE *)fopen("out_analyze_new1.dat","w");
	    file_out_new2            =(FILE *)fopen("out_analyze_new2.dat","w");
	}


	start = clock();
    if (read_pos==1) {
        printf("///////////////////////////////////////////////////////////////////////////////\n");
        printf("/////         initialize MD from POSITIONs %d steps\n",zeitschritte);
        printf("///////////////////////////////////////////////////////////////////////////////\n");
        file_in_positions           =(FILE *)fopen(filename_in_positions,"rb");
	    file_out_dudl               =(FILE *)fopen("out_dudl.dat","w");
	    file_out_dudl_av            =(FILE *)fopen("out_dudlav.dat","w");
        write_dudl_head(file_out_dudl,file_out_dudl_av);
	    if (write_analyze_forces==1) {
	        file_forces              =(FILE *)fopen("out_forcesdiff.dat","w");
	        file_forces_av           =(FILE *)fopen("out_forcesstd.dat","w");
	        file_forces_vs_forces_dft=(FILE *)fopen("out_forces_vs_dft_prl2015_Fig3b.dat","w");
            write_forces_head(file_forces,file_forces_av);
            write_analyze_forces_vs_dft_head(file_forces_vs_forces_dft);
        };




        get_pos0_from_externalfile_if_readpos(read_pos,filename_in_positions,columns);
        print_pos_0_to_screen(12, "POS_0 from externalfile");
        get_1NN(310, "POS_0 from get_1NN");
        //exit(1);
        check_pos0_for_duplicates(faktor,filename_in_positions);
        check_l1nn_for_pos0_and_resort_1nn();
        check_pos0_for_duplicates(faktor,filename_in_positions);
        check_pos_vs_pos0_if_same_positions_and_check_if_correct_alat();
        if (read_hesse!=0) {
            check_if_file_exists(filename_in_hesse);
            read_hessematrix(hessemat,read_hesse,filename_in_hesse,verbose);
            }

        //print_pos_0_to_screen(12, "POS_0");
        //print_pos_to_screen(12, "POS",i);

        ///////////////////////////////////////////////////////////////////////////////
        // TO GET TEMPERATURE AND VELOCITIES IS NOT NECESSARY WHEN RUNNING ON POSITIONS
        ///////////////////////////////////////////////////////////////////////////////
	    //file_out_temp               =(FILE *)fopen("out_temperature.dat","w");
	    //file_out_temp_av            =(FILE *)fopen("out_temperatureav.dat","w");


        if (verbose>2) {printf("STARTING MD  \n");};
        printf("STARTING MD  \n");
        //int l = 100;
        for (i=0;i<(zeitschritte);i++) {
            //double progress = (i / zeitschritte) * 100;
            // this strats from 0 to have right numbering with dudl;
            // is i=0 also necessary to get the correct first positioins? -> no
            // for the new way we average it might be also necessary to start with 0
            if (verbose>1) {printf("(a) Step %d   \b",i);}
            //if (i%l==l-1) {if (verbose>0) {printf("%d\b%",i);}};
            u_dft                  = arr_u_dft[i];
            get_pos_from_externalfile(file_in_positions,columns);
            //calculate_forces_energy_la(dt,verbose);
            // HIER !!!!!!!!
            calculate_forces_energy_la_from_simplelist(dt,verbose);
            //printf("hierxxz\n");
            //exit(1);
            if (read_hesse==1) {calculate_forces_energy_hesse(hessemat);};

            // analyze_forces is only of interest if DFT forces available
            analyze_forces(i,file_forces,file_forces_av,  file_forces_vs_forces_dft,write_analyze_forces,read_hesse);

            if (write_analyze==1) {
                write_analyze_positions(i);}; //,file_out_check_dist_xyz,file_out_check_dist_r,file_out_check_dist_nn,file_out_check_dist_nn_proj,file_out_prl15_2a,file_out_prl15_2au);};

            //forces_to_velocities(dt); // only necessary when volocities are necessary


            if (verbose>=2) {
                print_pos_forcesdft_to_screen(4,"======> DFT POS FORCES ",faktor,i);
                print_pos_forcestmp_to_screen(12,"======> TMP POS FORCES ",faktor,i);
                print_forces_to_screen(4,"======> LA FORCES ",faktor,i);
                print_forcesharm_to_screen(4,"======> Harmonic FORCES ",faktor,i);
                //print_velocities_to_screen(4,"======> VELOCITIES ",faktor,i);
                // pos          : pos[j].x*faktor
                //
                // forcedft     : forcedft[j].x
                // forcela      : force[j].x*faktor_force
                //              : force[j].x/1.602176e-9
                // forceharm    : forceharm[j].x*faktor
                //
                // vel dft      : DONT have those currently but could be checked from forces
                // vel la       : forcela*dtm  (dtm = dt/m_element;)
                // vel harm     : DONT have those currently
                //
                // temperature  : from velocities
                // for verlet   : from velocities
                //
                // u_dft
                // u_la
                // u_harm
                }
            if (verbose>=2) {
                printf("Step %d   \n",i);
                printf("           DFT    %10.3f\n",arr_u_dft[i]);
                printf("           LA LON %10.3f\n",u_la-u_la_tox);
                printf("           LA TOX %10.3f\n",u_la_tox);
                printf("           LA SUM %10.3f (long+tox)\n",u_la);
                printf("           HARM   %10.3f\n",u_harm);
                printf("\n");
                printf("\n");
                }

            if (write_positions_forces==1) {
                if (i%l==l-1) {write_positions_forces_tofile(file_out_positions);};}
            if (write_positions==1) {
                if (i%l==l-1) {write_positions_tofile(file_out_positions,write_positions_rel);};}

            ///////////////////////////////////////////////////////////////////////////////
            // TO GET TEMPERATURE AND VELOCITIES IS NOT NECESSARY WHEN RUNNING ON POSITIONS
            ///////////////////////////////////////////////////////////////////////////////
            //get_temperature_from_velocities();
            //write_temperature(file_out_temp,i);  // write it for every step?
            //write_temperature_average(file_out_temp_av);  // write it for every step?
            //printf("%d %10.3f\n",i,temperature);
            write_dudl(i,file_out_dudl,file_out_dudl_av,read_hesse,read_uoutcar);
            }
        printf("\n");
        printf("///////////////////////////////////////////////////////////////////////////////\n");
        printf("/////      MD from POSITIONs DONE last step was %d \n",i);
        if (michael_poly_yes_no == 1) {
            printf("/////      MODEL: POLYNOM MICHAEL\n");
            printf("/////      MODEL: poly_aa %3.8f \n",poly_aa);
            printf("/////      MODEL: poly_bb %3.8f \n",poly_bb);
            printf("/////      MODEL: poly_cc %3.8f \n",poly_cc);
            printf("/////      MODEL: poly_dd %3.8f \n",poly_dd);
            printf("/////      MODEL: poly_ee %3.8f \n",poly_ee);
            printf("/////      MODEL:   ktr_tox %3.8f \n",ktr_tox);
            printf("/////      MODEL: poly_rcut %3.8f \n",poly_rcut);
        };
        if (michael_poly_yes_no == 0) {
            if (ktr_tox == 0.0)   {printf("/////      MODEL:LA       \n");};
            if (ktr_tox != 0.0)   {printf("/////      MODEL:LA + TOX \n");};
        };
        printf("///////////////////////////////////////////////////////////////////////////////\n");

        if (write_analyze==1){
        write_out_analyze_positions_distances_from_equilibrium_xyzmean();};

        printf("Energy std (DFT-H)   : %6.4f meV/atom (works in general )\n",dudl[zeitschritte-1].dudl_dft_harm_std);
        printf("Energy std (DFT-LA)  : %6.6f meV/atom (works in general lon and tox, eq. to michael)\n",dudl[zeitschritte-1].dudl_dft_la_std);
        printf("\n");
        //printf("Forces std <DFT-H>   : %6.3f (eV/A) delta angle degree: \n",sqrt(fdudl_dft_harm_av*1000.));
        printf("Forces std <DFT-H>   : %6.5f (eV/A)\n",sqrt(fdiff_dft_harm_av));
        printf("Forces std <DFT-LA>  : %6.6f (eV/A) (ddof = 0)\n",sqrt(fdiff_dft_la_av));
        printf("Forces std <DFT-LA>  : %6.6f (eV/A) (ddof = 1)\n",sqrt(fdiff_dft_la_av));
        printf("Forces std <DFT>     : %6.5f (eV/A) (this inly considers DFT forces)\n",sqrt(fstd_dft));
        printf("Forces std <LA>      : %6.5f (eV/A) (this only considers LA  forces)\n",sqrt(fstd_la));
        printf("pearson correlation F: %6.6f (between LA and DFT)\n",fcov_dft_la/(sqrt(fstd_dft)*sqrt(fstd_la)));
        printf("\n");
        printf("Energy dudl (DFT-H)  : %6.4f meV/atom (works in general )\n",dudl_dft_harm_av);
        printf("Energy dudl (DFT-LA) : %6.4f meV/atom (works in general lon and tox)\n",dudl_dft_la_av);
        printf("\n");
        printf("MAKE FOR EVERY STEP: rmin, rmax, deltaF_max(between LA and DFT)\n");
        printf("\n");
        printf("rmin           %5.3f Angstrom; %5.3f relaive; rmin_distmax: %5.3f Angstrom\n",rmin*da0_alat_lattice,rmin*da0_alat_lattice/alat_lattice,alat_lattice/sqrt(2)-rmin*da0_alat_lattice);
        printf("nndist         %5.3f Angstrom; %5.3f relaive;\n",alat_lattice/sqrt(2),alat_lattice/sqrt(2)/alat_lattice);
        printf("rmax           %5.3f Angstrom; %5.3f relaive; rmax_distmax: %5.3f Angstrom\n",rmax*da0_alat_lattice,rmax*da0_alat_lattice/alat_lattice,-1.*(alat_lattice/sqrt(2)-rmax*da0_alat_lattice));


        if (write_analyze==1){
            fclose(file_forces_vs_forces_dft);
        };
	    fclose(file_out_dudl);
	    fclose(file_out_dudl_av);
	    if (write_analyze_forces==1) {
	        fclose(file_forces);
	        fclose(file_forces_av);
	        fclose(file_forces_vs_forces_dft);}
        if (write_positions==1) {fclose(file_out_positions);}
        if (write_positions_forces==1) {fclose(file_out_positions);}
    }

    //////////////////////////////////
    // MD using VERLET (NO POSITIONs)
    //////////////////////////////////
    if (read_pos==0) {
        printf("///////////////////////////////////////////////////////////////////////////////\n");
        printf("/////      MD using VERLET %d steps\n",zeitschritte);
        printf("///////////////////////////////////////////////////////////////////////////////\n");
	    file_out_temp     =(FILE *)fopen("out_temperature.dat","w");
	    file_out_temp_av  =(FILE *)fopen("out_temperatureav.dat","w");




	    if (read_hesse==1) {
            check_if_file_exists(filename_in_hesse);
            read_hessematrix(hessemat,read_hesse,filename_in_hesse,verbose);
            get_pos0_from_externalfile_if_readpos(1,filename_in_positions0,columns0);  // get pos0
            }

	    if (read_pos0==1) {
            get_pos0_from_externalfile_if_readpos(1,filename_in_positions0,columns0);}

        check_pos0_for_duplicates(faktor,filename_in_positions);
        check_l1nn_for_pos0_and_resort_1nn();
        check_pos0_for_duplicates(faktor,filename_in_positions);
        //print_pos_0_to_screen(12, "POS_0");
        //print_pos_to_screen(12, "POS",i);
        check_pos_vs_pos0_if_same_positions_and_check_if_correct_alat();




        preequilibration_if_no_readpos(dt,read_pos,file_out_temp,hessemat,faktor_verlet,read_hesse,evolve_md_on_hesse);

        printf("... MD start ... for %d steps\n",zeitschritte);
        for (i=0;i<(zeitschritte+1);i++) {
            //printf("\n");
            //printf("----> i: %d %d\n",i,zeitschritte);

            calculate_forces_energy_la_from_simplelist(dt,verbose);
            //printf("----> j:\n");
            //if (read_hesse==1) {calculate_forces_energy_hesse(hessemat);}
            if (evolve_md_on_hesse==1) {forcesharm_to_velocities(dt);}
            else {forces_to_velocities(dt);}
            //printf("----> k:\n");


            if (verbose>=2) {
                print_pos_forcesdft_to_screen(4,"======> POS FORCES DFT=0 ",faktor,i);
                print_forces_to_screen(4,"======> LA FORCES ",faktor,i);
                print_forcesharm_to_screen(4,"======> Harmonic FORCES ",faktor,i);
                print_velocities_to_screen(4,"======> LA VELOCITIES ",faktor,i);
                // pos          : pos[j].x*faktor
                //
                // forcedft     : forcedft[j].x
                // forcela      : force[j].x*faktor_force // faktor_force = 1/1.602176e-9
                //              : force[j].x/1.602176e-9
                // forceharm    : forceharm[j].x*faktor
                //
                // vel dft      : DONT have those currently but could be checked from forces
                // vel la       : forcela*dtm  (dtm = dt/m_element;)
                // vel harm     : DONT have those currently
                //
                // temperature  : from velocities
                // for verlet   : from velocities
                }
            if (verbose>=2) {
                printf("Step %d   \n",i);
                printf("           DFT   %10.3f\n",arr_u_dft[i]);
                printf("           LA    %10.3f\n",u_la);
                printf("           HARM  %10.3f\n",u_harm);
                printf("\n");
                printf("\n");
                }

            //printf("----> l:\n");
            if (write_positions_forces==1) {
                if (i%l==l-1) {write_positions_forces_tofile(file_out_positions);};}
            //printf("----> m:\n");
            if (write_positions==1) {
                if (i%l==l-1) {write_positions_tofile(file_out_positions,write_positions_rel);};}

            if (write_analyze==1) {
                write_analyze_positions(i);
                };


            //printf("----> n:\n");
            get_temperature_from_velocities();                 // write it for every step?
            //printf("----> i: %d %4.10f\n",i,temperature);
            write_temperature(file_out_temp,i);  // write it for every step?
            //printf("----> o:\n");
            write_temperature_average(file_out_temp_av);  // write it for every step?

            //printf("----> p: ater that problem with write dudl\n");
            //printf("----> q:\n");
            //if (i%l==0) { // this happens only every tenth step or so
            //    //printf("----> icc: %d\n",i);
            //    get_temperature_from_velocities();
            //    //u_la_sum+=u_la; //*faktor_energy_atom;
            //    //u_harm_sum+=u_harm;
            //    write_temperature(file_out_temp); // or write only every 10th step?
            //};
            //printf("----> r: %5.10f\n",u_la_sum);
            //printf("----> s: %5.10f\n",u_harm_sum);

            verletstep_from_velocities(dt,file_out_temp,faktor_verlet);
            //printf("----> t:\n");
        }
    printf("... MD done !\n");
    }


    if (write_positions==1) {fclose(file_out_positions);}

	fclose(file_out_temp);
	fclose(file_out_temp_av);
	if (write_analyze==1){
	    fclose(file_out_check_dist_xyz);
	    fclose(file_out_check_dist_r);
	    fclose(file_out_check_dist_nn);
	    fclose(file_out_check_dist_nn_proj);
	    fclose(file_out_new1);
	    fclose(file_out_new2);
	    fclose(file_out_prl15_2a);
	    fclose(file_out_prl15_2au);
    }
    if (read_pos == 1) {fclose(file_in_positions);}

	end = clock();
    write_out_info_add_timing(start,end,"STOP");

    printf("... return 0 !\n");
	return 0;
}
