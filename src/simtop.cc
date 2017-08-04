/*
  simtop.cc
 
  Simulation of a fluid of molecules each consisting of four
  Lennard-Jones sites.
 
  This program is meant for the purpose of illustrating the use of the
  implementation of the exact rotation of a rigid body in the context
  of the symplectic integration scheme using the exact rotation
  matrix.  The calculations of the forces and the random number
  generator are therefore implemented only in the most straightforward
  but computationally costly way.
 
  To compile type:
    c++ simtop.cc tops.cc -lm -o simtop
  
  To run type:
    simtop <inifile>

  Input: The input file <inifile> should contain a list of numbers
  which signify
 
     number of particles
     linear size of the cubic simulation box
     run-time of the simulation
     integration step-size 
 
  Output: The program writes a line per time step containing five
  columns representing time, total energy, translational kinetic
  energy, rotational kinetic energy and potential energy. 

  Ramses van Zon
  Last modified: May 29, 2009
*/

#include <iostream>
#include <iomanip>
#include <fstream>
#include <cstdlib>
#include <cmath>
#include <string.h>
#include "tops.h"
#include "emdee.h"

#define TRUE (_Bool)1
#define FALSE (_Bool)0

using namespace std; 

// To obtain a random number between 0 and 1
double rnd() 
{
  return rand()/(double)RAND_MAX; 
}

// The Molecule structure contains the properties of one molecule.
struct Molecule 
{
  Vector r;  // the position 
  Vector p;  // the momentum 
  Matrix A;  // the attitude 
  Matrix At; // the attitude transpose
  Vector L;  // the angular momentum 
  Vector F;  // the force 
  Vector t;  // the torque 
};

Matrix quat_to_mat(double quat[4])
{
  double w2 = quat[0]*quat[0];
  double i2 = quat[1]*quat[1];
  double j2 = quat[2]*quat[2];
  double k2 = quat[3]*quat[3];
  double twoij = 2.0*quat[1]*quat[2];
  double twoik = 2.0*quat[1]*quat[3];
  double twojk = 2.0*quat[2]*quat[3];
  double twoiw = 2.0*quat[1]*quat[0];
  double twojw = 2.0*quat[2]*quat[0];
  double twokw = 2.0*quat[3]*quat[0];
  return Matrix(w2+i2-j2-k2, twoij-twokw, twojw+twoik,
                twoij+twokw, w2-i2+j2-k2, twojk-twoiw,
                twoik-twojw, twojk+twoiw, w2-i2-j2+k2);
}

double cbox( double r, double L )
{
  return r - L*floor(r/L);
}

//==================================================================================================

// The System class contains all the molecules and the methods to
// initialize and perform the simulation.
class System 
{
 public:
  char *base;     // Base for file names
  int N;          // Number of molecules
  double L;       // Simulation box side length
  double Rc;      // Cutoff distance
  double Rs;      // Neighbor list skin
  double smooth;  // Skin for potential smoothing
  double Temp;    // Temperature
  int seed;       // Seed for random numbers
  int Nsteps;     // Number of steps
  double Dt;      // Integration time step
  int Nprop;      // Interval for printing properties
  double mvv2e;   // Energy convertion factor
  double Pconv;   // Pressure convertion factor
  double kB;      // Boltzmann's constant
  double kCoul;   // Coulomb's constant
  double alpha;   // Electrostatic damping parameter

  tEmDee md;      // EmDee system

  static const int    nSites = 3;       // the number of sites in each molecule
  static const int    nTypes = 2;       // the number of site types

 private:
  static const Vector site[nSites];     // the sites in the molecule
  static const int    type[nSites];     // the type of each site
  static const char   element[nSites];  // the chemical element at each site
  static const double charge[nSites];   // the chemical element at each site

  static const double mass[nTypes];     // the mass of each site
  static const double epsilon[nTypes];
  static const double sigma[nTypes];

  static const Vector I;            // the moments of inertia of a molecule
  static const double m;            // the mass of a molecule
  static TopRecur top;              // see tops.h and doctops.pdf

  double Epot;

  FILE *xyz;

  double RUNTIME; // the runtime 
  Molecule*mol;   // the molecules 

  //------------------------------------------------------------------------------------------------
  // To take a single time-step 
  void timeStep() {
    md.Energy.Compute = TRUE;
    // force propagatation by half a step 
    for (int i=0; i<N; i++) {
      mol[i].p += 0.5*Dt*mol[i].F;
      mol[i].L += 0.5*Dt*mol[i].t;
    }
    // free propagation by a full step 
    for (int i=0; i<N; i++ ) {
      mol[i].r += Dt*mol[i].p/m;
      Vector omega = mol[i].A*mol[i].L; 
      omega.x /= I.x;
      omega.y /= I.y; 
      omega.z /= I.z; 
      top.Propagation(Dt, omega, mol[i].A);
      mol[i].At = Transpose(mol[i].A);
    }
    // recompute forces 
    compute();
    // force propagatation by another half step 
    for (int i=0; i<N; i++) {
      mol[i].p += 0.5*Dt*mol[i].F;
      mol[i].L += 0.5*Dt*mol[i].t;
    }
  }

  //------------------------------------------------------------------------------------------------
  // To report energetic properties to standard output 
  void report(double t) const {
    double Elin = 0;
    double Erot = 0;
    for (int i=0; i<N; i++) {
      Vector Lb = mol[i].A*mol[i].L;
      Elin += 0.5*mol[i].p.nrm2()/m;
      Erot += 0.5*(Lb.x*Lb.x/I.x+Lb.y*Lb.y/I.y+Lb.z*Lb.z/I.z);
    }
    cout<<t<<' '<<mvv2e*(Elin+Erot+Epot)<<' '<<mvv2e*Elin<<' '<<mvv2e*Erot<<' '<<mvv2e*Epot<<'\n';
  }

  //------------------------------------------------------------------------------------------------
  // To compute forces and torques.
  void compute() {
    int natoms = N*nSites;
    double R[natoms][3], F[natoms][3];
    for (int i = 0; i < N; i++) {
      for (int j = 0; j < nSites; j++) {
        int k = nSites*i + j;
        Vector r = mol[i].r + mol[i].At*site[j];
        R[k][0] = r.x;
        R[k][1] = r.y;
        R[k][2] = r.z;
      }
    }
    EmDee_upload( &md, (char*)"coordinates", &R[0][0] );
    Epot = md.Energy.Potential;
    EmDee_download( md, (char*)"forces", &F[0][0] );
    for (int i = 0; i < N; i++) {
      mol[i].F.zero();
      mol[i].t.zero();
      for (int j = 0; j < nSites; j++) {
        int k = i*nSites + j;
        Vector f = Vector(F[k][0],F[k][1],F[k][2]);
        mol[i].F += f;
        mol[i].t += (mol[i].At*site[j])^f;
      }
    }
  }

 public:

  //------------------------------------------------------------------------------------------------
  // To construct the system.
  System() {
    mol = 0;
    // cause accurate output
    cout<<setprecision(14);   
  }

  //------------------------------------------------------------------------------------------------
  // To deconstruct the system.
  ~System() {
    // release memory
    delete[] mol;
    fclose(xyz);
  }

  //------------------------------------------------------------------------------------------------
  // To initialize the system.
  void initialize(string filebase, int threads) {
    string name = filebase + ".inp";
    FILE *file = fopen(name.c_str(),"r");
    if (!file) {
      cerr << "Error: file " << filebase << ".inp not found\n";
      exit(1);
    }
    char line[256];
    #define readline \
      if (!fgets(line, sizeof(line), file)) { \
        cerr << "ERROR: could not read data.\n"; \
        exit(0); \
      }
    readline; readline; sscanf( line, "%d",  &N );
    readline; readline; sscanf( line, "%lf", &L );
    readline; readline; sscanf( line, "%lf", &Dt );
    readline; readline; sscanf( line, "%lf", &Rc );
    readline; readline; sscanf( line, "%lf", &Rs );
    readline; readline; sscanf( line, "%lf", &smooth );
    readline; readline; sscanf( line, "%d",  &seed );
    readline; readline; sscanf( line, "%d",  &Nsteps );
    readline; readline; sscanf( line, "%d",  &Nprop );
    readline; readline; sscanf( line, "%lf", &Temp );
    readline; readline; sscanf( line, "%lf", &mvv2e );
    readline; readline; sscanf( line, "%lf", &Pconv );
    readline; readline; sscanf( line, "%lf", &kB );
    readline; readline; sscanf( line, "%lf", &kCoul );
    readline; readline; sscanf( line, "%lf", &alpha );
    readline;
    mol = new Molecule[N];
    int natoms = N*nSites;
    int body[natoms], atom_type[natoms];
    double Q[natoms];
    for (int i = 0; i < N; i++) {
      double r[3], q[4];
      readline; sscanf( line, "%lf %lf %lf %lf %lf %lf %lf",
                        &r[0], &r[1], &r[2], &q[0], &q[1], &q[2], &q[3] );
      double qnorm = 0.0;
      for (int j = 0; j < 4; j++) qnorm += q[j]*q[j];
      mol[i].r = Vector(cbox(r[0],L),cbox(r[1],L),cbox(r[2],L));
      mol[i].A = quat_to_mat(q);
      mol[i].At = Transpose(mol[i].A);
      for (int j = 0; j < nSites; j++) {
        int k = nSites*i + j;
        body[k] = i+1;
        atom_type[k] = type[j];
        Q[k] = charge[j];
      }
    }
    #undef readline
    name = filebase + ".xyz";
    xyz = fopen(name.c_str(),"w");

    // Initialize EmDee system:
    md = EmDee_system( threads, 1, Rc, Rs, natoms, &atom_type[0], (double*)&mass[0], &body[0] );
    md.Options.AutoForceCompute = TRUE;
    md.Options.AutoBodyUpdate = FALSE;

    for (int i = 0; i < nTypes; i++) {
      void *pair;
      if (epsilon[i] == 0.0)
        pair = EmDee_pair_none();
      else
        pair = EmDee_smoothed( EmDee_pair_lj_cut( epsilon[i]/mvv2e, sigma[i] ), Rc-smooth );
      EmDee_set_pair_model( md, i+1, i+1, pair, kCoul );
    }
    EmDee_set_coul_model( md, EmDee_coul_damped_smoothed( alpha, smooth ) );
    EmDee_upload( &md, (char*)"charges", &Q[0] );
    EmDee_upload( &md, (char*)"box", &L );
    compute();
    EmDee_random_momenta( &md, kB*Temp, TRUE, seed );
    double p[natoms][3];
    EmDee_download( md, (char*)"momenta", &p[0][0] );
    for (int i = 0; i < N; i++) {
      mol[i].p.zero();
      mol[i].L.zero();
      Matrix At = Transpose(mol[i].A);
      for (int j = 0; j < nSites; j++) {
        int k = i*nSites + j;
        mol[i].p += Vector(p[k][0],p[k][1],p[k][2]);
        mol[i].L += (At*site[j])^Vector(p[k][0],p[k][1],p[k][2]);
      }
    }
  }

  //------------------------------------------------------------------------------------------------
  // To write down an xyz file:
  void write_xyz() {
    fprintf(xyz,"%d\n\n",N*nSites);
    for (int i = 0; i < N; i++) {
      Vector Ri = mol[i].r;
      Ri.x = cbox(Ri.x,L);
      Ri.y = cbox(Ri.y,L);
      Ri.z = cbox(Ri.z,L);
      for (int j = 0; j < nSites; j++) {
        Vector r = Ri + Transpose(mol[i].A)*site[j];
        fprintf(xyz,"%c %f %f %f\n",element[j],r.x,r.y,r.z);
      }
    }
  }

  //------------------------------------------------------------------------------------------------
  // To run the simulation while outputting.
  void run() {
    report(0);
    for (int step=0; step < Nsteps; step++) {
      timeStep();
      report(step*Dt);
    }
  }
};

//==================================================================================================

// Define the system's parameters for TIP3P Water:
const Vector System::site[nSites] = {Vector(-0.75695,-0.52031970, 0.0),  // A
                                     Vector( 0.00000, 0.06556274, 0.0),  // A
                                     Vector( 0.75695,-0.52031970, 0.0)}; // A

const int    System::type[nSites] = {1,2,1};
const char   System::element[nSites] = {'H','O','H'};
const double System::charge[nSites] = {0.417,-0.834,0.417};

const double System::mass[nTypes] =    {1.008, 15.9994};  // Da
const double System::epsilon[nTypes] = {  0.0,  0.1520};  // kcal/mol
const double System::sigma[nTypes] =   {  0.0,  3.1507};  // A

const Vector System::I(0.61457,1.15511,1.76968); // Da.A^2
const double System::m  = 18.0154; // Da
TopRecur System::top(System::I.x,System::I.y,System::I.z);

//--------------------------------------------------------------------------------------------------
// Main function called at start up.

int main(int argc, char**argv) 
{
  System system;                  // the system object

  if (argc == 2)
    system.initialize((string)argv[1],1);
  else if (argc == 3)
    system.initialize((string)argv[2],atoi(argv[1]));
  else {
    cerr<<"Usage: simtop [nthreads] filebase\n";
    exit(1);
  }

  system.write_xyz();

  // run the simulation
  system.run();
}
