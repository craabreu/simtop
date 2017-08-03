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
#include "tops.h"
#include "emdee.h"

using namespace std; 

// To obtain a random number between 0 and 1
double rnd() 
{
  return rand()/(double)RAND_MAX; 
}

// The Molecule structure contains the properties of one molecule.
struct Molecule 
{
  Vector q; // the position 
  Vector p; // the momentum 
  Matrix A; // the attitude 
  Vector L; // the angular momentum 
  Vector F; // the force 
  Vector t; // the torque 
  double V; // the potential energy 
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

  Matrix mat;

  mat.xx = w2+i2-j2-k2;
  mat.xy = twoij-twokw;
  mat.xz = twojw+twoik;

  mat.yx = twoij+twokw;
  mat.yy = w2-i2+j2-k2;
  mat.yz = twojk-twoiw;

  mat.zx = twoik-twojw;
  mat.zy = twojk+twoiw;
  mat.zz = w2-i2-j2+k2;

  return mat;
}


// The System class contains all the molecules and the methods to
// initialize and perform the simulation.
class System 
{
 public:
  int N;          // Number of molecules
  double L;       // Simulation box side length
  double Rc;      // Cutoff distance
  double Rs;      // Neighbor list skin
  double Temp;    // Temperature
  int seed;       // Seed for random numbers
  int Nsteps;     // Number of steps
  double Dt;      // Integration time step
  int Nprop;      // Interval for printing properties
  double mvv2e;   // Energy convertion factor
  double Pconv;   // Pressure convertion factor
  double kB;      // Boltzmann's constant
  double kCoul;   // Coulomb's constant

  tEmDee md;      // EmDee system

  static const int nSites = 3;      // the number of sites in each molecule
  static const int nTypes = 2;      // the number of site types
  static const Vector site[nSites]; // the sites in the molecule
  static const int type[nSites];    // the type of each site

 private:
  static const Vector I;            // the moments of inertia of a molecule
  static const double m;            // the mass of a molecule
  static const double epsilon[nTypes];
  static const double sigma[nTypes];
  static double cubic[4];           // the interpolation coefficients
  static TopRecur top;              // see tops.h and doctops.pdf

  double RUNTIME; // the runtime 
  Molecule*mol;   // the molecules 

  //------------------------------------------------------------------------------------------------
  // To pre-compute interpolation coefficients.
  static const double RCA  = 2.5;
  static const double RCB  = 4;
  static const double RCA2 = 6.25;
  static const double RCB2 = 16;
  void computeCubicCoefficients() {
    double ir2 = 1.0/(double)RCA2;
    double ir6 = ir2*ir2*ir2;
    double V   = ir6*(4*ir6-4);
    double F   = ir6*(48*ir6-24)/RCA;
    double a   = RCB-RCA;
    double b   = 3*V-a*F;
    double c   = -(2*V-a*F);
    double ia2 = 1/(a*a);
    double ia3 = ia2/a;
    cubic[0] = b*RCB2*ia2+c*RCB*RCB2*ia3;
    cubic[1] = -2*b*RCB*ia2-3*c*RCB2*ia3;
    cubic[2] = b*ia2+3*c*RCB*ia3;
    cubic[3] = -c*ia3;
  }

  //------------------------------------------------------------------------------------------------
  // To compute the site-site potential with a smooth cut-off.
  // For r < RCA, full Lennard-Jones is used;
  // for r > RCB, the potential and force are both zero;
  // for RCA < r < RCB, an interpolating cubic is used.
  void siteForce(const Vector&r, Vector&f, double&V) const {
    double r2 = r.nrm2();
    if (r2>RCB2) { 
      V = 0.0; 
      f.zero(); 
    } else if (r2>RCA2) { 
      // smooth interpolation 
      double rnrm = sqrt(r2);
      V = cubic[0]+rnrm*cubic[1]+r2*(cubic[2]+rnrm*cubic[3]);
      f = -(cubic[1]/rnrm+2*cubic[2]+3*rnrm*cubic[3])*r;
    } else { 
      double ir2 = 1/r2;
      double ir6 = ir2*ir2*ir2;
      V = ir6*(4*ir6-4);
      f = ir6*ir2*(48*ir6-24)*r;
    } 
  }

  //------------------------------------------------------------------------------------------------
  // To compute the force F and torque t on molecule B due to A.
  void molForce(const Molecule&A, const Molecule&B, const Vector&r,
                 Vector&F, Vector&t, double&V) const {
     V = 0;          
     F.zero();
     t.zero();   
     for (int a=0; a<nSites; a++) 
       for (int b=0; b<nSites; b++) { 
         Vector ra  = Transpose(A.A)*site[a];
         Vector rb  = Transpose(B.A)*site[b];
         Vector rba = r+rb-ra;
         Vector fba;
         double vba;
         siteForce(rba, fba, vba);         
         V += vba;
         F += fba;
         t += rb^fba;
      }
  }

  //------------------------------------------------------------------------------------------------
  // To compute the forces between the parts of all molecules.
  void computeForces() {
    for (int i=0; i<N; i++){
      mol[i].F.zero();
      mol[i].t.zero();
      mol[i].V = 0;
    }    
    for (int i=0; i<N; i++) 
      for (int j=0; j<i; j++) {
        Vector r=mol[j].q-mol[i].q;  	
          // apply periodic boundary conditions 
          r.x -= L*rint(r.x/L);
          r.y -= L*rint(r.y/L);
          r.z -= L*rint(r.z/L);
          // compute force and torque on j due to i 
          double V;    
          Vector F;
          Vector t;
          molForce(mol[i], mol[j], r, F, t, V);            
          // add to total forces and torques
          mol[j].F += F;
          mol[i].F -= F;
          mol[j].t += t;
          // torque on i is -torque on j, - r x Fj: 
          mol[i].t -= t+(r^F); 
          mol[i].V += 0.5*V;
          mol[j].V += 0.5*V;
        }
  }

  //------------------------------------------------------------------------------------------------
  // To take a single time-step 
  void timeStep() {
    // force propagatation by half a step 
    for (int i=0; i<N; i++) {
      mol[i].p += 0.5*Dt*mol[i].F;
      mol[i].L += 0.5*Dt*mol[i].t;
    }
    // free propagation by a full step 
    for (int i=0; i<N; i++ ) {
      mol[i].q += Dt*mol[i].p/m;
      Vector omega = mol[i].A*mol[i].L; 
      omega.x /= I.x;
      omega.y /= I.y; 
      omega.z /= I.z; 
      top.Propagation(Dt, omega, mol[i].A); 
    }
    // recompute forces 
    computeForces();
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
    double Epot = 0;
    for (int i=0; i<N; i++) {
      Vector Lb = mol[i].A*mol[i].L;
      Epot += mol[i].V;
      Elin += 0.5*mol[i].p.nrm2()/m;
      Erot += 0.5*(Lb.x*Lb.x/I.x+Lb.y*Lb.y/I.y+Lb.z*Lb.z/I.z);
    }
    cout<<t<<' '<<Elin+Erot+Epot<<' '<<Elin<<' '<<Erot<<' '<<Epot<<'\n';
  }

 public:

  //------------------------------------------------------------------------------------------------
  // To construct the system.
  System() {
    mol = 0;
    // cause accurate output
    cout<<setprecision(14);   
    // precompute coefficients
    computeCubicCoefficients(); 
  }

  // To deconstruct the system.
  ~System() {
    // release memory
    delete[] mol;
  }

  // To initialize the system.
  void initialize(const char*filename, int threads) {
    FILE *file;
    file = fopen(filename,"r");
    #define readline \
      if (!fgets(line, sizeof(line), file)) { \
        cerr << "ERROR: could not read data.\n"; \
        exit(0); \
      }
    char line[256];
    if (file) {
      readline; readline; sscanf( line, "%d",  &N );
      readline; readline; sscanf( line, "%lf", &L );
      readline; readline; sscanf( line, "%lf", &Dt );
      readline; readline; sscanf( line, "%lf", &Rc );
      readline; readline; sscanf( line, "%lf", &Rs );
      readline; readline; sscanf( line, "%d",  &seed );
      readline; readline; sscanf( line, "%d",  &Nsteps );
      readline; readline; sscanf( line, "%d",  &Nprop );
      readline; readline; sscanf( line, "%lf", &Temp );
      readline; readline; sscanf( line, "%lf", &mvv2e );
      readline; readline; sscanf( line, "%lf", &Pconv );
      readline; readline; sscanf( line, "%lf", &kB );
      readline; readline; sscanf( line, "%lf", &kCoul );
      mol = new Molecule[N];
      readline;
      for (int i = 0; i < N; i++) {
        double quat[4];
        readline; sscanf( line, "%lf %lf %lf %lf %lf %lf %lf", 
                                &mol[i].q.x, &mol[i].q.y, &mol[i].q.z,
                                &quat[0], &quat[1], &quat[2], &quat[3] );
        mol[i].A = quat_to_mat(quat);
      }
    }
    #undef readline
    RUNTIME = Nsteps*Dt;

    int natoms = N*nSites;
    int body[natoms], atom_type[natoms];
    for (int i = 0; i < N; i++)
      for (int j = 0; j < nSites; j++) {
        body[nSites*i + j] = i+1;
        atom_type[nSites*i + j] = type[j];
      }

    md = EmDee_system( threads, 1, Rc, Rs, natoms, &atom_type[0], NULL, &body[0] );
  }
  
  //------------------------------------------------------------------------------------------------
  // To run the simulation while outputting.
  void run() {
    report(0);
    for (double t=0; t<RUNTIME; t+=Dt) {
      timeStep();
      report(t);
    }
  }
};

// Define the system's parameters for TIP3P Water:
const Vector System::site[nSites] = {Vector(-0.75695,-0.52031970, 0.0),  // A
                                     Vector( 0.00000, 0.06556274, 0.0),  // A
                                     Vector( 0.75695,-0.52031970, 0.0)}; // A
const int    System::type[nSites] = {1,2,1};
const double System::epsilon[nTypes] = {0.0};
const double System::sigma[nTypes] = {0.0};
const Vector System::I(0.61457,1.15511,1.76968); // Da.A^2
const double System::m  = 18.0154; // Da
TopRecur System::top(System::I.x,System::I.y,System::I.z);
double System::cubic[4];


// Main function called at start up.

int main(int argc, char**argv) 
{
  System system;                  // the system object

  if (argc == 2)
    system.initialize(argv[1],1);
  else if (argc == 3)
    system.initialize(argv[2],atoi(argv[1]));
  else {
    cerr<<"Usage: simtop [nthreads] filename\n";
    exit(1);
  }


//  char def_ini[] = "simtop.ini";  // the default <inifile>'s name

//  // initialize
//  if (argc>1) system.initialize(argv[1]);
//      else    system.initialize(def_ini);

  // run the simulation
//  system.run();
}
