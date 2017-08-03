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

// The System class contains all the molecules and the methods to
// initialize and perform the simulation.
class System 
{
 private:
  static const int nSites = 4;      // the number of site in each molecule
  static const Vector site[nSites]; // the sites in the molecule 
  static const Vector I;            // the moments of inertia of a molecule
  static const double m;            // the mass of a molecule
  static double cubic[4];           // the interpolation coefficients
  static TopRecur top;              // see tops.h and doctops.pdf

  int N;          // the number of particles 
  double L;       // the linear length of the (cubic) system 
  double RUNTIME; // the runtime 
  double DT;      // the integration time step 
  Molecule*mol;   // the molecules 

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

  // To take a single time-step 
  void timeStep() {
    // force propagatation by half a step 
    for (int i=0; i<N; i++) {
      mol[i].p += 0.5*DT*mol[i].F;
      mol[i].L += 0.5*DT*mol[i].t;
    }
    // free propagation by a full step 
    for (int i=0; i<N; i++ ) {
      mol[i].q += DT*mol[i].p/m;
      Vector omega = mol[i].A*mol[i].L; 
      omega.x /= I.x;
      omega.y /= I.y; 
      omega.z /= I.z; 
      top.Propagation(DT, omega, mol[i].A); 
    }
    // recompute forces 
    computeForces();
    // force propagatation by another half step 
    for (int i=0; i<N; i++) {
      mol[i].p += 0.5*DT*mol[i].F;
      mol[i].L += 0.5*DT*mol[i].t;
    }
  }

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
  void initialize(const char*filename) {
    // read in parameters 
    ifstream f(filename); 
    if (!f) {
      cerr<<"Cannot find file "<<filename<<"\nExiting.\n";
      exit(1);
    }
    f >> N;
    f >> L;
    f >> RUNTIME;
    f >> DT;
    // initialize the system 
    mol = new Molecule[N];
    // generate random initial conditions
    for (int i=0; i<N; i++) {
      // random positions 
      mol[i].q.x = rnd()-0.5; 
      mol[i].q.y = rnd()-0.5;
      mol[i].q.z = rnd()-0.5;
      // make positions range from -L/2 to L/2 
      mol[i].q *= L;          
      // put each molecule initially in upright attitude 
      mol[i].A.one();         
      mol[i].p.x = rnd()-0.5;
      // random linear momenta 
      mol[i].p.y = rnd()-0.5;
      mol[i].p.z = rnd()-0.5;
      // scale with mass
      mol[i].p *= sqrt(m);    
      // draw random angular momenta scaled with respective
      // inertial moments
      mol[i].L.x = sqrt(I.x)*(rnd()-0.5); 
      mol[i].L.y = sqrt(I.y)*(rnd()-0.5);
      mol[i].L.z = sqrt(I.z)*(rnd()-0.5);
    }
    computeForces();
  }
  
  // To run the simulation while outputting.
  void run() {
    report(0);
    for (double t=0; t<RUNTIME; t+=DT) {
      timeStep();
      report(t);
    }
  }
};

// Define the system's parameters:

const Vector System::site[nSites] = {Vector(0,.075,.05), Vector(0,-.075,.05),
                                     Vector(.075,0,-.1), Vector(-.075,0,-.1)};
const Vector System::I(0.3375,0.45,0.525);
const double System::m  = 60.0;            
TopRecur System::top(System::I.x,System::I.y,System::I.z);
double System::cubic[4];

// Main function called at start up.

int main(int argc, char**argv) 
{
  char def_ini[] = "simtop.ini";  // the default <inifile>'s name
  System system;                  // the system object

  // initialize
  if (argc>1) system.initialize(argv[1]);
      else    system.initialize(def_ini);

  // run the simulation
  system.run();
}
