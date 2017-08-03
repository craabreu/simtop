/*
   tops.cc
 
   Implementation of the 'Top' classes defined in tops.h.
 
   Ramses van Zon
 
   version 2, May 28, 2009

   - Includes computation of elliptic functions/integrals, so the gsl library is not needed.
   - Contains optional drift-correction. Define REFINE to enable (at a computational cost!)
   - An error in the TopProlate was corrected.
*/
#include <limits>
#include <cstdio>
#include <cstdlib>
#include "tops.h"

/*
 *  The AgmScale class uses the technique of the arithmetic-geometric
 *  mean (agm) scale to compute elliptic functions and integrals used
 *  by the TopAsymmetric and TopRecur classes.
 */
class AgmScale 
{
 public:
  AgmScale();
  void precompute( DOUBLE m );
  void sncndn( DOUBLE x, DOUBLE & sn, DOUBLE & cn, DOUBLE & dn ) const;
  DOUBLE F() const;
  DOUBLE F( DOUBLE phi ) const;
  DOUBLE F_from_sin( DOUBLE sinphi ) const;
  DOUBLE F_from_tan( DOUBLE tanphi ) const;    
 private:
  enum { MAXN = 13 }; 
  int     N;    
  DOUBLE  m;    
  DOUBLE  a[MAXN+1], b[MAXN+1], d[MAXN+1];
  DOUBLE  K;                              
  DOUBLE  Kscale;                         
};

/* Useful macros: */
#define odd(i)      ((bool)((i)&1))
#define sqr(x)      ((x)*(x))
#define maxNTerms   18

/* Compute the sin and cosine of an angle phi. */
inline void mysincos(const float phi, float & sinphi, float & cosphi) {
  sincosf(phi, &sinphi, &cosphi);
}
inline void mysincos(double phi, double & sinphi, double & cosphi) {
  sincos(phi, &sinphi, &cosphi);
}
inline void mysincos(const long double phi, long double & sinphi, long double & cosphi) {
  sinphi = sin(phi);
  cosphi = cos(phi);
}

#ifdef REFINE
/* Enforce that the squares of two numbers add up to one.
 * x is the first number (a reference to a DOUBLE)
 * y is the second number (a reference to a DOUBLE)
 */
inline void refine(DOUBLE & x, DOUBLE & y) {
  if (fabs(x)>fabs(y))
    x = copysign(sqrt(1-sqr(y)),x); 
  else
    y = copysign(sqrt(1-sqr(x)),y); 
}
/* 
 *  Enforce that the energy is equal to a target value.  The energy
 *  deviates from the target value due to round-off effects in the
 *  least significant bits of the numbers, leading to drift.  The
 *  energy drift is minimized in this routine by flipping the least
 *  significant bits, up to a maximum of maxBits bits.  Note that the
 *  energy is defined here as the squares of three angular velocities
 *  times their respective inertial moments (this is really twice the
 *  kinetic energy).
 *
 *  omega1 is the first number (a reference to a DOUBLE) 
 *  omega2 is the second number (a reference to a DOUBLE) 
 *  omega3 is the third number (a reference to a DOUBLE) 
 *  I1 is the first inertial moment
 *  I2 is the second inertial moment
 *  I3 is the third inertial moment
 *  Etarget the target energy
 */
inline void refineBits(DOUBLE & omega1, DOUBLE & omega2, DOUBLE & omega3, 
                       DOUBLE I1, DOUBLE I2, DOUBLE I3, DOUBLE Etarget) {
  /* On the choice of maxBits: 
   * - maxBits may range from zero (which does nothing) to nine.
   * - Too large a value can destroy conservation of angular momentum.
   * - Too small a value does not significantly reduce the drift.  
   * - From experience, we found maxBits = 4 to be reasonable overall.
   * - The true optimal value depends on the machine and the compiler.
   */
  static const unsigned char maxBits = 4;
  static const unsigned char indexCombination[220][3] = {
  { 0,0,0 }, { 1,0,0 }, { 0,1,0 }, { 0,0,1 }, { 2,0,0 }, { 1,1,0 }, { 1,0,1 }, { 0,2,0 }, { 0,1,1 }, { 0,0,2 },
  { 3,0,0 }, { 2,1,0 }, { 2,0,1 }, { 1,2,0 }, { 1,1,1 }, { 1,0,2 }, { 0,3,0 }, { 0,2,1 }, { 0,1,2 }, { 0,0,3 },
  { 4,0,0 }, { 3,1,0 }, { 3,0,1 }, { 2,2,0 }, { 2,1,1 }, { 2,0,2 }, { 1,3,0 }, { 1,2,1 }, { 1,1,2 }, { 1,0,3 },
  { 0,4,0 }, { 0,3,1 }, { 0,2,2 }, { 0,1,3 }, { 0,0,4 }, { 5,0,0 }, { 4,1,0 }, { 4,0,1 }, { 3,2,0 }, { 3,1,1 },
  { 3,0,2 }, { 2,3,0 }, { 2,2,1 }, { 2,1,2 }, { 2,0,3 }, { 1,4,0 }, { 1,3,1 }, { 1,2,2 }, { 1,1,3 }, { 1,0,4 },
  { 0,5,0 }, { 0,4,1 }, { 0,3,2 }, { 0,2,3 }, { 0,1,4 }, { 0,0,5 }, { 6,0,0 }, { 5,1,0 }, { 5,0,1 }, { 4,2,0 },
  { 4,1,1 }, { 4,0,2 }, { 3,3,0 }, { 3,2,1 }, { 3,1,2 }, { 3,0,3 }, { 2,4,0 }, { 2,3,1 }, { 2,2,2 }, { 2,1,3 },
  { 2,0,4 }, { 1,5,0 }, { 1,4,1 }, { 1,3,2 }, { 1,2,3 }, { 1,1,4 }, { 1,0,5 }, { 0,6,0 }, { 0,5,1 }, { 0,4,2 },
  { 0,3,3 }, { 0,2,4 }, { 0,1,5 }, { 0,0,6 }, { 7,0,0 }, { 6,1,0 }, { 6,0,1 }, { 5,2,0 }, { 5,1,1 }, { 5,0,2 },
  { 4,3,0 }, { 4,2,1 }, { 4,1,2 }, { 4,0,3 }, { 3,4,0 }, { 3,3,1 }, { 3,2,2 }, { 3,1,3 }, { 3,0,4 }, { 2,5,0 },
  { 2,4,1 }, { 2,3,2 }, { 2,2,3 }, { 2,1,4 }, { 2,0,5 }, { 1,6,0 }, { 1,5,1 }, { 1,4,2 }, { 1,3,3 }, { 1,2,4 },
  { 1,1,5 }, { 1,0,6 }, { 0,7,0 }, { 0,6,1 }, { 0,5,2 }, { 0,4,3 }, { 0,3,4 }, { 0,2,5 }, { 0,1,6 }, { 0,0,7 },
  { 8,0,0 }, { 7,1,0 }, { 7,0,1 }, { 6,2,0 }, { 6,1,1 }, { 6,0,2 }, { 5,3,0 }, { 5,2,1 }, { 5,1,2 }, { 5,0,3 },
  { 4,4,0 }, { 4,3,1 }, { 4,2,2 }, { 4,1,3 }, { 4,0,4 }, { 3,5,0 }, { 3,4,1 }, { 3,3,2 }, { 3,2,3 }, { 3,1,4 },
  { 3,0,5 }, { 2,6,0 }, { 2,5,1 }, { 2,4,2 }, { 2,3,3 }, { 2,2,4 }, { 2,1,5 }, { 2,0,6 }, { 1,7,0 }, { 1,6,1 },
  { 1,5,2 }, { 1,4,3 }, { 1,3,4 }, { 1,2,5 }, { 1,1,6 }, { 1,0,7 }, { 0,8,0 }, { 0,7,1 }, { 0,6,2 }, { 0,5,3 }, 
  { 0,4,4 }, { 0,3,5 }, { 0,2,6 }, { 0,1,7 }, { 0,0,8 }, { 9,0,0 }, { 8,1,0 }, { 8,0,1 }, { 7,2,0 }, { 7,1,1 },
  { 7,0,2 }, { 6,3,0 }, { 6,2,1 }, { 6,1,2 }, { 6,0,3 }, { 5,4,0 }, { 5,3,1 }, { 5,2,2 }, { 5,1,3 }, { 5,0,4 },
  { 4,5,0 }, { 4,4,1 }, { 4,3,2 }, { 4,2,3 }, { 4,1,4 }, { 4,0,5 }, { 3,6,0 }, { 3,5,1 }, { 3,4,2 }, { 3,3,3 },
  { 3,2,4 }, { 3,1,5 }, { 3,0,6 }, { 2,7,0 }, { 2,6,1 }, { 2,5,2 }, { 2,4,3 }, { 2,3,4 }, { 2,2,5 }, { 2,1,6 }, 
  { 2,0,7 }, { 1,8,0 }, { 1,7,1 }, { 1,6,2 }, { 1,5,3 }, { 1,4,4 }, { 1,3,5 }, { 1,2,6 }, { 1,1,7 }, { 1,0,8 },
  { 0,9,0 }, { 0,8,1 }, { 0,7,2 }, { 0,6,3 }, { 0,5,4 }, { 0,4,5 }, { 0,3,6 }, { 0,2,7 }, { 0,1,8 }, { 0,0,9 } };
  static const unsigned char nCombinations[10] = { 1, 4, 10, 20, 35, 56, 84, 120, 165, 220 };

  DOUBLE Efinal = I1*sqr(omega1)+I2*sqr(omega2)+I3*sqr(omega3);

  if (Efinal == Etarget) 
    return;

  DOUBLE E1[maxBits+1];
  DOUBLE E2[maxBits+1];
  DOUBLE E3[maxBits+1];
  DOUBLE newOmega1[maxBits+1];
  DOUBLE newOmega2[maxBits+1];
  DOUBLE newOmega3[maxBits+1];

  for (unsigned char mask = 0; mask <= maxBits; mask++) {
    newOmega1[mask] = omega1;
    newOmega2[mask] = omega2;
    newOmega3[mask] = omega3;
    *(unsigned char*)&(newOmega1[mask]) ^= mask;
    *(unsigned char*)&(newOmega2[mask]) ^= mask;
    *(unsigned char*)&(newOmega3[mask]) ^= mask;
    E1[mask] = I1*sqr(newOmega1[mask]);
    E2[mask] = I2*sqr(newOmega2[mask]); 
    E3[mask] = I3*sqr(newOmega3[mask]);
  }

  for (unsigned char p = 1; p < nCombinations[maxBits]; p++) {
    const unsigned char i = indexCombination[p][0];
    const unsigned char j = indexCombination[p][1];
    const unsigned char k = indexCombination[p][2];
    DOUBLE Etrial = E1[i] + E2[j] + E3[k];
    
    if ( fabs(Etrial-Etarget) < fabs(Efinal-Etarget) ) {
      Efinal = Etrial;
      omega1 = newOmega1[i];
      omega2 = newOmega2[j];
      omega3 = newOmega3[k];
    }
        
    if (Efinal == Etarget) 
        break;
  }
}
#else
#define refine(x,y)
#define refineBits(omega1,omega2,omega3,DOUBLEI1,I2,I3,Etarget)
#endif

TopSpherical::TopSpherical(DOUBLE _I) : 
  I(_I) 
{}

void TopSpherical::Initialization(const Vector & omega, const Matrix & A) {
  omegab = omega;
  A0     = A;
}

void TopSpherical::Evolution(DOUBLE t, Vector & omega, Matrix & A) {
  omega = omegab;
  A     = Rodrigues( (-t)*omegab ) * A0; 
}


TopProlate::TopProlate(DOUBLE _Ix, DOUBLE _Iz) : 
  I1(_Iz), 
  I3(_Ix) 
{}

void TopProlate::Initialization(const Vector & omega, const Matrix & A) {
  omegab0    = omega;
  A0         = A;
  LboverI1.x = I3*omegab0.x/I1;
  LboverI1.y = omegab0.y;
  LboverI1.z = omegab0.z;
  omegap     = (1 - I3/I1)*omegab0.x;
}

void TopProlate::Evolution(DOUBLE t, Vector & omega, Matrix & A) {
  DOUBLE phi = omegap*t;
  DOUBLE s;
  DOUBLE c;

  mysincos(phi, s, c);

  omega.x = omegab0.x;
  omega.y =  c * omegab0.y + s * omegab0.z;
  omega.z = -s * omegab0.y + c * omegab0.z;

  A = Rodrigues( (-t)*LboverI1 ) * A0;

  DOUBLE d;
  d = A.zz;
  A.zz *= c;
  A.zz -= A.yz*s;
  A.yz *= c;
  A.yz += d*s;
  d =  A.zy;
  A.zy *= c;
  A.zy -= A.yy*s;
  A.yy *= c;
  A.yy += d*s;
  d =  A.zx;
  A.zx *= c;
  A.zx -= A.yx*s;
  A.yx *= c;
  A.yx += d*s;
}


TopOblate::TopOblate(DOUBLE _Ix, DOUBLE _Iz) : 
  I1(_Ix), 
  I3(_Iz)
{}

void TopOblate::Initialization(const Vector& omega, const Matrix& A) {
  omegab0    = omega;
  A0         = A;
  LboverI1.x = omegab0.x;
  LboverI1.y = omegab0.y;
  LboverI1.z = I3*omegab0.z/I1;
  omegap     = (1 - I3/I1)*omegab0.z;
}

void TopOblate::Evolution(DOUBLE t, Vector & omega, Matrix & A) {
  DOUBLE phi = omegap*t;
  DOUBLE s;
  DOUBLE c;
 
  mysincos(phi, s, c);

  omega.x =  c * omegab0.x + s * omegab0.y;
  omega.y = -s * omegab0.x + c * omegab0.y; 
  omega.z = omegab0.z;

  A = Rodrigues( (-t)*LboverI1 ) * A0;

  DOUBLE d;  // d is a place holder
  d = A.xx;
  A.xx *= c;
  A.xx += A.yx*s;
  A.yx *= c;
  A.yx -= d*s;
  d    =  A.xy;
  A.xy *= c;
  A.xy += A.yy*s;
  A.yy *= c;
  A.yy -= d*s;
  d    =  A.xz;
  A.xz *= c;
  A.xz += A.yz*s;
  A.yz *= c;
  A.yz -= d*s;
}

TopAsymmetricCommonPart::TopAsymmetricCommonPart(DOUBLE _Ix, DOUBLE _Iy, DOUBLE _Iz) :
  Ix(_Ix), 
  Iy(_Iy), 
  Iz(_Iz),
  agm(new AgmScale)
{}

TopAsymmetricCommonPart::~TopAsymmetricCommonPart() {
  delete agm;
}

TopAsymmetric::TopAsymmetric(DOUBLE _Ix, DOUBLE _Iy, DOUBLE _Iz) :
  TopAsymmetricCommonPart(_Ix, _Iy, _Iz),
  agmComplement(new AgmScale),
  maxNumTerms(3), 
  Cqr(new DOUBLE[4]), 
  Cqi(new DOUBLE[4])
{}

TopAsymmetric::~TopAsymmetric() {
  delete[] Cqr;
  delete[] Cqi;
  delete agmComplement;
}

#define COMMON_INITIALIZATION \
    DOUBLE omega1;            \
    DOUBLE omega2;            \
    DOUBLE omega3;            \
    DOUBLE Lx = Ix*omega.x;   \
    DOUBLE Ly = Iy*omega.y;   \
    DOUBLE Lz = Iz*omega.z;            \
    DOUBLE L2 = sqr(Lx)+sqr(Ly)+sqr(Lz);     \
    E2 = Lx*omega.x + Ly*omega.y + Lz*omega.z;  \
    L  = sqrt(L2);                                              \
    if ( (E2 > L2/Iy and Ix < Iz) or (E2 < L2/Iy and Ix > Iz) ) {       \
        orderFlag = true;                                               \
        I1 = Iz;                                                        \
        I2 = Iy;                                                        \
        I3 = Ix;                                                        \
        omega1 = omega.z;                                               \
        omega2 = -omega.y;                                              \
        omega3 = omega.x;                                               \
        DOUBLE Lp = hypot(Ly, Lz);                                      \
        DOUBLE LpOverL = Lp/L;                                          \
        DOUBLE LxOverL = Lx/L;                                          \
        refine(LpOverL, LxOverL);                                       \
        DOUBLE LyOverLp = Ly/Lp;                                        \
        DOUBLE LzOverLp = Lz/Lp;                                        \
        refine(LyOverLp, LzOverLp);                                     \
        B.xx = -LpOverL; B.xy = LyOverLp*LxOverL; B.xz = LxOverL*LzOverLp; \
        B.yx = 0;        B.yy = -LzOverLp;        B.yz = LyOverLp;      \
        B.zx = LxOverL;  B.zy = LyOverLp*LpOverL; B.zz = LzOverLp*LpOverL; \
    } else {                                                            \
        orderFlag = false;                                              \
        I1 = Ix;                                                        \
        I2 = Iy;                                                        \
        I3 = Iz;                                                        \
        omega1 = omega.x;                                               \
        omega2 = omega.y;                                               \
        omega3 = omega.z;                                               \
        DOUBLE Lp = hypot(Lx, Ly);                                      \
        DOUBLE aoverf = Lx/Lp;                                          \
        DOUBLE boverf = Ly/Lp;                                          \
        refine (aoverf, boverf);                                        \
        DOUBLE coverL = Lz/L;                                           \
        DOUBLE foverL = Lp/L;                                           \
        refine (coverL, foverL);                                        \
        B.xx = aoverf*coverL; B.xy = boverf*coverL; B.xz = -foverL;     \
        B.yx = -boverf;       B.yy = aoverf;        B.yz = 0;           \
        B.zx = aoverf*foverL; B.zy = boverf*foverL; B.zz = coverL;      \
    }                                                                   \
    B *= A;                                                             \
    DOUBLE L12 = L2 - E2*I3;                                            \
    DOUBLE L23 = L2 - E2*I1;                                            \
    DOUBLE I13 = I1 - I3;                                               \
    DOUBLE I23 = I2 - I3;                                               \
    omega1m = copysign(sqrt(L12/I1/I13),omega1);                        \
    omega2m = -copysign(sqrt(omega2*omega2+I1*I13*omega1*omega1/I2/I23),omega1); \
    omega3m = copysign(sqrt(-L23/I3/I13),omega3);                       \
    omegap  = I23*copysign(sqrt(L23/(-I23)/I1/I2/I3), omega3);          \
    m       = L12*(I2-I1)/(L23*I23);                                    \
    snepsilon = omega2/omega2m;                                         \
    cnepsilon = omega1/omega1m;                                         \
    dnepsilon = omega3/omega3m;                                         \
    agm->precompute(m);

void TopAsymmetric::Initialization(const Vector& omega, const Matrix& A) {
  COMMON_INITIALIZATION
  DOUBLE sinam = omega2/omega2m;
  if (fabs(sinam) > 1e-3) {
    DOUBLE cosam = omega1/omega1m;
    epsilon = agm->F_from_tan(sinam/cosam);
  } else 
    epsilon = agm->F_from_sin(sinam);  
  InitializationLastAngle(omega1, omega2, omega3);
}

void TopRecur::Initialization(const Vector& omega, const Matrix& A){
  COMMON_INITIALIZATION
  InitializationLastAngle(omega1, omega2, omega3);
}

void TopAsymmetric::InitializationLastAngle(DOUBLE omega1, DOUBLE omega2, DOUBLE omega3) {
  // calculation of constants A_1, A_2 & coefficients theta function series:  
  agmComplement->precompute(1-m);

  DOUBLE Kp = agmComplement->F();     // complementary quarterperiod
  DOUBLE F  = agmComplement->F_from_sin(I3*omega3m/L);// F(I3 omega3m/L | 1-m)

  relfreq = M_PI_2/agm->F();       // freq of angular velocity (rel. to omegap)

  DOUBLE q = exp(-2.0*relfreq*Kp);      // elliptic nome

  // calculate exp{ pi/2K [K' - F(I3 omega3m/L | 1-m)]}                  
  DOUBLE sqrtxi = exp(relfreq*(copysign(Kp, omega3m)-F));
  DOUBLE xi     = sqrtxi*sqrtxi;

  A2 = L/I1+relfreq*omegap*(xi+1)/(xi-1);      // first term in A2 series.
  DOUBLE q2 = q*q;
  DOUBLE q2n = 1.0;                            // will be q^{2n}           
  DOUBLE xin = 1.0;                            // will be xi^n
  for (int n = 1; n < 10000; n++) {            // no more than 10000 terms
    q2n *= q2;                               // update recursively
    xin *= xi;
    // compute the next term in A2
    DOUBLE dA2 = -2.0*relfreq*omegap*q2n/(1.0-q2n)*(xin-1.0/xin); 
    // add to series
    A2 += dA2;
    // stop upon convergence
    if (fabs(dA2/A2) <= std::numeric_limits<DOUBLE>::epsilon())
        break;
  }
   
  // determine upper bound on number of terms needed in theta function series
  numTerms = (int)(log(std::numeric_limits<DOUBLE>::epsilon())/log(q)+.5);
  // allocate enough memory for the coefficients in these series
  if (numTerms > maxNumTerms) {
    delete[] Cqi;
    delete[] Cqr;
    Cqr = new DOUBLE[numTerms+1];
    Cqi = new DOUBLE[numTerms+1];
    maxNumTerms = numTerms;
  }

  DOUBLE a = 1.0;                              // a will be q^2n.
  DOUBLE b = 1.0;                              // b will be (-1)^n q^(n(n+1)). 
  DOUBLE e = sqrtxi;
  DOUBLE s, c;
  Cqr[0] =  (e+1.0/e);                         // zeroth term in the series
  Cqi[0] = -(e-1.0/e);                         // for real and imag part.
  DOUBLE u = relfreq*epsilon;
  mysincos(u, s, c);                           
  DOUBLE g = 2.0*c*s;                    // = sin(2u).
  DOUBLE h = 2.0*c*c-1.0;                // = cos(2u). 
  DOUBLE r = Cqr[0]*s;                         // Re(theta_1), zeroth term 
  DOUBLE i = Cqi[0]*c;                         // Im(theta_1), zeroth term
  for (int n = 1; n <= numTerms; n++) {
    e *= xi;                                 // e = xi^(n+1/2). 
    a *= q2;                                 // update a and b recursively
    b *= -a;
    Cqr[n] = b*(e+1.0/e);                    // compute next coef Re(theta1)
    Cqi[n] = -b*(e-1.0/e);                   // compute next coef Im(theta1)
    DOUBLE d = s;                      // compute sin&cos recursively.
    s = h*s+g*c;
    c = h*c-g*d;
    r += Cqr[n]*s;                          // add next term Re(theta_1).
    i += Cqi[n]*c;                          // add next term Im(theta_1).
    if ( (fabs(Cqr[n]) < std::numeric_limits<DOUBLE>::epsilon())
         and (fabs(Cqi[n]) < std::numeric_limits<DOUBLE>::epsilon()) )
        numTerms = n-1;                     // if converged, adjust numTerms
  }
  A1 = atan2(i,r);                            // compute arg(r+iI).
}

#define COMMON_EVOLUTION \
    DOUBLE omega1, omega2, omega3;              \
    DOUBLE snwpt, cnwpt, dnwpt;                 \
    agm->sncndn(omegap*t, snwpt, cnwpt, dnwpt); \
    DOUBLE snsn = snwpt*snepsilon;              \
    DOUBLE cncn = cnwpt*cnepsilon;              \
    DOUBLE dndn = dnwpt*dnepsilon;              \
    DOUBLE den = 1./(1-m*snsn*snsn);            \
    omega1 = omega1m*(cncn-snsn*dndn)*den;                              \
    omega2 = omega2m*(snwpt*cnepsilon*dnepsilon+cnwpt*dnwpt*snepsilon)*den; \
    omega3 = omega3m*(dndn-m*snsn*cncn)*den;                            \
    refineBits(omega1, omega2, omega3, I1, I2, I3, E2);                 \
    if (orderFlag) {                                                    \
        omega.x = omega3;                                               \
        omega.y = -omega2;                                              \
        omega.z = omega1;                                               \
        DOUBLE o1 = I1*omega1;                                          \
        DOUBLE o2 = I2*omega2;                                          \
        DOUBLE o3 = I3/L*omega3;                                        \
        DOUBLE f = hypot(o1, o2);                                       \
        o1 /= f;                                                        \
        o2 /= f;                                                        \
        refine (o1, o2);                                                \
        f /= L;                                                         \
        refine (o3, f);                                                 \
        A.xx = -f;     A.xy = 0;   A.xz = o3;                           \
        A.yx = -o2*o3; A.yy = -o1; A.yz = -o2*f;                        \
        A.zx = o1*o3;  A.zy = -o2; A.zz = o1*f;                         \
    } else {                                                            \
        omega.x = omega1;                                               \
        omega.y = omega2;                                               \
        omega.z = omega3;                                               \
        DOUBLE o1 = omega1*I1;                                          \
        DOUBLE o2 = omega2*I2;                                          \
        DOUBLE o3 = omega3*I3/L;                                        \
        DOUBLE f = hypot(o1, o2);                                       \
        o1 /= f;                                                        \
        o2 /= f;                                                        \
        refine (o1, o2);                                                \
        f /= L;                                                         \
        refine (o3, f);                                                 \
        A.xx = o1*o3; A.xy = -o2; A.xz = o1*f;                          \
        A.yx = o2*o3; A.yy = o1;  A.yz = o2*f;                          \
        A.zx = -f;    A.zy = 0;   A.zz = o3;                            \
    }                                                                   \
    DOUBLE s, c;                                                        \
    EvolutionLastAngle(omega1, omega2, omega3, t, s, c);                \
    DOUBLE d;                                                           \
    d = A.xx;                                                           \
    A.xx *= c;                                                          \
    A.xx -= A.xy*s;                                                     \
    A.xy *= c;                                                          \
    A.xy += d*s;                                                        \
    d = A.yx;                                                           \
    A.yx *= c;                                                          \
    A.yx -= A.yy*s;                                                     \
    A.yy *= c;                                                          \
    A.yy += d*s;                                                        \
    d = A.zx;                                                           \
    A.zx *= c;                                                          \
    A.zx -= A.zy*s;                                                     \
    A.zy *= c;                                                          \
    A.zy += d*s;                                                        \
    A *= B;

void TopAsymmetric::Evolution(DOUBLE t, Vector& omega, Matrix& A) {
  COMMON_EVOLUTION
}

void TopRecur::Evolution(DOUBLE t, Vector& omega, Matrix& A) {
  COMMON_EVOLUTION
}

void TopAsymmetric::EvolutionLastAngle(DOUBLE omega1, DOUBLE omega2, DOUBLE omega3,
                                       DOUBLE t, DOUBLE & s, DOUBLE & c)
{  
  DOUBLE u = omegap*t+epsilon;      // compute argument ell fncts   

  mysincos(relfreq*u, s, c);              

  DOUBLE g = 2.0*c*s;               // g= sin(2x), h= cos(2x), used
  DOUBLE h = 2.0*c*c-1.0;           // in recursion.                
  DOUBLE r = Cqr[0]*s;                    // zeroth term Re(theta1)       
  DOUBLE i = Cqi[0]*c;                    // zeroth term Im(theta1)       
  
  for (int n = 1; n <= numTerms; n++) {   // compute series
    DOUBLE d = s;
    s = h*s+g*c;                         // computes sin((2n+1)x) and    
    c = h*c-g*d;                         // cos((2n+1)x) recursively.    
    r += Cqr[n]*s;                       // next term in series for      
    i += Cqi[n]*c;                       // r= Re(theta1) & i= Im(theta1). 
  }
  
  mysincos(A1+A2*t, s, c);
  
  // use addition formula to compute s = sin(psi) and c = cos(psi)
  DOUBLE d = s;                           // dummy variable
  s = s*r-c*i;                                           
  c = c*r+d*i;                            
  d = hypot(r,i);                         // where psi= A1+A2 t+arg(r+iI)  
  s /= d;
  c /= d;
  refine (s, c);
}

TopRecur::TopRecur(DOUBLE _Ix, DOUBLE _Iy, DOUBLE _Iz) :  
  TopAsymmetricCommonPart(_Ix, _Iy, _Iz) 
{
  term  = (DOUBLE *)  malloc( sizeof(DOUBLE)  * maxNTerms );
  coeff = (DOUBLE **) malloc( sizeof(DOUBLE*) * maxNTerms );

  coeff[0] = new DOUBLE[1];
  coeff[0][0] = 1/2.;

  for (int a = 1; a < maxNTerms; a++) {
    coeff[a] = new DOUBLE[a+1];
    coeff[a][0] = 0.5;
    coeff[a][a] = coeff[a-1][a-1]/(4.0*a+2.0);
    DOUBLE den = 1.0/((2.0*a+2.0)*(2.0*a+1.0));
    for (int b = 1; b < a; b++)
        coeff[a][b] = (2.0*a-b+1.0)*den*(coeff[a-1][b-1]+(2.0*a-b)*coeff[a-1][b]);
  }
  for (int order=0; order<=1; order++) {
    q[order] = (DOUBLE***)malloc(maxNTerms*sizeof(DOUBLE**));
    for (int j=0;j<maxNTerms;j++) {
      if (odd(j)) {
        q[order][j] = (DOUBLE**)malloc(j*sizeof(DOUBLE*));
        for (int k=0; k<j; k++)
          q[order][j][k] = (DOUBLE*)calloc(j, sizeof(DOUBLE)); // coefficients set to zero initially
      } else { // even j:
        q[order][j] = (DOUBLE**)malloc((j+2)*sizeof(DOUBLE*));
        for (int k=0; k <= j+1; k++)
          q[order][j][k] = (DOUBLE*)calloc(j+1, sizeof(DOUBLE)); // coefficients set to zero initially
      }
    }
  }

  DOUBLE a = Iz/Ix - 1.0;
  DOUBLE b = Iz/Iy - 1.0;
  X[0] = a*b;       // Note: [0] = Jacobi ordering satisfied, [1] = reversed.
  Y[0] = a+b;

  a = Ix/Iz - 1.0;
  b = Ix/Iy - 1.0;
  X[1] = a*b; 
  Y[1] = a+b;

  // set up the recursion up the maxNTerms:
  for (int order=0; order<=1; order++) {
    q[order][0][0][0] = 1.0;
    // go to maxNTerms levels:
    for (int j=1;j<maxNTerms;j++) {
      if (odd(j)) {
        for (int k=0; k<j; k++)
          for (int l=0; l<j; l++)
              q[order][j][k][l] = 2*(j-k)*q[order][j-1][k][l];
      } else { // even j:
        for (int k=0; k<= j+1; k++) {
          int x = 2*(j-k);
          for (int l=0; l<= j; l++) {
            if (l>1 && l<=j) {
              if (k>0 && k<j)
                q[order][j][k][l] += q[order][j-1][k-1][l-2]*(x+1);
              if (k<(j-1))
                q[order][j][k][l] -= q[order][j-1][k][l-2]*x;
            }   
            if (l>0 && l<j) {
              if (k>1 && k<=j)
                q[order][j][k][l] += Y[order]*q[order][j-1][k-2][l-1]*(x+2);
              if (k>0 && k<j)
                q[order][j][k][l] -= Y[order]*q[order][j-1][k-1][l-1]*(x+1);
            }
            if (l<(j-1)) {
              if (k>2 && k<(j+2))
                q[order][j][k][l] += X[order]*q[order][j-1][k-3][l]*(x+3);
              if (k>1 && k<=j)
                q[order][j][k][l] -= X[order]*q[order][j-1][k-2][l]*(x+2);
            }
          }                                                      
        }
      }
    }
  }
}

TopRecur::~TopRecur() {
  for (int a=0;a<maxNTerms;a++)
      delete[] coeff[a];
  free(coeff);
  free(term);
  for (int order=0;order<=1;order++) {
    for (int j=0;j<maxNTerms;j++) {
      if (odd(j))
        for (int k=0;k<j;k++)
          free(q[order][j][k]);
      else
        for (int k=0;k<=(j+1);k++)
          free(q[order][j][k]);
      free(q[order][j]);
    }
    free(q[order]);
  }
}

static DOUBLE invwi;
static DOUBLE invwf;
static DOUBLE invwi2;
static DOUBLE invwf2;
static DOUBLE scaledL1L2L3i;
static DOUBLE scaledL1L2L3f;
static DOUBLE eps;
static DOUBLE LoverI3;
static DOUBLE minusepsilonLoverI3nplusone;
static DOUBLE epspower[maxNTerms];

// jth derivative of Omega:
void TopRecur::firstOmega(DOUBLE omega1i, DOUBLE omega2i, DOUBLE omega3i,
			  DOUBLE omega1f, DOUBLE omega2f, DOUBLE omega3f,
			  DOUBLE &Omegai, DOUBLE &Omegaf)
{
  DOUBLE L1 = I1*omega1i;
  DOUBLE L2 = I2*omega2i;
  DOUBLE L3 = I3*omega3i;
  DOUBLE Lsqr = L*L;
  eps = 1.0 - I3*E2/Lsqr;
  minusepsilonLoverI3nplusone = -eps;
  LoverI3 = L/I3;
  invwi = Lsqr/(Lsqr-L3*L3);
  invwi2 = invwi*invwi;
  scaledL1L2L3i = L1*L2*L3/Lsqr/L*I3*(I1-I2)/(I1*I2);
  L1 = I1*omega1f;
  L2 = I2*omega2f;
  L3 = I3*omega3f;     
  invwf = Lsqr/(Lsqr-L3*L3);
  invwf2 = invwf*invwf;
  scaledL1L2L3f = L1*L2*L3/Lsqr/L*I3*(I1-I2)/(I1*I2);
  // construct powers of eps also recursively
  epspower[0] = 1.0;

  Omegai = q[orderFlag][0][0][0]*invwi;
  Omegaf = q[orderFlag][0][0][0]*invwf;

  minusepsilonLoverI3nplusone *= LoverI3;
  Omegai *= minusepsilonLoverI3nplusone;
  Omegaf *= minusepsilonLoverI3nplusone;    

  Omegai += LoverI3;
  Omegaf += LoverI3;
}

// jth derivative of Omega:
inline void TopRecur::nextOmega(int j, DOUBLE & Omegai, DOUBLE & Omegaf) {

  if (j == maxNTerms) {
    // cannot reach convergence before set number of recursions
    Omegai = Omegaf = 0;  // for lack of anything better
  }

  epspower[j] = eps*epspower[j-1];
  Omegai      = 0.0;
  Omegaf      = 0.0;

  DOUBLE** qj = q[orderFlag][j];

  if ( j & 1 ) { 
    for (int k = 0; k < j; k++) {
      int     jminusone = j - 1;
      DOUBLE* qjk       = qj[k];
      int     largestl  = jminusone;
      int     dl        = k - j/2;
      if (dl > 0) 
          largestl -= dl;
      DOUBLE ck  = qjk[largestl];
      int smallestl = jminusone - k;
      while (largestl > smallestl) {
          ck *= eps;
          ck += qjk[--largestl];
      } 
      ck    *= epspower[smallestl];
      Omegai = Omegai*invwi+ck;
      Omegaf = Omegaf*invwf+ck;
    }   
    minusepsilonLoverI3nplusone *= LoverI3;
    Omegai *= minusepsilonLoverI3nplusone*invwi2*scaledL1L2L3i;
    Omegaf *= minusepsilonLoverI3nplusone*invwf2*scaledL1L2L3f;  
  } else {
    for (int k = 0; k < (j+2); k++) {        
      DOUBLE* qjk      = qj[k];
      int     largestl = j;
      int     dl       = k - j/2;     
      if (dl > 0)  
        largestl -= dl;
      if (largestl < 0) 
        largestl = 0;
      int smallestl = j - k;
      if (smallestl < 0) 
        smallestl = 0;
      DOUBLE ck = qjk[largestl];
      while (largestl > smallestl) {
        ck *= eps;
        ck += qjk[--largestl];
      }
      ck    *= epspower[smallestl];
      Omegai = Omegai*invwi+ck;
      Omegaf = Omegaf*invwf+ck;
    }

    minusepsilonLoverI3nplusone *= LoverI3;
    Omegai *= minusepsilonLoverI3nplusone;
    Omegaf *= minusepsilonLoverI3nplusone;    
  }
}

void TopRecur::InitializationLastAngle(DOUBLE omega1, DOUBLE omega2, DOUBLE omega3) {
  initialomega1 = omega1;
  initialomega2 = omega2;
  initialomega3 = omega3;
}

void TopRecur::EvolutionLastAngle(DOUBLE omega1, DOUBLE omega2, DOUBLE omega3,
                                  DOUBLE t, DOUBLE &s, DOUBLE &c) 
{
  DOUBLE finalomega1 = omega1;
  DOUBLE finalomega2 = omega2;
  DOUBLE finalomega3 = omega3;
  static DOUBLE oldt;
  static DOUBLE tpower[maxNTerms]; // tpower[n] = t^(n-1)
  if (t!=oldt) {
    oldt = t;
    tpower[0] = t;
    for (int i = 1; i<maxNTerms; i++)
      tpower[i] = tpower[i-1]*t;
  }
  int sign = +1;
  int n = 0;
  // first approximation:
  DOUBLE Omegai;
  DOUBLE Omegaf;
  firstOmega(initialomega1, initialomega2, initialomega3, 
             finalomega1, finalomega2, finalomega3,
             Omegai, Omegaf);
  term[n] = tpower[n]*(Omegai + Omegaf);
  DOUBLE psi0 = coeff[0][0]*term[0];
  DOUBLE dpsi, dpsi0;
  // higher order approximations:
  static const int nInitialTerms = 3;
  while (n < nInitialTerms) {
    n++;
    sign *= -1;
    nextOmega(n, Omegai, Omegaf);
    term[n] = tpower[n]*(Omegai+sign*Omegaf);
  }
  dpsi = 0;
  for (int k = 1; k <= n; k++) 
    dpsi += coeff[n][k]*term[k];
  do {
    n++;
    sign *= -1;
    nextOmega(n, Omegai, Omegaf);
    term[n] = tpower[n]*(Omegai+sign*Omegaf);
    dpsi0 = dpsi;
    dpsi = 0.0;
    for (int k = 1; k <= n; k++) 
      dpsi += coeff[n][k]*term[k];
  } while ( (n < (maxNTerms-1)) and 
            (fabs((dpsi-dpsi0)/psi0) > std::numeric_limits<DOUBLE>::epsilon()) );
  mysincos(psi0 + dpsi, s, c);
}

/* IMPLEMENTATION OF THE AGM SCALE */
AgmScale::AgmScale() :
  m(-1.0)
{}

static const DOUBLE ACC = sqrt(std::numeric_limits<DOUBLE>::epsilon());

void AgmScale::precompute(DOUBLE m_) 
{
  if (m_ < 0.0) {
      fprintf(stderr, "ERROR: m < 0.0 in AgmScale (m = %f): FIX IT!\n", m_);
      exit(1);
  }
  if (m != m_) {
    m    = m_;
    N    = 0;
    a[N] = 1.0;
    b[N] = sqrt(1.0 - m);
    d[N] = b[N]/a[N];   
    while (N < MAXN) {
      N++;
      a[N] = (a[N-1] + b[N-1])/2;
      b[N] = sqrt(a[N-1] * b[N-1]);
      d[N] = b[N]/a[N];     
      if (fabs(a[N-1]-b[N-1]) <= 2*ACC*a[N]) 
        break;
    }  
    K      = M_PI_2/a[N];
    Kscale = 1.0/((1<<N)*a[N]);
  }
}

DOUBLE AgmScale::F() const {
  return K;
};

DOUBLE AgmScale::F_from_tan(DOUBLE t) const {
  int s = 0;
  for (int i = 0; i < N; i++) {
    DOUBLE c = 1.0 - d[i]*t*t;
    s *= 2;
    if (c < 0.0) s += (t>0.0?1:-1);
    t = (1.0 + d[i])*t/c;
  }
  return Kscale*(atan(t)+s*M_PI);
};

DOUBLE AgmScale::F_from_sin(DOUBLE sinphi) const {
  if (sinphi <= -1.0)   // Better fix may be needed 
    return -F();
  else if (sinphi >= 1.0) 
    return F();
  else
    return F_from_tan(sinphi/sqrt(1-sinphi*sinphi));
}

DOUBLE AgmScale::F(DOUBLE phi) const {
  return F_from_tan(tan(phi));
}

void AgmScale::sncndn(DOUBLE x, DOUBLE &sn, DOUBLE &cn, DOUBLE &dn) const {
  if (m != 0.0) {
    mysincos(a[N]*x, sn, cn);
    dn = 1.0;
    if (sn) {
      DOUBLE d = cn/sn;
      DOUBLE c = a[N]*d;
      for (int i = N; i >= 0 ; i--) {
        d *= c;
        c *= dn;
        dn = (b[i] + d)/(a[i] + d);
        d  = c/a[i];
      }
      sn = copysign(1.0/sqrt(c*c + 1.0),sn);
      cn = c*sn;
    }
  } else {
    sn = tanh(x);
    cn = dn = 1.0/cosh(x);
  }
}
