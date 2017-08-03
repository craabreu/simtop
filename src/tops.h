/*
   tops.h
  
   Header file for the top classes for free rotation of general tops.
 
   Ramses van Zon
   May  9, 2007 (v1)
   May 29, 2009 (v2)
 
   From an abstract class "Top", the following classes are derived:
                   
   TopSpherical : for spherical tops,  for which Ix = Iy = Iz. See [1,2].
   TopProlate   : for prolate tops,    for which Ix < Iy = Iz. See [1,2].
   TopOblate    : for spherical tops,  for which Ix = Iy < Iz. See [1,2].
   TopAsymmetric: for asymmetric tops, for which Ix < Iy < Iz. See [2].
   TopRecur     : also for asymmetric tops, but computed using
                  a recursive scheme. See [3].
 
   [1] Lisandro Hernandez de la Pena, Ramses van Zon, Jeremy Schofield
        & Sheldon B. Opps, "Discontinuous molecular dynamics for
        semiflexible and rigid bodies" Journal of Chemical Physics
        126, 074105 (2007).
   [2] Ramses van Zon and Jeremy Schofield, "Numerical implementation
        of the exact dynamics of free rigid bodies" Journal of
        Computational Physics, 223, 145 (2007).
   [3] Ramses van Zon and Jeremy Schofield, "Symplectic algorithms for
        simulations of rigid-body systems using the exact solution of
        free motion", Physical Review E 75, 056701 (2007).
*/

#ifndef _TOPSH_
#define _TOPSH_

#include "vecmat3.h"

/*
 * Abstract base class Top, derivatives of which compute torque-free rotation.
 */
class Top 
{
 public:
  // To initialize.
  virtual void Initialization( const Vector & omega, const Matrix & A ) = 0;
  
  // To compute the angular velocity vector and the orientation matrix at time t.
  virtual void Evolution( DOUBLE t, Vector & omega, Matrix & A) = 0;
  
  // To propagate over a given time interval.
  void Propagation( DOUBLE dt, Vector & omega, Matrix & A ) {
    Initialization(omega, A);
    Evolution(dt, omega, A);
  }  

  // To deconstruct.
  virtual ~Top() {}
};

/*
 * Class to compute the torque-free rotation of a spherical top.
 */
class TopSpherical: public Top 
{
 public:
  // To set the moment of inertia. 
  explicit TopSpherical( DOUBLE I );

  // To initialize the spherical top.
  void Initialization( const Vector & omega, const Matrix & A );

  // To compute the angular velocity vector and the orientation matrix of the spherical top at time t.
  void Evolution( DOUBLE t, Vector & omega, Matrix & A );

 private:
  DOUBLE I;
  Vector omegab;
  Matrix A0;
};

/*
 *  Class to compute the torque-free rotation of a prolate top.
 */
class TopProlate: public Top {
  public:
    // To set the moments of inertia for the prolate top.
    TopProlate( DOUBLE Ix, DOUBLE Iz );

    // To initialize the prolate top.
    void Initialization( const Vector & omega, const Matrix & A );

    // Compute the angular velocity vector and orientation matrix of the prolate top at time t.
    void Evolution( DOUBLE t, Vector & omega, Matrix & A );

  private:
    DOUBLE  I1;
    DOUBLE  I3;
    Vector  omegab0;
    Vector  LboverI1;
    Matrix  A0;
    DOUBLE  omegap;
};

/*
 *  Class to compute the torque-free rotation of an oblate top.
 */
class TopOblate: public Top 
{
 public:
  // To initialize the moments of inertia for the oblate top.
  TopOblate( DOUBLE Ix, DOUBLE Iz );

  // To initialize the oblate top.
  void Initialization( const Vector & omega, const Matrix & A );

  // To compute the angular velocity vector and orientation matrix of the oblate top at time t.
  void Evolution( DOUBLE t, Vector & omega, Matrix & A );

 private:
  DOUBLE  I1;
  DOUBLE  I3;
  Vector  omegab0;
  Vector  LboverI1;
  Matrix  A0;
  DOUBLE  omegap;
};

class AgmScale; // forward declaration so pointers to this helper class can be used.

/* 
 * Parent class to compute torque-free rotation of a general asymmetric top. See [2].
 * This class contains common parts between TopAsymmetric and TopRecur.
 */
class TopAsymmetricCommonPart: public Top 
{
 public:
  // To construct the moments of inertia for the asymmetric top. 
  TopAsymmetricCommonPart( DOUBLE Ix, DOUBLE Iy, DOUBLE Iz );

  // To deconstruct.
  ~TopAsymmetricCommonPart();
    
 protected:
  bool      orderFlag;
  DOUBLE    Ix;
  DOUBLE    Iy;
  DOUBLE    Iz;
  DOUBLE    I1;
  DOUBLE    I2;
  DOUBLE    I3;
  DOUBLE    omega1m;
  DOUBLE    omega2m;
  DOUBLE    omega3m;
  DOUBLE    omegap;
  DOUBLE    snepsilon;
  DOUBLE    cnepsilon;
  DOUBLE    dnepsilon;
  DOUBLE    L;
  DOUBLE    E2;
  DOUBLE    m;
  Matrix    B;
  AgmScale* agm;
};

/* 
 *  Class to compute the torque-free rotation of a general asymmetric top. See [2].
 */
class TopAsymmetric: public TopAsymmetricCommonPart 
{
 public:
  // To construct the asymmetric top with give moments of inertia. 
  TopAsymmetric(DOUBLE Ix, DOUBLE Iy, DOUBLE Iz);

  //  To deconstruct.
  ~TopAsymmetric();

  // To initialization routine for the asymmetric top.
  void Initialization( const Vector & omega, const Matrix & A );

  // To compute the angular velocity vector and orientation matrix of the asymmetric top at time t.
  void Evolution( DOUBLE t, Vector & omega, Matrix & A );

 private:
  DOUBLE    epsilon;
  AgmScale* agmComplement;
  int       numTerms;
  DOUBLE    relfreq;
  DOUBLE    A1;
  DOUBLE    A2;
  int       maxNumTerms;
  DOUBLE*   Cqr;
  DOUBLE*   Cqi;

  // To initialize the computation of the last angle for the asymmetric top.
  void InitializationLastAngle(DOUBLE omega1, DOUBLE omega2, DOUBLE omega3);

  // To compute the last angle needed for the orientation matrix.
  void EvolutionLastAngle(DOUBLE omega1, DOUBLE omega2, DOUBLE omega3,
                            DOUBLE t, DOUBLE & s, DOUBLE & c);
};

/*
 *  Class to compute the torque-free rotation of a general asymmetric
 *  top using a recursive formula for the 'last angle'.  See [3].
 */
class TopRecur: public TopAsymmetricCommonPart 
{
 public:
  //  To construct the recursive top with given moments of inertia. 
  TopRecur( DOUBLE Ix, DOUBLE Iy, DOUBLE Iz );

  //  To deconstruct.
  ~TopRecur();

  // To initialize the asymmetric top.
  void Initialization( const Vector & omega, const Matrix & A );

  // To compute the angular velocity vector and orientation matrix of the asymmetric top at time t.
  void Evolution( DOUBLE t, Vector & omega, Matrix & A );

 private:
  DOUBLE**  coeff;
  DOUBLE*   term;
  DOUBLE    initialomega1;
  DOUBLE    initialomega2;
  DOUBLE    initialomega3;
  DOUBLE    X[2];    // two because of two ordering cases
  DOUBLE    Y[2];
  DOUBLE*** q[2];    // two three-index arrays for the recursion
  

  // To initialize the computation of the last angle for the asymmetric top.
  void InitializationLastAngle(DOUBLE omega1, DOUBLE omega2, DOUBLE omega3);

  // To compute the last angle needed for the orientation matrix.
  void EvolutionLastAngle(DOUBLE omega1, DOUBLE omega2, DOUBLE omega3, 
                            DOUBLE t, DOUBLE & s, DOUBLE & c);  

  /// Internal routines
  void firstOmega(DOUBLE omega1i, DOUBLE omega2i, DOUBLE omega3i,
                  DOUBLE omega1f, DOUBLE omega2f, DOUBLE omega3f,
                  DOUBLE &Omegai, DOUBLE &Omegaf);
  void nextOmega(int j, DOUBLE &Omegai, DOUBLE &Omegaf);
};

#endif
