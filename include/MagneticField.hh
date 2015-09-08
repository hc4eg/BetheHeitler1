//  
//    *************************************
//    *                                   *
//    *    MagneticField.hh               *
//    *                                   *
//    *************************************

#ifndef MagneticField_h
#define Magneticfield_h 1

#include "globals.hh"
#include "G4SystemOfUnits.hh"
#include "G4MagneticField.hh"
#include "G4ios.hh"

#include <fstream>
#include <vector>
#include <cmath>

using namespace std;


class MagneticField : public G4MagneticField
{
public:
  MagneticField( );
  ~MagneticField( );
  void  GetFieldValue( const  G4double Point[4],
		       G4double *Bfield          ) const;
  inline void SetFieldB0( G4double field) {B0 = field;};
  inline G4double GetFieldB0() {return B0;};
  
private:
  G4double B0; // central field
  G4double D; // gap
  G4double C[6]; // coefficients defining the fringe field
  G4double EFB; // radius of EFB of circular dipole magnet
	// dipole centered at (0, 0, 0)
	// central field is (0, 0, B0) 
  G4double R_max; // maximum distance from center for finite field.
  G4double R_min; // minimum distance from center for constant field B0.
  G4double DG; // parameter for grid for calculating off axis fields.

  G4double CalFieldZ(G4double, G4double) const;
};

#endif
