#include "MagneticField.hh"
using namespace std;

MagneticField::MagneticField()
	{
	// set the defaults
	B0 = -0.4536*tesla;
	EFB = 22.26*cm;
	D = 5.08*cm;
	C[0] = 0.35;
	C[1] = 2.71;
	C[2] = -0.94;
	C[3] = 0.27;
	C[4] = -0.027;
	C[5] = 0.0009;
	R_max = EFB + 4.*D;
	R_min = EFB - 3.*D;
	if(R_min < 0.)R_min = 0.;
	DG = 0.3*D;
	}
MagneticField::~MagneticField()
	{}

void
MagneticField::GetFieldValue(const G4double point[4],
				      G4double *Bfield ) const
	{
	G4double x,y,z;
	x = point[0]; y = point[1]; z = point[2];
	// note geant coordinate z up.
	G4double R = sqrt(x*x + y*y);
	if(R > R_max)
		{
		// field is zero
		Bfield[0] = 0;
		Bfield[1] = 0; 
		Bfield[2] = 0;
		return;
		}
	if(R < R_min)
		{
		// field is central value
		Bfield[0] = 0;
		Bfield[1] = 0; 
		Bfield[2] = B0;
		return;
		}
	// first calculate field at (x,y,z = 0)
	G4double B00 = CalFieldZ(x, y);
	if(z == 0.) 
		{
		Bfield[0] = 0;
		Bfield[1] = 0; 
		Bfield[2] = B00;
		return;
		}
	// else we need to estimate components of field
	// off the z = 0 plane
	// first calculate the fields in the 13 point grid
	// see RAYTRACE manual
	G4double B10   = CalFieldZ( x+DG,    y       );
	G4double B_10  = CalFieldZ( x-DG,    y       );
	G4double B20   = CalFieldZ( x+2.*DG, y       );
	G4double B_20  = CalFieldZ( x-2.*DG, y       );
	G4double B01   = CalFieldZ( x,       y+DG    );
	G4double B11   = CalFieldZ( x+DG,    y+DG    );
	G4double B_11  = CalFieldZ( x-DG,    y+DG    );
	G4double B02   = CalFieldZ( x,       y+2.*DG );
	G4double B0_1  = CalFieldZ( x,       y-DG    );
	G4double B1_1  = CalFieldZ( x+DG,    y-DG    );
	G4double B_1_1 = CalFieldZ( x-DG,    y-DG    );
	G4double B0_2  = CalFieldZ( x,       y-2.*DG );
	// now calculate the field components
	G4double Bz = B00;
	G4double s1 = B10 + B_10 + B01 + B0_1 - 4.*B00;
	G4double s2 = B20 + B_20 + B02 + B0_2 - 4.*B00;
	Bz = Bz - z*z/(DG*DG) * (2./3.*s1 - 1./24.*s2);
	G4double s3 = B11 + B_11 + B1_1 + B_1_1 
			-2.*B10 - 2.*B_10 - 2.*B01 -2.*B0_1 + 4.*B00;
	Bz = Bz + z*z*z*z/(DG*DG*DG*DG) * (-1./6.*s1 + 1./24.*s2 + 1./12.*s3);

	s1 = B10 - B_10;
	s2 = B20 - B_20;
	G4double Bx = z/DG*(2./3.*s1 - 1./12.*s2);
	s3 = B11 + B1_1 - B_11 - B_1_1 - 2.*B10 + 2.*B_10;
	Bx = Bx + z*z*z/(DG*DG*DG)*( 1./6.*s1 - 1./12.*s2 - 1./12.*s3);

	s1 = B01 - B0_1;
	s2 = B02 - B0_2;
	G4double By = z/DG*(2./3.*s1 - 1./12.*s2);
	s3 = B11 + B_11 - B1_1 - B_1_1 - 2.*B01 + 2.*B0_1;
	By = By + z*z*z/(DG*DG*DG)*( 1./6.*s1 - 1./12.*s2 - 1./12.*s3);

	Bfield[0] = Bx;
	Bfield[1] = By; 
	Bfield[2] = Bz;
	return;
	}

G4double MagneticField::CalFieldZ( G4double x, G4double y) const
	{
	// private member function to calculate the field
	// at (x, y, 0) Field is only Bz
	G4double R = sqrt(x*x + y*y);
	G4double s0 = (R - EFB)/D;
	G4double ss = 1.;
	G4double S = C[0];
	for(int i = 1; i < 6; i++)
		{
		ss *= s0;
		S += C[i]*ss;
		}
	G4double Bz = B0/(1. + exp(S));
	return(Bz);
	}
