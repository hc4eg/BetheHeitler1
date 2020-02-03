//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifndef PrimaryGeneratorAction_h
#define PrimaryGeneratorAction_h 1

#include "G4VUserPrimaryGeneratorAction.hh"
#include "globals.hh"
#include "G4SystemOfUnits.hh"
#include "OutputFile.hh"
#include "GeneratorMessenger.hh"
#include "DetectorConstruction.hh"
#include "G4Event.hh"
#include "G4ParticleGun.hh"

#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <cmath>

class GeneratorMessenger;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class PrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction
{
public:
  PrimaryGeneratorAction(DetectorConstruction*);    
  virtual ~PrimaryGeneratorAction();

  void GeneratePrimaries(G4Event*);

  inline void Set_central_energy(G4double val) { central_energy = val;};
  inline void Set_delta_min(G4double val) { delta_min = val;};
  inline void Set_delta_max(G4double val) { delta_max = val;};
  inline void Set_radius_max(G4double val) { radius_max = val;};
  inline void Set_x_min(G4double val) { x_min = val;};
  inline void Set_x_max(G4double val) { x_max = val;};
  inline void Set_y_min(G4double val) { y_min = val;};
  inline void Set_y_max(G4double val) { y_max = val;};
  inline void Set_angle_max(G4double val) { angle_max = val;};
  inline void Set_theta_min(G4double val) { theta_min = val;};
  inline void Set_theta_max(G4double val) { theta_max = val;};
  inline void Set_phi_min(G4double val) { phi_min = val;};
  inline void Set_phi_max(G4double val) { phi_max = val;};
  inline void Set_gamma_mode(G4bool mode) { gamma_mode = mode; }
  inline void Set_pair_mode(G4bool mode) { pair_mode = mode;}
  //inline G4double Get_central_energy() { return central_energy;};

  // New functions declaration that read BH_cross_sections generated 
  // pair data into primary vertex and run G4 simulation
  void OpenFile(string);
  int Convert(int linenumber);
  int ConvertNext();
  inline void Clear(){
	Ee = 0; Ep = 0;
	KEe = 0; KEp = 0;
	Pe = 0; Pp = 0;
	Thetae = 0; Thetap = 0;
	Phie = 0; Phip = 0;
	Pex = 0; Pey = 0; Pez = 0;
	Ppx = 0; Ppy = 0; Ppz = 0;
  }
  inline G4double Get_E(int n) { if(n == 0) return Ep; else return Ee; };
  inline G4double Get_KE(int n) { if(n == 0) return KEp; else return KEe; };
  inline G4double Get_P(int n) { if(n == 0) return Pp; else return Pe; };
  inline G4double Get_P_com( int n, int m) {
	if( n == 0 && m == 0) return Ppx;
	else if( n == 0 && m == 1) return Ppy;
	else if( n == 0 && m == 2) return Ppz;
	else if( n == 1 && m == 0) return Pex;
	else if( n == 1 && m == 1) return Pey;
	else if( n == 1 && m == 2) return Pez;
	else { 
		cerr << "Wrong function parameters, return 0." << endl;
		return 0;
	}
  };
  inline G4double Get_Theta(int n) { if(n == 0) return Thetap; else return Thetae; };
  inline G4double Get_Phi(int n) { if(n == 0) return Phip; else return Phie; };

  void PrintVertex();
  void PrintParameters();
  void PrintParametersToFile();

private:
  G4ParticleGun*           particleGun;	 //pointer a to G4  class
  DetectorConstruction*    Detector;     //pointer to the geometry
  GeneratorMessenger* gunMessenger;   //messenger of this class

  G4double delta_min, delta_max;
  G4double x_min, x_max, y_min, y_max, radius_max;
  G4double angle_max, theta_min, theta_max, phi_min, phi_max;
  G4double central_energy, target_position;
  G4double target_thickness;
  G4bool gamma_mode;

  G4bool pair_mode;
  G4double E_gamma,c;
  G4double Ee,KEe, Pe,Thetae, Phie, Pex, Pey, Pez, Me;
  G4double Ep,KEp, Pp,Thetap, Phip, Ppx, Ppy, Ppz;
  // new members related to reading data file
  int maxline;
  ifstream infile;

  OutputFile * pOutputFile;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif


