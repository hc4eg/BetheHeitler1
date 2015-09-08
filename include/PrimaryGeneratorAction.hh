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
  //inline G4double Get_central_energy() { return central_energy;};
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

  OutputFile * pOutputFile;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif


