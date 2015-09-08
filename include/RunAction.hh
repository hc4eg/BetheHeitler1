//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifndef RunAction_h
#define RunAction_h 1

#include "G4UserRunAction.hh"
#include "globals.hh"
#include "G4SystemOfUnits.hh"
#include <vector>
#include <fstream>
#include <iostream>
using namespace std;
#include "G4ThreeVector.hh"
#include "OutputFile.hh"
#include "PrimaryGeneratorAction.hh"
#include "DetectorConstruction.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class G4Run;

class RunAction : public G4UserRunAction
{
public:
  RunAction( PrimaryGeneratorAction * prim, DetectorConstruction * det);
  virtual ~RunAction();

  void BeginOfRunAction(const G4Run*);
  void   EndOfRunAction(const G4Run*);
    
  //void FillPerEvent(G4double E_i, G4double t_i, G4double phi_i, G4ThreeVector r_i, G4double t_f, G4double phi_f, G4ThreeVector r_f); 

private:
  //vector < vector <G4double> > fRunData;
  OutputFile * pOutputFile;
  PrimaryGeneratorAction * pPrimaryGeneratorAction;
  DetectorConstruction * pDetectorConstruction;
 
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif

