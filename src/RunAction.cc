//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "RunAction.hh"
//CLHEP includes
#include <CLHEP/Random/Randomize.h>
#include "G4Run.hh"
#include "G4RunManager.hh"
#include "G4UnitsTable.hh"
#include <iostream>
using namespace std;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

RunAction::RunAction(PrimaryGeneratorAction * prim, DetectorConstruction * det)
{
	// Set up the event-by-event output file characteristics
	// File name is <prefix>.<runno>.dat
	G4String prefix = "BH";
	pOutputFile = OutputFile::GetOutputFilePointer();
	pOutputFile->Set_filename_prefix(prefix);
	pOutputFile->Set_use_unique_filename(false); // filename set by G4 run number
	//pOutputFile->Set_use_unique_filename(true); // filename uses run number from file
	pPrimaryGeneratorAction = prim;
	pDetectorConstruction = det;

	G4int seed = time(NULL);
	CLHEP::HepRandom::setTheSeed(seed);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

RunAction::~RunAction()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunAction::BeginOfRunAction(const G4Run* aRun)
{ 
  G4int G4runno = aRun->GetRunID();
  G4cout << "### Begin of Run Action: G4Run " << G4runno << " start." << G4endl;

  //inform the runManager to save random number seed
  G4RunManager::GetRunManager()->SetRandomNumberStore(true);

  // Open the event-by-event output file.
  pOutputFile->OpenFile(G4runno);

  // Output the run parameters to that file.
  pPrimaryGeneratorAction->PrintParametersToFile();
  G4double thick = pDetectorConstruction->GetTargetThickness();
  pOutputFile->WriteComment("----Detector Parameters----");
  char line[100];
  sprintf(line, "  Target Thickness = %.4g mm", thick/mm);
  pOutputFile->WriteComment( G4String(line));
  G4String mat = "  Target Material = " + pDetectorConstruction->GetTargetMaterialName();
  pOutputFile->WriteComment( mat);
  G4double B0 = pDetectorConstruction->GetMagneticField();
  sprintf(line, "  Magnetic Field = %.4g Tesla", B0/tesla);
  pOutputFile->WriteComment( G4String(line));
  G4cout << "End Begin of Run Action." << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunAction::EndOfRunAction(const G4Run* aRun)
{
  G4int nEvents = aRun->GetNumberOfEvent();
  G4cout << "End Of Run Action: Events = " << nEvents << G4endl;

  // Close the event-by-event output file.
  pOutputFile->CloseFile();
  return;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
