//
//Glen Pridham, 2014
//This is a frankenstein program with lots of bits and pieces from other places, mainly Geant4 examples.  There are; consequently, numerous instances of obsolete code floating around which serve no purpose.

#include <iostream>
using namespace std;

#include "G4RunManager.hh"
#include "G4UImanager.hh"

#include "Randomize.hh"

#include "DetectorConstruction.hh"
//#include "PhysicsList.hh"
#include "PrimaryGeneratorAction.hh"
#include "RunAction.hh"
#include "EventAction.hh"
//#include "SteppingAction.hh"
//#include "SteppingVerbose.hh"

//#include "G4EmLivermorePhysics.hh" //valid for photons and electrons from 250 eV - 100 GeV
//G4EmLivermorePolarizedPhysics.hh for polarized reactions
#include "QGSP_BERT_HP.hh" //http://geant4.web.cern.ch/geant4/physics_lists/geant4.4.1/tracker/index.html //I add BERT_HP just in case I get some nuclear photoprocesses of significance

#ifdef G4VIS_USE
#include "G4VisExecutive.hh"
#endif

#ifdef G4UI_USE
#include "G4UIExecutive.hh"
#endif

#define useDefault (0)

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

int main(int argc,char** argv)
{
  // Choose the Random engine
  //
  CLHEP::HepRandom::setTheEngine(new CLHEP::RanecuEngine);
  
  // User Verbose output class
  //
  //G4VSteppingVerbose::SetInstance(new SteppingVerbose);
     
  // Construct the default run manager
  //
  G4RunManager * runManager = new G4RunManager; //note: this can access the world volume

  // Set mandatory initialization classes
  //
  DetectorConstruction* detector = new DetectorConstruction();
  runManager->SetUserInitialization(detector);
  //

  //#if useDefault 
  //PhysicsList* physics = new PhysicsList;
  //#else
  QGSP_BERT_HP * physics = new QGSP_BERT_HP();
  //#endif
  runManager->SetUserInitialization(physics);    

  // Set user action classes
  //
  PrimaryGeneratorAction* gen_action = 
                          new PrimaryGeneratorAction(detector);
  runManager->SetUserAction(gen_action);
  //
  RunAction* run_action = new RunAction(gen_action, detector);  
  runManager->SetUserAction(run_action);
  //
  EventAction* event_action = new EventAction(detector);
  runManager->SetUserAction(event_action);
  //
  //SteppingAction* stepping_action =
  //                  new SteppingAction(detector, event_action);
  //runManager->SetUserAction(stepping_action);
  
  // Initialize G4 kernel
  //
  runManager->Initialize();
  
#ifdef G4VIS_USE
  // Initialize visualization
  G4VisManager* visManager = new G4VisExecutive;
  // G4VisExecutive can take a verbosity argument - see /vis/verbose guidance.
  // G4VisManager* visManager = new G4VisExecutive("Quiet");
  visManager->Initialize();
#endif

  // Get the pointer to the User Interface manager
  G4UImanager* UImanager = G4UImanager::GetUIpointer();

  if (argc!=1)   // batch mode
    {
      G4String command = "/control/execute ";
      G4String fileName = argv[1];
      UImanager->ApplyCommand(command+fileName);
    }
  else
    {  // interactive mode : define UI session
#ifdef G4UI_USE
      G4UIExecutive* ui = new G4UIExecutive(argc, argv);
#ifdef G4VIS_USE
      UImanager->ApplyCommand("/control/execute vis.mac"); 
#endif
      if (ui->IsGUI())
	UImanager->ApplyCommand("/control/execute gui.mac");
      ui->SessionStart();
      delete ui;
#endif
    }

  // Job termination
  // Free the store: user actions, physics_list and detector_description are
  //                 owned and deleted by the run manager, so they should not
  //                 be deleted in the main() program !
#ifdef G4VIS_USE
  delete visManager;
#endif
  delete runManager;

  return 0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....
