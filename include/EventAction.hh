#ifndef EventAction_h
#define EventAction_h 1

#include "G4UserEventAction.hh"
#include "globals.hh"
#include "G4SystemOfUnits.hh"
#include <fstream>
#include <iostream>
using namespace std;

#include "G4SDManager.hh"
#include "WireChamberHit.hh"
#include "PaddleHit.hh"
#include "MonitorHit.hh"
#include "G4TwoVector.hh"
#include "OutputFile.hh"
#include "DetectorConstruction.hh"

#define NUMPADDLES 29


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class EventAction : public G4UserEventAction
{
public:
  EventAction(DetectorConstruction*);
  virtual ~EventAction();

  void  BeginOfEventAction(const G4Event*);
  void    EndOfEventAction(const G4Event*);
    
private:
   G4ThreeVector fR_i;
                    
   G4int     printModulo;
                             
   //EventActionMessenger*  eventMessenger;
   OutputFile * pOutputFile;
   G4double fVDCSpacing;
   DetectorConstruction* pDetector;
   G4double fVDCSigma; // std. dev. of position resolution
   G4bool fRequireHodoscopeHit;
   G4double fPaddleThreshold;
   G4double fsigmaLightStatistical;
   G4double fsigmaLightNoise;
   

};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif

    
