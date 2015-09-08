
#ifndef MonitorSD_h
#define MonitorSD_h 1

#include "globals.hh"
#include "G4SystemOfUnits.hh"
#include "G4VSensitiveDetector.hh"
#include "MonitorHit.hh"


//#include <vector>

class G4Step;
class G4HCofThisEvent;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

/// Hodoscope paddle sensitive detector class

class MonitorSD : public G4VSensitiveDetector
{
  public:
    MonitorSD(G4String name, G4String hitsCollectionName);
    virtual ~MonitorSD();
  
    // methods from base class
    virtual void   Initialize(G4HCofThisEvent* hitCollection);
    virtual G4bool ProcessHits(G4Step* step, G4TouchableHistory* history);
    virtual void   EndOfEvent(G4HCofThisEvent* hitCollection);

  private:
	MonitorHitCollection* fHitsCollection;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
