
#ifndef WireChamberSD_h
#define WireChamberSD_h 1

#include "globals.hh"
#include "G4SystemOfUnits.hh"
#include "G4VSensitiveDetector.hh"
#include "WireChamberHit.hh"


//#include <vector>

class G4Step;
class G4HCofThisEvent;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

/// WireChamber sensitive detector class

class WireChamberSD : public G4VSensitiveDetector
{
  public:
    WireChamberSD(G4String name, G4String hitsCollectionName);
    virtual ~WireChamberSD();
  
    // methods from base class
    virtual void   Initialize(G4HCofThisEvent* hitCollection);
    virtual G4bool ProcessHits(G4Step* step, G4TouchableHistory* history);
    virtual void   EndOfEvent(G4HCofThisEvent* hitCollection);

  private:
	WireChamberHitsCollection* fHitsCollection;
	G4int fparentageLayer;
	G4int fparentageVDC;
	G4int fparentagePackage;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
