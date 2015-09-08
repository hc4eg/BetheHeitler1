
#ifndef PaddleSD_h
#define PaddleSD_h 1

#include "globals.hh"
#include "G4SystemOfUnits.hh"
#include "G4VSensitiveDetector.hh"
#include "PaddleHit.hh"


//#include <vector>

class G4Step;
class G4HCofThisEvent;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

/// Hodoscope paddle sensitive detector class

class PaddleSD : public G4VSensitiveDetector
{
  public:
    PaddleSD(G4String name, G4String hitsCollectionName);
    virtual ~PaddleSD();
  
    // methods from base class
    virtual void   Initialize(G4HCofThisEvent* hitCollection);
    virtual G4bool ProcessHits(G4Step* step, G4TouchableHistory* history);
    virtual void   EndOfEvent(G4HCofThisEvent* hitCollection);

  private:
	PaddleHitCollection* fHitsCollection;
	G4int fparentageHodoscope;
	G4int fparentagePaddle;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
