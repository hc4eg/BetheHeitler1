//

#include "PaddleSD.hh"
#include "G4HCofThisEvent.hh"
#include "G4Step.hh"
#include "G4ThreeVector.hh"
#include "G4TwoVector.hh"
#include "G4SDManager.hh"
#include "G4ios.hh"
#include "G4VSolid.hh"
#include "G4DisplacedSolid.hh"
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PaddleSD::PaddleSD(G4String name, G4String hitsCollectionName) 
 : G4VSensitiveDetector(name)
{
  collectionName.insert(hitsCollectionName);
  verboseLevel = 0;
  fparentageHodoscope = 2;
  fparentagePaddle = 1;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PaddleSD::~PaddleSD() 
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PaddleSD::Initialize(G4HCofThisEvent* hce)
{
  // Create hits collection

  fHitsCollection 
    = new PaddleHitCollection(SensitiveDetectorName, collectionName[0]); 

  // Add this collection in hce

  G4int hcID = G4SDManager::GetSDMpointer()->GetCollectionID(fHitsCollection);
  hce->AddHitsCollection( hcID, fHitsCollection ); 

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4bool PaddleSD::ProcessHits(G4Step* aStep, G4TouchableHistory*)
{  
  // energy deposit
  G4double edep = aStep->GetTotalEnergyDeposit();
  if(edep == 0.) return false;
  G4StepPoint* preStepPoint = aStep->GetPreStepPoint();
  G4double time = preStepPoint->GetGlobalTime();
  G4double dtime = aStep->GetDeltaTime();
  time += dtime/2.;
  //G4double preEng = preStepPoint->GetKineticEnergy();

  //G4TouchableHistory* theTouchable
  // = (G4TouchableHistory*)(preStepPoint->GetTouchable());
  //G4VPhysicalVolume* thePhysical = theTouchable->GetVolume();
  G4TouchableHandle theTouchHandle = preStepPoint->GetTouchableHandle();

  G4int Paddle = theTouchHandle->GetReplicaNumber(fparentagePaddle);
  G4int Hodoscope = theTouchHandle->GetReplicaNumber(fparentageHodoscope);
  
//  G4cout << "======= Begin Process Hits =========" << G4endl;
#if 0
//    This is just for verifying what parentage is appropriate 
  G4cout << "======= Begin Process Hits =========" << G4endl;
for(G4int fparentage = 0; fparentage < 3; fparentage++)
  {
  G4int number = aStep->GetPreStepPoint()->GetTouchable()->GetReplicaNumber(fparentage);
  G4String name = aStep->GetPreStepPoint()->GetTouchable()->GetVolume(fparentage)->GetName();
  G4cout << "Parentage = "<<fparentage<<" Volume name = " << name << "  copy number = " << number << G4endl;
  }
#endif
// -----------------------------------
  PaddleHit* newHit = new PaddleHit();
  fHitsCollection->insert( newHit );

  // get the other stuff we want to record in the hit

   newHit->SetTrackID  (aStep->GetTrack()->GetTrackID());
   newHit->SetPaddleNB(Paddle);
   newHit->SetHodoscopeNB(Hodoscope);
   newHit->SetEdep(edep);
	// for now not implementing light output response
   newHit->SetLight(edep);
   newHit->SetTime(time);

   //newHit->Print();

  return true;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PaddleSD::EndOfEvent(G4HCofThisEvent*)
{
  if ( verboseLevel>1 ) { 
     G4int nofHits = fHitsCollection->entries();
     G4cout << "\n-------->Hits Collection: in this event they are " << nofHits 
            << " hits in the VDC chambers: " << G4endl;
     for ( G4int i=0; i<nofHits; i++ ) (*fHitsCollection)[i]->Print();
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
