//

#include "MonitorSD.hh"
#include "G4HCofThisEvent.hh"
#include "G4Step.hh"
#include "G4ThreeVector.hh"
#include "G4TwoVector.hh"
#include "G4SDManager.hh"
#include "G4ios.hh"
#include "G4VSolid.hh"
#include "G4DisplacedSolid.hh"
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

MonitorSD::MonitorSD(G4String name, G4String hitsCollectionName) 
 : G4VSensitiveDetector(name)
{
  collectionName.insert(hitsCollectionName);
  verboseLevel = 0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

MonitorSD::~MonitorSD() 
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void MonitorSD::Initialize(G4HCofThisEvent* hce)
{
  // Create hits collection

  fHitsCollection 
    = new MonitorHitCollection(SensitiveDetectorName, collectionName[0]); 

  // Add this collection in hce

  G4int hcID = G4SDManager::GetSDMpointer()->GetCollectionID(fHitsCollection);
  hce->AddHitsCollection( hcID, fHitsCollection ); 

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4bool MonitorSD::ProcessHits(G4Step* aStep, G4TouchableHistory*)
{  
  // energy deposit
  G4StepPoint* preStepPoint = aStep->GetPreStepPoint();
  //G4double time = preStepPoint->GetGlobalTime();
  //G4double preEng = preStepPoint->GetKineticEnergy();

  //G4TouchableHistory* theTouchable
  // = (G4TouchableHistory*)(preStepPoint->GetTouchable());
  //G4VPhysicalVolume* thePhysical = theTouchable->GetVolume();
  //G4TouchableHandle theTouchHandle = preStepPoint->GetTouchableHandle();

  //G4int Monitor = theTouchHandle->GetReplicaNumber(fparentageMonitor);
  //G4int Hodoscope = theTouchHandle->GetReplicaNumber(fparentageHodoscope);
  
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
  MonitorHit* newHit = new MonitorHit();
  fHitsCollection->insert( newHit );

  // get the other stuff we want to record in the hit

   newHit->SetTrackID  (aStep->GetTrack()->GetTrackID());
   newHit->SetTime(preStepPoint->GetGlobalTime());
   newHit->SetGlobalPrePosition(preStepPoint->GetPosition());
   newHit->SetMomentumDirection(preStepPoint->GetMomentumDirection());
   newHit->SetKineticEnergy(preStepPoint->GetKineticEnergy());
   //particle's name
   G4Track * aTrack = aStep->GetTrack();
   G4ParticleDefinition * aParticle = aTrack->GetDefinition();
   G4String part = aParticle->GetParticleName();
   newHit->SetParticle(part);
	// for now not implementing light output response

   //newHit->Print();

  return true;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void MonitorSD::EndOfEvent(G4HCofThisEvent*)
{
  if ( verboseLevel>1 ) { 
     G4int nofHits = fHitsCollection->entries();
     G4cout << "\n--->Hits Collection: in this event they are " << nofHits 
            << " hits in the Monitor: " << G4endl;
     for ( G4int i=0; i<nofHits; i++ ) (*fHitsCollection)[i]->Print();
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
