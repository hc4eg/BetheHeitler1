//

#include "WireChamberSD.hh"
#include "G4HCofThisEvent.hh"
#include "G4Step.hh"
#include "G4ThreeVector.hh"
#include "G4TwoVector.hh"
#include "G4SDManager.hh"
#include "G4ios.hh"
#include "G4VSolid.hh"
#include "G4DisplacedSolid.hh"
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

WireChamberSD::WireChamberSD(G4String name, G4String hitsCollectionName) 
 : G4VSensitiveDetector(name)
{
  collectionName.insert(hitsCollectionName);
  verboseLevel = 0;
  fparentageLayer = 0;
  fparentageVDC = 1;
  fparentagePackage = 2;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

WireChamberSD::~WireChamberSD() 
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void WireChamberSD::Initialize(G4HCofThisEvent* hce)
{
  // Create hits collection

  fHitsCollection 
    = new WireChamberHitsCollection(SensitiveDetectorName, collectionName[0]); 

  // Add this collection in hce

  G4int hcID = G4SDManager::GetSDMpointer()->GetCollectionID(fHitsCollection);
  hce->AddHitsCollection( hcID, fHitsCollection ); 

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4bool WireChamberSD::ProcessHits(G4Step* aStep, 
                                     G4TouchableHistory*)
{  
  // energy deposit
  G4double edep = aStep->GetTotalEnergyDeposit();
  if(edep == 0.) return false;
  G4StepPoint* preStepPoint = aStep->GetPreStepPoint();

  //G4TouchableHistory* theTouchable
  // = (G4TouchableHistory*)(preStepPoint->GetTouchable());
  //G4VPhysicalVolume* thePhysical = theTouchable->GetVolume();
  G4TouchableHandle theTouchHandle = preStepPoint->GetTouchableHandle();

  G4int VDClayer = theTouchHandle->GetReplicaNumber(fparentageLayer);
  G4int VDCnumber = theTouchHandle->GetReplicaNumber(fparentageVDC);
  G4int VDCpackage = theTouchHandle->GetReplicaNumber(fparentagePackage);
  
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
  WireChamberHit* newHit = new WireChamberHit();
  fHitsCollection->insert( newHit );

  // get the other stuff we want to record in the hit
  // what is important to us is the position of the track
  // since that is what we will get from the the VDC
  G4double time = preStepPoint->GetGlobalTime();
  const G4AffineTransform transformation = preStepPoint->GetTouchable()->GetHistory()->GetTopTransform();
   G4ThreeVector globalPrePosition = preStepPoint->GetPosition();
   G4ThreeVector localPrePosition = transformation.TransformPoint(globalPrePosition);
   G4ThreeVector globalPostPosition = aStep->GetPostStepPoint()->GetPosition();
   G4ThreeVector localPostPosition = transformation.TransformPoint(globalPostPosition);

   newHit->SetTrackID  (aStep->GetTrack()->GetTrackID());
   newHit->SetVDCnumber(VDCnumber);
   newHit->SetVDClayer(VDClayer);
   newHit->SetVDCpackage(VDCpackage);
   newHit->SetEdep(edep);
   newHit->SetTime(time);
   newHit->SetGlobalPrePosition(globalPrePosition);
   newHit->SetLocalPrePosition(localPrePosition);
   newHit->SetGlobalPostPosition(globalPostPosition);
   newHit->SetLocalPostPosition(localPostPosition);

 // newHit->Print();

  return true;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void WireChamberSD::EndOfEvent(G4HCofThisEvent*)
{
  if ( verboseLevel>1 ) { 
     G4int nofHits = fHitsCollection->entries();
     G4cout << "\n-------->Hits Collection: in this event they are " << nofHits 
            << " hits in the VDC chambers: " << G4endl;
     for ( G4int i=0; i<nofHits; i++ ) (*fHitsCollection)[i]->Print();
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
