#include "PaddleHit.hh"
#include "G4UnitsTable.hh"
#include "G4VVisManager.hh"
#include "G4Circle.hh"
#include "G4Colour.hh"
#include "G4VisAttributes.hh"

#include <iomanip>

G4Allocator<PaddleHit> PaddleHitAllocator;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PaddleHit::PaddleHit()
 : G4VHit(),
   fTrackID(-1),
   fHodoscopeNB(-1),
   fPaddleNB(-1),
   fEdep(0.),
   fLight(0.),
   fTime(0.)
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PaddleHit::~PaddleHit() {}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PaddleHit::PaddleHit(const PaddleHit& right)
  : G4VHit()
{
  fTrackID   = right.fTrackID;
  fHodoscopeNB = right.fHodoscopeNB;
  fPaddleNB = right.fPaddleNB;
  fEdep      = right.fEdep;
  fLight      = right.fLight;
  fTime      = right.fTime;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

const PaddleHit& PaddleHit::operator=(const PaddleHit& right)
{
  fTrackID   = right.fTrackID;
  fHodoscopeNB = right.fHodoscopeNB;
  fPaddleNB = right.fPaddleNB;
  fEdep      = right.fEdep;
  fLight      = right.fLight;
  fTime      = right.fTime;
  return *this;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4int PaddleHit::operator==(const PaddleHit& right) const
{
  return ( this == &right ) ? 1 : 0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PaddleHit::Draw()
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PaddleHit::Print()
{
  G4cout
     << "  trackID: " << fTrackID << " Hodoscope: " << fHodoscopeNB
     << "  Paddle: " << fPaddleNB << G4endl;
  G4cout.precision(4);
  G4cout 
	<< " Edep: " << std::setw(7) << G4BestUnit(fEdep,"Energy") 
	<< " Light: " << std::setw(7) << G4BestUnit(fLight,"Energy") 
	<< " Time: " << std::setw(7) << G4BestUnit(fTime,"Time") 
	<< G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
