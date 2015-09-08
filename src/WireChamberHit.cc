#include "WireChamberHit.hh"
#include "G4UnitsTable.hh"
#include "G4VVisManager.hh"
#include "G4Circle.hh"
#include "G4Colour.hh"
#include "G4VisAttributes.hh"

#include <iomanip>

G4Allocator<WireChamberHit> WireChamberHitAllocator;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

WireChamberHit::WireChamberHit()
 : G4VHit(),
   fTrackID(-1),
   fVDCnumber(-1),
   fVDClayer(-1),
   fEdep(0.),
   fTime(0.),
   fGlobalPrePosition(G4ThreeVector()),
   fLocalPrePosition(G4ThreeVector()),
   fGlobalPostPosition(G4ThreeVector()),
   fLocalPostPosition(G4ThreeVector())
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

WireChamberHit::~WireChamberHit() {}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

WireChamberHit::WireChamberHit(const WireChamberHit& right)
  : G4VHit()
{
  fTrackID   = right.fTrackID;
  fVDCnumber = right.fVDCnumber;
  fVDClayer = right.fVDClayer;
  fEdep      = right.fEdep;
  fTime      = right.fTime;
  fGlobalPrePosition       = right.fGlobalPrePosition;
  fLocalPrePosition       = right.fLocalPrePosition;
  fGlobalPostPosition       = right.fGlobalPostPosition;
  fLocalPostPosition       = right.fLocalPostPosition;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

const WireChamberHit& WireChamberHit::operator=(const WireChamberHit& right)
{
  fTrackID   = right.fTrackID;
  fVDCnumber = right.fVDCnumber;
  fVDClayer = right.fVDClayer;
  fEdep      = right.fEdep;
  fTime      = right.fTime;
  fGlobalPrePosition       = right.fGlobalPrePosition;
  fLocalPrePosition       = right.fLocalPrePosition;
  fGlobalPostPosition       = right.fGlobalPostPosition;
  fLocalPostPosition       = right.fLocalPostPosition;

  return *this;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4int WireChamberHit::operator==(const WireChamberHit& right) const
{
  return ( this == &right ) ? 1 : 0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void WireChamberHit::Draw()
{
  G4VVisManager* pVVisManager = G4VVisManager::GetConcreteInstance();
  if(pVVisManager)
  {
    G4Circle circle(fGlobalPrePosition);
    circle.SetScreenSize(4.);
    circle.SetFillStyle(G4Circle::filled);
    G4Colour colour(1.,0.,0.);
    G4VisAttributes attribs(colour);
    circle.SetVisAttributes(attribs);
    pVVisManager->Draw(circle);
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void WireChamberHit::Print()
{
  G4cout
     << "  trackID: " << fTrackID << " VDCnumber: " << fVDCnumber
     << "  layer: " << fVDClayer << G4endl;
  G4cout.precision(4);
  G4cout  << "Edep: "
     << std::setw(7) << G4BestUnit(fEdep,"Energy") 
     << "   Time: "
     << std::setw(7) << G4BestUnit(fTime,"Time") 
	<< G4endl;
  G4cout << " Position Global Pre: "
     << std::setw(7) << G4BestUnit( fGlobalPrePosition,"Length") <<G4endl;
  G4cout << " Position Global Post: "
     << std::setw(7) << G4BestUnit( fGlobalPostPosition,"Length") <<G4endl;
  G4cout << " Position Local Pre: "
     << std::setw(7) << G4BestUnit( fLocalPrePosition,"Length") <<G4endl;
  G4cout << " Position Local Post: "
     << std::setw(7) << G4BestUnit( fLocalPostPosition,"Length") <<G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
