#include "MonitorHit.hh"
#include "G4UnitsTable.hh"
#include "G4VVisManager.hh"
#include "G4Circle.hh"
#include "G4Colour.hh"
#include "G4VisAttributes.hh"

#include <iomanip>

G4Allocator<MonitorHit> MonitorHitAllocator;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

MonitorHit::MonitorHit()
 : G4VHit(),
   fTrackID(-1),
   fParticle(""),
   fMomentumDirection(G4ThreeVector(1.,0.,0.)),
   fGlobalPrePosition(G4ThreeVector(0.,0.,0.)),
   fKineticEnergy(0.),
   fTime(0.)
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

MonitorHit::~MonitorHit() {}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

MonitorHit::MonitorHit(const MonitorHit& right)
  : G4VHit()
{
  fTrackID   = right.fTrackID;
  fMomentumDirection = right.fMomentumDirection;
  fGlobalPrePosition = right.fGlobalPrePosition;
  fKineticEnergy = right.fKineticEnergy;
  fParticle      = right.fParticle;
  fTime      = right.fTime;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

const MonitorHit& MonitorHit::operator=(const MonitorHit& right)
{
  fTrackID   = right.fTrackID;
  fMomentumDirection = right.fMomentumDirection;
  fGlobalPrePosition = right.fGlobalPrePosition;
  fKineticEnergy = right.fKineticEnergy;
  fParticle      = right.fParticle;
  fTime      = right.fTime;
  return *this;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4int MonitorHit::operator==(const MonitorHit& right) const
{
  return ( this == &right ) ? 1 : 0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void MonitorHit::Draw()
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void MonitorHit::Print()
{
  G4cout
     << " Monitor: TrackID: " << fTrackID 
     << " Particle = " << fParticle << G4endl;
  G4cout.precision(4);
  G4cout 
	<< " Kinetic Energy: " << std::setw(7) << G4BestUnit(fKineticEnergy,"Energy") 
	<< " Time: " << std::setw(7) << G4BestUnit(fTime,"Time") 
	<< G4endl;
  G4cout << " Position Global Pre: "
     << std::setw(7) << G4BestUnit( fGlobalPrePosition,"Length") <<G4endl;
  G4cout << " Momentum Direction: "
     << std::setw(7) << fMomentumDirection <<G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
