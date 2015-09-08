#ifndef MonitorHit_h
#define MonitorHit_h 1

#include "globals.hh"
#include "G4SystemOfUnits.hh"
#include "G4VHit.hh"
#include "G4THitsCollection.hh"
#include "G4Allocator.hh"
#include "G4ThreeVector.hh"

/// Hodoscope paddle hit class
///
/// It defines data members to store the trackID, hodoscopeNb, paddleNb energy deposit,
/// light output, hit time.
/// - fTrackID, fhodoscopeNB, fpaddleNB fEdep, flight, ftime

class MonitorHit : public G4VHit
{
  public:
    MonitorHit();
    MonitorHit(const MonitorHit&);
    virtual ~MonitorHit();

    // operators
    const MonitorHit& operator=(const MonitorHit&);
    G4int operator==(const MonitorHit&) const;

    inline void* operator new(size_t);
    inline void  operator delete(void*);

    // methods from base class
    virtual void Draw();
    virtual void Print();

    // Set methods
    void SetTrackID(G4int track) { fTrackID = track; };
    void SetMomentumDirection(G4ThreeVector vec) { fMomentumDirection = vec; };
    void SetGlobalPrePosition(G4ThreeVector vec) { fGlobalPrePosition = vec; };
    void SetParticle(G4String part) { fParticle = part; };
    void SetKineticEnergy(G4double e) { fKineticEnergy = e; };
    void SetTime(G4double time) { fTime = time; };


    // Get methods
    G4int GetTrackID() const { return fTrackID; };
    G4ThreeVector GetMomentumDirection() const { return fMomentumDirection; };
    G4ThreeVector GetGlobalPrePosition() const { return fGlobalPrePosition; };
    G4String GetParticle() const { return fParticle; };
    G4double GetKineticEnergy() const { return fKineticEnergy; };
    G4double GetTime() const { return fTime; };

  private:

      G4int         fTrackID;
      G4String	    fParticle;
      G4ThreeVector fMomentumDirection;
      G4ThreeVector fGlobalPrePosition;
      G4double      fKineticEnergy;
      G4double      fTime;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

typedef G4THitsCollection<MonitorHit> MonitorHitCollection;

extern G4Allocator<MonitorHit> MonitorHitAllocator;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

inline void* MonitorHit::operator new(size_t)
{
  void *hit;
  hit = (void *) MonitorHitAllocator.MallocSingle();
  return hit;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

inline void MonitorHit::operator delete(void *hit)
{
  MonitorHitAllocator.FreeSingle((MonitorHit*) hit);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
