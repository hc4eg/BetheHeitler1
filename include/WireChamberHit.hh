/// \brief Definition of the WireChamberHit

#ifndef WireChamberHit_h
#define WireChamberHit_h 1

#include "globals.hh"
#include "G4SystemOfUnits.hh"
#include "G4VHit.hh"
#include "G4THitsCollection.hh"
#include "G4Allocator.hh"
#include "G4ThreeVector.hh"

/// Tracker hit class
///
/// It defines data members to store the trackID, chamberNb, energy deposit,
/// and position of charged particles in a selected volume:
/// - fTrackID, fChamberNB, fEdep, fPos

class WireChamberHit : public G4VHit
{
  public:
    WireChamberHit();
    WireChamberHit(const WireChamberHit&);
    virtual ~WireChamberHit();

    // operators
    const WireChamberHit& operator=(const WireChamberHit&);
    G4int operator==(const WireChamberHit&) const;

    inline void* operator new(size_t);
    inline void  operator delete(void*);

    // methods from base class
    virtual void Draw();
    virtual void Print();

    // Set methods
    void SetTrackID(G4int track)        { fTrackID = track; };
    void SetVDCnumber(G4int chamb)      { fVDCnumber = chamb; };
    void SetVDClayer(G4int layer)       { fVDClayer = layer; };
    void SetVDCpackage(G4int package)   { fVDCpackage = package; };
    void SetEdep(G4double de)           { fEdep = de; };
    void SetTime(G4double time)         { fTime = time; };
    void SetGlobalPrePosition(G4ThreeVector pos)  { fGlobalPrePosition = pos; };
    void SetLocalPrePosition(G4ThreeVector pos)   { fLocalPrePosition = pos; };
    void SetGlobalPostPosition(G4ThreeVector pos) { fGlobalPostPosition = pos; };
    void SetLocalPostPosition(G4ThreeVector pos)  { fLocalPostPosition = pos; };


    // Get methods
    G4int GetTrackID() const     { return fTrackID; };
    G4int GetVDCnumber() const   { return fVDCnumber; };
    G4int GetVDClayer() const    { return fVDClayer; };
    G4int GetVDCpackage() const  { return fVDCpackage; };
    G4double GetEdep() const     { return fEdep; };
    G4double GetTime() const     { return fTime; };
    G4ThreeVector GetGlobalPrePosition() const  { return fGlobalPrePosition; };
    G4ThreeVector GetLocalPrePosition() const   { return fLocalPrePosition; };
    G4ThreeVector GetGlobalPostPosition() const { return fGlobalPostPosition; };
    G4ThreeVector GetLocalPostPosition() const  { return fLocalPostPosition; };

  private:

      G4int         fTrackID;
      G4int         fVDCnumber;
      G4int         fVDClayer;
      G4int         fVDCpackage;
      G4double      fEdep;
      G4double      fTime;
      G4ThreeVector fGlobalPrePosition;
      G4ThreeVector fLocalPrePosition;
      G4ThreeVector fGlobalPostPosition;
      G4ThreeVector fLocalPostPosition;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

typedef G4THitsCollection<WireChamberHit> WireChamberHitsCollection;

extern G4Allocator<WireChamberHit> WireChamberHitAllocator;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

inline void* WireChamberHit::operator new(size_t)
{
  void *hit;
  hit = (void *) WireChamberHitAllocator.MallocSingle();
  return hit;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

inline void WireChamberHit::operator delete(void *hit)
{
  WireChamberHitAllocator.FreeSingle((WireChamberHit*) hit);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
