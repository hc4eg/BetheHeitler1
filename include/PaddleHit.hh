#ifndef PaddleHit_h
#define PaddleHit_h 1

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

class PaddleHit : public G4VHit
{
  public:
    PaddleHit();
    PaddleHit(const PaddleHit&);
    virtual ~PaddleHit();

    // operators
    const PaddleHit& operator=(const PaddleHit&);
    G4int operator==(const PaddleHit&) const;

    inline void* operator new(size_t);
    inline void  operator delete(void*);

    // methods from base class
    virtual void Draw();
    virtual void Print();

    // Set methods
    void SetTrackID(G4int track)        { fTrackID = track; };
    void SetHodoscopeNB(G4int num)      { fHodoscopeNB = num; };
    void SetPaddleNB(G4int num)      { fPaddleNB = num; };
    void SetEdep(G4double de)           { fEdep = de; };
    void SetLight(G4double de)           { fLight = de; };
    void SetTime(G4double time)         { fTime = time; };


    // Get methods
    G4int GetTrackID() const     { return fTrackID; };
    G4int GetHodoscopeNB() const   { return fHodoscopeNB; };
    G4int GetPaddleNB() const   { return fPaddleNB; };
    G4double GetEdep() const     { return fEdep; };
    G4double GetLight() const     { return fLight; };
    G4double GetTime() const     { return fTime; };

  private:

      G4int         fTrackID;
      G4int         fHodoscopeNB;
      G4int         fPaddleNB;
      G4double      fEdep;
      G4double      fLight;
      G4double      fTime;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

typedef G4THitsCollection<PaddleHit> PaddleHitCollection;

extern G4Allocator<PaddleHit> PaddleHitAllocator;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

inline void* PaddleHit::operator new(size_t)
{
  void *hit;
  hit = (void *) PaddleHitAllocator.MallocSingle();
  return hit;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

inline void PaddleHit::operator delete(void *hit)
{
  PaddleHitAllocator.FreeSingle((PaddleHit*) hit);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
