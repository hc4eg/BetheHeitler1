#ifndef EventAction_h
#define EventAction_h 1

#include "G4UserEventAction.hh"
#include "globals.hh"
#include "G4SystemOfUnits.hh"
#include <fstream>
#include <iostream>
using namespace std;

#include "G4SDManager.hh"
#include "WireChamberHit.hh"
#include "PaddleHit.hh"
#include "MonitorHit.hh"
#include "G4TwoVector.hh"
#include "OutputFile.hh"
#include "DetectorConstruction.hh"

#define NUMPADDLES 29

// New class representing hit data of a single wire hit collection step used in EventAction.cc
class SingleWireHit{
public:
  G4double Get_KE() {return KE;}
  G4double Get_ToF() {return ToF;}
  G4double Get_Charge() {return Charge;}
  G4int    Get_Particle() {return Particle;}
  void Clear() { KE = 0.; ToF = 0.; Charge = 0.;}
  void Set_KE(G4double val) { KE = val;}
  void Set_ToF(G4double val) { ToF = val;}
  void Set_Charge(G4double val) { Charge = val;}
  void Set_Particle(G4int val) { Particle = val;}
private:
  G4double KE;
  G4double ToF;
  //PDGEncoding to specify particle type
  //See pdg.lbl.gov/2007/reviews/montecarlorpp.pdf
  //For reference: e+ -11, e- +11, gamma +22.
  G4int Particle;
  G4double Charge;
};

// New class representing hit data of a single wire used in EventAction.cc
class WireHit{
public:
  //G4double* Get_Edep() {return &Edep;}
  //G4ThreeVector* Get_Position() {return &Position;}
  
  void Clear() { 
    WireNum = -1; // The lowest indexed wire WireNum = 0;
    Edep = 0.;
    Position = G4ThreeVector(0.,0.,0.); 
    //if(SWH.size()>0) SWH.clear();
  }
  
  void Set_WireNum(G4int val) { WireNum = val;}
  void Set_Edep(G4double val) { Edep = val;}
  void Set_Position(G4ThreeVector val) { 
    Position.setX(val.getX()); 
    Position.setY(val.getY());
    Position.setZ(val.getZ());}

  //Temporary as public data for easy access.
  vector<SingleWireHit> SWH;
  G4int WireNum;
  G4double Edep;
  G4ThreeVector Position;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class EventAction : public G4UserEventAction
{
public:
  EventAction(DetectorConstruction*);
  virtual ~EventAction();

  void  BeginOfEventAction(const G4Event*);
  void    EndOfEventAction(const G4Event*);
    
private:
   G4ThreeVector fR_i;
                    
   G4int     printModulo;
                             
   //EventActionMessenger*  eventMessenger;
   OutputFile * pOutputFile;
   G4double fVDCSpacing;
   DetectorConstruction* pDetector;
   G4double fVDCSigma; // std. dev. of position resolution
   G4bool fRequireHodoscopeHit;
   G4double fPaddleThreshold;
  // New member, assigning threshold of energy deposite of VDC wires.
   G4double fWireThreshold;
   G4double fsigmaLightStatistical;
   G4double fsigmaLightNoise;

};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif

    
