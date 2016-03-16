#ifndef OutputFile_h
#define OutputFile_h 1

#include "globals.hh"
#include "G4SystemOfUnits.hh"
class OutputFileMessenger;

#include <stdio.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include "TFile.h"
#include "TTree.h"
#include "TObject.h"

using namespace std;

class Monitor{
public:
    Float_t Charge;
    Float_t Energy;
    Float_t X;
    Float_t Y;
    Float_t Theta;
    Float_t Phi;
    
  inline void SetMonitor(Monitor* M){
    Charge = M->Charge; Energy = M->Energy; X = M->X;
    Y = M->Y; Theta = M->Theta; Phi = M->Phi;};
  };

struct Input{
  Float_t Energy;
  Float_t Delta;
  Float_t X;
  Float_t Y;
  Float_t Theta;
  Float_t Phi;

  inline void SetInput(Input* I){
    Energy = I->Energy; Delta = I->Delta; X = I->X;
    Y = I->Y; Theta = I->Theta; Phi = I->Phi;};
};

struct Wire{
  Int_t   WireNum;
  Float_t X;
  Float_t Y;
  Float_t Edep;
  Float_t KE;
  Float_t ToF;
  Float_t Charge;
  Int_t   Particle;
};

struct Paddle{
  Int_t PNum;
  Float_t Edep;
  Float_t Light;
  Float_t Time;
  //Float_t KE;
};

struct Detector{
  //vector<Wire> V;
  // Now assign 4 wireplanes on each side detector.
  vector<Wire> V0U, V0V, V1U, V1V;
  vector<Paddle> P;
};

class BH_Event : public TObject{
public:
  Int_t ENum;
  Monitor M0;
  Monitor M1;
  Input I0;
  Input I1;
  Detector D0;
  Detector D1;
  vector<Double_t> HKE0;
  vector<Double_t> HKE1;

  private:
    ClassDef(BH_Event,1)
};

class OutputWire{
public:
  inline void Clear_Wire () {
    WireNum = 0;
    X = 0.;
    Y = 0.;
    Edep = 0.;
    KE = 0.;
    ToF = 0.;
    Charge = 0.;
    Particle = 0;
  };
  inline void Set_WireNum_f (G4int val) { WireNum = val; };
  inline void Set_X_f (G4double val) { X = val; };
  inline void Set_Y_f (G4double val) { Y = val; };
  inline void Set_Edep_f (G4double val) { Edep = val; };
  inline void Set_KE_f (G4double val) { KE = val; };
  inline void Set_ToF_f (G4double val) { ToF = val; };
  inline void Set_Charge_f (G4double val) { Charge = val; };
  inline void Set_Particle_f (G4int val) { Particle = val; };
  G4int Get_WireNum_f() { return WireNum; };
  G4double Get_X_f() { return X;};
  G4double Get_Y_f() { return Y;};
  G4double Get_Edep_f() { return Edep;};
  G4double Get_KE_f() { return KE;};
  G4double Get_ToF_f() { return ToF;};
  G4double Get_Charge_f() { return Charge;};
  G4int Get_Particle_f() { return Particle;};

private:
  G4int WireNum;
  G4double Edep;
  G4double X,Y;
  G4double KE, ToF, Charge;
  G4int Particle;
};


class OutputFile {
	public:
	static OutputFile *  GetOutputFilePointer(void);
	~OutputFile();
	G4bool OpenFile(G4int runno);
	void CloseFile(void);
	void WriteEvent(void);
	void WriteComment(G4String line);
	void PrintParameters();


        inline void Set_out_pair_mode(G4bool val){out_pair_mode = val;};

	inline void Set_event_number(G4int val) {fevent_number = val;};
  //Input
        inline void Set_energy_i(G4int n,G4double val) {fenergy_i[n] = val;};
        inline void Set_delta_i(G4int n, G4double val) {fdelta_i[n] = val;};
        inline void Set_x_i(G4int n,G4double val) {fx_i[n] = val;};
        inline void Set_y_i(G4int n,G4double val) {fy_i[n] = val;};
        inline void Set_theta_i(G4int n,G4double val) {ftheta_i[n] = val;};
        inline void Set_phi_i(G4int n,G4double val) {fphi_i[n] = val;};
  //vdc
	inline void Set_detector_package(G4int i, G4bool val) {fdetector_package[i] = val;};
	/*	inline void Set_x_f(G4int i, G4double val) {fx_f[i] = val;};
	inline void Set_y_f(G4int i, G4double val) {fy_f[i] = val;};
  //inline void Set_theta_f(G4int i, G4double val) {ftheta_f[i] = val;};
  //inline void Set_phi_f(G4int i, G4double val) {fphi_f[i] = val;};

        inline void Set_KE_f(G4int i, G4double val) {fKE_f[i] = val;};
  inline void Set_ToF_f(G4int i, G4double val) { fToF_f[i] = val;};
  inline void Set_Charge_f(G4int i, G4double val) {fCharge_f[i] = val;};
  */

  //Now introducing wires for each wireplane
  // Before storing new data in fVDC_f, Clear_Wire is needed
  /*
  inline void Clear_Wire () {
    fWireHit_f.WireNum = 0;
    fWireHit_f.X = 0.;
    fWireHit_f.Y = 0.;
    fWireHit_f.Edep = 0.;
    fWireHit_f.KE = 0.;
    fWireHit_f.ToF = 0.;
    fWireHit_f.Charge = 0.;
    //fWireHit_f.Particle = 0;
  };
  inline void Set_WireNum_f (G4int WireNum) { fWireHit_f.WireNum = WireNum; };
  inline void Set_X_f (G4double X) { fWireHit_f.X = X; };
  inline void Set_Y_f (G4double Y) { fWireHit_f.Y = Y; };
  inline void Set_Edep_f (G4double Edep) { fWireHit_f.Edep = Edep; };
  inline void Set_KE_f (G4double KE) { fWireHit_f.KE = KE; };
  inline void Set_ToF_f (G4double ToF) { fWireHit_f.ToF = ToF; };
  inline void Set_Charge_f (G4double Charge) { fWireHit_f.Charge = Charge; };
  //inline void Set_Particle_f (G4int Particle) { fWireHit_f.Particle = Particle; };
  */

  // Use fVDC_f objects to store temporary hit data extract from EventAction.cc on each of 4 wireplanes in a detector package
  // Clear fVDC_f objects
  inline void Clear_VDC_f(){
    	  for(G4int i = 0; i < 2; i++)
	    for(G4int j = 0; j < 2; j++)
	      for(G4int k = 0; k < 2; k++)
		fVDC_f[i][j][k].clear();}
  
  // Set a wire hit in fVDC_f, then move it to root file.
  inline void Set_VDC_f (G4int i, G4int j, G4int k, OutputWire Wire) {fVDC_f[i][j][k].push_back(Wire);}

  // Use Get_Wire_f to get Wire object from fVDC_f, used it to store data in root file 
  inline OutputWire Get_Wire_f (G4int i, G4int j, G4int k, G4int l) { return fVDC_f[i][j][k].at(l); };
  // Use this to get number of iteration when put fVDC_f[] data into root file
  inline G4int Get_NumHit_f (G4int i, G4int j, G4int k) { return fVDC_f[i][j][k].size(); };


  //Paddle
	inline void Set_pad_hit(G4int hod, G4int pad, G4bool val) { fPad_hit[hod][pad] = val;};
	inline void Set_hod_hit(G4int hod, G4bool val) { fHod_hit[hod] = val;};
	inline void Set_pad_energy(G4int hod, G4int pad, G4double val) { fPadEnergy[hod][pad] = val;};
	inline void Set_pad_light(G4int hod, G4int pad, G4double val) { fPadLight[hod][pad] = val;};
	inline void Set_pad_time(G4int hod, G4int pad, G4double val) { fPadTime[hod][pad] = val;};
        inline void Set_pad_KE(G4int hod, G4double val){ fPadKE[hod] = val;};

  //Monitor
	inline void SetUseMonitor(G4bool val) { fUseMonitor = val;};
	inline G4bool GetUseMonitor() { return fUseMonitor;};
	inline void SetMonitorKineticEnergy(G4int i, G4double val) { fMonitorKineticEnergy[i] = val;};
	inline void SetMonitorX(G4int i, G4double val) { fMonitorX[i] = val;};
	inline void SetMonitorY(G4int i, G4double val) { fMonitorY[i] = val;};
	inline void SetMonitorTheta(G4int i, G4double val) { fMonitorTheta[i] = val;};
	inline void SetMonitorPhi(G4int i, G4double val) { fMonitorPhi[i] = val;};
	inline void SetMonitorHit(G4int i, G4bool val) { fMonitorHit[i] = val;};

	inline void Set_filename_prefix(G4String val){fFilename_prefix = val;};
	inline G4String Get_filename_prefix(){return fFilename_prefix;};
	inline void Set_use_unique_filename(G4bool val) {fuse_unique_filename = val;};
	inline void Set_output_lines(G4int val) { foutput_lines = val;};
        inline void ClearKE (){fPadKE[0] = 0.; fPadKE[1] = 0.;};

	private:
	// private constructor
	OutputFile();

	// pointer OutputFile
	static OutputFile * pOutputFilePointer;

	OutputFileMessenger * fmessenger;

	G4int frunno;
	G4int foutput_lines;
	G4int flines_for_heading;
        G4bool out_pair_mode;
        //Root File Ouput
        TFile* bh_tree_file;
        TTree* tree;
        BH_Event* event;
        TBranch* branch;
	// data to be written per event
	G4int fevent_number;
	// input data
	G4double fx_i[2], fy_i[2], ftheta_i[2], fphi_i[2];
	G4double fenergy_i[2], fdelta_i[2];
	

  // which detector package got a hit
	G4bool fdetector_package[2];
	/*       G4int  fWire_f;
        G4double fx_f[2], fy_f[2];
        G4double fedep_f[2][2][2];
        G4double fKE_f[2], fToF_f[2], fCharge_f[2];
  */
  // Newly defined data structure including wirehit
  vector<OutputWire> fVDC_f[2][2][2];
  //OutputWire fWireHit_f;


  // which hodoscope and paddles got hit
	G4bool fHod_hit[2];
	G4bool fPad_hit[2][29];
	G4double fPadEnergy[2][29];
	G4double fPadLight[2][29];
	G4double fPadTime[2][29];
        G4double fPadKE[2];

	// Monitor output stuff
	G4bool fUseMonitor;
	G4bool fMonitorHit[2];
	G4double fMonitorKineticEnergy[2];
	G4double fMonitorX[2];
	G4double fMonitorY[2];
	G4double fMonitorTheta[2];
	G4double fMonitorPhi[2];


	// file handling stuff
	//FILE * fd;
	G4bool fFile_open;
	G4String fFilename_prefix;
	char filename[100];
	G4bool fuse_unique_filename;	
};	

#endif
