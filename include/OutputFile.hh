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

struct VDC{
    Float_t X;
    Float_t Y;
    Float_t Theta;
    Float_t Phi;
    Float_t E0u;
    Float_t E0v;
    Float_t E1u;
    Float_t E1v;
};

struct Paddle{
  Int_t PNum;
  Float_t Edep;
  Float_t Light;
  Float_t Time;
};

struct Detector{
  vector<VDC> V;
  vector<Paddle> P;
};

class BH_Event : public TObject{
//class BH_Event{
public:
  Int_t ENum;
  Monitor M1;
  Monitor M0;
  Input I0;
  Input I1;
  Detector D0;
  Detector D1;

  private:
    ClassDef(BH_Event,1)
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

        inline void Set_energy_i(G4int n,G4double val) {fenergy_i[n] = val;};
        inline void Set_delta_i(G4int n, G4double val) {fdelta_i[n] = val;};
        inline void Set_x_i(G4int n,G4double val) {fx_i[n] = val;};
        inline void Set_y_i(G4int n,G4double val) {fy_i[n] = val;};
        inline void Set_theta_i(G4int n,G4double val) {ftheta_i[n] = val;};
        inline void Set_phi_i(G4int n,G4double val) {fphi_i[n] = val;};

	inline void Set_detector_package(G4int i, G4bool val) {fdetector_package[i] = val;};
	inline void Set_x_f(G4int i, G4double val) {fx_f[i] = val;};
	inline void Set_y_f(G4int i, G4double val) {fy_f[i] = val;};
	inline void Set_theta_f(G4int i, G4double val) {ftheta_f[i] = val;};
	inline void Set_phi_f(G4int i, G4double val) {fphi_f[i] = val;};
        //Add edep for vdc
        inline void Set_edep_f(G4int i, G4int j, G4int k, G4double val){fedep_f[i][j][k] = val;};

	inline void Set_pad_hit(G4int hod, G4int pad, G4bool val) { fPad_hit[hod][pad] = val;};
	inline void Set_hod_hit(G4int hod, G4bool val) { fHod_hit[hod] = val;};
	inline void Set_pad_energy(G4int hod, G4int pad, G4double val) { fPadEnergy[hod][pad] = val;};
	inline void Set_pad_light(G4int hod, G4int pad, G4double val) { fPadLight[hod][pad] = val;};
	inline void Set_pad_time(G4int hod, G4int pad, G4double val) { fPadTime[hod][pad] = val;};

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
	G4double fx_f[2], fy_f[2], ftheta_f[2], fphi_f[2];

        G4double fedep_f[2][2][2];
	// which hodoscope and paddles got hit
	G4bool fHod_hit[2];
	G4bool fPad_hit[2][29];
	G4double fPadEnergy[2][29];
	G4double fPadLight[2][29];
	G4double fPadTime[2][29];
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
