#ifndef OutputFile_h
#define OutputFile_h 1

#include "globals.hh"
#include "G4SystemOfUnits.hh"
class OutputFileMessenger;

#include <stdio.h>

class OutputFile {
	public:
	static OutputFile *  GetOutputFilePointer(void);
	~OutputFile();
	G4bool OpenFile(G4int runno);
	void CloseFile(void);
	void WriteEvent(void);
	void WriteComment(G4String line);
	void PrintParameters();

	inline void Set_event_number(G4int val) {fevent_number = val;};
	inline void Set_energy_i(G4double val) {fenergy_i = val;};
	inline void Set_delta_i(G4double val) {fdelta_i = val;};
	inline void Set_x_i(G4double val) {fx_i = val;};
	inline void Set_y_i(G4double val) {fy_i = val;};
	inline void Set_theta_i(G4double val) {ftheta_i = val;};
	inline void Set_phi_i(G4double val) {fphi_i = val;};

	inline void Set_detector_package(G4int i, G4bool val) {fdetector_package[i] = val;};
	inline void Set_x_f(G4int i, G4double val) {fx_f[i] = val;};
	inline void Set_y_f(G4int i, G4double val) {fy_f[i] = val;};
	inline void Set_theta_f(G4int i, G4double val) {ftheta_f[i] = val;};
	inline void Set_phi_f(G4int i, G4double val) {fphi_f[i] = val;};

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

	// data to be written per event
	G4int fevent_number;
	// input data
	G4double fx_i, fy_i, ftheta_i, fphi_i;
	G4double fenergy_i, fdelta_i;
	// which detector package got a hit
	G4bool fdetector_package[2];
	G4double fx_f[2], fy_f[2], ftheta_f[2], fphi_f[2];
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
	FILE * fd;
	G4bool fFile_open;
	G4String fFilename_prefix;
	char filename[100];
	G4bool fuse_unique_filename;
	
};	

#endif
