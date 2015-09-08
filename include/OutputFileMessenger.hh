#ifndef OutputFileMessenger_hh
#define OutputFileMessenger_hh 1

#include "OutputFile.hh"
#include "G4UImessenger.hh"
#include "globals.hh"
#include "G4SystemOfUnits.hh"
#include "G4UIdirectory.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithoutParameter.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWithADouble.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "G4UIcmdWithABool.hh"

class OutputFileMessenger: public G4UImessenger
{
   public:
	OutputFileMessenger( OutputFile * out);
	~OutputFileMessenger();

   public:
	void SetNewValue(G4UIcommand * command, G4String newValues);
	//G4String GetCurrentValue(G4UIcommand * command);
   private:
	OutputFile * fout;
	G4UIdirectory * Dir1;

	G4UIcmdWithoutParameter * PrintCmd;
	G4UIcmdWithAString * Set_prefix_cmd;
	G4UIcmdWithAnInteger * Set_header_frequency_cmd;
	G4UIcmdWithABool * Set_unique_cmd;
	G4UIcmdWithABool * Set_use_monitor_cmd;
};
#endif
