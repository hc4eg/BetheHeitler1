#ifndef GeneratorMessenger_hh
#define GeneratorMessenger_hh 1

#include "PrimaryGeneratorAction.hh"
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
class PrimaryGeneratorAction;

class GeneratorMessenger: public G4UImessenger
{
   public:
	GeneratorMessenger( PrimaryGeneratorAction * gen);
	~GeneratorMessenger();

   public:
	void SetNewValue(G4UIcommand * command, G4String newValues);
	//G4String GetCurrentValue(G4UIcommand * command);
   private:
	PrimaryGeneratorAction * generator;
	G4UIdirectory * Dir1;

	G4UIcmdWithoutParameter * PrintCmd;
	G4UIcmdWithADoubleAndUnit * Set_central_energy_cmd;
	G4UIcmdWithADouble * Set_delta_min_cmd;
	G4UIcmdWithADouble * Set_delta_max_cmd;
	G4UIcmdWithADoubleAndUnit * Set_radius_max_cmd;
	G4UIcmdWithADoubleAndUnit * Set_x_min_cmd;
	G4UIcmdWithADoubleAndUnit * Set_x_max_cmd;
	G4UIcmdWithADoubleAndUnit * Set_y_min_cmd;
	G4UIcmdWithADoubleAndUnit * Set_y_max_cmd;
	G4UIcmdWithADoubleAndUnit * Set_theta_min_cmd;
	G4UIcmdWithADoubleAndUnit * Set_theta_max_cmd;
	G4UIcmdWithADoubleAndUnit * Set_phi_min_cmd;
	G4UIcmdWithADoubleAndUnit * Set_phi_max_cmd;
	G4UIcmdWithADoubleAndUnit * Set_angle_max_cmd;
	G4UIcmdWithABool * Set_gamma_mode_cmd;
};
#endif
