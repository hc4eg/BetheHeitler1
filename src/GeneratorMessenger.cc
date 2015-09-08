#include "GeneratorMessenger.hh"

GeneratorMessenger::GeneratorMessenger( PrimaryGeneratorAction* gen)
{
	generator = gen;
	Dir1 = new G4UIdirectory("/BH/gun/");
	Dir1->SetGuidance("Change BH Simulation gun Parameters");

	Set_central_energy_cmd = new G4UIcmdWithADoubleAndUnit("/BH/gun/set_central_energy", this);
	Set_central_energy_cmd->SetGuidance("Set the central energy.");
	Set_central_energy_cmd->SetGuidance("Usage: set_central_energy energy [unit]");
	Set_central_energy_cmd->SetDefaultUnit("MeV");
	Set_central_energy_cmd->SetUnitCandidates("MeV keV");
	Set_central_energy_cmd->SetParameterName("central_energy", false, true);
	Set_central_energy_cmd->AvailableForStates(G4State_Idle);

	Set_radius_max_cmd = new G4UIcmdWithADoubleAndUnit("/BH/gun/set_radius_max", this);
	Set_radius_max_cmd->SetGuidance("Set maximum position radius.");
	Set_radius_max_cmd->SetGuidance("Usage: set_radius_max radius_max [unit]");
	Set_radius_max_cmd->SetDefaultUnit("cm");
	Set_radius_max_cmd->SetUnitCandidates("mm cm m");
	Set_radius_max_cmd->SetParameterName("radius_max", false, true);
	Set_radius_max_cmd->AvailableForStates(G4State_Idle);

	Set_x_min_cmd = new G4UIcmdWithADoubleAndUnit("/BH/gun/set_x_min", this);
	Set_x_min_cmd->SetGuidance("Set minimum x position.");
	Set_x_min_cmd->SetGuidance("Usage: set_x_min x_min [unit]");
	Set_x_min_cmd->SetDefaultUnit("cm");
	Set_x_min_cmd->SetUnitCandidates("mm cm m");
	Set_x_min_cmd->SetParameterName("x_min", false, true);
	Set_x_min_cmd->AvailableForStates(G4State_Idle);

	Set_x_max_cmd = new G4UIcmdWithADoubleAndUnit("/BH/gun/set_x_max", this);
	Set_x_max_cmd->SetGuidance("Set maximum x position.");
	Set_x_max_cmd->SetGuidance("Usage: set_x_max x_max [unit]");
	Set_x_max_cmd->SetDefaultUnit("cm");
	Set_x_max_cmd->SetUnitCandidates("mm cm m");
	Set_x_max_cmd->SetParameterName("x_max", false, true);
	Set_x_max_cmd->AvailableForStates(G4State_Idle);

	Set_y_min_cmd = new G4UIcmdWithADoubleAndUnit("/BH/gun/set_y_min", this);
	Set_y_min_cmd->SetGuidance("Set minimum x position.");
	Set_y_min_cmd->SetGuidance("Usage: set_y_min y_min [unit]");
	Set_y_min_cmd->SetDefaultUnit("cm");
	Set_y_min_cmd->SetUnitCandidates("mm cm m");
	Set_y_min_cmd->SetParameterName("y_min", false, true);
	Set_y_min_cmd->AvailableForStates(G4State_Idle);

	Set_y_max_cmd = new G4UIcmdWithADoubleAndUnit("/BH/gun/set_y_max", this);
	Set_y_max_cmd->SetGuidance("Set maximum x position.");
	Set_y_max_cmd->SetGuidance("Usage: set_y_max y_max [unit]");
	Set_y_max_cmd->SetDefaultUnit("cm");
	Set_y_max_cmd->SetUnitCandidates("mm cm m");
	Set_y_max_cmd->SetParameterName("y_max", false, true);
	Set_y_max_cmd->AvailableForStates(G4State_Idle);

	Set_angle_max_cmd = new G4UIcmdWithADoubleAndUnit("/BH/gun/set_angle_max", this);
	Set_angle_max_cmd->SetGuidance("Set maximum angle.");
	Set_angle_max_cmd->SetGuidance("Usage: set_angle_max angle_max [unit]");
	Set_angle_max_cmd->SetDefaultUnit("mrad");
	Set_angle_max_cmd->SetUnitCandidates("mrad rad deg");
	Set_angle_max_cmd->SetParameterName("angle_max", false, true);
	Set_angle_max_cmd->AvailableForStates(G4State_Idle);

	Set_theta_min_cmd = new G4UIcmdWithADoubleAndUnit("/BH/gun/set_theta_min", this);
	Set_theta_min_cmd->SetGuidance("Set minimum theta angle.");
	Set_theta_min_cmd->SetGuidance("Usage: set_theta_min theta_min [unit]");
	Set_theta_min_cmd->SetDefaultUnit("mrad");
	Set_theta_min_cmd->SetUnitCandidates("mrad rad deg");
	Set_theta_min_cmd->SetParameterName("theta_min", false, true);
	Set_theta_min_cmd->AvailableForStates(G4State_Idle);

	Set_theta_max_cmd = new G4UIcmdWithADoubleAndUnit("/BH/gun/set_theta_max", this);
	Set_theta_max_cmd->SetGuidance("Set maximum theta angle.");
	Set_theta_max_cmd->SetGuidance("Usage: set_theta_max theta_max [unit]");
	Set_theta_max_cmd->SetDefaultUnit("mrad");
	Set_theta_max_cmd->SetUnitCandidates("mrad rad deg");
	Set_theta_max_cmd->SetParameterName("theta_max", false, true);
	Set_theta_max_cmd->AvailableForStates(G4State_Idle);

	Set_phi_min_cmd = new G4UIcmdWithADoubleAndUnit("/BH/gun/set_phi_min", this);
	Set_phi_min_cmd->SetGuidance("Set minimum phi angle.");
	Set_phi_min_cmd->SetGuidance("Usage: set_phi_min phi_min [unit]");
	Set_phi_min_cmd->SetDefaultUnit("mrad");
	Set_phi_min_cmd->SetUnitCandidates("mrad rad deg");
	Set_phi_min_cmd->SetParameterName("phi_min", false, true);
	Set_phi_min_cmd->AvailableForStates(G4State_Idle);

	Set_phi_max_cmd = new G4UIcmdWithADoubleAndUnit("/BH/gun/set_phi_max", this);
	Set_phi_max_cmd->SetGuidance("Set maximum phi angle.");
	Set_phi_max_cmd->SetGuidance("Usage: set_phi_max phi_max [unit]");
	Set_phi_max_cmd->SetDefaultUnit("mrad");
	Set_phi_max_cmd->SetUnitCandidates("mrad rad deg");
	Set_phi_max_cmd->SetParameterName("phi_max", false, true);
	Set_phi_max_cmd->AvailableForStates(G4State_Idle);

	Set_delta_min_cmd = new G4UIcmdWithADouble("/BH/gun/set_delta_min", this);
	Set_delta_min_cmd->SetGuidance("Set minimum % from central energy.");
	Set_delta_min_cmd->SetGuidance("Usage: set_delta_min percent");
	Set_delta_min_cmd->SetParameterName("delta_min", false, false);
	Set_delta_min_cmd->AvailableForStates(G4State_Idle);

	Set_delta_max_cmd = new G4UIcmdWithADouble("/BH/gun/set_delta_max", this);
	Set_delta_max_cmd->SetGuidance("Set maximum % from central energy.");
	Set_delta_max_cmd->SetGuidance("Usage: set_delta_max percent");
	Set_delta_max_cmd->SetParameterName("delta_max", false, false);
	Set_delta_max_cmd->AvailableForStates(G4State_Idle);

	Set_gamma_mode_cmd = new G4UIcmdWithABool("/BH/gun/set_gamma_mode", this);
	Set_gamma_mode_cmd->SetGuidance("Set Gamma Mode: Particles will be generated");
	Set_gamma_mode_cmd->SetGuidance("upstream of target.");
	Set_gamma_mode_cmd->AvailableForStates(G4State_Idle);

	PrintCmd = new G4UIcmdWithoutParameter("/BH/gun/print", this);
	PrintCmd->SetGuidance("Print current gun parameters.");
	PrintCmd->AvailableForStates(G4State_Idle);

	//G4cout << "End Generator Messenger Constructor" << G4endl;
}

GeneratorMessenger::~GeneratorMessenger()
{
	//delete PrintCmd;
	delete Set_central_energy_cmd;
	delete Set_delta_min_cmd;
	delete Set_delta_max_cmd;
	delete Set_radius_max_cmd;
	delete Set_x_min_cmd;
	delete Set_x_max_cmd;
	delete Set_y_min_cmd;
	delete Set_y_max_cmd;
	delete Set_theta_min_cmd;
	delete Set_theta_max_cmd;
	delete Set_phi_min_cmd;
	delete Set_phi_max_cmd;
	delete Set_angle_max_cmd;
	delete Set_gamma_mode_cmd;
	delete Dir1;
}

void
GeneratorMessenger::SetNewValue( G4UIcommand *command, G4String newValue)
{
	if(command == PrintCmd) { generator->PrintParameters(); }
	else if(command == Set_gamma_mode_cmd)
		{ generator->Set_gamma_mode( Set_gamma_mode_cmd->GetNewBoolValue(newValue)); }
	else if(command == Set_central_energy_cmd)
		{ generator->Set_central_energy( Set_central_energy_cmd->GetNewDoubleValue(newValue)); }
	else if(command == Set_delta_min_cmd)
		{ generator->Set_delta_min( Set_delta_min_cmd->GetNewDoubleValue(newValue)); }
	else if(command == Set_delta_max_cmd)
		{ generator->Set_delta_max( Set_delta_max_cmd->GetNewDoubleValue(newValue)); }
	else if(command == Set_radius_max_cmd)
		{ generator->Set_radius_max( Set_radius_max_cmd->GetNewDoubleValue(newValue)); }
	else if(command == Set_x_min_cmd)
		{ generator->Set_x_min( Set_x_min_cmd->GetNewDoubleValue(newValue)); }
	else if(command == Set_x_max_cmd)
		{ generator->Set_x_max( Set_x_max_cmd->GetNewDoubleValue(newValue)); }
	else if(command == Set_y_min_cmd)
		{ generator->Set_y_min( Set_y_min_cmd->GetNewDoubleValue(newValue)); }
	else if(command == Set_y_max_cmd)
		{ generator->Set_y_max( Set_y_max_cmd->GetNewDoubleValue(newValue)); }
	else if(command == Set_theta_min_cmd)
		{ generator->Set_theta_min( Set_theta_min_cmd->GetNewDoubleValue(newValue)); }
	else if(command == Set_theta_max_cmd)
		{ generator->Set_theta_max( Set_theta_max_cmd->GetNewDoubleValue(newValue)); }
	else if(command == Set_phi_min_cmd)
		{ generator->Set_phi_min( Set_phi_min_cmd->GetNewDoubleValue(newValue)); }
	else if(command == Set_phi_max_cmd)
		{ generator->Set_phi_max( Set_phi_max_cmd->GetNewDoubleValue(newValue)); }
	else if(command == Set_angle_max_cmd)
		{ generator->Set_angle_max( Set_angle_max_cmd->GetNewDoubleValue(newValue)); }
}

//G4String GeneratorMessenger::GetCurrentValue(G4UIcommand *command)
//{
//	return "Unused for now!";
//}
	
