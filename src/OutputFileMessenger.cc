#include "OutputFileMessenger.hh"

OutputFileMessenger::OutputFileMessenger( OutputFile* out)
{
	fout = out;
	Dir1 = new G4UIdirectory("/BH/outputfile/");
	Dir1->SetGuidance("Change BH Output File Parameters");

	Set_use_monitor_cmd = new G4UIcmdWithABool("/BH/outputfile/set_use_monitor", this);
	Set_use_monitor_cmd->SetGuidance("If true: Record and included monitor information in output file.");
	Set_use_monitor_cmd->SetGuidance("If false: Do not record monitor information.");
	Set_use_monitor_cmd->SetGuidance("Usage: set_use_monitor <bool>");
	Set_use_monitor_cmd->AvailableForStates(G4State_Idle);

	Set_unique_cmd = new G4UIcmdWithABool("/BH/outputfile/set_use_unique_filename", this);
	Set_unique_cmd->SetGuidance("Output filename is <prefix>.<runno>.dat");
	Set_unique_cmd->SetGuidance("If false: runno = G4 run number.");
	Set_unique_cmd->SetGuidance("If true: runno = run number read from file <prefix>.runno.");
	Set_unique_cmd->SetGuidance("Usage: set_unique_filename <bool>");
	Set_unique_cmd->AvailableForStates(G4State_Idle);

	Set_header_frequency_cmd = new G4UIcmdWithAnInteger("/BH/outputfile/set_header_frequency", this);
	Set_header_frequency_cmd->SetGuidance("Set the number of output files lines between when a header is printed.");
	Set_header_frequency_cmd->SetGuidance("Usage: set_header_frequency <lines>");
	Set_header_frequency_cmd->SetGuidance("lines = -1 means no headers after the first (default)");
	Set_header_frequency_cmd->AvailableForStates(G4State_Idle);

	Set_prefix_cmd = new G4UIcmdWithAString("/BH/outputfile/set_prefix", this);
	Set_prefix_cmd->SetGuidance("Set the prefix for the output file filename.");
	Set_prefix_cmd->SetGuidance("Usage: set_prefix <prefix>");
	Set_prefix_cmd->AvailableForStates(G4State_Idle);

	PrintCmd = new G4UIcmdWithoutParameter("/BH/outputfile/print", this);
	PrintCmd->SetGuidance("Print current output file parameters.");
	PrintCmd->AvailableForStates(G4State_Idle);

	//G4cout << "End Generator Messenger Constructor" << G4endl;
}

OutputFileMessenger::~OutputFileMessenger()
{
	delete PrintCmd;
	delete Set_prefix_cmd;
	delete Set_header_frequency_cmd;
	delete Set_unique_cmd;
	delete Set_use_monitor_cmd;
	delete Dir1;
}

void
OutputFileMessenger::SetNewValue( G4UIcommand *command, G4String newValue)
{
	if(command == PrintCmd) { fout->PrintParameters(); }
	else if(command == Set_header_frequency_cmd)
		{ fout->Set_output_lines( Set_header_frequency_cmd->GetNewIntValue(newValue)); }
	else if(command == Set_unique_cmd)
		{ fout->Set_use_unique_filename( Set_unique_cmd->GetNewBoolValue(newValue)); }
	else if(command == Set_use_monitor_cmd)
		{ fout->SetUseMonitor( Set_use_monitor_cmd->GetNewBoolValue(newValue)); }
	else if(command == Set_prefix_cmd)
		{ fout->Set_filename_prefix( newValue); }
}

//G4String GeneratorMessenger::GetCurrentValue(G4UIcommand *command)
//{
//	return "Unused for now!";
//}
	
