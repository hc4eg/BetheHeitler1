//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "DetectorMessenger.hh"

#include "DetectorConstruction.hh"
#include "G4UIdirectory.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWithoutParameter.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorMessenger::DetectorMessenger( DetectorConstruction* Det)
:Detector(Det)
{ 
  BHDir = new G4UIdirectory("/BH/");
  BHDir->SetGuidance("UI commands of this BetheHeitler Magnet");
  
  detDir = new G4UIdirectory("/BH/detector/");
  detDir->SetGuidance("detector control");

	// Get materials list
      	G4String matlist;
	const G4MaterialTable * matable = G4Material::GetMaterialTable();
	for(size_t i=0; i < G4Material::GetNumberOfMaterials(); i++)
		{
		matlist += (*matable)[i]->GetName();
		matlist += " ";
		}	
	//G4cout << "Detector Messenger Constructor...." << G4endl;
	//G4cout << matlist << G4endl;
 
  Update_cmd = new G4UIcmdWithoutParameter("/BH/detector/update",this);
  Update_cmd->SetGuidance("Update geometry.");
  Update_cmd->SetGuidance("This command MUST be applied before \"beamOn\" ");
  Update_cmd->SetGuidance("if you changed geometrical value(s).");
  Update_cmd->AvailableForStates(G4State_Idle);
      
  Print_cmd = new G4UIcmdWithoutParameter("/BH/detector/print",this);
  Print_cmd->SetGuidance("Print changeable geometry parameters.");
  Print_cmd->AvailableForStates(G4State_Idle);

  SetVerbose_cmd = new G4UIcmdWithABool("/BH/detector/verbose",this);
  SetVerbose_cmd->SetGuidance("Set verbose output during detector construction.");
  SetVerbose_cmd->SetGuidance("Usage: verbose true/false.");
  SetVerbose_cmd->SetParameterName("verbose", false, false);
  SetVerbose_cmd->AvailableForStates(G4State_Idle);
      
  Dump_cmd = new G4UIcmdWithAnInteger("/BH/detector/dump",this);
  Dump_cmd->SetGuidance("Print information about geometry.");
  Dump_cmd->SetGuidance("Usage: dump [depth].");
  Dump_cmd->SetGuidance("Default depth = 2");
  Dump_cmd->SetParameterName("depth", true, false);
  Dump_cmd->SetDefaultValue(2);
  Dump_cmd->AvailableForStates(G4State_Idle);

  SetB0_cmd = new G4UIcmdWithADoubleAndUnit("/BH/detector/set_B0", this);
  SetB0_cmd->SetGuidance("Set central magnetic field.");
  SetB0_cmd->SetGuidance("Usage: set_B0 B0 [unit]");
  SetB0_cmd->SetDefaultUnit("tesla");
  SetB0_cmd->SetUnitCandidates("tesla gauss kilogauss");
  SetB0_cmd->SetParameterName("B0", false, true);
  SetB0_cmd->AvailableForStates(G4State_Idle);

  SetTargetMat_cmd = new G4UIcmdWithAString("/BH/detector/set_target_material", this);
  SetTargetMat_cmd->SetGuidance("Select Target Material.");
  SetTargetMat_cmd->SetGuidance("Usage: set_target_material material");
  SetTargetMat_cmd->SetGuidance("MUST execute /BH/detector/update for changes to take effect.");
  SetTargetMat_cmd->SetParameterName("material", false, false);
  SetTargetMat_cmd->SetCandidates(matlist);
  SetTargetMat_cmd->AvailableForStates(G4State_Idle);

  SetTargetThick_cmd = new G4UIcmdWithADoubleAndUnit("/BH/detector/set_target_thickness", this);
  SetTargetThick_cmd->SetGuidance("Set thickness of target.");
  SetTargetThick_cmd->SetGuidance("Usage: set_target_thickness thick [unit]");
  SetTargetThick_cmd->SetGuidance("MUST execute /BH/detector/update for changes to take effect.");
  SetTargetThick_cmd->SetDefaultUnit("mm");
  SetTargetThick_cmd->SetUnitCandidates("mm cm microm");
  SetTargetThick_cmd->SetParameterName("thick", false, true);
  SetTargetThick_cmd->AvailableForStates(G4State_Idle);

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorMessenger::~DetectorMessenger()
{
	delete Update_cmd;
	delete Print_cmd;
	delete Dump_cmd;
	delete SetB0_cmd;
	delete detDir;
	delete BHDir;  
	delete SetTargetMat_cmd;
	delete SetTargetThick_cmd;
	delete SetVerbose_cmd;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorMessenger::SetNewValue(G4UIcommand* command,G4String newValue)
{ 
 
  if( command == Update_cmd )
	{
	Detector->UpdateGeometry();
	}
  else if(command == Dump_cmd)
	{
	G4int depth = Dump_cmd->GetNewIntValue(newValue);
	Detector->DumpGeometricalTree(Detector->GetPhysWorld(), depth);
	}
  else if(command == Print_cmd)
	{
	Detector->Print();
	}
  else if(command == SetVerbose_cmd)
	{
	Detector->SetVerbose( SetVerbose_cmd->GetNewBoolValue(newValue));
	}
  else if(command == SetB0_cmd)
	{
	Detector->SetB0( SetB0_cmd->GetNewDoubleValue(newValue));
	}
  else if(command == SetTargetThick_cmd)
	{
	Detector->SetTargetThickness( SetTargetThick_cmd->GetNewDoubleValue(newValue));
	}
  else if(command == SetTargetMat_cmd)
	{
	Detector->SetTargetMaterial( newValue);
	}


}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
