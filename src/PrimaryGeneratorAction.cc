//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "PrimaryGeneratorAction.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "Randomize.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PrimaryGeneratorAction::PrimaryGeneratorAction( DetectorConstruction* DC)
:Detector(DC)
{
  G4int n_particle = 1;
  particleGun  = new G4ParticleGun(n_particle);

  pOutputFile = OutputFile::GetOutputFilePointer();
  
  //create a messenger for this class
  gunMessenger = new GeneratorMessenger(this);

  // default particle kinematic

  G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
  G4String particleName;
  G4ParticleDefinition * p_electron = particleTable->FindParticle(particleName="e-");

  //some default values defining source
  central_energy = 30.0*MeV;
  //target_position = -42.26*cm;
  // get from geometry
  target_position = -Detector->GetTargetDistance();
  target_thickness = Detector->GetTargetThickness();
  delta_min = 0.; // percent from central_energy
  delta_max = 0.;
  radius_max = 0.65*cm;
  x_min = 0.*cm;
  x_max = 0.*cm;
  y_min = 0.*cm;
  y_max = 0.*cm;
  angle_max = 175.*mrad;
  theta_min = 0*mrad;
  theta_max = 0*mrad;
  phi_min = 0*mrad;
  phi_max = 0*mrad;

  particleGun->SetParticleDefinition(p_electron);
  particleGun->SetParticleMomentumDirection(G4ThreeVector(1.,0.,0.));
  particleGun->SetParticleEnergy(central_energy);
  particleGun->SetParticlePosition(G4ThreeVector(target_position, 0.*cm, 0.*cm));

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PrimaryGeneratorAction::~PrimaryGeneratorAction()
{
  delete particleGun;
  delete pOutputFile;
  delete gunMessenger;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{
	G4double x_in, y_in, radius, energy_in, theta_in, phi_in;
	G4double delta_in, v_x, v_y, v_z, targ_in;
  //this function is called at the begining of event
  //
  // get from geometry  - may have changed after constructor
  target_position = -Detector->GetTargetDistance();
  target_thickness = Detector->GetTargetThickness();

   	//G4cout << "Generate Primaries called." << G4endl;
  	// choose position
  	do	{
		if(x_max == x_min) { x_in = x_min; }
		else { x_in =  CLHEP::RandFlat::shoot(x_min, x_max); }
		if(y_max == y_min) { y_in = y_min; }
		else { y_in =  CLHEP::RandFlat::shoot(y_min, y_max); }
		radius = sqrt(x_in*x_in +y_in*y_in);
		} while (radius > radius_max);
	// choose angle
	do	{
		if(theta_max == theta_min) { theta_in = theta_min; }
		else { theta_in =  CLHEP::RandFlat::shoot(theta_min, theta_max); }
		if(phi_max == phi_min) { phi_in = phi_min; }
		else { phi_in =  CLHEP::RandFlat::shoot(phi_min, phi_max); }
		radius = sqrt(theta_in*theta_in +phi_in*phi_in);
		} while (radius > angle_max);
	// choose delta
	if(delta_max == delta_min) { delta_in = delta_min; }
	else { delta_in =  CLHEP::RandFlat::shoot(delta_min, delta_max); }
	energy_in = central_energy*(1. + delta_in/100.);
	particleGun->SetParticleEnergy(energy_in);
	if(gamma_mode)
		{
		targ_in = target_position - Detector->GetTargetContainerThick()/2.;
		}
	else
		{
		// choose random depth within target
		targ_in =  target_position + CLHEP::RandFlat::shoot(-target_thickness/2., target_thickness/2.);
		}
	particleGun->SetParticlePosition(G4ThreeVector(targ_in, x_in, y_in));
	v_x = 1.*cm;
	v_y = tan(theta_in)*v_x;
	v_z = tan(phi_in)*v_x;
 	particleGun->SetParticleMomentumDirection(G4ThreeVector(v_x, v_y, v_z));

	pOutputFile->Set_energy_i(energy_in);
	pOutputFile->Set_delta_i(delta_in);
	pOutputFile->Set_x_i(x_in);
	pOutputFile->Set_y_i(y_in);
	pOutputFile->Set_theta_i(theta_in);
	pOutputFile->Set_phi_i(phi_in);
	
  //G4cout << "Generate Primary...." << G4endl;
  particleGun->GeneratePrimaryVertex(anEvent);
  //G4cout << "Generated Primary...." << G4endl;

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void
PrimaryGeneratorAction::PrintParameters()
	{
	if(gamma_mode)
		{
		G4cout << "Gamma mode:" << G4endl;
		G4cout << "Photon Energy = " << central_energy/MeV << " MeV" << G4endl;
		}
	else 	{
		G4cout << "Normal mode: " << G4endl;
		G4cout << "Central Energy = " << central_energy/MeV << " MeV" << G4endl;
		}
	G4String particle_name = particleGun->GetParticleDefinition()->GetParticleName();
	G4cout << "Particle       = " << particle_name << G4endl;
	G4cout << "Delta min      = " << delta_min << " %" << G4endl;
	G4cout << "Delta max      = " << delta_max << " %" << G4endl;
	G4cout << "Radius max     = " << radius_max/cm << " cm" << G4endl;
	G4cout << "X min          = " << x_min/cm << " cm" << G4endl;
	G4cout << "X max          = " << x_max/cm << " cm" << G4endl;
	G4cout << "Y min          = " << y_min/cm << " cm" << G4endl;
	G4cout << "Y max          = " << y_max/cm << " cm" << G4endl;
	G4cout << "Angle max      = " << angle_max/mrad << " mrad" << G4endl;
	G4cout << "Theta min      = " << theta_min/mrad << " mrad" << G4endl;
	G4cout << "Theta max      = " << theta_max/mrad << " mrad" << G4endl;
	G4cout << "Phi min        = " << phi_min/mrad << " mrad" << G4endl;
	G4cout << "Phi max        = " << phi_max/mrad << " mrad" << G4endl;
	}
void
PrimaryGeneratorAction::PrintParametersToFile()
	{
	char line[100];
	pOutputFile->WriteComment("----Primary Generator Parameters----");
	if(gamma_mode)
		{
		sprintf(line, "Gamma mode:");
		pOutputFile->WriteComment(G4String(line));
		sprintf(line, "Photon Energy = %.3g MeV", central_energy/MeV);
		pOutputFile->WriteComment(G4String(line));
		}
	else 	{
		sprintf(line, "Normal mode:");
		pOutputFile->WriteComment(G4String(line));
		sprintf(line, "Central Energy = %.3g MeV", central_energy/MeV);
		pOutputFile->WriteComment(G4String(line));
		}
	G4String particle_name = particleGun->GetParticleDefinition()->GetParticleName();
	G4String part = "Particle: "+particle_name;
	pOutputFile->WriteComment(part);
	sprintf(line, "Delta min      = %.3g %%", delta_min);
	pOutputFile->WriteComment(G4String(line));
	sprintf(line, "Delta max      = %.3g %%", delta_max);
	pOutputFile->WriteComment(G4String(line));
	sprintf(line, "Radius max     = %.3g cm", radius_max/cm);
	pOutputFile->WriteComment(G4String(line));
	sprintf(line, "X min          = %.3g cm", x_min/cm);
	pOutputFile->WriteComment(G4String(line));
	sprintf(line, "X max          = %.3g cm", x_max/cm);
	pOutputFile->WriteComment(G4String(line));
	sprintf(line, "Y min          = %.3g cm", y_min/cm);
	pOutputFile->WriteComment(G4String(line));
	sprintf(line, "Y max          = %.3g cm", y_max/cm);
	pOutputFile->WriteComment(G4String(line));
	sprintf(line, "Angle max      = %.3g mrad", angle_max/mrad);
	pOutputFile->WriteComment(G4String(line));
	sprintf(line, "Theta min      = %.3g mrad", theta_min/mrad);
	pOutputFile->WriteComment(G4String(line));
	sprintf(line, "Theta max      = %.3g mrad", theta_max/mrad);
	pOutputFile->WriteComment(G4String(line));
	sprintf(line, "Phi min        = %.3g mrad", phi_min/mrad);
	pOutputFile->WriteComment(G4String(line));
	sprintf(line, "Phi max        = %.3g mrad", phi_max/mrad);
	pOutputFile->WriteComment(G4String(line));
	}
