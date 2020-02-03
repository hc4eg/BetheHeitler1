//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "PrimaryGeneratorAction.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "Randomize.hh"
#include "math.h"
#include <cmath>
#define PI 3.14159265
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

inline double randfloat(){
	double r = random()/ double (RAND_MAX);
	return r;
}

PrimaryGeneratorAction::PrimaryGeneratorAction( DetectorConstruction* DC)
:Detector(DC)
{
  pair_mode = true;
  //pair_mode = false;

  /*
  cerr << "Pair mode = " << (pair_mode==true?"True":"False") << ", continue?"<< endl;
  string s_pair_mode;
  cin >> s_pair_mode;
  */

  G4int n_particle;
  if(pair_mode)
  	n_particle = 2;
  else
	n_particle = 1;

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
  target_position = -42.26*cm;
  // get from geometry
  //target_position = -Detector->GetTargetDistance();
  target_thickness = Detector->GetTargetThickness();
  delta_min = 0.; // percent from central_energy
  delta_max = 0.;

  //Here primary particles are generated within rectangle bounded by
  // x = x_min, x = x_max, y = y_min, y = y_max.
  // with additional condition sqrt(x^2+y^2)<radius_max (circle with radius radius_max)
  /*
  x_min = 0.*cm;
  x_max = 0.*cm;
  y_min = 0.*cm;
  y_max = 0.*cm;
  */

  x_min = -0.5*cm;
  x_max = 0.5*cm;
  y_min = -0.5*cm;
  y_max = 0.5*cm;
  //radius_max = 0.65*cm;
  radius_max = 0.5*cm;

  angle_max = 175.*mrad;
  theta_min = 0*mrad;
  theta_max = 0*mrad;
  phi_min = 0*mrad;
  phi_max = 0*mrad;

  particleGun->SetParticleDefinition(p_electron);
  particleGun->SetParticleMomentumDirection(G4ThreeVector(1.,0.,0.));
  particleGun->SetParticleEnergy(central_energy);
  particleGun->SetParticlePosition(G4ThreeVector(target_position, 0.*cm, 0.*cm));

  E_gamma = 60*MeV;
  c = CLHEP::c_light;
  Me = CLHEP::electron_mass_c2;

  if(pair_mode){
	string filenumber;
	cerr << "Please input filenumber: (Format XXX)" << endl;
	cin >> filenumber;
	OpenFile(filenumber);

	maxline = 0;
	//  infile.seekg(0);
	while(!infile.eof()){
		infile.ignore(200,'\n');
		maxline++;
	//if(maxline%100 == 0 ) cerr << "Read line: " << maxline << endl;
	}

	maxline--;
	cerr << "Data file has Total line " << maxline << endl;
	infile.seekg(0);
  }

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PrimaryGeneratorAction::~PrimaryGeneratorAction()
{
  delete particleGun;
  delete pOutputFile;
  delete gunMessenger;
  infile.close();
}



//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void PrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{
	  ////////////////////////////////////////////////
	  // In branch LinearEnergy:
	//1. Energy uniformly  distributed from E_min to E_max
	//2. Based on energy primary vertex position follows energy by:
	// 		x = xmin + (E - Emin)*(xmax - xmin)/(Emax - Emin)
	static const G4double inch = 2.54 * cm;
	G4double x_out[2], y_out[2], z_out[2];

	G4double E_min = 15.*MeV;
	G4double E_max = 45.*MeV;

	G4double X_min = -40.*inch;
	G4double X_max = 40.*inch;

	G4double E[2];
	E[1] = CLHEP::RandFlat::shoot(E_min, E_max);
	E[0] = 60*MeV-E[1];

	for(int i = 0; i < 2; i++){
		x_out[i] = X_min + (E[i] - E_min)*(X_max - X_min)/(E_max - E_min);
		//x_out[i] = 0;
		y_out[i] = pow(-1., (double)i) * 20.0 * inch;
		z_out[i] = 0.;
		cerr << "Input particle " << i << " position: (" << x_out[i]/cm
		<< ", " << y_out[i]/cm << "," << z_out[i]/cm << ")" << endl;
		cerr << "With energy " << E[i] << endl;
	}


	//if (pair_mode)
	if(1)
		{
		  //ConvertNext();
		  G4int CurrentEvent = anEvent->GetEventID();
		  G4ParticleTable* ParticleTable = G4ParticleTable::GetParticleTable();
		  G4String ParticleName;

		  // Note Theta, Phi are projection angles than spherical angles
		  // FIXME: units for Theta and Phi?
		  // for output: index 0: positron, index 1:eletron.

		  // Positron generation
		  G4ParticleDefinition * P_positron = ParticleTable->FindParticle(ParticleName="e+");
		  particleGun->SetParticleDefinition(P_positron);
		  particleGun->SetParticlePosition(G4ThreeVector(x_out[0], y_out[0], z_out[0]));
		  particleGun->SetParticleEnergy(E[0]);
		  particleGun->SetParticleMomentumDirection(G4ThreeVector(0., 1., 0.));
			pOutputFile->Set_energy_i(0, E[0]);
			pOutputFile->Set_x_i(0, 0.);
			pOutputFile->Set_y_i(0, 0.);
			pOutputFile->Set_phi_i(0, 0.);
			pOutputFile->Set_theta_i(0, 0.);
		  particleGun->GeneratePrimaryVertex(anEvent);

		  // Electron generation
		  G4ParticleDefinition * P_electron = ParticleTable->FindParticle(ParticleName="e-");
		  particleGun->SetParticleDefinition(P_positron);
		  particleGun->SetParticlePosition(G4ThreeVector(x_out[1], y_out[1], z_out[1]));
		  particleGun->SetParticleEnergy(E[1]);
		  particleGun->SetParticleMomentumDirection(G4ThreeVector(0., -1., 0.));
			pOutputFile->Set_energy_i(1, E[1]);
			pOutputFile->Set_x_i(1, 0.);
			pOutputFile->Set_y_i(1, 0.);
			pOutputFile->Set_phi_i(1, 0.);
			//pOutputFile->Set_theta_i(0, 180.);
			pOutputFile->Set_theta_i(1, 180. * deg);
		  particleGun->GeneratePrimaryVertex(anEvent);
		}

}



//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
// OpenFile: open file that contains generator pair data
void PrimaryGeneratorAction::OpenFile(string filenumber){
	string filename;
	filename = "/home/haoyu/BHsim/BetheHeitler-1.3/src/events.run"+filenumber+".dat";
	cerr << "Primary vertex data filename is " << filename  << endl;
	infile.open(filename.c_str());
	if(!infile.is_open()) cerr << "File open failed."  << endl;
}

// Initialize and Convert lines of datafile parameters to primary vertex parameters.
int PrimaryGeneratorAction::Convert(int linenumber){
	if(!infile.is_open()) { cerr << "Primary vertex data file not open." << endl; Clear(); return -1;}
	if(linenumber == 0) { cerr <<"Line number 0." << endl; Clear(); return -1;}
	infile.seekg(0);
	if(linenumber > maxline) { cerr << "Line number exceed total line number: " << maxline << endl; Clear(); return -1; }
	
	for(int i = 0; i < linenumber-1; i++) { infile.ignore(200, '\n');}
	ConvertNext();
	return 0;
}

// Code for sphereical angles stored in pair data file
/*
int PrimaryGeneratorAction::ConvertNext(){
	if(!infile.eof()){
		infile >> KEe >> KEp >> Thetae >> Phie >> Thetap >> Phip;
		infile.ignore(200, '\n');
		KEe *=MeV;
		KEp *=MeV;
		Thetae *= deg;
		Thetap *= deg;
		Phie *= deg;
		Phip *= deg;

		// From the way Py computed: 
		// phi = 0, left , phi = 180 deg , rihgt
	        // Clockwise looking towards x direction
		Pe = sqrt((KEe+Me)*(KEe+Me) - Me*Me);
		Pp = sqrt((KEp+Me)*(KEp+Me) - Me*Me);
		Pex = Pe*cos(Thetae);
		Ppx = Pp*cos(Thetap);
		Pey = Pe*sin(Thetae)*cos(Phie);	
		Ppy = Pp*sin(Thetap)*cos(Phip);	
		Pez = Pe*sin(Thetae)*sin(Phie);
		Ppz = Pp*sin(Thetap)*sin(Phip);
	}
	else { cerr << "End of file, no conversion. Return -1" << endl; return -1;}
	return 0;
}
*/

// Code for projection angles stored in pair data file
int PrimaryGeneratorAction::ConvertNext(){
	if(!infile.eof()){
		infile >> KEe >> KEp >> Thetae >> Phie >> Thetap >> Phip;
		infile.ignore(200, '\n');

		// KE stored in pair data file in units MeV
		KEe *=MeV;
		KEp *=MeV;

		// Theta,Phi stored in data file in unit degree
		// deg = PI/180; so Theta(in deg)*deg turns Theta(in rad) for computation
		Thetae *= deg;
		Thetap *= deg;
		Phie *= deg;
		Phip *= deg;

		// From the way Py computed: 
		// phi = 0, left , phi = 180 deg , right
	        // Clockwise looking towards x direction
		Pe = sqrt((KEe+Me)*(KEe+Me) - Me*Me);
		Pp = sqrt((KEp+Me)*(KEp+Me) - Me*Me);


		// Check if data in EventAction.cc will follow pairs in primary generator.
		/*
		bool b_s_dir = 0;
		if(b_s_dir){
			Thetae = 6.*deg;
			Thetap = -6.*deg;
			Phie = 30*deg;
			Phip = -20*deg;
		}
		*/

		// Covert to projection angles:

		// x is beam direction
		// y is horrizontal, left positive
		// z is vertical, up positive
		// Makes (x,y,z) right hand
		Pex = Pe/sqrt(1+tan(Thetae)*tan(Thetae)+tan(Phie)*tan(Phie));
		Pey = Pex*tan(Thetae);
		Pez = Pex*tan(Phie);

		Ppx = Pp/sqrt(1+tan(Thetap)*tan(Thetap)+tan(Phip)*tan(Phip));
		Ppy = Ppx*tan(Thetap);
		Ppz = Ppx*tan(Phip);
	}
	else { cerr << "End of file, no conversion. Return -1" << endl; return -1;}
	return 0;
}

// PrintVertex(): Print e(electron) and p(ositron)'s Theta, Phi, P , Px, Py, Pz
void PrimaryGeneratorAction::PrintVertex(){
	cout << "KEe " << KEe/MeV << " , Thetae " << Thetae/deg << ", Phie " << Phie/deg << " , Pe " << Pe/MeV << " , Pex " << Pex/MeV << " , Pey " << Pey/MeV << " , Pez " << Pez/MeV << endl; 
	cout << "KEp " << KEp/MeV << " , Thetap " << Thetap/deg << ", Phip " << Phip/deg << " , Pp " << Pp/MeV << " , Ppx " << Ppx/MeV << " , Ppy " << Ppy/MeV << " , Ppz " << Ppz/MeV << endl; 
}

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
