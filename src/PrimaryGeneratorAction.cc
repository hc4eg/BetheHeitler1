//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "PrimaryGeneratorAction.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "Randomize.hh"
#include "math.h"

#define PI 3.14159265
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

  //Here primary particles are generated within rectangle bounded by
  // x = x_min, x = x_max, y = y_min, y = y_max.
  // with additional condition sqrt(x^2+y^2)<radius_max (circle with radius radius_max)

  //x_min = 0.*cm;
  //x_max = 0.*cm;
  //y_min = 0.*cm;
  //y_max = 0.*cm;

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

  E = 60*MeV;
  c = CLHEP::c_light;
  Me = CLHEP::electron_mass_c2;

  string filenumber;
  cerr << "Please input filenumber: (Format XXX)" << endl;
  cin >> filenumber;
  OpenFile(filenumber);
  maxline = 0;
//  infile.seekg(0);
  while(!infile.eof()){
	infile.ignore(200,'\n');
	maxline++;
//	if(maxline%100 == 0 ) cerr << "Read line: " << maxline << endl;
  }

  maxline--;
  cerr << "Data file has Total line " << maxline << endl;
  infile.seekg(0);
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
      
       
      //Code below: Generating Fixed Asymmetry distribution in 3D
      //Consider artificially generate e+ e- pair
      // some bug at BH/gun/set_pair_mode
      if (pair_mode)
	{

	  /*
	  G4double A,a, Temp, Randa, ThetaM;
	  //Assigning Asymmetry
	  A = 0.5;
	  a = (1-A)/(1+A);

	  //Compute electron energy, momentum components
	  // Generating Fixed Asymmetry distritution in 3D
	  {
	    do{
	      Randa = CLHEP::RandFlat::shoot(0.0 , 1+a);
	      if (Randa <= 1.0){
		Ee = CLHEP::RandFlat::shoot(15.0*MeV, 30.0*MeV);
	      }
	      else{
		Ee = CLHEP::RandFlat::shoot(30.0*MeV, 45.0*MeV);		  
	      }
	      
		do{
		  Temp = sqrt(((E-Ee)*(E-Ee)-Me*Me)/(Ee*Ee-Me*Me));	  
		  if(-1 < Temp && 1 > Temp){
		    ThetaM = asin(Temp)*rad;
		    Thetae = CLHEP::RandFlat::shoot(0.0, ThetaM);
		  }
		  else
		    Thetae = CLHEP::RandFlat::shoot(0.0*deg,90.0*deg);
	    
		  Phie = CLHEP::RandFlat::shoot(0.0*deg, 360.0*deg);

		  Pe = sqrt(Ee*Ee-Me*Me);
		  Pex = Pe*cos(Thetae);
		  Pey = Pe*sin(Thetae)*cos(Phie);
		  Pez = Pe*sin(Thetae)*sin(Phie);
		  KEe = Ee - Me;
	    
		  Ep = E - Ee;
		  Ppx = sqrt((E-Ee)*(E-Ee)-(Ee*Ee-Me*Me)*sin(Thetae)*sin(Thetae)-Me*Me);
		  Ppy = -Pey;
		  Ppz = -Pez;
		  KEp = Ep - Me;
		  
		  //Note: Thetap here is always in 0 to 90deg
		  Thetap = atan(+Pe*sin(Thetae)/Ppx);
		  //Make positron opposite direction than electron, and lies in 0 to 360deg
		  if (Phie >= 0.*deg && Phie < 180.0*deg)
		    Phip = Phie + 180.0*deg;
		  else if ( Phie >= 180.0*deg && Phie < 360.0*deg)
		    Phip = Phie - 180.0*deg;
	  
		}while(Thetae > 20.0*deg || Thetap > 20.0*deg || Thetae < -20.0*deg || Thetap < -20.0*deg);

	      }while (((E-Ee)*(E-Ee)-Pey*Pey-Me*Me) < 0.0 || Ee < Me);	      
	    //	cerr << "Thetae = " << Thetae/deg << " deg. Thetap = " << Thetap/deg << " deg." << endl;
	    //		cerr << "Ee = " << Ee/MeV << " MeV. Ep = "  << Ep/MeV << " MeV." << endl;
		      // && (Ee < 15.0*MeV) && (Ee > 45.0*MeV) && (abs(Thetae) > 20.0*deg) && (abs(Thetap) > 20.0*deg));
	  }
	  */
	  G4int CurrentEvent = anEvent->GetEventID(); 
	  if(CurrentEvent == 1) { Convert(1) ;}
	  else ConvertNext();
	  //cerr << "Primary data file converted." << endl;
	  //PrintVertex();

	  // Code above: Generate e+ e- pair with same energy and momentum distribution 
	  G4ParticleTable* ParticleTable = G4ParticleTable::GetParticleTable();
	  G4String ParticleName;

	  G4ParticleDefinition * P_electron = ParticleTable->FindParticle(ParticleName="e-");
	  particleGun->SetParticleDefinition(P_electron);
	  
	  // Primary particle are generated uniformly within depth of target.
	  targ_in =  target_position + CLHEP::RandFlat::shoot(-target_thickness/2., target_thickness/2.);
	  //targ_in = target_position;
	  particleGun->SetParticlePosition(G4ThreeVector(targ_in, x_in, y_in));
	  particleGun->SetParticleEnergy(KEe);
	  particleGun->SetParticleMomentumDirection(G4ThreeVector(Pex, Pey, Pez));
	  particleGun->GeneratePrimaryVertex(anEvent);

	  // for output: index 0: positron, index 1:eletron.
	  pOutputFile->Set_energy_i(1,KEe);
	  //pOutputFile->Set_delta_i(1,delta_in);
	  pOutputFile->Set_x_i(1,x_in);
	  pOutputFile->Set_y_i(1,y_in);
	  pOutputFile->Set_theta_i(1,Thetae);
	  pOutputFile->Set_phi_i(1,Phie);

      	  //cerr << "Ee = " << Ee/MeV << " MeV." << "Pex = " << Pex/MeV << " MeV/c."<< "Pey = " << Pey/MeV << " MeV/c."<< "Thetae = "<< Thetae/deg << endl;

	  G4ParticleDefinition * P_positron = ParticleTable->FindParticle(ParticleName="e+");
	  particleGun->SetParticleDefinition(P_positron);
	  particleGun->SetParticlePosition(G4ThreeVector(targ_in, x_in, y_in));
	  particleGun->SetParticleEnergy(KEp);
	  particleGun->SetParticleMomentumDirection(G4ThreeVector(Ppx, Ppy, Ppz));
	  particleGun->GeneratePrimaryVertex(anEvent);

	  pOutputFile->Set_energy_i(0,KEp);
	  //pOutputFile->Set_delta_i(1,delta_in);
	  pOutputFile->Set_x_i(0,x_in);
	  pOutputFile->Set_y_i(0,y_in);
	  pOutputFile->Set_theta_i(0,Thetap);
	  pOutputFile->Set_phi_i(0,Phip);

	}
      //Code Above:Fixed Asymmetry 3D */
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

// Code for sphereical angle
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

// Code for projection angle
int PrimaryGeneratorAction::ConvertNext(){
	if(!infile.eof()){
		infile >> KEe >> KEp >> Thetae >> Phie >> Thetap >> Phip;
		infile.ignore(200, '\n');
		KEe *=MeV;
		KEp *=MeV;
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
