//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "PrimaryGeneratorAction.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "Randomize.hh"
#include "math.h"
#define PI 3.14159265
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

inline double randfloat(){
	double r = random()/ double (RAND_MAX);
	return r;
}

PrimaryGeneratorAction::PrimaryGeneratorAction( DetectorConstruction* DC)
:Detector(DC)
{
  //pair_mode = true;
  pair_mode = false;

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

  E = 60*MeV;
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
      // Consider artificially generate e+ e- pair
      // FIXME:some bug at BH/gun/set_pair_mode

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
	  // Code above: Generate e+ e- pair with same energy and momentum distribution 
	  */



	  G4int CurrentEvent = anEvent->GetEventID(); 
	  //cerr << "Current Event = " << CurrentEvent << endl;
	  ConvertNext();
	  //PrintVertex();

	  G4ParticleTable* ParticleTable = G4ParticleTable::GetParticleTable();
	  G4String ParticleName;
	  // Primary particle are generated uniformly within depth of target.
	  targ_in =  target_position + CLHEP::RandFlat::shoot(-target_thickness/2., target_thickness/2.);
	  //targ_in = target_position;

	/*
	  // Throwing uniformly random solid angle at circumference of cone:
	  // Numbers in DetectorConstruction:
	  double inch = 2.54*cm;
	  double fTargetDistance = 42.26*cm;
	  double fMagnetX = 30.75*inch;
	  double fYokeInnerX = 22.75*inch;
	  double fYokeSideX = (fMagnetX - fYokeInnerX)/2.;

	  double deg = PI/180.;
	  double fConeAddZ = 5.82*cm ;
	  //double fConeAddZ = 0.;
	  double fConeAngle = 5 * deg ;

	  double fConeR = ( fTargetDistance - fYokeInnerX/2. - fConeAddZ ) * tan(fConeAngle) + 0.001*cm;
	  targ_in = -( fYokeInnerX/2. + fConeAddZ );

	  cerr << "fConeR = " << fConeR/cm << "cm, targ_in = " << targ_in/cm << "cm."<< endl;
	  cerr << "continue?" << endl;
	  string ss;
	  cin >> ss;

	  double ang = CLHEP::RandFlat::shoot(0., 2*PI);
	  //cerr << "angle = " << ang/deg << "(deg)" << endl;
	  x_in = fConeR * cos(ang);
	  y_in = fConeR * sin(ang);



		// Below: use uniformly random spherical angle (rather than reading pairs from data file)
		// with theta from 4 to 90 deg as pairticle angles at primary vertex
		double th_min = 4 * deg;
		double th_max = 90 * deg;
		double t_Thetae = acos( cos(th_min) - randfloat()*( cos(th_min) - cos(th_max) ) );
		double t_Thetap = acos( cos(th_min) - randfloat()*( cos(th_min) - cos(th_max) ) );
		//Phie = randfloat() * 360.*deg;
		//Phip = randfloat() * 360.*deg;
		double t_Phie =  90* deg + randfloat() * 180. * deg;
		double t_Phip = -90* deg + randfloat() * 180. * deg;

		// From spherical angle to projection angle as input
		Thetae = atan( tan(t_Thetae) * cos(t_Phie) );
		Phie   = atan( tan(t_Thetae) * sin(t_Phie) );
		Thetap = atan( tan(t_Thetap) * cos(t_Phip) );
		Phip   = atan( tan(t_Thetap) * sin(t_Phip) );

		// [18, 42] MeV/C momentum
		Pe = 18.*MeV + randfloat()*24.*MeV;
		Pp = 18.*MeV + randfloat()*24.*MeV;
		//Pe = 30.*MeV;
		//Pp = 30.*MeV;

		// Makes (x,y,z) right hand
		Pex = Pe/sqrt(1+tan(Thetae)*tan(Thetae)+tan(Phie)*tan(Phie));
		Pey = Pex*tan(Thetae);
		Pez = Pex*tan(Phie);

		Ppx = Pp/sqrt(1+tan(Thetap)*tan(Thetap)+tan(Phip)*tan(Phip));
		Ppy = Ppx*tan(Thetap);
		Ppz = Ppx*tan(Phip);

	KEe = sqrt(Pe*Pe - Me*Me) -Me;
	KEp = sqrt(Pp*Pp - Me*Me) -Me;
	*/





	  // for output: index 0: positron, index 1:eletron.
	  // Positron generation
	  G4ParticleDefinition * P_positron = ParticleTable->FindParticle(ParticleName="e+");
	  particleGun->SetParticleDefinition(P_positron);
	  particleGun->SetParticlePosition(G4ThreeVector(targ_in, x_in, y_in));
	  particleGun->SetParticleEnergy(KEp);
	  // Note Thetae, Phie are projection angles than spherical angles
	  int flag = 1;
	  // FIXME: if SetParticleMomentumDirection() swapped Pp and Pe,
	  // SWAP Thetae, Phie in Set_theta_i(), Set_phi_i()
	  // Otherwise in root file Input Phi_i/Theta_i[0/1] will be correspond to Monitor Phi_m/Theta_m[1/0]
	  // FIXME: Below two lines comment maybe wrong
	  // flag == 1: original magnetic field direction: B vertical down direction, e- goes right, e+ left
	  // flag == 0: reversed magnetic field direction: B vertical up direction,   e- goes left,  e+ right
	  if( flag == 1){
		particleGun->SetParticleMomentumDirection(G4ThreeVector(Ppx, Ppy, Ppz));
		pOutputFile->Set_theta_i(0,Thetap);
	  }
	  else{
		particleGun->SetParticleMomentumDirection(G4ThreeVector(Ppx, -Ppy, Ppz));
		pOutputFile->Set_theta_i(0, -Thetap);
	  }
	  particleGun->GeneratePrimaryVertex(anEvent);

	  pOutputFile->Set_energy_i(0,KEp);
	  //pOutputFile->Set_delta_i(1,delta_in);
	  pOutputFile->Set_x_i(0,x_in);
	  pOutputFile->Set_y_i(0,y_in);
	  pOutputFile->Set_phi_i(0,Phip);

	  // Electron generation
	  G4ParticleDefinition * P_electron = ParticleTable->FindParticle(ParticleName="e-");
	  particleGun->SetParticleDefinition(P_electron);
	  particleGun->SetParticlePosition(G4ThreeVector(targ_in, x_in, y_in));
	  particleGun->SetParticleEnergy(KEe);
	  // Note Thetae, Phie are projection angles than spherical angles
	  if(flag == 1){
		particleGun->SetParticleMomentumDirection(G4ThreeVector(Pex, Pey, Pez));
		pOutputFile->Set_theta_i(1,Thetae);
	  }
	  else{
		particleGun->SetParticleMomentumDirection(G4ThreeVector(Pex, -Pey, Pez));
		pOutputFile->Set_theta_i(1,-Thetae);
	  }
	  particleGun->GeneratePrimaryVertex(anEvent);
	  pOutputFile->Set_energy_i(1,KEe);
	  //pOutputFile->Set_delta_i(1,delta_in);
	  pOutputFile->Set_x_i(1,x_in);
	  pOutputFile->Set_y_i(1,y_in);
	  pOutputFile->Set_phi_i(1,Phie);

	/*
	cerr << endl << endl;
	// Note: deg = PI/180
	cerr << "Input 0: Energy = " << KEp << ", (Ppx, Ppy, Ppz) = " << "(" << Ppx <<  " , " << Ppy << " , " << Ppz << ")" << ", Ppy/Ppz = " << Ppy/Ppz << endl;
	cerr << "t_theta_p = " << atan2(Ppy, Ppx)/deg << ", t_phi_p  = " << atan2(Ppz, Ppx)/deg << endl;
	cerr << "Input 1: Energy = " << KEe << ", (Pex, Pey, Pez) = " << "(" << Pex <<  " , " << Pey << " , " << Pez << ")" << ", Pey/Pez = " << Pey/Pez << endl;
	cerr << "t_theta_e = " << atan2(Pey, Pex)/deg << ", t_phi_e  = " << atan2(Pez, Pex)/deg << endl;
	*/
	}
      //Code Above: Read pair data file from generator obeys theory calculation. */
	else{
	  G4int CurrentEvent = anEvent->GetEventID(); 
	  // Code above: Generate e+ e- pair with same energy and momentum distribution 
	  G4ParticleTable* ParticleTable = G4ParticleTable::GetParticleTable();
	  G4String ParticleName;

	  G4ParticleDefinition * P_gamma = ParticleTable->FindParticle(ParticleName="gamma");
	  particleGun->SetParticleDefinition(P_gamma);
	  // Primary particle are generated uniformly within depth of target.
	  //targ_in =  target_position + CLHEP::RandFlat::shoot(-target_thickness/2., -target_thickness/2.);
	  //targ_in =  target_position;

	  // For checking background, using gamma before target:
	  targ_in = target_position - target_thickness/2.;
	  G4double E_gamma = 60.; 
	  particleGun->SetParticlePosition(G4ThreeVector(targ_in, x_in, y_in));
	  particleGun->SetParticleEnergy(E_gamma);
	  particleGun->SetParticleMomentumDirection(G4ThreeVector(+1., 0., 0.));
	  particleGun->GeneratePrimaryVertex(anEvent);

	  pOutputFile->Set_energy_i(0,E_gamma);
	  pOutputFile->Set_x_i(0,x_in);
	  pOutputFile->Set_y_i(0,y_in);
	  pOutputFile->Set_theta_i(0,0.);
	  pOutputFile->Set_phi_i(0,0.);
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

// Code for projection angle stored in pair data file
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
		bool b_s_dir = 0;
		if(b_s_dir){
			Thetae = 6.*deg;
			Thetap = -6.*deg;
			Phie = 30*deg;
			Phip = -20*deg;
		}

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
