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
      
      /*
      //Code Below: fixed asymmetry 2D
      //Consider artificially generate e+ e- pair with fixed asymmetry in 2D
  
      if (pair_mode)
	{
	  
	  //Code below: generate e+ e- pair with same energy and theta distribution.
	  E = 60*MeV;
	  c = CLHEP::c_light;
	  Me = CLHEP::electron_mass_c2;
	  
	  G4double A,a, Temp, Randa, ThetaM;
	  //Assigning Asymmetry
	  A = 0.2;
	  a = (1-A)/(1+A);

	  //Compute electron energy, momentum components
 
	  //G4double Rand = CLHEP::RandFlat::shoot(0.0, 1.0);
	  //if(Rand <= 0.50){
	  // Generating Fixed Asymmetry distritution in 2D
	  if(1){
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
		    Thetae = CLHEP::RandFlat::shoot(-ThetaM, ThetaM);
		  }
		  else
		    Thetae = CLHEP::RandFlat::shoot(-90.0*deg,90.0*deg);
	    
		  Pe = sqrt(Ee*Ee-Me*Me);
		  Pex = Pe*cos(Thetae);
		  Pey = Pe*sin(Thetae);
		  KEe = Ee - Me;
	    
		  Ep = E - Ee;
		  Ppx = sqrt((E-Ee)*(E-Ee)-(Ee*Ee-Me*Me)*sin(Thetae)*sin(Thetae)-Me*Me);
		  Ppy = -sqrt(Ee*Ee-Me*Me)*sin(Thetae);
		  KEp = Ep - Me;
		  Thetap = atan2(Ppy,Ppx);
		}while(Thetae > 20.0*deg || Thetap > 20.0*deg || Thetae < -20.0*deg || Thetap < -20.0*deg);

	      }while (((E-Ee)*(E-Ee)-Pey*Pey-Me*Me) < 0.0 || Ee < Me);	      
		cerr << "Thetae = " << Thetae/deg << " deg. Thetap = " << Thetap/deg << " deg." << endl;
		cerr << "Ee = " << Ee/MeV << " MeV. Ep = "  << Ep/MeV << " MeV." << endl;
		      // && (Ee < 15.0*MeV) && (Ee > 45.0*MeV) && (abs(Thetae) > 20.0*deg) && (abs(Thetap) > 20.0*deg));
	  }
      */

	  //Fixed energy , x-only direction beam
	  /*
	  Ee = 45.0*MeV;
	  Pe = sqrt(Ee*Ee-Me*Me);
	  Pex = Pe;
	  Pey = 0;
	  KEe = Ee - Me;
	  
	  Ep = E - Ee;
	  Ppx = sqrt(Ep*Ep- Me*Me);
	  Ppy = 0;
	  KEp = Ep - Me;
	  */
	  
	  /*
	  else{
	      do{
		Randa = CLHEP::RandFlat::shoot(0.0 , 1+a);
		if (Randa <= 1.0){
		  Ep = CLHEP::RandFlat::shoot(0.0*MeV, 30.0*MeV);
		}
		else{
		  Ep = CLHEP::RandFlat::shoot(30.0*MeV, 60.0*MeV);		  
		}


		//Assigning angle theta
		Temp = sqrt(((E-Ep)*(E-Ep)-Me*Me)/(Ep*Ep-Me*Me));	  
		if(-1 < Temp && 1 > Temp){
		  ThetaM = asin(Temp);
		  Thetap = CLHEP::RandFlat::shoot(-ThetaM*rad, ThetaM*rad);
		}
		else
		  Thetap = CLHEP::RandFlat::shoot(-90*deg,90*deg);
		
		Pp = sqrt(Ep*Ep-Me*Me);
		Ppx = Pp*cos(Thetap);
		Ppy = Pp*sin(Thetap);
		KEp = Ep - Me;

		//Compute positron energy, momentum components
		Ee = E - Ep;
		Pex = sqrt((E-Ep)*(E-Ep)-Ppy*Ppy-Me*Me);
		Pey = -Ppy;
		Thetae = atan(Pey/Pex);
		KEe = Ee - Me;

	      }while (Ep < Me || ((E-Ep)*(E-Ep)-Ppy*Ppy-Me*Me) < 0.0);
	      }*/
	  /////////////////////////
	  // Code above: Generate e+ e- pair with same energy and momentum distribution
      /*	  
	  G4ParticleTable* ParticleTable = G4ParticleTable::GetParticleTable();
	  G4String ParticleName;

	  G4ParticleDefinition * P_electron = ParticleTable->FindParticle(ParticleName="e-");
	  particleGun->SetParticleDefinition(P_electron);
	  
	  targ_in = target_position;
	  //targ_in =  target_position + CLHEP::RandFlat::shoot(-target_thickness/2., target_thickness/2.);
	  particleGun->SetParticlePosition(G4ThreeVector(targ_in, x_in, y_in));
	  particleGun->SetParticleEnergy(KEe);
	  particleGun->SetParticleMomentumDirection(G4ThreeVector(Pex, Pey, 0.0));
	  particleGun->GeneratePrimaryVertex(anEvent);

	  pOutputFile->Set_energy_i(1,KEe);
	  //pOutputFile->Set_delta_i(1,delta_in);
	  pOutputFile->Set_x_i(1,x_in);
	  pOutputFile->Set_y_i(1,y_in);
	  pOutputFile->Set_theta_i(1,abs(Thetae));
	  if(Pey > 0.)pOutputFile->Set_phi_i(1, 0.0);
	  else pOutputFile->Set_phi_i(1, 180.0*deg);

	  //	  cerr << "Ee = " << Ee/MeV << " MeV." << "Pex = " << Pex/MeV << " MeV/c."<< "Pey = " << Pey/MeV << " MeV/c."<< "Thetae = "<< Thetae/deg << endl;

	  G4ParticleDefinition * P_positron = ParticleTable->FindParticle(ParticleName="e+");
	  particleGun->SetParticleDefinition(P_positron);
	  particleGun->SetParticlePosition(G4ThreeVector(targ_in, x_in, y_in));
	  particleGun->SetParticleEnergy(KEp);
	  particleGun->SetParticleMomentumDirection(G4ThreeVector(Ppx, Ppy, 0.0));
	  particleGun->GeneratePrimaryVertex(anEvent);

	  pOutputFile->Set_energy_i(0,KEp);
	  //pOutputFile->Set_delta_i(1,delta_in);
	  pOutputFile->Set_x_i(0,x_in);
	  pOutputFile->Set_y_i(0,y_in);
	  pOutputFile->Set_theta_i(0,abs(Thetap));
	  if( Ppy > 0.) pOutputFile->Set_phi_i(0, 0.0);
	  else pOutputFile->Set_phi_i(0, 180.0*deg);
	}

      */
      //Code Above: fixed asymmetry 2D



       
      //Code below: Generating Fixed Asymmetry distribution in 3D
      //Consider artificially generate e+ e- pair
      // some bug at BH/gun/set_pair_mode
      if (pair_mode)
	{
	  
	  //Code below: generate e+ e- pair with same energy and theta distribution.
	  E = 60*MeV;
	  c = CLHEP::c_light;
	  Me = CLHEP::electron_mass_c2;

	  G4double A,a, Temp, Randa, ThetaM;
	  //Assigning Asymmetry
	  A = 0.0;
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
	  G4ParticleTable* ParticleTable = G4ParticleTable::GetParticleTable();
	  G4String ParticleName;

	  G4ParticleDefinition * P_electron = ParticleTable->FindParticle(ParticleName="e-");
	  particleGun->SetParticleDefinition(P_electron);
	  
	  //targ_in =  target_position + CLHEP::RandFlat::shoot(-target_thickness/2., target_thickness/2.);
	  targ_in = target_position;
	  particleGun->SetParticlePosition(G4ThreeVector(targ_in, x_in, y_in));
	  particleGun->SetParticleEnergy(KEe);
	  particleGun->SetParticleMomentumDirection(G4ThreeVector(Pex, Pey, Pez));
	  particleGun->GeneratePrimaryVertex(anEvent);

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


      /*
      //Symmetric energy distribution 3D
      // Code below: Generte e+ e- pair with same energy and momentum distribution in 3D
      if (pair_mode)
	{	  
	  E = 60*MeV;
	  c = CLHEP::c_light;
	  Me = CLHEP::electron_mass_c2;

	  //Compute electron energy, momentum components
 
	  G4double Rand = CLHEP::RandFlat::shoot(0.0, 1.0);
	  if(Rand <= 0.50){
	      do{
		Ee = CLHEP::RandFlat::shoot(20.00*MeV, 40.0*MeV);
		//Ee = CLHEP::RandFlat::shoot(0.0*MeV, 60.0*MeV);
		Thetae = CLHEP::RandFlat::shoot( 0.0, 1);
		Phie = CLHEP::RandFlat::shoot(0.0*deg, 360.0*deg);
		Thetae = (acos(-2*Thetae + 1))*rad;

		Pe = sqrt(Ee*Ee-Me*Me);
		Pex = Pe*cos(Thetae);
		Pey = Pe*sin(Thetae)*cos(Phie);
		Pez = Pe*sin(Thetae)*sin(Phie);
		KEe = Ee - Me;

		//Compute positron energy, momentum components
		Ep = E - Ee;
		Ppx = sqrt((E-Ee)*(E-Ee)-Pey*Pey-Pez*Pez-Me*Me);
		Ppy = -Pey;
		Ppz = -Pez;
		KEp = Ep - Me;

	      }while (Ee < Me || ((E-Ee)*(E-Ee)-Pey*Pey-Pez*Pez-Me*Me) < 0.0 || Thetae > 20.0*deg);

	      Thetap = (atan2(sqrt(Ppy*Ppy+Ppz*Ppz),Ppx))*rad;
	  }
	  else {
	      do{
		Ep = CLHEP::RandFlat::shoot(20.0*MeV, 40.0*MeV);
		//Ep = CLHEP::RandFlat::shoot(0.0*MeV, 60.0*MeV);
		Thetap = CLHEP::RandFlat::shoot( 0.0, 1);
		Phip = CLHEP::RandFlat::shoot(0.0*deg, 360.0*deg);
		Thetap = (acos(-2*Thetap + 1))*rad;

		Pp = sqrt(Ep*Ep-Me*Me);
		Ppx = Pp*cos(Thetap);
		Ppy = Pp*sin(Thetap)*cos(Phip);
		Ppz = Pp*sin(Thetap)*sin(Phip);
		KEp = Ep - Me;

		//Compute positron energy, momentum components
		Ee = E - Ep;
		Pex = sqrt((E-Ep)*(E-Ep)-Ppy*Ppy-Ppz*Ppz-Me*Me);
		Pey = -Ppy;
		Pez = -Ppz;
		KEe = Ee - Me;

	      }while (Ep < Me || ((E-Ep)*(E-Ep)-Ppy*Ppy-Ppz*Ppz-Me*Me) < 0.0 || Thetap > 20.0*deg);

	      Thetae = (atan2(sqrt(Pey*Pey+Pez*Pez),Pex))*rad;
	  }




	  //Below just assigning data to particle gun.
	  G4ParticleTable* ParticleTable = G4ParticleTable::GetParticleTable();
	  G4String ParticleName;

	  G4ParticleDefinition * P_electron = ParticleTable->FindParticle(ParticleName="e-");
	  particleGun->SetParticleDefinition(P_electron);
	  
	  //targ_in =  target_position + CLHEP::RandFlat::shoot(-target_thickness/2., target_thickness/2.);
	  targ_in = target_position;
	  particleGun->SetParticlePosition(G4ThreeVector(targ_in, x_in, y_in));
	  particleGun->SetParticleEnergy(KEe);
	  particleGun->SetParticleMomentumDirection(G4ThreeVector(Pex, Pey, Pez));
	  particleGun->GeneratePrimaryVertex(anEvent);

	  pOutputFile->Set_energy_i(1,KEe);
	  //pOutputFile->Set_delta_i(1,delta_in);
	  pOutputFile->Set_x_i(1,x_in);
	  pOutputFile->Set_y_i(1,y_in);
	  pOutputFile->Set_theta_i(1,Thetae);
	  pOutputFile->Set_phi_i(1,Phie);

	  //	  cerr << "Ee = " << Ee/MeV << " MeV." << "Pex = " << Pex/MeV << " MeV/c."<< "Pey = " << Pey/MeV << " MeV/c."<< "Thetae = "<< Thetae/deg << endl;

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
	  //	  cerr << "Ep = " << Ep/MeV << " MeV/c." << "Ppx = " << Ppx/MeV << " MeV/c."<< "Ppy = " << Ppy/MeV << " MeV/c."<< "Thetap = "<< Thetap/deg << endl;
	  //	  cerr << "Ee + Ep = " << (Ee+Ep)/MeV << "MeV." << endl;	
      //Code Above:Symmetric energy distribution 3D
      */



    /*
      //Code Below: Normal case or gamma mode.

      else
	{
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

      pOutputFile->Set_energy_i(0,energy_in);
      pOutputFile->Set_delta_i(0,delta_in);
      pOutputFile->Set_x_i(0,x_in);
      pOutputFile->Set_y_i(0,y_in);
      pOutputFile->Set_theta_i(0,theta_in);
      pOutputFile->Set_phi_i(0,phi_in);
	
      //G4cout << "Generate Primary...." << G4endl;
      particleGun->GeneratePrimaryVertex(anEvent);
      //G4cout << "Generated Primary...." << G4endl;
      }
      */
      //Code Above: Normal case or gamma mode.
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
