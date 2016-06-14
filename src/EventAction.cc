#include "EventAction.hh"
#include "RunAction.hh"
#include "G4Event.hh"
#include "G4UnitsTable.hh"
#include "Randomize.hh"
#include <iomanip>
#include "OutputFile.hh"

class OutputWire;
class SingleWireHit;
class WireHit;
#define PI 3.14159265

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

EventAction::EventAction(DetectorConstruction* det)
:printModulo(100),pDetector(det)

{
  //eventMessenger = new EventActionMessenger(this);
  pOutputFile = OutputFile::GetOutputFilePointer();
  fVDCSigma = 0.1*mm; //std. dev. of position resolution.
  fRequireHodoscopeHit = true;
  fPaddleThreshold = 20.*keV;
  fWireThreshold = 0.2*keV;
	// these may need to be adjusted
  fsigmaLightStatistical = 77.0*keV;
  fsigmaLightNoise = 5.0*keV;
  //pOutputFile->ClearKE();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

EventAction::~EventAction()
{
  //delete eventMessenger;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void EventAction::BeginOfEventAction(const G4Event* evt)
{  
  //G4cout << "Begin of Event Action called." << G4endl;
  G4int event_number = evt->GetEventID();
  // get needed detector info here because it may have changed since constructor
  fVDCSpacing = pDetector->GetVDCSpacing();
  pOutputFile->ClearKE();
  if (event_number < 5 || event_number%printModulo == 0)
	{ 
	if(event_number >= 10*printModulo && printModulo < 100000) printModulo *= 10;
	G4cout << "---> Begin of event: " << event_number << G4endl;
    //CLHEP::HepRandom::showEngineStatus();
	}
 
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.... 

void EventAction::EndOfEventAction(const G4Event* evt)
{

  //G4cout << "End of event action called." << G4endl; //debug
  G4int event_number = evt->GetEventID();
  // --------------- Monitor -----------------
  // Skip this part if we are not going to output Monitor info
  G4bool use_monitor = pOutputFile->GetUseMonitor();

if(use_monitor)
  {
  G4int MONID = G4SDManager::GetSDMpointer()->GetCollectionID("Monitor/MonitorHitCollection");
  MonitorHitCollection * MonitorHC = (MonitorHitCollection *) (evt->GetHCofThisEvent()->GetHC(MONID));
  G4int totalMonHits = MonitorHC->entries();
  MonitorHit * mHit;
  // We just need to find earliest time for each of e- and e+ hits
  G4double min_time[2]; G4int min_time_hit[2]; G4int ipart;
  G4bool mon_hit[2];
  for(G4int i=0 ; i < 2; i++)
	{
  	min_time[i] = 0;
	min_time_hit[i] = -1;
	mon_hit[i] = false;
	}
  for(G4int jj = 0; jj < totalMonHits; jj++)
	{
	mHit = (*MonitorHC)[jj];
	G4double now = mHit->GetTime();
	if(mHit->GetParticle() == "e-") ipart = 1;
	else if(mHit->GetParticle() == "e+") ipart = 0;
	else continue;
	if( !mon_hit[ipart] || min_time[ipart] > now)
		{
		min_time[ipart] = now;
		min_time_hit[ipart] = jj;
		mon_hit[ipart] = true;
		}
	}

  //we need both monitor hit to get result of e+ /e- pair production
  if(mon_hit[0] && mon_hit[1])
    //  if(mon_hit[0] || mon_hit[1])
	{
  	for(ipart = 0; ipart < 2; ipart++)
		{
		pOutputFile->SetMonitorHit(ipart, mon_hit[ipart] );
		if(mon_hit[ipart])
			{
			mHit = (*MonitorHC)[min_time_hit[ipart]];
			//cerr << "Time when monitor hit: " << min_time[ipart]/ns << "(ns)" << endl;
			// Set into output file
			// Translating to RAYTRACE coordinates
			pOutputFile->SetMonitorKineticEnergy(ipart, mHit->GetKineticEnergy() );
			G4ThreeVector pos = mHit->GetGlobalPrePosition();
			pOutputFile->SetMonitorX(ipart, pos.y() );
			pOutputFile->SetMonitorY(ipart, pos.z() );
			G4ThreeVector dir = mHit->GetMomentumDirection();
			//G4double theta_m = atan2( dir.y(), dir.x() );
			//G4double phi_m = atan2( dir.z(), dir.x() );

			// theta_m , phi_m here are spherical angle about x-axis.
			//Since direction from target to monitorSD always makes dir.x positive, theta_m is always in (0,PI/2)
			G4double theta_m = atan2( sqrt(dir.y()*dir.y()+dir.z()*dir.z()), dir.x() );
			G4double phi_m = 0.0;
			if( (dir.z() >= 0. && dir.y() > 0.) )
			  phi_m = atan( dir.z()/dir.y() );
			else if( (dir.z() <= 0. && dir.y() < 0.) || (dir.z() >= 0. && dir.y() < 0.) )
			  phi_m = atan( dir.z()/dir.y() ) + PI;
			else if( dir.z() <= 0. && dir.y() > 0.)
			  phi_m = atan( dir.z()/dir.y() ) + 2.*PI;
			else if ( dir.z() > 0. && dir.y() == 0.)
			  phi_m = PI/2.;
			else if ( dir.z() < 0. && dir.y() == 0.)
			  phi_m = PI/2.;

			pOutputFile->SetMonitorTheta(ipart, theta_m );
			pOutputFile->SetMonitorPhi(ipart, phi_m );
			}
		}
	}
  else
	{ use_monitor = false; } // surpress writing event if no monitor hits
  }
	
			

  // --------------- VDC ---------------------
  G4int CHCID = G4SDManager::GetSDMpointer()->GetCollectionID("WireChamber/HitCollection");
  WireChamberHitsCollection * chamberHC = (WireChamberHitsCollection *) (evt->GetHCofThisEvent()->GetHC(CHCID));
  
  
  G4int totalHits = chamberHC->entries();
  //G4cout << " Event " << event_number << ": Number of chamber hits = " << totalHits  << G4endl;//debug
  
  WireChamberHit * aHit;
  // Temporary hit data
  vector<G4ThreeVector> position[2][2][2]; // position on [package][chamber][layer] and wire
  vector<G4double> energyTotal[2][2][2];
  //(unique) hit wire number in this event
  vector<G4int> WireNum[2][2][2];
  vector<G4double> kineticenergy[2][2][2];
  vector<G4double> charge[2][2][2];
  //vector<G4double> particle[2][2][2];
  vector<G4double> VDC_time[2][2][2];
  //vector<G4double> kineticenergy[2]; 
  //vector<G4double> KEavg[2];
  //G4int KENum[2];
  vector<WireHit> VDCWireHit[2][2][2];
  
  // Number of hits for each hit wire in this event:
  vector<G4int> WireHitCount[2][2][2];
  // Temporary: position of each step
  G4ThreeVector pos;  
  // Temporary: Number of steps in a wire hit
  G4int Num[2][2][2];

  // Initialize all data needed to be stored
  for(G4int i = 0; i < 2; i++)
    for(G4int j = 0; j < 2; j++)
      for(G4int k = 0; k < 2; k++) {
	position[i][j][k].clear();
	energyTotal[i][j][k].clear();
	WireNum[i][j][k].clear();
	kineticenergy[i][j][k].clear();
	charge[i][j][k].clear();
	//particle[i][j][k].clear();
	VDC_time[i][j][k].clear();
	WireHitCount[i][j][k].clear();
	//VDCWireHit[i][j][k].clear();
      }
  
  //Find out how many individual wires hit on each wire plane in this event:
  vector<G4int> wirehit[2][2][2];
  for(G4int jj = 0; jj < totalHits; jj++){
    aHit = (*chamberHC)[jj];
    G4int wire = aHit->GetVDCwire();
    G4int chamber = aHit->GetVDCnumber();
    G4int layer = aHit->GetVDClayer();
    G4int package = aHit->GetVDCpackage();

    //Attempt to generate wirehit object with unique hit wire number
    /*
    if(wirehit[package][chamber][layer].size() == 0) wirehit[package][chamber][layer].push_back(wire);
    else{
    
      for(unsigned int i = 0; i < wirehit[package][chamber][layer].size(); i++){
	G4int Temp = wirehit[package][chamber][layer].at(i);
	if( wire == Temp ) continue;
	else { if(i == (wirehit[package][chamber][layer].size()-1)) wirehit[package][chamber][layer].push_back(wire);
	  else continue;}
	  }}*/
    wirehit[package][chamber][layer].push_back(wire);
  }
  //Till now for a certain [package][chamber][layer] each vector element of wirehit will look like:
  // 100, 50, 51, 51, 51, 50




  // Sort wirehit[i][j][k] in decrement order:
  for(G4int i = 0; i < 2; i++)
      for(G4int j = 0; j < 2; j++)
	for(G4int k = 0; k < 2; k++){
	  if(wirehit[i][j][k].size() >= 2){
	  for(unsigned int l = 1; l < wirehit[i][j][k].size(); l++){
	    for(unsigned int m = 0; m < l; m++){ 
	      G4int TEMP =  wirehit[i][j][k].at(l); 
	      if( wirehit[i][j][k].at(m) < wirehit[i][j][k].at(l)){ 
		wirehit[i][j][k].at(l) =  wirehit[i][j][k].at(m);
		wirehit[i][j][k].at(m) =  TEMP;}
	    }}}
	}

  /*// Print out number of hits on each wire plane:
  for(G4int i = 0; i < 2; i++)
    for(G4int j = 0; j < 2; j++)
      for(G4int k = 0; k < 2; k++){
	cerr << "Detector[" << i << "] Chamber[" << j << "] Layer[" << k << "] has " << wirehit[i][j][k].size() << " hits" << endl;
	for (unsigned int l = 0; l < wirehit[i][j][k].size(); l++)
	  cerr << "    Wire number " << wirehit[i][j][k].at(l) << " hit." << endl;
	  }*/
  //Till now the reordered wirehit for the [package][chamber][layer] will look like:
  //100, 51, 51, 51, 50, 50






  //Find unique hit wire numbers, and associate number of hits in the wire
  for(G4int i = 0; i < 2; i++)
    for(G4int j = 0; j < 2; j++)
      for(G4int k = 0; k < 2; k++){
	//Find number of same hit wire number, store this number as number of hits in the wire, 
	//and store wire number in WireNum object
	while(wirehit[i][j][k].size() != 0){
	  // Find same hit wire number, count it and associate number of hits in the wire.
	  unsigned int Count = 0;
	  // Store the unique distinct hit wire number.
	  WireNum[i][j][k].push_back(wirehit[i][j][k].at(0));
	  // Find how many hits in this wire, store it in WireHitCount
	  for(unsigned int l = 0; l < wirehit[i][j][k].size(); l++){
	    G4int TEMP = wirehit[i][j][k].at(0);
	      if(wirehit[i][j][k].at(l) == TEMP) Count++;
	      else break;}
	  WireHitCount[i][j][k].push_back(Count);

	  // Remove these hit wire numbers in wirehit object, till it's empty
	  for(unsigned int l = Count; l < wirehit[i][j][k].size(); l++){
	    wirehit[i][j][k].at(l-Count) =  wirehit[i][j][k].at(l);}
	  for(unsigned int l = 0; l < Count; l++){ wirehit[i][j][k].pop_back(); }}
      }
  /*
   //Print WireNum and WireHitCount:
  for(G4int i = 0; i < 2; i++)
    for(G4int j = 0; j < 2; j++)
      for(G4int k = 0; k < 2; k++){
	cerr << "VDClayer: ["<<i<<"]["<<j<<"]["<<k<<"]" << endl;
	for(unsigned int l = 0; l < WireNum[i][j][k].size(); l++)
	  {cerr <<
	      " Wire [" << WireNum[i][j][k].at(l) <<"] hit, and number of hits is: " 
		<< WireHitCount[i][j][k].at(l) <<"." << endl;
	  }}
  */
  //Till now, for the [package][chamber][layer],
  //WireNum[package][chamber][layer].at(0), at(1), at(2) will be: 100, 51, 50
  //Corresponding WireHitCount[package][chamber][layer].at(0), at(1), at(2) will be: 1,3,2







  //Set WireHit object has same number of entries as WireNum for SWH
  //cerr << "VDCWireHit object:" << endl;  
  for(G4int i = 0; i < 2; i++)
    for(G4int j = 0; j < 2; j++)
      for(G4int k = 0; k < 2; k++){
	VDCWireHit[i][j][k].clear();
	for(unsigned int l = 0; l < WireNum[i][j][k].size(); l++){
	  WireHit NullWire;
	  NullWire.Set_WireNum(-1);
	  NullWire.Set_Edep(0.);
	  NullWire.Set_Position(G4ThreeVector(0.,0.,0.));
	  VDCWireHit[i][j][k].push_back(NullWire);
	}
	//cerr << "[" << i << "][" << j << "][" << k << "] have "<< VDCWireHit[i][j][k].size() << " different hit wire" << endl;
      }
  // Till now VDCWireHit[package][chamber][layer].at(0),(1),(2) will be NullWire (a WireHit object, store hit data for a single wire):
  // 3 elements for the vector here since 3 individual wire hit(Num 100, 51, 50 )
  // NullWire.WireNum = -1 (since the WireNum indexing is from 0 to MaxWireNum-1, so -1 means WireNum not valid)
  // NullWire.Edep = 0.; (will store total Edep in the wire)
  // NullWire.Position = (0.,0.,0); (will store avg. hit position of all hitcollection steps in the wire)
  // NullWire.SWH (vector<SingleWireHit> object), will be empty (will store all individual hitcollection steps in a certain wire)
  // In a SingleWireHit will have following data: KE, ToF, Particle(PDGEncoding) , Charge




  // Iterate of each hit in hit collection and then hit wires, from 0 to WireHitCount - 1 of the hit wire. 
  // Store following data in each non-empty WireHit type object in VDCWireHit[package][chamber][layer]:
  for(G4int jj = 0; jj < totalHits; jj++) //for all hits
	{
	  aHit = (*chamberHC)[jj];
	  //	aHit->Print();
	  // Get all indexing (copy number) of wire, layer, chamber, and package
	  G4int wire = aHit->GetVDCwire();
	  G4int chamber = aHit->GetVDCnumber();
	  G4int layer = aHit->GetVDClayer();
	  G4int package = aHit->GetVDCpackage();
	  // edep from each step
	  G4double edep = aHit->GetEdep();
	  // Time of Flight, KineticEnergy, Charge, Particle Type and Position of the step
	  G4double time = aHit->GetTime();
	  G4double KE = aHit->GetKE();
	  G4double Charge = aHit->GetCharge();
	  G4int Particle = aHit->GetParticle();
	  G4ThreeVector pos = ( aHit->GetLocalPrePosition() + aHit->GetLocalPostPosition() )/2.;
	  
	  // Compare hit wire number, pacakge, chamber, layer, and
	  // store associate  WireNumber, Position, Edep, Kineticenergy, Charge, ToF in VDCWireHit[package][chamber][layer]
	  for (unsigned int i = 0; i < WireNum[package][chamber][layer].size(); i++){
	    // Compare WireNum from this hit collection step to data in WireNum[package][chamber][layer].at()
	    // Till now at(0),(1),(2) are 100, 51, 50 respectively
	    if(wire ==  WireNum[package][chamber][layer].at(i)){
	      // Store WireNum the same index as WireNum[package][chamber][layer]
	      VDCWireHit[package][chamber][layer].at(i).WireNum = wire;
	      // Add edep from the step to total Edep for the wire
	      VDCWireHit[package][chamber][layer].at(i).Edep += edep;
	      // Add avg portion of Position of the step for the wire
	      VDCWireHit[package][chamber][layer].at(i).Position += pos/(WireHitCount[package][chamber][layer].at(i));

	      //cerr << "Layer [" << package << "][" << chamber << "][" << layer << "]"; cerr<< "Wire " << wire << " hit. ";
	      //cerr<< "Edep = "  << VDCWireHit[package][chamber][layer].at(i).Edep << endl;
	      
	      {
		// The KE, ToF, Charge, Particle info will be stored as a SingleWireHit object as a vector element of the wire
		// Thus when MAX KE is found, coresponding ToF, Charge, Particle of the hitcollection step can be found easily.
		SingleWireHit SingleHit;
		SingleHit.Set_KE(KE);
		SingleHit.Set_ToF(time);
		SingleHit.Set_Charge(Charge);
		SingleHit.Set_Particle(Particle);
		VDCWireHit[package][chamber][layer].at(i).SWH.push_back(SingleHit);}
	      break;
	    }} 
	}

  /*
  //Print WireHit objects:
	  for(G4int i = 0; i < 2; i++)
	    for(G4int j = 0; j < 2; j++)
	      for(G4int k = 0; k < 2; k++){
		for(unsigned int l = 0; l < VDCWireHit[i][j][k].size(); l++){
		  cerr << "Wire[" << i << "][" << j << "][" << k <<"] number [" << VDCWireHit[i][j][k].at(l).WireNum <<"]:"<<endl;
		  cerr << "Total Edep:" << VDCWireHit[i][j][k].at(l).Edep << "." << endl;
		  for(unsigned int m = 0; m < VDCWireHit[i][j][k].at(l).SWH.size(); m++){
		    cerr << "The " << m <<"th hit has KE: " << VDCWireHit[i][j][k].at(l).SWH.at(m).Get_KE()
			 << ". ToF: " << VDCWireHit[i][j][k].at(l).SWH.at(m).Get_ToF()
			 << ". Charge: " << VDCWireHit[i][j][k].at(l).SWH.at(m).Get_Charge()
			 << "." << endl;}}
	      }
  */
	  //Till Now the VDCWireHit[package][chamber][layer] may look like:
	  // VDCWireHit[package][chamber][layer] (vector of WireHit objects) have:
	  
	  // at(0): WireNum = 100, Edep = 10keV, Position = (1,2,3), 
	  // and SWH.at(0) have KE = 40 (MeV), with cooresponding ToF, Charge, Particle.
	   
	  // at(1): WireNum = 51, Edep = 4+3.5+2 = 9.5keV, Position = (()+()+())/3 .
	  // and SWH.at(0) KE = 0.6 MeV, with corresponding ToF, Charge, Particle
	  //     SWH.at(1) KE = 0.8 MeV, with....
	  //     SWH.at(2) KE = 0.75 MeV, with...
	  
	  // at(2): WireNum = 50, Edep = 3+4 = 7 keV, Position = (()+())/2 .
	  // and SWH.at(0) KE = 0.5 MeV, with coresponding ToF, Charge, Particle
	  //     SWH.at(1) KE = 0.6 MeV, with...



	  //If Edep of the wire is bigger than threshold:
	  //Find out MAX KE hit for the wire and associate X, Y, ToF, Edep, Charge, particle type
	  //Now check if total energy deposite > threshold for the wire:
	  //if is, find out largest KE of the hit, and move data to pOutputFile
	  G4bool det_hit[2] = {false, false};
	  // clear fVDC_f[i][j][k] objects before storing data:
	  pOutputFile->Clear_VDC_f();
	  for(G4int i = 0; i < 2; i++)
	    for(G4int j = 0; j < 2; j++)
	      for(G4int k = 0; k < 2; k++){
		for(unsigned int l = 0; l < VDCWireHit[i][j][k].size(); l++){
		  
		  
		  // Only when the wire edep exceed threshold (Now 0.2keV), data will be stored
		  if( VDCWireHit[i][j][k].at(l).Edep > fWireThreshold){
		    det_hit[i] = true;

		    
		    //Find out highest KE at Wire VDCWireHit[i][j][k].at(l), and its index
		    G4double TEMPKE = 0.;
		    G4int TEMP = -1, TEMPPart = 0;
		    for(unsigned int m = 0; m < VDCWireHit[i][j][k].at(l).SWH.size(); m++){
		      if (VDCWireHit[i][j][k].at(l).SWH.at(m).Get_KE() > TEMPKE){
			TEMPKE = VDCWireHit[i][j][k].at(l).SWH.at(m).Get_KE();
			TEMP = m;
		      }}

		    //Put highest KE and related data into fVDC_f in output.
		    //Found out in DetectorConstruction.cc: 
		    //U WirePlane (layer copy number 0) is back WirePlane, V (layer copy number 1) is front WirePlane  
		      OutputWire TEMPWire;
		      // Here k denotes layer copy number
		      if( j== 0 && k==0)TEMPWire.Set_WirePlane_f(1);
		      else if( j==0 && k==1)TEMPWire.Set_WirePlane_f(0);
		      else if( j==1 && k==0)TEMPWire.Set_WirePlane_f(3);
		      else if( j==1 && k==1)TEMPWire.Set_WirePlane_f(2);

		      TEMPWire.Set_WireNum_f(VDCWireHit[i][j][k].at(l).WireNum);
		      TEMPWire.Set_X_f(VDCWireHit[i][j][k].at(l).Position.getX());

		      // In DetectorConstruction cc:Z rotated 180 deg for V WirePlane(layer copy number 1), then need - sign:
		      if(k == 1) TEMPWire.Set_Y_f(-VDCWireHit[i][j][k].at(l).Position.getZ());
		      else TEMPWire.Set_Y_f(VDCWireHit[i][j][k].at(l).Position.getZ());
		      TEMPWire.Set_Edep_f(VDCWireHit[i][j][k].at(l).Edep);
		      TEMPWire.Set_KE_f(TEMPKE);
		      TEMPWire.Set_ToF_f(VDCWireHit[i][j][k].at(l).SWH.at(TEMP).Get_ToF());
		      TEMPWire.Set_Charge_f(VDCWireHit[i][j][k].at(l).SWH.at(TEMP).Get_Charge());
		      TEMPWire.Set_Particle_f(VDCWireHit[i][j][k].at(l).SWH.at(TEMP).Get_Particle());

		      // A particle check condition to store data if needed.
		      TEMPPart = VDCWireHit[i][j][k].at(l).SWH.at(TEMP).Get_Particle();
		      // Store data of each wire in fVDC_f[i][j][k] object when e+ at left, or e- at right
		      //if((i == 0 && TEMPPart == -11) || (i == 1 && TEMPPart == 11))
		      pOutputFile->Set_VDC_f(i,j,k,TEMPWire);
		  }}
		//Now put flag det_hit[i] into and fdetector_package[i], as a flag to write VDC hit data to root file in OutputFile.cc
		pOutputFile->Set_detector_package(i, det_hit[i]);		
	      }
	  // Till now the final data before they're put into root file in OutputFile.cc
	  // The data vector<OutputWire> fVDC_f[package][chamber][layer]may look like:
	  // at(0): WireNum = 100, Edep = 10keV, (X,Y) = (1,3),
	  // KE = 40 (MeV), with cooresponding ToF, Charge, Particle.
	   
	  // at(1): WireNum = 51, Edep = 4+3.5+2 = 9.5keV, (X,Y) = (()+()+())/3 .
	  // KE = 0.8 MeV, with corresponding ToF, Charge, Particle
	  // 2 lower KE steps data are discarded
	  
	  // at(2): WireNum = 50, Edep = 3+ 4 = 7 keV, (X,Y) = (()+())/2 .
	  // KE = 0.6 MeV, with coresponding ToF, Charge, Particle
	  // 1 lower KE step data is discarded

	  // If there is a wire with total Edep less than 0.2 keV, it will not be stored in fVDC_f.




	  // Below: old code without wires in wireplane
	    /*
	  if(edep > 0.)
	  {
	    // Storing edep for each wireplane
	    energyTotal[package][chamber][layer] += edep;
	    pos = ( aHit->GetLocalPrePosition() + aHit->GetLocalPostPosition() )/2.;
	    // Storing number of hit for each wireplane
	    Num[package][chamber][layer] ++;
	    // Edep weighted sum of position
	    //position[package][chamber][layer] += pos*edep;
	    // Sum of position
	    position[package][chamber][layer] += pos;
	    
	    //////////////////////////*
	    // Storing KE at each step in hit collection of 1st wireplane in each side
	    if(layer == 0 && chamber== 0 && aHit->GetKE() > 0.){
	      kineticenergy[package].push_back(aHit->GetKE());
	      KENum[package] ++;} 
	    ////////////////////////////////////////////////
	  }
	  else kineticenergy[package] = 0.;
	}

  
	//Sort kineticenergy in decrement order in each hit wire. While find out the associate ToF, Charge, Particle.
	for (G4int m = 0; m < 2; m++)
	  for (unsigned n = 0; n < kineticenergy[m].size(); n++)
	    for (unsigned l = n; l < kineticenergy[m].size() ; l++ ){
	      if(kineticenergy[m].at(l) > kineticenergy[m].at(n) )
		{ G4double Temp; Temp = kineticenergy[m].at(l); kineticenergy[m].at(l) = kineticenergy[m].at(n); kineticenergy[m].at(n) = Temp;}
	    }
*/


	/*
	// Storing data in OutPutFile object
  G4bool det_hit[2];
  for(G4int i = 0; i < 2; i++)
	{
	det_hit[i] = true;
  	
	for(G4int j = 0; j < 2; j++)
  	for(G4int k = 0; k < 2; k++)
	  {
	    //change following line to deal with vacuum setting.
	    if(energyTotal[i][j][k] < fWireThreshold)
	    //if (kineticenergy[i]/KENum[i] < 2.0*MeV)
	    //if (position[i][j][k].getX() == 0. && position[i][j][k].getY() == 0. && position[i][j][k].getZ() == 0.)
	      { det_hit[i] = false;}
	    else
	      { 
		// Storing largest group KEavg
		//if( KEavg[i].size() > 0 && KEavg[i].at(0) >= 2.*MeV)  pOutputFile->Set_KE_f(i,KEavg[i].at(0));
		    
		    if(energyTotal[i][j][k] > fWireThreshold){
		      
		      if ( kineticenergy[i] > 0 && j == 0 && k ==0) {
			pOutputFile->Set_KE_f(i, kineticenergy[i]); 
			pOutputFile->Set_ToF_f(i, VDC_time[i]);
			pOutputFile->Set_Charge_f(i, charge[i]);
			//cerr << 			  
			//  "KE[" << i << "] = " << kineticenergy[i] << 
			//  " , and ToF = " << VDC_time[i]/ns << "(ns)."  << 
			//  " Particle is: " <<  particle[i] << " ." << 
			//  " Charge is: " << charge[i] << endl;
		      }
		      
		      // Storing largest KE of the event
		      //if ( kineticenergy[i].size() > 0) pOutputFile->Set_KE_f(i, kineticenergy[i].at(0));
		      // Storing largest group KEavg
		      //if( KEavg[i].size() > 0)  pOutputFile->Set_KE_f(i,KEavg[i].at(0));

		      pOutputFile->Set_edep_f(i,j,k,energyTotal[i][j][k]);
		      //position[i][j][k] /= energyTotal[i][j][k];
		      // Storing avg position
		      if (position[i][j][k].getX() != 0. && position[i][j][k].getY()!= 0. && position[i][j][k].getZ() != 0.){
			      position[i][j][k] /=Num[i][j][k];
			      position[i][j][k].setY(0.);}
		    }
	      }	
	  }
	}
	*/

  // ----------------- Hodoscope -----------------------

  G4int PADID = G4SDManager::GetSDMpointer()->GetCollectionID("Paddle/PaddleHitCollection");
  PaddleHitCollection * paddleHC = (PaddleHitCollection *) (evt->GetHCofThisEvent()->GetHC(PADID));
  G4int totalPadHits = paddleHC->entries();
  PaddleHit * pHit;
  // sort the hits in order of time
  // temporary arrays to store time and hit number
	G4double * hit_time = new G4double[totalPadHits];
	G4int * hit_num = new G4int[totalPadHits];
	//G4cout << "Hodoscope total hits = " << totalPadHits << G4endl;
	G4int valid_hits = 0;
	for(G4int jj = 0; jj < totalPadHits; jj++) //for all hits
		{
		pHit = (*paddleHC)[jj];
		//pHit->Print();
		if( pHit->GetEdep() > 0.)
			{
			G4double now = pHit->GetTime();
			hit_time[valid_hits] = now;
			hit_num[valid_hits] = jj;
			for(G4int kk = valid_hits-1; kk >= 0; kk--)
				{
				if(now < hit_time[kk])
					{
					hit_time[kk+1] = hit_time[kk];
					hit_time[kk] = now;
					hit_num[kk+1] = hit_num[kk];
					hit_num[kk] = jj;
					} 
				}
			valid_hits++;
			}
		}
	// set hits to zero
	//G4cout << "Valid Hits = " << valid_hits << G4endl;
	G4double PadEnergy[2][NUMPADDLES], PadLight[2][NUMPADDLES], PadTime[2][NUMPADDLES];
	G4bool pad_hit[2][NUMPADDLES], hod_hit[2];



	//Store paddle KE (PadKE), and counts of paddle KE(PadKENum)
	G4double PadKE[2], PadKENum[2];
	for(G4int i = 0; i < 2; i++)
		{
		hod_hit[i] = false;




		PadKE[i] = 0.; PadKENum[i] = 0.;
		for(G4int j = 0; j < NUMPADDLES; j++)
			{
			PadEnergy[i][j] = PadLight[i][j] = PadTime[i][j] = 0;
			pad_hit[i][j] = false;
			}
		}
	for(G4int jj = 0; jj < valid_hits; jj++) //for all valid hits sorted in time
		{
		G4int hit_jj = hit_num[jj];
		pHit = (*paddleHC)[hit_jj];
		//G4cout << "Valid Hit Number "<<jj<<" Hit Number " << hit_jj << G4endl;
		//pHit->Print();
		G4int hod = pHit->GetHodoscopeNB();
		G4int pad = pHit->GetPaddleNB();
		PadEnergy[hod][pad] += pHit->GetEdep();
		PadLight[hod][pad] += pHit->GetLight();




		// May change this value: only when KE at each step is bigger than 2MeV, KE of the step will be stored
		if( (pHit->GetKE()) >= 2.0*MeV){
		  PadKE[hod] += pHit->GetKE();
		  PadKENum[hod]++;}
		
		// Below invalidate in vacuum
		if(PadLight[hod][pad] >= fPaddleThreshold && PadTime[hod][pad] == 0.)
		  //if (pHit->GetKE() >= 2.0*MeV)
			{
			  PadTime[hod][pad] = pHit->GetTime();
			  pad_hit[hod][pad] = true;
			}
		}
	delete[] hit_time;
	delete[] hit_num;

  // add some statistical noise to light output and see if we still have a hit
  G4double lightOutSmear;
  for(G4int hod = 0; hod < 2; hod++)
	for(G4int pad = 0; pad < NUMPADDLES; pad++)
		{
		  if(pad_hit[hod][pad])
		  {
		    //Units are carried by sigmaLightStatistical and sigmaLightNoise
		    //Statistical -- proportional to sqrt(light_output)
		    lightOutSmear = (sqrt(PadLight[hod][pad]/MeV)) * G4RandGauss::shoot(0.0, fsigmaLightStatistical);
		    //Noise -- constant contribution
		    lightOutSmear += G4RandGauss::shoot(0.0, fsigmaLightNoise);
		    PadLight[hod][pad] += lightOutSmear;


		    //use below line when in vacuum
		    //if( PadTime[hod][pad] == 0.)
		    if(PadLight[hod][pad] < fPaddleThreshold)
		      pad_hit[hod][pad] = false;
		    else hod_hit[hod] = true;
		  }
		}

  for(G4int hod = 0; hod < 2; hod++)
	{
	pOutputFile->Set_hod_hit(hod, hod_hit[hod]);
	if(hod_hit[hod])
		{
		  //Store avg. kinetic energy in every hit collection step that
		  // the step with KE > 2.0MeV
		  if( PadKE[hod]/PadKENum[hod]>2.0*MeV &&  PadKE[hod]/PadKENum[hod] < 60.0*MeV ){
		    pOutputFile->Set_pad_KE(hod, PadKE[hod]/PadKENum[hod]);
		    //cerr << "Number " << hod << " hit. With KE " <<  PadKE[hod]/PadKENum[hod] << endl;
		  }
		for(G4int pad = 0; pad < NUMPADDLES; pad++)
			{
			pOutputFile->Set_pad_hit(hod, pad, pad_hit[hod][pad]);
			if(pad_hit[hod][pad])
				{
				pOutputFile->Set_pad_energy(hod, pad, PadEnergy[hod][pad]);
				pOutputFile->Set_pad_light(hod, pad, PadLight[hod][pad]);
				//cerr << "PadLight[" << hod << "][" << pad << "] = " << PadLight[hod][pad]/MeV << endl;
				pOutputFile->Set_pad_time(hod, pad, PadTime[hod][pad]);
				}
			}
		}
	}

  

  // ----------------------------------------------------------
  // VDC position, angle information
  // Now since all 4 wire plane may have multiple wire hit data
  // x,y in each wire plane now is ambiguous, temporary commented
  pOutputFile->Set_event_number(event_number);
  /*
  for(G4int det = 0; det < 2; det++)
	{
  	pOutputFile->Set_detector_package(det, det_hit[det]);
  	if(det_hit[det])
		{
	  	G4ThreeVector VDC0 = (position[det][0][0] + position[det][0][1])/2.; 
	  	G4ThreeVector VDC1 = (position[det][1][0] + position[det][1][1])/2.; 
	  	// G4cout.precision(4);
	  	// G4cout << " Position VDC0 : "
	  	//   << std::setw(7) << G4BestUnit( VDC0,"Length") <<G4endl; 
	  	// G4cout << " Position VDC1 : "
	  	//   << std::setw(7) << G4BestUnit( VDC1,"Length") <<G4endl; 
		
	  	// change to RAYTRACE coordinate system for output
	  	// and add a smearing of the output consistent with the VDC resolution
	  	G4double x_out = -VDC0.getX() + CLHEP::RandGauss::shoot(0., fVDCSigma);
	  	G4double y_out = VDC0.getZ() + CLHEP::RandGauss::shoot(0., fVDCSigma);
	  	G4double x2_out = -VDC1.getX() + CLHEP::RandGauss::shoot(0., fVDCSigma);
	  	G4double y2_out = VDC1.getZ() + CLHEP::RandGauss::shoot(0., fVDCSigma);
	  	G4double theta_out = atan2((x2_out - x_out), fVDCSpacing)*rad;
	  	G4double phi_out = atan2((y2_out - fVDCSpacing - y_out), fVDCSpacing)*rad;
	  	//G4cout << "Begin setting variables" << G4endl;
	  	//G4cout << "x_out = " << x_out/cm << " cm" << G4endl;
	  	//G4cout << "y_out = " << y_out/cm << " cm" << G4endl;
	  	//G4cout << "theta_out = " << theta_out/mrad << " mrad" << G4endl;
	  	//G4cout << "phi_out = " << phi_out/mrad << " mrad" << G4endl;
		

		///////////////////////////////////////

                // new choose for data coordinate for output.
                // theta and phi are angle between trajectory of electrons (positrons) and global coordinate system 
                // where theta and phi here are pi/2 minus ordinary spherical coordinate theta or phi in global coordinate system (+x downstream direction, +z up direction)
                // x,y are also measured same direction as global coordinate system, but relative to center of 1st VDC for each detector package.
		G4double x_f = x_out;
		G4double y_f = y_out/sqrt(2);
		G4double theta_f = atan2((y_out-y2_out)/sqrt(2),sqrt((x_out-x2_out)*(x_out-x2_out)+2*(fVDCSpacing+(y_out-y2_out)/2)*(fVDCSpacing+(y_out-y2_out)/2)))*rad;
		G4double phi_f = (atan2((x_out-x2_out),sqrt(2)*(fVDCSpacing+(y_out-y2_out)/2)))*rad;
					 
					 x_out = x_f;
					 y_out = y_f;
					 theta_out = theta_f;
					 phi_out = phi_f;
					 //////////////////////

					 pOutputFile->Set_x_f(det, x_out);
					 pOutputFile->Set_y_f(det, y_out);
					 pOutputFile->Set_theta_f(det, theta_out);
					 pOutputFile->Set_phi_f(det, phi_out);
					 
		}
	}*/

  
  // If we are using the monitor - always write event
  // else write the event if there is a hit in a VDC
  // but not if we are requiring a hit in the hodoscope and there is no
  // hit in the corresponding hodoscope
  if(use_monitor)
    {
      //If the following line is added, monitors will write data when all detectors are hit.
      /*
      if(det_hit[0] && hod_hit[0] && det_hit[1] && hod_hit[1])
	{
	  cerr<< det_hit[0] << hod_hit[0] << det_hit[1] << hod_hit[1] << endl;
	  pOutputFile->WriteEvent();
	}
      */
      //If following lines are used, such data will be stored if at least one of the SD is hit.
      //if(det_hit[0] && hod_hit[0] && det_hit[1] && hod_hit[1])
      if(det_hit[0] || hod_hit[0] || det_hit[1] || hod_hit[1])
	//Below is cut closer to real experiment, only trigger data acquisition by scintillater signal
	//if(hod_hit[0] || hod_hit[1])
	{
	  //cerr<< det_hit[0] << hod_hit[0] << det_hit[1] << hod_hit[1] << endl;
	  //cerr << hod_hit[0] << hod_hit[1] << endl;
	  pOutputFile->WriteEvent();
	}
    }

  /*
  else
	{
	if( fRequireHodoscopeHit )
		{
		if((det_hit[0] && hod_hit[0]) || (det_hit[1] && hod_hit[1]) )
  			pOutputFile->WriteEvent();
		}
	else
		{
		if( det_hit[0] || det_hit[1])
  			pOutputFile->WriteEvent();
		}
	}
  
 

  //Rob's code
   if(use_monitor)
    {
	  pOutputFile->WriteEvent();
    }

  
  else
	{
	if( fRequireHodoscopeHit )
		{
		if((det_hit[0] && hod_hit[0]) || (det_hit[1] && hod_hit[1]) )
  			pOutputFile->WriteEvent();
		}
	else
		{
		if( det_hit[0] || det_hit[1])
  			pOutputFile->WriteEvent();
		}
	}

  */
  return;

	}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
