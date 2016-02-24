#include "EventAction.hh"
#include "RunAction.hh"
#include "G4Event.hh"
#include "G4UnitsTable.hh"
#include "Randomize.hh"
#include <iomanip>

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

  //G4bool mon_hit[2];
  
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
  G4ThreeVector position[2][2][2]; // position on [package][chamber][layer]
  G4double energyTotal[2][2][2];
  vector<G4double> kineticenergy[2]; 
  //vector<G4double> KEavg[2];
  G4double KEavg[2];
  G4int KENum[2], Num[2][2][2];
  //G4double kineticenergy[2];
  //G4double charge[2];
  //G4String particle[2];
  //G4double VDC_time[2];
  G4ThreeVector pos;
  
  for(G4int i = 0; i < 2; i++)
  for(G4int j = 0; j < 2; j++)
  for(G4int k = 0; k < 2; k++)
    { position[i][j][k] = G4ThreeVector(0.,0.,0.); energyTotal[i][j][k] = 0.; 
      kineticenergy[i].clear(); KENum[i] = 0; Num[i][j][k] = 0;
      KEavg[i] = 0.;
      //KEavg[i].clear();
      //kineticenergy[i] = 0.;
      //charge[i] = 0.;
      //VDC_time[i] = 0;
    }
  
  for(G4int jj = 0; jj < totalHits; jj++) //for all hits
	{
	  aHit = (*chamberHC)[jj];
	  //	aHit->Print();
	  G4int chamber = aHit->GetVDCnumber();
	  G4int layer = aHit->GetVDClayer();
	  G4int package = aHit->GetVDCpackage();
	  G4double edep = aHit->GetEdep();
	  //G4double KE[2] = {0.,0.};


	  if(edep > 0.)
	  {
	    // Storing edep for each wireplane
	    energyTotal[package][chamber][layer] += edep;
	    pos = (aHit->GetLocalPrePosition() + aHit->GetLocalPostPosition() )/2.;
	    // Storing number of hit for each wireplane
	    Num[package][chamber][layer] ++;
	    // Edep weighted sum of position
	    //position[package][chamber][layer] += pos*edep;
	    // Sum of position
	    position[package][chamber][layer] += pos;
	    
	    // Storing KE at each step in hit collection of 1st wireplane in each side
	    if(layer == 0 && chamber== 0 && aHit->GetKE() > 0.){
	      kineticenergy[package].push_back(aHit->GetKE());
	      KENum[package] ++;}
	  }
	  //else kineticenergy[package] = 0.;
	}

  
	//sort kineticenergy in decrement order:
	for (G4int m = 0; m < 2; m++)
	  for (unsigned n = 0; n < kineticenergy[m].size(); n++)
	    for (unsigned l = n; l < kineticenergy[m].size() ; l++ ){
	      if(kineticenergy[m].at(l) > kineticenergy[m].at(n) )
		{ G4double Temp; Temp = kineticenergy[m].at(l); kineticenergy[m].at(l) = kineticenergy[m].at(n); kineticenergy[m].at(n) = Temp;}
	    }
 
	
	//find out the largest group of kineticenergy, epsilon = 0.5*MeV
	for (G4int m = 0; m < 2; m++){
	  G4double avg = 0., epsilon = 0.5*MeV, num = 0.;
	  
	  if ( kineticenergy[m].size() > 0){
	    avg = kineticenergy[m].at(0), num = 0.;
	    if ( kineticenergy[m].size() > 1){
	      for(unsigned n = 1; n < kineticenergy[m].size() ; n++){
		if ( abs( kineticenergy[m].at(n) - avg ) < epsilon ){num ++;avg = ( avg * num + kineticenergy[m].at(n) )/( num + 1.0 );}}}}
	  
	  KEavg[m] = avg;
	  if ( m == 1 && KEavg[1] != 0.) cerr << "D1.V.KE is " << KEavg[1] << " Step Count " << num+1 << ". Out of " << kineticenergy[1].size() << endl; 
	}

	
	/*	
	}
	*/




  G4bool det_hit[2];
  for(G4int i = 0; i < 2; i++)
	{
	det_hit[i] = true;
  	
	for(G4int j = 0; j < 2; j++)
  	for(G4int k = 0; k < 2; k++)
	  {
	    //change following line to deal with vacuum setting.
	    if(energyTotal[i][j][k] < 0.2*keV)
	    //if (kineticenergy[i]/KENum[i] < 2.0*MeV)
	    //if (position[i][j][k].getX() == 0. && position[i][j][k].getY() == 0. && position[i][j][k].getZ() == 0.)
	      { det_hit[i] = false;}
	    else
	      { 
		    
		    if(energyTotal[i][j][k] > 0.2*keV){
		      
		      // Storing largest group KEavg
		      //if( KEavg[i].size() > 0)  pOutputFile->Set_KE_f(i,KEavg[i].at(0));	      
		      if( KEavg[i] > 0)  pOutputFile->Set_KE_f(i,KEavg[i]);

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
		  if( PadKE[hod]/PadKENum[hod]>2.0*MeV &&  PadKE[hod]/PadKENum[hod] < 60.0*MeV){
		    pOutputFile->Set_pad_KE(hod, PadKE[hod]/PadKENum[hod]);
		    //cerr << "Number " << hod << " hit. With KE " <<  PadKE[hod]/PadKENum[hod] << endl;}
		  }
		for(G4int pad = 0; pad < NUMPADDLES; pad++)
			{
			pOutputFile->Set_pad_hit(hod, pad, pad_hit[hod][pad]);
			if(pad_hit[hod][pad])
				{
				pOutputFile->Set_pad_energy(hod, pad, PadEnergy[hod][pad]);
				pOutputFile->Set_pad_light(hod, pad, PadLight[hod][pad]);
				pOutputFile->Set_pad_time(hod, pad, PadTime[hod][pad]);
				}
			}
		}
	}

  // ----------------------------------------------------------
  pOutputFile->Set_event_number(event_number);
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
		


/*

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
*/

					 pOutputFile->Set_x_f(det, x_out);
					 pOutputFile->Set_y_f(det, y_out);
					 pOutputFile->Set_theta_f(det, theta_out);
					 pOutputFile->Set_phi_f(det, phi_out);
					 
		}
	}

  
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
	{
	  cerr<< det_hit[0] << hod_hit[0] << det_hit[1] << hod_hit[1] << endl;
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
