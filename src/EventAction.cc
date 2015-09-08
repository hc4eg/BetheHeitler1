#include "EventAction.hh"

#include "RunAction.hh"
//#include "EventActionMessenger.hh"

#include "G4Event.hh"
#include "G4UnitsTable.hh"
#include "Randomize.hh"
#include <iomanip>


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
  for(G4int i; i < 2; i++)
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
	
  if(mon_hit[0] || mon_hit[1])
	{
  	for(ipart = 0; ipart < 2; ipart++)
		{
		pOutputFile->SetMonitorHit(ipart, mon_hit[ipart] );
		if(mon_hit[ipart])
			{
			mHit = (*MonitorHC)[min_time_hit[ipart]];
			// Set into output file
			// Translating to RAYTRACE coordinates
			pOutputFile->SetMonitorKineticEnergy(ipart, mHit->GetKineticEnergy() );
			G4ThreeVector pos = mHit->GetGlobalPrePosition();
			pOutputFile->SetMonitorX(ipart, pos.y() );
			pOutputFile->SetMonitorY(ipart, pos.z() );
			G4ThreeVector dir = mHit->GetMomentumDirection();
			G4double theta_m = atan2( dir.y(), dir.x() );
			G4double phi_m = atan2( dir.z(), dir.x() );
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
  G4ThreeVector pos;
  for(G4int i = 0; i < 2; i++)
  for(G4int j = 0; j < 2; j++)
  for(G4int k = 0; k < 2; k++)
	{ position[i][j][k] = G4ThreeVector(0.,0.,0.); energyTotal[i][j][k] = 0.;}
  for(G4int jj = 0; jj < totalHits; jj++) //for all hits
	{
	aHit = (*chamberHC)[jj];
//	aHit->Print();
	G4int chamber = aHit->GetVDCnumber();
	G4int layer = aHit->GetVDClayer();
	G4int package = aHit->GetVDCpackage();
	G4double edep = aHit->GetEdep();
	if(edep > 0.)
		{
		energyTotal[package][chamber][layer] += edep;
		pos = (aHit->GetLocalPrePosition() + aHit->GetLocalPostPosition() )/2.;
		position[package][chamber][layer] += pos*edep;
		}
	}
  G4bool det_hit[2];
  for(G4int i = 0; i < 2; i++)
	{
	det_hit[i] = true;
  	for(G4int j = 0; j < 2; j++)
  	for(G4int k = 0; k < 2; k++)
		{
		if(energyTotal[i][j][k] < 0.2*keV)
			{ det_hit[i] = false;}
		else
			{
			position[i][j][k] /= energyTotal[i][j][k];
			position[i][j][k].setY(0.);
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
	for(G4int i = 0; i < 2; i++)
		{
		hod_hit[i] = false;
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
		if(PadLight[hod][pad] >= fPaddleThreshold && PadTime[hod][pad] == 0.)
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
			if(PadLight[hod][pad] < fPaddleThreshold) pad_hit[hod][pad] = false;
			else hod_hit[hod] = true;
			}
		}

  for(G4int hod = 0; hod < 2; hod++)
	{
	pOutputFile->Set_hod_hit(hod, hod_hit[hod]);
	if(hod_hit[hod])
		{
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
  	pOutputFile->WriteEvent();
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
  return;

}  

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
