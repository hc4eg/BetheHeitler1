#include "OutputFile.hh"
#include "OutputFileMessenger.hh"
#include "PrimaryGeneratorAction.hh"
#include "G4ios.hh"
#include <fstream>

using namespace std;
class OutputFile;
#include "TBranch.h"
ClassImp(BH_Event)

// initialize pointer
OutputFile * OutputFile::pOutputFilePointer = 0;

// private constructor
OutputFile::OutputFile()
	{
	// just some default values
	fFilename_prefix = "output";
	fFile_open = false;
	fuse_unique_filename = false;
	flines_for_heading = -1; // -1 = no additional headings
	fUseMonitor = false;
	out_pair_mode = false;
	// Start a messenger
	fmessenger = new OutputFileMessenger( this);
	}
OutputFile::~OutputFile()
	{
	if(fFile_open)
		{
		  //fclose(fd);
		fFile_open = false;
		}
	pOutputFilePointer = 0;
	}
OutputFile * OutputFile::GetOutputFilePointer()
	{
	if(pOutputFilePointer == 0)
		{
		pOutputFilePointer = new OutputFile();
		}
	return pOutputFilePointer;
	}
// Open a file
G4bool
OutputFile::OpenFile(G4int runno)
	{
	if(fFile_open)
		{
		G4cout << "WARNING: OutputFile::OpenFile():"
		<<" File " << G4String(filename) << " already open." << G4endl;
		return true;
		}
	frunno = runno;
	if(fuse_unique_filename)
		{
		ifstream runin;
		ofstream runout;
		// if the file "<prefix>.runno" exists - use its contents as the
		// output file run number
		//G4cout << "Using unique filename." << G4endl;
		G4String infile = fFilename_prefix+".runno";
		runin.open(infile);
		if(runin.good())
			{
			runin>>frunno;
			//G4cout << "Reading file " << infile << " runno = " << frunno << G4endl;
			runin.close();
			}
		else
			{
			frunno = 0;
			}
		frunno++;
		runout.open(infile, ios::trunc);
		if(runout.good())
			{
			runout << frunno << G4endl;
			runout.close();
			}
		else
			{
			G4cout << "Cannot open " << infile << " for output." << G4endl;
			}
		}
	//sprintf(filename, "%s.%.4d.dat", fFilename_prefix.c_str(), frunno);
	//sprintf(filename, "%s.%.4d.root", fFilename_prefix.c_str(), frunno);
	//fd = fopen(filename, "w");
	/*
	if(fd == NULL)
		{
		G4cout << "ERROR: OutputFile::Openfile(): Cannot open file " << G4String(filename)
			<< " for output." << G4endl;
		return false;
		}
	*/
	//cerr << "got here 1" << endl;
	TString fname = Form("%s.%.4d.root", fFilename_prefix.c_str(), frunno);
	
	bh_tree_file = new TFile(fname,"recreate");
	if(bh_tree_file == NULL) {
	  G4cerr << "ERROR: OutputFile::Openfile(): Cannot open file " << fname
			<< " for output." << G4endl;
		return false;
	}
	cerr << "The root file: " << fname << " is opened." << endl;

	tree = new TTree("T","root_tree");

	event = new BH_Event();

	branch = (TBranch*)tree->Branch("B","BH_Event",&event);

	foutput_lines = 0;
	fFile_open = true;
	G4cout << "File " << fname << " opened for output." << G4endl;
	return true;
	
	}
// Close a file
void
OutputFile::CloseFile()
	{	  
	 
	if(fFile_open)
		{
		  //fclose(fd);
	
        	bh_tree_file->Write();
		delete tree;
		fFile_open = false;
		//G4cout << "Output File " << fname << " closed." <<G4endl;
		}
	else
		{
		G4cout << "WARNING: OutputFile::CloseFile(): File not open." << G4endl;
		}
	}
	
// Write an event
void 
OutputFile::WriteEvent()
	{
	  //event = new BH_Event;
	  //event->event_num = (Int_t)fevent_number;
	  //cerr << "Event number " << event->event_num << "has been read" << endl;
          
	  //branch->SetAddress(&event);

	  //tree->Fill();
	  
	  /*
	//G4cout << "OutputFile::WriteEvent(): Called." << G4endl;
	if(!fFile_open)
		{
		G4cout << "ERROR: OutputFile::WriteEvent(): Cannot write event. No file open." << G4endl;
		return;
		}
	if(foutput_lines == 0 || ( flines_for_heading > 0 && (foutput_lines%flines_for_heading == 0) ))
		{
		fprintf(fd,"# ---- Output Format ----\n");
		fprintf(fd,"#  Event: event_number\n");
		fprintf(fd,"#  Input: ");
		fprintf(fd," energy_i(MeV)");      
		fprintf(fd,"    delta_i(%%)");      
		fprintf(fd,"       x_i(cm)");      
		fprintf(fd,"       y_i(cm)");      
		fprintf(fd," theta_i(mrad)");      
		fprintf(fd,"   phi_i(mrad)");      
		fprintf(fd, "\n");
		if(fUseMonitor)
			{
			fprintf(fd,"#  Monitor: <e-/e+>");
			fprintf(fd," Energy(MeV)");      
			fprintf(fd,"       x_m(cm)");      
			fprintf(fd,"       y_m(cm)");      
			fprintf(fd," theta_m(mrad)");      
			fprintf(fd,"   phi_m(mrad)");      
			fprintf(fd, "\n");
			}
		fprintf(fd,"#  Detector: detector_package (0 = beam left, 1 = beam right)\n");
		fprintf(fd,"#  VDC:");
		fprintf(fd,"       x_f(cm)");      
		fprintf(fd,"       y_f(cm)");      
		fprintf(fd," theta_f(mrad)");      
		fprintf(fd,"   phi_f(mrad)");      
		fprintf(fd,"\n");
		fprintf(fd,"#  Paddle: paddle_number");
		fprintf(fd,"       Edep(MeV)");      
		fprintf(fd,"  Light_out(MeV)");      
		fprintf(fd,"    Hit_time(ns)");      
		fprintf(fd,"\n");
		fprintf(fd,"#    .... [possible more Paddle lines]\n");      
		fprintf(fd,"#    .... [possible repeat for other detector package]\n");      
		fprintf(fd,"# --------\n");
		}
	  */

	  //fprintf(fd,"Event: %8d\n",fevent_number); 
	  event->ENum = (Int_t)fevent_number;
	  //Input data
	  //After Input change to vector to store e-,e+ data for pair_mode
	  if (out_pair_mode){
	        
		event->I0.Energy = (Float_t)(fenergy_i[0]/MeV);  
		event->I0.Delta = (Float_t)(fdelta_i[0]);
		event->I0.X = (Float_t)(fx_i[0]/cm);
		event->I0.Y = (Float_t)(fy_i[0]/cm);
		event->I0.Theta = (Float_t)(ftheta_i[0]/mrad);
		event->I0.Phi = (Float_t)(fphi_i[0]/mrad);

		event->I1.Energy = (Float_t)(fenergy_i[1]/MeV);  
		event->I1.Delta = (Float_t)(fdelta_i[1]);
		event->I1.X = (Float_t)(fx_i[1]/cm);
		event->I1.Y = (Float_t)(fy_i[1]/cm);
		event->I1.Theta = (Float_t)(ftheta_i[1]/mrad);
		event->I1.Phi = (Float_t)(fphi_i[1]/mrad);
	  }
	  else{
	    event->I0.Energy = (Float_t)(fenergy_i[0]/MeV);
	    //fprintf(fd," %13.4g",fdelta_i);   
	    event->I0.Delta = (Float_t)(fdelta_i[0]);	 
	    //fprintf(fd," %13.4g",fx_i/cm);
	    event->I0.X = (Float_t)(fx_i[0]/cm);
	    //fprintf(fd," %13.4g",fy_i/cm);
	    event->I0.Y = (Float_t)(fy_i[0]/cm);
	    //fprintf(fd," %13.4g",ftheta_i/mrad);
	    event->I0.Theta = (Float_t)(ftheta_i[0]/mrad);
	    //fprintf(fd," %13.4g\n",fphi_i/mrad);
	    event->I0.Phi = (Float_t)(fphi_i[0]/mrad);
	  }



	  // Monitor data
	if(fUseMonitor)
		{
		  //Check if there is single e+ or e- produced, or e+/e- pair produced
		  // Only monitor 0 is hit, set data for monitor 1 to be all 0
		  if(fMonitorHit[0] && !fMonitorHit[1])
		    {
		      		  Monitor* moni = new Monitor;
				  moni->Charge = 1.0;
				  moni->Energy = (Float_t)(fMonitorKineticEnergy[0]/MeV);
				  //cerr << "Ep = " << fMonitorKineticEnergy[0] << endl; 
				  moni->X = (Float_t)(fMonitorX[0]/cm);
				  moni->Y = (Float_t)(fMonitorY[0]/cm);     
				  moni->Theta = (Float_t)(fMonitorTheta[0]/mrad);
				  moni->Phi = (Float_t)(fMonitorPhi[0]/mrad);
				  event->M0.SetMonitor(moni);
				  delete moni;
				  
				  Monitor* moni1 = new Monitor;
				  moni1->Charge = 0.0;
				  moni1->Energy = (Float_t)(0.0/MeV); 
				  moni1->X = (Float_t)(0.0/cm);
				  moni1->Y = (Float_t)(0.0/cm);     
				  moni1->Theta = (Float_t)(0.0/mrad);
				  moni1->Phi = (Float_t)(0.0/mrad);
				  event->M1.SetMonitor(moni1);
				  delete moni1;
		    }
		  // Only monitor 1 is hit, set data for monitor 0 to be all 0
		  else if(!fMonitorHit[0] && fMonitorHit[1])
		    {
		      		  Monitor* moni = new Monitor;
				  moni->Charge = -1.0;
				  moni->Energy = (Float_t)(fMonitorKineticEnergy[1]/MeV);
				  //cerr << "Ee = " << fMonitorKineticEnergy[1] << endl; 
				  moni->X = (Float_t)(fMonitorX[1]/cm);
				  moni->Y = (Float_t)(fMonitorY[1]/cm);     
				  moni->Theta = (Float_t)(fMonitorTheta[1]/mrad);
				  moni->Phi = (Float_t)(fMonitorPhi[1]/mrad);
				  event->M1.SetMonitor(moni);
				  delete moni;

				  Monitor* moni1 = new Monitor;
				  moni1->Charge = 0.0;
				  moni1->Energy = (Float_t)(0.0/MeV); 
				  moni1->X = (Float_t)(0.0/cm);
				  moni1->Y = (Float_t)(0.0/cm);     
				  moni1->Theta = (Float_t)(0.0/mrad);
				  moni1->Phi = (Float_t)(0.0/mrad);
				  event->M0.SetMonitor(moni1);
				  delete moni;				  
		    }
		  else
		    for(G4int i = 0; i < 2; i++)
		      {
			//if(fMonitorHit[i])
				{
				  Monitor* moni = new Monitor;
				  //fprintf(fd,"Monitor:");
				  //if(i == 1) fprintf(fd," e-");
				  if(i==1) moni->Charge = -1.0;
				  //else fprintf(fd," e+");
				  else moni->Charge = 1.0;
				  //fprintf(fd," %13.4g",fMonitorKineticEnergy[i]/MeV);
				  moni->Energy = (Float_t)(fMonitorKineticEnergy[i]/MeV);   
				  moni->X = (Float_t)(fMonitorX[i]/cm);
				  //fprintf(fd," %13.4g",fMonitorY[i]/cm); 
				  moni->Y = (Float_t)(fMonitorY[i]/cm);
				  //fprintf(fd," %13.4g",fMonitorTheta[i]/mrad);      
				  moni->Theta = (Float_t)(fMonitorTheta[i]/mrad);
				  //fprintf(fd," %13.4g\n",fMonitorPhi[i]/mrad);
				  moni->Phi = (Float_t)(fMonitorPhi[i]/mrad);
				  if (i == 0 )
				    event->M0.SetMonitor(moni);
				  else
				    event->M1.SetMonitor(moni);
				  delete moni;
				}
			}
		}



	// Detector data
	/*
	event->D0.V0U.clear();
	event->D0.V0V.clear();
	event->D0.V1U.clear();
	event->D0.V1V.clear();
	event->D1.V0U.clear();
	event->D1.V0V.clear();
	event->D1.V1U.clear();
	event->D1.V1V.clear();*/
	/*
	for(G4int n = 0; n < 4; n++){
	  event->D0.WP[n].clear();
	  event->D1.WP[n].clear();}*/

	event->D0.W.clear();
	event->D1.W.clear();
	event->D0.P.clear();
	event->D1.P.clear();
	event->HKE1.clear();
	event->HKE0.clear();
	for(G4int i = 0; i < 2; i++)
		{

		  //Wire Hit data
		  //cerr << "Event Number" << event->ENum << endl;
		  if(fdetector_package[i]){
		    for(G4int j = 0; j < 2; j++)
		      for(G4int k = 0; k < 2; k++){
			for(G4int l = 0; l < Get_NumHit_f(i,j,k); l++){
			  OutputWire OutWire = Get_Wire_f(i,j,k,l);
			  Wire Wire;
			  // Wire hit data for a single hit wire.
			  Wire.WirePlane = (Int_t)OutWire.Get_WirePlane_f();
			  Wire.WireNum = (Int_t)OutWire.Get_WireNum_f();
			  Wire.X = (Float_t)(OutWire.Get_X_f()/cm);
			  Wire.Y = (Float_t)(OutWire.Get_Y_f()/cm);
			  Wire.KE = (Float_t)(OutWire.Get_KE_f()/MeV);
			  Wire.Edep = (Float_t)(OutWire.Get_Edep_f()/MeV);
			  Wire.ToF = (Float_t)(OutWire.Get_ToF_f()/ns);
			  Wire.Charge = (Float_t)OutWire.Get_Charge_f();
			  Wire.Particle = (Int_t)OutWire.Get_Particle_f();
			  //cerr << "[" << i << "][" << j << "][" << k << "] have " << Get_NumHit_f(i,j,k) << " hits." << endl;
			  
			  //In order to use TTree::MakeClass() Method, change different naming
			  //Of wireplane into wireplane array:
			  /*
			  // Assign wire to corresponding wireplane in root file
			  if( i == 0 && j == 0 && k == 0)event->D0.V0U.push_back(Wire);
			  else if( i == 0 && j == 0 && k == 1)  event->D0.V0V.push_back(Wire);
			  else if( i == 0 && j == 1 && k == 0)  event->D0.V1U.push_back(Wire);
			  else if( i == 0 && j == 1 && k == 1)  event->D0.V1V.push_back(Wire);
			  else if( i == 1 && j == 0 && k == 0)  event->D1.V0U.push_back(Wire);
			  else if( i == 1 && j == 0 && k == 1)  event->D1.V0V.push_back(Wire);
			  else if( i == 1 && j == 1 && k == 0)  event->D1.V1U.push_back(Wire);
			  else if( i == 1 && j == 0 && k == 1)  event->D1.V1V.push_back(Wire);*/
			  //It seems root file can't store and show in TBrowser with array of vector data structure
			  /*
			  if( i == 0 && j == 0 && k == 0)event->D0.WP[0].push_back(Wire);
			  else if( i == 0 && j == 0 && k == 1)  event->D0.WP[1].push_back(Wire);
			  else if( i == 0 && j == 1 && k == 0)  event->D0.WP[2].push_back(Wire);
			  else if( i == 0 && j == 1 && k == 1)  event->D0.WP[3].push_back(Wire);
			  else if( i == 1 && j == 0 && k == 0)  event->D1.WP[0].push_back(Wire);
			  else if( i == 1 && j == 0 && k == 1)  event->D1.WP[1].push_back(Wire);
			  else if( i == 1 && j == 1 && k == 0)  event->D1.WP[2].push_back(Wire);
			  else if( i == 1 && j == 1 && k == 1)  event->D1.WP[3].push_back(Wire);*/
			  
			  //Assign WirePlane information inside hit wires now:
			  if( i == 0) event->D0.W.push_back(Wire);
			  else if( i == 1) event->D1.W.push_back(Wire);
			}
		      }}

		// Hodoscope data
			if(fHod_hit[i])
			   {
			     if (i == 1 && fPadKE[1] > 0.*MeV && fPadKE[1] < 60.0*MeV) {
			       event->HKE1.push_back((Double_t)fPadKE[1]); 
			       // cerr << "HKE1 = " << event->HKE1.at(0) << endl;
			     }
			     else if(i == 0 && fPadKE[0] > 0.*MeV && fPadKE[0] < 60.0*MeV) {
			       event->HKE0.push_back((Double_t)fPadKE[0]); 
			       //  cerr << "HKE0 = " << event->HKE0.at(0) << endl;
			     }

			   for(G4int j = 0; j < 29; j++)
				{
				if(fPad_hit[i][j])
					{
					  Paddle paddle;
					  //fprintf(fd,"Paddle: %d ", j);
					  paddle.PNum = (Int_t)j;
					  //fprintf(fd," %13.4g",fPadEnergy[i][j]/MeV);  
					  paddle.Edep = (Float_t)(fPadEnergy[i][j]/MeV);
					  //fprintf(fd," %13.4g",fPadLight[i][j]/MeV);
					  paddle.Light = (Float_t)(fPadLight[i][j]/MeV);
					  //cerr << "paddle.Light[" << i << "][" << j << "] =" << fPadLight[i][j]/MeV << endl;
					  //fprintf(fd," %13.4g",fPadTime[i][j]/ns);      
					  paddle.Time = (Float_t)(fPadTime[i][j]/ns);
					  //fprintf(fd,"\n");
					  if(i==0) event->D0.P.push_back(paddle);
					  else event->D1.P.push_back(paddle);
					}
				}
			   }
		}

			//G4cout << "OutputFile::WriteEvent(): done." << G4endl;
			//foutput_lines++;
			tree->Fill();
	}
	


void
OutputFile::WriteComment(G4String line)
	{
	  //fprintf(fd,"# %s\n",line.c_str());
	}
	
void
OutputFile::PrintParameters()
	{
	G4cout << "   ----- Output File Parameters -----" << G4endl;
	G4cout << "Filename prefix: " << fFilename_prefix << G4endl;
	G4cout << "Filename will be " << fFilename_prefix << ".<runno>.dat" << G4endl;
	if(fuse_unique_filename)
		G4cout << "Filename uses runno = number read from file: " << fFilename_prefix << ".runno" << G4endl;
	else
		G4cout << "Filename uses runno = G4 run number" << G4endl;
	if(fUseMonitor)
		G4cout << "Monitor information will be recorded and included in output file." << G4endl;
	else
		G4cout << "Monitor information will not be recorded." << G4endl;	
	}
