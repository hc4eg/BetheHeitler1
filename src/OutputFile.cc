#include "OutputFile.hh"
#include "OutputFileMessenger.hh"
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
cerr << "got here 1" << endl;
	TString fname = Form("%s.%.4d.root", fFilename_prefix.c_str(), frunno);
	
	bh_tree_file = new TFile(fname,"recreate");
	if(bh_tree_file == NULL) {
	  G4cerr << "ERROR: OutputFile::Openfile(): Cannot open file " << fname
			<< " for output." << G4endl;
		return false;
	}
cerr << "The root file: " << fname << " is opened." << endl;
cerr << "got here 2" << endl;
//     TTree* T = new TTree("T","root tree");	
      //        tree = (TTree*)bh_tree_file->Get("T");
      tree = new TTree("tree","root_tree");


cerr << "got here 3" << endl; 
 event = new BH_Event();
	//        event = 0;
cerr << "got here 3.1" << endl;

	branch = (TBranch*)tree->Branch("Event_Branch","BH_Event",&event);

//branch = (TBranch*)tree->FindBranch("Event_Branch");

//branch = tree->GetBranch("Event_Branch");	
cerr << "got here 4" << endl;
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
		cerr << "The root file has been written" << endl;
		
		fFile_open = false;
		G4cout << "Output File " << G4String(filename) << " closed." <<G4endl;
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
			fprintf(fd," energy_m(MeV)");      
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
	  event->event_num = (Int_t)fevent_number;
	  cerr << "fevent_number is: " << fevent_number << endl;
	  //cerr << "get event number: " << event->event_num << endl;
	  //fprintf(fd,"Input:");
	  //fprintf(fd," %13.4g",fenergy_i/MeV);
	  //event->input->energy = 30.0;
	  //event->input->energy = (Float_t)30.0;
	    //(Float_t)(fenergy_i);
	  event->input.energy = (Float_t)(fenergy_i/MeV);
	  //cerr << "fenergy_i is: " << (fenergy_i/MeV) << " [MeV]." << endl;
	  //cerr << "fenergy_i is: " << (fenergy_i/keV) << " [keV]." << endl;
	  
	//fprintf(fd," %13.4g",fdelta_i);   
	  event->input.delta = (Float_t)(fdelta_i);
	  //cerr << "fdelta_i is: " << fdelta_i << endl;
   	  //cerr << "get input delta: " << event->input.delta << endl;
	 
	//fprintf(fd," %13.4g",fx_i/cm);      
	  event->input.x = (Float_t)(fx_i/cm);
	  //cerr << "fx_i is: " << fx_i << endl;
	//fprintf(fd," %13.4g",fy_i/cm);
	  event->input.y = (Float_t)(fy_i/cm);
	  //cerr << "fy_i is: " << fy_i << endl;
	//fprintf(fd," %13.4g",ftheta_i/mrad);
	  event->input.theta = (Float_t)(ftheta_i/mrad);
	  //cerr << "ftheta_i is: "<< ftheta_i << endl;
	//fprintf(fd," %13.4g\n",fphi_i/mrad);
	  event->input.phi = (Float_t)(fphi_i/mrad);
	  // cerr << "fphi_i is: "<< fphi_i << endl;

	  
	  //	  fUseMonitor = 1;
	  //cerr << "fUseMonitor is: " << fUseMonitor << endl;
	  event->monitors.clear();
	if(fUseMonitor)
		{
		for(G4int i = 0; i < 2; i++)
			{
			if(fMonitorHit[i])
				{
				  Monitor moni;
				  //fprintf(fd,"Monitor:");
				  //if(i == 1) fprintf(fd," e-");
				  if(i==1) moni.charge = -1;
				  //else fprintf(fd," e+");
				  else moni.charge = 1;
				  //fprintf(fd," %13.4g",fMonitorKineticEnergy[i]/MeV);
				  moni.energy_m = (Float_t)(fMonitorKineticEnergy[i]/MeV);
				  cerr << "fMonitorKineticEnergy is: " << fMonitorKineticEnergy[i];
				  //fprintf(fd," %13.4g",fMonitorX[i]/cm);   
				  moni.x_m = (Float_t)(fMonitorX[i]/cm);
				  //fprintf(fd," %13.4g",fMonitorY[i]/cm); 
				  moni.y_m = (Float_t)(fMonitorY[i]/cm);
				  //fprintf(fd," %13.4g",fMonitorTheta[i]/mrad);      
				  moni.theta_m = (Float_t)(fMonitorTheta[i]/mrad);
				  //fprintf(fd," %13.4g\n",fMonitorPhi[i]/mrad);
				  moni.phi_m = (Float_t)(fMonitorPhi[i]/mrad);
				  event->monitors.push_back(moni);
				}
			}
		}
	  
	event->detector0.vdc.clear();
	event->detector1.vdc.clear();
	event->detector0.paddles.clear();
	event->detector1.paddles.clear();
	for(G4int i = 0; i < 2; i++)
		{
		if(fdetector_package[i])
			{
			  Detector detector;
			  //fprintf(fd,"Detector: %d\n",i);//if (i==0) to detector0
			  VDC dc;
			  //fprintf(fd,"VDC:");
			  //fprintf(fd," %13.4g",fx_f[i]/cm);
			  dc.x_f = (Float_t)(fx_f[i]/cm);
			  //fprintf(fd," %13.4g",fy_f[i]/cm);
			  dc.y_f = (Float_t)(fy_f[i]/cm);
			  //fprintf(fd," %13.4g",ftheta_f[i]/mrad);
			  dc.theta_f = (Float_t)(ftheta_f[i]/mrad);
			  //fprintf(fd," %13.4g",fphi_f[i]/mrad);  
			  dc.phi_f = (Float_t)(fphi_f[i]/mrad);
			  //fprintf(fd,"\n");
			  //cerr << "i: " << i << endl
			  //     << "vdc.x_f: " << dc.x_f <<endl;
			  //<< "y_f: " << dc.y_f <<endl
			  //<< "thta_f:" << dc.theta_f <<endl
			  //<< "phi_f: " << dc.phi_f <<endl;
			  if(i==0) {event->detector0.vdc.push_back(dc);
			    // cerr << "detector0 vdc x_f: " << event->detector0.vdc.at(0).x_f << endl;}
			  }
			  else {event->detector1.vdc.push_back(dc);
			    // cerr << "detector1 vdc x_f: " << event->detector1.vdc.at(0).x_f << endl;}
			  }
			if(fHod_hit[i])
			   {
			   for(G4int j = 0; j < 29; j++)
				{
				if(fPad_hit[i][j])
					{
					  Paddle paddle;
					  //fprintf(fd,"Paddle: %d ", j);
					  paddle.paddle_num = (Int_t)j;
					  //fprintf(fd," %13.4g",fPadEnergy[i][j]/MeV);  
					  paddle.Edep = (Float_t)(fPadEnergy[i][j]/MeV);
					  //fprintf(fd," %13.4g",fPadLight[i][j]/MeV);
					  paddle.Light_out = (Float_t)(fPadLight[i][j]/MeV);
					  //fprintf(fd," %13.4g",fPadTime[i][j]/ns);      
					  paddle.Hit_time = (Float_t)(fPadTime[i][j]/ns);
					  //fprintf(fd,"\n");
					  if(i==0) event->detector0.paddles.push_back(paddle);
					  else event->detector1.paddles.push_back(paddle);
					  //else detector1->paddles->push_back(paddle);
					}
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
	  /*
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
	  */
	}


