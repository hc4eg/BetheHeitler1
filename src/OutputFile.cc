#include "OutputFile.hh"
#include "OutputFileMessenger.hh"
//#include "G4ios.hh"
#include <fstream>
using namespace std;
class OutputFile;

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
		fclose(fd);
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
	sprintf(filename, "%s.%.4d.dat", fFilename_prefix.c_str(), frunno);
	fd = fopen(filename, "w");
	if(fd == NULL)
		{
		G4cout << "ERROR: OutputFile::Openfile(): Cannot open file " << G4String(filename)
			<< " for output." << G4endl;
		return false;
		}

	foutput_lines = 0;
	fFile_open = true;
	G4cout << "File " << G4String(filename) << " opened for output." << G4endl;
	return true;
	}
// Close a file
void
OutputFile::CloseFile()
	{
	if(fFile_open)
		{
		fclose(fd);
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

	fprintf(fd,"Event: %8d\n",fevent_number);
	fprintf(fd,"Input:");
	fprintf(fd," %13.4g",fenergy_i/MeV);      
	fprintf(fd," %13.4g",fdelta_i);      
	fprintf(fd," %13.4g",fx_i/cm);      
	fprintf(fd," %13.4g",fy_i/cm);      
	fprintf(fd," %13.4g",ftheta_i/mrad);      
	fprintf(fd," %13.4g\n",fphi_i/mrad);      
	if(fUseMonitor)
		{
		for(G4int i = 0; i < 2; i++)
			{
			if(fMonitorHit[i])
				{
				fprintf(fd,"Monitor:");
				if(i == 1) fprintf(fd," e-");
				else fprintf(fd," e+");
				fprintf(fd," %13.4g",fMonitorKineticEnergy[i]/MeV);      
				fprintf(fd," %13.4g",fMonitorX[i]/cm);      
				fprintf(fd," %13.4g",fMonitorY[i]/cm);      
				fprintf(fd," %13.4g",fMonitorTheta[i]/mrad);      
				fprintf(fd," %13.4g\n",fMonitorPhi[i]/mrad);      
				}
			}
		}
	for(G4int i = 0; i < 2; i++)
		{
		if(fdetector_package[i])
			{
			fprintf(fd,"Detector: %d\n",i);
			fprintf(fd,"VDC:");
			fprintf(fd," %13.4g",fx_f[i]/cm);      
			fprintf(fd," %13.4g",fy_f[i]/cm);      
			fprintf(fd," %13.4g",ftheta_f[i]/mrad);      
			fprintf(fd," %13.4g",fphi_f[i]/mrad);      
			fprintf(fd,"\n");
			if(fHod_hit[i])
			   {
			   for(G4int j = 0; j < 29; j++)
				{
				if(fPad_hit[i][j])
					{
					fprintf(fd,"Paddle: %d ", j);
					fprintf(fd," %13.4g",fPadEnergy[i][j]/MeV);      
					fprintf(fd," %13.4g",fPadLight[i][j]/MeV);      
					fprintf(fd," %13.4g",fPadTime[i][j]/ns);      
					fprintf(fd,"\n");
					}
				}
			   }
			}
		}
	//G4cout << "OutputFile::WriteEvent(): done." << G4endl;
	foutput_lines++;
	}
void
OutputFile::WriteComment(G4String line)
	{
	fprintf(fd,"# %s\n",line.c_str());
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
