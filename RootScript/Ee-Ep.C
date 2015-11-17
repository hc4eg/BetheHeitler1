#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include "TFile.h"
#include "TTree.h"
#include "TObject.h"   
#include "TClonesArray.h"
#include "Riostream.h"
#include "TCanvas.h"
#include "TH1.h"

using namespace std;

//This root script is used for:
//1.Read Ep, Ee from monitors branch of .root file to make Ee - Ep 1D histogram.
//2.Read Ep, Ee and calculate Theta(opening angle) from monitors branch of .root file to make Ee Ep vs. Theta 2D histogram.

//Declare structures used to store data from a single event
class Monitor: public TObject {
public:
  Monitor() {};
  ~Monitor() {};
    Int_t enumber;
    Int_t charge;
    Float_t energy_m;
    Float_t x_m;
    Float_t y_m;
    Float_t theta_m;
    Float_t phi_m;

  inline Float_t GetEnergy(){return energy_m;}
  inline void ClearMonitor(){
    enumber = 0; charge = 0; energy_m = 0.; x_m = 0.; y_m = 0.; theta_m = 0.; phi_m = 0.;
  }
protected:
  ClassDef(Monitor,2);
  };

class Input: public TObject{
public:
  Input() {};
  ~Input() {};
  Float_t energy;
  Float_t delta;
  Float_t x;
  Float_t y;
  Float_t theta;
  Float_t phi;
  
protected:
  ClassDef(Input,2);
};

class VDC: public TObject{
public:
  VDC() {};
  ~VDC() {};
    Float_t x_f;
    Float_t y_f;
    Float_t theta_f;
    Float_t phi_f;

protected:
  ClassDef(VDC,2);
};

class Paddle: public TObject{
public:
  Paddle() {};
  ~Paddle() {};
  Int_t paddle_num;
  Float_t Edep;
  Float_t Light_out;
  Float_t Hit_time;

protected:
  ClassDef(Paddle,2);
};

//The number of VDC is varying for each Detector package (0 or 1)
//Also number of paddle hit is varying for each Detector package (0 to 29)

class Detector: public TObject{
public:
  //Int_t detector_num;

  Detector() {};
  ~Detector() {};
  vector<VDC> vdc;
  vector<Paddle> paddles;
protected:
  ClassDef(Detector,2);
};

//The number of monitors are varying (1 to 2), since one of e+ or e- or both is recorded in a single event

class BH_Event: public TObject{
public:
  BH_Event() {};
  ~BH_Event() {};
  UInt_t event_num;
  Monitor monitor0;
  Monitor monitor1;
  Input input;
  Detector detector0;
  Detector detector1;

  Monitor* GetMonitor(Int_t i){
		if(i == 1) return &monitor1;
		else return &monitor0; }
protected:
  ClassDef(BH_Event,2);
};

  ClassImp(BH_Event)
  ClassImp(Input)
  ClassImp(VDC)
  ClassImp(Monitor)
  ClassImp(Detector)
  ClassImp(Paddle)

//Main function used to read ttree and plot histograms.
void Histo(void){
  //Declare data to be transformed to ttree
  //For each BH_Event, there are 1 event_num,1 or 2 Monitors, 1 Input, and 1 or 2 Detectors
  //For each Detector, there are 1 VDC, and possible multiple Paddles
  //open the .dat file needed
  
  string filename;

  //Type in name of the root file
  cout<< "Input datafile name: (Format: BH.####.root)" << endl;
  getline(cin,filename);
  cout << endl;

  TFile * histofile = TFile::Open(filename.c_str());
  /*
  while(datafile.fail()){
    cout << "File not found: " << filename << endl
         <<"Please enter new name: ";
    getline(cin, filename);
    //cin.ignore(1, '\n');
    datafile.open(filename.c_str());
  }
  */
  cerr<<"BH data file " << filename << " is successfully opened." << endl;

  // Find tree in the TFile
  TTree* tree = new TTree();
  tree = (TTree * ) histofile->Get("tree");
  
  
  //SetBranchAddress to branch Event and Monitor
  BH_Event* event = new BH_Event();
  //vector<Monitor>* moni0 = new vector<Monitor>;
  //vector<Monitor>* moni1 = new vector<Monitor>;
  
  //  BH_Event* event = 0;
  Monitor* moni0 = new Monitor();
  Monitor* moni1 = new Monitor();
  moni0->ClearMonitor();
  moni1->ClearMonitor();

  //TBranch* bmoni0 = 0;
  //TBranch* bmoni1 = 0;

    tree -> SetBranchAddress("Event_Branch",&event);
    cerr << "Event_Branch is set." << endl;
  //  bmoni0 = tree -> GetBranch("monitor0");
  //  bmoni1 = tree -> GetBranch("monitor1");
  //cerr << "monitors branch is set, with address: " << &moni0 << endl;

    //    bmoni0 -> SetAddress(&moni0);
    //    bmoni1 -> SetAddress(&moni1);
  
  //Int_t NumEntry = tree -> GetEntries();
  Int_t NumEntry = tree -> GetEntries();
  cerr << "Data entries " << NumEntry << endl;
  Float_t energy_max = 30.;

  //Temporary storage of charge, energy data of monitor in a single event
  //Int_t mcharge[2];
  Float_t mE[2];

  //Canvas
  TCanvas * c_energy = new TCanvas("c_energy", "Energy Diff Distribution", 800, 400);
	c_energy->ToggleEventStatus();
	c_energy->Iconify();
 
  //Histograms
  TH1F * h_energy = new TH1F("h_energy", "Energy difference",
				120, -energy_max, energy_max);
  h_energy->SetXTitle("Energy Difference Delta(MeV)");
  h_energy->SetYTitle("Counts");

  //Read data from TTree and fill in the histograms
  Int_t j = 0;
  Int_t i = 0;

      cerr << "mE0 = " << mE[0] << endl;
      cerr << "mE1 = " << mE[1] << endl;
      cerr << "moni0E = " << moni0->energy_m << endl;
      cerr << "moni1E = " << moni1->energy_m << endl;

  //  for (Int_t i = 0; i < NumEntry; i++)
  
  for (; i < NumEntry ; ) {
    //      delete moni0; moni0 = new Monitor;
    //      delete moni1; moni1 = new Monitor;

    //    moni0->ClearMonitor();
    //    moni1->ClearMonitor();
        tree -> GetEntry(i);
    //          bmoni0 -> GetEvent(i);
    //          bmoni1 -> GetEvent(j);
      //cerr << "Reading the " << i << "th entry." << endl;
      //cerr << "Size of monitor at this event is " << event->monitors.size() << endl;
      //      mcharge[0] = 0; mcharge[1] = 0;
      mE[0] = 0.0; mE[1] = 0.0;
      /*
	  if(moni0->enumber == moni1->enumber) 
	    {
	      mE[0] = moni0->energy_m;
	      mE[1] = moni1->energy_m;
	    }
	  else if(moni0->enumber > moni1->enumber) 
	    do{
	        j++;
		bmoni1 -> GetEntry(j);
	    }while(moni0->enumber != moni1->enumber);
	  else 
	    do{
	        i++;
		bmoni0 -> GetEntry(i);
	    }while(moni0->enumber != moni1->enumber);
      */
      mE[0] = event->GetMonitor(0)->GetEnergy();
      mE[1] = event->GetMonitor(1)->GetEnergy();
      
      //      cerr << "Processing " << i+1 << "th event." << endl;
      //     cerr << "Ep = " << mE[0] << endl;
      //     cerr << "Ee = " << mE[1] << endl;
	  i++;
	  j++;
      //cerr << << "Charge in monitor is " << monit0->at(0).charge << endl;
      
      /*  Int_t j = 0;
      for (vector<Monitor>::iterator it = moni0->begin() ; it != moni0->end(); ++it)
    {
      mcharge[j] = (*it)->charge;
      mE[j] = (*it)->energy_m;
      j++;
    }
      
      */

      /*
      for (unsigned j = 0 ; j < monitor->size() ; j++)
    {
      mcharge[j] = monitor->at(j).charge;
      mE[j] = monitor->at(j).energy_m;
      cerr << "Size of monitor at this event is " << monitor->size() << endl;
      cerr << "Charge in monitor is " << mcharge[j] << endl;
      cerr << "Energy is " << mE[j] << endl;
      j++;
    }
      */
      
      //if((mcharge[0] != 0) && (mcharge[1] != 0))
	h_energy->Fill((mE[1]-mE[0]));
   }

  c_energy -> cd(); h_energy -> Draw();
  
  }
