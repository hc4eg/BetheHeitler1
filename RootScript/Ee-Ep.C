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
#include "TGraph.h"
#include "TH1.h"
#include "TH2.h"
#include <cmath>

using namespace std;

//This root script is used for:
//1.Read Ep, Ee from monitors branch of .root file to make Ee - Ep 1D histogram.
//2.Read Ep, Ee and calculate Theta(opening angle) from monitors branch of .root file to make Ee Ep vs. Theta 2D histogram.

#define PI 3.1415926
Float_t MradToDeg(Float_t mrad){return mrad*180/(1000*PI);}

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

  // Get functions are used for reading data from root tree.
  inline Float_t GetTheta(){return theta_m;}
  inline Float_t GetPhi(){return phi_m;}
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

  //Get name of the root file
  cout<< "Input datafile name: (Format: ****.root)" << endl;
  getline(cin,filename);
  cout << endl;

  TFile * histofile = TFile::Open(filename.c_str());
  //Uncomment below for reread name from standard in if there is no name match in 1st type.
  /*
  while(datafile.fail()){
    cout << "File not found: " << filename << endl
         <<"Please enter new name: ";
    getline(cin, filename);
    //cin.ignore(1, '\n');
    datafile.open(filename.c_str());
  }
  
  cerr<<"BH data file " << filename << " is successfully opened." << endl;*/

  // Find tree in the TFile
  TTree* tree = new TTree();
  tree = (TTree * ) histofile->Get("tree");
  
  //SetBranchAddress to branch Event
  BH_Event* event = new BH_Event();
 
  Monitor* moni0 = new Monitor();
  Monitor* moni1 = new Monitor();
  moni0->ClearMonitor();
  moni1->ClearMonitor();

  tree -> SetBranchAddress("Event_Branch",&event);

  Int_t NumEntry = tree -> GetEntries();
  cerr << "Data entries " << NumEntry << endl;
  Float_t energy_max = 30.;

  //Canvas
   TCanvas * c_energy = new TCanvas("c_energy", "Energy Diff (Ee-Ep)/2 Distribution", 800, 400);
	c_energy->ToggleEventStatus();
	c_energy->Iconify();

   TCanvas * c_totenergy = new TCanvas("c_totenergy", "Energy Sum (Ee+Ep) Distribution", 800, 400);
	c_totenergy->ToggleEventStatus();
	c_totenergy->Iconify();
     
   TCanvas * c_theta = new TCanvas("c_theta", "Openning angle Ditribution", 800, 400);
	c_theta->ToggleEventStatus();
	c_theta->Iconify();

	/*
   TCanvas * c_thetae = new TCanvas("c_thetae", "Theta Electron Ditribution", 800, 400);
	c_theta->ToggleEventStatus();
	c_theta->Iconify();

   TCanvas * c_thetap = new TCanvas("c_thetap", "Theta Positron Ditribution", 800, 400);
	c_theta->ToggleEventStatus();
	c_theta->Iconify();
	*/
   TCanvas * c_etheta = new TCanvas("c_etheta", "EnergyDifference  (Ee-Ep)/2 vs Theta Ditribution", 800, 400);
	c_etheta->ToggleEventStatus();
	c_etheta->Iconify();

	TCanvas * c_EpDel = new TCanvas("c_EpDel", "Asymmetry Epsilon vs Delta from electron",800,400);
	c_EpDel->ToggleEventStatus();
	c_EpDel->Iconify();
	
	TCanvas * c_N = new TCanvas("c_N", "Energy difference (Ee-Ep)/2 counts",800,400);
	c_EpDel->ToggleEventStatus();
	c_EpDel->Iconify();
	
  //Histograms
  TH1F * h_energy = new TH1F("h_energy", "Energy difference (Ee-Ep)/2",
				600, -30.0, 30.0);
  h_energy->SetXTitle("Energy Difference Delta(MeV)");
  h_energy->SetYTitle("Counts");

  TH1F * h_N = new TH1F("h_N", "Energy difference (Ee-Ep)/2 count",
				120, -30.0, 30.0);
  h_energy->SetXTitle("Energy Difference Delta(MeV)");
  h_energy->SetYTitle("Counts");

  TH1F * h_totenergy = new TH1F("h_totenergy", "Ee+Ep",
				200, 58, 59.5);
  h_totenergy->SetXTitle("Total Energy(MeV)");
  h_totenergy->SetYTitle("Counts");

  TH1F * h_theta = new TH1F("h_theta", "Opening angle",
				300, 0, 40);
  h_totenergy->SetXTitle("Opening Angle(Degree)");
  h_totenergy->SetYTitle("Counts");

  TH1F * h_thetae = new TH1F("h_thetae", "Thetae",
				300, -45, 45);
  h_totenergy->SetXTitle("Electron Angle(Degree)");
  h_totenergy->SetYTitle("Counts");

  TH1F * h_thetap = new TH1F("h_thetap", "Thetap",
				300, -45, 45);
  h_totenergy->SetXTitle("Positron Angle(Degree)");
  h_totenergy->SetYTitle("Counts");

  TH2F * h_etheta = new TH2F("h_etheta", "Delta vs Theta", 100, 0, 40, 600, -energy_max, energy_max);
  h_etheta->SetXTitle("Opening Angle(Degree)");
  h_etheta->SetYTitle("Energy Difference Delta(MeV)");

  //Temporary storage of charge, energy data of monitor in a single event
  //Int_t mcharge[2];
  // mE[0] positron engergy, mE[1] electron energy.
  Float_t mE[2];
  // theta: opening angle, thetae: electron theta, thetap: positron theta, phie: eletron phi, phip: eletron phi
  Float_t theta = 0., thetae = 0., thetap = 0., phie = 0., phip = 0.;
  // AvEnergy: average energy for pairs with designate energy cut. NAv: counts of pairs for counting average energy
  Float_t AvEnergy = 0, NAv = 0;
  // Ne[]: array to store electron counts, N[]: array to store N counts(therefore to compute asymmetry epsilon)
  // Del[]: array to store delta (therefore to draw epsilon vs delta TGraph)
  Float_t Ne[120], N[120], Del[120];
  for(int i = 0; i < 120; i++){ Ne[i] = 0.; N[i] = 0.; Del[i] = 0.;}


  //Read data from TTree and fill in the histograms
  Int_t j = 0;
  Int_t i = 0;
  //  for (Int_t i = 0; i < NumEntry; i++)  
  for (; i < NumEntry ; ) {
    // Turn on flag to count valid e+ e- pairs
    Bool_t flag = kFALSE;
    tree -> GetEntry(i);

      mE[0] = 0.0; mE[1] = 0.0;

      theta = 0., thetae = 0., thetap = 0., phie = 0., phip = 0.;
      
      thetae = event->GetMonitor(1)->GetTheta();
      thetap = event->GetMonitor(0)->GetTheta();
      phie = event->GetMonitor(1)->GetPhi();
      phip = event->GetMonitor(0)->GetPhi();

      // Following code computes openging angle theta from thetae, thetap, phie, phip
      // /1000 converts mrad to rad, *180/PI converts rad to degree
      theta = 1 + tan(thetae/1000)*tan(thetap/1000) + tan(phie)*tan(phip/1000);
      theta /= sqrt(1 + tan(thetae/1000)*tan(thetae/1000) + tan(phie/1000)*tan(phie/1000));
      theta /= sqrt(1 + tan(thetap/1000)*tan(thetap/1000) + tan(phip/1000)*tan(phip/1000));      
      theta = acos(theta);
      theta = theta*180/PI;

      // Read energy of e+ mE[0], and energy of e- me[1]
      mE[0] = event->GetMonitor(0)->GetEnergy();
      mE[1] = event->GetMonitor(1)->GetEnergy();
      // Print out total energy and event number from root file 
      cerr << "Ee- = " << mE[1] << ". " << "Ee+ = " << mE[0] << ". " << "Etot = "<< mE[1]+mE[0] << ". " << "Event number = " << i << endl;
      //Apply following cut: both e+ and e- with energy larger than 0.1 MeV, 
      //and sum of e+ e- energy between 58.00 and 60.00 MeV will flag kTRUE.
      if(mE[0] > 0.1 && mE[1] > 0.1)
	{
	  for (Int_t n = -60; n < 60; n++)
	    {
	      if (((Float_t)(n)< (mE[1]-mE[0])) && ((mE[1]-mE[0]) < (Float_t)((n+1))))
		if (((mE[1]+mE[0]) > 58.00) && ((mE[1]+mE[0]) < 60.0))
		  {
		    N[n+60]++;
		    AvEnergy += (mE[0]+mE[1]);
		    NAv++;
		    flag = kTRUE; //Only Ee>0.1MeV, Ep>0.1MeV, total energy from 58.00 to 60.5MeV fill histogram.
		    // cerr << "Number of effective event " << NAv << endl;
		  }
	    }
	}
	  i++;
	  j++;

	  // Filling histograms with cutted data
	  if(flag == kTRUE)
	    {
	      h_totenergy->Fill((mE[1]+mE[0]));
	      h_energy->Fill((mE[1]-mE[0])/2);
	      h_theta->Fill(theta);
	      h_thetae->Fill(thetae);
	      h_thetap->Fill(thetap);
	      h_etheta->Fill(theta,((mE[1]-mE[0])/2));
	    }
   }

  AvEnergy /= NAv;
  cerr << "Average energy = " << AvEnergy << " Effective event entry "<< NAv << endl;
  //  for(Int_t n = 0; n< 120 ; n++) cerr << "N["<< (n-60.)/2. << "] = " << N[n] << endl;

  //  for (Int_t n = 0; n < 60; n++)
  for (Int_t n = 0; n <120; n++)
  {      
      //N[n] = (N[60+n]-N[60-n])/(N[60+n]+N[60-n]);
      //     cerr << "Asymmetry = " << N[n] << endl;
      //if(abs(N[n]) > 1.0) N[n] = 0.0;
      //Del[n] = (-29.5) + n;
      Del[n] = (-29.75)+n/2.0;
      //   cerr << "Ne 30+Del = " << Ne[60+n] << "Ne 30-Del" << Ne[60-n];
      //cerr << "Epsilon = " << N[n] << " ." << "Del = " << Del[n] << endl;
    }
  /*
  for (Int_t n = 0; n < 120; n++)
    {
      Int Nmax = N[n];
      for (Int_t i = 0; i < Nmax ; i++)
	{
	  h_N->Fill((-29.5)+n);
	}
    }
  */
  
  TGraph * g_EpDel = new TGraph(120, Del, N);
  c_EpDel->cd(); g_EpDel->Draw();
  //  c_N -> cd();h_N -> Draw();
  c_theta -> cd(); h_theta -> Draw();
  //    c_thetae -> cd(); h_thetae -> Draw();
  //    c_thetap -> cd(); h_thetap -> Draw();

  c_totenergy -> cd (); h_totenergy -> Draw();
  c_energy -> cd(); h_energy -> Draw();
  c_etheta -> cd(); h_etheta -> Draw();  
  }
