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

  Float_t GetX(){return x_f;}
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

  vector<VDC>* GetVDC(){return &vdc;}
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

  Detector* GetDetector(Int_t i){
		if(i == 1) return &detector1;
		else return &detector0; }
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

  /*
   TCanvas * c_energy = new TCanvas("c_energy", "Energy Diff (Ee-Ep)/2 Distribution", 800, 400);
	c_energy->ToggleEventStatus();
	c_energy->Iconify();

   TCanvas * c_totenergy = new TCanvas("c_totenergy", "Energy Sum (Ee+Ep) Distribution", 800, 400);
	c_totenergy->ToggleEventStatus();
	c_totenergy->Iconify();
     
   TCanvas * c_theta = new TCanvas("c_theta", "Openning angle Ditribution", 800, 400);
	c_theta->ToggleEventStatus();
	c_theta->Iconify();

	
   TCanvas * c_thetae = new TCanvas("c_thetae", "Theta Electron Ditribution", 800, 400);
	c_theta->ToggleEventStatus();
	c_theta->Iconify();

   TCanvas * c_thetap = new TCanvas("c_thetap", "Theta Positron Ditribution", 800, 400);
	c_theta->ToggleEventStatus();
	c_theta->Iconify();
	
   TCanvas * c_etheta = new TCanvas("c_etheta", "EnergyDifference  (Ee-Ep)/2 vs Theta Ditribution", 800, 400);
	c_etheta->ToggleEventStatus();
	c_etheta->Iconify();

   TCanvas * c_EpDel = new TCanvas("c_EpDel", "Asymmetry Epsilon vs Delta from electron",800,400);
	c_EpDel->ToggleEventStatus();
	c_EpDel->Iconify();

   TCanvas * c_EpDel10 = new TCanvas("c_EpDel10", "Asymmetry Epsilon vs Delta from electron at opening angle around 10 degree",800,400);
	c_EpDel->ToggleEventStatus();
	c_EpDel->Iconify();
	
	TCanvas * c_N = new TCanvas("c_N", "Energy difference (Ee-Ep)/2 counts",800,400);
	c_EpDel->ToggleEventStatus();
	c_EpDel->Iconify();
	*/
  //Canvas to Test those abnormal data points |Ee-Ep|>20MeV and data in vdc and paddles which @size is 0

   TCanvas * c_BadETheta1 = new TCanvas("c_BadETheta1", "EnergyDifference  (Ee-Ep)/2 (from -60MeV to 20MeV) vs Opening Angle Ditribution", 800, 400);
	c_BadETheta1->ToggleEventStatus();
	c_BadETheta1->Iconify();
   
   TCanvas * c_BadEX0 = new TCanvas("c_BadEX0", "Positron Energy vs Detector0 VDC X Ditribution", 800, 400);
	c_BadEX0->ToggleEventStatus();
	c_BadEX0->Iconify();

   TCanvas * c_BadEX1 = new TCanvas("c_BadEX1", "Electron Energy vs Detector1 VDC X Ditribution", 800, 400);
	c_BadEX1->ToggleEventStatus();
	c_BadEX1->Iconify();

   TCanvas * c_PaddleCheck = new TCanvas("c_PaddleCheck", "Checking if hit paddles are adjacent", 800, 400);
	c_PaddleCheck->ToggleEventStatus();
	c_PaddleCheck->Iconify();
  //Histograms
  /*
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
  h_etheta->SetOption("COLZ");
  */

  //Canvas to Test those abnormal data points |Ee-Ep|>20MeV and data in vdc and paddles which @size is 0
  TH2F * h_BadETheta1 = new TH2F("h_BadETheta1", "Delta vs Theta", 100, 0, 40, 200, -60.0, -20);
  h_BadETheta1->SetXTitle("Opening Angle(Degree)");
  h_BadETheta1->SetYTitle("Energy Difference Delta(MeV)");
  h_BadETheta1->SetOption("COLZ");
  
  TH2F * h_BadEX0 = new TH2F("h_BadEX0", "Positron Energy vs Detector0 VDC X", 100, -100.0, 100.0, 300, 0.0, 60.0);
  h_BadEX0->SetXTitle("Detector0 VDC X(cm)");
  h_BadEX0->SetYTitle("Positron Energy (MeV)");
  h_BadEX0->SetOption("COLZ");
  
  TH2F * h_BadEX1 = new TH2F("h_BadEX1", "Electron Energy vs Detector1 VDC X", 100, -100.0, 100.0, 300, 0.0, 60.0);
  h_BadEX1->SetXTitle("Detctor1 VDC X(cm)");
  h_BadEX1->SetYTitle("Electron Energy(MeV)");
  h_BadEX1->SetOption("COLZ");
  
  TH1F * h_PaddleCheck = new TH1F("h_PaddleCheck", "PaddleCheck",
			     300, 0,  30);
  h_PaddleCheck->SetXTitle("|HitPaddleNum1 - HitPaddleNum0|");
  h_PaddleCheck->SetYTitle("Counts");
  //Temporary storage of charge, energy data of monitor in a single event
  //Int_t mcharge[2];
  // mE[0] positron engergy, mE[1] electron energy.
  Float_t mE[2];
  // theta: opening angle, thetae: electron theta, thetap: positron theta, phie: eletron phi, phip: eletron phi
  Float_t theta = 0., thetae = 0., thetap = 0., phie = 0., phip = 0.;
  // AvEnergy: average energy for pairs with designate energy cut. NAv: counts of pairs for counting average energy
  Float_t AvEnergy = 0, NAv = 0;

  //Below are method to compute (Epsilon - Delta) by counting number of Ee-Ep
  //DelInt represents the interval of Delta (or (Ee-Ep)/2)
  //NPt is the number of Pt in TGraph
  Float_t DelInt = 2.0;
  Int_t NPt = (Int_t)(60/DelInt);

  // Ne[]: array to store electron counts, N[]: array to store N counts(therefore to compute asymmetry epsilon)
  // Del[]: array to store delta (therefore to draw epsilon vs delta TGraph)
  Float_t Ne[(Int_t)(NPt/2)],Np[(Int_t)(NPt/2)], N[NPt], Del[NPt];
  for(int i = 0; i < NPt; i++){ N[i] = 0.; Del[i] = 0.;}
  for(int i = 0; i < (Int_t)(NPt/2); i++){ Ne[i] = 0.; Np[i] = 0.;}

  Float_t Ne10[(Int_t)(NPt/2)],Np10[(Int_t)(NPt/2)], N10[NPt], Del10[NPt];
  for(int i = 0; i < NPt; i++){ N10[i] = 0.; Del10[i] = 0.;}
  for(int i = 0; i < (Int_t)(NPt/2); i++){ Ne10[i] = 0.; Np10[i] = 0.;}



  //Tesing variable
  Float_t BadE1 = 0.0;
  Float_t X0 = 0.0, X1 = 0.0;
  Int_t PNumDiff0 = 0, PNumDiff1 = 0;

  //Read data from TTree and fill in the histograms
  Int_t j = 0;
  Int_t i = 0;
 
  for (; i < NumEntry ; ) {
    // Turn on flag to count valid e+ e- pairs
      Bool_t flag = kFALSE;
      tree -> GetEntry(i);

      mE[0] = 0.0; mE[1] = 0.0;
      theta = 0., thetae = 0., thetap = 0., phie = 0., phip = 0.;

      BadE1 = 0.0;
      X0 = 0.0, X1 = 0.0;
      PNumDiff0 = 0, PNumDiff1 = 0;
      
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

      //cerr << "Size of detector0 vdc is " << event->detector0.vdc.size() << "." << endl;
      //cerr << "Size of detector1 vdc is " << event->detector1.vdc.size() << "." << endl;
      //cerr << "Size of detector0 paddles is " << event->detector0.paddles.size() << "." << endl;
      
      if(event->detector0.vdc.size() == 1 && event->detector1.vdc.size() == 1 )
	{
	  //cerr << "Detector0 vdc x = " << event->detector0.vdc.at(0).x_f  << 
	  //  ". Detector1 vdc x = " <<  event->detector1.vdc.at(0).x_f<< "." << endl;
	  
	  //cerr << "Size of detector0 vdc is " << event->detector0.vdc.size() << "." << endl;
	  //cerr << "Size of detector1 vdc is " << event->detector1.vdc.size() << "." << endl;
	  X0 = event->detector0.vdc.at(0).x_f;
	  X1 = event->detector1.vdc.at(0).x_f;
	}
      

      
     if(event -> detector0.paddles.size() == 2)
	{
	  //cerr << "Hit paddle 0 number = " << event->detector0.paddles.at(0).paddle_num << endl;
	  // cerr << " Hit paddle 1 number = " << event->detector0.paddles.at(1).paddle_num<< "." << endl;
	  PNumDiff0 = abs(event->detector0.paddles.at(0).paddle_num - event->detector0.paddles.at(1).paddle_num);
	}
      
      // Print out total energy and event number from root file 
      //cerr << "Ee- = " << mE[1] << ". " << "Ee+ = " << mE[0] << ". " << "Etot = "<< mE[1]+mE[0] << ". " << "Event number = " << i << endl;
      //Apply following cut: both e+ and e- with energy larger than 0.1 MeV, 
      //and sum of e+ e- energy between 58.00 and 60.00 MeV will flag kTRUE.
      if(mE[0] > 0.1 && mE[1] > 0.1)
	{
		if (((mE[1]+mE[0]) > 58.00) && ((mE[1]+mE[0]) < 59.2))
		  {
		    for(Int_t n = 0; n < NPt; n++)
		      {
			if (((Float_t)(n*2.0*DelInt-60.0)< (mE[1]-mE[0])) && ((mE[1]-mE[0]) < (Float_t)((n+1)*2.0*DelInt-60.0)))
			  {
			    flag = kTRUE;
			    N[n]++;
			    if (theta < 11.00 && theta >= 9.00) N10[n]++;
			  }
		      }
		  }

	}

      if(mE[1]-mE[0] >= -60.0 && mE[1]-mE[0] <= -20.0)
	BadE1 = mE[1]-mE[0];

	  i++;
	  j++;

	  // Filling histograms with cutted data
	  if(flag == kTRUE)
	    {
	      /*
	      h_totenergy->Fill((mE[1]+mE[0]));
	      h_energy->Fill((mE[1]-mE[0])/2);
	      h_theta->Fill(theta);
	      h_thetae->Fill(thetae);
	      h_thetap->Fill(thetap);
	      h_etheta->Fill(theta,((mE[1]-mE[0])/2));
	      */
	      if(event -> detector0.paddles.size() == 2)
		h_PaddleCheck->Fill(PNumDiff0);
	      
	      //X filter preventing meaningless events get into histogram
	      if(BadE1 != 0.0 && abs(X0) >= 0.001 && abs(X1) >= 0.001)
	      //if(BadE1 != 0.0)
		{

		  //cerr << "mE[0] = " << mE[0] << " , " << "mE[1] = " << mE[1] << " ." << endl;
		  //cerr << "X0 = " << X0 << " , " << "X1 = " << X1 << " ." << endl;
		  h_BadETheta1->Fill(theta,BadE1);
		  //h_BadEX0->Fill(X0,mE[0]);
		  //h_BadEX1->Fill(X1,mE[1]);
		}
	      //X filter preventing meaningless events get into gistogram
	      if (abs(X0) >= 0.001 && abs(X1) >= 0.001)
		{
		  if(event -> detector0.paddles.size() == 2 && PNumDiff0 != 1) ;
		  else{
		    h_BadEX0->Fill(X0,mE[0]);
		    h_BadEX1->Fill(X1,mE[1]);
		  }
		}
	    }
   }

    for (Int_t n = 0 ; n < 30000; n++)
    {
      Float_t Count = h_BadEX0->GetBinContent(n);
      //cerr << "Bin content at bin #" << n << " is " << Count << endl;
      if (Count <= 200) h_BadEX0->SetBinContent(n,0.0);
    }
  
  AvEnergy /= NAv;
  cerr << "Average energy = " << AvEnergy << " Effective event entry "<< NAv << endl;
  for(Int_t n = 0; n< NPt ; n++) cerr << "N["<< (n*2*DelInt-60.)/2. << "] = " << N[n] << endl;
  //for(Int_t n = 0; n< 120 ; n++) cerr << "N10["<< (n-60.)/2. << "] = " << N10[n] << endl;

  for (Int_t n = 0; n < (Int_t)(NPt/2); n++)
    {
      Ne[n] = N[n+(Int_t)(NPt/2)];
      Np[n] = N[n];
    }
  
  for (Int_t n = 0; n < (Int_t)(NPt/2); n++)
  {   
    //Need to use temporary array to store N[n]      
    N[n] = (Np[n]-Ne[(NPt/2)-n-1])/(Np[n]+Ne[(NPt/2)-n-1]);
    if (Np[n]==Ne[(NPt/2)-n-1]) N[n]=0.0;
    N[(NPt/2)+n] = (Ne[n]-Np[(NPt/2)-n-1])/(Ne[n]+Np[(NPt/2)-n-1]);
    if (Ne[n]==Np[(NPt/2)-n-1]) N[(NPt/2)+n]=0.0;
  }

  for (Int_t n = 0; n < NPt; n++)
    {
      //  cerr << "Asymmetry[" << n <<"] = " << N[n] << endl;
      Del[n] = -30.0+ (DelInt/2) + DelInt*n;
    }
  
  // Ploting Epsilon vs Delta graph near opening angle 10 degree
  for (Int_t n = 0; n < (Int_t)(NPt/2); n++)
    {
      Ne10[n] = N10[n+(Int_t)(NPt/2)];
      Np10[n] = N10[n];
    }
  
  for (Int_t n = 0; n < (Int_t)(NPt/2); n++)
  {   
      //Need to use temporary array to store N[n]      
    N10[n] = (Np10[n]-Ne10[(NPt/2)-n-1])/(Np10[n]+Ne10[(NPt/2)-n-1]);
    if (Np10[n]==Ne10[(NPt/2)-n-1]) N10[n]=0.0;
    N10[(NPt/2)+n] = (Ne10[n]-Np10[(NPt/2)-n-1])/(Ne10[n]+Np10[(NPt/2)-n-1]);
    if (Ne10[n]==Np10[(NPt/2)-n-1]) N10[(NPt/2)+n]=0.0;
  }

  for (Int_t n = 0; n < NPt; n++)
    {
      //  cerr << "Asymmetry[" << n <<"] = " << N[n] << endl;
      Del10[n] = -30.0+ (DelInt/2) + DelInt*n;
    }
  

  /*
  //Graph and histograms in normal analyze
  TGraph * g_EpDel = new TGraph(NPt, Del, N);
  c_EpDel->cd(); g_EpDel->Draw();

  TGraph * g_EpDel10 = new TGraph(NPt, Del10, N10);
  c_EpDel10->cd(); g_EpDel10->Draw();
  //  c_N -> cd();h_N -> Draw();
  c_theta -> cd(); h_theta -> Draw();
  //    c_thetae -> cd(); h_thetae -> Draw();
  //    c_thetap -> cd(); h_thetap -> Draw();

  c_totenergy -> cd (); h_totenergy -> Draw();
  c_energy -> cd(); h_energy -> Draw();
  c_etheta -> cd(); h_etheta -> Draw();
  */

  //Graphes and histograms in abnormal analyze
  c_BadETheta1 -> cd (); h_BadETheta1 -> Draw();
  c_BadEX0 -> cd (); h_BadEX0 -> Draw();
  c_BadEX1 -> cd (); h_BadEX1 -> Draw();
  c_PaddleCheck -> cd(); h_PaddleCheck -> Draw();
  }
