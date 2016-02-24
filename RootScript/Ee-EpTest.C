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
    Int_t Charge;

    Float_t Energy;
    Float_t X;
    Float_t Y;
    Float_t Theta;
    Float_t Phi;

  // Get functions are used for reading data from root tree.
  inline Float_t GetTheta(){return Theta;}
  inline Float_t GetPhi(){return Phi;}
  inline Float_t GetEnergy(){return Energy;}
  inline void ClearMonitor(){
    Charge = 0; Energy = 0.; X = 0.; Y = 0.; Theta = 0.; Phi = 0.;
  }
protected:
  ClassDef(Monitor,2);
  };

class Input: public TObject{
public:
  Input() {};
  ~Input() {};
  Float_t Energy;
  Float_t Delta;
  Float_t X;
  Float_t Y;
  Float_t Theta;
  Float_t Phi;
  
protected:
  ClassDef(Input,2);
};

class VDC: public TObject{
public:
  VDC() {};
  ~VDC() {};
    Float_t X;
    Float_t Y;
    Float_t Theta;
    Float_t Phi;
  Float_t E0u;
  Float_t E0v;
  Float_t E1u;
  Float_t E1v;
    Float_t GetX(){return X;}
protected:
  ClassDef(VDC,2);
};

class Paddle: public TObject{
public:
  Paddle() {};
  ~Paddle() {};
  Int_t PNum;
  Float_t Edep;
  Float_t Light;
  Float_t Time;

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
  vector<VDC> V;
  vector<Paddle> P;

  vector<VDC>* GetVDC(){return &V;}
protected:
  ClassDef(Detector,2);
};

//The number of monitors are varying (1 to 2), since one of e+ or e- or both is recorded in a single event

class BH_Event: public TObject{
public:
  BH_Event() {};
  ~BH_Event() {};
  UInt_t ENum;
  Monitor M0;
  Monitor M1;
  Input I;
  Detector D0;
  Detector D1;

  Monitor* GetMonitor(Int_t i){
		if(i == 1) return &M1;
		else return &M0; }

  Detector* GetDetector(Int_t i){
		if(i == 1) return &D1;
		else return &D0; }
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
  tree = (TTree * ) histofile->Get("T");
  
  //SetBranchAddress to branch Event
  BH_Event* event = new BH_Event();
 
  Monitor* moni0 = new Monitor();
  Monitor* moni1 = new Monitor();
  moni0->ClearMonitor();
  moni1->ClearMonitor();

  tree -> SetBranchAddress("B",&event);

  Int_t NumEntry = tree -> GetEntries();
  cerr << "Data entries " << NumEntry << endl;
  Float_t energy_max = 30.;

  //Canvas

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

      //cerr << "Size of detector0 vdc is " << event->D0.vdc.size() << "." << endl;
      //cerr << "Size of detector1 vdc is " << event->detector1.vdc.size() << "." << endl;
      //cerr << "Size of detector0 paddles is " << event->D0.paddles.size() << "." << endl;
      
      if(event->D0.V.size() == 1 && event->D1.V.size() == 1 )
	{
	  //cerr << "Detector0 vdc x = " << event->D0.vdc.at(0).x_f  << 
	  //  ". Detector1 vdc x = " <<  event->detector1.vdc.at(0).x_f<< "." << endl;
	  
	  //cerr << "Size of detector0 vdc is " << event->D0.vdc.size() << "." << endl;
	  //cerr << "Size of detector1 vdc is " << event->detector1.vdc.size() << "." << endl;
	  X0 = event->D0.V.at(0).X;
	  X1 = event->D1.V.at(0).X;
	}      

     if(event -> D0.P.size() == 2)
	{
	  //cerr << "Hit paddle 0 number = " << event->D0.paddles.at(0).paddle_num << endl;
	  // cerr << " Hit paddle 1 number = " << event->D0.paddles.at(1).paddle_num<< "." << endl;
	  PNumDiff0 = abs(event->D0.P.at(0).PNum - event->D0.P.at(1).PNum);
	}
      
      // Print out total energy and event number from root file 
      //Apply following cut: both e+ and e- with energy larger than 0.1 MeV, 
      //and sum of e+ e- energy between 58.00 and 60.00 MeV will flag kTRUE.
      if(mE[0] > 0.1 && mE[1] > 0.1)
	{
	  // Only evnets of total e- and e+ energy in monitor volume ranges from 58.00 to 59.20 will be analyzed
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
	      if(event -> D0.P.size() == 2)
		h_PaddleCheck->Fill(PNumDiff0);
	      
	      //X filter preventing meaningless events get into histogram
	      if(BadE1 != 0.0 && abs(X0) >= 0.0001 && abs(X1) >= 0.0001)
	      //if(BadE1 != 0.0)
		{

		  //cerr << "mE[0] = " << mE[0] << " , " << "mE[1] = " << mE[1] << " ." << endl;
		  //cerr << "X0 = " << X0 << " , " << "X1 = " << X1 << " ." << endl;
		  h_BadETheta1->Fill(theta,BadE1);
		  //h_BadEX0->Fill(X0,mE[0]);
		  //h_BadEX1->Fill(X1,mE[1]);
		}
	      //X filter preventing meaningless events get into gistogram
	      if (abs(X0) >= 0.0001 && abs(X1) >= 0.0001)
		{
		  if(event -> D0.P.size() == 2 && PNumDiff0 != 1) ;
		  else{
		    h_BadEX0->Fill(X0,mE[0]);
		    h_BadEX1->Fill(X1,mE[1]);
		  }
		}
	    }
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

  //Graphes and histograms in abnormal analyze
  c_BadETheta1 -> cd (); h_BadETheta1 -> Draw();
  c_BadEX0 -> cd (); h_BadEX0 -> Draw();
  c_BadEX1 -> cd (); h_BadEX1 -> Draw();
  c_PaddleCheck -> cd(); h_PaddleCheck -> Draw();
  }
