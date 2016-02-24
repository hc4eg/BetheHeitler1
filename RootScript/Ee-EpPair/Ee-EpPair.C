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
// For the new version of the root files of the simulation. 
// (1. Name of data are abbreviated. 2. Change the way monitor data stored: leave e+, e- data to be 0 if not hit)
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

  inline Float_t GetEnergy(){return Energy;}
  inline Float_t GetDelta(){return Delta;}
  inline Float_t GetX(){return X;}
  inline Float_t GetY(){return Y;}
  inline Float_t GetTheta(){return Theta;}
  inline Float_t GetPhi(){return Phi;}
  
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
  Input I0;
  Input I1;
  Detector D0;
  Detector D1;

  Monitor* GetMonitor(Int_t i){
    if(i == 1) return &M1;
    else return &M0; }
  Input* GetInput(Int_t i){
    if(i == 1) return &I1;
    else return &I0;}
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
void HistoPair(void){
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
  /*
   TCanvas * c_energy = new TCanvas("c_energy", "Energy Diff (Ee-Ep) Distribution", 800, 400);
	c_energy->ToggleEventStatus();
	c_energy->Iconify();

   TCanvas * c_totenergy = new TCanvas("c_totenergy", "Energy Sum (Ee+Ep) Distribution", 800, 400);
	c_totenergy->ToggleEventStatus();
	c_totenergy->Iconify();
     
   TCanvas * c_theta = new TCanvas("c_theta", "Openning angle Ditribution", 800, 400);
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
  */


  //Canvas for Both Input and Monitor, draw correspondding histograms in the same canvas pad
  // Energy diffrence and energy sum
	TCanvas * c_energy = new TCanvas ("c_energy", "Energy Diff and Sum Distribution", 1800, 1000);
	c_energy->ToggleEventStatus();
	c_energy->Iconify();
	c_energy->Divide(2,1);

	// Opening angle ditribution and opening vs energy difference distribution
	TCanvas * c_theta = new TCanvas("c_theta", "Openning angle and Delta vs Theta Ditribution", 1800, 1000);
	c_theta->ToggleEventStatus();
	c_theta->Iconify();
	c_theta->Divide(2,1);
	
	// Asymmetry and asymmetry near 10 degree.
	TCanvas * c_EpDel = new TCanvas("c_EpDel", "Asymmetry Epsilon vs Delta from electron (And one near 10degree)",1800,1000);
	c_EpDel->ToggleEventStatus();
	c_EpDel->Iconify();
	c_EpDel->Divide(2,1);
	
  //Histograms
	TH1F * h_energy = new TH1F("h_energy", "Energy difference (Ee-Ep)",
				600, -60.0, 60.0);
	h_energy->SetXTitle("Energy Difference Delta(MeV)");
	h_energy->SetYTitle("Counts");

	TH1F * h_totenergy = new TH1F("h_totenergy", "Ee+Ep",
				200, 58, 59.5);
	h_totenergy->SetXTitle("Total Energy(MeV)");
	h_totenergy->SetYTitle("Counts");

	TH1F * h_theta = new TH1F("h_theta", "Opening angle",
				300, 0, 40);
	h_theta->SetXTitle("Opening Angle(Degree)");
	h_theta->SetYTitle("Counts");

	TH2F * h_etheta = new TH2F("h_etheta", "Delta vs Theta", 100, 0, 40, 600, -energy_max, energy_max);
	h_etheta->SetXTitle("Opening Angle(Degree)");
	h_etheta->SetYTitle("Energy Difference Delta(MeV)");
	h_etheta->SetOption("COLZ");

	//Histograms for Input. All Input data are set to red line while plotting histogram
	TH1F * h_Ienergy = new TH1F("h_Ienergy", "Energy difference (Ee-Ep)",
				600, -60.0, 60.0);
	h_Ienergy->SetLineColor(kRed);

	TH1F * h_Itotenergy = new TH1F("h_Itotenergy", "Ee+Ep",
				200, 58, 59.5);
	h_Itotenergy->SetLineColor(kRed);

	TH1F * h_Itheta = new TH1F("h_Itheta", "Opening angle",
				300, 0, 40);
	h_Itheta->SetLineColor(kRed);

	TH2F * h_Ietheta = new TH2F("h_Ietheta", "Delta vs Theta", 100, 0, 40, 600, -energy_max, energy_max);
	h_Ietheta->SetOption("COLZ");

  //Temporary storage of charge, energy data of monitor in a single event

  // ME[0] positron engergy, ME[1] electron energy.
  Float_t ME[2],IE[2];
  // theta: opening angle, thetae: electron theta, thetap: positron theta, phie: eletron phi, phip: eletron phi
  Float_t theta[2] = {0.,0.}, thetae[2] = {0.,0.}, thetap[2] = {0.,0.}, phie[2] = {0.,0.}, phip[2] = {0.,0.};

  //Below are method to compute (Epsilon vs Delta) by counting number of Ee-Ep
  //DelInt represents the interval of Delta (or (Ee-Ep)/2)
  //NPt is the number of Point in TGraph
  //Make sure 60/DelInt is an even number!
  Float_t DelInt = 2.0;
  Int_t NPt = (Int_t)(60/DelInt);

  // Ne[]: array to store electron counts, N[]: array to store N counts(therefore to compute asymmetry epsilon)
  // Del[]: array to store delta (therefore to draw epsilon vs delta TGraph)
  Float_t Ne[2][(Int_t)(NPt/2)],Np[2][(Int_t)(NPt/2)], N[2][NPt], Del[2][NPt];
  for(int i = 0; i < 2; i++)
    for(int j = 0; j < NPt; j++)
      { N[i][j] = 0.; Del[i][j] = 0.;}
  for (int i = 0.; i < 2; i++)
    for(int j = 0; j < (Int_t)(NPt/2); j++)
      { Ne[i][j] = 0.; Np[i][j] = 0.;}


  //Ne10 ... Version are only counting points where near angle = 10 degrees.
  Float_t Ne10[2][(Int_t)(NPt/2)],Np10[2][(Int_t)(NPt/2)], N10[2][NPt], Del10[2][NPt];
  for(int i = 0; i < 2; i++)
    for(int j = 0; j < NPt; j++)
      { N10[i][j] = 0.; Del10[i][j] = 0.;}
  for(int i = 0; i < 2; i++)
    for(int j = 0; j < (Int_t)(NPt/2); j++)
      { Ne10[i][j] = 0.; Np10[i][j] = 0.;}

  //Read data from TTree and fill in the histograms
  Int_t i = 0;
  for (; i < NumEntry ; ) {
    // Turn on flag to kTRUE only when counting valid e+ e- pairs
    Bool_t flag[2] = {kFALSE, kFALSE};
    tree -> GetEntry(i);

    ME[0] = 0.0; ME[1] = 0.0;
    IE[0] = 0.0; IE[1] = 0.0;
    for (int n = 0; n < 2; n++){ 
      theta[n] = 0., thetae[n] = 0., thetap[n] = 0., phie[n] = 0., phip[n] = 0.;}
   
    // Note these thetas from root files are using mrads as unit.
    thetae[1] = event->GetMonitor(1)->GetTheta();
    thetap[1] = event->GetMonitor(0)->GetTheta();
    phie[1] = event->GetMonitor(1)->GetPhi();
    phip[1] = event->GetMonitor(0)->GetPhi();

    // Do the same for Input
    thetae[0] = event->GetInput(1)->GetTheta();
    thetap[0] = event->GetInput(0)->GetTheta();
    phie[0] = event->GetInput(1)->GetPhi();
    phip[0] = event->GetInput(0)->GetPhi();     

    // theta[0] = abs(thetae[0]-thetap[0]);
    //theta[0] = theta[0]*180/(1000*PI);
    for(int n = 0; n < 2; n++)
      {
	thetae[n]=thetae[n]/1000;
	thetap[n]=thetap[n]/1000;
	phie[n]=phie[n]/1000;
	phip[n]=phip[n]/1000;
	Double_t SP = cos(thetae[n])*cos(thetap[n])+sin(thetae[n])*sin(thetap[n])*cos(phie[n]-phip[n]);
	theta[n] = acos(SP);
	theta[n] = theta[n]*180.0/PI;
      }

      // Following code computes openging angle theta from thetae, thetap, phie, phip
      // /1000 converts mrad to rad, *180/PI converts rad to degree, the final unit is degree
    /*
    for (int n = 1; n < 2; n++){
      theta[n] = 1 + tan(thetae[n]/1000)*tan(thetap[n]/1000) + tan(phie[n])*tan(phip[n]/1000);
      theta[n] /= sqrt(1 + tan(thetae[n]/1000)*tan(thetae[n]/1000) + tan(phie[n]/1000)*tan(phie[n]/1000));
      theta[n] /= sqrt(1 + tan(thetap[n]/1000)*tan(thetap[n]/1000) + tan(phip[n]/1000)*tan(phip[n]/1000));      
      theta[n] = acos(theta[n]);
      theta[n] = theta[n]*180/PI;
      }*/
    // Read energy of e+ ME[0], and energy of e- me[1]
    ME[0] = event->GetMonitor(0)->GetEnergy();
    ME[1] = event->GetMonitor(1)->GetEnergy();
    // Same for Input
    IE[0] = event->GetInput(0)->GetEnergy();
    IE[1] = event->GetInput(1)->GetEnergy();
    // Print out total energy and event number from root file 
    //cerr << "Ee- = " << ME[1] << ". " << "Ee+ = " << ME[0] << ". " << "Etot = "<< ME[1]+ME[0] << ". " << "Event number = " << i << endl;
  
    //Apply following cut: both e+ and e- with energy larger than 0.1 MeV, 
    //and sum of e+ e- energy between 58.00 and 60.00 MeV will flag kTRUE.
    if(ME[0] > 0.1 && ME[1] > 0.1)
      {
	if (((ME[1]+ME[0]) > 58.00) && ((ME[1]+ME[0]) < 60.0))
	  {
	    for(Int_t n = 0; n < NPt; n++)
	      {
		if (((Float_t)(n*2.0*DelInt-60.0)< (ME[1]-ME[0])) && ((ME[1]-ME[0]) < (Float_t)((n+1)*2.0*DelInt-60.0)))
		  {
		    // Counting N[n]
		    N[1][n]++;
		    // Couting N near 10 degree (points with opening angle from 9.00 degree to 11.00 degree)
		    if (theta[1] <= 11.00 && theta[1] > 9.00) N10[1][n]++;
		  }
	      }	      
	    flag[1] = kTRUE;
	  }
      }
    if(IE[0] > 0.1 && IE[1] > 0.1)
      {
	if (((IE[1]+IE[0]) > 58.00) && ((IE[1]+IE[0]) < 60.0))
	  {
	    for(Int_t n = 0; n < NPt; n++)
	      {
		if (((Float_t)(n*2.0*DelInt-60.0)< (IE[1]-IE[0])) && ((IE[1]-IE[0]) < (Float_t)((n+1)*2.0*DelInt-60.0)))
		  {
		    // Counting N[n]
		    N[0][n]++;
		    // Couting N near 10 degree (points with opening angle from 9.00 degree to 11.00 degree)
		    if (theta[0] <= 11.00 && theta[0] > 9.00) N10[0][n]++;
		  }
	      }	      
	    flag[0] = kTRUE;
	  }
      }

    i++;
    // Filling histograms with cutted data
    if(flag[1] == kTRUE)
      {
	h_totenergy->Fill((ME[1]+ME[0]));
	h_energy->Fill((ME[1]-ME[0]));
	h_theta->Fill(theta[1]);
	h_etheta->Fill(theta[1],((ME[1]-ME[0])/2));
      }
    if(flag[0] == kTRUE)
      {
	h_Itotenergy->Fill((IE[1]+IE[0]));
	h_Ienergy->Fill((IE[1]-IE[0]));
	h_Itheta->Fill(theta[0]);
	h_Ietheta->Fill(theta[0],((IE[1]-IE[0])/2));
      }   	  
  }
  
  for(Int_t i = 0; i < 2; i++)
    for(Int_t j = 0; j < NPt ; j++) 
      cerr << "N[" << i <<"]["<< (j*2*DelInt-60.)/2. << "] = " << N[i][j] << endl;
  //for(Int_t n = 0; n< 120 ; n++) cerr << "N10["<< (n-60.)/2. << "] = " << N10[n] << endl;
  
  //Below are procedures to compute Epsilon for Epsilon (Asymmetry) vs Delta ((Ee-Ep)/2)
  for (Int_t i = 0; i < 2; i++)
    for (Int_t j = 0; j < (Int_t)(NPt/2); j++)
      {
	Ne[i][j] = N[i][j+(Int_t)(NPt/2)];
	Np[i][j] = N[i][j];
      }
  
  for (Int_t i = 0; i < 2; i++)
    for (Int_t j = 0; j < (Int_t)(NPt/2); j++)
      {   
	//Need to use temporary array to store N[n]      
	N[i][j] = (Np[i][j]-Ne[i][(NPt/2)-j-1])/(Np[i][j]+Ne[i][(NPt/2)-j-1]);
	if (Np[i][j]==Ne[i][(NPt/2)-j-1]) N[i][j]=0.0;
	N[i][(NPt/2)+j] = (Ne[i][j]-Np[i][(NPt/2)-j-1])/(Ne[i][j]+Np[i][(NPt/2)-j-1]);
	if (Ne[i][j]==Np[i][(NPt/2)-j-1]) N[i][(NPt/2)+j]=0.0;
      }

  for (Int_t i = 0; i < 2; i++)
    for (Int_t j = 0; j < NPt; j++)
      {
	//  cerr << "Asymmetry[" << n <<"] = " << N[n] << endl;
	Del[i][j] = -30.0+ (DelInt/2) + DelInt*j;
      }
  
  // Ploting Epsilon vs Delta graph near opening angle 10 degree
  for (Int_t i = 0; i < 2; i++ )
    for (Int_t j = 0; j < (Int_t)(NPt/2); j++)
      {
	Ne10[i][j] = N10[i][j+(Int_t)(NPt/2)];
	Np10[i][j] = N10[i][j];
      }

  for (Int_t i = 0; i < 2; i++)
    for(Int_t j = 0; j< NPt ; j++) 
      cerr << "N10[" << i <<"]["<< (j*2*DelInt-60.)/2. << "] = " << N10[i][j] << endl;

  for (Int_t i = 0; i < 2; i++)
    for (Int_t j = 0; j < (Int_t)(NPt/2); j++)
      {   
	//Need to use temporary array to store N[n]      
	N10[i][j] = (Np10[i][j]-Ne10[i][(NPt/2)-j-1])/(Np10[i][j]+Ne10[i][(NPt/2)-j-1]);
	if (Np10[i][j]==Ne10[i][(NPt/2)-j-1]) N10[i][j]=0.0;
	N10[i][(NPt/2)+j] = (Ne10[i][j]-Np10[i][(NPt/2)-j-1])/(Ne10[i][j]+Np10[i][(NPt/2)-j-1]);
	if (Ne10[i][j]==Np10[i][(NPt/2)-j-1]) N10[i][(NPt/2)+j]=0.0;
      }

  for (Int_t i = 0; i < 2; i++)
    for (Int_t j = 0; j < NPt; j++)
      {
	Del10[i][j] = -30.0+ (DelInt/2) + DelInt*j;
      }


  //
  TGraph * g_EpDel = new TGraph(NPt, Del[1], N[1]);
  TGraph * g_IEpDel = new TGraph(NPt, Del[0], N[0]);    
  g_EpDel->SetLineColor(kRed);    
  c_EpDel->cd(1); g_EpDel->Draw(); g_IEpDel->Draw("same");
  
  TGraph * g_EpDel10 = new TGraph(NPt, Del10[1], N10[1]);
  TGraph * g_IEpDel10 = new TGraph(NPt, Del10[0],N10[0]);
  g_IEpDel10->SetLineColor(kRed);
  c_EpDel->cd(2); g_EpDel10->Draw();g_IEpDel10->Draw("same");
 
  //
  c_theta -> cd(1); h_Itheta -> Draw();h_theta -> Draw("same");
  c_theta -> cd(2); h_etheta -> Draw(); h_Ietheta -> Draw("same");
  
  c_energy -> cd(2); h_totenergy -> Draw(); h_Itotenergy -> Draw("same");
  c_energy -> cd(1); h_energy -> Draw(); h_Ienergy -> Draw("same");
  }
