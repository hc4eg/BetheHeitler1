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
  Double_t KE;
  Double_t ToF;
  Double_t Charge;

  inline Double_t GetKE(){return KE;}
  inline Double_t GetCharge(){return Charge;}
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

  vector<VDC>* GetVDC(){ return &V;}
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
  vector<Double_t> HKE0;
  vector<Double_t> HKE1;

  vector<Double_t>* GetHKE(Int_t i){
    if(i == 1) return &HKE1;
    else return &HKE0; }

  Monitor* GetMonitor(Int_t i){
    if(i == 1) return &M1;
    else return &M0; }
  Input* GetInput(Int_t i){
    if(i == 1) return &I1;
    else return &I0;}
  Detector* GetDetector(Int_t i){
    if(i == 1) return &D1;
    else return &D0;}
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
void HistoKE(void){
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
	TCanvas * c_EpDel = new TCanvas("c_EpDel", "Total Asymmetry vs Delta: Input Red, Monitor Black, VDC Blue, Hodoscope Green",1800,1000);
	c_EpDel->ToggleEventStatus();
	c_EpDel->Iconify();
	c_EpDel->Divide(2,1);

	TCanvas * c_IMEpDel = new TCanvas("c_IMEpDel", "A vs Delta near some opening angle, Left Input, Right Monitor",1800,1000);
	c_IMEpDel->ToggleEventStatus();
	c_IMEpDel->Iconify();
	c_IMEpDel->Divide(2,1);	
	
  //Histograms
	TH1F * h_energy = new TH1F("h_energy", "Energy difference (Ee-Ep)",
				600, -60.0, 60.0);
	h_energy->SetXTitle("Energy Difference Delta(MeV)");
	h_energy->SetYTitle("Counts");
	h_energy->SetLineColor(kBlack);

	TH1F * h_totenergy = new TH1F("h_totenergy", "Ee+Ep",
				200, 30.0, 80.0);
	h_totenergy->SetXTitle("Total Energy(MeV)");
	h_totenergy->SetYTitle("Counts");
	h_totenergy->SetLineColor(kBlack);

	TH1F * h_theta = new TH1F("h_theta", "Opening angle",
				300, 0, 40);
	h_theta->SetXTitle("Opening Angle(Degree)");
	h_theta->SetYTitle("Counts");
	h_theta->SetLineColor(kBlack);

	TH2F * h_etheta = new TH2F("h_etheta", "Delta vs Theta", 100, 0, 40, 600, -energy_max, energy_max);
	h_etheta->SetXTitle("Opening Angle(Degree)");
	h_etheta->SetYTitle("Energy Difference Delta(MeV)");
	h_etheta->SetOption("COLZ");

	//Histograms for Input. All Input data are set to red line while plotting histogram
	TH1F * h_Ienergy = new TH1F("h_Ienergy", "Energy difference (Ee-Ep)",
				600, -60.0, 60.0);
	h_Ienergy->SetLineColor(kRed);

	TH1F * h_Itotenergy = new TH1F("h_Itotenergy", "Ee+Ep",
				200, 30.0, 80.0);
	h_Itotenergy->SetLineColor(kRed);

	TH1F * h_Itheta = new TH1F("h_Itheta", "Opening angle",
				300, 0, 40);
	h_Itheta->SetLineColor(kRed);

	TH2F * h_Ietheta = new TH2F("h_Ietheta", "Delta vs Theta", 100, 0, 40, 600, -energy_max, energy_max);
	h_Ietheta->SetOption("COLZ");

	//Histograms for VDC and Hodoscope. VDC use Blue, Hod use Green lines.
	TH1F * h_Venergy = new TH1F("h_Venergy", "VDC Energy difference (Ee-Ep)",
				600, -60.0, 60.0);
	h_Venergy->SetLineColor(kBlue);

	TH1F * h_Vtotenergy = new TH1F("h_Vtotenergy", "Ee+Ep",
				200, 30.0, 80.0);
	h_Vtotenergy->SetLineColor(kBlue);	
	//Hod
	TH1F * h_Henergy = new TH1F("h_Henergy", "Hodoscope Energy difference (Ee-Ep)",
				600, -60.0, 60.0);
	h_Henergy->SetLineColor(kGreen);

	TH1F * h_Htotenergy = new TH1F("h_Htotenergy", "Ee+Ep",
				200, 30.0, 80.0);
	h_Htotenergy->SetLineColor(kGreen);
  //Temporary storage of charge, energy data of monitor in a single event

  // ME[0] positron engergy, ME[1] electron energy. And keep track of Charge in VDC. 
	Float_t ME[2],IE[2], VE[2], HE[2]; Float_t VCharge[2] = {100, 100};
  // theta: opening angle, thetae: electron theta, thetap: positron theta, phie: eletron phi, phip: eletron phi
  Float_t theta[2] = {0.,0.}, thetae[2] = {0.,0.}, thetap[2] = {0.,0.}, phie[2] = {0.,0.}, phip[2] = {0.,0.};

  //Below are method to compute (Epsilon vs Delta) by counting number of Ee-Ep
  //DelInt represents the interval of Delta (or (Ee-Ep)/2)
  //NPt is the number of Point in TGraph
  //Make sure 60/DelInt is an even number!
  Float_t DelInt = 2.0;
  Int_t NPt = (Int_t)(60/DelInt);

  // Ne[]: array to store electron counts, N[]: array to store N counts(therefore to compute asymmetry epsilon)
  // N[0] Input data; N[1] Monitor data; N[2] VDC data; N[3] Hod data;
  // Del: array to store delta (therefore to draw epsilon vs delta TGraph)
  Float_t Ne[4][(Int_t)(NPt/2)],Np[4][(Int_t)(NPt/2)], N[4][NPt], Del[NPt];
  for(int i = 0; i < 4; i++)
    for(int j = 0; j < NPt; j++)
      { N[i][j] = 0.; Del[j] = 0.;}
  for (int i = 0.; i < 4; i++)
    for(int j = 0; j < (Int_t)(NPt/2); j++)
      { Ne[i][j] = 0.; Np[i][j] = 0.;}


  //NeOA ... Version are only counting points where near some certain opening angle (determined by OAmax)
  Float_t NeOA[2][8][(Int_t)(NPt/2)],NpOA[2][8][(Int_t)(NPt/2)], NOA[2][8][NPt], DelOA[NPt]; Float_t OAmax = 40.0;
  cerr << "Asymmetry plot near following opening angle is plot:";
  for (int n = 0; n <8; n++) cerr << " " << OAmax/8.0*(n+1) << " " ; cerr << endl; 
  for(int i = 0; i < 2; i++)
    for(int j = 0; j < 8; j++)
      for(int k = 0; k < NPt; k++)
	{ NOA[i][j][k] = 0.; DelOA[k] = 0.;}
  for(int i = 0; i < 2; i++)
    for (int j = 0; j < 2; j++)
      for(int k = 0; k < (Int_t)(NPt/2); k++)
      { NeOA[i][j][k] = 0.; NpOA[i][j][k] = 0.;}





  //Read data from TTree and fill in the histograms
  Int_t i = 0;
  for (; i < NumEntry ; ) {
    // Turn on flag to kTRUE only when counting valid e+ e- pairs
    Bool_t flag[4] = {kFALSE, kFALSE, kFALSE, kFALSE};
    tree -> GetEntry(i);

    for (int n = 0; n < 2; n++){
      ME[n] = 0.0; IE[n] = 0.0; VE[n] = 0.0; HE[n] = 0.0; VCharge[n] = 100.0;}
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





    // Read energy of e+ ME[0], and energy of e- me[1]
    ME[0] = event->GetMonitor(0)->GetEnergy();
    ME[1] = event->GetMonitor(1)->GetEnergy();
    // Same for Input
    IE[0] = event->GetInput(0)->GetEnergy();
    IE[1] = event->GetInput(1)->GetEnergy();
    // Print out total energy and event number from root file 
    //cerr << "Ee- = " << ME[1] << ". " << "Ee+ = " << ME[0] << ". " << "Etot = "<< ME[1]+ME[0] << ". " << "Event number = " << i << endl;
    if ( event->GetDetector(0)->GetVDC()->size() == 1 && event->GetDetector(1)->GetVDC()->size() == 1){
      VE[0] = event->GetDetector(0)->GetVDC()->at(0).GetKE(); VE[1] = event->GetDetector(1)->GetVDC()->at(0).GetKE();}
      //cerr << "D0 VDC KE = " << VE[0] << ". D1 VDC KE = " << VE[1] << endl; }
    if ( event->GetHKE(0)->size() == 1 && event->GetHKE(1)->size() == 1){
      HE[0] = event->GetHKE(0)->at(0); HE[1] = event->GetHKE(1)->at(0);}
    // Find out charge of particles in VDC
    if (event->GetDetector(0)->GetVDC()->size() == 1) VCharge[0] = event->GetDetector(0)->GetVDC()->at(0).GetCharge();
    if (event->GetDetector(1)->GetVDC()->size() == 1) VCharge[1] = event->GetDetector(1)->GetVDC()->at(0).GetCharge();
  




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
		    // Couting N near some Opening Angle +-1.0 degree
		    for (Int_t m = 0; m < 8; m++){
		      if (theta[1] <= (OAmax/8.0)*(m+1)+1.0 && (OAmax/8.0)*(m+1)-1.0) NOA[1][m][n]++;}
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
		    for (Int_t m = 0; m < 8; m++){
		      if (theta[0] <= (OAmax/8.0)*(m+1)+1.0 && (OAmax/8.0)*(m+1)-1.0) NOA[0][m][n]++;}		    
		  }
	      }	      
	    flag[0] = kTRUE;
	  }
      }

    //VDC kinetic energy
    if(VE[0] > 0.1 && VE[1] > 0.1 && (VCharge[0] == 1. && VCharge[1] == -1.) && abs((VE[1] - VE[0])- (IE[1] - IE[0])) < 5.0 )
      //if(VE[0] > 0.1 && VE[1] > 0.1)
      {
      //if ( VE[1]+VE[0] < 80.0){
	{
	  for(Int_t n = 0; n < NPt; n++){
	    if (((Float_t)(n*2.0*DelInt-60.0)< (VE[1]-VE[0])) && ((VE[1]-VE[0]) < (Float_t)((n+1)*2.0*DelInt-60.0))) N[2][n]++;}
	}
	flag[2] = kTRUE;
      }
    //Hodoscope kinetic energy
    if(HE[0] > 0.1 && HE[1] > 0.1 && (VCharge[0] == 1. && VCharge[1] == -1.) ){
      //if ( HE[1]+HE[0] < 80.0){
      {      
	for(Int_t n = 0; n < NPt; n++){
	  if (((Float_t)(n*2.0*DelInt-60.0)< (HE[1]-HE[0])) && ((HE[1]-HE[0]) < (Float_t)((n+1)*2.0*DelInt-60.0))) N[3][n]++;}
      }
      flag[3] = kTRUE;
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
    if(flag[2] == kTRUE){
	h_Vtotenergy->Fill((VE[1]+VE[0]));
	h_Venergy->Fill((VE[1]-VE[0]));}
    if(flag[3] == kTRUE){
	h_Htotenergy->Fill((HE[1]+HE[0]));
	h_Henergy->Fill((HE[1]-HE[0]));}    
  }
  
  
  for(Int_t i = 0; i < 4; i++)
    for(Int_t j = 0; j < NPt ; j++) 
      cerr << "N[" << i <<"]["<< (j*2*DelInt-60.)/2. << "] = " << N[i][j] << endl;
  //for(Int_t n = 0; n< 120 ; n++) cerr << "NOA["<< (n-60.)/2. << "] = " << NOA[n] << endl;


  //Below are procedures to compute Epsilon for Epsilon (Asymmetry) vs Delta ((Ee-Ep)/2)
  for (Int_t i = 0; i < 4; i++)
    for (Int_t j = 0; j < (Int_t)(NPt/2); j++)
      {
	Ne[i][j] = N[i][j+(Int_t)(NPt/2)];
	Np[i][j] = N[i][j];
      }

  for (Int_t i = 0; i < 4; i++)
    for (Int_t j = 0; j < (Int_t)(NPt/2); j++)
      {   
	//Need to use temporary array to store N[n]      
	N[i][j] = (Np[i][j]-Ne[i][(NPt/2)-j-1])/(Np[i][j]+Ne[i][(NPt/2)-j-1]);
	if (Np[i][j]==Ne[i][(NPt/2)-j-1]) N[i][j]=0.0;
	N[i][(NPt/2)+j] = (Ne[i][j]-Np[i][(NPt/2)-j-1])/(Ne[i][j]+Np[i][(NPt/2)-j-1]);
	if (Ne[i][j]==Np[i][(NPt/2)-j-1]) N[i][(NPt/2)+j]=0.0;
      }

  for (Int_t i = 0; i < NPt; i++){
    //  cerr << "Asymmetry[" << n <<"] = " << N[n] << endl;
    Del[i] = -30.0+ (DelInt/2) + DelInt*i;}

  // Ploting Epsilon vs Delta graph near opening angle 10 degree
  for (Int_t i = 0; i < 2; i++ )
    for (Int_t j = 0; j < 8; j++)
      for (Int_t k = 0; k < (Int_t)(NPt/2); k++){
	NeOA[i][j][k] = NOA[i][j][k+(Int_t)(NPt/2)];
	NpOA[i][j][k] = NOA[i][j][k];}

  /*
  for (Int_t i = 0; i < 2; i++)
    for(Int_t k = 0; k < NPt ; k++) 
      cerr << "NOA[" << i <<"]["<< (k*2*DelInt-60.)/2. << "] = " << NOA[i][4][k] << endl;
  */

  for (Int_t i = 0; i < 2; i++)
    for (Int_t j = 0; j < 8; j++)
      for (Int_t k = 0; k < (Int_t)(NPt/2); k++)
	{   
	  //Need to use temporary array to store N[n]      
	  NOA[i][j][k] = (NpOA[i][j][k]-NeOA[i][j][(NPt/2)-k-1])/(NpOA[i][j][k]+NeOA[i][j][(NPt/2)-k-1]);
	  if (NpOA[i][j][k]==NeOA[i][j][(NPt/2)-k-1]) NOA[i][j][k]=0.0;
	  NOA[i][j][(NPt/2)+k] = (NeOA[i][j][k]-NpOA[i][j][(NPt/2)-k-1])/(NeOA[i][j][k]+NpOA[i][j][(NPt/2)-k-1]);
	  if (NeOA[i][j][k]==NpOA[i][j][(NPt/2)-k-1]) NOA[i][j][(NPt/2)+k]=0.0;
	}

    for (Int_t i = 0; i < NPt; i++) DelOA[i] = -30.0+ (DelInt/2) + DelInt*i;


  //
  TGraph * g_MTEpDel = new TGraph(NPt, Del, N[1]);
  TGraph * g_ITEpDel = new TGraph(NPt, Del, N[0]);
  TGraph * g_VTEpDel = new TGraph(NPt, Del, N[2]);
  TGraph * g_HTEpDel = new TGraph(NPt, Del, N[3]);
  g_ITEpDel->SetLineColor(kRed);
  g_VTEpDel->SetLineColor(kBlue);
  g_HTEpDel->SetLineColor(kGreen);
  c_EpDel->cd(1);g_HTEpDel->Draw(); g_MTEpDel->Draw("same"); g_ITEpDel->Draw("same"); g_VTEpDel->Draw("same");

   vector<TGraph*> g_IEpDel;
   vector<TGraph*> g_MEpDel;
  for (Int_t n = 0; n < 8; n++){
    g_IEpDel.push_back(new TGraph(NPt, DelOA, NOA[0][n]));
    g_IEpDel.at(n)->SetLineColor(n+1);
    g_MEpDel.push_back(new TGraph(NPt, DelOA, NOA[1][n]));
    g_MEpDel.at(n)->SetLineColor(n+1);    
  }

  c_IMEpDel->cd(1); g_IEpDel.at(0)->Draw();
    for (Int_t i = 0; i < 8; i++) g_IEpDel.at(i)->Draw("same");
    c_IMEpDel->cd(2);g_MEpDel.at(0)->Draw();
    for (Int_t i = 0; i < 8; i++) g_MEpDel.at(i)->Draw("same");
 
  //
  c_theta -> cd(1); h_Itheta -> Draw();h_theta -> Draw("same");
  c_theta -> cd(2); h_etheta -> Draw(); h_Ietheta -> Draw("same");
  // Ee + Ep
  c_energy -> cd(2);h_Vtotenergy -> Draw(); h_Itotenergy -> Draw("same");  h_totenergy -> Draw("same");h_Htotenergy -> Draw("same");
  // Ee - Ep
  c_energy -> cd(1);h_Ienergy -> Draw(); h_energy -> Draw("same"); h_Venergy -> Draw("same"); h_Henergy -> Draw("same");
  }
