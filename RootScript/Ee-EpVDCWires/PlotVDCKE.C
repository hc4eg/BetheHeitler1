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
#include <cmath>

using namespace std;

//Usage of Plot.C:
// Ananlyzing simulation root files, extract data by using AND cut (Only take events with all 4 detector hit)
// Start to use at Pair.0152.root
// New Histograms: VDC: KE , ToF, Charge. Hodoscope KE.

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
  
  inline Int_t GetPNum(){ return PNum;}
  inline Float_t GetLight(){ return Light;}

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
  vector<Paddle>* GetPad(){ return &P;}
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

void PlotInput(){
  //Open the root file
  string filename;
  cout<< "Input datafile name: (Format: ****.root)" << endl;
  getline(cin,filename);
  cout << endl;
  TFile * histofile = TFile::Open(filename.c_str());

  string OoN;
  cout << "Old or New file? (O/N)" << endl;
  getline(cin,OoN);
  cout << endl;

  TTree* T = new TTree();
  T = (TTree * ) histofile->Get("T");

  //TString Cut("B.M1.Charge !=0 && B.D0.V.E0u > 0.0002 && B.D0.V.E0v > 0.0002 && B.D0.V.E1u > 0.0002 && B.D0.V.E1v > 0.0002 && B.D0.P.Edep > 0.02 && B.M0.Charge !=0 && B.D1.V.E0u > 0.0002 && B.D1.V.E0v > 0.0002 && B.D1.V.E1u > 0.0002 && B.D1.V.E1v > 0.0002 && B.D1.P.Edep > 0.02");
  TString Cut(
    "B.M1.Charge !=0 && "  
    "B.D0.V.E0u > 0.0002 && "  
    "B.D0.V.E0v > 0.0002 && "  
    "B.D0.V.E1u > 0.0002 && "  
    "B.D0.V.E1v > 0.0002 && "  
    "B.D0.P.Edep > 0.02 && "  
    "B.M0.Charge !=0 && "  
    "B.D1.V.E0u > 0.0002 && "  
    "B.D1.V.E0v > 0.0002 && "  
    "B.D1.V.E1u > 0.0002 && "  
    "B.D1.V.E1v > 0.0002 && "  
    "B.D1.P.Edep > 0.02"
  );
  if (OoN == "N"){

    cerr << "OoN is " << OoN << endl;
    TCanvas* CENum = new TCanvas("CEnum", "ENum",800,400);
    CENum->cd();           T->Draw("B.ENum",Cut);

    TCanvas* CI0 = new TCanvas ("CI0", "Input0 data", 1600, 1500);
    CI0->Divide(3,2);

    CI0->cd(1);            T->Draw("B.I0.Energy",Cut);
    CI0->cd(2);            T->Draw("B.I0.Delta",Cut);
    CI0->cd(3);            T->Draw("B.I0.X",Cut);
    CI0->cd(4);            T->Draw("B.I0.Y",Cut);
    CI0->cd(5);            T->Draw("B.I0.Theta",Cut);
    CI0->cd(6);            T->Draw("B.I0.Phi",Cut);

    TCanvas* CI1 = new TCanvas ("CI1", "Input1 data", 1600, 1500);
    CI1->Divide(3,2);

    CI1->cd(1);            T->Draw("B.I1.Energy",Cut);
    CI1->cd(2);            T->Draw("B.I1.Delta",Cut);
    CI1->cd(3);            T->Draw("B.I1.X",Cut);
    CI1->cd(4);            T->Draw("B.I1.Y",Cut);
    CI1->cd(5);            T->Draw("B.I1.Theta",Cut);
    CI1->cd(6);            T->Draw("B.I1.Phi",Cut);
  }
  else{
    TCanvas* CENum = new TCanvas("CEnum", "ENum",800,400);
    CENum->cd();           T->Draw("B.ENum",Cut);
 
    TCanvas* CI = new TCanvas ("CI", "Input data", 1600, 1500);
    CI->Divide(3,2);

    CI->cd(1);            T->Draw("B.I.Energy",Cut);
    CI->cd(2);            T->Draw("B.I.Delta",Cut);
    CI->cd(3);            T->Draw("B.I.X",Cut);
    CI->cd(4);            T->Draw("B.I.Y",Cut);
    CI->cd(5);            T->Draw("B.I.Theta",Cut);
    CI->cd(6);            T->Draw("B.I.Phi",Cut);
  }    
}

void PlotMonitor(){
  //Open the root file
  string filename;
  cout<< "Input datafile name: (Format: ****.root)" << endl;
  getline(cin,filename);
  cout << endl;
  TFile * histofile = TFile::Open(filename.c_str());

  TTree* T = new TTree();
  T = (TTree * ) histofile->Get("T");

  TString Cut("B.M1.Charge !=0 && B.D0.V.E0u > 0.0002 && B.D0.V.E0v > 0.0002 && B.D0.V.E1u > 0.0002 && B.D0.V.E1v > 0.0002 && B.D0.P.Edep > 0.02 && B.M0.Charge !=0 && B.D1.V.E0u > 0.0002 && B.D1.V.E0v > 0.0002 && B.D1.V.E1u > 0.0002 && B.D1.V.E1v > 0.0002 && B.D1.P.Edep > 0.02");
    
    TCanvas* CENum = new TCanvas("CEnum", "ENum",800,400);

    TCanvas* CM0 = new TCanvas ("CM0", "Monitor0 data", 1600, 1500);
    CM0->Divide(3,2);

    TCanvas* CM1 = new TCanvas ("CM1", "Monitor1 data", 1600, 1500);
    CM1->Divide(3,2);

    CENum->cd();           T->Draw("B.ENum",Cut);
    
    CM1->cd(1);            T->Draw("B.M1.Charge",Cut);
    CM1->cd(2);            T->Draw("B.M1.Energy",Cut);
    CM1->cd(3);            T->Draw("B.M1.X",Cut);
    CM1->cd(4);            T->Draw("B.M1.Y",Cut);
    CM1->cd(5);            T->Draw("B.M1.Theta",Cut);
    CM1->cd(6);            T->Draw("B.M1.Phi",Cut);

    CM0->cd(1);            T->Draw("B.M0.Charge",Cut);
    CM0->cd(2);            T->Draw("B.M0.Energy",Cut);
    CM0->cd(3);            T->Draw("B.M0.X",Cut);
    CM0->cd(4);            T->Draw("B.M0.Y",Cut);
    CM0->cd(5);            T->Draw("B.M0.Theta",Cut);
    CM0->cd(6);            T->Draw("B.M0.Phi",Cut);
}

void PlotVDC(){
  //Open the root file
  string filename;
  cout<< "Input datafile name: (Format: ****.root)" << endl;
  getline(cin,filename);
  cout << endl;
  TFile * histofile = TFile::Open(filename.c_str());

  TTree* T = new TTree();
  T = (TTree * ) histofile->Get("T");

  TString Cut("B.M1.Charge !=0 && B.D0.V.E0u > 0.0002 && B.D0.V.E0v > 0.0002 && B.D0.V.E1u > 0.0002 && B.D0.V.E1v > 0.0002 && B.D0.P.Edep > 0.02 && B.M0.Charge !=0 && B.D1.V.E0u > 0.0002 && B.D1.V.E0v > 0.0002 && B.D1.V.E1u > 0.0002 && B.D1.V.E1v > 0.0002 && B.D1.P.Edep > 0.02");


    TCanvas* CENum = new TCanvas("CEnum", "ENum",800,400);

    TCanvas* CD1V = new TCanvas ("D1VDC", "Detector1 VDC data", 1600, 1000);
    CD1V->Divide(2,2);
    TCanvas* CD1VKE = new TCanvas ("D1VDC_KE", "Detector1 VDC KE ToF Charge data", 1600,1000);
    CD1VKE->Divide(2,2);

    CENum->cd();             T->Draw("B.ENum",Cut);
    
    CD1V->cd(1);             T->Draw("B.D1.V.X",Cut);
    CD1V->cd(2);             T->Draw("B.D1.V.Y",Cut);
    CD1V->cd(3);             T->Draw("B.D1.V.Theta",Cut);
    CD1V->cd(4);             T->Draw("B.D1.V.Phi",Cut);

    CD1VKE->cd(1);           T->Draw("B.D1.V.KE", Cut);
    CD1VKE->cd(2);           T->Draw("B.D1.V.ToF", Cut);
    CD1VKE->cd(3);           T->Draw("B.D1.V.Charge",Cut);

    TCanvas* CD0V = new TCanvas ("D0VDC", "Detector0 VDC data", 1600, 1000);
    CD0V->Divide(2,2);    
    TCanvas* CD0VKE = new TCanvas ("D0VDC_KE", "Detector0 VDC: KE ToF Charge data", 1600,1000);
    CD0VKE->Divide(2,2);

    CD0V->cd(1);             T->Draw("B.D0.V.X",Cut);
    CD0V->cd(2);             T->Draw("B.D0.V.Y",Cut);
    CD0V->cd(3);             T->Draw("B.D0.V.Theta",Cut);
    CD0V->cd(4);             T->Draw("B.D0.V.Phi",Cut);

    CD0VKE->cd(1);           T->Draw("B.D0.V.KE", Cut);
    CD0VKE->cd(2);           T->Draw("B.D0.V.ToF", Cut);
    CD0VKE->cd(3);           T->Draw("B.D0.V.Charge",Cut);
}

void PlotPaddle(){
    //Open the root file
    string filename;
    cout<< "Input datafile name: (Format: ****.root)" << endl;
    getline(cin,filename);
    cout << endl;
    TFile * histofile = TFile::Open(filename.c_str());
  
    TTree* T = new TTree();
    T = (TTree * ) histofile->Get("T");
  
    //SetBranchAddress to branch Event
    BH_Event* E = new BH_Event();
    T -> SetBranchAddress("B",&E);

    Int_t NumEntry = T -> GetEntries();
    cerr << "Data entries " << NumEntry << endl;

    //TString Cut("B.M1.Charge !=0 && B.D0.V.E0u > 0.0002 && B.D0.V.E0v > 0.0002 && B.D0.V.E1u > 0.0002 && B.D0.V.E1v > 0.0002 && B.D0.P.Edep > 0.02 && B.M0.Charge !=0 && B.D1.V.E0u > 0.0002 && B.D1.V.E0v > 0.0002 && B.D1.V.E1u > 0.0002 && B.D1.V.E1v > 0.0002 && B.D1.P.Edep > 0.02");
    TString Cut(
		//"B.M1.Charge !=0 && "  
		"B.D0.V.E0u > 0.0002 && "  
		"B.D0.V.E0v > 0.0002 && "  
		"B.D0.V.E1u > 0.0002 && "  
		"B.D0.V.E1v > 0.0002 && "  
		//"B.D0.P.Edep > 0.02 && "
		"B.D0.P.Light > 0.02 &&"
		//"B.M0.Charge !=0 && "  
		"B.D1.V.E0u > 0.0002 && "  
		"B.D1.V.E0v > 0.0002 && "  
		"B.D1.V.E1u > 0.0002 && "  
		"B.D1.V.E1v > 0.0002 && "  
		//"B.D1.P.Edep > 0.02"
		"B.D1.P.Light > 0.02"
		);

    TCanvas* CENum = new TCanvas("CEnum", "ENum",800,400);

    TCanvas* CD1P = new TCanvas("D1P","Detector1 Paddle data",1600,1200);
    CD1P->Divide(2,2);

    TH1F* H_D1PNum = new TH1F("H_D1PNum", "D1 Paddle Num distribution.", 29,0.0,29.0);
    H_D1PNum->SetLineColor(kRed);
    TH1F* H_D0PNum = new TH1F("H_D0PNum", "D0 Paddle Num distribution.", 29,0.0,29.0);
    H_D0PNum->SetLineColor(kRed);

    for(Int_t i = 0; i < NumEntry; i ++){
      // This is extremely important!
      T ->GetEntry(i);
      cerr << "Entry Num = " << i << "." << endl;
      cerr << "D0 VDC size = " << E->GetDetector(0)-> GetVDC()->size() << endl;
      cerr << "D0P size = " << E-> GetDetector(0)-> GetPad()->size() << endl;
      
      for(UInt_t j = 0; j <  E-> GetDetector(0)-> GetPad()->size();j++){
	
	if(E-> GetDetector(0) -> GetPad() -> size() > 0 &&  
	   E-> GetDetector(1) -> GetPad() -> size() > 0){
	  H_D0PNum -> Fill(E->GetDetector(0) -> GetPad() -> at(j).PNum);
	  cerr << "D0 PNum = " << E->GetDetector(0) -> GetPad() -> at(j).PNum << endl;
	}}
     
      cerr << "D1P size = " << E->GetDetector(1)-> GetPad()->size() << endl;
      for(UInt_t j = 0; j < E->GetDetector(1)-> GetPad()->size() ;j++){
	if(E-> GetDetector(0) -> GetPad() -> size() > 0 &&  
	   E-> GetDetector(1) -> GetPad() -> size() > 0){
	  H_D1PNum -> Fill(E->GetDetector(1) -> GetPad() -> at(j).PNum);
	  cerr << "D1 PNum = " << E->GetDetector(1) -> GetPad() -> at(j).PNum<< endl;
	}}      
    }
    CENum->cd();              T->Draw("B.ENum",Cut);

    CD1P->cd(1);              T->Draw("B.D1.P.PNum>>hd1(28,0,28)",Cut);
    //CD1P->cd(1);              H_D1PNum->Draw();
                              T->Draw("B.D1.P.PNum",Cut,"same");
    CD1P->cd(2);              T->Draw("B.D1.P.Edep",Cut);
    CD1P->cd(3);              T->Draw("B.D1.P.Light",Cut);
    CD1P->cd(4);              T->Draw("B.D1.P.Time",Cut);

    TCanvas* CD0P = new TCanvas("D0P","Detector0 Paddle data",1600,1200);
    CD0P->Divide(2,2);

    CD0P->cd(1);              T->Draw("B.D0.P.PNum >> (28,0,28)",Cut);
    //CD0P->cd(1);              H_D0PNum->Draw();
                              T->Draw("B.D0.P.PNum",Cut,"same");
    //CD0P->cd(1);              T->Draw("B.D0.P.PNum >> (28,0,28)",Cut);
    CD0P->cd(2);              T->Draw("B.D0.P.Edep",Cut);
    CD0P->cd(3);              T->Draw("B.D0.P.Light",Cut);
    CD0P->cd(4);              T->Draw("B.D0.P.Time",Cut);

    TCanvas* CHKE = new TCanvas("HKE","Hodoscope KE",1600,600);
    CHKE->Divide(2,1);
    
    CHKE->cd(1);              T->Draw("B.HKE0",Cut);
    CHKE->cd(2);              T->Draw("B.HKE1",Cut);
}
