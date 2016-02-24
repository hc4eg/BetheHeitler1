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
#include <cmath>

using namespace std;

//Usage of Plot.C:
// Ananlyzing simulation root files, extract data by using AND cut (Only take events with all 4 detector hit)
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

  TString Cut("B.M1.Charge !=0 && B.D0.V.E0u > 0.0002 && B.D0.V.E0v > 0.0002 && B.D0.V.E1u > 0.0002 && B.D0.V.E1v > 0.0002 && B.D0.P.Edep > 0.02 && B.M0.Charge !=0 && B.D1.V.E0u > 0.0002 && B.D1.V.E0v > 0.0002 && B.D1.V.E1u > 0.0002 && B.D1.V.E1v > 0.0002 && B.D1.P.Edep > 0.02");
    
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

    CENum->cd();             T->Draw("B.ENum",Cut);
    
    CD1V->cd(1);             T->Draw("B.D1.V.X",Cut);
    CD1V->cd(2);             T->Draw("B.D1.V.Y",Cut);
    CD1V->cd(3);             T->Draw("B.D1.V.Theta",Cut);
    CD1V->cd(4);             T->Draw("B.D1.V.Phi",Cut);

    TCanvas* CD0V = new TCanvas ("D0VDC", "Detector0 VDC data", 1600, 1000);
    CD0V->Divide(2,2);    

    CD0V->cd(1);             T->Draw("B.D0.V.X",Cut);
    CD0V->cd(2);             T->Draw("B.D0.V.Y",Cut);
    CD0V->cd(3);             T->Draw("B.D0.V.Theta",Cut);
    CD0V->cd(4);             T->Draw("B.D0.V.Phi",Cut);
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

  TString Cut(
	      //"B.M1.Charge !=0 &&" 
	      "(B.D0.V.E0u > 0.0002 || "
	      "B.D0.V.E0v > 0.0002 || "
	      "B.D0.V.E1u > 0.0002 || "
	      "B.D0.V.E1v > 0.0002  ||"
	      "B.D0.P.Edep > 0.02) ||" 
	      //"B.M0.Charge !=0 && "
	      "(B.D1.V.E0u > 0.0002 || "
	      "B.D1.V.E0v > 0.0002 || "
	      "B.D1.V.E1u > 0.0002 || "
	      "B.D1.V.E1v > 0.0002 ||" 
	      "B.D1.P.Edep > 0.02)");

    TCanvas* CENum = new TCanvas("CEnum", "ENum",800,400);

    TCanvas* CD1P = new TCanvas("D1P","Detector1 Paddle data",1600,1200);
    CD1P->Divide(2,2);

    CENum->cd();              T->Draw("B.ENum",Cut);

    CD1P->cd(1);              //T->Draw("B.D1.P.PNum",Cut);
    //T->Draw("B.D1.P.PNum >> hd1");
    //T->Draw("B.D1.P.PNum","B.D0.P.Edep > 0.02 || B.D1.P.Edep > 0.02");
    T->Draw("B.D1.P.PNum","B.D0.P.Edep > 0.02 && B.D1.P.Edep > 0.02");
    CD1P->cd(2);              T->Draw("B.D1.P.Edep",Cut);
    CD1P->cd(3);              T->Draw("B.D1.P.Light",Cut);
    CD1P->cd(4);              T->Draw("B.D1.P.Time",Cut);

    TCanvas* CD0P = new TCanvas("D0P","Detector0 Paddle data",1600,1200);
    CD0P->Divide(2,2);

    CD0P->cd(1);              //T->Draw("B.D0.P.PNum >> (29,0,29)",Cut);
    //T->Draw("B.D0.P.PNum >> hd2");
    //T->Draw("B.D0.P.PNum","B.D0.P.Edep > 0.02 || B.D1.P.Edep > 0.02");
    T->Draw("B.D0.P.PNum","B.D0.P.Edep > 0.02 && B.D1.P.Edep > 0.02");
    CD0P->cd(2);              T->Draw("B.D0.P.Edep",Cut);
    CD0P->cd(3);              T->Draw("B.D0.P.Light",Cut);
    CD0P->cd(4);              T->Draw("B.D0.P.Time",Cut);
}
