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
    TCanvas* CM1Energy = new TCanvas("CM1Energy", "M1.Energy",800,400);
    TCanvas* CM1Charge = new TCanvas("CM1Charge", "M1.Charge",800,400);
    TCanvas* CM1X = new TCanvas("CM1X", "M1.X",800,400);
    TCanvas* CM1Y = new TCanvas("CM1Y", "M1.Y",800,400);
    TCanvas* CM1Theta = new TCanvas("CM1Theta", "M1.Theta",800,400);
    TCanvas* CM1Phi = new TCanvas("CM1Phi", "M1.Phi",800,400);

    TCanvas* CM0Energy = new TCanvas("CM0Energy", "M0.Energy",800,400);
    TCanvas* CM0Charge = new TCanvas("CM0Charge", "M0.Charge",800,400);
    TCanvas* CM0X = new TCanvas("CM0X", "M0.X",800,400);
    TCanvas* CM0Y = new TCanvas("CM0Y", "M0.Y",800,400);
    TCanvas* CM0Theta = new TCanvas("CM0Theta", "M0.Theta",800,400);
    TCanvas* CM0Phi = new TCanvas("CM0Phi", "M0.Phi",800,400);

    TCanvas* CM0 = new TCanvas ("CM0", "Monitor0 data", 1600, 1200);

    CENum->cd();                T->Draw("B.ENum",Cut);
    CM1Charge->cd();            T->Draw("B.M1.Charge",Cut);
    CM1Energy->cd();            T->Draw("B.M1.Energy",Cut);
    CM1X->cd();                 T->Draw("B.M1.X",Cut);
    CM1Y->cd();                 T->Draw("B.M1.Y",Cut);
    CM1Theta->cd();             T->Draw("B.M1.Theta",Cut);
    CM1Phi->cd();               T->Draw("B.M1.Phi",Cut);

    CM0Charge->cd();            T->Draw("B.M0.Charge",Cut);
    CM0Energy->cd();            T->Draw("B.M0.Energy",Cut);
    CM0X->cd();                 T->Draw("B.M0.X",Cut);
    CM0Y->cd();                 T->Draw("B.M0.Y",Cut);
    CM0Theta->cd();             T->Draw("B.M0.Theta",Cut);
    CM0Phi->cd();               T->Draw("B.M0.Phi",Cut);
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
    TCanvas* CD1VX = new TCanvas("CD1VX", "D1.V.X",800,400);
    TCanvas* CD1VY = new TCanvas("CD1VY", "D1.V.Y",800,400);
    TCanvas* CD1VTheta = new TCanvas("CD1VTheta", "D1.V.Theta",800,400);
    TCanvas* CD1VPhi = new TCanvas("CD1VPhi", "D1.V.Phi",800,400);



    CENum->cd();                 T->Draw("B.ENum",Cut);
    CD1VX->cd();                 T->Draw("B.D1.V.X",Cut);
    CD1VY->cd();                 T->Draw("B.D1.V.Y",Cut);
    CD1VTheta->cd();             T->Draw("B.D1.V.Theta",Cut);
    CD1VPhi->cd();               T->Draw("B.D1.V.Phi",Cut);



    TCanvas* CD0VX = new TCanvas("CD0VX", "D0.V.X",800,400);
    TCanvas* CD0VY = new TCanvas("CD0VY", "D0.V.Y",800,400);
    TCanvas* CD0VTheta = new TCanvas("CD0VTheta", "D0.V.Theta",800,400);
    TCanvas* CD0VPhi = new TCanvas("CD0VPhi", "D0.V.Phi",800,400);

    CD0VX->cd();                 T->Draw("B.D0.V.X",Cut);
    CD0VY->cd();                 T->Draw("B.D0.V.Y",Cut);
    CD0VTheta->cd();             T->Draw("B.D0.V.Theta",Cut);
    CD0VPhi->cd();               T->Draw("B.D0.V.Phi",Cut);
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

  TString Cut("B.M1.Charge !=0 && B.D0.V.E0u > 0.0002 && B.D0.V.E0v > 0.0002 && B.D0.V.E1u > 0.0002 && B.D0.V.E1v > 0.0002 && B.D0.P.Edep > 0.02 && B.M0.Charge !=0 && B.D1.V.E0u > 0.0002 && B.D1.V.E0v > 0.0002 && B.D1.V.E1u > 0.0002 && B.D1.V.E1v > 0.0002 && B.D1.P.Edep > 0.02");

    TCanvas* CENum = new TCanvas("CEnum", "ENum",800,400);
    TCanvas* CD1PPNum = new TCanvas("CD1PPNum", "D1.P.PNum",800,400);
    TCanvas* CD1PEdep = new TCanvas("CD1PEdep", "D1.P.Edep",800,400);
    TCanvas* CD1PLight = new TCanvas("CD1PLight", "D1.P.Light",800,400);
    TCanvas* CD1PTime = new TCanvas("CD1PTime", "D1.P.Time",800,400);

    CENum->cd();                  T->Draw("B.ENum",Cut);
    CD1PPNum->cd();               T->Draw("B.D1.P.PNum",Cut);
    CD1PEdep->cd();               T->Draw("B.D1.P.Edep",Cut);
    CD1PLight->cd();              T->Draw("B.D1.P.Light",Cut);
    CD1PTime->cd();               T->Draw("B.D1.P.Time",Cut);


    TCanvas* CD0PPNum = new TCanvas("CD0PPNum", "D0.P.PNum",800,400);
    TCanvas* CD0PEdep = new TCanvas("CD0PEdep", "D0.P.Edep",800,400);
    TCanvas* CD0PLight = new TCanvas("CD0PLight", "D0.P.Light",800,400);
    TCanvas* CD0PTime = new TCanvas("CD0PTime", "D0.P.Time",800,400);


    CD0PPNum->cd();               T->Draw("B.D0.P.PNum",Cut);
    CD0PEdep->cd();               T->Draw("B.D0.P.Edep",Cut);
    CD0PLight->cd();              T->Draw("B.D0.P.Light",Cut);
    CD0PTime->cd();               T->Draw("B.D0.P.Time",Cut);
}
