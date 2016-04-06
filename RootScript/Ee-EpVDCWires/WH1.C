#define WH1_cxx
#include "WH.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <cstdio>
#include <iostream>

class WH1: public WH{
public:
  //void PlotVDC();
  //void Loop();
  //WH1(){}
  //WH1(TTree* tree):WH(tree) { cerr << "WH1 Tree Constructed." << endl;}

  //void WH1::Loop()
  void Plot()
{
   if (fChain == 0) return;

   //Canvases
   TCanvas* CENum = new TCanvas("CENum", "Event Number", 800, 400);
   TCanvas* CD0KE = new TCanvas("CD0KE", "Detector0 VDC WirePlane KE", 1600, 1200);
   TCanvas* CD0X = new TCanvas("CD0X", "Detector0 VDC WirePlane X",1600, 1200);
   TCanvas* CD0Y = new TCanvas("CD0Y", "Detector0 VDC WirePlane Y",1600, 1200);
   TCanvas* CD0ToF = new TCanvas("CD0ToF", "Detector0 VDC WirePlane ToF",1600, 1200);
   TCanvas* CD0Charge = new TCanvas("CD0Charge", "Detector0 VDC WirePlane Charge",1600, 1200);
   TCanvas* CD0Particle = new TCanvas("CD0Particle", "Detector0 VDC WirePlane Particle Type",1600, 1200);

   CD0KE->Divide(2,2);
   CD0X->Divide(2,2);
   CD0Y->Divide(2,2);
   CD0ToF->Divide(2,2);
   CD0Charge->Divide(2,2);
   CD0Particle->Divide(2,2);

   TCanvas* CD1KE = new TCanvas("CD1KE", "Detector1 VDC WirePlane KE", 1600, 1200);
   TCanvas* CD1X = new TCanvas("CD1X", "Detector1 VDC WirePlane X",1600, 1200);
   TCanvas* CD1Y = new TCanvas("CD1Y", "Detector1 VDC WirePlane Y",1600, 1200);
   TCanvas* CD1ToF = new TCanvas("CD1ToF", "Detector1 VDC WirePlane ToF",1600, 1200);
   TCanvas* CD1Charge = new TCanvas("CD1Charge", "Detector1 VDC WirePlane Charge",1600, 1200);
   TCanvas* CD1Particle = new TCanvas("CD1Particle", "Detector1 VDC WirePlane Particle Type",1600, 1200);
   
   CD1KE->Divide(2,2);
   CD1X->Divide(2,2);
   CD1Y->Divide(2,2);
   CD1ToF->Divide(2,2);
   CD1Charge->Divide(2,2);
   CD1Particle->Divide(2,2);   

   //Histograms
   TH1D* HENum = new TH1D("HENum", "Event Number", 100, 0, 100000);
   HENum->SetBit(TH1::kCanRebin);
   
   TH1F* H0KE[4];
   TH1F* H1KE[4];
   TH1F* H0X[4]; 
   TH1F* H1X[4];
   TH1F* H0Y[4]; 
   TH1F* H1Y[4];
   TH1F* H0ToF[4]; 
   TH1F* H1ToF[4];
   TH1F* H0Charge[4]; 
   TH1F* H1Charge[4];
   TH1F* H0Particle[4];
   TH1F* H1Particle[4];
   
   for (int i = 0; i< 4; i++){
     H0KE[i] = new TH1F(Form("H0KE WirePlane%d", i) , Form("Detector0 VDC WirePlane%d KE", i), 600, 0,60.0);
     H0X[i] = new TH1F(Form("H0X WirePlane%d", i) , Form("Detector0 VDC WirePlane%d X", i), 400, -50., 50.);
     H0X[i]->SetBit(TH1::kCanRebin);
     H0Y[i] = new TH1F(Form("H0Y WirePlane%d", i) , Form("Detector0 VDC WirePlane%d Y", i), 100, -25., 25.);
     H0Y[i]->SetBit(TH1::kCanRebin);
     H0ToF[i] = new TH1F(Form("H0ToF WirePlane%d", i), Form("Detector0 VDC WirePlane%d ToF", i), 100, 0., 10.);
     H0Charge[i] = new TH1F(Form("H0Charge WirePlane%d", i), Form("Detector0 VDC WirePlane%d Charge", i), 4, -2., 2.);
     H0Particle[i] = new TH1F(Form("H0Particle WirePlane%d", i), Form("Detector0 VDC WirePlane%d Particle",i), 50, -25., 25.);
     H0Particle[i]->SetBit(TH1::kCanRebin);

     H1KE[i] = new TH1F(Form("H1KE WirePlane%d", i) , Form("Detector1 VDC WirePlane%d KE", i), 600, 0,60.0);
     H1X[i] = new TH1F(Form("H1X WirePlane%d", i) , Form("Detector1 VDC WirePlane%d X", i), 400, -50., 50.);
     H1X[i]->SetBit(TH1::kCanRebin);
     H1Y[i] = new TH1F(Form("H1Y WirePlane%d", i) , Form("Detector1 VDC WirePlane%d Y", i), 100, -25., 25.);
     H1Y[i]->SetBit(TH1::kCanRebin);
     H1ToF[i] = new TH1F(Form("H1ToF WirePlane%d", i), Form("Detector1 VDC WirePlane%d ToF", i), 100, 0., 10.);
     H1Charge[i] = new TH1F(Form("H1Charge WirePlane%d", i), Form("Detector1 VDC WirePlane%d Charge", i), 4, -2., 2.);
     H1Particle[i] = new TH1F(Form("H1Particle WirePlane%d", i), Form("Detector1 VDC WirePlane%d Particle",i), 50, -25., 25.);
     H1Particle[i]->SetBit(TH1::kCanRebin);
   }


   // Original code
   Long64_t nentries = fChain->GetEntriesFast();

   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;

      //Fill CENum:
      HENum->Fill(ENum);

      /*
      //Find out how many data in each wireplane on each detector
      //Still all data with different wireplane is mixed
      int Index0[4]={0,0,0,0};     
      for(int n = 0; n < D0_W_; n++){
	for(int m = 0; m < 4; m++){
	  if(D0_W_WirePlane[n] == m){
	    Index0[m]++;
	  }}}

      // Do the same for D1:
      int Index1[4]={0,0,0,0};
      double TEMPKE1[4]={0.,0.,0.,0.};     
      for(int n = 0; n < D1_W_; n++){
	for(int m = 0; m < 4; m++){
	  if(D1_W_WirePlane[n] == m){
	    Index1[m]++;
	    }}}*/

      // Use Index information above and "AND" condition cuts to fill the histograms
      if( D0_W_>0 && D1_W_>0 && D0_P_>0 && D1_P_> 0){
	// Fill all data with "AND" condition in D0
	for (int i = 0; i < D0_W_; i++){
	  H0KE[i]->Fill(D0_W_KE[i]);
	  H0X[i]->Fill(D0_W_X[i]);
	  H0Y[i]->Fill(D0_W_Y[i]);
	  H0ToF[i]->Fill(D0_W_ToF[i]);
	  H0Charge[i]->Fill(D0_W_Charge[i]);
	  H0Particle[i]->Fill(D0_W_Particle[i]);}
	// Same for D1
	for(int i= 0; i < D1_W_; i++){
	  H1KE[i]->Fill(D1_W_KE[i]);
	  H1X[i]->Fill(D1_W_X[i]);
	  H1Y[i]->Fill(D1_W_Y[i]);
	  H1ToF[i]->Fill(D1_W_ToF[i]);
	  H1Charge[i]->Fill(D1_W_Charge[i]);
	  H1Particle[i]->Fill(D1_W_Particle[i]);}
      }
   }


   //Plot those histograms: these histograms shows the highest KE in each wire plane
   // and associate X,Y,ToF,Charge,Particle
   // Also note, only "AND" condition data are drawn.
   CENum->cd(); HENum->Draw();
   for(int i = 0; i < 4; i++){
     CD0KE->cd(i+1);    H0KE[i]->Draw();
     CD0X->cd(i+1);     H0X[i]->Draw();
     CD0Y->cd(i+1);     H0Y[i]->Draw();
     CD0ToF->cd(i+1);   H0ToF[i]->Draw();
     CD0Charge->cd(i+1);H0Charge[i]->Draw();
     CD0Particle->cd(i+1);H0Particle[i]->Draw();

     CD1KE->cd(i+1);    H1KE[i]->Draw();
     CD1X->cd(i+1);     H1X[i]->Draw();
     CD1Y->cd(i+1);     H1Y[i]->Draw();
     CD1ToF->cd(i+1);   H1ToF[i]->Draw();
     CD1Charge->cd(i+1);H1Charge[i]->Draw(); 
     CD1Particle->cd(i+1);H1Particle[i]->Draw();
   }
}




private:
};
