#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include "TFile.h"
#include "TTree.h"
#include "TRandom.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "TH1.h"
#include "TH2.h"
#include <cmath>

using namespace std;

#define PI 3.14159265

// Usage of the file: Computing following situation: A gamma photon in x direction hit a rest nucleus in lab reference frame.
// Then a e+ e- pair is generated with recoiled nucleus moving in photon direction.
// The electron/positron is generated with uniform random Energy and Theta (Angle relative to x-axis in xy plane).
// Then compute the distribution of energy and theta, plot histogram.

// Plots: 1. Thetap , Thetae distribution
//        2. Opening Angle distribution
//        3. Ee, Ep distribution
//        4. Asymmetry: Total, near 10 degree, near 20 degree curve


void TpTeFiedxASym (){
  Double_t Ee, Ep;
  Double_t Thetae, Thetap, Phie, Phip;
  Double_t E,Me;
  Double_t Pp,Ppx,Ppy,Ppz,Pe, Pex,Pey,Pez;
  Double_t ThetaM, Temp;

  Double_t A,a, Tempa;
  // Genertate e+, e- pairs with fixed asymmetry A
  A = 0.2;
  a = (1-A)/(1+A);

  Int_t N = 50000000;

  E = 60.0;
  Me = 0.511;
  //Me = 10.0;

  /*
  TString fname = "TeTpSym.root"
  TFile* F = new TFile(fname,"recreate");
  TTree* T = new TTree("T","tree");
  TBranch* B = T->Branch("Thetap", &Thetap,"Thetap/D");
  */


  //Canvas
  TCanvas * c_thetape = new TCanvas("c_thetap", "Theta Positron Distribution", 800, 400);
	c_thetape->ToggleEventStatus();
	c_thetape->Iconify();

  TCanvas * c_theta = new TCanvas("c_thetae", "Theta Electron Distribution", 800, 400);
	c_theta->ToggleEventStatus();
	c_theta->Iconify();

  TCanvas * c_Ee = new TCanvas("c_Ee", "Electron Energy Distribution", 800, 400);
	c_Ee->ToggleEventStatus();
	c_Ee->Iconify();

  TCanvas * c_Ep = new TCanvas("c_Ep", "Positron Energy Distribution", 800, 400);
	c_Ep->ToggleEventStatus();
	c_Ep->Iconify();

  //Histogram
  TH1F * h_thetap = new TH1F("h_thetap", "Thetap",
				300, -90.0, 90.0);
  h_thetap->SetXTitle("Positron Angle(Degree)");
  h_thetap->SetYTitle("Counts");
  
  TH1F * h_thetae = new TH1F("h_thetae", "Thetae",
				300, -90.0, 90.0);
  h_thetae->SetXTitle("Electron Angle(Degree)");
  h_thetae->SetYTitle("Counts");
  h_thetae->SetLineColor(kRed);

  TH1F * h_theta = new TH1F("h_theta", "Theta",
				300, 0.0, 180.0);
  h_theta->SetXTitle("Opening Angle(Degree)");
  h_theta->SetYTitle("Counts");
  h_theta->SetLineColor(kRed);

  TH1F * h_Ep = new TH1F("h_Ep", "Ep",
				300, 0.0, 60.0);
  h_Ep->SetXTitle("Ep(MeV)");
  h_Ep->SetYTitle("Counts");

  TH1F * h_Ee = new TH1F("h_Ee", "Ee",
				300, 0.0, 60.0);
  h_Ee->SetXTitle("Ee(MeV)");
  h_Ee->SetYTitle("Counts");


  //Computing asymmetry from histogram
  Double_t Del = 2.0; Int_t Nm = (Int_t)(120.0/Del);
  Double_t Num[Nm], Asym[Nm], del[Nm], N10[Nm], A10[Nm], N20[Nm], A20[Nm], AConst[Nm];
  for (Int_t n = 0; n < Nm; n++){
    Num[n] = 0.0; Asym[n] = 0.0;
    N10[n] = 0.0; A10[n] = 0.0;
    N20[n] = 0.0; A20[n] = 0.0;
    AConst[n] = 0.0;
    del[n] = 0.0;
  }
  
  for (Int_t i = 0; i < N; i++)
    {      
      //Following code is to generate e+ e- pair with same distribution.
      //if( i%2 == 0){
	do{
	  Tempa = gRandom -> Uniform(0,1+a);
	  

	  // Here the range of Ee distribution can be altered as long as inside 0 to 60 MeV range.
	  if(Tempa <= 1){
	    Ee = gRandom->Uniform(0.0, 30.0);
	  }
	  else{
	    Ee = gRandom->Uniform(30.0, 60.0);
	  }
	  //Thetae = gRandom->Uniform(-PI/4, PI/4);

	  //Assigning angle theta
	  //do{
	    Temp = sqrt(((E-Ee)*(E-Ee)-Me*Me)/(Ee*Ee-Me*Me));	  
	    if( 1 > abs(Temp)){
	      ThetaM = asin(Temp);
	      Thetae = gRandom->Uniform(-ThetaM, ThetaM);
	    }
	    else
	      Thetae = gRandom->Uniform(-PI/2,PI/2);
	    
	    Pe = sqrt(Ee*Ee-Me*Me);
	    Pex = Pe*cos(Thetae);
	    Pey = Pe*sin(Thetae);
	    
	    Ep = E - Ee;
	    Ppx = sqrt((E-Ee)*(E-Ee)-(Ee*Ee-Me*Me)*sin(Thetae)*sin(Thetae)-Me*Me);
	    Ppy = -sqrt(Ee*Ee-Me*Me)*sin(Thetae);

	    Thetap = atan2(Ppy,Ppx);
	    Thetap = Thetap*180/PI;
	    Thetae = Thetae*180/PI;
	    //}while(abs(Thetae*1.0) > 20.0 || abs(Thetap*1.0) > 20.0);
	}while(  ((E-Ee)*(E-Ee)-Pey*Pey-Me*Me < 0.0) || Ee < Me); 

	for (Int_t n = 0; n < Nm; n++)
	  {
	    if ( (Ee-Ep) > (-60.0 + (Double_t)(n*Del/1.0)) && (Ee-Ep) <= (-60.0 + (Double_t)((n+1)*Del/1.0)))
	      {
		Num[n]++;
		if (abs(Thetap-Thetae)<=11.0 && abs(Thetap-Thetae)>9.0)
		  N10[n]++;
		if (abs(Thetap-Thetae)<=32.0 && abs(Thetap-Thetae)>28.0)
		  N20[n]++;
	      }
	  }
	
		 //(Ee < 15.0) || (Ee > 45.0) || Thetae > 20.0 || Thetap > 20.0 || Thetae < -20.0 || Thetap < -20.0  );
	    //     }
      
      /*
      else{
	do{

	  Tempa = gRandom -> Uniform(0,1+a);
	  
	  if(Tempa <= 1){
	    Ep = gRandom->Uniform(0.0, 30.0);
	  }
	  else{
	    Ep = gRandom->Uniform(30.0, 60.0);
	  } 
	  //Thetae = gRandom->Uniform(-PI/4, PI/4);
	  Temp = sqrt(((E-Ep)*(E-Ep)-Me*Me)/(Ep*Ep-Me*Me));
	  if(-1 < Temp && 1 > Temp){
	    ThetaM = asin(Temp);
	    Thetap = gRandom->Uniform(-ThetaM, ThetaM);
	  }
	  else
	    Thetap = gRandom->Uniform(-90,90);
	  
	  Pp = sqrt(Ep*Ep-Me*Me);
	  Ppx = Pp*cos(Thetap);
	  Ppy = Pp*sin(Thetap);
	  
	  Ee = E - Ep;
	  Pex = sqrt((E-Ep)*(E-Ep)-(Ep*Ep-Me*Me)*sin(Thetap)*sin(Thetap)-Me*Me);
	  Pey = -sqrt(Ep*Ep-Me*Me)*sin(Thetap);
	}while (Ep < Me || ((E-Ep)*(E-Ep)-Ppy*Ppy-Me*Me <0) ||(E-Ep) < Me);
	    
	    Thetae = atan2(Pey,Pex);
	    Thetae = Thetae*180/PI;
	    Thetap = Thetap*180/PI;
      }
      */

      // Condition more like simulation with OR contition (require a range of energy to go through magnetic field, one side detector hit)
      // if((abs(Thetae) < 15.0 && (Ee < 40.0) && (Ee > 20.0)) || (abs(Thetap) < 15.0 && (Ep < 40.0) && (Ep > 20.0))){
      
      // Condition more like simulation with AND condition (require both side detector hit)
      //if((abs(Thetae) < 20.0 && (Ee < 45.0) && (Ee > 15.0)) && (abs(Thetap) < 20.0 && (Ep < 45.0) && (Ep > 15.0))){	
	
      // Looser OR requirement, simulate events go through yoke of magnet
      //if(abs(Thetae) < 15.0 || abs(Thetap) < 15.0){

      // Looser AND requirement, simulate events go through yoke of magnet
      //if(abs(Thetae*1.0) < 20.0 && abs(Thetap*1.0) < 20.0){
	h_thetap->Fill(Thetap);
	h_thetae->Fill(Thetae);
	h_theta->Fill(abs(Thetap-Thetae));
	h_Ee->Fill(Ee);
	h_Ep->Fill(Ep);
	//}
    }

  //Compute asymmetry
  for (Int_t n = 0; n < Nm; n++){
    del[n] = -30.0 + Del/4.0 + (Double_t)(n*Del/2.0);
    Asym[n] = (Num[n]-Num[Nm-n-1])/(Num[n]+Num[Nm-n-1]);
    if (Num[n] == Num[Nm-n-1]) Asym[n] = 0.0;
    if( Num[n] <= (N/(Nm*5.0)) || Num[Nm-n-1] <= (N/(Nm*5.0))) Asym[n] = 0.0;
    
    A10[n] = (N10[n]-N10[Nm-n-1])/(N10[n]+N10[Nm-n-1]);
    if (N10[n] == N10[Nm-n-1]) A10[n] = 0.0;    
    if (N10[n] == 0.0 || N10[Nm-n-1] == 0.0) A10[n] = 0.0;

    A20[n] = (N20[n]-N20[Nm-n-1])/(N20[n]+N20[Nm-n-1]);
    if (N20[n] == N20[Nm-n-1]) A20[n] = 0.0;    
    if (N20[n] == 0.0 || N20[Nm-n-1] == 0.0) A20[n] = 0.0;

    if(n < Nm/2.0) AConst[n] = A;
    else AConst[n] = -A;
  }
  
  TCanvas * c_EpDel = new TCanvas("c_EpDel", "Asymmetry Epsilon vs Delta from electron",1600,900);
  c_EpDel->ToggleEventStatus();
  c_EpDel->Iconify();
  TGraph* g_EpDel = new TGraph(Nm, del, Asym);
  
  TGraph* g_EpDel10 = new TGraph(Nm, del, A10);
  g_EpDel10->SetLineColor(kRed);
  
  TGraph* g_EpDel20 = new TGraph(Nm, del, A20);
  g_EpDel20->SetLineColor(kBlue);

  TGraph* g_Const = new TGraph(Nm, del, AConst);
  g_Const->SetLineColor(kGreen);
  
  c_EpDel->cd(); g_EpDel->Draw(); g_EpDel10->Draw("same");g_EpDel20->Draw("same"); g_Const->Draw("same");

  //Ploting histograms
  c_thetape->cd(); h_thetap->Draw(); h_thetae->Draw("same"); 
  c_theta->cd(); h_theta->Draw();
  c_Ee->cd();h_Ee->Draw();
  c_Ep->cd();h_Ep->Draw();
}
