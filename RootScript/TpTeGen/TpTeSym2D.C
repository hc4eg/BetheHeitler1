#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include "TFile.h"
#include "TTree.h"
#include "TRandom.h"
#include "TCanvas.h"
#include "TH1.h"
#include "TH2.h"
#include <cmath>

using namespace std;

#define PI 3.1415926

// Usage of the file: Computing following situation: A gamma photon in x direction hit a rest nucleus in lab reference frame.
// Then a e+ e- pair is generated with recoiled nucleus moving in photon direction.
// The electron/positron is generated with uniform random Energy and Theta (Angle relative to x-axis in xy plane).
// Then compute the distribution of energy and theta, plot histogram.


void TpTeSym (){
  Double_t Ee, Ep;
  Double_t Thetae, Thetap, Phie, Phip;
  Double_t E,Me;
  Double_t Pp,Ppx,Ppy,Ppz,Pe, Pex,Pey,Pez;

  Int_t N = 500000;

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
  TCanvas * c_thetap = new TCanvas("c_thetap", "Theta Positron Distribution", 800, 400);
	c_thetap->ToggleEventStatus();
	c_thetap->Iconify();

  TCanvas * c_thetae = new TCanvas("c_thetae", "Theta Electron Distribution", 800, 400);
	c_thetae->ToggleEventStatus();
	c_thetae->Iconify();

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

  TH1F * h_Ep = new TH1F("h_Ep", "Ep",
				300, 0.0, 60.0);
  h_Ep->SetXTitle("Ep(MeV)");
  h_Ep->SetYTitle("Counts");

  TH1F * h_Ee = new TH1F("h_Ee", "Ee",
				300, 0.0, 60.0);
  h_Ee->SetXTitle("Ee(MeV)");
  h_Ee->SetYTitle("Counts");


  
  for (Int_t i = 0; i < N; i++)
    {
      // Following is a try of generate e+/e- pair in 3D solid angle with uniform solid angle distribution.
      /*
      do{
	Ee = gRandom->Uniform(0.0, 60.0);
	Phie = gRandom->Uniform(0, 2*PI);
	Thetae = gRandom->Uniform();
	Thetae = acos(-2*Thetae+1);
	
	Pe = sqrt(Ee*Ee-Me*Me);
	Pex = Pe*cos(Thetae);
	Pey = Pe*sin(Thetae)*cos(Phie);
	Pez = Pe*sin(Thetae)*sin(Phie);

	Ep = E-Ee;
	Ppy = -Pey;
	Ppz = -Pez;
	Ppx = sqrt((E-Ee)*(E-Ee)-Pey*Pey-Pez*Pez-Me*Me);
      }while (Ee < Me || ((E-Ee)*(E-Ee)-Pey*Pey-Pez*Pez-Me*Me) < 0.0 || Thetae>(PI/2));
      
      //RMS for Thetae is 42.3 deg
      Thetap = atan2(sqrt(Ppy*Ppy+Ppz*Ppz),Ppx);
      Thetap = Thetap*180/PI;
      Thetae = Thetae*180/PI;
      */

      //Compare to Monte Carlo to see the Thetae is generated correctly
      //RMS~39.2deg
      /*
	Thetae = gRandom->Uniform();
	Thetae = acos(-2*Thetae+1);
	Thetae = Thetae*180/PI;
	Double_t Ran;
	do{
	  Thetap = gRandom->Uniform(0,PI);
	  Ran = gRandom->Uniform(0,1);
	}while((sin(Thetap))<Ran);
	Thetap = Thetap*180/PI;
      */

      
      //Following code is to generate e+ e- pair with same distribution.
      //if(i%2 == 0)
	  {
	    do{
	      Ee = gRandom->Uniform(0.0, 60.0);
	      Thetae = gRandom->Uniform(-PI/2, PI/2);
	      //Thetae = gRandom->Uniform(-PI/4, PI/4);

	      Pe = sqrt(Ee*Ee-Me*Me);
	      Pex = Pe*cos(Thetae);
	      Pey = Pe*sin(Thetae);

	      Ep = E - Ee;
	      Ppx = sqrt((E-Ee)*(E-Ee)-(Ee*Ee-Me*Me)*sin(Thetae)*sin(Thetae)-Me*Me);
	      Ppy = -sqrt(Ee*Ee-Me*Me)*sin(Thetae);
	    }while (Ee < Me || ((E-Ee)*(E-Ee)-Pey*Pey-Me*Me) < 0.0 || (E-Ee) < Me);
	    
	    Thetap = atan2(Ppy,Ppx);
	    Thetap = Thetap*180/PI;
	    Thetae = Thetae*180/PI;
	  }
	  /*
       	else
	  {
	    do{
	      Ep = gRandom->Uniform(0.0, 60.0);
	      Thetap = gRandom->Uniform(-PI/2, PI/2);
	      //Thetap = gRandom->Uniform(-PI/4, PI/4);

	      Pp = sqrt(Ep*Ep-Me*Me);
	      Ppx = Pp*cos(Thetap);
	      Ppy = Pp*sin(Thetap);

	      Ee = E - Ep;
	      Pex = sqrt((E-Ep)*(E-Ep)-(Ep*Ep-Me*Me)*sin(Thetap)*sin(Thetap)-Me*Me);
	      Pey = -sqrt(Ep*Ep-Me*Me)*sin(Thetap);
	    }while (Ep < Me || ((E-Ep)*(E-Ep)-Ppy*Ppy-Me*Me) < 0.0 || (E-Ep) < Me);
	    
	    Thetae = atan2(Pey,Pex);
	    Thetae = Thetae*180/PI;
	    Thetap = Thetap*180/PI;

	  }
	  */

      // Condition more like simulation with OR contition (require a range of energy to go through magnetic field, one side detector hit)
      // if((abs(Thetae) < 15.0 && (Ee < 40.0) && (Ee > 20.0)) || (abs(Thetap) < 15.0 && (Ep < 40.0) && (Ep > 20.0))){
      
      // Condition more like simulation with AND condition (require both side detector hit)
      // if((abs(Thetae) < 15.0 && (Ee < 40.0) && (Ee > 20.0)) && (abs(Thetap) < 15.0 && (Ep < 40.0) && (Ep > 20.0))){	
	
      // Looser OR requirement, simulate events go through yoke of magnet
      //if(abs(Thetae) < 15.0 || abs(Thetap) < 15.0){

      // Looser AND requirement, simulate events go through yoke of magnet
      //if(abs(Thetae) < 15.0 && abs(Thetap) < 15.0){
      h_thetap->Fill(Thetap);
      h_thetae->Fill(Thetae);
    
	h_Ee->Fill(Ee);
	h_Ep->Fill(Ep);
	// }

    }

  c_thetap->cd();h_thetap->Draw();
  c_thetae->cd();h_thetae->Draw();
  c_Ee->cd();h_Ee->Draw();
  c_Ep->cd();h_Ep->Draw();
}
