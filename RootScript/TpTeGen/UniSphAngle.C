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

void UniSphAngle (){

  Double_t Theta, R;
  Int_t N = 50000000;

  //canvas
  TCanvas * c_solid = new TCanvas("c_solid","Opening solid angle distribution",800, 400);
  c_solid->ToggleEventStatus();
  c_solid->Iconify();

  TCanvas * c_theta = new TCanvas("c_theta","Angle distribution",800, 400);
  c_theta->ToggleEventStatus();
  c_theta->Iconify();

  //histogram
  TH1F * h_theta = new TH1F("h_theta", "Theta",
				100, 0.0, PI);
  h_theta->SetXTitle("Angle");
  h_theta->SetYTitle("Counts");
  
  TH1F * h_solid = new TH1F("h_solid", "Opening Solid Angle",
			    300, 0.0, 1.);//4*PI);
  h_solid->SetXTitle("Opening Solid Angle (4PI^=1)");
  h_solid->SetYTitle("Counts");

  for (Int_t n = 0; n < N; n++){
    //R = gRandom -> Uniform(0,1);
    //Theta = acos(-2*R+1);
    Theta = gRandom -> Uniform(0, PI);
    h_theta->Fill(Theta);
    //h_solid->Fill(R);
    h_solid->Fill(0.5*(1-cos(Theta)));
  }
  c_solid->cd(); h_solid->Draw();
  c_theta->cd(); h_theta->Draw();
}
