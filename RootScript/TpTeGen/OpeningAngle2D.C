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

// Testing what distribution of Opening Angle will look like in 2D case
// In such condition: Thetae Thetap uniformly distributed from -90 to 90 degree, independent.

void OpeningAngle(){

  Double_t Thetae, Thetap, Theta;
  Int_t N = 5000000;

  TCanvas* c_theta = new TCanvas("c_theta", "Opening Angle Distribution", 800, 400);
	c_theta->ToggleEventStatus();
	c_theta->Iconify();

  TH1F * h_theta = new TH1F("h_theta", "Opening",
				300, 0.0, 180.0);
  h_theta->SetXTitle("Opening Angle(Degree)");
  h_theta->SetYTitle("Counts");


  for (Int_t n = 0; n < N; n++){
    Thetae = gRandom->Uniform(-90,90);
    Thetap = gRandom->Uniform(-90,90);
    Theta = abs(Thetae - Thetap);
    h_theta->Fill(Theta);
  }

  c_theta->cd(); h_theta->Draw();
}
