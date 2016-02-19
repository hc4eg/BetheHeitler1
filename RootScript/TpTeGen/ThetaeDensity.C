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

// Compute and draw density function of Thetae distribution Correspond to TpTeFixedAsym2.C using Monte Carlo integral

void ThetaeDensity(){
  Double_t A = 0.3, a = (1-A)/(1+A), C = 1/(30*(1+a)), E = 60.0, Me = 0.511, Area;
  cerr << "a = " << a << " . " << "C = " << C << endl;
  Int_t N = 150, Nmont = 1000000, nmont = 0;
  Double_t Thetae[N], ThetaeRad[N], fThetae[N], Eemax;
  Double_t EeTheta, fEeTheta, fmax;

  for(Int_t n = 0; n < N; n++)
    {
      Thetae[n] = 0.3+0.6*n;
      ThetaeRad[n] = Thetae[n]*PI/180.0;
      fThetae[n] = 0.0;
    }

  for(Int_t n = 0; n < N; n++)
    {
      Eemax = (E - sqrt(pow(sin(ThetaeRad[n]),2.0)*E*E+pow(cos(ThetaeRad[n]),4.0)*Me*Me))/(pow(cos(ThetaeRad[n]),2.0));
      fmax = 1.0/(2.0*asin(sqrt(((E-Eemax)*(E-Eemax)-Me*Me)/(Eemax*Eemax-Me*Me))));
      //cerr << "ThetaMax = " << asin(sqrt(((E-Eemax)*(E-Eemax)-Me*Me)/(Eemax*Eemax-Me*Me)))*180/PI << " . " << "Eemax = "<< Eemax << " . " << "fmax[" << n << "] = " << fmax << " . ";
      Area = (Eemax-E/2.0)*fmax;
      //cerr << "Area = " << Area << " . ";
	
      nmont = 0;

      for(Int_t m = 0; m < Nmont; m++){
	EeTheta = gRandom->Uniform(E/2.0, Eemax);
	fEeTheta = 1.0/(2.0*asin(sqrt(((E-EeTheta)*(E-EeTheta)-Me*Me)/(EeTheta*EeTheta-Me*Me))));
	//if (m == Nmont/2) cerr << "EeTheta = " << EeTheta << " . "<< "fEeTheta = " << fEeTheta << " . ";
	Double_t Rand = gRandom->Uniform(0.0,fmax);
	if (Rand < fEeTheta) nmont++;
      }
      
      //cerr << "nomt = "  << nmont << ". ";

      fThetae[n] = ((Double_t)nmont/(Double_t)Nmont)*Area*a*C + C*30.0/PI;

      //cerr << "fThetae["<< n << "] = " << fThetae[n] << "." << endl;
    }
  
  Double_t Sum = 0.0;
  for(Int_t n = 0; n < N ; n++){
    Sum += 0.6*PI/180*fThetae[n];
  }
  cerr << "Sum = " << Sum << endl;

  TCanvas* C_fThetae = new TCanvas("C_fThetae","Thetae distribution", 1600, 900);
  C_fThetae->ToggleEventStatus();
  C_fThetae->Iconify();
  TGraph* G_fThetae = new TGraph(N, Thetae, fThetae);
  C_fThetae->cd(); G_fThetae->Draw();
}
