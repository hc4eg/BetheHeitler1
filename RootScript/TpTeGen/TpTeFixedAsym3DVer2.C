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

//Ver2: means the condition when generating e+,e- pairs are using following range:
// 15.0 MeV < Ee, Ep < 45.0 MeV,  0 deg < abs(Thetae), abs(Thetap)< 20 deg

// Usage of the file: Computing following situation: A gamma photon in x direction hit a rest nucleus in lab reference frame.
// Then a e+ e- pair is generated with recoiled nucleus moving in photon direction.
// The electron/positron is generated with uniform random Energy and Theta (Angle relative to x-axis in xy plane).
// Then compute the distribution of energy and theta, plot histogram.

// 3D case for TpTeFixedAsym2D.C

void TpTeFixedASym3DVer2 (){
  Double_t Ee, Ep;
  Double_t Thetae, Thetap, Phie, Phip, Theta;
  Double_t E, Me;
  Double_t Pp, Ppx, Ppy, Ppz, Pe, Pex, Pey, Pez;
  Double_t ThetaM, Temp;

  Double_t A,a, Tempa;
  // Genertate e+, e- pairs with fixed asymmetry A
  A = 0.5;
  a = (1-A)/(1+A);
  cerr << "Asymmetry = " << A << " . " << endl;
  
  //Number of total events generated
  Int_t N = 50000000;

  E = 60.0;
  Me = 0.511;
  //Me = 10.0;

  //Canvas
  TCanvas * c_thetape = new TCanvas("c_thetape", "Positron and Electron Theta Distribution", 800, 400);
	c_thetape->ToggleEventStatus();
	c_thetape->Iconify();

	
  TCanvas * c_theta = new TCanvas("c_theta", "Opening Angle Distribution", 800, 400);
	c_theta->ToggleEventStatus();
	c_theta->Iconify();
	
	TCanvas * c_solid = new TCanvas("c_solid","Opening solid angle distribution",800, 400);
	c_solid->ToggleEventStatus();
	c_solid->Iconify();

  TCanvas * c_E = new TCanvas("c_E", "Energy Distribution", 800, 400);
	c_E->ToggleEventStatus();
	c_E->Iconify();

       
  TCanvas * c_DelE = new TCanvas("c_Ep", "Energy Difference Distribution", 800, 400);
	c_DelE->ToggleEventStatus();
	c_DelE->Iconify();

	TCanvas * c_2DTheta = new TCanvas("c_2DTheta", "(Thetae, Thetap) distribution; openning angle vs delta: (Ee-Ep)/2 distribution", 1200, 600);
	c_2DTheta->ToggleEventStatus();
	c_2DTheta->Iconify();
	c_2DTheta->Divide(2,1);

  //Histogram
  TH1F * h_thetap = new TH1F("h_thetap", "Thetap Blue",
				300, 0.0, 90.0);
  h_thetap->SetXTitle("Positron Angle(Degree)");
  h_thetap->SetYTitle("Counts");
  
  TH1F * h_thetae = new TH1F("h_thetae", "Thetae",
				300, 0.0, 90.0);
  h_thetae->SetXTitle("Electron Angle(Degree)");
  h_thetae->SetYTitle("Counts");
  h_thetae->SetLineColor(kRed);

  TH1F * h_theta = new TH1F("h_theta", "Opening",
				300, 0.0, 180.0);
  h_theta->SetXTitle("Opening Angle(Degree)");
  h_theta->SetYTitle("Counts");
  h_theta->SetLineColor(kRed);
  
  TH1F * h_solid = new TH1F("h_solid", "Opening Solid Angle",
				300, 0.0, 4*PI);
  h_solid->SetXTitle("Opening Solid Angle (4PI^=1)");
  h_solid->SetYTitle("Counts"); 

  TH1F * h_Ep = new TH1F("h_Ep", "Ep Blue",
				300, 0.0, 60.0);
  h_Ep->SetXTitle("Ep(MeV)");
  h_Ep->SetYTitle("Counts");

  TH1F * h_Ee = new TH1F("h_Ee", "Ee",
				300, 0.0, 60.0);
  h_Ee->SetXTitle("Ee(MeV)");
  h_Ee->SetYTitle("Counts");
  h_Ee->SetLineColor(kRed);

  TH1F * h_DelE = new TH1F("h_DelE", "(Ee-Ep)",
				300, -60.0, 60.0);
  h_DelE->SetXTitle("Ee-Ep(MeV)");
  h_DelE->SetYTitle("Counts");

  //(Thetae, Thetap)
  TH2F* h_TeTp = new TH2F("h_TeTp", "(Thetae, Thetap)", 
			  200, 0., 20., 200, 0., 20.);
  h_TeTp->SetXTitle("Thetae");
  h_TeTp->SetYTitle("Thetap");
  h_TeTp->SetOption("COLZ");

  TH2F* h_DelTheta = new TH2F("h_DelTheta", "Delta (Ee-Ep)/2 vs Opening Angle", 
			      200, 0., 40., 300, -15., 15.);
  h_DelTheta->SetXTitle("Opening Angle");
  h_DelTheta->SetYTitle("Delta");
  h_DelTheta->SetOption("COLZ");

  //Computing asymmetry from histogram
  Double_t Del = 2.0; Int_t Nm = (Int_t)(120.0/Del);
  Double_t del[Nm], N0[8][Nm], A0[8][Nm], AConst[Nm], ThetaP = 40.0;
  
  // Print out which opening angles are plotted in Asymmetry plot
  cerr << "Asymmetry near following angles are plotted:" << endl;
  for (Int_t i =1; i < 8; i++) cerr << (Double_t)(i)*ThetaP/7 << "degrees. ";
  cerr << endl;
  
  for (Int_t j = 0; j < Nm; j++){
    for (Int_t i = 0; i < 8; i++)
      {N0[i][j] = 0.0; A0[i][j] = 0.0;}
    AConst[j] = 0.0;
    del[j] = 0.0;
  }
  
  for (Int_t i = 0; i < N; i++)
    {
      if(i%10000 == 0) cerr << "Starting " << i << "th event." << endl;
      //Following code is to generate e+ e- pair with same distribution.
	do{
	  Tempa = gRandom -> Uniform(0,1+a);
	  

	  // Here the range of Ee distribution can be altered as long as inside 0 to 60 MeV range.
	  if(Tempa <= 1){
	    Ee = gRandom->Uniform(15.0, 30.0);
	  }
	  else{
	    Ee = gRandom->Uniform(30.0, 45.0);
	  }
	  //Assigning angle theta
	  do{
	    Temp = sqrt(((E-Ee)*(E-Ee)-Me*Me)/(Ee*Ee-Me*Me));
	    if( 1 > abs(Temp)){
	      ThetaM = asin(Temp);
	      Thetae = gRandom->Uniform(0.0, ThetaM);
	    }
	    else
	      Thetae = gRandom->Uniform(0.0,PI/2);

	    Phie = gRandom->Uniform(0, 2*PI);
	    
	    Pe = sqrt(Ee*Ee-Me*Me);
	    Pex = Pe*cos(Thetae);
	    Pey = Pe*sin(Thetae)*cos(Phie);
	    Pez = Pe*sin(Thetae)*sin(Phie);
	    
	    Ep = E - Ee;
	    Ppx = sqrt((E-Ee)*(E-Ee)-(Ee*Ee-Me*Me)*sin(Thetae)*sin(Thetae)-Me*Me);
	    //Ppy = -sqrt(Ee*Ee-Me*Me)*sin(Thetae)*cos(Phie);
	    //Ppz = -sqrt(Ee*Ee-Me*Me)*sin(Thetae)*sin(Phie);
	    Ppy = -Pey;
	    Ppz = -Pez;

	    Thetap = atan((sqrt(Ee*Ee-Me*Me)*sin(Thetae))/Ppx);
	    if(Phie < PI) Phip = Phie + PI;
	    else if(Phie > PI) Phip = Phie - PI;

	    // Method 1, computing Opening Angle. SP : scalar product
	    Double_t SP = cos(Thetae)*cos(Thetap)+sin(Thetae)*sin(Thetap)*(cos(Phie)*cos(Phip)+sin(Phie)*sin(Phip));
	    // Openging angle in unit rad
	    Theta = (acos(SP));
	    // Method 1, since Phip = Phie + PI or Phip = Phie - PI, then SP = cos(Thetae + Thetap). And Opening angle Theta = Thetae + Thetap.
	    //Theta = (Thetae + Thetap);

	    // Convert rad to deg unit.
	    Thetap = Thetap*180/PI;
	    Thetae = Thetae*180/PI;
	    Theta = Theta*180/PI;
	  }while(abs(Thetae*1.0) > 20.0 || abs(Thetap*1.0) > 20.0);
	}while(  ((E-Ee)*(E-Ee)-Pey*Pey-Me*Me < 0.0) || Ee < Me);
	
	//(Ee < 15.0) || (Ee > 45.0) || Thetae > 20.0 || Thetap > 20.0 || Thetae < -20.0 || Thetap < -20.0  );
	//     }
      
      // Condition more like simulation with OR contition (require a range of energy to go through magnetic field, one side detector hit)
	//if((abs(Thetae) < 15.0 && (Ee < 40.0) && (Ee > 20.0)) || (abs(Thetap) < 15.0 && (Ep < 40.0) && (Ep > 20.0))){
      
      // Condition more like simulation with AND condition (require both side detector hit)
      //if((abs(Thetae*1.0) < 20.0 && (Ee < 45.0) && (Ee > 15.0)) && (abs(Thetap*1.0) < 20.0 && (Ep < 45.0) && (Ep > 15.0))){	
	
      // Looser OR requirement, simulate events go through yoke of magnet
      //if(abs(Thetae) < 15.0 || abs(Thetap) < 15.0){

      // Looser AND requirement, simulate events go through yoke of magnet
      //if(abs(Thetae*1.0) < 20.0 && abs(Thetap*1.0) < 20.0){

	for (Int_t n = 0; n < Nm; n++)
	  {

	    // Counting total N(Ee-Ep) and differential N(Ee-Ep, theta):
	    if ( (Ee-Ep) > (-60.0 + (Double_t)(n*Del/1.0)) && (Ee-Ep) <= (-60.0 + (Double_t)((n+1)*Del/1.0)))
	      {
		N0[0][n]++;
		for(Int_t m = 1; m < 8; m++){
		  if ( Theta <= (m*ThetaP/7.0+1.0) && Theta > (m*ThetaP/7.0-1.0))
		    N0[m][n]++;}
	      }
	  }// comment out if any of the 'if' condtion above is commented

	

	//if(abs(Thetae*1.0) < 20.0 && abs(Thetap*1.0) < 20.0){	  
	h_thetap->Fill(Thetap);
	h_thetae->Fill(Thetae);
	h_TeTp->Fill(Thetae, Thetap);
	h_DelTheta->Fill(Theta, (Ee-Ep)/2.);
	//h_theta->Fill(abs(Thetap-Thetae));
	h_theta->Fill(Theta);
	//h_solid->Fill(2*PI*sin(abs(Thetae)*PI/180));
	h_solid->Fill(2*PI*(1-cos(abs(Thetae)*PI/180.)));
	h_Ee->Fill(Ee);
	h_Ep->Fill(Ep);
	h_DelE->Fill(Ee-Ep);
	//}// Comment with any of 4 if statement above for(Int_t n = 0; n < Nm; n++)
    }

  //Compute asymmetry
  for (Int_t n = 0; n < Nm; n++){
    del[n] = -30.0 + Del/4.0 + (Double_t)(n*Del/2.0);
    A0[0][n] = (N0[0][n]-N0[0][Nm-n-1])/(N0[0][n]+N0[0][Nm-n-1]);
    if (N0[0][n] == N0[0][Nm-n-1]) A0[0][n] = 0.0;
    if( N0[0][n] <= (N/(Nm*5.0)) || N0[0][Nm-n-1] <= (N/(Nm*5.0))) A0[0][n] = 0.0;
    
    for(Int_t m = 1; m < 8; m++){
      A0[m][n] = (N0[m][n]-N0[m][Nm-n-1])/(N0[m][n]+N0[m][Nm-n-1]);
      if (N0[m][n] == N0[m][Nm-n-1]) A0[m][n] = 0.0;    
      if (N0[m][n] == 0.0 || N0[m][Nm-n-1] == 0.0) A0[m][n] = 0.0;}

    if(n < Nm/2.0) AConst[n] = A;
    else AConst[n] = -A;
  }
  
  TCanvas * c_EpDel = new TCanvas("c_EpDel", "Asymmetry Epsilon vs Delta from electron",1600,900);
  c_EpDel->ToggleEventStatus();
  c_EpDel->Iconify();

  //Histograms of asymmetry
  vector<TGraph*> g_EpDel;
  for (Int_t n = 0; n < 8; n++){
    g_EpDel.push_back(new TGraph(Nm, del, A0[n]));
    g_EpDel.at(n)->SetLineColor(n+1);
  }
  
  TGraph* g_Const = new TGraph(Nm, del, AConst);
  g_Const->SetLineColor(kBlack);
  
  c_EpDel->cd(); 
  g_EpDel.at(5)->Draw();
  //g_EpDel.at(0)->Draw();
  for (unsigned n = 0; n < g_EpDel.size() ; n++){
    g_EpDel.at(n)->Draw("same");}
  g_Const->Draw("same");

  //Ploting histograms
  c_thetape->cd(); h_thetap->Draw(); h_thetae->Draw("same");
  c_theta->cd(); h_theta->Draw();
  c_solid->cd(); h_solid->Draw();
  c_E->cd();h_Ep->Draw();h_Ee->Draw("same");
  c_DelE->cd();h_DelE->Draw();

  c_2DTheta->cd(1); h_TeTp->Draw(); c_2DTheta->cd(2); h_DelTheta->Draw();
}
