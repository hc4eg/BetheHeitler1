#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include "TFile.h"
#include "TTree.h"
#include "TRandom.h"
#include "TCanvas.h"
#include "TGraphErrors.h"
#include "TH1.h"
#include "TH2.h"
#include <cmath>

using namespace std;

#define PI 3.14159265

// Usage of the file: Computing following situation: A gamma photon in x direction hit a rest nucleus in lab reference frame.
// Then a e+ e- pair is generated with recoiled nucleus moving in photon direction.
// The electron/positron is generated with uniform random Energy and Theta (Angle relative to x-axis in xy plane).
// Then compute the distribution of energy and theta, plot histogram.

// 3D case for TpTeFixedAsym2D.C

void Accecptance (){
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
  Int_t N = 20000000;
  //Int_t N = 50000;

  E = 60.0;
  Me = 0.511;
  //Me = 10.0;

  //Canvas
  /*
  TCanvas * c_thetape = new TCanvas("c_thetape", "Positron and Electron Theta Distribution", 800, 400);
	c_thetape->ToggleEventStatus();
	c_thetape->Iconify();

	
  TCanvas * c_theta = new TCanvas("c_theta", "Opening Angle Distribution", 800, 400);
	c_theta->ToggleEventStatus();
	c_theta->Iconify();
	
	TCanvas * c_solid = new TCanvas("c_solid","Opening solid angle distribution",800, 400);
	c_solid->ToggleEventStatus();
	c_solid->Iconify();
    
  TCanvas * c_DelE = new TCanvas("c_Ep", "Energy Difference Distribution", 800, 400);
	c_DelE->ToggleEventStatus();
	c_DelE->Iconify();
	*/

  // Canvas to plot Ee distributions after sequential geometry constraints
  TCanvas * c_Ee = new TCanvas("c_Ee", "e- Energy Distribution", 1200, 800);
  c_Ee->ToggleEventStatus();
  c_Ee->Iconify();
  c_Ee->Divide(3,2);
  
  // Same for Ep
  TCanvas * c_Ep = new TCanvas("c_Ep", "e+ Energy Distribution", 1200, 800);
  c_Ep->ToggleEventStatus();
  c_Ep->Iconify();
  c_Ep->Divide(3,2);

  //'AND' condition Ee-Ep plots
  TCanvas * c_EeEp = new TCanvas("c_EeEp", "Ee-Ep distribution", 1200, 800);
  c_EeEp->ToggleEventStatus();
  c_EeEp->Iconify();
  c_EeEp->Divide(2,2);

  // (Thetae, Phie) and (Thetap, Phip) distribution
  TCanvas * c_ThetaPhi = new TCanvas("c_ThetaPhi", "(Theta,Phi) distribution, left e+, right e-", 1200, 800);
  c_ThetaPhi->ToggleEventStatus();
  c_ThetaPhi->Iconify();
  c_ThetaPhi->Divide(2,2);
  // X vs Y distribution in VDC0 and VDC1 in either AND or OR condition
  TCanvas * c_eXY = new TCanvas("c_eXY", "(X,Y) distribution, up OR condition down AND condition, left VDC0, right VDC1", 1200, 800);
  c_eXY->ToggleEventStatus();
  c_eXY->Iconify();
  c_eXY->Divide(2,2);

  TCanvas * c_pXY = new TCanvas("c_pXY", "(X,Y) distribution, up OR condition down AND condition, left VDC0, right VDC1", 1200, 800);
  c_pXY->ToggleEventStatus();
  c_pXY->Iconify();
  c_pXY->Divide(2,2);

  TCanvas * c_Asym = new TCanvas("c_Asym", "Asymmetry, [0] use OR condition E, [1] use VDC1 E, [2] use Hodo E", 1200, 800);
  c_Asym->ToggleEventStatus();
  c_Asym->Iconify();
  c_Asym->Divide(2,2);
 
  //Histogram
	/*
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
				300, 0.0, 2*PI);
  h_solid->SetXTitle("Opening Solid Angle (4PI^=1)");
  h_solid->SetYTitle("Counts"); 

  TH1F * h_Ep = new TH1F("h_Ep", "Ep Blue",
				300, 0.0, 60.0);
  h_Ep->SetXTitle("Ep(MeV)");
  h_Ep->SetYTitle("Counts");

  TH1F* h_Ee[6];
  for (Int_t i = 0; i < 6; i++){
    h_Ee[i] = new TH1F("h_Ee", "Ee",300, 0.0, 60.0);
    h_Ee[i]->SetXTitle("Ee(MeV)");
    h_Ee[i]->SetYTitle("Counts");
    h_Ee[i]->SetLineColor(kRed);}

  TH1F * h_DelE = new TH1F("h_DelE", "(Ee-Ep)",
				300, -60.0, 60.0);
  h_DelE->SetXTitle("Ee-Ep(MeV)");
  h_DelE->SetYTitle("Counts");*/

  // Histograms to plot Ee after all geometry constraints
  TH1F* h_Ee[6];
  for (Int_t i = 0; i < 6; i++){
    h_Ee[i] = new TH1F("h_Ee", "Ee",300, 0.0, 60.0);
    h_Ee[i]->SetXTitle("Ee(MeV)");
    h_Ee[i]->SetYTitle("Counts");
    //  h_Ee[i]->SetLineColor(kRed);
  }
  TH1F* h_Ep[6];
  for (Int_t i = 0; i < 6; i++){
    h_Ep[i] = new TH1F("h_Ep", "Ep",300, 0.0, 60.0);
    h_Ep[i]->SetXTitle("Ep(MeV)");
    h_Ep[i]->SetYTitle("Counts");
  }
  TH1F* h_EeEp[3];
  for (Int_t i = 0; i < 3; i++){
    h_EeEp[i] = new TH1F("h_EeEp", Form("Ee-Ep[%d]",i),400, -40.0, 40.0);
    h_EeEp[i]->SetXTitle("Ee-Ep(MeV)");
    h_EeEp[i]->SetYTitle("Counts");
  }

  // Data used to draw asymmetry curves by TGraphErrors
  Int_t Num = 60;
  Float_t Count[4][Num]; Float_t Asym[4][Num], AsymErr[4][Num]; Float_t Del[Num]; Float_t delta = 120.0/(Float_t)Num;
  for(Int_t i = 0; i < 4; i++)
    for(Int_t j = 0; j < Num; j++){
      Count[i][j] = 0.; Asym[i][j] = 0.; AsymErr[i][j] = 0.;
    }
  for(Int_t i = 0; i < Num; i++) Del[i] = -30.+delta/4.+i*delta/2.;
  
  TH2F* h_ThetaPhiOR[2];
  for (Int_t i = 0; i < 2; i++){
    h_ThetaPhiOR[i] = new TH2F("h_ThetaPhiOR", Form("( Theta, Phi)[%d]",i), 200, 0., 20. ,360, 0., 360.0);
    h_ThetaPhiOR[i]->SetXTitle("Theta");
    h_ThetaPhiOR[i]->SetYTitle("Phi");
    h_ThetaPhiOR[i]->SetOption("COLZ");
  }
  TH2F* h_ThetaPhiAND[2];
  for (Int_t i = 0; i < 2; i++){
    h_ThetaPhiAND[i] = new TH2F("h_ThetaPhiAND", Form("(Theta, Phi)[%d]",i), 200, 0., 20. ,360, 0., 360.0);
    h_ThetaPhiAND[i]->SetXTitle("Theta");
    h_ThetaPhiAND[i]->SetYTitle("Phi");
    h_ThetaPhiAND[i]->SetOption("COLZ");
  }

  
  // 1st []: 0 OR condition, 1 AND condition; 2nd []: 0 VDC0, 1 VDC1
  TH2F* h_eXY[2][2];
  for(Int_t i = 0; i < 2; i++)
    for(Int_t j = 0; j < 2; j++){
      h_eXY[i][j] = new TH2F("h_eXY", Form("(X, Y)[%d]",i), 205, -1.025, 1.025, 54, -0.27, 0.27);
      h_eXY[i][j]->SetXTitle("X");
      h_eXY[i][j]->SetYTitle("Y");
      h_eXY[i][j]->SetOption("COLZ");
    }

  TH2F* h_pXY[2][2];
  for(Int_t i = 0; i < 2; i++)
    for(Int_t j = 0; j < 2; j++){
      h_pXY[i][j] = new TH2F("h_pXY", Form("(X, Y)[%d]",i), 205, -1.025, 1.025, 54, -0.27, 0.27);
      h_pXY[i][j]->SetXTitle("X");
      h_pXY[i][j]->SetYTitle("Y");
      h_pXY[i][j]->SetOption("COLZ");
    }
  // Loop over each single "Pair production"
  for (Int_t i = 0; i < N; i++)
    {
      if( i%10000 == 0) cerr << "Start the " << i << "th event." << endl; 
      //Following code is to generate e+ e- pair with same distribution.
      do{
	Tempa = gRandom -> Uniform(0,1+a);
	  

	// Here the range of Ee distribution can be altered as long as inside 0 to 60 MeV range.
	if(Tempa <= 1){
	  Ee = gRandom->Uniform(15., 30.0);
	}
	else{
	  Ee = gRandom->Uniform(30.0, 45.);
	}
	  //Assigning angle theta
	//do ... while condition that all initial angle is less than 20 degrees.
	//  do{

		do{
		  Temp = sqrt(((E-Ee)*(E-Ee)-Me*Me)/(Ee*Ee-Me*Me));	  
		  if(-1 < Temp && 1 > Temp){
		    ThetaM = asin(Temp);
		    Thetae = gRandom->Uniform(0.0, ThetaM);
		  }
		  else
		    Thetae = gRandom->Uniform(0.0,PI/2);
	    
		  Phie = gRandom->Uniform(0.0, 2*PI);

		  Pe = sqrt(Ee*Ee-Me*Me);
		  Pex = Pe*cos(Thetae);
		  Pey = Pe*sin(Thetae)*cos(Phie);
		  Pez = Pe*sin(Thetae)*sin(Phie);
	    
		  Ep = E - Ee;
		  Ppx = sqrt((E-Ee)*(E-Ee)-(Ee*Ee-Me*Me)*sin(Thetae)*sin(Thetae)-Me*Me);
		  Ppy = -Pey;
		  Ppz = -Pez;
		  
		  //Note: Thetap here is always in 0 to 90deg
		  Thetap = atan(+Pe*sin(Thetae)/Ppx);
		  //Make positron opposite direction than electron, and lies in 0 to 360deg
		  if (Phie >= 0. && Phie < PI)
		    Phip = Phie + PI;
		  else if ( Phie >= PI && Phie < 2*PI)
		    Phip = Phie - PI;

	// Method 1, computing Opening Angle
		  Double_t SP = cos(Thetae)*cos(Thetap)+sin(Thetae)*sin(Thetap)*(cos(Phie)*cos(Phip)+sin(Phie)*sin(Phip));
		  Theta = (acos(SP));
		  Thetap = Thetap*180/PI;
		  Thetae = Thetae*180/PI;
		  //cerr << "Thetap = " << Thetap << ", Thetae = " << Thetae << endl;
		  Theta = Theta*180/PI;
	  
		}while(abs(Thetae) > 20.0 || abs(Thetap) > 20.0);

	    
	// Method 1, computing Opening Angle
	    
	// Method 1, since Phip = Phie + PI, then SP = cos(Thetae + Thetap). And Opening angle Theta = Thetae + Thetap.
	//Theta = e(Thetae + Thetap);

	//Convert into degrees.
	    // }while(abs(Thetae*1.0) > 20.0 || abs(Thetap*1.0) > 20.0);
      }while( ((E-Ee)*(E-Ee)-Pey*Pey-Me*Me < 0.0) || Ee < Me);
	     //(Ee < 15.0) || (Ee > 45.0) || Thetae > 20.0 || Thetap > 20.0 || Thetae < -20.0 || Thetap < -20.0  );
	    //     }

      // Using various flags to indicate if e- goes through a specific geometry constraint
      Int_t Flage[6] = {0,0,0,0,0,0};
      // Now check detectors acceptance from magnets, VDC till Paddles.
      // Looser OR requirement, simulate events go through yoke of magnet
      // 1st cut: the hole at magnet
      if(abs(Thetae) < 14.67 || abs(Thetap) < 14.67){
	Flage[0] = 1;
	//Check only e- Acceptance and only Ee is used
	h_Ee[0]->Fill(Ee);
      }//if(abs(Thetae) < 15.0 || abs(Thetap) < 15.0)

      //2nd Cut, check if particle make out magnetic field
      //list of constants for 2nd cut computation, most use SI units
      Float_t in = 2.54;
      // R: magnetic field radius, D: distance between center of magnet and primary vertex
      // C: speed of light, B: center magnetic field
      Float_t R = 18.5*in/200., D = 42.26/100., C = 299792458, B = 0.4536;
      Thetae = Thetae*PI/180.;
      // Compute Position and velocity at Primary vertex
      // Note: all array: [0] for e+, [1] for e-
      Float_t XA[2], YA[2], ZA[2];
      XA[1] = -D, YA[1] = 0., ZA[1] = 0.;
      Float_t VXA[2], VYA[2], VZA[2];
      VXA[1] = (-1)*C*cos(Thetae), VYA[1] = C*sin(Thetae)*cos(Phie), VZA[1] = C*sin(Thetae)*sin(Phie);
      // Compute Position and velocity at Point when particle enters B field
      Float_t K[2];
      K[1]= tan(Thetae)*cos(Phie);
      //cerr << "K[1] = " << K[1] << " sqrt((R*R*(1+K[1]*K[1])-K[1]*K[1]*D*D) " << sqrt((R*R*(1+K[1]*K[1])-K[1]*K[1]*D*D)) <<  endl;
      Float_t XB[2], YB[2], ZB[2];
      XB[1] = (D*K[1]*K[1]+sqrt((R*R*(1.+K[1]*K[1])-K[1]*K[1]*D*D)))/(1.+K[1]*K[1]);
      YB[1] = (D-XB[1])*K[1];
      ZB[1] = YB[1]*tan(Phie);
      //cerr << "Thetae = " << Thetae*180./PI << "deg. Phie = " << Phie*180./PI << "deg." << endl; 
      //cerr << " (XB[1], YB[1], ZB[1]) = (" << XB[1] << " , " << YB[1] << " , " << ZB[1] << " ) meters." << endl;   
      
      //Center of the Helix(Trajectory of e- inside uniform B field)
      // r: radius of projection of helix in XY plane, Alpha: angle between initial velocity in XY plane projection and X direction.
      Float_t r[2], Alpha[2]; 
      r[1]= Ee*1000000.*sqrt(1-sin(Thetae)*sin(Phie)*sin(Thetae)*sin(Phie))/(C*B); Alpha[1] = atan(tan(Thetae)*cos(Phie));
      Float_t XH[2], YH[2], ZH[2];
      XH[1] = XB[1] + r[1]*sin(Alpha[1]); YH[1] = YB[1] + r[1]*cos(Alpha[1]); ZH[1] = ZB[1];

      //Compute Time, Position and Velocity at Point when particle leaves B field
      Float_t TC[2];
      TC[1] = 2*r[1]*acos((XH[1]*XH[1]+YH[1]*YH[1]+r[1]*r[1]-R*R)/(2*r[1]*sqrt(XH[1]*XH[1]+YH[1]*YH[1])))/(C*sqrt(1-sin(Thetae)*sin(Phie)*sin(Thetae)*sin(Phie)));
      Float_t Omega[2], OmegaTC[2];
      Omega[1] = C*C*B/(Ee*1000000); OmegaTC[1] = 2.*acos((XH[1]*XH[1]+YH[1]*YH[1]+r[1]*r[1]-R*R)/(2*r[1]*sqrt(XH[1]*XH[1]+YH[1]*YH[1])));
      Float_t XC[2], YC[2], ZC[2];
      XC[1] = XH[1] - r[1]*sin(OmegaTC[1]+Alpha[1]); YC[1] = YH[1] - r[1]*cos(OmegaTC[1]+Alpha[1]); ZC[1] = ZH[1] + C*TC[1]*sin(Thetae)*sin(Phie);
      //cerr << " r[1] = " << r[1] << "meters. Alpha[1] = " << Alpha[1]*180./PI << "degrees." << endl; 
      //cerr << " (XC[1], YC[1], ZC[1]) = (" << XC[1] << " , " << YC[1] << " , " << ZC[1] << " ) meters." << endl;
      Float_t PoleGap = 0.0508;

      // Check if e- hits poles of the magnet
      if( Flage[0]==1 && abs(ZC[1]) <= (PoleGap/2.) ){
	Flage[1] = 1;
	//Check e- Acceptance and only Ee is used
	h_Ee[1]->Fill(Ee);}

      // 3rd cut: Check if e- hits poles of the magnet
      // Compute e- velocity and trajectory when it exits magnet
      // VP: velocity parallel to horizontal plane XY
      Float_t VP[2];
      VP[1] = C*sqrt(1-sin(Thetae)*sin(Phie)*sin(Thetae)*sin(Phie)) ;
      Float_t VXC[2], VYC[2], VZC[2];
      VXC[1] = (-1.)*VP[1]*cos(OmegaTC[1]+Alpha[1]); VYC[1] = VP[1]*sin(OmegaTC[1]+Alpha[1]); VZC[1] = C*sin(Thetae)*sin(Phie);
      //Compute the point when e- reaches the boundary of side of magnet
      Float_t MagInnerX = 22.75*in/100. , MagY = 22.75*in/100. , MagInnerZ = 11.50*in/100.;
      // Time from exit of B field to side boundary of Magnet
      Float_t TD[2];
      TD[1] = ((MagY/2.) - YC[1])/VYC[1];
      Float_t XD[2], YD[2], ZD[2];
      XD[1] = XC[1]+VXC[1]*TD[1]; YD[1] = MagY/2. ; ZD[1] = ZC[1]+VZC[1]*TD[1];
      //Check if e- hits poles of the magnet
      if( Flage[1] == 1 && abs(ZD[1]) <= MagInnerZ/2. && abs(XD[1]) <= MagInnerX/2. ){
	Flage[2] = 1;
	h_Ee[2]->Fill(Ee);}

      // Note : cut 4,5,6 simplified that hits considered at center plane of detectors
      // 4th cut: Check if e- hits effective area of VDC0
      Float_t VDC0X = 80.5*in/100. , VDC0Y = 0.85, VDC0Z = 21.*in/(100.*sqrt(2));
      // Time from point C to VDC0
      Float_t TE[2];
      TE[1]= (VDC0Y-(YC[1]+ZC[1]))/(VYC[1]+VZC[1]);
      Float_t XE[2], YE[2], ZE[2];
      XE[1] = XC[1]+VXC[1]*TE[1]; YE[1] = YC[1]+VYC[1]*TE[1] ; ZE[1] = ZC[1]+VZC[1]*TE[1];
      if ( Flage[2] == 1 && abs(ZE[1]) <= (VDC0Z/2.) && abs(XE[1]) <= (VDC0X/2.)){
	Flage[3] = 1;
	h_Ee[3]->Fill(Ee);}

      // 5th cut: Check if e- hits effective area of VDC1
      Float_t VDC1X = VDC0X, VDC1Y = VDC0Y + 15.5*in*sqrt(2)/100., VDC1Z = VDC0Z;
      // Time from point C to VDC0
      Float_t TF[2];
      TF[1] = (VDC1Y-(YC[1]+ZC[1]))/(VYC[1]+VZC[1]);
      Float_t XF[2], YF[2], ZF[2];
      XF[1] = XC[1]+VXC[1]*TF[1]; YF[1] = YC[1]+VYC[1]*TF[1] ; ZF[1] = ZC[1]+VZC[1]*TF[1];
      if ( Flage[3] == 1 && abs(ZF[1]) <= (VDC1Z/2.0)&& abs(XF[1]) <= (VDC1X/2.) ){
	Flage[4] = 1;
	h_Ee[4]->Fill(Ee);}
      
      // 6th cut: Check if e- hits effective area of Hodoscope
      Float_t HodoX = 1.17 ,HodoY = 1.85, HodoZ = 0.35/sqrt(2);
      // Time from point C to VDC0
      Float_t TG[2];
      TG[1] = (HodoY-(YC[1]+ZC[1]))/(VYC[1]+VZC[1]);
      Float_t XG[2], YG[2], ZG[2];
      XG[1] = XC[1]+VXC[1]*TG[1]; YG[1] = YC[1]+VYC[1]*TG[1] ; ZG[1] = ZC[1]+VZC[1]*TG[1];
      if ( Flage[4] == 1 && abs(ZG[1]) <= (HodoZ/2.) && abs(XG[1]) <= (HodoX/2.) ){
	Flage[5] = 1;
	h_Ee[5]->Fill(Ee);}















      // Now consider Ep distribution
      // By symmetry of the geometry, only change to make is Thetae->Thetap, Phie->-Phie from e- to e+.
      // Using various flags to indicate if e+ goes through a specific geometry constraint
      Int_t Flagp[6] = {0,0,0,0,0,0};
      // Now check detectors acceptance from magnets, VDC till Paddles.
      // Looser OR requirement, simulate events go through yoke of magnet
      // 1st cut: the hole at magnet
      Thetae = Thetae*180/PI;
      //cerr << " Check point 2: Thetap = " << Thetap << ", Thetae = " << Thetae << endl;
      if(abs(Thetae) < 14.67 || abs(Thetap) < 14.67){
	Flagp[0] = 1;
	//Check only e+ Acceptance and only Ep is used
	h_Ep[0]->Fill(Ep);
      }//if(abs(Thetae) < 15.0 || abs(Thetap) < 15.0)

      //2nd Cut, check if particle make out magnetic field
      //list of constants for 2nd cut computation, most use SI units
      Thetap = Thetap*PI/180.;
      // Compute Position and velocity at Primary vertex
      XA[0] = -D, YA[0] = 0., ZA[0] = 0.;
      VXA[0] = (-1)*C*cos(Thetap), VYA[0] = C*sin(Thetap)*cos(-Phie), VZA[0] = C*sin(Thetap)*sin(-Phie);
      // Compute Position and velocity at Point when particle enters B field
      K[0]= tan(Thetap)*cos(-Phie);
      //cerr << "K[0] = " << K[0] << " sqrt((R*R*(1+K[0]*K[0])-K[0]*K[0]*D*D) " << sqrt((R*R*(1+K[0]*K[0])-K[0]*K[0]*D*D)) <<  endl;
      XB[0] = (D*K[0]*K[0]+sqrt((R*R*(1.+K[0]*K[0])-K[0]*K[0]*D*D)))/(1.+K[0]*K[0]);
      YB[0] = (D-XB[0])*K[0];
      ZB[0] = YB[0]*tan(-Phie);
      //cerr << "Thetap = " << Thetap*180./PI << "deg. -Phie = " << -Phie*180./PI << "deg." << endl; 
      //cerr << " (XB[0], YB[0], ZB[0]) = (" << XB[0] << " , " << YB[0] << " , " << ZB[0] << " ) meters." << endl;   
      
      //Center of the Helix(Trajectory of e+ inside uniform B field)
      // r: radius of projection of helix in XY plane, Alpha: angle between initial velocity in XY plane projection and X direction.
      r[0]= Ep*1000000.*sqrt(1-sin(Thetap)*sin(-Phie)*sin(Thetap)*sin(-Phie))/(C*B); Alpha[0] = atan(tan(Thetap)*cos(-Phie));
      XH[0] = XB[0] + r[0]*sin(Alpha[0]); YH[0] = YB[0] + r[0]*cos(Alpha[0]); ZH[0] = ZB[0];

      //Compute Time, Position and Velocity at Point when particle leaves B field
      TC[0] = 2*r[0]*acos((XH[0]*XH[0]+YH[0]*YH[0]+r[0]*r[0]-R*R)/(2*r[0]*sqrt(XH[0]*XH[0]+YH[0]*YH[0])))/(C*sqrt(1-sin(Thetap)*sin(-Phie)*sin(Thetap)*sin(-Phie)));
      Omega[0] = C*C*B/(Ep*1000000); OmegaTC[0] = 2.*acos((XH[0]*XH[0]+YH[0]*YH[0]+r[0]*r[0]-R*R)/(2*r[0]*sqrt(XH[0]*XH[0]+YH[0]*YH[0])));
      XC[0] = XH[0] - r[0]*sin(OmegaTC[0]+Alpha[0]); YC[0] = YH[0] - r[0]*cos(OmegaTC[0]+Alpha[0]); ZC[0] = ZH[0] + C*TC[0]*sin(Thetap)*sin(-Phie);
      //cerr << " r[0] = " << r[0] << "meters. Alpha[0] = " << Alpha[0]*180./PI << "degrees." << endl; 
      //cerr << " (XC[0], YC[0], ZC[0]) = (" << XC[0] << " , " << YC[0] << " , " << ZC[0] << " ) meters." << endl;



      // Check if e+ hits poles of the magnet
      if( Flagp[0]==1 && abs(ZC[0]) <= (PoleGap/2.) ){
	Flagp[1] = 1;
	//Check e+ Acceptance and only Ep is used
	h_Ep[1]->Fill(Ep);}

      // 3rd cut: Check if e+ hits poles of the magnet
      // Compute e+ velocity and trajectory when it exits magnet
      VP[0] = C*sqrt(1-sin(Thetap)*sin(-Phie)*sin(Thetap)*sin(-Phie)) ;
      VXC[0] = (-1.)*VP[0]*cos(OmegaTC[0]+Alpha[0]); VYC[0] = VP[0]*sin(OmegaTC[0]+Alpha[0]); VZC[0] = C*sin(Thetap)*sin(-Phie);
      //Compute the point when e+ reaches the boundary of side of magnet
      // Time from exit of B field to side boundary of Magnet
      TD[0] = ((MagY/2.) - YC[0])/VYC[0];
      XD[0] = XC[0]+VXC[0]*TD[0]; YD[0] = MagY/2. ; ZD[0] = ZC[0]+VZC[0]*TD[0];
      //Check if e+ hits poles of the magnet
      if( Flagp[1] == 1 && abs(ZD[0]) <= MagInnerZ/2. && abs(XD[0]) <= MagInnerX/2. ){
	Flagp[2] = 1;
	h_Ep[2]->Fill(Ep);}

      // Note : cut 4,5,6 simplified that hits considered at center plane of detectors
      // 4th cut: Check if e+ hits effective area of VDC0
      // Time from point C to VDC0
      TE[0]= (VDC0Y-(YC[0]+ZC[0]))/(VYC[0]+VZC[0]);
      XE[0] = XC[0]+VXC[0]*TE[0]; YE[0] = YC[0]+VYC[0]*TE[0] ; ZE[0] = ZC[0]+VZC[0]*TE[0];
      if ( Flagp[2] == 1 && abs(ZE[0]) <= (VDC0Z/2.) && abs(XE[0]) <= (VDC0X/2.)){
	Flagp[3] = 1;
	h_Ep[3]->Fill(Ep);}

      // 5th cut: Check if e+ hits effective area of VDC1
      // Time from point C to VDC0
      TF[0] = (VDC1Y-(YC[0]+ZC[0]))/(VYC[0]+VZC[0]);
      XF[0] = XC[0]+VXC[0]*TF[0]; YF[0] = YC[0]+VYC[0]*TF[0] ; ZF[0] = ZC[0]+VZC[0]*TF[0];
      if ( Flagp[3] == 1 && abs(ZF[0]) <= (VDC1Z/2.0)&& abs(XF[0]) <= (VDC1X/2.) ){
	Flagp[4] = 1;
	h_Ep[4]->Fill(Ep);}
      
      // 6th cut: Check if e+ hits effective area of Hodoscope
      // Time from point C to VDC0
      TG[0] = (HodoY-(YC[0]+ZC[0]))/(VYC[0]+VZC[0]);
      XG[0] = XC[0]+VXC[0]*TG[0]; YG[0] = YC[0]+VYC[0]*TG[0] ; ZG[0] = ZC[0]+VZC[0]*TG[0];
      if ( Flagp[4] == 1 && abs(ZG[0]) <= (HodoZ/2.) && abs(XG[0]) <= (HodoX/2.) ){
	Flagp[5] = 1;
	h_Ep[5]->Fill(Ep);}
      






      //Now add 'AND' cut condition.
      // PairFlag 1st Index 0 e+, 1 e+; 2nd Index 0 VDC package hit flag, 1 Hodoscope hit flag
      Int_t PairFlag[2][2];
      for (Int_t i = 0; i < 2; i++)
	for(Int_t j = 0; j < 2; j++){
	  PairFlag[i][j] = 0;
	}
      if(Flage[4]==1) PairFlag[1][0] = 1;
      if(Flagp[4]==1) PairFlag[0][0] = 1;
      if(Flage[2]==1 &&  abs(ZG[1]) <= (HodoZ/2.) && abs(XG[1]) <= (HodoX/2.))PairFlag[1][1] = 1;
      if(Flagp[2]==1 &&  abs(ZG[0]) <= (HodoZ/2.) && abs(XG[0]) <= (HodoX/2.))PairFlag[0][1] = 1;
      if(PairFlag[0][0]==1 && PairFlag[1][0]==1){
	h_EeEp[1]->Fill(Ee-Ep);
	// Use this condition to plot asymmetry curve on VDC energy
	for(Int_t i = 0; i < Num ; i++){
	  if( (-60.+ (Float_t)i*delta <= Ee-Ep) && (Ee-Ep < -60.+(Float_t)(i+1)*delta))
	    Count[1][i]++;
	}	
      }
      if(PairFlag[0][1]==1 && PairFlag[1][1]==1){
	h_EeEp[2]->Fill(Ee-Ep);
	// Use this condition to plot asymmetry curve on Hodoscope energy
	for(Int_t i = 0; i < Num ; i++){
	  if( (-60.+ (Float_t)i*delta <= Ee-Ep) && (Ee-Ep < -60.+(Float_t)(i+1)*delta))
	    Count[2][i]++;
	}
      }


      // Check Theta,Phi acceptance for OR condition.
      // If all PairFlag in OR condition Plot (Thetae, Phie) and (Thetap, Phip) distribution
      Float_t TEMPThetae, TEMPThetap, TEMPPhie, TEMPPhip;
      TEMPThetae = Thetae;
      TEMPThetap = Thetap*180/PI;
      TEMPPhie = Phie*180/PI;
      TEMPPhip = Phip*180/PI;
      // if(...) impose OR condition
      if(PairFlag[0][0] || PairFlag[0][1] || PairFlag[1][0] || PairFlag[1][1]){
	h_ThetaPhiOR[0]->Fill(TEMPThetap, TEMPPhip);
	h_ThetaPhiOR[1]->Fill(TEMPThetae, TEMPPhie);
	// Ee-Ep histogram
	h_EeEp[0]->Fill(Ee-Ep);
	// Counts of Ee-Ep for OR condition asymmetry graph.
	for(Int_t i = 0; i < Num ; i++){
	  if( (-60.+ (Float_t)i*delta <= Ee-Ep) && (Ee-Ep < -60.+(Float_t)(i+1)*delta))
	    Count[0][i]++;
	}
	h_eXY[0][0]->Fill(XE[1], ZE[1]*sqrt(2.));
	h_eXY[0][1]->Fill(XF[1], ZF[1]*sqrt(2.));
	h_pXY[0][0]->Fill(XE[0], ZE[0]*sqrt(2.));
	h_pXY[0][1]->Fill(XF[0], ZF[0]*sqrt(2.));
	//cerr << "Thetae = " << TEMPThetae << " , Phie = " << TEMPPhie << endl;
      }

      // Check Theta,Phi acceptance for AND condition, which is near to experiment condition
      // if(...) impose AND condition
      if(PairFlag[0][0] && PairFlag[0][1] && PairFlag[1][0] && PairFlag[1][1]){
	h_ThetaPhiAND[0]->Fill(TEMPThetap, TEMPPhip);
	h_ThetaPhiAND[1]->Fill(TEMPThetae, TEMPPhie);
	// Counts of Ee-Ep AND condition for asymmetry graph
	for(Int_t i = 0; i < Num ; i++){
	  if( (-60.+ (Float_t)i*delta <= Ee-Ep) && (Ee-Ep < -60.+(Float_t)(i+1)*delta))
	    Count[3][i]++;
	}	
	h_eXY[1][0]->Fill(XE[1], ZE[1]*sqrt(2.));
	h_eXY[1][1]->Fill(XF[1], ZF[1]*sqrt(2.));
	h_pXY[1][0]->Fill(XE[0], ZE[0]*sqrt(2.));
	h_pXY[1][1]->Fill(XF[0], ZF[0]*sqrt(2.));
	//cerr << "Thetae = " << TEMPThetae << " , Phie = " << TEMPPhie << endl;
      }
      
    }//for (Int_t i = 0; i < N; i++)


  // Compute asymmetry and add data to TGraphErrors:
  for(Int_t i = 0; i < 4; i++)
    for(Int_t j = 0; j < Num/2; j++){
      if(Count[i][j] != 0. && Count[i][Num-j] != 0.){
	Int_t N1 = Count[i][j], N2 = Count[i][Num-j];
	Asym[i][j] = (Count[i][j]-Count[i][Num-j])/(Count[i][j]+Count[i][Num-j]);
	AsymErr[i][j] = 2.*sqrt(N1*N2)/((N1+N2)*sqrt(N1+N2));
	Asym[i][Num-j] = -Asym[i][j];
	AsymErr[i][Num-j] = AsymErr[i][j];
	cerr << "N[" << i << "][" << j << "] = " << Count[i][j]
	     << ", " << "N[" << i <<"]["<< Num-j << "] = " << Count[i][Num-j] << endl;
      }
      else{
	//Asym[i][j] = 0.; Asym[i][Num-j] = 0.;
      }
    }
  
  TGraphErrors* h_Asym[4];
  for (Int_t i = 0; i < 4; i++){
    h_Asym[i] = new TGraphErrors(Num , Del, Asym[i], 0, AsymErr[i]);
    h_Asym[i]->SetTitle(Form("Asymmetry[%d]",i));
  }
      
  //Plot Histograms
  c_Ee->cd(1); h_Ee[0]->Draw();  c_Ee->cd(2); h_Ee[1]->Draw(); c_Ee->cd(3); h_Ee[2]->Draw();
  c_Ee->cd(4); h_Ee[3]->Draw();  c_Ee->cd(5); h_Ee[4]->Draw(); c_Ee->cd(6); h_Ee[5]->Draw();
  
  c_Ep->cd(1); h_Ep[0]->Draw();  c_Ep->cd(2); h_Ep[1]->Draw(); c_Ep->cd(3); h_Ep[2]->Draw();
  c_Ep->cd(4); h_Ep[3]->Draw();  c_Ep->cd(5); h_Ep[4]->Draw(); c_Ep->cd(6); h_Ep[5]->Draw();
  
  c_EeEp->cd(1); h_EeEp[0]->Draw();  c_EeEp->cd(2); h_EeEp[1]->Draw();
  c_EeEp->cd(3); h_EeEp[2]->Draw();

  c_ThetaPhi->cd(1); h_ThetaPhiOR[0]->Draw(); c_ThetaPhi->cd(2), h_ThetaPhiOR[1]->Draw();
  c_ThetaPhi->cd(3); h_ThetaPhiAND[0]->Draw(); c_ThetaPhi->cd(4), h_ThetaPhiAND[1]->Draw();

  c_eXY->cd(1); h_eXY[0][0]->Draw(); c_eXY->cd(2); h_eXY[0][1]->Draw();
  c_eXY->cd(3); h_eXY[1][0]->Draw(); c_eXY->cd(4); h_eXY[1][1]->Draw();

  c_pXY->cd(1); h_pXY[0][0]->Draw(); c_pXY->cd(2); h_pXY[0][1]->Draw();
  c_pXY->cd(3); h_pXY[1][0]->Draw(); c_pXY->cd(4); h_pXY[1][1]->Draw();

  for(Int_t i = 0; i < 4; i++){c_Asym->cd(i+1);  h_Asym[i]->Draw();}
}
