//#define WH_cxx
#include "WH.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <iostream>
#include <string>
#include "Riostream.h"
#include <cmath>
#include <vector>

using namespace std;

// Note all below functions are used to plot data stored in VDC objects
// Since other data type are unchanged, they're plotted by PlotVDCKE.C
#define PI 3.1415926
const Double_t Wire_angle = 26.45*PI/180.;

//IConvert: Since Input directions now have been changed to (t_theta,t_phi)
//Need to convert (t_theta,t_phi) back to spherical angle (theta,phi)
//NOTE: need input in unit rad
void IConvert(Float_t& theta,Float_t&  phi){
	Float_t t = theta;
	Float_t p = phi;
	theta = atan(sqrt(pow(tan(t),2)+pow(tan(p),2)));
	phi = atan2(tan(p),tan(t));
	phi = phi<0 ? 2*PI+phi : phi;
}

// IConvertAngle: convert input projection angle (t_theta0, t_phi0) and (t_theta1, t_phi1) in Mrad
// to (theta0,phi0) and (theta1,phi1) in degrees
inline void IConvertAngle(Float_t& theta0, Float_t& phi0, Float_t& theta1, Float_t& phi1){
			theta0 /= 1000.;
			phi0 /= 1000.;
			theta1 /= 1000.;
			phi1 /= 1000.;

			IConvert(theta0, phi0);
			IConvert(theta1, phi1);

			// rad to degrees
			theta0 *= 180/PI;
			phi0 *= 180/PI;
			theta1 *= 180/PI;
			phi1 *= 180/PI;
}

inline void MradDeg(Double_t& angle){
	angle *= (180./(PI*1000.));
}

void WH::PlotAsym(){

	/*
	// Constraints , Bin size, etc.
	Double_t D_theta_tol = 1.;
	Double_t D_phi_tol = 1.;
	Double_t D_theta_bin = 0.5;
	//Double_t D_theta_bin = 15.;
	Double_t D_delta_bin = 0.25;
	*/

	Double_t D_theta_tol = 1.;
	Double_t D_phi_tol = 5.;
	//Double_t D_theta_bin = 10;
	//Double_t D_theta_bin = 2;
	//Double_t D_theta_bin = 0.5;
	//Double_t D_delta_bin = 2;

	Double_t D_theta_bin = 1;
	Double_t D_delta_bin = 0.5;

	//Int_t NPts = 50;
	Int_t NPts = 41;
	//Int_t NPts = 1001;
	//Double_t theta_bin[NPts];
	Double_t theta_bin;
	cerr << "Which theta for asymmetry? (4~15 degrees)" << endl;
	cin >> theta_bin;
	Double_t delta_bin[NPts];
	Double_t NPos[NPts];
	Double_t NNeg[NPts];
	Double_t Asym[NPts];
	//Uncertainty of Asymmetry
	Double_t Sig_A[NPts];
	Double_t Sig_Del[NPts];
	// Initialize:
	for(Int_t i = 0; i < NPts; i++){
		NPos[i] = 0;
		NNeg[i] = 0;
		Asym[i] = 0;
		Sig_A[i] = 0;
		Sig_Del[i] = 0;
	}
	for(Int_t i = 0; i < NPts; i++){
		//delta_bin[i] = 24.*((double)i+1.)/((double)NPts);
		// bins centered at delta = {-20,-19,...,20} (MeV)
		delta_bin[i] = 40.*((double)i/(double)(NPts-1))-20.;
	}

	Double_t Del_Phi = 50.;
	//TH1F* HPhi = new TH1F("HPhi", "Phi distribution", 200, 100., 300.);

	/*
	TH1F* HTheta = new TH1F("HTheta", "Theta distribution", 150, 0., 15.);
	TH1F* HPhi = new TH1F("HPhi", "Phi distribution", 50, 180.-Del_Phi , 180.+Del_Phi);
	TH1F* HDelta = new TH1F("HDelta", "Delta distribution", 500, -25, 25.);
	*/

	TH1F* HTheta = new TH1F("HTheta", "Theta distribution", 200, 0., 20.);
	TH1F* HPhi = new TH1F("HPhi", "Phi distribution", 400, 180.-Del_Phi , 180.+Del_Phi);
	TH1F* HDelta = new TH1F("HDelta", "Delta distribution", 500, -30, 30.);

	//TH1F* HTheta = new TH1F("HTheta", "Theta distribution", 100, theta_bin-D_theta_bin/2., theta_bin+D_theta_bin/2.);
	//TH1F* HPhi = new TH1F("HPhi", "Phi distribution", 200, 100., 300.);
	//TH1F* HPhi = new TH1F("HPhi", "Phi distribution", 50, 175., 185.);
	//TH1F* HDelta = new TH1F("HDelta", "Delta distribution", 100, delta_bin[(NPts-1)/2+6]-D_delta_bin/2.,delta_bin[(NPts-1)/2+6]+D_delta_bin/2.);
	
	TH2F* HBin = new TH2F("HBin","Bin distribution",
				100, delta_bin[(NPts-1)/2+6]-D_delta_bin/2., delta_bin[(NPts-1)/2+6]+D_delta_bin/2.,
				100, theta_bin-D_theta_bin/2., theta_bin+D_theta_bin/2.);

	// Check data with same bin in HBin:
	vector<vector<double> > BinData;
	vector<double> TData(7);
	long TIndex = 0;

	cerr << "Use Bin from Theta " << theta_bin-D_theta_bin/2. << " to " << theta_bin+D_theta_bin/2.;
	cerr << "Bin from Delta  " << delta_bin[(NPts-1)/2+6]-D_delta_bin/2. << " to " << delta_bin[(NPts-1)/2+6]+D_delta_bin/2. << endl;

	// Conversion constant from mrad to degrees
	Double_t DegMrad = (0.18/PI);

	if (fChain == 0) return;
   	Long64_t nentries = fChain->GetEntriesFast();

   	Long64_t nbytes = 0, nb = 0;
   	for (Long64_t jentry=0; jentry<nentries;jentry++) {
      		Long64_t ientry = LoadTree(jentry);
      		if (ientry < 0) break;
      		nb = fChain->GetEntry(jentry);   nbytes += nb;
      		// if (Cut(ientry) < 0) continue;
	
		TIndex++;

							IConvertAngle(I0_Theta,I0_Phi,I1_Theta,I1_Phi);
							if(D0_W_>0 && D1_W_>0){
								HTheta->Fill((I0_Theta+I1_Theta)/2.);
								//HPhi->Fill(180+abs(abs(I0_Phi-I1_Phi)-180));
								HPhi->Fill(abs(I0_Phi-I1_Phi));
								HDelta->Fill(I1_Energy-I0_Energy);
							}	

		if( D0_W_>0 && D1_W_>0 && D0_P_>0 && D1_P_> 0){

			//IConvertAngle(I0_Theta,I0_Phi,I1_Theta,I1_Phi);

			//cerr << "Input0 spherical angles (" << I0_Theta<< "," << I0_Phi << ") (degrees)" << endl;
			//cerr << "Input1 spherical angles (" << I1_Theta<< "," << I1_Phi << ") (degrees)" << endl;

			/*
			if(abs(I0_Theta-I1_Theta) < D_theta_tol && abs(abs(I0_Phi-I1_Phi)-180) < D_phi_tol){
			//if(abs(I0_Theta-I1_Theta) < D_theta_tol)
			//if(abs(abs(I0_Phi-I1_Phi)-180) < D_phi_tol)
				//cerr << "Input (Theta,Phi) difference: ( " << abs(I0_Theta-I1_Theta) << ", " << abs(abs(I0_Phi-I1_Phi)-180) <<") (degrees)" << endl;
				//cerr << "Theta = " << (I0_Theta+I1_Theta)/2. << endl;
			}
			*/

			//cerr << "Input Phi difference " << abs(I0_Phi-I1_Phi) << " (degrees)" << endl;

			//if(abs((I0_Theta-I1_Theta)*DegMrad) < D_theta_tol && abs((I0_Phi-I1_Phi)*DegMrad) < D_phi_tol){
			if(abs(I0_Theta-I1_Theta) < D_theta_tol && abs(abs(I0_Phi-I1_Phi)-180) < D_phi_tol){
				// Following code was to check issue on randomness of pair generation, and it's found that 
				// there're seeding issues. (By checking if a small bin with multiple data having all data the same)
				// Also note: it seems jentry and ientry cannot be accessed, reason?
				if( abs((I0_Theta+I1_Theta)/2.-6.325) < 0.005 && abs((I1_Energy-I0_Energy)-6.1125)<0.0025 ){
					TData[0] = I0_Theta;
					TData[1] = I1_Theta;
					TData[2] = I0_Phi;
					TData[3] = I1_Phi;
					TData[4] = I0_Energy;
					TData[5] = I1_Energy;
					//TData[7] = jentry;
					//TData[7] = ientry;
					TData[6] = TIndex;
					BinData.push_back(TData);
				}
				
				/*
				// Condition: abs(I0_Theta-I1_Theta) < D_theta_tol && abs(abs(I0_Phi-I1_Phi)-180) < D_phi_tol)
							HTheta->Fill((I0_Theta+I1_Theta)/2.);
							HPhi->Fill(abs(I0_Phi-I1_Phi));
							HDelta->Fill(I1_Energy-I0_Energy);
				*/

				//cerr << "Phi difference: " << (I0_Phi-I1_Phi)*DegMrad << "(Degrees)" << endl;
				for(int i = 0; i < NPts; i++){
					if(abs(theta_bin-(I0_Theta+I1_Theta)/2.) < (D_theta_bin/2.)
					&& abs(delta_bin[i]-(I1_Energy-I0_Energy)) < (D_delta_bin/2.)){ 
						NPos[i]++;
						//if(i == (NPts-1)/2+6) HBin->Fill((I0_Theta+I1_Theta)/2., (I1_Energy-I0_Energy));
						if(i == (NPts-1)/2+6) {
							//cerr << "i = " << i << endl;
							//HBin->Fill( 6,6 );
							//HBin->Fill( (I0_Theta+I1_Theta)/2., (I1_Energy-I0_Energy) );
							HBin->Fill( (I1_Energy-I0_Energy), (I0_Theta+I1_Theta)/2. );
							//cerr << "Theta bin = " << theta_bin << ", Theta = " << (I0_Theta+I1_Theta)/2. << endl;
							//cerr << "delta bin = " << delta_bin[i] << ", delta =  " << (I1_Energy-I0_Energy) << endl << endl;
						}
					}
					if(abs(theta_bin-(I0_Theta+I1_Theta)/2.) < (D_theta_bin/2.)
					&& abs(delta_bin[i]+(I1_Energy-I0_Energy)) < (D_delta_bin/2.)) NNeg[i]++;
				}
			}
   		}
	}

	// Compute Asymmetry and its standard deviation by error propagation
	for(int i = 0; i < NPts; i++) {
		if(NPos[i] == 0. && NNeg[i] == 0.) Asym[i] = 0;
		else{
			Asym[i] = (NPos[i]-NNeg[i])/(NPos[i]+NNeg[i]); 
			cerr << "Npos[i] = " << NPos[i] << ", Nneg[i] = " << NNeg[i] << endl;

			//Compute statistical uncertainties
			Sig_A[i] = 4*NPos[i]*NNeg[i]/(pow((NPos[i]+NNeg[i]),3.));
			Sig_A[i] = sqrt(Sig_A[i]);
		}
	}

	// Canvas
	TCanvas* CAsym = new TCanvas("CAsym", "Asymmetry graph", 1600,900);
	TCanvas* CBin =  new TCanvas("CBin","Bin distribution",1600,900);
	TCanvas* CTheta = new TCanvas("CTheta","Theta distribution",1600,900);
	TCanvas* CPhi= new TCanvas("CPhi","Phi distribution",1600,900);
	TCanvas* CDelta = new TCanvas("CDelta","Delta distribution",1600,900);
	//TGraph* GAsym = new TGraph(NPts, delta_bin, Asym);
	TGraphErrors* GAsym = new TGraphErrors(NPts, delta_bin, Asym, Sig_Del, Sig_A);
	GAsym->SetTitle(Form("Asymmetry for fixed theta = %f",theta_bin));
	GAsym->GetXaxis()->SetTitle("Delta");
	GAsym->GetYaxis()->SetTitle("Asymmetry");
	GAsym->GetXaxis()->CenterTitle();
	GAsym->GetYaxis()->CenterTitle();
	CAsym->cd(); GAsym->Draw();

	//HBin: Draw hit distribution in a (Theta,Delta) bin
	HBin->SetTitle(Form("Bin distribution for fixed theta = %f",theta_bin));
	HBin->GetXaxis()->SetTitle("Delta");
	HBin->GetYaxis()->SetTitle("Theta");
	HBin->GetXaxis()->CenterTitle();
	HBin->GetYaxis()->CenterTitle();
	CBin->cd(); HBin->Draw("COLZ");

	//Print Bin data
	for(unsigned i = 0; i < BinData.size(); i++){
		for(unsigned j = 0; j < TData.size() ; j++){
			cerr << " " << BinData[i][j];
		}
		cerr << endl;
	}

	CTheta->cd(); HTheta->Draw();
	CPhi->cd(); HPhi->Draw();
	CDelta->cd(); HDelta->Draw();
}

void WH::PlotInput(){
	//Canvas
	TCanvas* CIth1th2 = new TCanvas("CIth1th2","Input theta1 vs theta2",1600,900);
	TCanvas* CIDiff = new TCanvas("CIDiff","Input differences",800,900);
	CIDiff->Divide(1,2);
	TCanvas* CI[2];
	for(int i = 0; i < 2; i++){
		CI[i] = new TCanvas(Form("Input %d", i),Form("Input %d data with AND condition", i),1600,900);
		CI[i]->Divide(3,2);
	}
	//Histograms
	TH1F* HIEnergy[2];
	TH1F* HIDelta[2];
	TH1F* HIX[2];
	TH1F* HIY[2];
	TH1F* HITheta[2];
	TH1F* HIPhi[2];
	TH1F* HIThetaDiff = new TH1F("HIThetaDiff","Input Theta difference",220,-11.,11.);
	TH1F* HIPhiDiff = new TH1F("HIPhiDiff","Input Phi difference",3600,-360.,360);
	TH2F* HIth1th2 = new TH2F("HIth1th2", "Input theta_p vs theta_q", 150, 0., 15., 150, 0., 15.);
	HIth1th2->SetOption("COLZ");
	for(int i = 0; i < 2; i++){
		HIEnergy[i] = new TH1F(Form("HIEnergy %d",i), Form("Detector %d input energy",i),600,0, 60.0 );	
		HIDelta[i] = new TH1F(Form("HIDelta %d",i), Form("Detector %d input delta",i),600,-30.,30. );	
		HIX[i] = new TH1F(Form("HIX %d",i), Form("Detector %d input X",i),80,-4.,4. );	
		HIY[i] = new TH1F(Form("HIY %d",i), Form("Detector %d input Y",i),80,-4.,4. );	
		HITheta[i] = new TH1F(Form("HITheta %d",i), Form("Detetor %d input theta",i),300,0., 15.);	
		HIPhi[i] = new TH1F(Form("HIPhi %d",i), Form("Detector %d phi",i),3600,0., 360.);	
		//HIPhi[i] = new TH1F(Form("HIPhi %d",i), Form("Detector %d phi",i),3600,-180., 180.);	
	}
	Double_t DegMrad = (0.18/PI);
	//cerr << "1 Degree = " << DegMrad << "mrad." << endl;

   if (fChain == 0) return;
   Long64_t nentries = fChain->GetEntriesFast();

   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
	Long64_t ientry = LoadTree(jentry);
      	if (ientry < 0) break;
      	nb = fChain->GetEntry(jentry);   nbytes += nb;
      	// if (Cut(ientry) < 0) continue;

	IConvertAngle(I0_Theta,I0_Phi,I1_Theta,I1_Phi);

	if( D0_W_>0 && D1_W_>0 && D0_P_>0 && D1_P_> 0){
		for (int i = 0; i < 2; i++){
			if(i == 0){
				HIEnergy[i]->Fill(I0_Energy);
				HIDelta[i]->Fill(I0_Delta);
				HIX[i]->Fill(I0_X);
				HIY[i]->Fill(I0_Y);
				HITheta[i]->Fill(I0_Theta);
				HIPhi[i]->Fill(I0_Phi);
			}
			else{
				HIEnergy[i]->Fill(I1_Energy);
				HIDelta[i]->Fill(I1_Delta);
				HIX[i]->Fill(I1_X);
				HIY[i]->Fill(I1_Y);
				HITheta[i]->Fill(I1_Theta);
				HIPhi[i]->Fill(I1_Phi);
			}
		}
		HIThetaDiff->Fill((I1_Theta-I0_Theta));
		//Maybe not reflect real case
		//HIPhiDiff->Fill(abs(I1_Phi-I0_Phi));
		Double_t Temp_I0_Phi = I0_Phi;
		Double_t Temp_I1_Phi = I1_Phi;
		if(I0_Phi < 90) Temp_I0_Phi += 360;
		if(I1_Phi < 90) Temp_I1_Phi += 360;
		HIPhiDiff->Fill(Temp_I1_Phi-Temp_I0_Phi);
		HIth1th2->Fill(I0_Theta, I1_Theta);
   	}
   }
   
   // Plot histograms
   for(int i = 0; i < 2; i++){
	CI[i]->cd(1); HIEnergy[i]->Draw();
	CI[i]->cd(2); HIDelta[i]->Draw();
	CI[i]->cd(3); HIX[i]->Draw();
	CI[i]->cd(4); HIY[i]->Draw();
	CI[i]->cd(5); HITheta[i]->Draw();
	CI[i]->cd(6); HIPhi[i]->Draw();
   }
   CIDiff->cd(1); HIThetaDiff->Draw();
   CIDiff->cd(2); HIPhiDiff->Draw();
   CIth1th2->cd(); HIth1th2->Draw();
}

void WH::PlotMonitor(){
	//Canvas
	TCanvas* CMDiff = new TCanvas("CMDiff","Monitor differences",800,900);
	CMDiff->Divide(1,2);
	TCanvas* CM[2];
	for(int i = 0; i < 2; i++){
		CM[i] = new TCanvas(Form("Monitor %d", i),Form("Monitor %d data with AND condition", i),1600,900);
		CM[i]->Divide(3,2);
	}
	//Histograms
	TH1F* HMEnergy[2];
	TH1F* HMCharge[2];
	TH1F* HMX[2];
	TH1F* HMY[2];
	TH1F* HMTheta[2];
	TH1F* HMPhi[2];
	TH1F* HMThetaDiff = new TH1F("HMThetaDiff","Monitor Theta difference",220,-11.,11.);
	TH1F* HMPhiDiff = new TH1F("HMPhiDiff","Monitor Phi difference",3600,-360.,360);
	for(int i = 0; i < 2; i++){
		HMEnergy[i] = new TH1F(Form("HMEnergy %d",i), Form("Detector %d monitor energy",i),600,0, 60.0 );	
		HMCharge[i] = new TH1F(Form("HMCharge %d",i), Form("Detector %d monitor delta",i),600,-30.,30. );	
		HMX[i] = new TH1F(Form("HMX %d",i), Form("Detector %d monitor X",i),80,-4.,4. );	
		HMY[i] = new TH1F(Form("HMY %d",i), Form("Detector %d monitor Y",i),80,-4.,4. );	
		HMTheta[i] = new TH1F(Form("HMTheta %d",i), Form("Detetor %d monitor theta",i),300,0., 15.);	
		HMPhi[i] = new TH1F(Form("HMPhi %d",i), Form("Detector %d phi",i),3600,0., 360.);	
	}
	Double_t DegMrad = (0.18/PI);
	//cerr << "1 Degree = " << DegMrad << "mrad." << endl;

   if (fChain == 0) return;
   Long64_t nentries = fChain->GetEntriesFast();

   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
	Long64_t ientry = LoadTree(jentry);
      	if (ientry < 0) break;
      	nb = fChain->GetEntry(jentry);   nbytes += nb;
      	// if (Cut(ientry) < 0) continue;
	if( D0_W_>0 && D1_W_>0 && D0_P_>0 && D1_P_> 0){
		for (int i = 0; i < 2; i++){
			if(i == 0){
				HMEnergy[i]->Fill(M0_Energy);
				HMCharge[i]->Fill(M0_Charge);
				HMX[i]->Fill(M0_X);
				HMY[i]->Fill(M0_Y);
				HMTheta[i]->Fill(M0_Theta*DegMrad);
				HMPhi[i]->Fill(M0_Phi*DegMrad);
			}
			else{
				HMEnergy[i]->Fill(M1_Energy);
				HMCharge[i]->Fill(M1_Charge);
				HMX[i]->Fill(M1_X);
				HMY[i]->Fill(M1_Y);
				HMTheta[i]->Fill(M1_Theta*DegMrad);
				HMPhi[i]->Fill(M1_Phi*DegMrad);
			}
		}
		HMThetaDiff->Fill((M1_Theta-M0_Theta)*DegMrad);
		HMPhiDiff->Fill((M1_Phi-M0_Theta)*DegMrad);
   	}
   }
   
   // Plot histograms
   for(int i = 0; i < 2; i++){
	CM[i]->cd(1); HMEnergy[i]->Draw();
	CM[i]->cd(2); HMCharge[i]->Draw();
	CM[i]->cd(3); HMX[i]->Draw();
	CM[i]->cd(4); HMY[i]->Draw();
	CM[i]->cd(5); HMTheta[i]->Draw();
	CM[i]->cd(6); HMPhi[i]->Draw();
   }
   CMDiff->cd(1); HMThetaDiff->Draw();
   CMDiff->cd(2); HMPhiDiff->Draw();
}


void WH::PlotPaddle(){
	//Canvas
	TCanvas* CP[2];
	for(int i = 0; i < 2; i++){
		CP[i] = new TCanvas(Form("Paddle %d", i),Form("Paddle %d data with AND condition", i),1600,900);
		CP[i]->Divide(2,2);
	}
	//Histograms
	TH1I* HPNum[2];
	TH1F* HPEdep[2];
	TH1F* HPLight[2];
	TH1F* HPTime[2];
	for(int i = 0; i < 2; i++){
		HPNum[i] = new TH1I(Form("HPNum %d",i), Form("Detector %d Paddle Number",i),29,0,29);
		HPEdep[i] = new TH1F(Form("HPEdep %d",i), Form("Detector %d Paddle Energy deposite",i),100,0.,10.);
		HPLight[i] = new TH1F(Form("HPLight %d",i), Form("Detctor %d Paddle Light", i),100,0.,10.);
		HPTime[i] = new TH1F(Form("HPTime %d",i), Form("Detector %d Paddle Time",i),100,5,9.);
		//HPTime[i]->SetBit(TH1::kCanRebin);
	}

	if (fChain == 0) return;
   	Long64_t nentries = fChain->GetEntriesFast();

   	Long64_t nbytes = 0, nb = 0;
   	for (Long64_t jentry=0; jentry<nentries;jentry++) {
		Long64_t ientry = LoadTree(jentry);
      		if (ientry < 0) break;
      		nb = fChain->GetEntry(jentry);   nbytes += nb;
      		// if (Cut(ientry) < 0) continue;
		// Add AND condition for plotting
		if( D0_W_>0 && D1_W_>0 && D0_P_>0 && D1_P_> 0){
			for(Int_t i = 0; i < D0_P_; i++){
				HPNum[0]->Fill(D0_P_PNum[i]);
				HPEdep[0]->Fill(D0_P_Edep[i]);
				HPLight[0]->Fill(D0_P_Light[i]);
				HPTime[0]->Fill(D0_P_Time[i]);
			}
			for(Int_t i = 0; i < D1_P_; i++){
				HPNum[1]->Fill(D1_P_PNum[i]);
				HPEdep[1]->Fill(D1_P_Edep[i]);
				HPLight[1]->Fill(D1_P_Light[i]);
				HPTime[1]->Fill(D1_P_Time[i]);
			}
		}
   	}
	//Plot histograms
	for(int i = 0; i < 2; i++){
		CP[i]->cd(1); HPNum[i]->Draw();
		CP[i]->cd(2); HPEdep[i]->Draw();
		CP[i]->cd(3); HPLight[i]->Draw();
		CP[i]->cd(4); HPTime[i]->Draw();
	}
}
/*
// Below is ROOT Make Class Loop() example for copy
void WH::Loop(){
   if (fChain == 0) return;
   Long64_t nentries = fChain->GetEntriesFast();

   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;
   }
}
*/
//Loop: Only hit with MAX KE in each wireplane with AND condition  associated data stored.
void WH::PlotVDC()
{
//   In a ROOT session, you can do:
//      Root > .L WH.C
//      Root > WH t
//      Root > t.GetEntry(12); // Fill t data members with entry number 12
//      Root > t.Show();       // Show values of entry 12
//      Root > t.Show(16);     // Read and show values of entry 16
//      Root > t.Loop();       // Loop on all entries
//

//     This is the loop skeleton where:
//    jentry is the global entry number in the chain
//    ientry is the entry number in the current Tree
//  Note that the argument to GetEntry must be:
//    jentry for TChain::GetEntry
//    ientry for TTree::GetEntry and TBranch::GetEntry
//
//       To read only selected branches, Insert statements like:
// METHOD1:
//    fChain->SetBranchStatus("*",0);  // disable all branches
//    fChain->SetBranchStatus("branchname",1);  // activate branchname
// METHOD2: replace line
//    fChain->GetEntry(jentry);       //read all branches
//by  b_branchname->GetEntry(ientry); //read only this branch
   if (fChain == 0) return;

      // Fill all histograms with different conditions
      bool AND, MAXKE;
      string SAND, SMAXKE;
      cerr << "AND condition ? (Y/N)" << endl;
      //getline(cin, SAND);	
      cin >> SAND;
      cerr << "MAXKE ? (Y/N) "<< endl;
      //getline(cin, SMAXKE);
      cin >> SMAXKE;

      //cerr << "AND is " << SAND << ", MAXKE is " << SMAXKE << endl;

      if(SAND == "Y") AND = true;
      else AND  = false;
      if(SMAXKE == "Y") MAXKE = true;
      else MAXKE  = false;
   
   //Canvases
   TCanvas* CENum = new TCanvas("CENum", "Event Number", 800, 400);
   TCanvas* CD0Wire = new TCanvas("CD0Wire", "Detector0 Wire Number", 1600, 1200);
   TCanvas* CD0KE = new TCanvas("CD0KE", "Detector0 VDC WirePlane KE", 1600, 1200);
   TCanvas* CD0X = new TCanvas("CD0X", "Detector0 VDC WirePlane X",1600, 1200);
   TCanvas* CD0Y = new TCanvas("CD0Y", "Detector0 VDC WirePlane Y",1600, 1200);
   TCanvas* CD0ToF = new TCanvas("CD0ToF", "Detector0 VDC WirePlane ToF",1600, 1200);
   TCanvas* CD0Charge = new TCanvas("CD0Charge", "Detector0 VDC WirePlane Charge",1600, 1200);
   TCanvas* CD0Particle = new TCanvas("CD0Particle", "Detector0 VDC WirePlane Particle Type",1600, 1200);

   CD0Wire->Divide(2,2);
   CD0KE->Divide(2,2);
   CD0X->Divide(2,2);
   CD0Y->Divide(2,2);
   CD0ToF->Divide(2,2);
   CD0Charge->Divide(2,2);
   CD0Particle->Divide(2,2);

   TCanvas* CD1Wire = new TCanvas("CD1Wire", "Detector1 Wire Number", 1600, 1200);
   TCanvas* CD1KE = new TCanvas("CD1KE", "Detector1 VDC WirePlane KE", 1600, 1200);
   TCanvas* CD1X = new TCanvas("CD1X", "Detector1 VDC WirePlane X",1600, 1200);
   TCanvas* CD1Y = new TCanvas("CD1Y", "Detector1 VDC WirePlane Y",1600, 1200);
   TCanvas* CD1ToF = new TCanvas("CD1ToF", "Detector1 VDC WirePlane ToF",1600, 1200);
   TCanvas* CD1Charge = new TCanvas("CD1Charge", "Detector1 VDC WirePlane Charge",1600, 1200);
   TCanvas* CD1Particle = new TCanvas("CD1Particle", "Detector1 VDC WirePlane Particle Type",1600, 1200);

   CD1Wire->Divide(2,2);
   CD1KE->Divide(2,2);
   CD1X->Divide(2,2);
   CD1Y->Divide(2,2);
   CD1ToF->Divide(2,2);
   CD1Charge->Divide(2,2);
   CD1Particle->Divide(2,2);   

   //To find out abnormality of 2nd wireplanes of each wire chamber, plot (X,Y) distribution of each wire plane
   TCanvas* CD0XY = new TCanvas("CD0XY", "Detector0 (X,Y)", 1600, 1200);
   CD0XY->Divide(2,2);
   TCanvas* CD1XY = new TCanvas("CD1XY", "Detector1 (X,Y)", 1600, 1200);
   CD1XY->Divide(2,2);
   //Also plot hit coordinate parallel to wire direction
   TCanvas* CD0Xp = new TCanvas("CD0Xp", "Detector0 Xp", 1600, 1200);
   CD0Xp -> Divide(2,2);
   TCanvas* CD1Xp = new TCanvas("CD1Xp", "Detector1 Xp", 1600, 1200);
   CD1Xp -> Divide(2,2);

   //Histograms
   TH1D* HENum = new TH1D("HENum", "Event Number", 100, 0, 100000);
   HENum->SetBit(TH1::kCanRebin);

   TH1I* H0Wire[4];
   TH1I* H1Wire[4];
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

   // Additional (X,Y) and Xp check
   TH2F* H0XY[4];
   TH2F* H1XY[4];
   TH1F* H0Xp[4];
   TH1F* H1Xp[4];
   
   for (int i = 0; i< 4; i++){
     H0Wire[i] = new TH1I(Form("H0Wire WirePlane%d", i), Form("Detector0 VDC WirePlane%d Wire Number", i), 279, 0, 279);
     H0KE[i] = new TH1F(Form("H0KE WirePlane%d", i) , Form("Detector0 VDC WirePlane%d KE", i), 600, 0,60.0);


     H0X[i] = new TH1F(Form("H0X WirePlane%d", i) , Form("Detector0 VDC WirePlane%d X", i), 4000, -50., 50.);
     //  H0X[i]->SetBit(TH1::kCanRebin);
     H0Y[i] = new TH1F(Form("H0Y WirePlane%d", i) , Form("Detector0 VDC WirePlane%d Y", i), 1000, -25., 25.);
     //  H0Y[i]->SetBit(TH1::kCanRebin);


     H0ToF[i] = new TH1F(Form("H0ToF WirePlane%d", i), Form("Detector0 VDC WirePlane%d ToF", i), 100, 0., 10.);
     H0Charge[i] = new TH1F(Form("H0Charge WirePlane%d", i), Form("Detector0 VDC WirePlane%d Charge", i), 4, -2., 2.);
     H0Particle[i] = new TH1F(Form("H0Particle WirePlane%d", i), Form("Detector0 VDC WirePlane%d Particle",i), 50, -25., 25.);
     //  H0Particle[i]->SetBit(TH1::kCanRebin);

     H1Wire[i] = new TH1I(Form("H1Wire WirePlane%d", i), Form("Detector1 VDC WirePlane%d Wire Number", i), 279, 0, 279);
     H1KE[i] = new TH1F(Form("H1KE WirePlane%d", i) , Form("Detector1 VDC WirePlane%d KE", i), 600, 0,60.0);


     H1X[i] = new TH1F(Form("H1X WirePlane%d", i) , Form("Detector1 VDC WirePlane%d X", i), 4000, -50., 50.);
     //  H1X[i]->SetBit(TH1::kCanRebin);
     H1Y[i] = new TH1F(Form("H1Y WirePlane%d", i) , Form("Detector1 VDC WirePlane%d Y", i), 1000, -25., 25.);
     //  H1Y[i]->SetBit(TH1::kCanRebin);


     H1ToF[i] = new TH1F(Form("H1ToF WirePlane%d", i), Form("Detector1 VDC WirePlane%d ToF", i), 100, 0., 10.);
     H1Charge[i] = new TH1F(Form("H1Charge WirePlane%d", i), Form("Detector1 VDC WirePlane%d Charge", i), 4, -2., 2.);
     H1Particle[i] = new TH1F(Form("H1Particle WirePlane%d", i), Form("Detector1 VDC WirePlane%d Particle",i), 50, -25., 25.);
     //  H1Particle[i]->SetBit(TH1::kCanRebin);

	//(X,Y) distribution
	H0XY[i] = new TH2F(Form("H0 XY WirePlane %d",i), Form("Detector0 VDC WirePlane%d (X,Y)",i), 500, -50., 50., 250, -25., 25.);
	H1XY[i] = new TH2F(Form("H1 XY WirePlane %d",i), Form("Detector1 VDC WirePlane%d (X,Y)",i), 500, -50., 50., 250, -25., 25.);

	H0Xp[i] = new TH1F(Form("H0 Xp WirePlane %d",i), Form("Detector0 VDC WirePlane%d Xp",i), 5000, -50, 50);
	H1Xp[i] = new TH1F(Form("H1 Xp WirePlane %d",i), Form("Detector1 VDC WirePlane%d Xp",i), 5000, -50, 50);
   }

	Double_t Xp0[D0_W_];
	Double_t Xp1[D1_W_];
	for(int i = 0; i < D0_W_; i++) Xp0[i] = 0.;
	for(int i = 0; i < D1_W_; i++) Xp1[i] = 0.;
	Double_t Rp[2] = {tan(Wire_angle), -1.};
	for(int i = 0; i < 2; i++) { Rp[i] /= sqrt(1+pow(tan(Wire_angle),2.)); }
	//cerr << "Rp = (" << Rp[0] << "," << Rp[1] << ")" << "   ,   |Rp|" <<  sqrt(1+pow(tan(Wire_angle),2.)) << endl;

   Long64_t nentries = fChain->GetEntriesFast();

   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;

      //Fill CENum:
      HENum->Fill(ENum);
      
	// Loop over data from a single hit, find MAX KE in each wireplane and its index
	// flag0 flag1 to make histograms fill with MAX KE condition
	int flag0[4]={0,0,0,0};
	int flag1[4]={0,0,0,0};
	int Index0[4]={0,0,0,0};
      	double TEMPKE0[4]={0.,0.,0.,0.};     
      	int Index1[4]={0,0,0,0};
      	double TEMPKE1[4]={0.,0.,0.,0.};     
	if(MAXKE == true){
      		for(int n = 0; n < D0_W_; n++){
			for(int m = 0; m < 4; m++){
		  	if(D0_W_WirePlane[n] == m && TEMPKE0[m] < D0_W_KE[n]){
				TEMPKE0[m] = D0_W_KE[n];
		    		Index0[m] = n;
				flag0[m] = 1;
		}}}

		// Do the same for D1:
      		for(int n = 0; n < D1_W_; n++){
			for(int m = 0; m < 4; m++){
		  	if(D1_W_WirePlane[n] == m && TEMPKE1[m] < D1_W_KE[n]){
	   	 		TEMPKE1[m] = D1_W_KE[n];
	           		Index1[m] = n;
				flag1[m] = 1;
		}}}
	}

	// Use MAX KE in each wireplane associate data and "AND" condition cuts to fill the histograms
      	// "AND condition"
      	if(AND == true && MAXKE == true){
		if( D0_W_>0 && D1_W_>0 && D0_P_>0 && D1_P_> 0){
			for (int i = 0; i < 4; i++){
				if(flag0[i] == 1){
					H0Wire[i]->Fill(D0_W_WireNum[Index0[i]]);
				  	H0KE[i]->Fill(D0_W_KE[Index0[i]]);
				  	H0X[i]->Fill(D0_W_X[Index0[i]]);
				  	H0Y[i]->Fill(D0_W_Y[Index0[i]]);
				  	H0ToF[i]->Fill(D0_W_ToF[Index0[i]]);
				  	H0Charge[i]->Fill(D0_W_Charge[Index0[i]]);
				  	H0Particle[i]->Fill(D0_W_Particle[Index0[i]]);

					H0XY[i]->Fill(D0_W_X[Index0[i]], D0_W_Y[Index0[i]]);
					H0Xp[i]->Fill( D0_W_X[Index0[i]]*Rp[0] + D0_W_Y[Index0[i]]*Rp[1] );
	  			}
				if(flag1[i] == 1){
				  	H1Wire[i]->Fill(D1_W_WireNum[Index1[i]]);
				  	H1KE[i]->Fill(D1_W_KE[Index1[i]]);
					H1X[i]->Fill(D1_W_X[Index1[i]]);
				  	H1Y[i]->Fill(D1_W_Y[Index1[i]]);
				  	H1ToF[i]->Fill(D1_W_ToF[Index1[i]]);
				  	H1Charge[i]->Fill(D1_W_Charge[Index1[i]]);
				  	H1Particle[i]->Fill(D1_W_Particle[Index1[i]]);

					H1XY[i]->Fill(D1_W_X[Index1[i]], D1_W_Y[Index1[i]]);
					H1Xp[i]->Fill( D1_W_X[Index1[i]]*Rp[0] + D1_W_Y[Index1[i]]*Rp[1] );
				}
			}
		}
       	}
       	else if (AND == true && MAXKE == false){
		// Use "AND" condition cuts to fill the histograms
      		if( D0_W_>0 && D1_W_>0 && D0_P_>0 && D1_P_> 0){
			// Fill all data with "AND" condition in D0
			for (int i = 0; i < D0_W_; i++){
				for (int j = 0; j < 4; j++){
					if(D0_W_WirePlane[i] == j){
						H0Wire[j]->Fill(D0_W_WireNum[i]);
					      	H0KE[j]->Fill(D0_W_KE[i]);
					      	H0X[j]->Fill(D0_W_X[i]);
					  	H0Y[j]->Fill(D0_W_Y[i]);
					      	H0ToF[j]->Fill(D0_W_ToF[i]);
						H0Charge[j]->Fill(D0_W_Charge[i]);
					      	H0Particle[j]->Fill(D0_W_Particle[i]);

						H0XY[j]->Fill(D0_W_X[i], D0_W_Y[i]);
						H0Xp[i]->Fill( D0_W_X[i]*Rp[0] + D0_W_Y[i]*Rp[1] );
			}}}
				// Same for D1
			for(int i = 0; i < D1_W_; i++){
				for(int j = 0; j < 4; j++){
					if(D1_W_WirePlane[i] == j){
						H1Wire[j]->Fill(D1_W_WireNum[i]);
			 		     	H1KE[j]->Fill(D1_W_KE[i]);
		     	   		     	H1X[j]->Fill(D1_W_X[i]);
		           		     	H1Y[j]->Fill(D1_W_Y[i]);
		     	     		     	H1ToF[j]->Fill(D1_W_ToF[i]);
		      	      		     	H1Charge[j]->Fill(D1_W_Charge[i]);
			     	 	      	H1Particle[j]->Fill(D1_W_Particle[i]);

						H1XY[j]->Fill(D1_W_X[i], D1_W_Y[i]);
						H1Xp[i]->Fill( D1_W_X[i]*Rp[0] + D1_W_Y[i]*Rp[1] );
			}}}
      		}
   	}
	else if (AND == false && MAXKE == true){
		for (int i = 0; i < 4; i++){
			if(flag0[i] == 1){
				H0Wire[i]->Fill(D0_W_WireNum[Index0[i]]);
			  	H0KE[i]->Fill(D0_W_KE[Index0[i]]);
			  	H0X[i]->Fill(D0_W_X[Index0[i]]);
			  	H0Y[i]->Fill(D0_W_Y[Index0[i]]);
			  	H0ToF[i]->Fill(D0_W_ToF[Index0[i]]);
			  	H0Charge[i]->Fill(D0_W_Charge[Index0[i]]);
			  	H0Particle[i]->Fill(D0_W_Particle[Index0[i]]);

					H0XY[i]->Fill(D0_W_X[Index0[i]], D0_W_Y[Index0[i]]);
					H0Xp[i]->Fill( D0_W_X[Index0[i]]*Rp[0] + D0_W_Y[Index0[i]]*Rp[1] );
  			}
			if(flag1[i] == 1){
			  	H1Wire[i]->Fill(D1_W_WireNum[Index1[i]]);
			  	H1KE[i]->Fill(D1_W_KE[Index1[i]]);
				H1X[i]->Fill(D1_W_X[Index1[i]]);
			  	H1Y[i]->Fill(D1_W_Y[Index1[i]]);
			  	H1ToF[i]->Fill(D1_W_ToF[Index1[i]]);
			  	H1Charge[i]->Fill(D1_W_Charge[Index1[i]]);
			  	H1Particle[i]->Fill(D1_W_Particle[Index1[i]]);

					H1XY[i]->Fill(D1_W_X[Index1[i]], D1_W_Y[Index1[i]]);
					H1Xp[i]->Fill( D1_W_X[Index1[i]]*Rp[0] + D1_W_Y[Index1[i]]*Rp[1] );
			}
		}
	}
	else if (AND == false && MAXKE == false){
			for (int i = 0; i < D0_W_; i++){
				for (int j = 0; j < 4; j++){
					if(D0_W_WirePlane[i] == j){
						H0Wire[j]->Fill(D0_W_WireNum[i]);
					      	H0KE[j]->Fill(D0_W_KE[i]);
					      	H0X[j]->Fill(D0_W_X[i]);
					  	H0Y[j]->Fill(D0_W_Y[i]);
					      	H0ToF[j]->Fill(D0_W_ToF[i]);
						H0Charge[j]->Fill(D0_W_Charge[i]);
					      	H0Particle[j]->Fill(D0_W_Particle[i]);

						H0XY[j]->Fill(D0_W_X[i], D0_W_Y[i]);
						//cerr << "D0_W_X " << D0_W_X[i] << ", D0_W_Y " << D0_W_Y[i] << endl;
						//cerr << "Rp[0] " << Rp[0] << ", Rp[1] " << Rp[1] << endl;
						H0Xp[i]->Fill( D0_W_X[i]*Rp[0] + D0_W_Y[i]*Rp[1] );
}}}
				// Same for D1
			for(int i = 0; i < D1_W_; i++){
				for(int j = 0; j < 4; j++){
					if(D1_W_WirePlane[i] == j){
						H1Wire[j]->Fill(D1_W_WireNum[i]);
			 		     	H1KE[j]->Fill(D1_W_KE[i]);
		     	   		     	H1X[j]->Fill(D1_W_X[i]);
		           		     	H1Y[j]->Fill(D1_W_Y[i]);
		     	     		     	H1ToF[j]->Fill(D1_W_ToF[i]);
		      	      		     	H1Charge[j]->Fill(D1_W_Charge[i]);
			     	 	      	H1Particle[j]->Fill(D1_W_Particle[i]);

						H1XY[j]->Fill(D1_W_X[i], D1_W_Y[i]);
						H1Xp[i]->Fill( D1_W_X[i]*Rp[0] + D1_W_Y[i]*Rp[1] );
}}}
	}
}
   
   //Plot those histograms: these histograms shows the highest KE in each wire plane
   // and associate X,Y,ToF,Charge,Particle
   // Also note, only "AND" condition data are drawn.
   CENum->cd(); HENum->Draw();
   for(int i = 0; i < 4; i++){
     CD0Wire->cd(i+1);  H0Wire[i]->Draw();
     CD0KE->cd(i+1);    H0KE[i]->Draw();
     CD0X->cd(i+1);     H0X[i]->Draw();
     CD0Y->cd(i+1);     H0Y[i]->Draw();
     CD0ToF->cd(i+1);   H0ToF[i]->Draw();
     CD0Charge->cd(i+1);H0Charge[i]->Draw();
     CD0Particle->cd(i+1);H0Particle[i]->Draw();
	CD0XY->cd(i+1); H0XY[i]->Draw("COLZ");
	CD0Xp->cd(i+1); H0Xp[i]->Draw();

     CD1Wire->cd(i+1);  H1Wire[i]->Draw();
     CD1KE->cd(i+1);    H1KE[i]->Draw();
     CD1X->cd(i+1);     H1X[i]->Draw();
     CD1Y->cd(i+1);     H1Y[i]->Draw();
     CD1ToF->cd(i+1);   H1ToF[i]->Draw();
     CD1Charge->cd(i+1);H1Charge[i]->Draw(); 
     CD1Particle->cd(i+1);H1Particle[i]->Draw();
	CD1XY->cd(i+1); H1XY[i]->Draw("COLZ");
	CD1Xp->cd(i+1); H1Xp[i]->Draw();
     }
}


// All below are integrated in PlotVDC()

// Compute and plot asymmetry curves using MAX KE ("OR" condition)  in each wireplane
// So there are 4 MAX KE computed asymmetry curve plots
// Also WhreNum vs X, KE vs X.
void WH::PlotVDCMAXKE(){
   if (fChain == 0) return;
   // "AND" flag
   bool AND;
   string SAND;
   cerr << "Apply AND condition? (y/n)" << endl;
   cin >> SAND;
   if(SAND == "y") AND = true;
   else if(SAND == "n") AND = false;
   else { cerr << "Incorrect input AND : " << SAND << endl;  return; }

   //Canvas
   TCanvas* CDN = new TCanvas("CDN", "Number counts in (Ee-Ep)", 1600, 1200);
   CDN -> Divide(2,2);
   TCanvas* CDA = new TCanvas("CDA", "Total Asymmetry ", 1600, 1200);
   CDA -> Divide(2,2);
   //WireNum vs X
   TCanvas* C0WX = new TCanvas("C0WX", "Detector0 WireNum vs X", 1600, 1200);
   C0WX -> Divide(2,2);
   TCanvas* C1WX = new TCanvas("C1WX", "Detector1 WireNum vs X", 1600, 1200);
   C1WX -> Divide(2,2);
   //KE vs X
   TCanvas* C0EX = new TCanvas("C0EX","Detector0 KE vs X", 1600, 1200);
   C0EX -> Divide(2,2);
   TCanvas* C1EX = new TCanvas("C1EX","Detector1 KE vs X", 1600, 1200);
   C1EX -> Divide(2,2);

   //Histograms:
   TH2F* H0WX[4];
   TH2F* H1WX[4];
   TH2F* H0EX[4];
   TH2F* H1EX[4];

   for (int i = 0; i< 4; i++){
     H0WX[i] = new TH2F(Form("H0 WireNum vs X WirePlane%d", i), Form("Detector0 VDC WirePlane%d Wire vs X", i), 200, -100., 100., 279, 0., 279.);
     H0EX[i] = new TH2F(Form("H0 KE vs X WirePlane%d", i), Form("Detector0 VDC WirePlane%d KE vs X", i), 200, -100., 100., 240, 0., 60.);
     H0WX[i] -> SetOption("COLZ");
     H0EX[i] -> SetOption("COLZ");

     H1WX[i] = new TH2F(Form("H1 WireNum vs X WirePlane%d", i), Form("Detector1 VDC WirePlane%d Wire vs X", i), 200, -100., 100., 279, 0., 279.);
     H1EX[i] = new TH2F(Form("H1 KE vs X WirePlane%d", i), Form("Detector1 VDC WirePlane%d KE vs X", i), 200, -100., 100., 240, 0., 60.);
     H1WX[i] -> SetOption("COLZ");
     H1EX[i] -> SetOption("COLZ");
   }
  
   // Data for Graphs
   Int_t NPt = 60;
   Double_t del = 120./(Double_t)NPt;
   Int_t N[4][NPt];
   Double_t A[4][NPt],  Del[NPt];

   for(int n = 0 ; n < 4; n++)
     for(int m = 0; m < NPt; m++){
       N[n][m] = 0; Del[m] = -30.0 + m*del/2. + del/4.; A[n][m] = 0.;}
 
   Long64_t nentries = fChain->GetEntriesFast();

   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;

       // Loop over data from a single hit, find MAX KE in each wireplane and its index
      int flag0[4] = {0,0,0,0}; 
      int Index0[4]={0,0,0,0};
      double TEMPKE0[4]={0.,0.,0.,0.};
      double TEMPX0[4]= {0.,0.,0.,0.};
      for(int n = 0; n < D0_W_; n++){
	for(int m = 0; m < 4; m++){
	  if(D0_W_WirePlane[n] == m && TEMPKE0[m] < D0_W_KE[n]){
	    TEMPKE0[m] = D0_W_KE[n];
	    TEMPX0[m] = D0_W_X[n];
	    Index0[m] = n;
	    flag0[m] = 1;
	  }}}
/*
      cerr << "Event number " << jentry << endl;
      for(int i = 0; i < 4; i++){
      cerr << "D0 WirePlane " << i << ", KE " << TEMPKE0[i] << ", X " << TEMPX0[i] << endl << endl;}
*/
      // Do the same for D1:
      int flag1[4] = {0,0,0,0}; 
      int Index1[4]={0,0,0,0};
      double TEMPKE1[4]={0.,0.,0.,0.};
      double TEMPX1[4]={0.,0.,0.,0.};
      for(int n = 0; n < D1_W_; n++){
	for(int m = 0; m < 4; m++){
	  if(D1_W_WirePlane[n] == m && TEMPKE1[m] < D1_W_KE[n]){
	    TEMPKE1[m] = D1_W_KE[n];
	    TEMPX1[m] = D1_W_X[n];
	    Index1[m] = n;
	    flag1[m] = 1;
	  }}}

      // Find out number of Ee-Ep counts for each wireplane MAX KE data, 
      // also need MAX KE on both side exist and both larger than 2 MeV
      for(int n = 0; n < 4; n++){
	if(TEMPKE0[n] > 2.0 && TEMPKE1[n] > 2.0){
	  for(int m = 0; m < NPt; m++){
	    if(TEMPKE1[n] - TEMPKE0[n] > -60.0 + (double)m*del && TEMPKE1[n]- TEMPKE0[n] < -60.0 + ((double)m+1.0)*del)
	      N[n][m]++;}}}
 
      //Fill Histograms with MAX KE on each WirePlane condition
	if(AND == false){
      		for(int i = 0; i < 4; i++){
	 		if(flag0[i]==1) H0WX[i]->Fill(D0_W_X[Index0[i]], D0_W_WireNum[Index0[i]]);
	  		if(flag1[i]==1) H1WX[i]->Fill(D1_W_X[Index1[i]], D1_W_WireNum[Index1[i]]);
	  		if(flag0[i]==1) H0EX[i]->Fill(D0_W_X[Index0[i]], D0_W_KE[Index0[i]]);
	  		if(flag1[i]==1) H1EX[i]->Fill(D1_W_X[Index1[i]], D1_W_KE[Index1[i]]);	
      		}
	}
	else if(AND == true){
		if( D0_W_>0 && D1_W_>0 && D0_P_>0 && D1_P_> 0){
      			for(int i = 0; i < 4; i++){
	 			if(flag0[i]==1) H0WX[i]->Fill(D0_W_X[Index0[i]], D0_W_WireNum[Index0[i]]);
	  			if(flag1[i]==1) H1WX[i]->Fill(D1_W_X[Index1[i]], D1_W_WireNum[Index1[i]]);
	  			if(flag0[i]==1) H0EX[i]->Fill(D0_W_X[Index0[i]], D0_W_KE[Index0[i]]);
	  			if(flag1[i]==1) H1EX[i]->Fill(D1_W_X[Index1[i]], D1_W_KE[Index1[i]]);	
      			}
		}	
	}
}
   for(int i = 0; i < NPt; i++){
     cerr << "N[0][" << i << "] = " << N[0][i] << endl;}
   
   //Number of counts graph and asymmetry graph
   TGraph* GDN[4];
   TGraph* GDA[4];
   
   Double_t ND[4][NPt];
   for(Int_t i = 0; i < 4; i++)
     for(Int_t j = 0; j < NPt; j++){
       ND[i][j] = (Double_t)N[i][j];}

   for(int i = 0; i < NPt; i++){
     cerr << "ND[0][" << i << "] = " << ND[0][i] << endl;}
   
   for(Int_t i = 0; i < 4; i++){
     GDN[i] = new TGraph (NPt, Del, ND[i]);
     for(Int_t j = 0; j < NPt/2; j++){
       if(N[i][NPt/2+j] == 0 && N[i][NPt/2-j-1] == 0) A[i][NPt/2+j] = 0.;
       else A[i][NPt/2+j] = ((Double_t)N[i][NPt/2+j] - (Double_t)N[i][NPt/2-j-1])/((Double_t)N[i][NPt/2+j] + (Double_t)N[i][NPt/2-j-1]);
       A[i][NPt/2-j-1] = -A[i][NPt/2+j];}
     GDA[i] = new TGraph (NPt, Del, A[i]);}

   for(Int_t i = 0; i < 4; i++){
     CDN->cd(i+1); GDN[i] -> Draw();
     CDA->cd(i+1); GDA[i] -> Draw();
     
     GDN[i]->SetTitle(Form("WirePlane%d Ee-Ep distribution",i));
     GDN[i]->GetXaxis()->SetTitle("Delta: (Ee-Ep)/2");
     GDN[i]->GetYaxis()->SetTitle("Counts");
     GDN[i]->GetXaxis()->CenterTitle();
     GDN[i]->GetYaxis()->CenterTitle();
     
     GDA[i]->SetTitle(Form("WirePlane%d total asymmetry",i));
     GDA[i]->GetXaxis()->SetTitle("Delta: (Ee-Ep)/2");
     GDA[i]->GetYaxis()->SetTitle("Asymmetry");
     GDA[i]->GetXaxis()->CenterTitle();
     GDA[i]->GetYaxis()->CenterTitle();

     C0WX->cd(i+1); H0WX[i]->Draw();
     C1WX->cd(i+1); H1WX[i]->Draw();
     C0EX->cd(i+1); H0EX[i]->Draw();
     C1EX->cd(i+1); H1EX[i]->Draw();
   }
}
