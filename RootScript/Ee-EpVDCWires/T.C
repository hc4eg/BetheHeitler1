#define T_cxx
#include "T.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

void T::Loop()
{
//   In a ROOT session, you can do:
//      Root > .L T.C
//      Root > T t
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

   TCanvas* CENum = new TCanvas("CENum", "ENum", 800, 400);
   TCanvas* C1Pad = new TCanvas("C1Pad", "PNum and ToF", 1600, 1200);
   C1Pad->Divide(2,2);
   TCanvas* C2Pad = new TCanvas("C2Pad", "Edep and Light", 1600, 1200);
   C2Pad->Divide(2,2);
   TCanvas* C3Pad = new TCanvas("C3Pad", "Hit Counts and Kinetic Energy", 1600, 1200);
   C3Pad->Divide(2,2);

   TH1F* H_ENum = new TH1F("H_ENum", "Event Number", 100, 0, 100000);
   H_ENum->SetBit(TH1::kCanRebin);
   TH1F* H_PNum[2];
   TH1F* H_Time[2];
   TH1F* H_Edep[2];
   TH1F* H_Light[2];
   TH1F* H_Count[2];
   TH1F* H_KE[2];

   for(int i = 0; i<2; i++){
     H_PNum[i] = new TH1F(Form("H_D%dPNum",i), Form("D%d Paddle Number distribution",i), 29, 0.0, 29.0);
     H_Time[i] = new TH1F(Form("H_D%dTime",i), Form("D%d Time of Flight",i), 120, 0., 12.);
     H_Edep[i] = new TH1F(Form("H_D%dEdep",i), Form("D%d Energy deposite", i), 100, 0., 10.);
     H_Light[i] = new TH1F(Form("H_D%dLight",i), Form("D%d Paddle Light", i), 100, 0., 10.);
     H_Count[i] = new TH1F(Form("H_D%dCount",i), Form("D%d Hit Paddle Count",i), 10, 0, 10);
     H_Count[i]->SetBit(TH1::kCanRebin);
     H_KE[i] = new TH1F(Form("H_D%dKE",i), Form("D%d Hodoscope avg. kinetic energy",i), 120, 0., 60.);
   }

   Long64_t nentries = fChain->GetEntriesFast();

   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;
      if( D0_V_>0 && D1_V_>0 && D0_P_>0 && D1_P_>0){
	H_ENum->Fill(ENum);
	//Detector0 Hodoscope
	for(int n = 0; n< D0_P_; n++){
	  H_PNum[0]->Fill(D0_P_PNum[n]);
	  H_Time[0]->Fill(D0_P_Time[n]);
	  H_Edep[0]->Fill(D0_P_Edep[n]);
	  H_Light[0]->Fill(D0_P_Light[n]);
	}
	//Detector1 Hodoscope
	for(int n = 0; n< D1_P_; n++){
	  H_PNum[1]->Fill(D1_P_PNum[n]);
	  H_Time[1]->Fill(D1_P_Time[n]);
	  H_Edep[1]->Fill(D1_P_Edep[n]);
	  H_Light[1]->Fill(D1_P_Light[n]);
	}
	H_Count[0]->Fill(D0_P_);
	H_Count[1]->Fill(D1_P_);
	if(HKE0.size()>0 && HKE1.size()>0){
	  H_KE[0]->Fill(HKE0.at(0));
	  H_KE[1]->Fill(HKE1.at(0));}
      }
   }
   //Draw the histograms
   CENum->cd(); H_ENum->Draw();
   C1Pad->cd(1); H_PNum[0]->Draw();
   C1Pad->cd(2); H_PNum[1]->Draw();
   C1Pad->cd(3); H_Time[0]->Draw();
   C1Pad->cd(4); H_Time[1]->Draw();

   C2Pad->cd(1); H_Edep[0]->Draw();
   C2Pad->cd(2); H_Edep[1]->Draw();
   C2Pad->cd(3); H_Light[0]->Draw();
   C2Pad->cd(4); H_Light[1]->Draw();

   C3Pad->cd(1); H_Count[0]->Draw();
   C3Pad->cd(2); H_Count[1]->Draw();
   C3Pad->cd(3); H_KE[0]->Draw();
   C3Pad->cd(4); H_KE[1]->Draw();
}
