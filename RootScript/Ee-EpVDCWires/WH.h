//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Sat Mar 19 01:07:42 2016 by ROOT version 5.34/14
// from TTree T/root_tree
// found on file: Pair.0255.root
//////////////////////////////////////////////////////////

#ifndef _WH_H_
#define _WH_H_

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <vector>

// Header file for the classes stored in the TTree if any.
#include <TObject.h>

// Fixed size dimensions of array or collections stored in the TTree if any.
// Changed to larger value to suit multi root files.
const Int_t kMaxD0_W = 200;
const Int_t kMaxD0_P = 20;
const Int_t kMaxD1_W = 200;
const Int_t kMaxD1_P = 20;

class WH {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   // Declaration of leaf types
 //BH_Event        *B;
   UInt_t          fUniqueID;
   UInt_t          fBits;
   Int_t           ENum;
   Float_t         M0_Charge;
   Float_t         M0_Energy;
   Float_t         M0_X;
   Float_t         M0_Y;
   Float_t         M0_Theta;
   Float_t         M0_Phi;
   Float_t         M1_Charge;
   Float_t         M1_Energy;
   Float_t         M1_X;
   Float_t         M1_Y;
   Float_t         M1_Theta;
   Float_t         M1_Phi;
   Float_t         I0_Energy;
   Float_t         I0_Delta;
   Float_t         I0_X;
   Float_t         I0_Y;
   Float_t         I0_Theta;
   Float_t         I0_Phi;
   Float_t         I1_Energy;
   Float_t         I1_Delta;
   Float_t         I1_X;
   Float_t         I1_Y;
   Float_t         I1_Theta;
   Float_t         I1_Phi;
   Int_t           D0_W_;
   Int_t           D0_W_WirePlane[kMaxD0_W];   //[D0.W_]
   Int_t           D0_W_WireNum[kMaxD0_W];   //[D0.W_]
   Float_t         D0_W_X[kMaxD0_W];   //[D0.W_]
   Float_t         D0_W_Y[kMaxD0_W];   //[D0.W_]
   Float_t         D0_W_Edep[kMaxD0_W];   //[D0.W_]
   Float_t         D0_W_KE[kMaxD0_W];   //[D0.W_]
   Float_t         D0_W_ToF[kMaxD0_W];   //[D0.W_]
   Float_t         D0_W_Charge[kMaxD0_W];   //[D0.W_]
   Int_t           D0_W_Particle[kMaxD0_W];   //[D0.W_]
   Int_t           D0_P_;
   Int_t           D0_P_PNum[kMaxD0_P];   //[D0.P_]
   Float_t         D0_P_Edep[kMaxD0_P];   //[D0.P_]
   Float_t         D0_P_Light[kMaxD0_P];   //[D0.P_]
   Float_t         D0_P_Time[kMaxD0_P];   //[D0.P_]
   Int_t           D1_W_;
   Int_t           D1_W_WirePlane[kMaxD1_W];   //[D1.W_]
   Int_t           D1_W_WireNum[kMaxD1_W];   //[D1.W_]
   Float_t         D1_W_X[kMaxD1_W];   //[D1.W_]
   Float_t         D1_W_Y[kMaxD1_W];   //[D1.W_]
   Float_t         D1_W_Edep[kMaxD1_W];   //[D1.W_]
   Float_t         D1_W_KE[kMaxD1_W];   //[D1.W_]
   Float_t         D1_W_ToF[kMaxD1_W];   //[D1.W_]
   Float_t         D1_W_Charge[kMaxD1_W];   //[D1.W_]
   Int_t           D1_W_Particle[kMaxD1_W];   //[D1.W_]
   Int_t           D1_P_;
   Int_t           D1_P_PNum[kMaxD1_P];   //[D1.P_]
   Float_t         D1_P_Edep[kMaxD1_P];   //[D1.P_]
   Float_t         D1_P_Light[kMaxD1_P];   //[D1.P_]
   Float_t         D1_P_Time[kMaxD1_P];   //[D1.P_]
   std::vector<Double_t> HKE0;
   std::vector<Double_t> HKE1;

   // List of branches
   TBranch        *b_B_fUniqueID;   //!
   TBranch        *b_B_fBits;   //!
   TBranch        *b_B_ENum;   //!
   TBranch        *b_B_M0_Charge;   //!
   TBranch        *b_B_M0_Energy;   //!
   TBranch        *b_B_M0_X;   //!
   TBranch        *b_B_M0_Y;   //!
   TBranch        *b_B_M0_Theta;   //!
   TBranch        *b_B_M0_Phi;   //!
   TBranch        *b_B_M1_Charge;   //!
   TBranch        *b_B_M1_Energy;   //!
   TBranch        *b_B_M1_X;   //!
   TBranch        *b_B_M1_Y;   //!
   TBranch        *b_B_M1_Theta;   //!
   TBranch        *b_B_M1_Phi;   //!
   TBranch        *b_B_I0_Energy;   //!
   TBranch        *b_B_I0_Delta;   //!
   TBranch        *b_B_I0_X;   //!
   TBranch        *b_B_I0_Y;   //!
   TBranch        *b_B_I0_Theta;   //!
   TBranch        *b_B_I0_Phi;   //!
   TBranch        *b_B_I1_Energy;   //!
   TBranch        *b_B_I1_Delta;   //!
   TBranch        *b_B_I1_X;   //!
   TBranch        *b_B_I1_Y;   //!
   TBranch        *b_B_I1_Theta;   //!
   TBranch        *b_B_I1_Phi;   //!
   TBranch        *b_B_D0_W_;   //!
   TBranch        *b_D0_W_WirePlane;   //!
   TBranch        *b_D0_W_WireNum;   //!
   TBranch        *b_D0_W_X;   //!
   TBranch        *b_D0_W_Y;   //!
   TBranch        *b_D0_W_Edep;   //!
   TBranch        *b_D0_W_KE;   //!
   TBranch        *b_D0_W_ToF;   //!
   TBranch        *b_D0_W_Charge;   //!
   TBranch        *b_D0_W_Particle;   //!
   TBranch        *b_B_D0_P_;   //!
   TBranch        *b_D0_P_PNum;   //!
   TBranch        *b_D0_P_Edep;   //!
   TBranch        *b_D0_P_Light;   //!
   TBranch        *b_D0_P_Time;   //!
   TBranch        *b_B_D1_W_;   //!
   TBranch        *b_D1_W_WirePlane;   //!
   TBranch        *b_D1_W_WireNum;   //!
   TBranch        *b_D1_W_X;   //!
   TBranch        *b_D1_W_Y;   //!
   TBranch        *b_D1_W_Edep;   //!
   TBranch        *b_D1_W_KE;   //!
   TBranch        *b_D1_W_ToF;   //!
   TBranch        *b_D1_W_Charge;   //!
   TBranch        *b_D1_W_Particle;   //!
   TBranch        *b_B_D1_P_;   //!
   TBranch        *b_D1_P_PNum;   //!
   TBranch        *b_D1_P_Edep;   //!
   TBranch        *b_D1_P_Light;   //!
   TBranch        *b_D1_P_Time;   //!
   TBranch        *b_B_HKE0;   //!
   TBranch        *b_B_HKE1;   //!

   WH(TTree *tree=0);
   WH(Int_t runNumber, TTree* tree=0);
   virtual ~WH();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     PlotAsym();
   virtual void     PlotInput();
   virtual void     PlotMonitor();
   virtual void     PlotPaddle();
   virtual void     PlotVDC();
   virtual void     PlotVDCMAXKE();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

WH::WH(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("Pair.0255.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("Pair.0255.root");
      }
      f->GetObject("T",tree);

   }
   Init(tree);
}

WH::WH(Int_t runNumber, TTree* tree)
{
     if (tree == 0) {
       // original Pair.%4d.root files
       TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject(Form("Pair.%04d.root",runNumber));
       // Theoretic cross section random generated pair as primary vertex
       //TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject(Form("TheoXsec.%03d.root",runNumber));
       //TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject(Form("TheoXsecLead.%03d.root",runNumber));
      if (!f || !f->IsOpen()) {
	f = new TFile(Form("Pair.%04d.root",runNumber));
	//f = new TFile(Form("TheoXsec.%03d.root",runNumber));
	//f = new TFile(Form("TheoXsecLead.%03d.root",runNumber));
      }
      //f->GetObject("T",tree);
      tree = (TTree*)gDirectory->Get("T");
   }
   Init(tree);
}

WH::~WH()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t WH::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t WH::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void WH::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("fUniqueID", &fUniqueID, &b_B_fUniqueID);
   fChain->SetBranchAddress("fBits", &fBits, &b_B_fBits);
   fChain->SetBranchAddress("ENum", &ENum, &b_B_ENum);
   fChain->SetBranchAddress("M0.Charge", &M0_Charge, &b_B_M0_Charge);
   fChain->SetBranchAddress("M0.Energy", &M0_Energy, &b_B_M0_Energy);
   fChain->SetBranchAddress("M0.X", &M0_X, &b_B_M0_X);
   fChain->SetBranchAddress("M0.Y", &M0_Y, &b_B_M0_Y);
   fChain->SetBranchAddress("M0.Theta", &M0_Theta, &b_B_M0_Theta);
   fChain->SetBranchAddress("M0.Phi", &M0_Phi, &b_B_M0_Phi);
   fChain->SetBranchAddress("M1.Charge", &M1_Charge, &b_B_M1_Charge);
   fChain->SetBranchAddress("M1.Energy", &M1_Energy, &b_B_M1_Energy);
   fChain->SetBranchAddress("M1.X", &M1_X, &b_B_M1_X);
   fChain->SetBranchAddress("M1.Y", &M1_Y, &b_B_M1_Y);
   fChain->SetBranchAddress("M1.Theta", &M1_Theta, &b_B_M1_Theta);
   fChain->SetBranchAddress("M1.Phi", &M1_Phi, &b_B_M1_Phi);
   fChain->SetBranchAddress("I0.Energy", &I0_Energy, &b_B_I0_Energy);
   fChain->SetBranchAddress("I0.Delta", &I0_Delta, &b_B_I0_Delta);
   fChain->SetBranchAddress("I0.X", &I0_X, &b_B_I0_X);
   fChain->SetBranchAddress("I0.Y", &I0_Y, &b_B_I0_Y);
   fChain->SetBranchAddress("I0.Theta", &I0_Theta, &b_B_I0_Theta);
   fChain->SetBranchAddress("I0.Phi", &I0_Phi, &b_B_I0_Phi);
   fChain->SetBranchAddress("I1.Energy", &I1_Energy, &b_B_I1_Energy);
   fChain->SetBranchAddress("I1.Delta", &I1_Delta, &b_B_I1_Delta);
   fChain->SetBranchAddress("I1.X", &I1_X, &b_B_I1_X);
   fChain->SetBranchAddress("I1.Y", &I1_Y, &b_B_I1_Y);
   fChain->SetBranchAddress("I1.Theta", &I1_Theta, &b_B_I1_Theta);
   fChain->SetBranchAddress("I1.Phi", &I1_Phi, &b_B_I1_Phi);
   fChain->SetBranchAddress("D0.W", &D0_W_, &b_B_D0_W_);
   fChain->SetBranchAddress("D0.W.WirePlane", D0_W_WirePlane, &b_D0_W_WirePlane);
   fChain->SetBranchAddress("D0.W.WireNum", D0_W_WireNum, &b_D0_W_WireNum);
   fChain->SetBranchAddress("D0.W.X", D0_W_X, &b_D0_W_X);
   fChain->SetBranchAddress("D0.W.Y", D0_W_Y, &b_D0_W_Y);
   fChain->SetBranchAddress("D0.W.Edep", D0_W_Edep, &b_D0_W_Edep);
   fChain->SetBranchAddress("D0.W.KE", D0_W_KE, &b_D0_W_KE);
   fChain->SetBranchAddress("D0.W.ToF", D0_W_ToF, &b_D0_W_ToF);
   fChain->SetBranchAddress("D0.W.Charge", D0_W_Charge, &b_D0_W_Charge);
   fChain->SetBranchAddress("D0.W.Particle", D0_W_Particle, &b_D0_W_Particle);
   fChain->SetBranchAddress("D0.P", &D0_P_, &b_B_D0_P_);
   fChain->SetBranchAddress("D0.P.PNum", D0_P_PNum, &b_D0_P_PNum);
   fChain->SetBranchAddress("D0.P.Edep", D0_P_Edep, &b_D0_P_Edep);
   fChain->SetBranchAddress("D0.P.Light", D0_P_Light, &b_D0_P_Light);
   fChain->SetBranchAddress("D0.P.Time", D0_P_Time, &b_D0_P_Time);
   fChain->SetBranchAddress("D1.W", &D1_W_, &b_B_D1_W_);
   fChain->SetBranchAddress("D1.W.WirePlane", D1_W_WirePlane, &b_D1_W_WirePlane);
   fChain->SetBranchAddress("D1.W.WireNum", D1_W_WireNum, &b_D1_W_WireNum);
   fChain->SetBranchAddress("D1.W.X", D1_W_X, &b_D1_W_X);
   fChain->SetBranchAddress("D1.W.Y", D1_W_Y, &b_D1_W_Y);
   fChain->SetBranchAddress("D1.W.Edep", D1_W_Edep, &b_D1_W_Edep);
   fChain->SetBranchAddress("D1.W.KE", D1_W_KE, &b_D1_W_KE);
   fChain->SetBranchAddress("D1.W.ToF", D1_W_ToF, &b_D1_W_ToF);
   fChain->SetBranchAddress("D1.W.Charge", D1_W_Charge, &b_D1_W_Charge);
   fChain->SetBranchAddress("D1.W.Particle", D1_W_Particle, &b_D1_W_Particle);
   fChain->SetBranchAddress("D1.P", &D1_P_, &b_B_D1_P_);
   fChain->SetBranchAddress("D1.P.PNum", D1_P_PNum, &b_D1_P_PNum);
   fChain->SetBranchAddress("D1.P.Edep", D1_P_Edep, &b_D1_P_Edep);
   fChain->SetBranchAddress("D1.P.Light", D1_P_Light, &b_D1_P_Light);
   fChain->SetBranchAddress("D1.P.Time", D1_P_Time, &b_D1_P_Time);
   fChain->SetBranchAddress("HKE0", &HKE0, &b_B_HKE0);
   fChain->SetBranchAddress("HKE1", &HKE1, &b_B_HKE1);
   Notify();
}

Bool_t WH::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void WH::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t WH::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif
