//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Thu Mar  3 16:08:34 2016 by ROOT version 5.34/14
// from TTree T/root_tree
// found on file: Asym0PerAndPadNoOver.root
//////////////////////////////////////////////////////////

#ifndef T_h
#define T_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.
#include <TObject.h>

// Fixed size dimensions of array or collections stored in the TTree if any.
const Int_t kMaxD0_V = 1;
const Int_t kMaxD0_P = 5;
const Int_t kMaxD1_V = 1;
const Int_t kMaxD1_P = 4;

class T {
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
   Int_t           D0_V_;
   Float_t         D0_V_X[kMaxD0_V];   //[D0.V_]
   Float_t         D0_V_Y[kMaxD0_V];   //[D0.V_]
   Float_t         D0_V_Theta[kMaxD0_V];   //[D0.V_]
   Float_t         D0_V_Phi[kMaxD0_V];   //[D0.V_]
   Float_t         D0_V_E0u[kMaxD0_V];   //[D0.V_]
   Float_t         D0_V_E0v[kMaxD0_V];   //[D0.V_]
   Float_t         D0_V_E1u[kMaxD0_V];   //[D0.V_]
   Float_t         D0_V_E1v[kMaxD0_V];   //[D0.V_]
   Float_t         D0_V_KE[kMaxD0_V];   //[D0.V_]
   Float_t         D0_V_ToF[kMaxD0_V];   //[D0.V_]
   Float_t         D0_V_Charge[kMaxD0_V];   //[D0.V_]
   Int_t           D0_P_;
   Int_t           D0_P_PNum[kMaxD0_P];   //[D0.P_]
   Float_t         D0_P_Edep[kMaxD0_P];   //[D0.P_]
   Float_t         D0_P_Light[kMaxD0_P];   //[D0.P_]
   Float_t         D0_P_Time[kMaxD0_P];   //[D0.P_]
   Int_t           D1_V_;
   Float_t         D1_V_X[kMaxD1_V];   //[D1.V_]
   Float_t         D1_V_Y[kMaxD1_V];   //[D1.V_]
   Float_t         D1_V_Theta[kMaxD1_V];   //[D1.V_]
   Float_t         D1_V_Phi[kMaxD1_V];   //[D1.V_]
   Float_t         D1_V_E0u[kMaxD1_V];   //[D1.V_]
   Float_t         D1_V_E0v[kMaxD1_V];   //[D1.V_]
   Float_t         D1_V_E1u[kMaxD1_V];   //[D1.V_]
   Float_t         D1_V_E1v[kMaxD1_V];   //[D1.V_]
   Float_t         D1_V_KE[kMaxD1_V];   //[D1.V_]
   Float_t         D1_V_ToF[kMaxD1_V];   //[D1.V_]
   Float_t         D1_V_Charge[kMaxD1_V];   //[D1.V_]
   Int_t           D1_P_;
   Int_t           D1_P_PNum[kMaxD1_P];   //[D1.P_]
   Float_t         D1_P_Edep[kMaxD1_P];   //[D1.P_]
   Float_t         D1_P_Light[kMaxD1_P];   //[D1.P_]
   Float_t         D1_P_Time[kMaxD1_P];   //[D1.P_]
   vector<Double_t> HKE0;
   vector<Double_t> HKE1;

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
   TBranch        *b_B_D0_V_;   //!
   TBranch        *b_D0_V_X;   //!
   TBranch        *b_D0_V_Y;   //!
   TBranch        *b_D0_V_Theta;   //!
   TBranch        *b_D0_V_Phi;   //!
   TBranch        *b_D0_V_E0u;   //!
   TBranch        *b_D0_V_E0v;   //!
   TBranch        *b_D0_V_E1u;   //!
   TBranch        *b_D0_V_E1v;   //!
   TBranch        *b_D0_V_KE;   //!
   TBranch        *b_D0_V_ToF;   //!
   TBranch        *b_D0_V_Charge;   //!
   TBranch        *b_B_D0_P_;   //!
   TBranch        *b_D0_P_PNum;   //!
   TBranch        *b_D0_P_Edep;   //!
   TBranch        *b_D0_P_Light;   //!
   TBranch        *b_D0_P_Time;   //!
   TBranch        *b_B_D1_V_;   //!
   TBranch        *b_D1_V_X;   //!
   TBranch        *b_D1_V_Y;   //!
   TBranch        *b_D1_V_Theta;   //!
   TBranch        *b_D1_V_Phi;   //!
   TBranch        *b_D1_V_E0u;   //!
   TBranch        *b_D1_V_E0v;   //!
   TBranch        *b_D1_V_E1u;   //!
   TBranch        *b_D1_V_E1v;   //!
   TBranch        *b_D1_V_KE;   //!
   TBranch        *b_D1_V_ToF;   //!
   TBranch        *b_D1_V_Charge;   //!
   TBranch        *b_B_D1_P_;   //!
   TBranch        *b_D1_P_PNum;   //!
   TBranch        *b_D1_P_Edep;   //!
   TBranch        *b_D1_P_Light;   //!
   TBranch        *b_D1_P_Time;   //!
   TBranch        *b_B_HKE0;   //!
   TBranch        *b_B_HKE1;   //!

   T(TTree *tree=0);
   virtual ~T();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef T_cxx
T::T(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("Asym0PerAndPadNoOver.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("Asym0PerAndPadNoOver.root");
      }
      f->GetObject("T",tree);

   }
   Init(tree);
}

T::~T()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t T::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t T::LoadTree(Long64_t entry)
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

void T::Init(TTree *tree)
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
   fChain->SetBranchAddress("D0.V", &D0_V_, &b_B_D0_V_);
   fChain->SetBranchAddress("D0.V.X", D0_V_X, &b_D0_V_X);
   fChain->SetBranchAddress("D0.V.Y", D0_V_Y, &b_D0_V_Y);
   fChain->SetBranchAddress("D0.V.Theta", D0_V_Theta, &b_D0_V_Theta);
   fChain->SetBranchAddress("D0.V.Phi", D0_V_Phi, &b_D0_V_Phi);
   fChain->SetBranchAddress("D0.V.E0u", D0_V_E0u, &b_D0_V_E0u);
   fChain->SetBranchAddress("D0.V.E0v", D0_V_E0v, &b_D0_V_E0v);
   fChain->SetBranchAddress("D0.V.E1u", D0_V_E1u, &b_D0_V_E1u);
   fChain->SetBranchAddress("D0.V.E1v", D0_V_E1v, &b_D0_V_E1v);
   fChain->SetBranchAddress("D0.V.KE", D0_V_KE, &b_D0_V_KE);
   fChain->SetBranchAddress("D0.V.ToF", D0_V_ToF, &b_D0_V_ToF);
   fChain->SetBranchAddress("D0.V.Charge", D0_V_Charge, &b_D0_V_Charge);
   fChain->SetBranchAddress("D0.P", &D0_P_, &b_B_D0_P_);
   fChain->SetBranchAddress("D0.P.PNum", D0_P_PNum, &b_D0_P_PNum);
   fChain->SetBranchAddress("D0.P.Edep", D0_P_Edep, &b_D0_P_Edep);
   fChain->SetBranchAddress("D0.P.Light", D0_P_Light, &b_D0_P_Light);
   fChain->SetBranchAddress("D0.P.Time", D0_P_Time, &b_D0_P_Time);
   fChain->SetBranchAddress("D1.V", &D1_V_, &b_B_D1_V_);
   fChain->SetBranchAddress("D1.V.X", D1_V_X, &b_D1_V_X);
   fChain->SetBranchAddress("D1.V.Y", D1_V_Y, &b_D1_V_Y);
   fChain->SetBranchAddress("D1.V.Theta", D1_V_Theta, &b_D1_V_Theta);
   fChain->SetBranchAddress("D1.V.Phi", D1_V_Phi, &b_D1_V_Phi);
   fChain->SetBranchAddress("D1.V.E0u", D1_V_E0u, &b_D1_V_E0u);
   fChain->SetBranchAddress("D1.V.E0v", D1_V_E0v, &b_D1_V_E0v);
   fChain->SetBranchAddress("D1.V.E1u", D1_V_E1u, &b_D1_V_E1u);
   fChain->SetBranchAddress("D1.V.E1v", D1_V_E1v, &b_D1_V_E1v);
   fChain->SetBranchAddress("D1.V.KE", D1_V_KE, &b_D1_V_KE);
   fChain->SetBranchAddress("D1.V.ToF", D1_V_ToF, &b_D1_V_ToF);
   fChain->SetBranchAddress("D1.V.Charge", D1_V_Charge, &b_D1_V_Charge);
   fChain->SetBranchAddress("D1.P", &D1_P_, &b_B_D1_P_);
   fChain->SetBranchAddress("D1.P.PNum", D1_P_PNum, &b_D1_P_PNum);
   fChain->SetBranchAddress("D1.P.Edep", D1_P_Edep, &b_D1_P_Edep);
   fChain->SetBranchAddress("D1.P.Light", D1_P_Light, &b_D1_P_Light);
   fChain->SetBranchAddress("D1.P.Time", D1_P_Time, &b_D1_P_Time);
   fChain->SetBranchAddress("HKE0", &HKE0, &b_B_HKE0);
   fChain->SetBranchAddress("HKE1", &HKE1, &b_B_HKE1);
   Notify();
}

Bool_t T::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void T::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t T::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef T_cxx
