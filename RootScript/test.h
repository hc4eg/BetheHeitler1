//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Thu Nov 19 16:04:49 2015 by ROOT version 5.34/14
// from TTree tree/root_tree
// found on file: Gamma10M.root
//////////////////////////////////////////////////////////

#ifndef test_h
#define test_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.
#include <TObject.h>

// Fixed size dimensions of array or collections stored in the TTree if any.
const Int_t kMaxdetector0_vdc = 1;
const Int_t kMaxdetector0_paddles = 1;
const Int_t kMaxdetector1_vdc = 1;
const Int_t kMaxdetector1_paddles = 1;

class test {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   // Declaration of leaf types
 //BH_Event        *Event_Branch;
   UInt_t          fUniqueID;
   UInt_t          fBits;
   Int_t           event_num;
   Int_t           monitor1_enumber;
   Float_t         monitor1_charge;
   Float_t         monitor1_energy_m;
   Float_t         monitor1_x_m;
   Float_t         monitor1_y_m;
   Float_t         monitor1_theta_m;
   Float_t         monitor1_phi_m;
   Int_t           monitor0_enumber;
   Float_t         monitor0_charge;
   Float_t         monitor0_energy_m;
   Float_t         monitor0_x_m;
   Float_t         monitor0_y_m;
   Float_t         monitor0_theta_m;
   Float_t         monitor0_phi_m;
   Float_t         input_energy;
   Float_t         input_delta;
   Float_t         input_x;
   Float_t         input_y;
   Float_t         input_theta;
   Float_t         input_phi;
   Int_t           detector0_vdc_;
   Float_t         detector0_vdc_x_f[kMaxdetector0_vdc];   //[detector0.vdc_]
   Float_t         detector0_vdc_y_f[kMaxdetector0_vdc];   //[detector0.vdc_]
   Float_t         detector0_vdc_theta_f[kMaxdetector0_vdc];   //[detector0.vdc_]
   Float_t         detector0_vdc_phi_f[kMaxdetector0_vdc];   //[detector0.vdc_]
   Int_t           detector0_paddles_;
   Int_t           detector0_paddles_paddle_num[kMaxdetector0_paddles];   //[detector0.paddles_]
   Float_t         detector0_paddles_Edep[kMaxdetector0_paddles];   //[detector0.paddles_]
   Float_t         detector0_paddles_Light_out[kMaxdetector0_paddles];   //[detector0.paddles_]
   Float_t         detector0_paddles_Hit_time[kMaxdetector0_paddles];   //[detector0.paddles_]
   Int_t           detector1_vdc_;
   Float_t         detector1_vdc_x_f[kMaxdetector1_vdc];   //[detector1.vdc_]
   Float_t         detector1_vdc_y_f[kMaxdetector1_vdc];   //[detector1.vdc_]
   Float_t         detector1_vdc_theta_f[kMaxdetector1_vdc];   //[detector1.vdc_]
   Float_t         detector1_vdc_phi_f[kMaxdetector1_vdc];   //[detector1.vdc_]
   Int_t           detector1_paddles_;
   Int_t           detector1_paddles_paddle_num[kMaxdetector1_paddles];   //[detector1.paddles_]
   Float_t         detector1_paddles_Edep[kMaxdetector1_paddles];   //[detector1.paddles_]
   Float_t         detector1_paddles_Light_out[kMaxdetector1_paddles];   //[detector1.paddles_]
   Float_t         detector1_paddles_Hit_time[kMaxdetector1_paddles];   //[detector1.paddles_]

   // List of branches
   TBranch        *b_Event_Branch_fUniqueID;   //!
   TBranch        *b_Event_Branch_fBits;   //!
   TBranch        *b_Event_Branch_event_num;   //!
   TBranch        *b_Event_Branch_monitor1_enumber;   //!
   TBranch        *b_Event_Branch_monitor1_charge;   //!
   TBranch        *b_Event_Branch_monitor1_energy_m;   //!
   TBranch        *b_Event_Branch_monitor1_x_m;   //!
   TBranch        *b_Event_Branch_monitor1_y_m;   //!
   TBranch        *b_Event_Branch_monitor1_theta_m;   //!
   TBranch        *b_Event_Branch_monitor1_phi_m;   //!
   TBranch        *b_Event_Branch_monitor0_enumber;   //!
   TBranch        *b_Event_Branch_monitor0_charge;   //!
   TBranch        *b_Event_Branch_monitor0_energy_m;   //!
   TBranch        *b_Event_Branch_monitor0_x_m;   //!
   TBranch        *b_Event_Branch_monitor0_y_m;   //!
   TBranch        *b_Event_Branch_monitor0_theta_m;   //!
   TBranch        *b_Event_Branch_monitor0_phi_m;   //!
   TBranch        *b_Event_Branch_input_energy;   //!
   TBranch        *b_Event_Branch_input_delta;   //!
   TBranch        *b_Event_Branch_input_x;   //!
   TBranch        *b_Event_Branch_input_y;   //!
   TBranch        *b_Event_Branch_input_theta;   //!
   TBranch        *b_Event_Branch_input_phi;   //!
   TBranch        *b_Event_Branch_detector0_vdc_;   //!
   TBranch        *b_detector0_vdc_x_f;   //!
   TBranch        *b_detector0_vdc_y_f;   //!
   TBranch        *b_detector0_vdc_theta_f;   //!
   TBranch        *b_detector0_vdc_phi_f;   //!
   TBranch        *b_Event_Branch_detector0_paddles_;   //!
   TBranch        *b_detector0_paddles_paddle_num;   //!
   TBranch        *b_detector0_paddles_Edep;   //!
   TBranch        *b_detector0_paddles_Light_out;   //!
   TBranch        *b_detector0_paddles_Hit_time;   //!
   TBranch        *b_Event_Branch_detector1_vdc_;   //!
   TBranch        *b_detector1_vdc_x_f;   //!
   TBranch        *b_detector1_vdc_y_f;   //!
   TBranch        *b_detector1_vdc_theta_f;   //!
   TBranch        *b_detector1_vdc_phi_f;   //!
   TBranch        *b_Event_Branch_detector1_paddles_;   //!
   TBranch        *b_detector1_paddles_paddle_num;   //!
   TBranch        *b_detector1_paddles_Edep;   //!
   TBranch        *b_detector1_paddles_Light_out;   //!
   TBranch        *b_detector1_paddles_Hit_time;   //!

   test(TTree *tree=0);
   virtual ~test();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef test_cxx
test::test(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("Gamma10M.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("Gamma10M.root");
      }
      f->GetObject("tree",tree);

   }
   Init(tree);
}

test::~test()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t test::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t test::LoadTree(Long64_t entry)
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

void test::Init(TTree *tree)
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

   fChain->SetBranchAddress("fUniqueID", &fUniqueID, &b_Event_Branch_fUniqueID);
   fChain->SetBranchAddress("fBits", &fBits, &b_Event_Branch_fBits);
   fChain->SetBranchAddress("event_num", &event_num, &b_Event_Branch_event_num);
   fChain->SetBranchAddress("monitor1.enumber", &monitor1_enumber, &b_Event_Branch_monitor1_enumber);
   fChain->SetBranchAddress("monitor1.charge", &monitor1_charge, &b_Event_Branch_monitor1_charge);
   fChain->SetBranchAddress("monitor1.energy_m", &monitor1_energy_m, &b_Event_Branch_monitor1_energy_m);
   fChain->SetBranchAddress("monitor1.x_m", &monitor1_x_m, &b_Event_Branch_monitor1_x_m);
   fChain->SetBranchAddress("monitor1.y_m", &monitor1_y_m, &b_Event_Branch_monitor1_y_m);
   fChain->SetBranchAddress("monitor1.theta_m", &monitor1_theta_m, &b_Event_Branch_monitor1_theta_m);
   fChain->SetBranchAddress("monitor1.phi_m", &monitor1_phi_m, &b_Event_Branch_monitor1_phi_m);
   fChain->SetBranchAddress("monitor0.enumber", &monitor0_enumber, &b_Event_Branch_monitor0_enumber);
   fChain->SetBranchAddress("monitor0.charge", &monitor0_charge, &b_Event_Branch_monitor0_charge);
   fChain->SetBranchAddress("monitor0.energy_m", &monitor0_energy_m, &b_Event_Branch_monitor0_energy_m);
   fChain->SetBranchAddress("monitor0.x_m", &monitor0_x_m, &b_Event_Branch_monitor0_x_m);
   fChain->SetBranchAddress("monitor0.y_m", &monitor0_y_m, &b_Event_Branch_monitor0_y_m);
   fChain->SetBranchAddress("monitor0.theta_m", &monitor0_theta_m, &b_Event_Branch_monitor0_theta_m);
   fChain->SetBranchAddress("monitor0.phi_m", &monitor0_phi_m, &b_Event_Branch_monitor0_phi_m);
   fChain->SetBranchAddress("input.energy", &input_energy, &b_Event_Branch_input_energy);
   fChain->SetBranchAddress("input.delta", &input_delta, &b_Event_Branch_input_delta);
   fChain->SetBranchAddress("input.x", &input_x, &b_Event_Branch_input_x);
   fChain->SetBranchAddress("input.y", &input_y, &b_Event_Branch_input_y);
   fChain->SetBranchAddress("input.theta", &input_theta, &b_Event_Branch_input_theta);
   fChain->SetBranchAddress("input.phi", &input_phi, &b_Event_Branch_input_phi);
   fChain->SetBranchAddress("detector0.vdc", &detector0_vdc_, &b_Event_Branch_detector0_vdc_);
   fChain->SetBranchAddress("detector0.vdc.x_f", &detector0_vdc_x_f, &b_detector0_vdc_x_f);
   fChain->SetBranchAddress("detector0.vdc.y_f", &detector0_vdc_y_f, &b_detector0_vdc_y_f);
   fChain->SetBranchAddress("detector0.vdc.theta_f", &detector0_vdc_theta_f, &b_detector0_vdc_theta_f);
   fChain->SetBranchAddress("detector0.vdc.phi_f", &detector0_vdc_phi_f, &b_detector0_vdc_phi_f);
   fChain->SetBranchAddress("detector0.paddles", &detector0_paddles_, &b_Event_Branch_detector0_paddles_);
   fChain->SetBranchAddress("detector0.paddles.paddle_num", &detector0_paddles_paddle_num, &b_detector0_paddles_paddle_num);
   fChain->SetBranchAddress("detector0.paddles.Edep", &detector0_paddles_Edep, &b_detector0_paddles_Edep);
   fChain->SetBranchAddress("detector0.paddles.Light_out", &detector0_paddles_Light_out, &b_detector0_paddles_Light_out);
   fChain->SetBranchAddress("detector0.paddles.Hit_time", &detector0_paddles_Hit_time, &b_detector0_paddles_Hit_time);
   fChain->SetBranchAddress("detector1.vdc", &detector1_vdc_, &b_Event_Branch_detector1_vdc_);
   fChain->SetBranchAddress("detector1.vdc.x_f", &detector1_vdc_x_f, &b_detector1_vdc_x_f);
   fChain->SetBranchAddress("detector1.vdc.y_f", &detector1_vdc_y_f, &b_detector1_vdc_y_f);
   fChain->SetBranchAddress("detector1.vdc.theta_f", &detector1_vdc_theta_f, &b_detector1_vdc_theta_f);
   fChain->SetBranchAddress("detector1.vdc.phi_f", &detector1_vdc_phi_f, &b_detector1_vdc_phi_f);
   fChain->SetBranchAddress("detector1.paddles", &detector1_paddles_, &b_Event_Branch_detector1_paddles_);
   fChain->SetBranchAddress("detector1.paddles.paddle_num", &detector1_paddles_paddle_num, &b_detector1_paddles_paddle_num);
   fChain->SetBranchAddress("detector1.paddles.Edep", &detector1_paddles_Edep, &b_detector1_paddles_Edep);
   fChain->SetBranchAddress("detector1.paddles.Light_out", &detector1_paddles_Light_out, &b_detector1_paddles_Light_out);
   fChain->SetBranchAddress("detector1.paddles.Hit_time", &detector1_paddles_Hit_time, &b_detector1_paddles_Hit_time);
   Notify();
}

Bool_t test::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void test::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t test::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef test_cxx
