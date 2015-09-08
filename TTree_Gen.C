#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include "TFile.h"
#include "TTree.h"

using namespace std;

  enum Charge{e_p,e_n};  //positron or electron
  
  struct Monitor{
    Charge charge;
    Float_t energy_m;
    Float_t x_m;
    Float_t y_m;
    Float_t theta_m;
    Float_t phi_m;    
  };

  struct Input{
    Float_t energy;
    Float_t delta;
    Float_t x;
    Float_t y;
    Float_t theta;
    Float_t phi;
  };
  
  struct VDC{
    Float_t x_f;
    Float_t y_f;
    Float_t theta_f;
    Float_t phi_f;
  };
 
  struct Paddle{
    Int_t paddle_num;
    Float_t Edep;
    Float_t Light_out;
    Float_t Hit_time;
  };

  struct Detector{
    VDC vdc;
    vector<Paddle> paddles;
  };

  struct BH_Event{
    Int_t event_num;
    vector<Monitor> monitors;
    Input input;
    vector<Detector> detectors;    
  };

void MakeTree(void){
  //declare data to be transformed to ttree
  //For each BH_Event, there are 1 event_num, 1 or 2 Monitors, 1 Input, and 1 or 2 Detectors
  //For each Detector, there are 1 VDC, and multiple Paddles

  BH_Event event;

  //open the .dat file needed
  ifstream datafile;
  string filename;

  cout << "Input data file name: (Format: BH.####.dat)";
  getline(cin,filename);
  cout << endl;
  datafile.open(filename.c_str());

  while(datafile.fail()){
    datafile.clear();
    cout << "File not found: " << filename << endl
         << "Please enter new name: ";
    getline(cin,filename);
    cin.ignore(1,'\n');
    datafile.open(filename.c_str());
  }

  cout << "BH data file" << filename <<" is successfully opened."<< endl;
  
  //open ROOT file to save TTree in
  TFile *bh_ftree= new TFile("BHdata.root","recreate");

  //create TTree
  TTree *tree= new TTree("tree","BH data TTree");

  //create branch in TTree
  tree->Branch("Event_Branch",&event, "");

  //read the whole data file, and fill TTree

  //close datafile
    datafile.close();
  //write and close the ROOT file

}
