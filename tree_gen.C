#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <TFile.h>
#include <TTree.h>

using namespace std;

  struct Monitor{
    Int_t charge;
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
  //Int_t detector_num;
  vector<VDC> vdc;
  vector<Paddle> paddles;
};

struct BH_Event{
public:
  UInt_t event_num;
  vector<Monitor> monitors;
  Input input;
  Detector detector0;
  Detector detector1;
};

//ClassImp(BH_Event)
void InInput(ifstream& datafile, Input& input){ 
  datafile >> input.energy >> input.delta >> input.x >>  input.y >> input.theta >> input.phi;
  datafile.ignore(200, '\n');
  //cout << "An input has been read" << endl;
}

void InMonitor(ifstream& datafile, vector<Monitor>& monitors){
  //Monitor tempmonitors;
  Monitor temp;
  string temps;
  datafile >> temps >> temp.energy_m >> temp.x_m >> temp.y_m 
           >> temp.theta_m >> temp.phi_m;
  if(temps=="e+") temp.charge=1;
  else if(temps=="e-") temp.charge=-1;
  monitors.push_back(temp);

  datafile.ignore(200,'\n');
  //cout << "A monitor has been read" << endl;
}

void InVDC(ifstream& datafile, vector<VDC>& vdc){
  VDC temp;
  datafile >> temp.x_f >> temp.y_f >> temp.theta_f >> temp.phi_f;
  vdc.push_back(temp);
  datafile.ignore(200,'\n');
  //cout << "A VDC has been read" << endl;
}

void InPaddle(ifstream& datafile, vector<Paddle>& paddles){
  Paddle temp;
  datafile >> temp.paddle_num >> temp.Edep >> temp.Light_out >> temp.Hit_time;
  paddles.push_back(temp);
  datafile.ignore(200,'\n');
  cout << "A Paddle has been read"<< endl;
}
/*
void InDetectors(ifstream& datafile, Detector& detector, Int_t& detector_flag){
  Detector temp;
  datafile >> temp.detector_num;
  detector_flag = temp.detector_num;
  // only add detector_num here, info of VDC and Paddles are in InVDC and InPaddle function.
  detector.push_back(temp);
  datafile.ignore(200,'\n');
  cout << "A dectector array number "<< temp.detector_num <<" has been read" << endl;
  }*/

void SetDetector(ifstream& datafile, Int_t& detector_flag){
  datafile >> detector_flag;
  datafile.ignore(200,'\n');
}

void ProcessEvent(ifstream& datafile, BH_Event* event, bool& newevent, string& linetitle, Int_t& detector_flag)
{
  while(datafile && newevent == false){
  if(linetitle=="Event:") {
    datafile >> event->event_num;
    datafile.ignore(200,'\n');
    cout << "A new event " << event->event_num <<" has been read" << endl;
  }
  
  datafile >> linetitle;  
  //cout << "Linetitle is " << linetitle << endl; 
  if(linetitle=="Input:"){        
    InInput(datafile, event->input);
  }
  else if(linetitle=="Monitor:"){
    InMonitor(datafile, event->monitors);
  }
  else if(linetitle=="Detector:"){
    SetDetector(datafile, detector_flag);
    //InDetectors(datafile, event->detectors, detector_flag);
  }
  else if(linetitle=="VDC:"){
    
    //InVDC(datafile, event->detectors.at(detector_index).vdc);
    if(detector_flag == 0){
      InVDC(datafile, event->detector0.vdc);
    }
    else if(detector_flag == 1){
      InVDC(datafile,event->detector1.vdc);
    }
  }
  else if(linetitle=="Paddle:"){
    cout << "Adding paddles for detector number: " << detector_flag << endl;
    if(detector_flag == 0){
      InPaddle(datafile, event->detector0.paddles);
    }
    else if(detector_flag == 1){
      InPaddle(datafile,event->detector1.paddles);
    }
    //datafile.ignore(200,'\n');
  }
  else if(linetitle=="Event:"){
    newevent=true;
    //if(newevent==true)
      //cout << "newevent has been changed to true." << endl;
  }
  
  }
  //else datafile.ignore(200,'\n');
}

void MakeTree(void){
  //Declare data to be transformed to ttree
  //For each BH_Event, there are 1 event_num,1 or 2 Monitors, 1 Input, and 1 or 2 Detectors
  //For each Detector, there are 1 VDC, and possible multiple Paddles
  //open the .dat file needed
  ifstream datafile;
  string filename;

  cout<< "Input datafile name: (Format: BH.####.dat)" << endl;
  getline(cin,filename);
  cout << endl;
  datafile.open(filename.c_str());

  while(datafile.fail()){
    cout << "File not found: " << filename << endl
         <<"Please enter new name: ";
    getline(cin, filename);
    //cin.ignore(1, '\n');
    datafile.open(filename.c_str());
  }

  cout<<"BH data file " << filename << " is successfully opened." << endl;

  // open ROOT file to save TTree
  TFile* bh_ftree= new TFile("BHdata.root", "recreate");

  // create TTree
  TTree* tree= new TTree("tree", "BH data tree");
  BH_Event* event= new BH_Event;
  
  //create branch in TTree
  tree->Branch("Event_Branch",&event);

  //read and parse the commented lines infront
  char ch;
  ch=datafile.peek();

  while (datafile && ch == '#'){
    datafile.ignore(INT_MAX, '\n');
    ch=datafile.peek();
  }
  //cout << "Comment lines are ignored." <<endl;
  //string next; 
  //getline(datafile,next);
  //datafile >> next;
  //cout << "Next line is: "<< next << endl; 

  //read data and fill ttree
  bool newevent=false;
  string linetitle;
  datafile >> linetitle;
  
  while(datafile){
    //data fill for a single evnet:
    event->monitors.clear();
    //event->detectors.clear();
    event->detector1.paddles.clear();
    event->detector0.paddles.clear();
    event->detector1.vdc.clear();
    event->detector0.vdc.clear();
    
    Int_t detector_flag=-1;

    ProcessEvent(datafile, event, newevent,linetitle, detector_flag);

    tree->Fill();
    newevent=false;
    if(!datafile) break;
    //if (datafile.eof()) newevent=true;
    //cout << newevent << endl;  
  }
  

  //close datafile
  datafile.close();
  //write and close the ROOT file
  bh_ftree->Write();
  delete tree; 
}

