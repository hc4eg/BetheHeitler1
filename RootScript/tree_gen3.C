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
  Int_t detector_num;
  VDC vdc;
  vector<Paddle> paddles;
};

struct BH_Event{
public:
  UInt_t event_num;
  vector<Monitor> monitors;
  Input input;
  vector<Detector> detectors;
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

void InVDC(ifstream& datafile, VDC& vdc){
  datafile >> vdc.x_f >> vdc.y_f >> vdc.theta_f >> vdc.phi_f;
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

void InDetectors(ifstream& datafile, vector<Detector>& detectors/*,Int_t& detector_index */){
  Detector temp;
  datafile >> temp.detector_num;
  // only add detector_num here, info of VDC and Paddles are in InVDC and InPaddle function.
  detectors.push_back(temp);
  datafile.ignore(200,'\n');
  cout << "A dectector array number "<< temp.detector_num <<" has been read" << endl;
}

void ProcessEvent(ifstream& datafile, BH_Event* event, bool& newevent, string& linetitle, Int_t& detector_index)
{
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
    InDetectors(datafile, event->detectors);
    detector_index++;
  }
  else if(linetitle=="VDC:"){ 
    InVDC(datafile, event->detectors.at(detector_index).vdc);
  }
  else if(linetitle=="Paddle:"){
    cout << "Adding paddles for detector array index number: " << detector_index << endl;
    InPaddle(datafile, event->detectors.at(detector_index).paddles);
    //datafile.ignore(200,'\n');
  }
  else if(linetitle=="Event:"){
    newevent=true;
    //if(newevent==true)
      //cout << "newevent has been changed to true." << endl;
  }
  //else datafile.ignore(200,'\n');
}

void MakeTree3(void){
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
    event->detectors.clear();
    Int_t detector_index=-1;

    while(datafile && newevent==false){
      ProcessEvent(datafile, event, newevent,linetitle, detector_index);
      if(!datafile) break;
    }
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
